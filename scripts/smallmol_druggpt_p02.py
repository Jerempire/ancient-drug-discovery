"""
DrugGPT candidate generation conditioned on the ERAP2 P02 allosteric pocket.

Strategy: feed the P02 lining-residue sequence context (with the divergent-vs-both
selectivity handles flagged) as the protein-side prompt to liyuesen/druggpt, generate
N SMILES, RDKit-validate, Lipinski filter, save ranked candidates.

This is NOT 3D-pocket-conditioned generation (DrugGPT is sequence-conditioned).
3D rescoring with AutoDock-Vina or Boltz-2 affinity happens downstream.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# Windows SSL fix: load OS trust store before any HTTPS lib is imported
try:
    import truststore
    truststore.inject_into_ssl()
except ImportError:
    pass

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
POCKET_JSON = ROOT / "data/results/smallmol/pocket_inventory.json"
ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
OUT = ROOT / "data/results/smallmol/druggpt_p02_candidates.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_erap2_seq() -> dict:
    """Return {resnum: aa1} from the AF ERAP2 PDB."""
    seq = {}
    with open(ERAP2_PDB) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
                rn = int(line[22:26])
                resname = line[17:20].strip()
                if resname in AA3TO1:
                    seq[rn] = AA3TO1[resname]
    return seq


def main():
    with open(POCKET_JSON) as f:
        inv = json.load(f)
    top = inv["top_allosteric"]
    lining = top["lining_residues"]
    handles = top["divergent_vs_both"]
    pocket_id = top["pocket_id"]
    centroid = top["centroid_xyz"]
    print(f"Target pocket: {pocket_id}  |  lining residues: {lining}")
    print(f"Selectivity handles (vs ERAP1 AND IRAP): {handles}")
    print(f"Centroid (for downstream docking): {centroid}")

    # DrugGPT — sequence-conditioned SMILES generator
    print("\nLoading liyuesen/druggpt ...")
    tok = AutoTokenizer.from_pretrained("liyuesen/druggpt")
    model = AutoModelForCausalLM.from_pretrained("liyuesen/druggpt")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device).eval()
    print(f"Loaded on {device}")

    # Build pocket-focused prompt from real PDB sequence: 35-residue window
    # around the dominant selectivity-handle cluster (residues 701-712)
    erap2_seq = load_erap2_seq()
    window_residues = list(range(686, 721))
    pocket_window = "".join(erap2_seq.get(r, "X") for r in window_residues)
    prompt = f"<|startoftext|><P>{pocket_window}<L>"
    print(f"Prompt: {prompt[:80]}... ({len(prompt)} chars)")

    n_target = 200
    print(f"\nGenerating {n_target} candidates (temp=0.9, top_k=50)...")

    all_smiles = []
    batch = 20
    while len(all_smiles) < n_target:
        inputs = tok([prompt] * batch, return_tensors="pt", padding=True).to(device)
        with torch.no_grad():
            out = model.generate(
                **inputs,
                max_new_tokens=128,
                do_sample=True,
                top_k=50,
                temperature=0.9,
                pad_token_id=tok.eos_token_id,
            )
        for ids in out:
            text = tok.decode(ids, skip_special_tokens=False)
            # Extract substring after <L> tag, up to <|endoftext|>
            m = re.search(r"<L>([^<]*)", text)
            if not m:
                continue
            raw = m.group(1).strip()
            all_smiles.append(raw)
        print(f"  generated {len(all_smiles)}")

    print(f"\nRaw candidates: {len(all_smiles)}")

    # RDKit validation + drug-likeness
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    results = []
    seen = set()
    for smi in all_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        canon = Chem.MolToSmiles(mol, canonical=True)
        if canon in seen:
            continue
        seen.add(canon)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        lipinski = hbd <= 5 and hba <= 10 and mw <= 500 and logp <= 5
        veber = rotb <= 10 and tpsa <= 140
        results.append({
            "smiles": canon,
            "mw": round(mw, 2),
            "logp": round(logp, 2),
            "hbd": hbd, "hba": hba, "tpsa": round(tpsa, 2),
            "rotb": rotb, "rings": rings,
            "lipinski": lipinski,
            "veber": veber,
            "lead_like": lipinski and veber and 150 <= mw <= 450 and rings >= 1,
        })

    valid = len(results)
    lead = sum(1 for r in results if r["lead_like"])
    print(f"\nValid + unique: {valid}/{len(all_smiles)} ({100*valid/len(all_smiles):.0f}%)")
    print(f"Lead-like (Lipinski+Veber+MW 150-450+ring): {lead}")

    # Rank by lead-like first, then MW closeness to 350
    results.sort(key=lambda r: (not r["lead_like"], abs(r["mw"] - 350)))
    print("\nTop 10 lead-like candidates:")
    for r in [r for r in results if r["lead_like"]][:10]:
        print(f"  MW={r['mw']:>6.1f}  logP={r['logp']:>5.2f}  rings={r['rings']}  {r['smiles']}")

    out = {
        "target_pocket": pocket_id,
        "pocket_centroid": centroid,
        "selectivity_handles_vs_both": handles,
        "prompt_window_residues": "ERAP2[685-720]",
        "n_generated": len(all_smiles),
        "n_valid_unique": valid,
        "n_lead_like": lead,
        "candidates": results,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
