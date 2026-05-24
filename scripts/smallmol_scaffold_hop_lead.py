"""
Scaffold-hop around the Pocket2Mol-Vina lead for ERAP2 P02:
  CC(=O)N1CC2CC(=O)C3C2(C)CCC31C(=O)O   (margin -1.89, ~25x selective)

Strategy:
  - Apply RDKit bioisosteric replacements at the COOH and acetamide handles.
  - Generate ~50 analogs, score drug-like, dock vs ERAP2 P02 + counterscreen vs ERAP1+IRAP.

Bioisosteres:
  - COOH -> tetrazole, sulfonamide, hydroxamate, acyl-sulfonamide
  - acetamide -> trifluoroacetamide, methylcarbamate, methanesulfonamide
  - Ring sizes 5<->6 not tried (would change the rigid core)
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Lipinski, Crippen, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/results/smallmol"
WORK = DATA / "vina_work_hop"
WORK.mkdir(parents=True, exist_ok=True)
OUT = DATA / "scaffold_hop_results.json"

VINA = Path(r"C:\Users\Jeremy\vina.exe")
OBABEL = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts\obabel.exe")

LEAD = "CC(=O)N1CC2CC(=O)C3C2(C)CCC31C(=O)O"
P02 = [-19.05, 9.19, 11.49]
E1 = [-8.45, 3.69, -12.46]
IR = [6.56, -28.94, -23.53]
BOX = 22.0
EXH = 8
N_POSES = 3

ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = ROOT / "data/structures/erap1_wt_alphafold.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"

# Bioisosteric replacements: (SMARTS pattern in lead, list of replacement SMILES fragments)
# Each fragment uses [*] for the attachment point to the rest of the molecule.
COOH_ISOSTERES = [
    ("C(=O)O", "C(=O)O", "COOH (parent)"),
    ("C(=O)O", "c1nnn[nH]1", "tetrazole"),
    ("C(=O)O", "S(=O)(=O)N", "sulfonamide"),
    ("C(=O)O", "C(=O)NO", "hydroxamate"),
    ("C(=O)O", "C(=O)NS(C)(=O)=O", "acyl-methanesulfonamide"),
    ("C(=O)O", "C(=O)NCC#N", "cyanoacetamide"),
    ("C(=O)O", "C#N", "nitrile"),
    ("C(=O)O", "C(N)=O", "primary amide"),
]
ACETYL_REPLACEMENTS = [
    ("CC(=O)", "CC(=O)", "acetyl (parent)"),
    ("CC(=O)", "FC(F)(F)C(=O)", "trifluoroacetyl"),
    ("CC(=O)", "CC(=O)OC", "methylcarbamate (rev)"),  # actually a carbamate-style mod
    ("CC(=O)", "CS(=O)(=O)", "methanesulfonyl"),
    ("CC(=O)", "C(=O)", "formyl"),
    ("CC(=O)", "CCC(=O)", "propionyl"),
]


def enumerate_analogs(lead_smi: str) -> list[tuple[str, str]]:
    """Apply each (COOH replacement) x (acetyl replacement) combination via simple string sub.
    Returns (description, smiles) pairs, unique."""
    out = []
    seen = set()
    # The lead string contains "C(=O)O" once (the COOH) and "CC(=O)N" once (the acetamide -> N).
    # We do simple substring replacement at one site at a time, then combine.
    for old_a, new_a, name_a in ACETYL_REPLACEMENTS:
        for old_c, new_c, name_c in COOH_ISOSTERES:
            modified = lead_smi
            # Replace acetyl portion first (left side: "CC(=O)N" stem)
            # The acetyl appears as "CC(=O)N1" in lead — replace "CC(=O)" with new_a only at that position
            if old_a != new_a:
                modified = modified.replace(old_a + "N1", new_a + "N1", 1)
            if old_c != new_c:
                # COOH at the end: "C(=O)O" terminal
                modified = re.sub(re.escape(old_c) + r"$", new_c, modified)
            mol = Chem.MolFromSmiles(modified)
            if mol is None:
                continue
            canon = Chem.MolToSmiles(mol, canonical=True)
            if canon in seen:
                continue
            seen.add(canon)
            label = f"{name_a} + {name_c}"
            out.append((label, canon))
    return out


def score(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol); hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol); rotb = Descriptors.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    return {"smiles": smi, "mw": round(mw, 2), "logp": round(logp, 2),
            "hbd": hbd, "hba": hba, "tpsa": round(tpsa, 2),
            "rotb": rotb, "rings": rings,
            "lipinski": hbd <= 5 and hba <= 10 and mw <= 500 and logp <= 5,
            "veber": rotb <= 10 and tpsa <= 140}


def smi_to_pdbqt(smi: str, out_pdbqt: Path) -> bool:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return False
    mol = Chem.AddHs(mol)
    try:
        if AllChem.EmbedMolecule(mol, randomSeed=42) < 0: return False
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        setups = MoleculePreparation().prepare(mol)
        if not setups: return False
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok: return False
        out_pdbqt.write_text(pdbqt)
        return out_pdbqt.stat().st_size > 200
    except Exception:
        return False


def prep_receptor(input_pdb: Path, out_pdbqt: Path, keep_chain="A"):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000:
        return
    clean = WORK / (input_pdb.stem + "_clean.pdb")
    with open(input_pdb) as f, open(clean, "w") as o:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                if len(line) >= 22 and line[21] == keep_chain and line[16] in (" ", "A"):
                    o.write(line)
            elif line.startswith("END"):
                o.write(line)
    cmd = [str(OBABEL), "-ipdb", str(clean), "-opdbqt", "-O", str(out_pdbqt),
           "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError(f"obabel fail: {r.stderr[-200:]}")


def dock(rec: Path, lig: Path, center, out: Path) -> float | None:
    cmd = [str(VINA), "--receptor", str(rec), "--ligand", str(lig),
           "--center_x", f"{center[0]}", "--center_y", f"{center[1]}", "--center_z", f"{center[2]}",
           "--size_x", f"{BOX}", "--size_y", f"{BOX}", "--size_z", f"{BOX}",
           "--exhaustiveness", str(EXH), "--num_modes", str(N_POSES), "--cpu", "4",
           "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0: return None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m: return float(m.group(1))
    return None


def main():
    print(f"Lead: {LEAD}")
    analogs = enumerate_analogs(LEAD)
    print(f"Enumerated {len(analogs)} unique analogs")
    print(f"Including parent? {'yes' if any(a[1] == Chem.MolToSmiles(Chem.MolFromSmiles(LEAD)) for a in analogs) else 'no'}")

    scored = []
    for label, smi in analogs:
        s = score(smi)
        if not s: continue
        s["label"] = label
        scored.append(s)
    print(f"Valid: {len(scored)}")

    print("\n=== Receptor prep ===")
    e2 = WORK / "erap2_receptor.pdbqt"
    e1 = WORK / "erap1_receptor.pdbqt"
    ir = WORK / "irap_receptor.pdbqt"
    prep_receptor(ERAP2_PDB, e2); prep_receptor(ERAP1_PDB, e1); prep_receptor(IRAP_PDB, ir)
    print("  receptors ready")

    print(f"\n=== Triple dock {len(scored)} analogs (ERAP2 P02 + ERAP1 + IRAP) ===")
    results = []
    for i, c in enumerate(scored, 1):
        lig = WORK / f"hop_{i:03d}.pdbqt"
        if not smi_to_pdbqt(c["smiles"], lig):
            print(f"  [{i:>2}] prep fail   {c['label']:<40}  {c['smiles']}")
            continue
        d2 = dock(e2, lig, P02, WORK / f"o2_{i:03d}.pdbqt")
        d1 = dock(e1, lig, E1, WORK / f"o1_{i:03d}.pdbqt")
        dr = dock(ir, lig, IR, WORK / f"or_{i:03d}.pdbqt")
        if d2 is None or d1 is None or dr is None:
            print(f"  [{i:>2}] dock fail   {c['label']}")
            continue
        margin = d2 - max(d1, dr)
        sel = margin <= -1.4
        flag = "**SELECTIVE**" if sel else ""
        results.append({**c, "dG_erap2_p02": d2, "dG_erap1_active": d1,
                       "dG_irap_active": dr, "selectivity_margin": round(margin, 2),
                       "selective": sel})
        print(f"  [{i:>2}]  E2={d2:>6.2f}  E1={d1:>6.2f}  IR={dr:>6.2f}  "
              f"margin={margin:>+6.2f}  {flag}  {c['label']}")

    results.sort(key=lambda r: r["selectivity_margin"])
    selective = [r for r in results if r["selective"]]

    print()
    print("=" * 80)
    print(f"SELECTIVE ANALOGS: {len(selective)} / {len(results)}")
    print("=" * 80)
    if selective:
        for r in selective:
            print(f"  margin={r['selectivity_margin']:>+6.2f}  E2={r['dG_erap2_p02']:>6.2f}  "
                  f"{r['label']:<40}  {r['smiles']}")

    # SAR insight: did any modification beat the parent margin?
    parent = next((r for r in results if r["label"] == "acetyl (parent) + COOH (parent)"), None)
    if parent:
        print(f"\nParent margin: {parent['selectivity_margin']:>+6.2f}")
        better = [r for r in results if r["selectivity_margin"] < parent["selectivity_margin"]]
        print(f"Analogs beating parent: {len(better)}")
        for r in better[:5]:
            print(f"  {r['label']:<40}  margin={r['selectivity_margin']:>+6.2f}  "
                  f"(parent {parent['selectivity_margin']:>+6.2f})  {r['smiles']}")

    out_payload = {
        "lead_parent": LEAD,
        "n_analogs": len(scored),
        "n_docked": len(results),
        "n_selective": len(selective),
        "selective_analogs": selective,
        "all_results_by_margin": results,
    }
    with open(OUT, "w") as f:
        json.dump(out_payload, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
