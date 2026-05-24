"""
Full pipeline for Pocket2Mol P02 candidates:
  1. Extract SMILES from samples_*.pt files (rdkit Mols saved by Pocket2Mol)
  2. RDKit lead-like filter (Lipinski + Veber + MW 150-450 + ring)
  3. Vina rescore vs ERAP2 P02
  4. Counterscreen vs ERAP1 + IRAP active sites
  5. Rank by selectivity margin
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import pickle
from pathlib import Path

import torch
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
P2M_DIR = next((DATA / "pocket2mol").glob("sample_for_pdb_*"))
WORK = DATA / "vina_work_p2m"
WORK.mkdir(parents=True, exist_ok=True)
OUT_CANDIDATES = DATA / "pocket2mol_candidates.json"
OUT_RESCORE = DATA / "pocket2mol_vina_rescore.json"
OUT_COUNTER = DATA / "pocket2mol_counterscreen.json"

VINA_EXE = Path(r"C:\Users\Jeremy\vina.exe")
OBABEL = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts\obabel.exe")

P02_CENTER = [-19.05, 9.19, 11.49]
ERAP1_CENTER = [-8.45, 3.69, -12.46]
IRAP_CENTER = [6.56, -28.94, -23.53]
BOX_SIZE = 22.0
EXHAUSTIVENESS = 8
N_POSES = 3
SELECTIVITY_THRESHOLD = -1.4  # ~10x

ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = ROOT / "data/structures/erap1_wt_alphafold.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"


def extract_smiles_from_pt() -> list[str]:
    seen = set()
    smiles = []
    for pt in sorted(P2M_DIR.glob("samples_*.pt")):
        try:
            data = torch.load(pt, map_location="cpu", weights_only=False)
        except Exception as e:
            print(f"  load fail {pt.name}: {e}")
            continue
        # Pocket2Mol saves a dict or list; inspect
        if isinstance(data, dict):
            # Try common keys
            candidates = data.get("mol") or data.get("smiles") or data.get("rdmol") or []
            if isinstance(candidates, str):
                candidates = [candidates]
        elif isinstance(data, list):
            candidates = data
        else:
            candidates = [data]
        for c in candidates if hasattr(candidates, "__iter__") else [candidates]:
            smi = None
            if isinstance(c, str):
                smi = c
            elif hasattr(c, "rdmol"):
                try: smi = Chem.MolToSmiles(c.rdmol)
                except: pass
            elif isinstance(c, dict):
                smi = c.get("smiles") or c.get("smi")
                if not smi and "rdmol" in c:
                    try: smi = Chem.MolToSmiles(c["rdmol"])
                    except: pass
            else:
                try: smi = Chem.MolToSmiles(c)
                except: pass
            if smi:
                m = Chem.MolFromSmiles(smi)
                if m is None: continue
                canon = Chem.MolToSmiles(m, canonical=True)
                if canon in seen: continue
                seen.add(canon)
                smiles.append(canon)
    return smiles


def extract_smiles_from_log() -> list[str]:
    """Fallback: parse the run log for 'Success: <smiles>' lines."""
    log = P2M_DIR / "log.txt"
    if not log.exists(): return []
    seen, smiles = set(), []
    with open(log, encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = re.search(r"Success:\s+(\S+)", line)
            if not m: continue
            smi = m.group(1)
            mol = Chem.MolFromSmiles(smi)
            if mol is None: continue
            canon = Chem.MolToSmiles(mol, canonical=True)
            if canon in seen: continue
            seen.add(canon)
            smiles.append(canon)
    return smiles


def score_drug_like(smi: str) -> dict | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    lipinski = hbd <= 5 and hba <= 10 and mw <= 500 and logp <= 5
    veber = rotb <= 10 and tpsa <= 140
    return {
        "smiles": smi, "mw": round(mw, 2), "logp": round(logp, 2),
        "hbd": hbd, "hba": hba, "tpsa": round(tpsa, 2),
        "rotb": rotb, "rings": rings,
        "lipinski": lipinski, "veber": veber,
        "lead_like": lipinski and veber and 150 <= mw <= 450 and rings >= 1,
    }


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


def prep_receptor(input_pdb: Path, out_pdbqt: Path, clean_chain="A"):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000:
        return
    clean_pdb = WORK / (input_pdb.stem + "_clean.pdb")
    with open(input_pdb) as f, open(clean_pdb, "w") as out:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                if len(line) >= 22 and line[21] == clean_chain and line[16] in (" ", "A"):
                    out.write(line)
            elif line.startswith("END"):
                out.write(line)
    cmd = [str(OBABEL), "-ipdb", str(clean_pdb), "-opdbqt", "-O", str(out_pdbqt),
           "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0 or out_pdbqt.stat().st_size < 5000:
        raise RuntimeError(f"obabel failed on {input_pdb}: {r.stderr[-300:]}")


def vina_dock(receptor: Path, ligand: Path, center, out: Path) -> float | None:
    cmd = [str(VINA_EXE),
           "--receptor", str(receptor), "--ligand", str(ligand),
           "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}",
           "--center_z", f"{center[2]:.3f}",
           "--size_x", f"{BOX_SIZE:.1f}", "--size_y", f"{BOX_SIZE:.1f}",
           "--size_z", f"{BOX_SIZE:.1f}",
           "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(N_POSES),
           "--cpu", "4", "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0: return None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m: return float(m.group(1))
    return None


def main():
    print(f"P2M dir: {P2M_DIR}")
    smis = extract_smiles_from_pt()
    if not smis:
        print("  .pt extraction yielded nothing, falling back to log parser")
        smis = extract_smiles_from_log()
    print(f"Unique SMILES extracted: {len(smis)}")

    scored = [s for s in (score_drug_like(smi) for smi in smis) if s]
    lead = [s for s in scored if s["lead_like"]]
    print(f"Valid: {len(scored)}  |  Lead-like: {len(lead)}")
    with open(OUT_CANDIDATES, "w") as f:
        json.dump({"all": scored, "lead_like": lead, "pocket_centroid": P02_CENTER}, f, indent=2)
    print(f"Saved: {OUT_CANDIDATES}")

    print("\n=== Receptor prep ===")
    e2 = WORK / "erap2_receptor.pdbqt"
    e1 = WORK / "erap1_receptor.pdbqt"
    irp = WORK / "irap_receptor.pdbqt"
    prep_receptor(ERAP2_PDB, e2)
    prep_receptor(ERAP1_PDB, e1)
    prep_receptor(IRAP_PDB, irp)
    print(f"  receptors ready")

    print(f"\n=== Vina rescore ERAP2 P02 ===  ({len(lead)} lead-like)")
    results = []
    for i, c in enumerate(lead):
        lig = WORK / f"lig_{i:03d}.pdbqt"
        if not smi_to_pdbqt(c["smiles"], lig):
            print(f"  [{i+1:>3}/{len(lead)}] prep fail  {c['smiles'][:60]}")
            continue
        dg = vina_dock(e2, lig, P02_CENTER, WORK / f"out_p02_{i:03d}.pdbqt")
        if dg is None:
            print(f"  [{i+1:>3}/{len(lead)}] dock fail  {c['smiles'][:60]}")
            continue
        results.append({**c, "dG_erap2_p02": dg, "_lig_idx": i})
        print(f"  [{i+1:>3}/{len(lead)}]  ERAP2 dG={dg:>6.2f}  MW={c['mw']:>5.0f}  {c['smiles'][:60]}")

    results.sort(key=lambda r: r["dG_erap2_p02"])
    with open(OUT_RESCORE, "w") as f:
        json.dump({"target": "ERAP2 P02", "n_docked": len(results),
                  "all_results": results, "top10": results[:10]}, f, indent=2)
    print(f"Saved: {OUT_RESCORE}")

    print(f"\n=== Counterscreen TOP 20 vs ERAP1 + IRAP ===")
    top = results[:20]
    counter = []
    for j, c in enumerate(top, 1):
        lig = WORK / f"lig_{c['_lig_idx']:03d}.pdbqt"
        e1dg = vina_dock(e1, lig, ERAP1_CENTER, WORK / f"out_e1_{c['_lig_idx']:03d}.pdbqt")
        irdg = vina_dock(irp, lig, IRAP_CENTER, WORK / f"out_ir_{c['_lig_idx']:03d}.pdbqt")
        if e1dg is None or irdg is None:
            print(f"  [{j:>2}] counterscreen fail e1={e1dg} ir={irdg}")
            continue
        worst = max(e1dg, irdg)
        margin = c["dG_erap2_p02"] - worst
        sel = margin <= SELECTIVITY_THRESHOLD
        flag = "**SELECTIVE**" if sel else ""
        print(f"  [{j:>2}]  E2={c['dG_erap2_p02']:>6.2f}  E1={e1dg:>6.2f}  IR={irdg:>6.2f}  "
              f"margin={margin:>+6.2f}  {flag}  {c['smiles'][:50]}")
        counter.append({**c, "dG_erap1_active": e1dg, "dG_irap_active": irdg,
                       "selectivity_margin": round(margin, 2), "selective": sel})

    counter.sort(key=lambda r: r["selectivity_margin"])
    selective = [r for r in counter if r["selective"]]
    print(f"\nSELECTIVE HITS (margin <= {SELECTIVITY_THRESHOLD} kcal/mol): "
          f"{len(selective)} / {len(counter)}")
    if selective:
        print("Recommended for wet-lab triage:")
        for r in selective:
            print(f"  margin={r['selectivity_margin']:>+6.2f}  E2={r['dG_erap2_p02']:>6.2f} "
                  f"E1={r['dG_erap1_active']:>6.2f} IR={r['dG_irap_active']:>6.2f}  {r['smiles']}")
    else:
        print("No selective hits at this threshold. Top 5 by margin:")
        for r in counter[:5]:
            print(f"  margin={r['selectivity_margin']:>+6.2f}  {r['smiles']}")

    with open(OUT_COUNTER, "w") as f:
        json.dump({"selectivity_threshold_kcal_mol": SELECTIVITY_THRESHOLD,
                  "n_selective": len(selective),
                  "selective_hits": selective,
                  "all_results_by_margin": counter}, f, indent=2)
    print(f"Saved: {OUT_COUNTER}")


if __name__ == "__main__":
    main()
