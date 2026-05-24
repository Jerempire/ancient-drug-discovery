"""Run the dock+counterscreen pipeline on the 561-molecule Pocket2Mol-large output.

Pulls SMILES from SMILES.txt directly (skips the heavy .pt files). Then:
  1. RDKit lead-like filter
  2. Vina rescore vs ERAP2 P02
  3. Counterscreen top-N vs ERAP1 + IRAP
  4. Selectivity gate at margin <= -1.4 kcal/mol
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
P2M_DIR = next((DATA / "pocket2mol_large").glob("sample_p02_large_*"))
SMILES_TXT = P2M_DIR / "SMILES.txt"
WORK = DATA / "vina_work_large"
WORK.mkdir(parents=True, exist_ok=True)
OUT_CANDIDATES = DATA / "pocket2mol_large_candidates.json"
OUT_RESCORE = DATA / "pocket2mol_large_vina_rescore.json"
OUT_COUNTER = DATA / "pocket2mol_large_counterscreen.json"

VINA = Path(r"C:\Users\Jeremy\vina.exe")
OBABEL = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts\obabel.exe")

P02 = [-19.05, 9.19, 11.49]
E1 = [-8.45, 3.69, -12.46]
IR = [6.56, -28.94, -23.53]
BOX = 22.0
EXH = 8
N_POSES = 3
SEL = -1.4
N_COUNTER = 40  # counterscreen top N

ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = ROOT / "data/structures/erap1_wt_alphafold.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"


def score(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol); hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol); rotb = Descriptors.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    lipinski = hbd <= 5 and hba <= 10 and mw <= 500 and logp <= 5
    veber = rotb <= 10 and tpsa <= 140
    return {"smiles": smi, "mw": round(mw, 2), "logp": round(logp, 2),
            "hbd": hbd, "hba": hba, "tpsa": round(tpsa, 2),
            "rotb": rotb, "rings": rings,
            "lipinski": lipinski, "veber": veber,
            "lead_like": lipinski and veber and 150 <= mw <= 450 and rings >= 1}


def smi_to_pdbqt(smi, out_pdbqt):
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


def prep_receptor(input_pdb, out_pdbqt, keep_chain="A"):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000: return
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
    subprocess.run(cmd, capture_output=True, text=True, timeout=180, check=True)


def dock(rec, lig, center, out):
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
    print(f"Pocket2Mol-large SMILES.txt: {SMILES_TXT}")
    raw = [s.strip() for s in SMILES_TXT.read_text().splitlines() if s.strip()]
    seen = set(); smis = []
    for s in raw:
        m = Chem.MolFromSmiles(s)
        if m is None: continue
        c = Chem.MolToSmiles(m, canonical=True)
        if c in seen: continue
        seen.add(c); smis.append(c)
    print(f"Unique parseable: {len(smis)}  (raw: {len(raw)})")

    scored = [s for s in (score(x) for x in smis) if s]
    lead = [s for s in scored if s["lead_like"]]
    print(f"Valid: {len(scored)}  |  Lead-like: {len(lead)}")
    with open(OUT_CANDIDATES, "w") as f:
        json.dump({"all": scored, "lead_like": lead, "pocket_centroid": P02}, f, indent=2)

    print("\n=== Receptor prep ===")
    e2_r = WORK / "erap2_receptor.pdbqt"
    e1_r = WORK / "erap1_receptor.pdbqt"
    ir_r = WORK / "irap_receptor.pdbqt"
    prep_receptor(ERAP2_PDB, e2_r); prep_receptor(ERAP1_PDB, e1_r); prep_receptor(IRAP_PDB, ir_r)
    print("  ready")

    print(f"\n=== Vina rescore ERAP2 P02 ({len(lead)} lead-like) ===")
    results = []
    for i, c in enumerate(lead, 1):
        lig = WORK / f"lig_{i:04d}.pdbqt"
        if not smi_to_pdbqt(c["smiles"], lig):
            continue
        dg = dock(e2_r, lig, P02, WORK / f"o2_{i:04d}.pdbqt")
        if dg is None: continue
        results.append({**c, "dG_erap2_p02": dg, "_idx": i})
        if i % 20 == 0:
            print(f"  [{i}/{len(lead)}]  dG={dg:>6.2f}  best so far: {min(r['dG_erap2_p02'] for r in results):>6.2f}")

    results.sort(key=lambda r: r["dG_erap2_p02"])
    print(f"\nDocked: {len(results)}/{len(lead)}  best dG: {results[0]['dG_erap2_p02']:.2f}")
    with open(OUT_RESCORE, "w") as f:
        json.dump({"n_docked": len(results), "all_results": results, "top20": results[:20]}, f, indent=2)

    n_cs = min(N_COUNTER, len(results))
    print(f"\n=== Counterscreen top {n_cs} vs ERAP1 + IRAP ===")
    counter = []
    for j, c in enumerate(results[:n_cs], 1):
        lig = WORK / f"lig_{c['_idx']:04d}.pdbqt"
        d1 = dock(e1_r, lig, E1, WORK / f"o1_{c['_idx']:04d}.pdbqt")
        dr = dock(ir_r, lig, IR, WORK / f"or_{c['_idx']:04d}.pdbqt")
        if d1 is None or dr is None:
            continue
        margin = c["dG_erap2_p02"] - max(d1, dr)
        sel = margin <= SEL
        flag = "**SEL**" if sel else ""
        counter.append({**c, "dG_erap1_active": d1, "dG_irap_active": dr,
                       "selectivity_margin": round(margin, 2), "selective": sel})
        print(f"  [{j:>2}]  E2={c['dG_erap2_p02']:>6.2f}  E1={d1:>6.2f}  IR={dr:>6.2f}  "
              f"margin={margin:>+6.2f}  {flag}  {c['smiles'][:55]}")

    counter.sort(key=lambda r: r["selectivity_margin"])
    selective = [r for r in counter if r["selective"]]

    print()
    print("=" * 80)
    print(f"SELECTIVE HITS: {len(selective)} / {len(counter)}")
    print("=" * 80)
    for r in selective:
        print(f"  margin={r['selectivity_margin']:>+6.2f}  E2={r['dG_erap2_p02']:>6.2f}  "
              f"E1={r['dG_erap1_active']:>6.2f}  IR={r['dG_irap_active']:>6.2f}  "
              f"MW={r['mw']:>5.0f}  {r['smiles']}")

    with open(OUT_COUNTER, "w") as f:
        json.dump({"selectivity_threshold_kcal_mol": SEL,
                  "n_counterscreened": len(counter),
                  "n_selective": len(selective),
                  "selective_hits": selective,
                  "all_by_margin": counter}, f, indent=2)
    print(f"\nSaved: {OUT_COUNTER}")


if __name__ == "__main__":
    main()
