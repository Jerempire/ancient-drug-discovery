"""
Vina rescore of DrugGPT lead-like candidates against ERAP2 P02 3D pocket.

Tooling on this Windows box:
  - vina binary: C:\\Users\\Jeremy\\vina.exe (AutoDock Vina v1.2.5)
  - meeko CLI: mk_prepare_receptor.exe + mk_prepare_ligand.exe (in env Scripts/)
  - PDBFixer for H-addition on receptor

Pipeline per ligand:
  SMILES -> RDKit 3D embed (ETKDG + UFF) -> SDF -> meeko ligand PDBQT
  -> vina.exe dock against pre-prepared ERAP2 P02 PDBQT
  -> parse best-pose affinity
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
CANDIDATES = ROOT / "data/results/smallmol/druggpt_p02_candidates.json"
OUT = ROOT / "data/results/smallmol/vina_rescore_p02.json"
WORK = ROOT / "data/results/smallmol/vina_work"
WORK.mkdir(parents=True, exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)

VINA_EXE = Path(r"C:\Users\Jeremy\vina.exe")
ENV_SCRIPTS = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts")
OBABEL = ENV_SCRIPTS / "obabel.exe"

BOX_SIZE = 22.0
EXHAUSTIVENESS = 8
N_POSES = 3
CPU = 4


def prep_receptor(input_pdb: Path, out_pdbqt: Path):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000:
        print(f"  receptor cached: {out_pdbqt}")
        return out_pdbqt
    # obabel: PDB -> PDBQT rigid receptor, add Hs at pH 7.4, compute Gasteiger charges
    print("  receptor prep: obabel PDB -> PDBQT (rigid, pH 7.4) ...")
    cmd = [str(OBABEL), "-ipdb", str(input_pdb), "-opdbqt", "-O", str(out_pdbqt),
           "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0 or not out_pdbqt.exists() or out_pdbqt.stat().st_size < 5000:
        print(r.stdout); print(r.stderr)
        raise RuntimeError("obabel receptor prep failed")
    print(f"  receptor PDBQT: {out_pdbqt} ({out_pdbqt.stat().st_size:,} bytes)")
    return out_pdbqt


def smiles_to_pdbqt(smi: str, out_pdbqt: Path) -> bool:
    """SMILES -> RDKit 3D embed -> meeko PDBQT (vina-compatible torsion tree)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    try:
        if AllChem.EmbedMolecule(mol, randomSeed=42) < 0:
            return False
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        return False
    try:
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            return False
        pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return False
        out_pdbqt.write_text(pdbqt_str)
        return out_pdbqt.stat().st_size > 200
    except Exception:
        return False


def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path, center, out_pdbqt: Path,
            box_size=BOX_SIZE, exh=EXHAUSTIVENESS, n_poses=N_POSES, cpu=CPU):
    cmd = [
        str(VINA_EXE),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{box_size:.1f}",
        "--size_y", f"{box_size:.1f}",
        "--size_z", f"{box_size:.1f}",
        "--exhaustiveness", str(exh),
        "--num_modes", str(n_poses),
        "--cpu", str(cpu),
        "--out", str(out_pdbqt),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        return None, r.stderr[-300:] if r.stderr else "vina failed"
    # Parse best affinity from stdout
    # Vina prints a table:   mode |   affinity | ...
    #                          1     -7.123
    best = None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            best = float(m.group(1))
            break
    return best, None


def main():
    if not VINA_EXE.exists():
        sys.exit(f"vina not found: {VINA_EXE}")
    if not OBABEL.exists():
        sys.exit(f"obabel not found: {OBABEL}")

    with open(CANDIDATES) as f:
        cands_data = json.load(f)
    p02_center = cands_data["pocket_centroid"]
    cands = [c for c in cands_data["candidates"] if c["lead_like"]]
    print(f"P02 centroid: {p02_center}")
    print(f"Lead-like candidates: {len(cands)}")

    rec_pdbqt = WORK / "erap2_receptor.pdbqt"
    print("\n=== Receptor prep ===")
    prep_receptor(ERAP2_PDB, rec_pdbqt)

    print("\n=== Docking ===")
    results = []
    fails = []
    for i, c in enumerate(cands):
        smi = c["smiles"]
        lig_pdbqt = WORK / f"lig_{i:03d}.pdbqt"
        out_pdbqt = WORK / f"out_{i:03d}.pdbqt"
        if not smiles_to_pdbqt(smi, lig_pdbqt):
            print(f"  [{i+1:>2}/{len(cands)}] LIG-PREP FAIL  {smi[:60]}")
            fails.append({"smiles": smi, "stage": "lig_prep"})
            continue
        try:
            dg, err = run_vina(rec_pdbqt, lig_pdbqt, p02_center, out_pdbqt)
        except subprocess.TimeoutExpired:
            print(f"  [{i+1:>2}/{len(cands)}] DOCK TIMEOUT   {smi[:60]}")
            fails.append({"smiles": smi, "stage": "dock_timeout"})
            continue
        if dg is None:
            print(f"  [{i+1:>2}/{len(cands)}] DOCK FAIL      {smi[:60]}  ({err[:100] if err else ''})")
            fails.append({"smiles": smi, "stage": "dock", "err": err[:200] if err else ""})
            continue
        results.append({**c, "vina_dG_kcal_mol": round(dg, 2)})
        print(f"  [{i+1:>2}/{len(cands)}]  dG={dg:>6.2f}  MW={c['mw']:>5.0f}  {smi[:60]}")

    results.sort(key=lambda r: r["vina_dG_kcal_mol"])
    print(f"\nDocked {len(results)} / {len(cands)} successfully ({len(fails)} fails)")

    print("\n=== TOP 10 by Vina dG ===")
    print(f"{'Rank':<5}{'dG':>7}{'MW':>7}{'logP':>7}{'Rings':>7}  SMILES")
    print("-" * 90)
    for i, r in enumerate(results[:10], 1):
        print(f"{i:<5}{r['vina_dG_kcal_mol']:>6.2f}  {r['mw']:>6.1f}{r['logp']:>7.2f}{r['rings']:>7}  {r['smiles']}")

    out = {
        "target": "ERAP2 P02 allosteric pocket",
        "receptor": str(rec_pdbqt),
        "center": p02_center,
        "box_size_A": BOX_SIZE,
        "exhaustiveness": EXHAUSTIVENESS,
        "n_candidates_in": len(cands),
        "n_docked": len(results),
        "n_failed": len(fails),
        "top10": results[:10],
        "all_results": results,
        "failures": fails,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
