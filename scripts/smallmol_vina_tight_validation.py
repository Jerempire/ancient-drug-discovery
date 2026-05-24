"""
Tight-box multi-seed validation of top 5 Pocket2Mol-Vina hits.

What changed vs previous counterscreens:
  - ERAP1 receptor: erap1_2YD0_experimental.pdb (with Zn HETATM) instead of AF model
  - ERAP1 box center: actual Zn coords (6.512, 71.495, 21.017) instead of CA centroid
  - IRAP box center: actual Zn coords (9.347, -31.132, -22.851) instead of CA centroid
  - Smaller box: 18 A instead of 22 A (tighter focus on the Zn pocket)
  - Higher exhaustiveness: 16 instead of 8
  - 5 random seeds per dock; report median + IQR
  - ERAP2 box: P02 centroid unchanged (we trust it; it's the target by design)

Goal: see whether the -2.39 margin on top hit survives a tighter, less-noisy comparison.
"""
from __future__ import annotations

import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/results/smallmol"
WORK = DATA / "vina_tight_validation"
WORK.mkdir(parents=True, exist_ok=True)
OUT = DATA / "vina_tight_validation.json"

VINA = Path(r"C:\Users\Jeremy\vina.exe")
OBABEL = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts\obabel.exe")

ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = ROOT / "data/structures/erap1_2YD0_experimental.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"

# Centers
P02_CENTER = [-19.05, 9.19, 11.49]           # ERAP2 allosteric pocket (unchanged)
ERAP1_ZN = [6.512, 71.495, 21.017]            # ERAP1 catalytic Zn from 2YD0
IRAP_ZN = [9.347, -31.132, -22.851]           # IRAP catalytic Zn from 5MJ6, chain A

BOX_TIGHT = 18.0
EXH = 16
N_POSES = 3
SEEDS = [1, 7, 42, 123, 2026]
SEL_THRESHOLD = -1.4  # ~10x

# Top-5 candidates from the large-run counterscreen, by margin
CANDIDATES = [
    ("pyridine-azabicyclic (champion)",      "CC12Cc3ccncc3C(C1)C(=O)C1CCCC=C12"),
    ("decalone-pyrrolidinone",                "CC1C2CCC(O)C1N1C(=O)CCCC1C2"),
    ("salicylate diacid-amide",               "CCC(C(=O)NC(C)=O)c1cccc(O)c1C(=O)O"),
    ("diphenyl-lactam",                       "CCC1(c2ccccc2)NC(=O)Cc2ccccc21"),
    ("bridged ene-lactone",                   "CC1=C2CCC(=O)C3CC(=O)OC(C(O)C1O)C23C"),
    ("tetrahydronaphtho-pyranone-A",          "CC(C)C1CCC2(C)c3cccc(O)c3OC(=O)C12"),
    ("polycyclic acetamide (original lead)",  "CC(=O)N1CC2CC(=O)C3C2(C)CCC31C(=O)O"),
]


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


def clean_pdb(input_pdb, out_pdb, keep_chain="A"):
    with open(input_pdb) as f, open(out_pdb, "w") as o:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                if len(line) >= 22 and line[21] == keep_chain and line[16] in (" ", "A"):
                    o.write(line)
            elif line.startswith("END"):
                o.write(line)


def prep_receptor(input_pdb, out_pdbqt, keep_chain="A"):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000:
        return
    clean = WORK / (input_pdb.stem + "_clean.pdb")
    clean_pdb(input_pdb, clean, keep_chain)
    cmd = [str(OBABEL), "-ipdb", str(clean), "-opdbqt", "-O", str(out_pdbqt),
           "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0 or out_pdbqt.stat().st_size < 5000:
        raise RuntimeError(f"obabel failed on {input_pdb}: {r.stderr[-300:]}")


def dock(rec, lig, center, out, seed):
    cmd = [str(VINA), "--receptor", str(rec), "--ligand", str(lig),
           "--center_x", f"{center[0]}", "--center_y", f"{center[1]}", "--center_z", f"{center[2]}",
           "--size_x", f"{BOX_TIGHT}", "--size_y", f"{BOX_TIGHT}", "--size_z", f"{BOX_TIGHT}",
           "--exhaustiveness", str(EXH), "--num_modes", str(N_POSES), "--cpu", "4",
           "--seed", str(seed), "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        return None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m: return float(m.group(1))
    return None


def median_iqr(vals):
    if not vals: return None, None, None
    s = sorted(vals)
    n = len(s)
    med = statistics.median(s)
    q1 = s[max(0, n // 4)]
    q3 = s[min(n - 1, (3 * n) // 4)]
    return round(med, 2), round(q1, 2), round(q3, 2)


def main():
    print("Tight-box multi-seed validation of top candidates")
    print(f"  ERAP2 P02 box:  center={P02_CENTER}  size={BOX_TIGHT}")
    print(f"  ERAP1 Zn  box:  center={ERAP1_ZN}  size={BOX_TIGHT}  (2YD0)")
    print(f"  IRAP  Zn  box:  center={IRAP_ZN}  size={BOX_TIGHT}  (5MJ6)")
    print(f"  exh={EXH}, {len(SEEDS)} seeds: {SEEDS}")

    e2_r = WORK / "erap2_receptor.pdbqt"
    e1_r = WORK / "erap1_2YD0_receptor.pdbqt"
    ir_r = WORK / "irap_5MJ6_receptor.pdbqt"
    print("\n=== Receptor prep ===")
    prep_receptor(ERAP2_PDB, e2_r)
    prep_receptor(ERAP1_PDB, e1_r)
    prep_receptor(IRAP_PDB, ir_r)
    print("  ready")

    out_rows = []
    for i, (label, smi) in enumerate(CANDIDATES, 1):
        print(f"\n--- [{i}/{len(CANDIDATES)}] {label} ---")
        print(f"  SMILES: {smi}")
        lig = WORK / f"lig_{i:02d}.pdbqt"
        if not smi_to_pdbqt(smi, lig):
            print("  PREP FAIL"); continue

        rows = {"label": label, "smiles": smi, "by_seed": {}}
        for target_name, rec, center in [
            ("erap2_p02", e2_r, P02_CENTER),
            ("erap1_zn", e1_r, ERAP1_ZN),
            ("irap_zn", ir_r, IRAP_ZN),
        ]:
            vals = []
            for s in SEEDS:
                out_path = WORK / f"out_{i:02d}_{target_name}_s{s}.pdbqt"
                dg = dock(rec, lig, center, out_path, s)
                if dg is not None:
                    vals.append(dg)
            med, q1, q3 = median_iqr(vals)
            print(f"  {target_name:>10}: median={med}  IQR=[{q1}, {q3}]  n={len(vals)}/{len(SEEDS)}  vals={[round(v,2) for v in vals]}")
            rows["by_seed"][target_name] = {"median": med, "q1": q1, "q3": q3,
                                            "values": [round(v, 2) for v in vals],
                                            "n": len(vals)}

        # Selectivity using medians; worst-case (using q1 of off-target, q3 of ERAP2) as stress test
        e2_med = rows["by_seed"]["erap2_p02"]["median"]
        e1_med = rows["by_seed"]["erap1_zn"]["median"]
        ir_med = rows["by_seed"]["irap_zn"]["median"]
        if None not in (e2_med, e1_med, ir_med):
            worst_off = max(e1_med, ir_med)
            margin = round(e2_med - worst_off, 2)
            # Worst case (95% CI-ish): ERAP2 is at its weakest (q3, less negative), off-target at its strongest (q1, more negative)
            e2_q3 = rows["by_seed"]["erap2_p02"]["q3"]
            e1_q1 = rows["by_seed"]["erap1_zn"]["q1"]
            ir_q1 = rows["by_seed"]["irap_zn"]["q1"]
            worst_off_strong = min(e1_q1, ir_q1)
            margin_worst = round(e2_q3 - worst_off_strong, 2)
            rows["margin_median"] = margin
            rows["margin_worst_case"] = margin_worst
            rows["selective_median"] = margin <= SEL_THRESHOLD
            rows["selective_worst_case"] = margin_worst <= SEL_THRESHOLD
            print(f"  ==> margin_median = {margin:+.2f}  ({'SEL' if rows['selective_median'] else 'NOT SEL'})")
            print(f"  ==> margin_worst  = {margin_worst:+.2f}  ({'SEL' if rows['selective_worst_case'] else 'NOT SEL'})")
        out_rows.append(rows)

    # Summary
    print()
    print("=" * 100)
    print("SUMMARY — tight-box multi-seed validation")
    print("=" * 100)
    print(f"{'Label':<40}{'E2med':>8}{'E1med':>8}{'IRmed':>8}{'mar_med':>10}{'mar_wc':>10}  Result")
    print("-" * 100)
    for r in out_rows:
        if "margin_median" not in r:
            print(f"{r['label']:<40}  --skipped--")
            continue
        e2 = r["by_seed"]["erap2_p02"]["median"]
        e1 = r["by_seed"]["erap1_zn"]["median"]
        ir = r["by_seed"]["irap_zn"]["median"]
        sm = r["margin_median"]; sw = r["margin_worst_case"]
        verdict = "**HOLDS**" if r["selective_median"] else "lost"
        wc = " (worst-case still SEL)" if r["selective_worst_case"] else ""
        print(f"{r['label']:<40}{e2:>8.2f}{e1:>8.2f}{ir:>8.2f}{sm:>+10.2f}{sw:>+10.2f}  {verdict}{wc}")

    holds_medium = sum(1 for r in out_rows if r.get("selective_median"))
    holds_strict = sum(1 for r in out_rows if r.get("selective_worst_case"))
    print()
    print(f"Hits surviving by median margin:      {holds_medium}/{len(out_rows)}")
    print(f"Hits surviving worst-case (q3 vs q1): {holds_strict}/{len(out_rows)}")

    payload = {
        "validation_protocol": {
            "exhaustiveness": EXH, "seeds": SEEDS, "box_A": BOX_TIGHT,
            "erap2_box_center": P02_CENTER, "erap1_box_center": ERAP1_ZN, "irap_box_center": IRAP_ZN,
            "erap1_receptor": "2YD0 experimental (with Zn)",
            "irap_receptor": "5MJ6 experimental (with Zn)",
            "selectivity_threshold_kcal_mol": SEL_THRESHOLD,
        },
        "n_candidates": len(out_rows),
        "n_holds_by_median": holds_medium,
        "n_holds_by_worst_case": holds_strict,
        "results": out_rows,
    }
    with open(OUT, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
