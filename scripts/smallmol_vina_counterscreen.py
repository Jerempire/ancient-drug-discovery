"""
Counterscreen the top Vina hits against ERAP1 and IRAP active sites.

Pipeline:
  1. Load top-N (default 10) from vina_rescore_p02.json
  2. Re-dock the same already-prepared ligand PDBQTs against:
       - ERAP1 active site (H353/H357/E376 centroid)
       - IRAP active site (H464/H468/E487 centroid)
  3. Compute selectivity margin: ERAP2_dG - max(ERAP1_dG, IRAP_dG)
  4. A negative margin means selective for ERAP2 (more negative dG = tighter binding).
  5. Flag candidates with margin <= -1.4 kcal/mol (~10x selectivity).
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/results/smallmol"
WORK = DATA / "vina_work"
RESCORE = DATA / "vina_rescore_p02.json"
OUT = DATA / "vina_counterscreen.json"

VINA_EXE = Path(r"C:\Users\Jeremy\vina.exe")
OBABEL = Path(r"C:\Users\Jeremy\anaconda3\envs\gpu\Scripts\obabel.exe")

ERAP1_PDB = ROOT / "data/structures/erap1_wt_alphafold.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"

# Active-site centroids (mean of catalytic His/Glu Ca coords, computed from PDBs)
ERAP1_CENTER = [-8.45, 3.69, -12.46]   # H353, H357, E376
IRAP_CENTER = [6.56, -28.94, -23.53]   # H464, H468, E487

BOX_SIZE = 22.0
EXHAUSTIVENESS = 8
N_POSES = 3
CPU = 4
N_TOP = 10  # how many ERAP2 top hits to counterscreen

SELECTIVITY_THRESHOLD = -1.4  # ~10x affinity margin (kcal/mol)


def clean_pdb_for_obabel(input_pdb: Path, out_pdb: Path, keep_chain: str = "A"):
    """Strip HETATM, alt-locs, and non-target chains. obabel chokes on dirty crystal PDBs."""
    with open(input_pdb) as f, open(out_pdb, "w") as out:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                if len(line) >= 22 and line[21] == keep_chain:
                    # Skip alt-loc B/C
                    if len(line) >= 17 and line[16] in (" ", "A"):
                        out.write(line)
            elif line.startswith("END"):
                out.write(line)


def prep_receptor(input_pdb: Path, out_pdbqt: Path):
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 5000:
        print(f"  cached: {out_pdbqt.name}")
        return
    clean_pdb = WORK / (input_pdb.stem + "_clean.pdb")
    clean_pdb_for_obabel(input_pdb, clean_pdb)
    print(f"  obabel pdb -> pdbqt: {clean_pdb.name} -> {out_pdbqt.name}")
    cmd = [str(OBABEL), "-ipdb", str(clean_pdb), "-opdbqt", "-O", str(out_pdbqt),
           "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0 or not out_pdbqt.exists() or out_pdbqt.stat().st_size < 5000:
        print(r.stdout); print(r.stderr)
        raise RuntimeError(f"obabel failed on {input_pdb}")


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
        return None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def main():
    with open(RESCORE) as f:
        data = json.load(f)
    top = data["top10"][:N_TOP]
    print(f"Counterscreening top {len(top)} ERAP2 hits against ERAP1 + IRAP active sites")

    print("\n=== Receptor prep ===")
    erap1_pdbqt = WORK / "erap1_receptor.pdbqt"
    irap_pdbqt = WORK / "irap_receptor.pdbqt"
    prep_receptor(ERAP1_PDB, erap1_pdbqt)
    prep_receptor(IRAP_PDB, irap_pdbqt)

    print("\n=== Counter-docking ===")
    out_rows = []
    for rank, c in enumerate(top, 1):
        # Find the existing ligand PDBQT from the rescore run by index
        # Strategy: re-locate via SMILES match to all_results
        smi = c["smiles"]
        # Find original index in the full lead-like list
        with open(DATA / "druggpt_p02_candidates.json") as f:
            dg = json.load(f)
        leadlike = [x for x in dg["candidates"] if x["lead_like"]]
        orig_idx = next((i for i, x in enumerate(leadlike) if x["smiles"] == smi), None)
        if orig_idx is None:
            print(f"  [{rank}] no ligand pdbqt found for {smi[:60]}")
            continue
        lig_pdbqt = WORK / f"lig_{orig_idx:03d}.pdbqt"
        if not lig_pdbqt.exists():
            print(f"  [{rank}] missing {lig_pdbqt.name}")
            continue

        e2 = c["vina_dG_kcal_mol"]
        e1_out = WORK / f"cs_erap1_{orig_idx:03d}.pdbqt"
        irap_out = WORK / f"cs_irap_{orig_idx:03d}.pdbqt"
        e1 = run_vina(erap1_pdbqt, lig_pdbqt, ERAP1_CENTER, e1_out)
        irap = run_vina(irap_pdbqt, lig_pdbqt, IRAP_CENTER, irap_out)
        if e1 is None or irap is None:
            print(f"  [{rank}] dock fail e1={e1} irap={irap}")
            continue

        worst_off = max(e1, irap)  # less-negative = weaker binding to off-target
        margin = e2 - worst_off    # negative = ERAP2-selective
        selective = margin <= SELECTIVITY_THRESHOLD
        flag = "**SELECTIVE**" if selective else ""
        print(f"  [{rank:>2}] ERAP2={e2:>6.2f}  ERAP1={e1:>6.2f}  IRAP={irap:>6.2f}  "
              f"margin={margin:>+6.2f}  {flag}  {smi[:55]}")
        out_rows.append({
            **c,
            "dG_erap2_p02": e2,
            "dG_erap1_active": e1,
            "dG_irap_active": irap,
            "selectivity_margin": round(margin, 2),
            "selective": selective,
        })

    out_rows.sort(key=lambda r: r["selectivity_margin"])
    selective = [r for r in out_rows if r["selective"]]

    print()
    print("=" * 80)
    print(f"SELECTIVE HITS (margin <= {SELECTIVITY_THRESHOLD} kcal/mol, ~10x): "
          f"{len(selective)} / {len(out_rows)}")
    print("=" * 80)
    if selective:
        print("Recommended for wet-lab triage:")
        for r in selective:
            print(f"  margin={r['selectivity_margin']:>+6.2f}  "
                  f"ERAP2={r['dG_erap2_p02']:>6.2f}  ERAP1={r['dG_erap1_active']:>6.2f}  "
                  f"IRAP={r['dG_irap_active']:>6.2f}  {r['smiles']}")
    else:
        print("No candidates clear the selectivity gate at this exhaustiveness.")
        print("Top 3 by margin (weak signal — DO NOT order yet):")
        for r in out_rows[:3]:
            print(f"  margin={r['selectivity_margin']:>+6.2f}  {r['smiles']}")

    summary = {
        "n_counterscreened": len(out_rows),
        "n_selective": len(selective),
        "selectivity_threshold_kcal_mol": SELECTIVITY_THRESHOLD,
        "erap1_center": ERAP1_CENTER,
        "irap_center": IRAP_CENTER,
        "selective_hits": selective,
        "all_results_by_margin": out_rows,
    }
    with open(OUT, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
