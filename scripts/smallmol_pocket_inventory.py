"""
Small-molecule pocket inventory for ERAP2.

Finds druggable pockets on AF-ERAP2, scores each by:
  - geometric druggability (LIGSITE-style grid burial)
  - hydrophobicity (Kyte-Doolittle of lining residues)
  - divergence vs ERAP1 (selectivity handle for ERAP2-selective inhibitors)
  - proximity to catalytic Zn triad (active site H370/H374/E393)

Outputs ranked pockets, flags top-1 allosteric (non-active-site) target.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.Align import PairwiseAligner
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
ERAP2 = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1 = ROOT / "data/structures/erap1_wt_alphafold.pdb"
IRAP = ROOT / "data/structures/irap_5MJ6_experimental.pdb"
OUT = ROOT / "data/results/smallmol/pocket_inventory.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

ZN_TRIAD = (370, 374, 393)
KD = {  # Kyte-Doolittle hydrophobicity
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_residues(pdb_path: Path, chain_id: str | None = None):
    parser = PDB.PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))
    model = s[0]
    # If chain not specified, use the longest protein chain
    if chain_id is None:
        best, best_len = None, 0
        for ch in model:
            n = sum(1 for r in ch if r.get_resname() in AA3TO1 and "CA" in r)
            if n > best_len:
                best, best_len = ch.id, n
        chain_id = best
    residues = []
    for r in model[chain_id]:
        if r.get_resname() in AA3TO1 and "CA" in r:
            residues.append(r)
    return residues


def residue_seq(residues):
    return "".join(AA3TO1[r.get_resname()] for r in residues)


def divergent_resnums(e2_res, other_res):
    s2, so = residue_seq(e2_res), residue_seq(other_res)
    aln = PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -2
    aln.extend_gap_score = -0.5
    a = aln.align(s2, so)[0]
    s2a, soa = str(a[0]), str(a[1])
    divergent = set()
    i2 = 0
    for c2, co in zip(s2a, soa):
        if c2 != "-":
            if co == "-" or c2 != co:
                divergent.add(e2_res[i2].get_id()[1])
            i2 += 1
    return divergent


def heavy_atom_coords(residues):
    coords, owner = [], []
    for r in residues:
        for a in r:
            if a.element != "H":
                coords.append(a.coord)
                owner.append(r.get_id()[1])
    return np.array(coords), np.array(owner)


def ligsite_grid(atom_coords, spacing=1.5, probe=1.4, vdw=1.8):
    """Build a grid; return interior grid points buried in pockets (PSP count >=4)."""
    mn = atom_coords.min(axis=0) - 5
    mx = atom_coords.max(axis=0) + 5
    xs = np.arange(mn[0], mx[0], spacing)
    ys = np.arange(mn[1], mx[1], spacing)
    zs = np.arange(mn[2], mx[2], spacing)
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    tree = cKDTree(atom_coords)
    nearest, _ = tree.query(grid, k=1)
    # Solvent-accessible: not inside protein
    solvent = (nearest > (vdw + probe)) & (nearest < 8.0)
    grid_s = grid[solvent]

    # 7 scan directions (3 axes + 4 cube diagonals)
    dirs = np.array([
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1],
    ], dtype=float)
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    step = spacing
    max_reach = 12.0
    n_steps = int(max_reach / step)

    psp = np.zeros(len(grid_s), dtype=int)
    for d in dirs:
        # forward
        hit_fwd = np.zeros(len(grid_s), dtype=bool)
        for k in range(1, n_steps + 1):
            probe_pts = grid_s + d * (k * step)
            dist, _ = tree.query(probe_pts, k=1)
            hit_fwd |= dist < vdw
        # backward
        hit_bwd = np.zeros(len(grid_s), dtype=bool)
        for k in range(1, n_steps + 1):
            probe_pts = grid_s - d * (k * step)
            dist, _ = tree.query(probe_pts, k=1)
            hit_bwd |= dist < vdw
        psp += (hit_fwd & hit_bwd).astype(int)

    pocket_pts = grid_s[psp >= 5]
    return pocket_pts


def cluster_pockets(pocket_pts, eps=2.5, min_samples=8):
    if len(pocket_pts) == 0:
        return []
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(pocket_pts)
    clusters = []
    for cid in sorted(set(db.labels_) - {-1}):
        clusters.append(pocket_pts[db.labels_ == cid])
    clusters.sort(key=len, reverse=True)
    return clusters


def lining_residues(cluster_pts, atom_coords, owner, cutoff=5.0):
    tree = cKDTree(atom_coords)
    near = tree.query_ball_point(cluster_pts, r=cutoff)
    res_set = set()
    for idx_list in near:
        for i in idx_list:
            res_set.add(int(owner[i]))
    return sorted(res_set)


def score_pocket(cluster_pts, lining, div_e1, div_irap, resname_by_num, triad_centroid):
    n_res = len(lining)
    n_div_e1 = len([r for r in lining if r in div_e1])
    n_div_irap = len([r for r in lining if r in div_irap])
    n_div_both = len([r for r in lining if r in div_e1 and r in div_irap])
    df_e1 = n_div_e1 / n_res if n_res else 0.0
    df_irap = n_div_irap / n_res if n_res else 0.0
    df_both = n_div_both / n_res if n_res else 0.0
    hphob = float(np.mean([KD.get(resname_by_num.get(r, "ALA"), 0.0) for r in lining])) if lining else 0.0
    centroid = cluster_pts.mean(axis=0)
    dist_active = float(np.linalg.norm(centroid - triad_centroid))
    is_active = dist_active < 10.0
    volume = float(len(cluster_pts) * (1.5 ** 3))
    # Druggability composite — selectivity now requires divergence vs BOTH ERAP1 and IRAP
    size_score = min(volume / 300.0, 2.0)
    drug_score = size_score + 0.5 * hphob / 4.5 + 1.5 * df_both + 0.5 * (df_e1 + df_irap) / 2
    return {
        "n_lining_residues": n_res,
        "n_divergent_vs_erap1": n_div_e1,
        "n_divergent_vs_irap": n_div_irap,
        "n_divergent_vs_both": n_div_both,
        "div_frac_erap1": round(df_e1, 3),
        "div_frac_irap": round(df_irap, 3),
        "div_frac_both": round(df_both, 3),
        "mean_hydrophobicity_kd": round(hphob, 2),
        "volume_A3_est": round(volume, 1),
        "centroid_xyz": [round(float(x), 2) for x in centroid],
        "dist_to_zn_triad_A": round(dist_active, 2),
        "is_active_site": bool(is_active),
        "druggability_score": round(drug_score, 3),
    }


def main():
    print("=" * 70)
    print("ERAP2 SMALL-MOLECULE POCKET INVENTORY")
    print("=" * 70)
    print(f"Structure: {ERAP2}")

    e2_res = load_residues(ERAP2)
    e1_res = load_residues(ERAP1)
    irap_res = load_residues(IRAP)
    print(f"ERAP2: {len(e2_res)} residues | ERAP1: {len(e1_res)} residues | IRAP: {len(irap_res)} residues")

    div_e1 = divergent_resnums(e2_res, e1_res)
    div_irap = divergent_resnums(e2_res, irap_res)
    div_both = div_e1 & div_irap
    print(f"Divergent vs ERAP1: {len(div_e1)}  |  vs IRAP: {len(div_irap)}  |  vs BOTH: {len(div_both)}")

    resname_by_num = {r.get_id()[1]: r.get_resname() for r in e2_res}
    atoms, owner = heavy_atom_coords(e2_res)
    print(f"Heavy atoms: {len(atoms)}")

    triad_coords = []
    for rn in ZN_TRIAD:
        for r in e2_res:
            if r.get_id()[1] == rn and "CA" in r:
                triad_coords.append(r["CA"].coord)
    triad_centroid = np.mean(triad_coords, axis=0)
    print(f"Zn triad centroid: {[round(float(x),2) for x in triad_centroid]}")

    print("\nScanning grid (LIGSITE-style PSP, this takes ~1-2 min)...")
    pts = ligsite_grid(atoms, spacing=1.5)
    print(f"Pocket grid points: {len(pts)}")

    clusters = cluster_pockets(pts, eps=2.5, min_samples=8)
    print(f"Pocket clusters: {len(clusters)}")

    pockets = []
    for i, cl in enumerate(clusters[:20]):
        lining = lining_residues(cl, atoms, owner)
        if len(lining) < 6:
            continue
        s = score_pocket(cl, lining, div_e1, div_irap, resname_by_num, triad_centroid)
        s["pocket_id"] = f"P{i+1:02d}"
        s["lining_residues"] = lining
        s["divergent_vs_erap1"] = [r for r in lining if r in div_e1]
        s["divergent_vs_irap"] = [r for r in lining if r in div_irap]
        s["divergent_vs_both"] = [r for r in lining if r in div_e1 and r in div_irap]
        pockets.append(s)

    pockets.sort(key=lambda p: p["druggability_score"], reverse=True)

    print()
    print("=" * 110)
    print(f"{'ID':<5}{'Vol':>6}{'Res':>4}{'/E1':>6}{'/IRAP':>7}{'/BOTH':>7}{'KD':>6}{'d(Zn)':>7}{'Site':>11}{'DrugScore':>11}")
    print("-" * 110)
    for p in pockets[:15]:
        site = "ACTIVE" if p["is_active_site"] else "allosteric"
        print(f"{p['pocket_id']:<5}{p['volume_A3_est']:>6.0f}{p['n_lining_residues']:>4}"
              f"{p['div_frac_erap1']*100:>5.0f}%{p['div_frac_irap']*100:>6.0f}%{p['div_frac_both']*100:>6.0f}%"
              f"{p['mean_hydrophobicity_kd']:>6.2f}{p['dist_to_zn_triad_A']:>6.1f}A{site:>11}{p['druggability_score']:>11.2f}")

    # Top allosteric pocket
    allosteric = [p for p in pockets if not p["is_active_site"]]
    top_allo = allosteric[0] if allosteric else None

    print()
    print("=" * 70)
    if top_allo:
        print(f"TOP ALLOSTERIC POCKET: {top_allo['pocket_id']}")
        print(f"  Volume: {top_allo['volume_A3_est']:.0f} A^3   "
              f"Lining: {top_allo['n_lining_residues']} res")
        print(f"  Selectivity vs ERAP1: {top_allo['n_divergent_vs_erap1']}/{top_allo['n_lining_residues']} "
              f"({top_allo['div_frac_erap1']*100:.0f}%)")
        print(f"  Selectivity vs IRAP:  {top_allo['n_divergent_vs_irap']}/{top_allo['n_lining_residues']} "
              f"({top_allo['div_frac_irap']*100:.0f}%)")
        print(f"  Divergent vs BOTH:    {top_allo['n_divergent_vs_both']}/{top_allo['n_lining_residues']} "
              f"({top_allo['div_frac_both']*100:.0f}%)  <-- true ERAP2-unique handles")
        print(f"  Distance to Zn triad: {top_allo['dist_to_zn_triad_A']:.1f} A")
        print(f"  Druggability score: {top_allo['druggability_score']:.2f}")
        print(f"  Centroid (for docking grid): {top_allo['centroid_xyz']}")
        print(f"  Divergent-vs-BOTH residues (use these for selectivity SAR):")
        print(f"    {top_allo['divergent_vs_both']}")
    else:
        print("No allosteric pocket found above thresholds.")
    print("=" * 70)

    out = {
        "structure": str(ERAP2),
        "n_divergent_vs_erap1_total": len(div_e1),
        "n_divergent_vs_irap_total": len(div_irap),
        "n_divergent_vs_both_total": len(div_both),
        "zn_triad_centroid": [round(float(x), 2) for x in triad_centroid],
        "pockets": pockets,
        "top_allosteric": top_allo,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
