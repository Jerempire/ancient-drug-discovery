"""Zoomed-in view of the P02 pocket with all 4 candidate poses + labeled residues."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
WORK = ROOT / "data/results/smallmol/vina_tight_validation"
FIG_DIR = ROOT / "data/results/smallmol/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

P02_LINING = {426, 568, 667, 668, 671, 701, 702, 705, 706, 709, 712, 940, 941, 944}
P02_HANDLES = {426, 568, 701, 705, 706, 709, 712, 940}
ZN_TRIAD = {370, 374, 393}
P02_CENTER = np.array([-19.05, 9.19, 11.49])

CANDIDATES = [
    (1, "pyridine-azabicyclic (champion)",   "out_01_erap2_p02_s1.pdbqt", "limegreen"),
    (3, "salicylate diacid-amide",           "out_03_erap2_p02_s1.pdbqt", "deepskyblue"),
    (6, "tetrahydronaphtho-pyranone-A",      "out_06_erap2_p02_s1.pdbqt", "magenta"),
    (7, "polycyclic acetamide (orig lead)",  "out_07_erap2_p02_s1.pdbqt", "darkorange"),
]

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def load_ca(pdb, chain="A"):
    p = PDBParser(QUIET=True)
    s = p.get_structure("x", str(pdb))[0]
    resnums, coords = [], []
    for r in s[chain]:
        if r.get_resname() in AA3 and "CA" in r:
            resnums.append(r.get_id()[1])
            coords.append(r["CA"].coord)
    return np.array(resnums), np.array(coords)


def load_pose(pdbqt):
    coords = []; in_m = False; seen = False
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                if seen: break
                seen = True; in_m = True
            elif in_m and line.startswith(("ATOM","HETATM")):
                elem = line[76:78].strip() if len(line) >= 78 else ""
                if elem and elem.startswith("H"): continue
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError: pass
            elif line.startswith("ENDMDL"):
                break
    return np.array(coords) if coords else np.zeros((0,3))


def main():
    resnums, coords = load_ca(ERAP2_PDB)
    res_idx = {r: i for i, r in enumerate(resnums)}

    # Determine zoom bounds from P02 lining residues + ligand poses
    nearby = []
    for r in P02_LINING:
        if r in res_idx:
            nearby.append(coords[res_idx[r]])
    poses = []
    for idx, label, fname, color in CANDIDATES:
        pose = load_pose(WORK / fname)
        poses.append((label, color, pose))
        nearby.extend(pose.tolist())
    nearby = np.array(nearby)
    center = nearby.mean(axis=0)
    # Box: 18 A radius from centroid (full pocket + margin)
    R = 18.0
    xlim = (center[0]-R, center[0]+R)
    ylim = (center[1]-R, center[1]+R)
    zlim = (center[2]-R, center[2]+R)

    fig = plt.figure(figsize=(14, 11))
    ax = fig.add_subplot(111, projection="3d")

    # Backbone trace — only show Ca atoms within the zoom box (plus a small buffer)
    buffer = 6.0
    mask = np.all((coords > [xlim[0]-buffer, ylim[0]-buffer, zlim[0]-buffer]) &
                  (coords < [xlim[1]+buffer, ylim[1]+buffer, zlim[1]+buffer]), axis=1)
    visible_idx = np.where(mask)[0]

    # Draw backbone segments locally
    segs = []
    for i in visible_idx[:-1]:
        if i+1 < len(coords) and resnums[i+1] - resnums[i] <= 2:
            segs.append([coords[i], coords[i+1]])
    if segs:
        lc = Line3DCollection(segs, colors="#888888", linewidths=1.2, alpha=0.5)
        ax.add_collection3d(lc)

    # Backbone Ca dots (visible only) for spatial reference
    vc = coords[visible_idx]
    ax.scatter(vc[:,0], vc[:,1], vc[:,2], s=8, c="#aaaaaa", alpha=0.5, zorder=2)

    # Zn triad if in view
    for r in ZN_TRIAD:
        if r in res_idx:
            c = coords[res_idx[r]]
            if (xlim[0] <= c[0] <= xlim[1] and ylim[0] <= c[1] <= ylim[1]
                and zlim[0] <= c[2] <= zlim[1]):
                ax.scatter(*c, s=500, c="red", edgecolors="black", linewidths=1.5,
                          label="Zn triad" if r == 370 else None, zorder=10, alpha=0.95)
                ax.text(c[0]+1, c[1]+1, c[2]+1, f"H{r}" if r != 393 else f"E{r}",
                       fontsize=11, color="darkred", weight="bold", zorder=15)

    # P02 lining (shared)
    for r in P02_LINING - P02_HANDLES:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.scatter(*c, s=250, c="steelblue", edgecolors="black", linewidths=0.8,
                      label="P02 lining (shared)" if r == 667 else None, zorder=8, alpha=0.85)
            ax.text(c[0]+1, c[1]+1, c[2]+1, str(r), fontsize=9, color="#1a4f80", zorder=15)

    # P02 selectivity handles
    for r in P02_HANDLES:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.scatter(*c, s=500, c="gold", edgecolors="black", linewidths=1.5,
                      label="P02 selectivity handles" if r == 426 else None, zorder=9, alpha=0.95)
            ax.text(c[0]+1, c[1]+1, c[2]+1, str(r), fontsize=11, color="darkgoldenrod",
                   weight="bold", zorder=16)

    # P02 box center
    ax.scatter(*P02_CENTER, s=300, c="black", marker="D",
              label="P02 box center", zorder=11)

    # Ligand poses
    for label, color, pose in poses:
        # Draw atoms as spheres
        ax.scatter(pose[:,0], pose[:,1], pose[:,2], s=80, c=color,
                  edgecolors="black", linewidths=0.5, alpha=0.85,
                  label=label, zorder=12)
        # Draw simple "bonds" between consecutive atoms
        for i in range(len(pose)-1):
            d = np.linalg.norm(pose[i+1] - pose[i])
            if d < 2.5:
                ax.plot([pose[i,0], pose[i+1,0]],
                       [pose[i,1], pose[i+1,1]],
                       [pose[i,2], pose[i+1,2]],
                       color=color, alpha=0.5, linewidth=1.5, zorder=11)

    handles, labels = ax.get_legend_handles_labels()
    seen = set(); uniq = []
    for h, l in zip(handles, labels):
        if l in seen: continue
        seen.add(l); uniq.append((h, l))
    ax.legend([h for h, l in uniq], [l for h, l in uniq],
             loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=10, framealpha=0.95)

    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlim)
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
    ax.set_title("ERAP2 P02 pocket — zoom (gold = ERAP2-unique selectivity handles)\n"
                 f"4 docked candidates engage 5-6 handle residues each", fontsize=12)

    # View 1: looking into the lid from the catalytic chamber side
    ax.view_init(elev=15, azim=-70)
    out1 = FIG_DIR / "erap2_p02_zoom_view1.png"
    plt.tight_layout(); plt.savefig(out1, dpi=160, bbox_inches="tight")
    print(f"Saved: {out1}")

    ax.view_init(elev=15, azim=110)
    out2 = FIG_DIR / "erap2_p02_zoom_view2.png"
    plt.savefig(out2, dpi=160, bbox_inches="tight")
    print(f"Saved: {out2}")

    ax.view_init(elev=70, azim=-90)
    out3 = FIG_DIR / "erap2_p02_zoom_top.png"
    plt.savefig(out3, dpi=160, bbox_inches="tight")
    print(f"Saved: {out3}")
    plt.close()


if __name__ == "__main__":
    main()
