"""
Visualize where the top candidates bind on ERAP2.

Outputs:
  1. data/results/smallmol/figures/erap2_p02_binding.png — 3D matplotlib figure
      - ERAP2 Ca backbone trace (light gray), colored by domain
      - Zn triad residues (red spheres)
      - P02 lining residues (blue spheres)
      - P02 selectivity handles vs both ERAP1+IRAP (gold larger spheres)
      - Champion ligand pose (green spheres)
      - View labeled with anatomical orientation
  2. data/results/smallmol/figures/p02_scene.pdb — minimal PDB with backbone + ligand,
      drop into ChimeraX / VMD / PyMOL / etc.
  3. Clean contact list (standard amino acids only).
"""
from __future__ import annotations

import json
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
OUT_PNG = FIG_DIR / "erap2_p02_binding.png"
OUT_SCENE = FIG_DIR / "p02_scene.pdb"
OUT_CONTACTS = ROOT / "data/results/smallmol/pose_contacts_clean.json"

P02_LINING = {426, 568, 667, 668, 671, 701, 702, 705, 706, 709, 712, 940, 941, 944}
P02_HANDLES = {426, 568, 701, 705, 706, 709, 712, 940}
ZN_TRIAD = {370, 374, 393}
P02_CENTER = np.array([-19.05, 9.19, 11.49])

# ERAP2 domain boundaries (UniProt Q6P179 / ERAP family canonical)
DOMAINS = [
    ("Domain I (N-saddle, β-sheet)",     1,    254, "#cfe2f3"),
    ("Domain II (catalytic, HEXXH+Zn)",   255,  528, "#f4cccc"),
    ("Domain III (β-sheet)",              529,  614, "#d9ead3"),
    ("Domain IV (helical CAP / lid)",     615,  960, "#fff2cc"),
]
DOMAIN_OUTLINE = {
    1: "#1f77b4",
    255: "#d62728",
    529: "#2ca02c",
    615: "#bf8a00",
}

CANDIDATES = [
    (1, "pyridine-azabicyclic (champion)",   "out_01_erap2_p02_s1.pdbqt"),
    (3, "salicylate diacid-amide",           "out_03_erap2_p02_s1.pdbqt"),
    (6, "tetrahydronaphtho-pyranone-A",      "out_06_erap2_p02_s1.pdbqt"),
    (7, "polycyclic acetamide (orig lead)",  "out_07_erap2_p02_s1.pdbqt"),
]

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def load_ca_chain(pdb_path: Path, chain="A"):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))[0]
    resnums, coords = [], []
    for r in s[chain]:
        if r.get_resname() in AA3 and "CA" in r:
            resnums.append(r.get_id()[1])
            coords.append(r["CA"].coord)
    return np.array(resnums), np.array(coords)


def load_pose(pdbqt: Path):
    coords = []
    in_model = False
    seen = False
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                if seen: break
                seen = True; in_model = True
            elif in_model and line.startswith(("ATOM","HETATM")):
                elem = line[76:78].strip() if len(line) >= 78 else ""
                if elem == "H" or (elem and elem.startswith("H")):
                    continue
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    pass
            elif line.startswith("ENDMDL"):
                break
    return np.array(coords) if coords else np.zeros((0, 3))


def domain_for(resnum):
    for name, lo, hi, color in DOMAINS:
        if lo <= resnum <= hi:
            return name, lo, color
    return "?", 0, "#cccccc"


def fig_3d(resnums, coords, ligand_poses, output: Path):
    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(111, projection="3d")

    # Backbone trace colored by domain
    segments = []
    seg_colors = []
    for i in range(len(coords) - 1):
        if resnums[i + 1] - resnums[i] <= 2:
            _, lo, color = domain_for(resnums[i])
            segments.append([coords[i], coords[i + 1]])
            seg_colors.append(DOMAIN_OUTLINE.get(lo, "#888888"))
    lc = Line3DCollection(segments, colors=seg_colors, linewidths=0.8, alpha=0.4)
    ax.add_collection3d(lc)

    # Zn triad — red spheres
    res_idx = {r: i for i, r in enumerate(resnums)}
    for r in ZN_TRIAD:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.scatter(*c, s=400, c="red", edgecolors="black", linewidths=1.5,
                      label="Zn triad" if r == 370 else None, zorder=10, alpha=0.9)

    # P02 lining (not handles) — light blue
    for r in P02_LINING - P02_HANDLES:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.scatter(*c, s=180, c="steelblue", edgecolors="black", linewidths=0.5,
                      label="P02 lining (shared)" if r == 667 else None, zorder=8, alpha=0.85)

    # P02 selectivity handles vs BOTH ERAP1+IRAP — gold/yellow
    for r in P02_HANDLES:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.scatter(*c, s=320, c="gold", edgecolors="black", linewidths=1.2,
                      label="P02 selectivity handles" if r == 426 else None, zorder=9, alpha=0.95)

    # P02 box center — black diamond
    ax.scatter(*P02_CENTER, s=200, c="black", marker="D",
              label="P02 box center", zorder=11)

    # Ligand poses
    ligand_colors = ["limegreen", "deepskyblue", "magenta", "darkorange"]
    for i, (label, pose) in enumerate(ligand_poses):
        color = ligand_colors[i % len(ligand_colors)]
        ax.scatter(pose[:, 0], pose[:, 1], pose[:, 2], s=40, c=color,
                  edgecolors="black", linewidths=0.3, alpha=0.7,
                  label=label, zorder=12)

    # Annotate selectivity handles with residue numbers
    for r in P02_HANDLES:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.text(c[0] + 1.5, c[1] + 1.5, c[2] + 1.5, str(r), fontsize=9,
                   color="darkgoldenrod", weight="bold", zorder=13)
    # Annotate Zn triad too
    for r in ZN_TRIAD:
        if r in res_idx:
            c = coords[res_idx[r]]
            ax.text(c[0] + 1.5, c[1] - 2, c[2] - 2, str(r), fontsize=9,
                   color="darkred", weight="bold", zorder=13)

    # Domain legend
    domain_handles = [
        Line2D([0], [0], color=DOMAIN_OUTLINE[1], linewidth=3,
               label="Domain I (β-saddle, res 1-254)"),
        Line2D([0], [0], color=DOMAIN_OUTLINE[255], linewidth=3,
               label="Domain II (CATALYTIC, res 255-528)"),
        Line2D([0], [0], color=DOMAIN_OUTLINE[529], linewidth=3,
               label="Domain III (β-sheet, res 529-614)"),
        Line2D([0], [0], color=DOMAIN_OUTLINE[615], linewidth=3,
               label="Domain IV (HELICAL CAP / LID, res 615-960)"),
    ]

    handles, labels = ax.get_legend_handles_labels()
    seen = set(); uniq = []
    for h, l in zip(handles, labels):
        if l in seen: continue
        seen.add(l); uniq.append((h, l))
    ax.legend(handles=domain_handles + [h for h, l in uniq],
             labels=[h.get_label() for h in domain_handles] + [l for h, l in uniq],
             loc="upper left", bbox_to_anchor=(1.05, 1), fontsize=9, framealpha=0.95)

    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
    ax.set_title(
        "ERAP2 with P02 pocket binding (4 candidates)\n"
        f"Zn triad (red) is ~{np.linalg.norm(P02_CENTER - coords[res_idx[370]]):.1f} Å from P02 box center (black ◆)",
        fontsize=11)
    # 30° elevation, view aimed at P02
    ax.view_init(elev=20, azim=-60)
    plt.tight_layout()
    plt.savefig(output, dpi=140, bbox_inches="tight")
    print(f"Saved figure: {output}")

    # Second angle
    ax.view_init(elev=20, azim=120)
    output2 = output.with_stem(output.stem + "_view2")
    plt.savefig(output2, dpi=140, bbox_inches="tight")
    print(f"Saved figure: {output2}")
    plt.close()


def write_scene_pdb(resnums, coords, ligand_poses, output: Path,
                   lining_set, handles_set, zn_set):
    """Write a minimal PDB with full Ca trace + ligand poses, residue numbers preserved."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(ERAP2_PDB))[0]
    atom_serial = 1
    with open(output, "w") as f:
        f.write("HEADER    ERAP2 P02 binding scene (Ca backbone + lining + ligand poses)\n")
        f.write("REMARK 100 Zn triad: H370 H374 E393 (red)\n")
        f.write("REMARK 100 P02 lining: 426 568 667 668 671 701 702 705 706 709 712 940 941 944\n")
        f.write("REMARK 100 Selectivity handles (vs ERAP1+IRAP): 426 568 701 705 706 709 712 940\n")
        # Write all heavy atoms of P02 lining + Zn triad with original residue info
        important = lining_set | zn_set
        for r in s["A"]:
            if r.get_resname() not in AA3:
                continue
            rn = r.get_id()[1]
            if rn in important:
                for a in r:
                    if a.element == "H": continue
                    f.write(f"ATOM  {atom_serial:>5} {a.get_name():<4} {r.get_resname():<3} A"
                            f"{rn:>4}    "
                            f"{a.coord[0]:>8.3f}{a.coord[1]:>8.3f}{a.coord[2]:>8.3f}"
                            f"  1.00 50.00           {a.element:>2}\n")
                    atom_serial += 1
        # Ca trace for all other protein residues
        for r in s["A"]:
            if r.get_resname() not in AA3 or "CA" not in r:
                continue
            rn = r.get_id()[1]
            if rn in important:
                continue
            ca = r["CA"]
            f.write(f"ATOM  {atom_serial:>5}  CA  {r.get_resname():<3} A"
                    f"{rn:>4}    "
                    f"{ca.coord[0]:>8.3f}{ca.coord[1]:>8.3f}{ca.coord[2]:>8.3f}"
                    f"  1.00 30.00            C\n")
            atom_serial += 1
        f.write("TER\n")
        # Ligand poses as separate chains
        chain_ids = "LMNO"
        for i, (label, pose) in enumerate(ligand_poses):
            ch = chain_ids[i % len(chain_ids)]
            for j, c in enumerate(pose):
                f.write(f"HETATM{atom_serial:>5}  C{j+1:<2} LIG {ch}{i+1:>4}    "
                        f"{c[0]:>8.3f}{c[1]:>8.3f}{c[2]:>8.3f}"
                        f"  1.00 50.00            C\n")
                atom_serial += 1
            f.write("TER\n")
        f.write("END\n")
    print(f"Saved scene PDB: {output}")


def clean_contacts(ligand_coords, protein_pdb, cutoff=4.0):
    """Recompute contacts filtering to STANDARD AMINO ACIDS ONLY."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(protein_pdb))[0]
    aa_atoms = []
    for chain in s:
        for r in chain:
            if r.get_resname() in AA3:
                for a in r:
                    if a.element != "H":
                        aa_atoms.append((a, r))
    if not aa_atoms or len(ligand_coords) == 0:
        return []
    coords = np.array([a.coord for a, _ in aa_atoms])
    contact_res = set()
    for lc in ligand_coords:
        d = np.linalg.norm(coords - lc, axis=1)
        for i, dist in enumerate(d):
            if dist <= cutoff:
                r = aa_atoms[i][1]
                contact_res.add((r.get_id()[1], r.get_resname(), r.get_parent().id))
    return sorted(contact_res)


def main():
    print("Loading ERAP2 backbone ...")
    resnums, coords = load_ca_chain(ERAP2_PDB)
    print(f"  {len(resnums)} Ca atoms")

    ligand_poses = []
    for idx, label, fname in CANDIDATES:
        pose = load_pose(WORK / fname)
        print(f"  loaded {label}: {len(pose)} heavy atoms")
        ligand_poses.append((label, pose))

    # Recompute clean contacts for all top hits, all 3 targets
    print("\n--- Clean contact analysis (amino acids only) ---")
    targets = [
        ("ERAP2 P02", ERAP2_PDB, P02_CENTER, P02_LINING, ZN_TRIAD),
        ("ERAP1 Zn",  ROOT / "data/structures/erap1_2YD0_experimental.pdb",
         np.array([6.512, 71.495, 21.017]), set(), {353, 357, 376}),
        ("IRAP Zn",   ROOT / "data/structures/irap_5MJ6_experimental.pdb",
         np.array([9.347, -31.132, -22.851]), set(), {464, 468, 487}),
    ]
    out_payload = []
    pose_files = {
        "ERAP2 P02": "erap2_p02",
        "ERAP1 Zn":  "erap1_zn",
        "IRAP Zn":   "irap_zn",
    }
    for idx, label, _ in CANDIDATES:
        print(f"\n[{idx}] {label}")
        row = {"idx": idx, "label": label}
        for tname, tpdb, tcenter, lining, zn_res in targets:
            pose = load_pose(WORK / f"out_{idx:02d}_{pose_files[tname]}_s1.pdbqt")
            contacts = clean_contacts(pose, tpdb)
            nums = [r[0] for r in contacts]
            centroid = pose.mean(axis=0) if len(pose) else np.zeros(3)
            offset = float(np.linalg.norm(centroid - tcenter))
            zn_overlap = [n for n in nums if n in zn_res]
            print(f"  {tname:>10}: offset={offset:.1f} Å, "
                  f"{len(contacts)} AA contacts, Zn-triad overlap: {zn_overlap}")
            print(f"             contacts: {nums}")
            if lining:
                lo = [n for n in nums if n in lining]
                ho = [n for n in nums if n in P02_HANDLES]
                print(f"             P02 lining overlap: {lo} ({100*len(lo)/max(len(contacts),1):.0f}%)")
                print(f"             selectivity-handle overlap: {ho}")
            row[tname] = {
                "offset_from_box_A": round(offset, 2),
                "n_contacts": len(contacts),
                "contacts": [{"resnum": n, "resname": rn, "chain": ch} for n, rn, ch in contacts],
                "zn_triad_overlap": zn_overlap,
                "p02_lining_overlap": [n for n in nums if n in lining] if lining else [],
                "p02_handle_overlap": [n for n in nums if n in P02_HANDLES] if lining else [],
            }
        out_payload.append(row)
    with open(OUT_CONTACTS, "w") as f:
        json.dump(out_payload, f, indent=2)
    print(f"\nSaved clean contacts: {OUT_CONTACTS}")

    print("\n--- Generating 3D figure ---")
    fig_3d(resnums, coords, ligand_poses, OUT_PNG)
    write_scene_pdb(resnums, coords, ligand_poses, OUT_SCENE,
                    P02_LINING, P02_HANDLES, ZN_TRIAD)


if __name__ == "__main__":
    main()
