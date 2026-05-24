"""
For each top candidate, parse the Vina output pose and report which protein
residues the ligand is actually contacting (any protein heavy atom within 4 A
of any ligand heavy atom).

Compares: declared box center, actual ligand centroid, contacted residues.
Cross-references against P02 lining residues / Zn-binding residues.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "data/results/smallmol/vina_tight_validation"
ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = ROOT / "data/structures/erap1_2YD0_experimental.pdb"
IRAP_PDB = ROOT / "data/structures/irap_5MJ6_experimental.pdb"
OUT = ROOT / "data/results/smallmol/pose_contacts.json"

P02_LINING = {426, 568, 667, 668, 671, 701, 702, 705, 706, 709, 712, 940, 941, 944}
P02_HANDLES_BOTH = {426, 568, 701, 705, 706, 709, 712, 940}
P02_CENTROID = np.array([-19.05, 9.19, 11.49])
ERAP1_ZN = np.array([6.512, 71.495, 21.017])
IRAP_ZN = np.array([9.347, -31.132, -22.851])

# HEXXH zinc-binding residues (catalytic triad)
ERAP2_ZN_RES = {370, 374, 393}
ERAP1_ZN_RES = {353, 357, 376}
IRAP_ZN_RES = {464, 468, 487}

CANDIDATES = [
    (1, "pyridine-azabicyclic (champion)"),
    (2, "decalone-pyrrolidinone"),
    (3, "salicylate diacid-amide"),
    (6, "tetrahydronaphtho-pyranone-A"),
    (7, "polycyclic acetamide (orig lead)"),
]


def parse_pdbqt_ligand(pdbqt: Path) -> np.ndarray:
    """Extract heavy-atom coords of the first MODEL (best pose) from a Vina output PDBQT."""
    coords = []
    in_model_1 = False
    seen_model = False
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                if seen_model:
                    break  # only want MODEL 1
                seen_model = True
                in_model_1 = True
            elif in_model_1 and line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip() if len(line) >= 78 else line[12:14].strip()
                if elem == "H" or (elem and elem[0] == "H"):
                    continue
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
            elif line.startswith("ENDMDL"):
                in_model_1 = False
                break
    return np.array(coords) if coords else np.zeros((0, 3))


def get_contacts(protein_pdb: Path, ligand_coords: np.ndarray, cutoff=4.0, chain="A"):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(protein_pdb))[0]
    atoms = [a for a in s[chain].get_atoms() if a.element != "H"]
    if len(ligand_coords) == 0 or len(atoms) == 0:
        return [], None
    ns = NeighborSearch(atoms)
    contact_res = set()
    for lc in ligand_coords:
        for a in ns.search(lc, cutoff, level="A"):
            r = a.get_parent()
            try:
                resnum = r.get_id()[1]
                resname = r.get_resname()
                contact_res.add((resnum, resname))
            except Exception:
                pass
    centroid = ligand_coords.mean(axis=0)
    return sorted(contact_res), centroid


def fmt_res_set(contacts, highlight_set, label):
    nums = [r[0] for r in contacts]
    overlap = [n for n in nums if n in highlight_set]
    return f"  {label}: {len(contacts)} residues   overlap with target set: {overlap}"


def main():
    print("=" * 80)
    print("POSE CONTACT ANALYSIS — where do the candidates actually bind?")
    print("=" * 80)
    print(f"\nP02 centroid (ERAP2 target): {P02_CENTROID.tolist()}")
    print(f"P02 lining residues: {sorted(P02_LINING)}")
    print(f"P02 selectivity handles (vs both ERAP1+IRAP): {sorted(P02_HANDLES_BOTH)}")
    print(f"\nERAP1 catalytic Zn coord: {ERAP1_ZN.tolist()}   triad: {sorted(ERAP1_ZN_RES)}")
    print(f"IRAP catalytic Zn coord:  {IRAP_ZN.tolist()}   triad: {sorted(IRAP_ZN_RES)}")

    results = []
    for idx, label in CANDIDATES:
        print(f"\n{'=' * 80}")
        print(f"[{idx}] {label}")
        print(f"{'=' * 80}")
        row = {"idx": idx, "label": label}
        for target_name, prot_pdb, box_center, lining_set, zn_res in [
            ("ERAP2 P02", ERAP2_PDB, P02_CENTROID, P02_LINING, ERAP2_ZN_RES),
            ("ERAP1 Zn", ERAP1_PDB, ERAP1_ZN, set(), ERAP1_ZN_RES),
            ("IRAP Zn", IRAP_PDB, IRAP_ZN, set(), IRAP_ZN_RES),
        ]:
            lig_pdbqt = WORK / f"out_{idx:02d}_{target_name.lower().replace(' ', '_').replace('p02', 'p02').replace('zn', 'zn')}_s1.pdbqt"
            # Fix the filename mapping
            if "ERAP2" in target_name:
                lig_pdbqt = WORK / f"out_{idx:02d}_erap2_p02_s1.pdbqt"
            elif "ERAP1" in target_name:
                lig_pdbqt = WORK / f"out_{idx:02d}_erap1_zn_s1.pdbqt"
            else:
                lig_pdbqt = WORK / f"out_{idx:02d}_irap_zn_s1.pdbqt"

            if not lig_pdbqt.exists():
                print(f"\n  {target_name}: pose file missing: {lig_pdbqt.name}")
                continue

            lig_coords = parse_pdbqt_ligand(lig_pdbqt)
            if len(lig_coords) == 0:
                print(f"\n  {target_name}: empty pose")
                continue

            contacts, centroid = get_contacts(prot_pdb, lig_coords)
            offset = np.linalg.norm(centroid - box_center)
            nums = [r[0] for r in contacts]
            zn_overlap = [n for n in nums if n in zn_res]
            print(f"\n  {target_name}")
            print(f"    ligand centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})")
            print(f"    distance from box center: {offset:.2f} A")
            print(f"    contacts ({len(contacts)} residues): {nums}")
            if zn_res:
                print(f"    overlap with Zn triad ({sorted(zn_res)}): {zn_overlap}  "
                      f"{'*** binds catalytic Zn site ***' if len(zn_overlap) >= 2 else ('partial' if zn_overlap else 'NOT at Zn site')}")
            if lining_set:
                lining_overlap = [n for n in nums if n in lining_set]
                handle_overlap = [n for n in nums if n in P02_HANDLES_BOTH]
                print(f"    overlap with P02 lining ({len(lining_set)} res): {lining_overlap}  "
                      f"({100*len(lining_overlap)/max(len(contacts),1):.0f}% of contacts)")
                print(f"    overlap with P02 selectivity handles: {handle_overlap}  "
                      f"{'*** engages selectivity handles ***' if len(handle_overlap) >= 2 else 'weak engagement'}")
            row[target_name] = {
                "ligand_centroid": [round(float(x), 2) for x in centroid],
                "offset_from_box": round(float(offset), 2),
                "contacts": [{"resnum": n, "resname": rn} for n, rn in contacts],
                "zn_triad_overlap": zn_overlap,
                "lining_overlap": [n for n in nums if n in lining_set] if lining_set else [],
                "handle_overlap": [n for n in nums if n in P02_HANDLES_BOTH] if lining_set else [],
            }
        results.append(row)

    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
