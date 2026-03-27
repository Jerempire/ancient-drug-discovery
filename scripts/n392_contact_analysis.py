"""
N392 Contact Analysis — Phase 1 of N392 selectivity enhancement.

Identifies which Y87A_Y89A binder residues are within contact distance of
ERAP2 position 392 and nearby divergent residues (Y398, A403, D414).

Uses the Boltz-2 predicted complex CIF file.
Outputs ranked list of binder positions for mutation design.

Usage: python scripts/n392_contact_analysis.py
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import numpy as np
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
CIF_FILE = PROJECT / "data" / "results" / "y87a_cif_files" / "n248_trim_c5_Y87A_Y89A_erap2.cif"
OUTPUT_DIR = PROJECT / "data" / "results" / "n392_selectivity"

OFFSET = 349  # cropped residue number + OFFSET = full ERAP2 numbering

# Key ERAP2 residues to analyze proximity to
TARGET_RESIDUES = {
    392: "K392 (selectivity polymorphism)",
    398: "Y398 (ERAP2-unique Tyr OH)",
    403: "A403 (ERAP2-unique hydrophobic)",
    406: "A406 (ERAP2 vs IRAP-K500)",
    414: "D414 (ERAP2-unique negative charge)",
}

CONTACT_CUTOFF = 6.0  # Angstroms
CLOSE_CUTOFF = 4.5    # close contact

# Y87A_Y89A binder sequence (92aa)
BINDER_SEQ = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_cif_atoms(cif_path):
    """Parse CIF file into list of atom dicts with coordinates."""
    atoms = []
    headers = []
    with open(cif_path) as f:
        in_atom = False
        for line in f:
            line = line.strip()
            if line.startswith('_atom_site.'):
                headers.append(line.split('.')[1])
                in_atom = True
                continue
            if in_atom and not line.startswith('_') and not line.startswith('#') and line:
                parts = line.split()
                if len(parts) >= len(headers):
                    record = dict(zip(headers, parts))
                    record['x'] = float(record['Cartn_x'])
                    record['y'] = float(record['Cartn_y'])
                    record['z'] = float(record['Cartn_z'])
                    record['chain'] = record.get('auth_asym_id', record.get('label_asym_id', '?'))
                    record['resnum'] = int(record.get('auth_seq_id', record.get('label_seq_id', '0')))
                    record['resname'] = record.get('label_comp_id', 'UNK')
                    record['atomname'] = record.get('auth_atom_id', record.get('label_atom_id', '?'))
                    element = record.get('type_symbol', record['atomname'][0])
                    record['element'] = element.upper()
                    atoms.append(record)
            elif in_atom and (line.startswith('#') or line.startswith('loop_')):
                if atoms:
                    break
    return atoms


def get_chain_atoms(atoms, chain_id, skip_h=True):
    """Filter atoms by chain, optionally skipping hydrogen."""
    result = []
    for a in atoms:
        if a['chain'] == chain_id:
            if skip_h and a['element'] == 'H':
                continue
            result.append(a)
    return result


def atoms_to_coords(atom_list):
    """Extract Nx3 coordinate array from atom list."""
    return np.array([[a['x'], a['y'], a['z']] for a in atom_list])


def find_proximal_residues(binder_atoms, target_atoms, target_resnum_full, cutoff):
    """Find binder residues within cutoff of a target residue (full ERAP2 numbering)."""
    # Target residue in cropped numbering
    target_resnum_crop = target_resnum_full - OFFSET

    # Get target residue atoms
    t_atoms = [a for a in target_atoms if a['resnum'] == target_resnum_crop]
    if not t_atoms:
        # Try full numbering directly (some CIFs use full numbering)
        t_atoms = [a for a in target_atoms if a['resnum'] == target_resnum_full]
    if not t_atoms:
        return []

    t_coords = atoms_to_coords(t_atoms)

    # Group binder atoms by residue
    binder_by_res = defaultdict(list)
    for a in binder_atoms:
        binder_by_res[a['resnum']].append(a)

    contacts = []
    for bres_num, b_atoms_list in binder_by_res.items():
        b_coords = atoms_to_coords(b_atoms_list)

        # Compute all pairwise distances
        min_dist = float('inf')
        best_pair = None
        for i, bc in enumerate(b_coords):
            for j, tc in enumerate(t_coords):
                d = np.linalg.norm(bc - tc)
                if d < min_dist:
                    min_dist = d
                    best_pair = (b_atoms_list[i], t_atoms[j])

        if min_dist <= cutoff:
            b_atom, t_atom = best_pair
            binder_aa = AA3TO1.get(b_atom['resname'], '?')
            # 0-indexed position in binder sequence
            binder_pos_0 = bres_num - 1  # CIF residues are 1-indexed for chain B
            if 0 <= binder_pos_0 < len(BINDER_SEQ):
                expected_aa = BINDER_SEQ[binder_pos_0]
            else:
                expected_aa = '?'

            contacts.append({
                "binder_resnum": bres_num,
                "binder_pos_0idx": binder_pos_0,
                "binder_resname": b_atom['resname'],
                "binder_aa": binder_aa,
                "expected_aa": expected_aa,
                "binder_atom": b_atom['atomname'],
                "target_resnum_full": target_resnum_full,
                "target_atom": t_atom['atomname'],
                "target_resname": t_atom['resname'],
                "distance_A": round(min_dist, 2),
                "is_close": bool(min_dist <= CLOSE_CUTOFF),
            })

    return sorted(contacts, key=lambda x: x['distance_A'])


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Parsing CIF: {CIF_FILE}")
    all_atoms = parse_cif_atoms(str(CIF_FILE))
    print(f"  Total atoms: {len(all_atoms)}")

    # Identify chains
    chains = defaultdict(int)
    for a in all_atoms:
        chains[a['chain']] += 1
    print(f"  Chains: {dict(chains)}")

    # Determine which chain is ERAP2 (larger) and which is binder (smaller)
    chain_ids = sorted(chains.keys(), key=lambda c: chains[c], reverse=True)
    target_chain = chain_ids[0]  # larger = ERAP2
    binder_chain = chain_ids[1]  # smaller = binder
    print(f"  Target chain: {target_chain} ({chains[target_chain]} atoms)")
    print(f"  Binder chain: {binder_chain} ({chains[binder_chain]} atoms)")

    target_atoms = get_chain_atoms(all_atoms, target_chain)
    binder_atoms = get_chain_atoms(all_atoms, binder_chain)

    # Check residue numbering scheme
    target_resnums = sorted(set(a['resnum'] for a in target_atoms))
    binder_resnums = sorted(set(a['resnum'] for a in binder_atoms))
    print(f"  Target residue range: {target_resnums[0]}-{target_resnums[-1]} ({len(target_resnums)} residues)")
    print(f"  Binder residue range: {binder_resnums[0]}-{binder_resnums[-1]} ({len(binder_resnums)} residues)")

    # Determine if target uses cropped (1-151) or full (350-500) numbering
    if target_resnums[0] < 100:
        numbering = "cropped"
        effective_offset = OFFSET
        print(f"  Numbering: cropped (res 1 = ERAP2 {1 + OFFSET})")
    else:
        numbering = "full"
        effective_offset = 0
        print(f"  Numbering: full ERAP2")

    # Analyze proximity to each key residue
    all_results = {}
    print(f"\n{'='*70}")
    print("BINDER RESIDUES PROXIMAL TO KEY ERAP2 POSITIONS")
    print(f"Cutoff: {CONTACT_CUTOFF} A | Close: {CLOSE_CUTOFF} A")
    print(f"{'='*70}")

    for target_res, description in TARGET_RESIDUES.items():
        print(f"\n--- {description} ---")

        if numbering == "cropped":
            contacts = find_proximal_residues(binder_atoms, target_atoms, target_res, CONTACT_CUTOFF)
        else:
            # Full numbering — search directly
            t_atoms = [a for a in target_atoms if a['resnum'] == target_res]
            if not t_atoms:
                print(f"  WARNING: residue {target_res} not found in target chain")
                continue
            # Use a modified approach
            contacts = []
            t_coords = atoms_to_coords(t_atoms)
            binder_by_res = defaultdict(list)
            for a in binder_atoms:
                binder_by_res[a['resnum']].append(a)
            for bres_num, b_atoms_list in binder_by_res.items():
                b_coords = atoms_to_coords(b_atoms_list)
                min_dist = float('inf')
                best_b = best_t = None
                for i, bc in enumerate(b_coords):
                    for j, tc in enumerate(t_coords):
                        d = np.linalg.norm(bc - tc)
                        if d < min_dist:
                            min_dist = d
                            best_b = b_atoms_list[i]
                            best_t = t_atoms[j]
                if min_dist <= CONTACT_CUTOFF:
                    binder_aa = AA3TO1.get(best_b['resname'], '?')
                    binder_pos_0 = bres_num - 1
                    expected_aa = BINDER_SEQ[binder_pos_0] if 0 <= binder_pos_0 < len(BINDER_SEQ) else '?'
                    contacts.append({
                        "binder_resnum": bres_num,
                        "binder_pos_0idx": binder_pos_0,
                        "binder_resname": best_b['resname'],
                        "binder_aa": binder_aa,
                        "expected_aa": expected_aa,
                        "binder_atom": best_b['atomname'],
                        "target_resnum_full": target_res,
                        "target_atom": best_t['atomname'],
                        "target_resname": best_t['resname'],
                        "distance_A": round(min_dist, 2),
                        "is_close": bool(min_dist <= CLOSE_CUTOFF),
                    })
            contacts.sort(key=lambda x: x['distance_A'])

        if not contacts:
            print(f"  No binder residues within {CONTACT_CUTOFF} A")
        else:
            print(f"  {'Pos':>4} {'AA':>3} {'Binder Atom':>12} {'Tgt Atom':>10} {'Dist':>6} {'Close':>6}")
            print(f"  {'-'*48}")
            for c in contacts:
                close_mark = " *" if c['is_close'] else ""
                print(f"  {c['binder_pos_0idx']:>4} {c['binder_aa']:>3} {c['binder_atom']:>12} "
                      f"{c['target_atom']:>10} {c['distance_A']:>6.2f}{close_mark}")

        all_results[str(target_res)] = contacts

    # Summary: which binder positions are most relevant for 392
    print(f"\n{'='*70}")
    print("MUTATION CANDIDATES (binder residues near position 392)")
    print(f"{'='*70}")

    res392_contacts = all_results.get("392", [])
    if res392_contacts:
        print(f"\nRanked by distance to residue 392:")
        print(f"{'Pos':>4} {'AA':>3} {'Seq':>4} {'Dist':>6} {'Also near':>30}")
        print(f"{'-'*55}")
        for c in res392_contacts:
            pos = c['binder_pos_0idx']
            # Check what other key residues this binder position is near
            also_near = []
            for tgt, contacts in all_results.items():
                if tgt == "392":
                    continue
                for c2 in contacts:
                    if c2['binder_pos_0idx'] == pos:
                        also_near.append(f"{tgt}({c2['distance_A']:.1f}A)")
                        break
            near_str = ", ".join(also_near) if also_near else "-"
            print(f"{pos:>4} {c['binder_aa']:>3} {BINDER_SEQ[pos] if pos < len(BINDER_SEQ) else '?':>4} "
                  f"{c['distance_A']:>6.2f} {near_str:>30}")

        # Recommend mutation positions
        print(f"\n--- RECOMMENDED MUTATION POSITIONS ---")
        # Prefer positions that are: (a) close to 392, (b) not near other key divergent residues
        # (mutating near Y398/D414 could break ERAP2 selectivity)
        for c in res392_contacts[:5]:
            pos = c['binder_pos_0idx']
            aa = BINDER_SEQ[pos] if pos < len(BINDER_SEQ) else '?'
            print(f"\nPosition {pos} ({aa}, {c['binder_resname']}):")
            print(f"  Distance to K392: {c['distance_A']:.2f} A via {c['binder_atom']}<->{c['target_atom']}")
            print(f"  Current side chain: {aa}")
            if aa in "GAVLI":
                print(f"  Small/hydrophobic -> good candidate for bulky/charged replacement")
            elif aa in "DE":
                print(f"  Negative charge -> likely forms salt bridge with K392 (!)")
                print(f"  Removing this salt bridge would DISFAVOR K392 -> N392 preference")
            elif aa in "KR":
                print(f"  Positive charge -> already repels K392, forms H-bonds with N392")
            elif aa in "NQ":
                print(f"  Amide -> may already H-bond with both K/N392")
            elif aa in "FYW":
                print(f"  Aromatic -> contributes to hydrophobic packing, mutation is risky")
    else:
        print("  No binder residues found near position 392!")
        print("  Check CIF numbering scheme and offset.")

    # Save results
    output = {
        "cif_file": str(CIF_FILE),
        "contact_cutoff_A": CONTACT_CUTOFF,
        "close_cutoff_A": CLOSE_CUTOFF,
        "binder_sequence": BINDER_SEQ,
        "numbering": numbering,
        "contacts_by_target": {k: v for k, v in all_results.items()},
        "mutation_candidates_392": res392_contacts[:5] if res392_contacts else [],
    }
    out_path = OUTPUT_DIR / "n392_contact_analysis.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
