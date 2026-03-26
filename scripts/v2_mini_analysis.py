"""V2-Mini Analysis: Find the minimal SPPS-synthesizable fragment of the V2 binder."""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import numpy as np
from scipy.spatial.distance import cdist

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(PROJECT, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
OUTPUT = os.path.join(PROJECT, "data/results/v43_validation/terminal_e/v2_mini_analysis.json")
OFFSET = 349

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def parse_cif(path, chain):
    atoms = []
    headers = []
    with open(path) as f:
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
                    r = dict(zip(headers, parts))
                    ch = r.get('auth_asym_id', r.get('label_asym_id'))
                    if ch == chain:
                        atoms.append({
                            'resnum': int(r.get('auth_seq_id', '0')),
                            'resname': r.get('label_comp_id', 'UNK'),
                            'atomname': r.get('auth_atom_id', r.get('label_atom_id', '?')),
                            'x': float(r['Cartn_x']),
                            'y': float(r['Cartn_y']),
                            'z': float(r['Cartn_z']),
                        })
            elif in_atom and (line.startswith('#') or line.startswith('loop_')):
                if atoms:
                    break
    return atoms


def main():
    erap2 = parse_cif(CIF, 'A')
    binder = parse_cif(CIF, 'B')

    coords_e = np.array([[a['x'], a['y'], a['z']] for a in erap2])
    coords_b = np.array([[a['x'], a['y'], a['z']] for a in binder])
    dists = cdist(coords_b, coords_e)

    binder_residues = sorted(set(a['resnum'] for a in binder))
    dp_set = set(range(353, 368)) | set(range(400, 415))

    print("=" * 80)
    print("V2-MINI ANALYSIS: Which binder residues drive the interface?")
    print("=" * 80)

    # Per-residue interface analysis
    residue_data = {}
    print(f"\n{'Res#':>5} {'AA':>3} {'MinDist':>8} {'Cnt<5A':>7} {'Cnt<4A':>7} {'Role':>10}")
    print("-" * 50)

    for resnum in binder_residues:
        b_idx = [i for i, a in enumerate(binder) if a['resnum'] == resnum]
        resname = binder[b_idx[0]]['resname']
        aa = AA3TO1.get(resname, 'X')

        min_dist = float(dists[b_idx].min())
        contacts_5 = int((dists[b_idx] < 5.0).any(axis=0).sum())
        contacts_4 = int((dists[b_idx] < 4.0).any(axis=0).sum())

        erap2_contacts = set()
        for bi in b_idx:
            for ei in np.where(dists[bi] < 5.0)[0]:
                erap2_contacts.add(erap2[ei]['resnum'] + OFFSET)

        k392 = 392 in erap2_contacts
        dp = bool(erap2_contacts & dp_set)
        role = 'K392!' if k392 else 'DP' if dp else 'iface' if contacts_5 > 0 else ''

        residue_data[resnum] = {
            'aa': aa, 'min_dist': round(min_dist, 2),
            'contacts_5': contacts_5, 'contacts_4': contacts_4,
            'erap2_contacts': sorted(erap2_contacts), 'role': role,
        }

        if contacts_5 > 0:
            sig = '***' if contacts_4 > 3 else '**' if contacts_4 > 0 else '*'
            print(f"{resnum:>5} {aa:>3} {min_dist:>7.2f}A {contacts_5:>7} {contacts_4:>7} {role:>10} {sig}")

    # Identify contiguous interface segments
    interface_res = [r for r in binder_residues if residue_data[r]['contacts_5'] > 0]
    total_contacts_5 = sum(residue_data[r]['contacts_5'] for r in binder_residues)
    total_contacts_4 = sum(residue_data[r]['contacts_4'] for r in binder_residues)

    print(f"\nInterface: {len(interface_res)}/{len(binder_residues)} residues")
    print(f"Total contacts (<5A): {total_contacts_5}, (<4A): {total_contacts_4}")

    segments = []
    current = [interface_res[0]]
    for r in interface_res[1:]:
        if r <= current[-1] + 2:
            current.append(r)
        else:
            segments.append(current)
            current = [r]
    segments.append(current)

    print(f"\nContiguous interface segments:")
    for seg in segments:
        seq = ''.join(residue_data.get(r, {'aa': '?'})['aa'] for r in range(seg[0], seg[-1] + 1))
        erap2_touched = set()
        for r in seg:
            if r in residue_data:
                erap2_touched.update(residue_data[r]['erap2_contacts'])
        dp_n = len(erap2_touched & dp_set)
        k = 392 in erap2_touched
        print(f"  {seg[0]}-{seg[-1]} ({seg[-1]-seg[0]+1}aa): {seq}")
        print(f"    ERAP2 contacts: {len(erap2_touched)}, DP: {dp_n}, K392: {k}")

    # Sliding window analysis
    print("\n" + "=" * 80)
    print("SLIDING WINDOW: Best contiguous fragment by size")
    print("=" * 80)

    results = []
    for ws in [20, 25, 30, 35, 40, 45, 50]:
        best = {'score': 0}
        for start in range(1, 93 - ws + 1):
            end = start + ws - 1
            c5 = sum(residue_data.get(r, {'contacts_5': 0})['contacts_5'] for r in range(start, end + 1))
            c4 = sum(residue_data.get(r, {'contacts_4': 0})['contacts_4'] for r in range(start, end + 1))
            if c5 > best['score']:
                best = {'score': c5, 'c4': c4, 'start': start, 'end': end}

        cov = best['score'] / total_contacts_5 * 100 if total_contacts_5 else 0
        seq = ''.join(residue_data.get(r, {'aa': '?'})['aa'] for r in range(best['start'], best['end'] + 1))

        erap2_w = set()
        for r in range(best['start'], best['end'] + 1):
            if r in residue_data:
                erap2_w.update(residue_data[r]['erap2_contacts'])
        dp_n = len(erap2_w & dp_set)
        k = 392 in erap2_w

        cost = f"~${100 + ws * 3}" if ws <= 50 else "recombinant"
        spps = "YES" if ws <= 50 else "NO"

        print(f"\n  {ws}aa: residues {best['start']}-{best['end']}")
        print(f"    Coverage: {cov:.0f}% of contacts (<5A)")
        print(f"    Close contacts (<4A): {best['c4']}")
        print(f"    K392: {k}, DP residues: {dp_n}")
        print(f"    SPPS: {spps} ({cost})")
        print(f"    Seq: {seq}")

        results.append({
            'window': ws, 'start': best['start'], 'end': best['end'],
            'coverage_pct': round(cov, 1), 'contacts_5': best['score'],
            'contacts_4': best['c4'], 'k392': k, 'dp_contacts': dp_n,
            'spps_ok': ws <= 50, 'est_cost': cost, 'sequence': seq,
        })

    # Two-fragment analysis: can two short peptides capture the interface?
    print("\n" + "=" * 80)
    print("TWO-FRAGMENT APPROACH: Best pair of short peptides")
    print("=" * 80)

    best_pair = {'score': 0}
    for s1 in range(1, 70):
        for l1 in [12, 15, 18, 20]:
            e1 = min(s1 + l1 - 1, 92)
            c1 = sum(residue_data.get(r, {'contacts_5': 0})['contacts_5'] for r in range(s1, e1 + 1))
            for s2 in range(e1 + 3, 80):  # at least 2 residue gap
                for l2 in [12, 15, 18, 20]:
                    e2 = min(s2 + l2 - 1, 92)
                    c2 = sum(residue_data.get(r, {'contacts_5': 0})['contacts_5'] for r in range(s2, e2 + 1))
                    total = c1 + c2
                    if total > best_pair['score']:
                        best_pair = {
                            'score': total, 'frag1': (s1, e1, l1, c1),
                            'frag2': (s2, e2, l2, c2),
                        }

    if best_pair['score'] > 0:
        f1 = best_pair['frag1']
        f2 = best_pair['frag2']
        cov = best_pair['score'] / total_contacts_5 * 100
        seq1 = ''.join(residue_data.get(r, {'aa': '?'})['aa'] for r in range(f1[0], f1[1] + 1))
        seq2 = ''.join(residue_data.get(r, {'aa': '?'})['aa'] for r in range(f2[0], f2[1] + 1))

        # Check what each fragment contacts
        e2_f1 = set()
        for r in range(f1[0], f1[1] + 1):
            if r in residue_data:
                e2_f1.update(residue_data[r]['erap2_contacts'])
        e2_f2 = set()
        for r in range(f2[0], f2[1] + 1):
            if r in residue_data:
                e2_f2.update(residue_data[r]['erap2_contacts'])

        total_len = f1[2] + f2[2]
        linked_len = f1[2] + 4 + f2[2]  # GSGS linker

        print(f"  Best pair: {cov:.0f}% coverage")
        print(f"  Fragment 1: res {f1[0]}-{f1[1]} ({f1[2]}aa): {seq1}")
        print(f"    Contacts: {f1[3]}, ERAP2 res: {sorted(e2_f1)}")
        print(f"  Fragment 2: res {f2[0]}-{f2[1]} ({f2[2]}aa): {seq2}")
        print(f"    Contacts: {f2[3]}, ERAP2 res: {sorted(e2_f2)}")
        print(f"  As linked peptide: {seq1}-GSGS-{seq2} ({linked_len}aa)")
        print(f"  SPPS: {'YES' if linked_len <= 50 else 'NO'} (~${100 + linked_len * 3})")

    # Full sequence reference
    full = ''.join(residue_data.get(r, {'aa': '?'})['aa'] for r in range(1, 93))
    print(f"\nFull V2 (92aa): {full}")

    # Save
    output = {
        'binder_sequence': full,
        'total_residues': len(binder_residues),
        'interface_residues': len(interface_res),
        'total_contacts_5A': total_contacts_5,
        'sliding_windows': results,
        'segments': [{'start': s[0], 'end': s[-1], 'length': s[-1] - s[0] + 1} for s in segments],
    }
    with open(OUTPUT, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTPUT}")


if __name__ == '__main__':
    main()
