"""V2-Mini stapling analysis: find helix positions and optimal staple sites."""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import numpy as np
from scipy.spatial.distance import cdist
import json, os

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(PROJECT, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
OUTPUT = os.path.join(PROJECT, "data/results/v43_validation/terminal_e/v2_staple_analysis.json")

SEQ = "FKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNIL"  # res 6-40


def parse_cif(path, chain):
    atoms, headers = [], []
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


def get_coord(atoms, resnum, atomname):
    for a in atoms:
        if a['resnum'] == resnum and a['atomname'] == atomname:
            return np.array([a['x'], a['y'], a['z']])
    return None


def calc_dihedral(p1, p2, p3, p4):
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    norm1 = np.linalg.norm(n1)
    norm2 = np.linalg.norm(n2)
    if norm1 < 1e-10 or norm2 < 1e-10:
        return None
    n1 /= norm1
    n2 /= norm2
    m = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2)))


def main():
    binder = parse_cif(CIF, 'B')
    erap2 = parse_cif(CIF, 'A')

    coords_b = np.array([[a['x'], a['y'], a['z']] for a in binder])
    coords_e = np.array([[a['x'], a['y'], a['z']] for a in erap2])
    dists = cdist(coords_b, coords_e)

    print("=" * 85)
    print("V2-MINI-35 STAPLING ANALYSIS")
    print("Sequence: " + SEQ)
    print("Binder residues 6-40 (35aa)")
    print("=" * 85)

    # Compute phi/psi and contacts for each residue
    print("\n--- Secondary Structure + Interface Map ---")
    print(f"{'Pos':>4} {'Res':>4} {'AA':>3} {'Phi':>8} {'Psi':>8} {'SS':>7} {'Cnt':>5} {'Face':>8}")
    print("-" * 55)

    ss_data = []
    for idx in range(35):
        resnum = idx + 6
        pos = idx + 1
        aa = SEQ[idx]

        C_prev = get_coord(binder, resnum - 1, 'C')
        N_i = get_coord(binder, resnum, 'N')
        CA_i = get_coord(binder, resnum, 'CA')
        C_i = get_coord(binder, resnum, 'C')
        N_next = get_coord(binder, resnum + 1, 'N')

        phi = calc_dihedral(C_prev, N_i, CA_i, C_i) if all(
            x is not None for x in [C_prev, N_i, CA_i, C_i]) else None
        psi = calc_dihedral(N_i, CA_i, C_i, N_next) if all(
            x is not None for x in [N_i, CA_i, C_i, N_next]) else None

        # SS classification
        ss = '?'
        if phi is not None and psi is not None:
            if -80 < phi < -40 and -60 < psi < -20:
                ss = 'helix'
            elif (-160 < phi < -100) and (100 < psi < 180 or -60 < psi < -10):
                ss = 'sheet'
            else:
                ss = 'coil'

        # Interface contacts
        b_idx = [i for i, a in enumerate(binder) if a['resnum'] == resnum]
        contacts = int((dists[b_idx] < 5.0).any(axis=0).sum()) if b_idx else 0

        # Determine face: interface (buried) vs solvent
        face = 'BURIED' if contacts >= 5 else 'surface' if contacts > 0 else 'solvent'

        ss_data.append({
            'pos': pos, 'resnum': resnum, 'aa': aa,
            'phi': round(phi, 1) if phi else None,
            'psi': round(psi, 1) if psi else None,
            'ss': ss, 'contacts': contacts, 'face': face,
        })

        phi_s = f"{phi:>7.1f}" if phi else "    N/A"
        psi_s = f"{psi:>7.1f}" if psi else "    N/A"
        marker = " <-helix" if ss == 'helix' else ""
        print(f"{pos:>4} {resnum:>4} {aa:>3} {phi_s} {psi_s} {ss:>7} {contacts:>5} {face:>8}{marker}")

    # Identify helical stretches
    print("\n--- Helical Segments ---")
    helix_runs = []
    current_run = []
    for d in ss_data:
        if d['ss'] == 'helix':
            current_run.append(d)
        else:
            if len(current_run) >= 3:
                helix_runs.append(current_run)
            current_run = []
    if len(current_run) >= 3:
        helix_runs.append(current_run)

    if not helix_runs:
        print("  No helical segments >= 3 residues found.")
        print("  The V2-Mini-35 may be a beta-sheet or mixed coil structure.")
        print("  Stapling is NOT appropriate for non-helical peptides.")
    else:
        for run in helix_runs:
            s = run[0]['pos']
            e = run[-1]['pos']
            seg_seq = SEQ[s - 1:e]
            total_cnt = sum(d['contacts'] for d in run)
            solvent_pos = [d['pos'] for d in run if d['contacts'] <= 2]
            print(f"  Helix: pos {s}-{e} ({e - s + 1}aa): {seg_seq}")
            print(f"    Total contacts: {total_cnt}")
            print(f"    Solvent-facing positions (low contact): {solvent_pos}")

    # Stapling candidates: i, i+4 pairs
    print("\n--- Staple Candidates (i, i+4) ---")
    print("  Hydrocarbon staple replaces both positions with (S)-pentenyl-alanine (S5)")
    print("  IDEAL: both helical, both solvent-facing (contacts <= 2)")
    print("  GOOD: both helical, at least one low-contact")
    print()
    print(f"{'Pair':>8} {'AA':>6} {'SS':>12} {'Cnt':>10} {'Grade':>8}")
    print("-" * 50)

    candidates = []
    for i in range(len(ss_data) - 4):
        d_i = ss_data[i]
        d_i4 = ss_data[i + 4]
        both_h = d_i['ss'] == 'helix' and d_i4['ss'] == 'helix'
        low_both = d_i['contacts'] <= 2 and d_i4['contacts'] <= 2
        low_one = d_i['contacts'] <= 3 or d_i4['contacts'] <= 3

        if both_h and low_both:
            grade = 'IDEAL'
        elif both_h and low_one:
            grade = 'GOOD'
        elif both_h:
            grade = 'OK'
        else:
            grade = 'SKIP'

        if grade in ('IDEAL', 'GOOD', 'OK'):
            p1 = d_i['pos']
            p2 = d_i4['pos']
            candidates.append({
                'pos_i': p1, 'pos_i4': p2,
                'aa_i': d_i['aa'], 'aa_i4': d_i4['aa'],
                'contacts_lost': d_i['contacts'] + d_i4['contacts'],
                'grade': grade,
            })
            print(f"{p1:>3}-{p2:<3} {d_i['aa']}+{d_i4['aa']:>3} "
                  f"{'H+H':>12} {d_i['contacts']:>4}+{d_i4['contacts']:<4} {grade:>8}")

    # Also check i, i+7
    print(f"\n--- Extended Staple Candidates (i, i+7) ---")
    for i in range(len(ss_data) - 7):
        d_i = ss_data[i]
        d_i7 = ss_data[i + 7]
        both_h = d_i['ss'] == 'helix' and d_i7['ss'] == 'helix'
        low_both = d_i['contacts'] <= 2 and d_i7['contacts'] <= 2

        if both_h and (low_both or d_i['contacts'] <= 3 or d_i7['contacts'] <= 3):
            grade = 'IDEAL' if low_both else 'GOOD'
            p1 = d_i['pos']
            p2 = d_i7['pos']
            print(f"  {p1:>3}-{p2:<3} {d_i['aa']}+{d_i7['aa']} "
                  f"cnt={d_i['contacts']}+{d_i7['contacts']} {grade}")

    # Final recommendation
    print("\n" + "=" * 85)
    print("SYNTHESIS RECOMMENDATIONS")
    print("=" * 85)

    ideal = [c for c in candidates if c['grade'] == 'IDEAL']
    good = [c for c in candidates if c['grade'] == 'GOOD']
    best = ideal if ideal else good

    if best:
        b = best[0]
        stapled = list(SEQ)
        stapled[b['pos_i'] - 1] = 'S5'
        stapled[b['pos_i4'] - 1] = 'S5'

        print(f"\n  BEST STAPLE: positions {b['pos_i']} and {b['pos_i4']}")
        print(f"  Replace: {b['aa_i']}{b['pos_i']} and {b['aa_i4']}{b['pos_i4']} with (S)-pentenyl-Ala")
        print(f"  Contacts lost: {b['contacts_lost']} (out of ~290 total)")
        print(f"  Original:  Ac-{SEQ}-NH2")

        # Show with markers
        marked = list(SEQ)
        marked[b['pos_i'] - 1] = f"[S5/{b['aa_i']}]"
        marked[b['pos_i4'] - 1] = f"[S5/{b['aa_i4']}]"
        print(f"  Stapled:   Ac-{''.join(marked)}-NH2")
    else:
        print("\n  No suitable stapling positions found.")
        print("  The peptide may not have a long enough helix for stapling.")

    print(f"\n  --- ORDER LIST (Tier 1, all SPPS) ---")
    print(f"  1. Ac-VAGSAF-NH2                              (6aa,  ~$120, V4 lead)")
    print(f"  2. Ac-IAFSAF-NH2                              (6aa,  ~$120, V4 selectivity lead)")
    print(f"  3. FASGAV                                     (6aa,  ~$100, scrambled control)")
    print(f"  4. Ac-{SEQ}-NH2  (35aa, ~$250, V2-Mini)")
    if best:
        b = best[0]
        print(f"  5. Ac-V2-Mini-STAPLED(pos{b['pos_i']},{b['pos_i4']})-NH2  (35aa, ~$350, stapled V2-Mini)")
    print(f"  6. Ac-VAGSAF-GSGSGS-ALAPYIP-NH2               (19aa, ~$180, peptide-PROTAC)")
    print(f"  7. LEEYLKNLPKVVDMLVDLYS-GSGS-TNILVKDDKFYAIDFGSAYI-NH2  (44aa, ~$280, two-frag linked)")
    total = 120 + 120 + 100 + 250 + 350 + 180 + 280
    print(f"\n  TOTAL (all 7): ~${total}")
    print(f"  Budget pick (items 1-4,6): ~${120+120+100+250+180}")

    # Save
    output = {
        'sequence': SEQ,
        'ss_data': ss_data,
        'helix_segments': [{'start': r[0]['pos'], 'end': r[-1]['pos'],
                            'seq': SEQ[r[0]['pos']-1:r[-1]['pos']]}
                           for r in helix_runs],
        'staple_candidates': candidates,
        'best_staple': best[0] if best else None,
    }
    with open(OUTPUT, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
