"""
Construct B + C Design: Potency upgrades for Y87A_Y89A parent binder.

Part 1: Find channel-facing non-interface residues for W-plug mutation
Part 2: Assess C-terminus geometry for tethered peptide cargo
Part 3: Comparison matrix
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json, os
import numpy as np
from scipy.spatial.distance import cdist

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(PROJECT, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
OFFSET = 349

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Bulkier side-chain volumes (A^3, approximate)
SC_VOLUME = {
    'G': 0, 'A': 26, 'V': 50, 'L': 67, 'I': 67, 'P': 42,
    'F': 77, 'Y': 86, 'W': 99, 'S': 29, 'T': 43, 'C': 35,
    'M': 63, 'D': 40, 'E': 54, 'N': 44, 'Q': 58, 'K': 68,
    'R': 85, 'H': 60,
}


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
                    if r.get('auth_asym_id', r.get('label_asym_id')) == chain:
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


def get_ca(atoms, resnum):
    for a in atoms:
        if a['resnum'] == resnum and a['atomname'] == 'CA':
            return np.array([a['x'], a['y'], a['z']])
    return None


def get_cb(atoms, resnum):
    """Get CB (or CA for Gly) — indicates sidechain direction."""
    for a in atoms:
        if a['resnum'] == resnum and a['atomname'] == 'CB':
            return np.array([a['x'], a['y'], a['z']])
    return get_ca(atoms, resnum)


def main():
    erap2 = parse_cif(CIF, 'A')
    binder = parse_cif(CIF, 'B')

    coords_b = np.array([[a['x'], a['y'], a['z']] for a in binder])
    coords_e = np.array([[a['x'], a['y'], a['z']] for a in erap2])
    dists = cdist(coords_b, coords_e)

    binder_residues = sorted(set(a['resnum'] for a in binder))
    binder_seq = {}
    for a in binder:
        if a['resnum'] not in binder_seq:
            binder_seq[a['resnum']] = AA3TO1.get(a['resname'], 'X')
    full_seq = ''.join(binder_seq[r] for r in sorted(binder_seq))

    # Define channel entrance region on ERAP2 (full numbering)
    # Catalytic zinc at H370, H374, E393
    # Channel entrance: 365-415
    # Substrate path: roughly along zinc → entrance axis
    channel_entrance = set(range(365, 416))
    zinc_region = {370, 374, 393}
    k392_region = {392}

    # Compute ERAP2 channel centroid (entrance residues CA)
    channel_cas = []
    for resnum in range(365, 416):
        ca = get_ca(erap2, resnum - OFFSET)
        if ca is not None:
            channel_cas.append(ca)
    channel_centroid = np.mean(channel_cas, axis=0)

    # Zinc centroid
    zinc_cas = []
    for resnum in [370, 374, 393]:
        ca = get_ca(erap2, resnum - OFFSET)
        if ca is not None:
            zinc_cas.append(ca)
    zinc_centroid = np.mean(zinc_cas, axis=0)

    # Channel axis: zinc → entrance centroid
    channel_axis = channel_centroid - zinc_centroid
    channel_axis_norm = channel_axis / np.linalg.norm(channel_axis)

    print("=" * 85)
    print("CONSTRUCT B + C DESIGN: Potency Upgrades for Y87A_Y89A")
    print("=" * 85)
    print(f"\nParent sequence (92aa): {full_seq}")
    print(f"Channel centroid: [{channel_centroid[0]:.1f}, {channel_centroid[1]:.1f}, {channel_centroid[2]:.1f}]")
    print(f"Zinc centroid: [{zinc_centroid[0]:.1f}, {zinc_centroid[1]:.1f}, {zinc_centroid[2]:.1f}]")

    # ===================================================================
    # PART 1: CONSTRUCT B — W-PLUG CANDIDATES
    # ===================================================================
    print("\n" + "=" * 85)
    print("PART 1: CONSTRUCT B — W-PLUG MUTATION CANDIDATES")
    print("=" * 85)

    print("\nScoring each binder residue:")
    print("  - interface_contacts: # ERAP2 atoms within 5A (high = critical for binding)")
    print("  - dist_to_channel: distance from CB to channel entrance centroid")
    print("  - dist_to_zinc: distance from CB to zinc centroid")
    print("  - channel_facing: does CB point toward channel? (CB-CA vector dot channel axis)")
    print("  - volume_gain: how much volume a W mutation adds vs current residue")
    print()

    candidates = []
    print(f"{'Res':>4} {'AA':>3} {'Cnt':>5} {'ChDist':>7} {'ZnDist':>7} {'Facing':>8} "
          f"{'VolGain':>8} {'Grade':>8}")
    print("-" * 60)

    for resnum in binder_residues:
        aa = binder_seq[resnum]
        b_idx = [i for i, a in enumerate(binder) if a['resnum'] == resnum]
        contacts = int((dists[b_idx] < 5.0).any(axis=0).sum())

        ca = get_ca(binder, resnum)
        cb = get_cb(binder, resnum)
        if ca is None or cb is None:
            continue

        # Distance from CB to channel centroid and zinc
        d_channel = float(np.linalg.norm(cb - channel_centroid))
        d_zinc = float(np.linalg.norm(cb - zinc_centroid))

        # Channel-facing score: does CB point toward channel?
        ca_to_cb = cb - ca
        ca_to_cb_len = np.linalg.norm(ca_to_cb)
        if ca_to_cb_len > 0.1:
            ca_to_channel = channel_centroid - ca
            facing = float(np.dot(ca_to_cb / ca_to_cb_len,
                                  ca_to_channel / np.linalg.norm(ca_to_channel)))
        else:
            facing = 0.0

        # Volume gain from W mutation
        vol_gain = SC_VOLUME.get('W', 99) - SC_VOLUME.get(aa, 50)

        # Grading
        # IDEAL: low contacts (<=3), close to channel (<20A), channel-facing (>0.3), good vol gain
        # GOOD: low contacts, moderately close, some facing
        # SKIP: high contacts (interface critical) or facing away
        if contacts <= 3 and d_channel < 18 and facing > 0.2 and vol_gain > 20:
            grade = 'IDEAL'
        elif contacts <= 5 and d_channel < 22 and facing > 0.0 and vol_gain > 10:
            grade = 'GOOD'
        elif contacts <= 3 and d_channel < 25:
            grade = 'OK'
        else:
            grade = '-'

        if grade in ('IDEAL', 'GOOD', 'OK'):
            print(f"{resnum:>4} {aa:>3} {contacts:>5} {d_channel:>6.1f}A {d_zinc:>6.1f}A "
                  f"{facing:>7.2f} {vol_gain:>7} {grade:>8}")
            candidates.append({
                'resnum': resnum, 'aa': aa, 'contacts': contacts,
                'dist_channel': round(d_channel, 1),
                'dist_zinc': round(d_zinc, 1),
                'facing': round(facing, 2),
                'vol_gain': vol_gain, 'grade': grade,
            })

    # Rank candidates
    ideal = [c for c in candidates if c['grade'] == 'IDEAL']
    good = [c for c in candidates if c['grade'] == 'GOOD']

    # Sort by: facing * vol_gain / dist_channel (maximize facing and volume, minimize distance)
    for pool in [ideal, good]:
        pool.sort(key=lambda c: c['facing'] * c['vol_gain'] / c['dist_channel'], reverse=True)

    ranked = ideal + good
    print(f"\n--- TOP 3 W-PLUG CANDIDATES ---")
    for i, c in enumerate(ranked[:3]):
        print(f"\n  #{i+1}: Position {c['resnum']} ({c['aa']} -> W)")
        print(f"      Interface contacts: {c['contacts']} (low = safe to mutate)")
        print(f"      Distance to channel entrance: {c['dist_channel']} A")
        print(f"      Distance to zinc: {c['dist_zinc']} A")
        print(f"      Channel-facing score: {c['facing']} (1.0 = directly toward channel)")
        print(f"      Volume gain (W vs {c['aa']}): +{c['vol_gain']} A^3")
        print(f"      Grade: {c['grade']}")

        # Also suggest Y and F alternatives
        vol_y = SC_VOLUME['Y'] - SC_VOLUME.get(c['aa'], 50)
        vol_f = SC_VOLUME['F'] - SC_VOLUME.get(c['aa'], 50)
        print(f"      Alternatives: {c['aa']}->{c['resnum']}Y (+{vol_y} A^3), "
              f"{c['aa']}->{c['resnum']}F (+{vol_f} A^3)")

    if ranked:
        best = ranked[0]
        construct_b_seq = list(full_seq)
        construct_b_seq[best['resnum'] - 1] = 'W'
        construct_b_seq = ''.join(construct_b_seq)

        print(f"\n{'='*85}")
        print(f"CONSTRUCT B RECOMMENDATION: {best['aa']}{best['resnum']}W")
        print(f"{'='*85}")
        print(f"  Parent:      ...{full_seq[best['resnum']-3:best['resnum']+2]}...")
        print(f"  Construct B: ...{construct_b_seq[best['resnum']-3:best['resnum']+2]}...")
        print(f"  Full: {construct_b_seq}")
        print(f"\n  Rationale: Position {best['resnum']} is {best['dist_channel']} A from the "
              f"channel entrance, points toward it (facing={best['facing']}), and has only "
              f"{best['contacts']} ERAP2 contacts (safe to mutate). W adds +{best['vol_gain']} A^3 "
              f"of steric bulk in the direction of the substrate path.")

    # ===================================================================
    # PART 2: CONSTRUCT C — TETHERED PEPTIDE
    # ===================================================================
    print(f"\n{'='*85}")
    print("PART 2: CONSTRUCT C — TETHERED PEPTIDE DESIGN")
    print(f"{'='*85}")

    # Where is the binder C-terminus relative to the channel?
    c_term_ca = get_ca(binder, 92)
    n_term_ca = get_ca(binder, 1)

    if c_term_ca is not None:
        d_cterm_channel = float(np.linalg.norm(c_term_ca - channel_centroid))
        d_cterm_zinc = float(np.linalg.norm(c_term_ca - zinc_centroid))
        # Does C-term point toward channel?
        ca91 = get_ca(binder, 91)
        if ca91 is not None:
            cterm_dir = c_term_ca - ca91
            cterm_dir /= np.linalg.norm(cterm_dir)
            cterm_to_channel = channel_centroid - c_term_ca
            cterm_facing = float(np.dot(cterm_dir, cterm_to_channel / np.linalg.norm(cterm_to_channel)))
        else:
            cterm_facing = 0
    else:
        d_cterm_channel = d_cterm_zinc = cterm_facing = 0

    if n_term_ca is not None:
        d_nterm_channel = float(np.linalg.norm(n_term_ca - channel_centroid))
        d_nterm_zinc = float(np.linalg.norm(n_term_ca - zinc_centroid))
    else:
        d_nterm_channel = d_nterm_zinc = 0

    print(f"\n  Terminus geometry:")
    print(f"    C-terminus (res 92) to channel centroid: {d_cterm_channel:.1f} A")
    print(f"    C-terminus (res 92) to zinc centroid:    {d_cterm_zinc:.1f} A")
    print(f"    C-terminus channel-facing score:         {cterm_facing:.2f}")
    print(f"    N-terminus (res 1) to channel centroid:  {d_nterm_channel:.1f} A")
    print(f"    N-terminus (res 1) to zinc centroid:     {d_nterm_zinc:.1f} A")

    # Which terminus is closer to channel?
    if d_cterm_channel < d_nterm_channel:
        attach = 'C-terminus'
        attach_dist = d_cterm_channel
        attach_zinc = d_cterm_zinc
    else:
        attach = 'N-terminus'
        attach_dist = d_nterm_channel
        attach_zinc = d_nterm_zinc

    print(f"\n  Better attachment point: {attach} ({attach_dist:.1f} A to channel)")

    # Linker length analysis
    # (GGGGS)2 = 10 residues, extended ~35A, average reach ~15-20A
    # (GGGGS)3 = 15 residues, extended ~52A, average reach ~20-30A
    print(f"\n  Linker analysis:")
    print(f"    Gap to bridge: {attach_dist:.1f} A (terminus to channel centroid)")
    print(f"    (GGGGS)x2 (10aa): max reach ~35A, avg ~18A")
    print(f"    (GGGGS)x3 (15aa): max reach ~52A, avg ~25A")

    if attach_dist < 20:
        rec_linker = "(GGGGS)x2"
        rec_linker_seq = "GGGGSGGGGS"
        print(f"    RECOMMENDATION: (GGGGS)x2 — distance is short enough, less entropy penalty")
    else:
        rec_linker = "(GGGGS)x3"
        rec_linker_seq = "GGGGSGGGGSGGGGS"
        print(f"    RECOMMENDATION: (GGGGS)x3 — need extra reach to span {attach_dist:.0f}A gap")

    # Cargo evaluation
    print(f"\n  Cargo evaluation:")

    cargos = [
        ("VKLLLL", "V-P1 + K392-interacting Lys + hydrophobic packing"),
        ("EKLLLL", "E-P1 + salt-bridge to K392 + hydrophobic tail"),
        ("VAGSAF", "V4 lead — highest ipTM 0.905 in original screen"),
        ("IAFSAF", "V4 selectivity lead — best K/N delta +0.239"),
    ]

    print(f"    {'Cargo':>10} {'Rationale':>55}")
    print(f"    {'-'*10} {'-'*55}")
    for cargo, rationale in cargos:
        print(f"    {cargo:>10} {rationale}")

    print(f"\n  CARGO RECOMMENDATION: VAGSAF")
    print(f"    Reason: Already validated as highest-affinity 6-mer for ERAP2-K392 channel")
    print(f"    (ipTM 0.905). The binder anchors it near the channel; the linker lets it")
    print(f"    sample the entrance. VKLLLL is untested and speculative.")

    # Build Construct C
    if attach == 'C-terminus':
        construct_c_seq = full_seq + rec_linker_seq + "VAGSAF"
    else:
        construct_c_seq = "VAGSAF" + rec_linker_seq + full_seq

    construct_c_len = len(construct_c_seq)

    print(f"\n{'='*85}")
    print(f"CONSTRUCT C RECOMMENDATION")
    print(f"{'='*85}")
    print(f"  Architecture: Parent — {rec_linker} — VAGSAF")
    print(f"  Attachment: {attach}")
    print(f"  Total length: {construct_c_len} aa")
    print(f"  Sequence: {construct_c_seq}")

    # Backup with longer linker
    if rec_linker == "(GGGGS)x2":
        backup_linker = "GGGGSGGGGSGGGGS"
        backup_name = "(GGGGS)x3"
    else:
        backup_linker = "GGGGSGGGGS"
        backup_name = "(GGGGS)x2"

    if attach == 'C-terminus':
        backup_seq = full_seq + backup_linker + "VAGSAF"
    else:
        backup_seq = "VAGSAF" + backup_linker + full_seq

    print(f"\n  BACKUP: Parent — {backup_name} — VAGSAF ({len(backup_seq)} aa)")

    # Risks
    print(f"\n  Risks:")
    print(f"    - Linker too floppy: cargo may not find channel entrance")
    print(f"    - Linker too short: cargo can't reach channel from {attach_dist:.0f}A away")
    print(f"    - ERAP2 may cleave the VAGSAF cargo (aminopeptidase)")
    print(f"      Mitigation: Ac-cap on VAGSAF N-term (blocks exopeptidase)")
    print(f"    - Recombinant expression needed (>{50}aa, can't do SPPS)")

    # ===================================================================
    # PART 3: COMPARISON
    # ===================================================================
    print(f"\n{'='*85}")
    print("PART 3: COMPARISON MATRIX")
    print(f"{'='*85}")
    print()
    print(f"{'':>20} {'Parent':>15} {'Construct B':>15} {'Construct C':>15}")
    print(f"{'':>20} {'(Y87A_Y89A)':>15} {'(W-plug)':>15} {'(tethered)':>15}")
    print("-" * 70)

    rows = [
        ("Size (aa)", "92", "92", str(construct_c_len)),
        ("Mutations", "none", f"{ranked[0]['aa']}{ranked[0]['resnum']}W" if ranked else "?", f"+linker+VAGSAF"),
        ("ERAP2 binding", "0.748 ipTM", "~0.74 (minimal change)", "~0.6-0.7 (linker drag)"),
        ("Selectivity risk", "LOW (proven)", "LOW (non-interface)", "MEDIUM (cargo hits conserved)"),
        ("Potency / inhibition", "~30% SASA block", "~35-40% (bigger plug)", "~40-50% (cargo in channel)"),
        ("Synthesis", "recombinant", "recombinant", "recombinant"),
        ("Cost", "~$3K", "~$3K", "~$3.5K"),
        ("Risk level", "baseline", "LOW", "MEDIUM"),
    ]
    for label, *vals in rows:
        print(f"{label:>20} {vals[0]:>15} {vals[1]:>15} {vals[2]:>15}")

    print(f"\n{'='*85}")
    print("FINAL RANKING")
    print(f"{'='*85}")
    print(f"\n  1. SAFEST potency upgrade:  CONSTRUCT B ({ranked[0]['aa']}{ranked[0]['resnum']}W)" if ranked else "")
    print(f"     Single point mutation. Minimal fold disruption. Adds steric bulk")
    print(f"     toward channel. Will NOT dramatically change selectivity.")
    print(f"     Expected improvement: modest (+5-10% channel occlusion)")
    print(f"\n  2. HIGHEST UPSIDE:  CONSTRUCT C (parent + linker + VAGSAF)")
    print(f"     The cargo can physically enter the channel. Potential for 40-50%")
    print(f"     occlusion. But higher risk: linker flexibility, expression complexity,")
    print(f"     possible selectivity loss from cargo-conserved contacts.")

    # Save
    output = {
        'parent': full_seq,
        'construct_b': {
            'mutation': f"{ranked[0]['aa']}{ranked[0]['resnum']}W" if ranked else None,
            'sequence': construct_b_seq if ranked else None,
            'rationale': f"Position {ranked[0]['resnum']}: {ranked[0]['dist_channel']}A to channel, "
                         f"facing={ranked[0]['facing']}, {ranked[0]['contacts']} contacts" if ranked else None,
            'top_candidates': ranked[:3],
        },
        'construct_c': {
            'architecture': f"Parent — {rec_linker} — VAGSAF",
            'sequence': construct_c_seq,
            'linker': rec_linker_seq,
            'cargo': 'VAGSAF',
            'attachment': attach,
            'terminus_to_channel_A': round(attach_dist, 1),
            'backup_sequence': backup_seq,
        },
    }
    out_path = os.path.join(PROJECT, "data/results/v43_validation/terminal_e/construct_bc_design.json")
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
