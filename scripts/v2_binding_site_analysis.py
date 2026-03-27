"""
V2 Lead Binder Binding Site Analysis

Parses the Boltz-2 CIF complex for Y87A_Y89A vs ERAP2 and confirms
the binder contacts the channel mouth near the active site.

Output:
  1. Full contact map (all ERAP2 residues within 5A of binder)
  2. Key residue check (K392, channel hotspots, zinc site)
  3. Binding site classification
  4. Contact table

Usage: python scripts/v2_binding_site_analysis.py
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import numpy as np
from pathlib import Path
from collections import defaultdict

PROJECT = Path(__file__).resolve().parent.parent
CIF_FILE = PROJECT / "data" / "results" / "y87a_cif_files" / "n248_trim_c5_Y87A_Y89A_erap2.cif"

CONTACT_CUTOFF = 5.0  # Angstroms
CLOSE_CUTOFF = 3.5    # close contact / H-bond range

# Cropped residue numbering offset
OFFSET = 349  # cropped res 1 = full ERAP2 res 350

# Key residue sets (FULL ERAP2 numbering)
VARIANT_SITE = {392}
CHANNEL_HOTSPOTS = {393, 401, 402, 403, 404, 405, 406}
ZINC_SITE = {383, 387, 455}  # zinc coordination sphere
CHANNEL_MOUTH = set(range(365, 420))  # broad channel entrance
DIVERGENT_PATCH_1 = set(range(353, 368))
DIVERGENT_PATCH_2 = set(range(400, 415))
DEEP_CHANNEL = set(range(440, 500))
DOMAIN_IV = set(range(500, 700))

# Binder sequence for reference
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


def main():
    print(f"Parsing: {CIF_FILE.name}")
    all_atoms = parse_cif_atoms(str(CIF_FILE))
    print(f"Total atoms: {len(all_atoms)}")

    # Identify chains
    chains = defaultdict(int)
    for a in all_atoms:
        chains[a['chain']] += 1
    chain_ids = sorted(chains.keys(), key=lambda c: chains[c], reverse=True)
    target_chain = chain_ids[0]
    binder_chain = chain_ids[1]
    print(f"Target (ERAP2): chain {target_chain} ({chains[target_chain]} atoms)")
    print(f"Binder: chain {binder_chain} ({chains[binder_chain]} atoms)")

    # Get heavy atoms only
    target_atoms = [a for a in all_atoms if a['chain'] == target_chain and a['element'] != 'H']
    binder_atoms = [a for a in all_atoms if a['chain'] == binder_chain and a['element'] != 'H']

    # Determine numbering scheme
    target_resnums = sorted(set(a['resnum'] for a in target_atoms))
    if target_resnums[0] < 100:
        numbering = "cropped"
        print(f"Numbering: cropped (res {target_resnums[0]}-{target_resnums[-1]}, +{OFFSET} for full)")
    else:
        numbering = "full"
        print(f"Numbering: full ERAP2 ({target_resnums[0]}-{target_resnums[-1]})")

    def to_full(resnum):
        return resnum + OFFSET if numbering == "cropped" else resnum

    # Build binder coordinate array for fast distance computation
    b_coords = np.array([[a['x'], a['y'], a['z']] for a in binder_atoms])

    # For each target residue, find minimum distance to any binder atom
    target_by_res = defaultdict(list)
    for a in target_atoms:
        target_by_res[a['resnum']].append(a)

    contacts = []
    for resnum in sorted(target_by_res.keys()):
        res_atoms = target_by_res[resnum]
        t_coords = np.array([[a['x'], a['y'], a['z']] for a in res_atoms])

        # Compute all pairwise distances
        min_dist = float('inf')
        closest_t_atom = None
        closest_b_atom = None
        closest_b_res = None
        for i, tc in enumerate(t_coords):
            dists = np.linalg.norm(b_coords - tc, axis=1)
            idx = np.argmin(dists)
            if dists[idx] < min_dist:
                min_dist = dists[idx]
                closest_t_atom = res_atoms[i]
                closest_b_atom = binder_atoms[idx]

        full_resnum = to_full(resnum)
        resname = res_atoms[0]['resname']
        aa1 = AA3TO1.get(resname, '?')

        # Classify region
        is_variant = full_resnum in VARIANT_SITE
        is_hotspot = full_resnum in CHANNEL_HOTSPOTS
        is_zinc = full_resnum in ZINC_SITE
        is_channel_mouth = full_resnum in CHANNEL_MOUTH
        is_dp1 = full_resnum in DIVERGENT_PATCH_1
        is_dp2 = full_resnum in DIVERGENT_PATCH_2
        is_deep = full_resnum in DEEP_CHANNEL
        is_domiv = full_resnum in DOMAIN_IV

        if is_variant:
            region = "VARIANT_SITE"
        elif is_zinc:
            region = "ZINC_SITE"
        elif is_hotspot:
            region = "CHANNEL_HOTSPOT"
        elif is_dp1:
            region = "DIVERGENT_PATCH_1"
        elif is_dp2:
            region = "DIVERGENT_PATCH_2"
        elif is_channel_mouth:
            region = "CHANNEL_MOUTH"
        elif is_deep:
            region = "DEEP_CHANNEL"
        elif is_domiv:
            region = "DOMAIN_IV"
        else:
            region = "OTHER"

        b_resnum = closest_b_atom['resnum']
        b_aa = AA3TO1.get(closest_b_atom['resname'], '?')

        contacts.append({
            "erap2_resnum": full_resnum,
            "erap2_crop_resnum": resnum,
            "erap2_aa": aa1,
            "erap2_resname": resname,
            "min_dist": round(min_dist, 2),
            "is_contact": min_dist <= CONTACT_CUTOFF,
            "is_close": min_dist <= CLOSE_CUTOFF,
            "region": region,
            "erap2_atom": closest_t_atom['atomname'],
            "binder_resnum": b_resnum,
            "binder_aa": b_aa,
            "binder_atom": closest_b_atom['atomname'],
        })

    # =====================================================
    # REPORT
    # =====================================================
    SEP = "=" * 85
    DASH = "-" * 85

    print(f"\n{SEP}")
    print("1. CONTACT MAP — ERAP2 residues within 5A of binder")
    print(SEP)

    contact_residues = [c for c in contacts if c['is_contact']]
    close_residues = [c for c in contacts if c['is_close']]
    print(f"Total ERAP2 residues: {len(contacts)}")
    print(f"Contacted (<{CONTACT_CUTOFF}A): {len(contact_residues)}")
    print(f"Close contacts (<{CLOSE_CUTOFF}A): {len(close_residues)}")

    print(f"\n{'ERAP2':>6} {'AA':>3} {'Dist':>6} {'Region':>20} {'Binder':>8} {'AA':>3} {'Atoms':>20}")
    print(DASH)
    for c in contact_residues:
        close_mark = " *" if c['is_close'] else ""
        atoms = f"{c['erap2_atom']}<->{c['binder_atom']}"
        print(f"{c['erap2_resnum']:>6} {c['erap2_aa']:>3} {c['min_dist']:>5.2f}{close_mark:2s} "
              f"{c['region']:>20} {c['binder_resnum']:>8} {c['binder_aa']:>3} {atoms:>20}")

    # =====================================================
    print(f"\n{SEP}")
    print("2. KEY RESIDUE CHECK")
    print(SEP)

    key_residues = {
        392: "K/N variant site",
        393: "E393 zinc coordinator",
        401: "Channel hotspot",
        402: "Channel hotspot",
        403: "A403 (ERAP2-unique)",
        404: "Channel hotspot",
        405: "Channel hotspot",
        406: "A406 (ERAP2 vs IRAP-K500)",
        383: "Zinc binding (H383)",
        387: "Zinc binding",
        455: "Zinc binding",
    }

    print(f"\n{'Residue':>8} {'Description':>30} {'Dist':>7} {'Contact':>8} {'Close':>6}")
    print(DASH)
    for resnum, desc in sorted(key_residues.items()):
        match = next((c for c in contacts if c['erap2_resnum'] == resnum), None)
        if match:
            is_c = "YES" if match['is_contact'] else "no"
            is_cl = "YES" if match['is_close'] else "no"
            print(f"{resnum:>8} {desc:>30} {match['min_dist']:>7.2f} {is_c:>8} {is_cl:>6}")
        else:
            print(f"{resnum:>8} {desc:>30} {'N/A':>7} {'N/A':>8} {'N/A':>6}")

    # =====================================================
    print(f"\n{SEP}")
    print("3. BINDING SITE CLASSIFICATION")
    print(SEP)

    # Count contacts by region
    region_counts = defaultdict(int)
    region_close = defaultdict(int)
    for c in contact_residues:
        region_counts[c['region']] += 1
        if c['is_close']:
            region_close[c['region']] += 1

    print(f"\n{'Region':>25} {'Contacts':>10} {'Close':>8}")
    print(DASH[:50])
    for region in ["VARIANT_SITE", "CHANNEL_HOTSPOT", "ZINC_SITE", "CHANNEL_MOUTH",
                   "DIVERGENT_PATCH_1", "DIVERGENT_PATCH_2", "DEEP_CHANNEL", "DOMAIN_IV", "OTHER"]:
        if region_counts[region] > 0 or region in ["VARIANT_SITE", "CHANNEL_HOTSPOT", "ZINC_SITE"]:
            print(f"{region:>25} {region_counts[region]:>10} {region_close[region]:>8}")

    # Classification logic
    has_392 = any(c['erap2_resnum'] == 392 and c['is_contact'] for c in contacts)
    has_hotspots = sum(1 for c in contact_residues if c['region'] == "CHANNEL_HOTSPOT")
    has_dp1 = sum(1 for c in contact_residues if c['region'] == "DIVERGENT_PATCH_1")
    has_dp2 = sum(1 for c in contact_residues if c['region'] == "DIVERGENT_PATCH_2")
    has_deep = sum(1 for c in contact_residues if c['region'] == "DEEP_CHANNEL")
    has_domiv = sum(1 for c in contact_residues if c['region'] == "DOMAIN_IV")
    total_contacts = len(contact_residues)

    print(f"\n--- VERDICT ---")
    if has_392 and has_hotspots >= 2:
        print(f"CHANNEL MOUTH (GOOD)")
        print(f"  Contacts K392: YES")
        print(f"  Channel hotspots (401-406): {has_hotspots} residues contacted")
        print(f"  Divergent patch 1 (353-367): {has_dp1} residues")
        print(f"  Divergent patch 2 (400-414): {has_dp2} residues")
        if has_deep > 5:
            print(f"  Also extends into deep channel: {has_deep} residues (bridge topology)")
        verdict = "CHANNEL_MOUTH"
    elif has_domiv > total_contacts * 0.5:
        print(f"WRONG SITE (BAD) — majority of contacts on domain IV")
        verdict = "WRONG_SITE"
    elif total_contacts > 0 and has_392 == 0 and has_hotspots == 0:
        print(f"NONSPECIFIC (BAD) — no channel residues contacted")
        verdict = "NONSPECIFIC"
    else:
        print(f"PARTIAL — some channel contacts but incomplete")
        verdict = "PARTIAL"

    # =====================================================
    print(f"\n{SEP}")
    print("4. FULL CONTACT TABLE")
    print(SEP)
    print(f"\n{'ERAP2':>6} {'AA':>3} {'MinDist':>8} {'Channel':>8} {'Region':>20}")
    print(DASH)
    for c in sorted(contacts, key=lambda x: x['min_dist']):
        if c['min_dist'] > 8.0:
            continue  # skip very distant
        is_ch = "Y" if c['region'] in ("VARIANT_SITE", "CHANNEL_HOTSPOT", "ZINC_SITE",
                                        "CHANNEL_MOUTH", "DIVERGENT_PATCH_1", "DIVERGENT_PATCH_2") else "N"
        marker = " <--" if c['is_contact'] else ""
        print(f"{c['erap2_resnum']:>6} {c['erap2_aa']:>3} {c['min_dist']:>8.2f} {is_ch:>8} {c['region']:>20}{marker}")

    # Summary stats
    print(f"\n{SEP}")
    print("SUMMARY")
    print(SEP)
    print(f"Verdict: {verdict}")
    print(f"Total contacts (<5A): {len(contact_residues)} ERAP2 residues")
    print(f"Close contacts (<3.5A): {len(close_residues)} ERAP2 residues")
    print(f"K392 contacted: {'YES (%.2fA)' % next((c['min_dist'] for c in contacts if c['erap2_resnum'] == 392), 99) if has_392 else 'NO'}")
    print(f"Channel hotspots (401-406): {has_hotspots} contacted")
    print(f"Divergent patches: {has_dp1 + has_dp2} residues ({has_dp1} DP1 + {has_dp2} DP2)")
    print(f"Deep channel: {has_deep} residues")
    print(f"Interface spans: {min(c['erap2_resnum'] for c in contact_residues)}-{max(c['erap2_resnum'] for c in contact_residues)} (full ERAP2 numbering)")


if __name__ == "__main__":
    main()
