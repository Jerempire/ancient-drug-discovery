"""
V2 Protein Binder Channel-Capping Analysis (Terminal E)

Analyzes whether the Y87A_Y89A binder covers the ERAP2 substrate channel
entrance by computing:
1. Contact footprint mapping (cropped → full ERAP2 numbering)
2. SASA of channel entrance residues ± binder
3. SASA of zinc coordination sphere ± binder
4. Spatial coverage: does the binder physically occlude the channel opening?

Dependencies: numpy, scipy, BioPython (for SASA via Shrake-Rupley)
Output: data/results/v43_validation/terminal_e/

Usage: python scripts/v2_channel_cap_analysis.py
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import numpy as np
from datetime import datetime, timezone

# === CONFIG ===
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF_FILE = os.path.join(PROJECT_DIR, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
ERAP2_PDB = os.path.join(PROJECT_DIR, "data/structures/erap2_wt_alphafold.pdb")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "data/results/v43_validation/terminal_e")
os.makedirs(OUTPUT_DIR, exist_ok=True)

OFFSET = 349  # cropped residue + 349 = full ERAP2 residue number

# Key structural regions (full ERAP2 numbering)
ZINC_COORD = {370, 374, 393}  # H370, H374, E393
CHANNEL_ENTRANCE = set(range(350, 420))  # broad entrance
DIVERGENT_PATCH_1 = set(range(353, 368))  # 353-367
DIVERGENT_PATCH_2 = set(range(400, 415))  # 400-414
K392_REGION = {392, 398}
CHANNEL_FLOOR = {392, 398}  # K392 selectivity handle
CHANNEL_WALL = {403, 406}   # IRAP/ERAP1 evasion
CHANNEL_CEILING = {412, 414}  # channel cap


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
                    atoms.append(record)
            elif in_atom and (line.startswith('#') or line.startswith('loop_')):
                if atoms:
                    break
    return atoms


def compute_sasa_shrake_rupley(coords, radii, probe_radius=1.4, n_points=100):
    """
    Shrake-Rupley SASA algorithm.
    coords: Nx3 array of atom positions
    radii: N array of van der Waals radii
    Returns: N array of per-atom SASA values
    """
    n_atoms = len(coords)
    sasa = np.zeros(n_atoms)

    # Generate points on unit sphere (Fibonacci lattice)
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    sphere_points = np.column_stack([
        np.sin(phi) * np.cos(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(phi)
    ])

    for i in range(n_atoms):
        ri = radii[i] + probe_radius
        test_points = coords[i] + ri * sphere_points

        # Check each test point against all other atoms
        accessible = np.ones(n_points, dtype=bool)
        for j in range(n_atoms):
            if i == j:
                continue
            rj = radii[j] + probe_radius
            dist_sq = np.sum((test_points - coords[j])**2, axis=1)
            accessible &= dist_sq > rj**2

        sasa[i] = 4 * np.pi * ri**2 * np.sum(accessible) / n_points

    return sasa


# Van der Waals radii by element
VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20,
    'ZN': 1.39, 'FE': 1.40, 'SE': 1.90
}


def get_vdw_radius(atom_record):
    """Get vdw radius from atom type."""
    element = atom_record.get('type_symbol', atom_record['atomname'][0])
    return VDW_RADII.get(element.upper(), 1.70)


def run_analysis():
    print("=" * 60)
    print("V2 PROTEIN BINDER CHANNEL-CAPPING ANALYSIS")
    print("Terminal E — Ancient Drug Discovery V4.3 Validation")
    print("=" * 60)

    # Parse CIF
    print("\n[1/5] Parsing CIF file...")
    all_atoms = parse_cif_atoms(CIF_FILE)
    erap2_atoms = [a for a in all_atoms if a['chain'] == 'A']
    binder_atoms = [a for a in all_atoms if a['chain'] == 'B']
    print(f"  ERAP2: {len(erap2_atoms)} atoms, Binder: {len(binder_atoms)} atoms")

    # === CONTACT ANALYSIS ===
    print("\n[2/5] Computing interface contacts (5A cutoff)...")
    from scipy.spatial.distance import cdist

    coords_e = np.array([[a['x'], a['y'], a['z']] for a in erap2_atoms])
    coords_b = np.array([[a['x'], a['y'], a['z']] for a in binder_atoms])

    dists = cdist(coords_e, coords_b)
    contact_mask = dists < 5.0

    # ERAP2 contact residues (full numbering)
    contact_erap2_idx = np.where(contact_mask.any(axis=1))[0]
    contact_erap2_full = sorted(set(erap2_atoms[i]['resnum'] + OFFSET for i in contact_erap2_idx))

    # Closest approach distances per ERAP2 residue
    erap2_res_min_dist = {}
    for i in range(len(erap2_atoms)):
        full_res = erap2_atoms[i]['resnum'] + OFFSET
        min_d = dists[i].min()
        if full_res not in erap2_res_min_dist or min_d < erap2_res_min_dist[full_res]:
            erap2_res_min_dist[full_res] = float(min_d)

    # Region overlap
    dp1_contacts = set(contact_erap2_full) & DIVERGENT_PATCH_1
    dp2_contacts = set(contact_erap2_full) & DIVERGENT_PATCH_2
    zinc_contacts = set(contact_erap2_full) & ZINC_COORD
    k392_contacts = set(contact_erap2_full) & K392_REGION
    entrance_contacts = set(contact_erap2_full) & CHANNEL_ENTRANCE

    print(f"  ERAP2 residues in contact: {len(contact_erap2_full)}")
    print(f"  Divergent patch 1 (353-367): {len(dp1_contacts)}/15 -> {sorted(dp1_contacts)}")
    print(f"  Divergent patch 2 (400-414): {len(dp2_contacts)}/15 -> {sorted(dp2_contacts)}")
    print(f"  Zinc coordination: {len(zinc_contacts)}/3 -> {sorted(zinc_contacts)}")
    print(f"  K392 region: {len(k392_contacts)}/2 -> {sorted(k392_contacts)}")
    print(f"  Channel entrance (350-419): {len(entrance_contacts)}/70")

    # === SASA ANALYSIS ===
    print("\n[3/5] Computing SASA of channel entrance residues...")

    # Select entrance residue atoms from ERAP2
    entrance_atom_idx = [i for i, a in enumerate(erap2_atoms)
                         if (a['resnum'] + OFFSET) in CHANNEL_ENTRANCE]

    # SASA with binder present (all atoms)
    all_coords = np.vstack([coords_e, coords_b])
    all_radii = np.array([get_vdw_radius(a) for a in erap2_atoms + binder_atoms])

    # SASA without binder (ERAP2 only)
    erap2_radii = np.array([get_vdw_radius(a) for a in erap2_atoms])

    print("  Computing apo SASA (ERAP2 only)...")
    sasa_apo = compute_sasa_shrake_rupley(coords_e, erap2_radii, n_points=92)

    print("  Computing holo SASA (ERAP2 + binder)...")
    sasa_holo = compute_sasa_shrake_rupley(all_coords, all_radii, n_points=92)
    sasa_holo_erap2 = sasa_holo[:len(erap2_atoms)]  # only ERAP2 atoms

    # Per-region SASA
    def region_sasa(region_set, label):
        idx = [i for i, a in enumerate(erap2_atoms) if (a['resnum'] + OFFSET) in region_set]
        if not idx:
            return {'apo': 0, 'holo': 0, 'reduction_pct': 0, 'label': label}
        apo_sum = sum(sasa_apo[i] for i in idx)
        holo_sum = sum(sasa_holo_erap2[i] for i in idx)
        reduction = (1 - holo_sum / apo_sum) * 100 if apo_sum > 0 else 0
        return {
            'apo_sasa': round(float(apo_sum), 1),
            'holo_sasa': round(float(holo_sum), 1),
            'reduction_pct': round(float(reduction), 1),
            'label': label
        }

    sasa_results = {
        'channel_entrance_350_419': region_sasa(CHANNEL_ENTRANCE, 'Channel entrance (350-419)'),
        'divergent_patch_1': region_sasa(DIVERGENT_PATCH_1, 'Divergent patch 1 (353-367)'),
        'divergent_patch_2': region_sasa(DIVERGENT_PATCH_2, 'Divergent patch 2 (400-414)'),
        'zinc_coordination': region_sasa(ZINC_COORD, 'Zinc coord (H370, H374, E393)'),
        'k392_region': region_sasa(K392_REGION, 'K392 selectivity (392, 398)'),
        'channel_floor': region_sasa(CHANNEL_FLOOR, 'Floor (K392, Y398)'),
        'channel_wall': region_sasa(CHANNEL_WALL, 'Wall (A403, A406)'),
        'channel_ceiling': region_sasa(CHANNEL_CEILING, 'Ceiling (Q412, D414)'),
    }

    print("\n  SASA Results (Angstrom^2):")
    print(f"  {'Region':<35} {'Apo':>8} {'Holo':>8} {'Reduction':>10}")
    print(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*10}")
    for key, val in sasa_results.items():
        print(f"  {val['label']:<35} {val['apo_sasa']:>8.1f} {val['holo_sasa']:>8.1f} {val['reduction_pct']:>9.1f}%")

    # === SPATIAL COVERAGE: Convex hull projection ===
    print("\n[4/5] Computing spatial channel occlusion...")

    # Define channel axis: roughly along the vector from zinc center to entrance
    # H370 is cropped res 21, entrance is around cropped res 1-5
    zinc_center_idx = [i for i, a in enumerate(erap2_atoms)
                       if a['resnum'] + OFFSET in {370, 374, 393}]
    entrance_idx = [i for i, a in enumerate(erap2_atoms)
                    if (a['resnum'] + OFFSET) in {350, 351, 352, 353}]

    if zinc_center_idx and entrance_idx:
        zinc_center = coords_e[zinc_center_idx].mean(axis=0)
        entrance_center = coords_e[entrance_idx].mean(axis=0)
        channel_axis = entrance_center - zinc_center
        channel_axis /= np.linalg.norm(channel_axis)

        # Project binder atoms onto plane perpendicular to channel axis at entrance
        binder_projections = []
        for bc in coords_b:
            # Project onto plane at entrance_center
            v = bc - entrance_center
            proj_along = np.dot(v, channel_axis)
            proj_on_plane = v - proj_along * channel_axis
            binder_projections.append(proj_on_plane)

        binder_proj = np.array(binder_projections)

        # Similarly project entrance ERAP2 atoms
        entrance_erap2_idx = [i for i, a in enumerate(erap2_atoms)
                              if (a['resnum'] + OFFSET) in CHANNEL_ENTRANCE]
        erap2_entrance_proj = []
        for ei in entrance_erap2_idx:
            v = coords_e[ei] - entrance_center
            proj_along = np.dot(v, channel_axis)
            proj_on_plane = v - proj_along * channel_axis
            erap2_entrance_proj.append(proj_on_plane)

        # Compute approximate channel opening radius from ERAP2 entrance residues
        erap2_proj = np.array(erap2_entrance_proj)
        channel_radius = np.percentile(np.linalg.norm(erap2_proj, axis=1), 50)

        # How much of the channel cross-section does the binder cover?
        # Count binder atoms within the channel radius
        binder_in_channel = np.sum(np.linalg.norm(binder_proj, axis=1) < channel_radius)
        binder_near_entrance = np.sum(np.abs(np.dot(coords_b - entrance_center,
                                                      channel_axis.reshape(1, -1).T).flatten()) < 10)

        print(f"  Channel axis direction: [{channel_axis[0]:.3f}, {channel_axis[1]:.3f}, {channel_axis[2]:.3f}]")
        print(f"  Channel radius (median): {channel_radius:.1f} A")
        print(f"  Binder atoms within channel cross-section: {binder_in_channel}/{len(coords_b)}")
        print(f"  Binder atoms within 10A of entrance plane: {binder_near_entrance}/{len(coords_b)}")

        spatial_coverage = {
            'channel_radius_A': round(float(channel_radius), 1),
            'binder_atoms_in_crosssection': int(binder_in_channel),
            'binder_atoms_total': len(coords_b),
            'crosssection_fraction': round(float(binder_in_channel) / len(coords_b), 3),
            'binder_atoms_near_entrance': int(binder_near_entrance),
            'entrance_proximity_fraction': round(float(binder_near_entrance) / len(coords_b), 3)
        }
    else:
        spatial_coverage = {'error': 'Could not compute channel axis'}

    # === VERDICT ===
    print("\n[5/5] Computing verdict...")

    entrance_reduction = sasa_results['channel_entrance_350_419']['reduction_pct']
    zinc_reduction = sasa_results['zinc_coordination']['reduction_pct']
    bridges_both = len(dp1_contacts) > 0 and len(dp2_contacts) > 0
    contacts_zinc = len(zinc_contacts) > 0
    contacts_k392 = len(k392_contacts) > 0

    # Deep contacts suggest binder extends into channel, not just sits on surface
    deep_contacts = set(contact_erap2_full) & set(range(440, 500))
    extends_deep = len(deep_contacts) > 3

    if entrance_reduction >= 80:
        cap_verdict = "FULL_CAP"
        cap_explanation = "Binder reduces channel entrance SASA by >80%. Acts as a physical lid."
    elif entrance_reduction >= 50:
        cap_verdict = "PARTIAL_CAP"
        cap_explanation = f"Binder reduces entrance SASA by {entrance_reduction:.0f}%. Partial blockage — substrates may still enter through gaps."
    elif entrance_reduction >= 20 and bridges_both:
        cap_verdict = "BRIDGE_NOT_CAP"
        cap_explanation = f"Binder bridges both divergent patches ({entrance_reduction:.0f}% SASA reduction) but doesn't fully occlude the entrance. It spans across the channel mouth but leaves gaps."
    else:
        cap_verdict = "SURFACE_BINDER"
        cap_explanation = f"Binder sits on the surface ({entrance_reduction:.0f}% SASA reduction). Does not meaningfully occlude channel entrance."

    # Functional assessment
    if contacts_zinc and entrance_reduction >= 50:
        functional_verdict = "LIKELY_INHIBITORY"
        functional_explanation = "Binder contacts zinc coordination sphere AND reduces entrance SASA substantially. Should impede substrate access."
    elif contacts_k392 and entrance_reduction >= 30:
        functional_verdict = "POSSIBLY_INHIBITORY"
        functional_explanation = "Binder contacts K392 selectivity handle and partially blocks entrance. May reduce activity."
    else:
        functional_verdict = "UNLIKELY_INHIBITORY"
        functional_explanation = "Insufficient evidence that binder blocks substrate access to catalytic zinc."

    # Redesign recommendation
    if cap_verdict in ("FULL_CAP",):
        redesign_needed = False
        redesign_recommendation = "No redesign needed. Binder already caps channel entrance."
    elif cap_verdict == "PARTIAL_CAP":
        redesign_needed = True
        redesign_recommendation = (
            "Partial coverage. Recommend extending binder loops toward uncovered regions "
            f"(especially {sorted(CHANNEL_ENTRANCE - set(contact_erap2_full))[:10]}...). "
            "Could use Boltz-2 with extended hotspot residues or add loop insertions."
        )
    else:
        redesign_needed = True
        redesign_recommendation = (
            "Binder does not adequately cap channel. Options:\n"
            "  1. Reposition binder to center over channel opening (re-dock with entrance hotspots)\n"
            "  2. Extend binder with loop insertions that reach across uncovered gaps\n"
            "  3. Add a channel-plugging peptide tail to the binder (hybrid V2+V4 approach)\n"
            f"  Key gaps to fill: residues {sorted(set(range(388,400)) - set(contact_erap2_full))}"
        )

    results = {
        'analysis': 'V2 Protein Binder Channel-Capping Assessment',
        'binder': 'n248_trim_c5_Y87A_Y89A',
        'binder_size_aa': 92,
        'target': 'ERAP2 (cropped 350-500)',
        'timestamp': datetime.now(timezone.utc).isoformat(),
        'contact_analysis': {
            'total_erap2_contacts': len(contact_erap2_full),
            'erap2_contact_residues_full': contact_erap2_full,
            'divergent_patch_1_contacts': sorted(dp1_contacts),
            'divergent_patch_2_contacts': sorted(dp2_contacts),
            'zinc_coord_contacts': sorted(zinc_contacts),
            'k392_region_contacts': sorted(k392_contacts),
            'entrance_contacts_count': len(entrance_contacts),
            'deep_contacts_440_500': sorted(deep_contacts),
            'bridges_both_patches': bridges_both,
            'extends_deep_into_channel': extends_deep,
            'closest_approach_to_zinc': {
                str(r): round(erap2_res_min_dist.get(r, 99.9), 2) for r in sorted(ZINC_COORD)
            }
        },
        'sasa_analysis': sasa_results,
        'spatial_coverage': spatial_coverage,
        'verdicts': {
            'channel_capping': cap_verdict,
            'channel_capping_explanation': cap_explanation,
            'functional_inhibition': functional_verdict,
            'functional_explanation': functional_explanation,
            'redesign_needed': redesign_needed,
            'redesign_recommendation': redesign_recommendation
        }
    }

    # Write results
    output_file = os.path.join(OUTPUT_DIR, 'v2_channel_cap_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results written to: {output_file}")

    # Print verdict
    print(f"\n{'='*60}")
    print(f"CHANNEL CAPPING VERDICT: {cap_verdict}")
    print(f"  {cap_explanation}")
    print(f"\nFUNCTIONAL VERDICT: {functional_verdict}")
    print(f"  {functional_explanation}")
    print(f"\nREDESIGN NEEDED: {'YES' if redesign_needed else 'NO'}")
    if redesign_needed:
        print(f"  {redesign_recommendation}")
    print(f"{'='*60}")

    return results


if __name__ == '__main__':
    results = run_analysis()
