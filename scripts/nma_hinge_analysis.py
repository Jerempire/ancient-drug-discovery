"""
NMA Hinge Jammer Analysis: Does V2 binder jam ERAP2's open/closed transition?

Computes Anisotropic Network Model (ANM) for apo vs V2-bound ERAP2 and checks
whether the V2 binder suppresses the primary domain hinge mode.

KEY IMPROVEMENT: Superposes V2 binder onto FULL apo ERAP2 to build a complete
holo model (960 + 92 = 1052 CA atoms), enabling direct NMA comparison.

Uses SciPy directly (no ProDy dependency). ANM is a Kirchhoff/Hessian matrix
eigendecomposition on CA atoms with a distance cutoff.

Output: data/results/hinge_jammer/nma_comparison.json
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import eigh
import json
import os

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APO_PDB = os.path.join(PROJECT, "data/structures/erap2_wt_alphafold.pdb")
COMPLEX_CIF = os.path.join(PROJECT, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
OUTPUT_DIR = os.path.join(PROJECT, "data/results/hinge_jammer")
OUTPUT_JSON = os.path.join(OUTPUT_DIR, "nma_comparison.json")

# CIF residue 1 = apo residue 350; verified by sequence alignment (0 mismatches)
CIF_OFFSET = 349

# ERAP2 domain boundaries (UniProt Q6P179 numbering)
# Domain I:   1-254   (N-terminal beta-sandwich)
# Domain II:  255-528  (catalytic thermolysin-like, contains zinc site)
# Domain III: 529-614  (beta-sandwich)
# Domain IV:  615-960  (C-terminal armadillo-like, forms the "lid")
# The open->closed hinge is primarily between Domain II and Domain IV
DOMAIN_BOUNDARIES = {
    'I':   (1, 254),
    'II':  (255, 528),
    'III': (529, 614),
    'IV':  (615, 960),
}
# For hinge analysis: "body" = domains I+II, "lid" = domains III+IV
BODY_RANGE = (1, 528)
LID_RANGE = (529, 960)

ANM_CUTOFF = 15.0  # Angstroms, standard for ANM


def parse_pdb_ca(path):
    """Extract CA atoms from PDB file. Returns list of (resnum, x, y, z)."""
    cas = []
    with open(path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                cas.append((resnum, x, y, z))
    return cas


def parse_cif_ca(path, chain):
    """Extract CA atoms from mmCIF file for a given chain."""
    cas = []
    headers = []
    in_atom = False
    with open(path) as f:
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
                    atomname = r.get('auth_atom_id', r.get('label_atom_id', '?'))
                    if ch == chain and atomname == 'CA':
                        resnum = int(r.get('auth_seq_id', '0'))
                        x = float(r['Cartn_x'])
                        y = float(r['Cartn_y'])
                        z = float(r['Cartn_z'])
                        cas.append((resnum, x, y, z))
            elif in_atom and (line.startswith('#') or line.startswith('loop_')):
                if cas:
                    break
    return cas


def superpose_binder_onto_full_apo(apo_cas, cif_erap2_cas, cif_binder_cas):
    """Superpose CIF crop+binder onto full apo ERAP2 using Kabsch algorithm.

    Returns transformed binder coordinates and RMSD.
    """
    # Build lookup: apo resnum -> coords
    apo_dict = {c[0]: np.array([c[1], c[2], c[3]]) for c in apo_cas}

    # Match CIF ERAP2 CA atoms to apo (CIF resnum + CIF_OFFSET = apo resnum)
    matched_cif, matched_apo = [], []
    for c in cif_erap2_cas:
        apo_resnum = c[0] + CIF_OFFSET
        if apo_resnum in apo_dict:
            matched_cif.append(np.array([c[1], c[2], c[3]]))
            matched_apo.append(apo_dict[apo_resnum])

    P = np.array(matched_cif)
    Q = np.array(matched_apo)

    # Kabsch: find rotation R and translation t such that R*P + t ~ Q
    p_center = P.mean(axis=0)
    q_center = Q.mean(axis=0)
    P_centered = P - p_center
    Q_centered = Q - q_center

    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T

    # RMSD
    P_aligned = (P_centered @ R.T) + q_center
    rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q) ** 2, axis=1)))

    # Transform binder coordinates
    binder_coords_raw = np.array([[c[1], c[2], c[3]] for c in cif_binder_cas])
    binder_transformed = ((binder_coords_raw - p_center) @ R.T) + q_center

    return binder_transformed, rmsd, len(matched_cif)


def build_anm_hessian(coords, cutoff=ANM_CUTOFF):
    """Build the 3Nx3N ANM Hessian matrix."""
    n = len(coords)
    H = np.zeros((3 * n, 3 * n))
    dists = cdist(coords, coords)

    for i in range(n):
        for j in range(i + 1, n):
            if dists[i, j] < cutoff:
                diff = coords[j] - coords[i]
                d2 = dists[i, j] ** 2
                block = -np.outer(diff, diff) / d2  # gamma=1
                H[3*i:3*i+3, 3*j:3*j+3] = block
                H[3*j:3*j+3, 3*i:3*i+3] = block
                H[3*i:3*i+3, 3*i:3*i+3] -= block
                H[3*j:3*j+3, 3*j:3*j+3] -= block

    return H


def compute_anm_modes(coords, cutoff=ANM_CUTOFF, n_modes=20):
    """Compute lowest n_modes non-trivial ANM modes."""
    H = build_anm_hessian(coords, cutoff)
    eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, n_modes + 6 - 1])
    # Skip first 6 trivial modes (translation + rotation)
    return eigenvalues[6:], eigenvectors[:, 6:]


def compute_hinge_vector(coords, resnums, body_range, lid_range):
    """Compute reference hinge vector: body stays, lid moves."""
    body_mask = np.array([(body_range[0] <= r <= body_range[1]) for r in resnums])
    lid_mask = np.array([(lid_range[0] <= r <= lid_range[1]) for r in resnums])

    body_com = coords[body_mask].mean(axis=0)
    lid_com = coords[lid_mask].mean(axis=0)

    hinge_dir = lid_com - body_com
    hinge_dir /= np.linalg.norm(hinge_dir)

    n = len(coords)
    hinge_vec = np.zeros(3 * n)
    for i in range(n):
        if lid_mask[i]:
            hinge_vec[3*i:3*i+3] = hinge_dir

    norm = np.linalg.norm(hinge_vec)
    if norm > 0:
        hinge_vec /= norm
    return hinge_vec


def mode_overlap(mode_vec, ref_vec):
    """Overlap (cosine similarity) between a mode vector and a reference vector."""
    return abs(np.dot(mode_vec, ref_vec))


def compute_collectivity(mode_vec, n_atoms):
    """Compute collectivity of a mode (1/N=localized, 1.0=all atoms)."""
    disp = mode_vec.reshape(n_atoms, 3)
    sq_disp = np.sum(disp**2, axis=1)
    sq_disp /= sq_disp.sum()
    entropy = -np.sum(sq_disp * np.log(sq_disp + 1e-20))
    return np.exp(entropy) / n_atoms


def domain_mobility(mode_vec, resnums, domain_boundaries):
    """Compute per-domain mean squared displacement for a mode."""
    n = len(resnums)
    disp = mode_vec.reshape(n, 3)
    sq_disp = np.sum(disp**2, axis=1)

    result = {}
    for name, (start, end) in domain_boundaries.items():
        mask = np.array([(start <= r <= end) for r in resnums])
        if mask.any():
            result[name] = float(sq_disp[mask].mean())
    return result


def inter_domain_angle(coords, resnums, body_range, lid_range):
    """Compute angle between body and lid domain centroids relative to hinge."""
    body_mask = np.array([(body_range[0] <= r <= body_range[1]) for r in resnums])
    lid_mask = np.array([(lid_range[0] <= r <= lid_range[1]) for r in resnums])

    body_com = coords[body_mask].mean(axis=0)
    lid_com = coords[lid_mask].mean(axis=0)

    hinge_residues = [i for i, r in enumerate(resnums)
                      if body_range[1] - 10 <= r <= lid_range[0] + 10]
    if hinge_residues:
        hinge_point = coords[hinge_residues].mean(axis=0)
    else:
        hinge_point = (body_com + lid_com) / 2

    v1 = body_com - hinge_point
    v2 = lid_com - hinge_point
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
    return float(np.degrees(np.arccos(np.clip(cos_angle, -1, 1))))


def analyze_modes(coords, resnums, n_modes=20, label=""):
    """Run NMA and analyze modes. Returns eigenvalues and mode analysis list."""
    print(f"\n   Computing ANM ({len(coords)} CA atoms)...")
    eigenvalues, eigenvectors = compute_anm_modes(coords, n_modes=n_modes)
    print(f"   Eigenvalue range: {eigenvalues[0]:.6f} - {eigenvalues[-1]:.6f}")

    hinge_vec = compute_hinge_vector(coords, resnums, BODY_RANGE, LID_RANGE)

    mode_results = []
    hinge_mode_idx = None
    max_overlap = 0

    for i in range(min(10, n_modes)):
        mode = eigenvectors[:, i]
        ov = mode_overlap(mode, hinge_vec)
        coll = compute_collectivity(mode, len(coords))
        dom_mob = domain_mobility(mode, resnums, DOMAIN_BOUNDARIES)

        lid_mob = dom_mob.get('III', 0) + dom_mob.get('IV', 0)
        body_mob = dom_mob.get('I', 0) + dom_mob.get('II', 0)
        lid_ratio = lid_mob / (body_mob + 1e-10)

        assessment = ""
        if ov > 0.3 and lid_ratio > 1.5:
            assessment = "HINGE MODE"
            if ov > max_overlap:
                max_overlap = ov
                hinge_mode_idx = i
        elif ov > 0.2:
            assessment = "partial hinge"
        elif coll > 0.3:
            assessment = "collective"

        mode_results.append({
            'mode': i + 1,
            'eigenvalue': float(eigenvalues[i]),
            'frequency': float(np.sqrt(abs(eigenvalues[i]))),
            'hinge_overlap': float(ov),
            'collectivity': float(coll),
            'domain_mobility': dom_mob,
            'lid_body_ratio': float(lid_ratio),
            'assessment': assessment,
        })

    return eigenvalues, eigenvectors, mode_results, hinge_mode_idx, max_overlap


def save_ca_pdb(path, apo_coords, apo_resnums, binder_coords=None, binder_resnums=None):
    """Write CA-only PDB for visualization."""
    with open(path, 'w') as f:
        idx = 1
        for coord, resnum in zip(apo_coords, apo_resnums):
            f.write(f"ATOM  {idx:5d}  CA  ALA A{resnum:4d}    "
                    f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00\n")
            idx += 1
        if binder_coords is not None:
            f.write("TER\n")
            for coord, resnum in zip(binder_coords, binder_resnums):
                f.write(f"ATOM  {idx:5d}  CA  ALA B{resnum:4d}    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00\n")
                idx += 1
        f.write("TER\nEND\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 80)
    print("ERAP2 HINGE JAMMER ANALYSIS -- NMA")
    print("Does V2 binder suppress the open/closed domain hinge mode?")
    print("=" * 80)

    # --- Load structures ---
    print("\n1. Loading structures...")
    apo_cas = parse_pdb_ca(APO_PDB)
    cif_erap2_cas = parse_cif_ca(COMPLEX_CIF, 'A')
    cif_binder_cas = parse_cif_ca(COMPLEX_CIF, 'B')
    print(f"   Apo ERAP2: {len(apo_cas)} CA atoms (residues {apo_cas[0][0]}-{apo_cas[-1][0]})")
    print(f"   CIF ERAP2 crop: {len(cif_erap2_cas)} CA atoms (CIF res {cif_erap2_cas[0][0]}-{cif_erap2_cas[-1][0]})")
    print(f"   CIF V2 binder: {len(cif_binder_cas)} CA atoms")
    print(f"   CIF offset: +{CIF_OFFSET} (CIF res 1 = apo res {1 + CIF_OFFSET})")

    apo_resnums = np.array([c[0] for c in apo_cas])
    apo_coords = np.array([[c[1], c[2], c[3]] for c in apo_cas])

    # --- Superpose V2 onto full ERAP2 ---
    print("\n2. Superposing V2 binder onto full apo ERAP2...")
    binder_coords, sup_rmsd, n_matched = superpose_binder_onto_full_apo(
        apo_cas, cif_erap2_cas, cif_binder_cas)
    print(f"   Matched {n_matched} CA atoms, superposition RMSD: {sup_rmsd:.3f} A")

    binder_resnums = np.array([c[0] for c in cif_binder_cas])

    # Build holo: full apo ERAP2 + superposed V2 binder
    holo_coords = np.vstack([apo_coords, binder_coords])
    # Binder uses resnums 961+ to avoid collision with ERAP2
    holo_resnums = np.concatenate([apo_resnums, binder_resnums + 960])
    n_erap2 = len(apo_coords)
    n_binder = len(binder_coords)
    print(f"   Holo complex: {len(holo_coords)} CA atoms ({n_erap2} ERAP2 + {n_binder} binder)")

    # Save complex PDB
    complex_pdb = os.path.join(OUTPUT_DIR, "full_erap2_v2_complex_ca.pdb")
    save_ca_pdb(complex_pdb, apo_coords, apo_resnums, binder_coords, binder_resnums)
    print(f"   Saved: {complex_pdb}")

    # --- Binder contact analysis ---
    print("\n3. Binder contact analysis (hinge region)...")
    cross_dists = cdist(binder_coords, apo_coords)
    min_dist_per_erap2 = cross_dists.min(axis=0)

    binder_contacts = {}
    for i, resnum in enumerate(apo_resnums):
        d = min_dist_per_erap2[i]
        if d < ANM_CUTOFF:
            domain = "?"
            for name, (s, e) in DOMAIN_BOUNDARIES.items():
                if s <= resnum <= e:
                    domain = name
                    break
            binder_contacts[int(resnum)] = {
                'min_dist': float(d),
                'domain': domain,
                'is_hinge_region': bool(518 <= resnum <= 540)
            }

    close_contacts = {k: v for k, v in binder_contacts.items() if v['min_dist'] < 8.0}
    hinge_contacts = {k: v for k, v in close_contacts.items() if v['is_hinge_region']}
    print(f"   ERAP2 residues within {ANM_CUTOFF}A of binder: {len(binder_contacts)}")
    print(f"   Within 8A (close contacts): {len(close_contacts)}")
    print(f"   Hinge-region contacts (518-540): {len(hinge_contacts)}")

    top_close = sorted(close_contacts.items(), key=lambda x: x[1]['min_dist'])[:15]
    for resnum, info in top_close:
        tag = " <-- HINGE" if info['is_hinge_region'] else ""
        print(f"     Res {resnum} (D{info['domain']}): {info['min_dist']:.1f}A{tag}")

    # --- APO NMA ---
    print("\n4. APO NMA (full-length ERAP2)...")
    apo_evals, apo_evecs, apo_modes, apo_hinge_idx, apo_max_ov = analyze_modes(
        apo_coords, apo_resnums, label="APO")

    print(f"\n   {'Mode':>5} {'Freq':>10} {'Overlap':>10} {'Collect':>10} {'Lid/Body':>10} {'Assessment':>15}")
    print("   " + "-" * 65)
    for mr in apo_modes:
        print(f"   {mr['mode']:>5} {mr['frequency']:>10.6f} {mr['hinge_overlap']:>10.4f} "
              f"{mr['collectivity']:>10.4f} {mr['lid_body_ratio']:>10.2f} {mr['assessment']:>15}")

    if apo_hinge_idx is not None:
        print(f"\n   APO primary hinge: Mode {apo_hinge_idx + 1} (overlap {apo_max_ov:.4f})")
    else:
        print("\n   No single dominant hinge mode in apo. Motion may be distributed.")

    # --- HOLO NMA ---
    print("\n5. HOLO NMA (full ERAP2 + V2 binder)...")

    # For holo, the hinge vector should only reference ERAP2 residues
    # We compute NMA on full holo, but extract ERAP2-only displacements for overlap
    print(f"   Building Hessian for {len(holo_coords)} atoms...")
    holo_evals, holo_evecs = compute_anm_modes(holo_coords, n_modes=20)

    # Compute hinge overlap using ERAP2-only portion of eigenvectors
    hinge_vec_erap2 = compute_hinge_vector(apo_coords, apo_resnums, BODY_RANGE, LID_RANGE)
    # Extend to full holo size (binder portion = zero, doesn't contribute to hinge overlap)
    hinge_vec_holo = np.zeros(3 * len(holo_coords))
    hinge_vec_holo[:3 * n_erap2] = hinge_vec_erap2
    hinge_vec_holo /= np.linalg.norm(hinge_vec_holo)

    holo_modes = []
    holo_hinge_idx = None
    holo_max_ov = 0

    for i in range(min(10, 20)):
        mode = holo_evecs[:, i]
        ov = mode_overlap(mode, hinge_vec_holo)
        coll = compute_collectivity(mode, len(holo_coords))

        # Domain mobility using ERAP2 portion only
        erap2_mode = mode[:3 * n_erap2]
        dom_mob = domain_mobility(erap2_mode, apo_resnums, DOMAIN_BOUNDARIES)

        lid_mob = dom_mob.get('III', 0) + dom_mob.get('IV', 0)
        body_mob = dom_mob.get('I', 0) + dom_mob.get('II', 0)
        lid_ratio = lid_mob / (body_mob + 1e-10)

        assessment = ""
        if ov > 0.3 and lid_ratio > 1.5:
            assessment = "HINGE MODE"
            if ov > holo_max_ov:
                holo_max_ov = ov
                holo_hinge_idx = i
        elif ov > 0.2:
            assessment = "partial hinge"
        elif coll > 0.3:
            assessment = "collective"

        holo_modes.append({
            'mode': i + 1,
            'eigenvalue': float(holo_evals[i]),
            'frequency': float(np.sqrt(abs(holo_evals[i]))),
            'hinge_overlap': float(ov),
            'collectivity': float(coll),
            'domain_mobility': dom_mob,
            'lid_body_ratio': float(lid_ratio),
            'assessment': assessment,
        })

    print(f"\n   {'Mode':>5} {'Freq':>10} {'Overlap':>10} {'Collect':>10} {'Lid/Body':>10} {'Assessment':>15}")
    print("   " + "-" * 65)
    for mr in holo_modes:
        print(f"   {mr['mode']:>5} {mr['frequency']:>10.6f} {mr['hinge_overlap']:>10.4f} "
              f"{mr['collectivity']:>10.4f} {mr['lid_body_ratio']:>10.2f} {mr['assessment']:>15}")

    if holo_hinge_idx is not None:
        print(f"\n   HOLO primary hinge: Mode {holo_hinge_idx + 1} (overlap {holo_max_ov:.4f})")
    else:
        print("\n   No single dominant hinge mode in holo.")

    # --- Direct comparison ---
    print("\n6. DIRECT COMPARISON: Apo vs Holo hinge modes")
    print("-" * 70)

    if apo_hinge_idx is not None and holo_hinge_idx is not None:
        apo_hinge = apo_modes[apo_hinge_idx]
        holo_hinge = holo_modes[holo_hinge_idx]

        freq_ratio = holo_hinge['frequency'] / apo_hinge['frequency']
        eval_ratio = holo_hinge['eigenvalue'] / apo_hinge['eigenvalue']

        print(f"   Apo hinge:  Mode {apo_hinge['mode']}, freq={apo_hinge['frequency']:.6f}, "
              f"overlap={apo_hinge['hinge_overlap']:.4f}")
        print(f"   Holo hinge: Mode {holo_hinge['mode']}, freq={holo_hinge['frequency']:.6f}, "
              f"overlap={holo_hinge['hinge_overlap']:.4f}")
        print(f"\n   Frequency ratio (holo/apo): {freq_ratio:.4f}")
        print(f"   Eigenvalue ratio (holo/apo): {eval_ratio:.4f}")

        if freq_ratio > 1.5:
            freq_verdict = "STRONG suppression (>1.5x stiffer)"
        elif freq_ratio > 1.2:
            freq_verdict = "MODERATE suppression (1.2-1.5x stiffer)"
        elif freq_ratio > 1.05:
            freq_verdict = "WEAK suppression (1.05-1.2x stiffer)"
        elif freq_ratio > 0.95:
            freq_verdict = "NO CHANGE (~1.0x)"
        else:
            freq_verdict = "SOFTENING (<0.95x, binder makes hinge MORE flexible)"

        print(f"   Frequency verdict: {freq_verdict}")

        # Overlap change (higher in holo = hinge character redistributed)
        ov_change = holo_hinge['hinge_overlap'] - apo_hinge['hinge_overlap']
        print(f"   Overlap change: {ov_change:+.4f} ({'more' if ov_change > 0 else 'less'} hinge-like)")

    else:
        freq_ratio = None
        eval_ratio = None
        freq_verdict = "CANNOT COMPARE (hinge mode not identified in one or both)"
        print(f"   {freq_verdict}")

    # Also compare first 5 modes side-by-side
    print(f"\n   Mode-by-mode comparison (first 10):")
    print(f"   {'Mode':>5} | {'Apo freq':>10} {'Holo freq':>10} {'Ratio':>8} | "
          f"{'Apo ov':>8} {'Holo ov':>8}")
    print("   " + "-" * 65)
    for i in range(min(10, len(apo_modes), len(holo_modes))):
        ratio = holo_modes[i]['frequency'] / apo_modes[i]['frequency'] if apo_modes[i]['frequency'] > 0 else 0
        print(f"   {i+1:>5} | {apo_modes[i]['frequency']:>10.6f} {holo_modes[i]['frequency']:>10.6f} "
              f"{ratio:>8.4f} | {apo_modes[i]['hinge_overlap']:>8.4f} {holo_modes[i]['hinge_overlap']:>8.4f}")

    # --- Inter-domain angle ---
    print("\n7. Inter-domain angle...")
    apo_angle = inter_domain_angle(apo_coords, apo_resnums, BODY_RANGE, LID_RANGE)
    print(f"   Apo body-hinge-lid angle: {apo_angle:.1f} degrees")
    print("   Note: Static angle is identical in holo (built by superposition).")
    print("   The NMA frequency shift IS the key metric for hinge suppression.")

    # --- KILL GATE ---
    print("\n" + "=" * 80)
    print("KILL GATE ASSESSMENT")
    print("=" * 80)

    score = 0
    reasons = []

    # 1. Frequency ratio
    if freq_ratio is not None:
        if freq_ratio > 1.5:
            score += 40
            reasons.append(f"STRONG: Hinge mode {freq_ratio:.2f}x stiffer in holo")
        elif freq_ratio > 1.2:
            score += 25
            reasons.append(f"MODERATE: Hinge mode {freq_ratio:.2f}x stiffer in holo")
        elif freq_ratio > 1.05:
            score += 10
            reasons.append(f"WEAK: Hinge mode only {freq_ratio:.2f}x stiffer")
        else:
            reasons.append(f"NO SUPPRESSION: Ratio {freq_ratio:.2f}x")
    else:
        score -= 10
        reasons.append("Could not identify clear hinge mode for comparison")

    # 2. Hinge mode character
    if apo_hinge_idx is not None:
        score += 15
        reasons.append(f"Clear apo hinge mode (Mode {apo_hinge_idx+1}, overlap {apo_max_ov:.3f})")
    else:
        reasons.append("No dominant apo hinge mode (motion distributed)")

    # 3. Binder contacts at hinge
    if hinge_contacts:
        score += 20
        reasons.append(f"V2 directly contacts {len(hinge_contacts)} hinge-region residues (518-540)")
    elif close_contacts:
        near_hinge = {k: v for k, v in close_contacts.items() if 490 <= k <= 530}
        if near_hinge:
            score += 10
            reasons.append(f"V2 contacts {len(near_hinge)} residues near hinge (490-530)")
        else:
            reasons.append("V2 does not contact hinge region directly")

    # 4. Network effect (binder adds cross-domain contacts)
    n_body_contacts = sum(1 for k, v in binder_contacts.items()
                         if BODY_RANGE[0] <= k <= BODY_RANGE[1] and v['min_dist'] < 8.0)
    n_lid_contacts = sum(1 for k, v in binder_contacts.items()
                        if LID_RANGE[0] <= k <= LID_RANGE[1] and v['min_dist'] < 8.0)
    if n_body_contacts > 0 and n_lid_contacts > 0:
        score += 15
        reasons.append(f"V2 bridges body ({n_body_contacts} contacts) and lid ({n_lid_contacts} contacts)")
    elif n_body_contacts > 0 or n_lid_contacts > 0:
        score += 5
        which = "body" if n_body_contacts > 0 else "lid"
        reasons.append(f"V2 contacts {which} domain only — no cross-domain bridging")

    # Verdict
    print(f"\n   Hinge Suppression Score: {score}/100")
    for r in reasons:
        print(f"   - {r}")

    if score >= 50:
        verdict = "PASS -- V2 likely jams the hinge. Allosteric inhibition is plausible."
        action = "Priority shift: present to user, consider allosteric paper. Cyclic peptide track becomes independent."
    elif score >= 25:
        verdict = "INCONCLUSIVE -- Some evidence but ambiguous. Run targeted MD to confirm."
        action = "Run 10ns MD on full ERAP2+V2 complex (Vast.ai, ~$0.50) to measure hinge angle dynamics."
    else:
        verdict = "FAIL -- V2 does not significantly suppress the hinge mode."
        action = "Proceed full-speed on Track 2 (cyclic peptides) and Track 3 (HLA-A29)."

    print(f"\n   VERDICT: {verdict}")
    print(f"   ACTION:  {action}")

    # --- Save results ---
    results = {
        'analysis': 'NMA Hinge Jammer (v2 - superposition method)',
        'structures': {
            'apo': APO_PDB,
            'complex_cif': COMPLEX_CIF,
            'complex_pdb_built': complex_pdb,
            'superposition_rmsd': float(sup_rmsd),
            'n_matched_atoms': n_matched,
            'cif_offset': CIF_OFFSET,
        },
        'apo_nma': {
            'n_atoms': n_erap2,
            'modes': apo_modes,
            'primary_hinge_mode': apo_hinge_idx + 1 if apo_hinge_idx is not None else None,
            'max_hinge_overlap': float(apo_max_ov),
        },
        'holo_nma': {
            'n_atoms': len(holo_coords),
            'n_erap2': n_erap2,
            'n_binder': n_binder,
            'modes': holo_modes,
            'primary_hinge_mode': holo_hinge_idx + 1 if holo_hinge_idx is not None else None,
            'max_hinge_overlap': float(holo_max_ov),
        },
        'comparison': {
            'frequency_ratio': float(freq_ratio) if freq_ratio else None,
            'eigenvalue_ratio': float(eval_ratio) if eval_ratio else None,
            'frequency_verdict': freq_verdict,
        },
        'contact_analysis': {
            'total_within_cutoff': len(binder_contacts),
            'close_contacts_8A': len(close_contacts),
            'hinge_region_contacts': len(hinge_contacts),
            'body_contacts': n_body_contacts,
            'lid_contacts': n_lid_contacts,
            'top_contacts': {str(k): v for k, v in sorted(close_contacts.items(),
                                                           key=lambda x: x[1]['min_dist'])[:20]}
        },
        'inter_domain_angle_deg': float(apo_angle),
        'kill_gate': {
            'score': score,
            'threshold_pass': 50,
            'threshold_inconclusive': 25,
            'verdict': verdict,
            'action': action,
            'reasons': reasons,
        }
    }

    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n   Results saved: {OUTPUT_JSON}")

    # Also save hinge geometry separately
    hinge_json = os.path.join(OUTPUT_DIR, "hinge_geometry.json")
    with open(hinge_json, 'w') as f:
        json.dump({
            'apo_angle_deg': float(apo_angle),
            'domain_boundaries': DOMAIN_BOUNDARIES,
            'body_range': list(BODY_RANGE),
            'lid_range': list(LID_RANGE),
            'binder_contact_residues': {str(k): v for k, v in close_contacts.items()},
        }, f, indent=2)
    print(f"   Hinge geometry saved: {hinge_json}")

    return results


if __name__ == '__main__':
    main()
