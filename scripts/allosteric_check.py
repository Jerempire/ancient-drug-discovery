"""
Allosteric Check: Does V2 binder binding perturb ERAP2's zinc coordination?

Compares apo ERAP2 (AlphaFold) vs V2-bound ERAP2 (Boltz-2 complex).
After Kabsch superposition, measures shifts in the catalytic triad (H370, H374, E393).

Output: data/results/v43_validation/terminal_e/allosteric_check.json
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import numpy as np

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APO_PDB = os.path.join(PROJECT, "data/structures/erap2_wt_alphafold.pdb")
HOLO_CIF = os.path.join(PROJECT, "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif")
OUTPUT = os.path.join(PROJECT, "data/results/v43_validation/terminal_e/allosteric_check.json")
OFFSET = 349


def parse_pdb_atoms(pdb_file, chain='A'):
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('ATOM') and line[21] == chain:
                atoms.append({
                    'resnum': int(line[22:26]),
                    'resname': line[17:20].strip(),
                    'atomname': line[12:16].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                })
    return atoms


def parse_cif_atoms(cif_file, chain='A'):
    atoms = []
    headers = []
    with open(cif_file) as f:
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
                    ch = record.get('auth_asym_id', record.get('label_asym_id', '?'))
                    if ch == chain:
                        atoms.append({
                            'resnum': int(record.get('auth_seq_id', '0')),
                            'resname': record.get('label_comp_id', 'UNK'),
                            'atomname': record.get('auth_atom_id', record.get('label_atom_id', '?')),
                            'x': float(record['Cartn_x']),
                            'y': float(record['Cartn_y']),
                            'z': float(record['Cartn_z']),
                        })
            elif in_atom and (line.startswith('#') or line.startswith('loop_')):
                if atoms:
                    break
    return atoms


def get_atom_coord(atoms, resnum, atomname):
    for a in atoms:
        if a['resnum'] == resnum and a['atomname'] == atomname:
            return np.array([a['x'], a['y'], a['z']])
    return None


def get_ca(atoms, resnum):
    return get_atom_coord(atoms, resnum, 'CA')


def kabsch_align(P, Q):
    p_center = P.mean(axis=0)
    q_center = Q.mean(axis=0)
    P_c = P - p_center
    Q_c = Q - q_center
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T
    P_aligned = (P_c @ R.T) + q_center
    rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q) ** 2, axis=1)))
    return P_aligned, rmsd


def main():
    print("=" * 70)
    print("ALLOSTERIC CHECK: Zinc Coordination — Apo vs V2-Bound ERAP2")
    print("=" * 70)

    apo_atoms = parse_pdb_atoms(APO_PDB, chain='A')
    holo_atoms_raw = parse_cif_atoms(HOLO_CIF, chain='A')
    for a in holo_atoms_raw:
        a['resnum'] += OFFSET

    print(f"Apo: {len(apo_atoms)} atoms | Holo: {len(holo_atoms_raw)} atoms (offset +{OFFSET})")

    # Collect CA atoms for alignment (residues 350-500)
    common = []
    apo_ca = []
    holo_ca = []
    for r in range(350, 501):
        ac = get_ca(apo_atoms, r)
        hc = get_ca(holo_atoms_raw, r)
        if ac is not None and hc is not None:
            common.append(r)
            apo_ca.append(ac)
            holo_ca.append(hc)

    apo_ca = np.array(apo_ca)
    holo_ca = np.array(holo_ca)
    apo_aligned, global_rmsd = kabsch_align(apo_ca, holo_ca)

    print(f"Global CA RMSD ({len(common)} residues): {global_rmsd:.2f} A")

    apo_dict = {r: apo_aligned[i] for i, r in enumerate(common)}
    holo_dict = {r: holo_ca[i] for i, r in enumerate(common)}

    # === 1. Per-residue shifts in catalytic zone ===
    print("\n--- Per-Residue CA Shifts (368-395) ---")
    roles = {
        370: "H370 ZINC COORD", 374: "H374 ZINC COORD", 393: "E393 ZINC COORD",
        392: "K392 SELECTIVITY", 398: "Y398",
        371: "channel", 372: "channel", 373: "channel", 375: "HEXXH",
        376: "channel", 377: "channel", 378: "channel", 379: "channel",
        380: "channel", 381: "channel", 382: "channel", 383: "channel",
        384: "channel", 385: "channel", 386: "channel",
        387: "near K392", 388: "near K392", 389: "near K392",
        390: "near K392", 391: "adj K392", 394: "post-zinc", 395: "post-zinc",
        368: "entrance", 369: "entrance",
    }

    per_residue = {}
    print(f"{'Res':>6} {'Shift':>8} {'Role':>20} {'Signal':>10}")
    print("-" * 50)
    for r in range(368, 400):
        if r in apo_dict and r in holo_dict:
            shift = float(np.linalg.norm(apo_dict[r] - holo_dict[r]))
            role = roles.get(r, "")
            sig = "***" if shift > 1.5 else "**" if shift > 1.0 else "*" if shift > 0.5 else ""
            print(f"{r:>6} {shift:>7.2f}A {role:>20} {sig:>10}")
            per_residue[r] = {"shift_A": round(shift, 3), "role": role}

    # === 2. Zinc triad geometry ===
    print("\n--- Zinc Triad Geometry ---")
    triad = [370, 374, 393]
    triad_shifts = []
    for r in triad:
        if r in apo_dict and r in holo_dict:
            s = float(np.linalg.norm(apo_dict[r] - holo_dict[r]))
            triad_shifts.append(s)
            print(f"  {r} CA shift: {s:.2f} A")

    pairs = [(370, 374), (370, 393), (374, 393)]
    pair_results = {}
    print()
    print(f"{'Pair':>15} {'Apo':>8} {'Holo':>8} {'Change':>8}")
    print("-" * 45)
    for r1, r2 in pairs:
        if r1 in apo_dict and r2 in apo_dict and r1 in holo_dict and r2 in holo_dict:
            d_apo = float(np.linalg.norm(apo_dict[r1] - apo_dict[r2]))
            d_holo = float(np.linalg.norm(holo_dict[r1] - holo_dict[r2]))
            change = d_holo - d_apo
            label = f"H{r1}-{'H' if r2 != 393 else 'E'}{r2}"
            print(f"{label:>15} {d_apo:>7.2f}A {d_holo:>7.2f}A {change:>+7.2f}A")
            pair_results[label] = {
                "apo_A": round(d_apo, 2),
                "holo_A": round(d_holo, 2),
                "change_A": round(change, 2),
            }

    # === 3. Triangle area (zinc triad planarity) ===
    if all(r in apo_dict for r in triad) and all(r in holo_dict for r in triad):
        def triangle_area(p1, p2, p3):
            return 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))

        area_apo = triangle_area(apo_dict[370], apo_dict[374], apo_dict[393])
        area_holo = triangle_area(holo_dict[370], holo_dict[374], holo_dict[393])
        area_change = ((area_holo - area_apo) / area_apo) * 100
        print(f"\n  Triad triangle area: apo={area_apo:.1f} A^2, holo={area_holo:.1f} A^2, change={area_change:+.1f}%")

    # === 4. Regional RMSD comparison ===
    print("\n--- Regional CA RMSD ---")
    regions = {
        "Zinc triad (370,374,393)": [370, 374, 393],
        "Catalytic zone (368-395)": list(range(368, 396)),
        "Divergent patch 1 (353-367)": list(range(353, 368)),
        "Divergent patch 2 (400-414)": list(range(400, 415)),
        "Channel ceiling (409-419)": list(range(409, 420)),
        "Deep channel (440-490)": list(range(440, 491)),
        "Far from binder (460-500)": list(range(460, 501)),
    }

    region_rmsds = {}
    print(f"{'Region':>35} {'N':>4} {'RMSD':>8} {'Max':>8} {'Assessment':>15}")
    print("-" * 75)
    for label, residues in regions.items():
        matched = [r for r in residues if r in apo_dict and r in holo_dict]
        if len(matched) < 2:
            continue
        apo_r = np.array([apo_dict[r] for r in matched])
        holo_r = np.array([holo_dict[r] for r in matched])
        diffs = np.linalg.norm(apo_r - holo_r, axis=1)
        rmsd = float(np.sqrt(np.mean(diffs ** 2)))
        mx = float(np.max(diffs))
        assess = "RIGID" if rmsd < 0.5 else "MINOR" if rmsd < 1.0 else "NOTABLE" if rmsd < 2.0 else "MAJOR"
        print(f"{label:>35} {len(matched):>4} {rmsd:>7.2f}A {mx:>7.2f}A {assess:>15}")
        region_rmsds[label] = {"n": len(matched), "rmsd": round(rmsd, 3), "max_shift": round(mx, 3)}

    # === VERDICT ===
    avg_zinc = np.mean(triad_shifts) if triad_shifts else 0
    max_zinc = np.max(triad_shifts) if triad_shifts else 0

    if avg_zinc > 1.5:
        verdict = "STRONG_ALLOSTERIC_SIGNAL"
        explanation = f"Avg zinc triad shift = {avg_zinc:.2f} A. V2 binding significantly distorts catalytic geometry. Strong candidate for allosteric inhibition."
    elif avg_zinc > 0.75:
        verdict = "MODERATE_ALLOSTERIC_SIGNAL"
        explanation = f"Avg zinc triad shift = {avg_zinc:.2f} A. Detectable perturbation — worth investigating with MD. May reduce catalytic efficiency."
    elif avg_zinc > 0.3:
        verdict = "WEAK_ALLOSTERIC_SIGNAL"
        explanation = f"Avg zinc triad shift = {avg_zinc:.2f} A. Within prediction noise for Boltz-2 vs AlphaFold comparison. Unlikely functional."
    else:
        verdict = "NO_ALLOSTERIC_SIGNAL"
        explanation = f"Avg zinc triad shift = {avg_zinc:.2f} A. Zinc coordination is rigid despite V2 binding."

    print()
    print("=" * 70)
    print(f"VERDICT: {verdict}")
    print(f"  {explanation}")
    print("=" * 70)

    # Save
    output = {
        "analysis": "Allosteric check: apo vs V2-bound ERAP2 zinc coordination",
        "structures": {"apo": "AlphaFold erap2_wt", "holo": "Boltz-2 Y87A_Y89A complex"},
        "global_ca_rmsd": round(global_rmsd, 3),
        "n_aligned_residues": len(common),
        "zinc_triad_shifts": {str(r): round(s, 3) for r, s in zip(triad, triad_shifts)},
        "zinc_triad_avg_shift": round(avg_zinc, 3),
        "zinc_triad_max_shift": round(max_zinc, 3),
        "triad_distances": pair_results,
        "triad_triangle_area": {
            "apo_A2": round(float(area_apo), 1),
            "holo_A2": round(float(area_holo), 1),
            "change_pct": round(float(area_change), 1),
        },
        "per_residue_shifts": per_residue,
        "regional_rmsds": region_rmsds,
        "verdict": verdict,
        "explanation": explanation,
    }
    with open(OUTPUT, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
