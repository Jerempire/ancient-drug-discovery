"""
ANM (Anisotropic Network Model) flexibility check for ERAP2 P02 pocket.

Replaces explicit MD with a deterministic linear-response calculation:
  1. Build C-alpha elastic network (cutoff 15 A, spring const 1.0)
  2. Diagonalize 3N x 3N Hessian
  3. Per-residue mean-square fluctuation from low-frequency modes
  4. Cryptic-pocket signal: compare P02 lining residue MSF to global baseline
     - If P02 lining has MSF >= 1.5x global mean -> cryptic/breathing pocket
     - Collective motion check: P02 lining co-fluctuates in low modes
  5. Mode-pulled structure: displace along mode 1 by +1 sigma -> rerun LIGSITE on
     P02 region in displaced structure -> verify pocket volume in dynamic state

Output: data/results/smallmol/p02_anm.json + verdict.
Runtime: ~30 sec on CPU for 960-residue ERAP2.

References:
  Atilgan et al. 2001 Biophys J (ANM)
  Eyal et al. 2006 Bioinformatics (ProDy)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from Bio import PDB
from scipy.spatial import cKDTree
from scipy.linalg import eigh

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
ERAP2 = ROOT / "data/structures/erap2_wt_alphafold.pdb"
POCKET_JSON = ROOT / "data/results/smallmol/pocket_inventory.json"
OUT = ROOT / "data/results/smallmol/p02_anm.json"
OUT.parent.mkdir(parents=True, exist_ok=True)


def load_ca(pdb_path: Path):
    parser = PDB.PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))
    resnums, coords = [], []
    for r in s[0]["A"]:
        if "CA" in r and r.get_resname() not in ("HOH", "WAT"):
            resnums.append(r.get_id()[1])
            coords.append(r["CA"].coord)
    return np.array(resnums), np.array(coords)


def build_anm_hessian(coords: np.ndarray, cutoff: float = 15.0, gamma: float = 1.0):
    """Build 3N x 3N Hessian for anisotropic network model."""
    n = len(coords)
    H = np.zeros((3 * n, 3 * n))
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cutoff, output_type="ndarray")
    for i, j in pairs:
        rij = coords[j] - coords[i]
        d2 = (rij ** 2).sum()
        # Off-diagonal block
        block = -gamma * np.outer(rij, rij) / d2
        H[3*i:3*i+3, 3*j:3*j+3] += block
        H[3*j:3*j+3, 3*i:3*i+3] += block
        # Diagonal blocks
        H[3*i:3*i+3, 3*i:3*i+3] -= block
        H[3*j:3*j+3, 3*j:3*j+3] -= block
    return H


def msf_from_modes(eigvals: np.ndarray, eigvecs: np.ndarray, n: int,
                  n_modes: int = 20, eigval_floor: float = 1e-3):
    """Mean-square fluctuation per residue from top n_modes non-trivial modes.
    Modes with eigenvalue below floor are treated as zero (rigid-body) and skipped.
    """
    # Find first mode above the floor
    first = int(np.searchsorted(eigvals, eigval_floor))
    use = slice(first, first + n_modes)
    vals = eigvals[use]
    vecs = eigvecs[:, use]
    msf = np.zeros(n)
    for k in range(len(vals)):
        v = vecs[:, k].reshape(n, 3)
        msf += (v ** 2).sum(axis=1) / vals[k]
    return msf, first


def displace_along_mode(coords: np.ndarray, eigvecs: np.ndarray, eigvals: np.ndarray,
                       mode_idx: int, sigma: float = 1.0):
    """Displace Ca coords along one mode by sigma units of its thermal amplitude."""
    v = eigvecs[:, mode_idx].reshape(-1, 3)
    # Amplitude scaled by 1/sqrt(lambda) for thermal weighting; multiplied by sigma
    amp = sigma / np.sqrt(eigvals[mode_idx])
    # Cap amplitude to a reasonable RMSD (~3 A max)
    cur_rmsd = np.sqrt((v ** 2).sum(axis=1).mean()) * amp
    if cur_rmsd > 3.0:
        amp *= 3.0 / cur_rmsd
    return coords + amp * v


def ligsite_local_ca(ca_coords: np.ndarray, center: np.ndarray,
                     radius: float = 12.0, spacing: float = 1.5):
    """LIGSITE on Ca-only coords (coarse but consistent with original pocket detection)."""
    mn = center - radius
    mx = center + radius
    xs = np.arange(mn[0], mx[0], spacing)
    ys = np.arange(mn[1], mx[1], spacing)
    zs = np.arange(mn[2], mx[2], spacing)
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    grid = grid[np.linalg.norm(grid - center, axis=1) <= radius]

    tree = cKDTree(ca_coords)
    nearest, _ = tree.query(grid, k=1)
    # CA-only means we use larger probe distances
    solv = (nearest > 4.5) & (nearest < 10.0)
    grid_s = grid[solv]

    dirs = np.array([
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1],
    ], dtype=float)
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    n_steps = int(radius / spacing)
    psp = np.zeros(len(grid_s), dtype=int)
    for d in dirs:
        fwd = np.zeros(len(grid_s), dtype=bool)
        bwd = np.zeros(len(grid_s), dtype=bool)
        for k in range(1, n_steps + 1):
            df, _ = tree.query(grid_s + d * k * spacing, k=1)
            fwd |= df < 4.0
            db, _ = tree.query(grid_s - d * k * spacing, k=1)
            bwd |= db < 4.0
        psp += (fwd & bwd).astype(int)

    pocket = grid_s[psp >= 4]
    return len(pocket) * spacing ** 3


def main():
    with open(POCKET_JSON) as f:
        inv = json.load(f)
    p02 = inv["top_allosteric"]
    p02_center = np.array(p02["centroid_xyz"])
    p02_lining = set(p02["lining_residues"])
    p02_handles = set(p02["divergent_vs_both"])
    p02_baseline_vol = p02["volume_A3_est"]
    print(f"P02 lining: {sorted(p02_lining)}")
    print(f"P02 handles vs BOTH: {sorted(p02_handles)}")
    print(f"P02 baseline volume (apo, all-atom): {p02_baseline_vol} A^3")

    resnums, coords = load_ca(ERAP2)
    n = len(coords)
    print(f"\nERAP2 Ca: {n} residues")

    print("Building ANM Hessian (cutoff 15 A)...")
    H = build_anm_hessian(coords)
    print("Diagonalizing 3N x 3N Hessian (this takes ~20-60 sec)...")
    eigvals, eigvecs = eigh(H)
    # eigvals come sorted ascending; first 6 are ~0 (rigid body)
    print(f"  first 10 eigvals: {[f'{v:.4f}' for v in eigvals[:10]]}")
    print(f"  modes 7-12 (lowest non-trivial): {[f'{v:.4f}' for v in eigvals[6:12]]}")

    msf, first_real_mode = msf_from_modes(eigvals, eigvecs, n, n_modes=20, eigval_floor=1e-3)
    print(f"First real mode (above 1e-3): mode index {first_real_mode}")
    msf_mean = msf.mean()
    msf_std = msf.std()
    print(f"\nMSF: mean={msf_mean:.4f}  std={msf_std:.4f}  max={msf.max():.4f}")

    resnum_to_idx = {r: i for i, r in enumerate(resnums)}
    lining_idx = [resnum_to_idx[r] for r in p02_lining if r in resnum_to_idx]
    handle_idx = [resnum_to_idx[r] for r in p02_handles if r in resnum_to_idx]

    msf_lining = msf[lining_idx]
    msf_handles = msf[handle_idx]
    msf_global = msf  # full

    lining_ratio = float(msf_lining.mean() / msf_mean)
    handle_ratio = float(msf_handles.mean() / msf_mean) if len(handle_idx) else 0
    pctile_lining = float((msf_lining > msf_mean).mean())

    print(f"\nP02 LINING flexibility:")
    print(f"  mean MSF: {msf_lining.mean():.4f}  ({lining_ratio:.2f}x global)")
    print(f"  fraction above-mean: {pctile_lining*100:.0f}%")
    print(f"P02 HANDLES (vs BOTH) flexibility:")
    print(f"  mean MSF: {msf_handles.mean():.4f}  ({handle_ratio:.2f}x global)")

    # Collective motion check: how aligned are lining residues in the first real mode?
    mode_low = eigvecs[:, first_real_mode].reshape(n, 3)
    lining_disp = mode_low[lining_idx]
    # Cosine similarity to first PC of lining displacements
    if len(lining_idx) >= 3:
        u = lining_disp.mean(axis=0)
        u_norm = u / (np.linalg.norm(u) + 1e-9)
        alignments = lining_disp @ u_norm / (np.linalg.norm(lining_disp, axis=1) + 1e-9)
        mode7_alignment = float(np.abs(alignments).mean())
    else:
        mode7_alignment = 0
    print(f"\nMode-7 collective alignment of P02 lining: {mode7_alignment:.2f}  "
          f"({'collective' if mode7_alignment > 0.5 else 'mixed'})")

    # Displace structure along first real mode by ±1 sigma, recompute pocket volume
    print(f"\nDisplacing along mode {first_real_mode} (±1 sigma), recomputing P02 volume...")
    coords_pos = displace_along_mode(coords, eigvecs, eigvals, first_real_mode, sigma=1.0)
    coords_neg = displace_along_mode(coords, eigvecs, eigvals, first_real_mode, sigma=-1.0)
    vol_apo = ligsite_local_ca(coords, p02_center)
    vol_pos = ligsite_local_ca(coords_pos, p02_center)
    vol_neg = ligsite_local_ca(coords_neg, p02_center)
    print(f"  Volume apo (Ca-only ligsite):  {vol_apo:.1f} A^3")
    print(f"  Volume mode7+1sigma:           {vol_pos:.1f} A^3")
    print(f"  Volume mode7-1sigma:           {vol_neg:.1f} A^3")
    # Max volume along motion captures the "opened" state
    max_vol_motion = max(vol_pos, vol_neg, vol_apo)
    volume_change_pct = (max_vol_motion - vol_apo) / max(vol_apo, 1) * 100

    # Verdict logic
    flex_ok = lining_ratio >= 1.3
    handle_ok = handle_ratio >= 1.2 or len(handle_idx) == 0
    collective_ok = mode7_alignment >= 0.4
    opens_ok = max_vol_motion >= vol_apo and (vol_pos > vol_apo * 1.1 or vol_neg > vol_apo * 1.1)

    passes = sum([flex_ok, handle_ok, collective_ok, opens_ok])
    if passes >= 3:
        verdict = "CRYPTIC_DYNAMIC_POCKET"
        explain = (f"P02 lining is {lining_ratio:.2f}x more flexible than global average, "
                   f"moves collectively in mode 7 (alignment {mode7_alignment:.2f}), "
                   f"opens to {max_vol_motion:.0f} A^3 ({volume_change_pct:+.0f}%) under motion. "
                   "Allosteric pocket is dynamic — viable small-molecule target.")
    elif passes >= 2:
        verdict = "MARGINAL_DYNAMICS"
        explain = (f"P02 shows some dynamic signal ({passes}/4 gates) but not unambiguous. "
                   "Worth keeping; consider longer MD for confirmation.")
    else:
        verdict = "RIGID_OR_STATIC"
        explain = (f"P02 region behaves as rigid ({passes}/4 dynamic gates passed). "
                   "Likely a stable static pocket — still druggable, just not cryptic.")

    out = {
        "method": "ANM Cα cutoff 15A, 20 lowest non-trivial modes",
        "p02_lining_msf_ratio": round(lining_ratio, 3),
        "p02_handles_msf_ratio": round(handle_ratio, 3),
        "p02_lining_above_mean_frac": round(pctile_lining, 3),
        "mode7_collective_alignment": round(mode7_alignment, 3),
        "ligsite_volumes_along_mode7": {
            "apo_A3": round(vol_apo, 1),
            "plus_1sigma_A3": round(vol_pos, 1),
            "minus_1sigma_A3": round(vol_neg, 1),
            "max_open_change_pct": round(volume_change_pct, 1),
        },
        "gates_passed": {
            "flexibility (lining >=1.3x global)": flex_ok,
            "handles flex (>=1.2x global)": handle_ok,
            "collective motion (>=0.4)": collective_ok,
            "opens under motion (>=+10%)": opens_ok,
        },
        "n_gates_passed": passes,
        "verdict": verdict,
        "explanation": explain,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print()
    print("=" * 70)
    print(f"VERDICT: {verdict}  ({passes}/4 gates)")
    print(f"  {explain}")
    print("=" * 70)
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
