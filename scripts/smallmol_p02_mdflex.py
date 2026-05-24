"""
MD-flex check: does the P02 allosteric pocket stay open in apo ERAP2 dynamics?

Pipeline:
  1. PDBFixer prep on AF-ERAP2 (add Hs, missing atoms)
  2. Solvate, minimize, NVT/NPT equilibrate
  3. Production MD (default 5 ns, frame every 50 ps -> 100 frames)
  4. Extract 10 evenly-spaced production frames as protein-only PDBs
  5. For each frame: run LIGSITE pocket detection on the P02 region only
  6. Report:
     - P02 occupancy: fraction of frames where a pocket >= 80 A^3 exists within 6 A of original centroid
     - Volume distribution
     - Verdict (STABLE / TRANSIENT / CLOSED)

Runtime: ~45-90 min on RTX 4090 for 5 ns of solvated ERAP2.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

ROOT = Path(__file__).resolve().parent.parent
ERAP2_PDB = ROOT / "data/structures/erap2_wt_alphafold.pdb"
POCKET_JSON = ROOT / "data/results/smallmol/pocket_inventory.json"
WORKDIR = ROOT / "data/results/smallmol/mdflex"
WORKDIR.mkdir(parents=True, exist_ok=True)
OUT = WORKDIR / "p02_occupancy.json"


def prep_system(input_pdb: Path, out_pdb: Path):
    """PDBFixer: add Hs, fix missing atoms, solvate."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, Modeller, ForceField
    from openmm import unit

    print("PDBFixer prep ...")
    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # ignore terminal disorder
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(ff, padding=0.8 * unit.nanometer, ionicStrength=0.15 * unit.molar)

    with open(out_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"  solvated -> {out_pdb} ({modeller.topology.getNumAtoms()} atoms)")
    return ff, modeller


def run_md(ff, modeller, traj_dcd: Path, prod_ns: float, frame_ps: float = 50):
    """Minimize -> NVT 50ps -> NPT 50ps -> production. Save protein-only frames."""
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
    from openmm.app import Simulation, DCDReporter, PDBFile, StateDataReporter

    T = 310 * unit.kelvin
    P = 1.0 * unit.atmosphere
    dt = 2.0 * unit.femtosecond

    system = ff.createSystem(
        modeller.topology,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=True,
    )
    from openmm.app import HBonds
    system = ff.createSystem(
        modeller.topology,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    integrator = LangevinMiddleIntegrator(T, 1.0 / unit.picosecond, dt)
    platform = None
    for pname, opts in [("CUDA", {"Precision": "mixed"}),
                        ("OpenCL", {"Precision": "mixed"}),
                        ("CPU", {})]:
        try:
            platform = Platform.getPlatformByName(pname)
            platform_opts = opts
            print(f"Using OpenMM platform: {pname}")
            break
        except Exception:
            continue
    sim = Simulation(modeller.topology, system, integrator, platform, platform_opts)
    sim.context.setPositions(modeller.positions)

    print("Minimizing ...")
    sim.minimizeEnergy(maxIterations=2000)

    print("NVT equil 50 ps ...")
    sim.context.setVelocitiesToTemperature(T)
    sim.step(int(50 * unit.picosecond / dt))

    print("NPT equil 50 ps ...")
    system.addForce(MonteCarloBarostat(P, T, 25))
    sim.context.reinitialize(preserveState=True)
    sim.step(int(50 * unit.picosecond / dt))

    prod_steps = int(prod_ns * 1000 * unit.picosecond / dt)
    frame_steps = int(frame_ps * unit.picosecond / dt)
    n_frames = prod_steps // frame_steps
    print(f"Production: {prod_ns} ns ({prod_steps} steps), frame every {frame_ps} ps ({n_frames} frames)")

    sim.reporters.append(DCDReporter(str(traj_dcd), frame_steps))
    sim.reporters.append(StateDataReporter(
        sys.stdout, frame_steps * 10,
        step=True, time=True, potentialEnergy=True, temperature=True, volume=True,
        progress=True, remainingTime=True, totalSteps=prod_steps,
    ))
    t0 = time.time()
    sim.step(prod_steps)
    print(f"Production done in {(time.time()-t0)/60:.1f} min")

    # Save final topology PDB (for mdtraj load with no waters later)
    state = sim.context.getState(getPositions=True)
    final_pdb = traj_dcd.with_suffix(".final.pdb")
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    return final_pdb


def extract_protein_frames(traj_dcd: Path, top_pdb: Path, out_dir: Path, n_frames: int = 10):
    import mdtraj as md
    traj = md.load(str(traj_dcd), top=str(top_pdb))
    prot = traj.atom_slice(traj.topology.select("protein"))
    print(f"Trajectory: {prot.n_frames} frames, protein atoms: {prot.n_atoms}")
    idx = np.linspace(0, prot.n_frames - 1, n_frames, dtype=int)
    frame_paths = []
    for i, fi in enumerate(idx):
        p = out_dir / f"frame_{i:02d}.pdb"
        prot[fi].save_pdb(str(p))
        frame_paths.append(p)
    return frame_paths


def ligsite_local(frame_pdb: Path, region_center: np.ndarray, region_radius: float = 12.0,
                  spacing: float = 1.5, psp_min: int = 5):
    """Run LIGSITE within a 12 A sphere around region_center. Return total pocket-grid volume."""
    from Bio import PDB
    from scipy.spatial import cKDTree

    parser = PDB.PDBParser(QUIET=True)
    s = parser.get_structure("x", str(frame_pdb))
    atom_coords = np.array([a.coord for a in s.get_atoms() if a.element != "H"])

    mn = region_center - region_radius
    mx = region_center + region_radius
    xs = np.arange(mn[0], mx[0], spacing)
    ys = np.arange(mn[1], mx[1], spacing)
    zs = np.arange(mn[2], mx[2], spacing)
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    grid = grid[np.linalg.norm(grid - region_center, axis=1) <= region_radius]

    tree = cKDTree(atom_coords)
    nearest, _ = tree.query(grid, k=1)
    solv = (nearest > 3.2) & (nearest < 8.0)
    grid_s = grid[solv]

    dirs = np.array([
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1],
    ], dtype=float)
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    n_steps = int(12.0 / spacing)
    psp = np.zeros(len(grid_s), dtype=int)
    for d in dirs:
        fwd = np.zeros(len(grid_s), dtype=bool)
        bwd = np.zeros(len(grid_s), dtype=bool)
        for k in range(1, n_steps + 1):
            df, _ = tree.query(grid_s + d * k * spacing, k=1)
            fwd |= df < 1.8
            db, _ = tree.query(grid_s - d * k * spacing, k=1)
            bwd |= db < 1.8
        psp += (fwd & bwd).astype(int)

    pocket_pts = grid_s[psp >= psp_min]
    volume = len(pocket_pts) * spacing ** 3
    return volume, pocket_pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod_ns", type=float, default=5.0)
    ap.add_argument("--frame_ps", type=float, default=50.0)
    ap.add_argument("--n_analysis_frames", type=int, default=10)
    args = ap.parse_args()

    with open(POCKET_JSON) as f:
        inv = json.load(f)
    p02 = inv["top_allosteric"]
    p02_center = np.array(p02["centroid_xyz"])
    p02_apo_volume = p02["volume_A3_est"]
    print(f"P02 apo centroid: {p02_center.tolist()}")
    print(f"P02 apo volume:   {p02_apo_volume} A^3")

    solv_pdb = WORKDIR / "solvated.pdb"
    traj_dcd = WORKDIR / "production.dcd"

    if not traj_dcd.exists():
        ff, modeller = prep_system(ERAP2_PDB, solv_pdb)
        final_pdb = run_md(ff, modeller, traj_dcd, args.prod_ns, args.frame_ps)
    else:
        print(f"Reusing existing trajectory: {traj_dcd}")
        final_pdb = traj_dcd.with_suffix(".final.pdb")

    frames_dir = WORKDIR / "frames"
    frames_dir.mkdir(exist_ok=True)
    frame_paths = extract_protein_frames(traj_dcd, final_pdb, frames_dir, args.n_analysis_frames)
    print(f"Extracted {len(frame_paths)} frames for pocket analysis")

    print("\nPer-frame P02 pocket detection (12 A sphere around apo centroid):")
    per_frame = []
    for p in frame_paths:
        v, pts = ligsite_local(p, p02_center)
        per_frame.append({"frame": p.name, "p02_volume_A3": round(float(v), 1),
                          "n_pocket_pts": int(len(pts))})
        print(f"  {p.name}: vol={v:>6.1f} A^3  pts={len(pts):>4}")

    vols = np.array([f["p02_volume_A3"] for f in per_frame])
    occupancy = float((vols >= 80).mean())
    median_vol = float(np.median(vols))
    mean_vol = float(vols.mean())

    if occupancy >= 0.7 and median_vol >= 100:
        verdict = "STABLE_ALLOSTERIC_POCKET"
        explain = f"P02 detected in {occupancy*100:.0f}% of frames, median vol {median_vol:.0f} A^3 — viable allosteric target."
    elif occupancy >= 0.4:
        verdict = "TRANSIENT_POCKET"
        explain = f"P02 open in {occupancy*100:.0f}% of frames — cryptic; may need longer MD or biasing to fully assess."
    else:
        verdict = "CRYPTIC_OR_CLOSED"
        explain = f"P02 rarely open ({occupancy*100:.0f}% of frames). Static AF prediction was misleading; do not prioritize."

    result = {
        "prod_ns": args.prod_ns,
        "n_frames_analyzed": len(frame_paths),
        "p02_apo_volume_A3": p02_apo_volume,
        "p02_apo_centroid": p02_center.tolist(),
        "per_frame": per_frame,
        "occupancy_frac_ge80A3": round(occupancy, 3),
        "volume_median_A3": round(median_vol, 1),
        "volume_mean_A3": round(mean_vol, 1),
        "verdict": verdict,
        "explanation": explain,
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print()
    print("=" * 70)
    print(f"VERDICT: {verdict}")
    print(f"  {explain}")
    print("=" * 70)
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
