"""MD Retention Test — GPU version with explicit CUDA platform.

Runs 2ns MD on best docked pose in ERAP2-K392 and IRAP.
Measures: does the peptide stay (stable RMSD) or drift out (rising RMSD)?

Fix: OpenMM CUDA platform + GLY CB atom cleanup + PDBFixer.
"""
import os, sys, json, glob
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np

def find_best_pose(dock_dir, target, peptide):
    pep_path = os.path.join(dock_dir, target, peptide)
    pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))[:16]
    best_pdb, best_c = None, 0
    for pdb in pdbs:
        pep, rec = [], []
        with open(pdb) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if line[21] == "A":
                    pep.append(xyz)
                elif line[21] == "B":
                    rec.append(xyz)
        if not pep or not rec:
            continue
        diff = np.array(pep)[:, None, :] - np.array(rec)[None, :, :]
        c = int(np.sum(np.sqrt(np.sum(diff**2, axis=-1)) < 4.5))
        if c > best_c:
            best_c, best_pdb = c, pdb
    if best_pdb:
        print(f"  Best: {os.path.basename(best_pdb)} ({best_c} contacts)")
    return best_pdb


def fix_pdb(input_path, output_path):
    """Remove CB atoms from GLY residues (DiffPepDock artifact)."""
    lines = []
    with open(input_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                if resname == "GLY" and atom == "CB":
                    continue
            lines.append(line)
    with open(output_path, "w") as f:
        f.writelines(lines)


def run_md_gpu(pdb_path, label, n_steps=1000000):
    """Run 2ns MD with explicit CUDA platform."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit

    # Fix GLY artifacts
    fixed_raw = f"/workspace/{label}_raw_fixed.pdb"
    fix_pdb(pdb_path, fixed_raw)

    # PDBFixer: add missing atoms and hydrogens
    print(f"  Fixing PDB...")
    fixer = PDBFixer(filename=fixed_raw)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    fixed_path = f"/workspace/{label}_fixed.pdb"
    with open(fixed_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # Load and solvate
    pdb = PDBFile(fixed_path)
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(pdb.topology, pdb.positions)

    print(f"  Adding solvent...")
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers)
    n_atoms = modeller.topology.getNumAtoms()
    print(f"  System: {n_atoms} atoms")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
    )

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
    )
    system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))

    # EXPLICIT CUDA PLATFORM
    platform = Platform.getPlatformByName("CUDA")
    properties = {"CudaPrecision": "mixed"}
    print(f"  Platform: CUDA (mixed precision)")

    simulation = Simulation(
        modeller.topology, system, integrator, platform, properties
    )
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print(f"  Minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)

    # Equilibrate 100ps
    print(f"  Equilibrating (100ps)...")
    simulation.step(50000)

    # Find peptide and receptor CA atoms
    chains = list(modeller.topology.chains())
    pep_ca = [a.index for a in chains[0].atoms() if a.name == "CA"]
    rec_ca = []
    if len(chains) > 1:
        rec_ca = [a.index for a in chains[1].atoms() if a.name == "CA"]
    print(f"  Peptide CAs: {len(pep_ca)}, Receptor CAs: {len(rec_ca)}")

    # Record initial positions
    state = simulation.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
    init_pep = pos[pep_ca].copy()
    init_rec_centroid = np.mean(pos[rec_ca], axis=0) if rec_ca else np.zeros(3)

    # Production: 2ns with snapshots every 50ps
    step_size = 0.002  # ps
    interval = 25000  # 50ps
    total_snapshots = n_steps // interval

    print(f"  Production ({n_steps * step_size / 1000:.0f}ns, {total_snapshots} snapshots)...")

    rmsds = []
    dists = []
    times = []

    for i in range(total_snapshots):
        simulation.step(interval)
        state = simulation.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

        cur_pep = pos[pep_ca]
        rmsd = float(np.sqrt(np.mean(np.sum((cur_pep - init_pep) ** 2, axis=1))))

        if rec_ca:
            pep_cent = np.mean(cur_pep, axis=0)
            rec_cent = np.mean(pos[rec_ca], axis=0)
            dist = float(np.linalg.norm(pep_cent - rec_cent))
        else:
            dist = 0.0

        t = (i + 1) * interval * step_size / 1000  # ns

        rmsds.append(round(rmsd, 2))
        dists.append(round(dist, 2))
        times.append(round(t, 3))

        if (i + 1) % 8 == 0:
            print(f"    t={t:.2f}ns  RMSD={rmsd:.1f}A  dist={dist:.1f}A")

    result = {
        "label": label,
        "pdb": os.path.basename(pdb_path),
        "n_atoms": n_atoms,
        "n_steps": n_steps,
        "time_ns": times,
        "rmsd_A": rmsds,
        "distance_A": dists,
        "final_rmsd": rmsds[-1],
        "final_dist": dists[-1],
        "avg_rmsd_last_1ns": round(float(np.mean(rmsds[-20:])), 2),
        "avg_dist_last_1ns": round(float(np.mean(dists[-20:])), 2),
        "rmsd_drift": round(rmsds[-1] - rmsds[0], 2),
        "dist_drift": round(dists[-1] - dists[0], 2),
        "max_rmsd": round(float(max(rmsds)), 2),
        "max_dist": round(float(max(dists)), 2),
    }

    return result


def main():
    print("=" * 70)
    print("MD RETENTION TEST (GPU) — Does peptide stay or drift?")
    print("=" * 70)

    # Check CUDA
    from openmm import Platform
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        print(f"  Platform {i}: {p.getName()}")

    DOCK_DIR = "/workspace/DiffPepBuilder/runs/docking"
    subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
    if subdirs:
        DOCK_DIR = subdirs[-1]

    # Find best poses for parent peptide
    peptide = "v41c_parent"
    print(f"\nFinding best poses for {peptide}...")

    print("\n  ERAP2-K392:")
    k392_pdb = find_best_pose(DOCK_DIR, "erap2_k392", peptide)
    print("  IRAP:")
    irap_pdb = find_best_pose(DOCK_DIR, "irap", peptide)

    if not k392_pdb or not irap_pdb:
        print("ERROR: Could not find poses")
        return

    os.makedirs("/workspace/results", exist_ok=True)

    # Run 2ns MD on each
    print("\n" + "=" * 70)
    print("ERAP2-K392 MD")
    print("=" * 70)
    k392_result = run_md_gpu(k392_pdb, "erap2_k392")

    print("\n" + "=" * 70)
    print("IRAP MD")
    print("=" * 70)
    irap_result = run_md_gpu(irap_pdb, "irap")

    # Comparison
    print("\n" + "=" * 70)
    print("COMPARISON")
    print("-" * 70)
    print(f"{'Metric':>25s}  {'ERAP2-K392':>12s}  {'IRAP':>12s}  {'Verdict':>20s}")
    print("-" * 70)

    metrics = [
        ("final_rmsd", "Final RMSD (A)"),
        ("avg_rmsd_last_1ns", "Avg RMSD last 1ns"),
        ("max_rmsd", "Max RMSD (A)"),
        ("rmsd_drift", "RMSD drift (A)"),
        ("final_dist", "Final distance (A)"),
        ("dist_drift", "Distance drift (A)"),
    ]

    for key, label in metrics:
        k = k392_result[key]
        i = irap_result[key]
        if "drift" in key or "max" in key:
            verdict = "IRAP WORSE" if i > k * 1.3 else "SIMILAR" if i < k * 1.3 else ""
        else:
            verdict = ""
        print(f"{label:>25s}  {k:>12.2f}  {i:>12.2f}  {verdict:>20s}")

    # Overall verdict
    k_drift = k392_result["rmsd_drift"]
    i_drift = irap_result["rmsd_drift"]
    if i_drift > k_drift * 2:
        verdict = "IRAP PEPTIDE DRIFTS — binding is transient (GOOD NEWS)"
    elif i_drift > k_drift * 1.5:
        verdict = "IRAP less stable — some evidence of weaker binding"
    else:
        verdict = "SIMILAR STABILITY — IRAP binding appears persistent"

    print(f"\nVERDICT: {verdict}")
    print(f"  ERAP2 drift: {k_drift:.2f}A")
    print(f"  IRAP drift:  {i_drift:.2f}A")

    # Save
    combined = {
        "erap2_k392": k392_result,
        "irap": irap_result,
        "verdict": verdict,
    }
    with open("/workspace/results/md_retention.json", "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nSaved: /workspace/results/md_retention.json")


if __name__ == "__main__":
    main()
