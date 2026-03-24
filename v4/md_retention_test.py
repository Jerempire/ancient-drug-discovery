"""MD Retention Test: Does the peptide stay in the channel or drift out?

Takes the best docked pose from DiffPepDock, runs 10ns MD with OpenMM,
measures peptide RMSD and distance from binding site over time.

If peptide stays (low RMSD, stable distance) = real binding.
If peptide drifts (rising RMSD, increasing distance) = transient/non-specific.

Run on Vast.ai GPU after docking.
"""
import os, sys, json, glob
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np

def find_best_pose(dock_dir, target, peptide):
    """Find the docked PDB with most contacts (best pose)."""
    pep_path = os.path.join(dock_dir, target, peptide)
    if not os.path.isdir(pep_path):
        print(f"  Not found: {pep_path}")
        return None

    pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))
    if not pdbs:
        return None

    best_pdb = None
    best_contacts = 0

    for pdb in pdbs[:16]:  # Check first 16
        try:
            pep_atoms = []
            rec_atoms = []
            with open(pdb) as f:
                for line in f:
                    if not line.startswith("ATOM"):
                        continue
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    if chain == "A":
                        pep_atoms.append(np.array([x, y, z]))
                    elif chain == "B":
                        rec_atoms.append(np.array([x, y, z]))

            if not pep_atoms or not rec_atoms:
                continue

            pep_xyz = np.array(pep_atoms)
            rec_xyz = np.array(rec_atoms)
            diff = pep_xyz[:, None, :] - rec_xyz[None, :, :]
            dists = np.sqrt(np.sum(diff**2, axis=-1))
            contacts = int(np.sum(dists < 4.5))

            if contacts > best_contacts:
                best_contacts = contacts
                best_pdb = pdb
        except:
            continue

    if best_pdb:
        print(f"  Best pose for {target}/{peptide}: {os.path.basename(best_pdb)} ({best_contacts} contacts)")
    return best_pdb


def run_md(pdb_path, output_prefix, n_steps=5000000, step_size=0.002):
    """Run MD simulation using OpenMM.

    n_steps=5000000 at 2fs = 10ns
    """
    try:
        from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation, StateDataReporter, DCDReporter
        from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
    except ImportError:
        from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation, StateDataReporter, DCDReporter
        from simtk.openmm import LangevinMiddleIntegrator, MonteCarloBarostat
        from simtk import unit

    print(f"\n  Loading {pdb_path}")
    pdb = PDBFile(pdb_path)

    # Use AMBER force field
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    print("  Adding solvent (10A padding)...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, padding=1.0*unit.nanometers)

    print(f"  System: {modeller.topology.getNumAtoms()} atoms")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=HBonds,
    )

    # Langevin integrator at 300K
    integrator = LangevinMiddleIntegrator(
        300*unit.kelvin,
        1.0/unit.picoseconds,
        step_size*unit.picoseconds,
    )

    # Add barostat for NPT
    system.addForce(MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print("  Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    # Equilibrate 100ps
    print("  Equilibrating (100ps)...")
    simulation.step(50000)

    # Get initial peptide positions (chain A CA atoms)
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

    # Find peptide CA atom indices
    pep_ca_indices = []
    rec_ca_indices = []
    for atom in modeller.topology.atoms():
        if atom.name == "CA":
            if atom.residue.chain.id == "A" or atom.residue.chain.index == 0:
                pep_ca_indices.append(atom.index)
            else:
                rec_ca_indices.append(atom.index)

    if not pep_ca_indices:
        # Try by residue index — peptide is first chain
        chains = list(modeller.topology.chains())
        if len(chains) >= 2:
            for atom in chains[0].atoms():
                if atom.name == "CA":
                    pep_ca_indices.append(atom.index)
            for atom in chains[1].atoms():
                if atom.name == "CA":
                    rec_ca_indices.append(atom.index)

    print(f"  Peptide CA atoms: {len(pep_ca_indices)}, Receptor CA atoms: {len(rec_ca_indices)}")

    initial_pep_pos = positions[pep_ca_indices]
    initial_pep_centroid = np.mean(initial_pep_pos, axis=0)

    # Production run with snapshots every 100ps
    print(f"  Production run ({n_steps * step_size / 1000:.0f}ns)...")

    report_interval = 50000  # Every 100ps
    n_reports = n_steps // report_interval

    rmsd_trajectory = []
    distance_trajectory = []
    time_points = []

    for i in range(n_reports):
        simulation.step(report_interval)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

        # Peptide RMSD from initial pose
        current_pep = positions[pep_ca_indices]
        diff = current_pep - initial_pep_pos
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        # Distance between peptide and receptor centroids
        pep_centroid = np.mean(current_pep, axis=0)
        rec_centroid = np.mean(positions[rec_ca_indices], axis=0) if rec_ca_indices else initial_pep_centroid
        dist = np.linalg.norm(pep_centroid - rec_centroid)

        time_ns = (i + 1) * report_interval * step_size / 1000
        rmsd_trajectory.append(float(rmsd))
        distance_trajectory.append(float(dist))
        time_points.append(float(time_ns))

        if (i + 1) % 10 == 0:
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"    t={time_ns:.1f}ns  RMSD={rmsd:.1f}A  dist={dist:.1f}A  PE={pe:.0f}kJ/mol")

    result = {
        "pdb": os.path.basename(pdb_path),
        "n_steps": n_steps,
        "time_ns": time_points,
        "rmsd_A": rmsd_trajectory,
        "centroid_distance_A": distance_trajectory,
        "final_rmsd": rmsd_trajectory[-1] if rmsd_trajectory else 0,
        "final_distance": distance_trajectory[-1] if distance_trajectory else 0,
        "avg_rmsd_last_2ns": float(np.mean(rmsd_trajectory[-20:])) if len(rmsd_trajectory) >= 20 else 0,
        "avg_distance_last_2ns": float(np.mean(distance_trajectory[-20:])) if len(distance_trajectory) >= 20 else 0,
        "rmsd_drift": rmsd_trajectory[-1] - rmsd_trajectory[0] if len(rmsd_trajectory) > 1 else 0,
        "distance_drift": distance_trajectory[-1] - distance_trajectory[0] if len(distance_trajectory) > 1 else 0,
    }

    out_path = f"/workspace/results/{output_prefix}_md.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")

    return result


def main():
    print("=" * 80)
    print("MD RETENTION TEST: Does peptide stay or drift?")
    print("=" * 80)

    DOCK_DIR = "/workspace/DiffPepBuilder/runs/docking"
    subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
    if subdirs:
        DOCK_DIR = subdirs[-1]

    # Test our lead peptide against ERAP2-K392 and IRAP
    peptide = "v41c_parent"  # or whatever the parent is named
    # Try multiple possible names
    for pep_name in ["v41c_parent", "v41_parent_glu_long", "pep_glu_long_01"]:
        k392_pdb = find_best_pose(DOCK_DIR, "erap2_k392", pep_name)
        irap_pdb = find_best_pose(DOCK_DIR, "irap", pep_name)
        if k392_pdb and irap_pdb:
            peptide = pep_name
            break

    if not k392_pdb or not irap_pdb:
        print("ERROR: Could not find docked poses. Check dock dir.")
        print(f"Dock dir: {DOCK_DIR}")
        print(f"Contents: {os.listdir(DOCK_DIR) if os.path.exists(DOCK_DIR) else 'NOT FOUND'}")
        return

    os.makedirs("/workspace/results", exist_ok=True)

    # Run shorter MD (2ns instead of 10ns for time/cost)
    # 2ns = 1,000,000 steps at 2fs
    N_STEPS = 1000000  # 2ns — enough to see if peptide drifts

    print(f"\n--- ERAP2-K392 ---")
    k392_result = run_md(k392_pdb, "erap2_k392", n_steps=N_STEPS)

    print(f"\n--- IRAP ---")
    irap_result = run_md(irap_pdb, "irap", n_steps=N_STEPS)

    # Compare
    print("\n" + "=" * 80)
    print("COMPARISON")
    print("-" * 80)
    print(f"{'Metric':>25s}  {'ERAP2-K392':>12s}  {'IRAP':>12s}  {'Verdict':>20s}")
    print("-" * 80)

    for metric, label in [
        ("final_rmsd", "Final RMSD (A)"),
        ("avg_rmsd_last_2ns", "Avg RMSD last phase"),
        ("rmsd_drift", "RMSD drift (A)"),
        ("final_distance", "Final centroid dist"),
        ("distance_drift", "Distance drift (A)"),
    ]:
        k = k392_result.get(metric, 0)
        i = irap_result.get(metric, 0)
        if "drift" in metric:
            verdict = "IRAP drifts more" if abs(i) > abs(k) * 1.5 else "Similar"
        elif "rmsd" in metric:
            verdict = "IRAP less stable" if i > k * 1.5 else "Similar"
        else:
            verdict = ""
        print(f"{label:>25s}  {k:>12.2f}  {i:>12.2f}  {verdict:>20s}")

    # Save combined
    combined = {
        "erap2_k392": k392_result,
        "irap": irap_result,
        "interpretation": {
            "erap2_stable": k392_result.get("rmsd_drift", 0) < 3.0,
            "irap_drifts": irap_result.get("rmsd_drift", 0) > k392_result.get("rmsd_drift", 0) * 1.5,
        }
    }
    with open("/workspace/results/md_retention_comparison.json", "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nSaved: /workspace/results/md_retention_comparison.json")


if __name__ == "__main__":
    main()
