"""Step 2: OpenMM CUDA — Run MD retention test on refined complex.

Takes PyRosetta-refined PDB, adds solvent, runs 2ns MD on GPU,
tracks peptide RMSD and centroid distance from receptor.

Usage: python simulate.py <input.pdb> <label> <output.json>
"""
import sys
import os
import json
import numpy as np

def run_md(input_pdb, label, output_json, n_steps=1000000):
    """Run 2ns MD with OpenMM CUDA platform."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit

    # Verify CUDA
    platform = Platform.getPlatformByName("CUDA")
    print(f"  Platform: {platform.getName()} (speed={platform.getSpeed()})")

    # PDBFixer to add missing atoms + hydrogens
    print(f"  Fixing {input_pdb}...")
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    # Solvate
    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    print("  Adding solvent...")
    modeller.addSolvent(ff, padding=0.8 * unit.nanometers)
    n_atoms = modeller.topology.getNumAtoms()
    print(f"  System: {n_atoms} atoms")

    # Build system
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
    )
    system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))

    # CUDA simulation
    sim = Simulation(
        modeller.topology, system, integrator,
        platform, {"CudaPrecision": "mixed"}
    )
    sim.context.setPositions(modeller.positions)

    # Minimize
    print("  Minimizing...")
    sim.minimizeEnergy(maxIterations=1000)

    # Equilibrate 100ps
    print("  Equilibrating (100ps)...")
    sim.step(50000)

    # Find peptide (chain A) and receptor (chain B) CA atoms
    chains = list(modeller.topology.chains())
    pep_ca = [a.index for a in chains[0].atoms() if a.name == "CA"]
    rec_ca = [a.index for a in chains[1].atoms() if a.name == "CA"] if len(chains) > 1 else []
    print(f"  Peptide CAs: {len(pep_ca)}, Receptor CAs: {len(rec_ca)}")

    # Record initial positions
    state = sim.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
    init_pep = pos[pep_ca].copy()

    # Production: 2ns, snapshots every 50ps
    interval = 25000  # 50ps
    total_snaps = n_steps // interval
    print(f"  Production ({n_steps * 0.002 / 1000:.0f}ns, {total_snaps} snapshots)...")

    rmsds, dists, times = [], [], []

    for i in range(total_snaps):
        sim.step(interval)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

        cur_pep = pos[pep_ca]
        rmsd = float(np.sqrt(np.mean(np.sum((cur_pep - init_pep) ** 2, axis=1))))

        if rec_ca:
            dist = float(np.linalg.norm(
                np.mean(cur_pep, axis=0) - np.mean(pos[rec_ca], axis=0)
            ))
        else:
            dist = 0.0

        t = (i + 1) * interval * 0.002 / 1000
        rmsds.append(round(rmsd, 2))
        dists.append(round(dist, 2))
        times.append(round(t, 3))

        if (i + 1) % 8 == 0:
            print(f"    t={t:.2f}ns  RMSD={rmsd:.1f}A  dist={dist:.1f}A")

    result = {
        "label": label,
        "pdb": os.path.basename(input_pdb),
        "n_atoms": n_atoms,
        "pep_residues": len(pep_ca),
        "rec_residues": len(rec_ca),
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
    }

    with open(output_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {output_json}")

    return result


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python simulate.py <input.pdb> <label> <output.json>")
        sys.exit(1)
    run_md(sys.argv[1], sys.argv[2], sys.argv[3])
