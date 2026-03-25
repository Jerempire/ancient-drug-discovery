"""
OpenMM 10ns MD Protocol for ERAP2 Short-Stack Peptide Validation.

Runs a full MD simulation: minimize → NVT equil → NPT equil → 10ns production.
Computes RMSD, COM drift, contact fraction, zinc distance, H-bond occupancy,
and 3-region contact distribution during production.

Usage:
  python openmm_md_protocol.py <complex.pdb> <output_dir> [--production_ns 10]

Example:
  python openmm_md_protocol.py /workspace/structures/VAGSAF_vs_erap2k392.pdb /workspace/results/md_VAGSAF_erap2k392/

Dependencies: openmm, mdtraj, pdbfixer, numpy
Install: conda install -c conda-forge openmm mdtraj pdbfixer numpy
  OR: pip install openmm mdtraj pdbfixer numpy

GPU: Required. RTX 4090 does ~150 ns/day for ~10K atom system.
10ns ≈ 1.5-2 hours per simulation on full ERAP2 (~960 residues + peptide + water).

Output:
  trajectory.dcd          — full trajectory (frames every 10 ps)
  trajectory_prot.dcd     — protein+peptide only (no water, smaller)
  energies.csv            — potential energy, temperature, volume per ps
  rmsd.csv                — peptide RMSD per frame
  com_drift.csv           — peptide center-of-mass drift per frame
  contacts.csv            — native contact fraction per frame
  zinc_distance.csv       — min peptide-zinc distance per frame
  hbond_occupancy.json    — key H-bond occupancies
  region_contacts.csv     — floor/wall/ceiling contact counts per frame
  summary.json            — final verdict + all metrics
"""
import sys
import os
import json
import time
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
TEMPERATURE = 310  # K (physiological)
PRESSURE = 1.0     # atm
TIMESTEP_FS = 2.0  # femtoseconds (production)
PADDING_NM = 1.0   # water box padding
IONIC_STRENGTH = 0.15  # M NaCl
MINIMIZE_STEPS = 1000
NVT_EQUIL_PS = 100
NPT_EQUIL_PS = 100
REPORT_INTERVAL_PS = 10  # save coordinates every 10 ps
ENERGY_INTERVAL_PS = 1   # save energies every 1 ps

# ERAP2 catalytic zinc coordinating residues (1-indexed)
ZINC_COORD_RESIDUES = [370, 374, 393]  # H370, H374, E393 (HEXXH motif)

# 3-Region contact definitions (ERAP2 residue numbers, 1-indexed)
FLOOR_RESIDUES = [392, 398]       # K392 selectivity handle, Y398
WALL_RESIDUES = [403, 406]        # A403, A406 (IRAP/ERAP1 evasion)
CEILING_RESIDUES = [412, 414]     # Q412, D414 (channel cap)

CONTACT_CUTOFF_NM = 0.45  # 4.5 Angstrom


def setup_system(pdb_path):
    """Load PDB, fix missing atoms/residues, solvate, add ions."""
    from pdbfixer import PDBFixer
    from openmm.app import Modeller, ForceField, PME, HBonds
    from openmm import unit

    print("Loading and fixing PDB...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    print("Setting up force field...")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    print("Solvating...")
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        forcefield,
        padding=PADDING_NM * unit.nanometers,
        ionicStrength=IONIC_STRENGTH * unit.molar,
    )

    print("Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
    )

    return system, modeller.topology, modeller.positions, forcefield


def find_peptide_indices(topology):
    """Find atom indices for the peptide chain (chain B, or last chain)."""
    chains = list(topology.chains())
    if len(chains) < 2:
        print("WARNING: Only 1 chain found. Assuming last residues are peptide.")
        # Fallback: last 6-7 residues
        all_res = list(topology.residues())
        peptide_res = all_res[-7:]
    else:
        # Peptide is chain B (second chain)
        peptide_res = list(chains[-1].residues())

    peptide_atoms = []
    peptide_ca = []
    peptide_heavy = []
    for res in peptide_res:
        for atom in res.atoms():
            peptide_atoms.append(atom.index)
            if atom.element.symbol != "H":
                peptide_heavy.append(atom.index)
            if atom.name == "CA":
                peptide_ca.append(atom.index)

    return {
        "all": peptide_atoms,
        "heavy": peptide_heavy,
        "ca": peptide_ca,
        "residues": peptide_res,
    }


def find_protein_indices(topology):
    """Find protein (non-water, non-ion) atom indices."""
    protein_atoms = []
    for chain in topology.chains():
        for res in chain.residues():
            if res.name not in ("HOH", "WAT", "NA", "CL", "SOL"):
                for atom in res.atoms():
                    protein_atoms.append(atom.index)
    return protein_atoms


def find_zinc_index(topology):
    """Find zinc atom index if present."""
    for atom in topology.atoms():
        if atom.element and atom.element.symbol == "ZN":
            return atom.index
    # No zinc found — look for zinc-coordinating residues' metal-binding atoms
    return None


def run_simulation(system, topology, positions, output_dir, production_ns):
    """Run minimization + equilibration + production."""
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
    from openmm.app import Simulation, DCDReporter, StateDataReporter, PDBFile
    from openmm import unit

    os.makedirs(output_dir, exist_ok=True)

    # Choose platform
    try:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"DeviceName": "0", "Precision": "mixed"}
        print("Using CUDA platform")
    except Exception:
        try:
            platform = Platform.getPlatformByName("OpenCL")
            properties = {}
            print("Using OpenCL platform")
        except Exception:
            platform = Platform.getPlatformByName("CPU")
            properties = {}
            print("WARNING: Using CPU platform — will be very slow!")

    # Integrator
    integrator = LangevinMiddleIntegrator(
        TEMPERATURE * unit.kelvin,
        1.0 / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds,
    )

    # Simulation
    sim = Simulation(topology, system, integrator, platform, properties)
    sim.context.setPositions(positions)

    # --- Phase 1: Minimize ---
    print("Minimizing (%d steps)..." % MINIMIZE_STEPS)
    sim.minimizeEnergy(maxIterations=MINIMIZE_STEPS)

    # Save minimized structure
    state = sim.context.getState(getPositions=True)
    with open(os.path.join(output_dir, "minimized.pdb"), "w") as f:
        PDBFile.writeFile(topology, state.getPositions(), f)

    # --- Phase 2: NVT Equilibration ---
    nvt_steps = int(NVT_EQUIL_PS / (TIMESTEP_FS / 1000))
    print("NVT equilibration (%d ps, %d steps)..." % (NVT_EQUIL_PS, nvt_steps))
    sim.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
    sim.step(nvt_steps)

    # --- Phase 3: NPT Equilibration ---
    system.addForce(MonteCarloBarostat(PRESSURE * unit.atmospheres, TEMPERATURE * unit.kelvin))
    sim.context.reinitialize(preserveState=True)
    npt_steps = int(NPT_EQUIL_PS / (TIMESTEP_FS / 1000))
    print("NPT equilibration (%d ps, %d steps)..." % (NPT_EQUIL_PS, npt_steps))
    sim.step(npt_steps)

    # --- Phase 4: Production ---
    production_steps = int(production_ns * 1000 / (TIMESTEP_FS / 1000))
    report_interval = int(REPORT_INTERVAL_PS / (TIMESTEP_FS / 1000))
    energy_interval = int(ENERGY_INTERVAL_PS / (TIMESTEP_FS / 1000))

    print("Production MD (%.1f ns, %d steps, saving every %d ps)..." % (
        production_ns, production_steps, REPORT_INTERVAL_PS))

    # Reporters
    traj_path = os.path.join(output_dir, "trajectory.dcd")
    energy_path = os.path.join(output_dir, "energies.csv")
    sim.reporters.append(DCDReporter(traj_path, report_interval))
    sim.reporters.append(StateDataReporter(
        energy_path, energy_interval,
        step=True, time=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True,
    ))
    sim.reporters.append(StateDataReporter(
        sys.stdout, report_interval * 10,
        step=True, time=True, speed=True, remainingTime=True,
        totalSteps=production_steps,
    ))

    t0 = time.time()
    sim.step(production_steps)
    wall_time = (time.time() - t0) / 3600

    # Save final state
    state = sim.context.getState(getPositions=True, getVelocities=True)
    with open(os.path.join(output_dir, "final_state.xml"), "w") as f:
        f.write(sim.context.createCheckpoint().hex() if hasattr(sim.context, "createCheckpoint") else "")
    with open(os.path.join(output_dir, "final.pdb"), "w") as f:
        PDBFile.writeFile(topology, state.getPositions(), f)

    print("Production complete. Wall time: %.2f hours" % wall_time)
    return traj_path, wall_time


def analyze_trajectory(traj_path, pdb_path, output_dir, peptide_name, target_name):
    """Compute all analysis metrics from the trajectory."""
    import mdtraj as md

    print("Loading trajectory for analysis...")
    traj = md.load(traj_path, top=pdb_path)
    n_frames = traj.n_frames
    print("Loaded %d frames" % n_frames)

    # Identify peptide and protein
    top = traj.topology
    chains = list(top.chains)
    if len(chains) >= 2:
        peptide_chain = chains[-1]
        protein_chain = chains[0]
    else:
        print("WARNING: Single chain — skipping detailed analysis")
        return {}

    peptide_atoms = top.select("chainid %d" % peptide_chain.index)
    peptide_heavy = top.select("chainid %d and not element H" % peptide_chain.index)
    peptide_ca = top.select("chainid %d and name CA" % peptide_chain.index)
    protein_ca = top.select("chainid %d and name CA" % protein_chain.index)

    # --- A. RMSD ---
    print("Computing RMSD...")
    # Align on protein CA, measure peptide RMSD
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=protein_ca)
    if len(peptide_heavy) > 0:
        rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0, atom_indices=peptide_heavy) * 10  # nm to Å
    else:
        rmsd = np.zeros(n_frames)

    np.savetxt(os.path.join(output_dir, "rmsd.csv"),
               np.column_stack([np.arange(n_frames) * REPORT_INTERVAL_PS / 1000, rmsd]),
               header="time_ns,rmsd_angstrom", delimiter=",", fmt="%.4f")

    # --- B. COM Drift ---
    print("Computing COM drift...")
    peptide_xyz = traj_aligned.xyz[:, peptide_heavy, :]
    com = peptide_xyz.mean(axis=1)  # (n_frames, 3)
    com_drift = np.sqrt(((com - com[0]) ** 2).sum(axis=1)) * 10  # nm to Å

    np.savetxt(os.path.join(output_dir, "com_drift.csv"),
               np.column_stack([np.arange(n_frames) * REPORT_INTERVAL_PS / 1000, com_drift]),
               header="time_ns,com_drift_angstrom", delimiter=",", fmt="%.4f")

    # --- C. Contact Fraction ---
    print("Computing contact fraction...")
    # Native contacts: protein-peptide contacts in frame 0
    pairs = []
    for pa in peptide_heavy:
        for prot_a in top.select("chainid %d and not element H" % protein_chain.index):
            pairs.append((pa, prot_a))

    if pairs:
        distances_0 = md.compute_distances(traj[0], pairs)[0]
        native_mask = distances_0 < CONTACT_CUTOFF_NM
        native_pairs = [pairs[i] for i in range(len(pairs)) if native_mask[i]]
        n_native = len(native_pairs)

        if n_native > 0 and native_pairs:
            all_distances = md.compute_distances(traj_aligned, native_pairs)
            contact_fraction = (all_distances < CONTACT_CUTOFF_NM).sum(axis=1) / n_native
        else:
            contact_fraction = np.zeros(n_frames)
    else:
        contact_fraction = np.zeros(n_frames)
        n_native = 0

    np.savetxt(os.path.join(output_dir, "contacts.csv"),
               np.column_stack([np.arange(n_frames) * REPORT_INTERVAL_PS / 1000, contact_fraction]),
               header="time_ns,contact_fraction", delimiter=",", fmt="%.4f")

    # --- D. Zinc Distance ---
    print("Computing zinc distance...")
    zn_atoms = top.select("element Zn or element ZN")
    if len(zn_atoms) > 0:
        zn_pairs = [(pa, zn_atoms[0]) for pa in peptide_heavy]
        zn_distances = md.compute_distances(traj_aligned, zn_pairs)
        min_zn_dist = zn_distances.min(axis=1) * 10  # nm to Å
    else:
        # No zinc — use approximate position from coordinating residues
        coord_atoms = top.select("resid %s and name CA" % " ".join(str(r-1) for r in ZINC_COORD_RESIDUES))
        if len(coord_atoms) > 0:
            zn_approx = traj_aligned.xyz[:, coord_atoms, :].mean(axis=1)  # approximate zinc position
            pep_com = peptide_xyz.mean(axis=1)
            min_zn_dist = np.sqrt(((pep_com - zn_approx) ** 2).sum(axis=1)) * 10
        else:
            min_zn_dist = np.full(n_frames, 99.0)

    np.savetxt(os.path.join(output_dir, "zinc_distance.csv"),
               np.column_stack([np.arange(n_frames) * REPORT_INTERVAL_PS / 1000, min_zn_dist]),
               header="time_ns,min_zinc_dist_angstrom", delimiter=",", fmt="%.4f")

    # --- E. 3-Region Contact Distribution ---
    print("Computing 3-region contacts...")
    region_data = {"floor": [], "wall": [], "ceiling": []}

    for region_name, residue_ids in [("floor", FLOOR_RESIDUES), ("wall", WALL_RESIDUES), ("ceiling", CEILING_RESIDUES)]:
        # Select protein atoms in this region (0-indexed residue IDs)
        region_sel = " or ".join(["resid %d" % (r - 1) for r in residue_ids])
        region_atoms = top.select("chainid %d and (%s) and not element H" % (protein_chain.index, region_sel))

        if len(region_atoms) > 0 and len(peptide_heavy) > 0:
            rpairs = [(pa, ra) for pa in peptide_heavy for ra in region_atoms]
            if rpairs:
                rdist = md.compute_distances(traj_aligned, rpairs)
                region_contacts = (rdist < CONTACT_CUTOFF_NM).sum(axis=1)
                region_data[region_name] = region_contacts.tolist()
            else:
                region_data[region_name] = [0] * n_frames
        else:
            region_data[region_name] = [0] * n_frames

    region_array = np.column_stack([
        np.arange(n_frames) * REPORT_INTERVAL_PS / 1000,
        region_data["floor"],
        region_data["wall"],
        region_data["ceiling"],
    ])
    np.savetxt(os.path.join(output_dir, "region_contacts.csv"),
               region_array, header="time_ns,floor_contacts,wall_contacts,ceiling_contacts",
               delimiter=",", fmt="%.4f")

    # Region verdict
    avg_floor = np.mean(region_data["floor"])
    avg_wall = np.mean(region_data["wall"])
    avg_ceiling = np.mean(region_data["ceiling"])
    regions_contacted = sum(1 for avg in [avg_floor, avg_wall, avg_ceiling] if avg > 0.5)

    if regions_contacted >= 3:
        region_verdict = "VOLUME_FILLER"
    elif regions_contacted == 2:
        region_verdict = "PARTIAL"
    else:
        region_verdict = "LOCAL_ANCHOR"

    # --- F. H-bond Occupancy ---
    print("Computing H-bond occupancy...")
    hbonds = md.baker_hubbard(traj_aligned, freq=0.1)  # bonds present in >10% of frames
    hbond_info = []
    for hb in hbonds:
        donor = top.atom(hb[0])
        acceptor = top.atom(hb[2])
        # Check if involves peptide and K392 region
        is_peptide_donor = donor.residue.chain.index == peptide_chain.index
        is_peptide_acceptor = acceptor.residue.chain.index == peptide_chain.index
        if is_peptide_donor or is_peptide_acceptor:
            # Compute occupancy
            distances = md.compute_distances(traj_aligned, [[hb[0], hb[2]]])
            occupancy = (distances < 0.35).mean()  # 3.5 Å
            hbond_info.append({
                "donor": "%s-%s-%s" % (donor.residue.name, donor.residue.resSeq, donor.name),
                "acceptor": "%s-%s-%s" % (acceptor.residue.name, acceptor.residue.resSeq, acceptor.name),
                "occupancy": float(occupancy),
                "involves_k392": any(r.resSeq == 392 for r in [donor.residue, acceptor.residue]),
            })

    hbond_info.sort(key=lambda x: -x["occupancy"])
    with open(os.path.join(output_dir, "hbond_occupancy.json"), "w") as f:
        json.dump(hbond_info[:20], f, indent=2)  # top 20

    # K392 H-bond occupancy
    k392_hbonds = [h for h in hbond_info if h["involves_k392"]]
    k392_occ = max([h["occupancy"] for h in k392_hbonds]) if k392_hbonds else 0

    # --- Compile Summary ---
    rmsd_final = float(rmsd[-1])
    rmsd_mean = float(rmsd.mean())
    rmsd_max = float(rmsd.max())
    com_final = float(com_drift[-1])
    cf_final = float(contact_fraction[-1])
    zn_mean = float(min_zn_dist.mean())
    zn_min = float(min_zn_dist.min())

    # Verdict
    if rmsd_mean < 2.0:
        verdict = "LOCKED"
    elif rmsd_mean < 5.0:
        verdict = "SHIFTING"
    else:
        verdict = "DRIFTING"

    summary = {
        "peptide": peptide_name,
        "target": target_name,
        "rmsd_final": round(rmsd_final, 3),
        "rmsd_mean": round(rmsd_mean, 3),
        "rmsd_max": round(rmsd_max, 3),
        "com_drift_final": round(com_final, 3),
        "contact_fraction_final": round(cf_final, 3),
        "n_native_contacts": n_native,
        "zinc_distance_mean": round(zn_mean, 3),
        "zinc_distance_min": round(zn_min, 3),
        "floor_contacts_avg": round(avg_floor, 3),
        "wall_contacts_avg": round(avg_wall, 3),
        "ceiling_contacts_avg": round(avg_ceiling, 3),
        "region_verdict": region_verdict,
        "hbond_k392_occupancy": round(k392_occ, 3),
        "n_peptide_hbonds": len(hbond_info),
        "verdict": verdict,
        "simulation_ns": float(n_frames * REPORT_INTERVAL_PS / 1000),
        "n_frames": n_frames,
    }

    with open(os.path.join(output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY: %s vs %s" % (peptide_name, target_name))
    print("=" * 60)
    print("RMSD (mean/max): %.1f / %.1f Å" % (rmsd_mean, rmsd_max))
    print("COM drift (final): %.1f Å" % com_final)
    print("Contact fraction (final): %.2f" % cf_final)
    print("Zinc distance (mean/min): %.1f / %.1f Å" % (zn_mean, zn_min))
    print("Region contacts (F/W/C): %.1f / %.1f / %.1f → %s" % (avg_floor, avg_wall, avg_ceiling, region_verdict))
    print("K392 H-bond occupancy: %.1f%%" % (k392_occ * 100))
    print("VERDICT: %s" % verdict)
    print("=" * 60)

    return summary


def main():
    parser = argparse.ArgumentParser(description="OpenMM 10ns MD for ERAP2 peptide validation")
    parser.add_argument("pdb", help="Input PDB file (protein + peptide complex)")
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument("--production_ns", type=float, default=10.0, help="Production MD length in ns")
    args = parser.parse_args()

    # Extract peptide and target names from filename
    basename = os.path.splitext(os.path.basename(args.pdb))[0]
    parts = basename.split("_vs_")
    peptide_name = parts[0] if len(parts) > 1 else "unknown"
    target_name = parts[1] if len(parts) > 1 else "unknown"

    print("=" * 60)
    print("OPENMM MD PROTOCOL")
    print("Peptide: %s | Target: %s | Duration: %.1f ns" % (peptide_name, target_name, args.production_ns))
    print("PDB: %s" % args.pdb)
    print("Output: %s" % args.output_dir)
    print("=" * 60)

    # Setup
    system, topology, positions, forcefield = setup_system(args.pdb)

    # Save solvated structure for analysis reference
    os.makedirs(args.output_dir, exist_ok=True)
    from openmm.app import PDBFile
    ref_path = os.path.join(args.output_dir, "solvated.pdb")
    with open(ref_path, "w") as f:
        PDBFile.writeFile(topology, positions, f)

    # Run simulation
    t_start = time.time()
    traj_path, wall_time = run_simulation(system, topology, positions, args.output_dir, args.production_ns)

    # Analyze
    summary = analyze_trajectory(traj_path, ref_path, args.output_dir, peptide_name, target_name)
    summary["wall_time_hours"] = round(wall_time, 2)

    # Update summary with wall time
    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    total_time = (time.time() - t_start) / 3600
    print("\nTotal time (setup + MD + analysis): %.2f hours" % total_time)


if __name__ == "__main__":
    main()
