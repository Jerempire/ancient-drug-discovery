"""OpenMM MD: Place peptides near ERAP2 channel and simulate.
Uses AlphaFold PDB as receptor. Builds peptide, places near res 392.
No Boltz-2 needed.

2 peptides x 2 variants = 4 sims, 10ns each.
"""
import sys, os, json, time
import numpy as np

try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer
    print("OpenMM + PDBFixer OK")
except ImportError as e:
    print(f"Missing: {e}")
    sys.exit(1)

try:
    import mdtraj
    HAS_MDTRAJ = True
    print("mdtraj OK")
except Exception:
    HAS_MDTRAJ = False
    print("mdtraj not available")

RESULTS_DIR = "/workspace/md/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

TEMPERATURE = 300 * unit.kelvin
PRESSURE = 1 * unit.atmosphere
TIMESTEP = 2 * unit.femtoseconds
EQUIL_STEPS = 25000       # 50 ps
PROD_STEPS = 5000000      # 10 ns
REPORT_INTERVAL = 2500    # every 5 ps

PEPTIDES = {
    "V_pep": "VKLLLLSI",
    "D_pep": "DKLLLLSI",
}


def run_md(receptor_pdb, peptide_seq, variant_name, case_name):
    print(f"\n{'='*60}")
    print(f"  MD: {case_name} ({peptide_seq} + ERAP2 {variant_name})")
    print(f"{'='*60}")
    t_start = time.time()

    print("  Loading receptor...")
    fixer = PDBFixer(receptor_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    if "N392" in variant_name:
        print("  Mutating K392 -> N392...")
        fixer.applyMutations(["LYS-392-ASN"], "A")
        # Re-fix after mutation
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

    # Find residue 392 CA position
    target_pos = None
    for atom in fixer.topology.atoms():
        if atom.residue.index == 391 and atom.name == "CA":
            target_pos = fixer.positions[atom.index]
            break

    if target_pos is None:
        print("  WARNING: residue 392 CA not found, using center")
        coords = np.array([[p.x, p.y, p.z] for p in fixer.positions])
        center = np.mean(coords, axis=0)
        target_pos = openmm.Vec3(center[0], center[1], center[2]) * unit.nanometers

    print(f"  Res 392 CA: ({target_pos.x:.2f}, {target_pos.y:.2f}, {target_pos.z:.2f}) nm")

    # Build peptide as extended chain using PDBFixer
    aa_map = {
        "A": "ALA", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
        "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
        "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR",
        "V": "VAL", "W": "TRP", "Y": "TYR",
    }

    print("  Building peptide...")
    pdb_lines = ["HEADER PEPTIDE\n"]
    atom_idx = 1
    for i, aa in enumerate(peptide_seq):
        resname = aa_map[aa]
        x = i * 3.8  # 3.8A per residue in extended chain
        # N
        pdb_lines.append(
            f"ATOM  {atom_idx:5d}  N   {resname} B{i+1:4d}    "
            f"{x - 1.2:8.3f}{0.7:8.3f}{0.0:8.3f}  1.00  0.00           N\n"
        )
        atom_idx += 1
        # CA
        pdb_lines.append(
            f"ATOM  {atom_idx:5d}  CA  {resname} B{i+1:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
        )
        atom_idx += 1
        # C
        pdb_lines.append(
            f"ATOM  {atom_idx:5d}  C   {resname} B{i+1:4d}    "
            f"{x + 1.2:8.3f}{-0.7:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
        )
        atom_idx += 1
        # O
        pdb_lines.append(
            f"ATOM  {atom_idx:5d}  O   {resname} B{i+1:4d}    "
            f"{x + 1.2:8.3f}{-1.9:8.3f}{0.0:8.3f}  1.00  0.00           O\n"
        )
        atom_idx += 1
    pdb_lines.append("TER\n")
    pdb_lines.append("END\n")

    pep_tmp = f"/workspace/md/structures/{case_name}_peptide_tmp.pdb"
    with open(pep_tmp, "w") as f:
        f.writelines(pdb_lines)

    pep_fixer = PDBFixer(pep_tmp)
    pep_fixer.findMissingResidues()
    pep_fixer.findMissingAtoms()
    pep_fixer.addMissingAtoms()
    pep_fixer.addMissingHydrogens(7.4)

    # Place peptide 2nm from res 392
    pep_positions = []
    offset = openmm.Vec3(2.0, 0.5, 0.0) * unit.nanometers
    for pos in pep_fixer.positions:
        new_pos = pos + target_pos + offset
        pep_positions.append(new_pos)

    # Combine
    print("  Combining receptor + peptide...")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.add(pep_fixer.topology, pep_positions)

    # Solvate
    print("  Solvating...")
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers,
                        ionicStrength=0.15 * unit.molar)

    n_atoms = modeller.topology.getNumAtoms()
    print(f"  System: {n_atoms} atoms")

    # Build system
    print("  Building system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )
    system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

    integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, 1 / unit.picosecond, TIMESTEP)
    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision": "mixed"}
        print("  Using CUDA platform")
    except Exception:
        platform = openmm.Platform.getPlatformByName("CPU")
        properties = {}
        print(f"  Using CPU platform ({os.cpu_count()} cores)")

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print("  Minimizing...")
    simulation.minimizeEnergy(maxIterations=2000)

    # Equilibrate
    print(f"  Equilibrating ({EQUIL_STEPS * 2 / 1000:.0f} ps)...")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    simulation.step(EQUIL_STEPS)

    # Save start
    out_pdb = f"{RESULTS_DIR}/{case_name}_start.pdb"
    state = simulation.context.getState(getPositions=True)
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Production
    out_dcd = f"{RESULTS_DIR}/{case_name}_traj.dcd"
    out_log = f"{RESULTS_DIR}/{case_name}_energy.csv"

    simulation.reporters.append(app.DCDReporter(out_dcd, REPORT_INTERVAL))
    simulation.reporters.append(app.StateDataReporter(
        out_log, REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, temperature=True,
        speed=True, separator=","
    ))

    print(f"  Production ({PROD_STEPS * 2 / 1e6:.0f} ns)...")
    prod_start = time.time()
    simulation.step(PROD_STEPS)
    prod_time = time.time() - prod_start
    total_time = time.time() - t_start

    ns_per_day = (PROD_STEPS * 2e-6) / (prod_time / 86400)
    print(f"  Done: {total_time:.0f}s total, {ns_per_day:.0f} ns/day")

    result = {
        "case": case_name,
        "peptide": peptide_seq,
        "variant": variant_name,
        "n_atoms": n_atoms,
        "prod_ns": PROD_STEPS * 2e-6,
        "wall_time_s": total_time,
        "ns_per_day": ns_per_day,
    }

    # Analysis
    if HAS_MDTRAJ and os.path.exists(out_dcd):
        print("  Analyzing trajectory...")
        try:
            traj = mdtraj.load(out_dcd, top=out_pdb)
            chains = list(traj.topology.chains)
            protein_chains = [c for c in chains if c.n_residues > 5]

            if len(protein_chains) >= 2:
                rec_ca = traj.topology.select("chainid 0 and name CA")
                pep_ca = traj.topology.select(f"chainid {protein_chains[1].index} and name CA")

                if len(pep_ca) > 0 and len(rec_ca) > 0:
                    traj_aligned = traj.superpose(traj, atom_indices=rec_ca)
                    pep_rmsd = mdtraj.rmsd(traj_aligned, traj_aligned, atom_indices=pep_ca) * 10
                    result["peptide_rmsd_mean_A"] = float(np.mean(pep_rmsd))
                    result["peptide_rmsd_std_A"] = float(np.std(pep_rmsd))
                    result["peptide_rmsd_max_A"] = float(np.max(pep_rmsd))

                    # Min distance to channel region (res 385-415)
                    rec_channel = traj.topology.select("chainid 0 and resid 385 to 415 and name CA")
                    if len(rec_channel) > 0:
                        pairs = [(p, r) for p in pep_ca for r in rec_channel]
                        dists = mdtraj.compute_distances(traj, pairs)
                        min_per_frame = np.min(dists, axis=1) * 10
                        result["min_dist_mean_A"] = float(np.mean(min_per_frame))
                        result["min_dist_max_A"] = float(np.max(min_per_frame))

                        dissoc = np.where(min_per_frame > 10.0)[0]
                        result["first_dissociation_ns"] = float(dissoc[0]) * REPORT_INTERVAL * 2e-6 if len(dissoc) > 0 else None

                    print(f"    RMSD: {result.get('peptide_rmsd_mean_A', 0):.1f} +/- {result.get('peptide_rmsd_std_A', 0):.1f} A")
                    print(f"    Min dist: {result.get('min_dist_mean_A', 0):.1f} A")
                    d = result.get("first_dissociation_ns")
                    print(f"    {'STABLE' if d is None else f'DISSOCIATED at {d:.1f} ns'}")
        except Exception as e:
            print(f"    Analysis error: {e}")

    return result


if __name__ == "__main__":
    receptor = "/workspace/md/structures/erap2_wt_alphafold.pdb"
    all_results = {}

    for pep_name, pep_seq in PEPTIDES.items():
        for variant in ["K392", "N392"]:
            case = f"{pep_name}_{variant}"
            try:
                result = run_md(receptor, pep_seq, variant, case)
                all_results[case] = result
            except Exception as e:
                print(f"\nFAILED {case}: {e}")
                import traceback
                traceback.print_exc()
                all_results[case] = {"case": case, "error": str(e)}

    # Summary
    print(f"\n{'='*70}")
    print(f"  MD SUMMARY")
    print(f"{'='*70}")
    print(f"{'Case':>15s}  {'Pep':>10s}  {'Var':>5s}  {'RMSD(A)':>12s}  {'Dist(A)':>8s}  {'Dissoc':>10s}")
    print("-" * 65)
    for name, r in all_results.items():
        if "error" in r:
            print(f"{name:>15s}  ERROR: {r['error'][:40]}")
            continue
        rmsd = f"{r.get('peptide_rmsd_mean_A', 0):.1f}+/-{r.get('peptide_rmsd_std_A', 0):.1f}" if "peptide_rmsd_mean_A" in r else "N/A"
        mdist = f"{r.get('min_dist_mean_A', 0):.1f}" if "min_dist_mean_A" in r else "N/A"
        d = r.get("first_dissociation_ns")
        dissoc = "STABLE" if d is None else f"{d:.1f}ns"
        print(f"{name:>15s}  {r['peptide']:>10s}  {r['variant']:>5s}  {rmsd:>12s}  {mdist:>8s}  {dissoc:>10s}")

    with open(f"{RESULTS_DIR}/md_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {RESULTS_DIR}/md_results.json")
