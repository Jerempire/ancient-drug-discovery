"""MD + MMGBSA binding free energy for peptide-protein complexes.

Runs 20ns explicit-solvent MD on GPU, then computes MMGBSA binding
free energy using OBC2 implicit solvent on trajectory snapshots.

Usage: python md_mmgbsa.py <input.pdb> <label> <output.json> [--ns 20]

Output JSON contains:
  - RMSD trajectory (peptide stability)
  - Contact count trajectory
  - MMGBSA ΔG_bind (kcal/mol) with std
  - Per-frame energy decomposition
"""
import sys
import os
import json
import argparse
import time
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass


def get_chain_indices(topology):
    """Get atom indices for peptide (chain 0) and receptor (chain 1+)."""
    chains = list(topology.chains())
    if len(chains) < 2:
        raise ValueError(f"Need >= 2 chains, got {len(chains)}")

    pep_atoms = [a.index for a in chains[0].atoms()]
    rec_atoms = []
    for ch in chains[1:]:
        rec_atoms.extend(a.index for a in ch.atoms())

    pep_ca = [a.index for a in chains[0].atoms() if a.name == "CA"]
    rec_ca = []
    for ch in chains[1:]:
        rec_ca.extend(a.index for a in ch.atoms() if a.name == "CA")

    # Protein-only atoms (no water/ions) for MMGBSA
    protein_atoms = sorted(pep_atoms + rec_atoms)

    return {
        "pep_all": pep_atoms,
        "rec_all": rec_atoms,
        "pep_ca": pep_ca,
        "rec_ca": rec_ca,
        "protein": protein_atoms,
    }


def run_md(input_pdb, label, sim_ns=20):
    """Run explicit-solvent MD, return topology + trajectory frames."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
    from openmm import unit

    platform = Platform.getPlatformByName("CUDA")
    print(f"  Platform: {platform.getName()}")

    # Fix structure
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
    modeller.addSolvent(ff, padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)
    n_atoms = modeller.topology.getNumAtoms()
    print(f"  System: {n_atoms} atoms")

    indices = get_chain_indices(modeller.topology)
    print(f"  Peptide: {len(indices['pep_ca'])} res, Receptor: {len(indices['rec_ca'])} res")

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

    sim = Simulation(
        modeller.topology, system, integrator,
        platform, {"CudaPrecision": "mixed"}
    )
    sim.context.setPositions(modeller.positions)

    # Minimize
    print("  Minimizing...")
    sim.minimizeEnergy(maxIterations=2000)

    # Equilibrate 200ps
    print("  Equilibrating (200ps)...")
    sim.step(100000)

    # Get initial reference
    state = sim.context.getState(getPositions=True)
    pos0 = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
    init_pep_ca = pos0[indices["pep_ca"]].copy()

    # Production
    n_steps = int(sim_ns * 1e6 / 2)  # 2fs timestep
    snap_interval = 50000  # 100ps between snapshots
    n_snaps = n_steps // snap_interval
    print(f"  Production: {sim_ns}ns, {n_snaps} snapshots...")

    frames = []
    rmsds = []
    dists = []
    contacts_list = []
    times = []

    t_start = time.time()
    for i in range(n_snaps):
        sim.step(snap_interval)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

        # RMSD of peptide CA
        cur_pep_ca = pos[indices["pep_ca"]]
        rmsd = float(np.sqrt(np.mean(np.sum((cur_pep_ca - init_pep_ca) ** 2, axis=1))))

        # Centroid distance
        pep_cent = np.mean(cur_pep_ca, axis=0)
        rec_cent = np.mean(pos[indices["rec_ca"]], axis=0)
        dist = float(np.linalg.norm(pep_cent - rec_cent))

        # Contact count (<4.5A between any pep-rec heavy atom pair)
        pep_pos = pos[indices["pep_all"]]
        rec_pos = pos[indices["rec_all"]]
        # Sample contacts (full pairwise is expensive for large receptors)
        n_contacts = 0
        for pp in pep_pos:
            d = np.sqrt(np.sum((rec_pos - pp) ** 2, axis=1))
            n_contacts += int(np.sum(d < 4.5))

        t = (i + 1) * snap_interval * 0.002 / 1000  # ns
        rmsds.append(round(rmsd, 2))
        dists.append(round(dist, 2))
        contacts_list.append(n_contacts)
        times.append(round(t, 3))

        # Save protein-only coordinates for MMGBSA (every 500ps = every 5th frame)
        if (i + 1) % 5 == 0:
            protein_pos = pos[indices["protein"]].copy()
            frames.append(protein_pos)

        if (i + 1) % 20 == 0:
            elapsed = time.time() - t_start
            rate = t / elapsed * 60 if elapsed > 0 else 0
            print(f"    t={t:.1f}ns  RMSD={rmsd:.1f}A  dist={dist:.1f}A  contacts={n_contacts}  ({rate:.1f} ns/min)")

    md_time = time.time() - t_start
    print(f"  MD complete in {md_time/60:.1f} min")

    md_data = {
        "label": label,
        "sim_ns": sim_ns,
        "n_atoms": n_atoms,
        "pep_residues": len(indices["pep_ca"]),
        "rec_residues": len(indices["rec_ca"]),
        "time_ns": times,
        "rmsd_A": rmsds,
        "distance_A": dists,
        "contacts": contacts_list,
        "final_rmsd": rmsds[-1],
        "avg_rmsd_last_5ns": round(float(np.mean(rmsds[-50:])), 2),
        "max_rmsd": round(float(max(rmsds)), 2),
        "avg_contacts_last_5ns": round(float(np.mean(contacts_list[-50:])), 1),
        "md_time_min": round(md_time / 60, 1),
    }

    return modeller.topology, indices, frames, md_data


def compute_mmgbsa(topology, indices, frames):
    """Compute MMGBSA binding free energy from MD trajectory frames.

    Uses OBC2 (GBn2) implicit solvent for solvation energy.
    1-trajectory approach: complex, receptor, ligand energies from same frames.
    """
    import parmed
    from openmm.app import ForceField, Simulation, NoCutoff
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm import unit

    print(f"\n  Computing MMGBSA on {len(frames)} frames...")

    # Build parmed structure from topology
    # We need the protein-only topology
    all_chains = list(topology.chains())
    pep_chain = all_chains[0]
    rec_chains = all_chains[1:]

    # Get protein-only residue indices for subsetting
    pep_local = list(range(len(indices["pep_all"])))
    rec_local = list(range(len(indices["pep_all"]), len(indices["protein"])))

    # Force field with implicit solvent (OBC2/GBn2)
    ff_gb = ForceField("amber14-all.xml", "implicit/obc2.xml")

    # We need to create 3 systems: complex, receptor, peptide
    # Strategy: build from the protein-only subset of topology
    # Use PDBFixer-like approach to create clean topologies

    # Create protein-only topology
    from openmm.app import Topology, Modeller
    protein_top = Topology()
    protein_top.setPeriodicBoxVectors(None)

    # Rebuild topology with only protein atoms
    chain_map = {}
    atom_map = {}
    idx = 0
    for chain in all_chains:
        new_chain = protein_top.addChain(chain.id)
        chain_map[chain] = new_chain
        for res in chain.residues():
            # Skip water/ions
            if res.name in ("HOH", "WAT", "NA", "CL", "K", "MG"):
                continue
            new_res = protein_top.addResidue(res.name, new_chain, res.id)
            for atom in res.atoms():
                if atom.index in indices["protein"]:
                    new_atom = protein_top.addAtom(atom.name, atom.element, new_res)
                    atom_map[atom.index] = idx
                    idx += 1

    # Add bonds for protein atoms
    for bond in topology.bonds():
        a1, a2 = bond[0].index, bond[1].index
        if a1 in atom_map and a2 in atom_map:
            atoms_list = list(protein_top.atoms())
            protein_top.addBond(atoms_list[atom_map[a1]], atoms_list[atom_map[a2]])

    # Create complex system
    try:
        complex_system = ff_gb.createSystem(protein_top, nonbondedMethod=NoCutoff)
    except Exception as e:
        print(f"  WARNING: Could not create GB system: {e}")
        print("  Falling back to interaction energy only.")
        return None

    platform = Platform.getPlatformByName("CUDA")
    integrator_c = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    sim_complex = Simulation(protein_top, complex_system, integrator_c, platform, {"CudaPrecision": "mixed"})

    # Now we need receptor-only and peptide-only systems
    # Build receptor topology
    rec_top = Topology()
    rec_atoms_list = []
    r_atom_map = {}
    r_idx = 0
    for chain in rec_chains:
        new_chain = rec_top.addChain(chain.id)
        for res in chain.residues():
            if res.name in ("HOH", "WAT", "NA", "CL", "K", "MG"):
                continue
            new_res = rec_top.addResidue(res.name, new_chain, res.id)
            for atom in res.atoms():
                if atom.index in indices["rec_all"]:
                    new_atom = rec_top.addAtom(atom.name, atom.element, new_res)
                    r_atom_map[atom.index] = r_idx
                    r_idx += 1

    for bond in topology.bonds():
        a1, a2 = bond[0].index, bond[1].index
        if a1 in r_atom_map and a2 in r_atom_map:
            atoms_l = list(rec_top.atoms())
            rec_top.addBond(atoms_l[r_atom_map[a1]], atoms_l[r_atom_map[a2]])

    rec_system = ff_gb.createSystem(rec_top, nonbondedMethod=NoCutoff)
    integrator_r = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    sim_rec = Simulation(rec_top, rec_system, integrator_r, platform, {"CudaPrecision": "mixed"})

    # Build peptide topology
    pep_top = Topology()
    p_atom_map = {}
    p_idx = 0
    pep_chain_obj = all_chains[0]
    new_chain = pep_top.addChain(pep_chain_obj.id)
    for res in pep_chain_obj.residues():
        if res.name in ("HOH", "WAT", "NA", "CL", "K", "MG"):
            continue
        new_res = pep_top.addResidue(res.name, new_chain, res.id)
        for atom in res.atoms():
            if atom.index in indices["pep_all"]:
                new_atom = pep_top.addAtom(atom.name, atom.element, new_res)
                p_atom_map[atom.index] = p_idx
                p_idx += 1

    for bond in topology.bonds():
        a1, a2 = bond[0].index, bond[1].index
        if a1 in p_atom_map and a2 in p_atom_map:
            atoms_l = list(pep_top.atoms())
            pep_top.addBond(atoms_l[p_atom_map[a1]], atoms_l[p_atom_map[a2]])

    pep_system = ff_gb.createSystem(pep_top, nonbondedMethod=NoCutoff)
    integrator_p = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    sim_pep = Simulation(pep_top, pep_system, integrator_p, platform, {"CudaPrecision": "mixed"})

    # Map protein-only indices back to local subsystem indices
    protein_list = sorted(indices["protein"])
    pep_set = set(indices["pep_all"])
    rec_set = set(indices["rec_all"])

    pep_frame_indices = [i for i, orig in enumerate(protein_list) if orig in pep_set]
    rec_frame_indices = [i for i, orig in enumerate(protein_list) if orig in rec_set]

    # Compute energies per frame
    dg_values = []
    kj_to_kcal = 0.239006

    for fi, frame_pos in enumerate(frames):
        pos_nm = frame_pos * 0.1  # Angstrom -> nm
        from openmm import Vec3
        complex_positions = [Vec3(pos_nm[i][0], pos_nm[i][1], pos_nm[i][2]) for i in range(len(pos_nm))]

        # Complex energy
        sim_complex.context.setPositions(complex_positions * unit.nanometers)
        e_complex = sim_complex.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

        # Receptor energy
        rec_positions = [complex_positions[i] for i in rec_frame_indices]
        sim_rec.context.setPositions(rec_positions * unit.nanometers)
        e_rec = sim_rec.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

        # Peptide energy
        pep_positions = [complex_positions[i] for i in pep_frame_indices]
        sim_pep.context.setPositions(pep_positions * unit.nanometers)
        e_pep = sim_pep.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

        dg = (e_complex - e_rec - e_pep) * kj_to_kcal
        dg_values.append(dg)

        if (fi + 1) % 10 == 0:
            print(f"    Frame {fi+1}/{len(frames)}: dG={dg:.1f} kcal/mol")

    dg_arr = np.array(dg_values)
    # Discard first 20% as equilibration
    equil_cutoff = max(1, len(dg_arr) // 5)
    dg_prod = dg_arr[equil_cutoff:]

    result = {
        "dg_bind_mean_kcal": round(float(np.mean(dg_prod)), 2),
        "dg_bind_std_kcal": round(float(np.std(dg_prod)), 2),
        "dg_bind_median_kcal": round(float(np.median(dg_prod)), 2),
        "n_frames_total": len(dg_values),
        "n_frames_production": len(dg_prod),
        "dg_per_frame": [round(float(x), 2) for x in dg_values],
    }
    print(f"  MMGBSA: dG_bind = {result['dg_bind_mean_kcal']:.1f} +/- {result['dg_bind_std_kcal']:.1f} kcal/mol")
    return result


def main():
    parser = argparse.ArgumentParser(description="MD + MMGBSA validation")
    parser.add_argument("input_pdb", help="Input PDB (peptide=chain A, receptor=chain B)")
    parser.add_argument("label", help="Label for this run (e.g., VAGSAF_vs_erap2k392)")
    parser.add_argument("output_json", help="Output JSON path")
    parser.add_argument("--ns", type=float, default=20, help="Simulation length in ns (default: 20)")
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"  MD + MMGBSA: {args.label}")
    print(f"  Input: {args.input_pdb}")
    print(f"  Duration: {args.ns}ns")
    print(f"{'='*60}\n")

    t0 = time.time()

    # Step 1: MD
    topology, indices, frames, md_data = run_md(args.input_pdb, args.label, sim_ns=args.ns)

    # Step 2: MMGBSA
    mmgbsa_data = compute_mmgbsa(topology, indices, frames)

    # Combine results
    result = {
        "label": args.label,
        "input_pdb": os.path.basename(args.input_pdb),
        "md": md_data,
        "mmgbsa": mmgbsa_data,
        "total_time_min": round((time.time() - t0) / 60, 1),
    }

    # Verdict
    if mmgbsa_data:
        dg = mmgbsa_data["dg_bind_mean_kcal"]
        if dg < -15:
            result["binding_verdict"] = "STRONG"
        elif dg < -5:
            result["binding_verdict"] = "MODERATE"
        elif dg < 0:
            result["binding_verdict"] = "WEAK"
        else:
            result["binding_verdict"] = "NO_BINDING"
    else:
        result["binding_verdict"] = "MMGBSA_FAILED"

    rmsd = md_data["avg_rmsd_last_5ns"]
    if rmsd < 3.0:
        result["retention_verdict"] = "STABLE"
    elif rmsd < 5.0:
        result["retention_verdict"] = "MODERATE"
    else:
        result["retention_verdict"] = "UNSTABLE"

    with open(args.output_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved: {args.output_json}")
    print(f"  Binding: {result['binding_verdict']}, Retention: {result['retention_verdict']}")
    print(f"  Total time: {result['total_time_min']:.1f} min")

    return result


if __name__ == "__main__":
    main()
