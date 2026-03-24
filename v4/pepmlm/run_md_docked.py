"""Corrected OpenMM MD from docked channel-bound poses.

Starts from Boltz-2 predicted complexes (peptide already in substrate channel).
4 systems x 3 replicates x 10ns = 120ns total.

Reports: RMSD, P1-392 distance, channel contacts, H-bonds/salt bridges,
         bound fraction, variant-dependent differences.
"""
import sys, os, json, time, glob
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import openmm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

try:
    import mdtraj
    HAS_MDTRAJ = True
except Exception:
    HAS_MDTRAJ = False
    print("WARNING: mdtraj not available")

RESULTS_DIR = "/workspace/md_docked/results"
STRUCT_DIR = "/workspace/md_docked/structures"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(STRUCT_DIR, exist_ok=True)

# ── MD Parameters ─────────────────────────────────────────────────────────────
TEMPERATURE = 300 * unit.kelvin
PRESSURE = 1 * unit.atmosphere
TIMESTEP = 2 * unit.femtoseconds
EQUIL_STEPS = 50000       # 100 ps
PROD_STEPS = 5000000      # 10 ns
REPORT_INTERVAL = 2500    # every 5 ps
N_REPLICATES = 3

# Channel residues (0-indexed) for contact analysis
CHANNEL_RES_START = 384   # res 385 in 1-indexed
CHANNEL_RES_END = 414     # res 415 in 1-indexed
RES_392_IDX = 391         # 0-indexed

SYSTEMS = {
    "V_K392": {"peptide": "VKLLLLSI", "p1": "V", "variant": "K392"},
    "V_N392": {"peptide": "VKLLLLSI", "p1": "V", "variant": "N392"},
    "D_K392": {"peptide": "DKLLLLSI", "p1": "D", "variant": "K392"},
    "D_N392": {"peptide": "DKLLLLSI", "p1": "D", "variant": "N392"},
}

ERAP2_K392 = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPS"
    "VVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYP"
    "AHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLF"
    "KANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASP"
    "DKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLW"
    "VTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPI"
    "SKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTS"
    "GGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQER"
    "YLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYF"
    "KPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAG"
    "WNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRE"
    "NWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLR"
    "TWLMVNT"
)
ERAP2_N392 = ERAP2_K392[:391] + "N" + ERAP2_K392[392:]


def generate_boltz2_structures():
    """Generate docked structures with Boltz-2."""
    print("=" * 60)
    print("  Generating Boltz-2 docked structures")
    print("=" * 60)

    for sys_name, info in SYSTEMS.items():
        erap_seq = ERAP2_K392 if info["variant"] == "K392" else ERAP2_N392
        yaml_path = f"{STRUCT_DIR}/{sys_name}.yaml"
        with open(yaml_path, "w") as f:
            f.write(f"sequences:\n")
            f.write(f"  - protein:\n")
            f.write(f"      id: A\n")
            f.write(f"      sequence: {erap_seq}\n")
            f.write(f"      msa: empty\n")
            f.write(f"  - protein:\n")
            f.write(f"      id: B\n")
            f.write(f"      sequence: {info['peptide']}\n")
            f.write(f"      msa: empty\n")

        print(f"\n  Boltz-2: {sys_name} ({info['peptide']} + ERAP2 {info['variant']})...")
        out_dir = f"{STRUCT_DIR}/{sys_name}_boltz"
        ret = os.system(
            f"boltz predict {yaml_path} "
            f"--diffusion_samples 1 --seed 42 "
            f"--out_dir {out_dir} 2>&1 | tail -3"
        )

        # Find CIF
        cifs = glob.glob(f"{out_dir}/**/*.cif", recursive=True)
        if cifs:
            os.system(f"cp {cifs[0]} {STRUCT_DIR}/{sys_name}.cif")
            print(f"  Structure saved: {sys_name}.cif")
        else:
            print(f"  ERROR: No CIF for {sys_name}")


def cif_to_fixed_pdb(cif_path, pdb_path, mutate_n392=False):
    """Convert CIF to fixed PDB with all atoms."""
    fixer = PDBFixer(cif_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    with open(pdb_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return pdb_path


def run_single_md(pdb_path, case_name, seed):
    """Run a single MD replicate."""
    print(f"\n  --- {case_name} (seed {seed}) ---")
    t_start = time.time()

    # Load structure
    pdb = app.PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)

    # Solvate
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers,
                        ionicStrength=0.15 * unit.molar)

    n_atoms = modeller.topology.getNumAtoms()
    print(f"    {n_atoms} atoms")

    # System
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )
    system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

    integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, 1 / unit.picosecond, TIMESTEP)
    integrator.setRandomNumberSeed(seed)

    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision": "mixed"}
        print("    CUDA platform")
    except Exception:
        platform = openmm.Platform.getPlatformByName("CPU")
        properties = {}
        print("    CPU platform (slow!)")

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print("    Minimizing...")
    simulation.minimizeEnergy(maxIterations=2000)

    # Equilibrate
    print(f"    Equilibrating ({EQUIL_STEPS * 2 / 1000:.0f} ps)...")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, seed)
    simulation.step(EQUIL_STEPS)

    # Save equilibrated structure
    prefix = f"{RESULTS_DIR}/{case_name}_s{seed}"
    state = simulation.context.getState(getPositions=True)
    with open(f"{prefix}_start.pdb", "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Production reporters
    simulation.reporters.append(app.DCDReporter(f"{prefix}_traj.dcd", REPORT_INTERVAL))
    simulation.reporters.append(app.StateDataReporter(
        f"{prefix}_energy.csv", REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, temperature=True,
        speed=True, separator=","
    ))

    print(f"    Production ({PROD_STEPS * 2 / 1e6:.0f} ns)...")
    prod_start = time.time()
    simulation.step(PROD_STEPS)
    prod_time = time.time() - prod_start
    total_time = time.time() - t_start
    ns_per_day = (PROD_STEPS * 2e-6) / (prod_time / 86400)
    print(f"    Done: {total_time:.0f}s, {ns_per_day:.0f} ns/day")

    return {
        "case": case_name, "seed": seed, "n_atoms": n_atoms,
        "prod_ns": PROD_STEPS * 2e-6,
        "wall_time_s": total_time, "ns_per_day": ns_per_day,
    }


def analyze_replicate(case_name, seed, sys_info):
    """Analyze a single MD trajectory."""
    prefix = f"{RESULTS_DIR}/{case_name}_s{seed}"
    pdb_path = f"{prefix}_start.pdb"
    dcd_path = f"{prefix}_traj.dcd"

    if not os.path.exists(dcd_path) or not HAS_MDTRAJ:
        return {}

    print(f"    Analyzing {case_name} seed {seed}...")
    traj = mdtraj.load(dcd_path, top=pdb_path)
    topology = traj.topology
    result = {}

    # Find chains
    chains = list(topology.chains)
    protein_chains = [c for c in chains if c.n_residues > 3]
    if len(protein_chains) < 2:
        print("    WARNING: <2 protein chains found")
        return result

    receptor_chain = protein_chains[0]
    peptide_chain = protein_chains[1]
    pep_chain_idx = peptide_chain.index

    # Atom selections
    rec_ca = topology.select(f"chainid {receptor_chain.index} and name CA")
    pep_ca = topology.select(f"chainid {pep_chain_idx} and name CA")
    pep_heavy = topology.select(f"chainid {pep_chain_idx} and not element H")

    if len(pep_ca) == 0 or len(rec_ca) == 0:
        return result

    # 1. Peptide RMSD (aligned to receptor backbone)
    traj_aligned = traj.superpose(traj, atom_indices=rec_ca)
    pep_rmsd = mdtraj.rmsd(traj_aligned, traj_aligned, atom_indices=pep_ca) * 10  # A
    result["peptide_rmsd_mean"] = float(np.mean(pep_rmsd))
    result["peptide_rmsd_std"] = float(np.std(pep_rmsd))
    result["peptide_rmsd_max"] = float(np.max(pep_rmsd))

    # 2. P1 to residue 392 distance
    # P1 = first residue of peptide chain
    p1_atoms = topology.select(f"chainid {pep_chain_idx} and resid 0 and name CA")
    # Res 392 of receptor (0-indexed = 391)
    res392_atoms = topology.select(f"chainid {receptor_chain.index} and resid {RES_392_IDX} and name CA")

    if len(p1_atoms) > 0 and len(res392_atoms) > 0:
        p1_392_pairs = [(p1_atoms[0], res392_atoms[0])]
        p1_392_dist = mdtraj.compute_distances(traj, p1_392_pairs)[:, 0] * 10  # A
        result["p1_392_dist_mean"] = float(np.mean(p1_392_dist))
        result["p1_392_dist_std"] = float(np.std(p1_392_dist))
        result["p1_392_dist_min"] = float(np.min(p1_392_dist))
        result["p1_392_dist_max"] = float(np.max(p1_392_dist))

    # 3. Contacts with channel residues 385-415
    channel_atoms = topology.select(
        f"chainid {receptor_chain.index} and "
        f"resid {CHANNEL_RES_START} to {CHANNEL_RES_END} and name CA"
    )
    if len(channel_atoms) > 0 and len(pep_ca) > 0:
        contact_pairs = [(p, c) for p in pep_ca for c in channel_atoms]
        contact_dists = mdtraj.compute_distances(traj, contact_pairs)
        # Contact = <5A (0.5 nm)
        contacts_per_frame = np.sum(contact_dists < 0.5, axis=1)
        result["channel_contacts_mean"] = float(np.mean(contacts_per_frame))
        result["channel_contacts_std"] = float(np.std(contacts_per_frame))
        result["channel_contacts_min"] = int(np.min(contacts_per_frame))

        # Bound fraction (>= 3 contacts)
        bound_frames = np.sum(contacts_per_frame >= 3)
        result["bound_fraction"] = float(bound_frames / traj.n_frames)

        # Time to first unbinding (contacts < 2)
        unbound = np.where(contacts_per_frame < 2)[0]
        if len(unbound) > 0:
            result["first_unbinding_ns"] = float(unbound[0]) * REPORT_INTERVAL * 2e-6
        else:
            result["first_unbinding_ns"] = None  # never unbound

    # 4. Hydrogen bonds between peptide and channel
    pep_residues = list(peptide_chain.residues)
    channel_residues = [r for r in receptor_chain.residues
                        if CHANNEL_RES_START <= r.index <= CHANNEL_RES_END]

    try:
        hbonds = mdtraj.baker_hubbard(traj, freq=0.1)  # present >10% of time
        # Filter for peptide-channel H-bonds
        pep_atom_set = set(topology.select(f"chainid {pep_chain_idx}"))
        channel_atom_set = set()
        for r in channel_residues:
            for a in r.atoms:
                channel_atom_set.add(a.index)

        pc_hbonds = 0
        for hb in hbonds:
            donor, _, acceptor = hb
            if (donor in pep_atom_set and acceptor in channel_atom_set) or \
               (donor in channel_atom_set and acceptor in pep_atom_set):
                pc_hbonds += 1
        result["peptide_channel_hbonds"] = pc_hbonds
    except Exception as e:
        result["peptide_channel_hbonds"] = f"error: {e}"

    # 5. Salt bridge: P1 charged atoms to res 392 sidechain
    p1_res = list(peptide_chain.residues)[0]
    res392 = list(receptor_chain.residues)[RES_392_IDX]

    # Get charged/polar heavy atoms
    p1_charged = [a.index for a in p1_res.atoms
                  if a.element.symbol in ("O", "N") and a.name not in ("N", "C", "O", "CA")]
    r392_charged = [a.index for a in res392.atoms
                    if a.element.symbol in ("O", "N") and a.name not in ("N", "C", "O", "CA")]

    if p1_charged and r392_charged:
        sb_pairs = [(p, r) for p in p1_charged for r in r392_charged]
        sb_dists = mdtraj.compute_distances(traj, sb_pairs)
        min_sb_dist = np.min(sb_dists, axis=1) * 10  # A
        result["salt_bridge_dist_mean"] = float(np.mean(min_sb_dist))
        result["salt_bridge_dist_min"] = float(np.min(min_sb_dist))
        # Salt bridge occupancy (<4A)
        sb_occupancy = float(np.mean(min_sb_dist < 4.0))
        result["salt_bridge_occupancy"] = sb_occupancy
    else:
        result["salt_bridge_dist_mean"] = None
        result["salt_bridge_occupancy"] = 0.0

    return result


def main():
    # Step 1: Generate Boltz-2 structures
    generate_boltz2_structures()

    # Verify structures exist
    for sys_name in SYSTEMS:
        cif_path = f"{STRUCT_DIR}/{sys_name}.cif"
        if not os.path.exists(cif_path):
            print(f"FATAL: No structure for {sys_name}. Cannot proceed.")
            return

    # Step 2: Fix structures
    print("\n" + "=" * 60)
    print("  Preparing PDB files")
    print("=" * 60)
    for sys_name in SYSTEMS:
        cif_path = f"{STRUCT_DIR}/{sys_name}.cif"
        pdb_path = f"{STRUCT_DIR}/{sys_name}.pdb"
        print(f"  Fixing {sys_name}...")
        cif_to_fixed_pdb(cif_path, pdb_path)
        n = sum(1 for _ in app.PDBFile(pdb_path).topology.atoms())
        print(f"    {n} atoms")

    # Step 3: Run MD
    print("\n" + "=" * 60)
    print(f"  Running MD: {len(SYSTEMS)} systems x {N_REPLICATES} replicates")
    print("=" * 60)

    all_results = {}
    seeds = [42, 123, 7]

    for sys_name, info in SYSTEMS.items():
        pdb_path = f"{STRUCT_DIR}/{sys_name}.pdb"
        for rep_idx in range(N_REPLICATES):
            seed = seeds[rep_idx]
            case = f"{sys_name}_rep{rep_idx}"
            try:
                md_result = run_single_md(pdb_path, case, seed)
                # Analyze
                analysis = analyze_replicate(case, seed, info)
                md_result.update(analysis)
                md_result["system"] = sys_name
                md_result["peptide"] = info["peptide"]
                md_result["p1"] = info["p1"]
                md_result["variant"] = info["variant"]
                all_results[case] = md_result
            except Exception as e:
                print(f"    FAILED: {e}")
                import traceback
                traceback.print_exc()
                all_results[case] = {"case": case, "error": str(e)}

    # Step 4: Summary
    print("\n" + "=" * 70)
    print("  MD SUMMARY (per-replicate)")
    print("=" * 70)
    print(f"{'Case':>20s}  {'RMSD':>8s}  {'P1-392':>8s}  {'Contacts':>9s}  {'Bound%':>7s}  {'SB occ':>7s}")
    print("-" * 65)
    for name, r in sorted(all_results.items()):
        if "error" in r:
            print(f"{name:>20s}  ERROR: {r['error'][:40]}")
            continue
        rmsd = f"{r.get('peptide_rmsd_mean', 0):.1f}" if "peptide_rmsd_mean" in r else "N/A"
        p1d = f"{r.get('p1_392_dist_mean', 0):.1f}" if "p1_392_dist_mean" in r else "N/A"
        cont = f"{r.get('channel_contacts_mean', 0):.1f}" if "channel_contacts_mean" in r else "N/A"
        bound = f"{r.get('bound_fraction', 0)*100:.0f}%" if "bound_fraction" in r else "N/A"
        sb = f"{r.get('salt_bridge_occupancy', 0)*100:.0f}%" if r.get("salt_bridge_occupancy") is not None else "N/A"
        print(f"{name:>20s}  {rmsd:>8s}  {p1d:>8s}  {cont:>9s}  {bound:>7s}  {sb:>7s}")

    # Aggregate by system (average across replicates)
    print("\n" + "=" * 70)
    print("  AGGREGATED (mean across replicates)")
    print("=" * 70)
    print(f"{'System':>10s}  {'Pep':>8s}  {'Var':>5s}  {'RMSD':>8s}  {'P1-392':>8s}  {'Contacts':>9s}  {'Bound%':>7s}  {'SB occ':>7s}")
    print("-" * 75)

    for sys_name, info in SYSTEMS.items():
        reps = [r for k, r in all_results.items() if r.get("system") == sys_name and "error" not in r]
        if not reps:
            print(f"{sys_name:>10s}  NO DATA")
            continue

        def avg(key):
            vals = [r[key] for r in reps if key in r and r[key] is not None]
            return np.mean(vals) if vals else None

        rmsd = avg("peptide_rmsd_mean")
        p1d = avg("p1_392_dist_mean")
        cont = avg("channel_contacts_mean")
        bound = avg("bound_fraction")
        sb = avg("salt_bridge_occupancy")

        r_str = f"{rmsd:.1f}" if rmsd is not None else "N/A"
        p_str = f"{p1d:.1f}" if p1d is not None else "N/A"
        c_str = f"{cont:.1f}" if cont is not None else "N/A"
        b_str = f"{bound*100:.0f}%" if bound is not None else "N/A"
        s_str = f"{sb*100:.0f}%" if sb is not None else "N/A"
        print(f"{sys_name:>10s}  {info['peptide']:>8s}  {info['variant']:>5s}  "
              f"{r_str:>8s}  {p_str:>8s}  {c_str:>9s}  {b_str:>7s}  {s_str:>7s}")

    # Save
    with open(f"{RESULTS_DIR}/md_docked_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {RESULTS_DIR}/md_docked_results.json")


if __name__ == "__main__":
    main()
