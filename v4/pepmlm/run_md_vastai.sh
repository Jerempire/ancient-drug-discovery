#!/bin/bash
# OpenMM MD on 2 SAR peptides x 2 ERAP2 variants = 4 simulations
# Peptides: VKLLLLSI (V at P1, K392-sel) vs DKLLLLSI (D at P1, N392-sel)
# Each: minimize → equilibrate (100ps) → production (20ns)
# ~2-3 hours on A100 for all 4
set -e

echo "============================================"
echo "  OpenMM MD: SAR Validation — $(date)"
echo "  2 peptides x 2 variants = 4 sims"
echo "============================================"

# ── Install deps ──────────────────────────────────────────────────────────────
echo "=== Installing dependencies ==="
pip install -q openmm pdbfixer biopython mdtraj 2>&1 | tail -3
pip install -q boltz 2>&1 | tail -3
pip uninstall -y cuequivariance-torch cuequivariance-ops-torch-cu12 2>/dev/null || true
echo "Dependencies ready"

# ── Step 1: Generate starting structures with Boltz-2 ─────────────────────────
echo ""
echo "=== Generating Boltz-2 starting structures ==="
mkdir -p /workspace/md/structures /workspace/md/results

# Write ERAP2 sequences
python3 << 'SEQPY'
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

cases = {
    "V_K392": ("VKLLLLSI", ERAP2_K392),
    "V_N392": ("VKLLLLSI", ERAP2_N392),
    "D_K392": ("DKLLLLSI", ERAP2_K392),
    "D_N392": ("DKLLLLSI", ERAP2_N392),
}

for name, (pep, erap) in cases.items():
    with open(f"/workspace/md/structures/{name}.yaml", "w") as f:
        f.write(f"sequences:\n")
        f.write(f"  - protein:\n")
        f.write(f"      id: A\n")
        f.write(f"      sequence: {erap}\n")
        f.write(f"      msa: empty\n")
        f.write(f"  - protein:\n")
        f.write(f"      id: B\n")
        f.write(f"      sequence: {pep}\n")
        f.write(f"      msa: empty\n")
    print(f"Created {name}.yaml: {pep} + ERAP2 {'K392' if 'K392' in name else 'N392'}")
SEQPY

# Run Boltz-2 for each (single sample is fine for MD starting structure)
for yaml in /workspace/md/structures/*.yaml; do
    name=$(basename "$yaml" .yaml)
    echo "Boltz-2: $name ..."
    boltz predict "$yaml" \
        --diffusion_samples 1 \
        --seed 42 \
        --out_dir "/workspace/md/structures/${name}_boltz" \
        --no_kernels 2>&1 | tail -3

    # Find and copy CIF
    cif=$(find "/workspace/md/structures/${name}_boltz" -name "*.cif" | head -1)
    if [ -n "$cif" ]; then
        cp "$cif" "/workspace/md/structures/${name}.cif"
        echo "  Structure: ${name}.cif"
    else
        echo "  WARNING: No CIF found for $name"
    fi
done

# ── Step 2: Run OpenMM MD ────────────────────────────────────────────────────
echo ""
echo "=== Running OpenMM MD ==="

python3 << 'MDPY'
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import os
import json
import time
import numpy as np
from pathlib import Path

try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer
except ImportError:
    print("ERROR: openmm or pdbfixer not installed")
    sys.exit(1)

try:
    import mdtraj
    HAS_MDTRAJ = True
except ImportError:
    HAS_MDTRAJ = False
    print("WARNING: mdtraj not available, skipping trajectory analysis")

STRUCT_DIR = Path("/workspace/md/structures")
RESULTS_DIR = Path("/workspace/md/results")
RESULTS_DIR.mkdir(exist_ok=True)

# MD parameters
TEMPERATURE = 300 * unit.kelvin
PRESSURE = 1 * unit.atmosphere
FRICTION = 1 / unit.picosecond
TIMESTEP = 2 * unit.femtoseconds
EQUIL_STEPS = 50000       # 100 ps equilibration
PROD_STEPS = 10000000     # 20 ns production
REPORT_INTERVAL = 5000    # save every 10 ps
ENERGY_INTERVAL = 5000    # log energy every 10 ps

cases = ["V_K392", "V_N392", "D_K392", "D_N392"]
all_results = {}

for case_name in cases:
    cif_path = STRUCT_DIR / f"{case_name}.cif"
    if not cif_path.exists():
        print(f"\nSKIPPING {case_name}: no CIF file")
        continue

    print(f"\n{'='*60}")
    print(f"  MD: {case_name}")
    print(f"{'='*60}")
    t_start = time.time()

    # Fix structure
    print("  Fixing structure with PDBFixer...")
    fixer = PDBFixer(str(cif_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    # Create system
    print("  Building system...")
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # Solvate in a box with 10A padding
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers,
                        ionicStrength=0.15 * unit.molar)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )
    system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

    integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Use GPU
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"CudaPrecision": "mixed"}

    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print("  Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    # Equilibration
    print(f"  Equilibrating ({EQUIL_STEPS * 2 / 1000:.0f} ps)...")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    simulation.step(EQUIL_STEPS)

    # Production
    out_dcd = RESULTS_DIR / f"{case_name}_traj.dcd"
    out_log = RESULTS_DIR / f"{case_name}_energy.csv"
    out_pdb = RESULTS_DIR / f"{case_name}_start.pdb"

    # Save starting structure
    state = simulation.context.getState(getPositions=True)
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    simulation.reporters.append(app.DCDReporter(str(out_dcd), REPORT_INTERVAL))
    simulation.reporters.append(app.StateDataReporter(
        str(out_log), ENERGY_INTERVAL,
        step=True, time=True, potentialEnergy=True, temperature=True,
        speed=True, separator=","
    ))

    print(f"  Production run ({PROD_STEPS * 2 / 1e6:.0f} ns)...")
    prod_start = time.time()
    simulation.step(PROD_STEPS)
    prod_time = time.time() - prod_start

    total_time = time.time() - t_start
    print(f"  Done in {total_time:.0f}s ({prod_time:.0f}s production)")

    # Basic analysis
    result = {
        "case": case_name,
        "peptide": "VKLLLLSI" if "V_" in case_name else "DKLLLLSI",
        "variant": "K392" if "K392" in case_name else "N392",
        "prod_ns": PROD_STEPS * 2e-6,
        "wall_time_s": total_time,
        "ns_per_day": (PROD_STEPS * 2e-6) / (prod_time / 86400),
    }

    # Trajectory analysis with mdtraj
    if HAS_MDTRAJ and out_dcd.exists():
        print("  Analyzing trajectory...")
        traj = mdtraj.load(str(out_dcd), top=str(out_pdb))

        # Find peptide chain (chain B, last chain)
        topology = traj.topology
        chains = list(topology.chains)
        protein_chains = [c for c in chains if any(r.name in mdtraj.core.residue_names._PROTEIN_RESIDUES for r in c.residues)]

        if len(protein_chains) >= 2:
            receptor_atoms = topology.select(f"chainid 0 and name CA")
            peptide_atoms = topology.select(f"chainid {protein_chains[-1].index} and name CA")

            if len(peptide_atoms) > 0 and len(receptor_atoms) > 0:
                # Peptide RMSD (relative to receptor)
                traj_aligned = traj.superpose(traj, atom_indices=receptor_atoms)
                pep_rmsd = mdtraj.rmsd(traj_aligned, traj_aligned, atom_indices=peptide_atoms) * 10  # nm -> A
                result["peptide_rmsd_mean_A"] = float(np.mean(pep_rmsd))
                result["peptide_rmsd_std_A"] = float(np.std(pep_rmsd))
                result["peptide_rmsd_max_A"] = float(np.max(pep_rmsd))

                # Interface contacts (4.5A cutoff)
                contacts_per_frame = []
                for frame in range(traj.n_frames):
                    pairs = []
                    for pa in peptide_atoms:
                        for ra in receptor_atoms[:50]:  # sample for speed
                            pairs.append((pa, ra))
                    if pairs:
                        distances = mdtraj.compute_distances(traj[frame], pairs)[0]
                        n_contacts = int(np.sum(distances < 0.45))  # 4.5A in nm
                        contacts_per_frame.append(n_contacts)
                if contacts_per_frame:
                    result["contacts_mean"] = float(np.mean(contacts_per_frame))
                    result["contacts_std"] = float(np.std(contacts_per_frame))
                    result["contacts_min"] = int(np.min(contacts_per_frame))
                    # "Dissociation" = frame where contacts drop below 2
                    dissociated = [i for i, c in enumerate(contacts_per_frame) if c < 2]
                    if dissociated:
                        result["first_dissociation_ns"] = dissociated[0] * REPORT_INTERVAL * 2e-6
                    else:
                        result["first_dissociation_ns"] = None  # never dissociated

                print(f"  RMSD: {result.get('peptide_rmsd_mean_A', '?'):.1f} +/- {result.get('peptide_rmsd_std_A', '?'):.1f} A")
                print(f"  Contacts: {result.get('contacts_mean', '?'):.1f} +/- {result.get('contacts_std', '?'):.1f}")
                if result.get("first_dissociation_ns") is None:
                    print(f"  Peptide stayed bound for full {result['prod_ns']} ns")
                else:
                    print(f"  First dissociation at {result['first_dissociation_ns']:.1f} ns")

    all_results[case_name] = result

# Save results
with open(RESULTS_DIR / "md_results.json", "w") as f:
    json.dump(all_results, f, indent=2)

# Summary table
print(f"\n{'='*70}")
print(f"  MD SUMMARY")
print(f"{'='*70}")
print(f"{'Case':>12s}  {'Peptide':>10s}  {'Variant':>7s}  {'RMSD(A)':>10s}  {'Contacts':>10s}  {'Dissoc(ns)':>12s}")
print("-" * 70)
for name, r in all_results.items():
    rmsd = f"{r.get('peptide_rmsd_mean_A', 0):.1f}+/-{r.get('peptide_rmsd_std_A', 0):.1f}" if "peptide_rmsd_mean_A" in r else "N/A"
    contacts = f"{r.get('contacts_mean', 0):.1f}" if "contacts_mean" in r else "N/A"
    dissoc = f"{r['first_dissociation_ns']:.1f}" if r.get("first_dissociation_ns") is not None else "STABLE"
    print(f"{name:>12s}  {r['peptide']:>10s}  {r['variant']:>7s}  {rmsd:>10s}  {contacts:>10s}  {dissoc:>12s}")

print(f"\nResults saved to {RESULTS_DIR / 'md_results.json'}")
MDPY

echo ""
echo "============================================"
echo "  MD Complete — $(date)"
echo "============================================"
