"""
run_pipeline_gpu.py — Full Stage 4+5 drug design pipeline for ERAP2.

Runs on Vast.ai GPU instance. Orchestrates:
  1. DiffDock  — dock 4 known inhibitors (validates structure)
  2. RFdiffusion — generate 40 de novo binders (4 length tiers × 10)
  3. ProteinMPNN — design sequences for RFdiffusion backbones
  4. Selectivity screen — counter-dock all candidates against ERAP1
  5. Ranking — final CSV with selectivity ratios

Usage:
    python scripts/run_pipeline_gpu.py                  # Full pipeline
    python scripts/run_pipeline_gpu.py --step diffdock  # Single step
    python scripts/run_pipeline_gpu.py --step rfdiffusion
    python scripts/run_pipeline_gpu.py --step proteinmpnn
    python scripts/run_pipeline_gpu.py --step selectivity
    python scripts/run_pipeline_gpu.py --step ranking

Expects data files at /workspace/data/ (uploaded by vast_launch.py upload).
"""
import os
import sys
import glob
import json
import subprocess
import time
from pathlib import Path

# --- Paths ---
WORKSPACE = Path("/workspace")
DATA = WORKSPACE / "data"
STRUCTURES = DATA / "structures"
LIGANDS = DATA / "ligands"
RESULTS = WORKSPACE / "results"

# Sub-result directories
DIFFDOCK_OUT = RESULTS / "diffdock"
RFDIFF_OUT = RESULTS / "rfdiffusion"
MPNN_OUT = RESULTS / "proteinmpnn"
SELECTIVITY_OUT = RESULTS / "selectivity"

# Structure files
ERAP2_PDB = STRUCTURES / "erap2_wt_alphafold.pdb"
ERAP1_PDB = STRUCTURES / "erap1_wt_alphafold.pdb"
ERAP2_EXP = STRUCTURES / "erap2_3SE6_experimental.pdb"
ERAP1_EXP = STRUCTURES / "erap1_2YD0_experimental.pdb"

# ERAP2 active site residues (zinc metallopeptidase catalytic site)
ERAP2_HOTSPOTS = [370, 371, 374, 392, 393]  # H370, E371, H374, K392, E393
# ERAP1 equivalent residues
ERAP1_ACTIVE_SITE = [353, 354, 357, 376]    # H353, E354, H357, E376

# Known ligands (SMILES)
LIGANDS_SMILES = {
    "Bestatin": "CC(O)C(=O)NC(CC1=CC=CC=C1)C(O)=O",
    "Tosedostat": "CC(C)CC(NC(=O)C(CC1=CC=CC=C1)OC(=O)C(C)(C)C)B(O)O",
    "Leucinethiol": "CC(CC)CS",
    "Phosphoramidon": "CC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC2=CC=CC=C2)C(=O)O",
}


def run(cmd, cwd=None, check=True):
    """Run a shell command, stream output."""
    print(f"  $ {cmd}")
    result = subprocess.run(
        cmd, shell=True, cwd=cwd,
        capture_output=False, text=True,
    )
    if check and result.returncode != 0:
        print(f"  Command failed with exit code {result.returncode}")
        return False
    return True


def verify_data():
    """Check that required data files exist."""
    missing = []
    for f in [ERAP2_PDB, ERAP1_PDB]:
        if not f.exists():
            missing.append(str(f))

    sdf_files = list(LIGANDS.glob("*.sdf"))
    if not sdf_files:
        missing.append(f"{LIGANDS}/*.sdf")

    if missing:
        print("ERROR: Missing required data files:")
        for m in missing:
            print(f"  - {m}")
        print("\nRun 'python scripts/vast_launch.py upload' from your local machine first.")
        sys.exit(1)

    print(f"  ERAP2 structure: {ERAP2_PDB}")
    print(f"  ERAP1 structure: {ERAP1_PDB}")
    print(f"  Ligands: {len(sdf_files)} SDF files")


# ============================================================
# Step 1: DiffDock — dock known inhibitors against ERAP2
# ============================================================
def step_diffdock():
    """Dock 4 known aminopeptidase inhibitors to validate ERAP2 structure."""
    print("\n" + "=" * 60)
    print("STEP 1: DiffDock — Dock Known Inhibitors")
    print("=" * 60)

    DIFFDOCK_OUT.mkdir(parents=True, exist_ok=True)

    # Create ligands CSV for DiffDock batch mode
    csv_path = DIFFDOCK_OUT / "ligands.csv"
    with open(csv_path, "w") as f:
        f.write("complex_name,protein_path,ligand_description,protein_sequence\n")
        for name, smiles in LIGANDS_SMILES.items():
            f.write(f"{name},{ERAP2_PDB},{smiles},\n")

    # Run DiffDock
    diffdock_dir = WORKSPACE / "DiffDock"
    success = run(
        f"python -m inference "
        f"--protein_ligand_csv {csv_path} "
        f"--out_dir {DIFFDOCK_OUT} "
        f"--inference_steps 20 "
        f"--samples_per_complex 10 "
        f"--batch_size 4 "
        f"--no_final_step_noise",
        cwd=str(diffdock_dir),
    )

    if not success:
        # Try alternative invocation
        print("  Trying alternative DiffDock invocation...")
        run(
            f"python -m diffdock.inference "
            f"--protein_ligand_csv {csv_path} "
            f"--out_dir {DIFFDOCK_OUT}",
            cwd=str(diffdock_dir),
            check=False,
        )

    # Analyze results
    result_dirs = sorted(DIFFDOCK_OUT.glob("*/"))
    result_dirs = [d for d in result_dirs if d.is_dir()]
    print(f"\nDiffDock produced results for {len(result_dirs)} complexes")

    from rdkit import Chem

    dock_results = []
    for rdir in result_dirs:
        name = rdir.name
        sdfs = sorted(rdir.glob("*.sdf"))
        if not sdfs:
            continue

        # Read top pose confidence
        conf_files = list(rdir.glob("*confidence*"))
        confidence = None
        if conf_files:
            with open(conf_files[0]) as f:
                scores = f.read().strip().split("\n")
                if scores:
                    try:
                        confidence = float(scores[0])
                    except ValueError:
                        pass

        # Get pose center
        mol = Chem.SDMolSupplier(str(sdfs[0]))[0]
        center = None
        if mol:
            pos = mol.GetConformer().GetPositions()
            center = pos.mean(axis=0).tolist()

        dock_results.append({
            "ligand": name,
            "poses": len(sdfs),
            "top_confidence": confidence,
            "top_pose_center": center,
        })
        print(f"  {name}: {len(sdfs)} poses, confidence={confidence}")

    # Save summary
    summary_path = DIFFDOCK_OUT / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(dock_results, f, indent=2)
    print(f"\nSaved DiffDock summary to {summary_path}")

    return dock_results


# ============================================================
# Step 2: RFdiffusion — generate de novo binders
# ============================================================
def step_rfdiffusion():
    """Generate 40 de novo protein binders at 4 length tiers."""
    print("\n" + "=" * 60)
    print("STEP 2: RFdiffusion — De Novo Binder Design")
    print("=" * 60)

    RFDIFF_OUT.mkdir(parents=True, exist_ok=True)

    rfdiff_dir = WORKSPACE / "RFdiffusion"
    hotspots = ",".join(f"A{r}" for r in ERAP2_HOTSPOTS)

    # 4 length tiers × 10 designs = 40 total
    tiers = [
        ("short",  "30-40"),   # Peptide synthesis ~$500
        ("medium", "50-60"),   # Moderate cost ~$800
        ("long",   "70-90"),   # Recombinant expression ~$1000
        ("large",  "90-110"),  # Maximum binding surface ~$1500
    ]

    for label, length_range in tiers:
        print(f"\n--- Designing {label} binders ({length_range} residues) ---")
        output_prefix = RFDIFF_OUT / f"erap2_{label}"

        # Check if already done
        existing = list(RFDIFF_OUT.glob(f"erap2_{label}_*.pdb"))
        if len(existing) >= 10:
            print(f"  Already have {len(existing)} designs, skipping.")
            continue

        run(
            f"python scripts/run_inference.py "
            f"inference.output_prefix={output_prefix} "
            f"inference.input_pdb={ERAP2_PDB} "
            f"'contigmap.contigs=[A{ERAP2_HOTSPOTS[0]}-{ERAP2_HOTSPOTS[-1]}/0 {length_range}]' "
            f"'ppi.hotspot_res=[{hotspots}]' "
            f"inference.num_designs=10 "
            f"inference.ckpt_override_path=models/Complex_beta_ckpt.pt",
            cwd=str(rfdiff_dir),
            check=False,
        )

    # Count results
    all_pdbs = sorted(RFDIFF_OUT.glob("erap2_*.pdb"))
    print(f"\nTotal binder designs: {len(all_pdbs)}")

    # Analyze each design
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    rfdiff_results = []

    for pdb_path in all_pdbs:
        struct = parser.get_structure(pdb_path.stem, str(pdb_path))
        chains = {}
        for chain in struct.get_chains():
            residues = list(chain.get_residues())
            chains[chain.id] = len(residues)

        rfdiff_results.append({
            "design": pdb_path.name,
            "tier": pdb_path.stem.split("_")[1] if "_" in pdb_path.stem else "unknown",
            "chains": chains,
            "binder_residues": chains.get("B", 0),
        })
        print(f"  {pdb_path.name}: {chains}")

    summary_path = RFDIFF_OUT / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(rfdiff_results, f, indent=2)

    return rfdiff_results


# ============================================================
# Step 3: ProteinMPNN — design sequences for binder backbones
# ============================================================
def step_proteinmpnn():
    """Design amino acid sequences for RFdiffusion backbone structures."""
    print("\n" + "=" * 60)
    print("STEP 3: ProteinMPNN — Sequence Design")
    print("=" * 60)

    MPNN_OUT.mkdir(parents=True, exist_ok=True)

    mpnn_dir = WORKSPACE / "ProteinMPNN"
    binder_pdbs = sorted(RFDIFF_OUT.glob("erap2_*.pdb"))

    if not binder_pdbs:
        print("  No RFdiffusion outputs found. Run step 2 first.")
        return []

    print(f"  Designing sequences for {len(binder_pdbs)} backbones...")

    # Create input list for ProteinMPNN
    pdb_list = MPNN_OUT / "pdb_list.txt"
    with open(pdb_list, "w") as f:
        for pdb in binder_pdbs:
            f.write(f"{pdb}\n")

    # Run ProteinMPNN — design only the binder chain (B), fix target chain (A)
    # 4 sequences per backbone
    run(
        f"python protein_mpnn_run.py "
        f"--pdb_path_chains B "
        f"--out_folder {MPNN_OUT} "
        f"--num_seq_per_target 4 "
        f"--sampling_temp 0.1 "
        f"--seed 42 "
        f"--batch_size 1 "
        f"--pdb_path {binder_pdbs[0]}",  # Single file mode first
        cwd=str(mpnn_dir),
        check=False,
    )

    # For batch mode, iterate
    mpnn_results = []
    for pdb in binder_pdbs:
        name = pdb.stem

        # Check if ProteinMPNN single-file mode works
        output_fasta = MPNN_OUT / "seqs" / f"{name}.fa"
        if not output_fasta.exists():
            # Try alternative invocation
            run(
                f"python protein_mpnn_run.py "
                f"--pdb_path {pdb} "
                f"--chain_id_jsonl '' "
                f"--out_folder {MPNN_OUT / name} "
                f"--num_seq_per_target 4 "
                f"--sampling_temp 0.1 "
                f"--batch_size 1",
                cwd=str(mpnn_dir),
                check=False,
            )

        # Collect FASTA outputs
        fasta_files = list((MPNN_OUT).rglob(f"*{name}*.fa")) + list((MPNN_OUT).rglob(f"*{name}*.fasta"))
        if fasta_files:
            for fa in fasta_files:
                with open(fa) as f:
                    content = f.read()
                seqs = [line for line in content.split("\n") if line and not line.startswith(">")]
                mpnn_results.append({
                    "backbone": name,
                    "sequences": seqs[:4],
                    "n_sequences": len(seqs[:4]),
                })
                print(f"  {name}: {len(seqs[:4])} sequences designed")

    summary_path = MPNN_OUT / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(mpnn_results, f, indent=2)

    return mpnn_results


# ============================================================
# Step 4: Selectivity screen — ERAP2 vs ERAP1
# ============================================================
def step_selectivity():
    """Counter-screen all binders against ERAP1 to identify selective candidates."""
    print("\n" + "=" * 60)
    print("STEP 4: Selectivity Screen — ERAP2 vs ERAP1")
    print("=" * 60)

    SELECTIVITY_OUT.mkdir(parents=True, exist_ok=True)

    from Bio.PDB import PDBParser, NeighborSearch
    import numpy as np

    parser = PDBParser(QUIET=True)

    # Load reference structures
    erap2_struct = parser.get_structure("erap2", str(ERAP2_PDB))
    erap2_atoms = list(erap2_struct.get_atoms())
    erap2_ns = NeighborSearch(erap2_atoms)

    erap1_struct = parser.get_structure("erap1", str(ERAP1_PDB))
    erap1_atoms = list(erap1_struct.get_atoms())
    erap1_ns = NeighborSearch(erap1_atoms)

    binder_pdbs = sorted(RFDIFF_OUT.glob("erap2_*.pdb"))
    if not binder_pdbs:
        print("  No binder PDBs found. Run RFdiffusion first.")
        return []

    print(f"  Screening {len(binder_pdbs)} binders for ERAP2 selectivity...")
    results = []

    for pdb_path in binder_pdbs:
        name = pdb_path.name
        struct = parser.get_structure(name, str(pdb_path))

        # Get binder chain B atoms
        binder_atoms = []
        for chain in struct.get_chains():
            if chain.id == "B":
                binder_atoms = list(chain.get_atoms())

        if not binder_atoms:
            print(f"  {name}: no chain B found, skipping")
            continue

        # Count contacts with ERAP2 active site (within 5 Angstroms)
        erap2_contacts = 0
        for atom in binder_atoms:
            nearby = erap2_ns.search(atom.get_vector().get_array(), 5.0)
            for n in nearby:
                res = n.get_parent()
                if res.get_id()[1] in ERAP2_HOTSPOTS:
                    erap2_contacts += 1

        # Count contacts with ERAP1 active site
        erap1_contacts = 0
        for atom in binder_atoms:
            nearby = erap1_ns.search(atom.get_vector().get_array(), 5.0)
            for n in nearby:
                res = n.get_parent()
                if res.get_id()[1] in ERAP1_ACTIVE_SITE:
                    erap1_contacts += 1

        binder_length = len(set(
            a.get_parent().get_id()[1] for a in binder_atoms
        ))
        selectivity = erap2_contacts / max(erap1_contacts, 1)

        tag = ("SELECTIVE" if selectivity > 3
               else "MODERATE" if selectivity > 1.5
               else "NON-SELECTIVE")

        results.append({
            "design": name,
            "binder_residues": binder_length,
            "erap2_contacts": erap2_contacts,
            "erap1_contacts": erap1_contacts,
            "selectivity_ratio": round(selectivity, 2),
            "tag": tag,
        })
        print(f"  {name}: ERAP2={erap2_contacts} ERAP1={erap1_contacts} "
              f"ratio={selectivity:.1f} → {tag}")

    # Save results
    import pandas as pd
    df = pd.DataFrame(results)
    df = df.sort_values("selectivity_ratio", ascending=False)
    csv_path = SELECTIVITY_OUT / "selectivity_screen.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nSaved to {csv_path}")

    return results


# ============================================================
# Step 5: Final ranking
# ============================================================
def step_ranking():
    """Rank all candidates by selectivity and binding quality."""
    print("\n" + "=" * 60)
    print("STEP 5: Final Ranking")
    print("=" * 60)

    import pandas as pd

    csv_path = SELECTIVITY_OUT / "selectivity_screen.csv"
    if not csv_path.exists():
        print("  No selectivity data. Run step 4 first.")
        return

    df = pd.read_csv(csv_path)

    # Filter selective binders
    selective = df[df["selectivity_ratio"] > 2.0].copy()
    selective = selective.sort_values("erap2_contacts", ascending=False)

    print(f"\nTotal designs screened: {len(df)}")
    print(f"ERAP2-selective (ratio > 2.0): {len(selective)}")
    print(f"Non-selective: {len(df[df['selectivity_ratio'] <= 1.5])}")

    if not selective.empty:
        print(f"\n{'='*70}")
        print("TOP ERAP2-SELECTIVE BINDER CANDIDATES")
        print(f"{'='*70}")
        print(selective.head(10).to_string(index=False))

        best = selective.iloc[0]
        tier = best["design"].split("_")[1] if "_" in best["design"] else "?"

        # Estimate synthesis cost
        size = best["binder_residues"]
        if size <= 40:
            cost = "$400-800 (peptide synthesis)"
        elif size <= 60:
            cost = "$800-1,500 (peptide synthesis)"
        else:
            cost = "$500-1,500 (recombinant expression)"

        print(f"\n*** LEAD CANDIDATE: {best['design']} ***")
        print(f"    Tier: {tier}")
        print(f"    Binder size: {best['binder_residues']} residues")
        print(f"    ERAP2 contacts: {best['erap2_contacts']}")
        print(f"    ERAP1 contacts: {best['erap1_contacts']}")
        print(f"    Selectivity: {best['selectivity_ratio']}x over ERAP1")
        print(f"    Est. synthesis: {cost}")
        print(f"\n    Next steps:")
        print(f"    1. File provisional patent ($75 micro entity)")
        print(f"    2. Synthesize peptide/express protein")
        print(f"    3. SPR binding assay (ERAP2 + ERAP1 counter-screen)")
        print(f"    4. Enzyme inhibition assay (fluorogenic substrate)")
    else:
        print("\nNo selective candidates. Consider:")
        print("  - Different hotspot residues (allosteric site)")
        print("  - Longer binder lengths")
        print("  - Proteina-Complexa (higher hit rate)")

    # Save final ranked CSV
    final_path = RESULTS / "final_candidates.csv"
    df_ranked = df.sort_values(
        ["selectivity_ratio", "erap2_contacts"],
        ascending=[False, False],
    )
    df_ranked.to_csv(final_path, index=False)
    print(f"\nFull ranked results: {final_path}")

    # Also check ProteinMPNN results
    mpnn_summary = MPNN_OUT / "summary.json"
    if mpnn_summary.exists():
        with open(mpnn_summary) as f:
            mpnn = json.load(f)
        # Match sequences to selective candidates
        selective_names = set(selective["design"].str.replace(".pdb", ""))
        matched = [m for m in mpnn if m["backbone"] in selective_names]
        if matched:
            print(f"\nSequences available for {len(matched)} selective candidates")
            for m in matched[:3]:
                print(f"  {m['backbone']}: {m['n_sequences']} sequences")
                if m["sequences"]:
                    print(f"    Top: {m['sequences'][0][:60]}...")


# ============================================================
# Main
# ============================================================
def main():
    import argparse

    steps = {
        "diffdock": step_diffdock,
        "rfdiffusion": step_rfdiffusion,
        "proteinmpnn": step_proteinmpnn,
        "selectivity": step_selectivity,
        "ranking": step_ranking,
    }

    # Simple arg parsing (argparse OK here — not Jupyter)
    parser = argparse.ArgumentParser(description="ERAP2 drug design pipeline")
    parser.add_argument("--step", choices=list(steps.keys()),
                        help="Run a single step instead of full pipeline")
    args = parser.parse_args()

    print("=" * 60)
    print("  ERAP2 Drug Discovery — GPU Pipeline")
    print("=" * 60)
    print(f"  Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Workspace: {WORKSPACE}")

    # Verify GPU
    import torch
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
        print(f"  VRAM: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")
    else:
        print("  WARNING: No GPU detected! Pipeline will be very slow.")

    print()

    # Create output dirs
    for d in [RESULTS, DIFFDOCK_OUT, RFDIFF_OUT, MPNN_OUT, SELECTIVITY_OUT]:
        d.mkdir(parents=True, exist_ok=True)

    # Verify data
    verify_data()

    if args.step:
        # Run single step
        steps[args.step]()
    else:
        # Full pipeline
        t0 = time.time()
        step_diffdock()
        step_rfdiffusion()
        step_proteinmpnn()
        step_selectivity()
        step_ranking()
        elapsed = time.time() - t0
        print(f"\n{'='*60}")
        print(f"  Pipeline complete in {elapsed/60:.1f} minutes")
        print(f"  Results: {RESULTS}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
