#!/bin/bash
# Master orchestrator: MD + MMGBSA validation of elite four peptides.
#
# Runs 12 complexes (4 peptides x 3 targets) through:
#   1. CIF -> PDB conversion
#   2. 20ns MD + MMGBSA binding free energy
#   3. Selectivity analysis
#
# Usage: bash /workspace/md_pipeline/run_elite_four.sh [--ns 20]
#
# Expected runtime: ~4-6 hrs on RTX 4090 (~$1-2)
set -e

NS=${1:-20}
ELITE=/workspace/results/elite_four
OUTDIR=/workspace/results/md_mmgbsa
PIPELINE=/workspace/md_pipeline

mkdir -p $OUTDIR/pdbs $OUTDIR/results

echo "============================================"
echo "  Elite Four MD + MMGBSA Validation"
echo "  Sim length: ${NS}ns per complex"
echo "  Complexes: 12 (4 peptides x 3 targets)"
echo "  Started: $(date)"
echo "============================================"

# Define the matrix
PEPTIDES="VAGSAF IAFSAF VAWSAF FASGAV"
TARGETS="erap2k392 erap1 irap"

# Step 1: Convert all CIFs to PDB
echo ""
echo "=== Step 1: CIF -> PDB conversion ==="
for pep in $PEPTIDES; do
    for tgt in $TARGETS; do
        label="${pep}_vs_${tgt}"
        cif_dir="$ELITE/${label}/boltz_results_${label}/predictions/${label}"
        cif="$cif_dir/${label}_model_0.cif"
        pdb="$OUTDIR/pdbs/${label}.pdb"

        if [ ! -f "$cif" ]; then
            echo "  WARNING: Missing $cif — skipping"
            continue
        fi

        if [ -f "$pdb" ]; then
            echo "  Already converted: $label"
            continue
        fi

        python3 $PIPELINE/cif_to_pdb.py "$cif" "$pdb"
    done
done

echo ""
echo "=== Step 2: MD + MMGBSA (${NS}ns each) ==="
COMPLETED=0
TOTAL=12
for pep in $PEPTIDES; do
    for tgt in $TARGETS; do
        label="${pep}_vs_${tgt}"
        pdb="$OUTDIR/pdbs/${label}.pdb"
        result="$OUTDIR/results/${label}.json"

        if [ ! -f "$pdb" ]; then
            echo "  SKIP: No PDB for $label"
            continue
        fi

        if [ -f "$result" ]; then
            echo "  Already completed: $label"
            COMPLETED=$((COMPLETED + 1))
            continue
        fi

        echo ""
        echo "--- [$((COMPLETED + 1))/$TOTAL] $label ---"
        python3 $PIPELINE/md_mmgbsa.py "$pdb" "$label" "$result" --ns $NS

        COMPLETED=$((COMPLETED + 1))
        echo "  Progress: $COMPLETED/$TOTAL complete"
    done
done

# Step 3: Analyze selectivity
echo ""
echo "=== Step 3: Selectivity Analysis ==="
python3 $PIPELINE/analyze_selectivity.py $OUTDIR/results $OUTDIR/selectivity_report.json

echo ""
echo "============================================"
echo "  COMPLETE: $(date)"
echo "  Results: $OUTDIR/"
echo "============================================"
