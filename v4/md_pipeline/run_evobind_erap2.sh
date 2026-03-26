#!/bin/bash
# EvoBind: Design + validate 6-mer peptides against ERAP2 K392 active site.
#
# Three experiments:
#   1. PREDICT_ONLY: Score VAGSAF in ERAP2 — does AF2 agree with Boltz-2?
#   2. DESIGN: Evolve new 6-mers targeting K392 pocket (1000 iterations)
#   3. DESIGN_BROAD: Evolve 6-mers targeting full active site (1000 iterations)
#
# Usage: bash /workspace/md_pipeline/run_evobind_erap2.sh
#
# Prerequisites: bash /workspace/md_pipeline/gpu_setup_validation.sh
# Expected runtime: ~2-3 hrs on RTX 4090
set -e

EVOBIND=/workspace/evobind
PIPELINE=/workspace/md_pipeline
OUTDIR=/workspace/results/evobind
AF2_PARAMS=$EVOBIND/src/AF2/params
ERAP2_FASTA=$PIPELINE/erap2_k392.fasta
ERAP2_MSA=$EVOBIND/data/erap2_msa/erap2_k392.a3m

# ERAP2 K392 active site residues (0-indexed)
# Key selectivity residues from structural equivalence table:
#   K392=391, Y398=397, A403=402, A406=405, Q412=411, D414=413
# Extended pocket (±3 residues around each key position):
K392_POCKET="388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416"

# Tight core (just the 6 key selectivity residues)
K392_CORE="391,397,402,405,411,413"

echo "============================================"
echo "  EvoBind ERAP2 Peptide Design"
echo "  $(date)"
echo "============================================"

# Verify setup
if [ ! -f "$AF2_PARAMS/params_model_1.npz" ]; then
    echo "ERROR: AF2 params not found. Run gpu_setup_validation.sh first."
    exit 1
fi
if [ ! -f "$ERAP2_MSA" ]; then
    echo "ERROR: ERAP2 MSA not found. Run gpu_setup_validation.sh first."
    exit 1
fi

mkdir -p $OUTDIR

# ============================================================
# Experiment 1: Predict VAGSAF binding (no optimization)
# Does AF2 predict VAGSAF binds the K392 pocket?
# ============================================================
echo ""
echo "=== Experiment 1: VAGSAF prediction (no optimization) ==="
EXP1=$OUTDIR/exp1_vagsaf_predict
mkdir -p $EXP1

if [ ! -f "$EXP1/metrics.csv" ]; then
    cd $EVOBIND
    python3 src/mc_design.py \
        --receptor_fasta_path=$ERAP2_FASTA \
        --peptide_length=6 \
        --peptide_sequence=VAGSAF \
        --receptor_if_residues=$K392_POCKET \
        --msas=$ERAP2_MSA \
        --output_dir=$EXP1/ \
        --model_names=model_1 \
        --data_dir=$AF2_PARAMS \
        --max_recycles=8 \
        --num_iterations=1 \
        --predict_only=True
    echo "  VAGSAF prediction complete"
else
    echo "  Already completed"
fi

# ============================================================
# Experiment 2: Design 6-mers targeting K392 pocket
# 1000 iterations of directed evolution
# ============================================================
echo ""
echo "=== Experiment 2: Design 6-mers (K392 pocket, 1000 iter) ==="
EXP2=$OUTDIR/exp2_design_pocket
mkdir -p $EXP2

if [ ! -f "$EXP2/metrics.csv" ]; then
    cd $EVOBIND
    python3 src/mc_design.py \
        --receptor_fasta_path=$ERAP2_FASTA \
        --peptide_length=6 \
        --receptor_if_residues=$K392_POCKET \
        --msas=$ERAP2_MSA \
        --output_dir=$EXP2/ \
        --model_names=model_1 \
        --data_dir=$AF2_PARAMS \
        --max_recycles=8 \
        --num_iterations=1000
    echo "  Pocket design complete"
else
    echo "  Already completed"
fi

# ============================================================
# Experiment 3: Design 6-mers targeting core selectivity residues
# Tighter focus on the 6 ERAP2-unique positions
# ============================================================
echo ""
echo "=== Experiment 3: Design 6-mers (core selectivity, 1000 iter) ==="
EXP3=$OUTDIR/exp3_design_core
mkdir -p $EXP3

if [ ! -f "$EXP3/metrics.csv" ]; then
    cd $EVOBIND
    python3 src/mc_design.py \
        --receptor_fasta_path=$ERAP2_FASTA \
        --peptide_length=6 \
        --receptor_if_residues=$K392_CORE \
        --msas=$ERAP2_MSA \
        --output_dir=$EXP3/ \
        --model_names=model_1 \
        --data_dir=$AF2_PARAMS \
        --max_recycles=8 \
        --num_iterations=1000
    echo "  Core design complete"
else
    echo "  Already completed"
fi

# ============================================================
# Analyze results
# ============================================================
echo ""
echo "=== Results Summary ==="
python3 << 'ANALYZE'
import os, csv, json, sys

OUTDIR = "/workspace/results/evobind"
experiments = {
    "exp1_vagsaf_predict": "VAGSAF Prediction",
    "exp2_design_pocket": "Pocket Design (1000 iter)",
    "exp3_design_core": "Core Design (1000 iter)",
}

results = {}
for exp_dir, name in experiments.items():
    metrics_path = os.path.join(OUTDIR, exp_dir, "metrics.csv")
    if not os.path.exists(metrics_path):
        print(f"  {name}: NO RESULTS")
        continue

    with open(metrics_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print(f"  {name}: EMPTY")
        continue

    # Get best iteration (lowest loss)
    best = min(rows, key=lambda r: float(r.get("loss", 999)))
    last = rows[-1]

    results[exp_dir] = {
        "name": name,
        "n_iterations": len(rows),
        "best_iteration": int(best.get("iteration", 0)),
        "best_sequence": best.get("sequence", "?"),
        "best_loss": float(best.get("loss", 0)),
        "best_if_dist": float(best.get("if_dist_peptide", 0)),
        "best_plddt": float(best.get("plddt", 0)),
        "final_sequence": last.get("sequence", "?"),
        "final_loss": float(last.get("loss", 0)),
    }

    print(f"\n  {name}:")
    print(f"    Best sequence: {results[exp_dir]['best_sequence']}")
    print(f"    Best loss: {results[exp_dir]['best_loss']:.4f}")
    print(f"    Interface dist: {results[exp_dir]['best_if_dist']:.2f} A")
    print(f"    pLDDT: {results[exp_dir]['best_plddt']:.2f}")
    print(f"    Iteration: {results[exp_dir]['best_iteration']}/{results[exp_dir]['n_iterations']}")

# Save summary
summary_path = os.path.join(OUTDIR, "evobind_summary.json")
with open(summary_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  Summary saved: {summary_path}")

# Collect unique designed sequences for Boltz-2 cross-screening
designed_seqs = set()
for exp_dir, data in results.items():
    if "predict" not in exp_dir:
        designed_seqs.add(data["best_sequence"])
        designed_seqs.add(data["final_sequence"])

if designed_seqs:
    print(f"\n  Designed sequences for Boltz-2 cross-screen:")
    for seq in sorted(designed_seqs):
        print(f"    {seq}")
    # Save for Boltz-2 pipeline
    with open(os.path.join(OUTDIR, "designed_sequences.txt"), "w") as f:
        for seq in sorted(designed_seqs):
            f.write(seq + "\n")
    print(f"  Saved to: {OUTDIR}/designed_sequences.txt")
    print(f"\n  NEXT STEP: Run these through Boltz-2 against ERAP2, ERAP1, and IRAP")
    print(f"  to check if EvoBind designs are also triple-selective.")
ANALYZE

echo ""
echo "============================================"
echo "  EvoBind complete: $(date)"
echo "  Results: $OUTDIR/"
echo "============================================"
