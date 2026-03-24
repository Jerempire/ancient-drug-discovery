#!/bin/bash
# V4 Validation Experiment: 18 peptides x 2 alleles x 5 diffusion samples
# Expected: ~3-4 hours on RTX 4090, ~$1.50
set -e

echo "============================================"
echo "  V4 Validation Experiment — $(date)"
echo "  18 peptides x 2 alleles x 5 seeds"
echo "============================================"

# ── Install deps ──────────────────────────────────────────────────────────────
echo "=== Installing dependencies ==="
pip install -q boltz biopython 2>&1 | tail -3

# Uninstall broken cuequivariance if present (torch 2.5.1 incompatible)
pip uninstall -y cuequivariance-torch cuequivariance-ops-torch-cu12 2>/dev/null || true

echo "Dependencies ready"

# ── Generate experiment ───────────────────────────────────────────────────────
echo ""
echo "=== Generating Experiment ==="
cd /workspace/validation
python validation_experiment.py

# ── Run Boltz-2 ───────────────────────────────────────────────────────────────
echo ""
echo "=== Running Boltz-2 (36 complexes x 5 diffusion samples) ==="

YAML_DIR="boltz2_inputs"
OUT_DIR="boltz2_out"
mkdir -p "$OUT_DIR"

TOTAL=$(ls "$YAML_DIR"/*.yaml | wc -l)
DONE=0
FAILED=0
SKIPPED=0

for yaml in "$YAML_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    DONE=$((DONE + 1))

    # Skip if already completed (has confidence files)
    if [ -d "$OUT_DIR/$name" ] && find "$OUT_DIR/$name" -name "confidence_*.json" -print -quit 2>/dev/null | grep -q .; then
        echo "[$DONE/$TOTAL] $name — already done, skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "[$DONE/$TOTAL] $name ..."
    if boltz predict "$yaml" \
        --diffusion_samples 5 \
        --seed 42 \
        --out_dir "$OUT_DIR/$name" \
        --no_kernels 2>&1 | tail -2; then
        echo "  OK"
    else
        echo "  FAILED"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Boltz-2 Complete ==="
echo "  Total: $TOTAL"
echo "  Skipped: $SKIPPED"
echo "  Failed: $FAILED"

# ── Analyze ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Validation Analysis ==="
python validation_experiment.py --analyze --results-dir "$OUT_DIR"

echo ""
echo "============================================"
echo "  Validation Complete — $(date)"
echo "  Results: boltz2_out/"
echo "  Analysis: boltz2_out/validation_analysis.json"
echo "============================================"
