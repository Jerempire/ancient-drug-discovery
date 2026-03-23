#!/bin/bash
# SAR Grid: Structured mapping of scaffold x P1 x length x variant
# 40 peptides x 2 variants = 80 Boltz-2 runs x 3 diffusion samples
# Expected: ~3-4 hours on RTX 4090, ~$1.00
set -e

echo "============================================"
echo "  SAR Grid Campaign — $(date)"
echo "  40 peptides x 2 variants x 3 seeds"
echo "============================================"

# ── Install deps ──────────────────────────────────────────────────────────────
echo "=== Installing dependencies ==="
pip install -q boltz biopython 2>&1 | tail -3

# Uninstall broken cuequivariance if present (torch 2.5.1 incompatible)
pip uninstall -y cuequivariance-torch cuequivariance-ops-torch-cu12 2>/dev/null || true

echo "Dependencies ready"

# ── Generate grid ─────────────────────────────────────────────────────────────
echo ""
echo "=== Generating SAR Grid ==="
cd /workspace/pepmlm
python sar_grid.py

# ── Run Boltz-2 ───────────────────────────────────────────────────────────────
echo ""
echo "=== Running Boltz-2 (80 complexes x 3 diffusion samples) ==="
echo "This will take ~3-4 hours on RTX 4090"

YAML_DIR="sar_results/boltz2_inputs"
OUT_DIR="sar_results/boltz2_out"
mkdir -p "$OUT_DIR"

# Run each YAML individually for progress tracking and crash recovery
TOTAL=$(ls "$YAML_DIR"/*.yaml | wc -l)
DONE=0
FAILED=0

for yaml in "$YAML_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    DONE=$((DONE + 1))

    # Skip if already completed
    if [ -d "$OUT_DIR/$name" ] && find "$OUT_DIR/$name" -name "confidence_*.json" -print -quit 2>/dev/null | grep -q .; then
        echo "[$DONE/$TOTAL] $name — already done, skipping"
        continue
    fi

    echo "[$DONE/$TOTAL] $name ..."
    if boltz predict "$yaml" \
        --diffusion_samples 3 \
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
echo "  Failed: $FAILED"

# ── Analyze ───────────────────────────────────────────────────────────────────
echo ""
echo "=== SAR Analysis ==="
python sar_grid.py --analyze --results-dir "$OUT_DIR"

echo ""
echo "============================================"
echo "  SAR Grid Complete — $(date)"
echo "  Results: sar_results/boltz2_out/"
echo "  Analysis: sar_results/boltz2_out/sar_analysis.json"
echo "============================================"
