#!/bin/bash
# V4.2 Multi-Handle Campaign: 20 peptides x 4 targets x 5 diffusion samples
# Expected: ~75 min on H200, ~$3-4
set -e

echo "============================================"
echo "  V4.2 Multi-Handle Campaign — $(date)"
echo "  20 peptides x 4 targets x 5 seeds"
echo "============================================"

# ── Install deps ──────────────────────────────────────────────────────────────
echo "=== Installing dependencies ==="
pip3 install -q boltz 2>&1 | tail -3
pip3 install -q cuequivariance-torch cuequivariance-ops-torch-cu12 2>&1 | tail -3
echo "Dependencies ready"

# ── Generate experiment ───────────────────────────────────────────────────────
echo ""
echo "=== Generating V4.2 Campaign ==="
cd /workspace/validation
python3 v42_campaign.py

# ── Run Boltz-2 ───────────────────────────────────────────────────────────────
echo ""
YAML_DIR="v42_boltz2_inputs"
OUT_DIR="v42_boltz2_out"
mkdir -p "$OUT_DIR"

TOTAL=$(ls "$YAML_DIR"/*.yaml | wc -l)
DONE=0
FAILED=0
SKIPPED=0

echo "=== Running Boltz-2 ($TOTAL complexes x 5 diffusion samples) ==="

for yaml in "$YAML_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    DONE=$((DONE + 1))

    if [ -d "$OUT_DIR/$name" ] && find "$OUT_DIR/$name" -name "confidence_*.json" -print -quit 2>/dev/null | grep -q .; then
        echo "[$DONE/$TOTAL] $name — skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "[$DONE/$TOTAL] $name ..."
    if boltz predict "$yaml" \
        --diffusion_samples 5 \
        --seed 42 \
        --out_dir "$OUT_DIR/$name" 2>&1 | tail -3; then
        echo "  OK"
    else
        echo "  FAILED"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Boltz-2 Complete ==="
echo "  Total: $TOTAL | Skipped: $SKIPPED | Failed: $FAILED"

# ── Analyze ───────────────────────────────────────────────────────────────────
echo ""
echo "=== V4.2 Analysis ==="
python3 v42_campaign.py --analyze --results-dir "$OUT_DIR"

echo ""
echo "============================================"
echo "  V4.2 Campaign Complete — $(date)"
echo "============================================"
