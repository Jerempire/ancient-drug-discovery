#!/bin/bash
# n392_vastai_run.sh — Run N392 selectivity mutation screen on Vast.ai
#
# Upload this + boltz_yamls/ to /workspace, then run:
#   bash /workspace/n392_vastai_run.sh
#
# Expected: ~34 min on RTX 4090, 68 predictions (single sample each)

set -euo pipefail

echo "============================================"
echo "  N392 Selectivity Mutation Screen"
echo "  $(date)"
echo "============================================"

# Install Boltz-2
echo "[1/3] Installing Boltz-2..."
pip install -q boltz biopython 2>&1 | tail -3

# Verify GPU
python3 -c "
import torch
print('PyTorch:', torch.__version__)
print('CUDA:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU:', torch.cuda.get_device_name(0))
    print('VRAM: %.1f GB' % (torch.cuda.get_device_properties(0).total_memory / 1e9))
"

# Verify boltz
which boltz || (echo "ERROR: boltz not found after pip install"; exit 1)

# Create results directory
mkdir -p /workspace/results

# Count YAMLs
N_YAMLS=$(ls /workspace/yamls/*.yaml 2>/dev/null | wc -l)
echo ""
echo "[2/3] Running Boltz-2 on $N_YAMLS YAMLs..."
echo "Estimated time: ~$((N_YAMLS * 30 / 60)) minutes"
echo ""

# Run Boltz-2 on each YAML
DONE=0
FAILED=0
START_TIME=$(date +%s)

for yaml_file in /workspace/yamls/*.yaml; do
    name=$(basename "$yaml_file" .yaml)
    out_dir="/workspace/results/$name"

    if [ -d "$out_dir" ] && ls "$out_dir"/*/confidence_*.json >/dev/null 2>&1; then
        echo "  SKIP $name (already done)"
        DONE=$((DONE + 1))
        continue
    fi

    echo -n "  $name... "
    if boltz predict "$yaml_file" --out_dir "$out_dir" --diffusion_samples 1 > /dev/null 2>&1; then
        # Extract ipTM from confidence file
        IPTM=$(python3 -c "
import json, glob
files = glob.glob('$out_dir/*/confidence_*.json')
if files:
    with open(files[0]) as f:
        c = json.load(f)
    print('%.3f' % c.get('iptm', c.get('i_ptm', 0)))
else:
    print('N/A')
" 2>/dev/null || echo "N/A")
        echo "ipTM=$IPTM"
        DONE=$((DONE + 1))
    else
        echo "FAILED"
        FAILED=$((FAILED + 1))
    fi
done

END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 60 ))

echo ""
echo "[3/3] Summary"
echo "============================================"
echo "  Completed: $DONE / $N_YAMLS"
echo "  Failed: $FAILED"
echo "  Time: ${ELAPSED}m"
echo "  Results: /workspace/results/"
echo "============================================"
echo ""
echo "Download results with:"
echo "  scp -P <PORT> -r root@<HOST>:/workspace/results/ ./boltz_results/"
