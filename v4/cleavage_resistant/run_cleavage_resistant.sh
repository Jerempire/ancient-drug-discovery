#!/bin/bash
# V4 Cleavage-Resistant Peptide Screen
# Runs 12 Boltz-2 predictions (4 modifications × 3 targets)
# Expected runtime: ~20-30 min on V100-32GB

set -e

INPUT_DIR="/workspace/v4_cleavage_resistant"
OUTPUT_DIR="/workspace/v4_cr_results"
LOG_FILE="/workspace/v4_cr_run.log"

echo "=== V4 Cleavage-Resistant Screen ===" | tee "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Run each YAML individually to isolate failures
# (if CCD codes fail for one modification, others still run)
for yaml in "$INPUT_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    out="$OUTPUT_DIR/$name"

    echo "--- Running: $name ---" | tee -a "$LOG_FILE"
    echo "Start: $(date)" | tee -a "$LOG_FILE"

    if boltz predict "$yaml" \
        --diffusion_samples 3 \
        --seed 42 \
        --out_dir "$out" \
        --num_workers 0 \
        2>&1 | tee -a "$LOG_FILE"; then
        echo "SUCCESS: $name" | tee -a "$LOG_FILE"
    else
        echo "FAILED: $name (exit code $?)" | tee -a "$LOG_FILE"
    fi

    echo "" | tee -a "$LOG_FILE"
done

echo "=== All predictions complete ===" | tee -a "$LOG_FILE"
echo "Finished: $(date)" | tee -a "$LOG_FILE"

# Quick summary: extract ipTM from results
echo "" | tee -a "$LOG_FILE"
echo "=== ipTM Summary ===" | tee -a "$LOG_FILE"
for result_dir in "$OUTPUT_DIR"/*/; do
    name=$(basename "$result_dir")
    # Find confidence JSON files
    json_files=$(find "$result_dir" -name "confidence_*_model_0.json" 2>/dev/null)
    if [ -n "$json_files" ]; then
        for jf in $json_files; do
            iptm=$(python3 -c "import json; d=json.load(open('$jf')); print(f'{d.get(\"iptm\", \"N/A\")}')" 2>/dev/null || echo "parse_error")
            echo "$name: ipTM=$iptm" | tee -a "$LOG_FILE"
        done
    else
        echo "$name: NO RESULTS" | tee -a "$LOG_FILE"
    fi
done
