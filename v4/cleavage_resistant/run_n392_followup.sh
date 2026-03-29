#!/bin/bash
# N392 follow-up: test cleavage-resistant peptides against N392 allele
set -e

INPUT_DIR="/workspace/v4_cleavage_resistant"
OUTPUT_DIR="/workspace/v4_cr_results"
LOG_FILE="/workspace/v4_cr_n392.log"

echo "=== N392 Follow-up Screen ===" | tee "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"

for yaml in "$INPUT_DIR"/*n392*.yaml; do
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

echo "=== N392 screen complete ===" | tee -a "$LOG_FILE"
echo "Finished: $(date)" | tee -a "$LOG_FILE"
