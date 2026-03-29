#!/bin/bash
# Kill everything, clean restart
pkill -9 -f "boltz predict" 2>/dev/null
pkill -9 -f "run_cleavage" 2>/dev/null
pkill -9 -f "run_kn" 2>/dev/null
sleep 3

# Verify clean
remaining=$(ps aux | grep -E "boltz predict|run_cleavage" | grep -v grep | wc -l)
echo "Remaining processes: $remaining"

if [ "$remaining" -gt 0 ]; then
    echo "FORCE KILLING ALL PYTHON"
    killall -9 python3 2>/dev/null
    sleep 2
fi

# Clean old results
rm -rf /workspace/v4_cr_results /workspace/v4_cr_run.log /workspace/v4_cr_nohup.log

# GPU check
nvidia-smi --query-gpu=memory.used --format=csv,noheader

# Start fresh - single process
echo "Starting fresh run..."
nohup bash /workspace/run_cleavage_resistant.sh > /workspace/v4_cr_nohup.log 2>&1 &
sleep 5

# Verify
ps aux | grep "boltz predict" | grep -v grep | wc -l
echo "Launched."
