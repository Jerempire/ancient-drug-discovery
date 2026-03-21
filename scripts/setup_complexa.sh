#!/bin/bash
# setup_complexa.sh — Install Proteina-Complexa on Vast.ai
# Usage: ssh root@host "bash /workspace/scripts/setup_complexa.sh"
set -euo pipefail

echo "============================================"
echo "  Proteina-Complexa Setup"
echo "============================================"

cd /workspace

# 1. Clone repo
if [ ! -d "Proteina-Complexa" ]; then
    echo "[1/4] Cloning Proteina-Complexa..."
    git clone --depth 1 https://github.com/NVIDIA-Digital-Bio/Proteina-Complexa.git
else
    echo "[1/4] Already cloned."
fi

cd Proteina-Complexa

# 2. Install via UV (recommended by repo)
echo "[2/4] Building UV environment..."
if [ ! -d ".venv" ]; then
    pip install -q uv 2>/dev/null || true
    ./env/build_uv_env.sh 2>&1 | tail -20
else
    echo "  .venv already exists"
fi

# Activate
source .venv/bin/activate 2>/dev/null || true
export PATH="/workspace/Proteina-Complexa/.venv/bin:$PATH"

# 3. Download model weights
echo "[3/4] Downloading model weights..."
complexa init 2>&1 | tail -5
complexa download 2>&1 | tail -10

# 4. Setup target data
echo "[4/4] Setting up ERAP2 target..."
mkdir -p assets/target_data/erap2
cp /workspace/data/structures/erap2_wt_alphafold.pdb assets/target_data/erap2/

# Add ERAP2 target to targets dict
cat >> configs/targets/targets_dict.yaml << 'TARGEOF'

  # ERAP2 divergent substrate channel — selective binder design
  erap2_channel_short:
    source: erap2
    target_filename: erap2_wt_alphafold
    target_input: A350-500
    hotspot_residues: [A353, A355, A360, A363, A367, A401, A403, A406, A408, A412]
    binder_length: [35, 55]
    pdb_id: 3SE6

  erap2_channel_medium:
    source: erap2
    target_filename: erap2_wt_alphafold
    target_input: A350-500
    hotspot_residues: [A353, A355, A360, A363, A367, A401, A403, A406, A408, A412]
    binder_length: [55, 75]
    pdb_id: 3SE6

  erap2_channel_long:
    source: erap2
    target_filename: erap2_wt_alphafold
    target_input: A350-500
    hotspot_residues: [A353, A355, A360, A363, A367, A401, A403, A406, A408, A412]
    binder_length: [75, 100]
    pdb_id: 3SE6
TARGEOF

echo ""
echo "============================================"
echo "  Setup complete!"
echo "  Run designs with:"
echo "    cd /workspace/Proteina-Complexa"
echo "    source .venv/bin/activate"
echo "    complexa design configs/search_binder_local_pipeline.yaml \\"
echo "      ++run_name=erap2_short ++generation.task_name=erap2_channel_short"
echo "============================================"
