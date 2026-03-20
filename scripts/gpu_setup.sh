#!/bin/bash
# gpu_setup.sh — Setup for Vast.ai GPU instance (rosettacommons/rfdiffusion image)
#
# IMPORTANT: Launch Vast.ai instances with image "rosettacommons/rfdiffusion:latest"
# This image ships with RFdiffusion + all model weights pre-installed.
# This script only installs the additional tools (ProteinMPNN, BioPython, etc).
#
# Usage:
#   ssh -p <PORT> root@<HOST> "cd /workspace && bash scripts/gpu_setup.sh"
#
# Expected runtime: ~2-3 min

set -euo pipefail

echo "============================================"
echo "  ERAP2 Drug Discovery — GPU Setup"
echo "  Image: rosettacommons/rfdiffusion:latest"
echo "============================================"
echo "Started: $(date)"
echo ""

# --- Verify RFdiffusion is pre-installed ---
echo "[1/4] Verifying RFdiffusion..."
if [ -d "/app/RFdiffusion" ]; then
    echo "  RFdiffusion found at /app/RFdiffusion"
    ls -lh /app/RFdiffusion/models/*.pt 2>/dev/null | head -3
    echo "  ($(ls /app/RFdiffusion/models/*.pt 2>/dev/null | wc -l) model weights)"
else
    echo "  ERROR: /app/RFdiffusion not found!"
    echo "  Did you launch with image 'rosettacommons/rfdiffusion:latest'?"
    exit 1
fi

# --- Install additional Python packages ---
echo ""
echo "[2/4] Installing additional packages (BioPython, pandas)..."
python3 -m ensurepip 2>/dev/null || true
python3 -m pip install -q biopython pandas 2>&1 | tail -3
echo "  Done."

# --- Clone ProteinMPNN ---
echo ""
echo "[3/4] Installing ProteinMPNN..."
if [ ! -d "/workspace/ProteinMPNN" ]; then
    git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git /workspace/ProteinMPNN 2>&1 | tail -2
    echo "  ProteinMPNN installed."
else
    echo "  ProteinMPNN already present."
fi

# --- Verify ---
echo ""
echo "[4/4] Verification..."
python3 -c "
import torch
print('PyTorch: %s' % torch.__version__)
print('CUDA: %s' % torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU: %s' % torch.cuda.get_device_name(0))
    print('VRAM: %.1f GB' % (torch.cuda.get_device_properties(0).total_memory / 1e9))
"

cd /app/RFdiffusion
python3 -c "from rfdiffusion.inference.utils import parse_pdb; print('RFdiffusion: OK')" 2>&1 | tail -1
cd /workspace

python3 -c "from Bio.PDB import PDBParser; print('BioPython: OK')" 2>/dev/null || echo "BioPython: FAILED"
ls /workspace/ProteinMPNN/protein_mpnn_run.py > /dev/null 2>&1 && echo "ProteinMPNN: OK" || echo "ProteinMPNN: FAILED"

echo ""
echo "Data files:"
ls /workspace/data/structures/*.pdb 2>/dev/null || echo "  No PDBs — run 'vast_launch.py upload' first"
ls /workspace/data/ligands/*.sdf 2>/dev/null || echo "  No SDFs — run 'vast_launch.py upload' first"

echo ""
echo "============================================"
echo "  Setup complete! $(date)"
echo "============================================"
echo ""
echo "Run pipeline:"
echo "  cd /app/RFdiffusion && python scripts/run_inference.py \\"
echo "    inference.output_prefix=/workspace/results/rfdiffusion/erap2_short \\"
echo "    inference.input_pdb=/workspace/data/structures/erap2_wt_alphafold.pdb \\"
echo "    'contigmap.contigs=[A370-393/0 30-40]' \\"
echo "    'ppi.hotspot_res=[A370,A371,A374,A392,A393]' \\"
echo "    inference.num_designs=10 \\"
echo "    inference.ckpt_override_path=/app/RFdiffusion/models/Complex_beta_ckpt.pt"
