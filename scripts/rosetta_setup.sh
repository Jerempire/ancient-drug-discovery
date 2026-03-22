#!/bin/bash
# rosetta_setup.sh — Setup for Vast.ai: Boltz-2 + PyRosetta
#
# Image: pytorch/pytorch:2.5.1-cuda12.4-cudnn9-devel
#
# Usage:
#   ssh -p <PORT> root@<HOST> "cd /workspace && bash scripts/rosetta_setup.sh"
#
# Expected runtime: ~5-8 min (PyRosetta is ~2 GB)

set -euo pipefail

echo "============================================"
echo "  Rosetta Interface Analysis — GPU Setup"
echo "  Boltz-2 + PyRosetta"
echo "============================================"
echo "Started: $(date)"
echo ""

# --- Install Boltz-2 ---
echo "[1/4] Installing Boltz-2..."
pip install -q boltz biopython pandas 2>&1 | tail -5
echo "  Done."

# --- Install PyRosetta via conda ---
echo ""
echo "[2/4] Installing PyRosetta (this takes a few minutes)..."
conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta 2>&1 | tail -10
echo "  Done."

# --- Install gemmi (CIF conversion fallback) ---
echo ""
echo "[3/4] Installing gemmi..."
pip install -q gemmi 2>&1 | tail -3
echo "  Done."

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

python3 -c "from pyrosetta import init; init('-mute all'); print('PyRosetta: OK')" 2>&1 | tail -1
python3 -c "import boltz; print('Boltz-2: OK')" 2>&1 | tail -1
python3 -c "import gemmi; print('gemmi: OK')" 2>&1 | tail -1
python3 -c "from Bio.PDB import PDBParser; print('BioPython: OK')" 2>&1 | tail -1

echo ""
echo "============================================"
echo "  Setup complete! $(date)"
echo "============================================"
echo ""
echo "Next: python3 /workspace/scripts/run_rosetta_pipeline.py"
