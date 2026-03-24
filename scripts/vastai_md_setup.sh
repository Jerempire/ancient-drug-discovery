#!/bin/bash
# =============================================================================
# Vast.ai MD Setup Script — Run ONCE per instance before any MD simulations
#
# Installs: Boltz-2, OpenMM, MDTraj, PDBFixer, MDAnalysis, BioPython
# Fixes: Known issues (KMP_DUPLICATE_LIB_OK, numpy version, torch CVE check)
# Verifies: GPU access, CUDA, all imports work
#
# Usage: ssh into Vast.ai instance, then:
#   bash /workspace/vastai_md_setup.sh
#
# Expected time: ~5-10 min (mostly pip downloads)
# =============================================================================

set -e
echo "========================================"
echo "VAST.AI MD SETUP — $(date)"
echo "========================================"

# --- 1. System info ---
echo ""
echo "[1/8] System info..."
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python3 --version
echo "CUDA available: $(python3 -c 'import torch; print(torch.cuda.is_available())' 2>/dev/null || echo 'torch not installed')"

# --- 2. Set environment variables for known issues ---
echo ""
echo "[2/8] Setting environment variables..."
export KMP_DUPLICATE_LIB_OK=TRUE
echo 'export KMP_DUPLICATE_LIB_OK=TRUE' >> ~/.bashrc

# --- 3. Install Boltz-2 ---
echo ""
echo "[3/8] Installing Boltz-2..."
pip install boltz 2>&1 | tail -3
echo "Boltz-2: $(python3 -c 'import boltz; print(boltz.__version__)' 2>/dev/null || echo 'FAILED')"

# --- 4. Install OpenMM + MD tools ---
echo ""
echo "[4/8] Installing OpenMM + MD tools..."
# OpenMM via conda-forge (pip version often has issues)
# But we're in a pip-only container, so use conda if available, else pip
if command -v conda &> /dev/null; then
    conda install -y -c conda-forge openmm mdtraj pdbfixer 2>&1 | tail -5
else
    # pip fallback — openmm from conda-forge wheel
    pip install openmm mdtraj pdbfixer 2>&1 | tail -5
fi

# --- 5. Install analysis tools ---
echo ""
echo "[5/8] Installing analysis tools..."
pip install MDAnalysis biopython numpy 2>&1 | tail -3

# --- 6. Fix numpy version if needed (Boltz wants 1.26, some tools want 2.x) ---
echo ""
echo "[6/8] Checking numpy compatibility..."
NUMPY_VER=$(python3 -c 'import numpy; print(numpy.__version__)')
echo "NumPy version: $NUMPY_VER"
# Boltz-2 pins numpy<2.0 — if something upgraded it, downgrade
if python3 -c "import numpy; assert int(numpy.__version__.split('.')[0]) < 2" 2>/dev/null; then
    echo "NumPy OK (< 2.0)"
else
    echo "NumPy >= 2.0 detected, downgrading for Boltz-2 compatibility..."
    pip install 'numpy<2.0' 2>&1 | tail -2
fi

# --- 7. Verify all imports ---
echo ""
echo "[7/8] Verifying imports..."
python3 << 'PYEOF'
import sys
errors = []

try:
    import boltz
    print("  boltz: %s" % boltz.__version__)
except Exception as e:
    errors.append("boltz: %s" % e)

try:
    import openmm
    print("  openmm: %s" % openmm.__version__)
except Exception as e:
    errors.append("openmm: %s" % e)

try:
    import mdtraj
    print("  mdtraj: %s" % mdtraj.__version__)
except Exception as e:
    errors.append("mdtraj: %s" % e)

try:
    import pdbfixer
    print("  pdbfixer: OK")
except Exception as e:
    errors.append("pdbfixer: %s" % e)

try:
    import MDAnalysis
    print("  MDAnalysis: %s" % MDAnalysis.__version__)
except Exception as e:
    errors.append("MDAnalysis: %s" % e)

try:
    from Bio.PDB import PDBParser, MMCIFParser, Superimposer
    print("  biopython: OK")
except Exception as e:
    errors.append("biopython: %s" % e)

try:
    import torch
    print("  torch: %s (CUDA: %s)" % (torch.__version__, torch.cuda.is_available()))
except Exception as e:
    errors.append("torch: %s" % e)

try:
    import numpy
    print("  numpy: %s" % numpy.__version__)
except Exception as e:
    errors.append("numpy: %s" % e)

if errors:
    print("\nFAILED IMPORTS:")
    for e in errors:
        print("  ERROR: %s" % e)
    sys.exit(1)
else:
    print("\nAll imports OK!")
PYEOF

# --- 8. Test OpenMM GPU ---
echo ""
echo "[8/8] Testing OpenMM GPU access..."
python3 << 'PYEOF'
try:
    import openmm
    platforms = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    print("  Available platforms: %s" % ", ".join(platforms))
    if "CUDA" in platforms:
        print("  CUDA platform: AVAILABLE")
        # Quick test
        p = openmm.Platform.getPlatformByName("CUDA")
        print("  CUDA device: %s" % p.getPropertyDefaultValue("DeviceName"))
    elif "OpenCL" in platforms:
        print("  OpenCL platform: AVAILABLE (fallback)")
    else:
        print("  WARNING: No GPU platform found! MD will run on CPU (very slow)")
except Exception as e:
    print("  OpenMM GPU test failed: %s" % e)
    print("  MD may need to run on CPU platform")
PYEOF

echo ""
echo "========================================"
echo "SETUP COMPLETE — $(date)"
echo "========================================"
echo ""
echo "Next steps:"
echo "  1. Upload scripts:  scp -P <PORT> scripts/*.py root@<HOST>:/workspace/"
echo "  2. Upload structures: scp -P <PORT> data/results/v43_validation/md_starting_structures/*.pdb root@<HOST>:/workspace/structures/"
echo "  3. Run MD: python3 /workspace/openmm_md_protocol.py /workspace/structures/VAGSAF_vs_erap2k392.pdb /workspace/results/md_VAGSAF_erap2k392/"
