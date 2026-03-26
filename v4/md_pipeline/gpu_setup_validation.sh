#!/bin/bash
# Setup script for MD + MMGBSA + EvoBind validation on Vast.ai
# Run after instance is up: bash /workspace/md_pipeline/gpu_setup_validation.sh
set -e

echo "============================================"
echo "  Validation Pipeline Setup"
echo "  $(date)"
echo "============================================"

# ---- Part 1: MD + MMGBSA deps ----
echo ""
echo "=== [1/4] Installing MD/MMGBSA dependencies ==="
pip install gemmi parmed mdtraj 2>/dev/null

python3 -c "
from openmm import Platform
p = Platform.getPlatformByName('CUDA')
print(f'  OpenMM CUDA OK (speed={p.getSpeed()})')
"
python3 -c "import parmed; print(f'  parmed {parmed.__version__} OK')"
python3 -c "import gemmi; print(f'  gemmi OK')"

# ---- Part 2: EvoBind clone + deps ----
echo ""
echo "=== [2/4] Setting up EvoBind ==="
if [ ! -d /workspace/evobind ]; then
    cd /workspace
    git clone https://github.com/patrickbryant1/EvoBind.git evobind
    cd /workspace
fi

# Install EvoBind deps into existing env (avoid full conda env — conflicts with PyTorch)
# Core deps that aren't already in the PyTorch image:
pip install \
    dm-haiku==0.0.12 \
    dm-tree==0.1.8 \
    chex==0.1.86 \
    ml-collections==0.1.1 \
    immutabledict==4.2.0 \
    jmp==0.0.4 \
    biopython \
    tabulate \
    2>/dev/null

# JAX with CUDA 12 (EvoBind needs JAX, not just PyTorch)
pip install "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html 2>/dev/null

python3 -c "import jax; print(f'  JAX {jax.__version__} OK, devices: {jax.devices()}')"

# ---- Part 3: AlphaFold2 parameters ----
echo ""
echo "=== [3/4] Downloading AlphaFold2 parameters ==="
AF2_PARAMS=/workspace/evobind/src/AF2/params
mkdir -p $AF2_PARAMS

if [ ! -f "$AF2_PARAMS/params_model_1.npz" ]; then
    echo "  Downloading AF2 params (~3.5GB)..."
    cd $AF2_PARAMS
    wget -q https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
    tar xf alphafold_params_2021-07-14.tar
    rm alphafold_params_2021-07-14.tar
    echo "  AF2 params installed"
else
    echo "  AF2 params already present"
fi

# ---- Part 4: ERAP2 MSA via ColabFold API (replaces 74GB Uniclust30) ----
echo ""
echo "=== [4/4] Generating ERAP2 MSA via ColabFold API ==="
ERAP2_MSA=/workspace/evobind/data/erap2_msa
mkdir -p $ERAP2_MSA

if [ ! -f "$ERAP2_MSA/erap2_k392.a3m" ]; then
    python3 /workspace/md_pipeline/evobind_msa.py \
        /workspace/md_pipeline/erap2_k392.fasta \
        $ERAP2_MSA/erap2_k392.a3m
else
    echo "  ERAP2 MSA already generated"
fi

echo ""
echo "============================================"
echo "  Setup complete! $(date)"
echo ""
echo "  Run MD+MMGBSA:  bash /workspace/md_pipeline/run_elite_four.sh"
echo "  Run EvoBind:    bash /workspace/md_pipeline/run_evobind_erap2.sh"
echo "============================================"
