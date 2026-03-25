#!/bin/bash
# Short-Stack Grind: 6-7mer N392-selective peptides
# Full pipeline: DiffPepDock → PyRosetta refine → OpenMM 5ns MD → Residence time
set -e

echo "============================================"
echo "  Short-Stack Grind: 6-7mer N392 Selectivity"
echo "  $(date)"
echo "============================================"

# --- STEP 1: Install deps ---
echo "=== Installing dependencies ==="
pip install -q biopython pandas pyrootutils hydra-core omegaconf fair-esm easydict \
  tmtools biotite mdtraj svgwrite ml-collections dm-tree wget pdbfixer wandb \
  scikit-learn scipy tqdm GPUtil 2>&1 | tail -1

conda install -y -c https://conda.rosettacommons.org pyrosetta 2>&1 | tail -2
conda install -y -c conda-forge openmm cuda-version=12.4 2>&1 | tail -2
pip install -q 'numpy<2' pdbfixer 2>&1 | tail -1

# Verify
python3 -c 'import pyrosetta; print("PyRosetta OK")'
python3 -c '
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    if p.getName() == "CUDA": print(f"OpenMM CUDA OK (speed={p.getSpeed()})")
'
echo "DEPS DONE"

# --- STEP 2: Prepare structures ---
echo ""
echo "=== Preparing structures ==="
python3 /workspace/prep_structures.py

# --- STEP 3: DiffPepDock ---
echo ""
echo "=== Setting up DiffPepDock ==="
cd /workspace
git clone -q https://github.com/YuzheWangPKU/DiffPepBuilder.git 2>/dev/null || true
cd DiffPepBuilder && pip install -q -e . 2>&1 | tail -1
mkdir -p experiments/checkpoints
[ -f experiments/checkpoints/diffpepdock_v1.pth ] || \
  wget -q https://zenodo.org/records/15398020/files/diffpepdock_v1.pth \
  -O experiments/checkpoints/diffpepdock_v1.pth

# Preprocess
rm -rf /workspace/dock_processed
python experiments/process_batch_dock.py \
  --pdb_dir /workspace/docking_data \
  --write_dir /workspace/dock_processed \
  --receptor_info_path /workspace/docking_data/docking_cases.json \
  --peptide_seq_path /workspace/docking_data/peptide_library.fasta \
  2>&1 | tail -5

# Dock
echo ""
echo "=== Docking ==="
export BASE_PATH=/workspace/DiffPepBuilder WANDB_MODE=disabled HYDRA_FULL_ERROR=1
export LOCAL_RANK=0 RANK=0 WORLD_SIZE=1 MASTER_ADDR=localhost MASTER_PORT=29500
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python experiments/run_docking.py \
  data.val_csv_path=/workspace/dock_processed/metadata_test.csv \
  data.filtering.max_len=1200 \
  experiment.num_gpus=1 \
  experiment.use_ddp=False \
  experiment.eval_batch_size=8 \
  2>&1 | tail -5

DOCK_DIR=$(ls -d /workspace/DiffPepBuilder/runs/docking/*D_*M_*Y_* | tail -1)
echo "Docked PDBs: $(find $DOCK_DIR -name '*.pdb' | wc -l)"

# --- STEP 4: Contact scoring ---
echo ""
echo "=== Contact scoring ==="
python3 /workspace/contact_score.py

# --- STEP 5: MD retention on top hits ---
echo ""
echo "=== MD Retention Test (subprocess switch) ==="
python3 /workspace/select_and_md.py

echo ""
echo "============================================"
echo "  COMPLETE: $(date)"
echo "============================================"
