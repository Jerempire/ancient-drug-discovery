#!/bin/bash
# V4 Full Counter-Screen: DiffPepDock + Pose Quality Analysis
# Top 3 peptides x 4 targets (K392, N392, ERAP1, IRAP)
# Run on Vast.ai RTX 3090
set -e

echo "============================================"
echo "  V4 Full Counter-Screen + Pose Analysis"
echo "  $(date)"
echo "============================================"

# --- STEP 1: Install everything ---
echo "=== Installing dependencies ==="
pip install -q biopython pandas pyrootutils hydra-core omegaconf fair-esm easydict \
  tmtools biotite mdtraj svgwrite ml-collections dm-tree wget pdbfixer wandb \
  scikit-learn scipy tqdm openmm 2>&1 | tail -3

conda install -y -c https://conda.rosettacommons.org pyrosetta 2>&1 | tail -3

# --- STEP 2: Clone DiffPepDock ---
echo ""
echo "=== Setting up DiffPepDock ==="
cd /workspace
if [ ! -d "DiffPepBuilder" ]; then
    git clone -q https://github.com/YuzheWangPKU/DiffPepBuilder.git
    cd DiffPepBuilder
    pip install -q -e . 2>&1 | tail -2
    mkdir -p experiments/checkpoints
    wget -q https://zenodo.org/records/15398020/files/diffpepdock_v1.pth \
      -O experiments/checkpoints/diffpepdock_v1.pth
    echo "DiffPepDock ready"
else
    cd DiffPepBuilder
    echo "DiffPepDock already installed"
fi

# --- STEP 3: Prepare structures ---
echo ""
echo "=== Preparing structures ==="
cd /workspace
python3 /workspace/counterscreen_setup.py

# Also create ERAP2 K392 and N392 cropped structures
python3 << 'PREP'
from Bio.PDB import PDBParser, PDBIO, Select

class ChannelSelect(Select):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def accept_residue(self, residue):
        resnum = residue.get_id()[1]
        return self.start <= resnum <= self.end

parser = PDBParser(QUIET=True)

# K392 (wildtype)
s = parser.get_structure("erap2", "/workspace/data/erap2_wt_alphafold.pdb")
io = PDBIO()
io.set_structure(s)
io.save("/workspace/docking_data/erap2_k392.pdb", ChannelSelect(350, 450))
print("ERAP2-K392 cropped")

# N392 variant
with open("/workspace/docking_data/erap2_k392.pdb") as f:
    lines = f.readlines()
n392_lines = []
for line in lines:
    if line.startswith("ATOM") and line[22:26].strip() == "392":
        atom_name = line[12:16].strip()
        if atom_name in ("N", "CA", "C", "O", "CB"):
            n392_lines.append(line[:17] + "ASN" + line[20:])
        elif atom_name == "CG":
            n392_lines.append(line[:17] + "ASN" + line[20:])
        elif atom_name == "CD":
            n392_lines.append(line[:12] + " OD1" + line[16:17] + "ASN" + line[20:])
        elif atom_name == "CE":
            n392_lines.append(line[:12] + " ND2" + line[16:17] + "ASN" + line[20:])
    else:
        n392_lines.append(line)
with open("/workspace/docking_data/erap2_n392.pdb", "w") as f:
    f.writelines(n392_lines)
print("ERAP2-N392 cropped")
PREP

# Update docking config for all 4 targets
python3 << 'CONFIG'
import json
config = {
    "erap2_k392": {
        "pdb": "erap2_k392.pdb",
        "binding_site_residues": "388-410",
        "description": "ERAP2 K392 substrate channel"
    },
    "erap2_n392": {
        "pdb": "erap2_n392.pdb",
        "binding_site_residues": "388-410",
        "description": "ERAP2 N392 substrate channel"
    },
    "erap1": {
        "pdb": "erap1.pdb",
        "binding_site_residues": "388-410",
        "description": "ERAP1 substrate channel"
    },
    "irap": {
        "pdb": "irap.pdb",
        "binding_site_residues": "388-410",
        "description": "IRAP substrate channel"
    }
}
# Top 3 peptides only
fasta = """>pep_glu_long_01
EALVAAGLAGLA
>pep_asp_long_01
DALVAAGLAGLA
>pep_glu_e3_01
EAELAAGLAA
"""
with open("/workspace/docking_data/docking_cases.json", "w") as f:
    json.dump(config, f, indent=2)
with open("/workspace/docking_data/peptide_seq.fasta", "w") as f:
    f.write(fasta)
print("Config: 4 targets x 3 peptides = 12 cases")
CONFIG

# --- STEP 4: Preprocess ---
echo ""
echo "=== Preprocessing ==="
cd /workspace/DiffPepBuilder
python experiments/process_batch_dock.py \
    --pdb_dir /workspace/docking_data \
    --write_dir /workspace/dock_processed \
    --receptor_info_path /workspace/docking_data/docking_cases.json \
    --peptide_seq_path /workspace/docking_data/peptide_seq.fasta \
    2>&1 | tail -5

# --- STEP 5: Dock ---
echo ""
echo "=== Docking (12 cases x 32 samples = 384 PDBs) ==="
export BASE_PATH=/workspace/DiffPepBuilder
export WANDB_MODE=disabled
export HYDRA_FULL_ERROR=1
export LOCAL_RANK=0 RANK=0 WORLD_SIZE=1 MASTER_ADDR=localhost MASTER_PORT=29500
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python experiments/run_docking.py \
    data.val_csv_path=/workspace/dock_processed/metadata_test.csv \
    data.filtering.max_len=1200 \
    experiment.num_gpus=1 \
    experiment.use_ddp=False \
    experiment.eval_batch_size=8 \
    2>&1 | tail -5

echo ""
echo "=== Docking complete ==="
find /workspace/DiffPepBuilder/runs/docking/ -name "*.pdb" | wc -l
echo "PDB files generated"

# --- STEP 6: Pose quality analysis ---
echo ""
echo "=== Pose Quality Analysis ==="
python /workspace/pose_quality_analysis.py

echo ""
echo "============================================"
echo "  COMPLETE: $(date)"
echo "============================================"
