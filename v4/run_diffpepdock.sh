#!/bin/bash
# V4: DiffPepDock peptide-channel docking for K392/N392 selectivity
# Runs on Vast.ai with PyTorch image
set -e

echo "============================================"
echo "  V4: DiffPepDock K392/N392 Selectivity"
echo "  $(date)"
echo "============================================"

# --- STEP 1: Install DiffPepDock ---
echo "=== Installing DiffPepDock ==="
cd /workspace
if [ ! -d "DiffPepBuilder" ]; then
    git clone https://github.com/YuzheWangPKU/DiffPepBuilder.git
    cd DiffPepBuilder
    pip install -q -e . 2>&1 | tail -5
    # Download model weights
    mkdir -p experiments/checkpoints
    wget -q https://zenodo.org/records/15398020/files/diffpepdock_v1.pth -O experiments/checkpoints/diffpepdock_v1.pth
    echo "DiffPepDock installed, weights downloaded"
else
    cd DiffPepBuilder
    echo "DiffPepDock already installed"
fi

# Install additional deps
pip install -q biopython 2>/dev/null

# --- STEP 2: Prepare ERAP2 K392 and N392 structures ---
echo ""
echo "=== Preparing structures ==="
cd /workspace

# Create N392 variant from wildtype K392 PDB
python3 << 'PREP'
import os

# Read wildtype PDB
with open("/workspace/data/erap2_wt_alphafold.pdb") as f:
    lines = f.readlines()

# K392 version (wildtype) - just copy
with open("/workspace/docking_data/erap2_k392.pdb", "w") as f:
    f.writelines(lines)
print("K392 PDB ready")

# N392 version - mutate K392 to N
# In AlphaFold PDB, residue 392 is LYS chain A
# Replace LYS with ASN at position 392 (keep backbone, truncate sidechain)
n392_lines = []
for line in lines:
    if line.startswith("ATOM") and int(line[22:26].strip()) == 392:
        # Keep backbone atoms (N, CA, C, O, CB)
        atom_name = line[12:16].strip()
        if atom_name in ("N", "CA", "C", "O", "CB"):
            # Change residue name to ASN
            new_line = line[:17] + "ASN" + line[20:]
            n392_lines.append(new_line)
        elif atom_name == "CG":
            new_line = line[:17] + "ASN" + line[20:]
            n392_lines.append(new_line)
        elif atom_name == "CD":
            # Rename to OD1 for ASN
            new_line = line[:12] + " OD1" + line[16:17] + "ASN" + line[20:]
            n392_lines.append(new_line)
        elif atom_name == "CE":
            # Rename to ND2 for ASN
            new_line = line[:12] + " ND2" + line[16:17] + "ASN" + line[20:]
            n392_lines.append(new_line)
        # Skip NZ (lysine only)
    else:
        n392_lines.append(line)

with open("/workspace/docking_data/erap2_n392.pdb", "w") as f:
    f.writelines(n392_lines)
print("N392 PDB ready")
PREP

# --- STEP 3: Preprocess for DiffPepDock ---
echo ""
echo "=== Preprocessing ==="
cd /workspace/DiffPepBuilder

python experiments/process_batch_dock.py \
    --pdb_dir /workspace/docking_data \
    --write_dir /workspace/dock_processed \
    --receptor_info_path /workspace/docking_data/docking_cases.json \
    --peptide_seq_path /workspace/docking_data/peptide_seq.fasta \
    2>&1 | tail -10

# --- STEP 4: Run docking ---
echo ""
echo "=== Docking (this takes a while) ==="

torchrun --nproc-per-node=1 experiments/run_docking.py \
    data.val_csv_path=/workspace/dock_processed/metadata_test.csv \
    2>&1 | tail -20

echo ""
echo "=== Docking complete ==="
ls runs/docking/*.pdb 2>/dev/null | wc -l
echo "PDB files generated"

# --- STEP 5: Score and rank ---
echo ""
echo "=== Scoring ==="

python3 << 'SCORE'
import glob
import json
import os
from pathlib import Path

results = []
dock_dir = "runs/docking"

# Check for postprocess results
csv_path = os.path.join(dock_dir, "postprocess_results.csv")
if os.path.exists(csv_path):
    import csv
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            results.append(row)
    print(f"Loaded {len(results)} results from postprocess CSV")
else:
    # Parse PDB filenames for scores
    pdbs = glob.glob(os.path.join(dock_dir, "*.pdb"))
    print(f"Found {len(pdbs)} docked PDBs (no postprocess CSV)")

# Group by peptide and target
if results:
    k392_scores = {}
    n392_scores = {}
    for r in results:
        pep = r.get("peptide_name", r.get("name", "?"))
        target = r.get("target", "?")
        score = float(r.get("ddG", r.get("score", 0)))
        if "k392" in target.lower():
            k392_scores[pep] = score
        elif "n392" in target.lower():
            n392_scores[pep] = score

    print(f"\n{'Peptide':>20s}  {'K392':>8s}  {'N392':>8s}  {'Delta':>8s}  Verdict")
    print("-" * 60)
    combined = []
    for pep in sorted(set(list(k392_scores.keys()) + list(n392_scores.keys()))):
        k = k392_scores.get(pep, 0)
        n = n392_scores.get(pep, 0)
        delta = k - n
        verdict = "K392-SEL" if delta < -1.0 else "N392-SEL" if delta > 1.0 else "NEUTRAL"
        combined.append({"peptide": pep, "k392_ddG": k, "n392_ddG": n, "delta": delta, "verdict": verdict})
        print(f"{pep:>20s}  {k:>8.2f}  {n:>8.2f}  {delta:>+8.2f}  {verdict}")

    combined.sort(key=lambda x: x["delta"])
    with open("/workspace/results/diffpepdock_results.json", "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nSaved: /workspace/results/diffpepdock_results.json")

print("\nDone!")
SCORE

echo ""
echo "============================================"
echo "  V4 COMPLETE: $(date)"
echo "============================================"
