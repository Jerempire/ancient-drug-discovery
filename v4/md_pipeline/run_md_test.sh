#!/bin/bash
# MD Retention Test: Subprocess Switch Pattern
# Step 1: PyRosetta refines DiffPepDock poses (locks CUDA context)
# Step 2: OpenMM runs MD on GPU (fresh CUDA context)
# The PDB file is the handoff between processes.
set -e

echo "============================================"
echo "  MD Retention Test — Subprocess Switch"
echo "  $(date)"
echo "============================================"

DOCK_DIR=$(ls -d /workspace/DiffPepBuilder/runs/docking/*D_*M_*Y_* 2>/dev/null | tail -1)
RESULTS=/workspace/results
mkdir -p $RESULTS

# Find best poses for parent peptide
echo ""
echo "=== Finding best docked poses ==="
python3 << 'FINDPOSE'
import os, glob, json
import numpy as np

DOCK_DIR = os.environ.get("DOCK_DIR", "")
if not DOCK_DIR:
    dirs = sorted(glob.glob("/workspace/DiffPepBuilder/runs/docking/*D_*M_*Y_*"))
    DOCK_DIR = dirs[-1] if dirs else ""

def find_best(target, peptide):
    pep_path = os.path.join(DOCK_DIR, target, peptide)
    pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))[:16]
    best_pdb, best_c = None, 0
    for pdb in pdbs:
        pep, rec = [], []
        with open(pdb) as f:
            for line in f:
                if not line.startswith("ATOM"): continue
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if line[21] == "A": pep.append(xyz)
                elif line[21] == "B": rec.append(xyz)
        if not pep or not rec: continue
        diff = np.array(pep)[:, None, :] - np.array(rec)[None, :, :]
        c = int(np.sum(np.sqrt(np.sum(diff**2, axis=-1)) < 4.5))
        if c > best_c: best_c, best_pdb = c, pdb
    return best_pdb, best_c

# Try multiple peptide names
for pep in ["v41c_parent", "pep_glu_long_01"]:
    k392, kc = find_best("erap2_k392", pep)
    irap, ic = find_best("irap", pep)
    if k392 and irap:
        print(f"K392: {k392} ({kc} contacts)")
        print(f"IRAP: {irap} ({ic} contacts)")
        with open("/workspace/results/best_poses.json", "w") as f:
            json.dump({"k392_pdb": k392, "irap_pdb": irap, "peptide": pep}, f)
        break
else:
    print("ERROR: No poses found")
    exit(1)
FINDPOSE

# Read best poses
K392_PDB=$(python3 -c "import json; d=json.load(open('/workspace/results/best_poses.json')); print(d['k392_pdb'])")
IRAP_PDB=$(python3 -c "import json; d=json.load(open('/workspace/results/best_poses.json')); print(d['irap_pdb'])")

echo "K392 pose: $K392_PDB"
echo "IRAP pose: $IRAP_PDB"

# --- STEP 1: PyRosetta refine (separate process) ---
echo ""
echo "=== Step 1: PyRosetta Refine ==="
echo "  Refining K392 complex..."
python3 /workspace/md_pipeline/refine.py "$K392_PDB" "$RESULTS/k392_refined.pdb"

echo "  Refining IRAP complex..."
python3 /workspace/md_pipeline/refine.py "$IRAP_PDB" "$RESULTS/irap_refined.pdb"

echo "  Refinement complete. PyRosetta process ended."

# --- STEP 2: OpenMM MD (fresh process, clean CUDA context) ---
echo ""
echo "=== Step 2: OpenMM MD (CUDA) ==="
echo "  Running K392 MD (2ns)..."
python3 /workspace/md_pipeline/simulate.py "$RESULTS/k392_refined.pdb" "erap2_k392" "$RESULTS/md_k392.json"

echo "  Running IRAP MD (2ns)..."
python3 /workspace/md_pipeline/simulate.py "$RESULTS/irap_refined.pdb" "irap" "$RESULTS/md_irap.json"

# --- STEP 3: Compare ---
echo ""
echo "=== Comparison ==="
python3 << 'COMPARE'
import json

k = json.load(open("/workspace/results/md_k392.json"))
i = json.load(open("/workspace/results/md_irap.json"))

print(f"{'Metric':>25s}  {'ERAP2-K392':>12s}  {'IRAP':>12s}")
print("-" * 55)
for m in ["final_rmsd", "avg_rmsd_last_1ns", "max_rmsd", "rmsd_drift", "final_dist", "dist_drift"]:
    print(f"{m:>25s}  {k[m]:>12.2f}  {i[m]:>12.2f}")

kd = k["rmsd_drift"]
id_ = i["rmsd_drift"]
if id_ > kd * 2:
    v = "IRAP DRIFTS — binding is transient (GOOD)"
elif id_ > kd * 1.5:
    v = "IRAP less stable — weaker binding"
else:
    v = "Similar stability — IRAP binding persistent"
print(f"\nVERDICT: {v}")

combined = {"erap2_k392": k, "irap": i, "verdict": v}
with open("/workspace/results/md_retention_final.json", "w") as f:
    json.dump(combined, f, indent=2)
print("Saved: /workspace/results/md_retention_final.json")
COMPARE

echo ""
echo "============================================"
echo "  COMPLETE: $(date)"
echo "============================================"
