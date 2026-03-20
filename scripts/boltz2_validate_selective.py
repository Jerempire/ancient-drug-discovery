"""
boltz2_validate_selective.py — Validate top ERAP2-selective binders with Boltz-2.

For each top ProteinMPNN sequence:
  1. Predict ERAP2 + binder complex structure (confidence = binding quality)
  2. Predict ERAP1 + binder complex structure (cross-reactivity check)
  3. Compare ipTM/pTM scores — higher ERAP2 vs ERAP1 = more selective

Uses Boltz-2 via CLI: `boltz predict input.yaml`

Outputs:
  /workspace/results/boltz2_validation/
    erap2_<design>/  — Boltz-2 output for ERAP2 complex
    erap1_<design>/  — Boltz-2 output for ERAP1 complex
    validation_summary.json — ranked results
"""
import json
import os
import subprocess
import sys
import glob
import time

RESULTS_DIR = "/workspace/results/boltz2_validation"
MPNN_SUMMARY = "/workspace/results/mpnn/mpnn_summary.json"

# ERAP2 and ERAP1 sequences (first 200 residues of substrate channel region)
# We use truncated target sequences around the binding site for speed
# Full sequence would take too long per prediction

# Read full sequences from structures
from Bio.PDB import PDBParser

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

parser = PDBParser(QUIET=True)


def get_region_sequence(pdb_path, start, end):
    """Extract sequence for a residue range from PDB."""
    struct = parser.get_structure("s", pdb_path)
    seq = []
    for r in struct.get_residues():
        resnum = r.get_id()[1]
        resname = r.get_resname()
        if resname in AA3TO1 and start <= resnum <= end:
            seq.append(AA3TO1[resname])
    return "".join(seq)


def get_full_sequence(pdb_path):
    """Get full protein sequence from PDB."""
    struct = parser.get_structure("s", pdb_path)
    seq = []
    for r in struct.get_residues():
        if r.get_resname() in AA3TO1:
            seq.append(AA3TO1[r.get_resname()])
    return "".join(seq)


# Get target region sequences (channel region 370-500 to cover binding context)
erap2_region = get_region_sequence(
    "/workspace/data/structures/erap2_wt_alphafold.pdb", 370, 500)
erap1_region = get_region_sequence(
    "/workspace/data/structures/erap1_wt_alphafold.pdb", 370, 500)

print("ERAP2 region (370-500): %d aa" % len(erap2_region))
print("ERAP1 region (370-500): %d aa" % len(erap1_region))


def write_boltz_yaml(target_seq, binder_seq, output_path, name="complex"):
    """Write Boltz-2 input YAML for a protein complex."""
    yaml_content = """version: 1
sequences:
  - protein:
      id: A
      sequence: %s
  - protein:
      id: B
      sequence: %s
""" % (target_seq, binder_seq)

    with open(output_path, "w") as f:
        f.write(yaml_content)


def run_boltz(yaml_path, output_dir):
    """Run Boltz-2 prediction."""
    cmd = [
        "/opt/conda/envs/bioreason/bin/boltz", "predict",
        yaml_path,
        "--out_dir", output_dir,
        "--recycling_steps", "3",
        "--diffusion_samples", "1",
        "--accelerator", "gpu",
        "--devices", "1",
        "--use_msa_server",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return result


def parse_boltz_scores(output_dir):
    """Parse Boltz-2 confidence scores from output."""
    # Look for confidence JSON
    json_files = glob.glob(os.path.join(output_dir, "**", "confidence_*.json"), recursive=True)
    if not json_files:
        json_files = glob.glob(os.path.join(output_dir, "**", "*.json"), recursive=True)

    scores = {"ptm": None, "iptm": None, "ranking_score": None}

    for jf in json_files:
        try:
            with open(jf) as f:
                data = json.load(f)
            if isinstance(data, dict):
                for key in ["ptm", "iptm", "ranking_score", "complex_plddt",
                            "pair_chains_iptm", "aggregate_score"]:
                    if key in data:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue

    # Also check for scores in CIF/PDB filenames or log
    log_files = glob.glob(os.path.join(output_dir, "**", "*.log"), recursive=True)
    for lf in log_files:
        try:
            with open(lf) as f:
                for line in f:
                    if "iptm" in line.lower() or "ptm" in line.lower():
                        scores["log_line"] = line.strip()
        except IOError:
            pass

    return scores


# --- Main ---
os.makedirs(RESULTS_DIR, exist_ok=True)

# Load top MPNN sequences
with open(MPNN_SUMMARY) as f:
    mpnn_results = json.load(f)

# Take top 5 unique designs (best score per design)
seen_designs = set()
top_sequences = []
for r in mpnn_results:
    if r["design"] not in seen_designs and len(top_sequences) < 5:
        seen_designs.add(r["design"])
        top_sequences.append(r)

print("\nValidating %d top binder sequences..." % len(top_sequences))

validation_results = []

for i, seq_info in enumerate(top_sequences):
    design = seq_info["design"]
    binder_seq = seq_info["sequence"]
    mpnn_score = seq_info.get("score", "?")

    print("\n--- [%d/%d] %s (MPNN=%.3f, %d aa) ---" % (
        i + 1, len(top_sequences), design,
        mpnn_score if isinstance(mpnn_score, float) else 0,
        len(binder_seq)))

    result = {
        "design": design,
        "binder_sequence": binder_seq,
        "binder_length": len(binder_seq),
        "mpnn_score": mpnn_score,
    }

    # --- ERAP2 complex ---
    e2_dir = os.path.join(RESULTS_DIR, "erap2_%s" % design)
    e2_yaml = os.path.join(RESULTS_DIR, "erap2_%s.yaml" % design)
    write_boltz_yaml(erap2_region, binder_seq, e2_yaml)

    print("  ERAP2 complex...", end=" ", flush=True)
    t0 = time.time()
    r = run_boltz(e2_yaml, e2_dir)
    elapsed = time.time() - t0

    if r.returncode == 0:
        e2_scores = parse_boltz_scores(e2_dir)
        result["erap2_scores"] = e2_scores
        print("OK (%.0fs) ptm=%s iptm=%s" % (
            elapsed, e2_scores.get("ptm"), e2_scores.get("iptm")))
    else:
        print("FAILED (%.0fs)" % elapsed)
        result["erap2_error"] = r.stderr[-200:] if r.stderr else "unknown"
        # Print last bit of output for debugging
        if r.stdout:
            print("  stdout: %s" % r.stdout[-300:])
        if r.stderr:
            print("  stderr: %s" % r.stderr[-300:])

    # --- ERAP1 complex (cross-reactivity) ---
    e1_dir = os.path.join(RESULTS_DIR, "erap1_%s" % design)
    e1_yaml = os.path.join(RESULTS_DIR, "erap1_%s.yaml" % design)
    write_boltz_yaml(erap1_region, binder_seq, e1_yaml)

    print("  ERAP1 complex...", end=" ", flush=True)
    t0 = time.time()
    r = run_boltz(e1_yaml, e1_dir)
    elapsed = time.time() - t0

    if r.returncode == 0:
        e1_scores = parse_boltz_scores(e1_dir)
        result["erap1_scores"] = e1_scores
        print("OK (%.0fs) ptm=%s iptm=%s" % (
            elapsed, e1_scores.get("ptm"), e1_scores.get("iptm")))
    else:
        print("FAILED (%.0fs)" % elapsed)
        result["erap1_error"] = r.stderr[-200:] if r.stderr else "unknown"
        if r.stdout:
            print("  stdout: %s" % r.stdout[-300:])
        if r.stderr:
            print("  stderr: %s" % r.stderr[-300:])

    # --- Selectivity ---
    e2_iptm = result.get("erap2_scores", {}).get("iptm")
    e1_iptm = result.get("erap1_scores", {}).get("iptm")
    if e2_iptm is not None and e1_iptm is not None and e1_iptm > 0:
        result["iptm_selectivity"] = round(e2_iptm / e1_iptm, 3)
        result["iptm_delta"] = round(e2_iptm - e1_iptm, 4)
    else:
        result["iptm_selectivity"] = None
        result["iptm_delta"] = None

    validation_results.append(result)

# Sort by ERAP2 iptm (higher = better binding)
validation_results.sort(
    key=lambda x: x.get("erap2_scores", {}).get("iptm") or 0, reverse=True)

# Save
summary_path = os.path.join(RESULTS_DIR, "validation_summary.json")
with open(summary_path, "w") as f:
    json.dump(validation_results, f, indent=2)

# Report
print("\n" + "=" * 80)
print("BOLTZ-2 VALIDATION — ERAP2 Selective Binders")
print("=" * 80)

print("\n%-30s %6s %6s %6s %8s %8s" % (
    "Design", "MPNN", "E2_iPTM", "E1_iPTM", "Sel_Ratio", "Tag"))
print("-" * 90)

for r in validation_results:
    e2i = r.get("erap2_scores", {}).get("iptm")
    e1i = r.get("erap1_scores", {}).get("iptm")
    sel = r.get("iptm_selectivity")
    mpnn = r.get("mpnn_score")

    e2_str = "%.3f" % e2i if e2i is not None else "  N/A"
    e1_str = "%.3f" % e1i if e1i is not None else "  N/A"
    sel_str = "%.2fx" % sel if sel is not None else "  N/A"
    mpnn_str = "%.3f" % mpnn if mpnn is not None else "N/A"

    tag = ""
    if sel is not None:
        if sel > 1.3:
            tag = "SELECTIVE"
        elif sel > 1.1:
            tag = "MODERATE"
        else:
            tag = "NON-SEL"

    print("%-30s %6s %6s %6s %8s %8s" % (
        r["design"][:30], mpnn_str, e2_str, e1_str, sel_str, tag))

print("\nSaved: %s" % summary_path)
