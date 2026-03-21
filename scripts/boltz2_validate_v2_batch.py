"""
boltz2_validate_v2_batch.py — Batch Boltz-2 validation for all V2 binder designs.

Validates all 30 V2 designs (short/medium/long × 10) through Boltz-2:
  1. Predict ERAP2 + binder complex (binding quality)
  2. Predict ERAP1 + binder complex (cross-reactivity)
  3. Compute selectivity metrics (ipTM delta, ratio)

Outputs validation_summary_v2.json in same format as validation_summary.json,
compatible with candidate_scorer.py.

Usage (on Vast.ai GPU instance):
    python3 /workspace/scripts/boltz2_validate_v2_batch.py
    python3 /workspace/scripts/boltz2_validate_v2_batch.py --designs short  # only short tier
    python3 /workspace/scripts/boltz2_validate_v2_batch.py --resume          # skip already done

Prereqs:
    - Boltz-2 installed (pip install boltz)
    - MPNN FASTA files at /workspace/results/mpnn_v2/erap2v2_*.fa
    - Target PDBs at /workspace/data/structures/
"""
import glob
import json
import os
import subprocess
import sys
import time

from Bio.PDB import PDBParser

# --- Config ---
WORKSPACE = "/workspace"
MPNN_DIR = os.path.join(WORKSPACE, "results", "mpnn_v2")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
OUTPUT_DIR = os.path.join(WORKSPACE, "results", "boltz2_validation_v2")

ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")

# Region around divergent channel for Boltz-2 (residues 350-500)
TARGET_REGION = (350, 500)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_region_sequence(pdb_path, start, end):
    """Extract sequence for a residue range from PDB."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    seq = []
    for r in struct.get_residues():
        resnum = r.get_id()[1]
        resname = r.get_resname()
        if resname in AA3TO1 and start <= resnum <= end:
            seq.append(AA3TO1[resname])
    return "".join(seq)


def parse_mpnn_fasta(fa_path):
    """Parse MPNN FASTA file. Returns list of (header, sequence, score) tuples.

    Takes the best-scoring designed sequence (lowest score = most designable).
    """
    entries = []
    current_header = None
    current_seq = []

    with open(fa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None and current_seq:
                    entries.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
    if current_header and current_seq:
        entries.append((current_header, "".join(current_seq)))

    # First entry is the wild-type backbone, skip it
    # Designed sequences start from index 1
    designed = []
    for header, seq in entries[1:]:
        score = None
        for part in header.split(","):
            part = part.strip()
            if part.startswith("score="):
                try:
                    score = float(part.split("=")[1])
                except ValueError:
                    pass
        designed.append({"header": header, "sequence": seq, "score": score})

    # Sort by score (lower = better)
    designed.sort(key=lambda x: x["score"] if x["score"] is not None else 999)
    return designed


def write_boltz_yaml(target_seq, binder_seq, output_path):
    """Write Boltz-2 input YAML for a protein complex (empty MSA for speed)."""
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {target_seq}
      msa: empty
  - protein:
      id: B
      sequence: {binder_seq}
      msa: empty
"""
    with open(output_path, "w") as f:
        f.write(yaml_content)


def run_boltz(yaml_path, output_dir):
    """Run Boltz-2 prediction. Returns subprocess result."""
    # Try common boltz install locations
    boltz_cmd = None
    for candidate in [
        "/opt/conda/envs/boltz/bin/boltz",
        "/opt/conda/envs/bioreason/bin/boltz",
        "/opt/conda/bin/boltz",
        "boltz",
    ]:
        try:
            subprocess.run([candidate, "--help"], capture_output=True, timeout=10)
            boltz_cmd = candidate
            break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue

    if boltz_cmd is None:
        print("ERROR: boltz not found. Install with: pip install boltz")
        sys.exit(1)

    cmd = [
        boltz_cmd, "predict",
        yaml_path,
        "--out_dir", output_dir,
        "--recycling_steps", "3",
        "--diffusion_samples", "1",
        "--accelerator", "gpu",
        "--devices", "1",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=900)


def parse_boltz_scores(output_dir):
    """Parse Boltz-2 confidence scores from output.

    Boltz-2 v2.x stores scores in:
      - confidence_*.json files
      - predictions/<name>/confidence_model_0.json
      - or as .npz files that need numpy
    """
    scores = {}
    json_files = glob.glob(os.path.join(output_dir, "**", "*.json"), recursive=True)

    for jf in json_files:
        if "manifest" in jf:
            continue
        try:
            with open(jf) as f:
                data = json.load(f)
            if isinstance(data, dict):
                for key in ["ptm", "iptm", "ranking_score", "complex_plddt",
                            "pair_chains_iptm", "confidence_score"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue

    # Boltz-2 v2.x also writes .npz confidence files
    if not scores:
        import numpy as np
        npz_files = glob.glob(os.path.join(output_dir, "**", "confidence_*.npz"), recursive=True)
        for npz_path in npz_files:
            try:
                data = dict(np.load(npz_path, allow_pickle=True))
                for key in ["ptm", "iptm", "complex_plddt", "pair_chains_iptm"]:
                    if key in data and key not in scores:
                        val = data[key]
                        if hasattr(val, 'item'):
                            val = val.item()
                        elif hasattr(val, 'tolist'):
                            val = val.tolist()
                        scores[key] = val
            except Exception:
                continue

    return scores


def main():
    # Parse args
    design_filter = None
    resume = False
    for arg in sys.argv[1:]:
        if arg == "--resume":
            resume = True
        elif arg.startswith("--designs="):
            design_filter = arg.split("=")[1]
        elif arg in ("short", "medium", "long"):
            design_filter = arg

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load target region sequences
    print("Loading target sequences...")
    erap2_region = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap1_region = get_region_sequence(ERAP1_PDB, *TARGET_REGION)
    print(f"  ERAP2 region ({TARGET_REGION[0]}-{TARGET_REGION[1]}): {len(erap2_region)} aa")
    print(f"  ERAP1 region ({TARGET_REGION[0]}-{TARGET_REGION[1]}): {len(erap1_region)} aa")

    # Find all V2 MPNN FASTA files
    fa_pattern = os.path.join(MPNN_DIR, "erap2v2_*.fa")
    fa_files = sorted(glob.glob(fa_pattern))
    if not fa_files:
        print(f"ERROR: No FASTA files found at {fa_pattern}")
        print("Upload MPNN results first: vast_launch.py upload")
        sys.exit(1)

    # Filter by tier if requested
    if design_filter:
        fa_files = [f for f in fa_files if design_filter in os.path.basename(f)]

    print(f"\nFound {len(fa_files)} V2 designs to validate")

    # Load existing results for resume
    summary_path = os.path.join(OUTPUT_DIR, "validation_summary_v2.json")
    existing_results = {}
    if resume and os.path.exists(summary_path):
        with open(summary_path) as f:
            for r in json.load(f):
                existing_results[r["design"]] = r
        print(f"Resume mode: {len(existing_results)} already done")

    validation_results = list(existing_results.values())
    total = len(fa_files)
    t_start = time.time()

    for idx, fa_path in enumerate(fa_files):
        design_name = os.path.splitext(os.path.basename(fa_path))[0]

        if design_name in existing_results:
            print(f"[{idx+1}/{total}] {design_name} — skipped (already done)")
            continue

        # Parse MPNN sequences, take the best one
        designed_seqs = parse_mpnn_fasta(fa_path)
        if not designed_seqs:
            print(f"[{idx+1}/{total}] {design_name} — no designed sequences, skipping")
            continue

        best = designed_seqs[0]
        binder_seq = best["sequence"]
        mpnn_score = best["score"]

        print(f"\n[{idx+1}/{total}] {design_name} (MPNN={mpnn_score:.3f}, {len(binder_seq)} aa)")

        result = {
            "design": design_name,
            "binder_sequence": binder_seq,
            "binder_length": len(binder_seq),
            "mpnn_score": mpnn_score,
        }

        # --- ERAP2 complex ---
        e2_dir = os.path.join(OUTPUT_DIR, f"erap2_{design_name}")
        e2_yaml = os.path.join(OUTPUT_DIR, f"erap2_{design_name}.yaml")
        write_boltz_yaml(erap2_region, binder_seq, e2_yaml)

        print(f"  ERAP2 complex...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e2_yaml, e2_dir)
        elapsed = time.time() - t0

        if r.returncode == 0:
            e2_scores = parse_boltz_scores(e2_dir)
            result["erap2_scores"] = e2_scores
            iptm = e2_scores.get("iptm", "?")
            print(f"OK ({elapsed:.0f}s) iptm={iptm}")
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["erap2_scores"] = {}
            result["erap2_error"] = (r.stderr or "")[-200:]

        # --- ERAP1 complex ---
        e1_dir = os.path.join(OUTPUT_DIR, f"erap1_{design_name}")
        e1_yaml = os.path.join(OUTPUT_DIR, f"erap1_{design_name}.yaml")
        write_boltz_yaml(erap1_region, binder_seq, e1_yaml)

        print(f"  ERAP1 complex...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e1_yaml, e1_dir)
        elapsed = time.time() - t0

        if r.returncode == 0:
            e1_scores = parse_boltz_scores(e1_dir)
            result["erap1_scores"] = e1_scores
            iptm = e1_scores.get("iptm", "?")
            print(f"OK ({elapsed:.0f}s) iptm={iptm}")
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["erap1_scores"] = {}
            result["erap1_error"] = (r.stderr or "")[-200:]

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

        # Save incrementally (in case of crash)
        validation_results.sort(
            key=lambda x: x.get("erap2_scores", {}).get("iptm") or 0,
            reverse=True,
        )
        with open(summary_path, "w") as f:
            json.dump(validation_results, f, indent=2)

    # Final save + report
    validation_results.sort(
        key=lambda x: x.get("erap2_scores", {}).get("iptm") or 0,
        reverse=True,
    )
    with open(summary_path, "w") as f:
        json.dump(validation_results, f, indent=2)

    elapsed_total = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"BOLTZ-2 BATCH VALIDATION COMPLETE — {len(validation_results)} designs")
    print(f"Total time: {elapsed_total/60:.1f} min")
    print("=" * 80)

    print(f"\n{'Design':<28} {'MPNN':>6} {'E2_iPTM':>8} {'E1_iPTM':>8} {'Delta':>8} {'Tag'}")
    print("-" * 80)

    for r in validation_results:
        e2i = r.get("erap2_scores", {}).get("iptm")
        e1i = r.get("erap1_scores", {}).get("iptm")
        delta = r.get("iptm_delta")
        mpnn = r.get("mpnn_score")

        tag = ""
        if delta is not None:
            tag = "SELECTIVE" if delta > 0.05 else ("MODERATE" if delta > 0 else "NON-SEL")

        print(f"{r['design']:<28} {mpnn or 0:>6.3f} "
              f"{e2i or 0:>8.4f} {e1i or 0:>8.4f} "
              f"{delta or 0:>+8.4f} {tag}")

    selective = [r for r in validation_results if (r.get("iptm_delta") or 0) > 0.05]
    print(f"\nSelective (delta > 0.05): {len(selective)}/{len(validation_results)}")
    print(f"\nOutput: {summary_path}")
    print("Run candidate scorer: python scoring/candidate_scorer.py "
          f"{summary_path}")


if __name__ == "__main__":
    main()
