"""
run_mpnn_selective.py — Run ProteinMPNN on top selective ERAP2 binders.

For each binder PDB:
  1. Fix chain IDs (ERAP2=A, binder=B)
  2. Run ProteinMPNN to redesign ONLY the binder chain (fix ERAP2)
  3. Generate 8 sequence candidates per design
  4. Score and rank by ProteinMPNN log-likelihood

Outputs:
  /workspace/results/mpnn/seqs/<design>/*.fa      — designed sequences
  /workspace/results/mpnn/scores/<design>.json     — per-sequence scores
  /workspace/results/mpnn/mpnn_summary.json        — ranked summary
"""
import glob
import json
import os
import subprocess
import sys

from Bio.PDB import PDBParser, PDBIO, Select

MPNN_DIR = "/workspace/ProteinMPNN"
BINDER_DIR = "/workspace/results/rfdiffusion"
OUTPUT_DIR = "/workspace/results/mpnn"
NUM_SEQS = 8  # sequences per design
TEMPERATURE = 0.1  # lower = more confident sequences

# Top designs from selectivity v3 (channel campaign)
TOP_DESIGNS = [
    "erap2_sel_channel_13.pdb",
    "erap2_sel_channel_8.pdb",
    "erap2_sel_channel_2.pdb",
    "erap2_sel_channel_4.pdb",
    "erap2_sel_channel_3.pdb",
    "erap2_sel_channel_1.pdb",
    "erap2_sel_channel_0.pdb",
    "erap2_sel_channel_14.pdb",
    "erap2_sel_channel_7.pdb",
    "erap2_sel_channel_19.pdb",
]

parser = PDBParser(QUIET=True)


class ChainSelect(Select):
    """Select specific chains."""
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def prepare_pdb(pdb_path, out_path):
    """Ensure PDB has chain A (target) and chain B (binder)."""
    struct = parser.get_structure("s", pdb_path)
    chains = list(struct.get_chains())

    if len(chains) < 2:
        print("  WARNING: %s has only %d chain(s), skipping" % (pdb_path, len(chains)))
        return False

    # Write both chains
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path)
    return True


def make_chain_json(pdb_path, design_name):
    """Create ProteinMPNN chain definition JSON.

    Fix chain A (ERAP2), redesign chain B (binder).
    """
    struct = parser.get_structure("s", pdb_path)
    chains = list(struct.get_chains())
    chain_ids = [c.id for c in chains]

    # Chain to redesign = non-A chain
    binder_chain = [c for c in chain_ids if c != "A"]
    if not binder_chain:
        return None

    # ProteinMPNN format: specify which chains to design
    chain_spec = {
        pdb_path: {
            "fixed_chains": ["A"],
            "designed_chains": binder_chain,
        }
    }

    json_path = os.path.join(OUTPUT_DIR, "chain_defs", "%s.json" % design_name)
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(chain_spec, f, indent=2)
    return json_path


def run_mpnn(pdb_path, design_name):
    """Run ProteinMPNN on a single design."""
    out_seqs = os.path.join(OUTPUT_DIR, "seqs")
    os.makedirs(out_seqs, exist_ok=True)

    # Create jsonl for input
    pdb_dir = os.path.dirname(pdb_path)
    jsonl_path = os.path.join(OUTPUT_DIR, "parsed", "%s.jsonl" % design_name)
    os.makedirs(os.path.dirname(jsonl_path), exist_ok=True)

    # Step 1: Parse PDB to jsonl
    parse_cmd = [
        "python3", os.path.join(MPNN_DIR, "helper_scripts", "parse_multiple_chains.py"),
        "--input_path", pdb_dir,
        "--output_path", jsonl_path,
    ]
    subprocess.run(parse_cmd, capture_output=True, text=True)

    # Step 2: Create chains to design JSON
    chains_json = os.path.join(OUTPUT_DIR, "chain_defs", "%s_chains.jsonl" % design_name)
    os.makedirs(os.path.dirname(chains_json), exist_ok=True)

    # Use the assign_fixed_chains helper
    assign_cmd = [
        "python3", os.path.join(MPNN_DIR, "helper_scripts", "assign_fixed_chains.py"),
        "--input_path", jsonl_path,
        "--output_path", chains_json,
        "--chain_list", "A",  # Fix chain A (ERAP2)
    ]
    subprocess.run(assign_cmd, capture_output=True, text=True)

    # Step 3: Run ProteinMPNN
    mpnn_cmd = [
        "python3", os.path.join(MPNN_DIR, "protein_mpnn_run.py"),
        "--jsonl_path", jsonl_path,
        "--chain_id_jsonl", chains_json,
        "--out_folder", out_seqs,
        "--num_seq_per_target", str(NUM_SEQS),
        "--sampling_temp", str(TEMPERATURE),
        "--batch_size", "1",
        "--seed", "42",
    ]

    result = subprocess.run(mpnn_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("  MPNN ERROR: %s" % result.stderr[-200:] if result.stderr else "unknown")
        return None

    return out_seqs


def parse_mpnn_output(design_name, seqs_dir):
    """Parse ProteinMPNN FASTA output and extract scores."""
    # Find the output fasta
    fasta_files = glob.glob(os.path.join(seqs_dir, "seqs", "*.fa"))
    if not fasta_files:
        fasta_files = glob.glob(os.path.join(seqs_dir, "**", "*.fa"), recursive=True)

    results = []
    for fasta_path in fasta_files:
        if design_name.replace(".pdb", "") not in fasta_path:
            continue
        with open(fasta_path) as f:
            lines = f.readlines()

        for i in range(0, len(lines), 2):
            if i + 1 >= len(lines):
                break
            header = lines[i].strip()
            sequence = lines[i + 1].strip()

            # Parse score from header
            # Format: >design_name, score=X.XXX, ...
            score = None
            global_score = None
            seq_recovery = None
            for part in header.split(","):
                part = part.strip()
                if "score=" in part:
                    try:
                        score = float(part.split("=")[1])
                    except (ValueError, IndexError):
                        pass
                if "global_score=" in part:
                    try:
                        global_score = float(part.split("=")[1])
                    except (ValueError, IndexError):
                        pass
                if "seq_recovery=" in part:
                    try:
                        seq_recovery = float(part.split("=")[1])
                    except (ValueError, IndexError):
                        pass

            results.append({
                "header": header[:100],
                "sequence": sequence,
                "length": len(sequence),
                "score": score,
                "global_score": global_score,
                "seq_recovery": seq_recovery,
            })

    return results


# --- Main ---
print("=" * 70)
print("ProteinMPNN Sequence Design — Top Selective ERAP2 Binders")
print("=" * 70)

# Prepare input directory with only top designs
prep_dir = os.path.join(OUTPUT_DIR, "input_pdbs")
os.makedirs(prep_dir, exist_ok=True)

for design in TOP_DESIGNS:
    src = os.path.join(BINDER_DIR, design)
    dst = os.path.join(prep_dir, design)
    if os.path.exists(src):
        prepare_pdb(src, dst)

# Run MPNN on all designs at once (more efficient)
print("\nParsing PDBs...")
jsonl_path = os.path.join(OUTPUT_DIR, "parsed", "all_selective.jsonl")
os.makedirs(os.path.dirname(jsonl_path), exist_ok=True)

parse_cmd = [
    "python3", os.path.join(MPNN_DIR, "helper_scripts", "parse_multiple_chains.py"),
    "--input_path", prep_dir,
    "--output_path", jsonl_path,
]
r = subprocess.run(parse_cmd, capture_output=True, text=True)
print("  Parse: %s" % ("OK" if r.returncode == 0 else r.stderr[-100:]))

# Assign fixed chains
chains_json = os.path.join(OUTPUT_DIR, "parsed", "fixed_chains.jsonl")
assign_cmd = [
    "python3", os.path.join(MPNN_DIR, "helper_scripts", "assign_fixed_chains.py"),
    "--input_path", jsonl_path,
    "--output_path", chains_json,
    "--chain_list", "A",
]
r = subprocess.run(assign_cmd, capture_output=True, text=True)
print("  Chains: %s" % ("OK" if r.returncode == 0 else r.stderr[-100:]))

# Run MPNN
print("\nRunning ProteinMPNN (%d designs x %d sequences = %d total)..." % (
    len(TOP_DESIGNS), NUM_SEQS, len(TOP_DESIGNS) * NUM_SEQS))

seqs_dir = os.path.join(OUTPUT_DIR, "output")
mpnn_cmd = [
    "python3", os.path.join(MPNN_DIR, "protein_mpnn_run.py"),
    "--jsonl_path", jsonl_path,
    "--chain_id_jsonl", chains_json,
    "--out_folder", seqs_dir,
    "--num_seq_per_target", str(NUM_SEQS),
    "--sampling_temp", str(TEMPERATURE),
    "--batch_size", "1",
    "--seed", "42",
]
r = subprocess.run(mpnn_cmd, capture_output=True, text=True)
print(r.stdout[-500:] if r.stdout else "")
if r.returncode != 0:
    print("MPNN STDERR: %s" % r.stderr[-500:])
    sys.exit(1)

# Parse results
print("\nParsing MPNN output...")
fasta_files = sorted(glob.glob(os.path.join(seqs_dir, "seqs", "*.fa")))
print("  Found %d FASTA files" % len(fasta_files))

all_results = []
for fasta_path in fasta_files:
    design_name = os.path.basename(fasta_path).replace(".fa", "")
    with open(fasta_path) as f:
        lines = f.readlines()

    for i in range(0, len(lines), 2):
        if i + 1 >= len(lines):
            break
        header = lines[i].strip()
        sequence = lines[i + 1].strip()

        score = None
        global_score = None
        seq_recovery = None
        for part in header.split(","):
            part = part.strip()
            if part.startswith("score="):
                try:
                    score = float(part.split("=")[1])
                except (ValueError, IndexError):
                    pass
            elif part.startswith("global_score="):
                try:
                    global_score = float(part.split("=")[1])
                except (ValueError, IndexError):
                    pass
            elif part.startswith("seq_recovery="):
                try:
                    seq_recovery = float(part.split("=")[1])
                except (ValueError, IndexError):
                    pass

        all_results.append({
            "design": design_name,
            "header": header[:120],
            "sequence": sequence,
            "length": len(sequence),
            "score": score,
            "global_score": global_score,
            "seq_recovery": seq_recovery,
        })

# Sort by score (lower = better for MPNN)
all_results.sort(key=lambda x: x.get("score") or 999)

# Save
summary_path = os.path.join(OUTPUT_DIR, "mpnn_summary.json")
with open(summary_path, "w") as f:
    json.dump(all_results, f, indent=2)

# Report
print("\n" + "=" * 70)
print("PROTEINMPNN RESULTS — ERAP2 Selective Binder Sequences")
print("=" * 70)
print("\nTotal sequences designed: %d" % len(all_results))

if all_results:
    scores = [r["score"] for r in all_results if r["score"] is not None]
    if scores:
        print("Score range: %.3f to %.3f (lower = better)" % (min(scores), max(scores)))
        print("Mean score: %.3f" % (sum(scores) / len(scores)))

    print("\nTOP 15 SEQUENCES:")
    print("%-40s %6s %4s %s" % ("Design", "Score", "Len", "Sequence (first 50)"))
    print("-" * 105)
    for r in all_results[:15]:
        s = r.get("score")
        score_str = "%.3f" % s if s is not None else "  N/A"
        print("%-40s %6s %4d %s..." % (
            r["design"][:40], score_str, r["length"], r["sequence"][:50]))

    print("\n*** BEST SEQUENCE ***")
    best = all_results[0]
    print("  Design: %s" % best["design"])
    print("  Score: %s" % best.get("score"))
    print("  Length: %d aa" % best["length"])
    print("  Sequence: %s" % best["sequence"])

print("\nSaved: %s" % summary_path)
