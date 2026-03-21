"""
trimc5_optimization.py — Expanded panel + N-core de-aromatization for n248_trim_c5.

1. Re-run ERAP2/ERAP1/IRAP/ANPEP panel on trim_c5
2. Test N-core aromatic reversals on trim_c5 backbone:
   - Y4A only (single most suspicious residue)
   - Y4A + F5A (pair)
   - Y4A + F5A + Y11A + Y23A (full N-core de-aromatization)

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/trimc5_optimization.py
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "trimc5_optimization")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")
TARGET_REGION = (350, 500)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# n248 wild-type (97aa) and trim_c5 (92aa)
N248_WT = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK"
N248_TRIM_C5 = N248_WT[:-5]  # Remove YFFEK

# Targets for expanded panel
IRAP_UNIPROT = "Q9UIQ6"
ANPEP_UNIPROT = "P15144"


def get_region_sequence(pdb_path, start, end):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    seq = []
    for chain in struct.get_chains():
        if chain.id == "A":
            for r in chain.get_residues():
                if r.get_resname() in AA3TO1 and start <= r.get_id()[1] <= end:
                    seq.append(AA3TO1[r.get_resname()])
            break
    return "".join(seq)


def fetch_uniprot_region(uniprot_id, start, end):
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp = urllib.request.urlopen(url, timeout=30)
    lines = resp.read().decode().strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    return seq[start-1:end]


def mutate(seq, positions):
    seq_list = list(seq)
    for pos in positions:
        if 0 <= pos < len(seq_list):
            seq_list[pos] = "A"
    return "".join(seq_list)


def write_boltz_yaml(target_seq, binder_seq, output_path):
    with open(output_path, "w") as f:
        f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {target_seq}
      msa: empty
  - protein:
      id: B
      sequence: {binder_seq}
      msa: empty
""")


def run_boltz(yaml_path, output_dir):
    boltz_cmd = "/opt/conda/bin/boltz"
    env = os.environ.copy()
    return subprocess.run(
        [boltz_cmd, "predict", yaml_path, "--out_dir", output_dir,
         "--recycling_steps", "3", "--diffusion_samples", "1",
         "--accelerator", "gpu", "--devices", "1"],
        capture_output=True, text=True, timeout=900, env=env
    )


def parse_boltz_scores(output_dir):
    scores = {}
    for jf in glob.glob(os.path.join(output_dir, "**", "*.json"), recursive=True):
        if "manifest" in jf:
            continue
        try:
            with open(jf) as f:
                data = json.load(f)
            if isinstance(data, dict):
                for key in ["ptm", "iptm", "complex_plddt", "pair_chains_iptm"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue
    return scores


def screen_variant(name, seq, target_seqs, results_dir):
    """Screen a variant against all targets. Returns dict of results."""
    result = {"name": name, "sequence": seq, "length": len(seq), "screens": {}}

    for target_name, target_seq in target_seqs.items():
        yaml_path = os.path.join(results_dir, f"{name}_{target_name}.yaml")
        pred_dir = os.path.join(results_dir, f"{name}_{target_name}")
        write_boltz_yaml(target_seq, seq, yaml_path)

        print(f"  vs {target_name}...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(yaml_path, pred_dir)
        elapsed = time.time() - t0

        if r.returncode == 0:
            scores = parse_boltz_scores(pred_dir)
            iptm = scores.get("iptm", 0)
            print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}")
            result["screens"][target_name] = {"iptm": iptm, "ptm": scores.get("ptm", 0)}
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["screens"][target_name] = {"iptm": 0, "error": True}

    return result


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    # Load target sequences
    print("Loading target sequences...")
    target_seqs = {
        "erap2": get_region_sequence(ERAP2_PDB, *TARGET_REGION),
        "erap1": get_region_sequence(ERAP1_PDB, *TARGET_REGION),
    }

    # Fetch IRAP and ANPEP
    try:
        target_seqs["irap"] = fetch_uniprot_region(IRAP_UNIPROT, 350, 550)
        print(f"  irap: {len(target_seqs['irap'])} aa")
    except Exception as e:
        print(f"  irap: FAILED ({e})")
    try:
        target_seqs["anpep"] = fetch_uniprot_region(ANPEP_UNIPROT, 350, 550)
        print(f"  anpep: {len(target_seqs['anpep'])} aa")
    except Exception as e:
        print(f"  anpep: FAILED ({e})")

    for name, seq in target_seqs.items():
        print(f"  {name}: {len(seq)} aa")

    # Define variants
    variants = [
        ("trim_c5", N248_TRIM_C5,
         "Primary lead (92aa) — expanded panel recheck"),
        ("trim_c5_Y4A", mutate(N248_TRIM_C5, [4]),
         "Trim + Y4A — single N-core de-aromatization"),
        ("trim_c5_Y4A_F5A", mutate(N248_TRIM_C5, [4, 5]),
         "Trim + Y4A + F5A — pair de-aromatization"),
        ("trim_c5_ncore_dearom", mutate(N248_TRIM_C5, [4, 5, 11, 23]),
         "Trim + full N-core de-aromatization (Y4A/F5A/Y11A/Y23A)"),
    ]

    print(f"\n{'='*70}")
    print(f"TRIM_C5 OPTIMIZATION — {len(variants)} variants x {len(target_seqs)} targets")
    print(f"{'='*70}")

    for name, seq, desc in variants:
        print(f"\n  {name}: {desc}")
        # Show mutations vs trim_c5
        if name != "trim_c5":
            diffs = [f"{N248_TRIM_C5[i]}{i+1}{seq[i]}" for i in range(len(seq)) if seq[i] != N248_TRIM_C5[i]]
            print(f"  Mutations: {', '.join(diffs)}")

    all_results = []
    for name, seq, desc in variants:
        print(f"\n{'='*60}")
        print(f"{name} ({len(seq)} aa) — {desc}")
        print(f"{'='*60}")
        result = screen_variant(name, seq, target_seqs, RESULTS_DIR)
        result["description"] = desc
        all_results.append(result)

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n\n{'='*100}")
    print(f"OPTIMIZATION RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*100}")

    targets = list(target_seqs.keys())
    header = f"{'Variant':<30} {'Len':>4}"
    for t in targets:
        header += f" {t:>8}"
    header += f" {'E2-E1':>8} {'E2-IRAP':>8} {'E2-ANPEP':>8} {'Clean?'}"
    print(header)
    print("-" * len(header))

    for r in all_results:
        row = f"{r['name']:<30} {r['length']:>4}"
        iptms = {}
        for t in targets:
            iptm = r["screens"].get(t, {}).get("iptm", 0)
            iptms[t] = iptm
            row += f" {iptm:>8.4f}"

        e2 = iptms.get("erap2", 0)
        deltas = {}
        for compare in ["erap1", "irap", "anpep"]:
            if compare in iptms:
                d = e2 - iptms[compare]
                deltas[compare] = d
                row += f" {d:>+8.4f}"
            else:
                row += f" {'N/A':>8}"

        # Clean = no off-target > 0.3
        concerns = [t for t in ["erap1", "irap", "anpep"] if iptms.get(t, 0) > 0.3]
        clean = "YES" if not concerns else f"NO ({','.join(concerns)})"
        row += f" {clean}"
        print(row)

    # Best variant
    print(f"\n--- RECOMMENDATION ---")
    best = max(all_results, key=lambda r: (
        r["screens"].get("erap2", {}).get("iptm", 0) -
        max(r["screens"].get("erap1", {}).get("iptm", 0),
            r["screens"].get("irap", {}).get("iptm", 0),
            r["screens"].get("anpep", {}).get("iptm", 0))
    ))
    e2 = best["screens"].get("erap2", {}).get("iptm", 0)
    worst_off = max(
        best["screens"].get("erap1", {}).get("iptm", 0),
        best["screens"].get("irap", {}).get("iptm", 0),
        best["screens"].get("anpep", {}).get("iptm", 0),
    )
    print(f"Best: {best['name']} ({best['length']} aa)")
    print(f"  ERAP2 ipTM: {e2:.4f}")
    print(f"  Worst off-target: {worst_off:.4f}")
    print(f"  Selectivity gap: {e2 - worst_off:+.4f}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "optimization_results.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
