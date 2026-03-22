"""
v3_variant_counterscreen.py — Screen V3 candidates against K392 vs N392 ERAP2 variants.

For each candidate:
  1. Boltz-2 vs ERAP2-K392 (wildtype)
  2. Boltz-2 vs ERAP2-N392 (plague variant)
  3. Boltz-2 vs ERAP1 (paralog counter-screen)
  4. Compute variant delta (K392 - N392)

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/v3_variant_counterscreen.py
"""
import csv
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v3_counterscreen")
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


def make_k392n_variant(seq_k392, region_start, mutation_pos=392):
    """Create N392 variant by mutating K→N at position 392 within the region."""
    idx = mutation_pos - region_start
    if 0 <= idx < len(seq_k392):
        seq_list = list(seq_k392)
        original = seq_list[idx]
        seq_list[idx] = "N"
        return "".join(seq_list), original
    return seq_k392, "?"


def extract_binder_sequence(pdb_path):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    chain_seqs = {}
    for chain in struct.get_chains():
        seq = []
        for r in chain.get_residues():
            if r.get_resname() in AA3TO1:
                seq.append(AA3TO1[r.get_resname()])
        chain_seqs[chain.id] = "".join(seq)
    if len(chain_seqs) < 2:
        return list(chain_seqs.values())[0] if chain_seqs else ""
    sorted_chains = sorted(chain_seqs.items(), key=lambda x: len(x[1]))
    return sorted_chains[0][1]


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
                for key in ["ptm", "iptm", "complex_plddt"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue
    return scores


def find_top_v3_candidates(n=15):
    """Find top V3 candidates from Complexa reward CSVs."""
    candidates = []
    base = os.path.join(WORKSPACE, "Proteina-Complexa", "inference")

    for d in sorted(glob.glob(os.path.join(base, "*v3*"))):
        csvs = glob.glob(os.path.join(d, "rewards_*.csv")) or glob.glob(os.path.join(d, "all_rewards_*.csv"))
        if not csvs:
            continue
        with open(csvs[0]) as f:
            for r in csv.DictReader(f):
                candidates.append({
                    "pdb_path": r["pdb_path"],
                    "af2_iptm": float(r["af2folding_i_ptm_log"]),
                    "af2_reward": float(r["total_reward"]),
                })

    candidates.sort(key=lambda x: x["af2_reward"], reverse=True)
    return candidates[:n]


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    # Get ERAP2 region and create both variants
    print("Loading target sequences...")
    erap2_k392 = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap2_n392, original_aa = make_k392n_variant(erap2_k392, TARGET_REGION[0])
    erap1_region = get_region_sequence(ERAP1_PDB, *TARGET_REGION)

    k392_idx = 392 - TARGET_REGION[0]
    print(f"  ERAP2 K392 (wildtype): {len(erap2_k392)} aa, pos {k392_idx} = {erap2_k392[k392_idx]}")
    print(f"  ERAP2 N392 (variant):  {len(erap2_n392)} aa, pos {k392_idx} = {erap2_n392[k392_idx]}")
    print(f"  ERAP1 region: {len(erap1_region)} aa")

    # Also fetch IRAP/ANPEP
    target_seqs = {
        "erap2_k392": erap2_k392,
        "erap2_n392": erap2_n392,
        "erap1": erap1_region,
    }
    try:
        import urllib.request
        for uid, name in [("Q9UIQ6", "irap"), ("P15144", "anpep")]:
            resp = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            lines = resp.read().decode().strip().split("\n")
            seq = "".join(l for l in lines if not l.startswith(">"))
            target_seqs[name] = seq[349:550]
            print(f"  {name}: {len(target_seqs[name])} aa")
    except Exception as e:
        print(f"  WARNING: Could not fetch IRAP/ANPEP: {e}")

    # Find top V3 candidates
    top_candidates = find_top_v3_candidates(15)
    if not top_candidates:
        # If no Complexa V3 output yet, use PDB files directly
        pdb_dir = os.path.join(WORKSPACE, "v3_pdbs")
        if os.path.isdir(pdb_dir):
            for pdb in sorted(glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True))[:15]:
                seq = extract_binder_sequence(pdb)
                if seq and len(seq) > 20:
                    top_candidates.append({"pdb_path": pdb, "sequence": seq, "af2_iptm": 0, "af2_reward": 0})

    print(f"\nScreening {len(top_candidates)} V3 candidates")

    # Extract sequences from PDBs
    for c in top_candidates:
        if "sequence" not in c:
            c["sequence"] = extract_binder_sequence(c["pdb_path"])
            c["name"] = os.path.splitext(os.path.basename(c["pdb_path"]))[0]

    # Screen each candidate
    results = []
    for i, c in enumerate(top_candidates):
        name = c.get("name", f"v3_{i}")
        seq = c["sequence"]
        if not seq or len(seq) < 20:
            continue

        print(f"\n[{i+1}/{len(top_candidates)}] {name} ({len(seq)} aa)")
        result = {"name": name, "sequence": seq, "length": len(seq),
                  "af2_iptm": c.get("af2_iptm", 0), "screens": {}}

        for target_name, target_seq in target_seqs.items():
            yaml_path = os.path.join(RESULTS_DIR, f"{name}_{target_name}.yaml")
            pred_dir = os.path.join(RESULTS_DIR, f"{name}_{target_name}")
            write_boltz_yaml(target_seq, seq, yaml_path)

            print(f"  vs {target_name}...", end=" ", flush=True)
            t0 = time.time()
            r = run_boltz(yaml_path, pred_dir)
            elapsed = time.time() - t0
            if r.returncode == 0:
                scores = parse_boltz_scores(pred_dir)
                iptm = scores.get("iptm", 0)
                print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}")
                result["screens"][target_name] = iptm
            else:
                print(f"FAILED ({elapsed:.0f}s)")
                result["screens"][target_name] = 0

        # Compute variant delta
        k = result["screens"].get("erap2_k392", 0)
        n = result["screens"].get("erap2_n392", 0)
        result["variant_delta"] = round(k - n, 4)
        result["best_variant_iptm"] = max(k, n)
        result["prefers"] = "K392" if k > n else "N392" if n > k else "equal"

        results.append(result)

    # Sort by variant delta magnitude
    results.sort(key=lambda x: abs(x["variant_delta"]), reverse=True)

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n{'='*100}")
    print(f"V3 VARIANT-SELECTIVE SCREEN ({elapsed_total/60:.1f} min)")
    print(f"{'='*100}")
    print(f"{'Name':<35} {'Len':>4} {'K392':>7} {'N392':>7} {'Delta':>7} {'Pref':>6} {'E1':>7} {'IRAP':>7} {'ANPEP':>7}")
    print("-" * 100)

    for r in results:
        s = r["screens"]
        print(f"{r['name']:<35} {r['length']:>4} "
              f"{s.get('erap2_k392',0):>7.4f} {s.get('erap2_n392',0):>7.4f} "
              f"{r['variant_delta']:>+7.4f} {r['prefers']:>6} "
              f"{s.get('erap1',0):>7.4f} {s.get('irap',0):>7.4f} {s.get('anpep',0):>7.4f}")

    # Highlight interesting candidates
    interesting = [r for r in results if abs(r["variant_delta"]) >= 0.10 and r["best_variant_iptm"] >= 0.50]
    print(f"\nInteresting (|delta| >= 0.10 AND best ipTM >= 0.50): {len(interesting)}")
    for r in interesting:
        print(f"  {r['name']}: delta={r['variant_delta']:+.4f}, prefers {r['prefers']}, "
              f"K={r['screens'].get('erap2_k392',0):.3f}, N={r['screens'].get('erap2_n392',0):.3f}")

    # Save
    with open(os.path.join(RESULTS_DIR, "v3_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {os.path.join(RESULTS_DIR, 'v3_results.json')}")


if __name__ == "__main__":
    main()
