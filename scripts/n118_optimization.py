"""
n118_optimization.py — Optimize long_n118 backup candidate for selectivity.

Strategy: Same as n248 Y87A/Y89A — mutate C-terminal aromatics that likely
drive IRAP/ANPEP cross-reactivity. n118 has F87, Y89, Y90, Y93, F94 cluster.

Tests:
  1. n118_Y89A_Y90A — conservative double (mirrors n248 approach)
  2. n118_F87A_Y89A_Y90A — triple, remove aromatic cluster core
  3. n118_trim_c4 — remove AAAA tail (96aa), expose aromatics differently
  4. n118_trim_c10 — remove last 10aa including aromatic cluster (90aa)

Usage (on Vast.ai/RunPod GPU):
    python3 /workspace/scripts/n118_optimization.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "n118_optimization")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")

DIFFUSION_SAMPLES = 3

# Original n118 sequence (100aa)
N118_WT = "KTYFLGKDETTGIEYYLMIITDPELFKKLTGSANTLAFAILDIGKKTAYVVLGPNLLKKRPEIKNPKIKPSKLPGEEYILEVSKKLFSYYKEYFELAAAA"

# Variants
CONSTRUCTS = [
    {
        "name": "n118_wt",
        "sequence": N118_WT,
        "description": "Original (ERAP2=0.877, ERAP1=0.303, IRAP=0.351, ANPEP=0.638)",
    },
    {
        "name": "n118_Y89A_Y90A",
        "sequence": N118_WT[:88] + "AA" + N118_WT[90:],
        "description": "Double Ala at C-terminal Y89,Y90 (mirrors n248 approach)",
    },
    {
        "name": "n118_F87A_Y89A_Y90A",
        "sequence": N118_WT[:86] + "A" + N118_WT[87:88] + "AA" + N118_WT[90:],
        "description": "Triple: F87A + Y89A + Y90A",
    },
    {
        "name": "n118_5aro_cleanup",
        "sequence": N118_WT[:86] + "A" + N118_WT[87:88] + "AA" + N118_WT[90:92] + "AA" + N118_WT[94:],
        "description": "All 5 C-terminal aromatics → Ala (F87A,Y89A,Y90A,Y93A,F94A)",
    },
    {
        "name": "n118_trim_c4",
        "sequence": N118_WT[:96],
        "description": "Remove AAAA tail (96aa)",
    },
    {
        "name": "n118_trim_c10",
        "sequence": N118_WT[:90],
        "description": "Remove last 10aa including aromatic cluster (90aa)",
    },
]

TARGETS = {
    "erap2": {"pdb": os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb"), "region": (350, 500)},
    "erap1": {"pdb": os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb"), "region": (350, 500)},
    "irap":  {"uniprot": "Q9UIQ6", "region": (350, 550)},
    "anpep": {"uniprot": "P15144", "region": (350, 550)},
}

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


def fetch_uniprot_sequence(uniprot_id, start, end):
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp = urllib.request.urlopen(url, timeout=30)
    lines = resp.read().decode().strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    return seq[start-1:end]


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


def find_boltz_cmd():
    for cmd in ["/usr/local/bin/boltz", "/opt/conda/bin/boltz", "boltz"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    print("ERROR: boltz not found")
    sys.exit(1)


def run_boltz(boltz_cmd, yaml_path, output_dir):
    env = os.environ.copy()
    return subprocess.run(
        [boltz_cmd, "predict", yaml_path,
         "--out_dir", output_dir,
         "--recycling_steps", "3",
         "--diffusion_samples", str(DIFFUSION_SAMPLES),
         "--accelerator", "gpu", "--devices", "1"],
        capture_output=True, text=True, timeout=600, env=env
    )


def parse_boltz_scores(output_dir):
    for jf in glob.glob(os.path.join(output_dir, "**", "*.json"), recursive=True):
        if "manifest" in jf:
            continue
        try:
            with open(jf) as f:
                data = json.load(f)
            if isinstance(data, dict) and "iptm" in data:
                return {"iptm": data["iptm"], "ptm": data.get("ptm", 0)}
        except:
            continue
    return {"iptm": 0, "ptm": 0}


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()
    boltz_cmd = find_boltz_cmd()

    # Verify mutations
    print("Construct verification:")
    for c in CONSTRUCTS:
        print(f"  {c['name']} ({len(c['sequence'])}aa): ...{c['sequence'][-15:]}")

    # Load targets
    print("\nLoading targets...")
    target_seqs = {}
    for tname, info in TARGETS.items():
        if "pdb" in info and os.path.exists(info["pdb"]):
            target_seqs[tname] = get_region_sequence(info["pdb"], *info["region"])
        elif "uniprot" in info:
            try:
                target_seqs[tname] = fetch_uniprot_sequence(info["uniprot"], *info["region"])
            except:
                print(f"  WARNING: Failed to fetch {tname}")
        print(f"  {tname}: {len(target_seqs.get(tname, ''))} aa")

    # Run predictions
    all_results = []
    for c in CONSTRUCTS:
        print(f"\n{'='*60}")
        print(f"{c['name']} ({len(c['sequence'])}aa) — {c['description']}")
        print(f"{'='*60}")

        result = {
            "name": c["name"],
            "sequence": c["sequence"],
            "length": len(c["sequence"]),
            "description": c["description"],
            "screens": {},
        }

        for tname, tseq in target_seqs.items():
            yaml_path = os.path.join(RESULTS_DIR, f"{c['name']}_{tname}.yaml")
            pred_dir = os.path.join(RESULTS_DIR, f"{c['name']}_{tname}")
            write_boltz_yaml(tseq, c["sequence"], yaml_path)

            t0 = time.time()
            r = run_boltz(boltz_cmd, yaml_path, pred_dir)
            elapsed = time.time() - t0

            if r.returncode == 0:
                scores = parse_boltz_scores(pred_dir)
                print(f"  {tname}: ipTM={scores['iptm']:.4f} ({elapsed:.0f}s)")
                result["screens"][tname] = scores["iptm"]
            else:
                print(f"  {tname}: FAILED ({elapsed:.0f}s)")
                result["screens"][tname] = 0

        all_results.append(result)

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n\n{'='*90}")
    print(f"n118 OPTIMIZATION RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*90}")

    targets_list = list(target_seqs.keys())
    print(f"{'Construct':<25} {'Len':>4}", end="")
    for t in targets_list:
        print(f" {t:>8}", end="")
    print(f" {'Clean':>6}")
    print("-" * 80)

    for r in all_results:
        print(f"{r['name']:<25} {r['length']:>4}", end="")
        clean = True
        for t in targets_list:
            val = r["screens"].get(t, 0)
            print(f" {val:>8.3f}", end="")
            if t != "erap2" and val > 0.3:
                clean = False
        print(f" {'YES' if clean else 'no':>6}")

    # Benchmark
    print(f"\n{'Y87A/Y89A benchmark':<25} {'92':>4}   {'0.748':>8}   {'0.112':>8}   {'0.186':>8}   {'0.183':>8} {'YES':>6}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "n118_optimization_results.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved: {summary_path}")


if __name__ == "__main__":
    main()
