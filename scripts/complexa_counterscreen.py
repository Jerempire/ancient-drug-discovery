"""
complexa_counterscreen.py — Boltz-2 ERAP1 counter-screen for Proteina-Complexa candidates.

Extracts binder sequences from Complexa PDB outputs, runs Boltz-2 against
both ERAP2 and ERAP1, computes selectivity metrics.

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/complexa_counterscreen.py
"""
import glob
import json
import os
import subprocess
import sys
import time

# --- Config ---
WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "complexa_counterscreen")
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

# Top 10 Complexa candidates (sequences extracted from PDBs on instance)
CANDIDATES = []  # Will be populated from PDB files


def get_region_sequence(pdb_path, start, end):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    seq = []
    for chain in struct.get_chains():
        if chain.id == "A":
            for r in chain.get_residues():
                resnum = r.get_id()[1]
                resname = r.get_resname()
                if resname in AA3TO1 and start <= resnum <= end:
                    seq.append(AA3TO1[resname])
            break
    return "".join(seq)


def extract_binder_sequence(pdb_path):
    """Extract binder chain sequence from Complexa PDB output."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    chains = list(struct.get_chains())

    # Complexa outputs: target is chain A (longer), binder is chain B (shorter)
    # Or sometimes just two chains — binder is the shorter one
    chain_seqs = {}
    for chain in chains:
        seq = []
        for r in chain.get_residues():
            resname = r.get_resname()
            if resname in AA3TO1:
                seq.append(AA3TO1[resname])
        chain_seqs[chain.id] = "".join(seq)

    if len(chain_seqs) < 2:
        # Single chain — try to split by chain ID
        if "B" in chain_seqs:
            return chain_seqs["B"]
        return list(chain_seqs.values())[0] if chain_seqs else ""

    # Return the shorter chain as binder
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
    for boltz_cmd in ["/opt/conda/envs/boltz/bin/boltz", "/opt/conda/bin/boltz", "boltz"]:
        try:
            subprocess.run([boltz_cmd, "--help"], capture_output=True, timeout=10)
            break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    else:
        print("ERROR: boltz not found")
        sys.exit(1)

    env = os.environ.copy()
    # Add CUDA nvrtc to LD_LIBRARY_PATH
    for path in glob.glob("/opt/conda/**/libnvrtc.so*", recursive=True):
        env["LD_LIBRARY_PATH"] = os.path.dirname(path) + ":" + env.get("LD_LIBRARY_PATH", "")
        break

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
                for key in ["ptm", "iptm", "complex_plddt", "pair_chains_iptm", "confidence_score"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue
    return scores


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load target region sequences
    print("Loading target sequences...")
    erap2_region = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap1_region = get_region_sequence(ERAP1_PDB, *TARGET_REGION)
    print(f"  ERAP2 region: {len(erap2_region)} aa")
    print(f"  ERAP1 region: {len(erap1_region)} aa")

    # Find Complexa PDB files
    pdb_dir = os.path.join(WORKSPACE, "complexa_pdbs")
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True))
    print(f"\nFound {len(pdb_files)} Complexa PDB files")

    # Extract binder sequences
    candidates = []
    for pdb_path in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        seq = extract_binder_sequence(pdb_path)
        if seq and len(seq) > 20:
            candidates.append({"name": name, "sequence": seq, "length": len(seq), "pdb": pdb_path})
            print(f"  {name}: {len(seq)} aa")
        else:
            print(f"  {name}: skipped (seq too short: {len(seq)})")

    if not candidates:
        print("ERROR: No valid candidates found")
        return

    print(f"\n{'='*70}")
    print(f"COUNTER-SCREENING {len(candidates)} CANDIDATES")
    print(f"{'='*70}")

    results = []
    t_start = time.time()

    for i, c in enumerate(candidates):
        name = c["name"]
        seq = c["sequence"]
        print(f"\n[{i+1}/{len(candidates)}] {name} ({c['length']} aa)")

        result = {"design": name, "sequence": seq, "length": c["length"]}

        # ERAP2
        e2_yaml = os.path.join(RESULTS_DIR, f"e2_{name}.yaml")
        e2_dir = os.path.join(RESULTS_DIR, f"e2_{name}")
        write_boltz_yaml(erap2_region, seq, e2_yaml)
        print(f"  ERAP2...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e2_yaml, e2_dir)
        elapsed = time.time() - t0
        if r.returncode == 0:
            e2 = parse_boltz_scores(e2_dir)
            result["erap2_iptm"] = e2.get("iptm", 0)
            result["erap2_ptm"] = e2.get("ptm", 0)
            result["erap2_plddt"] = e2.get("complex_plddt", 0)
            print(f"OK ({elapsed:.0f}s) iptm={result['erap2_iptm']:.4f}")
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["erap2_iptm"] = 0

        # ERAP1
        e1_yaml = os.path.join(RESULTS_DIR, f"e1_{name}.yaml")
        e1_dir = os.path.join(RESULTS_DIR, f"e1_{name}")
        write_boltz_yaml(erap1_region, seq, e1_yaml)
        print(f"  ERAP1...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e1_yaml, e1_dir)
        elapsed = time.time() - t0
        if r.returncode == 0:
            e1 = parse_boltz_scores(e1_dir)
            result["erap1_iptm"] = e1.get("iptm", 0)
            result["erap1_ptm"] = e1.get("ptm", 0)
            print(f"OK ({elapsed:.0f}s) iptm={result['erap1_iptm']:.4f}")
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["erap1_iptm"] = 0

        # Selectivity
        e2i = result.get("erap2_iptm", 0)
        e1i = result.get("erap1_iptm", 0)
        result["iptm_delta"] = round(e2i - e1i, 4)
        result["iptm_selectivity"] = round(e2i / max(e1i, 0.001), 3)
        result["selective"] = result["iptm_delta"] > 0.05

        results.append(result)

    # Save results
    elapsed_total = time.time() - t_start
    results.sort(key=lambda x: x.get("iptm_delta", 0), reverse=True)

    summary_path = os.path.join(RESULTS_DIR, "counterscreen_results.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    # Report
    print(f"\n{'='*80}")
    print(f"COUNTER-SCREEN RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*80}")
    print(f"{'Design':<50} {'E2_ipTM':>8} {'E1_ipTM':>8} {'Delta':>8} {'Sel?'}")
    print("-" * 80)
    for r in results:
        sel = "YES" if r["selective"] else "no"
        print(f"{r['design']:<50} {r['erap2_iptm']:>8.4f} {r['erap1_iptm']:>8.4f} {r['iptm_delta']:>+8.4f} {sel}")

    selective = [r for r in results if r["selective"]]
    print(f"\nSelective (delta > 0.05): {len(selective)}/{len(results)}")
    if selective:
        best = selective[0]
        print(f"Best: {best['design']} (delta={best['iptm_delta']:+.4f}, E2={best['erap2_iptm']:.4f})")

    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
