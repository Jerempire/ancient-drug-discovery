"""
expanded_counterscreen.py — Screen top Complexa candidates against full selectivity panel.

Targets: ERAP2 (primary), ERAP1, IRAP/LNPEP, Aminopeptidase N (ANPEP/CD13)
Uses AlphaFold DB structures for IRAP and ANPEP.

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/expanded_counterscreen.py
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "expanded_counterscreen")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")

# Target region: substrate channel / active site area
# Each target has its own numbering — we use equivalent functional regions
TARGETS = {
    "erap2": {
        "pdb": os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb"),
        "region": (350, 500),
        "description": "Primary target — divergent substrate channel",
    },
    "erap1": {
        "pdb": os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb"),
        "region": (350, 500),
        "description": "Closest paralog — must NOT bind",
    },
    "irap": {
        "sequence": None,  # Will be fetched from AlphaFold DB
        "uniprot": "Q9UIQ6",  # Human IRAP/LNPEP
        "region": (350, 550),  # Aminopeptidase domain
        "description": "IRAP/LNPEP — related aminopeptidase, safety check",
    },
    "anpep": {
        "sequence": None,
        "uniprot": "P15144",  # Human Aminopeptidase N / CD13
        "region": (350, 550),  # Catalytic domain
        "description": "Aminopeptidase N/CD13 — broad aminopeptidase, safety check",
    },
}

# Top 3 Complexa candidates
CANDIDATES = [
    {
        "name": "complexa_lead_n240",
        "sequence": "IYTAKNASEAIELANQIELKLKPGDLAVIVVKDNKSFPVGGMEFPGVLVLSNNFPTVPIPKTVSTDETKVIIEKLIKYFEEYFELPYPF",
        "length": 89,
        "af2_iptm": 0.846,
        "boltz_erap2_iptm": 0.868,
        "boltz_erap1_iptm": 0.166,
        "boltz_delta": 0.701,
    },
    {
        "name": "complexa_runner_n248",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "length": 97,
        "af2_iptm": 0.866,
        "boltz_erap2_iptm": 0.796,
        "boltz_erap1_iptm": 0.312,
        "boltz_delta": 0.484,
    },
    {
        "name": "complexa_third_n249",
        "sequence": "NIKVYEYDSKTGSLKDITKESKNLINLTKELSNILGLTYPYPSLNLVLVPNFGAGAMENRGLVYLNESLKNAKKIDPFSLEIAKKLFNYFKNYFELDN",
        "length": 98,
        "af2_iptm": 0.807,
        "boltz_erap2_iptm": 0.721,
        "boltz_erap1_iptm": 0.228,
        "boltz_delta": 0.493,
    },
]

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


def fetch_alphafold_sequence(uniprot_id, start, end):
    """Fetch sequence region from AlphaFold DB via UniProt."""
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = urllib.request.urlopen(url, timeout=30)
        lines = resp.read().decode().strip().split("\n")
        seq = "".join(l for l in lines if not l.startswith(">"))
        # Return region (0-indexed)
        region = seq[start-1:end]
        return region
    except Exception as e:
        print(f"  WARNING: Failed to fetch {uniprot_id}: {e}")
        return None


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
    for boltz_cmd in ["/opt/conda/bin/boltz", "boltz"]:
        try:
            subprocess.run([boltz_cmd, "--help"], capture_output=True, timeout=10)
            break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    else:
        print("ERROR: boltz not found")
        sys.exit(1)

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


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    # 1. Load/fetch target sequences
    print("=" * 70)
    print("EXPANDED SELECTIVITY PANEL")
    print("=" * 70)
    print("\nLoading target sequences...")

    target_seqs = {}
    for target_name, info in TARGETS.items():
        if "pdb" in info and info["pdb"] and os.path.exists(info["pdb"]):
            seq = get_region_sequence(info["pdb"], *info["region"])
            target_seqs[target_name] = seq
            print(f"  {target_name}: {len(seq)} aa (from PDB, region {info['region']})")
        elif "uniprot" in info:
            seq = fetch_alphafold_sequence(info["uniprot"], *info["region"])
            if seq:
                target_seqs[target_name] = seq
                print(f"  {target_name}: {len(seq)} aa (from UniProt {info['uniprot']})")
            else:
                print(f"  {target_name}: FAILED to fetch")

    if len(target_seqs) < 2:
        print("ERROR: Need at least 2 target sequences")
        return

    # 2. Screen each candidate against each target
    all_results = []

    for c in CANDIDATES:
        print(f"\n{'='*60}")
        print(f"Candidate: {c['name']} ({c['length']} aa)")
        print(f"{'='*60}")

        candidate_result = {
            "name": c["name"],
            "sequence": c["sequence"],
            "length": c["length"],
            "af2_iptm": c["af2_iptm"],
            "screens": {},
        }

        for target_name, target_seq in target_seqs.items():
            print(f"\n  vs {target_name} ({TARGETS[target_name]['description']})...", end=" ", flush=True)

            yaml_path = os.path.join(RESULTS_DIR, f"{c['name']}_{target_name}.yaml")
            pred_dir = os.path.join(RESULTS_DIR, f"{c['name']}_{target_name}")
            write_boltz_yaml(target_seq, c["sequence"], yaml_path)

            t0 = time.time()
            r = run_boltz(yaml_path, pred_dir)
            elapsed = time.time() - t0

            if r.returncode == 0:
                scores = parse_boltz_scores(pred_dir)
                iptm = scores.get("iptm", 0)
                ptm = scores.get("ptm", 0)
                plddt = scores.get("complex_plddt", 0)
                print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}")
                candidate_result["screens"][target_name] = {
                    "iptm": iptm, "ptm": ptm, "plddt": plddt,
                }
            else:
                print(f"FAILED ({elapsed:.0f}s)")
                candidate_result["screens"][target_name] = {"iptm": 0, "error": True}

        all_results.append(candidate_result)

    # 3. Summary
    elapsed_total = time.time() - t_start
    print(f"\n\n{'='*90}")
    print(f"SELECTIVITY PANEL RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*90}")

    # Header
    targets = list(target_seqs.keys())
    header = f"{'Candidate':<30}"
    for t in targets:
        header += f" {t:>10}"
    header += f" {'E2-E1':>10} {'E2-IRAP':>10} {'E2-ANPEP':>10}"
    print(header)
    print("-" * len(header))

    for r in all_results:
        row = f"{r['name']:<30}"
        iptms = {}
        for t in targets:
            iptm = r["screens"].get(t, {}).get("iptm", 0)
            iptms[t] = iptm
            row += f" {iptm:>10.4f}"

        # Deltas
        e2 = iptms.get("erap2", 0)
        for compare in ["erap1", "irap", "anpep"]:
            if compare in iptms:
                delta = e2 - iptms[compare]
                row += f" {delta:>+10.4f}"
            else:
                row += f" {'N/A':>10}"
        print(row)

    # Selectivity summary
    print(f"\nSelectivity Assessment:")
    for r in all_results:
        e2 = r["screens"].get("erap2", {}).get("iptm", 0)
        screens = r["screens"]
        issues = []
        for t in ["erap1", "irap", "anpep"]:
            if t in screens and screens[t].get("iptm", 0) > 0.3:
                issues.append(f"{t}={screens[t]['iptm']:.3f}")

        if not issues:
            verdict = "CLEAN — selective across entire panel"
        elif len(issues) == 1 and "erap1" in issues[0]:
            verdict = f"MODERATE — some ERAP1 signal ({issues[0]})"
        else:
            verdict = f"CONCERN — cross-reactivity: {', '.join(issues)}"

        print(f"  {r['name']}: {verdict}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "expanded_panel_results.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
