"""
af2m_crossvalidation.py — AF2-Multimer cross-validation of lead constructs.

Runs ColabFold (LocalColabFold) in single-sequence mode to get independent
ipTM predictions for ERAP2 binding. Compares to Boltz-2 scores.

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/af2m_crossvalidation.py
"""
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "af2m_crossval")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")

# Lead constructs (3-sample validated scores for reference)
CONSTRUCTS = {
    "n248_trim_c5_Y87A_Y89A": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
        "boltz2_erap2": 0.748,
        "boltz2_erap1": 0.112,
        "boltz2_irap": 0.186,
    },
    "n248_trim_c5": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "boltz2_erap2": 0.691,
        "boltz2_erap1": 0.362,
        "boltz2_irap": 0.601,
    },
    "n248_wt": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "boltz2_erap2": 0.475,
        "boltz2_erap1": 0.182,
        "boltz2_irap": 0.491,
    },
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


def write_colabfold_input(target_seq, binder_seq, name, output_dir):
    """Write A3M-format input for ColabFold multimer prediction."""
    fasta_path = os.path.join(output_dir, f"{name}.fasta")
    # ColabFold multimer: separate chains with ':'
    combined = f"{target_seq}:{binder_seq}"
    with open(fasta_path, "w") as f:
        f.write(f">{name}\n{combined}\n")
    return fasta_path


def parse_af2m_scores(result_dir, name):
    """Parse ColabFold output for ipTM and pTM."""
    scores = {}
    # ColabFold saves scores in JSON
    for jf in sorted(os.listdir(result_dir)):
        if jf.endswith("_scores_rank_001.json") or jf.endswith("scores.json"):
            path = os.path.join(result_dir, jf)
            with open(path) as f:
                data = json.load(f)
            scores["iptm"] = data.get("iptm", data.get("ipTM", 0))
            scores["ptm"] = data.get("ptm", data.get("pTM", 0))
            scores["plddt"] = data.get("mean_plddt", data.get("pLDDT", 0))
            break

    # Also try ranking_debug.json
    ranking_path = os.path.join(result_dir, "ranking_debug.json")
    if os.path.exists(ranking_path) and not scores:
        with open(ranking_path) as f:
            data = json.load(f)
        if "iptm+ptm" in data:
            best_model = list(data["iptm+ptm"].keys())[0] if data["iptm+ptm"] else None
            if best_model:
                scores["iptm_ptm"] = data["iptm+ptm"][best_model]

    return scores


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    # Check ColabFold installation
    colabfold_cmd = None
    for cmd in ["colabfold_batch", "/usr/local/bin/colabfold_batch",
                "/opt/conda/bin/colabfold_batch"]:
        try:
            r = subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            if r.returncode == 0:
                colabfold_cmd = cmd
                break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue

    if not colabfold_cmd:
        print("ERROR: colabfold_batch not found. Install with:")
        print("  pip install colabfold[alphafold]")
        print("  OR: pip install -q 'colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold'")
        sys.exit(1)

    print("=" * 70)
    print("AF2-MULTIMER CROSS-VALIDATION")
    print(f"ColabFold: {colabfold_cmd}")
    print("=" * 70)

    # Load target sequences
    targets = {}

    # ERAP2 (primary)
    erap2_pdb = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
    if os.path.exists(erap2_pdb):
        targets["erap2"] = get_region_sequence(erap2_pdb, 350, 500)
        print(f"  erap2: {len(targets['erap2'])} aa (from PDB)")

    # ERAP1
    erap1_pdb = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")
    if os.path.exists(erap1_pdb):
        targets["erap1"] = get_region_sequence(erap1_pdb, 350, 500)
        print(f"  erap1: {len(targets['erap1'])} aa (from PDB)")

    # IRAP
    try:
        targets["irap"] = fetch_uniprot_sequence("Q9UIQ6", 350, 550)
        print(f"  irap: {len(targets['irap'])} aa (from UniProt)")
    except Exception as e:
        print(f"  irap: FAILED ({e})")

    all_results = []

    for cname, cinfo in CONSTRUCTS.items():
        print(f"\n{'='*60}")
        print(f"{cname} ({len(cinfo['sequence'])} aa)")
        print(f"{'='*60}")

        result = {
            "name": cname,
            "length": len(cinfo["sequence"]),
            "boltz2_scores": {k: v for k, v in cinfo.items() if k.startswith("boltz2_")},
            "af2m_scores": {},
        }

        for tname, tseq in targets.items():
            pred_name = f"{cname}_{tname}"
            input_dir = os.path.join(RESULTS_DIR, "inputs")
            output_dir = os.path.join(RESULTS_DIR, pred_name)
            os.makedirs(input_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            fasta = write_colabfold_input(tseq, cinfo["sequence"], pred_name, input_dir)

            print(f"\n  vs {tname}...", end=" ", flush=True)
            t0 = time.time()

            try:
                r = subprocess.run(
                    [colabfold_cmd, fasta, output_dir,
                     "--msa-mode", "single_sequence",
                     "--num-recycle", "3",
                     "--num-models", "1",
                     "--model-type", "alphafold2_multimer_v3",
                     "--rank", "iptm"],
                    capture_output=True, text=True, timeout=600,
                )
                elapsed = time.time() - t0

                if r.returncode == 0:
                    scores = parse_af2m_scores(output_dir, pred_name)
                    iptm = scores.get("iptm", 0)
                    ptm = scores.get("ptm", 0)
                    print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}, pTM={ptm:.4f}")
                    result["af2m_scores"][tname] = scores

                    # Compare to Boltz-2
                    b2_key = f"boltz2_{tname}"
                    if b2_key in cinfo:
                        b2_val = cinfo[b2_key]
                        delta = iptm - b2_val
                        agree = "AGREE" if abs(delta) < 0.15 else "DISAGREE"
                        print(f"    Boltz-2: {b2_val:.3f}, AF2M: {iptm:.3f}, delta: {delta:+.3f} [{agree}]")
                else:
                    print(f"FAILED ({time.time()-t0:.0f}s)")
                    if r.stderr:
                        print(f"    {r.stderr[:200]}")
                    result["af2m_scores"][tname] = {"error": True}

            except subprocess.TimeoutExpired:
                print(f"TIMEOUT (>600s)")
                result["af2m_scores"][tname] = {"error": "timeout"}

        all_results.append(result)

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n\n{'='*90}")
    print(f"CROSS-VALIDATION SUMMARY ({elapsed_total/60:.1f} min)")
    print(f"{'='*90}")
    print(f"\n{'Construct':<30} {'Target':<8} {'Boltz-2':>8} {'AF2M':>8} {'Delta':>8} {'Verdict':>10}")
    print("-" * 80)

    for r in all_results:
        for tname in targets:
            b2_key = f"boltz2_{tname}"
            b2_val = r["boltz2_scores"].get(b2_key, 0)
            af2m_val = r["af2m_scores"].get(tname, {}).get("iptm", 0)
            if af2m_val:
                delta = af2m_val - b2_val
                agree = "AGREE" if abs(delta) < 0.15 else "DISAGREE"
                print(f"{r['name']:<30} {tname:<8} {b2_val:>8.3f} {af2m_val:>8.3f} {delta:>+8.3f} {agree:>10}")
            else:
                print(f"{r['name']:<30} {tname:<8} {b2_val:>8.3f} {'N/A':>8} {'N/A':>8} {'N/A':>10}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "af2m_crossval_results.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved: {summary_path}")
    print(f"Total time: {elapsed_total/60:.1f} min")


if __name__ == "__main__":
    main()
