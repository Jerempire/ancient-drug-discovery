"""
y87a_y89a_validation.py — Validate Y87A/Y89A IRAP-reduction mutant + multi-sample re-validation.

Two goals:
  1. Test n248_trim_c5_Y87A_Y89A against full 4-target panel (ERAP2, ERAP1, IRAP, ANPEP)
  2. Re-validate original 3 constructs (trim_c5, Y4A, wt) with 3 diffusion samples for variance

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/y87a_y89a_validation.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "y87a_validation")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
CIF_DIR = os.path.join(WORKSPACE, "results", "y87a_cif_files")

DIFFUSION_SAMPLES = 3
SEED = 42

# ======================================================================
# Sequences
# ======================================================================

# Original n248_trim_c5 (92aa): positions 87=Y, 89=Y
TRIM_C5 = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN"

# Y87A/Y89A double mutant — targeting IRAP-only interface residues
#   Position 87: Y -> A
#   Position 89: Y -> A
TRIM_C5_Y87A_Y89A = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

# Y4A variant (existing)
TRIM_C5_Y4A = "DIRHAFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN"

# Wildtype parent (97aa)
N248_WT = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK"

CONSTRUCTS = [
    {"name": "n248_trim_c5_Y87A_Y89A", "sequence": TRIM_C5_Y87A_Y89A, "new_mutant": True},
    {"name": "n248_trim_c5", "sequence": TRIM_C5, "new_mutant": False},
    {"name": "n248_trim_c5_Y4A", "sequence": TRIM_C5_Y4A, "new_mutant": False},
    {"name": "n248_wt", "sequence": N248_WT, "new_mutant": False},
]

# Targets — same as expanded_counterscreen.py
TARGETS = {
    "erap2": {
        "pdb": os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb"),
        "region": (350, 500),
    },
    "erap1": {
        "pdb": os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb"),
        "region": (350, 500),
    },
    "irap": {
        "uniprot": "Q9UIQ6",
        "region": (350, 550),
    },
    "anpep": {
        "uniprot": "P15144",
        "region": (350, 550),
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


def fetch_alphafold_sequence(uniprot_id, start, end):
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = urllib.request.urlopen(url, timeout=30)
        lines = resp.read().decode().strip().split("\n")
        seq = "".join(l for l in lines if not l.startswith(">"))
        return seq[start-1:end]
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


def find_boltz_cmd():
    for cmd in ["/opt/conda/bin/boltz", "boltz"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    print("ERROR: boltz not found")
    sys.exit(1)


def run_boltz(boltz_cmd, yaml_path, output_dir, diffusion_samples=1):
    env = os.environ.copy()
    return subprocess.run(
        [boltz_cmd, "predict", yaml_path,
         "--out_dir", output_dir,
         "--recycling_steps", "3",
         "--diffusion_samples", str(diffusion_samples),
         "--accelerator", "gpu",
         "--devices", "1"],
        capture_output=True, text=True, timeout=1200, env=env
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


def collect_cif_files(pred_dir, construct_name, target_name):
    """Copy CIF files to a central directory for PyRosetta analysis."""
    os.makedirs(CIF_DIR, exist_ok=True)
    for cif in glob.glob(os.path.join(pred_dir, "**", "*.cif"), recursive=True):
        dest = os.path.join(CIF_DIR, f"{construct_name}_{target_name}.cif")
        shutil.copy2(cif, dest)
        return dest
    return None


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(CIF_DIR, exist_ok=True)
    t_start = time.time()

    boltz_cmd = find_boltz_cmd()

    # Verify Y87A/Y89A mutations
    print("Sequence verification:")
    print(f"  trim_c5:          ...{TRIM_C5[84:92]}")
    print(f"  Y87A/Y89A mutant: ...{TRIM_C5_Y87A_Y89A[84:92]}")
    assert TRIM_C5[86] == "Y" and TRIM_C5_Y87A_Y89A[86] == "A", "Y87A mutation failed"
    assert TRIM_C5[88] == "Y" and TRIM_C5_Y87A_Y89A[88] == "A", "Y89A mutation failed"
    assert TRIM_C5[:86] == TRIM_C5_Y87A_Y89A[:86], "Prefix mismatch"
    print("  Mutations verified: Y87->A, Y89->A\n")

    # Load target sequences
    print("=" * 70)
    print("MULTI-SAMPLE VALIDATION + Y87A/Y89A IRAP-REDUCTION TEST")
    print(f"Diffusion samples: {DIFFUSION_SAMPLES}, Seed: {SEED}")
    print("=" * 70)

    target_seqs = {}
    for name, info in TARGETS.items():
        if "pdb" in info and os.path.exists(info["pdb"]):
            seq = get_region_sequence(info["pdb"], *info["region"])
            target_seqs[name] = seq
            print(f"  {name}: {len(seq)} aa (from PDB)")
        elif "uniprot" in info:
            seq = fetch_alphafold_sequence(info["uniprot"], *info["region"])
            if seq:
                target_seqs[name] = seq
                print(f"  {name}: {len(seq)} aa (from UniProt)")
            else:
                print(f"  {name}: FAILED")

    # Run predictions
    all_results = []

    for construct in CONSTRUCTS:
        cname = construct["name"]
        cseq = construct["sequence"]
        is_new = construct["new_mutant"]

        print(f"\n{'='*70}")
        label = "NEW MUTANT" if is_new else "RE-VALIDATION"
        print(f"[{label}] {cname} ({len(cseq)} aa)")
        print(f"{'='*70}")

        result = {
            "name": cname,
            "sequence": cseq,
            "length": len(cseq),
            "is_new_mutant": is_new,
            "diffusion_samples": DIFFUSION_SAMPLES,
            "screens": {},
        }

        for tname, tseq in target_seqs.items():
            print(f"\n  vs {tname}...", end=" ", flush=True)

            yaml_path = os.path.join(RESULTS_DIR, f"{cname}_{tname}.yaml")
            pred_dir = os.path.join(RESULTS_DIR, f"{cname}_{tname}")
            write_boltz_yaml(tseq, cseq, yaml_path)

            t0 = time.time()
            r = run_boltz(boltz_cmd, yaml_path, pred_dir, DIFFUSION_SAMPLES)
            elapsed = time.time() - t0

            if r.returncode == 0:
                scores = parse_boltz_scores(pred_dir)
                iptm = scores.get("iptm", 0)
                ptm = scores.get("ptm", 0)
                plddt = scores.get("complex_plddt", 0)
                print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}, ptm={ptm:.4f}")

                cif_path = collect_cif_files(pred_dir, cname, tname)
                result["screens"][tname] = {
                    "iptm": iptm, "ptm": ptm, "plddt": plddt,
                    "cif_saved": cif_path is not None,
                }
            else:
                print(f"FAILED ({elapsed:.0f}s)")
                if r.stderr:
                    print(f"    stderr: {r.stderr[:200]}")
                result["screens"][tname] = {"iptm": 0, "error": True}

        all_results.append(result)

    # ======================================================================
    # Summary
    # ======================================================================
    elapsed_total = time.time() - t_start
    print(f"\n\n{'='*100}")
    print(f"RESULTS SUMMARY ({elapsed_total/60:.1f} min, {DIFFUSION_SAMPLES} diffusion samples)")
    print(f"{'='*100}")

    targets = list(target_seqs.keys())
    header = f"{'Construct':<30}"
    for t in targets:
        header += f" {t:>10}"
    header += f" {'E2-E1':>10} {'E2-IRAP':>10} {'E2-ANPEP':>10}"
    print(header)
    print("-" * len(header))

    for r in all_results:
        tag = " *NEW*" if r["is_new_mutant"] else ""
        row = f"{r['name']:<30}"
        iptms = {}
        for t in targets:
            iptm = r["screens"].get(t, {}).get("iptm", 0)
            iptms[t] = iptm
            row += f" {iptm:>10.4f}"

        e2 = iptms.get("erap2", 0)
        for compare in ["erap1", "irap", "anpep"]:
            if compare in iptms:
                delta = e2 - iptms[compare]
                row += f" {delta:>+10.4f}"
            else:
                row += f" {'N/A':>10}"
        print(row + tag)

    # Y87A/Y89A vs original comparison
    print(f"\n{'='*70}")
    print("Y87A/Y89A IMPACT ASSESSMENT")
    print(f"{'='*70}")

    orig = next((r for r in all_results if r["name"] == "n248_trim_c5"), None)
    mutant = next((r for r in all_results if r["name"] == "n248_trim_c5_Y87A_Y89A"), None)

    if orig and mutant:
        for t in targets:
            o_iptm = orig["screens"].get(t, {}).get("iptm", 0)
            m_iptm = mutant["screens"].get(t, {}).get("iptm", 0)
            delta = m_iptm - o_iptm
            direction = "UP" if delta > 0.02 else "DOWN" if delta < -0.02 else "SAME"
            flag = ""
            if t == "irap" and delta < -0.05:
                flag = " << SUCCESS: IRAP reduced"
            elif t == "erap2" and delta < -0.10:
                flag = " << WARNING: ERAP2 binding lost"
            print(f"  {t:>8}: {o_iptm:.4f} -> {m_iptm:.4f} ({delta:+.4f}) [{direction}]{flag}")

    # Multi-sample variance (compare to original single-sample values)
    print(f"\n{'='*70}")
    print("MULTI-SAMPLE vs ORIGINAL SINGLE-SAMPLE")
    print(f"{'='*70}")

    original_scores = {
        "n248_trim_c5": {"erap2": 0.772, "erap1": 0.069, "irap": 0.576, "anpep": 0.215},
        "n248_trim_c5_Y4A": {"erap2": 0.838, "erap1": 0.195, "irap": 0.506, "anpep": 0.188},
        "n248_wt": {"erap2": 0.640, "erap1": 0.361, "irap": 0.176, "anpep": 0.238},
    }

    for r in all_results:
        if r["name"] in original_scores:
            print(f"\n  {r['name']}:")
            for t in targets:
                if t in original_scores[r["name"]]:
                    orig_val = original_scores[r["name"]][t]
                    new_val = r["screens"].get(t, {}).get("iptm", 0)
                    shift = new_val - orig_val
                    flag = " << UNSTABLE" if abs(shift) > 0.15 else ""
                    print(f"    {t:>8}: {orig_val:.3f} -> {new_val:.3f} ({shift:+.3f}){flag}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "y87a_validation_results.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved: {summary_path}")
    print(f"CIF files saved: {CIF_DIR}/")
    print(f"\nTotal time: {elapsed_total/60:.1f} min")


if __name__ == "__main__":
    main()
