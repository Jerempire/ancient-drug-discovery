"""
v2_k392_n392_test.py — Test V2 synthesis candidates for K392 vs N392 variant selectivity.

Runs all 4 V2 constructs against both ERAP2-K392 (wildtype) and ERAP2-N392
(plague-selected variant) with 3 diffusion samples for reliability.

The K392 data already exists (from y87a_y89a_validation.py). This script
fills in the missing N392 column and computes the variant selectivity delta.

Key question: Is the V2 lead (Y87A_Y89A, ipTM 0.748 on K392) variant-selective?

Possible outcomes:
  - Delta > +0.30: Strong K392 selectivity → precision medicine, strongest patent
  - Delta +0.10-0.25: Moderate selectivity → genotype-guided, broad market
  - Delta ~0: Non-selective → general ERAP2 inhibitor, simplest clinical path

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/v2_k392_n392_test.py

Cost: ~$0.30 on RTX 4090 (~10 min)
"""
import glob
import json
import math
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v2_variant_selectivity")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
CIF_DIR = os.path.join(RESULTS_DIR, "cif_files")

DIFFUSION_SAMPLES = 3
SEED = 42

# Target region — ERAP2 divergent channel where position 392 sits
TARGET_REGION = (350, 500)

# ======================================================================
# V2 Synthesis Constructs (from v2/v2_reference.md)
# ======================================================================

CONSTRUCTS = [
    {
        "name": "n248_trim_c5_Y87A_Y89A",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
        "known_k392": 0.748,
        "role": "PRIMARY LEAD",
    },
    {
        "name": "n248_trim_c5",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "known_k392": 0.691,
        "role": "Parent comparator",
    },
    {
        "name": "n248_wt",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "known_k392": 0.475,
        "role": "Full-length comparator",
    },
    {
        "name": "n248_ko_all_aromatics",
        "sequence": "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
        "known_k392": 0.315,
        "role": "Negative control",
    },
]

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_region_sequence(pdb_path, start, end):
    """Extract amino acid sequence from PDB for a given residue range."""
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


def make_n392_variant(k392_seq):
    """Mutate K392 -> N392 in the extracted ERAP2 sequence.

    Position 392 in full ERAP2 maps to index (392 - TARGET_REGION[0]) in our fragment.
    """
    idx = 392 - TARGET_REGION[0]
    seq_list = list(k392_seq)
    assert seq_list[idx] == "K", f"Expected K at position 392 (index {idx}), got {seq_list[idx]}"
    seq_list[idx] = "N"
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


def find_boltz_cmd():
    for cmd in ["/opt/conda/bin/boltz", "boltz"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    print("ERROR: boltz not found")
    sys.exit(1)


def get_env():
    env = os.environ.copy()
    try:
        import nvidia.cuda_nvrtc
        nvrtc_lib = os.path.join(os.path.dirname(nvidia.cuda_nvrtc.__file__), "lib")
        env["LD_LIBRARY_PATH"] = nvrtc_lib + ":" + env.get("LD_LIBRARY_PATH", "")
    except ImportError:
        pass
    env["CC"] = "/usr/bin/gcc"
    return env


def run_boltz(boltz_cmd, yaml_path, output_dir, diffusion_samples=3):
    env = get_env()
    return subprocess.run(
        [boltz_cmd, "predict", yaml_path,
         "--out_dir", output_dir,
         "--recycling_steps", "3",
         "--diffusion_samples", str(diffusion_samples),
         "--seed", str(SEED),
         "--accelerator", "gpu",
         "--devices", "1"],
        capture_output=True, text=True, timeout=1800, env=env
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


def collect_cif_files(pred_dir, construct_name, variant_name):
    """Copy CIF files for downstream PyRosetta analysis."""
    os.makedirs(CIF_DIR, exist_ok=True)
    for cif in glob.glob(os.path.join(pred_dir, "**", "*.cif"), recursive=True):
        dest = os.path.join(CIF_DIR, f"{construct_name}_{variant_name}.cif")
        import shutil
        shutil.copy2(cif, dest)
        return dest
    return None


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(CIF_DIR, exist_ok=True)
    t_start = time.time()

    boltz_cmd = find_boltz_cmd()

    # ======================================================================
    # Build K392 and N392 target sequences
    # ======================================================================
    erap2_pdb = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
    if not os.path.exists(erap2_pdb):
        print(f"ERROR: ERAP2 structure not found at {erap2_pdb}")
        print("Upload via: python scripts/vast_launch.py upload")
        sys.exit(1)

    k392_seq = get_region_sequence(erap2_pdb, *TARGET_REGION)
    n392_seq = make_n392_variant(k392_seq)

    print("=" * 90)
    print("V2 K392 vs N392 VARIANT SELECTIVITY TEST")
    print(f"4 constructs x 2 variants x {DIFFUSION_SAMPLES} diffusion samples (seed {SEED})")
    print("=" * 90)
    print(f"\n  K392 target: {len(k392_seq)} aa (position 392 = K)")
    print(f"  N392 target: {len(n392_seq)} aa (position 392 = N)")
    idx = 392 - TARGET_REGION[0]
    print(f"  Mutation site: index {idx} — K392: ...{k392_seq[idx-3:idx+4]}... -> N392: ...{n392_seq[idx-3:idx+4]}...")

    variants = {"k392": k392_seq, "n392": n392_seq}

    # ======================================================================
    # Run predictions
    # ======================================================================
    all_results = []

    for ci, construct in enumerate(CONSTRUCTS):
        cname = construct["name"]
        cseq = construct["sequence"]
        role = construct["role"]

        print(f"\n{'='*90}")
        print(f"[{ci+1}/{len(CONSTRUCTS)}] {cname} ({len(cseq)} aa) — {role}")
        print(f"{'='*90}")

        result = {
            "name": cname,
            "sequence": cseq,
            "length": len(cseq),
            "role": role,
            "known_k392": construct["known_k392"],
            "variants": {},
        }

        for vname, vseq in variants.items():
            print(f"\n  vs ERAP2-{vname.upper()}...", end=" ", flush=True)

            yaml_path = os.path.join(RESULTS_DIR, f"{cname}_{vname}.yaml")
            pred_dir = os.path.join(RESULTS_DIR, f"{cname}_{vname}")
            write_boltz_yaml(vseq, cseq, yaml_path)

            t0 = time.time()
            r = run_boltz(boltz_cmd, yaml_path, pred_dir, DIFFUSION_SAMPLES)
            elapsed = time.time() - t0

            if r.returncode == 0:
                scores = parse_boltz_scores(pred_dir)
                iptm = scores.get("iptm", 0)
                ptm = scores.get("ptm", 0)
                plddt = scores.get("complex_plddt", 0)
                print(f"OK ({elapsed:.0f}s) ipTM={iptm:.4f}")

                cif_path = collect_cif_files(pred_dir, cname, vname)
                result["variants"][vname] = {
                    "iptm": iptm, "ptm": ptm, "plddt": plddt,
                    "cif_saved": cif_path is not None,
                }
            else:
                print(f"FAILED ({elapsed:.0f}s)")
                if r.stderr:
                    print(f"    stderr: {r.stderr[:300]}")
                result["variants"][vname] = {"iptm": 0, "error": True}

        all_results.append(result)

    # ======================================================================
    # Results Table
    # ======================================================================
    elapsed_total = time.time() - t_start

    print(f"\n\n{'='*100}")
    print(f"V2 VARIANT SELECTIVITY RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*100}")
    print(f"{'Construct':<30} {'K392 ipTM':>10} {'N392 ipTM':>10} {'Delta(K-N)':>12} {'Selective?':>12}")
    print("-" * 100)

    for r in all_results:
        k_iptm = r["variants"].get("k392", {}).get("iptm", 0)
        n_iptm = r["variants"].get("n392", {}).get("iptm", 0)
        delta = k_iptm - n_iptm

        if abs(delta) >= 0.30:
            sel = "STRONG K392" if delta > 0 else "STRONG N392"
        elif abs(delta) >= 0.10:
            sel = "MODERATE K" if delta > 0 else "MODERATE N"
        elif abs(delta) >= 0.05:
            sel = "WEAK"
        else:
            sel = "NONE"

        print(f"{r['name']:<30} {k_iptm:>10.4f} {n_iptm:>10.4f} {delta:>+12.4f} {sel:>12}")

    # ======================================================================
    # Comparison with previously known K392 values
    # ======================================================================
    print(f"\n{'='*100}")
    print("K392 REPRODUCIBILITY CHECK (new 3-sample vs previous 3-sample)")
    print(f"{'='*100}")
    print(f"{'Construct':<30} {'Previous K392':>14} {'New K392':>10} {'Shift':>10}")
    print("-" * 100)

    for r in all_results:
        prev = r["known_k392"]
        new = r["variants"].get("k392", {}).get("iptm", 0)
        shift = new - prev
        flag = " << UNSTABLE" if abs(shift) > 0.10 else ""
        print(f"{r['name']:<30} {prev:>14.4f} {new:>10.4f} {shift:>+10.4f}{flag}")

    # ======================================================================
    # Interpretation
    # ======================================================================
    lead = next((r for r in all_results if r["name"] == "n248_trim_c5_Y87A_Y89A"), None)
    neg = next((r for r in all_results if r["name"] == "n248_ko_all_aromatics"), None)

    if lead:
        k = lead["variants"].get("k392", {}).get("iptm", 0)
        n = lead["variants"].get("n392", {}).get("iptm", 0)
        d = k - n

        print(f"\n{'='*100}")
        print("LEAD CONSTRUCT INTERPRETATION: n248_trim_c5_Y87A_Y89A")
        print(f"{'='*100}")
        print(f"  K392 ipTM:  {k:.4f}")
        print(f"  N392 ipTM:  {n:.4f}")
        print(f"  Delta:      {d:+.4f}")

        if d >= 0.30:
            print(f"\n  OUTCOME 1: STRONG K392 SELECTIVITY")
            print(f"  → Precision medicine drug for K392 homozygous (~30% European)")
            print(f"  → Strongest patent position: 'variant-selective ERAP2 inhibitor'")
            print(f"  → SPR needs both K392 + N392 recombinant protein")
        elif d >= 0.10:
            print(f"\n  OUTCOME 2: MODERATE SELECTIVITY")
            print(f"  → Genotype-guided therapy, works best in K392 carriers (~80% European)")
            print(f"  → Strong patent: 'genotype-guided ERAP2 modulator'")
            print(f"  → SPR needs both variants to quantify selectivity ratio")
        else:
            print(f"\n  OUTCOME 3: NON-SELECTIVE")
            print(f"  → General ERAP2 inhibitor, works regardless of genotype")
            print(f"  → Simpler clinical path, no companion diagnostic needed")
            print(f"  → Competes with ERAP1 inhibitors in development")

    if neg:
        neg_k = neg["variants"].get("k392", {}).get("iptm", 0)
        neg_n = neg["variants"].get("n392", {}).get("iptm", 0)
        print(f"\n  Negative control (ko_all_aromatics): K392={neg_k:.4f}, N392={neg_n:.4f}")
        if max(neg_k, neg_n) < 0.40:
            print(f"  → Control validates: aromatic knockout kills binding on both variants")

    # ======================================================================
    # Save results
    # ======================================================================
    output = {
        "experiment": "V2 K392 vs N392 variant selectivity",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": {
            "diffusion_samples": DIFFUSION_SAMPLES,
            "seed": SEED,
            "target_region": list(TARGET_REGION),
            "k392_position_index": 392 - TARGET_REGION[0],
        },
        "results": all_results,
    }

    results_path = os.path.join(RESULTS_DIR, "v2_variant_selectivity_results.json")
    with open(results_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved: {results_path}")
    print(f"CIF files saved: {CIF_DIR}/")
    print(f"Total time: {elapsed_total/60:.1f} min")


if __name__ == "__main__":
    main()
