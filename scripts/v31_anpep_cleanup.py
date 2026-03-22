"""
v31_anpep_cleanup.py — V3.1 rescue optimization for n243_id_12.

Goal: Reduce ANPEP compatibility while preserving K392 binding and variant delta.
Strategy: Mutate ANPEP-likely contact residues (periphery aromatics/hydrophobics),
          preserve core 392/393 discrimination interface.

10 variants: 6 single mutants + 4 paired mutants.
Screen each against: ERAP2(K392), ERAP2(N392), ERAP1, IRAP, ANPEP.

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/v31_anpep_cleanup.py
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v31_cleanup")
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

# n243_id_12 parent sequence (92aa)
PARENT = "SKALKEFLSNLNKAEDYKNKGNLAFNNGNYSDAISFYKKSLSELNKAKTIINNDKNLKKMLDNKTYLGKIYQNLEKLVTNNLKAAQYYKNNP"


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


def make_k392n_variant(seq, region_start, pos=392):
    idx = pos - region_start
    if 0 <= idx < len(seq):
        s = list(seq)
        s[idx] = "N"
        return "".join(s)
    return seq


def mutate(seq, mutations):
    """Apply mutations. mutations = list of (0-indexed pos, new_aa)."""
    s = list(seq)
    for pos, new_aa in mutations:
        if 0 <= pos < len(s):
            s[pos] = new_aa
    return "".join(s)


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
    # Try multiple boltz locations
    for cmd in ["/workspace/Proteina-Complexa/.venv/bin/boltz", "/opt/conda/bin/boltz", "boltz"]:
        if os.path.exists(cmd):
            boltz_cmd = cmd
            break
    else:
        boltz_cmd = "boltz"

    env = os.environ.copy()
    return subprocess.run(
        [boltz_cmd, "predict", yaml_path, "--out_dir", output_dir,
         "--recycling_steps", "3", "--diffusion_samples", "1",
         "--accelerator", "gpu", "--devices", "1"],
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


def design_variants():
    """Design 10 variants targeting ANPEP reduction."""
    variants = []

    # Parent control
    variants.append({
        "name": "v31_parent",
        "sequence": PARENT,
        "mutations": [],
        "rationale": "Parent n243_id_12 — baseline",
        "category": "control",
    })

    # === TIER 1: Single mutants (6 variants) ===

    # Target C-terminal aromatics (HIGH risk, far from 392/393 core)
    # Y87S, Y88S — C-terminal YY doublet, likely ANPEP surface contact
    variants.append({
        "name": "v31_Y87S",
        "sequence": mutate(PARENT, [(86, "S")]),
        "mutations": ["Y87S"],
        "rationale": "C-term aromatic, likely ANPEP contact",
        "category": "single",
    })
    variants.append({
        "name": "v31_Y88S",
        "sequence": mutate(PARENT, [(87, "S")]),
        "mutations": ["Y88S"],
        "rationale": "C-term aromatic, likely ANPEP contact",
        "category": "single",
    })

    # Target N-terminal aromatics
    # F7Q — N-terminal phenylalanine in hydrophobic patch pos 4-8
    variants.append({
        "name": "v31_F7Q",
        "sequence": mutate(PARENT, [(6, "Q")]),
        "mutations": ["F7Q"],
        "rationale": "N-term aromatic in hydrophobic patch",
        "category": "single",
    })

    # Target mid-C-term aromatics (pos 66-71 patch)
    # Y71S — in the 66-71 hydrophobic cluster
    variants.append({
        "name": "v31_Y71S",
        "sequence": mutate(PARENT, [(70, "S")]),
        "mutations": ["Y71S"],
        "rationale": "C-term hydrophobic cluster 66-71",
        "category": "single",
    })

    # F25Q — N-terminal region aromatic
    variants.append({
        "name": "v31_F25Q",
        "sequence": mutate(PARENT, [(24, "Q")]),
        "mutations": ["F25Q"],
        "rationale": "N-term aromatic, potential off-target contact",
        "category": "single",
    })

    # L74T — C-terminal hydrophobic in patch 74-78
    variants.append({
        "name": "v31_L74T",
        "sequence": mutate(PARENT, [(73, "T")]),
        "mutations": ["L74T"],
        "rationale": "C-term hydrophobic in LEKLV patch",
        "category": "single",
    })

    # === TIER 2: Paired mutants (4 variants) ===

    # Double C-term aromatic removal
    variants.append({
        "name": "v31_Y87S_Y88S",
        "sequence": mutate(PARENT, [(86, "S"), (87, "S")]),
        "mutations": ["Y87S", "Y88S"],
        "rationale": "Remove C-terminal YY doublet entirely",
        "category": "double",
    })

    # N-term + C-term aromatic pair
    variants.append({
        "name": "v31_F7Q_Y71S",
        "sequence": mutate(PARENT, [(6, "Q"), (70, "S")]),
        "mutations": ["F7Q", "Y71S"],
        "rationale": "Remove one N-term + one C-term aromatic",
        "category": "double",
    })

    # C-term cluster cleanup
    variants.append({
        "name": "v31_Y71S_L74T",
        "sequence": mutate(PARENT, [(70, "S"), (73, "T")]),
        "mutations": ["Y71S", "L74T"],
        "rationale": "Clean up C-term hydrophobic cluster 66-78",
        "category": "double",
    })

    # Conservative: just reduce the two most peripheral aromatics
    variants.append({
        "name": "v31_F25Q_Y87S",
        "sequence": mutate(PARENT, [(24, "Q"), (86, "S")]),
        "mutations": ["F25Q", "Y87S"],
        "rationale": "Remove most peripheral N + C aromatics",
        "category": "double",
    })

    return variants


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    # Load target sequences
    print("Loading target sequences...")
    erap2_k392 = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap2_n392 = make_k392n_variant(erap2_k392, TARGET_REGION[0])
    erap1 = get_region_sequence(ERAP1_PDB, *TARGET_REGION)

    target_seqs = {"erap2_k392": erap2_k392, "erap2_n392": erap2_n392, "erap1": erap1}

    # Fetch IRAP/ANPEP
    try:
        import urllib.request
        for uid, name in [("Q9UIQ6", "irap"), ("P15144", "anpep")]:
            resp = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            lines = resp.read().decode().strip().split("\n")
            seq = "".join(l for l in lines if not l.startswith(">"))
            target_seqs[name] = seq[349:550]
    except Exception as e:
        print(f"  WARNING: {e}")

    for name, seq in target_seqs.items():
        print(f"  {name}: {len(seq)} aa")

    # Design variants
    variants = design_variants()
    print(f"\nDesigned {len(variants)} variants:")
    for v in variants:
        muts = ", ".join(v["mutations"]) if v["mutations"] else "none"
        print(f"  {v['name']:<25} [{v['category']}] {muts} — {v['rationale']}")

    # Screen
    print(f"\n{'='*100}")
    print(f"V3.1 ANPEP CLEANUP — {len(variants)} variants x {len(target_seqs)} targets")
    print(f"{'='*100}")

    results = []
    for i, v in enumerate(variants):
        name = v["name"]
        seq = v["sequence"]
        print(f"\n[{i+1}/{len(variants)}] {name} ({len(seq)} aa, {v['category']})")

        result = {**v, "length": len(seq), "screens": {}}

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

        # Compute metrics
        k = result["screens"].get("erap2_k392", 0)
        n = result["screens"].get("erap2_n392", 0)
        result["variant_delta"] = round(k - n, 4)
        result["anpep"] = result["screens"].get("anpep", 0)

        results.append(result)

    # Summary
    elapsed_total = time.time() - t_start
    parent = [r for r in results if r["category"] == "control"][0]
    p_k = parent["screens"].get("erap2_k392", 0)
    p_anpep = parent["screens"].get("anpep", 0)

    print(f"\n{'='*110}")
    print(f"V3.1 RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*110}")
    print(f"{'Name':<25} {'Cat':<8} {'K392':>7} {'N392':>7} {'Delta':>7} {'E1':>7} {'IRAP':>7} {'ANPEP':>7} {'ANPEP chg':>9} {'Verdict'}")
    print("-" * 110)

    for r in sorted(results, key=lambda x: x.get("anpep", 1)):
        s = r["screens"]
        k = s.get("erap2_k392", 0)
        anpep = s.get("anpep", 0)
        anpep_change = anpep - p_anpep
        delta = r["variant_delta"]

        # Verdict
        if k >= 0.55 and delta >= 0.10 and anpep < p_anpep - 0.05:
            verdict = "IMPROVED"
        elif k >= 0.55 and delta >= 0.10:
            verdict = "OK"
        elif k < 0.55:
            verdict = "K392 lost"
        elif delta < 0.10:
            verdict = "delta lost"
        else:
            verdict = "no gain"

        print(f"{r['name']:<25} {r['category']:<8} {k:>7.4f} {s.get('erap2_n392',0):>7.4f} "
              f"{delta:>+7.4f} {s.get('erap1',0):>7.4f} {s.get('irap',0):>7.4f} "
              f"{anpep:>7.4f} {anpep_change:>+9.4f} {verdict}")

    # Save
    with open(os.path.join(RESULTS_DIR, "v31_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {os.path.join(RESULTS_DIR, 'v31_results.json')}")


if __name__ == "__main__":
    main()
