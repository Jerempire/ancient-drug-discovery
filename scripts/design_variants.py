"""
design_variants.py — Generate interface knockouts + trimmed variants for n248.

A. Interface knockouts: mutate key interface residues to alanine → should kill binding
B. Rational trimming: remove N/C-terminal residues not contributing to interface

All variants screened against ERAP2 via Boltz-2 to confirm:
  - Knockouts lose binding (validates mechanism)
  - Trimmed versions retain binding (identifies minimal binder)

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/design_variants.py
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "n248_variants")
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

# n248 lead sequence
N248_SEQ = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK"
N248_LEN = len(N248_SEQ)

# Residues likely at the interface (aromatic/charged residues that could
# drive binding — these are the ones we'll knock out to test mechanism)
# Identified by: large side chains, charged/aromatic, positioned in regions
# that would face ERAP2's divergent channel
INTERFACE_CANDIDATES = {
    # (0-indexed position, residue, rationale)
    3: ("H", "Histidine — potential metal/H-bond coordination"),
    4: ("Y", "Tyrosine — large aromatic, likely interface"),
    5: ("F", "Phenylalanine — aromatic stacking"),
    11: ("Y", "Tyrosine — aromatic"),
    22: ("L", "Leucine — hydrophobic core/interface"),
    23: ("Y", "Tyrosine — aromatic"),
    31: ("F", "Phenylalanine — aromatic"),
    44: ("F", "Phenylalanine — potential key contact"),
    45: ("Y", "Tyrosine — aromatic"),
    51: ("Y", "Tyrosine — aromatic interface"),
    87: ("F", "Phenylalanine — C-terminal aromatic cluster"),
    88: ("F", "Phenylalanine — C-terminal aromatic cluster"),
    93: ("F", "Phenylalanine — near C-term"),
    94: ("E", "Glutamate — charged, potential salt bridge"),
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
        sys.exit("boltz not found")

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


def mutate(seq, positions, replacement="A"):
    """Replace residues at given positions with replacement (alanine scan)."""
    seq_list = list(seq)
    for pos in positions:
        if 0 <= pos < len(seq_list):
            seq_list[pos] = replacement
    return "".join(seq_list)


def generate_variants():
    """Generate all variants to test."""
    variants = []

    # 0. Wild-type control
    variants.append({
        "name": "n248_wildtype",
        "sequence": N248_SEQ,
        "category": "control",
        "description": "Original lead — baseline",
    })

    # A. Interface knockouts (alanine scanning)
    # Knockout 1: Core aromatics (Y4, F5, Y11, Y23)
    variants.append({
        "name": "n248_ko_core_aromatics",
        "sequence": mutate(N248_SEQ, [4, 5, 11, 23]),
        "category": "knockout",
        "description": "Mutate Y4A, F5A, Y11A, Y23A — core aromatic cluster",
    })

    # Knockout 2: C-terminal aromatics (F87, F88, F93)
    variants.append({
        "name": "n248_ko_cterm_aromatics",
        "sequence": mutate(N248_SEQ, [87, 88, 93]),
        "category": "knockout",
        "description": "Mutate F87A, F88A, F93A — C-terminal aromatic cluster",
    })

    # Knockout 3: Mid-region contacts (F31, F44, Y45, Y51)
    variants.append({
        "name": "n248_ko_mid_contacts",
        "sequence": mutate(N248_SEQ, [31, 44, 45, 51]),
        "category": "knockout",
        "description": "Mutate F31A, F44A, Y45A, Y51A — mid-region contacts",
    })

    # Knockout 4: All aromatics at once (should devastate binding)
    all_aromatic_pos = [i for i, aa in enumerate(N248_SEQ) if aa in "FYW"]
    variants.append({
        "name": "n248_ko_all_aromatics",
        "sequence": mutate(N248_SEQ, all_aromatic_pos),
        "category": "knockout",
        "description": f"Mutate all {len(all_aromatic_pos)} F/Y/W to alanine — nuclear option",
    })

    # C. Rational trimming
    # Trim 1: Remove first 5 residues (DIRHY → these may be floppy)
    variants.append({
        "name": "n248_trim_n5",
        "sequence": N248_SEQ[5:],
        "category": "trimmed",
        "description": f"Remove N-terminal 5 aa (DIRHY) — {N248_LEN-5} aa",
    })

    # Trim 2: Remove last 5 residues (YFFEK → C-terminal tail)
    variants.append({
        "name": "n248_trim_c5",
        "sequence": N248_SEQ[:-5],
        "category": "trimmed",
        "description": f"Remove C-terminal 5 aa (YFFEK) — {N248_LEN-5} aa",
    })

    # Trim 3: Remove both termini (5+5)
    variants.append({
        "name": "n248_trim_both5",
        "sequence": N248_SEQ[5:-5],
        "category": "trimmed",
        "description": f"Remove both termini 5+5 — {N248_LEN-10} aa",
    })

    # Trim 4: Aggressive trim (10+10)
    variants.append({
        "name": "n248_trim_both10",
        "sequence": N248_SEQ[10:-10],
        "category": "trimmed",
        "description": f"Remove both termini 10+10 — {N248_LEN-20} aa",
    })

    # Trim 5: Core only (remove 15+15)
    variants.append({
        "name": "n248_core_only",
        "sequence": N248_SEQ[15:-15],
        "category": "trimmed",
        "description": f"Core only, remove 15+15 — {N248_LEN-30} aa",
    })

    return variants


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    erap2_region = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap1_region = get_region_sequence(ERAP1_PDB, *TARGET_REGION)
    print(f"ERAP2 region: {len(erap2_region)} aa")
    print(f"ERAP1 region: {len(erap1_region)} aa")

    variants = generate_variants()
    print(f"\nGenerated {len(variants)} variants:")
    for v in variants:
        print(f"  {v['name']:<30} {len(v['sequence']):>3} aa  [{v['category']}] {v['description']}")

    print(f"\n{'='*80}")
    print(f"SCREENING {len(variants)} VARIANTS")
    print(f"{'='*80}")

    results = []
    for i, v in enumerate(variants):
        name = v["name"]
        seq = v["sequence"]
        print(f"\n[{i+1}/{len(variants)}] {name} ({len(seq)} aa, {v['category']})")

        result = {**v, "length": len(seq)}

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
            print(f"OK ({elapsed:.0f}s) ipTM={result['erap2_iptm']:.4f}")
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            result["erap2_iptm"] = 0

        # ERAP1 (only for wild-type and trimmed — knockouts don't need it)
        if v["category"] in ("control", "trimmed"):
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
                result["iptm_delta"] = round(result["erap2_iptm"] - result["erap1_iptm"], 4)
                print(f"OK ({elapsed:.0f}s) ipTM={result['erap1_iptm']:.4f} delta={result['iptm_delta']:+.4f}")
            else:
                print(f"FAILED")
                result["erap1_iptm"] = 0

        results.append(result)

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n{'='*90}")
    print(f"VARIANT SCREEN RESULTS ({elapsed_total/60:.1f} min)")
    print(f"{'='*90}")

    # Knockouts
    print(f"\n--- INTERFACE KNOCKOUTS (should LOSE binding) ---")
    wt_iptm = [r for r in results if r["name"] == "n248_wildtype"][0]["erap2_iptm"]
    print(f"{'Name':<30} {'Len':>4} {'E2_ipTM':>8} {'Change':>8} {'Binding killed?'}")
    for r in results:
        if r["category"] != "knockout":
            continue
        change = r["erap2_iptm"] - wt_iptm
        killed = "YES" if change < -0.15 else "partial" if change < -0.05 else "NO"
        print(f"{r['name']:<30} {r['length']:>4} {r['erap2_iptm']:>8.4f} {change:>+8.4f} {killed}")

    # Trimmed
    print(f"\n--- TRIMMED VARIANTS (should RETAIN binding) ---")
    print(f"{'Name':<30} {'Len':>4} {'E2_ipTM':>8} {'E1_ipTM':>8} {'Delta':>8} {'Viable?'}")
    for r in results:
        if r["category"] not in ("control", "trimmed"):
            continue
        e1 = r.get("erap1_iptm", 0)
        delta = r.get("iptm_delta", 0)
        viable = "YES" if r["erap2_iptm"] > 0.4 and delta > 0.1 else "maybe" if r["erap2_iptm"] > 0.3 else "NO"
        print(f"{r['name']:<30} {r['length']:>4} {r['erap2_iptm']:>8.4f} {e1:>8.4f} {delta:>+8.4f} {viable}")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "variant_results.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
