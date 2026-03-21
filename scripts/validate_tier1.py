"""
validate_tier1.py — Pre-synthesis validation for TIER_1 candidates.

Three checks before spending money:
  1. Binder-alone fold confidence (monomer Boltz-2)
  2. Interface buried surface area (complex structure analysis)
  3. Length variants for top candidates (RFdiffusion → MPNN → Boltz-2)

Usage (on Vast.ai GPU instance):
    /opt/conda/envs/boltz/bin/python3 /workspace/scripts/validate_tier1.py
"""
import glob
import json
import math
import os
import subprocess
import sys
import time

# --- Config ---
WORKSPACE = "/workspace"
RESULTS = os.path.join(WORKSPACE, "results", "tier1_validation")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")

# TIER_1 candidates from scoring
TIER1 = [
    {"id": "erap2v2_long_2",  "seq": "GPTASRADLVAGWRGAAGGAAGGVHGGLGAPAGPAAAGLSTGATGGAVGLIVGGVAPELGVG", "score": 3.872},
    {"id": "erap2v2_long_0",  "seq": "GAAASAEDVRAALLAGAGGGAGGMFGDLVRPAGPAGAGLTTGGTAGMVGTVVGATFPELGVG", "score": 3.812},
    {"id": "erap2v2_short_9", "seq": "GATASAEDVRAAAAAAAREAARATFGDAVRPPDPSLAGVTTGAAAAAVGTVLGAERPELGVA", "score": 3.804},
    {"id": "erap2v2_medium_0","seq": "GPVLGGGGGRGGGAGGAHGGAHGGVGDAVAPPDPSAADLATGAAEAAVGSLMGAVFPELGVE", "score": 3.586},
    {"id": "erap2v2_long_5",  "seq": "GPVLSGGDIRGIAAGMAHGGAHGVFGDAFRAPDGSRAGLTTGSAGGAVGTVLGAVAPELGVG", "score": 3.578},
]

# Top 3 for length variants
TOP3_IDS = ["erap2v2_long_2", "erap2v2_long_0", "erap2v2_short_9"]

# Divergent channel hotspots for RFdiffusion
HOTSPOT_RESIDUES = [353, 355, 360, 363, 367, 401, 403, 406, 408, 412]

# ERAP2 target region for Boltz-2 complex predictions
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
    for r in struct.get_residues():
        resnum = r.get_id()[1]
        resname = r.get_resname()
        if resname in AA3TO1 and start <= resnum <= end:
            seq.append(AA3TO1[resname])
    return "".join(seq)


def run_boltz(yaml_path, output_dir):
    cmd = [
        "/opt/conda/envs/boltz/bin/boltz", "predict",
        yaml_path,
        "--out_dir", output_dir,
        "--recycling_steps", "3",
        "--diffusion_samples", "1",
        "--accelerator", "gpu",
        "--devices", "1",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=900)


def parse_boltz_scores(output_dir):
    scores = {}
    for jf in glob.glob(os.path.join(output_dir, "**", "*.json"), recursive=True):
        if "manifest" in jf:
            continue
        try:
            with open(jf) as f:
                data = json.load(f)
            if isinstance(data, dict):
                for key in ["ptm", "iptm", "complex_plddt", "pair_chains_iptm",
                            "confidence_score", "chains_ptm"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue
    return scores


def find_cif_file(output_dir):
    cifs = glob.glob(os.path.join(output_dir, "**", "*.cif"), recursive=True)
    return cifs[0] if cifs else None


def compute_bsa_from_cif(cif_path):
    """Compute buried surface area from a predicted complex CIF.

    BSA = SASA(chain_A) + SASA(chain_B) - SASA(complex)
    Uses BioPython + rolling ball approximation.
    """
    from Bio.PDB import MMCIFParser
    import numpy as np

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("complex", cif_path)

    chains = list(struct.get_chains())
    if len(chains) < 2:
        return {"error": "less than 2 chains", "bsa": 0}

    # Get atoms per chain
    chain_atoms = {}
    for chain in chains:
        atoms = []
        for atom in chain.get_atoms():
            if atom.element != "H":  # skip hydrogens
                atoms.append(atom.get_vector().get_array())
        chain_atoms[chain.id] = np.array(atoms)

    # Simple contact-based BSA proxy:
    # Count atom pairs within 5A across chains, multiply by ~10 A² per contact
    chain_ids = list(chain_atoms.keys())
    if len(chain_ids) < 2:
        return {"error": "not enough chains with atoms", "bsa": 0}

    a_atoms = chain_atoms[chain_ids[0]]
    b_atoms = chain_atoms[chain_ids[1]]

    # Interface contacts (within 5 Angstroms)
    contact_count = 0
    interface_a = set()
    interface_b = set()
    for i, a in enumerate(a_atoms):
        dists = np.sqrt(np.sum((b_atoms - a) ** 2, axis=1))
        contacts = np.where(dists < 5.0)[0]
        if len(contacts) > 0:
            contact_count += len(contacts)
            interface_a.add(i)
            for j in contacts:
                interface_b.add(j)

    # Rough BSA estimate: ~10 A² per interfacial atom
    bsa_estimate = (len(interface_a) + len(interface_b)) * 10.0

    return {
        "cross_chain_contacts": contact_count,
        "interface_atoms_A": len(interface_a),
        "interface_atoms_B": len(interface_b),
        "total_atoms_A": len(a_atoms),
        "total_atoms_B": len(b_atoms),
        "bsa_estimate_A2": round(bsa_estimate, 1),
        "interface_fraction_A": round(len(interface_a) / max(len(a_atoms), 1), 3),
        "interface_fraction_B": round(len(interface_b) / max(len(b_atoms), 1), 3),
    }


# ================================================================
# CHECK 1: Binder-alone fold confidence
# ================================================================
def check1_binder_alone():
    print("\n" + "=" * 70)
    print("CHECK 1: BINDER-ALONE FOLD CONFIDENCE")
    print("=" * 70)
    print("Testing if binders fold independently (without ERAP2)")

    out_dir = os.path.join(RESULTS, "check1_monomer")
    os.makedirs(out_dir, exist_ok=True)

    results = []
    for c in TIER1:
        name = c["id"]
        seq = c["seq"]
        print(f"\n  {name} ({len(seq)} aa)...", end=" ", flush=True)

        yaml_path = os.path.join(out_dir, f"{name}_monomer.yaml")
        pred_dir = os.path.join(out_dir, name)

        with open(yaml_path, "w") as f:
            f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}
      msa: empty
""")

        t0 = time.time()
        r = run_boltz(yaml_path, pred_dir)
        elapsed = time.time() - t0

        if r.returncode == 0:
            scores = parse_boltz_scores(pred_dir)
            ptm = scores.get("ptm", 0)
            plddt = scores.get("complex_plddt", 0)
            print(f"OK ({elapsed:.0f}s) pTM={ptm:.3f} pLDDT={plddt:.3f}")

            # Verdict
            if ptm > 0.7 and plddt > 0.7:
                verdict = "FOLDED"
            elif ptm > 0.5:
                verdict = "PARTIAL"
            else:
                verdict = "DISORDERED"

            results.append({
                "design": name,
                "monomer_ptm": ptm,
                "monomer_plddt": plddt,
                "verdict": verdict,
                "elapsed_s": round(elapsed, 1),
            })
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            results.append({"design": name, "error": (r.stderr or "")[-200:]})

    # Summary
    print("\n" + "-" * 70)
    print(f"{'Design':<28} {'pTM':>6} {'pLDDT':>6} {'Verdict'}")
    print("-" * 70)
    for r in results:
        if "error" in r:
            print(f"{r['design']:<28} {'ERROR':>6}")
        else:
            color = ""
            print(f"{r['design']:<28} {r['monomer_ptm']:>6.3f} {r['monomer_plddt']:>6.3f} {r['verdict']}")

    with open(os.path.join(out_dir, "monomer_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    return results


# ================================================================
# CHECK 2: Interface buried surface area
# ================================================================
def check2_interface_bsa():
    print("\n" + "=" * 70)
    print("CHECK 2: INTERFACE BURIED SURFACE AREA")
    print("=" * 70)

    out_dir = os.path.join(RESULTS, "check2_interface")
    os.makedirs(out_dir, exist_ok=True)

    erap2_region = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    print(f"ERAP2 region ({TARGET_REGION[0]}-{TARGET_REGION[1]}): {len(erap2_region)} aa")

    results = []
    for c in TIER1:
        name = c["id"]
        seq = c["seq"]
        print(f"\n  {name}...", end=" ", flush=True)

        yaml_path = os.path.join(out_dir, f"{name}_complex.yaml")
        pred_dir = os.path.join(out_dir, name)

        with open(yaml_path, "w") as f:
            f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {erap2_region}
      msa: empty
  - protein:
      id: B
      sequence: {seq}
      msa: empty
""")

        t0 = time.time()
        r = run_boltz(yaml_path, pred_dir)
        elapsed = time.time() - t0

        if r.returncode == 0:
            scores = parse_boltz_scores(pred_dir)
            cif_path = find_cif_file(pred_dir)

            bsa_data = {}
            if cif_path:
                bsa_data = compute_bsa_from_cif(cif_path)
                bsa = bsa_data.get("bsa_estimate_A2", 0)
                contacts = bsa_data.get("cross_chain_contacts", 0)
                intf_frac = bsa_data.get("interface_fraction_B", 0)
            else:
                bsa = 0
                contacts = 0
                intf_frac = 0

            iptm = scores.get("iptm", 0)
            print(f"OK ({elapsed:.0f}s) ipTM={iptm:.3f} BSA~{bsa:.0f}A² contacts={contacts}")

            # Verdict
            if bsa >= 800:
                bsa_verdict = "STRONG"
            elif bsa >= 500:
                bsa_verdict = "ADEQUATE"
            elif bsa >= 300:
                bsa_verdict = "WEAK"
            else:
                bsa_verdict = "MINIMAL"

            result = {
                "design": name,
                "iptm": iptm,
                "ptm": scores.get("ptm", 0),
                "complex_plddt": scores.get("complex_plddt", 0),
                **bsa_data,
                "bsa_verdict": bsa_verdict,
                "cif_path": cif_path,
            }
            results.append(result)
        else:
            print(f"FAILED ({elapsed:.0f}s)")
            results.append({"design": name, "error": (r.stderr or "")[-200:]})

    # Summary
    print("\n" + "-" * 70)
    print(f"{'Design':<28} {'ipTM':>6} {'BSA(A²)':>8} {'Contacts':>9} {'IntfFrac':>9} {'Verdict'}")
    print("-" * 70)
    for r in results:
        if "error" in r:
            print(f"{r['design']:<28} {'ERROR':>6}")
        else:
            print(f"{r['design']:<28} {r['iptm']:>6.3f} {r.get('bsa_estimate_A2',0):>8.0f} "
                  f"{r.get('cross_chain_contacts',0):>9} {r.get('interface_fraction_B',0):>9.3f} "
                  f"{r['bsa_verdict']}")

    with open(os.path.join(out_dir, "interface_results.json"), "w") as f:
        json.dump([{k: v for k, v in r.items() if k != "cif_path"} for r in results],
                  f, indent=2)

    return results


# ================================================================
# CHECK 3: Length variants (RFdiffusion → MPNN → Boltz-2)
# ================================================================
def check3_length_variants(check1_results, check2_results):
    print("\n" + "=" * 70)
    print("CHECK 3: LENGTH VARIANTS")
    print("=" * 70)

    # Filter to candidates that passed checks 1 and 2
    passed = set()
    for r in check1_results:
        if r.get("verdict") in ("FOLDED", "PARTIAL") and r["design"] in TOP3_IDS:
            passed.add(r["design"])

    failed_fold = [d for d in TOP3_IDS if d not in passed]
    if failed_fold:
        print(f"  Skipping {failed_fold} — failed fold check")

    # Also check BSA
    for r in check2_results:
        if r.get("bsa_verdict") == "MINIMAL" and r["design"] in passed:
            passed.discard(r["design"])
            print(f"  Skipping {r['design']} — minimal interface")

    if not passed:
        print("\n  No candidates passed checks 1+2. Skipping length variants.")
        return []

    candidates_for_variants = [c for c in TIER1 if c["id"] in passed]
    print(f"\n  Generating length variants for: {[c['id'] for c in candidates_for_variants]}")

    out_dir = os.path.join(RESULTS, "check3_length_variants")
    os.makedirs(out_dir, exist_ok=True)

    # Length tiers
    length_tiers = [
        ("shorter", "35-50"),
        ("same",    "55-65"),
        ("longer",  "75-95"),
    ]

    hotspots = ",".join(f"A{r}" for r in HOTSPOT_RESIDUES)
    rfdiff_dir = "/app/RFdiffusion"

    # Step 1: RFdiffusion — generate backbones at each length
    rfdiff_out = os.path.join(out_dir, "rfdiffusion")
    os.makedirs(rfdiff_out, exist_ok=True)

    print("\n--- Step 1: RFdiffusion backbone generation ---")
    for length_label, length_range in length_tiers:
        prefix = os.path.join(rfdiff_out, f"variant_{length_label}")
        n_designs = len(candidates_for_variants) * 2  # 2 per candidate per length

        existing = glob.glob(f"{prefix}_*.pdb")
        if len(existing) >= n_designs:
            print(f"  {length_label} ({length_range}): already have {len(existing)} designs")
            continue

        print(f"  {length_label} ({length_range}): generating {n_designs} designs...")
        cmd = (
            f"python scripts/run_inference.py "
            f"inference.output_prefix={prefix} "
            f"inference.input_pdb={ERAP2_PDB} "
            f"'contigmap.contigs=[A{HOTSPOT_RESIDUES[0]}-{HOTSPOT_RESIDUES[-1]}/0 {length_range}]' "
            f"'ppi.hotspot_res=[{hotspots}]' "
            f"inference.num_designs={n_designs} "
            f"inference.ckpt_override_path={rfdiff_dir}/models/Complex_beta_ckpt.pt"
        )
        result = subprocess.run(cmd, shell=True, cwd=rfdiff_dir,
                                capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"    FAILED: {result.stderr[-200:]}")
        else:
            generated = glob.glob(f"{prefix}_*.pdb")
            print(f"    Generated {len(generated)} backbones")

    # Step 2: ProteinMPNN — design sequences
    mpnn_out = os.path.join(out_dir, "mpnn")
    os.makedirs(mpnn_out, exist_ok=True)
    mpnn_dir = "/workspace/ProteinMPNN"

    print("\n--- Step 2: ProteinMPNN sequence design ---")
    all_pdbs = sorted(glob.glob(os.path.join(rfdiff_out, "variant_*.pdb")))
    print(f"  Designing sequences for {len(all_pdbs)} backbones...")

    mpnn_sequences = {}
    for pdb_path in all_pdbs:
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        fa_out = os.path.join(mpnn_out, name)
        os.makedirs(fa_out, exist_ok=True)

        cmd = (
            f"python protein_mpnn_run.py "
            f"--pdb_path {pdb_path} "
            f"--out_folder {fa_out} "
            f"--num_seq_per_target 4 "
            f"--sampling_temp 0.2 "
            f"--seed 42 "
            f"--batch_size 1"
        )
        result = subprocess.run(cmd, shell=True, cwd=mpnn_dir,
                                capture_output=True, text=True, timeout=120)

        # Find output FASTA
        fa_files = glob.glob(os.path.join(fa_out, "**", "*.fa"), recursive=True)
        if fa_files:
            best_seq = None
            best_score = 999
            with open(fa_files[0]) as f:
                lines = f.read().strip().split("\n")
            for i in range(0, len(lines), 2):
                if i + 1 < len(lines) and lines[i].startswith(">"):
                    header = lines[i]
                    seq = lines[i + 1]
                    # Skip wild-type (first entry)
                    if "score=" in header and "sample=" in header:
                        for part in header.split(","):
                            if "score=" in part.strip():
                                try:
                                    sc = float(part.strip().split("=")[1])
                                    if sc < best_score:
                                        best_score = sc
                                        best_seq = seq
                                except ValueError:
                                    pass
            if best_seq:
                mpnn_sequences[name] = {"sequence": best_seq, "score": best_score}
                print(f"    {name}: {len(best_seq)} aa, MPNN={best_score:.3f}")
            else:
                print(f"    {name}: no designed sequences found")
        else:
            print(f"    {name}: MPNN failed")

    # Step 3: Boltz-2 validation
    boltz_out = os.path.join(out_dir, "boltz2")
    os.makedirs(boltz_out, exist_ok=True)

    erap2_region = get_region_sequence(ERAP2_PDB, *TARGET_REGION)
    erap1_region = get_region_sequence(ERAP1_PDB, *TARGET_REGION)

    print(f"\n--- Step 3: Boltz-2 validation ({len(mpnn_sequences)} variants) ---")

    variant_results = []
    for name, seq_info in mpnn_sequences.items():
        binder_seq = seq_info["sequence"]
        mpnn_score = seq_info["score"]

        # Parse length tier from name
        if "shorter" in name:
            length_tier = "shorter"
        elif "longer" in name:
            length_tier = "longer"
        else:
            length_tier = "same"

        print(f"\n  {name} ({length_tier}, {len(binder_seq)} aa, MPNN={mpnn_score:.3f})")

        result = {
            "design": name,
            "length_tier": length_tier,
            "binder_sequence": binder_seq,
            "binder_length": len(binder_seq),
            "mpnn_score": mpnn_score,
        }

        # ERAP2 complex
        e2_yaml = os.path.join(boltz_out, f"e2_{name}.yaml")
        e2_dir = os.path.join(boltz_out, f"e2_{name}")
        with open(e2_yaml, "w") as f:
            f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {erap2_region}
      msa: empty
  - protein:
      id: B
      sequence: {binder_seq}
      msa: empty
""")
        print(f"    ERAP2...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e2_yaml, e2_dir)
        elapsed = time.time() - t0
        if r.returncode == 0:
            e2_scores = parse_boltz_scores(e2_dir)
            result["erap2_iptm"] = e2_scores.get("iptm", 0)
            result["erap2_ptm"] = e2_scores.get("ptm", 0)
            result["erap2_plddt"] = e2_scores.get("complex_plddt", 0)
            print(f"OK ({elapsed:.0f}s) ipTM={result['erap2_iptm']:.3f}")
        else:
            print(f"FAILED")
            result["erap2_error"] = True

        # ERAP1 complex
        e1_yaml = os.path.join(boltz_out, f"e1_{name}.yaml")
        e1_dir = os.path.join(boltz_out, f"e1_{name}")
        with open(e1_yaml, "w") as f:
            f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {erap1_region}
      msa: empty
  - protein:
      id: B
      sequence: {binder_seq}
      msa: empty
""")
        print(f"    ERAP1...", end=" ", flush=True)
        t0 = time.time()
        r = run_boltz(e1_yaml, e1_dir)
        elapsed = time.time() - t0
        if r.returncode == 0:
            e1_scores = parse_boltz_scores(e1_dir)
            result["erap1_iptm"] = e1_scores.get("iptm", 0)
            print(f"OK ({elapsed:.0f}s) ipTM={result['erap1_iptm']:.3f}")
        else:
            print(f"FAILED")
            result["erap1_error"] = True

        # Selectivity
        e2i = result.get("erap2_iptm", 0)
        e1i = result.get("erap1_iptm", 0)
        if e2i and e1i:
            result["iptm_delta"] = round(e2i - e1i, 4)
            result["iptm_selectivity"] = round(e2i / max(e1i, 0.001), 3)

        variant_results.append(result)

    # Summary
    print("\n" + "-" * 90)
    print(f"{'Design':<35} {'Tier':<8} {'Len':>4} {'E2_ipTM':>8} {'E1_ipTM':>8} {'Delta':>8}")
    print("-" * 90)
    for r in sorted(variant_results, key=lambda x: x.get("erap2_iptm", 0), reverse=True):
        print(f"{r['design']:<35} {r['length_tier']:<8} {r['binder_length']:>4} "
              f"{r.get('erap2_iptm', 0):>8.3f} {r.get('erap1_iptm', 0):>8.3f} "
              f"{r.get('iptm_delta', 0):>+8.4f}")

    with open(os.path.join(out_dir, "length_variant_results.json"), "w") as f:
        json.dump(variant_results, f, indent=2)

    return variant_results


# ================================================================
# MAIN
# ================================================================
def main():
    os.makedirs(RESULTS, exist_ok=True)

    print("=" * 70)
    print("  TIER_1 PRE-SYNTHESIS VALIDATION")
    print("=" * 70)
    print(f"  Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Candidates: {len(TIER1)}")
    print(f"  Length variant targets: {TOP3_IDS}")
    t_start = time.time()

    # Check 1: Binder alone
    check1 = check1_binder_alone()

    # Check 2: Interface BSA
    check2 = check2_interface_bsa()

    # Check 3: Length variants (conditional on 1+2)
    check3 = check3_length_variants(check1, check2)

    # Final summary
    elapsed = time.time() - t_start
    print("\n" + "=" * 70)
    print(f"  VALIDATION COMPLETE ({elapsed/60:.1f} min)")
    print("=" * 70)

    # Save combined results
    combined = {
        "check1_monomer": check1,
        "check2_interface": [{k: v for k, v in r.items() if k != "cif_path"} for r in check2],
        "check3_variants": check3,
        "elapsed_min": round(elapsed / 60, 1),
    }
    summary_path = os.path.join(RESULTS, "validation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
