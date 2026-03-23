"""
run_rosetta_pipeline.py — Combined Boltz-2 prediction + PyRosetta analysis on Vast.ai.

Predicts ERAP2/IRAP complexes for locked synthesis candidates,
then runs full Rosetta interface analysis on each CIF.

Usage (on Vast.ai GPU):
    python3 /workspace/scripts/run_rosetta_pipeline.py

Prerequisites:
    bash /workspace/scripts/rosetta_setup.sh
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "rosetta_pipeline")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")

# ERAP2 target region (substrate channel)
ERAP2_REGION = (350, 500)

# 4 locked synthesis candidates + their sequences
CANDIDATES = [
    {
        "name": "n248_trim_c5",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "length": 92,
        "description": "Primary lead (92aa)",
        "synthesis_priority": 1,
    },
    {
        "name": "n248_trim_c5_Y4A",
        "sequence": "DIRHAFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "length": 92,
        "description": "Trim + Y4A (92aa) — highest ERAP2 binding",
        "synthesis_priority": 2,
    },
    {
        "name": "n248_wt",
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "length": 97,
        "description": "Wild-type parent (97aa) — comparator",
        "synthesis_priority": 3,
    },
    {
        "name": "n248_ko_all_aromatics",
        "sequence": "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
        "length": 97,
        "description": "All 16 F/Y/W→A (97aa) — negative control",
        "synthesis_priority": 4,
    },
]

# Targets to screen each candidate against
TARGETS = {
    "erap2": {
        "pdb": os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb"),
        "region": (350, 500),
    },
    "irap": {
        "uniprot": "Q9UIQ6",
        "region": (350, 550),
    },
}

# Only screen IRAP for the 2 candidates with cross-reactivity concern
IRAP_CANDIDATES = ["n248_trim_c5", "n248_trim_c5_Y4A"]

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# ================================================================
# Helper functions (matching existing pipeline patterns)
# ================================================================

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
        return seq[start - 1:end]
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

    return subprocess.run(
        [boltz_cmd, "predict", yaml_path, "--out_dir", output_dir,
         "--recycling_steps", "3", "--diffusion_samples", "3",
         "--seed", "42",
         "--accelerator", "gpu", "--devices", "1",
         "--no_kernels",
         ],
        capture_output=True, text=True, timeout=2400,
        env=os.environ.copy()
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
                for key in ["ptm", "iptm", "complex_plddt", "pair_chains_iptm",
                            "confidence_score"]:
                    if key in data and key not in scores:
                        scores[key] = data[key]
        except (json.JSONDecodeError, IOError):
            continue
    return scores


def find_cif_file(output_dir):
    cifs = glob.glob(os.path.join(output_dir, "**", "*.cif"), recursive=True)
    return cifs[0] if cifs else None


# ================================================================
# Rosetta Analysis (inline — avoids import path issues on Vast.ai)
# ================================================================

def run_rosetta_analysis(cif_path, candidate_name, target_name):
    """Run FastRelax + InterfaceAnalyzer + alanine scan on a CIF file."""
    from pyrosetta import init as pyrosetta_init, Pose, pose_from_pdb
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.numeric import xyzVector_double_t

    # Convert CIF to PDB if needed
    ext = os.path.splitext(cif_path)[1].lower()
    if ext == ".cif":
        try:
            import gemmi
            doc = gemmi.cif.read(cif_path)
            st = gemmi.make_structure_from_block(doc[0])
            for model in st:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.occ == 0.0:
                                atom.occ = 1.0
            pdb_path = cif_path.replace(".cif", "_converted.pdb")
            st.write_pdb(pdb_path)
        except Exception as e:
            print(f"  gemmi conversion failed ({e}), trying BioPython...")
            from Bio.PDB import MMCIFParser, PDBIO
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("c", cif_path)
            for atom in structure.get_atoms():
                if atom.get_occupancy() == 0.0:
                    atom.set_occupancy(1.0)
            pdb_path = cif_path.replace(".cif", "_converted.pdb")
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_path)
    else:
        pdb_path = cif_path

    # Load into PyRosetta
    pose = pose_from_pdb(pdb_path)
    print(f"  Loaded: {pose.total_residue()} residues, {pose.num_chains()} chains")

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")

    # --- FastRelax (1 cycle) ---
    print("  Running FastRelax (1 cycle)...")
    score_before = sfxn(pose)
    relaxed = Pose(pose)
    relax = FastRelax()
    relax.set_scorefxn(sfxn)
    relax.max_iter(200)
    relax.apply(relaxed)
    score_after = sfxn(relaxed)
    print(f"    Energy: {score_before:.1f} -> {score_after:.1f} REU")

    # --- InterfaceAnalyzer ---
    print("  Running InterfaceAnalyzer...")
    sfxn(relaxed)
    iam = InterfaceAnalyzerMover()
    iam.set_interface("A_B")
    iam.set_scorefunction(sfxn)
    iam.set_pack_separated(True)
    iam.set_compute_interface_sc(True)
    iam.set_compute_packstat(True)
    iam.apply(relaxed)

    interface = {
        "dG_separated": round(iam.get_interface_dG(), 2),
        "dSASA_int": round(iam.get_interface_delta_sasa(), 1),
        "packstat": round(iam.get_interface_packstat(), 4),
        "shape_complementarity": round(iam.get_interface_sc(), 4),
        "hbonds_int": iam.get_interface_delta_hbond_unsat(),
        "nres_int": iam.get_num_interface_residues(),
    }
    interface["per_residue_energy_int"] = round(
        interface["dG_separated"] / max(interface["nres_int"], 1), 2
    )

    print(f"    dG={interface['dG_separated']:.1f} REU, "
          f"BSA={interface['dSASA_int']:.0f} A^2, "
          f"packstat={interface['packstat']:.3f}, "
          f"Sc={interface['shape_complementarity']:.3f}")

    # --- Alanine Scan (binder chain B only) ---
    print("  Running alanine scan (chain B)...")

    def calc_binding_energy(p):
        bound = Pose(p)
        bound_score = sfxn(bound)
        unbound = Pose(p)
        trans = unbound.jump(1)
        shift = xyzVector_double_t(500.0, 0.0, 0.0)
        trans.set_translation(shift)
        unbound.set_jump(1, trans)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        packer = PackRotamersMover(sfxn)
        packer.task_factory(tf)
        packer.apply(unbound)
        return bound_score - sfxn(unbound)

    wt_dG = calc_binding_energy(relaxed)
    ala_scan = []
    pdb_info = relaxed.pdb_info()

    for i in range(1, relaxed.total_residue() + 1):
        if pdb_info and pdb_info.chain(i) != "B":
            continue
        resname = relaxed.residue(i).name3().strip()
        if resname in ("GLY", "PRO", "ALA"):
            continue

        # Quick interface check (CA distance to any chain A residue)
        is_interface = False
        try:
            ca_i = relaxed.residue(i).xyz("CA")
            for j in range(1, relaxed.total_residue() + 1):
                if pdb_info and pdb_info.chain(j) == "B":
                    continue
                ca_j = relaxed.residue(j).xyz("CA")
                if ca_i.distance(ca_j) < 12.0:
                    is_interface = True
                    break
        except Exception:
            continue

        if not is_interface:
            continue

        mutant = Pose(relaxed)
        MutateResidue(i, "ALA").apply(mutant)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        packer = PackRotamersMover(sfxn)
        packer.task_factory(tf)
        packer.apply(mutant)

        mut_dG = calc_binding_energy(mutant)
        ddG = round(mut_dG - wt_dG, 2)

        ala_scan.append({
            "resnum": pdb_info.number(i) if pdb_info else i,
            "resname": resname,
            "chain": "B",
            "ddG": ddG,
            "is_hotspot": ddG > 1.0,
        })

    hotspots = [r for r in ala_scan if r["is_hotspot"]]
    print(f"    Scanned {len(ala_scan)} interface residues, "
          f"{len(hotspots)} hotspots (ddG > 1.0)")

    # --- Flags ---
    flags = []
    dG = interface["dG_separated"]
    if dG > 0:
        flags.append(f"CRITICAL: dG_separated={dG:.1f} (positive = unfavorable)")
    elif dG > -5:
        flags.append(f"WEAK: dG_separated={dG:.1f} (> -5 REU)")
    if interface["dSASA_int"] < 500:
        flags.append(f"WEAK: BSA={interface['dSASA_int']:.0f} A^2 (< 500)")
    if interface["packstat"] < 0.50:
        flags.append(f"WEAK: packstat={interface['packstat']:.3f} (< 0.50)")

    if flags:
        print("    FLAGS:")
        for f in flags:
            print(f"      - {f}")

    return {
        "candidate": candidate_name,
        "target": target_name,
        "relax": {"score_before": round(score_before, 2), "score_after": round(score_after, 2)},
        "interface": interface,
        "ala_scan": ala_scan,
        "n_hotspots": len(hotspots),
        "hotspot_residues": hotspots,
        "flags": flags,
    }


# ================================================================
# Main Pipeline
# ================================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    print("=" * 70)
    print("ROSETTA INTERFACE ANALYSIS PIPELINE")
    print("Boltz-2 prediction + PyRosetta analysis")
    print("=" * 70)

    # --- Step 1: Load target sequences ---
    print("\n[STEP 1] Loading target sequences...")
    target_seqs = {}

    # ERAP2 from PDB
    erap2_pdb = TARGETS["erap2"]["pdb"]
    start, end = TARGETS["erap2"]["region"]
    target_seqs["erap2"] = get_region_sequence(erap2_pdb, start, end)
    print(f"  ERAP2: {len(target_seqs['erap2'])} aa (residues {start}-{end})")

    # IRAP from UniProt
    start, end = TARGETS["irap"]["region"]
    target_seqs["irap"] = fetch_alphafold_sequence(
        TARGETS["irap"]["uniprot"], start, end
    )
    if target_seqs["irap"]:
        print(f"  IRAP: {len(target_seqs['irap'])} aa (residues {start}-{end})")
    else:
        print("  IRAP: FAILED to fetch — skipping IRAP screens")

    # --- Step 2: Run Boltz-2 predictions ---
    print(f"\n[STEP 2] Running Boltz-2 predictions...")
    boltz_results = {}

    for cand in CANDIDATES:
        for target_name in ["erap2"]:
            run_name = f"{cand['name']}_vs_{target_name}"
            out_dir = os.path.join(RESULTS_DIR, "boltz2", run_name)
            yaml_path = os.path.join(RESULTS_DIR, "boltz2", f"{run_name}.yaml")
            os.makedirs(os.path.dirname(yaml_path), exist_ok=True)

            print(f"\n  Predicting: {run_name}")
            write_boltz_yaml(target_seqs[target_name], cand["sequence"], yaml_path)

            result = run_boltz(yaml_path, out_dir)
            if result.returncode != 0:
                print(f"    ERROR: {result.stderr[:200]}")
                continue

            scores = parse_boltz_scores(out_dir)
            cif_path = find_cif_file(out_dir)
            boltz_results[run_name] = {
                "scores": scores,
                "cif_path": cif_path,
                "candidate": cand["name"],
                "target": target_name,
            }
            print(f"    ipTM={scores.get('iptm', 'N/A')}, CIF={'found' if cif_path else 'MISSING'}")

        # IRAP only for specific candidates
        if cand["name"] in IRAP_CANDIDATES and target_seqs.get("irap"):
            run_name = f"{cand['name']}_vs_irap"
            out_dir = os.path.join(RESULTS_DIR, "boltz2", run_name)
            yaml_path = os.path.join(RESULTS_DIR, "boltz2", f"{run_name}.yaml")
            os.makedirs(os.path.dirname(yaml_path), exist_ok=True)

            print(f"\n  Predicting: {run_name}")
            write_boltz_yaml(target_seqs["irap"], cand["sequence"], yaml_path)

            result = run_boltz(yaml_path, out_dir)
            if result.returncode != 0:
                print(f"    ERROR: {result.stderr[:200]}")
                continue

            scores = parse_boltz_scores(out_dir)
            cif_path = find_cif_file(out_dir)
            boltz_results[run_name] = {
                "scores": scores,
                "cif_path": cif_path,
                "candidate": cand["name"],
                "target": "irap",
            }
            print(f"    ipTM={scores.get('iptm', 'N/A')}, CIF={'found' if cif_path else 'MISSING'}")

    # Save Boltz-2 results
    boltz_summary_path = os.path.join(RESULTS_DIR, "boltz2_summary.json")
    with open(boltz_summary_path, "w") as f:
        # Strip non-serializable cif_path for summary
        summary = {}
        for k, v in boltz_results.items():
            summary[k] = {
                "scores": v["scores"],
                "cif_path": v["cif_path"],
                "candidate": v["candidate"],
                "target": v["target"],
            }
        json.dump(summary, f, indent=2)
    print(f"\n  Boltz-2 summary saved to {boltz_summary_path}")

    # --- Step 3: PyRosetta analysis ---
    print(f"\n[STEP 3] Running PyRosetta interface analysis...")

    # Initialize PyRosetta once
    from pyrosetta import init as pyrosetta_init
    pyrosetta_init("-mute all -ignore_unrecognized_res -use_input_sc")

    rosetta_results = []

    for run_name, br in boltz_results.items():
        cif_path = br["cif_path"]
        if not cif_path:
            print(f"\n  SKIPPING {run_name}: no CIF file")
            continue

        print(f"\n{'='*60}")
        print(f"  Analyzing: {run_name}")
        print(f"  CIF: {cif_path}")

        try:
            result = run_rosetta_analysis(
                cif_path, br["candidate"], br["target"]
            )
            # Add Boltz-2 scores for comparison
            result["boltz2_scores"] = br["scores"]
            rosetta_results.append(result)

            # Save individual result
            result_path = os.path.join(RESULTS_DIR, f"rosetta_{run_name}.json")
            with open(result_path, "w") as f:
                json.dump(result, f, indent=2)

        except Exception as e:
            print(f"    ERROR in Rosetta analysis: {e}")
            import traceback
            traceback.print_exc()

    # --- Step 4: Summary ---
    print(f"\n\n{'='*100}")
    print("ROSETTA INTERFACE ANALYSIS — FINAL SUMMARY")
    print(f"{'='*100}")
    print(f"{'Candidate':<24} {'Target':<8} {'Boltz ipTM':>10} {'dG_sep':>8} "
          f"{'BSA':>8} {'Pack':>6} {'Sc':>6} {'Hotspots':>9} {'Flags'}")
    print("-" * 100)

    for r in rosetta_results:
        iface = r["interface"]
        boltz_iptm = r["boltz2_scores"].get("iptm", "N/A")
        if isinstance(boltz_iptm, float):
            boltz_str = f"{boltz_iptm:.3f}"
        else:
            boltz_str = str(boltz_iptm)
        flags_str = "; ".join(f[:25] for f in r.get("flags", [])[:2])

        print(f"{r['candidate']:<24} {r['target']:<8} {boltz_str:>10} "
              f"{iface['dG_separated']:>8.1f} "
              f"{iface['dSASA_int']:>8.0f} "
              f"{iface['packstat']:>6.3f} "
              f"{iface['shape_complementarity']:>6.3f} "
              f"{r['n_hotspots']:>9} "
              f"{flags_str}")

    print("-" * 100)

    # IRAP verdict
    erap2_dGs = {}
    irap_dGs = {}
    for r in rosetta_results:
        if r["target"] == "erap2":
            erap2_dGs[r["candidate"]] = r["interface"]["dG_separated"]
        elif r["target"] == "irap":
            irap_dGs[r["candidate"]] = r["interface"]["dG_separated"]

    if irap_dGs:
        print(f"\nIRAP CROSS-REACTIVITY VERDICT:")
        for cand, irap_dG in irap_dGs.items():
            erap2_dG = erap2_dGs.get(cand, 0)
            print(f"  {cand}: ERAP2 dG={erap2_dG:.1f}, IRAP dG={irap_dG:.1f} REU")
            if irap_dG < -15:
                print(f"    -> REAL CONCERN (IRAP dG < -15 REU)")
            elif irap_dG < -5:
                print(f"    -> MARGINAL (IRAP dG between -5 and -15 REU)")
            else:
                print(f"    -> LIKELY ARTIFACT (IRAP dG > -5 REU)")

    # Save combined results
    combined_path = os.path.join(RESULTS_DIR, "rosetta_combined_results.json")
    with open(combined_path, "w") as f:
        json.dump(rosetta_results, f, indent=2)

    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed / 60:.1f} min")
    print(f"Results saved to {RESULTS_DIR}")
    print(f"\nDownload results with:")
    print(f"  scp -r -P <PORT> root@<HOST>:{RESULTS_DIR} .")


if __name__ == "__main__":
    main()
