"""
complexa_diversify.py — Generate diverse Complexa families beyond n248.

Runs 3 campaigns with different length ranges and increased sampling,
then screens top hits through 4-target Boltz-2 panel (3 diffusion samples).

Usage (on Vast.ai A6000):
    # After setup_complexa.sh
    source /workspace/Proteina-Complexa/.venv/bin/activate
    python3 /workspace/scripts/complexa_diversify.py
"""
import glob
import json
import os
import subprocess
import sys
import time

WORKSPACE = "/workspace"
COMPLEXA_DIR = os.path.join(WORKSPACE, "Proteina-Complexa")
RESULTS_DIR = os.path.join(WORKSPACE, "results", "complexa_diversify")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")

# 3 campaigns with different length ranges for diversity
CAMPAIGNS = [
    {
        "name": "div_short",
        "task": "erap2_channel_short",   # 35-55aa
        "nsamples": 40,
        "replicas": 4,
    },
    {
        "name": "div_medium",
        "task": "erap2_channel_medium",  # 55-75aa
        "nsamples": 40,
        "replicas": 4,
    },
    {
        "name": "div_long",
        "task": "erap2_channel_long",    # 75-100aa (n248 was 97aa)
        "nsamples": 40,
        "replicas": 4,
    },
]

# Boltz-2 panel targets
TARGETS = {
    "erap2": {"pdb": os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb"), "region": (350, 500)},
    "erap1": {"pdb": os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb"), "region": (350, 500)},
    "irap":  {"uniprot": "Q9UIQ6", "region": (350, 550)},
    "anpep": {"uniprot": "P15144", "region": (350, 550)},
}

DIFFUSION_SAMPLES = 3

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


def run_complexa_campaign(campaign):
    """Run one Complexa generation campaign."""
    name = campaign["name"]
    task = campaign["task"]
    nsamples = campaign["nsamples"]
    replicas = campaign["replicas"]

    print(f"\n{'='*60}")
    print(f"CAMPAIGN: {name} (task={task}, {nsamples} samples x {replicas} replicas)")
    print(f"{'='*60}")

    output_dir = os.path.join(RESULTS_DIR, name)
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "complexa", "generate",
        "configs/search_binder_local_pipeline.yaml",
        f"++run_name={name}",
        f"++generation.task_name={task}",
        "++gen_njobs=1",
        f"++generation.reward_model.reward_models.af2folding.af_params_dir={COMPLEXA_DIR}/community_models/ckpts/AF2",
        f"++generation.dataloader.dataset.nres.nsamples={nsamples}",
        f"++generation.search.best_of_n.replicas={replicas}",
    ]

    t0 = time.time()
    r = subprocess.run(
        cmd,
        capture_output=True, text=True, timeout=3600,
        cwd=COMPLEXA_DIR,
    )
    elapsed = time.time() - t0

    if r.returncode == 0:
        print(f"  Completed in {elapsed/60:.1f} min")
    else:
        print(f"  FAILED after {elapsed/60:.1f} min")
        if r.stderr:
            print(f"  stderr: {r.stderr[:500]}")
        return []

    # Parse Complexa output — find reward CSVs
    candidates = []
    import csv
    for csv_path in glob.glob(os.path.join(COMPLEXA_DIR, "outputs", name, "**", "rewards*.csv"), recursive=True):
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq = row.get("sequence", "")
                iptm = float(row.get("af2_iptm", row.get("iptm", 0)))
                if seq and iptm > 0.3:  # Only keep promising candidates
                    candidates.append({
                        "name": f"{name}_{len(candidates)}",
                        "sequence": seq,
                        "length": len(seq),
                        "af2_iptm": iptm,
                        "campaign": name,
                    })

    # Also check workspace outputs
    for csv_path in glob.glob(os.path.join(WORKSPACE, "results", name, "**", "rewards*.csv"), recursive=True):
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq = row.get("sequence", "")
                iptm = float(row.get("af2_iptm", row.get("iptm", 0)))
                if seq and iptm > 0.3:
                    candidates.append({
                        "name": f"{name}_{len(candidates)}",
                        "sequence": seq,
                        "length": len(seq),
                        "af2_iptm": iptm,
                        "campaign": name,
                    })

    print(f"  Found {len(candidates)} candidates with AF2 ipTM > 0.3")

    # Sort by AF2 ipTM and take top 5 per campaign
    candidates.sort(key=lambda x: x["af2_iptm"], reverse=True)
    top = candidates[:5]
    for c in top:
        print(f"    {c['name']}: {c['length']}aa, AF2 ipTM={c['af2_iptm']:.3f}")

    return top


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
    for cmd in ["/opt/conda/bin/boltz", "boltz",
                "/workspace/Proteina-Complexa/.venv/bin/boltz"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
            return cmd
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return None


def run_boltz_panel(candidates, target_seqs, boltz_cmd):
    """Screen candidates through 4-target Boltz-2 panel with 3 diffusion samples."""
    panel_dir = os.path.join(RESULTS_DIR, "boltz_panel")
    os.makedirs(panel_dir, exist_ok=True)

    results = []
    for c in candidates:
        print(f"\n  Screening {c['name']} ({c['length']}aa)...")
        c_result = dict(c)
        c_result["boltz_screens"] = {}

        for tname, tseq in target_seqs.items():
            yaml_path = os.path.join(panel_dir, f"{c['name']}_{tname}.yaml")
            pred_dir = os.path.join(panel_dir, f"{c['name']}_{tname}")
            write_boltz_yaml(tseq, c["sequence"], yaml_path)

            env = os.environ.copy()
            r = subprocess.run(
                [boltz_cmd, "predict", yaml_path,
                 "--out_dir", pred_dir,
                 "--recycling_steps", "3",
                 "--diffusion_samples", str(DIFFUSION_SAMPLES),
                 "--accelerator", "gpu", "--devices", "1"],
                capture_output=True, text=True, timeout=600, env=env,
            )

            if r.returncode == 0:
                # Parse scores
                for jf in glob.glob(os.path.join(pred_dir, "**", "*.json"), recursive=True):
                    if "manifest" in jf:
                        continue
                    try:
                        with open(jf) as f:
                            data = json.load(f)
                        if isinstance(data, dict) and "iptm" in data:
                            c_result["boltz_screens"][tname] = data["iptm"]
                            print(f"    {tname}: ipTM={data['iptm']:.4f}")
                            break
                    except:
                        continue
            else:
                print(f"    {tname}: FAILED")
                c_result["boltz_screens"][tname] = 0

        # Compute selectivity
        e2 = c_result["boltz_screens"].get("erap2", 0)
        e1 = c_result["boltz_screens"].get("erap1", 0)
        irap = c_result["boltz_screens"].get("irap", 0)
        anpep = c_result["boltz_screens"].get("anpep", 0)
        c_result["selectivity"] = {
            "e2_e1": e2 - e1,
            "e2_irap": e2 - irap,
            "e2_anpep": e2 - anpep,
            "all_clean": e1 < 0.3 and irap < 0.3 and anpep < 0.3,
        }
        results.append(c_result)

    return results


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    print("=" * 70)
    print("COMPLEXA DIVERSIFICATION — GENERATING NEW FAMILIES")
    print("=" * 70)

    # Phase 1: Generate candidates
    all_candidates = []
    for campaign in CAMPAIGNS:
        top = run_complexa_campaign(campaign)
        all_candidates.extend(top)

    if not all_candidates:
        print("\nERROR: No candidates generated. Check Complexa setup.")
        sys.exit(1)

    print(f"\n\nTotal candidates for screening: {len(all_candidates)}")

    # Phase 2: Boltz-2 panel
    print(f"\n{'='*70}")
    print("BOLTZ-2 4-TARGET PANEL (3 diffusion samples)")
    print(f"{'='*70}")

    # Load target sequences
    target_seqs = {}
    for tname, info in TARGETS.items():
        if "pdb" in info and os.path.exists(info["pdb"]):
            target_seqs[tname] = get_region_sequence(info["pdb"], *info["region"])
        elif "uniprot" in info:
            try:
                target_seqs[tname] = fetch_uniprot_sequence(info["uniprot"], *info["region"])
            except:
                print(f"  WARNING: Failed to fetch {tname}")
    print(f"  Loaded {len(target_seqs)} targets")

    boltz_cmd = find_boltz_cmd()
    if not boltz_cmd:
        print("WARNING: Boltz-2 not found. Install with: pip install boltz")
        print("Saving Complexa results without Boltz-2 screening.")
        results = all_candidates
    else:
        results = run_boltz_panel(all_candidates, target_seqs, boltz_cmd)

    # Phase 3: Summary
    elapsed = time.time() - t_start
    print(f"\n\n{'='*90}")
    print(f"DIVERSIFICATION RESULTS ({elapsed/60:.1f} min)")
    print(f"{'='*90}")

    header = f"{'Name':<25} {'Len':>4} {'AF2':>6} {'ERAP2':>7} {'ERAP1':>7} {'IRAP':>7} {'ANPEP':>7} {'Clean':>6}"
    print(header)
    print("-" * len(header))

    clean_hits = []
    for r in sorted(results, key=lambda x: x.get("boltz_screens", {}).get("erap2", 0), reverse=True):
        bs = r.get("boltz_screens", {})
        sel = r.get("selectivity", {})
        clean = "YES" if sel.get("all_clean", False) else "no"
        print(f"{r['name']:<25} {r['length']:>4} {r['af2_iptm']:>6.3f} "
              f"{bs.get('erap2', 0):>7.3f} {bs.get('erap1', 0):>7.3f} "
              f"{bs.get('irap', 0):>7.3f} {bs.get('anpep', 0):>7.3f} {clean:>6}")
        if sel.get("all_clean", False) and bs.get("erap2", 0) > 0.4:
            clean_hits.append(r)

    print(f"\nClean hits (ERAP2 > 0.4, all off-targets < 0.3): {len(clean_hits)}")
    for h in clean_hits:
        bs = h["boltz_screens"]
        print(f"  {h['name']}: ERAP2={bs['erap2']:.3f}, {h['length']}aa, seq={h['sequence'][:30]}...")

    # Compare to n248_trim_c5_Y87A_Y89A benchmark
    print(f"\nBenchmark: n248_trim_c5_Y87A_Y89A = ERAP2 0.748, ERAP1 0.112, IRAP 0.186, ANPEP 0.183")

    # Save
    summary_path = os.path.join(RESULTS_DIR, "diversify_results.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {summary_path}")
    print(f"Total time: {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
