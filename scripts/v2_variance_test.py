"""
v2_variance_test.py — Measure Boltz-2 ipTM variance across seeds for V2 TIER_1 candidates.

Runs each construct against ERAP2(K392), ERAP1, IRAP with --diffusion_samples 1
and seeds 1-5. Reports mean, std, min, max per construct×target.

Usage: python3 /workspace/scripts/v2_variance_test.py
"""
import glob, json, math, os, subprocess, sys, time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v2_variance_test")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")
TARGET_REGION = (350, 500)
SEEDS = [1, 2, 3, 4, 5]

AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def get_seq(pdb, s, e):
    from Bio.PDB import PDBParser
    p = PDBParser(QUIET=True); st = p.get_structure("s", pdb)
    return "".join(AA3TO1[r.get_resname()] for c in st.get_chains() if c.id=="A"
                   for r in c.get_residues() if r.get_resname() in AA3TO1 and s <= r.get_id()[1] <= e)

def write_yaml(t, b, path):
    with open(path, "w") as f:
        f.write(f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {t}\n      msa: empty\n  - protein:\n      id: B\n      sequence: {b}\n      msa: empty\n")

def get_env():
    env = os.environ.copy()
    try:
        import nvidia.cuda_nvrtc
        nvrtc_lib = os.path.join(os.path.dirname(nvidia.cuda_nvrtc.__file__), "lib")
        env["LD_LIBRARY_PATH"] = nvrtc_lib + ":" + env.get("LD_LIBRARY_PATH", "")
    except: pass
    env["CC"] = "/usr/bin/gcc"
    return env

def run_boltz(yaml_path, out_dir, seed):
    env = get_env()
    for cmd in ["/opt/conda/bin/boltz", "boltz"]:
        if os.path.exists(cmd) or cmd == "boltz":
            return subprocess.run([cmd, "predict", yaml_path, "--out_dir", out_dir,
                "--recycling_steps", "3", "--diffusion_samples", "1",
                "--seed", str(seed),
                "--accelerator", "gpu", "--devices", "1"],
                capture_output=True, text=True, timeout=1800, env=env)

def parse_scores(d):
    for jf in glob.glob(os.path.join(d, "**", "*.json"), recursive=True):
        if "manifest" in jf: continue
        try:
            data = json.load(open(jf))
            if "iptm" in data: return data["iptm"]
        except: pass
    return 0

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t0 = time.time()

    e2k = get_seq(ERAP2_PDB, *TARGET_REGION)
    e1 = get_seq(ERAP1_PDB, *TARGET_REGION)
    e2n = list(e2k); e2n[392-350] = "N"; e2n = "".join(e2n)
    targets = {"erap2_k392": e2k, "erap2_n392": e2n, "erap1": e1}

    # Fetch IRAP and ANPEP
    try:
        import urllib.request
        for uid, nm in [("Q9UIQ6", "irap"), ("P15144", "anpep")]:
            r = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            s = "".join(l for l in r.read().decode().strip().split("\n") if not l.startswith(">"))
            targets[nm] = s[349:550]
    except: pass

    for n, s in targets.items():
        print(f"  {n}: {len(s)} aa")

    # V2 TIER_1 constructs + n243_id_12
    constructs = [
        ("erap2v2_long_2",   "GPTASRADLVAGWRGAAGGAAGGVHGGLGAPAGPAAAGLSTGATGGAVGLIVGGVAPELGVG"),
        ("erap2v2_long_0",   "GAAASAEDVRAALLAGAGGGAGGMFGDLVRPAGPAGAGLTTGGTAGMVGTVVGATFPELGVG"),
        ("erap2v2_short_9",  "GATASAEDVRAAAAAAAREAARATFGDAVRPPDPSLAGVTTGAAAAAVGTVLGAERPELGVA"),
        ("erap2v2_medium_0", "GPVLGGGGGRGGGAGGAHGGAHGGVGDAVAPPDPSAADLATGAAEAAVGSLMGAVFPELGVE"),
        ("erap2v2_long_5",   "GPVLSGGDIRGIAAGMAHGGAHGVFGDAFRAPDGSRAGLTTGSAGGAVGTVLGAVAPELGVG"),
        ("n243_id_12",       "SKALKEFLSNLNKAEDYKNKGNLAFNNGNYSDAISFYKKSLSELNKAKTIINNDKNLKKMLDNKTYLGKIYQNLEKLVTNNLKAAQSYKNNP"),
    ]

    # Check for duplicate sequences
    seen = {}
    deduped = []
    for name, seq in constructs:
        if seq in seen:
            print(f"  NOTE: {name} has same sequence as {seen[seq]}, skipping")
            continue
        seen[seq] = name
        deduped.append((name, seq))
    constructs = deduped

    total_preds = len(constructs) * len(targets) * len(SEEDS)
    print(f"\n{'='*120}")
    print(f"V2 VARIANCE TEST — {len(constructs)} constructs x {len(targets)} targets x {len(SEEDS)} seeds = {total_preds} predictions")
    print(f"{'='*120}")

    all_results = {}
    pred_count = 0
    for ci, (cname, cseq) in enumerate(constructs):
        print(f"\n[{ci+1}/{len(constructs)}] {cname} ({len(cseq)} aa)")
        all_results[cname] = {"sequence": cseq, "length": len(cseq), "targets": {}}

        for tname, tseq in targets.items():
            scores = []
            print(f"  vs {tname}:", end=" ", flush=True)
            for seed in SEEDS:
                tag = f"{cname}_{tname}_s{seed}"
                yp = os.path.join(RESULTS_DIR, f"{tag}.yaml")
                od = os.path.join(RESULTS_DIR, tag)
                write_yaml(tseq, cseq, yp)
                t = time.time()
                br = run_boltz(yp, od, seed)
                el = time.time() - t
                pred_count += 1
                if br and br.returncode == 0:
                    ip = parse_scores(od)
                    scores.append(ip)
                    print(f"s{seed}={ip:.4f}({el:.0f}s)", end=" ", flush=True)
                else:
                    print(f"s{seed}=FAIL({el:.0f}s)", end=" ", flush=True)
            print()

            if scores:
                mean = sum(scores) / len(scores)
                std = math.sqrt(sum((x - mean)**2 for x in scores) / len(scores))
                all_results[cname]["targets"][tname] = {
                    "scores": scores,
                    "mean": round(mean, 4),
                    "std": round(std, 4),
                    "min": round(min(scores), 4),
                    "max": round(max(scores), 4),
                    "range": round(max(scores) - min(scores), 4),
                }
                print(f"    mean={mean:.4f}  std={std:.4f}  min={min(scores):.4f}  max={max(scores):.4f}  range={max(scores)-min(scores):.4f}")

    # Summary table
    elapsed = (time.time() - t0) / 60
    print(f"\n{'='*140}")
    print(f"V2 VARIANCE SUMMARY ({elapsed:.1f} min, {pred_count} predictions)")
    print(f"{'='*140}")
    print(f"{'Construct':<25} {'Target':<15} {'Seed1':>7} {'Seed2':>7} {'Seed3':>7} {'Seed4':>7} {'Seed5':>7} {'Mean':>7} {'Std':>7} {'Range':>7}")
    print("-"*140)

    for cname, cdata in all_results.items():
        for tname, tdata in cdata["targets"].items():
            scores = tdata["scores"]
            row = f"{cname:<25} {tname:<15}"
            for s in scores:
                row += f" {s:>7.4f}"
            for _ in range(5 - len(scores)):
                row += f" {'FAIL':>7}"
            row += f" {tdata['mean']:>7.4f} {tdata['std']:>7.4f} {tdata['range']:>7.4f}"
            print(row)
        print()

    # Delta analysis — per-seed deltas for each construct
    print(f"\n{'='*140}")
    print(f"DELTA ANALYSIS (per-seed K392 minus off-target)")
    print(f"{'='*140}")
    print(f"{'Construct':<25} {'Delta type':<20} {'S1':>7} {'S2':>7} {'S3':>7} {'S4':>7} {'S5':>7} {'Mean':>7} {'Std':>7}")
    print("-"*140)

    for cname, cdata in all_results.items():
        k392_scores = cdata["targets"].get("erap2_k392", {}).get("scores", [])
        if not k392_scores:
            continue
        for off_name in ["erap2_n392", "erap1", "irap", "anpep"]:
            off_scores = cdata["targets"].get(off_name, {}).get("scores", [])
            if not off_scores or len(off_scores) != len(k392_scores):
                continue
            deltas = [k - o for k, o in zip(k392_scores, off_scores)]
            dmean = sum(deltas) / len(deltas)
            dstd = math.sqrt(sum((x - dmean)**2 for x in deltas) / len(deltas))
            label = f"K392 - {off_name.replace('erap2_','')}"
            row = f"{cname:<25} {label:<20}"
            for d in deltas:
                row += f" {d:>+7.4f}"
            row += f" {dmean:>+7.4f} {dstd:>7.4f}"
            print(row)
        print()

    # Overall stats
    all_stds = [t["std"] for c in all_results.values() for t in c["targets"].values()]
    all_ranges = [t["range"] for c in all_results.values() for t in c["targets"].values()]
    if all_stds:
        print(f"Overall: mean std={sum(all_stds)/len(all_stds):.4f}, mean range={sum(all_ranges)/len(all_ranges):.4f}")
        print(f"         max std={max(all_stds):.4f}, max range={max(all_ranges):.4f}")

    with open(os.path.join(RESULTS_DIR, "v2_variance_results.json"), "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {os.path.join(RESULTS_DIR, 'v2_variance_results.json')}")

if __name__ == "__main__":
    main()
