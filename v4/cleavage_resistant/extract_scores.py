import json, glob, os, statistics

results_dir = "/workspace/v4_cr_results/"
for d in sorted(os.listdir(results_dir)):
    pattern = os.path.join(results_dir, d, f"boltz_results_{d}", "predictions", d, f"confidence_{d}_model_*.json")
    files = sorted(glob.glob(pattern))
    if files:
        iptms = [json.load(open(f)).get("iptm", 0) for f in files]
        avg = statistics.mean(iptms)
        vals = " ".join([f"{x:.3f}" for x in iptms])
        print(f"{d}: avg_ipTM={avg:.3f}  samples=[{vals}]")
    else:
        print(f"{d}: RUNNING")
