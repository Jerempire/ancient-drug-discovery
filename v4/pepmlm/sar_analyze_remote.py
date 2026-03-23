"""Analyze SAR Boltz-2 results on remote instance."""
import json, glob

scores = {}
base = "/workspace/pepmlm/sar_results/boltz2_out"
for conf in sorted(glob.glob(f"{base}/**/confidence_*.json", recursive=True)):
    data = json.load(open(conf))
    iptm = data.get("iptm")
    if iptm is None:
        continue
    parent = conf.split("/boltz_results_")[0].split("/")[-1]
    scores.setdefault(parent, []).append(float(iptm))

avg = {k: sum(v)/len(v) for k, v in scores.items()}

pairs = {}
for name, iptm in avg.items():
    parts = name.split("_")
    variant = parts[-1]
    pep_name = "_".join(parts[:-1])
    scaffold = parts[1]
    p1 = parts[2]
    length = int(parts[3].replace("aa", ""))
    if pep_name not in pairs:
        pairs[pep_name] = {"scaffold": scaffold, "p1": p1, "length": length}
    pairs[pep_name][variant] = iptm

for scaffold in ["VKLLLL", "ATSKSK"]:
    print(f"\n### Scaffold: {scaffold}")
    header = f"{'P1':>3s}  {'Len':>4s}  {'K392':>8s}  {'N392':>8s}  {'Delta':>8s}  {'Verdict':>10s}"
    print(header)
    print("-" * len(header))
    for key in sorted(pairs, key=lambda k: (pairs[k]["p1"], pairs[k]["length"])):
        p = pairs[key]
        if p["scaffold"] != scaffold:
            continue
        k = p.get("K392", 0)
        n = p.get("N392", 0)
        d = k - n
        v = "K392-SEL" if d > 0.03 else ("N392-SEL" if d < -0.03 else "NEUTRAL")
        print(f"{p['p1']:>3s}  {p['length']:>4d}  {k:>8.3f}  {n:>8.3f}  {d:>+8.3f}  {v:>10s}")

print("\n### P1 Effect (mean delta across all scaffolds and lengths)")
p1d = {}
for p in pairs.values():
    d = p.get("K392", 0) - p.get("N392", 0)
    p1d.setdefault(p["p1"], []).append(d)
for aa in ["E", "D", "A", "L", "V"]:
    if aa in p1d:
        m = sum(p1d[aa]) / len(p1d[aa])
        direction = "K392" if m > 0.01 else ("N392" if m < -0.01 else "NEUTRAL")
        print(f"  {aa}: {m:+.4f} ({direction})")

print("\n### Length Effect (mean delta across all scaffolds and P1)")
ld = {}
for p in pairs.values():
    d = p.get("K392", 0) - p.get("N392", 0)
    ld.setdefault(p["length"], []).append(d)
for l in sorted(ld):
    m = sum(ld[l]) / len(ld[l])
    direction = "K392" if m > 0.01 else ("N392" if m < -0.01 else "NEUTRAL")
    print(f"  {l}aa: {m:+.4f} ({direction})")

# P1 x Length interaction
print("\n### P1 x Length Interaction")
print(f"{'P1':>3s}  {'8aa':>8s}  {'9aa':>8s}  {'10aa':>8s}  {'11aa':>8s}")
print("-" * 40)
pxl = {}
for p in pairs.values():
    d = p.get("K392", 0) - p.get("N392", 0)
    key = (p["p1"], p["length"])
    pxl.setdefault(key, []).append(d)
for aa in ["E", "D", "A", "L", "V"]:
    row = f"{aa:>3s}"
    for l in [8, 9, 10, 11]:
        vals = pxl.get((aa, l), [])
        if vals:
            m = sum(vals) / len(vals)
            row += f"  {m:>+8.4f}"
        else:
            row += f"  {'N/A':>8s}"
    print(row)

with open(f"{base}/sar_analysis.json", "w") as f:
    json.dump({"pairs": pairs, "averages": avg}, f, indent=2)
print(f"\nSaved. {len(pairs)} pairs, {len(avg)} scores")
