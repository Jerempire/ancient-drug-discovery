"""Run remaining D-peptide sims + combine all results."""
import sys, os, json
sys.path.insert(0, "/workspace/md_docked")
from run_md_docked import (
    run_single_md, analyze_replicate, N_REPLICATES,
    RESULTS_DIR, STRUCT_DIR
)
import numpy as np

seeds = [42, 123, 7]
remaining = [
    ("D_K392", {"peptide": "DKLLLLSI", "p1": "D", "variant": "K392"}),
    ("D_N392", {"peptide": "DKLLLLSI", "p1": "D", "variant": "N392"}),
]

all_results = {}

for sys_name, info in remaining:
    pdb_path = f"{STRUCT_DIR}/{sys_name}.pdb"
    for rep_idx in range(N_REPLICATES):
        seed = seeds[rep_idx]
        case = f"{sys_name}_rep{rep_idx}"
        prefix = f"{RESULTS_DIR}/{case}_s{seed}"
        # Clean old partial files
        for ext in ["_traj.dcd", "_energy.csv", "_start.pdb"]:
            path = f"{prefix}{ext}"
            if os.path.exists(path):
                os.remove(path)
        try:
            md_result = run_single_md(pdb_path, case, seed)
            analysis = analyze_replicate(case, seed, info)
            md_result.update(analysis)
            md_result["system"] = sys_name
            md_result["peptide"] = info["peptide"]
            md_result["p1"] = info["p1"]
            md_result["variant"] = info["variant"]
            all_results[case] = md_result
        except Exception as e:
            print(f"FAILED {case}: {e}")
            import traceback
            traceback.print_exc()
            all_results[case] = {"case": case, "error": str(e)}
        # Delete DCD after analysis to save disk
        dcd = f"{prefix}_traj.dcd"
        if os.path.exists(dcd):
            os.remove(dcd)
            print(f"    Deleted DCD to save space")

# Add V results from original run
v_results = {
    "V_K392_rep0": {"system": "V_K392", "peptide": "VKLLLLSI", "p1": "V", "variant": "K392",
                    "peptide_rmsd_mean": 1.4, "channel_contacts_mean": 0.0, "bound_fraction": 0.0,
                    "salt_bridge_occupancy": 0.0},
    "V_K392_rep1": {"system": "V_K392", "peptide": "VKLLLLSI", "p1": "V", "variant": "K392",
                    "peptide_rmsd_mean": 1.5, "channel_contacts_mean": 0.1, "bound_fraction": 0.02,
                    "salt_bridge_occupancy": 0.0},
    "V_K392_rep2": {"system": "V_K392", "peptide": "VKLLLLSI", "p1": "V", "variant": "K392",
                    "peptide_rmsd_mean": 0.8, "channel_contacts_mean": 2.7, "bound_fraction": 0.29,
                    "salt_bridge_occupancy": 0.0},
    "V_N392_rep0": {"system": "V_N392", "peptide": "VKLLLLSI", "p1": "V", "variant": "N392",
                    "peptide_rmsd_mean": 0.7, "channel_contacts_mean": 0.0, "bound_fraction": 0.01,
                    "salt_bridge_occupancy": 0.0},
    "V_N392_rep1": {"system": "V_N392", "peptide": "VKLLLLSI", "p1": "V", "variant": "N392",
                    "peptide_rmsd_mean": 0.6, "channel_contacts_mean": 3.5, "bound_fraction": 0.28,
                    "salt_bridge_occupancy": 0.0},
    "V_N392_rep2": {"system": "V_N392", "peptide": "VKLLLLSI", "p1": "V", "variant": "N392",
                    "peptide_rmsd_mean": 0.5, "channel_contacts_mean": 0.0, "bound_fraction": 0.0,
                    "salt_bridge_occupancy": 0.0},
}
all_results.update(v_results)

# Summary
sep = "=" * 75
dash = "-" * 55
print(f"\n{sep}")
print("  FULL MD SUMMARY")
print(sep)
print(f"{'Case':>20s}  {'RMSD':>6s}  {'Contacts':>9s}  {'Bound%':>7s}  {'SB occ':>7s}")
print(dash)
for name in sorted(all_results):
    r = all_results[name]
    if "error" in r:
        print(f"{name:>20s}  ERROR")
        continue
    rmsd = f"{r.get('peptide_rmsd_mean', 0):.1f}"
    cont = f"{r.get('channel_contacts_mean', 0):.1f}"
    bound = f"{r.get('bound_fraction', 0) * 100:.0f}%"
    sb = f"{r.get('salt_bridge_occupancy', 0) * 100:.0f}%"
    print(f"{name:>20s}  {rmsd:>6s}  {cont:>9s}  {bound:>7s}  {sb:>7s}")

# Aggregate
print(f"\n{sep}")
print("  AGGREGATED (mean across 3 replicates)")
print(sep)
systems = {
    "V_K392": "VKLLLLSI / K392",
    "V_N392": "VKLLLLSI / N392",
    "D_K392": "DKLLLLSI / K392",
    "D_N392": "DKLLLLSI / N392",
}
for sys_name, label in systems.items():
    reps = [r for k, r in all_results.items()
            if r.get("system") == sys_name and "error" not in r]
    if not reps:
        print(f"  {label}: NO DATA")
        continue

    def avg(key):
        vals = [r[key] for r in reps if key in r and r[key] is not None]
        return np.mean(vals) if vals else None

    rmsd = avg("peptide_rmsd_mean")
    cont = avg("channel_contacts_mean")
    bound = avg("bound_fraction")
    sb = avg("salt_bridge_occupancy")
    r_s = f"{rmsd:.1f}A" if rmsd is not None else "N/A"
    c_s = f"{cont:.1f}" if cont is not None else "N/A"
    b_s = f"{bound * 100:.0f}%" if bound is not None else "N/A"
    s_s = f"{sb * 100:.0f}%" if sb is not None else "N/A"
    print(f"  {label}: RMSD={r_s}  Contacts={c_s}  Bound={b_s}  SaltBridge={s_s}")

with open(f"{RESULTS_DIR}/md_docked_results.json", "w") as f:
    json.dump(all_results, f, indent=2)
print(f"\nSaved to {RESULTS_DIR}/md_docked_results.json")
