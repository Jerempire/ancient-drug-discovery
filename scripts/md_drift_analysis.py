"""
MD Drift Analysis — quick post-hoc check on a trajectory.

Can run on Vast.ai or locally. Reads summary.json from each MD run
and compiles into a selectivity matrix.

Usage:
  # Analyze single run:
  python md_drift_analysis.py /path/to/md_run/

  # Compile all runs into matrix:
  python md_drift_analysis.py --compile /path/to/v43_validation/
"""
import sys
import os
import json
import glob
import argparse


def analyze_single(run_dir):
    """Print analysis of a single MD run."""
    summary_path = os.path.join(run_dir, "summary.json")
    if not os.path.exists(summary_path):
        print("ERROR: No summary.json found in %s" % run_dir)
        return None

    with open(summary_path) as f:
        s = json.load(f)

    print("=" * 60)
    print("%s vs %s" % (s.get("peptide", "?"), s.get("target", "?")))
    print("=" * 60)
    print("Verdict:          %s" % s.get("verdict", "?"))
    print("RMSD (mean/max):  %.1f / %.1f A" % (s.get("rmsd_mean", 0), s.get("rmsd_max", 0)))
    print("COM drift:        %.1f A" % s.get("com_drift_final", 0))
    print("Contact fraction: %.2f" % s.get("contact_fraction_final", 0))
    print("Zinc distance:    %.1f A (min %.1f)" % (s.get("zinc_distance_mean", 99), s.get("zinc_distance_min", 99)))
    print("Region:           F=%.1f W=%.1f C=%.1f -> %s" % (
        s.get("floor_contacts_avg", 0), s.get("wall_contacts_avg", 0),
        s.get("ceiling_contacts_avg", 0), s.get("region_verdict", "?")))
    print("K392 H-bond:      %.0f%%" % (s.get("hbond_k392_occupancy", 0) * 100))
    print("Wall time:        %.2f hrs" % s.get("wall_time_hours", 0))
    return s


def compile_all(base_dir):
    """Compile all terminal results into a selectivity matrix."""
    summaries = []

    for terminal in ["terminal_b", "terminal_c", "terminal_d"]:
        term_dir = os.path.join(base_dir, terminal)
        if not os.path.exists(term_dir):
            continue
        for run_dir in sorted(glob.glob(os.path.join(term_dir, "md_*"))):
            summary_path = os.path.join(run_dir, "summary.json")
            if os.path.exists(summary_path):
                with open(summary_path) as f:
                    s = json.load(f)
                s["terminal"] = terminal
                summaries.append(s)

    if not summaries:
        print("No summary.json files found in %s/terminal_*/md_*/" % base_dir)
        return

    # Build matrix
    peptides = sorted(set(s["peptide"] for s in summaries))
    targets = sorted(set(s["target"] for s in summaries))

    SEP = "=" * 90
    DASH = "-" * 90

    print(SEP)
    print("V4.3 VALIDATION — MD SELECTIVITY MATRIX")
    print(SEP)
    print("%-10s" % "Peptide", end="")
    for t in targets:
        print(" | %15s" % t, end="")
    print()
    print(DASH)

    for pep in peptides:
        print("%-10s" % pep, end="")
        for tgt in targets:
            match = [s for s in summaries if s["peptide"] == pep and s["target"] == tgt]
            if match:
                s = match[0]
                verdict = s.get("verdict", "?")
                rmsd = s.get("rmsd_mean", 99)
                print(" | %6.1fA %7s" % (rmsd, verdict), end="")
            else:
                print(" | %15s" % "—", end="")
        print()

    # Gate check
    print("\n" + SEP)
    print("GATE CHECK: FASGAV vs ERAP2")
    print(SEP)
    fasgav_erap2 = [s for s in summaries if s["peptide"] == "FASGAV" and "erap2" in s["target"]]
    if fasgav_erap2:
        s = fasgav_erap2[0]
        if s.get("verdict") == "DRIFTING" or s.get("rmsd_mean", 0) > 5.0:
            print("GATE PASSED — scrambled control drifts out (RMSD %.1f A)" % s["rmsd_mean"])
            print("Selectivity results are trustworthy.")
        elif s.get("verdict") == "LOCKED" or s.get("rmsd_mean", 0) < 2.0:
            print("*** GATE FAILED *** — scrambled control LOCKS in ERAP2 (RMSD %.1f A)" % s["rmsd_mean"])
            print("Pocket may be artificially sticky. All results suspect.")
        else:
            print("GATE AMBIGUOUS — scrambled control shifts (RMSD %.1f A)" % s["rmsd_mean"])
    else:
        print("FASGAV vs ERAP2 not found — gate not tested yet")

    # Save
    compiled_dir = os.path.join(base_dir, "compiled")
    os.makedirs(compiled_dir, exist_ok=True)

    matrix_path = os.path.join(compiled_dir, "selectivity_matrix.json")
    with open(matrix_path, "w") as f:
        json.dump(summaries, f, indent=2)

    # CSV version
    csv_path = os.path.join(compiled_dir, "all_rmsd.csv")
    with open(csv_path, "w") as f:
        f.write("peptide,target,terminal,rmsd_mean,rmsd_max,com_drift,contact_frac,zinc_dist,region_verdict,verdict\n")
        for s in summaries:
            f.write("%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%s,%s\n" % (
                s.get("peptide", ""), s.get("target", ""), s.get("terminal", ""),
                s.get("rmsd_mean", 0), s.get("rmsd_max", 0), s.get("com_drift_final", 0),
                s.get("contact_fraction_final", 0), s.get("zinc_distance_mean", 99),
                s.get("region_verdict", ""), s.get("verdict", "")))

    print("\nSaved to %s and %s" % (matrix_path, csv_path))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="MD run directory or v43_validation base dir")
    parser.add_argument("--compile", action="store_true", help="Compile all terminals into matrix")
    args = parser.parse_args()

    if args.compile:
        compile_all(args.path)
    else:
        analyze_single(args.path)


if __name__ == "__main__":
    main()
