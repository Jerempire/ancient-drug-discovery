"""
N392 Mutation Analysis — Phase 3 of N392 selectivity enhancement.

Parses Boltz-2 results from the mutation panel screen and ranks constructs
by N392 selectivity (ipTM_N392 - ipTM_K392) while filtering for ERAP1/IRAP
counter-selectivity.

Usage:
  python scripts/n392_mutation_analysis.py [results_dir]

If results_dir not provided, defaults to:
  data/results/n392_selectivity/boltz_results/
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
from pathlib import Path
from collections import defaultdict

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_RESULTS = PROJECT / "data" / "results" / "n392_selectivity" / "boltz_results"
MANIFEST_PATH = PROJECT / "data" / "results" / "n392_selectivity" / "mutation_manifest.json"
OUTPUT_DIR = PROJECT / "data" / "results" / "n392_selectivity"

# Success criteria
MIN_N392_IPTM = 0.70
MAX_ERAP1_IPTM = 0.35
MAX_IRAP_IPTM = 0.50
MIN_N392_DELTA = 0.10  # N392 - K392


def parse_boltz_results(results_dir):
    """Parse confidence files from Boltz-2 output directories."""
    scores = {}  # {yaml_name: iptm}

    results_path = Path(results_dir)
    for d in sorted(results_path.iterdir()):
        if not d.is_dir():
            continue
        # Look for confidence JSON files
        conf_files = list(d.rglob("confidence_*.json"))
        if not conf_files:
            # Try predictions subdirectory
            conf_files = list(d.rglob("**/confidence*.json"))
        if not conf_files:
            continue

        iptms = []
        for cf in conf_files:
            with open(cf) as f:
                c = json.load(f)
            iptm = c.get("iptm", c.get("i_ptm", c.get("ipTM", 0)))
            if iptm:
                iptms.append(float(iptm))

        if iptms:
            scores[d.name] = sum(iptms) / len(iptms)

    return scores


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else str(DEFAULT_RESULTS)
    print(f"Parsing results from: {results_dir}")

    scores = parse_boltz_results(results_dir)
    if not scores:
        print("ERROR: No results found. Check the results directory path.")
        print(f"Expected Boltz-2 output directories in: {results_dir}")
        sys.exit(1)

    print(f"Found {len(scores)} scored predictions")

    # Load manifest for mutation metadata
    with open(MANIFEST_PATH) as f:
        manifest_data = json.load(f)
    manifest = {m["name"]: m for m in manifest_data["manifest"]}

    # Group scores by construct (mutation)
    constructs = defaultdict(dict)
    for name, iptm in scores.items():
        meta = manifest.get(name)
        if not meta:
            # Try matching without extension
            for mname, mmeta in manifest.items():
                if mname in name or name in mname:
                    meta = mmeta
                    break
        if not meta:
            print(f"  WARNING: {name} not in manifest, skipping")
            continue

        mut = meta.get("mutation", "unknown")
        target = meta.get("target", "unknown")
        constructs[mut][target] = iptm

    # Build results table
    results = []
    for mut, targets in sorted(constructs.items()):
        k392 = targets.get("erap2_k392", 0)
        n392 = targets.get("erap2_n392", 0)
        erap1 = targets.get("erap1", 0)
        irap = targets.get("irap", 0)
        delta = n392 - k392

        # Check pass criteria
        passes_n392 = n392 >= MIN_N392_IPTM
        passes_erap1 = erap1 <= MAX_ERAP1_IPTM
        passes_irap = irap <= MAX_IRAP_IPTM
        passes_delta = delta >= MIN_N392_DELTA
        passes_all = passes_n392 and passes_erap1 and passes_irap and passes_delta

        results.append({
            "mutation": mut,
            "erap2_k392": k392,
            "erap2_n392": n392,
            "erap1": erap1,
            "irap": irap,
            "delta_n392_k392": delta,
            "passes_n392": passes_n392,
            "passes_erap1": passes_erap1,
            "passes_irap": passes_irap,
            "passes_delta": passes_delta,
            "passes_all": passes_all,
        })

    # Sort by N392 delta (descending)
    results.sort(key=lambda x: -x["delta_n392_k392"])

    # Print report
    SEP = "=" * 85
    DASH = "-" * 85

    print(f"\n{SEP}")
    print("N392 SELECTIVITY MUTATION SCREEN — RESULTS")
    print(f"Success: N392>={MIN_N392_IPTM}, ERAP1<={MAX_ERAP1_IPTM}, IRAP<={MAX_IRAP_IPTM}, delta>={MIN_N392_DELTA}")
    print(SEP)

    # Wildtype baseline
    wt = next((r for r in results if r["mutation"] == "wt"), None)
    if wt:
        print(f"\nWILDTYPE BASELINE:")
        print(f"  K392={wt['erap2_k392']:.3f}  N392={wt['erap2_n392']:.3f}  "
              f"ERAP1={wt['erap1']:.3f}  IRAP={wt['irap']:.3f}  "
              f"delta={wt['delta_n392_k392']:+.3f}")

    # Full table
    print(f"\n{'Mutation':<16} {'K392':>6} {'N392':>6} {'ERAP1':>6} {'IRAP':>6} {'Delta':>7} {'Pass':>5}")
    print(DASH)
    for r in results:
        flag = " HIT" if r["passes_all"] else ""
        print(f"{r['mutation']:<16} {r['erap2_k392']:>6.3f} {r['erap2_n392']:>6.3f} "
              f"{r['erap1']:>6.3f} {r['irap']:>6.3f} {r['delta_n392_k392']:>+7.3f}{flag}")

    # Hits
    hits = [r for r in results if r["passes_all"]]
    print(f"\n{SEP}")
    print(f"HITS: {len(hits)} / {len(results)} constructs pass all criteria")
    print(SEP)

    if hits:
        print("\nRANKED HITS (by N392 delta):")
        for i, h in enumerate(hits, 1):
            print(f"  {i}. {h['mutation']}: N392={h['erap2_n392']:.3f}, "
                  f"delta={h['delta_n392_k392']:+.3f}, "
                  f"ERAP1={h['erap1']:.3f}, IRAP={h['irap']:.3f}")
        print(f"\nRECOMMENDATION: Re-run top {min(3, len(hits))} hits with --diffusion_samples 3")
    else:
        # Look for near-misses
        near = [r for r in results if r["delta_n392_k392"] > 0 and r["mutation"] != "wt"]
        if near:
            print("\nNo full hits, but N392-leaning constructs found:")
            for r in near[:5]:
                fails = []
                if not r["passes_n392"]: fails.append(f"N392={r['erap2_n392']:.3f}<{MIN_N392_IPTM}")
                if not r["passes_erap1"]: fails.append(f"ERAP1={r['erap1']:.3f}>{MAX_ERAP1_IPTM}")
                if not r["passes_irap"]: fails.append(f"IRAP={r['irap']:.3f}>{MAX_IRAP_IPTM}")
                if not r["passes_delta"]: fails.append(f"delta={r['delta_n392_k392']:+.3f}<{MIN_N392_DELTA}")
                print(f"  {r['mutation']}: delta={r['delta_n392_k392']:+.3f} | fails: {', '.join(fails)}")
            print("\nConsider: lower delta threshold to +0.05 or proceed to Approach B (ProteinMPNN)")
        else:
            print("\nNo N392-leaning constructs found.")
            print("All mutations favor K392 equally or more than wildtype.")
            print("\nRECOMMENDATION: Proceed to Approach B (ProteinMPNN conditioned on N392 complex)")

    # Per-position analysis
    print(f"\n{SEP}")
    print("PER-POSITION ANALYSIS")
    print(SEP)

    wt_delta = wt["delta_n392_k392"] if wt else 0

    for pos_str in ["37", "50"]:
        pos_results = [r for r in results if str(r.get("mutation", "")).startswith(("N38", "F51")[int(pos_str == "50")])]
        if not pos_results:
            continue
        print(f"\nPosition {pos_str} mutations (WT delta = {wt_delta:+.3f}):")
        for r in sorted(pos_results, key=lambda x: -x["delta_n392_k392"]):
            improvement = r["delta_n392_k392"] - wt_delta
            print(f"  {r['mutation']:<10}: delta={r['delta_n392_k392']:+.3f} "
                  f"(vs WT: {improvement:+.3f})")

    # Save results
    output = {
        "analysis_date": __import__("datetime").datetime.now().isoformat(),
        "results_dir": str(results_dir),
        "criteria": {
            "min_n392_iptm": MIN_N392_IPTM,
            "max_erap1_iptm": MAX_ERAP1_IPTM,
            "max_irap_iptm": MAX_IRAP_IPTM,
            "min_n392_delta": MIN_N392_DELTA,
        },
        "wildtype": wt,
        "hits": hits,
        "all_results": results,
    }
    out_path = OUTPUT_DIR / "n392_mutation_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
