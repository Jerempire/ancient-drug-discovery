"""Analyze MD + MMGBSA results across elite four peptides.

Computes selectivity ratios (ddG between targets) and produces
a go/no-go verdict for each peptide before synthesis.

Usage: python analyze_selectivity.py <results_dir> <output.json>
"""
import sys
import os
import json
import glob

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass

PEPTIDES = ["VAGSAF", "IAFSAF", "VAWSAF", "FASGAV"]
TARGETS = ["erap2k392", "erap1", "irap"]
# From structure_manifest.json — ipTM values for comparison
IPTM_REF = {
    "VAGSAF": {"erap2k392": 0.910, "erap1": 0.486, "irap": 0.307},
    "IAFSAF": {"erap2k392": 0.838, "erap1": 0.192, "irap": 0.393},
    "VAWSAF": {"erap2k392": 0.897, "erap1": 0.162, "irap": 0.229},
    "FASGAV": {"erap2k392": 0.813, "erap1": 0.142, "irap": 0.391},
}


def load_results(results_dir):
    """Load all JSON result files."""
    data = {}
    for pep in PEPTIDES:
        data[pep] = {}
        for tgt in TARGETS:
            path = os.path.join(results_dir, f"{pep}_vs_{tgt}.json")
            if os.path.exists(path):
                with open(path) as f:
                    data[pep][tgt] = json.load(f)
            else:
                print(f"  WARNING: Missing {pep}_vs_{tgt}.json")
    return data


def analyze(data):
    """Compute selectivity metrics and verdicts."""
    report = {"peptides": {}, "ranking": [], "recommendation": ""}

    for pep in PEPTIDES:
        pep_data = data.get(pep, {})
        if not pep_data:
            continue

        entry = {"targets": {}, "selectivity": {}, "verdict": ""}

        for tgt in TARGETS:
            if tgt not in pep_data:
                continue
            r = pep_data[tgt]
            md = r.get("md", {})
            mmgbsa = r.get("mmgbsa", {})
            entry["targets"][tgt] = {
                "dg_bind": mmgbsa.get("dg_bind_mean_kcal") if mmgbsa else None,
                "dg_std": mmgbsa.get("dg_bind_std_kcal") if mmgbsa else None,
                "avg_rmsd_last_5ns": md.get("avg_rmsd_last_5ns"),
                "avg_contacts_last_5ns": md.get("avg_contacts_last_5ns"),
                "retention_verdict": r.get("retention_verdict"),
                "binding_verdict": r.get("binding_verdict"),
                "iptm": IPTM_REF.get(pep, {}).get(tgt),
            }

        # Compute selectivity (ddG = dG_offtarget - dG_ontarget)
        # Positive ddG = selective (binds off-target less strongly)
        e2 = entry["targets"].get("erap2k392", {})
        e1 = entry["targets"].get("erap1", {})
        ir = entry["targets"].get("irap", {})

        dg_e2 = e2.get("dg_bind")
        dg_e1 = e1.get("dg_bind")
        dg_ir = ir.get("dg_bind")

        if dg_e2 is not None and dg_e1 is not None:
            ddg_e1 = round(dg_e1 - dg_e2, 2)
            entry["selectivity"]["ddG_vs_erap1"] = ddg_e1
            entry["selectivity"]["erap1_selective"] = ddg_e1 > 3.0
        if dg_e2 is not None and dg_ir is not None:
            ddg_ir = round(dg_ir - dg_e2, 2)
            entry["selectivity"]["ddG_vs_irap"] = ddg_ir
            entry["selectivity"]["irap_selective"] = ddg_ir > 3.0

        # Verdict
        triple_sel = (
            entry["selectivity"].get("erap1_selective", False)
            and entry["selectivity"].get("irap_selective", False)
        )
        strong_bind = e2.get("binding_verdict") in ("STRONG", "MODERATE")
        stable = e2.get("retention_verdict") == "STABLE"

        if triple_sel and strong_bind and stable:
            entry["verdict"] = "GO_FOR_SYNTHESIS"
        elif triple_sel and (strong_bind or stable):
            entry["verdict"] = "PROMISING_NEEDS_REVIEW"
        elif not triple_sel:
            entry["verdict"] = "SELECTIVITY_CONCERN"
        else:
            entry["verdict"] = "WEAK_BINDING"

        report["peptides"][pep] = entry

    # Rank by ERAP2 binding strength
    ranked = []
    for pep in PEPTIDES:
        e = report["peptides"].get(pep, {})
        dg = e.get("targets", {}).get("erap2k392", {}).get("dg_bind")
        if dg is not None:
            ranked.append((pep, dg, e.get("verdict", "N/A")))
    ranked.sort(key=lambda x: x[1])  # Most negative = strongest
    report["ranking"] = [{"peptide": r[0], "dg_erap2": r[1], "verdict": r[2]} for r in ranked]

    # Overall recommendation
    go_peptides = [p for p, d in report["peptides"].items() if d.get("verdict") == "GO_FOR_SYNTHESIS"]
    if go_peptides:
        report["recommendation"] = f"SYNTHESIZE: {', '.join(go_peptides)}"
    else:
        promising = [p for p, d in report["peptides"].items() if "PROMISING" in d.get("verdict", "")]
        if promising:
            report["recommendation"] = f"REVIEW BEFORE SYNTHESIS: {', '.join(promising)}"
        else:
            report["recommendation"] = "CAUTION: No peptides passed all validation gates"

    return report


def print_report(report):
    """Print human-readable summary."""
    print("\n" + "=" * 70)
    print("  MD + MMGBSA SELECTIVITY REPORT")
    print("=" * 70)

    # Summary table
    print(f"\n{'Peptide':>8s}  {'dG(E2)':>8s}  {'dG(E1)':>8s}  {'dG(IR)':>8s}  {'ddG(E1)':>8s}  {'ddG(IR)':>8s}  {'RMSD':>6s}  {'Verdict'}")
    print("-" * 85)

    for pep in PEPTIDES:
        e = report["peptides"].get(pep, {})
        t = e.get("targets", {})
        s = e.get("selectivity", {})
        dg_e2 = t.get("erap2k392", {}).get("dg_bind", "N/A")
        dg_e1 = t.get("erap1", {}).get("dg_bind", "N/A")
        dg_ir = t.get("irap", {}).get("dg_bind", "N/A")
        ddg_e1 = s.get("ddG_vs_erap1", "N/A")
        ddg_ir = s.get("ddG_vs_irap", "N/A")
        rmsd = t.get("erap2k392", {}).get("avg_rmsd_last_5ns", "N/A")

        fmt = lambda x: f"{x:>8.1f}" if isinstance(x, (int, float)) else f"{str(x):>8s}"
        print(f"{pep:>8s}  {fmt(dg_e2)}  {fmt(dg_e1)}  {fmt(dg_ir)}  {fmt(ddg_e1)}  {fmt(ddg_ir)}  {fmt(rmsd)}  {e.get('verdict', 'N/A')}")

    # ipTM comparison
    print(f"\n  ipTM comparison (Boltz-2 structural confidence):")
    for pep in PEPTIDES:
        ref = IPTM_REF.get(pep, {})
        t = report["peptides"].get(pep, {}).get("targets", {})
        dg_e2 = t.get("erap2k392", {}).get("dg_bind", "?")
        print(f"    {pep}: ipTM={ref.get('erap2k392','?'):.3f}, dG={dg_e2} -> {'CONSISTENT' if isinstance(dg_e2, (int,float)) and dg_e2 < -5 else 'CHECK'}")

    # Ranking
    print(f"\n  Ranking by ERAP2 binding strength (most negative = strongest):")
    for i, r in enumerate(report.get("ranking", []), 1):
        print(f"    {i}. {r['peptide']}: dG = {r['dg_erap2']:.1f} kcal/mol ({r['verdict']})")

    print(f"\n  RECOMMENDATION: {report.get('recommendation', 'N/A')}")

    # Selectivity interpretation
    print(f"\n  Selectivity guide:")
    print(f"    ddG > +5 kcal/mol  = STRONG selectivity (>4000x preference)")
    print(f"    ddG > +3 kcal/mol  = GOOD selectivity (>100x preference)")
    print(f"    ddG > +1 kcal/mol  = WEAK selectivity (~5x preference)")
    print(f"    ddG < +1 kcal/mol  = NO meaningful selectivity")
    print("=" * 70)


def main():
    if len(sys.argv) < 3:
        print("Usage: python analyze_selectivity.py <results_dir> <output.json>")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_path = sys.argv[2]

    data = load_results(results_dir)
    report = analyze(data)
    print_report(report)

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  Full report saved: {output_path}")


if __name__ == "__main__":
    main()
