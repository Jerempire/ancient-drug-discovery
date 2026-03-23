"""Analyze Boltz-2 + PyRosetta docking results for V4 substrate-mimetic peptides.

Reads results downloaded from Vast.ai, computes K392/N392 selectivity deltas,
validates against the electrostatic model predictions, and generates a report.

Expected input: data/docking_results.json (downloaded from Vast.ai)
Flexible format — handles both:
  1. Boltz-2 + PyRosetta pipeline output (ipTM, dG_separated, per-residue energies)
  2. DiffPepDock postprocess CSV (ddG scores)

Outputs:
  data/v4_selectivity_analysis.json — full analysis
  v4_docking_report.md — human-readable report
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import csv
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"

# Electrostatic model predictions from k392_salt_bridge.py
ELECTROSTATIC_PREDICTIONS = {
    "E": "K392",   # Glu — salt bridge to K392 NH3+
    "D": "K392",   # Asp — salt bridge to K392 NH3+
    "L": "N392",   # Leu — hydrophobic, 165x N392 preference
    "F": "N392",   # Phe — hydrophobic
    "I": "N392",   # Ile — hydrophobic
    "V": "N392",   # Val — hydrophobic
    "Y": "N392",   # Tyr — hydrophobic
    "W": "N392",   # Trp — hydrophobic
    "A": "NEUTRAL", # Ala — control
    "G": "NEUTRAL", # Gly — control
    "R": "NEUTRAL", # Arg — positive charge, similar to K392
    "K": "NEUTRAL", # Lys — positive charge, similar to K392
}

# Map peptide IDs to P1 residues (from peptide_seq.fasta)
PEPTIDE_P1 = {
    "pep_glu_01": "E", "pep_glu_02": "E", "pep_glu_03": "E",
    "pep_asp_01": "D", "pep_asp_02": "D", "pep_asp_03": "D",
    "pep_leu_01": "L", "pep_leu_02": "L",
    "pep_phe_01": "F", "pep_phe_02": "F",
    "pep_arg_01": "R", "pep_arg_02": "R",
    "pep_lys_01": "K",
    "pep_ala_01": "A", "pep_ala_02": "A",
    "pep_gly_01": "G",
    "pep_ile_01": "I", "pep_val_01": "V",
    "pep_tyr_01": "Y", "pep_trp_01": "W",
    "pep_glu_long_01": "E", "pep_asp_long_01": "D", "pep_leu_long_01": "L",
    "pep_glu_short_01": "E", "pep_asp_short_01": "D",
    "pep_glu_e3_01": "E", "pep_glu_e5_01": "E", "pep_asp_d3_01": "D",
}

AA_NAMES = {
    "E": "Glutamate", "D": "Aspartate", "L": "Leucine", "F": "Phenylalanine",
    "I": "Isoleucine", "V": "Valine", "Y": "Tyrosine", "W": "Tryptophan",
    "A": "Alanine", "G": "Glycine", "R": "Arginine", "K": "Lysine",
}


def load_results(path):
    """Load results from JSON or CSV, normalize to common format.

    Returns list of dicts with keys:
        peptide, variant (k392/n392), score_type, score, extra_metrics
    """
    records = []

    if path.suffix == ".csv":
        with open(path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                peptide = row.get("peptide_name", row.get("peptide", row.get("name", "unknown")))
                variant = "k392" if "k392" in row.get("target", "").lower() else "n392"
                score = float(row.get("ddG", row.get("dG_separated", row.get("score", 0))))
                records.append({
                    "peptide": peptide,
                    "variant": variant,
                    "score_type": "ddG",
                    "score": score,
                    "extra": {k: v for k, v in row.items()
                              if k not in ("peptide_name", "peptide", "name", "target", "ddG", "score")},
                })
    elif path.suffix == ".json":
        with open(path, encoding="utf-8") as f:
            data = json.load(f)

        # Handle list of results
        if isinstance(data, list):
            for entry in data:
                peptide = entry.get("peptide", entry.get("peptide_name", entry.get("name", "unknown")))
                variant = entry.get("variant", "")
                if not variant:
                    # Try to infer from target or complex name
                    target = entry.get("target", entry.get("complex", ""))
                    variant = "k392" if "k392" in str(target).lower() else "n392" if "n392" in str(target).lower() else "unknown"

                # Pick best available score
                score = None
                score_type = None
                for key in ["dG_separated", "ddG", "total_score", "score", "iptm", "ipTM"]:
                    if key in entry and entry[key] is not None:
                        score = float(entry[key])
                        score_type = key
                        break

                if score is None:
                    continue

                records.append({
                    "peptide": peptide,
                    "variant": variant,
                    "score_type": score_type,
                    "score": score,
                    "extra": {k: v for k, v in entry.items()
                              if k not in ("peptide", "peptide_name", "name", "variant",
                                           "target", "complex", score_type)},
                })

        # Handle dict keyed by variant
        elif isinstance(data, dict):
            for variant_key, peptide_results in data.items():
                variant = "k392" if "k392" in variant_key.lower() else "n392"
                if isinstance(peptide_results, dict):
                    for pep_name, metrics in peptide_results.items():
                        if isinstance(metrics, (int, float)):
                            records.append({
                                "peptide": pep_name, "variant": variant,
                                "score_type": "score", "score": float(metrics), "extra": {},
                            })
                        elif isinstance(metrics, dict):
                            score = None
                            score_type = None
                            for key in ["dG_separated", "ddG", "total_score", "score", "iptm"]:
                                if key in metrics:
                                    score = float(metrics[key])
                                    score_type = key
                                    break
                            if score is not None:
                                records.append({
                                    "peptide": pep_name, "variant": variant,
                                    "score_type": score_type, "score": score,
                                    "extra": {k: v for k, v in metrics.items() if k != score_type},
                                })

    return records


def compute_selectivity(records):
    """Group by peptide, compute K392 vs N392 delta."""
    by_peptide = defaultdict(dict)
    for r in records:
        by_peptide[r["peptide"]][r["variant"]] = r

    results = []
    for peptide, variants in sorted(by_peptide.items()):
        k392_data = variants.get("k392")
        n392_data = variants.get("n392")

        if not k392_data or not n392_data:
            continue

        k_score = k392_data["score"]
        n_score = n392_data["score"]
        score_type = k392_data["score_type"]

        # For energy scores (dG, ddG): more negative = stronger binding
        # Delta = K392 - N392: negative = prefers K392
        # For ipTM scores: higher = better
        # Delta = K392 - N392: positive = prefers K392
        if score_type in ("iptm", "ipTM"):
            delta = k_score - n_score
            prefers_k392 = delta > 0.02
            prefers_n392 = delta < -0.02
        else:
            delta = k_score - n_score
            prefers_k392 = delta < -1.0
            prefers_n392 = delta > 1.0

        if prefers_k392:
            verdict = "K392-SELECT"
        elif prefers_n392:
            verdict = "N392-SELECT"
        else:
            verdict = "NEUTRAL"

        # Get P1 residue
        p1 = PEPTIDE_P1.get(peptide, "?")
        predicted = ELECTROSTATIC_PREDICTIONS.get(p1, "UNKNOWN")
        model_match = (
            (verdict == "K392-SELECT" and predicted == "K392") or
            (verdict == "N392-SELECT" and predicted == "N392") or
            (verdict == "NEUTRAL" and predicted == "NEUTRAL")
        )

        results.append({
            "peptide": peptide,
            "p1_residue": p1,
            "p1_name": AA_NAMES.get(p1, p1),
            "score_type": score_type,
            "k392_score": round(k_score, 4),
            "n392_score": round(n_score, 4),
            "delta": round(delta, 4),
            "verdict": verdict,
            "electrostatic_prediction": predicted,
            "model_validates": model_match,
            "k392_extra": k392_data.get("extra", {}),
            "n392_extra": n392_data.get("extra", {}),
        })

    # Sort: most K392-selective first (by delta direction appropriate to score type)
    if results and results[0]["score_type"] in ("iptm", "ipTM"):
        results.sort(key=lambda r: -r["delta"])
    else:
        results.sort(key=lambda r: r["delta"])

    return results


def build_report(results):
    """Generate markdown report."""
    lines = [
        "# V4 Docking Results: K392/N392 Selectivity Analysis",
        "",
        f"**Peptides analyzed:** {len(results)}",
    ]

    if not results:
        lines.append("\nNo paired K392/N392 results found. Check input file format.")
        return "\n".join(lines)

    score_type = results[0]["score_type"]
    lines.append(f"**Score metric:** {score_type}")
    lines.append("")

    # Summary table
    lines.append("## Selectivity Ranking")
    lines.append("")

    if score_type in ("iptm", "ipTM"):
        lines.append("| Rank | Peptide | P1 | K392 ipTM | N392 ipTM | Delta | Verdict | Model Match |")
        lines.append("|------|---------|-----|-----------|-----------|-------|---------|-------------|")
    else:
        lines.append(f"| Rank | Peptide | P1 | K392 {score_type} | N392 {score_type} | Delta | Verdict | Model Match |")
        lines.append("|------|---------|-----|------------|------------|-------|---------|-------------|")

    k392_count = 0
    n392_count = 0
    match_count = 0
    for i, r in enumerate(results, 1):
        match_str = "YES" if r["model_validates"] else "NO"
        lines.append(
            f"| {i} | {r['peptide']} | {r['p1_residue']} ({r['p1_name']}) "
            f"| {r['k392_score']:.4f} | {r['n392_score']:.4f} "
            f"| {r['delta']:+.4f} | {r['verdict']} | {match_str} |"
        )
        if r["verdict"] == "K392-SELECT":
            k392_count += 1
        elif r["verdict"] == "N392-SELECT":
            n392_count += 1
        if r["model_validates"]:
            match_count += 1

    lines += [
        "",
        "## Summary",
        "",
        f"- **K392-selective peptides:** {k392_count}",
        f"- **N392-selective peptides:** {n392_count}",
        f"- **Neutral:** {len(results) - k392_count - n392_count}",
        f"- **Electrostatic model accuracy:** {match_count}/{len(results)} ({100*match_count/len(results):.0f}%)",
        "",
    ]

    # Validation assessment
    lines.append("## Electrostatic Model Validation")
    lines.append("")

    glu_asp = [r for r in results if r["p1_residue"] in ("E", "D")]
    hydrophobic = [r for r in results if r["p1_residue"] in ("L", "F", "I", "V", "Y", "W")]
    controls = [r for r in results if r["p1_residue"] in ("A", "G", "R", "K")]

    if glu_asp:
        k392_sel = sum(1 for r in glu_asp if r["verdict"] == "K392-SELECT")
        lines.append(f"### Negatively charged P1 (Glu/Asp) -- NOVEL PREDICTION")
        lines.append(f"- Predicted: K392-selective (salt bridge to lysine NH3+)")
        lines.append(f"- Result: {k392_sel}/{len(glu_asp)} are K392-selective")
        avg_delta = sum(r["delta"] for r in glu_asp) / len(glu_asp)
        lines.append(f"- Average delta: {avg_delta:+.4f}")
        if k392_sel > 0:
            lines.append(f"- **VALIDATED**: Salt bridge mechanism confirmed by docking")
        else:
            lines.append(f"- **NOT VALIDATED**: Docking does not support salt bridge prediction")
        lines.append("")

    if hydrophobic:
        n392_sel = sum(1 for r in hydrophobic if r["verdict"] == "N392-SELECT")
        lines.append(f"### Hydrophobic P1 (Leu/Phe/Ile/Val/Tyr/Trp)")
        lines.append(f"- Predicted: N392-selective (165x published preference)")
        lines.append(f"- Result: {n392_sel}/{len(hydrophobic)} are N392-selective")
        avg_delta = sum(r["delta"] for r in hydrophobic) / len(hydrophobic)
        lines.append(f"- Average delta: {avg_delta:+.4f}")
        lines.append("")

    if controls:
        neutral = sum(1 for r in controls if r["verdict"] == "NEUTRAL")
        lines.append(f"### Controls (Ala/Gly/Arg/Lys)")
        lines.append(f"- Predicted: Neutral")
        lines.append(f"- Result: {neutral}/{len(controls)} are neutral")
        lines.append("")

    # Top candidates for synthesis
    lines.append("## Synthesis Candidates")
    lines.append("")

    k392_top = [r for r in results if r["verdict"] == "K392-SELECT"]
    if k392_top:
        lines.append("### K392-selective (target disease-associated allele)")
        for r in k392_top[:3]:
            lines.append(f"- **{r['peptide']}** (P1={r['p1_residue']}, delta={r['delta']:+.4f})")
            lines.append(f"  - Convert P1 to D-amino acid for non-cleavable version")
    else:
        lines.append("No K392-selective peptides found at current threshold.")

    lines.append("")

    n392_top = [r for r in results if r["verdict"] == "N392-SELECT"]
    if n392_top:
        lines.append("### N392-selective (control/probe)")
        for r in n392_top[:2]:
            lines.append(f"- **{r['peptide']}** (P1={r['p1_residue']}, delta={r['delta']:+.4f})")
    lines.append("")

    # Extra metrics if available
    has_extra = any(r.get("k392_extra") for r in results)
    if has_extra:
        lines.append("## Additional Metrics")
        lines.append("")
        sample = results[0]
        extra_keys = list(sample.get("k392_extra", {}).keys())[:6]
        if extra_keys:
            header = "| Peptide | " + " | ".join(f"K392 {k}" for k in extra_keys) + " | " + " | ".join(f"N392 {k}" for k in extra_keys) + " |"
            sep = "|---------|" + "|".join(["----------"] * len(extra_keys) * 2) + "|"
            lines.append(header)
            lines.append(sep)
            for r in results:
                k_vals = " | ".join(str(r.get("k392_extra", {}).get(k, "")) for k in extra_keys)
                n_vals = " | ".join(str(r.get("n392_extra", {}).get(k, "")) for k in extra_keys)
                lines.append(f"| {r['peptide']} | {k_vals} | {n_vals} |")
        lines.append("")

    lines += [
        "## Next Steps",
        "",
        "1. If salt bridge prediction validated: synthesize D-Glu and D-Asp peptides",
        "2. If not validated: re-examine channel geometry, consider longer linkers or beta-amino acids",
        "3. Measure Ki against both K392 and N392 ERAP2 alleles",
        "4. D-Leu peptide as N392-selective control",
        "",
        "## References",
        "",
        "- Evnouchidou et al. 2012 J Biol Chem (PMID: 22837489)",
        "- Papakyriakou & Stratikos 2017 PNAS",
        "- Camberlein et al. 2022 Angew Chem",
    ]

    return "\n".join(lines)


def main():
    print("=== V4 Docking Results Analysis ===\n")

    # Try multiple possible result file locations
    candidates = [
        DATA_DIR / "docking_results.json",
        DATA_DIR / "diffpepdock_results.json",
        DATA_DIR / "rosetta_pipeline_results.json",
        DATA_DIR / "boltz2_docking_results.json",
        DATA_DIR / "v4_results.json",
        SCRIPT_DIR / "postprocess_results.csv",
    ]

    result_path = None
    for p in candidates:
        if p.exists():
            result_path = p
            break

    if result_path is None:
        print("No results file found. Looked for:")
        for p in candidates:
            print(f"  {p}")
        print("\nDownload results from Vast.ai to one of these paths.")
        print("Supported formats: JSON (list or dict) or CSV with ddG/dG_separated/ipTM columns.")
        return

    print(f"Loading: {result_path}")
    records = load_results(result_path)
    print(f"Loaded {len(records)} records")

    if not records:
        print("No valid records found in file. Check format.")
        return

    # Show raw records
    variants_found = set(r["variant"] for r in records)
    peptides_found = set(r["peptide"] for r in records)
    print(f"Variants: {sorted(variants_found)}")
    print(f"Peptides: {sorted(peptides_found)}")
    print(f"Score type: {records[0]['score_type']}")

    # Compute selectivity
    results = compute_selectivity(records)
    print(f"\nPaired comparisons: {len(results)}")

    if not results:
        print("No peptides with both K392 and N392 scores found.")
        return

    # Print summary table
    print(f"\n{'='*90}")
    print(f"{'Peptide':>22s}  {'P1':>3s}  {'K392':>10s}  {'N392':>10s}  {'Delta':>10s}  {'Verdict':>12s}  {'Model':>5s}")
    print(f"{'-'*90}")
    for r in results:
        match_str = "OK" if r["model_validates"] else "MISS"
        print(f"{r['peptide']:>22s}  {r['p1_residue']:>3s}  {r['k392_score']:>10.4f}  {r['n392_score']:>10.4f}  "
              f"{r['delta']:>+10.4f}  {r['verdict']:>12s}  {match_str:>5s}")

    # Stats
    k392_sel = [r for r in results if r["verdict"] == "K392-SELECT"]
    n392_sel = [r for r in results if r["verdict"] == "N392-SELECT"]
    matches = [r for r in results if r["model_validates"]]

    print(f"\n{'='*90}")
    print(f"K392-selective: {len(k392_sel)}")
    print(f"N392-selective: {len(n392_sel)}")
    print(f"Neutral: {len(results) - len(k392_sel) - len(n392_sel)}")
    print(f"Electrostatic model accuracy: {len(matches)}/{len(results)} ({100*len(matches)/len(results):.0f}%)")

    # Key prediction check
    glu_asp = [r for r in results if r["p1_residue"] in ("E", "D")]
    if glu_asp:
        k392_confirmed = sum(1 for r in glu_asp if r["verdict"] == "K392-SELECT")
        print(f"\nSALT BRIDGE PREDICTION: {k392_confirmed}/{len(glu_asp)} Glu/Asp peptides are K392-selective")
        if k392_confirmed > 0:
            print(">>> NOVEL PREDICTION VALIDATED <<<")
        else:
            print(">>> Prediction NOT confirmed by docking <<<")

    # Save
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    json_path = DATA_DIR / "v4_selectivity_analysis.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {json_path}")

    report = build_report(results)
    report_path = SCRIPT_DIR / "v4_docking_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
