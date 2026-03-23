"""Analyze PepMLM + Boltz-2 results and compare with hand-designed V4 peptides.

Reads PepMLM selectivity CSVs and optional Boltz-2 scores, generates a
comparison report against the hand-designed 28-peptide library.

Usage (local, after downloading results from Vast.ai):
    python analyze_pepmlm_results.py

Outputs:
    results/pepmlm_analysis_report.md
    results/pepmlm_vs_handdesigned.csv
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import csv
from pathlib import Path
from collections import Counter

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
V4_DIR = SCRIPT_DIR.parent  # v4/
HANDDESIGNED_FASTA = V4_DIR / "inputs" / "peptide_seq.fasta"

# Electrostatic model predictions (from main V4 analyze_docking_results.py)
SALT_BRIDGE_P1 = {
    "E": "K392", "D": "K392",       # negatively charged → K392 salt bridge
    "L": "N392", "F": "N392",       # hydrophobic → N392 preference
    "I": "N392", "V": "N392",
    "Y": "N392", "W": "N392",
    "A": "NEUTRAL", "G": "NEUTRAL",
    "R": "NEUTRAL", "K": "NEUTRAL",
}


def load_handdesigned():
    """Load the 28 hand-designed peptides from V4 inputs."""
    peptides = {}
    if not HANDDESIGNED_FASTA.exists():
        return peptides
    lines = HANDDESIGNED_FASTA.read_text().strip().split("\n")
    for i in range(0, len(lines), 2):
        name = lines[i][1:].strip()
        seq = lines[i + 1].strip()
        peptides[name] = seq
    return peptides


def load_boltz2_scores(boltz2_dir):
    """Load Boltz-2 ipTM scores from output directory."""
    scores = {}
    if not boltz2_dir.exists():
        return scores
    for json_file in boltz2_dir.glob("**/*scores*.json"):
        try:
            data = json.loads(json_file.read_text())
            name = json_file.stem
            iptm = None
            if isinstance(data, dict):
                iptm = data.get("iptm", data.get("ipTM"))
            if iptm is not None:
                scores[name] = float(iptm)
        except (json.JSONDecodeError, ValueError):
            continue
    return scores


def analyze():
    # Load PepMLM results
    narrow_path = RESULTS_DIR / "pepmlm_selective_narrow.csv"
    wide_path = RESULTS_DIR / "pepmlm_selective_wide.csv"

    dfs = []
    for path, label in [(narrow_path, "narrow"), (wide_path, "wide")]:
        if path.exists():
            df = pd.read_csv(path)
            df["context"] = label
            dfs.append(df)
            print(f"Loaded {len(df)} peptides from {label} context")

    if not dfs:
        print("ERROR: No PepMLM results found in results/. Run generate_peptides.py first.")
        return

    combined = pd.concat(dfs, ignore_index=True)

    # Deduplicate keeping best PPL ratio per peptide
    combined = combined.sort_values("ppl_ratio", ascending=False).drop_duplicates("peptide", keep="first")
    print(f"Total unique scored peptides: {len(combined)}")

    # Load hand-designed
    handdesigned = load_handdesigned()
    print(f"Hand-designed peptides: {len(handdesigned)}")

    # Load Boltz-2 if available
    boltz2_scores = load_boltz2_scores(RESULTS_DIR / "boltz2_out")
    if boltz2_scores:
        print(f"Boltz-2 scores loaded: {len(boltz2_scores)}")

    # ── Analysis ──────────────────────────────────────────────────────────────
    report = []
    report.append("# PepMLM vs Hand-Designed Peptide Analysis\n")

    # 1. Selectivity distribution
    sel_counts = combined["selectivity"].value_counts()
    report.append("## Selectivity Distribution\n")
    for sel, count in sel_counts.items():
        report.append(f"- **{sel}**: {count} ({100*count/len(combined):.1f}%)")
    report.append("")

    # 2. P1 residue enrichment in K392-selective
    k392_sel = combined[combined["selectivity"] == "K392-SEL"]
    if len(k392_sel) > 0:
        p1_counts = Counter(pep[0] for pep in k392_sel["peptide"])
        report.append("## P1 Residue Enrichment (K392-selective hits)\n")
        report.append("| P1 | Count | % | Salt Bridge Prediction |")
        report.append("|---|---|---|---|")
        for aa, count in p1_counts.most_common():
            pct = 100 * count / len(k392_sel)
            pred = SALT_BRIDGE_P1.get(aa, "?")
            match = "CONFIRMED" if pred == "K392" else ("unexpected" if pred == "N392" else "neutral")
            report.append(f"| {aa} | {count} | {pct:.0f}% | {pred} ({match}) |")
        report.append("")

        # Key question: are E/D enriched?
        ed_count = p1_counts.get("E", 0) + p1_counts.get("D", 0)
        ed_pct = 100 * ed_count / len(k392_sel) if len(k392_sel) > 0 else 0
        report.append(f"**Salt bridge validation**: E+D account for {ed_count}/{len(k392_sel)} "
                      f"({ed_pct:.0f}%) of K392-selective hits.")
        if ed_pct > 30:
            report.append("**STRONG support** for electrostatic model.")
        elif ed_pct > 15:
            report.append("**Moderate support** for electrostatic model.")
        else:
            report.append("**Weak/no support** — PepMLM may have found alternative selectivity mechanisms.")
        report.append("")

    # 3. Top 15 K392-selective peptides
    report.append("## Top 15 K392-Selective Peptides\n")
    report.append("| Rank | Peptide | Length | PPL_K392 | PPL_N392 | Ratio | P1 |")
    report.append("|---|---|---|---|---|---|---|")
    for i, (_, row) in enumerate(k392_sel.head(15).iterrows()):
        report.append(
            f"| {i+1} | `{row['peptide']}` | {row.get('length_k392', len(row['peptide']))} "
            f"| {row['ppl_k392']:.2f} | {row['ppl_n392']:.2f} | {row['ppl_ratio']:.2f} "
            f"| {row['peptide'][0]} |"
        )
    report.append("")

    # 4. Overlap with hand-designed
    if handdesigned:
        report.append("## Overlap with Hand-Designed Library\n")
        hd_seqs = set(handdesigned.values())
        pepmlm_seqs = set(combined["peptide"])
        overlap = hd_seqs & pepmlm_seqs
        report.append(f"- Hand-designed: {len(hd_seqs)} peptides")
        report.append(f"- PepMLM unique: {len(pepmlm_seqs)} peptides")
        report.append(f"- Exact overlap: {len(overlap)}")
        if overlap:
            for seq in overlap:
                names = [k for k, v in handdesigned.items() if v == seq]
                report.append(f"  - `{seq}` ({', '.join(names)})")
        report.append("")

        # Score hand-designed peptides in PepMLM space
        report.append("### Hand-Designed Peptides in PepMLM Score Space\n")
        report.append("| Name | Sequence | P1 | In PepMLM? | PPL Ratio |")
        report.append("|---|---|---|---|---|")
        for name, seq in sorted(handdesigned.items()):
            in_pepmlm = seq in pepmlm_seqs
            ratio = ""
            if in_pepmlm:
                match = combined[combined["peptide"] == seq]
                if len(match) > 0:
                    ratio = f"{match.iloc[0]['ppl_ratio']:.2f}"
            report.append(f"| {name} | `{seq}` | {seq[0]} | {'Yes' if in_pepmlm else 'No'} | {ratio} |")
        report.append("")

    # 5. Boltz-2 results if available
    if boltz2_scores:
        report.append("## Boltz-2 Validation (ipTM scores)\n")
        report.append("(Higher ipTM = better predicted binding)\n")
        # Parse scores into K392/N392 pairs
        k392_iptm = {}
        n392_iptm = {}
        for name, score in boltz2_scores.items():
            if "k392" in name.lower():
                pep_name = name.replace("_k392", "").replace("_scores", "")
                k392_iptm[pep_name] = score
            elif "n392" in name.lower():
                pep_name = name.replace("_n392", "").replace("_scores", "")
                n392_iptm[pep_name] = score

        if k392_iptm:
            report.append("| Peptide | ipTM_K392 | ipTM_N392 | Delta | Consistent? |")
            report.append("|---|---|---|---|---|")
            for pep in sorted(k392_iptm.keys()):
                k = k392_iptm.get(pep, 0)
                n = n392_iptm.get(pep, 0)
                delta = k - n
                consistent = "Yes" if delta > 0.05 else ("No" if delta < -0.05 else "~")
                report.append(f"| {pep} | {k:.3f} | {n:.3f} | {delta:+.3f} | {consistent} |")
            report.append("")

    # Write report
    report_text = "\n".join(report)
    report_path = RESULTS_DIR / "pepmlm_analysis_report.md"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"\nReport written to {report_path}")
    print()
    print(report_text)


if __name__ == "__main__":
    analyze()
