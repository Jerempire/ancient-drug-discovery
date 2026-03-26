"""
V4.3 Results Compiler — Run after all terminals report complete.

Reads:
  - terminal_a/sasa/*.json
  - terminal_a/contact_maps/*.json
  - terminal_b/md_*/summary.json
  - terminal_c/md_*/summary.json
  - terminal_d/md_*/summary.json
  - coordination.json

Writes:
  - compiled/all_rmsd.csv
  - compiled/all_drift.csv
  - compiled/selectivity_matrix.json
  - compiled/gate_result.json
  - compiled/final_verdict.md
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import csv
from pathlib import Path
from datetime import datetime, timezone

PROJECT = Path(r"C:\Users\jmj2z\Projects\medical\ancient-drug-discovery")
V43 = PROJECT / "data" / "results" / "v43_validation"
COMPILED = V43 / "compiled"

PEPTIDES = ["VAGSAF", "IAFSAF", "VAWSAF", "FASGAV"]
TARGETS = ["erap2k392", "irap", "erap1"]
TERMINAL_MAP = {
    ("VAGSAF", "erap2k392"): "terminal_b",
    ("VAGSAF", "irap"): "terminal_b",
    ("IAFSAF", "erap2k392"): "terminal_b",
    ("IAFSAF", "irap"): "terminal_b",
    ("VAWSAF", "erap2k392"): "terminal_c",
    ("VAWSAF", "irap"): "terminal_c",
    ("FASGAV", "erap2k392"): "terminal_c",
    ("FASGAV", "irap"): "terminal_c",
    ("VAGSAF", "erap1"): "terminal_d",
    ("IAFSAF", "erap1"): "terminal_d",
    ("VAWSAF", "erap1"): "terminal_d",
    ("FASGAV", "erap1"): "terminal_d",
}


def load_md_summary(peptide, target):
    """Load summary.json for a given MD run."""
    terminal = TERMINAL_MAP.get((peptide, target))
    if not terminal:
        return None
    summary_path = V43 / terminal / f"md_{peptide}_{target}" / "summary.json"
    if not summary_path.exists():
        return None
    with open(summary_path) as f:
        return json.load(f)


def load_sasa_results():
    """Load SASA summary."""
    path = V43 / "terminal_a" / "sasa" / "sasa_summary.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


def load_contact_results():
    """Load contact map summary."""
    path = V43 / "terminal_a" / "contact_maps" / "contact_map_summary.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


def check_gate():
    """Check FASGAV gate result."""
    coord_path = V43 / "coordination.json"
    with open(coord_path) as f:
        coord = json.load(f)

    gate = coord.get("gate_passed")
    tc = coord.get("terminal_c", {})
    gate_result = tc.get("gate_result")

    # Also check the actual MD summary if available
    fasgav_md = load_md_summary("FASGAV", "erap2k392")

    result = {
        "coordination_gate_passed": gate,
        "coordination_gate_result": gate_result,
        "fasgav_md_available": fasgav_md is not None,
    }

    if fasgav_md:
        result["fasgav_rmsd_final"] = fasgav_md.get("rmsd_final")
        result["fasgav_rmsd_mean"] = fasgav_md.get("rmsd_mean")
        result["fasgav_verdict"] = fasgav_md.get("verdict")
        # Gate passes if FASGAV drifts out (RMSD > 3A)
        if fasgav_md.get("rmsd_final", 0) > 3.0:
            result["gate_verdict"] = "PASS"
        elif fasgav_md.get("rmsd_final", 0) < 2.0:
            result["gate_verdict"] = "FAIL"
        else:
            result["gate_verdict"] = "AMBIGUOUS"
    else:
        result["gate_verdict"] = "PENDING"

    return result


def main():
    COMPILED.mkdir(parents=True, exist_ok=True)

    # Check coordination status
    coord_path = V43 / "coordination.json"
    with open(coord_path) as f:
        coord = json.load(f)

    print("=== Terminal Status ===")
    for term in ["terminal_a", "terminal_b", "terminal_c", "terminal_d"]:
        t = coord.get(term, {})
        print(f"  {term}: {t.get('status', 'unknown')} | "
              f"done: {len(t.get('runs_completed', t.get('tasks_completed', [])))} | "
              f"remaining: {len(t.get('runs_remaining', t.get('tasks_remaining', [])))}")

    # Collect all MD summaries
    all_md = []
    missing = []
    for peptide in PEPTIDES:
        for target in TARGETS:
            summary = load_md_summary(peptide, target)
            if summary:
                all_md.append(summary)
            else:
                missing.append(f"{peptide}_{target}")

    print(f"\n=== MD Results: {len(all_md)}/12 available ===")
    if missing:
        print(f"  Missing: {', '.join(missing)}")

    # Write all_rmsd.csv
    if all_md:
        rmsd_path = COMPILED / "all_rmsd.csv"
        with open(rmsd_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "peptide", "target", "terminal", "rmsd_final", "rmsd_mean",
                "rmsd_max", "com_drift_final", "contact_fraction_final", "verdict"
            ])
            writer.writeheader()
            for md in sorted(all_md, key=lambda x: (x.get("peptide", ""), x.get("target", ""))):
                writer.writerow({k: md.get(k, "") for k in writer.fieldnames})
        print(f"  Wrote: {rmsd_path.name}")

    # Build selectivity matrix
    matrix = {}
    for peptide in PEPTIDES:
        matrix[peptide] = {}
        for target in TARGETS:
            md = load_md_summary(peptide, target)
            if md:
                matrix[peptide][target] = {
                    "rmsd_final": md.get("rmsd_final"),
                    "verdict": md.get("verdict"),
                    "com_drift_final": md.get("com_drift_final"),
                }
            else:
                matrix[peptide][target] = {"status": "pending"}

    matrix_path = COMPILED / "selectivity_matrix.json"
    with open(matrix_path, "w") as f:
        json.dump(matrix, f, indent=2)
    print(f"  Wrote: {matrix_path.name}")

    # Gate result
    gate = check_gate()
    gate_path = COMPILED / "gate_result.json"
    with open(gate_path, "w") as f:
        json.dump(gate, f, indent=2)
    print(f"\n=== Gate Result: {gate['gate_verdict']} ===")

    # Load SASA and contact data
    sasa = load_sasa_results()
    contacts = load_contact_results()

    # Write final verdict
    verdict_path = COMPILED / "final_verdict.md"
    with open(verdict_path, "w") as f:
        f.write("# V4.3 Wet Lab Validation — Final Results\n\n")
        f.write(f"**Compiled:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n")
        f.write(f"**MD Results Available:** {len(all_md)}/12\n")
        f.write(f"**Gate Status:** {gate['gate_verdict']}\n\n")

        # SASA section
        f.write("## 1. SASA Zinc Capping Test (Terminal A)\n\n")
        if sasa:
            f.write("| Peptide | Occlusion % | Min Zn Distance | Verdict |\n")
            f.write("|---------|------------|-----------------|----------|\n")
            for r in sasa.get("results", []):
                f.write(f"| {r['peptide']} | {r['occlusion_pct']}% | {r['min_zinc_dist_A']} A | {r['verdict']} |\n")
            f.write(f"\n**Overall:** {sasa.get('overall_verdict', 'N/A')}\n\n")
        else:
            f.write("*SASA results not yet available.*\n\n")

        # Contact maps section
        f.write("## 2. Contact Map Analysis (Terminal A)\n\n")
        if contacts:
            f.write("| Peptide | Total Contacts | H-bonds | K392 | Floor | Wall | Ceiling | Verdict |\n")
            f.write("|---------|---------------|---------|------|-------|------|---------|----------|\n")
            for r in contacts.get("results", []):
                f.write(f"| {r['peptide']} | {r['total_contacts']} | {r['hbond_candidates']} | "
                        f"{r['k392_contacts']} | {r['floor']} | {r['wall']} | {r['ceiling']} | "
                        f"{r['region_verdict']} |\n")
            f.write("\n**Note:** All peptides classify as LOCAL_ANCHOR — they sit at the channel "
                    "entrance (residues 324-371) rather than deep near K392 (392-414). "
                    "This is still a valid inhibition mechanism (entrance plugging blocks substrate access). "
                    "MD will reveal if peptides migrate deeper during dynamics.\n\n")
        else:
            f.write("*Contact map results not yet available.*\n\n")

        # MD section
        f.write("## 3. MD Stability Results\n\n")
        if all_md:
            f.write("| Peptide | Target | RMSD Final | RMSD Mean | COM Drift | Verdict |\n")
            f.write("|---------|--------|-----------|-----------|-----------|----------|\n")
            for md in sorted(all_md, key=lambda x: (x.get("peptide", ""), x.get("target", ""))):
                f.write(f"| {md.get('peptide', '?')} | {md.get('target', '?')} | "
                        f"{md.get('rmsd_final', '?')} | {md.get('rmsd_mean', '?')} | "
                        f"{md.get('com_drift_final', '?')} | {md.get('verdict', '?')} |\n")
            f.write("\n")
        else:
            f.write("*MD results not yet available. Awaiting Terminals B, C, D.*\n\n")

        if missing:
            f.write(f"**Pending runs:** {', '.join(missing)}\n\n")

        # Gate section
        f.write("## 4. FASGAV Gate Experiment\n\n")
        if gate["gate_verdict"] == "PASS":
            f.write("**GATE PASSED** — Scrambled control drifted out of ERAP2-K392. "
                    "Lead results are trustworthy.\n\n")
        elif gate["gate_verdict"] == "FAIL":
            f.write("**GATE FAILED** — Scrambled control locked in ERAP2-K392. "
                    "Pocket is artificially sticky. ALL selectivity claims are suspect.\n\n")
        elif gate["gate_verdict"] == "AMBIGUOUS":
            f.write("**GATE AMBIGUOUS** — FASGAV RMSD between 2-3 A. Inconclusive. "
                    "May need longer simulation or additional controls.\n\n")
        else:
            f.write("**GATE PENDING** — Awaiting Terminal C FASGAV results.\n\n")

        # Synthesis recommendation
        f.write("## 5. Synthesis Recommendation\n\n")
        if len(all_md) == 12 and gate["gate_verdict"] == "PASS":
            locked_leads = [md for md in all_md
                           if md.get("target") == "erap2k392"
                           and md.get("verdict") == "LOCKED"
                           and md.get("peptide") != "FASGAV"]
            if locked_leads:
                f.write("**PROCEED WITH SYNTHESIS** for:\n")
                for lead in locked_leads:
                    f.write(f"- {lead['peptide']} (RMSD: {lead.get('rmsd_final', '?')} A)\n")
                f.write("\nSee `docs/CRO_SYNTHESIS_SPECS.md` for order details.\n")
                f.write("See `docs/PROVISIONAL_PATENT_APPLICATION.md` for IP filing.\n")
            else:
                f.write("**NO LEADS PASSED** — All peptides drifted out of ERAP2-K392. "
                        "Boltz-2 predictions did not translate to MD stability.\n")
        elif len(all_md) < 12:
            f.write("**WAITING** — MD simulations not yet complete.\n")
        elif gate["gate_verdict"] == "FAIL":
            f.write("**DO NOT SYNTHESIZE** — Gate experiment failed. Re-evaluate starting poses.\n")
        else:
            f.write("**REVIEW REQUIRED** — Check individual results before proceeding.\n")

    print(f"\n  Wrote: {verdict_path.name}")
    print(f"\nCompilation complete. {len(all_md)}/12 MD runs available.")


if __name__ == "__main__":
    main()
