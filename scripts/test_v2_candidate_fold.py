#!/usr/bin/env python3
"""Multi-algorithm fold assessment for v2 synthesis candidates.

Two modes:
  1. Manual entry — paste scores from free web servers into fold_scores.yaml
  2. Auto-parse  — reads foldism output directory (if Modal was used)

Web servers (free, paste sequence):
  - Chai Discovery:    https://lab.chaidiscovery.com
  - Boltz-2 Tamarind:  https://tamarind.bio/tools/boltz
  - Boltz-2 Neurosnap: https://neurosnap.ai/service/Boltz-2%20(AlphaFold3)
  - AlphaFold Server:  https://alphafoldserver.com
  - Boltz-2 Rowan:     https://rowansci.com/features/free-online-boltz-2

Usage:
  # 1. Initialize empty scores template
  python scripts/test_v2_candidate_fold.py --init

  # 2. Fill in fold_scores.yaml with results from web servers

  # 3. Generate report
  python scripts/test_v2_candidate_fold.py

  # Or auto-parse foldism output (if you ran Modal)
  python scripts/test_v2_candidate_fold.py --foldism-dir data/results/foldism/fold_candidates
"""

import argparse
import json
import sys
from pathlib import Path

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

SCORES_FILE = Path("v2/fold_scores.yaml")
CANDIDATES = [
    "n248_trim_c5_Y87A_Y89A",
    "n248_trim_c5",
    "n248_wt",
    "n248_ko_all_aromatics",
]
SEQUENCES = {
    "n248_trim_c5_Y87A_Y89A": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
    "n248_trim_c5":           "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
    "n248_wt":                "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
    "n248_ko_all_aromatics":  "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
}
ESMFOLD_PLDDT = {
    "n248_trim_c5_Y87A_Y89A": 62.0,
    "n248_trim_c5": 62.0,
    "n248_wt": 58.0,
    "n248_ko_all_aromatics": None,
}

TEMPLATE = {
    cand: {
        "sequence": SEQUENCES[cand],
        "length": len(SEQUENCES[cand]),
        "esmfold": {"plddt": ESMFOLD_PLDDT[cand]},
        "boltz2": {"plddt": None, "ptm": None, "confidence": None, "source": "tamarind / neurosnap / rowan"},
        "chai1": {"aggregate_score": None, "ptm": None, "source": "lab.chaidiscovery.com"},
        "alphafold3": {"plddt": None, "ptm": None, "source": "alphafoldserver.com"},
    }
    for cand in CANDIDATES
}


def init_scores():
    """Write template fold_scores.yaml for manual entry."""
    SCORES_FILE.parent.mkdir(parents=True, exist_ok=True)
    if HAS_YAML:
        with open(SCORES_FILE, "w") as f:
            yaml.dump(TEMPLATE, f, default_flow_style=False, sort_keys=False, width=120)
    else:
        # Fallback: write as JSON
        scores_json = SCORES_FILE.with_suffix(".json")
        with open(scores_json, "w") as f:
            json.dump(TEMPLATE, f, indent=2)
        print(f"PyYAML not installed, wrote JSON instead: {scores_json}")
        print("  pip install pyyaml  # then re-run --init for YAML format")
        return

    print(f"Wrote template: {SCORES_FILE}")
    print()
    print("Fill in scores from web servers:")
    print()
    for cand in CANDIDATES:
        print(f"  {cand} ({len(SEQUENCES[cand])}aa):")
        print(f"    {SEQUENCES[cand]}")
        print()
    print("Web servers to use:")
    print("  1. Chai Discovery    — https://lab.chaidiscovery.com")
    print("  2. Boltz-2 (Tamarind)— https://tamarind.bio/tools/boltz")
    print("  3. AlphaFold Server  — https://alphafoldserver.com")
    print()
    print(f"Then re-run: python scripts/test_v2_candidate_fold.py")


def load_scores() -> dict:
    """Load scores from YAML or JSON."""
    if SCORES_FILE.exists() and HAS_YAML:
        with open(SCORES_FILE) as f:
            return yaml.safe_load(f)
    json_path = SCORES_FILE.with_suffix(".json")
    if json_path.exists():
        with open(json_path) as f:
            return json.load(f)
    return None


def parse_foldism_dir(foldism_dir: Path) -> dict:
    """Auto-parse foldism output into scores dict."""
    scores = {}
    for cand in CANDIDATES:
        scores[cand] = {
            "esmfold": {"plddt": ESMFOLD_PLDDT[cand]},
            "boltz2": {},
            "chai1": {},
            "alphafold3": {},
        }

    # Boltz-2
    for f in foldism_dir.glob("**/boltz2/**/confidence*.json"):
        with open(f) as fh:
            data = json.load(fh)
            # Apply to first candidate found (single-sequence run)
            for cand in CANDIDATES:
                if not scores[cand]["boltz2"]:
                    scores[cand]["boltz2"] = {
                        "plddt": data.get("plddt"),
                        "ptm": data.get("ptm"),
                        "confidence": data.get("confidence_score"),
                    }
                    break

    # Chai-1
    try:
        import numpy as np
        for f in foldism_dir.glob("**/chai1/**/*.npz"):
            data = np.load(str(f))
            for cand in CANDIDATES:
                if not scores[cand]["chai1"]:
                    scores[cand]["chai1"] = {
                        "aggregate_score": float(data.get("aggregate_score", [0])[0]),
                        "ptm": float(data["ptm"][0]) if "ptm" in data else None,
                    }
                    break
    except ImportError:
        pass

    return scores


def get_val(entry: dict, *keys) -> float | None:
    """Safely extract a numeric value from nested dict."""
    for key in keys:
        if isinstance(entry, dict) and key in entry:
            entry = entry[key]
        else:
            return None
    if entry is None:
        return None
    try:
        return float(entry)
    except (ValueError, TypeError):
        return None


def fmt(val: float | None, decimals: int = 1) -> str:
    if val is None:
        return "—"
    return f"{val:.{decimals}f}"


def assess_fold(plddt_vals: list[float]) -> str:
    """Assess fold confidence from available pLDDT/confidence values."""
    valid = [v for v in plddt_vals if v is not None]
    if not valid:
        return "NO DATA"
    avg = sum(valid) / len(valid)
    n = len(valid)
    if avg > 70:
        return f"FOLDS ({n} methods agree)"
    elif avg > 50:
        return f"PARTIAL ({n} methods, avg {avg:.0f})"
    else:
        return f"DISORDERED ({n} methods, avg {avg:.0f})"


def report(scores: dict):
    """Print fold assessment report."""
    print("=" * 90)
    print("  MULTI-ALGORITHM FOLD ASSESSMENT — V2 SYNTHESIS CANDIDATES")
    print("=" * 90)
    print()

    # Detailed per-candidate
    for cand in CANDIDATES:
        entry = scores.get(cand, {})
        seq_len = len(SEQUENCES[cand])
        print(f"  {cand} ({seq_len}aa)")
        print(f"  {'─' * 60}")

        esm = get_val(entry, "esmfold", "plddt")
        b2_plddt = get_val(entry, "boltz2", "plddt")
        b2_ptm = get_val(entry, "boltz2", "ptm")
        b2_conf = get_val(entry, "boltz2", "confidence")
        c1_agg = get_val(entry, "chai1", "aggregate_score")
        c1_ptm = get_val(entry, "chai1", "ptm")
        af3_plddt = get_val(entry, "alphafold3", "plddt")
        af3_ptm = get_val(entry, "alphafold3", "ptm")

        print(f"    ESMFold     pLDDT: {fmt(esm)}")
        print(f"    Boltz-2     pLDDT: {fmt(b2_plddt)}   pTM: {fmt(b2_ptm, 3)}   confidence: {fmt(b2_conf, 3)}")
        print(f"    Chai-1      aggr:  {fmt(c1_agg, 3)}   pTM: {fmt(c1_ptm, 3)}")
        print(f"    AlphaFold3  pLDDT: {fmt(af3_plddt)}   pTM: {fmt(af3_ptm, 3)}")

        # Collect all pLDDT-like values for consensus
        plddt_vals = [esm, b2_plddt, af3_plddt]
        # Chai aggregate_score is 0-1, scale to 0-100 for comparison
        if c1_agg is not None and c1_agg <= 1.0:
            plddt_vals.append(c1_agg * 100)
        elif c1_agg is not None:
            plddt_vals.append(c1_agg)

        verdict = assess_fold(plddt_vals)
        print(f"    >>> VERDICT: {verdict}")
        print()

    # Summary table
    print("─" * 90)
    print(f"  {'Candidate':<30} {'ESMFold':>8} {'Boltz2':>8} {'Chai1':>8} {'AF3':>8} {'Verdict':>20}")
    print("─" * 90)

    for cand in CANDIDATES:
        entry = scores.get(cand, {})
        esm = get_val(entry, "esmfold", "plddt")
        b2 = get_val(entry, "boltz2", "plddt")
        c1 = get_val(entry, "chai1", "aggregate_score")
        af3 = get_val(entry, "alphafold3", "plddt")

        plddt_vals = [esm, b2, af3]
        if c1 is not None and c1 <= 1.0:
            plddt_vals.append(c1 * 100)
        elif c1 is not None:
            plddt_vals.append(c1)

        verdict = assess_fold(plddt_vals)
        print(f"  {cand:<30} {fmt(esm):>8} {fmt(b2):>8} {fmt(c1, 3):>8} {fmt(af3):>8} {verdict:>20}")

    print("─" * 90)
    print()
    print("  SCALE:  pLDDT > 70 = confident fold  |  50-70 = partial  |  < 50 = disordered")
    print("  Chai-1 aggregate_score is 0-1 (scaled to 0-100 for comparison)")
    print()

    # Synthesis decision
    lead = scores.get("n248_trim_c5_Y87A_Y89A", {})
    neg = scores.get("n248_ko_all_aromatics", {})
    lead_b2 = get_val(lead, "boltz2", "plddt")
    neg_b2 = get_val(neg, "boltz2", "plddt")

    print("  SYNTHESIS DECISION:")
    if lead_b2 is not None and lead_b2 > 70:
        print("    ✓ Lead candidate folds confidently — PROCEED with synthesis")
    elif lead_b2 is not None and lead_b2 > 50:
        print("    ~ Lead candidate partially folds — may be conditionally structured")
        print("      (folds upon binding ERAP2, common for peptide binders)")
        print("      Consider: synthesis still viable, add circular dichroism validation")
    elif lead_b2 is not None:
        print("    ✗ Lead candidate likely disordered — RECONSIDER synthesis")
        print("      Options: add disulfide staple, cyclize, or redesign scaffold")
    else:
        print("    ? Insufficient data — fill in more scores and re-run")

    if lead_b2 is not None and neg_b2 is not None:
        if neg_b2 < lead_b2:
            print(f"    ✓ Negative control ({fmt(neg_b2)}) scores below lead ({fmt(lead_b2)}) — good")
        else:
            print(f"    ✗ WARNING: Negative control ({fmt(neg_b2)}) >= lead ({fmt(lead_b2)}) — investigate")
    print()


def main():
    parser = argparse.ArgumentParser(description="V2 candidate fold assessment")
    parser.add_argument("--init", action="store_true", help="Create template fold_scores.yaml")
    parser.add_argument("--foldism-dir", type=Path, help="Auto-parse foldism output directory")
    args = parser.parse_args()

    if args.init:
        init_scores()
        return

    if args.foldism_dir:
        scores = parse_foldism_dir(args.foldism_dir)
        report(scores)
        return

    scores = load_scores()
    if scores is None:
        print("No scores found. Initialize template first:")
        print("  python scripts/test_v2_candidate_fold.py --init")
        print()
        print("Then fill in v2/fold_scores.yaml with web server results.")
        sys.exit(1)

    report(scores)


if __name__ == "__main__":
    main()
