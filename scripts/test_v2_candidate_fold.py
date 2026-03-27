#!/usr/bin/env python3
"""Run foldism multi-algorithm fold prediction on v2 synthesis candidates.

Usage (run locally with Modal configured):
    # All 4 algorithms (Boltz-2, Chai-1, Protenix, AlphaFold2) — ~$1-2
    uvx modal run ../foldism/foldism.py --input-faa v2/fold_candidates.faa --out-dir data/results/foldism

    # Fast mode: skip MSA for quicker results (lower quality)
    uvx modal run ../foldism/foldism.py --input-faa v2/fold_candidates.faa --no-use-msa --out-dir data/results/foldism

    # Just Boltz-2 + Chai-1 (cheapest, ~$0.50)
    uvx modal run ../foldism/foldism.py --input-faa v2/fold_candidates.faa --algorithms boltz2,chai1 --out-dir data/results/foldism

After running, use this script to parse and summarize results:
    python scripts/test_v2_candidate_fold.py
"""

import json
import sys
from pathlib import Path

RESULTS_DIR = Path("data/results/foldism/fold_candidates")
CANDIDATES = [
    "n248_trim_c5_Y87A_Y89A",
    "n248_trim_c5",
    "n248_wt",
    "n248_ko_all_aromatics",
]
ALGORITHMS = ["boltz2", "chai1", "protenix", "alphafold2"]


def parse_boltz2_scores(result_dir: Path) -> dict:
    """Extract confidence scores from Boltz-2 output."""
    for f in result_dir.glob("**/confidence*.json"):
        with open(f) as fh:
            scores = json.load(fh)
            return {
                "confidence": scores.get("confidence_score"),
                "plddt": scores.get("plddt"),
                "ptm": scores.get("ptm"),
                "ipTM": scores.get("iptm"),
            }
    return {}


def parse_chai1_scores(result_dir: Path) -> dict:
    """Extract aggregate score from Chai-1 NPZ output."""
    try:
        import numpy as np
        for f in result_dir.glob("**/*.npz"):
            data = np.load(str(f))
            return {
                "aggregate_score": float(data.get("aggregate_score", [0])[0]),
                "ptm": float(data.get("ptm", [0])[0]) if "ptm" in data else None,
            }
    except ImportError:
        print("  (numpy not available, skipping Chai-1 NPZ parsing)")
    return {}


def parse_protenix_scores(result_dir: Path) -> dict:
    """Extract scores from Protenix output."""
    for f in result_dir.glob("**/*score*.json"):
        with open(f) as fh:
            return json.load(fh)
    for f in result_dir.glob("**/*.json"):
        with open(f) as fh:
            data = json.load(fh)
            if "plddt" in data or "confidence" in data or "ptm" in data:
                return data
    return {}


def parse_alphafold2_scores(result_dir: Path) -> dict:
    """Extract pLDDT/pTM from AlphaFold2 output."""
    for f in result_dir.glob("**/*ranking*.json"):
        with open(f) as fh:
            return json.load(fh)
    for f in result_dir.glob("**/*.json"):
        with open(f) as fh:
            data = json.load(fh)
            if "plddt" in data or "ptm" in data:
                return data
    return {}


PARSERS = {
    "boltz2": parse_boltz2_scores,
    "chai1": parse_chai1_scores,
    "protenix": parse_protenix_scores,
    "alphafold2": parse_alphafold2_scores,
}


def main():
    if not RESULTS_DIR.exists():
        print(f"Results directory not found: {RESULTS_DIR}")
        print()
        print("Run foldism first:")
        print("  cd /path/to/ancient-drug-discovery")
        print("  uvx modal run foldism/foldism.py --input-faa v2/fold_candidates.faa --out-dir data/results/foldism")
        sys.exit(1)

    print("=" * 80)
    print("FOLDISM MULTI-ALGORITHM FOLD ASSESSMENT — V2 CANDIDATES")
    print("=" * 80)

    # Discover what's available
    algo_dirs = [d for d in RESULTS_DIR.iterdir() if d.is_dir()]
    print(f"\nFound algorithm results: {[d.name for d in algo_dirs]}")
    print()

    # Build results table
    results = {}
    for algo_dir in sorted(algo_dirs):
        algo = algo_dir.name
        parser = PARSERS.get(algo)
        if not parser:
            print(f"  No parser for {algo}, listing files:")
            for f in algo_dir.rglob("*"):
                if f.is_file():
                    print(f"    {f.relative_to(RESULTS_DIR)}")
            continue

        scores = parser(algo_dir)
        results[algo] = scores
        print(f"  {algo}: {json.dumps(scores, indent=2, default=str)}")

    # Summary table
    print()
    print("-" * 80)
    print(f"{'Candidate':<30} {'ESMFold':>8} {'Boltz2':>8} {'Chai1':>8} {'AF2':>8} {'Protenix':>8}")
    print("-" * 80)

    esmfold_plddt = {
        "n248_trim_c5_Y87A_Y89A": "60-64",
        "n248_trim_c5": "~62",
        "n248_wt": "~58",
        "n248_ko_all_aromatics": "N/A",
    }

    for cand in CANDIDATES:
        esm = esmfold_plddt.get(cand, "?")
        boltz = results.get("boltz2", {}).get("plddt", "—")
        chai = results.get("chai1", {}).get("aggregate_score", "—")
        af2 = results.get("alphafold2", {}).get("plddt", "—")
        prot = results.get("protenix", {}).get("plddt", "—")

        if isinstance(boltz, float):
            boltz = f"{boltz:.1f}"
        if isinstance(chai, float):
            chai = f"{chai:.3f}"
        if isinstance(af2, float):
            af2 = f"{af2:.1f}"
        if isinstance(prot, float):
            prot = f"{prot:.1f}"

        print(f"{cand:<30} {esm:>8} {str(boltz):>8} {str(chai):>8} {str(af2):>8} {str(prot):>8}")

    print("-" * 80)
    print()
    print("INTERPRETATION:")
    print("  pLDDT > 70: Confident fold  |  50-70: Partial/uncertain  |  < 50: Likely disordered")
    print("  If 3+/4 methods agree on fold → high confidence for synthesis")
    print("  If methods disagree → investigate; may be conditionally folded (folds on binding)")
    print()


if __name__ == "__main__":
    main()
