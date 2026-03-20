#!/usr/bin/env python
"""
08_bioreason_targets.py — Run all ancient-drug-discovery targets through BioReason-Pro GO-GPT.

Fetches protein sequences from UniProt, runs GO-GPT predictions, and saves
structured results per target to data/processed/<gene>/bioreason/.

Usage:
    python scripts/08_bioreason_targets.py                    # all targets
    python scripts/08_bioreason_targets.py --targets ERAP2 G6PD  # specific targets
    python scripts/08_bioreason_targets.py --device cpu       # force CPU
    python scripts/08_bioreason_targets.py --model wanglab/gogpt  # HuggingFace model

Outputs per target:
    data/processed/<gene>/bioreason/go_predictions.json
    data/processed/<gene>/bioreason/summary.txt
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

# Fix OpenMP duplicate library issue on Windows conda
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import requests
import yaml

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
BIOREASON_DIR = PROJECT_ROOT / "BioReason-Pro"
TARGETS_YAML = PROJECT_ROOT / "targets.yaml"
DATA_DIR = PROJECT_ROOT / "data" / "processed"

# Add BioReason-Pro to path
sys.path.insert(0, str(BIOREASON_DIR))
sys.path.insert(0, str(BIOREASON_DIR / "gogpt" / "src"))


def load_targets(yaml_path: Path, filter_targets: Optional[List[str]] = None) -> dict:
    """Load target definitions from targets.yaml."""
    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    targets = config["targets"]
    if filter_targets:
        targets = {k: v for k, v in targets.items() if k in filter_targets}

    return targets


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    # Skip header line (starts with >)
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence


def run_gogpt(predictor, sequence: str, organism: str = "Homo sapiens") -> dict:
    """Run GO-GPT prediction on a single sequence."""
    predictions = predictor.predict(sequence=sequence, organism=organism)
    return predictions


def format_go_terms(predictions: dict, go_info: Optional[dict] = None) -> str:
    """Format GO predictions as human-readable text."""
    aspect_names = {
        "MF": "Molecular Function",
        "BP": "Biological Process",
        "CC": "Cellular Component",
    }

    lines = []
    for aspect in ["MF", "BP", "CC"]:
        terms = predictions.get(aspect, [])
        lines.append(f"\n{aspect_names[aspect]} ({aspect}):")
        if not terms:
            lines.append("  (none predicted)")
        else:
            for go_id in terms:
                if go_info:
                    name, _ = go_info.get(go_id, ("Unknown", ""))
                    lines.append(f"  {go_id} — {name}")
                else:
                    lines.append(f"  {go_id}")

    return "\n".join(lines)


def build_summary(gene: str, target: dict, predictions: dict, sequence: str) -> str:
    """Build a full summary report for a target."""
    lines = [
        f"{'='*70}",
        f"BioReason-Pro GO-GPT Analysis: {gene}",
        f"{'='*70}",
        f"",
        f"Protein: {target['protein_name']}",
        f"UniProt: {target['uniprot_id']}",
        f"Sequence length: {len(sequence)} aa",
        f"Ancient pathogen: {target['ancient_pathogen']}",
        f"Immune step: {target['immune_step']}",
        f"Drug status: {target['drug_status']}",
        f"",
        f"--- GO Term Predictions ---",
        format_go_terms(predictions),
        f"",
        f"--- Cancer Context ---",
        f"{target['cancer_mechanism']}",
        f"Cancer types: {', '.join(target.get('cancer_types', []))}",
        f"",
        f"--- Pipeline Notes ---",
        f"These GO predictions can be cross-referenced with:",
        f"  - Script 07 cancer evidence matrix (OpenTargets)",
        f"  - Script 05 ESM-2 variant effect scores",
        f"  - RFdiffusion binding site selection",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Run ancient-drug-discovery targets through BioReason-Pro GO-GPT."
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help="Target gene symbols to run (default: all from targets.yaml)",
    )
    parser.add_argument(
        "--model",
        default="wanglab/gogpt",
        help="HuggingFace model ID (default: wanglab/gogpt)",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Device: 'cuda', 'cpu', or None for auto-detect",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Cache directory for model weights",
    )
    args = parser.parse_args()

    # Load targets
    targets = load_targets(TARGETS_YAML, args.targets)
    print(f"Loaded {len(targets)} targets: {', '.join(targets.keys())}")

    # Load GO-GPT model
    print(f"\nLoading GO-GPT model: {args.model}")
    from gogpt_api import load_predictor

    predictor = load_predictor(args.model, cache_dir=args.cache_dir)

    # Try to load GO term name lookup
    go_info = None
    try:
        from bioreason2.dataset.cafa5.processor import _GO_INFO
        go_info = _GO_INFO
    except ImportError:
        print("(GO term names not available — will show IDs only)")

    # Process each target
    results_all = {}

    for gene, target in targets.items():
        print(f"\n{'─'*50}")
        print(f"Processing: {gene} ({target['protein_name']})")
        print(f"{'─'*50}")

        # Fetch sequence from UniProt
        uniprot_id = target["uniprot_id"]
        print(f"  Fetching sequence from UniProt ({uniprot_id})...", end=" ", flush=True)
        try:
            sequence = fetch_uniprot_sequence(uniprot_id)
            print(f"{len(sequence)} aa")
        except Exception as e:
            print(f"FAILED: {e}")
            continue

        # Run GO-GPT
        print(f"  Running GO-GPT prediction...", end=" ", flush=True)
        t0 = time.time()
        try:
            predictions = run_gogpt(predictor, sequence)
            elapsed = time.time() - t0
            total_terms = sum(len(v) for v in predictions.values())
            print(f"{total_terms} GO terms ({elapsed:.1f}s)")
        except Exception as e:
            print(f"FAILED: {e}")
            continue

        # Print predictions
        print(format_go_terms(predictions, go_info))

        # Save results
        out_dir = DATA_DIR / gene / "bioreason"
        out_dir.mkdir(parents=True, exist_ok=True)

        # JSON with full predictions + metadata
        result = {
            "gene": gene,
            "uniprot_id": uniprot_id,
            "protein_name": target["protein_name"],
            "sequence_length": len(sequence),
            "model": args.model,
            "predictions": predictions,
            "ancient_pathogen": target["ancient_pathogen"],
            "immune_step": target["immune_step"],
            "cancer_types": target.get("cancer_types", []),
        }

        json_path = out_dir / "go_predictions.json"
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  Saved: {json_path}")

        # Human-readable summary
        summary = build_summary(gene, target, predictions, sequence)
        summary_path = out_dir / "summary.txt"
        with open(summary_path, "w") as f:
            f.write(summary)
        print(f"  Saved: {summary_path}")

        results_all[gene] = result

    # Final summary
    print(f"\n{'='*50}")
    print(f"COMPLETE: {len(results_all)}/{len(targets)} targets processed")
    print(f"{'='*50}")

    for gene, result in results_all.items():
        preds = result["predictions"]
        counts = {k: len(v) for k, v in preds.items()}
        print(f"  {gene}: MF={counts['MF']} BP={counts['BP']} CC={counts['CC']}")

    print(f"\nResults saved to: {DATA_DIR}/<gene>/bioreason/")


if __name__ == "__main__":
    main()
