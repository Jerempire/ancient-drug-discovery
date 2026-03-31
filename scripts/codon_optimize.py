"""
Codon Optimization Module for Ancient Drug Discovery Pipeline.

Integrates CodonTransformer (BigBird, 164 organisms, 1M+ DNA-protein pairs)
into the pipeline as Stage 6b: peptide → codon-optimized mRNA.

Pipeline position:
  Stage 4 (Boltz-2 binder design) → Stage 5 (validation) →
  Stage 6a (synthesis planning) → **Stage 6b (codon optimization)** →
  Stage 7 (mRNA construct design)

Usage:
  # Single peptide
  python codon_optimize.py --peptide VAGSAF --organism "Homo sapiens"

  # From candidates.db (all validated peptides)
  python codon_optimize.py --from-db --organism "Homo sapiens"

  # Multi-sequence diversity (for library screening)
  python codon_optimize.py --peptide VAGSAF --organism "Homo sapiens" --n 5 --temperature 0.4

  # Evaluate existing DNA sequence
  python codon_optimize.py --evaluate GTGGCCGGCTCTGCCTTC --organism "Homo sapiens"

Requires: pip install CodonTransformer
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass  # Jupyter — already UTF-8

import json
import sqlite3
import argparse
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional

# ── Config ──────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "data" / "results" / "codon_optimized"
DB_PATH = PROJECT_ROOT / "candidates.db"

# Default organism for human therapeutic peptides
DEFAULT_ORGANISM = "Homo sapiens"

# Evaluation thresholds
GC_CONTENT_MIN = 30.0   # Too low = poor stability
GC_CONTENT_MAX = 70.0   # Too high = secondary structures
RARE_CODON_THRESHOLD = 0.10  # Flag if >10% rare codons


@dataclass
class CodonResult:
    """Result of codon optimization for a single peptide."""
    peptide: str
    organism: str
    dna_sequence: str
    gc_content: float
    csi_score: float           # Codon Similarity Index (like CAI)
    rare_codon_fraction: float
    sequence_complexity: float
    n_codons: int
    warnings: list


def _load_model(device_str: str = "cpu"):
    """Load CodonTransformer model and tokenizer."""
    import torch
    from CodonTransformer.CodonPrediction import load_model, load_tokenizer

    device = torch.device(device_str)
    tokenizer = load_tokenizer("adibvafa/CodonTransformer")
    model = load_model(model_path="adibvafa/CodonTransformer", device=device)
    return model, tokenizer, device


def optimize_peptide(
    peptide: str,
    organism: str = DEFAULT_ORGANISM,
    model=None,
    tokenizer=None,
    device=None,
    deterministic: bool = True,
    temperature: float = 0.2,
    top_p: float = 0.95,
    num_sequences: int = 1,
) -> list[CodonResult]:
    """
    Optimize codon usage for a peptide sequence.

    Args:
        peptide: Amino acid sequence (e.g., "VAGSAF")
        organism: Target organism name (164 supported)
        model/tokenizer/device: Pre-loaded model (optional, loaded if None)
        deterministic: If True, pick highest-probability codons
        temperature: Sampling temperature (lower = more conservative)
        top_p: Nucleus sampling threshold
        num_sequences: Number of diverse sequences to generate

    Returns:
        List of CodonResult objects
    """
    from CodonTransformer.CodonPrediction import predict_dna_sequence
    from CodonTransformer.CodonEvaluation import (
        get_GC_content,
        get_CSI_weights,
        get_CSI_value,
        get_cfd,
        get_sequence_complexity,
    )
    from CodonTransformer.CodonUtils import ORGANISM2ID
    import torch

    # Load model if not provided
    if model is None or tokenizer is None or device is None:
        model, tokenizer, device = _load_model(
            "cuda" if torch.cuda.is_available() else "cpu"
        )

    # Generate DNA sequence(s)
    if num_sequences > 1:
        deterministic = False

    predictions = predict_dna_sequence(
        protein=peptide,
        organism=organism,
        device=device,
        tokenizer=tokenizer,
        model=model,
        deterministic=deterministic,
        temperature=temperature,
        top_p=top_p,
        num_sequences=num_sequences,
        match_protein=True,  # Ensure valid codons for each amino acid
    )

    if not isinstance(predictions, list):
        predictions = [predictions]

    # Evaluate each prediction
    results = []
    for pred in predictions:
        dna = pred.predicted_dna
        warnings = []

        # GC content
        gc = get_GC_content(dna)
        if gc < GC_CONTENT_MIN:
            warnings.append(f"Low GC content ({gc:.1f}%) — may have poor mRNA stability")
        if gc > GC_CONTENT_MAX:
            warnings.append(f"High GC content ({gc:.1f}%) — risk of secondary structures")

        # Codon Similarity Index (organism-specific)
        try:
            csi_weights = get_CSI_weights([dna])
            csi = get_CSI_value(dna, csi_weights)
        except Exception:
            csi = -1.0
            warnings.append("CSI calculation failed — may need reference sequences")

        # Rare codon fraction
        try:
            # Build codon frequency table for the organism
            from CodonTransformer.CodonUtils import get_codon_frequencies
            codon_freq = get_codon_frequencies(organism)
            rare_frac = get_cfd(dna, codon_freq, threshold=0.1)
        except Exception:
            rare_frac = -1.0

        if rare_frac > RARE_CODON_THRESHOLD:
            warnings.append(f"High rare codon fraction ({rare_frac:.2f}) — may reduce expression")

        # Sequence complexity
        complexity = get_sequence_complexity(dna)

        results.append(CodonResult(
            peptide=peptide,
            organism=organism,
            dna_sequence=dna,
            gc_content=round(gc, 2),
            csi_score=round(csi, 4) if csi >= 0 else csi,
            rare_codon_fraction=round(rare_frac, 4) if rare_frac >= 0 else rare_frac,
            sequence_complexity=round(complexity, 4),
            n_codons=len(peptide),
            warnings=warnings,
        ))

    return results


def optimize_from_db(
    organism: str = DEFAULT_ORGANISM,
    min_iptm: float = 0.7,
) -> list[CodonResult]:
    """
    Optimize all validated peptides from candidates.db.

    Reads peptides with ipTM >= min_iptm from the database,
    runs codon optimization, and saves results.
    """
    import torch

    if not DB_PATH.exists():
        print(f"candidates.db not found at {DB_PATH}")
        return []

    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row

    # Try common table/column names for peptide candidates
    try:
        rows = conn.execute(
            "SELECT DISTINCT peptide_sequence FROM candidates WHERE iptm >= ? ORDER BY iptm DESC",
            (min_iptm,)
        ).fetchall()
        peptides = [r['peptide_sequence'] for r in rows]
    except sqlite3.OperationalError:
        # Fallback: try other common schemas
        try:
            rows = conn.execute(
                "SELECT DISTINCT sequence FROM peptide_candidates WHERE score >= ? ORDER BY score DESC",
                (min_iptm,)
            ).fetchall()
            peptides = [r['sequence'] for r in rows]
        except sqlite3.OperationalError:
            print("Could not find peptide table in candidates.db. Expected 'candidates' or 'peptide_candidates'.")
            conn.close()
            return []

    conn.close()

    if not peptides:
        print(f"No peptides found with ipTM >= {min_iptm}")
        return []

    print(f"Optimizing {len(peptides)} peptides for {organism}...")

    # Load model once, reuse for all peptides
    model, tokenizer, device = _load_model(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    all_results = []
    for pep in peptides:
        results = optimize_peptide(
            pep, organism=organism,
            model=model, tokenizer=tokenizer, device=device,
        )
        all_results.extend(results)
        print(f"  {pep} → {results[0].dna_sequence} (GC={results[0].gc_content}%)")

    return all_results


def save_results(results: list[CodonResult], output_dir: Optional[Path] = None):
    """Save codon optimization results to JSON and summary table."""
    out = output_dir or RESULTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    # JSON (full detail)
    json_path = out / "codon_optimization_results.json"
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump([asdict(r) for r in results], f, indent=2, ensure_ascii=False)
    print(f"Results saved to {json_path}")

    # Summary table (quick reference)
    summary_path = out / "codon_optimization_summary.md"
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("# Codon Optimization Results\n\n")
        f.write(f"Organism: {results[0].organism}\n")
        f.write(f"Peptides: {len(results)}\n\n")
        f.write("| Peptide | DNA Sequence | GC% | CSI | Rare Codons | Warnings |\n")
        f.write("|---------|-------------|-----|-----|-------------|----------|\n")
        for r in results:
            warns = "; ".join(r.warnings) if r.warnings else "None"
            f.write(f"| {r.peptide} | {r.dna_sequence} | {r.gc_content} | {r.csi_score} | {r.rare_codon_fraction} | {warns} |\n")

    print(f"Summary saved to {summary_path}")


def main():
    parser = argparse.ArgumentParser(description="Codon optimization for peptide candidates")
    parser.add_argument("--peptide", type=str, help="Amino acid sequence to optimize")
    parser.add_argument("--organism", type=str, default=DEFAULT_ORGANISM, help="Target organism")
    parser.add_argument("--from-db", action="store_true", help="Optimize all validated peptides from candidates.db")
    parser.add_argument("--min-iptm", type=float, default=0.7, help="Minimum ipTM for --from-db")
    parser.add_argument("--n", type=int, default=1, help="Number of diverse sequences to generate")
    parser.add_argument("--temperature", type=float, default=0.2, help="Sampling temperature")
    parser.add_argument("--evaluate", type=str, help="Evaluate an existing DNA sequence")
    parser.add_argument("--output", type=str, help="Output directory (default: data/results/codon_optimized/)")

    args = parser.parse_args()
    output_dir = Path(args.output) if args.output else None

    if args.evaluate:
        from CodonTransformer.CodonEvaluation import get_GC_content, get_sequence_complexity
        dna = args.evaluate.upper()
        print(f"DNA: {dna}")
        print(f"Length: {len(dna)} bp ({len(dna)//3} codons)")
        print(f"GC content: {get_GC_content(dna):.1f}%")
        print(f"Sequence complexity: {get_sequence_complexity(dna):.4f}")
        return

    if args.from_db:
        results = optimize_from_db(organism=args.organism, min_iptm=args.min_iptm)
    elif args.peptide:
        results = optimize_peptide(
            args.peptide,
            organism=args.organism,
            num_sequences=args.n,
            temperature=args.temperature,
        )
    else:
        parser.print_help()
        return

    if results:
        for r in results:
            print(f"\n{'='*60}")
            print(f"Peptide:    {r.peptide}")
            print(f"DNA:        {r.dna_sequence}")
            print(f"Organism:   {r.organism}")
            print(f"GC content: {r.gc_content}%")
            print(f"CSI score:  {r.csi_score}")
            print(f"Rare codons:{r.rare_codon_fraction}")
            print(f"Complexity: {r.sequence_complexity}")
            if r.warnings:
                print(f"⚠ Warnings: {'; '.join(r.warnings)}")

        save_results(results, output_dir)


if __name__ == "__main__":
    main()
