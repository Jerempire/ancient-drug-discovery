#!/usr/bin/env python3
"""
Foldseek structure similarity search wrapper.
Extracted from ScienceClaw (Apache 2.0).

Requires foldseek binary:
    conda install -c conda-forge -c bioconda foldseek

Or download from: https://mmseqs.com/foldseek/
"""

import subprocess
import sys
from pathlib import Path
from typing import Optional

try:
    import pandas as pd
except ImportError:
    print("Error: pandas is required. Install with: pip install pandas")
    sys.exit(1)


def foldseek_search(
    query_pdb: str,
    database: str,
    tmp_dir: str = "tmp/",
    e_value: float = 1e-3,
    max_hits: int = 100,
    exhaustive: bool = False,
    output_file: Optional[str] = None,
) -> "pd.DataFrame":
    """
    Run Foldseek structure similarity search.

    Args:
        query_pdb: Path to query PDB file
        database: Path to Foldseek database (e.g., from `foldseek databases PDB pdb_db tmp/`)
        tmp_dir: Temporary directory for Foldseek
        e_value: E-value threshold
        max_hits: Maximum number of hits
        exhaustive: Use exhaustive search (slower, more sensitive)
        output_file: Path for results TSV (default: auto-generated)

    Returns:
        DataFrame with columns: query, target, pident, alnlen, evalue, bits, prob, alntmscore, taxname
    """
    if output_file is None:
        output_file = str(Path(query_pdb).stem) + "_foldseek_results.tsv"

    cols = "query,target,pident,alnlen,evalue,bits,prob,alntmscore,taxname"

    cmd = [
        "foldseek", "easy-search",
        query_pdb, database, output_file, tmp_dir,
        "-e", str(e_value),
        "--max-seqs", str(max_hits),
        "--format-output", cols,
    ]

    if exhaustive:
        cmd.extend(["--exhaustive-search", "1"])

    print(f"Running Foldseek: {Path(query_pdb).name} vs {Path(database).name}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Foldseek failed: {result.stderr}")

    df = pd.read_csv(output_file, sep="\t", names=cols.split(","))
    df = df.sort_values("alntmscore", ascending=False).reset_index(drop=True)

    print(f"Found {len(df)} structural hits")
    return df


def check_design_novelty(query_pdb: str, database: str, tm_threshold: float = 0.5) -> dict:
    """
    Check if a designed protein has structural homologs in PDB.

    Returns:
        dict with 'novel' (bool), 'top_hits' (DataFrame), 'num_homologs' (int)
    """
    results = foldseek_search(query_pdb, database)
    homologs = results[results["alntmscore"] > tm_threshold]

    return {
        "novel": len(homologs) == 0,
        "num_homologs": len(homologs),
        "top_hits": results.head(10),
    }


def batch_search(query_dir: str, database: str, tmp_dir: str = "tmp/") -> "pd.DataFrame":
    """
    Search all PDB files in a directory against a database.

    Args:
        query_dir: Directory containing PDB files
        database: Foldseek database path

    Returns:
        Combined DataFrame of all results
    """
    # Create query database from directory
    query_db = str(Path(tmp_dir) / "query_db")
    result_db = str(Path(tmp_dir) / "result_db")
    output_file = str(Path(tmp_dir) / "batch_results.tsv")

    cols = "query,target,pident,alnlen,evalue,bits,prob,alntmscore,taxname"

    subprocess.run(["foldseek", "createdb", query_dir, query_db], check=True)
    subprocess.run([
        "foldseek", "search", query_db, database, result_db, tmp_dir,
        "-e", "1e-3"
    ], check=True)
    subprocess.run([
        "foldseek", "convertalis", query_db, database, result_db, output_file,
        "--format-output", cols
    ], check=True)

    df = pd.read_csv(output_file, sep="\t", names=cols.split(","))
    return df.sort_values(["query", "alntmscore"], ascending=[True, False]).reset_index(drop=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Foldseek structure similarity search")
    parser.add_argument("query", help="Query PDB file or directory")
    parser.add_argument("database", help="Foldseek database path")
    parser.add_argument("--max-hits", type=int, default=100)
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--exhaustive", action="store_true")
    parser.add_argument("--novelty-check", action="store_true", help="Check design novelty")
    args = parser.parse_args()

    if args.novelty_check:
        result = check_design_novelty(args.query, args.database)
        if result["novel"]:
            print("NOVEL — no structural homologs found in PDB")
        else:
            print(f"Found {result['num_homologs']} structural homologs (TM-score > 0.5)")
        print(result["top_hits"].to_string())
    else:
        df = foldseek_search(
            args.query, args.database,
            max_hits=args.max_hits,
            e_value=args.evalue,
            exhaustive=args.exhaustive,
        )
        print(df.to_string())
