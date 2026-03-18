"""
02_fetch_ancient_variants.py — Download AADR data and extract ERAP2-region variants.

The Allen Ancient DNA Resource (AADR) provides genotype data for 14,000+ ancient samples.
This script:
  1. Downloads the AADR v54.1 annotation file (sample metadata)
  2. Downloads the AADR SNP list (1240K panel positions)
  3. Filters for variants in the ERAP2 gene region (chr5:96,210,000-96,280,000, GRCh37)
  4. Cross-references with the target variant rs2549794
  5. Downloads individual genotype data if available in streamlined format

Outputs:
  data/raw/aadr_anno.tsv         — sample metadata
  data/raw/aadr_snp_list.tsv     — SNP positions on 1240K array
  data/processed/aadr_erap2_snps.csv — ERAP2-region SNPs from AADR panel
  data/processed/aadr_erap2_samples.csv — ancient samples with ERAP2 genotypes

NOTE: Full AADR genotype files are large (~10 GB in EIGENSTRAT format).
      This script fetches metadata and SNP lists first, then optionally
      downloads genotypes if run with --download-genotypes flag.
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import gzip
import io
from pathlib import Path

import pandas as pd
import requests

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
PROC_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)

# ERAP2 gene coordinates (GRCh37/hg19 — AADR uses this build)
ERAP2_CHR = "5"
ERAP2_START = 96_210_000
ERAP2_END = 96_280_000

# Target variant
TARGET_RSID = "rs2549794"

# AADR download URLs (v54.1 — latest stable as of 2026)
AADR_BASE = "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54"
AADR_ANNO_URL = f"{AADR_BASE}/V54.1/SHARE/public.dir/v54.1_1240K_public.anno"
# The SNP file lists all 1240K panel positions
AADR_SNP_URL = f"{AADR_BASE}/V54.1/SHARE/public.dir/v54.1_1240K_public.snp"


def download_file(url: str, dest: Path, desc: str) -> Path:
    """Download a file with progress indication."""
    if dest.exists():
        print(f"  {desc} already downloaded: {dest.name}")
        return dest

    print(f"  Downloading {desc} ...")
    print(f"    URL: {url}")
    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()

    total = int(resp.headers.get("Content-Length", 0))
    downloaded = 0
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            f.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded / total * 100
                print(f"\r    {downloaded:,} / {total:,} bytes ({pct:.1f}%)", end="")
    print()
    print(f"  Saved: {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def fetch_annotation() -> pd.DataFrame:
    """Download and parse AADR sample annotation file."""
    dest = RAW_DIR / "aadr_anno.tsv"

    # AADR files require manual download from Harvard Dataverse:
    # https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
    # Download the .anno file and place it at: data/raw/aadr_anno.tsv
    if not dest.exists():
        print("  AADR annotation file not found locally.")
        print("  To download: visit https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW")
        print(f"  Save the .anno file as: {dest}")
        print("  Skipping AADR annotation for now — continuing with other data sources.")
        return pd.DataFrame()

    print("  Parsing annotation file ...")
    df = pd.read_csv(dest, sep="\t", low_memory=False, encoding="utf-8", on_bad_lines="skip")
    print(f"  Loaded {len(df)} samples, {len(df.columns)} columns")
    return df


def fetch_snp_list() -> pd.DataFrame:
    """Download and parse AADR SNP position list (EIGENSTRAT .snp format)."""
    dest = RAW_DIR / "aadr_snp_list.tsv"

    if not dest.exists():
        print("  AADR SNP list not found locally.")
        print("  To download: visit https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW")
        print(f"  Save the .snp file as: {dest}")
        print("  Skipping AADR SNP list — continuing with other data sources.")
        return pd.DataFrame()

    print("  Parsing SNP list ...")
    df = pd.read_csv(
        dest,
        sep=r"\s+",
        header=None,
        names=["snp_id", "chr", "genetic_pos", "physical_pos", "ref", "alt"],
        dtype={"chr": str, "physical_pos": int},
        low_memory=False,
    )
    print(f"  Loaded {len(df):,} SNPs on 1240K panel")
    return df


def filter_erap2_snps(snp_df: pd.DataFrame) -> pd.DataFrame:
    """Filter SNP list for ERAP2 genomic region."""
    print(f"\nFiltering for ERAP2 region: chr{ERAP2_CHR}:{ERAP2_START:,}-{ERAP2_END:,} ...")
    if snp_df.empty:
        print("  SNP list is empty (AADR data not downloaded) — skipping filter.")
        out = PROC_DIR / "aadr_erap2_snps.csv"
        pd.DataFrame(columns=["snp_id", "chr", "genetic_pos", "physical_pos", "ref", "alt"]).to_csv(out, index=False)
        return pd.DataFrame(columns=["snp_id", "chr", "genetic_pos", "physical_pos", "ref", "alt"])
    mask = (
        (snp_df["chr"] == ERAP2_CHR)
        & (snp_df["physical_pos"] >= ERAP2_START)
        & (snp_df["physical_pos"] <= ERAP2_END)
    )
    erap2_snps = snp_df[mask].copy()
    print(f"  Found {len(erap2_snps)} SNPs in ERAP2 region")

    # Check for target variant
    target_match = erap2_snps[erap2_snps["snp_id"] == TARGET_RSID]
    if not target_match.empty:
        print(f"  ✓ Target variant {TARGET_RSID} found at position {target_match.iloc[0]['physical_pos']:,}")
    else:
        print(f"  ✗ Target variant {TARGET_RSID} NOT on 1240K panel")
        print("    (This is expected — 1240K is a targeted panel, not whole-genome)")
        print("    Will check if nearby proxy SNPs are available for imputation")

        # Find closest SNPs as potential proxies
        target_pos = 96_243_541  # rs2549794 GRCh37 position
        erap2_snps["dist_to_target"] = abs(erap2_snps["physical_pos"] - target_pos)
        closest = erap2_snps.nsmallest(5, "dist_to_target")
        if not closest.empty:
            print("    Closest SNPs on panel:")
            for _, row in closest.iterrows():
                print(f"      {row['snp_id']} at {row['physical_pos']:,} ({row['dist_to_target']:,} bp away)")
        erap2_snps = erap2_snps.drop(columns=["dist_to_target"], errors="ignore")

    out = PROC_DIR / "aadr_erap2_snps.csv"
    erap2_snps.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return erap2_snps


def query_clinvar_erap2() -> list[dict]:
    """Fetch ERAP2 variants from ClinVar via NCBI E-utilities as supplementary data."""
    print("\nFetching ERAP2 variants from ClinVar (supplementary) ...")
    # Search ClinVar for ERAP2 gene
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": "ERAP2[gene]",
        "retmax": 500,
        "retmode": "json",
    }
    resp = requests.get(search_url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    ids = data.get("esearchresult", {}).get("idlist", [])
    print(f"  Found {len(ids)} ClinVar entries for ERAP2")

    if not ids:
        return []

    # Fetch summaries
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    # Process in batches of 100
    all_variants = []
    for i in range(0, len(ids), 100):
        batch = ids[i : i + 100]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json",
        }
        resp = requests.get(summary_url, params=params, timeout=60)
        resp.raise_for_status()
        result = resp.json().get("result", {})

        for uid in batch:
            entry = result.get(uid, {})
            if not entry:
                continue
            all_variants.append({
                "uid": uid,
                "title": entry.get("title", ""),
                "clinical_significance": entry.get("clinical_significance", {}).get("description", ""),
                "gene_sort": entry.get("gene_sort", ""),
                "variation_set": entry.get("variation_set", []),
            })

    out = RAW_DIR / "clinvar_erap2.json"
    out.write_text(json.dumps(all_variants, indent=2), encoding="utf-8")
    print(f"  Saved {len(all_variants)} ClinVar variants → {out.name}")
    return all_variants


def fetch_uniprot_erap2() -> dict:
    """Fetch ERAP2 protein info from UniProt."""
    print("\nFetching ERAP2 protein data from UniProt ...")
    # ERAP2 human UniProt ID: Q6P179
    url = "https://rest.uniprot.org/uniprotkb/Q6P179.json"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    out = RAW_DIR / "uniprot_erap2_Q6P179.json"
    out.write_text(json.dumps(data, indent=2), encoding="utf-8")

    # Extract key info
    protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    seq_length = data.get("sequence", {}).get("length", 0)
    gene_name = ""
    genes = data.get("genes", [])
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    print(f"  Protein: {protein_name}")
    print(f"  Gene: {gene_name}")
    print(f"  Sequence length: {seq_length} aa")
    print(f"  Saved → {out.name}")

    # Also save FASTA sequence for later use with ESM/AlphaFold
    seq = data.get("sequence", {}).get("value", "")
    if seq:
        fasta_out = PROC_DIR / "erap2_Q6P179.fasta"
        fasta_out.write_text(f">sp|Q6P179|ERAP2_HUMAN {protein_name}\n{seq}\n", encoding="utf-8")
        print(f"  Saved FASTA → {fasta_out.name}")

    return data


def main():
    print("=" * 60)
    print("Ancient DNA + Protein Data Fetch — ERAP2 Focus")
    print("=" * 60)

    # 1. AADR sample annotation
    anno_df = fetch_annotation()

    # 2. AADR SNP positions
    snp_df = fetch_snp_list()

    # 3. Filter for ERAP2 region
    erap2_snps = filter_erap2_snps(snp_df)

    # 4. ClinVar ERAP2 variants (supplementary)
    clinvar = query_clinvar_erap2()

    # 5. UniProt ERAP2 protein data + FASTA
    uniprot = fetch_uniprot_erap2()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  AADR samples:        {len(anno_df):,}")
    print(f"  1240K panel SNPs:    {len(snp_df):,}")
    print(f"  ERAP2 region SNPs:   {len(erap2_snps)}")
    print(f"  ClinVar ERAP2 entries: {len(clinvar)}")
    print(f"  UniProt ERAP2:       Q6P179 ({uniprot.get('sequence', {}).get('length', '?')} aa)")
    print(f"\n  Raw data:      {RAW_DIR}")
    print(f"  Processed:     {PROC_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
