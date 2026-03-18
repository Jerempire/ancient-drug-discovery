"""
01_fetch_gwas.py — Fetch ERAP2 / Crohn's disease associations from GWAS Catalog REST API.

Queries:
  1. All GWAS associations for ERAP2 gene
  2. All GWAS associations for Crohn's disease (EFO_0000384)
  3. Specific variant rs2549794 (the Black Death → Crohn's link)

Outputs:
  data/raw/gwas_erap2_associations.json
  data/raw/gwas_crohns_associations.json
  data/raw/gwas_rs2549794.json
  data/processed/erap2_crohns_overlap.csv
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import time
from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
PROC_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)

HEADERS = {"Accept": "application/json"}

# Rate-limit helper
def _get(url: str, params: dict | None = None) -> dict:
    """GET with retry on 429."""
    for attempt in range(3):
        resp = requests.get(url, headers=HEADERS, params=params, timeout=60)
        if resp.status_code == 429:
            wait = int(resp.headers.get("Retry-After", 5))
            print(f"  Rate-limited, waiting {wait}s ...")
            time.sleep(wait)
            continue
        resp.raise_for_status()
        return resp.json()
    raise RuntimeError(f"Failed after 3 attempts: {url}")


def _paginate_associations(url: str, params: dict | None = None, max_pages: int = 20) -> list[dict]:
    """Follow GWAS Catalog pagination links to collect all associations."""
    all_assocs = []
    page = 0
    while url and page < max_pages:
        data = _get(url, params if page == 0 else None)
        embedded = data.get("_embedded", {})
        assocs = embedded.get("associations", [])
        all_assocs.extend(assocs)
        # Next page link
        links = data.get("_links", {})
        next_link = links.get("next", {}).get("href")
        url = next_link
        params = None  # params are baked into next URL
        page += 1
        if next_link:
            print(f"  Page {page}, {len(all_assocs)} associations so far ...")
    return all_assocs


def fetch_erap2_associations() -> list[dict]:
    """ERAP2-related associations via known ERAP2 GWAS SNPs."""
    print("Fetching ERAP2 associations from GWAS Catalog ...")
    # GWAS Catalog has no findByGene endpoint — search by known ERAP2 rsIDs instead
    erap2_snps = [
        "rs2549794",   # Black Death → Crohn's link
        "rs2248374",   # ERAP2 splice variant
        "rs2910686",   # ERAP2 regulatory variant
        "rs2287987",   # ERAP2 missense
        "rs26653",     # ERAP2 missense
        "rs17408150",  # ERAP2 region
        "rs2927608",   # ERAP2 region
    ]
    all_assocs = []
    for rsid in erap2_snps:
        print(f"  Querying {rsid} ...")
        url = f"{BASE_URL}/associations/search/findByRsId"
        params = {"rsId": rsid}
        try:
            assocs = _paginate_associations(url, params, max_pages=5)
            for a in assocs:
                a["_query_rsid"] = rsid
            all_assocs.extend(assocs)
            print(f"    Found {len(assocs)} associations")
        except Exception as e:
            print(f"    No results or error: {e}")
        time.sleep(0.5)  # Be polite to the API

    out = RAW_DIR / "gwas_erap2_associations.json"
    out.write_text(json.dumps(all_assocs, indent=2), encoding="utf-8")
    print(f"  Saved {len(all_assocs)} ERAP2 associations → {out.name}")
    return all_assocs


def fetch_crohns_associations() -> list[dict]:
    """GWAS associations for Crohn's disease (EFO_0000384)."""
    print("Fetching Crohn's disease associations from GWAS Catalog ...")
    url = f"{BASE_URL}/efoTraits/EFO_0000384/associations"
    params = None
    assocs = _paginate_associations(url, params, max_pages=50)
    out = RAW_DIR / "gwas_crohns_associations.json"
    out.write_text(json.dumps(assocs, indent=2), encoding="utf-8")
    print(f"  Saved {len(assocs)} Crohn's associations → {out.name}")
    return assocs


def fetch_rs2549794() -> dict:
    """Fetch data for the specific Black Death → Crohn's variant."""
    print("Fetching rs2549794 variant details ...")
    url = f"{BASE_URL}/singleNucleotidePolymorphisms/rs2549794"
    data = _get(url)
    out = RAW_DIR / "gwas_rs2549794.json"
    out.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"  Saved rs2549794 details → {out.name}")
    return data


def extract_association_rows(assocs: list[dict], source_label: str) -> list[dict]:
    """Flatten GWAS association JSON into tabular rows."""
    rows = []
    for a in assocs:
        p_value = a.get("pvalue")
        risk_allele_info = a.get("riskFrequency", "")
        # Extract SNPs
        snps = a.get("snps", [])
        if not snps:
            # Try loci → strongestRiskAlleles
            loci = a.get("loci", [])
            for locus in loci:
                for sra in locus.get("strongestRiskAlleles", []):
                    risk_name = sra.get("riskAlleleName", "")
                    rsid = risk_name.split("-")[0] if "-" in risk_name else risk_name
                    rows.append({
                        "source": source_label,
                        "rsid": rsid,
                        "risk_allele": risk_name,
                        "p_value": p_value,
                        "risk_frequency": risk_allele_info,
                    })
        else:
            for snp in snps:
                rsid = snp.get("rsId", "")
                rows.append({
                    "source": source_label,
                    "rsid": rsid,
                    "risk_allele": "",
                    "p_value": p_value,
                    "risk_frequency": risk_allele_info,
                })
    return rows


def build_overlap(erap2_assocs: list[dict], crohns_assocs: list[dict]) -> pd.DataFrame:
    """Find variants that appear in BOTH ERAP2 and Crohn's association sets."""
    print("Building ERAP2 ∩ Crohn's overlap ...")
    erap2_rows = extract_association_rows(erap2_assocs, "ERAP2")
    crohns_rows = extract_association_rows(crohns_assocs, "Crohns")

    df_erap2 = pd.DataFrame(erap2_rows)
    df_crohns = pd.DataFrame(crohns_rows)

    print(f"  ERAP2 rows: {len(df_erap2)}, Crohn's rows: {len(df_crohns)}")

    if df_erap2.empty or df_crohns.empty:
        print("  WARNING: One dataset is empty — overlap will be empty.")
        overlap = pd.DataFrame(columns=["rsid", "erap2_p", "crohns_p"])
    else:
        erap2_snps = set(df_erap2["rsid"].dropna().unique())
        crohns_snps = set(df_crohns["rsid"].dropna().unique())
        common = erap2_snps & crohns_snps
        print(f"  Overlapping SNPs: {len(common)}")

        # Build overlap table
        overlap_rows = []
        for rsid in sorted(common):
            e_p = df_erap2.loc[df_erap2["rsid"] == rsid, "p_value"].min()
            c_p = df_crohns.loc[df_crohns["rsid"] == rsid, "p_value"].min()
            overlap_rows.append({"rsid": rsid, "erap2_p": e_p, "crohns_p": c_p})
        overlap = pd.DataFrame(overlap_rows)

    # Also check if rs2549794 is in either set
    target = "rs2549794"
    in_erap2 = target in set(df_erap2["rsid"]) if not df_erap2.empty else False
    in_crohns = target in set(df_crohns["rsid"]) if not df_crohns.empty else False
    print(f"  rs2549794 in ERAP2 set: {in_erap2}")
    print(f"  rs2549794 in Crohn's set: {in_crohns}")

    out = PROC_DIR / "erap2_crohns_overlap.csv"
    overlap.to_csv(out, index=False)
    print(f"  Saved overlap → {out.name} ({len(overlap)} SNPs)")

    # Save full tables too for reference
    if not df_erap2.empty:
        df_erap2.to_csv(PROC_DIR / "gwas_erap2_flat.csv", index=False)
    if not df_crohns.empty:
        df_crohns.to_csv(PROC_DIR / "gwas_crohns_flat.csv", index=False)

    return overlap


def main():
    print("=" * 60)
    print("GWAS Catalog Data Fetch — ERAP2 / Crohn's Disease")
    print("=" * 60)

    erap2 = fetch_erap2_associations()
    crohns = fetch_crohns_associations()
    _ = fetch_rs2549794()
    overlap = build_overlap(erap2, crohns)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  ERAP2 associations:  {len(erap2)}")
    print(f"  Crohn's associations: {len(crohns)}")
    print(f"  Overlapping SNPs:     {len(overlap)}")
    print(f"  Output directory:     {PROC_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
