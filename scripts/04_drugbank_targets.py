"""
04_drugbank_targets.py — Query local DrugBank SQLite for ERAP2-related data.

Uses the existing DrugBank database at pharmacy-ai-project/data/drugbank/drugbank.db
to find:
  1. Drugs targeting ERAP2 or related aminopeptidases
  2. SNP effects related to ERAP2 gene region
  3. Known drugs for Crohn's disease (therapeutic context)
  4. Aminopeptidase inhibitors (structural analogs for drug design)

Outputs:
  data/processed/drugbank_erap2_targets.csv
  data/processed/drugbank_crohns_drugs.csv
  data/processed/drugbank_aminopeptidase_inhibitors.csv
  data/processed/drugbank_snp_effects_erap2.csv
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import sqlite3
from pathlib import Path

import pandas as pd

DRUGBANK_DB = Path("C:/Users/jmj2z/Projects/medical/pharmacy-ai-project/data/drugbank/drugbank.db")
PROC_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)


def get_conn() -> sqlite3.Connection:
    if not DRUGBANK_DB.exists():
        raise FileNotFoundError(f"DrugBank DB not found: {DRUGBANK_DB}")
    conn = sqlite3.connect(str(DRUGBANK_DB))
    conn.row_factory = sqlite3.Row
    return conn


def query_erap2_targets(conn: sqlite3.Connection) -> pd.DataFrame:
    """Find drugs that target ERAP2 or related aminopeptidases."""
    print("Querying DrugBank for ERAP2/aminopeptidase targets ...")

    query = """
    SELECT DISTINCT
        d.drugbank_id,
        d.name AS drug_name,
        d.indication,
        d.mechanism_of_action,
        dt.target_name,
        dt.actions,
        dt.known_action,
        dt.organism
    FROM drugs d
    JOIN drug_targets dt ON d.drugbank_id = dt.drugbank_id
    WHERE (
        dt.target_name LIKE '%aminopeptidase%'
        OR dt.target_name LIKE '%ERAP%'
        OR dt.target_name LIKE '%endoplasmic reticulum aminopeptidase%'
        OR dt.target_name LIKE '%leucyl aminopeptidase%'
        OR dt.target_name LIKE '%puromycin%'
    )
    AND dt.organism = 'Humans'
    ORDER BY d.name
    """
    df = pd.read_sql_query(query, conn)
    print(f"  Found {len(df)} drug-target relationships for aminopeptidases")

    if not df.empty:
        unique_drugs = df["drug_name"].nunique()
        unique_targets = df["target_name"].nunique()
        print(f"  Unique drugs: {unique_drugs}")
        print(f"  Unique targets: {unique_targets}")
        print(f"  Target names: {df['target_name'].unique().tolist()[:10]}")

    out = PROC_DIR / "drugbank_erap2_targets.csv"
    df.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return df


def query_crohns_drugs(conn: sqlite3.Connection) -> pd.DataFrame:
    """Find approved drugs indicated for Crohn's disease."""
    print("\nQuerying DrugBank for Crohn's disease drugs ...")

    query = """
    SELECT DISTINCT
        d.drugbank_id,
        d.name AS drug_name,
        d.indication,
        d.mechanism_of_action,
        d.pharmacodynamics,
        d.half_life,
        d.route_of_elimination
    FROM drugs d
    JOIN drug_groups dg ON d.drugbank_id = dg.drugbank_id
    WHERE dg.group_name = 'approved'
    AND (
        d.indication LIKE '%Crohn%'
        OR d.indication LIKE '%inflammatory bowel%'
        OR d.indication LIKE '%IBD%'
    )
    ORDER BY d.name
    """
    df = pd.read_sql_query(query, conn)
    print(f"  Found {len(df)} approved Crohn's/IBD drugs")

    if not df.empty:
        print(f"  Drug names: {df['drug_name'].tolist()[:15]}")

    out = PROC_DIR / "drugbank_crohns_drugs.csv"
    df.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return df


def query_aminopeptidase_inhibitors(conn: sqlite3.Connection) -> pd.DataFrame:
    """Find known aminopeptidase inhibitors — structural references for drug design."""
    print("\nQuerying DrugBank for aminopeptidase inhibitors ...")

    query = """
    SELECT DISTINCT
        d.drugbank_id,
        d.name AS drug_name,
        d.description,
        d.mechanism_of_action,
        dt.target_name,
        dt.actions
    FROM drugs d
    JOIN drug_targets dt ON d.drugbank_id = dt.drugbank_id
    WHERE dt.target_name LIKE '%aminopeptidase%'
    AND dt.actions LIKE '%inhibitor%'
    AND dt.organism = 'Humans'
    ORDER BY d.name
    """
    df = pd.read_sql_query(query, conn)
    print(f"  Found {len(df)} aminopeptidase inhibitor relationships")

    out = PROC_DIR / "drugbank_aminopeptidase_inhibitors.csv"
    df.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return df


def query_snp_effects(conn: sqlite3.Connection) -> pd.DataFrame:
    """Find pharmacogenomic SNP effects related to ERAP2 or immune/IBD pathways."""
    print("\nQuerying DrugBank for relevant SNP effects ...")

    # Direct ERAP2 SNP effects
    query_erap2 = """
    SELECT
        d.drugbank_id,
        d.name AS drug_name,
        se.gene_symbol,
        se.rs_id,
        se.allele,
        se.defining_change,
        se.description,
        se.protein_name,
        se.uniprot_id
    FROM drugs d
    JOIN snp_effects se ON d.drugbank_id = se.drugbank_id
    WHERE se.gene_symbol IN ('ERAP2', 'ERAP1', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'NOD2', 'ATG16L1', 'IL23R')
    ORDER BY se.gene_symbol, d.name
    """
    df = pd.read_sql_query(query_erap2, conn)
    print(f"  Found {len(df)} SNP effects for ERAP2/immune pathway genes")

    if not df.empty:
        gene_counts = df["gene_symbol"].value_counts()
        print(f"  By gene: {gene_counts.to_dict()}")

    out = PROC_DIR / "drugbank_snp_effects_erap2.csv"
    df.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return df


def query_external_ids_erap2(conn: sqlite3.Connection) -> dict:
    """Find cross-references for aminopeptidase-targeting drugs (PubChem, ChEBI, etc.)."""
    print("\nQuerying external identifiers for ERAP2-targeting drugs ...")

    query = """
    SELECT DISTINCT
        d.drugbank_id,
        d.name,
        ei.resource,
        ei.identifier
    FROM drugs d
    JOIN drug_targets dt ON d.drugbank_id = dt.drugbank_id
    JOIN external_identifiers ei ON d.drugbank_id = ei.drugbank_id
    WHERE dt.target_name LIKE '%aminopeptidase%'
    AND dt.organism = 'Humans'
    AND ei.resource IN ('PubChem Compound', 'PubChem Substance', 'ChEBI', 'KEGG Drug', 'ChEMBL')
    ORDER BY d.name, ei.resource
    """
    df = pd.read_sql_query(query, conn)
    print(f"  Found {len(df)} external cross-references")

    if not df.empty:
        by_resource = df.groupby("resource").size()
        print(f"  By resource: {by_resource.to_dict()}")

    # Return as dict for easy lookup
    xrefs = {}
    for _, row in df.iterrows():
        drug_id = row["drugbank_id"]
        if drug_id not in xrefs:
            xrefs[drug_id] = {"name": row["name"], "xrefs": {}}
        xrefs[drug_id]["xrefs"][row["resource"]] = row["identifier"]

    return xrefs


def main():
    print("=" * 60)
    print("DrugBank Query — ERAP2 / Crohn's Disease Context")
    print("=" * 60)

    conn = get_conn()

    try:
        erap2_targets = query_erap2_targets(conn)
        crohns_drugs = query_crohns_drugs(conn)
        inhibitors = query_aminopeptidase_inhibitors(conn)
        snp_effects = query_snp_effects(conn)
        xrefs = query_external_ids_erap2(conn)

        # Summary
        print("\n" + "=" * 60)
        print("DRUGBANK SUMMARY")
        print("=" * 60)
        print(f"  Aminopeptidase-targeting drugs:  {erap2_targets['drug_name'].nunique() if not erap2_targets.empty else 0}")
        print(f"  Approved Crohn's/IBD drugs:      {len(crohns_drugs)}")
        print(f"  Aminopeptidase inhibitors:        {inhibitors['drug_name'].nunique() if not inhibitors.empty else 0}")
        print(f"  Relevant SNP effects:             {len(snp_effects)}")
        print(f"  External cross-references:        {len(xrefs)} drugs")
        print(f"\n  Output directory: {PROC_DIR}")
        print("Done.")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
