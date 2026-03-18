"""
03_variant_analysis.py — Cross-reference ancient vs. modern variants at the ERAP2 locus.

This script:
  1. Loads GWAS ERAP2/Crohn's overlap data (from script 01)
  2. Loads AADR ERAP2-region SNPs (from script 02)
  3. Loads ClinVar ERAP2 data (from script 02)
  4. Queries OpenTargets for ERAP2 target-disease evidence
  5. Builds a unified variant annotation table
  6. Identifies the most promising variants for the pipeline PoC

Outputs:
  data/processed/erap2_variant_master.csv  — unified variant table
  data/processed/erap2_opentargets.json    — OpenTargets evidence
  data/processed/poc_target_summary.json   — summary for downstream pipeline stages
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
from pathlib import Path

import pandas as pd
import requests

PROC_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"

# ERAP2 identifiers
ERAP2_ENSEMBL = "ENSG00000164308"
ERAP2_UNIPROT = "Q6P179"

# Key variants from the Black Death / Crohn's literature
LITERATURE_VARIANTS = {
    "rs2549794": {
        "context": "Black Death survival variant → Crohn's risk",
        "source": "Klunk et al. 2022 Nature",
        "chr": "5",
        "pos_grch37": 96_243_541,
        "functional_impact": "Affects ERAP2 expression (eQTL), alters antigen presentation",
        "ancient_selection": "Positive selection during Black Death (1346-1353)",
    },
    "rs2248374": {
        "context": "ERAP2 splice variant — determines transcript A vs. B",
        "source": "Andrés et al. 2010 Nature Genetics",
        "chr": "5",
        "pos_grch37": 96_244_549,
        "functional_impact": "A allele → full-length ERAP2; G allele → truncated/NMD",
        "ancient_selection": "Balancing selection maintained both alleles for >700,000 years",
    },
    "rs2910686": {
        "context": "ERAP2 regulatory variant linked to Crohn's",
        "source": "GWAS Catalog / IBD genetics consortium",
        "chr": "5",
        "pos_grch37": 96_252_589,
        "functional_impact": "Regulatory — affects ERAP2 expression level",
        "ancient_selection": "Part of ERAP2 selective sweep haplotype",
    },
}


def load_gwas_overlap() -> pd.DataFrame:
    """Load the ERAP2 ∩ Crohn's overlap from script 01."""
    path = PROC_DIR / "erap2_crohns_overlap.csv"
    if path.exists():
        df = pd.read_csv(path)
        print(f"  Loaded GWAS overlap: {len(df)} SNPs")
        return df
    print("  WARNING: GWAS overlap file not found — run script 01 first")
    return pd.DataFrame()


def load_aadr_erap2_snps() -> pd.DataFrame:
    """Load AADR ERAP2-region SNPs from script 02."""
    path = PROC_DIR / "aadr_erap2_snps.csv"
    if path.exists():
        df = pd.read_csv(path)
        print(f"  Loaded AADR ERAP2 SNPs: {len(df)}")
        return df
    print("  WARNING: AADR ERAP2 SNPs not found — run script 02 first")
    return pd.DataFrame()


def load_clinvar() -> list[dict]:
    """Load ClinVar ERAP2 data from script 02."""
    path = RAW_DIR / "clinvar_erap2.json"
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        print(f"  Loaded ClinVar ERAP2: {len(data)} entries")
        return data
    print("  WARNING: ClinVar data not found — run script 02 first")
    return []


def query_opentargets_erap2() -> dict:
    """Query OpenTargets GraphQL API for ERAP2 target evidence."""
    print("\nQuerying OpenTargets for ERAP2 (ENSG00000164308) ...")
    url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Query 1: Target info
    query_target = """
    query {
      target(ensemblId: "ENSG00000164308") {
        id
        approvedSymbol
        approvedName
        biotype
        proteinAnnotations {
          subcellularLocations
          functions
        }
        tractability {
          label
          modality
          value
        }
      }
    }
    """

    # Query 2: Disease associations
    query_diseases = """
    query {
      target(ensemblId: "ENSG00000164308") {
        associatedDiseases(page: {index: 0, size: 25}) {
          count
          rows {
            disease {
              id
              name
            }
            score
            datatypeScores {
              componentId: id
              score
            }
          }
        }
      }
    }
    """

    results = {}

    # Fetch target info
    resp = requests.post(url, json={"query": query_target}, timeout=60)
    if resp.ok:
        target_data = resp.json().get("data", {}).get("target", {})
        results["target_info"] = target_data
        print(f"  Target: {target_data.get('approvedSymbol')} — {target_data.get('approvedName')}")

        # Tractability
        tract = target_data.get("tractability", [])
        if tract:
            print("  Tractability:")
            for t in tract:
                print(f"    {t.get('modality', '?')}: {t.get('label', '?')} = {t.get('value', '?')}")
    else:
        print(f"  WARNING: Target query failed ({resp.status_code})")

    # Fetch disease associations
    resp = requests.post(url, json={"query": query_diseases}, timeout=60)
    if resp.ok:
        assoc_data = resp.json().get("data", {}).get("target", {}).get("associatedDiseases", {})
        results["disease_associations"] = assoc_data
        count = assoc_data.get("count", 0)
        print(f"\n  Disease associations: {count}")
        for row in assoc_data.get("rows", [])[:10]:
            disease = row.get("disease", {})
            print(f"    {disease.get('name', '?'):40s} score={row.get('score', 0):.3f}  ({disease.get('id', '')})")
    else:
        print(f"  WARNING: Disease query failed ({resp.status_code})")

    out = PROC_DIR / "erap2_opentargets.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"  Saved → {out.name}")
    return results


def build_master_table(
    gwas_overlap: pd.DataFrame,
    aadr_snps: pd.DataFrame,
    clinvar: list[dict],
) -> pd.DataFrame:
    """Build unified variant annotation table."""
    print("\nBuilding master variant table ...")

    rows = []

    # Start with literature variants (known important)
    for rsid, info in LITERATURE_VARIANTS.items():
        row = {
            "rsid": rsid,
            "chr": info["chr"],
            "pos_grch37": info["pos_grch37"],
            "source": "literature",
            "context": info["context"],
            "literature_source": info["source"],
            "functional_impact": info["functional_impact"],
            "ancient_selection": info["ancient_selection"],
            "in_gwas_overlap": False,
            "in_aadr_panel": False,
            "clinvar_significance": "",
        }

        # Check if in GWAS overlap
        if not gwas_overlap.empty and rsid in gwas_overlap["rsid"].values:
            row["in_gwas_overlap"] = True
            match = gwas_overlap[gwas_overlap["rsid"] == rsid].iloc[0]
            row["gwas_erap2_p"] = match.get("erap2_p")
            row["gwas_crohns_p"] = match.get("crohns_p")

        # Check if in AADR panel
        if not aadr_snps.empty and rsid in aadr_snps["snp_id"].values:
            row["in_aadr_panel"] = True

        rows.append(row)

    # Add GWAS overlap variants not in literature
    if not gwas_overlap.empty:
        for _, grow in gwas_overlap.iterrows():
            rsid = grow["rsid"]
            if rsid in LITERATURE_VARIANTS:
                continue
            rows.append({
                "rsid": rsid,
                "source": "gwas_overlap",
                "in_gwas_overlap": True,
                "gwas_erap2_p": grow.get("erap2_p"),
                "gwas_crohns_p": grow.get("crohns_p"),
                "in_aadr_panel": rsid in aadr_snps["snp_id"].values if not aadr_snps.empty else False,
            })

    # Add AADR panel SNPs in ERAP2 region not already listed
    if not aadr_snps.empty:
        existing_rsids = {r["rsid"] for r in rows}
        for _, arow in aadr_snps.iterrows():
            rsid = arow["snp_id"]
            if rsid in existing_rsids or not rsid.startswith("rs"):
                continue
            rows.append({
                "rsid": rsid,
                "chr": str(arow["chr"]),
                "pos_grch37": arow["physical_pos"],
                "source": "aadr_panel",
                "in_gwas_overlap": False,
                "in_aadr_panel": True,
            })

    df = pd.DataFrame(rows)
    print(f"  Master table: {len(df)} variants")
    print(f"    From literature: {len([r for r in rows if r.get('source') == 'literature'])}")
    print(f"    From GWAS overlap: {len([r for r in rows if r.get('source') == 'gwas_overlap'])}")
    print(f"    From AADR panel: {len([r for r in rows if r.get('source') == 'aadr_panel'])}")

    out = PROC_DIR / "erap2_variant_master.csv"
    df.to_csv(out, index=False)
    print(f"  Saved → {out.name}")
    return df


def build_poc_summary(master_df: pd.DataFrame, opentargets: dict) -> dict:
    """Build a summary JSON for downstream pipeline stages."""
    print("\nBuilding PoC target summary ...")

    # Primary target variant
    primary = LITERATURE_VARIANTS["rs2549794"]

    # ERAP2 protein info
    uniprot_path = RAW_DIR / "uniprot_erap2_Q6P179.json"
    protein_info = {}
    if uniprot_path.exists():
        up_data = json.loads(uniprot_path.read_text(encoding="utf-8"))
        protein_info = {
            "uniprot_id": ERAP2_UNIPROT,
            "ensembl_id": ERAP2_ENSEMBL,
            "sequence_length": up_data.get("sequence", {}).get("length"),
            "function": "Aminopeptidase that trims peptides for MHC class I antigen presentation",
        }

    # Disease links from OpenTargets
    disease_links = []
    if opentargets.get("disease_associations"):
        for row in opentargets["disease_associations"].get("rows", [])[:5]:
            disease_links.append({
                "disease": row["disease"]["name"],
                "efo_id": row["disease"]["id"],
                "score": row["score"],
            })

    # Tractability
    tractability = {}
    if opentargets.get("target_info", {}).get("tractability"):
        for t in opentargets["target_info"]["tractability"]:
            tractability[t.get("label", "")] = {
                "modality": t.get("modality"),
                "value": t.get("value"),
            }

    summary = {
        "project": "Ancient Drug Discovery PoC",
        "target_gene": "ERAP2",
        "target_protein": protein_info,
        "primary_variant": {
            "rsid": "rs2549794",
            **primary,
        },
        "secondary_variants": [
            {"rsid": k, **v}
            for k, v in LITERATURE_VARIANTS.items()
            if k != "rs2549794"
        ],
        "disease_links": disease_links,
        "tractability": tractability,
        "total_variants_cataloged": len(master_df),
        "pipeline_next_steps": {
            "stage_2": "Run ESM Cambrian on ERAP2 wild-type vs. rs2549794 mutant to score functional impact",
            "stage_3": "Predict ERAP2 wild-type and mutant 3D structures via AlphaFold3/OpenFold3",
            "stage_4": "Dock known inhibitors (e.g., DG013A) + generate novel candidates via LigandForge/DrugGPT",
            "stage_5": "Score candidates with Boltz-2 + OpenADMET filtering",
        },
        "known_inhibitors": [
            {
                "name": "DG013A",
                "type": "Small molecule",
                "mechanism": "Active site inhibitor of ERAP2 aminopeptidase",
                "reference": "Zervoudi et al. 2013 PNAS",
                "use": "Validation — compare our predictions against known binder",
            },
        ],
        "fasta_file": "data/processed/erap2_Q6P179.fasta",
    }

    out = PROC_DIR / "poc_target_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"  Saved → {out.name}")
    return summary


def main():
    print("=" * 60)
    print("Variant Analysis — ERAP2 Ancient-to-Modern Bridge")
    print("=" * 60)

    # Load data from prior scripts
    print("\nLoading data from prior stages ...")
    gwas_overlap = load_gwas_overlap()
    aadr_snps = load_aadr_erap2_snps()
    clinvar = load_clinvar()

    # Query OpenTargets
    opentargets = query_opentargets_erap2()

    # Build unified table
    master = build_master_table(gwas_overlap, aadr_snps, clinvar)

    # Build PoC summary
    summary = build_poc_summary(master, opentargets)

    # Final report
    print("\n" + "=" * 60)
    print("VARIANT ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"  Total variants cataloged:  {len(master)}")
    print(f"  Literature key variants:   {len(LITERATURE_VARIANTS)}")
    print(f"  In GWAS overlap:           {master['in_gwas_overlap'].sum() if 'in_gwas_overlap' in master else 0}")
    print(f"  In AADR panel:             {master['in_aadr_panel'].sum() if 'in_aadr_panel' in master else 0}")
    print(f"  OpenTargets diseases:      {len(summary.get('disease_links', []))}")
    print(f"\n  Primary target: ERAP2 / rs2549794")
    print(f"  Known inhibitor: DG013A (validation reference)")
    print(f"\n  Next: Run notebooks/04_esm_variant_effects.ipynb on Colab")
    print("Done.")


if __name__ == "__main__":
    main()
