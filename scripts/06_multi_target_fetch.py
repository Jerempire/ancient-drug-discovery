"""
06_multi_target_fetch.py — Generalized data fetcher for all pathogen-selected immune targets.

Reads target definitions from targets.yaml and fetches:
  1. UniProt protein data + FASTA sequence
  2. OpenTargets disease associations + tractability
  3. GWAS associations for linked diseases
  4. PDB structure availability check

Works for any target in targets.yaml — not just ERAP2.

Outputs per target (in data/processed/<gene_symbol>/):
  <gene>_uniprot.json
  <gene>_<uniprot_id>.fasta
  <gene>_opentargets.json
  <gene>_gwas_associations.json
  <gene>_structure_info.json
  <gene>_target_summary.json
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
import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent
TARGETS_PATH = PROJECT_DIR / "targets.yaml"
DATA_DIR = PROJECT_DIR / "data"

GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"
OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v1/query"


def load_targets() -> dict:
    """Load all target definitions from targets.yaml."""
    data = yaml.safe_load(TARGETS_PATH.read_text(encoding="utf-8"))
    return data["targets"]


def _get(url: str, params: dict | None = None, headers: dict | None = None) -> dict:
    """GET with retry on 429."""
    if headers is None:
        headers = {"Accept": "application/json"}
    for attempt in range(3):
        resp = requests.get(url, headers=headers, params=params, timeout=60)
        if resp.status_code == 429:
            wait = int(resp.headers.get("Retry-After", 5))
            print(f"  Rate-limited, waiting {wait}s ...")
            time.sleep(wait)
            continue
        resp.raise_for_status()
        return resp.json()
    raise RuntimeError(f"Failed after 3 attempts: {url}")


def _gql(query: str) -> dict:
    """Execute OpenTargets GraphQL query."""
    resp = requests.post(OPENTARGETS_URL, json={"query": query}, timeout=60)
    resp.raise_for_status()
    return resp.json()


# ──────────────────────────────────────────────────────────────
# Data fetchers
# ──────────────────────────────────────────────────────────────

def fetch_uniprot(gene: str, uniprot_id: str, out_dir: Path) -> dict:
    """Fetch protein data from UniProt."""
    print(f"\n  [{gene}] Fetching UniProt {uniprot_id} ...")
    url = f"{UNIPROT_BASE}/{uniprot_id}.json"
    try:
        data = _get(url)
    except Exception as e:
        print(f"    WARNING: UniProt fetch failed: {e}")
        return {}

    # Save full JSON
    out = out_dir / f"{gene.lower()}_uniprot.json"
    out.write_text(json.dumps(data, indent=2), encoding="utf-8")

    # Extract key info
    protein_name = (
        data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "unknown")
    )
    seq = data.get("sequence", {}).get("value", "")
    seq_len = data.get("sequence", {}).get("length", 0)

    print(f"    Protein: {protein_name}")
    print(f"    Length: {seq_len} aa")

    # Save FASTA
    if seq:
        fasta_out = out_dir / f"{gene.lower()}_{uniprot_id}.fasta"
        fasta_out.write_text(
            f">sp|{uniprot_id}|{gene}_HUMAN {protein_name}\n{seq}\n",
            encoding="utf-8",
        )
        print(f"    Saved FASTA: {fasta_out.name}")

    # Extract functional annotations
    features = data.get("features", [])
    active_sites = [f for f in features if f.get("type") == "Active site"]
    binding_sites = [f for f in features if f.get("type") == "Binding site"]
    variants = [f for f in features if f.get("type") == "Natural variant"]

    print(f"    Active sites: {len(active_sites)}")
    print(f"    Binding sites: {len(binding_sites)}")
    print(f"    Natural variants: {len(variants)}")

    return {
        "uniprot_id": uniprot_id,
        "protein_name": protein_name,
        "sequence_length": seq_len,
        "active_sites": len(active_sites),
        "binding_sites": len(binding_sites),
        "natural_variants": len(variants),
    }


def fetch_opentargets(gene: str, ensembl_id: str, out_dir: Path) -> dict:
    """Fetch OpenTargets disease associations and tractability."""
    print(f"\n  [{gene}] Querying OpenTargets {ensembl_id} ...")

    query = f"""
    query {{
      target(ensemblId: "{ensembl_id}") {{
        id
        approvedSymbol
        approvedName
        biotype
        proteinAnnotations {{
          subcellularLocations
          functions
        }}
        tractability {{
          label
          modality
          value
        }}
        associatedDiseases(page: {{index: 0, size: 25}}) {{
          count
          rows {{
            disease {{
              id
              name
              therapeuticAreas {{
                id
                name
              }}
            }}
            score
            datatypeScores {{
              componentId: id
              score
            }}
          }}
        }}
      }}
    }}
    """

    try:
        result = _gql(query)
    except Exception as e:
        print(f"    WARNING: OpenTargets query failed: {e}")
        return {}

    target_data = result.get("data", {}).get("target", {})
    if not target_data:
        print(f"    WARNING: No OpenTargets data for {ensembl_id}")
        return {}

    out = out_dir / f"{gene.lower()}_opentargets.json"
    out.write_text(json.dumps(target_data, indent=2), encoding="utf-8")

    # Print summary
    assoc = target_data.get("associatedDiseases", {})
    count = assoc.get("count", 0)
    print(f"    Disease associations: {count}")

    # Find cancer-related associations
    cancer_assocs = []
    for row in assoc.get("rows", []):
        disease = row.get("disease", {})
        areas = disease.get("therapeuticAreas", [])
        area_ids = [a.get("id", "") for a in areas]
        if any("neoplasm" in aid.lower() or "cancer" in aid.lower() or "EFO_0000616" in aid for aid in area_ids):
            cancer_assocs.append({
                "disease": disease.get("name"),
                "efo_id": disease.get("id"),
                "score": row.get("score"),
            })

    if cancer_assocs:
        print(f"    Cancer-specific associations: {len(cancer_assocs)}")
        for ca in cancer_assocs[:5]:
            print(f"      {ca['disease']:40s} score={ca['score']:.3f}")

    # Tractability
    tract = target_data.get("tractability", [])
    if tract:
        print(f"    Tractability assessments: {len(tract)}")
        for t in tract[:5]:
            print(f"      {t.get('modality', '?')}: {t.get('label', '?')} = {t.get('value', '?')}")

    return {
        "total_associations": count,
        "cancer_associations": cancer_assocs,
        "tractability": tract,
    }


def fetch_gwas_for_diseases(gene: str, efo_ids: list[str], out_dir: Path) -> dict:
    """Fetch GWAS associations for diseases linked to this target."""
    if not efo_ids:
        print(f"\n  [{gene}] No GWAS EFO IDs configured — skipping GWAS fetch")
        return {}

    print(f"\n  [{gene}] Fetching GWAS associations for {len(efo_ids)} diseases ...")
    all_data = {}

    for efo_id in efo_ids:
        print(f"    Querying {efo_id} ...")
        url = f"{GWAS_BASE}/efoTraits/{efo_id}/associations"
        try:
            data = _get(url)
            assocs = data.get("_embedded", {}).get("associations", [])
            all_data[efo_id] = {
                "count": len(assocs),
                "sample_rsids": [],
            }
            # Extract sample rsIDs
            for a in assocs[:20]:
                loci = a.get("loci", [])
                for locus in loci:
                    for sra in locus.get("strongestRiskAlleles", []):
                        risk_name = sra.get("riskAlleleName", "")
                        rsid = risk_name.split("-")[0] if "-" in risk_name else risk_name
                        if rsid.startswith("rs"):
                            all_data[efo_id]["sample_rsids"].append(rsid)
            print(f"      Found {len(assocs)} associations")
        except Exception as e:
            print(f"      WARNING: GWAS fetch failed for {efo_id}: {e}")
            all_data[efo_id] = {"count": 0, "error": str(e)}
        time.sleep(0.5)

    out = out_dir / f"{gene.lower()}_gwas_associations.json"
    out.write_text(json.dumps(all_data, indent=2), encoding="utf-8")
    print(f"    Saved: {out.name}")
    return all_data


def check_pdb_structures(gene: str, target_cfg: dict, out_dir: Path) -> dict:
    """Check PDB for available experimental structures."""
    print(f"\n  [{gene}] Checking PDB structures ...")
    uniprot_id = target_cfg.get("uniprot_id", "")
    known_pdb = target_cfg.get("structure", {}).get("pdb_id")

    info = {
        "known_pdb": known_pdb,
        "alphafold_id": target_cfg.get("structure", {}).get("alphafold_id"),
        "alphafold_plddt": target_cfg.get("structure", {}).get("alphafold_plddt"),
        "pdb_search_results": [],
    }

    # Search PDB by UniProt ID
    if uniprot_id:
        search_query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "return_type": "entry",
            "request_options": {"results_content_type": ["experimental"]},
        }
        try:
            resp = requests.post(
                PDB_SEARCH_URL,
                json=search_query,
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            if resp.ok:
                results = resp.json().get("result_set", [])
                info["pdb_search_results"] = [r.get("identifier") for r in results[:10]]
                print(f"    PDB entries found: {len(results)} (showing first 10)")
                for pdb_id in info["pdb_search_results"][:5]:
                    print(f"      {pdb_id}")
            else:
                print(f"    PDB search returned {resp.status_code}")
        except Exception as e:
            print(f"    WARNING: PDB search failed: {e}")

    if known_pdb:
        print(f"    Known PDB: {known_pdb}")
    if info["alphafold_id"]:
        print(f"    AlphaFold: {info['alphafold_id']}")

    out = out_dir / f"{gene.lower()}_structure_info.json"
    out.write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info


def build_target_summary(
    gene: str,
    target_cfg: dict,
    uniprot_info: dict,
    opentargets_info: dict,
    gwas_info: dict,
    structure_info: dict,
    out_dir: Path,
) -> dict:
    """Build a comprehensive target summary for downstream pipeline stages."""
    summary = {
        "gene_symbol": gene,
        "ensembl_id": target_cfg.get("ensembl_id"),
        "uniprot_id": target_cfg.get("uniprot_id"),
        "protein": uniprot_info,
        "immune_step": target_cfg.get("immune_step"),
        "ancient_pathogen": target_cfg.get("ancient_pathogen"),
        "selection_event": target_cfg.get("selection_event"),
        "selection_era": target_cfg.get("selection_era"),
        "cancer_mechanism": target_cfg.get("cancer_mechanism", "").strip(),
        "cancer_types": target_cfg.get("cancer_types", []),
        "autoimmune_diseases": target_cfg.get("autoimmune_diseases", []),
        "drug_status": target_cfg.get("drug_status"),
        "known_compounds": target_cfg.get("known_compounds", []),
        "key_variants": target_cfg.get("key_variants", []),
        "opentargets": opentargets_info,
        "gwas": gwas_info,
        "structure": structure_info,
        "key_references": target_cfg.get("key_references", []),
        "pipeline_status": {
            "data_collection": "complete",
            "variant_analysis": "pending",
            "structure_prediction": "pending",
            "drug_design": "pending",
            "validation": "pending",
        },
    }

    out = out_dir / f"{gene.lower()}_target_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\n  [{gene}] Saved target summary: {out.name}")
    return summary


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def process_target(gene: str, target_cfg: dict) -> dict:
    """Run the full data collection pipeline for one target."""
    print(f"\n{'='*60}")
    print(f"Processing target: {gene}")
    print(f"  Pathogen: {target_cfg.get('ancient_pathogen')}")
    print(f"  Immune step: {target_cfg.get('immune_step')}")
    print(f"  Drug status: {target_cfg.get('drug_status', 'unknown')}")
    print(f"{'='*60}")

    # Create target-specific output directory
    out_dir = DATA_DIR / "processed" / gene.lower()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Fetch data
    uniprot_id = target_cfg.get("uniprot_id", "")
    ensembl_id = target_cfg.get("ensembl_id", "")
    efo_ids = target_cfg.get("gwas_efo_ids", [])

    uniprot_info = fetch_uniprot(gene, uniprot_id, out_dir) if uniprot_id else {}
    opentargets_info = fetch_opentargets(gene, ensembl_id, out_dir) if ensembl_id else {}
    gwas_info = fetch_gwas_for_diseases(gene, efo_ids, out_dir)
    structure_info = check_pdb_structures(gene, target_cfg, out_dir)

    # Build summary
    summary = build_target_summary(
        gene, target_cfg, uniprot_info, opentargets_info, gwas_info, structure_info, out_dir
    )
    return summary


def main():
    print("=" * 70)
    print("Multi-Target Data Fetch — Ancient Pathogen Selection & Cancer")
    print("=" * 70)

    targets = load_targets()
    print(f"\nLoaded {len(targets)} targets from targets.yaml")

    # Skip ERAP2 if it already has data (scripts 01-05 cover it)
    skip_existing = True

    summaries = {}
    for gene, cfg in targets.items():
        if skip_existing and gene == "ERAP2":
            erap2_summary = DATA_DIR / "processed" / "poc_target_summary.json"
            if erap2_summary.exists():
                print(f"\n  Skipping ERAP2 — existing data from scripts 01-05")
                continue

        try:
            summaries[gene] = process_target(gene, cfg)
        except Exception as e:
            print(f"\n  ERROR processing {gene}: {e}")
            summaries[gene] = {"error": str(e)}

    # Final summary
    print("\n" + "=" * 70)
    print("MULTI-TARGET FETCH SUMMARY")
    print("=" * 70)
    for gene, summary in summaries.items():
        if "error" in summary:
            print(f"  {gene:12s} ERROR: {summary['error']}")
        else:
            n_cancer = len(summary.get("opentargets", {}).get("cancer_associations", []))
            n_struct = len(summary.get("structure", {}).get("pdb_search_results", []))
            print(f"  {gene:12s} cancer_assocs={n_cancer:2d}  pdb_structures={n_struct:2d}  drug_status={summary.get('drug_status', '?')[:50]}")

    print(f"\n  Output: {DATA_DIR / 'processed'}")
    print("Done.")


if __name__ == "__main__":
    main()
