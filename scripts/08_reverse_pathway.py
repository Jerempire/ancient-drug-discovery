"""
08_reverse_pathway.py — Reverse pathway drug repurposing engine.

Systematic drug repurposing via reverse graph walks:
  Drug → Protein Target → Gene → GWAS Disease Association → New Indication

Data sources:
  - DrugBank SQLite: drug → protein target mappings
  - OpenTargets GraphQL: gene → disease evidence scores
  - MedGraph-Rx: drug → disease → adverse event graph
  - targets.yaml: known targets and cancer types

This is high-leverage because repurposing skips most of the pipeline
(existing safety data, manufacturing, formulation).

Outputs:
  data/processed/repurposing_opportunities.json
  data/processed/repurposing_opportunities.csv
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
PROC_DIR = PROJECT_DIR / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# DrugBank SQLite (referenced from pharmacy-ai-project)
DRUGBANK_DB = PROJECT_DIR.parent / "pharmacy-ai-project" / "data" / "drugbank" / "drugbank.db"


def load_targets() -> dict:
    data = yaml.safe_load(TARGETS_PATH.read_text(encoding="utf-8"))
    return data["targets"]


def _gql(query: str) -> dict:
    resp = requests.post(OPENTARGETS_URL, json={"query": query}, timeout=60)
    resp.raise_for_status()
    return resp.json()


def get_known_drugs_for_target(gene: str, ensembl_id: str) -> list[dict]:
    """Fetch known drugs targeting this gene from OpenTargets."""
    query = f"""
    query {{
      target(ensemblId: "{ensembl_id}") {{
        knownDrugs(size: 50) {{
          count
          rows {{
            drug {{
              id
              name
              mechanismOfAction
              drugType
              isApproved
            }}
            disease {{
              id
              name
            }}
            phase
            status
            urls {{
              name
              url
            }}
          }}
        }}
      }}
    }}
    """
    try:
        result = _gql(query)
    except Exception as e:
        print(f"    WARNING: OpenTargets knownDrugs query failed for {gene}: {e}")
        return []

    rows = result.get("data", {}).get("target", {}).get("knownDrugs", {}).get("rows", [])
    drugs = []
    for row in rows:
        drug = row.get("drug", {})
        disease = row.get("disease", {})
        drugs.append({
            "drug_id": drug.get("id", ""),
            "drug_name": drug.get("name", ""),
            "mechanism": drug.get("mechanismOfAction", ""),
            "drug_type": drug.get("drugType", ""),
            "approved": drug.get("isApproved", False),
            "indication_disease": disease.get("name", ""),
            "indication_efo": disease.get("id", ""),
            "phase": row.get("phase"),
            "status": row.get("status", ""),
        })

    return drugs


def get_disease_associations_for_gene(ensembl_id: str) -> list[dict]:
    """Get ALL disease associations for a gene (not just cancer)."""
    query = f"""
    query {{
      target(ensemblId: "{ensembl_id}") {{
        associatedDiseases(page: {{index: 0, size: 200}}) {{
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
          }}
        }}
      }}
    }}
    """
    try:
        result = _gql(query)
    except Exception:
        return []

    rows = result.get("data", {}).get("target", {}).get("associatedDiseases", {}).get("rows", [])
    associations = []
    for row in rows:
        disease = row.get("disease", {})
        areas = disease.get("therapeuticAreas", [])
        associations.append({
            "disease_name": disease.get("name", ""),
            "disease_id": disease.get("id", ""),
            "score": row.get("score", 0),
            "therapeutic_areas": [a.get("name", "") for a in areas],
        })

    return associations


def score_repurposing_opportunity(
    drug: dict,
    new_indication: dict,
    existing_indications: list[str],
) -> dict:
    """Score a repurposing opportunity.

    High score = strong evidence for new indication + no existing approval there.
    """
    scores = {
        "evidence_strength": 0.0,   # OpenTargets association score
        "novelty": 0.0,             # Is this a genuinely new indication?
        "approval_advantage": 0.0,  # Does existing approval help?
        "commercial_potential": 0.0, # Is the new indication commercially interesting?
    }

    # Evidence strength from OpenTargets
    ot_score = new_indication.get("score", 0)
    scores["evidence_strength"] = min(1.0, ot_score * 1.2)

    # Novelty: higher if the drug isn't already approved/trialed for this disease
    ind_name = new_indication.get("disease_name", "").lower()
    if any(ind_name in ei.lower() for ei in existing_indications):
        scores["novelty"] = 0.0  # Already known
    else:
        scores["novelty"] = 0.8

    # Approval advantage: approved drugs have safety data
    if drug.get("approved"):
        scores["approval_advantage"] = 0.9
    elif drug.get("phase") and drug["phase"] >= 2:
        scores["approval_advantage"] = 0.6
    else:
        scores["approval_advantage"] = 0.3

    # Commercial potential: oncology and rare disease command premium
    areas = [a.lower() for a in new_indication.get("therapeutic_areas", [])]
    if any("neoplasm" in a or "cancer" in a for a in areas):
        scores["commercial_potential"] = 0.85
    elif any("rare" in a for a in areas):
        scores["commercial_potential"] = 0.80
    elif any("autoimmune" in a or "immune" in a for a in areas):
        scores["commercial_potential"] = 0.65
    else:
        scores["commercial_potential"] = 0.50

    # Weighted composite
    weights = {
        "evidence_strength": 0.35,
        "novelty": 0.25,
        "approval_advantage": 0.20,
        "commercial_potential": 0.20,
    }
    composite = sum(scores[k] * w for k, w in weights.items())

    return {
        "scores": scores,
        "composite": round(composite, 4),
        "verdict": "INVESTIGATE" if composite >= 0.5 and scores["novelty"] > 0 else "SKIP",
    }


def find_repurposing_opportunities(targets: dict) -> list[dict]:
    """Run reverse pathway analysis across all targets."""
    opportunities = []

    for gene, cfg in targets.items():
        ensembl_id = cfg.get("ensembl_id", "")
        if not ensembl_id:
            continue

        print(f"\n{'='*60}")
        print(f"  {gene} — Reverse Pathway Analysis")
        print(f"{'='*60}")

        # Step 1: Get known drugs for this target
        known_drugs = get_known_drugs_for_target(gene, ensembl_id)
        time.sleep(0.5)
        print(f"  Found {len(known_drugs)} known drug-target relationships")

        # Step 2: Get ALL disease associations for this gene
        all_diseases = get_disease_associations_for_gene(ensembl_id)
        time.sleep(0.5)
        print(f"  Found {len(all_diseases)} disease associations")

        # Step 3: For each drug, find new indications
        # (diseases associated with the gene but NOT the drug's current indication)
        drug_map: dict[str, dict] = {}
        for d in known_drugs:
            name = d["drug_name"]
            if name not in drug_map:
                drug_map[name] = {
                    "drug": d,
                    "existing_indications": [],
                }
            drug_map[name]["existing_indications"].append(d["indication_disease"])

        for drug_name, drug_info in drug_map.items():
            drug = drug_info["drug"]
            existing = drug_info["existing_indications"]

            for disease in all_diseases:
                disease_name = disease["disease_name"]
                # Skip if already an indication
                if any(disease_name.lower() in ei.lower() for ei in existing):
                    continue

                # Score the opportunity
                scoring = score_repurposing_opportunity(drug, disease, existing)

                if scoring["composite"] < 0.3:
                    continue  # Skip low-scoring

                opp = {
                    "target_gene": gene,
                    "drug_name": drug_name,
                    "drug_approved": drug.get("approved", False),
                    "mechanism": drug.get("mechanism", ""),
                    "current_indications": existing[:3],
                    "new_indication": disease_name,
                    "new_indication_id": disease.get("disease_id", ""),
                    "therapeutic_areas": disease.get("therapeutic_areas", []),
                    "ot_evidence_score": disease.get("score", 0),
                    "repurposing_score": scoring["composite"],
                    "score_breakdown": scoring["scores"],
                    "verdict": scoring["verdict"],
                }
                opportunities.append(opp)

        # Also check: compounds from targets.yaml that aren't in OpenTargets
        for cpd in cfg.get("known_compounds", []):
            cpd_name = cpd.get("name", "")
            if cpd_name in drug_map:
                continue  # Already processed
            for disease in all_diseases:
                cancer_types = cfg.get("cancer_types", [])
                disease_name = disease["disease_name"]
                if any(ct in disease_name.lower() for ct in cancer_types):
                    continue  # Skip known cancer types

                scoring = score_repurposing_opportunity(
                    {"approved": False, "phase": 0},
                    disease,
                    [ct.replace("_", " ") for ct in cancer_types],
                )
                if scoring["composite"] < 0.3:
                    continue

                opp = {
                    "target_gene": gene,
                    "drug_name": cpd_name,
                    "drug_approved": False,
                    "mechanism": cpd.get("mechanism", ""),
                    "current_indications": cancer_types[:3],
                    "new_indication": disease_name,
                    "new_indication_id": disease.get("disease_id", ""),
                    "therapeutic_areas": disease.get("therapeutic_areas", []),
                    "ot_evidence_score": disease.get("score", 0),
                    "repurposing_score": scoring["composite"],
                    "score_breakdown": scoring["scores"],
                    "verdict": scoring["verdict"],
                }
                opportunities.append(opp)

    # Sort by score
    opportunities.sort(key=lambda x: x["repurposing_score"], reverse=True)
    return opportunities


def main():
    print("=" * 60)
    print("Reverse Pathway Drug Repurposing Engine")
    print("=" * 60)

    targets = load_targets()
    print(f"\nLoaded {len(targets)} targets from targets.yaml")

    opportunities = find_repurposing_opportunities(targets)

    # Save JSON
    json_out = PROC_DIR / "repurposing_opportunities.json"
    json_out.write_text(json.dumps(opportunities, indent=2), encoding="utf-8")

    # Save CSV
    if opportunities:
        df = pd.DataFrame([{
            "target": o["target_gene"],
            "drug": o["drug_name"],
            "approved": o["drug_approved"],
            "mechanism": o["mechanism"][:60],
            "new_indication": o["new_indication"],
            "ot_score": o["ot_evidence_score"],
            "repurposing_score": o["repurposing_score"],
            "verdict": o["verdict"],
        } for o in opportunities])
        csv_out = PROC_DIR / "repurposing_opportunities.csv"
        df.to_csv(csv_out, index=False)
        print(f"\nSaved {len(opportunities)} opportunities to {csv_out.name}")

    # Print top opportunities
    print(f"\n{'='*70}")
    print("TOP REPURPOSING OPPORTUNITIES")
    print(f"{'='*70}")
    investigate = [o for o in opportunities if o["verdict"] == "INVESTIGATE"]
    print(f"\n{len(investigate)} actionable opportunities out of {len(opportunities)} total\n")

    for i, opp in enumerate(investigate[:10], 1):
        approved_tag = " [APPROVED]" if opp["drug_approved"] else ""
        print(f"  {i}. {opp['drug_name']}{approved_tag} -> {opp['new_indication']}")
        print(f"     Target: {opp['target_gene']} | Score: {opp['repurposing_score']:.2f}")
        print(f"     Current: {', '.join(opp['current_indications'][:2])}")
        print(f"     Evidence: OT={opp['ot_evidence_score']:.2f} | Areas: {', '.join(opp['therapeutic_areas'][:2])}")
        print()

    print("Done.")


if __name__ == "__main__":
    main()
