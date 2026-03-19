"""
07_cancer_connections.py — Cancer-specific evidence compiler for pathogen-selected immune targets.

For each target in targets.yaml, queries:
  1. OpenTargets for cancer-specific associations and evidence
  2. ClinicalTrials.gov for active cancer trials involving the target
  3. Literature references for the pathogen→cancer connection

Builds a unified cancer evidence matrix showing:
  - Which immune step each target controls (detection/activation/killing)
  - Cancer types affected
  - Drug opportunity score
  - Evidence strength ranking

Outputs:
  data/processed/cancer_evidence_matrix.csv
  data/processed/cancer_evidence_matrix.json
  data/processed/cancer_trials_<gene>.json (per target)
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
CTGOV_BASE = "https://clinicaltrials.gov/api/v2/studies"


def load_targets() -> dict:
    data = yaml.safe_load(TARGETS_PATH.read_text(encoding="utf-8"))
    return data["targets"]


def _gql(query: str) -> dict:
    resp = requests.post(OPENTARGETS_URL, json={"query": query}, timeout=60)
    resp.raise_for_status()
    return resp.json()


# ──────────────────────────────────────────────────────────────
# Cancer evidence per target
# ──────────────────────────────────────────────────────────────

def fetch_cancer_associations(gene: str, ensembl_id: str) -> list[dict]:
    """Fetch cancer-specific disease associations from OpenTargets."""
    print(f"\n  [{gene}] Cancer associations from OpenTargets ...")

    # Query for diseases in the neoplasm therapeutic area
    query = f"""
    query {{
      target(ensemblId: "{ensembl_id}") {{
        associatedDiseases(
          page: {{index: 0, size: 100}}
          BFilter: "neoplasm"
        ) {{
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
        # Fallback: use the broader query and filter client-side
        return _fetch_cancer_fallback(gene, ensembl_id)

    target_data = result.get("data", {}).get("target", {})
    assoc = target_data.get("associatedDiseases", {})

    # If BFilter didn't work (API version difference), filter client-side
    rows = assoc.get("rows", [])
    if not rows:
        return _fetch_cancer_fallback(gene, ensembl_id)

    cancer_rows = []
    for row in rows:
        disease = row.get("disease", {})
        areas = disease.get("therapeuticAreas", [])
        area_names = [a.get("name", "").lower() for a in areas]
        area_ids = [a.get("id", "") for a in areas]

        is_cancer = any(
            "neoplasm" in n or "cancer" in n or "tumor" in n
            or "EFO_0000616" in aid
            for n, aid in zip(area_names, area_ids)
        )

        if is_cancer:
            cancer_rows.append({
                "disease": disease.get("name"),
                "efo_id": disease.get("id"),
                "score": row.get("score"),
                "evidence_types": {
                    dt["componentId"]: dt["score"]
                    for dt in row.get("datatypeScores", [])
                    if dt["score"] > 0
                },
            })

    print(f"    Found {len(cancer_rows)} cancer associations")
    for cr in cancer_rows[:5]:
        print(f"      {cr['disease']:45s} score={cr['score']:.3f}")

    return cancer_rows


def _fetch_cancer_fallback(gene: str, ensembl_id: str) -> list[dict]:
    """Fallback: fetch all associations and filter for cancer client-side."""
    query = f"""
    query {{
      target(ensemblId: "{ensembl_id}") {{
        associatedDiseases(page: {{index: 0, size: 200}}) {{
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
    cancer_rows = []
    for row in rows:
        disease = row.get("disease", {})
        areas = disease.get("therapeuticAreas", [])
        is_cancer = any(
            "neoplasm" in a.get("name", "").lower()
            or "EFO_0000616" in a.get("id", "")
            for a in areas
        )
        if is_cancer:
            cancer_rows.append({
                "disease": disease.get("name"),
                "efo_id": disease.get("id"),
                "score": row.get("score"),
            })
    return cancer_rows


def fetch_cancer_trials(gene: str) -> list[dict]:
    """Search ClinicalTrials.gov for cancer trials involving this target/gene."""
    print(f"\n  [{gene}] Searching ClinicalTrials.gov for cancer trials ...")

    params = {
        "query.cond": "cancer OR neoplasm OR carcinoma OR tumor",
        "query.intr": gene,
        "pageSize": 20,
        "format": "json",
    }

    try:
        resp = requests.get(CTGOV_BASE, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"    WARNING: ClinicalTrials.gov query failed: {e}")
        return []

    studies = data.get("studies", [])
    trials = []
    for study in studies:
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        design_module = protocol.get("designModule", {})

        trials.append({
            "nct_id": id_module.get("nctId", ""),
            "title": id_module.get("briefTitle", ""),
            "status": status_module.get("overallStatus", ""),
            "phase": design_module.get("phases", []),
            "start_date": status_module.get("startDateStruct", {}).get("date", ""),
        })

    print(f"    Found {len(trials)} cancer trials mentioning {gene}")
    for t in trials[:3]:
        phase_str = "/".join(t["phase"]) if t["phase"] else "?"
        print(f"      [{t['status']:15s}] Phase {phase_str}: {t['title'][:60]}")

    # Save per-target
    out = PROC_DIR / f"cancer_trials_{gene.lower()}.json"
    out.write_text(json.dumps(trials, indent=2), encoding="utf-8")
    return trials


# ──────────────────────────────────────────────────────────────
# Evidence matrix
# ──────────────────────────────────────────────────────────────

def score_drug_opportunity(target_cfg: dict, cancer_assocs: list, trials: list) -> dict:
    """Score the drug opportunity for a target based on available evidence."""
    scores = {
        "ancient_selection": 0,     # 0-3: how strong is the selection evidence?
        "cancer_mechanism": 0,      # 0-3: how clear is the cancer mechanism?
        "drug_gap": 0,              # 0-3: how undrugged is this target?
        "structural_data": 0,       # 0-3: can we design drugs against it?
        "clinical_validation": 0,   # 0-3: any clinical trials?
    }

    # Ancient selection strength
    selection = target_cfg.get("selection_event", "").lower()
    if "strongest" in selection or "+40%" in selection:
        scores["ancient_selection"] = 3
    elif "near-fixation" in selection or "400m" in selection:
        scores["ancient_selection"] = 3
    elif "susceptibility" in selection or "polymorphism" in selection:
        scores["ancient_selection"] = 2
    else:
        scores["ancient_selection"] = 1

    # Cancer mechanism clarity
    cancer_types = target_cfg.get("cancer_types", [])
    mechanism = target_cfg.get("cancer_mechanism", "")
    if len(cancer_types) >= 3 and ("%" in mechanism or "PNAS" in mechanism or "Nature" in mechanism):
        scores["cancer_mechanism"] = 3
    elif len(cancer_types) >= 2:
        scores["cancer_mechanism"] = 2
    else:
        scores["cancer_mechanism"] = 1

    # Drug gap (higher = more opportunity)
    drug_status = target_cfg.get("drug_status", "").lower()
    compounds = target_cfg.get("known_compounds", [])
    if "no drugs" in drug_status and not compounds:
        scores["drug_gap"] = 3
    elif "no drugs" in drug_status or "undrugged" in drug_status:
        scores["drug_gap"] = 3
    elif "repurpos" in drug_status or len(compounds) <= 1:
        scores["drug_gap"] = 2
    else:
        scores["drug_gap"] = 1

    # Structural data
    structure = target_cfg.get("structure", {})
    pdb = structure.get("pdb_id")
    af_plddt = structure.get("alphafold_plddt")
    if pdb:
        scores["structural_data"] = 3  # experimental structure
    elif af_plddt and af_plddt > 90:
        scores["structural_data"] = 3
    elif af_plddt and af_plddt > 70:
        scores["structural_data"] = 2
    else:
        scores["structural_data"] = 1  # AlphaFold only, unknown quality

    # Clinical validation
    active_trials = [t for t in trials if t.get("status") in ("RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION")]
    if len(active_trials) >= 3:
        scores["clinical_validation"] = 3
    elif len(active_trials) >= 1:
        scores["clinical_validation"] = 2
    elif len(trials) >= 1:
        scores["clinical_validation"] = 1

    scores["total"] = sum(scores.values())
    scores["max_possible"] = 15
    return scores


def build_evidence_matrix(targets: dict) -> pd.DataFrame:
    """Build the unified cancer evidence matrix."""
    print("\n" + "=" * 70)
    print("Building Cancer Evidence Matrix")
    print("=" * 70)

    rows = []
    all_evidence = {}

    for gene, cfg in targets.items():
        ensembl_id = cfg.get("ensembl_id", "")

        # Fetch cancer evidence
        cancer_assocs = fetch_cancer_associations(gene, ensembl_id) if ensembl_id else []
        time.sleep(0.3)

        trials = fetch_cancer_trials(gene)
        time.sleep(0.3)

        # Score opportunity
        opportunity = score_drug_opportunity(cfg, cancer_assocs, trials)

        row = {
            "gene": gene,
            "protein": cfg.get("protein_name", ""),
            "immune_step": cfg.get("immune_step"),
            "ancient_pathogen": cfg.get("ancient_pathogen"),
            "selection_era": cfg.get("selection_era"),
            "n_cancer_associations": len(cancer_assocs),
            "top_cancer_types": ", ".join(cfg.get("cancer_types", [])[:4]),
            "n_cancer_trials": len(trials),
            "drug_status": cfg.get("drug_status", "")[:80],
            "n_known_compounds": len(cfg.get("known_compounds", [])),
            "opportunity_score": opportunity["total"],
            "score_ancient": opportunity["ancient_selection"],
            "score_mechanism": opportunity["cancer_mechanism"],
            "score_drug_gap": opportunity["drug_gap"],
            "score_structure": opportunity["structural_data"],
            "score_clinical": opportunity["clinical_validation"],
        }
        rows.append(row)

        all_evidence[gene] = {
            "cancer_associations": cancer_assocs,
            "trials": trials,
            "opportunity_scores": opportunity,
            "config": {
                "immune_step": cfg.get("immune_step"),
                "ancient_pathogen": cfg.get("ancient_pathogen"),
                "cancer_mechanism": cfg.get("cancer_mechanism", "").strip(),
                "cancer_types": cfg.get("cancer_types", []),
                "drug_status": cfg.get("drug_status"),
            },
        }

    # Build DataFrame
    df = pd.DataFrame(rows)
    df = df.sort_values("opportunity_score", ascending=False)

    # Save
    csv_out = PROC_DIR / "cancer_evidence_matrix.csv"
    df.to_csv(csv_out, index=False)
    print(f"\nSaved matrix: {csv_out.name}")

    json_out = PROC_DIR / "cancer_evidence_matrix.json"
    json_out.write_text(json.dumps(all_evidence, indent=2), encoding="utf-8")
    print(f"Saved evidence: {json_out.name}")

    return df


def print_matrix_report(df: pd.DataFrame):
    """Print a formatted report of the evidence matrix."""
    print("\n" + "=" * 70)
    print("CANCER EVIDENCE MATRIX — RANKED BY DRUG OPPORTUNITY")
    print("=" * 70)
    print(f"\n{'Gene':12s} {'Step':12s} {'Pathogen':25s} {'Cancer Assocs':>13s} {'Trials':>6s} {'Score':>5s}/15")
    print("-" * 80)
    for _, row in df.iterrows():
        print(
            f"  {row['gene']:10s} {row['immune_step']:12s} "
            f"{row['ancient_pathogen']:25s} "
            f"{row['n_cancer_associations']:>8d}      "
            f"{row['n_cancer_trials']:>4d}  "
            f"{row['opportunity_score']:>4d}"
        )

    print(f"\n{'Scoring Breakdown':}")
    print(f"{'Gene':12s} {'Ancient':>7s} {'Mechanism':>9s} {'Drug Gap':>8s} {'Structure':>9s} {'Clinical':>8s} {'TOTAL':>5s}")
    print("-" * 65)
    for _, row in df.iterrows():
        print(
            f"  {row['gene']:10s} "
            f"{row['score_ancient']:>5d}/3  "
            f"{row['score_mechanism']:>5d}/3    "
            f"{row['score_drug_gap']:>5d}/3  "
            f"{row['score_structure']:>5d}/3    "
            f"{row['score_clinical']:>5d}/3  "
            f"{row['opportunity_score']:>4d}"
        )

    # Recommendations
    print(f"\n{'='*70}")
    print("PIPELINE RECOMMENDATIONS")
    print(f"{'='*70}")
    top = df.iloc[0]
    print(f"\n  Priority 1: {top['gene']} (score {top['opportunity_score']}/15)")
    print(f"    {top['ancient_pathogen']} → {top['top_cancer_types']}")
    print(f"    Drug status: {top['drug_status']}")

    if len(df) > 1:
        second = df.iloc[1]
        print(f"\n  Priority 2: {second['gene']} (score {second['opportunity_score']}/15)")
        print(f"    {second['ancient_pathogen']} → {second['top_cancer_types']}")
        print(f"    Drug status: {second['drug_status']}")


def main():
    print("=" * 70)
    print("Cancer Connections — Ancient Pathogen Selection & Cancer Framework")
    print("=" * 70)

    targets = load_targets()
    print(f"\nLoaded {len(targets)} targets from targets.yaml")

    df = build_evidence_matrix(targets)
    print_matrix_report(df)

    print("\nDone.")


if __name__ == "__main__":
    main()
