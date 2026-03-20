"""
07_cancer_connections.py — Cancer-specific evidence compiler for pathogen-selected immune targets.

For each target in targets.yaml, queries:
  1. OpenTargets for cancer-specific associations and evidence
  2. ClinicalTrials.gov for active cancer trials involving the target
  3. Literature references for the pathogen→cancer connection

Builds a unified cancer evidence matrix showing:
  - Which immune step each target controls (detection/activation/killing)
  - Cancer types affected
  - Drug opportunity score (legacy 15-point + Bayesian BetaEstimator)
  - Evidence strength ranking with Crawford-Sobel tiers
  - Commercial viability assessment

Outputs:
  data/processed/cancer_evidence_matrix.csv
  data/processed/cancer_evidence_matrix.json
  data/processed/cancer_trials_<gene>.json (per target)
  data/target_priors.json (Bayesian posterior state)
  data/processed/bayesian_rankings.json
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

# Add scoring module to path
sys.path.insert(0, str(PROJECT_DIR))
from scoring.bayes_target import TargetScorer, score_targets_from_yaml
from scoring.commercial import score_commercial, CommercialScore
from scoring.evidence_tiers import classify_evidence, TIERS

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


def build_evidence_matrix(targets: dict) -> tuple[pd.DataFrame, TargetScorer]:
    """Build the unified cancer evidence matrix with Bayesian scoring."""
    print("\n" + "=" * 70)
    print("Building Cancer Evidence Matrix")
    print("=" * 70)

    # Initialize Bayesian scorer from targets.yaml (seeds priors)
    scorer = score_targets_from_yaml(targets)

    rows = []
    all_evidence = {}
    commercial_scores: dict[str, CommercialScore] = {}

    for gene, cfg in targets.items():
        ensembl_id = cfg.get("ensembl_id", "")

        # Fetch cancer evidence
        cancer_assocs = fetch_cancer_associations(gene, ensembl_id) if ensembl_id else []
        time.sleep(0.3)

        trials = fetch_cancer_trials(gene)
        time.sleep(0.3)

        # --- Layer OpenTargets evidence into Bayesian scorer ---
        for assoc in cancer_assocs:
            ot_score = assoc.get("score", 0)
            if ot_score > 0:
                tier = classify_evidence(source_type="cohort", text="OpenTargets disease association")
                scorer.update(
                    gene, "tumor_relevance", "opentargets",
                    probability=min(0.95, 0.4 + ot_score * 0.5),
                    confidence=min(0.8, ot_score),
                    evidence_tier=tier,
                )

        # --- Layer clinical trials into Bayesian scorer ---
        active_statuses = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"}
        active_trials = [t for t in trials if t.get("status") in active_statuses]
        if active_trials:
            for trial in active_trials[:5]:  # Cap updates
                phases = trial.get("phase", [])
                phase_str = " ".join(phases).lower()
                tier = classify_evidence(text=phase_str)
                prob = 0.55
                if "phase 3" in phase_str or "phase iii" in phase_str:
                    prob = 0.80
                elif "phase 2" in phase_str or "phase ii" in phase_str:
                    prob = 0.70
                scorer.update(
                    gene, "clinical_validation", "clinical_trial",
                    probability=prob, confidence=0.60, evidence_tier=tier,
                )

        # --- Commercial scoring ---
        commercial = score_commercial(cfg, trials)
        commercial_scores[gene] = commercial

        # Feed commercial composite into Bayesian scorer
        scorer.update(
            gene, "commercial_viability", "expert_curation",
            probability=commercial.composite,
            confidence=0.55,
        )

        # Legacy scoring (kept for comparison)
        opportunity = score_drug_opportunity(cfg, cancer_assocs, trials)

        # Bayesian score
        bayes_result = scorer.score_target(gene)

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
            # Legacy score
            "opportunity_score": opportunity["total"],
            "score_ancient": opportunity["ancient_selection"],
            "score_mechanism": opportunity["cancer_mechanism"],
            "score_drug_gap": opportunity["drug_gap"],
            "score_structure": opportunity["structural_data"],
            "score_clinical": opportunity["clinical_validation"],
            # Bayesian score
            "bayes_composite": bayes_result["composite_score"],
            "bayes_confidence": bayes_result["confidence"],
            "bayes_strength": bayes_result["total_strength"],
            "bayes_verdict": bayes_result["verdict"],
            # Commercial score
            "commercial_composite": commercial.composite,
            "commercial_unmet_need": commercial.unmet_need,
            "commercial_market_size": commercial.market_size,
            "commercial_ip_room": commercial.ip_room,
        }
        rows.append(row)

        all_evidence[gene] = {
            "cancer_associations": cancer_assocs,
            "trials": trials,
            "opportunity_scores": opportunity,
            "bayesian_scores": bayes_result,
            "commercial_scores": commercial.to_dict(),
            "config": {
                "immune_step": cfg.get("immune_step"),
                "ancient_pathogen": cfg.get("ancient_pathogen"),
                "cancer_mechanism": cfg.get("cancer_mechanism", "").strip(),
                "cancer_types": cfg.get("cancer_types", []),
                "drug_status": cfg.get("drug_status"),
            },
        }

    # Build DataFrame — sort by Bayesian composite (primary), legacy (secondary)
    df = pd.DataFrame(rows)
    df = df.sort_values(["bayes_composite", "opportunity_score"], ascending=[False, False])

    # Save
    csv_out = PROC_DIR / "cancer_evidence_matrix.csv"
    df.to_csv(csv_out, index=False)
    print(f"\nSaved matrix: {csv_out.name}")

    json_out = PROC_DIR / "cancer_evidence_matrix.json"
    json_out.write_text(json.dumps(all_evidence, indent=2), encoding="utf-8")
    print(f"Saved evidence: {json_out.name}")

    # Save Bayesian priors for persistence
    scorer.save()
    print(f"Saved Bayesian priors: target_priors.json")

    # Save ranked Bayesian results
    rankings = scorer.rank_targets()
    rank_out = PROC_DIR / "bayesian_rankings.json"
    rank_out.write_text(json.dumps(rankings, indent=2), encoding="utf-8")
    print(f"Saved Bayesian rankings: {rank_out.name}")

    return df, scorer


def print_matrix_report(df: pd.DataFrame, scorer: TargetScorer | None = None):
    """Print a formatted report of the evidence matrix."""
    # --- Legacy scoring ---
    print("\n" + "=" * 70)
    print("CANCER EVIDENCE MATRIX — LEGACY SCORING (15-point)")
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

    # --- Bayesian scoring ---
    print(f"\n{'='*70}")
    print("BAYESIAN TARGET RANKING (BetaEstimator + Crawford-Sobel tiers)")
    print(f"{'='*70}")
    print(f"\n{'Gene':12s} {'Composite':>9s} {'Confidence':>10s} {'Strength':>8s}  {'Verdict'}")
    print("-" * 65)
    for _, row in df.iterrows():
        print(
            f"  {row['gene']:10s} "
            f"{row['bayes_composite']:>8.1%}  "
            f"{row['bayes_confidence']:>8.1%}   "
            f"{row['bayes_strength']:>7.1f}   "
            f"{row['bayes_verdict']}"
        )

    # --- Per-target Bayesian detail ---
    if scorer:
        print(f"\n{'='*70}")
        print("PER-DIMENSION BAYESIAN BREAKDOWN")
        print(f"{'='*70}")
        for _, row in df.iterrows():
            gene = row["gene"]
            result = scorer.score_target(gene)
            print(f"\n  {gene}:")
            for dim, est in result["dimensions"].items():
                ci_w = est["ci_high"] - est["ci_low"]
                bar = "#" * int(est["mean"] * 20)
                print(
                    f"    {dim:25s} {est['mean']:5.1%} "
                    f"[{est['ci_low']:.0%}-{est['ci_high']:.0%}] "
                    f"str={est['strength']:4.1f}  {bar}"
                )

    # --- Commercial scoring ---
    print(f"\n{'='*70}")
    print("COMMERCIAL VIABILITY")
    print(f"{'='*70}")
    print(f"\n{'Gene':12s} {'Composite':>9s} {'Unmet':>6s} {'Market':>7s} {'IP':>5s}")
    print("-" * 50)
    for _, row in df.iterrows():
        print(
            f"  {row['gene']:10s} "
            f"{row['commercial_composite']:>8.1%}  "
            f"{row['commercial_unmet_need']:>5.0%} "
            f"{row['commercial_market_size']:>6.0%}  "
            f"{row['commercial_ip_room']:>4.0%}"
        )

    # --- Recommendations ---
    print(f"\n{'='*70}")
    print("PIPELINE RECOMMENDATIONS (Bayesian-ranked)")
    print(f"{'='*70}")
    top = df.iloc[0]
    print(f"\n  Priority 1: {top['gene']} "
          f"(Bayes={top['bayes_composite']:.0%}, Legacy={top['opportunity_score']}/15)")
    print(f"    {top['ancient_pathogen']} -> {top['top_cancer_types']}")
    print(f"    Verdict: {top['bayes_verdict']}")
    print(f"    Drug status: {top['drug_status']}")

    if len(df) > 1:
        second = df.iloc[1]
        print(f"\n  Priority 2: {second['gene']} "
              f"(Bayes={second['bayes_composite']:.0%}, Legacy={second['opportunity_score']}/15)")
        print(f"    {second['ancient_pathogen']} -> {second['top_cancer_types']}")
        print(f"    Verdict: {second['bayes_verdict']}")
        print(f"    Drug status: {second['drug_status']}")


def main():
    print("=" * 70)
    print("Cancer Connections — Ancient Pathogen Selection & Cancer Framework")
    print("=" * 70)

    targets = load_targets()
    print(f"\nLoaded {len(targets)} targets from targets.yaml")

    df, scorer = build_evidence_matrix(targets)
    print_matrix_report(df, scorer)

    print("\nDone.")


if __name__ == "__main__":
    main()
