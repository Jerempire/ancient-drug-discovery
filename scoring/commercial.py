"""Commercial viability scoring for drug targets.

Evaluates the business case for pursuing a target, independent of biology.
This is the user's unique edge — pharmacy/reimbursement expertise that
typical AI drug discovery teams lack.

Six sub-dimensions:
1. Unmet need — existing treatments, approval history
2. Market size — prevalence x price ceiling
3. IP room — patent landscape density
4. Competitive pipeline — active trials in same target space
5. Reimbursement potential — payer willingness, indication tier
6. Biomarker strategy — companion diagnostic feasibility
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class CommercialScore:
    unmet_need: float = 0.0        # 0-1
    market_size: float = 0.0       # 0-1
    ip_room: float = 0.0           # 0-1
    competitive_pipeline: float = 0.0  # 0-1 (higher = less competition = better)
    reimbursement: float = 0.0     # 0-1
    biomarker_strategy: float = 0.0  # 0-1
    notes: list[str] = field(default_factory=list)

    @property
    def composite(self) -> float:
        """Weighted composite commercial score."""
        weights = {
            "unmet_need": 0.25,
            "market_size": 0.20,
            "ip_room": 0.20,
            "competitive_pipeline": 0.15,
            "reimbursement": 0.10,
            "biomarker_strategy": 0.10,
        }
        total = sum(
            getattr(self, dim) * w for dim, w in weights.items()
        )
        return round(total, 4)

    def to_dict(self) -> dict[str, Any]:
        return {
            "unmet_need": round(self.unmet_need, 3),
            "market_size": round(self.market_size, 3),
            "ip_room": round(self.ip_room, 3),
            "competitive_pipeline": round(self.competitive_pipeline, 3),
            "reimbursement": round(self.reimbursement, 3),
            "biomarker_strategy": round(self.biomarker_strategy, 3),
            "composite": self.composite,
            "notes": self.notes,
        }


# Reimbursement tier mapping
# Oncology biologics command premium pricing; generics face price compression
INDICATION_REIMBURSEMENT_TIERS: dict[str, float] = {
    "oncology_biologic": 0.95,
    "oncology_small_molecule": 0.85,
    "rare_disease": 0.90,
    "autoimmune_biologic": 0.80,
    "autoimmune_small_molecule": 0.65,
    "infectious_disease": 0.50,
    "generic_small_molecule": 0.30,
}


def score_commercial(target_cfg: dict, trials: list[dict] | None = None) -> CommercialScore:
    """Score commercial viability from targets.yaml config + trial data.

    Args:
        target_cfg: Single target entry from targets.yaml
        trials: ClinicalTrials.gov results (from 07_cancer_connections.py)
    """
    trials = trials or []
    score = CommercialScore()

    # --- 1. Unmet need ---
    drug_status = target_cfg.get("drug_status", "").lower()
    compounds = target_cfg.get("known_compounds", [])

    if "no drugs" in drug_status and not compounds:
        score.unmet_need = 0.90
        score.notes.append("Completely undrugged target — maximum unmet need")
    elif "no drugs" in drug_status or "undrugged" in drug_status:
        score.unmet_need = 0.85
    elif "repurpos" in drug_status:
        score.unmet_need = 0.60
        score.notes.append("Existing drug being repurposed — moderate unmet need")
    elif len(compounds) == 1:
        score.unmet_need = 0.55
    else:
        score.unmet_need = 0.35
        score.notes.append("Multiple known compounds — limited unmet need")

    # --- 2. Market size (prevalence proxy from cancer types) ---
    cancer_types = target_cfg.get("cancer_types", [])
    # Large-market cancers
    large_market = {"breast", "colorectal", "lung", "prostate", "melanoma", "hepatocellular"}
    n_large = sum(1 for ct in cancer_types if ct in large_market)
    n_total = len(cancer_types)

    if n_large >= 2:
        score.market_size = 0.85
    elif n_large >= 1:
        score.market_size = 0.70
    elif n_total >= 2:
        score.market_size = 0.55
    elif n_total >= 1:
        score.market_size = 0.40
    else:
        score.market_size = 0.20

    # Autoimmune indications add market breadth
    autoimmune = target_cfg.get("autoimmune_diseases", [])
    if autoimmune:
        score.market_size = min(0.95, score.market_size + 0.10 * len(autoimmune))
        score.notes.append(f"Dual oncology+autoimmune indication ({len(autoimmune)} autoimmune)")

    # --- 3. IP room ---
    if "no drugs" in drug_status and not compounds:
        score.ip_room = 0.90
        score.notes.append("Green field — no existing IP to navigate")
    elif not compounds:
        score.ip_room = 0.80
    elif len(compounds) == 1:
        # One compound = need to design around it
        score.ip_room = 0.55
    else:
        # Crowded space
        score.ip_room = 0.30

    # Check for recent patents (2022+ compounds suggest active IP landscape)
    for cpd in compounds:
        ref = cpd.get("reference", "")
        if any(str(y) in ref for y in range(2022, 2027)):
            score.ip_room = max(0.20, score.ip_room - 0.15)
            score.notes.append(f"Recent compound: {cpd.get('name', '?')} — active IP landscape")
            break

    # --- 4. Competitive pipeline ---
    active_statuses = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"}
    active_trials = [t for t in trials if t.get("status") in active_statuses]

    if not trials:
        # No trials = either greenfield (good) or no validation (bad)
        score.competitive_pipeline = 0.70
        score.notes.append("No competing clinical trials found")
    elif len(active_trials) == 0:
        score.competitive_pipeline = 0.75
    elif len(active_trials) <= 2:
        score.competitive_pipeline = 0.55
        score.notes.append(f"{len(active_trials)} active competing trials")
    elif len(active_trials) <= 5:
        score.competitive_pipeline = 0.40
    else:
        score.competitive_pipeline = 0.20
        score.notes.append(f"{len(active_trials)} active trials — crowded space")

    # --- 5. Reimbursement potential ---
    immune_step = target_cfg.get("immune_step", "")
    mechanism = target_cfg.get("cancer_mechanism", "").lower()

    # Oncology biologics / novel mechanisms get premium reimbursement
    if "checkpoint" in mechanism or "immunotherapy" in mechanism or immune_step == "detection":
        score.reimbursement = 0.85
        score.notes.append("Immuno-oncology mechanism — premium reimbursement tier")
    elif cancer_types:
        score.reimbursement = 0.75  # Any oncology indication
    elif autoimmune:
        score.reimbursement = 0.65
    else:
        score.reimbursement = 0.50

    # --- 6. Biomarker strategy ---
    key_variants = target_cfg.get("key_variants", [])
    if key_variants:
        # Known variants = companion diagnostic potential
        score.biomarker_strategy = min(0.85, 0.40 + 0.15 * len(key_variants))
        score.notes.append(f"{len(key_variants)} known variants — companion Dx feasible")
    else:
        score.biomarker_strategy = 0.30

    # Selection era as population-level biomarker
    if target_cfg.get("selection_era") in ("medieval", "neolithic_onward"):
        score.biomarker_strategy = min(0.90, score.biomarker_strategy + 0.10)

    return score
