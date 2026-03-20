"""Bayesian target prioritization using Beta distribution updating.

Adapted from polymarket-edge/shared/bayes.py for drug discovery.

Instead of market probabilities, we estimate P(target is viable) across
multiple scoring dimensions. Each piece of evidence (GWAS hit, clinical trial,
structural data, etc.) updates the posterior via source-weighted pseudo-observations.

Key differences from the market version:
- Source weights tuned for biomedical evidence types (GWAS, literature, in-silico)
- Time decay in years, not hours (drug development moves slower)
- Evidence tiers (Crawford-Sobel) modulate update strength
- Multi-axis scoring: each target has separate estimators per dimension
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .evidence_tiers import EvidenceTier, classify_evidence, tier_weight

PRIORS_FILE = Path(__file__).parent.parent / "data" / "target_priors.json"

# Source-specific reliability weights for drug discovery evidence
DEFAULT_SOURCE_WEIGHTS: dict[str, float] = {
    # Tier 1: validated experimental data
    "gwas": 0.80,
    "rct": 0.90,
    "crispr_screen": 0.85,
    "clinical_trial": 0.75,
    # Tier 2: computational / curated
    "opentargets": 0.60,
    "literature": 0.50,
    "drugbank": 0.55,
    "in_silico_docking": 0.30,
    "esm2_score": 0.35,
    "alphafold": 0.40,
    # Tier 3: manual / heuristic
    "expert_curation": 0.45,
    "targets_yaml": 0.40,
    "press_release": 0.20,
}

# Scoring dimensions for drug target viability
SCORING_DIMENSIONS = [
    "tumor_relevance",       # Is the target involved in cancer?
    "ancient_selection",     # How strong is the evolutionary evidence?
    "druggability",          # Can we make a drug against it?
    "clinical_validation",   # Any clinical evidence (trials, biomarkers)?
    "ip_room",               # Patent landscape — room for new IP?
    "commercial_viability",  # Market size, reimbursement, competitive landscape
]


class TargetBetaEstimator:
    """Beta distribution estimator for a single target scoring dimension.

    Starts with Beta(1,1) = uniform prior (no information).
    Each evidence update shifts the distribution based on:
    - Source reliability weight (GWAS > in-silico docking)
    - Evidence tier (Crawford-Sobel: Phase 3 > press release)
    - Signal confidence (how directly relevant is this evidence?)
    - Time decay (old studies count less)
    """

    def __init__(
        self,
        alpha: float = 1.0,
        beta: float = 1.0,
        source_weights: dict[str, float] | None = None,
        decay_years: float = 10.0,
    ):
        self.alpha = alpha
        self.beta = beta
        self.source_weights = source_weights or DEFAULT_SOURCE_WEIGHTS
        self.decay_years = decay_years
        self.signal_history: list[dict] = []
        self.last_updated: str | None = None

    def update(
        self,
        source: str,
        probability: float,
        confidence: float = 0.5,
        evidence_tier: EvidenceTier | None = None,
        year: int | None = None,
    ) -> None:
        """Bayesian soft update from an evidence signal.

        Args:
            source: Evidence source type (key in source_weights)
            probability: How strongly this evidence supports viability (0-1)
            confidence: How directly relevant this evidence is (0-1)
            evidence_tier: Crawford-Sobel tier (modulates update strength)
            year: Publication year (for time decay)
        """
        probability = max(0.001, min(0.999, probability))
        confidence = max(0.0, min(1.0, confidence))

        source_weight = self.source_weights.get(source, 0.3)

        # Apply evidence tier weight if provided
        tier_w = 1.0
        if evidence_tier:
            tier_w = tier_weight(evidence_tier, year=year)

        # Scale factor: controls pseudo-observations per update
        # 3.0 is more conservative than polymarket's 5.0 — drug evidence is sparser
        scale = 3.0
        n_obs = source_weight * confidence * tier_w * scale

        self.alpha += n_obs * probability
        self.beta += n_obs * (1 - probability)

        now = datetime.now(timezone.utc).isoformat()
        self.last_updated = now
        self.signal_history.append({
            "source": source,
            "probability": round(probability, 4),
            "confidence": round(confidence, 4),
            "tier": evidence_tier.tier if evidence_tier else None,
            "tier_weight": round(tier_w, 3),
            "n_obs": round(n_obs, 3),
            "year": year,
            "timestamp": now,
        })

    def apply_decay(self) -> None:
        """Apply time-decay to move posterior toward prior.

        Uses decay_years instead of decay_hours — drug development timescales.
        """
        if not self.last_updated or self.decay_years <= 0:
            return
        try:
            last_dt = datetime.fromisoformat(self.last_updated)
        except (ValueError, TypeError):
            return

        years_elapsed = (datetime.now(timezone.utc) - last_dt).total_seconds() / (3600 * 24 * 365.25)
        if years_elapsed < 0.1:  # No decay within ~5 weeks
            return

        decay_factor = math.exp(-years_elapsed / self.decay_years)
        prior_a, prior_b = 1.0, 1.0
        self.alpha = prior_a + (self.alpha - prior_a) * decay_factor
        self.beta = prior_b + (self.beta - prior_b) * decay_factor

    def get_estimate(self) -> dict[str, Any]:
        """Return posterior estimate with credible interval."""
        a, b = self.alpha, self.beta
        mean = a / (a + b)

        if a > 1 and b > 1:
            mode = (a - 1) / (a + b - 2)
        else:
            mode = mean

        var = (a * b) / ((a + b) ** 2 * (a + b + 1))
        std = math.sqrt(var)
        ci_low = max(0.0, mean - 1.96 * std)
        ci_high = min(1.0, mean + 1.96 * std)

        return {
            "mean": round(mean, 4),
            "mode": round(mode, 4),
            "ci_low": round(ci_low, 4),
            "ci_high": round(ci_high, 4),
            "alpha": round(a, 3),
            "beta": round(b, 3),
            "n_updates": len(self.signal_history),
            "strength": round(a + b, 1),
        }

    def to_dict(self) -> dict:
        return {
            "alpha": round(self.alpha, 4),
            "beta": round(self.beta, 4),
            "last_updated": self.last_updated,
            "signal_history": self.signal_history[-30:],
        }

    @classmethod
    def from_dict(cls, data: dict, **kwargs) -> TargetBetaEstimator:
        est = cls(
            alpha=data.get("alpha", 1.0),
            beta=data.get("beta", 1.0),
            **kwargs,
        )
        est.last_updated = data.get("last_updated")
        est.signal_history = data.get("signal_history", [])
        return est


class TargetScorer:
    """Multi-dimensional Bayesian scorer for drug targets.

    Each target gets one TargetBetaEstimator per scoring dimension.
    Evidence updates flow to the appropriate dimension(s).
    Final ranking uses posterior means weighted by dimension importance.
    """

    # Dimension importance weights (sum to 1.0)
    DIMENSION_WEIGHTS: dict[str, float] = {
        "tumor_relevance": 0.25,
        "ancient_selection": 0.15,
        "druggability": 0.20,
        "clinical_validation": 0.20,
        "ip_room": 0.10,
        "commercial_viability": 0.10,
    }

    def __init__(self, source_weights: dict[str, float] | None = None):
        self.source_weights = source_weights or DEFAULT_SOURCE_WEIGHTS
        self._targets: dict[str, dict[str, TargetBetaEstimator]] = {}

    def get_estimator(self, target: str, dimension: str) -> TargetBetaEstimator:
        """Get or create estimator for a target+dimension pair."""
        if target not in self._targets:
            self._targets[target] = {}
        if dimension not in self._targets[target]:
            self._targets[target][dimension] = TargetBetaEstimator(
                source_weights=self.source_weights,
            )
        return self._targets[target][dimension]

    def update(
        self,
        target: str,
        dimension: str,
        source: str,
        probability: float,
        confidence: float = 0.5,
        evidence_tier: EvidenceTier | None = None,
        year: int | None = None,
    ) -> None:
        """Update a single dimension for a target."""
        est = self.get_estimator(target, dimension)
        est.update(source, probability, confidence, evidence_tier, year)

    def score_target(self, target: str, registry: "SignalRegistry | None" = None) -> dict[str, Any]:
        """Compute weighted composite score for a target.

        If registry is provided, applies a correlation penalty to avoid
        double-counting correlated dimensions.
        """
        dimensions = {}
        weighted_sum = 0.0
        total_weight = 0.0
        total_strength = 0.0

        for dim in SCORING_DIMENSIONS:
            est = self.get_estimator(target, dim)
            estimate = est.get_estimate()
            dimensions[dim] = estimate
            w = self.DIMENSION_WEIGHTS.get(dim, 0.1)
            weighted_sum += estimate["mean"] * w
            total_weight += w
            total_strength += estimate["strength"]

        composite = weighted_sum / total_weight if total_weight > 0 else 0.5

        # Correlation penalty: subtract pairwise rho * w_a * w_b * mean_a * mean_b
        correlation_penalty = 0.0
        if registry is not None:
            seen = set()
            for (dim_a, dim_b), rho in registry.correlation_matrix().items():
                pair = tuple(sorted([dim_a, dim_b]))
                if pair in seen:
                    continue
                seen.add(pair)
                if dim_a in dimensions and dim_b in dimensions:
                    w_a = self.DIMENSION_WEIGHTS.get(dim_a, 0.1) / total_weight if total_weight > 0 else 0
                    w_b = self.DIMENSION_WEIGHTS.get(dim_b, 0.1) / total_weight if total_weight > 0 else 0
                    mean_a = dimensions[dim_a]["mean"]
                    mean_b = dimensions[dim_b]["mean"]
                    correlation_penalty += rho * w_a * w_b * mean_a * mean_b
            composite = max(0.0, min(1.0, composite - correlation_penalty))

        # Confidence: average CI width across dimensions
        ci_widths = [d["ci_high"] - d["ci_low"] for d in dimensions.values()]
        avg_ci = sum(ci_widths) / len(ci_widths) if ci_widths else 1.0

        result = {
            "target": target,
            "composite_score": round(composite, 4),
            "confidence": round(1.0 - avg_ci, 4),  # Higher = more confident
            "total_strength": round(total_strength, 1),
            "dimensions": dimensions,
            "verdict": _verdict(composite, avg_ci),
        }
        if registry is not None:
            result["correlation_penalty"] = round(correlation_penalty, 6)
        return result

    def rank_targets(self) -> list[dict[str, Any]]:
        """Rank all targets by composite score."""
        scores = [self.score_target(t) for t in self._targets]
        return sorted(scores, key=lambda x: x["composite_score"], reverse=True)

    def save(self, path: Path | None = None) -> None:
        """Persist all estimators to JSON."""
        path = path or PRIORS_FILE
        store = {}
        for target, dims in self._targets.items():
            store[target] = {dim: est.to_dict() for dim, est in dims.items()}
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(store, indent=2), encoding="utf-8")

    def load(self, path: Path | None = None) -> None:
        """Load estimators from JSON."""
        path = path or PRIORS_FILE
        if not path.exists():
            return
        try:
            store = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return
        for target, dims in store.items():
            for dim, data in dims.items():
                est = TargetBetaEstimator.from_dict(
                    data, source_weights=self.source_weights,
                )
                if target not in self._targets:
                    self._targets[target] = {}
                self._targets[target][dim] = est


def _verdict(composite: float, ci_width: float) -> str:
    """Human-readable verdict based on score and confidence."""
    if composite >= 0.7 and ci_width < 0.3:
        return "STRONG GO"
    elif composite >= 0.6:
        return "GO" if ci_width < 0.4 else "LEAN GO (low confidence)"
    elif composite >= 0.4:
        return "MAYBE" if ci_width < 0.4 else "INSUFFICIENT EVIDENCE"
    else:
        return "NO-GO" if ci_width < 0.4 else "INSUFFICIENT EVIDENCE"


def score_targets_from_yaml(targets: dict) -> TargetScorer:
    """Initialize a TargetScorer and populate from targets.yaml data.

    This converts the existing heuristic evidence in targets.yaml into
    Bayesian priors, then layers on OpenTargets + ClinicalTrials.gov data.
    """
    scorer = TargetScorer()

    for gene, cfg in targets.items():
        # --- Ancient selection ---
        selection = cfg.get("selection_event", "").lower()
        if "+40%" in selection or "strongest" in selection:
            prob, conf = 0.90, 0.85
        elif "near-fixation" in selection or "400m" in selection:
            prob, conf = 0.85, 0.80
        elif "susceptibility" in selection or "polymorphism" in selection:
            prob, conf = 0.65, 0.60
        else:
            prob, conf = 0.50, 0.40

        tier = classify_evidence(text=selection)
        scorer.update(gene, "ancient_selection", "targets_yaml", prob, conf, tier)

        # --- Tumor relevance ---
        cancer_types = cfg.get("cancer_types", [])
        mechanism = cfg.get("cancer_mechanism", "")
        n_types = len(cancer_types)
        if n_types >= 3:
            prob = 0.80
        elif n_types >= 2:
            prob = 0.65
        else:
            prob = 0.45
        conf = min(0.7, 0.3 + 0.1 * n_types)
        tier = classify_evidence(text=mechanism)
        scorer.update(gene, "tumor_relevance", "targets_yaml", prob, conf, tier)

        # --- Druggability ---
        drug_status = cfg.get("drug_status", "").lower()
        compounds = cfg.get("known_compounds", [])
        structure = cfg.get("structure", {})
        pdb = structure.get("pdb_id")
        plddt = structure.get("alphafold_plddt")

        # Structure-based druggability
        if pdb:
            scorer.update(gene, "druggability", "targets_yaml", 0.75, 0.70,
                          classify_evidence(source_type="rct"))
        elif plddt and plddt > 90:
            scorer.update(gene, "druggability", "alphafold", 0.70, 0.55)
        elif plddt and plddt > 70:
            scorer.update(gene, "druggability", "alphafold", 0.55, 0.40)
        else:
            scorer.update(gene, "druggability", "alphafold", 0.40, 0.30)

        # Existing compounds boost druggability
        if compounds:
            for cpd in compounds:
                ref = cpd.get("reference", "")
                tier = classify_evidence(text=ref)
                scorer.update(gene, "druggability", "drugbank", 0.70, 0.60, tier)

        # --- Clinical validation ---
        # Will be updated later with ClinicalTrials.gov data
        # For now, seed from compounds with clinical references
        if any("phase" in c.get("reference", "").lower() for c in compounds):
            scorer.update(gene, "clinical_validation", "targets_yaml", 0.65, 0.50,
                          classify_evidence(text="phase II trial"))

        # --- IP room ---
        if "no drugs" in drug_status and not compounds:
            scorer.update(gene, "ip_room", "targets_yaml", 0.85, 0.60)
        elif "no drugs" in drug_status or "undrugged" in drug_status:
            scorer.update(gene, "ip_room", "targets_yaml", 0.80, 0.55)
        elif "repurpos" in drug_status:
            scorer.update(gene, "ip_room", "targets_yaml", 0.50, 0.50)
        elif len(compounds) >= 2:
            scorer.update(gene, "ip_room", "targets_yaml", 0.35, 0.45)
        else:
            scorer.update(gene, "ip_room", "targets_yaml", 0.55, 0.40)

        # --- Commercial viability (basic seed, A5 will extend) ---
        scorer.update(gene, "commercial_viability", "targets_yaml", 0.50, 0.30)

    return scorer
