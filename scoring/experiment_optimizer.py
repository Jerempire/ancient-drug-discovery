"""CRO experiment optimizer — rank experiments by information gain per dollar.

For each candidate at a decision point, lists possible experiments,
estimates cost and information gain, and ranks by information_gain/cost.

The key insight: use the BetaEstimator's `strength` field to identify
which dimensions have the most uncertainty (low strength) on the
highest-weight dimensions. Cheap experiments that resolve these
uncertainties have maximum expected value.

Usage:
    from scoring.experiment_optimizer import recommend_experiments
    from scoring.bayes_target import score_targets_from_yaml

    targets = load_targets()
    scorer = score_targets_from_yaml(targets)
    recommendations = recommend_experiments(scorer, targets)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .bayes_target import TargetScorer, SCORING_DIMENSIONS


@dataclass
class Experiment:
    name: str
    description: str
    estimated_cost_usd: int
    duration_weeks: int
    dimensions_updated: list[str]
    expected_confidence_gain: float  # 0-1: how much this reduces CI width
    requires: list[str] = field(default_factory=list)  # prerequisite experiments

    @property
    def cost_per_dimension(self) -> float:
        n = len(self.dimensions_updated)
        return self.estimated_cost_usd / n if n > 0 else float("inf")


# Experiment catalog — typical CRO costs for drug discovery
EXPERIMENT_CATALOG: list[Experiment] = [
    Experiment(
        name="SPR binding assay",
        description="Surface plasmon resonance — confirms direct binding to protein target",
        estimated_cost_usd=5_000,
        duration_weeks=2,
        dimensions_updated=["druggability"],
        expected_confidence_gain=0.35,
    ),
    Experiment(
        name="Cell viability assay (IC50)",
        description="Dose-response in cancer cell lines — confirms functional activity",
        estimated_cost_usd=8_000,
        duration_weeks=3,
        dimensions_updated=["tumor_relevance", "druggability"],
        expected_confidence_gain=0.40,
    ),
    Experiment(
        name="Selectivity panel (10 off-targets)",
        description="Counter-screen against related proteins — confirms target specificity",
        estimated_cost_usd=12_000,
        duration_weeks=3,
        dimensions_updated=["druggability"],
        expected_confidence_gain=0.25,
    ),
    Experiment(
        name="CRISPR KO validation",
        description="Gene knockout in cancer cell line — confirms target is essential",
        estimated_cost_usd=15_000,
        duration_weeks=6,
        dimensions_updated=["tumor_relevance", "clinical_validation"],
        expected_confidence_gain=0.50,
    ),
    Experiment(
        name="Mouse xenograft (single agent)",
        description="Tumor growth in vivo — confirms efficacy in animal model",
        estimated_cost_usd=35_000,
        duration_weeks=10,
        dimensions_updated=["tumor_relevance", "clinical_validation"],
        expected_confidence_gain=0.45,
        requires=["Cell viability assay (IC50)"],
    ),
    Experiment(
        name="Microsomal stability + CYP inhibition",
        description="Metabolic stability and drug-drug interaction risk — basic ADMET",
        estimated_cost_usd=6_000,
        duration_weeks=2,
        dimensions_updated=["druggability"],
        expected_confidence_gain=0.20,
    ),
    Experiment(
        name="Crystallography / Cryo-EM",
        description="Co-crystal structure of compound bound to target — guides optimization",
        estimated_cost_usd=25_000,
        duration_weeks=8,
        dimensions_updated=["druggability", "ip_room"],
        expected_confidence_gain=0.55,
    ),
    Experiment(
        name="Patent landscape search",
        description="Freedom-to-operate analysis by IP attorney",
        estimated_cost_usd=3_000,
        duration_weeks=2,
        dimensions_updated=["ip_room", "commercial_viability"],
        expected_confidence_gain=0.40,
    ),
    Experiment(
        name="Biomarker correlation analysis",
        description="Retrospective analysis of target expression vs clinical outcome in public datasets",
        estimated_cost_usd=2_000,
        duration_weeks=2,
        dimensions_updated=["clinical_validation", "commercial_viability"],
        expected_confidence_gain=0.30,
    ),
    Experiment(
        name="In silico ADMET prediction",
        description="Computational prediction of absorption, distribution, metabolism, excretion, toxicity",
        estimated_cost_usd=500,
        duration_weeks=1,
        dimensions_updated=["druggability"],
        expected_confidence_gain=0.15,
    ),
    Experiment(
        name="Market sizing analysis",
        description="Epidemiology + pricing analysis for target indication",
        estimated_cost_usd=1_500,
        duration_weeks=1,
        dimensions_updated=["commercial_viability"],
        expected_confidence_gain=0.35,
    ),
]


def _information_value(
    scorer: TargetScorer,
    target: str,
    experiment: Experiment,
) -> float:
    """Estimate the information value of an experiment for a target.

    Higher value when:
    - The experiment updates high-weight dimensions
    - Those dimensions currently have low strength (high uncertainty)
    - The expected confidence gain is large
    """
    total_value = 0.0

    for dim in experiment.dimensions_updated:
        est = scorer.get_estimator(target, dim)
        estimate = est.get_estimate()

        # Dimension weight (importance)
        dim_weight = scorer.DIMENSION_WEIGHTS.get(dim, 0.1)

        # Current uncertainty (inverse of strength)
        strength = estimate["strength"]
        uncertainty = 1.0 / (1.0 + strength)  # 0→1 as strength→0

        # CI width as direct measure of remaining uncertainty
        ci_width = estimate["ci_high"] - estimate["ci_low"]

        # Information value = dimension importance * current uncertainty * expected gain
        value = dim_weight * uncertainty * ci_width * experiment.expected_confidence_gain
        total_value += value

    return total_value


def recommend_experiments(
    scorer: TargetScorer,
    targets: dict,
    budget_usd: int | None = None,
    max_recommendations: int = 10,
) -> list[dict[str, Any]]:
    """Recommend experiments ranked by information gain per dollar.

    Args:
        scorer: Initialized TargetScorer with current posteriors
        targets: targets.yaml dict
        budget_usd: Optional budget constraint
        max_recommendations: Max number of recommendations

    Returns:
        List of {target, experiment, info_value, cost, ratio, reasoning}
    """
    recommendations = []

    for gene in targets:
        target_score = scorer.score_target(gene)

        for exp in EXPERIMENT_CATALOG:
            info_value = _information_value(scorer, gene, exp)
            if info_value < 0.001:
                continue  # No meaningful information gain

            ratio = info_value / exp.estimated_cost_usd * 100_000  # Normalize

            # Build reasoning
            weak_dims = []
            for dim in exp.dimensions_updated:
                est = scorer.get_estimator(gene, dim).get_estimate()
                if est["strength"] < 5:
                    weak_dims.append(f"{dim} (str={est['strength']:.1f})")

            reasoning = ""
            if weak_dims:
                reasoning = f"Resolves uncertainty in: {', '.join(weak_dims)}"

            recommendations.append({
                "target": gene,
                "experiment": exp.name,
                "description": exp.description,
                "cost_usd": exp.estimated_cost_usd,
                "duration_weeks": exp.duration_weeks,
                "info_value": round(info_value, 4),
                "value_per_dollar": round(ratio, 2),
                "dimensions": exp.dimensions_updated,
                "reasoning": reasoning,
                "target_verdict": target_score["verdict"],
                "requires": exp.requires,
            })

    # Sort by value-per-dollar ratio
    recommendations.sort(key=lambda x: x["value_per_dollar"], reverse=True)

    # Apply budget constraint if specified
    if budget_usd is not None:
        filtered = []
        spent = 0
        for rec in recommendations:
            if spent + rec["cost_usd"] <= budget_usd:
                filtered.append(rec)
                spent += rec["cost_usd"]
        recommendations = filtered

    return recommendations[:max_recommendations]


def print_experiment_report(recommendations: list[dict], budget_usd: int | None = None) -> None:
    """Print a formatted experiment recommendation report."""
    print("\n" + "=" * 70)
    print("CRO EXPERIMENT OPTIMIZER — RANKED BY INFORMATION GAIN PER DOLLAR")
    print("=" * 70)

    if budget_usd:
        print(f"\nBudget: ${budget_usd:,}")

    total_cost = sum(r["cost_usd"] for r in recommendations)
    print(f"Recommended experiments: {len(recommendations)} (total: ${total_cost:,})\n")

    for i, rec in enumerate(recommendations, 1):
        print(f"  {i}. [{rec['target']}] {rec['experiment']} — ${rec['cost_usd']:,}")
        print(f"     {rec['description']}")
        print(f"     Value/$ ratio: {rec['value_per_dollar']:.1f} | "
              f"Duration: {rec['duration_weeks']}w | "
              f"Dims: {', '.join(rec['dimensions'])}")
        if rec["reasoning"]:
            print(f"     {rec['reasoning']}")
        if rec["requires"]:
            print(f"     Requires: {', '.join(rec['requires'])}")
        print()
