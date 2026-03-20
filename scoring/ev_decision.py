"""Expected Value decision layer and portfolio selection.

Combines pipeline approval probability, commercial scoring, and market
estimates into a single EV number per target. Portfolio selection picks
the best *set* of targets within a budget constraint.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .bayes_target import TargetScorer
from .commercial import CommercialScore
from .markov_pipeline import CandidateState

# Rough oncology market caps by indication type (USD)
MARKET_CAP_ESTIMATES: dict[str, int] = {
    "large_oncology": 2_000_000_000,    # $2B (breast, lung, colorectal)
    "medium_oncology": 500_000_000,     # $500M (prostate, melanoma)
    "rare_oncology": 200_000_000,       # $200M (niche indications)
    "autoimmune": 1_000_000_000,        # $1B
}

# Per-state experiment cost estimates (USD) from current state to approval
STATE_COSTS: dict[str, int] = {
    "discovery": 5_000_000,
    "lead_optimization": 3_000_000,
    "IND_enabling": 2_000_000,
    "phase1": 10_000_000,
    "phase2": 30_000_000,
    "phase3": 100_000_000,
    "NDA": 5_000_000,
    "approved": 0,
    "rejected": 0,
}


@dataclass
class EVEstimate:
    target: str
    p_approval: float
    payoff_usd: float
    total_experiment_cost: float
    ev_usd: float
    roi: float
    risk_tier: str
    components: dict[str, Any] = field(default_factory=dict)


@dataclass
class Portfolio:
    targets: list[EVEstimate]
    total_ev: float
    total_cost: float
    portfolio_roi: float
    risk_distribution: dict[str, int] = field(default_factory=dict)


def _estimate_remaining_cost(current_state: str) -> float:
    """Sum costs from current state through approval."""
    from .markov_pipeline import STATES
    try:
        idx = STATES.index(current_state)
    except ValueError:
        idx = 0
    return float(sum(
        STATE_COSTS.get(s, 0)
        for s in STATES[idx:]
        if s not in ("approved", "rejected")
    ))


def _assign_risk_tier(p_approval: float, has_compounds: bool) -> str:
    if p_approval >= 0.05 and has_compounds:
        return "safe"
    elif p_approval >= 0.01:
        return "medium"
    else:
        return "high_risk"


def compute_ev(
    candidate: CandidateState,
    scorer: TargetScorer,
    commercial: CommercialScore,
    market_type: str = "medium_oncology",
    market_estimates: dict[str, int] | None = None,
) -> EVEstimate:
    """Compute expected value for a single candidate.

    Args:
        candidate: Pipeline state with p_approval calculation
        scorer: TargetScorer with Bayesian estimates (used to check compound existence)
        commercial: CommercialScore with composite payoff proxy
        market_type: Key into MARKET_CAP_ESTIMATES
        market_estimates: Override default market cap dict
    """
    estimates = market_estimates or MARKET_CAP_ESTIMATES
    market_cap = estimates.get(market_type, 500_000_000)

    p_approval = candidate.p_approval()
    payoff_usd = market_cap * commercial.composite
    total_cost = _estimate_remaining_cost(candidate.current_state)
    ev_usd = p_approval * payoff_usd - total_cost

    has_compounds = candidate.current_state not in ("discovery",)
    risk_tier = _assign_risk_tier(p_approval, has_compounds)
    roi = ev_usd / total_cost if total_cost > 0 else 0.0

    return EVEstimate(
        target=candidate.target,
        p_approval=p_approval,
        payoff_usd=round(payoff_usd, 2),
        total_experiment_cost=total_cost,
        ev_usd=round(ev_usd, 2),
        roi=round(roi, 4),
        risk_tier=risk_tier,
        components={
            "candidate": candidate.candidate_name,
            "current_state": candidate.current_state,
            "market_type": market_type,
            "market_cap": market_cap,
            "commercial_composite": commercial.composite,
        },
    )


def rank_targets_by_ev(
    candidates: list[CandidateState],
    scorer: TargetScorer,
    commercial_scores: dict[str, CommercialScore],
    market_type: str = "medium_oncology",
    market_estimates: dict[str, int] | None = None,
) -> list[EVEstimate]:
    """Rank all candidates by EV descending."""
    evs = []
    for c in candidates:
        commercial = commercial_scores.get(c.target)
        if commercial is None:
            continue
        evs.append(compute_ev(c, scorer, commercial, market_type, market_estimates))
    return sorted(evs, key=lambda e: e.ev_usd, reverse=True)


def select_portfolio(
    ev_estimates: list[EVEstimate],
    budget: float,
    min_targets: int = 3,
    max_targets: int = 6,
) -> Portfolio:
    """Select best portfolio of targets within budget.

    Greedy algorithm with diversification constraint:
    1. Must include at least 1 target from each available risk tier
    2. Greedily add highest-ROI targets that fit within budget
    """
    by_tier: dict[str, list[EVEstimate]] = {"safe": [], "medium": [], "high_risk": []}
    for ev in ev_estimates:
        by_tier.setdefault(ev.risk_tier, []).append(ev)
    for tier_list in by_tier.values():
        tier_list.sort(key=lambda e: e.roi, reverse=True)

    selected: list[EVEstimate] = []
    selected_targets: set[str] = set()
    remaining_budget = budget

    # Phase 1: one from each available tier (diversification)
    for tier in ("safe", "medium", "high_risk"):
        for ev in by_tier.get(tier, []):
            if ev.total_experiment_cost <= remaining_budget and ev.target not in selected_targets:
                selected.append(ev)
                selected_targets.add(ev.target)
                remaining_budget -= ev.total_experiment_cost
                break

    # Phase 2: greedily fill by ROI
    all_by_roi = sorted(ev_estimates, key=lambda e: e.roi, reverse=True)
    for ev in all_by_roi:
        if len(selected) >= max_targets:
            break
        if ev.target in selected_targets:
            continue
        if ev.total_experiment_cost <= remaining_budget:
            selected.append(ev)
            selected_targets.add(ev.target)
            remaining_budget -= ev.total_experiment_cost

    total_ev = sum(e.ev_usd for e in selected)
    total_cost = sum(e.total_experiment_cost for e in selected)
    risk_dist = {}
    for e in selected:
        risk_dist[e.risk_tier] = risk_dist.get(e.risk_tier, 0) + 1

    return Portfolio(
        targets=selected,
        total_ev=round(total_ev, 2),
        total_cost=round(total_cost, 2),
        portfolio_roi=round(total_ev / total_cost, 4) if total_cost > 0 else 0.0,
        risk_distribution=risk_dist,
    )
