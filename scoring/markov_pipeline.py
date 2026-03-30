"""Markov state machine for drug development pipeline progression.

Ported from scenario-forecaster/models/engine.py.

Each drug candidate tracks its current state in the development pipeline.
Transition probabilities are initialized from BIO/Informa Phase Transition
data and adjusted by target-specific evidence.

States: discovery → lead_optimization → IND_enabling → phase1 → phase2 →
        phase3 → NDA → approved
        (with absorbing state: rejected)

Base rates source: BIO/Informa/QLS "Clinical Development Success Rates
2011-2020" and Thomas et al. Nature Reviews Drug Discovery (2016).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import yaml


# BIO/Informa historical phase transition success rates (oncology-specific where available)
# These are P(advance to next phase | entered current phase)
BASE_TRANSITIONS: dict[str, dict[str, float]] = {
    "discovery": {
        "lead_optimization": 0.30,
        "rejected": 0.70,
    },
    "lead_optimization": {
        "IND_enabling": 0.50,
        "rejected": 0.50,
    },
    "IND_enabling": {
        "phase1": 0.69,  # IND clearance rate
        "rejected": 0.31,
    },
    "phase1": {
        "phase2": 0.52,  # Oncology: 52% (lower than all-indication 66%)
        "rejected": 0.48,
    },
    "phase2": {
        "phase3": 0.29,  # Oncology: 29% (lower than all-indication 39%)
        "rejected": 0.71,
    },
    "phase3": {
        "NDA": 0.58,     # Oncology: 58% (vs all-indication 68%)
        "rejected": 0.42,
    },
    "NDA": {
        "approved": 0.85,
        "rejected": 0.15,
    },
    "approved": {
        "approved": 1.0,  # Absorbing state
    },
    "rejected": {
        "rejected": 1.0,  # Absorbing state
    },
}

STATES = list(BASE_TRANSITIONS.keys())

# Expected time in each state (years) — for timeline estimation
STATE_DURATIONS: dict[str, float] = {
    "discovery": 2.0,
    "lead_optimization": 1.5,
    "IND_enabling": 1.0,
    "phase1": 1.5,
    "phase2": 2.5,
    "phase3": 3.0,
    "NDA": 1.5,
    "approved": 0.0,
    "rejected": 0.0,
}


@dataclass
class CandidateState:
    """Track a drug candidate's position in the development pipeline."""
    target: str
    candidate_name: str
    current_state: str = "discovery"
    transition_adjustments: dict[str, dict[str, float]] = field(default_factory=dict)
    history: list[dict] = field(default_factory=list)

    def get_transitions(self) -> dict[str, float]:
        """Get current transition probabilities (base + adjustments)."""
        base = dict(BASE_TRANSITIONS.get(self.current_state, {}))

        # Apply target-specific adjustments
        adj = self.transition_adjustments.get(self.current_state, {})
        for next_state, delta in adj.items():
            if next_state in base:
                base[next_state] = max(0.01, min(0.99, base[next_state] + delta))

        # Renormalize
        total = sum(base.values())
        if total > 0:
            base = {k: v / total for k, v in base.items()}

        return base

    def adjust_transition(self, from_state: str, to_state: str, delta: float, reason: str = "") -> None:
        """Adjust a transition probability based on evidence.

        Positive delta = more likely to advance.
        Negative delta = less likely to advance.
        """
        if from_state not in self.transition_adjustments:
            self.transition_adjustments[from_state] = {}
        current = self.transition_adjustments[from_state].get(to_state, 0.0)
        self.transition_adjustments[from_state][to_state] = current + delta
        self.history.append({
            "action": "adjust",
            "from": from_state,
            "to": to_state,
            "delta": delta,
            "reason": reason,
        })

    def p_approval(self) -> float:
        """Calculate cumulative probability of reaching approval from current state.

        Multiplies forward transition probabilities through remaining states.
        """
        prob = 1.0
        state = self.current_state
        while state not in ("approved", "rejected"):
            transitions = dict(BASE_TRANSITIONS.get(state, {}))
            adj = self.transition_adjustments.get(state, {})
            for ns, delta in adj.items():
                if ns in transitions:
                    transitions[ns] = max(0.01, min(0.99, transitions[ns] + delta))
            # Renormalize
            total = sum(transitions.values())
            transitions = {k: v / total for k, v in transitions.items()} if total > 0 else transitions

            # Find the "advance" transition (not rejected)
            advance_prob = sum(v for k, v in transitions.items() if k != "rejected")
            prob *= advance_prob

            # Find next non-rejected state
            next_states = [k for k in transitions if k != "rejected" and transitions[k] > 0]
            if not next_states:
                break
            state = next_states[0]

        return round(prob, 4)

    def p_approval_mc(self, n_sims: int = 10_000, seed: int = 42) -> dict[str, Any]:
        """Monte Carlo estimate of P(approval) with confidence intervals.

        Instead of multiplying point-estimate transition probabilities, samples
        each phase transition from a Beta distribution centered on the adjusted
        base rate. Returns median, mean, and 90% CI.

        The Beta shape parameters are derived from an effective sample size of 50
        (moderate confidence in BIO/Informa base rates). Higher N = tighter prior.
        """
        rng = np.random.default_rng(seed)
        effective_n = 50  # pseudo-observations backing each base rate

        # Build the ordered forward path from current state
        forward_states = []
        state = self.current_state
        while state not in ("approved", "rejected"):
            transitions = dict(BASE_TRANSITIONS.get(state, {}))
            adj = self.transition_adjustments.get(state, {})
            for ns, delta in adj.items():
                if ns in transitions:
                    transitions[ns] = max(0.01, min(0.99, transitions[ns] + delta))
            total = sum(transitions.values())
            transitions = {k: v / total for k, v in transitions.items()} if total > 0 else transitions

            advance_prob = sum(v for k, v in transitions.items() if k != "rejected")
            forward_states.append((state, advance_prob))

            next_states = [k for k in transitions if k != "rejected" and transitions[k] > 0]
            if not next_states:
                break
            state = next_states[0]

        if not forward_states:
            return {
                "mean": 1.0 if self.current_state == "approved" else 0.0,
                "median": 1.0 if self.current_state == "approved" else 0.0,
                "ci_5": 1.0 if self.current_state == "approved" else 0.0,
                "ci_95": 1.0 if self.current_state == "approved" else 0.0,
                "n_sims": n_sims,
                "per_phase": {},
            }

        # Sample each phase transition independently from Beta(alpha, beta)
        # where alpha = p * N, beta = (1-p) * N
        phase_samples = {}
        for state_name, p in forward_states:
            alpha = p * effective_n
            beta = (1.0 - p) * effective_n
            samples = rng.beta(alpha, beta, size=n_sims)
            phase_samples[state_name] = samples

        # Cumulative P(approval) = product of all phase advances
        cumulative = np.ones(n_sims)
        for state_name, _ in forward_states:
            cumulative *= phase_samples[state_name]

        # Per-phase summary stats
        per_phase = {}
        for state_name, p in forward_states:
            s = phase_samples[state_name]
            per_phase[state_name] = {
                "point_estimate": round(p, 4),
                "mc_mean": round(float(np.mean(s)), 4),
                "mc_ci_5": round(float(np.percentile(s, 5)), 4),
                "mc_ci_95": round(float(np.percentile(s, 95)), 4),
            }

        return {
            "mean": round(float(np.mean(cumulative)), 6),
            "median": round(float(np.median(cumulative)), 6),
            "ci_5": round(float(np.percentile(cumulative, 5)), 6),
            "ci_95": round(float(np.percentile(cumulative, 95)), 6),
            "std": round(float(np.std(cumulative)), 6),
            "n_sims": n_sims,
            "point_estimate": self.p_approval(),
            "per_phase": per_phase,
        }

    def expected_timeline_years(self) -> float:
        """Estimate years to approval from current state."""
        total = 0.0
        idx = STATES.index(self.current_state) if self.current_state in STATES else 0
        for state in STATES[idx:]:
            if state in ("approved", "rejected"):
                break
            total += STATE_DURATIONS.get(state, 1.0)
        return round(total, 1)

    def to_dict(self, include_mc: bool = True) -> dict[str, Any]:
        result = {
            "target": self.target,
            "candidate": self.candidate_name,
            "current_state": self.current_state,
            "p_approval": self.p_approval(),
            "expected_years": self.expected_timeline_years(),
            "transitions": self.get_transitions(),
            "adjustments": self.transition_adjustments,
        }
        if include_mc:
            result["p_approval_mc"] = self.p_approval_mc()
        return result


def create_candidates_from_targets(targets: dict) -> list[CandidateState]:
    """Initialize CandidateState objects from targets.yaml.

    Infers current pipeline state from drug_status and known_compounds.
    """
    candidates = []

    for gene, cfg in targets.items():
        drug_status = cfg.get("drug_status", "").lower()
        compounds = cfg.get("known_compounds", [])

        if not compounds:
            # No compounds — still in discovery
            candidates.append(CandidateState(
                target=gene,
                candidate_name=f"{gene}_novel",
                current_state="discovery",
            ))
        else:
            for cpd in compounds:
                name = cpd.get("name", "unknown")
                ref = cpd.get("reference", "").lower()
                repurposing = cpd.get("repurposing", "").lower()

                # Infer state from references
                if "phase iii" in ref or "phase 3" in ref or "phase iii" in repurposing:
                    state = "phase3"
                elif "phase ii" in ref or "phase 2" in ref or "phase ii" in repurposing:
                    state = "phase2"
                elif "phase i" in ref or "phase 1" in ref:
                    state = "phase1"
                elif "fda" in ref and "approved" in ref:
                    state = "approved"
                elif "fda" in ref:
                    state = "NDA"
                elif "clinical" in ref or "trial" in ref:
                    state = "phase1"
                elif "nanomolar" in ref or "inhibitor" in ref:
                    state = "lead_optimization"
                else:
                    state = "discovery"

                # Repurposing candidates skip some early stages
                if "repurpos" in drug_status or "repurpos" in repurposing:
                    if state in ("discovery", "lead_optimization"):
                        state = "IND_enabling"

                candidate = CandidateState(
                    target=gene,
                    candidate_name=name,
                    current_state=state,
                )

                # Adjust probabilities based on known attributes
                selectivity = cpd.get("selectivity", "").lower()
                if "selective" in selectivity and "non" not in selectivity:
                    candidate.adjust_transition(
                        "phase1", "phase2", 0.05,
                        reason=f"Selective compound ({selectivity})"
                    )
                if "repurpos" in repurposing:
                    # Repurposed drugs have better safety profile
                    candidate.adjust_transition(
                        "phase1", "phase2", 0.10,
                        reason="Existing safety data from approved indication"
                    )

                candidates.append(candidate)

    return candidates
