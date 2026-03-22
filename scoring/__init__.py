"""Scoring infrastructure for ancient drug discovery target prioritization."""

from .bayes_target import TargetScorer, TargetBetaEstimator, score_targets_from_yaml
from .signal_registry import SignalRegistry
from .commercial import CommercialScore, score_commercial
from .markov_pipeline import CandidateState, create_candidates_from_targets
from .ev_decision import EVEstimate, Portfolio, compute_ev, rank_targets_by_ev, select_portfolio
from .calibration import CalibrationTracker, Prediction
from .experiment_optimizer import recommend_experiments
from .candidate_scorer import CandidateScore, score_candidates, export_results
from .developability import DevelopabilityProfile, assess_developability

# PyRosetta is optional — only import if installed
try:
    from .rosetta_interface_analysis import analyze_complex, run_interface_analyzer
    _HAS_PYROSETTA = True
except ImportError:
    _HAS_PYROSETTA = False

__all__ = [
    "TargetScorer", "TargetBetaEstimator", "score_targets_from_yaml",
    "SignalRegistry",
    "CommercialScore", "score_commercial",
    "CandidateState", "create_candidates_from_targets",
    "EVEstimate", "Portfolio", "compute_ev", "rank_targets_by_ev", "select_portfolio",
    "CalibrationTracker", "Prediction",
    "recommend_experiments",
    "CandidateScore", "score_candidates", "export_results",
    "DevelopabilityProfile", "assess_developability",
]

if _HAS_PYROSETTA:
    __all__.extend(["analyze_complex", "run_interface_analyzer"])
