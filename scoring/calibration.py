"""Calibration tracking — record predictions vs outcomes over time.

Tracks Bayesian estimates against ground truth as it arrives:
- Experiment outcomes (IC50 assay confirms tumor_relevance?)
- Literature validation (new paper confirms/contradicts druggability?)
- Trial readouts (Phase 2 result matches clinical_validation prediction?)

Storage: JSON file at data/calibration_log.json.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .bayes_target import TargetBetaEstimator

DEFAULT_CALIBRATION_PATH = Path(__file__).parent.parent / "data" / "calibration_log.json"


@dataclass
class Prediction:
    target: str
    dimension: str
    predicted_mean: float
    predicted_ci_low: float
    predicted_ci_high: float
    timestamp: str
    outcome: float | None = None
    outcome_timestamp: str | None = None


class CalibrationTracker:
    """Track predictions and outcomes for calibration analysis."""

    def __init__(self, path: Path | None = None):
        self.path = path or DEFAULT_CALIBRATION_PATH
        self._predictions: list[Prediction] = []
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            return
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return
        for entry in data:
            self._predictions.append(Prediction(**entry))

    def _save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps([asdict(p) for p in self._predictions], indent=2),
            encoding="utf-8",
        )

    def record_prediction(
        self,
        target: str,
        dimension: str,
        estimator: TargetBetaEstimator,
    ) -> None:
        """Record a prediction from a BetaEstimator's current state."""
        est = estimator.get_estimate()
        pred = Prediction(
            target=target,
            dimension=dimension,
            predicted_mean=est["mean"],
            predicted_ci_low=est["ci_low"],
            predicted_ci_high=est["ci_high"],
            timestamp=datetime.now(timezone.utc).isoformat(),
        )
        self._predictions.append(pred)
        self._save()

    def record_outcome(
        self,
        target: str,
        dimension: str,
        outcome: float,
    ) -> None:
        """Attach an outcome to the most recent unresolved prediction for target+dimension."""
        for pred in reversed(self._predictions):
            if pred.target == target and pred.dimension == dimension and pred.outcome is None:
                pred.outcome = outcome
                pred.outcome_timestamp = datetime.now(timezone.utc).isoformat()
                self._save()
                return

    def calibration_report(self) -> dict[str, Any]:
        """Generate calibration statistics."""
        resolved = [p for p in self._predictions if p.outcome is not None]
        n_total = len(self._predictions)
        n_resolved = len(resolved)

        if not resolved:
            return {
                "n_predictions": n_total,
                "n_resolved": 0,
                "coverage_90": None,
                "mean_absolute_error": None,
                "brier_score": None,
                "by_dimension": {},
                "by_source": {},
            }

        # Coverage: fraction of outcomes within 90% CI
        # The stored CI is 95% (±1.96σ). Scale to 90% (±1.645σ) by shrinking.
        in_ci = 0
        abs_errors = []
        brier_terms = []
        by_dim: dict[str, list[float]] = {}
        for p in resolved:
            err = abs(p.predicted_mean - p.outcome)
            abs_errors.append(err)

            # 90% CI: shrink the 95% CI by ratio 1.645/1.96
            ci_half_95 = (p.predicted_ci_high - p.predicted_ci_low) / 2
            ci_half_90 = ci_half_95 * (1.645 / 1.96)
            ci_center = (p.predicted_ci_high + p.predicted_ci_low) / 2
            if ci_center - ci_half_90 <= p.outcome <= ci_center + ci_half_90:
                in_ci += 1

            # Brier score for binary interpretation (outcome >= 0.5 = "confirmed")
            binary_outcome = 1.0 if p.outcome >= 0.5 else 0.0
            brier_terms.append((p.predicted_mean - binary_outcome) ** 2)

            by_dim.setdefault(p.dimension, []).append(err)

        coverage_90 = in_ci / n_resolved
        mae = sum(abs_errors) / n_resolved
        brier = sum(brier_terms) / n_resolved

        dim_report = {
            dim: {"n": len(errs), "mae": round(sum(errs) / len(errs), 4)}
            for dim, errs in by_dim.items()
        }

        return {
            "n_predictions": n_total,
            "n_resolved": n_resolved,
            "coverage_90": round(coverage_90, 4),
            "mean_absolute_error": round(mae, 4),
            "brier_score": round(brier, 4),
            "by_dimension": dim_report,
            "by_source": {},  # Populated when source tracking is added to Prediction
        }
