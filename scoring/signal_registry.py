"""Signal registry for drug target scoring dimensions.

Ported from alpha-engine/engine/signal_loader.py.

Each scoring dimension is a 'signal' with:
- Status tracking (active, candidate, decaying, retired)
- Source list (which evidence types feed into it)
- Decay configuration
- Baseline prior

The registry is the single source of truth for what dimensions exist
and how they're configured. The TargetScorer in bayes_target.py reads
dimension weights from here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

REGISTRY_PATH = Path(__file__).parent / "registry.yaml"

VALID_STATUSES = {"active", "candidate", "decaying", "retired"}
VALID_DOMAINS = {"oncology", "genomics", "chemistry", "clinical", "commercial", "safety"}


@dataclass
class SignalSpec:
    name: str
    description: str
    domain: str
    status: str
    weight: float
    sources: list[str] = field(default_factory=list)
    decay: dict[str, Any] = field(default_factory=dict)
    baseline: float = 0.5
    notes: str = ""

    def __post_init__(self):
        if self.status not in VALID_STATUSES:
            raise ValueError(f"Invalid status '{self.status}'. Must be one of {VALID_STATUSES}")
        if self.domain not in VALID_DOMAINS:
            raise ValueError(f"Invalid domain '{self.domain}'. Must be one of {VALID_DOMAINS}")


class SignalRegistry:
    """Load and query the signal registry."""

    def __init__(self, path: Path | None = None):
        self.path = path or REGISTRY_PATH
        self._signals: dict[str, SignalSpec] = {}
        self._correlations: dict[tuple[str, str], float] = {}
        self._load()

    def _load(self) -> None:
        data = yaml.safe_load(self.path.read_text(encoding="utf-8"))
        for key, cfg in data.get("signals", {}).items():
            self._signals[key] = SignalSpec(
                name=cfg["name"],
                description=cfg["description"],
                domain=cfg["domain"],
                status=cfg["status"],
                weight=cfg["weight"],
                sources=cfg.get("sources", []),
                decay=cfg.get("decay", {}),
                baseline=cfg.get("baseline", 0.5),
                notes=cfg.get("notes", ""),
            )
        for entry in data.get("correlations", []):
            dim_a, dim_b, rho = entry[0], entry[1], float(entry[2])
            self._correlations[(dim_a, dim_b)] = rho
            self._correlations[(dim_b, dim_a)] = rho

    def active(self) -> dict[str, SignalSpec]:
        return {k: v for k, v in self._signals.items() if v.status == "active"}

    def by_domain(self, domain: str) -> dict[str, SignalSpec]:
        return {k: v for k, v in self._signals.items() if v.domain == domain}

    def by_status(self, status: str) -> dict[str, SignalSpec]:
        return {k: v for k, v in self._signals.items() if v.status == status}

    def get(self, key: str) -> SignalSpec | None:
        return self._signals.get(key)

    def active_weights(self) -> dict[str, float]:
        """Return {dimension: weight} for all active signals. Normalized to sum=1."""
        active = self.active()
        raw = {k: v.weight for k, v in active.items()}
        total = sum(raw.values())
        if total == 0:
            return raw
        return {k: round(w / total, 4) for k, w in raw.items()}

    def correlation(self, dim_a: str, dim_b: str) -> float:
        """Return pairwise correlation between two dimensions (0.0 if not defined)."""
        return self._correlations.get((dim_a, dim_b), 0.0)

    def correlation_matrix(self) -> dict[tuple[str, str], float]:
        """Return all defined pairwise correlations."""
        return dict(self._correlations)

    @property
    def all_signals(self) -> dict[str, SignalSpec]:
        return dict(self._signals)
