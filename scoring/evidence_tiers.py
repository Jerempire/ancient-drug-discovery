"""Crawford-Sobel evidence tier classification for drug discovery.

Adapted from prediction-market-models/models/geopolitical_bayes.py.

Evidence is classified into tiers based on costly-signal theory (Crawford & Sobel 1982):
signals that are expensive to fabricate carry more information than cheap talk.

In drug discovery, the cost hierarchy is:
  Tier 1: Irreproducible wet-lab results (Phase 3, CRISPR KO) — hardest to fake
  Tier 2a: Regulated/reviewed data (FDA, published RCT) — institutional reputation at stake
  Tier 2b: Cheap-to-produce announcements (press releases, abstracts) — no replication yet
  Tier 3a: Expert opinion with conflicts (KOL with competing program) — biased but informed
  Tier 3b: Unvetted sources (preprints, anonymous analysts) — no accountability
  Tier 4: Pure speculation (social media, rumors) — zero information content
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class EvidenceTier:
    tier: str
    label: str
    weight: float
    description: str


TIERS: dict[str, EvidenceTier] = {
    "1": EvidenceTier(
        tier="1",
        label="Hard causal",
        weight=1.0,
        description="Phase 3 readout, confirmed CRISPR knockout, published RCT with replication",
    ),
    "2a": EvidenceTier(
        tier="2a",
        label="Costly soft",
        weight=0.95,
        description="FDA designation, published peer-reviewed RCT, IND approval",
    ),
    "2b": EvidenceTier(
        tier="2b",
        label="Cheap talk",
        weight=0.70,
        description="Company press release, conference abstract, poster presentation",
    ),
    "3a": EvidenceTier(
        tier="3a",
        label="Costly commentary",
        weight=0.65,
        description="KOL opinion with competing program, editorial in high-impact journal",
    ),
    "3b": EvidenceTier(
        tier="3b",
        label="Unvetted",
        weight=0.30,
        description="Anonymous analyst, preprint without peer review, blog post",
    ),
    "4": EvidenceTier(
        tier="4",
        label="Speculation",
        weight=0.0,
        description="Social media rumor, AI-generated hypothesis without validation",
    ),
}


# Keywords that trigger each tier classification
_TIER_KEYWORDS: dict[str, list[str]] = {
    "1": [
        "phase 3", "phase iii", "pivotal trial", "crispr knockout",
        "replicated rct", "confirmed in vivo", "clinical endpoint met",
        "survival benefit", "overall survival",
    ],
    "2a": [
        "phase 2", "phase ii", "fda", "ema", "breakthrough therapy",
        "orphan drug", "fast track", "published rct", "peer-reviewed",
        "nature", "science", "cell", "nejm", "lancet", "jama",
        "pnas", "immunity", "j med chem",
    ],
    "2b": [
        "press release", "conference abstract", "poster", "aacr",
        "asco", "esmo", "preclinical", "in vitro", "cell line",
        "phase 1", "phase i",
    ],
    "3a": [
        "kol", "editorial", "review article", "expert opinion",
        "advisory board", "competing program",
    ],
    "3b": [
        "preprint", "biorxiv", "medrxiv", "anonymous", "analyst note",
        "blog", "unreviewed",
    ],
    "4": [
        "rumor", "speculation", "social media", "twitter", "reddit",
        "ai-generated", "unvalidated hypothesis",
    ],
}


def classify_evidence(
    source_type: str = "",
    source_ref: str = "",
    text: str = "",
    year: int | None = None,
) -> EvidenceTier:
    """Classify a piece of evidence into a Crawford-Sobel tier.

    Uses keyword matching against source_type, source_ref, and free text.
    Falls back to tier 2b (cheap talk) when ambiguous — conservative default.
    """
    combined = f"{source_type} {source_ref} {text}".lower()

    # Check tiers from highest to lowest — first match wins
    for tier_key in ["1", "2a", "2b", "3a", "3b", "4"]:
        keywords = _TIER_KEYWORDS[tier_key]
        if any(kw in combined for kw in keywords):
            return TIERS[tier_key]

    # Default: tier 2b (cheap talk) — most evidence falls here
    return TIERS["2b"]


def classify_source_type(source_type: str) -> EvidenceTier:
    """Classify by MedGraph-Rx source_type enum directly."""
    mapping = {
        "rct": TIERS["2a"],
        "meta_analysis": TIERS["1"],
        "fda_label": TIERS["2a"],
        "guideline": TIERS["2a"],
        "cohort": TIERS["2b"],
        "seed_curation": TIERS["3a"],
        "faers": TIERS["2b"],
        "expert_opinion": TIERS["3a"],
    }
    return mapping.get(source_type, TIERS["2b"])


def tier_weight(tier: EvidenceTier, year: int | None = None, current_year: int = 2026) -> float:
    """Get effective weight with optional time decay.

    Evidence older than 10 years gets a 20% discount.
    Evidence older than 20 years gets a 40% discount.
    """
    w = tier.weight
    if year and current_year:
        age = current_year - year
        if age > 20:
            w *= 0.6
        elif age > 10:
            w *= 0.8
    return w
