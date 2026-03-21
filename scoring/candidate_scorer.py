"""Candidate-level binder scorer for ranking individual designs.

Scores each binder across 5 categories (interface quality, stability,
selectivity, developability, mechanistic value), applies hard filters,
and tiers candidates as TIER_1 / TIER_2 / TIER_3 / KILL.
"""

import json
import math
import csv
import sqlite3
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from collections import Counter
from typing import Optional

try:
    from .developability import assess_developability, DevelopabilityProfile
except ImportError:
    from developability import assess_developability, DevelopabilityProfile

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Divergent hotspot residues (ERAP2 numbering, chain A)
HOTSPOT_RESIDUES = {353, 355, 360, 363, 367, 401, 403, 406, 408, 412}

# Category weights
WEIGHTS = {
    "interface": 0.30,
    "stability": 0.20,
    "selectivity": 0.30,
    "developability": 0.10,
    "mechanism": 0.10,
}


@dataclass
class CandidateScore:
    candidate_id: str
    target: str
    site: str
    family: str
    binder_length: int
    sequence: str
    interface_score: float
    stability_score: float
    selectivity_score: float
    developability_score: float
    mechanism_score: float
    total_score: float
    tier: str
    hard_filter_pass: bool
    red_flags: list = field(default_factory=list)
    raw_metrics: dict = field(default_factory=dict)


def load_validation_data(path: Path) -> list[dict]:
    """Read validation_summary.json."""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def load_selectivity_data(path: Path) -> dict[str, dict]:
    """Read selectivity_v3.json, keyed by design name (without .pdb)."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    result = {}
    for entry in data:
        key = entry["design"].replace(".pdb", "")
        result[key] = entry
    return result


def compute_seq_complexity(seq: str) -> float:
    """Shannon entropy of amino acid composition (bits)."""
    counts = Counter(seq)
    n = len(seq)
    if n == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def _normalize(value: float, min_val: float, max_val: float) -> float:
    """Min-max normalize to 0-5 scale. Clamp to [0, 5]."""
    if max_val == min_val:
        return 2.5
    score = 5.0 * (value - min_val) / (max_val - min_val)
    return max(0.0, min(5.0, score))


def _normalize_inverted(value: float, min_val: float, max_val: float) -> float:
    """Inverted min-max normalize (lower raw = higher score)."""
    if max_val == min_val:
        return 2.5
    score = 5.0 * (max_val - value) / (max_val - min_val)
    return max(0.0, min(5.0, score))


def _length_penalty(length: int) -> float:
    """Score 0-5 for binder length. Optimal 30-120, penalize extremes."""
    if 30 <= length <= 120:
        return 5.0
    if length < 30:
        return max(0.0, 5.0 * length / 30)
    # > 120
    return max(0.0, 5.0 * (1.0 - (length - 120) / 120))


def _parse_site_from_design(design: str) -> str:
    """Extract site tag from design name."""
    if "channel" in design:
        return "divergent_channel"
    if "cterm" in design:
        return "conserved_cterm"
    if "domainIV" in design:
        return "domainIV"
    return "unknown"


def _parse_family_from_design(design: str) -> str:
    """Extract family (short/medium/long) if present in name."""
    for fam in ("short", "medium", "long"):
        if fam in design:
            return fam
    return "unknown"


def _collect_raw_metrics(candidates: list[dict], selectivity_data: dict) -> list[dict]:
    """Build raw metrics dict for each candidate."""
    results = []
    for c in candidates:
        design = c["design"]
        seq = c["binder_sequence"]
        sel = selectivity_data.get(design, {})

        # Binder pLDDT from pair_chains_iptm diagonal (chain B = index 1)
        binder_plddt = c["erap2_scores"]["pair_chains_iptm"]["1"]["1"]

        # Full developability assessment
        dev = assess_developability(seq)

        metrics = {
            "design": design,
            "sequence": seq,
            "binder_length": c["binder_length"],
            # Interface
            "iptm_erap2": c["erap2_scores"]["iptm"],
            "complex_plddt": c["erap2_scores"]["complex_plddt"],
            "ptm_erap2": c["erap2_scores"]["ptm"],
            # Stability
            "binder_plddt": binder_plddt,
            "mpnn_score": c["mpnn_score"],
            # Selectivity (Boltz-2)
            "iptm_delta": c["iptm_delta"],
            "iptm_selectivity": c["iptm_selectivity"],
            # Selectivity (contact-based)
            "divergent_frac": sel.get("divergent_frac"),
            "interface_selectivity": sel.get("interface_selectivity"),
            # Developability (from developability module)
            "seq_complexity": dev.shannon_entropy,
            "charge": dev.net_charge_ph7,
            "isoelectric_point": dev.isoelectric_point,
            "hydrophobic_fraction": dev.hydrophobic_fraction,
            "max_hydrophobic_patch": dev.max_hydrophobic_patch,
            "mean_hydrophobicity": dev.mean_hydrophobicity,
            "ptm_liability_count": dev.ptm_liability_count,
            "deamidation_sites": len(dev.deamidation_sites),
            "oxidation_sites": len(dev.oxidation_sites),
            "isomerization_sites": len(dev.isomerization_sites),
            "dev_flags": dev.flags,
            # Mechanism
            "divergent_residues_hit": sel.get("divergent_residues_hit", []),
            "site": _parse_site_from_design(design),
            "family": _parse_family_from_design(design),
        }
        results.append(metrics)
    return results


def _check_hard_filters(m: dict) -> tuple[bool, list[str]]:
    """Apply hard filters. Returns (pass, list_of_reasons)."""
    flags = []
    if m["iptm_erap2"] < 0.05:
        flags.append(f"iptm_erap2={m['iptm_erap2']:.4f} < 0.05 (no meaningful binding)")
    if m["iptm_delta"] < -0.10:
        flags.append(f"iptm_delta={m['iptm_delta']:.4f} < -0.10 (prefers ERAP1)")
    if m["seq_complexity"] < 1.5:
        flags.append(f"seq_complexity={m['seq_complexity']:.2f} < 1.5 (likely poly-alanine)")
    if m["binder_length"] < 20:
        flags.append(f"binder_length={m['binder_length']} < 20 (too short)")
    return (len(flags) == 0, flags)


def score_candidates(candidates: list[dict], selectivity_data: dict) -> list[CandidateScore]:
    """Score all candidates and return sorted list."""
    all_metrics = _collect_raw_metrics(candidates, selectivity_data)

    if not all_metrics:
        return []

    # Compute dataset-wide min/max for normalization
    def _minmax(key):
        vals = [m[key] for m in all_metrics if m[key] is not None]
        if not vals:
            return 0.0, 1.0
        return min(vals), max(vals)

    iptm_min, iptm_max = _minmax("iptm_erap2")
    cplddt_min, cplddt_max = _minmax("complex_plddt")
    ptm_min, ptm_max = _minmax("ptm_erap2")
    bplddt_min, bplddt_max = _minmax("binder_plddt")
    mpnn_min, mpnn_max = _minmax("mpnn_score")
    delta_min, delta_max = _minmax("iptm_delta")
    sel_min, sel_max = _minmax("iptm_selectivity")
    dfrac_min, dfrac_max = _minmax("divergent_frac")
    isel_min, isel_max = _minmax("interface_selectivity")

    scores = []
    for m in all_metrics:
        passed, kill_reasons = _check_hard_filters(m)

        # 1. Interface Quality (0.30)
        s_iptm = _normalize(m["iptm_erap2"], iptm_min, iptm_max)
        s_cplddt = _normalize(m["complex_plddt"], cplddt_min, cplddt_max)
        s_ptm = _normalize(m["ptm_erap2"], ptm_min, ptm_max)
        interface = (s_iptm * 0.5 + s_cplddt * 0.25 + s_ptm * 0.25)

        # 2. Stability (0.20)
        s_bplddt = _normalize(m["binder_plddt"], bplddt_min, bplddt_max)
        s_mpnn = _normalize_inverted(m["mpnn_score"], mpnn_min, mpnn_max)  # lower = better
        s_length = _length_penalty(m["binder_length"])
        stability = (s_bplddt * 0.4 + s_mpnn * 0.4 + s_length * 0.2)

        # 3. Selectivity (0.30)
        s_delta = _normalize(m["iptm_delta"], delta_min, delta_max)
        s_sel = _normalize(m["iptm_selectivity"], sel_min, sel_max)
        if m["divergent_frac"] is not None:
            s_dfrac = _normalize(m["divergent_frac"], dfrac_min, dfrac_max)
            s_isel = _normalize(m["interface_selectivity"], isel_min, isel_max)
            selectivity = (s_delta * 0.3 + s_sel * 0.3 + s_dfrac * 0.2 + s_isel * 0.2)
        else:
            selectivity = (s_delta * 0.5 + s_sel * 0.5)

        # 4. Developability (0.10) — uses full developability profile
        complexity = m["seq_complexity"]
        if complexity >= 3.5:
            s_complexity = 5.0
        elif complexity >= 2.5:
            s_complexity = 5.0 * (complexity - 1.5) / 2.0
        else:
            s_complexity = max(0.0, 5.0 * (complexity - 1.0) / 1.5)

        charge = abs(m["charge"])
        if charge <= 5:
            s_charge = 5.0
        elif charge <= 10:
            s_charge = 5.0 * (1.0 - (charge - 5) / 10)
        else:
            s_charge = max(0.0, 5.0 * (1.0 - (charge - 5) / 20))

        # Aggregation risk: penalize large hydrophobic patches
        patch = m["max_hydrophobic_patch"]
        if patch <= 4:
            s_aggregation = 5.0
        elif patch <= 7:
            s_aggregation = 5.0 * (1.0 - (patch - 4) / 6)
        else:
            s_aggregation = max(0.0, 5.0 * (1.0 - (patch - 4) / 10))

        # PTM liability: penalize high counts
        ptm_count = m["ptm_liability_count"]
        if ptm_count <= 2:
            s_ptm = 5.0
        elif ptm_count <= 5:
            s_ptm = 5.0 * (1.0 - (ptm_count - 2) / 6)
        else:
            s_ptm = max(0.0, 5.0 * (1.0 - (ptm_count - 2) / 10))

        developability = (
            s_complexity * 0.25
            + s_charge * 0.25
            + s_aggregation * 0.25
            + s_ptm * 0.25
        )

        # 5. Mechanistic Value (0.10)
        hotspots_hit = set(m["divergent_residues_hit"]) & HOTSPOT_RESIDUES
        hotspot_coverage = len(hotspots_hit) / len(HOTSPOT_RESIDUES) * 5.0
        site_bonus = 5.0 if m["site"] == "divergent_channel" else 2.5
        mechanism = (hotspot_coverage * 0.6 + site_bonus * 0.4)

        # Total
        total = (
            interface * WEIGHTS["interface"]
            + stability * WEIGHTS["stability"]
            + selectivity * WEIGHTS["selectivity"]
            + developability * WEIGHTS["developability"]
            + mechanism * WEIGHTS["mechanism"]
        )

        # Red flags (non-fatal warnings)
        red_flags = []
        if m["seq_complexity"] < 2.5:
            red_flags.append(f"low_complexity ({m['seq_complexity']:.2f})")
        if abs(m["charge"]) > 8:
            red_flags.append(f"extreme_charge ({m['charge']:+.1f})")
        if m["iptm_delta"] < 0:
            red_flags.append(f"negative_selectivity (delta={m['iptm_delta']:.4f})")
        red_flags.extend(m.get("dev_flags", []))

        # Tier
        if not passed:
            tier = "KILL"
        elif total >= 3.5:
            tier = "TIER_1"
        elif total >= 2.5:
            tier = "TIER_2"
        else:
            tier = "TIER_3"

        scores.append(CandidateScore(
            candidate_id=m["design"],
            target="ERAP2",
            site=m["site"],
            family=m["family"],
            binder_length=m["binder_length"],
            sequence=m["sequence"],
            interface_score=round(interface, 3),
            stability_score=round(stability, 3),
            selectivity_score=round(selectivity, 3),
            developability_score=round(developability, 3),
            mechanism_score=round(mechanism, 3),
            total_score=round(total, 3),
            tier=tier,
            hard_filter_pass=passed,
            red_flags=red_flags if passed else kill_reasons,
            raw_metrics={k: v for k, v in m.items()
                         if k not in ("sequence", "divergent_residues_hit", "dev_flags")},
        ))

    # Sort: KILL last, then by total descending
    tier_order = {"TIER_1": 0, "TIER_2": 1, "TIER_3": 2, "KILL": 3}
    scores.sort(key=lambda s: (tier_order[s.tier], -s.total_score))
    return scores


def export_results(scores: list[CandidateScore], out_dir: Path) -> None:
    """Write JSON and CSV outputs."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = out_dir / "candidate_ranking.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump([asdict(s) for s in scores], f, indent=2)
    print(f"Wrote {json_path}")

    # CSV
    csv_path = out_dir / "candidate_ranking.csv"
    csv_fields = [
        "candidate_id", "tier", "total_score",
        "interface_score", "stability_score", "selectivity_score",
        "developability_score", "mechanism_score",
        "binder_length", "site", "family", "hard_filter_pass", "red_flags",
    ]
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for s in scores:
            row = asdict(s)
            row["red_flags"] = "; ".join(row["red_flags"])
            writer.writerow({k: row[k] for k in csv_fields})
    print(f"Wrote {csv_path}")


DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS design_rounds (
    round_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    source_file TEXT,
    scored_at TEXT NOT NULL,
    num_candidates INTEGER,
    num_tier1 INTEGER,
    num_tier2 INTEGER,
    num_kill INTEGER,
    notes TEXT
);

CREATE TABLE IF NOT EXISTS candidates (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    round_id INTEGER NOT NULL REFERENCES design_rounds(round_id),
    candidate_id TEXT NOT NULL,
    target TEXT NOT NULL DEFAULT 'ERAP2',
    site TEXT,
    family TEXT,
    sequence TEXT NOT NULL,
    binder_length INTEGER,
    -- Category scores (0-5)
    interface_score REAL,
    stability_score REAL,
    selectivity_score REAL,
    developability_score REAL,
    mechanism_score REAL,
    total_score REAL,
    tier TEXT,
    hard_filter_pass INTEGER,
    red_flags TEXT,
    -- Raw Boltz-2 metrics
    iptm_erap2 REAL,
    iptm_erap1 REAL,
    iptm_delta REAL,
    iptm_selectivity REAL,
    complex_plddt REAL,
    ptm_erap2 REAL,
    binder_plddt REAL,
    mpnn_score REAL,
    -- Contact-based selectivity
    divergent_frac REAL,
    interface_selectivity REAL,
    -- Developability
    seq_complexity REAL,
    charge REAL,
    isoelectric_point REAL,
    hydrophobic_fraction REAL,
    max_hydrophobic_patch INTEGER,
    mean_hydrophobicity REAL,
    ptm_liability_count INTEGER,
    deamidation_sites INTEGER,
    oxidation_sites INTEGER,
    isomerization_sites INTEGER,
    UNIQUE(round_id, candidate_id)
);

CREATE INDEX IF NOT EXISTS idx_candidates_tier ON candidates(tier);
CREATE INDEX IF NOT EXISTS idx_candidates_total ON candidates(total_score DESC);
CREATE INDEX IF NOT EXISTS idx_candidates_round ON candidates(round_id);
"""


def init_db(db_path: Path) -> sqlite3.Connection:
    """Initialize SQLite database with schema."""
    conn = sqlite3.connect(str(db_path))
    conn.executescript(DB_SCHEMA)
    return conn


def export_to_db(scores: list[CandidateScore], db_path: Path,
                 round_name: str, source_file: str = "") -> None:
    """Write scored candidates to SQLite database."""
    conn = init_db(db_path)
    tier_counts = Counter(s.tier for s in scores)

    # Insert design round
    conn.execute(
        """INSERT INTO design_rounds (name, source_file, scored_at,
           num_candidates, num_tier1, num_tier2, num_kill, notes)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        (round_name, source_file, datetime.now().isoformat(),
         len(scores), tier_counts.get("TIER_1", 0),
         tier_counts.get("TIER_2", 0), tier_counts.get("KILL", 0), None)
    )
    round_id = conn.execute("SELECT last_insert_rowid()").fetchone()[0]

    # Insert candidates
    for s in scores:
        m = s.raw_metrics
        conn.execute(
            """INSERT OR REPLACE INTO candidates (
                round_id, candidate_id, target, site, family, sequence, binder_length,
                interface_score, stability_score, selectivity_score,
                developability_score, mechanism_score, total_score,
                tier, hard_filter_pass, red_flags,
                iptm_erap2, iptm_delta, iptm_selectivity, complex_plddt,
                ptm_erap2, binder_plddt, mpnn_score,
                divergent_frac, interface_selectivity,
                seq_complexity, charge, isoelectric_point,
                hydrophobic_fraction, max_hydrophobic_patch, mean_hydrophobicity,
                ptm_liability_count, deamidation_sites, oxidation_sites, isomerization_sites
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (round_id, s.candidate_id, s.target, s.site, s.family,
             s.sequence, s.binder_length,
             s.interface_score, s.stability_score, s.selectivity_score,
             s.developability_score, s.mechanism_score, s.total_score,
             s.tier, int(s.hard_filter_pass), "; ".join(s.red_flags),
             m.get("iptm_erap2"), m.get("iptm_delta"), m.get("iptm_selectivity"),
             m.get("complex_plddt"), m.get("ptm_erap2"), m.get("binder_plddt"),
             m.get("mpnn_score"),
             m.get("divergent_frac"), m.get("interface_selectivity"),
             m.get("seq_complexity"), m.get("charge"), m.get("isoelectric_point"),
             m.get("hydrophobic_fraction"), m.get("max_hydrophobic_patch"),
             m.get("mean_hydrophobicity"),
             m.get("ptm_liability_count"), m.get("deamidation_sites"),
             m.get("oxidation_sites"), m.get("isomerization_sites"))
        )

    conn.commit()
    conn.close()
    print(f"Wrote {len(scores)} candidates to {db_path} (round: {round_name})")


def print_summary(scores: list[CandidateScore]) -> None:
    """Print a summary table to stdout."""
    print("\n" + "=" * 90)
    print("CANDIDATE RANKING")
    print("=" * 90)
    print(f"{'Design':<28} {'Tier':<8} {'Total':>6} {'Intf':>6} {'Stab':>6} "
          f"{'Sel':>6} {'Dev':>6} {'Mech':>6} {'Flags'}")
    print("-" * 90)
    for s in scores:
        flags = "; ".join(s.red_flags[:2]) if s.red_flags else ""
        print(f"{s.candidate_id:<28} {s.tier:<8} {s.total_score:>6.3f} "
              f"{s.interface_score:>6.3f} {s.stability_score:>6.3f} "
              f"{s.selectivity_score:>6.3f} {s.developability_score:>6.3f} "
              f"{s.mechanism_score:>6.3f} {flags}")
    print("-" * 90)

    # Tier summary
    tier_counts = Counter(s.tier for s in scores)
    print(f"\nTier Summary: ", end="")
    for tier in ["TIER_1", "TIER_2", "TIER_3", "KILL"]:
        if tier_counts.get(tier, 0) > 0:
            print(f"{tier}={tier_counts[tier]}  ", end="")
    print()


def print_developability_report(scores: list[CandidateScore]) -> None:
    """Print developability details for non-KILL candidates."""
    viable = [s for s in scores if s.tier != "KILL"]
    if not viable:
        return
    print("\n" + "=" * 90)
    print("DEVELOPABILITY REPORT")
    print("=" * 90)
    print(f"{'Design':<28} {'pI':>5} {'Charge':>7} {'HydroFrac':>9} "
          f"{'MaxPatch':>8} {'PTMLiab':>7} {'Flags'}")
    print("-" * 90)
    for s in viable:
        m = s.raw_metrics
        flags = "; ".join(s.red_flags[:3]) if s.red_flags else ""
        print(f"{s.candidate_id:<28} {m.get('isoelectric_point', 0):>5.1f} "
              f"{m.get('charge', 0):>+7.1f} {m.get('hydrophobic_fraction', 0):>9.3f} "
              f"{m.get('max_hydrophobic_patch', 0):>8d} {m.get('ptm_liability_count', 0):>7d} "
              f"{flags}")
    print()


def main():
    # Accept optional custom validation JSON path as first arg
    validation_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        PROJECT_ROOT / "data" / "results" / "boltz2_validation" / "validation_summary.json"
    )
    selectivity_path = PROJECT_ROOT / "data" / "results" / "selectivity" / "selectivity_v3.json"
    out_dir = PROJECT_ROOT / "data" / "results"

    if not validation_path.exists():
        print(f"ERROR: {validation_path} not found")
        return

    candidates = load_validation_data(validation_path)
    print(f"Loaded {len(candidates)} candidates from {validation_path.name}")

    selectivity_data = {}
    if selectivity_path.exists():
        selectivity_data = load_selectivity_data(selectivity_path)
        print(f"Loaded {len(selectivity_data)} entries from selectivity_v3.json")
    else:
        print("WARNING: selectivity_v3.json not found, scoring without contact data")

    scores = score_candidates(candidates, selectivity_data)
    print_summary(scores)
    print_developability_report(scores)
    export_results(scores, out_dir)

    # SQLite export
    db_path = out_dir / "candidates.db"
    round_name = validation_path.stem  # e.g. "validation_summary_v2"
    export_to_db(scores, db_path, round_name, str(validation_path))


if __name__ == "__main__":
    main()
