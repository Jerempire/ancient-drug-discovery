"""Sequence-based developability assessment for protein binder candidates.

Computes biologic-relevant ADME-equivalent metrics from sequence alone:
  - Aggregation propensity (hydrophobic patch scanning)
  - PTM liabilities (deamidation, oxidation, isomerization)
  - Isoelectric point + solubility proxy
  - Hydrophobicity profile (Kyte-Doolittle)
  - Sequence complexity (Shannon entropy)
  - Charge balance

These metrics feed into the candidate_scorer's developability category.
"""

import math
from collections import Counter
from dataclasses import dataclass, field


# Kyte-Doolittle hydrophobicity scale
KD_HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "D": -3.5,
    "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}

# pKa values for ionizable residues
PKA = {
    "D": 3.65, "E": 4.25,  # acidic
    "H": 6.00,              # His
    "C": 8.18,              # Cys (usually not charged)
    "Y": 10.07,             # Tyr (rarely charged)
    "K": 10.53, "R": 12.48, # basic
}

# N-terminal pKa ~9.69, C-terminal pKa ~2.34
NTERM_PKA = 9.69
CTERM_PKA = 2.34

# Deamidation-prone motifs (Asn followed by small residues)
DEAMIDATION_MOTIFS = {"NG", "NS", "NT", "NA", "NH", "ND"}

# Oxidation-liable residues
OXIDATION_RESIDUES = {"M", "W", "C"}

# Isomerization-prone motifs (Asp followed by small residues)
ISOMERIZATION_MOTIFS = {"DG", "DS", "DT", "DA"}


@dataclass
class DevelopabilityProfile:
    """Full developability assessment for a binder sequence."""
    sequence: str
    length: int

    # Aggregation
    max_hydrophobic_patch: int  # longest consecutive hydrophobic stretch
    mean_hydrophobicity: float
    hydrophobic_fraction: float  # fraction of I/V/L/F/A/M/W

    # PTM liabilities
    deamidation_sites: list = field(default_factory=list)  # positions of NG/NS etc
    oxidation_sites: list = field(default_factory=list)     # positions of M/W/C
    isomerization_sites: list = field(default_factory=list)  # positions of DG/DS etc
    ptm_liability_count: int = 0

    # Charge / solubility
    net_charge_ph7: float = 0.0
    isoelectric_point: float = 0.0
    charge_density: float = 0.0  # |net charge| / length
    pos_charge_count: int = 0
    neg_charge_count: int = 0

    # Complexity
    shannon_entropy: float = 0.0
    unique_aa_count: int = 0

    # Summary flags
    flags: list = field(default_factory=list)

    @property
    def ptm_liability_density(self) -> float:
        return self.ptm_liability_count / max(self.length, 1)


def compute_net_charge(seq: str, ph: float = 7.0) -> float:
    """Compute net charge at given pH using Henderson-Hasselbalch."""
    charge = 0.0
    # N-terminal (positive when protonated)
    charge += 1.0 / (1.0 + 10 ** (ph - NTERM_PKA))
    # C-terminal (negative when deprotonated)
    charge -= 1.0 / (1.0 + 10 ** (CTERM_PKA - ph))

    for aa in seq:
        if aa in ("K", "R", "H"):
            # Basic: positive when protonated
            charge += 1.0 / (1.0 + 10 ** (ph - PKA[aa]))
        elif aa in ("D", "E"):
            # Acidic: negative when deprotonated
            charge -= 1.0 / (1.0 + 10 ** (PKA[aa] - ph))
    return charge


def compute_isoelectric_point(seq: str) -> float:
    """Binary search for isoelectric point (pH where net charge = 0)."""
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        charge = compute_net_charge(seq, mid)
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def find_hydrophobic_patches(seq: str, threshold: float = 1.5) -> list[tuple[int, int, int]]:
    """Find stretches of consecutive hydrophobic residues.

    Returns list of (start, end, length) for patches where average
    hydrophobicity exceeds threshold.
    """
    patches = []
    start = None
    run_sum = 0.0
    run_len = 0

    for i, aa in enumerate(seq):
        h = KD_HYDROPHOBICITY.get(aa, 0)
        if h >= threshold:
            if start is None:
                start = i
                run_sum = 0
                run_len = 0
            run_sum += h
            run_len += 1
        else:
            if start is not None and run_len >= 5:
                patches.append((start, i - 1, run_len))
            start = None
            run_len = 0
            run_sum = 0

    if start is not None and run_len >= 5:
        patches.append((start, len(seq) - 1, run_len))

    return patches


def find_ptm_liabilities(seq: str) -> tuple[list, list, list]:
    """Scan for post-translational modification liability motifs."""
    deamidation = []
    oxidation = []
    isomerization = []

    for i in range(len(seq) - 1):
        dipeptide = seq[i:i + 2]
        if dipeptide in DEAMIDATION_MOTIFS:
            deamidation.append(i)
        if dipeptide in ISOMERIZATION_MOTIFS:
            isomerization.append(i)

    for i, aa in enumerate(seq):
        if aa in OXIDATION_RESIDUES:
            oxidation.append(i)

    return deamidation, oxidation, isomerization


def compute_shannon_entropy(seq: str) -> float:
    """Shannon entropy of amino acid composition."""
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


def assess_developability(seq: str) -> DevelopabilityProfile:
    """Run full developability assessment on a binder sequence."""
    seq = seq.upper().strip()
    n = len(seq)

    # Hydrophobicity
    hydrophobic_set = set("IVLFAMW")
    hydrophobic_fraction = sum(1 for aa in seq if aa in hydrophobic_set) / max(n, 1)
    mean_hydro = sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in seq) / max(n, 1)
    patches = find_hydrophobic_patches(seq)
    max_patch = max((p[2] for p in patches), default=0)

    # PTM liabilities
    deamidation, oxidation, isomerization = find_ptm_liabilities(seq)
    ptm_total = len(deamidation) + len(oxidation) + len(isomerization)

    # Charge
    charge_ph7 = compute_net_charge(seq, 7.0)
    pi = compute_isoelectric_point(seq)
    pos_count = sum(1 for aa in seq if aa in "KRH")
    neg_count = sum(1 for aa in seq if aa in "DE")

    # Complexity
    entropy = compute_shannon_entropy(seq)
    unique_aa = len(set(seq))

    # Build flags
    flags = []
    if entropy < 2.0:
        flags.append("low_complexity")
    if max_patch >= 7:
        flags.append(f"hydrophobic_patch_{max_patch}aa")
    if hydrophobic_fraction > 0.5:
        flags.append("high_hydrophobicity")
    if abs(charge_ph7) > 10:
        flags.append(f"extreme_charge_{charge_ph7:+.1f}")
    if pi < 4.0 or pi > 10.0:
        flags.append(f"extreme_pI_{pi:.1f}")
    if ptm_total > 5:
        flags.append(f"high_ptm_liability_{ptm_total}")
    if len(deamidation) >= 3:
        flags.append(f"deamidation_prone_{len(deamidation)}_sites")
    if unique_aa <= 5:
        flags.append("very_low_aa_diversity")

    return DevelopabilityProfile(
        sequence=seq,
        length=n,
        max_hydrophobic_patch=max_patch,
        mean_hydrophobicity=round(mean_hydro, 3),
        hydrophobic_fraction=round(hydrophobic_fraction, 3),
        deamidation_sites=deamidation,
        oxidation_sites=oxidation,
        isomerization_sites=isomerization,
        ptm_liability_count=ptm_total,
        net_charge_ph7=round(charge_ph7, 2),
        isoelectric_point=pi,
        charge_density=round(abs(charge_ph7) / max(n, 1), 4),
        pos_charge_count=pos_count,
        neg_charge_count=neg_count,
        shannon_entropy=round(entropy, 3),
        unique_aa_count=unique_aa,
        flags=flags,
    )
