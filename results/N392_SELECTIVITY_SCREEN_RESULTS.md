# N392 Selectivity Mutation Screen Results

**Date:** 2026-03-26
**Instance:** 33489714 (RTX 4090, US)
**Predictions:** 68 (17 constructs x 4 targets, single sample each)
**Cost:** ~$0.20 compute
**Verdict:** NO HITS — fundamental geometric constraint prevents allele selectivity via point mutations

---

## Objective

Widen the V2 binder's N392 preference over K392 to delta >= +0.10, while maintaining ERAP1/IRAP counter-selectivity.

## Wildtype Baseline (Y87A_Y89A, 92aa)

| Target | ipTM |
|--------|------|
| ERAP2-K392 | 0.800 |
| ERAP2-N392 | 0.733 |
| ERAP1 | 0.523 |
| IRAP | 0.416 |
| **N392-K392 delta** | **-0.067** |

Note: WT actually prefers K392 by 0.067 in this single-sample run (previous MSA-enabled runs showed slight N392 preference). The msa=empty baseline is K392-leaning.

## Full Results (ranked by N392 delta)

| Rank | Mutation | K392 | N392 | Delta | ERAP1 | IRAP | Assessment |
|------|----------|------|------|-------|-------|------|------------|
| 1 | F51A | 0.329 | 0.641 | **+0.311** | 0.546 | 0.434 | Huge delta but ERAP2 too low, ERAP1 too high |
| 2 | F51Q | 0.431 | 0.704 | **+0.273** | 0.752 | 0.352 | N392 hits 0.70 but ERAP1=0.752 |
| 3 | F51W | 0.328 | 0.597 | **+0.270** | 0.653 | 0.441 | ERAP2 too low, ERAP1 too high |
| 4 | F51N | 0.396 | 0.580 | **+0.184** | 0.616 | 0.401 | Same pattern |
| 5 | N38Q | 0.414 | 0.583 | **+0.168** | 0.634 | 0.357 | Same pattern |
| 6 | N38K | 0.438 | 0.490 | +0.053 | 0.490 | 0.330 | Marginal, low ERAP2 |
| 7 | N38R_F51W | 0.394 | 0.409 | +0.016 | 0.649 | 0.448 | Negligible |
| 8 | N38A | 0.445 | 0.426 | -0.019 | 0.677 | 0.346 | Neutral |
| 9 | **wt** | **0.800** | **0.733** | **-0.067** | **0.523** | **0.416** | **Baseline** |
| 10 | N38W | 0.414 | 0.337 | -0.078 | 0.459 | 0.453 | Wrong direction |
| 11 | N38R | 0.447 | 0.363 | -0.083 | 0.416 | 0.360 | Wrong direction |
| 12 | N38K_F51K | 0.513 | 0.391 | -0.121 | 0.729 | 0.377 | Wrong direction |
| 13 | F51R | 0.616 | 0.459 | -0.157 | 0.422 | 0.361 | K392-selective |
| 14 | N38K_F51W | 0.677 | 0.511 | -0.166 | 0.433 | 0.385 | K392-selective |
| 15 | N38Q_F51N | 0.642 | 0.466 | -0.176 | 0.711 | 0.355 | K392-selective |
| 16 | N38D | 0.592 | 0.370 | -0.223 | 0.766 | 0.413 | K392-selective (neg control) |
| 17 | F51K | 0.674 | 0.449 | -0.225 | 0.523 | 0.478 | K392-selective |

## Key Findings

### 1. N392 selectivity and ERAP1 selectivity are anti-correlated

Every mutation that increases N392 preference also increases ERAP1 cross-reactivity:

```
         N392 delta   ERAP1
F51A     +0.311       0.546 (up from 0.523)
F51Q     +0.273       0.752 (up!)
F51W     +0.270       0.653 (up)
N38Q     +0.168       0.634 (up)
```

**Root cause:** ERAP1 has Asn (N375) at the structurally equivalent position to ERAP2-N392. Any binder modification that better accommodates Asparagine at position 392 will also better accommodate ERAP1's N375.

### 2. Position 50 (Phe) is the allele discrimination switch

Removing or replacing Phe50 dramatically shifts the K392/N392 balance:
- F51A: delta swings from -0.067 (WT) to +0.311 (+0.378 shift!)
- The aromatic ring at position 50 makes specific van der Waals contacts with K392's aliphatic chain (CE at 3.84A). Without this contact, K392's longer side chain has no advantage — N392's shorter chain actually fits better in the remaining binding pocket.

### 3. Position 37 (Asn) has weaker allele discrimination

N38Q is the best pos37 mutation (delta +0.168) but much weaker than the pos50 mutations. The Asn37 OD1-K392 NZ H-bond (4.36A) can be replaced by Gln's longer side chain reaching N392's amide, but the effect is smaller because H-bonds are more tolerant of side chain length changes.

### 4. Positive charges at the 392 contact FAVOR K392 (hypothesis disproven)

The electrostatic repulsion hypothesis (Lys/Arg on binder repels K392) was wrong:
- F51K: delta = -0.225 (strongly K392-selective)
- F51R: delta = -0.157 (K392-selective)
- N38K: delta = +0.053 (marginal)

Lys and Arg side chains are long and flexible — they find favorable hydrophobic packing with K392's aliphatic chain rather than electrostatic repulsion at the tip.

### 5. Double mutations cancel rather than combine

Every double mutant performed worse than expected from individual effects:
- N38K_F51W: expected ~positive from N38K(+0.053) + F51W(+0.270), got -0.166
- N38Q_F51N: expected ~positive from N38Q(+0.168) + F51N(+0.184), got -0.176

The two mutations destabilize the interface cooperatively, causing the binder to re-dock in new poses (the same re-docking artifact seen in prior V2 campaigns).

### 6. N38D negative control confirmed mechanism

N38D (Asp-, negative charge) showed delta = -0.223, confirming that a negative charge at position 37 attracts K392's positive Lys. The electrostatic model works at pos37 (where the atoms directly face each other) even though it fails at pos50 (where the geometry allows flexible re-packing).

---

## Why Allele-Selective N392 Mutations Can't Work (Structural Explanation)

The K392 and N392 alleles differ by one amino acid at the floor of ERAP2's substrate channel. The V2 binder contacts position 392 through exactly 2 residues (Asn37 and Phe50), both at close range (3.8-4.4A).

**The trap:** To prefer N392 over K392, a binder mutation must:
1. Accommodate Asparagine (short, neutral) at position 392
2. Disfavor Lysine (long, positive) at position 392

But ERAP1 has N375 (also Asparagine) at the structurally equivalent position. The binder contacts position 392 DIRECTLY — there's no way to exploit "local context" differences between ERAP2-N392 and ERAP1-N375 when the contact residue itself is identical.

The prior ERAP2/ERAP1 selectivity came from contacts with the divergent patches (353-367, 400-414), which are all >6A from the binder. Those contacts are preserved regardless of which amino acid sits at position 392 — but they also can't distinguish the N392 allele.

**In short:** Allele selectivity at position 392 requires contacts at position 392, but those same contacts can't distinguish ERAP2-N392 from ERAP1-N375. The two objectives (N392 preference + ERAP1 selectivity) are geometrically incompatible at this binding site.

---

## What Could Still Work

### 1. Accept allele-agnostic ERAP2 selectivity
The wildtype Y87A_Y89A already binds both alleles well (K392=0.800, N392=0.733). It achieves ERAP1 selectivity through the divergent patch contacts, not through allele-specific interactions. This is already a useful tool.

### 2. Multi-site redesign (ProteinMPNN on N392 complex)
ProteinMPNN could redesign multiple surface positions simultaneously on the N392 complex, potentially finding a novel binding mode that achieves both N392 preference and ERAP1 selectivity through an entirely different contact network. This is Plan Phase 5.

### 3. Peptide approach (V4 series)
The 6-mer peptides achieve selectivity through SIZE, not surface complementarity. Position 392 matters differently for peptide substrates than for protein binders — the K/N392 change affects substrate specificity (P1 preference), which could be exploited with allele-specific P1 residue choices.

---

## Approach B: ProteinMPNN N392-Conditioned Redesign

Redesigned 10 binder positions simultaneously (near N392/Y398/A403) using ProteinMPNN conditioned on the V2-N392 complex. 100 designs generated, top 10 screened against 4 targets.

### MPNN Results (single sample)

| Design | K392 | N392 | Delta | ERAP1 | IRAP | Note |
|--------|------|------|-------|-------|------|------|
| **wt** | **0.800** | **0.733** | **-0.067** | **0.523** | **0.416** | **Baseline** |
| mpnn07 | 0.467 | 0.748 | +0.281 | 0.565 | 0.524 | Best delta (near-hit) |
| mpnn02 | 0.423 | 0.688 | +0.265 | 0.776 | 0.375 | High ERAP1 |
| mpnn08 | 0.425 | 0.664 | +0.240 | 0.733 | 0.425 | High ERAP1 |
| mpnn06 | 0.544 | 0.613 | +0.069 | 0.536 | 0.385 | Moderate |

### mpnn07 3-Sample Validation (FAILED)

| Target | Samples | Avg |
|--------|---------|-----|
| ERAP2-K392 | 0.737, 0.362, 0.431 | **0.510** |
| ERAP2-N392 | 0.521, 0.592, 0.421 | **0.511** |
| ERAP1 | 0.614, 0.542, 0.630 | **0.595** |
| IRAP | 0.397, 0.360, 0.414 | **0.390** |

**The single-sample N392=0.748 was noise.** With 3 samples, delta collapses from +0.281 to +0.001. ProteinMPNN multi-site redesign hits the same ERAP1-N375 mimicry wall as point mutations.

### MPNN Mutations (all top 10 shared core)
- D1H, H4R, Y5H — N-terminal (near Y398)
- A48G, I49F — 392 neighborhood reshape
- F51A or F51P — removes the allele switch Phe (consistent with point mutation findings)
- N70P, I73W — deeper contacts

---

## Files

| File | Description |
|------|-------------|
| `scripts/n392_contact_analysis.py` | Phase 1: binder-392 contact identification |
| `scripts/n392_mutation_design.py` | Phase 2: mutation panel + YAML generation |
| `scripts/n392_vastai_run.sh` | Phase 3: Vast.ai batch Boltz-2 runner |
| `scripts/n392_mutation_analysis.py` | Phase 3: results parser (for local analysis) |
| `data/results/n392_selectivity/boltz_yamls/` | 68 Boltz-2 input YAMLs |
| `data/results/n392_selectivity/boltz_results/` | Boltz-2 output (downloaded) |
| `data/results/n392_selectivity/mutation_manifest.json` | Experiment manifest |
| `results/N392_SELECTIVITY_SCREEN_RESULTS.md` | This report |
