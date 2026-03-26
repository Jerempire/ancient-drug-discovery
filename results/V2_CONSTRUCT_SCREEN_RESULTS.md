# V2 Construct Screen Results — 2026-03-25

## Objective
Design potency upgrades for the 92aa ERAP2-selective binder (n248_trim_c5_Y87A_Y89A) to achieve 40-50% ERAP2 inhibition.

## Round 1: Binder Modifications (Constructs B/C/D)

### Constructs Tested
- **Parent**: n248_trim_c5_Y87A_Y89A (92aa, unmodified control)
- **B** (V19W+I49W): Double tryptophan plug at channel entrance + near zinc
- **C** (GGGGS cargo): N-term VAGSAF tethered via (GGGGS)x3 flexible linker
- **D** (combined): V19W+I49W + N-term VAGSAF via (GGGGS)x3

### Results

| Construct | ERAP2 avg ipTM | ERAP1 avg ipTM | IRAP avg ipTM | Selective? |
|-----------|---------------|---------------|--------------|------------|
| Parent    | 0.572         | 0.604         | 0.407        | NO         |
| B         | 0.507         | 0.623         | 0.370        | NO         |
| C         | 0.507         | 0.518         | 0.424        | NO         |
| D         | 0.455         | 0.616         | 0.376        | NO         |

### Key Finding
**The parent binder itself is not ERAP2-selective** (ERAP1 > ERAP2). All modifications worsened this. The binder binds a conserved surface region.

## Round 2: Partial Results (Instance died mid-run)

- **B2 (V19W only)**: ERAP2 0.519, ERAP1 0.690, IRAP 0.415 — worse than parent
- Confirmed: V19W mutation at channel entrance increases ERAP1 cross-reactivity
- C2 (EAAAK linker) results lost — rerun in VAGSAF screen below

## VAGSAF Stabilization Screen

### Variants Tested
1. **vagsaf_ctrl** (VAGSAF, 6aa): Original peptide control
2. **vagsaf_ss_cyc** (CVAGSAFC, 8aa): Disulfide-bracketed cyclization attempt
3. **vagsaf_ht_cyc** (GVAGSAFG, 8aa): Head-to-tail cyclization via bond constraint
4. **vagsaf_pro** (PVAGSAFP, 8aa): Proline-capped termini for protease resistance
5. **vagsaf_ext** (VAGSAFLL, 8aa): C-terminal leucine extension
6. **C2 EAAAK** (113aa): V19W binder + VAGSAF via rigid (EAAAK)x3 helical linker

### Results

| Variant | ERAP2 avg | ERAP1 avg | IRAP avg | E2>E1? | E2>IR? |
|---------|-----------|-----------|----------|--------|--------|
| **C2 EAAAK** | **0.759** | **0.193** | **0.609** | **YES (3.9x)** | **YES (1.2x)** |
| vagsaf_ext | 0.817 | 0.599 | 0.804 | YES | barely |
| vagsaf_ss_cyc | 0.772 | 0.755 | 0.811 | barely | NO |
| vagsaf_ctrl | 0.763 | 0.565 | 0.782 | YES | NO |
| vagsaf_pro | 0.742 | 0.589 | 0.806 | YES | NO |
| vagsaf_ht_cyc | 0.711 | 0.502 | 0.897 | YES | NO |

### Key Finding
**C2 EAAAK is the only construct with dual selectivity** (ERAP2 >> ERAP1, ERAP2 > IRAP).

## Structural Analysis

| Variant | Shape | Rg (Å) | N-C dist (Å) | S-S dist (Å) |
|---------|-------|--------|--------------|--------------|
| vagsaf_ctrl | Extended/linear | 5.5 | 17.8 | - |
| vagsaf_ext | Extended/linear | 6.3-7.4 | 19.7-23.6 | - |
| vagsaf_pro | Extended/linear | 5.1-5.3 | 16.4-16.8 | - |
| vagsaf_ht_cyc | **CYCLIC (ring)** | 4.0 | **1.3-1.4** | - |
| vagsaf_ss_cyc | Mostly open | 4.5-6.5 | 5.7-20.6 | 2.0-22.1 |
| C2 EAAAK | Compact/globular | 16-20 | 31-62 | - |

- Head-to-tail cyclization worked (bond constraint) but killed selectivity
- Disulfide cyclization failed in 2/3 models (no explicit constraint)
- Extended/linear peptides bind substrate groove but aren't IRAP-selective
- Globular C2 binder provides selective anchoring

## Critical Insights

1. **VAGSAF alone is NOT selective against IRAP** — all free peptide variants show IRAP binding >= ERAP2. The original triple selectivity (ipTM 0.905/0.335/0.236) was under different conditions.

2. **The binder provides selectivity, not potency** — the binder's ERAP2-specific surface anchoring prevents the tethered VAGSAF from docking into IRAP's conserved groove.

3. **Rigid linker is essential** — GGGGS (flexible) allows cargo to sample conserved surfaces on all enzymes. EAAAK (rigid helix) constrains cargo to the ERAP2 channel only.

4. **Architecture that works**: Globular ERAP2-selective binder → rigid EAAAK helical arm → VAGSAF substrate mimic at tip

## Winner: C2 EAAAK

```
Sequence (113aa):
VAGSAFEAAAKEAAAKEAAAKDIRHYFKSLEEYLKNLPKWVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN
```

- ERAP2 ipTM: 0.759 avg (best: 0.810)
- ERAP1 ipTM: 0.193 avg (3.9x selectivity)
- IRAP ipTM: 0.609 avg (1.2x selectivity)
- 113 amino acids, single chain, standard expression

## Next Steps
1. Validate C2 EAAAK with higher diffusion samples (5-10)
2. Test linker length variants: (EAAAK)x2 and (EAAAK)x4
3. Molecular dynamics on C2-ERAP2 complex (stability, residence time)
4. Expression feasibility assessment (the EAAAK helix should express well in E. coli)

## Files
- Round 1 results: `results/construct_bc_screen/`
- VAGSAF screen results: `results/vagsaf_screen/`
- Design script: `scripts/construct_bc_design.py`
- Vast.ai instances: 33420112 (died), 33548811 (destroyed after completion)
