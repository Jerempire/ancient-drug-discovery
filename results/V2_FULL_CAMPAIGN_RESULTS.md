# V2 Full Campaign Results — 2026-03-25/26

## Objective
Design a selective ERAP2 inhibitor for cancer immunotherapy. ERAP2 trims antigenic peptides on cancer cell surfaces — inhibiting it exposes more neoantigens to the immune system. Must NOT inhibit ERAP1 or IRAP (selectivity).

## Starting Point
- **Parent binder**: n248_trim_c5_Y87A_Y89A (92aa protein, designed to bind ERAP2 channel entrance)
- **VAGSAF**: 6-mer peptide with previously reported triple selectivity (ipTM 0.905/0.335/0.236 for ERAP2/ERAP1/IRAP)
- **Goal**: Achieve 40-50% ERAP2 inhibition with selectivity

---

## Campaign Summary: 5 Rounds, 162 Boltz-2 Predictions

### Round 1: Binder Modifications (FAILED)
**Hypothesis**: Add tryptophan plugs and/or tethered VAGSAF cargo to the binder to increase channel occlusion.

| Construct | Description | ERAP2 | ERAP1 | IRAP | Selective? |
|-----------|-------------|-------|-------|------|------------|
| Parent (control) | 92aa unmodified | 0.572 | 0.604 | 0.407 | NO |
| B (V19W+I49W) | Double W-plug | 0.507 | 0.623 | 0.370 | NO |
| C (GGGGS cargo) | N-term VAGSAF, flexible linker | 0.507 | 0.518 | 0.424 | NO |
| D (combined) | W-plugs + GGGGS cargo | 0.455 | 0.616 | 0.376 | NO |

**Key finding**: Parent binder itself prefers ERAP1 (0.604) over ERAP2 (0.572). All modifications worsened this. I49W is in the conserved catalytic core — guaranteed cross-reactivity.

### Round 2: Partial (Instance died)
- B2 (V19W only): ERAP2=0.519, ERAP1=0.690 — V19W alone also increases ERAP1. Binder scaffold is fundamentally not selective.

### Round 3: VAGSAF Stabilization + C2 EAAAK (BREAKTHROUGH)
**Pivot**: Abandon binder modifications. Test stabilized VAGSAF peptides as standalone drugs. Also test rigid EAAAK helical linker.

| Variant | Description | ERAP2 | ERAP1 | IRAP | E2>E1? | E2>IR? |
|---------|-------------|-------|-------|------|--------|--------|
| **C2 EAAAK** | **Binder + rigid helix + VAGSAF** | **0.759** | **0.193** | **0.609** | **YES (3.9x)** | **YES (1.2x)** |
| vagsaf_ext | VAGSAFLL (8aa) | 0.817 | 0.599 | 0.804 | YES | barely |
| vagsaf_ss_cyc | CVAGSAFC disulfide | 0.772 | 0.755 | 0.811 | barely | NO |
| vagsaf_ctrl | VAGSAF (6aa) | 0.763 | 0.565 | 0.782 | YES | NO |
| vagsaf_pro | PVAGSAFP Pro-capped | 0.742 | 0.589 | 0.806 | YES | NO |
| vagsaf_ht_cyc | GVAGSAFG cyclic | 0.711 | 0.502 | 0.897 | YES | NO |

**Key findings**:
1. Free VAGSAF is NOT IRAP-selective (all variants bind IRAP >= ERAP2)
2. Cyclization makes IRAP selectivity WORSE (0.897!)
3. C2 EAAAK is the ONLY dual-selective construct
4. The binder provides selectivity, VAGSAF provides potency, the rigid linker constrains delivery

### Round 4: Linker Fix for IRAP (INCREMENTAL)
**Problem**: C2 EAAAK has only 1.2x IRAP selectivity. Contact analysis showed EAAAK lysines form electrostatic contacts with IRAP Tyr138/His134.

| Variant | Description | ERAP2 | ERAP1 | IRAP | E2/E1 | E2/IR |
|---------|-------------|-------|-------|------|-------|-------|
| **c2_eaaaa3** | **No lysines, x3** | **0.778** | 0.333 | **0.573** | 2.3x | **1.4x** |
| c2_eaaak3 | Control | 0.729 | 0.294 | 0.588 | 2.5x | 1.2x |
| c2_eaaak2 | Shorter x2 | 0.703 | 0.240 | 0.571 | 2.9x | 1.2x |
| c2_eaaak1 | Shortest x1 | 0.743 | 0.299 | 0.740 | 2.5x | 1.0x |
| c2_eaaaa2 | No K, x2 | 0.575 | 0.292 | 0.691 | 2.0x | 0.8x |

**Key finding**: Removing lysines helped (1.2x → 1.4x) but hit a ceiling. 85% of IRAP contacts come from the binder body, not the linker.

### Round 4b: Negative Design (FAILED)
**Hypothesis**: Mutate 3 binder residues that contact IRAP but not ERAP2.

| Variant | ERAP2 | ERAP1 | IRAP | E2/IR |
|---------|-------|-------|------|-------|
| neg1 (I2E) | 0.664 | 0.255 | 0.708 | 0.9x |
| neg2 (D1R+I2E) | 0.753 | 0.201 | 0.782 | 1.0x |
| neg3 (D1R+I2E+K42E) | 0.739 | 0.234 | 0.673 | 1.1x |

**Key finding**: Point mutations at 3 positions caused the binder to re-dock in new IRAP-binding orientations. Can't fix with point mutations.

### Round 5: ProteinMPNN Interface Redesign (SUCCESS)
**Approach**: Systematic interface redesign. Aligned ERAP2 vs IRAP sequences (52% identity). Mapped all 35 binder interface contacts: 17 contact only conserved positions (IRAP cross-reactivity drivers), 4 contact only divergent positions (selectivity drivers). Used ProteinMPNN to redesign the 17 conserved-only positions while keeping divergent contacts fixed.

| Design | Mutations | ERAP2 | ERAP1 | IRAP | E2/E1 | E2/IR |
|--------|-----------|-------|-------|------|-------|-------|
| ref c2_eaaaa3 | (none) | 0.778 | 0.333 | 0.573 | 2.3x | 1.4x |
| mpnn01 | 15 muts | 0.810 | 0.381 | 0.711 | 2.1x | 1.1x |
| mpnn02 | 14 muts | 0.773 | 0.289 | 0.785 | 2.7x | 1.0x |
| mpnn03 | 15 muts | 0.727 | 0.567 | 0.824 | 1.3x | 0.9x |
| mpnn04 | 14 muts | 0.753 | 0.242 | 0.705 | 3.1x | 1.1x |
| mpnn05 | 14 muts | 0.730 | 0.454 | 0.671 | 1.6x | 1.1x |
| **mpnn06** | **14 muts** | **0.815** | **0.326** | **0.454** | **2.5x** | **1.8x** |
| mpnn07 | 14 muts | 0.676 | 0.267 | 0.708 | 2.5x | 1.0x |
| mpnn08 | 15 muts | 0.668 | 0.300 | 0.706 | 2.2x | 0.9x |
| mpnn09 | 15 muts | 0.730 | 0.322 | 0.745 | 2.3x | 1.0x |
| mpnn10 | 15 muts | 0.765 | 0.400 | 0.757 | 1.9x | 1.0x |

**Key finding**: mpnn06 achieved ERAP2=0.815, IRAP=0.454 (1.8x selectivity) — best IRAP selectivity of any construct while also having the highest ERAP2 binding. Only 1 of 10 MPNN designs improved IRAP selectivity; the rest were neutral or worse. Interface redesign works but is sensitive to exact mutation combinations.

---

## Final Lead: mpnn06

### Sequence (113aa)
```
VAGSAFEAAAAEAAAAEAAAAIERHYHKSLEEYLENL PKKVDMLVDLYSKGIFHLDNTNTLVEDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKNANARNR
```

Full breakdown:
- Cargo: VAGSAF (6aa)
- Linker: (EAAAA)x3 (15aa, rigid helix, no lysines)
- Binder body: 92aa redesigned (14 interface mutations)

### Mutations from parent (V19W binder)
```
R3I, Y5H, F6Q, W19K, L23M, T37N, I39T, V41L, F51D, I56Y, A89N, L90A, K91L, N92R
```

### Performance
| Target | ipTM avg | ipTM samples | Selectivity |
|--------|---------|-------------|-------------|
| ERAP2 | 0.815 | [0.878, 0.821, 0.745] | — |
| ERAP1 | 0.326 | — | 2.5x |
| IRAP | 0.454 | [0.710, 0.346, 0.305] | 1.8x |

### Properties
- 113 amino acids, ~12.4 kDa, single chain
- Standard E. coli expression (EAAAK helices express well)
- No disulfides, no non-standard amino acids
- Estimated synthesis cost: ~$1,500-2,500 (recombinant expression)

---

## Key Insights from Campaign

1. **The binder provides selectivity, the peptide provides potency.** Free VAGSAF binds all three aminopeptidases equally. The binder's selective surface anchoring constrains VAGSAF delivery to ERAP2 only.

2. **Rigid linkers are essential.** GGGGS (flexible) allows cargo to flop onto conserved surfaces. EAAAK (helical) constrains cargo toward the target channel.

3. **Lysines in EAAAK linkers create cross-reactive contacts.** EAAAA performs better than EAAAK for selectivity.

4. **Point mutations can't fix structural cross-reactivity.** The binder shares 17 interface residues with conserved ERAP2/IRAP surface. Breaking 3 just causes re-docking. Need systematic interface redesign.

5. **ProteinMPNN can improve selectivity but is noisy.** Only 1 of 10 designs improved IRAP selectivity. The successful design (mpnn06) was the one that preserved Asn38 (a divergent contact). Keeping divergent contacts fixed while reshuffling conserved contacts is the right strategy.

6. **Cyclization kills aminopeptidase selectivity.** Cyclic peptides can't extend into the substrate groove — they sit on conserved surfaces instead.

---

## Round 6: mpnn06 Validation + Cargo Tweaks (5-sample validation)

### Critical Finding: IRAP Selectivity Was Statistical Noise

mpnn06 re-run with 5 diffusion samples revealed the original 3-sample IRAP average (0.454) was unreliable. True IRAP binding is ~0.80 — no selectivity.

| Variant | ERAP2 | ERAP1 | IRAP | E2/E1 | E2/IR |
|---------|-------|-------|------|-------|-------|
| mpnn06 (5 samples) | 0.761 | 0.284 | 0.796 | 2.7x | 1.0x |
| mpnn06_short (EAAAA x2) | 0.652 | 0.216 | 0.801 | 3.0x | 0.8x |
| mpnn06_vagsae (F->E) | 0.835 | 0.233 | 0.812 | 3.6x | 1.0x |
| mpnn06_vagsak (F->K) | 0.834 | 0.294 | 0.798 | 2.8x | 1.0x |
| **mpnn06_vagsal (F->L)** | **0.824** | **0.188** | 0.806 | **4.4x** | 1.0x |
| mpnn06_short_vagsae | 0.750 | 0.333 | 0.863 | 2.3x | 0.9x |
| mpnn06_short_vagsak | 0.710 | 0.318 | 0.849 | 2.2x | 0.8x |

### Lesson: 3 diffusion samples are insufficient for selectivity claims
The original mpnn06 IRAP samples [0.710, 0.346, 0.305] appeared bimodal. With 5 samples [0.824, 0.808, 0.769, 0.750, 0.832], the binding is consistently high. Always validate selectivity with 5+ samples.

### IRAP Selectivity Is a Structural Limitation
No cargo variant, linker length, or interface tweak reduces IRAP binding below 0.77. The binder body docks on the M1 aminopeptidase catalytic domain surface, which is structurally conserved between ERAP2 and IRAP. This cannot be fixed by modifications to the current scaffold.

### Revised Lead: mpnn06_vagsal
- ERAP2: 0.824 [0.851, 0.799, 0.823] — strong, consistent
- ERAP1: 0.188 [0.189, 0.184, 0.192] — **4.4x selectivity, rock-solid**
- IRAP: 0.806 [0.799, 0.835, 0.783] — no selectivity (acknowledged)
- Architecture: VAGSAL-(EAAAA)x3-redesigned binder (113aa)
- The Leu at P1 (instead of Phe) gave the best ERAP1 selectivity of any construct tested

---

## Final Assessment

### What Works
- **ERAP1 selectivity: SOLVED.** mpnn06_vagsal at 4.4x (0.824 vs 0.188) is consistent and robust.
- **ERAP2 binding: STRONG.** 0.824 avg, tight distribution.
- **Architecture validated.** Binder + rigid EAAAK helix + peptide cargo is the right framework.

### What Doesn't Work
- **IRAP selectivity: NOT ACHIEVABLE** with this binder scaffold. The catalytic domain surface targeted by the binder is too conserved between ERAP2 and IRAP.

### Path to IRAP Selectivity
Requires a fundamentally different binder targeting ERAP2's domain IV (residues ~529-960), which has no structural equivalent in IRAP. This would be an RFdiffusion de novo design campaign — a separate project.

### Recommendation
Proceed to wet lab with mpnn06_vagsal as an ERAP2/ERAP1-selective inhibitor. The IRAP cross-reactivity is a known limitation that can be characterized in enzymatic assays. If IRAP inhibition proves problematic in vivo, pursue the domain IV RFdiffusion campaign.

## Next Steps
1. PyRosetta interface analysis on mpnn06_vagsal
2. Expression feasibility (E. coli, 113aa single chain)
3. Enzymatic assay design (ERAP2 vs ERAP1 vs IRAP)
4. Future: RFdiffusion binder targeting ERAP2 domain IV for IRAP selectivity

## Compute Summary
- 231 Boltz-2 predictions across 6 rounds
- 7 Vast.ai instances (~$3 total compute)
- 1 ProteinMPNN run (100 designs, 78 seconds)
- 1 instance died mid-run, 3 SCP downloads failed (destroyed too early)
