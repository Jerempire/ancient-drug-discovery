# Peptide Synthesis & Assay Specifications — ERAP2 Short-Stack Leads

**Project:** Ancient Drug Discovery — ERAP2-K392 Selective Peptide Inhibitors
**Date:** 2026-03-24
**Contact:** Jeremy Johnson
**Status:** PENDING — awaiting MD validation results before placing order

---

## 1. Peptide Synthesis Order

### Vendor Options (ranked by preference)
1. **GenScript** — Custom peptide synthesis, fast turnaround, competitive pricing
2. **CPC Scientific** — Specializes in short peptides, bulk discounts
3. **AAPPTEC** — Good for non-standard modifications (D-amino acids, Phase 2)
4. **Bachem** — Premium quality, higher cost

### Order Specifications

| Parameter | Specification |
|---|---|
| Form | Free peptide (no modifications) |
| Salt form | TFA (trifluoroacetate) |
| Purity | >95% by RP-HPLC |
| Quantity | 5 mg per peptide |
| QC deliverables | HPLC chromatogram + MALDI-TOF mass spectrum |
| Turnaround | Standard (2-3 weeks) |

### Peptide Panel

| # | Code Name | Sequence | 1-Letter | MW (Da) | pI (est) | Purpose |
|---|---|---|---|---|---|---|
| 1 | The Hammer | Val-Ala-Gly-Ser-Ala-Phe | VAGSAF | ~536 | 5.5 | Lead: highest K392 ipTM (0.905) |
| 2 | The Scalpel | Ile-Ala-Phe-Ser-Ala-Phe | IAFSAF | ~624 | 5.5 | Lead: best allele-selective (K392-N392 delta +0.239) |
| 3 | The Trp Edge | Val-Ala-Trp-Ser-Ala-Phe | VAWSAF | ~651 | 5.5 | Probe: does Trp add residence time? |
| 4 | Scrambled Ctrl | Phe-Ala-Ser-Gly-Ala-Val | FASGAV | ~536 | 5.5 | Negative control: VAGSAF residues, scrambled order |

### Solubility Notes
- All peptides are relatively hydrophobic (Ala, Val, Phe, Trp)
- Recommend dissolving in DMSO at 10 mM stock, then diluting into assay buffer
- Final DMSO concentration in assay should not exceed 1%
- If solubility issues arise with VAWSAF (Trp), try 10% acetonitrile in buffer

### Estimated Cost
- 4 peptides x 5 mg x >95% purity: **$400-800 total**
- GenScript typical pricing: $100-150/peptide at this spec

---

## 2. Enzymatic Assay Protocol — IC50 Determination

### Reagents

| Reagent | Source | Catalog # | Notes |
|---|---|---|---|
| Recombinant human ERAP2 | R&D Systems | 3735-ZN | Includes both K392 and N392 variants |
| Recombinant human ERAP1 | R&D Systems | 2334-ZN | Counterscreen |
| Recombinant human IRAP | R&D Systems | 2680-ZN | Counterscreen |
| L-Leucine-AMC (Leu-AMC) | Sigma-Aldrich | L2145 | Fluorogenic substrate |
| Bestatin | Sigma-Aldrich | B8385 | Positive control inhibitor |
| DMSO | Sigma-Aldrich | D2650 | Vehicle control |

### Assay Buffer
- 50 mM Tris-HCl, pH 7.4
- 150 mM NaCl
- 0.01% Triton X-100
- Prepare fresh, filter sterilize

### Assay Conditions
- **Temperature:** 37C
- **Plate format:** 96-well black, flat-bottom (Corning #3694)
- **Final volume:** 100 uL per well
- **Pre-incubation:** 15 min enzyme + inhibitor at 37C before adding substrate
- **Substrate concentration:** 100 uM Leu-AMC (above Km to ensure sensitivity)
- **Enzyme concentration:** Titrate to give linear signal for 30 min; typically 0.5-2 nM final

### Inhibitor Concentrations (8-point dose-response)
| Point | Concentration |
|---|---|
| 1 | 0.1 uM |
| 2 | 0.3 uM |
| 3 | 1 uM |
| 4 | 3 uM |
| 5 | 10 uM |
| 6 | 30 uM |
| 7 | 100 uM |
| 8 | 300 uM |

### Controls (every plate)
- **Vehicle:** DMSO-only (matching highest DMSO concentration)
- **Positive control:** 100 uM bestatin (known aminopeptidase inhibitor)
- **No-enzyme blank:** Substrate only (background fluorescence)
- **All in triplicate**

### Readout
- **Instrument:** Fluorescence plate reader (e.g., BioTek Synergy, Tecan Infinite)
- **Excitation:** 380 nm
- **Emission:** 460 nm
- **Mode:** Kinetic read, every 30 seconds for 30 minutes
- **Analysis:** Calculate initial rates (slope of linear phase, typically first 10 min)

### Data Analysis
1. Plot initial rate vs log[inhibitor]
2. Fit to 4-parameter logistic (Hill equation):
   ```
   v = v_min + (v_max - v_min) / (1 + (IC50/[I])^n)
   ```
3. Report: IC50 +/- 95% CI, Hill coefficient (n)
4. If IC50 < 10 uM, determine Ki using Cheng-Prusoff:
   ```
   Ki = IC50 / (1 + [S]/Km)
   ```
   Literature Km values: ERAP2 ~150 uM for Leu-AMC

---

## 3. Selectivity Panel

Run the IC50 assay for each peptide against all 4 targets:

| Peptide | ERAP2-K392 | ERAP2-N392 | ERAP1 | IRAP |
|---|---|---|---|---|
| VAGSAF | PRIMARY | Counterscreen | Counterscreen | Counterscreen |
| IAFSAF | PRIMARY | Counterscreen | Counterscreen | Counterscreen |
| VAWSAF | PRIMARY | Counterscreen | Counterscreen | Counterscreen |
| FASGAV | Control | Control | Control | Control |

### Minimum Selectivity Criteria
- **10-fold selectivity**: IC50(ERAP1) / IC50(ERAP2-K392) >= 10
- **FASGAV negative control**: IC50 > 500 uM on all targets
- **Allele selectivity** (IAFSAF): IC50(N392) / IC50(K392) >= 3

### Expected Outcomes

| Peptide | ERAP2-K392 | ERAP2-N392 | ERAP1 | IRAP |
|---|---|---|---|---|
| VAGSAF | < 50 uM | < 100 uM | > 500 uM | > 500 uM |
| IAFSAF | < 50 uM | > 200 uM | > 500 uM | > 500 uM |
| VAWSAF | < 50 uM | < 100 uM | > 500 uM | > 500 uM |
| FASGAV | > 500 uM | > 500 uM | > 500 uM | > 500 uM |

### Decision Matrix

| IC50 Range | Interpretation | Next Step |
|---|---|---|
| < 10 uM | Strong lead | Proceed to D-amino acid optimization + SPR |
| 10-50 uM | Moderate lead | Optimize P1/P3 positions, test analogs |
| 50-100 uM | Weak but proves mechanism | Structure-guided optimization needed |
| > 100 uM | Boltz-2 predictions unreliable | Reassess computational pipeline |
| > 500 uM | No activity | Peptide does not inhibit this target |

---

## 4. Phase 2 — D-Amino Acid Variants (contingent on Phase 1 results)

Only order if L-amino acid IC50 < 50 uM for at least one lead.

| # | Name | Modification | Estimated Cost |
|---|---|---|---|
| 5 | D-VAGSAF | All D-amino acids | +$50-100 surcharge |
| 6 | D-IAFSAF | All D-amino acids | +$50-100 surcharge |
| 7 | RI-VAGSAF | Retro-inverso (D-FASGA-v) | +$50-100 surcharge |

Purpose: Protease resistance for in vivo stability. D-peptides resist serum proteases but must be re-validated for binding.

---

## 5. Phase 2 — SPR Binding Kinetics (contingent on IC50 confirmation)

| Parameter | Specification |
|---|---|
| Instrument | Biacore T200 or Cytiva 8K |
| Chip | CM5 (amine coupling) |
| Ligand (immobilized) | ERAP2 (amine coupling, target ~5000 RU) |
| Analyte | Peptides, 0.1 - 100 uM, 2-fold serial dilution |
| Buffer | HBS-EP+ (10 mM HEPES, 150 mM NaCl, 3 mM EDTA, 0.05% P20) |
| Temperature | 25C |
| Contact time | 120 s |
| Dissociation time | 300 s |
| Regeneration | 10 mM glycine pH 2.0 (if needed) |

### Key Metrics
- **KD** (equilibrium dissociation constant)
- **kon** (association rate)
- **koff** (dissociation rate)
- **Residence time** = 1/koff — critical for VAWSAF (Trp) vs VAGSAF comparison

---

## 6. Budget Summary

| Item | Cost |
|---|---|
| Peptide synthesis (4 x 5 mg, >95%) | $400-800 |
| Recombinant enzymes (ERAP2, ERAP1, IRAP) | $500-800 |
| Substrates + reagents | $100-200 |
| Plate reader time (if outsourced) | $0-500 |
| Provisional patent filing | $75 |
| **Total Phase 1** | **$1,075 - $2,375** |
| D-amino acid peptides (Phase 2) | +$150-300 |
| SPR (Phase 2, outsourced) | +$500-1,500 |
| **Total with Phase 2** | **$1,725 - $4,175** |

---

## 7. Timeline

| Week | Activity |
|---|---|
| Week 0 | Complete MD validation, place peptide order |
| Week 2-3 | Receive synthesized peptides |
| Week 3-4 | Run IC50 assays (ERAP2-K392 primary) |
| Week 4-5 | Run selectivity panel (ERAP1, IRAP, N392) |
| Week 5-6 | Analyze results, decide on D-amino acid order |
| Week 8-10 | SPR kinetics (if IC50 confirms activity) |

---

## Notes

1. **ERAP2 variant-specific enzyme availability is uncertain.** R&D Systems 3735-ZN may contain mixed K392/N392. If variant-specific recombinant enzyme is unavailable commercially, consider ordering from Addgene or custom expression.

2. **If outsourcing the assay to a CRO**, recommended: Reaction Biology, Eurofins DiscoverX, or BPS Bioscience. Provide this document as the assay specification.

3. **All peptides are L-amino acids in Phase 1.** Do not order D-amino acid versions until enzymatic data confirms the mechanism.

4. **Store peptides at -20C, dessicated.** Aliquot upon receipt to avoid freeze-thaw. Reconstitute in DMSO immediately before use.
