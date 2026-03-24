# ERAP2 Peptide Synthesis Panel

**Locked:** 2026-03-23
**Target:** ERAP2 K392/N392 allele-selective substrate-mimic peptides
**Budget:** ~$1,200 (4 peptides) or ~$900 (compact 3-peptide set)

## Final Panel

| # | Name | Sequence | Len | P1 Residue | Arm | Source | Role |
|---|------|----------|-----|------------|-----|--------|------|
| 1 | pep_glu_long_01 | EALVAAGLAGLA | 12 | E (Glutamate) | K392 | Hand-designed | Proof-of-concept K392 lead |
| 2 | hybrid_E_VKLLLL | EKLLLLSIGK | 10 | E (Glutamate) | K392 | PepMLM hybrid | Robust scaffold K392 |
| 3 | pep_ala_01 | AALVAAGLA | 9 | A (Alanine) | N392 | Hand-designed | N392 lead |
| 4 | pep_leu_01 | LALVAAGLA | 9 | L (Leucine) | N392 | Hand-designed | N392 comparator |

**Compact set (if only 3):** Drop #4 (pep_leu_01), keep 1-3. Saves ~$300; loses the N392 comparator but retains both arms.

## Design Rationale

- **K392 arm** uses negatively charged P1 (Glu) to form a salt bridge with K392's NH3+ — the core selectivity mechanism validated by DiffPepDock
- **N392 arm** uses small/hydrophobic P1 (Ala, Leu) matching the 165x published hydrophobic preference for the ancestral allele
- **Two scaffold families** (hand-designed ALVAAG backbone vs PepMLM-generated KLLLL backbone) provide independent confirmation if both bind

## Cost Estimates

| Item | Est. Cost |
|------|-----------|
| pep_glu_long_01 (12-mer) | ~$300-360 |
| hybrid_E_VKLLLL (10-mer) | ~$250-300 |
| pep_ala_01 (9-mer) | ~$225-270 |
| pep_leu_01 (9-mer) | ~$225-270 |
| **Total (4 peptides)** | **~$975-1,170** |

Pricing based on $3.50-4.50/residue for standard L-amino acid peptides at >95% purity, 5 mg scale.

## Vendor Notes

- **GenScript** — Standard peptide synthesis, competitive pricing at scale, reliable turnaround (~2-3 weeks). Get quote via web portal.
- **AAPPTEC** — Lower per-residue cost ($3.50-4.00), good for simple sequences without modifications. Longer turnaround.
- Both vendors handle 9-12mer peptides routinely. No modifications needed (all standard L-amino acids).
- D-amino acid versions (non-cleavable) are a future option at ~1.5-2x cost premium.

## Evidence Summary

### DiffPepDock (V4 docking, 28 peptides)
- **pep_glu_long_01**: Rank 1, delta = -64.9 (strongest K392 selectivity)
- **pep_leu_01**: Rank 17, delta = +43.1 (N392-selective, model match)
- **pep_ala_01**: Rank 20, delta = +44.4 (N392-selective)
- Salt bridge mechanism validated: 5/13 charged P1 peptides are K392-selective
- Full report: `v4/v4_docking_report.md`

### PepMLM Hybrid Boltz-2 Validation
- **hybrid_E_VKLLLL** (K392): mean ipTM = 0.795, best = 0.818 (3 diffusion samples)
- **hybrid_E_VKLLLL** (N392): mean ipTM = 0.709, best = 0.804
- K392 preference confirmed (higher mean ipTM on target allele)
- Top performer among all 14 PepMLM hybrids tested
- Full data: `data/results/pepmlm_hybrids/boltz_results/hybrid_summary.json`

### Bias Check
- Results: `data/results/bias_check_results.json`

### SAR Grid
- Systematic P1 substitution analysis across scaffold variants
- Results: `v4/pepmlm/sar_results/sar_analysis.json`

## Next Steps

1. **File provisional patent** ($75 micro entity) — 12+ claims, 5+ SEQ IDs, PepMLM novelty evidence, two independent scaffold confirmations
2. **Order peptide synthesis** (~$1,200 for 4 peptides) — get quotes from GenScript + AAPPTEC
3. **Optional:** PepINVENT (AstraZeneca, open source) for RL-based D-amino acid optimization on Vast.ai
4. **Optional:** V2 protein binder synthesis ($3K) — deprioritized vs peptides
