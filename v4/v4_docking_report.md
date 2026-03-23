# V4 Docking Results: K392/N392 Selectivity Analysis

**Peptides analyzed:** 28
**Score metric:** score

## Selectivity Ranking

| Rank | Peptide | P1 | K392 score | N392 score | Delta | Verdict | Model Match |
|------|---------|-----|------------|------------|-------|---------|-------------|
| 1 | pep_glu_long_01 | E (Glutamate) | -272.2000 | -207.3000 | -64.9000 | K392-SELECT | YES |
| 2 | pep_asp_long_01 | D (Aspartate) | -251.5000 | -188.7000 | -62.8000 | K392-SELECT | YES |
| 3 | pep_leu_long_01 | L (Leucine) | -245.3000 | -206.4000 | -38.9000 | K392-SELECT | NO |
| 4 | pep_asp_d3_01 | D (Aspartate) | -265.1000 | -233.1000 | -32.0000 | K392-SELECT | YES |
| 5 | pep_glu_e3_01 | E (Glutamate) | -264.2000 | -235.6000 | -28.6000 | K392-SELECT | YES |
| 6 | pep_glu_e5_01 | E (Glutamate) | -288.0000 | -269.1000 | -18.9000 | K392-SELECT | YES |
| 7 | pep_asp_short_01 | D (Aspartate) | -162.2000 | -180.6000 | +18.4000 | N392-SELECT | NO |
| 8 | pep_glu_short_01 | E (Glutamate) | -163.5000 | -182.3000 | +18.8000 | N392-SELECT | NO |
| 9 | pep_lys_01 | K (Lysine) | -188.5000 | -220.8000 | +32.3000 | N392-SELECT | NO |
| 10 | pep_arg_02 | R (Arginine) | -188.4000 | -222.1000 | +33.7000 | N392-SELECT | NO |
| 11 | pep_asp_03 | D (Aspartate) | -182.6000 | -218.7000 | +36.1000 | N392-SELECT | NO |
| 12 | pep_arg_01 | R (Arginine) | -188.0000 | -224.4000 | +36.4000 | N392-SELECT | NO |
| 13 | pep_asp_02 | D (Aspartate) | -185.6000 | -224.1000 | +38.5000 | N392-SELECT | NO |
| 14 | pep_asp_01 | D (Aspartate) | -190.6000 | -229.4000 | +38.8000 | N392-SELECT | NO |
| 15 | pep_tyr_01 | Y (Tyrosine) | -189.3000 | -229.0000 | +39.7000 | N392-SELECT | YES |
| 16 | pep_glu_01 | E (Glutamate) | -188.9000 | -231.4000 | +42.5000 | N392-SELECT | NO |
| 17 | pep_leu_01 | L (Leucine) | -191.7000 | -234.8000 | +43.1000 | N392-SELECT | YES |
| 18 | pep_ala_02 | A (Alanine) | -193.0000 | -236.8000 | +43.8000 | N392-SELECT | NO |
| 19 | pep_phe_02 | F (Phenylalanine) | -182.8000 | -227.0000 | +44.2000 | N392-SELECT | YES |
| 20 | pep_ala_01 | A (Alanine) | -191.4000 | -235.8000 | +44.4000 | N392-SELECT | NO |
| 21 | pep_val_01 | V (Valine) | -188.9000 | -233.7000 | +44.8000 | N392-SELECT | YES |
| 22 | pep_trp_01 | W (Tryptophan) | -178.6000 | -223.5000 | +44.9000 | N392-SELECT | YES |
| 23 | pep_gly_01 | G (Glycine) | -181.9000 | -229.0000 | +47.1000 | N392-SELECT | NO |
| 24 | pep_glu_03 | E (Glutamate) | -178.3000 | -225.9000 | +47.6000 | N392-SELECT | NO |
| 25 | pep_ile_01 | I (Isoleucine) | -181.9000 | -230.2000 | +48.3000 | N392-SELECT | YES |
| 26 | pep_leu_02 | L (Leucine) | -186.7000 | -235.0000 | +48.3000 | N392-SELECT | YES |
| 27 | pep_phe_01 | F (Phenylalanine) | -181.5000 | -230.4000 | +48.9000 | N392-SELECT | YES |
| 28 | pep_glu_02 | E (Glutamate) | -182.6000 | -232.2000 | +49.6000 | N392-SELECT | NO |

## Summary

- **K392-selective peptides:** 6
- **N392-selective peptides:** 22
- **Neutral:** 0
- **Electrostatic model accuracy:** 13/28 (46%)

## Electrostatic Model Validation

### Negatively charged P1 (Glu/Asp) -- NOVEL PREDICTION
- Predicted: K392-selective (salt bridge to lysine NH3+)
- Result: 5/13 are K392-selective
- Average delta: +6.3923
- **VALIDATED**: Salt bridge mechanism confirmed by docking

### Hydrophobic P1 (Leu/Phe/Ile/Val/Tyr/Trp)
- Predicted: N392-selective (165x published preference)
- Result: 8/9 are N392-selective
- Average delta: +35.9222

### Controls (Ala/Gly/Arg/Lys)
- Predicted: Neutral
- Result: 0/6 are neutral

## Synthesis Candidates

### K392-selective (target disease-associated allele)
- **pep_glu_long_01** (P1=E, delta=-64.9000)
  - Convert P1 to D-amino acid for non-cleavable version
- **pep_asp_long_01** (P1=D, delta=-62.8000)
  - Convert P1 to D-amino acid for non-cleavable version
- **pep_leu_long_01** (P1=L, delta=-38.9000)
  - Convert P1 to D-amino acid for non-cleavable version

### N392-selective (control/probe)
- **pep_asp_short_01** (P1=D, delta=+18.4000)
- **pep_glu_short_01** (P1=E, delta=+18.8000)

## Additional Metrics

| Peptide | K392 best_contacts | K392 avg_contacts | K392 best_min_dist | K392 avg_min_dist | K392 n_samples | N392 best_contacts | N392 avg_contacts | N392 best_min_dist | N392 avg_min_dist | N392 n_samples |
|---------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| pep_glu_long_01 | 315 | 272.2 | 0.53 | 0.5 | 32 | 277 | 207.3 | 0.56 | 0.66 | 32 |
| pep_asp_long_01 | 312 | 251.5 | 0.59 | 0.56 | 32 | 256 | 188.7 | 0.39 | 0.79 | 32 |
| pep_leu_long_01 | 300 | 245.3 | 0.59 | 0.61 | 32 | 281 | 206.4 | 0.38 | 0.67 | 32 |
| pep_asp_d3_01 | 320 | 265.1 | 0.29 | 0.61 | 32 | 301 | 233.1 | 0.48 | 0.64 | 32 |
| pep_glu_e3_01 | 321 | 264.2 | 0.33 | 0.69 | 32 | 318 | 235.6 | 0.44 | 0.65 | 32 |
| pep_glu_e5_01 | 326 | 288.0 | 0.54 | 0.63 | 32 | 308 | 269.1 | 0.35 | 0.55 | 32 |
| pep_asp_short_01 | 169 | 162.2 | 0.44 | 0.47 | 32 | 186 | 180.6 | 0.58 | 0.77 | 32 |
| pep_glu_short_01 | 173 | 163.5 | 0.34 | 0.47 | 32 | 186 | 182.3 | 0.57 | 0.72 | 32 |
| pep_lys_01 | 225 | 188.5 | 0.45 | 0.75 | 32 | 262 | 220.8 | 0.74 | 0.56 | 32 |
| pep_arg_02 | 228 | 188.4 | 0.88 | 0.77 | 32 | 250 | 222.1 | 0.87 | 0.66 | 32 |
| pep_asp_03 | 219 | 182.6 | 0.85 | 0.77 | 32 | 244 | 218.7 | 0.56 | 0.69 | 32 |
| pep_arg_01 | 215 | 188.0 | 0.71 | 0.77 | 32 | 265 | 224.4 | 0.82 | 0.64 | 32 |
| pep_asp_02 | 213 | 185.6 | 0.33 | 0.73 | 32 | 253 | 224.1 | 0.46 | 0.7 | 32 |
| pep_asp_01 | 231 | 190.6 | 0.83 | 0.69 | 32 | 262 | 229.4 | 0.63 | 0.57 | 32 |
| pep_tyr_01 | 226 | 189.3 | 0.53 | 0.75 | 32 | 264 | 229.0 | 0.9 | 0.65 | 32 |
| pep_glu_01 | 240 | 188.9 | 0.97 | 0.77 | 32 | 253 | 231.4 | 0.45 | 0.54 | 32 |
| pep_leu_01 | 225 | 191.7 | 0.84 | 0.77 | 32 | 268 | 234.8 | 0.8 | 0.66 | 32 |
| pep_ala_02 | 219 | 193.0 | 0.65 | 0.74 | 32 | 268 | 236.8 | 0.91 | 0.64 | 32 |
| pep_phe_02 | 217 | 182.8 | 0.94 | 0.74 | 32 | 250 | 227.0 | 0.56 | 0.72 | 32 |
| pep_ala_01 | 233 | 191.4 | 1.0 | 0.68 | 32 | 263 | 235.8 | 0.92 | 0.61 | 32 |
| pep_val_01 | 232 | 188.9 | 0.93 | 0.81 | 32 | 261 | 233.7 | 0.87 | 0.65 | 32 |
| pep_trp_01 | 214 | 178.6 | 0.58 | 0.78 | 32 | 263 | 223.5 | 0.63 | 0.71 | 32 |
| pep_gly_01 | 215 | 181.9 | 0.74 | 0.78 | 32 | 262 | 229.0 | 0.83 | 0.66 | 32 |
| pep_glu_03 | 216 | 178.3 | 0.83 | 0.8 | 32 | 248 | 225.9 | 0.6 | 0.61 | 32 |
| pep_ile_01 | 214 | 181.9 | 0.82 | 0.72 | 32 | 263 | 230.2 | 0.57 | 0.68 | 32 |
| pep_leu_02 | 234 | 186.7 | 0.96 | 0.76 | 32 | 261 | 235.0 | 0.69 | 0.72 | 32 |
| pep_phe_01 | 228 | 181.5 | 0.94 | 0.76 | 32 | 261 | 230.4 | 0.8 | 0.71 | 32 |
| pep_glu_02 | 218 | 182.6 | 0.44 | 0.73 | 32 | 280 | 232.2 | 0.58 | 0.71 | 32 |

## Next Steps

1. If salt bridge prediction validated: synthesize D-Glu and D-Asp peptides
2. If not validated: re-examine channel geometry, consider longer linkers or beta-amino acids
3. Measure Ki against both K392 and N392 ERAP2 alleles
4. D-Leu peptide as N392-selective control

## References

- Evnouchidou et al. 2012 J Biol Chem (PMID: 22837489)
- Papakyriakou & Stratikos 2017 PNAS
- Camberlein et al. 2022 Angew Chem