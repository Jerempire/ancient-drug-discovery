# V4 Cleavage-Resistant Peptide Screen

## Purpose
Test whether modified (non-cleavable) versions of V4 short-stack peptides still bind ERAP2's channel.
Original V4 peptides are selective substrates — ERAP2 eats them. These modifications should make them
resistant to hydrolysis while (hopefully) preserving binding.

## Modifications tested

| Prefix | Modification | Mechanism | Base peptide |
|--------|-------------|-----------|-------------|
| dVAGSAF | D-Valine at P1 | Wrong chirality — zinc can't orient for cleavage | VAGSAF |
| dIAFSAF | D-Isoleucine at P1 | Same, on best K392-selective peptide | IAFSAF |
| NMeVAGSAF | N-methyl-Valine at P1 | Blocks backbone NH needed for catalysis | VAGSAF |
| riVAGSAF | Retro-inverso (all-D, reversed) | Mimics L-peptide shape but fully D | FASGAV |

## CCD codes used
- DVA: D-Valine
- DIL: D-Isoleucine
- MVA: N-methyl-Valine
- DPN: D-Phenylalanine
- DAL: D-Alanine
- DSN: D-Serine

## Targets (same 100aa truncations as V4)
- erap2k392: ERAP2 K392 wild-type (want HIGH ipTM — must still bind)
- erap1: ERAP1 (want LOW ipTM — must stay selective)
- irap: IRAP (want LOW ipTM — must stay selective)

## Success criteria
- K392 ipTM stays above 0.7 (original VAGSAF = 0.905)
- ERAP1 and IRAP ipTM stay below 0.5
- Any modification maintaining selectivity + binding = candidate for synthesis

## Run command
```bash
boltz predict boltz2_inputs/ --diffusion_samples 3 --seed 42 --out_dir ../boltz2_results/
```

## Baselines (original L-peptides, from V43 validation)
| Peptide | K392 ipTM | ERAP1 ipTM | IRAP ipTM |
|---------|-----------|------------|-----------|
| VAGSAF  | 0.905     | 0.335      | 0.236     |
| IAFSAF  | 0.870     | —          | —         |
