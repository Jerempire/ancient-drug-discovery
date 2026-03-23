# V4/PepMLM — Parallel Branch (Extended PepMLM Exploration)

This folder was created in a **parallel Claude Code session** running alongside the main session that produced `scripts/pepmlm_generate.py` and `data/results/pepmlm_hybrids/`.

## What this branch did

**Extended PepMLM generation** — ran PepMLM-650M at larger scale (1,595 unique peptides vs 48 in the main branch) to stress-test whether the model would ever generate E/D at P1.

### Key finding: Zero E/D at P1 across 1,595 peptides

Two context windows (narrow res 388-410, wide res 350-430), four peptide lengths (10, 11, 12, 13), 400 peptides per context. The model overwhelmingly generates K (42-47%), S (22-36%), and N/I/V. Not a single Glu or Asp at P1.

This is the **novelty evidence** for the patent — a state-of-the-art protein language model trained on 150M+ sequences does not converge on the salt bridge design strategy.

## Files

```
v4/pepmlm/
├── generate_peptides.py          # PepMLM generation + cross-filter + FASTA output
├── run_pepmlm_vastai.sh          # Vast.ai runner (narrow + wide context)
├── analyze_pepmlm_results.py     # Comparison with hand-designed V4 library
├── boltz2_3hits.sh               # Boltz-2 validation of 3 raw K392-selective hits
├── results/
│   ├── pepmlm_k392_candidates.csv      # 798 peptides scored against K392
│   ├── pepmlm_n392_candidates.csv      # 798 peptides scored against N392
│   ├── pepmlm_selective_narrow.csv     # Cross-filtered (narrow context)
│   ├── pepmlm_selective_wide.csv       # Cross-filtered (wide context)
│   ├── pepmlm_boltz2_narrow.fasta      # 3 K392-selective from narrow
│   ├── pepmlm_boltz2_wide.fasta        # Top 30 by K392 PPL from wide
│   ├── pepmlm_boltz2_merged.fasta      # 33 merged candidates
│   └── boltz2_inputs/                  # 66 Boltz-2 YAMLs (33 peptides x 2 variants)
└── README.md                           # This file
```

## Relationship to main branch

```
Main session (scripts/ + data/results/):
  pepmlm_generate.py → 48 scaffolds → pepmlm_saltbridge_hybrids.py → 12 E/D hybrids
  → Boltz-2 on Vast.ai → 7/12 K392-selective (the actionable result)

This session (v4/pepmlm/):
  generate_peptides.py → 1,595 scaffolds → cross-filter → 3 weak K392-selective
  → Confirmed: PepMLM alone cannot find salt bridge (patent novelty evidence)
```

Both sessions ran on the same day (2026-03-23). The main session's hybrid approach (grafting E/D onto PepMLM scaffolds) was the productive path. This session's extended search confirms the novelty claim.

## Compute cost

- PepMLM generation + scoring: ~$0.05 (RTX 4090, ~8 min)
- Boltz-2 attempt on 3 hits: ~$0.15 (incomplete — cuequivariance issues on torch 2.5.1)
