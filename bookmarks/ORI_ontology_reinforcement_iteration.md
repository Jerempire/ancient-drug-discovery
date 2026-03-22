# ORI: Functional Protein Design with Ontology Reinforcement Iteration

**Paper:** He, B., Qin, C., Zhao, Y. et al. *Functional protein design and enhancement with ontology reinforcement iteration.* Nature Communications (2026).
**DOI:** https://doi.org/10.1038/s41467-026-69855-6
**Authors:** Tencent AI for Life Sciences Lab
**Status:** Article in Press (unedited manuscript)
**Saved PDF:** `~/Downloads/Nature_article.pdf`

## What it does

Closed-loop protein engineering framework with 3 components:
- **PDA** (Protein Design Agent) — LLM translates goals into ontology prompts (function, species, structure, thermostability)
- **PGM** (3B-param generative model) — generates protein sequences conditioned on ontology terms
- **USM** (Unified Sequence Model) — predicts properties (expression, stability, function)
- **RLWF** (Reinforcement Learning from Wet-lab Feedback) — iterative improvement using experimental results

Key results: 100x lysozyme activity over natural baseline, chitinase stable at 85C, dual-function enzymes.

## Why bookmarked

- **Not relevant now** — ORI generates standalone functional proteins, not binders/complexes. Can't express selectivity targets or interface hotspots in its ontology vocabulary. Not open source.
- **Relevant later** — When optimizing a lead binder for developability (expression, stability, aggregation resistance), ORI's RLWF concept could be adapted. Also useful if thermostability engineering is needed for therapeutic candidates.

## When to revisit

- After synthesis/wet-lab validation of lead candidates
- If code/weights become publicly available
- When optimizing developability properties of confirmed binders
