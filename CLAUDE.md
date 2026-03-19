# CLAUDE.md - Ancient Drug Discovery Pipeline

## Project Overview
AI-driven pipeline mining evolutionary immune selection for undrugged cancer/autoimmune therapeutic targets. Maps ancient pathogen selection → immune gene variants → cancer immune evasion mechanisms → protein structures → drug candidates.

**Core thesis**: Genes shaped by 10,000 years of plague, malaria, TB, and smallpox are the same genes cancers hijack to evade the immune system. See FRAMEWORK.md for full thesis and evidence.

## Pipeline Stages
1. Data collection (GWAS, ClinVar, UniProt, OpenTargets, DrugBank) — local scripts 01-06
2. Cancer evidence compilation (OpenTargets cancer, ClinicalTrials.gov) — local script 07
3. Variant analysis + ESM-2 scoring — local script 05
4. Protein structure prediction (AlphaFold3/Boltz-2) — Colab
5. Drug design (RFdiffusion, DiffDock, DrugGPT) — Colab
6. Validation (Boltz-2 affinity, RDKit ADMET) — Colab + local

## Multi-Target Architecture
All targets defined in `targets.yaml`. Three immune defense steps:

| Step | Gene | Pathogen | Cancer Link |
|---|---|---|---|
| 1. Detection | **ERAP2** (primary) | Black Death | Antigen presentation / immune evasion |
| 1. Detection | **NOD2** | Leprosy | Crohn's → colorectal cancer |
| 2. Activation | **CCR5** | Smallpox | Metastasis (breast, colorectal) |
| 3. Killing | **G6PD** | Malaria | Tumor metabolic vulnerability |
| 3. Killing | **SLC11A1** | TB | Macrophage M1/M2 polarization |
| 3. Killing | **DARC** | Malaria | Prostate cancer disparity |

## Key Files
- `targets.yaml` — multi-target config (variants, compounds, references)
- `FRAMEWORK.md` — unified thesis + historical timeline + evidence
- `NARRATIVE.md` — ERAP2 proof-of-concept narrative
- `scripts/01-05` — ERAP2-specific pipeline (original)
- `scripts/06_multi_target_fetch.py` — generalized data fetcher for all targets
- `scripts/07_cancer_connections.py` — cancer evidence matrix + clinical trials

## Environment
- Conda: `ancient-drug-discovery` (Python 3.11, rdkit, biopython, pandas, pyyaml)
- Colab: notebooks/ dir (synced to Google Drive for GPU inference)
- DrugBank: referenced from `pharmacy-ai-project/data/drugbank/drugbank.db`

## Rules
- Scripts 01-07 run locally (CPU, data collection/querying)
- Notebooks 04-09 run on Colab (GPU inference)
- Never commit data/ contents (large files, .gitignored)
- Cross-validate structure predictions: AlphaFold vs Boltz-2 vs experimental PDB
- Per-target data goes in `data/processed/<gene_symbol>/`
- Patent: one provisional per target class, $75-150 micro entity
