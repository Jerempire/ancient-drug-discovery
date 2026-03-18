# CLAUDE.md - Ancient Drug Discovery Pipeline

## Project Overview
AI-driven genomics-to-drug-discovery pipeline. Maps ancient + modern genetic variants → disease mechanisms → protein structures → drug candidates. PoC: Black Death ERAP2 variants → Crohn's disease → novel therapeutics.

## Pipeline Stages
1. Data collection (GWAS, ClinVar, AADR, UniProt, OpenTargets, DrugBank) — local scripts
2. Variant → disease mapping (V2P, ESM Cambrian) — Colab
3. Protein structure prediction (AlphaFold3/OpenFold3, Boltz-2) — Colab
4. Drug design (LigandForge, RFdiffusion, DiffDock, DrugGPT) — Colab
5. Validation (Boltz-2 affinity, DeltaForge, RDKit ADMET) — Colab + local

## Key Targets
- **Gene:** ERAP2 (ENSG00000164308)
- **Protein:** Q6P179 (960 aa, aminopeptidase)
- **Primary variant:** rs2549794 (Black Death → Crohn's link)
- **Known inhibitor:** DG013A (validation reference)
- **Disease:** Crohn's disease (EFO_0000384)

## Environment
- Conda: `ancient-drug-discovery` (Python 3.11, rdkit, biopython, pandas)
- Colab: notebooks/ dir (synced to Google Drive for GPU inference)
- DrugBank: referenced from `pharmacy-ai-project/data/drugbank/drugbank.db`

## Rules
- Scripts 01-04 run locally (CPU, data collection/querying)
- Notebooks 04-09 run on Colab (GPU inference)
- Never commit data/ contents (large files, .gitignored)
- Cross-validate structure predictions: AlphaFold vs Boltz-2 vs experimental PDB
