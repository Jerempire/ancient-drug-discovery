# Ancient Drug Discovery Pipeline

AI-driven genomics-to-drug-discovery pipeline: maps genetic variants from ancient and modern disease datasets → protein structures → functional drug candidates.

## Pipeline Architecture

1. **Disease Genomics Data Collection** — GWAS Catalog, ClinVar, AADR (ancient DNA), UniProt, OpenTargets, DrugBank
2. **Variant → Disease Mapping** — ESM-2 (650M), V2P, Tranception, OpenTargets
3. **Protein Structure** — AlphaFold DB (pre-computed, 93.31 pLDDT for ERAP2), Boltz-2, ESMFold
4. **Drug/Binder Design** — **Proteina-Complexa** (NVIDIA, 68% hit rate), Dyno Psi-1, DiffDock, DrugGPT, LigandForge, RDKit
5. **Validation & Filtering** — Boltz-2 binding affinity, DeltaForge, OpenADMET, AutoDock Vina

## Proof of Concept: Black Death → Crohn's Disease

ERAP2 gene variants selected during the Black Death (1346-1353) are linked to modern Crohn's disease risk. We use the full pipeline to:
- Map ancient selective pressures at the ERAP2 locus
- Predict how variants alter protein function
- Design novel drug candidates targeting mutant ERAP2
- Validate computationally against known inhibitor DG013A

## Project Structure

```
scripts/
  01_fetch_gwas.py          — GWAS Catalog data (ERAP2 + Crohn's)
  02_fetch_ancient_variants.py — AADR ancient DNA + ClinVar + UniProt
  03_variant_analysis.py    — Cross-reference ancient/modern + OpenTargets
  04_drugbank_targets.py    — Query local DrugBank for ERAP2-related drugs
notebooks/
  04_esm_variant_effects.ipynb   — ESM Cambrian variant scoring (Colab)
  05_alphafold_structure.ipynb   — Structure prediction (Colab)
  06_diffdock_screening.ipynb    — Molecular docking (Colab)
  07_ligandforge_peptides.ipynb  — Peptide binder generation (Colab)
  08_druggpt_smallmol.ipynb      — Small molecule generation (Colab)
  09_validation.ipynb            — Binding affinity + ADMET filtering (Colab)
data/
  raw/       — Downloaded source data
  processed/ — Cleaned/merged outputs
```

## Setup

```bash
conda create -n ancient-drug-discovery python=3.11 rdkit biopython requests pandas -c conda-forge
conda activate ancient-drug-discovery
pip install transformers torch
```

## What Makes This Novel

- **Ancient-to-modern disease bridge** — systematic mapping of paleogenomic selective pressures to modern drug targets
- **Pharmacy domain knowledge** — pharmacological plausibility evaluation at every stage
- **Full open-source stack** — every tool is free for research
- **DrugBank integration** — 19,830 drugs, 24K targets, 2.9M interactions as reference

## Key References

- Klunk et al. 2022 Nature — Black Death survivors carry ERAP2 variants → Crohn's risk
- Conyngham approach — non-traditional researcher designed mRNA cancer vaccine with AI tools
- Zervoudi et al. 2013 PNAS — DG013A as ERAP2 inhibitor (validation reference)

## License

Research use only. DrugBank data subject to DrugBank license terms.
