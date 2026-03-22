# Research Tools

Standalone scientific database query tools extracted from [ScienceClaw](https://github.com/lamm-mit/scienceclaw) (MIT, Apache 2.0).

## Tools

| Script | Database | Use Case | Dependency |
|--------|----------|----------|------------|
| `chembl_search.py` | ChEMBL (EBI) | Find existing ERAP2-targeting compounds, competitive landscape | `requests` |
| `pubmed_search.py` | PubMed (NCBI) | Prior art literature search for patent filing | `biopython` |
| `pubchem_search.py` | PubChem (NCBI) | Chemical properties, SMILES, cross-references | `requests` |
| `uniprot_fetch.py` | UniProt | ERAP2/ERAP1 protein info, sequences, annotations | `requests` |
| `blast_search.py` | NCBI BLAST | Sequence homology for designed binders | `biopython` |
| `foldseek_search.py` | PDB/AlphaFold DB | Structural homolog search for V2 divergent channel | `pandas`, foldseek binary |
| `tdc_predict.py` | TDC (HuggingFace) | ADMET predictions (BBB, hERG, CYP3A4) | `PyTDC`, `DeepPurpose`, `torch` |

## Install

```bash
# Core tools (ChEMBL, PubChem, UniProt, PubMed, BLAST)
pip install requests biopython

# Foldseek (structural search)
conda install -c conda-forge -c bioconda foldseek

# TDC (optional, heavy — only if predicting ADMET)
pip install PyTDC DeepPurpose torch dgl
```

## Usage — Python Import

```python
from tools import search_molecules, search_pubmed, fetch_articles, search_uniprot
from tools.foldseek_search import foldseek_search, check_design_novelty

# Find existing ERAP2 compounds
mols = search_molecules("ERAP2")

# Search literature for patent prior art
pmids = search_pubmed("ERAP2 selective binder peptide")
articles = fetch_articles(pmids)

# Check if designed binder has structural homologs
result = check_design_novelty("data/results/rfdiffusion_v2/binder_01.pdb", "pdb_db")
```

## Usage — CLI

```bash
python tools/chembl_search.py --query "ERAP2" --format detailed
python tools/pubmed_search.py --query "ERAP2 peptide binder" --year 2024 --format detailed
python tools/uniprot_fetch.py --search "ERAP2" --organism human --reviewed --format detailed
python tools/blast_search.py --query data/results/mpnn_v2/seq_01.fasta --program blastp
python tools/foldseek_search.py data/results/rfdiffusion_v2/binder_01.pdb pdb_db
```
