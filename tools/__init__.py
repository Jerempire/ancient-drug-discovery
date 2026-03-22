"""
Research tools extracted from ScienceClaw (MIT, Apache 2.0).
https://github.com/lamm-mit/scienceclaw

Each module is standalone — no framework dependencies.
Use as importable functions or CLI scripts.

Dependencies:
    pip install requests biopython

Optional (for TDC predictions):
    pip install PyTDC DeepPurpose torch dgl

Foldseek requires binary install:
    conda install -c conda-forge -c bioconda foldseek
"""

# ChEMBL — drug/compound lookup (needs: requests)
from .chembl_search import search_molecules, get_molecule

# PubMed — literature search (needs: biopython)
from .pubmed_search import search_pubmed, fetch_articles

# PubChem — chemical properties (needs: requests)
from .pubchem_search import search_compound, get_compound_properties

# UniProt — protein info & sequences (needs: requests)
from .uniprot_fetch import search_uniprot, fetch_protein, fetch_fasta

# BLAST — sequence homology search (needs: biopython)
from .blast_search import run_blast

# Foldseek — see SKILL.md in scienceclaw/skills/foldseek/ for CLI usage
# (binary tool, not a Python import — use subprocess or the function in foldseek SKILL.md)
