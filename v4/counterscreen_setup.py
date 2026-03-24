"""Prepare ERAP1 and IRAP structures for DiffPepDock counter-screen.

Crops equivalent channel regions and creates docking config.
Run on Vast.ai after uploading structures.
"""
import os
import json
from Bio.PDB import PDBParser, PDBIO, Select

class ChannelSelect(Select):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def accept_residue(self, residue):
        resnum = residue.get_id()[1]
        return self.start <= resnum <= self.end

parser = PDBParser(QUIET=True)
os.makedirs("/workspace/docking_data", exist_ok=True)

# ERAP1: equivalent channel region 350-450 (same numbering as ERAP2, homologous fold)
# Residue 392 in ERAP1 = PRO (no charge)
s = parser.get_structure("erap1", "/workspace/data/erap1_wt_alphafold.pdb")
io = PDBIO()
io.set_structure(s)
io.save("/workspace/docking_data/erap1.pdb", ChannelSelect(350, 450))
count = sum(1 for r in s.get_residues() if 350 <= r.get_id()[1] <= 450 and r.get_id()[0] == ' ')
print(f"ERAP1 cropped: {count} residues (350-450), pos 392 = PRO")

# IRAP: crystal structure 5MJ6, equivalent channel region
# IRAP is larger (~1025 res in crystal) but the aminopeptidase domain aligns to ERAP2
# Equivalent channel region is approximately 350-450 (same fold)
# Residue 392 in IRAP = TYR
s2 = parser.get_structure("irap", "/workspace/data/irap_5MJ6_experimental.pdb")
# Check which chains are available
chains = [c.get_id() for c in s2.get_chains()]
print(f"IRAP chains: {chains}")

# Use first chain only
class ChainChannelSelect(Select):
    def __init__(self, chain, start, end):
        self.chain = chain
        self.start = start
        self.end = end
    def accept_chain(self, chain):
        return chain.get_id() == self.chain
    def accept_residue(self, residue):
        resnum = residue.get_id()[1]
        return self.start <= resnum <= self.end

io2 = PDBIO()
io2.set_structure(s2)
target_chain = chains[0]
io2.save("/workspace/docking_data/irap.pdb", ChainChannelSelect(target_chain, 350, 450))
count2 = sum(1 for r in s2[0][target_chain].get_residues() if 350 <= r.get_id()[1] <= 450 and r.get_id()[0] == ' ')
print(f"IRAP cropped: {count2} residues (350-450) chain {target_chain}, pos 392 = TYR")

# Create docking config — only test top 5 peptides (the ones that matter)
config = {
    "erap1": {
        "pdb": "erap1.pdb",
        "binding_site_residues": "388-410",
        "description": "ERAP1 substrate channel (counter-screen)"
    },
    "irap": {
        "pdb": "irap.pdb",
        "binding_site_residues": "388-410",
        "description": "IRAP substrate channel (counter-screen)"
    }
}
with open("/workspace/docking_data/docking_cases.json", "w") as f:
    json.dump(config, f, indent=2)

# Use reduced peptide library — just the synthesis candidates
fasta = """>pep_glu_long_01
EALVAAGLAGLA
>pep_asp_long_01
DALVAAGLAGLA
>pep_glu_e3_01
EAELAAGLAA
>pep_leu_01
LALVAAGLA
>pep_ala_01
AALVAAGLA
"""
with open("/workspace/docking_data/peptide_seq.fasta", "w") as f:
    f.write(fasta)

print("\nDocking config written: 2 targets x 5 peptides = 10 cases")
print("Ready for DiffPepDock preprocessing")
