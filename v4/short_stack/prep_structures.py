"""Prepare all 4 receptor structures for short-stack docking."""
import os
from Bio.PDB import PDBParser, PDBIO, Select

class ChannelSelect(Select):
    def __init__(self, start, end):
        self.start, self.end = start, end
    def accept_residue(self, residue):
        return self.start <= residue.get_id()[1] <= self.end

class ChainChannelSelect(Select):
    def __init__(self, chain, start, end):
        self.chain, self.start, self.end = chain, start, end
    def accept_chain(self, chain):
        return chain.get_id() == self.chain
    def accept_residue(self, residue):
        return self.start <= residue.get_id()[1] <= self.end

parser = PDBParser(QUIET=True)
io = PDBIO()
os.makedirs("/workspace/docking_data", exist_ok=True)

# ERAP2 K392 (wildtype)
s = parser.get_structure("e2", "/workspace/data/erap2_wt_alphafold.pdb")
io.set_structure(s)
io.save("/workspace/docking_data/erap2_k392.pdb", ChannelSelect(350, 450))
print("ERAP2-K392 cropped")

# ERAP2 N392
with open("/workspace/docking_data/erap2_k392.pdb") as f:
    lines = f.readlines()
n392_lines = []
for line in lines:
    if line.startswith("ATOM") and line[22:26].strip() == "392":
        atom_name = line[12:16].strip()
        if atom_name in ("N", "CA", "C", "O", "CB"):
            n392_lines.append(line[:17] + "ASN" + line[20:])
        elif atom_name == "CG":
            n392_lines.append(line[:17] + "ASN" + line[20:])
        elif atom_name == "CD":
            n392_lines.append(line[:12] + " OD1" + line[16:17] + "ASN" + line[20:])
        elif atom_name == "CE":
            n392_lines.append(line[:12] + " ND2" + line[16:17] + "ASN" + line[20:])
    else:
        n392_lines.append(line)
with open("/workspace/docking_data/erap2_n392.pdb", "w") as f:
    f.writelines(n392_lines)
print("ERAP2-N392 cropped")

# ERAP1
s1 = parser.get_structure("e1", "/workspace/data/erap1_wt_alphafold.pdb")
io.set_structure(s1)
io.save("/workspace/docking_data/erap1.pdb", ChannelSelect(350, 450))
print("ERAP1 cropped")

# IRAP
s2 = parser.get_structure("ir", "/workspace/data/irap_5MJ6_experimental.pdb")
io.set_structure(s2)
io.save("/workspace/docking_data/irap.pdb", ChainChannelSelect("A", 350, 450))
print("IRAP cropped")

print("\n4 structures ready: erap2_k392, erap2_n392, erap1, irap")
