"""Find ERAP2 surface residues that diverge from ERAP1 — selectivity targets."""
from Bio.PDB import PDBParser, NeighborSearch
from Bio import pairwise2
import json

parser = PDBParser(QUIET=True)

e2 = parser.get_structure("erap2", "/workspace/data/structures/erap2_wt_alphafold.pdb")
e1 = parser.get_structure("erap1", "/workspace/data/structures/erap1_wt_alphafold.pdb")

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_seq_and_residues(struct):
    seq, residues = [], []
    for r in struct.get_residues():
        if r.get_resname() in AA3TO1:
            seq.append(AA3TO1[r.get_resname()])
            residues.append(r)
    return "".join(seq), residues


e2_seq, e2_res = get_seq_and_residues(e2)
e1_seq, e1_res = get_seq_and_residues(e1)
print("ERAP2: %d residues" % len(e2_seq))
print("ERAP1: %d residues" % len(e1_seq))

# Align
alignments = pairwise2.align.globalxx(e2_seq, e1_seq, one_alignment_only=True)
e2_aln, e1_aln = alignments[0].seqA, alignments[0].seqB

# Find divergent positions
divergent = []
e2_idx = 0
for i in range(len(e2_aln)):
    if e2_aln[i] != "-":
        if e1_aln[i] == "-" or e2_aln[i] != e1_aln[i]:
            if e2_idx < len(e2_res):
                divergent.append(e2_res[e2_idx].get_id()[1])
        e2_idx += 1

pct = 100 * len(divergent) / len(e2_seq)
print("\nDivergent residues: %d / %d (%.0f%%)" % (len(divergent), len(e2_seq), pct))

# Surface-exposed divergent residues
all_atoms = list(e2.get_atoms())
ns = NeighborSearch(all_atoms)
resnum_to_idx = {r.get_id()[1]: i for i, r in enumerate(e2_res)}

surface_divergent = []
for resnum in divergent:
    idx = resnum_to_idx.get(resnum)
    if idx is None:
        continue
    r = e2_res[idx]
    try:
        ca = r["CA"]
    except KeyError:
        continue
    neighbors = ns.search(ca.get_vector().get_array(), 8.0)
    n_neighbor_res = len(set(n.get_parent().get_id()[1] for n in neighbors)) - 1
    if n_neighbor_res < 18:
        surface_divergent.append(resnum)

print("Surface-exposed divergent: %d" % len(surface_divergent))

# Cluster nearby divergent surface residues
clusters = []
used = set()
for r1 in surface_divergent:
    if r1 in used:
        continue
    cluster = [r1]
    used.add(r1)
    idx1 = resnum_to_idx[r1]
    try:
        ca1 = e2_res[idx1]["CA"]
    except KeyError:
        continue
    for r2 in surface_divergent:
        if r2 in used:
            continue
        idx2 = resnum_to_idx.get(r2)
        if idx2 is None:
            continue
        try:
            ca2 = e2_res[idx2]["CA"]
        except KeyError:
            continue
        dist = ca1 - ca2
        if dist < 15.0:
            cluster.append(r2)
            used.add(r2)
    if len(cluster) >= 4:
        clusters.append(sorted(cluster))

clusters.sort(key=len, reverse=True)

print("\nDivergent surface clusters (>=4 residues within 15A):")
for i, c in enumerate(clusters[:8]):
    rng = "%d-%d" % (min(c), max(c))
    print("  Cluster %d: %d residues, range %s" % (i + 1, len(c), rng))
    aa_list = []
    for resnum in c:
        idx = resnum_to_idx[resnum]
        aa_list.append("%s%d" % (e2_seq[idx], resnum))
    print("    Residues: %s" % c)
    print("    AAs: %s" % ", ".join(aa_list))

# Save for use by RFdiffusion
output = {
    "divergent_count": len(divergent),
    "surface_divergent_count": len(surface_divergent),
    "surface_divergent": surface_divergent,
    "clusters": [{"residues": c, "range": "%d-%d" % (min(c), max(c))} for c in clusters[:8]],
}
with open("/workspace/results/erap2_divergent_surface.json", "w") as f:
    json.dump(output, f, indent=2)
print("\nSaved: /workspace/results/erap2_divergent_surface.json")
