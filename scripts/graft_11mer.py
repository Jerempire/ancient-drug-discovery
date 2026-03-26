"""Graft 11-mer AILKLYSSKKY from Boltz-2 channel prediction into full ERAP2."""
import sys
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer

pdb_parser = PDBParser(QUIET=True)
cif_parser = MMCIFParser(QUIET=True)

# Load full ERAP2
full_struct = pdb_parser.get_structure("erap2", "/workspace/structures/erap2_wt_alphafold.pdb")
full_model = full_struct[0]
full_chain_A = full_model["A"]
print("Full ERAP2: %d residues" % len(list(full_chain_A.get_residues())))

# Load Boltz-2 prediction (best model = model_1, ipTM 0.892)
cif_path = "/workspace/boltz_11mer/boltz_results_n392s_045_11aa_K392/predictions/n392s_045_11aa_K392/n392s_045_11aa_K392_model_1.cif"
boltz_struct = cif_parser.get_structure("boltz", cif_path)
boltz_model = boltz_struct[0]
boltz_chains = list(boltz_model.get_chains())
boltz_channel = boltz_chains[0]
boltz_peptide = boltz_chains[1]
print("Boltz channel: %d residues" % len(list(boltz_channel.get_residues())))
print("Boltz peptide: %d residues" % len(list(boltz_peptide.get_residues())))

# Channel = residues 350-450 of full ERAP2 (1-indexed)
channel_start = 350

boltz_cas = []
full_cas = []
boltz_channel_res = sorted(boltz_channel.get_residues(), key=lambda r: r.id[1])

for i, res in enumerate(boltz_channel_res):
    if "CA" not in res:
        continue
    full_resid = channel_start + i
    try:
        full_res = full_chain_A[(" ", full_resid, " ")]
        if "CA" in full_res:
            boltz_cas.append(res["CA"])
            full_cas.append(full_res["CA"])
    except KeyError:
        pass

print("Aligned %d CA atoms" % len(boltz_cas))

if len(boltz_cas) < 10:
    print("ERROR: Too few alignment atoms!")
    sys.exit(1)

# Superimpose Boltz onto full protein coordinate frame
sup = Superimposer()
sup.set_atoms(full_cas, boltz_cas)
print("Alignment RMSD: %.2f A" % sup.rms)

# Apply transformation to all Boltz atoms
sup.apply(list(boltz_model.get_atoms()))

# Add peptide chain to full structure
boltz_peptide.id = "B"
# Detach from boltz model first
boltz_peptide.detach_parent()
full_model.add(boltz_peptide)

# Save
io = PDBIO()
io.set_structure(full_struct)
outpath = "/workspace/structures/AILKLYSSKKY_vs_erap2k392_full.pdb"
io.save(outpath)

# Verify
check = pdb_parser.get_structure("check", outpath)
for chain in check[0]:
    n_res = len([r for r in chain if r.id[0] == " "])
    print("Output chain %s: %d residues" % (chain.id, n_res))

print("\nSaved: %s" % outpath)
