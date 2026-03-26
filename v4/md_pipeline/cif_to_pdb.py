"""Convert Boltz-2 CIF files to PDB for OpenMM.

Usage: python cif_to_pdb.py <input.cif> <output.pdb>

Handles Boltz-2 quirks:
- Multi-chain complexes (peptide=A, receptor=B)
- Non-standard residue naming
"""
import sys
import os

def cif_to_pdb(cif_path, pdb_path):
    try:
        import gemmi
    except ImportError:
        sys.exit("ERROR: pip install gemmi")

    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)

    if len(st) == 0:
        sys.exit(f"ERROR: No models in {cif_path}")

    st.setup_entities()
    st.assign_label_seq_id()

    st.write_pdb(pdb_path)
    model = st[0]
    chains = [ch.name for ch in model]
    n_res = sum(len(ch) for ch in model)
    print(f"  Converted: {os.path.basename(cif_path)} -> {os.path.basename(pdb_path)}")
    print(f"    Chains: {chains}, Residues: {n_res}")
    return pdb_path


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python cif_to_pdb.py <input.cif> <output.pdb>")
        sys.exit(1)
    cif_to_pdb(sys.argv[1], sys.argv[2])
