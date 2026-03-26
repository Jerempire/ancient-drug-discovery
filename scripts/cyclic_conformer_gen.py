"""
Cyclic Peptide Conformer Generation for ERAP2 Inhibitor Panel

Generates 3D conformers for 6 constrained peptide constructs using RDKit.
For disulfide-bridged variants: builds linear peptide, adds SS constraint,
generates conformers, minimizes with MMFF94, filters by geometry.

Design Panel:
  Tier 1 (Lead):
    1. Relaxed Lasso:    R-C-G-G-W-C-F  (7aa, SS C2-C6)
    2. VAGSAF-Lasso:     V-A-C-S-A-C-F  (7aa, SS C3-C6)
    3. Hydroxyl-Bridge:  R-S-G-G-W-F    (6aa, linear control)
  Tier 2 (Exploratory):
    4. Lasso-Bolt:       R-G-C-W-C-F    (6aa, SS C3-C5)
    5. Bicyclic Cage:    C-R-C-G-C-W-C  (7aa, SS C1-C5 + C3-C7)
    6. Zn-Thiol:         C-G-G-W-A-F    (6aa, linear + free thiol)

Output: data/results/cyclic_peptides/conformers/<name>_conf<N>.pdb
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import os
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolTransforms
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.ERROR)

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(PROJECT, "data/results/cyclic_peptides/conformers")

# --- Peptide constructs ---
CONSTRUCTS = [
    {
        'name': 'relaxed_lasso',
        'tier': 1,
        'sequence': 'RCGGWCF',
        'length': 7,
        'topology': 'monocyclic',
        'disulfides': [(2, 6)],  # 1-indexed residue positions
        'rationale': 'Arg S1 bolt + Trp S3\' plug + relaxed 5-residue ring',
    },
    {
        'name': 'vagsaf_lasso',
        'tier': 1,
        'sequence': 'VACSACF',
        'length': 7,
        'topology': 'monocyclic',
        'disulfides': [(3, 6)],
        'rationale': 'Cyclized VAGSAF core, preserves existing binding contacts',
    },
    {
        'name': 'hydroxyl_bridge',
        'tier': 1,
        'sequence': 'RSGGWF',
        'length': 6,
        'topology': 'linear',
        'disulfides': [],
        'rationale': 'Linear control; Ser2 H-bonds Y398 for drift stabilization',
    },
    {
        'name': 'lasso_bolt',
        'tier': 2,
        'sequence': 'RGCWCF',
        'length': 6,
        'topology': 'monocyclic',
        'disulfides': [(3, 5)],
        'rationale': 'Tight lasso; N-terminal Arg tail, C-terminal ring',
    },
    {
        'name': 'bicyclic_cage',
        'tier': 2,
        'sequence': 'CRCGCWC',
        'length': 7,
        'topology': 'bicyclic',
        'disulfides': [(1, 5), (3, 7)],
        'rationale': 'Zero-drift moonshot; maximum rigidity, high clash risk',
    },
    {
        'name': 'zn_thiol',
        'tier': 2,
        'sequence': 'CGGWAF',
        'length': 6,
        'topology': 'linear',
        'disulfides': [],
        'rationale': 'N-terminal Cys thiol for direct Zn coordination; speculative',
    },
]

# One-letter to three-letter code
AA_3LETTER = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
}

# Amino acid SMILES (L-amino acids, backbone: N-CA-C(=O))
AA_SMILES = {
    'G': 'NCC(=O)O',
    'A': 'N[C@@H](C)C(=O)O',
    'V': 'N[C@@H](CC(C)C)C(=O)O',
    'S': 'N[C@@H](CO)C(=O)O',
    'F': 'N[C@@H](Cc1ccccc1)C(=O)O',
    'W': 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
    'R': 'N[C@@H](CCCNC(=N)N)C(=O)O',
    'C': 'N[C@@H](CS)C(=O)O',
}


def build_peptide_smiles(sequence):
    """Build a linear peptide SMILES string from amino acid sequence."""
    # For RDKit, we build the peptide as amide-linked residues
    # Using a simplified HELM-like approach
    parts = []
    for i, aa in enumerate(sequence):
        if aa == 'G':
            if i == 0:
                parts.append('NCC(=O)')
            else:
                parts.append('NCC(=O)')
        elif aa == 'A':
            parts.append('N[C@@H](C)C(=O)')
        elif aa == 'V':
            parts.append('N[C@@H](CC(C)C)C(=O)')
        elif aa == 'S':
            parts.append('N[C@@H](CO)C(=O)')
        elif aa == 'F':
            parts.append('N[C@@H](Cc1ccccc1)C(=O)')
        elif aa == 'W':
            parts.append('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)')
        elif aa == 'R':
            parts.append('N[C@@H](CCCNC(=N)N)C(=O)')
        elif aa == 'C':
            parts.append('N[C@@H](CS)C(=O)')
        else:
            raise ValueError(f"Unknown amino acid: {aa}")

    # Connect peptide bonds: the C(=O) of residue i bonds to N of residue i+1
    # We'll use RDKit's peptide builder instead
    return None  # We'll use a different approach


def build_peptide_mol(sequence, disulfides=None):
    """Build a peptide molecule using RDKit's fragment joining.

    Returns an RDKit Mol object with 3D coordinates.
    """
    # Build peptide SMILES manually
    # Format: H-[NH]-[CaH(R)]-[C(=O)]-[NH]-[CaH(R)]-[C(=O)]-...-OH
    fragments = []
    for i, aa in enumerate(sequence):
        side_chain = {
            'G': '[H]',
            'A': 'C',
            'V': 'CC(C)C',
            'S': 'CO',
            'F': 'Cc1ccccc1',
            'W': 'Cc1c[nH]c2ccccc12',
            'R': 'CCCNC(=N)N',
            'C': 'CS',
        }[aa]

        if aa == 'G':
            frag = 'NCC(=O)'
        else:
            frag = f'N[C@@H]({side_chain})C(=O)'
        fragments.append(frag)

    # Join with peptide bonds
    # Simple approach: build the full SMILES
    peptide_smiles = fragments[0]
    for frag in fragments[1:]:
        # Remove the terminal O from previous C(=O)O and add peptide bond
        peptide_smiles += frag

    # Add terminal OH
    peptide_smiles += 'O'

    # Fix: we need to properly link the amide bonds
    # Let's use a proper peptide SMILES notation
    # H2N-CH(R1)-CO-NH-CH(R2)-CO-...-NH-CH(Rn)-COOH
    smiles_parts = []
    for i, aa in enumerate(sequence):
        sc = {
            'G': '',  # No side chain for Gly (H only)
            'A': '(C)',
            'V': '(CC(C)C)',
            'S': '(CO)',
            'F': '(Cc1ccccc1)',
            'W': '(Cc1c[nH]c2ccccc12)',
            'R': '(CCCNC(=N)N)',
            'C': '(CS)',
        }[aa]

        if i == 0:
            if aa == 'G':
                smiles_parts.append('NCC(=O)')
            else:
                smiles_parts.append(f'N[C@@H]{sc}C(=O)')
        elif i == len(sequence) - 1:
            if aa == 'G':
                smiles_parts.append('NCC(=O)O')
            else:
                smiles_parts.append(f'N[C@@H]{sc}C(=O)O')
        else:
            if aa == 'G':
                smiles_parts.append('NCC(=O)')
            else:
                smiles_parts.append(f'N[C@@H]{sc}C(=O)')

    peptide_smiles = ''.join(smiles_parts)

    mol = Chem.MolFromSmiles(peptide_smiles)
    if mol is None:
        print(f"   WARNING: Could not parse SMILES: {peptide_smiles}")
        return None, peptide_smiles

    mol = Chem.AddHs(mol)
    return mol, peptide_smiles


def find_sulfur_atoms(mol, sequence, disulfide_pairs):
    """Find the sulfur atom indices for disulfide bond formation.

    Args:
        disulfide_pairs: list of (res_i, res_j) 1-indexed residue positions
    Returns:
        list of (S_atom_idx_i, S_atom_idx_j) pairs
    """
    if not disulfide_pairs:
        return []

    # Find all sulfur atoms in CYS residues
    sulfur_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            sulfur_indices.append(atom.GetIdx())

    # Map CYS positions in sequence to sulfur indices
    cys_positions = [i for i, aa in enumerate(sequence) if aa == 'C']

    pairs = []
    for res_i, res_j in disulfide_pairs:
        # Convert 1-indexed to 0-indexed position in CYS list
        cys_i_pos = [p for p in cys_positions if p == res_i - 1]
        cys_j_pos = [p for p in cys_positions if p == res_j - 1]

        if not cys_i_pos or not cys_j_pos:
            print(f"   WARNING: Cannot find CYS at positions {res_i}, {res_j}")
            continue

        # Map to sulfur index (in order of appearance)
        cys_idx_i = cys_positions.index(res_i - 1)
        cys_idx_j = cys_positions.index(res_j - 1)

        if cys_idx_i < len(sulfur_indices) and cys_idx_j < len(sulfur_indices):
            pairs.append((sulfur_indices[cys_idx_i], sulfur_indices[cys_idx_j]))

    return pairs


def add_disulfide_bonds(mol, sulfur_pairs):
    """Add disulfide bonds between specified sulfur atoms.

    Handles multiple SS bonds (bicyclic) by tracking atom indices through
    H-removal and bond-addition steps. Uses atom map numbers to survive
    index shifts when atoms are removed.
    """
    if not sulfur_pairs:
        return mol

    emol = Chem.RWMol(mol)

    # Tag all sulfur atoms with unique map numbers so we can track them
    sulfur_map = {}  # map_num -> original_idx
    map_counter = 100
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 16:
            atom.SetAtomMapNum(map_counter)
            sulfur_map[map_counter] = atom.GetIdx()
            map_counter += 1

    # Build pair list using map numbers
    # sulfur_pairs contains original atom indices; map them to map numbers
    idx_to_map = {atom.GetIdx(): atom.GetAtomMapNum()
                  for atom in emol.GetAtoms() if atom.GetAtomMapNum() > 0}
    pair_maps = []
    for s1_idx, s2_idx in sulfur_pairs:
        m1 = idx_to_map.get(s1_idx)
        m2 = idx_to_map.get(s2_idx)
        if m1 and m2:
            pair_maps.append((m1, m2))

    # Remove H from ALL sulfur atoms that participate in SS bonds first
    participating_maps = set()
    for m1, m2 in pair_maps:
        participating_maps.add(m1)
        participating_maps.add(m2)

    h_to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomMapNum() in participating_maps:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    h_to_remove.append(nbr.GetIdx())

    # Remove H atoms in reverse order to preserve lower indices
    for h_idx in sorted(set(h_to_remove), reverse=True):
        emol.RemoveAtom(h_idx)

    # Now add SS bonds using map numbers to find current indices
    for m1, m2 in pair_maps:
        s1_new = None
        s2_new = None
        for atom in emol.GetAtoms():
            if atom.GetAtomMapNum() == m1:
                s1_new = atom.GetIdx()
            elif atom.GetAtomMapNum() == m2:
                s2_new = atom.GetIdx()

        if s1_new is not None and s2_new is not None:
            # Check bond doesn't already exist
            existing_bond = emol.GetBondBetweenAtoms(s1_new, s2_new)
            if existing_bond is None:
                emol.AddBond(s1_new, s2_new, Chem.BondType.SINGLE)

    # Clear map numbers
    for atom in emol.GetAtoms():
        atom.SetAtomMapNum(0)

    try:
        mol = emol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"   WARNING: Sanitization failed after disulfide: {e}")
        return emol.GetMol()


def generate_conformers(mol, n_confs=100, random_seed=42):
    """Generate conformers using ETKDG v3."""
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0  # Use all cores
    params.pruneRmsThresh = 0.5  # Remove very similar conformers

    n_generated = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    return n_generated


def minimize_conformers(mol, max_iters=500):
    """Minimize all conformers with MMFF94 force field.

    Returns list of (conf_id, energy, converged).
    """
    results = []

    # Try MMFF94 first, fall back to UFF
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    use_uff = mmff_props is None

    for conf_id in range(mol.GetNumConformers()):
        try:
            if use_uff:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            else:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)

            if ff is None:
                results.append((conf_id, float('inf'), False))
                continue

            converged = ff.Minimize(maxIts=max_iters)
            energy = ff.CalcEnergy()
            results.append((conf_id, energy, converged == 0))
        except Exception as e:
            results.append((conf_id, float('inf'), False))

    return results


def check_ss_geometry(mol, conf_id, sulfur_pairs_indices):
    """Check disulfide bond geometry for a conformer.

    Ideal SS bond: 2.03 A, dihedral ~90 degrees.
    Returns dict with measurements.
    """
    if not sulfur_pairs_indices:
        return {'valid': True, 'note': 'no disulfide'}

    conf = mol.GetConformer(conf_id)
    results = {}

    for i, (s1, s2) in enumerate(sulfur_pairs_indices):
        pos1 = conf.GetAtomPosition(s1)
        pos2 = conf.GetAtomPosition(s2)

        dist = pos1.Distance(pos2)
        results[f'ss_{i+1}_distance'] = float(dist)
        results[f'ss_{i+1}_valid'] = 1.8 < dist < 2.3

    results['all_valid'] = all(v for k, v in results.items() if k.endswith('_valid'))
    return results


def compute_radius_of_gyration(mol, conf_id):
    """Compute radius of gyration for a conformer."""
    conf = mol.GetConformer(conf_id)
    positions = []
    masses = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])
        masses.append(atom.GetMass())

    positions = np.array(positions)
    masses = np.array(masses)
    total_mass = masses.sum()

    com = (positions.T * masses).T.sum(axis=0) / total_mass
    diffs = positions - com
    rg = np.sqrt((masses * np.sum(diffs**2, axis=1)).sum() / total_mass)
    return float(rg)


def save_conformer_pdb(mol, conf_id, path, name, sequence):
    """Save a conformer as PDB file."""
    try:
        pdb_block = Chem.MolToPDBBlock(mol, confId=conf_id)
        with open(path, 'w') as f:
            f.write(f"REMARK  Construct: {name}\n")
            f.write(f"REMARK  Sequence: {sequence}\n")
            f.write(pdb_block)
    except Exception as e:
        print(f"   WARNING: Could not save PDB for {name} conf {conf_id}: {e}")


def process_construct(construct):
    """Full pipeline for one construct: build, conformers, minimize, filter."""
    name = construct['name']
    seq = construct['sequence']
    ss_pairs = construct.get('disulfides', [])

    print(f"\n{'='*60}")
    print(f"Construct: {name} ({seq}, {construct['length']}aa, {construct['topology']})")
    print(f"Disulfides: {ss_pairs if ss_pairs else 'none'}")
    print(f"Rationale: {construct['rationale']}")
    print(f"{'='*60}")

    # Build molecule
    print("\n  Building molecule...")
    mol, smiles = build_peptide_mol(seq, ss_pairs)
    if mol is None:
        return {'name': name, 'status': 'FAILED', 'error': f'SMILES parse failed: {smiles}'}

    print(f"  SMILES: {smiles}")
    print(f"  Heavy atoms: {mol.GetNumHeavyAtoms()}, Total atoms: {mol.GetNumAtoms()}")

    # Add disulfide bonds if needed
    ss_atom_pairs = []
    if ss_pairs:
        print("  Adding disulfide constraints...")
        ss_atom_pairs = find_sulfur_atoms(mol, seq, ss_pairs)
        if ss_atom_pairs:
            mol = add_disulfide_bonds(mol, ss_atom_pairs)
            print(f"  Disulfide bonds added: {len(ss_atom_pairs)}")
            # Re-find sulfur pairs after modification
            ss_atom_pairs = find_sulfur_atoms(mol, seq, ss_pairs)

    # Generate conformers
    print("  Generating conformers (ETKDG v3)...")
    n_gen = generate_conformers(mol, n_confs=100)
    print(f"  Generated: {n_gen} conformers")

    if n_gen == 0:
        print("  FAILED: No conformers generated")
        return {'name': name, 'status': 'FAILED', 'error': 'No conformers generated'}

    # Minimize
    print("  Minimizing with MMFF94...")
    min_results = minimize_conformers(mol)
    converged = [r for r in min_results if r[2]]
    print(f"  Converged: {len(converged)}/{len(min_results)}")

    if not converged:
        # Use all results even if not converged
        converged = min_results

    # Sort by energy
    converged.sort(key=lambda x: x[1])

    # Check SS geometry
    print("  Checking geometry...")
    valid_confs = []
    for conf_id, energy, conv in converged:
        ss_geom = check_ss_geometry(mol, conf_id, ss_atom_pairs)
        rg = compute_radius_of_gyration(mol, conf_id)

        valid_confs.append({
            'conf_id': conf_id,
            'energy': energy,
            'converged': conv,
            'ss_geometry': ss_geom,
            'radius_of_gyration': rg,
            'ss_valid': ss_geom.get('all_valid', True),
        })

    # Filter: valid SS geometry, lowest energy
    if ss_pairs:
        ss_valid = [c for c in valid_confs if c['ss_valid']]
        print(f"  SS-valid conformers: {len(ss_valid)}/{len(valid_confs)}")
        if ss_valid:
            valid_confs = ss_valid

    # Take top 5
    top_confs = valid_confs[:5]

    # Save PDBs
    construct_dir = os.path.join(OUTPUT_DIR, name)
    os.makedirs(construct_dir, exist_ok=True)

    print(f"\n  Top 5 conformers:")
    print(f"  {'Rank':>5} {'ConfID':>7} {'Energy':>12} {'Rg':>8} {'SS_Valid':>8}")
    print(f"  {'-'*45}")

    saved_files = []
    for rank, conf in enumerate(top_confs):
        print(f"  {rank+1:>5} {conf['conf_id']:>7} {conf['energy']:>12.2f} "
              f"{conf['radius_of_gyration']:>8.2f} {str(conf['ss_valid']):>8}")

        pdb_path = os.path.join(construct_dir, f"{name}_conf{rank+1}.pdb")
        save_conformer_pdb(mol, conf['conf_id'], pdb_path, name, seq)
        saved_files.append(pdb_path)

    # Summary
    result = {
        'name': name,
        'tier': construct['tier'],
        'sequence': seq,
        'length': construct['length'],
        'topology': construct['topology'],
        'disulfides': ss_pairs,
        'status': 'SUCCESS',
        'smiles': smiles,
        'n_heavy_atoms': mol.GetNumHeavyAtoms(),
        'n_conformers_generated': n_gen,
        'n_converged': len([r for r in min_results if r[2]]),
        'n_ss_valid': len([c for c in valid_confs if c['ss_valid']]),
        'top_conformers': top_confs,
        'saved_files': saved_files,
        'best_energy': top_confs[0]['energy'] if top_confs else None,
        'best_rg': top_confs[0]['radius_of_gyration'] if top_confs else None,
    }

    # Save per-construct JSON
    json_path = os.path.join(construct_dir, f"{name}_summary.json")
    with open(json_path, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\n  Summary saved: {json_path}")

    return result


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 80)
    print("CYCLIC PEPTIDE CONFORMER GENERATION")
    print("ERAP2 Inhibitor Design Panel — 6 Constructs")
    print("=" * 80)

    all_results = []
    for construct in CONSTRUCTS:
        result = process_construct(construct)
        all_results.append(result)

    # Summary table
    print("\n" + "=" * 80)
    print("PANEL SUMMARY")
    print("=" * 80)
    print(f"\n{'Name':<20} {'Tier':>4} {'Seq':<10} {'Status':<8} {'Confs':>6} "
          f"{'SS_OK':>6} {'BestE':>10} {'Rg':>8}")
    print("-" * 78)

    for r in all_results:
        if r['status'] == 'SUCCESS':
            print(f"{r['name']:<20} {r['tier']:>4} {r['sequence']:<10} "
                  f"{r['status']:<8} {r['n_conformers_generated']:>6} "
                  f"{r['n_ss_valid']:>6} {r['best_energy']:>10.1f} "
                  f"{r['best_rg']:>8.2f}")
        else:
            print(f"{r['name']:<20} {r.get('tier','?'):>4} {r.get('sequence','?'):<10} "
                  f"{r['status']:<8} {'--':>6} {'--':>6} {'--':>10} {'--':>8}")

    # Notes on limitations
    print("\n" + "-" * 80)
    print("LIMITATIONS & NEXT STEPS")
    print("-" * 80)
    print("""
  1. These are GAS-PHASE conformers, not docked into ERAP2.
     Next step: Boltz-2 docking (scripts/cyclic_boltz2_yamls/) or Rosetta CyclicPeptideDocking.

  2. Boltz-2 treats all peptides as LINEAR — it cannot model cyclic topology natively.
     Workaround: dock linear version, then check if Cys residues are close enough for ring closure.
     Alternative: Rosetta CyclicPeptideDocking (better for constrained peptides).

  3. Bicyclic cage (#5) has very tight geometric tolerance (~0.5A).
     If no SS-valid conformers, the topology is likely incompatible.

  4. Radius of gyration (Rg) indicates compactness.
     ERAP2 channel diameter is ~8-10A. Rg > 6A suggests the peptide is too extended.

  5. SAVE CHECKPOINT: All conformers saved to data/results/cyclic_peptides/conformers/
""")

    # Save master summary
    master_json = os.path.join(OUTPUT_DIR, "panel_summary.json")
    with open(master_json, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"  Master summary: {master_json}")


if __name__ == '__main__':
    main()
