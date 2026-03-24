"""
Structural alignment of ERAP2 (7SH0), ERAP1 (6RQX), IRAP (4PJ6)
to find spatially equivalent residues for ERAP2 K392 and nearby S1' pocket residues.

Uses BioPython PDB module for:
1. Fetch PDB structures
2. Structural superposition (CE-like via Superimposer on matched CAs)
3. Nearest-residue search by CA-CA distance after alignment
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import os
import warnings
warnings.filterwarnings('ignore')

from Bio.PDB import PDBParser, PDBList, Superimposer, Selection
from Bio.PDB.Polypeptide import is_aa
import numpy as np

# ── Config ──────────────────────────────────────────────────────────────
ERAP2_PDB = "7SH0"
ERAP1_PDB = "6RQX"
IRAP_PDB  = "4PJ6"

# ERAP2 residues of interest (chain A assumed)
ERAP2_RESIDUES = [392, 398, 403, 406, 412, 414]

WORK_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(WORK_DIR, "pdb_cache")
os.makedirs(PDB_DIR, exist_ok=True)

# ── Fetch PDBs ──────────────────────────────────────────────────────────
pdbl = PDBList()
parser = PDBParser(QUIET=True)

def fetch_and_parse(pdb_id):
    fn = pdbl.retrieve_pdb_file(pdb_id, pdir=PDB_DIR, file_format="pdb")
    structure = parser.get_structure(pdb_id, fn)
    return structure

print("Fetching structures...")
erap2_struct = fetch_and_parse(ERAP2_PDB)
erap1_struct = fetch_and_parse(ERAP1_PDB)
irap_struct  = fetch_and_parse(IRAP_PDB)

# ── Helper: get CA atoms for a chain ────────────────────────────────────
def get_ca_residues(structure, chain_id=None):
    """Return list of (residue, CA_atom) for standard amino acids."""
    model = structure[0]
    results = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for res in chain:
            if is_aa(res, standard=True) and 'CA' in res:
                results.append((res, res['CA']))
    return results

# ── Determine chains ────────────────────────────────────────────────────
def list_chains(structure):
    model = structure[0]
    chains = []
    for chain in model:
        n_aa = sum(1 for r in chain if is_aa(r, standard=True))
        if n_aa > 50:  # skip small peptides/ligands
            chains.append((chain.id, n_aa))
    return chains

print("\nChains in each structure:")
for name, struct in [("ERAP2 (7SH0)", erap2_struct), ("ERAP1 (6RQX)", erap1_struct), ("IRAP (4PJ6)", irap_struct)]:
    chains = list_chains(struct)
    print(f"  {name}: {chains}")

# We'll use chain A for all (or the longest chain)
ERAP2_CHAIN = "A"
ERAP1_CHAIN = "A"
IRAP_CHAIN = "A"

# ── Structural alignment via sequence-independent CA matching ───────────
# BioPython's Superimposer needs paired atoms. We'll use a two-step approach:
# 1. Get all CAs from both structures
# 2. Do iterative closest-point alignment (or use sequence alignment to pair CAs)
# For M1 aminopeptidases with ~50% identity, sequence-based pairing works well.

from Bio import pairwise2
from Bio.Data.IUPACData import protein_letters_3to1

def three_to_one(resname):
    """Convert 3-letter AA code to 1-letter."""
    key = resname[0] + resname[1:].lower()  # e.g. ALA -> Ala
    return protein_letters_3to1[key]

def get_sequence_and_cas(structure, chain_id):
    """Get sequence string and corresponding CA atoms."""
    model = structure[0]
    chain = model[chain_id]
    seq = []
    cas = []
    resis = []
    for res in chain:
        if is_aa(res, standard=True) and 'CA' in res:
            try:
                aa = three_to_one(res.get_resname())
            except KeyError:
                continue
            seq.append(aa)
            cas.append(res['CA'])
            resis.append(res)
    return ''.join(seq), cas, resis

print("\nExtracting sequences and CAs...")
erap2_seq, erap2_cas, erap2_resis = get_sequence_and_cas(erap2_struct, ERAP2_CHAIN)
erap1_seq, erap1_cas, erap1_resis = get_sequence_and_cas(erap1_struct, ERAP1_CHAIN)
irap_seq, irap_cas, irap_resis    = get_sequence_and_cas(irap_struct, IRAP_CHAIN)

print(f"  ERAP2: {len(erap2_seq)} residues")
print(f"  ERAP1: {len(erap1_seq)} residues")
print(f"  IRAP:  {len(irap_seq)} residues")

def align_and_superimpose(ref_seq, ref_cas, ref_resis, mob_seq, mob_cas, mob_resis, label=""):
    """Align sequences, pair CAs, superimpose mobile onto reference."""
    # Global alignment
    alignments = pairwise2.align.globalxx(ref_seq, mob_seq, one_alignment_only=True)
    aln = alignments[0]
    ref_aln, mob_aln = aln.seqA, aln.seqB

    # Pair matched positions
    ref_paired = []
    mob_paired = []
    ref_idx = 0
    mob_idx = 0

    for i in range(len(ref_aln)):
        ref_char = ref_aln[i]
        mob_char = mob_aln[i]

        if ref_char != '-' and mob_char != '-':
            ref_paired.append(ref_cas[ref_idx])
            mob_paired.append(mob_cas[mob_idx])

        if ref_char != '-':
            ref_idx += 1
        if mob_char != '-':
            mob_idx += 1

    print(f"  {label}: {len(ref_paired)} paired CAs")

    # Superimpose
    sup = Superimposer()
    ref_atoms = ref_paired
    mob_atoms = mob_paired

    # Get coordinates
    sup.set_atoms(ref_atoms, mob_atoms)

    # Apply rotation/translation to ALL atoms in mobile structure
    all_atoms = list(mob_resis[0].get_parent().get_parent().get_parent().get_atoms())
    sup.apply(all_atoms)

    print(f"  {label}: RMSD = {sup.rms:.2f} A over {len(ref_paired)} CAs")

    return sup.rms

# ── Align ERAP1 and IRAP onto ERAP2 ────────────────────────────────────
print("\nPerforming structural alignment...")
rmsd1 = align_and_superimpose(erap2_seq, erap2_cas, erap2_resis,
                               erap1_seq, erap1_cas, erap1_resis, "ERAP1->ERAP2")
rmsd2 = align_and_superimpose(erap2_seq, erap2_cas, erap2_resis,
                               irap_seq, irap_cas, irap_resis, "IRAP->ERAP2")

# ── Find spatially equivalent residues ──────────────────────────────────
# After alignment, find the residue in ERAP1/IRAP whose CA is closest to each ERAP2 target CA

AA_PROPERTIES = {
    'A': ('nonpolar', 'small', 0),
    'R': ('positive', 'large', +1),
    'N': ('polar', 'medium', 0),
    'D': ('negative', 'medium', -1),
    'C': ('polar', 'small', 0),
    'E': ('negative', 'large', -1),
    'Q': ('polar', 'large', 0),
    'G': ('nonpolar', 'tiny', 0),
    'H': ('positive*', 'large', 0),  # partially charged at physiological pH
    'I': ('nonpolar', 'large', 0),
    'L': ('nonpolar', 'large', 0),
    'K': ('positive', 'large', +1),
    'M': ('nonpolar', 'large', 0),
    'F': ('nonpolar', 'large', 0),
    'P': ('nonpolar', 'medium', 0),
    'S': ('polar', 'small', 0),
    'T': ('polar', 'medium', 0),
    'W': ('nonpolar', 'large', 0),
    'Y': ('polar', 'large', 0),
    'V': ('nonpolar', 'medium', 0),
}

def get_residue_info(res):
    """Get residue number, name, one-letter code."""
    resname = res.get_resname()
    resnum = res.id[1]
    try:
        aa1 = three_to_one(resname)
    except KeyError:
        aa1 = '?'
    return resnum, resname, aa1

def find_nearest_residue(target_res, search_resis, n_closest=3):
    """Find the nearest residue(s) in search_resis to target_res, by CA distance."""
    target_ca = target_res['CA'].get_vector().get_array()

    distances = []
    for res in search_resis:
        if 'CA' not in res:
            continue
        ca = res['CA'].get_vector().get_array()
        dist = np.linalg.norm(target_ca - ca)
        distances.append((dist, res))

    distances.sort(key=lambda x: x[0])
    return distances[:n_closest]

def find_nearest_sidechain(target_res, search_resis, n_closest=3):
    """Find nearest residue by minimum sidechain atom distance (more precise than CA-CA)."""
    # Get sidechain heavy atoms of target
    target_atoms = []
    for atom in target_res:
        if atom.name not in ('N', 'CA', 'C', 'O') and atom.element != 'H':
            target_atoms.append(atom.get_vector().get_array())
    if not target_atoms:
        # Glycine - use CA
        target_atoms = [target_res['CA'].get_vector().get_array()]

    distances = []
    for res in search_resis:
        sc_atoms = []
        for atom in res:
            if atom.name not in ('N', 'CA', 'C', 'O') and atom.element != 'H':
                sc_atoms.append(atom.get_vector().get_array())
        if not sc_atoms:
            if 'CA' in res:
                sc_atoms = [res['CA'].get_vector().get_array()]
            else:
                continue

        # Min distance between any sidechain atom pair
        min_dist = float('inf')
        for ta in target_atoms:
            for sa in sc_atoms:
                d = np.linalg.norm(ta - sa)
                if d < min_dist:
                    min_dist = d

        distances.append((min_dist, res))

    distances.sort(key=lambda x: x[0])
    return distances[:n_closest]

# ── Build the residue map for ERAP2 target residues ────────────────────
# Re-extract post-alignment coordinates
erap1_resis_post = []
for res in erap1_struct[0][ERAP1_CHAIN]:
    if is_aa(res, standard=True) and 'CA' in res:
        erap1_resis_post.append(res)

irap_resis_post = []
for res in irap_struct[0][IRAP_CHAIN]:
    if is_aa(res, standard=True) and 'CA' in res:
        irap_resis_post.append(res)

# Find ERAP2 target residues
erap2_target_resis = {}
for res in erap2_struct[0][ERAP2_CHAIN]:
    if is_aa(res, standard=True) and res.id[1] in ERAP2_RESIDUES:
        erap2_target_resis[res.id[1]] = res

print("\n" + "="*100)
print("STRUCTURAL EQUIVALENCE TABLE: ERAP2 S1' POCKET RESIDUES")
print("="*100)
print(f"{'ERAP2':^12} | {'ERAP2':^8} | {'ERAP1 equivalent':^25} | {'IRAP equivalent':^25} | {'Chemistry Notes'}")
print(f"{'position':^12} | {'residue':^8} | {'(resi + AA + dist)':^25} | {'(resi + AA + dist)':^25} |")
print("-"*100)

results = []

for resi_num in ERAP2_RESIDUES:
    if resi_num not in erap2_target_resis:
        print(f"  WARNING: ERAP2 residue {resi_num} not found in structure!")
        continue

    target = erap2_target_resis[resi_num]
    t_num, t_name, t_aa = get_residue_info(target)

    # Find nearest in ERAP1 (by CA and sidechain)
    erap1_ca_near = find_nearest_residue(target, erap1_resis_post, n_closest=5)
    erap1_sc_near = find_nearest_sidechain(target, erap1_resis_post, n_closest=5)

    # Find nearest in IRAP
    irap_ca_near = find_nearest_residue(target, irap_resis_post, n_closest=5)
    irap_sc_near = find_nearest_sidechain(target, irap_resis_post, n_closest=5)

    # Best match by CA distance
    e1_dist, e1_res = erap1_ca_near[0]
    e1_num, e1_name, e1_aa = get_residue_info(e1_res)

    ir_dist, ir_res = irap_ca_near[0]
    ir_num, ir_name, ir_aa = get_residue_info(ir_res)

    # Chemistry comparison
    t_props = AA_PROPERTIES.get(t_aa, ('?', '?', 0))
    e1_props = AA_PROPERTIES.get(e1_aa, ('?', '?', 0))
    ir_props = AA_PROPERTIES.get(ir_aa, ('?', '?', 0))

    notes = []
    if t_props[2] != 0:
        charge_str = "+" if t_props[2] > 0 else "-"
        if e1_props[2] != t_props[2]:
            notes.append(f"ERAP2 {charge_str} vs ERAP1 {'+' if e1_props[2]>0 else '-' if e1_props[2]<0 else '0'}")
        if ir_props[2] != t_props[2]:
            notes.append(f"ERAP2 {charge_str} vs IRAP {'+' if ir_props[2]>0 else '-' if ir_props[2]<0 else '0'}")
    if t_props[0] != e1_props[0] or t_props[0] != ir_props[0]:
        notes.append(f"polarity: E2={t_props[0]}, E1={e1_props[0]}, IR={ir_props[0]}")
    if not notes:
        notes.append("conserved chemistry")

    notes_str = "; ".join(notes)

    print(f"  {resi_num:^10} | {t_aa:^3}({t_name:3}) | {e1_aa}{e1_num} ({e1_dist:.1f}A CA-CA){' ':>8} | {ir_aa}{ir_num} ({ir_dist:.1f}A CA-CA){' ':>8} | {notes_str}")

    results.append({
        'erap2_resi': resi_num,
        'erap2_aa': t_aa,
        'erap2_name': t_name,
        'erap1_resi': e1_num,
        'erap1_aa': e1_aa,
        'erap1_name': e1_name,
        'erap1_ca_dist': e1_dist,
        'irap_resi': ir_num,
        'irap_aa': ir_aa,
        'irap_name': ir_name,
        'irap_ca_dist': ir_dist,
    })

# ── Detailed analysis for each position ─────────────────────────────────
print("\n" + "="*100)
print("DETAILED ANALYSIS PER POSITION")
print("="*100)

for resi_num in ERAP2_RESIDUES:
    if resi_num not in erap2_target_resis:
        continue

    target = erap2_target_resis[resi_num]
    t_num, t_name, t_aa = get_residue_info(target)

    print(f"\n--- ERAP2 {t_aa}{resi_num} ({t_name}) ---")

    # Top 3 nearest in ERAP1
    erap1_ca_near = find_nearest_residue(target, erap1_resis_post, n_closest=5)
    erap1_sc_near = find_nearest_sidechain(target, erap1_resis_post, n_closest=3)

    print(f"  ERAP1 nearest by CA-CA:")
    for dist, res in erap1_ca_near[:5]:
        rn, rname, raa = get_residue_info(res)
        props = AA_PROPERTIES.get(raa, ('?','?',0))
        print(f"    {raa}{rn:>4} ({rname}) - {dist:.2f}A  [{props[0]}, {props[1]}, charge={props[2]}]")

    print(f"  ERAP1 nearest by sidechain:")
    for dist, res in erap1_sc_near[:3]:
        rn, rname, raa = get_residue_info(res)
        print(f"    {raa}{rn:>4} ({rname}) - {dist:.2f}A sidechain")

    # Top 3 nearest in IRAP
    irap_ca_near = find_nearest_residue(target, irap_resis_post, n_closest=5)
    irap_sc_near = find_nearest_sidechain(target, irap_resis_post, n_closest=3)

    print(f"  IRAP nearest by CA-CA:")
    for dist, res in irap_ca_near[:5]:
        rn, rname, raa = get_residue_info(res)
        props = AA_PROPERTIES.get(raa, ('?','?',0))
        print(f"    {raa}{rn:>4} ({rname}) - {dist:.2f}A  [{props[0]}, {props[1]}, charge={props[2]}]")

    print(f"  IRAP nearest by sidechain:")
    for dist, res in irap_sc_near[:3]:
        rn, rname, raa = get_residue_info(res)
        print(f"    {raa}{rn:>4} ({rname}) - {dist:.2f}A sidechain")

# ── Key question: K392 salt bridge selectivity ──────────────────────────
print("\n" + "="*100)
print("KEY QUESTION: ERAP2 K392 UNIQUE POSITIVE CHARGE?")
print("="*100)

k392 = erap2_target_resis.get(392)
if k392:
    t_num, t_name, t_aa = get_residue_info(k392)
    print(f"\nERAP2 K392: Lysine (positive charge, +1, long flexible sidechain)")

    # ERAP1 equivalent
    erap1_near = find_nearest_residue(k392, erap1_resis_post, n_closest=3)
    print(f"\nERAP1 structurally equivalent residues (by CA proximity):")
    for dist, res in erap1_near:
        rn, rname, raa = get_residue_info(res)
        props = AA_PROPERTIES.get(raa, ('?','?',0))
        charge = "POSITIVE" if props[2] > 0 else "NEGATIVE" if props[2] < 0 else "NEUTRAL"
        print(f"  {raa}{rn} ({rname}): {dist:.2f}A, {props[0]}, charge={charge}")

    # IRAP equivalent
    irap_near = find_nearest_residue(k392, irap_resis_post, n_closest=3)
    print(f"\nIRAP structurally equivalent residues (by CA proximity):")
    for dist, res in irap_near:
        rn, rname, raa = get_residue_info(res)
        props = AA_PROPERTIES.get(raa, ('?','?',0))
        charge = "POSITIVE" if props[2] > 0 else "NEGATIVE" if props[2] < 0 else "NEUTRAL"
        print(f"  {raa}{rn} ({rname}): {dist:.2f}A, {props[0]}, charge={charge}")

    # Salt bridge analysis
    e1_best = erap1_near[0]
    ir_best = irap_near[0]
    e1_rn, e1_rname, e1_aa = get_residue_info(e1_best[1])
    ir_rn, ir_rname, ir_aa = get_residue_info(ir_best[1])

    e1_charge = AA_PROPERTIES.get(e1_aa, ('?','?',0))[2]
    ir_charge = AA_PROPERTIES.get(ir_aa, ('?','?',0))[2]

    print(f"\n*** SALT BRIDGE SELECTIVITY ASSESSMENT ***")
    print(f"  ERAP2 K392 (Lys, +1) — would form salt bridge with Glu/Asp peptide")
    print(f"  ERAP1 equivalent: {e1_aa}{e1_rn} (charge={'+1' if e1_charge>0 else '-1' if e1_charge<0 else '0'})")
    print(f"  IRAP  equivalent: {ir_aa}{ir_rn} (charge={'+1' if ir_charge>0 else '-1' if ir_charge<0 else '0'})")

    if e1_charge <= 0 and ir_charge <= 0:
        print(f"\n  >> YES: K392 positive charge is UNIQUE to ERAP2 at this position.")
        print(f"  >> A Glu/Asp in the peptide substrate WOULD form a selective salt bridge")
        print(f"     only with ERAP2, NOT with ERAP1 or IRAP.")
    elif e1_charge > 0 or ir_charge > 0:
        which = []
        if e1_charge > 0:
            which.append(f"ERAP1 ({e1_aa}{e1_rn})")
        if ir_charge > 0:
            which.append(f"IRAP ({ir_aa}{ir_rn})")
        print(f"\n  >> CAUTION: Positive charge also present in: {', '.join(which)}")
        print(f"  >> Salt bridge selectivity may be REDUCED at this position.")

# ── PyMOL reproduction commands ─────────────────────────────────────────
print("\n" + "="*100)
print("PyMOL COMMANDS TO REPRODUCE THIS ANALYSIS")
print("="*100)
print("""
# Fetch structures
fetch 7SH0, erap2
fetch 6RQX, erap1
fetch 4PJ6, irap

# Align to ERAP2
align erap1, erap2
align irap, erap2

# Select ERAP2 target residues
select erap2_targets, erap2 and chain A and resi 392+398+403+406+412+414

# Show target residues
show sticks, erap2_targets
color yellow, erap2_targets

# For each ERAP2 residue, find nearest in ERAP1/IRAP:
""")

# Print specific select/distance commands for each result
for r in results:
    e2r = r['erap2_resi']
    e1r = r['erap1_resi']
    irr = r['irap_resi']
    print(f"# ERAP2 {r['erap2_aa']}{e2r} -> ERAP1 {r['erap1_aa']}{e1r}, IRAP {r['irap_aa']}{irr}")
    print(f"select e2_{e2r}, erap2 and chain A and resi {e2r}")
    print(f"select e1_{e1r}, erap1 and chain A and resi {e1r}")
    print(f"select ir_{irr}, irap and chain A and resi {irr}")
    print(f"distance dist_e1_{e2r}, e2_{e2r} and name CA, e1_{e1r} and name CA")
    print(f"distance dist_ir_{e2r}, e2_{e2r} and name CA, ir_{irr} and name CA")
    print(f"show sticks, e1_{e1r} or ir_{irr}")
    print()

print("# Color by enzyme")
print("color cyan, erap2")
print("color green, erap1")
print("color salmon, irap")
print("color yellow, erap2_targets")
print()
print("# Zoom to pocket")
print("zoom erap2_targets, 12")

print("\n\nDone.")
