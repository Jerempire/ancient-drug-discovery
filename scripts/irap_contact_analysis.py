"""
IRAP Contact Analysis: Diagnose where IRAP binding comes from in C2 EAAAK.

Decomposes the C2 construct (113aa) into 3 regions:
  - Cargo: residues 1-6 (VAGSAF)
  - Linker: residues 7-21 (EAAAKEAAAKEAAAK)
  - Binder body: residues 22-113 (V19W parent binder)

Compares contact maps of C2 vs ERAP2 (best, model_1 ipTM=0.810)
and C2 vs IRAP (best, model_0 ipTM=0.709) to identify what drives
IRAP cross-reactivity.
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import math
from pathlib import Path

# Use gemmi if available, otherwise fallback to manual CIF parsing
try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False
    print("WARNING: gemmi not available, using basic CIF parser")


RESULTS = Path("C:/Users/jmj2z/Projects/medical/ancient-drug-discovery/results/vagsaf_screen/results")

# C2 regions (1-indexed residue numbers in chain B)
CARGO = range(1, 7)       # VAGSAF (res 1-6)
LINKER = range(7, 22)     # EAAAKEAAAKEAAAK (res 7-21)
BINDER = range(22, 114)   # parent binder body (res 22-113)

CONTACT_CUTOFF = 5.0  # Angstroms for heavy-atom contact


def region_label(resnum):
    if resnum in CARGO:
        return "cargo"
    elif resnum in LINKER:
        return "linker"
    elif resnum in BINDER:
        return "binder"
    return "unknown"


def get_contacts_gemmi(cif_path, cutoff=CONTACT_CUTOFF):
    """Extract inter-chain contacts between chain A (target) and chain B (C2 construct)."""
    st = gemmi.read_structure(str(cif_path))
    model = st[0]

    chain_a = model.find_chain("A")
    chain_b = model.find_chain("B")

    if not chain_a or not chain_b:
        print(f"  Missing chains in {cif_path}")
        return [], [], []

    # Get all heavy atoms for each chain
    a_atoms = []
    for res in chain_a:
        for atom in res:
            if atom.element.name != "H":
                a_atoms.append((res.seqid.num, res.name, atom.name,
                               atom.pos.x, atom.pos.y, atom.pos.z))

    b_atoms = []
    for res in chain_b:
        for atom in res:
            if atom.element.name != "H":
                b_atoms.append((res.seqid.num, res.name, atom.name,
                               atom.pos.x, atom.pos.y, atom.pos.z))

    # Find contacts
    contacts = []
    for b_resnum, b_resname, b_atom, bx, by, bz in b_atoms:
        for a_resnum, a_resname, a_atom, ax, ay, az in a_atoms:
            d = math.sqrt((bx-ax)**2 + (by-ay)**2 + (bz-az)**2)
            if d <= cutoff:
                contacts.append({
                    "b_res": b_resnum,
                    "b_resname": b_resname,
                    "b_atom": b_atom,
                    "a_res": a_resnum,
                    "a_resname": a_resname,
                    "a_atom": a_atom,
                    "dist": d,
                    "region": region_label(b_resnum),
                })

    return contacts, a_atoms, b_atoms


def analyze_contacts(contacts):
    """Summarize contacts by region."""
    from collections import defaultdict

    by_region = defaultdict(list)
    for c in contacts:
        by_region[c["region"]].append(c)

    # Unique residue-residue pairs per region
    unique_pairs = defaultdict(set)
    for c in contacts:
        unique_pairs[c["region"]].add((c["b_res"], c["a_res"]))

    # Unique target residues contacted per region
    target_residues = defaultdict(set)
    for c in contacts:
        target_residues[c["region"]].add(c["a_res"])

    # Unique binder residues making contacts per region
    binder_residues = defaultdict(set)
    for c in contacts:
        binder_residues[c["region"]].add(c["b_res"])

    return by_region, unique_pairs, target_residues, binder_residues


def print_analysis(name, contacts):
    by_region, unique_pairs, target_res, binder_res = analyze_contacts(contacts)

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"  Total atom-atom contacts (<{CONTACT_CUTOFF}A): {len(contacts)}")
    print(f"{'='*60}")

    for region in ["cargo", "linker", "binder"]:
        n_contacts = len(by_region[region])
        n_pairs = len(unique_pairs[region])
        n_target = len(target_res[region])
        n_binder = len(binder_res[region])
        pct = 100 * n_contacts / max(1, len(contacts))

        print(f"\n  {region.upper():8s}: {n_contacts:4d} contacts ({pct:5.1f}%)")
        print(f"           {n_pairs:4d} unique res-res pairs")
        print(f"           {n_binder:4d} C2 residues involved")
        print(f"           {n_target:4d} target residues contacted")

        # List the specific binder residues making contacts
        if n_binder > 0 and region in ["cargo", "linker"]:
            sorted_bres = sorted(binder_res[region])
            print(f"           C2 residues: {sorted_bres}")

        # Show closest contacts
        if by_region[region]:
            closest = sorted(by_region[region], key=lambda x: x["dist"])[:5]
            print(f"           Closest contacts:")
            for c in closest:
                print(f"             B:{c['b_resname']}{c['b_res']}.{c['b_atom']} - "
                      f"A:{c['a_resname']}{c['a_res']}.{c['a_atom']} = {c['dist']:.1f}A")


def compare_target_contacts(erap2_contacts, irap_contacts):
    """Find IRAP-specific contacts not present in ERAP2."""
    print(f"\n{'='*60}")
    print(f"  DIFFERENTIAL ANALYSIS: IRAP-specific contacts")
    print(f"{'='*60}")

    # Get residue-residue pairs for each
    erap2_pairs = set()
    for c in erap2_contacts:
        erap2_pairs.add((c["b_res"], c["region"]))

    irap_pairs = set()
    for c in irap_contacts:
        irap_pairs.add((c["b_res"], c["region"]))

    # C2 residues that contact IRAP but NOT ERAP2
    irap_only = set()
    for c in irap_contacts:
        key = (c["b_res"], c["region"])
        if c["b_res"] not in {p[0] for p in erap2_pairs}:
            irap_only.add((c["b_res"], c["b_resname"], c["region"]))

    if irap_only:
        print(f"\n  C2 residues that contact IRAP but NOT ERAP2:")
        for res, name, region in sorted(irap_only):
            print(f"    {name}{res} ({region})")

    # C2 residues that contact ERAP2 but NOT IRAP
    erap2_only = set()
    for c in erap2_contacts:
        if c["b_res"] not in {p[0] for p in irap_pairs}:
            erap2_only.add((c["b_res"], c["b_resname"], c["region"]))

    if erap2_only:
        print(f"\n  C2 residues that contact ERAP2 but NOT IRAP:")
        for res, name, region in sorted(erap2_only):
            print(f"    {name}{res} ({region})")

    # Shared contacts (same C2 residue contacts both targets)
    shared = set()
    erap2_bres = {c["b_res"] for c in erap2_contacts}
    irap_bres = {c["b_res"] for c in irap_contacts}
    both = erap2_bres & irap_bres
    print(f"\n  Shared C2 residues (contact both): {len(both)}")
    print(f"  ERAP2-only C2 residues: {len(erap2_bres - irap_bres)}")
    print(f"  IRAP-only C2 residues: {len(irap_bres - erap2_bres)}")


def main():
    # Best models: ERAP2 model_1 (ipTM=0.810), IRAP model_0 (ipTM=0.709), ERAP1 model_0 (ipTM=0.315)
    erap2_cif = RESULTS / "c2_eaaak_vs_erap2k392/boltz_results_c2_eaaak_vs_erap2k392/predictions/c2_eaaak_vs_erap2k392/c2_eaaak_vs_erap2k392_model_1.cif"
    irap_cif = RESULTS / "c2_eaaak_vs_irap/boltz_results_c2_eaaak_vs_irap/predictions/c2_eaaak_vs_irap/c2_eaaak_vs_irap_model_0.cif"
    erap1_cif = RESULTS / "c2_eaaak_vs_erap1/boltz_results_c2_eaaak_vs_erap1/predictions/c2_eaaak_vs_erap1/c2_eaaak_vs_erap1_model_0.cif"

    if not HAS_GEMMI:
        print("Need gemmi to parse CIF files. Install with: pip install gemmi")
        return

    print("IRAP CONTACT ANALYSIS FOR C2 EAAAK")
    print("Decomposing 113aa construct: cargo(1-6) | linker(7-21) | binder(22-113)")
    print(f"Contact cutoff: {CONTACT_CUTOFF} Angstroms")

    # Analyze each target
    erap2_contacts, _, _ = get_contacts_gemmi(erap2_cif)
    print_analysis("C2 vs ERAP2 (model_1, ipTM=0.810)", erap2_contacts)

    irap_contacts, _, _ = get_contacts_gemmi(irap_cif)
    print_analysis("C2 vs IRAP (model_0, ipTM=0.709)", irap_contacts)

    erap1_contacts, _, _ = get_contacts_gemmi(erap1_cif)
    print_analysis("C2 vs ERAP1 (model_0, ipTM=0.315)", erap1_contacts)

    # Differential analysis
    compare_target_contacts(erap2_contacts, irap_contacts)

    # Summary
    print(f"\n{'='*60}")
    print(f"  DIAGNOSIS SUMMARY")
    print(f"{'='*60}")

    for label, contacts in [("ERAP2", erap2_contacts), ("IRAP", irap_contacts), ("ERAP1", erap1_contacts)]:
        by_region, _, _, _ = analyze_contacts(contacts)
        total = len(contacts)
        cargo_pct = 100 * len(by_region["cargo"]) / max(1, total)
        linker_pct = 100 * len(by_region["linker"]) / max(1, total)
        binder_pct = 100 * len(by_region["binder"]) / max(1, total)
        print(f"\n  {label:6s}: cargo={cargo_pct:5.1f}%  linker={linker_pct:5.1f}%  binder={binder_pct:5.1f}%  (total={total})")


if __name__ == "__main__":
    main()
