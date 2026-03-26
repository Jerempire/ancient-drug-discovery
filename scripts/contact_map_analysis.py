"""
Contact Map Analysis — V4.3 Wet Lab Validation

For each Elite Four peptide vs ERAP2-K392:
1. Compute all peptide-protein contacts (< 4.5 A)
2. Classify contacts by 3-region model (Floor/Wall/Ceiling)
3. Identify key interactions (H-bonds, hydrophobic)
4. Output JSON contact maps for each peptide

Region definitions from V43_WETLAB_VALIDATION_PLAN.md:
  Floor:   K392, Y398 — variant selectivity
  Wall:    A403, A406 — paralog selectivity
  Ceiling: Q412, D414 — channel cap / volume filling
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import numpy as np
from pathlib import Path
from datetime import datetime, timezone
from Bio.PDB import PDBParser

PROJECT = Path(r"C:\Users\jmj2z\Projects\medical\ancient-drug-discovery")
STARTING_STRUCTS = PROJECT / "data" / "results" / "v43_validation" / "md_starting_structures"
OUTPUT_DIR = PROJECT / "data" / "results" / "v43_validation" / "terminal_a" / "contact_maps"

PEPTIDES = ["VAGSAF", "IAFSAF", "VAWSAF", "FASGAV"]
TARGETS = ["erap2k392"]  # contact maps for primary target

CONTACT_CUTOFF = 4.5  # Angstroms
HBOND_CUTOFF = 3.5    # Angstroms, donor-acceptor distance
PROTEIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"

# 3-region model residues
REGIONS = {
    "floor": {"residues": [392, 398], "description": "K392 selectivity handle"},
    "wall":  {"residues": [403, 406], "description": "IRAP/ERAP1 evasion handle"},
    "ceiling": {"residues": [412, 414], "description": "Channel cap / volume filling"}
}

# Expanded contact neighborhood (within ~8A of active site channel)
CHANNEL_RESIDUES = list(range(365, 420))

# H-bond capable atoms
HBOND_DONORS = {"N", "NE", "NE2", "ND1", "ND2", "NH1", "NH2", "NZ", "OG", "OG1", "OH", "NE1"}
HBOND_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD"}


def get_contacts(structure, cutoff=CONTACT_CUTOFF):
    """Find all peptide-protein atom pairs within cutoff distance."""
    model = structure[0]
    protein_atoms = []
    peptide_atoms = []

    for chain in model:
        if chain.id == PROTEIN_CHAIN:
            for res in chain:
                if res.id[0] != " ":
                    continue  # skip heteroatoms
                for atom in res:
                    if atom.element != "H":
                        protein_atoms.append((atom, res))
        elif chain.id == PEPTIDE_CHAIN:
            for res in chain:
                for atom in res:
                    if atom.element != "H":
                        peptide_atoms.append((atom, res))

    contacts = []
    for patom, pres in protein_atoms:
        pcoord = patom.get_vector().get_array()
        for qatom, qres in peptide_atoms:
            qcoord = qatom.get_vector().get_array()
            dist = np.linalg.norm(pcoord - qcoord)
            if dist <= cutoff:
                contacts.append({
                    "protein_chain": PROTEIN_CHAIN,
                    "protein_resname": pres.resname,
                    "protein_resid": pres.id[1],
                    "protein_atom": patom.name,
                    "peptide_chain": PEPTIDE_CHAIN,
                    "peptide_resname": qres.resname,
                    "peptide_resid": qres.id[1],
                    "peptide_atom": qatom.name,
                    "distance_A": round(dist, 2),
                    "is_hbond_candidate": (
                        (patom.name in HBOND_DONORS and qatom.name in HBOND_ACCEPTORS) or
                        (patom.name in HBOND_ACCEPTORS and qatom.name in HBOND_DONORS)
                    ) and dist <= HBOND_CUTOFF
                })

    return contacts


def classify_contacts(contacts):
    """Classify contacts by 3-region model."""
    region_contacts = {region: [] for region in REGIONS}
    region_contacts["other_channel"] = []
    region_contacts["outside_channel"] = []

    for contact in contacts:
        resid = contact["protein_resid"]
        assigned = False
        for region_name, region_def in REGIONS.items():
            if resid in region_def["residues"]:
                region_contacts[region_name].append(contact)
                assigned = True
                break
        if not assigned:
            if resid in CHANNEL_RESIDUES:
                region_contacts["other_channel"].append(contact)
            else:
                region_contacts["outside_channel"].append(contact)

    return region_contacts


def get_region_verdict(region_contacts):
    """Classify binding mode based on 3-region contacts."""
    n_regions = sum(1 for r in ["floor", "wall", "ceiling"] if len(region_contacts[r]) > 0)
    if n_regions == 3:
        return "VOLUME_FILLER"
    elif n_regions == 2:
        return "PARTIAL"
    else:
        return "LOCAL_ANCHOR"


def summarize_contacts(contacts):
    """Get unique protein residues contacted."""
    seen = set()
    unique = []
    for c in contacts:
        key = (c["protein_resid"], c["protein_resname"])
        if key not in seen:
            seen.add(key)
            unique.append({"resid": c["protein_resid"], "resname": c["protein_resname"]})
    return sorted(unique, key=lambda x: x["resid"])


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)
    all_results = []

    for peptide in PEPTIDES:
        for target in TARGETS:
            pdb_file = STARTING_STRUCTS / f"{peptide}_vs_{target}.pdb"
            if not pdb_file.exists():
                print(f"WARNING: {pdb_file} not found, skipping")
                continue

            print(f"\n{'='*60}")
            print(f"{peptide} vs {target}")
            print(f"{'='*60}")

            structure = parser.get_structure(f"{peptide}_{target}", str(pdb_file))
            contacts = get_contacts(structure)
            region_contacts = classify_contacts(contacts)

            # H-bond candidates
            hbonds = [c for c in contacts if c["is_hbond_candidate"]]

            # Region summary
            region_verdict = get_region_verdict(region_contacts)
            print(f"Total contacts (< {CONTACT_CUTOFF} A): {len(contacts)}")
            print(f"H-bond candidates: {len(hbonds)}")
            print(f"Region verdict: {region_verdict}")

            for region_name in ["floor", "wall", "ceiling", "other_channel", "outside_channel"]:
                rc = region_contacts[region_name]
                unique_res = summarize_contacts(rc)
                label = REGIONS.get(region_name, {}).get("description", region_name)
                print(f"  {region_name:>16} ({label}): {len(rc)} contacts, "
                      f"{len(unique_res)} residues")
                if unique_res:
                    for r in unique_res:
                        print(f"    {r['resname']}{r['resid']}")

            # Key interactions
            print(f"\nKey H-bond candidates:")
            for hb in sorted(hbonds, key=lambda x: x["distance_A"]):
                print(f"  {hb['protein_resname']}{hb['protein_resid']}.{hb['protein_atom']} "
                      f"<-> {hb['peptide_resname']}{hb['peptide_resid']}.{hb['peptide_atom']} "
                      f"= {hb['distance_A']:.2f} A")

            # K392 specific contacts
            k392_contacts = [c for c in contacts if c["protein_resid"] == 392]
            print(f"\nK392 contacts: {len(k392_contacts)}")
            for c in sorted(k392_contacts, key=lambda x: x["distance_A"])[:5]:
                print(f"  K392.{c['protein_atom']} <-> {c['peptide_resname']}{c['peptide_resid']}.{c['peptide_atom']} = {c['distance_A']:.2f} A")

            result = {
                "peptide": peptide,
                "target": target,
                "total_contacts": len(contacts),
                "hbond_candidates": len(hbonds),
                "region_verdict": region_verdict,
                "floor_contacts": len(region_contacts["floor"]),
                "wall_contacts": len(region_contacts["wall"]),
                "ceiling_contacts": len(region_contacts["ceiling"]),
                "other_channel_contacts": len(region_contacts["other_channel"]),
                "outside_channel_contacts": len(region_contacts["outside_channel"]),
                "floor_residues": summarize_contacts(region_contacts["floor"]),
                "wall_residues": summarize_contacts(region_contacts["wall"]),
                "ceiling_residues": summarize_contacts(region_contacts["ceiling"]),
                "k392_contacts": len(k392_contacts),
                "k392_closest_A": min((c["distance_A"] for c in k392_contacts), default=None),
                "hbonds": [{
                    "protein": f"{h['protein_resname']}{h['protein_resid']}.{h['protein_atom']}",
                    "peptide": f"{h['peptide_resname']}{h['peptide_resid']}.{h['peptide_atom']}",
                    "distance_A": h["distance_A"]
                } for h in sorted(hbonds, key=lambda x: x["distance_A"])],
                "unique_protein_residues_contacted": summarize_contacts(contacts)
            }

            # Write individual result
            out_file = OUTPUT_DIR / f"{peptide}_vs_{target}_contacts.json"
            with open(out_file, "w") as f:
                json.dump(result, f, indent=2)

            all_results.append(result)

    # Write summary
    summary = {
        "test": "Contact Map Analysis",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "contact_cutoff_A": CONTACT_CUTOFF,
        "hbond_cutoff_A": HBOND_CUTOFF,
        "regions": {k: v for k, v in REGIONS.items()},
        "results": [{
            "peptide": r["peptide"],
            "target": r["target"],
            "total_contacts": r["total_contacts"],
            "hbond_candidates": r["hbond_candidates"],
            "k392_contacts": r["k392_contacts"],
            "k392_closest_A": r["k392_closest_A"],
            "floor": r["floor_contacts"],
            "wall": r["wall_contacts"],
            "ceiling": r["ceiling_contacts"],
            "region_verdict": r["region_verdict"]
        } for r in all_results]
    }

    summary_file = OUTPUT_DIR / "contact_map_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print(f"CONTACT MAP SUMMARY")
    print(f"{'='*60}")
    print(f"{'Peptide':<10} {'Total':>6} {'HBond':>6} {'K392':>5} {'Floor':>6} {'Wall':>5} {'Ceil':>5} {'Verdict':<15}")
    print(f"{'-'*65}")
    for r in all_results:
        print(f"{r['peptide']:<10} {r['total_contacts']:>6} {r['hbond_candidates']:>6} "
              f"{r['k392_contacts']:>5} {r['floor_contacts']:>6} {r['wall_contacts']:>5} "
              f"{r['ceiling_contacts']:>5} {r['region_verdict']:<15}")


if __name__ == "__main__":
    main()
