"""
V2 Protein Binder Channel-Capping Analysis (Terminal E)

Analyzes whether the V2 binder (n248_trim_c5_Y87A_Y89A, 92aa) caps the
ERAP2 substrate channel entrance, or just binds the surface without blocking
substrate access.

Usage:
    python scripts/v2_channel_capping_analysis.py

Output: data/results/v43_validation/terminal_e/capping_analysis.json
"""
import json
import os
import sys
import warnings
from datetime import datetime, timezone

import numpy as np
from Bio.PDB import NeighborSearch
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

# === Config ===
PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF_FILE = os.path.join(
    PROJECT,
    "data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_erap2.cif",
)
OUTPUT_DIR = os.path.join(
    PROJECT, "data/results/v43_validation/terminal_e"
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ERAP2 residue offset: chain A residue 1 = ERAP2 residue 350
OFFSET = 350

# Key structural regions (original ERAP2 numbering)
CHANNEL_ENTRANCE = list(range(353, 368)) + list(range(400, 415))
ZINC_COORD = [370, 374, 393]
CHANNEL_INTERIOR = list(range(368, 400))
K392 = [392]

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def to_chain_ids(reslist):
    return set(r - OFFSET + 1 for r in reslist)


def compute_region_sasa(apo_sasa, holo_sasa, region_ids):
    apo_total = sum(apo_sasa.get(r, 0) for r in region_ids)
    holo_total = sum(holo_sasa.get(r, 0) for r in region_ids)
    reduction = (apo_total - holo_total) / apo_total * 100 if apo_total > 0 else 0
    if reduction > 80:
        verdict = "CAPPED"
    elif reduction > 50:
        verdict = "PARTIALLY_CAPPED"
    else:
        verdict = "NOT_CAPPED"
    return {
        "apo_sasa": round(apo_total, 1),
        "holo_sasa": round(holo_total, 1),
        "reduction_pct": round(reduction, 1),
        "verdict": verdict,
    }


def main():
    parser = MMCIFParser(QUIET=True)

    # Load complex (holo)
    struct_holo = parser.get_structure("holo", CIF_FILE)
    model_holo = struct_holo[0]
    chain_a = model_holo["A"]
    chain_b = model_holo["B"]

    # Compute SASA with binder
    sr = ShrakeRupley()
    sr.compute(model_holo, level="R")
    holo_sasa = {res.get_id()[1]: res.sasa for res in chain_a}

    # Load apo (remove binder)
    struct_apo = parser.get_structure("apo", CIF_FILE)
    struct_apo[0].detach_child("B")
    sr2 = ShrakeRupley()
    sr2.compute(struct_apo[0], level="R")
    apo_sasa = {res.get_id()[1]: res.sasa for res in struct_apo[0]["A"]}

    # Region definitions
    entrance_ids = to_chain_ids(CHANNEL_ENTRANCE)
    zinc_ids = to_chain_ids(ZINC_COORD)
    interior_ids = to_chain_ids(CHANNEL_INTERIOR)
    all_channel_ids = entrance_ids | zinc_ids | interior_ids

    # SASA by region
    results = {
        "binder": "n248_trim_c5_Y87A_Y89A",
        "binder_length": 92,
        "target": "ERAP2 (cropped 350-500)",
        "cif_file": os.path.basename(CIF_FILE),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "terminal": "e",
    }

    results["sasa"] = {
        "channel_entrance": compute_region_sasa(apo_sasa, holo_sasa, entrance_ids),
        "zinc_coordinating": compute_region_sasa(apo_sasa, holo_sasa, zinc_ids),
        "channel_interior": compute_region_sasa(apo_sasa, holo_sasa, interior_ids),
        "entire_channel": compute_region_sasa(apo_sasa, holo_sasa, all_channel_ids),
    }

    # Interface contacts
    erap2_atoms = list(chain_a.get_atoms())
    binder_atoms = list(chain_b.get_atoms())
    ns = NeighborSearch(erap2_atoms)

    contact_residues_a = set()
    contact_residues_b = set()
    for atom in binder_atoms:
        close = ns.search(atom.get_vector().get_array(), 4.5, "R")
        for res in close:
            contact_residues_a.add(res.get_id()[1])
            contact_residues_b.add(atom.get_parent().get_id()[1])

    contact_orig = sorted([r + OFFSET - 1 for r in contact_residues_a])
    results["interface"] = {
        "erap2_contact_residues": contact_orig,
        "binder_contact_residues": sorted(contact_residues_b),
        "n_erap2_contacts": len(contact_residues_a),
        "n_binder_contacts": len(contact_residues_b),
        "contacts_divergent_patch_1": [r for r in contact_orig if 353 <= r <= 367],
        "contacts_divergent_patch_2": [r for r in contact_orig if 400 <= r <= 414],
        "contacts_zinc_region": [r for r in contact_orig if r in [370, 374, 393]],
        "contacts_k392": 392 in contact_orig,
        "contacts_interior": [r for r in contact_orig if 368 <= r <= 399],
    }

    # Uncapped entrance residues
    uncapped = []
    for rid in sorted(entrance_ids):
        orig = rid + OFFSET - 1
        a = apo_sasa.get(rid, 0)
        h = holo_sasa.get(rid, 0)
        if a > 5 and h / a > 0.5:
            uncapped.append({
                "residue": orig,
                "apo_sasa": round(a, 1),
                "holo_sasa": round(h, 1),
                "pct_exposed": round(h / a * 100, 0),
            })
    results["uncapped_entrance_residues"] = uncapped

    # Geometric analysis
    entrance_coords = []
    for res in chain_a:
        if res.get_id()[1] in entrance_ids:
            for atom in res:
                entrance_coords.append(atom.get_vector().get_array())
    entrance_center = np.array(entrance_coords).mean(axis=0)

    binder_coords = np.array([a.get_vector().get_array() for a in binder_atoms])
    binder_com = binder_coords.mean(axis=0)

    interior_coords = []
    for res in chain_a:
        if res.get_id()[1] in interior_ids:
            for atom in res:
                interior_coords.append(atom.get_vector().get_array())
    interior_center = np.array(interior_coords).mean(axis=0)

    zinc_coords = []
    for res in chain_a:
        if res.get_id()[1] in zinc_ids:
            for atom in res:
                zinc_coords.append(atom.get_vector().get_array())
    zinc_center = np.array(zinc_coords).mean(axis=0)

    results["geometry"] = {
        "binder_com_to_entrance": round(float(np.linalg.norm(binder_com - entrance_center)), 1),
        "binder_com_to_interior": round(float(np.linalg.norm(binder_com - interior_center)), 1),
        "binder_com_to_zinc": round(float(np.linalg.norm(binder_com - zinc_center)), 1),
        "nterm_to_entrance": round(float(np.linalg.norm(
            chain_b[(" ", 1, " ")]["CA"].get_vector().get_array() - entrance_center
        )), 1),
        "cterm_to_entrance": round(float(np.linalg.norm(
            chain_b[(" ", 92, " ")]["CA"].get_vector().get_array() - entrance_center
        )), 1),
    }

    # Binder residues near entrance (extension candidates)
    extension_candidates = []
    entrance_atom_arr = np.array(entrance_coords)
    for res in chain_b:
        rid = res.get_id()[1]
        rname = res.get_resname()
        ca = res["CA"] if "CA" in res else list(res.get_atoms())[0]
        coord = ca.get_vector().get_array()
        dists = np.linalg.norm(entrance_atom_arr - coord, axis=1)
        min_d = float(dists.min())

        if min_d < 12.0:
            contacts_erap2 = False
            for atom in res:
                close = ns.search(atom.get_vector().get_array(), 4.5, "R")
                if close:
                    contacts_erap2 = True
                    break
            if not contacts_erap2:
                extension_candidates.append({
                    "binder_residue": rid,
                    "aa": AA3TO1.get(rname, "X"),
                    "distance_to_entrance": round(min_d, 1),
                    "is_solvent_exposed": True,
                })

    results["extension_candidates"] = extension_candidates

    # Binder sequence
    seq = "".join(
        AA3TO1.get(res.get_resname(), "X") for res in chain_b
    )
    results["binder_sequence"] = seq

    # Overall verdict
    entrance_red = results["sasa"]["channel_entrance"]["reduction_pct"]
    interior_red = results["sasa"]["channel_interior"]["reduction_pct"]

    if entrance_red > 80 and interior_red > 80:
        overall = "CHANNEL_CAP"
        recommendation = "Binder caps both entrance and interior — functional inhibitor"
    elif entrance_red > 80:
        overall = "ENTRANCE_CAP"
        recommendation = "Binder caps entrance but interior exposed — may still inhibit"
    elif interior_red > 50:
        overall = "INTERIOR_PLUG"
        recommendation = (
            "Binder occludes interior but entrance open — substrates can still enter. "
            "Need to extend binder to cover entrance residues 353-363 and 409-414."
        )
    else:
        overall = "SURFACE_STICKER"
        recommendation = "Binder sits on surface without capping — needs major redesign"

    results["verdict"] = {
        "classification": overall,
        "entrance_reduction_pct": entrance_red,
        "interior_reduction_pct": interior_red,
        "recommendation": recommendation,
    }

    # Write output
    out_path = os.path.join(OUTPUT_DIR, "capping_analysis.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    print("=" * 60)
    print("V2 BINDER CHANNEL-CAPPING ANALYSIS")
    print("=" * 60)
    print(f"Binder: {results['binder']} ({results['binder_length']}aa)")
    print()
    print("SASA Reduction by Region:")
    for region, data in results["sasa"].items():
        print(f"  {region:25s}: {data['reduction_pct']:5.1f}% -> {data['verdict']}")
    print()
    print(f"Verdict: {overall}")
    print(f"  {recommendation}")
    print()
    print(f"Uncapped entrance residues: {len(uncapped)}")
    print(f"Extension candidate positions: {[c['binder_residue'] for c in extension_candidates]}")
    print()
    print(f"Output: {out_path}")

    return results


if __name__ == "__main__":
    main()
