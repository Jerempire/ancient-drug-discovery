"""
selectivity_screen_v3.py — Fixed ERAP2 vs ERAP1 selectivity with structural alignment.

Key fix: Align ERAP1 onto ERAP2 using conserved CA atoms FIRST, then measure
binder cross-reactivity in the shared coordinate frame.

Scores all RFdiffusion binder designs (old active-site + new divergent campaigns).
"""
import csv
import glob
import json
import os

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Superimposer

# --- Config ---
ERAP2_PDB = "/workspace/data/structures/erap2_wt_alphafold.pdb"
ERAP1_PDB = "/workspace/data/structures/erap1_wt_alphafold.pdb"
DIVERGENT_JSON = "/workspace/results/erap2_divergent_surface.json"
BINDER_GLOBS = [
    "/workspace/results/rfdiffusion/erap2_sel_*.pdb",
    "/workspace/results/rfdiffusion/erap2_short_*.pdb",
    "/workspace/results/rfdiffusion/erap2_medium_*.pdb",
    "/workspace/results/rfdiffusion/erap2_long_*.pdb",
    "/workspace/results/rfdiffusion/erap2_large_*.pdb",
]
OUTPUT_DIR = "/workspace/results/selectivity"
CONTACT_CUTOFF = 5.0

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

parser = PDBParser(QUIET=True)


def get_sequence(struct):
    """Extract 1-letter sequence and CA coords from structure."""
    seq, cas, resnums = [], [], []
    for r in struct.get_residues():
        if r.get_resname() in AA3TO1:
            seq.append(AA3TO1[r.get_resname()])
            resnums.append(r.get_id()[1])
            try:
                cas.append(r["CA"])
            except KeyError:
                cas.append(None)
    return "".join(seq), cas, resnums


def align_erap1_onto_erap2(e2_struct, e1_struct):
    """Structurally align ERAP1 onto ERAP2 using conserved CA atoms.

    This puts both proteins in the same coordinate frame so that
    binder contacts can be meaningfully compared.
    """
    from Bio import pairwise2

    e2_seq, e2_cas, e2_nums = get_sequence(e2_struct)
    e1_seq, e1_cas, e1_nums = get_sequence(e1_struct)

    # Align sequences
    alignments = pairwise2.align.globalxx(e2_seq, e1_seq, one_alignment_only=True)
    e2_aln, e1_aln = alignments[0].seqA, alignments[0].seqB

    # Find conserved (identical) positions with valid CA atoms
    e2_matched_cas = []
    e1_matched_cas = []
    e2_idx, e1_idx = 0, 0

    for i in range(len(e2_aln)):
        is_e2 = e2_aln[i] != "-"
        is_e1 = e1_aln[i] != "-"

        if is_e2 and is_e1 and e2_aln[i] == e1_aln[i]:
            # Conserved position
            if (e2_idx < len(e2_cas) and e1_idx < len(e1_cas)
                    and e2_cas[e2_idx] is not None and e1_cas[e1_idx] is not None):
                e2_matched_cas.append(e2_cas[e2_idx])
                e1_matched_cas.append(e1_cas[e1_idx])

        if is_e2:
            e2_idx += 1
        if is_e1:
            e1_idx += 1

    print("Aligning ERAP1 onto ERAP2 using %d conserved CA atoms..." % len(e2_matched_cas))

    # Superimpose ERAP1 onto ERAP2
    sup = Superimposer()
    sup.set_atoms(e2_matched_cas, e1_matched_cas)
    sup.apply(list(e1_struct.get_atoms()))

    print("  RMSD after alignment: %.2f A" % sup.rms)
    return sup.rms


# --- Load and align structures ---
print("Loading structures...")
erap2 = parser.get_structure("erap2", ERAP2_PDB)
erap1 = parser.get_structure("erap1", ERAP1_PDB)

rmsd = align_erap1_onto_erap2(erap2, erap1)

# Build neighbor search on aligned structures
erap2_atoms = list(erap2.get_atoms())
erap1_atoms = list(erap1.get_atoms())
erap2_ns = NeighborSearch(erap2_atoms)
erap1_ns = NeighborSearch(erap1_atoms)

# Load divergent residues
with open(DIVERGENT_JSON) as f:
    div_data = json.load(f)
divergent_set = set(div_data["surface_divergent"])
print("Loaded %d divergent surface residues\n" % len(divergent_set))


def get_binder_atoms(pdb_path):
    """Extract binder chain atoms."""
    struct = parser.get_structure("b", pdb_path)
    chains = list(struct.get_chains())
    if len(chains) < 2:
        return [], 0
    for chain in chains:
        if chain.id != "A":
            atoms = list(chain.get_atoms())
            n_res = len(set(a.get_parent().get_id()[1] for a in atoms))
            return atoms, n_res
    return [], 0


def count_contacts(binder_atoms, target_ns, residue_filter=None):
    """Count contacts with optional residue filter."""
    contact_residues = set()
    n_contacts = 0
    for a in binder_atoms:
        neighbors = target_ns.search(a.get_vector().get_array(), CONTACT_CUTOFF)
        for n in neighbors:
            resnum = n.get_parent().get_id()[1]
            if residue_filter is None or resnum in residue_filter:
                contact_residues.add(resnum)
                n_contacts += 1
    return n_contacts, contact_residues


def score_binder(pdb_path):
    """Score a binder for ERAP2 selectivity vs aligned ERAP1."""
    name = os.path.basename(pdb_path)
    atoms, n_res = get_binder_atoms(pdb_path)
    if not atoms:
        return None

    # ERAP2 contacts (all + divergent-only)
    e2_total, e2_res = count_contacts(atoms, erap2_ns)
    e2_div, e2_div_res = count_contacts(atoms, erap2_ns, divergent_set)

    # ERAP1 contacts (aligned — same coordinate frame)
    e1_total, e1_res = count_contacts(atoms, erap1_ns)

    # --- Selectivity metrics ---
    raw_ratio = e2_total / max(e1_total, 1)
    div_frac = e2_div / max(e2_total, 1)

    # Interface residue selectivity (how many ERAP2 contact residues are NOT in ERAP1)
    e2_only_res = e2_res - e1_res
    interface_selectivity = len(e2_only_res) / max(len(e2_res), 1)

    # Composite score
    composite = (
        0.25 * min(raw_ratio / 5.0, 1.0) +
        0.30 * div_frac +
        0.25 * interface_selectivity +
        0.10 * min(n_res / 50.0, 1.0) +
        0.10 * min(len(e2_div_res) / 15.0, 1.0)
    )

    # Campaign type
    if "sel_channel" in name:
        campaign = "divergent_channel"
    elif "sel_domainIV" in name:
        campaign = "divergent_domIV"
    elif "sel_cterm" in name:
        campaign = "divergent_cterm"
    else:
        campaign = "active_site_v1"

    # Tag
    if composite >= 0.60:
        tag = "HIGHLY_SELECTIVE"
    elif composite >= 0.45:
        tag = "SELECTIVE"
    elif composite >= 0.30:
        tag = "MODERATE"
    else:
        tag = "NON_SELECTIVE"

    return {
        "design": name,
        "campaign": campaign,
        "binder_residues": n_res,
        "e2_contacts": e2_total,
        "e2_contact_residues": len(e2_res),
        "e2_divergent_contacts": e2_div,
        "e2_divergent_residues": len(e2_div_res),
        "e1_contacts": e1_total,
        "e1_contact_residues": len(e1_res),
        "e2_only_residues": len(e2_only_res),
        "raw_ratio": round(raw_ratio, 2),
        "divergent_frac": round(div_frac, 3),
        "interface_selectivity": round(interface_selectivity, 3),
        "composite_score": round(composite, 4),
        "tag": tag,
        "divergent_residues_hit": sorted(e2_div_res),
    }


# --- Main ---
all_pdbs = []
for pattern in BINDER_GLOBS:
    all_pdbs.extend(sorted(glob.glob(pattern)))
all_pdbs = sorted(set(all_pdbs))
print("Screening %d binder designs...\n" % len(all_pdbs))

results = []
for pdb_path in all_pdbs:
    r = score_binder(pdb_path)
    if r:
        results.append(r)
        print("  %-35s [%-20s] comp=%.3f div=%.0f%% raw=%.1fx iface_sel=%.0f%% -> %s" % (
            r["design"], r["campaign"],
            r["composite_score"], r["divergent_frac"] * 100,
            r["raw_ratio"], r["interface_selectivity"] * 100, r["tag"]))

results.sort(key=lambda x: -x["composite_score"])

# Save
os.makedirs(OUTPUT_DIR, exist_ok=True)

with open(os.path.join(OUTPUT_DIR, "selectivity_v3.json"), "w") as f:
    json.dump(results, f, indent=2)

csv_fields = [
    "design", "campaign", "binder_residues",
    "e2_contacts", "e2_divergent_contacts", "e2_divergent_residues",
    "e1_contacts", "e2_only_residues",
    "raw_ratio", "divergent_frac", "interface_selectivity",
    "composite_score", "tag",
]
with open(os.path.join(OUTPUT_DIR, "selectivity_v3.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(results)

# --- Report ---
print("\n" + "=" * 90)
print("SELECTIVITY SCREEN V3 — ERAP2 vs ERAP1 (structurally aligned, RMSD=%.2fA)" % rmsd)
print("=" * 90)

campaigns = {}
for r in results:
    c = r["campaign"]
    if c not in campaigns:
        campaigns[c] = {"total": 0, "selective": 0, "highly": 0, "scores": []}
    campaigns[c]["total"] += 1
    campaigns[c]["scores"].append(r["composite_score"])
    if r["tag"] == "SELECTIVE":
        campaigns[c]["selective"] += 1
    elif r["tag"] == "HIGHLY_SELECTIVE":
        campaigns[c]["highly"] += 1

print("\nBy campaign:")
for c, s in sorted(campaigns.items(), key=lambda x: -max(x[1]["scores"])):
    avg = sum(s["scores"]) / len(s["scores"])
    best = max(s["scores"])
    print("  %-25s: %2d designs | %d highly + %d selective | avg=%.3f best=%.3f" % (
        c, s["total"], s["highly"], s["selective"], avg, best))

tags = {}
for r in results:
    tags[r["tag"]] = tags.get(r["tag"], 0) + 1
print("\nOverall: %s" % " | ".join("%s=%d" % (k, v) for k, v in sorted(tags.items())))

print("\nTOP 10 CANDIDATES:")
print("%-35s %-20s %6s %5s %5s %6s %s" % (
    "Design", "Campaign", "Score", "Div%", "Raw", "ISel%", "Tag"))
print("-" * 105)
for r in results[:10]:
    print("%-35s %-20s %6.3f %4.0f%% %4.1fx %5.0f%% %s" % (
        r["design"], r["campaign"],
        r["composite_score"], r["divergent_frac"] * 100,
        r["raw_ratio"], r["interface_selectivity"] * 100, r["tag"]))

if results:
    lead = results[0]
    print("\n*** LEAD: %s ***" % lead["design"])
    print("  Campaign: %s" % lead["campaign"])
    print("  Score: %.3f (%s)" % (lead["composite_score"], lead["tag"]))
    print("  ERAP2: %d contacts, %d residues (%d divergent)" % (
        lead["e2_contacts"], lead["e2_contact_residues"], lead["e2_divergent_residues"]))
    print("  ERAP1: %d contacts, %d residues" % (lead["e1_contacts"], lead["e1_contact_residues"]))
    print("  ERAP2-only interface residues: %d (%.0f%% selective)" % (
        lead["e2_only_residues"], lead["interface_selectivity"] * 100))
    print("  Divergent fraction: %.0f%%" % (lead["divergent_frac"] * 100))
    print("  Raw selectivity: %.1fx" % lead["raw_ratio"])

    # Compare campaigns
    print("\n--- CAMPAIGN COMPARISON ---")
    for c in ["divergent_channel", "divergent_cterm", "divergent_domIV", "active_site_v1"]:
        camp_results = [r for r in results if r["campaign"] == c]
        if not camp_results:
            continue
        avg_div = sum(r["divergent_frac"] for r in camp_results) / len(camp_results)
        avg_isel = sum(r["interface_selectivity"] for r in camp_results) / len(camp_results)
        avg_raw = sum(r["raw_ratio"] for r in camp_results) / len(camp_results)
        print("  %-25s: avg_div=%.0f%% avg_iface_sel=%.0f%% avg_raw=%.1fx" % (
            c, avg_div * 100, avg_isel * 100, avg_raw))

print("\nSaved: selectivity_v3.json, selectivity_v3.csv")
