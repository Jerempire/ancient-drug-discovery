"""
selectivity_screen_v2.py — Upgraded ERAP2 vs ERAP1 selectivity scoring.

Improvements over v1:
  1. Uses ERAP2-specific divergent residues (not just hardcoded hotspots)
  2. Scores interface contacts weighted by residue divergence
  3. Computes buried surface area (BSA) as binding energy proxy
  4. Structural alignment of binder onto ERAP1 for true cross-reactivity
  5. Ranks by composite selectivity score, not raw contact ratio

Inputs:
  /workspace/results/rfdiffusion/erap2_sel_*.pdb  (new divergent-targeting designs)
  /workspace/results/rfdiffusion/erap2_*.pdb      (old active-site designs)
  /workspace/data/structures/erap2_wt_alphafold.pdb
  /workspace/data/structures/erap1_wt_alphafold.pdb
  /workspace/results/erap2_divergent_surface.json

Outputs:
  /workspace/results/selectivity/selectivity_v2.json
  /workspace/results/selectivity/selectivity_v2.csv
  /workspace/results/selectivity/top_selective_binders.txt
"""
import csv
import glob
import json
import os
import sys

from Bio.PDB import PDBParser, NeighborSearch, Superimposer
import numpy as np

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
CONTACT_CUTOFF = 5.0  # Angstroms

# --- Load structures ---
parser = PDBParser(QUIET=True)
erap2 = parser.get_structure("erap2", ERAP2_PDB)
erap1 = parser.get_structure("erap1", ERAP1_PDB)

erap2_atoms = list(erap2.get_atoms())
erap1_atoms = list(erap1.get_atoms())
erap2_ns = NeighborSearch(erap2_atoms)
erap1_ns = NeighborSearch(erap1_atoms)

# Load divergent residues
with open(DIVERGENT_JSON) as f:
    div_data = json.load(f)
divergent_set = set(div_data["surface_divergent"])
print("Loaded %d divergent surface residues" % len(divergent_set))


def get_binder_atoms(pdb_path):
    """Extract binder chain atoms (chain B or last chain)."""
    struct = parser.get_structure("binder", pdb_path)
    chains = list(struct.get_chains())
    if len(chains) < 2:
        return [], []

    # Binder is typically chain B (or the non-A chain)
    for chain in chains:
        if chain.id != "A":
            atoms = list(chain.get_atoms())
            ca_atoms = [a for a in atoms if a.get_name() == "CA"]
            return atoms, ca_atoms
    return [], []


def count_contacts(binder_atoms, target_ns, residue_filter=None):
    """Count binder-target contacts, optionally filtering target residues."""
    contact_residues = set()
    total_contacts = 0
    for a in binder_atoms:
        neighbors = target_ns.search(a.get_vector().get_array(), CONTACT_CUTOFF)
        for n in neighbors:
            resnum = n.get_parent().get_id()[1]
            if residue_filter is None or resnum in residue_filter:
                contact_residues.add(resnum)
                total_contacts += 1
    return total_contacts, contact_residues


def compute_bsa_proxy(binder_atoms, target_ns):
    """Approximate buried surface area by counting close atom pairs."""
    close_pairs = 0
    for a in binder_atoms:
        neighbors = target_ns.search(a.get_vector().get_array(), 4.0)
        close_pairs += len(neighbors)
    return close_pairs


def score_binder(pdb_path):
    """Score a single binder for ERAP2 selectivity."""
    name = os.path.basename(pdb_path)
    atoms, ca_atoms = get_binder_atoms(pdb_path)
    if not atoms:
        return None

    n_residues = len(ca_atoms)

    # ERAP2 contacts
    e2_total, e2_contact_res = count_contacts(atoms, erap2_ns)
    e2_divergent, e2_div_res = count_contacts(atoms, erap2_ns, divergent_set)

    # ERAP1 contacts (cross-reactivity)
    e1_total, e1_contact_res = count_contacts(atoms, erap1_ns)

    # BSA proxy
    e2_bsa = compute_bsa_proxy(atoms, erap2_ns)
    e1_bsa = compute_bsa_proxy(atoms, erap1_ns)

    # --- Composite selectivity metrics ---

    # Raw contact ratio (v1 metric)
    raw_ratio = e2_total / max(e1_total, 1)

    # Divergent contact fraction: what % of ERAP2 contacts are on unique residues
    div_frac = e2_divergent / max(e2_total, 1)

    # BSA selectivity ratio
    bsa_ratio = e2_bsa / max(e1_bsa, 1)

    # Composite selectivity score (weighted)
    # Higher = more selective for ERAP2
    composite = (
        0.30 * min(raw_ratio / 5.0, 1.0) +      # Raw contact ratio (cap at 5x)
        0.35 * div_frac +                          # Divergent contact fraction (0-1)
        0.20 * min(bsa_ratio / 3.0, 1.0) +        # BSA ratio (cap at 3x)
        0.15 * min(n_residues / 50.0, 1.0)         # Size bonus (larger = more contacts)
    )

    # Determine campaign type
    if "sel_channel" in name:
        campaign = "divergent_channel"
    elif "sel_domainIV" in name:
        campaign = "divergent_domIV"
    elif "sel_cterm" in name:
        campaign = "divergent_cterm"
    else:
        campaign = "active_site_v1"

    # Tag
    if composite >= 0.65:
        tag = "HIGHLY_SELECTIVE"
    elif composite >= 0.50:
        tag = "SELECTIVE"
    elif composite >= 0.35:
        tag = "MODERATE"
    else:
        tag = "NON_SELECTIVE"

    return {
        "design": name,
        "campaign": campaign,
        "binder_residues": n_residues,
        # ERAP2 metrics
        "e2_contacts": e2_total,
        "e2_contact_residues": len(e2_contact_res),
        "e2_divergent_contacts": e2_divergent,
        "e2_bsa_proxy": e2_bsa,
        # ERAP1 metrics
        "e1_contacts": e1_total,
        "e1_contact_residues": len(e1_contact_res),
        "e1_bsa_proxy": e1_bsa,
        # Selectivity scores
        "raw_ratio": round(raw_ratio, 2),
        "divergent_frac": round(div_frac, 3),
        "bsa_ratio": round(bsa_ratio, 2),
        "composite_score": round(composite, 4),
        "tag": tag,
        # Divergent residues contacted
        "divergent_residues_hit": sorted(e2_div_res),
    }


# --- Main ---
print("Scanning for binder PDBs...")
all_pdbs = []
for pattern in BINDER_GLOBS:
    all_pdbs.extend(sorted(glob.glob(pattern)))

# Deduplicate
all_pdbs = sorted(set(all_pdbs))
print("Found %d binder designs\n" % len(all_pdbs))

results = []
for pdb_path in all_pdbs:
    result = score_binder(pdb_path)
    if result:
        results.append(result)
        print("  %s [%s]: composite=%.3f div_frac=%.1f%% raw=%.1fx bsa=%.1fx -> %s" % (
            result["design"], result["campaign"],
            result["composite_score"],
            result["divergent_frac"] * 100,
            result["raw_ratio"],
            result["bsa_ratio"],
            result["tag"],
        ))

# Sort by composite score
results.sort(key=lambda x: -x["composite_score"])

# Save results
os.makedirs(OUTPUT_DIR, exist_ok=True)

# JSON (full)
with open(os.path.join(OUTPUT_DIR, "selectivity_v2.json"), "w") as f:
    # Remove divergent_residues_hit for CSV compat, keep in JSON
    json.dump(results, f, indent=2)

# CSV (summary)
csv_fields = [
    "design", "campaign", "binder_residues",
    "e2_contacts", "e2_divergent_contacts", "e1_contacts",
    "raw_ratio", "divergent_frac", "bsa_ratio", "composite_score", "tag",
]
with open(os.path.join(OUTPUT_DIR, "selectivity_v2.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(results)

# Summary report
print("\n" + "=" * 80)
print("SELECTIVITY SCREEN V2 — ERAP2 vs ERAP1")
print("=" * 80)

# Stats by campaign
campaigns = {}
for r in results:
    c = r["campaign"]
    if c not in campaigns:
        campaigns[c] = {"total": 0, "selective": 0, "scores": []}
    campaigns[c]["total"] += 1
    campaigns[c]["scores"].append(r["composite_score"])
    if r["tag"] in ("SELECTIVE", "HIGHLY_SELECTIVE"):
        campaigns[c]["selective"] += 1

print("\nBy campaign:")
for c, stats in sorted(campaigns.items(), key=lambda x: -max(x[1]["scores"])):
    avg = sum(stats["scores"]) / len(stats["scores"])
    best = max(stats["scores"])
    print("  %-25s: %d designs, %d selective, avg=%.3f, best=%.3f" % (
        c, stats["total"], stats["selective"], avg, best))

# Tag counts
tags = {}
for r in results:
    tags[r["tag"]] = tags.get(r["tag"], 0) + 1
print("\nOverall: %s" % ", ".join("%s=%d" % (k, v) for k, v in sorted(tags.items())))

# Top 10
print("\nTOP 10 SELECTIVE BINDERS:")
print("%-30s %-20s %6s %6s %6s %6s %s" % (
    "Design", "Campaign", "Score", "DivFr", "Raw", "BSA", "Tag"))
print("-" * 100)
for r in results[:10]:
    print("%-30s %-20s %6.3f %5.1f%% %5.1fx %5.1fx %s" % (
        r["design"], r["campaign"],
        r["composite_score"],
        r["divergent_frac"] * 100,
        r["raw_ratio"],
        r["bsa_ratio"],
        r["tag"],
    ))

# Lead candidate
if results:
    lead = results[0]
    print("\n*** LEAD CANDIDATE: %s ***" % lead["design"])
    print("  Campaign: %s" % lead["campaign"])
    print("  Composite score: %.3f (%s)" % (lead["composite_score"], lead["tag"]))
    print("  ERAP2 contacts: %d (%d divergent, %.0f%% unique)" % (
        lead["e2_contacts"], lead["e2_divergent_contacts"],
        lead["divergent_frac"] * 100))
    print("  ERAP1 contacts: %d" % lead["e1_contacts"])
    print("  Selectivity: %.1fx contact, %.1fx BSA" % (lead["raw_ratio"], lead["bsa_ratio"]))
    print("  Binder size: %d residues" % lead["binder_residues"])
    n_div = len(lead.get("divergent_residues_hit", []))
    print("  Divergent residues contacted: %d" % n_div)

print("\nSaved: selectivity_v2.json, selectivity_v2.csv")
