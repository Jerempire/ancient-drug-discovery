"""Selectivity screen: score all RFdiffusion binders against ERAP2 vs ERAP1."""
import os, glob, json, csv
from Bio.PDB import PDBParser, NeighborSearch

ERAP2_HOTSPOTS = [370, 371, 374, 392, 393]
ERAP1_ACTIVE = [353, 354, 357, 376]

parser = PDBParser(QUIET=True)
erap2 = parser.get_structure("erap2", "/workspace/data/structures/erap2_wt_alphafold.pdb")
erap2_ns = NeighborSearch(list(erap2.get_atoms()))
erap1 = parser.get_structure("erap1", "/workspace/data/structures/erap1_wt_alphafold.pdb")
erap1_ns = NeighborSearch(list(erap1.get_atoms()))

binder_pdbs = sorted(glob.glob("/workspace/results/rfdiffusion/erap2_*.pdb"))
print("Screening %d binders...\n" % len(binder_pdbs))

results = []
for pdb_path in binder_pdbs:
    name = os.path.basename(pdb_path)
    struct = parser.get_structure(name, pdb_path)
    binder_atoms = []
    for chain in struct.get_chains():
        if chain.id == "B":
            binder_atoms = list(chain.get_atoms())
    if not binder_atoms:
        continue

    e2c = 0
    for a in binder_atoms:
        for n in erap2_ns.search(a.get_vector().get_array(), 5.0):
            if n.get_parent().get_id()[1] in ERAP2_HOTSPOTS:
                e2c += 1

    e1c = 0
    for a in binder_atoms:
        for n in erap1_ns.search(a.get_vector().get_array(), 5.0):
            if n.get_parent().get_id()[1] in ERAP1_ACTIVE:
                e1c += 1

    binder_res = len(set(a.get_parent().get_id()[1] for a in binder_atoms))
    sel = e2c / max(e1c, 1)
    tier = name.split("_")[1]
    tag = "SELECTIVE" if sel > 3 else ("MODERATE" if sel > 1.5 else "NON-SELECTIVE")

    results.append(dict(
        design=name, tier=tier, binder_residues=binder_res,
        erap2_contacts=e2c, erap1_contacts=e1c,
        selectivity_ratio=round(sel, 2), tag=tag,
    ))
    print("  %s: E2=%d E1=%d ratio=%.1f %s" % (name, e2c, e1c, sel, tag))

results.sort(key=lambda x: (-x["selectivity_ratio"], -x["erap2_contacts"]))

os.makedirs("/workspace/results/selectivity", exist_ok=True)
with open("/workspace/results/selectivity/selectivity_screen.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)

with open("/workspace/results/selectivity/selectivity_screen.json", "w") as f:
    json.dump(results, f, indent=2)

sel_list = [r for r in results if r["selectivity_ratio"] > 2.0]
ns = sum(1 for r in results if r["tag"] == "SELECTIVE")
nm = sum(1 for r in results if r["tag"] == "MODERATE")
nn = sum(1 for r in results if r["tag"] == "NON-SELECTIVE")

print("\n" + "=" * 70)
print("Total: %d | SELECTIVE(>3): %d | MODERATE: %d | NON-SEL: %d" % (len(results), ns, nm, nn))

if sel_list:
    print("\nTOP CANDIDATES:")
    for i, r in enumerate(sel_list[:5]):
        sz = r["binder_residues"]
        cost = "$400-800" if sz <= 40 else ("$800-1500" if sz <= 60 else "$500-1500 recombinant")
        d = r["design"]
        print("  %d. %s: %d res, sel=%.1fx, E2=%d E1=%d | %s" % (
            i + 1, d, sz, r["selectivity_ratio"], r["erap2_contacts"], r["erap1_contacts"], cost))
    lead = sel_list[0]
    print("\n*** LEAD: %s (%.1fx selective) ***" % (lead["design"], lead["selectivity_ratio"]))
else:
    print("No selective candidates found")
