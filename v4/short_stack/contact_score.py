"""Contact scoring for short-stack peptides across all 4 targets."""
import os, json, glob, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np
from collections import defaultdict

DOCK_DIR = "/workspace/DiffPepBuilder/runs/docking"
subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
if subdirs:
    DOCK_DIR = subdirs[-1]

def score_pdb(pdb_path, cutoff=4.5):
    pep, rec = [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if line[21] == "A":
                pep.append(xyz)
            elif line[21] == "B":
                rec.append(xyz)
    if not pep or not rec:
        return None
    diff = np.array(pep)[:, None, :] - np.array(rec)[None, :, :]
    dists = np.sqrt(np.sum(diff**2, axis=-1))
    return int(np.sum(dists < cutoff))

def compute_hbonds(pdb_path, cutoff=3.5):
    donor_atoms = {"N", "NE", "NE1", "NE2", "ND1", "ND2", "NH1", "NH2", "NZ", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2"}
    pep, rec = [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            entry = {
                "atom": line[12:16].strip(), "res": line[17:20].strip(),
                "chain": line[21],
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            }
            if entry["chain"] == "A":
                pep.append(entry)
            elif entry["chain"] == "B":
                rec.append(entry)
    hb = 0
    for pa in pep:
        if pa["atom"] in donor_atoms:
            for ra in rec:
                if ra["atom"] in acceptor_atoms and np.linalg.norm(pa["xyz"] - ra["xyz"]) < cutoff:
                    hb += 1
    for ra in rec:
        if ra["atom"] in donor_atoms:
            for pa in pep:
                if pa["atom"] in acceptor_atoms and np.linalg.norm(ra["xyz"] - pa["xyz"]) < cutoff:
                    hb += 1
    return hb

results = []
for target in sorted(os.listdir(DOCK_DIR)):
    target_path = os.path.join(DOCK_DIR, target)
    if not os.path.isdir(target_path):
        continue
    for pep_dir in sorted(os.listdir(target_path)):
        pep_path = os.path.join(target_path, pep_dir)
        if not os.path.isdir(pep_path):
            continue
        pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))
        contacts_list, hbond_list = [], []
        for pdb in pdbs:
            c = score_pdb(pdb)
            if c is not None:
                contacts_list.append(c)
                hbond_list.append(compute_hbonds(pdb))
        if contacts_list:
            results.append({
                "target": target,
                "peptide": pep_dir,
                "avg_contacts": round(float(np.mean(contacts_list)), 1),
                "avg_hbonds": round(float(np.mean(hbond_list)), 1),
                "n_samples": len(contacts_list),
            })
            print(f"  {target}/{pep_dir}: contacts={np.mean(contacts_list):.0f} hb={np.mean(hbond_list):.1f}")

os.makedirs("/workspace/results", exist_ok=True)
with open("/workspace/results/short_stack_contacts.json", "w") as f:
    json.dump(results, f, indent=2)

# Print comparison table
print("\n" + "=" * 100)
print(f"{'Peptide':>15s}  {'N392 cont':>9s}  {'N392 hb':>7s}  {'K392 cont':>9s}  {'ERAP1 cont':>10s}  {'IRAP cont':>9s}  {'N392/IRAP':>9s}")
print("-" * 100)

by_pep = defaultdict(dict)
for r in results:
    by_pep[r["peptide"]][r["target"]] = r

for pep in sorted(by_pep.keys()):
    t = by_pep[pep]
    n392 = t.get("erap2_n392", {}).get("avg_contacts", 0)
    n392h = t.get("erap2_n392", {}).get("avg_hbonds", 0)
    k392 = t.get("erap2_k392", {}).get("avg_contacts", 0)
    e1 = t.get("erap1", {}).get("avg_contacts", 0)
    ir = t.get("irap", {}).get("avg_contacts", 0)
    ratio = f"{n392/ir:.2f}x" if ir > 0 else "n/a"
    print(f"{pep:>15s}  {n392:>9.0f}  {n392h:>7.1f}  {k392:>9.0f}  {e1:>10.0f}  {ir:>9.0f}  {ratio:>9s}")

print(f"\nSaved: /workspace/results/short_stack_contacts.json")
