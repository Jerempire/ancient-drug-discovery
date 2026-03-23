"""Fast contact-based scoring of DiffPepDock results — no PyRosetta needed."""
import os, json, glob, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np

DOCK_DIR = "/workspace/DiffPepBuilder/runs/docking/23D_03M_2026Y_02h_15m"

def parse_pdb_contacts(pdb_path, cutoff=4.5):
    pep_atoms = []
    rec_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            chain = line[21]
            if chain == "A":
                pep_atoms.append((x, y, z))
            elif chain == "B":
                rec_atoms.append((x, y, z))

    if not pep_atoms or not rec_atoms:
        return None

    pep_xyz = np.array(pep_atoms)
    rec_xyz = np.array(rec_atoms)
    diff = pep_xyz[:, None, :] - rec_xyz[None, :, :]
    dists = np.sqrt(np.sum(diff**2, axis=-1))

    return {
        "n_contacts": int(np.sum(dists < cutoff)),
        "min_dist": float(np.min(dists)),
    }

results = []
for variant in ["erap2_k392", "erap2_n392"]:
    variant_dir = os.path.join(DOCK_DIR, variant)
    if not os.path.isdir(variant_dir):
        continue
    for pep_dir in sorted(os.listdir(variant_dir)):
        pep_path = os.path.join(variant_dir, pep_dir)
        if not os.path.isdir(pep_path):
            continue
        pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))

        contacts_list = []
        for pdb in pdbs:
            try:
                c = parse_pdb_contacts(pdb)
                if c:
                    contacts_list.append(c)
            except:
                continue

        if contacts_list:
            best = max(contacts_list, key=lambda c: c["n_contacts"])
            avg_contacts = np.mean([c["n_contacts"] for c in contacts_list])
            avg_min_dist = np.mean([c["min_dist"] for c in contacts_list])

            results.append({
                "peptide": pep_dir,
                "variant": variant.replace("erap2_", ""),
                "target": variant,
                "best_contacts": best["n_contacts"],
                "avg_contacts": round(float(avg_contacts), 1),
                "best_min_dist": round(best["min_dist"], 2),
                "avg_min_dist": round(float(avg_min_dist), 2),
                "n_samples": len(contacts_list),
                "score": round(-float(avg_contacts), 1),
            })
            print(f"  {variant}/{pep_dir}: contacts={avg_contacts:.0f} min_d={avg_min_dist:.1f} ({len(contacts_list)} samples)")

os.makedirs("/workspace/results", exist_ok=True)
with open("/workspace/results/diffpepdock_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nDone: {len(results)} results -> /workspace/results/diffpepdock_results.json")
