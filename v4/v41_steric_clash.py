"""V4.1 Steric Clash Design — Anti-IRAP selectivity via A406/N406 handle.

ERAP2 has Ala406 (tiny), IRAP has Asn406 (bulky).
Strategy: place bulky residue on peptide at position facing 406.
Since we don't know exact peptide-receptor register, scan P3-P6
with Trp (biggest natural AA) and Phe substitutions.

Run on Vast.ai after DiffPepDock docking.
"""
import os, json, glob, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np
from collections import defaultdict

DOCK_DIR = os.environ.get("DOCK_DIR", "/workspace/DiffPepBuilder/runs/docking")
subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
if subdirs:
    DOCK_DIR = subdirs[-1]


def parse_pdb(pdb_path):
    pep_atoms = []
    rec_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            entry = {
                "atom": atom_name, "res": resname, "chain": chain,
                "resnum": resnum, "xyz": np.array([x, y, z]),
            }
            if chain == "A":
                pep_atoms.append(entry)
            elif chain == "B":
                rec_atoms.append(entry)
    return pep_atoms, rec_atoms


def analyze_contacts_by_residue(pep_atoms, rec_atoms, cutoff=4.5):
    """For each peptide residue, count contacts with each receptor residue."""
    contacts = defaultdict(lambda: defaultdict(int))
    for pa in pep_atoms:
        for ra in rec_atoms:
            dist = np.linalg.norm(pa["xyz"] - ra["xyz"])
            if dist < cutoff:
                contacts[pa["resnum"]][ra["resnum"]] += 1
    return contacts


def compute_hbonds(pep_atoms, rec_atoms, cutoff=3.5):
    donor_atoms = {"N", "NE", "NE1", "NE2", "ND1", "ND2", "NH1", "NH2", "NZ", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2"}
    hbonds = 0
    for pa in pep_atoms:
        if pa["atom"] not in donor_atoms:
            continue
        for ra in rec_atoms:
            if ra["atom"] not in acceptor_atoms:
                continue
            if np.linalg.norm(pa["xyz"] - ra["xyz"]) < cutoff:
                hbonds += 1
    for ra in rec_atoms:
        if ra["atom"] not in donor_atoms:
            continue
        for pa in pep_atoms:
            if pa["atom"] not in acceptor_atoms:
                continue
            if np.linalg.norm(ra["xyz"] - pa["xyz"]) < cutoff:
                hbonds += 1
    return hbonds


def main():
    print("=" * 90)
    print("V4.1 STERIC CLASH ANALYSIS")
    print("=" * 90)

    all_results = []

    for target in sorted(os.listdir(DOCK_DIR)):
        target_path = os.path.join(DOCK_DIR, target)
        if not os.path.isdir(target_path):
            continue

        for pep_dir in sorted(os.listdir(target_path)):
            pep_path = os.path.join(target_path, pep_dir)
            if not os.path.isdir(pep_path):
                continue
            pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))

            all_contacts = []
            all_hbonds = []
            all_total_contacts = []

            for pdb in pdbs:
                try:
                    pep, rec = parse_pdb(pdb)
                    if not pep or not rec:
                        continue

                    contacts = analyze_contacts_by_residue(pep, rec)
                    all_contacts.append(contacts)

                    hb = compute_hbonds(pep, rec)
                    all_hbonds.append(hb)

                    pep_xyz = np.array([a["xyz"] for a in pep])
                    rec_xyz = np.array([a["xyz"] for a in rec])
                    diff = pep_xyz[:, None, :] - rec_xyz[None, :, :]
                    dists = np.sqrt(np.sum(diff**2, axis=-1))
                    all_total_contacts.append(int(np.sum(dists < 4.5)))
                except:
                    continue

            if not all_contacts:
                continue

            avg_hb = np.mean(all_hbonds)
            avg_contacts = np.mean(all_total_contacts)

            # Aggregate: which receptor residue does each peptide position contact most?
            pep_rec_map = defaultdict(lambda: defaultdict(list))
            for contacts in all_contacts:
                for pep_res, rec_contacts in contacts.items():
                    for rec_res, count in rec_contacts.items():
                        pep_rec_map[pep_res][rec_res].append(count)

            # For position 406: which peptide residues contact it?
            pos406_contacts = {}
            for pep_res in sorted(pep_rec_map.keys()):
                if 406 in pep_rec_map[pep_res]:
                    avg = np.mean(pep_rec_map[pep_res][406])
                    pos406_contacts[pep_res] = avg

            result = {
                "target": target,
                "peptide": pep_dir,
                "avg_contacts": round(float(avg_contacts), 0),
                "avg_hbonds": round(float(avg_hb), 1),
                "n_samples": len(all_contacts),
                "pos406_contacts": {str(k): round(v, 1) for k, v in pos406_contacts.items()},
            }
            all_results.append(result)

            p406_str = ", ".join(f"P{k}={v:.1f}" for k, v in sorted(pos406_contacts.items()))
            print(f"  {target}/{pep_dir}: contacts={avg_contacts:.0f} hb={avg_hb:.1f} | 406-contacts: {p406_str}")

    # Save
    os.makedirs("/workspace/results", exist_ok=True)
    with open("/workspace/results/v41_steric_clash.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary table
    print("\n" + "=" * 90)
    print(f"{'Target':>12s}  {'Peptide':>25s}  {'Contacts':>9s}  {'H-bonds':>8s}  {'IRAP vs ERAP2':>15s}")
    print("-" * 90)

    # Group by peptide for comparison
    by_peptide = defaultdict(dict)
    for r in all_results:
        by_peptide[r["peptide"]][r["target"]] = r

    for pep in sorted(by_peptide.keys()):
        targets = by_peptide[pep]
        for t in ["erap2_k392", "erap2_n392", "erap1", "irap"]:
            if t in targets:
                r = targets[t]
                k392_c = targets.get("erap2_k392", {}).get("avg_contacts", 0)
                ratio = ""
                if t == "irap" and k392_c > 0:
                    ratio = f"{r['avg_contacts']/k392_c:.2f}x"
                print(f"{t:>12s}  {pep:>25s}  {r['avg_contacts']:>9.0f}  {r['avg_hbonds']:>8.1f}  {ratio:>15s}")
        print()

    print(f"\nSaved: /workspace/results/v41_steric_clash.json")


if __name__ == "__main__":
    main()
