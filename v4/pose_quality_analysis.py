"""Pose quality analysis: clustering consistency + H-bond/salt bridge counting.

Analyzes DiffPepDock output PDBs across ERAP2-K392, ERAP2-N392, ERAP1, IRAP.
Answers: is binding specific (convergent poses, directional H-bonds) or
non-specific (scattered poses, only van der Waals)?

Run on Vast.ai after docking completes.
"""
import os, json, glob, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np
from collections import defaultdict

DOCK_DIR = os.environ.get("DOCK_DIR", "/workspace/DiffPepBuilder/runs/docking")
# Find the actual timestamped subdirectory
subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
if subdirs:
    DOCK_DIR = subdirs[-1]
    print(f"Using dock dir: {DOCK_DIR}")


def parse_pdb(pdb_path):
    """Parse PDB into peptide (chain A) and receptor (chain B) atoms."""
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


def compute_hbonds(pep_atoms, rec_atoms, cutoff=3.5, angle_cutoff=None):
    """Count potential H-bonds between peptide and receptor.

    H-bond donors: N-H, O-H (backbone NH, Ser/Thr/Tyr OH, etc.)
    H-bond acceptors: O, N (backbone C=O, sidechain carboxyl, etc.)
    Distance cutoff: 3.5A between donor heavy atom and acceptor.
    """
    donor_atoms = {"N", "NE", "NE1", "NE2", "ND1", "ND2", "NH1", "NH2", "NZ", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2"}

    hbonds = 0

    # Peptide donors -> receptor acceptors
    for pa in pep_atoms:
        if pa["atom"] not in donor_atoms:
            continue
        for ra in rec_atoms:
            if ra["atom"] not in acceptor_atoms:
                continue
            dist = np.linalg.norm(pa["xyz"] - ra["xyz"])
            if dist < cutoff:
                hbonds += 1

    # Receptor donors -> peptide acceptors
    for ra in rec_atoms:
        if ra["atom"] not in donor_atoms:
            continue
        for pa in pep_atoms:
            if pa["atom"] not in acceptor_atoms:
                continue
            dist = np.linalg.norm(ra["xyz"] - pa["xyz"])
            if dist < cutoff:
                hbonds += 1

    return hbonds


def compute_salt_bridges(pep_atoms, rec_atoms, cutoff=4.0):
    """Count salt bridges (charged pairs within cutoff).

    Positive: LYS NZ, ARG NH1/NH2/NE, HIS ND1/NE2
    Negative: ASP OD1/OD2, GLU OE1/OE2
    """
    pos_atoms = {("LYS", "NZ"), ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE")}
    neg_atoms = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}

    bridges = 0

    # Peptide positive -> receptor negative
    for pa in pep_atoms:
        if (pa["res"], pa["atom"]) not in pos_atoms:
            continue
        for ra in rec_atoms:
            if (ra["res"], ra["atom"]) not in neg_atoms:
                continue
            if np.linalg.norm(pa["xyz"] - ra["xyz"]) < cutoff:
                bridges += 1

    # Peptide negative -> receptor positive
    for pa in pep_atoms:
        if (pa["res"], pa["atom"]) not in neg_atoms:
            continue
        for ra in rec_atoms:
            if (ra["res"], ra["atom"]) not in pos_atoms:
                continue
            if np.linalg.norm(pa["xyz"] - ra["xyz"]) < cutoff:
                bridges += 1

    return bridges


def compute_peptide_centroid(pep_atoms):
    """Get peptide center of mass (CA atoms only)."""
    ca_coords = [a["xyz"] for a in pep_atoms if a["atom"] == "CA"]
    if not ca_coords:
        return None
    return np.mean(ca_coords, axis=0)


def compute_peptide_ca_coords(pep_atoms):
    """Get ordered CA coordinates for RMSD calculation."""
    return np.array([a["xyz"] for a in pep_atoms if a["atom"] == "CA"])


def kabsch_rmsd(P, Q):
    """Compute RMSD after optimal superposition (Kabsch algorithm)."""
    if len(P) != len(Q) or len(P) == 0:
        return float("inf")

    # Center both
    p_center = P.mean(axis=0)
    q_center = Q.mean(axis=0)
    P_c = P - p_center
    Q_c = Q - q_center

    # Kabsch
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T

    P_rot = (R @ P_c.T).T
    diff = P_rot - Q_c
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


def analyze_target(target_dir, target_name):
    """Analyze all peptides docked to one target."""
    results = []

    for pep_dir in sorted(os.listdir(target_dir)):
        pep_path = os.path.join(target_dir, pep_dir)
        if not os.path.isdir(pep_path):
            continue

        pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))
        if not pdbs:
            continue

        # Analyze each sample
        centroids = []
        ca_coords_list = []
        hbond_counts = []
        salt_bridge_counts = []
        contact_counts = []

        for pdb in pdbs:
            try:
                pep, rec = parse_pdb(pdb)
                if not pep or not rec:
                    continue

                centroid = compute_peptide_centroid(pep)
                if centroid is not None:
                    centroids.append(centroid)

                ca = compute_peptide_ca_coords(pep)
                if len(ca) > 0:
                    ca_coords_list.append(ca)

                hb = compute_hbonds(pep, rec)
                hbond_counts.append(hb)

                sb = compute_salt_bridges(pep, rec)
                salt_bridge_counts.append(sb)

                # Quick contact count for comparison
                pep_xyz = np.array([a["xyz"] for a in pep])
                rec_xyz = np.array([a["xyz"] for a in rec])
                diff = pep_xyz[:, None, :] - rec_xyz[None, :, :]
                dists = np.sqrt(np.sum(diff**2, axis=-1))
                contact_counts.append(int(np.sum(dists < 4.5)))

            except Exception as e:
                continue

        if not centroids:
            continue

        # Pose clustering: centroid spread (lower = more convergent)
        centroids_arr = np.array(centroids)
        centroid_mean = centroids_arr.mean(axis=0)
        centroid_spread = np.sqrt(np.mean(np.sum((centroids_arr - centroid_mean)**2, axis=1)))

        # Pairwise RMSD between poses (sample up to 16 for speed)
        sample_cas = ca_coords_list[:16]
        pairwise_rmsds = []
        ref_len = min(len(c) for c in sample_cas) if sample_cas else 0
        if ref_len > 0:
            trimmed = [c[:ref_len] for c in sample_cas]
            for i in range(len(trimmed)):
                for j in range(i+1, len(trimmed)):
                    pairwise_rmsds.append(kabsch_rmsd(trimmed[i], trimmed[j]))

        avg_rmsd = np.mean(pairwise_rmsds) if pairwise_rmsds else 0

        result = {
            "target": target_name,
            "peptide": pep_dir,
            "n_samples": len(centroids),
            "centroid_spread_A": round(float(centroid_spread), 2),
            "avg_pairwise_rmsd_A": round(float(avg_rmsd), 2),
            "avg_hbonds": round(float(np.mean(hbond_counts)), 1),
            "std_hbonds": round(float(np.std(hbond_counts)), 1),
            "avg_salt_bridges": round(float(np.mean(salt_bridge_counts)), 1),
            "std_salt_bridges": round(float(np.std(salt_bridge_counts)), 1),
            "avg_contacts": round(float(np.mean(contact_counts)), 0),
            "hbond_per_contact": round(float(np.mean(hbond_counts)) / max(np.mean(contact_counts), 1) * 100, 1),
        }
        results.append(result)

        print(f"  {target_name}/{pep_dir}: spread={centroid_spread:.1f}A "
              f"RMSD={avg_rmsd:.1f}A hbonds={np.mean(hbond_counts):.1f} "
              f"salt={np.mean(salt_bridge_counts):.1f} contacts={np.mean(contact_counts):.0f}")

    return results


def main():
    print("=" * 80)
    print("POSE QUALITY ANALYSIS: Clustering + H-bonds + Salt Bridges")
    print("=" * 80)

    all_results = []

    for target in sorted(os.listdir(DOCK_DIR)):
        target_path = os.path.join(DOCK_DIR, target)
        if not os.path.isdir(target_path):
            continue
        print(f"\n--- {target} ---")
        results = analyze_target(target_path, target)
        all_results.extend(results)

    # Save raw results
    os.makedirs("/workspace/results", exist_ok=True)
    with open("/workspace/results/pose_quality.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Print comparison table
    print("\n" + "=" * 120)
    print(f"{'Target':>12s}  {'Peptide':>20s}  {'Spread(A)':>10s}  {'RMSD(A)':>8s}  "
          f"{'H-bonds':>8s}  {'Salt-br':>8s}  {'Contacts':>9s}  {'Hb/Cont%':>9s}")
    print("-" * 120)

    for r in all_results:
        print(f"{r['target']:>12s}  {r['peptide']:>20s}  {r['centroid_spread_A']:>10.1f}  "
              f"{r['avg_pairwise_rmsd_A']:>8.1f}  {r['avg_hbonds']:>8.1f}  "
              f"{r['avg_salt_bridges']:>8.1f}  {r['avg_contacts']:>9.0f}  "
              f"{r['hbond_per_contact']:>9.1f}")

    # Summary by target
    print("\n" + "=" * 80)
    print("SUMMARY BY TARGET")
    print("-" * 80)

    by_target = defaultdict(list)
    for r in all_results:
        by_target[r["target"]].append(r)

    for target in sorted(by_target.keys()):
        rs = by_target[target]
        avg_spread = np.mean([r["centroid_spread_A"] for r in rs])
        avg_rmsd = np.mean([r["avg_pairwise_rmsd_A"] for r in rs])
        avg_hb = np.mean([r["avg_hbonds"] for r in rs])
        avg_salt = np.mean([r["avg_salt_bridges"] for r in rs])
        avg_cont = np.mean([r["avg_contacts"] for r in rs])
        avg_hb_pct = np.mean([r["hbond_per_contact"] for r in rs])

        # Variance of metrics across peptides (does target discriminate?)
        spread_cv = np.std([r["centroid_spread_A"] for r in rs]) / max(avg_spread, 0.01) * 100
        hb_cv = np.std([r["avg_hbonds"] for r in rs]) / max(avg_hb, 0.01) * 100

        print(f"{target:>12s}: spread={avg_spread:.1f}A  RMSD={avg_rmsd:.1f}A  "
              f"hbonds={avg_hb:.1f}  salt={avg_salt:.1f}  contacts={avg_cont:.0f}  "
              f"hb%={avg_hb_pct:.1f}  peptide_CV={hb_cv:.0f}%")

    print("\nInterpretation:")
    print("  Low spread + low RMSD = convergent poses (specific binding)")
    print("  High spread + high RMSD = scattered poses (non-specific)")
    print("  High hbond% = directional/specific interactions")
    print("  Low hbond% = van der Waals only (non-specific)")
    print("  High peptide_CV = target discriminates between peptides (good)")
    print("  Low peptide_CV = target treats all peptides same (non-specific)")

    print(f"\nSaved: /workspace/results/pose_quality.json ({len(all_results)} entries)")


if __name__ == "__main__":
    main()
