"""Select top hits from contact scoring and run MD retention test.

Uses subprocess pattern: PyRosetta refine → OpenMM CUDA MD.
Measures residence time: fraction of frames where P1 atom is within
4.0A of the nearest receptor heavy atom (proxy for catalytic zinc distance).

Goal: ERAP2-N392 residence > 80%, IRAP residence < 20%.
"""
import os, sys, json, glob, subprocess
try:
    sys.stdout.reconfigure(encoding="utf-8")
except:
    pass
import numpy as np

DOCK_DIR = "/workspace/DiffPepBuilder/runs/docking"
subdirs = sorted(glob.glob(os.path.join(DOCK_DIR, "*D_*M_*Y_*")))
if subdirs:
    DOCK_DIR = subdirs[-1]


def find_best_pose(target, peptide):
    pep_path = os.path.join(DOCK_DIR, target, peptide)
    pdbs = sorted(glob.glob(os.path.join(pep_path, "*.pdb")))[:16]
    best_pdb, best_c = None, 0
    for pdb in pdbs:
        pep, rec = [], []
        with open(pdb) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if line[21] == "A":
                    pep.append(xyz)
                elif line[21] == "B":
                    rec.append(xyz)
        if not pep or not rec:
            continue
        diff = np.array(pep)[:, None, :] - np.array(rec)[None, :, :]
        c = int(np.sum(np.sqrt(np.sum(diff**2, axis=-1)) < 4.5))
        if c > best_c:
            best_c, best_pdb = c, pdb
    return best_pdb, best_c


def fix_gly_cb(input_path, output_path):
    lines = []
    with open(input_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[17:20].strip() == "GLY" and line[12:16].strip() == "CB":
                continue
            lines.append(line)
    with open(output_path, "w") as f:
        f.writelines(lines)


def run_refine(input_pdb, output_pdb):
    """Run PyRosetta refine as subprocess."""
    script = f'''
import pyrosetta
pyrosetta.init("-mute all -ignore_unrecognized_res")
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta import get_fa_scorefxn
pose = pyrosetta.pose_from_pdb("{input_pdb}")
sfxn = get_fa_scorefxn()
relax = FastRelax()
relax.set_scorefxn(sfxn)
relax.max_iter(100)
relax.apply(pose)
pose.dump_pdb("{output_pdb}")
print("Refined: {output_pdb}")
'''
    result = subprocess.run(
        ["python3", "-c", script],
        capture_output=True, text=True, timeout=300
    )
    if result.returncode != 0:
        print(f"  Refine error: {result.stderr[-200:]}")
        return False
    print(f"  {result.stdout.strip()}")
    return True


def run_md(input_pdb, label, output_json, n_steps=2500000):
    """Run OpenMM MD as subprocess. 5ns = 2,500,000 steps at 2fs."""
    script = f'''
import os, json
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds, Simulation
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit

fixer = PDBFixer(filename="{input_pdb}")
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
modeller = Modeller(fixer.topology, fixer.positions)
modeller.addSolvent(ff, padding=0.8*unit.nanometers)
print(f"Atoms: {{modeller.topology.getNumAtoms()}}")

system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
system.addForce(MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))

platform = Platform.getPlatformByName("CUDA")
sim = Simulation(modeller.topology, system, integrator, platform, {{"CudaPrecision": "mixed"}})
sim.context.setPositions(modeller.positions)
sim.minimizeEnergy(maxIterations=1000)
sim.step(50000)  # 100ps equilibration

chains = list(modeller.topology.chains())
pep_atoms = [a.index for a in chains[0].atoms()]
pep_ca = [a.index for a in chains[0].atoms() if a.name == "CA"]
rec_atoms = [a.index for a in chains[1].atoms()] if len(chains) > 1 else []

# P1 atoms (first residue of peptide)
p1_atoms = [a.index for a in chains[0].residues().__next__().atoms()]

state = sim.context.getState(getPositions=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
init_pep = pos[pep_ca].copy()

interval = 25000  # 50ps
n_snaps = {n_steps} // interval
rmsds, dists, residence = [], [], []

for i in range(n_snaps):
    sim.step(interval)
    state = sim.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)

    cur_pep = pos[pep_ca]
    rmsd = float(np.sqrt(np.mean(np.sum((cur_pep - init_pep)**2, axis=1))))

    # Centroid distance
    if rec_atoms:
        dist = float(np.linalg.norm(np.mean(cur_pep, 0) - np.mean(pos[[a for a in rec_atoms if a < len(pos)]], 0)))
    else:
        dist = 0.0

    # Residence: min distance from any P1 atom to any receptor atom < 4.0A
    p1_pos = pos[p1_atoms]
    rec_pos = pos[rec_atoms] if rec_atoms else pos[pep_atoms]
    min_p1_dist = float(np.min(np.sqrt(np.sum((p1_pos[:, None, :] - rec_pos[None, :, :])**2, axis=-1))))
    in_contact = min_p1_dist < 4.0

    rmsds.append(round(rmsd, 2))
    dists.append(round(dist, 2))
    residence.append(1 if in_contact else 0)

    t = (i+1) * interval * 0.002 / 1000
    if (i+1) % 20 == 0:
        res_pct = sum(residence) / len(residence) * 100
        print(f"  t={{t:.1f}}ns RMSD={{rmsd:.1f}}A dist={{dist:.1f}}A res={{res_pct:.0f}}%")

res_total = sum(residence) / len(residence) * 100
result = {{
    "label": "{label}",
    "n_steps": {n_steps},
    "final_rmsd": rmsds[-1],
    "avg_rmsd": round(float(np.mean(rmsds[-20:])), 2),
    "rmsd_drift": round(rmsds[-1] - rmsds[0], 2),
    "residence_pct": round(res_total, 1),
    "residence_last_2ns": round(sum(residence[-40:]) / max(len(residence[-40:]), 1) * 100, 1),
}}
with open("{output_json}", "w") as f:
    json.dump(result, f, indent=2)
print(f"Residence: {{res_total:.1f}}%")
print(f"Saved: {output_json}")
'''
    result = subprocess.run(
        ["python3", "-c", script],
        capture_output=True, text=True, timeout=1800
    )
    if result.returncode != 0:
        print(f"  MD error: {result.stderr[-300:]}")
        return None
    print(f"  {result.stdout.strip()}")
    return output_json


def main():
    print("=" * 70)
    print("SHORT-STACK MD: Select top hits → Refine → 5ns MD → Residence Time")
    print("=" * 70)

    # Load contact scores to find best peptides
    scores_path = "/workspace/results/short_stack_contacts.json"
    if not os.path.exists(scores_path):
        print("ERROR: Run contact_score.py first")
        return

    with open(scores_path) as f:
        scores = json.load(f)

    # Rank by N392 contacts, filter for ERAP1 selectivity
    from collections import defaultdict
    by_pep = defaultdict(dict)
    for r in scores:
        by_pep[r["peptide"]][r["target"]] = r

    candidates = []
    for pep, targets in by_pep.items():
        n392 = targets.get("erap2_n392", {}).get("avg_contacts", 0)
        e1 = targets.get("erap1", {}).get("avg_contacts", 0)
        ir = targets.get("irap", {}).get("avg_contacts", 0)
        selectivity = n392 / max(e1, 1)
        irap_ratio = ir / max(n392, 1)
        candidates.append({
            "peptide": pep,
            "n392_contacts": n392,
            "erap1_contacts": e1,
            "irap_contacts": ir,
            "selectivity": round(selectivity, 2),
            "irap_ratio": round(irap_ratio, 2),
        })

    # Sort by N392 contacts, take top 3 with best ERAP1 selectivity
    candidates.sort(key=lambda c: (-c["selectivity"], -c["n392_contacts"]))
    top3 = candidates[:3]

    print("\nTop 3 candidates for MD:")
    for c in top3:
        print(f"  {c['peptide']}: N392={c['n392_contacts']:.0f} ERAP1={c['erap1_contacts']:.0f} "
              f"IRAP={c['irap_contacts']:.0f} sel={c['selectivity']}x irap_ratio={c['irap_ratio']}x")

    os.makedirs("/workspace/results", exist_ok=True)
    all_md_results = []

    for candidate in top3:
        pep = candidate["peptide"]
        print(f"\n{'='*60}")
        print(f"Processing: {pep}")
        print(f"{'='*60}")

        for target in ["erap2_n392", "irap"]:
            print(f"\n  --- {target} ---")
            pose_pdb, contacts = find_best_pose(target, pep)
            if not pose_pdb:
                print(f"  No pose found for {target}/{pep}")
                continue

            print(f"  Best pose: {contacts} contacts")

            # Fix GLY CB
            fixed = f"/workspace/results/{pep}_{target}_fixed.pdb"
            fix_gly_cb(pose_pdb, fixed)

            # Refine (subprocess — PyRosetta)
            refined = f"/workspace/results/{pep}_{target}_refined.pdb"
            print(f"  Refining with PyRosetta...")
            ok = run_refine(fixed, refined)
            if not ok or not os.path.exists(refined):
                print(f"  Refine failed, using fixed PDB")
                refined = fixed

            # MD (subprocess — OpenMM CUDA)
            md_json = f"/workspace/results/md_{pep}_{target}.json"
            print(f"  Running 5ns MD...")
            run_md(refined, f"{pep}_{target}", md_json, n_steps=2500000)

            if os.path.exists(md_json):
                with open(md_json) as f:
                    md_result = json.load(f)
                md_result["peptide"] = pep
                md_result["target"] = target
                all_md_results.append(md_result)

    # Final comparison
    print("\n" + "=" * 70)
    print("RESIDENCE TIME COMPARISON")
    print("-" * 70)
    print(f"{'Peptide':>15s}  {'Target':>12s}  {'Residence%':>10s}  {'Last 2ns%':>10s}  {'Avg RMSD':>9s}")
    print("-" * 70)

    for r in all_md_results:
        print(f"{r['peptide']:>15s}  {r['target']:>12s}  {r['residence_pct']:>10.1f}  "
              f"{r.get('residence_last_2ns', 0):>10.1f}  {r['avg_rmsd']:>9.2f}")

    # Save combined
    with open("/workspace/results/short_stack_md_final.json", "w") as f:
        json.dump(all_md_results, f, indent=2)
    print(f"\nSaved: /workspace/results/short_stack_md_final.json")

    # Check if any meet criteria
    print("\n" + "=" * 70)
    print("PASS/FAIL CRITERIA: N392 residence > 80%, IRAP residence < 20%")
    print("-" * 70)
    from collections import defaultdict
    by_pep_md = defaultdict(dict)
    for r in all_md_results:
        by_pep_md[r["peptide"]][r["target"]] = r

    for pep, targets in by_pep_md.items():
        n392_res = targets.get("erap2_n392", {}).get("residence_pct", 0)
        irap_res = targets.get("irap", {}).get("residence_pct", 0)
        n392_pass = "PASS" if n392_res > 80 else "FAIL"
        irap_pass = "PASS" if irap_res < 20 else "FAIL"
        overall = "CANDIDATE" if n392_pass == "PASS" and irap_pass == "PASS" else "REJECT"
        print(f"  {pep}: N392={n392_res:.0f}% ({n392_pass}) IRAP={irap_res:.0f}% ({irap_pass}) → {overall}")


if __name__ == "__main__":
    main()
