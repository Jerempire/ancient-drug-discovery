"""Standalone PyRosetta analysis — reads existing Boltz-2 CIF files.

Usage (on Vast.ai after Boltz-2 predictions exist):
    python3 /workspace/scripts/run_rosetta_analysis_only.py
"""
import glob
import json
import os
import sys
import time

RESULTS_DIR = "/workspace/results/rosetta_pipeline"
BOLTZ_DIR = os.path.join(RESULTS_DIR, "boltz2")

COMPLEXES = [
    ("n248_trim_c5", "erap2"),
    ("n248_trim_c5", "irap"),
    ("n248_trim_c5_Y4A", "erap2"),
    ("n248_trim_c5_Y4A", "irap"),
    ("n248_wt", "erap2"),
    ("n248_ko_all_aromatics", "erap2"),
]


def find_cif(candidate, target):
    pattern = os.path.join(BOLTZ_DIR, f"{candidate}_vs_{target}", "**", "*.cif")
    cifs = glob.glob(pattern, recursive=True)
    return cifs[0] if cifs else None


def find_confidence_json(candidate, target):
    pattern = os.path.join(BOLTZ_DIR, f"{candidate}_vs_{target}", "**", "confidence_*.json")
    jsons = glob.glob(pattern, recursive=True)
    if jsons:
        with open(jsons[0]) as f:
            return json.load(f)
    return {}


def analyze(cif_path, candidate, target):
    from pyrosetta import Pose, pose_from_pdb
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.numeric import xyzVector_double_t

    # CIF -> PDB conversion
    import gemmi
    doc = gemmi.cif.read(cif_path)
    st = gemmi.make_structure_from_block(doc[0])
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.occ == 0.0:
                        atom.occ = 1.0
    pdb_path = cif_path.replace(".cif", "_conv.pdb")
    st.write_pdb(pdb_path)

    pose = pose_from_pdb(pdb_path)
    print(f"  Loaded: {pose.total_residue()} residues, {pose.num_chains()} chains")

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")

    # FastRelax
    print("  FastRelax...")
    score_before = sfxn(pose)
    relaxed = Pose(pose)
    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.max_iter(200)
    fr.apply(relaxed)
    score_after = sfxn(relaxed)
    print(f"    {score_before:.1f} -> {score_after:.1f} REU")

    # InterfaceAnalyzer
    print("  InterfaceAnalyzer...")
    sfxn(relaxed)
    iam = InterfaceAnalyzerMover()
    iam.set_interface("A_B")
    iam.set_scorefunction(sfxn)
    iam.set_pack_separated(True)
    iam.set_compute_packstat(True)
    iam.apply(relaxed)

    dG = iam.get_interface_dG()
    dSASA = iam.get_interface_delta_sasa()
    packstat = iam.get_interface_packstat()
    hbonds_unsat = iam.get_interface_delta_hbond_unsat()
    nres = iam.get_num_interface_residues()
    crossterm = iam.get_crossterm_interface_energy()
    hbond_E = iam.get_total_Hbond_E()

    print(f"    dG={dG:.1f}, BSA={dSASA:.0f}, pack={packstat:.3f}, "
          f"nres={nres}, unsat_hb={hbonds_unsat}")

    # Alanine scan (chain B)
    print("  Alanine scan...")

    def binding_energy(p):
        bound_score = sfxn(Pose(p))
        unbound = Pose(p)
        t = unbound.jump(1)
        t.set_translation(xyzVector_double_t(500.0, 0.0, 0.0))
        unbound.set_jump(1, t)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        pk = PackRotamersMover(sfxn)
        pk.task_factory(tf)
        pk.apply(unbound)
        return bound_score - sfxn(unbound)

    wt_dG = binding_energy(relaxed)
    ala_scan = []
    pdb_info = relaxed.pdb_info()

    for i in range(1, relaxed.total_residue() + 1):
        if pdb_info and pdb_info.chain(i) != "B":
            continue
        rn = relaxed.residue(i).name3().strip()
        if rn in ("GLY", "PRO", "ALA"):
            continue
        # Interface check
        is_iface = False
        try:
            ca_i = relaxed.residue(i).xyz("CA")
            for j in range(1, relaxed.total_residue() + 1):
                if pdb_info and pdb_info.chain(j) == "B":
                    continue
                if ca_i.distance(relaxed.residue(j).xyz("CA")) < 12.0:
                    is_iface = True
                    break
        except Exception:
            continue
        if not is_iface:
            continue

        mut = Pose(relaxed)
        MutateResidue(i, "ALA").apply(mut)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        pk = PackRotamersMover(sfxn)
        pk.task_factory(tf)
        pk.apply(mut)
        ddG = round(binding_energy(mut) - wt_dG, 2)
        ala_scan.append({
            "resnum": pdb_info.number(i) if pdb_info else i,
            "resname": rn, "chain": "B",
            "ddG": ddG, "is_hotspot": ddG > 1.0,
        })

    hotspots = [r for r in ala_scan if r["is_hotspot"]]
    print(f"    {len(ala_scan)} scanned, {len(hotspots)} hotspots")

    # Flags
    flags = []
    if dG > 0:
        flags.append(f"CRITICAL: dG={dG:.1f} (positive)")
    elif dG > -5:
        flags.append(f"WEAK: dG={dG:.1f}")
    if dSASA < 500:
        flags.append(f"WEAK: BSA={dSASA:.0f}")
    if packstat < 0.50:
        flags.append(f"WEAK: pack={packstat:.3f}")

    return {
        "candidate": candidate, "target": target,
        "score_before": round(score_before, 1),
        "score_after": round(score_after, 1),
        "dG_separated": round(dG, 2),
        "dSASA_int": round(dSASA, 1),
        "packstat": round(packstat, 4),
        "crossterm_energy": round(crossterm, 2),
        "hbond_energy": round(hbond_E, 2),
        "hbonds_unsat": hbonds_unsat,
        "nres_int": nres,
        "n_hotspots": len(hotspots),
        "hotspots": hotspots,
        "ala_scan": ala_scan,
        "flags": flags,
    }


def main():
    from pyrosetta import init
    init("-mute all -ignore_unrecognized_res -use_input_sc")

    t0 = time.time()
    all_results = []

    for candidate, target in COMPLEXES:
        cif = find_cif(candidate, target)
        if not cif:
            print(f"\nSKIP: {candidate} vs {target} — no CIF")
            continue

        conf = find_confidence_json(candidate, target)
        boltz_iptm = conf.get("iptm", "N/A")

        print(f"\n{'='*60}")
        print(f"  {candidate} vs {target} (Boltz ipTM={boltz_iptm})")
        print(f"  {cif}")

        try:
            r = analyze(cif, candidate, target)
            r["boltz_iptm"] = boltz_iptm
            all_results.append(r)

            # Save individual JSON
            out = os.path.join(RESULTS_DIR, f"rosetta_{candidate}_vs_{target}.json")
            with open(out, "w") as f:
                json.dump(r, f, indent=2)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()

    # Summary table
    print(f"\n\n{'='*110}")
    print("ROSETTA INTERFACE ANALYSIS — FINAL SUMMARY")
    print(f"{'='*110}")
    print(f"{'Candidate':<24} {'Target':<7} {'ipTM':>7} {'dG_sep':>8} "
          f"{'BSA':>7} {'Pack':>6} {'xterm':>7} {'HB_E':>7} "
          f"{'Unsat':>5} {'nRes':>5} {'Hot':>4} {'Flags'}")
    print("-" * 110)

    for r in all_results:
        iptm = r.get("boltz_iptm", "N/A")
        iptm_str = f"{iptm:.3f}" if isinstance(iptm, float) else str(iptm)
        flags_str = "; ".join(r.get("flags", [])[:2])
        print(f"{r['candidate']:<24} {r['target']:<7} {iptm_str:>7} "
              f"{r['dG_separated']:>8.1f} {r['dSASA_int']:>7.0f} "
              f"{r['packstat']:>6.3f} {r['crossterm_energy']:>7.1f} "
              f"{r['hbond_energy']:>7.1f} {r['hbonds_unsat']:>5} "
              f"{r['nres_int']:>5} {r['n_hotspots']:>4} {flags_str}")
    print("-" * 110)

    # IRAP verdict
    erap2 = {r["candidate"]: r for r in all_results if r["target"] == "erap2"}
    irap = {r["candidate"]: r for r in all_results if r["target"] == "irap"}

    if irap:
        print(f"\nIRAP CROSS-REACTIVITY VERDICT:")
        for c, ir in irap.items():
            e2 = erap2.get(c, {})
            e2_dG = e2.get("dG_separated", 0)
            ir_dG = ir["dG_separated"]
            print(f"  {c}:")
            print(f"    ERAP2 dG = {e2_dG:.1f} REU (ipTM={e2.get('boltz_iptm', 'N/A')})")
            print(f"    IRAP  dG = {ir_dG:.1f} REU (ipTM={ir.get('boltz_iptm', 'N/A')})")
            if ir_dG < -15:
                print("    -> REAL CONCERN")
            elif ir_dG < -5:
                print("    -> MARGINAL")
            else:
                print("    -> LIKELY ARTIFACT (weak physics binding)")

    # Save combined
    combined = os.path.join(RESULTS_DIR, "rosetta_combined_results.json")
    with open(combined, "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\nTotal time: {(time.time()-t0)/60:.1f} min")
    print(f"Results: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
