"""Step 1: PyRosetta — Fix DiffPepDock PDB and relax sidechains.

Takes raw DiffPepDock output, removes GLY CB artifacts, adds hydrogens,
runs FastRelax, saves clean PDB for OpenMM.

Usage: python refine.py <input.pdb> <output.pdb>
"""
import sys
import os

def fix_gly_cb(input_path, output_path):
    """Remove CB atoms from GLY residues (DiffPepDock artifact)."""
    lines = []
    with open(input_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                if resname == "GLY" and atom == "CB":
                    continue
            lines.append(line)
    with open(output_path, "w") as f:
        f.writelines(lines)
    return output_path


def refine(input_pdb, output_pdb):
    """Fix artifacts and relax with PyRosetta."""
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res")

    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta import get_fa_scorefxn

    # Fix GLY CB first
    temp_pdb = input_pdb + ".fixed.pdb"
    fix_gly_cb(input_pdb, temp_pdb)

    # Load into PyRosetta
    pose = pyrosetta.pose_from_pdb(temp_pdb)
    print(f"  Loaded: {pose.total_residue()} residues, {pose.num_chains()} chains")

    # Quick relax (1 cycle, not full 5)
    scorefxn = get_fa_scorefxn()
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(100)

    print("  Relaxing (1 cycle)...")
    relax.apply(pose)

    score = scorefxn(pose)
    print(f"  Score after relax: {score:.1f} REU")

    # Save
    pose.dump_pdb(output_pdb)
    print(f"  Saved: {output_pdb}")

    # Cleanup temp
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)

    return output_pdb


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python refine.py <input.pdb> <output.pdb>")
        sys.exit(1)
    refine(sys.argv[1], sys.argv[2])
