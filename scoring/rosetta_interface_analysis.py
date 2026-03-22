"""PyRosetta-based interface analysis for ERAP2 binder candidates.

Adds physics-based structural validation on top of Boltz-2 ML predictions.
Takes complex CIF/PDB files and runs: FastRelax, InterfaceAnalyzer,
alanine scanning, and per-residue energy breakdown.

Requires PyRosetta (free academic license from graylab.jhu.edu).

Usage:
    python -m scoring.rosetta_interface_analysis
    # or from project root:
    python scoring/rosetta_interface_analysis.py
"""

import glob
import json
import os
import sqlite3
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# --- Config ---
CONFIG = {
    "interface": "A_B",           # Rosetta chain designation: target=A, binder=B
    "relax_cycles": 1,            # 1 for retrospective, 3 for V3 forward
    "ala_scan_cutoff": 8.0,       # Angstroms for interface residue detection
    "ala_scan_chain": "B",        # Scan binder residues only by default
    "output_dir": str(PROJECT_ROOT / "data" / "results" / "rosetta"),
    "db_path": str(PROJECT_ROOT / "data" / "results" / "candidates.db"),
    "complexes_dir": str(PROJECT_ROOT / "data" / "results" / "boltz2_complexes"),
}

# Divergent hotspot residues (ERAP2 numbering, chain A)
HOTSPOT_RESIDUES = {353, 355, 360, 363, 367, 401, 403, 406, 408, 412}

# Interpretation thresholds (REU = Rosetta Energy Units)
THRESHOLDS = {
    "dG_separated": {"weak": -5.0, "adequate": -15.0},   # more negative = stronger
    "dSASA_int": {"weak": 500.0, "adequate": 800.0},     # Angstrom^2
    "packstat": {"weak": 0.50, "adequate": 0.65},         # 0-1
    "shape_complementarity": {"weak": 0.55, "adequate": 0.65},  # 0-1
    "delta_unsatHbonds": {"concern": 4},                  # more = worse
}


# ================================================================
# CIF/PDB Loading
# ================================================================

def find_cif_file(output_dir: str) -> Optional[str]:
    """Locate a CIF file in a directory tree."""
    cifs = glob.glob(os.path.join(output_dir, "**", "*.cif"), recursive=True)
    return cifs[0] if cifs else None


def _cif_to_pdb_gemmi(cif_path: str) -> str:
    """Convert CIF to PDB using gemmi. Returns path to temp PDB file."""
    import gemmi
    doc = gemmi.cif.read(cif_path)
    st = gemmi.make_structure_from_block(doc[0])
    # Fix occupancy=0.0 (Boltz-2 artifact)
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.occ == 0.0:
                        atom.occ = 1.0
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(tmp.name)
    return tmp.name


def _cif_to_pdb_biopython(cif_path: str) -> str:
    """Fallback CIF->PDB conversion using BioPython."""
    from Bio.PDB import MMCIFParser, PDBIO

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", cif_path)
    # Fix occupancy
    for atom in structure.get_atoms():
        if atom.get_occupancy() == 0.0:
            atom.set_occupancy(1.0)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp.name)
    return tmp.name


def load_complex_pose(path: str):
    """Load a CIF or PDB file into a PyRosetta Pose.

    Tries pose_from_mmcif first (PyRosetta native), falls back to
    gemmi conversion, then BioPython conversion.
    """
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.core.import_pose import pose_from_mmcif

    ext = os.path.splitext(path)[1].lower()

    if ext == ".pdb":
        return pose_from_pdb(path)

    # CIF: try native PyRosetta first
    try:
        pose = pose_from_mmcif(path)
        if pose.total_residue() > 0:
            return pose
    except Exception:
        pass

    # Fallback: gemmi
    try:
        pdb_path = _cif_to_pdb_gemmi(path)
        pose = pose_from_pdb(pdb_path)
        os.unlink(pdb_path)
        return pose
    except Exception:
        pass

    # Fallback: BioPython
    pdb_path = _cif_to_pdb_biopython(path)
    pose = pose_from_pdb(pdb_path)
    os.unlink(pdb_path)
    return pose


# ================================================================
# FastRelax
# ================================================================

def run_fast_relax(pose, n_cycles: int = 1):
    """Energy-minimize a complex with FastRelax (REF2015).

    Returns (relaxed_pose, score_before, score_after).
    """
    from pyrosetta import Pose
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.relax import FastRelax

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    score_before = sfxn(pose)

    relaxed = Pose(pose)  # deep copy
    relax = FastRelax()
    relax.set_scorefxn(sfxn)
    relax.max_iter(n_cycles * 200)  # default is 200 per cycle
    relax.apply(relaxed)

    score_after = sfxn(relaxed)
    print(f"  FastRelax: {score_before:.1f} -> {score_after:.1f} REU "
          f"(delta={score_after - score_before:.1f})")
    return relaxed, score_before, score_after


# ================================================================
# InterfaceAnalyzer
# ================================================================

def run_interface_analyzer(pose, interface: str = "A_B") -> dict:
    """Run InterfaceAnalyzerMover and return all metrics.

    Returns dict with: dG_separated, dSASA_int, packstat,
    shape_complementarity, hbonds_int, delta_unsatHbonds, etc.
    """
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    sfxn(pose)  # score first

    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    iam.set_scorefunction(sfxn)
    iam.set_pack_separated(True)
    iam.set_compute_interface_sc(True)
    iam.set_compute_packstat(True)
    iam.apply(pose)

    metrics = {
        "dG_separated": iam.get_interface_dG(),
        "dSASA_int": iam.get_interface_delta_sasa(),
        "packstat": iam.get_interface_packstat(),
        "shape_complementarity": iam.get_interface_sc(),
        "hbonds_int": iam.get_interface_delta_hbond_unsat(),
        "nres_int": iam.get_num_interface_residues(),
        "per_residue_energy_int": iam.get_interface_dG() / max(iam.get_num_interface_residues(), 1),
    }

    # delta_unsatHbonds: get from the pose's cached data
    try:
        metrics["delta_unsatHbonds"] = iam.get_interface_delta_hbond_unsat()
    except Exception:
        metrics["delta_unsatHbonds"] = None

    print(f"  InterfaceAnalyzer: dG={metrics['dG_separated']:.1f} REU, "
          f"BSA={metrics['dSASA_int']:.0f} A^2, "
          f"packstat={metrics['packstat']:.3f}, "
          f"Sc={metrics['shape_complementarity']:.3f}")
    return metrics


# ================================================================
# Alanine Scanning
# ================================================================

def run_alanine_scan(pose, interface: str = "A_B",
                     scan_chain: str = "B",
                     cutoff: float = 8.0) -> list[dict]:
    """Identify interface hotspot residues by alanine scanning.

    For each interface residue on scan_chain: mutate to Ala, repack,
    compute ddG. Positive ddG = residue is important for binding.

    Returns list of {resnum, resname, chain, ddG, is_hotspot}.
    """
    from pyrosetta import Pose
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking,
        PreventRepacking,
        OperateOnResidueSubset,
    )
    from pyrosetta.rosetta.core.select.residue_selector import (
        ChainSelector,
        NeighborhoodResidueSelector,
        ResidueIndexSelector,
    )
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")

    def calc_binding_energy(p):
        """Compute binding energy by separation method."""
        from pyrosetta import Pose
        bound = Pose(p)
        sfxn(bound)
        bound_score = sfxn(bound)

        # Separate chains by large translation
        unbound = Pose(p)
        jump_id = 1
        trans = unbound.jump(jump_id)
        from pyrosetta.rosetta.numeric import xyzVector_double_t
        shift = xyzVector_double_t(500.0, 0.0, 0.0)
        trans.set_translation(shift)
        unbound.set_jump(jump_id, trans)

        # Repack after separation
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        packer = PackRotamersMover(sfxn)
        packer.task_factory(tf)
        packer.apply(unbound)
        unbound_score = sfxn(unbound)

        return bound_score - unbound_score

    # Get wild-type binding energy
    wt_dG = calc_binding_energy(pose)

    # Find interface residues on scan_chain
    results = []
    chain_id = ord(scan_chain) - ord('A') + 1  # 1-indexed

    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if pdb_info and pdb_info.chain(i) != scan_chain:
            continue

        resname = pose.residue(i).name3().strip()

        # Skip Gly (already minimal), Pro (backbone constraints), Ala (no-op)
        if resname in ("GLY", "PRO", "ALA"):
            continue

        # Check if this residue is at the interface (any atom within cutoff of other chain)
        is_interface = False
        for j in range(1, pose.total_residue() + 1):
            if pdb_info and pdb_info.chain(j) == scan_chain:
                continue
            # Check CA-CA distance as quick proxy
            try:
                ca_i = pose.residue(i).xyz("CA")
                ca_j = pose.residue(j).xyz("CA")
                dist = ca_i.distance(ca_j)
                if dist < cutoff + 4.0:  # generous cutoff for CA-based check
                    is_interface = True
                    break
            except Exception:
                continue

        if not is_interface:
            continue

        # Mutate to Ala
        mutant = Pose(pose)
        mutate = MutateResidue(i, "ALA")
        mutate.apply(mutant)

        # Repack around mutation site
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        packer = PackRotamersMover(sfxn)
        packer.task_factory(tf)
        packer.apply(mutant)

        mut_dG = calc_binding_energy(mutant)
        ddG = mut_dG - wt_dG  # positive = destabilizing = important residue

        resnum = pdb_info.number(i) if pdb_info else i

        results.append({
            "resnum": resnum,
            "resname": resname,
            "chain": scan_chain,
            "ddG": round(ddG, 2),
            "is_hotspot": ddG > 1.0,
        })

    hotspots = [r for r in results if r["is_hotspot"]]
    print(f"  Alanine scan: {len(results)} interface residues scanned, "
          f"{len(hotspots)} hotspots (ddG > 1.0 REU)")
    return results


# ================================================================
# Residue Energy Breakdown
# ================================================================

def run_residue_energy_breakdown(pose) -> dict:
    """Per-residue total energy decomposition.

    Returns {pose_resnum: {chain, resnum, resname, total_energy}}.
    Useful for identifying which residues drive binding vs cross-reactivity.
    """
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    sfxn(pose)

    breakdown = {}
    pdb_info = pose.pdb_info()

    for i in range(1, pose.total_residue() + 1):
        energies = pose.energies()
        total_e = energies.residue_total_energy(i)

        chain = pdb_info.chain(i) if pdb_info else "?"
        resnum = pdb_info.number(i) if pdb_info else i
        resname = pose.residue(i).name3().strip()

        breakdown[i] = {
            "chain": chain,
            "resnum": resnum,
            "resname": resname,
            "total_energy": round(total_e, 3),
        }

    return breakdown


# ================================================================
# flex_ddg (V3 — variant discrimination)
# ================================================================

def run_flex_ddg(pose, position: int, from_aa: str, to_aa: str,
                 chain: str = "A") -> float:
    """Predict ddG for a point mutation at the interface.

    Used for V3 K392N variant discrimination: does the binder
    prefer K392 or N392?

    Returns ddG in REU (positive = mutation destabilizes binding).
    """
    from pyrosetta import Pose
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.relax import FastRelax

    sfxn = ScoreFunctionFactory.create_score_function("ref2015")

    # Find the pose residue number for this PDB position + chain
    pdb_info = pose.pdb_info()
    target_pose_num = None
    for i in range(1, pose.total_residue() + 1):
        if pdb_info and pdb_info.chain(i) == chain and pdb_info.number(i) == position:
            target_pose_num = i
            break

    if target_pose_num is None:
        print(f"  WARNING: Residue {chain}{position} not found in pose")
        return float("nan")

    actual_res = pose.residue(target_pose_num).name1()
    if actual_res != from_aa:
        print(f"  WARNING: Expected {from_aa} at {chain}{position}, found {actual_res}")

    # Wild-type: relax + score
    wt = Pose(pose)
    relax_wt = FastRelax()
    relax_wt.set_scorefxn(sfxn)
    relax_wt.max_iter(200)
    relax_wt.apply(wt)
    wt_score = sfxn(wt)

    # Mutant: mutate + relax + score
    mut = Pose(pose)
    # Convert 1-letter to 3-letter for MutateResidue
    aa_map = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    mutate = MutateResidue(target_pose_num, aa_map.get(to_aa, to_aa))
    mutate.apply(mut)

    relax_mut = FastRelax()
    relax_mut.set_scorefxn(sfxn)
    relax_mut.max_iter(200)
    relax_mut.apply(mut)
    mut_score = sfxn(mut)

    ddG = mut_score - wt_score
    print(f"  flex_ddg {chain}{position}{from_aa}->{to_aa}: "
          f"ddG = {ddG:.2f} REU (wt={wt_score:.1f}, mut={mut_score:.1f})")
    return round(ddG, 2)


# ================================================================
# Interpretation / Flagging
# ================================================================

def interpret_metrics(metrics: dict) -> list[str]:
    """Generate warning flags based on threshold interpretation."""
    flags = []

    dG = metrics.get("dG_separated")
    if dG is not None:
        if dG > 0:
            flags.append(f"CRITICAL: dG_separated={dG:.1f} (positive = unfavorable binding)")
        elif dG > THRESHOLDS["dG_separated"]["weak"]:
            flags.append(f"WEAK: dG_separated={dG:.1f} (> {THRESHOLDS['dG_separated']['weak']})")

    bsa = metrics.get("dSASA_int")
    if bsa is not None and bsa < THRESHOLDS["dSASA_int"]["weak"]:
        flags.append(f"WEAK: BSA={bsa:.0f} A^2 (< {THRESHOLDS['dSASA_int']['weak']})")

    ps = metrics.get("packstat")
    if ps is not None and ps < THRESHOLDS["packstat"]["weak"]:
        flags.append(f"WEAK: packstat={ps:.3f} (< {THRESHOLDS['packstat']['weak']})")

    unsat = metrics.get("delta_unsatHbonds")
    if unsat is not None and unsat > THRESHOLDS["delta_unsatHbonds"]["concern"]:
        flags.append(f"CONCERN: {unsat} buried unsatisfied H-bonds (> {THRESHOLDS['delta_unsatHbonds']['concern']})")

    return flags


# ================================================================
# Database Export
# ================================================================

ROSETTA_SCHEMA = """
CREATE TABLE IF NOT EXISTS rosetta_metrics (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id TEXT NOT NULL,
    target_protein TEXT NOT NULL DEFAULT 'ERAP2',
    analyzed_at TEXT NOT NULL,
    source_structure TEXT,
    relax_cycles INTEGER,
    score_pre_relax REAL,
    score_post_relax REAL,
    dG_separated REAL,
    dSASA_int REAL,
    packstat REAL,
    shape_complementarity REAL,
    hbonds_int INTEGER,
    delta_unsatHbonds INTEGER,
    nres_int INTEGER,
    per_residue_energy_int REAL,
    n_hotspot_residues INTEGER,
    hotspot_residues TEXT,
    energy_breakdown_json TEXT,
    flex_ddg_K392N REAL,
    rosetta_flags TEXT,
    UNIQUE(candidate_id, target_protein)
);

CREATE INDEX IF NOT EXISTS idx_rosetta_candidate
    ON rosetta_metrics(candidate_id);
"""


def init_rosetta_db(db_path: str) -> sqlite3.Connection:
    """Create rosetta_metrics table if it doesn't exist."""
    conn = sqlite3.connect(db_path)
    conn.executescript(ROSETTA_SCHEMA)
    return conn


def save_to_db(results: dict, candidate_id: str, target_protein: str,
               source_path: str, db_path: str) -> None:
    """Insert or update rosetta metrics in the database."""
    conn = init_rosetta_db(db_path)

    hotspots = [r for r in results.get("ala_scan", []) if r["is_hotspot"]]

    conn.execute(
        """INSERT OR REPLACE INTO rosetta_metrics (
            candidate_id, target_protein, analyzed_at, source_structure,
            relax_cycles, score_pre_relax, score_post_relax,
            dG_separated, dSASA_int, packstat, shape_complementarity,
            hbonds_int, delta_unsatHbonds, nres_int, per_residue_energy_int,
            n_hotspot_residues, hotspot_residues,
            energy_breakdown_json, flex_ddg_K392N, rosetta_flags
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        (
            candidate_id,
            target_protein,
            datetime.now().isoformat(),
            source_path,
            results.get("relax_cycles"),
            results.get("score_pre_relax"),
            results.get("score_post_relax"),
            results["interface"].get("dG_separated"),
            results["interface"].get("dSASA_int"),
            results["interface"].get("packstat"),
            results["interface"].get("shape_complementarity"),
            results["interface"].get("hbonds_int"),
            results["interface"].get("delta_unsatHbonds"),
            results["interface"].get("nres_int"),
            results["interface"].get("per_residue_energy_int"),
            len(hotspots),
            json.dumps(hotspots),
            json.dumps(results.get("energy_breakdown")),
            results.get("flex_ddg_K392N"),
            json.dumps(results.get("flags", [])),
        )
    )
    conn.commit()
    conn.close()


# ================================================================
# Export to JSON
# ================================================================

def export_results(results: dict, candidate_id: str, target_protein: str,
                   out_dir: str) -> None:
    """Write analysis results to JSON file."""
    os.makedirs(out_dir, exist_ok=True)
    filename = f"rosetta_{candidate_id}_{target_protein}.json"
    out_path = os.path.join(out_dir, filename)

    # Make JSON-serializable
    export = {
        "candidate_id": candidate_id,
        "target_protein": target_protein,
        "analyzed_at": datetime.now().isoformat(),
        "relax": {
            "cycles": results.get("relax_cycles"),
            "score_before": results.get("score_pre_relax"),
            "score_after": results.get("score_post_relax"),
        },
        "interface": results.get("interface", {}),
        "ala_scan": results.get("ala_scan", []),
        "flags": results.get("flags", []),
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(export, f, indent=2)
    print(f"  Wrote {out_path}")


# ================================================================
# Main Analysis Pipeline
# ================================================================

def analyze_complex(cif_path: str, candidate_id: str,
                    target_protein: str = "ERAP2",
                    relax_cycles: int = 1,
                    run_ala_scan: bool = True,
                    run_energy_breakdown: bool = True,
                    flex_ddg_mutation: Optional[tuple] = None) -> dict:
    """Run the full Rosetta analysis battery on a single complex.

    Args:
        cif_path: Path to CIF or PDB file of the complex
        candidate_id: Identifier (e.g. "n248_trim_c5")
        target_protein: Which target (ERAP2, IRAP, etc.)
        relax_cycles: Number of FastRelax cycles (1=fast, 3=thorough)
        run_ala_scan: Whether to run alanine scanning
        run_energy_breakdown: Whether to run per-residue energy decomposition
        flex_ddg_mutation: Optional (position, from_aa, to_aa, chain) for flex_ddg

    Returns:
        dict with all results
    """
    print(f"\n{'='*60}")
    print(f"Analyzing: {candidate_id} vs {target_protein}")
    print(f"Source: {cif_path}")
    print(f"{'='*60}")

    results = {}

    # Load structure
    print("Loading structure...")
    pose = load_complex_pose(cif_path)
    print(f"  Loaded: {pose.total_residue()} residues, "
          f"{pose.num_chains()} chains")

    # FastRelax
    print(f"Running FastRelax ({relax_cycles} cycle(s))...")
    pose, score_before, score_after = run_fast_relax(pose, relax_cycles)
    results["relax_cycles"] = relax_cycles
    results["score_pre_relax"] = round(score_before, 2)
    results["score_post_relax"] = round(score_after, 2)

    # InterfaceAnalyzer
    print("Running InterfaceAnalyzer...")
    interface_metrics = run_interface_analyzer(pose, CONFIG["interface"])
    results["interface"] = interface_metrics

    # Alanine scanning
    if run_ala_scan:
        print("Running alanine scan...")
        ala_results = run_alanine_scan(
            pose, CONFIG["interface"],
            CONFIG["ala_scan_chain"], CONFIG["ala_scan_cutoff"]
        )
        results["ala_scan"] = ala_results

    # Residue energy breakdown
    if run_energy_breakdown:
        print("Running residue energy breakdown...")
        breakdown = run_residue_energy_breakdown(pose)
        results["energy_breakdown"] = breakdown

    # flex_ddg (V3 only)
    if flex_ddg_mutation:
        pos, from_aa, to_aa, chain = flex_ddg_mutation
        print(f"Running flex_ddg: {chain}{pos}{from_aa}->{to_aa}...")
        ddG = run_flex_ddg(pose, pos, from_aa, to_aa, chain)
        results["flex_ddg_K392N"] = ddG

    # Interpretation flags
    results["flags"] = interpret_metrics(interface_metrics)
    if results["flags"]:
        print("  FLAGS:")
        for f in results["flags"]:
            print(f"    - {f}")

    return results


def print_summary_table(all_results: list[tuple[str, str, dict]]) -> None:
    """Print a comparison table of all analyzed complexes."""
    print("\n" + "=" * 100)
    print("ROSETTA INTERFACE ANALYSIS SUMMARY")
    print("=" * 100)
    print(f"{'Candidate':<24} {'Target':<8} {'dG_sep':>8} {'BSA':>8} "
          f"{'Packstat':>9} {'Sc':>6} {'Hbonds':>7} {'Unsat':>6} "
          f"{'Hotspots':>9} {'Flags'}")
    print("-" * 100)

    for candidate_id, target, results in all_results:
        iface = results.get("interface", {})
        hotspots = [r for r in results.get("ala_scan", []) if r["is_hotspot"]]
        flags = results.get("flags", [])
        flag_str = "; ".join(f[:20] for f in flags[:2]) if flags else ""

        print(f"{candidate_id:<24} {target:<8} "
              f"{iface.get('dG_separated', 0):>8.1f} "
              f"{iface.get('dSASA_int', 0):>8.0f} "
              f"{iface.get('packstat', 0):>9.3f} "
              f"{iface.get('shape_complementarity', 0):>6.3f} "
              f"{iface.get('hbonds_int', 0):>7} "
              f"{iface.get('delta_unsatHbonds', 0):>6} "
              f"{len(hotspots):>9} "
              f"{flag_str}")

    print("-" * 100)


# ================================================================
# Batch Runner
# ================================================================

# Candidates to analyze — update paths after Vast.ai CIF acquisition
CANDIDATES = [
    # (candidate_id, target_protein, cif_subdir_or_path)
    ("n248_trim_c5", "ERAP2", "n248_trim_c5_erap2"),
    ("n248_trim_c5_Y4A", "ERAP2", "n248_trim_c5_y4a_erap2"),
    ("n248_wt", "ERAP2", "n248_wt_erap2"),
    ("n248_ko_all_aromatics", "ERAP2", "n248_ko_aromatics_erap2"),
    ("n248_trim_c5", "IRAP", "n248_trim_c5_irap"),
    ("n248_trim_c5_Y4A", "IRAP", "n248_trim_c5_y4a_irap"),
]


def main():
    from pyrosetta import init
    init("-mute all -ignore_unrecognized_res -use_input_sc")

    out_dir = CONFIG["output_dir"]
    db_path = CONFIG["db_path"]
    complexes_dir = CONFIG["complexes_dir"]
    relax_cycles = CONFIG["relax_cycles"]

    os.makedirs(out_dir, exist_ok=True)

    all_results = []

    for candidate_id, target, cif_subdir in CANDIDATES:
        # Locate CIF file
        cif_dir = os.path.join(complexes_dir, cif_subdir)
        cif_path = find_cif_file(cif_dir)

        if not cif_path:
            # Try direct path
            direct = os.path.join(complexes_dir, f"{cif_subdir}.cif")
            if os.path.exists(direct):
                cif_path = direct
            else:
                print(f"\nSKIPPING {candidate_id} vs {target}: "
                      f"no CIF found in {cif_dir}")
                continue

        results = analyze_complex(
            cif_path=cif_path,
            candidate_id=candidate_id,
            target_protein=target,
            relax_cycles=relax_cycles,
            run_ala_scan=True,
            run_energy_breakdown=(target == "ERAP2"),  # skip for cross-reactivity targets
        )

        # Export
        export_results(results, candidate_id, target, out_dir)
        save_to_db(results, candidate_id, target, cif_path, db_path)

        all_results.append((candidate_id, target, results))

    # Summary
    if all_results:
        print_summary_table(all_results)

        # IRAP cross-reactivity verdict
        erap2_results = {c: r for c, t, r in all_results if t == "ERAP2"}
        irap_results = {c: r for c, t, r in all_results if t == "IRAP"}

        if irap_results:
            print("\n" + "=" * 60)
            print("IRAP CROSS-REACTIVITY VERDICT")
            print("=" * 60)
            for cand, irap_r in irap_results.items():
                irap_dG = irap_r["interface"].get("dG_separated", 0)
                erap2_dG = erap2_results.get(cand, {}).get("interface", {}).get("dG_separated", 0)
                print(f"  {cand}:")
                print(f"    ERAP2 dG = {erap2_dG:.1f} REU")
                print(f"    IRAP  dG = {irap_dG:.1f} REU")
                if irap_dG < -15:
                    print(f"    VERDICT: REAL CONCERN (IRAP dG < -15)")
                elif irap_dG < -5:
                    print(f"    VERDICT: MARGINAL (IRAP dG between -5 and -15)")
                else:
                    print(f"    VERDICT: LIKELY ARTIFACT (IRAP dG > -5)")

    print(f"\nDone. Results in {out_dir}")


if __name__ == "__main__":
    main()
