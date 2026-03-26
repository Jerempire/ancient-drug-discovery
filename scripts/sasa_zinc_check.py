"""
SASA Zinc Capping Test — V4.3 Wet Lab Validation

Tests whether 6-mer peptides functionally CAP the catalytic zinc site in ERAP2.
AlphaFold structures lack explicit Zn2+, so we measure SASA of the three
zinc-coordinating residues (H370, H374, E393) with and without peptide bound.

If peptide reduces their combined SASA by >80%, it's a functional plug.

Uses BioPython ShrakeRupley SASA calculator.
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
from pathlib import Path
from datetime import datetime, timezone

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

# --- Config ---
PROJECT = Path(r"C:\Users\jmj2z\Projects\medical\ancient-drug-discovery")
STARTING_STRUCTS = PROJECT / "data" / "results" / "v43_validation" / "md_starting_structures"
OUTPUT_DIR = PROJECT / "data" / "results" / "v43_validation" / "terminal_a" / "sasa"
APO_STRUCTURE = PROJECT / "data" / "structures" / "erap2_wt_alphafold.pdb"

PEPTIDES = ["VAGSAF", "IAFSAF", "VAWSAF", "FASGAV"]
ZINC_RESIDUES = [370, 374, 393]  # H370, H374, E393 — zinc-coordinating in HEXXH motif
PROTEIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"


class ProteinOnlySelect(Select):
    """Select only chain A (protein) atoms for apo SASA calculation."""
    def accept_chain(self, chain):
        return chain.id == PROTEIN_CHAIN


def get_zinc_residue_sasa(structure, model_idx=0):
    """Compute per-atom SASA, then sum for zinc-coordinating residues."""
    sr = ShrakeRupley()
    model = structure[model_idx]
    sr.compute(model, level="A")  # atom-level SASA

    residue_sasa = {}
    for chain in model:
        if chain.id != PROTEIN_CHAIN:
            continue
        for residue in chain:
            resid = residue.id[1]
            if resid in ZINC_RESIDUES:
                atom_sasa = sum(atom.sasa for atom in residue if atom.sasa is not None)
                residue_sasa[resid] = {
                    "resname": residue.resname,
                    "resid": resid,
                    "sasa_A2": round(atom_sasa, 2),
                    "n_atoms": len(list(residue.get_atoms()))
                }
    return residue_sasa


def compute_apo_sasa(parser):
    """Compute SASA for apo ERAP2 (no peptide) using the AlphaFold structure."""
    apo = parser.get_structure("apo", str(APO_STRUCTURE))
    return get_zinc_residue_sasa(apo)


def compute_complex_sasa(parser, pdb_path):
    """Compute SASA for the full complex (protein + peptide)."""
    complex_struct = parser.get_structure("complex", str(pdb_path))
    return get_zinc_residue_sasa(complex_struct)


def compute_protein_only_sasa(parser, pdb_path):
    """
    Compute SASA for protein chain only (peptide removed) from the complex PDB.
    This gives the 'holo-minus-peptide' reference — same protein conformation
    as the complex but without the peptide atoms blocking solvent.
    """
    complex_struct = parser.get_structure("prot_only", str(pdb_path))
    model = complex_struct[0]

    # Remove peptide chain
    chains_to_remove = [c for c in model if c.id != PROTEIN_CHAIN]
    for chain in chains_to_remove:
        model.detach_child(chain.id)

    sr = ShrakeRupley()
    sr.compute(model, level="A")

    residue_sasa = {}
    for chain in model:
        for residue in chain:
            resid = residue.id[1]
            if resid in ZINC_RESIDUES:
                atom_sasa = sum(atom.sasa for atom in residue if atom.sasa is not None)
                residue_sasa[resid] = {
                    "resname": residue.resname,
                    "resid": resid,
                    "sasa_A2": round(atom_sasa, 2)
                }
    return residue_sasa


def classify(reduction_pct):
    if reduction_pct >= 80:
        return "FUNCTIONAL_PLUG"
    elif reduction_pct >= 50:
        return "PARTIAL_PLUG"
    else:
        return "PASSENGER"


def peptide_zinc_distance(parser, pdb_path):
    """Compute minimum distance from any peptide heavy atom to zinc-coordinating residue atoms."""
    import numpy as np
    struct = parser.get_structure("dist", str(pdb_path))
    model = struct[0]

    zinc_atoms = []
    peptide_atoms = []

    for chain in model:
        if chain.id == PROTEIN_CHAIN:
            for res in chain:
                if res.id[1] in ZINC_RESIDUES:
                    for atom in res:
                        if atom.element != "H":
                            zinc_atoms.append(atom.get_vector().get_array())
        elif chain.id == PEPTIDE_CHAIN:
            for res in chain:
                for atom in res:
                    if atom.element != "H":
                        peptide_atoms.append(atom.get_vector().get_array())

    if not zinc_atoms or not peptide_atoms:
        return None

    zinc_coords = np.array(zinc_atoms)
    pep_coords = np.array(peptide_atoms)

    # Pairwise distance matrix
    diff = zinc_coords[:, None, :] - pep_coords[None, :, :]
    dists = np.sqrt(np.sum(diff**2, axis=-1))

    return round(float(np.min(dists)), 2)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)

    # 1. Compute apo SASA (reference)
    print("Computing apo ERAP2 SASA (AlphaFold structure)...")
    apo_sasa = compute_apo_sasa(parser)
    apo_total = sum(v["sasa_A2"] for v in apo_sasa.values())
    print(f"  Apo zinc-site SASA: {apo_total:.1f} A^2")
    for resid, data in sorted(apo_sasa.items()):
        print(f"    {data['resname']}{resid}: {data['sasa_A2']:.1f} A^2")

    results = {}
    summary_rows = []

    for peptide in PEPTIDES:
        pdb_file = STARTING_STRUCTS / f"{peptide}_vs_erap2k392.pdb"
        if not pdb_file.exists():
            print(f"\nWARNING: {pdb_file} not found, skipping {peptide}")
            continue

        print(f"\n--- {peptide} vs ERAP2-K392 ---")

        # Complex SASA (with peptide)
        complex_sasa = compute_complex_sasa(parser, pdb_file)
        complex_total = sum(v["sasa_A2"] for v in complex_sasa.values())

        # Protein-only SASA (peptide removed, same protein conformation)
        prot_only_sasa = compute_protein_only_sasa(parser, pdb_file)
        prot_only_total = sum(v["sasa_A2"] for v in prot_only_sasa.values())

        # The meaningful comparison: protein-only vs complex
        # (how much does the peptide occlude the zinc site?)
        reduction = prot_only_total - complex_total
        reduction_pct = (reduction / prot_only_total * 100) if prot_only_total > 0 else 0

        # Also compare to apo
        apo_vs_complex_reduction = apo_total - complex_total
        apo_vs_complex_pct = (apo_vs_complex_reduction / apo_total * 100) if apo_total > 0 else 0

        # Peptide-to-zinc-site distance
        min_dist = peptide_zinc_distance(parser, pdb_file)

        verdict = classify(reduction_pct)

        print(f"  Complex (with peptide):   {complex_total:.1f} A^2")
        print(f"  Protein-only (no peptide): {prot_only_total:.1f} A^2")
        print(f"  Apo (AlphaFold):           {apo_total:.1f} A^2")
        print(f"  Peptide occlusion:         {reduction:.1f} A^2 ({reduction_pct:.1f}%)")
        print(f"  vs Apo reduction:          {apo_vs_complex_reduction:.1f} A^2 ({apo_vs_complex_pct:.1f}%)")
        print(f"  Min peptide-zinc dist:     {min_dist} A")
        print(f"  Verdict:                   {verdict}")

        result = {
            "peptide": peptide,
            "target": "erap2k392",
            "apo_zinc_sasa_A2": round(apo_total, 2),
            "complex_zinc_sasa_A2": round(complex_total, 2),
            "protein_only_zinc_sasa_A2": round(prot_only_total, 2),
            "peptide_occlusion_A2": round(reduction, 2),
            "peptide_occlusion_pct": round(reduction_pct, 1),
            "apo_vs_complex_reduction_A2": round(apo_vs_complex_reduction, 2),
            "apo_vs_complex_reduction_pct": round(apo_vs_complex_pct, 1),
            "min_peptide_zinc_distance_A": min_dist,
            "per_residue": {
                str(resid): {
                    "resname": complex_sasa[resid]["resname"],
                    "apo_sasa": apo_sasa[resid]["sasa_A2"],
                    "complex_sasa": complex_sasa[resid]["sasa_A2"],
                    "prot_only_sasa": prot_only_sasa[resid]["sasa_A2"]
                }
                for resid in ZINC_RESIDUES if resid in complex_sasa
            },
            "verdict": verdict
        }

        results[peptide] = result

        # Write individual result
        out_file = OUTPUT_DIR / f"{peptide}_sasa.json"
        with open(out_file, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  Saved: {out_file.name}")

        summary_rows.append({
            "peptide": peptide,
            "zinc_sasa_complex": round(complex_total, 1),
            "zinc_sasa_prot_only": round(prot_only_total, 1),
            "occlusion_pct": round(reduction_pct, 1),
            "min_zinc_dist_A": min_dist,
            "verdict": verdict
        })

    # Write summary
    summary = {
        "test": "SASA Zinc Capping",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "apo_reference": {
            "source": str(APO_STRUCTURE.name),
            "zinc_residues": {str(k): v for k, v in apo_sasa.items()},
            "total_zinc_sasa_A2": round(apo_total, 2)
        },
        "thresholds": {
            "functional_plug": ">= 80% occlusion",
            "partial_plug": "50-80% occlusion",
            "passenger": "< 50% occlusion"
        },
        "results": summary_rows,
        "overall_verdict": "PASS" if any(r["verdict"] == "FUNCTIONAL_PLUG" for r in summary_rows) else "REVIEW"
    }

    summary_file = OUTPUT_DIR / "sasa_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"{'Peptide':<10} {'Occlusion%':>10} {'Zn Dist':>8} {'Verdict':<20}")
    print(f"{'-'*50}")
    for row in summary_rows:
        print(f"{row['peptide']:<10} {row['occlusion_pct']:>9.1f}% {row['min_zinc_dist_A']:>7.1f} {row['verdict']:<20}")
    print(f"\nSummary saved: {summary_file.name}")


if __name__ == "__main__":
    main()
