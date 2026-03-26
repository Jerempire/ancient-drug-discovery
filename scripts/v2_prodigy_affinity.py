"""
V2 Protein Binder — PRODIGY Binding Affinity Prediction

Runs PRODIGY (PROtein binDIng enerGY prediction) on V2 Boltz-2 CIF files
to estimate Kd for each enzyme target. Compares to known aminopeptidase
inhibitor affinities.

Usage: conda run -n ancient-drug-discovery python scripts/v2_prodigy_affinity.py
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import math
from datetime import datetime, timezone

from Bio.PDB import MMCIFParser
from prodigy_prot.cli import Prodigy

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF_DIR = os.path.join(PROJECT, "data", "results", "y87a_cif_files")
OUTPUT_DIR = os.path.join(PROJECT, "data", "results", "v2_inhibition")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# V2 constructs and their CIF files (chain A = target enzyme, chain B = binder)
CONSTRUCTS = {
    "Y87A_Y89A": {
        "erap2": "n248_trim_c5_Y87A_Y89A_erap2.cif",
        "erap1": "n248_trim_c5_Y87A_Y89A_erap1.cif",
        "irap":  "n248_trim_c5_Y87A_Y89A_irap.cif",
        "anpep": "n248_trim_c5_Y87A_Y89A_anpep.cif",
    },
    "parent_trim_c5": {
        "erap2": "n248_trim_c5_erap2.cif",
        "erap1": "n248_trim_c5_erap1.cif",
        "irap":  "n248_trim_c5_irap.cif",
        "anpep": "n248_trim_c5_anpep.cif",
    },
    "wt": {
        "erap2": "n248_wt_erap2.cif",
        "erap1": "n248_wt_erap1.cif",
        "irap":  "n248_wt_irap.cif",
        "anpep": "n248_wt_anpep.cif",
    },
}

# Known aminopeptidase inhibitor affinities for context
REFERENCE_INHIBITORS = {
    "bestatin (ERAP2)":     {"Kd_uM": 5.0,   "note": "Known aminopeptidase inhibitor, broad-spectrum"},
    "DG011A (ERAP2)":       {"Kd_uM": 0.05,  "note": "Selective ERAP2 inhibitor, published lead"},
    "leucinethiol (ERAP1)": {"Kd_uM": 0.3,   "note": "Potent ERAP1 inhibitor"},
    "ERAP2-IN-1":           {"Kd_uM": 1.2,   "note": "Literature compound"},
}

TEMP_C = 25.0  # Temperature for Kd calculation


def run_prodigy(cif_path, selection=None, temp=TEMP_C):
    """Run PRODIGY on a CIF file and return prediction dict."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", cif_path)
    model = structure[0]

    chains = [c.id for c in model.get_chains()]
    if len(chains) < 2:
        return {"error": f"Only {len(chains)} chain(s) found, need 2"}

    # PRODIGY selection format: chain pairs like ["A", "B"]
    if selection is None:
        selection = [chains[0], chains[1]]

    prod = Prodigy(model, name=os.path.basename(cif_path), selection=selection, temp=temp)
    prod.predict()
    result = prod.as_dict()

    return result


def kd_to_dg(kd_M, temp_K=298.15):
    """Convert Kd (M) to dG (kcal/mol)."""
    R = 1.987e-3  # kcal/(mol*K)
    return R * temp_K * math.log(kd_M)


def main():
    print("=" * 70)
    print("V2 PROTEIN BINDER — PRODIGY BINDING AFFINITY PREDICTION")
    print("=" * 70)

    all_results = {}

    for construct_name, targets in CONSTRUCTS.items():
        print(f"\n{'─' * 60}")
        print(f"Construct: {construct_name}")
        print(f"{'─' * 60}")
        all_results[construct_name] = {}

        for target, cif_file in targets.items():
            cif_path = os.path.join(CIF_DIR, cif_file)
            if not os.path.exists(cif_path):
                print(f"  {target:8s}: SKIP (file not found)")
                all_results[construct_name][target] = {"error": "file not found"}
                continue

            try:
                result = run_prodigy(cif_path)
                dg = result.get("ba_val", None)
                kd = result.get("kd_val", None)

                if kd is not None:
                    kd_uM = kd * 1e6
                    print(f"  {target:8s}: dG = {dg:+.2f} kcal/mol  |  Kd = {kd_uM:.3f} uM  |  Kd = {kd:.2e} M")
                else:
                    print(f"  {target:8s}: prediction failed — {result}")

                all_results[construct_name][target] = {
                    "dG_kcal_mol": round(dg, 3) if dg else None,
                    "Kd_M": kd,
                    "Kd_uM": round(kd * 1e6, 4) if kd else None,
                    "contacts": result.get("contacts", {}),
                    "nis_charged": result.get("nis_c", None),
                    "nis_polar": result.get("nis_p", None),
                    "nis_apolar": result.get("nis_a", None),
                }
            except Exception as e:
                print(f"  {target:8s}: ERROR — {e}")
                all_results[construct_name][target] = {"error": str(e)}

    # Selectivity analysis for lead construct
    print(f"\n{'=' * 70}")
    print("SELECTIVITY ANALYSIS — Y87A_Y89A")
    print(f"{'=' * 70}")

    lead = all_results.get("Y87A_Y89A", {})
    erap2_kd = lead.get("erap2", {}).get("Kd_uM")

    if erap2_kd:
        for target in ["erap1", "irap", "anpep"]:
            off_kd = lead.get(target, {}).get("Kd_uM")
            if off_kd:
                selectivity = off_kd / erap2_kd
                print(f"  ERAP2 vs {target:8s}: {selectivity:.1f}x selectivity  (ERAP2 Kd={erap2_kd:.3f} uM, {target} Kd={off_kd:.3f} uM)")
            else:
                print(f"  ERAP2 vs {target:8s}: cannot compute")

    # Compare to known inhibitors
    print(f"\n{'=' * 70}")
    print("COMPARISON TO KNOWN INHIBITORS")
    print(f"{'=' * 70}")
    print(f"  {'Compound':<30s} {'Kd (uM)':>10s}  Note")
    print(f"  {'─' * 30} {'─' * 10}  {'─' * 30}")

    if erap2_kd:
        print(f"  {'V2 Y87A_Y89A (predicted)':<30s} {erap2_kd:>10.3f}  This work (PRODIGY, Boltz-2 pose)")

    for name, data in REFERENCE_INHIBITORS.items():
        print(f"  {name:<30s} {data['Kd_uM']:>10.3f}  {data['note']}")

    # Verdict
    print(f"\n{'=' * 70}")
    print("VERDICT")
    print(f"{'=' * 70}")

    if erap2_kd:
        if erap2_kd < 0.1:
            print("  STRONG BINDER — predicted Kd in nanomolar range")
            print("  Caveat: PRODIGY on predicted structures overestimates affinity")
        elif erap2_kd < 1.0:
            print("  MODERATE-STRONG BINDER — sub-micromolar predicted Kd")
            print("  In range of known selective ERAP2 inhibitors")
        elif erap2_kd < 10.0:
            print("  MODERATE BINDER — low micromolar predicted Kd")
            print("  Comparable to bestatin class (broad-spectrum)")
        elif erap2_kd < 100.0:
            print("  WEAK BINDER — high micromolar predicted Kd")
            print("  Would likely need optimization for therapeutic use")
        else:
            print("  VERY WEAK / NON-BINDER — predicted Kd > 100 uM")

        print()
        print("  IMPORTANT CAVEATS:")
        print("  1. PRODIGY was trained on crystal structures, not Boltz-2 predictions")
        print("  2. PRODIGY predicts protein-protein Kd, not enzyme inhibition (IC50)")
        print("  3. Binding != inhibition — V2 may bind strongly but not block substrates")
        print("  4. Treat as order-of-magnitude estimate only")

    # Save results
    output = {
        "analysis": "PRODIGY binding affinity prediction for V2 binders",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "temperature_C": TEMP_C,
        "constructs": all_results,
        "reference_inhibitors": REFERENCE_INHIBITORS,
    }

    out_path = os.path.join(OUTPUT_DIR, "prodigy_affinity.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Results saved to: {out_path}")

    return output


if __name__ == "__main__":
    main()
