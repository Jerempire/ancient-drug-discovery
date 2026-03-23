"""Structured SAR grid for ERAP2 K392/N392 variant selectivity.

Systematically varies scaffold, P1 residue, and length to answer:
1. Where is the length switch (N392→K392 preference)?
2. Does P1 chemistry or scaffold drive selectivity?
3. Is the N392-selective arm real and general?

Grid: 2 scaffolds x 5 P1 residues x 4 lengths = 40 peptides
      x 2 variants (K392, N392) = 80 Boltz-2 runs
      x 3 diffusion samples each

Usage:
    python sar_grid.py                    # Generate YAMLs locally
    python sar_grid.py --run-boltz2       # Generate + run Boltz-2 (on Vast.ai)
    python sar_grid.py --analyze          # Analyze results after download
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import os
from pathlib import Path
from itertools import product

SCRIPT_DIR = Path(__file__).resolve().parent
SAR_DIR = SCRIPT_DIR / "sar_results"

# ── ERAP2 sequences ──────────────────────────────────────────────────────────
ERAP2_K392 = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPS"
    "VVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYP"
    "AHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLF"
    "KANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASP"
    "DKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLW"
    "VTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPI"
    "SKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTS"
    "GGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQER"
    "YLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYF"
    "KPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAG"
    "WNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRE"
    "NWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLR"
    "TWLMVNT"
)
ERAP2_N392 = ERAP2_K392[:391] + "N" + ERAP2_K392[392:]

# ── SAR Grid Definition ──────────────────────────────────────────────────────
# Scaffolds: the body after P1 (from PepMLM's top hits)
SCAFFOLDS = {
    "VKLLLL": "KLLLLSIGK",    # PepMLM top binder (original: VKLLLLSIGK, 10-mer)
    "ATSKSK": "TSKSKSSRSSTVKK",  # PepMLM best selectivity (original: ATSKSKSSRSSTVKK, 15-mer)
}

P1_RESIDUES = ["E", "D", "A", "L", "V"]

LENGTHS = [8, 9, 10, 11]


def build_peptide(scaffold_body, p1, target_length):
    """Build a peptide: P1 + scaffold_body, truncated or padded to target_length."""
    full = p1 + scaffold_body
    if len(full) >= target_length:
        return full[:target_length]
    else:
        # Extend with the scaffold pattern (repeat last few residues)
        while len(full) < target_length:
            full += scaffold_body[len(full) % len(scaffold_body)]
        return full[:target_length]


def generate_grid():
    """Generate all SAR grid peptides."""
    peptides = []
    for scaffold_name, scaffold_body in SCAFFOLDS.items():
        for p1 in P1_RESIDUES:
            for length in LENGTHS:
                seq = build_peptide(scaffold_body, p1, length)
                name = f"sar_{scaffold_name}_{p1}_{length}aa"
                peptides.append({
                    "name": name,
                    "sequence": seq,
                    "scaffold": scaffold_name,
                    "p1": p1,
                    "length": length,
                })
    return peptides


def write_boltz2_yamls(peptides, output_dir):
    """Write Boltz-2 input YAMLs for all peptides x both variants."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for pep in peptides:
        for variant_name, erap_seq in [("K392", ERAP2_K392), ("N392", ERAP2_N392)]:
            yaml_name = f"{pep['name']}_{variant_name}.yaml"
            yaml_content = (
                f"sequences:\n"
                f"  - protein:\n"
                f"      id: A\n"
                f"      sequence: {erap_seq}\n"
                f"      msa: empty\n"
                f"  - protein:\n"
                f"      id: B\n"
                f"      sequence: {pep['sequence']}\n"
                f"      msa: empty\n"
            )
            (output_dir / yaml_name).write_text(yaml_content)
            count += 1

    return count


def analyze_results(results_dir):
    """Analyze Boltz-2 SAR results and build the selectivity table."""
    results_dir = Path(results_dir)

    # Collect all confidence scores
    scores = {}
    for json_file in sorted(results_dir.glob("**/*confidence*.json")):
        try:
            data = json.loads(json_file.read_text())
            iptm = data.get("iptm", data.get("ipTM"))
            if iptm is None:
                continue

            # Parse name from path: sar_SCAFFOLD_P1_LENaa_VARIANT/...
            parts = json_file.parts
            for part in parts:
                if part.startswith("sar_"):
                    name = part
                    break
            else:
                continue

            if name not in scores:
                scores[name] = []
            scores[name].append(float(iptm))
        except (json.JSONDecodeError, ValueError):
            continue

    if not scores:
        print("No confidence scores found. Check results directory structure.")
        return

    # Average across diffusion samples
    avg_scores = {name: sum(vals) / len(vals) for name, vals in scores.items()}

    # Parse into structured data
    rows = []
    seen_peptides = set()
    for full_name, avg_iptm in avg_scores.items():
        # Parse: sar_SCAFFOLD_P1_LENaa_VARIANT
        # Example: sar_VKLLLL_E_10aa_K392
        parts = full_name.rsplit("_", 1)
        if len(parts) != 2:
            continue
        pep_name, variant = parts
        if pep_name not in seen_peptides:
            seen_peptides.add(pep_name)

        # Parse scaffold, P1, length from pep_name
        # sar_VKLLLL_E_10aa
        pep_parts = pep_name.split("_")
        if len(pep_parts) >= 4:
            scaffold = pep_parts[1]
            p1 = pep_parts[2]
            length = int(pep_parts[3].replace("aa", ""))
        else:
            continue

        rows.append({
            "name": pep_name,
            "scaffold": scaffold,
            "p1": p1,
            "length": length,
            "variant": variant,
            "iptm": avg_iptm,
            "n_samples": len(scores[full_name]),
        })

    # Build paired K392/N392 comparisons
    pairs = {}
    for row in rows:
        key = row["name"]
        if key not in pairs:
            pairs[key] = {"scaffold": row["scaffold"], "p1": row["p1"], "length": row["length"]}
        pairs[key][row["variant"]] = row["iptm"]
        pairs[key][f"{row['variant']}_n"] = row["n_samples"]

    # Calculate deltas and print SAR table
    print("=" * 90)
    print("SAR GRID RESULTS: ERAP2 K392/N392 Variant Selectivity")
    print("=" * 90)
    print()

    # Group by scaffold
    for scaffold_name in SCAFFOLDS:
        print(f"### Scaffold: {scaffold_name}")
        print(f"{'Name':>30s}  {'P1':>3s}  {'Len':>4s}  {'K392':>8s}  {'N392':>8s}  {'Delta':>8s}  {'Verdict':>10s}")
        print("-" * 80)

        scaffold_pairs = {k: v for k, v in pairs.items() if v["scaffold"] == scaffold_name}

        for key in sorted(scaffold_pairs, key=lambda k: (pairs[k]["p1"], pairs[k]["length"])):
            p = pairs[key]
            k392 = p.get("K392", 0)
            n392 = p.get("N392", 0)
            delta = k392 - n392

            if delta > 0.03:
                verdict = "K392-SEL"
            elif delta < -0.03:
                verdict = "N392-SEL"
            else:
                verdict = "NEUTRAL"

            print(f"{key:>30s}  {p['p1']:>3s}  {p['length']:>4d}  {k392:>8.3f}  {n392:>8.3f}  {delta:>+8.3f}  {verdict:>10s}")

        print()

    # Summary: P1 effect (averaged across scaffolds and lengths)
    print("### P1 Effect (mean delta across all scaffolds and lengths)")
    p1_deltas = {}
    for p in pairs.values():
        k392 = p.get("K392", 0)
        n392 = p.get("N392", 0)
        delta = k392 - n392
        p1 = p["p1"]
        if p1 not in p1_deltas:
            p1_deltas[p1] = []
        p1_deltas[p1].append(delta)

    print(f"{'P1':>5s}  {'Mean Delta':>12s}  {'N':>4s}  {'Direction':>12s}")
    print("-" * 40)
    for p1 in P1_RESIDUES:
        if p1 in p1_deltas:
            mean_d = sum(p1_deltas[p1]) / len(p1_deltas[p1])
            direction = "K392" if mean_d > 0.01 else ("N392" if mean_d < -0.01 else "NEUTRAL")
            print(f"{p1:>5s}  {mean_d:>+12.4f}  {len(p1_deltas[p1]):>4d}  {direction:>12s}")
    print()

    # Summary: Length effect (averaged across scaffolds and P1)
    print("### Length Effect (mean delta across all scaffolds and P1)")
    len_deltas = {}
    for p in pairs.values():
        k392 = p.get("K392", 0)
        n392 = p.get("N392", 0)
        delta = k392 - n392
        length = p["length"]
        if length not in len_deltas:
            len_deltas[length] = []
        len_deltas[length].append(delta)

    print(f"{'Length':>8s}  {'Mean Delta':>12s}  {'N':>4s}  {'Direction':>12s}")
    print("-" * 40)
    for length in sorted(len_deltas):
        mean_d = sum(len_deltas[length]) / len(len_deltas[length])
        direction = "K392" if mean_d > 0.01 else ("N392" if mean_d < -0.01 else "NEUTRAL")
        print(f"{length:>8d}  {mean_d:>+12.4f}  {len(len_deltas[length]):>4d}  {direction:>12s}")
    print()

    # Save results
    output = {
        "pairs": pairs,
        "p1_effect": {p1: {"mean_delta": sum(d)/len(d), "n": len(d)} for p1, d in p1_deltas.items()},
        "length_effect": {str(l): {"mean_delta": sum(d)/len(d), "n": len(d)} for l, d in len_deltas.items()},
    }
    output_path = results_dir / "sar_analysis.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved analysis to {output_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="SAR grid for ERAP2 variant selectivity")
    parser.add_argument("--analyze", action="store_true", help="Analyze existing results")
    parser.add_argument("--run-boltz2", action="store_true", help="Generate YAMLs and run Boltz-2")
    parser.add_argument("--results-dir", type=str, default=str(SAR_DIR / "boltz2_out"),
                        help="Results directory for analysis")
    args = parser.parse_args()

    if args.analyze:
        analyze_results(args.results_dir)
        return

    # Generate grid
    peptides = generate_grid()
    print(f"SAR Grid: {len(peptides)} peptides")
    print(f"  Scaffolds: {list(SCAFFOLDS.keys())}")
    print(f"  P1 residues: {P1_RESIDUES}")
    print(f"  Lengths: {LENGTHS}")
    print()

    # Print grid
    print(f"{'Name':>30s}  {'Sequence':>20s}  {'Scaffold':>8s}  {'P1':>3s}  {'Len':>4s}")
    print("-" * 75)
    for pep in peptides:
        print(f"{pep['name']:>30s}  {pep['sequence']:>20s}  {pep['scaffold']:>8s}  {pep['p1']:>3s}  {pep['length']:>4d}")
    print()

    # Write YAMLs
    yaml_dir = SAR_DIR / "boltz2_inputs"
    n = write_boltz2_yamls(peptides, yaml_dir)
    print(f"Wrote {n} Boltz-2 input YAMLs to {yaml_dir}/")

    # Save grid definition
    grid_path = SAR_DIR / "sar_grid_definition.json"
    SAR_DIR.mkdir(parents=True, exist_ok=True)
    with open(grid_path, "w") as f:
        json.dump(peptides, f, indent=2)
    print(f"Saved grid definition to {grid_path}")

    if args.run_boltz2:
        print("\nRunning Boltz-2...")
        os.system(f"boltz predict {yaml_dir}/ --diffusion_samples 3 --seed 42 "
                  f"--out_dir {SAR_DIR}/boltz2_out/ --no_kernels 2>&1 | tail -20")
        print("\nBoltz-2 complete. Run with --analyze to see results.")


if __name__ == "__main__":
    main()
