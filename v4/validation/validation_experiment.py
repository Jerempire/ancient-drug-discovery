"""V4 Validation Experiment: Proper multi-seed Boltz-2 with controls.

Tests whether the 4 locked synthesis candidates show real allele selectivity
by comparing them against hydrophobic controls, isostere controls, and
random peptides.

Experiment:
    4 leads + 2 hydrophobic + 2 isostere + 10 random = 18 peptides
    x 2 alleles (K392, N392) = 36 Boltz-2 runs
    x 5 diffusion samples each = 180 total samples

Usage:
    python validation_experiment.py                  # Generate YAMLs locally
    python validation_experiment.py --run-boltz2     # Generate + run on Vast.ai
    python validation_experiment.py --analyze        # Analyze after download
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import os
import random
import statistics
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
YAML_DIR = SCRIPT_DIR / "boltz2_inputs"
RESULTS_DIR = SCRIPT_DIR / "boltz2_out"

DIFFUSION_SAMPLES = 5
SEED = 42

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

# ── Experiment Definition ─────────────────────────────────────────────────────

LEADS = [
    {"name": "lead_pep_glu_long_01",    "sequence": "EALVAAGLAGLA", "p1": "E", "arm": "K392", "group": "lead", "source": "Hand-designed"},
    {"name": "lead_hybrid_E_VKLLLL",    "sequence": "EKLLLLSIGK",   "p1": "E", "arm": "K392", "group": "lead", "source": "PepMLM hybrid"},
    {"name": "lead_pep_ala_01",         "sequence": "AALVAAGLA",     "p1": "A", "arm": "N392", "group": "lead", "source": "Hand-designed"},
    {"name": "lead_pep_leu_01",         "sequence": "LALVAAGLA",     "p1": "L", "arm": "N392", "group": "lead", "source": "Hand-designed"},
]

# Hydrophobic controls: expect N392 preference (matches published 165x hydrophobic bias)
HYDROPHOBIC_CONTROLS = [
    {"name": "ctrl_hydro_leu_long",  "sequence": "LALVAAGLAGLA", "p1": "L", "arm": "N392-expected", "group": "hydrophobic", "source": "L at P1, same scaffold as lead #1"},
    {"name": "ctrl_hydro_val",       "sequence": "VALVAAGLA",    "p1": "V", "arm": "N392-expected", "group": "hydrophobic", "source": "V at P1, same scaffold as leads #3-4"},
]

# Isostere controls: E->Q removes charge but keeps size/polarity
# Tests: is the salt bridge (charge) driving K392 selectivity?
ISOSTERE_CONTROLS = [
    {"name": "ctrl_iso_gln_long",    "sequence": "QALVAAGLAGLA", "p1": "Q", "arm": "neutral-expected", "group": "isostere", "source": "E->Q on lead #1 (charge removal)"},
    {"name": "ctrl_iso_gln_hybrid",  "sequence": "QKLLLLSIGK",   "p1": "Q", "arm": "neutral-expected", "group": "isostere", "source": "E->Q on lead #2 (charge removal)"},
]

# Random peptides: background distribution (fixed seed for reproducibility)
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
RANDOM_LENGTHS = [9, 9, 9, 10, 10, 10, 11, 11, 11, 10]  # 10 peptides, mixed lengths


def generate_random_peptides(n=10, seed=12345):
    """Generate n random peptides with reproducible seed."""
    rng = random.Random(seed)
    peptides = []
    for i in range(n):
        length = RANDOM_LENGTHS[i]
        seq = "".join(rng.choice(STANDARD_AA) for _ in range(length))
        peptides.append({
            "name": f"ctrl_random_{i:02d}",
            "sequence": seq,
            "p1": seq[0],
            "arm": "unknown",
            "group": "random",
            "source": f"Random {length}-mer (seed={seed})",
        })
    return peptides


def get_all_peptides():
    """Return the full experiment peptide list."""
    randoms = generate_random_peptides()
    all_peps = LEADS + HYDROPHOBIC_CONTROLS + ISOSTERE_CONTROLS + randoms
    return all_peps


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


def print_experiment_table(peptides):
    """Print the full experiment table."""
    print(f"{'#':>3s}  {'Name':>28s}  {'Sequence':>14s}  {'Len':>3s}  {'P1':>2s}  {'Group':>12s}  {'Expected':>16s}")
    print("-" * 95)
    for i, pep in enumerate(peptides, 1):
        print(f"{i:>3d}  {pep['name']:>28s}  {pep['sequence']:>14s}  {len(pep['sequence']):>3d}  {pep['p1']:>2s}  {pep['group']:>12s}  {pep['arm']:>16s}")


# ── Analysis ──────────────────────────────────────────────────────────────────

def collect_scores(results_dir):
    """Collect all ipTM/pTM scores from Boltz-2 output."""
    results_dir = Path(results_dir)
    scores = {}  # key: "peptide_name_VARIANT" -> list of ipTM

    for json_file in sorted(results_dir.rglob("confidence_*.json")):
        try:
            data = json.loads(json_file.read_text())
            iptm = data.get("iptm", data.get("ipTM"))
            ptm = data.get("ptm", data.get("pTM"))
            if iptm is None:
                continue

            # Walk up to find the run directory name
            # Structure: boltz2_out/<peptide_VARIANT>/boltz_results_<...>/predictions/.../confidence_*.json
            # OR: boltz2_out/<peptide_VARIANT>/predictions/.../confidence_*.json
            run_name = None
            for parent in json_file.parents:
                if parent.parent == results_dir or parent == results_dir:
                    run_name = parent.name
                    break

            if run_name is None:
                continue

            if run_name not in scores:
                scores[run_name] = []
            scores[run_name].append({"iptm": float(iptm), "ptm": float(ptm)})
        except (json.JSONDecodeError, ValueError):
            continue

    return scores


def analyze_results(results_dir):
    """Full analysis: per-peptide stats, lead vs random comparison, decision."""
    results_dir = Path(results_dir)
    scores = collect_scores(results_dir)

    if not scores:
        print("No confidence scores found. Check results directory.")
        print(f"Searched: {results_dir}")
        return

    # Build paired comparisons
    peptides = get_all_peptides()
    pep_lookup = {p["name"]: p for p in peptides}

    pairs = {}  # peptide_name -> {K392: [iptm...], N392: [iptm...], meta: {...}}
    for run_name, sample_scores in scores.items():
        # Parse: lead_pep_glu_long_01_K392
        if run_name.endswith("_K392"):
            pep_name = run_name[:-5]
            variant = "K392"
        elif run_name.endswith("_N392"):
            pep_name = run_name[:-5]
            variant = "N392"
        else:
            continue

        if pep_name not in pairs:
            meta = pep_lookup.get(pep_name, {"group": "unknown", "arm": "unknown", "p1": "?"})
            pairs[pep_name] = {"K392": [], "N392": [], "meta": meta}

        pairs[pep_name][variant] = [s["iptm"] for s in sample_scores]

    # ── Per-peptide statistics ────────────────────────────────────────────────
    print("=" * 110)
    print("VALIDATION EXPERIMENT RESULTS")
    print("=" * 110)
    print()

    rows = []
    for pep_name in sorted(pairs, key=lambda k: (pairs[k]["meta"].get("group", "z"), k)):
        p = pairs[pep_name]
        k_scores = p["K392"]
        n_scores = p["N392"]
        meta = p["meta"]

        if not k_scores or not n_scores:
            continue

        k_mean = statistics.mean(k_scores)
        n_mean = statistics.mean(n_scores)
        delta = k_mean - n_mean

        # Success ratio: what fraction of seeds favor K392?
        # Compare seed-by-seed if same count, otherwise just use means
        if len(k_scores) == len(n_scores):
            k_wins = sum(1 for k, n in zip(k_scores, n_scores) if k > n)
            success_ratio = k_wins / len(k_scores)
        else:
            success_ratio = float(delta > 0)

        k_std = statistics.stdev(k_scores) if len(k_scores) > 1 else 0
        n_std = statistics.stdev(n_scores) if len(n_scores) > 1 else 0

        row = {
            "name": pep_name,
            "group": meta.get("group", "unknown"),
            "p1": meta.get("p1", "?"),
            "expected_arm": meta.get("arm", "?"),
            "sequence": meta.get("sequence", "?"),
            "k392_mean": k_mean,
            "k392_std": k_std,
            "k392_n": len(k_scores),
            "n392_mean": n_mean,
            "n392_std": n_std,
            "n392_n": len(n_scores),
            "delta": delta,
            "success_ratio": success_ratio,
        }
        rows.append(row)

    # Print by group
    for group_name in ["lead", "hydrophobic", "isostere", "random"]:
        group_rows = [r for r in rows if r["group"] == group_name]
        if not group_rows:
            continue

        label = {
            "lead": "LEADS (synthesis candidates)",
            "hydrophobic": "HYDROPHOBIC CONTROLS (expect N392)",
            "isostere": "ISOSTERE CONTROLS (E->Q, tests charge)",
            "random": "RANDOM PEPTIDES (background)",
        }[group_name]

        print(f"### {label}")
        print(f"{'Name':>28s}  {'P1':>2s}  {'K392 mean':>10s}  {'N392 mean':>10s}  {'Delta':>8s}  {'K wins':>7s}  {'Verdict':>12s}")
        print("-" * 90)

        for r in group_rows:
            if r["delta"] > 0.03:
                verdict = "K392-SEL"
            elif r["delta"] < -0.03:
                verdict = "N392-SEL"
            else:
                verdict = "NEUTRAL"

            sr_str = f"{r['success_ratio']:.0%}"
            print(f"{r['name']:>28s}  {r['p1']:>2s}  "
                  f"{r['k392_mean']:>7.3f}+/-{r['k392_std']:.3f}"
                  f"  {r['n392_mean']:>7.3f}+/-{r['n392_std']:.3f}"
                  f"  {r['delta']:>+8.3f}  {sr_str:>7s}  {verdict:>12s}")
        print()

    # ── Lead vs Random comparison ─────────────────────────────────────────────
    random_deltas = [r["delta"] for r in rows if r["group"] == "random"]
    lead_rows = [r for r in rows if r["group"] == "lead"]

    if random_deltas and lead_rows:
        rand_mean = statistics.mean(random_deltas)
        rand_std = statistics.stdev(random_deltas) if len(random_deltas) > 1 else 0.01

        print("### LEAD vs RANDOM COMPARISON")
        print(f"Random background: mean delta = {rand_mean:+.4f}, std = {rand_std:.4f} (n={len(random_deltas)})")
        print()
        print(f"{'Lead':>28s}  {'Delta':>8s}  {'Z-score':>8s}  {'Outlier?':>10s}")
        print("-" * 62)

        for r in lead_rows:
            z = (r["delta"] - rand_mean) / rand_std if rand_std > 0 else 0
            outlier = "YES" if abs(z) > 2.0 else ("MAYBE" if abs(z) > 1.5 else "NO")
            print(f"{r['name']:>28s}  {r['delta']:>+8.3f}  {z:>+8.2f}  {outlier:>10s}")
        print()

    # ── Isostere check ────────────────────────────────────────────────────────
    iso_rows = [r for r in rows if r["group"] == "isostere"]
    if iso_rows and lead_rows:
        print("### ISOSTERE CHECK (does charge drive selectivity?)")
        print("If E->Q kills K392 preference, charge (salt bridge) is the mechanism.")
        print()
        # Compare lead E-peptides vs their Q-isosteres
        for iso in iso_rows:
            # Find matching lead
            if "long" in iso["name"]:
                matching_lead = next((l for l in lead_rows if "glu_long" in l["name"]), None)
            elif "hybrid" in iso["name"]:
                matching_lead = next((l for l in lead_rows if "hybrid" in l["name"]), None)
            else:
                matching_lead = None

            if matching_lead:
                shift = iso["delta"] - matching_lead["delta"]
                print(f"  {matching_lead['name']} (E): delta = {matching_lead['delta']:+.3f}")
                print(f"  {iso['name']} (Q): delta = {iso['delta']:+.3f}")
                print(f"  Shift from E->Q: {shift:+.3f}")
                if shift < -0.03:
                    print(f"  -> Charge IS driving K392 selectivity (delta dropped)")
                elif abs(shift) < 0.03:
                    print(f"  -> Charge NOT the main driver (delta unchanged)")
                else:
                    print(f"  -> Unexpected: Q is MORE K392-selective than E")
                print()

    # ── Decision ──────────────────────────────────────────────────────────────
    print("=" * 80)
    print("DECISION")
    print("=" * 80)

    # Classify outcome
    strong_leads = [r for r in lead_rows if abs(r["delta"]) > 0.05 and r["success_ratio"] >= 0.6]
    weak_leads = [r for r in lead_rows if 0.03 < abs(r["delta"]) <= 0.05]

    if strong_leads:
        # Check if they beat randoms
        if random_deltas:
            rand_std_safe = rand_std if rand_std > 0 else 0.01
            outlier_leads = [r for r in strong_leads
                             if abs((r["delta"] - rand_mean) / rand_std_safe) > 2.0]
            if outlier_leads:
                print("CASE 1: STRONG SIGNAL")
                print(f"  {len(outlier_leads)} lead(s) are consistent, beat randoms:")
                for r in outlier_leads:
                    print(f"    {r['name']}: delta={r['delta']:+.3f}, K-wins={r['success_ratio']:.0%}")
                print("  -> Proceed to MD on top 1-2 candidates")
            else:
                print("CASE 2: WEAK BUT CONSISTENT")
                print("  Leads are consistent but don't clearly separate from randoms.")
                print("  -> Test top 1-2 in MD before deciding")
        else:
            print("CASE 1: STRONG SIGNAL (no random baseline to compare)")
            for r in strong_leads:
                print(f"    {r['name']}: delta={r['delta']:+.3f}, K-wins={r['success_ratio']:.0%}")
    elif weak_leads:
        print("CASE 2: WEAK BUT CONSISTENT")
        print(f"  {len(weak_leads)} lead(s) show small but stable signal:")
        for r in weak_leads:
            print(f"    {r['name']}: delta={r['delta']:+.3f}, K-wins={r['success_ratio']:.0%}")
        print("  -> Test 1-2 peptides in MD")
    else:
        print("CASE 3: NO SIGNAL")
        print("  No leads show consistent allele selectivity.")
        print("  -> Stop V4 or redesign")
    print()

    # ── Save results ──────────────────────────────────────────────────────────
    output = {
        "experiment": {
            "diffusion_samples": DIFFUSION_SAMPLES,
            "seed": SEED,
            "n_leads": len(lead_rows),
            "n_hydrophobic": len([r for r in rows if r["group"] == "hydrophobic"]),
            "n_isostere": len(iso_rows),
            "n_random": len(random_deltas),
        },
        "results": rows,
        "random_baseline": {
            "mean_delta": rand_mean if random_deltas else None,
            "std_delta": rand_std if random_deltas else None,
            "n": len(random_deltas),
        },
    }
    output_path = results_dir / "validation_analysis.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"Saved analysis to {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(description="V4 Validation Experiment")
    parser.add_argument("--analyze", action="store_true", help="Analyze existing results")
    parser.add_argument("--run-boltz2", action="store_true", help="Generate YAMLs and run Boltz-2")
    parser.add_argument("--results-dir", type=str, default=str(RESULTS_DIR),
                        help="Results directory for analysis")
    args = parser.parse_args()

    if args.analyze:
        analyze_results(Path(args.results_dir))
        return

    # Generate experiment
    peptides = get_all_peptides()

    print("=" * 95)
    print("V4 VALIDATION EXPERIMENT")
    print(f"  {len(peptides)} peptides x 2 alleles = {len(peptides) * 2} runs x {DIFFUSION_SAMPLES} samples")
    print("=" * 95)
    print()

    print_experiment_table(peptides)
    print()

    # Write YAMLs
    n = write_boltz2_yamls(peptides, YAML_DIR)
    print(f"Wrote {n} Boltz-2 input YAMLs to {YAML_DIR}/")

    # Save experiment definition
    defn_path = SCRIPT_DIR / "experiment_definition.json"
    with open(defn_path, "w") as f:
        json.dump(peptides, f, indent=2)
    print(f"Saved experiment definition to {defn_path}")

    if args.run_boltz2:
        print(f"\nRunning Boltz-2 ({DIFFUSION_SAMPLES} diffusion samples, seed={SEED})...")
        for yaml_file in sorted(YAML_DIR.glob("*.yaml")):
            name = yaml_file.stem
            out = RESULTS_DIR / name
            if out.exists() and list(out.rglob("confidence_*.json")):
                print(f"  {name} — already done, skipping")
                continue
            print(f"  {name} ...")
            os.system(
                f'boltz predict "{yaml_file}" '
                f'--diffusion_samples {DIFFUSION_SAMPLES} '
                f'--seed {SEED} '
                f'--out_dir "{out}" '
                f'--no_kernels 2>&1 | tail -2'
            )


if __name__ == "__main__":
    main()
