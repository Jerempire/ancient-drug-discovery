"""
Multi-seed confirmation for top V4 peptide candidates.

Tests whether K392/N392 selectivity survives across 5 random seeds.
This is the noise-vs-signal gate before committing to synthesis.

4 peptides x 2 variants x 5 seeds = 40 Boltz-2 predictions.
Expected runtime: ~30 min on RTX 4090.

Usage (on Vast.ai):
    python3 /workspace/multiseed_confirm.py
"""
import subprocess
import json
import os
import glob
import sys
from collections import defaultdict
from statistics import mean, stdev

# --- CONFIG ---
SEEDS = [42, 123, 456, 789, 1337]
DIFFUSION_SAMPLES = 3
OUT_BASE = "/workspace/results/multiseed"

# ERAP2 channel sequences (residues 350-450)
ERAP2_FULL = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
K392_CHANNEL = ERAP2_FULL[349:450]
N392_CHANNEL = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:450]

# 4 candidates to confirm: 2 K392 leaders + 2 N392 leaders
CANDIDATES = {
    "pep_glu_long_01": {"seq": "EALVAAGLAGLA", "arm": "K392", "source": "V4 original"},
    "hybrid_E_VKLLLL": {"seq": "EKLLLLSIGK", "arm": "K392", "source": "PepMLM hybrid"},
    "pep_leu_01":      {"seq": "LALVAAGLA",  "arm": "N392", "source": "V4 original"},
    "pep_ala_01":      {"seq": "AALVAAGLA",  "arm": "N392", "source": "V4 original"},
}


def write_yaml(name, variant, seed, pep_seq, channel_seq):
    """Write a Boltz-2 input YAML."""
    yaml_dir = os.path.join(OUT_BASE, "yamls")
    os.makedirs(yaml_dir, exist_ok=True)
    fname = f"{name}_{variant}_seed{seed}.yaml"
    path = os.path.join(yaml_dir, fname)
    with open(path, "w") as f:
        f.write(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {channel_seq}
      msa: empty
  - protein:
      id: B
      sequence: {pep_seq}
      msa: empty
""")
    return path


def run_boltz(yaml_path, out_dir, seed):
    """Run Boltz-2 prediction."""
    cmd = [
        "boltz", "predict", yaml_path,
        "--out_dir", out_dir,
        "--diffusion_samples", str(DIFFUSION_SAMPLES),
        "--seed", str(seed),
        "--no_kernels",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-200:]}")
    return result.returncode == 0


def extract_iptm(out_dir):
    """Extract average ipTM from confidence JSONs."""
    conf_files = glob.glob(os.path.join(out_dir, "**/confidence_*.json"), recursive=True)
    iptms = []
    for cf in conf_files:
        with open(cf) as f:
            c = json.load(f)
        iptm = c.get("iptm", c.get("i_ptm", 0))
        if iptm:
            iptms.append(iptm)
    return mean(iptms) if iptms else 0


def main():
    os.makedirs(OUT_BASE, exist_ok=True)
    total = len(CANDIDATES) * 2 * len(SEEDS)
    done = 0

    print(f"Multi-seed confirmation: {len(CANDIDATES)} peptides x 2 variants x {len(SEEDS)} seeds = {total} predictions")
    print(f"Seeds: {SEEDS}")
    print()

    # Run all predictions
    all_scores = {}  # {name: {variant: [iptm_per_seed]}}

    for name, info in CANDIDATES.items():
        all_scores[name] = {"K392": [], "N392": [], "info": info}

        for variant, channel_seq in [("K392", K392_CHANNEL), ("N392", N392_CHANNEL)]:
            for seed in SEEDS:
                done += 1
                tag = f"{name}_{variant}_seed{seed}"
                out_dir = os.path.join(OUT_BASE, tag)
                print(f"[{done}/{total}] {tag}...", end=" ", flush=True)

                yaml_path = write_yaml(name, variant, seed, info["seq"], channel_seq)
                ok = run_boltz(yaml_path, out_dir, seed)

                if ok:
                    iptm = extract_iptm(out_dir)
                    all_scores[name][variant].append(iptm)
                    print(f"ipTM={iptm:.3f}")
                else:
                    all_scores[name][variant].append(0)
                    print("FAILED")

    # --- Analysis ---
    SEP = "=" * 80
    DASH = "-" * 80

    print("\n" + SEP)
    print("MULTI-SEED CONFIRMATION RESULTS")
    print(SEP)
    print("%-22s %6s  %20s  %20s  %8s  %8s" % (
        "Peptide", "Arm", "K392 (mean +/- sd)", "N392 (mean +/- sd)", "Delta", "Stable?"))
    print(DASH)

    results_summary = []
    for name, data in all_scores.items():
        info = data["info"]
        k_scores = data["K392"]
        n_scores = data["N392"]

        k_mean = mean(k_scores) if k_scores else 0
        n_mean = mean(n_scores) if n_scores else 0
        k_sd = stdev(k_scores) if len(k_scores) > 1 else 0
        n_sd = stdev(n_scores) if len(n_scores) > 1 else 0
        delta = k_mean - n_mean

        # Signal is stable if delta direction is consistent and magnitude > 2x combined SD
        combined_sd = (k_sd**2 + n_sd**2)**0.5
        stable = abs(delta) > 2 * combined_sd and abs(delta) > 0.02

        direction = "K392" if delta > 0.02 else "N392" if delta < -0.02 else "neutral"

        print("%-22s %6s  %8.3f +/- %.3f    %8.3f +/- %.3f    %+8.3f  %8s" % (
            name, info["arm"], k_mean, k_sd, n_mean, n_sd, delta,
            "YES" if stable else "NO"))

        results_summary.append({
            "name": name,
            "sequence": info["seq"],
            "arm": info["arm"],
            "source": info["source"],
            "K392_scores": k_scores,
            "N392_scores": n_scores,
            "K392_mean": k_mean,
            "K392_sd": k_sd,
            "N392_mean": n_mean,
            "N392_sd": n_sd,
            "delta": delta,
            "combined_sd": combined_sd,
            "direction": direction,
            "stable": stable,
        })

    # Verdict
    print("\n" + SEP)
    print("VERDICT")
    print(SEP)
    stable_count = sum(1 for r in results_summary if r["stable"])
    print(f"{stable_count}/{len(results_summary)} candidates have stable selectivity across {len(SEEDS)} seeds")
    for r in results_summary:
        status = "CONFIRMED" if r["stable"] else "NOISY — selectivity may be artifact"
        print(f"  {r['name']}: {r['direction']} delta={r['delta']:+.3f} sd={r['combined_sd']:.3f} -> {status}")

    # Save
    out_path = os.path.join(OUT_BASE, "multiseed_results.json")
    with open(out_path, "w") as f:
        json.dump(results_summary, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
