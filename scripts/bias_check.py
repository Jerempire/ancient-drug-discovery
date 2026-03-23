"""
Boltz-2 K392 bias check + VKLLLL P1 scan.

Part 1: Do random/control peptides also score higher on K392?
Part 2: Does P1 chemistry drive selectivity on the VKLLLL scaffold?

26 predictions total. ~15 min on RTX 4090.
"""
import subprocess
import json
import os
import glob
import random
from statistics import mean
from collections import defaultdict

OUT_BASE = "/workspace/results/bias_check"
SEED = 42
DIFFUSION_SAMPLES = 3

ERAP2_FULL = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
K392_CHANNEL = ERAP2_FULL[349:450]
N392_CHANNEL = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:450]

# --- PART 1: BIAS CHECK CONTROLS ---
random.seed(42)
AA = "ACDEFGHIKLMNPQRSTVWY"

def scramble(seq):
    """Scramble a sequence (same AA composition, random order)."""
    chars = list(seq)
    random.shuffle(chars)
    return "".join(chars)

def random_peptide(length):
    """Generate random peptide."""
    return "".join(random.choice(AA) for _ in range(length))

BIAS_CONTROLS = {
    # Random peptides (matched lengths to our leads)
    "random_10mer": random_peptide(10),
    "random_11mer": random_peptide(11),
    "random_9mer": random_peptide(9),
    # Scrambled versions of our leads
    "scramble_glu_long": scramble("EALVAAGLAGLA"),   # scramble of pep_glu_long_01
    "scramble_VKLLLL": scramble("EKLLLLSIGK"),       # scramble of hybrid_E_VKLLLL
    "scramble_leu_9mer": scramble("LALVAAGLA"),       # scramble of pep_leu_01
    # Poly-Ala controls (zero chemistry)
    "polyala_10mer": "AAAAAAAAAA",
    "polyala_11mer": "AAAAAAAAAAA",
}

# --- PART 2: VKLLLL P1 SCAN ---
VKLLLL_BACKBONE = "KLLLLSIGK"  # 9 residues after P1
P1_SCAN = {
    "vkllll_P1_E": "E" + VKLLLL_BACKBONE,  # Glu (current lead)
    "vkllll_P1_D": "D" + VKLLLL_BACKBONE,  # Asp
    "vkllll_P1_A": "A" + VKLLLL_BACKBONE,  # Ala (neutral)
    "vkllll_P1_L": "L" + VKLLLL_BACKBONE,  # Leu (hydrophobic)
    "vkllll_P1_V": "V" + VKLLLL_BACKBONE,  # Val (hydrophobic)
}

ALL_PEPTIDES = {}
ALL_PEPTIDES.update(BIAS_CONTROLS)
ALL_PEPTIDES.update(P1_SCAN)


def write_yaml(name, variant, pep_seq, channel_seq):
    yaml_dir = os.path.join(OUT_BASE, "yamls")
    os.makedirs(yaml_dir, exist_ok=True)
    path = os.path.join(yaml_dir, f"{name}_{variant}.yaml")
    with open(path, "w") as f:
        f.write(f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {channel_seq}\n      msa: empty\n  - protein:\n      id: B\n      sequence: {pep_seq}\n      msa: empty\n")
    return path


def run_boltz(yaml_path, out_dir):
    cmd = ["boltz", "predict", yaml_path, "--out_dir", out_dir,
           "--diffusion_samples", str(DIFFUSION_SAMPLES), "--seed", str(SEED), "--no_kernels"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def extract_iptm(out_dir):
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
    total = len(ALL_PEPTIDES) * 2
    done = 0

    print("BOLTZ-2 K392 BIAS CHECK + VKLLLL P1 SCAN")
    print("=" * 60)

    # Print what we're testing
    print("\nPart 1 — Bias controls:")
    for name, seq in BIAS_CONTROLS.items():
        print(f"  {name}: {seq} ({len(seq)} aa)")
    print("\nPart 2 — VKLLLL P1 scan:")
    for name, seq in P1_SCAN.items():
        print(f"  {name}: {seq} ({len(seq)} aa)")
    print(f"\nTotal: {total} predictions\n")

    scores = {}
    for name, seq in ALL_PEPTIDES.items():
        scores[name] = {}
        for variant, channel in [("K392", K392_CHANNEL), ("N392", N392_CHANNEL)]:
            done += 1
            tag = f"{name}_{variant}"
            out_dir = os.path.join(OUT_BASE, tag)
            print(f"[{done}/{total}] {tag}...", end=" ", flush=True)
            yaml_path = write_yaml(name, variant, seq, channel)
            ok = run_boltz(yaml_path, out_dir)
            if ok:
                iptm = extract_iptm(out_dir)
                scores[name][variant] = iptm
                print(f"ipTM={iptm:.3f}")
            else:
                scores[name][variant] = 0
                print("FAILED")

    # --- ANALYSIS ---
    SEP = "=" * 70
    DASH = "-" * 70

    print("\n" + SEP)
    print("PART 1: BOLTZ-2 K392 BIAS CHECK")
    print(SEP)
    print("%-25s %8s %8s %8s %10s" % ("Control", "K392", "N392", "Delta", "Leans"))
    print(DASH)

    bias_deltas = []
    for name in BIAS_CONTROLS:
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        bias_deltas.append(d)
        lean = "K392" if d > 0.02 else "N392" if d < -0.02 else "neutral"
        print("%-25s %8.3f %8.3f %+8.3f %10s" % (name, k, n, d, lean))

    avg_bias = mean(bias_deltas)
    k392_count = sum(1 for d in bias_deltas if d > 0.02)
    n392_count = sum(1 for d in bias_deltas if d < -0.02)

    print(DASH)
    print("Average control delta: %+.3f" % avg_bias)
    print("Controls leaning K392: %d/%d" % (k392_count, len(bias_deltas)))
    print("Controls leaning N392: %d/%d" % (n392_count, len(bias_deltas)))

    if avg_bias > 0.02:
        print("\n>>> SYSTEMATIC K392 BIAS DETECTED (avg delta %+.3f)" % avg_bias)
        print(">>> Peptide selectivity deltas should be interpreted RELATIVE to this baseline")
        print(">>> Corrected pep_glu_long_01 delta: %+.3f (raw +0.097 minus bias %+.3f)" % (0.097 - avg_bias, avg_bias))
        print(">>> Corrected hybrid_E_VKLLLL delta: %+.3f (raw +0.075 minus bias %+.3f)" % (0.075 - avg_bias, avg_bias))
    elif avg_bias < -0.02:
        print("\n>>> SYSTEMATIC N392 BIAS DETECTED — K392 selectivity signal is STRONGER than it looks")
    else:
        print("\n>>> NO SYSTEMATIC BIAS — selectivity deltas are interpretable at face value")

    print("\n" + SEP)
    print("PART 2: VKLLLL P1 CHEMISTRY SCAN")
    print(SEP)
    print("%-25s %4s %8s %8s %8s %10s" % ("Peptide", "P1", "K392", "N392", "Delta", "Leans"))
    print(DASH)

    for name in P1_SCAN:
        p1 = P1_SCAN[name][0]
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        lean = "K392" if d > 0.02 else "N392" if d < -0.02 else "neutral"
        # Bias-corrected delta
        corrected = d - avg_bias
        print("%-25s %4s %8.3f %8.3f %+8.3f %10s  (corrected: %+.3f)" % (name, p1, k, n, d, lean, corrected))

    # Save everything
    output = {
        "bias_controls": {name: scores[name] for name in BIAS_CONTROLS},
        "p1_scan": {name: scores[name] for name in P1_SCAN},
        "average_bias": avg_bias,
        "bias_interpretation": "K392_BIAS" if avg_bias > 0.02 else "N392_BIAS" if avg_bias < -0.02 else "NO_BIAS",
        "all_sequences": ALL_PEPTIDES,
    }
    out_path = os.path.join(OUT_BASE, "bias_check_results.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
