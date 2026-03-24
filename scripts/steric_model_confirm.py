"""
Steric model confirmation: P1 scan + KKY motif tests.

Three experiments in one run:
1. KILKLYSSKKY P1 scan (7 P1s x 2 variants = 14 predictions)
   - Ile, Val, Leu, Phe, Trp (hydrophobic/aromatic — predicted best by steric model)
   - Glu (charged — should underperform if steric model is right)
   - Ala (baseline)

2. KKY motif graft onto 3 Bucket A scaffolds (3 x 2 variants = 6 predictions)
   - Tests if KKY is a transferable selectivity module

3. KKY knockout on KILKLYSSKKY (1 x 2 variants = 2 predictions)
   - KILKLYSSKAA — is the KKY motif load-bearing?

Total: 22 predictions, ~20 min on RTX 4090.

Usage (on Vast.ai):
    python3 /workspace/steric_model_confirm.py
"""
import subprocess
import json
import os
import glob
from statistics import mean
from collections import defaultdict

OUT_BASE = "/workspace/results/steric_confirm"
SEED = 42
DIFFUSION_SAMPLES = 3  # 3 samples for more reliable scores

ERAP2_FULL = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
K392_CHANNEL = ERAP2_FULL[349:450]
N392_CHANNEL = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:450]

# ============================================================
# EXPERIMENT 1: KILKLYSSKKY P1 scan
# Original scaffold: KILKLYSSKKY (P1=K, but screened with P1=A)
# Backbone after P1: ILKLYSSKKY
# ============================================================
KILK_BACKBONE = "ILKLYSSKKY"  # 10 residues after P1

P1_SCAN = {
    "kilk_P1_I": {"seq": "I" + KILK_BACKBONE, "p1": "Ile", "category": "branched_hydrophobic"},
    "kilk_P1_V": {"seq": "V" + KILK_BACKBONE, "p1": "Val", "category": "branched_hydrophobic"},
    "kilk_P1_L": {"seq": "L" + KILK_BACKBONE, "p1": "Leu", "category": "linear_hydrophobic"},
    "kilk_P1_F": {"seq": "F" + KILK_BACKBONE, "p1": "Phe", "category": "aromatic"},
    "kilk_P1_W": {"seq": "W" + KILK_BACKBONE, "p1": "Trp", "category": "aromatic"},
    "kilk_P1_E": {"seq": "E" + KILK_BACKBONE, "p1": "Glu", "category": "charged"},
    "kilk_P1_A": {"seq": "A" + KILK_BACKBONE, "p1": "Ala", "category": "baseline"},
}

# ============================================================
# EXPERIMENT 2: KKY motif graft onto Bucket A scaffolds
# These are strong K392 binders but LOW selectivity
# Graft SSKKY onto their C-terminus (replace last 5 residues)
# ============================================================

# Bucket A scaffolds (K392 >= 0.7 but delta < 0.05)
# TSYYWSPLHK (10, K392=0.891, delta=+0.019) → TSYYWSSKKY
# SISLYIDRRLPGK (13, K392=0.874, delta=+0.015) → SISLYIDRRSSKKY (14)
# VTLKLYSIFHG (11, K392=0.870, delta=-0.005) → VTLKLYSSKKY

KKY_GRAFTS = {
    "graft_TSYYW": {"seq": "TSYYWSSKKY", "original": "TSYYWSPLHK", "note": "replaced SPLHK with SSKKY"},
    "graft_VTLKL": {"seq": "VTLKLYSSKKY", "original": "VTLKLYSIFHG", "note": "replaced SIFHG with SSKKY"},
    "graft_SSTPW": {"seq": "SSTPWSSKKY", "original": "SSTPWLSPISG", "note": "replaced LSPISG with SSKKY"},
}

# ============================================================
# EXPERIMENT 3: KKY knockout on KILKLYSSKKY
# Replace KKY with AAA — is the motif load-bearing?
# ============================================================
KKY_KNOCKOUT = {
    "kilk_KKY_knockout": {"seq": "AILKLYSSKAA", "note": "KKY->KAA, P1=A (same as screen baseline)"},
}

# Combine all experiments
ALL_EXPERIMENTS = {}
for name, info in P1_SCAN.items():
    ALL_EXPERIMENTS[name] = info["seq"]
for name, info in KKY_GRAFTS.items():
    ALL_EXPERIMENTS[name] = info["seq"]
for name, info in KKY_KNOCKOUT.items():
    ALL_EXPERIMENTS[name] = info["seq"]


def write_yaml(name, variant, pep_seq, channel_seq):
    yaml_dir = os.path.join(OUT_BASE, "yamls")
    os.makedirs(yaml_dir, exist_ok=True)
    path = os.path.join(yaml_dir, "%s_%s.yaml" % (name, variant))
    with open(path, "w") as f:
        f.write("version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: %s\n      msa: empty\n  - protein:\n      id: B\n      sequence: %s\n      msa: empty\n" % (channel_seq, pep_seq))
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
    total = len(ALL_EXPERIMENTS) * 2
    done = 0

    print("STERIC MODEL CONFIRMATION")
    print("=" * 70)
    print("Exp 1: KILKLYSSKKY P1 scan (7 P1s)")
    print("Exp 2: KKY motif graft (3 Bucket A scaffolds)")
    print("Exp 3: KKY knockout (KILKLYSSKAA)")
    print("Total: %d predictions (3 diffusion samples each)" % total)
    print()

    scores = {}
    for name, seq in ALL_EXPERIMENTS.items():
        scores[name] = {}
        for variant, channel in [("K392", K392_CHANNEL), ("N392", N392_CHANNEL)]:
            done += 1
            tag = "%s_%s" % (name, variant)
            out_dir = os.path.join(OUT_BASE, tag)
            print("[%d/%d] %s (%s)..." % (done, total, tag, seq), end=" ", flush=True)
            yaml_path = write_yaml(name, variant, seq, channel)
            ok = run_boltz(yaml_path, out_dir)
            if ok:
                iptm = extract_iptm(out_dir)
                scores[name][variant] = iptm
                print("ipTM=%.3f" % iptm)
            else:
                scores[name][variant] = 0
                print("FAILED")

    # ============================================================
    # ANALYSIS
    # ============================================================
    SEP = "=" * 80
    DASH = "-" * 80

    # Experiment 1: P1 scan
    print("\n" + SEP)
    print("EXPERIMENT 1: KILKLYSSKKY P1 SCAN")
    print("Steric model predicts: branched hydrophobic > aromatic > charged > baseline")
    print(SEP)
    print("%-18s %10s %8s %8s %8s %20s" % ("Peptide", "P1", "K392", "N392", "Delta", "Category"))
    print(DASH)

    p1_results = []
    for name in P1_SCAN:
        info = P1_SCAN[name]
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        print("%-18s %10s %8.3f %8.3f %+8.3f %20s" % (name, info["p1"], k, n, d, info["category"]))
        p1_results.append({"name": name, "p1": info["p1"], "category": info["category"], "K392": k, "N392": n, "delta": d})

    # Check steric model prediction
    branched = [r for r in p1_results if r["category"] == "branched_hydrophobic"]
    aromatic = [r for r in p1_results if r["category"] == "aromatic"]
    charged = [r for r in p1_results if r["category"] == "charged"]
    baseline = [r for r in p1_results if r["category"] == "baseline"]

    avg_branched = mean([r["delta"] for r in branched]) if branched else 0
    avg_aromatic = mean([r["delta"] for r in aromatic]) if aromatic else 0
    avg_charged = mean([r["delta"] for r in charged]) if charged else 0
    avg_baseline = mean([r["delta"] for r in baseline]) if baseline else 0

    print("\nCategory averages:")
    print("  Branched hydrophobic (I,V): %+.3f" % avg_branched)
    print("  Aromatic (F,W):             %+.3f" % avg_aromatic)
    print("  Charged (E):                %+.3f" % avg_charged)
    print("  Baseline (A):               %+.3f" % avg_baseline)

    if avg_branched > avg_charged and avg_branched > avg_baseline:
        print("\n>>> STERIC MODEL CONFIRMED: branched hydrophobic > charged for selectivity")
    elif avg_charged > avg_branched:
        print("\n>>> STERIC MODEL REJECTED: charged > branched — salt bridge may matter more here")
    else:
        print("\n>>> INCONCLUSIVE: differences too small to call")

    # Experiment 2: KKY graft
    print("\n" + SEP)
    print("EXPERIMENT 2: KKY MOTIF GRAFT")
    print("Does grafting SSKKY onto Bucket A scaffolds boost selectivity?")
    print(SEP)
    print("%-18s %22s %8s %8s %8s %12s" % ("Graft", "Sequence", "K392", "N392", "Delta", "Original D"))
    print(DASH)

    # Original Bucket A deltas for comparison
    original_deltas = {"graft_TSYYW": +0.019, "graft_VTLKL": -0.005, "graft_SSTPW": +0.028}
    for name in KKY_GRAFTS:
        info = KKY_GRAFTS[name]
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        orig = original_deltas.get(name, 0)
        improvement = d - orig
        print("%-18s %22s %8.3f %8.3f %+8.3f %+12.3f  (%+.3f change)" % (name, info["seq"], k, n, d, orig, improvement))

    # Experiment 3: KKY knockout
    print("\n" + SEP)
    print("EXPERIMENT 3: KKY KNOCKOUT")
    print("KILKLYSSKKY (Ala P1) = 0.892 K392, +0.427 delta")
    print("KILKLYSSKAA (KKY->KAA) = ?")
    print(SEP)

    for name in KKY_KNOCKOUT:
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        k_drop = 0.892 - k
        d_drop = 0.427 - d
        print("%-18s  K392=%.3f  N392=%.3f  delta=%+.3f" % (name, k, n, d))
        print("  K392 drop: %.3f (%.1f%%)" % (k_drop, k_drop / 0.892 * 100))
        print("  Selectivity drop: %.3f (%.1f%%)" % (d_drop, d_drop / 0.427 * 100))
        if d_drop > 0.2:
            print("  >>> KKY MOTIF IS LOAD-BEARING — selectivity drops >50%% without it")
        elif d_drop > 0.1:
            print("  >>> KKY contributes significantly but scaffold also matters")
        else:
            print("  >>> KKY is NOT the primary selectivity driver — scaffold geometry dominates")

    # Save
    output = {
        "p1_scan": p1_results,
        "kky_grafts": {name: {"seq": KKY_GRAFTS[name]["seq"], "K392": scores[name].get("K392", 0), "N392": scores[name].get("N392", 0), "delta": scores[name].get("K392", 0) - scores[name].get("N392", 0)} for name in KKY_GRAFTS},
        "kky_knockout": {name: {"seq": KKY_KNOCKOUT[name]["seq"], "K392": scores[name].get("K392", 0), "N392": scores[name].get("N392", 0), "delta": scores[name].get("K392", 0) - scores[name].get("N392", 0)} for name in KKY_KNOCKOUT},
        "steric_model": {
            "avg_branched_delta": avg_branched,
            "avg_aromatic_delta": avg_aromatic,
            "avg_charged_delta": avg_charged,
            "avg_baseline_delta": avg_baseline,
        },
    }
    out_path = os.path.join(OUT_BASE, "steric_confirm_results.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print("\nSaved to %s" % out_path)


if __name__ == "__main__":
    main()
