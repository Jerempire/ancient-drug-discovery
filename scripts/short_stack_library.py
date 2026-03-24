"""
V4.2 Short-Stack Library: 6-7mer peptides for triple selectivity.

Hypothesis (Gemini):
- 6-7mers evade ERAP1 molecular ruler (needs >8mer to close)
- Branched hydrophobic P1 packs against K392 aliphatic stem (steric model)
- Short length reduces IRAP "stickiness" (marble-in-bathtub)
- IRAP Tyr392 clashes with branched P1 designed for Lys392 stem
- C-terminal aromatics stack at channel exit

Design rules:
  P1: V, I, L (branched/hydrophobic — steric K392 packing)
  Mid: A, G, S (small, flexible — sit deep in S1' pocket)
  C-term: F, Y, W (aromatic stacking at channel exit)
  Length: 6 or 7 residues (below ERAP1 8-residue floor)

Also tests:
  - KKY motif on short scaffold (does it still work at 7mer?)
  - Position 406 clash handle (bulky P3/P4 for IRAP exclusion)

Screen against: K392, N392 channels
Future: ERAP1, IRAP counterscreen on winners

Usage (on Vast.ai, chain after steric_model_confirm.py):
    python3 /workspace/short_stack_library.py
"""
import subprocess
import json
import os
import glob
from statistics import mean

OUT_BASE = "/workspace/results/short_stack"
SEED = 42
DIFFUSION_SAMPLES = 3

ERAP2_FULL = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
K392_CHANNEL = ERAP2_FULL[349:450]
N392_CHANNEL = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:450]

# ============================================================
# SHORT-STACK LIBRARY
# ============================================================
# Design: P1[branched] + mid[small] + C-term[aromatic]
# 6-mers and 7-mers

LIBRARY = {
    # --- 6-mers: P1-mid-mid-mid-mid-Cterm ---
    # Core steric designs
    "ss6_VAGSAF": {"seq": "VAGSAF", "len": 6, "logic": "V-P1 + small mid + F-Cterm"},
    "ss6_IAGSAY": {"seq": "IAGSAY", "len": 6, "logic": "I-P1 + small mid + Y-Cterm"},
    "ss6_VAGSYW": {"seq": "VAGSYW", "len": 6, "logic": "V-P1 + aromatic pair Cterm"},
    "ss6_IAGSAW": {"seq": "IAGSAW", "len": 6, "logic": "I-P1 + W-Cterm (large aromatic)"},
    "ss6_LAGSAF": {"seq": "LAGSAF", "len": 6, "logic": "L-P1 (linear) + F-Cterm"},
    "ss6_VASSKY": {"seq": "VASSKY", "len": 6, "logic": "V-P1 + mini KY motif (cation-pi)"},

    # --- 7-mers: one extra mid residue ---
    "ss7_VAGSAAF": {"seq": "VAGSAAF", "len": 7, "logic": "V-P1 + extended small mid + F-Cterm"},
    "ss7_IAGSAAY": {"seq": "IAGSAAY", "len": 7, "logic": "I-P1 + extended + Y-Cterm"},
    "ss7_VAGSKKY": {"seq": "VAGSKKY", "len": 7, "logic": "V-P1 + KKY motif on 7-mer"},
    "ss7_IAGSKKY": {"seq": "IAGSKKY", "len": 7, "logic": "I-P1 + KKY motif on 7-mer"},
    "ss7_VAFSAGY": {"seq": "VAFSAGY", "len": 7, "logic": "V-P1 + F-mid aromatic + Y-Cterm"},
    "ss7_IALSAFW": {"seq": "IALSAFW", "len": 7, "logic": "I-P1 + L-mid + FW aromatic pair"},

    # --- 406 clash handle designs ---
    # Position 406: ERAP2=Ala (tiny), IRAP=Asn (bulky)
    # Put bulky residue at P3/P4 to fit ERAP2 hole, clash with IRAP Asn
    "ss6_VAWSAF": {"seq": "VAWSAF", "len": 6, "logic": "W at P3 — fits ERAP2 Ala406 hole, clashes IRAP Asn406"},
    "ss6_IAFSAF": {"seq": "IAFSAF", "len": 6, "logic": "F at P3 — aromatic clash handle"},
    "ss7_VAWSAAF": {"seq": "VAWSAAF", "len": 7, "logic": "W at P3, 7-mer version"},
    "ss7_IAFSAAY": {"seq": "IAFSAAY", "len": 7, "logic": "F at P3, 7-mer, Y-Cterm"},

    # --- Controls ---
    "ss6_AAGSAF": {"seq": "AAGSAF", "len": 6, "logic": "A-P1 baseline (no steric P1)"},
    "ss6_EAGSAF": {"seq": "EAGSAF", "len": 6, "logic": "E-P1 charged (salt bridge vs steric test)"},
}


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
    total = len(LIBRARY) * 2
    done = 0

    print("V4.2 SHORT-STACK LIBRARY")
    print("=" * 70)
    print("6-7mer peptides for triple selectivity (ERAP2 > ERAP1 > IRAP)")
    print("Design: branched P1 + small mid + aromatic C-term")
    print("Total: %d predictions (%d peptides x 2 variants)" % (total, len(LIBRARY)))
    print()

    # Print library
    sixmers = {k: v for k, v in LIBRARY.items() if v["len"] == 6}
    sevenmers = {k: v for k, v in LIBRARY.items() if v["len"] == 7}
    print("6-mers (%d):" % len(sixmers))
    for name, info in sixmers.items():
        print("  %s: %s — %s" % (name, info["seq"], info["logic"]))
    print("7-mers (%d):" % len(sevenmers))
    for name, info in sevenmers.items():
        print("  %s: %s — %s" % (name, info["seq"], info["logic"]))
    print()

    scores = {}
    for name, info in LIBRARY.items():
        scores[name] = {}
        for variant, channel in [("K392", K392_CHANNEL), ("N392", N392_CHANNEL)]:
            done += 1
            tag = "%s_%s" % (name, variant)
            out_dir = os.path.join(OUT_BASE, tag)
            print("[%d/%d] %s (%s, %d-mer)..." % (done, total, tag, info["seq"], info["len"]), end=" ", flush=True)
            yaml_path = write_yaml(name, variant, info["seq"], channel)
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

    print("\n" + SEP)
    print("SHORT-STACK RESULTS — RANKED BY K392 SELECTIVITY")
    print(SEP)
    print("%-15s %8s %4s %8s %8s %8s %s" % ("Name", "Sequence", "Len", "K392", "N392", "Delta", "Logic"))
    print(DASH)

    results = []
    for name, info in LIBRARY.items():
        k = scores[name].get("K392", 0)
        n = scores[name].get("N392", 0)
        d = k - n
        results.append({"name": name, "seq": info["seq"], "len": info["len"], "K392": k, "N392": n, "delta": d, "logic": info["logic"]})

    for r in sorted(results, key=lambda x: -x["delta"]):
        sel = "K392" if r["delta"] > 0.02 else "N392" if r["delta"] < -0.02 else "neut"
        print("%-15s %8s %4d %8.3f %8.3f %+8.3f %s" % (r["name"], r["seq"], r["len"], r["K392"], r["N392"], r["delta"], r["logic"][:40]))

    # Key comparisons
    print("\n" + SEP)
    print("KEY COMPARISONS")
    print(SEP)

    # 1. Do 6-7mers show ANY K392 selectivity?
    k392_selective = [r for r in results if r["delta"] > 0.02]
    print("K392-selective short peptides: %d / %d" % (len(k392_selective), len(results)))
    if k392_selective:
        print(">>> SHORT-STACK HYPOTHESIS ALIVE — some 6-7mers show K392 preference")
        best = max(k392_selective, key=lambda x: x["delta"])
        print(">>> Best: %s (%s, delta=%+.3f)" % (best["name"], best["seq"], best["delta"]))
    else:
        print(">>> SHORT-STACK HYPOTHESIS FAILS — 6-7mers cannot reach K392 steric handle")

    # 2. Branched P1 vs charged P1 on 6-mers
    v_6mer = next((r for r in results if r["name"] == "ss6_VAGSAF"), None)
    e_6mer = next((r for r in results if r["name"] == "ss6_EAGSAF"), None)
    a_6mer = next((r for r in results if r["name"] == "ss6_AAGSAF"), None)
    if v_6mer and e_6mer and a_6mer:
        print("\n6-mer P1 comparison (steric vs charged vs baseline):")
        print("  Val P1: K392=%.3f, delta=%+.3f" % (v_6mer["K392"], v_6mer["delta"]))
        print("  Glu P1: K392=%.3f, delta=%+.3f" % (e_6mer["K392"], e_6mer["delta"]))
        print("  Ala P1: K392=%.3f, delta=%+.3f" % (a_6mer["K392"], a_6mer["delta"]))
        if v_6mer["delta"] > e_6mer["delta"]:
            print("  >>> Steric > charged confirmed at 6-mer length")

    # 3. KKY motif on 7-mer
    kky_7 = [r for r in results if "KKY" in r["logic"] or "kky" in r["name"].lower()]
    if kky_7:
        print("\nKKY motif on 7-mers:")
        for r in kky_7:
            print("  %s: K392=%.3f, delta=%+.3f" % (r["seq"], r["K392"], r["delta"]))

    # 4. 406 clash handle
    clash = [r for r in results if "406" in r["logic"] or "clash" in r["logic"]]
    if clash:
        print("\n406 clash handle designs:")
        for r in clash:
            print("  %s: K392=%.3f, delta=%+.3f — %s" % (r["seq"], r["K392"], r["delta"], r["logic"][:50]))

    # 5. Comparison to 11-mer leads
    print("\n" + SEP)
    print("COMPARISON TO 11-MER LEADS")
    print(SEP)
    print("KILKLYSSKKY (11-mer, Ala P1): K392=0.892, delta=+0.427")
    print("VKLLLLSIGK  (10-mer, Val P1): K392=0.801, delta=+0.137")
    if k392_selective:
        best = max(k392_selective, key=lambda x: x["K392"])
        print("Best short-stack:             K392=%.3f, delta=%+.3f (%s, %d-mer)" % (best["K392"], best["delta"], best["seq"], best["len"]))
        k_ratio = best["K392"] / 0.892
        print("Short-stack retains %.0f%% of 11-mer K392 binding" % (k_ratio * 100))
        if k_ratio > 0.7:
            print(">>> VIABLE: >70%% binding retained at 6-7mer — ERAP1 evasion + K392 selectivity")
        else:
            print(">>> MARGINAL: significant binding loss at short length")
    else:
        print("No K392-selective short peptides found — length floor is real")

    # Save
    output = {
        "library": {name: {"seq": LIBRARY[name]["seq"], "len": LIBRARY[name]["len"], "logic": LIBRARY[name]["logic"]} for name in LIBRARY},
        "results": sorted(results, key=lambda x: -x["delta"]),
        "k392_selective_count": len(k392_selective),
        "total_tested": len(results),
    }
    out_path = os.path.join(OUT_BASE, "short_stack_results.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print("\nSaved to %s" % out_path)


if __name__ == "__main__":
    main()
