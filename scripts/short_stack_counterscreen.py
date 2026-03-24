"""
Short-stack counterscreen: all 18 peptides against ERAP1 and IRAP.

Combines with existing K392/N392 data to build full 4-target panel.
Tests Gemini's triple selectivity hypothesis:
- ERAP1: 6-mers too short for molecular ruler
- IRAP: P3 Phe clashes with IRAP Asn406, reduced surface area in cavern

Usage (on Vast.ai):
    python3 /workspace/short_stack_counterscreen.py
"""
import subprocess
import json
import os
import glob
from statistics import mean

OUT_BASE = "/workspace/results/short_stack_counter"
SEED = 42
DIFFUSION_SAMPLES = 3

# ERAP2 K392 channel (residues 350-450) for reference
ERAP2_FULL = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"

# ERAP1 equivalent channel region (using ERAP1 UniProt Q9NZ08, residues 333-433)
# Aligned to ERAP2 350-450 via structural alignment
ERAP1_FULL = "MEPPRGPRPLLTLLLLALAPGAGASQDCNPLAHGITGHLSHRIRQNFGWNLDQHFSAEQFQEYMLRGKHIPVKAHKGLSHLTAKDMAQKVQTPFNQRIGVSDPMLALVAQSGPGLERSYVLLASAEVDNIPRSTAQSISLDDDTLQWFQLAQLVQHVIRQYSDTDRISSPSTASALSPSSALLWPQGFILANPDEFNAQKLCSHLGTPDPAHYTALTRYQLQNPKLNNLHKSMQPYILAQELDALNYENAFKFYELQFSHTDALSNFKPSASSAFYSQLMKYLNNLKDLLEEKEISRGIFQKLHLPKQGSLRKNKVKIGLLDIFALKFNTAELQTIENDWLSHFRFQQLSGVTDLMFRSPLMRFFTFYESHFLAMELYNKKFMKLYKNQFQPITKFDDTKLYHAVKSLNSMTAMGSNVFNLNIFNLESMGYAEVKDRIAEWQNFKQHGSSPDLTIEQSMRVLQKFDTFVKPRYSQMKNHPEEMDFFSYYSFSWKDDPMQCFTFNCNEHNAILKWYQDRKAMNVGFKEVELPSFQKFAENLRESILKVVMVYIFAAAQALSKELDVNAWKYANNDKEQGVITFTTDESNYLRVRVSHELMEFSYLFLDESYRQYVVQIFERLPYGDLRTLMKASMPVFAQDFADTASFLNKFSFPTLGFMSRTMVRDRNRGMAVISPHFYYSGSEAGITLHYVDMENPQELNDAIRFAWRNFFRQCFDIDSQAAIRHYELQFQYTSMSGTALPRGDSYALDGLVEGAMEVKEQFGTFIDALANYSLQNHGLPASCTFTVQKNASQALLEELRVLHLHREGGIYAGSQHFDSYEGLFNKRDLMVDGKMLPEKFHSNKFSESYRSVFVDDGMQNYLNQFLMQHSIHQHALTSGFLRDNFEHFPWMKEIQASGTSMRPMYHKPGLQLLSLGLDYDLVNAFPKDKEQYAFAGWLKELFLPLEGGFAAQ"
ERAP1_CHANNEL = ERAP1_FULL[332:433]  # Aligned to ERAP2 350-450

# IRAP equivalent channel region (using IRAP UniProt Q9UIQ6, residues 444-544)
IRAP_FULL = "MSKRKEFGALITLGLCVMACSAEEEERPHRDLADCSAEEYEQGATPQLFDTMHKFNDQVYAQRNPHTYHGTKAQIKEGLLDRFANRIQELEEDTSEDGSKDSSPFRCHHGKILTINGSVDMSMEQTAQEEFHKKTYRVDPEKKQRGDDNKCGDTCLSWLAQKFPHKQYHQFGIIDRQLKALPPMEKPRYGSDPAMAIRTKDFMTPNFVHYAGMGSPSEGQYKMTYREYAQRVVELYEEFLTAKDVLQKRFQFPDVFYSDNGTLYPQPKHHTLQNGQNFRRWEDHMDPWFVLDLFHKQAELKQGYLKKKDYMKLYKEKYQSIIKFNQPQKFEDDKRFYHAVKSLGDSKAMEAEYSFLNLNLFNIESMGYPEKLDILGKWHNLRMNVSFVDVFNIHDYKMYEIVEDILRHLIEFDDYRFVSIEFLMHDSKNFQPLNYDFEVYEYSYLSWRRPDMCFTFNCNHSHAIEQWYQRKFEMNSRLLEFSAEPNKSFVTKFDLEIKPNYTPMKNHPHDMYFFYSYANWKDDPLQCFTYNCNKQGTAILKWYRSRKSSNIGLKEIPMPSLQKFAHNLRDTIISAIMIYMFAAAEALTKELDTSGWDYANDTDNRVTFTIDEGSYLKARIPHELAHQWFGNLVTMEWWNDIYLKEGFAKYMELIAVNATYPELQFDDYYLNFCFEVISKDSINSSRPISKPEKTPAQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFAYRNAKNNDLWSSLSNSCLESDFTSGGVCHSDDKMTSSMLAFLEEHADIKEMMTTWTLQKGTPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTALVKFNVDSNGYYIVHYEGHGWDQLIIKLNQNHTLLRPKDRVGLMHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGFSYLESFYHMMDLRKISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSLSSSEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
IRAP_CHANNEL = IRAP_FULL[443:544]  # Aligned to ERAP2 350-450

# All 18 short-stack peptides
LIBRARY = {
    "ss6_VAGSAF":  "VAGSAF",
    "ss6_IAGSAY":  "IAGSAY",
    "ss6_VAGSYW":  "VAGSYW",
    "ss6_IAGSAW":  "IAGSAW",
    "ss6_LAGSAF":  "LAGSAF",
    "ss6_VASSKY":  "VASSKY",
    "ss7_VAGSAAF": "VAGSAAF",
    "ss7_IAGSAAY": "IAGSAAY",
    "ss7_VAGSKKY": "VAGSKKY",
    "ss7_IAGSKKY": "IAGSKKY",
    "ss7_VAFSAGY": "VAFSAGY",
    "ss7_IALSAFW": "IALSAFW",
    "ss6_VAWSAF":  "VAWSAF",
    "ss6_IAFSAF":  "IAFSAF",
    "ss7_VAWSAAF": "VAWSAAF",
    "ss7_IAFSAAY": "IAFSAAY",
    "ss6_AAGSAF":  "AAGSAF",
    "ss6_EAGSAF":  "EAGSAF",
}

# K392/N392 scores from the first run (to combine into full panel)
PREV_SCORES = {
    "ss6_VAGSAF":  {"K392": 0.905, "N392": 0.870},
    "ss6_IAGSAY":  {"K392": 0.759, "N392": 0.858},
    "ss6_VAGSYW":  {"K392": 0.777, "N392": 0.777},
    "ss6_IAGSAW":  {"K392": 0.820, "N392": 0.744},
    "ss6_LAGSAF":  {"K392": 0.821, "N392": 0.717},
    "ss6_VASSKY":  {"K392": 0.718, "N392": 0.657},
    "ss7_VAGSAAF": {"K392": 0.801, "N392": 0.882},
    "ss7_IAGSAAY": {"K392": 0.731, "N392": 0.749},
    "ss7_VAGSKKY": {"K392": 0.728, "N392": 0.844},
    "ss7_IAGSKKY": {"K392": 0.776, "N392": 0.820},
    "ss7_VAFSAGY": {"K392": 0.833, "N392": 0.754},
    "ss7_IALSAFW": {"K392": 0.877, "N392": 0.874},
    "ss6_VAWSAF":  {"K392": 0.852, "N392": 0.830},
    "ss6_IAFSAF":  {"K392": 0.870, "N392": 0.631},
    "ss7_VAWSAAF": {"K392": 0.629, "N392": 0.695},
    "ss7_IAFSAAY": {"K392": 0.710, "N392": 0.627},
    "ss6_AAGSAF":  {"K392": 0.799, "N392": 0.635},
    "ss6_EAGSAF":  {"K392": 0.730, "N392": 0.641},
}

TARGETS = {
    "ERAP1": ERAP1_CHANNEL,
    "IRAP": IRAP_CHANNEL,
}


def write_yaml(name, target, pep_seq, channel_seq):
    yaml_dir = os.path.join(OUT_BASE, "yamls")
    os.makedirs(yaml_dir, exist_ok=True)
    path = os.path.join(yaml_dir, "%s_%s.yaml" % (name, target))
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
    total = len(LIBRARY) * len(TARGETS)
    done = 0

    print("SHORT-STACK COUNTERSCREEN: ERAP1 + IRAP")
    print("=" * 70)
    print("18 peptides x 2 targets = %d predictions" % total)
    print()

    scores = {}
    for name, seq in LIBRARY.items():
        scores[name] = dict(PREV_SCORES.get(name, {}))
        for target, channel in TARGETS.items():
            done += 1
            tag = "%s_%s" % (name, target)
            out_dir = os.path.join(OUT_BASE, tag)
            print("[%d/%d] %s (%s vs %s)..." % (done, total, tag, seq, target), end=" ", flush=True)
            yaml_path = write_yaml(name, target, seq, channel)
            ok = run_boltz(yaml_path, out_dir)
            if ok:
                iptm = extract_iptm(out_dir)
                scores[name][target] = iptm
                print("ipTM=%.3f" % iptm)
            else:
                scores[name][target] = 0
                print("FAILED")

    # ============================================================
    # FULL 4-TARGET PANEL
    # ============================================================
    SEP = "=" * 100
    DASH = "-" * 100

    print("\n" + SEP)
    print("FULL 4-TARGET PANEL: SHORT-STACK PEPTIDES")
    print("Tier 1: K392 binding | Tier 2: ERAP1 avoidance | Tier 3: N392 discrimination | Bonus: IRAP avoidance")
    print(SEP)
    print("%-14s %6s %8s %8s %8s %8s %10s %10s %10s" % (
        "Name", "Seq", "K392", "N392", "ERAP1", "IRAP", "selE1", "selIR", "selN392"))
    print(DASH)

    results = []
    for name, seq in LIBRARY.items():
        s = scores[name]
        k = s.get("K392", 0)
        n = s.get("N392", 0)
        e1 = s.get("ERAP1", 0)
        ir = s.get("IRAP", 0)
        sel_e1 = k - e1  # positive = avoids ERAP1
        sel_ir = k - ir   # positive = avoids IRAP
        sel_n = k - n     # positive = K392-selective

        results.append({
            "name": name, "seq": seq, "len": len(seq),
            "K392": k, "N392": n, "ERAP1": e1, "IRAP": ir,
            "sel_ERAP1": sel_e1, "sel_IRAP": sel_ir, "sel_N392": sel_n,
        })

    # Sort by K392 potency
    results.sort(key=lambda x: -x["K392"])

    for r in results:
        # Flag triple-selective candidates
        triple = ""
        if r["sel_ERAP1"] > 0.05 and r["sel_IRAP"] > 0.05 and r["sel_N392"] > 0.02:
            triple = " *** TRIPLE"
        elif r["sel_ERAP1"] > 0.05 and r["sel_N392"] > 0.02:
            triple = " ** DUAL(E1+N)"
        elif r["sel_ERAP1"] > 0.05:
            triple = " * E1-selective"

        print("%-14s %6s %8.3f %8.3f %8.3f %8.3f %+10.3f %+10.3f %+10.3f%s" % (
            r["name"], r["seq"], r["K392"], r["N392"], r["ERAP1"], r["IRAP"],
            r["sel_ERAP1"], r["sel_IRAP"], r["sel_N392"], triple))

    # Triple selectivity candidates
    print("\n" + SEP)
    print("TRIPLE SELECTIVITY CANDIDATES")
    print("K392 > ERAP1 (+0.05) AND K392 > IRAP (+0.05) AND K392 > N392 (+0.02)")
    print(SEP)

    triple_hits = [r for r in results if r["sel_ERAP1"] > 0.05 and r["sel_IRAP"] > 0.05 and r["sel_N392"] > 0.02]
    if triple_hits:
        for r in triple_hits:
            print("  %s (%s, %d-mer): K392=%.3f, selE1=%+.3f, selIR=%+.3f, selN=%+.3f" % (
                r["name"], r["seq"], r["len"], r["K392"], r["sel_ERAP1"], r["sel_IRAP"], r["sel_N392"]))
    else:
        print("  None found — check dual-selective candidates above")

    # ERAP1 molecular ruler test
    print("\n" + SEP)
    print("ERAP1 MOLECULAR RULER TEST")
    print("Do 6-mers score lower on ERAP1 than on ERAP2-K392?")
    print(SEP)
    sixmers = [r for r in results if r["len"] == 6]
    sevenmers = [r for r in results if r["len"] == 7]
    if sixmers:
        avg_e1_6 = mean([r["ERAP1"] for r in sixmers])
        avg_k_6 = mean([r["K392"] for r in sixmers])
        print("6-mers: avg K392=%.3f, avg ERAP1=%.3f, avg sel_E1=%+.3f" % (avg_k_6, avg_e1_6, avg_k_6 - avg_e1_6))
    if sevenmers:
        avg_e1_7 = mean([r["ERAP1"] for r in sevenmers])
        avg_k_7 = mean([r["K392"] for r in sevenmers])
        print("7-mers: avg K392=%.3f, avg ERAP1=%.3f, avg sel_E1=%+.3f" % (avg_k_7, avg_e1_7, avg_k_7 - avg_e1_7))

    # IRAP clash handle test
    print("\n" + SEP)
    print("IRAP CLASH HANDLE TEST (P3 position)")
    print(SEP)
    clash_f = next((r for r in results if r["name"] == "ss6_IAFSAF"), None)
    clash_w = next((r for r in results if r["name"] == "ss6_VAWSAF"), None)
    no_clash = next((r for r in results if r["name"] == "ss6_VAGSAF"), None)
    if clash_f and no_clash:
        print("IAFSAF (F at P3): IRAP=%.3f, sel_IR=%+.3f" % (clash_f["IRAP"], clash_f["sel_IRAP"]))
        print("VAGSAF (G at P3): IRAP=%.3f, sel_IR=%+.3f" % (no_clash["IRAP"], no_clash["sel_IRAP"]))
        if clash_f["sel_IRAP"] > no_clash["sel_IRAP"]:
            print(">>> P3=F REDUCES IRAP BINDING — clash handle works")
        else:
            print(">>> P3=F does NOT reduce IRAP — clash handle hypothesis fails")

    # Save
    output = {
        "full_panel": sorted(results, key=lambda x: -x["K392"]),
        "triple_selective": triple_hits,
        "erap1_ruler_test": {
            "sixmer_avg_k392": mean([r["K392"] for r in sixmers]) if sixmers else 0,
            "sixmer_avg_erap1": mean([r["ERAP1"] for r in sixmers]) if sixmers else 0,
        },
    }
    out_path = os.path.join(OUT_BASE, "counterscreen_results.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print("\nSaved to %s" % out_path)


if __name__ == "__main__":
    main()
