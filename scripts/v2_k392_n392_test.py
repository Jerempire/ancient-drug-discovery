"""
V2 K392/N392 Variant Selectivity Test

Tests all 4 V2 synthesis candidates against ERAP2-N392 to fill in the
missing allele selectivity data. K392 data already exists in v2_reference.md.

Runs on Vast.ai with Boltz-2 (3 diffusion samples, seed 42).

Usage: python3 /workspace/scripts/v2_k392_n392_test.py
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import subprocess
import glob
from datetime import datetime

RESULTS_DIR = "/workspace/results/v2_variant_selectivity"
YAML_DIR = RESULTS_DIR + "/yamls"

# V2 synthesis candidates (from v2_reference.md)
CONSTRUCTS = {
    "Y87A_Y89A": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
        "length": 92,
        "k392_iptm": 0.748,
        "erap1_iptm": 0.112,
        "irap_iptm": 0.186,
    },
    "trim_c5": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "length": 92,
        "k392_iptm": 0.691,
        "erap1_iptm": 0.362,
        "irap_iptm": 0.601,
    },
    "wt_97aa": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "length": 97,
        "k392_iptm": 0.475,
        "erap1_iptm": 0.182,
        "irap_iptm": 0.491,
    },
    "ko_aromatics": {
        "sequence": "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
        "length": 97,
        "k392_iptm": 0.315,
        "erap1_iptm": None,
        "irap_iptm": None,
    },
}

# ERAP2 N392 channel crop (151aa, residues 350-500 with K392N)
ERAP2_FULL = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNG"
    "ERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITN"
    "ATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEG"
    "FYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKV"
    "KTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYA"
    "LQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSAS"
    "DKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFL"
    "NVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGI"
    "IQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMM"
    "TTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSS"
    "NVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNIS"
    "DISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQ"
    "WMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEK"
    "LLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRM"
    "IISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWL"
    "MVNT"
)
# K392 crop (original)
ERAP2_K392 = ERAP2_FULL[349:500]
# N392 crop (K->N at position 392, index 391 in full sequence)
ERAP2_N392 = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:500]

DIFFUSION_SAMPLES = 3
SEED = 42


def main():
    os.makedirs(YAML_DIR, exist_ok=True)

    print("=" * 70)
    print("V2 K392/N392 VARIANT SELECTIVITY TEST")
    print("=" * 70)
    print(f"Constructs: {len(CONSTRUCTS)}")
    print(f"Targets: ERAP2-K392, ERAP2-N392 (3 diffusion samples each)")
    print(f"Total predictions: {len(CONSTRUCTS) * 2}")
    print()

    # Generate YAMLs for both K392 and N392
    for cname, cdata in CONSTRUCTS.items():
        for variant, crop in [("k392", ERAP2_K392), ("n392", ERAP2_N392)]:
            name = f"{cname}_{variant}"
            yaml_path = f"{YAML_DIR}/{name}.yaml"
            with open(yaml_path, "w") as f:
                f.write(f"version: 1\nsequences:\n"
                        f"  - protein:\n      id: A\n      sequence: {crop}\n      msa: empty\n"
                        f"  - protein:\n      id: B\n      sequence: {cdata['sequence']}\n      msa: empty\n")

    # Run Boltz-2
    all_results = {}
    yamls = sorted(glob.glob(f"{YAML_DIR}/*.yaml"))
    total = len(yamls)

    for yi, yaml_file in enumerate(yamls):
        name = os.path.basename(yaml_file).replace(".yaml", "")
        out_dir = f"{RESULTS_DIR}/{name}"

        # Skip if already done
        existing = glob.glob(f"{out_dir}/**/confidence_*.json", recursive=True)
        if existing and len(existing) >= DIFFUSION_SAMPLES:
            iptms = []
            for cf in existing:
                with open(cf) as f:
                    c = json.load(f)
                iptms.append(c.get("iptm", c.get("i_ptm", 0)))
            avg = sum(iptms) / len(iptms)
            print(f"  [{yi+1}/{total}] {name}: CACHED avg={avg:.3f} samples={[round(x,3) for x in iptms]}")
            all_results[name] = {"avg": avg, "samples": iptms}
            continue

        print(f"  [{yi+1}/{total}] {name}...", end=" ", flush=True)
        result = subprocess.run(
            ["boltz", "predict", yaml_file, "--out_dir", out_dir,
             "--diffusion_samples", str(DIFFUSION_SAMPLES), "--seed", str(SEED)],
            capture_output=True, text=True, timeout=600
        )

        if result.returncode == 0:
            conf_files = sorted(glob.glob(f"{out_dir}/**/confidence_*.json", recursive=True))
            iptms = []
            for cf in conf_files:
                with open(cf) as f:
                    c = json.load(f)
                iptms.append(c.get("iptm", c.get("i_ptm", 0)))
            avg = sum(iptms) / len(iptms) if iptms else 0
            print(f"avg={avg:.3f} samples={[round(x,3) for x in iptms]}")
            all_results[name] = {"avg": avg, "samples": iptms}
        else:
            print("FAILED")
            print(f"    stderr: {result.stderr[-200:]}")
            all_results[name] = {"avg": 0, "samples": []}

    # Results table
    print()
    print("=" * 70)
    print("RESULTS: V2 VARIANT SELECTIVITY")
    print("=" * 70)
    print()
    print(f"{'Construct':<16} {'K392':>8} {'N392':>8} {'Delta':>8} {'ERAP1':>8} {'IRAP':>8}")
    print("-" * 60)

    summary = []
    for cname, cdata in CONSTRUCTS.items():
        k_key = f"{cname}_k392"
        n_key = f"{cname}_n392"
        k_avg = all_results.get(k_key, {}).get("avg", 0)
        n_avg = all_results.get(n_key, {}).get("avg", 0)
        delta = n_avg - k_avg if k_avg and n_avg else 0

        # Use existing K392/ERAP1/IRAP data from v2_reference.md if our fresh K392 is missing
        k_ref = cdata.get("k392_iptm", 0) or 0
        e1 = cdata.get("erap1_iptm") or "-"
        ir = cdata.get("irap_iptm") or "-"

        e1_str = f"{e1:.3f}" if isinstance(e1, (int, float)) else e1
        ir_str = f"{ir:.3f}" if isinstance(ir, (int, float)) else ir

        print(f"{cname:<16} {k_avg:>8.3f} {n_avg:>8.3f} {delta:>+8.3f} {e1_str:>8} {ir_str:>8}")

        summary.append({
            "construct": cname,
            "k392_fresh": k_avg,
            "k392_fresh_samples": all_results.get(k_key, {}).get("samples", []),
            "n392_fresh": n_avg,
            "n392_fresh_samples": all_results.get(n_key, {}).get("samples", []),
            "delta_n392_k392": delta,
            "k392_reference": k_ref,
            "erap1_reference": cdata.get("erap1_iptm"),
            "irap_reference": cdata.get("irap_iptm"),
        })

    print()
    print("Note: K392 values are FRESH 3-sample runs (may differ from v2_reference.md")
    print("which used different conditions). Delta is computed from matched conditions.")

    # Save
    output = {
        "test": "V2 K392/N392 Variant Selectivity",
        "date": datetime.now().isoformat(),
        "diffusion_samples": DIFFUSION_SAMPLES,
        "seed": SEED,
        "msa": "empty",
        "results": summary,
        "raw_scores": {k: v for k, v in all_results.items()},
    }
    out_path = f"{RESULTS_DIR}/v2_variant_selectivity_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
