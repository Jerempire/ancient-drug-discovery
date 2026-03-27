"""
V2 K392/N392 Variant Selectivity Test — WITH MSAs

Same as v2_k392_n392_test.py but uses Boltz-2's built-in MSA server for
the ERAP2 target chain. Binder remains msa=empty (de novo, no homologs).

This gives reference-quality absolute ipTM scores comparable to v2_reference.md.

Usage: python3 /workspace/scripts/v2_k392_n392_test_msa.py
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json, os, subprocess, glob
from datetime import datetime

RESULTS_DIR = "/workspace/results/v2_variant_msa"
YAML_DIR = RESULTS_DIR + "/yamls"

CONSTRUCTS = {
    "Y87A_Y89A": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
        "k392_ref": 0.748, "erap1_ref": 0.112, "irap_ref": 0.186,
    },
    "trim_c5": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
        "k392_ref": 0.691, "erap1_ref": 0.362, "irap_ref": 0.601,
    },
    "wt_97aa": {
        "sequence": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
        "k392_ref": 0.475, "erap1_ref": 0.182, "irap_ref": 0.491,
    },
    "ko_aromatics": {
        "sequence": "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
        "k392_ref": 0.315, "erap1_ref": None, "irap_ref": None,
    },
}

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
ERAP2_K392 = ERAP2_FULL[349:500]
ERAP2_N392 = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:500]

DIFFUSION_SAMPLES = 3
SEED = 42


def write_yaml(target_seq, binder_seq, path):
    """Write Boltz-2 YAML with MSA server for target, empty for binder."""
    with open(path, "w") as f:
        # Target chain A: use default MSA (Boltz-2 will query server)
        # Binder chain B: msa=empty (de novo sequence)
        f.write(f"version: 1\nsequences:\n"
                f"  - protein:\n      id: A\n      sequence: {target_seq}\n"
                f"  - protein:\n      id: B\n      sequence: {binder_seq}\n      msa: empty\n")


def main():
    os.makedirs(YAML_DIR, exist_ok=True)

    print("=" * 70)
    print("V2 K392/N392 VARIANT SELECTIVITY TEST (MSA-ENABLED)")
    print("=" * 70)
    print(f"Target ERAP2 chain: MSA from server (evolutionary covariance)")
    print(f"Binder chain: msa=empty (de novo)")
    print(f"Diffusion samples: {DIFFUSION_SAMPLES}, seed: {SEED}")
    print()

    # Generate YAMLs
    for cname, cdata in CONSTRUCTS.items():
        for variant, crop in [("k392", ERAP2_K392), ("n392", ERAP2_N392)]:
            name = f"{cname}_{variant}"
            write_yaml(crop, cdata["sequence"], f"{YAML_DIR}/{name}.yaml")

    # Run Boltz-2
    all_results = {}
    yamls = sorted(glob.glob(f"{YAML_DIR}/*.yaml"))
    total = len(yamls)

    for yi, yaml_file in enumerate(yamls):
        name = os.path.basename(yaml_file).replace(".yaml", "")
        out_dir = f"{RESULTS_DIR}/{name}"

        existing = glob.glob(f"{out_dir}/**/confidence_*.json", recursive=True)
        if existing and len(existing) >= DIFFUSION_SAMPLES:
            iptms = []
            for cf in existing:
                with open(cf) as f:
                    c = json.load(f)
                iptms.append(c.get("iptm", c.get("i_ptm", 0)))
            avg = sum(iptms) / len(iptms)
            print(f"  [{yi+1}/{total}] {name}: CACHED avg={avg:.3f}")
            all_results[name] = {"avg": avg, "samples": iptms}
            continue

        print(f"  [{yi+1}/{total}] {name}...", end=" ", flush=True)
        result = subprocess.run(
            ["boltz", "predict", yaml_file, "--out_dir", out_dir,
             "--diffusion_samples", str(DIFFUSION_SAMPLES), "--seed", str(SEED)],
            capture_output=True, text=True, timeout=900  # longer timeout for MSA fetch
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
            # Show error for debugging
            stderr = result.stderr[-300:] if result.stderr else "no stderr"
            print(f"    {stderr}")
            all_results[name] = {"avg": 0, "samples": [], "error": stderr}

    # Results
    print()
    print("=" * 70)
    print("RESULTS: V2 VARIANT SELECTIVITY (MSA-ENABLED)")
    print("=" * 70)
    print()
    print(f"{'Construct':<16} {'K392':>8} {'N392':>8} {'Delta':>8} {'K392ref':>8} {'ERAP1':>8} {'IRAP':>8}")
    print("-" * 70)

    summary = []
    for cname, cdata in CONSTRUCTS.items():
        k_key = f"{cname}_k392"
        n_key = f"{cname}_n392"
        k_avg = all_results.get(k_key, {}).get("avg", 0)
        n_avg = all_results.get(n_key, {}).get("avg", 0)
        delta = n_avg - k_avg if k_avg and n_avg else 0
        k_ref = cdata.get("k392_ref", 0) or 0
        e1 = cdata.get("erap1_ref")
        ir = cdata.get("irap_ref")
        e1_s = f"{e1:.3f}" if isinstance(e1, (int, float)) else "-"
        ir_s = f"{ir:.3f}" if isinstance(ir, (int, float)) else "-"

        print(f"{cname:<16} {k_avg:>8.3f} {n_avg:>8.3f} {delta:>+8.3f} {k_ref:>8.3f} {e1_s:>8} {ir_s:>8}")
        summary.append({
            "construct": cname, "k392_msa": k_avg, "n392_msa": n_avg,
            "delta": delta, "k392_ref": k_ref,
            "k392_samples": all_results.get(k_key, {}).get("samples", []),
            "n392_samples": all_results.get(n_key, {}).get("samples", []),
            "erap1_ref": e1, "irap_ref": ir,
        })

    print()
    print("K392ref = original v2_reference.md value (different run conditions)")

    output = {
        "test": "V2 K392/N392 Variant Selectivity (MSA-enabled)",
        "date": datetime.now().isoformat(),
        "diffusion_samples": DIFFUSION_SAMPLES, "seed": SEED,
        "msa_mode": "server for ERAP2, empty for binder",
        "results": summary,
    }
    out_path = f"{RESULTS_DIR}/v2_variant_msa_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
