"""
V2 Substrate Competition — Boltz-2 YAML Generator + Analysis

Generates Boltz-2 prediction configs for 3-chain (ERAP2 + V2 binder + substrate)
and 2-chain (ERAP2 + substrate) controls to test whether V2 functionally blocks
substrate access.

Usage:
  1. Generate YAMLs: python scripts/v2_substrate_competition.py generate
  2. After running Boltz-2 on Vast.ai, analyze: python scripts/v2_substrate_competition.py analyze

Expected Vast.ai workflow:
  - Upload YAMLs to /workspace/v2_competition/
  - Run: for f in /workspace/v2_competition/*.yaml; do
      boltz predict "$f" --out_dir /workspace/v2_competition/out --diffusion_samples 3 --recycling_steps 3
    done
  - Download /workspace/v2_competition/out/ and place in data/results/v2_inhibition/substrate_competition/
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import glob as glob_mod
from datetime import datetime, timezone

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(PROJECT, "data", "results", "v2_inhibition", "substrate_competition")
os.makedirs(OUTPUT_DIR, exist_ok=True)
YAML_DIR = os.path.join(OUTPUT_DIR, "boltz2_yamls")
os.makedirs(YAML_DIR, exist_ok=True)

# Full ERAP2 sequence (K392 variant = wild-type UniProt Q6P179)
ERAP2_SEQ = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPL"
    "HYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLV"
    "PEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESR"
    "HIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKL"
    "LDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTM"
    "EWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGA"
    "CILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVK"
    "EMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDT"
    "LDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYY"
    "LQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNH"
    "APCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLK"
    "LIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKL"
    "FFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
)

# V2 lead binder (Y87A_Y89A, 92aa)
V2_SEQ = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

# Substrate peptides
SUBSTRATES = {
    "VAGSAF":      {"seq": "VAGSAF",      "desc": "V4.3 lead 6-mer (triple-selective)"},
    "KILKLYSSKKY": {"seq": "KILKLYSSKKY", "desc": "N392 scaffold screen 11-mer"},
    "LLLLLLLLL":   {"seq": "LLLLLLLLL",   "desc": "Poly-Leu 9-mer (universal substrate)"},
}


def write_yaml(filepath, chains):
    """Write a Boltz-2 YAML config with given chains."""
    lines = ["version: 1", "sequences:"]
    for chain in chains:
        lines.append(f"  - protein:")
        lines.append(f"      id: {chain['id']}")
        lines.append(f"      sequence: {chain['seq']}")
        if chain.get("msa_empty", False):
            lines.append(f"      msa: empty")
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")


def generate():
    """Generate all Boltz-2 YAML configs."""
    print("=" * 60)
    print("GENERATING SUBSTRATE COMPETITION BOLTZ-2 YAMLS")
    print("=" * 60)

    configs = []

    for sub_name, sub_info in SUBSTRATES.items():
        sub_seq = sub_info["seq"]

        # Control: ERAP2 + substrate only (2 chains)
        ctrl_name = f"ctrl_{sub_name}"
        ctrl_path = os.path.join(YAML_DIR, f"{ctrl_name}.yaml")
        write_yaml(ctrl_path, [
            {"id": "A", "seq": ERAP2_SEQ, "msa_empty": False},
            {"id": "C", "seq": sub_seq,    "msa_empty": True},
        ])
        configs.append({"name": ctrl_name, "type": "control", "substrate": sub_name, "chains": 2})

        # Test: ERAP2 + V2 binder + substrate (3 chains)
        test_name = f"test_{sub_name}"
        test_path = os.path.join(YAML_DIR, f"{test_name}.yaml")
        write_yaml(test_path, [
            {"id": "A", "seq": ERAP2_SEQ, "msa_empty": False},
            {"id": "B", "seq": V2_SEQ,    "msa_empty": True},
            {"id": "C", "seq": sub_seq,   "msa_empty": True},
        ])
        configs.append({"name": test_name, "type": "test", "substrate": sub_name, "chains": 3})

    # Also: V2 binder alone (sanity check — should match prior ipTM ~0.748)
    sanity_name = "sanity_v2_erap2"
    sanity_path = os.path.join(YAML_DIR, f"{sanity_name}.yaml")
    write_yaml(sanity_path, [
        {"id": "A", "seq": ERAP2_SEQ, "msa_empty": False},
        {"id": "B", "seq": V2_SEQ,    "msa_empty": True},
    ])
    configs.append({"name": sanity_name, "type": "sanity", "substrate": None, "chains": 2})

    # Write manifest
    manifest = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "total_predictions": len(configs),
        "diffusion_samples": 3,
        "seed": 42,
        "configs": configs,
        "estimated_cost": f"~${len(configs) * 0.08:.2f} on RTX 4090 (~4.5 min each)",
        "total_chains_largest": f"ERAP2 (960aa) + V2 (92aa) + substrate (6-11aa) = ~1063aa",
    }

    manifest_path = os.path.join(OUTPUT_DIR, "competition_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\n  Generated {len(configs)} YAML configs in: {YAML_DIR}")
    for c in configs:
        chains_str = f"{c['chains']} chains"
        sub = c['substrate'] or 'N/A'
        print(f"    {c['name']:<25s}  {c['type']:<8s}  {chains_str}  substrate={sub}")

    print(f"\n  Manifest: {manifest_path}")
    print(f"\n  Estimated total cost: {manifest['estimated_cost']}")
    print(f"  Largest prediction: {manifest['total_chains_largest']}")

    print(f"\n  --- VAST.AI LAUNCH COMMANDS ---")
    print(f"  1. Upload: scp -P <PORT> -r {YAML_DIR} root@<HOST>:/workspace/v2_competition/")
    print(f"  2. Run:")
    print(f"     for f in /workspace/v2_competition/*.yaml; do")
    print(f'       boltz predict "$f" --out_dir /workspace/v2_competition/out \\')
    print(f"         --diffusion_samples 3 --recycling_steps 3")
    print(f"     done")
    print(f"  3. Download: scp -P <PORT> -r root@<HOST>:/workspace/v2_competition/out {OUTPUT_DIR}/")


def analyze():
    """Analyze Boltz-2 results for substrate competition."""
    print("=" * 60)
    print("SUBSTRATE COMPETITION — ANALYSIS")
    print("=" * 60)

    # Look for Boltz-2 output directories
    results_dir = os.path.join(OUTPUT_DIR, "out")
    if not os.path.exists(results_dir):
        print(f"\n  No results found at: {results_dir}")
        print("  Run Boltz-2 on Vast.ai first, then download results here.")
        return

    # Parse confidence JSON files from Boltz-2 output
    def find_confidence(pred_dir):
        """Find ipTM from Boltz-2 output directory."""
        # Boltz-2 outputs confidence_*.json or summary.json
        json_files = glob_mod.glob(os.path.join(pred_dir, "**", "confidence_*.json"), recursive=True)
        if not json_files:
            json_files = glob_mod.glob(os.path.join(pred_dir, "**", "*.json"), recursive=True)

        iptm_values = []
        for jf in json_files:
            try:
                with open(jf) as f:
                    data = json.load(f)
                if "iptm" in data:
                    iptm_values.append(data["iptm"])
                elif "ipTM" in data:
                    iptm_values.append(data["ipTM"])
            except (json.JSONDecodeError, KeyError):
                continue

        if iptm_values:
            return sum(iptm_values) / len(iptm_values)
        return None

    # Collect results
    controls = {}
    tests = {}
    sanity_iptm = None

    for entry in os.listdir(results_dir):
        pred_dir = os.path.join(results_dir, entry)
        if not os.path.isdir(pred_dir):
            continue

        iptm = find_confidence(pred_dir)
        if iptm is None:
            continue

        if entry.startswith("ctrl_"):
            sub = entry.replace("ctrl_", "")
            controls[sub] = iptm
        elif entry.startswith("test_"):
            sub = entry.replace("test_", "")
            tests[sub] = iptm
        elif entry.startswith("sanity_"):
            sanity_iptm = iptm

    if not controls and not tests:
        print("  Could not parse any results. Check output directory structure.")
        return

    # Report
    if sanity_iptm is not None:
        print(f"\n  Sanity check (ERAP2 + V2 only): ipTM = {sanity_iptm:.3f}")
        print(f"  (Expected: ~0.748 from prior runs)")
        if abs(sanity_iptm - 0.748) > 0.1:
            print("  WARNING: Sanity check deviates significantly from prior result!")

    print(f"\n  {'Substrate':<15} {'Control ipTM':>13} {'+ V2 ipTM':>12} {'Delta':>8} {'% Drop':>8} {'Verdict':>12}")
    print(f"  {'─'*15} {'─'*13} {'─'*12} {'─'*8} {'─'*8} {'─'*12}")

    analysis_results = {}
    for sub_name in SUBSTRATES:
        ctrl = controls.get(sub_name)
        test = tests.get(sub_name)

        if ctrl is not None and test is not None:
            delta = test - ctrl
            pct_drop = (ctrl - test) / ctrl * 100 if ctrl > 0 else 0

            if pct_drop > 50:
                verdict = "BLOCKED"
            elif pct_drop > 25:
                verdict = "IMPAIRED"
            elif pct_drop > 10:
                verdict = "MILD_EFFECT"
            else:
                verdict = "NO_EFFECT"

            print(f"  {sub_name:<15} {ctrl:>13.3f} {test:>12.3f} {delta:>+8.3f} {pct_drop:>7.1f}% {verdict:>12}")
            analysis_results[sub_name] = {
                "control_iptm": round(ctrl, 4),
                "test_iptm": round(test, 4),
                "delta": round(delta, 4),
                "pct_drop": round(pct_drop, 1),
                "verdict": verdict,
            }
        else:
            ctrl_str = f"{ctrl:.3f}" if ctrl else "MISSING"
            test_str = f"{test:.3f}" if test else "MISSING"
            print(f"  {sub_name:<15} {ctrl_str:>13} {test_str:>12} {'—':>8} {'—':>8} {'INCOMPLETE':>12}")

    # Overall verdict
    blocked_count = sum(1 for v in analysis_results.values() if v["verdict"] == "BLOCKED")
    impaired_count = sum(1 for v in analysis_results.values() if v["verdict"] in ("BLOCKED", "IMPAIRED"))
    total = len(analysis_results)

    print(f"\n  {'=' * 60}")
    if blocked_count == total and total > 0:
        print("  OVERALL: V2 FUNCTIONALLY BLOCKS substrate access")
        print("  All substrates show >50% ipTM drop — strong inhibitor candidate")
    elif impaired_count > total / 2:
        print("  OVERALL: V2 PARTIALLY IMPAIRS substrate access")
        print("  Majority of substrates affected — moderate inhibitor")
    elif impaired_count > 0:
        print("  OVERALL: V2 has MIXED effect on substrates")
        print("  Some substrates affected, others not — mechanism-dependent")
    else:
        print("  OVERALL: V2 does NOT block substrate access")
        print("  Substrates dock normally despite V2 presence — binder without inhibition")
    print(f"  {'=' * 60}")

    # Save
    out = {
        "analysis": "V2 substrate competition",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sanity_iptm": sanity_iptm,
        "substrates": analysis_results,
    }
    out_path = os.path.join(OUTPUT_DIR, "competition_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Results saved to: {out_path}")


if __name__ == "__main__":
    import sys as _sys
    if len(_sys.argv) > 1 and _sys.argv[1] == "analyze":
        analyze()
    else:
        generate()
