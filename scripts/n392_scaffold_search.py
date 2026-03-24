"""
N392 scaffold discovery pipeline.

Strategy:
1. Generate ~100 PepMLM scaffolds conditioned on N392 channel (more diversity)
2. Put Ala at P1 on each (neutral baseline test)
3. Screen against both K392 and N392 channels via Boltz-2
4. Identify scaffolds that are N392-leaning or neutral at baseline
5. Output: Boltz-2 YAMLs ready for Vast.ai

The VKLLLL scaffold has an inherent K392 bias — we need different scaffolds
for the N392 arm. This script finds them.

Usage:
  # Step 1: Generate scaffolds + YAMLs (local, CPU)
  python scripts/n392_scaffold_search.py generate

  # Step 2: After Boltz-2 screening on Vast.ai, analyze results
  python scripts/n392_scaffold_search.py analyze /path/to/results/
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import os
from pathlib import Path
from datetime import datetime
from collections import defaultdict

PROJECT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT / "data" / "results" / "n392_scaffold_search"

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
K392_CHANNEL = ERAP2_FULL[349:450]
N392_CHANNEL = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:450]

# PepMLM config — more samples, focused on N392 channel
NUM_SCAFFOLDS = 120  # generate more for diversity
PEPTIDE_LENGTHS = [10, 11, 12, 13, 15]  # skip 8-9 (too short)
TOP_K = 5       # more diversity in sampling
TEMPERATURE = 1.5  # higher temp = more diverse


def cmd_generate():
    """Generate N392-conditioned scaffolds and create Ala-baseline YAMLs."""
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    yaml_dir = OUT_DIR / "boltz_yamls"
    yaml_dir.mkdir(exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    model_id = "TianlaiChen/PepMLM-650M"
    print(f"Loading {model_id}...")
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForMaskedLM.from_pretrained(model_id).to(device)
    model.eval()

    # Generate scaffolds conditioned on N392 channel
    print(f"\nGenerating {NUM_SCAFFOLDS} scaffolds conditioned on N392 channel...")
    print(f"Lengths: {PEPTIDE_LENGTHS}, temp={TEMPERATURE}, top_k={TOP_K}\n")

    all_scaffolds = []
    seen = set()
    per_length = NUM_SCAFFOLDS // len(PEPTIDE_LENGTHS)

    for pep_len in PEPTIDE_LENGTHS:
        attempts = per_length * 5  # oversample for uniqueness
        count = 0
        print(f"  {pep_len}-mers: generating...", end=" ", flush=True)

        for _ in range(attempts):
            if count >= per_length:
                break

            mask_token = tokenizer.mask_token
            masked_seq = N392_CHANNEL + mask_token * pep_len
            inputs = tokenizer(masked_seq, return_tensors="pt").to(device)

            with torch.no_grad():
                outputs = model(**inputs)
            logits = outputs.logits[0]

            mask_positions = (inputs["input_ids"][0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0]
            peptide_tokens = []
            for pos in mask_positions:
                pos_logits = logits[pos] / TEMPERATURE
                probs = torch.softmax(pos_logits, dim=-1)
                top_k_probs, top_k_ids = torch.topk(probs, TOP_K)
                top_k_probs = top_k_probs / top_k_probs.sum()
                chosen = top_k_ids[torch.multinomial(top_k_probs, 1)]
                peptide_tokens.append(chosen.item())

            peptide = tokenizer.decode(peptide_tokens).replace(" ", "")
            if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in peptide) and len(peptide) == pep_len and peptide not in seen:
                seen.add(peptide)
                all_scaffolds.append({"sequence": peptide, "length": pep_len, "p1": peptide[0]})
                count += 1

        print(f"got {count}")

    print(f"\nTotal unique scaffolds: {len(all_scaffolds)}")

    # Create Ala-P1 baseline versions (swap P1 to A)
    ala_scaffolds = []
    for s in all_scaffolds:
        ala_seq = "A" + s["sequence"][1:]
        ala_scaffolds.append({
            "original": s["sequence"],
            "ala_version": ala_seq,
            "length": s["length"],
            "original_p1": s["p1"],
        })

    # Write Boltz-2 YAMLs: each Ala-scaffold x both variants
    yaml_count = 0
    for i, s in enumerate(ala_scaffolds):
        for variant, channel in [("K392", K392_CHANNEL), ("N392", N392_CHANNEL)]:
            name = f"n392s_{i:03d}_{s['length']}aa_{variant}"
            yaml_path = yaml_dir / f"{name}.yaml"
            yaml_path.write_text(
                f"version: 1\nsequences:\n"
                f"  - protein:\n      id: A\n      sequence: {channel}\n      msa: empty\n"
                f"  - protein:\n      id: B\n      sequence: {s['ala_version']}\n      msa: empty\n"
            )
            yaml_count += 1

    # Save scaffold data
    output = {
        "metadata": {
            "generated_at": datetime.now().isoformat(),
            "model": "TianlaiChen/PepMLM-650M",
            "conditioned_on": "N392_channel",
            "temperature": TEMPERATURE,
            "top_k": TOP_K,
            "num_scaffolds": len(all_scaffolds),
            "lengths": PEPTIDE_LENGTHS,
        },
        "scaffolds": all_scaffolds,
        "ala_baselines": ala_scaffolds,
    }
    json_path = OUT_DIR / "n392_scaffolds.json"
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nScaffolds saved to {json_path}")
    print(f"Boltz-2 YAMLs: {yaml_count} files in {yaml_dir}")
    print(f"\nNext: upload {yaml_dir} to Vast.ai and run Boltz-2")
    print(f"Then: python scripts/n392_scaffold_search.py analyze /path/to/results/")


def cmd_analyze(results_dir: str):
    """Analyze Boltz-2 results — rank by K392 potency first, selectivity second."""
    from statistics import mean, median

    results_path = Path(results_dir)
    scaffolds_json = OUT_DIR / "n392_scaffolds.json"

    with open(scaffolds_json) as f:
        data = json.load(f)

    # Parse Boltz-2 confidence files
    scores = {}  # {scaffold_idx: {variant: iptm}}
    for d in sorted(results_path.iterdir()):
        if not d.is_dir() or not d.name.startswith("n392s_"):
            continue
        parts = d.name.split("_")
        idx = int(parts[1])
        variant = parts[-1]  # K392 or N392

        conf_files = list(d.rglob("confidence_*.json"))
        iptms = []
        for cf in conf_files:
            with open(cf) as f:
                c = json.load(f)
            iptm = c.get("iptm", c.get("i_ptm", 0))
            if iptm:
                iptms.append(iptm)

        if iptms:
            if idx not in scores:
                scores[idx] = {}
            scores[idx][variant] = sum(iptms) / len(iptms)

    # Build full scaffold list with scores
    all_entries = []
    for idx, variants in sorted(scores.items()):
        k = variants.get("K392", 0)
        n = variants.get("N392", 0)
        delta = k - n
        scaffold = data["ala_baselines"][idx]

        # Bucket classification (per ChatGPT framework)
        if k >= 0.7 and delta > 0.05:
            bucket = "B_allele_selective"
        elif k >= 0.7:
            bucket = "A_erap2_inhibitor"
        elif k >= 0.5:
            bucket = "C_moderate"
        else:
            bucket = "D_weak"

        all_entries.append({
            "idx": idx,
            "scaffold": scaffold["original"],
            "ala_version": scaffold["ala_version"],
            "length": scaffold["length"],
            "K392": k,
            "N392": n,
            "delta": delta,
            "bucket": bucket,
        })

    # --- REPORT ---
    SEP = "=" * 80
    DASH = "-" * 80

    print(SEP)
    print("N392 SCAFFOLD SCREEN — RANKED BY K392 POTENCY")
    print("Primary: K392 binding strength | Secondary: K392 vs N392 selectivity")
    print(SEP)

    # Summary stats
    k392_scores = [e["K392"] for e in all_entries]
    n392_scores = [e["N392"] for e in all_entries]
    deltas = [e["delta"] for e in all_entries]

    print(f"\nTotal scaffolds screened: {len(all_entries)}")
    print(f"K392 ipTM: mean={mean(k392_scores):.3f}, median={median(k392_scores):.3f}, max={max(k392_scores):.3f}")
    print(f"N392 ipTM: mean={mean(n392_scores):.3f}, median={median(n392_scores):.3f}, max={max(n392_scores):.3f}")
    print(f"Delta (K-N): mean={mean(deltas):+.3f}, median={median(deltas):+.3f}")
    print(f"K392-leaning: {sum(1 for d in deltas if d > 0.02)}")
    print(f"Neutral: {sum(1 for d in deltas if -0.02 <= d <= 0.02)}")
    print(f"N392-leaning: {sum(1 for d in deltas if d < -0.02)}")

    # Bucket counts
    buckets = defaultdict(list)
    for e in all_entries:
        buckets[e["bucket"]].append(e)

    print(f"\n{'='*80}")
    print("BUCKET CLASSIFICATION")
    print(f"{'='*80}")
    print(f"Bucket B (K392 >= 0.7 AND selective > +0.05): {len(buckets['B_allele_selective'])} scaffolds")
    print(f"Bucket A (K392 >= 0.7, any selectivity):      {len(buckets['A_erap2_inhibitor'])} scaffolds")
    print(f"Bucket C (K392 0.5-0.7):                      {len(buckets['C_moderate'])} scaffolds")
    print(f"Bucket D (K392 < 0.5):                        {len(buckets['D_weak'])} scaffolds")

    # Top 15 by K392 potency
    by_k392 = sorted(all_entries, key=lambda x: -x["K392"])
    print(f"\n{'='*80}")
    print("TOP 15 BY K392 BINDING STRENGTH")
    print(f"{'='*80}")
    print("%-5s %-22s %4s %8s %8s %8s %s" % ("Rank", "Scaffold", "Len", "K392", "N392", "Delta", "Bucket"))
    print(DASH)
    for i, e in enumerate(by_k392[:15], 1):
        print("%-5d %-22s %4d %8.3f %8.3f %+8.3f %s" % (
            i, e["scaffold"], e["length"], e["K392"], e["N392"], e["delta"], e["bucket"]))

    # Top 10 Bucket B (allele-selective)
    bucket_b = sorted(buckets["B_allele_selective"], key=lambda x: -x["K392"])
    if bucket_b:
        print(f"\n{'='*80}")
        print("BUCKET B: ALLELE-SELECTIVE K392 INHIBITORS (best outcome)")
        print("K392 >= 0.7 AND delta > +0.05")
        print(f"{'='*80}")
        print("%-5s %-22s %4s %8s %8s %8s" % ("Rank", "Scaffold", "Len", "K392", "N392", "Delta"))
        print(DASH)
        for i, e in enumerate(bucket_b[:10], 1):
            print("%-5d %-22s %4d %8.3f %8.3f %+8.3f" % (
                i, e["scaffold"], e["length"], e["K392"], e["N392"], e["delta"]))

    # N392-leaning (for mechanistic interest)
    n392_leaning = sorted([e for e in all_entries if e["delta"] < -0.02], key=lambda x: x["delta"])
    if n392_leaning:
        print(f"\n{'='*80}")
        print("N392-LEANING SCAFFOLDS (counterscreen interest / mechanistic)")
        print(f"{'='*80}")
        for e in n392_leaning[:10]:
            print("  %-22s (%daa)  K392=%8.3f  N392=%8.3f  delta=%+.3f" % (
                e["scaffold"], e["length"], e["K392"], e["N392"], e["delta"]))

    # Comparison to VKLLLL
    print(f"\n{'='*80}")
    print("COMPARISON TO EXISTING LEADS")
    print(f"{'='*80}")
    print("VKLLLL (P1=V): K392=0.801, N392=0.665, delta=+0.137 (bias-corrected +0.125)")
    print("VKLLLL (P1=E): K392=0.796, N392=0.709, delta=+0.087 (bias-corrected +0.076)")
    if bucket_b:
        best = bucket_b[0]
        print(f"Best new Bucket B: {best['scaffold']} K392={best['K392']:.3f} delta={best['delta']:+.3f}")
        if best["K392"] > 0.801:
            print(">>> NEW SCAFFOLD BEATS VKLLLL ON K392 POTENCY")
        elif best["delta"] > 0.137:
            print(">>> NEW SCAFFOLD BEATS VKLLLL ON SELECTIVITY")
        else:
            print(">>> VKLLLL remains the lead — new scaffolds are backup/diversity")

    # P1 scan recommendations
    print(f"\n{'='*80}")
    print("RECOMMENDED NEXT STEPS")
    print(f"{'='*80}")
    top_for_p1 = by_k392[:5]
    print(f"P1 scan candidates (top 5 by K392, test V/E/L/A/D at P1):")
    for e in top_for_p1:
        print(f"  {e['scaffold']} ({e['length']}aa, K392={e['K392']:.3f}, delta={e['delta']:+.3f})")

    # Save
    result = {
        "summary": {
            "total_screened": len(all_entries),
            "k392_mean": mean(k392_scores),
            "n392_mean": mean(n392_scores),
            "delta_mean": mean(deltas),
            "bucket_b_count": len(buckets["B_allele_selective"]),
            "bucket_a_count": len(buckets["A_erap2_inhibitor"]),
            "n392_leaning_count": len(n392_leaning),
        },
        "all_scaffolds": sorted(all_entries, key=lambda x: -x["K392"]),
        "bucket_b": bucket_b,
        "n392_leaning": n392_leaning,
        "p1_scan_candidates": [e["scaffold"] for e in top_for_p1],
    }
    out_path = OUT_DIR / "scaffold_screen_results.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python n392_scaffold_search.py [generate|analyze <results_dir>]")
        sys.exit(1)

    cmd = sys.argv[1]
    if cmd == "generate":
        cmd_generate()
    elif cmd == "analyze" and len(sys.argv) > 2:
        cmd_analyze(sys.argv[2])
    else:
        print(f"Unknown command: {cmd}")
        sys.exit(1)
