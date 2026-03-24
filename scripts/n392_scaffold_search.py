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
    """Analyze Boltz-2 results to find N392-leaning scaffolds."""
    import glob

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

    # Classify scaffolds
    n392_leaning = []
    neutral = []
    k392_leaning = []

    for idx, variants in sorted(scores.items()):
        k = variants.get("K392", 0)
        n = variants.get("N392", 0)
        delta = k - n
        scaffold = data["ala_baselines"][idx]

        entry = {
            "idx": idx,
            "scaffold": scaffold["original"],
            "ala_version": scaffold["ala_version"],
            "length": scaffold["length"],
            "K392": k,
            "N392": n,
            "delta": delta,
        }

        if delta < -0.02:
            n392_leaning.append(entry)
        elif delta > 0.02:
            k392_leaning.append(entry)
        else:
            neutral.append(entry)

    # Report
    print("=" * 70)
    print("N392 SCAFFOLD DISCOVERY RESULTS")
    print("=" * 70)
    print(f"Total scaffolds screened: {len(scores)}")
    print(f"N392-leaning: {len(n392_leaning)}")
    print(f"Neutral: {len(neutral)}")
    print(f"K392-leaning: {len(k392_leaning)}")

    print(f"\n{'='*70}")
    print("N392-LEANING SCAFFOLDS (best candidates for N392 arm)")
    print("=" * 70)
    for s in sorted(n392_leaning, key=lambda x: x["delta"]):
        print(f"  {s['scaffold']:<20} ({s['length']}aa)  K392={s['K392']:.3f}  N392={s['N392']:.3f}  delta={s['delta']:+.3f}")

    print(f"\n{'='*70}")
    print("NEUTRAL SCAFFOLDS (potential for either arm)")
    print("=" * 70)
    for s in sorted(neutral, key=lambda x: x["delta"]):
        print(f"  {s['scaffold']:<20} ({s['length']}aa)  K392={s['K392']:.3f}  N392={s['N392']:.3f}  delta={s['delta']:+.3f}")

    # Save
    result = {
        "n392_leaning": n392_leaning,
        "neutral": neutral,
        "k392_leaning": k392_leaning,
        "total_screened": len(scores),
    }
    out_path = OUT_DIR / "scaffold_bias_results.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")

    if n392_leaning:
        print(f"\nNext: Run P1 scan (A, L, V, D, F, I) on top {min(5, len(n392_leaning))} N392-leaning scaffolds")


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
