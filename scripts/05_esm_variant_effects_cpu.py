"""
05_esm_variant_effects_cpu.py — ESM-2 variant effect scoring (CPU-friendly).

Uses ESM-2 (650M params) to score the functional impact of ERAP2 mutations.
Runs on CPU — slower than GPU but works locally without Colab.

For each variant position, masks that residue and computes:
  - Wild-type log-probability (how "expected" the normal amino acid is)
  - Mutant log-probability (how "expected" the mutation is)
  - Delta (negative = mutation is disfavored by evolution = likely damaging)

Key variants scored:
  - rs2549794 → K528R (the Black Death → Crohn's link)
  - rs2248374 → splice variant (affects ERAP2 expression, not amino acid)
  - Additional known ERAP2 missense variants from ClinVar/literature

Outputs:
  data/processed/esm_variant_scores.csv
  data/processed/esm_conservation_profile.csv (sampled)

Requires: transformers, torch (CPU)
Estimated runtime: ~5-15 min on CPU for key variants, ~1-2 hrs for full scan
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import time
from pathlib import Path

# IMPORTANT: import torch BEFORE pandas on Windows —
# conda pandas loads MKL DLLs that conflict with pip-installed torch's shm.dll
import torch
from transformers import AutoTokenizer, EsmForMaskedLM
import pandas as pd

PROC_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
FASTA_PATH = PROC_DIR / "erap2_Q6P179.fasta"

# ──────────────────────────────────────────────────────────────
# ERAP2 variant map: rsID → amino acid change
# Mapped from genomic coordinates via UniProt Q6P179 annotations
# ──────────────────────────────────────────────────────────────
VARIANTS = {
    "rs2549794_proxy_K528R": {
        "position": 527,  # 0-indexed (K528 in 1-indexed)
        "wt_aa": "K",
        "mt_aa": "R",
        "context": "Black Death survival variant — proxy effect via linkage with ERAP2 haplotype",
        "note": "rs2549794 is intronic/regulatory; K528R is a linked coding variant on the same haplotype",
    },
    "rs26653_N392K": {
        "position": 391,  # 0-indexed
        "wt_aa": "N",
        "mt_aa": "K",
        "context": "ERAP2 missense variant — near active site zinc-binding residues",
    },
    "D198N_active_site": {
        "position": 197,
        "wt_aa": "D",
        "mt_aa": "N",
        "context": "Synthetic — near substrate binding region, tests sensitivity",
    },
    "H370A_zinc_binding": {
        "position": 369,  # Zinc binding residue
        "wt_aa": "H",
        "mt_aa": "A",
        "context": "Zinc-binding residue — mutation should be highly damaging",
    },
    "E371A_catalytic": {
        "position": 370,  # Catalytic residue
        "wt_aa": "E",
        "mt_aa": "A",
        "context": "Catalytic glutamate — mutation should abolish activity",
    },
    "H374A_zinc_binding": {
        "position": 373,  # Zinc binding residue
        "wt_aa": "H",
        "mt_aa": "A",
        "context": "Zinc-binding residue — mutation should be highly damaging",
    },
    "E393A_zinc_binding": {
        "position": 392,  # Zinc binding residue (HEXXH+E motif)
        "wt_aa": "E",
        "mt_aa": "A",
        "context": "Third zinc ligand (HEXXH+E motif) — mutation should be devastating",
    },
}


def load_sequence() -> str:
    """Load ERAP2 protein sequence from FASTA."""
    if not FASTA_PATH.exists():
        raise FileNotFoundError(
            f"ERAP2 FASTA not found at {FASTA_PATH}\n"
            "Run script 02 first to download from UniProt."
        )
    lines = FASTA_PATH.read_text(encoding="utf-8").strip().split("\n")
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    print(f"Loaded ERAP2 sequence: {len(seq)} amino acids")
    return seq


def load_model():
    """Load ESM-2 (650M) for CPU inference."""
    print("Loading ESM-2 model (650M params) — this takes ~1-2 min on CPU ...")
    model_name = "facebook/esm2_t33_650M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name)
    model.eval()
    param_count = sum(p.numel() for p in model.parameters()) / 1e6
    print(f"Model loaded: {param_count:.0f}M parameters (CPU mode)")
    return tokenizer, model


def score_single_mutation(
    sequence: str,
    position: int,
    wt_aa: str,
    mt_aa: str,
    tokenizer,
    model,
) -> dict:
    """
    Score a single amino acid mutation using ESM-2 masked marginals.

    Masks the target position, gets log-probabilities for all amino acids,
    and compares wild-type vs mutant likelihood.

    Returns dict with scores and interpretation.
    """
    # Verify the sequence has the expected wild-type AA
    actual_aa = sequence[position]
    if actual_aa != wt_aa:
        print(f"  WARNING: Expected {wt_aa} at position {position+1}, found {actual_aa}")

    # Create masked sequence
    masked_seq = sequence[:position] + tokenizer.mask_token + sequence[position + 1 :]

    # Tokenize — ESM-2 can handle full 960aa sequence on CPU (just slower)
    inputs = tokenizer(masked_seq, return_tensors="pt", truncation=True, max_length=1024)

    with torch.no_grad():
        outputs = model(**inputs)

    # Find the mask position in tokenized sequence
    mask_positions = (inputs.input_ids == tokenizer.mask_token_id).nonzero(as_tuple=True)
    if len(mask_positions[1]) == 0:
        return {"error": "Mask token not found in tokenized sequence"}

    mask_idx = mask_positions[1][0].item()
    logits = outputs.logits[0, mask_idx]
    log_probs = torch.log_softmax(logits, dim=-1)

    # Get token IDs for amino acids
    wt_token = tokenizer.convert_tokens_to_ids(wt_aa)
    mt_token = tokenizer.convert_tokens_to_ids(mt_aa)

    wt_logprob = log_probs[wt_token].item()
    mt_logprob = log_probs[mt_token].item()
    delta = mt_logprob - wt_logprob

    # Interpret the delta
    if delta < -3.0:
        effect = "HIGHLY DAMAGING"
    elif delta < -1.5:
        effect = "DAMAGING"
    elif delta < -0.5:
        effect = "MILDLY DELETERIOUS"
    elif delta < 0.5:
        effect = "NEUTRAL"
    else:
        effect = "POTENTIALLY BENEFICIAL"

    # Get top 5 most likely amino acids at this position (conservation context)
    top5_indices = torch.topk(log_probs, 5).indices.tolist()
    top5_aas = [tokenizer.convert_ids_to_tokens(idx) for idx in top5_indices]
    top5_probs = [log_probs[idx].item() for idx in top5_indices]

    return {
        "position_1indexed": position + 1,
        "wt_aa": wt_aa,
        "mt_aa": mt_aa,
        "mutation": f"{wt_aa}{position + 1}{mt_aa}",
        "wt_logprob": round(wt_logprob, 4),
        "mt_logprob": round(mt_logprob, 4),
        "delta_logprob": round(delta, 4),
        "effect": effect,
        "top5_aa": top5_aas,
        "top5_logprob": [round(p, 4) for p in top5_probs],
        "wt_rank": (log_probs > wt_logprob).sum().item() + 1,  # Rank of WT among all AAs
    }


def score_all_variants(sequence: str, tokenizer, model) -> pd.DataFrame:
    """Score all defined ERAP2 variants."""
    print(f"\nScoring {len(VARIANTS)} variants ...")
    print("-" * 80)

    results = []
    for variant_id, info in VARIANTS.items():
        pos = info["position"]
        wt = info["wt_aa"]
        mt = info["mt_aa"]
        ctx = info["context"]

        print(f"\n  {variant_id}: {wt}{pos+1}{mt}")
        print(f"    Context: {ctx}")

        start = time.time()
        result = score_single_mutation(sequence, pos, wt, mt, tokenizer, model)
        elapsed = time.time() - start

        result["variant_id"] = variant_id
        result["context"] = ctx
        result["compute_time_sec"] = round(elapsed, 1)

        print(f"    WT log-prob: {result['wt_logprob']:.4f} (rank {result.get('wt_rank', '?')}/20)")
        print(f"    MT log-prob: {result['mt_logprob']:.4f}")
        print(f"    Delta:       {result['delta_logprob']:.4f}")
        print(f"    Effect:      {result['effect']}")
        print(f"    Top 5 AAs:   {result['top5_aa']}")
        print(f"    Time:        {elapsed:.1f}s")

        results.append(result)

    df = pd.DataFrame(results)
    return df


def conservation_scan(sequence: str, tokenizer, model, step: int = 20) -> pd.DataFrame:
    """
    Quick conservation scan — sample every Nth position to find
    highly conserved regions (functionally important).
    """
    n_positions = len(range(0, len(sequence), step))
    print(f"\nConservation scan: sampling every {step}th position ({n_positions} positions) ...")
    print("  (Full scan of all 960 positions would take ~30-60 min on CPU)")

    rows = []
    for i, pos in enumerate(range(0, len(sequence), step)):
        wt_aa = sequence[pos]

        masked_seq = sequence[:pos] + tokenizer.mask_token + sequence[pos + 1 :]
        inputs = tokenizer(masked_seq, return_tensors="pt", truncation=True, max_length=1024)

        with torch.no_grad():
            outputs = model(**inputs)

        mask_positions = (inputs.input_ids == tokenizer.mask_token_id).nonzero(as_tuple=True)
        if len(mask_positions[1]) == 0:
            continue

        mask_idx = mask_positions[1][0].item()
        logits = outputs.logits[0, mask_idx]
        log_probs = torch.log_softmax(logits, dim=-1)

        wt_token = tokenizer.convert_tokens_to_ids(wt_aa)
        wt_logprob = log_probs[wt_token].item()

        # Entropy as measure of position variability
        probs = torch.softmax(logits, dim=-1)
        entropy = -(probs * torch.log(probs + 1e-10)).sum().item()

        rows.append({
            "position": pos + 1,
            "aa": wt_aa,
            "wt_logprob": round(wt_logprob, 4),
            "entropy": round(entropy, 4),
            "conserved": wt_logprob > -0.5,  # High WT probability = conserved
        })

        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{n_positions} positions scanned ...")

    df = pd.DataFrame(rows)
    conserved_count = df["conserved"].sum()
    print(f"  Highly conserved positions: {conserved_count}/{len(df)} ({conserved_count/len(df)*100:.0f}%)")
    return df


def main():
    print("=" * 70)
    print("ESM-2 Variant Effect Scoring — ERAP2 (CPU Mode)")
    print("=" * 70)

    # Load sequence
    sequence = load_sequence()

    # Validate expected residues at variant positions
    print("\nValidating variant positions against sequence ...")
    for vid, info in VARIANTS.items():
        actual = sequence[info["position"]]
        expected = info["wt_aa"]
        status = "OK" if actual == expected else f"MISMATCH (found {actual})"
        print(f"  {vid}: position {info['position']+1} = {actual} (expected {expected}) — {status}")

    # Load model
    tokenizer, model = load_model()

    # Score all variants
    variant_df = score_all_variants(sequence, tokenizer, model)

    # Save variant scores
    out_variants = PROC_DIR / "esm_variant_scores.csv"
    # Drop list columns for CSV (save separately)
    csv_df = variant_df.drop(columns=["top5_aa", "top5_logprob"], errors="ignore")
    csv_df.to_csv(out_variants, index=False)
    print(f"\nSaved variant scores → {out_variants.name}")

    # Save full results as JSON (preserves lists)
    out_json = PROC_DIR / "esm_variant_scores.json"
    out_json.write_text(
        json.dumps(variant_df.to_dict(orient="records"), indent=2),
        encoding="utf-8",
    )

    # Conservation scan (sampled)
    conservation_df = conservation_scan(sequence, tokenizer, model, step=20)
    out_conservation = PROC_DIR / "esm_conservation_profile.csv"
    conservation_df.to_csv(out_conservation, index=False)
    print(f"Saved conservation profile → {out_conservation.name}")

    # Final summary
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"\n{'Mutation':<20} {'Delta':>8} {'Effect':<25} {'Context'}")
    print("-" * 90)
    for _, row in variant_df.iterrows():
        print(f"  {row['mutation']:<18} {row['delta_logprob']:>8.4f} {row['effect']:<25} {row['context'][:40]}")

    print(f"\nKey finding for PoC:")
    poc_row = variant_df[variant_df["variant_id"] == "rs2549794_proxy_K528R"]
    if not poc_row.empty:
        r = poc_row.iloc[0]
        print(f"  K528R (Black Death variant): delta = {r['delta_logprob']:.4f} → {r['effect']}")
        if r["delta_logprob"] < -0.5:
            print(f"  This confirms the mutation is evolutionarily disfavored,")
            print(f"  consistent with it being protective during plague but")
            print(f"  potentially deleterious for normal aminopeptidase function.")
        else:
            print(f"  The mutation appears relatively neutral by ESM-2 scoring.")
            print(f"  This is expected for a variant under balancing selection —")
            print(f"  it may alter specificity rather than fold stability.")

    print(f"\n  Validation check: zinc-binding mutations (H370A, H374A, E393A)")
    print(f"  should score as HIGHLY DAMAGING or DAMAGING.")
    zinc_rows = variant_df[variant_df["variant_id"].str.contains("zinc|catalytic")]
    if not zinc_rows.empty:
        all_damaging = all(
            "DAMAGING" in row["effect"] or "HIGHLY" in row["effect"]
            for _, row in zinc_rows.iterrows()
        )
        print(f"  Result: {'PASS — all scored as damaging' if all_damaging else 'CHECK — review scores above'}")

    print("\nDone. Next: Upload variant scores + FASTA to Colab for notebook 05 (structure prediction).")


if __name__ == "__main__":
    main()
