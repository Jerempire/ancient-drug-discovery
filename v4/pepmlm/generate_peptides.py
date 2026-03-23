"""PepMLM-based peptide generation for ERAP2 K392/N392 selectivity.

Generates candidate peptide binders conditioned on the ERAP2 substrate channel
(residues 388-410) for both K392 (wildtype) and N392 (plague variant).
Cross-filters for K392-selective hits, then ranks by pseudo-perplexity.

Usage (Vast.ai GPU or local CPU):
    python generate_peptides.py [--num_peptides 200] [--top_k 3] [--device cuda]

Outputs:
    results/pepmlm_k392_candidates.csv
    results/pepmlm_n392_candidates.csv
    results/pepmlm_selective_hits.csv   (K392-selective after cross-filter)
    results/pepmlm_all_for_boltz2.fasta (top candidates formatted for Boltz-2)
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import os
import json
import time
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from torch.distributions import Categorical
from transformers import AutoTokenizer, AutoModelForMaskedLM

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── ERAP2 sequences ──────────────────────────────────────────────────────────
# Full ERAP2 sequence (Q6P179), residues 1-960
ERAP2_FULL = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPS"
    "VVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYP"
    "AHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLF"
    "KANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASP"
    "DKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLW"
    "VTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPI"
    "SKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTS"
    "GGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQER"
    "YLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYF"
    "KPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAG"
    "WNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRE"
    "NWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLR"
    "TWLMVNT"
)

# Substrate channel region: residues 388-410 (0-indexed: 387-409)
# K392 = position 391 (0-indexed) within full sequence
CHANNEL_START = 387  # 0-indexed
CHANNEL_END = 410    # exclusive

CHANNEL_K392 = ERAP2_FULL[CHANNEL_START:CHANNEL_END]  # wildtype K at position 392
CHANNEL_N392 = CHANNEL_K392[:4] + "N" + CHANNEL_K392[5:]  # K→N at position 392

# Also try a wider context window (residues 350-430) for better conditioning
WIDE_START = 349
WIDE_END = 430
WIDE_K392 = ERAP2_FULL[WIDE_START:WIDE_END]
WIDE_N392 = WIDE_K392[: (391 - WIDE_START)] + "N" + WIDE_K392[(392 - WIDE_START) :]


def compute_pseudo_perplexity(model, tokenizer, protein_seq, binder_seq):
    """Compute pseudo-perplexity of a peptide given a target sequence.

    Lower PPL = model is more confident this peptide belongs with the target.
    """
    sequence = protein_seq + binder_seq
    original_input = tokenizer.encode(sequence, return_tensors="pt").to(model.device)
    length_of_binder = len(binder_seq)

    masked_inputs = original_input.repeat(length_of_binder, 1)
    positions_to_mask = torch.arange(-length_of_binder - 1, -1, device=model.device)
    masked_inputs[torch.arange(length_of_binder), positions_to_mask] = tokenizer.mask_token_id

    labels = torch.full_like(masked_inputs, -100)
    labels[torch.arange(length_of_binder), positions_to_mask] = original_input[0, positions_to_mask]

    with torch.no_grad():
        outputs = model(masked_inputs, labels=labels)

    return np.exp(outputs.loss.item())


def generate_peptides(model, tokenizer, target_seq, peptide_length=9, top_k=3, num_binders=200):
    """Generate peptide binders conditioned on target sequence."""
    results = []
    for i in range(num_binders):
        masked_peptide = "<mask>" * peptide_length
        input_sequence = target_seq + masked_peptide
        inputs = tokenizer(input_sequence, return_tensors="pt").to(model.device)

        with torch.no_grad():
            logits = model(**inputs).logits

        mask_indices = (inputs["input_ids"] == tokenizer.mask_token_id).nonzero(as_tuple=True)[1]
        logits_at_masks = logits[0, mask_indices]

        top_k_logits, top_k_indices = logits_at_masks.topk(top_k, dim=-1)
        probabilities = torch.nn.functional.softmax(top_k_logits, dim=-1)
        predicted_indices = Categorical(probabilities).sample()
        predicted_token_ids = top_k_indices.gather(-1, predicted_indices.unsqueeze(-1)).squeeze(-1)

        binder = tokenizer.decode(predicted_token_ids, skip_special_tokens=True).replace(" ", "")
        results.append(binder)

        if (i + 1) % 50 == 0:
            print(f"  Generated {i + 1}/{num_binders} peptides")

    return results


def score_peptides(model, tokenizer, target_seq, peptides):
    """Score a list of peptides against a target sequence."""
    scored = []
    for i, pep in enumerate(peptides):
        ppl = compute_pseudo_perplexity(model, tokenizer, target_seq, pep)
        scored.append({"peptide": pep, "length": len(pep), "ppl": ppl})
        if (i + 1) % 50 == 0:
            print(f"  Scored {i + 1}/{len(peptides)}")
    return scored


def cross_filter_selective(k392_df, n392_df, ppl_ratio_threshold=1.3):
    """Find peptides that bind K392 better than N392.

    A peptide is K392-selective if its PPL on K392 is lower than on N392,
    meaning the model is more confident it belongs with K392.

    ppl_ratio = ppl_n392 / ppl_k392. Higher ratio = more K392-selective.
    """
    merged = k392_df.merge(n392_df, on="peptide", suffixes=("_k392", "_n392"))
    merged["ppl_ratio"] = merged["ppl_n392"] / merged["ppl_k392"]
    merged["selectivity"] = merged["ppl_ratio"].apply(
        lambda r: "K392-SEL" if r > ppl_ratio_threshold else ("N392-SEL" if r < 1 / ppl_ratio_threshold else "NEUTRAL")
    )
    return merged.sort_values("ppl_ratio", ascending=False)


def write_boltz2_fasta(selective_df, output_path, top_n=30):
    """Write top selective peptides as FASTA for Boltz-2 docking."""
    top = selective_df.head(top_n)
    with open(output_path, "w") as f:
        for _, row in top.iterrows():
            tag = row["selectivity"]
            ratio = row["ppl_ratio"]
            f.write(f">pepmlm_{row.name:03d}|{tag}|ratio={ratio:.2f}\n")
            f.write(f"{row['peptide']}\n")
    print(f"Wrote {len(top)} peptides to {output_path}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="PepMLM peptide generation for ERAP2")
    parser.add_argument("--num_peptides", type=int, default=200, help="Peptides per target variant")
    parser.add_argument("--top_k", type=int, default=3, help="Top-k sampling parameter")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--peptide_lengths", type=str, default="10,11,12,13",
                        help="Comma-separated peptide lengths to generate. "
                        "V4 DiffPepDock showed 9-mers universally prefer N392 — "
                        "only 10+ mers reach K392 salt bridge.")
    parser.add_argument("--use_wide_context", action="store_true",
                        help="Use wider channel context (res 350-430) instead of minimal (388-410)")
    parser.add_argument("--ppl_ratio_threshold", type=float, default=1.3,
                        help="PPL ratio threshold for selectivity (default 1.3)")
    parser.add_argument("--top_n_boltz2", type=int, default=30,
                        help="Number of top candidates to write for Boltz-2")
    args = parser.parse_args()

    lengths = [int(x) for x in args.peptide_lengths.split(",")]
    per_length = args.num_peptides // len(lengths)

    # Select channel context
    if args.use_wide_context:
        ctx_k392, ctx_n392 = WIDE_K392, WIDE_N392
        print(f"Using WIDE context: residues 350-430 ({len(ctx_k392)} aa)")
    else:
        ctx_k392, ctx_n392 = CHANNEL_K392, CHANNEL_N392
        print(f"Using CHANNEL context: residues 388-410 ({len(ctx_k392)} aa)")

    print(f"K392 channel: {ctx_k392}")
    print(f"N392 channel: {ctx_n392}")
    print(f"Difference at pos {391 - (WIDE_START if args.use_wide_context else CHANNEL_START)}: "
          f"K -> N")
    print(f"Device: {args.device}")
    print(f"Generating {args.num_peptides} peptides per variant "
          f"({per_length} per length: {lengths})")
    print()

    # Load model
    print("Loading PepMLM-650M...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("TianlaiChen/PepMLM-650M")
    model = AutoModelForMaskedLM.from_pretrained("TianlaiChen/PepMLM-650M").to(args.device)
    model.eval()
    print(f"Model loaded in {time.time() - t0:.1f}s")
    print()

    # ── Generate for K392 ─────────────────────────────────────────────────────
    print("=== Generating peptides for K392 (wildtype) ===")
    k392_peptides = []
    for plen in lengths:
        print(f"  Length {plen}:")
        peps = generate_peptides(model, tokenizer, ctx_k392,
                                 peptide_length=plen, top_k=args.top_k,
                                 num_binders=per_length)
        k392_peptides.extend(peps)

    # Deduplicate
    k392_unique = list(set(k392_peptides))
    print(f"K392: {len(k392_peptides)} generated, {len(k392_unique)} unique")
    print()

    # ── Generate for N392 ─────────────────────────────────────────────────────
    print("=== Generating peptides for N392 (plague variant) ===")
    n392_peptides = []
    for plen in lengths:
        print(f"  Length {plen}:")
        peps = generate_peptides(model, tokenizer, ctx_n392,
                                 peptide_length=plen, top_k=args.top_k,
                                 num_binders=per_length)
        n392_peptides.extend(peps)

    n392_unique = list(set(n392_peptides))
    print(f"N392: {len(n392_peptides)} generated, {len(n392_unique)} unique")
    print()

    # Combine all unique peptides from both sets
    all_unique = list(set(k392_unique + n392_unique))
    print(f"Total unique peptides across both: {len(all_unique)}")
    print()

    # ── Score ALL peptides against BOTH targets ───────────────────────────────
    print("=== Scoring all peptides against K392 ===")
    k392_scored = score_peptides(model, tokenizer, ctx_k392, all_unique)
    k392_df = pd.DataFrame(k392_scored)
    k392_df.to_csv(RESULTS_DIR / "pepmlm_k392_candidates.csv", index=False)
    print(f"Saved {len(k392_df)} K392 scores")
    print()

    print("=== Scoring all peptides against N392 ===")
    n392_scored = score_peptides(model, tokenizer, ctx_n392, all_unique)
    n392_df = pd.DataFrame(n392_scored)
    n392_df.to_csv(RESULTS_DIR / "pepmlm_n392_candidates.csv", index=False)
    print(f"Saved {len(n392_df)} N392 scores")
    print()

    # ── Cross-filter for selectivity ──────────────────────────────────────────
    print("=== Cross-filtering for K392 selectivity ===")
    selective = cross_filter_selective(k392_df, n392_df, args.ppl_ratio_threshold)
    selective.to_csv(RESULTS_DIR / "pepmlm_selective_hits.csv", index=False)

    n_k392_sel = (selective["selectivity"] == "K392-SEL").sum()
    n_n392_sel = (selective["selectivity"] == "N392-SEL").sum()
    n_neutral = (selective["selectivity"] == "NEUTRAL").sum()
    print(f"Results: {n_k392_sel} K392-selective, {n_n392_sel} N392-selective, {n_neutral} neutral")
    print()

    # Show top 10 K392-selective
    top_k392 = selective[selective["selectivity"] == "K392-SEL"].head(10)
    if len(top_k392) > 0:
        print("Top 10 K392-selective peptides:")
        print(f"{'Peptide':>15s}  {'PPL_K392':>10s}  {'PPL_N392':>10s}  {'Ratio':>8s}")
        print("-" * 50)
        for _, row in top_k392.iterrows():
            print(f"{row['peptide']:>15s}  {row['ppl_k392']:>10.2f}  {row['ppl_n392']:>10.2f}  {row['ppl_ratio']:>8.2f}")
    print()

    # Check if salt bridge model is supported: are E/D P1 residues enriched?
    k392_sel_peptides = selective[selective["selectivity"] == "K392-SEL"]["peptide"].tolist()
    if k392_sel_peptides:
        p1_counts = {}
        for pep in k392_sel_peptides:
            p1 = pep[0]
            p1_counts[p1] = p1_counts.get(p1, 0) + 1
        print("P1 residue distribution in K392-selective hits:")
        for aa, count in sorted(p1_counts.items(), key=lambda x: -x[1]):
            pct = 100 * count / len(k392_sel_peptides)
            marker = " <-- salt bridge" if aa in ("E", "D") else ""
            print(f"  {aa}: {count} ({pct:.0f}%){marker}")
    print()

    # ── Write FASTA for Boltz-2 ───────────────────────────────────────────────
    k392_selective = selective[selective["selectivity"] == "K392-SEL"]
    if len(k392_selective) > 0:
        write_boltz2_fasta(k392_selective, RESULTS_DIR / "pepmlm_all_for_boltz2.fasta",
                           top_n=args.top_n_boltz2)
    else:
        print("WARNING: No K392-selective hits found. Writing top by K392 PPL instead.")
        fallback = selective.sort_values("ppl_k392").head(args.top_n_boltz2)
        write_boltz2_fasta(fallback, RESULTS_DIR / "pepmlm_all_for_boltz2.fasta",
                           top_n=args.top_n_boltz2)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("=" * 60)
    print("PEPMLM GENERATION COMPLETE")
    print(f"  Total unique peptides: {len(all_unique)}")
    print(f"  K392-selective:        {n_k392_sel}")
    print(f"  N392-selective:        {n_n392_sel}")
    print(f"  Neutral:               {n_neutral}")
    print(f"  Top candidates FASTA:  {RESULTS_DIR / 'pepmlm_all_for_boltz2.fasta'}")
    print(f"  Next: run Boltz-2 on pepmlm_all_for_boltz2.fasta")
    print("=" * 60)


if __name__ == "__main__":
    main()
