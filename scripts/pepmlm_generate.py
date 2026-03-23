"""
PepMLM: Generate peptide binder candidates for ERAP2 K392 channel region.

Generates linear peptide binders conditioned on the ERAP2 substrate channel
subsequence (residues 350-430), for both K392 (ancestral) and N392 variants.
Downstream: filter candidates by P1 position chemistry (Glu/Asp = salt bridge
hypothesis), then screen with Boltz-2 for structural validation.

Usage:
  # Local (CPU, ~5-10 min for 50 peptides):
  pip install transformers torch
  python scripts/pepmlm_generate.py

  # Colab (GPU, ~1-2 min):
  Copy this script into a Colab cell after !pip install transformers torch

Reference: Chen et al., Nature Biotechnology 2025
Model: https://huggingface.co/TianlaiChen/PepMLM-650M
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
from pathlib import Path
from datetime import datetime

import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
MODEL_ID = "TianlaiChen/PepMLM-650M"
NUM_PEPTIDES = 24          # candidates per variant (reduced for CPU)
PEPTIDE_LENGTHS = [8, 10, 12, 15]  # residues — V4 library uses 8-15mers
CHANNEL_START = 349        # 0-indexed, residue 350 (1-indexed)
CHANNEL_END = 430          # 0-indexed exclusive, residue 430 (1-indexed)
VARIANT_POS = 391          # 0-indexed, residue 392 (1-indexed)
TOP_K = 3                  # top-k sampling for diversity
TEMPERATURE = 1.2          # sampling temperature (>1 = more diverse)

# Full ERAP2 sequence (UniProt Q6P179)
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

# Output directory
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "results" / "pepmlm"

# Amino acids with salt-bridge potential (Glu, Asp) — the V4 hypothesis
SALT_BRIDGE_AA = {"E", "D"}
# Hydrophobic P1 residues (control — predicted N392-selective)
HYDROPHOBIC_P1 = {"L", "I", "V", "F", "W", "Y", "M"}


def make_variants(seq: str) -> dict[str, str]:
    """Create K392 and N392 channel subsequences."""
    # Verify position 392 is K in wild-type
    assert seq[VARIANT_POS] == "K", f"Expected K at pos {VARIANT_POS+1}, got {seq[VARIANT_POS]}"
    channel_k = seq[CHANNEL_START:CHANNEL_END]  # K392 (ancestral)
    seq_n = seq[:VARIANT_POS] + "N" + seq[VARIANT_POS + 1:]
    channel_n = seq_n[CHANNEL_START:CHANNEL_END]  # N392 variant
    return {"K392": channel_k, "N392": channel_n}


def generate_peptides(
    model, tokenizer, target_seq: str, pep_length: int, n_samples: int, device: str
) -> list[str]:
    """Generate peptide binders for a target sequence using PepMLM."""
    mask_token = tokenizer.mask_token
    # PepMLM format: [TARGET_SEQ][MASK_TOKENS] at C-terminus
    masked_seq = target_seq + mask_token * pep_length
    peptides = []

    for _ in range(n_samples):
        inputs = tokenizer(masked_seq, return_tensors="pt").to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        logits = outputs.logits[0]  # (seq_len, vocab_size)

        # Decode only the masked positions (last pep_length tokens before </s>)
        mask_positions = (inputs["input_ids"][0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0]

        peptide_tokens = []
        for pos in mask_positions:
            pos_logits = logits[pos] / TEMPERATURE
            probs = torch.softmax(pos_logits, dim=-1)
            # Top-k sampling
            top_k_probs, top_k_ids = torch.topk(probs, TOP_K)
            top_k_probs = top_k_probs / top_k_probs.sum()
            chosen = top_k_ids[torch.multinomial(top_k_probs, 1)]
            peptide_tokens.append(chosen.item())

        peptide = tokenizer.decode(peptide_tokens).replace(" ", "")
        # Filter: only keep valid single-letter AA sequences
        if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in peptide) and len(peptide) == pep_length:
            peptides.append(peptide)

    return list(set(peptides))  # deduplicate


def classify_peptide(seq: str) -> dict:
    """Classify peptide by P1 position and salt-bridge potential."""
    p1 = seq[0]  # N-terminal = P1 in substrate-mimic design
    return {
        "p1_residue": p1,
        "salt_bridge_candidate": p1 in SALT_BRIDGE_AA,
        "hydrophobic_p1": p1 in HYDROPHOBIC_P1,
        "charge_class": (
            "negative" if p1 in {"D", "E"} else
            "positive" if p1 in {"K", "R", "H"} else
            "hydrophobic" if p1 in HYDROPHOBIC_P1 else
            "polar"
        ),
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")
    print(f"Loading PepMLM ({MODEL_ID})...")

    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_ID).to(device)
    model.eval()
    print("Model loaded.\n")

    variants = make_variants(ERAP2_FULL)
    all_results = {}

    for variant_name, channel_seq in variants.items():
        print(f"=== Generating for {variant_name} (channel residues 350-430) ===")
        print(f"Channel sequence ({len(channel_seq)} aa): {channel_seq[:30]}...{channel_seq[-10:]}")
        print()

        variant_results = []
        for pep_len in PEPTIDE_LENGTHS:
            n_per_length = NUM_PEPTIDES // len(PEPTIDE_LENGTHS)
            print(f"  Generating {n_per_length} x {pep_len}-mers...", end=" ", flush=True)
            peptides = generate_peptides(
                model, tokenizer, channel_seq, pep_len, n_per_length * 3, device
            )
            # Keep top n_per_length unique
            peptides = peptides[:n_per_length]
            print(f"got {len(peptides)} unique")

            for pep in peptides:
                entry = {
                    "sequence": pep,
                    "length": len(pep),
                    "variant": variant_name,
                    **classify_peptide(pep),
                }
                variant_results.append(entry)

        all_results[variant_name] = variant_results
        print(f"  Total {variant_name}: {len(variant_results)} candidates\n")

    # --- Analysis ---
    print("=" * 60)
    print("SALT BRIDGE CANDIDATES (D/E at P1 — predicted K392-selective)")
    print("=" * 60)
    salt_bridge_hits = []
    for variant, results in all_results.items():
        hits = [r for r in results if r["salt_bridge_candidate"]]
        salt_bridge_hits.extend(hits)
        print(f"\n{variant}: {len(hits)} candidates with D/E at P1")
        for h in hits[:10]:
            print(f"  {h['sequence']}  (P1={h['p1_residue']}, {h['length']}aa)")

    print(f"\n{'=' * 60}")
    print("P1 DISTRIBUTION (both variants)")
    print("=" * 60)
    from collections import Counter
    all_peptides = [r for results in all_results.values() for r in results]
    p1_counts = Counter(r["p1_residue"] for r in all_peptides)
    for aa, count in p1_counts.most_common():
        marker = " <-- SALT BRIDGE" if aa in SALT_BRIDGE_AA else ""
        print(f"  {aa}: {count}{marker}")

    # --- Save ---
    output = {
        "metadata": {
            "model": MODEL_ID,
            "generated_at": datetime.now().isoformat(),
            "device": device,
            "channel_region": f"{CHANNEL_START+1}-{CHANNEL_END}",
            "variant_position": VARIANT_POS + 1,
            "num_peptides_requested": NUM_PEPTIDES,
            "peptide_lengths": PEPTIDE_LENGTHS,
            "temperature": TEMPERATURE,
            "top_k": TOP_K,
        },
        "variants": {
            name: {"channel_sequence": variants[name], "peptides": results}
            for name, results in all_results.items()
        },
        "salt_bridge_candidates": salt_bridge_hits,
        "p1_distribution": dict(p1_counts.most_common()),
    }

    out_path = OUT_DIR / "pepmlm_candidates.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")

    # Also save a simple FASTA for Boltz-2 screening
    fasta_path = OUT_DIR / "pepmlm_candidates.fasta"
    with open(fasta_path, "w") as f:
        for variant, results in all_results.items():
            for i, r in enumerate(results):
                f.write(f">pepmlm_{variant}_{i:03d}_{r['length']}aa_P1{r['p1_residue']}\n")
                f.write(f"{r['sequence']}\n")
    print(f"FASTA saved to {fasta_path}")

    print(f"\nNext steps:")
    print(f"  1. Review salt bridge candidates (D/E at P1)")
    print(f"  2. Add D-amino acid modifications (D-Glu, D-Asp) to top hits")
    print(f"  3. Screen with Boltz-2: boltz predict pepmlm_candidate + ERAP2_K392")
    print(f"  4. Compare binding to K392 vs N392 for selectivity")


if __name__ == "__main__":
    main()
