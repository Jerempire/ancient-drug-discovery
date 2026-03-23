"""
Create hybrid peptides: PepMLM scaffolds + salt bridge P1 substitutions.

Takes PepMLM-generated scaffolds and substitutes P1 (position 1) with:
  - E (Glu) → predicted K392-selective via salt bridge with K392 lysine
  - D (Asp) → predicted K392-selective (shorter side chain variant)
  - L (Leu) → predicted N392-selective control (hydrophobic)
  - A (Ala) → neutral control

This creates a combinatorial library merging:
  1. PepMLM's learned binding scaffolds (model knows what fits the channel)
  2. V4 salt bridge hypothesis (engineered P1 for variant selectivity)

Output: FASTA + YAML for Boltz-2 docking on Vast.ai
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
from pathlib import Path
from datetime import datetime

# Paths
PROJECT = Path(__file__).resolve().parent.parent
PEPMLM_RESULTS = PROJECT / "data" / "results" / "pepmlm" / "pepmlm_candidates.json"
OUT_DIR = PROJECT / "data" / "results" / "pepmlm_hybrids"

# P1 substitutions matching V4 library design
P1_SUBS = {
    "E": {"label": "D-Glu", "hypothesis": "K392-selective (salt bridge)", "priority": 1},
    "D": {"label": "D-Asp", "hypothesis": "K392-selective (shorter salt bridge)", "priority": 1},
    "L": {"label": "D-Leu", "hypothesis": "N392-selective control (hydrophobic)", "priority": 2},
    "A": {"label": "D-Ala", "hypothesis": "neutral control", "priority": 3},
}

# ERAP2 sequences for Boltz-2 YAML (channel + flanking context)
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

# N392 variant
ERAP2_N392 = ERAP2_FULL[:391] + "N" + ERAP2_FULL[392:]


def select_diverse_scaffolds(peptides: list[dict], max_per_length: int = 2) -> list[dict]:
    """Pick diverse scaffolds: top 2 per length, prefer different P1 classes."""
    from collections import defaultdict
    by_length = defaultdict(list)
    for p in peptides:
        by_length[p["length"]].append(p)

    selected = []
    for length, peps in sorted(by_length.items()):
        # Sort by diversity: prefer different charge classes
        seen_classes = set()
        for p in peps:
            if p["charge_class"] not in seen_classes and len(selected) < len(by_length) * max_per_length:
                selected.append(p)
                seen_classes.add(p["charge_class"])
            if len([s for s in selected if s["length"] == length]) >= max_per_length:
                break
        # Fill remaining slots if not enough diversity
        for p in peps:
            if p not in selected and len([s for s in selected if s["length"] == length]) < max_per_length:
                selected.append(p)

    return selected


def make_hybrid(scaffold: dict, new_p1: str) -> dict:
    """Replace P1 (first residue) with new amino acid."""
    old_seq = scaffold["sequence"]
    new_seq = new_p1 + old_seq[1:]
    return {
        "sequence": new_seq,
        "original_scaffold": old_seq,
        "original_p1": old_seq[0],
        "new_p1": new_p1,
        "length": len(new_seq),
        "variant_target": scaffold["variant"],
        **P1_SUBS[new_p1],
    }


def write_boltz_yaml(hybrid: dict, erap2_seq: str, variant_label: str, out_dir: Path) -> Path:
    """Write Boltz-2 input YAML for a peptide-ERAP2 complex."""
    name = f"hybrid_{hybrid['new_p1']}_{hybrid['original_scaffold'][:6]}_{variant_label}"
    yaml_path = out_dir / f"{name}.yaml"

    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {erap2_seq}
  - protein:
      id: B
      sequence: {hybrid['sequence']}
      msa: empty
"""
    yaml_path.write_text(yaml_content)
    return yaml_path


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(PEPMLM_RESULTS) as f:
        data = json.load(f)

    # Select diverse scaffolds from K392 results (primary target)
    k392_peptides = data["variants"]["K392"]["peptides"]
    scaffolds = select_diverse_scaffolds(k392_peptides, max_per_length=2)

    print(f"Selected {len(scaffolds)} diverse scaffolds from K392 PepMLM results:")
    for s in scaffolds:
        print(f"  {s['sequence']} ({s['length']}aa, P1={s['p1_residue']}, {s['charge_class']})")

    # Generate hybrids
    hybrids = []
    for scaffold in scaffolds:
        for p1_aa in P1_SUBS:
            hybrid = make_hybrid(scaffold, p1_aa)
            hybrids.append(hybrid)

    print(f"\nGenerated {len(hybrids)} hybrid peptides ({len(scaffolds)} scaffolds x {len(P1_SUBS)} P1 subs)")

    # Priority 1 (salt bridge candidates) for immediate Boltz-2
    priority_1 = [h for h in hybrids if h["priority"] == 1]
    print(f"Priority 1 (salt bridge, E/D at P1): {len(priority_1)} peptides")

    # Save all hybrids
    output = {
        "metadata": {
            "generated_at": datetime.now().isoformat(),
            "source": "pepmlm_candidates.json",
            "num_scaffolds": len(scaffolds),
            "p1_substitutions": {k: v["label"] for k, v in P1_SUBS.items()},
            "total_hybrids": len(hybrids),
        },
        "scaffolds": scaffolds,
        "hybrids": hybrids,
        "priority_1_salt_bridge": priority_1,
    }

    json_path = OUT_DIR / "hybrid_candidates.json"
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {json_path}")

    # Save FASTA (all hybrids)
    fasta_path = OUT_DIR / "hybrid_candidates.fasta"
    with open(fasta_path, "w") as f:
        for i, h in enumerate(hybrids):
            name = f"hybrid_{i:03d}_P1{h['new_p1']}_{h['label']}_{h['length']}aa"
            f.write(f">{name}\n{h['sequence']}\n")
    print(f"FASTA saved to {fasta_path}")

    # Write Boltz-2 YAMLs for priority 1 (salt bridge) x both variants
    boltz_dir = OUT_DIR / "boltz_inputs"
    boltz_dir.mkdir(exist_ok=True)
    yaml_count = 0
    for h in priority_1:
        write_boltz_yaml(h, ERAP2_FULL, "K392", boltz_dir)
        write_boltz_yaml(h, ERAP2_N392, "N392", boltz_dir)
        yaml_count += 2
    print(f"Boltz-2 YAMLs: {yaml_count} files in {boltz_dir}")

    # Summary table
    print(f"\n{'='*70}")
    print(f"HYBRID LIBRARY SUMMARY")
    print(f"{'='*70}")
    print(f"{'Scaffold':<20} {'P1→E(Glu)':<12} {'P1→D(Asp)':<12} {'P1→L(Leu)':<12} {'P1→A(Ala)':<12}")
    print(f"{'-'*70}")
    for scaffold in scaffolds:
        row = f"{scaffold['sequence']:<20}"
        for p1 in ["E", "D", "L", "A"]:
            hybrid_seq = p1 + scaffold['sequence'][1:]
            row += f" {hybrid_seq:<12}"
        print(row)

    print(f"\n{'='*70}")
    print(f"NEXT STEPS")
    print(f"{'='*70}")
    print(f"  1. Upload {boltz_dir} to Vast.ai")
    print(f"  2. Run Boltz-2: boltz predict <yaml> --diffusion_samples 3 --seed 42")
    print(f"  3. Compare K392 vs N392 ipTM for each hybrid")
    print(f"  4. Salt bridge hits (E/D at P1 + K392-selective) → synthesis candidates")
    print(f"  5. Merge with V4 docking results for final ranking")


if __name__ == "__main__":
    main()
