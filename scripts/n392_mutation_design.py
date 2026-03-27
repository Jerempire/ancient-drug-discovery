"""
N392 Mutation Design — Phase 2 of N392 selectivity enhancement.

Generates Boltz-2 YAMLs for point mutations at positions 37 (Asn) and 50 (Phe)
of the Y87A_Y89A binder, screening against 4 targets: ERAP2-K392, ERAP2-N392,
ERAP1, IRAP.

Contact analysis (Phase 1) found:
  - Pos 37 (N, Asn): 4.36A from K392 NZ via OD1 -> likely H-bond with K392 NH3+
  - Pos 50 (F, Phe): 3.84A from K392 CE via CD2 -> hydrophobic contact with Lys chain

Mutation rationale:
  K392 = Lys: +1 charge, ~6.5A side chain
  N392 = Asn: neutral, ~3.5A side chain (3A shorter)

  To disfavor K392 while favoring N392:
  - Bulky: W, F -> steric clash with K392's extended chain, fills N392 void
  - Positive: K, R -> +/+ repulsion with K392, H-bond with N392 amide C=O
  - H-bond complementary: Q -> reciprocal H-bonds with N392 amide (dual donor+acceptor)
  - Small: A -> remove any favorable K392 contact, neutral baseline

Usage:
  python scripts/n392_mutation_design.py

Output: data/results/n392_selectivity/boltz_yamls/ (ready for Vast.ai upload)
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
from pathlib import Path
from datetime import datetime

PROJECT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT / "data" / "results" / "n392_selectivity"
YAML_DIR = OUTPUT_DIR / "boltz_yamls"

# Y87A_Y89A binder sequence (92aa)
BINDER_WT = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

# ERAP2 full sequence (from n392_scaffold_search.py)
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

# Target crops (151aa each, residues 350-500 in full numbering)
ERAP2_K392_CROP = ERAP2_FULL[349:500]  # K at position 392 (index 391)
ERAP2_N392_CROP = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:500]

# ERAP1 and IRAP crops (from mpnn_redesign.py — validated in prior campaigns)
ERAP1_CROP = "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR"
IRAP_CROP = "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN"

TARGETS = {
    "erap2_k392": ERAP2_K392_CROP,
    "erap2_n392": ERAP2_N392_CROP,
    "erap1": ERAP1_CROP,
    "irap": IRAP_CROP,
}

# Mutation positions and candidates
# Pos 37 (0-indexed): Asn -> H-bonds with K392 NZ at 4.36A
# Pos 50 (0-indexed): Phe -> hydrophobic contact with K392 CE at 3.84A
MUTATIONS = {
    37: {
        "wt": "N",
        "rationale": "OD1 H-bonds K392 NZ at 4.36A",
        "mutants": {
            "W": "Bulky aromatic — fills N392 void, clashes K392 extended chain",
            "K": "+/+ repulsion with K392, H-bond donor to N392 C=O",
            "R": "+/+ repulsion with K392, guanidinium H-bonds N392",
            "Q": "Longer amide — better reach for N392 H-bonds",
            "D": "Negative charge — ATTRACTS K392 (negative control, should favor K)",
            "A": "Remove H-bond entirely — baseline for effect size",
        }
    },
    50: {
        "wt": "F",
        "rationale": "CD2 contacts K392 CE at 3.84A (hydrophobic)",
        "mutants": {
            "W": "Larger aromatic — more steric clash with K392 chain",
            "K": "Positive charge replaces aromatic — repels K392",
            "R": "Positive + long chain — maximum K392 repulsion",
            "N": "H-bond donor/acceptor — complements N392 amide",
            "Q": "Longer H-bond — reaches deeper toward N392",
            "A": "Remove hydrophobic contact — baseline",
        }
    },
}


def mutate_binder(seq, pos_0idx, new_aa):
    """Create mutant binder sequence."""
    return seq[:pos_0idx] + new_aa + seq[pos_0idx + 1:]


def write_yaml(target_seq, binder_seq, path):
    """Write Boltz-2 input YAML."""
    path.write_text(
        f"version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {target_seq}\n      msa: empty\n"
        f"  - protein:\n      id: B\n      sequence: {binder_seq}\n      msa: empty\n"
    )


def main():
    YAML_DIR.mkdir(parents=True, exist_ok=True)

    # Generate wildtype control YAMLs
    yaml_count = 0
    manifest = []

    print("=" * 70)
    print("N392 SELECTIVITY MUTATION PANEL")
    print("=" * 70)

    # Wildtype control
    print("\n--- WILDTYPE CONTROL ---")
    print(f"Binder: Y87A_Y89A (92aa)")
    for tgt_name, tgt_seq in TARGETS.items():
        name = f"wt_{tgt_name}"
        yaml_path = YAML_DIR / f"{name}.yaml"
        write_yaml(tgt_seq, BINDER_WT, yaml_path)
        manifest.append({"name": name, "position": "wt", "mutation": "wt", "target": tgt_name})
        yaml_count += 1
    print(f"  4 YAMLs (wildtype x 4 targets)")

    # Mutations
    for pos, info in MUTATIONS.items():
        wt_aa = info["wt"]
        print(f"\n--- POSITION {pos} ({wt_aa}, {info['rationale']}) ---")
        for mut_aa, rationale in info["mutants"].items():
            mut_seq = mutate_binder(BINDER_WT, pos, mut_aa)
            mut_label = f"{wt_aa}{pos + 1}{mut_aa}"  # 1-indexed convention (N38W, F51W, etc.)
            print(f"  {mut_label}: {rationale}")

            for tgt_name, tgt_seq in TARGETS.items():
                name = f"{mut_label}_{tgt_name}"
                yaml_path = YAML_DIR / f"{name}.yaml"
                write_yaml(tgt_seq, mut_seq, yaml_path)
                manifest.append({
                    "name": name,
                    "position": pos,
                    "mutation": mut_label,
                    "mutant_aa": mut_aa,
                    "target": tgt_name,
                })
                yaml_count += 1

    # Also test double mutants for the most promising combos
    print(f"\n--- DOUBLE MUTANTS (top combos) ---")
    double_mutants = [
        (37, "K", 50, "W", "K38+W51: positive repulsion + steric fill"),
        (37, "R", 50, "W", "R38+W51: guanidinium repulsion + steric fill"),
        (37, "K", 50, "K", "K38+K51: double positive — max K392 repulsion"),
        (37, "Q", 50, "N", "Q38+N51: double H-bond complementary to N392"),
    ]
    for pos1, aa1, pos2, aa2, rationale in double_mutants:
        mut_seq = mutate_binder(BINDER_WT, pos1, aa1)
        mut_seq = mutate_binder(mut_seq, pos2, aa2)
        wt1 = MUTATIONS[pos1]["wt"]
        wt2 = MUTATIONS[pos2]["wt"]
        mut_label = f"{wt1}{pos1+1}{aa1}_{wt2}{pos2+1}{aa2}"
        print(f"  {mut_label}: {rationale}")

        for tgt_name, tgt_seq in TARGETS.items():
            name = f"{mut_label}_{tgt_name}"
            yaml_path = YAML_DIR / f"{name}.yaml"
            write_yaml(tgt_seq, mut_seq, yaml_path)
            manifest.append({
                "name": name,
                "position": f"{pos1}+{pos2}",
                "mutation": mut_label,
                "target": tgt_name,
            })
            yaml_count += 1

    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    n_constructs = 1 + len(MUTATIONS[37]["mutants"]) + len(MUTATIONS[50]["mutants"]) + len(double_mutants)
    print(f"Constructs: {n_constructs} (1 WT + {len(MUTATIONS[37]['mutants'])} pos37 + {len(MUTATIONS[50]['mutants'])} pos50 + {len(double_mutants)} doubles)")
    print(f"Targets: {len(TARGETS)}")
    print(f"Total YAMLs: {yaml_count}")
    print(f"Estimated time on RTX 4090: ~{yaml_count * 30 // 60} min ({yaml_count} x ~30 sec)")
    print(f"Estimated cost: ~${yaml_count * 30 / 3600 * 0.37:.2f}")
    print(f"\nYAMLs saved to: {YAML_DIR}")

    # Save manifest
    manifest_data = {
        "generated_at": datetime.now().isoformat(),
        "binder_wt": BINDER_WT,
        "mutations": {str(k): v for k, v in MUTATIONS.items()},
        "double_mutants": [{"pos1": d[0], "aa1": d[1], "pos2": d[2], "aa2": d[3], "rationale": d[4]} for d in double_mutants],
        "targets": list(TARGETS.keys()),
        "total_yamls": yaml_count,
        "manifest": manifest,
    }
    manifest_path = OUTPUT_DIR / "mutation_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest_data, f, indent=2)
    print(f"Manifest saved to: {manifest_path}")

    # Print Vast.ai upload instructions
    print(f"\n{'='*70}")
    print("VAST.AI INSTRUCTIONS")
    print(f"{'='*70}")
    print(f"""
1. Launch RTX 4090 instance:
   python scripts/vast_launch.py search

2. Upload YAMLs:
   scp -P <PORT> -r {YAML_DIR}/* root@<HOST>:/workspace/yamls/

3. Run Boltz-2 batch:
   ssh -p <PORT> root@<HOST> 'cd /workspace && for y in yamls/*.yaml; do
     name=$(basename $y .yaml)
     boltz predict $y --out_dir results/$name --diffusion_samples 1
   done'

4. Download results:
   scp -P <PORT> -r root@<HOST>:/workspace/results/ {OUTPUT_DIR / 'boltz_results'}/

5. Analyze:
   python scripts/n392_mutation_analysis.py
""")


if __name__ == "__main__":
    main()
