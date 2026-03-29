#!/usr/bin/env python3
"""Generate Boltz-2 input YAMLs for cleavage-resistant V4 peptide variants.

Tests whether D-amino acid and N-methylated modifications at P1
preserve ERAP2 channel binding while resisting catalytic cleavage.

Modifications tested:
  1. D-Val at P1  (dVAGSAF)  — CCD: DVA — chirality flip blocks zinc orientation
  2. D-Ile at P1  (dIAFSAF)  — CCD: DIL — same for best K392-selective peptide
  3. D-Phe at P1  (dFAGSAF)  — CCD: DPN — aromatic D-amino acid control
  4. N-Me-Val P1  (NMeVAGSAF)— CCD: MVA — blocks backbone NH for catalysis
  5. N-Me-Ile P1  (NMeIAFSAF)— CCD: MEI — same for IAFSAF
  6. L-Val control (VAGSAF)   — no mod  — baseline (gets cleaved)
  7. L-Ile control (IAFSAF)   — no mod  — baseline (gets cleaved)

Each peptide is tested against 4 targets:
  - ERAP2 K392 (primary target)
  - ERAP2 N392 (variant)
  - ERAP1 (off-target — should NOT bind)
  - IRAP (off-target — should NOT bind)

Run on Vast.ai:
  boltz predict v4/cleavage_resistant/boltz2_inputs/ --out_dir v4/cleavage_resistant/boltz2_out/
"""

from pathlib import Path

OUTPUT_DIR = Path("v4/cleavage_resistant/boltz2_inputs")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Full ERAP2 K392 sequence
ERAP2_K392 = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"

# ERAP2 N392 — K→N at position 392 (1-indexed)
ERAP2_N392 = ERAP2_K392[:391] + "N" + ERAP2_K392[392:]

# ERAP1 full sequence (from UniProt Q9NZ08)
ERAP1 = "MGSRLFILFTLAHIFTGAFSHEDFEMTPKIFNHQIKGFYCLTAILPQICICSQFAVSSSYQFSEDPGAFPPATNADRFPWQEIQIPGVVIPLHNDLFVHPNLTSLDFVASEKIEVLVSNTTQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKASFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDNLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"

# IRAP full sequence (from UniProt Q9UIQ6, truncated to match prior runs if needed)
IRAP = "MKFSILTLFNLCFFSAQAQNSDCNRQKENASNIHSREFETKIDGRQKDWALPGTFQTFEGFLKFESFQDAFQHFFAQTCNETVPFLKDALIKLTDYYKVIGKFNLAALEKNCVDPAKITRDAIRILKKTSFTTTTLKQSSKLSTHFFDIVEGYDTSESSELKPNAETDKMECFQKYWDTLLDQLNTQFHKVYAEELKWEPLRSQELMPKMTDSMKTLNISSWAMDEAGLVEAYLKKNDISKPVTDDKVIFRFANDSYMATTQFGKTSISLTFHNDINYAQVYFSILQPSFSIEPLVPIGSVRVEIMGTHPAGSILEKDDKKGFCVLLDYNIQQKPFVWFLEGQFRHENYPTLHIEFDEKVWTMDKYNPHWSYIHCPDTKSLVERGDNQTLWEDKTIYIFSKDIMSNKGFGDRAELEDPFENYSYRTTLMVHVQDSHIKPLNVSWMTEDDQFPDNYEFGEQGISVLMALDKEYAGKDCFTNCPVLYVANEAPYNRKRSLQIKAYNDKIWYVKNALAALNLYQIKQNQVDYKLINQFLPNPHQFTNASPQRDFRQFQAAGLQNAGFKFGEVRTIKPILQVVHTFEHKKASDQSSSFWRDPELIQPSLIMYELGTLYHKWTEEEQAKWISRYITDTIHYLENPKNGFVLQIDATTTKELKFPNKFNMNNFQHPAGQFIVNPDEVNDPETFPNYTLHYNYNEQDKPNLHLSAEEFMKAFESYMDSAIASEEFIKLGKNPNKEKVRYWDYEKRYGKLTKQFEGFIKAIDKYKSIDADEHIEHAIFTGTHHALTDLLNRVYPMSGSFREKSTLDQIGPCLNKENPSEQFMSIYNKSIYWKNKSIQEGEQYPTIQTLALENYYGKTIEWYGHSNRWSELQPVDVLQFYQNRVTFLDLRFNPLSTHLISHTSYVYNLHRGLEPAWKGRVPAQEKSYEEQIIAQFQALKQFFNHEAYRQTFRSCQYICLDSHALRHFHFQTIDKENSAAYTAFQIVSELVNFYCDYAEETQDDLQLFAHVLRNLNGTYQILKTFTQLHAPDTFYSNVQAFPLSATMDLQKALTLSEALYIVSNEFKKDDYLQSLTALNVQNFVMKSYAFDEQIRGFRTTHLDIALQKLHDIYNKMTADESYRSILPGLQAMAPATDAQISYMHFQNVAPMAGQVHWIKTTFSERKEQFQFTLEYAEMAKKMALEEIKDHKYKQFEEELKNLKIDFNELGRIPQYGEQLLPRNNKYPNLFVLGKGIDYDKEQARQYLKLLVRSYLGQFNWEKCEAQNLPTKTAPAELPTSTHEKKLNVNLAAYMQAFCSTVNLPLTY"

TARGETS = {
    "K392": ERAP2_K392,
    "N392": ERAP2_N392,
    "ERAP1": ERAP1,
    "IRAP": IRAP,
}

# Peptide designs
# Format: (name, sequence, modifications_list)
# modifications_list: [(position_1indexed, ccd_code), ...]
PEPTIDES = [
    # Controls (unmodified, will be cleaved)
    ("ctrl_VAGSAF", "VAGSAF", []),
    ("ctrl_IAFSAF", "IAFSAF", []),

    # D-amino acid at P1 (position 1) — zinc can't orient D-residue for hydrolysis
    ("dVal_dVAGSAF", "VAGSAF", [(1, "DVA")]),      # D-Valine
    ("dIle_dIAFSAF", "IAFSAF", [(1, "DIL")]),      # D-Isoleucine
    ("dPhe_dFAGSAF", "FAGSAF", [(1, "DPN")]),      # D-Phenylalanine (aromatic control)

    # N-methylation at P1 — blocks backbone NH needed for zinc catalysis
    ("NMe_NMeVAGSAF", "VAGSAF", [(1, "MVA")]),     # N-methyl-Valine
    ("NMe_NMeIAFSAF", "IAFSAF", [(1, "MEI")]),     # N-methyl-Isoleucine

    # D-amino acid at P1 + P6 — double protection
    ("dVal_dPhe_dVAGSAdF", "VAGSAF", [(1, "DVA"), (6, "DPN")]),

    # All-D retro-inverso (reversed sequence, all D) — mimics L-peptide shape
    ("retroinv_FASgav", "FASGAV", [
        (1, "DPN"),  # D-Phe
        (2, "DAL"),  # D-Ala
        (3, "DSN"),  # D-Ser
        (4, "GL3"),  # D-Gly (glycine is achiral but CCD has GL3)
        (5, "DAL"),  # D-Ala
        (6, "DVA"),  # D-Val
    ]),
]


def write_yaml(name: str, target_name: str, target_seq: str,
               peptide_seq: str, modifications: list, output_dir: Path):
    """Write a Boltz-2 input YAML file."""
    fname = f"{name}_{target_name}.yaml"
    filepath = output_dir / fname

    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {target_seq}",
        "      msa: empty",
        "  - protein:",
        "      id: B",
        f"      sequence: {peptide_seq}",
    ]

    if modifications:
        lines.append("      modifications:")
        for pos, ccd in modifications:
            lines.append(f"        - position: {pos}")
            lines.append(f"          ccd: {ccd}")

    lines.append("      msa: empty")
    lines.append("")

    filepath.write_text("\n".join(lines))
    return filepath


def main():
    total = 0
    for pep_name, pep_seq, mods in PEPTIDES:
        for target_name, target_seq in TARGETS.items():
            write_yaml(pep_name, target_name, target_seq, pep_seq, mods, OUTPUT_DIR)
            total += 1

    print(f"Generated {total} Boltz-2 input files in {OUTPUT_DIR}/")
    print()
    print("Peptides:")
    for pep_name, pep_seq, mods in PEPTIDES:
        mod_str = ", ".join(f"pos{p}={c}" for p, c in mods) if mods else "none (control)"
        print(f"  {pep_name:30s} seq={pep_seq:8s} mods={mod_str}")
    print()
    print(f"Targets: {', '.join(TARGETS.keys())}")
    print(f"Total predictions: {total}")
    print()
    print("Run on Vast.ai:")
    print("  boltz predict v4/cleavage_resistant/boltz2_inputs/ --out_dir v4/cleavage_resistant/boltz2_out/ --recycling_steps 3 --diffusion_samples 3")
    print()
    print("Expected cost: ~$0.50-1.00 on RTX 4090 (~20 min)")
    print()
    print("Key comparison:")
    print("  If dVal_dVAGSAF_K392 ipTM ≈ ctrl_VAGSAF_K392 ipTM → modification preserves binding → SYNTHESIZE")
    print("  If dVal_dVAGSAF_K392 ipTM << ctrl_VAGSAF_K392 ipTM → modification kills binding → try next mod")


if __name__ == "__main__":
    main()
