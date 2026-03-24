"""
Generate Boltz-2 structures for Elite Four peptides against all 3 targets.
Then graft peptide poses into full protein structures for MD.

Step 1 (Vast.ai GPU): Boltz-2 predictions on cropped channels — saves CIF files
Step 2 (local CPU): Graft peptide coordinates into full protein PDBs

Usage:
  # On Vast.ai:
  python3 elite_four_structures.py predict

  # Locally after downloading CIFs:
  python3 scripts/elite_four_structures.py graft
"""
import subprocess
import json
import os
import sys
import glob

SEED = 42
DIFFUSION_SAMPLES = 3

ERAP2_FULL_SEQ = "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"

K392_CHANNEL = ERAP2_FULL_SEQ[349:450]
N392_CHANNEL = ERAP2_FULL_SEQ[349:391] + "N" + ERAP2_FULL_SEQ[392:450]

# ERAP1 channel (aligned residues 333-433)
ERAP1_FULL_SEQ = "MEPPRGPRPLLTLLLLALAPGAGASQDCNPLAHGITGHLSHRIRQNFGWNLDQHFSAEQFQEYMLRGKHIPVKAHKGLSHLTAKDMAQKVQTPFNQRIGVSDPMLALVAQSGPGLERSYVLLASAEVDNIPRSTAQSISLDDDTLQWFQLAQLVQHVIRQYSDTDRISSPSTASALSPSSALLWPQGFILANPDEFNAQKLCSHLGTPDPAHYTALTRYQLQNPKLNNLHKSMQPYILAQELDALNYENAFKFYELQFSHTDALSNFKPSASSAFYSQLMKYLNNLKDLLEEKEISRGIFQKLHLPKQGSLRKNKVKIGLLDIFALKFNTAELQTIENDWLSHFRFQQLSGVTDLMFRSPLMRFFTFYESHFLAMELYNKKFMKLYKNQFQPITKFDDTKLYHAVKSLNSMTAMGSNVFNLNIFNLESMGYAEVKDRIAEWQNFKQHGSSPDLTIEQSMRVLQKFDTFVKPRYSQMKNHPEEMDFFSYYSFSWKDDPMQCFTFNCNEHNAILKWYQDRKAMNVGFKEVELPSFQKFAENLRESILKVVMVYIFAAAQALSKELDVNAWKYANNDKEQGVITFTTDESNYLRVRVSHELMEFSYLFLDESYRQYVVQIFERLPYGDLRTLMKASMPVFAQDFADTASFLNKFSFPTLGFMSRTMVRDRNRGMAVISPHFYYSGSEAGITLHYVDMENPQELNDAIRFAWRNFFRQCFDIDSQAAIRHYELQFQYTSMSGTALPRGDSYALDGLVEGAMEVKEQFGTFIDALANYSLQNHGLPASCTFTVQKNASQALLEELRVLHLHREGGIYAGSQHFDSYEGLFNKRDLMVDGKMLPEKFHSNKFSESYRSVFVDDGMQNYLNQFLMQHSIHQHALTSGFLRDNFEHFPWMKEIQASGTSMRPMYHKPGLQLLSLGLDYDLVNAFPKDKEQYAFAGWLKELFLPLEGGFAAQ"
ERAP1_CHANNEL = ERAP1_FULL_SEQ[332:433]

# IRAP channel (aligned residues 444-544)
IRAP_FULL_SEQ = "MSKRKEFGALITLGLCVMACSAEEEERPHRDLADCSAEEYEQGATPQLFDTMHKFNDQVYAQRNPHTYHGTKAQIKEGLLDRFANRIQELEEDTSEDGSKDSSPFRCHHGKILTINGSVDMSMEQTAQEEFHKKTYRVDPEKKQRGDDNKCGDTCLSWLAQKFPHKQYHQFGIIDRQLKALPPMEKPRYGSDPAMAIRTKDFMTPNFVHYAGMGSPSEGQYKMTYREYAQRVVELYEEFLTAKDVLQKRFQFPDVFYSDNGTLYPQPKHHTLQNGQNFRRWEDHMDPWFVLDLFHKQAELKQGYLKKKDYMKLYKEKYQSIIKFNQPQKFEDDKRFYHAVKSLGDSKAMEAEYSFLNLNLFNIESMGYPEKLDILGKWHNLRMNVSFVDVFNIHDYKMYEIVEDILRHLIEFDDYRFVSIEFLMHDSKNFQPLNYDFEVYEYSYLSWRRPDMCFTFNCNHSHAIEQWYQRKFEMNSRLLEFSAEPNKSFVTKFDLEIKPNYTPMKNHPHDMYFFYSYANWKDDPLQCFTYNCNKQGTAILKWYRSRKSSNIGLKEIPMPSLQKFAHNLRDTIISAIMIYMFAAAEALTKELDTSGWDYANDTDNRVTFTIDEGSYLKARIPHELAHQWFGNLVTMEWWNDIYLKEGFAKYMELIAVNATYPELQFDDYYLNFCFEVISKDSINSSRPISKPEKTPAQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFAYRNAKNNDLWSSLSNSCLESDFTSGGVCHSDDKMTSSMLAFLEEHADIKEMMTTWTLQKGTPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTALVKFNVDSNGYYIVHYEGHGWDQLIIKLNQNHTLLRPKDRVGLMHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGFSYLESFYHMMDLRKISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSLSSSEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"
IRAP_CHANNEL = IRAP_FULL_SEQ[443:544]

# Elite Four peptides
ELITE_FOUR = {
    "VAGSAF": "VAGSAF",
    "IAFSAF": "IAFSAF",
    "VAWSAF": "VAWSAF",
    "FASGAV": "FASGAV",
}

# Targets
TARGETS = {
    "erap2k392": K392_CHANNEL,
    "erap1": ERAP1_CHANNEL,
    "irap": IRAP_CHANNEL,
}

OUT_BASE = "/workspace/results/elite_four"


def write_yaml(pep_name, target_name, pep_seq, channel_seq):
    yaml_dir = os.path.join(OUT_BASE, "yamls")
    os.makedirs(yaml_dir, exist_ok=True)
    name = "%s_vs_%s" % (pep_name, target_name)
    path = os.path.join(yaml_dir, "%s.yaml" % name)
    with open(path, "w") as f:
        f.write("version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: %s\n      msa: empty\n  - protein:\n      id: B\n      sequence: %s\n      msa: empty\n" % (channel_seq, pep_seq))
    return path, name


def run_boltz(yaml_path, out_dir):
    cmd = ["boltz", "predict", yaml_path, "--out_dir", out_dir,
           "--diffusion_samples", str(DIFFUSION_SAMPLES), "--seed", str(SEED), "--no_kernels"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def cmd_predict():
    """Run Boltz-2 predictions on Vast.ai — saves CIF structure files."""
    os.makedirs(OUT_BASE, exist_ok=True)
    total = len(ELITE_FOUR) * len(TARGETS)
    done = 0

    print("ELITE FOUR STRUCTURE GENERATION")
    print("=" * 60)
    print("%d peptides x %d targets = %d predictions" % (len(ELITE_FOUR), len(TARGETS), total))
    print("Saving CIF files for MD grafting\n")

    results = {}
    for pep_name, pep_seq in ELITE_FOUR.items():
        results[pep_name] = {}
        for target_name, channel_seq in TARGETS.items():
            done += 1
            yaml_path, run_name = write_yaml(pep_name, target_name, pep_seq, channel_seq)
            out_dir = os.path.join(OUT_BASE, run_name)
            print("[%d/%d] %s..." % (done, total, run_name), end=" ", flush=True)
            ok = run_boltz(yaml_path, out_dir)
            if ok:
                # Find best CIF
                cifs = sorted(glob.glob(os.path.join(out_dir, "**/*.cif"), recursive=True))
                confs = sorted(glob.glob(os.path.join(out_dir, "**/confidence_*.json"), recursive=True))
                iptm = 0
                if confs:
                    with open(confs[0]) as f:
                        c = json.load(f)
                    iptm = c.get("iptm", c.get("i_ptm", 0))
                results[pep_name][target_name] = {"iptm": iptm, "cif": cifs[0] if cifs else None, "n_cifs": len(cifs)}
                print("ipTM=%.3f, %d CIFs saved" % (iptm, len(cifs)))
            else:
                results[pep_name][target_name] = {"iptm": 0, "cif": None, "error": True}
                print("FAILED")

    # Save manifest
    manifest = os.path.join(OUT_BASE, "structure_manifest.json")
    with open(manifest, "w") as f:
        json.dump(results, f, indent=2)
    print("\nManifest saved to %s" % manifest)
    print("\nDownload with:")
    print("  scp -r -P <PORT> root@<HOST>:%s <local_path>" % OUT_BASE)


def cmd_graft():
    """Graft peptide poses from cropped-channel CIFs into full protein PDBs.

    Run locally after downloading CIFs from Vast.ai.

    Requires: BioPython
    """
    try:
        from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer, Selection
        from Bio.PDB.Chain import Chain
    except ImportError:
        print("ERROR: pip install biopython")
        sys.exit(1)

    import numpy as np

    # Paths (local)
    PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    STRUCT_DIR = os.path.join(PROJECT, "data", "structures")
    CIF_DIR = os.path.join(PROJECT, "data", "results", "v43_validation", "elite_four_cifs")
    OUT_DIR = os.path.join(PROJECT, "data", "results", "v43_validation", "md_starting_structures")
    os.makedirs(OUT_DIR, exist_ok=True)

    # Full protein structures
    FULL_STRUCTURES = {
        "erap2k392": os.path.join(STRUCT_DIR, "erap2_wt_alphafold.pdb"),
        "erap1": os.path.join(STRUCT_DIR, "erap1_wt_alphafold.pdb"),
        "irap": os.path.join(STRUCT_DIR, "irap_5MJ6_experimental.pdb"),
    }

    # Channel residue ranges in full structures (for alignment)
    CHANNEL_RANGES = {
        "erap2k392": (350, 450),  # 1-indexed
        "erap1": (333, 433),
        "irap": (444, 544),
    }

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    # Find all CIF files
    if not os.path.exists(CIF_DIR):
        print("ERROR: CIF directory not found: %s" % CIF_DIR)
        print("Download CIFs from Vast.ai first:")
        print("  scp -r -P <PORT> root@<HOST>:/workspace/results/elite_four %s" % CIF_DIR)
        sys.exit(1)

    count = 0
    for pep_name in ELITE_FOUR:
        for target_name in TARGETS:
            run_name = "%s_vs_%s" % (pep_name, target_name)
            # Find CIF
            cif_pattern = os.path.join(CIF_DIR, run_name, "**", "*_model_0.cif")
            cifs = glob.glob(cif_pattern, recursive=True)
            if not cifs:
                print("SKIP %s: no CIF found" % run_name)
                continue

            cif_path = cifs[0]
            full_pdb_path = FULL_STRUCTURES[target_name]

            if not os.path.exists(full_pdb_path):
                print("SKIP %s: full structure not found: %s" % (run_name, full_pdb_path))
                continue

            print("Grafting %s..." % run_name, end=" ")

            # Load structures
            cif_struct = cif_parser.get_structure("cif", cif_path)
            full_struct = pdb_parser.get_structure("full", full_pdb_path)

            # Get chains from CIF (A=channel, B=peptide)
            cif_model = cif_struct[0]
            cif_chains = list(cif_model.get_chains())
            if len(cif_chains) < 2:
                print("SKIP: CIF has < 2 chains")
                continue

            cif_channel = cif_chains[0]  # Chain A = cropped channel
            cif_peptide = cif_chains[1]  # Chain B = peptide

            # Get channel residues from full structure for alignment
            full_model = full_struct[0]
            full_chain = list(full_model.get_chains())[0]
            start, end = CHANNEL_RANGES[target_name]

            # Get CA atoms for alignment
            cif_cas = []
            full_cas = []
            for res in cif_channel.get_residues():
                if "CA" in res:
                    cif_cas.append(res["CA"])

            full_channel_res = [r for r in full_chain.get_residues()
                               if start <= r.get_id()[1] <= end and "CA" in r]

            # Match by count (both should be ~100 residues)
            n = min(len(cif_cas), len(full_channel_res))
            if n < 10:
                print("SKIP: too few matching residues (%d)" % n)
                continue

            cif_cas = cif_cas[:n]
            full_cas = [r["CA"] for r in full_channel_res[:n]]

            # Superimpose: align CIF channel onto full protein channel
            sup = Superimposer()
            sup.set_atoms(full_cas, cif_cas)
            sup.apply(cif_model.get_atoms())

            print("RMSD=%.2f" % sup.rms, end=" ")

            # Now the CIF peptide is in the full protein's coordinate frame
            # Add peptide as chain B to the full structure
            # First remove any existing chain B
            for chain in list(full_model.get_chains()):
                if chain.id == "B":
                    full_model.detach_child(chain.id)

            # Rename peptide chain to B and add
            cif_peptide.id = "B"
            full_model.add(cif_peptide)

            # Save
            out_path = os.path.join(OUT_DIR, "%s.pdb" % run_name)
            io.set_structure(full_struct)
            io.save(out_path)
            print("-> %s" % out_path)
            count += 1

    print("\nGrafted %d structures to %s" % (count, OUT_DIR))
    print("These are ready for OpenMM MD simulations.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python elite_four_structures.py [predict|graft]")
        sys.exit(1)

    if sys.argv[1] == "predict":
        cmd_predict()
    elif sys.argv[1] == "graft":
        cmd_graft()
    else:
        print("Unknown command: %s" % sys.argv[1])
