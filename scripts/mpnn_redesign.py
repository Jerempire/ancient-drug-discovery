"""
ProteinMPNN interface redesign for IRAP selectivity.

Takes the best ERAP2-binder complex, redesigns the 17 conserved-only
interface positions on the binder while keeping divergent contacts and
core residues fixed. Then builds c2_eaaaa3-format constructs and
generates Boltz-2 YAMLs for screening.

Usage: python mpnn_redesign.py (run on Vast.ai after uploading CIF)
"""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import json
import os
import subprocess
from pathlib import Path

# ============================================================
# Config
# ============================================================

CIF_PATH = "/workspace/complex.cif"  # uploaded from local
OUTPUT_DIR = "/workspace/mpnn_output"
YAML_DIR = "/workspace/redesign_yamls"
N_DESIGNS = 100  # ProteinMPNN samples
TEMPERATURE = 0.1  # low = conservative, close to native
TOP_K = 10  # top designs to screen with Boltz-2

# Binder = chain B, residues 22-113 (in C2 numbering)
# Parent binder positions (1-indexed in parent): 1-92
# C2 positions = parent + 21

# Conserved-only contact positions (parent numbering) — REDESIGNABLE
CONSERVED_ONLY = [3, 5, 6, 11, 19, 23, 37, 38, 39, 40, 41, 51, 56, 89, 90, 91, 92]

# Divergent-only contact positions (parent numbering) — KEEP FIXED
DIVERGENT_ONLY = [17, 21, 28, 35]

# All other binder positions — KEEP FIXED (core fold)

# Target sequences for Boltz-2 screening
ERAP2 = "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
ERAP1 = "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR"
IRAP = "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN"

LINKER = "EAAAA" * 3
CARGO = "VAGSAF"

# Parent binder sequence (V19W version)
PARENT_BINDER = "DIRHYFKSLEEYLKNLPKWVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"


def step1_convert_cif_to_pdb():
    """Convert CIF to PDB using gemmi."""
    print("Step 1: Converting CIF to PDB...")
    import gemmi
    st = gemmi.read_structure(CIF_PATH)
    pdb_path = "/workspace/complex.pdb"
    st.write_pdb(pdb_path)
    print(f"  Written: {pdb_path}")
    return pdb_path


def step2_prepare_mpnn_inputs(pdb_path):
    """Create ProteinMPNN input JSON with fixed/redesignable positions."""
    print("Step 2: Preparing ProteinMPNN inputs...")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Build fixed positions: everything EXCEPT the conserved-only contacts on chain B
    # In the PDB, chain B residues are numbered 1-113 (C2 numbering)
    # Conserved-only in C2 numbering = parent_pos + 21
    redesignable_c2 = [p + 21 for p in CONSERVED_ONLY]

    # All chain B positions
    all_b_positions = list(range(1, 114))  # 1-113

    # Fixed = everything NOT in redesignable
    fixed_b = [p for p in all_b_positions if p not in redesignable_c2]

    # Chain A is fully fixed (it's the target)
    # We need to get chain A residue numbers from the PDB
    import gemmi
    st = gemmi.read_structure(pdb_path)
    model = st[0]
    chain_a = model.find_chain("A")
    all_a_positions = [res.seqid.num for res in chain_a]

    # Write chains to process
    chains_path = f"{OUTPUT_DIR}/chains_to_design.json"
    with open(chains_path, "w") as f:
        json.dump({"complex": ["B"]}, f)  # only redesign chain B

    # Write fixed positions
    fixed_path = f"{OUTPUT_DIR}/fixed_positions.json"
    fixed_dict = {"complex": {"B": fixed_b}}
    with open(fixed_path, "w") as f:
        json.dump(fixed_dict, f)

    # Write PDB list
    pdb_list = f"{OUTPUT_DIR}/pdb_list.txt"
    with open(pdb_list, "w") as f:
        f.write(pdb_path + "\n")

    print(f"  Redesignable positions (C2 numbering): {redesignable_c2}")
    print(f"  Fixed chain B positions: {len(fixed_b)} of {len(all_b_positions)}")
    print(f"  Chain A positions (all fixed): {len(all_a_positions)}")

    return chains_path, fixed_path, pdb_list


def step3_run_proteinmpnn(pdb_path, chains_path, fixed_path):
    """Run ProteinMPNN to generate redesigned sequences."""
    print(f"Step 3: Running ProteinMPNN ({N_DESIGNS} designs, T={TEMPERATURE})...")

    cmd = [
        "python", "/workspace/ProteinMPNN/protein_mpnn_run.py",
        "--pdb_path", pdb_path,
        "--chain_id_jsonl", chains_path,
        "--fixed_positions_jsonl", fixed_path,
        "--out_folder", OUTPUT_DIR,
        "--num_seq_per_target", str(N_DESIGNS),
        "--sampling_temp", str(TEMPERATURE),
        "--seed", "42",
        "--batch_size", "1",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    print(f"  stdout: {result.stdout[-500:]}")
    if result.returncode != 0:
        print(f"  stderr: {result.stderr[-500:]}")
        return None

    # Find output FASTA
    fasta_files = list(Path(OUTPUT_DIR).rglob("*.fa"))
    if not fasta_files:
        fasta_files = list(Path(OUTPUT_DIR).rglob("*.fasta"))
    if not fasta_files:
        print("  ERROR: No FASTA output found")
        return None

    print(f"  Output: {fasta_files[0]}")
    return fasta_files[0]


def step4_parse_and_rank(fasta_path):
    """Parse ProteinMPNN output and rank by score."""
    print("Step 4: Parsing and ranking designs...")

    designs = []
    current_header = None
    current_seq = ""

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    # Parse score from header
                    score = 0.0
                    for part in current_header.split(","):
                        if "score" in part.lower():
                            try:
                                score = float(part.split("=")[-1].strip())
                            except:
                                pass
                    designs.append({
                        "header": current_header,
                        "sequence": current_seq,
                        "score": score,
                    })
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line

    # Don't forget last entry
    if current_header and current_seq:
        score = 0.0
        for part in current_header.split(","):
            if "score" in part.lower():
                try:
                    score = float(part.split("=")[-1].strip())
                except:
                    pass
        designs.append({
            "header": current_header,
            "sequence": current_seq,
            "score": score,
        })

    # Sort by score (lower = better for ProteinMPNN)
    designs.sort(key=lambda x: x["score"])

    print(f"  Total designs: {len(designs)}")
    print(f"  Score range: {designs[0]['score']:.3f} to {designs[-1]['score']:.3f}")

    # Show top designs with mutations
    print(f"\n  Top {TOP_K} designs:")
    for i, d in enumerate(designs[:TOP_K]):
        seq = d["sequence"]
        # Extract binder portion (chain B part of the sequence)
        # ProteinMPNN outputs the designed chain only
        if len(seq) >= 92:
            binder_seq = seq[-92:]  # last 92 chars = binder body
        else:
            binder_seq = seq

        # Count mutations vs parent
        mutations = []
        for j, (orig, new) in enumerate(zip(PARENT_BINDER, binder_seq)):
            if orig != new:
                mutations.append(f"{orig}{j+1}{new}")

        print(f"    Design {i+1}: score={d['score']:.3f}, {len(mutations)} mutations: {', '.join(mutations[:10])}")

    return designs[:TOP_K]


def step5_build_yamls(top_designs):
    """Build Boltz-2 YAMLs for top designs."""
    print(f"Step 5: Building Boltz-2 YAMLs for {len(top_designs)} designs...")

    os.makedirs(YAML_DIR, exist_ok=True)
    targets = {"erap2k392": ERAP2, "erap1": ERAP1, "irap": IRAP}

    count = 0
    for i, design in enumerate(top_designs):
        seq = design["sequence"]
        # Build full c2_eaaaa3 construct with redesigned binder
        if len(seq) >= 92:
            binder_body = seq[-92:]
        else:
            binder_body = seq

        full_seq = CARGO + LINKER + binder_body

        for tname, tseq in targets.items():
            fname = f"{YAML_DIR}/mpnn{i+1:02d}_vs_{tname}.yaml"
            yaml_str = f"""version: 1
sequences:
  - protein:
      id: A
      msa: empty
      sequence: {tseq}
  - protein:
      id: B
      msa: empty
      sequence: {full_seq}
"""
            with open(fname, "w") as f:
                f.write(yaml_str)
            count += 1

    print(f"  Generated {count} YAMLs ({len(top_designs)} designs x 3 targets)")
    return count


if __name__ == "__main__":
    pdb_path = step1_convert_cif_to_pdb()
    chains_path, fixed_path, pdb_list = step2_prepare_mpnn_inputs(pdb_path)
    print("\nReady for ProteinMPNN. Run step3 after installing ProteinMPNN.")
    print(f"  PDB: {pdb_path}")
    print(f"  Chains: {chains_path}")
    print(f"  Fixed: {fixed_path}")
