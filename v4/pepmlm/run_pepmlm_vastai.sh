#!/bin/bash
# V4/PepMLM: Generate + score + cross-filter peptides on Vast.ai GPU
# Then run Boltz-2 validation on top K392-selective candidates
#
# Usage:
#   1. Launch Vast.ai instance (RTX 3090, pytorch image)
#   2. Upload this folder: scp -P <PORT> -r v4/pepmlm/ root@<HOST>:/workspace/pepmlm/
#   3. Upload ERAP2 structures: scp -P <PORT> v4/inputs/*.pdb root@<HOST>:/workspace/pepmlm/
#   4. SSH in and run: bash /workspace/pepmlm/run_pepmlm_vastai.sh
set -e

echo "============================================"
echo "  V4/PepMLM Pipeline — $(date)"
echo "============================================"

# ── STEP 1: Install deps ─────────────────────────────────────────────────────
echo "=== Installing dependencies ==="
pip install -q transformers torch pandas numpy biopython 2>&1 | tail -3
echo "Dependencies installed"

# ── STEP 2: Generate peptides ─────────────────────────────────────────────────
echo ""
echo "=== PepMLM Generation ==="
cd /workspace/pepmlm

# Run with both narrow and wide context for comparison
# V4 DiffPepDock finding: 9-mers universally prefer N392.
# Only 10+ mers reach the K392 salt bridge at position 392.
# Focus generation on 10-13 mers to maximize K392-selective hits.
echo "--- Narrow context (res 388-410) ---"
python generate_peptides.py \
    --num_peptides 400 \
    --top_k 3 \
    --peptide_lengths 10,11,12,13 \
    --device cuda \
    --top_n_boltz2 30

# Rename narrow results
mv results/pepmlm_selective_hits.csv results/pepmlm_selective_narrow.csv
mv results/pepmlm_all_for_boltz2.fasta results/pepmlm_boltz2_narrow.fasta

echo ""
echo "--- Wide context (res 350-430) ---"
python generate_peptides.py \
    --num_peptides 400 \
    --top_k 3 \
    --peptide_lengths 10,11,12,13 \
    --device cuda \
    --use_wide_context \
    --top_n_boltz2 30

mv results/pepmlm_selective_hits.csv results/pepmlm_selective_wide.csv
mv results/pepmlm_all_for_boltz2.fasta results/pepmlm_boltz2_wide.fasta

# ── STEP 3: Merge and deduplicate top candidates ─────────────────────────────
echo ""
echo "=== Merging narrow + wide candidates ==="
python3 << 'MERGE'
from pathlib import Path

narrow = Path("results/pepmlm_boltz2_narrow.fasta")
wide = Path("results/pepmlm_boltz2_wide.fasta")
merged = Path("results/pepmlm_boltz2_merged.fasta")

seen = set()
entries = []
for fasta in [narrow, wide]:
    if not fasta.exists():
        continue
    lines = fasta.read_text().strip().split("\n")
    for i in range(0, len(lines), 2):
        header, seq = lines[i], lines[i + 1]
        if seq not in seen:
            seen.add(seq)
            entries.append((header, seq))

with open(merged, "w") as f:
    for i, (header, seq) in enumerate(entries):
        f.write(f">pepmlm_{i:03d}|{header.split('|', 1)[-1] if '|' in header else ''}\n")
        f.write(f"{seq}\n")

print(f"Merged: {len(entries)} unique peptides -> {merged}")
MERGE

# ── STEP 4: Boltz-2 validation (optional, if boltz2 is installed) ─────────
if command -v boltz &> /dev/null || python3 -c "import boltz" 2>/dev/null; then
    echo ""
    echo "=== Boltz-2 Validation ==="
    echo "Boltz-2 found — running validation on top candidates"
    echo "NOTE: This step is optional. You can also download results and"
    echo "      run Boltz-2 separately with your existing pipeline."

    # Create Boltz-2 input YAMLs for each peptide x {K392, N392}
    python3 << 'BOLTZ_PREP'
import os
from pathlib import Path

fasta_path = Path("results/pepmlm_boltz2_merged.fasta")
boltz_dir = Path("results/boltz2_inputs")
boltz_dir.mkdir(exist_ok=True)

lines = fasta_path.read_text().strip().split("\n")
peptides = []
for i in range(0, len(lines), 2):
    name = lines[i][1:].split("|")[0]
    seq = lines[i + 1]
    peptides.append((name, seq))

# Read ERAP2 full sequence for Boltz-2
erap2 = (
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
erap2_n392 = erap2[:391] + "N" + erap2[392:]

count = 0
for name, pep_seq in peptides:
    for variant, erap_seq in [("k392", erap2), ("n392", erap2_n392)]:
        yaml_name = f"{name}_{variant}.yaml"
        with open(boltz_dir / yaml_name, "w") as f:
            f.write(f"sequences:\n")
            f.write(f"  - protein:\n")
            f.write(f"      id: A\n")
            f.write(f"      sequence: {erap_seq}\n")
            f.write(f"      msa: empty\n")
            f.write(f"  - protein:\n")
            f.write(f"      id: B\n")
            f.write(f"      sequence: {pep_seq}\n")
            f.write(f"      msa: empty\n")
        count += 1

print(f"Created {count} Boltz-2 input YAMLs in {boltz_dir}/")
print(f"Run: boltz predict {boltz_dir}/ --diffusion_samples 3 --seed 42 --output_dir results/boltz2_out/")
BOLTZ_PREP

    echo "Running Boltz-2 (3 diffusion samples)..."
    boltz predict results/boltz2_inputs/ \
        --diffusion_samples 3 \
        --seed 42 \
        --output_dir results/boltz2_out/ \
        2>&1 | tail -20

    echo "Boltz-2 complete. Results in results/boltz2_out/"
else
    echo ""
    echo "=== Boltz-2 not installed ==="
    echo "Download results/ folder and run Boltz-2 separately:"
    echo "  boltz predict results/boltz2_inputs/ --diffusion_samples 3 --seed 42"
fi

echo ""
echo "============================================"
echo "  PepMLM Pipeline Complete — $(date)"
echo "============================================"
echo ""
echo "Output files:"
ls -la results/
echo ""
echo "Next steps:"
echo "  1. Download results/ folder"
echo "  2. Run analyze_pepmlm_results.py locally"
echo "  3. Compare with hand-designed V4 peptides"
