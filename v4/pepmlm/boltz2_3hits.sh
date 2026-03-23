#!/bin/bash
# Boltz-2 validation of 3 PepMLM K392-selective hits
# 3 peptides x 2 variants (K392, N392) = 6 complexes x 3 diffusion samples
set -e

echo "============================================"
echo "  Boltz-2: PepMLM K392-selective hits"
echo "  $(date)"
echo "============================================"

# Install Boltz-2
echo "=== Installing Boltz-2 ==="
pip install -q boltz 2>&1 | tail -3

# Full ERAP2 sequence (Q6P179, 960 aa)
ERAP2_K392="MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRPKDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWLMVNT"

# N392 variant (K→N at position 392, 1-indexed)
ERAP2_N392="${ERAP2_K392:0:391}N${ERAP2_K392:392}"

# 3 PepMLM K392-selective peptides
declare -A PEPTIDES
PEPTIDES[pepmlm_hit1]="KPLALFGLLPYY"
PEPTIDES[pepmlm_hit2]="NVLLPPGLLQLYL"
PEPTIDES[pepmlm_hit3]="KTALLPAGLLGIK"

# Create input YAMLs
mkdir -p /workspace/boltz_inputs
for name in "${!PEPTIDES[@]}"; do
    seq="${PEPTIDES[$name]}"
    for variant in k392 n392; do
        if [ "$variant" == "k392" ]; then
            erap_seq="$ERAP2_K392"
        else
            erap_seq="$ERAP2_N392"
        fi
        cat > /workspace/boltz_inputs/${name}_${variant}.yaml << YAML
sequences:
  - protein:
      id: A
      sequence: ${erap_seq}
      msa: empty
  - protein:
      id: B
      sequence: ${seq}
      msa: empty
YAML
        echo "Created ${name}_${variant}.yaml"
    done
done

echo ""
echo "=== Running Boltz-2 (3 diffusion samples) ==="
echo "6 complexes total"

# Run each individually for progress tracking
for yaml in /workspace/boltz_inputs/*.yaml; do
    name=$(basename "$yaml" .yaml)
    echo ""
    echo "--- $name ---"
    boltz predict "$yaml" \
        --diffusion_samples 3 \
        --seed 42 \
        --out_dir "/workspace/results/${name}" \
        2>&1 | tail -5
done

# Collect scores
echo ""
echo "=== Results ==="
python3 << 'SCORE'
import json, glob, os

results = []
for d in sorted(glob.glob("/workspace/results/pepmlm_*")):
    name = os.path.basename(d)
    # Find scores JSON
    for jf in glob.glob(os.path.join(d, "**/*confidence*scores*.json"), recursive=True) + \
               glob.glob(os.path.join(d, "**/*ranking*.json"), recursive=True):
        try:
            data = json.load(open(jf))
            if isinstance(data, list):
                for entry in data:
                    iptm = entry.get("iptm", entry.get("ipTM", None))
                    ptm = entry.get("ptm", entry.get("pTM", None))
                    if iptm is not None:
                        results.append({"name": name, "iptm": iptm, "ptm": ptm, "file": jf})
            elif isinstance(data, dict):
                iptm = data.get("iptm", data.get("ipTM", None))
                ptm = data.get("ptm", data.get("pTM", None))
                if iptm is not None:
                    results.append({"name": name, "iptm": iptm, "ptm": ptm, "file": jf})
        except:
            pass

if not results:
    # Try alternate score locations
    for d in sorted(glob.glob("/workspace/results/pepmlm_*")):
        for jf in glob.glob(os.path.join(d, "**/*.json"), recursive=True):
            try:
                data = json.load(open(jf))
                print(f"  Found JSON: {jf} keys={list(data.keys()) if isinstance(data, dict) else 'list'}")
            except:
                pass
    print("No ipTM scores found — check output structure manually")
else:
    print(f"{'Name':>30s}  {'ipTM':>8s}  {'pTM':>8s}")
    print("-" * 52)
    for r in sorted(results, key=lambda x: x["name"]):
        iptm = f"{r['iptm']:.4f}" if r['iptm'] is not None else "?"
        ptm = f"{r['ptm']:.4f}" if r['ptm'] is not None else "?"
        print(f"{r['name']:>30s}  {iptm:>8s}  {ptm:>8s}")

    # Save
    with open("/workspace/results/pepmlm_boltz2_scores.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to /workspace/results/pepmlm_boltz2_scores.json")
SCORE

echo ""
echo "============================================"
echo "  Boltz-2 validation complete — $(date)"
echo "============================================"
