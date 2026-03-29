#!/bin/bash
set -e

echo "============================================"
echo "BindCraft ERAP2 Channel-Cap Setup"
echo "============================================"

# ---- Install BindCraft ----
cd /workspace
if [ ! -d "BindCraft" ]; then
    echo "Cloning BindCraft..."
    git clone https://github.com/martinpacesa/BindCraft.git
    cd BindCraft
    bash install_bindcraft.sh
else
    echo "BindCraft already installed."
    cd BindCraft
fi

# ---- Prepare target ----
echo ""
echo "Setting up ERAP2 target..."
mkdir -p /workspace/structures /workspace/results/bindcraft

# Check structures exist
if [ ! -f /workspace/structures/erap2_3SE6_experimental.pdb ]; then
    echo "ERROR: ERAP2 structure not found at /workspace/structures/erap2_3SE6_experimental.pdb"
    echo "Upload with: scp -P PORT structures/*.pdb root@HOST:/workspace/structures/"
    exit 1
fi

# ---- Write BindCraft target config ----
cat > /workspace/settings_target/erap2_channel_cap.json << 'TARGETEOF'
{
    "target_pdb_path": "/workspace/structures/erap2_3SE6_experimental.pdb",
    "chains_to_design": "A",
    "target_hotspot_residues": ["A363", "A317", "A414", "A438", "A488", "A215", "A308", "A347", "A423", "A857", "A888", "A896"],
    "lengths": [25, 30, 35, 40, 45, 50],
    "number_of_final_designs": 20
}
TARGETEOF

echo "Target config written."

# ---- Run BindCraft ----
echo ""
echo "============================================"
echo "Running BindCraft (20 designs, 25-50 aa)"
echo "Target: ERAP2 channel mouth (12 unique residues)"
echo "Core hotspots: W363, Y317, D414, K438, R488"
echo "============================================"

python bindcraft.py \
    --settings settings_target/erap2_channel_cap.json \
    --filters settings_filters/default_filters.json \
    --advanced settings_advanced/default_advanced.json \
    --output /workspace/results/bindcraft/ \
    2>&1 | tee /workspace/results/bindcraft/bindcraft_run.log

echo ""
echo "BindCraft complete. Results in /workspace/results/bindcraft/"
echo "============================================"

# ---- Validate with Boltz-2 ----
echo ""
echo "Installing Boltz-2 for validation..."
pip install "boltz<2.2" -q

echo "Running Boltz-2 validation of BindCraft designs..."
mkdir -p /workspace/results/bindcraft_boltz2

# Generate Boltz-2 YAMLs for each design
for pdb in /workspace/results/bindcraft/design_*.pdb; do
    if [ -f "$pdb" ]; then
        name=$(basename "$pdb" .pdb)

        # Extract designed chain sequence
        python3 -c "
from Bio.PDB import PDBParser
import warnings
warnings.filterwarnings('ignore')
p = PDBParser(QUIET=True)
s = p.get_structure('s', '$pdb')
aa = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
for chain in s[0].get_chains():
    seq = ''.join(aa.get(r.get_resname(),'X') for r in chain.get_residues() if r.id[0]==' ')
    print(f'{chain.id}:{seq}')
" > /tmp/${name}_seqs.txt

        # Get chain B (designed binder) sequence
        binder_seq=$(grep "^B:" /tmp/${name}_seqs.txt | cut -d: -f2)

        if [ -n "$binder_seq" ]; then
            # ERAP2 channel region sequence (from 3SE6 chain A, residues near active site)
            # Use the full ERAP2 sequence from the PDB
            erap2_seq=$(grep "^A:" /tmp/${name}_seqs.txt | cut -d: -f2)

            # Write Boltz-2 YAML for ERAP2
            cat > /workspace/results/bindcraft_boltz2/${name}_vs_erap2.yaml << YAMLEOF
version: 1
sequences:
- protein:
    id: A
    sequence: ${erap2_seq}
    msa: empty
    templates:
    - pdb: /workspace/structures/erap2_3SE6_experimental.pdb
      chain_id: [A]
- protein:
    id: B
    sequence: ${binder_seq}
    msa: empty
YAMLEOF

            # Write Boltz-2 YAML for ERAP1 (counter-screen)
            cat > /workspace/results/bindcraft_boltz2/${name}_vs_erap1.yaml << YAMLEOF
version: 1
sequences:
- protein:
    id: A
    sequence: placeholder_erap1
    msa: empty
    templates:
    - pdb: /workspace/structures/erap1_2YD0_experimental.pdb
      chain_id: [A]
- protein:
    id: B
    sequence: ${binder_seq}
    msa: empty
YAMLEOF

            # Write Boltz-2 YAML for IRAP (counter-screen)
            cat > /workspace/results/bindcraft_boltz2/${name}_vs_irap.yaml << YAMLEOF
version: 1
sequences:
- protein:
    id: A
    sequence: placeholder_irap
    msa: empty
    templates:
    - pdb: /workspace/structures/irap_5MJ6_experimental.pdb
      chain_id: [A]
- protein:
    id: B
    sequence: ${binder_seq}
    msa: empty
YAMLEOF

            echo "  Generated YAMLs for $name (binder: ${#binder_seq} aa)"
        fi
    fi
done

echo ""
echo "Running Boltz-2 selectivity screen..."
for yaml in /workspace/results/bindcraft_boltz2/*.yaml; do
    echo "  Predicting: $(basename $yaml)"
    boltz predict "$yaml" \
        --diffusion_samples 3 \
        --seed 42 \
        --out_dir /workspace/results/bindcraft_boltz2/ \
        --no_kernels \
        2>&1 | tail -1
done

echo ""
echo "============================================"
echo "ALL DONE. Download results with:"
echo "  scp -r -P PORT root@HOST:/workspace/results/ ."
echo "============================================"
