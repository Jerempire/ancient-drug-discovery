"""
BindCraft channel-capping binder design for ERAP2.

Designs 20-50 residue protein binders targeting the ERAP2 channel mouth
surface using the 12 ERAP2-unique residues identified by structural alignment.

Run on Vast.ai GPU instance after BindCraft installation.

Target: ERAP2 channel mouth (surface capping, NOT channel interior)
Selectivity: Contacts 12 residues unique to ERAP2 vs ERAP1/IRAP
"""
import os
import json
import sys

# ============================================================
# ERAP2-unique channel mouth residues (from BLOSUM62 alignment)
# ============================================================
# These residues are surface-exposed, within 20A of the active site zinc,
# and differ from BOTH ERAP1 and IRAP at the aligned position.

TARGET_RESIDUES = {
    # resid: (ERAP2_aa, ERAP1_equiv, IRAP_equiv, notes)
    215: ('L', 'A198', 'A310', 'Size difference - Leu vs Ala'),
    308: ('D', 'E291', 'E402', 'Charge/shape shift'),
    317: ('Y', 'P300', 'Q411', 'AROMATIC vs Pro vs polar - HIGH VALUE'),
    347: ('T', 'S330', 'E441', 'Polar vs charged'),
    363: ('W', 'G346', 'L457', 'TRYPTOPHAN vs Gly - HIGHEST VALUE'),
    414: ('D', 'G397', 'Y508', 'Charged vs tiny vs aromatic'),
    423: ('E', 'D406', 'K517', 'Charge flip on IRAP'),
    438: ('K', 'T421', 'S532', 'Charged vs polar'),
    488: ('R', 'K471', 'A582', 'Size/charge change'),
    857: ('A', 'Q834', 'F926', 'Small vs large'),
    888: ('D', 'E865', 'P957', 'Charge/shape shift'),
    896: ('M', 'H873', 'N965', 'Hydrophobic vs imidazole vs polar'),
}

# Priority residues for core contacts (largest shape differences)
CORE_RESIDUES = [363, 317, 414, 438, 488]  # W, Y, D, K, R

# All target residues
ALL_RESIDUES = sorted(TARGET_RESIDUES.keys())

def write_bindcraft_target_settings(output_dir):
    """Write BindCraft target settings JSON."""
    os.makedirs(output_dir, exist_ok=True)

    settings = {
        "target_pdb": "/workspace/structures/erap2_3SE6_experimental.pdb",
        "target_chain": "A",
        "target_residues": ALL_RESIDUES,
        "design_chain": "B",
        "design_length_min": 25,
        "design_length_max": 50,
        "num_designs": 20,
        "hotspot_residues": CORE_RESIDUES,
        "notes": "ERAP2 channel mouth capping binder. 12 ERAP2-unique surface residues."
    }

    settings_path = os.path.join(output_dir, "erap2_cap_target.json")
    with open(settings_path, 'w') as f:
        json.dump(settings, f, indent=2)
    print(f"Target settings written to {settings_path}")
    return settings_path


def write_vastai_setup_script(output_dir):
    """Write the Vast.ai setup + run script."""
    os.makedirs(output_dir, exist_ok=True)

    script = """#!/bin/bash
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

python bindcraft.py \\
    --settings settings_target/erap2_channel_cap.json \\
    --filters settings_filters/default_filters.json \\
    --advanced settings_advanced/default_advanced.json \\
    --output /workspace/results/bindcraft/ \\
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
    boltz predict "$yaml" \\
        --diffusion_samples 3 \\
        --seed 42 \\
        --out_dir /workspace/results/bindcraft_boltz2/ \\
        --no_kernels \\
        2>&1 | tail -1
done

echo ""
echo "============================================"
echo "ALL DONE. Download results with:"
echo "  scp -r -P PORT root@HOST:/workspace/results/ ."
echo "============================================"
"""

    script_path = os.path.join(output_dir, "run_bindcraft_erap2.sh")
    with open(script_path, 'w', newline='\n') as f:
        f.write(script)
    print(f"Setup script written to {script_path}")
    return script_path


if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "scripts")
    write_bindcraft_target_settings(output_dir)
    write_vastai_setup_script(output_dir)

    print("\n" + "="*50)
    print("NEXT STEPS:")
    print("="*50)
    print("1. Launch Vast.ai: python scripts/vast_launch.py launch")
    print("2. Upload structures:")
    print("   scp -P PORT data/structures/*.pdb root@HOST:/workspace/structures/")
    print("3. Upload BindCraft repo:")
    print("   scp -r -P PORT D:/ProjectData/ancient-drug-discovery/tools-staging/BindCraft root@HOST:/workspace/")
    print("4. Upload + run setup script:")
    print("   scp -P PORT scripts/run_bindcraft_erap2.sh root@HOST:/workspace/")
    print("   ssh -p PORT root@HOST 'bash /workspace/run_bindcraft_erap2.sh'")
    print("5. Download results:")
    print("   python scripts/vast_launch.py download")
