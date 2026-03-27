"""
Setup and run Boltz-2 with PDB templates for accurate target structures.

Downloads ERAP2 (5AB0), ERAP1 (2YD0), IRAP (5MJ6) crystal structures
from RCSB PDB and uses them as templates for Boltz-2 predictions.

Run on Vast.ai after installing Boltz-2.
"""
import os
import json
import subprocess
import urllib.request

WORKSPACE = "/workspace"
TEMPLATE_DIR = f"{WORKSPACE}/templates"
YAML_DIR = f"{WORKSPACE}/yamls"
RESULTS_DIR = f"{WORKSPACE}/results"

# PDB IDs for target templates
PDB_TEMPLATES = {
    "erap2": "5AB0",  # ERAP2 closed form, 2.1A resolution
    "erap1": "2YD0",  # ERAP1, 2.7A resolution
    "irap": "5MJ6",   # IRAP closed form, 2.5A resolution
}

# Full target sequences (from UniProt, matching the K392 region we've been using)
ERAP2_SEQ = "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
ERAP1_SEQ = "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR"
IRAP_SEQ = "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN"

# Parent binder (Y87A/Y89A) - the one we've been testing all session
PARENT_BINDER = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

# mpnn06_vagsal binder (our current lead)
MPNN06_BINDER = "IERHYHKSLEEYLKNLPKKVDMLVDLYSKGIFHLDNTNTLVEDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKNANALNR"
MPNN06_VAGSAL = "VAGSAL" + "EAAAA" * 3 + MPNN06_BINDER

targets = {
    "erap2": ERAP2_SEQ,
    "erap1": ERAP1_SEQ,
    "irap": IRAP_SEQ,
}


def download_pdb(pdb_id, out_path):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, out_path)
    print(f"  Saved: {out_path}")


def setup_templates():
    """Download all PDB templates."""
    os.makedirs(TEMPLATE_DIR, exist_ok=True)
    for target, pdb_id in PDB_TEMPLATES.items():
        out_path = f"{TEMPLATE_DIR}/{pdb_id}.pdb"
        if not os.path.exists(out_path):
            download_pdb(pdb_id, out_path)
        else:
            print(f"  {pdb_id}.pdb already exists")
    return {t: f"{TEMPLATE_DIR}/{pid}.pdb" for t, pid in PDB_TEMPLATES.items()}


def build_yaml(name, target_seq, binder_seq, template_path, use_msa=False):
    """Build Boltz-2 YAML with template and optional MSA."""
    # For designed binder: always msa: empty (no homologs)
    # For natural target: use template + optionally MSA
    msa_target = "empty"  # We use templates instead of MSA server

    yaml = f"""version: 1
sequences:
  - protein:
      id: A
      msa: {msa_target}
      sequence: {target_seq}
  - protein:
      id: B
      msa: empty
      sequence: {binder_seq}
templates:
  - path: {template_path}
    ids: [A]
"""
    fname = f"{YAML_DIR}/{name}.yaml"
    with open(fname, "w") as f:
        f.write(yaml)
    return fname


def build_yaml_no_template(name, target_seq, binder_seq):
    """Build Boltz-2 YAML WITHOUT template (control)."""
    yaml = f"""version: 1
sequences:
  - protein:
      id: A
      msa: empty
      sequence: {target_seq}
  - protein:
      id: B
      msa: empty
      sequence: {binder_seq}
"""
    fname = f"{YAML_DIR}/{name}.yaml"
    with open(fname, "w") as f:
        f.write(yaml)
    return fname


def main():
    os.makedirs(YAML_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Step 1: Download templates
    print("=== Downloading PDB templates ===")
    template_paths = setup_templates()

    # Step 2: Build YAMLs
    print("\n=== Building YAMLs ===")

    constructs = {
        "parent": PARENT_BINDER,
        "mpnn06_vagsal": MPNN06_VAGSAL,
    }

    count = 0
    for cname, cseq in constructs.items():
        for tname, tseq in targets.items():
            # With template
            build_yaml(f"{cname}_tmpl_vs_{tname}", tseq, cseq, template_paths[tname])
            count += 1
            print(f"  {cname}_tmpl_vs_{tname} (with {PDB_TEMPLATES[tname]} template)")

            # Without template (control, matches our previous runs)
            build_yaml_no_template(f"{cname}_ctrl_vs_{tname}", tseq, cseq)
            count += 1
            print(f"  {cname}_ctrl_vs_{tname} (no template, control)")

    print(f"\nTotal: {count} YAMLs")
    print(f"  {len(constructs)} constructs x {len(targets)} targets x 2 (template/control)")

    # Step 3: Write run script
    print("\n=== Writing run script ===")
    lines = [
        "#!/bin/bash",
        'echo "=== TEMPLATE vs NO-TEMPLATE SCREEN ==="',
        'echo "Start: $(date)"',
        f"mkdir -p {RESULTS_DIR}",
        "",
    ]

    for yaml_file in sorted(os.listdir(YAML_DIR)):
        if not yaml_file.endswith(".yaml"):
            continue
        name = yaml_file.replace(".yaml", "")
        yaml_path = f"{YAML_DIR}/{yaml_file}"
        out_dir = f"{RESULTS_DIR}/{name}"
        lines.append(f'echo "--- {name} ---"')
        lines.append(f"boltz predict {yaml_path} --out_dir {out_dir} --diffusion_samples 5 --seed 42 --accelerator gpu --devices 1 2>&1 | tail -2")
        lines.append(f'for conf in {out_dir}/boltz_results_*/predictions/*/confidence_*.json; do')
        lines.append(f'  if [ -f "$conf" ]; then')
        lines.append("    python3 -c \"import json,sys;d=json.load(open(sys.argv[1]));print('  ipTM='+str(round(d.get('iptm',0),3)))\" \"$conf\" 2>/dev/null")
        lines.append("  fi")
        lines.append("done")
        lines.append("")

    lines.append('echo "=== DONE ==="')
    lines.append('echo "End: $(date)"')

    run_script = f"{WORKSPACE}/run.sh"
    with open(run_script, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(run_script, 0o755)
    print(f"Run script: {run_script}")
    print(f"Total predictions: {count} x 5 samples = {count * 5} diffusion samples")


if __name__ == "__main__":
    main()
