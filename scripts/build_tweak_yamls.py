"""Build mpnn06 tweak YAMLs and run script for Vast.ai."""
import os

ERAP2 = "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
ERAP1 = "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR"
IRAP = "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN"

targets = {"erap2k392": ERAP2, "erap1": ERAP1, "irap": IRAP}

MPNN06_BINDER = "IERHYHKSLEEYLKNLPKKVDMLVDLYSKGIFHLDNTNTLVEDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKNANALNR"

LINKER3 = "EAAAA" * 3
LINKER2 = "EAAAA" * 2

variants = {
    "mpnn06_5s": ("VAGSAF" + LINKER3 + MPNN06_BINDER, 5),
    "mpnn06_short": ("VAGSAF" + LINKER2 + MPNN06_BINDER, 3),
    "mpnn06_vagsae": ("VAGSAE" + LINKER3 + MPNN06_BINDER, 3),
    "mpnn06_vagsak": ("VAGSAK" + LINKER3 + MPNN06_BINDER, 3),
    "mpnn06_vagsal": ("VAGSAL" + LINKER3 + MPNN06_BINDER, 3),
    "mpnn06_short_vagsae": ("VAGSAE" + LINKER2 + MPNN06_BINDER, 3),
    "mpnn06_short_vagsak": ("VAGSAK" + LINKER2 + MPNN06_BINDER, 3),
}

os.makedirs("/workspace/yamls", exist_ok=True)

count = 0
for vname, (vseq, n_samples) in sorted(variants.items()):
    for tname, tseq in targets.items():
        fname = f"/workspace/yamls/{vname}_vs_{tname}.yaml"
        y = f"version: 1\nsequences:\n  - protein:\n      id: A\n      msa: empty\n      sequence: {tseq}\n  - protein:\n      id: B\n      msa: empty\n      sequence: {vseq}\n"
        with open(fname, "w") as f:
            f.write(y)
        count += 1
        print(f"  {vname}_vs_{tname} ({len(vseq)}aa, {n_samples} samples)")

print(f"\nTotal: {count} YAMLs")

# Write run script
lines = [
    "#!/bin/bash",
    'echo "=== MPNN06 TWEAK SCREEN ==="',
    'echo "Start: $(date)"',
    "mkdir -p /workspace/results",
    "",
]

for vname, (vseq, n_samples) in sorted(variants.items()):
    for tname in targets:
        fname = f"/workspace/yamls/{vname}_vs_{tname}.yaml"
        outdir = f"/workspace/results/{vname}_vs_{tname}"
        lines.append(f'echo "--- {vname}_vs_{tname} ---"')
        lines.append(f"boltz predict {fname} --out_dir {outdir} --diffusion_samples {n_samples} --seed 42 --accelerator gpu --devices 1 2>&1 | tail -2")
        lines.append(f'for conf in {outdir}/boltz_results_*/predictions/*/confidence_*.json; do')
        lines.append(f'  if [ -f "$conf" ]; then')
        lines.append('    python3 -c "import json,sys;d=json.load(open(sys.argv[1]));print(\'  ipTM=\'+str(round(d.get(\'iptm\',0),3)))" "$conf" 2>/dev/null')
        lines.append("  fi")
        lines.append("done")
        lines.append("")

lines.append('echo "=== DONE ==="')
lines.append('echo "End: $(date)"')

with open("/workspace/run.sh", "w") as f:
    f.write("\n".join(lines) + "\n")

os.chmod("/workspace/run.sh", 0o755)
print("Run script written")
