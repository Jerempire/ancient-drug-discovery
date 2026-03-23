"""Filter PepMLM hybrid YAMLs to 10+ residue peptides only."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json, shutil
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
hybrid_dir = PROJECT / "data" / "results" / "pepmlm_hybrids" / "boltz_inputs"
upload_dir = PROJECT / "data" / "results" / "pepmlm_hybrids" / "boltz_upload_filtered"
upload_dir.mkdir(parents=True, exist_ok=True)

with open(PROJECT / "data" / "results" / "pepmlm_hybrids" / "hybrid_candidates.json") as f:
    data = json.load(f)

long_seqs = {h["sequence"] for h in data["hybrids"] if h["length"] >= 10}

copied = skipped = 0
for yaml_file in hybrid_dir.glob("*.yaml"):
    content = yaml_file.read_text()
    if any(seq in content for seq in long_seqs):
        shutil.copy2(yaml_file, upload_dir / yaml_file.name)
        copied += 1
    else:
        skipped += 1

print(f"Filtered: {copied} YAMLs kept (10+ aa), {skipped} skipped (8-mers)")
for f in sorted(upload_dir.glob("*.yaml")):
    print(f"  {f.name}")
