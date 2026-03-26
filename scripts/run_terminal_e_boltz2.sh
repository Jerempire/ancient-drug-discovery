#!/bin/bash
# Terminal E: V2 Capping Redesign — Boltz-2 Selectivity Screen
# Run on Vast.ai RTX 4090 instance
# 9 predictions: 3 variants x 3 targets, 3 diffusion samples each
#
# Instance: 33509299, ssh -p 29298 root@ssh7.vast.ai
# Expected runtime: ~15 min
# Expected cost: ~$0.10

set -e
cd /workspace

YAML_DIR="/workspace/boltz2_yamls"
RESULTS_DIR="/workspace/results_terminal_e"
DIFFUSION_SAMPLES=3
SEED=42

mkdir -p "$RESULTS_DIR"

echo "=== Terminal E: V2 Capping Redesign Boltz-2 Screen ==="
echo "Start: $(date)"
echo ""

# Run all 9 YAMLs
TOTAL=0
PASS=0
for yaml in "$YAML_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    echo "--- Running: $name ---"

    boltz predict "$yaml" \
        --out_dir "$RESULTS_DIR/$name" \
        --diffusion_samples $DIFFUSION_SAMPLES \
        --seed $SEED \
        --accelerator gpu \
        --devices 1 \
        2>&1 | tail -3

    # Extract confidence scores
    for model_json in "$RESULTS_DIR/$name"/boltz_results_*/predictions/*/confidence_*.json; do
        if [ -f "$model_json" ]; then
            echo "  Confidence: $(python3 -c "
import json
with open('$model_json') as f:
    d = json.load(f)
print(f'ipTM={d.get(\"interface_ptm\", d.get(\"iptm\", \"?\")):.3f}, pTM={d.get(\"ptm\", \"?\"):.3f}')
" 2>/dev/null)"
        fi
    done

    TOTAL=$((TOTAL + 1))
    echo ""
done

echo "=== All $TOTAL predictions complete ==="
echo "End: $(date)"

# Compile results summary
python3 -c "
import json, os, glob, sys

results = {}
results_dir = '$RESULTS_DIR'

for yaml_dir in sorted(glob.glob(os.path.join(results_dir, '*'))):
    if not os.path.isdir(yaml_dir):
        continue
    name = os.path.basename(yaml_dir)

    # Find confidence JSONs
    conf_files = glob.glob(os.path.join(yaml_dir, '**/confidence_*.json'), recursive=True)

    scores = []
    for cf in sorted(conf_files):
        with open(cf) as f:
            d = json.load(f)
        scores.append({
            'file': os.path.basename(cf),
            'iptm': d.get('interface_ptm', d.get('iptm', None)),
            'ptm': d.get('ptm', None),
        })

    if scores:
        iptms = [s['iptm'] for s in scores if s['iptm'] is not None]
        avg_iptm = sum(iptms) / len(iptms) if iptms else None
        results[name] = {
            'samples': scores,
            'avg_iptm': avg_iptm,
            'n_samples': len(scores),
        }

# Parse variant and target from name
summary = []
for name, data in sorted(results.items()):
    parts = name.rsplit('_vs_', 1)
    variant = parts[0] if len(parts) == 2 else name
    target = parts[1] if len(parts) == 2 else '?'
    summary.append({
        'variant': variant,
        'target': target,
        'avg_iptm': round(data['avg_iptm'], 3) if data['avg_iptm'] else None,
        'n_samples': data['n_samples'],
        'all_iptm': [round(s['iptm'], 3) for s in data['samples'] if s['iptm']],
    })

# Print table
print()
print('=== SELECTIVITY MATRIX ===')
print(f'{\"Variant\":>20} {\"Target\":>12} {\"avg ipTM\":>10} {\"samples\":>30} {\"Verdict\":>12}')
print('-' * 90)
for row in summary:
    iptm = row['avg_iptm'] or 0
    if 'erap2' in row['target']:
        verdict = 'PASS' if iptm > 0.7 else 'FAIL (weak)'
    else:
        verdict = 'PASS' if iptm < 0.3 else 'FAIL (cross-react!)'
    print(f'{row[\"variant\"]:>20} {row[\"target\"]:>12} {iptm:>10.3f} {str(row[\"all_iptm\"]):>30} {verdict:>12}')

# Save JSON
out = os.path.join(results_dir, 'selectivity_summary.json')
with open(out, 'w') as f:
    json.dump({'results': summary}, f, indent=2)
print(f'\nSaved: {out}')
"

echo ""
echo "=== DONE. Download results with: ==="
echo "scp -P 29298 -r root@ssh7.vast.ai:/workspace/results_terminal_e C:/Users/jmj2z/Projects/medical/ancient-drug-discovery/data/results/v43_validation/terminal_e/"
