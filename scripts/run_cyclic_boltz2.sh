#!/bin/bash
# Cyclic Peptide Panel: Boltz-2 Selectivity Screen
# Run on Vast.ai RTX 4090 instance
# 15 predictions: 5 constructs x 3 targets, 3 diffusion samples each
#
# Expected runtime: ~20 min
# Expected cost: ~$0.10

set -e
cd /workspace

YAML_DIR="/workspace/cyclic_boltz2_yamls"
RESULTS_DIR="/workspace/results_cyclic"
DIFFUSION_SAMPLES=3
SEED=42

mkdir -p "$RESULTS_DIR"

echo "=== Cyclic Peptide Panel: Boltz-2 Selectivity Screen ==="
echo "Start: $(date)"
echo "Constructs: relaxed_lasso, vagsaf_lasso, hydroxyl_bridge, lasso_bolt, zn_thiol"
echo "Targets: erap2k392, erap1, irap"
echo "Diffusion samples: $DIFFUSION_SAMPLES, Seed: $SEED"
echo ""

# Install Boltz-2 if not present
if ! command -v boltz &> /dev/null; then
    echo "Installing Boltz-2..."
    pip install boltz 2>&1 | tail -3
fi

# Run all 15 YAMLs
TOTAL=0
for yaml in "$YAML_DIR"/*.yaml; do
    name=$(basename "$yaml" .yaml)
    echo "--- Running: $name ($((TOTAL + 1))/15) ---"

    boltz predict "$yaml" \
        --out_dir "$RESULTS_DIR/$name" \
        --diffusion_samples $DIFFUSION_SAMPLES \
        --seed $SEED \
        --accelerator gpu \
        --devices 1 \
        2>&1 | tail -5

    # Extract confidence scores
    for model_json in "$RESULTS_DIR/$name"/boltz_results_*/predictions/*/confidence_*.json; do
        if [ -f "$model_json" ]; then
            echo "  $(python3 -c "
import json
with open('$model_json') as f:
    d = json.load(f)
iptm = d.get('interface_ptm', d.get('iptm', '?'))
ptm = d.get('ptm', '?')
print(f'ipTM={iptm:.3f}, pTM={ptm:.3f}')
" 2>/dev/null)"
        fi
    done

    TOTAL=$((TOTAL + 1))
    echo ""
done

echo "=== All $TOTAL predictions complete ==="
echo "End: $(date)"

# Compile selectivity matrix
python3 -c "
import json, os, glob

results = {}
results_dir = '$RESULTS_DIR'

for yaml_dir in sorted(glob.glob(os.path.join(results_dir, '*'))):
    if not os.path.isdir(yaml_dir):
        continue
    name = os.path.basename(yaml_dir)

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

# Parse construct and target from name
summary = []
for name, data in sorted(results.items()):
    parts = name.rsplit('_vs_', 1)
    construct = parts[0] if len(parts) == 2 else name
    target = parts[1] if len(parts) == 2 else '?'
    summary.append({
        'construct': construct,
        'target': target,
        'avg_iptm': round(data['avg_iptm'], 3) if data['avg_iptm'] else None,
        'n_samples': data['n_samples'],
        'all_iptm': [round(s['iptm'], 3) for s in data['samples'] if s['iptm']],
    })

# Print selectivity matrix
print()
print('=== CYCLIC PEPTIDE SELECTIVITY MATRIX ===')
print(f'{\"Construct\":>20} {\"Target\":>12} {\"avg ipTM\":>10} {\"samples\":>30} {\"Verdict\":>15}')
print('-' * 95)
for row in summary:
    iptm = row['avg_iptm'] or 0
    if 'erap2' in row['target']:
        verdict = 'PASS (binds)' if iptm > 0.7 else ('WEAK' if iptm > 0.5 else 'FAIL (no bind)')
    else:
        verdict = 'PASS (evades)' if iptm < 0.3 else ('CONCERN' if iptm < 0.5 else 'FAIL (cross-react!)')
    print(f'{row[\"construct\"]:>20} {row[\"target\"]:>12} {iptm:>10.3f} {str(row[\"all_iptm\"]):>30} {verdict:>15}')

# Construct-level summary
print()
print('=== PER-CONSTRUCT SUMMARY ===')
constructs = {}
for row in summary:
    c = row['construct']
    if c not in constructs:
        constructs[c] = {}
    constructs[c][row['target']] = row['avg_iptm'] or 0

for c, targets in sorted(constructs.items()):
    e2 = targets.get('erap2k392', 0)
    e1 = targets.get('erap1', 0)
    ir = targets.get('irap', 0)
    delta_e1 = e2 - e1
    delta_ir = e2 - ir
    selective = e2 > 0.5 and e1 < 0.3 and ir < 0.3
    print(f'{c:>20}: ERAP2={e2:.3f}  ERAP1={e1:.3f}  IRAP={ir:.3f}  '
          f'delta_E1={delta_e1:+.3f}  delta_IR={delta_ir:+.3f}  '
          f'{\"SELECTIVE\" if selective else \"NOT_SELECTIVE\"}')

# Rank by selectivity
ranked = sorted(constructs.items(),
                key=lambda x: x[1].get('erap2k392', 0) - max(x[1].get('erap1', 0), x[1].get('irap', 0)),
                reverse=True)
print()
print('=== RANKING (by selectivity delta) ===')
for rank, (c, targets) in enumerate(ranked, 1):
    e2 = targets.get('erap2k392', 0)
    best_off = max(targets.get('erap1', 0), targets.get('irap', 0))
    print(f'  {rank}. {c}: delta={e2 - best_off:+.3f} (ERAP2={e2:.3f}, best_off={best_off:.3f})')

# Save JSON
out = os.path.join(results_dir, 'selectivity_summary.json')
with open(out, 'w') as f:
    json.dump({'results': summary, 'constructs': {c: t for c, t in constructs.items()}}, f, indent=2)
print(f'\nSaved: {out}')
"

echo ""
echo "=== CHECKPOINT: Results saved to $RESULTS_DIR ==="
echo "Download with:"
echo "  scp -P <PORT> -r root@<HOST>:/workspace/results_cyclic ."
