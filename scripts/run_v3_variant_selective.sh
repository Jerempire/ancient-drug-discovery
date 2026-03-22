#!/bin/bash
# run_v3_variant_selective.sh — V3 campaign: K392/N392 variant-selective binders
#
# Runs Proteina-Complexa with hotspots centered on 392/393 region,
# then counter-screens top candidates against both K392 and N392 variants.
#
# Usage (on Vast.ai A6000):
#   bash /workspace/scripts/run_v3_variant_selective.sh

set -e
cd /workspace/Proteina-Complexa
source .venv/bin/activate
source env.sh 2>/dev/null || true
export TERM=xterm

AF2_PATH="/workspace/Proteina-Complexa/community_models/ckpts/AF2"

echo "============================================"
echo "  V3: K392/N392 VARIANT-SELECTIVE CAMPAIGN"
echo "  $(date)"
echo "============================================"

# Generate across 3 length tiers with 392/393 hotspots
for TIER in short medium long; do
    echo ""
    echo "=== TIER: $TIER ==="
    echo "Start: $(date)"

    complexa generate configs/search_binder_local_pipeline.yaml \
        ++run_name=v3_${TIER} \
        ++generation.task_name=erap2_v3_${TIER} \
        ++gen_njobs=1 \
        ++generation.reward_model.reward_models.af2folding.af_params_dir=$AF2_PATH \
        ++generation.dataloader.dataset.nres.nsamples=20 \
        ++generation.search.best_of_n.replicas=4 \
        2>&1 | tail -20

    echo "Tier $TIER done: $(date)"
done

echo ""
echo "============================================"
echo "  V3 GENERATION COMPLETE: $(date)"
echo "============================================"
