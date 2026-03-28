# Vast.ai Session Plan: New Tools + Cyclic Peptide Docking

## Pre-Session (Local, do before launching instance)

### 1. Prepare RAPiDock virtual screening CSV
Create `data/results/cyclic_peptides/rapidock_screen.csv` with all peptide sequences vs ERAP2/ERAP1/IRAP:

```csv
complex_name,protein_description,peptide_description
relaxed_lasso_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,RCGGWCF
relaxed_lasso_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,RCGGWCF
relaxed_lasso_irap,/workspace/structures/irap_5MJ6_experimental.pdb,RCGGWCF
lasso_bolt_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,RGCWCF
lasso_bolt_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,RGCWCF
lasso_bolt_irap,/workspace/structures/irap_5MJ6_experimental.pdb,RGCWCF
vagsaf_lasso_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,VACSACF
vagsaf_lasso_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,VACSACF
vagsaf_lasso_irap,/workspace/structures/irap_5MJ6_experimental.pdb,VACSACF
hydroxyl_bridge_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,RSGGWF
hydroxyl_bridge_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,RSGGWF
hydroxyl_bridge_irap,/workspace/structures/irap_5MJ6_experimental.pdb,RSGGWF
zn_thiol_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,CGGWAF
zn_thiol_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,CGGWAF
zn_thiol_irap,/workspace/structures/irap_5MJ6_experimental.pdb,CGGWAF
vagsaf_linear_erap2,/workspace/structures/erap2_3SE6_experimental.pdb,VAGSAF
vagsaf_linear_erap1,/workspace/structures/erap1_2YD0_experimental.pdb,VAGSAF
vagsaf_linear_irap,/workspace/structures/irap_5MJ6_experimental.pdb,VAGSAF
```

### 2. Prepare pocket PDBs for RAPiDock (local docking)
RAPiDock works best with pocket-truncated PDBs. Use their pocket tool:
```bash
# On Vast.ai after setup:
cd /workspace/RAPiDock
python pocket_trunction.py --protein /workspace/structures/erap2_3SE6_experimental.pdb --output /workspace/structures/erap2_pocket.pdb
python pocket_trunction.py --protein /workspace/structures/erap1_2YD0_experimental.pdb --output /workspace/structures/erap1_pocket.pdb
python pocket_trunction.py --protein /workspace/structures/irap_5MJ6_experimental.pdb --output /workspace/structures/irap_pocket.pdb
```

---

## Session Order of Operations

### Phase 1: Boltz-2 Cyclic Peptide Docking (~30 min)
Already planned — run the 15 YAMLs that are staged.
```bash
# Install Boltz-2
pip install "boltz<2.2"

# Run all 15 YAMLs (5 peptides x 3 targets)
for yaml in /workspace/cyclic_yamls/*.yaml; do
    boltz predict "$yaml" --diffusion_samples 5 --seed 42 --out_dir /workspace/results/cyclic_boltz2/ --no_kernels
done
```

### Phase 2: RAPiDock Selectivity Screening (~10 min)
Run all cyclic + linear peptides through RAPiDock for fast docking comparison.

```bash
# Install RAPiDock
cd /workspace/RAPiDock
conda create -n RAPiDock python=3.9 pytorch torchvision torchaudio pytorch-cuda=12.4 MDAnalysis pyg pytorch-cluster pytorch-scatter pytorch-sparse pyyaml -c pytorch -c nvidia -c pyg -c conda-forge -y
conda activate RAPiDock
pip install e3nn rdkit-pypi fair-esm

# Generate pocket PDBs
python pocket_trunction.py --protein /workspace/structures/erap2_3SE6_experimental.pdb --output /workspace/structures/erap2_pocket.pdb
python pocket_trunction.py --protein /workspace/structures/erap1_2YD0_experimental.pdb --output /workspace/structures/erap1_pocket.pdb
python pocket_trunction.py --protein /workspace/structures/irap_5MJ6_experimental.pdb --output /workspace/structures/irap_pocket.pdb

# Run virtual screening (all peptides x 3 targets)
python inference.py \
    --config default_inference_args.yaml \
    --protein_peptide_csv /workspace/rapidock_screen.csv \
    --output_dir /workspace/results/rapidock/ \
    --N 10 \
    --model_dir train_models/CGTensorProductEquivariantModel \
    --ckpt rapidock_local.pt \
    --scoring_function ref2015 \
    --batch_size 4 \
    --no_final_step_noise \
    --inference_steps 16 \
    --actual_steps 16 \
    --cpu 10
```

At 0.35 sec/complex x 10 samples x 18 pairs = ~63 seconds total. Compare `ref2015_score.csv` across targets for selectivity.

### Phase 3: PepTune Generation (~20 min)
Generate NEW selective peptide candidates using ERAP2 as the optimization target.

```bash
# Install PepTune deps
cd /workspace/PepTune
pip install torch transformers rdkit-pypi

# Generate with MCTG multi-objective guidance
# Objectives: ERAP2 binding (maximize) + ERAP1 non-binding (minimize cross-reactivity)
python generate_mcts.py \
    --config configs/generate_config.yaml \
    --checkpoint peptune-pretrained.ckpt \
    --num_samples 50 \
    --output_dir /workspace/results/peptune/
```
NOTE: Config will need customization for ERAP2-specific objectives. Check `configs/` and `scoring/` directories for how to plug in custom reward functions.

### Phase 4: RFdiffusion3 Cyclic Binder Design (~30 min)
Design cyclic peptide backbones targeting ERAP2 active site.

```bash
cd /workspace/RFdiffusion3/models/rfd3

# Design cyclic peptide binders for ERAP2
python run.py \
    --input /workspace/structures/erap2_3SE6_experimental.pdb \
    --target_chain A \
    --target_residues 392,398,403,406,412,414 \
    --design_chain B \
    --design_length 8-15 \
    --cyclic_chains B \
    --num_designs 20 \
    --output_dir /workspace/results/rfd3_cyclic/
```
NOTE: Exact CLI flags TBD — check `docs/tutorials/ppi_design_tutorial.md` for syntax. The key params are `cyclic_chains` and targeting the 6 ERAP2-unique channel residues.

---

## Upload Checklist (scp to Vast.ai)

```bash
# Structures
scp -P <PORT> C:/Users/jmj2z/Projects/medical/ancient-drug-discovery/data/structures/*.pdb root@<HOST>:/workspace/structures/

# Cyclic peptide Boltz-2 YAMLs
scp -P <PORT> C:/Users/jmj2z/Projects/medical/ancient-drug-discovery/data/results/cyclic_peptides/boltz2_yamls/*.yaml root@<HOST>:/workspace/cyclic_yamls/

# RAPiDock screening CSV
scp -P <PORT> rapidock_screen.csv root@<HOST>:/workspace/

# Tools (from D:\ staging)
scp -r -P <PORT> D:/ProjectData/ancient-drug-discovery/tools-staging/RAPiDock root@<HOST>:/workspace/
scp -r -P <PORT> D:/ProjectData/ancient-drug-discovery/tools-staging/PepTune root@<HOST>:/workspace/
scp -r -P <PORT> D:/ProjectData/ancient-drug-discovery/tools-staging/RFdiffusion3 root@<HOST>:/workspace/
```

## Download Checklist (after session)
```bash
# All results
scp -r -P <PORT> root@<HOST>:/workspace/results/ D:/ProjectData/ancient-drug-discovery/results/session_newtools/

# Specifically:
# - results/cyclic_boltz2/   → Boltz-2 CIF + ipTM scores
# - results/rapidock/        → ref2015_score.csv + docked PDBs
# - results/peptune/         → generated SMILES + Pareto front
# - results/rfd3_cyclic/     → designed cyclic backbones (PDB)
```

## Expected Costs
- RTX 4090 @ ~$0.30/hr
- Total session: ~90 min = ~$0.45
- Boltz-2 is the bottleneck (15 predictions x 5 samples x ~2 min each = ~2.5 hrs if sequential)
- Run Boltz-2 first, then RAPiDock/PepTune/RFd3 in parallel while waiting

## Decision Matrix After Session
| Result | Action |
|--------|--------|
| RAPiDock confirms Boltz-2 selectivity rankings | Trust both → proceed to synthesis |
| RAPiDock disagrees with Boltz-2 | Run 3rd method (DiffPepDock) as tiebreaker |
| PepTune generates high-scoring novel sequences | Dock with RAPiDock → validate with Boltz-2 |
| RFd3 cyclic designs score > 0.8 ipTM on ERAP2 | Run selectivity screen → add to synthesis panel |
