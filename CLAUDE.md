# CLAUDE.md - Ancient Drug Discovery Pipeline

## Project Overview
AI-driven pipeline mining evolutionary immune selection for undrugged cancer/autoimmune therapeutic targets. Maps ancient pathogen selection → immune gene variants → cancer immune evasion mechanisms → protein structures → drug candidates.

**Core thesis**: Genes shaped by 10,000 years of plague, malaria, TB, and smallpox are the same genes cancers hijack to evade the immune system. See FRAMEWORK.md for full thesis and evidence.

## Pipeline Stages
1. Data collection (GWAS, ClinVar, UniProt, OpenTargets, DrugBank) — local scripts 01-06
2. Cancer evidence compilation (OpenTargets cancer, ClinicalTrials.gov) — local script 07
3. Variant analysis + ESM-2 scoring — local script 05
4. Protein structure prediction (AlphaFold3/Boltz-2) — Colab
5. Drug design (RFdiffusion, DiffDock, DrugGPT) — Colab
6. Validation (Boltz-2 affinity, RDKit ADMET) — Colab + local
7. Physics-based interface analysis (PyRosetta) — local CPU

## Multi-Target Architecture
All targets defined in `targets.yaml`. Three immune defense steps:

| Step | Gene | Pathogen | Cancer Link |
|---|---|---|---|
| 1. Detection | **ERAP2** (primary) | Black Death | Antigen presentation / immune evasion |
| 1. Detection | **NOD2** | Leprosy | Crohn's → colorectal cancer |
| 2. Activation | **CCR5** | Smallpox | Metastasis (breast, colorectal) |
| 3. Killing | **G6PD** | Malaria | Tumor metabolic vulnerability |
| 3. Killing | **SLC11A1** | TB | Macrophage M1/M2 polarization |
| 3. Killing | **DARC** | Malaria | Prostate cancer disparity |

## Key Files
- `targets.yaml` — multi-target config (variants, compounds, references)
- `FRAMEWORK.md` — unified thesis + historical timeline + evidence
- `NARRATIVE.md` — ERAP2 proof-of-concept narrative
- `scripts/01-05` — ERAP2-specific pipeline (original)
- `scripts/06_multi_target_fetch.py` — generalized data fetcher for all targets
- `scripts/07_cancer_connections.py` — cancer evidence matrix + clinical trials

## Environment
- Conda: `ancient-drug-discovery` (Python 3.11, rdkit, biopython, pandas, pyyaml)
- DrugBank: referenced from `pharmacy-ai-project/data/drugbank/drugbank.db`

### GPU Compute (Vast.ai)
Replaces Colab for GPU inference. Claude Code SSHs directly — no browser needed.

**Two workflows:**

#### Boltz-2 Predictions (primary) — `boltz_runner.py`
Uses custom Docker image (`jerempire/boltz2-ancient:latest`) with pre-baked PyTorch, Boltz-2, CCD, and model weights. Zero setup time.

```bash
# End-to-end (launch → run → download → stop)
python scripts/boltz_runner.py full --yamls path/to/yamls/ --samples 3 --keep

# Or step-by-step:
python scripts/boltz_runner.py warm                    # reuse/restart/launch
python scripts/boltz_runner.py run --yamls path/       # upload + predict
python scripts/boltz_runner.py download                # CIF + JSON (enforced)
python scripts/boltz_runner.py stop                    # pause billing (keep disk)
python scripts/boltz_runner.py start                   # resume later
python scripts/boltz_runner.py destroy                 # full teardown
```

**Instance state:** `.vast_instance_boltz.json` (gitignored). Supports stop/start (disk preserved).
**GPU tiers:** `--tier small` (RTX 4090, <$0.80/hr) for <30 YAMLs; `--tier large` (A100/H100) for big campaigns.
**Docker image source:** `docker/Dockerfile.boltz2`

#### RFdiffusion — `vast_launch.py` (legacy)
**CRITICAL:** Always use Docker image `rosettacommons/rfdiffusion:latest` for RFdiffusion jobs. Do NOT use for Boltz-2.

**Scripts:** `vast_launch.py`, `gpu_setup.sh`, `run_pipeline_gpu.py`, `selectivity_screen.py`
**Instance state:** `.vast_instance.json` (separate from Boltz runner).

**MANDATORY: Always download CIF files** from Boltz-2 predictions (not just JSON scores). CIF files are needed for PyRosetta interface analysis. Save to `data/results/boltz2_complexes/`.

**Cost:** Boltz-2 run (12 complexes x 3 samples): ~$0.10 on RTX 4090. RFdiffusion (40 binders): ~$0.50.
**Results:** `data/results/boltz2_complexes/` (CIF+JSON), `data/results/rfdiffusion/` (PDBs)

## BioReason-Pro (Protein Function Reasoning)

Cloned into `BioReason-Pro/` — multimodal reasoning LLM for protein function prediction (Bo Wang lab, University of Toronto).

**What it does**: Takes a protein sequence → outputs GO term predictions (molecular function, biological process, cellular component), natural language functional reasoning, per-residue attention maps, and binding partner predictions.

**Architecture**: ESM3 protein embeddings + Gene Ontology graph encoder + RL-optimized LLM. Trained on 133K+ proteins across 3,135 organisms. Preferred over curated UniProt entries 79% of the time by human experts.

**Why we use it**: Fills the gap between evolutionary signal detection (scripts 01-07) and structural drug design (Colab/Vast.ai). Validates target mechanisms, surfaces cancer-relevant GO terms, and identifies binding sites for RFdiffusion.

**Integration points**:
| Pipeline Stage | BioReason-Pro Use |
|---|---|
| Script 05 (ESM-2 variant scoring) | Second opinion on variant functional effects via ESM3 embeddings |
| Script 07 (cancer connections) | GO-based functional reasoning may surface mechanisms OpenTargets misses |
| RFdiffusion binder design | Per-residue attention maps → optimal binding site selection |
| Scoring (bayes_target.py) | BioReason-Pro confidence as additional Bayesian signal |

**Key entry points**:
- `BioReason-Pro/gogpt_api.py` — GO term prediction CLI/API
- `BioReason-Pro/interpro_api.py` — InterPro domain annotation
- HuggingFace models: `wanglab/gogpt`, `wanglab/bioreason-pro-sft`, `wanglab/bioreason-pro-rl`

**Wrapper script**: `scripts/08_bioreason_targets.py` — runs all 6 targets through GO-GPT, saves results to `data/processed/<gene>/bioreason/`

**Runtime**: GPU recommended (Vast.ai). Can run GO-GPT on CPU for inference (slower).

## PyRosetta (Physics-Based Interface Analysis)

Post-Boltz-2 validation step using Rosetta energy functions. Runs locally on CPU.

**License**: Academic (free via graylab.jhu.edu). Installed in `ancient-drug-discovery` conda env.

**Script**: `scoring/rosetta_interface_analysis.py` — takes Boltz-2 CIF files, runs:
- **FastRelax** — energy minimization (~90 sec/complex)
- **InterfaceAnalyzer** — real dG_separated, BSA, packstat, shape complementarity, H-bonds
- **Alanine scan** — hotspot residue identification (seconds on CPU vs minutes on GPU)
- **Residue energy breakdown** — per-residue interaction energy for cross-reactivity analysis
- **flex_ddg** — point mutation ΔΔG (V3 K392N variant discrimination)

**Output**: JSON in `data/results/rosetta/`, SQLite `rosetta_metrics` table in `candidates.db`

**When to run**: After every Boltz-2 validation batch, before synthesis decisions.

**Key thresholds**: dG_separated < -15 REU = strong, > -5 = weak, positive = do not synthesize.

## Rules
- Scripts 01-07 run locally (CPU, data collection/querying)
- Notebooks 04-09 run on Colab (GPU inference)
- Never commit data/ contents (large files, .gitignored)
- Cross-validate structure predictions: AlphaFold vs Boltz-2 vs experimental PDB
- Per-target data goes in `data/processed/<gene_symbol>/`
- Patent: one provisional per target class, $75-150 micro entity
