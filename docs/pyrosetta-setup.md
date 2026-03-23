# PyRosetta Setup on Vast.ai

PyRosetta runs on the remote Vast.ai Linux instance, not locally on Windows.

## Architecture

```
Your Windows PC (SSH client)
    | ssh -p <PORT> root@<HOST>
Remote Vast.ai RTX 3090 (Linux)
    |-- Boltz-2 (GPU) -> predicts complexes, outputs CIF files
    |-- PyRosetta (CPU) -> reads CIFs, runs interface analysis
    | scp download
Results back on your Windows PC
```

## Instance Setup

**Image**: `pytorch/pytorch:2.5.1-cuda12.4-cudnn9-devel`

PyRosetta installs via the open conda channel — no credentials or license key needed.

```bash
# Install PyRosetta + dependencies
conda install -y -c https://conda.rosettacommons.org pyrosetta
pip install gemmi biopython

# Verify
python3 -c "import pyrosetta; pyrosetta.init(); print('PyRosetta: OK')"
python3 -c "import boltz; print('Boltz-2: OK')"
python3 -c "import gemmi; print('gemmi: OK')"
python3 -c "from Bio import PDB; print('BioPython: OK')"
```

## Running the Pipeline

The combined pipeline runs Boltz-2 first (GPU, generates CIF files), then PyRosetta (CPU, analyzes interfaces).

```bash
# From the Vast.ai instance
python3 /workspace/scripts/run_rosetta_pipeline.py
```

Expects ~40-50 minutes for 6 complexes (Boltz-2 predictions + PyRosetta analysis).

## Downloading Results

```bash
# From your local Windows machine
scp -P <PORT> root@<HOST>:/workspace/results/rosetta/*.json data/results/rosetta/
scp -P <PORT> root@<HOST>:/workspace/results/boltz2_complexes/*.cif data/results/boltz2_complexes/
```

## Key Thresholds

| Metric | Weak | Adequate | Strong |
|---|---|---|---|
| dG_separated (REU) | > -5 | < -15 | < -25 |
| dSASA_int (A^2) | < 500 | > 800 | > 1200 |
| packstat | < 0.50 | > 0.65 | > 0.75 |
| shape_complementarity | < 0.55 | > 0.65 | > 0.75 |

**Kill gate**: dG_separated > -5 REU = interface is not real, do not synthesize.

## Notes

- PyRosetta never needs to run on Windows — all analysis happens on the rented Linux box
- Only JSON results and CIF files are downloaded locally
- The conda channel (`conda.rosettacommons.org`) provides native Linux builds freely
- Verified working: 2026-03-22, RTX 3090 (25.3GB VRAM), PyTorch 2.5.1+cu124
