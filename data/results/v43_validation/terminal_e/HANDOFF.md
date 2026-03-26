# Terminal E Handoff Sheet

**Last updated:** 2026-03-25
**Instance:** 33509299 | SSH: `ssh -p 29298 root@ssh7.vast.ai`
**GPU:** RTX 4090 | Driver 535.183.01

---

## What Terminal E Does

Evaluates whether the V2 protein binder (Y87A_Y89A, 92aa) caps the ERAP2 substrate channel entrance, and tests redesigned variants that improve coverage.

## Current Status

**Phase 1 COMPLETE:** SASA analysis shows binder is BRIDGE_NOT_CAP (30.7% entrance blockage, needs >80%).

**Phase 2 IN PROGRESS:** Running 15 Boltz-2 predictions (5 variants x 3 targets) on Vast.ai.

## Instance Setup Notes

**CRITICAL: CUDA compat fix was applied.** If the instance is restarted or a new one is created, you MUST run:
```bash
# Remove CUDA 12.4 compat library that conflicts with driver 535
rm /etc/ld.so.conf.d/cuda-12-4-compat.conf 2>/dev/null
sed -i '/cuda-12.4\/compat/d' /etc/ld.so.conf.d/*.conf 2>/dev/null
ldconfig

# Verify GPU works
python3 -c "import torch; print(torch.cuda.is_available())"
# Should print: True
```

PyTorch version: 2.5.1+cu121 (installed via `pip install torch==2.5.1+cu121 --index-url https://download.pytorch.org/whl/cu121`)

## Commands to Run

### If predictions are still running:
```bash
# Check progress
ssh -p 29298 root@ssh7.vast.ai "ls /workspace/results_terminal_e/"
```

### If predictions completed — download and analyze:
```bash
cd C:\Users\jmj2z\Projects\medical\ancient-drug-discovery

# Download results
scp -P 29298 -r root@ssh7.vast.ai:/workspace/results_terminal_e "data/results/v43_validation/terminal_e/boltz2_results"

# Destroy instance
vastai destroy instance 33509299

# Run SASA capping re-analysis on best variant's CIF
python scripts/v2_channel_cap_analysis.py
# (update CIF_FILE in the script to point to the best variant's CIF first)
```

### If predictions failed — re-run:
```bash
ssh -p 29298 root@ssh7.vast.ai "rm -rf /workspace/results_terminal_e && bash /workspace/run_terminal_e_boltz2.sh"
```

### If instance died — start fresh:
```bash
# 1. Launch new RTX 4090 instance
vastai search offers "gpu_name=RTX_4090 num_gpus=1 inet_down>100 disk_space>30 reliability>0.95 dph<0.40"
vastai create instance <ID> --image pytorch/pytorch:2.5.1-cuda12.1-cudnn9-devel --disk 30

# 2. Get SSH info
vastai show instances --raw | python -c "import sys,json; d=json.load(sys.stdin); [print(f'ssh -p {i[\"ssh_port\"]} root@{i[\"ssh_host\"]}') for i in d]"

# 3. Install Boltz-2
ssh -p <PORT> root@<HOST> "pip install boltz biopython"

# 4. Upload YAMLs and run script
scp -P <PORT> -r data/boltz2_inputs/*.yaml root@<HOST>:/workspace/boltz2_yamls/
scp -P <PORT> scripts/run_terminal_e_boltz2.sh root@<HOST>:/workspace/  # or recreate from coordination.json notes

# 5. Run predictions
ssh -p <PORT> root@<HOST> "bash /workspace/run_terminal_e_boltz2.sh"
```

## 5 Variants Being Tested

| Variant | Design | Size | Source |
|---------|--------|------|--------|
| northcap | +4aa loop at pos 26 (ENYGS, toward ceiling) | 96aa | Other terminal |
| southcap | +4aa loop at pos 14 (RNSFN, toward floor) | 96aa | Other terminal |
| dualcap | Both insertions | 100aa | Other terminal |
| hybrid_vagsaf | +GGS+VAGSAF C-terminal peptide tail | 101aa | This terminal |
| hybrid_iafsaf | +GGS+IAFSAF C-terminal peptide tail | 101aa | This terminal |

Each tested against: ERAP2-K392, ERAP1, IRAP (3 targets x 3 diffusion samples = 9 CIFs per variant)

## Pass Criteria

- ipTM(ERAP2) > 0.7 for at least one variant
- ipTM(ERAP1) < 0.3 for that same variant
- ipTM(IRAP) < 0.3 for that same variant
- Then re-run SASA analysis: channel entrance SASA reduction must be >80%

## Key Files

| File | Location | Description |
|------|----------|-------------|
| SASA analysis script | `scripts/v2_channel_cap_analysis.py` | Reusable — update CIF_FILE path |
| Phase 1 results | `data/results/v43_validation/terminal_e/v2_channel_cap_analysis.json` | Contact + SASA analysis |
| Redesign strategy | `data/results/v43_validation/terminal_e/CHANNEL_CAP_REDESIGN.md` | Full design rationale |
| Boltz-2 YAMLs (local) | `data/boltz2_inputs/v2_hybrid_*.yaml` | Local copies of hybrid configs |
| Boltz-2 YAMLs (remote) | `/workspace/boltz2_yamls/*.yaml` | All 15 on the instance |
| Run script (remote) | `/workspace/run_terminal_e_boltz2.sh` | Master Boltz-2 runner |
| Coordination | `data/results/v43_validation/coordination.json` | terminal_e section |

## After Results Come In

1. Check selectivity matrix in `/workspace/results_terminal_e/selectivity_summary.json`
2. Find the best variant: highest ERAP2 ipTM with ERAP1 < 0.3 and IRAP < 0.3
3. Download its CIF file
4. Update `CIF_FILE` in `scripts/v2_channel_cap_analysis.py` and re-run SASA analysis
5. If SASA entrance reduction > 80%: SUCCESS — add to MD validation queue
6. Update `data/results/v43_validation/coordination.json` terminal_e status
7. Destroy the Vast.ai instance

## Binder Sequences (for reference)

**Original V2 (92aa):**
`DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN`

**Hybrid VAGSAF (101aa):**
`DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKNGGSVAGSAF`

**Hybrid IAFSAF (101aa):**
`DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKNGGSIAFSAF`

**Northcap (96aa):**
`DIRHYFKSLEEYLKNLPKVVDMLVDENYGSLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN`

**Southcap (96aa):**
`DIRHYFKSLEEYLKRNSFNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN`

**Dualcap (100aa):**
`DIRHYFKSLEEYLKRNSFNLPKVVDMLVDENYGSLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN`
