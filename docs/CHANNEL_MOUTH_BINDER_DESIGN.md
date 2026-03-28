# ERAP2 Channel Mouth Binder Design

## Date: 2026-03-28
## Status: Backbone validated, sequence optimization pending

## Problem
Linear peptides (VAGSAF) achieve triple selectivity but eject from ERAP2 channel in MD.
Cyclic peptides stay in pocket but lose selectivity (bind ERAP1/IRAP too).
Core tension: channel interior is conserved across M1 aminopeptidases.

## Solution: Channel-Mouth Capping
Instead of putting peptides INSIDE the channel, design a 25-50 residue protein binder
that sits ON TOP of the channel entrance. Selectivity comes from surface contacts
(divergent between enzymes), not channel interior contacts (conserved).

## ERAP2-Unique Channel Mouth Residues (12 total)
Identified via BLOSUM62 structural alignment of ERAP2 (3SE6) vs ERAP1 (2YD0) vs IRAP (5MJ6).
All surface-exposed (SASA > 40 A²), within 20A of active site zinc.

| ERAP2 | ERAP1 equiv | IRAP equiv | Notes |
|-------|-------------|------------|-------|
| **L215** | A198 | A310 | Size difference (Leu vs Ala) |
| **D308** | E291 | E402 | Charge/shape shift |
| **Y317** | P300 | Q411 | Aromatic vs Pro vs polar — HIGH VALUE |
| **T347** | S330 | E441 | Polar vs charged |
| **W363** | G346 | L457 | Trp vs Gly — HIGHEST VALUE |
| **D414** | G397 | Y508 | Charged vs tiny vs aromatic |
| **E423** | D406 | K517 | Charge flip on IRAP |
| **K438** | T421 | S532 | Charged vs polar |
| **R488** | K471 | A582 | Size/charge change |
| **A857** | Q834 | F926 | Small vs large |
| **D888** | E865 | P957 | Charge/shape shift |
| **M896** | H873 | N965 | Hydrophobic vs imidazole vs polar |

**Priority hotspots:** W363, Y317 (large aromatic side chains completely absent in ERAP1/IRAP)

**Target residue list:** `215,308,317,347,363,414,423,438,488,857,888,896`

## Sequence Identity
- ERAP2 vs ERAP1: 54.6% (869 vs 858 aa)
- ERAP2 vs IRAP: 46.2% (869 vs 870 aa)

## Interior Channel Residues (for reference — NOT target)
| ERAP2 | ERAP1 | IRAP | Note |
|-------|-------|------|------|
| K/N392 | N375 | N486 | Selectivity handle (K392 = Black Death variant) |
| Y398 | F381 | F492 | Unique |
| A403 | S386 | S497 | Unique |
| A406 | V389 | K500 | Unique |
| Q412 | K395 | S506 | Unique |
| D414 | G397 | Y508 | Unique (also surface-exposed = in both lists) |

## BindCraft Results (2026-03-28)

### Setup
- Instance: RTX 4090 ($0.315/hr), PyTorch Docker
- Target: 345-residue pocket PDB (25A around channel mouth)
- Full ERAP2 (869 res) OOMs at 24GB VRAM — pocket required

### Trajectory Results
| Design | Length | pLDDT | i_pTM | dG (REU) | MPNN Pass? |
|--------|--------|-------|-------|----------|------------|
| l50_s893817 | 50 | 0.77 | 0.86 | -73.4 | No (0/16) |
| l25_s27 | 25 | 0.85 | 0.92 | scoring OK | No (0/1) |

### Key Finding
**Backbones are excellent** — i_pTM 0.86-0.92, dG -73.4 REU.
The channel-mouth approach works geometrically.
**ProteinMPNN sequence redesign is the bottleneck** — sequences fail AF2 re-validation.
Even with relaxed filters (pLDDT 0.80→0.65, i_pTM 0.50→0.40), no designs passed.

### Root Cause
The target surface is flat and scattered (12 residues across 20A) — hard for MPNN
to find sequences maintaining both binder fold AND all interface contacts simultaneously.
BindCraft's two-step approach (hallucinate backbone → redesign sequence) fails here.

### Lesson
Need joint backbone+sequence co-design (Proteina-Complexa) instead of
sequential backbone-then-sequence (BindCraft).

## Next Steps

### 1. Proteina-Complexa (Docker approach)
- Build Docker image from `env/docker/Dockerfile` in the repo
- Requires Python 3.12+ and `atomworks` (NVIDIA internal framework)
- Push to Docker Hub, launch Vast.ai with that image
- Joint sequence+backbone design avoids MPNN bottleneck
- Target: same 12 channel-mouth residues, 25-50 aa binders

### 2. Boltz-2 Selectivity Screen (after getting designs)
- Run each design vs ERAP2 K392, ERAP1 (2YD0), IRAP (5MJ6)
- Use `boltz_runner.py` with PDB templates
- `--diffusion_samples 5` for selectivity claims

### 3. Fallback: BindCraft overnight run
- If Proteina-Complexa also fails, BindCraft with relaxed filters
  running overnight (8-10 hrs) would likely produce designs
- Cost: ~$2.50-3.00

## Files
- Pocket PDB: `D:/ProjectData/ancient-drug-discovery/results/bindcraft_session_20260328/erap2_pocket.pdb`
- BindCraft logs: `D:/ProjectData/ancient-drug-discovery/results/bindcraft_session_20260328/results/`
- Relaxed filters: `scripts/relaxed_filters.json`
- Analysis scripts (temp, deleted): `tmp_channel_analysis.py`, `tmp_proper_alignment.py`
- RAPiDock screening CSV: `data/results/cyclic_peptides/rapidock_screen.csv`

## Compute Cost
- Vast.ai RTX 4090: ~4 hours at $0.315/hr = **~$1.26 total**
