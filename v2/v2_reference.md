# V2 Protein Binder Synthesis Candidates

**Status**: READY TO ORDER (pending patent filing)
**Last Updated**: 2026-03-23
**Estimated cost**: ~$2,500-4,000 for 4 constructs
**Origin**: Proteina-Complexa (NVIDIA) n248 family, targeting ERAP2 divergent channel (res 353-412)

---

## Constructs (REVISED 2026-03-22)

All scores are 3-sample Boltz-2 averages (seed 42).

### 1. n248_trim_c5_Y87A_Y89A — PRIMARY LEAD

- **Length**: 92 aa
- **Sequence**: `DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN`
- **Mutations**: Y87A + Y89A (C-terminal IRAP-only interface residues, identified by PyRosetta alanine scan)
- **ERAP2 ipTM**: 0.748 | **ERAP1**: 0.112 | **IRAP**: 0.186 | **ANPEP**: 0.183
- **Rationale**: Cleanest selectivity panel of any construct. IRAP dropped from 0.601 (parent) to 0.186. ERAP2 actually improved (+0.057 vs parent).

### 2. n248_trim_c5 — Parent Comparator

- **Length**: 92 aa
- **Sequence**: `DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN`
- **ERAP2 ipTM**: 0.691 | **ERAP1**: 0.362 | **IRAP**: 0.601 | **ANPEP**: 0.277
- **Rationale**: Shows Y87A/Y89A improvement. If lead fails in SPR, this tells you whether the mutation or the whole family is the issue.

### 3. n248_wt — Full-Length Comparator

- **Length**: 97 aa
- **Sequence**: `DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK`
- **ERAP2 ipTM**: 0.475 | **ERAP1**: 0.182 | **IRAP**: 0.491 | **ANPEP**: 0.310
- **Rationale**: Anchor reference. If trim or Y87A/Y89A behave unexpectedly, WT tells you whether trimming/mutations caused it.

### 4. n248_ko_all_aromatics — Negative Control

- **Length**: 97 aa
- **Sequence**: `DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK`
- **Mutations**: All 16 F/Y/W → A
- **ERAP2 ipTM**: 0.315 (binding killed)
- **Rationale**: Mechanism validation. If this doesn't bind in SPR, proves binding is aromatic-network-dependent.

### DROPPED: n248_trim_c5_Y4A

- Multi-sample re-validation revealed ANPEP cross-reactivity at 0.812 (was 0.188 in single-sample). Dealbreaker.

---

## Structure Files

| Construct | ESMFold Monomer | Boltz-2 Complex CIFs | Rosetta Analysis |
|---|---|---|---|
| Y87A_Y89A | `data/results/esmfold/n248_trim_c5_Y87A_Y89A.pdb` | `data/results/y87a_cif_files/n248_trim_c5_Y87A_Y89A_*.cif` | — |
| trim_c5 | `data/results/esmfold/n248_trim_c5.pdb` | `data/results/y87a_cif_files/n248_trim_c5_*.cif` | `data/results/rosetta_pipeline/rosetta_n248_trim_c5_*.json` |
| wt | `data/results/esmfold/n248_wt.pdb` | `data/results/y87a_cif_files/n248_wt_*.cif` | `data/results/rosetta_pipeline/rosetta_n248_wt_*.json` |
| ko_all_aromatics | `data/results/esmfold/n248_ko_all_aromatics.pdb` | — | `data/results/rosetta_pipeline/rosetta_n248_ko_all_aromatics_*.json` |

Additional Boltz-2 full predictions (3 diffusion samples each): `data/results/boltz2_complexes/boltz2/`

---

## SPR Panel (in order)

| Target | Protein | Source | Purpose |
|---|---|---|---|
| ERAP2 | Human recombinant | Commercial (R&D Systems, Sino Bio) | Primary binding |
| IRAP/LNPEP | Human recombinant | Commercial (R&D Systems) | Safety — IRAP cross-reactivity confirmed by PyRosetta |
| ERAP1 | Human recombinant | Commercial (R&D Systems, Sino Bio) | Selectivity confirmation |

## Go/No-Go Gates

| Gate | Pass | Caution | Fail |
|---|---|---|---|
| ERAP2 KD | < 1 uM | 1-10 uM | > 10 uM |
| ERAP1 selectivity | > 10x | 5-10x | < 5x |
| IRAP selectivity | > 10x | 5-10x | < 5x |
| Negative control | No binding | Weak binding | Same as lead |

## Pre-Order Checklist

- [ ] File provisional patent ($75 micro entity)
- [ ] Confirm synthesis vendor (GenScript, Twist, or similar)
- [ ] Confirm recombinant expression system (E. coli for 92-97aa range)
- [ ] Order ERAP2/ERAP1/IRAP recombinant proteins for SPR
- [ ] Book SPR instrument time

## Key Validation Notes

- PyRosetta: all interfaces STRONG (dG -88 to -120 REU, BSA 3150-4100 A²)
- PyRosetta: IRAP cross-reactivity is real — physics says IRAP binding equals/exceeds ERAP2 for parent. Y87A/Y89A fixes this.
- PyRosetta: ko_all_aromatics still binds at dG -95.9 REU — backbone contributes, not just aromatics
- ESMFold monomer: leads 60-64 pLDDT (partial fold, moderate risk). ko_all_aromatics 44.1 (disordered).
- All sequences are novel — no natural protein has >30% identity
- 92-97aa is within E. coli expression range
