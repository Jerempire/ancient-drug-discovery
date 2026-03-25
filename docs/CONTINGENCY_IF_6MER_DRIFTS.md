# Contingency Plan: If 6-mer Peptides Drift in MD

**Trigger:** IAFSAF (or all Elite Four) show RMSD > 5 Å / COM drift > 5 Å in ERAP2-K392 MD
**Date prepared:** 2026-03-25

---

## Option 1: Fix the Method — Restrained Equilibration (TRY FIRST)

**Cost:** ~$0.50 (re-run 1 simulation) | **Time:** ~3 hours | **Probability of fixing it:** High

The current protocol throws the peptide into a fully solvated box with no initial restraints. A 6-mer is small enough that water molecules can jostle it out before it settles.

**Fix:** Add 1-2 ns of restrained equilibration:
- Hold peptide heavy atoms in place (10 kJ/mol/nm² restraints)
- Let water + protein relax around the peptide
- Gradually release restraints over 500 ps
- THEN run unrestrained production MD

This is standard practice for small ligand MD. We skipped it for speed. If the peptide stays after proper equilibration, the binding is real — we just needed to let the system settle first.

**How to implement:** Add to `openmm_md_protocol.py` between NPT equilibration and production:
```python
# Restrained equilibration (1 ns)
# Add position restraints on peptide heavy atoms
# Gradually reduce: 10 → 5 → 2 → 0 kJ/mol/nm² over 4 × 250 ps
```

## Option 2: Go Back to 10-11 mers

**Cost:** ~$2 (4 MD runs) | **Time:** ~8 hours | **Probability of working:** High

KILKLYSSKKY (11-mer) scored 0.892 ipTM with +0.427 selectivity and the KKY motif is confirmed load-bearing (81% selectivity drop without it). Longer peptides have more surface area to grip the channel walls.

**Trade-off:** Lose ERAP1 molecular ruler evasion (11-mers ARE long enough for ERAP1 to close around). But gain binding stability. ERAP1 selectivity would need to come from sequence design, not length.

**Candidates:**
- KILKLYSSKKY (11-mer, K392=0.892, delta=+0.427)
- VVLVWIFPKKK (11-mer, K392=0.882, delta=+0.200)
- VRLPWVSSKKY (11-mer, K392=0.874, delta=+0.201)

**Structures already available** from the N392 scaffold screen Boltz-2 runs. Would need grafting into full proteins.

## Option 3: Cyclize the 6-mer

**Cost:** ~$500-800 per peptide (synthesis) + ~$1 compute | **Time:** 1-2 weeks (synthesis) | **Probability:** Medium

A linear 6-mer is floppy — high entropy penalty for binding. A head-to-tail cyclic peptide is pre-organized. Benefits:
- Pre-pays conformational entropy cost → better binding free energy
- More drug-like (oral bioavailability, protease resistance)
- Smaller conformational ensemble → more predictable behavior

**Approach:**
- Use RDKit/OpenMM to model cyclic IAFSAF
- Re-dock cyclic version with Boltz-2
- If ipTM holds, synthesize cyclic version

**Caveat:** Boltz-2 may not handle cyclic peptides well (trained on linear). May need specialized docking.

## Option 4: D-amino Acid Backbone

**Cost:** ~$200 per peptide (synthesis surcharge) + ~$0.50 compute | **Time:** ~1 week | **Probability:** Medium

Switch to all-D amino acids. Backbone geometry flips, sidechain chemistry stays. Benefits:
- Protease-resistant (enzymes can't cleave D-peptides)
- Different backbone conformations may fit the channel better
- Can be modeled computationally before synthesis

**Approach:**
- PepINVENT (AstraZeneca, open source) can optimize D-amino acid peptides via RL
- Or manually build D-IAFSAF and re-dock

## Option 5: Non-Natural Anchor Residue

**Cost:** ~$1,000+ per peptide | **Time:** 2-4 weeks | **Probability:** Medium-Low

Replace one position with a non-natural amino acid that forms a stronger interaction:
- Boronic acid at P1 → reversible covalent bond with catalytic zinc
- N-methyl amino acid → reduced backbone flexibility
- Beta-amino acid insertion → altered backbone geometry

This is Phase 2 optimization — only pursue after confirming the binding site is correct.

## Option 6: Run All 18 Short-Stack with Restrained Equilibration

**Cost:** ~$6 (18 × 2 hrs on RTX 4090) | **Time:** ~36 hours | **Probability:** Medium

Before giving up on 6-mers entirely, re-run the full 18-peptide panel with proper restrained equilibration. The current MD failure might be a protocol issue, not a peptide issue. Some of the 18 may be stable when properly equilibrated even if IAFSAF isn't.

---

## Decision Tree

```
IAFSAF drifts in MD
    │
    ├─→ Try restrained equilibration on IAFSAF (Option 1)
    │       │
    │       ├─→ IAFSAF stays → 6-mers work, protocol was the issue
    │       │       └─→ Run all 18 with restrained equil (Option 6)
    │       │
    │       └─→ IAFSAF still drifts → 6-mers genuinely don't bind stably
    │               │
    │               ├─→ Run KILKLYSSKKY 11-mer through MD (Option 2)
    │               │       │
    │               │       ├─→ 11-mer stays → pivot to 11-mers, accept ERAP1 issue
    │               │       └─→ 11-mer drifts → deeper problem with starting poses
    │               │
    │               └─→ Cyclize IAFSAF (Option 3) — if 6-mer length is critical for ERAP1 evasion
    │
    └─→ Cost to reach answer: ~$0.50-3.00 depending on path
```

## Priority Order
1. Restrained equilibration (cheapest, fastest, most likely to work)
2. 11-mer fallback (proven scaffold, known selectivity)
3. Cyclization (if 6-mer length is non-negotiable for ERAP1 evasion)
4. D-amino acids (Phase 2 optimization)
5. Non-natural amino acids (Phase 2+)

---

## File Locations
- This document: `docs/CONTINGENCY_IF_6MER_DRIFTS.md`
- MD protocol: `scripts/openmm_md_protocol.py`
- Starting structures: `data/results/v43_validation/md_starting_structures/`
- 11-mer scaffold data: `data/results/n392_scaffold_search/scaffold_screen_results.json`
- Full validation plan: `docs/V43_WETLAB_VALIDATION_PLAN.md`
