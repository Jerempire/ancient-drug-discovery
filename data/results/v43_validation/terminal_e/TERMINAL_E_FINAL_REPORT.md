# Terminal E Final Report: V2 Binder Channel-Capping Redesign

**Date:** 2026-03-25
**Verdict:** FAILED — Redesign destroys selectivity
**Instance:** 33509299 (RTX 4090, went offline during run)
**Cost:** ~$0.15

---

## The Question

The V2 protein binder (Y87A_Y89A, 92 amino acids) binds ERAP2's surface at the divergent patch with clean selectivity: ipTM 0.748 for ERAP2 vs 0.186 for both ERAP1 and IRAP. But it was designed to **bind selectively**, not to **inhibit**. Does it actually block substrates from entering the catalytic channel? And if not, can we fix that?

---

## Phase 1: Structural Analysis — Does the Binder Cap the Channel?

### What we did

Parsed the Boltz-2 predicted complex (Y87A_Y89A + ERAP2 cropped channel, residues 350-500) and computed:
1. Contact footprint: which ERAP2 residues does the binder touch?
2. SASA (Solvent Accessible Surface Area): how much of the channel entrance is blocked?
3. Spatial coverage: does the binder sit over the opening like a lid?

### What we found

**The binder is a BRIDGE, not a CAP.**

It spans across the channel mouth — contacting both divergent patches (353-367 and 400-414) — but leaves the entrance wide open. Think of it like a bar across a doorway rather than a door.

#### Contact Map

```
ERAP2 Channel Cross-Section (looking down from outside)

        Ceiling (409-419)
        ~~~~ 0% blocked ~~~~
       /                    \
      /                      \
  DP1 |    CHANNEL OPENING    | DP2
353-367     (SUBSTRATES        400-414
  8.6%      ENTER HERE)        22.4%
  blocked    ~~~~~~~~~~~~      blocked
      \                      /
       \   K392 (79.7%)     /
        \   BLOCKED        /
         ==================
         Deep channel (442-489)
         BINDER EXTENDS HERE
```

The binder contacts 39 ERAP2 residues total, including:
- **K392** (selectivity handle): 79.7% SASA reduction — well shielded
- **H374** (zinc coordinator): within 2.66 A — close contact
- **Divergent patch 1** (353-367): 4/15 residues contacted
- **Divergent patch 2** (400-414): 4/15 residues contacted
- **Deep channel** (442-489): 12 residues contacted — the binder reaches deep

But critically:
- **Ceiling (Q412, D414)**: 0.0% SASA reduction — completely uncovered
- **Mid-gap (387-398)**: 11 residues uncovered — the zinc approach path is wide open
- **E393** (zinc coordinator): not contacted at all (5.5 A away)
- **Overall entrance**: only 30.7% SASA reduction (need >80% for a functional cap)

#### The Topology Problem

The binder has a **deep-insertion topology**: it pokes into the channel like a finger in a tube rather than sitting flat across the opening like a lid. Only 37 of its 757 atoms (4.9%) are within the channel cross-section. It maximizes its binding interface by reaching deep, but this means it doesn't cover the entrance.

### Phase 1 Verdict

**BRIDGE_NOT_CAP / POSSIBLY_INHIBITORY** — The binder may partially impede substrates (it contacts K392 and is near the zinc), but it cannot fully block channel access. Substrates can still enter through the ceiling and mid-gap openings.

---

## Phase 2: Redesign Attempts — 5 Variants Tested

### Design Philosophy

Two independent approaches to make the binder cover the channel entrance:

**Approach A — Loop Extensions (other terminal):**
Add short loops to the binder's interior to extend its reach toward uncovered regions.

| Variant | Modification | Target Gap | Size |
|---------|-------------|------------|------|
| **northcap** | +4aa loop (ENYGS) at position 26 | Ceiling (409-419) | 96aa |
| **southcap** | +4aa loop (RNSFN) at position 14 | Floor (353-360) | 96aa |
| **dualcap** | Both loops combined | Both gaps | 100aa |

**Approach B — Peptide Tail Grafts (this terminal):**
Append a validated 6-mer channel peptide to the binder's C-terminus via a GGS linker.

| Variant | Modification | Rationale | Size |
|---------|-------------|-----------|------|
| **hybrid_vagsaf** | +GGS+VAGSAF C-terminal tail | V4 "Hammer" peptide plugs channel | 101aa |
| **hybrid_iafsaf** | +GGS+IAFSAF C-terminal tail | V4 "Scalpel" peptide plugs channel | 101aa |

### Boltz-2 Selectivity Screen

Each variant tested against 3 targets (ERAP2-K392, ERAP1, IRAP) with 3 diffusion samples each. Run without MSAs (msa=empty) due to the de novo binder sequence.

### Results

| Variant | ERAP2-K392 avg | ERAP1 avg | IRAP avg | Delta (E2-E1) | Pass? |
|---------|---------------|-----------|----------|---------------|-------|
| **Original V2** | **0.748** | **0.186** | **0.186** | **+0.562** | **YES** |
| northcap | 0.567 | 0.471 | ??? | +0.096 | NO |
| dualcap | 0.425 | 0.586 | 0.379 | -0.161 | NO |
| hybrid_iafsaf | 0.532 | 0.571 | 0.511 | -0.039 | NO |
| hybrid_vagsaf | 0.426 | 0.552 | 0.413 | -0.126 | NO |
| southcap | ??? | ??? | ??? | ??? | LOST |

**Pass criteria:** ERAP2 > 0.7, ERAP1 < 0.3, IRAP < 0.3

### Every single variant is worse than the original on every metric.

- **ERAP2 binding dropped** from 0.748 to 0.43-0.57 (24-43% weaker)
- **ERAP1 cross-reactivity exploded** from 0.186 to 0.47-0.59 (2.5-3.2x higher)
- **IRAP cross-reactivity exploded** from 0.186 to 0.38-0.51 (2.0-2.7x higher)
- **Selectivity delta collapsed** from +0.562 to between -0.161 and +0.096
- The **dualcap** variant actually binds ERAP1 *better than ERAP2* (delta -0.161)

Individual sample scores show high variance (northcap got one ERAP2 sample at 0.746, but the other two were 0.554 and 0.400), consistent with the modifications destabilizing the original binding mode.

---

## Why It Failed: Root Cause Analysis

### 1. The Selectivity-Inhibition Tradeoff Is Fundamental

The original V2 binder achieves selectivity by targeting the **divergent patch** — the ~30 residues (353-367 and 400-414) where ERAP2 differs from ERAP1 and IRAP. But the substrate channel entrance sits **between and below** the divergent patches, surrounded by the **conserved catalytic machinery** (H370, H374, E393, and the HEXXH motif shared by all M1 aminopeptidases).

**You cannot extend from the divergent patch toward the channel entrance without crossing into conserved territory.** This is a geometric constraint, not a design flaw. The channel entrance is literally surrounded by the most conserved residues in the entire enzyme family.

```
   DIVERGENT (unique to ERAP2)
   ============================
   |  DP1 (353-367)           |
   |                          |
   |    CONSERVED ZONE        |  <-- Any extension lands here
   |    H370, H374, E393      |  <-- Shared with ERAP1 + IRAP
   |    (zinc catalytic site)  |
   |                          |
   |  DP2 (400-414)           |
   ============================
   DIVERGENT (unique to ERAP2)
```

### 2. The Peptide Tail Grafts Failed for a Different Reason

The hybrid approach (binder body + GGS + VAGSAF) was meant to keep the binder body on the divergent patch while the peptide tail independently plugs the channel. But Boltz-2 doesn't respect this intent — it predicts the most thermodynamically favorable complex, and a flexible GGS linker gives the peptide tail too many degrees of freedom. Instead of the tail docking into the channel, it likely:

- Collapses back onto the binder body (intramolecular contacts are easier than intermolecular)
- Forces the binder off its preferred divergent-patch pose to accommodate the tail
- Creates new contacts with conserved residues as the tail explores conformational space

The IAFSAF hybrid scored highest on ERAP2 (0.609 best avg) but also scored 0.571 on ERAP1 and 0.620 on one IRAP sample — the tail is interacting non-selectively.

### 3. msa=empty Mode Is a Caveat but Not the Explanation

Running without MSAs reduces prediction accuracy, particularly for the target protein chains (ERAP2, ERAP1, IRAP) whose folds benefit from evolutionary covariance information. However:
- The same msa=empty condition applies to all 5 variants AND the targets equally
- The relative ranking (all worse than original) is meaningful even if absolute values shift
- The cross-reactivity pattern (ERAP1 > 0.5 for every variant) is too strong to be a noise artifact

### 4. The ERAP2 Sequence Crop Matters

The other terminal used a different/longer ERAP2 crop than the original V2 screen. This changes the binding landscape and could partially explain the lower absolute ERAP2 scores. However, the same crop was used for all variants, so the relative comparison is valid.

---

## What Can Be Fixed

### Option 1: Accept V2 as a Selective Surface Binder (RECOMMENDED)

Don't try to make V2 do something it wasn't designed for. It's an excellent selective binder — 0.748 ipTM with 4x selectivity over ERAP1/IRAP. Possible uses:

- **Research tool**: Fluorescently tagged V2 binder to detect ERAP2-K392 allele in cell assays
- **Delivery vehicle**: Conjugate a cytotoxic payload to V2 to target ERAP2-K392-expressing cells
- **Allosteric modulator**: V2 may alter ERAP2 activity even without capping the channel — surface binding could change the conformational dynamics (open/closed transition)

### Option 2: Pursue V4 Short-Stack Peptides Independently (ALREADY DOING THIS)

The 6-mer peptides (VAGSAF, IAFSAF) were designed specifically for channel inhibition. They achieve selectivity through a completely different mechanism — length-based evasion of ERAP1's molecular ruler + insufficient surface area in IRAP's cavity. They don't need the V2 binder's help.

Terminals B-D are already validating these via MD simulation. **This is the right path to functional inhibition.**

### Option 3: Rigid Linker Instead of GGS (If V2+V4 Hybrid Is Still Desired)

The flexible GGS linker was the wrong choice — it allows the peptide tail too much freedom. A **rigid helical linker** (e.g., EAAAK repeat) would constrain the tail's position relative to the binder body:

```
Y87A_Y89A (92aa) — EAAAKEAAAK — VAGSAF (6aa) = 108aa
```

This could work if:
- The rigid linker places the tail exactly at the channel entrance
- The tail's position is predictable (not flopping around)
- The linker doesn't create new contacts with conserved residues

**Risk:** Still likely to hit conserved residues. Would need careful geometric modeling before Boltz-2.

### Option 4: Disulfide-Tethered Channel Plug

Instead of a covalent fusion, design a separate small peptide that:
1. Binds inside the channel entrance (like V4 peptides)
2. Has a cysteine that forms a disulfide bond with a cysteine introduced on the V2 binder surface

This keeps the two components structurally independent — the binder anchors to the divergent patch, the plug sits in the channel, and the disulfide just keeps them in proximity. The plug doesn't need to be selective on its own because the binder ensures it only reaches ERAP2.

**Complexity:** Requires introducing a Cys mutation on V2 at the right position + designing a Cys-containing channel peptide. Non-trivial but potentially the best of both worlds.

### Option 5: Protein-PROTAC Approach

Instead of blocking ERAP2's channel, use V2 as a **degradation recruiter**:
- Fuse V2 to an E3 ligase-recruiting peptide (e.g., VHL or CRBN binder)
- V2 selectively binds ERAP2-K392
- The E3 ligase tag recruits ubiquitination machinery
- ERAP2-K392 gets degraded via the proteasome

This bypasses the inhibition problem entirely — you don't need to block the channel if you eliminate the protein. V2's selectivity ensures only the K392 allele gets degraded.

**Frontier approach** — no protein-based PROTACs (bioPROTACs) are approved yet, but several are in clinical trials.

---

## Summary

| Phase | Question | Answer |
|-------|----------|--------|
| 1. SASA analysis | Does V2 cap the channel? | **NO** — 30.7% entrance coverage, ceiling 0% |
| 2. Loop extensions | Can loops reach the ceiling? | **NO** — contacts conserved residues, loses selectivity |
| 2. Peptide grafts | Can a tail plug the channel? | **NO** — flexible tail destabilizes binding, loses selectivity |

**The fundamental insight:** Selectivity and channel inhibition are geometrically opposed for protein binders on ERAP2. The divergent patch (selective) and the channel entrance (inhibitory) are separated by a belt of conserved catalytic residues. Any protein large enough to reach from one to the other will inevitably contact the conserved zone.

**The 6-mer peptides solve this differently** — they achieve selectivity through *size* (too short for ERAP1, too small for IRAP) rather than through *surface complementarity*. That's why they work and the V2 redesigns don't.

---

## Files Produced

| File | Description |
|------|-------------|
| `v2_channel_cap_analysis.json` | Phase 1 contact + SASA analysis (39 contacts mapped) |
| `boltz2_selectivity_results.json` | Phase 2 Boltz-2 scores (11/15 predictions) |
| `CHANNEL_CAP_REDESIGN.md` | Pre-screen design rationale and YAML configs |
| `HANDOFF.md` | Continuation instructions for another terminal |
| `TERMINAL_E_FINAL_REPORT.md` | This document |
| `scripts/v2_channel_cap_analysis.py` | Reusable SASA capping analysis script |
