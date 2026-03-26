# V2 Binder Channel-Capping Redesign — Terminal E

**Date:** 2026-03-25
**Status:** ANALYSIS COMPLETE, REDESIGN PROPOSALS READY

---

## Current State: BRIDGE_NOT_CAP

The Y87A_Y89A binder (92aa) binds ERAP2 at residues 350-489 but **does not cap the substrate channel entrance**. It bridges both divergent patches while leaving the channel opening exposed.

### SASA Analysis Results

| Region | Apo SASA | Holo SASA | Reduction |
|--------|----------|-----------|-----------|
| Channel entrance (350-419) | 3786 A² | 2625 A² | **30.7%** |
| Divergent patch 1 (353-367) | 1043 A² | 954 A² | 8.6% |
| Divergent patch 2 (400-414) | 814 A² | 631 A² | 22.4% |
| Zinc coord (H370,H374,E393) | 30 A² | 16 A² | 44.5% |
| K392 selectivity | 32 A² | 6 A² | **79.7%** |
| Ceiling (Q412, D414) | 189 A² | 189 A² | **0.0%** |

### Diagnosis

The binder has a **deep-insertion topology**: it contacts K392 and extends deep into the channel (442-489) like a finger poking into a tube, rather than sitting across the entrance like a lid. This leaves three critical gaps:

1. **Ceiling gap (409-419)**: 11 contiguous residues, 0% blocked. This is the largest opening.
2. **Mid-channel gap (387-398)**: The zinc approach path. E393 (zinc coord) is not contacted.
3. **N-terminal patch gap (353-360)**: 8 residues of divergent patch 1 uncovered.

**Why it doesn't cap:** The binder was designed for *selective binding* to the divergent patch, not for *functional inhibition*. It found a thermodynamically favorable pose that maximizes interface area by reaching deep, but this leaves the entrance mouth open.

---

## Redesign Strategy Options

### Option A: Hybrid V2+V4 — Graft a Peptide Plug (RECOMMENDED)

**Concept:** Attach a 6-mer peptide tail (VAGSAF or IAFSAF) to the V2 binder's C-terminus via a short flexible linker. The protein body anchors to the divergent patch (providing selectivity), while the peptide tail plugs the channel entrance (providing inhibition).

**Rationale:**
- V2 already contacts K392 and bridges both patches — selectivity is solved
- V4 short-stack peptides (VAGSAF, IAFSAF) dock into the channel entrance at the zinc region
- Combining them: protein anchor + peptide plug = selective channel cap
- The C-terminus of the binder is near residue 92, which contacts ERAP2 around 399-408 — close to the channel entrance

**Design:**
```
Y87A_Y89A binder (92aa) — GGS linker — VAGSAF (6aa)
Total: ~101 aa
```

**Boltz-2 validation:** Predict the chimeric construct bound to ERAP2 and check:
1. Does the binder body maintain its pose on the divergent patch?
2. Does the peptide tail insert into the channel entrance?
3. Does SASA of channel entrance + zinc drop to >80%?

**Advantages:**
- Builds on two validated components (V2 selectivity + V4 peptide binding)
- Short linker keeps peptide in the right neighborhood
- Doesn't require redesigning the binder core

**Risks:**
- Linker may be too flexible → peptide doesn't find the entrance
- Steric clash between binder body and peptide tail
- May lose selectivity if the peptide forces the binder off its preferred pose

### Option B: Loop Extension — Ceiling Coverage

**Concept:** Insert a 5-8 residue loop into the binder at a position that faces the ceiling gap (409-419). The loop would extend from the binder surface to cover the uncovered ceiling residues.

**Target:** The binder contacts 399, 402, 403, 407, 408 — the edge of the ceiling gap. A loop inserted between binder residues 50-60 (which face the 400-408 region of ERAP2) could reach residues 409-419.

**Approach:**
1. Identify binder residues closest to the ceiling gap
2. Insert a loop (GGSGGS or designed sequence) at that position
3. Re-predict with Boltz-2 using ceiling residues (409-419) as hotspot constraints

**Advantages:**
- Preserves the binder's core fold and binding mode
- Targeted — only adds coverage where needed

**Risks:**
- Loop insertion may destabilize the binder fold
- Boltz-2 may not reliably place the inserted loop at the ceiling

### Option C: Rotational Repositioning

**Concept:** The binder's deep contacts (442-489) waste surface area on buried residues that don't need blocking. Rotating the binder ~30-45 degrees around the channel axis would shift those contacts from deep-interior toward the ceiling.

**Approach:**
1. Re-dock with Boltz-2 using modified hotspot constraints:
   - ADD hotspots: 409, 412, 414 (ceiling)
   - ADD hotspots: 388, 393, 395, 398 (mid-gap)
   - REMOVE deep hotspots (no constraint on 442-489)
2. Use `pocket_ids` or equivalent constraint to bias the binder toward the entrance plane

**Advantages:**
- No sequence changes to the binder
- Leverages Boltz-2's ability to find new poses

**Risks:**
- May find a completely different pose rather than a rotated version
- Could lose K392 contact (selectivity loss)
- Boltz-2 hotspot constraints don't guarantee the binder will respect them

### Option D: De Novo Channel Lid (via Proteina-Complexa)

**Concept:** Use Proteina-Complexa to design a new binder specifically targeting the channel entrance plane. Specify the entrance residues (350-419) as the binding surface and constrain the binder to sit across the opening.

**Advantages:**
- Purpose-built for capping, not retrofitted
- Can specify exact interface residues

**Risks:**
- Starts from scratch — loses V2's proven selectivity
- Proteina-Complexa may not handle the concave entrance geometry well (same issue as RFdiffusion V3)

---

## Recommended Execution Order

1. **Option A (Hybrid V2+V4)** — fastest, lowest risk, builds on two validated components
2. **Option C (Repositioning)** — if hybrid doesn't work, try biased re-docking
3. **Option B (Loop extension)** — more invasive, try if A and C fail
4. **Option D (De novo)** — last resort

---

## Boltz-2 YAML Configs for Option A

### Config 1: Chimeric binder (Y87A_Y89A + GGS + VAGSAF) vs ERAP2

```yaml
# Save as: data/boltz2_inputs/v2_hybrid_vagsaf.yaml
sequences:
  - protein:
      id: "erap2_channel"
      sequence: "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
      # ERAP2 residues 350-500 (cropped channel)
  - protein:
      id: "v2_hybrid_vagsaf"
      sequence: "[V2_BINDER_92AA_SEQUENCE]GGSVAGSAF"
      # Y87A_Y89A binder + GGS linker + VAGSAF peptide
```

### Config 2: Same with IAFSAF tail

```yaml
sequences:
  - protein:
      id: "erap2_channel"
      sequence: "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
  - protein:
      id: "v2_hybrid_iafsaf"
      sequence: "[V2_BINDER_92AA_SEQUENCE]GGSIAFSAF"
```

**NOTE:** Replace `[V2_BINDER_92AA_SEQUENCE]` with the actual binder sequence extracted from the CIF file (chain B).

### Config 3: Repositioned binder (Option C) with ceiling hotspots

```yaml
sequences:
  - protein:
      id: "erap2_channel"
      sequence: "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
  - protein:
      id: "v2_binder_reposition"
      sequence: "[V2_BINDER_92AA_SEQUENCE]"
# Use --diffusion_samples 5 --seed 42 for multiple poses
# Score all 5 and pick the one with best ceiling coverage
```

---

## Success Criteria for Redesigned Binder

| Metric | Target | Why |
|--------|--------|-----|
| Channel entrance SASA reduction | >80% | Functional plug — substrates can't enter |
| Zinc coordination SASA reduction | >80% | Catalytic site blocked |
| Ceiling (Q412, D414) SASA reduction | >50% | Currently 0% — must improve |
| K392 SASA reduction | >50% | Maintain selectivity handle |
| Bridges both divergent patches | YES | Required for selectivity |
| ERAP2 ipTM | >0.6 | Confident binding |
| ERAP1 ipTM | <0.3 | Selectivity maintained |
| IRAP ipTM | <0.3 | Selectivity maintained |

---

## Next Steps

1. Extract V2 binder sequence from CIF file (chain B)
2. Generate hybrid constructs (binder + linker + peptide)
3. Run Boltz-2 predictions on Vast.ai (GPU, ~$0.50)
4. Re-run SASA analysis on new poses
5. If >80% entrance SASA reduction → proceed to MD validation
6. Run selectivity counterscreen (ERAP1, IRAP) on best hybrid
