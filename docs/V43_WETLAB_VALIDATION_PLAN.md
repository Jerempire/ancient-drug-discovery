# V4.3 Wet Lab Validation Plan — Short-Stack Peptide Leads

**Status:** READY FOR EXECUTION
**Date:** 2026-03-24
**Author:** Claude Code + Jeremy
**Goal:** Prepare 4 peptide leads for $2,000 wet lab validation + provisional patent filing

---

## Background

### The Discovery
6-mer peptides achieve triple selectivity against ERAP2-K392:
- **ERAP1 evasion**: Below 8-residue molecular ruler floor (avg ERAP1 ipTM: 0.266 vs K392: 0.805)
- **IRAP evasion**: Insufficient surface area in IRAP's cavernous active site (avg IRAP ipTM: 0.298)
- **K392 selectivity**: Steric packing against K392 lysine's 4-carbon aliphatic stem

### The "Elite Four" Test Panel

| Code Name | Sequence | Length | K392 | N392 | ERAP1 | IRAP | Purpose |
|---|---|---|---|---|---|---|---|
| **The Hammer** | VAGSAF | 6 | 0.905 | 0.870 | 0.335 | 0.236 | Highest K392 binder |
| **The Scalpel** | IAFSAF | 6 | 0.870 | 0.631 | 0.192 | 0.251 | Best allele-selective (+0.239 delta) |
| **The Trp Edge** | VAWSAF | 6 | 0.852 | 0.830 | 0.169 | 0.226 | Does Trp add residence time? |
| **Scrambled Ctrl** | FASGAV | 6 | TBD | TBD | TBD | TBD | Catches sticky-pocket bias |

### Previous Versions That Failed (DO NOT SUGGEST THESE)
- V1: Active-site binders (90%+ conserved with ERAP1, zero selectivity)
- V2: De novo protein binders (Y87A/Y89A works but hard to manufacture)
- V3: RFdiffusion can't design for concave channels
- V4.1: 11-mer Trp/Phe wedges failed to kill IRAP (cavity too permissive)

---

## TASK 1: Functional Check — SASA Zinc Capping Test

### Question
Does a 6-mer peptide actually CAP the catalytic zinc in ERAP2's active site, or does it just sit passively in the channel without blocking enzyme function?

### Method
Use MDAnalysis to compute Solvent Accessible Surface Area (SASA) of the catalytic zinc (Zn2+) and its coordinating residues with and without the peptide bound.

### ERAP2 Catalytic Zinc Coordination
- **Zinc position**: Coordinated by H370, H374, E393 (HEXXH motif)
- **Catalytic mechanism**: Zinc activates a water molecule for peptide bond hydrolysis
- **If SASA of zinc drops >80% when peptide is bound** → peptide is a functional plug
- **If SASA barely changes** → peptide is a passenger, not an inhibitor

### Input Files Needed
- Boltz-2 predicted CIF/PDB of each peptide + ERAP2-K392 complex
- Must use FULL ERAP2 structure (not cropped channel) for accurate zinc environment
- AlphaFold ERAP2 structure: `data/structures/erap2_wt_alphafold.pdb`

### Script Location
`scripts/sasa_zinc_check.py`

### Script Specification

```python
"""
SASA Zinc Capping Test

For each Elite Four peptide:
1. Load Boltz-2 predicted complex (ERAP2 + peptide)
2. Load apo ERAP2 (no peptide)
3. Compute SASA of zinc + coordinating residues (H370, H374, E393)
4. Compare: if peptide reduces zinc SASA >80%, it's a functional plug

Dependencies: MDAnalysis, numpy
Install: pip install MDAnalysis numpy

Input: PDB/CIF files from Boltz-2 predictions
Output: JSON with SASA values + plug/passenger verdict
"""
```

### Expected Runtime
~1 min per structure, CPU only. No GPU needed.

### Success Criteria
- SASA reduction >80%: **FUNCTIONAL PLUG** — peptide blocks substrate access to zinc
- SASA reduction 50-80%: **PARTIAL PLUG** — may reduce activity but not fully inhibit
- SASA reduction <50%: **PASSENGER** — peptide sits in channel but doesn't block catalysis

---

## TASK 2: MD Stability Protocol — OpenMM 10ns Simulations

### Question
Does each peptide STAY locked in ERAP2-K392 for 10ns while DRIFTING out of IRAP?

### The 3-Target Stress Test
For each of the 4 Elite Four peptides, run MD against:
1. **ERAP2-K392** — confirm persistence (low RMSD, peptide stays)
2. **IRAP** — confirm evasion (high RMSD, peptide drifts out)
3. **ERAP1** — confirm invisibility (peptide never finds stable pose)

**Total: 12 MD simulations** (4 peptides × 3 targets)

### System Preparation

#### Starting Structures
Use Boltz-2 predicted complexes. For ERAP1 runs, dock peptides into ERAP1 channel using the Boltz-2 YAML approach first (already done in counterscreen), then use those poses as MD starting points.

#### Force Field
- Protein: `amber14-all.xml`
- Water: `amber14/tip3pfb.xml`
- Solvent: explicit water box, 1.0 nm padding
- Ions: 150 mM NaCl (physiological)
- No zinc parameters needed for cropped channel (zinc is outside the 350-450 crop)

#### If Using Full ERAP2 (recommended for zinc capping):
- Use ZAFF (Zinc Amber Force Field) parameters for catalytic zinc
- Or use `openmmforcefields` with `gaff-2.11` for non-standard

### Simulation Protocol

```
Phase 1: Energy Minimization
  - 1000 steps steepest descent
  - Tolerance: 10 kJ/mol

Phase 2: NVT Equilibration (100 ps)
  - Temperature: 310 K (physiological)
  - Langevin integrator, friction 1.0 /ps
  - Timestep: 2 fs
  - Hydrogen mass repartitioning: YES (allows 4 fs timestep in production)
  - Restraints: 10 kJ/mol/nm² on protein backbone heavy atoms

Phase 3: NPT Equilibration (100 ps)
  - Same as NVT + MonteCarloBarostat at 1 atm
  - Restraints: 5 kJ/mol/nm² on protein backbone

Phase 4: Production MD (10 ns)
  - NPT ensemble, 310 K, 1 atm
  - Timestep: 4 fs (with HMR)
  - No restraints
  - Save coordinates every 10 ps (1000 frames total)
  - Save energies every 1 ps
  - PME for electrostatics, 1.0 nm cutoff
```

### Script Location
`scripts/openmm_md_protocol.py`

### Script Specification

```python
"""
OpenMM 10ns MD Protocol for Short-Stack Validation

Usage:
  python openmm_md_protocol.py <complex.pdb> <output_dir> [--target erap2|irap|erap1]

For each run:
1. Load PDB, add hydrogens, solvate
2. Minimize → NVT equil → NPT equil → 10ns production
3. Save trajectory (DCD), energies (CSV), final state (XML)
4. Compute RMSD of peptide relative to starting pose
5. Compute center-of-mass drift of peptide
6. Output: rmsd.csv, com_drift.csv, summary.json

Dependencies: openmm, openmmtools, mdtraj, pdbfixer
Install: conda install -c conda-forge openmm mdtraj pdbfixer

GPU: Required. RTX 4090 does ~150 ns/day for ~10K atom system.
10ns ≈ 1.5-2 hours per simulation.
"""
```

### Analysis: Drift Detection + Functional Validation

After each 10ns run, compute ALL of the following:

#### A. Core Stability Metrics

1. **Peptide RMSD** (backbone atoms, aligned to protein):
   - < 2.0 Å: LOCKED — peptide is stable
   - 2.0 - 5.0 Å: SHIFTING — peptide moves but stays in pocket
   - > 5.0 Å: DRIFTING — peptide is leaving

2. **Center-of-Mass (COM) Drift**:
   - Track peptide COM distance from starting position over time
   - Monotonically increasing COM = peptide exiting
   - Stable COM with fluctuations = peptide bound

3. **Contact Fraction**:
   - Count native contacts (contacts present in frame 1) maintained over time
   - > 70% maintained: STABLE complex
   - < 30% maintained: DISSOCIATED

#### B. Zinc Distance (CRITICAL — proves functional inhibition)

4. **Peptide-to-Zinc Distance**:
   - Track minimum distance from any peptide heavy atom to catalytic Zn2+ over time
   - < 5 Å: Peptide is near catalytic center — potential functional plug
   - 5-10 Å: Peptide is in channel but not blocking zinc
   - > 10 Å: Peptide is peripheral, not inhibitory
   - **Must use FULL ERAP2 structure** (not cropped channel) for this metric

#### C. 3-Region Contact Distribution (proves NOT a local anchor)

5. **Contact Distribution Across 3 Channel Neighborhoods**:
   A 6-mer must touch 3 distinct regions to be a real binder, not a local anchor:

   | Region | Name | Key Residues | What It Proves |
   |---|---|---|---|
   | **Floor** | K392 selectivity handle | K392, Y398 | Variant selectivity |
   | **Wall** | IRAP/ERAP1 evasion handle | A403, A406 | Paralog selectivity |
   | **Ceiling** | Channel cap | Q412, D414 | Peptide fills volume, not laying flat |

   For each frame, count contacts (< 4.5 Å) with each region:
   - All 3 regions contacted: **VOLUME FILLER** — real functional binder
   - 2 regions: **PARTIAL** — may still work but less confident
   - 1 region only: **LOCAL ANCHOR** — likely not a functional inhibitor

#### D. Hydrogen Bond Occupancy (proves specific interactions)

6. **Key H-bond Occupancy Over Trajectory**:
   - P1 ↔ K392 sidechain: Is there a persistent interaction? (>50% occupancy)
   - C-terminal Phe backbone ↔ channel residues: Does the aromatic anchor hold?
   - Any peptide ↔ Y398 hydroxyl: Does the ERAP2-unique Tyr contribute?
   - Report: H-bond donor, acceptor, occupancy %, avg distance

#### E. Per-Run Summary JSON (REQUIRED — all terminals must output this)

Each MD run's `summary.json` in the terminal output directory must contain:
```json
{
  "peptide": "VAGSAF",
  "target": "erap2k392",
  "terminal": "b",
  "rmsd_final": 1.45,
  "rmsd_mean": 1.23,
  "rmsd_max": 2.10,
  "com_drift_final": 0.82,
  "contact_fraction_final": 0.85,
  "zinc_distance_mean": 4.2,
  "zinc_distance_min": 3.1,
  "floor_contacts_pct": 0.92,
  "wall_contacts_pct": 0.78,
  "ceiling_contacts_pct": 0.65,
  "region_verdict": "VOLUME_FILLER",
  "hbond_p1_k392_occupancy": 0.73,
  "hbond_cterm_occupancy": 0.81,
  "verdict": "LOCKED",
  "simulation_ns": 10.0,
  "wall_time_hours": 1.8,
  "gpu": "RTX 4090",
  "instance_id": "33XXXXXX"
}
```

### Script Location for Analysis
`scripts/md_drift_analysis.py`

### FUTURE: Steered MD Competition Test (Phase 2, after initial MD)

If leads pass the 10ns stability test, run Steered Molecular Dynamics (SMD):
- Place a natural 10-mer substrate at the channel entrance while the 6-mer is bound
- Apply a pulling force to drag the substrate into the channel
- **If the 6-mer doesn't budge**: it has higher binding affinity than substrate → true competitive inhibitor
- **If the 6-mer gets displaced**: it's a weaker binder that won't survive in vivo substrate competition
- This is the computational equivalent of a Ki assay
- Requires ~5ns per pull, 3 replicas = ~15ns per peptide. Do AFTER confirming stability.

### Expected Results

| Peptide | ERAP2-K392 | IRAP | ERAP1 |
|---|---|---|---|
| VAGSAF | LOCKED (<2Å) | DRIFTING (>5Å) | DRIFTING (>5Å) |
| IAFSAF | LOCKED (<2Å) | DRIFTING (>5Å) | DRIFTING (>5Å) |
| VAWSAF | LOCKED (<2Å) | DRIFTING (>5Å) | DRIFTING (>5Å) |
| FASGAV | DRIFTING (>5Å) | DRIFTING (>5Å) | DRIFTING (>5Å) |

If FASGAV (scrambled control) also locks in ERAP2, the pocket is artificially sticky and the results are unreliable.

### Runtime Estimate
- 12 simulations × 2 hrs each = ~24 hrs on single RTX 4090
- Can parallelize across 2-3 instances: ~8-12 hrs
- Cost: ~$6-12 on Vast.ai

---

## TASK 3: CRO Readiness — Peptide Synthesis Specifications

### Order Specifications

**Vendor Options:** GenScript, AAPPTEC, Bachem, CPC Scientific

**For each peptide, order:**

| Parameter | Specification |
|---|---|
| Form | Naked peptide (no modifications unless noted) |
| Salt form | TFA (trifluoroacetic acid) salt |
| Purity | >95% by HPLC |
| Quantity | 5 mg per peptide (enough for IC50 + Ki + SPR) |
| QC included | HPLC chromatogram + MALDI-TOF MS |
| Turnaround | Standard (2-3 weeks) |

**Peptide Order Table:**

| # | Name | Sequence | MW (approx) | Notes |
|---|---|---|---|---|
| 1 | VAGSAF | Val-Ala-Gly-Ser-Ala-Phe | ~536 Da | Lead - highest K392 binding |
| 2 | IAFSAF | Ile-Ala-Phe-Ser-Ala-Phe | ~624 Da | Lead - best allele-selectivity |
| 3 | VAWSAF | Val-Ala-Trp-Ser-Ala-Phe | ~651 Da | Trp residence time test |
| 4 | FASGAV | Phe-Ala-Ser-Gly-Ala-Val | ~536 Da | Scrambled control (VAGSAF shuffled) |

**Estimated cost:** $100-200 per peptide at 5 mg, >95% purity = **$400-800 total**

### D-Amino Acid Variants (Phase 2, after initial validation)

If L-amino acid versions show activity, order D-amino acid versions for protease resistance:

| # | Name | Modifications | Notes |
|---|---|---|---|
| 5 | D-VAGSAF | All D-amino acids | Protease-resistant version |
| 6 | D-IAFSAF | All D-amino acids | Protease-resistant version |
| 7 | Retro-inverso VAGSAF | D-amino acids, reversed sequence (FASGA-v) | Maintains binding geometry |

**D-amino acid surcharge:** ~$50-100 per peptide additional

### Enzymatic Assay Specifications

**IC50 Assay:**
- **Enzyme:** Recombinant human ERAP2 (R&D Systems #3735-ZN or Sigma #SRP0330)
- **Substrate:** L-Leucine-AMC (Leu-AMC) — standard fluorogenic aminopeptidase substrate
- **Buffer:** 50 mM Tris-HCl pH 7.4, 150 mM NaCl, 0.01% Triton X-100
- **Temperature:** 37°C
- **Peptide concentrations:** 0.1, 0.3, 1, 3, 10, 30, 100, 300 µM (8-point curve)
- **Controls:** DMSO vehicle, bestatin (known aminopeptidase inhibitor, positive control)
- **Readout:** Fluorescence (ex 380/em 460), kinetic read every 30 sec for 30 min
- **Analysis:** Initial rate vs [inhibitor], fit to 4-parameter logistic → IC50

**Selectivity Panel:**
- Same assay run against:
  - ERAP2-K392 (primary, if variant-specific enzyme available)
  - ERAP2-N392 (counterscreen)
  - ERAP1 (counterscreen — R&D Systems #2334-ZN)
  - IRAP (counterscreen — R&D Systems #2680-ZN)

**Ki Determination (for leads with IC50 < 10 µM):**
- Cheng-Prusoff correction: Ki = IC50 / (1 + [S]/Km)
- Requires Km determination for each enzyme (literature values available)

**Expected outcomes:**
- VAGSAF/IAFSAF: IC50 < 50 µM on ERAP2, IC50 > 500 µM on ERAP1/IRAP
- FASGAV (scrambled): IC50 > 500 µM on everything (negative control)
- If IC50 < 10 µM: strong lead, proceed to D-amino acid optimization
- If IC50 10-100 µM: moderate, needs optimization but proves mechanism
- If IC50 > 100 µM: weak, Boltz-2 predictions may not translate to activity

### SPR (Surface Plasmon Resonance) — Optional Phase 2

If enzymatic assay confirms activity:
- **Instrument:** Biacore T200 or Cytiva 8K
- **Immobilize:** ERAP2 on CM5 chip (amine coupling)
- **Analyte:** Peptides at 0.1-100 µM
- **Gives:** kon, koff, KD (binding kinetics + affinity)
- **Key metric:** Residence time (1/koff) — does the peptide STAY bound?

---

## TASK 4: IP Strategy — Patent Framing

### Core Invention
**Length-Dependent Triple Selectivity of Peptide Aminopeptidase Inhibitors**

### Non-Obvious Inventive Step Arguments

#### 1. "Size Floor" Evasion of ERAP1
- **Prior art teaches:** ERAP1 and ERAP2 share conserved active sites; selective inhibition requires targeting divergent surface patches (see V1/V2 failures)
- **Our invention teaches:** ERAP1 has a conformational "molecular ruler" requiring peptides ≥8 residues to trigger catalytic closure. Peptides of 6-7 residues bind ERAP2 (which has a more permissive active site) while being invisible to ERAP1.
- **Non-obvious because:** No prior publication identifies the 8-residue floor as a selectivity handle for drug design. The molecular ruler was characterized for substrate processing, not inhibitor design.

#### 2. "Surface Area Evasion" of IRAP
- **Prior art teaches:** IRAP has a large, cavernous active site that accommodates diverse substrates
- **Our invention teaches:** 6-mer peptides lack sufficient surface area to form stable contacts in IRAP's cavity (IRAP ipTM 0.236-0.411 vs ERAP2 ipTM 0.718-0.905). The same property (small size) that evades ERAP1's ruler also evades IRAP's cavity.
- **Non-obvious because:** Prior inhibitor design focused on steric clashes or charge complementarity for selectivity. No prior art teaches that simply being SMALL achieves selectivity against a larger enzyme.

#### 3. Steric (Not Electrostatic) K392 Variant Selectivity
- **Prior art teaches:** K392 has a positive charge (lysine); drug design would suggest negatively charged inhibitors (salt bridge)
- **Our invention teaches:** The 4-carbon aliphatic STEM of K392's lysine, not its charge, is the primary selectivity handle. Branched hydrophobic P1 residues (Val, Ile) outperform charged residues (Glu) for K392 selectivity.
- **Non-obvious because:** Targeting a charged residue with a hydrophobic interaction is counterintuitive. Our data shows Val at P1 achieves +0.125 corrected selectivity vs Glu at +0.076.

#### 4. Length-Dependent Motif Inversion
- **Our data shows:** The KKY C-terminal motif is K392-selective on 11-mers (delta +0.427) but N392-selective on 7-mers (delta -0.116). The SAME motif produces OPPOSITE selectivity depending on peptide length.
- **Non-obvious because:** This is a novel structure-activity relationship with no precedent in aminopeptidase literature. It demonstrates that the channel microenvironment changes at different depths.

### Claim Structure (Draft)

**Independent Claims:**
1. A peptide of 5-7 amino acid residues that selectively inhibits ERAP2 over ERAP1 and IRAP, wherein the peptide is below the conformational gating threshold of ERAP1.
2. A peptide according to claim 1, wherein position 392 of the target ERAP2 is lysine (K392 variant).
3. A method of selectively inhibiting ERAP2-K392 aminopeptidase activity comprising administering a peptide of 5-7 residues having a C-terminal aromatic residue.

**Dependent Claims:**
4. The peptide of claim 1, wherein P1 is selected from Val, Ile, Leu, Ala.
5. The peptide of claim 1, wherein the C-terminal residue is Phe, Tyr, or Trp.
6. The peptide of claim 1, wherein P3 is an aromatic residue (Phe, Tyr, Trp).
7. The peptide of claim 1, comprising D-amino acids for protease resistance.
8. The peptide of claim 7, wherein all amino acids are D-configured.
9. The peptide of claim 7, in retro-inverso configuration.
10. A pharmaceutical composition comprising the peptide of claim 1 and a pharmaceutically acceptable carrier.
11. A method of treating cancer characterized by ERAP2-K392-mediated immune evasion, comprising administering the composition of claim 10.
12. A method of treating autoimmune disease associated with the ERAP2-K392 risk allele.

**Sequence Claims:**
13. A peptide selected from: VAGSAF, IAFSAF, VAWSAF, LAGSAF, IAGSAW, AAGSAF, or variants thereof having 1-2 conservative substitutions.

### Filing Details
- **Type:** US Provisional Patent Application
- **Cost:** $75 (micro entity)
- **Filing method:** USPTO EFS-Web
- **Priority date:** Filing date (12 months to file non-provisional)
- **Key evidence to include:** Boltz-2 selectivity data (4-target panel), bias check controls, structural equivalence table, PepMLM novelty evidence (zero E/D at P1)

---

## TASK ASSIGNMENT: Terminal Split (Parallelized)

### Terminal A: SASA + Local Analysis + CRO + Patent (CPU only, no GPU)
1. Generate Boltz-2 structures of Elite Four vs FULL ERAP2 (not cropped — need zinc)
2. Run SASA zinc capping test (`scripts/sasa_zinc_check.py`)
3. Analyze existing Boltz-2 CIF files for contact maps
4. Draft synthesis order for GenScript/AAPPTEC
5. Draft provisional patent application
6. Prepare assay protocol document for CRO
7. **Time:** ~3 hours
8. **Cost:** $0 (local CPU)

### Terminal B: MD — Leads vs ERAP2-K392 + IRAP (GPU instance 1)
**Peptides:** VAGSAF (Hammer) + IAFSAF (Scalpel)
**Targets:** ERAP2-K392 + IRAP
**Runs:** 4 simulations (2 peptides × 2 targets × 10ns each)
1. Launch RTX 4090 on Vast.ai
2. Install OpenMM: `conda install -c conda-forge openmm mdtraj pdbfixer`
3. Upload PDB structures + `scripts/openmm_md_protocol.py`
4. Run sequentially: VAGSAF→ERAP2, VAGSAF→IRAP, IAFSAF→ERAP2, IAFSAF→IRAP
5. Download trajectories + run `scripts/md_drift_analysis.py`
6. Destroy instance
7. **Time:** ~8 hours on RTX 4090
8. **Cost:** ~$2-3

### Terminal C: MD — Controls vs ERAP2-K392 + IRAP (GPU instance 2)
**Peptides:** VAWSAF (Trp Edge) + FASGAV (Scrambled Control)
**Targets:** ERAP2-K392 + IRAP
**Runs:** 4 simulations (2 peptides × 2 targets × 10ns each)
1. Same setup as Terminal B
2. **CRITICAL: FASGAV is the gate experiment.** If scrambled control LOCKS in ERAP2 (RMSD <2Å), the pocket is artificially sticky and ALL results from Terminals B and D are suspect. Monitor FASGAV first — if it drifts out of ERAP2 within 2-3ns, the control passes and you can trust the lead data.
3. **Time:** ~8 hours on RTX 4090
4. **Cost:** ~$2-3

### Terminal D: MD — All 4 peptides vs ERAP1 (GPU instance 3)
**Peptides:** VAGSAF, IAFSAF, VAWSAF, FASGAV
**Targets:** ERAP1 only
**Runs:** 4 simulations (4 peptides × 1 target × 10ns each)
1. Same setup as Terminal B
2. All 4 peptides should drift out of ERAP1 — confirms molecular ruler evasion in MD
3. **Time:** ~8 hours on RTX 4090
4. **Cost:** ~$2-3

### Execution Order
```
Hour 0:  Launch all 4 terminals simultaneously
Hour 1:  Terminal A delivers SASA results
Hour 2:  Terminal C — check FASGAV vs ERAP2 (the gate)
         If FASGAV RMSD > 3Å by 2ns → CONTROL PASSES, continue all
         If FASGAV RMSD < 2Å by 2ns → STOP, pocket is sticky, rethink
Hour 8:  Terminals B, C, D finish MD runs
Hour 9:  Run drift analysis on all 12 trajectories
Hour 10: Compile results, update synthesis panel
```

### Early Termination Rules
- **If FASGAV locks in ERAP2:** STOP terminals B and D. Destroy instances. The sticky pocket means Boltz-2 starting poses are biased. Need to re-equilibrate with randomized peptide placement.
- **If any lead drifts out of ERAP2-K392 within 2ns:** That lead is deprioritized. Don't waste the remaining 8ns — kill that run and reallocate GPU time.
- **If all leads drift out of IRAP within 2ns:** IRAP evasion confirmed early. Can kill IRAP runs and save GPU hours.

### Total Estimated Cost
- 3 GPU instances × 8 hrs × ~$0.27/hr = ~$6-8
- Peptide synthesis: ~$400-800
- Patent filing: $75
- Enzyme assay reagents: ~$500-800
- **Total: ~$1,000-1,700**

---

## OUTPUT DIRECTORY STRUCTURE

**Every terminal writes to its own subdirectory. No exceptions.**

```
data/results/v43_validation/
├── terminal_a/                    ← Terminal A (CPU, local)
│   ├── sasa/
│   │   ├── VAGSAF_sasa.json
│   │   ├── IAFSAF_sasa.json
│   │   ├── VAWSAF_sasa.json
│   │   ├── FASGAV_sasa.json
│   │   └── sasa_summary.json     ← zinc capping verdicts
│   ├── contact_maps/
│   │   └── *.json
│   └── terminal_a_status.json    ← updated after each task completes
│
├── terminal_b/                    ← Terminal B (GPU — leads vs ERAP2+IRAP)
│   ├── md_VAGSAF_erap2k392/
│   │   ├── trajectory.dcd
│   │   ├── rmsd.csv
│   │   ├── com_drift.csv
│   │   └── summary.json
│   ├── md_VAGSAF_irap/
│   │   ├── trajectory.dcd
│   │   ├── rmsd.csv
│   │   ├── com_drift.csv
│   │   └── summary.json
│   ├── md_IAFSAF_erap2k392/
│   │   └── (same structure)
│   ├── md_IAFSAF_irap/
│   │   └── (same structure)
│   └── terminal_b_status.json    ← updated after each run completes
│
├── terminal_c/                    ← Terminal C (GPU — controls vs ERAP2+IRAP)
│   ├── md_VAWSAF_erap2k392/
│   │   └── (same structure)
│   ├── md_VAWSAF_irap/
│   │   └── (same structure)
│   ├── md_FASGAV_erap2k392/      ← THE GATE EXPERIMENT
│   │   └── (same structure)
│   ├── md_FASGAV_irap/
│   │   └── (same structure)
│   └── terminal_c_status.json
│
├── terminal_d/                    ← Terminal D (GPU — all 4 vs ERAP1)
│   ├── md_VAGSAF_erap1/
│   │   └── (same structure)
│   ├── md_IAFSAF_erap1/
│   │   └── (same structure)
│   ├── md_VAWSAF_erap1/
│   │   └── (same structure)
│   ├── md_FASGAV_erap1/
│   │   └── (same structure)
│   └── terminal_d_status.json
│
├── compiled/                      ← Final compiled results (written AFTER all terminals finish)
│   ├── all_rmsd.csv              ← 12 runs, one row per run
│   ├── all_drift.csv
│   ├── selectivity_matrix.json   ← full 4×3 peptide×target matrix
│   ├── gate_result.json          ← FASGAV pass/fail verdict
│   └── final_verdict.md          ← human-readable summary
│
└── coordination.json              ← COORDINATION FILE (see below)
```

### Coordination File: `coordination.json`

Each terminal reads and updates this file to track progress. **Use atomic file writes (Python json.load/dump) to avoid corruption.**

```json
{
  "created": "2026-03-24T...",
  "terminal_a": {
    "status": "not_started",
    "tasks_completed": [],
    "tasks_remaining": ["sasa", "contact_maps"],
    "last_updated": null
  },
  "terminal_b": {
    "status": "not_started",
    "instance_id": null,
    "ssh": null,
    "runs_completed": [],
    "runs_remaining": ["VAGSAF_erap2k392", "VAGSAF_irap", "IAFSAF_erap2k392", "IAFSAF_irap"],
    "last_updated": null
  },
  "terminal_c": {
    "status": "not_started",
    "instance_id": null,
    "ssh": null,
    "runs_completed": [],
    "runs_remaining": ["FASGAV_erap2k392", "FASGAV_irap", "VAWSAF_erap2k392", "VAWSAF_irap"],
    "gate_result": null,
    "last_updated": null
  },
  "terminal_d": {
    "status": "not_started",
    "instance_id": null,
    "ssh": null,
    "runs_completed": [],
    "runs_remaining": ["VAGSAF_erap1", "IAFSAF_erap1", "VAWSAF_erap1", "FASGAV_erap1"],
    "last_updated": null
  },
  "gate_passed": null,
  "all_complete": false
}
```

### Rules for All Terminals

1. **Write ONLY to your own `terminal_X/` directory.** Never write to another terminal's directory.
2. **Update `coordination.json`** after completing each run. Move the run name from `runs_remaining` to `runs_completed`. Update `last_updated` timestamp.
3. **Check `coordination.json` before starting a new run.** If `gate_passed` is `false`, STOP and destroy your GPU instance.
4. **Terminal C sets the gate.** After FASGAV vs ERAP2-K392 finishes, Terminal C writes `gate_result` ("pass" or "fail") and sets `gate_passed` to `true` or `false`.
5. **Each MD run's `summary.json` must contain:**
   ```json
   {
     "peptide": "VAGSAF",
     "target": "erap2k392",
     "terminal": "b",
     "rmsd_final": 1.45,
     "rmsd_mean": 1.23,
     "rmsd_max": 2.10,
     "com_drift_final": 0.82,
     "contact_fraction_final": 0.85,
     "verdict": "LOCKED",
     "simulation_ns": 10.0,
     "wall_time_hours": 1.8,
     "gpu": "RTX 4090",
     "instance_id": "33XXXXXX"
   }
   ```
6. **Naming convention for trajectories on Vast.ai:** `/workspace/md_{PEPTIDE}_{TARGET}/` — matches the local directory name exactly. Download with `scp -r` into your terminal's directory.
7. **Do NOT run `git push` from GPU terminals.** Only push from the local machine after downloading results.

### Post-Completion: Compiling Results

After all terminals report `status: "complete"` in `coordination.json`, one terminal (preferably Terminal A since it's local) runs the compiler:

```bash
python scripts/compile_v43_results.py
```

This reads all 12 `summary.json` files, builds the selectivity matrix, and writes `compiled/final_verdict.md`.

---

## SUCCESS CRITERIA

### Computational (this phase)
- [ ] SASA shows >80% zinc capping for VAGSAF and IAFSAF
- [ ] MD shows VAGSAF/IAFSAF RMSD < 2Å in ERAP2-K392 at 10ns
- [ ] MD shows VAGSAF/IAFSAF RMSD > 5Å in IRAP (drifts out)
- [ ] MD shows FASGAV (scrambled) RMSD > 5Å in ERAP2 (not sticky pocket)
- [ ] ERAP1 simulations show no stable binding for any 6-mer

### Wet Lab (next phase)
- [ ] IC50 < 50 µM for VAGSAF or IAFSAF on ERAP2
- [ ] IC50 > 500 µM for same peptide on ERAP1
- [ ] FASGAV scrambled control shows no activity
- [ ] At least 10-fold selectivity ratio (ERAP2 IC50 / ERAP1 IC50)

### Patent (parallel)
- [ ] Provisional filed within 30 days
- [ ] Claims cover length-dependent selectivity mechanism
- [ ] Sequence claims cover Elite Four + variants
- [ ] Non-obvious arguments include steric model + molecular ruler + surface area evasion

---

## FILE LOCATIONS

| File | Path | Description |
|---|---|---|
| This plan | `docs/V43_WETLAB_VALIDATION_PLAN.md` | Master document |
| ERAP2 structure | `data/structures/erap2_wt_alphafold.pdb` | AlphaFold ERAP2 |
| ERAP2 sequence | `data/processed/erap2_Q6P179.fasta` | UniProt Q6P179 |
| Boltz-2 results | `data/results/` | All screening data |
| Short-stack results | `data/results/short_stack_results.json` | 18 peptide K392/N392 |
| Counterscreen | `data/results/counterscreen_results.json` | 18 peptide 4-target panel |
| Scaffold screen | `data/results/n392_scaffold_search/` | 120 scaffold screen |
| Structural equiv | `docs/structural_equivalence_table.md` | 6 ERAP2-unique residues |
| OpenMM script | `scripts/openmm_md_protocol.py` | TO BE WRITTEN |
| SASA script | `scripts/sasa_zinc_check.py` | TO BE WRITTEN |
| MD analysis | `scripts/md_drift_analysis.py` | TO BE WRITTEN |

---

## IMPORTANT NOTES FOR ALL TERMINALS

1. **Boltz-2 scores are single-seed for the short-stack panel.** The V4.3 campaign (160 predictions, running on instance 33430272) will provide multi-seed confirmation. Do NOT make synthesis decisions until multi-seed data is in.

2. **The cropped channel (residues 350-450) does NOT include the catalytic zinc.** For SASA zinc capping, you MUST use full ERAP2 structure and re-dock peptides into the full protein.

3. **FASGAV (scrambled control) is the most important experiment.** If it also locks in ERAP2, all our selectivity claims are suspect. Run it first or in parallel with leads.

4. **D-amino acid modifications are Phase 2.** Don't order D-peptides until L-amino acid IC50 data confirms activity. L-peptides are cheaper and faster to synthesize.

5. **The steric model (not salt bridge) should be the primary patent narrative.** Val > Glu for selectivity is counterintuitive and non-obvious. The salt bridge is a secondary contribution.
