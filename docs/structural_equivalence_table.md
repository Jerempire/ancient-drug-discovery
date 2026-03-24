# ERAP2 / ERAP1 / IRAP Structural Residue Equivalence

**Date:** 2026-03-23
**Method:** BioPython structural superposition (Superimposer) on PDB crystal structures
**Validated by:** UniProt sequence cross-check, HEXXH anchor offset analysis, local domain alignment, literature review

## Structures Used

| Enzyme | PDB ID | Chain | Residues | Notes |
|--------|--------|-------|----------|-------|
| ERAP2 | 7SH0 | A | 896 | N392 allele (K392 is alternate) |
| ERAP1 | 6RQX | A | 862 | |
| IRAP | 4PJ6 | A | 861 | |

## Alignment Quality

| Alignment | Global RMSD | Local (catalytic domain) RMSD | Paired CAs |
|-----------|-------------|-------------------------------|------------|
| ERAP1 -> ERAP2 | 4.44 A | 2.71 A | 499 / 109 |
| IRAP -> ERAP2 | 4.70 A | 3.79 A | 464 / 94 |

## Validated Equivalence Table

| ERAP2 | ERAP1 equiv | IRAP equiv | HEXXH offset (E2/E1/IR) | CA-CA dist (local) | Chemistry |
|-------|-------------|------------|--------------------------|---------------------|-----------|
| **K392** (Lys, +) | N375 (Asn, 0) | N486 (Asn, 0) | +22/+21/+22 | 0.3 / 1.3 A | **ERAP2-unique positive charge** |
| **Y398** (Tyr, OH) | F381 (Phe) | F492 (Phe) | +28/+27/+28 | 0.9 / 0.6 A | **ERAP2-unique H-bond donor** |
| **A403** (Ala, small) | S386 (Ser) | S497 (Ser) | +33/+32/+33 | 1.0 / 0.8 A | ERAP2 hydrophobic vs polar |
| **A406** (Ala, small) | V389 (Val) | K500 (Lys, +) | +36/+35/+36 | 1.4 / 0.8 A | **IRAP-unique positive charge** |
| **Q412** (Gln, 0) | K395 (Lys, +) | S506 (Ser, 0) | +42/+41/+42 | 2.0 / 1.0 A | **ERAP1-unique positive charge** |
| **D414** (Asp, -) | G397 (Gly) | Y508 (Tyr) | +44/+43/+44 | 1.5 / 1.6 A | **ERAP2-unique negative charge** |

Notes:
- CA-CA distances are post local alignment (catalytic domain only, tighter RMSD)
- Positions 398-414 fall in a loop with indels — sequence alignment shows GAPs but spatial alignment is consistent
- ERAP1 is systematically -1 offset from ERAP2/IRAP due to a single indel between HEXXH and this loop
- Position 392 in 7SH0 is crystallized as N (Asn); UniProt canonical Q6P179 is K (Lys) — rs2549782 K/N polymorphism

## Selectivity Assessment

### Best positions for ERAP2-selective peptide design (ranked):

1. **D414** — Only negative charge at this position. Peptide Arg/Lys here forms salt bridge exclusive to ERAP2 (ERAP1 has Gly, IRAP has Tyr — neither can reciprocate).

2. **K392** (K allele) — Only positive charge. Peptide Glu/Asp forms salt bridge exclusive to ERAP2 (ERAP1=N375, IRAP=N486, both neutral Asn). Confirmed by Evnouchidou 2012 (PMID 22837489): K392 directly affects S1 pocket specificity and peptide N-terminus stabilization.

3. **Y398** — Tyr hydroxyl enables H-bonds that ERAP1 F381 and IRAP F492 cannot form (Phe lacks OH). Design peptides with H-bond acceptors (Asn, Gln, Ser) at this contact.

4. **A406** — Small hydrophobic pocket in ERAP2. IRAP has K500 (bulky, positive) — steric + charge clash provides built-in IRAP counter-selection. ERAP1 has V389 (similar size, no charge).

5. **A403** — Hydrophobic in ERAP2 vs polar Ser in both ERAP1 and IRAP. Moderate selectivity via hydrophobic contacts.

6. **Q412** — ERAP1 has K395 (positive charge). Avoid placing negative charge at this peptide position or it will cross-react with ERAP1. Keep neutral.

### Key answers:

- **Does ERAP2 K392 have a unique positive charge?** YES (for K allele). ERAP1=Asn, IRAP=Asn at structurally equivalent positions.
- **Would Glu/Asp form a selective salt bridge?** YES, exclusively with ERAP2 K392. Neither ERAP1 N375 nor IRAP N486 can form a salt bridge.
- **D414 is also ERAP2-unique** — a negative charge not present in ERAP1 (Gly) or IRAP (Tyr). Peptide Arg/Lys at this contact is a second orthogonal selectivity handle.

## Validation Summary

| Check | Result |
|-------|--------|
| UniProt sequence identity (18 residues) | 18/18 match |
| Local domain alignment stability | All 6 equivalences held |
| HEXXH anchor offset consistency | ERAP2/IRAP identical; ERAP1 systematic -1 (known indel) |
| Literature (Evnouchidou 2012, Zervoudi 2011, Mpakali 2022) | Partial confirmation; no prior residue correspondence table for 398-414 |

## Literature References

- Evnouchidou et al. (2012) PMID 22837489 — K392N specificity switch, S1 pocket
- Zervoudi et al. (2011) PMID 21314638 — S1 pocket non-conserved residues across ERAP1/ERAP2/IRAP
- Mpakali et al. (2015) PMID 26381406 — ERAP2 crystal structure, exon10 loop divergence
- Mpakali et al. (2022) PMID 35178178 — ERAP2 S1' pocket, ERAP1/IRAP superposition (Fig. 4)
- Birtley et al. (2012) PMID 22106953 — ERAP2 3.08A structure (PDB 3SE6)

## Reproduction

Script: `scripts/structural_alignment_analysis.py`
Local alignment check: ran as inline script (BioPython Superimposer on residues 340-520)

### PyMOL commands:

```
fetch 7SH0, erap2
fetch 6RQX, erap1
fetch 4PJ6, irap

align erap1, erap2
align irap, erap2

select erap2_targets, erap2 and chain A and resi 392+398+403+406+412+414
show sticks, erap2_targets
color yellow, erap2_targets

# K392 equivalents
select e2_392, erap2 and chain A and resi 392
select e1_375, erap1 and chain A and resi 375
select ir_486, irap and chain A and resi 486
distance dist_e1_392, e2_392 and name CA, e1_375 and name CA
distance dist_ir_392, e2_392 and name CA, ir_486 and name CA

# Y398 equivalents
select e2_398, erap2 and chain A and resi 398
select e1_381, erap1 and chain A and resi 381
select ir_492, irap and chain A and resi 492
distance dist_e1_398, e2_398 and name CA, e1_381 and name CA
distance dist_ir_398, e2_398 and name CA, ir_492 and name CA

# A403 equivalents
select e2_403, erap2 and chain A and resi 403
select e1_386, erap1 and chain A and resi 386
select ir_497, irap and chain A and resi 497

# A406 equivalents
select e2_406, erap2 and chain A and resi 406
select e1_389, erap1 and chain A and resi 389
select ir_500, irap and chain A and resi 500

# Q412 equivalents
select e2_412, erap2 and chain A and resi 412
select e1_395, erap1 and chain A and resi 395
select ir_506, irap and chain A and resi 506

# D414 equivalents
select e2_414, erap2 and chain A and resi 414
select e1_397, erap1 and chain A and resi 397
select ir_508, irap and chain A and resi 508

show sticks, e1_375 or ir_486 or e1_381 or ir_492 or e1_386 or ir_497 or e1_389 or ir_500 or e1_395 or ir_506 or e1_397 or ir_508

color cyan, erap2
color green, erap1
color salmon, irap
color yellow, erap2_targets
zoom erap2_targets, 12
```
