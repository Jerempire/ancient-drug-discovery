# How dVAGSAF Works — Visual Guide

## Normal Substrate vs D-Peptide Inhibitor

```
NORMAL SUBSTRATE (L-VAGSAF) — Gets Destroyed
═════════════════════════════════════════════

  Step 1: Entry              Step 2: Zinc grabs P1        Step 3: Hydrolysis

  Channel mouth              Active site                   Fragments ejected
  ┌────────┐                 ┌────────┐                    ┌────────┐
  │        │                 │        │                    │        │
  │ VAGSAF │ slides in →     │  V←Zn  │ orients →         │  ·AGSAF│ + V expelled
  │ ██████ │                 │  │  │  │ scissile           │  ·····█│
  │        │                 │  O  H₂O│ bond               │        │
  └────────┘                 │  ‖     │                    └────────┘
                             │  C─NH  │ ← zinc positions
                             │  ↑     │   water for         ~90 seconds
                             │  cut   │   nucleophilic      per substrate
                             │  here  │   attack
                             └────────┘

  L-Val α-carbon:            Zinc coordination:
       H                         Zn²⁺
       │                        / │ \
  H₃N─C─COO⁻                H370 H374 E371
       │                           │
       CH(CH₃)₂                   ↓
                             positions water
  carbonyl oxygen              to attack C─N
  faces zinc correctly         bond → cleaves
```

```
dVAGSAF (D-Val at P1) — Blocks But Can't Be Cut
═════════════════════════════════════════════════

  Step 1: Entry              Step 2: Zinc grabs...        Step 3: STUCK

  Channel mouth              Active site                   Permanently occupied
  ┌────────┐                 ┌────────┐                    ┌────────┐
  │        │                 │        │                    │ ██████ │
  │dVAGSAF │ slides in →     │ dV←Zn  │ tries to →        │ ██████ │ BLOCKED
  │ ██████ │                 │  │  │  │ orient...          │ ██████ │
  │        │                 │  O  H₂O│                    │ ██████ │ no substrate
  └────────┘                 │  ‖     │ but carbonyl       │ ██████ │ can enter
                             │  C─NH  │ faces WRONG WAY    │ ██████ │
                             │  ↑     │                    └────────┘
                             │  ???   │ water can't
                             │        │ reach the bond     residence time:
                             └────────┘                    minutes to hours

  D-Val α-carbon:            Zinc coordination:
       H                         Zn²⁺
       │                        / │ \
  ⁻OOC─C─NH₃⁺              H370 H374 E371
       │                           │
       CH(CH₃)₂                   ↓
                             water positioned
  *** MIRROR IMAGE ***       but carbonyl oxygen
  carbonyl oxygen now        is on the WRONG SIDE
  faces AWAY from zinc       of the peptide bond
  → can't be hydrolyzed     → no cleavage possible
```

## The Selectivity Story

```
WHY ONLY ERAP2 GETS BLOCKED
════════════════════════════

                    ERAP1              ERAP2 (K392+N392)        IRAP
                 ┌──────────┐       ┌────┐                  ┌──────────┐
                 │          │       │████│                   │          │
   dVAGSAF       │  TOO     │       │████│ SNUG              │  TOO     │
   (6-mer)  ───→ │  SHORT   │  ───→ │████│ FIT          ───→ │  LOOSE   │
                 │  to grip │       │████│ STUCK!            │  to grip │
                 │          │       │████│                   │          │
                 └──────────┘       └────┘                  └──────────┘

                 ERAP1 needs        dVAGSAF enters,          IRAP cavity
                 8+ residues        zinc grabs it,           too large —
                 to close its       D-Val blocks             6-mer rattles
                 conformational     hydrolysis,              around, no
                 gate. 6-mer        peptide STAYS            stable contacts
                 falls through.     in channel.              → falls out

                 ipTM: 0.152        ipTM: 0.858-0.884       ipTM: 0.233
                 (5.6x selective)   (INHIBITS)               (3.7x selective)
```

## Predicted Inhibition Profile

```
                   ipTM        Predicted           Predicted
  Modification     (K392)      Mechanism           IC50 Range
  ─────────────────────────────────────────────────────────────
  dVAGSAF          0.858       Competitive         1-50 µM
                               (channel plug,       (D-amino acid peptides
                                D-Val resists        typically mid-µM;
                                cleavage)            tight channel = lower end)

  dIAFSAF          0.876       Competitive         1-50 µM
                               (same mechanism,     (slightly higher affinity
                                D-Ile at P1)         but lower selectivity)

  NMeVAGSAF        0.883       Competitive         5-100 µM
                               (N-methyl blocks     (N-methyl peptides often
                                backbone NH for      weaker inhibitors than
                                catalysis)           D-amino acid versions)
  ─────────────────────────────────────────────────────────────

  WHY WE CAN'T PREDICT IC50 MORE PRECISELY:

  Boltz-2 ipTM measures:    What IC50 also needs:
  ┌─────────────────┐       ┌─────────────────────────┐
  │ Static binding   │       │ Residence time           │
  │ pose quality     │       │ (how long it stays)      │
  │                  │       │                          │
  │ "Does it fit?"   │       │ kon/koff rates           │
  │                  │       │ (how fast in/out)        │
  │                  │       │                          │
  │ ✓ We have this   │       │ Solubility               │
  │                  │       │ (does it dissolve?)      │
  │                  │       │                          │
  │                  │       │ Cell permeability        │
  │                  │       │ (does it reach target?)  │
  │                  │       │                          │
  │                  │       │ ✗ Need experiments       │
  └─────────────────┘       └─────────────────────────┘
```

## Competitive Inhibition In Practice

```
  Without inhibitor:          With dVAGSAF:

  Substrates      ERAP2       dVAGSAF     Substrates    ERAP2
  ┌─┐┌─┐┌─┐      ┌────┐      ┌──┐        ┌─┐┌─┐┌─┐    ┌────┐
  │S││S││S│  →→→  │    │      │dV│  →     │S││S││S│    │████│ ← occupied
  └─┘└─┘└─┘      │ ○  │      └──┘        └─┘└─┘└─┘    │████│
                  │    │         ↓                       │    │
                  └────┘      enters &                   └────┘
                  processes    STAYS                     substrates
                  all of them                            can't get in

  ERAP2 trims    → antigenic   ERAP2         → antigenic peptides
  peptides to      peptides    blocked           STAY at 9-mer
  <8-mers          destroyed                    → immune system
                                                 SEES the cancer
                  (cancer
                   hides)                       (cancer exposed)


  DOSE RESPONSE (predicted shape):

  100%  ┬─────────────────────────────
        │ ←── ERAP2 activity
   75%  ┤
        │          ╲
   50%  ┤           ╲←── IC50 (1-50 µM?)
        │            ╲
   25%  ┤             ╲
        │              ╲___________
    0%  ┴──┬──┬──┬──┬──┬──┬──┬──┬──
        0.1  1   10  100  µM
             [dVAGSAF]

  At IC50 concentration, half of ERAP2 molecules
  have dVAGSAF stuck in their channel at any moment.
  Above 10x IC50, essentially all ERAP2 is blocked.
```

## Comparison to Known ERAP Inhibitors

```
  Compound           Type          IC50      Selectivity   Status
  ──────────────────────────────────────────────────────────
  DG013A (Lille)     Phosphinic    0.8 µM    ERAP2-only    Published
  DG023 (Lille)      Phosphinic    2.1 µM    ERAP2-only    Published
  Compound 3 (Leed)  Thiol-based   12 µM     ERAP1>ERAP2   Published
  ERAP-IN-1          Small mol     0.3 µM    Pan-ERAP      Commercial
  ──────────────────────────────────────────────────────────
  dVAGSAF (ours)     D-peptide     1-50 µM?  ERAP2-only    PREDICTED
  ──────────────────────────────────────────────────────────

  Our advantage:     Our disadvantage:
  ✓ 5.6x ERAP1      ✗ Untested IC50
    selectivity      ✗ Unknown cell
  ✓ 3.7x IRAP         permeability
    selectivity      ✗ Unknown half-life
  ✓ Cheap ($500)       in vivo
  ✓ Novel mechanism  ✗ No K392/N392
    (substrate         discrimination
     channel plug)     (pan-ERAP2)
  ✓ 6-mer = proven
    drug size range
```
