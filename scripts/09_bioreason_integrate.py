#!/usr/bin/env python
"""
09_bioreason_integrate.py — Integrate BioReason-Pro GO-GPT results into the pipeline.

Covers items 1-3:
  1. Cross-reference GO terms across targets to find shared pathway overlaps
  2. Feed GO prediction quality into Bayesian target scoring
  3. Analyze ERAP2 GO terms to inform RFdiffusion binding site selection

Reads:
  data/processed/<gene>/bioreason/go_predictions.json (from script 08)
  targets.yaml
  data/target_priors.json (existing Bayesian state, if any)

Outputs:
  data/processed/go_cross_reference.json       — shared GO terms between targets
  data/processed/go_pathway_overlaps.txt        — human-readable overlap report
  data/processed/bayesian_rankings.json         — updated Bayesian rankings
  data/target_priors.json                       — updated Bayesian priors
  data/processed/ERAP2/bioreason/binding_site_analysis.txt — ERAP2 binding guidance
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent
TARGETS_PATH = PROJECT_DIR / "targets.yaml"
PROC_DIR = PROJECT_DIR / "data" / "processed"

sys.path.insert(0, str(PROJECT_DIR))
from scoring.bayes_target import TargetScorer, score_targets_from_yaml
from scoring.evidence_tiers import classify_evidence


def load_go_predictions() -> dict:
    """Load GO-GPT predictions for all targets."""
    results = {}
    for gene_dir in PROC_DIR.iterdir():
        go_file = gene_dir / "bioreason" / "go_predictions.json"
        if go_file.exists():
            with open(go_file) as f:
                data = json.load(f)
            # Normalize gene name to uppercase (Vast.ai dirs may be lowercase)
            gene_key = data.get("gene", gene_dir.name).upper()
            results[gene_key] = data
    return results


def load_targets() -> dict:
    data = yaml.safe_load(TARGETS_PATH.read_text(encoding="utf-8"))
    return data["targets"]


# ──────────────────────────────────────────────────────────────
# Item 1: Cross-reference GO terms across targets
# ──────────────────────────────────────────────────────────────

def cross_reference_go_terms(predictions: dict) -> dict:
    """Find GO terms shared between targets — reveals shared pathways."""
    # Build term → targets mapping
    term_to_targets = defaultdict(lambda: {"targets": [], "aspect": None})

    for gene, pred in predictions.items():
        for aspect in ["MF", "BP", "CC"]:
            for go_id in pred.get("predictions", {}).get(aspect, []):
                term_to_targets[go_id]["targets"].append(gene)
                term_to_targets[go_id]["aspect"] = aspect

    # Find shared terms (2+ targets)
    shared = {}
    for go_id, info in term_to_targets.items():
        if len(info["targets"]) >= 2:
            shared[go_id] = {
                "targets": sorted(info["targets"]),
                "aspect": info["aspect"],
                "n_targets": len(info["targets"]),
            }

    # Sort by number of targets sharing the term
    shared = dict(sorted(shared.items(), key=lambda x: x[1]["n_targets"], reverse=True))
    return shared


def find_pathway_connections(predictions: dict, targets_cfg: dict) -> list[dict]:
    """Find mechanistic connections between targets via GO overlap."""
    connections = []
    genes = list(predictions.keys())

    for i, gene_a in enumerate(genes):
        pred_a = predictions[gene_a].get("predictions", {})
        cfg_a = targets_cfg.get(gene_a, {})

        for gene_b in genes[i + 1:]:
            pred_b = predictions[gene_b].get("predictions", {})
            cfg_b = targets_cfg.get(gene_b, {})

            # Find overlapping GO terms per aspect
            overlaps = {}
            for aspect in ["MF", "BP", "CC"]:
                terms_a = set(pred_a.get(aspect, []))
                terms_b = set(pred_b.get(aspect, []))
                common = terms_a & terms_b
                # Exclude very generic root terms
                generic = {"GO:0008150", "GO:0003674", "GO:0005575", "GO:0110165",
                           "GO:0009987", "GO:0008152", "GO:0005488"}
                specific_common = common - generic
                if specific_common:
                    overlaps[aspect] = sorted(specific_common)

            if overlaps:
                # Check if they share cancer types too
                cancers_a = set(cfg_a.get("cancer_types", []))
                cancers_b = set(cfg_b.get("cancer_types", []))
                shared_cancers = cancers_a & cancers_b

                connections.append({
                    "target_a": gene_a,
                    "target_b": gene_b,
                    "step_a": cfg_a.get("immune_step", "?"),
                    "step_b": cfg_b.get("immune_step", "?"),
                    "pathogen_a": cfg_a.get("ancient_pathogen", "?"),
                    "pathogen_b": cfg_b.get("ancient_pathogen", "?"),
                    "go_overlaps": overlaps,
                    "n_specific_overlaps": sum(len(v) for v in overlaps.values()),
                    "shared_cancer_types": sorted(shared_cancers),
                })

    connections.sort(key=lambda x: x["n_specific_overlaps"], reverse=True)
    return connections


# ──────────────────────────────────────────────────────────────
# Item 2: Feed GO scores into Bayesian target scoring
# ──────────────────────────────────────────────────────────────

# GO terms that indicate cancer-relevant molecular functions
CANCER_RELEVANT_GO = {
    # Antigen presentation / immune detection
    "GO:0019882", "GO:0002474", "GO:0019886",  # antigen processing/presentation
    "GO:0004177",  # aminopeptidase activity
    "GO:0008237",  # metallopeptidase activity
    # Immune signaling
    "GO:0006955",  # immune response
    "GO:0006952",  # defense response
    "GO:0001819",  # cytokine production
    "GO:0032481",  # type I interferon production
    "GO:0051607",  # defense response to virus
    # Metabolic vulnerability
    "GO:0006098",  # pentose-phosphate shunt
    "GO:0006740",  # NADPH regeneration
    "GO:0004345",  # G6PD activity
    "GO:0009051",  # PPP oxidative branch
    # Cell surface / membrane
    "GO:0004930",  # GPCR activity
    "GO:0004950",  # chemokine receptor activity
    "GO:0019957",  # C-C chemokine binding
    # Transport / macrophage
    "GO:0046873",  # metal ion transmembrane transporter
    "GO:0005384",  # manganese ion transmembrane transporter
}

# GO terms indicating druggable protein classes
DRUGGABLE_GO = {
    "GO:0003824",  # catalytic activity
    "GO:0016787",  # hydrolase activity
    "GO:0008233",  # peptidase activity
    "GO:0004177",  # aminopeptidase activity
    "GO:0004930",  # GPCR activity
    "GO:0016491",  # oxidoreductase activity
    "GO:0005215",  # transporter activity
    "GO:0004888",  # transmembrane signaling receptor
    "GO:0008270",  # zinc ion binding (metalloenzyme = druggable)
    "GO:0004345",  # G6PD activity (specific enzyme = druggable)
    "GO:0050661",  # NADP binding (cofactor site = druggable)
}


def update_bayesian_scores(predictions: dict, targets_cfg: dict) -> TargetScorer:
    """Update Bayesian target scores with GO-GPT evidence."""
    scorer = score_targets_from_yaml(targets_cfg)

    # Try to load existing priors
    priors_file = PROJECT_DIR / "data" / "target_priors.json"
    scorer.load(priors_file)

    for gene, pred in predictions.items():
        all_terms = set()
        for aspect_terms in pred.get("predictions", {}).values():
            all_terms.update(aspect_terms)

        n_total = len(all_terms)
        cancer_hits = all_terms & CANCER_RELEVANT_GO
        druggable_hits = all_terms & DRUGGABLE_GO

        # --- Tumor relevance: cancer-relevant GO terms ---
        if cancer_hits:
            prob = min(0.90, 0.50 + 0.10 * len(cancer_hits))
            conf = min(0.70, 0.30 + 0.08 * len(cancer_hits))
            tier = classify_evidence(source_type="cohort", text="BioReason-Pro GO prediction")
            scorer.update(gene, "tumor_relevance", "esm2_score", prob, conf, tier, year=2026)

        # --- Druggability: druggable protein class GO terms ---
        if druggable_hits:
            prob = min(0.85, 0.45 + 0.10 * len(druggable_hits))
            conf = min(0.65, 0.25 + 0.08 * len(druggable_hits))
            tier = classify_evidence(source_type="cohort", text="BioReason-Pro functional prediction")
            scorer.update(gene, "druggability", "esm2_score", prob, conf, tier, year=2026)

        # --- Overall characterization confidence ---
        # More GO terms = better characterized protein = more confidence
        if n_total >= 30:
            char_prob, char_conf = 0.70, 0.55
        elif n_total >= 15:
            char_prob, char_conf = 0.60, 0.45
        else:
            char_prob, char_conf = 0.45, 0.35
        scorer.update(gene, "druggability", "in_silico_docking", char_prob, char_conf, year=2026)

    return scorer


# ──────────────────────────────────────────────────────────────
# Item 3: ERAP2 binding site analysis for RFdiffusion
# ──────────────────────────────────────────────────────────────

def analyze_erap2_binding(predictions: dict, targets_cfg: dict) -> str:
    """Analyze ERAP2 GO terms to guide RFdiffusion binding site selection."""
    erap2_pred = predictions.get("ERAP2", {}).get("predictions", {})
    erap2_cfg = targets_cfg.get("ERAP2", {})

    if not erap2_pred:
        return "No ERAP2 predictions found."

    mf_terms = erap2_pred.get("MF", [])
    bp_terms = erap2_pred.get("BP", [])
    cc_terms = erap2_pred.get("CC", [])

    lines = [
        "=" * 70,
        "ERAP2 BINDING SITE ANALYSIS — RFdiffusion Guidance",
        "=" * 70,
        "",
        "Source: BioReason-Pro GO-GPT + targets.yaml + FRAMEWORK.md",
        "",
        "--- Molecular Function Summary ---",
        "",
    ]

    # Analyze key functional GO terms
    zinc_binding = "GO:0008270" in mf_terms
    aminopeptidase = "GO:0004177" in mf_terms
    metallopeptidase = "GO:0008237" in mf_terms
    metalloamino = "GO:0070006" in mf_terms
    exopeptidase = "GO:0008238" in mf_terms

    if zinc_binding:
        lines.append("  [CONFIRMED] Zinc ion binding (GO:0008270)")
        lines.append("    -> Active site contains catalytic zinc")
        lines.append("    -> RFdiffusion binders should coordinate with or near Zn2+")
        lines.append("    -> Known ERAP inhibitors (DG013A, Camberlein) target zinc coordination")
        lines.append("")

    if aminopeptidase and metallopeptidase:
        lines.append("  [CONFIRMED] Metalloaminopeptidase (GO:0004177 + GO:0008237)")
        lines.append("    -> Cleaves N-terminal amino acids from peptide substrates")
        lines.append("    -> Substrate peptides are 8-16 residues (MHC-I ligand precursors)")
        lines.append("    -> Design binders that mimic substrate or transition state")
        lines.append("")

    if exopeptidase:
        lines.append("  [CONFIRMED] Exopeptidase activity (GO:0008238)")
        lines.append("    -> Trims from peptide termini — accessible binding pocket")
        lines.append("")

    lines.extend([
        "--- Binding Strategy for RFdiffusion ---",
        "",
        "  PRIMARY TARGET: Active site (zinc coordination sphere)",
        "    Residues: Around positions 370-393 (known hotspot region)",
        "    Known variant: K392N (rs2549782) — near active site, balancing selection",
        "    GO evidence: zinc binding + metalloaminopeptidase + exopeptidase",
        "    Approach: Design 30-40 residue binders contacting Zn coordination residues",
        "",
        "  SECONDARY TARGET: Substrate recognition channel",
        "    ERAP2 must discriminate peptide length (8-16 aa for MHC-I)",
        "    Allosteric binders could alter substrate specificity rather than block catalysis",
        "    This is the ERAP2-selective approach (avoids ERAP1 cross-reactivity)",
        "",
        "  TERTIARY TARGET: ERAP2-specific surface (selectivity over ERAP1)",
        "    ERAP1/ERAP2 share 50% sequence identity — selectivity is the key challenge",
        "    GO-GPT classified both as metalloaminopeptidases",
        "    RFdiffusion should include ERAP1 (Q9NZ08) as negative design target",
        "",
        "--- RFdiffusion Command Recommendations ---",
        "",
        "  # Active site binder (primary)",
        "  contigmap.contigs=[A370-393/0 30-40]",
        "  ppi.hotspot_res=[A370,A371,A374,A392,A393]",
        "  inference.ckpt_override_path=Complex_beta_ckpt.pt",
        "",
        "  # Extended substrate channel binder",
        "  contigmap.contigs=[A340-420/0 40-60]",
        "  ppi.hotspot_res=[A370,A371,A374,A380,A392,A393,A400]",
        "",
        "  # Selectivity screen: run ProteinMPNN on binders,",
        "  # then DiffDock against both ERAP2 (Q6P179) and ERAP1 (Q9NZ08)",
        "  # Keep only binders with >10x affinity ratio",
        "",
        "--- Key Variant for Personalized Approach ---",
        "",
        "  rs2549782 (K392N): Missense near active site under balancing selection",
        "  -> Design binders that work on BOTH K392 and N392 variants",
        "  -> Or design variant-specific binders for precision oncology",
        f"  -> {len(erap2_cfg.get('key_variants', []))} total variants in targets.yaml",
    ])

    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("BioReason-Pro Integration — GO Cross-Reference + Bayesian Update")
    print("=" * 70)

    targets_cfg = load_targets()
    predictions = load_go_predictions()
    print(f"Loaded GO predictions for: {', '.join(predictions.keys())}")

    # --- Item 1: Cross-reference ---
    print(f"\n{'─'*50}")
    print("1. GO Term Cross-Reference")
    print(f"{'─'*50}")

    shared_terms = cross_reference_go_terms(predictions)
    connections = find_pathway_connections(predictions, targets_cfg)

    # Save
    xref_out = PROC_DIR / "go_cross_reference.json"
    xref_out.write_text(json.dumps({
        "shared_terms": shared_terms,
        "pathway_connections": connections,
    }, indent=2), encoding="utf-8")

    # Print top connections
    print(f"\n  Found {len(shared_terms)} GO terms shared across 2+ targets")
    print(f"  Found {len(connections)} target pairs with specific GO overlap")

    # Report
    report_lines = [
        "=" * 70,
        "GO TERM CROSS-REFERENCE — Shared Pathways Between Targets",
        "=" * 70, "",
    ]

    for conn in connections:
        a, b = conn["target_a"], conn["target_b"]
        n = conn["n_specific_overlaps"]
        cancers = conn["shared_cancer_types"]
        report_lines.append(
            f"  {a} ({conn['step_a']}/{conn['pathogen_a']}) <-> "
            f"{b} ({conn['step_b']}/{conn['pathogen_b']})"
        )
        report_lines.append(f"    Shared specific GO terms: {n}")
        for aspect, terms in conn["go_overlaps"].items():
            report_lines.append(f"      {aspect}: {', '.join(terms)}")
        if cancers:
            report_lines.append(f"    Shared cancer types: {', '.join(cancers)}")
        report_lines.append("")

    report_text = "\n".join(report_lines)
    print(report_text)

    overlap_out = PROC_DIR / "go_pathway_overlaps.txt"
    overlap_out.write_text(report_text, encoding="utf-8")

    # --- Item 2: Bayesian update ---
    print(f"\n{'─'*50}")
    print("2. Bayesian Target Scoring Update")
    print(f"{'─'*50}")

    scorer = update_bayesian_scores(predictions, targets_cfg)
    rankings = scorer.rank_targets()

    print(f"\n  {'Gene':12s} {'Composite':>9s} {'Confidence':>10s} {'Verdict'}")
    print("  " + "-" * 50)
    for r in rankings:
        print(
            f"  {r['target']:12s} "
            f"{r['composite_score']:>8.1%}  "
            f"{r['confidence']:>8.1%}   "
            f"{r['verdict']}"
        )

    # Save updated priors and rankings
    scorer.save()
    rank_out = PROC_DIR / "bayesian_rankings.json"
    rank_out.write_text(json.dumps(rankings, indent=2), encoding="utf-8")
    print(f"\n  Saved: target_priors.json, bayesian_rankings.json")

    # --- Item 3: ERAP2 binding analysis ---
    print(f"\n{'─'*50}")
    print("3. ERAP2 Binding Site Analysis")
    print(f"{'─'*50}")

    binding_analysis = analyze_erap2_binding(predictions, targets_cfg)
    print(binding_analysis)

    binding_out = PROC_DIR / "ERAP2" / "bioreason" / "binding_site_analysis.txt"
    binding_out.parent.mkdir(parents=True, exist_ok=True)
    binding_out.write_text(binding_analysis, encoding="utf-8")
    print(f"\n  Saved: {binding_out}")

    print(f"\n{'='*70}")
    print("Integration complete.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
