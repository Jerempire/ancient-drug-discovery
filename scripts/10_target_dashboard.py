#!/usr/bin/env python
"""
10_target_dashboard.py — Generate unified HTML dashboard for all targets.

Pulls together:
  - Evolutionary evidence (targets.yaml)
  - Cancer connections (script 07 outputs)
  - BioReason-Pro GO predictions (script 08)
  - Bayesian rankings (script 09)
  - RFdiffusion results (if available)

Output:
  data/results/target_dashboard.html
"""

import json
from datetime import datetime
from pathlib import Path

import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent
TARGETS_PATH = PROJECT_DIR / "targets.yaml"
PROC_DIR = PROJECT_DIR / "data" / "processed"
RESULTS_DIR = PROJECT_DIR / "data" / "results"
OUTPUT = RESULTS_DIR / "target_dashboard.html"


def load_targets() -> dict:
    data = yaml.safe_load(TARGETS_PATH.read_text(encoding="utf-8"))
    return data["targets"]


def load_json(path: Path) -> dict | list | None:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return None


def load_go_predictions(gene: str) -> dict | None:
    return load_json(PROC_DIR / gene / "bioreason" / "go_predictions.json")


def count_rfdiffusion_results(gene: str) -> int:
    """Count PDB files from RFdiffusion for a gene."""
    rfd_dir = RESULTS_DIR / "rfdiffusion"
    if not rfd_dir.exists():
        return 0
    return len([f for f in rfd_dir.glob(f"{gene.lower()}*.pdb")] +
               [f for f in rfd_dir.glob(f"erap2*.pdb")] if gene == "ERAP2" else
               [f for f in rfd_dir.glob(f"{gene.lower()}*.pdb")])


def build_html(targets: dict) -> str:
    """Build the full dashboard HTML."""
    # Load cross-reference and rankings
    xref = load_json(PROC_DIR / "go_cross_reference.json")
    rankings = load_json(PROC_DIR / "bayesian_rankings.json")
    evidence_matrix = load_json(PROC_DIR / "cancer_evidence_matrix.json")

    # Build rankings lookup
    rank_lookup = {}
    if rankings:
        for r in rankings:
            rank_lookup[r["target"]] = r

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Color coding
    step_colors = {
        "detection": "#4a90d9",
        "activation": "#e6a817",
        "killing": "#d94a4a",
    }

    verdict_colors = {
        "STRONG GO": "#27ae60",
        "GO": "#2ecc71",
        "LEAN GO (low confidence)": "#f39c12",
        "MAYBE": "#e67e22",
        "INSUFFICIENT EVIDENCE": "#95a5a6",
        "NO-GO": "#e74c3c",
    }

    # --- Build target cards ---
    target_cards = []
    priority = 0
    for gene in ["ERAP2", "G6PD", "NOD2", "SLC11A1", "CCR5", "DARC"]:
        cfg = targets.get(gene, {})
        if not cfg:
            continue
        priority += 1

        go_pred = load_go_predictions(gene)
        rank = rank_lookup.get(gene, {})
        ev = evidence_matrix.get(gene, {}) if evidence_matrix else {}
        step = cfg.get("immune_step", "unknown")
        step_color = step_colors.get(step, "#666")
        verdict = rank.get("verdict", "?")
        verdict_color = verdict_colors.get(verdict, "#666")
        composite = rank.get("composite_score", 0)
        confidence = rank.get("confidence", 0)

        # GO term counts
        go_mf = len(go_pred["predictions"].get("MF", [])) if go_pred else 0
        go_bp = len(go_pred["predictions"].get("BP", [])) if go_pred else 0
        go_cc = len(go_pred["predictions"].get("CC", [])) if go_pred else 0

        # Cancer associations
        n_cancer = len(ev.get("cancer_associations", []))
        n_trials = len(ev.get("trials", []))

        # Key GO terms (top 5 most specific per aspect)
        go_highlights = []
        if go_pred:
            for aspect in ["MF", "BP", "CC"]:
                terms = go_pred["predictions"].get(aspect, [])
                # Skip generic root terms
                generic = {"GO:0008150", "GO:0003674", "GO:0005575", "GO:0110165",
                           "GO:0009987", "GO:0008152", "GO:0005488", "GO:0065007",
                           "GO:0050896", "GO:0044238", "GO:0071704", "GO:0006807"}
                specific = [t for t in terms if t not in generic]
                go_highlights.extend(specific[:3])

        # Compounds
        compounds = cfg.get("known_compounds", [])
        compound_html = ""
        if compounds:
            compound_html = "<ul>" + "".join(
                f"<li><b>{c.get('name', '?')}</b> — {c.get('mechanism', '')}</li>"
                for c in compounds
            ) + "</ul>"
        else:
            compound_html = "<em>No known compounds — greenfield target</em>"

        # RFdiffusion results
        rfd_dir = RESULTS_DIR / "rfdiffusion"
        rfd_count = 0
        if rfd_dir.exists():
            rfd_count = len(list(rfd_dir.glob(f"{gene.lower()}*")) +
                          (list(rfd_dir.glob("erap2*")) if gene == "ERAP2" else []))

        # Dimension breakdown
        dims_html = ""
        if rank.get("dimensions"):
            dims_html = "<table class='dim-table'>"
            for dim, est in rank["dimensions"].items():
                bar_width = int(est["mean"] * 100)
                dims_html += (
                    f"<tr><td>{dim.replace('_', ' ').title()}</td>"
                    f"<td><div class='bar-bg'><div class='bar-fill' style='width:{bar_width}%'></div></div></td>"
                    f"<td>{est['mean']:.0%}</td></tr>"
                )
            dims_html += "</table>"

        card = f"""
        <div class="target-card">
            <div class="card-header" style="border-left: 5px solid {step_color}">
                <div class="card-title">
                    <span class="priority">#{priority}</span>
                    <h2>{gene}</h2>
                    <span class="step-badge" style="background:{step_color}">{step}</span>
                </div>
                <div class="verdict" style="color:{verdict_color}">{verdict}</div>
            </div>
            <div class="card-subtitle">{cfg.get('protein_name', '')}</div>

            <div class="card-grid">
                <div class="card-section">
                    <h3>Evolutionary Signal</h3>
                    <p><b>Pathogen:</b> {cfg.get('ancient_pathogen', '?')}</p>
                    <p><b>Era:</b> {cfg.get('selection_era', '?')}</p>
                    <p><b>Selection:</b> {cfg.get('selection_event', '?')[:100]}</p>
                    <p><b>Key variants:</b> {len(cfg.get('key_variants', []))}</p>
                </div>

                <div class="card-section">
                    <h3>Cancer Connection</h3>
                    <p><b>Types:</b> {', '.join(cfg.get('cancer_types', []))}</p>
                    <p><b>OpenTargets associations:</b> {n_cancer}</p>
                    <p><b>Clinical trials:</b> {n_trials}</p>
                    <p><b>Drug status:</b> {cfg.get('drug_status', '?')[:80]}</p>
                </div>

                <div class="card-section">
                    <h3>BioReason-Pro GO Terms</h3>
                    <p><b>MF:</b> {go_mf} | <b>BP:</b> {go_bp} | <b>CC:</b> {go_cc} | <b>Total:</b> {go_mf+go_bp+go_cc}</p>
                    <p class="go-terms">{' '.join(go_highlights[:6])}</p>
                </div>

                <div class="card-section">
                    <h3>Compounds</h3>
                    {compound_html}
                </div>
            </div>

            <div class="card-section">
                <h3>Bayesian Score: {composite:.0%} (confidence: {confidence:.0%})</h3>
                {dims_html}
            </div>

            <div class="card-section">
                <h3>Pipeline Status</h3>
                <div class="pipeline">
                    <span class="stage done">Data Collection</span>
                    <span class="arrow">&rarr;</span>
                    <span class="stage done">Cancer Evidence</span>
                    <span class="arrow">&rarr;</span>
                    <span class="stage done">GO-GPT</span>
                    <span class="arrow">&rarr;</span>
                    <span class="stage {'done' if gene == 'ERAP2' else 'todo'}">RFdiffusion</span>
                    <span class="arrow">&rarr;</span>
                    <span class="stage todo">Validation</span>
                    <span class="arrow">&rarr;</span>
                    <span class="stage todo">Patent</span>
                </div>
            </div>
        </div>
        """
        target_cards.append(card)

    # --- Pathway connections ---
    connections_html = ""
    if xref and xref.get("pathway_connections"):
        connections_html = "<div class='connections'><h2>Cross-Target Pathway Connections</h2>"
        for conn in xref["pathway_connections"]:
            a, b = conn["target_a"], conn["target_b"]
            n = conn["n_specific_overlaps"]
            cancers = conn.get("shared_cancer_types", [])
            connections_html += f"""
            <div class="connection-card">
                <b>{a}</b> ({conn['step_a']}) &harr; <b>{b}</b> ({conn['step_b']})
                — {n} shared GO terms
                {f'<br><small>Shared cancers: {", ".join(cancers)}</small>' if cancers else ''}
            </div>
            """
        connections_html += "</div>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Ancient Drug Discovery — Target Dashboard</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           background: #0d1117; color: #c9d1d9; padding: 20px; }}
    h1 {{ color: #58a6ff; margin-bottom: 5px; }}
    .subtitle {{ color: #8b949e; margin-bottom: 20px; }}
    .target-card {{ background: #161b22; border: 1px solid #30363d; border-radius: 8px;
                    padding: 20px; margin-bottom: 20px; }}
    .card-header {{ display: flex; justify-content: space-between; align-items: center;
                    padding-left: 10px; margin-bottom: 8px; }}
    .card-title {{ display: flex; align-items: center; gap: 12px; }}
    .priority {{ color: #8b949e; font-size: 1.2em; }}
    h2 {{ color: #f0f6fc; font-size: 1.5em; }}
    .step-badge {{ color: white; padding: 2px 10px; border-radius: 12px; font-size: 0.8em; }}
    .verdict {{ font-size: 1.1em; font-weight: bold; }}
    .card-subtitle {{ color: #8b949e; margin-bottom: 15px; padding-left: 15px; }}
    .card-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 15px; }}
    .card-section {{ background: #0d1117; padding: 12px; border-radius: 6px; }}
    .card-section h3 {{ color: #58a6ff; margin-bottom: 8px; font-size: 0.9em; }}
    .card-section p {{ margin-bottom: 4px; font-size: 0.85em; }}
    .go-terms {{ color: #7ee787; font-family: monospace; font-size: 0.75em; word-break: break-all; }}
    .dim-table {{ width: 100%; font-size: 0.8em; }}
    .dim-table td {{ padding: 3px 8px; }}
    .dim-table td:first-child {{ width: 35%; color: #8b949e; }}
    .dim-table td:last-child {{ width: 40px; text-align: right; }}
    .bar-bg {{ background: #21262d; border-radius: 4px; height: 12px; width: 100%; }}
    .bar-fill {{ background: #58a6ff; border-radius: 4px; height: 12px; min-width: 2px; }}
    .pipeline {{ display: flex; align-items: center; gap: 5px; flex-wrap: wrap; }}
    .stage {{ padding: 4px 10px; border-radius: 4px; font-size: 0.8em; }}
    .stage.done {{ background: #238636; color: white; }}
    .stage.todo {{ background: #21262d; color: #8b949e; }}
    .arrow {{ color: #484f58; }}
    .connections {{ margin-top: 20px; }}
    .connection-card {{ background: #161b22; border: 1px solid #30363d; border-radius: 6px;
                        padding: 10px 15px; margin-bottom: 8px; font-size: 0.9em; }}
    ul {{ margin-left: 16px; font-size: 0.85em; }}
    li {{ margin-bottom: 3px; }}
    .summary-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; margin-bottom: 25px; }}
    .summary-card {{ background: #161b22; border: 1px solid #30363d; border-radius: 8px; padding: 15px; text-align: center; }}
    .summary-card .big {{ font-size: 2em; color: #58a6ff; }}
    .summary-card .label {{ color: #8b949e; font-size: 0.85em; }}
</style>
</head>
<body>

<h1>Ancient Drug Discovery — Target Dashboard</h1>
<p class="subtitle">Mining evolutionary immune selection for undrugged cancer targets | Generated {now}</p>

<div class="summary-grid">
    <div class="summary-card">
        <div class="big">{len(targets)}</div>
        <div class="label">Targets</div>
    </div>
    <div class="summary-card">
        <div class="big">{sum(len(t.get('cancer_types',[])) for t in targets.values())}</div>
        <div class="label">Cancer Type Links</div>
    </div>
    <div class="summary-card">
        <div class="big">{sum(1 for t in targets.values() if 'no drugs' in t.get('drug_status','').lower() or 'undrugged' in t.get('drug_status','').lower())}</div>
        <div class="label">Undrugged Targets</div>
    </div>
</div>

{''.join(target_cards)}

{connections_html}

<div style="margin-top: 30px; color: #484f58; font-size: 0.8em; text-align: center;">
    Ancient Drug Discovery Pipeline | BioReason-Pro + Bayesian Scoring + RFdiffusion<br>
    Core thesis: Genes shaped by 10,000 years of plague, malaria, TB, and smallpox are the same genes cancers hijack.
</div>

</body>
</html>"""
    return html


def main():
    print("Building target dashboard...")
    targets = load_targets()
    html = build_html(targets)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(html, encoding="utf-8")
    print(f"Dashboard saved: {OUTPUT}")


if __name__ == "__main__":
    main()
