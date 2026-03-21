"""
generate_dashboard.py — Build interactive SeqDesign-style binder explorer dashboard.

Auto-discovers validation_summary*.json files in data/results/, merges with
candidate_ranking.json and selectivity data, embeds PDB structures inline,
and outputs a self-contained index.html.

Usage:
    cd ancient-drug-discovery
    python dashboard/generate_dashboard.py
"""

import glob
import json
import os
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "data" / "results"
OUTPUT_HTML = Path(__file__).resolve().parent / "index.html"

# Hotspot residues (ERAP2 numbering)
HOTSPOT_RESIDUES = [353, 355, 360, 363, 367, 401, 403, 406, 408, 412]


def detect_version(json_path: str) -> str:
    """Infer version from path (v2 if 'v2' in path, else v1)."""
    if "v2" in json_path.lower():
        return "v2"
    return "v1"


def detect_family(design_name: str) -> str:
    for fam in ("short", "medium", "long", "large"):
        if fam in design_name:
            return fam
    return "unknown"


def load_validation_summaries() -> list[dict]:
    """Find all validation_summary*.json and merge into unified list."""
    candidates = []
    seen_designs = set()

    # Find all validation summary files (skip tier1 which has different format)
    patterns = [
        str(RESULTS_DIR / "boltz2_validation_v2" / "validation_summary_v2.json"),
        str(RESULTS_DIR / "boltz2_validation" / "validation_summary.json"),
    ]

    # Also glob for any future ones
    for p in glob.glob(str(RESULTS_DIR / "**" / "validation_summary*.json"), recursive=True):
        if "tier1" not in p and p not in patterns:
            patterns.append(p)

    for path in patterns:
        if not os.path.exists(path):
            print(f"  Skipping (not found): {path}")
            continue

        version = detect_version(path)
        with open(path, encoding="utf-8") as f:
            data = json.load(f)

        if isinstance(data, dict):
            # tier1 format or other dict format — skip
            continue

        for entry in data:
            design = entry.get("design", "")
            if design in seen_designs:
                continue
            seen_designs.add(design)

            seq = entry.get("binder_sequence", entry.get("sequence", ""))
            e2_scores = entry.get("erap2_scores", {})
            e1_scores = entry.get("erap1_scores", {})

            # Extract binder pLDDT from pair_chains_iptm if available
            binder_plddt = None
            pair = e2_scores.get("pair_chains_iptm", {})
            if "1" in pair and "1" in pair["1"]:
                binder_plddt = pair["1"]["1"]

            candidates.append({
                "design_id": design,
                "version": version,
                "family": detect_family(design),
                "sequence": seq,
                "binder_length": entry.get("binder_length", len(seq)),
                "iptm_erap2": e2_scores.get("iptm", 0),
                "iptm_erap1": e1_scores.get("iptm", 0),
                "iptm_delta": entry.get("iptm_delta", 0),
                "iptm_selectivity": entry.get("iptm_selectivity", 0),
                "complex_plddt": e2_scores.get("complex_plddt", 0),
                "ptm": e2_scores.get("ptm", 0),
                "mpnn_score": entry.get("mpnn_score", 0),
                "binder_plddt": binder_plddt,
            })

        print(f"  Loaded {len(data)} candidates from {path} (version={version})")

    return candidates


def enrich_with_ranking(candidates: list[dict]) -> list[dict]:
    """Merge candidate_ranking.json data (scores, tiers, developability)."""
    ranking_path = RESULTS_DIR / "candidate_ranking.json"
    if not ranking_path.exists():
        print("  candidate_ranking.json not found, skipping enrichment")
        return candidates

    with open(ranking_path, encoding="utf-8") as f:
        ranking = json.load(f)

    ranking_map = {r["candidate_id"]: r for r in ranking}

    for c in candidates:
        r = ranking_map.get(c["design_id"])
        if r:
            c["tier"] = r.get("tier", "UNSCORED")
            c["total_score"] = r.get("total_score", 0)
            c["interface_score"] = r.get("interface_score", 0)
            c["stability_score"] = r.get("stability_score", 0)
            c["selectivity_score"] = r.get("selectivity_score", 0)
            c["developability_score"] = r.get("developability_score", 0)
            c["mechanism_score"] = r.get("mechanism_score", 0)
            c["hard_filter_pass"] = r.get("hard_filter_pass", True)
            c["red_flags"] = r.get("red_flags", [])
            rm = r.get("raw_metrics", {})
            c["seq_complexity"] = rm.get("seq_complexity", 0)
            c["charge"] = rm.get("charge", 0)
            c["hydrophobic_fraction"] = rm.get("hydrophobic_fraction", 0)
        else:
            c["tier"] = "UNSCORED"
            c["total_score"] = 0
            c["interface_score"] = 0
            c["stability_score"] = 0
            c["selectivity_score"] = 0
            c["developability_score"] = 0
            c["mechanism_score"] = 0
            c["hard_filter_pass"] = True
            c["red_flags"] = []
            c["seq_complexity"] = 0
            c["charge"] = 0
            c["hydrophobic_fraction"] = 0

    return candidates


def load_pdb_data(candidates: list[dict]) -> dict[str, str]:
    """Load PDB file contents for each candidate, keyed by design_id."""
    pdb_map = {}
    rfd_v2 = RESULTS_DIR / "rfdiffusion_v2"
    rfd_v1 = RESULTS_DIR / "rfdiffusion"

    for c in candidates:
        design = c["design_id"]
        # Try V2 directory first, then V1
        pdb_path = rfd_v2 / f"{design}.pdb"
        if not pdb_path.exists():
            pdb_path = rfd_v1 / f"{design}.pdb"
        if not pdb_path.exists():
            # Try without version prefix for V1 designs
            for p in glob.glob(str(rfd_v1 / "*.pdb")):
                if design in os.path.basename(p):
                    pdb_path = Path(p)
                    break

        if pdb_path.exists():
            with open(pdb_path, encoding="utf-8", errors="replace") as f:
                pdb_map[design] = f.read()
        else:
            print(f"  WARNING: No PDB found for {design}")

    print(f"  Loaded {len(pdb_map)} PDB structures")
    return pdb_map


def generate_html(candidates: list[dict], pdb_data: dict[str, str]) -> str:
    """Generate the full HTML dashboard."""
    candidates_json = json.dumps(candidates, indent=None)
    pdb_json = json.dumps(pdb_data, indent=None)
    hotspots_json = json.dumps(HOTSPOT_RESIDUES)

    # Read ERAP2 target sequence for the region display
    erap2_seq_path = RESULTS_DIR / "erap2_full_sequence.txt"
    erap2_seq = ""
    if erap2_seq_path.exists():
        with open(erap2_seq_path) as f:
            erap2_seq = f.read().strip()

    return HTML_TEMPLATE.replace("__CANDIDATES_JSON__", candidates_json) \
                        .replace("__PDB_DATA_JSON__", pdb_json) \
                        .replace("__HOTSPOTS_JSON__", hotspots_json) \
                        .replace("__ERAP2_SEQUENCE__", erap2_seq) \
                        .replace("__TOTAL_COUNT__", str(len(candidates)))


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Binder Explorer - Ancient Drug Discovery</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f6fa; color: #2d3436; }

.header {
    background: linear-gradient(135deg, #2d3436 0%, #636e72 100%);
    color: white; padding: 12px 24px; display: flex; align-items: center; justify-content: space-between;
}
.header h1 { font-size: 18px; font-weight: 600; }
.header .count { font-size: 14px; opacity: 0.8; }

.container { display: grid; grid-template-columns: 55fr 45fr; height: calc(100vh - 48px); gap: 0; }

/* Left Panel */
.left-panel { display: flex; flex-direction: column; overflow: hidden; border-right: 1px solid #dfe6e9; }
.toolbar {
    display: flex; align-items: center; gap: 8px; padding: 8px 12px;
    background: white; border-bottom: 1px solid #dfe6e9; flex-wrap: wrap;
}
.toolbar label { font-size: 12px; color: #636e72; font-weight: 500; }
.toolbar select {
    font-size: 12px; padding: 4px 8px; border: 1px solid #dfe6e9;
    border-radius: 4px; background: white; cursor: pointer;
}
.toolbar button {
    font-size: 12px; padding: 4px 12px; border: 1px solid #dfe6e9;
    border-radius: 4px; background: white; cursor: pointer; color: #2d3436;
}
.toolbar button:hover { background: #f0f0f0; }

.plot-container { flex: 1; min-height: 300px; }
#scatter-plot { width: 100%; height: 100%; }

.table-container {
    height: 240px; overflow-y: auto; border-top: 1px solid #dfe6e9; background: white;
}
.table-container table { width: 100%; border-collapse: collapse; font-size: 12px; }
.table-container th {
    position: sticky; top: 0; background: #f8f9fa; padding: 6px 8px;
    text-align: left; font-weight: 600; cursor: pointer; white-space: nowrap;
    border-bottom: 2px solid #dfe6e9; user-select: none;
}
.table-container th:hover { background: #e9ecef; }
.table-container th .sort-arrow { margin-left: 4px; opacity: 0.3; }
.table-container th.sorted .sort-arrow { opacity: 1; }
.table-container td { padding: 5px 8px; border-bottom: 1px solid #f1f2f6; white-space: nowrap; }
.table-container tr { cursor: pointer; }
.table-container tr:hover { background: #f8f9fa; }
.table-container tr.selected { background: #dfe6fd !important; }

.tier-badge {
    display: inline-block; padding: 1px 6px; border-radius: 3px;
    font-size: 10px; font-weight: 700; letter-spacing: 0.5px;
}
.tier-TIER_1 { background: #00b894; color: white; }
.tier-TIER_2 { background: #fdcb6e; color: #2d3436; }
.tier-TIER_3 { background: #b2bec3; color: white; }
.tier-KILL { background: #d63031; color: white; }
.tier-UNSCORED { background: #636e72; color: white; }

/* Right Panel */
.right-panel { display: flex; flex-direction: column; overflow: hidden; background: white; }

.sequence-panel {
    padding: 10px 12px; border-bottom: 1px solid #dfe6e9;
    font-family: 'Courier New', monospace; font-size: 11px; overflow-x: auto;
    max-height: 120px; overflow-y: auto; background: #fafbfc;
}
.seq-label { font-weight: 700; margin-right: 8px; display: inline-block; width: 70px; }
.seq-target { color: #d63031; }
.seq-binder { color: #00b894; }
.seq-hotspot { background: #ffeaa7; font-weight: 700; }
.seq-numbers { color: #b2bec3; font-size: 10px; }

.viewer-container { flex: 1; position: relative; min-height: 300px; }
#mol-viewer { width: 100%; height: 100%; position: relative; }
.viewer-placeholder {
    display: flex; align-items: center; justify-content: center;
    height: 100%; color: #b2bec3; font-size: 14px;
}

.info-card {
    padding: 10px 12px; border-top: 1px solid #dfe6e9;
    background: #fafbfc; font-size: 12px; min-height: 80px;
}
.info-card .info-grid { display: grid; grid-template-columns: repeat(4, 1fr); gap: 8px; }
.info-card .info-item { }
.info-card .info-label { color: #636e72; font-size: 10px; text-transform: uppercase; font-weight: 600; }
.info-card .info-value { font-size: 14px; font-weight: 600; margin-top: 2px; }
.info-card .info-flags { margin-top: 6px; color: #d63031; font-size: 11px; }
</style>
</head>
<body>

<div class="header">
    <h1>Binder Explorer &mdash; Ancient Drug Discovery</h1>
    <div class="count"><span id="shown-count">__TOTAL_COUNT__</span> of <span id="total-count">__TOTAL_COUNT__</span> results</div>
</div>

<div class="container">
    <div class="left-panel">
        <div class="toolbar">
            <label>X:</label>
            <select id="x-axis">
                <option value="iptm_erap2" selected>ipTM (ERAP2)</option>
                <option value="iptm_erap1">ipTM (ERAP1)</option>
                <option value="iptm_delta">ipTM Delta</option>
                <option value="iptm_selectivity">Selectivity Ratio</option>
                <option value="complex_plddt">Complex pLDDT</option>
                <option value="ptm">pTM</option>
                <option value="mpnn_score">MPNN Score</option>
                <option value="binder_length">Binder Length</option>
                <option value="total_score">Total Score</option>
                <option value="seq_complexity">Seq Complexity</option>
                <option value="hydrophobic_fraction">Hydrophobic Frac</option>
            </select>
            <label>Y:</label>
            <select id="y-axis">
                <option value="iptm_erap2">ipTM (ERAP2)</option>
                <option value="iptm_erap1">ipTM (ERAP1)</option>
                <option value="iptm_delta" selected>ipTM Delta</option>
                <option value="iptm_selectivity">Selectivity Ratio</option>
                <option value="complex_plddt">Complex pLDDT</option>
                <option value="ptm">pTM</option>
                <option value="mpnn_score">MPNN Score</option>
                <option value="binder_length">Binder Length</option>
                <option value="total_score">Total Score</option>
                <option value="seq_complexity">Seq Complexity</option>
                <option value="hydrophobic_fraction">Hydrophobic Frac</option>
            </select>
            <label>Color:</label>
            <select id="color-axis">
                <option value="total_score" selected>Total Score</option>
                <option value="iptm_erap2">ipTM (ERAP2)</option>
                <option value="iptm_delta">ipTM Delta</option>
                <option value="iptm_selectivity">Selectivity Ratio</option>
                <option value="complex_plddt">Complex pLDDT</option>
                <option value="mpnn_score">MPNN Score</option>
                <option value="seq_complexity">Seq Complexity</option>
            </select>
            <button onclick="downloadCSV()">Download CSV</button>
            <button onclick="resetPlot()">Reset</button>
            <button onclick="exportSVG()">Export SVG</button>
        </div>
        <div class="plot-container">
            <div id="scatter-plot"></div>
        </div>
        <div class="table-container">
            <table id="data-table">
                <thead><tr>
                    <th data-col="design_id">Design <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="version">Ver <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="family">Family <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="tier">Tier <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="iptm_erap2">ipTM(E2) <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="iptm_erap1">ipTM(E1) <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="iptm_delta">Delta <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="iptm_selectivity">Sel <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="complex_plddt">pLDDT <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="mpnn_score">MPNN <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="binder_length">Len <span class="sort-arrow">&#9650;</span></th>
                    <th data-col="total_score">Score <span class="sort-arrow">&#9650;</span></th>
                </tr></thead>
                <tbody></tbody>
            </table>
        </div>
    </div>

    <div class="right-panel">
        <div class="sequence-panel" id="sequence-panel">
            <div style="color: #b2bec3; text-align: center; padding: 20px;">
                Select a candidate to view sequence alignment
            </div>
        </div>
        <div class="viewer-container">
            <div id="mol-viewer">
                <div class="viewer-placeholder" id="viewer-placeholder">
                    Select a candidate to view 3D structure
                </div>
            </div>
        </div>
        <div class="info-card" id="info-card">
            <div style="color: #b2bec3; text-align: center;">Select a candidate to view details</div>
        </div>
    </div>
</div>

<script>
// === Data ===
const CANDIDATES = __CANDIDATES_JSON__;
const PDB_DATA = __PDB_DATA_JSON__;
const HOTSPOTS = __HOTSPOTS_JSON__;
const ERAP2_SEQ = "__ERAP2_SEQUENCE__";

let selectedIdx = null;
let sortCol = "total_score";
let sortAsc = false;
let sortedIndices = CANDIDATES.map((_, i) => i);
let viewer = null;

// === Metric display names ===
const METRIC_NAMES = {
    iptm_erap2: "ipTM (ERAP2)", iptm_erap1: "ipTM (ERAP1)",
    iptm_delta: "ipTM Delta", iptm_selectivity: "Selectivity Ratio",
    complex_plddt: "Complex pLDDT", ptm: "pTM", mpnn_score: "MPNN Score",
    binder_length: "Binder Length", total_score: "Total Score",
    seq_complexity: "Seq Complexity", hydrophobic_fraction: "Hydrophobic Frac",
};

// === Version marker shapes ===
const VERSION_SYMBOLS = { v1: "diamond", v2: "circle" };

// === Scatter Plot ===
function buildPlot() {
    const xKey = document.getElementById("x-axis").value;
    const yKey = document.getElementById("y-axis").value;
    const cKey = document.getElementById("color-axis").value;

    // Group by version for different markers
    const versions = [...new Set(CANDIDATES.map(c => c.version))];
    const traces = versions.map(ver => {
        const indices = [];
        const x = [], y = [], colors = [], texts = [], ids = [];
        CANDIDATES.forEach((c, i) => {
            if (c.version !== ver) return;
            indices.push(i);
            x.push(c[xKey] || 0);
            y.push(c[yKey] || 0);
            colors.push(c[cKey] || 0);
            ids.push(i);
            texts.push(
                `<b>${c.design_id}</b><br>` +
                `Tier: ${c.tier}<br>` +
                `ipTM(E2): ${fmt(c.iptm_erap2)}<br>` +
                `ipTM(E1): ${fmt(c.iptm_erap1)}<br>` +
                `Delta: ${fmt(c.iptm_delta)}<br>` +
                `Score: ${fmt(c.total_score)}`
            );
        });
        return {
            x, y, text: texts, customdata: ids, name: ver.toUpperCase(),
            type: "scatter", mode: "markers",
            marker: {
                size: 12, symbol: VERSION_SYMBOLS[ver] || "circle",
                color: colors,
                colorscale: [[0, "#ffeaa7"], [0.5, "#e17055"], [1, "#6c5ce7"]],
                showscale: ver === versions[0],
                colorbar: ver === versions[0] ? { title: METRIC_NAMES[cKey] || cKey, thickness: 15 } : undefined,
                line: { width: 1, color: "#2d3436" },
                opacity: 0.85,
            },
            hovertemplate: "%{text}<extra></extra>",
        };
    });

    // Highlight selected point
    if (selectedIdx !== null) {
        const c = CANDIDATES[selectedIdx];
        traces.push({
            x: [c[xKey] || 0], y: [c[yKey] || 0],
            type: "scatter", mode: "markers",
            marker: { size: 18, color: "rgba(0,0,0,0)", line: { width: 3, color: "#0984e3" } },
            showlegend: false, hoverinfo: "skip",
        });
    }

    const layout = {
        xaxis: { title: METRIC_NAMES[xKey] || xKey, gridcolor: "#f1f2f6" },
        yaxis: { title: METRIC_NAMES[yKey] || yKey, gridcolor: "#f1f2f6" },
        plot_bgcolor: "white", paper_bgcolor: "white",
        margin: { t: 10, r: 10, b: 50, l: 60 },
        legend: { x: 0.01, y: 0.99, bgcolor: "rgba(255,255,255,0.8)" },
        hovermode: "closest",
    };

    Plotly.react("scatter-plot", traces, layout, { responsive: true, displayModeBar: false });

    // Click handler
    document.getElementById("scatter-plot").removeAllListeners?.("plotly_click");
    document.getElementById("scatter-plot").on("plotly_click", function(data) {
        if (data.points.length > 0) {
            const idx = data.points[0].customdata;
            if (idx !== undefined) selectCandidate(idx);
        }
    });
}

function fmt(v) {
    if (v === null || v === undefined) return "N/A";
    return typeof v === "number" ? v.toFixed(4) : v;
}

function fmtShort(v) {
    if (v === null || v === undefined) return "N/A";
    return typeof v === "number" ? v.toFixed(3) : v;
}

// === Data Table ===
function buildTable() {
    // Sort indices
    sortedIndices = CANDIDATES.map((_, i) => i);
    sortedIndices.sort((a, b) => {
        let va = CANDIDATES[a][sortCol], vb = CANDIDATES[b][sortCol];
        if (va === null || va === undefined) va = sortAsc ? Infinity : -Infinity;
        if (vb === null || vb === undefined) vb = sortAsc ? Infinity : -Infinity;
        if (typeof va === "string") return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
        return sortAsc ? va - vb : vb - va;
    });

    const tbody = document.querySelector("#data-table tbody");
    tbody.innerHTML = "";
    sortedIndices.forEach(idx => {
        const c = CANDIDATES[idx];
        const tr = document.createElement("tr");
        tr.dataset.idx = idx;
        if (idx === selectedIdx) tr.classList.add("selected");
        tr.onclick = () => selectCandidate(idx);
        tr.innerHTML = `
            <td>${c.design_id}</td>
            <td>${c.version}</td>
            <td>${c.family}</td>
            <td><span class="tier-badge tier-${c.tier}">${c.tier}</span></td>
            <td>${fmtShort(c.iptm_erap2)}</td>
            <td>${fmtShort(c.iptm_erap1)}</td>
            <td>${fmtShort(c.iptm_delta)}</td>
            <td>${fmtShort(c.iptm_selectivity)}</td>
            <td>${fmtShort(c.complex_plddt)}</td>
            <td>${fmtShort(c.mpnn_score)}</td>
            <td>${c.binder_length}</td>
            <td>${fmtShort(c.total_score)}</td>
        `;
        tbody.appendChild(tr);
    });

    // Update sort arrows
    document.querySelectorAll("#data-table th").forEach(th => {
        th.classList.toggle("sorted", th.dataset.col === sortCol);
        const arrow = th.querySelector(".sort-arrow");
        if (th.dataset.col === sortCol) {
            arrow.innerHTML = sortAsc ? "&#9650;" : "&#9660;";
        } else {
            arrow.innerHTML = "&#9650;";
        }
    });
}

// Table header click sorting
document.querySelectorAll("#data-table th").forEach(th => {
    th.addEventListener("click", () => {
        const col = th.dataset.col;
        if (sortCol === col) {
            sortAsc = !sortAsc;
        } else {
            sortCol = col;
            sortAsc = false;
        }
        buildTable();
    });
});

// === Selection ===
function selectCandidate(idx) {
    selectedIdx = idx;
    const c = CANDIDATES[idx];

    // Update table highlight
    document.querySelectorAll("#data-table tr.selected").forEach(tr => tr.classList.remove("selected"));
    const row = document.querySelector(`#data-table tr[data-idx="${idx}"]`);
    if (row) {
        row.classList.add("selected");
        row.scrollIntoView({ block: "nearest" });
    }

    // Update plot (add selection ring)
    buildPlot();

    // Update sequence panel
    updateSequencePanel(c);

    // Update 3D viewer
    updateViewer(c);

    // Update info card
    updateInfoCard(c);
}

// === Sequence Panel ===
function updateSequencePanel(c) {
    const panel = document.getElementById("sequence-panel");
    const targetSeq = ERAP2_SEQ;
    const binderSeq = c.sequence;

    // Format target with hotspot highlighting (show region 350-500 if full seq available)
    let targetHtml = "";
    if (targetSeq.length > 0) {
        // Show residues 340-420 region (0-indexed: 339-419)
        const start = 339;
        const end = Math.min(419, targetSeq.length - 1);
        const region = targetSeq.substring(start, end + 1);
        targetHtml = '<div><span class="seq-label seq-target">Target</span><span class="seq-numbers">' + (start + 1) + ' </span>';
        for (let i = 0; i < region.length; i++) {
            const resNum = start + 1 + i;
            const isHotspot = HOTSPOTS.includes(resNum);
            targetHtml += isHotspot
                ? `<span class="seq-hotspot" title="Hotspot ${resNum}">${region[i]}</span>`
                : `<span class="seq-target">${region[i]}</span>`;
            if ((i + 1) % 10 === 0 && i < region.length - 1) targetHtml += ' ';
        }
        targetHtml += '<span class="seq-numbers"> ' + (end + 1) + '</span></div>';
    }

    // Format binder sequence
    let binderHtml = '<div style="margin-top: 4px;"><span class="seq-label seq-binder">Binder</span><span class="seq-numbers">1 </span>';
    for (let i = 0; i < binderSeq.length; i++) {
        binderHtml += `<span class="seq-binder">${binderSeq[i]}</span>`;
        if ((i + 1) % 10 === 0 && i < binderSeq.length - 1) binderHtml += ' ';
    }
    binderHtml += '<span class="seq-numbers"> ' + binderSeq.length + '</span></div>';

    // Hotspot legend
    let legend = '<div style="margin-top: 4px; font-size: 10px; color: #636e72;">' +
        'Hotspots: ' + HOTSPOTS.join(', ') +
        ' | <span class="seq-hotspot">highlighted</span> in target' +
        '</div>';

    panel.innerHTML = targetHtml + binderHtml + legend;
}

// === 3D Viewer ===
function updateViewer(c) {
    const pdb = PDB_DATA[c.design_id];
    const placeholder = document.getElementById("viewer-placeholder");

    if (!pdb) {
        if (placeholder) placeholder.style.display = "flex";
        if (viewer) { viewer.clear(); viewer.render(); }
        if (placeholder) placeholder.textContent = "No PDB structure available for " + c.design_id;
        return;
    }

    if (placeholder) placeholder.style.display = "none";

    if (!viewer) {
        const elem = document.getElementById("mol-viewer");
        viewer = $3Dmol.createViewer(elem, {
            backgroundColor: "white",
            antialias: true,
        });
    }

    viewer.clear();
    viewer.addModel(pdb, "pdb");

    // Target chain (A) in cartoon, red/maroon
    viewer.setStyle({ chain: "A" }, {
        cartoon: { color: "#b71c1c", opacity: 0.9 }
    });

    // Binder chain (B) in cartoon, green
    viewer.setStyle({ chain: "B" }, {
        cartoon: { color: "#00b894", opacity: 0.9 }
    });

    // Hotspot residues as sticks, yellow
    viewer.addStyle(
        { chain: "A", resi: HOTSPOTS },
        { stick: { color: "#f1c40f", radius: 0.15 }, cartoon: { color: "#f1c40f" } }
    );

    viewer.zoomTo();
    viewer.render();

    // Center on interface (zoom to binder + nearby target)
    try {
        viewer.zoomTo({ chain: "B" });
        viewer.zoom(0.8);
        viewer.render();
    } catch(e) {}
}

// === Info Card ===
function updateInfoCard(c) {
    const card = document.getElementById("info-card");
    const flags = (c.red_flags && c.red_flags.length > 0)
        ? `<div class="info-flags">Red flags: ${c.red_flags.join("; ")}</div>` : "";

    card.innerHTML = `
        <div style="display: flex; align-items: center; gap: 12px; margin-bottom: 6px;">
            <b>${c.design_id}</b>
            <span class="tier-badge tier-${c.tier}">${c.tier}</span>
            <span style="color: #636e72; font-size: 11px;">${c.version.toUpperCase()} | ${c.family} | ${c.binder_length} aa</span>
        </div>
        <div class="info-grid">
            <div class="info-item"><div class="info-label">ipTM (ERAP2)</div><div class="info-value">${fmtShort(c.iptm_erap2)}</div></div>
            <div class="info-item"><div class="info-label">ipTM (ERAP1)</div><div class="info-value">${fmtShort(c.iptm_erap1)}</div></div>
            <div class="info-item"><div class="info-label">Delta</div><div class="info-value" style="color: ${c.iptm_delta > 0 ? '#00b894' : '#d63031'}">${fmtShort(c.iptm_delta)}</div></div>
            <div class="info-item"><div class="info-label">Total Score</div><div class="info-value">${fmtShort(c.total_score)}</div></div>
            <div class="info-item"><div class="info-label">Interface</div><div class="info-value">${fmtShort(c.interface_score)}</div></div>
            <div class="info-item"><div class="info-label">Stability</div><div class="info-value">${fmtShort(c.stability_score)}</div></div>
            <div class="info-item"><div class="info-label">Selectivity</div><div class="info-value">${fmtShort(c.selectivity_score)}</div></div>
            <div class="info-item"><div class="info-label">Developability</div><div class="info-value">${fmtShort(c.developability_score)}</div></div>
        </div>
        ${flags}
    `;
}

// === Toolbar Actions ===
function downloadCSV() {
    const headers = ["design_id","version","family","tier","iptm_erap2","iptm_erap1",
        "iptm_delta","iptm_selectivity","complex_plddt","ptm","mpnn_score",
        "binder_length","total_score","sequence"];
    const rows = CANDIDATES.map(c => headers.map(h => c[h] ?? "").join(","));
    const csv = headers.join(",") + "\n" + rows.join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "binder_candidates.csv";
    a.click();
}

function resetPlot() {
    document.getElementById("x-axis").value = "iptm_erap2";
    document.getElementById("y-axis").value = "iptm_delta";
    document.getElementById("color-axis").value = "total_score";
    selectedIdx = null;
    buildPlot();
    buildTable();
    document.getElementById("sequence-panel").innerHTML = '<div style="color: #b2bec3; text-align: center; padding: 20px;">Select a candidate to view sequence alignment</div>';
    document.getElementById("info-card").innerHTML = '<div style="color: #b2bec3; text-align: center;">Select a candidate to view details</div>';
    if (viewer) { viewer.clear(); viewer.render(); }
    document.getElementById("viewer-placeholder").style.display = "flex";
    document.getElementById("viewer-placeholder").textContent = "Select a candidate to view 3D structure";
}

function exportSVG() {
    Plotly.downloadImage("scatter-plot", { format: "svg", filename: "binder_scatter" });
}

// Axis change listeners
["x-axis", "y-axis", "color-axis"].forEach(id => {
    document.getElementById(id).addEventListener("change", buildPlot);
});

// === Initialize ===
buildPlot();
buildTable();

// Auto-select first TIER_1 candidate
const firstTier1 = CANDIDATES.findIndex(c => c.tier === "TIER_1");
if (firstTier1 >= 0) selectCandidate(firstTier1);
</script>
</body>
</html>"""


def main():
    print("=== Binder Explorer Dashboard Generator ===")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Results dir: {RESULTS_DIR}")
    print()

    print("1. Loading validation summaries...")
    candidates = load_validation_summaries()
    print(f"   Total candidates: {len(candidates)}")
    print()

    print("2. Enriching with scoring data...")
    candidates = enrich_with_ranking(candidates)
    tiers = {}
    for c in candidates:
        tiers[c["tier"]] = tiers.get(c["tier"], 0) + 1
    print(f"   Tiers: {tiers}")
    print()

    print("3. Loading PDB structures...")
    pdb_data = load_pdb_data(candidates)
    print()

    print("4. Generating HTML...")
    html = generate_html(candidates, pdb_data)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    size_mb = os.path.getsize(OUTPUT_HTML) / 1024 / 1024
    print(f"   Written: {OUTPUT_HTML} ({size_mb:.1f} MB)")
    print()
    print("Done! Open dashboard/index.html in your browser.")


if __name__ == "__main__":
    main()
