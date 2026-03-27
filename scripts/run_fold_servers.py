#!/usr/bin/env python3
"""Submit v2 candidates to free fold prediction servers and collect scores.

Run on any machine with internet access:
    pip install requests pyyaml biopython
    python scripts/run_fold_servers.py

Hits ESMFold API (free, no account needed) for all 4 candidates, extracts
pLDDT scores, and updates v2/fold_scores.yaml automatically.

For Chai-1 and Boltz-2, opens browser tabs to the web servers since they
don't have free public APIs — paste the sequences manually there.
"""

import json
import os
import sys
import time
import webbrowser
from pathlib import Path

import requests

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

CANDIDATES = {
    "n248_trim_c5_Y87A_Y89A": "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN",
    "n248_trim_c5":           "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKN",
    "n248_wt":                "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIYNYLKNYFFEK",
    "n248_ko_all_aromatics":  "DIRHAAKSLEEALKNLPKVVDMLVDLASKGIAHLDNTNILVKDDKAAAIDAGSAAINEKKSTDATLKIKNDQISSEEAVKSVSEKIANALKNAAAEK",
}

ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
OUTPUT_DIR = Path("data/results/esmfold_v2")
SCORES_FILE = Path("v2/fold_scores.yaml")


def extract_plddt_from_pdb(pdb_text: str) -> float:
    """Extract mean pLDDT from B-factor column of PDB ATOM records."""
    bfactors = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                bfactors.append(bfactor)
            except ValueError:
                continue
    if not bfactors:
        return 0.0
    return sum(bfactors) / len(bfactors)


def extract_ca_plddt(pdb_text: str) -> list[float]:
    """Extract per-residue pLDDT from CA atoms only."""
    ca_bfactors = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                ca_bfactors.append(bfactor)
            except ValueError:
                continue
    return ca_bfactors


def run_esmfold(name: str, sequence: str) -> dict:
    """Submit sequence to ESMFold API, return scores."""
    print(f"  Submitting {name} ({len(sequence)}aa) to ESMFold API...")

    try:
        resp = requests.post(
            ESMFOLD_API,
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=120,
        )
        resp.raise_for_status()
    except requests.exceptions.ConnectionError:
        print(f"    ERROR: Cannot reach ESMFold API. Server may be down.")
        print(f"    Try: https://esmatlas.com/resources?action=fold")
        return {"plddt": None, "error": "connection_failed"}
    except requests.exceptions.Timeout:
        print(f"    ERROR: ESMFold API timed out (>120s)")
        return {"plddt": None, "error": "timeout"}
    except requests.exceptions.HTTPError as e:
        print(f"    ERROR: ESMFold API returned {e.response.status_code}")
        return {"plddt": None, "error": str(e)}

    pdb_text = resp.text
    mean_plddt = extract_plddt_from_pdb(pdb_text)
    ca_plddt = extract_ca_plddt(pdb_text)

    # Save PDB
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdb_path = OUTPUT_DIR / f"{name}.pdb"
    pdb_path.write_text(pdb_text)

    # Per-residue stats
    low_conf = sum(1 for v in ca_plddt if v < 50)
    high_conf = sum(1 for v in ca_plddt if v > 70)
    n_res = len(ca_plddt)

    result = {
        "plddt": round(mean_plddt, 1),
        "residues": n_res,
        "plddt_above_70": high_conf,
        "plddt_below_50": low_conf,
        "pdb_file": str(pdb_path),
    }

    print(f"    pLDDT: {mean_plddt:.1f}  ({high_conf}/{n_res} residues >70, {low_conf}/{n_res} <50)")
    print(f"    Saved: {pdb_path}")
    return result


def update_scores_yaml(esmfold_results: dict):
    """Update fold_scores.yaml with ESMFold results."""
    if not SCORES_FILE.exists():
        print(f"\n  {SCORES_FILE} not found — run: python scripts/test_v2_candidate_fold.py --init")
        return

    if HAS_YAML:
        with open(SCORES_FILE) as f:
            scores = yaml.safe_load(f)

        for name, result in esmfold_results.items():
            if name in scores and result.get("plddt") is not None:
                scores[name]["esmfold"]["plddt"] = result["plddt"]

        with open(SCORES_FILE, "w") as f:
            yaml.dump(scores, f, default_flow_style=False, sort_keys=False, width=120)

        print(f"\n  Updated {SCORES_FILE} with ESMFold pLDDT scores")
    else:
        print(f"\n  Install pyyaml to auto-update: pip install pyyaml")


def open_web_servers():
    """Open browser tabs for servers that need manual submission."""
    print("\n" + "=" * 70)
    print("MANUAL SUBMISSION — paste sequences into these web servers:")
    print("=" * 70)

    servers = [
        ("Chai Discovery", "https://lab.chaidiscovery.com"),
        ("Boltz-2 (Tamarind)", "https://tamarind.bio/tools/boltz"),
        ("AlphaFold Server", "https://alphafoldserver.com"),
    ]

    for server_name, url in servers:
        print(f"\n  {server_name}: {url}")

    print("\nSequences to paste:")
    for name, seq in CANDIDATES.items():
        print(f"\n  {name} ({len(seq)}aa):")
        print(f"  {seq}")

    answer = input("\nOpen browser tabs for these servers? [y/N] ")
    if answer.lower() in ("y", "yes"):
        for _, url in servers:
            webbrowser.open(url)
            time.sleep(0.5)
        print("  Opened browser tabs. Paste sequences and record scores in v2/fold_scores.yaml")


def main():
    print("=" * 70)
    print("  V2 CANDIDATE FOLD PREDICTION")
    print("=" * 70)
    print()

    # Phase 1: ESMFold API (automatic)
    print("Phase 1: ESMFold API (automatic, free, no account)")
    print("-" * 50)

    esmfold_results = {}
    for name, seq in CANDIDATES.items():
        result = run_esmfold(name, seq)
        esmfold_results[name] = result
        time.sleep(1)  # Be polite to the API

    # Summary
    print("\n" + "-" * 50)
    print("ESMFold Summary:")
    print(f"  {'Candidate':<30} {'pLDDT':>8} {'High(>70)':>10} {'Low(<50)':>10}")
    print("-" * 60)
    for name in CANDIDATES:
        r = esmfold_results[name]
        plddt = r.get("plddt")
        if plddt is not None:
            high = r.get("plddt_above_70", "?")
            low = r.get("plddt_below_50", "?")
            nres = r.get("residues", "?")
            print(f"  {name:<30} {plddt:>8.1f} {high:>5}/{nres:<4} {low:>5}/{nres:<4}")
        else:
            print(f"  {name:<30} {'FAILED':>8}")

    # Save raw results
    results_json = OUTPUT_DIR / "esmfold_results.json"
    with open(results_json, "w") as f:
        json.dump(esmfold_results, f, indent=2)
    print(f"\n  Raw results: {results_json}")

    # Update YAML
    update_scores_yaml(esmfold_results)

    # Phase 2: Manual web servers
    open_web_servers()

    print("\n" + "=" * 70)
    print("  After filling in all scores, run:")
    print("  python scripts/test_v2_candidate_fold.py")
    print("=" * 70)


if __name__ == "__main__":
    main()
