"""
Check commercial availability of the top 3 candidates via PubChem REST API.

For each SMILES:
  1. Try exact-match lookup (canonical SMILES against PubChem)
  2. If no exact match, do 2D similarity search at 85% Tanimoto threshold
  3. For each hit CID, fetch vendor list (PubChem aggregates ~200 vendors)
  4. Flag commercial vendors (Enamine, Mcule, Sigma, Aldrich, MolPort, ChemBridge, etc.)
"""
from __future__ import annotations

# Windows SSL fix
try:
    import truststore
    truststore.inject_into_ssl()
except ImportError:
    pass

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/results/smallmol/commercial_availability.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

CANDIDATES = [
    ("pyridine-azabicyclic (champion)",      "CC12Cc3ccncc3C(C1)C(=O)C1CCCC=C12"),
    ("polycyclic acetamide (orig lead)",     "CC(=O)N1CC2CC(=O)C3C2(C)CCC31C(=O)O"),
    ("salicylate diacid-amide",              "CCC(C(=O)NC(C)=O)c1cccc(O)c1C(=O)O"),
]

VENDOR_KEYWORDS = [
    "enamine", "mcule", "molport", "sigma-aldrich", "sigma aldrich",
    "aldrich", "ambeed", "chembridge", "asinex", "vitas", "specs",
    "chemspace", "lifechem", "bld pharm", "tcil", "tci", "alfa aesar",
    "fluorochem", "matrix scientific", "key organics", "otava",
    "uorsy", "indofine", "ark pharm", "combi-blocks", "combi blocks",
    "carbosynth", "biosynth", "akos", "amerigen", "selleck"
]

UA = {"User-Agent": "smallmol-availability-check/1.0"}


def fetch_json(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def fetch_text(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="ignore")


def exact_match(smi: str):
    """PubChem exact-match by canonical SMILES."""
    esc = urllib.parse.quote(smi, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{esc}/cids/JSON"
    try:
        d = fetch_json(url)
        return d.get("IdentifierList", {}).get("CID", [])
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise


def similarity_search(smi: str, threshold=85, max_records=20):
    """PubChem 2D similarity search (async, polls until done)."""
    esc = urllib.parse.quote(smi, safe="")
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/"
           f"smiles/{esc}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}")
    try:
        d = fetch_json(url, timeout=120)
        return d.get("IdentifierList", {}).get("CID", [])
    except urllib.error.HTTPError as e:
        print(f"    HTTP {e.code}: {e.reason}")
        return []


def get_vendors(cid: int):
    """Fetch vendor/source names for a CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/SourceName/JSON"
    try:
        d = fetch_json(url)
        info = d.get("InformationList", {}).get("Information", [{}])[0]
        return info.get("SourceName", [])
    except Exception:
        return []


def get_iupac(cid: int):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"
    try:
        d = fetch_json(url)
        return d.get("PropertyTable", {}).get("Properties", [{}])[0]
    except Exception:
        return {}


def commercial_vendors(sources):
    matches = []
    for s in sources:
        sl = s.lower()
        for kw in VENDOR_KEYWORDS:
            if kw in sl:
                matches.append(s)
                break
    return matches


def main():
    print("Commercial availability check via PubChem")
    print("=" * 80)
    payload = []
    for label, smi in CANDIDATES:
        print(f"\n[{label}]")
        print(f"  SMILES: {smi}")
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
        if canon != smi:
            print(f"  Canonical: {canon}")

        row = {"label": label, "smiles_input": smi, "smiles_canonical": canon}

        # 1. Exact match
        print("  -> exact-match search ...")
        cids = [c for c in exact_match(canon) if c != 0]
        if cids:
            print(f"    EXACT MATCH! PubChem CID(s): {cids[:5]}")
            row["exact_match_cids"] = cids
            for cid in cids[:3]:
                src = get_vendors(cid)
                comm = commercial_vendors(src)
                prop = get_iupac(cid)
                print(f"    CID {cid}: {prop.get('MolecularFormula', '?')}  MW={prop.get('MolecularWeight', '?')}")
                print(f"      IUPAC: {prop.get('IUPACName', '(none)')}")
                print(f"      {len(src)} sources; {len(comm)} commercial vendors:")
                for v in comm[:10]:
                    print(f"        - {v}")
                row.setdefault("exact_hits", []).append({
                    "cid": cid, "properties": prop, "all_sources": src,
                    "commercial_vendors": comm
                })
            time.sleep(1)
            continue
        else:
            print("    no exact match")

        # 2. Similarity search at 85%
        print("  -> 85% similarity search ...")
        sim = similarity_search(canon, threshold=85, max_records=10)
        if not sim:
            print("    no 85% similar compounds")
            row["similar_85"] = []
            payload.append(row)
            continue
        print(f"    {len(sim)} similar CIDs: {sim[:5]}...")
        row["similar_85_cids"] = sim
        # Check first 5 for vendors
        for cid in sim[:5]:
            src = get_vendors(cid)
            comm = commercial_vendors(src)
            prop = get_iupac(cid)
            print(f"    CID {cid}: {prop.get('MolecularFormula','?')}  MW={prop.get('MolecularWeight','?')}")
            print(f"      SMILES: {prop.get('CanonicalSMILES', '?')}")
            print(f"      {len(src)} sources; {len(comm)} commercial vendors:")
            for v in comm[:8]:
                print(f"        - {v}")
            row.setdefault("similar_hits", []).append({
                "cid": cid, "properties": prop, "all_sources": src,
                "commercial_vendors": comm, "tanimoto_min": 0.85
            })
            time.sleep(0.5)

        payload.append(row)
        time.sleep(1)

    with open(OUT, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\n\nSaved: {OUT}")

    # Summary
    print("\n=== SUMMARY ===")
    for r in payload:
        if r.get("exact_match_cids"):
            best_comm = sum(len(h["commercial_vendors"]) for h in r.get("exact_hits", []))
            print(f"  {r['label']}: EXACT MATCH, {best_comm} commercial vendors aggregate")
        elif r.get("similar_85_cids"):
            n_with_vendors = sum(1 for h in r.get("similar_hits", []) if h["commercial_vendors"])
            print(f"  {r['label']}: {len(r['similar_85_cids'])} similar, {n_with_vendors} have commercial vendors")
        else:
            print(f"  {r['label']}: NO HITS — needs custom synthesis")


if __name__ == "__main__":
    main()
