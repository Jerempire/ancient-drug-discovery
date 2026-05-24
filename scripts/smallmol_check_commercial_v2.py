"""Follow-up: for each candidate, find the closest commercially-available analog
and compute the actual Tanimoto similarity (Morgan fingerprints, radius 2)."""
from __future__ import annotations

try:
    import truststore; truststore.inject_into_ssl()
except ImportError:
    pass

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/results/smallmol/commercial_availability_v2.json"

CANDIDATES = [
    ("pyridine-azabicyclic (champion)",  "CC12Cc3ccncc3C(C1)C(=O)C1CCCC=C12"),
    ("polycyclic acetamide (orig lead)", "CC(=O)N1CC2CC(=O)C3C2(C)CCC31C(=O)O"),
    ("salicylate diacid-amide",          "CCC(C(=O)NC(C)=O)c1cccc(O)c1C(=O)O"),
]

VENDOR_KW = ["enamine", "mcule", "molport", "sigma", "aldrich", "ambeed",
             "chembridge", "asinex", "vitas", "specs", "chemspace", "lifechem",
             "bld pharm", "tcil", " tci", "alfa aesar", "fluorochem",
             "matrix scientific", "key organics", "otava", "uorsy", "indofine",
             "ark pharm", "combi-blocks", "combi blocks", "carbosynth",
             "biosynth", "akos", "amerigen", "selleck"]

UA = {"User-Agent": "smallmol-check-v2/1.0"}


def fetch_json(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def similarity_cids(smi, threshold, max_records=25):
    esc = urllib.parse.quote(smi, safe="")
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/"
           f"smiles/{esc}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}")
    try:
        d = fetch_json(url, timeout=120)
        return [c for c in d.get("IdentifierList", {}).get("CID", []) if c != 0]
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise


def get_props(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight/JSON"
    try:
        d = fetch_json(url)
        return d.get("PropertyTable", {}).get("Properties", [{}])[0]
    except Exception:
        return {}


def get_vendors(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/SourceName/JSON"
    try:
        d = fetch_json(url)
        return d.get("InformationList", {}).get("Information", [{}])[0].get("SourceName", [])
    except Exception:
        return []


def commercial(srcs):
    out = []
    for s in srcs:
        sl = s.lower()
        for kw in VENDOR_KW:
            if kw in sl:
                out.append(s); break
    return out


def tanimoto(smi_a, smi_b):
    a = Chem.MolFromSmiles(smi_a); b = Chem.MolFromSmiles(smi_b)
    if a is None or b is None: return None
    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, 2048)
    fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, 2048)
    return round(DataStructs.TanimotoSimilarity(fa, fb), 3)


def main():
    payload = []
    for label, smi in CANDIDATES:
        print(f"\n[{label}]")
        print(f"  Query SMILES: {smi}")
        target_mol = Chem.MolFromSmiles(smi)
        target_mw = round(Descriptors.MolWt(target_mol), 1)
        print(f"  Target MW: {target_mw}")

        row = {"label": label, "query_smiles": smi, "target_mw": target_mw,
              "thresholds": {}}

        for threshold in [90, 85, 75, 65]:
            print(f"\n  -- {threshold}% similarity --")
            cids = similarity_cids(smi, threshold, max_records=25)
            print(f"    {len(cids)} CIDs returned")
            if not cids:
                row["thresholds"][threshold] = {"n_cids": 0, "commercial_hits": []}
                time.sleep(0.5)
                continue

            commercial_hits = []
            for cid in cids[:15]:  # check top 15
                props = get_props(cid)
                vendors = get_vendors(cid)
                comm = commercial(vendors)
                if not comm:
                    continue
                hit_smi = props.get("IsomericSMILES") or props.get("CanonicalSMILES")
                if not hit_smi:
                    continue
                tan = tanimoto(smi, hit_smi)
                commercial_hits.append({
                    "cid": cid,
                    "smiles": hit_smi,
                    "formula": props.get("MolecularFormula"),
                    "mw": props.get("MolecularWeight"),
                    "tanimoto_morgan_r2": tan,
                    "mw_delta": round(float(props.get("MolecularWeight", 0)) - target_mw, 1) if props.get("MolecularWeight") else None,
                    "commercial_vendors": comm,
                    "n_total_sources": len(vendors),
                })
                time.sleep(0.3)

            commercial_hits.sort(key=lambda h: (-(h["tanimoto_morgan_r2"] or 0), abs(h["mw_delta"] or 0)))
            row["thresholds"][threshold] = {
                "n_cids": len(cids), "n_checked": min(15, len(cids)),
                "commercial_hits": commercial_hits[:5]
            }
            print(f"    {len(commercial_hits)} commercial-available analogs found")
            for h in commercial_hits[:5]:
                print(f"      Tanimoto={h['tanimoto_morgan_r2']}  MW={h['mw']} (Δ{h['mw_delta']:+})")
                print(f"        SMILES: {h['smiles']}")
                print(f"        Vendors: {', '.join(set(h['commercial_vendors']))}")
            if commercial_hits:
                # Stop at first threshold that gave us decent hits
                if any(h["tanimoto_morgan_r2"] and h["tanimoto_morgan_r2"] >= 0.6 for h in commercial_hits):
                    break
            time.sleep(0.5)

        payload.append(row)
        time.sleep(1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved: {OUT}")

    # Final summary
    print("\n" + "=" * 80)
    print("SUMMARY - closest commercially-available analog per candidate")
    print("=" * 80)
    for r in payload:
        best = None
        for t in [90, 85, 75, 65]:
            hits = r["thresholds"].get(t, {}).get("commercial_hits", [])
            if hits:
                cand = max(hits, key=lambda h: h["tanimoto_morgan_r2"] or 0)
                if best is None or (cand["tanimoto_morgan_r2"] or 0) > (best["tanimoto_morgan_r2"] or 0):
                    best = cand
                if best and best["tanimoto_morgan_r2"] and best["tanimoto_morgan_r2"] >= 0.6:
                    break
        print(f"\n{r['label']}:")
        if not best:
            print("  NO commercial analog at ≥65% similarity - CUSTOM SYNTHESIS REQUIRED")
        else:
            print(f"  Closest: Tanimoto {best['tanimoto_morgan_r2']}, MW Δ{best['mw_delta']:+}")
            print(f"  SMILES:  {best['smiles']}")
            print(f"  Vendors: {', '.join(set(best['commercial_vendors']))}")
            print(f"  PubChem CID: {best['cid']} (https://pubchem.ncbi.nlm.nih.gov/compound/{best['cid']})")


if __name__ == "__main__":
    main()
