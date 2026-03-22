"""
v311_erap1_cleanup.py — Polish Y87S: reduce ERAP1 without breaking K392/delta/ANPEP.

9 variants (6 single + 3 double) targeting peripheral aromatics and hydrophobics
that likely contact ERAP1 via conserved aminopeptidase surface.

Usage: python3 /workspace/scripts/v311_erap1_cleanup.py
"""
import glob, json, os, subprocess, sys, time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v311_erap1_cleanup")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")
TARGET_REGION = (350, 500)

AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

# Y87S parent (n243_id_12 + Y87S)
Y87S = "SKALKEFLSNLNKAEDYKNKGNLAFNNGNYSDAISFYKKSLSELNKAKTIINNDKNLKKMLDNKTYLGKIYQNLEKLVTNNLKAAQSYKNNP"

def get_seq(pdb, s, e):
    from Bio.PDB import PDBParser
    p = PDBParser(QUIET=True); st = p.get_structure("s", pdb)
    return "".join(AA3TO1[r.get_resname()] for c in st.get_chains() if c.id=="A" for r in c.get_residues() if r.get_resname() in AA3TO1 and s <= r.get_id()[1] <= e)

def mut(seq, muts):
    s = list(seq)
    for pos, aa in muts:
        s[pos] = aa
    return "".join(s)

def write_yaml(t, b, path):
    with open(path, "w") as f:
        f.write(f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {t}\n      msa: empty\n  - protein:\n      id: B\n      sequence: {b}\n      msa: empty\n")

def run_boltz(yaml_path, out_dir):
    for cmd in ["/opt/conda/bin/boltz", "boltz"]:
        if os.path.exists(cmd) or cmd == "boltz":
            return subprocess.run([cmd, "predict", yaml_path, "--out_dir", out_dir, "--recycling_steps", "3", "--diffusion_samples", "1", "--accelerator", "gpu", "--devices", "1"], capture_output=True, text=True, timeout=1800, env=os.environ.copy())

def parse_scores(d):
    for jf in glob.glob(os.path.join(d, "**", "*.json"), recursive=True):
        if "manifest" in jf: continue
        try:
            data = json.load(open(jf))
            if "iptm" in data: return data.get("iptm", 0)
        except: pass
    return 0

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t0 = time.time()

    e2k = get_seq(ERAP2_PDB, *TARGET_REGION)
    e2n = list(e2k); e2n[392-350] = "N"; e2n = "".join(e2n)
    e1 = get_seq(ERAP1_PDB, *TARGET_REGION)
    targets = {"k392": e2k, "n392": e2n, "erap1": e1}

    # Fetch IRAP/ANPEP
    try:
        import urllib.request
        for uid, nm in [("Q9UIQ6","irap"),("P15144","anpep")]:
            r = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            s = "".join(l for l in r.read().decode().strip().split("\n") if not l.startswith(">"))
            targets[nm] = s[349:550]
    except: pass

    for n, s in targets.items():
        print(f"  {n}: {len(s)} aa")

    # 9 variants + parent control
    variants = [
        ("y87s_parent", Y87S, "Y87S parent baseline"),
        ("y87s_Y17S", mut(Y87S, [(16,"S")]), "Y17S — N-term aromatic"),
        ("y87s_Y37S", mut(Y87S, [(36,"S")]), "Y37S — patch 33-38 aromatic"),
        ("y87s_L8T", mut(Y87S, [(7,"T")]), "L8T — N-term hydrophobic"),
        ("y87s_L4T", mut(Y87S, [(3,"T")]), "L4T — N-terminal"),
        ("y87s_Y66S", mut(Y87S, [(65,"S")]), "Y66S — C-term near Y87"),
        ("y87s_F36Q", mut(Y87S, [(35,"Q")]), "F36Q — patch 33-38"),
        ("y87s_Y17S_Y37S", mut(Y87S, [(16,"S"),(36,"S")]), "Y17S+Y37S — two non-core aromatics"),
        ("y87s_Y17S_L8T", mut(Y87S, [(16,"S"),(7,"T")]), "Y17S+L8T — clean N-term"),
        ("y87s_F36Q_Y66S", mut(Y87S, [(35,"Q"),(65,"S")]), "F36Q+Y66S — patch + C-term"),
    ]

    print(f"\n{'='*110}")
    print(f"V3.1.1 ERAP1 CLEANUP — {len(variants)} variants x {len(targets)} targets")
    print(f"{'='*110}")

    results = []
    for i, (name, seq, desc) in enumerate(variants):
        print(f"\n[{i+1}/{len(variants)}] {name} — {desc}")
        r = {"name": name, "desc": desc, "length": len(seq), "sequence": seq, "screens": {}}
        for tn, ts in targets.items():
            yp = os.path.join(RESULTS_DIR, f"{name}_{tn}.yaml")
            od = os.path.join(RESULTS_DIR, f"{name}_{tn}")
            write_yaml(ts, seq, yp)
            print(f"  vs {tn}...", end=" ", flush=True)
            t = time.time()
            br = run_boltz(yp, od)
            el = time.time() - t
            if br and br.returncode == 0:
                ip = parse_scores(od)
                print(f"OK ({el:.0f}s) ipTM={ip:.4f}")
                r["screens"][tn] = ip
            else:
                print(f"FAILED ({el:.0f}s)")
                r["screens"][tn] = 0
        r["delta"] = round(r["screens"].get("k392",0) - r["screens"].get("n392",0), 4)
        results.append(r)

    # Summary
    parent = results[0]
    pe1 = parent["screens"].get("erap1", 0)
    print(f"\n{'='*120}")
    print(f"V3.1.1 RESULTS ({(time.time()-t0)/60:.1f} min)")
    print(f"{'='*120}")
    print(f"{'Name':<25} {'K392':>7} {'N392':>7} {'Delta':>7} {'E1':>7} {'E1 chg':>7} {'IRAP':>7} {'ANPEP':>7} {'Verdict'}")
    print("-"*120)
    for r in sorted(results, key=lambda x: x["screens"].get("erap1", 1)):
        s = r["screens"]
        k = s.get("k392",0); e1 = s.get("erap1",0)
        e1c = e1 - pe1; d = r["delta"]
        v = "IMPROVED" if k>=0.65 and d>=0.15 and e1<pe1-0.03 else "OK" if k>=0.65 and d>=0.15 else "K392 lost" if k<0.55 else "delta lost" if d<0.10 else "no gain"
        print(f"{r['name']:<25} {k:>7.4f} {s.get('n392',0):>7.4f} {d:>+7.4f} {e1:>7.4f} {e1c:>+7.4f} {s.get('irap',0):>7.4f} {s.get('anpep',0):>7.4f} {v}")

    with open(os.path.join(RESULTS_DIR, "v311_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {os.path.join(RESULTS_DIR, 'v311_results.json')}")

if __name__ == "__main__":
    main()
