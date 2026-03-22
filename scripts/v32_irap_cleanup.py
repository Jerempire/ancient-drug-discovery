"""
v32_irap_cleanup.py — IRAP-targeted cleanup on Y87S+Y17S+L8T.

Tri-paralog divergent surface analysis (archaeogenetics/erap2_tri_divergent.py)
found that 10 of 30 V3 binder contact residues sit on ERAP2 positions that are
ERAP1-divergent but IRAP-conserved. That's why fixing ERAP1 caused IRAP to spike.

Strategy: mutate binder residues in the C-terminal region (where IRAP contacts
cluster) to break hydrophobic complementarity with IRAP-conserved ERAP2 patches
(297, 305, 332, 360, 401, 402, 822, 864, 876, 892), while preserving contacts
to tri-divergent positions (389, 392, 403, 405).

Parent: Y87S+Y17S+L8T (K392=0.779, delta=+0.282, ERAP1=0.268, IRAP=0.555, ANPEP=0.189)

Usage: python3 /workspace/scripts/v32_irap_cleanup.py
"""
import glob, json, os, subprocess, sys, time

WORKSPACE = "/workspace"
RESULTS_DIR = os.path.join(WORKSPACE, "results", "v32_irap_cleanup")
STRUCTURES = os.path.join(WORKSPACE, "data", "structures")
ERAP2_PDB = os.path.join(STRUCTURES, "erap2_wt_alphafold.pdb")
ERAP1_PDB = os.path.join(STRUCTURES, "erap1_wt_alphafold.pdb")
TARGET_REGION = (350, 500)

AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

# Y87S+Y17S+L8T parent (n243_id_12 with 3 mutations)
# Original: SKALKEFLSNLNKAEDYKNKGNLAFNNGNYSDAISFYKKSLSELNKAKTIINNDKNLKKMLDNKTYLGKIYQNLEKLVTNNLKAAQSYKNNP
# Y87S: pos 87 (idx 86) Y->S
# Y17S: pos 17 (idx 16) Y->S
# L8T:  pos 8  (idx 7)  L->T
PARENT = "SKALKEFTSNLNKAEDSKNKGNLAFNNGNYSDAISFYKKSLSELNKAKTIINNDKNLKKMLDNKTYLGKIYQNLEKLVTNNLKAAQSYKNNP"
# Verify: pos 8=T, pos 17=S, pos 87=S (check below)

def get_seq(pdb, s, e):
    from Bio.PDB import PDBParser
    p = PDBParser(QUIET=True); st = p.get_structure("s", pdb)
    return "".join(AA3TO1[r.get_resname()] for c in st.get_chains() if c.id=="A"
                   for r in c.get_residues() if r.get_resname() in AA3TO1 and s <= r.get_id()[1] <= e)

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
            return subprocess.run([cmd, "predict", yaml_path, "--out_dir", out_dir,
                "--recycling_steps", "3", "--diffusion_samples", "3", "--seed", "42",
                "--accelerator", "gpu", "--devices", "1"],
                capture_output=True, text=True, timeout=1800, env=os.environ.copy())

def parse_scores(d):
    scores = []
    for jf in glob.glob(os.path.join(d, "**", "*.json"), recursive=True):
        if "manifest" in jf: continue
        try:
            data = json.load(open(jf))
            if "iptm" in data: scores.append(data["iptm"])
        except: pass
    return sum(scores) / len(scores) if scores else 0

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t0 = time.time()

    # Verify parent sequence
    assert PARENT[7] == "T", f"L8T not applied: pos 8 is {PARENT[7]}"
    assert PARENT[16] == "S", f"Y17S not applied: pos 17 is {PARENT[16]}"
    assert PARENT[86] == "S", f"Y87S not applied: pos 87 is {PARENT[86]}"
    assert len(PARENT) == 92, f"Length {len(PARENT)} != 92"
    print(f"Parent Y87S+Y17S+L8T verified: {len(PARENT)} aa")
    print(f"  Mutations: L8T={PARENT[7]}, Y17S={PARENT[16]}, Y87S={PARENT[86]}")

    e2k = get_seq(ERAP2_PDB, *TARGET_REGION)
    e2n = list(e2k); e2n[392-350] = "N"; e2n = "".join(e2n)
    e1 = get_seq(ERAP1_PDB, *TARGET_REGION)
    targets = {"k392": e2k, "n392": e2n, "erap1": e1}

    try:
        import urllib.request
        for uid, nm in [("Q9UIQ6","irap"),("P15144","anpep")]:
            r = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            s = "".join(l for l in r.read().decode().strip().split("\n") if not l.startswith(">"))
            targets[nm] = s[349:550]
    except: pass

    for n, s in targets.items():
        print(f"  {n}: {len(s)} aa")

    # V3.2 IRAP cleanup variants
    #
    # The IRAP spike (0.555 -> 0.789 with Y17S alone) came from removing Y17
    # which was a contact suppressing IRAP. On the parent Y87S+Y17S+L8T, IRAP=0.555.
    #
    # IRAP-conserved ERAP2 positions the binder likely contacts: 401, 402
    # (right next to K392), plus 297, 305, 332, 360 (further out).
    #
    # The binder's C-terminal region (residues ~60-92) docks near the K392 channel.
    # Residues in the hydrophobic cluster 66-78 (TYLGKIYQNLEKLV) likely make
    # contacts with IRAP-conserved 401/402. The N-terminal region contacts
    # positions further from K392 including some IRAP-conserved patches.
    #
    # Targets for IRAP reduction:
    # - I69 (C-term hydrophobic, likely contacts 401/402 area)
    # - K70 (positive charge near 401/402, IRAP has similar pocket)
    # - Q72 (mid C-term, may bridge to IRAP-conserved patch)
    # - L75 (in LEKLV patch, hydrophobic contact to 401/402 region)
    # - V78 (terminal of LEKLV, anchors into channel)
    # - K47 (mid-region positive, may contact 360/332 area)
    #
    # Strategy: disrupt hydrophobic/charge complementarity at these positions
    # while leaving the tri-divergent contacts (389, 392, 403, 405) untouched.
    # The K392 discrimination depends on contacts to 392 itself — those are safe.

    variants = [
        ("v32_parent", PARENT, "Y87S+Y17S+L8T parent baseline"),
        # Single mutants — target C-term IRAP contacts
        ("v32_I69T",  mut(PARENT, [(68,"T")]),  "I69T — C-term hydrophobic near 401/402"),
        ("v32_K70E",  mut(PARENT, [(69,"E")]),  "K70E — flip charge at IRAP pocket"),
        ("v32_Q72N",  mut(PARENT, [(71,"N")]),  "Q72N — shorten sidechain, reduce IRAP fit"),
        ("v32_L75T",  mut(PARENT, [(74,"T")]),  "L75T — LEKLV patch hydrophobic"),
        ("v32_V78T",  mut(PARENT, [(77,"T")]),  "V78T — LEKLV anchor into channel"),
        ("v32_K47E",  mut(PARENT, [(46,"E")]),  "K47E — mid-region, may contact 360/332"),
        # Double mutants — combine best single targets
        ("v32_I69T_L75T",  mut(PARENT, [(68,"T"),(74,"T")]),  "I69T+L75T — clean C-term hydrophobics"),
        ("v32_K70E_V78T",  mut(PARENT, [(69,"E"),(77,"T")]),  "K70E+V78T — charge flip + anchor"),
        ("v32_I69T_K47E",  mut(PARENT, [(68,"T"),(46,"E")]),  "I69T+K47E — C-term + mid-region"),
    ]

    print(f"\n{'='*120}")
    print(f"V3.2 IRAP CLEANUP — {len(variants)} variants x {len(targets)} targets")
    print(f"{'='*120}")

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

    # Summary — sort by IRAP (ascending)
    parent = results[0]
    p_irap = parent["screens"].get("irap", 0)
    print(f"\n{'='*130}")
    print(f"V3.2 RESULTS ({(time.time()-t0)/60:.1f} min)")
    print(f"{'='*130}")
    print(f"{'Name':<25} {'K392':>7} {'N392':>7} {'Delta':>7} {'E1':>7} {'IRAP':>7} {'IRAP chg':>9} {'ANPEP':>7} {'Verdict'}")
    print("-"*130)

    for r in sorted(results, key=lambda x: x["screens"].get("irap", 1)):
        s = r["screens"]
        k = s.get("k392",0); irap = s.get("irap",0)
        ic = irap - p_irap; d = r["delta"]
        # Verdict: need K392>=0.65, delta>=0.15, IRAP lower than parent
        if k >= 0.65 and d >= 0.15 and irap < p_irap - 0.05:
            v = "IMPROVED"
        elif k >= 0.65 and d >= 0.15:
            v = "OK"
        elif k < 0.55:
            v = "K392 lost"
        elif d < 0.10:
            v = "delta lost"
        else:
            v = "no gain"
        print(f"{r['name']:<25} {k:>7.4f} {s.get('n392',0):>7.4f} {d:>+7.4f} {s.get('erap1',0):>7.4f} {irap:>7.4f} {ic:>+9.4f} {s.get('anpep',0):>7.4f} {v}")

    # Check if any variant clears ALL gates
    print(f"\n--- ALL-GATE CHECK (K392>=0.65, delta>=0.15, E1<0.35, IRAP<0.35, ANPEP<0.30) ---")
    for r in results:
        s = r["screens"]
        k=s.get("k392",0); d=r["delta"]; e1=s.get("erap1",0); ir=s.get("irap",0); an=s.get("anpep",0)
        gates = k>=0.65 and d>=0.15 and e1<0.35 and ir<0.35 and an<0.30
        if gates:
            print(f"  PASS: {r['name']} — K392={k:.3f} delta={d:+.3f} E1={e1:.3f} IRAP={ir:.3f} ANPEP={an:.3f}")
    else:
        # Check if any pass was printed
        any_pass = any(
            s.get("k392",0)>=0.65 and r["delta"]>=0.15 and s.get("erap1",0)<0.35 and s.get("irap",0)<0.35 and s.get("anpep",0)<0.30
            for r in results for s in [r["screens"]]
        )
        if not any_pass:
            print("  No variant passes all 5 gates. Closest candidates above.")

    with open(os.path.join(RESULTS_DIR, "v32_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {os.path.join(RESULTS_DIR, 'v32_results.json')}")

if __name__ == "__main__":
    main()
