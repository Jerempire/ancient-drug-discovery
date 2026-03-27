"""
ProteinMPNN N392 redesign pipeline — runs entirely on Vast.ai.
1. Convert N392 complex CIF to PDB
2. Find binder positions near Y398/D414/392 (could extend contact network)
3. Run ProteinMPNN with those positions redesignable
4. Build Boltz-2 YAMLs for top 10 designs x 4 targets
5. Run Boltz-2 screen + analyze

Usage: python3 /workspace/n392_mpnn_pipeline.py
"""
import json, os, subprocess, sys, glob
import gemmi
from collections import defaultdict

CIF_PATH = "/workspace/results/wt_erap2_n392/boltz_results_wt_erap2_n392/predictions/wt_erap2_n392/wt_erap2_n392_model_0.cif"
PDB_PATH = "/workspace/n392_complex.pdb"
MPNN_DIR = "/workspace/ProteinMPNN"
OUTPUT_DIR = "/workspace/mpnn_n392"
YAML_DIR = "/workspace/mpnn_yamls"
RESULTS_DIR = "/workspace/mpnn_results"
N_DESIGNS = 100
TEMPERATURE = 0.1
TOP_K = 10

BINDER_WT = "DIRHYFKSLEEYLKNLPKVVDMLVDLYSKGIFHLDNTNILVKDDKFYAIDFGSAYINEKKSTDYTLKIKNDQISSEEYVKSVSEKIANALKN"

ERAP2_FULL = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNG"
    "ERFPWQELRLPSVVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITN"
    "ATLQSEEDSRYMKPGKELKVLSYPAHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEG"
    "FYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLFKANFSIKIRRESRHIALSNMPKV"
    "KTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASPDKRNQTHYA"
    "LQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSAS"
    "DKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFL"
    "NVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGI"
    "IQYLKKFSYRNAKNDDLWSSLSNSCLESDFTSGGVCHSDPKMTSNMLAFLGENAEVKEMM"
    "TTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQERYLWHIPLTYSTSSS"
    "NVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNIS"
    "DISENLKRYLLQYFKPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQ"
    "WMESSGKLNIPTDVLKIVYSVGAQTTAGWNYLLEQYELSMSSAEQNKILYALSTSKHQEK"
    "LLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRENWTHLLKKFDLGSYDIRM"
    "IISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLRTWL"
    "MVNT"
)
ERAP2_K392 = ERAP2_FULL[349:500]
ERAP2_N392 = ERAP2_FULL[349:391] + "N" + ERAP2_FULL[392:500]
ERAP1 = "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR"
IRAP_SEQ = "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN"
TARGETS = {"erap2_k392": ERAP2_K392, "erap2_n392": ERAP2_N392, "erap1": ERAP1, "irap": IRAP_SEQ}


def step1_convert():
    print("[1/5] Converting CIF to PDB...")
    st = gemmi.read_structure(CIF_PATH)
    st.write_pdb(PDB_PATH)
    model = st[0]
    for chain in model:
        res_count = sum(1 for r in chain if r.name != "UNK")
        print(f"  Chain {chain.name}: {res_count} residues")
    return PDB_PATH


def step2_find_redesignable(pdb_path):
    print("[2/5] Finding redesignable positions...")
    st = gemmi.read_structure(pdb_path)
    model = st[0]
    chains = {}
    for chain in model:
        res_list = [r for r in chain if r.name != "UNK"]
        chains[chain.name] = len(res_list)
    sorted_chains = sorted(chains.items(), key=lambda x: -x[1])
    target_chain_id = sorted_chains[0][0]
    binder_chain_id = sorted_chains[1][0]
    print(f"  Target: chain {target_chain_id} ({chains[target_chain_id]} res)")
    print(f"  Binder: chain {binder_chain_id} ({chains[binder_chain_id]} res)")

    target_chain = model.find_chain(target_chain_id)
    binder_chain = model.find_chain(binder_chain_id)

    # Key ERAP2 residues (cropped numbering: full - 349)
    key_resnums = {43: "N392", 49: "Y398", 54: "A403", 65: "D414"}
    target_coords = {}
    for res in target_chain:
        rnum = res.seqid.num
        if rnum in key_resnums:
            atoms = [(a.name, a.pos) for a in res if a.element.name != "H"]
            if atoms:
                target_coords[rnum] = atoms

    binder_positions = []
    for res in binder_chain:
        rnum = res.seqid.num
        b_atoms = [(a.name, a.pos) for a in res if a.element.name != "H"]
        if not b_atoms:
            continue
        dists = {}
        for tnum, t_atoms in target_coords.items():
            min_d = float("inf")
            for _, bp in b_atoms:
                for _, tp in t_atoms:
                    d = bp.dist(tp)
                    if d < min_d:
                        min_d = d
            dists[key_resnums[tnum]] = min_d
        binder_positions.append({
            "resnum": rnum, "resname": res.name,
            "dist_N392": dists.get("N392", 99), "dist_Y398": dists.get("Y398", 99),
            "dist_A403": dists.get("A403", 99), "dist_D414": dists.get("D414", 99),
        })

    # Redesignable: within 10A of Y398 or D414, or within 8A of N392
    redesignable = []
    for bp in binder_positions:
        if bp["dist_Y398"] < 10.0 or bp["dist_D414"] < 10.0 or bp["dist_N392"] < 8.0:
            redesignable.append(bp["resnum"])

    print(f"  Redesignable: {len(redesignable)} positions")
    for bp in binder_positions:
        if bp["resnum"] in redesignable:
            print(f"    Pos {bp['resnum']:3d} ({bp['resname']:3s}): "
                  f"N392={bp['dist_N392']:.1f}A Y398={bp['dist_Y398']:.1f}A "
                  f"A403={bp['dist_A403']:.1f}A D414={bp['dist_D414']:.1f}A")

    return binder_chain_id, target_chain_id, redesignable


def step3_run_mpnn(pdb_path, binder_chain, target_chain, redesignable):
    print(f"[3/5] Running ProteinMPNN ({N_DESIGNS} designs, T={TEMPERATURE})...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Install if needed
    if not os.path.exists(f"{MPNN_DIR}/protein_mpnn_run.py"):
        print("  Installing ProteinMPNN...")
        subprocess.run(["git", "clone", "--depth", "1",
                       "https://github.com/dauparas/ProteinMPNN.git", MPNN_DIR],
                      capture_output=True, timeout=60)

    st = gemmi.read_structure(pdb_path)
    model = st[0]
    bc = model.find_chain(binder_chain)
    all_binder_res = [r.seqid.num for r in bc if r.name != "UNK"]
    fixed_binder = [r for r in all_binder_res if r not in redesignable]

    # Write JSONL files (ProteinMPNN format)
    jsonl_chains = f"{OUTPUT_DIR}/chains.jsonl"
    with open(jsonl_chains, "w") as f:
        f.write(json.dumps({os.path.basename(pdb_path).replace(".pdb", ""): [binder_chain]}) + "\n")

    jsonl_fixed = f"{OUTPUT_DIR}/fixed.jsonl"
    with open(jsonl_fixed, "w") as f:
        fixed_str = " ".join(str(r) for r in fixed_binder)
        f.write(json.dumps({os.path.basename(pdb_path).replace(".pdb", ""): {binder_chain: fixed_str}}) + "\n")

    print(f"  Redesignable: {len(redesignable)}, Fixed: {len(fixed_binder)}")

    cmd = [
        sys.executable, f"{MPNN_DIR}/protein_mpnn_run.py",
        "--out_folder", OUTPUT_DIR,
        "--num_seq_per_target", str(N_DESIGNS),
        "--sampling_temp", str(TEMPERATURE),
        "--batch_size", "1",
        "--pdb_path", pdb_path,
        "--chain_id_jsonl", jsonl_chains,
        "--fixed_positions_jsonl", jsonl_fixed,
    ]
    print(f"  Running: {' '.join(cmd[-6:])}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    print(f"  Return code: {result.returncode}")
    if result.stdout:
        print(f"  STDOUT (last 300): {result.stdout[-300:]}")
    if result.returncode != 0:
        print(f"  STDERR (last 500): {result.stderr[-500:]}")

    # Parse output
    fasta_files = glob.glob(f"{OUTPUT_DIR}/seqs/*.fa") + glob.glob(f"{OUTPUT_DIR}/seqs/*.fasta")
    if not fasta_files:
        print(f"  No FASTA output! Contents of {OUTPUT_DIR}:")
        for root, dirs, files in os.walk(OUTPUT_DIR):
            for fn in files:
                print(f"    {os.path.join(root, fn)}")
        return []

    designs = []
    for fa in fasta_files:
        with open(fa) as f:
            lines = f.readlines()
        for i in range(0, len(lines), 2):
            if i >= len(lines) or not lines[i].startswith(">"):
                continue
            header = lines[i].strip()
            seq_line = lines[i+1].strip() if i+1 < len(lines) else ""

            # Skip the first entry (native sequence)
            if "T=0.0000" in header or ", score=0.0000," in header:
                continue

            # ProteinMPNN concatenates chains with /
            parts = seq_line.split("/")
            binder_seq = None
            for p in parts:
                if len(p) == 92:
                    binder_seq = p
                    break
            if not binder_seq:
                # Try the part matching binder length
                for p in parts:
                    if len(p) == len(BINDER_WT):
                        binder_seq = p
                        break
            if not binder_seq:
                continue

            score = 0
            for tag in ["score=", "global_score="]:
                if tag in header:
                    try:
                        score = float(header.split(tag)[1].split(",")[0].split()[0])
                    except:
                        pass
                    break

            n_mutations = sum(1 for a, b in zip(BINDER_WT, binder_seq) if a != b)
            designs.append({"header": header, "sequence": binder_seq, "score": score, "n_mutations": n_mutations})

    designs.sort(key=lambda x: x["score"])
    top = designs[:TOP_K]
    print(f"  {len(designs)} valid designs, taking top {len(top)}")
    for i, d in enumerate(top):
        muts = []
        for j, (a, b) in enumerate(zip(BINDER_WT, d["sequence"])):
            if a != b:
                muts.append(f"{a}{j+1}{b}")
        mut_str = ", ".join(muts[:6])
        if len(muts) > 6:
            mut_str += f"... ({len(muts)} total)"
        print(f"    mpnn{i+1:02d}: {d['n_mutations']} muts, score={d['score']:.3f} | {mut_str}")
    return top


def step4_build_yamls(designs):
    print(f"[4/5] Building Boltz-2 YAMLs...")
    os.makedirs(YAML_DIR, exist_ok=True)
    count = 0
    for i, d in enumerate(designs):
        for tgt_name, tgt_seq in TARGETS.items():
            name = f"mpnn{i+1:02d}_{tgt_name}"
            with open(f"{YAML_DIR}/{name}.yaml", "w") as f:
                f.write(f"version: 1\nsequences:\n"
                        f"  - protein:\n      id: A\n      sequence: {tgt_seq}\n      msa: empty\n"
                        f"  - protein:\n      id: B\n      sequence: {d['sequence']}\n      msa: empty\n")
            count += 1
    print(f"  {count} YAMLs -> {YAML_DIR}")
    return count


def step5_run_boltz():
    print(f"[5/5] Running Boltz-2 screen...")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    yamls = sorted(glob.glob(f"{YAML_DIR}/*.yaml"))
    total = len(yamls)
    done = 0
    for yi, yaml_file in enumerate(yamls):
        name = os.path.basename(yaml_file).replace(".yaml", "")
        out_dir = f"{RESULTS_DIR}/{name}"
        if os.path.exists(out_dir) and glob.glob(f"{out_dir}/**/confidence_*.json", recursive=True):
            done += 1
            continue
        sys.stdout.write(f"\r  [{yi+1}/{total}] {name}...")
        sys.stdout.flush()
        result = subprocess.run(
            ["boltz", "predict", yaml_file, "--out_dir", out_dir, "--diffusion_samples", "1"],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode == 0:
            done += 1
        else:
            print(f"\n  FAILED: {name}")
    print(f"\n  Completed: {done}/{total}")


def analyze():
    print("\n" + "=" * 70)
    print("MPNN N392 REDESIGN RESULTS")
    print("=" * 70)
    scores = {}
    for d in sorted(glob.glob(f"{RESULTS_DIR}/*/")):
        name = os.path.basename(d.rstrip("/"))
        files = glob.glob(d + "**/confidence_*.json", recursive=True)
        if files:
            with open(files[0]) as f:
                c = json.load(f)
            scores[name] = c.get("iptm", c.get("i_ptm", 0))

    constructs = defaultdict(dict)
    for name, iptm in scores.items():
        for suffix in ("erap2_k392", "erap2_n392", "erap1", "irap"):
            if name.endswith(suffix):
                mut = name[:-(len(suffix)+1)]
                constructs[mut][suffix] = iptm
                break

    print("%-10s %6s %6s %7s %6s %6s  %s" % ("Design", "K392", "N392", "Delta", "ERAP1", "IRAP", ""))
    print("-" * 65)
    print("%-10s %6.3f %6.3f %+7.3f %6.3f %6.3f  %s" % ("wt", 0.800, 0.733, -0.067, 0.523, 0.416, "BASELINE"))

    for mut in sorted(constructs.keys()):
        t = constructs[mut]
        k = t.get("erap2_k392", 0)
        n = t.get("erap2_n392", 0)
        e = t.get("erap1", 0)
        i = t.get("irap", 0)
        delta = n - k if k and n else 0
        flag = ""
        if delta >= 0.10 and n >= 0.70 and e <= 0.35 and i <= 0.50:
            flag = "** HIT **"
        elif delta >= 0.10 and n >= 0.70:
            flag = "near"
        elif delta >= 0.05 and n >= 0.65:
            flag = "promising"
        print("%-10s %6.3f %6.3f %+7.3f %6.3f %6.3f  %s" % (mut, k, n, delta, e, i, flag))

    with open(f"{RESULTS_DIR}/mpnn_n392_summary.json", "w") as f:
        json.dump({"constructs": dict(constructs), "scores": scores}, f, indent=2)
    print(f"\nSaved to {RESULTS_DIR}/mpnn_n392_summary.json")


if __name__ == "__main__":
    pdb = step1_convert()
    binder_chain, target_chain, redesignable = step2_find_redesignable(pdb)
    designs = step3_run_mpnn(pdb, binder_chain, target_chain, redesignable)
    if designs:
        step4_build_yamls(designs)
        step5_run_boltz()
        analyze()
    else:
        print("\nNo designs generated. Check ProteinMPNN output above.")
