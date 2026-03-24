"""V4.3 Short-Stack + KILKLYSSKKY Validation Campaign.

5 experiments in one Vast.ai run, ALL scored against 4 targets:
  ERAP2 K392, ERAP2 N392, ERAP1, IRAP

Experiments:
  1. KILKLYSSKKY multi-seed confirmation (5 seeds x 4 targets = 20)
  2. KILKLYSSKKY P1 scan (8 P1s x 4 targets = 32)
  3. KKY motif tests: graft + knockout (5 designs x 4 targets = 20)
  4. Short-Stack 6-7mer library (13 designs x 4 targets = 52)
  5. v42_multi_TVL multi-seed confirmation (5 seeds x 4 targets = 20)

Total: ~144 Boltz-2 predictions (multi-seed counted as separate runs)

Usage:
    python v43_shortstack_campaign.py                    # Generate YAMLs
    python v43_shortstack_campaign.py --analyze          # Analyze results
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import json
import os
import statistics
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
YAML_DIR = SCRIPT_DIR / "v43_boltz2_inputs"
RESULTS_DIR = SCRIPT_DIR / "v43_boltz2_out"

DIFFUSION_SAMPLES = 5
BASE_SEED = 42
MULTI_SEEDS = [42, 123, 456, 789, 1024]

# ── Target sequences (same as V4.2) ──────────────────────────────────────────

ERAP2_K392 = (
    "MFHSSAMVNSHRKPMFNIHRGFYCLTAILPQICICSQFSVPSSYHFTEDPGAFPVATNGERFPWQELRLPS"
    "VVIPLHYDLFVHPNLTSLDFVASEKIEVLVSNATQFIILHSKDLEITNATLQSEEDSRYMKPGKELKVLSYP"
    "AHEQIALLVPEKLTPHLKYYVAMDFQAKLGDGFEGFYKSTYRTLGGETRILAVTDFEPTQARMAFPCFDEPLF"
    "KANFSIKIRRESRHIALSNMPKVKTIELEGGLLEDHFETTVKMSTYLVAYIVCDFHSLSGFTSSGVKVSIYASP"
    "DKRNQTHYALQASLKLLDFYEKYFDIYYPLSKLDLIAIPDFAPGAMENWGLITYRETSLLFDPKTSSASDKLW"
    "VTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPI"
    "SKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLSNSCLESDFTS"
    "GGVCHSDPKMTSNMLAFLGENAEVKEMMTTWTLQKGIPLLVVKQDGCSLRLQQERFLQGVFQEDPEWRALQER"
    "YLWHIPLTYSTSSSNVIHRHILKSKTDTLDLPEKTSWVKFNVDSNGYYIVHYEGHGWDQLITQLNQNHTLLRP"
    "KDRVGLIHDVFQLVGAGRLTLDKALDMTYYLQHETSSPALLEGLSYLESFYHMMDRRNISDISENLKRYLLQYF"
    "KPVIDRQSWSDKGSVWDRMLRSALLKLACDLNHAPCIQKAAELFSQWMESSGKLNIPTDVLKIVYSVGAQTTAG"
    "WNYLLEQYELSMSSAEQNKILYALSTSKHQEKLLKLIELGMEGKVIKTQNLAALLHAIARRPKGQQLAWDFVRE"
    "NWTHLLKKFDLGSYDIRMIISGTTAHFSSKDKLQEVKLFFESLEAQGSHLDIFQTVLETITKNIKWLEKNLPTLR"
    "TWLMVNT"
)
ERAP2_N392 = ERAP2_K392[:391] + "N" + ERAP2_K392[392:]

ERAP1 = (
    "MVFLPLKWSLATMSFLLSSLLALLTVSTPSWCQSTEASPKRSDGTPFPWNKIRLPEYVIPVHYDLLIHANLT"
    "TLTFWGTTKVEITASQPTSTIILHSHHLQISRATLRKGAGERLSEEEPLQVLEHPRQEQIALLAPEPLLVGLP"
    "YTVVIHYAGNLSETFHGFYKSTYRTKEGELRILASTQFEPTAARMAFPCFDEPAFKASFSIKIRREPRHLAISN"
    "MPLVKSVTVAEGLIEDHFDVTVKMSTYLVAFIISDFESVSKITKSGVKVSVYAVPDKINQADYALDAAVTLLEF"
    "YEDYFSIPYPLPKQDLAAIPDFQSGAMENWGLTTYRESALLFDAEKSSASSKLGITMTVAHELAHQWFGNLVTM"
    "EWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDK"
    "GACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSRSQHSSSSSHWHQEG"
    "VDVKTMMNTWTLQKGFPLITITVRGRNVHMKQEHYMKGSDGAPDTGYLWHVPLTFITSKSDMVHRFLLKTKTDV"
    "LILPEEVEWIKFNVGMNGYYIVHYEDDGWDSLTGLLKGTHTAVSSNDRASLINNAFQLVSIGKLSIEKALDLSLY"
    "LKHETEIMPVFQGLNELIPMYKLMEKRDMNEVETQFKAFLIRLLRDLIDKQTWTDEGSVSERMLRSQLLLLACVH"
    "NYQPCVQRAEGYFRKWKESNGNLSLPVDVTLAVFAVGAQSTEGHWDFLYSKYQFSLSSTEKSQIEFALCRTQNKE"
    "KLQWLLDESFKGDKIKTQEFPQILTLIGRNPVGYPLAWQFLRKNWNKLVQKFELGSSSIAHMVMGTTNQFSTRT"
    "RLEEVKGFFSSLKENGSSQLRCVQQTIETIEENIGWMDKNFDKIRVWLQSEKLERM"
)

IRAP = (
    "MEPFTNDRLQLPRNMIENSMFEEEPDVVDLAKEPCLHPLEPDEVEYEPRGSRLLVRGLGEHEMEEDEEDYESS"
    "AKLLGMSFMNRSSGLRNSATGYRQSPDGACSVPSARTMVVCAFVIVVAVSVIMVIYLLPRCTFTKEGCHKKNQS"
    "IGLIQPFATNGKLFPWAQIRLPTAVVPLRYELSSLHPNLTSMTFRGSVTISVQALQVTWNIILHSTGHNISRVT"
    "FMSAVSSQEKQAEILEYAYIHGQIAIVAPEALLAGHNYTLKIEYSANISSSYYGFYGFSYTDESNEKKYFAATQF"
    "EPLAARSSAFPCFDEPAFKATFIIKIIRDEQYTALSNMPKKSSVVLDDGLVQDEFSESVKMSTYLVAFIVGEMKNL"
    "SQAVNGTLVSIYAVPEKIGQVHYALETTVKLLEFFQNYFEIQYPLKKLDLVAIPDFEAGAMENWGLLTFREETML"
    "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKT"
    "MKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDS"
    "FNEVTNQTLDVKRMMKTWTLQKGFPLVTVQKKGKELFIQQERFFLNMKPEIQPSDTSYLWHIPLSYVTEGRNYSK"
    "YQSVSLLDKKSGVINLTEEVLWVKVNINMNGYYIVHYADDDWEALIHQLKINPYVLSDKDRANLINNIFELAGLGK"
    "VPLKRAFDLINYLGNENHTAPITEALFQTDLIYNLLEKLGYMDLASRLVTRVFKLLQNQIQQQTWTDEGTPSMREL"
    "RSALLEFACTHNLGNCSTTAMKLFDDWMASNGTQSLPTDVMTTVFKVGAKTDKGWSFLLGKYISIGSEAEKNKIL"
    "EALASSEDVRKLYWLMKSSLNGDNFRTQKLSFIIRTVGRHFPGHLLAWDFVKENWNKLVQKFPLGSYTIQNIVAG"
    "STYLFSTKTHLSEVQAFFENQSEATFRLRCVQEALEVIQLNIQWMEKNLKSLTWWL"
)

TARGETS = {
    "ERAP2_K392": ERAP2_K392,
    "ERAP2_N392": ERAP2_N392,
    "ERAP1": ERAP1,
    "IRAP": IRAP,
}

# ── Experiment 1: KILKLYSSKKY Multi-Seed ──────────────────────────────────────
# Single-sample showed +0.427 delta — need multi-seed to confirm

EXP1_MULTISEED = [
    {"name": "kilk_multiseed", "seq": "KILKLYSSKKY", "experiment": "E1_multiseed",
     "rationale": "KILKLYSSKKY multi-seed confirmation (single-sample: K392=0.892, delta=+0.427)"},
]

# ── Experiment 2: KILKLYSSKKY P1 Scan ────────────────────────────────────────
# Tests steric model: V/I/F should beat E/D if steric > electrostatic

EXP2_P1SCAN = [
    {"name": "kilk_p1_V", "seq": "VILKLYSSKKY", "experiment": "E2_p1scan", "p1": "V",
     "rationale": "P1=Val, branched hydrophobic — predicted best by steric model"},
    {"name": "kilk_p1_I", "seq": "IILKLYSSKKY", "experiment": "E2_p1scan", "p1": "I",
     "rationale": "P1=Ile, branched hydrophobic — should match Val if steric model holds"},
    {"name": "kilk_p1_E", "seq": "EILKLYSSKKY", "experiment": "E2_p1scan", "p1": "E",
     "rationale": "P1=Glu, salt bridge candidate — should underperform V/I if steric > electrostatic"},
    {"name": "kilk_p1_D", "seq": "DILKLYSSKKY", "experiment": "E2_p1scan", "p1": "D",
     "rationale": "P1=Asp, shorter salt bridge — contrast to E"},
    {"name": "kilk_p1_L", "seq": "LILKLYSSKKY", "experiment": "E2_p1scan", "p1": "L",
     "rationale": "P1=Leu, linear hydrophobic — should be slightly worse than branched V/I"},
    {"name": "kilk_p1_A", "seq": "AILKLYSSKKY", "experiment": "E2_p1scan", "p1": "A",
     "rationale": "P1=Ala, baseline neutral — this is the original PepMLM scaffold"},
    {"name": "kilk_p1_F", "seq": "FILKLYSSKKY", "experiment": "E2_p1scan", "p1": "F",
     "rationale": "P1=Phe, aromatic — can stack with K392 NH3+ tip (cation-pi)"},
    {"name": "kilk_p1_W", "seq": "WILKLYSSKKY", "experiment": "E2_p1scan", "p1": "W",
     "rationale": "P1=Trp, large aromatic — strongest cation-pi potential"},
]

# ── Experiment 3: KKY Motif Tests ─────────────────────────────────────────────
# Tests whether C-terminal KKY is a modular selectivity handle

# Top 3 Bucket A scaffolds (high K392 ipTM but low selectivity) for KKY grafting
# From N392 scaffold screen results:
#   TSYYWSPLHK  (K392=0.891, delta=+0.019) — replace HK tail with SSKKY
#   SVLLLLKFKY  (K392=0.881, delta=-0.041) — replace KFKY with SSKKY
#   VTLKLYSIFHG (K392=0.870, delta=-0.005) — replace IFHG with SSKKY

EXP3_KKY = [
    {"name": "kky_graft_TSYYW", "seq": "TSYYWSPSSKKY", "experiment": "E3_kky",
     "rationale": "TSYYWSPLHK Bucket A scaffold with SSKKY graft (replace LHK)"},
    {"name": "kky_graft_SVLLL", "seq": "SVLLLLSSKKY", "experiment": "E3_kky",
     "rationale": "SVLLLLKFKY Bucket A scaffold with SSKKY graft (replace KFKY)"},
    {"name": "kky_graft_VTLKL", "seq": "VTLKLYSSKKY", "experiment": "E3_kky",
     "rationale": "VTLKLYSIFHG Bucket A scaffold with SSKKY graft (replace IFHG)"},
    {"name": "kky_knockout_AAA", "seq": "KILKLYSSAAA", "experiment": "E3_kky",
     "rationale": "KILKLYSSKKY with KKY->AAA knockout — tests KKY necessity"},
    {"name": "kky_knockout_GGG", "seq": "KILKLYSSGGG", "experiment": "E3_kky",
     "rationale": "KILKLYSSKKY with KKY->GGG knockout — flexible linker control"},
]

# ── Experiment 4: Short-Stack 6-7mer Library ──────────────────────────────────
# Two scaffold families compressed to 6-7 residues
# Family A: TVL-derived from V4.2 winner EKTVLLSIGK (handles 392+398+403)
# Family B: KKY-derived from KILKLYSSKKY (KKY C-terminal motif)

EXP4_SHORTSTACK = [
    # Family A: TVL-derived (compress 392+398+403 handles into 6-7 residues)
    {"name": "ss_tvl_EKTVKY", "seq": "EKTVKY", "experiment": "E4_shortstack", "family": "A_TVL",
     "rationale": "6-mer core: P1=E(392) + T(398) + V(403) + KY(anchor)"},
    {"name": "ss_tvl_EKTVLKY", "seq": "EKTVLKY", "experiment": "E4_shortstack", "family": "A_TVL",
     "rationale": "7-mer: adds Leu spacer between V(403) and KY anchor"},
    {"name": "ss_tvl_DKTVKY", "seq": "DKTVKY", "experiment": "E4_shortstack", "family": "A_TVL",
     "rationale": "6-mer Asp variant: shorter salt bridge at P1"},
    {"name": "ss_tvl_VKTVKY", "seq": "VKTVKY", "experiment": "E4_shortstack", "family": "A_TVL",
     "rationale": "6-mer steric P1: Val packs K392 aliphatic stem"},
    {"name": "ss_tvl_IKTVKY", "seq": "IKTVKY", "experiment": "E4_shortstack", "family": "A_TVL",
     "rationale": "6-mer branched steric: Ile for maximum K392 stem packing"},

    # Family B: KKY-derived (preserve C-terminal KKY motif in 6-7 residues)
    {"name": "ss_kky_KILKKY", "seq": "KILKKY", "experiment": "E4_shortstack", "family": "B_KKY",
     "rationale": "6-mer: compress KILKLYSSKKY, keep KIL- head + KKY tail"},
    {"name": "ss_kky_KILKKYA", "seq": "KILKKYA", "experiment": "E4_shortstack", "family": "B_KKY",
     "rationale": "7-mer: adds Ala spacer for flexibility"},
    {"name": "ss_kky_VILKKY", "seq": "VILKKY", "experiment": "E4_shortstack", "family": "B_KKY",
     "rationale": "6-mer steric P1: Val replaces Lys at P1"},
    {"name": "ss_kky_IILKKY", "seq": "IILKKY", "experiment": "E4_shortstack", "family": "B_KKY",
     "rationale": "6-mer double branched: Ile at P1 + Ile at P2"},
    {"name": "ss_kky_EILKKY", "seq": "EILKKY", "experiment": "E4_shortstack", "family": "B_KKY",
     "rationale": "6-mer charged P1: Glu salt bridge + KKY motif"},

    # Controls
    {"name": "ss_ctrl_polyA", "seq": "AAAAAK", "experiment": "E4_shortstack", "family": "ctrl",
     "rationale": "6-mer poly-Ala negative control — no specific interactions"},
    {"name": "ss_ctrl_scramble", "seq": "KVTKEY", "experiment": "E4_shortstack", "family": "ctrl",
     "rationale": "6-mer scrambled EKTVKY — same composition, wrong order"},
    {"name": "ss_ctrl_trunc", "seq": "EKLLLL", "experiment": "E4_shortstack", "family": "ctrl",
     "rationale": "6-mer truncated VKLLLL scaffold — no handles, no anchor"},
]

# ── Experiment 5: v42_multi_TVL Multi-Seed ────────────────────────────────────
# Best V4.2 candidate needs multi-seed validation

EXP5_TVL_MULTISEED = [
    {"name": "tvl_multiseed", "seq": "EKTVLLSIGK", "experiment": "E5_multiseed",
     "rationale": "v42_multi_TVL multi-seed: V4.2 winner (K392=0.721, selE1=+0.150, selIR=+0.154)"},
]


def get_all_experiments():
    """Return all experiments organized by type."""
    return {
        "E1_multiseed": EXP1_MULTISEED,
        "E2_p1scan": EXP2_P1SCAN,
        "E3_kky": EXP3_KKY,
        "E4_shortstack": EXP4_SHORTSTACK,
        "E5_multiseed": EXP5_TVL_MULTISEED,
    }


def get_all_peptides():
    """Flat list of all unique peptides."""
    experiments = get_all_experiments()
    seen = set()
    peptides = []
    for exp_list in experiments.values():
        for p in exp_list:
            if p["name"] not in seen:
                seen.add(p["name"])
                peptides.append(p)
    return peptides


def write_boltz2_yamls(output_dir):
    """Generate Boltz-2 YAML input files for all experiments."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    experiments = get_all_experiments()

    for exp_name, pep_list in experiments.items():
        is_multiseed = exp_name in ("E1_multiseed", "E5_multiseed")

        for pep in pep_list:
            for target_name, target_seq in TARGETS.items():
                if is_multiseed:
                    # Multi-seed: one YAML per seed
                    for seed_idx, seed in enumerate(MULTI_SEEDS):
                        yaml_name = f"{pep['name']}_s{seed}_{target_name}.yaml"
                        yaml_content = (
                            f"sequences:\n"
                            f"  - protein:\n"
                            f"      id: A\n"
                            f"      sequence: {target_seq}\n"
                            f"      msa: empty\n"
                            f"  - protein:\n"
                            f"      id: B\n"
                            f"      sequence: {pep['seq']}\n"
                            f"      msa: empty\n"
                        )
                        (output_dir / yaml_name).write_text(yaml_content)
                        count += 1
                else:
                    # Standard: single YAML per peptide-target pair
                    yaml_name = f"{pep['name']}_{target_name}.yaml"
                    yaml_content = (
                        f"sequences:\n"
                        f"  - protein:\n"
                        f"      id: A\n"
                        f"      sequence: {target_seq}\n"
                        f"      msa: empty\n"
                        f"  - protein:\n"
                        f"      id: B\n"
                        f"      sequence: {pep['seq']}\n"
                        f"      msa: empty\n"
                    )
                    (output_dir / yaml_name).write_text(yaml_content)
                    count += 1

    return count


def print_experiment_summary():
    experiments = get_all_experiments()
    total = 0
    print("=" * 100)
    print("V4.3 SHORT-STACK + KILKLYSSKKY VALIDATION CAMPAIGN")
    print("=" * 100)
    print()

    for exp_name, pep_list in experiments.items():
        is_multiseed = exp_name in ("E1_multiseed", "E5_multiseed")
        n_targets = len(TARGETS)

        if is_multiseed:
            n_runs = len(pep_list) * n_targets * len(MULTI_SEEDS)
        else:
            n_runs = len(pep_list) * n_targets

        labels = {
            "E1_multiseed": "Exp 1: KILKLYSSKKY Multi-Seed Confirmation",
            "E2_p1scan": "Exp 2: KILKLYSSKKY P1 Scan (steric model test)",
            "E3_kky": "Exp 3: KKY Motif Graft + Knockout",
            "E4_shortstack": "Exp 4: Short-Stack 6-7mer Library",
            "E5_multiseed": "Exp 5: v42_multi_TVL Multi-Seed Confirmation",
        }
        print(f"### {labels[exp_name]}")
        print(f"    {len(pep_list)} peptides x {n_targets} targets"
              + (f" x {len(MULTI_SEEDS)} seeds" if is_multiseed else "")
              + f" = {n_runs} runs")

        for p in pep_list:
            print(f"    {p['name']:>25}  {p['seq']:>15}  ({len(p['seq'])}-mer)  {p['rationale']}")
        print()
        total += n_runs

    print(f"TOTAL: {total} Boltz-2 predictions")
    print(f"Estimated time on H100: ~{total}min  |  RTX 4090: ~{total * 1.5:.0f}min")
    print(f"Estimated cost: H100 ~${total/60*3.5:.2f}  |  RTX 4090 ~${total/60*0.30:.2f}")
    return total


# ── Analysis ──────────────────────────────────────────────────────────────────

def collect_scores(results_dir):
    """Collect ipTM scores from Boltz-2 confidence JSON files."""
    results_dir = Path(results_dir)
    scores = {}
    for json_file in sorted(results_dir.rglob("confidence_*.json")):
        try:
            data = json.loads(json_file.read_text())
            iptm = data.get("iptm", data.get("ipTM"))
            if iptm is None:
                continue
            run_name = None
            for parent in json_file.parents:
                if parent.parent == results_dir or parent == results_dir:
                    run_name = parent.name
                    break
            if run_name is None:
                continue
            if run_name not in scores:
                scores[run_name] = []
            scores[run_name].append(float(iptm))
        except (json.JSONDecodeError, ValueError):
            continue
    return scores


def parse_run_name(run_name):
    """Parse a run name into (peptide_name, target_name, seed_or_None)."""
    # Multi-seed format: kilk_multiseed_s42_ERAP2_K392
    # Standard format: kilk_p1_V_ERAP2_K392
    for target_name in TARGETS:
        suffix = f"_{target_name}"
        if run_name.endswith(suffix):
            prefix = run_name[:-len(suffix)]
            # Check for seed suffix: _s42, _s123, etc.
            seed = None
            for s in MULTI_SEEDS:
                seed_suffix = f"_s{s}"
                if prefix.endswith(seed_suffix):
                    pep_name = prefix[:-len(seed_suffix)]
                    seed = s
                    return pep_name, target_name, seed
            return prefix, target_name, None
    return None, None, None


def analyze_results(results_dir):
    """Analyze V4.3 campaign results."""
    results_dir = Path(results_dir)
    scores = collect_scores(results_dir)

    if not scores:
        print("No confidence scores found in", results_dir)
        return

    pep_lookup = {p["name"]: p for p in get_all_peptides()}

    # ── Aggregate scores per peptide per target ──────────────────────────
    # For multi-seed: combine all seeds
    # For standard: use diffusion sample average
    agg = {}  # pep_name -> target -> [iptm values]
    for run_name, iptm_list in scores.items():
        pep_name, target_name, seed = parse_run_name(run_name)
        if pep_name is None:
            continue
        if pep_name not in agg:
            agg[pep_name] = {}
        if target_name not in agg[pep_name]:
            agg[pep_name][target_name] = []
        # For multi-seed, each seed has multiple diffusion samples
        # Take the mean of diffusion samples per seed, then we'll average across seeds
        agg[pep_name][target_name].append(statistics.mean(iptm_list))

    # ── Build results rows ───────────────────────────────────────────────
    rows = []
    for pep_name, target_scores in agg.items():
        meta = pep_lookup.get(pep_name, {"seq": "?", "experiment": "?"})

        k392_vals = target_scores.get("ERAP2_K392", [])
        n392_vals = target_scores.get("ERAP2_N392", [])
        erap1_vals = target_scores.get("ERAP1", [])
        irap_vals = target_scores.get("IRAP", [])

        k392 = statistics.mean(k392_vals) if k392_vals else 0
        n392 = statistics.mean(n392_vals) if n392_vals else 0
        erap1 = statistics.mean(erap1_vals) if erap1_vals else 0
        irap = statistics.mean(irap_vals) if irap_vals else 0

        k392_sd = statistics.stdev(k392_vals) if len(k392_vals) > 1 else 0
        n392_sd = statistics.stdev(n392_vals) if len(n392_vals) > 1 else 0
        erap1_sd = statistics.stdev(erap1_vals) if len(erap1_vals) > 1 else 0
        irap_sd = statistics.stdev(irap_vals) if len(irap_vals) > 1 else 0

        delta_kn = k392 - n392
        erap2_best = max(k392, n392)
        sel_erap1 = erap2_best - erap1
        sel_irap = erap2_best - irap
        composite = erap2_best + sel_erap1 + sel_irap

        rows.append({
            "name": pep_name,
            "seq": meta.get("seq", "?"),
            "length": len(meta.get("seq", "")),
            "experiment": meta.get("experiment", "?"),
            "family": meta.get("family", ""),
            "p1": meta.get("p1", meta.get("seq", "?")[0] if meta.get("seq") else "?"),
            "k392": k392, "k392_sd": k392_sd, "n_k392": len(k392_vals),
            "n392": n392, "n392_sd": n392_sd,
            "erap1": erap1, "erap1_sd": erap1_sd,
            "irap": irap, "irap_sd": irap_sd,
            "delta_kn": delta_kn,
            "erap2_best": erap2_best,
            "sel_erap1": sel_erap1,
            "sel_irap": sel_irap,
            "composite": composite,
            "rationale": meta.get("rationale", ""),
        })

    # ── Print results by experiment ──────────────────────────────────────

    exp_labels = {
        "E1_multiseed": "EXPERIMENT 1: KILKLYSSKKY Multi-Seed Confirmation",
        "E2_p1scan": "EXPERIMENT 2: KILKLYSSKKY P1 Scan (Steric Model Test)",
        "E3_kky": "EXPERIMENT 3: KKY Motif Graft + Knockout",
        "E4_shortstack": "EXPERIMENT 4: Short-Stack 6-7mer Library",
        "E5_multiseed": "EXPERIMENT 5: v42_multi_TVL Multi-Seed Confirmation",
    }

    for exp_name, label in exp_labels.items():
        exp_rows = [r for r in rows if r["experiment"] == exp_name]
        if not exp_rows:
            continue

        print()
        print("=" * 140)
        print(label)
        print("=" * 140)

        is_multiseed = exp_name in ("E1_multiseed", "E5_multiseed")

        if is_multiseed:
            for r in exp_rows:
                print(f"\n  {r['name']} ({r['seq']}, {r['length']}-mer)")
                print(f"    K392: {r['k392']:.3f} +/- {r['k392_sd']:.3f}  (n={r['n_k392']})")
                print(f"    N392: {r['n392']:.3f} +/- {r['n392_sd']:.3f}")
                print(f"    ERAP1: {r['erap1']:.3f} +/- {r['erap1_sd']:.3f}")
                print(f"    IRAP:  {r['irap']:.3f} +/- {r['irap_sd']:.3f}")
                print(f"    Delta K-N: {r['delta_kn']:+.3f}")
                print(f"    Selectivity: vs ERAP1 {r['sel_erap1']:+.3f}  |  vs IRAP {r['sel_irap']:+.3f}")
                print(f"    Composite: {r['composite']:.3f}")
        else:
            # Sort by composite score
            exp_rows.sort(key=lambda r: r["composite"], reverse=True)
            header = (f"{'Rank':>4}  {'Name':>25}  {'Seq':>15}  {'Len':>3}  "
                      f"{'K392':>6}  {'N392':>6}  {'dKN':>6}  "
                      f"{'ERAP1':>6}  {'IRAP':>6}  {'selE1':>6}  {'selIR':>6}  {'Score':>6}")
            print(header)
            print("-" * len(header))
            for i, r in enumerate(exp_rows, 1):
                print(f"{i:>4}  {r['name']:>25}  {r['seq']:>15}  {r['length']:>3}  "
                      f"{r['k392']:.3f}  {r['n392']:.3f}  {r['delta_kn']:>+.3f}  "
                      f"{r['erap1']:.3f}  {r['irap']:.3f}  "
                      f"{r['sel_erap1']:>+.3f}  {r['sel_irap']:>+.3f}  {r['composite']:.3f}")

    # ── Key Questions ────────────────────────────────────────────────────
    print()
    print("=" * 100)
    print("KEY QUESTIONS")
    print("=" * 100)

    # Q1: Does KILKLYSSKKY hold up?
    kilk = next((r for r in rows if r["name"] == "kilk_multiseed"), None)
    if kilk:
        hold = "YES" if kilk["delta_kn"] > 0.15 else ("MARGINAL" if kilk["delta_kn"] > 0.05 else "NO")
        print(f"\n  Q1: KILKLYSSKKY multi-seed delta = {kilk['delta_kn']:+.3f} (SD={kilk['k392_sd']:.3f}) -> {hold}")
        print(f"      ERAP1 sel = {kilk['sel_erap1']:+.3f}  |  IRAP sel = {kilk['sel_irap']:+.3f}")

    # Q2: Steric model on KILKLYSSKKY?
    p1_rows = [r for r in rows if r["experiment"] == "E2_p1scan"]
    if p1_rows:
        p1_rows.sort(key=lambda r: r["delta_kn"], reverse=True)
        print(f"\n  Q2: P1 scan selectivity ranking (delta K-N):")
        for r in p1_rows:
            print(f"      {r['p1']:>3} ({r['seq'][:1]}): dKN={r['delta_kn']:+.3f}  K392={r['k392']:.3f}  "
                  f"ERAP1={r['erap1']:.3f}  IRAP={r['irap']:.3f}")
        hydrophobic_avg = statistics.mean([r["delta_kn"] for r in p1_rows if r.get("p1") in ("V", "I", "F")])
        charged_avg = statistics.mean([r["delta_kn"] for r in p1_rows if r.get("p1") in ("E", "D")])
        print(f"      Hydrophobic (V/I/F) avg delta: {hydrophobic_avg:+.3f}")
        print(f"      Charged (E/D) avg delta: {charged_avg:+.3f}")
        steric_confirmed = "YES" if hydrophobic_avg > charged_avg else "NO"
        print(f"      Steric > Electrostatic on KILKLYSSKKY? {steric_confirmed}")

    # Q3: KKY motif necessary?
    kky_ko = [r for r in rows if "knockout" in r["name"]]
    kky_parent = next((r for r in rows if r["name"] == "kilk_multiseed"), None)
    if kky_ko and kky_parent:
        print(f"\n  Q3: KKY knockout effect:")
        print(f"      Parent (KILKLYSSKKY): dKN={kky_parent['delta_kn']:+.3f}  IRAP={kky_parent['irap']:.3f}")
        for r in kky_ko:
            print(f"      {r['name']}: dKN={r['delta_kn']:+.3f}  IRAP={r['irap']:.3f}  "
                  f"(delta drop: {r['delta_kn'] - kky_parent['delta_kn']:+.3f})")

    # Q4: Can 6-mers bind ERAP2?
    ss_rows = [r for r in rows if r["experiment"] == "E4_shortstack" and r["family"] != "ctrl"]
    if ss_rows:
        best_6mer = max([r for r in ss_rows if r["length"] <= 7], key=lambda r: r["erap2_best"], default=None)
        print(f"\n  Q4: Best 6-7mer ERAP2 binding:")
        if best_6mer:
            binds = "YES" if best_6mer["erap2_best"] > 0.5 else "NO"
            print(f"      {best_6mer['name']} ({best_6mer['seq']}): ipTM={best_6mer['erap2_best']:.3f} -> {binds}")
            print(f"      ERAP1={best_6mer['erap1']:.3f}  IRAP={best_6mer['irap']:.3f}")

    # Q5: Do 6-mers evade ERAP1?
    if ss_rows:
        print(f"\n  Q5: 6-7mer ERAP1 evasion:")
        for r in sorted(ss_rows, key=lambda r: r["erap1"]):
            evades = "EVADES" if r["erap1"] < 0.4 else ("MARGINAL" if r["erap1"] < 0.6 else "BINDS")
            print(f"      {r['name']:>20} ({r['length']}-mer): ERAP1={r['erap1']:.3f} IRAP={r['irap']:.3f} -> {evades}")

    # ── Save results ─────────────────────────────────────────────────────
    output = {
        "experiment": "V4.3 Short-Stack + KILKLYSSKKY Validation Campaign",
        "n_peptides": len(get_all_peptides()),
        "targets": list(TARGETS.keys()),
        "diffusion_samples": DIFFUSION_SAMPLES,
        "multi_seeds": MULTI_SEEDS,
        "results": rows,
    }
    output_path = results_dir / "v43_analysis.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved analysis to {output_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="V4.3 Short-Stack + KILKLYSSKKY Campaign")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--results-dir", type=str, default=str(RESULTS_DIR))
    args = parser.parse_args()

    if args.analyze:
        analyze_results(Path(args.results_dir))
        return

    total = print_experiment_summary()
    print()

    n = write_boltz2_yamls(YAML_DIR)
    print(f"\nWrote {n} Boltz-2 input YAMLs to {YAML_DIR}/")

    # Save experiment definition
    defn_path = SCRIPT_DIR / "v43_experiment_definition.json"
    with open(defn_path, "w") as f:
        json.dump(get_all_peptides(), f, indent=2)
    print(f"Saved experiment definition to {defn_path}")


if __name__ == "__main__":
    main()
