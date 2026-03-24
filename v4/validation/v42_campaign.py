"""V4.2 Rational Peptide Design Campaign: Multi-handle ERAP2 channel selectivity.

Designs substrate-mimetic peptides that exploit multiple ERAP2-specific channel
features beyond just K392, scoring against ERAP2 K392, ERAP2 N392, ERAP1, and IRAP.

ERAP2-unique channel residues exploited:
  392: K vs N (variant)       — P1 position
  398: Y vs F (ERAP1/IRAP)    — P3 position (H-bond to hydroxyl)
  403: A vs S (ERAP1/IRAP)    — P4 position (small hydrophobic fit)
  406: A vs V/K (ERAP1/IRAP)  — P5 position (anti-IRAP charge clash)
  412: Q vs K/S (ERAP1/IRAP)  — P6 position
  414: D vs G/Y (ERAP1/IRAP)  — P7 position (salt-bridge)

Scaffold: VKLLLL family (EKLLLLSIGK, 10-mer)
Families: P1 scan, Y398 handle, anti-IRAP, multi-handle

Usage:
    python v42_campaign.py                    # Generate YAMLs
    python v42_campaign.py --analyze          # Analyze results
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
YAML_DIR = SCRIPT_DIR / "v42_boltz2_inputs"
RESULTS_DIR = SCRIPT_DIR / "v42_boltz2_out"

DIFFUSION_SAMPLES = 5
SEED = 42

# ── Target sequences ─────────────────────────────────────────────────────────

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

# ── Peptide Design ────────────────────────────────────────────────────────────
# Base scaffold: EKLLLLSIGK (10-mer)
# Positions: P1-P2-P3-P4-P5-P6-P7-P8-P9-P10
#             E  K  L  L  L  L  S  I  G  K
# Channel mapping (approximate):
#   P1  → faces 392 (catalytic site)
#   P3  → faces 398 region
#   P4  → faces 403 region
#   P5  → faces 406 region
#   P6  → faces 412 region
#   P7  → faces 414 region

PEPTIDES = [
    # --- Family 1: P1 scan on base scaffold (baseline) ---
    {"name": "v42_base_E",      "seq": "EKLLLLSIGK", "family": "F1_P1scan", "p1": "E",
     "handles": "392 only",
     "rationale": "Original VKLLLL scaffold, P1=E (K392-leaning from V4)"},
    {"name": "v42_base_V",      "seq": "VKLLLLSIGK", "family": "F1_P1scan", "p1": "V",
     "handles": "392 only",
     "rationale": "P1=V, hydrophobic, K392-leaning from published data"},
    {"name": "v42_base_D",      "seq": "DKLLLLSIGK", "family": "F1_P1scan", "p1": "D",
     "handles": "392 only",
     "rationale": "P1=D, negative charge, contrast to E"},
    {"name": "v42_base_A",      "seq": "AKLLLLSIGK", "family": "F1_P1scan", "p1": "A",
     "handles": "392 only",
     "rationale": "P1=A, small/neutral, N392-leaning control"},
    {"name": "v42_base_L",      "seq": "LKLLLLSIGK", "family": "F1_P1scan", "p1": "L",
     "handles": "392 only",
     "rationale": "P1=L, hydrophobic, N392-leaning control"},

    # --- Family 2: Y398 handle (P3 substitution) ---
    # ERAP2 Y398 has hydroxyl for H-bonding; ERAP1/IRAP F398 cannot
    {"name": "v42_y398_S",      "seq": "EKSLLLSIGK", "family": "F2_Y398", "p1": "E",
     "handles": "392+398",
     "rationale": "P3=S, hydroxyl H-bonds Y398; penalized by F398 in off-targets"},
    {"name": "v42_y398_T",      "seq": "EKTLLLSIGK", "family": "F2_Y398", "p1": "E",
     "handles": "392+398",
     "rationale": "P3=T, hydroxyl H-bonds Y398; methyl adds van der Waals"},
    {"name": "v42_y398_T_V",    "seq": "VKTLLLSIGK", "family": "F2_Y398", "p1": "V",
     "handles": "392+398",
     "rationale": "P1=V (K392) + P3=T (Y398), two-handle design"},
    {"name": "v42_y398_T_D",    "seq": "DKTLLLSIGK", "family": "F2_Y398", "p1": "D",
     "handles": "392+398",
     "rationale": "P1=D (contrast) + P3=T (Y398)"},

    # --- Family 3: Anti-IRAP K406 (P5 substitution) ---
    # ERAP2 A406 is small; IRAP K406 is bulky+charged → acidic P5 clashes with K406
    {"name": "v42_k406_E",      "seq": "EKLLELSIGK", "family": "F3_antiIRAP", "p1": "E",
     "handles": "392+406",
     "rationale": "P5=E, charge-charge repulsion with IRAP K406, neutral to ERAP2 A406"},
    {"name": "v42_k406_D",      "seq": "EKLLDLSIGK", "family": "F3_antiIRAP", "p1": "E",
     "handles": "392+406",
     "rationale": "P5=D, charge-charge repulsion with IRAP K406"},
    {"name": "v42_k406_E_V",    "seq": "VKLLELSIGK", "family": "F3_antiIRAP", "p1": "V",
     "handles": "392+406",
     "rationale": "P1=V (K392) + P5=E (anti-IRAP K406)"},

    # --- Family 4: Multi-handle combinations ---
    # Stack 2-3 ERAP2-specific handles for maximum selectivity
    {"name": "v42_multi_TE",    "seq": "EKTLELSIGK", "family": "F4_multi", "p1": "E",
     "handles": "392+398+406",
     "rationale": "P3=T (Y398 H-bond) + P5=E (anti-IRAP K406), two ERAP2-specific handles"},
    {"name": "v42_multi_TVL",   "seq": "EKTVLLSIGK", "family": "F4_multi", "p1": "E",
     "handles": "392+398+403",
     "rationale": "P3=T (Y398) + P4=V (fits A403, penalized by S403 in off-targets)"},
    {"name": "v42_multi_TVE",   "seq": "EKTVELSIGK", "family": "F4_multi", "p1": "E",
     "handles": "392+398+403+406",
     "rationale": "Triple handle: P3=T (Y398) + P4=V (A403) + P5=E (anti-K406)"},
    {"name": "v42_multi_TVE_V", "seq": "VKTVELSIGK", "family": "F4_multi", "p1": "V",
     "handles": "392+398+403+406",
     "rationale": "Quad handle: P1=V (K392) + P3=T + P4=V + P5=E"},
    {"name": "v42_multi_TVE_D", "seq": "DKTVELSIGK", "family": "F4_multi", "p1": "D",
     "handles": "392+398+403+406",
     "rationale": "Triple+contrast: P1=D + P3=T + P4=V + P5=E"},
    {"name": "v42_multi_TVEK",  "seq": "EKTVELKIGK", "family": "F4_multi", "p1": "E",
     "handles": "392+398+403+406+414",
     "rationale": "Five handles: P3=T + P4=V + P5=E + P7=K (salt-bridge D414)"},
    {"name": "v42_multi_TVER",  "seq": "EKTVELRIGK", "family": "F4_multi", "p1": "E",
     "handles": "392+398+403+406+414",
     "rationale": "Five handles: P3=T + P4=V + P5=E + P7=R (salt-bridge D414, longer reach)"},
    {"name": "v42_max_VTVEK",   "seq": "VKTVELKIGK", "family": "F4_multi", "p1": "V",
     "handles": "392+398+403+406+414",
     "rationale": "Maximum stack: P1=V + P3=T + P4=V + P5=E + P7=K, 5 ERAP2-specific handles"},
]


def get_all_peptides():
    return PEPTIDES


def write_boltz2_yamls(peptides, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for pep in peptides:
        for target_name, target_seq in TARGETS.items():
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


def print_design_table(peptides):
    print(f"{'#':>2}  {'Name':>22}  {'Sequence':>12}  {'Len':>3}  {'P1':>2}  {'Family':>14}  {'Handles':>20}")
    print("-" * 95)
    for i, p in enumerate(peptides, 1):
        print(f"{i:>2}  {p['name']:>22}  {p['seq']:>12}  {len(p['seq']):>3}  {p['p1']:>2}  {p['family']:>14}  {p['handles']:>20}")


# ── Analysis ──────────────────────────────────────────────────────────────────

def collect_scores(results_dir):
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


def analyze_results(results_dir):
    results_dir = Path(results_dir)
    scores = collect_scores(results_dir)

    if not scores:
        print("No confidence scores found.")
        return

    pep_lookup = {p["name"]: p for p in PEPTIDES}

    # Build per-peptide, per-target averages
    data = {}  # pep_name -> {target: mean_iptm, ...}
    for run_name, iptm_list in scores.items():
        # Parse: v42_base_E_ERAP2_K392
        for target_name in TARGETS:
            suffix = f"_{target_name}"
            if run_name.endswith(suffix):
                pep_name = run_name[:-len(suffix)]
                break
        else:
            continue

        if pep_name not in data:
            data[pep_name] = {}
        data[pep_name][target_name] = {
            "mean": statistics.mean(iptm_list),
            "std": statistics.stdev(iptm_list) if len(iptm_list) > 1 else 0,
            "n": len(iptm_list),
        }

    # ── Results table ─────────────────────────────────────────────────────
    print("=" * 140)
    print("V4.2 MULTI-HANDLE CAMPAIGN RESULTS")
    print("=" * 140)
    print()

    # Get baseline (v42_base_E) ERAP1/IRAP scores for comparison
    baseline_erap1 = data.get("v42_base_E", {}).get("ERAP1", {}).get("mean", 0)
    baseline_irap = data.get("v42_base_E", {}).get("IRAP", {}).get("mean", 0)

    rows = []
    for pep_name in sorted(data.keys(), key=lambda k: (pep_lookup.get(k, {}).get("family", "z"), k)):
        d = data[pep_name]
        meta = pep_lookup.get(pep_name, {})

        k392 = d.get("ERAP2_K392", {}).get("mean", 0)
        n392 = d.get("ERAP2_N392", {}).get("mean", 0)
        erap1 = d.get("ERAP1", {}).get("mean", 0)
        irap = d.get("IRAP", {}).get("mean", 0)
        delta_kn = k392 - n392
        erap2_best = max(k392, n392)

        # Selectivity: ERAP2 binding vs off-target binding
        selectivity_erap1 = erap2_best - erap1
        selectivity_irap = erap2_best - irap

        rows.append({
            "name": pep_name,
            "seq": meta.get("seq", "?"),
            "family": meta.get("family", "?"),
            "p1": meta.get("p1", "?"),
            "handles": meta.get("handles", "?"),
            "k392": k392,
            "n392": n392,
            "erap1": erap1,
            "irap": irap,
            "delta_kn": delta_kn,
            "erap2_best": erap2_best,
            "sel_erap1": selectivity_erap1,
            "sel_irap": selectivity_irap,
            "rationale": meta.get("rationale", ""),
        })

    # Print by family
    for fam in ["F1_P1scan", "F2_Y398", "F3_antiIRAP", "F4_multi"]:
        fam_rows = [r for r in rows if r["family"] == fam]
        if not fam_rows:
            continue

        labels = {
            "F1_P1scan": "FAMILY 1: P1 Scan (baseline)",
            "F2_Y398": "FAMILY 2: Y398 Handle (P3 substitution)",
            "F3_antiIRAP": "FAMILY 3: Anti-IRAP K406 (P5 substitution)",
            "F4_multi": "FAMILY 4: Multi-Handle Stack",
        }
        print(f"### {labels[fam]}")
        print(f"{'Name':>22}  {'Seq':>12}  {'P1':>2}  {'K392':>6}  {'N392':>6}  {'dKN':>6}  {'ERAP1':>6}  {'IRAP':>6}  {'selE1':>6}  {'selIR':>6}  {'Handles':>20}")
        print("-" * 130)

        for r in fam_rows:
            print(f"{r['name']:>22}  {r['seq']:>12}  {r['p1']:>2}  "
                  f"{r['k392']:>.3f}  {r['n392']:>.3f}  {r['delta_kn']:>+.3f}  "
                  f"{r['erap1']:>.3f}  {r['irap']:>.3f}  "
                  f"{r['sel_erap1']:>+.3f}  {r['sel_irap']:>+.3f}  "
                  f"{r['handles']:>20}")
        print()

    # ── Ranking ───────────────────────────────────────────────────────────
    print("=" * 100)
    print("COMPOSITE RANKING (ERAP2 channel engagement + selectivity over off-targets)")
    print("=" * 100)
    print()

    # Composite score: ERAP2_best + selectivity_erap1 + selectivity_irap
    # Higher = better ERAP2 binding AND lower off-target binding
    for r in rows:
        r["composite"] = r["erap2_best"] + r["sel_erap1"] + r["sel_irap"]

    ranked = sorted(rows, key=lambda r: r["composite"], reverse=True)

    print(f"{'Rank':>4}  {'Name':>22}  {'Seq':>12}  {'ERAP2':>6}  {'ERAP1':>6}  {'IRAP':>6}  {'selE1':>6}  {'selIR':>6}  {'Score':>6}  {'Handles':>20}")
    print("-" * 120)
    for i, r in enumerate(ranked, 1):
        marker = " ***" if i <= 4 else ""
        print(f"{i:>4}  {r['name']:>22}  {r['seq']:>12}  "
              f"{r['erap2_best']:>.3f}  {r['erap1']:>.3f}  {r['irap']:>.3f}  "
              f"{r['sel_erap1']:>+.3f}  {r['sel_irap']:>+.3f}  "
              f"{r['composite']:>.3f}  {r['handles']:>20}{marker}")
    print()

    # ── Top 4 shortlist ───────────────────────────────────────────────────
    print("=" * 100)
    print("TOP 4 SHORTLIST FOR SYNTHESIS")
    print("=" * 100)
    for i, r in enumerate(ranked[:4], 1):
        print(f"\n  #{i}: {r['name']}")
        print(f"      Sequence: {r['seq']} ({len(r['seq'])}-mer, P1={r['p1']})")
        print(f"      ERAP2 K392: {r['k392']:.3f}  |  N392: {r['n392']:.3f}  |  delta: {r['delta_kn']:+.3f}")
        print(f"      ERAP1: {r['erap1']:.3f}  |  IRAP: {r['irap']:.3f}")
        print(f"      Selectivity: vs ERAP1 {r['sel_erap1']:+.3f}  |  vs IRAP {r['sel_irap']:+.3f}")
        print(f"      Handles: {r['handles']}")
        print(f"      Rationale: {r['rationale']}")

    # ── Save ──────────────────────────────────────────────────────────────
    output = {
        "experiment": "V4.2 multi-handle ERAP2 channel selectivity campaign",
        "n_peptides": len(PEPTIDES),
        "targets": list(TARGETS.keys()),
        "diffusion_samples": DIFFUSION_SAMPLES,
        "ranked_results": ranked,
        "top4": ranked[:4],
        "baseline_erap1": baseline_erap1,
        "baseline_irap": baseline_irap,
    }
    output_path = results_dir / "v42_analysis.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved analysis to {output_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="V4.2 Multi-Handle Campaign")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--results-dir", type=str, default=str(RESULTS_DIR))
    args = parser.parse_args()

    if args.analyze:
        analyze_results(Path(args.results_dir))
        return

    peptides = get_all_peptides()
    print("=" * 95)
    print(f"V4.2 MULTI-HANDLE CAMPAIGN: {len(peptides)} peptides x {len(TARGETS)} targets = {len(peptides)*len(TARGETS)} runs")
    print("=" * 95)
    print()

    print_design_table(peptides)
    print()

    n = write_boltz2_yamls(peptides, YAML_DIR)
    print(f"Wrote {n} Boltz-2 input YAMLs to {YAML_DIR}/")

    defn_path = SCRIPT_DIR / "v42_experiment_definition.json"
    with open(defn_path, "w") as f:
        json.dump(peptides, f, indent=2)
    print(f"Saved experiment definition to {defn_path}")


if __name__ == "__main__":
    main()
