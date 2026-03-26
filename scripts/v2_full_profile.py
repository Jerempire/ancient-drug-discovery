"""Pull together all V2 variant data and identify best candidates for dual K/N binding."""
import sys, json, os
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load variant results
with open(os.path.join(PROJECT, "data/results/n248_variants/variant_results.json")) as f:
    variants = json.load(f)

print("=" * 90)
print("V2 BINDER FULL VARIANT PROFILE")
print("=" * 90)

# Print all variants with available scores
print(f"\n{'Name':>35} {'E2 ipTM':>9} {'E1 ipTM':>9} {'Delta':>8} {'Len':>5}")
print("-" * 70)
for v in variants:
    e2 = v.get('erap2_iptm', 0)
    e1 = v.get('erap1_iptm')
    delta = v.get('iptm_delta')
    name = v['name']
    length = v['length']
    e1_s = f"{e1:.3f}" if e1 is not None else "?"
    d_s = f"{delta:+.4f}" if delta is not None else "?"
    print(f"{name:>35} {e2:>8.3f} {e1_s:>9} {d_s:>8} {length:>5}")

# Find the best: highest ERAP2, lowest ERAP1
print("\n" + "=" * 90)
print("RANKING BY SELECTIVITY (ERAP2 - ERAP1 delta)")
print("=" * 90)
with_delta = [v for v in variants if v.get('iptm_delta') is not None]
with_delta.sort(key=lambda v: v['iptm_delta'], reverse=True)
print(f"\n{'Name':>35} {'E2':>8} {'E1':>8} {'Delta':>8}")
print("-" * 65)
for v in with_delta:
    print(f"{v['name']:>35} {v['erap2_iptm']:>7.3f} {v['erap1_iptm']:>7.3f} {v['iptm_delta']:>+7.4f}")

# Check what ERAP2 sequence was used (K392 or N392?)
print("\n" + "=" * 90)
print("KEY QUESTION: Was the original screen done with K392 or N392 ERAP2?")
print("=" * 90)

# Check if there are separate K/N results
k_results = [v for v in variants if 'k392' in v.get('name', '').lower()]
n_results = [v for v in variants if 'n392' in v.get('name', '').lower()]
print(f"  K392-specific variants found: {len(k_results)}")
print(f"  N392-specific variants found: {len(n_results)}")

# Look for the ERAP2 sequence used
yaml_dir = os.path.join(PROJECT, "data/boltz2_inputs")
if os.path.exists(yaml_dir):
    for yf in os.listdir(yaml_dir):
        if 'erap2' in yf.lower() and yf.endswith('.yaml'):
            print(f"  Found YAML: {yf}")

# The cropped ERAP2 sequence from our CIF analysis
print(f"\n  ERAP2 channel sequence used (from CIF, residues 350-500):")
print(f"  LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDS...")
print(f"  Position 392 (cropped pos 43) = K (Lysine)")
print(f"  This IS the K392 allele. N392 would have N at that position.")

# Identify the best V2 variant
print("\n" + "=" * 90)
print("BEST V2 VARIANT FOR DUAL K/N OPTIMIZATION")
print("=" * 90)

best = max(with_delta, key=lambda v: v['iptm_delta'])
print(f"\n  Best selectivity: {best['name']}")
print(f"    ERAP2: {best['erap2_iptm']:.3f}")
print(f"    ERAP1: {best['erap1_iptm']:.3f}")
print(f"    Delta: {best['iptm_delta']:+.4f}")
print(f"    Length: {best['length']}aa")
print(f"    Sequence: {best['sequence']}")

# What do we need to test?
print("\n" + "=" * 90)
print("WHAT WE NEED TO RUN ON BOLTZ-2")
print("=" * 90)
print("""
  For each of the top 3 V2 variants, test against 4 targets:

  1. ERAP2-K392 (current screen target)
  2. ERAP2-N392 (K->N mutation at position 392)
  3. ERAP1 (selectivity counterscreen)
  4. IRAP (selectivity counterscreen)

  The N392 ERAP2 sequence is identical to K392 except position 392:
  K392: ...KEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDS...
                K <-- position 43 in cropped (392 in full)
  N392: ...KEGFANYMELIAVNATYPELQFDDYFLNVCFEVITKDS...
                N <-- Asparagine instead of Lysine

  Top 3 candidates to test:
""")

# Rank by ERAP2 binding * delta (balance of binding + selectivity)
scored = [(v, v['erap2_iptm'] * v['iptm_delta']) for v in with_delta if v['erap2_iptm'] > 0.5]
scored.sort(key=lambda x: x[1], reverse=True)

for i, (v, score) in enumerate(scored[:3]):
    print(f"  #{i+1}: {v['name']} (E2={v['erap2_iptm']:.3f}, E1={v['erap1_iptm']:.3f}, "
          f"delta={v['iptm_delta']:+.3f}, len={v['length']})")

# Build the N392 ERAP2 sequence
erap2_k392 = "LFDPKTSSASDKLWVTRVIAHELAHQWFGNLVTMEWWNDIWLKEGFAKYMELIAVNATYPELQFDDYFLNVCFEVITKDSLNSSRPISKPAETPTQIQEMFDEVSYNKGACILNMLKDFLGEEKFQKGIIQYLKKFSYRNAKNDDLWSSLS"
# Position 392 = index 42 (0-based) in the cropped sequence (350+42=392)
pos_392 = 392 - 350  # = 42
erap2_n392 = erap2_k392[:pos_392] + 'N' + erap2_k392[pos_392+1:]

print(f"\n  K392 sequence (pos {pos_392}): ...{erap2_k392[pos_392-3:pos_392+4]}...")
print(f"  N392 sequence (pos {pos_392}): ...{erap2_n392[pos_392-3:pos_392+4]}...")
print(f"  Differ at position {pos_392}: K->N")

# Verify
assert erap2_k392[pos_392] == 'K', f"Expected K at pos {pos_392}, got {erap2_k392[pos_392]}"
assert erap2_n392[pos_392] == 'N'
assert len(erap2_k392) == len(erap2_n392)
diff_count = sum(1 for a, b in zip(erap2_k392, erap2_n392) if a != b)
print(f"  Total differences: {diff_count} (should be 1)")

# Save sequences for YAML generation
output = {
    "erap2_k392": erap2_k392,
    "erap2_n392": erap2_n392,
    "erap1": "TVAHELAHQWFGNLVTMEWWNDLWLNEGFAKFMEFVSVSVTHPELKVGDYFFGKCFDAMEVDALNSSHPVSTPVENPAQIREMFDDVSYDKGACILNMLREYLSADAFKSGIVQYLQKHSYKNTKNEDLWDSMASICPTDGVKGMDGFCSR",
    "irap": "LYDSNTSSMADRKLVTKIIAHELAHQWFGNLVTMKWWNDLWLNEGFATFMEYFSLEKIFKELSSYEDFLDARFKTMKKDSLNSSHPISSSVQSSEQIEEMFDSLSYFKGSSLLLMLKTYLSEDVFQHAVVLYLHNHSYASIQSDDLWDSFN",
    "top_candidates": [
        {"name": v['name'], "sequence": v['sequence'], "length": v['length'],
         "erap2_iptm": v['erap2_iptm'], "erap1_iptm": v.get('erap1_iptm'),
         "delta": v.get('iptm_delta')}
        for v, _ in scored[:3]
    ],
}
out_path = os.path.join(PROJECT, "data/results/v43_validation/terminal_e/v2_dual_binding_setup.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: {out_path}")
