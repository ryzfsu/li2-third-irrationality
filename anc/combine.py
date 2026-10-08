#!/usr/bin/env python3
"""Combine the per-item outputs into results.json with a PASS/FAIL per item."""
import json, sys
e = json.load(open('energy_results.json'))
a = json.load(open('arith_results.json'))
c = json.load(open('caps_results.json'))
p = json.load(open('parsed_tex.json'))
x = json.load(open('crosscheck_indep_caps_results.json')) if len(sys.argv) > 1 else None
# Compare each computed cap endpoint with the printed ``tab:caps`` bound,
# rounding the computed upper end upward at 1e-6.
from fractions import Fraction as Fr
import math
tab = {'x': 0, 'v': 1, 'z': 2, 'circ': 3}
cmp = {}
for s, i in tab.items():
    row = p['caps_table'][i]
    w = Fr(c['supports'][s]['worst_upper_end'])
    up6 = Fr(math.ceil(w * 10 ** 6), 10 ** 6)
    cmp[s] = {'computed_worst_upper_end': c['supports'][s]['worst_upper_end'],
              'rounded_up_1e-6': str(float(up6)), 'printed_upper_bound': row['upper_bound'],
              'printed_matches': up6 == Fr(row['upper_bound']),
              'leaves': c['supports'][s]['leaves'], 'printed_leaves': row['leaves'],
              'printed_cap': row['cap'], 'cap_minus_worst': c['supports'][s]['cap_minus_worst']}
least = min(cmp, key=lambda s: Fr(c['supports'][s]['cap']) - Fr(c['supports'][s]['worst_upper_end']))
out = {
 'source': {'tex': 'input_main.tex', 'tex_sha256': p['tex_sha256']},
 'item1_reference': {'status': 'PASS' if e['item1_reference']['pass'] else 'FAIL', 'detail': e['item1_reference']},
 'item2_energy': {'status': 'PASS' if e['item2_energy']['pass'] else 'FAIL', 'detail': e['item2_energy']},
 'item3_caps': {'status': 'PASS' if c['pass'] else 'FAIL', 'table_B2_comparison': cmp,
                'least_margin_support': least,
                'least_margin_exceeds_2e-6': Fr(c['supports'][least]['cap']) - Fr(c['supports'][least]['worst_upper_end']) > Fr(2, 10 ** 6),
                'total_leaves': c['total_leaves'], 'detail': c},
 'item4_arithmetic': {'status': 'PASS' if a['item4_arithmetic']['pass'] else 'FAIL', 'detail': a['item4_arithmetic']},
 'item5_alpha0': {'status': 'PASS' if a['item5_alpha0']['pass'] else 'FAIL', 'detail': a['item5_alpha0']},
}
if x is not None:
    out['crosscheck_reference_impl'] = {'status': x.get('status'), 'supports': {k: {kk: v[kk] for kk in
        ('base_leaves', 'accepted_leaves', 'bisected_failed_nodes', 'worst_upper_bound', 'cap_minus_worst')}
        for k, v in x['supports'].items()}, 'alpha': x['alpha']}
json.dump(out, open('results.json', 'w'), indent=1)
print(json.dumps({k: (v.get('status') if isinstance(v, dict) else v) for k, v in out.items()}, indent=1))
print(json.dumps(cmp, indent=1), least, out['item3_caps']['least_margin_exceeds_2e-6'])
