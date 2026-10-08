"""Rebuild the reference measure from the parsed TeX tables (item 1).

Exact rational data. The lower-half Poisson centres are produced by
conjugation of the printed rows k>=1; row k=0 (real centre) is single.
"""
import json
from fractions import Fraction as Fr

R = Fr(2, 5)


def load(path):
    d = json.load(open(path))
    chains = {}
    for key, ch in d['chains'].items():
        b0 = Fr(ch['b0'])
        T = [b0]
        W = []
        ks = []
        for k, bk, wk in ch['rows']:
            ks.append(k)
            T.append(Fr(bk))
            W.append(Fr(wk))
        chains[key] = {'T': T, 'W': W, 'ks': ks}
    rows = [(k, Fr(a), Fr(b), Fr(w)) for k, a, b, w in d['poisson_rows']]
    poi = []
    for k, a, b, w in rows:
        if k == 0:
            poi.append((a, b, w))
        else:
            poi.append((a, b, w))
            poi.append((a, -b, w))
    return d, chains, rows, poi


def checks(d, chains, rows, poi):
    res = {}
    ok = True
    exp_n = {'x': 131, 'v': 143, 'z_int': 105}
    exp_end = {'x': Fr(1), 'v': Fr(1), 'z_int': Fr(3)}
    for key, ch in chains.items():
        T, W, ks = ch['T'], ch['W'], ch['ks']
        c = {
            'rows': len(W),
            'rows_expected': exp_n[key],
            'indices_consecutive_1..n': ks == list(range(1, len(W) + 1)),
            'strictly_increasing_contiguous': all(T[i + 1] > T[i] for i in range(len(W))),
            'all_weights_positive': all(w > 0 for w in W),
            'start': str(T[0]), 'end': str(T[-1]), 'end_expected': str(exp_end[key]),
            'mass': str(sum(W)),
            'max_density_H': float(max(W[i] / (T[i + 1] - T[i]) for i in range(len(W)))),
        }
        c['pass'] = (c['rows'] == exp_n[key] and c['indices_consecutive_1..n'] and
                     c['strictly_increasing_contiguous'] and c['all_weights_positive'] and
                     T[-1] == exp_end[key])
        res[key] = c
        ok &= c['pass']
    mx = sum(chains['x']['W']); mv = sum(chains['v']['W']); mzi = sum(chains['z_int']['W'])
    mp = sum(w for _, _, w in poi)
    res['mass_x'] = str(mx); res['mass_v'] = str(mv); res['mass_z_int'] = str(mzi)
    res['mass_z_poisson'] = str(mp); res['mass_z_total'] = str(mzi + mp)
    res['mass_z_int_matches_text'] = mzi == Fr(170340748021, 500000000000)
    res['mass_z_poisson_matches_text'] = mp == Fr(329659251979, 500000000000)
    res['masses_(1,5/4,1)'] = (mx == 1 and mv == Fr(5, 4) and mzi + mp == 1)
    ok &= res['masses_(1,5/4,1)'] and res['mass_z_int_matches_text'] and res['mass_z_poisson_matches_text']
    res['v_chain_start_equals_vzero_macro'] = chains['v']['T'][0] == Fr(d['macros']['vzero'])
    res['x_chain_start'] = float(chains['x']['T'][0])
    res['z_chain_start'] = float(chains['z_int']['T'][0])
    res['z_chain_start_in_(1,3)'] = 1 < chains['z_int']['T'][0] < 3
    # chains inside the declared supports
    res['supp_x_in_[0,1]'] = chains['x']['T'][0] >= 0 and chains['x']['T'][-1] <= 1
    res['supp_v_in_[v0,1]_v0>0'] = chains['v']['T'][0] > 0 and chains['v']['T'][-1] <= 1
    res['supp_zint_in_[1,3]'] = chains['z_int']['T'][0] >= 1 and chains['z_int']['T'][-1] <= 3
    ok &= res['supp_x_in_[0,1]'] and res['supp_v_in_[v0,1]_v0>0'] and res['supp_zint_in_[1,3]'] and res['v_chain_start_equals_vzero_macro']
    # r not an endpoint of the v chain
    res['r_not_v_endpoint'] = R not in chains['v']['T']
    ok &= res['r_not_v_endpoint']
    # Poisson
    ks = [k for k, _, _, _ in rows]
    res['poisson_rows'] = len(rows)
    res['poisson_measures'] = len(poi)
    res['poisson_row_indices_0..83'] = ks == list(range(84))
    res['poisson_row0_real'] = rows[0][2] == 0
    res['poisson_rows_k>=1_Im>0'] = all(b > 0 for k, a, b, w in rows if k >= 1)
    res['poisson_weights_positive'] = all(w > 0 for _, _, _, w in rows)
    res['poisson_|c|<r_all'] = all(a * a + b * b < R * R for _, a, b, _ in rows)
    res['poisson_max_|c|^2/r^2'] = float(max((a * a + b * b) / (R * R) for _, a, b, _ in rows))
    # increasing order of arg c for k>=1 (Im>0, so arg increasing <=> a/|c| decreasing;
    # compare via cross product sign: arg(c_{k+1})>arg(c_k) <=> Im(conj(c_k) c_{k+1}) > 0)
    inc = True
    for (k1, a1, b1, _), (k2, a2, b2, _) in zip(rows[1:], rows[2:]):
        if not (a1 * b2 - b1 * a2 > 0):
            inc = False
    res['poisson_rows_increasing_arg'] = inc
    S = {}
    for a, b, w in poi:
        S[(a, b, w)] = S.get((a, b, w), 0) + 1
    res['poisson_conjugation_invariant'] = all(S.get((a, -b, w), 0) == n for (a, b, w), n in S.items())
    ok &= (len(poi) == 167 and res['poisson_row_indices_0..83'] and res['poisson_row0_real'] and
           res['poisson_rows_k>=1_Im>0'] and res['poisson_weights_positive'] and
           res['poisson_|c|<r_all'] and res['poisson_conjugation_invariant'])
    res['pass'] = bool(ok)
    return res
