#!/usr/bin/env python3
"""Items 4 and 5: exact rational recomputation from the printed slab table and the
definitions (3.1)-(3.2), plus alpha_0. Pure Fraction arithmetic; data from parsed_tex.json only.
"""
import sys, json, math, time
from fractions import Fraction as Fr

d = json.load(open(sys.argv[1]))
OUT = sys.argv[2]
M = {k: v for k, v in d['macros'].items()}
R = int(M['Rseam'])
assert R == 8
lo_r, hi_r = Fr(1, R), Fr(5, 4)
T0 = time.time()
res = {}

# ---------------- definitions (3.1)-(3.2) ----------------
H, Q4, Q2, T2, T4 = Fr(1, 2), Fr(1, 4), Fr(1, 2), Fr(3, 2), Fr(5, 2)


def cells(x, rho):
    K = [k for k in range(0, 64) if Fr(1, 2) <= x + k * rho < Fr(3, 2)]
    L = [l for l in range(0, 64) if Fr(1, 4) <= x + l * rho < Fr(5, 2)]
    W = [w for w in range(1, 64) if 0 < w * rho - x <= Fr(1, 4)]
    return K, L, W


def rr(x, rho, with_ind=True):
    K, L, W = cells(x, rho)
    v = min(2 * len(K) + len(W), len(L))
    special = (K == [0] and L == [0, 1] and W == [])
    if with_ind and special:
        v -= 1
    return v, special


# line families, from the index ranges 4z0/5-1 < k < z0 R and 1 <= w < R/4+1
def krange(z0):
    return [k for k in range(0, 200) if Fr(4, 5) * z0 - 1 < k < z0 * R]


fam = {
    'L_lower(1/4)': krange(Fr(1, 4)), 'K_lower(1/2)': krange(Fr(1, 2)),
    'K_upper(3/2)': krange(Fr(3, 2)), 'L_upper(5/2)': krange(Fr(5, 2)),
    'W': [w for w in range(1, 100) if 1 <= w < Fr(R, 4) + 1],
}
printed_fam = {'L_lower(1/4)': list(range(0, (R - 1) // 4 + 1)), 'K_lower(1/2)': list(range(0, (R - 1) // 2 + 1)),
               'K_upper(3/2)': list(range(1, (3 * R - 1) // 2 + 1)), 'L_upper(5/2)': list(range(2, (5 * R - 1) // 2 + 1)),
               'W': list(range(1, (R + 3) // 4 + 1))}
# printed R=8 line-family table (label ``tab:line-families``):
# 0..1, 0..3, 1..11, 2..19, 1..2; total 39.
printed_rendered = {'L_lower(1/4)': (0, 1), 'K_lower(1/2)': (0, 3), 'K_upper(3/2)': (1, 11),
                    'L_upper(5/2)': (2, 19), 'W': (1, 2)}
fam_ok = all(fam[k] == printed_fam[k] == list(range(printed_rendered[k][0], printed_rendered[k][1] + 1)) for k in fam)
nlines = 2 + sum(len(v) for v in fam.values())
res['line_families'] = {k: [v[0], v[-1], len(v)] for k, v in fam.items()}
res['line_families_match_printed'] = fam_ok
res['nlines'] = nlines
res['nlines_is_39'] = nlines == 39

# lines x = alpha + beta*rho
lines = [(Fr(0), Fr(0)), (Fr(0), Fr(1))]
for z0, key in ((Fr(1, 4), 'L_lower(1/4)'), (Fr(1, 2), 'K_lower(1/2)'), (Fr(3, 2), 'K_upper(3/2)'), (Fr(5, 2), 'L_upper(5/2)')):
    for k in fam[key]:
        lines.append((z0, Fr(-k)))
for w in fam['W']:
    lines.append((Fr(-1, 4), Fr(w)))

# the cut set C
C = set()
for i in range(len(lines)):
    for j in range(i + 1, len(lines)):
        a1, b1 = lines[i]; a2, b2 = lines[j]
        if b1 == b2:
            continue
        rho = (a2 - a1) / (b1 - b2)
        x = a1 + b1 * rho
        if lo_r < rho < hi_r and 0 <= x <= rho:
            C.add(rho)
C = sorted(C)
res['cut_points_in_(1/R,5/4)'] = len(C)

slabs = [[Fr(s) for s in row] for row in d['slabs']]
res['slab_rows'] = len(slabs)
consec = slabs[0][0] == lo_r and slabs[-1][1] == hi_r and all(slabs[i][1] == slabs[i + 1][0] for i in range(len(slabs) - 1))
res['slabs_consecutive_from_1/R_to_5/4'] = consec
internal = [s[1] for s in slabs[:-1]]
res['slab_internal_endpoints_equal_cut_set'] = internal == C
if internal != C:
    res['cut_minus_slab'] = [str(x) for x in sorted(set(C) - set(internal))]
    res['slab_minus_cut'] = [str(x) for x in sorted(set(internal) - set(C))]


def c_of(rho, with_ind=True):
    """c(rho) = 1 + int_0^rho r(x,rho) dx, exact; also the measure of the special cell and #lines strictly inside."""
    xs = {Fr(0), rho}
    inside = 0
    for a, b in lines[2:]:
        x = a + b * rho
        if 0 < x < rho:
            xs.add(x); inside += 1
    xs = sorted(xs)
    tot = Fr(0); sp = Fr(0)
    for u0, u1 in zip(xs, xs[1:]):
        v, s = rr((u0 + u1) / 2, rho, with_ind)
        tot += (u1 - u0) * v
        if s:
            sp += u1 - u0
    return 1 + tot, sp, inside


bad = []
five = True
integral = Fr(0)
special_int = Fr(0)
for i, (r0, r1, c0, c1) in enumerate(slabs):
    pts = [r0 + (r1 - r0) * Fr(k, 4) for k in (1, 2, 3)]
    vals = [c_of(p) for p in pts]
    for p, (cv, sp, ins) in zip(pts, vals):
        if cv != c0 + c1 * p:
            bad.append((i, str(p), str(cv), str(c0 + c1 * p)))
        if ins != 5:
            five = False
    # independent slope/intercept recovered from two interior points
    (ca, _, _), (cb, _, _) = vals[0], vals[2]
    slope = (cb - ca) / (pts[2] - pts[0]); icpt = ca - slope * pts[0]
    if slope != c1 or icpt != c0:
        bad.append((i, 'affine-fit', str(icpt), str(slope)))
    integral += (r1 - r0) * (c0 + c1 * (r0 + r1) / 2)
    special_int += (r1 - r0) * vals[1][1]       # special-cell measure is affine in rho on a slab
res['slab_rows_mismatching_definitions'] = bad
res['every_slab_has_exactly_5_interior_lines'] = five
inner_printed = Fr(M['innervalue'])
res['inner_integral'] = str(integral)
res['inner_integral_matches_printed'] = integral == inner_printed
res['rank_one_block_contribution'] = str(-special_int)
res['rank_one_block_is_-5/48'] = -special_int == Fr(-5, 48)
# direct closed form of the -5/48 computation
cf = -(Fr(3, 2) * (1 - Fr(11, 12) ** 2) - Fr(11, 4) * (1 - Fr(11, 12))) - (Fr(1, 2) * (Fr(25, 16) - 1) - Fr(3, 4) * Fr(1, 4))
res['rank_one_closed_form'] = str(cf)
outer = (Fr(19, 4) * Fr(3, 4) - (4 - Fr(25, 16))) + (Fr(11, 4) * Fr(3, 4) - (Fr(121, 16) - 4) / 2)
res['outer_integral'] = str(outer)
res['outer_is_45/32'] = outer == Fr(45, 32)
A = Fr(13, 4 * R) + integral + Fr(45, 32)
res['A'] = str(A)
res['A_matches_printed'] = A == Fr(M['Avalue'])

# lattice-sum identity (3.6) for small m (illustrative finite test)
lat_ok = True; lat_n = 0
for n in range(1, 41):
    m = 4 * n; U = 11 * n - 1
    for p in range(2, 5 * m // 4 + 1):
        if all(p % t for t in range(2, int(p ** .5) + 1)) and Fr(m, R) < p < Fr(5 * m, 4) and p * p > U:
            s = sum(rr(Fr(u, m), Fr(p, m))[0] for u in range(p))
            cv = c_of(Fr(p, m))[0]
            if s != m * (cv - 1):
                lat_ok = False
            lat_n += 1
res['lattice_sum_identity_checked_pairs'] = lat_n
res['lattice_sum_identity_ok'] = lat_ok

# ---------------- tau with the paper's log 3 bounds ----------------
L3lo = 2 * sum(Fr(1, (2 * k + 1) * 2 ** (2 * k + 1)) for k in range(100))
L3hi = L3lo + Fr(8, 3 * 201 * 2 ** 201)
assert L3lo < L3hi
# sanity vs float
res['log3_bounds_float'] = [float(L3lo), float(L3hi), math.log(3)]
tau_lo = A - Fr(53, 48) * L3hi
tau_hi = A - Fr(53, 48) * L3lo


def dec(x, n=22):
    s = '-' if x < 0 else ''
    x = abs(x)
    ip = x.numerator // x.denominator
    fr = x - ip
    return s + str(ip) + '.' + str((fr * 10 ** n).numerator // (fr * 10 ** n).denominator).zfill(n)


res['tau_lower'] = dec(tau_lo); res['tau_upper'] = dec(tau_hi)
tp = Fr(M['tauvalue_decimal'])
res['tau_printed'] = M['tauvalue_decimal']
res['tau_printed_digits_consistent'] = dec(tau_lo, 18) == M['tauvalue_decimal'] == dec(tau_hi, 18)
res['tau_upper_lt_taubound_3.87273'] = tau_hi < Fr(M['taubound'])

# ---------------- item 5: alpha_0 ----------------
qm, Kx, Kv, Kz, Kzc = (Fr(M[k]) for k in ('qminus', 'Kx', 'Kv', 'Kz', 'Kzcirc'))
C0 = Fr(1647919, 10 ** 6)
a0 = qm - C0 - Kx - Fr(5, 4) * Kv - Kz
it5 = {'q_minus': str(qm), 'K_x': str(Kx), 'K_v': str(Kv), 'K_z': str(Kz), 'K_zcirc': str(Kzc),
       'C0': str(C0), 'alpha0_exact': str(a0), 'alpha0_decimal': dec(a0, 8),
       'alpha0_matches_printed_3.875784': a0 == Fr(M['alphazero']),
       'three_halves_log3_upper_lt_C0': Fr(3, 2) * L3hi < C0,
       'three_halves_log3_bounds': [dec(Fr(3, 2) * L3lo, 12), dec(Fr(3, 2) * L3hi, 12)],
       'Kzcirc_lt_Kz': Kzc < Kz,
       'alpha0_gt_alphabound_3.8757': a0 > Fr(M['alphabound']),
       'alpha0_minus_tau_lower': dec(a0 - tau_hi, 12), 'alpha0_minus_tau_upper': dec(a0 - tau_lo, 12),
       'alpha0_minus_tau_gt_0.00305': a0 - tau_hi > Fr(M['marginvalue']),
       'taubound_lt_alphabound': Fr(M['taubound']) < Fr(M['alphabound'])}
it5['pass'] = all(it5[k] for k in ('alpha0_matches_printed_3.875784', 'three_halves_log3_upper_lt_C0', 'Kzcirc_lt_Kz',
                                    'alpha0_gt_alphabound_3.8757', 'alpha0_minus_tau_gt_0.00305'))
res['pass'] = bool(fam_ok and res['nlines_is_39'] and consec and res['slab_internal_endpoints_equal_cut_set'] and not bad and
                   five and res['inner_integral_matches_printed'] and res['rank_one_block_is_-5/48'] and
                   res['outer_is_45/32'] and res['A_matches_printed'] and res['tau_printed_digits_consistent'] and
                   res['tau_upper_lt_taubound_3.87273'] and lat_ok)
res['seconds'] = time.time() - T0
json.dump({'item4_arithmetic': res, 'item5_alpha0': it5}, open(OUT, 'w'), indent=1)
print(json.dumps({'item4_arithmetic': res, 'item5_alpha0': it5}, indent=1))
