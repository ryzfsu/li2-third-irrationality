#!/usr/bin/env python3
"""Item 1 + item 2: rebuild nu from the TeX tables and enclose Q(nu) in Arb.

Closed forms (paper, proof of Lemma 4.2, 'Energy'):
  interval-interval: G(u) = u^2/2 log|u| - 3u^2/4, G'' = log|u|
  interval-Poisson : primitive of log|t-c| (|t|>=r) and of log|r^2 - t cbar| - log r (|t|<=r), split at t=r
  Poisson-Poisson  : w w' (log|r^2 - c cbar'| - log r)
  D(nu_x,nu_v)     : (5/4)log3 - sum_{j<=N} M_j^x M_j^v/(j 3^j) - E_N, 0<=E_N<=(5/4)/(2(N+1)3^N), N=150
All inputs are exact rationals; all arithmetic is Arb balls (directed enclosures) at 320 bits.
"""
import sys, json, time
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx
import refmeasure

PREC = 320
ctx.prec = PREC
NE = 150
T0 = time.time()

d, chains, rows, poi = refmeasure.load(sys.argv[1])
OUT = sys.argv[2]
item1 = refmeasure.checks(d, chains, rows, poi)
print('item1', json.dumps(item1, indent=1), flush=True)


def q(x):
    return fmpq(x.numerator, x.denominator)


R = q(Fr(2, 5)); rA = arb(R); logr = rA.log()


def intervals(key):
    T, W = chains[key]['T'], chains[key]['W']
    return [(q(T[i]), q(T[i + 1]), q(W[i])) for i in range(len(W))]


IX, IV, IZ = intervals('x'), intervals('v'), intervals('z_int')
POI = [(q(a), q(b), q(w)) for a, b, w in poi]
POIA = [(acb(arb(a), arb(b)), arb(w)) for a, b, w in POI]

_gc = {}


def G(u):
    if u == 0:
        return arb(0)
    if u in _gc:
        return _gc[u]
    ua = arb(u)
    v = ua * ua * (abs(ua).log() / 2 - arb(fmpq(3, 4)))
    _gc[u] = v
    return v


def I_int_int(A, B):
    s = arb(0)
    for a0, a1, wa in A:
        Ha = arb(wa / (a1 - a0))
        for b0, b1, wb in B:
            Hb = arb(wb / (b1 - b0))
            s += Ha * Hb * (G(a1 - b0) - G(a0 - b0) - G(a1 - b1) + G(a0 - b1))
    return s


def fprim(u):  # complex u log u - u, principal branch
    return u * u.log() - u


def int_log_abs(t0, t1, s):
    """int_{t0}^{t1} log|t - s| dt for real t0<t1 and complex s not on [t0,t1] (segment avoids branch cut)."""
    return (fprim(acb(arb(t1)) - s) - fprim(acb(arb(t0)) - s)).real


def I_int_poi(A):
    s = arb(0)
    for t0, t1, w in A:
        pieces = []
        if t1 <= R:
            pieces.append((t0, t1, True))
        elif t0 >= R:
            pieces.append((t0, t1, False))
        else:
            pieces += [(t0, R, True), (R, t1, False)]
        H = arb(w / (t1 - t0))
        for p0, p1, inside in pieces:
            for (a, b, wc), (c, wca) in zip(POI, POIA):
                if inside:
                    # log|r^2 - t cbar| - log r = log|c| + log|t - r^2/cbar| - log r
                    sc = (rA * rA) / c.conjugate()
                    val = arb(p1 - p0) * (abs(c).log() - logr) + int_log_abs(p0, p1, sc)
                else:
                    val = int_log_abs(p0, p1, c)
                s += H * wca * val
    return s


def I_poi_poi():
    s = arb(0)
    r2 = acb(rA * rA)
    for c, w in POIA:
        for c2, w2 in POIA:
            s += w * w2 * (abs(r2 - c * c2.conjugate()).log() - logr)
    return s


def moments(A, N):
    M = []
    for j in range(1, N + 1):
        m = fmpq(0)
        for t0, t1, w in A:
            m += w * (t1 ** (j + 1) - t0 ** (j + 1)) / ((j + 1) * (t1 - t0))
        M.append(m)
    return M


def D_xv():
    Mx = moments(IX, NE); Mv = moments(IV, NE)
    s = fmpq(0)
    for j in range(1, NE + 1):
        s += Mx[j - 1] * Mv[j - 1] / (j * fmpq(3) ** j)
    tail = fmpq(5, 4) / (2 * (NE + 1) * fmpq(3) ** NE)
    base = arb(fmpq(5, 4)) * arb(3).log() - arb(s)
    return base - arb(0).union(arb(tail)), str(tail)


res = {'precision_bits': PREC, 'series_N': NE}
comp = {}
comp['I_xx'] = I_int_int(IX, IX); print('I_xx', comp['I_xx'], flush=True)
comp['I_vv'] = I_int_int(IV, IV); print('I_vv', comp['I_vv'], flush=True)
comp['I_int'] = I_int_int(IZ, IZ); print('I_int', comp['I_int'], flush=True)
comp['I_mix'] = I_int_poi(IZ); print('I_mix', comp['I_mix'], flush=True)
comp['I_Poi'] = I_poi_poi(); print('I_Poi', comp['I_Poi'], flush=True)
comp['I_vzint'] = I_int_int(IV, IZ); print('I_vzint', comp['I_vzint'], flush=True)
comp['I_vzPoi'] = I_int_poi(IV); print('I_vzPoi', comp['I_vzPoi'], flush=True)
comp['D_xv'], tail = D_xv(); print('D_xv', comp['D_xv'], flush=True)
Qv = (comp['I_xx'] + comp['I_vv'] + comp['I_int'] + 2 * comp['I_mix'] + comp['I_Poi'] +
      comp['I_vzint'] + comp['I_vzPoi'] - comp['D_xv'])
comp['Q'] = Qv
print('Q', Qv, flush=True)

tab = d['energy_table']
rows_out = {}
allc = True
for k, v in comp.items():
    lo, hi = arb(fmpq(Fr(tab[k][0]).numerator, Fr(tab[k][0]).denominator)), \
        arb(fmpq(Fr(tab[k][1]).numerator, Fr(tab[k][1]).denominator))
    contained = bool(v.lower() >= lo.lower() and v.upper() <= hi.upper())
    overlaps = bool(v.overlaps(lo.union(hi)))
    rows_out[k] = {'computed': v.str(25, radius=True), 'printed': tab[k],
                   'computed_inside_printed': contained, 'overlaps_printed': overlaps}
    allc &= contained
qm = Fr(d['macros']['qminus'])
qmA = arb(fmpq(qm.numerator, qm.denominator))
res['components'] = rows_out
res['D_tail_bound'] = tail
res['all_components_inside_printed_enclosures'] = bool(allc)
res['Q_lower'] = Qv.lower().str(25, radius=False)
res['q_minus'] = str(qm)
res['Q_gt_q_minus'] = bool(Qv > qmA)
res['Q_gt_q_minus_plus_5e-7'] = bool(Qv > qmA + arb(fmpq(5, 10 ** 7)))
res['Q_minus_q_minus'] = (Qv - qmA).str(15)
res['pass'] = bool(res['Q_gt_q_minus'] and res['Q_gt_q_minus_plus_5e-7'])
res['seconds'] = time.time() - T0
json.dump({'item1_reference': item1, 'item2_energy': res}, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
