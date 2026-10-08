#!/usr/bin/env python3
"""Item 3: the caps (4.3) of Lemma 4.2, with the leaf bound exactly as described in the
proof of Lemma 4.2 of the manuscript, on the initial pieces and uniform depths
of the printed cap table (no adaptive refinement: a single failing leaf is a FAIL).

Real leaf [l,u], centre c0, half-width h0:
  F = F_C + F_N, F_C = the (coefficient-2) interval component whose chain interval contains [l,u]
  max F <= max_{q in {l,u}} ( F_C(q) + F_N(c0) + F_N'(c0)(q-c0) ) + h0^2 B / 2,
  B >= max(0, sup F_N''), assembled from: other interval components (concave, 0);
  (3/4)log q, (1/4)log q (concave, 0); -(5/4)log q (<= (5/4)/l^2);
  -D_v, -D_x (<= M/4); Poisson terms (ball enclosure of the exact second derivative).
  -D series with N' = 90, value tail [0, M/(2(N'+1)3^N')], derivative tail [0, M/(2*3^N')].
Circle arc theta = pi t, t in [t0,t1] subset [0,1] (upper half; lower half by conjugation):
  bound 1: F(th0) + h0|F'(th0)| + h0^2 B/2, B >= max(0, sup F''), B from Poisson ball
           enclosures of Re(q c/(q-c)^2) and, per interval component, the smaller upper end of
           H Re(q log((q-t1)/(q-t0)) + q^2/(q-t1) - q^2/(q-t0)) and of w Re(q s/(q-s)^2), s in [t0,t1];
  bound 2: real-interval terms at the right end of the arc + bound 1 for the Poisson terms alone;
  the smaller of the two; on the arc starting at theta=0 only bound 2.
All in python-flint Arb at 192 bits; the test is the certain comparison  bound < cap.
"""
import sys, json, time
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx, arb_poly
import refmeasure

ctx.prec = 192
NF = 90
T0 = time.time()
d, chains, rows, poi = refmeasure.load(sys.argv[1])
OUT = sys.argv[2]
ONLY = sys.argv[3].split(',') if len(sys.argv) > 3 else ['x', 'v', 'z', 'circ']
M = d['macros']


def q_(x):
    x = Fr(x)
    return fmpq(x.numerator, x.denominator)


CAP = {'x': q_(M['Kx']), 'v': q_(M['Kv']), 'z': q_(M['Kz']), 'circ': q_(M['Kzcirc'])}
DEPTH = {'x': 9, 'v': 8, 'z': 7, 'circ': 8}
NARCS = int(M['Narcs'])
EPS0 = fmpq(1, 10 ** 12)
r = fmpq(2, 5); rA = arb(r)


class Chain:
    def __init__(self, key):
        T, W = chains[key]['T'], chains[key]['W']
        self.T = [q_(t) for t in T]
        self.W = [q_(w) for w in W]
        self.H = [self.W[k] / (self.T[k + 1] - self.T[k]) for k in range(len(self.W))]
        self.TA = [arb(t) for t in self.T]
        self.HA = [arb(h) for h in self.H]
        self.mass = sum(self.W, fmpq(0))
        self.n = len(self.W)

    def find(self, l, u):
        """index k with [l,u] inside [T_k,T_{k+1}], None if [l,u] is outside the chain; error if straddling."""
        hit = None
        for k in range(self.n):
            a, b = self.T[k], self.T[k + 1]
            if l >= a and u <= b:
                hit = k
            elif not (u <= a or l >= b):
                raise RuntimeError('leaf straddles a chain breakpoint')
        return hit

    def comp_real(self, k, qx):
        """k-th component potential at exact real qx: H[g(qx-T_k) - g(qx-T_{k+1})], g(u) = u log|u| - u."""
        def g(u):
            if u == 0:
                return arb(0)
            ua = arb(u)
            return ua * abs(ua).log() - ua
        return self.HA[k] * (g(qx - self.T[k]) - g(qx - self.T[k + 1]))

    def real_at(self, qx):
        """sum potential, derivative and per-component (value, derivative) at exact real qx (not a breakpoint)."""
        logs = []; gs = []
        for t in self.T:
            u = qx - t
            ua = arb(u); L = abs(ua).log()
            logs.append(L); gs.append(ua * L - ua)
        vals = [self.HA[k] * (gs[k] - gs[k + 1]) for k in range(self.n)]
        ders = [self.HA[k] * (logs[k] - logs[k + 1]) for k in range(self.n)]
        return sum(vals, arb(0)), sum(ders, arb(0)), vals, ders

    def circ_at(self, qc):
        """value and theta-derivative at complex qc in the closed upper half plane (qc not real-negative-crossing)."""
        Ls = []; fs = []
        for t in self.TA:
            u = qc - t; L = u.log()
            Ls.append(L); fs.append(u * L - u)
        iq = acb(0, 1) * qc
        val = arb(0); der = arb(0)
        for k in range(self.n):
            val += self.HA[k] * (fs[k] - fs[k + 1]).real
            der += self.HA[k] * (iq * (Ls[k] - Ls[k + 1])).real
        return val, der

    def circ_B(self, qb):
        """upper bound of sum_k d^2/dtheta^2 (component k) over the arc ball qb, min of the two enclosures."""
        tot = arb(0)
        for k in range(self.n):
            cand = []
            u0 = qb - self.TA[k]; u1 = qb - self.TA[k + 1]
            e1 = self.HA[k] * (qb * (u1.log() - u0.log()) + qb * qb / u1 - qb * qb / u0).real
            if e1.is_finite():
                cand.append(e1.upper())
            s = acb(self.TA[k].union(self.TA[k + 1]))
            dd = qb - s
            e2 = arb(self.W[k]) * (qb * s / (dd * dd)).real
            if e2.is_finite():
                cand.append(e2.upper())
            if not cand:
                return None
            tot += min(cand)
        return tot


CX, CV, CZ = Chain('x'), Chain('v'), Chain('z_int')
POI = [(acb(arb(q_(a)), arb(q_(b))), arb(q_(w))) for a, b, w in poi]


def negD_poly(ch):
    co = [arb(0)]; dco = []
    for j in range(1, NF + 1):
        mj = fmpq(0)
        for k in range(ch.n):
            mj += ch.W[k] * (ch.T[k + 1] ** (j + 1) - ch.T[k] ** (j + 1)) / ((j + 1) * (ch.T[k + 1] - ch.T[k]))
        co.append(arb(mj / (j * fmpq(3) ** j)))
        dco.append(arb(mj / fmpq(3) ** j))
    m = ch.mass
    return {'p': arb_poly(co), 'dp': arb_poly(dco), 'const': -arb(m) * arb(3).log(),
            'vt': arb(0).union(arb(m / (2 * (NF + 1) * fmpq(3) ** NF))),
            'dt': arb(0).union(arb(m / (2 * fmpq(3) ** NF))), 'B': arb(m / 4)}


ND_V = negD_poly(CV)   # -D_v, enters F_x
ND_X = negD_poly(CX)   # -D_x, enters F_v


def negD(nd, qa):
    return nd['const'] + nd['p'](qa) + nd['vt'], nd['dp'](qa) + nd['dt']


def poisson_real(qx, inside):
    """sum_c w * potential and derivative at exact real qx; inside: |qx|<=r formula."""
    qa = acb(arb(qx)); v = arb(0); dv = arb(0)
    r2 = acb(rA * rA)
    for c, w in POI:
        if inside:
            den = r2 - qa * c.conjugate()
            v += w * (abs(den).log() - rA.log())
            dv += w * (-c.conjugate() / den).real
        else:
            den = qa - c
            v += w * abs(den).log()
            dv += w * (1 / den).real
    return v, dv


def poisson_real_B(qball, inside):
    qa = acb(qball); B = arb(0); r2 = acb(rA * rA)
    for c, w in POI:
        if inside:
            den = r2 - qa * c.conjugate()
            B += w * (-(c.conjugate() ** 2) / (den * den)).real
        else:
            den = qa - c
            B += w * (-1 / (den * den)).real
    return B


def pos(B):
    up = B.upper()
    return arb(up) if up > 0 else arb(0)


def leaf_real(sup, l, u):
    c0 = (l + u) / 2; h0 = arb((u - l) / 2); ca = arb(c0)
    qball = arb(l).union(arb(u))
    conv = None
    if sup == 'x':
        k = CX.find(l, u)
        U, Ud, vals, ders = CX.real_at(c0)
        nd, ndd = negD(ND_V, ca)
        FN = 2 * U + nd + arb(fmpq(3, 4)) * ca.log()
        FNd = 2 * Ud + ndd + arb(fmpq(3, 4)) / ca
        if k is not None:
            FN -= 2 * vals[k]; FNd -= 2 * ders[k]; conv = (CX, k)
        B = ND_V['B']
    elif sup == 'v':
        inside = u <= r
        if not inside and l < r:
            raise RuntimeError('v leaf straddles r')
        k = CV.find(l, u)
        if CZ.find(l, u) is not None:
            raise RuntimeError
        U, Ud, vals, ders = CV.real_at(c0)
        Uz, Uzd, _, _ = CZ.real_at(c0)
        P, Pd = poisson_real(c0, inside)
        nd, ndd = negD(ND_X, ca)
        FN = 2 * U + Uz + P + nd + arb(fmpq(1, 4)) * ca.log()
        FNd = 2 * Ud + Uzd + Pd + ndd + arb(fmpq(1, 4)) / ca
        if k is not None:
            FN -= 2 * vals[k]; FNd -= 2 * ders[k]; conv = (CV, k)
        B = ND_X['B'] + poisson_real_B(qball, inside)
    else:  # z on [1,3]
        k = CZ.find(l, u)
        if CV.find(l, u) is not None:
            raise RuntimeError
        U, Ud, vals, ders = CZ.real_at(c0)
        Uv, Uvd, _, _ = CV.real_at(c0)
        P, Pd = poisson_real(c0, False)
        FN = 2 * U + 2 * P + Uv - arb(fmpq(5, 4)) * ca.log()
        FNd = 2 * Ud + 2 * Pd + Uvd - arb(fmpq(5, 4)) / ca
        if k is not None:
            FN -= 2 * vals[k]; FNd -= 2 * ders[k]; conv = (CZ, k)
        B = 2 * poisson_real_B(qball, False) + arb(fmpq(5, 4)) / (arb(l) ** 2)
    ends = []
    for qe, sg in ((l, -1), (u, 1)):
        FC = 2 * conv[0].comp_real(conv[1], qe) if conv else arb(0)
        ends.append(FC + FN + sg * FNd * h0)
    return ends[0].max(ends[1]) + pos(B) * h0 * h0 / 2


CONST = -arb(fmpq(5, 4)) * rA.log()


def cpoint(t):
    s, c = arb.sin_cos_pi_fmpq(t)
    return acb(rA * c, rA * s)


def poisson_circ(qc):
    v = arb(0); dv = arb(0); iq = acb(0, 1) * qc
    for c, w in POI:
        dd = qc - c
        v += w * abs(dd).log(); dv += w * (iq / dd).real
    return 2 * v, 2 * dv


def poisson_circ_B(qb):
    B = arb(0)
    for c, w in POI:
        dd = qb - c
        B += w * (qb * c / (dd * dd)).real
    return 2 * B


def realterms_at(t):
    """2 U_zint + U_v at r e^{i pi t}."""
    if t == 1:
        qx = -r
        vz = sum((CZ.comp_real(k, qx) for k in range(CZ.n)), arb(0))
        vv = sum((CV.comp_real(k, qx) for k in range(CV.n)), arb(0))
        return 2 * vz + vv
    qc = cpoint(t)
    return 2 * CZ.circ_at(qc)[0] + CV.circ_at(qc)[0]


def leaf_circ(t0, t1):
    tc = (t0 + t1) / 2
    h0 = arb.pi() * arb((t1 - t0) / 2)
    qc = cpoint(tc)
    th = arb.pi() * arb(t0).union(arb(t1))
    sn, cs = th.sin_cos()
    qb = acb(rA * cs, rA * sn)
    P, Pd = poisson_circ(qc)
    BP = poisson_circ_B(qb)
    cands = []
    # bound 2
    b2 = realterms_at(t1) + CONST + P + abs(Pd) * h0 + pos(BP) * h0 * h0 / 2
    cands.append(('bound2', b2))
    if t0 > 0:
        vz, dz = CZ.circ_at(qc); vv, dv = CV.circ_at(qc)
        Bz = CZ.circ_B(qb); Bv = CV.circ_B(qb)
        if Bz is not None and Bv is not None:
            F = 2 * vz + vv + P + CONST
            Fd = 2 * dz + dv + Pd
            b1 = F + abs(Fd) * h0 + pos(2 * Bz + Bv + BP) * h0 * h0 / 2
            cands.append(('bound1', b1))
    cands = [c for c in cands if c[1].is_finite()]
    best = min(cands, key=lambda c: c[1].upper())
    return best[1], best[0]


def pieces(sup):
    if sup == 'x':
        pts = sorted(set([EPS0, fmpq(1)] + [t for t in CX.T if EPS0 < t < 1]))
    elif sup == 'v':
        pts = sorted(set([EPS0, fmpq(1), r] + [t for t in CV.T if EPS0 < t < 1]))
    elif sup == 'z':
        pts = sorted(set([fmpq(1), fmpq(3)] + [t for t in CZ.T if 1 < t < 3]))
    else:
        pts = [fmpq(k, NARCS) for k in range(NARCS + 1)]
    return list(zip(pts, pts[1:]))


def run(sup, log):
    t0 = time.time()
    P = pieces(sup); dpt = DEPTH[sup]; n = 2 ** dpt
    cap = arb(CAP[sup])
    worst = None; wleaf = None; nfail = 0; nleaf = 0; used = {}
    fails = []
    for A, Bp in P:
        for j in range(1, n + 1):
            l = A + (Bp - A) * (j - 1) / n
            u = A + (Bp - A) * j / n
            if sup == 'circ':
                bd, how = leaf_circ(l, u)
                used[how] = used.get(how, 0) + 1
            else:
                bd = leaf_real(sup, l, u)
            nleaf += 1
            if not (bd < cap):
                nfail += 1
                if len(fails) < 20:
                    fails.append([str(l), str(u), bd.str(15)])
            up = bd.upper()
            if worst is None or up > worst:
                worst = up; wleaf = [str(l), str(u), bd.str(25)]
            if nleaf % 5000 == 0:
                log('%s %d leaves, worst %s, fails %d, %.0fs' % (sup, nleaf, worst.str(12), nfail, time.time() - t0))
    wq = worst.upper()
    return {'initial_pieces': len(P), 'depth': dpt, 'leaves': nleaf, 'failing_leaves': nfail,
            'fail_examples': fails, 'cap': str(CAP[sup]), 'worst_upper_end': wq.str(20, radius=False),
            'worst_leaf': wleaf, 'cap_minus_worst': (cap - worst).str(12), 'bound_used': used,
            'pass': nfail == 0, 'seconds': time.time() - t0}


def small():
    e = arb(EPS0)
    fx = -arb(fmpq(5, 4)) * (3 - e).log() + arb(fmpq(3, 4)) * e.log()
    fv = (3 + e).log() - (3 - e).log() + arb(fmpq(1, 4)) * e.log()
    # premises: |q - t| <= 1 on x,v supports; |q-t| <= 3+e0 for t in [1,3]; Poisson inside value <= log(r+e0) <= log(3+e0)
    return {'Fx_bound': fx.str(15), 'Fv_bound': fv.str(15),
            'Fx_bound_lt_-22.09': bool(fx < arb(fmpq(-2209, 100))), 'Fv_bound_lt_-6.90': bool(fv < arb(fmpq(-69, 10))),
            'Fx_below_Kx': bool(fx < arb(CAP['x'])), 'Fv_below_Kv': bool(fv < arb(CAP['v'])),
            'poisson_inside_bound_log(r+e0)_lt_log(3+e0)': bool((rA + e).log() < (3 + e).log())}


if __name__ == '__main__':
    def log(s):
        print('[%6.0fs] %s' % (time.time() - T0, s), flush=True)
    out = {'tex_sha256': d['tex_sha256'], 'prec_bits': 192, 'NF': NF, 'depths': DEPTH,
           'caps': {k: str(v) for k, v in CAP.items()}, 'small_arguments': small(), 'supports': {}}
    log(json.dumps(out['small_arguments']))
    for s in ONLY:
        out['supports'][s] = run(s, log)
        log('%s DONE %s' % (s, json.dumps(out['supports'][s])))
        json.dump(out, open(OUT, 'w'), indent=1)
    out['total_leaves'] = sum(v['leaves'] for v in out['supports'].values())
    out['pass'] = all(v['pass'] for v in out['supports'].values()) and out['small_arguments']['Fx_below_Kx'] and \
        out['small_arguments']['Fv_below_Kv']
    out['seconds'] = time.time() - T0
    json.dump(out, open(OUT, 'w'), indent=1)
    log('STATUS %s' % ('PASS' if out['pass'] else 'FAIL'))
