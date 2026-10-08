#!/usr/bin/env python3
"""Independent directed-ball verification of the printed Li2(1/3) caps.

Written from the paper formulas only (no reuse of any earlier replay code).

Potentials (nu = frozen reference measure, r = 2/5):
  U_i(q) = int log|q-t| dnu_i(t),  D_i(q) = int log|3-qt| dnu_i(t)
  F_x = 2U_x - D_v + (3/4) log q              on (0,1]
  F_v = 2U_v + U_z - D_x + (1/4) log q         on (0,1]
  F_z = 2U_z + U_v - (5/4) log|q|              on [1,3] and |q| = r
  U_z = U_zreal + P,  P = sum of Poisson potentials.

Leaf bound (real supports).  On a leaf [a,b] with centre c0, half width h,
split F = C + N, where C is the sum of the interval components whose
reference interval contains the leaf (these are convex on the leaf, since
U'' = H(t1-t0)/((q-t0)(t1-q)) > 0 inside), and N is everything else.  Then
for t in [a,b]
  C(t) <= chord of C,   N(t) <= N(c0) + N'(c0)(t-c0) + (1/2) B (t-c0)^2,
with B >= sup N''.  Hence
  sup F <= max_{t in {a,b}} [C(t) + N(c0) + N'(c0)(t-c0)] + (1/2) max(B,0) h^2.
Upper bounds for N'' are assembled per component:
  interval components outside the leaf: U'' <= 0 (concave);
  -D (series, all coefficients positive): -D'' = int t^2/(3-qt)^2 <= M/4;
  (3/4)log q, (1/4)log q: concave (<= 0);  -(5/4)log q: <= (5/4)/l^2;
  Poisson components: ball enclosure of the exact second derivative.

Circle |q|=r, q = r e^{i theta}, theta in [0,pi] (lower half by conjugation).
Method T: F(c0) + h|F'(c0)| + (1/2) max(B,0) h^2, B >= sup d^2F/dtheta^2,
  with d^2/dtheta^2 log|q-s| = Re(q s/(q-s)^2) and, for a uniform interval
  of density H on [t0,t1], the exact antiderivative
  H Re[ q(log(q-t1)-log(q-t0)) + q^2(1/(q-t1)-1/(q-t0)) ].
Method M (any leaf, used near theta = 0): every real-interval component
  (supports in [0,3], s >= 0) is increasing in theta on [0,pi] because
  |r e^{i theta}-s|^2 = r^2+s^2-2rs cos(theta); so it is bounded by its
  value at the right end, and only the Poisson part gets method T.
"""
import json, sys, time, hashlib, os
from flint import arb, acb, fmpq, ctx, arb_poly

PREC = int(os.environ.get('INDEP_PREC', '192'))
ctx.prec = PREC
NSER = 90
MAXEXTRA = 24

REF = sys.argv[1] if len(sys.argv) > 1 else 'reference.json'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'results.json'
ONLY = sys.argv[3].split(',') if len(sys.argv) > 3 else ['x', 'v', 'z_real', 'z_circle']

# Printed parameters (the macros and uniform depths are read from the manuscript).
CAPS = {'x': fmpq(-483148, 100000), 'v': fmpq(-474318, 100000),
        'z_real': fmpq(4318, 100000), 'z_circle': fmpq(4283, 100000)}
EPS0 = fmpq(1, 10**12)
DSHIFT = int(os.environ.get('INDEP_DEPTH_SHIFT', '0'))
DEPTH = {k: v + DSHIFT for k, v in {'x': 9, 'v': 7+1, 'z_real': 7, 'z_circle': 8}.items()}
NARCS = 32

def Q(s):
    s = str(s)
    if '/' in s:
        a, b = s.split('/')
        return fmpq(int(a), int(b))
    return fmpq(int(s))

raw = open(REF, 'rb').read()
SHA = hashlib.sha256(raw).hexdigest()
ref = json.loads(raw)
R = Q(ref['radius'])
assert R == fmpq(2, 5)
rA = arb(R)

def chain(lst):
    L = sorted([(Q(i['a']), Q(i['b']), Q(i['w'])) for i in lst], key=lambda t: t[0])
    T = [L[0][0]]
    W = []
    for a, b, w in L:
        assert a == T[-1] and b > a and w > 0
        T.append(b)
        W.append(w)
    H = [W[k] / (T[k + 1] - T[k]) for k in range(len(W))]
    return {'T': T, 'W': W, 'H': H, 'HA': [arb(h) for h in H], 'TA': [arb(t) for t in T],
            'M': sum(W, fmpq(0))}

MX = chain(ref['x']); MV = chain(ref['y']); MZ = chain(ref['z_real'])
POI = [(Q(c['re']), Q(c['im']), Q(c['w'])) for c in ref['z_poisson']]
POIA = [(acb(arb(a), arb(b)), arb(w), a, b, w) for a, b, w in POI]

# ---------- reference checks (exact) ----------
checks = {}
checks['mass_x'] = str(MX['M']); checks['mass_v'] = str(MV['M'])
checks['mass_z_real'] = str(MZ['M'])
mp = sum((w for _, _, w in POI), fmpq(0))
checks['mass_z_circle'] = str(mp); checks['mass_z_total'] = str(MZ['M'] + mp)
assert MX['M'] == 1 and MV['M'] == fmpq(5, 4) and MZ['M'] + mp == 1
assert MX['T'][0] >= 0 and MX['T'][-1] <= 1 and MV['T'][0] > 0 and MV['T'][-1] <= 1
assert MZ['T'][0] >= 1 and MZ['T'][-1] <= 3
assert all(a * a + b * b < R * R and w > 0 for a, b, w in POI)
S = set((a, b, w) for a, b, w in POI)
assert len(S) == len(POI) and all((a, -b, w) in S for a, b, w in POI)
checks['poisson_conjugate_symmetric'] = True
checks['poisson_centres_strictly_inside'] = True
checks['counts'] = {'x': len(MX['W']), 'v': len(MV['W']), 'z_real': len(MZ['W']), 'poisson': len(POI)}
checks['max_abs_c_over_r_sq'] = str(max((a * a + b * b) / (R * R) for a, b, _ in POI))

# ---------- -D series ----------
def neg_D_polys(m):
    # -D(q) = -M log3 + sum_j M_j q^j / (j 3^j) + tail, tail in [0, M/(2(N+1)3^N)]
    T, H = m['TA'], m['HA']
    coeffs = [arb(0)]
    dco = []
    for j in range(1, NSER + 1):
        Mj = arb(0)
        for k in range(len(H)):
            Mj += H[k] * (T[k + 1] ** (j + 1) - T[k] ** (j + 1))
        Mj = Mj / (j + 1)
        coeffs.append(Mj / (j * arb(3) ** j))
        dco.append(Mj / arb(3) ** j)        # coefficient of q^{j-1} in -D'
    M = arb(m['M'])
    t0 = M / (2 * (NSER + 1) * arb(3) ** NSER)
    t1 = M / (2 * arb(3) ** NSER)
    return {'p': arb_poly(coeffs), 'dp': arb_poly(dco), 'c': -M * arb(3).log(),
            'tail': arb(0).union(t0), 'dtail': arb(0).union(t1), 'B': M / 4}
NDV = neg_D_polys(MV); NDX = neg_D_polys(MX)

def negD(nd, q):
    return nd['c'] + nd['p'](q) + nd['tail'], nd['dp'](q) + nd['dtail']

# ---------- real-line primitives ----------
def gval(u):  # u log|u| - u for exact fmpq u
    if u == 0:
        return arb(0)
    ua = arb(u)
    return ua * abs(ua).log() - ua

def chain_center(m, q):
    """U(q), U'(q), per-interval value list at exact q (q not a breakpoint)."""
    T, H = m['T'], m['HA']
    g = []; L = []
    for t in T:
        u = q - t
        assert u != 0
        ua = arb(u); la = abs(ua).log()
        L.append(la); g.append(ua * la - ua)
    vals = [H[k] * (g[k] - g[k + 1]) for k in range(len(H))]
    ders = [H[k] * (L[k] - L[k + 1]) for k in range(len(H))]
    return sum(vals, arb(0)), sum(ders, arb(0)), vals, ders

def comp_val(m, k, q):
    return m['HA'][k] * (gval(q - m['T'][k]) - gval(q - m['T'][k + 1]))

def classify(m, a, b):
    """index of the reference interval containing [a,b], or None; assert no straddling."""
    T = m['T']; kin = None
    for k in range(len(T) - 1):
        t0, t1 = T[k], T[k + 1]
        if a >= t0 and b <= t1:
            kin = k
        elif b <= t0 or a >= t1:
            pass
        else:
            raise AssertionError('leaf straddles a reference interval')
    return kin

def poisson_real(q, inside):
    """sum of Poisson potentials and derivatives at exact real q."""
    val = arb(0); der = fmpq(0)
    R2 = R * R
    for a, b, w in POI:
        if inside:
            re = R2 - q * a; im = q * b
            d2 = re * re + im * im
            val += arb(w) * (arb(d2).log() / 2 - rA.log())
            der += w * (-(a * R2 - q * (a * a + b * b)) / d2)
        else:
            d2 = (q - a) ** 2 + b * b
            val += arb(w) * arb(d2).log() / 2
            der += w * (q - a) / d2
    return val, arb(der)

def poisson_real_B(qball, inside):
    """upper bound of sum w * second derivative of the Poisson potentials over a real ball."""
    q = acb(qball, 0); B = arb(0)
    r2 = acb(rA * rA, 0)
    for c, wa, a, b, w in POIA:
        if inside:
            cb = c.conjugate(); den = r2 - q * cb
            B += wa * (-(cb * cb) / (den * den)).real
        else:
            den = q - c
            B += wa * (-1 / (den * den)).real
    return B

# ---------- real-support leaf bound ----------
def leaf_real(sup, a, b):
    c0 = (a + b) / 2; h = (b - a) / 2
    ha = arb(h); qball = arb(a).union(arb(b))
    cap = arb(CAPS[sup])
    if sup == 'x':
        kx = classify(MX, a, b)
        U, Ud, vals, ders = chain_center(MX, c0)
        nd, ndd = negD(NDV, arb(c0))
        lg = arb(c0).log()
        N = 2 * U + nd + arb(fmpq(3, 4)) * lg
        Nd = 2 * Ud + ndd + arb(fmpq(3, 4)) / arb(c0)
        convex = []
        if kx is not None:
            N -= 2 * vals[kx]; Nd -= 2 * ders[kx]
            convex.append((2, MX, kx))
        B = NDV['B']
    elif sup == 'v':
        inside = b <= R
        assert inside or a >= R
        kv = classify(MV, a, b); kz = classify(MZ, a, b)
        assert kz is None
        Uv, Uvd, vals, ders = chain_center(MV, c0)
        Uz, Uzd, _, _ = chain_center(MZ, c0)
        P, Pd = poisson_real(c0, inside)
        nd, ndd = negD(NDX, arb(c0))
        N = 2 * Uv + Uz + P + nd + arb(fmpq(1, 4)) * arb(c0).log()
        Nd = 2 * Uvd + Uzd + Pd + ndd + arb(fmpq(1, 4)) / arb(c0)
        convex = []
        if kv is not None:
            N -= 2 * vals[kv]; Nd -= 2 * ders[kv]
            convex.append((2, MV, kv))
        B = NDX['B'] + poisson_real_B(qball, inside)
    elif sup == 'z_real':
        kz = classify(MZ, a, b); kv = classify(MV, a, b)
        assert kv is None
        Uz, Uzd, vals, ders = chain_center(MZ, c0)
        Uv, Uvd, _, _ = chain_center(MV, c0)
        P, Pd = poisson_real(c0, False)
        N = 2 * Uz + 2 * P + Uv - arb(fmpq(5, 4)) * arb(c0).log()
        Nd = 2 * Uzd + 2 * Pd + Uvd - arb(fmpq(5, 4)) / arb(c0)
        convex = []
        if kz is not None:
            N -= 2 * vals[kz]; Nd -= 2 * ders[kz]
            convex.append((2, MZ, kz))
        B = 2 * poisson_real_B(qball, False) + arb(fmpq(5, 4)) / (arb(a) * arb(a))
    else:
        raise ValueError
    ends = []
    for t, sgn in ((a, -1), (b, 1)):
        Cv = arb(0)
        for coef, m, k in convex:
            Cv += coef * comp_val(m, k, t)
        ends.append(Cv + N + sgn * Nd * ha)
    Bup = B.upper()
    Bpos = arb(Bup) if Bup > 0 else arb(0)
    if not Bpos.is_finite():
        return None
    bound = ends[0].max(ends[1]) + Bpos * ha * ha / 2
    return bound

# ---------- circle ----------
CONST_C = -arb(fmpq(5, 4)) * rA.log()

def circ_point(phi):
    """q = r e^{i pi phi} for exact rational phi (ball)."""
    s, c = arb.sin_cos_pi_fmpq(phi)
    return acb(rA * c, rA * s)

def interval_cplx(m, q, coef):
    """coef * sum_k H_k Re[f(q-t0)-f(q-t1)] and theta-derivative, q in open upper half plane."""
    T, H = m['TA'], m['HA']
    Ls = []; fs = []
    for t in T:
        u = q - t; L = u.log()
        Ls.append(L); fs.append(u * L - u)
    val = arb(0); der = arb(0)
    iq = acb(0, 1) * q
    for k in range(len(H)):
        val += H[k] * (fs[k] - fs[k + 1]).real
        der += H[k] * (iq * (Ls[k] - Ls[k + 1])).real
    return coef * val, coef * der

def interval_real_exact(m, q, coef):
    v = arb(0)
    for k in range(len(m['H'])):
        v += comp_val(m, k, q)
    return coef * v

def interval_B(m, qb, coef):
    """upper bound of coef * sum_k d^2/dtheta^2 (interval potential) over ball qb."""
    T, H, W = m['TA'], m['HA'], m['W']
    tot = arb(0)
    for k in range(len(H)):
        best = None
        try:
            u0 = qb - T[k]; u1 = qb - T[k + 1]
            e = (qb * (u1.log() - u0.log()) + qb * qb * (1 / u1 - 1 / u0)).real * H[k]
            if e.is_finite():
                best = e.upper()
        except Exception:
            pass
        try:
            s = acb(T[k].union(T[k + 1]), 0)
            d = qb - s
            e2 = (qb * s / (d * d)).real * arb(W[k])
            if e2.is_finite():
                u = e2.upper()
                best = u if best is None or u < best else best
        except Exception:
            pass
        if best is None:
            return None
        tot += arb(best)
    return coef * tot

def poisson_circ(q):
    val = arb(0); der = arb(0)
    iq = acb(0, 1) * q
    for c, wa, a, b, w in POIA:
        d = q - c
        val += wa * abs(d).log()
        der += wa * (iq / d).real
    return 2 * val, 2 * der

def poisson_circ_B(qb):
    tot = arb(0)
    for c, wa, a, b, w in POIA:
        d = qb - c
        tot += wa * (qb * c / (d * d)).real
    return 2 * tot

def realpart_at(phi):
    """real-interval part 2U_zreal + U_v at q = r e^{i pi phi}, phi in (0,1]."""
    if phi == 1:
        q = -R
        return interval_real_exact(MZ, q, 2) + interval_real_exact(MV, q, 1)
    q = circ_point(phi)
    vz, _ = interval_cplx(MZ, q, 2); vv, _ = interval_cplx(MV, q, 1)
    return vz + vv

def leaf_circle(pa, pb):
    pc = (pa + pb) / 2
    ha = arb.pi() * arb((pb - pa) / 2)
    q = circ_point(pc)
    vz, dz = interval_cplx(MZ, q, 2); vv, dv = interval_cplx(MV, q, 1)
    P, Pd = poisson_circ(q)
    th = arb.pi() * arb(pa).union(arb(pb))
    s, c = th.sin_cos()
    qb = acb(rA * c, rA * s)
    BP = poisson_circ_B(qb)
    res = []
    # method M
    if BP.is_finite():
        BPu = BP.upper(); BPp = arb(BPu) if BPu > 0 else arb(0)
        res.append(realpart_at(pb) + P + CONST_C + abs(Pd) * ha + BPp * ha * ha / 2)
    # method T
    if pa > 0:
        Bz = interval_B(MZ, qb, 2); Bv = interval_B(MV, qb, 1)
        if Bz is not None and Bv is not None and BP.is_finite():
            Bt = (Bz + Bv + BP).upper(); Btp = arb(Bt) if Bt > 0 else arb(0)
            F = vz + vv + P + CONST_C; Fd = dz + dv + Pd
            res.append(F + abs(Fd) * ha + Btp * ha * ha / 2)
    res = [x for x in res if x.is_finite()]
    if not res:
        return None
    best = res[0]
    for x in res[1:]:
        if x.upper() < best.upper():
            best = x
    return best

# ---------- driver ----------
def initial_pieces(sup):
    if sup == 'x':
        pts = sorted(set([EPS0] + MX['T']))
    elif sup == 'v':
        pts = sorted(set([EPS0] + MV['T'] + [R]))
    elif sup == 'z_real':
        pts = sorted(set([fmpq(1)] + MZ['T']))
    else:
        pts = [fmpq(k, NARCS) for k in range(NARCS + 1)]
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]

def run_support(sup, log):
    t0 = time.time()
    cap = arb(CAPS[sup])
    pieces = initial_pieces(sup)
    d = DEPTH[sup]; n = 2 ** d
    stack = []
    for (a, b) in pieces:
        for k in range(n):
            stack.append((a + (b - a) * k / n, a + (b - a) * (k + 1) / n, d))
    stack.reverse()
    nleaf = 0; nextra = 0; worst = None; worst_leaf = None; maxdepth = d
    min_margin = None; extra_leaves = []
    fn = leaf_circle if sup == 'z_circle' else (lambda a, b: leaf_real(sup, a, b))
    while stack:
        a, b, dep = stack.pop()
        bd = fn(a, b)
        ok = bd is not None and bd < cap
        if not ok:
            if dep - d >= MAXEXTRA:
                raise SystemExit('FAIL %s leaf [%s,%s] bound %s' % (sup, a, b, bd))
            m = (a + b) / 2
            stack.append((m, b, dep + 1)); stack.append((a, m, dep + 1))
            nextra += 1
            continue
        nleaf += 1
        if dep > d:
            extra_leaves.append(dep)
        maxdepth = max(maxdepth, dep)
        up = bd.upper()
        if worst is None or up > worst:
            worst = up; worst_leaf = (str(a), str(b), dep, bd.str(30))
        if nleaf % 2000 == 0:
            log('%s: %d leaves, worst so far %s, %.0fs' % (sup, nleaf, worst.str(20), time.time() - t0))
    margin = cap - worst
    assert margin > 0
    return {'initial_pieces': len(pieces), 'base_depth': d, 'base_leaves': len(pieces) * n,
            'accepted_leaves': nleaf, 'bisected_failed_nodes': nextra,
            'leaves_below_base_depth': len(extra_leaves), 'max_depth': maxdepth,
            'cap': str(CAPS[sup]), 'worst_upper_bound': worst.str(30, radius=False),
            'worst_leaf': worst_leaf, 'cap_minus_worst': margin.str(20),
            'pass': True, 'seconds': time.time() - t0}

def small_args():
    e = arb(EPS0)
    fx = -arb(fmpq(5, 4)) * (3 - e).log() + arb(fmpq(3, 4)) * e.log()
    fv = (3 + e).log() - (3 - e).log() + arb(fmpq(1, 4)) * e.log()
    return {'Fx_bound_on_(0,eps0]': fx.str(20), 'Fv_bound_on_(0,eps0]': fv.str(20),
            'Fx_below_cap': bool(fx < arb(CAPS['x'])), 'Fv_below_cap': bool(fv < arb(CAPS['v']))}

def alpha():
    qm = fmpq(-1298393, 250000)
    Kz = max(CAPS['z_real'], CAPS['z_circle'])
    a0 = qm - fmpq(1647919, 10**6) - CAPS['x'] - fmpq(5, 4) * CAPS['v'] - Kz
    l3 = arb(fmpq(3, 2)) * arb(3).log()
    return {'q_minus_input': str(qm), 'K_z_used': str(Kz), 'alpha0_exact': str(a0),
            'alpha0_decimal': arb(a0).str(12),
            'three_halves_log3': l3.str(25),
            'three_halves_log3_lt_1647919e-6': bool(l3 < arb(fmpq(1647919, 10**6))),
            'matches_claim_484473/125000': a0 == fmpq(484473, 125000)}

if __name__ == '__main__':
    T0 = time.time()
    def log(s):
        print('[%7.0fs] %s' % (time.time() - T0, s), flush=True)
    log('reference sha256 %s' % SHA)
    out = {'reference_sha256': SHA, 'precision_bits': PREC, 'depth_shift': DSHIFT, 'series_terms': NSER,
           'reference_checks': checks, 'small_arguments': small_args(), 'alpha': alpha(),
           'supports': {}}
    log(json.dumps(out['small_arguments'])); log(json.dumps(out['alpha']))
    for sup in ONLY:
        out['supports'][sup] = run_support(sup, log)
        log('%s DONE %s' % (sup, json.dumps(out['supports'][sup])))
        json.dump(out, open(OUT, 'w'), indent=1)
    out['elapsed_seconds'] = time.time() - T0
    out['status'] = 'PASS' if all(out['supports'][s]['pass'] for s in out['supports']) and \
        out['small_arguments']['Fx_below_cap'] and out['small_arguments']['Fv_below_cap'] else 'FAIL'
    json.dump(out, open(OUT, 'w'), indent=1)
    log('STATUS %s' % out['status'])
