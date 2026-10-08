#!/usr/bin/env python3
"""Exact finite replay of the source determinant and the merged moments."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial, gcd, lcm, prod
import mpmath as mp
from pathlib import Path
import platform
import re
import sys
import time

from tex_tables import compact, load_tex, parse_fraction, table_rows


HERE = Path(__file__).resolve().parent


def det(matrix):
    matrix = [[F(x) for x in row] for row in matrix]
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("determinant matrix must be square")
    value = F(1)
    for k in range(n):
        pivot = next((i for i in range(k, n) if matrix[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            value = -value
        p = matrix[k][k]
        value *= p
        for i in range(k + 1, n):
            q = matrix[i][k] / p
            for j in range(k + 1, n):
                matrix[i][j] -= q * matrix[k][j]
            matrix[i][k] = F(0)
    return value


def seed_matrix(T: F):
    a = 3
    I, J, lam = range(2, 6), range(1, 10), range(1, 2)
    t, b = [F(0)], [F(0)]
    for j in range(1, 10):
        t.append(a * t[-1] + F(1, j))
        b.append(a * b[-1] + F(1, j * j))
    C = [[F(a**i - a**j, j - i) if i != j else F(0) for j in J] for i in I]
    B = [[(t[i] - t[j]) / (j - i) if i != j else b[i] for j in J] for i in I]
    A = [[a**i if i == j else 0 for j in J] for i in I]
    V = [[F(1, j + ell) for j in J] for ell in lam]
    matrix = C + [[B[r][c] - T * A[r][c] for c in range(9)] for r in range(4)] + V
    return matrix, t, b


def interpolate_at_integers(values):
    """Newton interpolation in the binomial basis, returned in monomials."""
    differences = [F(x) for x in values]
    falling = [F(1)]
    coefficients = [F(0)] * len(values)
    for k in range(len(values)):
        scale = differences[0] / factorial(k)
        for j, item in enumerate(falling):
            coefficients[j] += scale * item
        differences = [differences[j + 1] - differences[j] for j in range(len(differences) - 1)]
        nxt = [F(0)] * (len(falling) + 1)
        for j, item in enumerate(falling):
            nxt[j] -= k * item
            nxt[j + 1] += item
        falling = nxt
    return coefficients


def vandermonde(points):
    return prod((points[j] - points[i] for i in range(len(points)) for j in range(i + 1, len(points))), start=F(1))


def determinant_factor_test(m, s, k, h):
    a, L, e = F(3), m + h, m * (m - 1) // 2
    I, J = range(s + k, s + k + m), range(s, s + 2 * m + h)
    x = [F(j + 1, m + 2) for j in range(m)]
    v = [F(j + 1, L + 2) for j in range(L)]
    z = [F(3) + F(7 * (j + 1), 5 * (m + 1)) for j in range(m)]
    dc = det([[a**i * q ** (-i - 1) for q in z] for i in I])
    dx = det([[q**i for q in x] for i in I])
    dk = det([[F(1, q * t - 3) for t in v] for q in x]
             + [[t**ell for t in v] for ell in range(h)])
    dp = det([[q ** (s + p) for q in z + v] for p in range(m + L)])
    dx_expected = vandermonde(x) * prod((q ** (s + k) for q in x), start=F(1))
    dc_expected = ((-1) ** e * a ** (m * (s + k) + e) * vandermonde(z)
                   * prod((q ** (-s - k - m) for q in z), start=F(1)))
    denominator = prod((3 - q * t for q in x for t in v), start=F(1))
    cross = prod((q - t for q in z for t in v), start=F(1))
    dk_expected = ((-1) ** (m * L) * 3**e * vandermonde(x) * vandermonde(v)
                   * prod((q**h for q in x), start=F(1)) / denominator)
    dp_expected = ((-1) ** (m * L) * vandermonde(z) * vandermonde(v) * cross
                   * prod((q**s for q in z), start=F(1))
                   * prod((q**s for q in v), start=F(1)))
    if (dc, dx, dk, dp) != (dc_expected, dx_expected, dk_expected, dp_expected):
        raise AssertionError((m, s, k, h, "factor determinant identity"))
    integrand = (vandermonde(x)**2 * vandermonde(v)**2 * vandermonde(z)**2 * cross
                 * prod((q**(s + k + h) for q in x), start=F(1))
                 * prod((q**s for q in v), start=F(1))
                 * prod((q**(-m - k) for q in z), start=F(1)) / denominator)
    scalar = (-1)**e * 3**(m * (s + k) + 2 * e)
    if dc * dx * dk * dp != scalar * integrand:
        raise AssertionError((m, s, k, h, "merged integrand factor"))


def check_g4_table(tex, scalar, primitive):
    rows = table_rows(tex, "tab:g4-coefficients")
    parsed = []
    for row in rows[1:]:
        if len(row) != 2:
            raise AssertionError("G4 coefficient table must have two columns")
        j = int(parse_fraction(row[0]))
        qj = int(parse_fraction(row[1]))
        parsed.append((j, qj))
    expected = list(enumerate(primitive))
    if parsed != expected:
        raise AssertionError("printed G4 coefficient table differs from the exact determinant")
    match = re.search(
        r"G_4\(T\)\s*=\s*\\frac\{(\d+)\}\{(\d+)\}\s*\\sum_\{j=0\}\^4q_jT\^j",
        tex,
    )
    if not match or F(int(match[1]), int(match[2])) != scalar:
        raise AssertionError("printed G4 common scalar differs from the exact determinant")


def direct_chain_moment(d: int) -> mp.mpc:
    """Numerically integrate z^(d-1) over [3,1]+(log(3)/(2*pi*i))*gamma."""
    radius = mp.mpf(2) / 5
    coefficient = mp.log(3) / (2 * mp.pi * 1j)
    segment = mp.quad(lambda t: (3 - 2 * t) ** (d - 1) * (-2), [0, 1])
    circle = coefficient * mp.quad(
        lambda theta: (radius * mp.exp(1j * theta)) ** (d - 1)
        * (1j * radius * mp.exp(1j * theta)),
        [0, 2 * mp.pi],
    )
    return segment + circle


def run(tex_path: Path | None):
    started = time.monotonic()
    mp.mp.dps = 80
    values = [det(seed_matrix(F(T))[0]) for T in range(5)]
    coefficients = interpolate_at_integers(values)
    denominator = lcm(*(q.denominator for q in coefficients))
    numerators = [int(q * denominator) for q in coefficients]
    content = 0
    for q in numerators:
        content = gcd(content, abs(q))
    primitive = [q // content for q in numerators]
    scalar = F(content, denominator)
    if coefficients != [scalar * q for q in primitive]:
        raise AssertionError("primitive polynomial normalization")

    # The chain moment at d=0 is the exact cancellation -log(3)+log(3)=0;
    # for d != 0 the circle moment vanishes and the segment gives (1-3^d)/d.
    moment_tests = []
    max_direct_error = mp.mpf(0)
    for d in range(-16, 17):
        if d == 0:
            moment = F(0)
            expected = F(0)
            cancellation = "-log(3)+log(3)=0"
        else:
            moment = (1 - F(3)**d) / d
            expected = moment
            cancellation = "circle moment is zero for nonzero integer exponent"
        if moment != expected:
            raise AssertionError((d, "chain moment"))
        for i in range(0, 17):
            j = i + d
            if j < 0:
                continue
            c = F(0) if i == j else F(3**i - 3**j, j - i)
            from_moment = F(3**i) * moment
            if c != from_moment:
                raise AssertionError((i, j, "C row from chain moment"))
        direct = direct_chain_moment(d)
        target = mp.mpf(0) if d == 0 else (1 - mp.mpf(3) ** d) / d
        direct_error = abs(direct - target)
        max_direct_error = max(max_direct_error, direct_error)
        if direct_error >= mp.mpf("1e-65"):
            raise AssertionError((d, "direct chain quadrature", mp.nstr(direct_error, 12)))
        moment_tests.append({"d": d, "segment_plus_circle_moment": str(moment),
                             "direct_chain_error": mp.nstr(direct_error, 12),
                             "zero_case": cancellation if d == 0 else None})

    tests = [(m, s, k, h) for m in (1, 2, 3)
             for (s, k, h) in ((0, 2, 0), (1, 0, 1), (2, 1, 2))]
    for test in tests:
        determinant_factor_test(*test)

    if tex_path:
        tex = load_tex(tex_path)
        check_g4_table(tex, scalar, primitive)
    else:
        tex = ""
    exact = {
        "definition": "G_4(T;3)=det[C;B-TA;V], I=2..5, J=1..9, Lambda={1}",
        "matrix_size": 9,
        "degree": 4,
        "t_j": [str(seed_matrix(F(0))[1][j]) for j in range(10)],
        "b_j": [str(seed_matrix(F(0))[2][j]) for j in range(10)],
        "coefficients_ascending": [str(x) for x in coefficients],
        "rational_scalar": str(scalar),
        "primitive_integer_coefficients_ascending": primitive,
    }
    result = {
        "ok": True,
        "scope": "Exact seed determinant, chain moments and finite Cauchy--Vandermonde factor tests; uniform proofs remain in the paper.",
        "G4_scalar": str(scalar),
        "G4_primitive_coefficients_ascending": primitive,
        "G4_rational_coefficients_ascending": [str(x) for x in coefficients],
        "chain_moment_tests": moment_tests,
        "chain_quadrature_dps": 80,
        "chain_quadrature_max_error": mp.nstr(max_direct_error, 12),
        "factor_identity_tests": len(tests),
        "tex_sha256": hashlib.sha256(tex.encode()).hexdigest() if tex_path else None,
        "raw_tex_sha256": hashlib.sha256(tex_path.read_bytes()).hexdigest() if tex_path else None,
        "seconds": round(time.monotonic() - started, 6),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    return result, exact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "verification" / "determinant.json")
    args = parser.parse_args()
    result, seed = run(args.tex)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (HERE / "verification" / "seed_exact.json").write_text(json.dumps(seed, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
