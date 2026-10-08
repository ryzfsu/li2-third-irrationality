#!/usr/bin/env python3
"""Recompute the 153-by-153 mod-17 tower checks and read the printed table."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import re
import sys
import time

from mod17_replay import check_tower, native_seed, qmod, rational_det
from tex_tables import compact, load_tex, parse_fraction, table_rows


HERE = Path(__file__).resolve().parent
P = 17


def verify_printed_seed(tex: str) -> None:
    seed_label = r"\label{eq:g-mod-17}"
    start = tex.index(seed_label)
    end = tex.index(r"\end{equation}", start)
    statement = compact(tex[start:end])
    if not re.search(r"G_4\(T\)\\equivg\(T\):=11\+5T\^2\+4T\^4\\pmod\{17\}", statement):
        raise AssertionError("printed seed congruence differs from 11+5T^2+4T^4")

    start = tex.index(r"\label{eq:g4-mod-matrix}")
    end = tex.index(r"\end{equation}", start)
    block = tex[start:end]
    match = re.search(r"\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", block, re.S)
    if not match:
        raise AssertionError("printed mod-17 seed matrix is missing")
    rows = [row.strip() for row in match.group(1).split(r"\\") if row.strip()]
    printed = []
    for row in rows:
        entries = []
        for cell in row.split("&"):
            value = compact(cell)
            item = re.fullmatch(r"(\d+)(?:\+(\d+)T)?", value)
            if not item:
                raise AssertionError(f"unsupported mod-17 matrix entry {cell!r}")
            entries.append((int(item[1]) % P, int(item[2] or 0) % P))
        printed.append(entries)
    base = native_seed(F(0))[0]
    shifted = native_seed(F(1))[0]
    expected = [[(qmod(base[i][j], P), (qmod(shifted[i][j], P) - qmod(base[i][j], P)) % P)
                 for j in range(9)] for i in range(9)]
    if printed != expected:
        raise AssertionError("printed mod-17 seed matrix differs from the exact source matrix")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "verification" / "mod17.json")
    args = parser.parse_args()
    start = time.monotonic()
    # Degree at most four and values at five points determine the polynomial.
    # The defining determinant is evaluated independently by Fraction Gaussian elimination.
    values = [qmod(rational_det(native_seed(F(T))[0]), P) for T in range(5)]
    # Convert the Newton values at 0,...,4 to the expected sparse coefficients
    # by direct interpolation over F_17 at every field point.
    polynomial = [0] * 5
    for w, val in enumerate(values):
        basis = [1]
        denominator = 1
        for u in range(5):
            if u == w:
                continue
            next_basis = [0] * (len(basis) + 1)
            for j, c in enumerate(basis):
                next_basis[j] = (next_basis[j] - u * c) % P
                next_basis[j + 1] = (next_basis[j + 1] + c) % P
            basis = next_basis
            denominator = denominator * (w - u) % P
        scale = val * pow(denominator, -1, P) % P
        for j, c in enumerate(basis):
            polynomial[j] = (polynomial[j] + scale * c) % P
    expected_seed = [11, 0, 5, 0, 4]
    if polynomial != expected_seed:
        raise AssertionError(("seed polynomial mod 17", polynomial))

    replay = check_tower(17, 1, True, [0, 1, 2])
    checks = replay["native_dense_modular_determinant_checks"]
    if [(q["W"], q["native_scaled_determinant_mod_p"], q["seed_value_mod_p"]) for q in checks] != [
        (0, 11, 11), (1, 3, 3), (2, 10, 10)
    ]:
        raise AssertionError("direct tower determinant values")

    if args.tex:
        tex = load_tex(args.tex)
        verify_printed_seed(tex)
        if r"\label{tab:mod17-seed}" in tex:
            rows = table_rows(tex, "tab:mod17-seed")
            actual = [[int(parse_fraction(c)) for c in row] for row in rows[1:]]
            if actual != [[0, 11, 11], [1, 3, 3], [2, 10, 10]]:
                raise AssertionError("printed mod-17 seed block table")

    result = {
        "status": "PASS",
        "scope": "Exact finite seed and direct 153x153 modular determinants; the all-t block identity is proved in the manuscript.",
        "prime": P,
        "seed_polynomial_mod17": "11+5W^2+4W^4",
        "seed_coefficients_ascending_mod17": polynomial,
        "tower_determinants": [
            {"W": row["W"], "seed_value_mod17": row["seed_value_mod_p"],
             "tower_det_mod17": row["native_scaled_determinant_mod_p"]}
            for row in checks
        ],
        "matrix_order": 153,
        "coefficient_pair_checks": replay["coefficient_pair_checks"],
        "scaled_sequence_checks": replay["scaled_sequence_checks"],
        "tex_sha256": hashlib.sha256(load_tex(args.tex).encode()).hexdigest() if args.tex else None,
        "raw_tex_sha256": hashlib.sha256(args.tex.read_bytes()).hexdigest() if args.tex else None,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "seconds": round(time.monotonic() - start, 6),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
