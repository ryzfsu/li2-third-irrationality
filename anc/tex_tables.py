#!/usr/bin/env python3
"""Generate the literal TeX tables from the exact frozen checker outputs."""
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent


def decimal_lower_endpoint(ball: str, significant_digits: int = 18) -> str:
    match = re.fullmatch(r"\[([^ ]+) \+/- ([^]]+)\]", ball)
    if not match:
        raise ValueError(f"unexpected Arb ball format: {ball!r}")
    with localcontext() as ctx:
        ctx.prec = 120
        lower = Decimal(match.group(1)) - Decimal(match.group(2))
        quantum = Decimal(1).scaleb(lower.adjusted() - significant_digits + 1)
        rounded = (lower / quantum).to_integral_value(rounding=ROUND_FLOOR) * quantum
        return format(rounded, "f")


def tex_fraction(value: F | str | int) -> str:
    value = F(value)
    sign = "-" if value < 0 else ""
    value = abs(value)
    if value.denominator == 1:
        return sign + str(value.numerator)
    return rf"{sign}\frac{{{value.numerator}}}{{{value.denominator}}}"


def exact_table(caption: str, label: str, columns: str, headers: list[str], rows: list[list[str]]) -> str:
    lines = [r"\begin{table}[p]", r"\centering", r"\scriptsize", r"\caption{" + caption + "}",
             r"\label{" + label + "}", r"\begin{tabular}{" + columns + "}", r"\hline",
             " & ".join(headers) + r" \\", r"\hline"]
    for row in rows:
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_tables(data_root: Path, output: Path) -> list[Path]:
    ver = data_root / "verification"
    arith = load(ver / "arithmetic.json")
    det = load(ver / "determinant.json")
    mod17 = load(ver / "mod17.json")
    analytic = load(ver / "analytic" / "result.json")
    audit = load(data_root / "analytic" / "audit.json")
    reference = load(data_root / "analytic" / "reference.json")

    summary = arith["summary"]
    inner = F(summary["inner_integral"])
    cap_inner = F(summary["cap_only_inner_integral"])
    outer = F(summary["outer_integral"])
    low = F(summary["native_low_strip"])
    total = F(summary["A_total"])
    integral_rows = [
        ["native strip $0<\\rho\\le 1/40$", tex_fraction(low)],
        ["inner rank integral $1/40\\le\\rho\\le5/4$", tex_fraction(inner)],
        ["inner dimension-cap integral", tex_fraction(cap_inner)],
        ["single rank-one saving", tex_fraction(cap_inner-inner)],
        ["outer integral $5/4<\\rho\\le11/4$", tex_fraction(outer)],
        ["$A_{40}$", tex_fraction(total)],
        ["$D_m/m^2$ coefficient", r"\frac{53}{48}"],
    ]
    r40_summary = exact_table("Exact R40 price components.", "tab:r40-summary", "lr",
                              ["component", "exact value"], integral_rows)
    outer_rows = [
        [r"$5/4<\rho<2$", r"$19/4-2\rho$", r"$9/8$"],
        [r"$2\le\rho\le11/4$", r"$11/4-\rho$", r"$9/32$"],
        [r"$5/4<\rho\le11/4$", r"$c_{\rm out}(\rho)$", tex_fraction(outer)],
    ]
    r40_outer = exact_table("Outer-prime exponent densities and integrals.", "tab:r40-outer", "lll",
                            [r"$\rho$ range", "density", "integral"], outer_rows)

    scalar = F(det["G4_scalar"])
    g4rows = [[str(j), str(value)]
              for j, value in enumerate(det["G4_primitive_coefficients_ascending"])]
    g4 = exact_table(r"Integer coefficients in the displayed formula for $G_4(T;3)$.",
                     "tab:g4-coefficients", "rl", ["degree $j$", "integer coefficient $q_j$"], g4rows)

    seed_poly = mod17["seed_polynomial_mod17"]
    mod_rows = [[str(row["W"]), str(row["seed_value_mod17"]), str(row["tower_det_mod17"])]
                for row in mod17["tower_determinants"]]
    mod17_table = exact_table(
        r"Seed polynomial $11+5W^2+4W^4$ and direct $153\times153$ tower determinants.",
        "tab:mod17-seed", "rrr", ["$W$", "$g(W)\bmod17$", r"$H_{17}(W)\bmod17$"], mod_rows)

    reference_parts = []
    for key, label, caption in [
        ("x", "tab:reference-x", r"Rational interval data for $\nu_x$ (weights are masses)."),
        ("y", "tab:reference-v", r"Rational interval data for $\nu_v$ (weights are masses)."),
        ("z_real", "tab:reference-z-real", r"Rational interval data for the real part of $\nu_z$."),
    ]:
        source_rows = reference[key]
        if key == "z_real":
            source_rows = sorted(source_rows, key=lambda row: F(row["a"]))
        rows = [[str(i), tex_fraction(F(row["a"])), tex_fraction(F(row["b"])), tex_fraction(F(row["w"]))]
                for i, row in enumerate(source_rows, 1)]
        reference_parts.append(exact_table(caption, label, "rlll",
                                           ["row", "$a$", "$b$", "mass"], rows))
    rows = [[str(i), tex_fraction(F(row["re"])), tex_fraction(F(row["im"])), tex_fraction(F(row["w"]))]
            for i, row in enumerate(reference["z_poisson"], 1)]
    reference_parts.append(exact_table(r"Rational centers and masses of the Poisson part of $\nu_z$.",
                                       "tab:reference-poisson", "rlll",
                                       ["row", r"$\Re c$", r"$\Im c$", "mass"], rows))
    reference_measure = "\n".join(reference_parts)

    interval_rows = []
    for i in range(max(len(reference["x"]), len(reference["y"]))):
        if i < len(reference["x"]):
            xrow = reference["x"][i]
            left = [str(i + 1), tex_fraction(F(xrow["a"])), tex_fraction(F(xrow["w"]))]
        else:
            left = ["", "", ""]
        if i < len(reference["y"]):
            vrow = reference["y"][i]
            right = [str(151 + i), tex_fraction(F(vrow["a"])), tex_fraction(F(vrow["w"]))]
        else:
            right = ["", "", ""]
        interval_rows.append(" & ".join(left + right) + r" \\")
    reference_intervals_rows = "\n".join(interval_rows) + "\n"
    z_real_rows = sorted(reference["z_real"], key=lambda row: F(row["a"]))
    z_rows = [f"{313 + i} & {tex_fraction(F(row['a']))} & {tex_fraction(F(row['w']))} " + r"\\"
              for i, row in enumerate(z_real_rows)]
    reference_zintervals_rows = "\n".join(z_rows) + "\n"
    poisson_rows = [
        f"{i} & {tex_fraction(F(row['re']))} & {tex_fraction(F(row['im']))} & {tex_fraction(F(row['w']))} " + r"\\"
        for i, row in enumerate(reference["z_poisson"], 1)
    ]
    reference_poisson_rows = "\n".join(poisson_rows) + "\n"

    coverage = audit["coverage"]
    order = [
        ("x", "[0,1]", "x"),
        ("v", "[0,1]", "v"),
        ("z_real", "[1,3]", "z real"),
        ("z_circle", "upper semicircle", "$|z|=2/5$"),
    ]
    cap_rows = []
    for key, support, shown in order:
        item = coverage[key]
        slack = item["minimum_margin_ball"]
        slack_decimal = decimal_lower_endpoint(slack)
        cap = {"x": "-9661/2000", "v": "-47421/10000", "z_real": "11/250", "z_circle": "11/250"}[key]
        cap_rows.append([shown, support, str(item["accepted_leaves"]), str(item["max_depth"]),
                         tex_fraction(F(cap)), f"${slack_decimal}$"])
    cap_table = exact_table("Full-support directed Arb cap certificate summary.", "tab:cap-leaf-stats",
                            "llrrrr", ["potential", "support", "leaves", "max depth", "cap", "minimum slack"],
                            cap_rows)

    # All irrational constants are enclosed elsewhere; the table prints the rational inputs.
    numeric_rows = [
        ["$Q(\\nu)$ lower bound", r"$-\frac{1298393}{250000}$"],
        [r"$C_0=(3/2)\log3$ upper bound", r"$\frac{1647919}{10^6}$"],
        [r"$K_x$", r"$-\frac{9661}{2000}$"],
        [r"$K_v$", r"$-\frac{47421}{10000}$"],
        [r"$K_z$", r"$\frac{11}{250}$"],
        [r"$\alpha_0$", r"$\frac{1936317}{500000}=3.872634$"],
        [r"$C$ error-budget bound", r"$10+3S<14514019$"],
        [r"$C$ used in the all-$m$ passage", r"$1.5\cdot10^7$"],
        [r"$\tau_{40}$", r"$3.86994289574041356037456595637820333081106094\ldots$"],
        [r"$\alpha-\tau_{40}$", r"$>0.00269$"],
    ]
    numeric = exact_table("Certified constants for the analytic and arithmetic comparison.",
                          "tab:numeric-constants", "lr", ["quantity", "value"], numeric_rows)

    generated = {
        "r40_integral_summary.tex": r40_summary,
        "r40_outer.tex": r40_outer,
        "g4_coefficients.tex": g4,
        "mod17_seed_block.tex": mod17_table,
        "reference_measure.tex": reference_measure,
        "reference_intervals_rows.tex": reference_intervals_rows,
        "reference_zintervals_rows.tex": reference_zintervals_rows,
        "reference_poisson_rows.tex": reference_poisson_rows,
        "cap_leaf_stats.tex": cap_table,
        "numeric_constants.tex": numeric,
    }
    output.mkdir(parents=True, exist_ok=True)
    written = []
    for name, contents in generated.items():
        path = output / name
        path.write_text(contents, encoding="utf-8")
        written.append(path)
    return written


def strip_comments(source: str) -> str:
    return re.sub(r"(?<!\\)%[^\n]*", "", source)


def read_tex(path: Path, active=()) -> str:
    path = path.resolve()
    if path in active:
        raise ValueError(f"cyclic TeX input: {path}")
    source = strip_comments(path.read_text(encoding="utf-8"))

    def include(match):
        child = path.parent / match.group(1)
        if not child.suffix:
            child = child.with_suffix(".tex")
        return read_tex(child, active + (path,))

    return re.sub(r"\\(?:input|include)\s*\{([^}]+)\}", include, source)


def load_tex(path: Path) -> str:
    """Expand a TeX source and its local input/include dependencies."""
    return read_tex(path)


def compact(cell: str) -> str:
    value = cell.replace("$", "")
    value = value.replace(r"\(", "").replace(r"\)", "")
    value = re.sub(r"\\[,;! ]", "", value)
    return re.sub(r"\s+", "", value)


def parse_fraction(cell: str) -> F:
    value = compact(cell)
    match = re.fullmatch(r"(-?)\\(?:t?frac)\{(-?\d+)\}\{(\d+)\}", value)
    if match:
        sign = -1 if match[1] == "-" else 1
        return F(sign * int(match[2]), int(match[3]))
    if re.fullmatch(r"[+-]?\d+(?:/\d+|\.\d+)?", value):
        return F(value)
    raise ValueError(f"not a literal rational cell: {cell!r}")


def table_rows(source: str, label: str) -> list[list[str]]:
    pattern = re.compile(r"\\begin\{(table|longtable)\}(.*?)\\end\{\1\}", re.S)
    matches = [m for m in pattern.finditer(source) if rf"\label{{{label}}}" in m[2]]
    if len(matches) != 1:
        raise ValueError(f"expected one table {label}, found {len(matches)}")
    body = matches[0][2]
    tabular = re.search(r"\\begin\{tabular\}(?:\[[^]]*\])?", body)
    if tabular:
        body = body[tabular.end():]
        depth = 0
        end = None
        for i, char in enumerate(body):
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        if end is None:
            raise ValueError(f"missing tabular column specification in {label}")
        body = body[end:].split(r"\end{tabular}", 1)[0]
    body = re.sub(r"\\(?:hline|toprule|midrule|bottomrule)\b", "", body)
    body = re.sub(r"\\(?:endfirsthead|endhead|endfoot|endlastfoot)\b", "", body)
    rows = []
    for raw in re.split(r"\\\\(?:\[[^]]*\])?", body):
        raw = raw.strip()
        if "&" in raw or r"\multicolumn" in raw:
            rows.append([cell.strip() for cell in raw.split("&")])
    return rows


def _longtable_body(source: str, label: str) -> str:
    """Return the unique longtable body carrying ``label``."""
    pattern = re.compile(
        rf"\\begin\{{longtable\}}(?:(?!\\end\{{longtable\}}).)*?"
        rf"\\label\{{{re.escape(label)}\}}(.*?)(?=\\end\{{longtable\}})", re.S)
    matches = [m.group(1) for m in pattern.finditer(source)]
    if len(matches) != 1:
        raise ValueError(f"expected one longtable {label}, found {len(matches)}")
    return matches[0]


def _row_chunks(body: str) -> list[list[str]]:
    """Split a longtable body into literal cell rows."""
    body = re.sub(r"\\(?:hline|endfirsthead|endhead|endfoot|endlastfoot)\b", "", body)
    rows = []
    for raw in re.split(r"\\\\(?:\[[^]]*\])?", body):
        raw = raw.strip()
        if "&" in raw or r"\multicolumn" in raw:
            rows.append([cell.strip() for cell in raw.split("&")])
    return rows


def parse_reference(source: str) -> dict:
    """Reconstruct the reference measure from the printed TeX tables.

    This is deliberately the only source used by the analytic checker.  The
    resulting dictionary is written only to a temporary directory for the
    closed-form Arb helpers; no JSON input is required in the ancillary tree.
    """
    ref = {"radius": "2/5", "source_a": 3, "x": [], "y": [],
           "z_real": [], "z_poisson": []}
    body = _longtable_body(source, "tab:reference-intervals")
    current = None
    for cells in _row_chunks(body):
        line = " ".join(cells)
        if "\\nu_x" in line or "\\nu_v" in line or "\\nu_z" in line:
            if "\\nu_x" in line:
                key = "x"
            elif "\\nu_v" in line:
                key = "y"
            else:
                key = "z_real"
            match = re.search(r"b_0\s*=\s*\$?([^$}]*)", line)
            if not match:
                raise ValueError(f"missing printed initial endpoint: {line}")
            current = {"key": key, "last": parse_fraction(match.group(1))}
            continue
        if current is None or not cells or not re.fullmatch(r"\$?[0-9]+\$?", compact(cells[0])):
            continue
        index = int(parse_fraction(cells[0]))
        entries = []
        for j in range(1, len(cells) - 1, 2):
            if not cells[j] or r"\multicolumn" in cells[j] or r"\multicolumn" in cells[j + 1]:
                break
            try:
                endpoint = parse_fraction(cells[j])
                weight = parse_fraction(cells[j + 1])
            except ValueError:
                break
            entries.append((index + len(entries), endpoint, weight))
        previous = current["last"]
        for _idx, endpoint, weight in entries:
            current["key"]
            ref[current["key"]].append({"a": str(previous), "b": str(endpoint),
                                        "w": str(weight)})
            previous = endpoint
        if entries:
            current["last"] = previous

    if [len(ref[k]) for k in ("x", "y", "z_real")] != [131, 143, 105]:
        raise ValueError("interval table did not reconstruct 131/143/105 rows")
    if any(F(row["b"]) != F(ref[k][i + 1]["a"])
           for k in ("x", "y", "z_real") for i, row in enumerate(ref[k][:-1])):
        raise ValueError("interval chains are not adjacent")
    body = _longtable_body(source, "tab:reference-poisson")
    displayed = []
    for cells in _row_chunks(body):
        if not cells:
            continue
        for j in range(0, len(cells) - 3, 4):
            if not re.fullmatch(r"\$?[0-9]+\$?", compact(cells[j])):
                continue
            try:
                displayed.append({"re": str(parse_fraction(cells[j + 1])),
                                  "im": str(parse_fraction(cells[j + 2])),
                                  "w": str(parse_fraction(cells[j + 3]))})
            except ValueError:
                continue
    positive = [row for row in displayed if F(row["im"]) > 0]
    real = [row for row in displayed if F(row["im"]) == 0]
    if len(positive) != 83 or len(real) != 1:
        raise ValueError("Poisson table must contain 83 upper rows and one real row")
    # The paper prints only the upper half.  Use a deterministic conjugate
    # expansion; energy and replay are invariant under this ordering.
    zpoisson = []
    for row in sorted(positive, key=lambda r: (F(r["re"]), F(r["im"]))):
        zpoisson.append(row)
        zpoisson.append({"re": row["re"], "im": str(-F(row["im"])), "w": row["w"]})
    zpoisson.append(real[0])
    ref["z_poisson"] = zpoisson
    return ref


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=HERE)
    parser.add_argument("--output", type=Path, default=HERE.parent.parent / "work" / "tables")
    args = parser.parse_args()
    paths = write_tables(args.data_root, args.output)
    for path in paths:
        print(path.name, path.stat().st_size)


if __name__ == "__main__":
    main()
