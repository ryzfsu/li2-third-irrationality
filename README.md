# Irrationality of the dilogarithm at 1/3

Ruiyi Zhang

The manuscript proves that `Li_2(1/3)` is irrational.

Irrationality of `Li_2(1/s)` was known for integers `s <= -5` and `s >= 6` (Hata; Rhin and
Viola), and Prellberg recently proved it for `s = -4, 5, -3`. The proof defines polynomials
`G_m(T)` with rational coefficients by an interpolation determinant of size `9n`
(`m = 4n`). The rank of this determinant modulo primes `p` with `m/8 < p < 5m/4` bounds the
denominators of the coefficients:

    limsup log K_m / m^2 <= tau < 3.87273.

A concavity inequality for a logarithmic energy, combined with the inversion `v -> 3/v`,
gives

    log |G_m(Li_2(1/3))| <= -3.8757 m^2 + C m log m.

A congruence modulo 17 shows that `G_m(xi) != 0` for every rational `xi` and
`m = 4 * 17^t`, `t` large. Interval arithmetic is used only for the energy and the
potentials of a printed reference measure (Lemma 4.2), on a uniform subdivision with
126,464 leaves. The other tables can be checked by hand.

## Contents

- `main.tex`, `main.pdf` (24 pages)
- `anc/` - the four checking programs and their helpers (Appendix D of the paper), with
  `anc/SCRIPTS_SHA256SUMS`, whose SHA-256 is printed in the paper. See `anc/README.md`.
- `anc/verification/` - the reports of the four checks run on this `main.tex`, with
  `VERIFICATION_SHA256SUMS`.
- `SHA256SUMS`

## Build

    pdflatex main.tex
    pdflatex main.tex

## Checks

Python 3 with `python-flint` (Arb) for the analytic check:

    cd anc
    python3 check_determinant.py --tex ../main.tex --output /tmp/determinant.json
    python3 check_mod17.py --tex ../main.tex --output /tmp/mod17.json
    python3 check_arith.py --tex ../main.tex --output /tmp/arithmetic.json
    python3 check_analytic.py --tex ../main.tex --replay --output /tmp/analytic.json
