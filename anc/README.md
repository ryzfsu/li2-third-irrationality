# Li2(1/3) ancillary verification

This directory is a paper-driven finite verification bundle for
the manuscript TeX supplied to each checker. The PDF is not a mathematical
input to these checks. No manuscript hash is embedded in the scripts: every
report records the raw SHA-256 of the TeX file actually checked. The same
bundle can be used for the v8 manuscript and for a later v9 manuscript
without changing the script sources.

The checkers parse the printed tables and macros from the TeX itself.  Any
JSON files created by the scripts are outputs in `verification/` or temporary
intermediates; no JSON file is required as a source of mathematical data.
All high-volume work is performed in a local temporary directory and the
temporary directory is removed at process exit.

## Entrypoints

* `check_determinant.py --tex MAIN --output verification/determinant.json`
  recomputes the finite G4 determinant, chain moments, and factor identities.
* `check_mod17.py --tex MAIN --output verification/mod17.json` recomputes the
  finite mod-17 seed and the three displayed tower determinants.
* `check_arith.py --tex MAIN --output verification/arithmetic.json` parses the
  R=8 line and slab tables and recomputes their exact rational arrangement,
  integral, tau bounds, divisor samples, and alpha assembly.
* `check_analytic.py --tex MAIN --replay --output verification/analytic.json`
  parses the reference measure, recomputes the energy with directed Arb, runs
  the printed uniform depths 9/8/7/8 with the Lemma 4.2 chord-plus-Taylor
  bound, and runs `indep_caps_params.py` as an independent directed-Arb
  cross-check.  Both implementations are required to have zero failed
  leaves.  Replay outputs never overwrite the TeX or this directory.

The helper scripts (`parse_tex.py`, `refmeasure.py`, `energy.py`, `arith.py`,
`caps.py`, `tex2refjson.py`, and `indep_caps_params.py`) are implementation
sources for those entrypoints.  The finite checks do not replace the written
uniform valuation, limiting, or nonvanishing proofs in the manuscript.

## Frozen results

The stable script manifest is `SCRIPTS_SHA256SUMS`; it contains only this
README and the Python checker/helper sources, so it is independent of the
manuscript version. The reports in `verification/` are a separate artifact:
`verification/VERIFICATION_SHA256SUMS` covers those reports and is regenerated
for the final TeX. The paper should cite the stable script-manifest hash, not
the verification-manifest hash.

For the current data, the paper-only caps run has 126464 leaves and no failures;
the independent caps run has the same accepted counts. The exact arithmetic
output records `A_8=902042411/177365760`, `D_m=(14,63,147,267)` for
`m=4,8,12,16`, and the stated positive alpha--tau gap.

The files in this directory are an ancillary audit artifact, not a public
claim of human refereeing.
