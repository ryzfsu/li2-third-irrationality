# manuscript verification reports

These reports were generated in `/tmp` from the exact manuscript TeX whose
raw SHA-256 is recorded in each JSON report.  The four entry-point runs all
passed: the determinant report has `ok=true`, and the arithmetic, mod-17, and
analytic reports have `status=PASS`.  The analytic replay checked the printed
uniform subdivision (126464 leaves) and the independent directed-Arb caps
implementation; replay intermediates were kept outside this directory.

`VERIFICATION_SHA256SUMS` covers the reports and this README.  It is a
version-specific audit manifest and is deliberately separate from the stable
`../SCRIPTS_SHA256SUMS` manifest cited by the paper. The reports bind the
exact TeX hash that was checked; the script manifest does not.
