# GD and Nesterov comparison

This directory supplies the code, fixed configuration and accepted records
for the four-method comparison in the working manuscript. CHB and HB are
read from the original records in `data/accepted/`; the additional methods
are gradient descent (GD) and quadratically tuned Nesterov (NAG-Q).

The completed study has 24 full trajectories and four two-step prefix checks:
12,288 main gradient updates plus eight prefix updates. These are recorded
calculations. Publishing and checking this distribution does not rerun them.

## Findings and design

See [the results](RESULTS.md), [protocol](protocol.md),
[configuration](config.json), [endpoint table](reference/endpoint_metrics.csv)
and [paired comparisons](reference/CHB_comparisons.csv).

Two synthetic quadratics (R = 1/2 and 1/8) and the existing raw Yacht
least-squares problem share their original starts and 512-query horizon.
Both exact gradients and the scalar-half relative-error policy are retained.
The 120-digit results are primary; the 80-digit repetition checks arithmetic.

CHB has smaller terminal distances on the synthetic examples. On Yacht,
NAG-Q has smaller objective gaps than CHB at every positive query under both
shared policies. A CHB distance advantage at query 512 reverses at query 511;
the maximum distance in the adjacent pair favors NAG-Q. The comparison does
not establish general empirical superiority or a worst-case error threshold.

## Verify saved records

Use an environment with `mpmath==1.3.0`, as described in the root README.
Every command block starts from the repository's parent directory.

```bash
cd periodic-momentum
.venv/bin/python -B comparisons/gd-nag/run_pilot.py verify
```

This checks the distribution manifest, original frozen baseline inputs,
28 accepted record receipts, all 24 full work-unit identities, 513 rows per
trajectory, the 12,296-update ledger and the 80/120-digit agreement. Matrix
and relative-error acceptance values are checked from the stored results.
No optimizer update, prefix, QR solve or eigensolver is called.

## Recreate tables and plots from saved records

```bash
cd periodic-momentum
.venv/bin/python -B comparisons/gd-nag/report.py --out outputs/gd_nag_reports
Rscript --vanilla comparisons/gd-nag/plot.R outputs/gd_nag_reports outputs/gd_nag_figures
```

The report uses saved records only. All four reporting CSVs were checked
against `reference/` and matched byte for byte after portable-path adaptation.
Outputs must be new directories. Rendering requires base R with Cairo and
creates SVG/PNG figures. The accepted figures are also in [figures/](figures/).
`D_pair` is the larger normalized distance in each adjacent odd/even pair.

Synthetic CHB/HB records retain full norm logs but objective values only at
saved states; missing objective traces are left blank. Yacht retains all
four methods' full objective traces. The synthetic norm logs were stored
as float64, and their comparisons use the recorded 1e-12 log tie tolerance;
the Yacht tolerance is 1e-30.

## Manuscript figures

The current four-method [Figure 2](manuscript_figures/Fig2.pdf) and
[Figure 3](manuscript_figures/Fig3.pdf) are supplied as the original vector
PDFs with embedded fonts. Figure 1 is unchanged in the original material.
Their plotting commands use the same saved reporting CSVs:

```bash
cd periodic-momentum
Rscript --vanilla comparisons/gd-nag/render_manuscript_figures.R outputs/gd_nag_reports outputs/gd_nag_manuscript_figures
```

This redraws only graphics. The stored PDFs include a structure-only removal
of redundant page transparency groups from opaque Cairo output; drawing
commands are unchanged. A redraw can differ in PDF metadata while retaining
the plotted results. No numerical method is called by the rendering script.

## Numerical implementation and replay

`adapter.py` retains the numerical functions used in the accepted study.
Publication changes are limited to repository-relative paths, separate
accepted/replay directories and manifest verification. The original 28
record files are unchanged. Configuration values are unchanged except for
the non-scientific repository-path field. See [provenance.json](provenance.json).

GD uses alpha = 2/(L + mu), beta = 0. NAG-Q uses alpha = 4/(3L + mu)
and beta = (sqrt(3 kappa + 1) - 2)/(sqrt(3 kappa + 1) + 2), with
kappa = L/mu. NAG-Q evaluates the perturbed gradient at its look-ahead point.
Coefficients are fixed before the comparison and are not retuned for errors.
These quadratic tunings are from [Lessard, Recht and Packard (2016),
Proposition 1](https://doi.org/10.1137/15M1009597).

The following command performs a fresh numerical replay; it is distinct
from the saved-record commands above:

```bash
cd periodic-momentum
.venv/bin/python -B comparisons/gd-nag/run_pilot.py run
```

It writes only `comparisons/gd-nag/outputs/run_001/`, refuses overwriting,
and enforces the fixed 24-trajectory/12,296-update and 120-second budgets.
The POSIX alarm boundary used here requires a POSIX Python environment.
`run --resume` validates signatures and completed receipts and skips those
units. An attempted but uncommitted unit stops automatic continuation.
The relocated replay runner was inspected and parsed; no fresh numerical
replay was performed for this publication. Read-only record validation and
saved-output reproduction were performed.

The code uses Python's standard library and mpmath. Accepted calculations
used Python 3.13.5, mpmath 1.3.0, and R 4.6.0 for figures. No package wrapper
is required. The repository's GPL-3.0-or-later license applies to this code.
