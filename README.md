# A fourth-power relative-error threshold for two-periodic heavy ball

Research materials for the working manuscript by **Bahadır Yüzbaşı and
Mohammad Arashi**, prepared for *Optimization Letters*.

[Research companion](https://byuzbasi.github.io/periodic-momentum/) ·
[Istanbul results](applications/istanbul/reference/RESULTS.md) ·
[Synthetic comparisons](comparisons/gd-nag/reference/results.md)

![Fourth-power threshold and two spectral bands](docs/assets/mark.svg)

The contribution is a sharp fourth-power order for the relative-gradient-error
threshold as a spectral gap closes at fixed condition number. The cyclical
heavy-ball coefficients are established in the literature; this repository
does not claim a new optimization algorithm. The paper proves a uniform
vector-error lower bound and a scalar obstruction of matching order.

## Numerical illustrations

- Controlled four-dimensional two-band quadratics illustrate the scalar
  upper obstruction, transients and comparisons with stationary heavy ball
  (HB), gradient descent (GD) and quadratically tuned Nesterov (NAG-Q).
- The **Istanbul Stock Exchange** application uses all 536 observations,
  seven raw international returns and the USD-based ISE response. It is
  no-intercept, unpenalized least squares without scaling or exclusions.
  CHB has smaller distance and objective gap than all three comparators
  at queries 62–512 under exact and certified mixing policies. Its early
  objective gap peaks at about 62.815 times initialization.

Both phases, adjacent-pair maxima, early losses and numerical checks are
reported. One selected Hessian and a tiny certified error budget do not
establish broad superiority, practical noise tolerance, forecasting accuracy
or wall-clock speedups. Arithmetic precision is not measurement precision.

## Contents

| Path | Purpose |
|---|---|
| `code/periodic_momentum_numerics_v06/code/numerics.py` | Original synthetic CHB/HB kernel, unchanged |
| `configs/synthetic.json` | Fixed synthetic protocol, no random draws |
| `data/accepted/synthetic` | Five scalar units and eight paired CHB/HB units |
| `comparisons/gd-nag` | Both synthetic GD/NAG-Q cases, 16 full records and two prefixes |
| `applications/istanbul` | Portable seven-dimensional four-method adapter, fixed protocol, input tokens and accepted records |
| `data/raw/data_akbilgic.xlsx` | Unmodified UCI source workbook |
| `reference/synthetic` | Original synthetic reporting CSVs |
| `docs` | Research website and source download |

The current distribution follows the manuscript's synthetic and Istanbul
scope. Previous distributions remain recoverable in Git history. The
synthetic comparator records form an explicitly documented subset of their
historical study, selected by case identity, with every retained record
unchanged. No outcomes, coefficients or trajectories were recomputed to
prepare this publication.

## Environment and saved-output reproduction

Validated with Python 3.13.5 and **mpmath 1.3.0**; base R 4.6.0 with Cairo
is used only for plotting. Python's ZIP/XML reader verifies the workbook
without pandas or a spreadsheet package. This is research code, not an
installable Python package.

```sh
cd periodic-momentum
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B scripts/reproduce.py check --out outputs/check_01
.venv/bin/python -B scripts/reproduce.py export --out outputs/synthetic_01
.venv/bin/python -B scripts/reproduce.py figures --analysis outputs/synthetic_01 --out outputs/synthetic_figures_01
.venv/bin/python -B comparisons/gd-nag/run_pilot.py verify
.venv/bin/python -B comparisons/gd-nag/report.py --out outputs/comparisons_01
.venv/bin/python -B applications/istanbul/verify_saved.py
.venv/bin/python -B applications/istanbul/report.py --out outputs/istanbul_01
Rscript --vanilla applications/istanbul/plot.R applications/istanbul/outputs/istanbul_01 outputs/istanbul_figures_01
```

Use new output names for repeated reporting. The `check`, `verify` and
reporting commands read saved results; they do not perform optimizer updates,
QR solves or eigensolves. The manuscript's 120/160-digit Istanbul records
contain 16 full trajectories plus four prefixes (8,200 historical updates).
Synthetic records use 80/120 digits. Precision pairs are arithmetic checks,
not independent statistical replications. Synthetic CHB/HB objective gaps
are reported only where saved states permit their evaluation.

## Article figures from saved outputs

The original CHB/HB diagnostic rendering remains available through the root
`figures` command. The current article's four-method figures are reproduced by:

```sh
cd periodic-momentum
Rscript --vanilla comparisons/gd-nag/render_article_figure.R comparisons/gd-nag/reference outputs/article_synthetic_01
Rscript --vanilla applications/istanbul/render_article_figure.R applications/istanbul/reference outputs/article_istanbul_01
```

Both renderers read accepted CSV reports only and refuse existing PDF outputs.

## Optional numerical replay

Saved-output checks are sufficient to inspect the reported results. Explicit
numerical commands are separate:

```sh
cd periodic-momentum
.venv/bin/python -B scripts/reproduce.py recompute --confirm-frozen-replay --out outputs/synthetic_replay_01
.venv/bin/python -B comparisons/gd-nag/run_pilot.py run
.venv/bin/python -B applications/istanbul/run_pilot.py run
```

These commands were **not run after relocation**. They retain fixed designs,
write fresh outputs and support `--resume` only for complete, verified units
with matching signatures. Interrupted/failed units require inspection;
automatic retry is disabled. See each module's README and protocol. Do not
confuse accepted historical records with a newly validated relocated replay.

## Data, license and citation

Akbilgic, O. (2013). *Istanbul Stock Exchange*. UCI Machine Learning
Repository. [DOI: 10.24432/C54P4J](https://doi.org/10.24432/C54P4J),
CC BY 4.0. See [data attribution](data/ATTRIBUTION.txt). The source workbook
is unchanged and its selected numeric cells are preserved as text tokens.
Code is GPL version 3 or later; see [LICENSE](LICENSE.md).

Use [CITATION.cff](CITATION.cff) for both authors and the working manuscript.
No publication DOI is claimed. Both authors declare no competing interests;
no funding was received for the study.

Questions or reproducibility issues can be reported through GitHub Issues.
Include the command, software versions and the relevant unit identifier;
do not upload private machine paths or credentials.
