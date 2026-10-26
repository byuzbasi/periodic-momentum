# Istanbul Stock Exchange application

All 536 rows, seven raw international returns and the USD-based ISE response
from UCI are retained. No intercept, scaling, centering, penalty, exclusions,
train/test split or parameter retuning is used. The fixed quadratic concerns
contemporaneous returns; it is not a forecasting evaluation.

The accepted records contain 16 full 512-query trajectories plus four two-step
prefix checks: 8,200 historical updates, two QR references and two eigensolves.
160 digits are primary; 120 digits check arithmetic. These are deterministic
arithmetic checks, not statistical replications. Data measurement precision
is distinct from working arithmetic precision.

From the repository root, with mpmath 1.3.0:

```sh
cd periodic-momentum
python -B applications/istanbul/verify_saved.py
python -B applications/istanbul/report.py --out outputs/istanbul_report_01
Rscript --vanilla applications/istanbul/plot.R applications/istanbul/outputs/istanbul_report_01 outputs/istanbul_figures_01
```

Reporting reads saved records and does not fit or optimize a model.
See [reference/RESULTS.md](reference/RESULTS.md), [protocol.md](protocol.md)
and [config.json](config.json). Reports are checked against accepted records.
The numerical record files and all feature/response tokens are unchanged; only the source-file metadata path is made portable.

CHB has smaller distance and objective gap than HB, GD and NAG-Q at each
recorded query 62 through 512 under both policies, after substantial early
overshoot: its objective gap peaks at about 62.815 times initialization.
Both adjacent terminal phases and their maxima are reported. A single
selected Hessian and tiny certified error budget do not establish broad
empirical superiority, financial prediction accuracy or practical noise tolerance.

Optional numerical replay (not executed after relocation):

```sh
cd periodic-momentum
python -B applications/istanbul/run_pilot.py design
python -B applications/istanbul/run_pilot.py run
```

Replay writes applications/istanbul/outputs/run_001 and refuses overwrites.
Use the same run command with --resume only for completed, signature-matching
units. Interrupted/failed units need inspection and are not automatically
retried. The replay keeps the original 8,200-update and 120-second limits.
Portability checks cover saved-output verification/reporting, static parsing
and mathematical function identity. No fresh relocated replay was performed.

The article figure is rendered from reference/trajectories.csv by
render_article_figure.R, using base R/Cairo and a new output directory.

The [candidate screening](selection/README.md) records all four eligible
datasets and the feature-only selection rule. The [reference-accuracy note](REFERENCE_ACCURACY.md)
defines the saved quadratic-energy measure and its residual correction.
