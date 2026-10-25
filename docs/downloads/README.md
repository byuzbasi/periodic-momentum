# A fourth-power relative-error threshold for two-periodic heavy ball

Companion code for the working manuscript by Bahadır Yüzbaşı and Mohammad
Arashi, prepared for *Optimization Letters*.

The numerical examples study established two-periodic cyclical heavy-ball
coefficients on quadratics with a gapped spectrum. They illustrate the scalar
upper bound's fourth-power scaling, finite-horizon transients, and alternating
phases on the Yacht Hydrodynamics least-squares problem. The uniform vector-error
guarantee is a theorem in the manuscript; these selected computations do not
establish an adversarial optimum or a general performance advantage.

This is a collection of scripts and frozen research inputs, not an installable
Python or R package.

[GitHub repository](https://github.com/byuzbasi/periodic-momentum) · [Research website](https://byuzbasi.github.io/periodic-momentum/) · [Citation](CITATION.cff)

## Research companion website

The [research website](https://byuzbasi.github.io/periodic-momentum/) presents the theorem,
the three existing figures, their interpretation and reproduction commands.
Its source is in `docs/index.html`. It can be previewed with a static web server
and requires no build or external
web dependencies. The download on the site contains the validated research
source and records. It does not start a calculation.

## Requirements

- Python 3.13.5 and `mpmath==1.3.0` were used for the accepted calculations and
  the limited checks described below. The scripts use the Python standard
  library and mpmath. The command-line interface requires Python 3.10 or later;
  other Python versions have not been validated here.
- Base R with Cairo support is needed only to redraw the figures. The validated
  rendering used R 4.6.0 and Arial. Different font installations may change
  typography; they do not change the numerical plotting data.
- CVXPY, SDP solvers, NumPy and SymPy are not dependencies of this workflow.

After placing this directory on your computer, an optional isolated environment
can be created with the following commands. An existing environment with the
pinned dependency can also be used.

Each command block starts from the parent of this directory. If you are already
inside it, omit the first `cd` command.

Commands use paths relative to the repository. No personal directory or machine
configuration is required. If your checkout has a different name, adjust only
the first `cd` command.

```bash
cd periodic-momentum
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Check the accepted records

```bash
cd periodic-momentum
.venv/bin/python scripts/reproduce.py check --out outputs/check_local
```

This verifies the frozen files, raw-data dimensions, 13 synthetic work units,
16 Yacht records, trace dimensions, and the original numerical acceptance
criteria. The Yacht precision comparison is recomputed from the saved 80/120-digit
logs. No optimizer iteration, QR solve or eigensolver is called. A completion
marker alone is not used as evidence of validity.

## Recreate the reporting data and figures

```bash
cd periodic-momentum
.venv/bin/python scripts/reproduce.py export --out outputs/report_local
.venv/bin/python scripts/reproduce.py figures --analysis outputs/report_local --out outputs/figures_local
```

The exporter recreates five reporting CSV files and the noise-effect and phase
diagnostics directly from accepted records. It compares them with `reference/`:
CSV files must match byte for byte, and JSON values must match exactly. These
comparisons passed in the validated environment. Small platform-dependent
floating-point formatting differences are not silently accepted or used to
change the frozen reference values.

`adjacent_queries.csv` additionally reports the two exact-gradient Yacht rows
at queries 511 and 512. `reporting_check.json` records the comparisons. The
existing directory names in the generated reports preserve the original plot
script's input contract.

The figure command reads these exported CSV files only. The resulting three
vector PDFs are in `outputs/figures_local/figures/`. Rendering never starts a
trajectory. In the limited verification, their text and drawing commands
matched all three manuscript figures, and their fonts were embedded.

All output commands require a new directory and refuse unintended overwrites.
Use a new output name for another reporting check.

## Optional eight-update prefix check

```bash
cd periodic-momentum
.venv/bin/python scripts/reproduce.py smoke --out outputs/smoke_local
```

This explicitly requested numerical check uses only existing exact-gradient
cases at 120 decimal digits:

| Existing case | Method | Work performed |
|---|---|---|
| Synthetic T01, R = 1/2 | CHB | Two independent state-transition updates |
| Synthetic T01, R = 1/2 | HB | Two independent state-transition updates |
| Yacht Y02 | CHB | Two gradient updates from the saved setup |
| Yacht Y04 | HB | Two gradient updates from the saved setup |

The synthetic transition implementation is compared against the saved direct
iteration states. The Yacht check exercises the unchanged gradient loop and
compares its first two checkpoints with accepted records. The original frozen
matrix-check tolerance is used. The Yacht Hessian and reference solution are
loaded from the accepted 120-digit setup; no new QR or spectrum computation is
performed. A two-step runtime prefix does not modify either 512-query
configuration file. It is not a new scientific scenario or a new full trajectory.

One such check was performed for this handoff: four state transitions and four
gradient updates, eight updates in total. It passed. See
`validation/limited_check.json`. Do not interpret this limited check as a fresh
recomputation of every 512-query result.

The complete verification scope and remaining limits are recorded in
[`validation/summary.json`](validation/summary.json).

## Frozen design

| Component | Fixed choices |
|---|---|
| Synthetic | mu = 1/16, L = 1, condition number 16, four endpoint modes, rest initialization |
| Scalar scaling | R = 1/2, 1/4, 1/8, 1/16, 1/32 |
| Synthetic trajectories | R = 1/2 and 1/8; exact, scalar-half, scalar-above and mixing-half policies |
| Yacht | All 308 rows, six raw inputs, seventh column as target |
| Objective | No-intercept, unpenalized least squares; no centering, scaling, exclusions or split |
| Yacht methods/policies | CHB and fixed noiseless HB; four original error policies |
| Precision/horizon | 80 and 120 decimal digits; 512 queries |
| Randomness | None; no RNG or seed is used |

The two mathematical kernels and the original synthetic and Yacht
configurations are retained byte for byte. The new adapters handle portable
file access, validation, reporting and explicitly requested execution.

The scalar upper obstruction is a threshold for retaining the benchmark rate,
not a divergence threshold. The fixed HB comparison is not the optimal noisy
stationary method. Yacht's even and odd phases and its large objective
transients are retained. Distance and objective metrics should not be conflated.

## Complete numerical replay: prepared, not run for this handoff

The following commands are provided for a separately authorized replay of the
same frozen design. They were **not executed** during the limited validation.
They perform fresh QR/spectral preparation and all existing trajectories;
there is no parameter search or additional dataset. The newly adapted complete
runner has only static/control-flow validation, not a fresh end-to-end replay.

```bash
cd periodic-momentum
.venv/bin/python scripts/reproduce.py recompute --out outputs/replay_local --confirm-frozen-replay
.venv/bin/python scripts/reproduce.py export --records outputs/replay_local --out outputs/replay_reports
```

The complete replay comprises two Yacht setup records, five scalar records,
eight synthetic paired records and sixteen Yacht records. It executes 48
trajectories and 24,576 gradient updates, plus two QR solves and the prescribed
spectral checks. These are the original work units, not new replications.
Runtime for this adapted full runner has not been measured. Storage is expected
to be of the same order as the small accepted archive, but a fresh measurement
has not been made.

Each accepted work unit is saved atomically with an input/software signature
and a receipt. Progress counts durably saved units (31 including setup), not
elapsed time. Inspect it with:

```bash
cd periodic-momentum
cat outputs/replay_local/progress.json
```

A compatible interrupted run can resume completed, verified work:

```bash
cd periodic-momentum
.venv/bin/python scripts/reproduce.py recompute --out outputs/replay_local --confirm-frozen-replay --resume
```

An attempted but uncommitted unit stops automatic continuation for inspection.
Changed signatures and invalid receipts are rejected. After a complete replay,
the exporter validates its records and compares the reporting values with the
accepted reference. Neither receipt counts nor `progress.json` replace the
scientific validation.

## Directory guide

| Path | Contents |
|---|---|
| `code/` | Unchanged numerical kernels, retaining their relative import layout |
| `configs/` | Original scientific configurations, spectral bounds and validation budget |
| `data/raw/` | Unchanged Yacht numeric file |
| `data/accepted/` | Accepted scalar, trajectory and setup records |
| `reference/` | Original reporting values and manuscript figures |
| `scripts/` | Portable command-line and reporting adapters |
| `validation/` | Results and scope of the limited verification |
| `docs/` | Static research website, figure assets and source download |
| `outputs/` | Local generated results, excluded by `.gitignore` |

## Data attribution, code licence and manuscript status

Gerritsma, J., Onnink, R., and Versluis, A. (1981), *Yacht Hydrodynamics*,
UCI Machine Learning Repository, [DOI: 10.24432/C5XG7R](https://doi.org/10.24432/C5XG7R).
The numeric file is distributed unchanged under CC BY 4.0; see
`data/ATTRIBUTION.txt`. The data licence does not assign a licence to the code.
The source code in `code/` and `scripts/`, together with the website's HTML,
CSS and JavaScript source, is distributed under the GNU General
Public License, version 3 or (at your option) any later version. See
[`LICENSE.md`](LICENSE.md). This code licence does not change the Yacht data
licence or assign a software licence to the manuscript figures.

The manuscript remains a working draft; no journal acceptance or published
paper DOI is claimed. The code is publicly available in this repository. The present checks cover the numerical
illustrations and reporting; they are not a machine-checked proof of the theorem
and do not replay the earlier SDP exploration.
