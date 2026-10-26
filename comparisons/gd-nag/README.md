# Synthetic GD and NAG-Q comparisons

Run from the repository root:

```sh
cd periodic-momentum
python -B comparisons/gd-nag/run_pilot.py verify
python -B comparisons/gd-nag/report.py --out outputs/synthetic_comparison_01
```

These commands check/export saved records only. Results are in
[reference/results.md](reference/results.md). Synthetic CHB/HB objective
values are available at saved state points; absent points are not interpolated.

An optional numerical replay (not performed after relocation) uses
`python -B comparisons/gd-nag/run_pilot.py run`, writes outputs/run_001,
and can continue completed, signature-matching units with `--resume`.
Interrupted units require inspection and are not automatically retried.

The article figure is rendered from reference/trajectories.csv by
render_article_figure.R, using base R/Cairo and a new output directory.
