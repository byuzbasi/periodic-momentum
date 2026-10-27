# Manuscript V32

**A Fourth-Power Relative-Error Threshold for Two-Periodic Heavy Ball**

Bahadır Yüzbaşı and Mohammad Arashi

Working manuscript prepared for *Optimization Letters*, 27 October 2026.

[Read the LaTeX source](main.tex) or download the
[self-contained source archive](manuscript_source_v32.zip).
The archive contains the publisher class, supporting styles and the three
unchanged vector figure assets. It contains no compiled manuscript PDF.

After extracting the archive, compile with an existing TeX installation:

```sh
cd manuscript_source_v32
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

## Analytic refinements

The fourth-power threshold retains the benchmark rate q for fixed CHB
coefficients. It does not describe the full rate-robustness tradeoff, the
optimal noise-dependent tuning or an absolute stability margin.

For 0 < R <= 1/2, separate phase bounds use
A_0 = 8 kappa/(kappa+3), A_1 = 4 kappa/(kappa+1),
B_0 = 6(kappa-1)/(kappa+3), and B_1 = 6(kappa-1)/(3 kappa+1).
Here A_r bounds alpha_r L, B_r bounds |1+m-alpha_r lambda| and
T_r = sqrt(B_r^2+2) bounds the phase-r state transition.
The improved energy and perturbation bounds are
K = sqrt((T_0 T_1+1)^2+2) and
G = T_1 A_0+T_0 A_1+A_0 A_1.
They give c_kappa = min(1,d_minus^2/(4 K G)), with
d_minus = (3/2) q^2 sqrt(kappa)/(kappa+1).
The odd-iterate factor is updated to T_0+A_0 epsilon.

Using m >= 3q^2/4 also improves the coarse upper coefficient to
C_kappa = 4q^2/[(kappa-1)(1+3q^2/4)^2]. The exact scalar upper
obstruction and the fourth-power order are unchanged. No optimal
vector-error coefficient is claimed.

The comparison with stationary methods is Corollary 1. Its proof uses
r_S*(epsilon) >= r_S*(0) = q because zero error is admissible.
It does not establish strict deterioration of the best stationary rate
at a positive error budget.

## Unchanged numerical study

The fixed configurations, code interfaces, data, accepted trajectories,
tables and figures are unchanged. The Istanbul budget remains
epsilon = c_kappa^p R^4/2, approximately 1.60364486e-15, where the
historical protocol constant c_kappa^p is explicitly defined in the
article and [application protocol](../../applications/istanbul/protocol.md).
The improved c_kappa is at least this earlier constant; no experiment
was retuned. The low-dimensional examples illustrate the theory and
transients rather than scalability or broad empirical superiority.

The article now defines the spectral enclosure and distinguishes rounded
display values from the full stored tuning values. Its numerical findings
do not include a new rounding-sensitivity experiment.

## Verification

The source archive compiled independently to the same 16-page article.
All three tables, all three figures, the references and the generative-AI
statement are preserved. Checks covered the phase bounds, noncommuting
two-step error product, repeated-root energy identity, rest initialization
and both iterate parities. Twenty exact symbolic/sign checks passed.
No optimizer updates, QR solves or eigensolves were performed for V32.
The symbolic audit is an internal verification, not an additional numerical
benchmark or a claim that the full theorem was machine-formalized.
