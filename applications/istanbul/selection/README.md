# Candidate selection for the Istanbul illustration

This table and the accompanying spectral records reproduce the completed,
feature-only screening preceding the Istanbul optimizer runs. They are not
a new experiment or a comparative optimization benchmark. The four candidates
were a targeted search for a Hessian satisfying the theorem; they are not a
random sample of applications. Earlier distributions remain in Git history.

| Dataset | n | p | Raw condition number | Maximum central gap | Used gap | Lower/upper band counts |
|---|---:|---:|---:|---:|---:|---:|
| Istanbul Stock Exchange | 536 | 7 | 163.584 | 0.501455 | 0.496441 | 6 / 1 |
| Combined Cycle Power Plant | 9568 | 4 | 113,096.356 | 0.999461 | 0.500000 | 3 / 1 |
| Concrete Slump Test | 103 | 7 | 211,234.538 | 0.983343 | 0.500000 | 6 / 1 |
| Computer Hardware | 209 | 6 | 12,524,622.586 | 0.959091 | 0.500000 | 5 / 1 |

All four candidates passed the recorded full-rank, positive-definiteness,
SVD/Gram agreement and central-band enclosure checks. Selection ranked eligible
candidates by the smallest raw condition number, then the largest used gap,
then the smallest feature and observation counts. Optimizer outcomes did not
enter this ranking. All features retained their original scale; no intercept,
centering, standardization, penalty or row removal was applied. Response and
identifier columns were excluded by their recorded roles. For Power, only
Sheet1 was used, without stacking the supplied shuffled copies.

The floating-point guard is tau = 64 eps lambda_max. Acceptance requires full
rank at max(n,p) eps s_max, lambda_min > tau, SVD/Gram eigenvalue differences
at most tau, and relative symmetry/eigenpair residuals at most tau/lambda_max.
The enclosing endpoints are lambda_min - tau and lambda_max + tau. The used
gap is min(1/2, R_kappa, 0.99 R_max_safe), with R_kappa defined in the theorem.
Both bands must be occupied. These are floating-point diagnostics, not an
interval-arithmetic certificate.

The unchanged JSON records contain the Gram matrices, all eigenvalues, SVD
cross-checks, feature means, enclosure endpoints and acceptance diagnostics.
The CSV retains the original row order and full stored numerical values;
the display above is rounded and sorted by the selection rule. config.json
documents the exact source locations, feature columns and safeguards. The
source links and dataset DOIs identify the original UCI data; only the chosen
Istanbul workbook is distributed with the application. No additional source
datasets were downloaded or processed for this publication update.

A selected favorable geometry cannot establish general empirical superiority.
The Istanbul application reports its early overshoot as well as the later
advantage, and does not evaluate financial forecasting.
