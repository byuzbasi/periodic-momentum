# Fixed four-method comparison

The question is how finite-horizon GD and Nesterov behavior compares with
the existing CHB/HB illustrations. The manuscript's contribution remains
the fourth-power relative-error threshold for fixed CHB coefficients.

## Problems and parameters

- Synthetic: mu = 1/16, L = 1, H = diag(mu, a, b, L), R = 1/2 or 1/8,
  optimum zero, and rest initialization x[-1] = x[0] = (1,1,1,1).
- Yacht: all 308 observations and six raw features, no-intercept,
  unpenalized least squares. There is no centering, standardization,
  exclusion or train/test split. Starts, Hessian, reference solution and
  safe bounds are reused from the accepted 80/120-digit setup files.
- New methods: GD with alpha = 2/(L+mu), and NAG-Q with
  alpha = 4/(3L+mu), beta = (sqrt(3 kappa+1)-2)/(sqrt(3 kappa+1)+2).
  Here kappa = L/mu. Coefficients are fixed and are not tuned under noise.
- The update is y = x + beta (x-x_previous), g = H y - h,
  x_next = y - alpha (g+e). Each update uses one gradient query.
- Profiles: exact (e = 0), and scalar_half (e = epsilon g), where
  epsilon is half the CHB scalar upper obstruction for that problem.
  This is not asserted to be GD/Nesterov's own critical threshold.

CHB/HB are read from accepted records. No new baseline trajectory, QR
solution or eigendecomposition is needed for these comparisons.

## Fixed budget and ordering

The order is configuration case, profile, method, precision. Three problems,
two profiles, two methods and two precisions give 24 full trajectories.
Each uses 512 updates. Four two-update exact-gradient prefixes cover GD and
NAG-Q on synthetic R=1/2 and Yacht at 120 digits. The total is exactly
12,288 + 8 = 12,296 updates. There is no RNG or seed, parameter search,
horizon extension, automatic retry or additional dataset.

120 decimal digits supply primary results; 80 digits check arithmetic.
This precision repetition is not a statistical replication.

## Acceptance and persistence

Store normalized distance D, objective gap F, and two-step state-norm
logarithms at k=0..512. Positive objective gaps are evaluated by the
Hessian quadratic form. Gradients remain H y-h, using the original
reference; the independent affine-state matrix retains its fixed residual.
Checkpoints are k=1,2,16,64,512.

- Maximum 80/120 log difference: at most 1e-45.
- Independent affine-state relative error: at most 1e-60.
- Positive relative-error-budget violation: at most 1e-60.
- Exact GD objective gaps must not increase beyond 1e-60 log roundoff.

One writer reserves each unit's query budget before execution, commits
accepted records atomically and maintains a receipt ledger. Attempted but
uncommitted units are not retried silently. Resume requires matching
software/input signatures and valid committed receipts. The cumulative
numerical-stage POSIX alarm budget is 120 seconds, not a runtime estimate.

## Reporting

Compare distances throughout the horizon, both queries 511 and 512, and
the larger error in each adjacent odd/even pair. Report both favorable
and unfavorable CHB comparisons. Do not extend the horizon after seeing
the results.

Synthetic CHB/HB objective values are available only at stored states:
0,1,2,15,16,63,64,511,512. Missing values are neither interpolated nor
recomputed. Yacht has complete objective traces for all four methods.
Synthetic baseline norm logs are stored as float64; their log tie tolerance
is 1e-12. The Yacht tie tolerance is 1e-30. These are storage/roundoff
tolerances, not scientific success criteria.

All 24 full trajectories and four prefixes passed in the recorded study.
The interpretation is to report the expanded comparison and distinguish
asymptotic guarantees from finite-horizon distance and objective behavior.
Further numerical work is not implied by arithmetic acceptance.
