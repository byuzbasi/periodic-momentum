# QR reference and the recorded energy gap

Let w_star minimize f(w) = ||Xw-y||^2/(2n), and let w_hat denote the
finite-precision QR reference. With H = X^T X/n, h = X^T y/n and
r = H w_hat - h, the exact quadratic expansion is

    f(w) - f(w_hat) = E(w) + r^T (w-w_hat),
    E(w) = (w-w_hat)^T H (w-w_hat)/2.

The saved F_k is E(w_k)/E(w_0), and D_k is
||w_k-w_hat||/||w_0-w_hat||. Thus F is a QR-referenced energy measure,
computed without cancellation, rather than the exact objective difference
relative to a finite-precision reference. E equals the true objective gap
when w_hat = w_star. The numerical adapter's historical `star` variable name
is retained; no trajectory, metric computation or accepted result was changed.

For H >= mu I with mu > 0, w_hat - w_star = H^(-1) r gives the bounds

    ||w_hat - w_star|| <= ||r||/mu,
    |[f(w)-f(w_star)] - E(w)|
        <= ||r|| ||w-w_hat|| + ||r||^2/(2 mu).

The first follows from the norm of H^(-1). The second follows by expanding
the exact-minimizer energy: the correction is
r^T(w-w_hat) + r^T H^(-1)r/2, and applying Cauchy-Schwarz and H >= mu I.
These formulas bound reference error without computing another solution.
For normalized quantities, if a_k and a_0 bound relative numerator and
denominator errors and a_0 < 1, their ratio differs relatively by at most
(a_k+a_0)/(1-a_0). The diagnostic applies this to both D and F.

## Check from existing records

The diagnostic evaluates the saved H, h and w_hat at 220 decimal digits,
then recomputes their residual. It separately records the original residual
used in the affine-state checks, because decimal serialization and the
original arithmetic give slightly different tiny residuals. Both are
reported; the larger reconstruction-based values are used below.

Each precision has eight accepted trajectories. Their saved current/previous
states cover k = 0, 1, 2, 15, 16, 63, 64, 511 and 512, giving 72 states per
precision. No interpolation or new optimizer iteration is used.

| Largest value over saved states | 120 digits | 160 digits |
|---|---:|---:|---:|
| Absolute linear correction / E | 1.480e-83 | 4.943e-123 |
| Relative energy error bound against the exact minimizer | 7.525e-83 | 1.099e-122 |
| Relative distance error bound against the exact minimizer | 5.774e-81 | 8.428e-121 |

The normalized-ratio bounds have the same displayed values. These
diagnostics are far below the precision shown in the article tables.
They are floating-point checks on stored quantities, not interval-arithmetic
certificates of the original real-valued data or checks of every unsaved state.

Reproduce only this saved-output diagnostic from the repository root:

```sh
cd periodic-momentum
python -B applications/istanbul/check_reference_accuracy.py --out outputs/reference_accuracy_01.json
```

The command refuses an existing output. It uses mpmath 1.3.0, performs no
optimizer updates, QR solves or eigensolves, and does not modify any input.
The published result is [reference/reference_accuracy.json](reference/reference_accuracy.json).
