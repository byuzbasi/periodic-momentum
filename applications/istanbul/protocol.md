# Istanbul Stock Exchange: fixed application protocol

## Contribution and selection

The fixed-Hessian real-data example tests finite-query CHB behavior under
the existing spectral-gap theorem's assumptions. The established CHB method
and its coefficients are not presented as a new algorithm. A single Hessian
cannot verify the fourth-power scaling law, universal empirical superiority,
financial prediction accuracy or wall-clock speedups.

Istanbul was selected by feature-only spectral geometry before its optimizer outcomes.
The [four-candidate screening records and selection rule](selection/README.md)
make the targeted selection explicit.
All early losses and overshoots are reported. This is a selected illustration,
not an unbiased population benchmark.

## Exact input/model definitions

Source: UCI data_akbilgic.xlsx, original worksheet `orjinal data`, all 536 rows
3–538 in source order. Seven inputs are SP, DAX, FTSE, NIKKEI, BOVESPA, EU, EM
(zero-based columns 3–9). The response is the USD BASED ISE series in C
(zero-based column 2), identified by the original first-row currency label.
Date and the TL-based ISE column are not inputs. input_tokens.json retains
every source numeric token before mpmath parsing; raw source bytes stay in
the repository data/raw directory. Source hashes and column identity are checked.

For X in R^(536x7), y in R^536 and w in R^7, the objective is
f(w)=||Xw-y||²/(2*536), H=X^T X/536, h=X^T y/536, and g(w)=Hw-h.
This fixed quadratic illustrates optimization behavior on contemporaneous
returns. No intercept, normalization, centering, regularization, observation
selection, lagged variable, imputation or train/test split is introduced.
No out-of-sample or causal claim is made.

At each working precision, original decimal tokens are read directly, H and
h are formed once, and one QR least-squares reference w_hat and one symmetric
eigensolve of H are computed. Numerical precision does not change the
measurement precision of the supplied returns. References and matrices are
stored at their full working precision; serialization retains the complete working precision.

The frozen bounds are mu=5.235251356211395e-6, L=0.0008564028670237406 and
R=0.49644090689845105, taken from the recorded safe enclosure. All seven
high-precision eigenvalues must lie within these endpoints and the theorem's
equal-width two bands [mu,a] and [b,L]. Here a=(L+mu-R*(L-mu))/2 and
b=(L+mu+R*(L-mu))/2. Failed enclosure stops the pilot; no retuning occurs.
Each method starts at w[-1]=w[0]=0; CHB uses alpha0 on its first query.

## Four methods, shared policies and independent state check

CHB/HB coefficients and the norm/oracle helpers reuse the unchanged shared numerical kernel.
The GD/NAG-Q coefficient function body is the same as the accepted GD/NAG
adapter. HB uses alpha=4/(sqrt(L)+sqrt(mu))², beta=q², eta=0, with
q=(sqrt(L/mu)-1)/(sqrt(L/mu)+1). GD uses alpha=2/(L+mu), beta=eta=0.
NAG-Q uses alpha=4/(3L+mu), beta=eta=(sqrt(3L/mu+1)-2)/(sqrt(3L/mu+1)+2).
CHB uses the manuscript's alternating alpha and momentum m, eta=0. Each
update is w_next=w+beta*(w-w_previous)-alpha*(g(z)+e), where
z=w+eta*(w-w_previous). The reported iterate is w, not NAG's query point z.

The two policies are exact (e=0) and certified_mixing. The latter uses
epsilon=c_kappa*R^4/2, evaluated from the same theorem formulas at each
precision, and e_k=epsilon*Q_k*g(z_k). For P(v1,...,v7)=(v7,v1,...,v6),
Q_k cycles P,P^T,-P beginning at k=0. The same policy and relative budget
apply to each method's own gradient; realized error vectors can differ.
This is controlled computational error, not measurement noise.

Independent checks use a centered affine state [w-w_hat;w_previous-w_hat;1].
They retain r=H*w_hat-h, the finite-precision reference residual. Writing
S=I+epsilon*Q_k, the top blocks are (1+beta)I-alpha*(1+eta)*S*H,
-beta*I+alpha*eta*S*H, and -alpha*S*r. The bottom blocks shift the
previous iterate and retain the constant 1. Their product over a period
of two for exact gradients and six for mixing gives the independent state
at k=1,2,16,64,512. Matrix powers are validation work, not gradient queries;
their time is included in the numerical budget.

## Frozen workload and acceptance

Run order is profile, method, working precision, in config.json order:
two policies x four methods x two precisions = 16 full 512-query trajectories,
8,192 updates. Before these, run each method for two exact-gradient updates
at 160 digits: four prefixes and eight updates. The exact total limit is
8,200 updates; prefixes cannot be repeated. Preparation is limited to two
QR references and two eigensolves. There is no RNG or statistical repetition.

The QR normalized gradient residual must be <=1e-60. Every full trace must
contain k=0,...,512; every independent state error is relative to the actual
two-iterate state and must be <=1e-60. Positive relative oracle violation
must be <=1e-60. Across 120/160 digits, all distance/objective/state log
differences must be <=1e-45. Exact GD's objective log must be nonincreasing,
allowing 1e-60 rounding slack. A failed gate stops the pilot without relaxing
thresholds or automatically retrying.

Persistent reservations precede reference calls and trajectory attempts.
The actual gradient counter is separate from the reserved budget. Atomic
JSON writes, receipts, a single-process file lock and a frozen scientific/
software signature support continuation only after completed, valid units.
Interrupted or failed units require inspection and are not automatically
rerun. Raw data and accepted outputs are read-only inputs to replay.

One local process has a cumulative 120-second numerical ceiling including
reference preparation, prefixes, trajectories and independent state checks.
This is a stop limit. The completed original run took about 6.27 seconds;
runtime on other environments is unknown. Publication checks do not execute
a new numerical replay.

## Reporting

Primary results are the 160-digit records, with 120 digits providing an
arithmetic check. D_k is the Euclidean distance to the QR reference divided
by its initial value. The exact minimizer w_star and the finite-precision QR reference w_hat are
distinct. The saved F_k is E(w_k)/E(w_0), where
E(w) = (w-w_hat)^T H(w-w_hat)/2. With r = H w_hat - h,
f(w)-f(w_hat) = E(w) + r^T(w-w_hat); the linear residual term is not
part of the stored energy. The adapter's historical variable name `star`
continues to denote w_hat. See [reference accuracy](REFERENCE_ACCURACY.md)
for checks from saved states and bounds relative to the exact minimizer. Both denominators
are positive. The normalized two-iterate norm is also retained.

Report both adjacent terminal queries, all-query curves, peaks and maxima
within odd/even pairs. Retain early losses. A single selected Hessian does
not establish broad superiority or the fourth-power scaling law by itself.
