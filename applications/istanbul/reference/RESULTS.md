# Istanbul: accepted four-method pilot

All four methods use the same 536 observations, seven raw inputs, USD-based ISE response, rest initialization and 512 gradient queries. Lower D and F are better. The 160-digit records are primary; 120 digits are an arithmetic check.
F is normalized quadratic energy about the computed QR reference. It approximates the exact objective gap; the residual correction and reference-error bounds are documented in ../REFERENCE_ACCURACY.md. The accepted metric values are unchanged.

## Adjacent terminal queries and peaks

| profile | method | D511 | D512 | F511 | F512 | pair_max_D | pair_max_F | peak_D | peak_F |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| exact | CHB | 6.220265314353432e-38 | 1.81623014783989e-38 | 2.907940972227531e-74 | 2.34298845519701e-75 | 6.220265314353432e-38 | 2.907940972227531e-74 | 2.922099971280495 | 62.81487844491687 |
| exact | HB | 5.250788425082081e-33 | 4.498024337421709e-33 | 2.033481127889818e-64 | 1.492228859515649e-64 | 5.250788425082081e-33 | 2.033481127889818e-64 | 1.621991062409944 | 17.90877247894399 |
| exact | GD | 0.00144177262152723 | 0.001424252407620241 | 3.061462442848426e-6 | 2.987509617539161e-6 | 0.00144177262152723 | 3.061462442848426e-6 | 1.0 | 1.0 |
| exact | NAG-Q | 3.337019161472048e-20 | 3.041869379431494e-20 | 5.178568925074923e-41 | 4.302944556296563e-41 | 3.337019161472048e-20 | 5.178568925074923e-41 | 1.0 | 1.0 |
| certified_mixing | CHB | 6.220265315071612e-38 | 1.816230148038767e-38 | 2.907940972904048e-74 | 2.34298845574404e-75 | 6.220265315071612e-38 | 2.907940972904048e-74 | 2.922099971280583 | 62.81487844492074 |
| certified_mixing | HB | 5.250788425529969e-33 | 4.498024337806908e-33 | 2.033481128245972e-64 | 1.49222885977804e-64 | 5.250788425529969e-33 | 2.033481128245972e-64 | 1.621991062409964 | 17.90877247894448 |
| certified_mixing | GD | 0.00144177262152737 | 0.001424252407620379 | 3.061462442851499e-6 | 2.987509617542177e-6 | 0.00144177262152737 | 3.061462442851499e-6 | 1.0 | 1.0 |
| certified_mixing | NAG-Q | 3.337019161471969e-20 | 3.041869379431422e-20 | 5.178568925074718e-41 | 4.302944556296392e-41 | 3.337019161471969e-20 | 5.178568925074718e-41 | 1.0 | 1.0 |

## CHB / comparator ratios

| profile | comparator | CHB_to_other_D512 | CHB_to_other_F512 | CHB_to_other_pair_max_D | CHB_to_other_pair_max_F |
| --- | --- | --- | --- | --- | --- |
| exact | HB | 4.037839752732336e-6 | 1.570126753852959e-11 | 1.184634536908844e-5 | 1.430030961361887e-10 |
| exact | GD | 1.275216484186675e-35 | 7.842613933162668e-70 | 4.314317820631436e-35 | 9.498535508807169e-69 |
| exact | NAG-Q | 5.970769685644201e-19 | 5.445081674985749e-35 | 1.864018458800071e-18 | 5.615337005841742e-34 |
| certified_mixing | HB | 4.037839752828689e-6 | 1.570126753943457e-11 | 1.184634536944572e-5 | 1.430030961444113e-10 |
| certified_mixing | GD | 1.275216484326187e-35 | 7.842613934985808e-70 | 4.31431782112914e-35 | 9.498535511007418e-69 |
| certified_mixing | NAG-Q | 5.970769686298141e-19 | 5.445081676257259e-35 | 1.864018459015331e-18 | 5.615337007148343e-34 |

## Counts across all 512 positive queries

| profile | comparator | D_all_CHB_wins | D_all_CHB_losses | D_all_ties | F_all_CHB_wins | F_all_CHB_losses | F_all_ties |
| --- | --- | --- | --- | --- | --- | --- | --- |
| exact | HB | 497 | 15 | 0 | 497 | 15 | 0 |
| exact | GD | 496 | 16 | 0 | 489 | 23 | 0 |
| exact | NAG-Q | 488 | 24 | 0 | 457 | 55 | 0 |
| certified_mixing | HB | 497 | 15 | 0 | 497 | 15 | 0 |
| certified_mixing | GD | 496 | 16 | 0 | 489 | 23 | 0 |
| certified_mixing | NAG-Q | 488 | 24 | 0 | 457 | 55 | 0 |

## Counts across 256 odd/even pairs

| profile | comparator | D_pair_CHB_wins | D_pair_CHB_losses | F_pair_CHB_wins | F_pair_CHB_losses | D_CHB_winning_suffix_start | F_CHB_winning_suffix_start |
| --- | --- | --- | --- | --- | --- | --- | --- |
| exact | HB | 241 | 15 | 241 | 15 | 30 | 30 |
| exact | GD | 245 | 11 | 242 | 14 | 22 | 28 |
| exact | NAG-Q | 240 | 16 | 225 | 31 | 32 | 62 |
| certified_mixing | HB | 241 | 15 | 241 | 15 | 30 | 30 |
| certified_mixing | GD | 245 | 11 | 242 | 14 | 22 | 28 |
| certified_mixing | NAG-Q | 240 | 16 | 225 | 31 | 32 | 62 |

## Scope

16 full trajectories plus four two-step prefixes: 8,200 total updates. Two QR references and two high-precision eigensolves were performed once. No retuning, additional dataset, longer horizon or automatic retry occurred.
Both metrics and both phases are reported, including peaks and early losses. This one fixed quadratic is an application illustration. It establishes neither broad empirical superiority, financial prediction accuracy, wall-clock speedups nor the fourth-power law by itself.
