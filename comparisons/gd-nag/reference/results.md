# GD/NAG bounded comparison: saved results

NAG-Q is the constant-parameter Nesterov tuning for strongly convex quadratics. All trajectories use 512 gradient queries. Lower D and F are better.

## Distance and objective at the final adjacent pair

| case | profile | method | D511 | D512 | F511 | F512 | pair_max_D | pair_max_F |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic_R_half | exact | GD | 1.182277855928544e-28 | 1.043186343466363e-28 | 1.397780928618996e-56 | 1.08823774719472e-56 | 1.182277855928544e-28 | 1.397780928618996e-56 |
| synthetic_R_half | exact | NAG-Q | 1.566261918409035e-73 | 1.12093293807525e-73 | 2.88674403858601e-147 | 1.47855805101879e-147 | 1.566261918409035e-73 | 2.88674403858601e-147 |
| synthetic_R_half | scalar_half | GD | 1.587781898717029e-28 | 1.401956724613397e-28 | 3.585544798679268e-56 | 2.796754571630554e-56 | 1.587781898717029e-28 | 3.585544798679268e-56 |
| synthetic_R_half | scalar_half | NAG-Q | 5.508903316295368e-75 | 4.143506882030191e-75 | 6.418066381994604e-150 | 3.476942909867063e-150 | 5.508903316295368e-75 | 6.418066381994604e-150 |
| synthetic_R_eighth | exact | GD | 1.182277855928544e-28 | 1.043186343466363e-28 | 1.397780928618996e-56 | 1.08823774719472e-56 | 1.182277855928544e-28 | 1.397780928618996e-56 |
| synthetic_R_eighth | exact | NAG-Q | 1.566261918409035e-73 | 1.12093293807525e-73 | 2.88674403858601e-147 | 1.47855805101879e-147 | 1.566261918409035e-73 | 2.88674403858601e-147 |
| synthetic_R_eighth | scalar_half | GD | 1.183242054358311e-28 | 1.044038774109431e-28 | 1.402342062378774e-56 | 1.091795764272298e-56 | 1.183242054358311e-28 | 1.402342062378774e-56 |
| synthetic_R_eighth | scalar_half | NAG-Q | 1.556436656939074e-73 | 1.113874019198915e-73 | 2.850651774967306e-147 | 1.460000597494756e-147 | 1.556436656939074e-73 | 2.850651774967306e-147 |
| synthetic_R_half | exact | CHB | 5.653420562436065e-127 | 1.094463731650922e-127 | 5.634290975611547e-253 | 1.918114033156575e-254 | 5.653420562436065e-127 | 5.634290975611547e-253 |
| synthetic_R_half | exact | HB | 1.822173658150321e-111 | 1.095439117563931e-111 | 5.902952666552966e-222 | 2.133372912840261e-222 | 1.822173658150321e-111 | 5.902952666552966e-222 |
| synthetic_R_half | scalar_half | CHB | 1.671809021416789e-117 | 3.53377133412548e-118 | 5.261073701817338e-234 | 2.350595734943319e-235 | 1.671809021416789e-117 | 5.261073701817338e-234 |
| synthetic_R_half | scalar_half | HB | 2.14773849086117e-103 | 1.346365161247107e-103 | 8.682881176708588e-206 | 3.412139571614072e-206 | 2.14773849086117e-103 | 8.682881176708588e-206 |
| synthetic_R_eighth | exact | CHB | 3.183668177582331e-112 | 1.515227812655893e-112 | 1.818075499291535e-223 | 4.0339603353333e-224 | 3.183668177582331e-112 | 1.818075499291535e-223 |
| synthetic_R_eighth | exact | HB | 1.822172753337059e-111 | 1.095439114570118e-111 | 5.902944863084136e-222 | 2.133372858332647e-222 | 1.822172753337059e-111 | 5.902944863084136e-222 |
| synthetic_R_eighth | scalar_half | CHB | 4.173117037302799e-112 | 1.982904882129787e-112 | 3.191741232821148e-223 | 7.122923110504624e-224 | 4.173117037302799e-112 | 3.191741232821148e-223 |
| synthetic_R_eighth | scalar_half | HB | 2.380478941953006e-111 | 1.432517430686921e-111 | 1.033240522510966e-221 | 3.741995244643128e-222 | 2.380478941953006e-111 | 1.033240522510966e-221 |

## CHB / comparator ratios

| case | profile | comparator | CHB_to_other_D511 | CHB_to_other_D512 | CHB_to_other_pair_max_D | CHB_to_other_F512 |
| --- | --- | --- | --- | --- | --- | --- |
| synthetic_R_half | exact | HB | 3.102569580648435e-16 | 9.991095936804068e-17 | 3.102569580648435e-16 | 8.990992721487679e-33 |
| synthetic_R_half | exact | GD | 4.781803646314553e-99 | 1.049154581543094e-99 | 4.781803646314553e-99 | 1.762587300524289e-198 |
| synthetic_R_half | exact | NAG-Q | 3.609498830296949e-54 | 9.763864496035078e-55 | 3.609498830296949e-54 | 1.29728692886621e-107 |
| synthetic_R_half | scalar_half | HB | 7.784043674453358e-15 | 2.624675263323235e-15 | 7.784043674453358e-15 | 6.888920237900461e-30 |
| synthetic_R_half | scalar_half | GD | 1.052921073585519e-89 | 2.520599439401348e-90 | 1.052921073585519e-89 | 8.404726531198204e-180 |
| synthetic_R_half | scalar_half | NAG-Q | 3.034740175002833e-43 | 8.528455327179377e-44 | 3.034740175002833e-43 | 6.760524391334315e-86 |
| synthetic_R_eighth | exact | HB | 0.1747182407239862 | 0.1383214997987828 | 0.1747182407239862 | 0.01890883874132565 |
| synthetic_R_eighth | exact | GD | 2.692825685280152e-84 | 1.452499663311353e-84 | 2.692825685280152e-84 | 3.706874114348743e-168 |
| synthetic_R_eighth | exact | NAG-Q | 2.032653759989398e-39 | 1.351755989308053e-39 | 2.032653759989398e-39 | 2.728307037085035e-77 |
| synthetic_R_eighth | scalar_half | HB | 0.1753057741346397 | 0.1384209950715185 | 0.1753057741346397 | 0.01903509396678545 |
| synthetic_R_eighth | scalar_half | GD | 3.526849829188957e-84 | 1.899263639725653e-84 | 3.526849829188957e-84 | 6.524043546965199e-168 |
| synthetic_R_eighth | scalar_half | NAG-Q | 2.681199404227443e-39 | 1.780187748301974e-39 | 2.681199404227443e-39 | 4.8787124626709e-77 |

## Odd/even distance counts (256 per phase)

| case | profile | comparator | odd_CHB_wins | odd_CHB_losses | odd_ties | even_CHB_wins | even_CHB_losses | even_ties |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic_R_half | exact | HB | 250 | 6 | 0 | 256 | 0 | 0 |
| synthetic_R_half | exact | GD | 253 | 3 | 0 | 255 | 1 | 0 |
| synthetic_R_half | exact | NAG-Q | 251 | 5 | 0 | 254 | 2 | 0 |
| synthetic_R_half | scalar_half | HB | 250 | 6 | 0 | 256 | 0 | 0 |
| synthetic_R_half | scalar_half | GD | 253 | 3 | 0 | 255 | 1 | 0 |
| synthetic_R_half | scalar_half | NAG-Q | 251 | 5 | 0 | 254 | 2 | 0 |
| synthetic_R_eighth | exact | HB | 235 | 21 | 0 | 256 | 0 | 0 |
| synthetic_R_eighth | exact | GD | 253 | 3 | 0 | 254 | 2 | 0 |
| synthetic_R_eighth | exact | NAG-Q | 251 | 5 | 0 | 253 | 3 | 0 |
| synthetic_R_eighth | scalar_half | HB | 235 | 21 | 0 | 256 | 0 | 0 |
| synthetic_R_eighth | scalar_half | GD | 253 | 3 | 0 | 254 | 2 | 0 |
| synthetic_R_eighth | scalar_half | NAG-Q | 251 | 5 | 0 | 253 | 3 | 0 |

## Validation and scope

The published synthetic subset contains 16 accepted full trajectories and two prefixes: 8,192 + 4 = 8,196 historical updates. No new CHB/HB trajectory, QR solve or eigensolve.
All new trajectories passed independent affine-state and relative-error checks, and 80/120-digit log comparisons. See numerical_quality.json.
Synthetic CHB/HB objective values exist only at saved state points. Their full objective peaks are intentionally missing.
Synthetic baseline norm logs were stored as float64: comparisons use a 1e-12 log tie tolerance. All reported distances refer to the same reference solution and initialization.
No claim about arbitrary gradient-error paths, parameter-optimal noisy methods, out-of-sample prediction or general empirical superiority follows from this pilot.
