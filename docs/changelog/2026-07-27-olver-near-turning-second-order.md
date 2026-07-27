# Olver map: second-order near-turning expansion, branch-cut placement, and cleanup

## Change

`fortran/hyperspherical_bessels_olver.f90`:

1. **Near-turning expansion carried to second order.** `compute_olver_z_amp` has a
   special branch for `|chi - chi_t| <= 1e-4*max(chi_t, z_t)`, because
   `analytic_amplitude` there is a ratio of two quantities that both vanish at the
   turning point. That branch previously used only the leading terms,
   `z = z_t + C_t^(1/3)*eps` with a *constant* `amp = C_t^(-1/6)`. It now carries the
   next order, at no cost in transcendental evaluations:

   ```
   c1 = C_t**(1/3)                                        ! C_t = sqrt(1 - K*z_t^2)
   c2 = (3*alpha^2*(c1^4 - C_t^2) - K) / (10*alpha*c1^2)
   z   = z_t + eps*(c1 + c2*eps)
   amp = 1/sqrt(c1 + 2*c2*eps)                            ! replaces C_t**(-1/6)
   ```

   `c1` and `amp` stay exactly mutually consistent (`amp = (dz/dchi)^(-1/2)`), and the
   `**(-1/6)` fractional power is replaced by a square root.

2. **Branch cuts in `invert_flat_action` re-placed.** `p = (3q)^(1/3)` cuts moved from
   `1.8`/`2.0` to `1.72`/`1.97`, where the fitted polynomial and the asymptotic form
   have equal error against an exact root solve.

3. **Redundant work removed.** `S_K(achi)` was evaluated twice per `phi_olver` call
   (once inside `compute_olver_z_amp`, once in `olver_value` to divide), and the
   recurrence fallback multiplied `phi_recurs` by `S_K` only for the caller to divide
   it straight back out. It is now evaluated once in `olver_value` and threaded down.
   The two cheap endpoint metrics are also tested before `open_alpha_cut`, which costs
   a real-exponent power; the accepted region is unchanged (both conditions were
   already ANDed).

4. **Dead code and duplication.** Removed the always-`.false.` `raw` argument, the
   uncalled public `olver_coordinate` and `phi_olver_smallchi`, `olver_smallchi_reduced`
   (inlined into `olver_reduced`, which removes a duplicated `bjl` tail), a duplicate
   local `CACHE_EPS` (now imported from `HypersphericalBesselUtils`), a dead `l > 0`
   test, and an unused `use MpiUtils`. `turn_scale` is now the single `turning_scale`
   helper used in all three places. `real(l*(l+1), dl)` in
   `compute_olver_z_amp_smallchi` became `real(l,dl)*real(l+1,dl)` (the integer product
   overflows for `l >= 46341`; `compute_olver_z_amp` already did it safely).

5. **Documented preconditions** on `phi_olver`/`u_olver`: `l >= 0`, `chi >= 0`, and for
   `K=1` an integer `nu > l` (hence `nu > sqrt(l(l+1))`, so `asin(sqrt(l(l+1))/nu)`
   exists and the closed gate denominator `2(nu-l)` is at least 2). These are validated
   by the Python wrapper; nothing in the Fortran enforces them.

`const_pi`/`const_twopi` from `constants.f90` are now used in
`fortran/hyperspherical_bessels_smallnu.f90` and `fortran/tests/phi_olver_gate_validation.f90`
in place of local `pi` parameters.

## Motivation / measurements

The near-turning branch was producing a localised error spike, 20-50x above the
surrounding error and right at the module's `1e-4` peak-normalized target, because the
window half-width `1e-4*max(chi_t, z_t)` mixes a chi-scale with a z-scale: in open space
at small `alpha = nu/l` we have `z_t >> chi_t`, so the window is wide enough that a
frozen amplitude is visibly wrong, and it sits inside the *first Airy lobe*, i.e. at the
peak. Measured against `phi_recurs` (max `|dphi|/peak` inside the window):

| l | alpha (K=-1) | before | after | just outside window |
|---|---|---|---|---|
| 3000 | 0.0955 | 9.14e-5 | **3.2e-6** | 5.1e-6 |
| 3000 | 0.11 | 7.58e-5 | **2.4e-6** | 4.3e-6 |
| 1000 | 0.11 | 9.32e-5 | **1.0e-5** | 1.3e-5 |
| 10000 | 0.0955 | 9.78e-5 | **6.5e-7** | 1.8e-6 |

It is live because `open_alpha_cut` floors at `0.095`, so `alpha` just above that takes
the raw-Olver path. Component errors at the window edge, against a 45-dps solve of the
map (`K=-1`, `alpha=0.0955`): amplitude `1.46e-4 -> 4.3e-8`, `z` `3.2e-8 -> 6.2e-12`.

Narrowing the window instead does *not* work: `analytic_amplitude` is catastrophically
ill-conditioned there, and dropping the branch entirely gives O(1) errors. The branch is
required; it just needed one more order.

Branch cuts, scored against an exact root solve of `q = t - tanh t` and
`q = tan(theta) - theta` at 40 dps, worst relative error in `u` over the whole line:

| branch | old cut | old worst | new cut | new worst |
|---|---|---|---|---|
| evanescent | p=1.8 (q=1.944) | 5.4e-6 | p=1.72 (q=1.696) | **1.5e-6** |
| oscillatory | p=2.0 (q=2.667) | 2.1e-8 | p=1.97 (q=2.548) | **8.5e-9** |

Speed, from the timing harness in `fortran/tests/phi_olver_gate_validation.f90`
(single thread, min of 4 runs):

| region | before | after |
|---|---|---|
| open high-L accepted | 1.58e-7 s/eval | **1.24e-7** (-21%) |
| closed high-L accepted | 1.22e-7 | **9.6e-8** (-21%) |
| closed mid-L accepted | 1.02e-7 | **9.4e-8** (-8%) |

## Accuracy

- `fortran/tests/phi_olver_gate_validation.f90` over 711530 points, before vs after:
  identical worst peak-normalized error (`1.157e-4`), identical count of warnings above
  `1e-4` (8), zero failures above the `2e-4` target. The calibrated grid does not sample
  the narrow near-turning windows densely enough to see the spike; the targeted scans
  above do.
- `python -m unittest camb.tests.camb_test`: 20 tests OK (117.8 s).
- `python -m unittest camb.tests.mathutils_test`: OK.
- `-Wall -Wextra` clean on the Olver module.

## Notes

- Independently verified while reviewing, and left unchanged because they are correct:
  `qintegral_exact` agrees with quadrature to `<1e-35` relative for all three `K`, both
  branches, over the full domain (including `chi -> 0`, `chi -> chi_t`, `chi -> pi/2`,
  `alpha` from `0.1` to `1e5`), with correct `atan2` branch selection; the `D` series in
  `compute_olver_z_amp_smallchi` is exactly `d(chi*F)/dchi`; `analytic_amplitude` is
  exactly `(dz/dchi)^(-1/2)`; and the `invert_flat_action` polynomial fits do reproduce
  the analytic small-`p` behaviour `1 -/+ p^2/2 + 0.075 p^4` to ~1e-5 in the coefficients.
- A loop-invariant memo for `ell`, `alpha`, `z_t`, `chi_t` keyed on `(l, K, nu)` was
  built and benchmarked out of tree, and **not** adopted: it gave at most ~5%, within
  run-to-run noise, and would have required OpenMP `threadprivate` mutable state plus
  dropping `pure` from `compute_olver_z_amp`.
- `docs/changelog/hyperspherical_bessels.tex` updated (eq. `turning-limit` extended to
  second order plus a new `turning-expansion`, and the `invert_flat_action` branch-cut
  placement); `docs/changelog/Olvers_hyperspherical_bessel.md` updated for the removed
  `raw`/`phi_olver_smallchi` symbols, the corrected module location of
  `invert_flat_action`, and to flag its historical gate constants as superseded.
