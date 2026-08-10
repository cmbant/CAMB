# Halofit non-linear scale solver and HMcode interpolation cleanup

## Summary

- **Standard halofit `R_nl`**: replaced the bisection for `sigma(R)=1` with a bracketed Newton
  iteration. `wint` already returns `d1 = dln sigma^2/dln R`, which is exactly the derivative
  needed (a log-derivative is base independent, so it is also `dlog10 sigma^2/dlog10 R`); the
  step falls back to bisection whenever it leaves the bracket or fails to halve the step size,
  so termination is still guaranteed. The convergence tolerance was tightened from
  `|sigma-1| <= 1e-3` to `1e-4` (`tol_r_nl`): the old tolerance left `R_nl` — and hence `P_nl`
  by up to ~0.3% — dependent on which iterates the solver happened to visit.
- The fixed `wint` quadrature wavenumbers are now a module parameter (`wint_k`) shared by
  `wint` and `wint_pk_table` instead of being rebuilt on every call, and `wint` accumulates the
  three moments from one `exp` per point.
- **One-halo term**: `p_1h` accumulates its trapezium sum as the halo windows are evaluated,
  dropping the `lut%n` temporary and the generic `inttab` call.
- **Growth ODE**: `ode_growth` (fixed-step solve, doubled, re-solved, compared, repeated; with
  unused Euler and midpoint branches) is replaced by `solve_growth_ODE`, which integrates
  directly onto the 64 output points with RK4. The intermediate `a_tab/d_tab/v_tab` arrays and
  the interpolation pass onto the growth grid are gone, and `gnorm` is now just the solution at
  `a=1`. The accumulated growth `G(a)=int_0^a g/a' da'` is now accumulated interval by interval
  instead of re-integrating from `a_min` for every output point (the old loop was `O(n^2)` in
  cubic fits).
- **NFW profile**: the scale radius `rv/c` and the profile normalisation
  `1/(log(1+c)-c/(1+c))` are cached per halo in `HM_tables` when the concentrations are filled,
  rather than recomputed at every `k`. `win`, `winnfw` and `mass` are collapsed into one
  function.
- **Interpolation**: the generic `find(x,xtab,ytab,n,iorder,ifind,imeth)` (order 1/2/3 x
  equal-spacing/bisection x coefficient/Lagrange, plus a reverse-the-arrays recovery path) is
  replaced by the three combinations actually configured: `interp_linear_uniform`,
  `interp_cubic_uniform` and `interp_cubic_sorted`. `find`, `table_integer`, `fit_line`,
  `fit_quadratic`, `Lagrange_polynomial`, `reverse` and ~20 `iorder_*`/`ifind_*`/`imeth_*`
  parameters are removed.
- **Feedback response**: `HMcode` no longer allocates `p_den(nk,nz)` and `p_num(nk,nz)`. The
  denominator is a `p_den(nk)` local to `HMcode_redshift` and the response is applied in the
  numerator pass, so the post-loop array expression and the allocatable 2-D dummy arguments go.
- Dead runtime selectors removed alongside: `inttab` keeps only the cubic rule and `integrate`
  only Simpson's rule (the only orders any caller ever requested).
- Smaller cleanups: the `timing_test` block in `HMcode` (dead debug scaffolding ending in a
  bare `stop`) is removed; `find_pk` no longer repeats its high-k extrapolation for each matter
  type; `init_wiggle` drops a redundant `nk` temporary; `calculate_nowiggle` drops an argument
  it only used for a size check; `int_split` bisects in integer arithmetic; the remaining bare
  `stop` is an `error stop`; `allocate_LUT` no longer zeroes every `HM_tables` array/scalar on
  every redshift, since each is either fully overwritten in the same call (by `halomod_tables`,
  `fill_conc`, or unconditionally in `halomod_init` right afterwards) or is a feedback-only
  field only ever read by the `imead` variant that fills it; `smooth_array_Gaussian` drops the
  runtime check for equally-spaced `x`, since its only caller (`init_wiggle`'s log-k grid, built
  by `fill_table`) is always exactly uniform; `rebin_pk` (a `.true.` compile-time parameter whose
  `.false.` branch was an immediate `error stop`, i.e. dead) and its branch are removed from
  `fill_plintab`, which now always builds the equal-log-spaced `k` table directly.

Net: 841 lines removed, 312 added. No public API or parameter changes.

## Numerical effect

Compared against the previous code over 6 cosmologies (including massive neutrinos and `w0-wa`),
10 halofit versions, `z = 0, 0.3, 1, 2.5, 6` and `k/h` in `[1e-3, 15]`:

| Change | Max fractional change in `P_nl` |
| --- | ---: |
| One-halo fusion, feedback arrays, NFW caching, smaller cleanups | `1.1e-14` |
| Growth ODE + interpolation rewrite | `1.9e-6` |
| `R_nl` solver (standard halofit only) | `3.2e-3` |

The `R_nl` shift is entirely the tolerance change, not the solver: building the new code with the
Newton step forced to bisection at a converged `1e-8` tolerance reproduces the Newton results to
`5.5e-8`, i.e. both find the same root. At the chosen `1e-4` the answer is within `2.5e-4` of that
converged root, against ~`3e-3` for the old `1e-3`. HMcode versions do not use this code path.

The growth-ODE change is a convergence improvement, not a re-tuning. Measured against the same
code run with `nsub_growth_ODE=64` (same `a_ini`), the truncation error of the growth solution
propagated into `P_nl` is:

| Cosmology | new (`nsub=8`) | previous `ode_growth` |
| --- | ---: | ---: |
| LCDM | `3.6e-10` | `2.2e-8` |
| massive neutrinos (0.06 eV / 0.3 eV) | `3.6e-10` | `2.4e-8` |
| `w0 = -0.9, wa = 0.2` | `3.7e-8` | `1.9e-6` |

`w0-wa` is where the two differ most, because the Dolag (2004) concentration correction takes
ratios of the growth in the model and in a matched LCDM one; almost all of the `1.9e-6` entry in
the first table is the old solver's error there. Changing `a_ini` from `1e-4` to `1e-5` moves both
codes by a common `~1e-7` (the growing-mode initial condition, which is unchanged).

`camb/tests/camb_test.py` `testPowers` has a hard-coded `halofit_version="takahashi"` reference
value, updated 57.723 -> 57.734 accordingly.

## Cost

Evaluations of `wint` for 64 redshifts: 655 before (to `1e-3`), 330 after (to `1e-4`) — half as
many, at a tolerance ten times tighter. `wint` is not the dominant cost of standard halofit, so
wall-clock timings for the non-linear correction alone are within noise (~0.4 s for 64 redshifts,
4 threads) for all versions.

`fill_growtab` (growth ODE plus accumulated growth, run once per HMcode call) takes 0.6 ms,
against ~270 ms for the non-linear correction at a single redshift; it was 1.2 ms before the
accumulated-growth loop was made linear. Running the ODE with 8x the steps is not measurable in
the total. A general adaptive integrator (`RungeKuttaDP45`) would therefore buy nothing here, and
could not be called without either putting the working `HM_cosmology` on `THalofit` — which the
parallel redshift loop relies on not happening — or type-punning it through the solver's opaque
first argument.

## Tests

```text
python setup.py make
python -m unittest camb.tests.camb_test    # 22 tests OK, 99 s
python -m camb.check_accuracy inifiles/planck_2018.ini   # PASS
pre-commit run --files fortran/halofit.f90 camb/tests/camb_test.py   # all passed
```
