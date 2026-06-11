# RECFAST Rosenbrock Integrator

## Summary

CAMB now has an optional ROS2 Rosenbrock integrator for the RECFAST H/He/Tm system.
It is intended as a fast semi-stiff replacement for the early part of the old
non-stiff evolution, with a handoff back to `RungeKuttaDP45` (the DP5 pair that
replaced DVERK) once hydrogen has recombined enough.

Relevant implementation points:

- Uses a two-stage ROS2 method with analytic Jacobian entries for `x_H`, `x_He`, and `a*T_m`.
- Includes the explicit non-autonomous `df/dz` term in the Rosenbrock stages.
- Integrates the full smooth H/He/Tm equations while Rosenbrock is active, with no Saha shortcut in that phase.
- Snaps the Rosenbrock-to-RK45 handoff to the higher-redshift RECFAST output node by redoing the crossing interval with `RungeKuttaDP45`.
- Uses an embedded Rosenbrock error estimate rather than full-step versus two-half-step differencing.
- Scales the effective Rosenbrock tolerance as `rosenbrock_tol / CP%Accuracy%IntTolBoost` before solving.
- Scales the `RungeKuttaDP45` tolerance as `rk45_tol / CP%Accuracy%IntTolBoost`; `rk45_tol = 3e-6` is used directly by the integrator (no internal `/ 5` rescaling, unlike the old DVERK).
- Uses the stored `RECFAST_nz` value as a baseline and converts it to an internal grid with `Nz * CP%Accuracy%BackgroundTimeStepBoost`.
- Uses a smaller ionization error-scale floor than the temperature floor so full-range Rosenbrock does not allow large relative tail errors when `x_e` is small.

## Current Parameters

| Parameter | Current value | Notes |
| --- | --- | --- |
| `use_rosenbrock` | `.true.` | Default RECFAST mode now uses the fast Rosenbrock-to-RK45 handoff. |
| `RECFAST_nz` | `2046` | Stored class default; the internal RECFAST grid uses `Nz * BackgroundTimeStepBoost`. |
| `rosenbrock_handoff_xH` | `0.976` | Fast-handoff default for the new default `Nz = 2046` baseline. |
| `rosenbrock_tol` | `3e-4` | Tuned for the fast handoff path, not for full-range Rosenbrock. |
| `rk45_tol` | `3e-6` | RECFAST passes this to `RungeKuttaDP45` after dividing by `IntTolBoost`; the integrator uses it directly, matching the old DVERK internal scale at `IntTolBoost = 1`. |
| `RECFAST_rosenbrock_ion_scale_floor` | `1e-3` | Protects low-`x_e` ionization components from effectively unconstrained relative errors. |

## How The Parameters Were Chosen

The tuning was done against a smooth internal reference built from full-range Rosenbrock
with `handoff = 0`, `rosenbrock_tol = 1e-7`, and `Nz = 160000`.
That reference was checked for grid convergence by comparing `Nz = 80000` and `Nz = 160000`.

The main fast-path settings were chosen as follows:

- `rosenbrock_handoff_xH` was scanned over snapped handoff candidates and chosen to minimize the recombination-era `x_e` mismatch while keeping `get_background` fast.
- For the default `Nz = 2046` baseline, `0.976` remains the best tested fast-handoff value and is used as the code default.
- For the nearby `Nz = 2048` regression grid, `0.976` gave the best tested recombination-era agreement against the smooth reference and is retained in the regression test.
- `rosenbrock_tol = 3e-4` was retained as the default because it is a good speed/accuracy tradeoff for the fast handoff mode.
- After the switch from DVERK to the `RungeKuttaDP45` DP5 pair (which uses its input tolerance directly), RECFAST passes `rk45_tol = 3e-6`, matching the old DVERK internal scale; RECFAST divides it by `IntTolBoost`, so the RK45 phase tightens together with the Rosenbrock phase.
- The solver now applies `rosenbrock_tol / IntTolBoost`, so increasing CAMB integration accuracy tightens the Rosenbrock substeps as well.
- The stored `Nz` value is treated as a baseline grid and is converted internally using `Nz * BackgroundTimeStepBoost` so RECFAST refines together with CAMB's background time-step boost.
- `RECFAST_rosenbrock_ion_scale_floor = 1e-3` was introduced after checking the low-`x_e` tail, where a unit floor was formally too weak for full-range Rosenbrock.

## Measured Behavior

Against the smooth Rosenbrock reference:

- The most useful accuracy summary is the recombination-era fractional `x_e` error, rather than the global max `x_e` error, because the latter is often dominated by the late low-`x_e` tail.
- Fast handoff, default baseline `Nz = 2046`, `handoff = 0.976`, `tol = 3e-4`: recombination-era fractional `x_e` error about `4.32e-5`, corresponding to absolute `x_e` error about `4.25e-5`, with recombination-era relative `T_b` error about `1.04e-6`.
- With `RungeKuttaDP45`, RECFAST passes `rk45_tol = 3e-6` directly, the same internal tolerance scale as the old DVERK at `IntTolBoost = 1`.
- Fast handoff, `Nz = 2048`, `handoff = 0.976`, `tol = 3e-4`: recombination-era fractional `x_e` error is about `4.30e-5`, corresponding to absolute `x_e` error about `4.23e-5`.
- Increasing `BackgroundTimeStepBoost` with the corrected internal scaling `Nz * BackgroundTimeStepBoost` gives smooth convergence to the same reference: the default RECFAST mode improves from about `4.32e-5` in recombination-era fractional `x_e` error at boost `1` to about `6.27e-6` at boost `2`, `5.14e-6` at boost `4`, and `2.29e-6` at boost `8`.
- Increasing `IntTolBoost` alone has only a tiny effect for this fast-handoff configuration, which is consistent with the early Rosenbrock plus snapped RK45 handoff already being dominated by the RECFAST output grid rather than by local substep tolerances.
- Full-range Rosenbrock with `tol = 3e-4` is too loose and was not adopted as the default mode.
- Full-range Rosenbrock becomes respectable around `tol = 1e-5`, but this is about `2.3x` to `2.6x` slower than the fast handoff configuration in the tested cases.

In short, the shipped defaults are tuned for the fast Rosenbrock-plus-RK45 handoff mode.
If full-range Rosenbrock is requested, it should normally be used with a tighter tolerance than the default.

## Review Fixes And Optimizations (2026-06-11)

A code review of `recfast.f90` led to two correctness fixes and a set of efficiency changes:

- Fixed a Jacobian bug where the singlet `AHcon` derivatives picked up the triplet KIV fit
  constants (`pb = 0.66`, `qb = 0.9` instead of `0.36` and `RECFAST_fudge_He`) left behind by the
  f-evaluation when `Heswitch = 6` (the default). The singlet/triplet fit shapes are now separate
  `pb_s/qb_s` and `pb_t/qb_t` constants. Only the Jacobian (step control) was affected; the fix
  shifts `x_e` by at most `1e-5` relative at `z ~ 1880`, within the Rosenbrock tolerance.
- Rosenbrock failure in `TRecfast_init` now calls `GlobalError` instead of returning silently
  with stale interpolation tables.
- Each ROS2 step now LU-factors the stage matrix once and back-substitutes for both stages,
  instead of running Gaussian elimination twice.
- `f`, the analytic Jacobian, and the finite-difference `df/dz` are computed once per (z, y)
  state and reused across step-size retries; `df/dz` uses a one-sided difference reusing the
  already-computed `f`. Each Rosenbrock attempt now costs 2 ODE evaluations instead of 4.
- Repeated `pow` evaluations in the Jacobian are cached, and the PPB/VF/triplet fit constants
  are now `parameter` declarations.

Net effect: results unchanged beyond the Jacobian fix above (the efficiency changes alone moved
`x_e` by `<1e-9` relative), the full default test suite passes, and `get_background` is roughly
25% faster (about `4.9 ms` to `3.5 ms` per call in the tested configuration).
