# RECFAST Rosenbrock Integrator

## Summary

CAMB now has an optional ROS2 Rosenbrock integrator for the RECFAST H/He/Tm system.
It is intended as a fast semi-stiff replacement for the early part of the old DVERK
evolution, with a handoff back to DVERK once hydrogen has recombined enough.

Relevant implementation points:

- Uses a two-stage ROS2 method with analytic Jacobian entries for `x_H`, `x_He`, and `a*T_m`.
- Includes the explicit non-autonomous `df/dz` term in the Rosenbrock stages.
- Integrates the full smooth H/He/Tm equations while Rosenbrock is active, with no Saha shortcut in that phase.
- Snaps the Rosenbrock-to-DVERK handoff to the higher-redshift RECFAST output node by redoing the crossing interval with DVERK.
- Uses an embedded Rosenbrock error estimate rather than full-step versus two-half-step differencing.
- Scales the effective Rosenbrock tolerance as `rosenbrock_tol / CP%Accuracy%IntTolBoost` before solving.
- Scales the RECFAST DVERK input tolerance as `dverk_tol / CP%Accuracy%IntTolBoost`, with `dverk` itself using an internal tolerance of `tol_in / 5`.
- Uses the stored `RECFAST_nz` value as a baseline and converts it to an internal grid with `Nz * CP%Accuracy%BackgroundTimeStepBoost`.
- Uses a smaller ionization error-scale floor than the temperature floor so full-range Rosenbrock does not allow large relative tail errors when `x_e` is small.

## Current Parameters

| Parameter | Current value | Notes |
| --- | --- | --- |
| `use_rosenbrock` | `.true.` | Default RECFAST mode now uses the fast Rosenbrock-to-DVERK handoff. |
| `RECFAST_nz` | `2046` | Stored class default; the internal RECFAST grid uses `Nz * BackgroundTimeStepBoost`. |
| `rosenbrock_handoff_xH` | `0.976` | Fast-handoff default for the new default `Nz = 2046` baseline. |
| `rosenbrock_tol` | `3e-4` | Tuned for the fast handoff path, not for full-range Rosenbrock. |
| `dverk_tol` input | `1.5e-5` | RECFAST passes this to DVERK, then divides by `IntTolBoost`; DVERK internally uses `tol_in / 5`, so this matches the old `3e-6` internal scale at `IntTolBoost = 1`. |
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
- The internal RECFAST `dverk` scale was re-tuned after switching `dverk` to a DP5 pair. The current RECFAST input value is `1.5e-5`, but `dverk` now converts that to an internal tolerance of `3e-6`; RECFAST also divides the input value by `IntTolBoost`, so the DVERK phase tightens together with the Rosenbrock phase.
- The solver now applies `rosenbrock_tol / IntTolBoost`, so increasing CAMB integration accuracy tightens the Rosenbrock substeps as well.
- The stored `Nz` value is treated as a baseline grid and is converted internally using `Nz * BackgroundTimeStepBoost` so RECFAST refines together with CAMB's background time-step boost.
- `RECFAST_rosenbrock_ion_scale_floor = 1e-3` was introduced after checking the low-`x_e` tail, where a unit floor was formally too weak for full-range Rosenbrock.

## Measured Behavior

Against the smooth Rosenbrock reference:

- The most useful accuracy summary is the recombination-era fractional `x_e` error, rather than the global max `x_e` error, because the latter is often dominated by the late low-`x_e` tail.
- Fast handoff, default baseline `Nz = 2046`, `handoff = 0.976`, `tol = 3e-4`: recombination-era fractional `x_e` error about `4.32e-5`, corresponding to absolute `x_e` error about `4.25e-5`, with recombination-era relative `T_b` error about `1.04e-6`.
- With the new DP5-based `dverk`, RECFAST now passes `dverk_tol = 1.5e-5`, which maps back to the same `3e-6` internal DVERK tolerance at `IntTolBoost = 1`.
- Fast handoff, `Nz = 2048`, `handoff = 0.976`, `tol = 3e-4`: recombination-era fractional `x_e` error is about `4.30e-5`, corresponding to absolute `x_e` error about `4.23e-5`.
- Increasing `BackgroundTimeStepBoost` with the corrected internal scaling `Nz * BackgroundTimeStepBoost` gives smooth convergence to the same reference: the default RECFAST mode improves from about `4.32e-5` in recombination-era fractional `x_e` error at boost `1` to about `6.27e-6` at boost `2`, `5.14e-6` at boost `4`, and `2.29e-6` at boost `8`.
- Increasing `IntTolBoost` alone has only a tiny effect for this fast-handoff configuration, which is consistent with the early Rosenbrock plus snapped DVERK handoff already being dominated by the RECFAST output grid rather than by local substep tolerances.
- Full-range Rosenbrock with `tol = 3e-4` is too loose and was not adopted as the default mode.
- Full-range Rosenbrock becomes respectable around `tol = 1e-5`, but this is about `2.3x` to `2.6x` slower than the fast handoff configuration in the tested cases.

In short, the shipped defaults are tuned for the fast Rosenbrock-plus-DVERK handoff mode.
If full-range Rosenbrock is requested, it should normally be used with a tighter tolerance than the default.