# RECFAST lag corrections

`scripts/recfast_lag_corrections.py` regenerates the numbers in this entry, apart from the timing medians, by comparing
RECFAST with `lag_corrections` off and on against a converged reference CAMB.

- Date: 2026-10-23
- Scope: optional RECFAST setting `lag_corrections` (ini key `RECFAST_lag_corrections`), off by default, and the four
  settings it reads, `Tm_handover_frac`, `hydrogen_saha_xH`, `Saha_He_end_xHe` and `rk45_tol`.

Here x_H and x_He are RECFAST's ionized fractions of hydrogen and helium, n_H+/n_H and n_He+/n_He. `lag_corrections`
also turns on `heIII_saha_continued` and `snap_saha_z`, described in the two entries of the same date. Either of those
can be set on its own, but `lag_corrections` needs `snap_saha_z`. With `use_rosenbrock = F` its hydrogen rows
integrate 0.985 < x_H < 0.99 with the full rate, where the right-hand side no longer records `Recombination_saha_z`,
so without `snap_saha_z` that redshift would stay at its starting value of 0.

## Symptom

With the setting off, RECFAST holds x_He at its value solved by the Saha equation until it falls to 0.99, and then
integrates the helium rate equation from there. However, equilibrium is not reached instantly as the equilibrium value
changes over time, thus helium is slightly more ionized, i.e. x_He lags behind the equilibrium value. RECFAST also
starts integrating with T_m = T_rad, after the gas has already cooled slightly below T_rad. With the setting off, the x_e error
reaches 5.8e-5, or 5.4e-5 as a fraction of x_e, about the size of RECFAST's remaining difference from CosmoRec in the
helium window (3.5e-5 RMS and 7.7e-5 at most, after the helium-rate fit in
`2026-06-18-recfast-helium-rate-correction-fit.md`, which was made with this setting off). These offsets set a floor
that making RECFAST's settings stricter does not get below, and two refinements make the error larger.

- Lowering `rosenbrock_tol` from 3e-4 to 3e-6 makes the lensed TT error 1.5 times larger (3.699e-6 to 5.696e-6).
- Refining the RECFAST grid 8 times makes it 1.58 times larger (3.705e-6 to 5.846e-6, both against a grid-converged
  reference CAMB at Nz 32736).
- Across the four settings the script tries (helium threshold 0.999, `rosenbrock_tol` 3e-6, `rosenbrock_handoff_xH`
  0.99 and Nz times 8), the smallest TT error is 3.252e-6, at the helium threshold of 0.999. Here
  `rosenbrock_handoff_xH` is the x_H at which RECFAST switches from the Rosenbrock solver back to `RungeKuttaDP45`.
  The 9.3e-5 HeIII step stays, since it comes from the cut at z = 5000.

The T_m = T_rad initial condition makes T_m too high by 2.2e-3 K all the way to z ≈ 1100. CAMB's tight-coupling
temperature equation (in `EvaluateRecfastODE`) has no term that pulls T_m back, as d(T_m − T_rad)/dz depends on
T_m only through the ionization rates. The script prints the offset at z ≈ 1100. Starting helium from its Saha
value at x_He = 0.99 contributes to the largest x_e error in CAMB below z = 3400, which is 5.790e-5 at z = 2810.

## Change

For an ionization fraction x with equilibrium value x_eq, the first-order lag is delta = (dx_eq/dz)/J, where
J = d(dx/dz)/dx at x_eq is the rate per unit redshift at which the species relaxes. With the setting on, RECFAST
adds, to first order, how far helium, hydrogen and T_m lag behind their equilibrium values.

- The helium Saha rows add the lag to x_He.
- With `use_rosenbrock = F` the hydrogen Saha rows count helium's electrons, follow the Saha slope through each step
  and add the lag. A row whose value is not above `hydrogen_saha_xH` is redone with the full hydrogen equation.
- With `use_rosenbrock = T` hydrogen is integrated in full from the end of the helium Saha rows.
- T_m is set to the first-order tight-coupling solution of its equation, T_m = T_rad (1 - H t_C/(1 + 2 H t_C)), with t_C the
  Compton cooling time, until H t_C passes `Tm_handover_frac`, and RECFAST then moves T_m to its full equation. The
  default 3.7e-4 hands over at z ≈ 855 in LCDM, close to the z ≈ 847 where RECFAST leaves tight coupling with the
  setting off.

This solution needs H but no dH/dz, and it takes H as CAMB computes it for the model being run, early dark energy
included. With the setting off, RECFAST evolves T_m during tight coupling with an equation that needs dH/dz (in
`EvaluateRecfastODE`), and it takes dH/dz from a fit to H(z) in matter, radiation and curvature, which leaves out
components such as early dark energy.

With `use_rosenbrock = F` the setting also retires two special cases, the band 0.985 < x_H < 0.99 where RECFAST drops
the Peebles factor and relaxes hydrogen at the bare case-B rate, and the frozen hydrogen rate above it.

With the setting on, a positive `Saha_He_end_xHe` below 0.99 is refused by `validate()` and by an ini run, because
the first-order helium lag overshoots there.

## Reference CAMB

The reference CAMB (the script's "ground truth") is a very high-accuracy version with `lag_corrections` on. It
switches to the full equations earlier, and its numerical error tolerance is lowered until the outputs stop changing.
There, x_He and x_H are computed from the Saha equation until they reach 0.999999 and 0.99999 respectively, and T_m is
solved with the full equation from the end of the first integrated step, where it starts from its first-order
tight-coupling value (`Tm_handover_frac = 1e-8`). Rosenbrock is turned off (`use_rosenbrock = F`), and `rk45_tol` is
set to 1e-10. The script prints three checks on it.

1. Loosening the three thresholds together, each by one step (to 0.99999, 0.9999 and 1e-7), moves x_e by 6.145e-9
   and lensed TT by 4.772e-10.
2. Rosenbrock on and off differ by 1.254e-8 in x_e and 1.263e-9 in TT.
3. Setting `Tm_handover_frac` to 0, a value `validate()` refuses but `camb.get_results` runs, leaves the history the
   same bit for bit, so at 1e-8 the full T_m equation already starts as early as the code allows, at the end of the
   first integrated step.

The high-accuracy run that `recfast_rosenbrock.md` tunes the Rosenbrock settings against (full-range Rosenbrock,
`rosenbrock_handoff_xH` 0, `rosenbrock_tol` 1e-7, Nz 160000, `lag_corrections` off) also computes helium from the Saha equation until x_He
falls to 0.99, and starts its Rosenbrock integration only there. Solved on the same grid, its x_e at z ≈ 2815 is 5.1e-5 lower than in the reference CAMB, a fractional
error of about 4.7e-5 and about twenty times the recombination-era fractional x_e error of 2.29e-6 that
`recfast_rosenbrock.md` reports for CAMB at `BackgroundTimeStepBoost = 8`. With only that choice
moved to 0.99999, the gap there falls to 5.2e-9 (the largest remaining difference, 1.6e-6, is at z ≈ 1305), so the
same convergence test, with that one change, would also capture this lag.

## Before and after

The model is a Planck-like LCDM (H0 = 67.4, ombh2 = 0.0224, omch2 = 0.12, mnu = 0.06 eV) at CAMB's default accuracy
(lmax 2500, Nz 2046), and each row is measured against the reference CAMB on the same grid.

| setting | thresholds (He, H, H t_C) | peak x_e error, max abs(Delta x_e), 400 < z < 3400 | theta* error, theta*/theta*_ref - 1 | lensed TT error, max abs(Delta C_l/C_l), 30 <= l <= 2500 | time |
| --- | --- | ---: | ---: | ---: | ---: |
| option off, CAMB's own settings | 0.99, 0.99, - | 5.8e-5 | 4.5e-7 | 3.7e-6 | 1.00 |
| option on, `use_rosenbrock = F` | 0.99, 0.99, 3.7e-4 | 1.7e-5 | -2.0e-8 | 8.9e-7 | ~0.94 |
| stricter, at equal cost | 0.995, 0.995, 1.25e-4 | 2.2e-6 | -7.0e-9 | 2.4e-7 | ~0.98 |
| stricter still, `rk45_tol = 1e-8` | 0.999, 0.997, 4.5e-5 | 9.0e-8 | -2.3e-9 | 6.7e-8 | ~1.23 |

The thresholds are the x_He and x_H at which helium and hydrogen leave their Saha phases, and the H t_C at which T_m
goes over to its full equation. Time is RECFAST plus background relative to the first row, the median of 40
interleaved timing runs on one machine. Single runs scatter by about 2% around it.

- At CAMB's own thresholds, with `use_rosenbrock = F`, the option cuts the theta* error 22.8 times, the peak x_e
  error 3.5 times (5.5 times counting the HeIII step) and the TT error 4.2 times, and it runs 4 to 8% faster.
- At equal cost, the x_e error is 26.3 times and the TT error 15.7 times smaller than with the option off.
- With the option, the same `rosenbrock_tol` change makes TT 4.3 times better (3.919e-6 to 9.082e-7, with
  `use_rosenbrock = T` at handoff 0.976), and the 8 times finer grid changes it by a factor of 1.04 (8.745e-7 to
  8.383e-7).
- With `use_rosenbrock = T`, the option has a TT error of 3.919e-6 at the default `rosenbrock_handoff_xH` of 0.976,
  slightly above the 3.699e-6 with the option off. At `rosenbrock_handoff_xH = 0.995` the TT error falls to 9.343e-7,
  close to the 8.912e-7 with `use_rosenbrock = F`, for about 5 to 6% more time (1.045 to 1.063 in five runs of the script
  on one machine, each timing the best of ten calls).

## Early dark energy

With the setting off, RECFAST's tight-coupling temperature equation takes dH/dz from a fit in matter, radiation and
curvature, which leaves an early dark energy fluid out. `python scripts/recfast_lag_corrections.py --ede` measures both
paths in 16 `AxionEffectiveFluid` backgrounds, with w_n = 1/3 and 1, f_de(z_c) = 0.05 and 0.1, and z_c = 1000, 1700,
3500 and 5000, each against a reference CAMB solved in that background. The scan covers that fluid family at those 16
points.

- With the setting off, the lensed TT error ranges from 2.444e-6 to 8.954e-6. It is above the LCDM value of 3.699e-6
  in 9 of the 16, and largest for w_n = 1/3 and f_de(z_c) = 0.1 at z_c = 1700, at 2.42 times the LCDM value.
- With the setting on, it stays between 8.419e-7 and 9.656e-7, against 8.912e-7 in LCDM.
- With the setting off, the theta* error ranges from 7.2e-8 to 4.4e-7, below the LCDM 4.54e-7 in each of the 16. T_m
  starts slightly too high, and the fit's missing term makes T_m lower than the full dH/dz would. We take the smaller
  theta* error to mean that these two errors in CAMB partly cancel, which the script does not test on its own. If they
  do, correcting the initial condition while keeping the fit would expose the error from the fit.
- For w_n = 1/3, f_de(z_c) = 0.1 and z_c = 1700, the reference CAMB solved with Rosenbrock on and off differs by 1.3e-8
  in x_e and 1.4e-9 in TT, as in LCDM.

## Validation

- With the new settings at their defaults, outputs match CAMB 2.0.5 bit for bit in 392,476 numbers from 18 runs, which
  vary the RECFAST and accuracy settings and cover LCDM, a curved w0wa model with massive neutrinos and an
  `AxionEffectiveFluid` model, with background histories and derived parameters in each run and CMB spectra and
  matter power in four of them. They also match in the 541 numerical files `fortran/tests/run_tests.py` writes.
- `testRecfastLagCorrections` covers the setting with both solvers, its thresholds, their validation and the ini
  round trip.
- `python scripts/recfast_lag_corrections.py --all` prints the table and three checks in about 30 s, with the
  supporting comparisons in its sections 1 to 10, and `--ede` adds the early dark energy scan in about 70 s more.
