# RECFAST explicit solver: first-order lag corrections

- Date: 2026-10-30
- Scope: RECFAST with `use_rosenbrock = False`, the path that integrates with `RungeKuttaDP45` alone. No new
  settings. With `use_rosenbrock = True`, the default, the history and the spectra are unchanged.

Here x_H and x_He are RECFAST's ionized fractions of hydrogen and helium, n_H+/n_H and n_He+/n_He.

## Symptom

RECFAST sets helium and hydrogen from their Saha equilibrium values on its Saha rows, and starts integrating each one
only once it falls below 0.99. The true fractions trail the falling equilibrium values. RECFAST also starts the
matter temperature at T_m = T_rad, after the gas has already cooled slightly below T_rad. With
`use_rosenbrock = False`, hydrogen is still held at Saha through helium recombination, and those rows leave helium's
electrons out of the hydrogen Saha relation.

Against the reference described below, the explicit solver's peak x_e error is 5.7e-4 at z ≈ 1900, against 5.8e-5
for the Rosenbrock solver.

## Change

For an ionized fraction x with equilibrium value x_eq, the first-order lag is delta = (dx_eq/dz)/J, where
J = d(dx/dz)/dx at x_eq is the rate per unit redshift at which x relaxes. With `use_rosenbrock = False`, RECFAST now
makes four changes.

- The helium Saha rows add the lag to x_He where the equilibrium value is below 0.999.
- The hydrogen Saha rows count helium's electrons, move x_H along the Saha slope through each step and add the lag. A
  row whose value is not above 0.99 is redone with the full hydrogen equation. Below x_H = 0.99 hydrogen uses the
  full Peebles rate, also in the band 0.985 < x_H < 0.99 where RECFAST used the bare recombination rate.
- T_m is the tightly coupled solution T_m = T_rad (1 - H t_C/(1 + 2 H t_C)), with t_C the Compton cooling time, on
  the Saha rows and on the integrated rows until H t_C reaches 3.7e-4, at z ≈ 855 in LCDM. RECFAST then integrates
  the full T_m equation. This solution needs H(z), which it takes from CAMB's background for the model being run, but
  no dH/dz. The tight-coupling equation that the Rosenbrock path keeps takes dH/dz from a fit to H(z) in matter,
  radiation and curvature only.
- `Recombination_saha_z`, the redshift below which `Evolve_delta_xe` evolves the perturbation to x_e, is recorded at
  the last output node with x_H > 0.985. The explicit solver needs this record: the bare-rate branch for
  0.985 < x_H < 0.99, which the full Peebles rate replaces, was where RECFAST recorded that redshift. The Rosenbrock
  path keeps its own records.

The thresholds are fixed values in `recfast.f90`: 0.99 for the end of the helium and the hydrogen Saha rows, as
before, and 3.7e-4 for H t_C.

## Before and after

The model is a Planck-like LCDM (H0 = 67.4, ombh2 = 0.0224, omch2 = 0.12, mnu = 0.06 eV) at CAMB's default accuracy
(lmax 2500, Nz 2046). Each row is measured against the reference.

| solver | peak abs(Delta x_e), 400 < z < 3400 | theta*/theta*_ref - 1 | lensed TT, max abs(Delta C_l/C_l), 30 <= l <= 2500 |
| --- | ---: | ---: | ---: |
| `use_rosenbrock = True`, unchanged | 5.8e-5 | 4.5e-7 | 3.7e-6 |
| `use_rosenbrock = False`, before | 5.7e-4 | 4.6e-7 | 6.3e-5 |
| `use_rosenbrock = False`, now | 1.7e-5 | -2.0e-8 | 1.1e-6 |

With the change, RECFAST and the background take 0.92 of the Rosenbrock path's time. That is the median of ten rounds
of interleaved timing on one machine, and the rounds range from 0.89 to 0.93.

In 16 `AxionEffectiveFluid` backgrounds, each measured against its own reference, the explicit solver's x_e error
stays between 1.67e-5 and 1.69e-5, and its theta* error is at most 2.0e-8. The backgrounds take w_n = 1/3 and 1,
f_de(z_c) = 0.05 and 0.1, and z_c = 1000, 1700, 3500 and 5000. Its lensed TT error is 7.2e-7 to 1.7e-6 in 13 of
them. It is 7.7e-6, 4.2e-6 and 5.5e-6 in the three with w_n = 1/3 and (f_de(z_c), z_c) = (0.05, 3500), (0.1, 3500)
and (0.1, 5000), for the reason in the next section. The Rosenbrock path's errors in the same backgrounds are 5.8e-5
in x_e, at most 4.4e-7 in theta*, and 2.4e-6 to 9.0e-6 in TT. In the (0.05, 3500) background its TT error, 7.0e-6,
is below the explicit solver's.

## Remaining HeIII step

Both solvers still set x_e = 1 + f_He for 3500 < z < 5000. At z = 5000 the He++/He+ Saha root still leaves about
1.1e-3 of the helium doubly ionized, so x_e drops by 9.3e-5 there. The reference continues that root down to
z = 3500. The step leaves x_e below z = 3400 and theta* unchanged.

It moves lensed TT by 3.3e-7 and lensed EE by 1.4e-6 in LCDM. In the three backgrounds named above it moves TT by
3.5e-6 to 6.9e-6, which accounts for most of the explicit solver's TT error there. The same numbers come out at `TimeStepBoost = 4`, so the
step, and not CAMB's time sampling, sets them.

## Reference

The reference is the explicit solver with these corrections at much stricter values:
- the helium and hydrogen Saha rows end at x_He = 0.999999 and x_H = 0.99999;
- the full T_m equation starts at the first integrated step (H t_C threshold 1e-8);
- the He++/He+ Saha root continues to z = 3500;
- the `RungeKuttaDP45` tolerance is 1e-10.

These values are constants in this version, so the reference was computed with an earlier version of this change,
commit 1f6d905 of pull request #214, which exposed them as settings. In that version, loosening every value one step
moved x_e by 6.1e-9 and TT by 4.8e-10. The same values with the Rosenbrock solver gave an x_e within 1.3e-8 of the
reference.

The reference that `recfast_rosenbrock.md` tunes against agrees with this one to 1.6e-6 in x_e once its helium Saha
rows also end at x_He = 0.99999. That reference uses full-range Rosenbrock with Nz 160000. With its helium Saha rows
ending at 0.99, it is 5.1e-5 lower at z ≈ 2815.

## Validation

- With `use_rosenbrock = True`, outputs match those of a CAMB 2.0.5 build without this change bit for bit in 562,586
  saved values from 20 runs, parameter names included. The runs vary the RECFAST and accuracy settings,
  `rosenbrock_handoff_xH` 0.995 and 0.996 included, and include `Evolve_delta_xe`. They cover LCDM, a curved w0wa
  model with massive neutrinos and an `AxionEffectiveFluid` model. Each run records background histories and derived
  parameters, and seven of them also record CMB spectra.
- With `rosenbrock_handoff_xH` at 11 values from 0 to 0.999, each with `Evolve_delta_xe` off and on, the
  `use_rosenbrock = True` outputs match the pull request's first version, whose new settings were all off by
  default, bit for bit in 1,381,116 numbers.
- `testRecfastRosenbrockAgreement` now requires the two solvers to agree in x_e to 1e-4, instead of 1e-3. It allows
  5e-6 in T_b, instead of 1e-6. Only the Rosenbrock path still uses the dH/dz fit, so near z = 850 the two T_b
  differ by up to about 1.5e-6. Near z = 825 a step error of `RungeKuttaDP45` in either solver raises that to about
  3.6e-6, and `IntTolBoost = 2` removes it.
- `testRecfastExplicitSolverSahaTemperature` checks that T_m stays below T_rad where the helium Saha rows end.
