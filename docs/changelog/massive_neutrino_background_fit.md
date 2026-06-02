# Massive Neutrino Background Fits

## Summary

Thermal massive-neutrino background density and pressure now use direct smooth
fits for the intermediate range `0.3 < a*m_nu < 70`. The small- and large-mass
series branches are unchanged outside this interval.

The fit variable maps the intermediate range to `-1 <= z <= 1`, with separate
low-order Horner polynomials for `rho`, `P`, and `D = rho - 3P`. `rho` and
`rho_P` share the same density fit, while `drho` evaluates the direct fit for
`D` used in `dot rho = H * D`.

## Cleanup

- Removed the cubic-Hermite thermal background table storage and initialization.
- Removed the legacy 100-step direct quadrature helper.
- Removed the optimized 8-point quadrature helper used only to fill the old
  table.
- Removed the now-unneeded table deallocation from `CAMB_FreeGlobalMemory`.
- Removed the old `ThermalNuBack` pointer workaround and the now-empty
  `ThermalNuBackground_init` hook.

## Accuracy

Against direct Fermi-Dirac quadrature on 801 logarithmically spaced points over
the fit range, the current low-order fits give approximately:

- `rho`: max relative error `2.9e-5`, RMS relative error `2.0e-5`
- `P`: max relative error `5.8e-5`, RMS relative error `4.1e-5`
- `D = rho - 3P`: max relative error `8.9e-5`, RMS relative error `6.3e-5`

For context, the removed 8-point quadrature table generation was much more
accurate for the table values themselves, with max relative errors around
`1.4e-7` in `rho` and `4.9e-7` in `P` against direct quadrature. The direct fits
trade this surplus table-generation accuracy for less setup, less global state,
and faster repeated background evaluations.

## Timing Notes

In local `lmax=4000` runs with standard massive neutrinos, after discarding the
first CAMB call, the direct-fit version was slightly faster than the spline
baseline in 20-run CPU timing:

- old spline baseline: mean `1.015s`, median `1.010s`, total `20.310s`
- direct low-order fits: mean `0.977s`, median `0.966s`, total `19.548s`

A focused `rho_P` microbenchmark showed a larger change in the hot path itself,
from about `0.419s` for the old spline interpolation to `0.350s` for the direct
fit over the same call pattern.
