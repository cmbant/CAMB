# Non-flat integration updates

This note summarizes the current non-flat integration changes in `fortran/cmbmain.f90`.

## Main algorithm changes

- Replaced the on-the-fly non-flat `u_{\nu l}` integration from RK4 stepping at every source sample with a Numerov-based update, using one RK4 bootstrap step followed by Numerov steps.
- Refactored both scalar and tensor non-flat line-of-sight integration so they compute `u_{\nu l}` values on the source time grid and then integrate against `Source_q` there.
- Removed interpolation of non-flat sources at internal Numerov substeps.
- For the scalar oscillatory-region sweep, added a high-`x` cutoff so the earliest, very oscillatory tail is skipped unless `full_bessel_integration` or bispectrum calculations require it.
- Kept the source-interval landing exact by splitting each source interval into an integer number of Numerov substeps, with
  `nSubSteps = ceiling(dchisource / dchimax)` and `delchi = dchisource / nSubSteps`.

## Code simplifications

- Removed the non-flat source second-derivative storage `ddSource_q`, which is no longer needed now that both scalar and tensor paths integrate on the source grid.
- Simplified the scalar inner accumulation loop by applying trapezoidal endpoint weights to the stored `u_{\nu l}` array before accumulation, and unrolled the common 3-source case.

## Accuracy and performance checks

- Against the stored RK4 scalar baseline over `Omega_k = \pm 0.02, \pm 0.002`, the current scalar implementation gives lensed `TT` max fractional differences about `8.4e-4 .. 9.6e-4` and `EE` about `3.9e-4 .. 4.4e-4` for `\ell >= 2`. For `\ell >= 30`, `TT` stays about `8.4e-4 .. 9.6e-4` while `EE` drops to about `1.62e-4 .. 1.78e-4`.
- In the same scalar runs, CPU time dropped from about `25 .. 32 s` for the stored RK4 baseline to about `11.5 .. 15.1 s` for the current code, roughly a `52% .. 55%` reduction.
- The scalar high-`x` cutoff itself gives an additional CPU reduction of about `2.6% .. 12.5%` versus the uncapped source-grid Numerov code over the same `Omega_k` cases.
- Relative to that uncapped source-grid Numerov code, the scalar cutoff changes spectra mainly at very low `\ell`: for `\ell >= 2`, max differences are about `TT = 4.4e-4 .. 6.5e-4`, `TE = 5.7e-7 .. 9.2e-7`, and `PP = 6.2e-6 .. 7.1e-6`. For `\ell >= 30`, the same comparison shrinks to about `TT = 1.8e-6 .. 3.7e-6`, `EE = 6.6e-11 .. 9.1e-11`, `TE = 5.9e-9 .. 1.3e-8`, and `PP = 1.0e-9 .. 3.0e-9`.
- Direct comparison of lensing potential `PP` against the pre-change code shows the largest fractional differences at very low `\ell`, peaking around `\ell = 3` at about `3.1e-3 .. 3.5e-3` over tested `Omega_k = \pm 0.02, \pm 0.002`.
- For `\ell >= 30`, the `PP` differences versus the pre-change code are smaller, with max fractional differences about `7.4e-4 .. 8.5e-4` and RMS about `2.4e-4 .. 2.7e-4`.
- In the same `PP` runs, the current code was faster than the pre-change code by roughly `20% .. 33%` in CPU time over the tested curvature values.

## Scope

These changes are localized to the non-flat on-the-fly hyperspherical Bessel integration in `cmbmain.f90`. Flat-space integration and precomputed flat Bessel interpolation are unchanged.
