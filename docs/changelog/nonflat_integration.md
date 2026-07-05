# Non-flat integration updates

This note summarizes the current non-flat integration changes in `fortran/cmbmain.f90`.

## Main algorithm changes

- Replaced the on-the-fly non-flat `u_{\nu l}` integration from RK4 stepping at every source sample with a Numerov-based update, using one RK4 bootstrap step followed by Numerov steps.
- Refactored both scalar and tensor non-flat line-of-sight integration so they compute `u_{\nu l}` values on the source time grid and then integrate against `Source_q` there.
- Removed interpolation of non-flat sources at internal Numerov substeps.
- For the scalar oscillatory-region sweep, added a high-`x` cutoff so the earliest, very oscillatory tail is skipped unless `full_bessel_integration` or bispectrum calculations require it.
- Kept the source-interval landing exact by splitting each source interval into an integer number of Numerov substeps, with
  `nSubSteps = ceiling(dchisource / dchimax)` and `delchi = dchisource / nSubSteps`.

## Near-flat shifted-ν approximation

- For scalar non-flat runs close to flat, `DoRangeInt` can now fill `u_{\nu l}` directly from
  `u_{\nu l}(\chi) \approx (\chi / S_K(\chi)) j_l(\nu_{\mathrm{eff}} \chi)` with
  `\nu_{\mathrm{eff}}^2 = \nu^2 - K_{\mathrm{sign}} l(l+1) / 3`.
- When the full scalar source range is inside that near-flat regime, `IntegrateSourcesBessels` now bypasses the dissipative/oscillatory split entirely and uses a dedicated `DoNearFlatIntegration` path that follows the same source-grid traversal order and the same high-`x` and Limber cutoffs as `DoFlatIntegration`, after replacing `q` by the shifted `\nu_{\mathrm{eff}} / R_K` and including the extra `\chi / S_K(\chi)` factor.
- This branch is only used when the shared near-flat runtime criteria are active in both `results.f90` and `cmbmain.f90`:
  `|State%scale - 1| <= 0.03`, `chi_max < near_flat_approx_chi_limit / sqrt(AccuracyBoost * NonFlatIntAccuracyBoost)`,
  `chiDisp < near_flat_approx_chidisp_limit`, and the shared local error estimate
  `l(l+1) chi_max^2 / (15 \nu_{\mathrm{eff}}^2) < 1.0e-3 / (AccuracyBoost * NonFlatIntAccuracyBoost)`.
- The near-flat scalar approximations are also guarded by the runtime switches `enable_do_near_flat_integration` and `enable_shifted_nu_scalar_approx`, which let the full-path and local shifted-ν branches be compared without rebuilding.
- In the same regime, `lSamples_init` keeps `Ascale = 1 / lSampleBoost` so the `l` sampling stays on the flat template, and the flat Bessel spline tables are initialized for the non-flat scalar integration as well.
- The flat Bessel precomputation now accepts the requested `k eta` range directly, and the near-flat non-flat caller uses a fixed analytic extra margin so all models that use the shifted-ν approximation with the same `lmax` and `max_eta_k` request the same cached table extent. The bound uses only `near_flat_scale_tol`, `near_flat_approx_chi_limit`, `lmax`, and `max_eta_k`; for `lmax = 4000`, `lens_potential_accuracy = 4`, and `max_eta_k = 72000`, this requests `x_max = 74249` for every approximation-eligible model. In the tested cases, this safely covers the observed worst shifted-ν range, which was about `72243.6` for `Omega_k = -0.002`.

## Code simplifications

- Removed the non-flat source second-derivative storage `ddSource_q`, which is no longer needed now that both scalar and tensor paths integrate on the source grid.
- Simplified the scalar inner accumulation loop by applying trapezoidal endpoint weights to the stored `u_{\nu l}` array before accumulation, and unrolled the common 3-source case.

## Accuracy and performance checks

- Fresh-process timings include the near-flat Bessel precomputation. In a same-process repeat check, `Omega_k = -0.002` dropped from about `1.25 s` on the first call to `1.08 s` on the second, while `Omega_k = 0.02` changed only from about `2.16 s` to `2.08 s`, consistent with the first-call-only spline setup cost appearing only in the near-flat branch.
- In the current fixed-`near_flat_approx_chidisp_limit = 0.1` mode-comparison scan over `|Omega_k| = 10^{-4} .. 10^{-2}`, the three runtime modes order cleanly in speed as `both_on < shifted_nu_only < both_off`. On the earlier `10^{-6} .. 6 \times 10^{-3}` scan used to establish the workflow, mean runtimes were about `1.04 s`, `1.21 s`, and `2.03 s` respectively.
- A direct A/B probe of the shared local near-flat gate with `err_tol = 2.0e-3 / boost_scale` in both `UseNearFlatScalarIntegration` and `UseShiftedNuScalarApprox` produced spectra identical to the current `1.0e-3 / boost_scale` setting over that comparison grid, with no measurable runtime improvement beyond run-to-run noise, so the baseline threshold was kept.
- Over the near-flat cases where the shifted-ν branch is active (`Omega_k = -0.002, -10^{-5}, +10^{-5}, +0.002`), fresh-process runtime ratios versus the current baseline non-flat implementation were about `0.71`, `0.81`, `0.67`, and `0.64` respectively. Over the same cases, max differences versus that baseline for `2 <= l <= 4000` were `TT \lesssim 1.57e-3`, `EE \lesssim 2.04e-3`, `\Delta TE / \sqrt{TT \cdot EE} \lesssim 1.37e-3`, and `PP \lesssim 1.42e-3`.
- For `|Omega_k| \ge 0.02`, the shifted-ν branch is disabled by the shared `State%scale` gate, the spectra are identical to the baseline default code in the tested runs, and timing stays within about `10%` of the baseline over `Omega_k = \pm 0.02, \pm 0.1`.
- Comparing the current non-flat results to the flat `Omega_k = 0` case gives a direct continuity check near the flat limit. For `|Omega_k| = 10^{-5}`, the max differences for `2 <= l <= 4000` remain about `TT = 2.0e-4 .. 3.3e-4`, `EE = 3.0e-4 .. 4.3e-4`, and `\Delta TE / \sqrt{TT \cdot EE} = 2.7e-4 .. 3.7e-4`; with the flat-ordered near-flat integration, the low-`l` `PP` continuity is much better, dropping to about `3.6e-5` over `2 <= l <= 80` for both signs of `Omega_k`.

- Against the stored RK4 scalar baseline over `Omega_k = \pm 0.02, \pm 0.002`, the current scalar implementation gives lensed `TT` max fractional differences about `8.4e-4 .. 9.6e-4` and `EE` about `3.9e-4 .. 4.4e-4` for `\ell >= 2`. For `\ell >= 30`, `TT` stays about `8.4e-4 .. 9.6e-4` while `EE` drops to about `1.62e-4 .. 1.78e-4`.
- In the same scalar runs, CPU time dropped from about `25 .. 32 s` for the stored RK4 baseline to about `11.5 .. 15.1 s` for the current code, roughly a `52% .. 55%` reduction.
- The scalar high-`x` cutoff itself gives an additional CPU reduction of about `2.6% .. 12.5%` versus the uncapped source-grid Numerov code over the same `Omega_k` cases.
- Relative to that uncapped source-grid Numerov code, the scalar cutoff changes spectra mainly at very low `\ell`: for `\ell >= 2`, max differences are about `TT = 4.4e-4 .. 6.5e-4`, `TE = 5.7e-7 .. 9.2e-7`, and `PP = 6.2e-6 .. 7.1e-6`. For `\ell >= 30`, the same comparison shrinks to about `TT = 1.8e-6 .. 3.7e-6`, `EE = 6.6e-11 .. 9.1e-11`, `TE = 5.9e-9 .. 1.3e-8`, and `PP = 1.0e-9 .. 3.0e-9`.
- Direct comparison of lensing potential `PP` against the pre-change code shows the largest fractional differences at very low `\ell`, peaking around `\ell = 3` at about `3.1e-3 .. 3.5e-3` over tested `Omega_k = \pm 0.02, \pm 0.002`.
- For `\ell >= 30`, the `PP` differences versus the pre-change code are smaller, with max fractional differences about `7.4e-4 .. 8.5e-4` and RMS about `2.4e-4 .. 2.7e-4`.
- In the same `PP` runs, the current code was faster than the pre-change code by roughly `20% .. 33%` in CPU time over the tested curvature values.

## Scope

These changes are localized to the non-flat hyperspherical Bessel integration path in `cmbmain.f90` and its reuse of the flat Bessel spline tables in `bessels.f90` for near-flat scalar runs. Flat-space integration itself is unchanged.

## Current lpa=4 edge-status follow-up

The current high-accuracy follow-up now has a checked-in repro and plotting
entry point, [../../compare_nonflat_accuracy_sweep.py](../../compare_nonflat_accuracy_sweep.py),
which writes per-case `.npz`, difference plots, and a `summary.json` under
`accuracy_plots/nonflat_accuracy_sweep/<tag>/`. The current default rerun is in
`accuracy_plots/nonflat_accuracy_sweep/default_ab1_ab2_lpa4/`.

The key result from the `lmax = 4000`, `lens_potential_accuracy = 4`
`AccuracyBoost = 1` versus `2` follow-up is that the remaining edge `PP` miss is
still controlled by the source-integration path, not only by downstream
lensing support:

1. At default settings, `Omega_k = \pm 0.03` each still carry a low-`L`
  `lens PP` row at `L = 2` in addition to the dominant source-side `EE`
  failures.
2. Turning on `enable_olver_source_integration` removes that `PP` row for both
  edge signs, while leaving the leading `EE` failures essentially unchanged.
3. The runtime cost is very large: `+0.03` rises from about `7.15 / 144.82 s`
  CPU for AB1 / AB2 to `28.45 / 517.86 s`, and `-0.03` rises from about
  `8.16 / 237.93 s` to `33.36 / 812.29 s`.

So the low-`L` `PP` edge miss is sensitive to the scalar source integration
scheme, but the current `olver` path is much slower than the default Numerov
path.

A temporary follow-up tightened the low-`L` Numerov rebootstrap interval in
`cmbmain.f90` from `100` to `50` while keeping `enable_olver_source_integration`
off. That did not recover the edge `PP` improvement: at `Omega_k = -0.03`, the
`L = 2` `PP` row remained and slightly worsened (`0.007451`, score `1.490`,
versus the default `0.007444`, score `1.489`). So the present evidence does not
support a simple rebootstrap-frequency substitution for the `olver` source
integration path.
