# lensing-template-tail-acoustic-shift: shift the high-L template tail by the peak-scale ratio

- Date / issue id: 2026-11-27. Remaining high-`L` lensed TE/EE `check_accuracy` failures in the generated
  test set (the "high-`L` EE/TE" group in `matter_power_transfer_accuracy_review.md`).
- Symptom & minimal reproducer: on master (`3d524fe`), with inis from
  `python fortran/tests/CAMB_test_files.py <dir> --make_ini --no_run_test`,
  `python -m camb.check_accuracy <dir>/params_omk_0.040.ini` fails `TE 600 <= L <= 2800` with
  `2.379e-3 > 1e-3` at `L = 2789`. Same pattern for `params_omch2_0.080.ini` (TE `2.065e-3` at 2788, EE
  `1.019e-3` at 2777), `params_omch2_0.150.ini` (TE `1.308e-3`) and `params_omk_-0.030.ini` (TE `1.697e-3`).
  The other inis in that group (`hubble_62/67/78`, `ombh2_0.025`, `omch2_0.100`, `w_-0.750/-1.200`) already pass
  on master, and the `omk_-0.030` interpolation-range error no longer occurs.
- Root cause: drill-down `AccuracyBoost` (only broad axis that restores it; `lAccuracyBoost=2` and
  `lSampleBoost=2` do not) -> `LensingBoost` -> `CorrFuncFullSkyApodized` -> `range_fac` (short angular range of the
  correlation-function integral). Unlensed spectra are converged near the top (`TE <= 2.3e-4`); the error appears
  only in the lensing step, oscillates in `L` with period ~`2*pi/theta_max`, and with identical unlensed input the
  full-range Gauss-Legendre method (`AccurateBB`) converges while the short-range one does not (`range_fac` 32, 16,
  8, 4 give TE `2.3e-3, 1.1e-3, 1.2e-3, 1.8e-3`; only the full range gives `2e-4`). The ringing source is the input:
  `PrepareLensedCLSpectra` (`fortran/lensing.f90`) extends TT/EE/TE above `CP%Max_l` with the Planck 2018 template,
  matching only the amplitude at `CP%Max_l`. For a model whose acoustic scale differs from the template, the TE tail
  is out of phase and jumps at `CP%Max_l` by `0.33, 0.28, 0.19, 0.24` of `sqrt(TT*EE)` for the four failing
  models (`0.01-0.03` for the passing ones). The short-range correlation integral smears that step into the top
  ~300 output multipoles. TT/EE are continuous at `CP%Max_l` but also out of phase above it (10-40% errors at
  `l = 3100-4500` against computed spectra).
- Controlling quantity: the multipole at which the template is evaluated in the tail,
  `highL_CL_template(l, ...)` in `PrepareLensedCLSpectra`.
- Fix: evaluate the TT, EE and TE template at `l/State%scale` (linear interpolation), where `State%scale` is the
  existing peak-scale ratio relative to the Planck 2018 drag angle, already used to shift the same template for the
  `use_cl_spline_template` Cl interpolation (`InterpolateClArrTemplated`). The lensing potential tail has no acoustic
  phase and is not shifted. The Python mirror `correlations.lensed_cls(..., use_lensing_template=True)` gets a
  `template_l_scale` argument (default 1, the old behaviour) so it can reproduce the Fortran convention with
  `CAMBdata.scale`.
- Full physical extent: all lensed CMB spectra that use the template tail, i.e. every Fortran lensing method
  (fast short-range, full-range `AccurateBB`, flat-sky, harmonic) through the shared `PrepareLensedCLSpectra`, and only
  output `L` within a few hundred of `lmax - lens_output_margin`. Unlensed spectra, lensing potential, matter power
  and derived parameters do not use this code. For the Planck 2018 fiducial `State%scale = 1` to `~1e-5` and results
  are unchanged.
- Locality / rejected alternatives:
  - `LensingBoost` / `LensRangeBoost = 2` passes `check_accuracy`, but only because the boosted reference uses the
    same setting; against the full-range reference it leaves TE at `1.1e-3`, and larger ranges are not monotone.
  - Wider apodization (`x2`, `x3`) and finer `interp_fac` or theta sampling: TE stays `1.5-2.8e-3`.
  - Smoothly blending the TE step away over 25-400 multipoles fixes the method ringing (`2.4e-4` at width >= 100) but
    biases the physical answer: against a run with unlensed spectra computed to `Max_l = 4500` TE gets worse
    (`8.0e-4 -> 1.9e-3` for `omk_0.040`, `8.7e-5 -> 1.8e-3` for `omch2_0.080`). Rejected.
  - Shifting only TE fixes the TE failure but leaves the EE tail phase error (`4.1e-3` lensed EE error near the top
    against `Max_l = 4500`), so TT/EE are shifted too.
  - No extra sampling or range, so there is no runtime trade-off to tune.
- Change: `fortran/lensing.f90` (`ShiftedHighLTemplate`, `PrepareLensedCLSpectra`), `camb/correlations.py`
  (`_extrapolate_lensing_inputs`, `lensed_cls`), `camb/tests/camb_test.py` (`testPowers` passes
  `template_l_scale=data.scale`; new `testLensingTemplateTailPhase`).
- Before/after (`check_accuracy`, standard vs default boosted reference, `600 <= L <= 2800`):

  | ini | TE before | TE after | EE before | EE after |
  | --- | ---: | ---: | ---: | ---: |
  | `omk_0.040` | `2.379e-3` | `2.18e-4` | `9.14e-4` | `5.04e-4` |
  | `omch2_0.080` | `2.065e-3` | `3.74e-4` | `1.019e-3` | `8.51e-4` |
  | `omch2_0.150` | `1.308e-3` | `3.14e-4` | `7.67e-4` | `7.67e-4` |
  | `omk_-0.030` | `1.697e-3` | `2.79e-4` | `6.22e-4` | `5.51e-4` |
  | `nu_masssplit` | `1.797e-3` | `3.02e-4` | `8.69e-4` | `7.23e-4` |
  | `omch2_0.100` | `5.75e-4` | `2.15e-4` | `7.63e-4` | `4.64e-4` |
  | `w_-0.750` | `8.45e-4` | `2.01e-4` | `5.98e-4` | `4.50e-4` |
  | `w_-1.200` | `6.12e-4` | `1.91e-4` | `6.99e-4` | `4.50e-4` |
  | `ombh2_0.025` | `2.65e-4` | `2.59e-4` | `6.01e-4` | `6.00e-4` |

  `params_nu_masssplit.ini` also failed TE on master (`1.797e-3`) and passes after; it was not in the review's list.

- Reference independence (same unlensed input, lensing-only reference `AccurateBB=T, LensingBoost=4`, standard fast
  method): `omk_0.040` TE `2.28e-3 -> 1.29e-4`, EE `8.0e-4 -> 3.1e-4`; `omch2_0.080` TE `1.88e-3 -> 1.31e-4`,
  EE `7.6e-4 -> 2.2e-4`.
- Physical accuracy of the tail (python `lensed_cls`, full range, identical numerics; truth uses unlensed spectra
  computed to `l = 5000` with `lSampleBoost = 2`; max over `600 <= L <= 2800`):

  | model | EE unshifted | EE shifted | TE unshifted | TE shifted | TT unshifted | TT shifted |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: |
  | `omk_0.040` | `4.0e-3` | `1.7e-4` | `8.2e-4` | `5.6e-5` | `1.4e-4` | `2.1e-4` |
  | `omk_-0.030` | `4.7e-3` | `9.8e-5` | `3.4e-4` | `1.1e-4` | `8.9e-4` | `4.7e-4` |
  | `omch2_0.150` | `7.0e-3` | `8.9e-5` | `2.7e-4` | `1.5e-4` | `1.5e-3` | `1.6e-4` |
  | `hubble_67` | `3.6e-3` | `5.9e-5` | `4.0e-4` | `7.0e-5` | `9.5e-4` | `2.3e-4` |
  | `omch2_0.080` | `8.5e-4` | `6.1e-4` | `9.7e-5` | `1.1e-4` | `9.6e-4` | `3.0e-4` |

  Directly, the shifted TT/EE tails match computed unlensed spectra at `l = 3100-4500` to 0.1-3% (unshifted: 4-40%).
  So the change also removes a lensed EE error near the top (`~4e-3` for typical models) that `check_accuracy`
  cannot see, because standard and boosted runs share the tail.
- Residual: lensed BB near the top is still limited by the (unshifted) lensing-potential template tail. With the
  CMB tails shifted, BB against the truth above is `+2e-3` at `L = 2100` and `+6e-3` at `L = 2500` (tolerance `2e-2`);
  using the true `C_L^{phi phi}` tail as well removes it (`<= 8e-4`). Before this change the unshifted EE phase error
  partly cancelled that bias in some models and doubled it in others (`hubble_67`: `+9.9e-3` at `L = 2500`). The
  PP tail is a separate issue and is not changed here. The unlensed baseline EE excess noted in the review
  (`~1.1e-3` in `600 <= L < 3000`, e.g. unlensed EE `1.37e-3` at `L = 2248` for `omch2_0.080`) does not involve
  lensing and is also not addressed.
- Acceptance matrix (main tree rebuilt with `python setup.py make`; "before" is a helper-made worktree of `3d524fe`):
  - `python -m unittest camb.tests.camb_test`: all 26 tests pass. `testPowers` (Python `lensed_cls` mirror with
    `template_l_scale=data.scale`) still agrees with Fortran to `1e-6`; the Planck 2018 high-`L` check is unchanged.
  - New `testLensingTemplateTailPhase`: fails on master (max TE `2.9e-3`), passes with the fix (`1.8e-4`, limit `5e-4`).
  - `pre-commit run --all-files`: clean.
  - `check_accuracy` before/after on 63 of the 136 generated inis (a sample, not the full set): the first 44
    alphabetically through `params_l_accuracy_boost2.ini` (`params_21cm_base2.ini` skipped for memory), plus 19
    selected for high-`L` CMB/lensing coverage (`lmax3000/6000`, `nearflat_omk_m0p03/p0p03`, `nu_masssplit`,
    `ombh2_0.025`, `omch2_0.080/0.100/0.150`, `omk_-0.001/-0.030/0.040`, both `ppf_w-0.7_wa0.5`/`ppf_w-1.2_wa-0.3`,
    `re_optical_depth_0.110`, `thetastar_lmax4000_lpa_auto`, `three_deg_omnuh2_0.03000`, `w_-0.750/-1.200`).
    Five go FAIL -> PASS (`nu_masssplit`, `omch2_0.080`, `omch2_0.150`, `omk_-0.030`, `omk_0.040`), the other 58 pass
    on both builds, and no row in any of the 63 grows by more than 10% where it is above `1e-4`.
- Runtime impact: no measurable change. Default `inifiles/params.ini` (lensing on, `l_max_scalar = 3000`, `State%scale = 1.010`), 24-thread CPU, 4 interleaved processes per build x 6 calls, first call of each discarded (20 warm samples each): `get_results` + lensed Cls median `81.1 ms` master vs `82.1 ms` fix (min `73.3` vs `72.4`); lensing step alone (`get_lensed_cls_with_spectrum`) median `2.38 ms` vs `2.26 ms`. The difference is within run-to-run noise; the change only adds a linear interpolation per tail multipole.
- Physical justification: the damping-tail TT/EE/TE of a model are, to percent level, the fiducial spectra with
  multipoles rescaled by the ratio of acoustic angular scales; `State%scale` is that ratio and is already used for
  the same template elsewhere. Matching the phase removes the spurious TE step and makes the extrapolated tail an
  actual approximation to the model.
- Checker verdict: no separate checker agent was run. Checked directly: the new test fails on master and passes with the fix; the TE/EE gains hold against an independent full-range reference (`AccurateBB=T`, `LensingBoost=4`) and against unlensed spectra computed to high `l`, so they are not an artefact of the boosted `check_accuracy` reference; the fiducial model is unchanged (`State%scale = 1` to `~1e-5`, Planck 2018 high-`L` check unchanged); the `check_accuracy` sample above shows no regressions; residuals (BB tail, unlensed EE) are stated, not hidden.
