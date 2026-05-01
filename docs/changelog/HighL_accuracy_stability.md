# Planck 2018 Accuracy Stability

## Summary

This investigation targeted

```bash
camb check_accuracy inifiles/planck_2018.ini --set-for-lmax 4000 --lens-potential-accuracy 4
```

with separate treatment of CMB power spectra and matter power. The final changes
avoid broad accuracy-boost floors and instead adjust the specific source-grid,
lensing-correlation, and late photon-hierarchy ranges that control the observed
instability.

## Main Conclusions

- The low-ell EE instability is controlled by low-`k` CMB source sampling, not by
  high-`L` lensing physics. A direct cap on the low-`k` logarithmic source step is
  preferable to increasing a global integration or source boost.
- The high-`L` lensed CMB instability at `lmax=4000` is controlled by
  correlation-function angular sampling. A separate `ThetaSampleBoost` floor is
  cleaner than changing `LensAccuracyBoost`, which also affects unrelated
  lensing work.
- The matter-power failure is not fixed by increasing output
  `transfer_k_per_logint`; sweeps up to very dense output sampling left the same
  interpolation-stable mismatch. The sensitive calculation is the shared CMB
  source `q` grid used by high-precision transfer output, plus the late photon
  hierarchy truncation time.
- The high-precision transfer source-grid spacing should follow an explicit
  dense `transfer_k_per_logint` request. For automatic sampling
  (`transfer_k_per_logint = 0`), 20 points per log interval is sufficient for the
  low-`q` source grid in this case.
- Extra refinement of the smooth high-`q` source tail via `dksmooth` was tested
  and found unnecessary after the low/intermediate source-grid fixes.
- At `lmax=6000` and `8000`, larger theta sampling did not resolve the remaining
  CMB differences, even with higher lensing accuracy. Those failures appear to be
  a separate high-`lmax` issue and are not addressed by this change set.

## Code Changes

- `fortran/cmbmain.f90`
  - Caps the low-`k` CMB source logarithmic step for accurate reionization
    polarization.
  - Uses a combined high-precision transfer/CMB-source condition for the source
    grid reused by transfer output.
  - Sets the low-`q` source log density to
    `max(20, transfer_k_per_logint)` for high-precision transfer output.
  - Refines the two linear source intervals, `dkn1` and `dkn2`, by a factor of
    `1.5` for this shared high-precision transfer source grid.
  - Caps `qmax_log` at the horizon-entry switch to avoid an unnecessary
    intermediate interval when the logarithmic range already reaches that scale.
  - Increases the low-`k` scalar integration grid minimum from 10 to 11 log
    samples for accurate reionization polarization.

- `fortran/equations.f90`
  - Isolates the high-precision transfer switch change to late photon multipole
    truncation, `tau_switch_no_phot_multpoles`.
  - Leaves late neutrino multipole truncation and massive-neutrino switch timing
    on the original `TimeSwitchBoost` behavior.

- `fortran/lensing.f90`
  - Adds `ThetaSampleBoost` as a separate angular-sampling factor for lensed CMB
    correlation functions.
  - Applies a floor of `1.8` only for `CP%Max_l > 3500`, avoiding a broad
    `LensAccuracyBoost` increase.

- `camb/tests/camb_test.py`
  - Updates the nonlinear matter-power regression check to compare at a fixed
    physical `k/h` by interpolation, rather than checking the grid-fragile
    `pk[0][160]` index.

## Tests

- Original target command passes with automatic transfer sampling:
  - matter power max error: `0.0008865`
  - standard run timing: about `8.33s` CPU, `1.08s` wall in the final check run
- Explicit dense transfer sampling also passes:
  - command uses `transfer_k_per_logint = 50`
  - matter power max error: `0.0007711`
  - standard run timing: about `11.59s` CPU, `1.50s` wall
- `python -m unittest camb.tests.camb_test` passes: 20 tests.
- `git diff --check` passes for the touched source and test files.

## Negative Tests

- Removing the `dkn2` refinement fails with matter-power max error `0.001343`
  near `k/h = 0.104`.
- Removing the `dkn1` refinement fails with matter-power max error `0.001279`
  near `k/h = 0.076`.
- Refining both `dkn1` and `dkn2` by only `1.25` passes but with too little
  margin: matter-power max error `0.0009885`.
- Refining `dksmooth` by `4` passes but is not needed and is slower; it was
  removed from the final patch.
- Neutrino-only late radiation truncation was not sufficient for the matter
  power stability target; photon-only late truncation was sufficient.
