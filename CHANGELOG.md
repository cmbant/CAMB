# CAMB v2 Change Summary

This is a high-level summary of important user-facing and result-facing changes
on the development branch that will become CAMB v2. It is not a complete commit
log; see the GitHub history and the detailed notes under `docs/changelog/` for
all tuning, validation, and implementation details.

## Accuracy Defaults And Diagnostics

- CAMB now enables targeted accuracy improvements by default with the internal
  `AccuracyTarget = 1`. These changes are aimed at modern high-precision CMB,
  lensing, and matter-power use cases. Set `AccuracyTarget = 0` in an `.ini`
  file, or `camb.config.AccuracyTarget = 0` from Python, for behavior closer to
  the CAMB 1.x numerical-error profile.
- The new `camb.check_accuracy` module and `camb check_accuracy` command compare
  a requested calculation to a higher-accuracy reference, report CMB, lensing
  potential, matter-power, and derived-parameter differences, optionally make
  plots, and can search for minimal accuracy boosts.
- `CAMBparams.set_for_lmax(..., lens_potential_accuracy=None)` is now the public
  default. `None` selects an automatic high-accuracy lensing-potential/kmax
  setting,
  `max(4, (lmax - 1500) / 500)`.
- Explicit `lens_potential_accuracy` values keep their old meaning. In
  particular, use `lens_potential_accuracy=0` to reproduce the old low-k default
  behavior. `set_params_cosmomc` keeps its historical default
  `lens_potential_accuracy=1`; pass `None` there to opt into the new automatic
  rule.
- The automatic lens-potential rule is calibrated for lensed CMB spectra and
  lensing-potential stability at the relevant accuracy target. At high
  multipoles, remaining numerical errors in lensed spectra can be much smaller
  than the uncertainty from non-linear matter modelling.
- `lens_output_margin` is now a first-class Python and `.ini` parameter. It
  consistently controls how far above the requested lensed output range CAMB
  calculates internally, including the Fortran lensing convolution support.

## Non-Flat Models And Hyperspherical Bessel Functions

- Non-flat scalar line-of-sight integration has been substantially refactored
  and sped up. The main changes are Numerov/source-grid integration, improved
  high-oscillation cutoffs, near-flat shifted-ν approximations, and direct Olver
  evaluation in high-substep ranges.
- Near-flat open and closed models can reuse flat Bessel table machinery where
  controlled local error estimates allow it. This improves speed near the flat
  limit while preserving continuity checks.
- The branch includes new hyperspherical Bessel implementations and validation
  paths: Olver-style approximations, small-chi/open-small-nu fallbacks, Airy
  utilities, and Python-accessible math utilities for testing.
- In the documented branch comparison, default non-flat/matter-power cases were
  about three times faster than the CAMB 1.x baseline in that test set. Focused
  non-flat scalar integration timings show larger speedups in some curvature
  regimes, with low-level numerical changes documented in `docs/changelog/`.

## Lensing Calculations

- CAMB has a new optimized lensing method selector. It keeps the long-standing
  fast curved-sky method for ordinary runs and uses a full Gauss-Legendre
  curved-sky correlation method when `AccurateBB=True`.
- The direct curved-sky lensing implementation can also be selected explicitly,
  and Python calls such as `get_lensed_cls_with_spectrum` can temporarily
  override the lensing method for comparisons.
- The full-sky correlation code was optimized substantially, including cached
  Gauss-Legendre nodes/weights, inlined accumulation, recurrence-based factors,
  and faster Legendre tables exposed through `camb.mathutils`.
- Low-l EE tapering and high-L template extension behavior have been clarified
  and made more consistent between the Python and Fortran lensing paths.

## Matter Power And Non-Linear Modelling

- Matter-power accuracy tuning was updated for massive neutrinos, photon and
  massless-neutrino hierarchy depths, and transfer-high-precision cases. The
  goal is better default agreement with boosted references without requiring
  broad global accuracy boosts.
- HMCode/Halofit evaluation was cleaned up and optimized. Cached HMCode
  redshift-local quantities give speedups of order 10-20% in the documented
  matter-power benchmarks, with only tiny changes from removing unintended
  single-precision round trips.
- CAMB now includes an `SPkNonLinear` model for the SP(k) baryon-suppression
  prescription, wrapping a base Halofit/HMCode model. It includes documented
  validity ranges, MCMC-friendly boundary behavior, and protections against
  double-counting baryonic feedback.
- New non-linear model hooks include `ExternalNonLinearRatio` for externally
  supplied non-linear ratios and `SecondOrderPK` for second-order perturbative
  matter-power ratios.

## Recombination, Reionization, And Backgrounds

- The default BBN consistency relation now uses the September 2024 PRIMAT
  helium and deuterium table, replacing the 2021 PRIMAT table. For typical
  Planck-like models this lowers the default helium mass fraction by about
  2e-4, with sub-per-mille effects on fixed-parameter CMB spectra.
- RECFAST now uses a fast Rosenbrock-to-RK45 handoff mode by default. The new
  path is tuned against high-accuracy internal references and scales with CAMB
  accuracy boosts. It is intended to improve the speed/accuracy tradeoff of the
  recombination background calculation.
- The default RECFAST approximation is now the `recfast_cosmorec` fit, including
  the helium-rate correction calibrated against direct CosmoRec histories.
  Planck-era RECFAST parameters remain available as `recfast_planck` and are
  explicitly used by Planck-specific compatibility inputs.
- The CosmoRec wrapper was updated for the newer CosmoRec vX interface and
  exposes the relevant CosmoRec controls through CAMB's recombination model.
- Reionization models now have an optional approximate heating switch that
  raises the baryon temperature and sound speed during reionization. It is off
  by default and is intended for order-of-magnitude low-redshift matter-power
  effects rather than precision thermal-history modelling.
- Thermal massive-neutrino background density and pressure now use direct smooth
  fits over the intermediate mass range, reducing setup/global state and modestly
  speeding repeated background evaluations.

## Python Interface And New Capabilities

- A Python bispectrum wrapper is now available as `camb.bispectrum`. It runs the
  existing Fortran CMB-lensing or local-primordial bispectrum calculation using
  normal `CAMBparams` objects, writes large tables directly to files, and
  returns small Fisher summaries when the library is built with Fisher support.
- Documentation now includes pages for the bispectrum wrapper, SP(k), nonlinear
  models, check-accuracy workflow, and math utilities.
- CAMB now targets Python 3.11+ and uses the ruff/pre-commit toolchain for
  Python formatting and linting.
- The development tree includes updated devcontainer and CI configuration, but
  those changes are primarily for contributors rather than result-facing users.

## Compatibility Notes

- Numerical outputs can change relative to CAMB 1.x because the v2 branch has a
  higher default accuracy target, different non-flat algorithms, updated
  lensing support, tuned matter-power accuracy settings, and RECFAST changes.
- For closer 1.x-style numerical behavior, start with `AccuracyTarget = 0`,
  explicit `lens_potential_accuracy=0` in `set_for_lmax`, and fixed legacy
  matter-power settings where comparing against older runs.
- Users comparing old and new results should use `camb check_accuracy` and
  compare at fixed physical output ranges and k ranges. Avoid judging changes
  only from sparse grid-index differences, especially for matter power and
  high-l lensing.
- Some new options are deliberately off by default because they change the
  physical model rather than only the numerical method, for example reionization
  heating.
