# Automatic lens-potential accuracy default

- Date: 2026-06-19
- Change: `CAMBparams.set_for_lmax(..., lens_potential_accuracy=None)` now selects
  `max(4, (lmax - 1500) / 500)` internally and uses the corresponding
  `lens_k_eta_reference * lens_potential_accuracy` floor for `max_eta_k`.
- Compatibility: pass `lens_potential_accuracy=0` to reproduce the previous
  low-accuracy default behavior. Explicit positive values keep their previous
  meaning.
- Diagnostics: `camb.check_accuracy` continues to use a fixed
  `lens_potential_accuracy=0` when it needs to call `set_for_lmax` and no lens
  potential accuracy was explicitly requested, so numerical diagnostics do not
  change the lensing physics by default.
- CosmoMC helper: `set_params_cosmomc` keeps its historical
  `lens_potential_accuracy=1` default; pass `None` there to request the new
  automatic high-accuracy rule.
- Motivation: upcoming high-precision lensing data need a higher default kmax
  than the historical Planck-era default. The rule follows the same-support
  convergence sweep summarized in `2026-06-07-lensing-lpa-sweeps-current.md`.
- Scope: the automatic rule keeps both the lensed CMB `C_l` and the lensing
  potential spectrum within `camb.check_accuracy` tolerances (checked both for
  numerical stability and against a `lens_potential_accuracy=30` reference).
- Floor rationale: in the floor-dominated range (`lmax < 3500`, where
  `(lmax - 1500) / 500 < 4`) the binding constraint is the lensing potential
  (`PP`) spectrum, not the lensed CMB. Numerical stability holds even at
  `lens_potential_accuracy=1`, and the lensed CMB converges against the
  `lens_potential_accuracy=30` reference by `lens_potential_accuracy≈2-3`, but
  keeping `PP` inside its check_accuracy tolerances needs `≈4` around
  `lmax≈2500-3000` (`≈3` near `lmax=2000`). The floor of 4 is therefore set by
  lensing-potential convergence and ensures the lensing spectrum is also good
  enough for CMB lensing reconstruction, with a small margin at lower `lmax`.
- At high `L`, the numerical lensed-spectrum residuals can be well below the
  non-linear modelling uncertainty.
- Docs/examples: the main demo now shows the automatic default for ordinary
  `set_for_lmax`/`set_params(lmax=...)` use, while convergence examples can still
  set fixed `lens_potential_accuracy` values explicitly.
