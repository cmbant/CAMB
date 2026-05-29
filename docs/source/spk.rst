.. _spk:

SP(k) baryon suppression model
==============================

The SP(k) model applies a multiplicative suppression to the non-linear matter
power spectrum to account for baryonic feedback effects.  It is calibrated
against hydrodynamical simulations and wrapped around a base non-linear
prescription (Halofit/HMCode by default).

Reference: `Salcido et al., MNRAS 523, 2247 (2023) <https://doi.org/10.1093/mnras/stad1474>`_

Calibrated validity domain
--------------------------

========  ==============================
Quantity  Range
========  ==============================
z         0 to 3
k         up to 12 h/Mpc
SO        200 or 500
========  ==============================

The baryon fraction f_b must also fall within per-(z, k, SO) calibrated limits
derived from the fitting coefficients.

Boundary behavior
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Condition
     - Action
     - Rationale
   * - z outside [0, 3]
     - Identity (no suppression)
     - Feedback sub-dominant at high z; model not calibrated there.
   * - k < k_min (~10⁻¹² h/Mpc)
     - Identity (no suppression)
     - Baryonic effects negligible on large scales.
   * - k > 12 h/Mpc
     - Clamp to k = 12
     - Suppression saturates at high k; avoids breaking P(k) integrals (σ₈, lensing).
   * - f_b outside calibrated limits
     - ``nonlin_ratio`` set to NaN
     - Model invalid — no reliable prediction possible.

When ``FeedbackLevel > 0``, a one-time warning is printed for each out-of-range
case encountered during a calculation.

MCMC / Cobaya considerations
----------------------------

The boundary choices are designed for stable sampling:

- **k-clamping** keeps P(k) integrals (σ₈, CMB lensing) finite for every
  sample. The clamp is physically motivated: SP(k) suppression saturates at
  high k.

- **z-skip** is safe because most cosmological observables (weak lensing,
  galaxy clustering, CMB lensing) probe z < 3.

- **f_b → NaN** acts as a hard prior boundary.  Any parameter combination
  pushing f_b outside the calibrated envelope produces NaN in the non-linear
  ratio, which propagates through
  :meth:`~camb.results.CAMBdata.get_matter_power_interpolator` and causes
  Cobaya to assign −∞ log-likelihood (sample rejected).

.. note::

   :meth:`~camb.results.CAMBdata.get_matter_power_spectrum` uses Fortran-side
   spline interpolation onto a regular log-k grid that does **not** preserve
   NaN (it produces tiny finite unphysical values ~10⁻²²).

   :meth:`~camb.results.CAMBdata.get_matter_power_interpolator` builds a 2D
   ``RectBivariateSpline`` over the internal (z, k) grid.  Any NaN on that
   grid propagates through the entire spline, so **all** evaluations return
   NaN — not just the invalid (z, k) cells.  This is the correct behaviour
   for MCMC: a single out-of-range point invalidates the whole sample.

   For diagnostics (locating *which* (z, k) cells are NaN), use
   :meth:`~camb.results.CAMBdata.get_linear_matter_power_spectrum` with
   ``nonlinear=True``, which returns the raw nonlinear P(k) on the internal
   transfer-function k-grid without any interpolation.  When comparing against
   other spectra, interpolate *those* onto this k-grid (not vice versa) to
   avoid spurious oscillations from resampling through BAO features.

**Prior guidance:** set priors on the SP(k) relation parameters to keep f_b
within calibrated limits for the bulk of your parameter space.  The NaN
boundary then acts only as a safety net for rare excursions, not a dominant
rejection mechanism.

Base non-linear model configuration
-------------------------------------

``SPkNonLinear`` wraps a ``Halofit`` instance as its base non-linear model.
The full set of Halofit/HMCode parameters — ``halofit_version``,
``HMCode_A_baryon``, ``HMCode_eta_baryon``, and ``HMCode_logT_AGN`` — can be
passed directly to :meth:`~camb.nonlinear.SPkNonLinear.set_params` and are
forwarded to the wrapped base model::

    spk = camb.SPkNonLinear()
    spk.set_params(
        halofit_version="mead2020_feedback",
        HMCode_logT_AGN=8.0,
        SPk_feedback=False,
    )

These parameters are also accepted as flat keyword arguments by
``camb.set_params()``, and as ``extra_args`` in Cobaya.

.. warning::

   ``SPk_feedback=True`` cannot be combined with
   ``halofit_version='mead2020_feedback'`` or non-default
   ``HMCode_A_baryon``/``HMCode_eta_baryon`` values (HMCode 2015/2016), as
   this would double-count baryonic corrections.  A ``CAMBValueError`` is
   raised if this combination is detected.

API reference
-------------

See :class:`camb.nonlinear.SPkNonLinear` for the full parameter list and
Cobaya YAML example.
