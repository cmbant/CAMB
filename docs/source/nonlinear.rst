Non-linear models
==================================

SP(k) model notes
-----------------

The ``SPkNonLinear`` model wraps a base non-linear prescription (Halofit by
default) and applies SP(k) suppression multiplicatively to CAMB's non-linear
ratio.

The implementation follows the calibrated SP(k) range. In practice:

- outside calibrated redshift range, SP(k) is not applied for that slice
- above calibrated ``k``, suppression is evaluated at calibrated ``k_max``

When ``FeedbackLevel > 0``, CAMB prints one-time warnings for these out-of-range
cases.

SP(k) should not be combined with HMCode baryon-feedback modes (for example
``halofit_version='mead2020_feedback'``).

.. autoclass:: camb.nonlinear.NonLinearModel
   :members:

.. autoclass:: camb.nonlinear.Halofit
   :show-inheritance:
   :members:

.. autoclass:: camb.nonlinear.SPkNonLinear
   :show-inheritance:
   :members:

.. autoclass:: camb.nonlinear.ExternalNonLinearRatio
   :show-inheritance:
   :members:

.. autoclass:: camb.nonlinear.SecondOrderPK
   :show-inheritance:
   :members:
