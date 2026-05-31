Non-linear models
==================================

.. autoclass:: camb.nonlinear.NonLinearModel
   :members:

.. autoclass:: camb.nonlinear.Halofit
   :show-inheritance:
   :members:

.. autoclass:: camb.nonlinear.ExternalNonLinearRatio
   :show-inheritance:
   :members:

.. autoclass:: camb.nonlinear.SecondOrderPK
   :show-inheritance:
   :members:

SP(k) model
------------

The ``SPkNonLinear`` model wraps a base non-linear prescription (Halofit by
default) and applies SP(k) baryon suppression multiplicatively.

See :doc:`spk` for calibrated validity domain, boundary behavior, and
MCMC/Cobaya usage guidance.

.. autoclass:: camb.nonlinear.SPkNonLinear
   :show-inheritance:
   :members:
