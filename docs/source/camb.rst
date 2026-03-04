Basic functions
==================================


.. automodule:: camb
   :members: get_results, get_background, get_transfer_functions, set_params, set_params_cosmomc, read_ini, get_matter_power_interpolator, get_age, get_zre_from_tau, set_feedback_level, run_ini, get_valid_numerical_params


Type aliases
------------

.. py:class:: camb.Array1D

   Type alias for 1D array-like inputs: ``Sequence[np.number | float | int] | NDArray[np.number]``.
   Functions accepting Array1D can take lists, tuples, or numpy arrays of numbers.
