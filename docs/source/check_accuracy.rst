.. _check-accuracy:

Accuracy stability checks
=========================

When developing new models, adding new approximations, changing numerical
methods, or updating the Fortran calculation, it is good practice to check that
the requested results are stable to tighter numerical accuracy settings. Small
changes in the code can otherwise move spectra in ways that are hard to spot
from a single run. The ``check_accuracy`` command compares a normal CAMB run
with a higher-accuracy reference run and reports where any differences exceed
configurable tolerances.

The same implementation is available from Python as the
:mod:`camb.check_accuracy` module, and from the command line as::

    camb check_accuracy inifiles/planck_2018.ini

This tests numerical stability with mostly fixed scale cuts. Also check
lens-potential-accuracy is set high enough (effectively increasing kmax), though
note that inaccuracies in the default non-linear halo model are often larger than
numerical errors.

What it compares
----------------

The checker loads an input ``.ini`` file using :func:`camb.read_ini`, runs it
once as requested, then runs a reference calculation with boosted accuracy
settings. By default the reference uses::

    AccuracyBoost = 2
    lSampleBoost = 2
    lAccuracyBoost = 2
    IntTolBoost = 2
    DoLateRadTruncation = True

It compares:

* derived parameters, reported as fractional changes;
* total lensed CMB spectra from ``get_total_cls()``, so tensor contributions are
  included when present;
* the lensing potential spectrum from ``get_lens_potential_cls()``;
* matter power spectra when transfer functions are requested by the input
  parameters.

For TT and EE the comparison is fractional in each multipole range. For TE it
uses ``Delta TE / sqrt(TT * EE)``, which avoids artificial problems at TE
zero-crossings. The reported tables include maximum and rms errors, the
tolerance, pass/fail status, and the multipole or grid point of the worst
sample. CPU and wall times are reported for the standard run, reference run, and
any candidate runs tested during boost searches.

Common uses
-----------

Basic stability check for a standard parameter file::

    camb check_accuracy inifiles/planck_2018.ini

Check spectra only up to a realistic analysis scale, while forcing CAMB to
calculate that far first::

    camb check_accuracy inifiles/planck_2018.ini --set-for-lmax 4000

Use higher lensing-potential accuracy, useful for high-accuracy lensed spectra::

    camb check_accuracy inifiles/planck_2018.ini --set-for-lmax 4000 --lens-potential-accuracy 4

Make diagnostic plots of the fractional differences::

    camb check_accuracy inifiles/planck_2018.ini --plot-dir accuracy_plots

Calculate a fiducial CMB delta chi-squared using a Simons Observatory-like
noise model::

    camb check_accuracy inifiles/planck_2018.ini --set-for-lmax 4000 --chi2 --chi2-config so

Search for the smallest top-level boost settings that pass against the
reference calculation::

    camb check_accuracy inifiles/planck_2018.ini --find-minimal-boosts

Then try to identify which underlying component accuracy parameters are most
important, keeping ``AccuracyBoost=1``::

    camb check_accuracy inifiles/planck_2018.ini --find-minimal-boosts --refine-accuracy-components

This component refinement varies ``IntTolBoost`` and the lower-level component
accuracy settings affected by ``AccuracyBoost``. It does not vary
``lSampleBoost`` or ``lAccuracyBoost``; those remain top-level boost parameters
handled by ``--find-minimal-boosts``.

Assess the physical effect of changing lensing reference settings without
changing the comparison run::

    camb check_accuracy inifiles/planck_2018.ini --reference-lens-potential-accuracy 8

Options
-------

The command-line option list below is generated from the current parser, so it
stays aligned with ``camb check_accuracy --help``.

.. program-output:: camb check_accuracy --help

Programmatic use
----------------

For scripts and tests that already have a :class:`camb.model.CAMBparams` object,
use :func:`camb.check_accuracy.compare_params_accuracy` directly rather than
writing a temporary ``.ini`` file. The returned
:class:`camb.check_accuracy.AccuracyCheckResult` contains the standard and
reference run outputs, comparison tables, timings, and optional chi-squared
summary.

For example::

    import camb
    from camb.check_accuracy import compare_params_accuracy

    params = camb.read_ini("inifiles/planck_2018.ini")
    result = compare_params_accuracy(params, set_for_lmax=4000)

    if not result.comparison.passed:
        print(result.comparison.worst_failure)

API reference
-------------

.. automodule:: camb.check_accuracy

.. autofunction:: camb.check_accuracy.compare_params_accuracy

.. autoclass:: camb.check_accuracy.AccuracyCheckResult
   :members:

.. autoclass:: camb.check_accuracy.SearchResult
   :members:

.. autoclass:: camb.check_accuracy.ComparisonResult
   :members:

.. autoclass:: camb.check_accuracy.RunOutput
   :members:

.. autoclass:: camb.check_accuracy.NoiseConfig
   :members:
