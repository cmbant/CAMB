Bispectrum
==================================

The bispectrum wrapper runs the existing Fortran calculation from Python, using normal
:class:`~camb.model.CAMBparams` objects rather than writing temporary ``.ini`` files. It can calculate CMB lensing
and local primordial bispectra, write large slice or full bispectrum tables directly to files, and return the small
in-memory Fisher summary when requested.

The module is not imported by default with ``import camb``. Import it explicitly when needed::

   import camb
   from camb import bispectrum

   pars = camb.set_params(lmax=600, lens_potential_accuracy=1)
   bpars = bispectrum.BispectrumParams(Slice_Base_L=10, deltas=[0, 2])
   result = bispectrum.get_bispectrum(pars, bpars, output_root="run1_")

   print(bpars.expected_output_files("run1_"))
   print(result.has_fisher)

The default :class:`~camb.bispectrum.BispectrumParams` calculates the CMB lensing bispectrum, so the CAMB parameters
must have lensing enabled. To calculate the local primordial bispectrum normalized to ``f_NL=1`` instead, use
``BispectrumParams(do_lensing_bispectrum=False, do_primordial_bispectrum=True)``.

Fisher Matrices
---------------

Fisher matrix output is disabled in normal builds to avoid requiring LAPACK. To use
``BispectrumParams(DoFisher=True)``, build the Fortran library from source with ``FISHER=Y`` and LAPACK/BLAS linked.
For the default gfortran makefile this is, for example::

   cd fortran
   make python FISHER=Y

The gfortran makefile uses ``-lblas -llapack`` by default for Fisher builds. Intel ``ifort`` builds use MKL when
``FISHER=Y`` is set.

.. automodule:: camb.bispectrum
   :members:
