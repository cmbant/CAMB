CAMB
=============

CAMB (Code for Anisotropies in the Microwave Background), a cosmology code for calculating CMB, lensing,
galaxy count, dark-age 21cm power spectra, matter power spectra and transfer functions.
There are also general utility function for cosmological calculations such as the background expansion, distances, etc.
The main code is Python with numerical calculations implemented efficiently in Python-wrapped modern Fortran.

See the `CAMB python example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_ for an
introductory set of examples of how to use the CAMB package. This is usually the fastest way to learn how to use it
and quickly see some of the capabilities.

For a standard non-editable installation use::

    pip install camb [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
If you want to work on the code from `GitHub <https://github.com/cmbant/camb>`_, you can also just install in place without copying anything using::

    git clone --recursive https://github.com/cmbant/CAMB.git
    pip install -e ./CAMB [--user]

You will need ifort or gfortran 6 or higher installed (and on your path) to compile; see :ref:`fortran-compilers` for
compiler installation details if needed. A compiled library for Windows is also provided, and is used if no
gfortran installation is found on Windows machines. If you have gfortran installed, "python setup.py make" will build
the Fortran library on all systems (including Windows without directly using a Makefile), and can be used to update
a source installation after changes or pulling an updated version.

Anaconda users can also install from conda-forge, best making a new clean environment using::

  conda create -n camb -c conda-forge python=3.9 camb
  activate camb

with no need for a Fortran compiler (unless you want to use custom sources/symbolic compilation features).
Check that conda installs the latest version, if not try installing
in a new clean conda environment as above.

After installation the camb python module can be loaded from your scripts using "import camb".

You can also run CAMB from the command line reading parameters from a .ini file, e.g.::

  camb inifiles/planck_2018.ini


Sample .ini files can be obtained from the `repository <https://github.com/cmbant/CAMB/tree/master/inifiles>`_. You can load parameters programmatically from an .ini file or URL using :func:`.camb.read_ini`.

Main high-level modules:

.. toctree::
   :maxdepth: 2

   camb
   model
   results
   symbolic

Other modules:

.. toctree::
   :maxdepth: 1

   bbn
   dark_energy
   initialpower
   nonlinear
   reionization
   recombination
   sources
   correlations
   postborn
   emission_angle

.. toctree::
   :maxdepth: 1

   transfer_variables
   fortran_compilers
   mathutils

* `Example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_
* :ref:`genindex`

===================

.. image:: https://cdn.cosmologist.info/antony/Sussex_white.svg
   :alt: University of Sussex
   :target: https://www.sussex.ac.uk/astronomy/
   :height: 170px
   :width: 170px

.. image:: https://cdn.cosmologist.info/antony/ERC_white.svg
   :alt: European Research Council
   :target: https://erc.europa.eu/
   :height: 170px
   :width: 170px

.. image:: https://cdn.cosmologist.info/antony/STFC_white.svg
   :alt: Science and Technology Facilities Council
   :target: https://stfc.ukri.org/
   :height: 170px
   :width: 170px
