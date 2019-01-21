Python CAMB
=============

CAMB (Code for Anisotropies in the Microwave Background), a cosmology code for calculating CMB, lensing,
galaxy count, dark-age 21cm power spectra, matter power spectra and transfer functions.
There are also general utility function for cosmological calculations. The main code is Python with numerical
calculations implemented efficiently in Python-wrapped modern Fortran.

See the `CAMB python example notebook <http://camb.readthedocs.org/en/latest/CAMBdemo.html>`_ for a quick
introductory set of examples of how to use the CAMB package.

For a standard non-editable installation use::

    pip install camb [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
To install from source, in the project source root directory use::

    python setup.py install [--user]

If you want to work on the code, you can also just install in place without copying anything using::

    pip install -e . [--user]

You will need gfortran 6 or higher installed to compile. Binary files for Windows are also provided, so these are used instead if no
gfortran installation is found on Windows machines. If you have gfortran installed, "python setup.py make" will build
the Fortran library on all systems (including Windows without directly using a Makefile).

After installation the camb python module can be loaded from your scripts using "import camb".
You can also run CAMB from the command line reading parameters from a .ini file, e.g.::

  camb inifiles/planck_2018.ini

You may need to check your python scripts directory is in your path for this to work.
Alternatively from the source package root directory (after make but without installation) use::

  python camb.py inifiles/planck_2018.ini

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

* :ref:`genindex`

