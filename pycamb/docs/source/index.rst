CAMB Python
==================================

CAMB is Code for Anisotropies in the Microwave Background, a cosmology code for calculating CMB and matter power spectra,
as well as general utility function for cosmological calculations. The main code is written in Fortran, which this
Python package encapulates for convenient python programming and plotting.
For an introduction to CAMB see the main `home page <http://camb.info/>`_ , and the
`CAMB python example notebook <http://camb.readthedocs.org/en/latest/CAMBdemo.html>`_ for a quick
introduction to how to use the CAMB Python package.

To install the CAMB python package on your computer, in the pycamb folder run::

    python setup.py install --user

You will need gfortran 4.9 or higher installed to compile. Binary files for Windows are also provided, so these are used instead if no
gfortran installation is found on Windows machines.


Main high-level modules:

.. toctree::
   :maxdepth: 2

   camb
   model
   symbolic

Other modules:

.. toctree::
   :maxdepth: 1

   bbn
   initialpower
   reionization
   recombination
   correlations
   postborn
   emission_angle

* :ref:`genindex`

