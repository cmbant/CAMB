===================
CAMB
===================
:CAMB: Code for Anisotropies in the Microwave Background
:Author: Antony Lewis and Anthony Challinor
:Homepage: https://camb.info/

.. image:: http://img.shields.io/pypi/v/camb.svg?style=flat
        :target: https://pypi.python.org/pypi/camb/
.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.org/en/latest


Description and installation
=============================

CAMB is a cosmology code for calculating cosmlogical observables, including
CMB, lensing, source count and 21cm angular power spectra, matter power spectra, transfer functions
and background evolution. The code is in Python and modern Fortran.

To install the CAMB python package on your computer run::

    pip install camb

or from the source

   python setup.py install

You will need gfortran 6 or higher installed to compile. Binary library builds for python on
Windows are also provided, so these are used instead if no gfortran installation
is found on Windows machines.

See the `CAMB python example notebook <https://camb.readthedocs.org/en/latest/CAMBdemo.html>`_ for a
quick introduction to how to use the CAMB Python package.

The python wrapper provides a module called "camb" documented in the `Python CAMB documentation <https://camb.readthedocs.io/en/latest/>`_.

.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.org/en/latest

To compile the Fortran command-line code run "Make" in the fortran directory. For full details
see the  `ReadMe <https://camb.info/readme.html>`_.

Branches
=============================

The master branch contains latest changes to the main release version.

.. image:: https://secure.travis-ci.org/cmbant/CAMB.png?branch=master
  :target: https://secure.travis-ci.org/cmbant/CAMB/builds
.. image:: https://mybinder.org/badge.svg
  :target: https://mybinder.org/v2/gh/cmbant/camb/master?filepath=docs%2FCAMBdemo.ipynb

The devel branch contains latest less-stable things in development.
The master and devel branches have an integrated test suite, which runs automatically on `Travis <http://travis-ci.org>`_  for new commits and pull requests.
Reference results and test outputs are stored in the `test outputs repository <https://github.com/cmbant/CAMB_test_outputs/>`_. Tests can also be run locally.

To reproduce legacy results, see these branches:

CAMB_sources is the old public `CAMB Sources <http://camb.info/sources/>`_ code.
CAMB_v1 is the old Fortran-oriented (gfortran 4.8-compatible) version as used by the Planck 2018 analysis.

=============

.. raw:: html

    <a href="http://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex.png" height="170px"></a>
    <a href="http://erc.europa.eu/"><img src="https://erc.europa.eu/sites/default/files/content/erc_banner-vertical.jpg" height="200px"></a>
