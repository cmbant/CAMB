===================
CAMB
===================
:CAMB: Code for Anisotropies in the Microwave Background
:Version: 0.1.1
:Author: Antony Lewis and Anthony Challinor
:Homepage: https://camb.info/

.. image:: http://img.shields.io/pypi/v/camb.svg?style=flat
        :target: https://pypi.python.org/pypi/camb/
.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.org/en/latest

Description
============

CAMB is a cosmology code for calculating CMB and matter power spectra,
as well as general utility function for cosmological calculations. The main code is written in Fortran, which this
Python package encapulates for convenient python programming and plotting.

For an introduction to CAMB see the main `home page <http://camb.info/>`_ , and the
`CAMB python example notebook <http://camb.readthedocs.org/en/latest/CAMBdemo.html>`_ for a quick
introduction to how to use the CAMB Python package.

To install the CAMB python package on your computer run::

    pip install --egg camb
  
You will need gfortran 4.9 or higher installed to compile. Binary files for Windows are also provided, so these are used instead if no
gfortran installation is found on Windows machines.
