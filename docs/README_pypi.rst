===================
CAMB
===================
:CAMB: Code for Anisotropies in the Microwave Background
:Author: Antony Lewis and Anthony Challinor
:Homepage: https://camb.info/

.. image:: https://img.shields.io/pypi/v/camb.svg?style=flat
        :target: https://pypi.python.org/pypi/camb/
.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.io/en/latest

Description
============


CAMB is a cosmology code for calculating cosmological observables, including
CMB, lensing, source count and 21cm angular power spectra, matter power spectra, transfer functions
and background evolution. The code is in Python, with numerical code implemented in fast modern Fortran.

See the `CAMB python example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_ for a
quick introduction to how to use the CAMB Python package.

To install the CAMB on your computer run::

    pip install camb [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
Binary wheels are provided for most systems; to compile from source you will need gfortran 6 or higher installed
(usually installed as part of gcc). gfortran is also required for just-in-time compilation of some functions like custom sources.

The python module is  "camb". see the `CAMB documentation <https://camb.readthedocs.io/en/latest/>`_.
