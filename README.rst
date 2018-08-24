===================
CAMB
===================
:CAMB:  Code for Anisotropies in the Microwave Background, Fortran 95 code and python module
:Homepage: http://camb.info/

.. image:: https://img.shields.io/pypi/v/camb.svg?style=flat
        :target: https://pypi.python.org/pypi/camb/
  
Description and installation
=============================

For full details of the Fortran code see the `ReadMe <http://camb.info/readme.html>`_.

The python wrapper provides a module called "camb", source in the "pycamb" folder and documented in the `Python CAMB documentation <https://camb.readthedocs.io/en/latest/>`_.

.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.org/en/latest

The master and devel branches have an integrated test suite, which runs automatically on `Travis <http://travis-ci.org>`_  for new commits and pull requests.
Reference results and test outputs are stored in the `test outputs repository <https://github.com/cmbant/CAMB_test_outputs/>`_. Tests can also be run locally.

Branches
=============================

The master branch contains latest changes to the main release version.

.. image:: https://secure.travis-ci.org/cmbant/CAMB.png?branch=master
  :target: https://secure.travis-ci.org/cmbant/CAMB/builds
.. image:: https://mybinder.org/badge.svg 
  :target: https://mybinder.org/v2/gh/cmbant/camb/master?filepath=pycamb%2Fdocs%2FCAMBdemo.ipynb

The devel branch is a development version, which integrates CAMB and CAMB sources, and uses Fortran 2008 (and hence requires ifort 14+ or gfortran 6+). It also allows runtime switching of the dark energy and halofit model, and easy access to lower-level accuracy parameters.

.. image:: https://secure.travis-ci.org/cmbant/CAMB.png?branch=devel
  :target: https://secure.travis-ci.org/cmbant/CAMB/builds
.. image:: https://mybinder.org/badge.svg
  :target: https://mybinder.org/v2/gh/cmbant/camb/devel?filepath=pycamb%2Fdocs%2FCAMBdemo.ipynb


CAMB_sources is the updated public `CAMB Sources <http://camb.info/sources/>`_ code.

=============

.. raw:: html

    <a href="http://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex.png" height="170px"></a>
    <a href="http://erc.europa.eu/"><img src="https://erc.europa.eu/sites/default/files/content/erc_banner-vertical.jpg" height="200px"></a>
