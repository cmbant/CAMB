===================
CAMB
===================
:CAMB:  Code for Anisotropies in the Microwave Background, Fortran 95 code
:Homepage: http://camb.info/

  
Description and installation
=============================

For full details see the `ReadMe <http://camb.info/readme.html>`_.

The master and devel branches have an integrated test suite, which runs automatically on `Travis <http://travis-ci.org>`_  for new commits and pull requests.
Reference results and test outputs are stored in the `test outputs repository <https://github.com/cmbant/CAMB_test_outputs/>`_. Tests can also be run locally.

Branches
=============================

The master branch contains latest changes to the main release version.

.. image:: https://secure.travis-ci.org/cmbant/CAMB.png?branch=master
  :target: https://secure.travis-ci.org/cmbant/CAMB

The devel branch is a developement version, which integrates CAMB and CAMB sources, and uses Fortran 2008 (and hence requires ifort 14+ or gfortran 6+).

.. image:: https://secure.travis-ci.org/cmbant/CAMB.png?branch=devel
  :target: https://secure.travis-ci.org/cmbant/CAMB


CAMB_sources is the updated public `CAMB Sources <http://camb.info/sources/>`_ code.

