.. _fortran-compilers:

Fortran compilers
=========================

CAMB internally uses modern (object-oriented) Fortran 2008 for most numerical calculations, and needs a recent
fortran compiler to build the numerical library. The recommended compilers are

- gfortran version 6.3 or higher
- Intel Fortran (ifort), version 18.0.1 or higher (some things may work with version 14+)

The gfortran compiler is part of the standard "gcc" compiler package, and may be pre-installed on recent unix systems.
Check the version using "gfortran --version".

If you do not have a suitable Fortran compiler, you can get one as follows:

:Mac:
    Download the `binary installation <https://gcc.gnu.org/wiki/GFortranBinaries>`_
:Windows:
    Download gfortran as part of `MinGW-w64 <https://sourceforge.net/projects/mingw-w64/files>`_ (select x86_64 option in the installation program)
:Linux:
    To install from the standard repository use:

     - "sudo apt-get update; sudo apt-get install gfortran"

    On Ubuntu systems where the default gfortran is too old, you can use this to install a later version

     - sudo add-apt-repository ppa:ubuntu-toolchain-r/test
     - sudo apt-get update
     - sudo apt install gfortran-8

    To make this the default gfortran then use

     - mkdir -p gfortran-symlinks
     - ln -s /usr/bin/gfortran-8 gfortran-symlinks/gfortran
     - export PATH=$PWD/gfortran-symlinks:$PATH

    To re-use next time, add gfortran-symlinks directory to your startup settings (.bashrc).

Alternatively you can compile and run in a container or virtual machine: e.g., see `CosmoBox <https://cosmologist.info/CosmoBox>`_.
For example, to run a configured shell in docker where you can install and run camb from the command line (after changing to the camb directory)::

    docker run -v /local/git/path/CAMB:/camb -i -t cmbant/cosmobox

Updating and modified Fortran code
===================================

In the main CAMB source root directory, to re-build the Fortran binary including any
pulled or local changes use::

    python setup.py make

This will also work on Windows as long as you have MinGW-w64 installed as described above.

Note that you will need to close all python instances using camb before you can re-load with an updated library.
This includes in Jupyter notebooks; just re-start the kernel or use::

    import IPython
    IPython.Application.instance().kernel.do_shutdown(True)

If you want to automamatically rebuild the library from Jupyter you can do something like this::

    import subprocess
    import sys
    import os
    src_dir = '/path/to/git/CAMB'
    try:
        subprocess.check_output(r'python "%s" make'%os.path.join(src_dir, 'setup.py'),
                                stderr=subprocess.STDOUT)
        sys.path.insert(0,src_dir)
        import camb
        print('Using CAMB %s installed at %s'%(camb.__version__,
                                        os.path.dirname(camb.__file__)))

    except subprocess.CalledProcessError as E:
        print(E.output.decode())
