#!/usr/bin/env python

import os
import sys
import os.path as osp
import platform
import subprocess
import io
import re
import os
from distutils.command.install import install
import struct

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

package_name = 'camb'
if platform.system() == "Windows":
    DLLNAME = 'cambdll.dll'
else:
    DLLNAME = 'camblib.so'
file_dir = os.path.abspath(os.path.dirname(__file__))

is32Bit = struct.calcsize("P") == 4


def find_version():
    version_file = io.open(os.path.join(file_dir, '%s/__init__.py' % package_name)).read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def has_win_gfortran():
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        if os.path.isfile(os.path.join(path, 'gfortran.exe')):
            return True
    return False


class SharedLibrary(install):
    def run(self):
        CAMBDIR = os.path.join(file_dir, '..')
        os.chdir(CAMBDIR)
        if platform.system() == "Windows":
            COMPILER = "gfortran"
            # note that TDM-GCC MingW 5.1 does not work due go general fortran bug.
            # This works: http://sourceforge.net/projects/mingw-w64/?source=typ_redirect
            # but need to use 32bit compiler to build 32 bit dll (contrary to what is implied)
            FFLAGS = "-shared -static -cpp -fopenmp -O3 -ffast-math -fmax-errors=4"
            if is32Bit: FFLAGS = "-m32 " + FFLAGS
            SOURCES = "constants.f90 utils.f90 subroutines.f90 inifile.f90 power_tilt.f90 recfast.f90 reionization.f90" \
                      " modules.f90 bessels.f90 equations.f90 halofit_ppf.f90 lensing.f90 SeparableBispectrum.f90 cmbmain.f90" \
                      " camb.f90 camb_python.f90"
            OUTPUT = r"-o pycamb\camb\%s" % DLLNAME
            scrs = os.listdir(os.getcwd())
            if not has_win_gfortran():
                print('WARNING: gfortran not in path (if you just installed you may need to log off and on again).')
                print('         You can get a Windows gfortan build from http://sourceforge.net/projects/mingw-w64/')
                print('         (get the %s version to match this python installation)'%(('x86_64','i686')[is32Bit]))
                print('Using pre-compiled binaries instead - any local changes will be ignored...')
                COPY = r'copy /Y pycamb\dlls\%s pycamb\camb\%s' % (('cambdll_x64.dll', DLLNAME)[is32Bit], DLLNAME)
                subprocess.call(COPY, shell=True)
            else:
                print(COMPILER + ' ' + FFLAGS + ' ' + SOURCES + ' ' + OUTPUT)
                subprocess.call(COMPILER + ' ' + FFLAGS + ' ' + SOURCES + ' ' + OUTPUT, shell=True)
            COPY = r"copy /Y HighLExtrapTemplate_lenspotentialCls.dat pycamb\camb"
            subprocess.call(COPY, shell=True)
            scrs.append(DLLNAME)
            if not osp.isfile('pycamb/camb/' + DLLNAME): sys.exit('Compilation failed')
            print("Removing temp files")
            nscrs = os.listdir(os.getcwd())
            for file in nscrs:
                if not file in scrs:
                    os.remove(file)

        else:
            print("Compiling source...")
            MAKE = "make camblib.so"
            subprocess.call(MAKE, shell=True)
            if not osp.isfile('pycamb/camb/camblib.so'): sys.exit('Compilation failed')
            CHOWN = "chmod 755 pycamb/camb/camblib.so"
            subprocess.call(CHOWN, shell=True)
            COPY = "cp HighLExtrapTemplate_lenspotentialCls.dat pycamb/camb"
            subprocess.call(COPY, shell=True)

        os.chdir(file_dir)
        install.run(self)
        print("Cleaning intermediate files...")
        if platform.system() == "Windows":
            DELETE = 'rmdir /s /q build'
        else:
            DELETE = 'rm -rf build'
        subprocess.call(DELETE, shell=True)


setup(name=package_name,
      version=find_version(),
      description='CAMB library for python',
      author='Antony Lewis',
      url="http://camb.info/",
      cmdclass={'install': SharedLibrary},
      packages=['camb', 'camb_tests'],
      package_data={'camb': [DLLNAME, 'HighLExtrapTemplate_lenspotentialCls.dat']},
      test_suite='camb_tests',
      classifiers=[
          "Programming Language :: Python :: 2",
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
      ],
      keywords=['cosmology', 'CAMB']
      )
