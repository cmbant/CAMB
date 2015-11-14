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

def find_version():
    version_file = io.open(os.path.join(file_dir, '%s/__init__.py'%package_name)).read()
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
        CAMBDIR = os.path.join(file_dir,'..')
        os.chdir(CAMBDIR)
        if platform.system() == "Windows":
            COMPILER = "gfortran"
            # note that TDM-GCC MingW 5.1 does not work due go general fortran bug.
            FFLAGS = "-shared -static -cpp -fopenmp -O3 -ffast-math -fmax-errors=4"
            SOURCES = "constants.f90 utils.f90 subroutines.f90 inifile.f90 power_tilt.f90 recfast.f90 reionization.f90"\
             " modules.f90 bessels.f90 equations.f90 halofit_ppf.f90 lensing.f90 SeparableBispectrum.f90 cmbmain.f90"\
             " camb.f90 camb_python.f90"
            OUTPUT = r"-o pycamb\camb\%s"%DLLNAME
            scrs = os.listdir(os.getcwd())
            if not has_win_gfortran():
                print('WARNING: gfortran not in path. If you just installed may need to log off and on again.')
            print COMPILER+' '+FFLAGS+' '+SOURCES+' '+OUTPUT
            subprocess.call(COMPILER+' '+FFLAGS+' '+SOURCES+' '+OUTPUT,shell=True)
            COPY = r"copy HighLExtrapTemplate_lenspotentialCls.dat cambpy\camb"
            subprocess.call(COPY,shell=True)
            scrs.append(DLLNAME)
            if not osp.isfile('pycamb/camb/'+DLLNAME) : sys.exit('Compilation failed')
            print "Removing temp files"
            nscrs = os.listdir(os.getcwd())
            for file in nscrs:
                if not file in scrs:
                    os.remove(file)

        else:
            print "Compiling source..."
            MAKE = "make camblib.so"
            subprocess.call(MAKE,shell=True)
            if not osp.isfile('pycamb/camb/camblib.so') : sys.exit('Compilation failed')
            CHOWN = "chmod 755 pycamb/camb/camblib.so"
            subprocess.call(CHOWN,shell=True)
            COPY = "cp HighLExtrapTemplate_lenspotentialCls.dat pycamb/camb"
            subprocess.call(COPY,shell=True)

        os.chdir(file_dir)
        install.run(self)
        print "Cleaning intermediate files..."
        if platform.system() == "Windows":              
            DELETE = 'rmdir /s /q build'
        else:
            DELETE = 'rm -rf build'
        subprocess.call(DELETE,shell=True)
        
        
setup(name= package_name,
      version=find_version(),
      description='CAMB library for python',
      author='Antony Lewis',
      url="http://camb.info/",
      cmdclass={'install': SharedLibrary},
      packages=['camb'],
      package_data = {'camb':[DLLNAME,'HighLExtrapTemplate_lenspotentialCls.dat']},
      classifiers=[
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
      ],
      keywords=['cosmology', 'CAMB']
     )
