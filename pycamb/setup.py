#!/usr/bin/env python

import sys
import platform
import subprocess
import io
import re
import os
from distutils.command.build import build
from distutils.command.install import install
import struct

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

package_name = 'camb'

is_windows = platform.system() == "Windows"
if is_windows:
    DLLNAME = 'cambdll.dll'
else:
    DLLNAME = 'camblib.so'
file_dir = os.path.abspath(os.path.dirname(__file__))

os.chdir(file_dir)

is32Bit = struct.calcsize("P") == 4

gfortran_min = '6'


def get_long_description():
    with open('README.rst') as f:
        return f.read()


def get_gfortran_version():
    try:
        return subprocess.check_output("gfortran -dumpversion", shell=True).decode().strip()
    except subprocess.CalledProcessError:
        return None


def check_gfortran(version=gfortran_min, msg=True, exit=False, import_fail_ok=True):
    gfortran_version = get_gfortran_version()
    version = str(version)
    if gfortran_version:
        try:
            from pkg_resources import parse_version
            ok = parse_version(version) <= parse_version(gfortran_version)
        except ImportError:
            ok = import_fail_ok
            pass
    else:
        ok = False
    if not ok and msg:
        raise Exception(
            'You need gfortran %s or higher to compile (found: %s).' % (
                version, gfortran_version))
    if exit:
        sys.exit(1 if ok else 0)
    return ok, gfortran_version


def find_version():
    version_file = io.open(os.path.join(file_dir, '%s/__init__.py' % package_name)).read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        version = version_match.group(1)
        commit = os.getenv('TRAVIS_BUILD_NUMBER')
        if commit and not os.getenv('TRAVIS_TAG'):
            version += '.' + commit
        return version
    raise RuntimeError("Unable to find version string.")


class SharedLibrary(build):
    def run(self):
        CAMBDIR = os.path.join(file_dir, '..')
        if not os.path.exists(os.path.join(CAMBDIR, 'lensing.f90')):
            CAMBDIR = os.path.join(file_dir, 'fortran')  # pypi install
            pycamb_path = '..'
        else:
            pycamb_path = 'pycamb'
        os.chdir(CAMBDIR)
        ok, gfortran_version = check_gfortran(msg=not is_windows)
        if is_windows:
            COMPILER = "gfortran"
            # note that TDM-GCC MingW 5.1 does not work due go general fortran bug.
            # This works: http://sourceforge.net/projects/mingw-w64/?source=typ_redirect
            # but need to use 32bit compiler to build 32 bit dll (contrary to what is implied)
            FFLAGS = "-shared -static -cpp -fopenmp -O3 -ffast-math -fmax-errors=4"
            if is32Bit: FFLAGS = "-m32 " + FFLAGS
            scrs = os.listdir(os.getcwd())
            if not ok:
                print(
                    'WARNING: gfortran %s or higher not in path (if you just installed you may need to log off and on again).' %
                    gfortran_min)
                print('         You can get a Windows gfortran build from http://sourceforge.net/projects/mingw-w64/')
                print('         (get the %s version to match this python installation)' % (('x86_64', 'i686')[is32Bit]))
                print('Using pre-compiled binaries instead - any local changes will be ignored...')
                subprocess.call(r'copy /Y %s\dlls\%s %s\camb\%s' % (
                    pycamb_path, ('cambdll_x64.dll', DLLNAME)[is32Bit], pycamb_path, DLLNAME), shell=True)
            else:
                FORUTILS = "MiscUtils.f90 StringUtils.f90 ArrayUtils.f90 MpiUtils.f90 FileUtils.f90 " \
                           "IniObjects.f90 RandUtils.f90 ObjectLists.f90 RangeUtils.f90 Interpolation.f90"
                SOURCES = " constants.f90 subroutines.f90 power_tilt.f90 recfast.f90 reionization.f90 DarkEnergyInterface.f90  modules.f90" \
                          " bessels.f90 equations.f90 DarkEnergyFluid.f90 DarkEnergyPPF.f90 halofit_ppf.f90 lensing.f90 SeparableBispectrum.f90" \
                          " cmbmain.f90 camb.f90 camb_python.f90"
                OUTPUT = r"-o %s\camb\%s" % (pycamb_path, DLLNAME)
                fpath = os.getenv('FORUTILSPATH')
                if not fpath:
                    dirs = ['', "..\\", "..\\..\\"]
                    for dir in dirs:
                        if os.path.isdir(dir + 'forutils'):
                            fpath = dir + 'forutils'
                            break
                    if not fpath:
                        raise Exception(
                            'First install forutils from https://github.com/cmbant/forutils; or set FORUTILSPATH variable')
                FORUTILS = " ".join([os.path.join(fpath, p) for p in FORUTILS.split()])
                print('Compiling sources...')
                cmd = COMPILER + ' ' + FFLAGS + ' ' + FORUTILS + ' ' + SOURCES + ' ' + OUTPUT
                print(cmd)
                subprocess.call(cmd, shell=True)
            subprocess.call(r"copy /Y HighLExtrapTemplate_lenspotentialCls.dat %s\camb" % pycamb_path, shell=True)
            scrs.append(DLLNAME)
            if not os.path.isfile(os.path.join(pycamb_path, 'camb', DLLNAME)): sys.exit('Compilation failed')
            print("Removing temp files")
            nscrs = os.listdir(os.getcwd())
            for file in nscrs:
                if not file in scrs:
                    os.remove(file)
        else:
            print("Compiling source...")
            subprocess.call("make camblib.so COMPILER=gfortran PYCAMB_OUTPUT_DIR=%s/camb/" % pycamb_path, shell=True)
            so_file = os.path.join(pycamb_path, 'camb', 'camblib.so')
            if not os.path.isfile(so_file): sys.exit('Compilation failed')
            subprocess.call("chmod 755 %s" % so_file, shell=True)
            subprocess.call("cp HighLExtrapTemplate_lenspotentialCls.dat %s/camb" % pycamb_path, shell=True)

        os.chdir(file_dir)
        build.run(self)


class CustomInstall(install):
    def run(self):
        self.run_command('build')
        install.run(self)


if __name__ == "__main__":
    setup(name=package_name,
          version=find_version(),
          description='Code for Anisotropies in the Microwave Background',
          long_description=get_long_description(),
          author='Antony Lewis',
          url="http://camb.info/",
          cmdclass={'build': SharedLibrary, 'install': CustomInstall},
          packages=['camb', 'camb_tests'],
          package_data={'camb': [DLLNAME, 'HighLExtrapTemplate_lenspotentialCls.dat',
                                 'PArthENoPE_880.2_marcucci.dat', 'PArthENoPE_880.2_standard.dat']},
          test_suite='camb_tests',
          classifiers=[
              "Programming Language :: Python :: 2",
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3',
              'Programming Language :: Python :: 3.4',
              'Programming Language :: Python :: 3.5',
              'Programming Language :: Python :: 3.6'
          ],
          keywords=['cosmology', 'CAMB']
          )
