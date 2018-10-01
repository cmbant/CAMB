#!/usr/bin/env python

import sys
import platform
import subprocess
import io
import re
import os
from distutils.command.build import build
from distutils.command.install import install
from distutils.command.sdist import sdist
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

gfortran_min = '4.9'


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
        try:
            ifort = subprocess.check_output("ifort -v", shell=True)
        except:
            ifort = False
        if not ifort:
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


class SharedLibrary(build, object):
    cluster = False

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
            SOURCES = "constants.f90 utils.f90 subroutines.f90 inifile.f90 power_tilt.f90 recfast.f90 reionization.f90" \
                      " modules.f90 bessels.f90 equations.f90 halofit_ppf.f90 lensing.f90 SeparableBispectrum.f90 cmbmain.f90" \
                      " camb.f90 camb_python.f90"
            OUTPUT = r"-o %s\camb\%s" % (pycamb_path, DLLNAME)
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
                print(COMPILER + ' ' + FFLAGS + ' ' + SOURCES + ' ' + OUTPUT)
                subprocess.call(COMPILER + ' ' + FFLAGS + ' ' + SOURCES + ' ' + OUTPUT, shell=True)
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
            subprocess.call("make camblib.so PYCAMB_OUTPUT_DIR=%s/camb/ CLUSTER_SAFE=%d" %
                            (pycamb_path, int(self.cluster)), shell=True)
            so_file = os.path.join(pycamb_path, 'camb', 'camblib.so')
            if not os.path.isfile(so_file): sys.exit('Compilation failed')
            subprocess.call("chmod 755 %s" % so_file, shell=True)
            subprocess.call("cp HighLExtrapTemplate_lenspotentialCls.dat %s/camb" % pycamb_path, shell=True)

        os.chdir(file_dir)
        build.run(self)


class SharedLibraryCluster(SharedLibrary):
    cluster = True

    def run(self):
        super(SharedLibraryCluster, self).run()


class CustomInstall(install):
    def run(self):
        self.run_command('build')
        install.run(self)


class CustomSdist(sdist):
    def read_template(self):
        sdist.read_template(self)
        self.filelist.process_template_line('recursive-include fortran Makefile* *.dat *.?90')

    def run(self):
        if not os.path.exists('fortran'):
            import shutil, glob
            os.mkdir('fortran')
            try:
                for file in glob.glob('..' + os.sep + '*.*90'):
                    shutil.copy(file, 'fortran')
                shutil.copy('..' + os.sep + 'Makefile', 'fortran')
                shutil.copy('..' + os.sep + 'Makefile_main', 'fortran')
                shutil.copy('..' + os.sep + 'HighLExtrapTemplate_lenspotentialCls.dat', 'fortran')
                sdist.run(self)
            finally:
                shutil.rmtree('fortran')
        else:
            sdist.run(self)


if __name__ == "__main__":
    setup(name=package_name,
          version=find_version(),
          description='Code for Anisotropies in the Microwave Background',
          long_description=get_long_description(),
          author='Antony Lewis',
          url="https://camb.info/",
          cmdclass={'build': SharedLibrary, 'build_cluster': SharedLibraryCluster,
                    'install': CustomInstall, 'sdist': CustomSdist},
          packages=['camb', 'camb_tests'],
          package_data={'camb': [DLLNAME, 'HighLExtrapTemplate_lenspotentialCls.dat',
                                 'PArthENoPE_880.2_marcucci.dat', 'PArthENoPE_880.2_standard.dat',
                                 'PRIMAT_Yp_DH_Error.dat']},
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
