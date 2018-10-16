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
    if ok and is_windows:
        version_str = str(subprocess.check_output("gfortran -dumpmachine", shell=True))
        ok = is32Bit and 'i686' in version_str or not is32Bit and 'x86_64' in version_str
    if not ok and msg:
        try:
            with open(os.devnull, 'w') as devnull:
                ifort = subprocess.check_output("ifort -v", shell=True, stderr=devnull)
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
    version_file = io.open(os.path.join(file_dir, 'camb', '__init__.py')).read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        version = version_match.group(1)
        commit = os.getenv('TRAVIS_BUILD_NUMBER')
        if commit and not os.getenv('TRAVIS_TAG'):
            version += '.' + commit
        return version
    raise RuntimeError("Unable to find version string.")


def get_forutils():
    fpath = os.getenv('FORUTILSPATH')
    if not fpath:
        dirs = ['.', '..', '../..']
        for dir in dirs:
            path = os.path.join(dir, 'forutils')
            if os.path.isdir(path):
                fpath = path
                break
        if not fpath:
            try:
                print('forutils not found, attempting to install using git')
                if subprocess.call("git clone --depth=1 https://github.com/cmbant/forutils", shell=True) == 0:
                    fpath = os.path.join('.', 'forutils')
            except Exception:
                print('Failed to install using git')
    if not path:
        raise Exception('Install forutils from https://github.com/cmbant/forutils; or set FORUTILSPATH variable')
    return fpath


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
        ifort = None
        if not is_windows:
            try:
                ifort = str(subprocess.check_output("ifort -v", shell=True))
            except Exception:
                pass
        if not ifort:
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
                SOURCES = " constants.f90 classes.f90 subroutines.f90 power_tilt.f90 recfast.f90 reionization.f90 DarkEnergyInterface.f90  modules.f90" \
                          " bessels.f90 equations.f90 DarkEnergyFluid.f90 DarkEnergyPPF.f90 halofit_ppf.f90 lensing.f90 SeparableBispectrum.f90" \
                          " cmbmain.f90 camb.f90 camb_python.f90"
                OUTPUT = r"-o %s\camb\%s" % (pycamb_path, DLLNAME)
                fpath = get_forutils()
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
            get_forutils()
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
            fpath = get_forutils()
            import shutil, glob
            try:
                for dir, root in zip(['fortran', 'fortran/forutils'], ['..', fpath]):
                    os.mkdir(dir)
                    for pat in ['*.*90', 'Makefile*']:
                        for file in glob.glob(os.path.join(root, pat)):
                            shutil.copy(file, dir)
                shutil.copy('..' + os.sep + 'HighLExtrapTemplate_lenspotentialCls.dat', 'fortran')
                sdist.run(self)
            finally:
                shutil.rmtree('fortran')
        else:
            sdist.run(self)


if __name__ == "__main__":
    setup(name=os.getenv('CAMB_PACKAGE_NAME', 'camb'),
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
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7'
          ],
          keywords=['cosmology', 'CAMB']
          )
