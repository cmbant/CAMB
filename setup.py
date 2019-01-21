#!/usr/bin/env python

import sys
import platform
import subprocess
import io
import re
import os
import shutil
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from distutils.core import Command
import struct

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
    with open(os.path.join('docs', 'README_pypi.rst')) as f:
        return f.read()


def get_gfortran_version():
    try:
        return subprocess.check_output("gfortran -dumpversion", shell=True).decode().strip()
    except subprocess.CalledProcessError:
        return None


def check_ifort():
    try:
        return subprocess.check_output("ifort -v", shell=True, stderr=subprocess.STDOUT)
    except:
        return False


def check_gfortran(version=gfortran_min, msg=False, import_fail_ok=True):
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
        raise Exception(
            'You need ifort or gfortran %s or higher to compile (found: %s).' % (
                version, gfortran_version))

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

    def git_install_forutils():
        try:
            print('forutils not found, attempting to install using git...')
            os.chdir('..')
            try:
                if subprocess.call("git clone --depth=1 https://github.com/cmbant/forutils", shell=True) == 0:
                    return os.path.join('..', 'forutils')
            finally:
                os.chdir('fortran')
        except Exception:
            print('Failed to install using git')

    if not fpath:
        dirs = ['..', '../..']
        for dir in dirs:
            path = os.path.join(dir, 'forutils')
            if os.path.isdir(path):
                fpath = path
                main_dir = dir
                break
        if not fpath:
            fpath = git_install_forutils()
        elif not os.path.exists(os.path.join(fpath, 'Makefile')):
            if os.path.isdir(os.path.join('..', '.git')):
                # submodule may not be pulled
                try:
                    try:
                        os.chdir('..')
                        print('forutils directory found but no Makefile. Attempting to clone submodule...')
                        if subprocess.call("git submodule update --init --recursive", shell=True, cwd=main_dir) != 0:
                            raise Exception()
                    finally:
                        os.chdir('fortran')
                except Exception:
                    fpath = None
                    print('Failed to install forutils using git')
            else:
                os.rmdir(fpath)
                if not os.path.exists(fpath):
                    fpath = git_install_forutils()
                else:
                    fpath = None

    if not fpath:
        raise Exception(
            'Install forutils from https://github.com/cmbant/forutils, pull the submodule, or set FORUTILSPATH variable')
    return fpath


def makefile_dict(filename):
    # this is very non-general, just enough for pulling source file names from Makefile
    with io.open(filename, 'r') as f:
        lines = f.readlines()
    vals = {}
    lastval = None
    append = False
    for line in lines:
        parts = line.split('\\')
        line = parts[0].strip()
        if '?=' in line:
            key, val = line.split('?=')
            env = os.environ.get(key.strip(), None)
            if env:
                vals[key] = env
                lastval = None
                append = False
                continue
            else:
                line = line.replace('?=', '=')
        if append and lastval:
            vals[lastval] += ' ' + line
        elif '=' in line:
            if len(line.split('=')) == 2 and ':' not in line:
                lastval, value = line.split('=')
                lastval = lastval.strip()
                vals[lastval] = value.strip()
        else:
            lastval = None
        append = len(parts) > 1

    def repl(groups):
        if groups.group(1) in vals:
            return vals[groups.group(1)]
        else:
            return groups.group(0)

    for key, value in vals.items():
        if '$' in value:
            vals[key] = re.sub(r'\$\((\w+)\)', repl, value)
    return vals


def make_library(cluster=False):
    CAMBDIR = os.path.join(file_dir, 'fortran')
    pycamb_path = '..'
    os.chdir(CAMBDIR)
    lib_file = os.path.join(pycamb_path, 'camb', DLLNAME)
    if is_windows or not check_ifort():
        ok, gfortran_version = check_gfortran(msg=not is_windows)
    if is_windows:
        COMPILER = "gfortran"
        # note that TDM-GCC MingW 5.1 does not work due go general fortran bug.
        # This works: http://sourceforge.net/projects/mingw-w64/?source=typ_redirect
        # but need to use 32bit compiler to build 32 bit dll (contrary to what is implied)
        FFLAGS = "-shared -static -cpp -fopenmp -O3 -fmax-errors=4"
        if is32Bit: FFLAGS = "-m32 " + FFLAGS
        if not ok:
            print(
                'WARNING: gfortran %s or higher not in path (if you just installed you may need to log off and on again).' %
                gfortran_min)
            print('         You can get a Windows gfortran build from http://sourceforge.net/projects/mingw-w64/')
            print('         (get the %s version to match this python installation)' % (('x86_64', 'i686')[is32Bit]))
            if is32Bit:
                raise IOError('No 32bit Windows DLL provided, you need to build or use 64 bit python')
            else:
                print('Using pre-compiled binary instead - any local changes will be ignored...')
        else:
            fpath = get_forutils()
            makefile = makefile_dict('Makefile_main')
            SOURCES = makefile['SOURCEFILES'].split() + [makefile['CAMBSO'].replace('.f90', '')]
            FORUTILS = [os.path.join(fpath, f.replace('.f90', '')) for f in
                        makefile_dict(os.path.join(fpath, 'Makefile'))['SRCS'].replace('MatrixUtils.f90',
                                                                                       '').split()]
            tmpdir = 'WinDLL' + ('', '32')[is32Bit]
            if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
            ofiles = []
            for source in FORUTILS + SOURCES:
                #manual Make using dependency files if available
                outroot = os.path.join(tmpdir, os.path.split(source)[1])
                fout = outroot + '.o'
                ofiles += [fout]
                modified = not os.path.exists(fout)
                if not modified:
                    o_time = os.path.getmtime(fout)
                    modified = o_time < os.path.getmtime(source + '.f90') or not os.path.exists(outroot + '.d')
                    if not modified:
                        with io.open(outroot + '.d', 'r') as f:
                            for dependence in " ".join(f.readlines()).replace("\\\n", "").split(':')[1].strip().split():
                                if os.path.getmtime(dependence) > o_time:
                                    modified = True
                                    break

                if modified:
                    cmd = COMPILER + ' ' + FFLAGS + ' ' + source + '.f90 -MMD -c -o %s -J%s' % (fout, tmpdir)
                    print(cmd)
                    if subprocess.call(cmd, shell=True) != 0:
                        raise IOError('Compilation failed')
            if os.path.exists(lib_file):
                # raise an exception if the file in use and cannot be deleted
                try:
                    os.remove(lib_file)
                except OSError:
                    raise IOError('dll file in use. Stop python codes and notebook kernels that are using camb.')
            print('Compiling sources...')
            cmd = COMPILER + ' ' + FFLAGS + ' ' + " ".join(ofiles) + ' -o %s -J%s' % (lib_file, tmpdir)
            print(cmd)
            if subprocess.call(cmd, shell=True) != 0:
                raise IOError('Compilation failed')
    else:
        get_forutils()
        print("Compiling source...")
        subprocess.call("make camblib.so PYCAMB_OUTPUT_DIR=%s/camb/ CLUSTER_SAFE=%d" %
                        (pycamb_path, int(cluster)), shell=True)
        subprocess.call("chmod 755 %s" % lib_file, shell=True)

    if not os.path.isfile(os.path.join(pycamb_path, 'camb', DLLNAME)):
        sys.exit('Compilation failed')
    tem_file = 'HighLExtrapTemplate_lenspotentialCls.dat'
    tem = os.path.join(pycamb_path, 'camb', tem_file)
    if not os.path.exists(tem) or os.path.getmtime(tem) < os.path.getmtime(tem_file):
        shutil.copy(tem_file, tem)

    os.chdir(file_dir)


class MakeLibrary(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        make_library(False)


class MakeLibraryCluster(MakeLibrary):

    def run(self):
        make_library(True)


class SharedLibrary(build_py):

    def run(self):
        make_library(False)
        build_py.run(self)


class SharedLibraryCluster(SharedLibrary):

    def run(self):
        make_library(True)
        build_py.run(self)


class DevelopLibrary(develop):

    def run(self):
        make_library(False)
        develop.run(self)


class DevelopLibraryCluster(develop):

    def run(self):
        make_library(True)
        develop.run(self)


if __name__ == "__main__":
    setup(name=os.getenv('CAMB_PACKAGE_NAME', 'camb'),
          version=find_version(),
          description='Code for Anisotropies in the Microwave Background',
          long_description=get_long_description(),
          author='Antony Lewis',
          url="https://camb.info/",
          zip_safe=False,
          cmdclass={'build_py': SharedLibrary, 'build_cluster': SharedLibraryCluster,
                    'make': MakeLibrary, 'make_cluster': MakeLibraryCluster,
                    'develop': DevelopLibrary, 'develop_cluster': DevelopLibraryCluster},
          packages=['camb', 'camb_tests'],
          package_data={'camb': [DLLNAME, 'HighLExtrapTemplate_lenspotentialCls.dat',
                                 'PArthENoPE_880.2_marcucci.dat', 'PArthENoPE_880.2_standard.dat',
                                 'PRIMAT_Yp_DH_Error.dat']},
          test_suite='camb_tests',
          entry_points={
              'console_scripts': [
                  'camb=camb._command_line:run_command_line',
              ]},
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Operating System :: OS Independent',
              'Intended Audience :: Science/Research',
              'Topic :: Scientific/Engineering :: Astronomy',
              "Programming Language :: Python :: 2",
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3',
              'Programming Language :: Python :: 3.5',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7'
          ],
          keywords=['cosmology', 'CAMB', 'CMB']
          )
