#!/usr/bin/env python

import sys
import subprocess
import re
import os
import shutil
from typing import Any
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from distutils.command.clean import clean
from distutils.core import Command

file_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(file_dir)

sys.path.insert(0, os.path.join(file_dir, 'camb'))
_compile: Any = __import__('_compilers')

if _compile.is_windows:
    DLLNAME = 'cambdll.dll'
else:
    DLLNAME = 'camblib.so'


def get_long_description():
    with open(os.path.join('docs', 'README_pypi.rst')) as f:
        return f.read()


def find_version():
    version_file = open(os.path.join(file_dir, 'camb', '__init__.py')).read()
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
            fbranch = os.getenv('FORUTILSBRANCH', '1.0.3' if os.environ.get('CONDA_BUILD') else 'master')
            try:
                if subprocess.call("git clone --branch %s --depth=1 https://github.com/cmbant/forutils" % fbranch,
                                   shell=True) == 0:
                    return os.path.join('..', 'forutils')
            finally:
                os.chdir('fortran')
        except Exception:
            print('Failed to install using git')

    if not fpath:
        dirs = ['..', '..' + os.sep + '..']
        for _dir in dirs:
            path = os.path.join(_dir, 'forutils')
            if os.path.isdir(path):
                fpath = path
                main_dir = _dir
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
        raise Exception('Install forutils from https://github.com/cmbant/forutils, '
                        'pull the forutils submodule, or set FORUTILSPATH variable.\n'
                        'If you are cloning with git, use "git clone --recursive"')
    return fpath


def clean_dir(path, rmdir=False):
    if os.path.isdir(path):
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
        if rmdir:
            os.rmdir(path)


def make_library(cluster=False):
    os.chdir(os.path.join(file_dir, 'fortran'))
    pycamb_path = '..'
    lib_file = os.path.join(pycamb_path, 'camb', DLLNAME)
    if _compile.is_windows or not _compile.check_ifort():
        ok, gfortran_version = _compile.check_gfortran(msg=not _compile.is_windows)
        if ok and '8.2.0' in gfortran_version:
            print('WARNING: gfortran 8.2.0 may be buggy and give unreliable results or crashes, upgrade gfortran.')
    if _compile.is_windows:
        COMPILER = "gfortran"
        FFLAGS = "-shared -static -cpp -fopenmp -O3 -fmax-errors=4"
        # FFLAGS = "-shared -static -cpp -fopenmp -g -fbacktrace -ffpe-trap=invalid,overflow,zero " \
        #         "-fbounds-check -fmax-errors=4"
        if _compile.is_32_bit:
            FFLAGS = "-m32 " + FFLAGS
        if not ok:
            print('WARNING: gfortran %s or higher not in path (if you just installed '
                  'you may need to log off and on again).' % _compile.gfortran_min)
            print('         You can get a Windows gfortran build from https://sourceforge.net/projects/mingw-w64/')
            print('         - go to Files, and download MinGW-W64 Online Installer.')
            if _compile.is_32_bit:
                raise IOError('No 32bit Windows DLL provided, you need to build or use 64 bit python')
            else:
                print('Using pre-compiled binary instead - any local changes will be ignored...')
        else:
            fpath = get_forutils()
            makefile = _compile.makefile_dict('Makefile_main')
            SOURCES = makefile['SOURCEFILES'].split()
            FORUTILS = [os.path.join(fpath, f.replace('.f90', '')) for f in
                        _compile.makefile_dict(os.path.join(fpath, 'Makefile'))['SRCS'].replace('MatrixUtils.f90',
                                                                                                '').split()]
            tmpdir = 'WinDLL'
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir)
            ofiles = []
            new_compiler = True
            ver_file = os.path.join(tmpdir, 'compiler.ver')
            if os.path.exists(ver_file):
                with open(ver_file, 'r') as f:
                    new_compiler = gfortran_version != f.readline().strip()
            if new_compiler:
                clean_dir(tmpdir)
                with open(ver_file, 'w') as f:
                    f.write(gfortran_version)

            need_compile = not os.path.exists(lib_file)
            if not need_compile:
                dll_time = os.path.getmtime(lib_file)
            for source in FORUTILS + SOURCES:
                # manual Make using dependency files if available
                outroot = os.path.join(tmpdir, os.path.split(source)[1])
                fout = outroot + '.o'
                ofiles += [fout]
                modified = new_compiler or not os.path.exists(fout)
                if not modified:
                    o_time = os.path.getmtime(fout)
                    modified = o_time < os.path.getmtime(source + '.f90') or not os.path.exists(outroot + '.d')
                    if not modified:
                        with open(outroot + '.d', 'r') as f:
                            for dependence in " ".join(f.readlines()).replace("\\\n", "").split(':')[1].strip().split():
                                if os.path.getmtime(dependence) > o_time:
                                    modified = True
                                    break

                if modified:
                    need_compile = True
                    cmd = COMPILER + ' ' + FFLAGS + ' ' + source + '.f90 -MMD -c -o %s -J%s' % (fout, tmpdir)
                    print(cmd)
                    if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
                        raise IOError('Compilation failed')
                elif not need_compile and dll_time < o_time:
                    need_compile = True

            if need_compile or not os.path.exists(lib_file):
                if os.path.exists(lib_file):
                    # raise an exception if the file in use and cannot be deleted
                    try:
                        os.remove(lib_file)
                    except OSError:
                        raise IOError('dll file in use. Stop python codes and notebook kernels that are using camb.')
                print('Compiling sources...')
                cmd = COMPILER + ' ' + FFLAGS + ' ' + " ".join(ofiles) + ' -o %s -J%s' % (lib_file, tmpdir)
                print(cmd)
                if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
                    raise IOError('Compilation failed')
            else:
                print('DLL up to date.')
    else:
        get_forutils()
        print("Compiling source...")
        subprocess.call("make python PYCAMB_OUTPUT_DIR=%s/camb/ CLUSTER_SAFE=%d" %
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


class CleanLibrary(clean):

    def run(self):
        if _compile.is_windows:
            clean_dir(os.path.join(file_dir, 'fortran', 'WinDLL'), rmdir=True)
        else:
            subprocess.call("make clean", shell=True, cwd=os.path.join(file_dir, 'fortran'))
        clean.run(self)


if __name__ == "__main__":
    setup(name=os.getenv('CAMB_PACKAGE_NAME', 'camb'),
          version=find_version(),
          description='Code for Anisotropies in the Microwave Background',
          long_description=get_long_description(),
          author='Antony Lewis',
          url="https://camb.info/",
          project_urls={
              'Documentation': 'https://camb.readthedocs.io',
              'Source': 'https://github.com/cmbant/camb',
              'Tracker': 'https://github.com/cmbant/camb/issues',
              'Reference': 'https://arxiv.org/abs/astro-ph/9911177',
              'Licensing': 'https://github.com/cmbant/CAMB/blob/master/LICENCE.txt'
          },
          zip_safe=False,
          cmdclass={'build_py': SharedLibrary, 'build_cluster': SharedLibraryCluster,
                    'make': MakeLibrary, 'make_cluster': MakeLibraryCluster, 'clean': CleanLibrary,
                    'develop': DevelopLibrary, 'develop_cluster': DevelopLibraryCluster},
          packages=['camb', 'camb.tests'],
          platforms="any",
          package_data={'camb': [DLLNAME, 'HighLExtrapTemplate_lenspotentialCls.dat',
                                 'PArthENoPE_880.2_marcucci.dat', 'PArthENoPE_880.2_standard.dat',
                                 'PRIMAT_Yp_DH_Error.dat', 'PRIMAT_Yp_DH_ErrorMC_2021.dat']},
          test_suite='camb.tests',
          entry_points={
              'console_scripts': [
                  'camb=camb._command_line:run_command_line',
              ]},
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Operating System :: OS Independent',
              'Intended Audience :: Science/Research',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Programming Language :: Python :: 3',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',
              'Programming Language :: Python :: 3.9'
          ],
          keywords=['cosmology', 'CAMB', 'CMB'],
          install_requires=['scipy>=1.0', 'sympy>=1.0'],
          python_requires='>=3.6'
          )
