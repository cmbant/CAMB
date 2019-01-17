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
    if not fpath:
        dirs = ['.', '..', '../..']
        for dir in dirs:
            path = os.path.join(dir, 'forutils')
            if os.path.isdir(path):
                fpath = path
                main_dir = dir
                break
        if not fpath:
            try:
                print('forutils not found, attempting to install using git...')
                if subprocess.call("git clone --depth=1 https://github.com/cmbant/forutils", shell=True) == 0:
                    fpath = os.path.join('.', 'forutils')
            except Exception:
                print('Failed to install using git')
        elif not os.path.exists(os.path.join(fpath, 'Makefilke')) and os.path.isdir(os.path.join(fpath, '.git')):
            # submodule may not be pulled
            try:
                print('forutils directory found but no Makefile. Attempting to clone submodule...')
                if subprocess.call("git submodule update --init --recursive", shell=True, cwd=main_dir) != 0:
                    raise Exception()
            except Exception:
                print('Failed to install using git')

    if not fpath:
        raise Exception('Install forutils from https://github.com/cmbant/forutils; or set FORUTILSPATH variable')
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
        if groups[1] in vals:
            return vals[groups[1]]
        else:
            return groups[0]

    for key, value in vals.items():
        if '$' in value:
            vals[key] = re.sub(r'\$\((\w+)\)', repl, value)
    return vals


class SharedLibrary(build, object):
    cluster = False

    def run(self):
        CAMBDIR = os.path.join(file_dir, 'fortran')
        pycamb_path = '..'
        os.chdir(CAMBDIR)
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
                    print('Using pre-compiled binaries instead - any local changes will be ignored...')
            else:
                fpath = get_forutils()
                makefile = makefile_dict('Makefile_main')
                SOURCES = makefile['SOURCEFILES'].split() + [makefile['CAMBSO'].replace('.f90', '')]
                FORUTILS = [os.path.join(fpath, f.replace('.f90', '')) for f in
                            makefile_dict(os.path.join(fpath, 'Makefile'))['SRCS'].replace('MatrixUtils.f90',
                                                                                           '').split()]
                tmpdir = 'WinDLL' + ('', '32')[is32Bit]
                if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
                modified = False
                ofiles = []
                for source in FORUTILS + SOURCES:
                    # simplest possible Makefile-free make without making full dependencies
                    fout = os.path.join(tmpdir, os.path.split(source)[1] + '.o')
                    ofiles += [fout]
                    if modified or not os.path.exists(fout) or os.path.getmtime(fout) < os.path.getmtime(
                            source + '.f90'):
                        modified = True
                        cmd = COMPILER + ' ' + FFLAGS + ' ' + source + '.f90 -c -o %s -J%s' % (fout, tmpdir)
                        print(cmd)
                        if subprocess.call(cmd, shell=True) != 0:
                            raise IOError('Compilation failed')
                OUTPUT = r"%s\camb\%s" % (pycamb_path, DLLNAME)
                if os.path.exists(OUTPUT):
                    # raise an exception if the file in use and cannot be deleted
                    try:
                        os.remove(OUTPUT)
                    except OSError:
                        raise IOError('dll file in use. Stop python codes and notebook kernels that are using camb.')
                print('Compiling sources...')
                cmd = COMPILER + ' ' + FFLAGS + ' ' + " ".join(ofiles) + ' -o %s -J%s' % (OUTPUT, tmpdir)
                print(cmd)
                if subprocess.call(cmd, shell=True) != 0:
                    raise IOError('Compilation failed')
            subprocess.call(r"copy /Y HighLExtrapTemplate_lenspotentialCls.dat %s\camb" % pycamb_path, shell=True)
            if not os.path.isfile(os.path.join(pycamb_path, 'camb', DLLNAME)):
                sys.exit('Compilation failed')
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
          scripts=['camb.py'],
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
          keywords=['cosmology', 'CAMB', 'CMB']
          )
