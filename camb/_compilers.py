import subprocess
import os
import struct
import platform
import re

is_windows = platform.system() == "Windows"

is_32_bit = struct.calcsize("P") == 4

compiler_environ = os.environ.copy()

gfortran_min = '6'
gfortran_bits = ('x86_64', 'i686')[is_32_bit]


def call_command(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, env=compiler_environ).decode().strip()
    except:
        return None


def get_ifort_version():
    return call_command("ifort -v")


def get_gfortran_version(command='gfortran'):
    ver = call_command(command + " -dumpversion")
    if ver and '.' not in ver:
        ver = call_command(command + " -dumpfullversion")
    return ver


def check_ifort():
    return get_ifort_version() or False


def check_gfortran(version=gfortran_min, msg=False, retry=False):
    from packaging.version import parse as parse_version
    global compiler_environ
    gfortran_version = get_gfortran_version()
    version = str(version)
    if gfortran_version:
        ok = parse_version(version) <= parse_version(gfortran_version)
    else:
        ok = False
    if not ok and is_windows and not retry:
        newpath = None
        best_version = gfortran_min
        for root in (os.path.expanduser('~'), os.environ["ProgramFiles"]):
            mingw = os.path.join(root, 'mingw-w64')
            if not os.path.isdir(mingw):
                mingw = os.path.join(root, 'mingw64')

            if os.path.isdir(mingw):
                # look for mingw installation
                dirs = [name for name in os.listdir(mingw) if
                        gfortran_bits in name and os.path.isdir(os.path.join(mingw, name))]
                for i, x in enumerate(dirs):
                    if '.' in x:
                        ver = x.split('-')[1]
                        if parse_version(best_version) <= parse_version(ver):
                            best_version = ver
                            newpath = os.path.join(mingw, x, 'mingw64', 'bin')
                bin = os.path.join(mingw, 'bin')
                if os.path.exists(bin):
                    ver = get_gfortran_version('"' + os.path.join(bin, 'gfortran') + '"')
                    if ver and parse_version(best_version) <= parse_version(ver):
                        best_version = ver
                        newpath = bin
        if newpath:
            if os.path.exists(os.path.join(newpath, 'gfortran.exe')):
                if not compiler_environ["PATH"].startswith(newpath):
                    compiler_environ["PATH"] = newpath + ';' + compiler_environ["PATH"]
                return check_gfortran(version, msg, retry=True)
    if ok and is_windows:
        version_str = str(subprocess.check_output("gfortran -dumpmachine", shell=True, env=compiler_environ))
        ok = gfortran_bits in version_str
    if not ok and msg:
        raise Exception(
            'You need ifort or gfortran %s or higher to compile (found: %s).\nSee %s' % (
                version, gfortran_version, 'https://camb.readthedocs.io/en/latest/fortran_compilers.html'))

    return ok, gfortran_version


def makefile_dict(filename):
    # this is very non-general, just enough for pulling source file names from Makefile
    with open(filename, 'r') as f:
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
