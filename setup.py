import os
import shutil
import subprocess
import sys
from typing import Any, ClassVar

from setuptools import Command, Extension, setup
from setuptools.command.build_ext import build_ext

try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install

file_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(file_dir)

sys.path.insert(0, os.path.join(file_dir, "camb"))
_compile: Any = __import__("_compilers")

if _compile.is_windows:
    DLLNAME = "cambdll.dll"
else:
    DLLNAME = "camblib.so"


def get_forutils():
    fpath = os.getenv("FORUTILSPATH")

    def git_install_forutils():
        try:
            print("forutils not found, attempting to install using git...")
            os.chdir("..")
            fbranch = os.getenv("FORUTILSBRANCH", "1.1.0" if os.environ.get("CONDA_BUILD") else "master")
            try:
                if (
                    subprocess.call(
                        f"git clone --branch {fbranch} --depth=1 https://github.com/cmbant/forutils", shell=True
                    )
                    == 0
                ):
                    return os.path.join("..", "forutils")
            finally:
                os.chdir("fortran")
        except OSError:
            print("Failed to install using git")

    if not fpath:
        dirs = ["..", ".." + os.sep + ".."]
        for _dir in dirs:
            path = os.path.join(_dir, "forutils")
            if os.path.isdir(path):
                fpath = path
                main_dir = os.path.abspath(_dir)
                break
        if not fpath:
            fpath = git_install_forutils()
        elif not os.path.exists(os.path.join(fpath, "Makefile")):
            if os.path.isdir(os.path.join("..", ".git")):
                # submodule may not be pulled
                try:
                    try:
                        os.chdir("..")
                        print("forutils directory found but no Makefile. Attempting to clone submodule...")
                        if subprocess.call("git submodule update --init --recursive", shell=True, cwd=main_dir) != 0:
                            raise RuntimeError("git submodule update failed")
                    finally:
                        os.chdir("fortran")
                except (OSError, RuntimeError):
                    fpath = None
                    print("Failed to install forutils using git")
            else:
                os.rmdir(fpath)
                if not os.path.exists(fpath):
                    fpath = git_install_forutils()
                else:
                    fpath = None

    if not fpath:
        raise RuntimeError(
            "Install forutils from https://github.com/cmbant/forutils, "
            "pull the forutils submodule, or set FORUTILSPATH variable.\n"
            'If you are cloning with git, use "git clone --recursive"'
        )
    return fpath


def clean_dir(path, rmdir=False):
    if os.path.isdir(path):
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
        if rmdir:
            os.rmdir(path)


def remove_stale_dependency_files(build_dir, expected_forutils_path):
    if not os.path.isdir(build_dir):
        return

    normalized_forutils_path = os.path.abspath(expected_forutils_path).replace("\\", "/")

    for name in os.listdir(build_dir):
        if not name.endswith(".d"):
            continue

        dep_path = os.path.join(build_dir, name)
        try:
            with open(dep_path, encoding="utf-8") as dep_file:
                dep_contents = dep_file.read().replace("\\", "/")
        except OSError:
            continue

        if "/forutils/" in dep_contents and normalized_forutils_path not in dep_contents:
            os.remove(dep_path)


def make_library(cluster=False):
    os.chdir(os.path.join(file_dir, "fortran"))
    pycamb_path = ".."
    lib_file = os.path.join(pycamb_path, "camb", DLLNAME)
    flang_command = None
    if _compile.is_windows:
        gfortran_ok, gfortran_version = _compile.check_gfortran()
        flang_command = _compile.find_flang_command()
        requested_compiler = os.environ.get("COMPILER")
        if requested_compiler == "flang" or (not gfortran_ok and flang_command):
            if not flang_command:
                raise RuntimeError("COMPILER=flang was requested, but LLVM flang was not found on PATH.")
            compiler = flang_command
            compiler_version = _compile.call_command(f'"{compiler}" --version')
            is_flang = True
            lld_link = shutil.which("lld-link")
            llvm_nm = shutil.which("llvm-nm")
            if not lld_link or not llvm_nm:
                raise RuntimeError(
                    "Windows builds with LLVM flang require lld-link and llvm-nm on PATH to export CAMB's Fortran symbols."
                )
        elif gfortran_ok:
            compiler = "gfortran"
            compiler_version = gfortran_version
            is_flang = False
        else:
            _compile.check_gfortran(msg=True)
    elif not _compile.check_ifort():
        ok, gfortran_version = _compile.check_gfortran(msg=True)
        if ok and "8.2.0" in gfortran_version:
            print("WARNING: gfortran 8.2.0 may be buggy and give unreliable results or crashes, upgrade gfortran.")
    if _compile.is_windows:
        if is_flang:
            FFLAGS = "-cpp -fopenmp -O3"
            link_flags = "-shared -fuse-ld=lld"
            module_dir_flag = "-module-dir"
            # -fopenmp's own implicit "-defaultlib:libomp.lib" does not reliably resolve every
            # symbol with lld-link (e.g. __kmpc_fork_call_if, used by "!$OMP ... IF(...)" clauses,
            # is left undefined even though it is present in libomp.lib) - passing the import
            # library explicitly on the link line avoids relying on that fallback path.
            flang_exe = shutil.which(compiler, path=_compile.compiler_environ.get("PATH"))
            flang_bin_dir = os.path.dirname(flang_exe) if flang_exe else None
            libomp_dll = os.path.join(flang_bin_dir, "libomp.dll") if flang_bin_dir else None
            libomp_lib = os.path.join(flang_bin_dir, "..", "lib", "libomp.lib") if flang_bin_dir else None
            if libomp_lib and os.path.exists(libomp_lib):
                link_flags += f' "{libomp_lib}"'
            else:
                print(
                    "WARNING: could not locate libomp.lib next to the flang executable; "
                    "linking may fail to resolve some OpenMP runtime symbols."
                )
        else:
            FFLAGS = "-shared -static -cpp -fopenmp -O3 -fmax-errors=4"
            link_flags = ""
            module_dir_flag = "-J"
        # FFLAGS = "-shared -static -cpp -fopenmp -g -fbacktrace -ffpe-trap=invalid,overflow,zero " \
        #         "-fbounds-check -fmax-errors=4"

        fpath = get_forutils()
        makefile = _compile.makefile_dict("Makefile_main")
        SOURCES = makefile["SOURCEFILES"].split()
        FORUTILS = [
            os.path.join(fpath, f.replace(".f90", ""))
            for f in _compile.makefile_dict(os.path.join(fpath, "Makefile"))["SRCS"]
            .replace("MatrixUtils.f90", "")
            .split()
        ]
        tmpdir = "WinDLL"
        if not os.path.isdir(tmpdir):
            os.mkdir(tmpdir)
        ofiles = []
        new_compiler = True
        ver_file = os.path.join(tmpdir, "compiler.ver")
        if os.path.exists(ver_file):
            with open(ver_file) as f:
                new_compiler = compiler_version != f.readline().strip()
        if new_compiler:
            clean_dir(tmpdir)
            with open(ver_file, "w") as f:
                f.write(compiler_version)

        if is_flang:
            # The conda-forge Windows flang package (unlike e.g. apt.llvm.org's Linux build) ships
            # libomp's C header but no Fortran omp_lib.mod, so "use omp_lib" fails to resolve. Only
            # omp_get_thread_num/omp_get_max_threads/omp_set_num_threads/omp_get_wtime are used
            # (via "use omp_lib, only: ...") anywhere in the codebase, so provide a minimal standin
            # module built fresh each time against whichever flang is active, binding directly to
            # libomp's C-interoperable exports (confirmed via llvm-objdump -p libomp.dll: plain,
            # unmangled "omp_get_thread_num" etc.) rather than depending on any particular Fortran
            # module packaging.
            omp_lib_mod = os.path.join(tmpdir, "omp_lib.mod")
            if not os.path.exists(omp_lib_mod):
                shim_source = os.path.join(tmpdir, "win_flang_omp_lib_shim.f90")
                with open(shim_source, "w") as f:
                    f.write(
                        "module omp_lib\n"
                        "  use iso_c_binding, only: c_int, c_double\n"
                        "  implicit none\n"
                        "  private\n"
                        "  public :: omp_get_thread_num, omp_get_max_threads, omp_set_num_threads, omp_get_wtime\n"
                        "  interface\n"
                        '    function omp_get_thread_num() bind(c, name="omp_get_thread_num")\n'
                        "      import :: c_int\n"
                        "      integer(c_int) :: omp_get_thread_num\n"
                        "    end function omp_get_thread_num\n"
                        '    function omp_get_max_threads() bind(c, name="omp_get_max_threads")\n'
                        "      import :: c_int\n"
                        "      integer(c_int) :: omp_get_max_threads\n"
                        "    end function omp_get_max_threads\n"
                        '    subroutine omp_set_num_threads(num_threads) bind(c, name="omp_set_num_threads")\n'
                        "      import :: c_int\n"
                        "      integer(c_int), value :: num_threads\n"
                        "    end subroutine omp_set_num_threads\n"
                        '    function omp_get_wtime() bind(c, name="omp_get_wtime")\n'
                        "      import :: c_double\n"
                        "      real(c_double) :: omp_get_wtime\n"
                        "    end function omp_get_wtime\n"
                        "  end interface\n"
                        "end module omp_lib\n"
                    )
                cmd = (
                    compiler
                    + " "
                    + FFLAGS
                    + f" {shim_source} -c -o {os.path.join(tmpdir, 'win_flang_omp_lib_shim.o')} {module_dir_flag}{tmpdir}"
                )
                print(cmd)
                if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
                    raise OSError("Compilation failed")

        need_compile = not os.path.exists(lib_file)
        if not need_compile:
            dll_time = os.path.getmtime(lib_file)
        for source in FORUTILS + SOURCES:
            # manual Make using dependency files if available
            outroot = os.path.join(tmpdir, os.path.split(source)[1])
            fout = outroot + ".o"
            ofiles += [fout]
            modified = new_compiler or not os.path.exists(fout)
            if not modified:
                o_time = os.path.getmtime(fout)
                modified = o_time < os.path.getmtime(source + ".f90")
                if not modified and not is_flang:
                    modified = not os.path.exists(outroot + ".d")
                if not modified and not is_flang:
                    with open(outroot + ".d") as f:
                        for dependence in " ".join(f.readlines()).replace("\\\n", "").split(":")[1].strip().split():
                            if os.path.getmtime(dependence) > o_time:
                                modified = True
                                break

            if modified:
                need_compile = True
                dependency_flag = "" if is_flang else " -MMD"
                cmd = (
                    compiler
                    + " "
                    + FFLAGS
                    + " "
                    + source
                    + f".f90{dependency_flag} -c -o {fout} {module_dir_flag}{tmpdir}"
                )
                print(cmd)
                if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
                    raise OSError("Compilation failed")
            elif not need_compile and dll_time < o_time:
                need_compile = True

        if need_compile or not os.path.exists(lib_file):
            if os.path.exists(lib_file):
                # raise an exception if the file in use and cannot be deleted
                try:
                    os.remove(lib_file)
                except OSError:
                    raise OSError("dll file in use. Stop python codes and notebook kernels that are using camb.")
            print("Compiling sources...")
            if is_flang:
                def_file = os.path.join(tmpdir, "cambdll.def")
                nm_output = subprocess.check_output(
                    [llvm_nm, "--extern-only", "--defined-only", *ofiles], text=True, env=_compile.compiler_environ
                )
                symbols = sorted(
                    {
                        symbol
                        for line in nm_output.splitlines()
                        if line and (symbol := line.rsplit(maxsplit=1)[-1]).startswith(("_QM", "camb_"))
                    }
                )
                with open(def_file, "w") as f:
                    f.write("EXPORTS\n" + "\n".join(symbols) + "\n")
                link_flags += f" -Wl,-def:{def_file}"
            cmd = (
                compiler
                + " "
                + FFLAGS
                + " "
                + link_flags
                + " "
                + " ".join(ofiles)
                + f" -o {lib_file} {module_dir_flag}{tmpdir}"
            )
            print(cmd)
            if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
                raise OSError("Compilation failed")
            if is_flang:
                # -fopenmp on Windows flang links dynamically against LLVM's libomp.dll (no static
                # option is available, unlike gfortran's -static libgomp); copy it alongside
                # cambdll.dll so the OpenMP runtime is found without needing PATH changes, since
                # Windows searches a loaded DLL's own directory when resolving its dependencies.
                if libomp_dll and os.path.exists(libomp_dll):
                    shutil.copy(libomp_dll, os.path.join(pycamb_path, "camb", "libomp.dll"))
                else:
                    print(
                        "WARNING: could not locate libomp.dll next to the flang executable; "
                        "the built cambdll.dll may fail to load unless libomp.dll is on PATH."
                    )
        else:
            print("DLL up to date.")
    else:
        if not _compile.call_command("make -v"):
            raise OSError(
                'Build failed - you must have "make" installed. '
                'E.g. on ubuntu install with "sudo apt install make" (or use build-essential package).'
            )
        fpath = get_forutils()
        remove_stale_dependency_files(os.path.join(file_dir, "fortran", "Releaselib"), fpath)
        remove_stale_dependency_files(os.path.join(file_dir, "fortran", "Debuglib"), fpath)
        if os.path.exists(lib_file) and not os.access(lib_file, os.W_OK):
            os.remove(lib_file)
        print("Compiling source...")
        cluster_safe = int(cluster if not os.getenv("GITHUB_ACTIONS") else 1)
        subprocess.call(
            f"make python PYCAMB_OUTPUT_DIR={pycamb_path}/camb/ CLUSTER_SAFE={cluster_safe}",
            shell=True,
        )

        if os.path.isfile(lib_file):
            subprocess.call(f"chmod 755 {lib_file}", shell=True)

    if not os.path.isfile(os.path.join(pycamb_path, "camb", DLLNAME)):
        sys.exit("Compilation failed")
    tem_file = "HighLExtrapTemplate_lenspotentialCls.dat"
    tem = os.path.join(pycamb_path, "camb", tem_file)
    if not os.path.exists(tem) or os.path.getmtime(tem) < os.path.getmtime(tem_file):
        shutil.copy(tem_file, tem)

    os.chdir(file_dir)


class MakeLibrary(Command):
    user_options: ClassVar[list] = []

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


class CleanLibrary(MakeLibrary):
    def run(self):
        if _compile.is_windows:
            clean_dir(os.path.join(file_dir, "fortran", "WinDLL"), rmdir=True)
        else:
            subprocess.call("make clean", shell=True, cwd=os.path.join(file_dir, "fortran"))


class BDistWheelNonPure(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, plat = super().get_tag()
        return "py3", "none", plat


class InstallPlatlib(install):
    def finalize_options(self):
        super().finalize_options()
        if self.distribution.has_ext_modules():
            self.install_lib = self.install_platlib


class BuildExtCommand(build_ext):
    """Ensure built extensions are added to the correct path in the wheel."""

    def run(self):
        pass


if __name__ == "__main__":
    setup(
        name=os.getenv("CAMB_PACKAGE_NAME", "camb"),
        zip_safe=False,
        cmdclass={
            "build_py": SharedLibrary,
            "build_cluster": SharedLibraryCluster,
            "make": MakeLibrary,
            "make_cluster": MakeLibraryCluster,
            "clean": CleanLibrary,
            "develop": DevelopLibrary,
            "develop_cluster": DevelopLibraryCluster,
            "bdist_wheel": BDistWheelNonPure,
            "install": InstallPlatlib,
            "build_ext": BuildExtCommand,
        },
        ext_modules=[Extension("camb.camblib", [])],
        packages=["camb", "camb.tests"],
        platforms="any",
        package_data={
            "camb": [
                DLLNAME,
                "libomp.dll",  # only present for Windows builds with COMPILER=flang; harmless if absent
                "HighLExtrapTemplate_lenspotentialCls.dat",
                "PArthENoPE_880.2_marcucci.dat",
                "PArthENoPE_880.2_standard.dat",
                "PRIMAT_Yp_DH_Error.dat",
                "PRIMAT_Yp_DH_ErrorMC_2021.dat",
                "PRIMAT_Yp_DH_ErrorMC_2024.dat",
            ]
        },
    )
