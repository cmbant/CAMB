import argparse
import ctypes
import os
import sys
from argparse import RawTextHelpFormatter

sys.path.insert(0, os.path.dirname(__file__))

from baseconfig import filepath_to_fortran, lib_import


def _get_version() -> str:
    init_path = os.path.join(os.path.dirname(__file__), "__init__.py")
    with open(init_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("__version__"):
                _, value = line.split("=", 1)
                return value.strip().strip("\"'")
    raise RuntimeError(f"Could not determine CAMB version from {init_path}")


def run_command_line():
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Python command line CAMB reading parameters from a .ini file."
        + "\n\nSample .ini files are provided in the source distribution, "
        "e.g. see inifiles/planck_2018.ini at "
        "https://github.com/cmbant/CAMB/tree/master/inifiles",
    )
    parser.add_argument("ini_file", help="text .ini file with parameter settings")
    parser.add_argument(
        "--validate", action="store_true", help="Just validate the .ini file, don't actually run anything"
    )
    parser.add_argument("-V", "--version", action="version", version=_get_version())
    args = parser.parse_args()

    if not os.path.exists(args.ini_file):
        sys.exit(f"File not found: {args.ini_file}")

    s, path_len = filepath_to_fortran(args.ini_file)

    # Import wrapper function round fortran command line program
    CAMB_RunCommandLine = lib_import("camb", "camb", "CommandLineValidate" if args.validate else "CommandLineRun")
    CAMB_RunCommandLine.argtypes = [ctypes.c_char_p, ctypes.c_long]
    CAMB_RunCommandLine(s, path_len)
    if args.validate:
        print("OK")


if __name__ == "__main__":
    run_command_line()
