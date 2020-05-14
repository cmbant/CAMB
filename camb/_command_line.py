import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import camb
from camb.baseconfig import lib_import
import argparse
from argparse import RawTextHelpFormatter
import ctypes


def run_command_line():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description='Python command line CAMB reading parameters from a .ini file.' +
                                                 '\n\nSample .ini files are provided in the source distribution, '
                                                 'e.g. see inifiles/planck_2018.ini at '
                                                 'https://github.com/cmbant/CAMB/tree/master/inifiles')
    parser.add_argument('ini_file', help='text .ini file with parameter settings')
    parser.add_argument('--validate', action='store_true',
                        help='Just validate the .ini file, don''t actually run anything')
    parser.add_argument('-V', '--version', action='version', version=camb.__version__)
    args = parser.parse_args()

    if not os.path.exists(args.ini_file):
        sys.exit('File not found: %s' % args.ini_file)

    s = ctypes.create_string_buffer(args.ini_file.encode("latin-1"))

    # Import wrapper function round fortran command line program
    CAMB_RunCommandLine = lib_import('camb', 'camb', 'CommandLineValidate' if args.validate else 'CommandLineRun')
    CAMB_RunCommandLine.argtypes = [ctypes.c_char_p, ctypes.c_long]
    CAMB_RunCommandLine(s, ctypes.c_long(len(args.ini_file)))
    if args.validate:
        print('OK')


if __name__ == "__main__":
    run_command_line()
