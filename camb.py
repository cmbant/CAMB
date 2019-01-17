#!/usr/bin/env python

# Python command line CAMB reading parameters from a .ini file
# an alternative to fortran binary compiled into fortran/camb using "make camb".
# To use .ini files from a python script use camb.run_ini or camb.read_ini

from camb._command_line import _run_command_line

_run_command_line()
