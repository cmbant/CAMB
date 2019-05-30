#!/usr/bin/env python

# Python command line CAMB reading parameters from a .ini file
# an alternative to fortran binary compiled into fortran/camb using "make camb".
# To use .ini files from a python script use camb.run_ini or camb.read_ini
# If you have installed the camb package, you can just use "camb params.ini" without using this script.

from camb._command_line import run_command_line

run_command_line()
