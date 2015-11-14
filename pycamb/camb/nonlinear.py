from ctypes import *
from baseconfig import camblib

# ---Parameters in halofit_ppf.f90

halofit_original = 1
halofit_bird = 2
halofit_peacock = 3
halofit_takahashi = 4
halofit_default = halofit_takahashi

halofit_version = c_int.in_dll(camblib, "__nonlinear_MOD_halofit_version")
# halofit_version.value = halofit_default
