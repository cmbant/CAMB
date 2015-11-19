from baseconfig import dll_import
from ctypes import c_int

# ---Parameters in halofit_ppf.f90

halofit_original = 1
halofit_bird = 2
halofit_peacock = 3
halofit_takahashi = 4
halofit_default = halofit_takahashi

halofit_version = dll_import(c_int, "nonlinear", "halofit_version")
# halofit_version.value = halofit_default
