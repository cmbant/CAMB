'''

In order to get access to functions in the shared library
you will need to use the following pattern by using the construction

camblib.__modulename_MOD_functionname(args)

set the data type for the functions
camblib.__modulename_MOD_functionname.restype = ctype

The arguments are passed by reference

'''

from ctypes import *
from baseconfig import camblib

# ---Variables in modules.f90
# To set the value please just put 
# variable_name.value = new_value

lensing_method_curv_corr = 1
lensing_method_flat_corr = 2
lensing_method_harmonic = 3

lensing_method = c_int.in_dll(camblib, "__lensing_MOD_lensing_method")
# lensing_method.value = lensing_method_curv_corr

ALens_Fiducial = c_double.in_dll(camblib, "__lensing_MOD_alens_fiducial")
# ALens_Fiducial.value = 0

lensing_includes_tensors = c_bool.in_dll(camblib, "__lensing_MOD_lensing_includes_tensors")
# lensing_includes_tensors.value = False

ALens = c_double.in_dll(camblib, "__cambmain_MOD_alens")
