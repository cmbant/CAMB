from ctypes import c_int, c_double, c_bool
from .baseconfig import dll_import

# ---Variables in modules.f90
# To set the value please just put 
# variable_name.value = new_value

lensing_method_curv_corr = 1
lensing_method_flat_corr = 2
lensing_method_harmonic = 3

lensing_method = dll_import(c_int, "lensing", "lensing_method")
# lensing_method.value = lensing_method_curv_corr

ALens_Fiducial = dll_import(c_double, "lensing", "alens_fiducial")
# ALens_Fiducial.value = 0

lensing_includes_tensors = dll_import(c_bool, "lensing", "lensing_includes_tensors")
# lensing_includes_tensors.value = False

ALens = dll_import(c_double, "cambmain", "alens")
