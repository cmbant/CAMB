from ctypes import *
from baseconfig import camblib, CAMB_Structure

# logical
Do21cm = c_double.in_dll(camblib, "__recombination_MOD_do21cm")
# Do21cm.value = False

# logical
doTmatTspin = c_bool.in_dll(camblib, "__recombination_MOD_dotmattspin")
# doTmatTspin.value = False

recombination_saha_z = c_double.in_dll(camblib, "__recombination_MOD_recombination_saha_z")

recombination_saha_tau = c_double.in_dll(camblib, "__recombination_MOD_recombination_saha_tau")


# ---Derived Types in recombination.f90

class RecombinationParams(CAMB_Structure):
    _fields_ = [
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Heswitch", c_int),
        ("RECFAST_Hswitch", c_int)  # logical
    ]
