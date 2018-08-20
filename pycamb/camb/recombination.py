from .baseconfig import CAMB_Structure, dll_import
from ctypes import c_bool, c_int, c_double

# logical
Do21cm = dll_import(c_double, "recombination", "do21cm")
# Do21cm.value = False

# logical
doTmatTspin = dll_import(c_bool, "recombination", "dotmattspin")
# doTmatTspin.value = False

recombination_saha_z = dll_import(c_double, "recombination", "recombination_saha_z")

recombination_saha_tau = dll_import(c_double, "recombination", "recombination_saha_tau")


# ---Derived Types in recombination.f90

class RecombinationParams(CAMB_Structure):
    """
    Holds parametes for the RECFAST recombination model (see recfast source for details).

    """
    _fields_ = [
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Heswitch", c_int),
        ("RECFAST_Hswitch", c_int)  # logical
    ]
