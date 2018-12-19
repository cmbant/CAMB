from .baseconfig import CAMB_Structure
from ctypes import c_int, c_double, c_bool


class RecombinationParams(CAMB_Structure):
    """
    Holds parametes for the RECFAST recombination model (see recfast source for details).

    """
    _fields_ = [
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Heswitch", c_int),
        ("RECFAST_Hswitch", c_bool),
        ("AGauss1", c_double),
        ("AGauss2", c_double),
        ("zGauss1", c_double),
        ("zGauss2", c_double),
        ("wGauss1", c_double),
        ("wGauss2", c_double)
    ]
