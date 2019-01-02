from ctypes import c_int, c_double, c_bool
from .baseconfig import F2003Class, fortran_class


class RecombinationModel(F2003Class):
    """
    Abstract base class for recombination models
    """
    _fields_ = [
        ("min_a_evolve_Tm", c_double,
         "minimum scale factor at which to solve matter temperature perturbation if evolving sound speed or ionization fraction perturbation")
    ]


@fortran_class
class Recfast(RecombinationModel):
    """
    Holds parameters for the RECFAST recombination model (see recfast source for details).

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

    _fortran_class_module_ = 'Recombination'
    _fortran_class_name_ = 'TRecfast'
