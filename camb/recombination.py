from ctypes import c_int, c_double, c_bool
from .baseconfig import F2003Class, fortran_class, optional_fortran_class


class RecombinationModel(F2003Class):
    """
    Abstract base class for recombination models
    """
    _fields_ = [
        ("min_a_evolve_Tm", c_double, "minimum scale factor at which to solve matter temperature "
                                      "perturbation if evolving sound speed or ionization fraction perturbations")
    ]


@fortran_class
class Recfast(RecombinationModel):
    """
    RECFAST recombination model (see recfast source for details).

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


@optional_fortran_class
class CosmoRec(RecombinationModel):
    """
    `CosmoRec <http://www.jb.man.ac.uk/~jchluba/Science/CosmoRec/CosmoRec.html>`_ recombination model.
    To use this, the library must be build with CosmoRec installed and RECOMBINATION_FILES including cosmorec
    in the Makefile.

    CosmoRec must be built with -fPIC added to the compiler flags.

    """
    _fortran_class_module_ = 'CosmoRec'
    _fortran_class_name_ = 'TCosmoRec'

    _fields_ = [
        ("runmode", c_int,
         "Default 0, with diffusion; 1: without diffusion; 2: RECFAST++, 3: RECFAST++ run with correction"),
        ("fdm", c_double, "Dark matter annihilation efficiency"),
        ("accuracy", c_double, "0-normal, 3-most accurate")
    ]


@optional_fortran_class
class HyRec(RecombinationModel):
    r"""
    `HyRec <https://github.com/nanoomlee/HYREC-2>`_ recombination model.
    To use this, the library must be build with HyRec installed and RECOMBINATION_FILES including hyrec in the Makefile.

    """
    _fortran_class_module_ = 'HyRec'
    _fortran_class_name_ = 'THyRec'
