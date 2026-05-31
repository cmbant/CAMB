from ctypes import c_bool, c_double, c_int

from .baseconfig import F2003Class, fortran_class, optional_fortran_class


class RecombinationModel(F2003Class):
    """
    Abstract base class for recombination models
    """

    _fields_ = (
        (
            "min_a_evolve_Tm",
            c_double,
            "minimum scale factor at which to solve matter temperature "
            "perturbation if evolving sound speed or ionization fraction perturbations",
        ),
    )

    def write_ini(self, state) -> None:
        state.set("recombination_model", self.__class__.__name__)


@fortran_class
class Recfast(RecombinationModel):
    """
    RECFAST recombination model (see recfast source for details).

    """

    _fields_ = (
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Heswitch", c_int),
        ("RECFAST_Hswitch", c_bool),
        ("AGauss1", c_double),
        ("AGauss2", c_double),
        ("zGauss1", c_double),
        ("zGauss2", c_double),
        ("wGauss1", c_double),
        ("wGauss2", c_double),
        ("Nz", c_int),
        ("use_rosenbrock", c_bool),
        ("rosenbrock_handoff_xH", c_double),
        ("rosenbrock_tol", c_double),
    )

    _fortran_class_module_ = "Recombination"
    _fortran_class_name_ = "TRecfast"

    def write_ini(self, state) -> None:
        super().write_ini(state)
        recfast_fudge = self.RECFAST_fudge
        if self.RECFAST_Hswitch:
            recfast_fudge += 1.14 - (1.105 + 0.02)
        state.set("RECFAST_fudge", recfast_fudge)
        state.write_fields(
            self,
            names=(
                "RECFAST_fudge_He",
                "RECFAST_Heswitch",
                "RECFAST_Hswitch",
                "AGauss1",
                "AGauss2",
                "zGauss1",
                "zGauss2",
                "wGauss1",
                "wGauss2",
                "Nz",
                "use_rosenbrock",
                "rosenbrock_handoff_xH",
                "rosenbrock_tol",
            ),
            rename={
                "Nz": "RECFAST_nz",
                "use_rosenbrock": "RECFAST_use_rosenbrock",
                "rosenbrock_handoff_xH": "RECFAST_rosenbrock_handoff_xH",
                "rosenbrock_tol": "RECFAST_rosenbrock_tol",
            },
        )


@optional_fortran_class
class CosmoRec(RecombinationModel):
    """
    `CosmoRec <https://www.jb.man.ac.uk/~jchluba/Science/CosmoRec/CosmoRec.html>`_ recombination model.
    To use this, the library must be built with CosmoRec installed and RECOMBINATION_FILES including cosmorec
    in the Makefile.

    CosmoRec must be built with -fPIC added to the compiler flags.

    """

    _fortran_class_module_ = "CosmoRec"
    _fortran_class_name_ = "TCosmoRec"

    _fields_ = (
        (
            "runmode",
            c_int,
            "Default 0, with diffusion; 1: without diffusion; 2: RECFAST++, 3: RECFAST++ run with correction",
        ),
        ("fdm", c_double, "Dark matter annihilation efficiency"),
        ("accuracy", c_double, "0-normal, 3-most accurate"),
    )

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.write_fields(
            self,
            names=("runmode", "accuracy", "fdm"),
            rename={
                "runmode": "cosmorec_runmode",
                "accuracy": "cosmorec_accuracy",
                "fdm": "cosmorec_fdm",
            },
        )


@optional_fortran_class
class HyRec(RecombinationModel):
    r"""
    `HyRec <https://github.com/nanoomlee/HYREC-2>`_ recombination model.
    To use this, the library must be built with HyRec installed and RECOMBINATION_FILES including hyrec in the Makefile.

    """

    _fortran_class_module_ = "HyRec"
    _fortran_class_name_ = "THyRec"
