from ctypes import c_bool, c_double, c_int

from .baseconfig import CAMBValueError, F2003Class, fortran_class, optional_fortran_class

recfast_planck = "planck"
recfast_cosmorec = "cosmorec"
recfast_hyrec = "hyrec"

recfast_default = recfast_cosmorec

recfast_approx_model_params = {
    recfast_planck: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.86,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": False,
        "AGauss1": -0.14,
        "AGauss2": 0.079,
        "zGauss1": 7.28,
        "zGauss2": 6.73,
        "wGauss1": 0.18,
        "wGauss2": 0.33,
    },
    recfast_cosmorec: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.8472367977,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": True,
        "AGauss1": -0.1395272483,
        "AGauss2": 0.0729891952,
        "zGauss1": 7.2813061282,
        "zGauss2": 6.7667038679,
        "wGauss1": 0.1638966410,
        "wGauss2": 0.2785834127,
    },
    recfast_hyrec: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.8658491370,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": False,
        "AGauss1": -0.1417072428,
        "AGauss2": 0.0733978462,
        "zGauss1": 7.2841809376,
        "zGauss2": 6.7316088397,
        "wGauss1": 0.1734480042,
        "wGauss2": 0.3373321220,
    },
}


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

    def set_params(self, **kwargs):
        pass


@fortran_class
class Recfast(RecombinationModel):
    """
    RECFAST recombination model (see recfast source for details).

    """

    _fields_ = (
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Hswitch", c_bool),
        ("RECFAST_He_rate_correction", c_bool),
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

    def set_params(self, recfast_approx_model=None):
        """
        Set RECFAST approximation parameters.

        :param recfast_approx_model: optional named approximation parameter set. One of

            - ``planck``: CAMB v1.x/RECFAST 1.5.2 Planck-era fudge parameters.
            - ``cosmorec``: RECFAST fit with an He rate correction to direct CosmoRec ``accuracy=6`` with
              10 H shells.
            - ``hyrec``: seven-parameter RECFAST fit to direct HyRec-2.

        :return: self
        """
        if recfast_approx_model is not None:
            try:
                params = recfast_approx_model_params[recfast_approx_model]
            except KeyError as err:
                allowed = ", ".join(sorted(recfast_approx_model_params))
                raise CAMBValueError(
                    f"Unknown RECFAST recfast_approx_model '{recfast_approx_model}'. Use one of: {allowed}"
                ) from err
            for name, value in params.items():
                setattr(self, name, value)

        return self

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("RECFAST_H_fudge", self.RECFAST_fudge)
        state.write_fields(
            self,
            names=(
                "RECFAST_fudge_He",
                "RECFAST_Hswitch",
                "RECFAST_He_rate_correction",
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
        ("accuracy", c_double, "CosmoRec accuracy preset; see CosmoRec source for supported values"),
        ("diff_iteration_max", c_int, "Override CosmoRec diffusion iteration count; negative uses batch default (2)"),
        ("n_shells", c_int, "Override CosmoRec hydrogen shells; negative uses selected accuracy preset"),
        ("n_shells_hei", c_int, "Override CosmoRec helium shells; negative uses selected accuracy preset"),
        ("flag_hi_absorption", c_int, "Override CosmoRec HI absorption flag; negative uses selected accuracy preset"),
        ("n_s_2gamma", c_int, "Override CosmoRec two-photon shell limit; negative uses selected accuracy preset"),
        ("n_s_raman", c_int, "Override CosmoRec Raman shell limit; negative uses selected accuracy preset"),
    )

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.write_fields(
            self,
            names=(
                "runmode",
                "accuracy",
                "fdm",
                "diff_iteration_max",
                "n_shells",
                "n_shells_hei",
                "flag_hi_absorption",
                "n_s_2gamma",
                "n_s_raman",
            ),
            rename={
                "runmode": "cosmorec_runmode",
                "accuracy": "cosmorec_accuracy",
                "fdm": "cosmorec_fdm",
                "diff_iteration_max": "cosmorec_diff_iteration_max",
                "n_shells": "cosmorec_n_shells",
                "n_shells_hei": "cosmorec_n_shells_hei",
                "flag_hi_absorption": "cosmorec_flag_hi_absorption",
                "n_s_2gamma": "cosmorec_n_s_2gamma",
                "n_s_raman": "cosmorec_n_s_raman",
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
