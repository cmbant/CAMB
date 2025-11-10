"""

CAMB, Code for Anisotropies in the Microwave Background (https://camb.info)
Computational modules are wrapped Fortran 2003, but can be used entirely from Python.

"""

__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__url__ = "https://camb.readthedocs.io"
__version__ = "1.6.5"

from . import baseconfig

baseconfig.check_fortran_version(__version__)
from . import dark_energy, initialpower, model, nonlinear, reionization
from ._config import config
from .baseconfig import CAMBError, CAMBFortranError, CAMBParamRangeError, CAMBUnknownArgumentError, CAMBValueError
from .camb import (
    free_global_memory,
    get_age,
    get_background,
    get_matter_power_interpolator,
    get_results,
    get_transfer_functions,
    get_valid_numerical_params,
    get_zre_from_tau,
    read_ini,
    run_ini,
    set_feedback_level,
    set_params,
    set_params_cosmomc,
)
from .dark_energy import DarkEnergyFluid, DarkEnergyPPF
from .initialpower import InitialPowerLaw, SplinedInitialPower
from .mathutils import threej
from .model import CAMBparams, TransferParams
from .nonlinear import Halofit
from .reionization import ExpReionization, TanhReionization
from .results import CAMBdata, ClTransferData, MatterTransferData
