# coding: utf8
"""

Python CAMB interface (https://camb.info)

"""
__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__version__ = "1.0.0"

from .camb import get_results, get_transfer_functions, get_background, \
    get_age, get_zre_from_tau, set_feedback_level, set_params, get_matter_power_interpolator, \
    set_params_cosmomc, read_ini, run_ini
from . import results
from . import model
from . import initialpower
from . import reionization
from .model import CAMBparams, TransferParams
from .results import CAMBdata, MatterTransferData, ClTransferData
from .reionization import TanhReionization
from .initialpower import InitialPowerLaw, SplinedInitialPower
from .mathutils import threej
from ._config import config
