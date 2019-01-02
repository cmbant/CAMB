# coding: utf8
"""

Python CAMB interface (http://camb.info)

"""
__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__status__ = "beta"
__version__ = "0.3.0"

from .camb import CAMBdata, MatterTransferData, ClTransferData, get_results, get_transfer_functions, get_background, \
    get_age, get_zre_from_tau, set_feedback_level, set_params, get_matter_power_interpolator, \
    set_params_cosmomc
from . import model
from . import initialpower
from . import reionization
from .model import CAMBparams, TransferParams
from .reionization import TanhReionization
from .initialpower import InitialPowerLaw, SplinedInitialPower
from .bispectrum import threej
from ._config import config
