# coding: utf8
"""

Python CAMB interface (http://camb.info)

"""
__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__status__ = "beta"
__version__ = "0.1.4"

from .baseconfig import dll_import
from .camb import CAMBdata, MatterTransferData, get_results, get_transfer_functions, get_background, \
    get_age, get_zre_from_tau, set_z_outputs, set_feedback_level, set_params, get_matter_power_interpolator
from . import model
from . import initialpower
from . import reionization
from .nonlinear import set_halofit_version
from .model import CAMBparams, TransferParams
from .reionization import ReionizationParams
from .initialpower import InitialPowerParams
from .bispectrum import threej
from ctypes import c_int, c_double, c_bool

ThreadNum = dll_import(c_int, "modelparams", "threadnum")
# ThreadNum.value = 0

# Variables from module GaugeInterface
DoTensorNeutrinos = dll_import(c_bool, "gaugeinterface", "dotensorneutrinos")
# DoTensorNeutrinos.value = True

Magnetic = dll_import(c_double, "gaugeinterface", "magnetic")
# Magnetic.value = 0.

vec_sig0 = dll_import(c_double, "gaugeinterface", "vec_sig0")
# vec_sig0.value = 1.
