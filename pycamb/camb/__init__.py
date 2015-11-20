# coding: utf8
"""

Python CAMB interface (http://camb.info)

"""
__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__status__ = "alpha"
__version__ = "0.1.0"

from .baseconfig import dll_import
from .camb import CAMBdata, MatterTransferData, get_results, get_transfer_functions, get_background, \
    get_age, get_zre_from_tau, set_z_outputs, set_feedback_level
from . import model
from . import initialpower
from . import reionization
from .model import CAMBparams, TransferParams
from .reionization import ReionizationParams
from .initialpower import InitialPowerParams
from ctypes import POINTER, c_int, c_double, c_float, c_bool

ThreadNum = dll_import(c_int, "modelparams", "threadnum")
# ThreadNum.value = 0

# logical
HighAccuracyDefault = dll_import(POINTER(c_bool), "modelparams", "highaccuracydefault")
HighAccuracyDefault.value = True

lSampleBoost = dll_import(c_double, "modelparams", "lsampleboost")
# lSampleBoost.value = 1.

AccuracyBoost = dll_import(c_double, "modelparams", "accuracyboost")
# AccuracyBoost.value = 1.

lAccuracyBoost = dll_import(c_float, "modelparams", "laccuracyboost")
# lAccuracyBoost.value = 1.

# Variables from module GaugeInterface
DoTensorNeutrinos = dll_import(c_bool, "gaugeinterface", "dotensorneutrinos")
# DoTensorNeutrinos.value = True

DoLateRadTruncation = dll_import(c_bool, "gaugeinterface", "dolateradtruncation")
# DoLateRadTruncation.value = True

Magnetic = dll_import(c_double, "gaugeinterface", "magnetic")
# Magnetic.value = 0.

vec_sig0 = dll_import(c_double, "gaugeinterface", "vec_sig0")
# vec_sig0.value = 1.
