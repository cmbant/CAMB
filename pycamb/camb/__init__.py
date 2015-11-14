#coding: utf8
"""

python CAMB interface (http://camb.info)

"""
__author__ = "Antony Lewis"
__contact__ = "antony at cosmologist dot info"
__status__ = "alpha"
__version__ = "0.1.0"

from camb import CAMBdata, get_results, get_age, get_zre_from_tau, set_z_outputs, set_feedback_level
from model import CAMBparams, TransferParams
from reionization import ReionizationParams
from initialpower import InitialPowerParams
from baseconfig import camblib
from ctypes import *

ThreadNum = c_int.in_dll(camblib, "__modelparams_MOD_threadnum")
# ThreadNum.value = 0

# logical
HighAccuracyDefault = POINTER(c_bool).in_dll(camblib, "__modelparams_MOD_highaccuracydefault")
HighAccuracyDefault.value = True

lSampleBoost = c_double.in_dll(camblib, "__modelparams_MOD_lsampleboost")
# lSampleBoost.value = 1.

AccuracyBoost = c_double.in_dll(camblib, "__modelparams_MOD_accuracyboost")
# AccuracyBoost.value = 1.

lAccuracyBoost = c_float.in_dll(camblib, "__modelparams_MOD_laccuracyboost")
# lAccuracyBoost.value = 1.

# Variables from module GaugeInterface
DoTensorNeutrinos = c_bool.in_dll(camblib, "__gaugeinterface_MOD_dotensorneutrinos")
# DoTensorNeutrinos.value = True

DoLateRadTruncation = c_bool.in_dll(camblib, "__gaugeinterface_MOD_dolateradtruncation")
# DoLateRadTruncation.value = True

Magnetic = c_double.in_dll(camblib, "__gaugeinterface_MOD_magnetic")
# Magnetic.value = 0.

vec_sig0 = c_double.in_dll(camblib, "__gaugeinterface_MOD_vec_sig0")
# vec_sig0.value = 1.
