from ctypes import c_int, c_double, c_char, POINTER
from .baseconfig import camblib
from numpy.ctypeslib import ndpointer
import numpy as np

# ---Parameters
Ini_max_string_len = 1024
max_bispectrum_deltas = 5


class TBispectrumParams:
    _fields_ = [
        ("do_lensing_bispectrum", c_int),  # logical
        ("do_primordial_bispectrum", c_int),  # logical
        ("nfields", c_int),
        ("Slice_Base_L", c_int),
        ("deltas", c_int * max_bispectrum_deltas),
        ("do_parity_odd", c_int),  # logical
        ("DoFisher", c_int),  # logical
        ("export_alpha_beta", c_int),  # logical
        ("FisherNoise", c_double),
        ("FisherNoisePol", c_double),
        ("FisherNoiseFwhmArcmin", c_double),
        ("FullOutputFile", c_char * Ini_max_string_len),
        ("SparseFullOutput", c_int),  # logical
    ]

int_arg = POINTER(c_int)
utils_3j = camblib.__amlutils_MOD_getthreejs
utils_3j.argtypes = [ndpointer(c_double, flags='C_CONTIGUOUS'), int_arg,int_arg,int_arg,int_arg]

def threej(l2,l3,m2,m3):
    """
    Convenience wrapper around standard 3j function, returning array for all allowed l1 values
    :param l2: L_2
    :param l3: L_3
    :param m2: M_2
    :param m3: M_3
    :return: array of 3j from  max(abs(l2-l3),abs(m2+m3)) .. l2+l3
    """
    l1min = np.max(np.abs(l2-l3),np.abs(m2+m3))
    result = np.zeros(int(l3+l2 -l1min+1))
    l2in, l3in, m2in, m3in = c_int(l2),c_int(l3),c_int(m2),c_int(m3)
    utils_3j(result, l2in, l3in, m2in, m3in)
    return result