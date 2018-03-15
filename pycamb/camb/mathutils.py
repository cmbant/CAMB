from ctypes import c_int, c_double, POINTER
from .baseconfig import camblib
import numpy as np
from numpy.ctypeslib import ndpointer

numpy_2d = ndpointer(c_double, flags='C_CONTIGUOUS', ndim=2)
numpy_1d = ndpointer(c_double, flags='C_CONTIGUOUS')
int_arg = POINTER(c_int)

_chi2 = camblib.__handles_MOD_utils_getchisquared
_chi2.argtypes = [numpy_2d, numpy_1d, POINTER(c_int)]
_chi2.restype = c_double


def chi_squared(covinv, x):
    """
    Utility function to efficiently calculate x^T covinv x

    :param covinv: symmetric inverse covariance matrix
    :param x: vector
    :return: covinv.dot(x).dot(x), but parallelized and using symmetry
    """
    if len(x) != covinv.shape[0] or covinv.shape[0] != covinv.shape[1]:
        raise ValueError('Wrong shape in chi_squared')
    return _chi2(covinv, x, c_int(len(x)))


_3j = camblib.__amlutils_MOD_getthreejs
_3j.argtypes = [numpy_1d, int_arg, int_arg, int_arg, int_arg]


def threej(l2, l3, m2, m3):
    """
    Convenience wrapper around standard 3j function, returning array for all allowed l1 values
    :param l2: L_2
    :param l3: L_3
    :param m2: M_2
    :param m3: M_3
    :return: array of 3j from  max(abs(l2-l3),abs(m2+m3)) .. l2+l3
    """
    l1min = max(np.abs(l2 - l3), np.abs(m2 + m3))
    result = np.zeros(int(l3 + l2 - l1min + 1))
    l2in, l3in, m2in, m3in = c_int(l2), c_int(l3), c_int(m2), c_int(m3)
    _3j(result, l2in, l3in, m2in, m3in)
    return result


_gauss_legendre = camblib.__gauss_legendre
_gauss_legendre.argtypes = [numpy_1d, numpy_1d, int_arg]


def gauss_legendre(xvals, weights, npoints):
    _gauss_legendre(xvals, weights, c_int(npoints))
