from ctypes import c_int, c_double, c_bool, POINTER
from .baseconfig import camblib, numpy_1d, numpy_2d, numpy_3d
import numpy as np

_chi2 = camblib.__mathutils_MOD_getchisquared
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


int_arg = POINTER(c_int)
_3j = camblib.__mathutils_MOD_getthreejs
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


# Utils_3j_integrate(W,lmax_w, n, dopol, M, lmax)
_coupling_3j = camblib.__mathutils_MOD_integrate_3j
_coupling_3j.argtypes = [numpy_2d, POINTER(c_int), POINTER(c_int), POINTER(c_bool), numpy_3d, POINTER(c_int)]


def threej_coupling(W, lmax, pol=False):
    """
    Calculate symmetric coupling matrix for given weights W (e.g. mask power)
    :param W: !d array of Weights for each L, or array of weights
    :param lmax: lmax for the output matrix (assumed symmetric, though not in principle)
    :param pol: if pol, produce TT, TE, EE, EB couplings for three input mask weights
    :return: coupling matrix or array of matrices
    """
    if not isinstance(W, (list, tuple)):
        W = [W]
    if pol:
        n = 4
        if len(W) == 1: W = W * 3
        assert len(W) == 3
    else:
        n = len(W)
    M = np.empty((n, lmax + 1, lmax + 1))
    nW = len(W)
    lmax_w = min(2 * lmax, len(W[0]) - 1)
    for m in W[1:]:
        assert lmax_w == min(2 * lmax, len(m) - 1)
    Wmat = np.empty((nW, lmax_w + 1))
    for i, m in enumerate(W):
        Wmat[i, :] = m[:lmax_w + 1]
    _coupling_3j(Wmat, c_int(lmax_w), c_int(nW), c_bool(pol), M, c_int(lmax))
    if n == 1:
        return M[0, :, :]
    else:
        return [M[i, :, :] for i in range(n)]


_gauss_legendre = camblib.__mathutils_MOD_gauss_legendre
_gauss_legendre.argtypes = [numpy_1d, numpy_1d, int_arg]


def gauss_legendre(xvals, weights, npoints):
    _gauss_legendre(xvals, weights, c_int(npoints))
