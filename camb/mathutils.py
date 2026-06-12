"""
This module contains some fast utility functions that are useful in the same contexts as camb. They are entirely
independent of the main camb code.

"""

from ctypes import POINTER, c_bool, c_double, c_int

import numpy as np

from .baseconfig import camblib, numpy_1d, numpy_2d, numpy_3d

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
        raise ValueError("Wrong shape in chi_squared")
    return _chi2(covinv, x, c_int(len(x)))


int_arg = POINTER(c_int)
double_arg = POINTER(c_double)
_3j = camblib.__mathutils_MOD_getthreejs
_3j.argtypes = [numpy_1d, int_arg, int_arg, int_arg, int_arg]

_phi_olver = camblib.camb_getphiolver
_phi_olver.argtypes = [c_int, c_int, c_double, c_double]
_phi_olver.restype = c_double

_phi_olver_array = camblib.camb_getphiolverarray
_phi_olver_array.argtypes = [numpy_1d, c_int, c_int, c_double, numpy_1d, c_int]

_phi_recurs = camblib.camb_getphirecurs
_phi_recurs.argtypes = [c_int, c_int, c_double, c_double]
_phi_recurs.restype = c_double

_phi_recurs_array = camblib.camb_getphirecursarray
_phi_recurs_array.argtypes = [numpy_1d, c_int, c_int, c_double, numpy_1d, c_int]

_phi_derivative = camblib.camb_getphiderivative
_phi_derivative.argtypes = [c_int, c_int, c_double, c_double]
_phi_derivative.restype = c_double

_phi_first_peak_chi = camblib.camb_getphifirstpeakchi
_phi_first_peak_chi.argtypes = [c_int, c_int, c_double]
_phi_first_peak_chi.restype = c_double

_phi_first_peak_no_peak_found = camblib.camb_getphifirstpeaknopeakfound
_phi_first_peak_no_peak_found.argtypes = [c_int, c_int, c_double]
_phi_first_peak_no_peak_found.restype = c_int

_phi_first_peak_amplitude = camblib.camb_getphifirstpeakamplitude
_phi_first_peak_amplitude.argtypes = [c_int, c_int, c_double]
_phi_first_peak_amplitude.restype = c_double

_airy_ai_fast_array = camblib.__mathutils_MOD_airyaifastarray
_airy_ai_fast_array.argtypes = [numpy_1d, numpy_1d, int_arg]

_airy_fast_array = camblib.__mathutils_MOD_airyfastarray
_airy_fast_array.argtypes = [numpy_1d, numpy_1d, numpy_1d, int_arg]


def airy_ai_fast(x):
    """
    Fast Airy :math:`Ai(x)` approximation.

    Uses a fitted Fortran implementation optimized for < 1e-7 absolute accuracy.

    :param x: scalar or array-like input values
    :return: Airy :math:`Ai(x)`, with scalar or array shape matching ``x``
    """
    x_array = np.asarray(x, dtype=np.float64)
    flat = np.ascontiguousarray(x_array.reshape(-1))
    result = np.empty_like(flat)
    _airy_ai_fast_array(result, flat, c_int(flat.size))
    if x_array.ndim == 0:
        return result.item()
    return result.reshape(x_array.shape)


def airy_fast(x):
    """
    Fast Airy :math:`Ai(x)` and derivative :math:`Ai'(x)` approximation.

    Uses a fitted Fortran implementation optimized for < 1e-7 absolute accuracy.

    :param x: scalar or array-like input values
    :return: tuple ``(ai, aip)`` for Airy :math:`Ai(x)` and :math:`Ai'(x)`, with scalar or array shapes matching ``x``
    """
    x_array = np.asarray(x, dtype=np.float64)
    flat = np.ascontiguousarray(x_array.reshape(-1))
    ai = np.empty_like(flat)
    aip = np.empty_like(flat)
    _airy_fast_array(ai, aip, flat, c_int(flat.size))
    if x_array.ndim == 0:
        return ai.item(), aip.item()
    return ai.reshape(x_array.shape), aip.reshape(x_array.shape)


def _hyperspherical_bessel_dispatch(function, vector_function, L, K, nu, chi):
    _validate_hyperspherical_bessel_inputs(L, K, nu)
    chi_array = np.asarray(chi, dtype=np.float64)
    if np.any(chi_array < 0):
        raise ValueError("chi must be non-negative")
    nu_in = c_double(nu)

    if chi_array.ndim == 0:
        return function(c_int(L), c_int(K), nu_in, c_double(float(chi_array)))

    if chi_array.ndim != 1:
        raise ValueError("chi must be a scalar or 1D array")

    chi_array = np.ascontiguousarray(chi_array)
    result = np.empty_like(chi_array)
    vector_function(result, c_int(L), c_int(K), nu_in, chi_array, c_int(len(chi_array)))
    return result


def _validate_hyperspherical_bessel_inputs(L, K, nu):
    if L < 0:
        raise ValueError("Bessel function index L must be non-negative")
    if K not in (-1, 0, 1):
        raise ValueError("K must be one of -1, 0 or 1")
    if nu < 0:
        raise ValueError("nu must be non-negative")
    if K == 1:
        inu = round(nu)
        if abs(nu - inu) > 100 * np.finfo(float).eps * max(1.0, abs(nu)):
            raise ValueError("nu must be an integer mode for K=1")
        if inu < 3:
            raise ValueError("nu must be >= 3 for K=1")
        if inu <= L:
            raise ValueError("nu must be > L for K=1")


def phi_olver(L, K, nu, chi):
    r"""
    Evaluate the regular hyperspherical Bessel function :math:`\phi_L^\nu(K,\chi)` using the
    leading-order Olver map to a flat spherical Bessel function.

    Fast with peak-relative accuracy around 1e-4; use :func:`phi_recurs` for a slower high-accuracy reference.
    For ``K=0`` this returns the spherical Bessel function :math:`j_L(\nu\chi)`.
    Falls back to the recursive result where the Olver approximation may be unreliable.

    :param L: multipole index
    :param K: dimensionless curvature sign, one of -1, 0, 1
    :param nu: dimensionless radial eigenvalue; for closed models (``K=1``), an integer mode with ``nu >= 3``
        and ``nu > L``
    :param chi: non-negative scalar dimensionless radial distance or 1D array of non-negative values
    :return: scalar value or 1D array matching chi
    """

    return _hyperspherical_bessel_dispatch(_phi_olver, _phi_olver_array, L, K, nu, chi)


def phi_recurs(L, K, nu, chi):
    r"""
    Evaluate the regular hyperspherical Bessel function :math:`\phi_L^\nu(K,\chi)` by recurrence.

    Uses upward recurrence in the safe oscillatory region and Miller backward recurrence elsewhere,
    normalized by the exact low-order solution. The recurrence follows Abbott and Schaefer (1986)
    for the seed and recurrence formulae. Miller starts use the continued-fraction construction of
    Tram (2017) and Lesgourgues and Tram (2014) for ``K=0,-1``; for ``K=1`` they use either the
    finite closed-spectrum endpoint or the closed-space Gegenbauer Miller start.

    :param L: multipole index
    :param K: dimensionless curvature sign, one of -1, 0, 1
    :param nu: dimensionless radial eigenvalue; for closed models (``K=1``), an integer mode with ``nu >= 3``
        and ``nu > L``
    :param chi: non-negative scalar dimensionless radial distance or 1D array of non-negative values
    :return: scalar value or 1D array matching chi
    """

    return _hyperspherical_bessel_dispatch(_phi_recurs, _phi_recurs_array, L, K, nu, chi)


def phi_derivative(L, K, nu, chi):
    r"""
    Evaluate ``d phi_L^nu(K, chi) / d chi`` using the adjacent-order recurrence.
    """

    _validate_hyperspherical_bessel_inputs(L, K, nu)
    if chi < 0:
        raise ValueError("chi must be non-negative")
    return _phi_derivative(c_int(L), c_int(K), c_double(nu), c_double(chi))


def phi_first_peak_chi(L, K, nu, return_status=False):
    r"""
    Return the first peak position at or after the hyperspherical Bessel turning point.

    If ``return_status`` is true, also return whether no stationary peak was found before
    the search boundary, in which case the returned position is that boundary.
    """

    _validate_hyperspherical_bessel_inputs(L, K, nu)
    chi = _phi_first_peak_chi(c_int(L), c_int(K), c_double(nu))
    if return_status:
        no_peak_found = bool(_phi_first_peak_no_peak_found(c_int(L), c_int(K), c_double(nu)))
        return chi, no_peak_found
    return chi


def phi_first_peak_amplitude(L, K, nu):
    r"""
    Return ``abs(phi_recurs(L, K, nu, phi_first_peak_chi(L, K, nu)))``.

    If ``phi_first_peak_chi(..., return_status=True)`` reports no peak found, this is the
    amplitude at the search boundary rather than at a stationary point.
    """

    _validate_hyperspherical_bessel_inputs(L, K, nu)
    return _phi_first_peak_amplitude(c_int(L), c_int(K), c_double(nu))


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


def threej_pt(l1, l2, l3, m1, m2, m3):
    """
    Convenience testing function to get 3j for specific arguments.
    Normally use threej to get an array at once for same cost.

    :param l1: L_1
    :param l2: L_2
    :param l3: L_3
    :param m1: M_1
    :param m2: M_2
    :param m3: M_3
    :return: Wigner 3j (integer zero if outside triangle constraints)
    """
    if m1 + m2 + m3:
        return 0
    l1min = max(np.abs(l2 - l3), np.abs(m1))
    if l1 < l1min or l1 > l2 + l3:
        return 0
    wigner = threej(l2, l3, m2, m3)
    return wigner[l1 - l1min]


# Utils_3j_integrate(W,lmax_w, n, dopol, M, lmax)
_coupling_3j = camblib.__mathutils_MOD_integrate_3j
_coupling_3j.argtypes = [
    numpy_2d,
    POINTER(c_int),
    POINTER(c_int),
    POINTER(c_bool),
    numpy_3d,
    POINTER(c_int),
]


def threej_coupling(W, lmax, pol=False):
    r"""
    Calculate symmetric coupling matrix :math`\Xi` for given weights :math:`W_{\ell}`,
    where :math:`\langle\tilde{C}_\ell\rangle = \Xi_{\ell \ell'} (2\ell'+1) C_\ell`.
    The weights are related to the power spectrum of the mask P
    by :math:`W_\ell = (2 \ell + 1) P_\ell / 4 \pi`.
    See e.g. Eq D16 of `arxiv:0801.0554 <http://arxiv.org/abs/0801.0554>`_.

    If pol is False and W is an array of weights, produces array of temperature couplings, otherwise for pol is True
    produces set of TT, TE, EE, EB couplings (and weights must have one spectrum - for same masks - or three).

    Use :func:`scalar_coupling_matrix` or :func:`pcl_coupling_matrix` to get the coupling matrix directly from the
    mask power spectrum.

    :param W: 1d array of Weights for each L, or list of arrays of weights (zero based)
    :param lmax: lmax for the output matrix (assumed symmetric, though not in principle)
    :param pol: if pol, produce TT, TE, EE, EB couplings for three input mask weights (or one if assuming same mask)
    :return: symmetric coupling matrix or array of matrices
    """
    if not isinstance(W, (list, tuple)):
        W = [W]
    if pol:
        n = 4
        if len(W) == 1:
            W = W * 3
        assert len(W) == 3
    else:
        n = len(W)
    M = np.zeros((n, lmax + 1, lmax + 1))
    nW = len(W)
    lmax_w = min(2 * lmax, len(W[0]) - 1)
    for m in W[1:]:
        assert lmax_w == min(2 * lmax, len(m) - 1)
    Wmat = np.empty((nW, lmax_w + 1))
    for i, m in enumerate(W):
        Wmat[i, :] = m[: lmax_w + 1]
    _coupling_3j(Wmat, c_int(lmax_w), c_int(nW), c_bool(pol), M, c_int(lmax))
    if n == 1:
        return M[0, :, :]
    else:
        return [M[i, :, :] for i in range(n)]


def scalar_coupling_matrix(P, lmax):
    """
    Get scalar Pseudo-Cl coupling matrix from power spectrum of mask, or array of power masks.
    Uses multiple threads. See Eq A31 of `astro-ph/0105302 <https://arxiv.org/abs/astro-ph/0105302>`_

    :param P: power spectrum of mask, or list of mask power spectra
    :param lmax: lmax for the matrix (assumed square)
    :return: coupling matrix (square but not symmetric), or list of couplings for different masks
    """

    if not isinstance(P, (list, tuple)):
        P = [P]
    elif any(x.size != P[0].size for x in P[1:]):
        raise ValueError("Mask power spectra must have same lmax")

    lmax_power = min(P[0].size - 1, 2 * lmax)
    if lmax_power < 2 * lmax:
        print("Warning: power spectrum lmax is less than 2*lmax")

    fac = (2 * np.arange(lmax_power + 1) + 1) / 4 / np.pi
    M = threej_coupling([fac * power for power in P], lmax)
    factor = 2 * np.arange(lmax + 1) + 1
    if len(P) == 1:
        return M * factor
    else:
        return [m * factor for m in M]


def pcl_coupling_matrix(P, lmax, pol=False):
    """
    Get Pseudo-Cl coupling matrix from power spectrum of mask.
    Uses multiple threads. See Eq A31 of `astro-ph/0105302 <https://arxiv.org/abs/astro-ph/0105302>`_

    :param P: power spectrum of mask
    :param lmax: lmax for the matrix
    :param pol: whether to calculate TE, EE, BB couplings
    :return: coupling matrix (square but not symmetric), or list of TT, TE, EE, BB if pol
    """

    lmax_power = min(P.size - 1, 2 * lmax)
    if lmax_power < 2 * lmax:
        print("Warning: power spectrum lmax is less than 2*lmax")

    W = (2 * np.arange(lmax_power + 1) + 1) * P / (4 * np.pi)
    M = threej_coupling(W, lmax, pol=pol)

    factor = 2 * np.arange(lmax + 1) + 1
    if pol:
        return [mat * factor for mat in M]
    else:
        return M * factor


_gauss_legendre = camblib.__mathutils_MOD_gauss_legendre
_gauss_legendre.argtypes = [numpy_1d, numpy_1d, int_arg]


def gauss_legendre(xvals, weights, npoints):
    _gauss_legendre(xvals, weights, c_int(npoints))


_legendre_table = camblib.__mathutils_MOD_legendre_table
_legendre_table.argtypes = [numpy_1d, numpy_2d, numpy_2d, int_arg, int_arg]


def legendre_polynomials(x, lmax):
    """
    Legendre polynomials :math:`P_\\ell(x)` and derivatives :math:`dP_\\ell/dx` for all
    :math:`0\\le \\ell \\le` lmax (requires :math:`|x| < 1`).

    :param x: scalar or 1D array of x values
    :param lmax: maximum :math:`\\ell`
    :return: P, dP arrays; shape (lmax+1,) for scalar x, else (len(x), lmax+1)
    """
    xarr = np.ascontiguousarray(np.atleast_1d(x), dtype=np.float64)
    P = np.empty((len(xarr), lmax + 1))
    dP = np.empty((len(xarr), lmax + 1))
    _legendre_table(xarr, P, dP, c_int(lmax), c_int(len(xarr)))
    if np.ndim(x) == 0:
        return P[0], dP[0]
    return P, dP
