r"""
Functions to transform CMB angular power spectra into correlation functions (cl2corr)
and vice versa (corr2cl), and calculate lensed power spectra from unlensed ones.

The lensed power spectrum functions are not intended to replace those calculated by default when getting CAMB results,
but may be useful for tests, e.g. using different lensing potential power spectra, partially-delensed
lensing power spectra, etc.

These functions are all pure python/scipy, and operate and return cls including
factors :math:`\ell(\ell+1)/2\pi` (for CMB) and :math:`[L(L+1)]^2/2\pi` (for lensing).

A. Lewis December 2016
"""

import numpy as np
import os
from scipy.special import lpn as legendrep

try:
    from .mathutils import gauss_legendre
except:
    # use np.polynomial.legendre if can't load fast native (so can use module without compiling camb)
    # Fortran version is much faster than current np.polynomial
    gauss_legendre = None

_gauss_legendre_cache = {}


def _cached_gauss_legendre(npoints, cache=True):
    if cache and npoints in _gauss_legendre_cache:
        return _gauss_legendre_cache[npoints]
    else:
        if gauss_legendre is not None:
            xvals = np.empty(npoints)
            weights = np.empty(npoints)
            gauss_legendre(xvals, weights, npoints)
            xvals.flags.writeable = False
            weights.flags.writeable = False
        else:
            xvals, weights = np.polynomial.legendre.leggauss(npoints)
        if cache:
            _gauss_legendre_cache[npoints] = xvals, weights
        return xvals, weights


def legendre_funcs(lmax, x, m=(0, 2), lfacs=None, lfacs2=None, lrootfacs=None):
    r"""
    Utility function to return array of Legendre and :math:`d_{mn}` functions for all :math:`\ell` up to lmax.
    Note that :math:`d_{mn}` arrays start at :math:`\ell_{\rm min} = \max(m,n)`, so returned arrays are different sizes

    :param lmax: maximum :math:`\ell`
    :param x: scalar value of :math:`\cos(\theta)` at which to evaluate
    :param m: m values to calculate :math:`d_{m,n}`, etc as relevant
    :param lfacs: optional pre-computed :math:`\ell(\ell+1)` float array
    :param lfacs2: optional pre-computed :math:`(\ell+2)*(\ell-1)` float array
    :param lrootfacs: optional pre-computed sqrt(lfacs*lfacs2) array
    :return: :math:`(P,P'),(d_{11},d_{-1,1}), (d_{20}, d_{22}, d_{2,-2})` as requested, where P starts
             at :math:`\ell=0`, but spin functions start at :math:`\ell=\ell_{\rm min}`
    """
    allP, alldP = legendrep(lmax, x)
    # Polarization functions all start at L=2
    fac1 = 1 - x
    fac2 = 1 + x
    res = []
    if 0 in m:
        res.append((allP, alldP))

    if 1 in m:
        lfacs1 = np.arange(1, lmax + 1, dtype=np.float64)
        lfacs1 *= (1 + lfacs1)
        d11 = fac1 * alldP[1:] / lfacs1 + allP[1:]
        dm11 = fac2 * alldP[1:] / lfacs1 - allP[1:]
        res.append((d11, dm11))

    if 2 in m:
        if lfacs is None:
            ls = np.arange(2, lmax + 1, dtype=np.float64)
            lfacs = ls * (ls + 1)
            lfacs2 = (ls + 2) * (ls - 1)
            lrootfacs = np.sqrt(lfacs * lfacs2)
        P = allP[2:]
        dP = alldP[2:]

        fac = fac1 / fac2
        d22 = (((4 * x - 8) / fac2 + lfacs) * P
               + 4 * fac * (fac2 + (x - 2) / lfacs) * dP) / lfacs2
        if x > 0.998:
            # for stability use series at small angles (thanks Pavel Motloch)
            d2m2 = np.empty(lmax - 1)
            indser = int(np.sqrt((400.0 + 3 / (1 - x ** 2)) / 150)) - 1
            d2m2[indser:] = ((lfacs[indser:] - (4 * x + 8) / fac1) * P[indser:]
                             + 4 / fac * (-fac1 + (x + 2) / lfacs[indser:]) * dP[indser:]) / lfacs2[indser:]
            sin2 = 1 - x ** 2
            d2m2[:indser] = lfacs[:indser] * lfacs2[:indser] * sin2 ** 2 / 7680 * (20 + sin2 * (16 - lfacs[:indser]))
        else:
            d2m2 = ((lfacs - (4 * x + 8) / fac1) * P
                    + 4 / fac * (-fac1 + (x + 2) / lfacs) * dP) / lfacs2
        d20 = (2 * x * dP - lfacs * P) / lrootfacs
        res.append((d20, d22, d2m2))

    return res


def cl2corr(cls, xvals, lmax=None):
    r"""
    Get the correlation function from the power spectra, evaluated at points cos(theta) = xvals.
    Use roots of Legendre polynomials (np.polynomial.legendre.leggauss) for accurate back integration with corr2cl.
    Note currently does not work at xvals=1 (can easily calculate that as special case!).

    :param cls: 2D array cls(L,ix), with L (:math:`\equiv \ell`) starting at zero and ix-0,1,2,3 in
                order TT, EE, BB, TE. cls should include :math:`\ell(\ell+1)/2\pi` factors.
    :param xvals: array of :math:`\cos(\theta)` values at which to calculate correlation function.
    :param lmax: optional maximum L to use from the cls arrays
    :return: 2D array of corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross
    """

    if lmax is None:
        lmax = cls.shape[0] - 1
    xvals = np.asarray(xvals)
    ls = np.arange(0, lmax + 1, dtype=np.float64)
    corrs = np.zeros((len(xvals), 4))
    lfacs = ls * (ls + 1)
    lfacs[0] = 1
    facs = (2 * ls + 1) / (4 * np.pi) * 2 * np.pi / lfacs

    ct = facs * cls[:lmax + 1, 0]
    # For polarization, all arrays start at 2
    cp = facs[2:] * (cls[2:lmax + 1, 1] + cls[2:lmax + 1, 2])
    cm = facs[2:] * (cls[2:lmax + 1, 1] - cls[2:lmax + 1, 2])
    cc = facs[2:] * cls[2:lmax + 1, 3]
    ls = ls[2:]
    lfacs = lfacs[2:]
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    for i, x in enumerate(xvals):
        (P, _), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 2], lfacs, lfacs2, lrootfacs)
        corrs[i, 0] = np.dot(ct, P)  # T
        corrs[i, 1] = np.dot(cp, d22)  # Q+U
        corrs[i, 2] = np.dot(cm, d2m2)  # Q-U
        corrs[i, 3] = np.dot(cc, d20)  # cross

    return corrs


def gauss_legendre_correlation(cls, lmax=None, sampling_factor=1):
    r"""
    Transform power specturm cls into correlation functions evaluated at the
    roots of the Legendre polynomials for Gauss-Legendre quadrature. Returns correlation function array,
    evaluation points and weights.
    Result can be passed to corr2cl for accurate back transform.

    :param cls: 2D array cls(L,ix), with L (:math:`\equiv \ell`) starting at zero and ix=0,1,2,3 in
                order TT, EE, BB, TE. Should include :math:`\ell(\ell+1)/2\pi` factors.
    :param lmax: optional maximum L to use
    :param sampling_factor: uses Gauss-Legendre with degree lmax*sampling_factor+1
    :return: corrs, xvals, weights; corrs[i, ix] is 2D array where ix=0,1,2,3 are T, Q+U, Q-U and cross
    """

    if lmax is None:
        lmax = cls.shape[0] - 1
    xvals, weights = _cached_gauss_legendre(int(sampling_factor * lmax) + 1)
    return cl2corr(cls, xvals, lmax), xvals, weights


def corr2cl(corrs, xvals, weights, lmax):
    r"""
    Transform from correlation functions to power spectra.
    Note that using cl2corr followed by corr2cl is generally very accurate (< 1e-5 relative error) if
    xvals, weights = np.polynomial.legendre.leggauss(lmax+1)

    :param corrs: 2D array, corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross
    :param xvals: values of :math:`\cos(\theta)` at which corrs stores values
    :param weights: weights for integrating each point in xvals. Typically from np.polynomial.legendre.leggauss
    :param lmax: maximum :math:`\ell` to calculate :math:`C_\ell`
    :return: array of power spectra, cl[L, ix], where L starts at zero and ix=0,1,2,3 in order TT, EE, BB, TE.
      They include :math:`\ell(\ell+1)/2\pi` factors.
    """
    # For polarization, all arrays start at 2
    ls = np.arange(2, lmax + 1, dtype=np.float64)
    lfacs = ls * (ls + 1)
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    cls = np.zeros((lmax + 1, 4))
    for i, (x, weight) in enumerate(zip(xvals, weights)):
        (P, _), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 2], lfacs, lfacs2, lrootfacs)
        cls[:, 0] += (weight * corrs[i, 0]) * P
        T2 = (corrs[i, 1] * weight / 2) * d22
        T4 = (corrs[i, 2] * weight / 2) * d2m2
        cls[2:, 1] += T2 + T4
        cls[2:, 2] += T2 - T4
        cls[2:, 3] += (weight * corrs[i, 3]) * d20

    cls[1, :] *= 2
    cls[2:, :] = (cls[2:, :].T * lfacs).T

    return cls


def lensing_correlations(clpp, xvals, lmax=None):
    r"""
    Get the :math:`\sigma^2(x)` and :math:`C_{{\rm gl},2}(x)` functions from the lensing power spectrum

    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
    :param xvals: array of :math:`\cos(\theta)` values at which to calculate correlation function.
    :param lmax: optional maximum L to use from the clpp array
    :return: array of :math:`\sigma^2(x)`, array of :math:`C_{{\rm gl},2}(x)`
    """
    if lmax is None:
        lmax = clpp.shape[0] - 1
    ls = np.arange(1, lmax + 1, dtype=np.float64)
    cldd = clpp[1:] / (ls * (ls + 1))
    cphil3 = (2 * ls + 1) * cldd / 2  # (2*L+1)L(L+1)/4pi C_phi_phi
    sigmasq = np.zeros(xvals.shape)
    Cg2 = np.zeros(xvals.shape)
    llp1 = ls * (ls + 1)
    for i, x in enumerate(xvals):
        P, dP = legendrep(lmax, x)
        fac1 = 1 - x
        fac2 = 1 + x
        d_11 = fac1 * dP[1:] / llp1 + P[1:]
        d_m11 = fac2 * dP[1:] / llp1 - P[1:]
        sigmasq[i] = np.dot(1 - d_11, cphil3)
        Cg2[i] = np.dot(d_m11, cphil3)
    return sigmasq, Cg2


def lensing_R(clpp, lmax=None):
    r"""
    Get :math:`R \equiv \frac{1}{2} \langle |\nabla \phi|^2\rangle`

    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum
    :param lmax: optional maximum L to use from the cls arrays
    :return: R
    """
    if lmax is None:
        lmax = clpp.shape[0] - 1
    ls = np.arange(1, lmax + 1, dtype=np.float64)
    cldd = clpp[1:] / (ls * (ls + 1))
    cphil3 = (2 * ls + 1) * cldd / 4
    return np.sum(cphil3)


def lensed_correlations(cls, clpp, xvals, weights=None, lmax=None, delta=False, theta_max=None,
                        apodize_point_width=10):
    r"""
    Get the lensed correlation function from the unlensed power spectra, evaluated at
    points :math:`\cos(\theta)` = xvals.
    Use roots of Legendre polynomials (np.polynomial.legendre.leggauss) for accurate back integration with corr2cl.
    Note currently does not work at xvals=1 (can easily calculate that as special case!).

    To get the lensed cls efficiently, set weights to the integral weights for each x value, then function returns
    lensed correlations and lensed cls.

    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of
    `astro-ph/0601594 <https://arxiv.org/abs/astro-ph/0601594>`_, to second order in :math:`C_{{\rm gl},2}`

    :param cls: 2D array of unlensed cls(L,ix), with L (:math:`\equiv\ell`) starting at zero and ix=0,1,2,3 in
                order TT, EE, BB, TE. cls should include :math:`\ell(\ell+1)/2\pi` factors.
    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
    :param xvals: array of :math:`\cos(\theta)` values at which to calculate correlation function.
    :param weights: if given also return lensed :math:`C_\ell`, otherwise just lensed correlations
    :param lmax: optional maximum L to use from the cls arrays
    :param delta: if true, calculate the difference between lensed and unlensed (default False)
    :param theta_max: maximum angle (in radians) to keep in the correlation functions
    :param apodize_point_width: smoothing scale for apodization at truncation of correlation function
    :return: 2D array of corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross;
        if weights is not None, then return corrs, lensed_cls
    """

    if lmax is None:
        lmax = cls.shape[0] - 1
    xvals = np.asarray(xvals)
    ls = np.arange(0, lmax + 1, dtype=np.float64)
    lfacs = ls * (ls + 1)
    lfacsall = lfacs.copy()
    lfacs[0] = 1
    cldd = clpp[1:lmax + 1] / lfacs[1:]
    cphil3 = (2 * ls[1:] + 1) * cldd / 2  # (2*L+1)L(L+1)/4pi C_phi_phi
    facs = (2 * ls + 1) / (4 * np.pi) * 2 * np.pi / lfacs

    ct = facs * cls[:lmax + 1, 0]
    # For polarization, all arrays start at 2
    cp = facs[2:] * (cls[2:lmax + 1, 1] + cls[2:lmax + 1, 2])
    cm = facs[2:] * (cls[2:lmax + 1, 1] - cls[2:lmax + 1, 2])
    cc = facs[2:] * cls[2:lmax + 1, 3]
    ls = ls[2:]
    lfacs = lfacs[2:]
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    rootfac1 = np.sqrt(lfacs2)
    rootfac2 = np.sqrt((ls[1:] + 3) * (ls[1:] - 2))
    rootrat = lfacs2[1:] / rootfac1[1:] / rootfac2
    rootfac3 = np.sqrt((ls[2:] - 3) * (ls[2:] + 4))
    if weights is not None:
        lensedcls = np.zeros((lmax + 1, 4))
    if delta:
        delta_diff = 1
    else:
        delta_diff = 0

    if theta_max is not None:
        xmin = np.cos(theta_max)
        imin = np.searchsorted(xvals, xmin)  # assume xvals sorted
    else:
        imin = 0

    corrs = np.zeros((len(xvals[imin:]), 4))

    for i, x in enumerate(xvals[imin:]):
        (P, dP), (d11, dm11), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 1, 2], lfacs, lfacs2,
                                                                lrootfacs)
        sigma2 = np.dot(1 - d11, cphil3)
        Cg2 = np.dot(dm11, cphil3)

        c2fac = lfacsall[1:] * Cg2 / 2
        c2fac2 = c2fac[1:] ** 2
        fac = np.exp(-lfacsall * sigma2 / 2)
        difffac = fac - delta_diff
        f = ct * fac
        # T (don't really need the term second order Cg2 here, but include for consistency)
        corrs[i, 0] = np.dot(ct * difffac, P) + np.dot(f[1:], c2fac * (dm11 + c2fac * P[1:] / 4)) \
                      + np.dot(f[2:], c2fac2 * d2m2) / 4
        sinth = np.sqrt(1 - x ** 2)
        sinfac = 4 / sinth
        fac1 = 1 - x
        fac2 = 1 + x
        d1m2 = sinth / rootfac1 * (dP[2:] - 2 / fac1 * dm11[1:])
        d12 = sinth / rootfac1 * (dP[2:] - 2 / fac2 * d11[1:])
        d1m3 = (-(x + 0.5) * sinfac * d1m2[1:] / rootfac2 - rootrat * dm11[2:])
        d2m3 = (-fac2 * d2m2[1:] * sinfac - rootfac1[1:] * d1m2[1:]) / rootfac2
        d3m3 = (-(x + 1.5) * d2m3 * sinfac - rootfac1[1:] * d1m3) / rootfac2
        d13 = ((x - 0.5) * sinfac * d12[1:] / rootfac2 - rootrat * d11[2:])
        d04 = ((-lfacs[2:] + (18 * x ** 2 + 6) / sinth ** 2) * d20[2:] -
               6 * x * lfacs2[2:] * dP[4:] / lrootfacs[2:]) / (rootfac2[1:] * rootfac3)
        d2m4 = (-(6 * x + 4) / sinth * d2m3[1:] - rootfac2[1:] * d2m2[2:]) / rootfac3
        d4m4 = (-7 / 5.0 * (lfacs2[2:] - 6) * d2m2[2:] +
                12 / 5.0 * (-lfacs2[2:] + (9 * x + 26) / fac1) * d3m3[1:]) / (lfacs2[2:] - 12)
        # + (second order Cg2 terms are needed for <1% accuracy on BB)
        f = cp * fac[2:]
        corrs[i, 1] = np.dot(cp * difffac[2:], d22) + np.dot(f[1:], c2fac[2:] * d13) \
                      + (np.dot(f, c2fac2 * d22) + np.dot(f[2:], c2fac2[2:] * d04)) / 4
        # -
        f = cm * fac[2:]
        corrs[i, 2] = np.dot(cm * difffac[2:], d2m2) + (np.dot(f, c2fac[1:] * dm11[1:])
                                                        + np.dot(f[1:], c2fac[2:] * d3m3)) / 2 \
                      + (np.dot(f, c2fac2 * (2 * d2m2 + P[2:])) + np.dot(f[2:], c2fac2[2:] * d4m4)) / 8

        # cross
        f = cc * fac[2:]
        corrs[i, 3] = np.dot(cc * difffac[2:], d20) + (np.dot(f, c2fac[1:] * d11[1:])
                                                       + np.dot(f[1:], c2fac[2:] * d1m3)) / 2 \
                      + (3 * np.dot(f, c2fac2 * d20) + np.dot(f[2:], c2fac2[2:] * d2m4)) / 8
        if weights is not None:
            weight = weights[i + imin]
            if theta_max is not None and i < apodize_point_width * 4:
                weight *= 1 - np.exp(-((i + 1.) / apodize_point_width) ** 2 / 2)

            lensedcls[:, 0] += (weight * corrs[i, 0]) * P
            T2 = (corrs[i, 1] * weight / 2) * d22
            T4 = (corrs[i, 2] * weight / 2) * d2m2
            lensedcls[2:, 1] += T2 + T4
            lensedcls[2:, 2] += T2 - T4
            lensedcls[2:, 3] += (weight * corrs[i, 3]) * d20

    if weights is not None:
        lensedcls[1, :] *= 2
        lensedcls[2:, :] = (lensedcls[2:, :].T * lfacs).T
        return corrs, lensedcls
    else:
        return corrs


def lensed_cls(cls, clpp, lmax=None, lmax_lensed=None, sampling_factor=1.4, delta_cls=False,
               theta_max=np.pi / 32, apodize_point_width=10, leggaus=True, cache=True):
    r"""
    Get the lensed power spectra from the unlensed power spectra and the lensing potential power.
    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of
    `astro-ph/0601594 <https://arxiv.org/abs/astro-ph/0601594>`_, to second order in :math:`C_{{\rm gl},2}`.

    Correlations are calculated for Gauss-Legendre integration if leggaus=True; this slows it by several seconds,
    but will be must faster on subsequent calls with the same lmax*sampling_factor.
    If Gauss-Legendre is not used, sampling_factor needs to be about 2 times larger for same accuracy.

    For a reference implementation with the full integral range and no apodization set theta_max=None.

    Note that this function does not pad high :math:`\ell` with a smooth fit (like CAMB's main functions); for
    accurate results should be called with lmax high enough that input cls are effectively band limited
    (lmax >= 2500, or higher for accurate BB to small scales).
    Usually lmax truncation errors are far larger than other numerical errors for lmax<4000.
    For a faster result use get_lensed_cls_with_spectrum.

    :param cls: 2D array of unlensed cls(L,ix), with L starting at zero and ix=0,1,2,3 in order TT, EE, BB, TE.
        cls should include :math:`\ell(\ell+1)/2\pi` factors.
    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
    :param lmax: optional maximum L to use from the cls arrays
    :param lmax_lensed: optional maximum L for the returned cl array (lmax_lensed <= lmax)
    :param sampling_factor: npoints = int(sampling_factor*lmax)+1
    :param delta_cls: if true, return the difference between lensed and unlensed (optional, default False)
    :param theta_max: maximum angle (in radians) to keep in the correlation functions; default: pi/32
    :param apodize_point_width: if theta_max is set, apodize around the cut using half Gaussian of approx
        width apodize_point_width/lmax*pi
    :param leggaus: whether to use Gauss-Legendre integration (default True)
    :param cache: if leggaus = True, set cache to save the x values and weights between calls (most of the time)
    :return: 2D array of cls[L, ix], with L starting at zero and ix=0,1,2,3 in order TT, EE, BB, TE.
        cls include :math:`\ell(\ell+1)/2\pi` factors.
    """
    if lmax is None:
        lmax = cls.shape[0] - 1
    npoints = int(sampling_factor * lmax) + 1
    if leggaus:
        xvals, weights = _cached_gauss_legendre(npoints, cache)
    else:
        theta = np.arange(1, npoints + 1) * np.pi / (npoints + 1)
        xvals = np.cos(theta[::-1])
        weights = np.pi / npoints * np.sin(theta)
    _, lensedcls = lensed_correlations(cls, clpp, xvals, weights, lmax, delta=True,
                                       theta_max=theta_max,
                                       apodize_point_width=int(apodize_point_width * sampling_factor))
    if not delta_cls:
        lensedcls += cls[:lmax + 1, :]
    if lmax_lensed is not None:
        return lensedcls[:lmax_lensed + 1, :]
    else:
        return lensedcls


def lensed_cl_derivatives(cls, clpp, lmax=None, theta_max=np.pi / 32,
                          apodize_point_width=10, sampling_factor=1.4):
    r"""
    Get derivative dcl of lensed :math:`D_\ell\equiv \ell(\ell+1)C_\ell/2\pi` with respect to :math:`\log(C^{\phi}_L)`.
    To leading order (and hence not actually accurate), the lensed correction to power spectrum is
    is given by dcl[ix,:,:].dot(np.ones(clpp.shape)).

    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of
    `astro-ph/0601594 <https://arxiv.org/abs/astro-ph/0601594>`_, to second order in :math:`C_{{\rm gl},2}`

    :param cls: 2D array of unlensed cls(L,ix), with L starting at zero and ix=0,1,2,3 in order TT, EE, BB, TE.
        cls should include :math:`\ell(\ell+1)/2\pi` factors.
    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
    :param lmax: optional maximum L to use from the cls arrays
    :param theta_max: maximum angle (in radians) to keep in the correlation functions
    :param apodize_point_width: if theta_max is set, apodize around the cut using half Gaussian of approx
        width apodize_point_width/lmax*pi
    :param sampling_factor: npoints = int(sampling_factor*lmax)+1
    :return: array dCL[ix, ell, L], where ix=0,1,2,3 are T, EE, BB, TE and result
               is :math:`d[D^{\rm ix}_\ell]/ d (\log C^{\phi}_L)`

    """

    if lmax is None:
        lmax = cls.shape[0] - 1
    npoints = int(sampling_factor * lmax) + 1
    xvals, weights = _cached_gauss_legendre(npoints)

    ls = np.arange(0, lmax + 1, dtype=np.float64)
    lfacs = ls * (ls + 1)
    lfacsall = lfacs.copy()
    lfacs[0] = 1
    cldd = clpp[1:lmax + 1] / lfacs[1:]
    cphil3 = (2 * ls[1:] + 1) * cldd / 2  # (2*L+1)L(L+1)/4pi C_phi_phi
    facs = (2 * ls + 1) / (4 * np.pi) * 2 * np.pi / lfacs

    ct = facs * cls[:lmax + 1, 0]
    # For polarization, all arrays start at 2
    cp = facs[2:] * (cls[2:lmax + 1, 1] + cls[2:lmax + 1, 2])
    cm = facs[2:] * (cls[2:lmax + 1, 1] - cls[2:lmax + 1, 2])
    cc = facs[2:] * cls[2:lmax + 1, 3]
    ls = ls[2:]
    lfacs = lfacs[2:]
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    rootfac1 = np.sqrt(lfacs2)
    rootfac2 = np.sqrt((ls[1:] + 3) * (ls[1:] - 2))
    rootrat = lfacs2[1:] / rootfac1[1:] / rootfac2
    rootfac3 = np.sqrt((ls[2:] - 3) * (ls[2:] + 4))
    dcl = np.zeros((4, lmax + 1, lmax + 1))

    if theta_max is not None:
        xmin = np.cos(theta_max)
        imin = np.searchsorted(xvals, xmin)  # assume xvals sorted
    else:
        imin = 0

    corrs = np.zeros((4, lmax + 1))

    for i, x in enumerate(xvals[imin:]):
        (P, dP), (d11, dm11), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 1, 2], lfacs, lfacs2,
                                                                lrootfacs)
        sigma2 = np.dot(1 - d11, cphil3)
        dsigma2 = 1 - d11
        Cg2 = np.dot(dm11, cphil3)
        dCg2 = dm11

        c2fac = lfacsall[1:] * Cg2 / 2
        c2fac2 = c2fac[1:] ** 2
        fac = np.exp(-lfacsall * sigma2 / 2)
        f = -lfacsall / 2 * ct * fac
        orderfac = 1  # set to zero to neglect second order in cg2 (doesn't make much difference)
        # T (don't really need the term second order Cg2 here, but include for consistency)
        corr = np.dot(f[1:], P[1:] + c2fac * (dm11 + orderfac * c2fac * P[1:] / 4)) \
               + orderfac * np.dot(f[2:], c2fac2 * d2m2) / 4
        f = -f
        corr2 = np.dot(f[1:], (dm11 + orderfac * 2 * c2fac * P[1:] / 4)) + orderfac * 2 * np.dot(f[2:],
                                                                                                 c2fac[1:] * d2m2) / 4
        corrs[0, 1:] = (dsigma2 * corr + dCg2 * corr2) * cphil3
        sinth = np.sqrt(1 - x ** 2)
        sinfac = 4 / sinth
        fac1 = 1 - x
        fac2 = 1 + x
        d1m2 = sinth / rootfac1 * (dP[2:] - 2 / fac1 * dm11[1:])
        d12 = sinth / rootfac1 * (dP[2:] - 2 / fac2 * d11[1:])
        d1m3 = (-(x + 0.5) * sinfac * d1m2[1:] / rootfac2 - rootrat * dm11[2:])
        d2m3 = (-fac2 * d2m2[1:] * sinfac - rootfac1[1:] * d1m2[1:]) / rootfac2
        d3m3 = (-(x + 1.5) * d2m3 * sinfac - rootfac1[1:] * d1m3) / rootfac2
        d13 = ((x - 0.5) * sinfac * d12[1:] / rootfac2 - rootrat * d11[2:])
        d04 = ((-lfacs[2:] + (18 * x ** 2 + 6) / sinth ** 2) * d20[2:] -
               6 * x * lfacs2[2:] * dP[4:] / lrootfacs[2:]) / (rootfac2[1:] * rootfac3)
        d2m4 = (-(6 * x + 4) / sinth * d2m3[1:] - rootfac2[1:] * d2m2[2:]) / rootfac3
        d4m4 = (-7 / 5.0 * (lfacs2[2:] - 6) * d2m2[2:] +
                12 / 5.0 * (-lfacs2[2:] + (9 * x + 26) / fac1) * d3m3[1:]) / (lfacs2[2:] - 12)
        # + (second order Cg2 terms are needed for <1% accuracy on BB)
        f = -lfacsall[2:] / 2 * cp * fac[2:]
        corr = np.dot(f, d22) + np.dot(f[1:], c2fac[2:] * d13) \
               + orderfac * (np.dot(f, c2fac2 * d22) + np.dot(f[2:], c2fac2[2:] * d04)) / 4
        f = lfacsall[2:] / 2 * cp * fac[2:]
        corr2 = np.dot(f[1:], d13) + orderfac * 2 * (np.dot(f, c2fac[1:] * d22) + np.dot(f[2:], c2fac[3:] * d04)) / 4
        corrs[1, 1:] = (dsigma2 * corr + dCg2 * corr2) * cphil3
        # -
        f = -lfacsall[2:] / 2 * cm * fac[2:]
        corr = np.dot(f, d2m2) + (np.dot(f, c2fac[1:] * dm11[1:])
                                  + np.dot(f[1:], c2fac[2:] * d3m3)) / 2 \
               + orderfac * (np.dot(f, c2fac2 * (2 * d2m2 + P[2:])) + np.dot(f[2:], c2fac2[2:] * d4m4)) / 8
        f = -f
        corr2 = (np.dot(f, dm11[1:]) + np.dot(f[1:], d3m3)) / 2 \
                + orderfac * 2 * (np.dot(f, c2fac[1:] * (2 * d2m2 + P[2:])) + np.dot(f[2:], c2fac[3:] * d4m4)) / 8

        corrs[2, 1:] = (dsigma2 * corr + dCg2 * corr2) * cphil3

        # cross
        f = -lfacsall[2:] / 2 * cc * fac[2:]
        corr = np.dot(f, d20) + (np.dot(f, c2fac[1:] * d11[1:])
                                 + np.dot(f[1:], c2fac[2:] * d1m3)) / 2 \
               + orderfac * (3 * np.dot(f, c2fac2 * d20) + np.dot(f[2:], c2fac2[2:] * d2m4)) / 8
        f = -f
        corr2 = (np.dot(f, d11[1:]) + np.dot(f[1:], d1m3)) / 2 \
                + orderfac * 2 * (3 * np.dot(f, c2fac[1:] * d20) + np.dot(f[2:], c2fac[3:] * d2m4)) / 8

        corrs[3, 1:] = (dsigma2 * corr + dCg2 * corr2) * cphil3

        weight = weights[i + imin]
        if theta_max is not None and i < apodize_point_width * 4:
            weight *= 1 - np.exp(-((i + 1.) / apodize_point_width) ** 2 / 2)

        dcl[0, :, :] += np.outer(P, (weight * corrs[0, :]))
        T2 = np.outer(d22, (corrs[1, :] * weight / 2))
        T4 = np.outer(d2m2, (corrs[2, :] * weight / 2))
        dcl[1, 2:, :] += T2 + T4
        dcl[2, 2:, :] += T2 - T4
        dcl[3, 2:, :] += np.outer(d20, (weight * corrs[3, :]))

    # put into ell(ell+1)C_ell/2pi units [two pi cancels from correlation integral]
    dcl[:, 1, :] *= 2
    for i in range(dcl.shape[0]):
        dcl[i, 2:, :] = (dcl[i, 2:, :].T * lfacs).T
    return dcl


def lensed_cl_derivative_unlensed(clpp, lmax=None, theta_max=np.pi / 32,
                                  apodize_point_width=10, sampling_factor=1.4):
    r"""
    Get derivative dcl of lensed minus unlensed power :math:`D_\ell \equiv \ell(\ell+1)\Delta C_\ell/2\pi` with respect
    to :math:`\ell(\ell+1)C^{\rm unlens}_\ell/2\pi`

    The difference in power in the lensed spectrum is given by dCL[ix, :, :].dot(cl),
    where cl is the appropriate :math:`\ell(\ell+1)C^{\rm unlens}_\ell/2\pi`.

    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of
    `astro-ph/0601594 <https://arxiv.org/abs/astro-ph/0601594>`_, to second order in :math:`C_{{\rm gl},2}`

    :param clpp: array of :math:`[L(L+1)]^2 C_L^{\phi\phi}/2\pi` lensing potential power spectrum (zero based)
    :param lmax: optional maximum L to use from the clpp array
    :param theta_max: maximum angle (in radians) to keep in the correlation functions
    :param apodize_point_width: if theta_max is set, apodize around the cut using half Gaussian of approx
        width apodize_point_width/lmax*pi
    :param sampling_factor: npoints = int(sampling_factor*lmax)+1
    :return: array dCL[ix, ell, L], where ix=0,1,2,3 are TT, EE, BB, TE and result is
         :math:`d\left(\Delta D^{\rm ix}_\ell\right) / d D^{{\rm unlens},j}_L` where j[ix] are TT, EE, EE, TE

    """

    if lmax is None:
        lmax = clpp.shape[0] - 1
    npoints = int(sampling_factor * lmax) + 1
    xvals, weights = _cached_gauss_legendre(npoints)

    ls = np.arange(0, lmax + 1, dtype=np.float64)
    lfacs = ls * (ls + 1)
    lfacsall = lfacs.copy()
    lfacs[0] = 1
    cldd = clpp[1:lmax + 1] / lfacs[1:]
    cphil3 = (2 * ls[1:] + 1) * cldd / 2  # (2*L+1)L(L+1)/4pi C_phi_phi
    facs = (2 * ls + 1) / (4 * np.pi) * 2 * np.pi / lfacs

    ct = facs
    # For polarization, all arrays start at 2
    cp = facs[2:]
    cm = facs[2:]
    cc = facs[2:]
    ls = ls[2:]
    lfacs = lfacs[2:]
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    rootfac1 = np.sqrt(lfacs2)
    rootfac2 = np.sqrt((ls[1:] + 3) * (ls[1:] - 2))
    rootrat = lfacs2[1:] / rootfac1[1:] / rootfac2
    rootfac3 = np.sqrt((ls[2:] - 3) * (ls[2:] + 4))
    dcl = np.zeros((4, lmax + 1, lmax + 1))

    if theta_max is not None:
        xmin = np.cos(theta_max)
        imin = np.searchsorted(xvals, xmin)  # assume xvals sorted
    else:
        imin = 0

    corr = np.zeros((4, lmax + 1))

    for i, x in enumerate(xvals[imin:]):
        (P, dP), (d11, dm11), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 1, 2], lfacs, lfacs2,
                                                                lrootfacs)
        sigma2 = np.dot(1 - d11, cphil3)
        Cg2 = np.dot(dm11, cphil3)

        c2fac = lfacsall[1:] * Cg2 / 2
        c2fac2 = c2fac[1:] ** 2
        fac = np.exp(-lfacsall * sigma2 / 2)
        difffac = fac - 1
        f = ct * fac
        # T (don't really need the term second order Cg2 here, but include for consistency)
        corr[0, :] = ct * difffac * P
        corr[0, 1:] += f[1:] * c2fac * (dm11 + c2fac * P[1:] / 4)
        corr[0, 2:] += f[2:] * c2fac2 * d2m2 / 4

        sinth = np.sqrt(1 - x ** 2)
        sinfac = 4 / sinth
        fac1 = 1 - x
        fac2 = 1 + x
        d1m2 = sinth / rootfac1 * (dP[2:] - 2 / fac1 * dm11[1:])
        d12 = sinth / rootfac1 * (dP[2:] - 2 / fac2 * d11[1:])
        d1m3 = (-(x + 0.5) * sinfac * d1m2[1:] / rootfac2 - rootrat * dm11[2:])
        d2m3 = (-fac2 * d2m2[1:] * sinfac - rootfac1[1:] * d1m2[1:]) / rootfac2
        d3m3 = (-(x + 1.5) * d2m3 * sinfac - rootfac1[1:] * d1m3) / rootfac2
        d13 = ((x - 0.5) * sinfac * d12[1:] / rootfac2 - rootrat * d11[2:])
        d04 = ((-lfacs[2:] + (18 * x ** 2 + 6) / sinth ** 2) * d20[2:] -
               6 * x * lfacs2[2:] * dP[4:] / lrootfacs[2:]) / (rootfac2[1:] * rootfac3)
        d2m4 = (-(6 * x + 4) / sinth * d2m3[1:] - rootfac2[1:] * d2m2[2:]) / rootfac3
        d4m4 = (-7 / 5.0 * (lfacs2[2:] - 6) * d2m2[2:] +
                12 / 5.0 * (-lfacs2[2:] + (9 * x + 26) / fac1) * d3m3[1:]) / (lfacs2[2:] - 12)
        # + (second order Cg2 terms are needed for <1% accuracy on BB)

        f = cp * fac[2:]
        corr[1, 2:] = cp * difffac[2:] * d22 + f * c2fac2 * d22 / 4
        corr[1, 3:] += f[1:] * c2fac[2:] * d13
        corr[1, 4:] += f[2:] * c2fac2[2:] * d04 / 4

        # -
        f = cm * fac[2:]
        corr[2, 2:] = cm * difffac[2:] * d2m2 + f * c2fac[1:] * dm11[1:] / 2 + f * c2fac2 * (2 * d2m2 + P[2:]) / 8
        corr[2, 3:] += f[1:] * c2fac[2:] * d3m3 / 2
        corr[2, 4:] += f[2:] * c2fac2[2:] * d4m4 / 8

        # cross
        f = cc * fac[2:]
        corr[3, 2:] = cc * difffac[2:] * d20 + f * c2fac[1:] * d11[1:] / 2 + 3 / 8. * f * c2fac2 * d20
        corr[3, 3:] += f[1:] * c2fac[2:] * d1m3 / 2
        corr[3, 4:] += f[2:] * c2fac2[2:] * d2m4 / 8

        weight = weights[i + imin]
        if theta_max is not None and i < apodize_point_width * 4:
            weight *= 1 - np.exp(-((i + 1.) / apodize_point_width) ** 2 / 2)

        dcl[0, :, :] += np.outer(P, (weight * corr[0, :]))
        T2 = np.outer(d22, (corr[1, :] * weight / 2))
        T4 = np.outer(d2m2, (corr[2, :] * weight / 2))
        dcl[1, 2:, :] += T2 + T4
        dcl[2, 2:, :] += T2 - T4
        dcl[3, 2:, :] += np.outer(d20, (weight * corr[3, :]))

    # put into ell(ell+1)C_ell/2pi units [two pi cancels from correlation integral]
    dcl[:, 1, :] *= 2
    for i in range(dcl.shape[0]):
        dcl[i, 2:, :] = (dcl[i, 2:, :].T * lfacs).T
    return dcl
