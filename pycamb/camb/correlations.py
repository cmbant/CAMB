"""
Functions to transform CMB angular power spectra into correlation functions (cl2corr)
and vice versa (corr2cl). These are all pure python/scipy.
"""

import numpy as np
import os

if not os.environ.get('READTHEDOCS', None):
    from scipy.special import lpn as legendreP


def legendre_funcs(lmax, x, lfacs=None, lfacs2=None, lrootfacs=None):
    """
    Utility function to return array of Legendre and d_{2m} functions for all L up to lmax.
    Note that d_{2m} arrays have size lmax-1, and start at L=2, Legendre P starts at 0 (size lmax+1)

    :param lmax: maximum L
    :param x: scalar value of cos(theta) at which to evaluate
    :param lfacs: optional pre-computed L(L+1) float array
    :param lfacs2: optional pre-computed (L+2)*(L-1) float array
    :param lrootfacs: optioanl pre-computed sqrt(lfacs*lfacs2) array
    :return: P, d20, d22, d2m2, where P starts at zero, but spin 2 functions start at L=2
    """
    if lfacs is None:
        ls = np.arange(2, lmax + 1, dtype=np.float64)
        lfacs = ls * (ls + 1)
        lfacs2 = (ls + 2) * (ls - 1)
        lrootfacs = np.sqrt(lfacs * lfacs2)
    allP, dP = legendreP(lmax, x)
    # Polarization functions all start at L=2
    P = allP[2:]
    dP = dP[2:]
    fac1 = 1 - x
    fac2 = 1 + x
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
    return allP, d20, d22, d2m2


def cl2corr(cls, xvals, lmax=None):
    """
    Get the correlation function from the power spectra, evaluated at points cos(theta) = xvals.
    Use roots of Legendre polynomials (np.polynomial.legendre.leggauss) for accurate back integration with corr2cl.
    Note currently does not work at xvals=1 (can easily calculate that as special case!).

    :param cls: 2D array cls(L,ix), with L starting at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
        cls should include l(l+1)/2pi factors.
    :param xvals: array of cos(theta) values at which to calculate correlation function.
    :param lmax: optional maximum L to use from the cls arrays
    :return: 2D array of corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross
    """

    if lmax is None: lmax = cls.shape[0] - 1
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
        P, d20, d22, d2m2 = legendre_funcs(lmax, x, lfacs, lfacs2, lrootfacs)
        corrs[i, 0] = np.dot(ct, P)  # T
        corrs[i, 1] = np.dot(cp, d22)  # Q+U
        corrs[i, 2] = np.dot(cm, d2m2)  # Q-U
        corrs[i, 3] = np.dot(cc, d20)  # cross

    return corrs


def gauss_legendre_correlation(cls, lmax=None, sampling_factor=1):
    """
    Transform power specturm cls into correlation functions evaluated at the
    roots of the Legendre polynomials for Gauss-Legendre quadrature. Returns correlation function array,
    evaluation points and weights.
    Result can be passed to corr2cl for accurate back transform.

    :param cls: 2D array cls(L,ix), with L starting at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
     Should include l*(l+1)/2pi factors.
    :param lmax: optional maximum L to use
    :param sampling_factor: uses Gauss-Legendre with degree lmax*sampling_factor+1
    :return: corrs, xvals, weights; corrs[i, ix] is 2D array where ix=0,1,2,3 are T, Q+U, Q-U and cross
    """

    if lmax is None: lmax = cls.shape[0] - 1
    xvals, weights = np.polynomial.legendre.leggauss(sampling_factor * lmax + 1)
    return cl2corr(cls, xvals, lmax), xvals, weights


def corr2cl(corrs, xvals, weights, lmax):
    """
    Transform from correlation functions to power spectra.
    Note that using cl2corr followed by corr2cl is generally very accurate (< 1e-5 relative error) if
    xvals, weights = np.polynomial.legendre.leggauss(lmax+1)

    :param corrs: 2D array, corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross
    :param xvals: values of cos(theta) at which corrs stores values
    :param weights: weights for integrating each point in xvals. Typically from np.polynomial.legendre.leggauss
    :param lmax: maximum L to calculate CL
    :return: array of power spectra, cl[L, ix], where L starts at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
      They include l(l+1)/2pi factors.
    """
    # For polarization, all arrays start at 2
    ls = np.arange(2, lmax + 1, dtype=np.float64)
    lfacs = ls * (ls + 1)
    lfacs2 = (ls + 2) * (ls - 1)
    lrootfacs = np.sqrt(lfacs * lfacs2)
    cls = np.zeros((lmax + 1, 4))
    for i, (x, weight) in enumerate(zip(xvals, weights)):
        P, d20, d22, d2m2 = legendre_funcs(lmax, x, lfacs, lfacs2, lrootfacs)
        cls[:, 0] += (weight * corrs[i, 0]) * P
        T2 = (corrs[i, 1] * weight / 2) * d22
        T4 = (corrs[i, 2] * weight / 2) * d2m2
        cls[2:, 1] += T2 + T4
        cls[2:, 2] += T2 - T4
        cls[2:, 3] += (weight * corrs[i, 3]) * d20

    cls[1, :] *= 2
    cls[2:,:] = (cls[2:, :].T * lfacs).T

    return cls
