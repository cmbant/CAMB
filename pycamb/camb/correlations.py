"""
Functions to transform CMB angular power spectra into correlation functions (cl2corr)
and vice versa (corr2cl), and calculate lensed power spectra from unlensed ones.

The lensed power spectrum functions are not intended to replace those calculated by default when getting CAMB results,
but may be useful for tests, e.g. using different lensing potential power spectra, partially-delensed
lensing power spectra, etc.

These functions are all pure python/scipy, and operate and return cls including L(L+1)/2pi and [L(L+1)]^2/2pi factors.

A. Lewis December 2016
"""

import numpy as np
import os

if not os.environ.get('READTHEDOCS', None):
    from scipy.special import lpn as legendreP


def legendre_funcs(lmax, x, m=[0, 2], lfacs=None, lfacs2=None, lrootfacs=None):
    """
    Utility function to return array of Legendre and d_{mn} functions for all L up to lmax.
    Note that d_{mn} arrays start at L_min = max(m,n), so returned arrays are different sizes

    :param lmax: maximum L
    :param x: scalar value of cos(theta) at which to evaluate
    :param m: m values to calculate d_{m,n}, etc as relevant
    :param lfacs: optional pre-computed L(L+1) float array
    :param lfacs2: optional pre-computed (L+2)*(L-1) float array
    :param lrootfacs: optional pre-computed sqrt(lfacs*lfacs2) array
    :return: (P,dP),(d11,dm11), (d20, d22, d2m2) as requested, where P starts at L=0, but spin functions start at L=Lmin
    """
    allP, alldP = legendreP(lmax, x)
    # Polarization functions all start at L=2
    fac1 = 1 - x
    fac2 = 1 + x
    res = []
    if 0 in m: res.append((allP, alldP))

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
        (P, _), (d20, d22, d2m2) = legendre_funcs(lmax, x, [0, 2], lfacs, lfacs2, lrootfacs)
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
    """
    Get the sigma^2(x) and C_{gl,2}(x) functions from the lensing power spectrum

    :param clpp: array of [l(l+1)]^2 C_phi_phi/2/pi lensing potential power spectrum
    :param xvals: array of cos(theta) values at which to calculate correlation function.
    :param lmax: optional maximum L to use from the cls arrays
    :return: sigmasq, Cg2
    """
    if lmax is None: lmax = clpp.shape[0] - 1
    ls = np.arange(1, lmax + 1, dtype=np.float64)
    cldd = clpp[1:] / (ls * (ls + 1))
    cphil3 = (2 * ls + 1) * cldd / 2  # (2*l+1)l(l+1)/4pi C_phi_phi
    sigmasq = np.zeros(xvals.shape)
    Cg2 = np.zeros(xvals.shape)
    llp1 = ls * (ls + 1)
    for i, x in enumerate(xvals):
        P, dP = legendreP(lmax, x)
        fac1 = 1 - x
        fac2 = 1 + x
        d_11 = fac1 * dP[1:] / llp1 + P[1:]
        d_m11 = fac2 * dP[1:] / llp1 - P[1:]
        sigmasq[i] = np.dot(1 - d_11, cphil3)
        Cg2[i] = np.dot(d_m11, cphil3)
    return sigmasq, Cg2


def lensed_correlations(cls, clpp, xvals, weights=None, lmax=None, delta=False):
    """
    Get the lensed correlation function from the unlensed power spectra, evaluated at points cos(theta) = xvals.
    Use roots of Legendre polynomials (np.polynomial.legendre.leggauss) for accurate back integration with corr2cl.
    Note currently does not work at xvals=1 (can easily calculate that as special case!).

    To get the lensed cls efficiently, set weights to the integral weights for each x value, then function returns
    lensed correlations and lensed cls.

    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of astro-ph/0601594, to second order in C_{gl,2}

    :param cls: 2D array of unlensed cls(L,ix), with L starting at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
        cls should include l(l+1)/2pi factors.
    :param clpp: array of [l(l+1)]^2 C_phi_phi/2/pi lensing potential power spectrum
    :param xvals: array of cos(theta) values at which to calculate correlation function.
    :param weights: if given also return lensed Cls, otherwise just lensed correlations
    :param lmax: optional maximum L to use from the cls arrays
    :param delta: if true, calculate the difference between lensed and unlensed (default False)
    :return: 2D array of corrs[i, ix], where ix=0,1,2,3 are T, Q+U, Q-U and cross;
        if weights is not None, then return corrs, lensed_cls
    """

    if lmax is None: lmax = cls.shape[0] - 1
    xvals = np.asarray(xvals)
    ls = np.arange(0, lmax + 1, dtype=np.float64)
    corrs = np.zeros((len(xvals), 4))
    lfacs = ls * (ls + 1)
    lfacsall = lfacs.copy()
    lfacs[0] = 1
    cldd = clpp[1:] / lfacs[1:]
    cphil3 = (2 * ls[1:] + 1) * cldd / 2  # (2*l+1)l(l+1)/4pi C_phi_phi
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

    for i, x in enumerate(xvals):
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
            weight = weights[i]
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


def lensed_cls(cls, clpp, lmax=None, lmax_lensed=None, sampling_factor=1, delta_cls=False):
    """
    Get the lensed power spectra from the unlensed power spectra and the lensing potential power.
    Uses the non-perturbative curved-sky results from Eqs 9.12 and 9.16-9.18 of astro-ph/0601594, to second order in C_{gl,2}.
    Correlations are calculated for Gauss-Legendre integration, function is reasonably fast but intended to be
    reference implementation rather than optimized (getting lensed cls directly from CAMB is much faster).

    :param cls: 2D array of unlensed cls(L,ix), with L starting at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
        cls should include l(l+1)/2pi factors.
    :param clpp: array of [l(l+1)]^2 C_phi_phi/2/pi lensing potential power spectrum (zero based)
    :param lmax: optional maximum L to use from the cls arrays
    :param lmax_lensed: optional maximum L for the returned cl array (lmax_lensed <= lmax)
    :param sampling_factor: factor by which to oversample Gauss-Legendre points, default 1 (which is very accurate)
    :param delta_cls: if true, return the difference between lensed and unlensed (optional, default False)
    :return: 2D array of cls[L, ix], with L starting at zero and ix-0,1,2,3 in order TT, EE, BB, TE.
        cls include l(l+1)/2pi factors.
    """
    if lmax is None: lmax = cls.shape[0] - 1
    xvals, weights = np.polynomial.legendre.leggauss(sampling_factor * lmax + 1)
    _, lensedcls = lensed_correlations(cls, clpp, xvals, weights, lmax, delta=delta_cls)
    if lmax_lensed is not None:
        return lensedcls[:lmax_lensed + 1, :]
    else:
        return lensedcls
