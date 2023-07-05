from . import camb, model
import numpy as np

from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline


def cl_kappa_limber(results, PK, ls, nz, chi_source, chi_source2=None):
    chi_source = np.float64(chi_source)
    if chi_source2 is None:
        chi_source2 = chi_source
    else:
        chi_source2 = np.float64(chi_source2)
        if chi_source2 < chi_source:
            chi_source, chi_source2 = chi_source2, chi_source
    chis = np.linspace(0, chi_source, nz, dtype=np.float64)
    zs = results.redshift_at_comoving_radial_distance(chis)
    dchis = (chis[2:] - chis[:-2]) / 2
    chis = chis[1:-1]
    zs = zs[1:-1]
    win = (1 / chis - 1 / chi_source) * (1 / chis - 1 / chi_source2) / chis ** 2
    cl = np.zeros(ls.shape)
    w = np.ones(chis.shape)
    for i, l in enumerate(ls):
        k = (l + 0.5) / chis
        w[:] = 1
        w[k < 1e-4] = 0
        w[k >= PK.kmax] = 0
        cl[i] = np.dot(dchis, w * PK.P(zs, k, grid=False) * win / k ** 4)
    cl *= (ls * (ls + 1)) ** 2
    return cl


def get_field_rotation_power(params, kmax=100, lmax=20000, non_linear=True, z_source=None,
                             k_per_logint=None, acc=1, lsamp=None):
    r"""
    Get field rotation power spectrum, :math:`C_L^{\omega\omega}`,
    following `arXiv:1605.05662 <https://arxiv.org/abs/1605.05662>`_. Uses the lowest Limber approximation.

    :param params: :class:`.model.CAMBparams` instance with cosmological parameters etc.
    :param kmax: maximum k (in :math:`{\rm Mpc}^{-1}` units)
    :param lmax: maximum L
    :param non_linear: include non-linear corrections
    :param z_source: redshift of source. If None, use peak of CMB visibility for CMB lensing
    :param k_per_logint: sampling to use in k
    :param acc: accuracy setting, increase to test stability
    :param lsamp: array of L values to compute output at. If not set, set to sampling good for interpolation
    :return: :math:`L`, :math:`C_L^{\omega\omega}`: the L sample values and corresponding rotation power
    """

    results = camb.get_background(params)
    if z_source:
        chi_source = results.comoving_radial_distance(z_source)
    else:
        chi_source = results.tau0 - results.tau_maxvis
        z_source = results.redshift_at_comoving_radial_distance(chi_source)

    PK = camb.get_matter_power_interpolator(params, nonlinear=non_linear,
                                            hubble_units=False, k_hunit=False, kmax=kmax, k_per_logint=k_per_logint,
                                            var1=model.Transfer_Weyl, var2=model.Transfer_Weyl, zmax=z_source)
    return get_field_rotation_power_from_PK(params, PK, chi_source, lmax, acc, lsamp)


def get_field_rotation_power_from_PK(params, PK, chi_source, lmax=20000, acc=1, lsamp=None):
    results = camb.get_background(params)
    nz = int(100 * acc)
    if lmax < 3000:
        raise ValueError('field rotation assumed lmax > 3000')
    ls = np.hstack((np.arange(2, 400, 1), np.arange(401, 2600, int(10. / acc)),
                    np.arange(2650, lmax, int(50. / acc)), np.arange(lmax, lmax + 1))).astype(np.float64)

    # get grid of C_L(chi_s,k) for different redshifts
    chimaxs = np.linspace(0, chi_source, nz)
    cls = np.zeros((nz, ls.size))
    for i, chimax in enumerate(chimaxs[1:]):
        cl = cl_kappa_limber(results, PK, ls, nz, chimax)
        cls[i + 1, :] = cl
    cls[0, :] = 0
    cl_chi = RectBivariateSpline(chimaxs, ls, cls)

    # Get M(L,L') matrix
    chis = np.linspace(0, chi_source, nz, dtype=np.float64)
    zs = results.redshift_at_comoving_radial_distance(chis)
    dchis = (chis[2:] - chis[:-2]) / 2
    chis = chis[1:-1]
    zs = zs[1:-1]
    win = (1 / chis - 1 / chi_source) ** 2 / chis ** 2
    w = np.ones(chis.shape)
    cchi = cl_chi(chis, ls, grid=True)
    M = np.zeros((ls.size, ls.size))
    for i, ell in enumerate(ls):
        k = (ell + 0.5) / chis
        w[:] = 1
        w[k < 1e-4] = 0
        w[k >= PK.kmax] = 0
        cl = np.dot(dchis * w * PK.P(zs, k, grid=False) * win / k ** 4, cchi)
        M[i, :] = cl * ell ** 4  # note we don't attempt to be accurate beyond lowest Limber
    Mf = RectBivariateSpline(ls, ls, np.log(M))

    # L sampling for output
    if lsamp is None:
        lsamp = np.hstack((np.arange(2, 20, 2), np.arange(25, 200, 10 // acc), np.arange(220, 1200, 30 // acc),
                           np.arange(1300, min(lmax // 2, 2600), 150 // acc),
                           np.arange(3000, lmax // 2 + 1, 1000 // acc)))

    # Get field rotation (curl) spectrum.
    diagm = np.diag(M)
    diagmsp = InterpolatedUnivariateSpline(ls, diagm)

    def high_curl_integrand(_ll, _lp):
        _lp = _lp.astype(int)
        r2 = (np.float64(_ll) / _lp) ** 2
        return _lp * r2 * diagmsp(_lp) / np.pi

    clcurl = np.zeros(lsamp.shape)
    lsall = np.arange(2, lmax + 1, dtype=np.float64)

    for i, ll in enumerate(lsamp):

        ell = np.float64(ll)
        lmin = lsall[0]
        lpmax = min(lmax, int(max(1000, ell * 2)))
        if ll < 500:
            lcalc = lsall[0:lpmax - 2]
        else:
            # sampling in L', with denser around L~L'
            lcalc = np.hstack((lsall[0:20:4],
                               lsall[29:ll - 200:35],
                               lsall[ll - 190:ll + 210:6],
                               lsall[ll + 220:lpmax + 60:60]))

        tmps = np.zeros(lcalc.shape)
        for ix, lp in enumerate(lcalc):
            llp = int(lp)
            lp = np.float64(lp)
            if abs(ll - llp) > 200 and lp > 200:
                nphi = 2 * int(min(lp / 10 * acc, 200)) + 1
            elif ll > 2000:
                nphi = 2 * int(lp / 10 * acc) + 1
            else:
                nphi = 2 * int(lp) + 1
            dphi = 2 * np.pi / nphi
            phi = np.linspace(dphi, (nphi - 1) / 2 * dphi, (nphi - 1) // 2)  # even and don't need zero
            w = 2 * np.ones(phi.size)
            cosphi = np.cos(phi)
            lrat = lp / ell
            lfact = np.sqrt(1 + lrat ** 2 - 2 * cosphi * lrat)
            lnorm = ell * lfact
            lfact[lfact <= 0] = 1
            w[lnorm < lmin] = 0
            w[lnorm > lmax] = 0

            lnorm = np.maximum(lmin, np.minimum(lmax, lnorm))
            tmps[ix] += lp * np.dot(w, (np.sin(phi) / lfact ** 2 * (cosphi - lrat)) ** 2 *
                                    np.exp(Mf(lnorm, lp, grid=False))) * dphi

        sp = InterpolatedUnivariateSpline(lcalc, tmps)
        clcurl[i] = sp.integral(2, lpmax - 1) * 4 / (2 * np.pi) ** 2

        if lpmax < lmax:
            tail = np.sum(high_curl_integrand(ll, lsall[lpmax - 2:]))
            clcurl[i] += tail

    return lsamp, clcurl


def get_field_rotation_BB(params, lmax=None, acc=1, CMB_unit='muK', raw_cl=False, spline=True):
    r"""
    Get the B-mode power spectrum from field post-born field rotation, based on perturbative and Limber approximations.
    See `arXiv:1605.05662 <https://arxiv.org/abs/1605.05662>`_.

    :param params: :class:`.model.CAMBparams` instance with cosmological parameters etc.
    :param lmax: maximum :math:`\ell`
    :param acc: accuracy
    :param CMB_unit: units for CMB output relative to dimensionless
    :param raw_cl: return :math:`C_\ell` rather than :math:`\ell(\ell+1)C_\ell/2\pi`
    :param spline: return InterpolatedUnivariateSpline, otherwise return tuple of lists of :math:`\ell`
                   and :math:`C_\ell`
    :return: InterpolatedUnivariateSpline (or arrays of sampled :math:`\ell` and) :math:`\ell^2 C_\ell^{BB}/(2 \pi)`
             (unless raw_cl, in which case just :math:`C_\ell^{BB}`)
    """

    par_CMB = params.copy()
    lmax = (lmax or 10000) * 2
    par_CMB.set_for_lmax(lmax)
    par_CMB.WantScalars = True
    par_CMB.WantCls = True
    results = camb.get_results(par_CMB)
    CE = results.get_unlensed_scalar_cls(lmax, CMB_unit=CMB_unit, raw_cl=True)[:, 1]
    CEsp = InterpolatedUnivariateSpline(np.arange(CE.shape[0]), CE)
    lsamp, clcurl = get_field_rotation_power(params, acc=acc)
    lsamp, BB = get_field_rotation_BB_integral(lsamp, clcurl, CEsp, lmax, acc=acc, raw_cl=raw_cl)
    if spline:
        return InterpolatedUnivariateSpline(lsamp, BB)
    else:
        return lsamp, BB


def get_field_rotation_BB_integral(lsamp, clcurl, cl_E_unlensed_sp, lmax=None, lsamp_out=None, acc=1, raw_cl=False):
    CurlSp = InterpolatedUnivariateSpline(lsamp, clcurl)
    lmax = lmax or lsamp[-1]
    if lsamp_out is None:
        lsamp_out = np.array([L for L in lsamp if L <= lmax // 2])
    Bcurl = np.zeros(lsamp_out.shape)

    for i, ll in enumerate(lsamp_out):
        ell = np.float64(ll)
        for llp in range(10, lmax):
            lp = np.float64(llp)
            if abs(ll - llp) > 200 and lp > 200:
                nphi = 2 * int(min(lp / 10 * acc, 200)) + 1
            elif ll > 2000:
                nphi = 2 * int(lp / 10 * acc) + 1
            else:
                nphi = 2 * int(lp) + 1
            dphi = 2 * np.pi / nphi
            phi = np.linspace(dphi, (nphi - 1) / 2 * dphi, (nphi - 1) // 2)
            w = 2 * np.ones(phi.size)
            cosphi = np.cos(phi)
            sinphi = np.sin(phi)
            sin2phi = np.sin(2 * phi)
            lpp = np.sqrt(lp ** 2 + ell ** 2 - 2 * cosphi * ell * lp)
            w[lpp < 2] = 0
            w[lpp > lmax] = 0
            curls = CurlSp(lpp)
            dCEs = cl_E_unlensed_sp(lp) * lp * dphi
            crossterm = sinphi * ell * lp / lpp ** 2
            Bcurl[i] += np.dot(w, curls * (crossterm * sin2phi) ** 2) * dCEs

    Bcurl *= 4 / (2 * np.pi) ** 2
    if not raw_cl:
        Bcurl *= lsamp_out * (lsamp_out + 1) / (2 * np.pi)
    return lsamp_out, Bcurl
