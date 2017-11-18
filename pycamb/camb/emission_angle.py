"""
This module calculates the corrections to the standard lensed CMB power spectra results due to time delay and
emission angle, following `arXiv:1706.02673 <http://arxiv.org/abs/1706.02673>`_. This can be combined with the result
from the postborn module to estimate the leading corrections to the standard lensing B modes.

Corrections to T and E are negligible, and not calculated. The result for BB includes approximately contributions
from reionization, but this can optionally be turned off.
"""

from .baseconfig import needs_scipy
from . import camb, model
import numpy as np
from .bispectrum import threej


def cl_deflection_limber(results, PK, ls, nz, chi_source, emit_pow=2, lens_pow=0):
    chi_source = np.float64(chi_source)
    chis = np.linspace(0, chi_source, nz, dtype=np.float64)
    zs = results.redshift_at_comoving_radial_distance(chis)
    dchis = (chis[2:] - chis[:-2]) / 2
    chis = chis[1:-1]
    zs = zs[1:-1]
    win = 1 / (chis * chi_source) ** emit_pow
    if lens_pow:
        win *= (-(1 - chis / chi_source) / chis ** 2) ** lens_pow
    cl = np.zeros(ls.shape)
    w = np.ones(chis.shape)
    for i, l in enumerate(ls):
        k = (l + 0.5) / chis
        w[:] = 1
        w[k < 1e-4] = 0
        w[k >= PK.kmax] = 0
        cl[i] = np.dot(dchis, w * PK.P(zs, k, grid=False) * win / k ** 4)
    cl *= 4 * (ls * (ls + 1))
    return cl


def get_emission_angle_powers(camb_background, PK, chi_source, lmax=3000, acc=1, lsamp=None):
    """
    Get the power spectrum of psi_d, the potential for the emission angle, and its cross with standard lensing.
    Uses the Limber approximation (and assumes flat universe).

    :param camb_background: a CAMB results object, used for calling background functions
    :param PK: a matter power spectrum interpolator (from camb.get_matter_power_interpolator)
    :param chi_source: comoving radial distance of source in Mpc
    :param lmax: maximum ell
    :param acc: accuracy parameter
    :param lsamp: ell sampling for the result
    :return: a UnivariateSpline object containing L(L+1) C_L
    """

    needs_scipy()
    from scipy.interpolate import UnivariateSpline

    assert (isinstance(camb_background, camb.CAMBdata) and not camb_background.Params.omegak)
    nz = int(100 * acc)
    ls = lsamp or np.hstack((np.arange(2, 60, 2), np.arange(60, 400, 10),
                             np.arange(400, lmax, int(50. / acc)),
                             np.arange(lmax, lmax + 1))).astype(np.float64)
    cl_psi_d = cl_deflection_limber(camb_background, PK, ls, nz, chi_source, emit_pow=2, lens_pow=0)
    cl_psi_d_cross = cl_deflection_limber(camb_background, PK, ls, nz, chi_source, emit_pow=1, lens_pow=1)

    return UnivariateSpline(ls, cl_psi_d, s=0), UnivariateSpline(ls, cl_psi_d_cross, s=0)


def get_emission_delay_BB(params, kmax=100, lmax=3000, non_linear=True, CMB_unit='muK', raw_cl=False,
                          acc=1, lsamp=None, return_terms=False, include_reionization=True):
    """
    Get B modes from emission angle and time delay effects.
    Uses full-sky result from appendix of `arXiv:1706.02673 <http://arxiv.org/abs/1706.02673>`_

    :param params: :class:`.model.CAMBparams` instance with cosmological parameters etc.
    :param kmax: maximum k (in Mpc^{-1} units)
    :param lmax: maximum L
    :param non_linear: include non-linear corrections
    :param CMB_unit: normalization for the result
    :param raw_cl: if true return C_L, else L(L+1)C_L/2pi
    :param acc: accuracy setting, increase to test stability
    :param lsamp: array of L values to compute output at. If not set, set to sampling good for interpolation
    :param return_terms: return the three sub-terms separately rather than the total
    :param include_reionization: approximately include reionization terms by second scattering surface
    :return: C_L^{BB}, the L sample values and corresponding rotation power
    """

    needs_scipy()
    from scipy.interpolate import UnivariateSpline

    assert (not params.omegak)
    camb_background = camb.get_background(params)
    tau0 = camb_background.conformal_time(0)
    chi_source = tau0 - model.tau_maxvis.value
    z_source = camb_background.redshift_at_comoving_radial_distance(chi_source)

    PK = camb.get_matter_power_interpolator(params, nonlinear=non_linear,
                                            hubble_units=False, k_hunit=False, kmax=kmax,
                                            var1=model.Transfer_Weyl, var2=model.Transfer_Weyl, zmax=z_source)

    assert (lmax > 250)
    lsampvelcl = np.hstack(
        (np.arange(2, 20, 2), np.arange(25, 200, 20 // acc), np.arange(220, lmax, 40 // acc)))

    lmax_e = max(1500, lmax * 2)
    pars = params.copy()
    pars.set_for_lmax(lmax_e, lens_potential_accuracy=1)
    cmb = get_source_cmb_cl(pars, CMB_unit=CMB_unit)

    totautoB = np.zeros(lsampvelcl.shape)
    totBEterm = np.zeros(lsampvelcl.shape)
    totBxterm = np.zeros(lsampvelcl.shape)

    for reion in [False, True]:
        if reion:
            if not include_reionization: break
            zreion = params.get_zre()
            chi_source = tau0 - camb_background.conformal_time(zreion)
            lmax_e = 300
            tag_E = 'E2'
            tag_zeta = 'emit2'
            lstep = 1
        else:
            tag_E = 'E1'
            tag_zeta = 'emit1'
            lstep = 5

        cl_psi_d_sp, cl_psi_d_x_lens_sp = get_emission_angle_powers(camb_background, PK, chi_source, lmax_e, acc, lsamp)

        lsarr = np.arange(2, lmax_e + 1, dtype=np.float64)
        llp1 = lsarr * (lsarr + 1.)
        cdd = cl_psi_d_sp(lsarr) * (lsarr + 2) * (lsarr - 1) * (2 * lsarr + 1)
        cd = cl_psi_d_sp(lsarr) / llp1 * (2 * lsarr + 1)
        cxdd = cl_psi_d_x_lens_sp(lsarr) * (2 * lsarr + 1)
        cxd = cxdd / llp1
        cEE = cmb['%sx%s' % (tag_E, tag_E)][2:lmax_e + 1] / llp1 * (2 * lsarr + 1)  # raw CL
        cEx = cmb['%sx%s' % (tag_E, tag_zeta)][2:lmax_e + 1] * (2 * lsarr + 1)
        cExx = cEx / llp1
        czeta = cmb['%sx%s' % (tag_zeta, tag_zeta)][2:lmax_e + 1] / llp1 * (2 * lsarr + 1)

        for i, ll in enumerate(lsampvelcl):
            if reion and ll > lmax_e: break
            for llp in range(2, lmax_e, lstep):
                lp = np.float64(llp)
                wig = threej(llp, ll, 2, -2)
                minl = np.abs(llp - ll)
                if minl < 2:
                    wig = wig[2 - minl:]

                wigx = threej(llp, ll, 0, -2)
                minl = max(2, np.abs(llp - ll))
                maxl = min(lmax_e, np.abs(llp + ll))
                off = 0
                if (minl + llp + ll) % 2 == 0: off = 1
                wig2 = wig[off:maxl - minl + 1:2] ** 2
                totautoB[i] += lstep * np.dot(wig2, czeta[minl + off - 2:maxl + 1 - 2:2]) * cdd[llp - 2]
                totBEterm[i] += lstep * np.dot(wig2,
                                               cd[minl + off - 2:maxl + 1 - 2:2] - cxdd[minl + off - 2:maxl + 1 - 2:2]
                                               - (lp * (lp + 1) - ll * (ll + 1)) * cxd[minl + off - 2:maxl + 1 - 2:2]) \
                                * cEE[llp - 2]
                wigx2 = wigx[off:maxl - minl + 1:2] * wig[off:maxl - minl + 1:2]
                totBxterm[i] += lstep * (np.dot(wigx2, cExx[minl + off - 2:maxl + 1 - 2:2]) \
                                         * (2 * cd[llp - 2] - ((lp * (lp + 1) - ll * (ll + 1)) * cxd[llp - 2]))
                                         - np.dot(wigx2, cEx[minl + off - 2:maxl + 1 - 2:2]) * cxd[llp - 2]) \
                                * np.sqrt(lp * (lp + 1) * (lp + 2) * (lp - 1))

    fac = 1 / 2.  # (4 * np.pi)  [CMB CL already have 1/2pi in ]
    if not raw_cl: fac *= lsampvelcl * (lsampvelcl + 1) / (2 * np.pi)
    totautoB *= fac
    totBEterm *= fac
    totBxterm *= fac
    if return_terms:
        return UnivariateSpline(lsampvelcl, totautoB, s=0), \
               UnivariateSpline(lsampvelcl, totBEterm, s=0), \
               UnivariateSpline(lsampvelcl, totBxterm, s=0)
    else:
        return UnivariateSpline(lsampvelcl, totBxterm + totBEterm + totautoB, s=0)


def get_source_cmb_cl(params, CMB_unit='muK'):
    """
    Get the angular power spectrum of emission angle and time delay sources psi_t, psi_zeta,
    as well as the perpendicular velocity and E polarization.
    All are returned with 1 and 2 versions, for recombination and reionization respectively.
    Note that this function destroys any custom sources currently configured.

    :param params: :class:`.model.CAMBparams` instance with cosmological parameters etc.
    :param CMB_unit: scale results from dimensionless, use 'muK' for muK^2 units
    :return: dictionary of power spectra, with L(L+1)/2pi factors.
    """

    import sympy
    from sympy import diff
    from . import symbolic as cs

    assert (not params.omegak)

    angdist = cs.tau0 - cs.t
    emission_sources = {
        'vperp': -(cs.sigma + cs.v_b) * cs.visibility / cs.k / angdist,
        'emit': 15 * diff(cs.polter * cs.visibility, cs.t) / 8 / cs.k ** 2 / angdist,
        'delay': 15 * diff(cs.polter * cs.visibility, cs.t) / 8 / cs.k ** 2 / angdist ** 2 * (cs.tau0 - cs.tau_maxvis),
        'E': cs.scalar_E_source}

    sources = {}
    for key, source in list(emission_sources.items()):
        sources[key + '1'] = sympy.Piecewise((source, 1 / cs.a - 1 > 30), (0, True))  # recombination
        sources[key + '2'] = sympy.Piecewise((source, 1 / cs.a - 1 <= 30), (0, True))  # reionization

    camb.set_custom_scalar_sources(sources, source_ell_scales={'E1': 2, 'E2': 2})
    try:
        cldic = camb.get_results(params).get_cmb_unlensed_scalar_array_dict(CMB_unit=CMB_unit)
    finally:
        camb.clear_custom_scalar_sources()
    return cldic
