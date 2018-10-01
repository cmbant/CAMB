from __future__ import absolute_import
import os
import sys
import unittest
import numpy as np
import logging

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    import camb
from camb import model, correlations, bbn
from camb.baseconfig import CAMBParamRangeError


class CambTest(unittest.TestCase):
    def testBackground(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=68.5, ombh2=0.022, omch2=0.122, YHe=0.2453, mnu=0.07, omk=0)

        zre = camb.get_zre_from_tau(pars, 0.06)
        age = camb.get_age(pars)
        self.assertAlmostEqual(zre, 8.39, 2)
        self.assertAlmostEqual(age, 13.65, 2)

        data = camb.CAMBdata()
        bao = data.get_BAO([0.57, 0.27], pars)

        data = camb.CAMBdata()
        data.calc_background(pars)
        DA = data.angular_diameter_distance(0.57)
        H = data.hubble_parameter(0.27)
        self.assertAlmostEqual(DA, bao[0][2], 3)
        self.assertAlmostEqual(H, bao[1][1], 3)

        age2 = data.physical_time(0)
        self.assertAlmostEqual(age, age2, 4)

        data.comoving_radial_distance(0.48)
        t0 = data.conformal_time(0)
        t1 = data.conformal_time(11.5)
        t2 = data.comoving_radial_distance(11.5)
        self.assertAlmostEqual(t2, t0 - t1, 2)
        self.assertAlmostEqual(t1, 4200.78, 2)
        chistar = data.conformal_time(0) - model.tau_maxvis.value
        chis = np.linspace(0, chistar, 197)
        zs = data.redshift_at_comoving_radial_distance(chis)
        chitest = data.comoving_radial_distance(zs)
        self.assertTrue(np.sum((chitest - chis) ** 2) < 1e-3)

        theta = data.cosmomc_theta()
        self.assertAlmostEqual(theta, 0.0104759965, 5)

        derived = data.get_derived_params()
        self.assertAlmostEqual(derived['age'], age, 2)
        self.assertAlmostEqual(derived['rdrag'], 146.976, 2)

        # Test BBN consistency, base_plikHM_TT_lowTEB best fit model
        pars.set_cosmology(H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.YHe, 0.245336, 5)
        data.calc_background(pars)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.01040862, 7)
        self.assertAlmostEqual(data.get_derived_params()['kd'], 0.14055, 4)

        pars.set_cosmology(H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0,
                           bbn_predictor=bbn.BBN_table_interpolator())
        self.assertAlmostEqual(pars.YHe, 0.2453469, 5)
        self.assertAlmostEqual(pars.get_Y_p(), bbn.BBN_table_interpolator().Y_p(0.022242, 0), 5)

        # test massive sterile models as in Planck papers
        pars.set_cosmology(H0=68.0, ombh2=0.022305, omch2=0.11873, mnu=0.06, nnu=3.073, omk=0, meffsterile=0.013)
        self.assertAlmostEqual(pars.omegan * (pars.H0 / 100) ** 2, 0.00078, 5)
        self.assertAlmostEqual(pars.YHe, 0.24573, 5)
        self.assertAlmostEqual(pars.N_eff(), 3.073, 4)

        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()['age'], 13.773, 2)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104103, 6)

        # test dark energy
        pars.set_cosmology(H0=68.26, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars.set_dark_energy(w=-1.0226, dark_energy_model='fluid')

        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()['age'], 13.789, 2)
        scal = data.luminosity_distance(1.4)
        vec = data.luminosity_distance([1.2, 1.4, 0.1, 1.9])
        self.assertAlmostEqual(scal, vec[1], 5)
        pars.set_dark_energy()  # re-set defaults

        # test theta
        pars.set_cosmology(cosmomc_theta=0.0104085, H0=None, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.H0, 67.5512, 2)
        with self.assertRaises(CAMBParamRangeError):
            pars.set_cosmology(cosmomc_theta=0.0204085, H0=None, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars = camb.set_params(cosmomc_theta=0.0104077, H0=None, ombh2=0.022, omch2=0.122, w=-0.95)
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)
        pars.set_dark_energy()

    def testEvolution(self):
        redshifts = [0.4, 31.5]
        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95,
                               redshifts=redshifts, kmax=0.1)
        pars.WantCls = False

        # check transfer function routines and evolution code agree
        # Note transfer function redshifts are re-sorted in outputs
        data = camb.get_transfer_functions(pars)
        mtrans = data.get_matter_transfer_data()
        transfer_k = mtrans.transfer_z('delta_cdm', z_index=1)
        transfer_k2 = mtrans.transfer_z('delta_baryon', z_index=0)

        kh = mtrans.transfer_z('k/h', z_index=1)
        ev = data.get_redshift_evolution(mtrans.q, redshifts, ['delta_baryon', 'delta_cdm', 'delta_photon'],
                                         lAccuracyBoost=1)
        self.assertTrue(np.all(np.abs(transfer_k * kh ** 2 * (pars.H0 / 100) ** 2 / ev[:, 0, 1] - 1) < 1e-3))
        ix = 1
        self.assertAlmostEqual(transfer_k2[ix] * kh[ix] ** 2 * (pars.H0 / 100) ** 2, ev[ix, 1, 0], 4)

    def testPowers(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_dark_energy()  # re-set defaults
        pars.InitPower.set_params(ns=0.965, As=2e-9)

        self.assertAlmostEqual(pars.scalar_power(1), 1.801e-9, 4)
        self.assertAlmostEqual(pars.scalar_power([1, 1.5])[0], 1.801e-9, 4)

        pars.set_matter_power(nonlinear=True)
        self.assertEqual(pars.NonLinear, model.NonLinear_pk)
        pars.set_matter_power(redshifts=[0., 0.17, 3.1], nonlinear=False)
        data = camb.get_results(pars)

        kh, z, pk = data.get_matter_power_spectrum(1e-4, 1, 20)

        kh2, z2, pk2 = data.get_linear_matter_power_spectrum()

        s8 = data.get_sigma8()
        self.assertAlmostEqual(s8[0], 0.24686, 3)
        self.assertAlmostEqual(s8[2], 0.80044, 3)

        pars.NonLinear = model.NonLinear_both
        data.calc_power_spectra(pars)
        kh3, z3, pk3 = data.get_matter_power_spectrum(1e-4, 1, 20)
        self.assertAlmostEqual(pk[-1][-3], 51.909, 2)
        self.assertAlmostEqual(pk3[-1][-3], 57.697, 2)
        self.assertAlmostEqual(pk2[-2][-4], 53.476, 2)

        camb.set_feedback_level(0)

        PKnonlin = camb.get_matter_power_interpolator(pars, nonlinear=True)
        pars.set_matter_power(redshifts=[0, 0.09, 0.15, 0.42, 0.76, 1.5, 2.3, 5.5, 8.9], kmax=10, k_per_logint=5)
        pars.NonLinear = model.NonLinear_both
        results = camb.get_results(pars)
        kh, z, pk = results.get_nonlinear_matter_power_spectrum()
        pk_interp = PKnonlin.P(z, kh)
        self.assertTrue(np.sum((pk / pk_interp - 1) ** 2) < 0.005)
        PKnonlin2 = results.get_matter_power_interpolator(nonlinear=True, extrap_kmax=500)
        pk_interp2 = PKnonlin2.P(z, kh)
        self.assertTrue(np.sum((pk_interp / pk_interp2 - 1) ** 2) < 0.005)

        camb.set_halofit_version('mead')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars, var1='delta_cdm', var2='delta_cdm')
        self.assertAlmostEqual(pk[0][160], 824.6, delta=0.5)

        lmax = 4000
        pars.set_for_lmax(lmax)
        cls = data.get_cmb_power_spectra(pars)
        data.get_total_cls(2000)
        data.get_unlensed_scalar_cls(2500)
        data.get_tensor_cls(2000)
        cls_lensed = data.get_lensed_scalar_cls(3000)
        data.get_lens_potential_cls(2000)

        # check lensed CL against python; will only agree well for high lmax as python has no extrapolation template
        cls_lensed2 = correlations.lensed_cls(cls['unlensed_scalar'], cls['lens_potential'][:, 0], delta_cls=False)
        self.assertTrue(np.all(np.abs(cls_lensed2[2:2000, 2] / cls_lensed[2:2000, 2] - 1) < 1e-3))
        self.assertTrue(np.all(np.abs(cls_lensed2[2:3000, 0] / cls_lensed[2:3000, 0] - 1) < 1e-3))
        self.assertTrue(np.all(np.abs(cls_lensed2[2:3000, 1] / cls_lensed[2:3000, 1] - 1) < 1e-3))
        self.assertTrue(np.all(np.abs((cls_lensed2[2:3000, 3] - cls_lensed[2:3000, 3]) /
                                      np.sqrt(cls_lensed2[2:3000, 0] * cls_lensed2[2:3000, 1])) < 1e-4))

        corr, xvals, weights = correlations.gauss_legendre_correlation(cls['lensed_scalar'])
        clout = correlations.corr2cl(corr, xvals, weights, 2500)
        self.assertTrue(np.all(np.abs(clout[2:2300, 2] / cls['lensed_scalar'][2:2300, 2] - 1) < 1e-3))

    def testSymbolic(self):
        import camb.symbolic as s

        monopole_source, ISW, doppler, quadrupole_source = s.get_scalar_temperature_sources()
        temp_source = monopole_source + ISW + doppler + quadrupole_source

        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, omk=0.1)
        data = camb.get_background(pars)
        tau = np.linspace(1, 1200, 300)
        ks = [0.001, 0.05, 1]
        monopole2 = s.make_frame_invariant(s.newtonian_gauge(monopole_source), 'Newtonian')
        Delta_c_N = s.make_frame_invariant(s.Delta_c, 'Newtonian')
        Delta_c_N2 = s.make_frame_invariant(s.synchronous_gauge(Delta_c_N), 'CDM')

        ev = data.get_time_evolution(ks, tau, ['delta_photon', s.Delta_g, Delta_c_N, Delta_c_N2,
                                               monopole_source, monopole2,
                                               temp_source, 'T_source'])
        self.assertTrue(np.allclose(ev[:, :, 0], ev[:, :, 1]))
        self.assertTrue(np.allclose(ev[:, :, 2], ev[:, :, 3]))
        self.assertTrue(np.allclose(ev[:, :, 4], ev[:, :, 5]))
        self.assertTrue(np.allclose(ev[:, :, 6], ev[:, :, 7]))

        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95)
        pars.set_accuracy(lSampleBoost=2)
        try:
            camb.set_custom_scalar_sources([monopole_source + ISW + doppler + quadrupole_source,
                                            s.scalar_E_source], source_names=['T2', 'E2'],
                                           source_ell_scales={'E2': 2})
            data = camb.get_results(pars)
            dic = data.get_cmb_unlensed_scalar_array_dict(CMB_unit='muK')
            self.assertTrue(np.all(np.abs(dic['T2xT2'][2:2000] / dic['TxT'][2:2000] - 1) < 1e-3))
            self.assertTrue(np.all(np.abs(dic['TxT2'][2:2000] / dic['TxT'][2:2000] - 1) < 1e-3))
            # default interpolation errors much worse for E
            self.assertTrue(np.all(np.abs(dic['E2xE2'][10:2000] / dic['ExE'][10:2000] - 1) < 2e-3))
            self.assertTrue(np.all(np.abs(dic['E2xE'][10:2000] / dic['ExE'][10:2000] - 1) < 2e-3))
            dic1 = data.get_cmb_power_spectra(CMB_unit='muK')
            self.assertTrue(np.allclose(dic1['unlensed_scalar'][2:2000, 1], dic['ExE'][2:2000]))
        finally:
            pars.set_accuracy(lSampleBoost=1)
            camb.clear_custom_scalar_sources()

        s.internal_consistency_checks()

    def test_extra_EmissionAnglePostBorn(self):
        from camb import emission_angle, postborn
        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, tau=0.055)
        BB = emission_angle.get_emission_delay_BB(pars, lmax=3500)
        self.assertAlmostEqual(BB(80) * 2 * np.pi / 80 / 81., 1.1e-10, delta=1e-11)

        Bom = postborn.get_field_rotation_BB(pars, lmax=3500)
        self.assertAlmostEqual(Bom(100) * 2 * np.pi / 100 / 101., 1.65e-11, delta=1e-12)
