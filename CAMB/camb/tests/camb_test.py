import os
import sys
import unittest
import platform
import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    import camb
from camb import model, correlations, bbn, dark_energy, initialpower
from camb.baseconfig import CAMBParamRangeError, CAMBValueError

fast = 'ci fast' in os.getenv("TRAVIS_COMMIT_MESSAGE", "")


class CambTest(unittest.TestCase):

    def testAssigments(self):
        ini = os.path.join(os.path.dirname(__file__), '..', 'inifiles', 'planck_2018.ini')
        if os.path.exists(ini):
            pars = camb.read_ini(ini)
            self.assertTrue(np.abs(camb.get_background(pars).cosmomc_theta() * 100 / 1.040909 - 1) < 2e-5)
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=68.5, ombh2=0.022, mnu=0, omch2=0.1)
        self.assertAlmostEqual(pars.omegam, (0.022 + 0.1) / 0.685 ** 2)
        with self.assertRaises(AttributeError):
            # noinspection PyPropertyAccess
            pars.omegam = 1
        pars.InitPower.set_params(ns=0.01)
        data = camb.CAMBdata()
        data.Params = pars
        self.assertEqual(data.Params.InitPower.ns, pars.InitPower.ns)
        d = dark_energy.DarkEnergyFluid(w=-0.95)
        pars.DarkEnergy = d
        self.assertEqual(pars.DarkEnergy.w, -0.95)
        pars.DarkEnergy = dark_energy.AxionEffectiveFluid(w_n=0.4)
        data.Params = pars
        self.assertEqual(pars.DarkEnergy.w_n, 0.4)
        pars.z_outputs = [0.1, 0.4]
        self.assertEqual(pars.z_outputs[1], 0.4)
        pars.z_outputs[0] = 0.3
        self.assertEqual(pars.z_outputs[0], 0.3)
        pars.z_outputs = pars.z_outputs
        pars.z_outputs = []
        pars.z_outputs = None
        # noinspection PyTypeChecker
        self.assertFalse(len(pars.z_outputs))
        with self.assertRaises(TypeError):
            pars.DarkEnergy = initialpower.InitialPowerLaw()
        pars.NonLinear = model.NonLinear_both
        printstr = str(pars)
        self.assertTrue('Want_CMB_lensing = True' in printstr and
                        "NonLinear = NonLinear_both" in printstr)
        pars.NonLinear = model.NonLinear_lens
        self.assertTrue(pars.NonLinear == model.NonLinear_lens)
        with self.assertRaises(ValueError):
            pars.NonLinear = 4
        pars.nu_mass_degeneracies = np.zeros(3)
        self.assertTrue(len(pars.nu_mass_degeneracies) == 3)
        pars.nu_mass_degeneracies = [1, 2, 3]
        self.assertTrue(pars.nu_mass_degeneracies[1] == 2)
        pars.nu_mass_degeneracies[1] = 5
        self.assertTrue(pars.nu_mass_degeneracies[1] == 5)
        with self.assertRaises(CAMBParamRangeError):
            pars.nu_mass_degeneracies = np.zeros(7)
        pars.nu_mass_eigenstates = 0
        self.assertFalse(len((pars.nu_mass_degeneracies[:1])))
        pars = camb.set_params(**{'InitPower.ns': 1.2, 'WantTransfer': True})
        self.assertEqual(pars.InitPower.ns, 1.2)
        self.assertTrue(pars.WantTransfer)
        pars.DarkEnergy = None
        pars = camb.set_params(**{'H0': 67, 'ombh2': 0.002, 'r': 0.1, 'Accuracy.AccurateBB': True})
        self.assertEqual(pars.Accuracy.AccurateBB, True)

        from camb.sources import GaussianSourceWindow
        pars = camb.CAMBparams()
        pars.SourceWindows = [GaussianSourceWindow(), GaussianSourceWindow(redshift=1)]
        self.assertEqual(pars.SourceWindows[1].redshift, 1)
        pars.SourceWindows[0].redshift = 2
        self.assertEqual(pars.SourceWindows[0].redshift, 2)
        self.assertTrue(len(pars.SourceWindows) == 2)
        pars.SourceWindows[0] = GaussianSourceWindow(redshift=3)
        self.assertEqual(pars.SourceWindows[0].redshift, 3)
        self.assertTrue('redshift = 3.0' in str(pars))
        pars.SourceWindows = pars.SourceWindows[0:1]
        self.assertTrue(len(pars.SourceWindows) == 1)
        pars.SourceWindows = []
        self.assertTrue(len(pars.SourceWindows) == 0)
        params = camb.get_valid_numerical_params()
        self.assertEqual(params, {'ombh2', 'deltazrei', 'omnuh2', 'tau', 'omk', 'zrei', 'thetastar', 'nrunrun',
                                  'meffsterile', 'nnu', 'ntrun', 'HMCode_A_baryon', 'HMCode_eta_baryon',
                                  'HMCode_logT_AGN', 'cosmomc_theta', 'YHe', 'wa', 'cs2', 'H0', 'mnu', 'Alens',
                                  'TCMB', 'ns', 'nrun', 'As', 'nt', 'r', 'w', 'omch2'})
        params2 = camb.get_valid_numerical_params(dark_energy_model='AxionEffectiveFluid')
        self.assertEqual(params2.difference(params), {'fde_zc', 'w_n', 'zc', 'theta_i'})

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
        self.assertAlmostEqual(t0, data.tau0)
        t1 = data.conformal_time(11.5)
        t2 = data.comoving_radial_distance(11.5)
        self.assertAlmostEqual(t2, t0 - t1, 2)
        self.assertAlmostEqual(t1, 4200.809, 2)
        chistar = data.conformal_time(0) - data.tau_maxvis
        chis = np.linspace(0, chistar, 197)
        zs = data.redshift_at_comoving_radial_distance(chis)
        chitest = data.comoving_radial_distance(zs)
        self.assertTrue(np.sum((chitest - chis) ** 2) < 1e-3)

        theta = data.cosmomc_theta()
        self.assertAlmostEqual(theta, 0.0104759965, 5)

        derived = data.get_derived_params()
        self.assertAlmostEqual(derived['age'], age, 2)
        self.assertAlmostEqual(derived['rdrag'], 146.986, 2)
        self.assertAlmostEqual(derived['rstar'], data.sound_horizon(derived['zstar']), 2)

        # Test BBN consistency, base_plikHM_TT_lowTEB best fit model
        pars.set_cosmology(H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.YHe, 0.2458, 5)
        data.calc_background(pars)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104090741, 7)
        self.assertAlmostEqual(data.get_derived_params()['kd'], 0.14055, 4)

        pars.set_cosmology(H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0,
                           bbn_predictor=bbn.BBN_table_interpolator())
        self.assertAlmostEqual(pars.YHe, 0.2458, 5)
        self.assertAlmostEqual(pars.get_Y_p(), bbn.BBN_table_interpolator().Y_p(0.022242, 0), 5)

        # test massive sterile models as in Planck papers
        pars.set_cosmology(H0=68.0, ombh2=0.022305, omch2=0.11873, mnu=0.06, nnu=3.073, omk=0, meffsterile=0.013)
        self.assertAlmostEqual(pars.omnuh2, 0.00078, 5)
        self.assertAlmostEqual(pars.YHe, 0.246218, 5)
        self.assertAlmostEqual(pars.N_eff, 3.073, 4)

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
        pars.set_cosmology(cosmomc_theta=0.0104085, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.H0, 67.537, 2)
        with self.assertRaises(CAMBParamRangeError):
            pars.set_cosmology(cosmomc_theta=0.0204085, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars = camb.set_params(cosmomc_theta=0.0104077, ombh2=0.022, omch2=0.122, w=-0.95)
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        pars = camb.set_params(thetastar=0.010311, ombh2=0.022, omch2=0.122)
        self.assertAlmostEqual(camb.get_background(pars).get_derived_params()['thetastar'] / 100, 0.010311, 7)
        pars = camb.set_params(thetastar=0.010311, ombh2=0.022, omch2=0.122, omk=-0.05)
        self.assertAlmostEqual(camb.get_background(pars).get_derived_params()['thetastar'] / 100, 0.010311, 7)
        self.assertAlmostEqual(pars.H0, 49.70624, places=3)

        pars = camb.set_params(cosmomc_theta=0.0104077, ombh2=0.022, omch2=0.122, w=-0.95, wa=0,
                               dark_energy_model='ppf')
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        pars = camb.set_params(cosmomc_theta=0.0104077, ombh2=0.022, omch2=0.122, w=-0.95,
                               dark_energy_model='DarkEnergyFluid', initial_power_model='InitialPowerLaw')
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        with self.assertRaises(CAMBValueError):
            camb.set_params(dark_energy_model='InitialPowerLaw')
        data.calc_background(pars)
        h2 = (data.Params.H0 / 100) ** 2
        self.assertAlmostEqual(data.get_Omega('baryon'), data.Params.ombh2 / h2, 7)
        self.assertAlmostEqual(data.get_Omega('nu'), data.Params.omnuh2 / h2, 7)
        self.assertAlmostEqual(data.get_Omega('photon') + data.get_Omega('neutrino') + data.get_Omega('de') +
                               (pars.ombh2 + pars.omch2 + pars.omnuh2) / h2 + pars.omk, 1, 8)
        pars.set_cosmology(H0=67, mnu=1, neutrino_hierarchy='normal')
        data.calc_background(pars)
        h2 = (pars.H0 / 100) ** 2
        self.assertAlmostEqual(data.get_Omega('photon') + data.get_Omega('neutrino') + data.get_Omega('de') +
                               (pars.ombh2 + pars.omch2 + pars.omnuh2) / h2 + pars.omk, 1, 8)
        redshifts = np.array([0.005, 0.01, 0.3, 0.9342, 4, 27, 321.5, 932, 1049, 1092, 2580, 1e4, 2.1e7])
        self.assertTrue(
            np.allclose(data.redshift_at_conformal_time(data.conformal_time(redshifts)), redshifts, rtol=1e-7))
        pars.set_dark_energy(w=-1.8)
        data.calc_background(pars)
        self.assertTrue(
            np.allclose(data.redshift_at_conformal_time(data.conformal_time(redshifts)), redshifts, rtol=1e-7))
        pars.set_cosmology(cosmomc_theta=0.0104085)
        data.calc_background(pars)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104085)
        derived = data.get_derived_params()
        pars.Accuracy.BackgroundTimeStepBoost = 2
        data.calc_background(pars)
        derived2 = data.get_derived_params()
        self.assertAlmostEqual(derived['thetastar'], derived2['thetastar'], places=5)
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.11, neutrino_hierarchy='inverted')
        self.assertEqual(pars.num_nu_massive, 3)
        self.assertEqual(pars.nu_mass_numbers[1], 1)
        self.assertEqual(pars.nu_mass_eigenstates, 2)
        self.assertAlmostEqual(pars.nu_mass_fractions[0], 0.915197, places=4)

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=68.5, ombh2=0.022, omch2=0.122, YHe=0.2453, mnu=0.07, omk=0, zrei=zre)
        results = camb.get_background(pars)
        self.assertEqual(results.Params.Reion.redshift, zre)

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=68.5, ombh2=0.022, omch2=0.122, YHe=0.2453, mnu=0.07, omk=-0.05)
        data = camb.get_background(pars)
        delta2 = data.curvature_radius / (1 + 0.25) * (
            np.sin((data.comoving_radial_distance(0.25) - data.comoving_radial_distance(0.05)) / data.curvature_radius))
        np.testing.assert_allclose(delta2, data.angular_diameter_distance2(0.05, 0.25))
        dists = data.angular_diameter_distance2([0.3, 0.05, 0.25], [1, 0.25, 0.05])
        self.assertAlmostEqual(delta2, dists[1])
        self.assertEqual(0, dists[2])

        self.assertEqual(data.physical_time(0.4), data.physical_time([0.2, 0.4])[1])
        d = data.conformal_time_a1_a2(0, 0.5) + data.conformal_time_a1_a2(0.5, 1)
        self.assertAlmostEqual(d, data.conformal_time_a1_a2(0, 1))
        self.assertAlmostEqual(d, sum(data.conformal_time_a1_a2([0, 0.5], [0.5, 1])))

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

    def testInstances(self):
        pars = camb.set_params(H0=69.1, ombh2=0.032, omch2=0.122, As=3e-9, ns=0.91, omk=0.013,
                               redshifts=[0.], kmax=0.5)
        data = camb.get_background(pars)
        res1 = data.angular_diameter_distance(0.7)
        drag1 = data.get_derived_params()['rdrag']
        pars2 = camb.set_params(H0=65, ombh2=0.022, omch2=0.122, As=3e-9, ns=0.91)
        data2 = camb.get_background(pars2)
        res2 = data2.angular_diameter_distance(1.7)
        drag2 = data2.get_derived_params()['rdrag']
        self.assertAlmostEqual(res1, data.angular_diameter_distance(0.7))
        self.assertAlmostEqual(res2, data2.angular_diameter_distance(1.7))
        self.assertAlmostEqual(drag1, data.get_derived_params()['rdrag'])
        self.assertEqual(pars2.InitPower.ns, data2.Params.InitPower.ns)
        data2.calc_background(pars)
        self.assertEqual(pars.InitPower.ns, data2.Params.InitPower.ns)
        self.assertAlmostEqual(res1, data2.angular_diameter_distance(0.7))
        data3 = camb.get_results(pars2)
        cl3 = data3.get_lensed_scalar_cls(1000)
        self.assertAlmostEqual(res2, data3.angular_diameter_distance(1.7))
        self.assertAlmostEqual(drag2, data3.get_derived_params()['rdrag'], places=3)
        self.assertAlmostEqual(drag1, data.get_derived_params()['rdrag'], places=3)
        pars.set_for_lmax(3000, lens_potential_accuracy=1)
        camb.get_results(pars)
        del data3
        data4 = camb.get_results(pars2)
        cl4 = data4.get_lensed_scalar_cls(1000)
        self.assertTrue(np.allclose(cl4, cl3))

    def testPowers(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_dark_energy()  # re-set defaults
        pars.InitPower.set_params(ns=0.965, As=2e-9)
        pars.NonLinearModel.set_params(halofit_version='takahashi')

        self.assertAlmostEqual(pars.scalar_power(1), 1.801e-9, 4)
        self.assertAlmostEqual(pars.scalar_power([1, 1.5])[0], 1.801e-9, 4)

        pars.set_matter_power(nonlinear=True)
        self.assertEqual(pars.NonLinear, model.NonLinear_pk)
        pars.set_matter_power(redshifts=[0., 0.17, 3.1], silent=True, nonlinear=False)
        data = camb.get_results(pars)

        kh, z, pk = data.get_matter_power_spectrum(1e-4, 1, 20)

        kh2, z2, pk2 = data.get_linear_matter_power_spectrum()

        s8 = data.get_sigma8()
        self.assertAlmostEqual(s8[0], 0.24686, 3)
        self.assertAlmostEqual(s8[2], 0.80044, 3)
        fs8 = data.get_fsigma8()
        self.assertAlmostEqual(fs8[0], 0.2431, 3)
        self.assertAlmostEqual(fs8[2], 0.424712, 3)

        pars.NonLinear = model.NonLinear_both

        data.calc_power_spectra(pars)
        kh3, z3, pk3 = data.get_matter_power_spectrum(1e-4, 1, 20)
        self.assertAlmostEqual(pk[-1][-3], 51.924, 2)
        self.assertAlmostEqual(pk3[-1][-3], 57.723, 2)
        self.assertAlmostEqual(pk2[-2][-4], 56.454, 2)
        camb.set_feedback_level(0)

        PKnonlin = camb.get_matter_power_interpolator(pars, nonlinear=True)
        pars.set_matter_power(redshifts=[0, 0.09, 0.15, 0.42, 0.76, 1.5, 2.3, 5.5, 8.9],
                              silent=True, kmax=10, k_per_logint=5)
        pars.NonLinear = model.NonLinear_both
        results = camb.get_results(pars)
        kh, z, pk = results.get_nonlinear_matter_power_spectrum()
        pk_interp = PKnonlin.P(z, kh)
        self.assertTrue(np.sum((pk / pk_interp - 1) ** 2) < 0.005)
        PKnonlin2 = results.get_matter_power_interpolator(nonlinear=True, extrap_kmax=500)
        pk_interp2 = PKnonlin2.P(z, kh)
        self.assertTrue(np.sum((pk_interp / pk_interp2 - 1) ** 2) < 0.005)

        pars.NonLinearModel.set_params(halofit_version='mead')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
        self.assertAlmostEqual(pk[0][160], 814.9, delta=0.5)

        pars.NonLinearModel.set_params(halofit_version='mead2016')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
        self.assertAlmostEqual(pk[0][160], 814.9, delta=0.5)

        pars.NonLinearModel.set_params(halofit_version='mead2015')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
        self.assertAlmostEqual(pk[0][160], 791.3, delta=0.5)

        pars.NonLinearModel.set_params(halofit_version='mead2020')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
        self.assertAlmostEqual(pk[0][160], 815.8, delta=0.5)

        pars.NonLinearModel.set_params(halofit_version='mead2020_feedback')
        _, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
        self.assertAlmostEqual(pk[0][160], 799.0, delta=0.5)

        lmax = 4000
        pars.set_for_lmax(lmax)
        cls = data.get_cmb_power_spectra(pars)
        data.get_total_cls(2000)
        cls_unlensed = data.get_unlensed_scalar_cls(2500)
        data.get_tensor_cls(2000)
        cls_lensed = data.get_lensed_scalar_cls(3000)
        data.get_lens_potential_cls(2000)

        cls_lensed2 = data.get_lensed_cls_with_spectrum(data.get_lens_potential_cls()[:, 0], lmax=3000)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_lensed[2:, :], rtol=1e-4)
        cls_lensed2 = data.get_partially_lensed_cls(1, lmax=3000)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_lensed[2:, :], rtol=1e-4)
        cls_lensed2 = data.get_partially_lensed_cls(0, lmax=2500)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_unlensed[2:, :], rtol=1e-4)

        # check lensed CL against python; will only agree well for high lmax as python has no extrapolation template
        cls_lensed2 = correlations.lensed_cls(cls['unlensed_scalar'], cls['lens_potential'][:, 0], delta_cls=False)
        np.testing.assert_allclose(cls_lensed2[2:2000, 2], cls_lensed[2:2000, 2], rtol=1e-3)
        np.testing.assert_allclose(cls_lensed2[2:2000, 1], cls_lensed[2:2000, 1], rtol=1e-3)
        np.testing.assert_allclose(cls_lensed2[2:2000, 0], cls_lensed[2:2000, 0], rtol=1e-3)
        self.assertTrue(np.all(np.abs((cls_lensed2[2:3000, 3] - cls_lensed[2:3000, 3]) /
                                      np.sqrt(cls_lensed2[2:3000, 0] * cls_lensed2[2:3000, 1])) < 1e-4))

        corr, xvals, weights = correlations.gauss_legendre_correlation(cls['lensed_scalar'])
        clout = correlations.corr2cl(corr, xvals, weights, 2500)
        self.assertTrue(np.all(np.abs(clout[2:2300, 2] / cls['lensed_scalar'][2:2300, 2] - 1) < 1e-3))

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=78, YHe=0.22)
        pars.set_for_lmax(2000, lens_potential_accuracy=1)
        pars.WantTensors = True
        results = camb.get_transfer_functions(pars)
        from camb import initialpower
        cls = []
        for r in [0, 0.2, 0.4]:
            inflation_params = initialpower.InitialPowerLaw()
            inflation_params.set_params(ns=0.96, r=r, nt=0)
            results.power_spectra_from_transfer(inflation_params, silent=True)
            cls += [results.get_total_cls(CMB_unit='muK')]
        self.assertTrue(np.allclose((cls[1] - cls[0])[2:300, 2] * 2, (cls[2] - cls[0])[2:300, 2], rtol=1e-3))

        # Check generating tensors and scalars together
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67)
        lmax = 2000
        pars.set_for_lmax(lmax, lens_potential_accuracy=1)
        pars.InitPower.set_params(ns=0.96, r=0)
        pars.WantTensors = False
        results = camb.get_results(pars)
        cl1 = results.get_total_cls(lmax, CMB_unit='muK')
        pars.InitPower.set_params(ns=0.96, r=0.1, nt=0)
        pars.WantTensors = True
        results = camb.get_results(pars)
        cl2 = results.get_lensed_scalar_cls(lmax, CMB_unit='muK')
        ctensor2 = results.get_tensor_cls(lmax, CMB_unit='muK')
        results = camb.get_transfer_functions(pars)
        results.Params.InitPower.set_params(ns=1.1, r=1)
        inflation_params = initialpower.InitialPowerLaw()
        inflation_params.set_params(ns=0.96, r=0.05, nt=0)
        results.power_spectra_from_transfer(inflation_params, silent=True)
        cl3 = results.get_lensed_scalar_cls(lmax, CMB_unit='muK')
        ctensor3 = results.get_tensor_cls(lmax, CMB_unit='muK')
        self.assertTrue(np.allclose(ctensor2, ctensor3 * 2, rtol=1e-4))
        self.assertTrue(np.allclose(cl1, cl2, rtol=1e-4))
        # These are identical because all scalar spectra were identical (non-linear corrections change it  otherwise)
        self.assertTrue(np.allclose(cl1, cl3, rtol=1e-4))

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_for_lmax(2500)
        pars.min_l = 2
        res = camb.get_results(pars)
        cls = res.get_lensed_scalar_cls(2000)
        pars.min_l = 1
        res = camb.get_results(pars)
        cls2 = res.get_lensed_scalar_cls(2000)
        np.testing.assert_allclose(cls[2:, 0:2], cls2[2:, 0:2], rtol=1e-4)
        self.assertAlmostEqual(cls2[1, 0], 1.30388e-10, places=13)
        self.assertAlmostEqual(cls[1, 0], 0)

    def testSigmaR(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.InitPower.set_params(ns=0.965, As=2e-9)
        pars.set_matter_power(nonlinear=False)
        results = camb.get_results(pars)
        sigma8 = results.get_sigma8_0()
        self.assertAlmostEqual(sigma8, results.get_sigmaR(8)[-1], places=3)
        self.assertAlmostEqual(sigma8, results.get_sigmaR(np.array([8]), z_indices=-1)[-1], places=3)
        self.assertAlmostEqual(results.get_sigmaR(8)[-1], results.get_sigmaR(8, z_indices=-1))
        pars.set_matter_power(nonlinear=False, k_per_logint=0, kmax=2)

        results = camb.get_results(pars)
        P, z, k = results.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False,
                                                        return_z_k=True, extrap_kmax=100, silent=True)
        truth = 0.800679  # from high kmax, high accuracy boost
        self.assertTrue(abs(results.get_sigmaR(8)[-1] / sigma8 - 1) < 1e-3)

        def get_sigma(_ks, dlogk):
            x = _ks * 8 / (pars.H0 / 100)
            w = (3 * (np.sin(x) - x * np.cos(x)) / x ** 3) ** 2
            w[x < 1e-2] = 1 - x[x < 1e-2] ** 2 / 2
            Ps = P.P(0, _ks) * _ks ** 3 / (2 * np.pi ** 2)
            return np.sqrt(np.dot(w, Ps * dlogk))

        logk = np.arange(np.log(1e-5), np.log(20.), 1. / 100)
        ks = np.exp(logk)
        py_sigma = get_sigma(ks, logk[1] - logk[0])
        self.assertAlmostEqual(py_sigma, truth, places=3)
        # no interpolation
        logk = np.log(k)
        diffs = (logk[2:] - logk[:-2]) / 2
        ks = k[1:-1]
        py_sigma2 = get_sigma(ks, diffs)
        self.assertAlmostEqual(py_sigma2, truth, places=3)
        self.assertTrue(abs(results.get_sigmaR(8)[-1] / truth - 1) < 1e-4)
        self.assertTrue(abs(results.get_sigmaR(np.array([8]), z_indices=-1)[-1] / truth - 1) < 1e-4)
        pars.set_matter_power(nonlinear=False, k_per_logint=0, kmax=1.2, redshifts=np.arange(0, 10, 2))
        results = camb.get_results(pars)
        sigmas = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        pars.Accuracy.AccuracyBoost = 2
        results = camb.get_results(pars)
        sigmas2 = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        self.assertTrue(np.all(np.abs(sigmas / sigmas2 - 1) < 1e-3))
        pars.Accuracy.AccuracyBoost = 1
        pars.set_matter_power(nonlinear=False, k_per_logint=100, kmax=10, redshifts=np.arange(0, 10, 2))
        results = camb.get_results(pars)
        sigmas2 = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        self.assertAlmostEqual(sigmas2[4, 2], 1.77346, places=3)
        self.assertTrue(np.all(np.abs(sigmas[:, 1:] / sigmas2[:, 1:] - 1) < 1e-3))
        self.assertTrue(np.all(np.abs(sigmas[:, 0] / sigmas2[:, 0] - 1) < 2e-3))

    def testTimeTransfers(self):
        from camb import initialpower
        pars = camb.set_params(H0=69, YHe=0.22, lmax=2000, lens_potential_accuracy=1, ns=0.96, As=2.5e-9)
        results1 = camb.get_results(pars)
        cl1 = results1.get_total_cls()

        pars = camb.set_params(H0=69, YHe=0.22, lmax=2000, lens_potential_accuracy=1)
        results = camb.get_transfer_functions(pars, only_time_sources=True)
        inflation_params = initialpower.InitialPowerLaw()
        inflation_params.set_params(ns=0.96, As=2.5e-9)
        results.power_spectra_from_transfer(inflation_params)
        cl2 = results.get_total_cls()
        np.testing.assert_allclose(cl1, cl2, rtol=1e-4)
        inflation_params.set_params(ns=0.96, As=1.9e-9)
        results.power_spectra_from_transfer(inflation_params)
        inflation_params.set_params(ns=0.96, As=2.5e-9)
        results.power_spectra_from_transfer(inflation_params)
        cl2 = results.get_total_cls()
        np.testing.assert_allclose(cl1, cl2, rtol=1e-4)

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=78, YHe=0.22)
        pars.set_for_lmax(2000, lens_potential_accuracy=1)
        pars.WantTensors = True
        results = camb.get_transfer_functions(pars, only_time_sources=True)
        cls = []
        for r in [0, 0.2, 0.4]:
            inflation_params = initialpower.InitialPowerLaw()
            inflation_params.set_params(ns=0.96, r=r, nt=0)
            results.power_spectra_from_transfer(inflation_params)
            cls += [results.get_total_cls(CMB_unit='muK')]
        self.assertTrue(np.allclose((cls[1] - cls[0])[2:300, 2] * 2, (cls[2] - cls[0])[2:300, 2], rtol=1e-3))

    def testDarkEnergy(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=71)
        pars.InitPower.set_params(ns=0.965, r=0)
        for m in ['fluid', 'ppf']:
            pars.set_dark_energy(w=-0.7, wa=0.2, dark_energy_model=m)
            C1 = camb.get_results(pars).get_cmb_power_spectra()
            a = np.logspace(-5, 0, 1000)
            w = -0.7 + 0.2 * (1 - a)
            pars2 = pars.copy()
            pars2.set_dark_energy_w_a(a, w, dark_energy_model=m)
            C2 = camb.get_results(pars2).get_cmb_power_spectra()
            for f in ['lens_potential', 'lensed_scalar']:
                self.assertTrue(np.allclose(C1[f][2:, 0], C2[f][2:, 0]))
            pars3 = pars2.copy()
            self.assertAlmostEqual(-0.7, pars3.DarkEnergy.w)

    def testInitialPower(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67)
        import ctypes
        P = camb.InitialPowerLaw()
        P2 = ctypes.pointer(P)
        self.assertEqual(P.As, pars.InitPower.As)
        As = 1.8e-9
        ns = 0.8
        P.set_params(As=As, ns=ns)
        self.assertEqual(P.As, As)
        self.assertEqual(P2.contents.As, As)

        pars2 = camb.CAMBparams()
        pars2.set_cosmology(H0=67)
        pars2.InitPower.set_params(As=1.7e-9, ns=ns)
        self.assertEqual(pars2.InitPower.As, 1.7e-9)
        pars.set_initial_power(pars2.InitPower)
        self.assertEqual(pars.InitPower.As, 1.7e-9)
        pars.set_initial_power(P)
        self.assertEqual(pars.InitPower.As, As)

        ks = np.logspace(-5.5, 2, 1000)
        pk = (ks / P.pivot_scalar) ** (ns - 1) * As
        pars2.set_initial_power_table(ks, pk)
        self.assertAlmostEqual(pars2.scalar_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)
        sp = camb.SplinedInitialPower(ks=ks, PK=pk)
        pars2.set_initial_power(sp)
        self.assertAlmostEqual(pars2.scalar_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)
        self.assertFalse(sp.has_tensors())
        self.assertFalse(pars2.InitPower.has_tensors())

        sp = camb.SplinedInitialPower()
        sp.set_scalar_log_regular(10 ** (-5.5), 10. ** 2, pk)
        pars2.set_initial_power(sp)
        self.assertAlmostEqual(pars2.scalar_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)

        sp.set_tensor_log_regular(10 ** (-5.5), 10. ** 2, pk)
        pars2.set_initial_power(sp)
        self.assertAlmostEqual(pars2.tensor_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)
        self.assertTrue(sp.has_tensors())
        sp.set_tensor_table([], [])
        self.assertFalse(sp.has_tensors())
        pars2.set_initial_power(sp)

        results = camb.get_results(pars2)
        cl = results.get_lensed_scalar_cls(CMB_unit='muK')
        pars.InitPower.set_params(As=As, ns=ns)
        results2 = camb.get_results(pars)
        cl2 = results2.get_lensed_scalar_cls(CMB_unit='muK')
        self.assertTrue(np.allclose(cl, cl2, rtol=1e-4))
        P = camb.InitialPowerLaw(As=2.1e-9, ns=0.9)
        pars2.set_initial_power(P)
        pars.InitPower.set_params(As=2.1e-9, ns=0.9)
        self.assertAlmostEqual(pars2.scalar_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)

        def PK(k, A, n):
            return A * (k / 0.05) ** (n - 1) * (1 + 0.1 * np.sin(10 * k))

        pars.set_initial_power_function(PK, args=(3e-9, 0.95))
        P = pars.scalar_power(ks)
        np.testing.assert_almost_equal(P, PK(ks, 3e-9, 0.95), decimal=4)

    # noinspection PyTypeChecker
    def testSources(self):
        from camb.sources import GaussianSourceWindow, SplinedSourceWindow
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=64, mnu=0)
        pars.set_for_lmax(1200)
        pars.Want_CMB = False
        pars.SourceWindows = [
            GaussianSourceWindow(redshift=0.17, source_type='counts', bias=1.2, sigma=0.04, dlog10Ndm=-0.2),
            GaussianSourceWindow(redshift=0.5, source_type='lensing', sigma=0.07, dlog10Ndm=0)]
        pars.SourceTerms.limber_windows = True
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        zs = np.arange(0, 0.5, 0.02)
        W = np.exp(-(zs - 0.17) ** 2 / 2 / 0.04 ** 2) / np.sqrt(2 * np.pi) / 0.04

        ks = np.logspace(-4, 3, 50)
        bias_kz = 1.2 * np.ones((len(ks), len(zs)))
        test_windows = [SplinedSourceWindow(bias=1.2, dlog10Ndm=-0.2, z=zs, W=W),
                        SplinedSourceWindow(bias_z=1.2 * np.ones_like(zs), dlog10Ndm=-0.2, z=zs, W=W),
                        SplinedSourceWindow(k_bias=ks, bias_kz=bias_kz, dlog10Ndm=-0.2, z=zs, W=W)]
        for window in test_windows:
            pars.SourceWindows[0] = window
            results = camb.get_results(pars)
            cls2 = results.get_source_cls_dict()
            self.assertTrue(np.allclose(cls2["W1xW1"][2:1200], cls["W1xW1"][2:1200], rtol=1e-3))

        pars.SourceWindows = [GaussianSourceWindow(redshift=1089, source_type='lensing', sigma=30)]
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        PP = cls['PxP']
        ls = np.arange(0, PP.shape[0])
        self.assertTrue(np.allclose(PP / 4 * (ls * (ls + 1)), cls['W1xW1'], rtol=1e-3))
        self.assertTrue(np.allclose(PP / 2 * np.sqrt(ls * (ls + 1)), cls['PxW1'], rtol=1e-3))
        # test something sharp with redshift distortions (tricky..)
        from scipy import signal
        zs = np.arange(1.9689, 2.1057, (2.1057 - 1.9689) / 2000)
        W = signal.windows.tukey(len(zs), alpha=0.1)
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
        pars.InitPower.set_params(As=2e-9, ns=0.965)
        pars.set_for_lmax(4000)
        pars.SourceWindows = [SplinedSourceWindow(z=zs, W=W, source_type='counts')]
        pars.SourceTerms.counts_redshift = True
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        self.assertAlmostEqual(np.sum(cls['PxW1'][10:3000:20]), 0.00020001, places=5)
        self.assertAlmostEqual(np.sum(cls['W1xW1'][10:3000:20]), 2.26413, places=3)
        self.assertAlmostEqual(np.sum(cls['W1xW1'][10]), 0.0001097, places=6)

    def testSymbolic(self):
        if fast:
            return
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
            pars.set_custom_scalar_sources([monopole_source + ISW + doppler + quadrupole_source,
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

        s.internal_consistency_checks()

    def test_mathutils(self):
        from camb.mathutils import chi_squared, threej_coupling, scalar_coupling_matrix, pcl_coupling_matrix
        cinv = np.linalg.inv(np.array([[1.2, 3], [3, 18.2]]))
        vec = np.array([0.5, 5.])
        self.assertAlmostEqual(chi_squared(cinv, vec), cinv.dot(vec).dot(vec))
        W = np.zeros(100)
        W[0] = 1
        lmax = len(W)
        Xi = threej_coupling(W, lmax)
        np.testing.assert_allclose(np.diag(Xi) * (2 * np.arange(lmax + 1) + 1), np.ones(lmax + 1))
        Xis = threej_coupling(W, lmax, pol=True)
        np.testing.assert_allclose(np.diag(Xis[0]) * (2 * np.arange(lmax + 1) + 1), np.ones(lmax + 1))
        P = W * 4 * np.pi
        M = scalar_coupling_matrix(P, lmax)
        np.testing.assert_allclose(M, np.eye(lmax + 1))
        M = pcl_coupling_matrix(P, lmax)
        np.testing.assert_allclose(M, np.eye(lmax + 1))

    def test_extra_EmissionAnglePostBorn(self):
        if fast:
            return
        from camb import emission_angle, postborn
        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, tau=0.055)
        BB = emission_angle.get_emission_delay_BB(pars, lmax=3500)
        self.assertAlmostEqual(BB(80) * 2 * np.pi / 80 / 81., 1.1e-10, delta=1e-11)

        Bom = postborn.get_field_rotation_BB(pars, lmax=3500)
        self.assertAlmostEqual(Bom(100) * 2 * np.pi / 100 / 101., 1.65e-11, delta=1e-12)

    def test_memory(self):
        if platform.system() != "Windows":
            import gc
            import resource
            last_usage = -1
            for i in range(3):
                pars = camb.CAMBparams()
                pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.12, mnu=0.06, omk=0, tau=0.17)
                results = camb.get_results(pars)  # noqa
                del pars, results
                gc.collect()
                usage = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1)
                if 0 < last_usage != usage:
                    print('Memory usage: %2.2f KB vs %2.2f KB' % (usage, last_usage))
                    raise Exception("Apparent memory leak")
                last_usage = usage

    def test_quintessence(self):
        n = 3
        # set zc and fde_zc
        pars = camb.set_params(ombh2=0.022, omch2=0.122, thetastar=0.01044341764253,
                               dark_energy_model='EarlyQuintessence',
                               m=8e-53, f=0.05, n=n, theta_i=3.1, use_zc=True, zc=1e4, fde_zc=0.1)
        camb.get_background(pars)
        results = camb.get_results(pars)
        self.assertAlmostEqual(results.get_derived_params()['thetastar'], 1.044341764253, delta=1e-5)
