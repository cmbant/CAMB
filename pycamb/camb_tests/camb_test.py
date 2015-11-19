from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import unittest

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    import camb
from camb import model


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

        r = data.comoving_radial_distance(0.48)
        t0 = data.conformal_time(0)
        t1 = data.conformal_time(11.5)
        t2 = data.comoving_radial_distance(11.5)
        self.assertAlmostEqual(t2, t0 - t1, 2)
        self.assertAlmostEqual(t1, 4200.78, 2)

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

        # test massive sterile models as in Planck papers
        pars.set_cosmology(H0=68.0, ombh2=0.022305, omch2=0.11873, mnu=0.06, nnu=3.073, omk=0, meffsterile=0.013)
        self.assertAlmostEquals(pars.omegan * (pars.H0 / 100) ** 2, 0.00078, 5)
        self.assertAlmostEquals(pars.YHe, 0.24573, 5)

        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()['age'], 13.773, 2)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104103, 6)

        # test dark energy
        pars.set_cosmology(H0=68.26, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars.set_dark_energy(w=-1.0226, dark_energy_model='fluid')
        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()['age'], 13.789, 2)
        scal= data.luminosity_distance(1.4)
        vec = data.luminosity_distance([1.2, 1.4, 0.1, 1.9])
        self.assertAlmostEqual(scal,vec[1],5)
        pars.set_dark_energy() #re-set defaults


    def testPowers(self):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_dark_energy() #re-set defaults
        pars.InitPower.set_params(ns=0.965)
        pars.set_matter_power(redshifts=[0., 0.17, 3.1])
        pars.NonLinear = model.NonLinear_none
        print(pars)
        data = camb.get_results(pars)
        print(data.Params)

        kh, z, pk = data.get_matter_power_spectrum(1e-4, 1, 20)

        kh2, z2, pk2 = data.get_linear_matter_power_spectrum()

        s8 = data.get_sigma8()
        self.assertAlmostEqual(s8[0], 0.247, 3)

        pars.NonLinear = model.NonLinear_both
        data.calc_power_spectra(pars)
        kh3, z3, pk3 = data.get_matter_power_spectrum(1e-4, 1, 20)
        self.assertAlmostEquals(pk[-1][-3], 51.909, 2)
        self.assertAlmostEquals(pk3[-1][-3], 57.697, 2)
        self.assertAlmostEquals(pk2[-2][-4], 53.47, 2)

        camb.set_feedback_level(0)
        cls = data.get_cmb_power_spectra(pars)

        MTrans = data.get_matter_transfer_data()

        cls_tot = data.get_total_cls(2000)
        cls_scal = data.get_unlensed_scalar_cls(2000)
        cls_tensor = data.get_tensor_cls(2000)
        cls_lensed = data.get_lensed_scalar_cls(2000)
        cls_phi = data.get_lens_potential_cls(2000)
