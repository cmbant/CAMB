import os
import pickle
import platform
import subprocess
import sys
import tempfile
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb
from camb import bbn, constants, correlations, dark_energy, initialpower, model, recombination
from camb.baseconfig import CAMBError, CAMBParamRangeError, CAMBValueError


def new_def_params(**kwargs):
    pars = camb.CAMBparams(**kwargs)
    pars.Recomb.set_params(recfast_approx_model=recombination.recfast_planck)
    return pars


def def_set_params(**kwargs):
    kwargs.setdefault("recfast_approx_model", recombination.recfast_planck)
    return camb.set_params(**kwargs)


class CambTest(unittest.TestCase):
    def testAssigments(self):
        ini = os.path.join(os.path.dirname(__file__), "..", "inifiles", "planck_2018.ini")
        if os.path.exists(ini):
            pars = camb.read_ini(ini)
            self.assertTrue(np.abs(camb.get_background(pars).cosmomc_theta() * 100 / 1.040909 - 1) < 2e-5)
        pars = new_def_params()
        pars.set_cosmology(H0=68.5, ombh2=0.022, mnu=0, omch2=0.1)
        self.assertAlmostEqual(pars.omegam, (0.022 + 0.1) / 0.685**2)
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
        self.assertTrue("Want_CMB_lensing = True" in printstr and "NonLinear = NonLinear_both" in printstr)
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
        self.assertFalse(len(pars.nu_mass_degeneracies[:1]))
        pars = def_set_params(**{"InitPower.ns": 1.2, "WantTransfer": True})
        self.assertEqual(pars.InitPower.ns, 1.2)
        self.assertTrue(pars.WantTransfer)
        pars.DarkEnergy = None
        pars = def_set_params(**{"H0": 67, "ombh2": 0.002, "r": 0.1, "Accuracy.AccurateBB": True})
        self.assertEqual(pars.Accuracy.AccurateBB, True)

        pars = new_def_params()
        pars.set_for_lmax(2500)
        self.assertEqual(pars.max_eta_k, 72000)
        self.assertEqual(pars.NonLinear, model.NonLinear_lens)
        pars = new_def_params()
        pars.set_for_lmax(2500, lens_potential_accuracy=0)
        self.assertEqual(pars.max_eta_k, (2500 + 200) * 2.5)
        self.assertEqual(pars.NonLinear, model.NonLinear_none)
        pars = def_set_params(lmax=4000)
        self.assertEqual(pars.max_eta_k, 90000)
        pars = def_set_params(lmax=4000, lens_potential_accuracy=0)
        self.assertEqual(pars.max_eta_k, (4000 + 200) * 2.5)
        cosmomc_params = {
            "H0": 67,
            "omegabh2": 0.022,
            "omegach2": 0.12,
            "tau": 0.054,
            "ns": 0.965,
            "A": 2.1,
        }
        pars = camb.set_params_cosmomc(cosmomc_params)
        self.assertEqual(pars.max_eta_k, 18000)
        pars = camb.set_params_cosmomc(cosmomc_params, lens_potential_accuracy=None)
        self.assertEqual(pars.max_eta_k, 72000)

        from camb import check_accuracy

        class FakeAccuracyParams:
            DoLensing = True
            max_l = 1000
            max_eta_k = 2500.0
            NonLinear = model.NonLinear_none

            def __init__(self):
                self.set_for_lmax_calls = []

            def set_for_lmax(self, lmax, **kwargs):
                self.set_for_lmax_calls.append((lmax, kwargs))
                self.max_l = lmax + kwargs.get("lens_output_margin", 200)
                self.max_eta_k = kwargs.get("max_eta_k") or self.max_l * kwargs.get("k_eta_fac", 2.5)

        fake = FakeAccuracyParams()
        check_accuracy.apply_lensing_settings(fake, set_for_lmax=4000)
        self.assertEqual(fake.set_for_lmax_calls[0][1]["lens_potential_accuracy"], 0.0)
        self.assertFalse(fake.set_for_lmax_calls[0][1]["nonlinear"])
        fake = FakeAccuracyParams()
        check_accuracy.apply_lensing_settings(fake, lens_output_margin=200)
        self.assertEqual(fake.set_for_lmax_calls[0][1]["lens_potential_accuracy"], 0.0)
        self.assertEqual(fake.set_for_lmax_calls[0][1]["max_eta_k"], 2500.0)

        from camb.sources import GaussianSourceWindow

        pars = new_def_params()
        pars.SourceWindows = [GaussianSourceWindow(), GaussianSourceWindow(redshift=1)]
        self.assertEqual(pars.SourceWindows[1].redshift, 1)
        pars.SourceWindows[0].redshift = 2
        self.assertEqual(pars.SourceWindows[0].redshift, 2)
        self.assertTrue(len(pars.SourceWindows) == 2)
        pars.SourceWindows[0] = GaussianSourceWindow(redshift=3)
        self.assertEqual(pars.SourceWindows[0].redshift, 3)
        self.assertTrue("redshift = 3.0" in str(pars))
        pars.SourceWindows = pars.SourceWindows[0:1]
        self.assertTrue(len(pars.SourceWindows) == 1)
        pars.SourceWindows = []
        self.assertTrue(len(pars.SourceWindows) == 0)
        params = camb.get_valid_numerical_params()
        self.assertEqual(
            params,
            {
                "ombh2",
                "deltazrei",
                "omnuh2",
                "tau",
                "omk",
                "zrei",
                "thetastar",
                "nrunrun",
                "meffsterile",
                "nnu",
                "ntrun",
                "HMCode_A_baryon",
                "HMCode_eta_baryon",
                "HMCode_logT_AGN",
                "cosmomc_theta",
                "YHe",
                "wa",
                "cs2",
                "H0",
                "mnu",
                "Alens",
                "TCMB",
                "ns",
                "nrun",
                "As",
                "nt",
                "r",
                "w",
                "omch2",
                "max_zrei",
                "wde_a_array",
                "wde_w_array",
            },
        )
        params2 = camb.get_valid_numerical_params(dark_energy_model="AxionEffectiveFluid")
        self.assertEqual(params2.difference(params), {"fde_zc", "w_n", "zc", "theta_i"})
        pars = def_set_params(
            H0=67,
            ombh2=0.022,
            omch2=0.12,
            dark_energy_model="AxionEffectiveFluid",
            w_n=0.4,
            fde_zc=0.05,
            zc=4000,
        )
        self.assertIsInstance(pars.DarkEnergy, dark_energy.AxionEffectiveFluid)
        self.assertEqual(pars.DarkEnergy.w_n, 0.4)
        self.assertEqual(pars.DarkEnergy.fde_zc, 0.05)
        self.assertEqual(pars.DarkEnergy.zc, 4000)

    def testWriteIniRoundTrip(self):
        script = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "fortran", "tests", "CAMB_test_files.py")
        )
        script_dir = os.path.dirname(script)
        repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
        fortran_dir = os.path.join(repo_root, "fortran")
        base_settings = os.path.join(repo_root, "inifiles", "params.ini")

        with tempfile.TemporaryDirectory() as ini_dir:
            subprocess.run(
                [
                    sys.executable,
                    script,
                    ini_dir,
                    "--make_ini",
                    "--no_run_test",
                    "--base_settings",
                    base_settings,
                ],
                check=True,
                cwd=fortran_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            ini_files = sorted(
                os.path.join(ini_dir, filename)
                for filename in os.listdir(ini_dir)
                if filename.endswith(".ini") and filename.startswith("params_")
            )
            self.assertTrue(ini_files)

            cwd = os.getcwd()
            os.chdir(fortran_dir)
            # Several test ini files (accurate_BB, share_delta_neff with nu_mass_degeneracies, ...)
            # intentionally trigger Fortran warnings on read; silence them to keep the test output clean.
            prior_print_fortran_warnings = camb.config.print_fortran_warnings
            camb.config.print_fortran_warnings = False
            try:
                for ini_file in ini_files:
                    with self.subTest(ini=os.path.basename(ini_file)):
                        pars = camb.read_ini(ini_file)
                        written_ini = os.path.join(ini_dir, "written_" + os.path.basename(ini_file))
                        camb.write_ini(pars, written_ini)
                        self.assertTrue(os.path.exists(written_ini))
            finally:
                camb.config.print_fortran_warnings = prior_print_fortran_warnings
                os.chdir(cwd)

    def testWriteIniFromPythonParams(self):
        pars = new_def_params()
        pars.set_cosmology(H0=67, ombh2=0.0224, omch2=0.119, tau=0.054, mnu=0.06)
        pars.set_dark_energy(w=-0.95, wa=0.15, dark_energy_model="ppf")
        pars.InitPower.set_params(As=2.1e-9, ns=0.965, nrun=0.01, r=0.03, nt=0.0)
        pars.set_matter_power(
            redshifts=[0.0, 0.5, 1.0], kmax=2.0, accurate_massive_neutrino_transfers=True, silent=True
        )
        pars.set_for_lmax(2200, lens_potential_accuracy=1)
        pars.WantTensors = True
        pars.NonLinearModel.set_params(halofit_version="mead2020_feedback", HMCode_logT_AGN=7.7)
        pars.Alens = 0.95

        with tempfile.TemporaryDirectory() as temp_dir:
            ini_file = os.path.join(temp_dir, "python_params.ini")
            pars.write_ini(ini_file)
            self.assertTrue(os.path.exists(ini_file))

    def testIniThetaInput(self):
        base_ini = os.path.join(os.path.dirname(__file__), "..", "..", "inifiles", "planck_2018.ini")
        base = camb.read_ini(base_ini)
        theta = camb.get_background(base).get_derived_params()["thetastar"] / 100

        py_pars = camb.read_ini(base_ini)
        py_pars.set_H0_for_theta(theta)

        with open(base_ini, encoding="utf-8") as handle:
            text = handle.read()

        theta_text = text.replace("hubble         = 67.32117", f"thetastar= {theta:.12f}")

        with tempfile.TemporaryDirectory() as temp_dir:
            theta_ini = os.path.join(temp_dir, "theta.ini")
            with open(theta_ini, "w", encoding="utf-8") as handle:
                handle.write(theta_text)

            ini_pars = camb.read_ini(theta_ini)
            ini_theta = camb.get_background(ini_pars).get_derived_params()["thetastar"] / 100
            py_theta = camb.get_background(py_pars).get_derived_params()["thetastar"] / 100
            self.assertAlmostEqual(ini_pars.H0, base.H0, places=4)
            self.assertAlmostEqual(ini_pars.H0, py_pars.H0, delta=1e-3)
            self.assertAlmostEqual(ini_theta, theta, places=7)
            self.assertAlmostEqual(ini_theta, py_theta, places=7)

            both_ini = os.path.join(temp_dir, "theta_and_hubble.ini")
            with open(both_ini, "w", encoding="utf-8") as handle:
                handle.write(theta_text + f"\nhubble         = {base.H0:.5f}\n")

            with self.assertRaises(CAMBValueError):
                camb.read_ini(both_ini)

    def testBackground(self):
        pars = new_def_params()
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
        self.assertAlmostEqual(derived["age"], age, 2)
        self.assertAlmostEqual(derived["rdrag"], 146.986, 2)
        self.assertAlmostEqual(derived["rstar"], data.sound_horizon(derived["zstar"]), 2)

        # Test BBN consistency, base_plikHM_TT_lowTEB best fit model
        pars.set_cosmology(H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.YHe, 0.24560055, 5)
        data.calc_background(pars)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104090733, 7)
        self.assertAlmostEqual(data.get_derived_params()["kd"], 0.14055499, 4)

        pars.set_cosmology(
            H0=67.31, ombh2=0.022242, omch2=0.11977, mnu=0.06, omk=0, bbn_predictor=bbn.BBN_table_interpolator()
        )
        self.assertAlmostEqual(pars.YHe, 0.24560055, 5)
        self.assertAlmostEqual(pars.get_Y_p(), bbn.BBN_table_interpolator().Y_p(0.022242, 0), 5)

        # test massive sterile models as in Planck papers
        pars.set_cosmology(H0=68.0, ombh2=0.022305, omch2=0.11873, mnu=0.06, nnu=3.073, omk=0, meffsterile=0.013)
        self.assertAlmostEqual(pars.omnuh2, 0.00078, 5)
        self.assertAlmostEqual(pars.YHe, 0.24601443, 5)
        self.assertAlmostEqual(pars.N_eff, 3.073, 4)

        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()["age"], 13.773, 2)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104102759, 6)

        # test dark energy
        pars.set_cosmology(H0=68.26, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars.set_dark_energy(w=-1.0226, dark_energy_model="fluid")

        data.calc_background(pars)
        self.assertAlmostEqual(data.get_derived_params()["age"], 13.789, 2)
        scal = data.luminosity_distance(1.4)
        vec = data.luminosity_distance([1.2, 1.4, 0.1, 1.9])
        self.assertAlmostEqual(scal, vec[1], 5)

        pars.set_dark_energy()  # re-set defaults

        # test theta
        pars.set_cosmology(cosmomc_theta=0.0104085, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        self.assertAlmostEqual(pars.H0, 67.537, 2)
        with self.assertRaises(CAMBParamRangeError):
            pars.set_cosmology(cosmomc_theta=0.0204085, ombh2=0.022271, omch2=0.11914, mnu=0.06, omk=0)
        pars = def_set_params(cosmomc_theta=0.0104077, ombh2=0.022, omch2=0.122, w=-0.95)
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        pars = def_set_params(thetastar=0.010311, ombh2=0.022, omch2=0.122)
        self.assertAlmostEqual(camb.get_background(pars).get_derived_params()["thetastar"] / 100, 0.010311, 7)
        pars = def_set_params(thetastar=0.010311, ombh2=0.022, omch2=0.122, omk=-0.05)
        self.assertAlmostEqual(camb.get_background(pars).get_derived_params()["thetastar"] / 100, 0.010311, 7)
        self.assertAlmostEqual(pars.H0, 49.70523, places=3)

        pars = def_set_params(cosmomc_theta=0.0104077, ombh2=0.022, omch2=0.122, w=-0.95, wa=0, dark_energy_model="ppf")
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        pars = def_set_params(
            cosmomc_theta=0.0104077,
            ombh2=0.022,
            omch2=0.122,
            w=-0.95,
            dark_energy_model="DarkEnergyFluid",
            initial_power_model="InitialPowerLaw",
        )
        self.assertAlmostEqual(camb.get_background(pars, no_thermo=True).cosmomc_theta(), 0.0104077, 7)

        with self.assertRaises(CAMBValueError):
            def_set_params(dark_energy_model="InitialPowerLaw")
        data.calc_background(pars)
        h2 = (data.Params.H0 / 100) ** 2
        self.assertAlmostEqual(data.get_Omega("baryon"), data.Params.ombh2 / h2, 7)
        self.assertAlmostEqual(data.get_Omega("nu"), data.Params.omnuh2 / h2, 7)
        self.assertAlmostEqual(
            data.get_Omega("photon")
            + data.get_Omega("neutrino")
            + data.get_Omega("de")
            + (pars.ombh2 + pars.omch2 + pars.omnuh2) / h2
            + pars.omk,
            1,
            8,
        )
        pars.set_cosmology(H0=67, mnu=1, neutrino_hierarchy="normal")
        data.calc_background(pars)
        h2 = (pars.H0 / 100) ** 2
        self.assertAlmostEqual(
            data.get_Omega("photon")
            + data.get_Omega("neutrino")
            + data.get_Omega("de")
            + (pars.ombh2 + pars.omch2 + pars.omnuh2) / h2
            + pars.omk,
            1,
            8,
        )
        redshifts = np.array([0.005, 0.01, 0.3, 0.9342, 4, 27, 321.5, 932, 1049, 1092, 2580, 1e4, 2.1e7])
        np.testing.assert_allclose(
            data.redshift_at_conformal_time(data.conformal_time(redshifts)), redshifts, rtol=1e-7
        )
        pars.set_dark_energy(w=-1.8)
        data.calc_background(pars)
        np.testing.assert_allclose(
            data.redshift_at_conformal_time(data.conformal_time(redshifts)), redshifts, rtol=1e-7
        )
        pars.set_cosmology(cosmomc_theta=0.0104085)
        data.calc_background(pars)
        self.assertAlmostEqual(data.cosmomc_theta(), 0.0104085)
        derived = data.get_derived_params()
        pars.Accuracy.BackgroundTimeStepBoost = 2
        data.calc_background(pars)
        derived2 = data.get_derived_params()
        self.assertAlmostEqual(derived["thetastar"], derived2["thetastar"], places=5)
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.11, neutrino_hierarchy="inverted")
        self.assertEqual(pars.num_nu_massive, 3)
        self.assertEqual(pars.nu_mass_numbers[1], 1)
        self.assertEqual(pars.nu_mass_eigenstates, 2)
        self.assertAlmostEqual(pars.nu_mass_fractions[0], 0.915197, places=4)

        pars = new_def_params()
        pars.set_cosmology(H0=68.5, ombh2=0.022, omch2=0.122, YHe=0.2453, mnu=0.07, omk=0, zrei=zre)
        results = camb.get_background(pars)
        self.assertEqual(results.Params.Reion.redshift, zre)

        pars = new_def_params()
        pars.set_cosmology(H0=68.5, ombh2=0.022, omch2=0.122, YHe=0.2453, mnu=0.07, omk=-0.05)
        data = camb.get_background(pars)
        delta2 = (
            data.curvature_radius
            / (1 + 0.25)
            * (
                np.sin(
                    (data.comoving_radial_distance(0.25) - data.comoving_radial_distance(0.05)) / data.curvature_radius
                )
            )
        )
        np.testing.assert_allclose(delta2, data.angular_diameter_distance2(0.05, 0.25))
        dists = data.angular_diameter_distance2([0.3, 0.05, 0.25], [1, 0.25, 0.05])
        self.assertAlmostEqual(delta2, dists[1])
        self.assertEqual(0, dists[2])

        self.assertEqual(data.physical_time(0.4), data.physical_time([0.2, 0.4])[1])
        d = data.conformal_time_a1_a2(0, 0.5) + data.conformal_time_a1_a2(0.5, 1)
        self.assertAlmostEqual(d, data.conformal_time_a1_a2(0, 1))
        self.assertAlmostEqual(d, sum(data.conformal_time_a1_a2([0, 0.5], [0.5, 1])))

    def testRecfastRosenbrockAgreement(self):
        redshifts = np.geomspace(1.0, 3001.0, 400) - 1.0

        def make_pars(nz, use_rosenbrock=False, handoff=0.985):
            pars = new_def_params()
            pars.set_cosmology(H0=67.4, ombh2=0.02237, omch2=0.12, mnu=0.06, tau=0.054, YHe=0.2453)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.Recomb.Nz = nz
            pars.Recomb.use_rosenbrock = use_rosenbrock
            pars.Recomb.rosenbrock_handoff_xH = handoff
            pars.Recomb.rosenbrock_tol = 3e-4
            return pars

        for nz, handoff in [(2048, 0.976), (10000, 0.985)]:
            with self.subTest(nz=nz, handoff=handoff):
                base = camb.get_background(make_pars(nz))
                ros = camb.get_background(make_pars(nz, use_rosenbrock=True, handoff=handoff))

                base_hist = base.get_background_redshift_evolution(redshifts, vars=["x_e", "T_b"], format="array")
                ros_hist = ros.get_background_redshift_evolution(redshifts, vars=["x_e", "T_b"], format="array")

                xe_denom = np.maximum(np.maximum(np.abs(base_hist[:, 0]), np.abs(ros_hist[:, 0])), 1e-12)
                tb_denom = np.maximum(np.maximum(np.abs(base_hist[:, 1]), np.abs(ros_hist[:, 1])), 1e-12)
                xe_rel = (base_hist[:, 0] - ros_hist[:, 0]) / xe_denom
                tb_rel = (base_hist[:, 1] - ros_hist[:, 1]) / tb_denom

                base_derived = base.get_derived_params()
                ros_derived = ros.get_derived_params()
                theta_rel = abs(ros_derived["thetastar"] / base_derived["thetastar"] - 1.0)
                zstar_rel = abs(ros_derived["zstar"] / base_derived["zstar"] - 1.0)

                np.testing.assert_allclose(xe_rel, 0, atol=1e-3)
                np.testing.assert_allclose(tb_rel, 0, atol=1e-6)
                self.assertLess(theta_rel, 1e-6)
                self.assertLess(zstar_rel, 1e-6)

        nz = 2048
        delta_z = 1.0e4 / nz
        nodes = 1.0e4 - np.arange(1, nz + 1, dtype=np.float64) * delta_z
        base_nodes = camb.get_background(make_pars(nz)).get_background_redshift_evolution(
            nodes, vars=["x_e"], format="array"
        )[:, 0]
        snap_index = np.flatnonzero(base_nodes < 0.976)[0]
        upper_x = base_nodes[snap_index - 1]
        lower_x = base_nodes[snap_index]
        handoff_a = lower_x + 0.25 * (upper_x - lower_x)
        handoff_b = lower_x + 0.75 * (upper_x - lower_x)
        ros_a = camb.get_background(make_pars(nz, use_rosenbrock=True, handoff=handoff_a))
        ros_b = camb.get_background(make_pars(nz, use_rosenbrock=True, handoff=handoff_b))
        hist_a = ros_a.get_background_redshift_evolution(redshifts, vars=["x_e", "T_b"], format="array")
        hist_b = ros_b.get_background_redshift_evolution(redshifts, vars=["x_e", "T_b"], format="array")
        np.testing.assert_allclose(hist_a, hist_b, rtol=0, atol=1e-15)

        with tempfile.TemporaryDirectory() as temp_dir:
            ini_path = os.path.join(temp_dir, "recfast_rosenbrock.ini")
            make_pars(2048, use_rosenbrock=True, handoff=0.985).write_ini(ini_path)
            roundtrip = camb.read_ini(ini_path)
            self.assertTrue(roundtrip.Recomb.use_rosenbrock)
            self.assertEqual(roundtrip.Recomb.Nz, 2048)
            self.assertAlmostEqual(roundtrip.Recomb.rosenbrock_handoff_xH, 0.985)
            self.assertAlmostEqual(roundtrip.Recomb.rosenbrock_tol, 3e-4)
            self.assertAlmostEqual(roundtrip.Recomb.RECFAST_fudge, 1.125)

            with open(ini_path) as ini_file:
                ini_lines = ini_file.readlines()
            self.assertTrue(any(line.startswith("RECFAST_H_fudge") for line in ini_lines))
            self.assertFalse(any(line.startswith("RECFAST_fudge =") for line in ini_lines))

            legacy_ini_path = os.path.join(temp_dir, "legacy_recfast_fudge.ini")
            with open(legacy_ini_path, "w") as legacy_ini:
                for line in ini_lines:
                    if line.startswith("RECFAST_H_fudge"):
                        legacy_ini.write("RECFAST_fudge = 1.14\n")
                    else:
                        legacy_ini.write(line)
            legacy_roundtrip = camb.read_ini(legacy_ini_path)
            self.assertAlmostEqual(legacy_roundtrip.Recomb.RECFAST_fudge, 1.125)

    def testTCMBRecombinationConsistency(self):
        default_tcmb = constants.COBE_CMBTemp
        low_tcmb = 1.7
        density_scale = (low_tcmb / default_tcmb) ** 3
        default_redshifts = np.geomspace(801.0, 3001.0, 80) - 1.0
        low_tcmb_redshifts = default_tcmb / low_tcmb * (1.0 + default_redshifts) - 1.0

        def make_pars(TCMB, scale, recfast_approx_model):
            pars = new_def_params()
            pars.Recomb.set_params(recfast_approx_model=recfast_approx_model)
            pars.set_cosmology(
                H0=67.4,
                ombh2=0.02237 * scale,
                omch2=0.12 * scale,
                mnu=0,
                YHe=0.2453,
                TCMB=TCMB,
            )
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            return pars

        for recfast_approx_model in [recombination.recfast_planck, recombination.recfast_cosmorec]:
            with self.subTest(recfast_approx_model=recfast_approx_model):
                default = camb.get_background(make_pars(default_tcmb, 1.0, recfast_approx_model))
                low_tcmb_model = camb.get_background(make_pars(low_tcmb, density_scale, recfast_approx_model))
                default_xe = default.get_background_redshift_evolution(default_redshifts, vars=["x_e"], format="array")[
                    :, 0
                ]
                low_tcmb_xe = low_tcmb_model.get_background_redshift_evolution(
                    low_tcmb_redshifts, vars=["x_e"], format="array"
                )[:, 0]

                np.testing.assert_allclose(low_tcmb_xe, default_xe, rtol=5e-5, atol=0)

    def testRecfastApproxModels(self):
        self.assertEqual(recombination.recfast_default, recombination.recfast_cosmorec)
        default_rec = recombination.Recfast()
        default_expected = recombination.recfast_approx_model_params[recombination.recfast_cosmorec]
        for name, value in default_expected.items():
            self.assertAlmostEqual(getattr(default_rec, name), value)

        rec = recombination.Recfast()
        rec.set_params(recfast_approx_model=recombination.recfast_cosmorec)
        expected = recombination.recfast_approx_model_params[recombination.recfast_cosmorec]
        for name, value in expected.items():
            self.assertAlmostEqual(getattr(rec, name), value)

        rec.set_params(recfast_approx_model=recombination.recfast_hyrec)
        rec.RECFAST_fudge_He = 0.85
        rec.Nz = 4096
        self.assertAlmostEqual(rec.RECFAST_fudge_He, 0.85)
        self.assertAlmostEqual(
            rec.AGauss1, recombination.recfast_approx_model_params[recombination.recfast_hyrec]["AGauss1"]
        )
        self.assertEqual(rec.Nz, 4096)

        pars = def_set_params(recfast_approx_model=recombination.recfast_planck)
        self.assertIsInstance(pars.Recomb, recombination.Recfast)
        pars.Recomb.RECFAST_fudge_He = 0.851
        self.assertAlmostEqual(pars.Recomb.RECFAST_fudge_He, 0.851)
        self.assertAlmostEqual(
            pars.Recomb.AGauss2, recombination.recfast_approx_model_params[recombination.recfast_planck]["AGauss2"]
        )

        with self.assertRaises(CAMBValueError):
            recombination.Recfast().set_params(recfast_approx_model="not_a_fit")

    def testErrors(self):
        redshifts = np.logspace(-1, np.log10(1089))
        pars = def_set_params(
            H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, redshifts=redshifts, kmax=0.1, silent=True
        )

        results = camb.get_background(pars)
        with self.assertRaises(CAMBError):
            results.get_matter_power_interpolator()

    def testEvolution(self):
        redshifts = [0.4, 31.5]
        pars = def_set_params(
            H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, redshifts=redshifts, kmax=0.1, silent=True
        )
        pars.WantCls = False

        # check transfer function routines and evolution code agree
        # Note transfer function redshifts are re-sorted in outputs
        data = camb.get_transfer_functions(pars)
        mtrans = data.get_matter_transfer_data()
        transfer_k = mtrans.transfer_z("delta_cdm", z_index=1)
        transfer_k2 = mtrans.transfer_z("delta_baryon", z_index=0)

        kh = mtrans.transfer_z("k/h", z_index=1)
        ev = data.get_redshift_evolution(
            mtrans.q, redshifts, ["delta_baryon", "delta_cdm", "delta_photon"], lAccuracyBoost=1
        )
        self.assertTrue(np.all(np.abs(transfer_k * kh**2 * (pars.H0 / 100) ** 2 / ev[:, 0, 1] - 1) < 1e-3))
        ix = 1
        self.assertAlmostEqual(transfer_k2[ix] * kh[ix] ** 2 * (pars.H0 / 100) ** 2, ev[ix, 1, 0], 4)

    def testInstances(self):
        pars = def_set_params(H0=69.1, ombh2=0.032, omch2=0.122, As=3e-9, ns=0.91, omk=0.013, redshifts=[0.0], kmax=0.5)
        data = camb.get_background(pars)
        res1 = data.angular_diameter_distance(0.7)
        drag1 = data.get_derived_params()["rdrag"]
        pars2 = def_set_params(H0=65, ombh2=0.022, omch2=0.122, As=3e-9, ns=0.91)
        data2 = camb.get_background(pars2)
        res2 = data2.angular_diameter_distance(1.7)
        drag2 = data2.get_derived_params()["rdrag"]
        self.assertAlmostEqual(res1, data.angular_diameter_distance(0.7))
        self.assertAlmostEqual(res2, data2.angular_diameter_distance(1.7))
        self.assertAlmostEqual(drag1, data.get_derived_params()["rdrag"])
        self.assertEqual(pars2.InitPower.ns, data2.Params.InitPower.ns)
        data2.calc_background(pars)
        self.assertEqual(pars.InitPower.ns, data2.Params.InitPower.ns)
        self.assertAlmostEqual(res1, data2.angular_diameter_distance(0.7))
        data3 = camb.get_results(pars2)
        cl3 = data3.get_lensed_scalar_cls(1000)
        self.assertAlmostEqual(res2, data3.angular_diameter_distance(1.7))
        self.assertAlmostEqual(drag2, data3.get_derived_params()["rdrag"], places=3)
        self.assertAlmostEqual(drag1, data.get_derived_params()["rdrag"], places=3)
        pars.set_for_lmax(3000, lens_potential_accuracy=1)
        camb.get_results(pars)
        del data3
        data4 = camb.get_results(pars2)
        cl4 = data4.get_lensed_scalar_cls(1000)
        np.testing.assert_allclose(cl4, cl3, atol=1e-20, rtol=1e-5)

    def testPowers(self):
        pars = new_def_params()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_dark_energy()  # re-set defaults
        pars.InitPower.set_params(ns=0.965, As=2e-9)
        pars.NonLinearModel.set_params(halofit_version="takahashi")

        self.assertAlmostEqual(pars.scalar_power(1), 1.801e-9, 4)
        self.assertAlmostEqual(pars.scalar_power([1, 1.5])[0], 1.801e-9, 4)

        pars.set_matter_power(nonlinear=True, silent=True)
        self.assertEqual(pars.NonLinear, model.NonLinear_pk)
        pars.set_matter_power(redshifts=[0.0, 0.17, 3.1], silent=True, nonlinear=False)
        data = camb.get_results(pars)

        kh, z, pk = data.get_matter_power_spectrum(1e-4, 1, 20)

        _kh2, _z2, pk2 = data.get_linear_matter_power_spectrum()

        s8 = data.get_sigma8()
        self.assertAlmostEqual(s8[0], 0.24686, 3)
        self.assertAlmostEqual(s8[2], 0.80044, 3)
        fs8 = data.get_fsigma8()
        self.assertAlmostEqual(fs8[0], 0.2431, 3)
        self.assertAlmostEqual(fs8[2], 0.424712, 3)

        pars.NonLinear = model.NonLinear_both

        data.calc_power_spectra(pars)
        _kh3, _z3, pk3 = data.get_matter_power_spectrum(1e-4, 1, 20)
        self.assertAlmostEqual(pk[-1][-3], 51.924, 2)
        self.assertAlmostEqual(pk3[-1][-3], 57.723, 2)
        self.assertAlmostEqual(pk2[-2][-4], 56.454, 2)
        camb.set_feedback_level(0)

        PKnonlin = camb.get_matter_power_interpolator(pars, nonlinear=True)
        pars.set_matter_power(
            redshifts=[0, 0.09, 0.15, 0.42, 0.76, 1.5, 2.3, 5.5, 8.9], silent=True, kmax=10, k_per_logint=5
        )
        pars.NonLinear = model.NonLinear_both
        results = camb.get_results(pars)
        kh, z, pk = results.get_nonlinear_matter_power_spectrum()
        pk_interp = PKnonlin.P(z, kh)
        self.assertTrue(np.sum((pk / pk_interp - 1) ** 2) < 0.005)
        PKnonlin2 = results.get_matter_power_interpolator(nonlinear=True, extrap_kmax=500)
        pk_interp2 = PKnonlin2.P(z, kh)
        self.assertTrue(np.sum((pk_interp / pk_interp2 - 1) ** 2) < 0.005)

        # The nonlinear matter grid can move when CMB source sampling changes; compare at a fixed physical k/h.
        def pk_at_fixed_kh(halofit_version, kh_value=0.44223076105117803):
            pars.NonLinearModel.set_params(halofit_version=halofit_version)
            kh, _, pk = results.get_nonlinear_matter_power_spectrum(params=pars)
            return np.exp(np.interp(np.log(kh_value), np.log(kh), np.log(pk[0])))

        self.assertAlmostEqual(pk_at_fixed_kh("mead"), 814.9, delta=0.5)

        self.assertAlmostEqual(pk_at_fixed_kh("mead2016"), 814.9, delta=0.5)

        self.assertAlmostEqual(pk_at_fixed_kh("mead2015"), 791.3, delta=0.5)

        self.assertAlmostEqual(pk_at_fixed_kh("mead2020"), 815.8, delta=0.5)

        self.assertAlmostEqual(pk_at_fixed_kh("mead2020_feedback"), 799.0, delta=0.5)

        lmax = 4000
        pars.set_for_lmax(lmax, lens_potential_accuracy=0)
        cls = data.get_cmb_power_spectra(pars)
        data.get_total_cls(2000)
        cls_unlensed = data.get_unlensed_scalar_cls(2500)
        data.get_tensor_cls(2000)
        cls_lensed = data.get_lensed_scalar_cls(3000)
        data.get_lens_potential_cls(2000)

        cls_lensed2 = data.get_lensed_cls_with_spectrum(data.get_lens_potential_cls()[:, 0], lmax=3000)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_lensed[2:, :], rtol=1e-4, atol=1e-18)
        cls_lensed2 = data.get_partially_lensed_cls(1, lmax=3000)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_lensed[2:, :], rtol=1e-4, atol=1e-18)
        cls_lensed2 = data.get_partially_lensed_cls(0, lmax=2500)
        np.testing.assert_allclose(cls_lensed2[2:, :], cls_unlensed[2:, :], rtol=1e-4, atol=1e-18)

        full_method = camb.lensing_method_curv_corr_full

        # Check lensed CL against python, including the same high-L template extension
        # used by the full Fortran curved-sky path.
        unlensed_scalar = data.get_unlensed_scalar_cls(data.Params.max_l)
        clpp = data.get_lens_potential_cls(data.Params.max_l)[:, 0]
        lens_lmax = data.Params.max_l - data.Params.lens_output_margin + 50 + 750
        if camb.config.AccuracyTarget > 0:
            boost = data.Params.Accuracy.AccuracyBoost * data.Params.Accuracy.LensingBoost
            lens_lmax = (
                data.Params.max_l
                - data.Params.lens_output_margin
                + 50
                + max(750, int(np.ceil((0.45 * (data.Params.max_l - data.Params.lens_output_margin) + 400) * boost)))
            )
        lens_lmax = min(8000, lens_lmax)
        cls_lensed2 = correlations.lensed_cls(
            unlensed_scalar,
            clpp,
            lmax=lens_lmax,
            lmax_lensed=3000,
            delta_cls=False,
            use_lensing_template=True,
            low_l_ee_taper=True,  # the Fortran short-range (AccurateBB=F) path tapers the kernel
        )
        cls_lensed_full = data.get_lensed_cls_with_spectrum(clpp, lmax=3000, lensing_method=full_method)
        np.testing.assert_allclose(cls_lensed2[2:3000, 2], cls_lensed_full[2:3000, 2], rtol=1e-6, atol=1e-18)
        np.testing.assert_allclose(cls_lensed2[2:3000, 1], cls_lensed_full[2:3000, 1], rtol=1e-6, atol=1e-18)
        np.testing.assert_allclose(cls_lensed2[2:3000, 0], cls_lensed_full[2:3000, 0], rtol=1e-6, atol=1e-18)
        self.assertTrue(
            np.all(
                np.abs(
                    (cls_lensed2[2:3000, 3] - cls_lensed_full[2:3000, 3])
                    / np.sqrt(cls_lensed2[2:3000, 0] * cls_lensed2[2:3000, 1])
                )
                < 1e-8
            )
        )

        optimized_method = camb.lensing_method_optimized
        pars = new_def_params()
        pars.set_cosmology(H0=67)
        pars.set_for_lmax(2500, lens_potential_accuracy=1)
        pars.Accuracy.AccurateBB = True
        data = camb.get_results(pars)
        clpp = data.get_lens_potential_cls()[:, 0]
        cls_lensed_full = data.get_lensed_cls_with_spectrum(clpp, lmax=2500, lensing_method=full_method)
        cls_lensed_optimized = data.get_lensed_cls_with_spectrum(clpp, lmax=2500, lensing_method=optimized_method)
        np.testing.assert_allclose(cls_lensed_optimized[2:, :], cls_lensed_full[2:, :], rtol=1e-12)

        pars = new_def_params()
        pars.set_cosmology(H0=67)
        pars.set_for_lmax(2500, lens_potential_accuracy=1)
        pars.Accuracy.AccurateBB = False
        data = camb.get_results(pars)
        clpp = data.get_lens_potential_cls()[:, 0]
        original_method = camb.config.lensing_method
        cls_lensed_curv = data.get_lensed_cls_with_spectrum(
            clpp, lmax=2500, lensing_method=camb.lensing_method_curv_corr
        )
        cls_lensed_optimized = data.get_lensed_cls_with_spectrum(clpp, lmax=2500, lensing_method=optimized_method)
        np.testing.assert_allclose(cls_lensed_optimized[2:, :], cls_lensed_curv[2:, :], rtol=1e-7)
        self.assertEqual(camb.config.lensing_method, original_method)

        corr, xvals, weights = correlations.gauss_legendre_correlation(cls["lensed_scalar"])
        clout = correlations.corr2cl(corr, xvals, weights, 2500)
        self.assertTrue(np.all(np.abs(clout[2:2300, 2] / cls["lensed_scalar"][2:2300, 2] - 1) < 1e-3))

        pars = new_def_params()
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
            cls += [results.get_total_cls(CMB_unit="muK")]
        np.testing.assert_allclose((cls[1] - cls[0])[2:300, 2] * 2, (cls[2] - cls[0])[2:300, 2], rtol=1e-3)

        # Check generating tensors and scalars together
        pars = new_def_params()
        pars.set_cosmology(H0=67)
        lmax = 2000
        pars.set_for_lmax(lmax, lens_potential_accuracy=1)
        pars.InitPower.set_params(ns=0.96, r=0)
        pars.WantTensors = False
        results = camb.get_results(pars)
        cl1 = results.get_total_cls(lmax, CMB_unit="muK")
        pars.InitPower.set_params(ns=0.96, r=0.1, nt=0)
        pars.WantTensors = True
        results = camb.get_results(pars)
        cl2 = results.get_lensed_scalar_cls(lmax, CMB_unit="muK")
        ctensor2 = results.get_tensor_cls(lmax, CMB_unit="muK")
        results = camb.get_transfer_functions(pars)
        results.Params.InitPower.set_params(ns=1.1, r=1)
        inflation_params = initialpower.InitialPowerLaw()
        inflation_params.set_params(ns=0.96, r=0.05, nt=0)
        results.power_spectra_from_transfer(inflation_params, silent=True)
        cl3 = results.get_lensed_scalar_cls(lmax, CMB_unit="muK")
        ctensor3 = results.get_tensor_cls(lmax, CMB_unit="muK")
        np.testing.assert_allclose(ctensor2, ctensor3 * 2, rtol=1e-4)
        np.testing.assert_allclose(cl1, cl2, rtol=1e-4)
        # These are identical because all scalar spectra were identical (non-linear corrections change it  otherwise)
        np.testing.assert_allclose(cl1, cl3, rtol=1e-4)

        pars = new_def_params()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.07, omk=0)
        pars.set_for_lmax(2500, lens_potential_accuracy=0)
        pars.min_l = 2
        res = camb.get_results(pars)
        cls = res.get_lensed_scalar_cls(2000)
        pars.min_l = 1
        res = camb.get_results(pars)
        cls2 = res.get_lensed_scalar_cls(2000)
        np.testing.assert_allclose(cls[2:, 0:2], cls2[2:, 0:2], rtol=1e-4)
        np.testing.assert_allclose(cls2[1, 0], 1.303942e-10, rtol=3e-3)
        self.assertAlmostEqual(cls[1, 0], 0)

    def testSave(self):
        pars = def_set_params(
            H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, redshifts=[0.4, 31.5], kmax=0.1, silent=True
        )
        pars.set_dark_energy(w=-0.7, wa=0.2, dark_energy_model="ppf")
        from camb.sources import GaussianSourceWindow

        pars.SourceWindows = [GaussianSourceWindow(), GaussianSourceWindow(redshift=1)]
        s = repr(pars)
        pars2 = eval(s)
        assert repr(pars2) == s
        assert "DarkEnergyPPF" in str(pars2)
        b = pickle.dumps(pars)
        pars2 = pickle.loads(b)
        assert repr(pars2) == s
        pars2.InitPower = initialpower.SplinedInitialPower()
        with self.assertRaises(TypeError):
            repr(pars2)

    def testSigmaR(self):
        pars = new_def_params()
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
        P, _z, k = results.get_matter_power_interpolator(
            nonlinear=False, hubble_units=False, k_hunit=False, return_z_k=True, extrap_kmax=100, silent=True
        )
        truth = 0.800679  # from high kmax, high accuracy boost
        self.assertTrue(abs(results.get_sigmaR(8)[-1] / sigma8 - 1) < 1e-3)

        def get_sigma(_ks, dlogk):
            x = _ks * 8 / (pars.H0 / 100)
            w = (3 * (np.sin(x) - x * np.cos(x)) / x**3) ** 2
            w[x < 1e-2] = 1 - x[x < 1e-2] ** 2 / 2
            Ps = P.P(0, _ks) * _ks**3 / (2 * np.pi**2)
            return np.sqrt(np.dot(w, Ps * dlogk))

        logk = np.arange(np.log(1e-5), np.log(20.0), 1.0 / 100)
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
        pars.set_matter_power(nonlinear=False, k_per_logint=0, kmax=1.2, redshifts=np.arange(0, 10, 2), silent=True)
        results = camb.get_results(pars)
        sigmas = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        pars.Accuracy.AccuracyBoost = 2
        results = camb.get_results(pars)
        sigmas2 = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        np.testing.assert_allclose(sigmas, sigmas2, rtol=1e-3)
        pars.Accuracy.AccuracyBoost = 1
        pars.set_matter_power(nonlinear=False, k_per_logint=100, kmax=10, redshifts=np.arange(0, 10, 2), silent=True)
        results = camb.get_results(pars)
        sigmas2 = results.get_sigmaR(np.arange(1, 20, 1), hubble_units=False, z_indices=None)
        self.assertAlmostEqual(sigmas2[4, 2], 1.77346, places=3)
        self.assertTrue(np.all(np.abs(sigmas[:, 1:] / sigmas2[:, 1:] - 1) < 1e-3))
        self.assertTrue(np.all(np.abs(sigmas[:, 0] / sigmas2[:, 0] - 1) < 2e-3))

    def testTimeTransfers(self):
        from camb import initialpower

        pars = def_set_params(H0=69, YHe=0.22, lmax=2000, lens_potential_accuracy=1, ns=0.96, As=2.5e-9)
        results1 = camb.get_results(pars)
        cl1 = results1.get_total_cls()

        pars = def_set_params(H0=69, YHe=0.22, lmax=2000, lens_potential_accuracy=1)
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

        pars = new_def_params()
        pars.set_cosmology(H0=78, YHe=0.22)
        pars.set_for_lmax(2000, lens_potential_accuracy=1)
        pars.WantTensors = True
        results = camb.get_transfer_functions(pars, only_time_sources=True)
        cls = []
        for r in [0, 0.2, 0.4]:
            inflation_params = initialpower.InitialPowerLaw()
            inflation_params.set_params(ns=0.96, r=r, nt=0)
            results.power_spectra_from_transfer(inflation_params)
            cls += [results.get_total_cls(CMB_unit="muK")]
        np.testing.assert_allclose((cls[1] - cls[0])[2:300, 2] * 2, (cls[2] - cls[0])[2:300, 2], rtol=1e-3)

    def testDarkEnergy(self):
        pars = new_def_params()
        pars.set_cosmology(H0=71)
        pars.InitPower.set_params(ns=0.965, r=0)
        for m in ["fluid", "ppf"]:
            pars.set_dark_energy(w=-0.7, wa=0.2, dark_energy_model=m)
            C1 = camb.get_results(pars).get_cmb_power_spectra()
            a = np.logspace(-5, 0, 1000)
            w = -0.7 + 0.2 * (1 - a)
            pars2 = pars.copy()
            pars2.set_dark_energy_w_a(a, w, dark_energy_model=m)
            C2 = camb.get_results(pars2).get_cmb_power_spectra()
            for f in ["lens_potential", "lensed_scalar"]:
                np.testing.assert_allclose(C1[f][2:, 0], C2[f][2:, 0], rtol=1e-4, atol=5e-14)
            pars3 = pars2.copy()
            self.assertAlmostEqual(-0.7, pars3.DarkEnergy.w)

    def testInitialPower(self):
        pars = new_def_params()
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

        pars2 = new_def_params()
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
        sp.set_scalar_log_regular(10 ** (-5.5), 10.0**2, pk)
        pars2.set_initial_power(sp)
        self.assertAlmostEqual(pars2.scalar_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)

        sp.set_tensor_log_regular(10 ** (-5.5), 10.0**2, pk)
        pars2.set_initial_power(sp)
        self.assertAlmostEqual(pars2.tensor_power(1.1), pars.scalar_power(1.1), delta=As * 1e-4)
        self.assertTrue(sp.has_tensors())
        sp.set_tensor_table([], [])
        self.assertFalse(sp.has_tensors())
        pars2.set_initial_power(sp)

        results = camb.get_results(pars2)
        cl = results.get_lensed_scalar_cls(CMB_unit="muK")
        pars.InitPower.set_params(As=As, ns=ns)
        results2 = camb.get_results(pars)
        cl2 = results2.get_lensed_scalar_cls(CMB_unit="muK")
        np.testing.assert_allclose(cl, cl2, rtol=1e-4)
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

        pars = new_def_params()
        pars.set_cosmology(H0=64, mnu=0)
        pars.set_for_lmax(1200, lens_potential_accuracy=0)
        pars.Want_CMB = False
        pars.SourceWindows = [
            GaussianSourceWindow(redshift=0.17, source_type="counts", bias=1.2, sigma=0.04, dlog10Ndm=-0.2),
            GaussianSourceWindow(redshift=0.5, source_type="lensing", sigma=0.07, dlog10Ndm=0),
        ]
        pars.SourceTerms.limber_windows = True
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        zs = np.arange(0, 0.5, 0.02)
        W = np.exp(-((zs - 0.17) ** 2) / 2 / 0.04**2) / np.sqrt(2 * np.pi) / 0.04

        ks = np.logspace(-4, 3, 50)
        bias_kz = 1.2 * np.ones((len(ks), len(zs)))
        test_windows = [
            SplinedSourceWindow(bias=1.2, dlog10Ndm=-0.2, z=zs, W=W),
            SplinedSourceWindow(bias_z=1.2 * np.ones_like(zs), dlog10Ndm=-0.2, z=zs, W=W),
            SplinedSourceWindow(k_bias=ks, bias_kz=bias_kz, dlog10Ndm=-0.2, z=zs, W=W),
        ]
        for window in test_windows:
            pars.SourceWindows[0] = window
            results = camb.get_results(pars)
            cls2 = results.get_source_cls_dict()
            np.testing.assert_allclose(cls2["W1xW1"][2:1200], cls["W1xW1"][2:1200], rtol=1e-3)

        pars.SourceWindows = [GaussianSourceWindow(redshift=1089, source_type="lensing", sigma=30)]
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        PP = cls["PxP"]
        ls = np.arange(0, PP.shape[0])
        np.testing.assert_allclose(PP / 4 * (ls * (ls + 1)), cls["W1xW1"], rtol=1e-3, atol=1e-10)
        np.testing.assert_allclose(PP / 2 * np.sqrt(ls * (ls + 1)), cls["PxW1"], rtol=1e-3, atol=1e-10)
        # test something sharp with redshift distortions (tricky..)
        from scipy import signal

        zs = np.arange(1.9689, 2.1057, (2.1057 - 1.9689) / 2000)
        W = signal.windows.tukey(len(zs), alpha=0.1)
        pars = new_def_params()
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
        pars.InitPower.set_params(As=2e-9, ns=0.965)
        pars.set_for_lmax(4000, lens_potential_accuracy=0)
        pars.SourceWindows = [SplinedSourceWindow(z=zs, W=W, source_type="counts")]
        pars.SourceTerms.counts_redshift = True
        results = camb.get_results(pars)
        cls = results.get_source_cls_dict()
        self.assertAlmostEqual(np.sum(cls["PxW1"][10:3000:20]), 0.00020001, places=5)
        self.assertAlmostEqual(np.sum(cls["W1xW1"][10:3000:20]), 2.26350, places=3)
        self.assertAlmostEqual(np.sum(cls["W1xW1"][10]), 0.0001097, places=6)

    def test_memory(self):
        if platform.system() != "Windows":
            import gc
            import resource

            last_usage = -1
            for i in range(3):
                pars = new_def_params()
                pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.12, mnu=0.06, omk=0, tau=0.17)
                results = camb.get_results(pars)
                del pars, results
                gc.collect()
                usage = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1)
                if 0 < last_usage != usage:
                    print(f"Memory usage: {usage:2.2f} KB vs {last_usage:2.2f} KB")
                    raise Exception("Apparent memory leak")
                last_usage = usage

        camb.free_global_memory()

    def test_quintessence(self):
        n = 3
        # set zc and fde_zc
        pars = def_set_params(
            ombh2=0.022,
            omch2=0.122,
            thetastar=0.01044341764253,
            dark_energy_model="EarlyQuintessence",
            m=8e-53,
            f=0.05,
            n=n,
            theta_i=3.1,
            use_zc=True,
            zc=1e4,
            fde_zc=0.1,
        )
        camb.get_background(pars)
        results = camb.get_results(pars)
        self.assertAlmostEqual(results.get_derived_params()["thetastar"], 1.044341764253, delta=1e-5)
