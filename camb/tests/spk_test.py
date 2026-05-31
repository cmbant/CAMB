import inspect
import os
import sys
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb

from camb.nonlinear import Halofit, SPkNonLinear  # type: ignore[attr-defined]


class SPkTest(unittest.TestCase):
    def _get_pk(self, model_obj, z=0.5, kmax=5.0):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=[z], kmax=kmax, k_per_logint=100)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel = model_obj
        results = camb.get_results(pars)
        kh, _z, pk = results.get_nonlinear_matter_power_spectrum()
        return kh, pk[0], results

    def test_spk_invalid_params(self):
        model = SPkNonLinear()
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_SO=300)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=99)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=1, SPk_fb_pivot=0.0)

    def test_spk_accepts_halofit_version(self):
        model = SPkNonLinear()
        model.set_params(
            halofit_version="mead2016",
            SPk_feedback=True,
            SPk_SO=200,
            SPk_relation_kind=1,
            SPk_fb_a=0.4,
            SPk_fb_pow=0.2,
            SPk_fb_pivot=1e14,
        )
        k, pk, _ = self._get_pk(model, z=0.5, kmax=3.0)
        self.assertTrue(np.all(np.isfinite(k)))
        self.assertTrue(np.all(np.isfinite(pk)))

    def test_spk_rejects_mead2020_feedback_via_set_params(self):
        model = SPkNonLinear()
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_feedback=True, halofit_version="mead2020_feedback")

    def test_spk_cobaya_friendly_set_params_signature(self):
        signature = inspect.signature(SPkNonLinear.set_params)
        self.assertIn("halofit_version", signature.parameters)
        self.assertIn("HMCode_A_baryon", signature.parameters)
        self.assertIn("HMCode_eta_baryon", signature.parameters)
        self.assertIn("HMCode_logT_AGN", signature.parameters)

    def test_spk_hmcode_params_passed_to_base(self):
        """HMCode parameters set via SPkNonLinear.set_params() are forwarded to BaseModel."""
        model = SPkNonLinear()
        model.set_params(
            halofit_version="mead2020_feedback",
            HMCode_A_baryon=3.5,
            HMCode_eta_baryon=0.7,
            HMCode_logT_AGN=8.0,
        )
        self.assertAlmostEqual(model.BaseModel.HMCode_A_baryon, 3.5)
        self.assertAlmostEqual(model.BaseModel.HMCode_eta_baryon, 0.7)
        self.assertAlmostEqual(model.BaseModel.HMCode_logT_AGN, 8.0)

    def test_spk_disabled_matches_base(self):
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k_base, pk_base, _ = self._get_pk(base)

        spk = SPkNonLinear()
        spk.set_params(halofit_version="mead2020", SPk_feedback=False)
        k_spk, pk_spk, _ = self._get_pk(spk)

        self.assertTrue(np.allclose(k_base, k_spk, rtol=0, atol=0))
        self.assertTrue(np.allclose(pk_base, pk_spk, rtol=2e-12, atol=1e-14))

    def test_spk_out_of_range_behaviour(self):
        """Verify suppression is skipped for z outside calibrated range and k is clamped."""
        # z=4 is beyond calibrated range [0, 3]: suppression should not be applied.
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k_base, pk_base_z4, _ = self._get_pk(base, z=4.0, kmax=20.0)

        spk = SPkNonLinear()
        spk.set_params(
            halofit_version="mead2020",
            SPk_feedback=True,
            SPk_SO=200,
            SPk_relation_kind=1,
            SPk_fb_a=0.4,
            SPk_fb_pow=0.3,
            SPk_fb_pivot=1e14,
        )
        k_spk, pk_spk_z4, _ = self._get_pk(spk, z=4.0, kmax=20.0)

        # At z=4 (out of range), SPk should be identity — P(k) unchanged.
        np.testing.assert_allclose(pk_spk_z4, pk_base_z4, rtol=1e-10)

        # At z=0.5 (in range) with k up to 20, k > 12 is clamped — suppression still applied.
        _, pk_base_z05, _ = self._get_pk(base, z=0.5, kmax=20.0)
        _, pk_spk_z05, _ = self._get_pk(spk, z=0.5, kmax=20.0)
        sup = pk_spk_z05 / pk_base_z05
        # Suppression should differ from 1 for k in calibrated range.
        mask = (k_spk > 0.1) & (k_spk <= 12.0)
        self.assertFalse(np.allclose(sup[mask], 1.0, atol=1e-4))

    def test_spk_off_node_boundary_regression(self):
        """Off-node z boundary cases should remain finite with Akima limit interpolation.

        At off-node redshifts the linear-interpolated min_fb is slightly higher
        than the Akima-interpolated value.  These test points sit in the gap:
        Akima says in-range (finite P(k)), linear would say out-of-range (NaN).

        Uses get_linear_matter_power_spectrum(nonlinear=True) which returns raw
        per-(z, k) nonlinear P(k) without spline smoothing — the only accessor
        that faithfully preserves per-cell NaN (see docs/source/spk.rst).
        """
        # (z, k_target, fb_a) tuples where fb_a sits between Akima and linear min_fb.
        cases = [
            (1.5, 0.5, 0.15399),
            (2.5, 11.5, 0.1805),
        ]

        # Limit coefficient arrays for SO=200 (from spk_model.f90).
        z_nodes = np.array([0.0, 0.125, 0.5, 1.0, 2.0, 3.0])
        min_x0_200 = np.array(
            [
                63.59373179416563,
                59.88726319810792,
                56.365373020207954,
                39.64033211739476,
                91.48777680660496,
                45.013496639467114,
            ]
        )
        min_x1_200 = np.array(
            [
                -9.731727022847117,
                -9.176876134517682,
                -8.677101127391419,
                -6.141984569473165,
                -14.545324008239655,
                -7.155837194116757,
            ]
        )
        min_x2_200 = np.array(
            [
                0.36717360571115487,
                0.34646698848026913,
                0.32901138538950075,
                0.23315608004243354,
                0.5737424941339003,
                0.280547072910215,
            ]
        )

        for z, k_target, fb_a in cases:
            with self.subTest(z=z, k=k_target, fb_a=fb_a):
                pars = camb.CAMBparams()
                pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
                pars.InitPower.set_params(As=2.1e-9, ns=0.965)
                pars.set_matter_power(redshifts=[z], kmax=12.0, k_per_logint=80)
                pars.NonLinear = camb.model.NonLinear_both
                pars.NonLinearModel = SPkNonLinear()
                pars.NonLinearModel.set_params(
                    halofit_version="mead2020",
                    SPk_feedback=True,
                    SPk_SO=200,
                    SPk_relation_kind=1,
                    SPk_fb_a=fb_a,
                    SPk_fb_pow=0.0,
                    SPk_fb_pivot=1.0e14,
                )

                data = camb.get_results(pars)
                kh, zs, pk = data.get_linear_matter_power_spectrum(nonlinear=True)
                self.assertEqual(len(zs), 1)

                idx = np.argmin(np.abs(kh - k_target))
                self.assertLess(abs(kh[idx] - k_target), 0.5)
                self.assertTrue(
                    np.isfinite(pk[0, idx]),
                    f"Expected finite P(k) at z={z}, k/h~{k_target} (Akima in-range)",
                )

                # Verify that the old linear limit would have rejected this point.
                x = 1.0 + z
                spk_a = 15.24311120000861 - 1.2436699435560352 * x + 0.14837558774401766 * x * x
                spk_b = 14.969187892657688 - 1.0993025612653198 * x + 0.12905587245129102 * x * x
                spk_g = 0.8000441576980428 - 0.01715621131893159 * x + 0.06131887249968379 * x * x
                best_mass = spk_a - (spk_a - spk_b) * (k_target**spk_g)
                m_opt = 10.0**best_mass
                logm = np.log10(m_opt)
                min_c0_lin = np.interp(z, z_nodes, min_x0_200)
                min_c1_lin = np.interp(z, z_nodes, min_x1_200)
                min_c2_lin = np.interp(z, z_nodes, min_x2_200)
                min_fb_lin = 0.8 * 10.0 ** (min_c0_lin + min_c1_lin * logm + min_c2_lin * logm * logm)
                self.assertLess(fb_a, min_fb_lin)

    def test_spk_fb_outside_calibrated_limits_produces_nan(self):
        """When fb is pushed far outside the calibrated fitting limits, P(k) should contain NaN."""
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=[0.5], kmax=3.0)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel = SPkNonLinear()
        pars.NonLinearModel.set_params(
            halofit_version="mead2020",
            SPk_feedback=True,
            SPk_SO=200,
            SPk_relation_kind=1,
            SPk_fb_a=100.0,
            SPk_fb_pow=0.0,
            SPk_fb_pivot=1.0,
        )
        data = camb.get_results(pars)
        interp = data.get_matter_power_interpolator(nonlinear=True)
        k_test = np.logspace(-1, 0.4, 10)
        pk = interp.P(0.5, k_test)  # type: ignore[union-attr]
        self.assertTrue(np.any(np.isnan(pk)), "Expected NaN for out-of-limits fb")

        # Verify valid params remain finite
        pars2 = camb.CAMBparams()
        pars2.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars2.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars2.set_matter_power(redshifts=[0.5], kmax=3.0)
        pars2.NonLinear = camb.model.NonLinear_both
        pars2.NonLinearModel = SPkNonLinear()
        pars2.NonLinearModel.set_params(
            halofit_version="mead2020",
            SPk_feedback=True,
            SPk_fb_a=0.4,
            SPk_fb_pow=0.2,
            SPk_fb_pivot=1e14,
        )
        data2 = camb.get_results(pars2)
        interp2 = data2.get_matter_power_interpolator(nonlinear=True)
        pk2 = interp2.P(0.5, k_test)  # type: ignore[union-attr]
        self.assertTrue(np.all(np.isfinite(pk2)), "Expected finite P(k) for valid fb")

    def test_spk_class_selection_via_set_classes(self):
        pars = camb.CAMBparams()
        pars.set_classes(non_linear_model="SPkNonLinear")
        self.assertEqual(pars.NonLinearModel.__class__.__name__, "SPkNonLinear")

        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=[0.5], kmax=3.0)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel.set_params(
            halofit_version="mead2020",
            SPk_feedback=True,
            SPk_SO=200,
            SPk_relation_kind=1,
            SPk_fb_a=0.4,
            SPk_fb_pow=0.2,
            SPk_fb_pivot=1e14,
        )

        data = camb.get_results(pars)
        k, z, pk = data.get_matter_power_spectrum(minkh=1e-2, maxkh=1.0, npoints=8)
        self.assertEqual(len(z), 1)
        self.assertTrue(np.all(np.isfinite(k)))
        self.assertTrue(np.all(np.isfinite(pk)))


if __name__ == "__main__":
    unittest.main()
