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
