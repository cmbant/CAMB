import inspect
import os
import sys
import types
import unittest
import warnings

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb

from camb.nonlinear import Halofit, SPkNonLinear  # type: ignore[attr-defined]

from _spk_helpers import get_pk, pyspk_kwargs

try:
    import pyspk
except ImportError:
    pyspk: types.ModuleType | None = None  # type: ignore[no-redef]


class SPkTest(unittest.TestCase):
    def _get_pk(self, model_obj, z=0.5, kmax=5.0):
        return get_pk(model_obj, z=z, kmax=kmax)

    def _assert_relation_match(self, relation_kind, z=0.5, so=200, rtol=1e-5, atol=1e-8, **kwargs):
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k, pk_base, results = self._get_pk(base, z=z)

        spk = SPkNonLinear()
        spk.set_params(halofit_version="mead2020", SPk_feedback=True, SPk_SO=so, SPk_relation_kind=relation_kind, **kwargs)
        _, pk_spk, _ = self._get_pk(spk, z=z)

        assert pyspk is not None
        py_kwargs = pyspk_kwargs(relation_kind, so, z, k, results, kwargs)
        k_ref, sup_ref = pyspk.sup_model(**py_kwargs)
        measured_sup = pk_spk / pk_base

        self.assertTrue(np.allclose(k_ref, k, rtol=0, atol=0))
        self.assertTrue(np.allclose(measured_sup, sup_ref, rtol=rtol, atol=atol))

    def test_spk_invalid_params(self):
        model = SPkNonLinear()
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_SO=300)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=99)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=1, SPk_fb_pivot=0.0)

    def test_spk_halofit_version_passthrough(self):
        model = SPkNonLinear()
        model.set_params(halofit_version="mead2016")
        self.assertEqual(model.BaseModel.get_halofit_version(), "mead2016")

    def test_spk_rejects_mead2020_feedback_via_set_params(self):
        model = SPkNonLinear()
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_feedback=True, halofit_version="mead2020_feedback")

    def test_spk_cobaya_friendly_set_params_signature(self):
        signature = inspect.signature(SPkNonLinear.set_params)
        self.assertIn("halofit_version", signature.parameters)
        self.assertNotIn("base_model", signature.parameters)

    def test_spk_disabled_matches_base(self):
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k_base, pk_base, _ = self._get_pk(base)

        spk = SPkNonLinear()
        spk.set_params(halofit_version="mead2020", SPk_feedback=False)
        k_spk, pk_spk, _ = self._get_pk(spk)

        self.assertTrue(np.allclose(k_base, k_spk, rtol=0, atol=0))
        self.assertTrue(np.allclose(pk_base, pk_spk, rtol=2e-12, atol=1e-14))

    @unittest.skipIf(pyspk is None, "pyspk not installed")
    def test_spk_power_law_matches_reference(self):
        self._assert_relation_match(
            relation_kind=1,
            z=0.5,
            so=200,
            rtol=1e-6,
            atol=1e-9,
            SPk_fb_a=0.4,
            SPk_fb_pow=0.2,
            SPk_fb_pivot=1e14,
        )

    @unittest.skipIf(pyspk is None, "pyspk not installed")
    def test_spk_cosmo_power_law_matches_reference(self):
        self._assert_relation_match(
            relation_kind=2,
            z=1.0,
            so=500,
            rtol=1e-6,
            atol=1e-9,
            SPk_alpha=3.4,
            SPk_beta=1.0,
            SPk_gamma=0.15,
        )

    @unittest.skipIf(pyspk is None, "pyspk not installed")
    def test_spk_double_power_law_matches_reference(self):
        self._assert_relation_match(
            relation_kind=3,
            z=1.5,
            so=200,
            rtol=1e-6,
            atol=1e-9,
            SPk_epsilon=0.55,
            SPk_alpha=0.2,
            SPk_beta=0.9,
            SPk_gamma=0.1,
            SPk_m_pivot=1e14,
        )

    @unittest.skipIf(pyspk is None, "pyspk not installed")
    def test_spk_high_k_boundary_continuity(self):
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k, pk_base, _ = self._get_pk(base, z=0.125, kmax=12.0)

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
        _, pk_spk, _ = self._get_pk(spk, z=0.125, kmax=12.0)

        # Mask to pySPK's calibrated range (CAMB's raw grid may exceed kmax slightly)
        mask = k < 12.0
        k, pk_base, pk_spk = k[mask], pk_base[mask], pk_spk[mask]

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=R"Scales with k_max > k_ny = 8\.0 \[h/Mpc\] may not be accurately reproduced by the model\.",
                category=UserWarning,
            )
            assert pyspk is not None
            _, sup_ref = pyspk.sup_model(
                SO=200,
                z=0.125,
                fb_a=0.4,
                fb_pow=0.3,
                fb_pivot=1e14,
                k_array=k,
                errors=False,
            )
        measured_sup = pk_spk / pk_base
        rel = np.abs(measured_sup / sup_ref - 1.0)

        self.assertLess(float(np.max(rel)), 1e-4)

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
