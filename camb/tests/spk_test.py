import os
import sys
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb

try:
    import pyspk
except ImportError:
    pyspk = None


class _CambEFunc:
    def __init__(self, results):
        self._results = results
        self._h0 = results.Params.H0

    def efunc(self, z):
        z_arr = np.asarray(z, dtype=float)
        flat = z_arr.reshape(-1)
        vals = np.array([self._results.hubble_parameter(float(zz)) / self._h0 for zz in flat], dtype=float)
        vals = vals.reshape(z_arr.shape)
        if np.isscalar(z):
            return float(vals)
        return vals


class SPkTest(unittest.TestCase):
    def _get_pk(self, model_obj, z=0.5):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        pars.set_matter_power(redshifts=[z], kmax=6.0)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel = model_obj
        results = camb.get_results(pars)
        k, _z, pk = results.get_matter_power_spectrum(minkh=1e-2, maxkh=5.0, npoints=28)
        return k, pk[0], results

    def _assert_relation_match(self, relation_kind, z=0.5, so=200, **kwargs):
        base = camb.Halofit()
        base.set_params(halofit_version="mead2020")
        k, pk_base, results = self._get_pk(base, z=z)

        spk = camb.SPkNonLinear()
        spk.BaseModel.set_params(halofit_version="mead2020")
        spk.set_params(SPk_feedback=True, SPk_SO=so, SPk_relation_kind=relation_kind, **kwargs)
        _, pk_spk, _ = self._get_pk(spk, z=z)

        cosmo = _CambEFunc(results)
        pyspk_kwargs = {
            "SO": so,
            "z": z,
            "k_array": k,
            "errors": False,
        }
        if relation_kind == 1:
            pyspk_kwargs.update(
                {
                    "fb_a": kwargs["SPk_fb_a"],
                    "fb_pow": kwargs["SPk_fb_pow"],
                    "fb_pivot": kwargs["SPk_fb_pivot"],
                }
            )
        elif relation_kind == 2:
            pyspk_kwargs.update(
                {
                    "alpha": kwargs["SPk_alpha"],
                    "beta": kwargs["SPk_beta"],
                    "gamma": kwargs["SPk_gamma"],
                    "cosmo": cosmo,
                }
            )
        elif relation_kind == 3:
            pyspk_kwargs.update(
                {
                    "epsilon": kwargs["SPk_epsilon"],
                    "alpha": kwargs["SPk_alpha"],
                    "beta": kwargs["SPk_beta"],
                    "gamma": kwargs["SPk_gamma"],
                    "m_pivot": kwargs["SPk_m_pivot"],
                    "cosmo": cosmo,
                }
            )
        else:
            raise ValueError("Unknown relation kind")

        k_ref, sup_ref = pyspk.sup_model(**pyspk_kwargs)
        measured_sup = pk_spk / pk_base

        self.assertTrue(np.allclose(k_ref, k, rtol=0, atol=0))
        self.assertTrue(np.allclose(measured_sup, sup_ref, rtol=1e-5, atol=1e-8))

    def test_spk_invalid_params(self):
        model = camb.SPkNonLinear()
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_SO=300)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=99)
        with self.assertRaises(camb.CAMBValueError):
            model.set_params(SPk_relation_kind=1, SPk_fb_pivot=0.0)

    def test_spk_disabled_matches_base(self):
        base = camb.Halofit()
        base.set_params(halofit_version="mead2020")
        k_base, pk_base, _ = self._get_pk(base)

        spk = camb.SPkNonLinear()
        spk.BaseModel.set_params(halofit_version="mead2020")
        spk.set_params(SPk_feedback=False)
        k_spk, pk_spk, _ = self._get_pk(spk)

        self.assertTrue(np.allclose(k_base, k_spk, rtol=0, atol=0))
        self.assertTrue(np.allclose(pk_base, pk_spk, rtol=2e-12, atol=1e-14))

    @unittest.skipIf(pyspk is None, "pyspk not installed")
    def test_spk_power_law_matches_reference(self):
        self._assert_relation_match(
            relation_kind=1,
            z=0.5,
            so=200,
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
            SPk_epsilon=0.55,
            SPk_alpha=0.2,
            SPk_beta=0.9,
            SPk_gamma=0.1,
            SPk_m_pivot=1e14,
        )


if __name__ == "__main__":
    unittest.main()
