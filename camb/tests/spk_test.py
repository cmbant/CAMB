import os
import sys
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb


def _poly2(x, c0, c1, c2):
    return c2 * x * x + c1 * x + c0


def _spk_params(so, z):
    x = 1.0 + z
    if so == 500:
        return {
            "a": _poly2(x, 14.783423122120318, -0.999062404857228, 0.12062854541689262),
            "b": _poly2(x, 14.620528368613265, -0.9136466201011957, 0.10835389086945699),
            "g": _poly2(x, 0.9671320682693298, -0.03185388045484575, 0.02650236152450093),
            "lambda_a": _poly2(x, 0.019349810078190303, -0.007410668383424459, 0.0008334762393555539),
            "lambda_b": _poly2(x, 2.9566773924238143, 0.6205340408676114, -0.001928273640110775),
            "mu_a": _poly2(x, 0.715853343781141, -0.19276613600825665, 0.04948240117059147),
            "mu_b": _poly2(x, 3.385355123440431, 0.9658906605139421, -0.06825861100375574),
            "mu_c": _poly2(x, 4.457257708010122, -2.191853871334233, 0.45457701107254733),
            "nu_a": _poly2(x, 478.86477329610375, 429.88795783439946, 249.25655627821902),
            "nu_b": _poly2(x, -11.227459319819815, -0.5581080204509223, 0.4489962047114509),
            "nu_c": _poly2(x, 3.499449440557995, -0.08488559389068073, -0.0923847866118189),
        }
    return {
        "a": _poly2(x, 15.24311120000861, -1.2436699435560352, 0.14837558774401766),
        "b": _poly2(x, 14.969187892657688, -1.0993025612653198, 0.12905587245129102),
        "g": _poly2(x, 0.8000441576980428, -0.01715621131893159, 0.06131887249968379),
        "lambda_a": _poly2(x, 0.02178116280689233, -0.0077564325654746955, 0.0007915576054589781),
        "lambda_b": _poly2(x, 3.0878286643613437, 0.4529677646796634, 0.001552571083240605),
        "mu_a": _poly2(x, 0.6930259177449359, -0.16913553700233935, 0.04263185199898842),
        "mu_b": _poly2(x, 3.161914061444856, 0.8616834297321924, 0.011346427353554053),
        "mu_c": _poly2(x, 5.532188503256583, -3.0864672185252537, 0.5083422518560442),
        "nu_a": _poly2(x, 413.00988701513904, 311.63957063032285, 37.89105940901369),
        "nu_b": _poly2(x, -11.243859405779181, -0.34421412616421965, 0.3343548325485801),
        "nu_c": _poly2(x, 3.476463891168505, -0.018333059687988575, -0.08276237963970698),
    }


def _spk_sup_power_law(kh, z, so=200, fb_a=0.4, fb_pow=0.2, fb_pivot=1e14):
    p = _spk_params(so, z)
    best_mass = p["a"] - (p["a"] - p["b"]) * (kh ** p["g"])
    m_opt = 10.0**best_mass
    fb = fb_a * ((m_opt / fb_pivot) ** fb_pow)

    x = np.log10(kh)
    x0 = 1.0 + p["lambda_a"] * np.exp(p["lambda_b"] * x)
    x1 = p["mu_a"] + ((1.0 - p["mu_a"]) / (1.0 + np.exp(p["mu_b"] * x + p["mu_c"])))
    x2 = p["nu_a"] * np.exp(-0.5 * ((x - p["nu_b"]) / p["nu_c"]) ** 2)
    sup = x0 - (x0 - x1) * np.exp(-x2 * fb)
    return np.maximum(sup, 1e-6)


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
        return k, pk[0]

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
        k_base, pk_base = self._get_pk(base)

        spk = camb.SPkNonLinear()
        spk.BaseModel.set_params(halofit_version="mead2020")
        spk.set_params(SPk_feedback=False)
        k_spk, pk_spk = self._get_pk(spk)

        self.assertTrue(np.allclose(k_base, k_spk, rtol=0, atol=0))
        self.assertTrue(np.allclose(pk_base, pk_spk, rtol=2e-12, atol=1e-14))

    def test_spk_power_law_matches_reference(self):
        fb_a = 0.4
        fb_pow = 0.2
        fb_pivot = 1e14
        z = 0.5

        base = camb.Halofit()
        base.set_params(halofit_version="mead2020")
        k, pk_base = self._get_pk(base, z=z)

        spk = camb.SPkNonLinear()
        spk.BaseModel.set_params(halofit_version="mead2020")
        spk.set_params(
            SPk_feedback=True,
            SPk_SO=200,
            SPk_relation_kind=1,
            SPk_fb_a=fb_a,
            SPk_fb_pow=fb_pow,
            SPk_fb_pivot=fb_pivot,
        )
        _, pk_spk = self._get_pk(spk, z=z)

        expected_sup = _spk_sup_power_law(k, z, so=200, fb_a=fb_a, fb_pow=fb_pow, fb_pivot=fb_pivot)
        measured_sup = pk_spk / pk_base
        self.assertTrue(np.allclose(measured_sup, expected_sup, rtol=1e-5, atol=1e-8))


if __name__ == "__main__":
    unittest.main()
