import os
import sys
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb


class CambSlowTest(unittest.TestCase):
    def testSymbolic(self):
        import camb.symbolic as s

        monopole_source, ISW, doppler, quadrupole_source = s.get_scalar_temperature_sources()
        temp_source = monopole_source + ISW + doppler + quadrupole_source

        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, omk=0.1)
        data = camb.get_background(pars)
        tau = np.linspace(1, 1200, 300)
        ks = [0.001, 0.05, 1]
        monopole2 = s.make_frame_invariant(s.newtonian_gauge(monopole_source), "Newtonian")
        Delta_c_N = s.make_frame_invariant(s.Delta_c, "Newtonian")
        Delta_c_N2 = s.make_frame_invariant(s.synchronous_gauge(Delta_c_N), "CDM")
        ev = data.get_time_evolution(
            ks,
            tau,
            ["delta_photon", s.Delta_g, Delta_c_N, Delta_c_N2, monopole_source, monopole2, temp_source, "T_source"],
        )
        self.assertTrue(np.allclose(ev[:, :, 0], ev[:, :, 1]))
        self.assertTrue(np.allclose(ev[:, :, 2], ev[:, :, 3]))
        self.assertTrue(np.allclose(ev[:, :, 4], ev[:, :, 5]))
        self.assertTrue(np.allclose(ev[:, :, 6], ev[:, :, 7]))

        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95)
        pars.set_accuracy(lSampleBoost=2)
        try:
            pars.set_custom_scalar_sources(
                [monopole_source + ISW + doppler + quadrupole_source, s.scalar_E_source],
                source_names=["T2", "E2"],
                source_ell_scales={"E2": 2},
            )
            data = camb.get_results(pars)
            dic = data.get_cmb_unlensed_scalar_array_dict(CMB_unit="muK")
            self.assertTrue(np.all(np.abs(dic["T2xT2"][2:2000] / dic["TxT"][2:2000] - 1) < 1e-3))
            self.assertTrue(np.all(np.abs(dic["TxT2"][2:2000] / dic["TxT"][2:2000] - 1) < 1e-3))
            # default interpolation errors much worse for E
            self.assertTrue(np.all(np.abs(dic["E2xE2"][10:2000] / dic["ExE"][10:2000] - 1) < 2e-3))
            self.assertTrue(np.all(np.abs(dic["E2xE"][10:2000] / dic["ExE"][10:2000] - 1) < 2e-3))
            dic1 = data.get_cmb_power_spectra(CMB_unit="muK")
            self.assertTrue(np.allclose(dic1["unlensed_scalar"][2:2000, 1], dic["ExE"][2:2000]))
        finally:
            pars.set_accuracy(lSampleBoost=1)

        s.internal_consistency_checks()

    def test_extra_EmissionAnglePostBorn(self):
        from camb import emission_angle, postborn

        pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95, tau=0.055)
        BB = emission_angle.get_emission_delay_BB(pars, lmax=3500)
        self.assertAlmostEqual(BB(80) * 2 * np.pi / 80 / 81.0, 1.1e-10, delta=1e-11)  # type: ignore

        Bom = postborn.get_field_rotation_BB(pars, lmax=3500)
        self.assertAlmostEqual(Bom(100) * 2 * np.pi / 100 / 101.0, 1.65e-11, delta=1e-12)  # type: ignore
