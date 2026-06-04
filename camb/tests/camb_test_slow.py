import os
import sys
import unittest

import numpy as np

try:
    import camb
    from camb import check_accuracy
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb
    from camb import check_accuracy


def _accuracy_failure_message(comparison):
    rows = check_accuracy.comparison_rows(comparison)
    failures = [row for row in rows if not row.passed]
    if not failures:
        return "accuracy comparison passed"
    return "\n".join(
        f"{row.quantity} {row.range_label}: max={row.max_abs:.4g}, tolerance={row.tolerance:.4g}, L={row.location}"
        for row in failures
    )


def _standard_model(omk=0, lmax=None, lens_potential_accuracy=None, matter_power=False):
    params = camb.set_params(
        H0=67.5,
        ombh2=0.022,
        omch2=0.122,
        omk=omk,
        mnu=0.06,
        tau=0.054,
        As=2.1e-9,
        ns=0.965,
    )
    params.set_dark_energy(w=-0.95, wa=0.15, dark_energy_model="ppf")
    if matter_power:
        params.set_matter_power(redshifts=[0.0, 0.5, 1.0], kmax=2.0, silent=True)
    if lmax is not None:
        lens_accuracy = 0 if lens_potential_accuracy is None else lens_potential_accuracy
        params.set_for_lmax(lmax, lens_potential_accuracy=lens_accuracy)
    return params


class CambSlowTest(unittest.TestCase):
    def test_high_l_lensing_accuracy_and_near_flat_continuity(self):
        lmax = 4000
        lens_potential_accuracy = 8

        accuracy_result = check_accuracy.compare_params_accuracy(
            _standard_model(matter_power=True),
            lmax=lmax,
            set_for_lmax=lmax,
            lens_potential_accuracy=lens_potential_accuracy,
        )
        self.assertTrue(
            accuracy_result.comparison.passed,
            _accuracy_failure_message(accuracy_result.comparison),
        )

        flat_cls = camb.get_results(
            _standard_model(lmax=lmax, lens_potential_accuracy=lens_potential_accuracy)
        ).get_total_cls(
            lmax,
        )
        ell = np.arange(lmax + 1)
        min_ell = 2
        spectra_slice = slice(min_ell, lmax + 1)
        high_tt_slice = slice(20, lmax + 1)

        for omk in (-1e-6, 1e-6):
            near_flat_cls = camb.get_results(
                _standard_model(omk=omk, lmax=lmax, lens_potential_accuracy=lens_potential_accuracy)
            ).get_total_cls(lmax)

            with np.errstate(divide="ignore", invalid="ignore"):
                low_tt_delta = (near_flat_cls[min_ell:20, 0] - flat_cls[min_ell:20, 0]) / flat_cls[min_ell:20, 0]
            self.assertLess(np.max(np.abs(low_tt_delta)), 1.1e-3)
            np.testing.assert_allclose(
                near_flat_cls[high_tt_slice, 0],
                flat_cls[high_tt_slice, 0],
                rtol=1e-3,
            )
            np.testing.assert_allclose(
                near_flat_cls[spectra_slice, 1],
                flat_cls[spectra_slice, 1],
                rtol=1e-3,
            )
            np.testing.assert_allclose(
                near_flat_cls[spectra_slice, 2],
                flat_cls[spectra_slice, 2],
                rtol=1e-3,
            )
            te_delta = check_accuracy.normalized_te_delta(near_flat_cls, flat_cls)
            self.assertLess(np.max(np.abs(te_delta[ell >= min_ell])), 1e-3)

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
