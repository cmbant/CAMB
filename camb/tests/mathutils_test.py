import contextlib
import io
import os
import sys
import unittest

import numpy as np
from scipy.optimize import brentq
from scipy.special import airy, spherical_jn

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb  # noqa: F401

from camb.mathutils import (
    airy_ai_fast,
    airy_fast,
    chi_squared,
    pcl_coupling_matrix,
    phi_derivative,
    phi_first_peak_amplitude,
    phi_first_peak_chi,
    phi_olver,
    phi_recurs,
    scalar_coupling_matrix,
    threej_coupling,
)


def hyperspherical_turning_point(ell, curvature, nu):
    turn_arg = np.sqrt(ell * (ell + 1)) / nu
    if curvature == -1:
        return np.arcsinh(turn_arg)
    if curvature == 1:
        return np.arcsin(min(turn_arg, 1))
    return turn_arg


def hyperspherical_derivative(ell, curvature, nu, chi):
    if curvature == -1:
        cot_k = 1 / np.tanh(chi)
        root_l = np.sqrt(nu**2 + ell**2)
    elif curvature == 1:
        cot_k = 1 / np.tan(chi)
        root_l = np.sqrt(nu**2 - ell**2)
    else:
        cot_k = 1 / chi
        root_l = nu
    return root_l * phi_recurs(ell - 1, curvature, nu, chi) - (ell + 1) * cot_k * phi_recurs(ell, curvature, nu, chi)


def first_peak_after_turn(ell, curvature, nu):
    turn = hyperspherical_turning_point(ell, curvature, nu)
    upper = np.pi / 2 if curvature == 1 else np.inf
    turn = min(max(turn, 1e-8), upper)
    if hyperspherical_derivative(ell, curvature, nu, turn) <= 0:
        return turn  # no rise after the turn (e.g. open models with nu << 1)

    if curvature == -1:
        cot_t = 1 / np.tanh(turn)
    elif curvature == 1:
        cot_t = 1 / np.tan(turn)
    else:
        cot_t = 1 / turn
    # Near the turn k^2 ~ 2 nu^2 cot_K(turn) (chi - turn), so the first peak is
    # close to the first maximum of Ai(-x) at x = 1.0188
    delta = 1.0188 / max(2 * nu**2 * cot_t, 1e-30) ** (1 / 3)
    hi = min(turn + 2 * delta, upper)
    while hyperspherical_derivative(ell, curvature, nu, hi) > 0:
        if hi >= upper:
            return upper  # search boundary reached before bracketing a stationary peak
        hi = min(turn + 2 * (hi - turn), upper)
    return brentq(lambda chi: hyperspherical_derivative(ell, curvature, nu, chi), turn, hi, xtol=1e-7)


def open_alpha_cut(ell):
    ell_float = float(ell)
    if ell_float < 500:
        alpha_cut = 0.12 * (500 / ell_float) ** 0.70
    else:
        alpha_cut = 0.12 * (500 / ell_float) ** 0.14
    return max(0.095, alpha_cut)


def open_smallnu_gate(ell):
    ell_float = float(ell)
    old_gate = max(1.2e-4 * ell_float**1.5, min(12 * (ell_float / 1000) ** 3, 0.032 * ell_float))
    if ell >= 3000:
        return max(old_gate, 0.04 * ell_float, min(8 * (ell_float / 1000) ** 2, 0.16 * ell_float))
    return old_gate


class MathutilsTest(unittest.TestCase):
    def assert_hyperspherical_peak_relative(self, ell, curvature, nu, chi_grid, max_peak_error=2e-4):
        recurs_grid = phi_recurs(ell, curvature, nu, chi_grid)
        olver_grid = phi_olver(ell, curvature, nu, chi_grid)
        peak_norm = phi_first_peak_amplitude(ell, curvature, nu)
        self.assertGreater(peak_norm, 0)
        np.testing.assert_allclose(
            (olver_grid - recurs_grid) / peak_norm,
            0,
            atol=max_peak_error,
            err_msg=f"peak-relative hyperspherical Bessel error for L={ell}, K={curvature}, nu={nu:g}",
        )

    def test_mathutils(self):
        airy_x = np.array([-20.0, -6.5, -4.4, -2.09, -1.0, 0.0, 2.09, 2.98, 6.5, 20.0, 30.0])
        ai_expected, aip_expected, _, _ = airy(airy_x)
        ai, aip = airy_fast(airy_x)
        np.testing.assert_allclose(ai, ai_expected, atol=8e-8, rtol=0)
        np.testing.assert_allclose(aip, aip_expected, atol=8e-8, rtol=0)
        np.testing.assert_allclose(airy_ai_fast(airy_x), ai, atol=0, rtol=1e-15)

        ai_scalar, aip_scalar = airy_fast(0.0)
        np.testing.assert_allclose(ai_scalar, airy(0.0)[0], atol=8e-8, rtol=0)
        np.testing.assert_allclose(aip_scalar, airy(0.0)[1], atol=8e-8, rtol=0)
        self.assertAlmostEqual(airy_ai_fast(0.0), ai_scalar)

        cinv = np.linalg.inv(np.array([[1.2, 3], [3, 18.2]]))
        vec = np.array([0.5, 5.0])
        self.assertAlmostEqual(chi_squared(cinv, vec), cinv.dot(vec).dot(vec))

        for function in (phi_recurs, phi_olver):
            with self.assertRaisesRegex(ValueError, "chi must be non-negative"):
                function(60, -1, 100.0, -0.1)
            with self.assertRaisesRegex(ValueError, "chi must be non-negative"):
                function(60, -1, 100.0, np.array([0.1, -0.1]))

        self.assertEqual(phi_recurs(2, 0, 0.0, 1.0), 0.0)
        open_zero_chi = 1.3
        open_zero_phi0 = open_zero_chi / np.sinh(open_zero_chi)
        open_zero_phi1 = (open_zero_chi / np.tanh(open_zero_chi) - 1) / np.sinh(open_zero_chi)
        open_zero_phi2 = (3 / np.tanh(open_zero_chi) * open_zero_phi1 - open_zero_phi0) / 2
        np.testing.assert_allclose(phi_recurs(0, -1, 0.0, open_zero_chi), open_zero_phi0, rtol=1e-12)
        np.testing.assert_allclose(phi_recurs(1, -1, 0.0, open_zero_chi), open_zero_phi1, rtol=1e-12)
        np.testing.assert_allclose(phi_recurs(2, -1, 0.0, open_zero_chi), open_zero_phi2, rtol=1e-12)
        self.assertEqual(phi_recurs(10, -1, 30.0, 1e-300), 0.0)
        self.assertEqual(phi_recurs(10, 0, 30.0, 1e-300), 0.0)
        self.assertEqual(phi_recurs(10, 1, 30.0, 1e-300), 0.0)
        np.testing.assert_allclose(
            phi_recurs(1, -1, 30.0, 1e-300),
            1e-300 * np.sqrt(30.0**2 + 1) / 3,
            rtol=1e-12,
            atol=0,
        )
        chis = np.array([0.0, 0.1, 1.0, 10.0])
        np.testing.assert_allclose(phi_recurs(25, 0, 3.5, chis), spherical_jn(25, 3.5 * chis), rtol=1e-12)
        np.testing.assert_allclose(phi_olver(25, 0, 3.5, chis), spherical_jn(25, 3.5 * chis), rtol=1e-12, atol=1e-24)
        for ell in (100, 1000):
            np.testing.assert_allclose(phi_recurs(ell, 0, 1.0, ell), spherical_jn(ell, ell), rtol=1e-12)
            np.testing.assert_allclose(phi_olver(ell, 0, 1.0, ell), spherical_jn(ell, ell), rtol=1e-9)
        for arg in (0.999e-4, 1.0e-4, 1.001e-4, 1e-3):
            expected = arg * (1 / 3 - arg**2 / 30 + arg**4 / 840)
            np.testing.assert_allclose(phi_recurs(1, 0, 2.0, arg / 2.0), expected, rtol=1e-13)

        small_chi = 5e-5
        for curvature in (-1, 1):
            nu = 200000.0
            sin_K = np.sinh(small_chi) if curvature == -1 else np.sin(small_chi)
            cot_K = 1 / np.tanh(small_chi) if curvature == -1 else 1 / np.tan(small_chi)
            root1 = np.sqrt(nu**2 - curvature)
            arg = nu * small_chi
            phi0 = np.sin(arg) / (nu * sin_K)
            phi1 = (np.sin(arg) * cot_K / (nu * sin_K) - np.cos(arg) / sin_K) / root1
            np.testing.assert_allclose(phi_recurs(0, curvature, nu, small_chi), phi0, rtol=1e-12)
            np.testing.assert_allclose(phi_recurs(1, curvature, nu, small_chi), phi1, rtol=1e-12)
        for curvature in (-1, 1):
            for arg in (0.999e-4, 1.0e-4, 1.001e-4):
                nu = 3.0 if curvature == 1 else 2.0
                small_chi = arg / nu
                root1 = np.sqrt(nu**2 - curvature)
                if abs(arg) < 1e-4:
                    chi2 = small_chi**2
                    phi1 = small_chi * root1 / 3 * (1 - (3 * nu**2 - 7 * curvature) * chi2 / 30)
                else:
                    chi2 = small_chi**2
                    arg2 = arg**2
                    arg4 = arg2**2
                    sinc = np.sin(arg) / arg
                    sinc_minus_cos = arg2 / 3 - arg4 / 30 + arg4 * arg2 / 840
                    chi_over_sin = 1 + curvature * chi2 / 6 + 7 * chi2**2 / 360
                    chi_cot_m1 = -curvature * chi2 / 3 - chi2**2 / 45
                    phi1 = (sinc_minus_cos + sinc * chi_cot_m1) * chi_over_sin / (root1 * small_chi)
                np.testing.assert_allclose(phi_recurs(1, curvature, nu, small_chi), phi1, rtol=1e-10)

        open_turning_chi = np.arctanh(2.0 / np.sqrt(2.0**2 + (2e-6) ** 2)) * (1 - 1e-4)
        self.assertAlmostEqual(phi_recurs(2, -1, 2e-6, open_turning_chi), 1.3026626711804369e-05)
        for args, expected in [
            ((1000, -1, 1500.0, 1.0), -5.798121096319023e-04),
            ((1000, 1, 1100.0, 1.2), -1.501736026757476e-03),
            ((1000, 1, 10000.0, 0.1), 1.427031695544760e-03),
            ((500, 1, 600.0, 1.0), 5.235950059980328e-03),
        ]:
            self.assertAlmostEqual(phi_recurs(*args), expected)
        for args, expected in [
            ((10, -1, 2.0, np.pi / 2), 1.241277248560563e-02),
            ((80, -1, 120.0, 0.55), 2.6546410483499028e-04),
            ((80, 1, 90.0, 1.0), 1.1532376890679929e-03),
        ]:
            np.testing.assert_allclose(phi_recurs(*args), expected, rtol=1e-11, atol=0)
        for args, expected in [
            ((4, -1, 0.008, 1.2191174089684585e-4), 5.610238786031219e-18),
            ((460, -1, 0.92, 6.82954817444722), 1.1752384591191022e-03),
            ((8000, -1, 16.0, 6.01704300224), 3.4785032603148615e-12),
            ((8000, -1, 16.0, np.arcsinh(0.999 * 8000 / 16)), 1.1060344911733589e-04),
            ((8000, -1, 16.0, 7.0), 1.496457080852643e-04),
            ((10000, -1, 20.0, 6.04202992999), 1.1204823640065975e-13),
            ((10000, -1, 100000.0, np.arcsinh(0.995 / 10)), 5.319664845847832e-06),
        ]:
            np.testing.assert_allclose(phi_recurs(*args), expected, rtol=1e-8, atol=0)

        for ell, curvature, nu in [
            (80, -1, 120.0),
            (500, -1, 25.0),
            (80, 1, 90.0),
            (500, 1, 700.0),
            (2000, 1, 2500.0),
        ]:
            peak_chi = first_peak_after_turn(ell, curvature, nu)
            fortran_peak_chi, no_peak_found = phi_first_peak_chi(ell, curvature, nu, return_status=True)
            self.assertFalse(no_peak_found)
            np.testing.assert_allclose(fortran_peak_chi, peak_chi, rtol=5e-7, atol=5e-8)
            np.testing.assert_allclose(
                phi_first_peak_amplitude(ell, curvature, nu),
                abs(phi_recurs(ell, curvature, nu, peak_chi)),
                rtol=5e-7,
                atol=0,
            )
            np.testing.assert_allclose(phi_derivative(ell, curvature, nu, fortran_peak_chi), 0, atol=1e-7)

        boundary_chi, no_peak_found = phi_first_peak_chi(80, 1, 81.0, return_status=True)
        self.assertTrue(no_peak_found)
        np.testing.assert_allclose(boundary_chi, np.pi / 2, rtol=0, atol=1e-12)

        for ell, curvature, nu, chi in [
            (80, -1, 120.0, 0.7),
            (500, -1, 25.0, 6.0),
            (80, 1, 90.0, 1.0),
            (500, 1, 700.0, 0.2),
        ]:
            np.testing.assert_allclose(
                phi_derivative(ell, curvature, nu, chi),
                hyperspherical_derivative(ell, curvature, nu, chi),
                rtol=1e-10,
                atol=1e-14,
            )

        for function in (phi_recurs, phi_olver):
            ell, curvature, nu, chi = 60, 1, 180.0, 0.4
            inu = round(nu)
            base = function(ell, curvature, nu, chi)
            for shifted_chi, sign in [
                (np.pi - chi, (-1) ** (inu - ell - 1)),
                (2 * np.pi - chi, (-1) ** ell),
                (np.pi + chi, (-1) ** (inu - 1)),
            ]:
                np.testing.assert_allclose(function(ell, curvature, nu, shifted_chi), sign * base, rtol=1e-12)

        for args in [
            (80, -1, 20.0, 6.0),
            (80, 1, 90.0, 1.0),
        ]:
            np.testing.assert_allclose(phi_olver(*args), phi_recurs(*args), rtol=1e-14, atol=0)

        nonflat_chis = np.geomspace(1e-5, np.pi, 192)
        for ell, curvature, nu in [
            (80, -1, 260.0),
            (80, 1, 260.0),
        ]:
            recurs_grid = phi_recurs(ell, curvature, nu, nonflat_chis)
            olver_grid = phi_olver(ell, curvature, nu, nonflat_chis)

            peak_norm = np.max(np.abs(recurs_grid))
            self.assertGreater(peak_norm, 0)
            np.testing.assert_allclose((olver_grid - recurs_grid) / peak_norm, 0, atol=1e-4)

        for args in [
            (100, -1, 130.0, 1.0),
            (500, -1, 700.0, 1.0),
            (100, 1, 500.0, 0.5),
            (500, 1, 700.0, 0.7),
        ]:
            np.testing.assert_allclose(phi_olver(*args), phi_recurs(*args), rtol=1e-4, atol=0)

        # Exercise both sides of the pointwise Olver, small-chi, Airy, and
        # small-nu fallback gates on log-spaced chi grids.
        for ell in (80, 500, 3000):
            for scale in (0.99, 1.01):
                nu = scale * open_alpha_cut(ell) * ell
                chi_grid = np.geomspace(1e-5, max(12, np.arcsinh(ell / nu) + 4), 48)
                self.assert_hyperspherical_peak_relative(ell, -1, nu, chi_grid)

        for ell in (80, 500, 3000):
            nu = 0.75 * open_alpha_cut(ell) * ell
            open_metric_chi = 2 * nu * 3e-3
            chi_grid = open_metric_chi * np.geomspace(0.1, 10, 48)
            self.assert_hyperspherical_peak_relative(ell, -1, nu, chi_grid)

        for ell in (100, 500, 3000):
            for scale in (0.99, 1.01):
                nu = scale * open_smallnu_gate(ell)
                chi_grid = np.geomspace(1e-6, max(12, np.arcsinh(ell / nu) + 4), 48)
                self.assert_hyperspherical_peak_relative(ell, -1, nu, chi_grid)

        for ell in (80, 500, 3000):
            nu = 3.0 * ell
            smallchi_edge = (5e-2 * nu / ell**2) ** (1 / 7)
            chi_grid = smallchi_edge * np.geomspace(0.1, 10, 48)
            for curvature in (-1, 1):
                self.assert_hyperspherical_peak_relative(ell, curvature, nu, chi_grid)

        for ell in (150, 500, 3000):
            for delta in (-1, 1):
                nu = ell + 100.0 + delta
                closed_metric_chi = 2 * (nu - ell) * 7e-3
                chi_grid = np.minimum(closed_metric_chi * np.geomspace(0.1, 10, 48), np.pi / 2)
                self.assert_hyperspherical_peak_relative(ell, 1, nu, chi_grid)

        for ell in (100, 500):
            for nu in (6.9, 7.1):
                chi_grid = np.geomspace(1e-4, max(10, np.arcsinh(ell / nu) + 4), 48)
                self.assert_hyperspherical_peak_relative(ell, -1, nu, chi_grid)
        for ell in (150, 500):
            lambda_ell = ell + 0.5
            boundary_nu = lambda_ell * np.sqrt(1 + 15 / lambda_ell)
            for nu in (np.floor(boundary_nu), np.ceil(boundary_nu)):
                chi_grid = np.geomspace(1e-4, np.pi / 2, 48)
                self.assert_hyperspherical_peak_relative(ell, 1, float(nu), chi_grid)

        W = np.zeros(100)
        W[0] = 1
        lmax = len(W)
        Xi = threej_coupling(W, lmax)
        np.testing.assert_allclose(np.diag(Xi) * (2 * np.arange(lmax + 1) + 1), np.ones(lmax + 1))
        Xis = threej_coupling(W, lmax, pol=True)
        np.testing.assert_allclose(np.diag(Xis[0]) * (2 * np.arange(lmax + 1) + 1), np.ones(lmax + 1))
        P = W * 4 * np.pi
        # mathutils prints a Python warning when the mask power lmax is < 2*lmax, which is expected here
        with contextlib.redirect_stdout(io.StringIO()):
            M = scalar_coupling_matrix(P, lmax)
            np.testing.assert_allclose(M, np.eye(lmax + 1))
            M = pcl_coupling_matrix(P, lmax)
            np.testing.assert_allclose(M, np.eye(lmax + 1))
