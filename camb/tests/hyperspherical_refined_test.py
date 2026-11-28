"""Independent reference values for the refined curved-space Olver fallback."""
import unittest

import numpy as np

from camb.mathutils import phi_olver


class HypersphericalRefinedTest(unittest.TestCase):
    def test_exact_representations(self):
        # 70-digit direct Mehler--Dirichlet integrals (open) and normalized
        # Gegenbauer polynomials (closed), at the stated double-precision inputs.
        # See docs/changelog/hyperspherical_refined.tex.
        cases = [
            (21, -1, 0.002775, 4.33147, 2.16801652908675084e-2),
            (32, -1, 0.0001979, 3.24647, 4.68404857833203386e-3),
            (10000, -1, 0.0, 12.0, 2.73360386187274689e-5),
            (1000, -1, 0.0, 6.0, 1.91190191127807087e-5),
            # Near the origin the outgoing-wave formula is ill-conditioned.
            (3, -1, 0.001, 1e-5, 5.71428960295216958e-17),
            (10000, 1, 10001.0, 1.565, 7.95767812981843480e-4),
            (10000, 1, 10002.0, 1.570795, 1.76628752881692199e-7),
            (10000, 1, 10006.0, 1.55, 9.71978847230654905e-5),
            (80, -1, 20.0, 6.0, -1.13135616748027954e-4),
        ]
        for ell, k, nu, chi, expected in cases:
            with self.subTest(ell=ell, k=k, nu=nu, chi=chi):
                np.testing.assert_allclose(phi_olver(ell, k, nu, chi), expected, rtol=5e-10, atol=0)
        # Closed odd parity about the equator, including the array wrapper.
        chi = 1.570795
        expected = 1.76628752881692199e-7
        np.testing.assert_allclose(
            phi_olver(10000, 1, 10002.0, np.array([chi, np.pi - chi])),
            [expected, -expected],
            rtol=2e-8,
            atol=1e-15,
        )

    def test_corrected_airy(self):
        # Both curvature signs and both sides of the effective-mu=4,16 cuts.
        cases = [
            (1000, -1, 6.9, 5.67, 7.72643215530129058e-4),
            (1000, -1, 20.0, 4.6055, 9.21450116711150708e-4),
            (500, 1, 520.0, 1.28, 1.74325858674533800e-3),
            (1000, 1, 1050.0, 1.3, -2.55823146406886125e-3),
            (1000, -1, 4.000031568404025, 6.215104077436099, 7.04711012643264917e-4),
            (1000, -1, 4.000032368410419, 6.215103877437697, 7.04711036447474497e-4),
            (1000, -1, 16.002044745844547, 4.828749770150602, 8.89204301945456261e-4),
            (1000, -1, 16.002047946253818, 4.8287495701761785, 8.89204331649054566e-4),
        ]
        for ell, k, nu, chi, expected in cases:
            with self.subTest(ell=ell, k=k, nu=nu, chi=chi):
                np.testing.assert_allclose(phi_olver(ell, k, nu, chi), expected, rtol=4e-6, atol=0)
