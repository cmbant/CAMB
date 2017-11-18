# Classes to produce BBN predictions for helium mass fraction and D/H given other parameters
# Fitting formula for Parthenelope via Julien Lesourges Dec 2014
# General interpolation tables, with defaults from Parthenope May 2017 (thanks Ofelia Pesanti)

import numpy as np
import os
import io
from .baseconfig import mock_load, needs_scipy

# Various useful constants
hbar = 1.05457e-34
c = 299792458.
kB = 1.380650e-23
MeV = 1.6021764e-13
eV = MeV / 1e6
G = 6.6738e-11
TCMB = 2.7255
mP = np.sqrt(hbar * c / G)
Mpc = 3.0856776e22

m_proton = 1.672621637e-27
m_H = 1.673575e-27
not4 = 3.9715
m_He = m_H * not4

zeta3 = 1.202056903

n_photon = (kB * TCMB / hbar / c) ** 3 * zeta3 * 2 / np.pi ** 2
omegafac = (1e5 / Mpc) ** 2 / (8 * np.pi * G) * 3


def yhe_to_ypBBN(Yp):
    return -4 * m_H * Yp / (Yp * m_He - 4 * Yp * m_H - m_He)


def ypBBN_to_yhe(YBBN):
    return -YBBN * m_He / (-YBBN * m_He + 4 * YBBN * m_H - 4 * m_H)


class BBNPredictor(object):
    def Y_p(self, ombh2, delta_neff=0.):
        """
        Get BBN helium nucleon fraction. Must be implemented by extensions.

        :param ombh2: Omega_b h^2
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return:  Y_p helium nucleon fraction predicted by BBN
        """
        raise Exception('Not implemented')

    def Y_He(self, ombh2, delta_neff=0.):
        """
        Get BBN helium mass fraction for CMB code.

        :param ombh2: Omega_b h^2
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return: Y_He helium mass fraction predicted by BBN
        """
        return ypBBN_to_yhe(self.Y_p(ombh2, delta_neff))


class BBN_table_interpolator(BBNPredictor):
    """
    BBN predictor based on interpolation on a table calculated from BBN code
    """

    def __init__(self, interpolation_table='PArthENoPE_880.2_standard.dat'):
        """
        Load table file and initialize interpolation

        :param interpolation_table: filename of interpolation table to use.
        """

        if not os.sep in interpolation_table:
            interpolation_table = os.path.join(os.path.dirname(__file__), interpolation_table)
        self.interpolation_table = interpolation_table

        comment = None
        with io.open(interpolation_table) as f:
            for line in f:
                line = line.strip()
                if line:
                    if line[0] == '#':
                        comment = line[1:]
                    else:
                        break
        assert comment
        columns = comment.split()
        ombh2_i = columns.index('ombh2')
        DeltaN_i = columns.index('DeltaN')
        Yp_i = columns.index('Yp^BBN')
        dh_i = columns.index('D/H')
        table = np.loadtxt(interpolation_table, usecols=[ombh2_i, DeltaN_i, Yp_i, dh_i])
        deltans = list(np.unique(table[:, 1]))
        ombh2s = list(np.unique(table[:, 0]))
        assert (table.shape[0] == len(ombh2s) * len(deltans))
        grid = np.zeros((len(ombh2s), len(deltans)))
        dh_grid = np.zeros(grid.shape)
        for i in range(table.shape[0]):
            grid[ombh2s.index(table[i, 0]), deltans.index(table[i, 1])] = table[i, 2]
            dh_grid[ombh2s.index(table[i, 0]), deltans.index(table[i, 1])] = table[i, 3]

        needs_scipy()
        from scipy.integrate import RectBivariateSpline
        
        self.interpolator_Yp = RectBivariateSpline(ombh2s, deltans, grid)
        self.interpolator_DH = RectBivariateSpline(ombh2s, deltans, dh_grid)

    def Y_p(self, ombh2, delta_neff=0.):
        """
        Get BBN helium nucleon fraction by intepolation in table.

        :param ombh2: Omega_b h^2
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return:  Y_p helium nucleon fraction predicted by BBN. Call Y_He() to get mass fraction instead.
        """
        return self.interpolator_Yp(ombh2, delta_neff)[0][0]

    def DH(self, ombh2, delta_neff=0.):
        """
        Get deuterium ration D/H by interpolation in table

        :param ombh2: Omega_b h^2
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return: D/H
        """
        return self.interpolator_DH(ombh2, delta_neff)[0][0]


class BBN_fitting_parthenope(BBNPredictor):
    """
    BBN predictions for Helium abundance using fitting formulae based on Parthenope (pre 2015)
    """

    def __init__(self, tau_neutron=None):
        """
        :param tau_neutron: fitting formula can use different neutron lifetime, defaults to 880.3s if not specified.
        """
        if tau_neutron is None:
            self.taun = 880.3  # Planck 2015 default
        else:
            self.taun = tau_neutron

    def Y_p(self, ombh2, delta_neff=0., tau_neutron=None):
        """
        Get BBN helium nucleon fraction.
        # Parthenope fits, as in Planck 2015 papers

        :param ombh2: Omega_b h^2
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :param tau_neutron: neutron lifetime
        :return:  Y_p helium nucleon fraction predicted by BBN
        """
        return (
                   0.2311 + 0.9502 * ombh2 - 11.27 * ombh2 * ombh2
                   + delta_neff * (0.01356 + 0.008581 * ombh2 - 0.1810 * ombh2 * ombh2)
                   + delta_neff * delta_neff * (-0.0009795 - 0.001370 * ombh2 + 0.01746 * ombh2 * ombh2)
               ) * pow((tau_neutron or self.taun) / 880.3, 0.728)

    def DH(self, ombh2, delta_neff, tau_neutron=None):
        return (
                   18.754 - 1534.4 * ombh2 + 48656. * ombh2 * ombh2 - 552670. * ombh2 ** 3
                   + delta_neff * (2.4914 - 208.11 * ombh2 + 6760.9 * ombh2 ** 2 - 78007. * ombh2 ** 3)
                   + delta_neff * delta_neff * (
                       0.012907 - 1.3653 * ombh2 + 37.388 * ombh2 ** 2 - 267.78 * ombh2 ** 3)
               ) * pow((tau_neutron or self.taun) / 880.3, 0.418) * 1e-5


_default_predictor = None


def get_default_predictor():
    """
    Get instance of default BBNPredictor class. Currently fitting formula to match Planck 2015 analysis.
    """
    global _default_predictor
    if _default_predictor is None:
        _default_predictor = BBN_fitting_parthenope()
    return _default_predictor


if __name__ == "__main__":
    print(BBN_table_interpolator().Y_He(0.01897, 1.2))
    print(BBN_fitting_parthenope().Y_He(0.01897, 1.2))
    print(BBN_table_interpolator().DH(0.02463, -0.6))
    print(BBN_fitting_parthenope().DH(0.02463, -0.6))
    print(BBN_table_interpolator('PArthENoPE_880.2_marcucci.dat').DH(0.02463, -0.6))
