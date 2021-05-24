# Classes to produce BBN predictions for helium mass fraction and D/H given other parameters
# Fitting formula for Parthenelope via Julien Lesgourgues Dec 2014
# General interpolation tables, with defaults from Parthenope May 2017 (thanks Ofelia Pesanti)
# Use PRIMAT_Yp_DH_Error.dat table for latest from the PRIMAT code (arXiv: 1801.08023, thanks Cyril Pitrou)

import numpy as np
import os
from scipy.interpolate import RectBivariateSpline

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

default_interpolation_table = 'PArthENoPE_880.2_standard.dat'


def yhe_to_ypBBN(Yp):
    return -4 * m_H * Yp / (Yp * m_He - 4 * Yp * m_H - m_He)


def ypBBN_to_yhe(YBBN):
    return -YBBN * m_He / (-YBBN * m_He + 4 * YBBN * m_H - 4 * m_H)


class BBNIterpolator(RectBivariateSpline):
    grid: np.ndarray


class BBNPredictor:
    """
    The base class for making BBN predictions for Helium abundance
    """

    def Y_p(self, ombh2, delta_neff=0.):
        r"""
        Get BBN helium nucleon fraction. Must be implemented by extensions.

        :param ombh2: :math:`\Omega_b h^2`
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return:  Y_p helium nucleon fraction predicted by BBN
        """
        raise Exception('Not implemented')

    def Y_He(self, ombh2, delta_neff=0.):
        r"""
        Get BBN helium mass fraction for CMB code.

        :param ombh2: :math:`\Omega_b h^2`
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :return: Y_He helium mass fraction predicted by BBN
        """
        return ypBBN_to_yhe(self.Y_p(ombh2, delta_neff))


class BBN_table_interpolator(BBNPredictor):
    """
    BBN predictor based on interpolation from a numerical table calculated by a BBN code.

    Tables are supplied for `Parthenope <http://parthenope.na.infn.it/>`_ 2017 (PArthENoPE_880.2_standard.dat, default),
    similar but with Marucci rates (PArthENoPE_880.2_marcucci.dat),
    `PRIMAT <http://www2.iap.fr/users/pitrou/primat.htm>`_ (PRIMAT_Yp_DH_Error.dat, PRIMAT_Yp_DH_ErrorMC_2021.dat).

    :param interpolation_table: filename of interpolation table to use.
    :param function_of: two variables that determine the interpolation grid (x,y) in the table,
         matching top column label comment. By default ombh2, DeltaN, and function argument names reflect that,
         but can also be used more generally.

    """

    def __init__(self, interpolation_table=default_interpolation_table, function_of=('ombh2', 'DeltaN')):

        if os.sep not in interpolation_table and '/' not in interpolation_table:
            interpolation_table = os.path.normpath(os.path.join(os.path.dirname(__file__), interpolation_table))
        self.interpolation_table = interpolation_table

        comment = None
        with open(interpolation_table) as f:
            for line in f:
                line = line.strip()
                if line:
                    if line[0] == '#':
                        comment = line[1:]
                    else:
                        break
        assert comment
        columns = comment.split()
        ombh2_i = columns.index(function_of[0])
        DeltaN_i = columns.index(function_of[1])

        table = np.loadtxt(interpolation_table)
        deltans = list(np.unique(table[:, DeltaN_i]))
        ombh2s = list(np.unique(table[:, ombh2_i]))
        assert (table.shape[0] == len(ombh2s) * len(deltans))
        self.interpolators = {}
        for i, col in enumerate(columns):
            if i != ombh2_i and i != DeltaN_i and np.count_nonzero(table[:, i]):
                grid = np.zeros((len(ombh2s), len(deltans)))
                for ix in range(table.shape[0]):
                    grid[ombh2s.index(table[ix, ombh2_i]), deltans.index(table[ix, DeltaN_i])] = table[ix, i]
                self.interpolators[col] = BBNIterpolator(ombh2s, deltans, grid)
                self.interpolators[col].grid = grid

        self.ombh2s = ombh2s
        self.deltans = deltans

    def Y_p(self, ombh2, delta_neff=0., grid=False):
        r"""
        Get BBN helium nucleon fraction by intepolation in table.

        :param ombh2: :math:`\Omega_b h^2` (or, more generally, value of function_of[0])
        :param delta_neff:  additional N_eff relative to standard value (of 3.046) (or value of function_of[1])
        :param grid: parameter for :class:`~scipy:scipy.interpolate.RectBivariateSpline` (whether to evaluate the
           results on a grid spanned by the input arrays, or at points specified by the input arrays)
        :return:  Y_p helium nucleon fraction predicted by BBN. Call Y_He() to get mass fraction instead.
        """
        return self.get('Yp^BBN', ombh2, delta_neff, grid)

    def DH(self, ombh2, delta_neff=0., grid=False):
        r"""
        Get deuterium ratio D/H by interpolation in table

        :param ombh2: :math:`\Omega_b h^2` (or, more generally, value of function_of[0])
        :param delta_neff:  additional N_eff relative to standard value (of 3.046) (or value of function_of[1])
        :param grid: parameter for :class:`~scipy:scipy.interpolate.RectBivariateSpline` (whether to evaluate the
           results on a grid spanned by the input arrays, or at points specified by the input arrays)
        :return: D/H
        """
        return self.get('D/H', ombh2, delta_neff, grid)

    def get(self, name, ombh2, delta_neff=0., grid=False):
        r"""
        Get value for variable "name" by intepolation from table (where name is given in the column header comment)
        For example get('sig(D/H)',0.0222,0) to get the error on D/H

        :param name: string name of the parameter, as given in header of interpolation table
        :param ombh2: :math:`\Omega_b h^2` (or, more generally, value of function_of[0])
        :param delta_neff:  additional N_eff relative to standard value (of 3.046) (or value of function_of[1])
        :param grid: parameter for :class:`~scipy:scipy.interpolate.RectBivariateSpline` (whether to evaluate the
           results on a grid spanned by the input arrays, or at points specified by the input arrays)
        :return:  Interpolated value (or grid)
        """
        if name not in self.interpolators:
            raise ValueError('Unknown BBN table column index "%s"' % name)
        res = self.interpolators[name](ombh2, delta_neff, grid=grid)
        if np.isscalar(ombh2) and np.isscalar(delta_neff):
            return np.float64(res)
        return res


class BBN_fitting_parthenope(BBNPredictor):
    """
    Old BBN predictions for Helium abundance using fitting formulae based on Parthenope (pre 2015).
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
        r"""
        Get BBN helium nucleon fraction.
        # Parthenope fits, as in Planck 2015 papers

        :param ombh2: :math:`\Omega_b h^2`
        :param delta_neff:  additional N_eff relative to standard value (of 3.046)
        :param tau_neutron: neutron lifetime
        :return:  :math:`Y_p^{\rm BBN}` helium nucleon fraction predicted by BBN
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


_predictors = {}


def get_predictor(predictor_name=None):
    """
    Get instance of default BBNPredictor class. Currently numerical table interpolation as Planck 2018 analysis.
    """
    global _predictors
    predictor_name = predictor_name or default_interpolation_table
    predictor = _predictors.get(predictor_name, None)
    if predictor is None:
        if predictor_name == 'BBN_fitting_parthenope':
            predictor = BBN_fitting_parthenope()
        else:
            predictor = BBN_table_interpolator(interpolation_table=predictor_name)
        _predictors[predictor_name] = predictor
    return predictor


if __name__ == "__main__":
    print(BBN_table_interpolator().Y_He(0.01897, 1.2))
    print(BBN_fitting_parthenope().Y_He(0.01897, 1.2))
    print(BBN_table_interpolator('PRIMAT_Yp_DH_Error.dat').Y_He(0.01897, 1.2))
    print(BBN_table_interpolator().DH(0.02463, -0.6))
    print(BBN_fitting_parthenope().DH(0.02463, -0.6))
    print(BBN_table_interpolator('PArthENoPE_880.2_marcucci.dat').DH(0.02463, -0.6))
    print(BBN_table_interpolator('PRIMAT_Yp_DH_Error.dat').DH(0.02463, -0.6))
    print(BBN_table_interpolator('PRIMAT_Yp_DH_ErrorMC_2021.dat').DH(0.02463, -0.6))
