module Constants

export const_pi, const_twopi, const_fourpi, const_eightpi,
       const_sqrt6, zeta3, zeta5, zeta7,
       c, h_P, hbar, G, sigma_thomson, sigma_boltz, k_B, eV,
       m_p, m_e, m_H, m_He4, mass_ratio_He_H,
       Gyr, Mpc, MPC_in_sec,
       barssc0, kappa, a_rad, Compton_CT,
       f_21cm, l_21cm, T_21cm, A10, B10, line21_const,
       COBE_CMBTemp, inv_neutrino_mass_fac, neutrino_mass_fac,
       default_nnu, delta_mnu21, delta_mnu31, mnu_min_normal,
       Ini_max_string_len

const_pi = 3.1415926535897932384626433832795
const_twopi = 2 * const_pi
const_fourpi = 4 * const_pi
const_eightpi = 8 * const_pi
const_sqrt6 = 2.4494897427831780981972840747059
zeta3 = 1.2020569031595942853997
zeta5 = 1.0369277551433699263313
zeta7 = 1.0083492773819228268397

c = 2.99792458e8
h_P = 6.62607015e-34
hbar = h_P / const_twopi

G = 6.67430e-11
sigma_thomson = 6.6524587321e-29
sigma_boltz = 5.670374419e-8
k_B = 1.380649e-23
eV = 1.602176634e-19

m_p = 1.67262192369e-27
m_e = 9.1093837015e-31
m_H = 1.673575e-27
m_He4 = 6.646479073e-27
mass_ratio_He_H = m_He4 / m_H

Gyr = 365.25e9 * 86400
Mpc = 3.085677581e22
MPC_in_sec = Mpc / c

barssc0 = k_B / m_p / c^2
kappa = 8 * const_pi * G
a_rad = 8 * const_pi^5 * k_B^4 / 15 / c^3 / h_P^3
Compton_CT = MPC_in_sec * (8 / 3) * (sigma_thomson / (m_e * c)) * a_rad

f_21cm = 1420.40575e6
l_21cm = c / f_21cm
T_21cm = h_P * f_21cm / k_B
A10 = 2.869e-15
B10 = l_21cm^3 / 2 / h_P / c * A10
line21_const = 3 * l_21cm^2 * c * h_P / 32 / const_pi / k_B * A10 * MPC_in_sec * 1000

COBE_CMBTemp = 2.7255
inv_neutrino_mass_fac = zeta3 * 3 / 2 / const_pi^2 * 4 / 11 * ((k_B * COBE_CMBTemp / hbar / c)^3 * kappa / 3 / (100 * 1000 / Mpc)^2 / (c^2 / eV))
neutrino_mass_fac = 1 / inv_neutrino_mass_fac

default_nnu = 3.044
delta_mnu21 = 7.54e-5
delta_mnu31 = 2.46e-3
mnu_min_normal = 0.06

Ini_max_string_len = 1024

end # module
