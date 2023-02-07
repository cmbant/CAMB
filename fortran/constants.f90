

    module Precision
    implicit none

    integer, parameter :: dl = KIND(1.d0)
    integer, parameter :: sp = KIND(1.0)
    end module Precision

    module constants
    use precision
    implicit none

    real(dl), parameter :: const_pi = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: const_twopi=2._dl*const_pi, const_fourpi=4._dl*const_pi
    real(dl), parameter :: const_eightpi=8._dl*const_pi
    real(dl), parameter :: const_sqrt6=2.4494897427831780981972840747059_dl
    real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
    real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
    real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

    !Constants from particle data book 2022

    real(dl), parameter :: c = 2.99792458e8_dl
    real(dl), parameter :: h_P = 6.62607015e-34_dl
    real(dl), parameter :: hbar = h_P/const_twopi

    real(dl), parameter :: G = 6.67430e-11_dl ! last two digit +/-15
    real(dl), parameter :: sigma_thomson = 6.6524587321e-29_dl
    real(dl), parameter :: sigma_boltz = 5.670374419e-8_dl
    real(dl), parameter :: k_B = 1.380649e-23_dl
    real(dl), parameter :: eV = 1.602176634e-19_dl

    real(dl), parameter :: m_p = 1.67262192369e-27_dl
    real(dl), parameter :: m_e = 9.1093837015e-31_dl

    !Isotopic masses from https://www.ciaaw.org/
    real(dl), parameter :: m_H = 1.673575e-27_dl !av. H atom for D/H=2.527e-5
    real(dl), parameter :: m_He4 = 6.646479073e-27_dl
    real(dl), parameter :: mass_ratio_He_H = m_He4/m_H

    real(dl), parameter :: Gyr = 365.25e9_dl*86400 !Julian gigayear
    real(dl), parameter :: Mpc = 3.085677581e22_dl
    real(dl), parameter :: MPC_in_sec = Mpc/c ! Mpc/c = 1.029271250e14 in SI units

    real(dl), parameter :: barssc0= k_B / m_p / c**2
    real(dl), parameter :: kappa=8._dl*const_pi*G
    real(dl), parameter :: a_rad = 8._dl*const_pi**5*k_B**4/15/c**3/h_p**3
    !7.565914e-16_dl !radiation constant for u=aT^4


    real(dl), parameter :: Compton_CT = MPC_in_sec*(8.d0/3.d0)*(sigma_thomson/(m_e*c))*a_rad
    !Compton_CT is CT in Mpc units, (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc
    !Used to get evolution of matter temperature

    !For 21cm
    real(dl), parameter :: f_21cm = 1420.40575e6_dl, l_21cm= c/f_21cm, T_21cm = h_P*f_21cm/k_B
    real(dl), parameter :: A10 = 2.869e-15, B10 = l_21cm**3/2/h_P/c*A10

    real(dl), parameter :: line21_const = 3*l_21cm**2*C*h_P/32/const_pi/k_B*A10 * Mpc_in_sec * 1000
    !1000 to get in MiliKelvin
    real(dl), parameter :: COBE_CMBTemp = 2.7255_dl !(Fixsen 2009) used as default value

    ! zeta3*3/2/pi^2*4/11*((k_B*COBE_CMBTemp/hbar/c)^3* 8*pi*G/3/(100*km/s/megaparsec)^2/(c^2/eV)
    real(dl), parameter :: inv_neutrino_mass_fac = zeta3*3._dl/2/const_pi**2*4._dl/11*((k_B*COBE_CMBTemp/hbar/c)**3* &
        kappa/3/(100*1000/Mpc)**2/(c**2/eV))
    !converts omnuh2 into sum m_nu in eV for non-relativistic but thermal neutrinos (no 0.046 factor); ~ 94.07
    real(dl), parameter :: neutrino_mass_fac= 1/inv_neutrino_mass_fac

    !Effective number of neutrinos, arXiv:2008.01074, arXiv:2012.02726
    real(dl), parameter :: default_nnu = 3.044_dl
    !Neutrino mass splittings
    real(dl), parameter :: delta_mnu21 = 7.54e-5_dl !eV^2 Particle Data Group 2015 (-0.22, + 0.26)
    real(dl), parameter :: delta_mnu31 = 2.46e-3_dl !eV^2 Particle Data Group 2015 (+- 0.06)
    !Round up to 0.06, so consistent with cosmomc's 1 neutrino default
    real(dl), parameter :: mnu_min_normal = 0.06_dl ! sqrt(delta_mnu31)+sqrt(delta_mnu21)

    !This constant initially was in inifile.f90.
    integer, parameter :: Ini_max_string_len = 1024
    end module constants




