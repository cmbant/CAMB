

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
    real(dl), parameter :: const_sqrt6=2.4494897427831780981972840747059_dl

    real(dl), parameter :: c = 2.99792458e8_dl
    real(dl), parameter :: h_P = 6.62606896e-34_dl
    real(dl), parameter :: hbar = h_P/const_twopi

    real(dl), parameter :: G=6.6738e-11_dl !data book 2012, last digit +/-8
    real(dl), parameter :: sigma_thomson = 6.6524616e-29_dl
    real(dl), parameter :: sigma_boltz = 5.6704e-8_dl
    real(dl), parameter :: k_B = 1.3806504e-23_dl
    real(dl), parameter :: eV = 1.60217646e-19_dl


    real(dl), parameter :: m_p = 1.672621637e-27_dl ! 1.672623e-27_dl
    real(dl), parameter :: m_H = 1.673575e-27_dl !av. H atom
    real(dl), parameter :: m_e = 9.10938215e-31_dl
    real(dl), parameter :: mass_ratio_He_H = 3.9715_dl


    real(dl), parameter :: Gyr=3.1556926e16_dl
    real(dl), parameter :: Mpc = 3.085678e22_dl !seem to be different definitions of this?
    real(dl), parameter :: MPC_in_sec = Mpc/c ! Mpc/c = 1.029272d14 in SI units

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
    real(dl), parameter :: default_nnu = 3.046_dl
    !Neutrino mass splittings
    real(dl), parameter :: delta_mnu21 = 7.54e-5_dl !eV^2 Particle Data Group 2015 (-0.22, + 0.26)
    real(dl), parameter :: delta_mnu31 = 2.46e-3_dl !eV^2 Particle Data Group 2015 (+- 0.06)
    !Round up to 0.06, so consistent with cosmomc's 1 neutrino default
    real(dl), parameter :: mnu_min_normal = 0.06_dl ! sqrt(delta_mnu31)+sqrt(delta_mnu21)

    end module constants


    module Errors
    implicit none

    integer :: global_error_flag=0
    character(LEN=1024) :: global_error_message = ''
    integer, parameter :: error_reionization=1
    integer, parameter :: error_recombination=2
    integer, parameter :: error_inital_power=3
    integer, parameter :: error_evolution=4
    integer, parameter :: error_unsupported_params=5

    contains

    subroutine GlobalError(message, id)
    character(LEN=*), intent(IN), optional :: message
    integer, intent(in), optional :: id

    if (present(message)) then
        global_error_message = message
    else
        global_error_message=''
    end if
    if (present(id)) then
        if (id==0) error stop 'Error id must be non-zero'
        global_error_flag=id
    else
        global_error_flag=-1
    end if

    end subroutine GlobalError

    end module Errors

