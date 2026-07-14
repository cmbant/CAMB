    !Recombination module for CAMB, using RECFAST

    !-------------------------------------------------------------------------------
    ! Integrator for Cosmic Recombination of Hydrogen and Helium,
    ! developed by Douglas Scott (dscott@astro.ubc.ca)
    ! based on calculations in the paper Seager, Sasselov & Scott
    ! (ApJ, 523, L1, 1999).
    !and "fudge" updates in Wong, Moss & Scott (2008).
    !
    ! Permission to use, copy, modify and distribute without fee or royalty at
    ! any tier, this software and its documentation, for any purpose and without
    ! fee or royalty is hereby granted, provided that you agree to comply with
    ! the following copyright notice and statements, including the disclaimer,
    ! and that the same appear on ALL copies of the software and documentation,
    ! including modifications that you make for internal use or for distribution:
    !
    ! Copyright 1999-2010 by University of British Columbia.  All rights reserved.
    !
    ! THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
    ! REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
    ! BY WAY OF EXAMPLE, BUT NOT LIMITATION,
    ! U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
    ! MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
    ! THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
    ! ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !-------------------------------------------------------------------------------
    !
    !CN     Name:        RECFAST
    !C
    !CP     Purpose:  Calculate ionised fraction as a function of redshift.
    !CP            Solves for H and He simultaneously, and includes
    !CP           H "fudge factor" for low z effect, as well as
    !CP           HeI fudge factor.
    !C
    !CD     Description: Solves for ionisation history since recombination
    !CD     using the equations in Seager, Sasselov & Scott (ApJ, 1999).
    !CD     The Cosmological model can be flat or open.
    !CD  The matter temperature is also followed, with an update from
    !CD  Scott & Scott (2009).
    !CD  The values for \alpha_B for H are from Hummer (1994).
    !CD  The singlet HeI coefficient is a fit from the full code.
    !CD  Additional He "fudge factors" are as described in Wong, Moss
    !CD  and Scott (2008).
    !CD  Extra fitting function included (in optical depth) to account
    !CD  for extra H physics described in Rubino-Martin et al. (2010).
    !CD  Care is taken to use the most accurate constants.
    !C
    !CA     Arguments:
    !CA     Name, Description
    !CA     real(dl) throughout
    !CA
    !CA     z is redshift - W is sqrt(1+z), like conformal time
    !CA     x is total ionised fraction, relative to H
    !CA     x_H is ionized fraction of H - y(1) in R-K routine
    !CA     x_He is ionized fraction of He - y(2) in R-K routine
    !CA       (note that x_He=n_He+/n_He here and not n_He+/n_H)
    !CA     Tmat is matter temperature
    !CA     aTmat=a*Tmat=Tmat/(1+z) is the stored temperature variable - y(3) in R-K routine
    !CA     f's are the derivatives of the Y's
    !CA     alphaB is case B recombination rate
    !CA     alpHe is the singlet only HeII recombination rate
    !CA     a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     a_VF is Verner and Ferland type fitting parameter for Helium
    !CA     b_VF is Verner and Ferland type fitting parameter for Helium
    !CA     T_0 is Verner and Ferland type fitting parameter for Helium
    !CA     T_1 is Verner and Ferland type fitting parameter for Helium
    !CA     Tnow is the observed CMB temperature today
    !CA     Yp is the primordial helium abundace
    !CA     fHe is He/H number ratio = Yp/4(1-Yp)
    !CA     Trad and Tmat are radiation and matter temperatures
    !CA     epsilon is the approximate difference (=Trad-Tmat) at high z
    !CA     OmegaB is Omega in baryons today
    !CA     H is Hubble constant in units of 100 km/s/Mpc
    !CA     HO is Hubble constant in SI units
    !CA     bigH is 100 km/s/Mpc in SI units
    !CA     Hz is the value of H at the specific z (in ION)
    !CA     G is grvitational constant
    !CA     n is number density of hydrogen
    !CA     Nnow is number density today
    !CA     x0 is initial ionized fraction
    !CA     x_H0 is initial ionized fraction of Hydrogen
    !CA     x_He0 is initial ionized fraction of Helium
    !CA     rhs is dummy for calculating x0
    !CA     zinitial and zfinal are starting and ending redshifts
    !CA     zeq is the redshift of matter-radiation equality
    !CA     zstart and zend are for each pass to the integrator
    !CA     C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
    !CA     m_e,m_H: electron mass and mass of H atom in SI
    !CA     not4: ratio of 4He atomic mass to 1H atomic mass
    !CA     sigma: Thomson cross-section
    !CA     a_rad: radiation constant for u=aT^4
    !CA     Lambda: 2s-1s two photon rate for Hydrogen
    !CA     Lambda_He: 2s-1s two photon rate for Helium
    !CA     DeltaB: energy of first excited state from continuum = 3.4eV
    !CA     DeltaB_He: energy of first excited state from cont. for He = 3.97eV
    !CA     L_H_ion: level for H ionization in m^-1
    !CA     L_H_alpha: level for H Ly alpha in m^-1
    !CA     L_He1_ion: level for HeI ionization
    !CA     L_He2_ion: level for HeII ionization
    !CA     L_He_2s: level for HeI 2s
    !CA     L_He_2p: level for HeI 2p (21P1-11S0) in m^-1
    !CA     Lalpha: Ly alpha wavelength in SI
    !CA     Lalpha_He: Helium I 2p-1s wavelength in SI
    !CA     mu_H,mu_T: mass per H atom and mass per particle
    !CA     H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
    !CA     CDB=DeltaB/k_B                     Constants derived from B1,B2,R
    !CA     CDB_He=DeltaB_He/k_B  n=2-infinity for He in Kelvin
    !CA     CB1=CDB*4.         Lalpha and sigma_Th, calculated
    !CA     CB1_He1: CB1 for HeI ionization potential
    !CA     CB1_He2: CB1 for HeII ionization potential
    !CA     CR=2*Pi*(m_e/h_P)*(k_B/h_P)  once and passed in a common block
    !CA     CK=Lalpha**3/(8.*Pi)
    !CA     CK_He=Lalpha_He**3/(8.*Pi)
    !CA     CL=C*h_P/(k_B*Lalpha)
    !CA     CL_He=C*h_P/(k_B*Lalpha_He)
    !CA     CT=(8./3.)*(sigma/(m_e*C))*a
    !CA     Bfact=exp((E_2p-E_2s)/kT)    Extra Boltzmann factor
    !CA b_He= "fudge factor" for HeI, to approximate higher z behaviour
    !CA Old Heswitch=6 HeI recombination physics is always used
    !CA Parameters and quantities to describe the extra triplet states
    !CA  and also the continuum opacity of H, with a fitting function
    !CA  suggested by KIV, astro-ph/0703438
    !CA a_trip: used to fit HeI triplet recombination rate
    !CA b_trip: used to fit HeI triplet recombination rate
    !CA L_He_2Pt: level for 23P012-11S0 in m^-1
    !CA L_He_2St: level for 23S1-11S0 in m^-1
    !CA L_He2St_ion: level for 23S1-continuum in m^-1
    !CA A2P_s: Einstein A coefficient for He 21P1-11S0
    !CA A2P_t: Einstein A coefficient for He 23P1-11S0
    !CA sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
    !CA sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
    !CA CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
    !CA CfHe_t: triplet statistical correction
    !CA Hswitch is an boolean for modifying the H recombination
    !CA AGauss1 is the amplitude of the 1st Gaussian for the H fudging
    !CA AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
    !CA zGauss1 is the ln(1+z) central value of the 1st Gaussian
    !CA zGauss2 is the ln(1+z) central value of the 2nd Gaussian
    !CA wGauss1 is the width of the 1st Gaussian
    !CA wGauss2 is the width of the 2nd Gaussian


    !CA     tol: tolerance for the integrator
    !CA     rk_settings,w(3,9): work space for RungeKuttaDP45
    !CA     Ndim: number of d.e.'s to solve (integer)
    !CA     Nz: number of output redshitf (integer)
    !CA     I: loop index (integer)
    !CA     ind,nw: work-space for RungeKuttaDP45 (integer)
    !C
    !CF     File & device access:
    !CF     Unit /I,IO,O  /Name (if known)
    !C
    !CM     Modules called:
    !CM     RungeKuttaDP45 (numerical integrator)
    !CM     GET_INIT (initial values for ionization fractions)
    !CM     ION (ionization and Temp derivatives)
    !C
    !CC     Comments:
    !CC     none
    !C
    !CH     History:
    !CH     CREATED            (simplest version) 19th March 1989
    !CH     RECREATED    11th January 1995
    !CH               includes variable Cosmology
    !CH               initial conditions are Saha
    !CH     MODIFIED     January 1995 (include Hummer's 1994 alpha table)
    !CH               January 1995 (include new value for 2s-1s rate)
    !CH               January 1995 (expand comments)
    !CH               March 1995 (add Saha for Helium)
    !CH               August 1997 (add HeII alpha table)
    !CH               July 1998 (include OmegaT correction and H fudge factor)
    !CH               Nov 1998 (change Trad to Tmat in Rup)
    !CH               Jan 1999 (tidied up for public consumption)
    !CH               Sept 1999 (switch to formula for alpha's, fix glitch)
    !CH                  Sept 1999 modified to CMBFAST by US & MZ
    !CH                     Nov 1999 modified for F90 and CAMB (AML)
    !CH                     Aug 2000 modified to prevent overflow erorr in He_Boltz (AML)
    !CH                     Feb 2001 corrected fix of Aug 2000 (AML)
    !CH                     Oct 2001 fixed error in hubble parameter, now uses global function (AML)
    !                       March 2003 fixed bugs reported by savita gahlaut
    !                       March 2005 added option for corrections from astro-ph/0501672.
    !                                  thanks to V.K.Dubrovich, S.I.Grachev
    !                       June 2006 defined RECFAST_fudge as free parameter (AML)
    !                       October 2006 (included new value for G)
    !                       October 2006 (improved m_He/m_H to be "not4")
    !                       October 2006 (fixed error, x for x_H in part of f(1))
    !CH              January 2008 (improved HeI recombination effects,
    !CH                       including HeI rec. fudge factor)
    !                Feb 2008   Recfast 1.4 changes above added (AML)
    !                           removed Dubrovich option (wrong anyway)
    !CH              Sept 2008 (added extra term to make transition, smoother for Tmat evolution)
    !                Sept 2008 Recfast 1.4.2 changes above added (AML)
    !                          General recombination module structure, fix to make He x_e smooth also in recfast (AML)
    !CH              Jan 2010 (added fitting function to modify K to match x_e(z) for new H physics)
    !AL              June 2012 updated fudge parameters to match HyRec and CosmoRec (AML)
    !AL              Sept 2012 changes now in public recfast, version number changed to match Recfast 1.5.2.
    !AL              Apr 2026 updated to use variable nZ and evolve/interpolate a*T_m with Rosenbrock while stiff,
    !AL                      and added optional HeI recombination rate correction fit. New CosmoRec fudges.


    module Recombination
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use constants
    use classes
    use DarkAge21cm
    use Interpolation, only : cubic_spline_regular_horner_coefficients, cubic_spline_regular_second_derivs
    use MathUtils
    use RungeKuttaDP45Module, only : RungeKuttaDP45Settings, TClassRungeKuttaDP45
    use Config, only : GlobalError, error_recombination
    use results
    use MpiUtils, only : MpiStop
    implicit none
    private

    real(dl), parameter ::  zinitial = 1e4_dl !highest redshift
    real(dl), parameter ::  zfinal=0._dl

    real(dl), parameter :: RECFAST_fudge_He_default = 0.8472367977_dl !Helium fudge parameter
    logical, parameter  :: RECFAST_Hswitch_default = .true. !include H corrections (v1.5, 2010)
    logical, parameter  :: RECFAST_He_rate_correction_default = .true.
    real(dl), parameter :: RECFAST_He_rate_zmin = 1500._dl
    real(dl), parameter :: RECFAST_He_rate_zmax = 3000._dl
    real(dl), parameter :: RECFAST_He_rate_a0 = 0.07805480599148856_dl
    real(dl), parameter :: RECFAST_He_rate_a1 = 0.899368710079757_dl
    real(dl), parameter :: RECFAST_He_rate_a2 = -3.063896868095767_dl
    real(dl), parameter :: RECFAST_He_rate_d2 = 3.752520187028219_dl
    real(dl), parameter :: RECFAST_He_rate_d4 = 18.126758791787946_dl
    real(dl), parameter :: RECFAST_He_rate_z0 = 1968.1134791219417_dl
    real(dl), parameter :: RECFAST_He_rate_width = 774.0312440582285_dl
    real(dl), parameter :: RECFAST_fudge_default = 1.14_dl !1.14_dl
    real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 + 0.02d0
    integer, parameter :: RECFAST_nz_default = 2046
    real(dl), parameter :: RECFAST_x_He_freeze_threshold = 1.e-8_dl
    logical, parameter :: RECFAST_use_rosenbrock_default = .true.
    real(dl), parameter :: RECFAST_rosenbrock_handoff_xH_default = 0.976_dl

    real(dl), parameter :: ROS2_gamma = 1._dl + 1._dl/sqrt(2._dl)
    real(dl), parameter :: ROS2_a21 = 1._dl/ROS2_gamma
    real(dl), parameter :: ROS2_c21 = -2._dl/ROS2_gamma
    real(dl), parameter :: ROS2_m1 = 3._dl/(2._dl*ROS2_gamma)
    real(dl), parameter :: ROS2_m2 = 1._dl/(2._dl*ROS2_gamma)
    real(dl), parameter :: ROS2_safety = 0.9_dl
    real(dl), parameter :: ROS2_min_scale = 0.2_dl
    real(dl), parameter :: ROS2_max_scale = 2.5_dl
    ! Keep the ionization fractions on a smaller absolute scale than temperature
    ! so full-range Rosenbrock does not tolerate O(1) relative tail errors.
    real(dl), parameter :: RECFAST_rosenbrock_ion_scale_floor = 1.e-3_dl
    ! Tuned for the fast handoff path; full-range Rosenbrock should use a tighter tol.
    real(dl), parameter :: RECFAST_rosenbrock_tol_default = 3.e-4_dl
    integer, parameter :: RECFAST_rosenbrock_max_steps = 4096

    Type RecombinationData
        real(dl) :: Recombination_saha_z !Redshift at which saha OK
        real(dl), private :: NNow, fHe
        integer, private :: nz = 0
        real(dl), private :: delta_z = 0._dl, minz = zfinal, maxz = zinitial
        real(dl), allocatable, private :: zrec(:), xrec(:), dxrec(:), Tsrec(:), dTsrec(:), tmrec(:), dtmrec(:), &
            xrec_horner(:, :), tsrec_horner(:, :), tmrec_horner(:, :)
        ! tmrec stores a*Tmat = Tmat/(1+z)
        ! tsrec stores a*Tspin = Tspin/(1+z)
        real(dl), private :: DeltaB,DeltaB_He,Lalpha,mu_H,mu_T

        real(dl), private :: HO, Tnow, fu
        integer, private :: n_eq = 3
        logical :: doTspin = .false.

        !The following only used for approximations where small effect
        real(dl) :: OmegaK, OmegaT, z_eq
        class(CAMBdata), pointer :: State
    end Type RecombinationData

    type, extends(TRecombinationModel) :: TRecfast
        real(dl) :: RECFAST_fudge  = RECFAST_fudge_default2
        real(dl) :: RECFAST_fudge_He = RECFAST_fudge_He_default
        logical  :: RECFAST_Hswitch = RECFAST_Hswitch_default
        logical  :: RECFAST_He_rate_correction = RECFAST_He_rate_correction_default
        ! HeI recombination always uses the old RECFAST_Heswitch=6 physics:
        ! singlet Sobolev escape, H continuum opacity in the singlet channel,
        ! triplet recombination, and H continuum opacity in the triplet channel.

        !fudge parameter if RECFAST_Hswitch
        !Gaussian fits for extra H physics (fit to CosmoRec accuracy6 H10, 2026)
        real(dl) :: AGauss1 =      -0.1395272483D0  !Amplitude of 1st Gaussian
        real(dl) :: AGauss2 =       0.0729891952D0  !Amplitude of 2nd Gaussian
        real(dl) :: zGauss1 =       7.2813061282D0  !ln(1+z) of 1st Gaussian
        real(dl) :: zGauss2=        6.7667038679D0  !ln(1+z) of 2nd Gaussian
        real(dl) :: wGauss1=        0.1638966410D0  !Width of 1st Gaussian
        real(dl) :: wGauss2=        0.2785834127D0  !Width of 2nd Gaussian
        integer  :: Nz = RECFAST_nz_default
        logical  :: use_rosenbrock = RECFAST_use_rosenbrock_default
        real(dl) :: rosenbrock_handoff_xH = RECFAST_rosenbrock_handoff_xH_default
        real(dl) :: rosenbrock_tol = RECFAST_rosenbrock_tol_default
        Type(RecombinationData), allocatable :: Calc
    contains
    procedure :: ReadParams => TRecfast_ReadParams
    procedure :: Validate => TRecfast_Validate
    procedure :: Init => TRecfast_init
    procedure :: x_e => TRecfast_xe
    procedure :: xe_Tm => TRecfast_xe_Tm !ionization fraction and baryon temperature
    procedure :: T_m => TRecfast_tm !baryon temperature
    procedure :: T_s => TRecfast_ts !Spin temperature
    procedure :: Version => TRecfast_version
    procedure :: dDeltaxe_dtau => TRecfast_dDeltaxe_dtau
    procedure :: get_Saha_z => TRecfast_Get_Saha_z
    procedure, nopass :: SelfPointer => TRecfast_SelfPointer

    end type TRecfast

    character(LEN=*), parameter :: Recfast_Version = 'Recfast_2.0.0'

    logical, parameter :: evolve_Ts = .false. !local equilibrium is very accurate
    real(dl), parameter :: Do21cm_minev = 1/(1+400.) !at which to evolve T_s

    real(dl), parameter :: bigH=100.0D3/Mpc !Ho in s-1
    real(dl), parameter :: sigma = sigma_thomson
    real(dl), parameter :: not4  = mass_ratio_He_H    !mass He/H atom

    real(dl), parameter :: B01 = 3*B10

    real(dl), parameter ::  Lambda = 8.2245809d0
    real(dl), parameter :: Lambda_He = 51.3d0    !new value from Dalgarno
    real(dl), parameter :: L_H_ion   = 1.096787737D7 !level for H ion. (in m^-1)
    real(dl), parameter :: L_H_alpha = 8.225916453D6 !averaged over 2 levels
    real(dl), parameter :: L_He1_ion = 1.98310772D7  !from Drake (1993)
    real(dl), parameter :: L_He2_ion = 4.389088863D7 !from JPhysChemRefData (1987)
    real(dl), parameter :: L_He_2s   = 1.66277434D7  !from Drake (1993)
    real(dl), parameter :: L_He_2p   = 1.71134891D7  !from Drake (1993)
    !   2 photon rates and atomic levels in SI units

    real(dl), parameter :: A2P_s     = 1.798287D9    !Morton, Wu & Drake (2006)
    real(dl), parameter :: A2P_t     = 177.58D0      !Lach & Pachuski (2001)
    real(dl), parameter :: L_He_2Pt  = 1.690871466D7 !Drake & Morton (2007)
    real(dl), parameter :: L_He_2St  = 1.5985597526D7 !Drake & Morton (2007)
    real(dl), parameter :: L_He2St_ion  =3.8454693845D6 !Drake & Morton (2007)
    real(dl), parameter :: sigma_He_2Ps  = 1.436289D-22  !Hummer & Storey (1998)
    real(dl), parameter :: sigma_He_2Pt  = 1.484872D-22  !Hummer & Storey (1998)
    !    Atomic data for HeI

    !       Set up some constants so they don't have to be calculated later
    real(dl), parameter :: Lalpha = 1.d0/L_H_alpha
    real(dl), parameter :: Lalpha_He = 1.d0/L_He_2p
    real(dl), parameter :: DeltaB = h_P*C*(L_H_ion-L_H_alpha)
    real(dl), parameter :: CDB = DeltaB/k_B
    real(dl), parameter :: DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)   !2s, not 2p
    real(dl), parameter :: CDB_He = DeltaB_He/k_B
    real(dl), parameter :: CB1 = h_P*C*L_H_ion/k_B
    real(dl), parameter :: CB1_He1 = h_P*C*L_He1_ion/k_B   !ionization for HeI
    real(dl), parameter :: CB1_He2 = h_P*C*L_He2_ion/k_B   !ionization for HeII
    real(dl), parameter :: CR = const_twopi*(m_e/h_P)*(k_B/h_P)
    real(dl), parameter :: CK = Lalpha**3/(const_eightpi)
    real(dl), parameter :: CK_He = Lalpha_He**3/(const_eightpi)
    real(dl), parameter :: CL = C*h_P/(k_B*Lalpha)
    real(dl), parameter :: CL_He = C*h_P/(k_B/L_He_2s) !comes from det.bal. of 2s-1s
    real(dl), parameter :: CT = Compton_CT / MPC_in_sec

    real(dl), parameter :: Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

    !       Matter departs from radiation when t(Th) > H_frac * t(H)
    !       choose some safely small number
    real(dl), parameter :: H_frac = 1D-3

    procedure(obj_function), private :: dtauda

    public TRecfast,  CB1

    contains

    subroutine TRecfast_ReadParams(this, Ini)
    use IniObjects
    class(TRecfast) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: legacy_heswitch

    this%RECFAST_fudge_He = Ini%Read_Double('RECFAST_fudge_He', RECFAST_fudge_He_default)
    if (Ini%HasKey('RECFAST_Heswitch')) then
        legacy_heswitch = Ini%Read_Int('RECFAST_Heswitch', 6)
        if (legacy_heswitch /= 6) then
            call GlobalError("RECFAST_Heswitch is no longer configurable; only RECFAST_Heswitch=6 is supported", &
                error_recombination)
            return
        end if
    end if
    this%RECFAST_Hswitch = Ini%Read_Logical('RECFAST_Hswitch', RECFAST_Hswitch_default)
    this%RECFAST_He_rate_correction = Ini%Read_Logical('RECFAST_He_rate_correction', &
        RECFAST_He_rate_correction_default)
    if (Ini%HasKey('RECFAST_H_fudge')) then
        this%RECFAST_fudge = Ini%Read_Double('RECFAST_H_fudge', RECFAST_fudge_default2)
    else
        this%RECFAST_fudge = Ini%Read_Double('RECFAST_fudge', RECFAST_fudge_default)
        if (this%RECFAST_Hswitch) then
            this%RECFAST_fudge = this%RECFAST_fudge - (RECFAST_fudge_default - RECFAST_fudge_default2)
        end if
    end if
    call Ini%Read('AGauss1',this%AGauss1)
    call Ini%Read('AGauss2',this%AGauss2)
    call Ini%Read('zGauss1',this%zGauss1)
    call Ini%Read('zGauss2',this%zGauss2)
    call Ini%Read('wGauss1',this%wGauss1)
    call Ini%Read('wGauss2',this%wGauss2)
    this%Nz = Ini%Read_Int("RECFAST_nz", this%Nz)
    this%use_rosenbrock = Ini%Read_Logical("RECFAST_use_rosenbrock", this%use_rosenbrock)
    this%rosenbrock_handoff_xH = Ini%Read_Double("RECFAST_rosenbrock_handoff_xH", this%rosenbrock_handoff_xH)
    this%rosenbrock_tol = Ini%Read_Double("RECFAST_rosenbrock_tol", this%rosenbrock_tol)
    end subroutine TRecfast_ReadParams

    subroutine TRecfast_Validate(this, OK)
    class(TRecfast),intent(in) :: this
    logical, intent(inout) :: OK

    if (this%Nz < 2) then
        OK = .false.
        write(*,*) "RECFAST_nz must be at least 2"
    end if
    if (this%rosenbrock_handoff_xH < 0._dl .or. this%rosenbrock_handoff_xH >= 1._dl) then
        OK = .false.
        write(*,*) "RECFAST_rosenbrock_handoff_xH must be in [0, 1)"
    end if
    if (this%rosenbrock_tol <= 0._dl) then
        OK = .false.
        write(*,*) "RECFAST_rosenbrock_tol must be > 0"
    end if
    end subroutine TRecfast_Validate


    function TRecfast_tm(this,a)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl) z, zst, TRecfast_tm, aTmat_stored, az
    integer ihi

    z=1/a-1
    associate( Calc => this%Calc)
        if (z >= Calc%maxz) then
            TRecfast_tm=Calc%Tnow/a
        else
            if (z <= Calc%minz) then
                TRecfast_tm=(1._dl + zfinal)*Calc%Tmrec(Calc%nz)
            else
                zst = (zinitial - z)/Calc%delta_z
                ihi = int(zst)
                az = zst - real(ihi, dl)
                aTmat_stored = Calc%tmrec_horner(1, ihi) + az*(Calc%tmrec_horner(2, ihi) &
                    + az*(Calc%tmrec_horner(3, ihi) +     az*Calc%tmrec_horner(4, ihi)))
                TRecfast_tm=(1._dl + z)*aTmat_stored
            endif
        endif
    end associate

    end function TRecfast_tm


    function TRecfast_ts(this,a)
    class(TRecfast) :: this
    !zrec(1) is zinitial-delta_z
    real(dl), intent(in) :: a
    real(dl) zst,z,az,TRecfast_ts, aTspin_stored
    integer ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z >= Calc%maxz) then
            TRecfast_ts=(1._dl + z)*Calc%tsrec(1)
        else
            if (z <= Calc%minz) then
                TRecfast_ts=(1._dl + zfinal)*Calc%tsrec(Calc%nz)
            else
                zst = (zinitial - z)/Calc%delta_z
                ihi = int(zst)
                az = zst - real(ihi, dl)
                aTspin_stored = Calc%tsrec_horner(1, ihi) + az*(Calc%tsrec_horner(2, ihi) &
                    + az*(Calc%tsrec_horner(3, ihi) + az*Calc%tsrec_horner(4, ihi)))
                TRecfast_ts = (1._dl + z)*aTspin_stored
            endif
        endif
    end associate
    end function TRecfast_ts

    function TRecfast_xe(this,a)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl) zst,z,az,TRecfast_xe
    integer ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z >= Calc%maxz) then
            TRecfast_xe=Calc%xrec(1)
        else
            if (z <= Calc%minz) then
                TRecfast_xe=Calc%xrec(Calc%nz)
            else
                zst = (zinitial - z)/Calc%delta_z
                ihi = int(zst)
                az = zst - real(ihi, dl)
                TRecfast_xe = Calc%xrec_horner(1, ihi) + az*(Calc%xrec_horner(2, ihi) + az*(Calc%xrec_horner(3, ihi) + &
                    az*Calc%xrec_horner(4, ihi)))
            endif
        endif
    end associate
    end function TRecfast_xe

    subroutine TRecfast_xe_Tm(this,a, xe, Tm)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: xe, Tm
    real(dl) z, zst, aTmat_stored, az
    integer ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z >= Calc%maxz) then
            xe=Calc%xrec(1)
            Tm = Calc%Tnow/a
        else
            if (z <= Calc%minz) then
                xe=Calc%xrec(Calc%nz)
                Tm = (1._dl + zfinal)*Calc%Tmrec(Calc%nz)
            else
                zst = (zinitial - z)/Calc%delta_z
                ihi = int(zst)
                az = zst - real(ihi, dl)
                xe = Calc%xrec_horner(1, ihi) + az*(Calc%xrec_horner(2, ihi) + az*(Calc%xrec_horner(3, ihi) + &
                    az*Calc%xrec_horner(4, ihi)))
                aTmat_stored = Calc%tmrec_horner(1, ihi) + az*(Calc%tmrec_horner(2, ihi) &
                    + az*(Calc%tmrec_horner(3, ihi) +  az*Calc%tmrec_horner(4, ihi)))
                Tm=(1._dl + z)*aTmat_stored
            endif
        endif
    end associate
    end subroutine TRecfast_xe_Tm

    subroutine EnsureRecfastStorage(Calc, target_nz, OK)
    type(RecombinationData), intent(inout) :: Calc
    integer, intent(in) :: target_nz
    logical, intent(out) :: OK
    logical :: needs_allocate

    OK = .false.
    if (target_nz < 2) then
        call GlobalError("recfast_nz must be at least 2", error_recombination)
        return
    end if

    needs_allocate = .not. allocated(Calc%zrec)
    if (.not. needs_allocate) needs_allocate = size(Calc%zrec) /= target_nz

    if (needs_allocate .and. allocated(Calc%zrec)) then
        deallocate(Calc%zrec, Calc%xrec, Calc%dxrec, Calc%tsrec, Calc%dtsrec, Calc%tmrec, Calc%dtmrec, &
            Calc%xrec_horner, Calc%tsrec_horner, Calc%tmrec_horner)
    end if

    Calc%nz = target_nz
    Calc%delta_z = (zinitial-zfinal)/real(Calc%nz, dl)
    Calc%minz = zfinal + 1.e-12_dl*Calc%delta_z
    Calc%maxz = zinitial - Calc%delta_z - 1.e-12_dl*Calc%delta_z
    if (needs_allocate) then
        allocate(Calc%zrec(Calc%nz), Calc%xrec(Calc%nz), Calc%dxrec(Calc%nz), Calc%tsrec(Calc%nz), &
            Calc%dtsrec(Calc%nz), Calc%tmrec(Calc%nz), Calc%dtmrec(Calc%nz), Calc%xrec_horner(4, Calc%nz - 1), &
            Calc%tsrec_horner(4, Calc%nz - 1), Calc%tmrec_horner(4, Calc%nz - 1))
    end if
    OK = .true.

    end subroutine EnsureRecfastStorage

    function TRecfast_version(this) result(this_version)
    class(TRecfast) :: this
    character(LEN=:), allocatable :: this_version

    this_version = Recfast_Version

    end function TRecfast_version

    subroutine TRecfast_init(this,State, WantTSpin)
    use MiscUtils
    implicit none
    class(TRecfast), target :: this
    class(TCAMBdata), target :: State
    real(dl) :: Trad,Tmat,Tspin, ainv
    integer :: I
    Type(RecombinationData), pointer :: Calc
    logical, intent(in), optional :: WantTSpin
    real(dl) :: z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0,H, Yp
    real(dl) :: zstart,zend,z_scale, rosenbrock_tol, rk45_tol_use, background_step_boost
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl), dimension(:,:), allocatable :: w
    real(dl) :: y(4), y_rosen_start(4)
    real(dl) :: C10, tau_21Ts
    integer :: ind, nw, internal_nz
    real(dl), parameter :: rk45_tol=3D-6            !Input tolerance for RungeKuttaDP45
    procedure(TClassRungeKuttaDP45) :: RungeKuttaDP45
    logical :: storage_ok, rosenbrock_handed_off, rosenbrock_ok


    if (.not. allocated(this%Calc)) allocate(this%Calc)
    Calc => this%Calc

    select type(State)
    class is (CAMBdata)
        Calc%State => State
        Calc%doTspin = DefaultFalse(WantTSpin)
        background_step_boost = max(State%CP%Accuracy%BackgroundTimeStepBoost, 1.e-12_dl)
        internal_nz = max(2, nint(this%Nz*background_step_boost))
        call EnsureRecfastStorage(Calc, internal_nz, storage_ok)
        if (.not. storage_ok) return

        Calc%n_eq = 3
        if (Evolve_Ts) Calc%n_eq=4
        allocate(w(Calc%n_eq,9))

        Calc%Recombination_saha_z=0.d0

        Calc%Tnow = State%CP%tcmb

        z = zinitial
        ainv = 1._dl + z

        H = State%CP%H0/100._dl

        !Not general, but only for approx
        Calc%OmegaT=(State%CP%omch2+State%CP%ombh2)/H**2        !total dark matter + baryons
        Calc%OmegaK=State%CP%omk       !curvature


        !       convert the Hubble constant units
        Calc%HO = H*bigH
        Yp = State%CP%Yhe

        !       sort out the helium abundance parameters
        Calc%mu_H = 1.d0/(1.d0-Yp)           !Mass per H atom
        Calc%mu_T = not4/(not4-(not4-1.d0)*Yp)   !Mass per atom
        Calc%fHe = Yp/(not4*(1.d0-Yp))       !n_He_tot / n_H_tot


        Calc%Nnow = 3._dl*bigH**2*State%CP%ombh2/(const_eightpi*G*Calc%mu_H*m_H)

        n = Calc%Nnow*ainv**3
        Calc%z_eq = State%z_eq

        !       Fudge factor to approximate for low z out of equilibrium effect
        Calc%fu=this%RECFAST_fudge
        rosenbrock_tol = this%rosenbrock_tol/State%CP%Accuracy%IntTolBoost
        rk45_tol_use = rk45_tol/State%CP%Accuracy%IntTolBoost

        !       Set initial matter temperature. y(3) stores a*Tmat directly.
        Tmat = Calc%Tnow*ainv                 !Initial rad. & mat. temperature
        y(3) = Calc%Tnow
        y(4) = Calc%Tnow
        Tspin = Tmat

        call get_init(Calc,z,x_H0,x_He0,x0)

        y(1) = x_H0
        y(2) = x_He0

        !       OK that's the initial conditions, now start writing output file

        !       Set up work-space stuff for RungeKuttaDP45
        ind  = 1
        nw   = Calc%n_eq
        rosenbrock_handed_off = .not. this%use_rosenbrock
        rk_settings%reuse_first_derivative = .true.

        do i = 1,Calc%nz
            !       calculate the start and end redshift for the interval at each z
            !       or just at each z
            zstart = zinitial  - real(i-1,dl)*Calc%delta_z
            zend   = zinitial  - real(i,dl)*Calc%delta_z

            ! Use Saha to get x_e, using the equation for x_e for ionized helium
            ! and for neutral helium.
            ! Everything ionized above z=8000.  First ionization over by z=5000.
            ! Assume He all singly ionized down to z=3500, then use He Saha until
            ! He is 99% singly ionized, and *then* switch to joint H/He recombination.

            z = zend
            ainv = 1._dl + z
            z_scale = Calc%Tnow/COBE_CMBTemp*ainv - 1
            if (.not. rosenbrock_handed_off .and. y(1) <= this%rosenbrock_handoff_xH) then
                rosenbrock_handed_off = .true.
            end if

            if (z_scale > 8000._dl) then

                x_H0 = 1._dl
                x_He0 = 1._dl
                x0 = 1._dl+2._dl*Calc%fHe
                y(1) = x_H0
                y(2) = x_He0
                Tmat = Calc%Tnow*ainv
                y(3) = Calc%Tnow
                y(4) = Calc%Tnow

            else if(z_scale > 5000._dl)then

                x_H0 = 1._dl
                x_He0 = 1._dl
                rhs = exp(1.5d0*log(CR*Calc%Tnow/ainv) - CB1_He2/(Calc%Tnow*ainv)) / Calc%Nnow
                rhs = rhs*1._dl            !ratio of g's is 1 for He++ <-> He+
                x0 = 0.5d0 * ( sqrt( (rhs-1._dl-Calc%fHe)**2 &
                    + 4._dl*(1._dl+2._dl*Calc%fHe)*rhs) - (rhs-1._dl-Calc%fHe) )
                y(1) = x_H0
                y(2) = x_He0
                Tmat = Calc%Tnow*ainv
                y(3) = Calc%Tnow
                y(4) = Calc%Tnow

            else if(z_scale > 3500._dl)then

                x_H0 = 1._dl
                x_He0 = 1._dl
                x0 = x_H0 + Calc%fHe*x_He0
                y(1) = x_H0
                y(2) = x_He0
                Tmat = Calc%Tnow*ainv
                y(3) = Calc%Tnow
                y(4) = Calc%Tnow

            else if(y(2) > 0.99)then

                x_H0 = 1._dl
                rhs = exp(1.5d0*log(CR*Calc%Tnow/ainv) - CB1_He1/(Calc%Tnow*ainv)) / Calc%Nnow
                rhs = rhs*4._dl            !ratio of g's is 4 for He+ <-> He0
                x_He0 = 0.5d0 * ( sqrt( (rhs-1._dl)**2 &
                    + 4._dl*(1._dl+Calc%fHe)*rhs )- (rhs-1._dl))
                x0 = x_He0
                x_He0 = (x0 - 1._dl)/Calc%fHe
                y(1) = x_H0
                y(2) = x_He0
                Tmat = Calc%Tnow*ainv
                y(3) = Calc%Tnow
                y(4) = Calc%Tnow

            else if (.not. rosenbrock_handed_off) then

                ! Integrate the full smooth H/He/Tm system with Rosenbrock until the
                ! hydrogen fraction has dropped enough to hand back to RungeKuttaDP45.
                y_rosen_start = y
                call RecfastRosenbrockAdvance(this, zstart, y(1:Calc%n_eq), zend, rosenbrock_tol, &
                    rosenbrock_ok)
                if (.not. rosenbrock_ok) then
                    call GlobalError("Recfast: Rosenbrock integration failed", error_recombination)
                    return
                end if
                if (.not. Evolve_Ts) y(4) = y(3)
                if (y(1) <= this%rosenbrock_handoff_xH) then
                    ! Snap the handoff to the higher-redshift grid node by redoing the
                    ! crossing interval with the original RungeKuttaDP45 evolution.
                    y = y_rosen_start
                    rosenbrock_handed_off = .true.
                    if (y(1) > 0.99d0) then
                        rhs = exp(1.5d0*log(CR*Calc%Tnow/ainv) - CB1/(Calc%Tnow*ainv)) / Calc%Nnow
                        x_H0 = 0.5d0 * (sqrt(rhs**2 + 4._dl*rhs) - rhs)

                        call RungeKuttaDP45(this, 3, ION, zstart, y, zend, rk45_tol_use, ind, rk_settings, nw, w)
                        y(1) = x_H0
                        x0 = y(1) + Calc%fHe*y(2)
                        y(4) = y(3)
                        rk_settings%first_derivative_valid = .false.
                    else
                        call RungeKuttaDP45(this, nw, ION, zstart, y, zend, rk45_tol_use, ind, rk_settings, nw, w)
                        x0 = y(1) + Calc%fHe*y(2)
                    end if
                else
                    x0 = y(1) + Calc%fHe*y(2)
                    if (y(1) > 0.985d0) Calc%Recombination_saha_z = zend
                end if

            else if (y(1) > 0.99d0) then

                rhs = exp(1.5d0*log(CR*Calc%Tnow/ainv) - CB1/(Calc%Tnow*ainv)) / Calc%Nnow
                x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )

                rk_settings%first_derivative_valid = .false.
                call RungeKuttaDP45(this, 3, ION, zstart, y, zend, rk45_tol_use, ind, rk_settings, nw, w)
                y(1) = x_H0
                x0 = y(1) + Calc%fHe*y(2)
                y(4)=y(3)
                rk_settings%first_derivative_valid = .false.
            else

                call RungeKuttaDP45(this, nw, ION, zstart, y, zend, rk45_tol_use, ind, rk_settings, nw, w)

                x0 = y(1) + Calc%fHe*y(2)

            end if

            ainv = 1._dl + zend
            Trad = Calc%Tnow*ainv
            Tmat = ainv*y(3)
            x_H = y(1)
            x_He = y(2)
            x = x0

            Calc%zrec(i) = zend
            Calc%xrec(i) = x
            Calc%tmrec(i) = y(3)


            if (Calc%doTspin) then
                if (Evolve_Ts .and. zend< 1/Do21cm_minev-1 ) then
                    Tspin = ainv*y(4)
                else
                    C10 = Calc%Nnow*ainv**3*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) &
                        + kappa_eH_21cm(Tmat,.false.)*x)
                    tau_21Ts = line21_const*Calc%NNow*ainv*dtauda(State,1/ainv)/1000

                    Tspin = Trad*( C10/Trad + A10/T_21cm)/(C10/Tmat + A10/T_21cm) + &
                        tau_21Ts/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) )

                    y(4) = Tspin/ainv
                end if

                Calc%tsrec(i) = y(4)

            end if

        end do

        call cubic_spline_regular_second_derivs(-Calc%delta_z,Calc%xrec,Calc%nz,Calc%dxrec)
        ! At low z, adiabatic cooling gives Tmat ~ (1+z)^2, so a*Tmat ~ (1+z).
        call cubic_spline_regular_second_derivs(-Calc%delta_z, Calc%tmrec, Calc%nz, Calc%dtmrec, &
            End2=Calc%tmrec(Calc%nz)/(1._dl + zfinal))
        call cubic_spline_regular_horner_coefficients(Calc%delta_z, Calc%xrec, Calc%dxrec, &
            Calc%xrec_horner, Calc%nz)
        call cubic_spline_regular_horner_coefficients(Calc%delta_z, Calc%tmrec, Calc%dtmrec, &
            Calc%tmrec_horner, Calc%nz)
        if (Calc%doTspin) then
            call cubic_spline_regular_second_derivs(-Calc%delta_z,Calc%tsrec,Calc%nz,Calc%dtsrec)
            call cubic_spline_regular_horner_coefficients(Calc%delta_z, Calc%tsrec, Calc%dtsrec, &
                Calc%tsrec_horner, Calc%nz)
        end if
    class default
        call MpiStop('Wrong state type')
    end select

    end subroutine TRecfast_init

    !       ===============================================================
    subroutine GET_INIT(Calc,z,x_H0,x_He0,x0)

    !       Set up the initial conditions so it will work for general,
    !       but not pathological choices of zstart
    !       Initial ionization fraction using Saha for relevant species
    Type(RecombinationData) :: Calc
    real(dl) z,x0,rhs,x_H0,x_He0, z_scale

    z_scale = Calc%Tnow/COBE_CMBTemp*(z+1)-1

    if(z_scale > 8000._dl)then
        x_H0 = 1._dl
        x_He0 = 1._dl
        x0 = 1._dl+2._dl*Calc%fHe

    else if(z_scale > 3500._dl)then

        x_H0 = 1._dl
        x_He0 = 1._dl
        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1_He2/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        rhs = rhs*1._dl    !ratio of g's is 1 for He++ <-> He+
        x0 = 0.5d0 * ( sqrt( (rhs-1._dl-Calc%fHe)**2 &
            + 4._dl*(1._dl+2._dl*Calc%fHe)*rhs) - (rhs-1._dl-Calc%fHe) )

    else if(z_scale > 2000._dl)then

        x_H0 = 1._dl
        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1_He1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        rhs = rhs*4._dl    !ratio of g's is 4 for He+ <-> He0
        x_He0 = 0.5d0  * ( sqrt( (rhs-1._dl)**2 + 4._dl*(1._dl+Calc%fHe)*rhs )- (rhs-1._dl))
        x0 = x_He0
        x_He0 = (x0 - 1._dl)/Calc%fHe

    else

        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )
        x_He0 = 0._dl
        x0 = x_H0
    end if

    end subroutine GET_INIT

    subroutine EscapeProbabilityAndDerivative(tau, p_escape, dp_dtau)
    real(dl), intent(in) :: tau
    real(dl), intent(out) :: p_escape, dp_dtau
    real(dl) :: tau2

    if (abs(tau) < 1.e-6_dl) then
        tau2 = tau*tau
        p_escape = 1._dl - tau/2._dl + tau2/6._dl
        dp_dtau = -0.5_dl + tau/3._dl - tau2/8._dl
    else
        p_escape = (1._dl - exp(-tau))/tau
        dp_dtau = (exp(-tau)*(tau + 1._dl) - 1._dl)/(tau*tau)
    end if

    end subroutine EscapeProbabilityAndDerivative

    subroutine FactorSmallLinearSystem(a, ipiv, ok)
    !In-place LU factorization with partial pivoting; multipliers stored below the diagonal
    real(dl), intent(inout) :: a(:, :)
    integer, intent(out) :: ipiv(:)
    logical, intent(out) :: ok
    real(dl) :: factor, maxabs, temp, pivot_scale
    integer :: i, j, k, n, pivot

    n = size(ipiv)
    ok = .false.
    if (.not. all(ieee_is_finite(a))) return

    do k = 1, n - 1
        pivot = k
        maxabs = abs(a(k, k))
        do i = k + 1, n
            if (abs(a(i, k)) > maxabs) then
                pivot = i
                maxabs = abs(a(i, k))
            end if
        end do
        pivot_scale = max(1._dl, maxval(abs(a(pivot, k:n))))
        if (maxabs <= 100._dl*epsilon(1._dl)*pivot_scale) return
        ipiv(k) = pivot
        if (pivot /= k) then
            do j = 1, n
                temp = a(k, j)
                a(k, j) = a(pivot, j)
                a(pivot, j) = temp
            end do
        end if
        do i = k + 1, n
            factor = a(i, k)/a(k, k)
            a(i, k) = factor
            do j = k + 1, n
                a(i, j) = a(i, j) - factor*a(k, j)
            end do
        end do
    end do
    ipiv(n) = n

    pivot_scale = max(1._dl, abs(a(n, n)))
    if (abs(a(n, n)) <= 100._dl*epsilon(1._dl)*pivot_scale) return
    ok = .true.

    end subroutine FactorSmallLinearSystem

    subroutine SolveFactoredSmallLinearSystem(a, ipiv, rhs, solution, ok)
    !Forward/back substitution using the output of FactorSmallLinearSystem
    real(dl), intent(in) :: a(:, :), rhs(:)
    integer, intent(in) :: ipiv(:)
    real(dl), intent(out) :: solution(:)
    logical, intent(out) :: ok
    real(dl) :: b(size(rhs)), temp
    integer :: i, j, k, n

    n = size(rhs)
    ok = .false.
    if (.not. all(ieee_is_finite(rhs))) return

    b = rhs
    do k = 1, n - 1
        if (ipiv(k) /= k) then
            temp = b(k)
            b(k) = b(ipiv(k))
            b(ipiv(k)) = temp
        end if
        do i = k + 1, n
            b(i) = b(i) - a(i, k)*b(k)
        end do
    end do

    solution(n) = b(n)/a(n, n)
    do i = n - 1, 1, -1
        solution(i) = b(i)
        do j = i + 1, n
            solution(i) = solution(i) - a(i, j)*solution(j)
        end do
        solution(i) = solution(i)/a(i, i)
    end do
    if (.not. all(ieee_is_finite(solution))) return
    ok = .true.

    end subroutine SolveFactoredSmallLinearSystem

    logical function RecfastRosenbrockStateOK(y)
    real(dl), intent(in) :: y(:)

    RecfastRosenbrockStateOK = all(ieee_is_finite(y))
    if (.not. RecfastRosenbrockStateOK) return
    if (size(y) >= 2) then
        RecfastRosenbrockStateOK = minval(y(1:2)) >= -1.e-8_dl
        if (.not. RecfastRosenbrockStateOK) return
    end if
    RecfastRosenbrockStateOK = y(3) > 0._dl

    end function RecfastRosenbrockStateOK

    subroutine EvaluateRecfastODETimeDerivative(this, Ndim, z, y, h, f0, f_t)
    !One-sided difference in z of the full smooth (force_full_hydrogen) system,
    !reusing the already-computed f0 = f(z, y)
    class(TRecfast), target :: this
    integer, intent(in) :: Ndim
    real(dl), intent(in) :: z, y(Ndim), h, f0(Ndim)
    real(dl), intent(out) :: f_t(Ndim)
    real(dl) :: delta_z, delta_z_floor, z_pert
    real(dl) :: f_pert(Ndim)

    delta_z = min(0.05_dl, max(1.e-4_dl, 1.e-5_dl*max(1._dl, abs(z))))
    if (abs(h) > 0._dl) then
        delta_z = min(delta_z, 0.5_dl*abs(h))
        delta_z_floor = 1.e-8_dl*max(1._dl, abs(z))
        delta_z = max(delta_z, delta_z_floor)
    end if
    !step in the direction of integration where possible
    z_pert = max(zfinal, min(zinitial, z + sign(delta_z, h)))
    if (z_pert == z) z_pert = max(zfinal, min(zinitial, z - sign(delta_z, h)))

    if (z_pert /= z) then
        call EvaluateRecfastODE(this, Ndim, z_pert, y, f_pert, .true.)
        f_t = (f_pert - f0)/(z_pert - z)
    else
        f_t = 0._dl
    end if

    end subroutine EvaluateRecfastODETimeDerivative

    subroutine RecfastROS2Step(this, z, y, h, f, f_t, jac, yout, yerr, ok)
    class(TRecfast), target :: this
    real(dl), intent(in) :: z, y(:), h, f(:), f_t(:), jac(:, :)
    real(dl), intent(out) :: yout(:)
    real(dl), intent(out) :: yerr(:)
    logical, intent(out) :: ok
    real(dl) :: f_stage(size(y)), matrix(size(y), size(y))
    real(dl) :: k1(size(y)), k2(size(y)), rhs(size(y)), y_stage(size(y))
    real(dl) :: gamma_h, gamma_h2
    integer :: i, n, ipiv(size(y))

    n = size(y)
    ok = .false.

    matrix = -ROS2_gamma*h*jac
    do i = 1, n
        matrix(i, i) = matrix(i, i) + 1._dl
    end do
    !factor once and reuse for both stage solves
    call FactorSmallLinearSystem(matrix, ipiv, ok)
    if (.not. ok) return
    gamma_h = ROS2_gamma*h
    gamma_h2 = gamma_h*gamma_h

    rhs = gamma_h*f + gamma_h2*f_t
    call SolveFactoredSmallLinearSystem(matrix, ipiv, rhs, k1, ok)
    if (.not. ok) return

    y_stage = y + ROS2_a21*k1
    call EvaluateRecfastODE(this, n, z + h, y_stage, f_stage, .true.)

    rhs = gamma_h*f_stage + ROS2_gamma*ROS2_c21*k1 - gamma_h2*f_t
    call SolveFactoredSmallLinearSystem(matrix, ipiv, rhs, k2, ok)
    if (.not. ok) return

    yout = y + ROS2_m1*k1 + ROS2_m2*k2
    yerr = (k1 + k2)/(2._dl*ROS2_gamma)
    ok = RecfastRosenbrockStateOK(yout) .and. all(ieee_is_finite(yerr))

    end subroutine RecfastROS2Step

    subroutine RecfastRosenbrockAdvance(this, zstart, y, zend, tol, ok)
    class(TRecfast), target :: this
    real(dl), intent(in) :: zstart, zend, tol
    real(dl), intent(inout) :: y(:)
    logical, intent(out) :: ok
    real(dl) :: direction, err, factor, min_step, scale, step, z
    real(dl) :: y_err(size(y)), y_trial(size(y))
    real(dl) :: f(size(y)), f_t(size(y)), jac(size(y), size(y))
    integer :: attempt, i
    logical :: step_ok, need_derivs

    ok = .false.
    if (zend == zstart) then
        ok = .true.
        return
    end if

    direction = sign(1._dl, zend - zstart)
    z = zstart
    step = zend - zstart
    min_step = max(abs(step)*1.e-8_dl, 1.e-10_dl)
    need_derivs = .true.

    do attempt = 1, RECFAST_rosenbrock_max_steps
        if (direction*(zend - z) <= 0._dl) then
            ok = .true.
            return
        end if
        if (direction*(zend - (z + step)) < 0._dl) step = zend - z

        if (need_derivs) then
            ! f, jacobian and f_t depend only on (z, y), so can be reused across
            ! step-size retries (f_t's stencil is sized from the first attempt's step, but
            ! it estimates an h-independent derivative on scales far below f's z-smoothness,
            ! so reuse at smaller h is fine). The Rosenbrock path always uses the full
            ! smooth hydrogen system rather than the old piecewise H switch.
            call EvaluateRecfastODE(this, size(y), z, y, f, .true., jac)
            call EvaluateRecfastODETimeDerivative(this, size(y), z, y, step, f, f_t)
            need_derivs = .false.
        end if
        call RecfastROS2Step(this, z, y, step, f, f_t, jac, y_trial, y_err, step_ok)

        if (step_ok) then
            err = 0._dl
            do i = 1, size(y)
                if (i <= 2) then
                    scale = max(RECFAST_rosenbrock_ion_scale_floor, abs(y(i)), abs(y_trial(i)))
                else
                    scale = max(1._dl, abs(y(i)), abs(y_trial(i)))
                end if
                err = max(err, abs(y_err(i))/scale)
            end do
            if (.not. ieee_is_finite(err)) step_ok = .false.
        end if

        if (step_ok .and. err <= tol) then
            y = y_trial
            y(1:min(2, size(y))) = max(y(1:min(2, size(y))), 0._dl)
            z = z + step
            need_derivs = .true.
            if (direction*(zend - z) <= 0._dl) then
                ok = .true.
                return
            end if
            if (err == 0._dl) then
                factor = ROS2_max_scale
            else
                factor = ROS2_safety*sqrt(tol/err)
                factor = min(ROS2_max_scale, max(ROS2_min_scale, factor))
            end if
            step = direction*min(abs(step)*factor, abs(zend - z))
        else
            if (step_ok) then
                factor = ROS2_safety*sqrt(tol/max(err, tiny(1._dl)))
                factor = min(0.9_dl, max(ROS2_min_scale, factor))
            else
                factor = ROS2_min_scale
            end if
            step = step*factor
            if (abs(step) < min_step) return
        end if
    end do

    end subroutine RecfastRosenbrockAdvance

    recursive subroutine EvaluateRecfastODE(this, Ndim, z, y, f, force_full_hydrogen, jacobian)
    class(TRecfast), target :: this
    integer, intent(in) :: Ndim
    real(dl), intent(in) :: z, y(Ndim)
    real(dl), intent(out) :: f(Ndim)
    logical, intent(in), optional :: force_full_hydrogen
    real(dl), intent(out), optional :: jacobian(Ndim, Ndim)
    real(dl) :: x, n, n_He, Trad, Tmat, Tspin, x_H, x_He, Hz, aTmat, ainv, aTs
    real(dl) :: Rup, Rdown, K, K_He, Rup_He, Rdown_He, He_Boltz
    real(dl) :: timeTh, timeH
    !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    real(dl), parameter :: a_PPB = 4.309d0, b_PPB = -0.6166d0, c_PPB = 0.6703d0, d_PPB = 0.5300d0
    !       the Verner and Ferland type fitting parameters for Helium
    !       fixed to match those in the SSS papers, and now correct
    real(dl), parameter :: a_VF = 10.d0**(-16.744d0), b_VF = 0.711d0
    real(dl), parameter :: T_0 = 10.d0**(0.477121d0), T_1 = 10.d0**(5.114d0)
    !      fitting parameters for HeI triplets
    !      (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)
    real(dl), parameter :: a_trip = 10.d0**(-16.306d0), b_trip = 0.761D0
    real(dl) :: ainv2, ainv3, inv_ainv, inv_Tmat, inv_Tmat2
    real(dl) :: Trad3, Trad4, theta, log_theta, theta_d, cr_tmat, cr_tmat_32
    real(dl) :: sq_0, sq_1, one_plus_sq_0, one_plus_sq_1, cpd, pgq, pgq_s, pgq_t, pgq_den, AH_fac
    real(dl) :: log_one_plus_sq_0, log_one_plus_sq_1, log_1pz, gauss_arg1, gauss_arg2
    real(dl) :: K_He2, C_H2, CHe2, x2, ion_sum, trip_denom, trip_denom2
    real(dl) :: tauHe_s, pHe_s, dpHe_s_dtau
    real(dl) :: Rdown_trip, Rup_trip
    real(dl) :: Doppler, gamma_2Ps, qb_s, AHcon
    !KIV (2007) H continuum opacity fit shapes; separate constants for the singlet
    !and triplet channels so the Jacobian cannot pick up stale values
    real(dl), parameter :: pb_s = 0.36d0, pb_t = 0.66d0, qb_t = 0.9d0
    real(dl) :: tauHe_t, pHe_t, dpHe_t_dtau, CL_PSt, CfHe_t, gamma_2Pt, AHcon_t
    real(dl) :: eps_tc, daTmat_dz, dTspin_dz
    real(dl) :: C10, dHdz, z_scale
    real(dl) :: He_rate_factor, He_rate_u, He_rate_u2
    real(dl) :: A_H, A_H_xH, A_H_xHe, A_H_T, B_H, C_H, C_H_xH, C_H_T
    real(dl) :: A_He, A_He_xH, A_He_xHe, A_He_T
    real(dl) :: A_trip_term, A_trip_xH, A_trip_xHe, A_trip_T
    real(dl) :: BHe, BHe_xH, BHe_xHe, BHe_T, CHe, CHe_xH, CHe_xHe, CHe_T
    real(dl) :: K_He_xH, K_He_xHe, K_He_T
    real(dl) :: AHcon_xH, AHcon_xHe, AHcon_dT
    real(dl) :: AHcon_t_xH, AHcon_t_xHe, AHcon_t_dT
    real(dl) :: dlnRdown, dRdown, dRup, dRupE
    real(dl) :: dlnRdown_He, dRdown_He, dRup_He, dHe_Boltz, dRupHeE
    real(dl) :: dlnRdown_trip, dRdown_trip, dRup_trip, dRupTripE, dEPSt
    real(dl) :: denH, EHe, EPSt, ETrip, LHe, LHe_xHe, LHe_T
    real(dl) :: MHe, MHe_xHe, MHe_T, pHe_s_xHe, pHe_t_xHe
    real(dl) :: RupE, RupHeE, RupTripE, S, S_T, eps_x, P, P_xH, P_xHe
    real(dl) :: Q, Q_xH, Q_xHe, Q_T, coupling_prefac, loose_prefac
    real(dl) :: Trip_source, Trip_source_xH, Trip_source_xHe, Trip_source_T
    real(dl) :: CfHe_t_xH, CfHe_t_xHe, CfHe_t_T, ypert_plus(Ndim), ypert_minus(Ndim), fpert_plus(Ndim), &
        fpert_minus(Ndim), delta
    real(dl) :: tauHe_s_const, tauHe_t_const
    integer :: col
    logical :: full_hydrogen, he_corrections
    type(RecombinationData), pointer :: Recomb

    Recomb => this%Calc
    full_hydrogen = .false.
    if (present(force_full_hydrogen)) full_hydrogen = force_full_hydrogen

    f = 0._dl
    if (present(jacobian)) jacobian = 0._dl

    x_H = y(1)
    x_He = y(2)
    he_corrections = x_He >= RECFAST_x_He_freeze_threshold .and. x_He <= 0.98d0
    x = x_H + Recomb%fHe*x_He
    ion_sum = 1.d0 + x + Recomb%fHe
    ainv = 1._dl + z
    inv_ainv = 1._dl/ainv
    ainv2 = ainv*ainv
    ainv3 = ainv2*ainv
    aTmat = y(3)
    Tmat = ainv*aTmat
    inv_Tmat = 1._dl/Tmat
    inv_Tmat2 = inv_Tmat*inv_Tmat

    n = Recomb%Nnow*ainv3
    n_He = Recomb%fHe*Recomb%Nnow*ainv3
    Trad = Recomb%Tnow*ainv
    Trad3 = Trad*Trad*Trad
    Trad4 = Trad3*Trad
    Hz = ainv2/dtauda(Recomb%State, inv_ainv)/MPC_in_sec
    denH = Hz*ainv
    cr_tmat = CR*Tmat
    cr_tmat_32 = cr_tmat*sqrt(cr_tmat)

    !       Get the radiative rates using PPQ fit, identical to Hummer's table
    theta = Tmat*1.d-4
    log_theta = log(theta)
    theta_d = exp(d_PPB*log_theta)
    cpd = c_PPB*theta_d
    Rdown = 1.d-19*a_PPB*exp(b_PPB*log_theta)/(1._dl + cpd)
    Rup = Rdown*cr_tmat_32*exp(-CDB*inv_Tmat)

    !       calculate He using a fit to a Verner & Ferland type formula
    sq_0 = sqrt(Tmat/T_0)
    sq_1 = sqrt(Tmat/T_1)
    one_plus_sq_0 = 1.d0 + sq_0
    one_plus_sq_1 = 1.d0 + sq_1
    log_one_plus_sq_0 = log(one_plus_sq_0)
    log_one_plus_sq_1 = log(one_plus_sq_1)
    !       typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF/sq_0*exp(-((1.d0 - b_VF)*log_one_plus_sq_0 &
        + (1.d0 + b_VF)*log_one_plus_sq_1))
    Rup_He = 4.d0*Rdown_He*cr_tmat_32*exp(-CDB_He*inv_Tmat)
    !       Avoid overflow (pointed out by Jacques Roland)
    if ((Bfact*inv_Tmat) > 680.d0) then
        He_Boltz = exp(680.d0)
    else
        He_Boltz = exp(Bfact*inv_Tmat)
    end if
    z_scale = this%Calc%Tnow/COBE_CMBTemp*ainv - 1

    !   now deal with H and its fudges
    if (.not. this%RECFAST_Hswitch) then
        K = CK/Hz
    else
        !  fit a double Gaussian correction function
        log_1pz = log(1.0d0 + z_scale)
        gauss_arg1 = (log_1pz - this%zGauss1)/this%wGauss1
        gauss_arg2 = (log_1pz - this%zGauss2)/this%wGauss2
        K = CK/Hz*(1.0d0 + this%AGauss1*exp(-gauss_arg1*gauss_arg1) &
            + this%AGauss2*exp(-gauss_arg2*gauss_arg2))
    end if

    !  add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip/sq_0*exp(-((1.d0 - b_trip)*log_one_plus_sq_0 &
        + (1.d0 + b_trip)*log_one_plus_sq_1))
    Rup_trip = Rdown_trip*exp(-h_P*C*L_He2St_ion/k_B*inv_Tmat)*cr_tmat_32*(4.d0/3.d0)
    !   last factor here is the statistical weight

    tauHe_s = 0._dl
    pHe_s = 1._dl
    dpHe_s_dtau = -0.5_dl
    K_He = CK_He/Hz
    AHcon = 0._dl
    AHcon_t = 0._dl
    CfHe_t = 0._dl
    CL_PSt = 0._dl
    EPSt = 0._dl
    pgq_s = 0._dl
    pgq_t = 0._dl
    He_rate_factor = 1._dl
    if (this%RECFAST_He_rate_correction .and. he_corrections .and. z_scale > RECFAST_He_rate_zmin &
        .and. z_scale < RECFAST_He_rate_zmax) then
        He_rate_u = (z_scale - RECFAST_He_rate_z0)/RECFAST_He_rate_width
        He_rate_u2 = He_rate_u*He_rate_u
        He_rate_factor = 1._dl + (RECFAST_He_rate_a0 + RECFAST_He_rate_a1*He_rate_u &
            + RECFAST_He_rate_a2*He_rate_u2)/(1._dl + RECFAST_He_rate_d2*He_rate_u2 &
            + RECFAST_He_rate_d4*He_rate_u2*He_rate_u2)
    end if

    ! Use Peebles coefficient for helium when He is frozen out or nearly fully ionized;
    ! otherwise apply the old RECFAST_Heswitch=6 helium corrections below.
    if (he_corrections) then
        ! Singlet 2^1P-1^1S Sobolev escape probability.
        tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0 - x_He)/Hz
        call EscapeProbabilityAndDerivative(tauHe_s, pHe_s, dpHe_s_dtau)
        K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0 - x_He))
        if (x_H < 0.9999999d0) then
            !AL changed July 08 to get smoother Helium

            ! H continuum opacity in the HeI singlet channel, using the fitting
            ! formula suggested by Kholupenko, Ivanchik & Varshalovich (2007).
            ! First get the Doppler width parameter.
            Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2p*sqrt(Doppler)
            gamma_2Ps = 3.d0*A2P_s*Recomb%fHe*(1.d0 - x_He)*C*C
            gamma_2Ps = gamma_2Ps/(sqrt(const_pi)*sigma_He_2Ps*const_eightpi*Doppler*(1.d0 - x_H))
            gamma_2Ps = gamma_2Ps/((C*L_He_2p)**2)
            qb_s = this%RECFAST_fudge_He
            !   calculate AHcon, the value of A*p_(con,H) for H continuum opacity
            pgq_s = pb_s*gamma_2Ps**qb_s
            AHcon = A2P_s/(1.d0 + pgq_s)
            K_He = 1.d0/((A2P_s*pHe_s + AHcon)*3.d0*n_He*(1.d0 - x_He))
        end if

        ! Triplet recombination through 2^3P-1^1S, with optional H continuum
        ! opacity in the triplet channel when neutral hydrogen is available.
        tauHe_t = A2P_t*n_He*(1.d0 - x_He)*3.d0
        tauHe_t = tauHe_t/(const_eightpi*Hz*L_He_2Pt**3)
        call EscapeProbabilityAndDerivative(tauHe_t, pHe_t, dpHe_t_dtau)
        CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
        EPSt = exp(-CL_PSt*inv_Tmat)
        if (x_H > 0.99999d0) then
            CfHe_t = A2P_t*pHe_t*EPSt
            CfHe_t = CfHe_t/(Rup_trip + CfHe_t)
        else
            Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2Pt*sqrt(Doppler)
            gamma_2Pt = 3.d0*A2P_t*Recomb%fHe*(1.d0 - x_He)*C*C
            gamma_2Pt = gamma_2Pt/(sqrt(const_pi)*sigma_He_2Pt*const_eightpi*Doppler*(1.d0 - x_H))
            gamma_2Pt = gamma_2Pt/((C*L_He_2Pt)**2)
            pgq_t = pb_t*gamma_2Pt**qb_t
            AHcon_t = A2P_t/(1.d0 + pgq_t)/3.d0
            CfHe_t = (A2P_t*pHe_t + AHcon_t)*EPSt
            CfHe_t = CfHe_t/(Rup_trip + CfHe_t)
        end if
    end if

    !       Estimates of Thomson scattering time and Hubble time
    timeTh = (1._dl/(CT*Trad4))*ion_sum/x
    timeH = 2._dl/(3._dl*Recomb%HO*ainv*sqrt(ainv))

    !       calculate the derivatives
    !       turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    !       (clunky, but seems to work)
    RupE = Rup*exp(-CL*inv_Tmat)
    A_H = x*x_H*n*Rdown - RupE*(1.d0 - x_H)
    if (.not. full_hydrogen) then
        if (x_H > 0.99d0) then
            f(1) = 0._dl
        else if (x_H > 0.985d0) then
            f(1) = A_H/denH
            Recomb%Recombination_saha_z = z
        else
            B_H = 1.d0 + K*Lambda*n*(1.d0 - x_H)
            C_H = 1.d0/Recomb%fu + K*Lambda*n*(1.d0 - x_H)/Recomb%fu + K*Rup*n*(1.d0 - x_H)
            f(1) = A_H*B_H/(denH*C_H)
        end if
    else
        B_H = 1.d0 + K*Lambda*n*(1.d0 - x_H)
        C_H = 1.d0/Recomb%fu + K*Lambda*n*(1.d0 - x_H)/Recomb%fu + K*Rup*n*(1.d0 - x_H)
        f(1) = A_H*B_H/(denH*C_H)
    end if

    !       turn off the He once it is small
    if (x_He < RECFAST_x_He_freeze_threshold) then
        f(2) = 0._dl
    else
        EHe = exp(-CL_He*inv_Tmat)
        RupHeE = Rup_He*EHe
        A_He = x*x_He*n*Rdown_He - RupHeE*(1.d0 - x_He)
        LHe = Lambda_He*n_He*(1.d0 - x_He)*He_Boltz
        MHe = (Lambda_He + Rup_He)*n_He*(1.d0 - x_He)*He_Boltz
        BHe = 1.d0 + K_He*LHe
        CHe = 1.d0 + K_He*MHe
        f(2) = A_He*BHe/(denH*CHe)

        if (he_corrections) then
            ETrip = exp(-h_P*C*L_He_2st/k_B*inv_Tmat)
            RupTripE = 3.d0*Rup_trip*ETrip
            A_trip_term = x*x_He*n*Rdown_trip - (1.d0 - x_He)*RupTripE
            f(2) = f(2) + A_trip_term*CfHe_t/denH
        end if
        f(2) = He_rate_factor*f(2)
    end if

    if (timeTh < H_frac*timeH) then
        ! Original RECFAST formula here is for dTmat/dz; written directly for aTmat.
        ! The first term is the exact tightly-coupled limit aTmat -> Tnow.
        dHdz = (Recomb%HO**2/2.d0/Hz)*(4.d0*ainv3/(1.d0 + Recomb%z_eq)*Recomb%OmegaT &
            + 3.d0*Recomb%OmegaT*ainv2 + 2.d0*Recomb%OmegaK*ainv)

        eps_tc = Hz*ion_sum/(CT*Trad3*x)
        daTmat_dz = (Recomb%Tnow - aTmat)*inv_ainv &
            + eps_tc*((1.d0 + Recomb%fHe)/ion_sum)*((f(1) + Recomb%fHe*f(2))/x)*inv_ainv &
            - eps_tc*dHdz/denH + 3.d0*eps_tc/ainv2
    else
        ! Original RECFAST formula is for dTmat/dz. Using Tmat=(1+z)*aTmat and Trad=(1+z)*Tnow gives:
        ! d(aTmat)/dz = [CT*Trad^4*x*(aTmat-Tnow)/(Hz*(1+x+fHe)) + aTmat]/(1+z).
        daTmat_dz = (CT*Trad4*x*(aTmat - Recomb%Tnow)/(Hz*ion_sum) + aTmat)*inv_ainv
    end if

    f(3) = daTmat_dz

    if (Evolve_Ts) then
        if (timeTh < H_frac*timeH) then
            f(4) = 0._dl
        else
            if (z < 1/Do21cm_minev - 1) then
                aTs = y(4)
                Tspin = ainv*aTs
                C10 = n*(kappa_HH_21cm(Tmat, .false.)*(1.d0 - x_H) + kappa_eH_21cm(Tmat, .false.)*x)
                dTspin_dz = 4*Tspin/(Hz*ainv)*((Tspin/Tmat - 1._dl)*C10 + Trad/T_21cm*(Tspin/Trad - 1._dl)*A10)
                dTspin_dz = dTspin_dz - f(1)*Tspin/(1.d0 - x_H)
                f(4) = (dTspin_dz - aTs)/ainv
            else
                f(4) = daTmat_dz
            end if
        end if
    end if

    if (.not. present(jacobian)) return

    dlnRdown = (b_PPB + cpd*(b_PPB - d_PPB))/(1.d0 + cpd)*inv_Tmat
    dRdown = Rdown*dlnRdown
    dRup = Rup*(dlnRdown + 1.5d0*inv_Tmat + CDB*inv_Tmat2)
    dRupE = RupE*(dlnRdown + 1.5d0*inv_Tmat + CB1*inv_Tmat2)

    A_H_xH = n*Rdown*(x + x_H) + RupE
    A_H_xHe = n*Rdown*Recomb%fHe*x_H
    A_H_T = n*x*x_H*dRdown - (1.d0 - x_H)*dRupE

    if (.not. full_hydrogen .and. x_H > 0.99d0) then
        jacobian(1, 1:3) = 0._dl
    else if (.not. full_hydrogen .and. x_H > 0.985d0) then
        jacobian(1, 1) = A_H_xH/denH
        jacobian(1, 2) = A_H_xHe/denH
        jacobian(1, 3) = ainv*A_H_T/denH
    else
        B_H = 1.d0 + K*Lambda*n*(1.d0 - x_H)
        C_H = 1.d0/Recomb%fu + K*Lambda*n*(1.d0 - x_H)/Recomb%fu + K*Rup*n*(1.d0 - x_H)
        C_H_xH = -K*n*(Lambda/Recomb%fu + Rup)
        C_H_T = K*n*(1.d0 - x_H)*dRup
        C_H2 = C_H*C_H
        jacobian(1, 1) = ((A_H_xH*B_H - A_H*K*Lambda*n)*C_H - A_H*B_H*C_H_xH)/(denH*C_H2)
        jacobian(1, 2) = A_H_xHe*B_H/(denH*C_H)
        jacobian(1, 3) = ainv*(A_H_T*B_H*C_H - A_H*B_H*C_H_T)/(denH*C_H2)
    end if

    dlnRdown_He = -(1.d0 + (1.d0 - b_VF)*sq_0/one_plus_sq_0 + (1.d0 + b_VF)*sq_1/one_plus_sq_1)
    dlnRdown_He = 0.5d0*dlnRdown_He*inv_Tmat
    dRdown_He = Rdown_He*dlnRdown_He
    dRup_He = Rup_He*(dlnRdown_He + 1.5d0*inv_Tmat + CDB_He*inv_Tmat2)
    if ((Bfact*inv_Tmat) > 680.d0) then
        dHe_Boltz = 0._dl
    else
        dHe_Boltz = -He_Boltz*Bfact*inv_Tmat2
    end if

    K_He_xH = 0._dl
    K_He_xHe = 0._dl
    K_He_T = 0._dl
    AHcon_xH = 0._dl
    AHcon_xHe = 0._dl
    AHcon_dT = 0._dl
    if (he_corrections) then
        tauHe_s_const = A2P_s*CK_He*3.d0*n_He/Hz
        pHe_s_xHe = dpHe_s_dtau*(-tauHe_s_const)
        if (x_H < 0.9999999d0) then
            pgq = pgq_s
            pgq_den = 1.d0 + pgq
            AH_fac = A2P_s*qb_s*pgq/(pgq_den*pgq_den)
            AHcon_xH = -AH_fac/(1.d0 - x_H)
            AHcon_xHe = AH_fac/(1.d0 - x_He)
            AHcon_dT = 0.5d0*AH_fac*inv_Tmat
        end if
        K_He2 = K_He*K_He
        K_He_xH = -K_He2*AHcon_xH*3.d0*n_He*(1.d0 - x_He)
        K_He_xHe = -K_He2*((A2P_s*pHe_s_xHe + AHcon_xHe)*3.d0*n_He*(1.d0 - x_He) &
            - (A2P_s*pHe_s + AHcon)*3.d0*n_He)
        K_He_T = -K_He2*AHcon_dT*3.d0*n_He*(1.d0 - x_He)
    end if

    if (x_He >= RECFAST_x_He_freeze_threshold) then
        dRupHeE = RupHeE*(dlnRdown_He + 1.5d0*inv_Tmat + CB1_He1*inv_Tmat2)
        A_He = x*x_He*n*Rdown_He - RupHeE*(1.d0 - x_He)
        A_He_xH = n*Rdown_He*x_He
        A_He_xHe = n*Rdown_He*(x + Recomb%fHe*x_He) + RupHeE
        A_He_T = n*x*x_He*dRdown_He - (1.d0 - x_He)*dRupHeE

        LHe = Lambda_He*n_He*(1.d0 - x_He)*He_Boltz
        LHe_xHe = -Lambda_He*n_He*He_Boltz
        LHe_T = Lambda_He*n_He*(1.d0 - x_He)*dHe_Boltz
        MHe = (Lambda_He + Rup_He)*n_He*(1.d0 - x_He)*He_Boltz
        MHe_xHe = -(Lambda_He + Rup_He)*n_He*He_Boltz
        MHe_T = n_He*(1.d0 - x_He)*(dRup_He*He_Boltz + (Lambda_He + Rup_He)*dHe_Boltz)

        BHe = 1.d0 + K_He*LHe
        CHe = 1.d0 + K_He*MHe
        BHe_xH = K_He_xH*LHe
        BHe_xHe = K_He_xHe*LHe + K_He*LHe_xHe
        BHe_T = K_He_T*LHe + K_He*LHe_T
        CHe_xH = K_He_xH*MHe
        CHe_xHe = K_He_xHe*MHe + K_He*MHe_xHe
        CHe_T = K_He_T*MHe + K_He*MHe_T

        CHe2 = CHe*CHe
        jacobian(2, 1) = ((A_He_xH*BHe + A_He*BHe_xH)*CHe - A_He*BHe*CHe_xH)/(denH*CHe2)
        jacobian(2, 2) = ((A_He_xHe*BHe + A_He*BHe_xHe)*CHe - A_He*BHe*CHe_xHe)/(denH*CHe2)
        jacobian(2, 3) = ainv*((A_He_T*BHe + A_He*BHe_T)*CHe - A_He*BHe*CHe_T)/(denH*CHe2)

        if (he_corrections) then
            dlnRdown_trip = -(1.d0 + (1.d0 - b_trip)*sq_0/one_plus_sq_0 &
                + (1.d0 + b_trip)*sq_1/one_plus_sq_1)
            dlnRdown_trip = 0.5d0*dlnRdown_trip*inv_Tmat
            dRdown_trip = Rdown_trip*dlnRdown_trip
            dRup_trip = Rup_trip*(dlnRdown_trip + 1.5d0*inv_Tmat + h_P*C*L_He2St_ion/k_B*inv_Tmat2)
            !exact exponent in RupTripE is h_P*C*(L_He2St_ion + L_He_2St)/k_B, which equals
            !CB1_He1 up to ~5e-6 relative (atomic-data rounding); good enough for the Jacobian
            dRupTripE = RupTripE*(dlnRdown_trip + 1.5d0*inv_Tmat + CB1_He1*inv_Tmat2)

            A_trip_term = x*x_He*n*Rdown_trip - (1.d0 - x_He)*RupTripE
            A_trip_xH = n*Rdown_trip*x_He
            A_trip_xHe = n*Rdown_trip*(x + Recomb%fHe*x_He) + RupTripE
            A_trip_T = n*x*x_He*dRdown_trip - (1.d0 - x_He)*dRupTripE

            AHcon_t_xH = 0._dl
            AHcon_t_xHe = 0._dl
            AHcon_t_dT = 0._dl
            pHe_t_xHe = dpHe_t_dtau*(-A2P_t*3.d0*n_He/(const_eightpi*Hz*L_He_2Pt**3))
            dEPSt = EPSt*CL_PSt*inv_Tmat2
            if (x_H <= 0.99999d0) then
                pgq = pgq_t
                pgq_den = 1.d0 + pgq
                AH_fac = A2P_t*qb_t*pgq/(3.d0*pgq_den*pgq_den)
                AHcon_t_xH = -AH_fac/(1.d0 - x_H)
                AHcon_t_xHe = AH_fac/(1.d0 - x_He)
                AHcon_t_dT = 0.5d0*AH_fac*inv_Tmat
            end if

            Trip_source = (A2P_t*pHe_t + AHcon_t)*EPSt
            Trip_source_xH = AHcon_t_xH*EPSt
            Trip_source_xHe = (A2P_t*pHe_t_xHe + AHcon_t_xHe)*EPSt
            Trip_source_T = AHcon_t_dT*EPSt + (A2P_t*pHe_t + AHcon_t)*dEPSt
            trip_denom = Rup_trip + Trip_source
            trip_denom2 = trip_denom*trip_denom
            CfHe_t_xH = Trip_source_xH*Rup_trip/trip_denom2
            CfHe_t_xHe = Trip_source_xHe*Rup_trip/trip_denom2
            CfHe_t_T = (Trip_source_T*Rup_trip - Trip_source*dRup_trip)/trip_denom2

            jacobian(2, 1) = jacobian(2, 1) + (A_trip_xH*CfHe_t + A_trip_term*CfHe_t_xH)/denH
            jacobian(2, 2) = jacobian(2, 2) + (A_trip_xHe*CfHe_t + A_trip_term*CfHe_t_xHe)/denH
            jacobian(2, 3) = jacobian(2, 3) + ainv*(A_trip_T*CfHe_t + A_trip_term*CfHe_t_T)/denH
        end if
        jacobian(2, 1:3) = He_rate_factor*jacobian(2, 1:3)
    else
        jacobian(2, 1:3) = 0._dl
    end if

    if (timeTh < H_frac*timeH) then
        S = f(1) + Recomb%fHe*f(2)
        S_T = (jacobian(1, 3) + Recomb%fHe*jacobian(2, 3))*inv_ainv
        x2 = x*x
        eps_x = -Hz*(1.d0 + Recomb%fHe)/(CT*Trad3*x2)
        P = (1.d0 + Recomb%fHe)/ion_sum
        P_xH = -P/ion_sum
        P_xHe = Recomb%fHe*P_xH
        Q = S/x
        Q_xH = ((jacobian(1, 1) + Recomb%fHe*jacobian(2, 1))*x - S)/x2
        Q_xHe = ((jacobian(1, 2) + Recomb%fHe*jacobian(2, 2))*x - S*Recomb%fHe)/x2
        Q_T = S_T/x
        coupling_prefac = -dHdz/(Hz*ainv) + 3.d0/ainv2

        jacobian(3, 1) = (eps_x*P*Q + eps_tc*P_xH*Q + eps_tc*P*Q_xH)*inv_ainv + eps_x*coupling_prefac
        jacobian(3, 2) = (Recomb%fHe*eps_x*P*Q + eps_tc*P_xHe*Q + eps_tc*P*Q_xHe)*inv_ainv
        jacobian(3, 2) = jacobian(3, 2) + Recomb%fHe*eps_x*coupling_prefac
        jacobian(3, 3) = -inv_ainv + eps_tc*P*Q_T
    else
        loose_prefac = CT*Trad4*(aTmat - Recomb%Tnow)/(Hz*ion_sum*ion_sum)
        jacobian(3, 1) = loose_prefac*(1.d0 + Recomb%fHe)*inv_ainv
        jacobian(3, 2) = loose_prefac*(1.d0 + Recomb%fHe)*Recomb%fHe*inv_ainv
        jacobian(3, 3) = (CT*Trad4*x/(Hz*ion_sum) + 1.d0)*inv_ainv
    end if

    if (Ndim > 3) then
        jacobian(1:3, 4) = 0._dl
        do col = 1, Ndim
            delta = 1.e-7_dl*max(1._dl, abs(y(col)))
            ypert_plus = y
            ypert_plus(col) = ypert_plus(col) + delta
            if (col <= 2 .and. y(col) - delta < 0._dl) then
                call EvaluateRecfastODE(this, Ndim, z, ypert_plus, fpert_plus, full_hydrogen)
                jacobian(4, col) = (fpert_plus(4) - f(4))/delta
            else
                ypert_minus = y
                ypert_minus(col) = ypert_minus(col) - delta
                call EvaluateRecfastODE(this, Ndim, z, ypert_plus, fpert_plus, full_hydrogen)
                call EvaluateRecfastODE(this, Ndim, z, ypert_minus, fpert_minus, full_hydrogen)
                jacobian(4, col) = (fpert_plus(4) - fpert_minus(4))/(2._dl*delta)
            end if
        end do
    end if

    end subroutine EvaluateRecfastODE

    subroutine ION(this,Ndim,z,Y,f)
    class(TRecfast), target :: this
    integer Ndim
    real(dl) :: z
    real(dl) :: y(Ndim), f(Ndim)

    call EvaluateRecfastODE(this, Ndim, z, y, f)

    end subroutine ION


    function TRecfast_dDeltaxe_dtau(this,a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa)
    !d x_e/d tau assuming Helium all neutral and temperature perturbations negligible
    !it is not accurate for x_e of order 1
    class(TRecfast) :: this
    real(dl) TRecfast_dDeltaxe_dtau
    real(dl), intent(in):: a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa
    real(dl) Delta_Tg
    real(dl) xedot,z,x,n,n_He,Trad,Tmat,x_H,Hz, C_r, dlnC_r
    real(dl) Rup,Rdown,K
    !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    real(dl), parameter :: a_PPB = 4.309d0, b_PPB = -0.6166d0, c_PPB = 0.6703d0, d_PPB = 0.5300d0
    real(dl) delta_alpha, delta_beta, delta_K, clh
    real(dl) xe

    associate(Calc=>this%Calc)
        Delta_tg =Delta_Tm
        call this%xe_Tm(a, xe, Tmat)
        x_H = min(1._dl,xe)

        z=1/a-1

        x = x_H

        n = Calc%Nnow /a**3
        n_He = Calc%fHe * n
        Trad = Calc%Tnow /a
        clh = adotoa !conformal time Hubble
        Hz = clh/a/MPC_in_sec !normal time in seconds

        !       Get the radiative rates using PPQ fit, identical to Hummer's table

        Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
            /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)   !alpha
        Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)

        K = CK/Hz              !Peebles coefficient K=lambda_a^3/8piH


        Rdown = Rdown*Calc%fu
        Rup = Rup*Calc%fu
        C_r =  a*(1.d0 + K*Lambda*n*(1.d0-x_H)) /( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) )*MPC_in_sec

        xedot = -(x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat))*C_r

        delta_alpha = (b_PPB + c_PPB*(Tmat/1d4)**d_PPB*(b_PPB-d_PPB))/(1+c_PPB*(Tmat/1d4)**d_PPB)*Delta_Tg
        delta_beta = delta_alpha + (3./2 + CDB/Tmat)*delta_Tg !(Rup = beta)
        delta_K = - hdot/clh - kvb/clh/3


        dlnC_r = -Rup*K*n*( (Delta_nH+Delta_K + Delta_beta*(1+K*Lambda*n*(1-x_H)))*(1-x_H) - x_H*Delta_xe) &
            / ( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) ) /(1.d0 + K*Lambda*n*(1.d0-x_H))

        TRecfast_dDeltaxe_dtau= xedot/x_H*(dlnC_r +Delta_alpha - Delta_xe) &
            - C_r*( (2*Delta_xe + Delta_nH)*x_H*n*Rdown + (Delta_xe - (3./2+ CB1/Tmat)*(1/x_H-1)*Delta_Tg)*Rup*exp(-CL/Tmat))
    end associate

    !Approximate form valid at late times
    !        dDeltaxe_dtau= xedot/x_H*(Delta_alpha + Delta_xe + Delta_nH)


    end function TRecfast_dDeltaxe_dtau

    real(dl) function TRecfast_Get_Saha_z(this)
    class(TRecfast) :: this
    TRecfast_Get_Saha_z =  this%Calc%recombination_saha_z
    end function


    subroutine TRecfast_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TRecfast), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TRecfast_SelfPointer

    end module Recombination
