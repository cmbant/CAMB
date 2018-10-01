    ! Modules used by cmbmain and other routines.

    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.
    !
    !     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
    !     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
    !     Original CMBFAST copyright and disclaimer:
    !
    !     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
    !     the Massachusetts Institute of Technology.  All rights reserved.
    !
    !     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
    !     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
    !     By way of example, but not limitation,
    !     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
    !     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
    !     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
    !     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !     portions of this software are based on the COSMICS package of
    !     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
    !     for restrictions on the modification and distribution of this software.


    module ModelParams
    use precision
    use Ranges
    use InitialPower
    use Reionization
    use Recombination
    use Errors

    implicit none
    public

    character(LEN=*), parameter :: version = 'Aug18'

    integer :: FeedbackLevel = 0 !if >0 print out useful information about the model

    logical :: output_file_headers = .true.


    logical, parameter :: DebugMsgs=.false. !Set to true to view progress and timing

    logical, parameter :: DebugEvolution = .false. !Set to true to do all the evolution for all k

    real(dl) :: DebugParam = 0._dl !not used but read in, useful for parameter-dependent tests

    logical ::  do_bispectrum  = .false.
    logical, parameter :: hard_bispectrum = .false. ! e.g. warm inflation where delicate cancellations

    logical, parameter :: full_bessel_integration = .false. !(go into the tails when calculating the sources)

    integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
    !For CAMBparams%MassiveNuMethod
    !Nu_int: always integrate distribution function
    !Nu_trunc: switch to expansion in velocity once non-relativistic
    !Nu_approx: approximate scheme - good for CMB, but not formally correct and no good for matter power
    !Nu_best: automatically use mixture which is fastest and most accurate

    integer, parameter :: max_Nu = 5 !Maximum number of neutrino species
    integer, parameter :: max_transfer_redshifts = 150
    integer, parameter :: fileio_unit = 13 !Any number not used elsewhere will do
    integer, parameter :: outNone=1

    integer :: max_bessels_l_index  = 1000000
    real(dl) :: max_bessels_etak = 1000000*2

    real(dl), parameter ::  OutputDenominator =twopi
    !When using outNone the output is l(l+1)Cl/OutputDenominator

    Type(Regions) :: TimeSteps

    type TransferParams
        logical     ::  high_precision
        logical     ::  accurate_massive_neutrinos
        integer     ::  num_redshifts
        real(dl)    ::  kmax         !these are acutally q values, but same as k for flat
        integer     ::  k_per_logint ! ..
        real(dl)    ::  redshifts(max_transfer_redshifts)
        !JD 08/13 Added so both NL lensing and PK can be run at the same time
        real(dl)    ::  PK_redshifts(max_transfer_redshifts)
        real(dl)    ::  NLL_redshifts(max_transfer_redshifts)
        integer     ::  PK_redshifts_index(max_transfer_redshifts)
        integer     ::  NLL_redshifts_index(max_transfer_redshifts)
        integer     ::  PK_num_redshifts
        integer     ::  NLL_num_redshifts

    end type TransferParams

    !other variables, options, derived variables, etc.

    integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2
    integer, parameter :: NonLinear_both=3  !JD 08/13 added so both can be done

    ! Main parameters type
    type CAMBparams

        logical   :: WantCls, WantTransfer
        logical   :: WantScalars, WantTensors, WantVectors
        logical   :: DoLensing
        logical   :: want_zstar, want_zdrag     !!JH for updated BAO likelihood.
        logical   :: PK_WantTransfer             !JD 08/13 Added so both NL lensing and PK can be run at the same time
        integer   :: NonLinear
        logical   :: Want_CMB

        integer   :: Max_l, Max_l_tensor
        real(dl)  :: Max_eta_k, Max_eta_k_tensor
        ! _tensor settings only used in initialization,
        !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

        real(dl)  :: omegab, omegac, omegav, omegan
        !Omega baryon, CDM, Lambda and massive neutrino
        real(dl)  :: H0,TCMB,yhe,Num_Nu_massless
        integer   :: Num_Nu_massive !sum of Nu_mass_numbers below
        integer   :: Nu_mass_eigenstates  !1 for degenerate masses
        logical   :: share_delta_neff !take fractional part to heat all eigenstates the same
        real(dl)  :: Nu_mass_degeneracies(max_nu)
        real(dl)  :: Nu_mass_fractions(max_nu) !The ratios of the total densities
        integer   :: Nu_mass_numbers(max_nu) !physical number per eigenstate

        integer   :: Scalar_initial_condition
        !must be one of the initial_xxx values defined in GaugeInterface

        integer   :: OutputNormalization
        !outNone, or C_OutputNormalization=1 if > 1

        logical   :: AccuratePolarization
        !Do you care about the accuracy of the polarization Cls?

        logical   :: AccurateBB
        !Do you care about BB accuracy (e.g. in lensing)

        !Reionization settings - used if Reion%Reionization=.true.
        logical   :: AccurateReionization
        !Do you care about pecent level accuracy on EE signal from reionization?

        integer   :: MassiveNuMethod

        type(InitialPowerParams) :: InitPower  !see power_tilt.f90 - you can change this
        type(ReionizationParams) :: Reion
        type(RecombinationParams):: Recomb
        type(TransferParams)     :: Transfer

        real(dl) ::  InitialConditionVector(1:10) !Allow up to 10 for future extensions
        !ignored unless Scalar_initial_condition == initial_vector

        logical OnlyTransfers !Don't use initial power spectrum data, instead get Delta_q_l array
        !If true, sigma_8 is not calculated either

        logical DerivedParameters !calculate various derived parameters  (ThermoDerivedParams)

        !Derived parameters, not set initially
        type(ReionizationHistory) :: ReionHist

        logical flat,closed,open
        real(dl) omegak
        real(dl) curv,r, Ksign !CP%r = 1/sqrt(|CP%curv|), CP%Ksign = 1,0 or -1
        real(dl) tau0,chi0 !time today and rofChi(CP%tau0/CP%r)

    end type CAMBparams

    type(CAMBparams), save :: CP  !Global collection of parameters

    real(dl) scale !relative to CP%flat. e.g. for scaling lSamp%l sampling.

    logical ::call_again = .false.
    !if being called again with same parameters to get different thing

    !     grhom =kappa*a^2*rho_m0
    !     grhornomass=grhor*number of massless neutrino species
    !     taurst,taurend - time at start/end of recombination
    !     dtaurec - dtau during recombination
    !     adotrad - a(tau) in radiation era

    real(dl) grhom,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
    real(dl) taurst,dtaurec,taurend, tau_maxvis,adotrad

    !Neutrinos
    real(dl) grhormass(max_nu)

    !     nu_masses=m_nu*c**2/(k_B*T_nu0)
    real(dl) :: nu_masses(max_nu)

    real(dl) akthom !sigma_T * (number density of protons now)
    real(dl) fHe !n_He_tot / n_H_tot
    real(dl) Nnow


    integer :: ThreadNum = 0
    !If zero assigned automatically, obviously only used if parallelised

    !Parameters for checking/changing overall accuracy
    !If HighAccuracyDefault=.false., the other parameters equal to 1 corresponds to ~0.3% scalar C_l accuracy
    !If HighAccuracyDefault=.true., the other parameters equal to 1 corresponds to ~0.1% scalar C_l accuracy (at L>600)
    logical :: HighAccuracyDefault = .true.

    real(dl) :: lSampleBoost=1._dl
    !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation

    real(dl) :: AccuracyBoost =1._dl

    !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
    !Can also be used to improve speed significantly if less accuracy is required.
    !or improving accuracy for extreme models.
    !Note this does not increase lSamp%l sampling or massive neutrino q-sampling

    real(sp) :: lAccuracyBoost=1.
    !Boost number of multipoles integrated in Boltzman heirarchy

    integer :: limber_phiphi = 0 !for l>limber_phiphi use limber approx for lensing potential
    integer :: num_redshiftwindows = 0
    integer :: num_extra_redshiftwindows = 0
    integer :: num_custom_sources = 0
    integer, allocatable :: custom_source_ell_scales(:)

    integer, parameter :: lmin = 2
    !must be either 1 or 2

    real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code

    real(dl),parameter :: tol=1.0d-4 !Base tolerance for integrations

    !     used as parameter for spline - tells it to use 'natural' end values
    real(dl), parameter :: spl_large=1.e40_dl

    integer, parameter:: l0max=4000

    !     lmax is max possible number of l's evaluated
    integer, parameter :: lmax_arr = l0max

    character(LEN=1024) :: highL_unlensed_cl_template = 'HighLExtrapTemplate_lenspotentialCls.dat'
    !fiducial high-accuracy high-L C_L used for making small cosmology-independent numerical corrections
    !to lensing and C_L interpolation. Ideally close to models of interest, but dependence is weak.
    logical :: use_spline_template = .true.
    integer, parameter :: lmax_extrap_highl = 8000
    real(dl), allocatable :: highL_CL_template(:,:)

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13
    integer, parameter :: nthermo_derived = 13

    real(dl) ThermoDerivedParams(nthermo_derived)

    Type TBackgroundOutputs
        real(dl), pointer :: z_outputs(:) => null()
        real(dl), allocatable :: H(:), DA(:), rs_by_D_v(:)
    end Type TBackgroundOutputs

    Type(TBackgroundOutputs), save :: BackgroundOutputs

    contains


    subroutine CAMBParams_Set(P, error, DoReion)
    use constants
    type(CAMBparams), intent(in) :: P
    real(dl) GetOmegak, fractional_number, conv
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical WantReion
    integer nu_i,actual_massless
    real(dl) nu_massless_degeneracy, neff_i
    external GetOmegak
    real(dl), save :: last_tau0
    !Constants in SI units

    global_error_flag = 0

    if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
        call GlobalError( 'Cannot generate tensor C_l and transfer without scalar C_l',error_unsupported_params)
    end if

    if (present(error)) error = global_error_flag
    if (global_error_flag/=0) return

    if (present(DoReion)) then
        WantReion = DoReion
    else
        WantReion = .true.
    end if

    CP=P
    if (call_again) CP%DerivedParameters = .false.

    CP%Max_eta_k = max(CP%Max_eta_k,CP%Max_eta_k_tensor)

    if (CP%WantTransfer) then
        CP%WantScalars=.true.
        if (.not. CP%WantCls) then
            CP%AccuratePolarization = .false.
            CP%Reion%Reionization = .false.
        end if
    else
        CP%transfer%num_redshifts=0
    end if

    if (CP%Num_Nu_Massive /= sum(CP%Nu_mass_numbers(1:CP%Nu_mass_eigenstates))) then
        if (sum(CP%Nu_mass_numbers(1:CP%Nu_mass_eigenstates))/=0) call MpiStop('Num_Nu_Massive is not sum of Nu_mass_numbers')
    end if
    if (CP%Omegan == 0 .and. CP%Num_Nu_Massive /=0) then
        if (CP%share_delta_neff) then
            CP%Num_Nu_Massless = CP%Num_Nu_Massless + CP%Num_Nu_Massive
        else
            CP%Num_Nu_Massless = CP%Num_Nu_Massless + sum(CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates))
        end if
        CP%Num_Nu_Massive  = 0
        CP%Nu_mass_numbers = 0
    end if

    nu_massless_degeneracy = CP%Num_Nu_massless !N_eff for massless neutrinos
    if (CP%Num_nu_massive > 0) then
        if (CP%Nu_mass_eigenstates==0) call MpiStop('Have Num_nu_massive>0 but no nu_mass_eigenstates')
        if (CP%Nu_mass_eigenstates==1 .and. CP%Nu_mass_numbers(1)==0) CP%Nu_mass_numbers(1) = CP%Num_Nu_Massive
        if (all(CP%Nu_mass_numbers(1:CP%Nu_mass_eigenstates)==0)) CP%Nu_mass_numbers=1 !just assume one for all
        if (CP%share_delta_neff) then
            !default case of equal heating of all neutrinos
            fractional_number = CP%Num_Nu_massless + CP%Num_Nu_massive
            actual_massless = int(CP%Num_Nu_massless + 1e-6_dl)
            neff_i = fractional_number/(actual_massless + CP%Num_Nu_massive)
            nu_massless_degeneracy = neff_i*actual_massless
            CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates) = CP%Nu_mass_numbers(1:CP%Nu_mass_eigenstates)*neff_i
        end if
        if (abs(sum(CP%Nu_mass_fractions(1:CP%Nu_mass_eigenstates))-1) > 1e-4) &
            call MpiStop('Nu_mass_fractions do not add up to 1')
    else
        CP%Nu_mass_eigenstates = 0
    end if

    if ((CP%WantTransfer).and. CP%MassiveNuMethod==Nu_approx) then
        CP%MassiveNuMethod = Nu_trunc
    end if

    CP%omegak = GetOmegak()

    CP%flat = (abs(CP%omegak) <= OmegaKFlat)
    CP%closed = CP%omegak < -OmegaKFlat

    CP%open = .not.CP%flat.and..not.CP%closed
    if (CP%flat) then
        CP%curv=0
        CP%Ksign=0
        CP%r=1._dl !so we can use tau/CP%r, etc, where CP%r's cancel
    else
        CP%curv=-CP%omegak/((c/1000)/CP%h0)**2
        CP%Ksign =sign(1._dl,CP%curv)
        CP%r=1._dl/sqrt(abs(CP%curv))
    end if
    !  grho gives the contribution to the expansion rate from: (g) photons,
    !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
    !  (m) nonrelativistic matter (for Omega=1).  grho is actually
    !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
    !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
    !  (Used only to set the initial conformal time.)

    !H0 is in km/s/Mpc

    grhom = 3*CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

    !grhom=3.3379d-11*h0*h0
    grhog = kappa/c**2*4*sigma_boltz/c**3*CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
    ! grhog=1.4952d-13*tcmb**4
    grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
    !grhor=3.3957d-14*tcmb**4

    !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
    !for massive Nu_mass_degeneracies parameters account for heating from grhor

    grhornomass=grhor*nu_massless_degeneracy
    grhormass=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass(nu_i)=grhor*CP%Nu_mass_degeneracies(nu_i)
    end do
    grhoc=grhom*CP%omegac
    grhob=grhom*CP%omegab
    grhov=grhom*CP%omegav
    grhok=grhom*CP%omegak
    !  adotrad gives the relation a(tau) in the radiation era:
    adotrad = sqrt((grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates)))/3)


    Nnow = CP%omegab*(1-CP%yhe)*grhom*c**2/kappa/m_H/Mpc**2

    akthom = sigma_thomson*Nnow*Mpc
    !sigma_T * (number density of protons now)

    fHe = CP%YHe/(mass_ratio_He_H*(1.d0-CP%YHe))  !n_He_tot / n_H_tot

    if (.not.call_again) then
        call init_massive_nu(CP%omegan /=0)
        call init_background
        if (global_error_flag==0) then
            CP%tau0=TimeOfz(0._dl)
            ! print *, 'chi = ',  (CP%tau0 - TimeOfz(0.15_dl)) * CP%h0/100
            last_tau0=CP%tau0
            if (WantReion) call Reionization_Init(CP%Reion,CP%ReionHist, CP%YHe, akthom, CP%tau0, FeedbackLevel)
        end if
    else
        CP%tau0=last_tau0
    end if

    !JD 08/13 Changes for nonlinear lensing of CMB + MPK compatibility
    !if ( CP%NonLinear==NonLinear_Lens) then
    if (CP%NonLinear==NonLinear_Lens .or. CP%NonLinear==NonLinear_both ) then
        CP%Transfer%kmax = max(CP%Transfer%kmax, CP%Max_eta_k/CP%tau0)
        if (FeedbackLevel > 0 .and. CP%Transfer%kmax== CP%Max_eta_k/CP%tau0) &
            write (*,*) 'max_eta_k changed to ', CP%Max_eta_k
    end if

    if (CP%closed .and. CP%tau0/CP%r >3.14) then
        call GlobalError('chi >= pi in closed model not supported',error_unsupported_params)
    end if

    if (global_error_flag/=0) then
        if (present(error)) error = global_error_flag
        return
    end if

    if (present(error)) then
        error = 0
    else if (FeedbackLevel > 0 .and. .not. call_again) then
        write(*,'("Om_b h^2             = ",f9.6)') CP%omegab*(CP%H0/100)**2
        write(*,'("Om_c h^2             = ",f9.6)') CP%omegac*(CP%H0/100)**2
        write(*,'("Om_nu h^2            = ",f9.6)') CP%omegan*(CP%H0/100)**2
        write(*,'("Om_Lambda            = ",f9.6)') CP%omegav
        write(*,'("Om_K                 = ",f9.6)') CP%omegak
        write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-CP%omegak-CP%omegav
        write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta()
        if (CP%Num_Nu_Massive > 0) then
            write(*,'("N_eff (total)        = ",f9.6)') nu_massless_degeneracy + &
                sum(CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates))
            do nu_i=1, CP%Nu_mass_eigenstates
                conv = k_B*(8*grhor/grhog/7)**0.25*CP%tcmb/eV * &
                    (CP%nu_mass_degeneracies(nu_i)/CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                write(*,'(I2, " nu, g=",f7.4," m_nu*c^2/k_B/T_nu0= ",f9.2," (m_nu= ",f6.3," eV)")') &
                    CP%nu_mass_numbers(nu_i), CP%nu_mass_degeneracies(nu_i), nu_masses(nu_i),conv*nu_masses(nu_i)
            end do
        end if
    end if
    CP%chi0=rofChi(CP%tau0/CP%r)
    scale= CP%chi0*CP%r/CP%tau0  !e.g. change l sampling depending on approx peak spacing

    end subroutine CAMBParams_Set


    function GetTestTime()
    real(sp) GetTestTime
    real(sp) atime

    !           GetTestTime = etime(tarray)
    !Can replace this if etime gives problems
    !Or just comment out - only used if DebugMsgs = .true.
    call cpu_time(atime)
    GetTestTime = atime

    end function GetTestTime


    function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
    real(dl) Chi,rofChi

    if (CP%closed) then
        rofChi=sin(chi)
    else if (CP%open) then
        rofChi=sinh(chi)
    else
        rofChi=chi
    endif
    end function rofChi


    function cosfunc (Chi)
    real(dl) Chi,cosfunc

    if (CP%closed) then
        cosfunc= cos(chi)
    else if (CP%open) then
        cosfunc=cosh(chi)
    else
        cosfunc = 1._dl
    endif
    end function cosfunc

    function tanfunc(Chi)
    real(dl) Chi,tanfunc
    if (CP%closed) then
        tanfunc=tan(Chi)
    else if (CP%open) then
        tanfunc=tanh(Chi)
    else
        tanfunc=Chi
    end if

    end  function tanfunc

    function invsinfunc(x)
    real(dl) invsinfunc,x

    if (CP%closed) then
        invsinfunc=asin(x)
    else if (CP%open) then
        invsinfunc=log((x+sqrt(1._dl+x**2)))
    else
        invsinfunc = x
    endif
    end function invsinfunc

    function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x
    f_K = CP%r*rofChi(x/CP%r)

    end function f_K


    function DeltaTime(a1,a2, in_tol)
    implicit none
    real(dl) DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) dtauda, rombint !diff of tau w.r.t a and integration
    external dtauda, rombint

    if (present(in_tol)) then
        atol = in_tol
    else
        atol = tol/1000/exp(AccuracyBoost-1)
    end if
    DeltaTime=rombint(dtauda,a1,a2,atol)

    end function DeltaTime

    function TimeOfz(z)
    implicit none
    real(dl) TimeOfz
    real(dl), intent(IN) :: z

    TimeOfz=DeltaTime(0._dl,1._dl/(z+1._dl))
    end function TimeOfz

    subroutine TimeOfzArr(nz, redshifts, outputs)
    integer, intent(in) :: nz
    real(dl), intent(in) :: redshifts(nz)
    real(dl), intent(out) :: outputs(nz)
    integer i

    !Dumb slow version
    !$OMP PARALLEL DO DEFAUlT(SHARED)
    do i=1, nz
        outputs(i) = timeOfZ(redshifts(i))
    end do
    !$OMP END PARALLEL DO

    end subroutine TimeOfzArr

    function DeltaPhysicalTimeGyr(a1,a2, in_tol)
    use constants
    real(dl), intent(in) :: a1, a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) rombint,DeltaPhysicalTimeGyr, atol
    external rombint

    if (present(in_tol)) then
        atol = in_tol
    else
        atol = 1d-4/exp(AccuracyBoost-1)
    end if
    DeltaPhysicalTimeGyr = rombint(dtda,a1,a2,atol)*Mpc/c/Gyr
    end function DeltaPhysicalTimeGyr

    function AngularDiameterDistance(z)
    !This is the physical (non-comoving) angular diameter distance in Mpc
    real(dl) AngularDiameterDistance
    real(dl), intent(in) :: z

    AngularDiameterDistance = CP%r/(1+z)*rofchi(ComovingRadialDistance(z) /CP%r)

    end function AngularDiameterDistance

    subroutine AngularDiameterDistanceArr(arr, z, n)
    !This is the physical (non-comoving) angular diameter distance in Mpc for array of z
    !z array must be monotonically increasing
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    call ComovingRadialDistanceArr(arr, z, n, 1e-4_dl)
    if (CP%flat) then
        arr = arr/(1+z)
    else
        do i=1, n
            arr(i) =  CP%r/(1+z(i))*rofchi(arr(i)/CP%r)
        end do
    end if

    end subroutine AngularDiameterDistanceArr


    function AngularDiameterDistance2(z1, z2) ! z1 < z2
    !From http://www.slac.stanford.edu/~amantz/work/fgas14/#cosmomc
    real(dl) AngularDiameterDistance2
    real(dl), intent(in) :: z1, z2

    AngularDiameterDistance2 = CP%r/(1+z2)*rofchi(ComovingRadialDistance(z2)/CP%r - ComovingRadialDistance(z1)/CP%r)

    end function AngularDiameterDistance2

    function LuminosityDistance(z)
    real(dl) LuminosityDistance
    real(dl), intent(in) :: z

    LuminosityDistance = AngularDiameterDistance(z)*(1+z)**2

    end function LuminosityDistance

    function ComovingRadialDistance(z)
    real(dl) ComovingRadialDistance
    real(dl), intent(in) :: z

    ComovingRadialDistance = DeltaTime(1/(1+z),1._dl)

    end function ComovingRadialDistance

    subroutine ComovingRadialDistanceArr(arr, z, n, tol)
    !z array must be monotonically increasing
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    real(dl), intent(in) :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        if (i==1) then
            if (z(i) < 1e-6_dl) then
                arr(i) = 0
            else
                arr(i) = DeltaTime(1/(1+z(i)),1._dl, tol)
            end if
        else
            if (z(i) < z(i-1)) error stop 'ComovingRadialDistanceArr redshifts out of order'
            arr(i) = DeltaTime(1/(1+z(i)),1/(1+z(i-1)),tol)
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine ComovingRadialDistanceArr

    function Hofz(z)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    real(dl) Hofz, dtauda,a
    real(dl), intent(in) :: z
    external dtauda

    a = 1/(1+z)
    Hofz = 1/(a**2*dtauda(a))

    end function Hofz

    subroutine HofzArr(arr, z, n)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    do i=1, n
        arr(i) = Hofz(z(i))
    end do

    end subroutine HofzArr

    real(dl) function BAO_D_v_from_DA_H(z, DA, Hz)
    real(dl), intent(in) :: z, DA, Hz
    real(dl) ADD

    ADD = DA*(1.d0+z)
    BAO_D_v_from_DA_H = ((ADD)**2.d0*z/Hz)**(1.d0/3.d0)

    end function BAO_D_v_from_DA_H

    real(dl) function BAO_D_v(z)
    real(dl), intent(IN) :: z

    BAO_D_v = BAO_D_v_from_DA_H(z,AngularDiameterDistance(z), Hofz(z))

    end function BAO_D_v

    function dsound_da_exact(a)
    implicit none
    real(dl) dsound_da_exact,dtauda,a,R,cs
    external dtauda

    R = 3*grhob*a / (4*grhog)
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_exact=dtauda(a)*cs

    end function dsound_da_exact


    function dsound_da(a)
    !approximate form used e.g. by CosmoMC for theta
    implicit none
    real(dl) dsound_da,dtauda,a,R,cs
    external dtauda

    R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
    !          R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and used for previous consistency
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da=dtauda(a)*cs

    end function dsound_da

    function dtda(a)
    real(dl) dtda,dtauda,a
    external dtauda
    dtda= dtauda(a)*a
    end function

    function CosmomcTheta()
    real(dl) zstar, astar, atol, rs, DA
    real(dl) CosmomcTheta
    real(dl) ombh2, omdmh2
    real(dl) rombint
    external rombint

    ombh2 = CP%omegab*(CP%h0/100.0d0)**2
    omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2

    !!From Hu & Sugiyama
    zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
        (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
        (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

    astar = 1/(1+zstar)
    atol = 1e-6
    rs = rombint(dsound_da,1d-8,astar,atol)
    DA = AngularDiameterDistance(zstar)/astar
    CosmomcTheta = rs/DA
    !       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

    end function CosmomcTheta

    end module ModelParams



    !ccccccccccccccccccccccccccccccccccccccccccccccccccc

    module lvalues
    use precision
    use ModelParams
    implicit none
    public

    Type lSamples
        integer l0
        integer l(lmax_arr)
    end Type lSamples

    Type(lSamples) :: lSamp
    !Sources
    logical :: Log_lvalues  = .false.

    contains

    function lvalues_indexOf(lSet,l)
    type(lSamples) :: lSet
    integer, intent(in) :: l
    integer lvalues_indexOf, i

    do i=2,lSet%l0
        if (l < lSet%l(i)) then
            lvalues_indexOf = i-1
            return
        end if
    end do
    lvalues_indexOf = lSet%l0

    end function  lvalues_indexOf

    subroutine initlval(lSet,max_l)

    ! This subroutines initializes lSet%l arrays. Other values will be interpolated.

    implicit none
    type(lSamples) :: lSet

    integer, intent(IN) :: max_l
    integer lind, lvar, step,top,bot,ls(lmax_arr)
    real(dl) AScale

    Ascale=scale/lSampleBoost

    if (lSampleBoost >=50) then
        !just do all of them
        lind=0
        do lvar=lmin, max_l
            lind=lind+1
            ls(lind)=lvar
        end do
        lSet%l0=lind
        lSet%l(1:lind) = ls(1:lind)
        return
    end if

    lind=0
    do lvar=lmin, 10
        lind=lind+1
        ls(lind)=lvar
    end do

    if (CP%AccurateReionization) then
        if (lSampleBoost > 1) then
            do lvar=11, 37,1
                lind=lind+1
                ls(lind)=lvar
            end do
        else
            do lvar=11, 37,2
                lind=lind+1
                ls(lind)=lvar
            end do
        end if

        step = max(nint(5*Ascale),2)
        bot=40
        top=bot + step*10
    else
        if (lSampleBoost >1) then
            do lvar=11, 15
                lind=lind+1
                ls(lind)=lvar
            end do
        else
            lind=lind+1
            ls(lind)=12
            lind=lind+1
            ls(lind)=15
        end if
        step = max(nint(10*Ascale),3)
        bot=15+max(step/2,2)
        top=bot + step*7
    end if

    do lvar=bot, top, step
        lind=lind+1
        ls(lind)=lvar
    end do

    !Sources
    if (Log_lvalues) then
        !Useful for generating smooth things like 21cm to high l
        step=max(nint(20*Ascale),4)
        do
            lvar = lvar + step
            if (lvar > max_l) exit
            lind=lind+1
            ls(lind)=lvar
            step = nint(step*1.2) !log spacing
        end do
    else
        step=max(nint(20*Ascale),4)
        bot=ls(lind)+step
        top=bot+step*2

        do lvar = bot,top,step
            lind=lind+1
            ls(lind)=lvar
        end do

        if (ls(lind)>=max_l) then
            do lvar=lind,1,-1
                if (ls(lvar)<=max_l) exit
            end do
            lind=lvar
            if (ls(lind)<max_l) then
                lind=lind+1
                ls(lind)=max_l
            end if
        else
            step=max(nint(25*Ascale),4)
            !Get EE right around l=200 by putting extra point at 175
            bot=ls(lind)+step
            top=bot+step

            do lvar = bot,top,step
                lind=lind+1
                ls(lind)=lvar
            end do

            if (ls(lind)>=max_l) then
                do lvar=lind,1,-1
                    if (ls(lvar)<=max_l) exit
                end do
                lind=lvar
                if (ls(lind)<max_l) then
                    lind=lind+1
                    ls(lind)=max_l
                end if
            else
                if (HighAccuracyDefault .and. .not. use_spline_template) then
                    step=max(nint(42*Ascale),7)
                else
                    step=max(nint(50*Ascale),7)
                end if
                bot=ls(lind)+step
                top=min(5000,max_l)

                do lvar = bot,top,step
                    lind=lind+1
                    ls(lind)=lvar
                end do

                if (max_l > 5000) then
                    !Should be pretty smooth or tiny out here
                    step=max(nint(400*Ascale),50)
                    lvar = ls(lind)
                    do
                        lvar = lvar + step
                        if (lvar > max_l) exit
                        lind=lind+1
                        ls(lind)=lvar
                        step = nint(step*1.5) !log spacing
                    end do
                end if
                !Sources
            end if !log_lvalues

            if (ls(lind) /=max_l) then
                lind=lind+1
                ls(lind)=max_l
            end if
            if (.not. CP%flat) ls(lind-1)=int(max_l+ls(lind-2))/2
            !Not in CP%flat case so interpolation table is the same when using lower l_max
        end if
    end if
    lSet%l0=lind
    lSet%l(1:lind) = ls(1:lind)

    end subroutine initlval

    subroutine InterpolateClArr(lSet,iCl, all_Cl, max_ind)
    type (lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lmin:*)
    integer, intent(in) :: max_ind
    integer il,llo,lhi, xi
    real(dl) ddCl(lSet%l0)
    real(dl) xl(lSet%l0)

    real(dl) a0,b0,ho
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

    if (max_ind > lSet%l0) call MpiStop('Wrong max_ind in InterpolateClArr')

    xl = real(lSet%l(1:lSet%l0),dl)
    call spline(xl,iCL(1),max_ind,cllo,clhi,ddCl(1))

    llo=1
    do il=lmin,lSet%l(max_ind)
        xi=il
        if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
            llo=llo+1
        end if
        lhi=llo+1
        ho=lSet%l(lhi)-lSet%l(llo)
        a0=(lSet%l(lhi)-xi)/ho
        b0=(xi-lSet%l(llo))/ho

        all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
            +(b0**3-b0)*ddCl(lhi))*ho**2/6
    end do

    end subroutine InterpolateClArr

    subroutine InterpolateClArrTemplated(lSet,iCl, all_Cl, max_ind, template_index)
    type (lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lmin:*)
    integer, intent(in) :: max_ind
    integer, intent(in), optional :: template_index
    integer maxdelta, il
    real(dl) DeltaCL(lSet%l0)
    real(dl), allocatable :: tmpall(:)

    if (max_ind > lSet%l0) call MpiStop('Wrong max_ind in InterpolateClArrTemplated')

    if (use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
            !interpolate only the difference between the C_l and an accurately interpolated template. Temp only for the mo.
            !Using unlensed for template, seems to be good enough
            maxdelta=max_ind
            do while (lSet%l(maxdelta) > lmax_extrap_highl)
                maxdelta=maxdelta-1
            end do
            DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

            call InterpolateClArr(lSet,DeltaCl, all_Cl, maxdelta)

            do il=lmin,lSet%l(maxdelta)
                all_Cl(il) = all_Cl(il) +  highL_CL_template(il,template_index)
            end do

            if (maxdelta < max_ind) then
                !directly interpolate high L where no t  emplate (doesn't effect lensing spectrum much anyway)
                allocate(tmpall(lmin:lSet%l(max_ind)))
                call InterpolateClArr(lSet,iCl, tmpall, max_ind)
                !overlap to reduce interpolation artefacts
                all_cl(lSet%l(maxdelta-2):lSet%l(max_ind) ) = tmpall(lSet%l(maxdelta-2):lSet%l(max_ind))
                deallocate(tmpall)
            end if
            return
        end if
    end if

    call InterpolateClArr(lSet,iCl, all_Cl, max_ind)


    end subroutine InterpolateClArrTemplated

    end module lvalues


    !ccccccccccccccccccccccccccccccccccccccccccccccccccc

    module ModelData
    use precision
    use ModelParams
    use InitialPower
    use lValues
    use Ranges
    use AMlUtils
    implicit none
    public

    Type LimberRec
        integer n1,n2 !corresponding time step array indices
        real(dl), dimension(:), pointer :: k  => NULL()
        real(dl), dimension(:), pointer :: Source  => NULL()
    end Type LimberRec

    Type ClTransferData
        !Cl transfer function variables
        !values of q for integration over q to get C_ls
        Type (lSamples) :: ls ! scalar and tensor l that are computed
        integer :: NumSources
        !Changes -scalars:  2 for just CMB, 3 for lensing
        !- tensors: T and E and phi (for lensing), and T, E, B respectively

        Type (Regions) :: q
        real(dl), dimension(:,:,:), pointer :: Delta_p_l_k => NULL()

        !The L index of the lowest L to use for Limber
        integer, dimension(:), pointer :: Limber_l_min => NULL()
        !For each l, the set of k in each limber window
        !indices LimberWindow(SourceNum,l)
        Type(LimberRec), dimension(:,:), pointer :: Limber_windows => NULL()

        !The maximum L needed for non-Limber
        integer max_index_nonlimber

    end Type ClTransferData

    Type(ClTransferData), save, target :: CTransScal, CTransTens, CTransVec

    !Computed output power spectra data

    integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4, C_PhiTemp = 5, C_PhiE=6
    integer :: C_last = C_PhiE
    integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4
    integer, parameter :: name_tag_len = 12
    character(LEN=name_tag_len), dimension(C_PhiE), parameter :: C_name_tags = ['TT','EE','TE','PP','TP','EP']
    character(LEN=name_tag_len), dimension(CT_Cross), parameter :: CT_name_tags = ['TT','EE','BB','TE']
    character(LEN=name_tag_len), dimension(7), parameter :: lens_pot_name_tags = ['TT','EE','BB','TE','PP','TP','EP']


    logical :: has_cl_2D_array = .false.

    real(dl), dimension (:,:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
    !Indices are Cl_xxx( l , intial_power_index, Cl_type)
    !where Cl_type is one of the above constants

    real(dl), dimension (:,:,:,:), allocatable :: Cl_Scalar_Array
    !Indices are Cl_xxx( l , intial_power_index, field1,field2)
    !where ordering of fields is T, E, \psi (CMB lensing potential), window_1, window_2...

    !The following are set only if doing lensing
    integer lmax_lensed !Only accurate to rather less than this
    real(dl) , dimension (:,:,:), allocatable :: Cl_lensed
    !Cl_lensed(l, power_index, Cl_type) are the interpolated Cls

    contains

    subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call Ranges_getArray(CTrans%q, .true.)

    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,&
        min(CTrans%max_index_nonlimber,CTrans%ls%l0), CTrans%q%npoints),  STAT = st)
    if (st /= 0) call MpiStop('Init_ClTransfer: Error allocating memory for transfer functions')
    CTrans%Delta_p_l_k = 0

    end subroutine Init_ClTransfer

    subroutine Init_Limber(CTrans)
    Type(ClTransferData) :: CTrans

    allocate(CTrans%Limber_l_min(CTrans%NumSources))
    CTrans%Limber_l_min = 0
    if (num_redshiftwindows>0 .or. limber_phiphi>0) then
        allocate(CTrans%Limber_windows(CTrans%NumSources,CTrans%ls%l0))
    end if

    end subroutine Init_Limber

    subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    nullify(CTrans%Delta_p_l_k)
    call Ranges_Free(CTrans%q)
    call Free_Limber(CTrans)

    end subroutine Free_ClTransfer

    subroutine Free_Limber(CTrans)
    Type(ClTransferData) :: CTrans
    integer st,i,j

    if (associated(CTrans%Limber_l_min)) then
        do i=1, CTrans%NumSources
            if (CTrans%Limber_l_min(i)/=0) then
                do j=CTrans%Limber_l_min(i), CTrans%ls%l0
                    deallocate(CTrans%Limber_windows(i, j)%k, STAT = st)
                    deallocate(CTrans%Limber_windows(i, j)%Source, STAT = st)
                end do
            end if
        end do
        deallocate(CTrans%Limber_l_min, STAT = st)
    end if
    deallocate(CTrans%Limber_windows, STAT = st)
    nullify(CTrans%Limber_l_min)
    nullify(CTrans%Limber_windows)

    end subroutine Free_Limber

    subroutine CheckLoadedHighLTemplate
    integer L
    real(dl) array(7)

    if (.not. allocated(highL_CL_template)) then
        allocate(highL_CL_template(lmin:lmax_extrap_highl, C_Temp:C_Phi))

        call OpenTxtFile(highL_unlensed_cl_template,fileio_unit)

        if (lmin==1) highL_CL_template(lmin,:)=0
        do
            read(fileio_unit,*, end=500) L , array
            if (L>lmax_extrap_highl) exit
            !  array = array * (2*l+1)/(4*pi) * 2*pi/(l*(l+1))
            highL_CL_template(L, C_Temp:C_E) =array(1:2)
            highL_CL_template(L, C_Cross) =array(4)
            highL_CL_template(L, C_Phi) =array(5)
        end do
500     if (L< lmax_extrap_highl) &
            call MpiStop('CheckLoadedHighLTemplate: template file does not go up to lmax_extrap_highl')
        close(fileio_unit)
    end if

    end subroutine CheckLoadedHighLTemplate


    subroutine Init_Cls

    call CheckLoadedHighLTemplate
    if (CP%WantScalars) then
        if (allocated(Cl_scalar)) deallocate(Cl_scalar)
        allocate(Cl_scalar(lmin:CP%Max_l, CP%InitPower%nn, C_Temp:C_last))
        Cl_scalar = 0
        if (has_cl_2D_array) then
            if (allocated(Cl_scalar_array)) deallocate(Cl_scalar_array)
            allocate(Cl_scalar_Array(lmin:CP%Max_l, CP%InitPower%nn, &
                3+num_redshiftwindows+num_custom_sources,3+num_redshiftwindows+num_custom_sources))
            Cl_scalar_array = 0
        end if
    end if

    if (CP%WantVectors) then
        if (allocated(Cl_vector)) deallocate(Cl_vector)
        allocate(Cl_vector(lmin:CP%Max_l, CP%InitPower%nn, CT_Temp:CT_Cross))
        Cl_vector = 0
    end if


    if (CP%WantTensors) then
        if (allocated(Cl_tensor)) deallocate(Cl_tensor)
        allocate(Cl_tensor(lmin:CP%Max_l_tensor, CP%InitPower%nn, CT_Temp:CT_Cross))
        Cl_tensor = 0
    end if

    end subroutine Init_Cls

    function open_file_header(filename, Col1, Columns, n) result(unit)
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: col1
    character(LEN=name_tag_len), intent(in) :: Columns(:)
    integer, intent(in), optional :: n
    integer :: unit, nn

    if (present(n)) then
        nn = n
    else
        nn = 6
    end if
    open(newunit=unit,file=filename,form='formatted',status='replace')
    if (output_file_headers) then
        write(unit,'("#",1A'//Trim(IntToStr(nn-1))//'," ",*(A15))') Col1,Columns
    end if

    end function open_file_header

    function scalar_fieldname(i)
    integer, intent(in) :: i
    character(LEN=5) :: scalar_fieldname
    character(LEN=3), parameter :: scalar_fieldnames = 'TEP'

    if (i<=3) then
        scalar_fieldname = scalar_fieldnames(i:i)
    else
        scalar_fieldname = 'W'//trim(IntToStr(i-3))
    end if

    end function scalar_fieldname

    subroutine output_cl_files(ScalFile,ScalCovFile,TensFile, TotFile, LensFile, LensTotFile, factor)
    implicit none
    integer in,il, i, j
    character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile,ScalCovfile
    real(dl), intent(in), optional :: factor
    real(dl) fact
    integer last_C
    real(dl), allocatable :: outarr(:,:)
    integer unit
    character(LEN=name_tag_len) :: cov_names((3+num_redshiftwindows)**2)


    if (present(factor)) then
        fact = factor
    else
        fact =1
    end if

    if (CP%WantScalars .and. ScalFile /= '') then
        last_C=min(C_PhiTemp,C_last)
        unit = open_file_header(ScalFile, 'L', C_name_tags(:last_C))
        do in=1,CP%InitPower%nn
            do il=lmin,min(10000,CP%Max_l)
                write(unit,trim(numcat('(1I6,',last_C))//'E15.6)')il ,fact*Cl_scalar(il,in,C_Temp:last_C)
            end do
            do il=10100,CP%Max_l, 100
                write(unit,trim(numcat('(1E15.6,',last_C))//'E15.6)') real(il),&
                    fact*Cl_scalar(il,in,C_Temp:last_C)
            end do
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. has_cl_2D_array .and. ScalCovFile /= '' .and. CTransScal%NumSources>2) then
        allocate(outarr(1:3+num_redshiftwindows,1:3+num_redshiftwindows))
        do i=1, 3+num_redshiftwindows
            do j=1, 3+num_redshiftwindows
                cov_names(j + (i-1)*(3+num_redshiftwindows)) = trim(scalar_fieldname(i))//'x'//trim(scalar_fieldname(j))
            end do
        end do
        unit = open_file_header(ScalCovFile, 'L', cov_names)

        do in=1,CP%InitPower%nn
            do il=lmin,min(10000,CP%Max_l)
                outarr=Cl_scalar_array(il,in,1:3+num_redshiftwindows,1:3+num_redshiftwindows)
                outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
                outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
                write(unit,trim(numcat('(1I6,',(3+num_redshiftwindows)**2))//'E15.6)') il, outarr
            end do
            do il=10100,CP%Max_l, 100
                outarr=Cl_scalar_array(il,in,1:3+num_redshiftwindows,1:3+num_redshiftwindows)
                outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
                outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
                write(unit,trim(numcat('(1E15.6,',(3+num_redshiftwindows)**2))//'E15.6)') real(il), outarr
            end do
        end do
        close(unit)
        deallocate(outarr)
    end if

    if (CP%WantTensors .and. TensFile /= '') then
        unit = open_file_header(TensFile, 'L', CT_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin,CP%Max_l_tensor
                write(unit,'(1I6,4E15.6)')il, fact*Cl_tensor(il, in, CT_Temp:CT_Cross)
            end do
        end do
        close(unit)
    end if

    if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
        unit = open_file_header(TotFile, 'L', CT_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin,CP%Max_l_tensor
                write(unit,'(1I6,4E15.6)')il, fact*(Cl_scalar(il, in, C_Temp:C_E)+ Cl_tensor(il,in, C_Temp:C_E)), &
                    fact*Cl_tensor(il,in, CT_B), fact*(Cl_scalar(il, in, C_Cross) + Cl_tensor(il, in, CT_Cross))
            end do
            do il=CP%Max_l_tensor+1,CP%Max_l
                write(unit,'(1I6,4E15.6)')il ,fact*Cl_scalar(il,in,C_Temp:C_E), 0._dl, fact*Cl_scalar(il,in,C_Cross)
            end do
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
        unit = open_file_header(LensFile, 'L', CT_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin, lmax_lensed
                write(unit,'(1I6,4E15.6)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
            end do
        end do
        close(unit)
    end if


    if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
        unit = open_file_header(LensTotFile, 'L', CT_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin,min(CP%Max_l_tensor,lmax_lensed)
                write(unit,'(1I6,4E15.6)')il, fact*(Cl_lensed(il, in, CT_Temp:CT_Cross)+ Cl_tensor(il,in, CT_Temp:CT_Cross))
            end do
            do il=min(CP%Max_l_tensor,lmax_lensed)+1,lmax_lensed
                write(unit,'(1I6,4E15.6)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
            end do
        end do
        close(unit)
    end if
    end subroutine output_cl_files

    subroutine output_lens_pot_files(LensPotFile, factor)
    !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed
    !This input supported by LensPix from 2010
    implicit none
    integer in,il
    real(dl), intent(in), optional :: factor
    real(dl) fact, scale, BB, TT, TE, EE
    character(LEN=*) LensPotFile
    integer unit
    !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi and [l(l+1)]^(3/2) C_phi_T/2pi
    !This is the format used by Planck_like but original LensPix uses scalar_output_file.

    !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi
    ! - for historical reasons)

    if (present(factor)) then
        fact = factor
    else
        fact =1
    end if

    if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then
        unit =  open_file_header(LensPotFile, 'L', lens_pot_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin,min(10000,CP%Max_l)
                TT = Cl_scalar(il, in, C_Temp)
                EE = Cl_scalar(il, in, C_E)
                TE = Cl_scalar(il, in, C_Cross)
                if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                    TT= TT+Cl_tensor(il,in, CT_Temp)
                    EE= EE+Cl_tensor(il,in, CT_E)
                    TE= TE+Cl_tensor(il,in, CT_Cross)
                    BB= Cl_tensor(il,in, CT_B)
                else
                    BB=0
                end if
                scale = (real(il+1)/il)**2/OutputDenominator !Factor to go from old l^4 factor to new

                write(unit,'(1I6,7E15.6)') il , fact*TT, fact*EE, fact*BB, fact*TE, scale*Cl_scalar(il,in,C_Phi),&
                    (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*Cl_scalar(il,in,C_PhiTemp:C_PhiE)
            end do
            do il=10100,CP%Max_l, 100
                scale = (real(il+1)/il)**2/OutputDenominator
                write(unit,'(1E15.6,7E15.6)') real(il), fact*Cl_scalar(il,in,C_Temp:C_E),0.,fact*Cl_scalar(il,in,C_Cross), &
                    scale*Cl_scalar(il,in,C_Phi),&
                    (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*Cl_scalar(il,in,C_PhiTemp:C_PhiE)
            end do
        end do
        close(unit)
    end if
    end subroutine output_lens_pot_files


    subroutine output_veccl_files(VecFile, factor)
    implicit none
    integer in,il
    character(LEN=*) VecFile
    real(dl), intent(in), optional :: factor
    real(dl) fact
    integer unit

    if (present(factor)) then
        fact = factor
    else
        fact =1
    end if


    if (CP%WantVectors .and. VecFile /= '') then
        unit =  open_file_header(VecFile, 'L', CT_name_tags)
        do in=1,CP%InitPower%nn
            do il=lmin,CP%Max_l
                write(unit,'(1I6,4E15.6)')il, fact*Cl_vector(il, in, CT_Temp:CT_Cross)
            end do
        end do
        close(unit)
    end if

    end subroutine output_veccl_files

    subroutine NormalizeClsAtL(lnorm)
    implicit none
    integer, intent(IN) :: lnorm
    integer in
    real(dl) Norm

    do in=1,CP%InitPower%nn
        if (CP%WantScalars) then
            Norm=1/Cl_scalar(lnorm,in, C_Temp)
            Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) = Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) * Norm
        end if

        if (CP%WantTensors) then
            if (.not.CP%WantScalars) Norm = 1/Cl_tensor(lnorm,in, C_Temp)
            !Otherwise Norm already set correctly
            Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) =  &
                Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) * Norm
        end if
    end do

    end  subroutine NormalizeClsAtL

    subroutine ModelData_Free

    call Free_ClTransfer(CTransScal)
    call Free_ClTransfer(CTransVec)
    call Free_ClTransfer(CTransTens)
    if (allocated(Cl_vector)) deallocate(Cl_vector)
    if (allocated(Cl_tensor)) deallocate(Cl_tensor)
    if (allocated(Cl_scalar)) deallocate(Cl_scalar)
    if (allocated(Cl_lensed)) deallocate(Cl_lensed)
    if (allocated(Cl_scalar_array)) deallocate(Cl_scalar_array)

    end subroutine ModelData_Free

    end module ModelData


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    module MassiveNu
    use precision
    use ModelParams
    implicit none
    private

    real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
    !const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
    real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
    real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
    real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

    ! zeta3*3/2/pi^2*4/11*((k_B*COBE_CMBTemp/hbar/c)^3* 8*pi*G/3/(100*km/s/megaparsec)^2/(c^2/eV)
    real(dl), parameter :: neutrino_mass_fac= 94.07_dl !converts omnuh2 into sum m_nu in eV

    integer, parameter  :: nrhopn=2000
    real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 600._dl
    !max a*m_nu to integrate

    real(dl),parameter  :: am_minp=am_min*1.1
    real(dl), parameter :: am_maxp=am_max*0.9

    real(dl) dlnam

    real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1

    !Sample for massive neutrino momentum
    !These settings appear to be OK for P_k accuate at 1e-3 level
    integer, parameter :: nqmax0=80 !maximum array size of q momentum samples
    real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0)

    integer nqmax !actual number of q modes evolves

    public const,Nu_Init,Nu_background, Nu_rho, Nu_drho,  nqmax0, nqmax, &
        nu_int_kernel, nu_q, sum_mnu_for_m1, neutrino_mass_fac
    contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sum_mnu_for_m1(summnu,dsummnu, m1, targ, sgn)
    use constants
    real(dl), intent(in) :: m1, targ, sgn
    real(dl), intent(out) :: summnu, dsummnu
    real(dl) :: m2,m3

    m2 = sqrt(m1**2 + delta_mnu21)
    m3 = sqrt(m1**2 + sgn*delta_mnu31)
    summnu = m1 + m2 + m3 - targ
    dsummnu = m1/m2+m1/m3 + 1

    end subroutine sum_mnu_for_m1

    subroutine Nu_init
    !  Initialize interpolation tables for massive neutrinos.
    !  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
    integer i
    real(dl) dq,dlfdlq, q, am, rhonu,pnu
    real(dl) spline_data(nrhopn)

    !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
    !  Get number density n of neutrinos from
    !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
    !  then m = Omega_nu/N_nu rho_crit /n
    !  Error due to velocity < 1e-5

    do i=1, CP%Nu_mass_eigenstates
        nu_masses(i)=const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
            /CP%Nu_mass_degeneracies(i)
    end do

    if (allocated(r1)) return
    allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))


    nqmax=3
    if (AccuracyBoost >1) nqmax=4
    if (AccuracyBoost >2) nqmax=5
    if (AccuracyBoost >3) nqmax=nint(AccuracyBoost*10)
    !note this may well be worse than the 5 optimized points

    if (nqmax > nqmax0) call MpiStop('Nu_Init: qmax > nqmax0')

    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes
    if (nqmax==3) then
        !Accurate at 2e-4 level
        nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
        nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
    else if (nqmax==4) then
        !This seems to be very accurate (limited by other numerics)
        nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
        nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
    else if (nqmax==5) then
        !exact for n=-4,-2..3
        !This seems to be very accurate (limited by other numerics)
        nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
        nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)
    else
        dq = (12 + nqmax/5)/real(nqmax)
        do i=1,nqmax
            q=(i-0.5d0)*dq
            nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
    end if
    nu_int_kernel=nu_int_kernel/const

    dlnam=-(log(am_min/am_max))/(nrhopn-1)


    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
    !$OMP & PRIVATE(am, rhonu,pnu)
    do i=1,nrhopn
        am=am_min*exp((i-1)*dlnam)
        call nuRhoPres(am,rhonu,pnu)
        r1(i)=log(rhonu)
        p1(i)=log(pnu)
    end do
    !$OMP END PARALLEL DO


    call splini(spline_data,nrhopn)
    call splder(r1,dr1,nrhopn,spline_data)
    call splder(p1,dp1,nrhopn,spline_data)
    call splder(dr1,ddr1,nrhopn,spline_data)


    end subroutine Nu_init

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine nuRhoPres(am,rhonu,pnu)
    !  Compute the density and pressure of one eigenstate of massive neutrinos,
    !  in units of the mean density of one flavor of massless neutrinos.

    real(dl),  parameter :: qmax=30._dl
    integer, parameter :: nq=100
    real(dl) dum1(nq+1),dum2(nq+1)
    real(dl), intent(in) :: am
    real(dl), intent(out) ::  rhonu,pnu
    integer i
    real(dl) q,aq,v,aqdn,adq


    !  q is the comoving momentum in units of k_B*T_nu0/c.
    !  Integrate up to qmax and then use asymptotic expansion for remainder.
    adq=qmax/nq
    dum1(1)=0._dl
    dum2(1)=0._dl
    do  i=1,nq
        q=i*adq
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        aqdn=adq*q*q*q/(exp(q)+1._dl)
        dum1(i+1)=aqdn/v
        dum2(i+1)=aqdn*v
    end do
    call splint(dum1,rhonu,nq+1)
    call splint(dum2,pnu,nq+1)
    !  Apply asymptotic corrrection for q>qmax and normalize by relativistic
    !  energy density.
    rhonu=(rhonu+dum1(nq+1)/adq)/const
    pnu=(pnu+dum2(nq+1)/adq)/const/3._dl

    end subroutine nuRhoPres

    !cccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_background(am,rhonu,pnu)
    use precision
    use ModelParams
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu

    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.

    real(dl) d
    integer i

    if (am <= am_minp) then
        rhonu=1._dl + const2*am**2
        pnu=(2-rhonu)/3._dl
        return
    else if (am >= am_maxp) then
        rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
        pnu = 900._dl/120._dl/const*(zeta5-63._dl/4*Zeta7/am**2)/am
        return
    end if


    d=log(am/am_min)/dlnam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
        -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
    pnu=p1(i)+d*(dp1(i)+d*(3._dl*(p1(i+1)-p1(i))-2._dl*dp1(i) &
        -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2._dl*(p1(i)-p1(i+1)))))
    rhonu=exp(rhonu)
    pnu=exp(pnu)

    end subroutine Nu_background

    !cccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_rho(am,rhonu)
    use precision
    use ModelParams
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu

    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.

    real(dl) d
    integer i

    if (am <= am_minp) then
        rhonu=1._dl + const2*am**2
        return
    else if (am >= am_maxp) then
        rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
        return
    end if

    d=log(am/am_min)/dlnam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
        -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
    rhonu=exp(rhonu)
    end subroutine Nu_rho

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    function Nu_drho(am,adotoa,rhonu) result (rhonudot)
    use precision
    use ModelParams

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) adotoa,rhonu,rhonudot
    real(dl) d
    real(dl), intent(IN) :: am
    integer i

    if (am< am_minp) then
        rhonudot = 2*const2*am**2*adotoa
    else if (am>am_maxp) then
        rhonudot = 3/(2*const)*(zeta3*am - (15*zeta5)/2/am)*adotoa
    else
        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i
        !  Cubic spline interpolation for rhonudot.
        rhonudot=dr1(i)+d*(ddr1(i)+d*(3._dl*(dr1(i+1)-dr1(i)) &
            -2._dl*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1) &
            +2._dl*(dr1(i)-dr1(i+1)))))

        rhonudot=rhonu*adotoa*rhonudot/dlnam
    end if

    end function Nu_drho

    end module MassiveNu

    ! wrapper function to avoid cirular module references
    subroutine init_massive_nu(has_massive_nu)
    use MassiveNu
    use ModelParams
    implicit none
    logical, intent(IN) :: has_massive_nu

    if (has_massive_nu) then
        call Nu_Init
    else
        nu_masses = 0
    end if
    end subroutine init_massive_nu


    !ccccccccccccccccccccccccccccccccccccccccccccccccccc

    module Transfer
    use ModelData
    use Errors
    implicit none
    public
    integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
        Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
        Transfer_tot=7, Transfer_nonu=8, Transfer_tot_de=9,  &
        ! total perturbations with and without neutrinos, with neutrinos+dark energy in the numerator
        Transfer_Weyl = 10, & ! the Weyl potential, for lensing and ISW
        Transfer_Newt_vel_cdm=11, Transfer_Newt_vel_baryon=12,   & ! -k v_Newtonian/H
        Transfer_vel_baryon_cdm = 13 !relative velocity of baryons and CDM

    integer, parameter :: Transfer_max = Transfer_vel_baryon_cdm
    character(LEN=name_tag_len) :: Transfer_name_tags(Transfer_max-1) = &
        ['CDM     ', 'baryon  ', 'photon  ', 'nu      ', 'mass_nu ', 'total   ', &
        'no_nu   ', 'total_de', 'Weyl    ', 'v_CDM   ', 'v_b     ', 'v_b-v_c ']

    logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
    !set to false to output calculated values for later interpolation

    integer :: transfer_power_var = Transfer_tot
    !What to use to calulcate the output matter power spectrum and sigma_8
    !Transfer_tot uses total matter perturbation

    logical :: get_growth_sigma8 = .true.
    !gets sigma_vdelta, like sigma8 but using velocity-density cross power,
    !in late LCDM f*sigma8 = sigma_vdelta^2/sigma8

    Type MatterTransferData
        !Computed data
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        real(dl), dimension (:), pointer :: q_trans => NULL()
        real(dl), dimension (:,:), pointer ::  sigma_8 => NULL()
        real(dl), dimension (:,:), pointer ::  sigma2_vdelta_8 => NULL() !growth from sigma_{v delta}
        real, dimension(:,:,:), pointer :: TransferData => NULL()
        !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
    end Type MatterTransferData

    Type MatterPowerData
        !everything is a function of k/h
        integer   ::  num_k, num_z
        real(dl), dimension(:), pointer :: log_kh => NULL(), redshifts => NULL()
        !matpower is log(P_k)
        real(dl), dimension(:,:), allocatable :: matpower, ddmat
        !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
        !function of k and redshift NonLinearScaling(k_index,z_index)
        real(dl), dimension(:,:), pointer :: nonlin_ratio => NULL()
    end Type MatterPowerData

    Type (MatterTransferData), save :: MT

    interface Transfer_GetMatterPower
    module procedure Transfer_GetMatterPowerD,Transfer_GetMatterPowerS
    end interface

    contains

    subroutine Transfer_GetUnsplinedPower(M,PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k)
    Type(MatterTransferData) :: M
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    real(dl) h, k
    integer nz, nk, zix, ik
    integer s1, s2
    logical hnorm

    s1 = transfer_power_var
    if (present(var1))  s1 = var1
    s2 = transfer_power_var
    if (present(var2))  s2 = var2
    hnorm = .true.
    if (present(hubble_units)) hnorm = hubble_units

    nk=M%num_q_trans
    nz=CP%Transfer%PK_num_redshifts
    if (nk/= size(PK,1) .or. nz/=size(PK,2)) call MpiStop('Trasfer_GetUnsplinedPower wrong size')

    h = CP%H0/100

    do ik=1,nk
        k = M%TransferData(Transfer_kh,ik,1)*h
        do zix=1,nz
            PK(ik,zix) = M%TransferData(s1,ik,CP%Transfer%PK_redshifts_index(nz-zix+1))*&
                M%TransferData(s2,ik,CP%Transfer%PK_redshifts_index(nz-zix+1))*k*pi*twopi*scalarPower(k,1)
        end do
    end do
    if (hnorm) PK=  PK * h**3

    end subroutine Transfer_GetUnsplinedPower

    subroutine Transfer_GetUnsplinedNonlinearPower(M,PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k) after re-scaling for non-linear evolution (if turned on)
    Type(MatterTransferData), intent(in) :: M
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    Type(MatterPowerData) :: PKdata
    integer zix

    call Transfer_GetUnsplinedPower(M,PK,var1,var2, hubble_units)
    do zix=1, CP%Transfer%PK_num_redshifts
        call Transfer_GetMatterPowerData(M, PKdata, 1, &
            CP%Transfer%PK_redshifts_index(CP%Transfer%PK_num_redshifts-zix+1))
        call NonLinear_GetRatios(PKdata)
        PK(:,zix) =  PK(:,zix) *PKdata%nonlin_ratio(:,1)**2
        call MatterPowerdata_Free(PKdata)
    end do

    end subroutine Transfer_GetUnsplinedNonlinearPower

    subroutine Transfer_GetMatterPowerData(MTrans, PK_data, power_ix, itf_only, var1, var2)
    !Does *NOT* include non-linear corrections
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !spectrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK_data
    integer, intent(in), optional :: power_ix
    integer, intent(in), optional :: itf_only
    integer, intent(in), optional :: var1, var2
    real(dl) h, kh, k, power
    integer ik
    integer nz,itf, itf_start, itf_end
    integer :: s1,s2, p_ix


    s1 = transfer_power_var
    if (present(var1))  s1 = var1
    s2 = transfer_power_var
    if (present(var2))  s2 = var2
    p_ix = 1
    if (present(power_ix)) p_ix = power_ix

    if (present(itf_only)) then
        itf_start=itf_only
        itf_end = itf_only
        nz = 1
    else
        itf_start=1
        nz= size(MTrans%TransferData,3)
        itf_end = nz
    end if
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))
    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = CP%Transfer%Redshifts(itf_start:itf_end)

    h = CP%H0/100

    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        k = kh*h
        PK_data%log_kh(ik) = log(kh)
        power = ScalarPower(k,p_ix)
        if (global_error_flag/=0) then
            call MatterPowerdata_Free(PK_data)
            return
        end if
        do itf = 1, nz
            PK_data%matpower(ik,itf) = &
                log(MTrans%TransferData(s1,ik,itf_start+itf-1)*&
                MTrans%TransferData(s2,ik,itf_start+itf-1)*k &
                *pi*twopi*h**3*power)
        end do
    end do

    call MatterPowerdata_getsplines(PK_data)

    end subroutine Transfer_GetMatterPowerData

    subroutine MatterPowerData_Load(PK_data,fname)
    !Loads in kh, P_k from file for one redshiftr and one initial power spectrum
    !Not redshift is not stored in file, so not set correctly
    !Also note that output _matterpower file is already interpolated, so re-interpolating is probs not a good idea

    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    use AmlUtils
    character(LEN=*) :: fname
    Type(MatterPowerData) :: PK_data
    real(dl)kh, Pk
    integer ik
    integer nz


    nz = 1
    call openTxtFile(fname, fileio_unit)

    PK_data%num_k = FileLines(fileio_unit)
    PK_Data%num_z = 1

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))

    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = 0

    do ik=1,PK_data%num_k
        read (fileio_unit,*) kh, Pk
        PK_data%matpower(ik,1) = log(Pk)
        PK_data%log_kh(ik) = log(kh)
    end do

    call MatterPowerdata_getsplines(PK_data)

    end subroutine MatterPowerData_Load


    subroutine MatterPowerdata_getsplines(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

    do i = 1,PK_Data%num_z
        call spline(PK_data%log_kh,PK_data%matpower(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddmat(1,i))
    end do

    end subroutine MatterPowerdata_getsplines

    subroutine MatterPowerdata_MakeNonlinear(PK_data)
    Type(MatterPowerData) :: PK_data

    call NonLinear_GetRatios(PK_data)
    PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
    call MatterPowerdata_getsplines(PK_data)

    end subroutine MatterPowerdata_MakeNonlinear

    subroutine MatterPowerdata_Free(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    deallocate(PK_data%log_kh,stat=i)
    deallocate(PK_data%matpower,stat=i)
    deallocate(PK_data%ddmat,stat=i)
    deallocate(PK_data%nonlin_ratio,stat=i)
    deallocate(PK_data%redshifts,stat=i)
    call MatterPowerdata_Nullify(PK_data)

    end subroutine MatterPowerdata_Free

    subroutine MatterPowerdata_Nullify(PK_data)
    Type(MatterPowerData) :: PK_data

    nullify(PK_data%log_kh)
    nullify(PK_data%nonlin_ratio)
    nullify(PK_data%redshifts)

    end subroutine MatterPowerdata_Nullify

    function MatterPowerData_k(PK,  kh, itf) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: kh
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho,a0,b0
    integer, save :: i_last = 1

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        dp = (PK%matpower(2,itf) -  PK%matpower(1,itf)) / &
            ( PK%log_kh(2)-PK%log_kh(1) )
        outpower = PK%matpower(1,itf) + dp*(logk - PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter

        dp = (PK%matpower(PK%num_k,itf) -  PK%matpower(PK%num_k-1,itf)) / &
            ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
        outpower = PK%matpower(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
    else
        llo=min(i_last,PK%num_k)
        do while (PK%log_kh(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_kh(llo+1)< logk)
            llo=llo+1
        end do
        i_last =llo
        lhi=llo+1
        ho=PK%log_kh(lhi)-PK%log_kh(llo)
        a0=(PK%log_kh(lhi)-logk)/ho
        b0=1-a0

        outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            ((a0**3-a0)* PK%ddmat(llo,itf) &
            +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6
    end if

    outpower = exp(outpower)

    end function MatterPowerData_k

    subroutine Transfer_GetMatterPowerS(MTrans,outpower, itf, in, minkh, dlnkh, npoints, var1, var2)
    Type(MatterTransferData), intent(in) :: MTrans
    integer, intent(in) :: itf, in, npoints
    integer, intent(in), optional :: var1, var2
    real, intent(out) :: outpower(*)
    real, intent(in) :: minkh, dlnkh
    real(dl) :: outpowerd(npoints)
    real(dl):: minkhd, dlnkhd

    minkhd = minkh; dlnkhd = dlnkh
    call Transfer_GetMatterPowerD(MTrans,outpowerd, itf, in, minkhd, dlnkhd, npoints,var1, var2)
    outpower(1:npoints) = outpowerd(1:npoints)

    end subroutine Transfer_GetMatterPowerS

    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    !Changed input variable from itf to itf_PK because we are looking for the itf_PK'th
    !redshift in the PK_redshifts array.  The position of this redshift in the master redshift
    !array, itf, is given by itf = CP%Transfer%Pk_redshifts_index(itf_PK)
    !Also changed (CP%NonLinear/=NonLinear_None) to
    !CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
    subroutine Transfer_GetMatterPowerD(MTrans,outpower, itf_PK, in, minkh, dlnkh, npoints, var1, var2)
    !Allows for non-smooth priordial spectra
    !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
    !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
    !in units of (h Mpc^{-1})^3.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK

    integer, intent(in) :: itf_PK, in, npoints
    real(dl), intent(out) :: outpower(npoints)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in), optional :: var1, var2

    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    integer ik, llo,il,lhi,lastix
    real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
    real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
    integer itf
    integer :: s1,s2

    s1 = transfer_power_var
    if (present(var1))  s1 = var1
    s2 = transfer_power_var
    if (present(var2))  s2 = var2

    itf = CP%Transfer%PK_redshifts_index(itf_PK)

    if (npoints < 2) call MpiStop('Need at least 2 points in Transfer_GetMatterPower')

    !         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
    !            stop 'Transfer_GetMatterPower: kh out of computed region'
    !          end if
    if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
        .and. FeedbackLevel > 0 ) &
        write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'


    if (CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens) then
        call Transfer_GetMatterPowerData(MTrans, PK, in, itf) ! Mar 16, changed to use default variable
        call NonLinear_GetRatios(PK)
    end if

    h = CP%H0/100
    logmink = log(minkh)
    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,itf)
        k = kh*h
        kvals(ik) = log(kh)
        atransfer=MTrans%TransferData(s1,ik,itf)*MTrans%TransferData(s2,ik,itf)
        if (CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens) &
            atransfer = atransfer* PK%nonlin_ratio(ik,1)**2 !only one element, this itf
        matpower(ik) = log(atransfer*k*pi*twopi*h**3)
        !Put in power spectrum later: transfer functions should be smooth, initial power may not be
    end do

    call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)

    llo=1
    lastix = npoints + 1
    do il=1, npoints
        xi=logmink + dlnkh*(il-1)
        if (xi < kvals(1)) then
            outpower(il)=-30.
            cycle
        end if
        do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
            llo=llo+1
            if (llo >= MTrans%num_q_trans) exit
        end do
        if (llo == MTrans%num_q_trans) then
            lastix = il
            exit
        end if
        lhi=llo+1
        ho=kvals(lhi)-kvals(llo)
        a0=(kvals(lhi)-xi)/ho
        b0=(xi-kvals(llo))/ho

        outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
            +(b0**3-b0)*ddmat(lhi))*ho**2/6
    end do

    do while (lastix <= npoints)
        !Do linear extrapolation in the log
        !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
        outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
        lastix = lastix+1
    end do

    outpower = exp(max(-30.d0,outpower))

    do il = 1, npoints
        k = exp(logmink + dlnkh*(il-1))*h
        outpower(il) = outpower(il) * ScalarPower(k,in)
        if (global_error_flag /= 0) exit
    end do

    if (CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens) call MatterPowerdata_Free(PK)

    end subroutine Transfer_GetMatterPowerD

    subroutine Transfer_Get_SigmaR(MTrans, R, outvals, var1, var2, power_ix, root)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc, for all requested redshifts
    !set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R
    integer, intent(in), optional :: var1, var2
    integer, intent(in), optional :: power_ix
    logical, intent(in), optional :: root !if true, give sigma8, otherwise sigma8^2
    real(dl), intent(out) :: outvals(:)
    integer ik
    real(dl) kh, k, h, x, win
    real(dl) lnk, dlnk, lnko
    real(dl), dimension(CP%Transfer%PK_num_redshifts) ::  dsig8, dsig8o, sig8, sig8o
    real(dl) powers
    integer s1,s2
    integer :: ix = 1

    s1 = transfer_power_var
    if (present(var1))  s1 = var1
    s2 = transfer_power_var
    if (present(var2))  s2 = var2
    if (present(power_ix)) ix =power_ix

    H=CP%h0/100._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    do ik=1, MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        if (kh==0) cycle
        k = kh*H

        dsig8 = MTrans%TransferData(s1,ik, &
            CP%Transfer%PK_redshifts_index(1:CP%Transfer%PK_num_redshifts))
        if (s1==s2) then
            dsig8 = dsig8**2
        else
            dsig8 = dsig8*MTrans%TransferData(s2,ik, &
                CP%Transfer%PK_redshifts_index(1:CP%Transfer%PK_num_redshifts))
        end if
        x= kh *R
        win =3*(sin(x)-x*cos(x))/x**3
        lnk=log(k)
        if (ik==1) then
            dlnk=0.5_dl
            !Approx for 2._dl/(CP%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
            !Contribution should be very small in any case
        else
            dlnk=lnk-lnko
        end if
        powers = ScalarPower(k,power_ix)
        dsig8=(win*k**2)**2*powers*dsig8
        sig8=sig8+(dsig8+dsig8o)*dlnk/2
        dsig8o=dsig8
        lnko=lnk
    end do

    if (present(root)) then
        if (root) sig8 =sqrt(sig8)
    else
        sig8 =sqrt(sig8)
    end if
    outvals(1:CP%Transfer%PK_num_redshifts) = sig8

    end subroutine Transfer_Get_SigmaR

    subroutine Transfer_GetSigmaRArray(MTrans, R, sigmaR, redshift_ix, var1, var2, power_ix)
    !Get array of SigmaR at (by default) redshift zero, for all values of R
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R(:)
    real(dl), intent(out) :: SigmaR(:)
    integer, intent(in), optional :: redshift_ix, var1, var2, power_ix
    integer red_ix, ik, subk
    real(dl) kh, k, h, dkh
    real(dl) lnk, dlnk, lnko, minR
    real(dl), dimension(size(R)) ::  x, win, dsig8, dsig8o, sig8, sig8o
    type(MatterPowerData) :: PKspline
    integer, parameter :: nsub = 5

    minR = minval(R)
    red_ix =  CP%Transfer%PK_redshifts_index(CP%Transfer%PK_num_redshifts)
    if (present(redshift_ix)) red_ix = redshift_ix

    call Transfer_GetMatterPowerData(MTrans, PKspline, power_ix, red_ix, var1, var2 )

    H=CP%h0/100._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    if (MTrans%TransferData(Transfer_kh,1,1)==0) call MpiStop('Transfer_GetSigmaRArray kh zero')
    do ik=1, MTrans%num_q_trans + 2
        if (ik < MTrans%num_q_trans) then
            dkh = (MTrans%TransferData(Transfer_kh,ik+1,1)- MTrans%TransferData(Transfer_kh,ik,1))/nsub
            !after last step just extrapolate a bit with previous size
        end if
        if (ik <= MTrans%num_q_trans) kh = MTrans%TransferData(Transfer_kh,ik,1)
        do subk = 1, nsub
            k = kh*H
            lnk=log(k)

            x= kh *R
            win =3*(sin(x)-x*cos(x))/x**3
            if (ik==1 .and. subk==1) then
                dlnk=0.5_dl
                !Approx for 2._dl/(CP%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
                !Contribution should be very small in any case
            else
                dlnk=lnk-lnko
            end if
            dsig8=win**2*(MatterPowerData_k(PKspline,  kh, 1)*k**3)
            sig8=sig8+(dsig8+dsig8o)*dlnk/2
            dsig8o=dsig8
            lnko=lnk
            kh = kh + dkh
        end do
    end do
    call MatterPowerdata_Free(PKspline)

    SigmaR=sqrt(sig8/(pi*twopi*h**3 ))

    end subroutine Transfer_GetSigmaRArray

    subroutine Transfer_Get_sigma8(MTrans, R, var1, var2)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc
    ! set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var1, var2
    integer ix
    real(dl) :: radius = 8._dl

    if (global_error_flag /= 0) return

    if (present(R)) radius = R

    do ix = 1, CP%InitPower%nn
        call Transfer_Get_SigmaR(MTrans, radius, MTrans%sigma_8(:,ix), var1,var2, ix)
    end do

    end subroutine Transfer_Get_sigma8

    subroutine Transfer_Get_sigmas(MTrans, R, var_delta, var_v)
    !Get sigma8 and sigma_{delta v} (for growth, like f sigma8 in LCDM)
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var_delta,var_v
    real(dl) :: radius = 8._dl
    integer s1, s2, ix

    if (global_error_flag /= 0) return

    if (present(R)) radius = R
    s1 = transfer_power_var
    if (present(var_delta))  s1 = var_delta
    s2 = Transfer_Newt_vel_cdm
    if (present(var_v))  s2 = var_v

    do ix = 1, CP%InitPower%nn
        call Transfer_Get_SigmaR(MTrans, radius, MTrans%sigma_8(:,ix), s1,s1, ix)
        if (get_growth_sigma8) call Transfer_Get_SigmaR(MTrans, radius, &
            MTrans%sigma2_vdelta_8(:,ix), s1, s2, ix, root=.false.)
    end do

    end subroutine Transfer_Get_sigmas

    subroutine Transfer_output_Sig8(MTrans)
    Type(MatterTransferData), intent(in) :: MTrans
    integer in, j
    !JD 08/13 Changes in here to PK arrays and variables
    integer j_PK

    do in=1, CP%InitPower%nn
        if (CP%InitPower%nn>1)  write(*,*) 'Power spectrum : ', in
        do j_PK=1, CP%Transfer%PK_num_redshifts
            j = CP%Transfer%PK_redshifts_index(j_PK)
            write(*,'("at z =",f7.3," sigma8 (all matter) = ",f7.4)') &
                CP%Transfer%redshifts(j), MTrans%sigma_8(j_PK,in)
        end do
        if (get_growth_sigma8) then
            do j_PK=1, CP%Transfer%PK_num_redshifts
                j = CP%Transfer%PK_redshifts_index(j_PK)
                write(*,'("at z =",f7.3," sigma8^2_vd/sigma8  = ",f7.4)') &
                    CP%Transfer%redshifts(j), MTrans%sigma2_vdelta_8(j_PK,in)/MTrans%sigma_8(j_PK,in)
            end do
        end if
    end do

    end subroutine Transfer_output_Sig8

    subroutine Transfer_Allocate(MTrans)
    Type(MatterTransferData) :: MTrans
    integer st

    deallocate(MTrans%q_trans, STAT = st)
    deallocate(MTrans%TransferData, STAT = st)
    deallocate(MTrans%sigma_8, STAT = st)
    if (get_growth_sigma8) deallocate(MTrans%sigma2_vdelta_8, STAT = st)
    allocate(MTrans%q_trans(MTrans%num_q_trans))
    allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans,CP%Transfer%num_redshifts))
    !JD 08/13 Changes in here to PK arrays and variables
    allocate(MTrans%sigma_8(CP%Transfer%PK_num_redshifts, CP%InitPower%nn))
    if (get_growth_sigma8) allocate(MTrans%sigma2_vdelta_8(CP%Transfer%PK_num_redshifts, CP%InitPower%nn))

    end  subroutine Transfer_Allocate

    subroutine Transfer_Nullify(Mtrans)
    Type(MatterTransferData):: MTrans

    Mtrans%num_q_trans = 0
    nullify(MTrans%q_trans)
    nullify(MTrans%TransferData)
    nullify(MTrans%sigma_8)
    nullify(MTrans%sigma2_vdelta_8)

    end subroutine Transfer_Nullify

    subroutine Transfer_Free(MTrans)
    Type(MatterTransferData):: MTrans
    integer st

    deallocate(MTrans%q_trans, STAT = st)
    deallocate(MTrans%TransferData, STAT = st)
    deallocate(MTrans%sigma_8, STAT = st)
    if (get_growth_sigma8) deallocate(MTrans%sigma2_vdelta_8, STAT = st)
    call Transfer_Nullify(MTrans)

    end subroutine Transfer_Free

    !JD 08/13 Changes for nonlinear lensing of CMB + MPK compatibility
    !Changed function below to write to only P%NLL_*redshifts* variables
    subroutine Transfer_SetForNonlinearLensing(P)
    Type(TransferParams) :: P
    integer i
    real maxRedshift

    P%kmax = max(P%kmax,5*AccuracyBoost)
    P%k_per_logint  = 0
    maxRedshift = 10
    P%NLL_num_redshifts =  nint(10*AccuracyBoost)
    if (HighAccuracyDefault .and. AccuracyBoost>=2) then
        !only notionally more accuracy, more stable for RS
        maxRedshift =15
    end if
    if (P%NLL_num_redshifts > max_transfer_redshifts) &
        call MpiStop('Transfer_SetForNonlinearLensing: Too many redshifts')
    do i=1,P%NLL_num_redshifts
        P%NLL_redshifts(i) = real(P%NLL_num_redshifts-i)/(P%NLL_num_redshifts/maxRedshift)
    end do

    end subroutine Transfer_SetForNonlinearLensing



    subroutine Transfer_SaveToFiles(MTrans,FileNames)
    use IniFile
    Type(MatterTransferData), intent(in) :: MTrans
    integer i,ik
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    !JD 08/13 Changes in here to PK arrays and variables
    integer i_PK
    integer unit

    do i_PK=1, CP%Transfer%PK_num_redshifts
        if (FileNames(i_PK) /= '') then
            i = CP%Transfer%PK_redshifts_index(i_PK)
            unit = open_file_header(FileNames(i_PK), 'k/h', transfer_name_tags, 14)
            do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                    write(unit,'(*(E15.6))') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
            end do
            close(unit)
        end if
    end do

    end subroutine Transfer_SaveToFiles

    subroutine Transfer_SaveMatterPower(MTrans, FileNames)
    use IniFile
    !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
    Type(MatterTransferData), intent(in) :: MTrans
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer itf,in,i
    integer points
    real, dimension(:,:,:), allocatable :: outpower
    real minkh,dlnkh
    Type(MatterPowerData) :: PK_data
    integer ncol
    !JD 08/13 Changes in here to PK arrays and variables
    integer itf_PK
    integer unit
    character(name_tag_len) :: columns(3)

    ncol=1
    if (CP%InitPower%nn>1 .and. output_file_headers) error stop 'InitPower%nn>1 deprecated'

    do itf=1, CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = CP%Transfer%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
                allocate(outpower(points,CP%InitPower%nn,ncol))

                do in = 1, CP%InitPower%nn
                    call Transfer_GetMatterPowerData(MTrans, PK_data, in, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)&
                        call MatterPowerdata_MakeNonlinear(PK_Data)

                    outpower(:,in,1) = exp(PK_data%matpower(:,1))
                    call MatterPowerdata_Free(PK_Data)
                end do
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData(Transfer_kh,i,1),outpower(i,1:CP%InitPower%nn,:)
                end do
                close(unit)
            else
                minkh = 1e-4
                dlnkh = 0.02
                points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
                !             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
                allocate(outpower(points,CP%InitPower%nn,1))
                do in = 1, CP%InitPower%nn
                    call Transfer_GetMatterPowerS(MTrans,outpower(1,in,1), itf, in, minkh,dlnkh, points)
                end do

                columns(1) = 'P'
                unit = open_file_header(FileNames(itf), 'k/h', columns(:1), 15)

                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1:CP%InitPower%nn,1)
                end do
                close(unit)
            end if

            deallocate(outpower)
        end if
    end do

    end subroutine Transfer_SaveMatterPower

    !JD 08/13 New function for nonlinear lensing of CMB + MPK compatibility
    !Build master redshift array from array of desired Nonlinear lensing (NLL)
    !redshifts and an array of desired Power spectrum (PK) redshifts.
    !At the same time fill arrays for NLL and PK that indicate indices
    !of their desired redshifts in the master redshift array.
    !Finally define number of redshifts in master array. This is usually given by:
    !P%num_redshifts = P%PK_num_redshifts + P%NLL_num_redshifts - 1.  The -1 comes
    !from the fact that z=0 is in both arrays (when non-linear is on)
    subroutine Transfer_SortAndIndexRedshifts(P)
    Type(TransferParams) :: P
    integer i, iPK, iNLL
    real(dl), parameter :: tol = 1.d-5

    i=0
    iPK=1
    iNLL=1
    do while (iPk<=P%PK_num_redshifts .or. iNLL<=P%NLL_num_redshifts)
        !JD write the next line like this to account for roundoff issues with ==. Preference given to PK_Redshift
        i=i+1
        if (i > max_transfer_redshifts) &
            call Mpistop('Transfer_SortAndIndexRedshifts: Too many redshifts')

        if(iNLL>P%NLL_num_redshifts .or. P%PK_redshifts(iPK)>P%NLL_redshifts(iNLL)+tol) then
            P%redshifts(i)=P%PK_redshifts(iPK)
            P%PK_redshifts_index(iPK)=i
            iPK=iPK+1
        else if(iPK>P%PK_num_redshifts .or. P%NLL_redshifts(iNLL)>P%PK_redshifts(iPK)+tol) then
            P%redshifts(i)=P%NLL_redshifts(iNLL)
            P%NLL_redshifts_index(iNLL)=i
            iNLL=iNLL+1
        else
            P%redshifts(i)=P%PK_redshifts(iPK)
            P%PK_redshifts_index(iPK)=i
            P%NLL_redshifts_index(iNLL)=i
            iPK=iPK+1
            iNLL=iNLL+1
        end if
    end do
    P%num_redshifts=i

    end subroutine Transfer_SortAndIndexRedshifts

    end module Transfer


    !ccccccccccccccccccccccccccccccccccccccccccccccccccc

    module ThermoData
    use ModelData
    implicit none
    private
    integer,parameter :: nthermo=20000

    real(dl) tb(nthermo),cs2(nthermo),xe(nthermo)
    real(dl) dcs2(nthermo)
    real(dl) dotmu(nthermo), ddotmu(nthermo)
    real(dl) sdotmu(nthermo),emmu(nthermo)
    real(dl) demmu(nthermo)
    real(dl) dddotmu(nthermo),ddddotmu(nthermo)
    real(dl) winlens(nthermo),dwinlens(nthermo), scalefactor(nthermo)
    real(dl) tauminn,dlntau,Maxtau
    logical, parameter :: dowinlens = .false.

    real(dl) :: tight_tau, actual_opt_depth
    !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
    real(dl) :: matter_verydom_tau
    real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.

    public thermo,inithermo, tight_tau, IonizationFunctionsAtTime, &
        Thermo_OpacityToTime,matter_verydom_tau, ThermoData_Free,&
        z_star, z_drag, GetBackgroundEvolution
    contains

    subroutine thermo(tau,cs2b,opacity, dopacity)
    !Compute unperturbed sound speed squared,
    !and ionization fraction by interpolating pre-computed tables.
    !If requested also get time derivative of opacity
    implicit none
    real(dl) tau,cs2b,opacity
    real(dl), intent(out), optional :: dopacity

    integer i
    real(dl) d

    d=log(tau/tauminn)/dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
        cs2b=cs2(1)+(d+i-1)*dcs2(1)
        opacity=dotmu(1)+(d-1)*ddotmu(1)
        call MpiStop('thermo out of bounds')
    else if (i >= nthermo) then
        cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
        opacity=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
        if (present(dopacity)) then
            dopacity = 0
            call MpiStop('thermo: shouldn''t happen')
        end if
    else
        !Cubic spline interpolation.
        cs2b=cs2(i)+d*(dcs2(i)+d*(3*(cs2(i+1)-cs2(i))  &
            -2*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)  &
            +2*(cs2(i)-cs2(i+1)))))
        opacity=dotmu(i)+d*(ddotmu(i)+d*(3*(dotmu(i+1)-dotmu(i)) &
            -2*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
            +2*(dotmu(i)-dotmu(i+1)))))

        if (present(dopacity)) then
            dopacity=(ddotmu(i)+d*(dddotmu(i)+d*(3*(ddotmu(i+1)  &
                -ddotmu(i))-2*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
                +dddotmu(i+1)+2*(ddotmu(i)-ddotmu(i+1))))))/(tau*dlntau)
        end if
    end if
    end subroutine thermo

    function Thermo_OpacityToTime(opacity)
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(dotmu(j)> opacity)
        j=j+1
    end do

    Thermo_OpacityToTime = exp((j-1)*dlntau)*tauminn

    end function Thermo_OpacityToTime

    subroutine inithermo(taumin,taumax)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use precision
    use ModelParams
    use MassiveNu
    real(dl) taumin,taumax


    real(dl) tau01,adot0,a0,a02,x1,x2,barssc,dtau
    real(dl) xe0,tau,a,a2
    real(dl) adot,tg0,ahalf,adothalf,fe,thomc,thomc0,etc,a2t
    real(dl) dtbdla,vfi,cf1,maxvis, vis
    integer ncount,i,j1,iv,ns
    real(dl) spline_data(nthermo)
    real(dl) last_dotmu
    real(dl) dtauda  !diff of tau w.CP%r.t a and integration
    external dtauda
    real(dl) a_verydom
    real(dl) awin_lens1p,awin_lens2p,dwing_lens, rs, DA
    real(dl) z_eq, a_eq
    real(dl) rombint
    integer noutput
    external rombint

    call Recombination_Init(CP%Recomb, CP%omegac, CP%omegab,CP%Omegan, CP%Omegav, &
        CP%h0,CP%tcmb,CP%yhe,CP%Num_Nu_massless + CP%Num_Nu_massive)
    !almost all the time spent here
    if (global_error_flag/=0) return
    Maxtau=taumax
    tight_tau = 0
    actual_opt_depth = 0
    ncount=0
    z_star=0.d0
    z_drag=0.d0
    thomc0= Compton_CT * CP%tcmb**4
    r_drag0 = 3.d0/4.d0*CP%omegab*grhom/grhog
    !thomc0=5.0577d-8*CP%tcmb**4

    tauminn=0.05d0*taumin
    dlntau=log(CP%tau0/tauminn)/(nthermo-1)
    last_dotmu = 0

    matter_verydom_tau = 0
    a_verydom = AccuracyBoost*5*(grhog+grhornomass)/(grhoc+grhob)

    !  Initial conditions: assume radiation-dominated universe.
    tau01=tauminn
    adot0=adotrad
    a0=adotrad*tauminn
    a02=a0*a0
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.
    tb(1)=CP%tcmb/a0
    xe0=1._dl
    x1=0._dl
    x2=1._dl
    xe(1)=xe0+0.25d0*CP%yhe/(1._dl-CP%yhe)*(x1+2*x2)
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(1))
    cs2(1)=4._dl/3._dl*barssc*tb(1)
    dotmu(1)=xe(1)*akthom/a02
    sdotmu(1)=0

    do i=2,nthermo
        tau=tauminn*exp((i-1)*dlntau)
        dtau=tau-tau01
        !  Integrate Friedmann equation using inverse trapezoidal rule.

        a=a0+adot0*dtau
        scaleFactor(i)=a
        a2=a*a

        adot=1/dtauda(a)

        if (matter_verydom_tau ==0 .and. a > a_verydom) then
            matter_verydom_tau = tau
        end if

        a=a0+2._dl*dtau/(1._dl/adot0+1._dl/adot)
        !  Baryon temperature evolution: adiabatic except for Thomson cooling.
        !  Use  quadrature solution.
        ! This is redundant as also calculated in REFCAST, but agrees well before reionization
        tg0=CP%tcmb/a0
        ahalf=0.5d0*(a0+a)
        adothalf=0.5d0*(adot0+adot)
        !  fe=number of free electrons divided by total number of free baryon
        !  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
        !  more accuracy is required (unlikely) then this can be iterated with
        !  the solution of the ionization equation.
        fe=(1._dl-CP%yhe)*xe(i-1)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i-1))
        thomc=thomc0*fe/adothalf/ahalf**3
        etc=exp(-thomc*(a-a0))
        a2t=a0*a0*(tb(i-1)-tg0)*etc-CP%tcmb/thomc*(1._dl-etc)
        tb(i)=CP%tcmb/a+a2t/(a*a)

        ! If there is re-ionization, smoothly increase xe to the
        ! requested value.
        if (CP%Reion%Reionization .and. tau > CP%ReionHist%tau_start) then
            if(ncount == 0) then
                ncount=i-1
            end if
            xe(i) = Reionization_xe(a, tau, xe(ncount))
            !print *,1/a-1,xe(i)
            if (CP%AccurateReionization .and. CP%DerivedParameters) then
                dotmu(i)=(Recombination_xe(a) - xe(i))*akthom/a2

                if (last_dotmu /=0) then
                    actual_opt_depth = actual_opt_depth - 2._dl*dtau/(1._dl/dotmu(i)+1._dl/last_dotmu)
                end if
                last_dotmu = dotmu(i)
            end if
        else
            xe(i)=Recombination_xe(a)
        end if

        !  Baryon sound speed squared (over c**2).
        dtbdla=-2._dl*tb(i)-thomc*adothalf/adot*(a*tb(i)-CP%tcmb)
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i))
        cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3._dl)


        ! Calculation of the visibility function
        dotmu(i)=xe(i)*akthom/a2

        if (tight_tau==0 .and. 1/(tau*dotmu(i)) > 0.005) tight_tau = tau !0.005
        !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)

        if (tau < 0.001) then
            sdotmu(i)=0
        else
            sdotmu(i)=sdotmu(i-1)+2._dl*dtau/(1._dl/dotmu(i)+1._dl/dotmu(i-1))
        end if

        a0=a
        tau01=tau
        adot0=adot
    end do !i

    if (CP%Reion%Reionization .and. (xe(nthermo) < 0.999d0)) then
        write(*,*)'Warning: xe at redshift zero is < 1'
        write(*,*) 'Check input parameters and Reionization_xe'
        write(*,*) 'function in the Reionization module'
    end if

    do j1=1,nthermo
        if (sdotmu(j1) - sdotmu(nthermo)< -69) then
            emmu(j1)=1.d-30
        else
            emmu(j1)=exp(sdotmu(j1)-sdotmu(nthermo))
            if (.not. CP%AccurateReionization .and. &
                actual_opt_depth==0 .and. xe(j1) < 1e-3) then
                actual_opt_depth = -sdotmu(j1)+sdotmu(nthermo)
            end if
            if (CP%AccurateReionization .and. CP%DerivedParameters .and. z_star==0.d0) then
                if (sdotmu(nthermo)-sdotmu(j1) - actual_opt_depth < 1) then
                    tau01=1-(sdotmu(nthermo)-sdotmu(j1) - actual_opt_depth)
                    tau01=tau01*(1._dl/dotmu(j1)+1._dl/dotmu(j1-1))/2
                    z_star = 1/(scaleFactor(j1)- tau01/dtauda(scaleFactor(j1))) -1
                end if
            end if
        end if
    end do

    if (CP%AccurateReionization .and. FeedbackLevel > 0 .and. CP%DerivedParameters) then
        write(*,'("Reion opt depth      = ",f7.4)') actual_opt_depth
    end if


    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time of maximum visibility
    if (ncount == 0) then
        cf1=1._dl
        ns=nthermo
    else
        cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
        ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
        vis = emmu(j1)*dotmu(j1)
        tau = tauminn*exp((j1-1)*dlntau)
        vfi=vfi+vis*cf1*dlntau*tau
        if ((iv == 0).and.(vfi > 1.0d-7/AccuracyBoost)) then
            taurst=9._dl/10._dl*tau
            iv=1
        elseif (iv == 1) then
            if (vis > maxvis) then
                maxvis=vis
                tau_maxvis = tau
            end if
            if (vfi > 0.995) then
                taurend=tau
                iv=2
                exit
            end if
        end if
    end do

    if (iv /= 2) then
        call GlobalError('inithermo: failed to find end of recombination',error_reionization)
        return
    end if

    if (dowinlens) then
        vfi=0
        awin_lens1p=0
        awin_lens2p=0
        winlens=0
        do j1=1,nthermo-1
            vis = emmu(j1)*dotmu(j1)
            tau = tauminn*exp((j1-1)*dlntau)
            vfi=vfi+vis*cf1*dlntau*tau
            if (vfi < 0.995) then
                dwing_lens =  vis*cf1*dlntau*tau / 0.995

                awin_lens1p = awin_lens1p + dwing_lens
                awin_lens2p = awin_lens2p + dwing_lens/(CP%tau0-tau)
            end if
            winlens(j1)= awin_lens1p/(CP%tau0-tau) - awin_lens2p
        end do
    end if

    ! Calculating the timesteps during recombination.

    if (CP%WantTensors) then
        dtaurec=min(dtaurec,taurst/160)/AccuracyBoost
    else
        dtaurec=min(dtaurec,taurst/40)/AccuracyBoost
        if (do_bispectrum .and. hard_bispectrum) dtaurec = dtaurec / 4
    end if

    if (CP%Reion%Reionization) taurend=min(taurend,CP%ReionHist%tau_start)

    if (DebugMsgs) then
        write (*,*) 'taurst, taurend = ', taurst, taurend
    end if

    call splini(spline_data,nthermo)
    call splder(cs2,dcs2,nthermo,spline_data)
    call splder(dotmu,ddotmu,nthermo,spline_data)
    call splder(ddotmu,dddotmu,nthermo,spline_data)
    call splder(dddotmu,ddddotmu,nthermo,spline_data)
    call splder(emmu,demmu,nthermo,spline_data)
    if (dowinlens) call splder(winlens,dwinlens,nthermo,spline_data)

    call SetTimeSteps

    if ((CP%want_zstar .or. CP%DerivedParameters) .and. z_star==0.d0) call find_z(optdepth,z_star)
    if (CP%want_zdrag .or. CP%DerivedParameters) call find_z(dragoptdepth,z_drag)

    if (CP%DerivedParameters) then
        rs =rombint(dsound_da_exact,1d-8,1/(z_star+1),1d-6)
        DA = AngularDiameterDistance(z_star)/(1/(z_star+1))

        ThermoDerivedParams( derived_Age ) = DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
        ThermoDerivedParams( derived_zstar ) = z_star
        ThermoDerivedParams( derived_rstar ) = rs
        ThermoDerivedParams( derived_thetastar ) = 100*rs/DA
        ThermoDerivedParams( derived_DAstar ) = DA/1000
        ThermoDerivedParams( derived_zdrag ) = z_drag
        rs =rombint(dsound_da_exact,1d-8,1/(z_drag+1),1d-6)
        ThermoDerivedParams( derived_rdrag ) = rs
        ThermoDerivedParams( derived_kD ) =  sqrt(1.d0/(rombint(ddamping_da, 1d-8, 1/(z_star+1), 1d-6)/6))
        ThermoDerivedParams( derived_thetaD ) =  100*pi/ThermoDerivedParams( derived_kD )/DA
        z_eq = (grhob+grhoc)/(grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates))) -1
        ThermoDerivedParams( derived_zEQ ) = z_eq
        a_eq = 1/(1+z_eq)
        ThermoDerivedParams( derived_kEQ ) = 1/(a_eq*dtauda(a_eq))
        ThermoDerivedParams( derived_thetaEQ ) = 100*timeOfz(z_eq)/DA
        ThermoDerivedParams( derived_theta_rs_EQ ) = 100*rombint(dsound_da_exact,1d-8,a_eq,1d-6)/DA

        if (associated(BackgroundOutputs%z_outputs)) then
            if (allocated(BackgroundOutputs%H)) &
                deallocate(BackgroundOutputs%H, BackgroundOutputs%DA, BackgroundOutputs%rs_by_D_v)
            noutput = size(BackgroundOutputs%z_outputs)
            allocate(BackgroundOutputs%H(noutput), BackgroundOutputs%DA(noutput), BackgroundOutputs%rs_by_D_v(noutput))
            do i=1,noutput
                BackgroundOutputs%H(i) = HofZ(BackgroundOutputs%z_outputs(i))
                BackgroundOutputs%DA(i) = AngularDiameterDistance(BackgroundOutputs%z_outputs(i))
                BackgroundOutputs%rs_by_D_v(i) = rs/BAO_D_v_from_DA_H(BackgroundOutputs%z_outputs(i), &
                    BackgroundOutputs%DA(i),BackgroundOutputs%H(i))
            end do
        end if

        if (FeedbackLevel > 0) then
            write(*,'("Age of universe/GYr  = ",f7.3)') ThermoDerivedParams( derived_Age )
            write(*,'("zstar                = ",f8.2)') ThermoDerivedParams( derived_zstar )
            write(*,'("r_s(zstar)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rstar )
            write(*,'("100*theta            = ",f9.6)') ThermoDerivedParams( derived_thetastar )
            write(*,'("DA(zstar)/Gpc        = ",f9.5)') ThermoDerivedParams( derived_DAstar )

            write(*,'("zdrag                = ",f8.2)') ThermoDerivedParams( derived_zdrag )
            write(*,'("r_s(zdrag)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rdrag )

            write(*,'("k_D(zstar) Mpc       = ",f7.4)') ThermoDerivedParams( derived_kD )
            write(*,'("100*theta_D          = ",f9.6)') ThermoDerivedParams( derived_thetaD )

            write(*,'("z_EQ (if v_nu=1)     = ",f8.2)') ThermoDerivedParams( derived_zEQ )
            write(*,'("k_EQ Mpc (if v_nu=1) = ",f9.6)') ThermoDerivedParams( derived_kEQ )
            write(*,'("100*theta_EQ         = ",f9.6)') ThermoDerivedParams( derived_thetaEQ )
            write(*,'("100*theta_rs_EQ      = ",f9.6)') ThermoDerivedParams( derived_theta_rs_EQ )

        end if
    end if

    end subroutine inithermo


    subroutine SetTimeSteps
    real(dl) dtau0
    integer nri0, nstep

    call Ranges_Init(TimeSteps)

    call Ranges_Add_delta(TimeSteps, taurst, taurend, dtaurec)

    ! Calculating the timesteps after recombination
    if (CP%WantTensors) then
        dtau0=max(taurst/40,Maxtau/2000._dl/AccuracyBoost)
    else
        dtau0=Maxtau/500._dl/AccuracyBoost
        if (do_bispectrum) dtau0 = dtau0/3
        !Don't need this since adding in Limber on small scales
        !  if (CP%DoLensing) dtau0=dtau0/2
        !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
    end if

    call Ranges_Add_delta(TimeSteps,taurend, CP%tau0, dtau0)

    if (CP%Reion%Reionization) then
        nri0=int(Reionization_timesteps(CP%ReionHist)*AccuracyBoost)
        !Steps while reionization going from zero to maximum
        call Ranges_Add(TimeSteps,CP%ReionHist%tau_start,CP%ReionHist%tau_complete,nri0)
    end if

    !Create arrays out of the region information.
    call Ranges_GetArray(TimeSteps)
    nstep = TimeSteps%npoints

    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Set ',nstep, ' time steps'

    end subroutine SetTimeSteps


    subroutine ThermoData_Free

    call Ranges_Free(TimeSteps)

    end subroutine ThermoData_Free


    subroutine IonizationFunctionsAtTime(tau, opac, dopac, ddopac, &
        vis, dvis, ddvis, expmmu, lenswin)
    real(dl), intent(in) :: tau
    real(dl), intent(out):: opac, dopac, ddopac, vis, dvis, ddvis, expmmu, lenswin
    real(dl) d
    integer i

    d=log(tau/tauminn)/dlntau+1._dl
    i=int(d)
    d=d-i

    if (i < nthermo) then
        opac=dotmu(i)+d*(ddotmu(i)+d*(3._dl*(dotmu(i+1)-dotmu(i)) &
            -2._dl*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
            +2._dl*(dotmu(i)-dotmu(i+1)))))
        dopac=(ddotmu(i)+d*(dddotmu(i)+d*(3._dl*(ddotmu(i+1)  &
            -ddotmu(i))-2._dl*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
            +dddotmu(i+1)+2._dl*(ddotmu(i)-ddotmu(i+1))))))/(tau &
            *dlntau)
        ddopac=(dddotmu(i)+d*(ddddotmu(i)+d*(3._dl*(dddotmu(i+1) &
            -dddotmu(i))-2._dl*ddddotmu(i)-ddddotmu(i+1)  &
            +d*(ddddotmu(i)+ddddotmu(i+1)+2._dl*(dddotmu(i) &
            -dddotmu(i+1)))))-(dlntau**2)*tau*dopac) &
            /(tau*dlntau)**2
        expmmu=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
            -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
            +2._dl*(emmu(i)-emmu(i+1)))))

        if (dowinlens) then
            lenswin=winlens(i)+d*(dwinlens(i)+d*(3._dl*(winlens(i+1)-winlens(i)) &
                -2._dl*dwinlens(i)-dwinlens(i+1)+d*(dwinlens(i)+dwinlens(i+1) &
                +2._dl*(winlens(i)-winlens(i+1)))))
        end if
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3*opac*dopac+ddopac)
    else
        opac=dotmu(nthermo)
        dopac=ddotmu(nthermo)
        ddopac=dddotmu(nthermo)
        expmmu=emmu(nthermo)
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3._dl*opac*dopac+ddopac)
    end if

    end subroutine IonizationFunctionsAtTime


    function ddamping_da(a)
    real(dl) :: ddamping_da
    real(dl), intent(in) :: a
    real(dl) :: R
    real(dl) :: dtauda
    external dtauda

    R=r_drag0*a
    !ignoring reionisation, not relevant for distance measures
    ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(a)*a**2/(Recombination_xe(a)*akthom)

    end function ddamping_da


    !!!!!!!!!!!!!!!!!!!
    !JH: functions and subroutines for calculating z_star and z_drag

    function doptdepth_dz(z)
    real(dl) :: doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a
    real(dl) :: dtauda
    external dtauda

    a = 1._dl/(1._dl+z)

    !ignoring reionisation, not relevant for distance measures
    doptdepth_dz = Recombination_xe(a)*akthom*dtauda(a)

    end function doptdepth_dz

    function optdepth(z)
    real(dl) :: rombint2
    external rombint2
    real(dl) optdepth
    real(dl),intent(in) :: z

    optdepth = rombint2(doptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function optdepth


    function ddragoptdepth_dz(z)
    real(dl) :: ddragoptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a
    real(dl) :: dtauda
    external dtauda

    a = 1._dl/(1._dl+z)
    ddragoptdepth_dz = doptdepth_dz(z)/r_drag0/a

    end function ddragoptdepth_dz


    function dragoptdepth(z)
    real(dl) :: rombint2
    external rombint2
    real(dl) dragoptdepth
    real(dl),intent(in) :: z

    dragoptdepth =  rombint2(ddragoptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function dragoptdepth


    subroutine find_z(func,zout)  !find redshift at which (photon/drag) optical depth = 1
    real(dl), external :: func
    real(dl), intent(out) :: zout
    real(dl) :: try1,try2,diff,avg
    integer :: i

    try1 = 0.d0
    try2 = 10000.d0

    i=0
    diff = 10.d0
    do while (diff .gt. 1d-3)
        i=i+1
        if (i .eq. 100) then
            call GlobalError('optical depth redshift finder did not converge',error_reionization)
            zout=0
            return
        end if

        diff = func(try2)-func(try1)
        avg = 0.5d0*(try2+try1)
        if (func(avg) .gt. 1.d0) then
            try2 = avg
        else
            try1 = avg
        end if
    end do

    zout = avg

    end subroutine find_z

    !!!!!!!!!!!!!!!!!!! end JH

    subroutine GetBackgroundEvolution(ntimes, times, outputs)
    integer, intent(in) :: ntimes
    real(dl), intent(in) :: times(ntimes)
    real(dl) :: outputs(5, ntimes)
    real(dl) spline_data(nthermo), ddxe(nthermo), ddTb(nthermo)
    real(dl) :: d, tau, cs2b, opacity, vis, Tbaryon
    integer i, ix

    call splini(spline_data,nthermo)
    call splder(xe,ddxe,nthermo,spline_data)
    call splder(Tb,ddTb,nthermo,spline_data)

    outputs = 0
    do ix = 1, ntimes
        tau = times(ix)
        if (tau < tauminn) cycle
        d=log(tau/tauminn)/dlntau+1._dl
        i=int(d)
        d=d-i
        call thermo(tau,cs2b, opacity)

        if (i < nthermo) then
            outputs(1,ix)=xe(i)+d*(ddxe(i)+d*(3._dl*(xe(i+1)-xe(i)) &
                -2._dl*ddxe(i)-ddxe(i+1)+d*(ddxe(i)+ddxe(i+1) &
                +2._dl*(xe(i)-xe(i+1)))))
            vis=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
                -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
                +2._dl*(emmu(i)-emmu(i+1)))))
            Tbaryon = tb(i)+d*(ddtb(i)+d*(3._dl*(tb(i+1)-tb(i)) &
                -2._dl*ddtb(i)-ddtb(i+1)+d*(ddtb(i)+ddtb(i+1) &
                +2._dl*(tb(i)-tb(i+1)))))
        else
            outputs(1,ix)=xe(nthermo)
            vis = emmu(nthermo)
            Tbaryon = Tb(nthermo)
        end if

        outputs(2, ix) = opacity
        outputs(3, ix) = opacity*vis
        outputs(4, ix) = cs2b
        outputs(5, ix) = Tbaryon
    end do

    end subroutine GetBackgroundEvolution

    end module ThermoData
