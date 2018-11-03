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

    module RedshiftSpaceData
    use precision
    implicit none

    integer, parameter :: window_21cm = 1, window_counts = 2, window_lensing = 3

    Type TRedWin
        integer kind
        real(dl) Redshift
        real(dl) sigma
        real(dl) a
        real(dl) bias, dlog10Ndm
        real(dl) tau, tau_start, tau_end
        real(dl) sigma_z ! estimate of of width in z
        real(dl) sigma_tau !width in conformal time (set by code)
        real(dl) chi0, chimin
        integer mag_index !The index into the extra sources using for adding magnification to counts
        real(dl), dimension(:), allocatable :: winF, wing,wing2,wingtau,dwing,dwing2,dwingtau,ddwing,ddwing2,ddwingtau,&
            winV,dwinV,ddwinV, win_lens, comoving_density_ev
        real(dl) Fq, optical_depth_21
        logical has_lensing_window
    end Type TRedWin


    integer, parameter :: max_redshiftwindows  = 20
    Type(TRedWin), target :: Redshift_W(max_redshiftwindows)

    Type SourceTermParams
        logical :: limber_windows = .false.
        logical :: do_counts_lensing = .false.
        logical :: line_phot_dipole = .false.
        logical :: line_phot_quadrupole= .false.
        logical :: line_basic = .true.
        logical :: line_extra = .false.
        logical :: line_distortions = .true.
        logical :: line_reionization = .false.
        logical :: counts_velocity  = .true.
        logical :: counts_radial = .false. !does not include time delay; subset of counts_velocity, just 1/(chi*H) term
        logical :: counts_density = .true.
        logical :: counts_redshift = .true.
        logical :: counts_timedelay = .true. !time delay terms * 1/(H*chi)
        logical :: counts_ISW = .true.
        logical :: counts_potential = .true. !terms in potentials at source
        logical :: counts_evolve = .false.
        logical :: use_21cm_mK = .true.
    end type SourceTermParams

    contains

    function count_obs_window_z(Win, z, winamp)
    !distribution function W(z) for the observed sources, used for lensing and number count spectrum
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    !note this is the total count distribution observed, not a fractional selection function on an underlying distribution
    Type(TRedWin) :: Win
    real(dl), intent(in) :: z
    real(dl) count_obs_window_z, dz,winamp
    real(dl), parameter :: root2pi = 2.506628274_dl

    dz = z-Win%Redshift
    winamp =  exp(-(dz/Win%sigma)**2/2)
    count_obs_window_z =winamp/Win%sigma/root2pi
    !!!
    !         winamp =  z**3*exp(-(z/3.35821)**13)/20
    !         count_obs_window_z = winamp*20/28.49261877

    end function count_obs_window_z

    function counts_background_z(Win, z)
    !if counts_evolve = T this function is used to get n(z) for the source population
    !(not the same as the number actually observed)
    Type(TRedWin) :: Win
    real(dl), intent(in) :: z
    real(dl) counts_background_z, winamp

    counts_background_z = count_obs_window_z(Win, z, winamp)
    !This is the special case where you observe all of the sources

    end function counts_background_z

    function Window_f_a(Win, a, winamp)
    !distribution function as function of scale factor a, used for 21cm
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    Type(TRedWin) :: Win
    real(dl), intent(in) :: a
    real(dl) Window_f_a, winamp
    real(dl), parameter :: root2pi = 2.506628274_dl
    !   real(dl) z

    if (Win%kind == window_21cm) then
        !Take W_T(S) = W_f(S)/S to be Gaussain, for frequency integration of T_b
        winamp =  exp(-((a-Win%a)/Win%sigma)**2/2)
        Window_f_a = a*winamp/Win%sigma/root2pi
    else
        Window_f_a = count_obs_window_z(Win, 1/a-1,winamp)/a**2
    end if

    end function Window_f_a



    subroutine InitRedshiftWindow(Win)
    Type(TRedWin) :: Win

    Win%bias = 1
    Win%mag_index = 0

    end subroutine InitRedshiftWindow

    end module RedshiftSpaceData

    module CambSettings
    use precision
    use Errors
    use classes
    use MassiveNu
    use InitialPower
    use Reionization
    use Recombination
    use DarkEnergyInterface
    use RedshiftSpaceData
    use constants, only : const_pi, const_twopi, const_fourpi, const_eightpi, COBE_CMBTemp, default_nnu
    use MiscUtils
    use RangeUtils
    use StringUtils
    use MpiUtils, only : MpiStop
    implicit none
    public

    character(LEN=*), parameter :: version = 'Oct18'

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
    integer, parameter :: outNone=1

    integer :: max_bessels_l_index  = 1000000
    real(dl) :: max_bessels_etak = 1000000*2

    real(dl), parameter ::  OutputDenominator = const_twopi
    !When using outNone the output is l(l+1)Cl/OutputDenominator

    type TransferParams
        logical     ::  high_precision = .false.
        logical     ::  accurate_massive_neutrinos = .false.
        integer     ::  num_redshifts = 1
        real(dl)    ::  kmax = 0.9_dl        !these are acutally q values, but same as k for flat
        integer     ::  k_per_logint =0
        real(dl)    ::  redshifts(max_transfer_redshifts)
        !JD 08/13 Added so both NL lensing and PK can be run at the same time
        real(dl)    ::  PK_redshifts(max_transfer_redshifts)
        real(dl)    ::  NLL_redshifts(max_transfer_redshifts)
        integer     ::  PK_redshifts_index(max_transfer_redshifts)
        integer     ::  NLL_redshifts_index(max_transfer_redshifts)
        integer     ::  PK_num_redshifts = 1
        integer     ::  NLL_num_redshifts = 0
    end type TransferParams

    type AccuracyParams
        !Parameters for checking/changing overall accuracy
        !parameters equal to 1 corresponds to ~0.1% scalar C_l accuracy (at L>600)

        real(dl) :: lSampleBoost=1._dl
        !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation

        real(dl) :: AccuracyBoost =1._dl

        !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
        !Can also be used to improve speed significantly if less accuracy is required.
        !or improving accuracy for extreme models.
        !Note this does not increase lSamp%l sampling or massive neutrino q-sampling

        real(dl) :: lAccuracyBoost=1._dl
        !Boost number of multipoles integrated in Boltzman heirarchy

        logical   :: AccuratePolarization = .true.
        !Do you care about the accuracy of the polarization Cls?

        logical   :: AccurateBB = .false.
        !Do you care about BB accuracy (e.g. in lensing)

        !Reionization settings - used if Reion%Reionization=.true.
        logical   :: AccurateReionization = .true.
        !Do you care about pecent level accuracy on EE signal from reionization?

        !The following allow separate tweaking (all also affected by AccuracyBoost above)

        real(dl) :: TimeStepBoost = 1._dl !sampling of timesteps during recombination

        real(dl) :: IntTolBoost = 1._dl !Tolerances for integrating differential equations

        real(dl) :: SourcekAccuracyBoost = 1._dl !Accuracy of k sampling for source time integration

        real(dl) :: IntkAccuracyBoost = 1._dl !Accuracy of k sampling for integration

        real(dl) :: TransferkBoost = 1._dl !Accuracy of k sampling for transfer functions

        real(dl) :: NonFlatIntAccuracyBoost = 1._dl !Accuracy of non-flat time integration

        real(dl) :: BessIntBoost = 1._dl !Accuracy of bessel integration truncation

        real(dl) :: BesselBoost = 1._dl !accuracy of bessel pre-computation sampling

        real(dl) :: LimberBoost = 1._dl !Accuracy of Limber approximation use

        real(dl) :: Kmax_Boost = 1._dl !Boost max k for source window functions

        real(dl) :: neutrino_q_boost = 1._dl !number of momenta integrated for neutrino perturbations

        real(dl) :: thermo_boost = 1._dl !number of time steps for background thermal history interpolation

    end type AccuracyParams

    !other variables, options, derived variables, etc.

    integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2
    integer, parameter :: NonLinear_both=3  !JD 08/13 added so both can be done

    ! Main parameters type
    type CAMBparams
        logical   :: WantCls  = .true.
        logical   :: WantTransfer = .false.
        logical   :: WantScalars = .true.
        logical   :: WantTensors = .false.
        logical   :: WantVectors = .false.
        logical   :: want_cl_2D_array = .true.
        logical   :: DoLensing  = .true.
        logical   :: Do21cm = .false.
        logical   :: want_zstar = .false.
        logical   :: want_zdrag = .false.     !!JH for updated BAO likelihood.
        logical   :: PK_WantTransfer = .false. !JD 08/13 Added so both NL lensing and PK can be run at the same time
        integer   :: NonLinear = NonLinear_none

        logical   :: Want_CMB = .true.
        logical   :: Want_CMB_lensing = .true.

        integer   :: Max_l = 2500
        integer   :: Max_l_tensor = 600
        real(dl)  :: Max_eta_k = 5000
        real(dl)  :: Max_eta_k_tensor = 1200
        ! _tensor settings only used in initialization,
        !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

        real(dl)  :: omegab = 0.045_dl
        real(dl)  :: omegac = 0.255_dl
        real(dl)  :: omegav = 0.7_dl
        real(dl)  :: omegan = 0._dl
        !Omega baryon, CDM, Lambda and massive neutrino
        real(dl)  :: H0 = 65._dl
        real(dl)  :: TCMB = COBE_CMBTemp
        real(dl)  :: Yhe = 0.24_dl
        real(dl)  :: Num_Nu_massless = default_nnu
        integer   :: Num_Nu_massive = 0 !sum of Nu_mass_numbers below
        integer   :: Nu_mass_eigenstates = 0  !1 for degenerate masses
        logical   :: share_delta_neff = .false. !take fractional part to heat all eigenstates the same
        real(dl)  :: Nu_mass_degeneracies(max_nu)
        real(dl)  :: Nu_mass_fractions(max_nu) !The ratios of the total densities
        integer   :: Nu_mass_numbers(max_nu) !physical number per eigenstate

        integer   :: Scalar_initial_condition = 1 !adiabatic
        !must be one of the initial_xxx values defined in GaugeInterface

        integer   :: OutputNormalization = outNone
        !outNone, or C_OutputNormalization=1 if > 1

        real(dl)  :: Alens = 1._dl

        integer   :: MassiveNuMethod = Nu_best

        type(ReionizationParams) :: Reion
        type(RecombinationParams):: Recomb
        type(TransferParams)     :: Transfer
        type(AccuracyParams)     :: Accuracy
        type(SourceTermParams)   :: SourceTerms

        logical :: DoLateRadTruncation = .true.
        !if true, use smooth approx to radition perturbations after decoupling on
        !small scales, saving evolution of irrelevant osciallatory multipole equations

        logical :: Evolve_baryon_cs = .false.
        !if true, evolves equation for Delta_{T_m} to get cs_2 = \delta p /\delta\rho for perfect gas

        logical :: Evolve_delta_xe = .false. !Include ionization fraction perturbations

        logical :: Evolve_delta_Ts =.false. !Equilibrium result agrees to sub-percent level

        real(dl) ::  InitialConditionVector(1:10) !Allow up to 10 for future extensions
        !ignored unless Scalar_initial_condition == initial_vector

        logical :: OnlyTransfers = .false. !Don't use initial power spectrum data, instead get Delta_q_l array
        !If true, sigma_8 is not calculated either

        !Sources
        logical :: transfer_21cm_cl = .false.

        logical :: DerivedParameters = .true.!calculate various derived parameters  (ThermoDerivedParams)

        !Sources
        logical :: Log_lvalues  = .false.

        class(TInitialPower), allocatable :: InitPower
        class(TDarkEnergyBase), allocatable :: DarkEnergy
        class(TNonLinearModel), allocatable :: NonLinearModel

    end type CAMBparams


    Type TBackgroundOutputs
        real(dl), pointer :: z_outputs(:) => null()
        real(dl), allocatable :: H(:), DA(:), rs_by_D_v(:)
    end Type TBackgroundOutputs

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13
    integer, parameter :: nthermo_derived = 13

    Type lSamples
        integer :: nl = 0
        integer, allocatable :: l(:)
    contains
    procedure :: Init => lSamples_init
    procedure :: IndexOf => lSamples_indexOf
    procedure :: InterpolateClArr
    procedure :: InterpolateClArrTemplated
    end Type lSamples

    logical, parameter :: dowinlens = .false. !not used, test getting CMB lensing using visibility
    integer, parameter :: thermal_history_def_timesteps = 40000

    Type TThermoData
        !Background thermal history, interpolated from precomputed tables
        integer :: nthermo !Number of table steps
        !baryon temperature, sound speed, ionization fractions, and e^(-tau)
        real(dl), dimension(:), allocatable :: tb, cs2, xe, emmu
        !derivatives for interpolation
        real(dl), dimension(:), allocatable :: dcs2, dotmu, ddotmu, sdotmu
        real(dl), dimension(:), allocatable :: demmu, dddotmu, ddddotmu
        real(dl), dimension(:), allocatable :: winlens, dwinlens
        real(dl) tauminn,dlntau,Maxtau
        real(dl) :: tight_tau, actual_opt_depth
        !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
        real(dl) :: recombination_saha_tau
        !sound horizon and recombination redshifts
        real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.
        real(dl), dimension(:), allocatable :: step_redshift, rhos_fac, drhos_fac
        real(dl) :: tau_start_redshiftwindows,tau_end_redshiftwindows
        logical :: has_lensing_windows = .false.
        real(dl) recombination_Tgas_tau
        !Mapping between redshift and time
        real(dl), private,dimension(:), allocatable :: redshift_time, dredshift_time
        real(dl), private, dimension(:), allocatable :: arhos_fac, darhos_fac, ddarhos_fac
    contains
    procedure :: Init => Thermo_Init
    procedure :: OpacityToTime => Thermo_OpacityToTime
    procedure :: values => Thermo_values
    procedure :: IonizationFunctionsAtTime
    procedure, private :: DoWindowSpline
    procedure, private :: SetTimeSteps
    procedure, private :: SetTimeStepWindows
    end type TThermoData

    !Sources
    Type CalWins
        real(dl), allocatable :: awin_lens(:),  dawin_lens(:)
    end Type CalWins

    Type LimberRec
        integer n1,n2 !corresponding time step array indices
        real(dl), dimension(:), pointer :: k  => NULL()
        real(dl), dimension(:), pointer :: Source  => NULL()
    end Type LimberRec

    Type ClTransferData
        !Cl transfer function variables
        !values of q for integration over q to get C_ls
        Type (lSamples) :: ls ! l that are computed
        integer :: NumSources
        !Changes -scalars:  2 for just CMB, 3 for lensing
        !- tensors: T and E and phi (for lensing), and T, E, B respectively

        type (TRanges) :: q
        real(dl), dimension(:,:,:), pointer :: Delta_p_l_k => NULL()

        !The L index of the lowest L to use for Limber
        integer, dimension(:), pointer :: Limber_l_min => NULL()
        !For each l, the set of k in each limber window
        !indices LimberWindow(SourceNum,l)
        Type(LimberRec), dimension(:,:), pointer :: Limber_windows => NULL()

        !The maximum L needed for non-Limber
        integer max_index_nonlimber

    end Type ClTransferData

    !Computed output power spectra data

    integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4, C_PhiTemp = 5, C_PhiE=6
    integer :: C_last = C_PhiE
    integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4
    integer, parameter :: name_tag_len = 12
    character(LEN=name_tag_len), dimension(C_PhiE), parameter :: C_name_tags = ['TT','EE','TE','PP','TP','EP']
    character(LEN=name_tag_len), dimension(CT_Cross), parameter :: CT_name_tags = ['TT','EE','BB','TE']
    character(LEN=name_tag_len), dimension(7), parameter :: lens_pot_name_tags = ['TT','EE','BB','TE','PP','TP','EP']

    type TCLdata
        Type(ClTransferData) :: CTransScal, CTransTens, CTransVec

        real(dl), dimension (:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , Cl_type)
        !where Cl_type is one of the above constants

        real(dl), dimension (:,:,:), allocatable :: Cl_Scalar_Array
        !Indices are Cl_xxx( l , field1,field2)
        !where ordering of fields is T, E, \psi (CMB lensing potential), window_1, window_2...

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this
        real(dl) , dimension (:,:), allocatable :: Cl_lensed
        !Cl_lensed(l, Cl_type) are the interpolated Cls
    contains
    procedure :: InitCls => TCLdata_InitCls
    procedure :: output_cl_files => TCLdata_output_cl_files
    procedure :: output_lens_pot_files => TCLdata_output_lens_pot_files
    procedure :: NormalizeClsAtL => TCLdata_NormalizeClsAtL
    procedure :: output_veccl_files => TCLdata_output_veccl_files
    end type TCLdata

    type CAMBstate

        type(CAMBparams) :: CP
        !Derived parameters, not set initially
        type(ReionizationHistory) :: ReionHist
        type(RecombinationData) :: Recombination

        logical flat,closed,open
        real(dl) omegak
        real(dl) curv, curvature_radius, Ksign !CP%r = 1/sqrt(|curv|), Ksign = 1,0 or -1
        real(dl) tau0,chi0 !time today and rofChi(CP%tau0/CP%r)
        real(dl) scale !relative to CP%flat. e.g. for scaling lSamp%l sampling.

        logical ::call_again = .false.
        !if being called again with same parameters to get different thing

        !     grhom =kappa*a^2*rho_m0
        !     grhornomass=grhor*number of massless neutrino species
        !     taurst,taurend - time at start/end of recombination
        !     dtaurec - dtau during recombination
        !     adotrad - a(tau) in radiation era

        real(dl) grhom,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
        real(dl) taurst,dtaurec,taurend,tau_maxvis,adotrad

        !Neutrinos
        real(dl) grhormass(max_nu)
        !     nu_masses=m_nu*c**2/(k_B*T_nu0)
        real(dl) nu_masses(max_nu)

        real(dl) akthom !sigma_T * (number density of protons now)
        real(dl) fHe !n_He_tot / n_H_tot
        real(dl) Nnow

        real(dl) ThermoDerivedParams(nthermo_derived)
        Type(TNuPerturbations) :: NuPerturbations

        Type(TBackgroundOutputs) :: BackgroundOutputs

        !Time steps for sampling sources
        Type(TRanges) :: TimeSteps
        !Background interpolation tables for thermal history etc.
        Type(TThermoData) :: ThermoData

        !Matter transfer data
        Type (MatterTransferData):: MT

        Type(TClData) :: CLdata

    contains
    procedure :: CAMBParams_Set
    procedure :: rofChi
    procedure :: cosfunc
    procedure :: tanfunc
    procedure :: invsinfunc
    end type CAMBstate

    type(CAMBstate), pointer :: State !Current state
    type(CAMBParams), pointer :: CP   !Current parameters (part of state)

    integer :: ThreadNum = 0
    !If zero assigned automatically, obviously only used if parallelised

    integer :: limber_phiphi = 0 !for l>limber_phiphi use limber approx for lensing potential
    integer :: num_redshiftwindows = 0
    integer :: num_extra_redshiftwindows = 0
    integer :: num_custom_sources = 0
    integer, allocatable :: custom_source_ell_scales(:)

    integer, parameter :: lmin = 2
    !must be either 1 or 2

    real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code

    real(dl), parameter :: tol=1.0d-4 !Base tolerance for integrations

    !     used as parameter for spline - tells it to use 'natural' end values
    real(dl), parameter :: spl_large=1.e40_dl

    character(LEN=1024) :: highL_unlensed_cl_template = 'HighLExtrapTemplate_lenspotentialCls.dat'
    !fiducial high-accuracy high-L C_L used for making small cosmology-independent numerical corrections
    !to lensing and C_L interpolation. Ideally close to models of interest, but dependence is weak.
    logical :: use_spline_template = .true.
    integer, parameter :: lmax_extrap_highl = 8000
    real(dl), allocatable :: highL_CL_template(:,:)

    real(dl), private, external :: dtauda, rombint, rombint2

    contains

    subroutine SetActiveState(P)
    class(CAMBstate), target :: P

    select type(P)
    type is (CAMBstate)
        State => P
        CP => P%CP
        class default
        error stop 'Unknown state type'
    end select
    end subroutine SetActiveState

    subroutine CAMBParams_Set(this, P, error, DoReion, call_again)
    !Initialize background variables; does not yet calculate thermal history
    use constants
    class(CAMBstate), target :: this
    type(CAMBparams), intent(in) :: P
    real(dl) GetOmegak, fractional_number, conv
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical, optional :: call_again
    logical WantReion, calling_again
    integer nu_i,actual_massless
    real(dl) nu_massless_degeneracy, neff_i, eta_k
    external GetOmegak
    !Constants in SI units

    global_error_flag = 0

    if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
        call GlobalError( 'Cannot generate tensor C_l and transfer without scalar C_l',error_unsupported_params)
    end if

    if (.not. allocated(P%DarkEnergy)) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if

    if (present(error)) error = global_error_flag
    if (global_error_flag/=0) return

    WantReion = DefaultTrue(DoReion)
    calling_again= DefaultFalse(call_again)

    if (calling_again) then
        this%CP%DerivedParameters = .false.
        this%CP%Accuracy = P%Accuracy
        this%CP%Transfer%high_precision = P%Transfer%high_precision
        this%CP%Reion%Reionization = P%Reion%Reionization
        this%CP%WantTransfer =P%WantTransfer
        this%CP%WantScalars =P%WantScalars
        this%CP%WantTensors =P%WantTensors
        this%CP%WantVectors =P%WantVectors
        this%CP%WantCls = P%WantCls
    else
        this%CP=P
        call SetActiveState(this)
        this%CP%Max_eta_k = max(this%CP%Max_eta_k,this%CP%Max_eta_k_tensor)
    end if

    if (P%WantTransfer) then
        this%CP%WantScalars=.true.
        if (.not. P%WantCls) then
            this%CP%Accuracy%AccuratePolarization = .false.
            !Sources
            this%CP%Reion%Reionization = this%CP%transfer_21cm_cl
        end if
    end if
    if (this%CP%WantTransfer.and. this%CP%MassiveNuMethod==Nu_approx) then
        this%CP%MassiveNuMethod = Nu_trunc
    end if

    if (.not. calling_again) then

        if (this%CP%Num_Nu_Massive /= sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))) then
            if (sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))/=0) &
                call MpiStop('Num_Nu_Massive is not sum of Nu_mass_numbers')
        end if
        if (this%CP%Omegan == 0 .and. this%CP%Num_Nu_Massive /=0) then
            if (this%CP%share_delta_neff) then
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + this%CP%Num_Nu_Massive
            else
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            end if
            this%CP%Num_Nu_Massive  = 0
            this%CP%Nu_mass_numbers = 0
        end if

        nu_massless_degeneracy = this%CP%Num_Nu_massless !N_eff for massless neutrinos
        if (this%CP%Num_nu_massive > 0) then
            if (this%CP%Nu_mass_eigenstates==0) call MpiStop('Have Num_nu_massive>0 but no nu_mass_eigenstates')
            if (this%CP%Nu_mass_eigenstates==1 .and. this%CP%Nu_mass_numbers(1)==0) &
                this%CP%Nu_mass_numbers(1) = this%CP%Num_Nu_Massive
            if (all(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)==0)) this%CP%Nu_mass_numbers=1 !just assume one for all
            if (this%CP%share_delta_neff) then
                !default case of equal heating of all neutrinos
                fractional_number = this%CP%Num_Nu_massless + this%CP%Num_Nu_massive
                actual_massless = int(this%CP%Num_Nu_massless + 1e-6_dl)
                neff_i = fractional_number/(actual_massless + this%CP%Num_Nu_massive)
                nu_massless_degeneracy = neff_i*actual_massless
                this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates) = &
                    this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)*neff_i
            end if
            if (abs(sum(this%CP%Nu_mass_fractions(1:this%CP%Nu_mass_eigenstates))-1) > 1e-4) &
                call MpiStop('Nu_mass_fractions do not add up to 1')
        else
            this%CP%Nu_mass_eigenstates = 0
        end if

        this%omegak = GetOmegak(this%CP)

        this%flat = (abs(this%omegak) <= OmegaKFlat)
        this%closed = this%omegak < -OmegaKFlat

        this%open = .not. this%flat .and. .not. this%closed
        if (this%flat) then
            this%curv=0
            this%Ksign=0
            this%curvature_radius=1._dl !so we can use tau/curvature_radius, etc, where r's cancel
        else
            this%curv=-this%omegak/((c/1000)/this%CP%h0)**2
            this%Ksign =sign(1._dl,this%curv)
            this%curvature_radius=1._dl/sqrt(abs(this%curv))
        end if
        !  grho gives the contribution to the expansion rate from: (g) photons,
        !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
        !  (m) nonrelativistic matter (for Omega=1).  grho is actually
        !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
        !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
        !  (Used only to set the initial conformal time.)

        !H0 is in km/s/Mpc

        this%grhom = 3*this%CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

        !grhom=3.3379d-11*h0*h0
        this%grhog = kappa/c**2*4*sigma_boltz/c**3*this%CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
        ! grhog=1.4952d-13*tcmb**4
        this%grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*this%grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        !grhor=3.3957d-14*tcmb**4

        !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
        !for massive Nu_mass_degeneracies parameters account for heating from grhor

        this%grhornomass=this%grhor*nu_massless_degeneracy
        this%grhormass=0
        do nu_i = 1, this%CP%Nu_mass_eigenstates
            this%grhormass(nu_i)=this%grhor*this%CP%Nu_mass_degeneracies(nu_i)
        end do
        this%grhoc=this%grhom*this%CP%omegac
        this%grhob=this%grhom*this%CP%omegab
        this%grhov=this%grhom*this%CP%omegav
        this%grhok=this%grhom*this%omegak
        !  adotrad gives the relation a(tau) in the radiation era:
        this%adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates)))/3)

        this%Nnow = this%CP%omegab*(1-this%CP%yhe)*this%grhom*c**2/kappa/m_H/Mpc**2

        this%akthom = sigma_thomson*this%Nnow*Mpc
        !sigma_T * (number density of protons now)

        this%fHe = this%CP%YHe/(mass_ratio_He_H*(1.d0-this%CP%YHe))  !n_He_tot / n_H_tot

        if (this%CP%omegan/=0) then
            !Initialize things for massive neutrinos
            call ThermalNuBackground%Init()
            call this%NuPerturbations%Init(P%Accuracy%AccuracyBoost*P%Accuracy%neutrino_q_boost)
            !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
            !  Get number density n of neutrinos from
            !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
            !  then m = Omega_nu/N_nu rho_crit /n
            !  Error due to velocity < 1e-5
            do nu_i=1, this%CP%Nu_mass_eigenstates
                this%nu_masses(nu_i)=fermi_dirac_const/(1.5d0*zeta3)*this%grhom/this%grhor* &
                    this%CP%omegan*this%CP%Nu_mass_fractions(nu_i)/this%CP%Nu_mass_degeneracies(nu_i)
            end do
        else
            this%nu_masses = 0
        end if
        call this%CP%DarkEnergy%Init()
        if (global_error_flag==0) then
            this%tau0=TimeOfz(this,0._dl)
            if (WantReion) call this%ReionHist%Init(this%CP%Reion,this%CP%YHe, this%akthom, this%tau0, FeedbackLevel)
            this%chi0=this%rofChi(this%tau0/this%curvature_radius)
            this%scale= this%chi0*this%curvature_radius/this%tau0  !e.g. change l sampling depending on approx peak spacing
            if (this%closed .and. this%tau0/this%curvature_radius >3.14) then
                call GlobalError('chi >= pi in closed model not supported',error_unsupported_params)
            end if
        end if
    end if

    !Sources
    if (this%CP%WantScalars .and. this%CP%WantCls .and. num_redshiftwindows>0) then
        eta_k = this%CP%Max_eta_k
        do nu_i = 1,num_redshiftwindows
            Redshift_w(nu_i)%tau = TimeOfz(this,Redshift_w(nu_i)%Redshift)
            Redshift_w(nu_i)%chi0 = this%tau0-Redshift_w(nu_i)%tau
            Redshift_w(nu_i)%chimin = min(Redshift_w(nu_i)%chi0,&
                this%tau0 - TimeOfz(this,max(0.05_dl,Redshift_w(nu_i)%Redshift - 1.5*Redshift_w(nu_i)%sigma)) )
            this%CP%Max_eta_k = max(this%CP%Max_eta_k, this%tau0*WindowKmaxForL(Redshift_w(nu_i),this%CP%max_l))
        end do
        if (eta_k /= this%CP%Max_eta_k .and. FeedbackLevel>0) &
            write (*,*) 'source max_eta_k: ', this%CP%Max_eta_k,'kmax = ', this%CP%Max_eta_k/this%tau0
    end if

    !JD 08/13 Changes for nonlinear lensing of CMB + MPK compatibility
    !if ( this%CP%NonLinear==NonLinear_Lens) then
    if (this%CP%NonLinear==NonLinear_Lens .or. this%CP%NonLinear==NonLinear_both ) then
        this%CP%Transfer%kmax = max(this%CP%Transfer%kmax, this%CP%Max_eta_k/this%tau0)
        if (FeedbackLevel > 0 .and. this%CP%Transfer%kmax== this%CP%Max_eta_k/this%tau0) &
            write (*,*) 'max_eta_k changed to ', this%CP%Max_eta_k
    end if

    if (global_error_flag/=0) then
        if (present(error)) error = global_error_flag
        return
    end if

    if (present(error)) then
        error = 0
    else if (FeedbackLevel > 0 .and. .not. calling_again) then
        write(*,'("Om_b h^2             = ",f9.6)') P%omegab*(P%H0/100)**2
        write(*,'("Om_c h^2             = ",f9.6)') P%omegac*(P%H0/100)**2
        write(*,'("Om_nu h^2            = ",f9.6)') P%omegan*(P%H0/100)**2
        write(*,'("Om_Lambda            = ",f9.6)') P%omegav
        write(*,'("Om_K                 = ",f9.6)') this%omegak
        write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-this%omegak-P%omegav
        write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta(this)
        if (this%CP%Num_Nu_Massive > 0) then
            write(*,'("N_eff (total)        = ",f9.6)') nu_massless_degeneracy + &
                sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            do nu_i=1, this%CP%Nu_mass_eigenstates
                conv = k_B*(8*this%grhor/this%grhog/7)**0.25*this%CP%tcmb/eV * &
                    (this%CP%nu_mass_degeneracies(nu_i)/this%CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                write(*,'(I2, " nu, g=",f7.4," m_nu*c^2/k_B/T_nu0= ",f9.2," (m_nu= ",f6.3," eV)")') &
                    this%CP%nu_mass_numbers(nu_i), this%CP%nu_mass_degeneracies(nu_i), &
                    this%nu_masses(nu_i),conv*this%nu_masses(nu_i)
            end do
        end if
    end if

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


    function rofChi(this,Chi) !sinh(chi) for open, sin(chi) for closed.
    class(CAMBstate) :: this
    real(dl) Chi,rofChi

    if (this%closed) then
        rofChi=sin(chi)
    else if (this%open) then
        rofChi=sinh(chi)
    else
        rofChi=chi
    endif
    end function rofChi


    function cosfunc (this,Chi)
    class(CAMBstate) :: this
    real(dl) Chi,cosfunc

    if (this%closed) then
        cosfunc= cos(chi)
    else if (this%open) then
        cosfunc=cosh(chi)
    else
        cosfunc = 1._dl
    endif
    end function cosfunc

    function tanfunc(this,Chi)
    class(CAMBstate) :: this
    real(dl) Chi,tanfunc
    if (this%closed) then
        tanfunc=tan(Chi)
    else if (this%open) then
        tanfunc=tanh(Chi)
    else
        tanfunc=Chi
    end if

    end function tanfunc

    function invsinfunc(this,x)
    class(CAMBstate) :: this
    real(dl) invsinfunc,x

    if (this%closed) then
        invsinfunc=asin(x)
    else if (this%open) then
        invsinfunc=log((x+sqrt(1._dl+x**2)))
    else
        invsinfunc = x
    endif
    end function invsinfunc

    function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x
    f_K = State%curvature_radius*State%rofChi(x/State%curvature_radius)

    end function f_K


    function DeltaTime(this, a1,a2, in_tol)
    Type(CAMBstate) :: this
    real(dl) DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol

    call SetActiveState(this)
    atol = PresentDefault(tol/1000/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1), in_tol)
    DeltaTime = rombint(dtauda,a1,a2,atol)

    end function DeltaTime

    function TimeOfz(this, z)
    Type(CAMBstate) :: this
    real(dl) TimeOfz
    real(dl), intent(IN) :: z

    TimeOfz=DeltaTime(this,0._dl,1._dl/(z+1._dl))
    end function TimeOfz

    subroutine TimeOfzArr(this, nz, redshifts, outputs)
    Type(CAMBstate) :: this
    integer, intent(in) :: nz
    real(dl), intent(in) :: redshifts(nz)
    real(dl), intent(out) :: outputs(nz)
    integer i

    !Dumb slow version
    !$OMP PARALLEL DO DEFAUlT(SHARED)
    do i=1, nz
        outputs(i) = timeOfZ(this,redshifts(i))
    end do
    !$OMP END PARALLEL DO

    end subroutine TimeOfzArr

    function DeltaPhysicalTimeGyr(this, a1,a2, in_tol)
    use constants
    Type(CAMBstate) :: this
    real(dl), intent(in) :: a1, a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) DeltaPhysicalTimeGyr, atol

    call SetActiveState(this)
    atol = PresentDefault(1d-4/exp(CP%Accuracy%AccuracyBoost-1), in_tol)
    DeltaPhysicalTimeGyr = rombint(dtda,a1,a2,atol)*Mpc/c/Gyr
    end function DeltaPhysicalTimeGyr

    function AngularDiameterDistance(this,z)
    Type(CAMBstate) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc
    real(dl) AngularDiameterDistance
    real(dl), intent(in) :: z

    AngularDiameterDistance = this%curvature_radius/(1+z)* &
        this%rofchi(ComovingRadialDistance(this,z) /this%curvature_radius)

    end function AngularDiameterDistance

    subroutine AngularDiameterDistanceArr(this, arr, z, n)
    Type(CAMBstate) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc for array of z
    !z array must be monotonically increasing
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    call ComovingRadialDistanceArr(this,arr, z, n, 1e-4_dl)
    if (this%flat) then
        arr = arr/(1+z)
    else
        do i=1, n
            arr(i) =  this%curvature_radius/(1+z(i))*this%rofchi(arr(i)/this%curvature_radius)
        end do
    end if

    end subroutine AngularDiameterDistanceArr


    function AngularDiameterDistance2(this,z1, z2) ! z1 < z2
    !From http://www.slac.stanford.edu/~amantz/work/fgas14/#cosmomc
    Type(CAMBstate) :: this
    real(dl) AngularDiameterDistance2
    real(dl), intent(in) :: z1, z2

    AngularDiameterDistance2 = this%curvature_radius/(1+z2)* &
        this%rofchi(ComovingRadialDistance(this,z2)/this%curvature_radius &
        - ComovingRadialDistance(this,z1)/this%curvature_radius)

    end function AngularDiameterDistance2

    function LuminosityDistance(this,z)
    Type(CAMBstate) :: this
    real(dl) LuminosityDistance
    real(dl), intent(in) :: z

    LuminosityDistance = AngularDiameterDistance(this,z)*(1+z)**2

    end function LuminosityDistance

    function ComovingRadialDistance(this, z)
    Type(CAMBstate) :: this
    real(dl) ComovingRadialDistance
    real(dl), intent(in) :: z

    ComovingRadialDistance = DeltaTime(this,1/(1+z),1._dl)

    end function ComovingRadialDistance

    subroutine ComovingRadialDistanceArr(this, arr, z, n, tol)
    !z array must be monotonically increasing
    Type(CAMBstate) :: this
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
                arr(i) = DeltaTime(this,1/(1+z(i)),1._dl, tol)
            end if
        else
            if (z(i) < z(i-1)) error stop 'ComovingRadialDistanceArr redshifts out of order'
            arr(i) = DeltaTime(this,1/(1+z(i)),1/(1+z(i-1)),tol)
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine ComovingRadialDistanceArr

    function Hofz(this,z)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    Type(CAMBstate) :: this
    real(dl) Hofz, a
    real(dl), intent(in) :: z

    call SetActiveState(this)
    a = 1/(1+z)
    Hofz = 1/(a**2*dtauda(a))

    end function Hofz

    subroutine HofzArr(this, arr, z, n)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    Type(CAMBstate) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i
    real(dl) :: a

    call SetActiveState(this)
    do i=1, n
        a = 1/(1+z(i))
        arr(i) = 1/(a**2*dtauda(a))
    end do

    end subroutine HofzArr

    real(dl) function BAO_D_v_from_DA_H(z, DA, Hz)
    real(dl), intent(in) :: z, DA, Hz
    real(dl) ADD

    ADD = DA*(1.d0+z)
    BAO_D_v_from_DA_H = ((ADD)**2.d0*z/Hz)**(1.d0/3.d0)

    end function BAO_D_v_from_DA_H

    real(dl) function BAO_D_v(this,z)
    Type(CAMBstate) :: this
    real(dl), intent(IN) :: z

    BAO_D_v = BAO_D_v_from_DA_H(z,AngularDiameterDistance(this,z), Hofz(this,z))

    end function BAO_D_v

    function dsound_da_exact(a)
    implicit none
    real(dl) dsound_da_exact,a,R,cs

    R = 3*state%grhob*a / (4*state%grhog)
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_exact=dtauda(a)*cs

    end function dsound_da_exact

    function dsound_da(a)
    !approximate form used e.g. by CosmoMC for theta
    implicit none
    real(dl) dsound_da,a,R,cs

    R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
    !          R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and used for previous consistency
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da=dtauda(a)*cs

    end function dsound_da

    function dtda(a)
    real(dl) dtda,a

    dtda= dtauda(a)*a
    end function

    function CosmomcTheta(this)
    Type(CAMBstate) :: this
    real(dl) zstar, astar, atol, rs, DA
    real(dl) CosmomcTheta
    real(dl) ombh2, omdmh2

    call SetActiveState(this)
    ombh2 = CP%omegab*(CP%h0/100.0d0)**2
    omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2

    !!From Hu & Sugiyama
    zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
        (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
        (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

    astar = 1/(1+zstar)
    atol = 1e-6
    rs = rombint(dsound_da,1d-8,astar,atol)
    DA = AngularDiameterDistance(this,zstar)/astar
    CosmomcTheta = rs/DA
    !       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

    end function CosmomcTheta

    !Sources
    function WindowKmaxForL(W,ell) result(res)
    Type(TRedWin), intent(in) :: W
    real(dl) res
    integer, intent(in)::  ell

    if (W%kind == window_lensing) then
        res = CP%Accuracy%AccuracyBoost*18*ell/W%chi0
    else
        !On large scales power can be aliased from smaller, so make sure k goes up until at least the turnover
        !in the matter power spectrum
        res = CP%Accuracy%AccuracyBoost*max(0.05_dl,2.5*ell/W%chimin)
    end if

    res = res* CP%Accuracy%Kmax_Boost
    end function WindowKmaxForL


    function lSamples_indexOf(lSet,l)
    class(lSamples) :: lSet
    integer, intent(in) :: l
    integer lSamples_indexOf, i

    do i=2,lSet%nl
        if (l < lSet%l(i)) then
            lSamples_indexOf = i-1
            return
        end if
    end do
    lSamples_indexOf = lSet%nl

    end function  lSamples_indexOf

    subroutine lSamples_init(this, State, max_l)
    ! This subroutines initializes lSet%l arrays. Other values will be interpolated.
    class(lSamples) :: this
    class(CAMBstate), target :: State
    integer, intent(IN) :: max_l
    integer lind, lvar, step, top, bot
    integer, allocatable :: ls(:)
    real(dl) AScale

    allocate(ls(max_l))
    if (allocated(this%l)) deallocate(this%l)
    associate(Accuracy => State%CP%Accuracy)
        Ascale=State%scale/Accuracy%lSampleBoost

        if (Accuracy%lSampleBoost >=50) then
            !just do all of them
            lind=0
            do lvar=lmin, max_l
                lind=lind+1
                ls(lind)=lvar
            end do
            this%nl=lind
            allocate(this%l, source = ls(1:lind))
            return
        end if

        lind=0
        do lvar=lmin, 10
            lind=lind+1
            ls(lind)=lvar
        end do

        if (Accuracy%AccurateReionization) then
            if (Accuracy%lSampleBoost > 1) then
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
            if (Accuracy%lSampleBoost >1) then
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
        if (State%CP%Log_lvalues) then
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
                    if (.not. use_spline_template) then
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
                if (.not. State%flat) ls(lind-1)=int(max_l+ls(lind-2))/2
                !Not in CP%flat case so interpolation table is the same when using lower l_max
            end if
        end if
    end associate
    this%nl=lind
    allocate(this%l, source=ls(1:lind))

    end subroutine lSamples_init

    subroutine InterpolateClArr(lSet,iCl, all_Cl, max_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lmin:*)
    integer, intent(in), optional :: max_index
    integer il,llo,lhi, xi
    real(dl) ddCl(lSet%nl)
    real(dl) xl(lSet%nl)
    real(dl) a0,b0,ho
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    integer max_ind

    max_ind = PresentDefault(lSet%nl, max_index)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArr')

    xl = real(lSet%l(1:lSet%nl),dl)
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
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lmin:*)
    integer, intent(in) :: max_ind
    integer, intent(in), optional :: template_index
    integer maxdelta, il
    real(dl) DeltaCL(lSet%nl)
    real(dl), allocatable :: tmpall(:)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArrTemplated')
    if (use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
            !interpolate only the difference between the C_l and an accurately interpolated template. Temp only for the mo.
            !Using unlensed for template, seems to be good enough
            maxdelta=max_ind
            do while (lSet%l(maxdelta) > lmax_extrap_highl)
                maxdelta=maxdelta-1
            end do
            DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

            call lSet%InterpolateClArr(DeltaCl, all_Cl, maxdelta)

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


    subroutine Thermo_values(this,tau, cs2b, opacity, dopacity)
    !Compute unperturbed sound speed squared,
    !and ionization fraction by interpolating pre-computed tables.
    !If requested also get time derivative of opacity
    class(TThermoData) :: this
    real(dl) :: tau, cs2b, opacity
    real(dl), intent(out), optional :: dopacity
    integer i
    real(dl) d

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
        cs2b=this%cs2(1)+(d+i-1)*this%dcs2(1)
        opacity=this%dotmu(1)+(d-1)*this%ddotmu(1)
        call MpiStop('thermo out of bounds')
    else if (i >= this%nthermo) then
        cs2b=this%cs2(this%nthermo)+(d+i-this%nthermo)*this%dcs2(this%nthermo)
        opacity=this%dotmu(this%nthermo)+(d-this%nthermo)*this%ddotmu(this%nthermo)
        if (present(dopacity)) then
            dopacity = 0
            call MpiStop('thermo: shouldn''t happen')
        end if
    else
        !Cubic spline interpolation.
        cs2b=this%cs2(i)+d*(this%dcs2(i)+d*(3*(this%cs2(i+1)-this%cs2(i))  &
            -2*this%dcs2(i)-this%dcs2(i+1)+d*(this%dcs2(i)+this%dcs2(i+1)  &
            +2*(this%cs2(i)-this%cs2(i+1)))))
        opacity=this%dotmu(i)+d*(this%ddotmu(i)+d*(3*(this%dotmu(i+1)-this%dotmu(i)) &
            -2*this%ddotmu(i)-this%ddotmu(i+1)+d*(this%ddotmu(i)+this%ddotmu(i+1) &
            +2*(this%dotmu(i)-this%dotmu(i+1)))))

        if (present(dopacity)) then
            dopacity=(this%ddotmu(i)+d*(this%dddotmu(i)+d*(3*(this%ddotmu(i+1)  &
                -this%ddotmu(i))-2*this%dddotmu(i)-this%dddotmu(i+1)+d*(this%dddotmu(i) &
                +this%dddotmu(i+1)+2*(this%ddotmu(i)-this%ddotmu(i+1))))))/(tau*this%dlntau)
        end if
    end if
    end subroutine Thermo_values

    function Thermo_OpacityToTime(this,opacity)
    class(TThermoData) :: this
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(this%dotmu(j)> opacity)
        j=j+1
    end do

    Thermo_OpacityToTime = exp((j-1)*this%dlntau)*this%tauminn

    end function Thermo_OpacityToTime

    subroutine Thermo_Init(this,taumin)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use StringUtils
    class(TThermoData) :: this
    real(dl), intent(in) :: taumin
    integer nthermo
    real(dl) tau01,adot0,a0,a02,x1,x2,barssc,dtau
    real(dl) xe0,tau,a,a2
    real(dl) adot,tg0,ahalf,adothalf,fe,thomc,thomc0,etc,a2t
    real(dl) dtbdla,vfi,cf1,maxvis, vis
    integer ncount,i,j1,iv,ns
    real(dl), allocatable :: spline_data(:), scaleFactor(:)
    real(dl) last_dotmu
    real(dl) a_verydom
    real(dl) awin_lens1p,awin_lens2p,dwing_lens, rs, DA
    real(dl) z_eq, a_eq
    integer noutput

    Type(CalWins), dimension(:), allocatable, target :: RW
    real(dl) awin_lens1(num_redshiftwindows),awin_lens2(num_redshiftwindows)
    real(dl) Tspin, Trad, rho_fac, window, tau_eps
    integer transfer_ix(max_transfer_redshifts)
    integer RW_i, j2
    real(dl) Tb21cm, winamp, z
    character(len=:), allocatable :: outstr

    nthermo = nint(thermal_history_def_timesteps*CP%Accuracy%thermo_boost)
    this%nthermo = nthermo
    call State%Recombination%Init(CP%Recomb, CP%omegac, CP%omegab,CP%Omegan, CP%Omegav, &
        CP%h0,CP%tcmb,CP%yhe,CP%Num_Nu_massless + CP%Num_Nu_massive, &
        CP%Evolve_baryon_cs .or. CP%Evolve_delta_xe .or. CP%Do21cm)
    !almost all the time spent here
    if (global_error_flag/=0) return

    if (allocated(this%tb) .and. this%nthermo/=size(this%tb)) then
        deallocate(this%tb, this%cs2, this%xe, this%emmu)
        deallocate(this%dcs2, this%dotmu, this%ddotmu, this%sdotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    endif
    if (.not. allocated(this%tb)) then
        allocate(this%tb(nthermo), this%cs2(nthermo), this%xe(nthermo), this%emmu(nthermo))
        allocate(this%dcs2(nthermo), this%dotmu(nthermo), this%ddotmu(nthermo), this%sdotmu(nthermo))
        allocate(this%demmu(nthermo), this%dddotmu(nthermo), this%ddddotmu(nthermo))
        if (dowinlens) allocate(this%winlens(nthermo), this%dwinlens(nthermo))
    end if

    if (num_redshiftwindows >0) then
        allocate(this%redshift_time(nthermo),this%dredshift_time(nthermo))
        allocate(this%arhos_fac(nthermo), this%darhos_fac(nthermo), this%ddarhos_fac(nthermo))
        allocate(RW(num_redshiftwindows))
    end if
    allocate(scaleFactor(nthermo), spline_data(nthermo))

    !Sources
    if (CP%Evolve_delta_xe) this%recombination_saha_tau  = TimeOfZ(State,State%Recombination%recombination_saha_z)
    if (CP%Evolve_baryon_cs .or. CP%Evolve_delta_xe .or. CP%Evolve_delta_Ts .or. CP%Do21cm) &
        this%recombination_Tgas_tau = TimeOfz(State,1/Do21cm_mina-1)
    transfer_ix =0

    this%Maxtau=State%tau0
    this%tight_tau = 0
    this%actual_opt_depth = 0
    ncount=0
    this%z_star=0.d0
    this%z_drag=0.d0
    thomc0= Compton_CT * CP%tcmb**4
    this%r_drag0 = 3.d0/4.d0*CP%omegab*State%grhom/State%grhog
    !thomc0=5.0577d-8*CP%tcmb**4

    this%tauminn=0.05d0*taumin
    this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)
    last_dotmu = 0

    this%matter_verydom_tau = 0
    a_verydom = CP%Accuracy%AccuracyBoost*5*(State%grhog+State%grhornomass)/(State%grhoc+State%grhob)

    !  Initial conditions: assume radiation-dominated universe.
    tau01=this%tauminn
    adot0=State%adotrad
    a0=State%adotrad*this%tauminn
    a02=a0*a0
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.
    this%tb(1)=CP%tcmb/a0
    xe0=1._dl
    x1=0._dl
    x2=1._dl
    this%xe(1)=xe0+0.25d0*CP%yhe/(1._dl-CP%yhe)*(x1+2*x2)
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(1))
    this%cs2(1)=4._dl/3._dl*barssc*this%tb(1)
    this%dotmu(1)=this%xe(1)*State%akthom/a02
    this%sdotmu(1)=0

    !Sources
    awin_lens1=0
    awin_lens2=0

    do RW_i = 1, num_redshiftwindows
        associate (RedWin => Redshift_w(RW_i))
            RedWin%tau_start = 0
            RedWin%tau_end = this%Maxtau
            if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts .and. CP%SourceTerms%do_counts_lensing) then
                allocate(RW(RW_i)%awin_lens(nthermo))
                allocate(RW(RW_i)%dawin_lens(nthermo))
            end if
        end associate
    end do

    do i=2,nthermo
        tau=this%tauminn*exp((i-1)*this%dlntau)
        dtau=tau-tau01
        !  Integrate Friedmann equation using inverse trapezoidal rule.

        a=a0+adot0*dtau
        scaleFactor(i)=a
        a2=a*a

        adot=1/dtauda(a)
        z= 1._dl/a-1._dl

        if (num_redshiftwindows>0) then
            this%redshift_time(i) = z
        end if
        if (a > 1d-4 .and. num_redshiftwindows>0) then
            if (CP%Do21cm ) then
                Tspin = State%Recombination%T_s(a)
                Trad = CP%TCMB/a
                rho_fac = line21_const*State%NNow/a**3
                tau_eps = a*line21_const*State%NNow/a**3/(adot/a)/Tspin/1000
                this%arhos_fac(i) = (1-exp(-tau_eps))/tau_eps*a*rho_fac*(1 - Trad/Tspin)/(adot/a)
                !         arhos_fac(i) = a*rho_fac*(1 - Trad/Tspin)/(adot/a)
            else
                rho_fac = State%grhoc/(a2*a)
                this%arhos_fac(i) = a*rho_fac/(adot/a)
            end if
        end if


        do RW_i = 1, num_redshiftwindows
            associate (Win => RW(RW_i), RedWin => Redshift_w(RW_i))
                if (a > 1d-4) then
                    window = Window_f_a(RedWin,a, winamp)

                    if  (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts  &
                        .and. CP%SourceTerms%do_counts_lensing) then
                        if (State%tau0 - tau > 2) then
                            dwing_lens =  adot * window *dtau
                            awin_lens1(RW_i) = awin_lens1(RW_i) + dwing_lens
                            awin_lens2(RW_i) = awin_lens2(RW_i) + dwing_lens/(State%tau0-tau)
                            Win%awin_lens(i) = awin_lens1(RW_i)/(State%tau0-tau) - awin_lens2(RW_i)
                        else
                            Win%awin_lens(i) = 0
                        end if
                    end if

                    if (RedWin%tau_start ==0 .and. winamp > 1e-8) then
                        RedWin%tau_start = tau01
                    else if (RedWin%tau_start /=0 .and. RedWin%tau_end==This%MaxTau .and. winamp < 1e-8) then
                        RedWin%tau_end = min(State%tau0,tau + dtau)
                        if (DebugMsgs) print *,'time window = ', RedWin%tau_start, RedWin%tau_end
                    end if
                else
                    if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts &
                        .and. CP%SourceTerms%do_counts_lensing) then
                        Win%awin_lens(i)=0
                    end if
                end if
            end associate
        end do
        !End sources

        if (this%matter_verydom_tau ==0 .and. a > a_verydom) then
            this%matter_verydom_tau = tau
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
        fe=(1._dl-CP%yhe)*this%xe(i-1)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i-1))
        thomc=thomc0*fe/adothalf/ahalf**3
        etc=exp(-thomc*(a-a0))
        a2t=a0*a0*(this%tb(i-1)-tg0)*etc-CP%tcmb/thomc*(1._dl-etc)
        this%tb(i)=CP%tcmb/a+a2t/(a*a)
        !Sources
        if (CP%WantTransfer) then
            do RW_i = 1, CP%Transfer%num_redshifts
                if (z< CP%Transfer%redshifts(RW_i) .and. transfer_ix(RW_i)==0) then
                    transfer_ix(RW_i) = i
                end if
            end do
        end if

        ! If there is re-ionization, smoothly increase xe to the
        ! requested value.
        if (CP%Reion%Reionization .and. tau > State%ReionHist%tau_start) then
            if(ncount == 0) then
                ncount=i-1
            end if
            this%xe(i) = State%ReionHist%x_e(a, tau, this%xe(ncount))
            !print *,1/a-1,xe(i)
            if (CP%Accuracy%AccurateReionization .and. CP%DerivedParameters) then
                this%dotmu(i)=(State%Recombination%x_e(a) - this%xe(i))*State%akthom/a2

                if (last_dotmu /=0) then
                    this%actual_opt_depth = this%actual_opt_depth - 2._dl*dtau/(1._dl/this%dotmu(i)+1._dl/last_dotmu)
                end if
                last_dotmu = this%dotmu(i)
            end if
        else
            this%xe(i)=State%Recombination%x_e(a)
        end if

        !  Baryon sound speed squared (over c**2).
        dtbdla=-2._dl*this%tb(i)-thomc*adothalf/adot*(a*this%tb(i)-CP%tcmb)
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        this%cs2(i)=barssc*this%tb(i)*(1-dtbdla/this%tb(i)/3._dl)


        ! Calculation of the visibility function
        this%dotmu(i)=this%xe(i)*State%akthom/a2

        if (this%tight_tau==0 .and. 1/(tau*this%dotmu(i)) > 0.005) this%tight_tau = tau !0.005
        !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)

        if (tau < 0.001) then
            this%sdotmu(i)=0
        else
            this%sdotmu(i)=this%sdotmu(i-1)+2._dl*dtau/(1._dl/this%dotmu(i)+1._dl/this%dotmu(i-1))
        end if

        a0=a
        tau01=tau
        adot0=adot
    end do !i

    if (CP%Reion%Reionization .and. (this%xe(nthermo) < 0.999d0)) then
        write(*,*)'Warning: xe at redshift zero is < 1'
        write(*,*) 'Check input parameters an Reionization_xe'
        write(*,*) 'function in the Reionization module'
    end if

    do j1=1,nthermo
        if (this%sdotmu(j1) - this%sdotmu(nthermo)< -69) then
            this%emmu(j1)=1.d-30
        else
            this%emmu(j1)=exp(this%sdotmu(j1)-this%sdotmu(nthermo))
            if (.not. CP%Accuracy%AccurateReionization .and. &
                this%actual_opt_depth==0 .and. this%xe(j1) < 1e-3) then
                this%actual_opt_depth = -this%sdotmu(j1)+this%sdotmu(nthermo)
            end if
            if (CP%ACcuracy%AccurateReionization .and. CP%DerivedParameters .and. this%z_star==0.d0) then
                if (this%sdotmu(nthermo)-this%sdotmu(j1) - this%actual_opt_depth < 1) then
                    tau01=1-(this%sdotmu(nthermo)-this%sdotmu(j1) - this%actual_opt_depth)
                    tau01=tau01*(1._dl/this%dotmu(j1)+1._dl/this%dotmu(j1-1))/2
                    this%z_star = 1/(scaleFactor(j1)- tau01/dtauda(scaleFactor(j1))) -1
                end if
            end if
        end if
    end do

    !Sources
    if (CP%WantTransfer) then
        if (associated(State%MT%optical_depths)) deallocate(State%MT%optical_depths, stat = RW_i)
        allocate(State%MT%optical_depths(CP%Transfer%num_redshifts))
        do RW_i = 1, CP%Transfer%num_redshifts
            if (CP%Transfer%Redshifts(RW_i) < 1e-3) then
                State%MT%optical_depths(RW_i) = 0 !zero may not be set correctly in transfer_ix
            else
                State%MT%optical_depths(RW_i) =  -this%sdotmu(transfer_ix(RW_i))+this%sdotmu(nthermo)
            end if
        end do
    end if

    if (CP%Accuracy%AccurateReionization .and. FeedbackLevel > 0 .and. CP%DerivedParameters) then
        write(*,'("Reion opt depth      = ",f7.4)') this%actual_opt_depth
    end if


    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time of maximum visibility
    if (ncount == 0) then
        cf1=1._dl
        ns=nthermo
    else
        cf1=exp(this%sdotmu(nthermo)-this%sdotmu(ncount))
        ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
        vis = this%emmu(j1)*this%dotmu(j1)
        tau = this%tauminn*exp((j1-1)*this%dlntau)
        vfi=vfi+vis*cf1*this%dlntau*tau
        if ((iv == 0).and.(vfi > 1.0d-7/CP%Accuracy%AccuracyBoost)) then
            State%taurst=9._dl/10._dl*tau
            iv=1
        elseif (iv == 1) then
            if (vis > maxvis) then
                maxvis=vis
                State%tau_maxvis = tau
            end if
            if (vfi > 0.995) then
                State%taurend=tau
                iv=2
                exit
            end if
        end if
    end do

    if (iv /= 2) then
        call GlobalError('ThemoData Init: failed to find end of recombination',error_reionization)
        return
    end if

    if (dowinlens) then
        vfi=0
        awin_lens1p=0
        awin_lens2p=0
        this%winlens=0
        do j1=1,nthermo-1
            vis = this%emmu(j1)*this%dotmu(j1)
            tau = this%tauminn*exp((j1-1)*this%dlntau)
            vfi=vfi+vis*cf1*this%dlntau*tau
            if (vfi < 0.995) then
                dwing_lens =  vis*cf1*this%dlntau*tau / 0.995

                awin_lens1p = awin_lens1p + dwing_lens
                awin_lens2p = awin_lens2p + dwing_lens/(State%tau0-tau)
            end if
            this%winlens(j1)= awin_lens1p/(State%tau0-tau) - awin_lens2p
        end do
    end if

    ! Calculating the timesteps during recombination.

    if (CP%WantTensors) then
        State%dtaurec=min(State%dtaurec,State%taurst/160)/CP%Accuracy%AccuracyBoost
    else
        State%dtaurec=min(State%dtaurec,State%taurst/40)/CP%Accuracy%AccuracyBoost
        if (do_bispectrum .and. hard_bispectrum) State%dtaurec = State%dtaurec / 4
    end if

    if (CP%Reion%Reionization) State%taurend=min(State%taurend,State%ReionHist%tau_start)

    if (DebugMsgs) then
        write (*,*) 'taurst, taurend = ', State%taurst, State%taurend
    end if

    !Sources
    do RW_i = 1, num_redshiftwindows
        associate(Win => RW(RW_i))
            if (Redshift_w(RW_i)%kind == window_lensing .or. &
                Redshift_w(RW_i)%kind == window_counts .and. CP%SourceTerms%do_counts_lensing) then
                this%has_lensing_windows = .true.
                Redshift_w(RW_i)%has_lensing_window = .true.
                if (FeedbackLevel>0) &
                    write(*,'(I1," Int W              = ",f9.6)') RW_i, awin_lens1(RW_i)

                Win%awin_lens=Win%awin_lens/awin_lens1(RW_i)
            else
                Redshift_w(RW_i)%has_lensing_window = .false.
            end if
        end associate
    end do

    call splini(spline_data,nthermo)
    call splder(this%cs2,this%dcs2,nthermo,spline_data)
    call splder(this%dotmu,this%ddotmu,nthermo,spline_data)
    call splder(this%ddotmu,this%dddotmu,nthermo,spline_data)
    call splder(this%dddotmu,this%ddddotmu,nthermo,spline_data)
    call splder(this%emmu,this%demmu,nthermo,spline_data)
    if (dowinlens) call splder(this%winlens,this%dwinlens,nthermo,spline_data)
    !Sources
    if (num_redshiftwindows >0) then
        call splder(this%redshift_time,this%dredshift_time,nthermo,spline_data)
        call splder(this%arhos_fac,this%darhos_fac,nthermo,spline_data)
        call splder(this%darhos_fac,this%ddarhos_fac,nthermo,spline_data)
        do j2 = 1, num_redshiftwindows
            if (Redshift_w(j2)%has_lensing_window) then
                call splder(RW(j2)%awin_lens,RW(j2)%dawin_lens,nthermo,spline_data)
            end if
        end do
    end if

    call this%SetTimeSteps(State%TimeSteps)

    if (num_redshiftwindows>0) then
        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
        do j2=1,State%TimeSteps%npoints
            call this%DoWindowSpline(j2,State%TimeSteps%points(j2), RW)
        end do
        !$OMP END PARALLEL DO
        call this%SetTimeStepWindows(State%TimeSteps)
    end if

    if ((CP%want_zstar .or. CP%DerivedParameters) .and. this%z_star==0.d0) call find_z(optdepth,this%z_star)
    if (CP%want_zdrag .or. CP%DerivedParameters) call find_z(dragoptdepth,this%z_drag)

    if (CP%DerivedParameters) then
        rs =rombint(dsound_da_exact,1d-8,1/(this%z_star+1),1d-6)
        DA = AngularDiameterDistance(State,this%z_star)/(1/(this%z_star+1))
        associate(ThermoDerivedParams => State%ThermoDerivedParams)
            ThermoDerivedParams( derived_Age ) = DeltaPhysicalTimeGyr(State,0.0_dl,1.0_dl)
            ThermoDerivedParams( derived_zstar ) = this%z_star
            ThermoDerivedParams( derived_rstar ) = rs
            ThermoDerivedParams( derived_thetastar ) = 100*rs/DA
            ThermoDerivedParams( derived_DAstar ) = DA/1000
            ThermoDerivedParams( derived_zdrag ) = this%z_drag
            rs =rombint(dsound_da_exact,1d-8,1/(this%z_drag+1),1d-6)
            ThermoDerivedParams( derived_rdrag ) = rs
            ThermoDerivedParams( derived_kD ) =  sqrt(1.d0/(rombint(ddamping_da, 1d-8, 1/(this%z_star+1), 1d-6)/6))
            ThermoDerivedParams( derived_thetaD ) =  100*const_pi/ThermoDerivedParams( derived_kD )/DA
            z_eq = (State%grhob+State%grhoc)/(State%grhog+State%grhornomass+sum(State%grhormass(1:CP%Nu_mass_eigenstates))) -1
            ThermoDerivedParams( derived_zEQ ) = z_eq
            a_eq = 1/(1+z_eq)
            ThermoDerivedParams( derived_kEQ ) = 1/(a_eq*dtauda(a_eq))
            ThermoDerivedParams( derived_thetaEQ ) = 100*timeOfz(State,z_eq)/DA
            ThermoDerivedParams( derived_theta_rs_EQ ) = 100*rombint(dsound_da_exact,1d-8,a_eq,1d-6)/DA

            associate(BackgroundOutputs => State%BackgroundOutputs)
                if (associated(BackgroundOutputs%z_outputs)) then
                    if (allocated(BackgroundOutputs%H)) &
                        deallocate(BackgroundOutputs%H, BackgroundOutputs%DA, BackgroundOutputs%rs_by_D_v)
                    noutput = size(BackgroundOutputs%z_outputs)
                    allocate(BackgroundOutputs%H(noutput), BackgroundOutputs%DA(noutput), BackgroundOutputs%rs_by_D_v(noutput))
                    do i=1,noutput
                        BackgroundOutputs%H(i) = HofZ(State,BackgroundOutputs%z_outputs(i))
                        BackgroundOutputs%DA(i) = AngularDiameterDistance(State,BackgroundOutputs%z_outputs(i))
                        BackgroundOutputs%rs_by_D_v(i) = rs/BAO_D_v_from_DA_H(BackgroundOutputs%z_outputs(i), &
                            BackgroundOutputs%DA(i),BackgroundOutputs%H(i))
                    end do
                end if
            end associate

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
        end associate
    end if

    !Sources
    if(num_redshiftwindows>0) then
        deallocate(RW,this%arhos_fac, this%darhos_fac, this%ddarhos_fac)
        deallocate(this%redshift_time, this%dredshift_time)
        do RW_i = 1, num_redshiftwindows
            associate (RedWin => Redshift_W(RW_i))
                if (RedWin%kind == window_21cm) then
                    outstr = 'z= '//trim(RealToStr(real(RedWin%Redshift),4))//': T_b = '//trim(RealToStr(real(RedWin%Fq),6))// &
                        'mK; tau21 = '//trim(RealToStr(real(RedWin%optical_depth_21),5))
                    write (*,*) RW_i,trim(outstr)
                end if
            end associate
        end do
    end if

    if (CP%Do21cm .and. CP%transfer_21cm_cl) then
        !           do RW_i=15,800
        !
        !             a=1._dl/(1+RW_i)
        !             Tspin = tsRecfast(a)
        !             Trad = CP%TCMB/a
        !             adot = 1/dtauda(a)
        !             tau_eps = a**2*line21_const*NNow/a**3/adot/Tspin/1000
        !             Tb21cm = 1000*(1-exp(-tau_eps))*a*(Tspin-Trad)
        !
        !             write (*,'(3E15.5)') 1/a-1, Tb21cm, tau_eps
        !            end do
        !           stop

        do RW_i=1,CP%Transfer%PK_num_redshifts
            a=1._dl/(1+CP%Transfer%PK_redshifts(RW_i))
            Tspin = State%Recombination%T_s(a)
            Trad = CP%TCMB/a
            adot = 1/dtauda(a)
            tau_eps = a**2*line21_const*State%NNow/a**3/adot/Tspin/1000
            Tb21cm = 1000*(1-exp(-tau_eps))*a*(Tspin-Trad)

            outstr = 'z= '//trim(RealToStr(real(CP%Transfer%PK_redshifts(RW_i))))// &
                ': tau_21cm = '//trim(RealToStr(real(tau_eps),5))//'; T_b = '//trim(RealToStr(real(Tb21cm),6))//'mK'
            write (*,*) trim(outstr)
        end do
    end if

    end subroutine Thermo_Init


    subroutine SetTimeSteps(this,TimeSteps)
    !Set time steps to use for sampling the source functions for the CMB power spectra
    class(TThermoData) :: this
    Type(TRanges) :: TimeSteps
    real(dl) dtau0
    integer nri0, nstep
    !Sources
    integer ix,i,nwindow, L_limb
    real(dl) keff, win_end

    call TimeSteps%Init()

    call TimeSteps%Add_delta(State%taurst, State%taurend, State%dtaurec)

    ! Calculating the timesteps after recombination
    if (CP%WantTensors) then
        dtau0=max(State%taurst/40,this%Maxtau/2000._dl/CP%Accuracy%AccuracyBoost)
    else
        dtau0=this%Maxtau/500._dl/CP%Accuracy%AccuracyBoost
        if (do_bispectrum) dtau0 = dtau0/3
        !Don't need this since adding in Limber on small scales
        !  if (CP%DoLensing) dtau0=dtau0/2
        !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
    end if

    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)

    !Sources

    this%tau_start_redshiftwindows = this%MaxTau
    this%tau_end_redshiftwindows = 0
    if (CP%WantScalars .or. CP%SourceTerms%line_reionization) then
        do ix=1, num_redshiftwindows
            associate (Win => Redshift_W(ix))

                Win%sigma_tau = Win%sigma_z*dtauda(1/(1+Win%Redshift))/(1+Win%Redshift)**2

                if (Win%tau_start==0) then
                    Win%tau_start = max(Win%tau - Win%sigma_tau*7,state%taurst)
                    Win%tau_end = min(State%tau0,Win%tau + Win%sigma_tau*7)
                end if

                this%tau_start_redshiftwindows = min(Win%tau_start,this%tau_start_redshiftwindows)

                if (Win%kind /= window_lensing) then
                    !Have to be careful to integrate dwinV as the window tails off
                    this%tau_end_redshiftwindows = max(Win%tau_end,this%tau_end_redshiftwindows)
                    nwindow = nint(150*CP%Accuracy%AccuracyBoost)
                    win_end = Win%tau_end
                else !lensing
                    nwindow = nint(CP%Accuracy%AccuracyBoost*Win%chi0/100)
                    win_end = State%tau0
                end if

                if (Win%kind == window_21cm .and. (CP%SourceTerms%line_phot_dipole .or. &
                    CP%SourceTerms%line_phot_quadrupole)) nwindow = nwindow *3

                L_limb = Win_limber_ell(Win,CP%max_l)
                keff = WindowKmaxForL(Win,L_limb)

                !Keep sampling in x better than Nyquist
                nwindow = max(nwindow, nint(CP%Accuracy%AccuracyBoost *(win_end- Win%tau_start)* keff/3))
                if (Feedbacklevel > 1) write (*,*) ix, 'nwindow =', nwindow

                call TimeSteps%Add(Win%tau_start, win_end, nwindow)
                !This should cover whole range where not tiny

                if (Win%kind /= window_lensing .and. Win%tau_end - Win%tau_start > Win%sigma_tau*7) then
                    call TimeSteps%Add(TimeOfZ(State,Win%Redshift+Win%sigma_z*3), &
                        max(0._dl,TimeOfZ(State,Win%Redshift-Win%sigma_z*3)), nwindow)
                    !This should be over peak
                end if
                !Make sure line of sight integral OK too
                ! if (dtau0 > Win%tau_end/300/AccuracyBoost) then
                !  call Ranges_Add_delta(TimeSteps, Win%tau_end, CP%tau0,  Win%tau_start/300/AccuracyBoost)
                ! end if
            end associate
        end do
    end if


    if (CP%Reion%Reionization) then
        nri0=int(Reionization_timesteps(State%ReionHist)*CP%Accuracy%AccuracyBoost)
        !Steps while reionization going from zero to maximum
        call TimeSteps%Add(State%ReionHist%tau_start,State%ReionHist%tau_complete,nri0)
    end if

    !Sources
    if (.not. CP%Want_CMB .and. CP%WantCls) then
        if (num_redshiftwindows==0) call MpiStop('Want_CMB=false, but not redshift windows either')
        call TimeSteps%Add_delta(this%tau_start_redshiftwindows, State%tau0, dtau0)
    end if

    !Create arrays out of the region information.
    call TimeSteps%GetArray()
    nstep = TimeSteps%npoints

    if (num_redshiftwindows > 0) then
        if (allocated(this%step_redshift)) deallocate(this%step_redshift, this%rhos_fac, this%drhos_fac)
        allocate(this%step_redshift(nstep), this%rhos_fac(nstep), this%drhos_fac(nstep))
        do i=1,num_redshiftwindows
            associate (Win => Redshift_W(i))
                allocate(Win%winF(nstep),Win%wing(nstep),Win%dwing(nstep),Win%ddwing(nstep), &
                    Win%winV(nstep),Win%dwinV(nstep),Win%ddwinV(nstep))
                allocate(Win%win_lens(nstep),Win%wing2(nstep),Win%dwing2(nstep),Win%ddwing2(nstep))
                allocate(Win%wingtau(nstep),Win%dwingtau(nstep),Win%ddwingtau(nstep))
                if (Win%kind == window_counts) then
                    allocate(Win%comoving_density_ev(nstep))
                end if
            end associate
        end do
    end if

    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Set ',nstep, ' time steps'

    end subroutine SetTimeSteps

    subroutine SetTimeStepWindows(this,TimeSteps)
    use Recombination
    use constants
    class(TThermoData) :: this
    Type(TRanges) :: TimeSteps
    integer i, j, jstart, ix
    real(dl) tau,  a, a2
    real(dl) Tspin, Trad, rho_fac, tau_eps
    real(dl) window, winamp
    real(dl) z,rhos, adot, exp_fac
    real(dl) tmp(TimeSteps%npoints), tmp2(TimeSteps%npoints), hubble_tmp(TimeSteps%npoints)
    real(dl), allocatable , dimension(:,:) :: int_tmp, back_count_tmp
    integer ninterp

    ! Prevent false positive warnings for uninitialized
    Tspin = 0._dl
    Trad = 0._dl
    tau_eps = 0._dl
    exp_fac = 0._dl

    jstart = TimeSteps%IndexOf(this%tau_start_redshiftwindows)
    ninterp = TimeSteps%npoints - jstart + 1

    do i = 1, num_redshiftwindows
        associate (RedWin => Redshift_W(i))
            RedWin%wing=0
            RedWin%winV=0
            RedWin%winF=0
            RedWin%wing2=0
            RedWin%dwing=0
            RedWin%dwinV=0
            RedWin%dwing2=0
            RedWin%ddwing=0
            RedWin%ddwinV=0
            RedWin%ddwing2=0
            RedWin%wingtau=0
            RedWin%dwingtau=0
            RedWin%ddwingtau=0
            RedWin%Fq = 0
            if (RedWin%kind == window_counts) then
                RedWin%comoving_density_ev  = 0
            end if
        end associate
    end do

    ! print *,'z = ',TimeSteps%points(jstart),step_redshift(jstart)
    allocate(int_tmp(jstart:TimeSteps%npoints,num_redshiftwindows))
    int_tmp = 0
    allocate(back_count_tmp(jstart:TimeSteps%npoints,num_redshiftwindows))
    back_count_tmp = 0

    do j=jstart, TimeSteps%npoints
        tau = TimeSteps%points(j)
        z = this%step_redshift(j)
        a = 1._dl/(1._dl+z)
        a2=a**2
        adot=1._dl/dtauda(a)


        if (CP%Do21cm) then
            Tspin = State%Recombination%T_s(a)
            Trad = CP%TCMB/a
            rho_fac = line21_const*State%NNow/a**3 !neglect ionization fraction
            tau_eps = a*rho_fac/(adot/a)/Tspin/1000
            exp_fac =  (1-exp(-tau_eps))/tau_eps
        else
            rho_fac = State%grhoc/a**3
        end if

        hubble_tmp(j) = adot/a

        do i = 1, num_redshiftwindows
            associate (RedWin => Redshift_W(i))
                if (tau < RedWin%tau_start) cycle

                window = Window_f_a(RedWin, a, winamp)

                if (RedWin%kind == window_21cm) then
                    rhos = rho_fac*(1 - Trad/Tspin)

                    !Want to integrate this...
                    int_tmp(j,i) = this%drhos_fac(j)*a*window

                    RedWin%WinV(j) = -exp(-tau_eps)*a2*rhos*window/(adot/a)

                    RedWin%wing(j) = exp_fac*a2*rhos*window

                    !The window that doesn't go to zero at T_s = T_gamma
                    RedWin%wing2(j) = exp_fac*a2*rho_fac*Trad/Tspin*window

                    !Window for tau_s for the self-absoption term
                    RedWin%wingtau(j) =  RedWin%wing(j)*(1 - exp(-tau_eps)/exp_fac)
                elseif ( RedWin%kind == window_counts) then

                    !window is n(a) where n is TOTAL not fractional number
                    !delta = int wing(eta) delta(eta) deta
                    RedWin%wing(j) = adot *window

                    !Window with 1/H in
                    RedWin%wing2(j) = RedWin%wing(j)/(adot/a)

                    !winv is g/chi for the ISW and time delay terms
                    RedWin%WinV(j) = 0
                    if (tau < State%tau0 -0.1) then
                        int_tmp(j,i) = RedWin%wing(j)/(State%tau0 - tau)
                    else
                        int_tmp(j,i)=0
                    end if

                    if (CP%SourceTerms%counts_evolve) then
                        back_count_tmp(j,i) =  counts_background_z(RedWin, 1/a-1)/a
                        if (tau < State%tau0 -0.1) then
                            RedWin%comoving_density_ev(j) = back_count_tmp(j,i)*(adot/a)/(State%tau0 - tau)**2
                        else
                            RedWin%comoving_density_ev(j) = 0
                        end if
                    end if
                end if
            end associate
            !Lensing windows should be fine from inithermo
        end do
    end do


    do i = 1, num_redshiftwindows
        associate (RedWin => Redshift_W(i))

            ! int (a*rho_s/H)' a W_f(a) d\eta, or for counts int g/chi deta
            call spline(TimeSteps%points(jstart),int_tmp(jstart,i),ninterp,spl_large,spl_large,tmp)
            call spline_integrate(TimeSteps%points(jstart),int_tmp(jstart,i),tmp, tmp2(jstart),ninterp)
            RedWin%WinV(jstart:TimeSteps%npoints) =  &
                RedWin%WinV(jstart:TimeSteps%npoints) + tmp2(jstart:TimeSteps%npoints)

            call spline(TimeSteps%points(jstart),RedWin%WinV(jstart),ninterp,spl_large,spl_large,RedWin%ddWinV(jstart))
            call spline_deriv(TimeSteps%points(jstart),RedWin%WinV(jstart),RedWin%ddWinV(jstart), RedWin%dWinV(jstart), ninterp)

            call spline(TimeSteps%points(jstart),RedWin%Wing(jstart),ninterp,spl_large,spl_large,RedWin%ddWing(jstart))
            call spline_deriv(TimeSteps%points(jstart),RedWin%Wing(jstart),RedWin%ddWing(jstart), RedWin%dWing(jstart), ninterp)

            call spline(TimeSteps%points(jstart),RedWin%Wing2(jstart),ninterp,spl_large,spl_large,RedWin%ddWing2(jstart))
            call spline_deriv(TimeSteps%points(jstart),RedWin%Wing2(jstart),RedWin%ddWing2(jstart), &
                RedWin%dWing2(jstart), ninterp)

            call spline_integrate(TimeSteps%points(jstart),RedWin%Wing(jstart),RedWin%ddWing(jstart), RedWin%WinF(jstart),ninterp)
            RedWin%Fq = RedWin%WinF(TimeSteps%npoints)

            if (RedWin%kind == window_21cm) then
                call spline_integrate(TimeSteps%points(jstart),RedWin%Wing2(jstart),&
                    RedWin%ddWing2(jstart), tmp(jstart),ninterp)
                RedWin%optical_depth_21 = tmp(TimeSteps%npoints) / (CP%TCMB*1000)
                !WinF not used.. replaced below

                call spline(TimeSteps%points(jstart),RedWin%Wingtau(jstart),ninterp,spl_large,spl_large,RedWin%ddWingtau(jstart))
                call spline_deriv(TimeSteps%points(jstart),RedWin%Wingtau(jstart),RedWin%ddWingtau(jstart), &
                    RedWin%dWingtau(jstart), ninterp)
            elseif (RedWin%kind == window_counts) then

                if (CP%SourceTerms%counts_evolve) then
                    call spline(TimeSteps%points(jstart),back_count_tmp(jstart,i),ninterp,spl_large,spl_large,tmp)
                    call spline_deriv(TimeSteps%points(jstart),back_count_tmp(jstart,i),tmp,tmp2(jstart),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%Wingtau(ix) = 0
                        else
                            !evo bias is computed with total derivative
                            RedWin%Wingtau(ix) =  -tmp2(ix) * RedWin%Wing(ix) / (back_count_tmp(ix,i)*hubble_tmp(ix)) &
                                !+ 5*RedWin%dlog10Ndm * ( RedWin%Wing(ix)- int_tmp(ix,i)/hubble_tmp(ix))
                                !The correction from total to partial derivative takes 1/adot(tau0-tau) cancels
                                + 10*RedWin%dlog10Ndm * RedWin%Wing(ix)
                        end if
                    end do

                    !comoving_density_ev is d log(a^3 n_s)/d eta * window
                    call spline(TimeSteps%points(jstart),RedWin%comoving_density_ev(jstart),ninterp,spl_large,spl_large,tmp)
                    call spline_deriv(TimeSteps%points(jstart),RedWin%comoving_density_ev(jstart),tmp,tmp2(jstart),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%comoving_density_ev(ix) = 0
                        elseif (RedWin%comoving_density_ev(ix)/=0._dl) then
                            !correction needs to be introduced from total derivative to partial derivative
                            RedWin%comoving_density_ev(ix) =   tmp2(ix) / RedWin%comoving_density_ev(ix) &
                                -5*RedWin%dlog10Ndm * ( hubble_tmp(ix) + int_tmp(ix,i)/RedWin%Wing(ix))
                        end if
                    end do
                else
                    RedWin%comoving_density_ev=0
                    call spline(TimeSteps%points(jstart),hubble_tmp(jstart),ninterp,spl_large,spl_large,tmp)
                    call spline_deriv(TimeSteps%points(jstart),hubble_tmp(jstart),tmp, tmp2(jstart), ninterp)

                    !assume d( a^3 n_s) of background population is zero, so remaining terms are
                    !wingtau =  g*(2/H\chi + Hdot/H^2)  when s=0; int_tmp = window/chi
                    RedWin%Wingtau(jstart:TimeSteps%npoints) = &
                        2*(1-2.5*RedWin%dlog10Ndm)*int_tmp(jstart:TimeSteps%npoints,i)/&
                        hubble_tmp(jstart:TimeSteps%npoints)&
                        + 5*RedWin%dlog10Ndm*RedWin%Wing(jstart:TimeSteps%npoints) &
                        + tmp2(jstart:TimeSteps%npoints)/hubble_tmp(jstart:TimeSteps%npoints)**2 &
                        *RedWin%Wing(jstart:TimeSteps%npoints)
                endif

                call spline(TimeSteps%points(jstart),RedWin%Wingtau(jstart),ninterp, &
                    spl_large,spl_large,RedWin%ddWingtau(jstart))
                call spline_deriv(TimeSteps%points(jstart),RedWin%Wingtau(jstart),RedWin%ddWingtau(jstart), &
                    RedWin%dWingtau(jstart), ninterp)

                !WinF is int[ g*(...)]
                call spline_integrate(TimeSteps%points(jstart),RedWin%Wingtau(jstart),&
                    RedWin%ddWingtau(jstart), RedWin%WinF(jstart),ninterp)
            end if
        end associate
    end do

    deallocate(int_tmp,back_count_tmp)

    end subroutine SetTimeStepWindows


    subroutine interp_window(TimeSteps,RedWin,tau,wing_t, wing2_t, winv_t)
    !for evolving sources for reionization we neglect wingtau self-absorption
    Type(TRanges) :: TimeSteps
    Type(TRedWin)  :: RedWin
    integer i
    real(dl) :: tau, wing_t, wing2_t,winv_t
    real(dl) a0,b0,ho

    i = TimeSteps%IndexOf(tau)
    if (i< TimeSteps%npoints) then
        ho=TimeSteps%points(i+1)-TimeSteps%points(i)
        a0=(TimeSteps%points(i+1)-tau)/ho
        b0=1-a0
        wing_t = a0*RedWin%wing(i)+ b0*RedWin%wing(i+1)+((a0**3-a0)* RedWin%ddwing(i) &
            +(b0**3-b0)*RedWin%ddwing(i+1))*ho**2/6
        wing2_t = a0*RedWin%wing2(i)+ b0*RedWin%wing2(i+1)+((a0**3-a0)* RedWin%ddwing2(i) &
            +(b0**3-b0)*RedWin%ddwing2(i+1))*ho**2/6
        winv_t = a0*RedWin%winv(i)+ b0*RedWin%winv(i+1)+((a0**3-a0)* RedWin%ddwinv(i) &
            +(b0**3-b0)*RedWin%ddwinv(i+1))*ho**2/6
    else
        wing_t = 0
        wing2_t = 0
        winv_t = 0
    end if
    end subroutine interp_window

    subroutine DoWindowSpline(this,j2,tau, RW)
    class(TThermoData) :: this
    integer j2, i, RW_i
    real(dl) d, tau
    Type(CalWins) :: RW(:)

    !     Cubic-spline interpolation.
    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)

    d=d-i
    if (i < this%nthermo) then

        this%step_redshift(j2) = this%redshift_time(i)+d*(this%dredshift_time(i)+ &
            d*(3._dl*(this%redshift_time(i+1)-this%redshift_time(i)) &
            -2._dl*this%dredshift_time(i)-this%dredshift_time(i+1)+d*(this%dredshift_time(i)+this%dredshift_time(i+1) &
            +2._dl*(this%redshift_time(i)-this%redshift_time(i+1)))))

        this%rhos_fac(j2) = this%arhos_fac(i)+d*(this%darhos_fac(i)+d*(3._dl*(this%arhos_fac(i+1)-this%arhos_fac(i)) &
            -2._dl*this%darhos_fac(i)-this%darhos_fac(i+1)+d*(this%darhos_fac(i)+this%darhos_fac(i+1) &
            +2._dl*(this%arhos_fac(i)-this%arhos_fac(i+1)))))
        this%drhos_fac(j2) = (this%darhos_fac(i)+d*(this%ddarhos_fac(i)+d*(3._dl*(this%darhos_fac(i+1)  &
            -this%darhos_fac(i))-2._dl*this%ddarhos_fac(i)-this%ddarhos_fac(i+1)+d*(this%ddarhos_fac(i) &
            +this%ddarhos_fac(i+1)+2._dl*(this%darhos_fac(i)-this%darhos_fac(i+1))))))/(tau &
            *this%dlntau)

        do RW_i=1, num_redshiftwindows
            if (Redshift_w(RW_i)%has_lensing_window) then
                associate(W => Redshift_W(RW_i), C=> RW(RW_i))

                    W%win_lens(j2) = C%awin_lens(i)+d*(C%dawin_lens(i)+d*(3._dl*(C%awin_lens(i+1)-C%awin_lens(i)) &
                        -2._dl*C%dawin_lens(i)-C%dawin_lens(i+1)+d*(C%dawin_lens(i)+C%dawin_lens(i+1) &
                        +2._dl*(C%awin_lens(i)-C%awin_lens(i+1)))))
                end associate
            end if
        end do

    else
        this%step_redshift(j2) = 0
        this%rhos_fac(j2)=0
        this%drhos_fac(j2)=0
        do RW_i=1, num_redshiftwindows
            associate (W => Redshift_W(RW_i))
                W%win_lens(j2)=0
            end associate
        end do
    end if

    end subroutine DoWindowSpline

    subroutine IonizationFunctionsAtTime(this,tau, opac, dopac, ddopac, &
        vis, dvis, ddvis, expmmu, lenswin)
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out):: opac, dopac, ddopac, vis, dvis, ddvis, expmmu, lenswin
    real(dl) d
    integer i

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i

    if (i < this%nthermo) then
        opac=this%dotmu(i)+d*(this%ddotmu(i)+d*(3._dl*(this%dotmu(i+1)-this%dotmu(i)) &
            -2._dl*this%ddotmu(i)-this%ddotmu(i+1)+d*(this%ddotmu(i)+this%ddotmu(i+1) &
            +2._dl*(this%dotmu(i)-this%dotmu(i+1)))))
        dopac=(this%ddotmu(i)+d*(this%dddotmu(i)+d*(3._dl*(this%ddotmu(i+1)  &
            -this%ddotmu(i))-2._dl*this%dddotmu(i)-this%dddotmu(i+1)+d*(this%dddotmu(i) &
            +this%dddotmu(i+1)+2._dl*(this%ddotmu(i)-this%ddotmu(i+1))))))/(tau &
            *this%dlntau)
        ddopac=(this%dddotmu(i)+d*(this%ddddotmu(i)+d*(3._dl*(this%dddotmu(i+1) &
            -this%dddotmu(i))-2._dl*this%ddddotmu(i)-this%ddddotmu(i+1)  &
            +d*(this%ddddotmu(i)+this%ddddotmu(i+1)+2._dl*(this%dddotmu(i) &
            -this%dddotmu(i+1)))))-(this%dlntau**2)*tau*dopac) &
            /(tau*this%dlntau)**2
        expmmu=this%emmu(i)+d*(this%demmu(i)+d*(3._dl*(this%emmu(i+1)-this%emmu(i)) &
            -2._dl*this%demmu(i)-this%demmu(i+1)+d*(this%demmu(i)+this%demmu(i+1) &
            +2._dl*(this%emmu(i)-this%emmu(i+1)))))

        if (dowinlens) then
            lenswin=this%winlens(i)+d*(this%dwinlens(i)+d*(3._dl*(this%winlens(i+1)-this%winlens(i)) &
                -2._dl*this%dwinlens(i)-this%dwinlens(i+1)+d*(this%dwinlens(i)+this%dwinlens(i+1) &
                +2._dl*(this%winlens(i)-this%winlens(i+1)))))
        end if
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3*opac*dopac+ddopac)
    else
        opac=this%dotmu(this%nthermo)
        dopac=this%ddotmu(this%nthermo)
        ddopac=this%dddotmu(this%nthermo)
        expmmu=this%emmu(this%nthermo)
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3._dl*opac*dopac+ddopac)
    end if

    end subroutine IonizationFunctionsAtTime

    function ddamping_da(a)
    real(dl) :: ddamping_da
    real(dl), intent(in) :: a
    real(dl) :: R

    R=State%ThermoData%r_drag0*a
    !ignoring reionisation, not relevant for distance measures
    ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(a)*a**2/(State%Recombination%x_e(a)*State%akthom)

    end function ddamping_da


    !!!!!!!!!!!!!!!!!!!
    !JH: functions and subroutines for calculating z_star and z_drag

    function doptdepth_dz(z)
    real(dl) :: doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)

    !ignoring reionisation, not relevant for distance measures
    doptdepth_dz = State%Recombination%x_e(a)*State%akthom*dtauda(a)

    end function doptdepth_dz


    function optdepth(z)
    real(dl) optdepth
    real(dl),intent(in) :: z

    optdepth = rombint2(doptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function optdepth


    function ddragoptdepth_dz(z)
    real(dl) :: ddragoptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)
    ddragoptdepth_dz = doptdepth_dz(z)/State%ThermoData%r_drag0/a

    end function ddragoptdepth_dz


    function dragoptdepth(z)
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

    subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call CTrans%q%getArray(.true.)

    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,&
        min(CTrans%max_index_nonlimber,CTrans%ls%nl), CTrans%q%npoints),  STAT = st)
    if (st /= 0) call MpiStop('Init_ClTransfer: Error allocating memory for transfer functions')
    CTrans%Delta_p_l_k = 0

    end subroutine Init_ClTransfer

    subroutine Init_Limber(CTrans)
    Type(ClTransferData) :: CTrans

    call Free_Limber(Ctrans)
    allocate(CTrans%Limber_l_min(CTrans%NumSources))
    CTrans%Limber_l_min = 0
    if (num_redshiftwindows>0 .or. limber_phiphi>0) then
        allocate(CTrans%Limber_windows(CTrans%NumSources,CTrans%ls%nl))
    end if

    end subroutine Init_Limber

    subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    nullify(CTrans%Delta_p_l_k)
    call CTrans%q%Free()
    call Free_Limber(CTrans)

    end subroutine Free_ClTransfer

    subroutine Free_Limber(CTrans)
    Type(ClTransferData) :: CTrans
    integer st,i,j

    if (associated(CTrans%Limber_l_min)) then
        do i=1, CTrans%NumSources
            if (CTrans%Limber_l_min(i)/=0) then
                do j=CTrans%Limber_l_min(i), CTrans%ls%nl
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


    function Win_Limber_ell(W,lmax) result(ell_limb)
    Type(TRedWin) :: W
    integer, intent(in) :: lmax
    integer ell_limb

    if (CP%SourceTerms%limber_windows) then
        !Turn on limber when k is a scale smaller than window width
        if (W%kind==window_lensing) then
            ell_limb = max(limber_phiphi,nint(50*CP%Accuracy%AccuracyBoost))
        else
            ell_limb = max(limber_phiphi, nint(CP%Accuracy%AccuracyBoost *6* W%chi0/W%sigma_tau))
        end if
    else
        ell_limb = lmax
    end if
    end function Win_Limber_ell


    subroutine CheckLoadedHighLTemplate
    use FileUtils
    integer :: L
    real(dl) :: array(7)
    type(TTextFile) :: infile

    if (.not. allocated(highL_CL_template)) then
        allocate(highL_CL_template(lmin:lmax_extrap_highl, C_Temp:C_Phi))
        call infile%Open(highL_unlensed_cl_template)
        if (lmin==1) highL_CL_template(lmin,:)=0
        do
            read(infile%unit, *, end=500) L, array
            if (L>lmax_extrap_highl) exit
            !  array = array * (2*l+1)/(4*pi) * 2*pi/(l*(l+1))
            highL_CL_template(L, C_Temp:C_E) =array(1:2)
            highL_CL_template(L, C_Cross) =array(4)
            highL_CL_template(L, C_Phi) =array(5)
        end do
500     if (L< lmax_extrap_highl) &
            call MpiStop('CheckLoadedHighLTemplate: template file does not go up to lmax_extrap_highl')
        call infile%Close()
    end if

    end subroutine CheckLoadedHighLTemplate


    subroutine TCLdata_InitCls(this, CP)
    class(TCLData) :: this
    class(CAMBParams) :: CP

    call CheckLoadedHighLTemplate
    if (CP%WantScalars) then
        if (allocated(this%Cl_scalar)) deallocate(this%Cl_scalar)
        allocate(this%Cl_scalar(lmin:CP%Max_l, C_Temp:C_last), source=0._dl)
        if (CP%want_cl_2D_array) then
            if (allocated(this%Cl_scalar_array)) deallocate(this%Cl_scalar_array)
            allocate(this%Cl_scalar_Array(lmin:CP%Max_l,  &
                3+num_redshiftwindows+num_custom_sources,3+num_redshiftwindows+num_custom_sources))
            this%Cl_scalar_array = 0
        end if
    end if

    if (CP%WantVectors) then
        if (allocated(this%Cl_vector)) deallocate(this%Cl_vector)
        allocate(this%Cl_vector(lmin:CP%Max_l, CT_Temp:CT_Cross), source=0._dl)
    end if

    if (CP%WantTensors) then
        if (allocated(this%Cl_tensor)) deallocate(this%Cl_tensor)
        allocate(this%Cl_tensor(lmin:CP%Max_l_tensor, CT_Temp:CT_Cross), source=0._dl)
    end if
    end subroutine TCLdata_InitCls

    function open_file_header(filename, Col1, Columns, n) result(unit)
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: col1
    character(LEN=name_tag_len), intent(in) :: Columns(:)
    integer, intent(in), optional :: n
    integer :: unit, nn

    nn = PresentDefault(6, n)

    open(newunit=unit,file=filename,form='formatted',status='replace')
    if (output_file_headers) then
        write(unit,'("#",1A'//Trim(IntToStr(nn-1))//'," ",*(A15))') Col1, Columns
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

    subroutine TCLdata_output_cl_files(this,ScalFile,ScalCovFile,TensFile, TotFile, LensFile, LensTotFile, factor)
    class(TCLData) :: this
    character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile,ScalCovfile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact
    integer :: last_C, il, i, j, unit
    real(dl), allocatable :: outarr(:,:)
    character(LEN=name_tag_len) :: cov_names((3+num_redshiftwindows)**2)

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. ScalFile /= '') then
        last_C=min(C_PhiTemp,C_last)
        unit = open_file_header(ScalFile, 'L', C_name_tags(:last_C))
        do il=lmin,min(10000,CP%Max_l)
            write(unit,trim(numcat('(1I6,',last_C))//'E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        do il=10100,CP%Max_l, 100
            write(unit,trim(numcat('(1E15.6,',last_C))//'E15.6)') real(il),&
                fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%want_cl_2D_array .and. ScalCovFile /= '' .and. this%CTransScal%NumSources>2) then
        allocate(outarr(1:3+num_redshiftwindows,1:3+num_redshiftwindows))

        do i=1, 3+num_redshiftwindows
            do j=1, 3+num_redshiftwindows
                cov_names(j + (i-1)*(3+num_redshiftwindows)) = trim(scalar_fieldname(i))//'x'//trim(scalar_fieldname(j))
            end do
        end do
        unit = open_file_header(ScalCovFile, 'L', cov_names)

        do il=lmin,min(10000,CP%Max_l)
            outarr=this%Cl_scalar_array(il,1:3+num_redshiftwindows,1:3+num_redshiftwindows)
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            write(unit,trim(numcat('(1I6,',(3+num_redshiftwindows)**2))//'E15.6)') il, real(outarr)
        end do
        do il=10100,CP%Max_l, 100
            outarr=this%Cl_scalar_array(il,1:3+num_redshiftwindows,1:3+num_redshiftwindows)
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            write(unit,trim(numcat('(1E15.5,',(3+num_redshiftwindows)**2))//'E15.6)') real(il), real(outarr)
        end do
        close(unit)
        deallocate(outarr)
    end if

    if (CP%WantTensors .and. TensFile /= '') then
        unit = open_file_header(TensFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_tensor(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if

    if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
        unit = open_file_header(TotFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_scalar(il, C_Temp:C_E)+ this%Cl_tensor(il,C_Temp:C_E)), &
                fact*this%Cl_tensor(il, CT_B), fact*(this%Cl_scalar(il, C_Cross) + this%Cl_tensor(il, CT_Cross))
        end do
        do il=CP%Max_l_tensor+1,CP%Max_l
            write(unit,'(1I6,4E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:C_E), 0._dl, fact*this%Cl_scalar(il,C_Cross)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
        unit = open_file_header(LensFile, 'L', CT_name_tags)
        do il=lmin, this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if


    if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
        unit = open_file_header(LensTotFile, 'L', CT_name_tags)
        do il=lmin,min(CP%Max_l_tensor,this%lmax_lensed)
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_lensed(il,CT_Temp:CT_Cross)+ &
                this%Cl_tensor(il, CT_Temp:CT_Cross))
        end do
        do il=min(CP%Max_l_tensor,this%lmax_lensed)+1,this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_cl_files

    subroutine TCLdata_output_lens_pot_files(this,LensPotFile, factor)
    !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed
    !This input supported by LensPix from 2010
    class(TCLdata) :: this
    integer :: il, unit
    real(dl), intent(in), optional :: factor
    real(dl) fact, scale, BB, TT, TE, EE
    character(LEN=*), intent(in) :: LensPotFile
    !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi and [l(l+1)]^(3/2) C_phi_T/2pi
    !This is the format used by Planck_like but original LensPix uses scalar_output_file.

    !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi
    ! - for historical reasons)

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then
        unit = open_file_header(LensPotFile, 'L', lens_pot_name_tags)
        do il=lmin, min(10000,CP%Max_l)
            TT = this%Cl_scalar(il, C_Temp)
            EE = this%Cl_scalar(il, C_E)
            TE = this%Cl_scalar(il, C_Cross)
            if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                TT= TT+this%Cl_tensor(il, CT_Temp)
                EE= EE+this%Cl_tensor(il, CT_E)
                TE= TE+this%Cl_tensor(il, CT_Cross)
                BB= this%Cl_tensor(il, CT_B)
            else
                BB=0
            end if
            scale = (real(il+1)/il)**2/OutputDenominator !Factor to go from old l^4 factor to new

            write(unit,'(1I6,7E15.6)') il , fact*TT, fact*EE, fact*BB, fact*TE, scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        do il=10100,CP%Max_l, 100
            scale = (real(il+1)/il)**2/OutputDenominator
            write(unit,'(1E15.5,7E15.6)') real(il), fact*this%Cl_scalar(il,C_Temp:C_E),0.,&
                fact*this%Cl_scalar(il,C_Cross), scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_lens_pot_files


    subroutine TCLdata_output_veccl_files(this,VecFile, factor)
    class(TCLData) :: this
    integer :: il, unit
    character(LEN=*), intent(in) :: VecFile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact

    fact = PresentDefault(1._dl, factor)

    if (CP%WantVectors .and. VecFile /= '') then
        unit = open_file_header(VecFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_vector(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if

    end subroutine TCLdata_output_veccl_files

    subroutine TCLdata_NormalizeClsAtL(this,lnorm)
    class(TCLData) :: this
    integer, intent(IN) :: lnorm
    real(dl) Norm

    if (CP%WantScalars) then
        Norm=1/this%Cl_scalar(lnorm, C_Temp)
        this%Cl_scalar(lmin:CP%Max_l,  C_Temp:C_Cross) = this%Cl_scalar(lmin:CP%Max_l, C_Temp:C_Cross) * Norm
    end if

    if (CP%WantTensors) then
        if (.not.CP%WantScalars) Norm = 1/this%Cl_tensor(lnorm, C_Temp)
        !Otherwise Norm already set correctly
        this%Cl_tensor(lmin:CP%Max_l_tensor, CT_Temp:CT_Cross) =  &
            this%Cl_tensor(lmin:CP%Max_l_tensor,  CT_Temp:CT_Cross) * Norm
    end if

    end subroutine TCLdata_NormalizeClsAtL


    end module CambSettings



    module Transfer
    use CAMBsettings
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
    !Sources
    !Alternatively for 21cm
    integer, parameter :: Transfer_monopole=4, Transfer_vnewt=5, Transfer_Tmat = 6

    integer, parameter :: Transfer_max = Transfer_vel_baryon_cdm
    character(LEN=name_tag_len) :: Transfer_name_tags(Transfer_max-1) = &
        ['CDM     ', 'baryon  ', 'photon  ', 'nu      ', 'mass_nu ', 'total   ', &
        'no_nu   ', 'total_de', 'Weyl    ', 'v_CDM   ', 'v_b     ', 'v_b-v_c ']
    character(LEN=name_tag_len) :: Transfer21cm_name_tags(Transfer_max-1) = &
        ['CDM      ', 'baryon   ', 'photon   ', 'monopole ', 'v_newt   ', 'delta_T_g', &
        'no_nu    ', 'total_de ', 'Weyl     ', 'v_CDM    ', 'v_b      ', 'v_b-v_c  ']

    logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
    !set to false to output calculated values for later interpolation

    integer :: transfer_power_var = Transfer_tot
    !What to use to calulcate the output matter power spectrum and sigma_8
    !Transfer_tot uses total matter perturbation

    logical :: get_growth_sigma8 = .true.
    !gets sigma_vdelta, like sigma8 but using velocity-density cross power,
    !in late LCDM f*sigma8 = sigma_vdelta^2/sigma8

    !Sources
    Type Cl21cmVars
        Type(MatterPowerData), pointer :: PK
        integer l, itf
        logical logs
        real(dl) chi
    end Type Cl21cmVars

    interface Transfer_GetMatterPower
    module procedure Transfer_GetMatterPowerD,Transfer_GetMatterPowerS
    end interface

    real(dl), private, external :: rombint_obj

    contains

    subroutine Transfer_GetUnsplinedPower(M, Params, PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k)
    Type(MatterTransferData) :: M
    Type(CambParams) :: Params
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    real(dl) :: h, k
    integer :: nz, nk, zix, ik, s1, s2
    logical :: hnorm

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)
    hnorm = DefaultTrue (hubble_units)

    nk=M%num_q_trans
    nz=CP%Transfer%PK_num_redshifts
    if (nk/= size(PK,1) .or. nz/=size(PK,2)) call MpiStop('Trasfer_GetUnsplinedPower wrong size')

    h = CP%H0/100

    do ik=1,nk
        k = M%TransferData(Transfer_kh,ik,1)*h
        do zix=1,nz
            PK(ik,zix) = M%TransferData(s1,ik,Params%Transfer%PK_redshifts_index(nz-zix+1))*&
                M%TransferData(s2,ik,Params%Transfer%PK_redshifts_index(nz-zix+1))*k*&
                const_pi*const_twopi*Params%InitPower%ScalarPower(k)
        end do
    end do
    if (hnorm) PK=  PK * h**3

    end subroutine Transfer_GetUnsplinedPower

    subroutine Transfer_GetUnsplinedNonlinearPower(M,Params,PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k) after re-scaling for non-linear evolution (if turned on)
    Type(MatterTransferData), intent(in) :: M
    Type(CambParams) :: Params
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    Type(MatterPowerData) :: PKdata
    integer zix

    call Transfer_GetUnsplinedPower(M,Params,PK,var1,var2, hubble_units)
    do zix=1, Params%Transfer%PK_num_redshifts
        call Transfer_GetMatterPowerData(M, Params, PKdata, &
            Params%Transfer%PK_redshifts_index(Params%Transfer%PK_num_redshifts-zix+1))
        call Params%NonLinearModel%GetNonLinRatios(PKdata)
        PK(:,zix) =  PK(:,zix) *PKdata%nonlin_ratio(:,1)**2
        call MatterPowerdata_Free(PKdata)
    end do

    end subroutine Transfer_GetUnsplinedNonlinearPower

    subroutine Transfer_GetMatterPowerData(MTrans, Params, PK_data, itf_only, var1, var2)
    !Does *NOT* include non-linear corrections
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !spectrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CambParams), intent(in) :: Params
    Type(MatterPowerData) :: PK_data
    integer, intent(in), optional :: itf_only
    integer, intent(in), optional :: var1, var2
    real(dl) :: h, kh, k, power
    integer :: ik, nz, itf, itf_start, itf_end, s1, s2

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

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
    PK_data%redshifts = Params%Transfer%Redshifts(itf_start:itf_end)

    h = Params%H0/100

    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        k = kh*h
        PK_data%log_kh(ik) = log(kh)
        power = Params%InitPower%ScalarPower(k)
        if (global_error_flag/=0) then
            call MatterPowerdata_Free(PK_data)
            return
        end if
        do itf = 1, nz
            PK_data%matpower(ik,itf) = &
                log(MTrans%TransferData(s1,ik,itf_start+itf-1)*&
                MTrans%TransferData(s2,ik,itf_start+itf-1)*k &
                *const_pi*const_twopi*h**3*power)
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
    use FileUtils
    character(LEN=*) :: fname
    type(MatterPowerData) :: PK_data
    type(TTextFile) :: F
    real(dl)kh, Pk
    integer ik
    integer nz

    nz = 1
    call F%Open(fname)

    PK_data%num_k = F%Lines()
    PK_Data%num_z = 1

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))

    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = 0

    do ik=1,PK_data%num_k
        read (F%unit, *) kh, Pk
        PK_data%matpower(ik,1) = log(Pk)
        PK_data%log_kh(ik) = log(kh)
    end do

    call MatterPowerdata_getsplines(PK_data)
    call F%Close()

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

    !Sources
    subroutine MatterPowerdata_getsplines21cm(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

    do i = 1,PK_Data%num_z
        call spline(PK_data%log_k,PK_data%matpower(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddmat(1,i))
        call spline(PK_data%log_k,PK_data%vvpower(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddvvpower(1,i))
        call spline(PK_data%log_k,PK_data%vdpower(1,i),PK_data%num_k,&
            cllo,clhi,PK_data%ddvdpower(1,i))
    end do

    end subroutine MatterPowerdata_getsplines21cm


    subroutine MatterPowerdata_MakeNonlinear(NonLinearModel,PK_data)
    class(TNonLinearModel) :: NonLinearModel
    Type(MatterPowerData) :: PK_data

    call NonLinearModel%GetNonLinRatios(PK_data)
    PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
    call MatterPowerdata_getsplines(PK_data)

    end subroutine MatterPowerdata_MakeNonlinear

    subroutine MatterPowerdata_Free(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    !this shouldn't be needed. But need to check for compiler bugs..
    deallocate(PK_data%log_kh,stat=i)
    deallocate(PK_data%matpower,stat=i)
    deallocate(PK_data%ddmat,stat=i)
    deallocate(PK_data%nonlin_ratio,stat=i)
    deallocate(PK_data%redshifts,stat=i)
    !Sources
    deallocate(PK_data%log_k,stat=i)
    deallocate(PK_data%nonlin_ratio_vv,stat=i)
    deallocate(PK_data%nonlin_ratio_vd,stat=i)
    deallocate(PK_data%vvpower,stat=i)
    deallocate(PK_data%ddvvpower,stat=i)
    deallocate(PK_data%vdpower,stat=i)
    deallocate(PK_data%ddvdpower,stat=i)

    call MatterPowerdata_Nullify(PK_data)
    end subroutine MatterPowerdata_Free

    subroutine MatterPowerdata_Nullify(PK_data)
    Type(MatterPowerData) :: PK_data

    nullify(PK_data%log_kh)
    nullify(PK_data%nonlin_ratio)
    nullify(PK_data%redshifts)
    !Sources
    nullify(PK_data%log_k)
    nullify(PK_data%nonlin_ratio_vv)
    nullify(PK_data%nonlin_ratio_vd)
    nullify(PK_data%vvpower)
    nullify(PK_data%ddvvpower)
    nullify(PK_data%vdpower)
    nullify(PK_data%ddvdpower)

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

    !Sources
    subroutine MatterPower21cm_k(PK,  k, itf, monopole, vv, vd)
    !Get monopole and velocity power at particular k by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: k
    real(dl), intent(out) :: monopole, vv, vd
    real(dl) :: logk
    integer llo,lhi
    real(dl) ho,a0,b0,f1,f2,f3
    integer, save :: i_last = 1

    logk = log(k)
    if (logk < PK%log_k(1)) then
        monopole = 0
        vv=0
        return
    end if
    if (logk > PK%log_k(PK%num_k)) then
        monopole=0
        vv=0
        return
        !stop 'MatterPower21cm_k: out of bounds'
    else
        llo=min(i_last,PK%num_k)
        do while (PK%log_k(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_k(llo+1)< logk)
            llo=llo+1
        end do
        i_last =llo
        lhi=llo+1
        ho=PK%log_k(lhi)-PK%log_k(llo)
        a0=(PK%log_k(lhi)-logk)/ho
        b0=1-a0
        f1= (a0**3-a0)
        f2= (b0**3-b0)
        f3= ho**2/6


        monopole = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            (f1* PK%ddmat(llo,itf) &
            +f2*PK%ddmat(lhi,itf))*f3
        vv = a0*PK%vvpower(llo,itf)+ b0*PK%vvpower(lhi,itf)+&
            (f1* PK%ddvvpower(llo,itf) &
            +f2*PK%ddvvpower(lhi,itf))*f3

        vd = a0*PK%vdpower(llo,itf)+ b0*PK%vdpower(lhi,itf)+&
            (f1* PK%ddvdpower(llo,itf) &
            +f2*PK%ddvdpower(lhi,itf))*f3
    end if

    monopole = exp(max(-30._dl,monopole))
    vv = exp(max(-30._dl,vv))
    vd = exp(max(-30._dl,vd))

    end subroutine MatterPower21cm_k


    subroutine Transfer_GetMatterPowerS(MTrans, Params, outpower, itf, minkh, dlnkh, npoints, var1, var2)
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CambParams) :: Params
    integer, intent(in) :: itf, npoints
    integer, intent(in), optional :: var1, var2
    real, intent(out) :: outpower(*)
    real, intent(in) :: minkh, dlnkh
    real(dl) :: outpowerd(npoints)
    real(dl):: minkhd, dlnkhd

    minkhd = minkh; dlnkhd = dlnkh
    call Transfer_GetMatterPowerD(MTrans, Params, outpowerd, itf, minkhd, dlnkhd, npoints,var1, var2)
    outpower(1:npoints) = outpowerd(1:npoints)

    end subroutine Transfer_GetMatterPowerS

    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    !Changed input variable from itf to itf_PK because we are looking for the itf_PK'th
    !redshift in the PK_redshifts array.  The position of this redshift in the master redshift
    !array, itf, is given by itf = CP%Transfer%Pk_redshifts_index(itf_PK)
    !Also changed (CP%NonLinear/=NonLinear_None) to
    !CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
    subroutine Transfer_GetMatterPowerD(MTrans, Params,outpower, itf_PK, minkh, dlnkh, npoints, var1, var2)
    !Allows for non-smooth priordial spectra
    !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
    !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
    !in units of (h Mpc^{-1})^3.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CambParams) :: Params
    Type(MatterPowerData) :: PK

    integer, intent(in) :: itf_PK, npoints
    real(dl), intent(out) :: outpower(npoints)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in), optional :: var1, var2

    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    integer ik, llo,il,lhi,lastix
    real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
    real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
    integer itf
    integer :: s1,s2

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

    itf = Params%Transfer%PK_redshifts_index(itf_PK)

    if (npoints < 2) call MpiStop('Need at least 2 points in Transfer_GetMatterPower')

    !         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
    !            stop 'Transfer_GetMatterPower: kh out of computed region'
    !          end if
    if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
        .and. FeedbackLevel > 0 ) &
        write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'


    if (Params%NonLinear/=NonLinear_none .and. Params%NonLinear/=NonLinear_Lens) then
        call Transfer_GetMatterPowerData(MTrans, Params, PK, itf) ! Mar 16, changed to use default variable
        call Params%NonLinearModel%GetNonLinRatios(PK)
    end if

    h = Params%H0/100
    logmink = log(minkh)
    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,itf)
        k = kh*h
        kvals(ik) = log(kh)
        atransfer=MTrans%TransferData(s1,ik,itf)*MTrans%TransferData(s2,ik,itf)
        if (Params%NonLinear/=NonLinear_none .and. Params%NonLinear/=NonLinear_Lens) &
            atransfer = atransfer* PK%nonlin_ratio(ik,1)**2 !only one element, this itf
        matpower(ik) = log(atransfer*k*const_pi*const_twopi*h**3)
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
        outpower(il) = outpower(il) * Params%InitPower%ScalarPower(k)
        if (global_error_flag /= 0) exit
    end do

    if (Params%NonLinear/=NonLinear_none .and. Params%NonLinear/=NonLinear_Lens) call MatterPowerdata_Free(PK)

    end subroutine Transfer_GetMatterPowerD

    subroutine Transfer_Get_SigmaR(MTrans, Params, R, outvals, var1, var2, root)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc, for all requested redshifts
    !set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    Type(MatterTransferData) :: MTrans
    Type(CambParams) :: Params
    real(dl), intent(in) :: R
    integer, intent(in), optional :: var1, var2
    logical, intent(in), optional :: root !if true, give sigma8, otherwise sigma8^2
    real(dl), intent(out) :: outvals(:)
    real(dl) :: kh, k, h, x, win, lnk, dlnk, lnko, powers
    real(dl), dimension(Params%Transfer%PK_num_redshifts) :: dsig8, dsig8o, sig8, sig8o
    integer :: s1, s2, ik

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)
    H=Params%h0/100._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    do ik=1, MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        if (kh==0) cycle
        k = kh*H

        dsig8 = MTrans%TransferData(s1,ik, &
            Params%Transfer%PK_redshifts_index(1:Params%Transfer%PK_num_redshifts))
        if (s1==s2) then
            dsig8 = dsig8**2
        else
            dsig8 = dsig8*MTrans%TransferData(s2,ik, &
                Params%Transfer%PK_redshifts_index(1:Params%Transfer%PK_num_redshifts))
        end if
        x= kh *R
        win =3*(sin(x)-x*cos(x))/x**3
        lnk=log(k)
        if (ik==1) then
            dlnk=0.5_dl
            !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
            !Contribution should be very small in any case
        else
            dlnk=lnk-lnko
        end if
        powers = Params%InitPower%ScalarPower(k)
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
    outvals(1:Params%Transfer%PK_num_redshifts) = sig8

    end subroutine Transfer_Get_SigmaR

    subroutine Transfer_GetSigmaRArray(MTrans, Params,R, sigmaR, redshift_ix, var1, var2)
    !Get array of SigmaR at (by default) redshift zero, for all values of R
    Type(MatterTransferData) :: MTrans
    Type(CambParams) :: Params
    real(dl), intent(in) :: R(:)
    real(dl), intent(out) :: SigmaR(:)
    integer, intent(in), optional :: redshift_ix, var1, var2
    integer red_ix, ik, subk
    real(dl) kh, k, h, dkh
    real(dl) lnk, dlnk, lnko, minR
    real(dl), dimension(size(R)) ::  x, win, dsig8, dsig8o, sig8, sig8o
    type(MatterPowerData) :: PKspline
    integer, parameter :: nsub = 5

    minR = minval(R)
    red_ix = PresentDefault (Params%Transfer%PK_redshifts_index(&
        Params%Transfer%PK_num_redshifts), redshift_ix)

    call Transfer_GetMatterPowerData(MTrans, Params, PKspline, red_ix, var1, var2 )

    H=Params%h0/100._dl
    dkh = 0._dl
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
                !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
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

    SigmaR=sqrt(sig8/(const_pi*const_twopi*h**3 ))

    end subroutine Transfer_GetSigmaRArray

    subroutine Transfer_Get_sigma8(MTrans, Params, R, var1, var2)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc
    ! set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    Type(MatterTransferData) :: MTrans
    Type(CambParams) :: Params
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var1, var2
    real(dl) :: radius

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)

    call Transfer_Get_SigmaR(MTrans, Params, radius, MTrans%sigma_8, var1,var2)

    end subroutine Transfer_Get_sigma8

    subroutine Transfer_Get_sigmas(MTrans, Params, R, var_delta, var_v)
    !Get sigma8 and sigma_{delta v} (for growth, like f sigma8 in LCDM)
    Type(MatterTransferData) :: MTrans
    Type(CambParams) :: Params
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var_delta, var_v
    real(dl) :: radius
    integer :: s1, s2

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)
    s1 = PresentDefault (transfer_power_var, var_delta)
    s2 = PresentDefault (Transfer_Newt_vel_cdm, var_v)

    call Transfer_Get_SigmaR(MTrans, Params, radius, MTrans%sigma_8, s1,s1)
    if (get_growth_sigma8) call Transfer_Get_SigmaR(MTrans, Params,radius, &
        MTrans%sigma2_vdelta_8(:), s1, s2, root=.false.)

    end subroutine Transfer_Get_sigmas

    subroutine Transfer_output_Sig8(MTrans)
    Type(MatterTransferData), intent(in) :: MTrans
    integer j
    !JD 08/13 Changes in here to PK arrays and variables
    integer j_PK

    do j_PK=1, CP%Transfer%PK_num_redshifts
        j = CP%Transfer%PK_redshifts_index(j_PK)
        write(*,'("at z =",f7.3," sigma8 (all matter) = ",f7.4)') &
            CP%Transfer%redshifts(j), MTrans%sigma_8(j_PK)
    end do
    if (get_growth_sigma8) then
        do j_PK=1, CP%Transfer%PK_num_redshifts
            j = CP%Transfer%PK_redshifts_index(j_PK)
            write(*,'("at z =",f7.3," sigma8^2_vd/sigma8  = ",f7.4)') &
                CP%Transfer%redshifts(j), MTrans%sigma2_vdelta_8(j_PK)/MTrans%sigma_8(j_PK)
        end do
    end if

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
    allocate(MTrans%sigma_8(CP%Transfer%PK_num_redshifts))
    if (get_growth_sigma8) allocate(MTrans%sigma2_vdelta_8(CP%Transfer%PK_num_redshifts))

    end subroutine Transfer_Allocate

    subroutine Transfer_Nullify(Mtrans)
    Type(MatterTransferData):: MTrans

    Mtrans%num_q_trans = 0
    nullify(MTrans%q_trans)
    nullify(MTrans%TransferData)
    nullify(MTrans%sigma_8)
    nullify(MTrans%sigma2_vdelta_8)
    nullify(MTrans%optical_depths)

    end subroutine Transfer_Nullify

    subroutine Transfer_Free(MTrans)
    Type(MatterTransferData):: MTrans
    integer st

    deallocate(MTrans%q_trans, STAT = st)
    deallocate(MTrans%TransferData, STAT = st)
    deallocate(MTrans%sigma_8, STAT = st)
    if (get_growth_sigma8) deallocate(MTrans%sigma2_vdelta_8, STAT = st)
    deallocate(MTrans%optical_depths, STAT = st)

    end subroutine Transfer_Free


    !JD 08/13 Changes for nonlinear lensing of CMB + MPK compatibility
    !Changed function below to write to only P%NLL_*redshifts* variables
    subroutine Transfer_SetForNonlinearLensing(P, eta_k_max)
    Type(TransferParams) :: P
    integer i
    real maxRedshift
    !Sources
    real(dl), intent(in), optional :: eta_k_max

    !Sources
    if (CP%Do21cm) then
        if (maxval(Redshift_w(1:num_redshiftwindows)%Redshift) /= minval(Redshift_w(1:num_redshiftwindows)%Redshift))  &
            stop 'Non-linear 21cm currently only for narrow window at one redshift'
        if (.not. present(eta_k_max)) stop 'bad call to Transfer_SetForNonlinearLensing'
        P%kmax = eta_k_max/10000.
        P%k_per_logint  = 0
        P%NLL_num_redshifts =  1
        P%NLL_redshifts(1) = Redshift_w(1)%Redshift
        return
    end if

    P%kmax = max(P%kmax,5*CP%Accuracy%AccuracyBoost)
    P%k_per_logint  = 0
    maxRedshift = 10
    P%NLL_num_redshifts =  nint(10*CP%Accuracy%AccuracyBoost)
    if (CP%Accuracy%AccuracyBoost>=2) then
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
    use constants
    Type(MatterTransferData), intent(in) :: MTrans
    integer i,ik
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    !JD 08/13 Changes in here to PK arrays and variables
    integer i_PK
    integer unit

    do i_PK=1, CP%Transfer%PK_num_redshifts
        if (FileNames(i_PK) /= '') then
            i = CP%Transfer%PK_redshifts_index(i_PK)
            if (CP%do21cm) then
                unit = open_file_header(FileNames(i_PK), 'k/h', Transfer21cm_name_tags, 14)
            else
                unit = open_file_header(FileNames(i_PK), 'k/h', transfer_name_tags, 14)
            end if
            do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                    write(unit,'(*(E15.6))') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
            end do
            close(unit)
        end if
    end do

    end subroutine Transfer_SaveToFiles

    subroutine Transfer_SaveMatterPower(MTrans, FileNames, all21cm)
    use constants
    !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
    Type(MatterTransferData), intent(in) :: MTrans
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    character(LEN=name_tag_len) :: columns(3)
    integer itf, i, unit
    integer points
    real, dimension(:,:), allocatable :: outpower
    real minkh,dlnkh
    Type(MatterPowerData) :: PK_data
    integer ncol
    logical, intent(in), optional :: all21cm
    logical all21
    !JD 08/13 Changes in here to PK arrays and variables
    integer itf_PK

    all21 = DefaultFalse(all21cm)
    if (all21) then
        ncol = 3
    else
        ncol = 1
    end if

    do itf=1, CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = CP%Transfer%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
                allocate(outpower(points,ncol))

                !Sources
                if (all21) then
                    call Transfer_Get21cmPowerData(MTrans, CP, PK_data, itf_PK)
                else
                    call Transfer_GetMatterPowerData(MTrans, CP, PK_data, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)&
                        call MatterPowerdata_MakeNonlinear(CP%NonLinearModel, PK_Data)
                end if

                outpower(:,1) = exp(PK_data%matpower(:,1))
                !Sources
                if (all21) then
                    outpower(:,3) = exp(PK_data%vvpower(:,1))
                    outpower(:,2) = exp(PK_data%vdpower(:,1))

                    outpower(:,1) = outpower(:,1)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,2) = outpower(:,2)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,3) = outpower(:,3)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                end if

                call MatterPowerdata_Free(PK_Data)
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData(Transfer_kh,i,1),outpower(i,:)
                end do
                close(unit)
            else
                if (all21) stop 'Transfer_SaveMatterPower: if output all assume not interpolated'
                minkh = 1e-4
                dlnkh = 0.02
                points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
                !             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
                allocate(outpower(points,1))
                call Transfer_GetMatterPowerS(MTrans, CP, outpower(1,1), itf,  minkh,dlnkh, points)

                columns(1) = 'P'
                unit = open_file_header(FileNames(itf), 'k/h', columns(:1), 15)

                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1)
                end do
                close(unit)
            end if

            deallocate(outpower)
        end if
    end do

    end subroutine Transfer_SaveMatterPower


    subroutine Transfer_Get21cmPowerData(MTrans, Params, PK_data, z_ix)
    !In terms of k, not k/h, and k^3 P_k /2pi rather than P_k
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CambParams) :: Params
    Type(MatterPowerData) :: PK_data, PK_cdm
    real(dl) h, k, pow
    integer ik
    integer z_ix,nz

    nz = 1
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%vvpower(PK_data%num_k,nz))
    allocate(PK_data%ddvvpower(PK_data%num_k,nz))
    allocate(PK_data%vdpower(PK_data%num_k,nz))
    allocate(PK_data%ddvdpower(PK_data%num_k,nz))
    allocate(PK_data%log_k(PK_data%num_k))
    allocate(PK_data%redshifts(nz))

    PK_data%redshifts = Params%Transfer%Redshifts(z_ix)

    h = Params%H0/100

    if (Params%NonLinear/=NonLinear_None .and. Params%NonLinear/=NonLinear_Lens) then
        if (z_ix>1) stop 'not tested more than one redshift with Nonlinear 21cm'
        call Transfer_GetMatterPowerData(MTrans, Params, PK_cdm, z_ix)
        call Params%NonLinearModel%GetNonLinRatios_All(PK_cdm)
    end if

    do ik=1,MTrans%num_q_trans
        k = MTrans%TransferData(Transfer_kh,ik,z_ix)*h
        PK_data%log_k(ik) = log(k)
        pow = Params%InitPower%ScalarPower(k)*1d10

        PK_data%matpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2)**2 * pow)
        PK_data%vvpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)**2 * pow)
        PK_data%vdpower(ik,1) = &
            log( abs((MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)*&
            (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2))* pow)

        if (Params%NonLinear/=NonLinear_None) then
            PK_data%matpower(ik,1) = PK_data%matpower(ik,1) + 2*log(PK_cdm%nonlin_ratio(ik,z_ix))
            PK_data%vvpower(ik,1) = PK_data%vvpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vv(ik,z_ix))
            PK_data%vdpower(ik,1) = PK_data%vdpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vd(ik,z_ix))
        end if

        !test              PK_data%matpower(ik,z_ix) = &
        !                    log( (MTrans%TransferData(Transfer_cdm,ik,z_ix)*k**2)**2 * pow)
        !                 PK_data%vvpower(ik,z_ix) = 0
    end do

    if (Params%NonLinear/=NonLinear_None)  call MatterPowerdata_Free(PK_cdm)


    call MatterPowerdata_getsplines21cm(PK_data)

    end subroutine Transfer_Get21cmPowerData

    function Get21cmCl_l(Vars,kin)
    !Direct integration with j^2/r, etc.
    Type(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,k, jlm1
    real(dl) Get21cmCl_l
    real(dl) monopole, vv , vd
    external BJL_EXTERNAL
    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK,  k, Vars%itf, monopole, vv, vd)
    call bjl_external(Vars%l, x, jl)
    call bjl_external(Vars%l-1, x, jlm1)
    ddjl = -( 2/x*jlm1-(Vars%l+2)*real(Vars%l+1,dl)/x**2*jl + jl)

    Get21cmCl_l = jl**2*monopole + ddjl**2*vv - 2._dl *ddjl*jl*vd
    if (.not. Vars%logs)  Get21cmCl_l =  Get21cmCl_l / k

    end function Get21cmCl_l


    function Get21cmCl_l_avg(Vars,kin)
    !Asymptotic results where we take <cos^2>=1/2 assuming smooth power spectrum
    Type(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,cross,k
    real(dl) Get21cmCl_l_avg
    real(dl) monopole, vv , vd,lphalf
    external BJL_EXTERNAL

    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK,  k, Vars%itf, monopole, vv, vd)
    lphalf=Vars%l+0.5_dl

    jl = 1/(2*x**2) /sqrt(1-(lphalf/x)**2)

    !  ddjl = (4/x**4+1)/(2*x**2)
    !
    ddjl = (x**4-2*x**2*lphalf**2+lphalf**4)/(x**4*sqrt(x**2-lphalf**2)*x)/2

    !    cross = (2-x**2)/(2*x**4)

    cross = (-x**2+lphalf**2)/(x**2*sqrt(x**2-lphalf**2)*x)/2

    Get21cmCl_l_avg = jl*monopole + ddjl*vv - 2._dl *cross*vd
    if (.not. Vars%logs)  Get21cmCl_l_avg =  Get21cmCl_l_avg / k

    !       Get21cmCl_l_avg=Get21cmCl_l_avg
    end function Get21cmCl_l_avg


    subroutine Transfer_Get21cmCls(MTrans, FileNames)
    use constants
    !Get 21cm C_l from sharp shell, using only monopole source and redshift distortions
    Type(MatterTransferData), intent(in) :: MTrans
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer itf,ik
    integer points
    character(LEN=name_tag_len), dimension(3), parameter :: Transfer_21cm_name_tags = &
        ['CL  ','P   ','P_vv']
    Type(MatterPowerData), target ::PK_data
    real(dl)  tol,atol, chi, Cl
    integer l, lastl, unit
    real(dl) k_min, k_max,k, avg_fac
    Type(Cl21cmVars) vars


    tol = 1e-5/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    do itf=1, CP%Transfer%num_redshifts
        if (FileNames(itf) /= '') then
            !print *, 'tau = ', MTrans%optical_depths(itf)
            chi = State%tau0-TimeOfz(State,CP%Transfer%redshifts(itf))

            points = MTrans%num_q_trans

            lastl=0

            call MatterPowerdata_Nullify(PK_data)
            call Transfer_Get21cmPowerData(MTrans, CP, PK_data, itf)

            unit = open_file_header(FileNames(itf), 'L', Transfer_21cm_name_tags, 8)

            do ik =1, points-1
                k =exp(PK_data%log_k(ik))
                l=nint(k*chi)
                !This is not an approximation, we are just chosing to sample l at values around (k samples)*chi

                if (l>1 .and. l/= lastl) then
                    lastl=l
                    Vars%l=l
                    Vars%chi = chi
                    Vars%PK => PK_data
                    Vars%itf = 1
                    Cl=0
                    atol = tol
                    avg_fac = 200
                    k_min = max(exp(PK_data%log_k(1)), k*(1-20*CP%Accuracy%AccuracyBoost/chi))
                    k_max = CP%Accuracy%AccuracyBoost*max(k*(1+avg_fac/chi), k*(1._dl+real(l,dl)**(-0.666_dl)))

                    if (k_max*chi < l+10) k_max = (l+10)/chi

                    Vars%logs = .false.
                    if (k_max < exp(PK_data%log_k(points))) then
                        !Integrate bessels properly
                        Cl = rombint_obj(Vars,Get21cmCl_l,k_min,k_max, atol, 25)

                        Vars%logs = .true.
                        k_min = log(k_max)


                        if (l>2e6) then
                            !In baryon damping
                            !                     Vars%logs = .false.
                            !                     atol = tol/10
                            !                     k_min = max(exp(PK_data%log_k(1)), k*(1-10/chi) )
                            !                     k_max = k*(1+100/chi)

                            k_max = min(log(5*k), PK_data%log_k(points))


                            !                    if (k_max < exp(PK_data%log_k(points))) then
                            !                      Cl = rombint_obj(Vars,Get21cmCl_l,k_min,k_max, atol, 25)
                            !                      Vars%logs = .true.
                            !                      k_min = log(k*(1+100/chi))
                            !                     k_max = min(log(3*k), PK_data%log_k(points))
                            !                    else
                            !                       k_max = exp(PK_data%log_k(points))
                            !                    end if
                        elseif (l>1e4) then
                            Vars%logs = .false.
                            k_min = k_max

                            k_max = min(k*35*CP%Accuracy%AccuracyBoost, exp(PK_data%log_k(points)))
                        else
                            !In white noise regime
                            k_max = min(log(max(0.3_dl,k)*18*CP%Accuracy%AccuracyBoost), PK_data%log_k(points))
                        end if

                        Cl = Cl+rombint_obj(Vars,Get21cmCl_l_avg,k_min,k_max, atol, 25)
                        !                   Cl = Cl+rombint_obj(Vars,Get21cmCl_l,k_min,k_max, atol, 25)
                    else
                        k_max = exp(PK_data%log_k(points))
                        Cl = rombint_obj(Vars,Get21cmCl_l,k_min,k_max, atol, 25)
                    end if


                    Cl=exp(-2*MTrans%optical_depths(itf))*const_fourpi*Cl* &
                        real(l,dl)*(l+1)/const_twopi/1d10

                    write (unit, '(1I8,3E15.5)') l, Cl, exp(PK_data%matpower(ik,1)/1d10), exp(PK_data%vvpower(ik,1)/1d10)
                end if
            end do

            close(unit)

            call MatterPowerdata_Free(PK_Data)
        end if
    end do

    end subroutine Transfer_Get21cmCls

    !JD 08/13 New function for nonlinear lensing of CMB + MPK compatibility
    !Build master redshift array from array of desired Nonlinear lensing (NLL)
    !redshifts and an array of desired Power spectrum (PK) redshifts.
    !At the same time fill arrays for NLL and PK that indicate indices
    !of their desired redshifts in the master redshift array.
    !Finally define number of redshifts in master array. This is usually given by:
    !P%num_redshifts = P%PK_num_redshifts + P%NLL_num_redshifts - 1.  The -1 comes
    !from the fact that z=0 is in both arrays (when non-linear is on)
    subroutine Transfer_SortAndIndexRedshifts(P)
    use MpiUtils, only : MpiStop
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
            call MpiStop('Transfer_SortAndIndexRedshifts: Too many redshifts')

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



