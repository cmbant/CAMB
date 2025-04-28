    module model
    use Precision
    use classes
    use SourceWindows
    use constants, only : COBE_CMBTemp, default_nnu
    use DarkEnergyInterface
    use MassiveNu
    use config
    use iso_c_binding
    implicit none

    integer, parameter :: outNone=1

    integer, parameter :: neutrino_hierarchy_normal = 1, neutrino_hierarchy_inverted = 2, neutrino_hierarchy_degenerate = 3

    integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
    !For CAMBparams%MassiveNuMethod
    !Nu_int: always integrate distribution function
    !Nu_trunc: switch to expansion in velocity once non-relativistic
    !Nu_approx: approximate scheme - good for CMB, but not formally correct and no good for matter power
    !Nu_best: automatically use mixture which is fastest and most accurate

    integer, parameter :: max_Nu = 5 !Maximum number of neutrino species
    integer, parameter :: max_transfer_redshifts = 256

    type TransferParams
        logical     ::  high_precision = .false.
        logical     ::  accurate_massive_neutrinos = .false.
        real(dl)    ::  kmax = 0.9_dl        !these are acutally q values, but same as k for flat
        integer     ::  k_per_logint =0
        integer     ::  PK_num_redshifts = 1
        real(dl)    ::  PK_redshifts(max_transfer_redshifts) = 0._dl
    end type TransferParams

    type AccuracyParams
        !Parameters for checking/changing overall accuracy
        !parameters equal to 1 corresponds to ~0.1% scalar C_l accuracy (at L>600)

        real(dl) :: AccuracyBoost =1._dl
        !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
        !Can also be used to improve speed significantly if less accuracy is required.
        !or improving accuracy for extreme models.
        !Note this does not increase lSamp%l sampling or massive neutrino q-sampling


        real(dl) :: lSampleBoost=1._dl
        !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation

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

        real(dl) :: TimeStepBoost = 1._dl !sampling timesteps

        real(dl) :: BackgroundTimeStepBoost = 1._dl !number of time steps for background thermal history interpolation

        real(dl) :: IntTolBoost = 1._dl !Tolerances for integrating differential equations

        real(dl) :: SourcekAccuracyBoost = 1._dl !Accuracy of k sampling for source time integration

        real(dl) :: IntkAccuracyBoost = 1._dl !Accuracy of k sampling for integration

        real(dl) :: TransferkBoost = 1._dl !Accuracy of k sampling for transfer functions

        real(dl) :: NonFlatIntAccuracyBoost = 1._dl !Accuracy of non-flat time integration

        real(dl) :: BessIntBoost = 1._dl !Accuracy of bessel integration truncation

        real(dl) :: LensingBoost = 1._dl !Accuracy for CMB lensing of CMB power spectra

        real(dl) :: NonlinSourceBoost = 1._dl !Steps and kmax for the non-linear correction

        real(dl) :: BesselBoost = 1._dl !accuracy of bessel pre-computation sampling

        real(dl) :: LimberBoost = 1._dl !Accuracy of Limber approximation use

        real(dl) :: SourceLimberBoost = 1._dl !Scales when to switch to Limber for source windows

        real(dl) :: KmaxBoost = 1._dl !Boost max k for source window functions

        real(dl) :: neutrino_q_boost = 1._dl !number of momenta integrated for neutrino perturbations

    end type AccuracyParams

    !For holding custom CMB source functions (use via python interface)
    Type TCustomSourceParams
        integer :: num_custom_sources = 0
        Type(C_FUNPTR) :: c_source_func = c_null_funptr
        integer, allocatable :: custom_source_ell_scales(:)
    end Type TCustomSourceParams

    !Non-linear corrections, either just P(k), or just CMB lensing/sources, or both
    integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2
    integer, parameter :: NonLinear_both=3

    ! Main parameters type
    type, extends (TCAMBParameters) :: CAMBparams
        logical   :: WantCls  = .true.
        logical   :: WantTransfer = .false.

        logical   :: WantScalars = .true.
        logical   :: WantTensors = .false.
        logical   :: WantVectors = .false.

        logical   :: WantDerivedParameters = .true.!calculate various derived parameters  (ThermoDerivedParams)
        logical   :: Want_cl_2D_Array = .true.
        logical   :: Want_CMB = .true.
        logical   :: Want_CMB_lensing = .true.
        logical   :: DoLensing  = .true.
        integer   :: NonLinear = NonLinear_none
        type(TransferParams) :: Transfer

        logical   :: want_zstar = .false.
        logical   :: want_zdrag = .false.     !!JH for updated BAO likelihood.

        integer   :: Min_l = 2 ! 1 or larger, usually 1 or 2
        integer   :: Max_l = 2500
        integer   :: Max_l_tensor = 600
        real(dl)  :: Max_eta_k = 5000
        real(dl)  :: Max_eta_k_tensor = 1200
        ! _tensor settings only used in initialization,
        !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

        real(dl)  :: ombh2 = 0._dl !baryon density Omega_b h^2
        real(dl)  :: omch2 = 0._dl !cold dark matter density Omega_c h^2
        real(dl)  :: omk = 0._dl !Omega_K
        real(dl)  :: omnuh2 = 0._dl !massive neutino Omega_nu h^2
        real(dl)  :: H0 = 67._dl !Hubble parameter in km/s/Mpc
        real(dl)  :: TCMB = COBE_CMBTemp
        real(dl)  :: Yhe = 0.24_dl
        real(dl)  :: Num_Nu_massless = default_nnu
        integer   :: Num_Nu_massive = 0 !sum of Nu_mass_numbers below
        integer   :: Nu_mass_eigenstates = 0  !1 for degenerate masses
        logical   :: share_delta_neff = .false. !take fractional part to heat all eigenstates the same
        real(dl)  :: Nu_mass_degeneracies(max_nu)
        real(dl)  :: Nu_mass_fractions(max_nu) !The ratios of the total densities
        integer   :: Nu_mass_numbers(max_nu) !physical number per eigenstate

        class(TInitialPower), allocatable :: InitPower
        class(TRecombinationModel), allocatable :: Recomb
        class(TReionizationModel), allocatable :: Reion
        class(TDarkEnergyModel), allocatable :: DarkEnergy
        class(TNonLinearModel), allocatable :: NonLinearModel
        type(AccuracyParams)     :: Accuracy
        type(SourceTermParams)   :: SourceTerms

        real(dl), allocatable :: z_outputs(:) !Redshifts to output background outputs

        integer   :: Scalar_initial_condition = 1 !adiabatic
        !must be one of the initial_xxx values defined in GaugeInterface
        real(dl), allocatable  :: InitialConditionVector(:) !ignored unless Scalar_initial_condition == initial_vector

        integer   :: OutputNormalization = outNone
        !outNone, or C_OutputNormalization=1 if > 1

        real(dl)  :: Alens = 1._dl !Unphysical rescaling parameter of the CMB lensing power

        integer   :: MassiveNuMethod = Nu_best

        logical :: DoLateRadTruncation = .true.
        !if true, use smooth approx to radition perturbations after decoupling on
        !small scales, saving evolution of irrelevant osciallatory multipole equations

        logical :: Evolve_baryon_cs = .false.
        !if true, evolves equation for Delta_{T_m} to get cs_2 = \delta p /\delta\rho for perfect gas

        logical :: Evolve_delta_xe = .false. !Include ionization fraction perturbations

        logical :: Evolve_delta_Ts =.false. !Equilibrium result agrees to sub-percent level

        logical :: Do21cm = .false.
        logical :: transfer_21cm_cl = .false.
        logical :: Log_lvalues  = .false. !useful for smooth results at very high L
        logical :: use_cl_spline_template = .true.    
        integer :: min_l_logl_sampling = 5000 ! increase to use linear sampling for longer
    
        Type(TSourceWindowHolder), allocatable :: SourceWindows(:)

        Type(TCustomSourceParams) :: CustomSources
    contains
    procedure, nopass :: PythonClass => CAMBparams_PythonClass
    procedure, nopass :: SelfPointer => CAMBparams_SelfPointer
    procedure :: Replace => CAMBParams_Replace
    procedure :: SetNeutrinoHierarchy => CAMBparams_SetNeutrinoHierarchy
    procedure :: Validate => CAMBparams_Validate
    procedure :: PrimordialPower => CAMBparams_PrimordialPower
    procedure :: SetCustomSourcesFunc => CAMBparams_SetCustomSourcesFunc
    procedure :: N_eff => CAMBparams_N_eff
    end type CAMBparams

    contains

    function CAMBparams_PythonClass()
    character(LEN=:), allocatable :: CAMBparams_PythonClass
    CAMBparams_PythonClass = 'CAMBparams'
    end function CAMBparams_PythonClass

    subroutine CAMBparams_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (CAMBparams), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine CAMBparams_SelfPointer

    subroutine CAMBparams_Replace(this, replace_with)
    class(CAMBParams), target :: this
    Type(CAMBParams), pointer :: p
    class(TPythonInterfacedClass) :: replace_with

    select type(this)
    type is (CAMBParams)
        p => this
        select type(replace_with)
        type is (CAMBParams)
            p = replace_with
        class default
            error stop 'Wrong assignment type'
        end select
    end select
    end subroutine CAMBparams_Replace

    subroutine CAMBparams_SetNeutrinoHierarchy(this, omnuh2, omnuh2_sterile, nnu, neutrino_hierarchy, num_massive_neutrinos)
    !Set neutrino hierarchy in the approximate two-eigenstate model (treating two as exactly degenerate, and assuming non-relativistic),
    !or use degenerate mass approximation.
    !omnuh2 is the massive total neutrino density today, omnuh2_sterile is the component of that due to steriles
    !omnuh2_sterile is interpreted as in the Planck parameter papers
    use MathUtils
    use constants
    class(CAMBparams), intent(inout) :: this
    real(dl), intent(in) :: omnuh2, omnuh2_sterile, nnu
    integer, intent(in) :: neutrino_hierarchy
    integer, intent(in), optional :: num_massive_neutrinos  !for degenerate hierarchy
    real(dl) normal_frac, m3, neff_massive_standard, mnu, m1

    this%omnuh2 = omnuh2
    if (omnuh2==0) then
        this%Num_Nu_Massless = nnu
        return
    end if
    this%Nu_mass_eigenstates=0
    this%share_delta_neff = .false.
    if (omnuh2 > omnuh2_sterile) then
        normal_frac =  (omnuh2-omnuh2_sterile)/omnuh2
        if (neutrino_hierarchy == neutrino_hierarchy_degenerate) then
            neff_massive_standard = num_massive_neutrinos*default_nnu/3
            this%Num_Nu_Massive = num_massive_neutrinos
            this%Nu_mass_eigenstates=this%Nu_mass_eigenstates+1
            if (nnu > neff_massive_standard) then
                this%Num_Nu_Massless = nnu - neff_massive_standard
            else
                this%Num_Nu_Massless = 0
                neff_massive_standard=nnu
            end if
            this%Nu_mass_numbers(this%Nu_mass_eigenstates) = num_massive_neutrinos
            this%Nu_mass_degeneracies(this%Nu_mass_eigenstates) = neff_massive_standard
            this%Nu_mass_fractions(this%Nu_mass_eigenstates) = normal_frac
        else
            !Use normal or inverted hierarchy, approximated as two eigenstates in physical regime, 1 at minimum and below
            mnu = (omnuh2 - omnuh2_sterile)*neutrino_mass_fac*(COBE_CMBTemp/this%TCMB)**3/ (default_nnu / 3) ** 0.75_dl
            if (neutrino_hierarchy == neutrino_hierarchy_normal) then
                if (mnu > mnu_min_normal + 1e-4_dl) then
                    !Two eigenstate approximation.
                    m1=Newton_Raphson2(0._dl, mnu, sum_mnu_for_m1, mnu, 1._dl)
                    this%Num_Nu_Massive = 3
                else
                    !One eigenstate
                    this%Num_Nu_Massive = 1
                end if
            else if (neutrino_hierarchy == neutrino_hierarchy_inverted) then
                if (mnu > sqrt(delta_mnu31)+sqrt(delta_mnu31+delta_mnu21) + 1e-4_dl ) then
                    !Valid case, two eigenstates
                    m1=Newton_Raphson2(sqrt(delta_mnu31)+1e-6_dl, mnu, sum_mnu_for_m1, mnu, -1._dl)
                    this%Num_Nu_Massive = 3
                else
                    !Unphysical low mass case: take one (2-degenerate) eigenstate
                    this%Num_Nu_Massive = 2
                end if
            else
                error stop 'Unknown neutrino_hierarchy setting'
            end if
            neff_massive_standard = this%Num_Nu_Massive *default_nnu/3
            if (nnu > neff_massive_standard) then
                this%Num_Nu_Massless = nnu - neff_massive_standard
            else
                this%Num_Nu_Massless = 0
                neff_massive_standard=nnu
            end if
            if (this%Num_Nu_Massive==3) then
                !two with mass m1, one with m3
                this%Nu_mass_eigenstates = 2
                this%Nu_mass_degeneracies(1) = neff_massive_standard*2/3._dl
                this%Nu_mass_degeneracies(2) = neff_massive_standard*1/3._dl
                m3 = mnu - 2*m1
                this%Nu_mass_fractions(1) = 2*m1/mnu*normal_frac
                this%Nu_mass_fractions(2) = m3/mnu*normal_frac
                this%Nu_mass_numbers(1) = 2
                this%Nu_mass_numbers(2) = 1
            else
                this%Nu_mass_degeneracies(1) = neff_massive_standard
                this%Nu_mass_numbers(1) = this%Num_Nu_Massive
                this%Nu_mass_eigenstates = 1
                this%Nu_mass_fractions(1) = normal_frac
            end if
        end if
    else
        neff_massive_standard=0
    end if
    if (omnuh2_sterile>0) then
        if (nnu<default_nnu) call MpiStop('nnu < 3.044 with massive sterile')
        this%Num_Nu_Massless = default_nnu - neff_massive_standard
        this%Num_Nu_Massive=this%Num_Nu_Massive+1
        this%Nu_mass_eigenstates=this%Nu_mass_eigenstates+1
        this%Nu_mass_numbers(this%Nu_mass_eigenstates) = 1
        this%Nu_mass_degeneracies(this%Nu_mass_eigenstates) = max(1d-6,nnu - default_nnu)
        this%Nu_mass_fractions(this%Nu_mass_eigenstates) = omnuh2_sterile/omnuh2
    end if
    end subroutine CAMBparams_SetNeutrinoHierarchy

    real(dl) function CAMBparams_N_eff(this)
    class(CAMBparams), intent(in) :: this

    CAMBparams_N_eff = this%Num_Nu_Massless+sum(this%Nu_mass_degeneracies(1:this%Nu_mass_eigenstates))

    end function CAMBparams_N_eff

    function CAMBparams_Validate(this) result(OK)
    class(CAMBparams), intent(in) :: this
    logical OK

    OK = .true.
    if (.not. this%WantTransfer .and. .not. this%WantCls) then
        OK = .false.
        write(*,*) 'There is nothing to do! Do transfer functions or Cls.'
    end if

    if (this%h0 < 20._dl.or.this%h0 > 100._dl) then
        OK = .false.
        write(*,*) '  Warning: H0 has units of km/s/Mpc. You have:', this%h0
    end if
    if (this%tcmb < 2.7d0.or.this%tcmb > 2.8d0) then
        write(*,*) '  Warning: Tcmb has units of K.  Your have:', this%tcmb
    end if

    if (this%yhe < 0.2d0.or.this%yhe > 0.8d0) then
        OK = .false.
        write(*,*) &
            '  Warning: YHe is the Helium fraction of baryons.', &
            '  Your have:', this%yhe
    end if
    if (this%Num_Nu_massive < 0) then
        OK = .false.
        write(*,*) &
            'Warning: Num_Nu_massive is strange:',this%Num_Nu_massive
    end if
    if (this%Num_Nu_massless < 0) then
        OK = .false.
        write(*,*) &
            'Warning: Num_nu_massless is strange:', this%Num_Nu_massless
    end if
    if (this%Num_Nu_massive < 1 .and. this%omnuh2 > 0.0) then
        OK = .false.
        write(*,*) &
            'Warning: You have omega_neutrino > 0, but no massive species'
    end if


    if (this%ombh2<0.0005 .or. this%omch2<0 .or. this%ombh2>0.5 .or. this%omch2>2) then
        OK = .false.
        write(*,*) 'Your matter densities are strange (may not have been set)'
    end if

    if (this%WantScalars .and. this%Max_eta_k < this%Max_l .or.  &
        this%WantTensors .and. this%Max_eta_k_tensor < this%Max_l_tensor) then
        OK = .false.
        write(*,*) 'You need Max_eta_k larger than Max_l to get good results'
    end if

    call this%Reion%Validate(OK)
    call this%Recomb%Validate(OK)

    if (this%WantTransfer) then
        if (this%transfer%PK_num_redshifts > max_transfer_redshifts .or. this%transfer%PK_num_redshifts<1) then
            OK = .false.
            write(*,*) 'Maximum ',  max_transfer_redshifts, &
                'redshifts. You have: ', this%transfer%PK_num_redshifts
        end if
        if (this%transfer%kmax < 0.01 .or. this%transfer%kmax > 50000 .or. &
            this%transfer%k_per_logint>0 .and.  this%transfer%k_per_logint <1) then
            !            OK = .false.
            write(*,*) 'Strange transfer function settings.'
        end if
    end if

    end function CAMBparams_Validate


    subroutine CAMBParams_SetCustomSourcesFunc(this, ncustomsources, c_source_func, ell_scales)
    class(CAMBParams) :: this
    integer, intent(in) :: ncustomsources
    integer, intent(in) :: ell_scales(ncustomsources)
    TYPE(C_FUNPTR), INTENT(IN) :: c_source_func

    this%CustomSources%num_custom_sources = ncustomsources
    if (allocated(this%CustomSources%custom_source_ell_scales)) deallocate(this%CustomSources%custom_source_ell_scales)
    if (ncustomsources > 0) then
        this%CustomSources%c_source_func = c_source_func
        allocate(this%CustomSources%custom_source_ell_scales(ncustomsources), source=ell_scales)
    else
        this%CustomSources%c_source_func = c_null_funptr
    end if

    end subroutine CAMBParams_SetCustomSourcesFunc

    function CAMBparams_PrimordialPower(this, k, powers, n,  i) result(err)
    class(CAMBparams) :: this
    integer, intent(in) :: i,n
    real(dl), intent(in) :: k(n)
    real(dl), intent(out) :: powers(n)
    integer err,ix

    global_error_flag = 0
    call this%InitPower%Init(this)
    if (global_error_flag==0) then
        do ix =1, n
            if (i==0) then
                powers(ix) = this%InitPower%ScalarPower(k(ix))
            elseif (i==2) then
                powers(ix) = this%InitPower%TensorPower(k(ix))
            else
                error stop 'Unknown power type index'
            end if
            if (global_error_flag /= 0) exit
        end do
    end if
    err= global_error_flag

    end function CAMBparams_PrimordialPower

    end module model
