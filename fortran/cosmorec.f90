    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using CosmoRec
    ! Author: Richard Shaw (CITA)
    !
    ! To use with the python wrapper add -fPIC to CCFLAGS in the CosmoRec Makefile (for gcc)
    !---------------------------------------------------------------------------------------------------
    ! 08.06.2012: added possibility to communicate Hubble (Jens Chluba)
    ! 12.06.2012: AL, changed interface to pass nnu directly; fixed spline extrapolation
    ! 09.01.2019: AL, updated for new class structure

    module CosmoRec
    use precision
    use constants
    use classes
    use MathUtils
    use results
    use config
    use Interpolation
    use MpiUtils, only : MpiStop
    implicit none
    private

    integer, parameter :: Nz = 10000
    real(dl), parameter :: zmax = 1d4
    real(dl), parameter :: zmin = 0._dl

    type, extends(TRecombinationModel) :: TCosmoRec

        integer :: runmode = 0
        real(dl) :: fdm  = 0._dl ! Dark matter annihilation efficiency

        ! Internal accuracy of CosmoRec (0 - normal, 3 - most accurate
        ! other values defined in CosmoRec.cpp source file)
        real(dl) :: accuracy = 0._dl
        ! Override CosmoRec diffusion iterations; negative uses CosmoRec batch default (2)
        integer :: diff_iteration_max = -1
        ! Overrides for selected CosmoRec atomic and radiative-transfer settings
        integer :: n_shells = -1
        integer :: n_shells_hei = -1
        integer :: flag_hi_absorption = -1
        integer :: n_s_2gamma = -1
        integer :: n_s_raman = -1
        !Internal data
        class(TRegularCubicSpline), allocatable :: xrec, tmrec
    contains
    procedure :: ReadParams => TCosmoRec_ReadParams
    procedure :: Validate => TCosmoRec_Validate
    procedure :: Init => TCosmoRec_init
    procedure :: x_e => TCosmoRec_xe
    procedure :: T_m => TCosmoRec_tm !baryon temperature
    procedure, nopass :: SelfPointer => TCosmoRec_SelfPointer
    end type TCosmoRec

    public TCosmoRec
    contains

    subroutine TCosmoRec_ReadParams(this, Ini)
    class(TCosmoRec) :: this
    class(TIniFile), intent(in) :: Ini

    call Ini%Read('cosmorec_runmode', this%runmode)
    call Ini%Read('cosmorec_accuracy', this%accuracy)
    call Ini%Read('cosmorec_fdm', this%fdm)
    this%diff_iteration_max = Ini%Read_Int('cosmorec_diff_iteration_max', this%diff_iteration_max)
    this%n_shells = Ini%Read_Int('cosmorec_n_shells', this%n_shells)
    this%n_shells_hei = Ini%Read_Int('cosmorec_n_shells_hei', this%n_shells_hei)
    this%flag_hi_absorption = Ini%Read_Int('cosmorec_flag_hi_absorption', this%flag_hi_absorption)
    this%n_s_2gamma = Ini%Read_Int('cosmorec_n_s_2gamma', this%n_s_2gamma)
    this%n_s_raman = Ini%Read_Int('cosmorec_n_s_raman', this%n_s_raman)

    end subroutine TCosmoRec_ReadParams


    subroutine TCosmoRec_Validate(this, OK)
    class(TCosmoRec),intent(in) :: this
    logical, intent(inout) :: OK

    if(this%runmode < 0 .or. this%runmode > 3) then
        write(*,*) "Invalid runmode for CosmoRec,"
        OK = .false.
    end if

    if(this%runmode < 2 .and. this%fdm > 1d-23) then
        write(*,*) "Dark matter annihilation rate too high. Will crash CosmoRec."
        OK = .false.
    end if

    if(this%accuracy < -1.0 .or. this%accuracy > 6.0) then
        write(*,*) "CosmoRec accuracy mode undefined."
        OK = .false.
    end if
    if(this%diff_iteration_max < -1) then
        write(*,*) "Invalid CosmoRec diffusion iteration override."
        OK = .false.
    end if
    if(this%n_shells < -1 .or. this%n_shells_hei < -1 .or. this%flag_hi_absorption < -1 .or. &
        this%n_s_2gamma < -1 .or. this%n_s_raman < -1) then
        write(*,*) "Invalid CosmoRec shell or radiative-transfer override."
        OK = .false.
    end if

    end subroutine TCosmoRec_Validate

    function TCosmoRec_tm(this,a)
    class(TCosmoRec) :: this
    real(dl), intent(in) :: a
    real(dl) TCosmoRec_tm
    real(dl) z

    z =1/a-1
    if (z >= this%tmrec%xmax) then
        TCosmoRec_tm = this%tmrec%F(Nz)*(1+z)/(1+this%tmrec%xmax)
    else if (z <= this%tmrec%xmin) then
        TCosmoRec_tm = this%tmrec%F(1)
    else
        TCosmoRec_tm = this%tmrec%Value(z)
    end if

    end function TCosmoRec_tm

    function TCosmoRec_xe(this,a)
    class(TCosmoRec) :: this
    real(dl), intent(in) :: a
    real(dl) TCosmoRec_xe
    real(dl) z

    z =1/a-1

    if (z >= this%xrec%xmax) then
        TCosmoRec_xe = this%xrec%F(Nz)
    else if (z <= this%xrec%xmin) then
        TCosmoRec_xe = this%xrec%F(1)
    else
        TCosmoRec_xe = this%xrec%Value(z)
    end if

    end function TCosmoRec_xe

    subroutine TCosmoRec_init(this, State, WantTSpin)
    use MiscUtils
    class(TCosmoRec), target :: this
    class(TCAMBdata), target :: State
    logical, intent(in), optional :: WantTSpin
    integer :: i, label
    real(dl), dimension(21) :: runpars
    procedure(obj_function) :: dtauda
    real(dl) OmegaB, OmegaC, OmegaK, h2
    real(dl), allocatable :: Hz(:), zrec(:), tmrec(:), xrec(:), tmp(:)
    external CosmoRec_calc_h_cpp

#ifdef MPI
    label = GetMpiRank()
#else
    label = 0
#endif
    ! Some feedback
    if (DefaultFalse(WantTSpin)) call MpiStop('CosmoRec does not support 21cm')

    select type(State)
    class is (CAMBdata)

        if (State%CP%Evolve_delta_xe) &
            call MpiStop('CosmoRec currently does not support evolving Delta x_e')
        h2 = (State%CP%H0/100)**2
        OmegaB = State%CP%ombh2/h2
        !These parameters are now redundant since using Hz array, just set to something
        Omegak = State%CP%omk
        OmegaC = State%CP%omch2/h2

        if (.not. allocated(this%xrec)) allocate(TRegularCubicSpline::this%xrec)
        if (.not. allocated(this%tmrec)) allocate(TRegularCubicSpline::this%tmrec)

        allocate(Hz(Nz), zrec(Nz), xrec(Nz), tmrec(Nz), tmp(Nz))

        ! Set runtime parameters using the CosmoRec.vX batch layout. The first 15 entries
        ! are the native vX API; the final entries are CAMB-specific overrides.
        runpars = -1._dl
        runpars(1) = this%fdm ! Set dark matter annihilation efficiency
        runpars(2) = this%accuracy
        runpars(3) = 0._dl ! Do not write CosmoRec outputs
        runpars(4) = 0._dl ! Use default A2s1s
        runpars(5) = 0._dl ! B0: no magnetic-field heating
        runpars(6) = -1._dl ! nB placeholder
        runpars(7) = 2._dl ! PMF mode placeholder; irrelevant for B0=0
        runpars(8) = 1._dl ! fine-structure constant scaling
        runpars(9) = 1._dl ! electron-mass scaling
        runpars(10) = 0._dl ! redshift-dependent scaling
        runpars(11) = 6._dl ! native vX rescaling mode placeholder
        runpars(12:15) = 0._dl ! Farhang et al. PCA terms
        runpars(16) = this%diff_iteration_max
        runpars(17) = this%n_shells
        runpars(18) = this%n_shells_hei
        runpars(19) = this%flag_hi_absorption
        runpars(20) = this%n_s_2gamma
        runpars(21) = this%n_s_raman

        ! Set redshifts to calculate at.
        do i=1,Nz
            zrec(i) = zmax - (i-1)*((zmax - zmin) / (Nz-1))
            Hz(i) = 1/dtauda(State,1/(1._dl+zrec(i))) &
                *(1._dl+zrec(i))**2/MPC_in_sec
        end do

        ! internal Hubble function of CosmoRec is used
        !call CosmoRec_calc_cpp(Recomb%runmode, runpars, &
        !     OmegaC, OmegaB, OmegaK, num_nu, h0inp, tcmb, yp, &
        !     zrec, xrec, tmrec, Nz, label)

        ! version which uses camb Hubble function
        call CosmoRec_calc_h_cpp(this%runmode, runpars, &
            OmegaC, OmegaB, OmegaK, State%CP%N_eff(), State%CP%H0, State%CP%tcmb, State%CP%yhe, &
            zrec, Hz, Nz, zrec, xrec, tmrec, Nz, label)

        !Init interpolation
        tmp =xrec(Nz:1:-1)
        call this%xrec%Init(zrec(Nz), zrec(1), Nz, tmp)
        tmp =tmrec(Nz:1:-1)
        call this%tmrec%Init(zrec(Nz), zrec(1), Nz, tmp)
    end select

    end subroutine TCosmoRec_init

    subroutine TCosmoRec_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TCosmoRec), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TCosmoRec_SelfPointer

    end module CosmoRec
