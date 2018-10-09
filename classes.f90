    module classes
    use precision
    use MpiUtils
    use interpolation, only : TSpline1D, TCubicSpline, TLogRegularCubicSpline
    use IniObjects
    implicit none

    Type MatterTransferData
        !Computed data
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        real(dl), dimension (:), pointer :: q_trans => NULL()
        real(dl), dimension (:), pointer ::  sigma_8 => NULL()
        real(dl), dimension (:), pointer ::  sigma2_vdelta_8 => NULL() !growth from sigma_{v delta}
        real, dimension(:,:,:), pointer :: TransferData => NULL()
        !Sources
        real(dl), dimension(:), pointer :: optical_depths => NULL()
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
        !Sources
        real(dl), dimension(:), pointer :: log_k => NULL()
        real(dl), dimension(:,:), pointer :: vvpower => NULL(), ddvvpower => NULL()
        real(dl), dimension(:,:), pointer :: vdpower => NULL(), ddvdpower => NULL()

        real(dl), dimension(:,:), pointer :: nonlin_ratio_vv => NULL()
        real(dl), dimension(:,:), pointer :: nonlin_ratio_vd => NULL()

    end Type MatterPowerData

    type TCambComponent
    contains
    procedure :: ReadParams => TCambComponent_ReadParams
    end type TCambComponent

    type, extends(TCambComponent) :: TNonLinearModel
        real(dl) :: Min_kh_nonlinear  = 0.005_dl
    contains
    procedure :: GetNonLinRatios => TNonLinearModel_GetNonLinRatios
    procedure :: GetNonLinRatios_All => TNonLinearModel_GetNonLinRatios_All
    end type TNonLinearModel

    type TInitialPower
        real(dl) :: curv = 0._dl !curvature parameter
    contains
    procedure :: ScalarPower => TInitialPower_ScalarPower
    procedure :: TensorPower => TInitialPower_TensorPower
    procedure :: PythonClass => TInitialPower_PythonClass
    procedure :: Init => TInitialPower_Init
    procedure :: ReadParams => TInitialPower_ReadParams
    procedure :: Effective_ns => TInitalPower_Effective_ns
    end type TInitialPower

    Type, extends(TInitialPower) :: TInitialPower2
        real(dl) :: curv2 = 0._dl
    contains
    end type

    Type, extends(TInitialPower) :: TSplinedInitialPower
        class(TSpline1D), allocatable :: Pscalar, Ptensor
    contains
    procedure :: SetScalarTable => TSplinedInitialPower_SetScalarTable
    procedure :: SetTensorTable => TSplinedInitialPower_SetTensorTable
    procedure :: SetScalarLogRegular => TSplinedInitialPower_SetScalarLogRegular
    procedure :: SetTensorLogRegular => TSplinedInitialPower_SetTensorLogRegular
    procedure :: ScalarPower => TSplinedInitialPower_ScalarPower
    procedure :: TensorPower => TSplinedInitialPower_TensorPower
    procedure :: PythonClass => TSplinedInitialPower_PythonClass
    end Type TSplinedInitialPower

    contains

    subroutine TCambComponent_ReadParams(this, Ini)
    class(TCambComponent), intent(inout) :: this
    class(TIniFile), intent(in) :: Ini

    end subroutine TCambComponent_ReadParams

    subroutine TNonLinearModel_GetNonLinRatios(this,CAMB_Pk)
    class(TNonLinearModel) :: this
    type(MatterPowerData) :: CAMB_Pk
    error stop 'GetNonLinRatios Not implemented'
    end subroutine TNonLinearModel_GetNonLinRatios

    subroutine TNonLinearModel_GetNonLinRatios_All(this,CAMB_Pk)
    class(TNonLinearModel) :: this
    type(MatterPowerData) :: CAMB_Pk
    error stop 'GetNonLinRatios_all  not supported (no non-linear velocities)'
    end subroutine TNonLinearModel_GetNonLinRatios_All


    function TInitialPower_ScalarPower(this, k)
    class(TInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TInitialPower_ScalarPower

    TInitialPower_ScalarPower = 0
    error stop 'ScalarPower not implemented'
    end function TInitialPower_ScalarPower

    function TInitialPower_TensorPower(this, k)
    class(TInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TInitialPower_TensorPower

    TInitialPower_TensorPower = 0
    error stop 'TensorPower not implemented'
    end function TInitialPower_TensorPower

    function TInitialPower_PythonClass(this)
    class(TInitialPower) :: this
    character(LEN=:), allocatable :: TInitialPower_PythonClass

    TInitialPower_PythonClass = ''

    end function TInitialPower_PythonClass

    subroutine TInitialPower_Init(this, acurv)
    class(TInitialPower) :: this
    real(dl), intent(in) :: acurv !Curvature parameter if non-flat

    this%curv = acurv
    end subroutine TInitialPower_Init


    subroutine TInitialPower_ReadParams(this, Ini, WantTensors)
    class(TInitialPower) :: this
    class(TIniFile), intent(in) :: Ini
    logical, intent(in) :: WantTensors

    end subroutine TInitialPower_ReadParams

    function TInitalPower_Effective_ns(this)
    class(TInitialPower) :: this
    real(dl) :: TInitalPower_Effective_ns

    TInitalPower_Effective_ns=1
    call MpiStop('Power spectrum does not support an effective n_s for halofit')

    end function TInitalPower_Effective_ns


    function TSplinedInitialPower_ScalarPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_ScalarPower

    TSplinedInitialPower_ScalarPower = this%Pscalar%Value(k)

    end function TSplinedInitialPower_ScalarPower

    function TSplinedInitialPower_TensorPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_TensorPower

    TSplinedInitialPower_TensorPower = this%Ptensor%Value(k)

    end function TSplinedInitialPower_TensorPower

    function TSplinedInitialPower_PythonClass(this)
    class(TSplinedInitialPower) :: this
    character(LEN=:), allocatable :: TSplinedInitialPower_PythonClass

    TSplinedInitialPower_PythonClass = 'SplinedInitialPower'

    end function TSplinedInitialPower_PythonClass

    subroutine TSplinedInitialPower_SetScalarTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
        allocate(TCubicSpline::this%Pscalar)
        select type (Sp => this%Pscalar)
        class is (TCubicSpline)
            call Sp%Init(k,PK)
        end select
    end if

    end subroutine TSplinedInitialPower_SetScalarTable
    
    
    subroutine TSplinedInitialPower_SetTensorTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%PTensor)) deallocate(this%PTensor)
    if (n>0) then
        allocate(TCubicSpline::this%PTensor)
        select type (Sp => this%PTensor)
        class is (TCubicSpline)
            call Sp%Init(k,PK)
        end select
    end if

    end subroutine TSplinedInitialPower_SetTensorTable

    subroutine TSplinedInitialPower_SetScalarLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
        allocate(TLogRegularCubicSpline::this%Pscalar)
        select type (Sp => this%Pscalar)
        class is (TLogRegularCubicSpline)
            call Sp%Init(kmin, kmax, n, PK)
        end select
    end if

    end subroutine TSplinedInitialPower_SetScalarLogRegular


    subroutine TSplinedInitialPower_SetTensorLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Ptensor)) deallocate(this%Ptensor)
    if (n>0) then
        allocate(TLogRegularCubicSpline::this%Ptensor)
        select type (Sp => this%Ptensor)
        class is (TLogRegularCubicSpline)
            call Sp%Init(kmin, kmax, n, PK)
        end select
    end if

    end subroutine TSplinedInitialPower_SetTensorLogRegular


    end module classes