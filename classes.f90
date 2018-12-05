    module classes
    use precision
    use MpiUtils
    use interpolation, only : TSpline1D, TCubicSpline, TLogRegularCubicSpline
    use IniObjects
    implicit none

    Type MatterTransferData
        !Computed data
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        real(dl), dimension (:), allocatable :: q_trans
        real(dl), dimension (:), allocatable ::  sigma_8
        real(dl), dimension (:), allocatable ::  sigma2_vdelta_8 !growth from sigma_{v delta}
        real, dimension(:,:,:), allocatable :: TransferData
        !Sources
        real(dl), dimension(:), allocatable :: optical_depths
        !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
    end Type MatterTransferData

    Type MatterPowerData
        !everything is a function of k/h
        integer   ::  num_k, num_z
        real(dl), dimension(:), allocatable :: log_kh, redshifts
        !matpower is log(P_k)
        real(dl), dimension(:,:), allocatable :: matpower, ddmat
        !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
        !function of k and redshift NonLinearScaling(k_index,z_index)
        real(dl), dimension(:,:), allocatable :: nonlin_ratio
        !Sources
        real(dl), dimension(:), allocatable :: log_k
        real(dl), dimension(:,:), allocatable :: vvpower, ddvvpower
        real(dl), dimension(:,:), allocatable :: vdpower, ddvdpower

        real(dl), dimension(:,:), allocatable :: nonlin_ratio_vv
        real(dl), dimension(:,:), allocatable :: nonlin_ratio_vd

    end Type MatterPowerData

    !Classes that can be accessed from python and contain allocatable objects and arrays or other classes
    !In Python inherited from F2003Class (defined in baseconfig)
    !All python-accessible inherited classes must define SelfPointer, and use @fortran_class decorator in python
    !Allocatable objects can only contain instances of python-accessible classes so they can be identified from python.
    !Class could be abstract and SelfPointer deferred, but Fortran then doesn't allow called inherited "super()" methods implemented in abstract classes
    Type TPythonInterfacedClass
        integer :: self !Not used, just need to be able to get reference to first entry if whole class is allocatable
    contains
    !SelfPointer msust be overriden for each class to be referenced from python.
    !Gets a class pointer from an untyped pointer.
    procedure, nopass :: SelfPointer
    !PythonClass gets string of python class name; not actually used internally
    procedure, nopass :: PythonClass
    !Replace is for copying over this type with an instance of the same type
    !It must be defined for any class used as a non-allocatable python structure field that can be assigned to
    procedure :: Replace
    end type TPythonInterfacedClass

    Type PythonClassPtr
        class(TPythonInterfacedClass), pointer :: Ref => null()
    end type

    Type PythonClassAllocatable
        class(TPythonInterfacedClass), allocatable :: P
    end Type PythonClassAllocatable

    type, extends(TPythonInterfacedClass) :: TCambComponent
    contains
    procedure :: ReadParams => TCambComponent_ReadParams
    end type TCambComponent

    type, extends(TPythonInterfacedClass) :: TCAMBParameters
        !Actual type defined in modules.f90
    end type TCAMBParameters

    Type, extends(TPythonInterfacedClass) :: TCAMBCalculation
        !Actual type defined in modules.f90
    end type TCAMBCalculation

    type, extends(TCambComponent) :: TNonLinearModel
        real(dl) :: Min_kh_nonlinear  = 0.005_dl
    contains
    procedure :: GetNonLinRatios => TNonLinearModel_GetNonLinRatios
    procedure :: GetNonLinRatios_All => TNonLinearModel_GetNonLinRatios_All
    end type TNonLinearModel

    type, extends(TPythonInterfacedClass) :: TInitialPower
    contains
    procedure :: ScalarPower => TInitialPower_ScalarPower
    procedure :: TensorPower => TInitialPower_TensorPower
    procedure :: Init => TInitialPower_Init
    procedure :: ReadParams => TInitialPower_ReadParams
    procedure :: Effective_ns => TInitalPower_Effective_ns
    end type TInitialPower

    Type, extends(TInitialPower) :: TSplinedInitialPower
        class(TSpline1D), allocatable :: Pscalar, Ptensor
    contains
    procedure :: SetScalarTable => TSplinedInitialPower_SetScalarTable
    procedure :: SetTensorTable => TSplinedInitialPower_SetTensorTable
    procedure :: SetScalarLogRegular => TSplinedInitialPower_SetScalarLogRegular
    procedure :: SetTensorLogRegular => TSplinedInitialPower_SetTensorLogRegular
    procedure :: ScalarPower => TSplinedInitialPower_ScalarPower
    procedure :: TensorPower => TSplinedInitialPower_TensorPower
    procedure :: HasTensors => TSplinedInitialPower_HasTensors
    procedure, nopass :: PythonClass => TSplinedInitialPower_PythonClass
    procedure, nopass :: SelfPointer => TSplinedInitialPower_SelfPointer
    end Type TSplinedInitialPower

    contains

    subroutine TCambComponent_ReadParams(this, Ini)
    class(TCambComponent) :: this
    class(TIniFile), intent(in) :: Ini

    end subroutine TCambComponent_ReadParams

    function PythonClass()
    character(LEN=:), allocatable :: PythonClass

    PythonClass = ''
    error stop 'PythonClass Not implemented'

    end function PythonClass

    subroutine SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TPythonInterfacedClass), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    error stop 'Class must define SelfPointer function returning pointer to actual type'

    end subroutine SelfPointer

    subroutine Replace(this, replace_with)
    class(TPythonInterfacedClass), target :: this
    class(TPythonInterfacedClass) :: replace_with

    error stop 'Assignment not implemented for this class'

    end subroutine Replace

    subroutine TNonLinearModel_GetNonLinRatios(this,Params,CAMB_Pk)
    class(TNonLinearModel) :: this
    class(TCAMBParameters) :: Params
    type(MatterPowerData) :: CAMB_Pk
    error stop 'GetNonLinRatios Not implemented'
    end subroutine TNonLinearModel_GetNonLinRatios

    subroutine TNonLinearModel_GetNonLinRatios_All(this,Params,CAMB_Pk)
    class(TNonLinearModel) :: this
    class(TCAMBParameters), intent(in) :: Params
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


    subroutine TInitialPower_Init(this, Params, Omegak)
    class(TInitialPower) :: this
    class(TCAMBParameters), intent(in) :: Params
    real(dl), intent(in) :: Omegak
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

    subroutine TSplinedInitialPower_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TSplinedInitialPower), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSplinedInitialPower_SelfPointer

    logical function TSplinedInitialPower_HasTensors(this)
    class(TSplinedInitialPower) :: this

    TSplinedInitialPower_HasTensors = allocated(this%Ptensor)

    end function TSplinedInitialPower_HasTensors

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

    function TSplinedInitialPower_PythonClass()
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