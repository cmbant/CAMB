    !This module provides the initial power spectra

    !TInitialPowerLaw is parameterized as an expansion in ln k
    !
    ! ln P_s = ln A_s + (n_s -1)*ln(k/k_0_scalar) + n_{run}/2 * ln(k/k_0_scalar)^2 + n_{runrun}/6 * ln(k/k_0_scalar)^3
    !
    ! so if n_{run} = 0, n_{runrun}=0
    !
    ! P_s = A_s (k/k_0_scalar)^(n_s-1)
    !
    !for the scalar spectrum, when n_s=an(in) is the in'th spectral index. k_0_scalar
    !is a pivot scale, fixed here to 0.05/Mpc (change it below as desired or via .ini file).
    !
    !The tensor spectrum has three different supported parameterizations giving
    !
    ! ln P_t = ln A_t + n_t*ln(k/k_0_tensor) + n_{t,run}/2 * ln(k/k_0_tensor)^2
    !
    ! tensor_parameterization==tensor_param_indeptilt (=1) (default, same as CAMB pre-April 2014)
    !
    ! A_t = r A_s
    !
    ! tensor_parameterization==tensor_param_rpivot (=2)
    !
    ! A_t = r P_s(k_0_tensor)
    !
    ! tensor_parameterization==tensor_param_AT (=3)
    !
    ! A_t =  tensor_amp
    !
    !The absolute normalization of the Cls is unimportant here, but the relative ratio
    !of the tensor and scalar Cls generated with this module will be correct for general models

    module InitialPower
    use Precision
    use MpiUtils, only : MpiStop
    use classes
    implicit none

    private

    integer, parameter, public :: tensor_param_indeptilt=1,  tensor_param_rpivot = 2, tensor_param_AT = 3

    Type, extends(TInitialPower) :: TInitialPowerLaw
        integer :: tensor_parameterization = tensor_param_indeptilt
        !For the default implementation return power spectra based on spectral indices
        real(dl) :: ns = 1._dl !scalar spectral indices
        real(dl) :: nrun = 0._dl !running of spectral index
        real(dl) :: nrunrun  = 0._dl !running of spectral index
        real(dl) :: nt  = 0._dl !tensor spectral indices
        real(dl) :: ntrun  = 0._dl !tensor spectral index running
        real(dl) :: r  = 0._dl !ratio of scalar to tensor initial power spectrum amplitudes
        real(dl) :: pivot_scalar = 0.05_dl !pivot scales in Mpc^{-1}
        real(dl) :: pivot_tensor = 0.05_dl
        real(dl) :: As = 1._dl
        real(dl) :: At = 1._dl !A_T at k_0_tensor if tensor_parameterization==tensor_param_AT
        real(dl), private :: curv = 0._dl !curvature parameter
    contains
    procedure :: Init => TInitialPowerLaw_Init
    procedure, nopass :: PythonClass => TInitialPowerLaw_PythonClass
    procedure, nopass :: SelfPointer => TInitialPowerLaw_SelfPointer
    procedure :: ScalarPower => TInitialPowerLaw_ScalarPower
    procedure :: TensorPower => TInitialPowerLaw_TensorPower
    procedure :: ReadParams => TInitialPowerLaw_ReadParams
    procedure :: Effective_ns => TInitalPowerLaw_Effective_ns
    end Type TInitialPowerLaw

    Type, extends(TInitialPower) :: TSplinedInitialPower
        real(dl) :: effective_ns_for_nonlinear = -1._dl !used for halofit
        real(dl) :: kmin_scalar, kmax_scalar
        real(dl) :: kmin_tensor, kmax_tensor
        class(TSpline1D), allocatable :: Pscalar, Ptensor
    contains
    procedure :: SetScalarTable => TSplinedInitialPower_SetScalarTable
    procedure :: SetTensorTable => TSplinedInitialPower_SetTensorTable
    procedure :: SetScalarLogRegular => TSplinedInitialPower_SetScalarLogRegular
    procedure :: SetTensorLogRegular => TSplinedInitialPower_SetTensorLogRegular
    procedure :: ScalarPower => TSplinedInitialPower_ScalarPower
    procedure :: TensorPower => TSplinedInitialPower_TensorPower
    procedure :: HasTensors => TSplinedInitialPower_HasTensors
    procedure :: Effective_ns => TSplinedInitialPower_Effective_ns
    procedure, nopass :: PythonClass => TSplinedInitialPower_PythonClass
    procedure, nopass :: SelfPointer => TSplinedInitialPower_SelfPointer
    end Type TSplinedInitialPower


    public TInitialPowerLaw, TSplinedInitialPower
    contains

    function TInitialPowerLaw_PythonClass()
    character(LEN=:), allocatable :: TInitialPowerLaw_PythonClass
    TInitialPowerLaw_PythonClass = 'InitialPowerLaw'
    end function TInitialPowerLaw_PythonClass

    subroutine TInitialPowerLaw_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TInitialPowerLaw), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TInitialPowerLaw_SelfPointer

    subroutine TInitialPowerLaw_Init(this, Params)
    use classes
    use results
    use constants, only : c
    class(TInitialPowerLaw) :: this
    class(TCAMBParameters), intent(in) :: Params

    select type(Params)
    class is (CAMBParams)
        !Curvature parameter if non-flat
        this%curv = -Params%Omk/((c/1000)/Params%H0)**2
    end select

    end subroutine TInitialPowerLaw_Init

    function TInitialPowerLaw_ScalarPower(this, k)
    class(TInitialPowerLaw) :: this
    real(dl), intent(in) :: k
    real(dl) TInitialPowerLaw_ScalarPower
    real(dl) lnrat
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + curv)/|curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.


    lnrat = log(k/this%pivot_scalar)
    TInitialPowerLaw_ScalarPower = this%As * exp(lnrat * (this%ns - 1 + &
        &             lnrat * (this%nrun / 2 + this%nrunrun / 6 * lnrat)))

    end function TInitialPowerLaw_ScalarPower


    function TInitialPowerLaw_TensorPower(this,k)
    use constants
    class(TInitialPowerLaw) :: this
    !TensorPower= const for scale invariant spectrum
    !The normalization is defined so that
    ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
    !for a closed model
    ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
    !for an open model
    !Here nu^2 = (k^2 + 3*curv)/|curv|
    real(dl), intent(in) :: k
    real(dl) TInitialPowerLaw_TensorPower
    real(dl), parameter :: PiByTwo=const_pi/2._dl
    real(dl) lnrat, k_dep

    lnrat = log(k/this%pivot_tensor)
    k_dep = exp(lnrat*(this%nt + this%ntrun/2*lnrat))
    if (this%tensor_parameterization==tensor_param_indeptilt) then
        TInitialPowerLaw_TensorPower = this%r*this%As*k_dep
    else if (this%tensor_parameterization==tensor_param_rpivot) then
        TInitialPowerLaw_TensorPower = this%r*this%ScalarPower(this%pivot_tensor) * k_dep
    else if (this%tensor_parameterization==tensor_param_At) then
        TInitialPowerLaw_TensorPower = this%At * k_dep
    end if
    if (this%curv < 0) TInitialPowerLaw_TensorPower= &
        TInitialPowerLaw_TensorPower*tanh(PiByTwo*sqrt(-k**2/this%curv-3))
    end function TInitialPowerLaw_TensorPower

    function CompatKey(Ini, name)
    class(TIniFile), intent(in) :: Ini
    character(LEN=*), intent(in) :: name
    character(LEN=:), allocatable :: CompatKey
    !Allow backwards compatibility with old .ini files where initial power parameters were arrays

    if (Ini%HasKey(name//'(1)')) then
        CompatKey = name//'(1)'
        if (Ini%HasKey(name)) call MpiStop('Must have one of '//trim(name)//' or '//trim(name)//'(1)')
    else
        CompatKey = name
    end if
    end function CompatKey

    subroutine TInitialPowerLaw_ReadParams(this, Ini)
    use IniObjects
    class(TInitialPowerLaw) :: this
    class(TIniFile), intent(in) :: Ini
    logical :: WantTensors

    WantTensors = Ini%Read_Logical('get_tensor_cls', .false.)

    call Ini%Read('pivot_scalar', this%pivot_scalar)
    call Ini%Read('pivot_tensor', this%pivot_tensor)
    if (Ini%Read_Int('initial_power_num', 1) /= 1) call MpiStop('initial_power_num>1 no longer supported')
    if (WantTensors) then
        this%tensor_parameterization =  Ini%Read_Int('tensor_parameterization', tensor_param_indeptilt)
        if (this%tensor_parameterization < tensor_param_indeptilt .or. &
            &   this%tensor_parameterization > tensor_param_AT) &
            &   call MpiStop('InitialPower: unknown tensor_parameterization')
    end if
    this%r = 1
    this%ns = Ini%Read_Double(CompatKey(Ini,'scalar_spectral_index'))
    call Ini%Read(CompatKey(Ini,'scalar_nrun'), this%nrun)
    call Ini%Read(CompatKey(Ini,'scalar_nrunrun'), this%nrunrun)

    if (WantTensors) then
        this%nt = Ini%Read_Double(CompatKey(Ini,'tensor_spectral_index'))
        call Ini%Read(CompatKey(Ini,'tensor_nrun'),this%ntrun)
        if (this%tensor_parameterization == tensor_param_AT) then
            this%At = Ini%Read_Double(CompatKey(Ini,'tensor_amp'))
        else
            this%r = Ini%Read_Double(CompatKey(Ini,'initial_ratio'))
        end if
    else
        this%r =0
        this%At=0
    end if

    call Ini%Read(CompatKey(Ini,'scalar_amp'),this%As)
    !Always need this as may want to set tensor amplitude even if scalars not computed

    end subroutine TInitialPowerLaw_ReadParams

    function TInitalPowerLaw_Effective_ns(this)
    class(TInitialPowerLaw) :: this
    real(dl) :: TInitalPowerLaw_Effective_ns

    TInitalPowerLaw_Effective_ns = this%ns

    end function TInitalPowerLaw_Effective_ns


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

    if (k <= this%kmin_scalar) then
        TSplinedInitialPower_ScalarPower = this%Pscalar%F(1)
    elseif (k >= this%kmax_scalar) then
        TSplinedInitialPower_ScalarPower = this%Pscalar%F(this%Pscalar%n)
    else
        TSplinedInitialPower_ScalarPower = this%Pscalar%Value(k)
    end if

    end function TSplinedInitialPower_ScalarPower

    function TSplinedInitialPower_TensorPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_TensorPower

    if (k <= this%kmin_tensor) then
        TSplinedInitialPower_TensorPower = this%Ptensor%F(1)
    elseif (k >= this%kmax_tensor) then
        TSplinedInitialPower_TensorPower = this%Ptensor%F(this%Ptensor%n)
    else
        TSplinedInitialPower_TensorPower = this%Ptensor%Value(k)
    end if

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
        this%kmin_scalar = k(1)
        this%kmax_scalar = k(n)
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
        this%kmin_tensor = k(1)
        this%kmax_tensor = k(n)
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
        this%kmin_scalar = kmin
        this%kmax_scalar = kmax
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
        this%kmin_tensor = kmin
        this%kmax_tensor = kmax
    end if

    end subroutine TSplinedInitialPower_SetTensorLogRegular

    function TSplinedInitialPower_Effective_ns(this)
    use config
    class(TSplinedInitialPower) :: this
    real(dl) :: TSplinedInitialPower_Effective_ns

    if (this%effective_ns_for_nonlinear==-1._dl) then
        call GlobalError('TSplinedInitialPower: effective_ns_for_nonlinear not set',error_inital_power)
    else
        TSplinedInitialPower_Effective_ns = this%effective_ns_for_nonlinear
    end if
    end function TSplinedInitialPower_Effective_ns



    end module InitialPower
