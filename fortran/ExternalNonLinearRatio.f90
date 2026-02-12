
    module ExternalNonLinearRatio
    ! Non-linear model that applies a user-supplied ratio sqrt(P_NL/P_L)
    ! from an external source (e.g. CCL, axionHMcode).
    ! The ratio is stored as a 2D interpolation grid over (k/h, z).
    use results
    use transfer
    use Interpolation, only : TInterpGrid2D
    implicit none
    private

    type, extends(TNonLinearModel) :: TExternalNonLinearRatio
        Type(TInterpGrid2D) :: Ratio
        logical :: ratio_set = .false.
    contains
    procedure :: GetNonLinRatios => TExternalNonLinearRatio_GetNonLinRatios
    procedure :: GetNonLinRatios_All => TExternalNonLinearRatio_GetNonLinRatios_All
    procedure :: SetRatio => TExternalNonLinearRatio_SetRatio
    procedure :: ClearRatio => TExternalNonLinearRatio_ClearRatio
    procedure, nopass :: SelfPointer => TExternalNonLinearRatio_SelfPointer
    end type TExternalNonLinearRatio

    public TExternalNonLinearRatio
    contains

    subroutine TExternalNonLinearRatio_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TExternalNonLinearRatio), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TExternalNonLinearRatio_SelfPointer

    subroutine TExternalNonLinearRatio_SetRatio(this, nk, nz, k_h_arr, z_arr, ratio_arr)
    ! Set the non-linear ratio grid sqrt(P_NL/P_L)
    class(TExternalNonLinearRatio) :: this
    integer, intent(in) :: nk, nz
    real(dl), intent(in) :: k_h_arr(nk)
    real(dl), intent(in) :: z_arr(nz)
    real(dl), intent(in) :: ratio_arr(nk, nz)

    call this%Ratio%Init(k_h_arr, z_arr, ratio_arr)
    this%ratio_set = .true.

    end subroutine TExternalNonLinearRatio_SetRatio

    subroutine TExternalNonLinearRatio_ClearRatio(this)
    ! Clear the stored ratio, releasing interpolation data
    class(TExternalNonLinearRatio) :: this

    call this%Ratio%Clear()
    this%ratio_set = .false.

    end subroutine TExternalNonLinearRatio_ClearRatio

    subroutine TExternalNonLinearRatio_GetNonLinRatios(this, State, CAMB_Pk)
    ! Fill CAMB_Pk%nonlin_ratio from the stored interpolation grid.
    ! Clamps k and z to the grid range to prevent unreliable polynomial extrapolation.
    class(TExternalNonLinearRatio) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    integer :: ik, iz
    real(dl) :: kh, z, kh_clamped, z_clamped
    real(dl) :: k_min, k_max, z_min, z_max

    if (.not. this%ratio_set) &
        error stop 'ExternalNonLinearRatio: ratio not set. Call SetRatio first.'

    ! Get bounds of the ratio grid to clamp values
    k_min = this%Ratio%x(1)
    k_max = this%Ratio%x(this%Ratio%nx)
    z_min = this%Ratio%y(1)
    z_max = this%Ratio%y(this%Ratio%ny)

    CAMB_Pk%nonlin_ratio = 1.0_dl
    do iz = 1, CAMB_Pk%num_z
        z = CAMB_Pk%redshifts(iz)
        z_clamped = max(z_min, min(z_max, z))
        do ik = 1, CAMB_Pk%num_k
            kh = exp(CAMB_Pk%log_kh(ik))
            kh_clamped = max(k_min, min(k_max, kh))
            CAMB_Pk%nonlin_ratio(ik, iz) = this%Ratio%Value(kh_clamped, z_clamped)
        end do
    end do

    end subroutine TExternalNonLinearRatio_GetNonLinRatios

    subroutine TExternalNonLinearRatio_GetNonLinRatios_All(this, State, CAMB_Pk)
    class(TExternalNonLinearRatio) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk

    error stop 'ExternalNonLinearRatio: GetNonLinRatios_All not supported (no velocity corrections)'

    end subroutine TExternalNonLinearRatio_GetNonLinRatios_All

    end module ExternalNonLinearRatio
