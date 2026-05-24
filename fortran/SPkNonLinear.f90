module SPkNonLinear
! Wrapper non-linear model: apply SP(k) suppression on top of a base non-linear model.
use results
use NonLinear, only: THalofit, halofit_mead, halofit_mead2015, halofit_mead2016, halofit_mead2020_feedback
use SPkModel
use config
use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
implicit none
private

type, extends(TNonLinearModel) :: TSPkNonLinear
    class(TNonLinearModel), allocatable :: BaseModel
    logical :: SPk_feedback = .false.
    integer :: SPk_SO = 200
    integer :: SPk_relation_kind = SPk_rel_power_law
    real(dl) :: SPk_fb_a = 1.0_dl
    real(dl) :: SPk_fb_pow = 0.0_dl
    real(dl) :: SPk_fb_pivot = 1.0_dl
    real(dl) :: SPk_alpha = 0.0_dl
    real(dl) :: SPk_beta = 0.0_dl
    real(dl) :: SPk_gamma = 0.0_dl
    real(dl) :: SPk_epsilon = 0.0_dl
    real(dl) :: SPk_m_pivot = 1.0_dl
contains
    procedure :: Init => TSPkNonLinear_Init
    procedure :: ReadParams => TSPkNonLinear_ReadParams
    procedure :: GetNonLinRatios => TSPkNonLinear_GetNonLinRatios
    procedure, nopass :: PythonClass => TSPkNonLinear_PythonClass
    procedure, nopass :: SelfPointer => TSPkNonLinear_SelfPointer
    procedure, private :: ValidateConfig => TSPkNonLinear_ValidateConfig
end type TSPkNonLinear

public TSPkNonLinear

contains

function TSPkNonLinear_PythonClass()
character(LEN=:), allocatable :: TSPkNonLinear_PythonClass

TSPkNonLinear_PythonClass = 'SPkNonLinear'

end function TSPkNonLinear_PythonClass

subroutine TSPkNonLinear_SelfPointer(cptr, P)
use iso_c_binding
Type(c_ptr) :: cptr
Type(TSPkNonLinear), pointer :: PType
class(TPythonInterfacedClass), pointer :: P

call c_f_pointer(cptr, PType)
P => PType

end subroutine TSPkNonLinear_SelfPointer

subroutine TSPkNonLinear_Init(this, State)
class(TSPkNonLinear) :: this
class(TCAMBdata), target :: State

if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)
this%Min_kh_nonlinear = this%BaseModel%Min_kh_nonlinear
call this%BaseModel%Init(State)

end subroutine TSPkNonLinear_Init

subroutine TSPkNonLinear_ReadParams(this, Ini)
use IniObjects
class(TSPkNonLinear) :: this
class(TIniFile), intent(in) :: Ini

if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)
call this%BaseModel%ReadParams(Ini)

this%SPk_feedback = Ini%Read_Logical('SPk_feedback', .false.)
this%SPk_SO = Ini%Read_Int('SPk_SO', 200)
this%SPk_relation_kind = Ini%Read_Int('SPk_relation_kind', SPk_rel_power_law)
this%SPk_fb_a = Ini%Read_Double('SPk_fb_a', 1.0_dl)
this%SPk_fb_pow = Ini%Read_Double('SPk_fb_pow', 0.0_dl)
this%SPk_fb_pivot = Ini%Read_Double('SPk_fb_pivot', 1.0_dl)
this%SPk_alpha = Ini%Read_Double('SPk_alpha', 0.0_dl)
this%SPk_beta = Ini%Read_Double('SPk_beta', 0.0_dl)
this%SPk_gamma = Ini%Read_Double('SPk_gamma', 0.0_dl)
this%SPk_epsilon = Ini%Read_Double('SPk_epsilon', 0.0_dl)
this%SPk_m_pivot = Ini%Read_Double('SPk_m_pivot', 1.0_dl)

end subroutine TSPkNonLinear_ReadParams

subroutine TSPkNonLinear_ValidateConfig(this)
class(TSPkNonLinear), intent(in) :: this

if (this%SPk_SO /= 200 .and. this%SPk_SO /= 500) then
    call MpiStop('SP(k): SPk_SO must be 200 or 500')
end if
if (this%SPk_relation_kind < SPk_rel_power_law .or. this%SPk_relation_kind > SPk_rel_double_power_law) then
    call MpiStop('SP(k): SPk_relation_kind must be 1 (power_law), 2 (cosmo_power_law), or 3 (double_power_law)')
end if

if (this%SPk_relation_kind == SPk_rel_power_law) then
    if (this%SPk_fb_pivot <= 0.0_dl) call MpiStop('SP(k): SPk_fb_pivot must be > 0 for power_law relation')
end if
if (this%SPk_relation_kind == SPk_rel_double_power_law) then
    if (this%SPk_m_pivot <= 0.0_dl) call MpiStop('SP(k): SPk_m_pivot must be > 0 for double_power_law relation')
end if

select type (base => this%BaseModel)
type is (THalofit)
    if (this%SPk_feedback .and. base%halofit_version == halofit_mead2020_feedback) then
        call MpiStop('SP(k) is not compatible with halofit_version=mead2020_feedback')
    end if
    if (this%SPk_feedback .and. &
        (base%halofit_version == halofit_mead .or. base%halofit_version == halofit_mead2015 .or. base%halofit_version == halofit_mead2016) .and. &
        (abs(base%HMcode_A_baryon - 3.13_dl) > 1e-12_dl .or. abs(base%HMcode_eta_baryon - 0.603_dl) > 1e-12_dl)) then
        call MpiStop('SP(k) cannot be combined with HMCode_A_baryon/HMCode_eta_baryon baryonic corrections in HMCode 2015/2016')
    end if
class default
    continue
end select

end subroutine TSPkNonLinear_ValidateConfig

subroutine TSPkNonLinear_GetNonLinRatios(this, State, CAMB_Pk)
class(TSPkNonLinear) :: this
class(TCAMBdata) :: State
type(MatterPowerData), target :: CAMB_Pk
integer :: itf, i
real(dl) :: rk, rk_eval, spk_sup, spk_href, spk_eratio
real(dl) :: spk_fb, spk_m_opt, spk_fb_min, spk_fb_max
logical, save :: warned_spk_z_outside = .false.
logical, save :: warned_spk_k_clamped = .false.
logical, save :: warned_spk_fb_outside = .false.

if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)

call this%ValidateConfig()
call this%BaseModel%GetNonLinRatios(State, CAMB_Pk)

if (.not. this%SPk_feedback) return

select type (State)
class is (CAMBdata)
    if (this%SPk_relation_kind == SPk_rel_cosmo_power_law .or. this%SPk_relation_kind == SPk_rel_double_power_law) then
        spk_href = State%Hofz(0.3_dl)
    else
        spk_href = 1.0_dl
    end if
class default
    call MpiStop('SP(k): unsupported state type for Hofz evaluation')
end select

do itf = 1, CAMB_Pk%num_z
    if (CAMB_Pk%redshifts(itf) < SPk_calibrated_z_min .or. CAMB_Pk%redshifts(itf) > SPk_calibrated_z_max) then
        if (FeedbackLevel > 0 .and. .not. warned_spk_z_outside) then
            write(*,'(A,F8.3,A,F6.2,A,F6.2,A)') 'WARNING: SP(k) skipped outside calibrated redshift range. z=', &
                CAMB_Pk%redshifts(itf), ', calibrated range=[', SPk_calibrated_z_min, ',', SPk_calibrated_z_max, '].'
            warned_spk_z_outside = .true.
        end if
        cycle
    end if
    select type (State)
    class is (CAMBdata)
        if (this%SPk_relation_kind == SPk_rel_cosmo_power_law .or. this%SPk_relation_kind == SPk_rel_double_power_law) then
            spk_eratio = State%Hofz(CAMB_Pk%redshifts(itf)) / spk_href
        else
            spk_eratio = 1.0_dl
        end if
    class default
        call MpiStop('SP(k): unsupported state type for Hofz evaluation')
    end select
    do i = 1, CAMB_Pk%num_k
        rk = exp(CAMB_Pk%log_kh(i))
        if (rk < SPk_calibrated_k_min) cycle
        if (rk > SPk_calibrated_k_max) then
            if (FeedbackLevel > 0 .and. .not. warned_spk_k_clamped) then
                write(*,'(A,F8.3,A,F6.2,A)') 'WARNING: SP(k) input k exceeds calibrated range; clamping to k_max=', &
                    rk, ' -> ', SPk_calibrated_k_max, ' h/Mpc.'
                warned_spk_k_clamped = .true.
            end if
        end if
        rk_eval = min(rk, SPk_calibrated_k_max)
        call SPk_ComputeFb(this%SPk_SO, rk_eval, CAMB_Pk%redshifts(itf), this%SPk_relation_kind, &
            this%SPk_fb_a, this%SPk_fb_pow, this%SPk_fb_pivot, this%SPk_alpha, this%SPk_beta, this%SPk_gamma, &
            this%SPk_epsilon, this%SPk_m_pivot, spk_eratio, spk_fb, spk_m_opt)
        call SPk_GetFbLimits(this%SPk_SO, CAMB_Pk%redshifts(itf), spk_m_opt, spk_fb_min, spk_fb_max)
        if (spk_fb < spk_fb_min .or. spk_fb > spk_fb_max) then
            if (FeedbackLevel > 0 .and. .not. warned_spk_fb_outside) then
                write(*,'(A)') 'WARNING: SP(k) baryon fraction outside calibrated fitting limits; setting out-of-range points to NaN.'
                warned_spk_fb_outside = .true.
            end if
            CAMB_Pk%nonlin_ratio(i, itf) = ieee_value(1.0_dl, ieee_quiet_nan)
            cycle
        end if
        spk_sup = SPk_Suppression(this%SPk_SO, rk_eval, CAMB_Pk%redshifts(itf), &
            this%SPk_relation_kind, this%SPk_fb_a, this%SPk_fb_pow, this%SPk_fb_pivot, &
            this%SPk_alpha, this%SPk_beta, this%SPk_gamma, &
            this%SPk_epsilon, this%SPk_m_pivot, spk_eratio)
        CAMB_Pk%nonlin_ratio(i, itf) = CAMB_Pk%nonlin_ratio(i, itf) * sqrt(spk_sup)
    end do
end do

end subroutine TSPkNonLinear_GetNonLinRatios

end module SPkNonLinear
