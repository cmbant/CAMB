    module RungeKuttaDP45Module
    use Precision
    implicit none
    private

    integer, parameter, public :: RK45ErrorMixed = 0
    integer, parameter, public :: RK45ErrorAbsolute = 1
    integer, parameter, public :: RK45ErrorRelative = 2
    integer, parameter, public :: RK45ErrorRelativeFloor = 3

    type, public :: RungeKuttaDP45Settings
        ! User-configurable controls. Defaults mainly preserves CAMB's historical usage
        ! unless the caller overrides fields and enters with ind = 2. The
        ! reuse_first_derivative request is also preserved across an ind = 1 reset
        ! so callers can enable FSAL without bypassing the normal state reset.
        ! error_control chooses the max-norm weighting:
        !   RK45ErrorMixed         1/max(1, abs(y(k))) (default).
        !   RK45ErrorAbsolute      1.
        !   RK45ErrorRelative      1/abs(y(k)).
        !   RK45ErrorRelativeFloor 1/max(error_floor, abs(y(k))).
        integer :: error_control = RK45ErrorMixed
        ! Scalar floor for RK45ErrorRelativeFloor. Default 0 is unused unless that
        ! mode is selected.
        real(dl) :: error_floor = 0._dl
        ! Explicit minimum trial-step magnitude. Default 0 uses the historical
        ! automatic estimate 10*max(tiny, eps*max(weighted_norm_y/tol, abs(x))).
        real(dl) :: min_step_size = 0._dl
        ! Explicit first-step magnitude. Default 0 uses computed_max_step_size * tol**(1/5).
        real(dl) :: initial_step_size = 0._dl
        ! Problem scale entering the acceptance test and optional hmax cap.
        ! Default 0 promotes to 1; if non-zero it also limits hmax to 2/abs(problem_scale).
        real(dl) :: problem_scale = 0._dl
        ! Explicit maximum trial-step magnitude. Default 0 uses CAMB default of 20 (used to be 2) unless
        ! problem_scale imposes a tighter 2/abs(scale) limit.
        real(dl) :: max_step_size = 0._dl
        ! Optional hard cap on derivative evaluations. Default 0 disables this check.
        integer :: max_function_evaluations = 0
        ! Return with ind = 4 after selecting the next step size but before the trial step.
        logical :: interrupt_before_trial_step = .false.
        ! Return with ind = 5 or 6 after estimating trial-step error and before
        ! applying the accept/reject update.
        logical :: interrupt_after_trial_step = .false.
        ! Reuse the first-same-as-last (FSAL) stage across calls: on reaching xend
        ! save the last stage-7 derivative (= fcn at the returned x, y) so the next
        ! ind = 3 continuation skips the initial derivative evaluation. Only valid
        ! if the caller does not change x, y, n or the equations between calls;
        ! callers that do must reset ind or set first_derivative_valid = .false.
        logical :: reuse_first_derivative = .false.
        ! Allow the step to stretch by up to this factor beyond the natural step
        ! size to reach xend in one step instead of two clipped ones. The stretched
        ! step must still pass the error test and respect the maximum step size.
        ! Values <= 1 recover the historical clip to half the remaining interval.
        real(dl) :: xend_stretch = 1.2_dl
        ! Persistent solver state below this point, maintained across re-entry.
        ! These correspond to the historical c(12:24) bookkeeping values.
        real(dl) :: weighted_norm_y = 0._dl
        real(dl) :: computed_min_step_size = 0._dl
        real(dl) :: current_step_size = 0._dl
        real(dl) :: computed_scale = 0._dl
        real(dl) :: computed_max_step_size = 0._dl
        real(dl) :: trial_x = 0._dl
        real(dl) :: trial_step = 0._dl
        real(dl) :: estimated_error = 0._dl
        real(dl) :: previous_xend = 0._dl
        logical :: xend_reached = .false.
        logical :: first_derivative_valid = .false.
        integer :: successful_steps = 0
        integer :: successive_failures = 0
        integer :: function_evaluations = 0
    end type RungeKuttaDP45Settings

    public :: TClassRungeKuttaDP45
    abstract interface
    subroutine TClassRungeKuttaDP45(this, n, fcn, x, y, xend, tol, ind, settings, nw, w)
    use Precision
    use classes, only : TCambComponent
    import :: RungeKuttaDP45Settings
    class(TCambComponent), target :: this
    integer, intent(in) :: n, nw
    integer, intent(inout) :: ind
    real(dl), intent(inout) :: x, y(n), w(nw, 9)
    real(dl), intent(in) :: xend, tol
    type(RungeKuttaDP45Settings), intent(inout) :: settings
    external fcn
    end subroutine TClassRungeKuttaDP45
    end interface

    end module RungeKuttaDP45Module


    subroutine RungeKuttaDP45(EV, n, fcn, x, y, xend, tol_in, ind, settings, nw, w)
    use Precision
    use Config, only : GlobalError, error_evolution
    use RungeKuttaDP45Module, only : RK45ErrorMixed, RK45ErrorAbsolute, RK45ErrorRelative, &
        RK45ErrorRelativeFloor, RungeKuttaDP45Settings
    implicit none
    integer, intent(in) :: n, nw
    integer, intent(inout) :: ind
    real(dl), intent(inout) :: x, y(n), w(nw, 9)
    real(dl), intent(in) :: xend, tol_in
    type(RungeKuttaDP45Settings), intent(inout) :: settings
    real(dl) :: tol, temp
    real(dl), parameter :: one_fifth = 1._dl / 5._dl
    real(dl), parameter :: default_max_step_size = 20._dl
    real(dl), parameter :: dp_a21 = 1._dl / 5._dl
    real(dl), parameter :: dp_a31 = 3._dl / 40._dl
    real(dl), parameter :: dp_a32 = 9._dl / 40._dl
    real(dl), parameter :: dp_a41 = 44._dl / 45._dl
    real(dl), parameter :: dp_a42 = -56._dl / 15._dl
    real(dl), parameter :: dp_a43 = 32._dl / 9._dl
    real(dl), parameter :: dp_a51 = 19372._dl / 6561._dl
    real(dl), parameter :: dp_a52 = -25360._dl / 2187._dl
    real(dl), parameter :: dp_a53 = 64448._dl / 6561._dl
    real(dl), parameter :: dp_a54 = -212._dl / 729._dl
    real(dl), parameter :: dp_a61 = 9017._dl / 3168._dl
    real(dl), parameter :: dp_a62 = -355._dl / 33._dl
    real(dl), parameter :: dp_a63 = 46732._dl / 5247._dl
    real(dl), parameter :: dp_a64 = 49._dl / 176._dl
    real(dl), parameter :: dp_a65 = -5103._dl / 18656._dl
    real(dl), parameter :: dp_b1 = 35._dl / 384._dl
    real(dl), parameter :: dp_b3 = 500._dl / 1113._dl
    real(dl), parameter :: dp_b4 = 125._dl / 192._dl
    real(dl), parameter :: dp_b5 = -2187._dl / 6784._dl
    real(dl), parameter :: dp_b6 = 11._dl / 84._dl
    real(dl), parameter :: dp_e1 = 71._dl / 57600._dl
    real(dl), parameter :: dp_e3 = -71._dl / 16695._dl
    real(dl), parameter :: dp_e4 = 71._dl / 1920._dl
    real(dl), parameter :: dp_e5 = -17253._dl / 339200._dl
    real(dl), parameter :: dp_e6 = 22._dl / 525._dl
    real(dl), parameter :: dp_e7 = -1._dl / 40._dl
    real(dl), parameter :: machine_roundoff = epsilon(1._dl)
    real(dl), parameter :: machine_tiny = tiny(1._dl)
    logical :: resume_after_interrupt1, resume_after_interrupt2
    real :: EV

! RungeKuttaDP45 preserves CAMB's historical control flow, interrupt
! semantics, and workspace layout while exposing the old c(*) control array as
! named fields on a settings/state object.
!
! Arguments:
!   EV, fcn: opaque object and derivative callback passed through unchanged.
!   n, x, y, xend, tol_in: standard ODE system size, state, target x, and
!       caller-provided tolerance. tol_in is used directly.
!   ind: solver entry/return status.
!       1 reset with default settings, preserving reuse_first_derivative,
!       2 reset using caller-populated settings,
!       3 continue normal re-entry, 4/5/6 resume after interrupts,
!       returns 3/4/5/6 or -1/-2/-3.
!   settings:
!       Configurable fields are error_control, error_floor, min_step_size,
!       initial_step_size, problem_scale, max_step_size,
!       max_function_evaluations, interrupt_before_trial_step,
!       interrupt_after_trial_step, reuse_first_derivative, and xend_stretch.
!       The remaining fields store persistent solver state between calls.
!   nw, w: work array dimensions retained for compatibility with existing
!       CAMB callers.
!
! Error control modes:
!   RK45ErrorMixed         default 1 / max(1, abs(y(k))) weights.
!   RK45ErrorAbsolute      absolute error control.
!   RK45ErrorRelative      relative error control.
!   RK45ErrorRelativeFloor relative control with scalar floor error_floor.

    external fcn
    tol = tol_in

    resume_after_interrupt1 = .false.
    resume_after_interrupt2 = .false.

    select case (ind)
    case (1, 2)
        if (n > nw .or. tol <= 0._dl) then
            call abort_runge_kutta_dp45()
            return
        end if

        if (ind == 1) then
            settings = RungeKuttaDP45Settings(reuse_first_derivative=settings%reuse_first_derivative)
        else
            call normalize_settings(settings)
        end if

        settings%previous_xend = x
        settings%xend_reached = .false.
        settings%first_derivative_valid = .false.
        settings%successful_steps = 0
        settings%successive_failures = 0
        settings%function_evaluations = 0

    case (3)
        if (settings%xend_reached .and. (x /= settings%previous_xend .or. xend == settings%previous_xend)) then
            call abort_runge_kutta_dp45()
            return
        end if
        settings%xend_reached = .false.

    case (4)
        resume_after_interrupt1 = .true.

    case (5, 6)
        resume_after_interrupt2 = .true.

    case default
        call abort_runge_kutta_dp45()
        return
    end select

    step_loop: do
        if (.not. resume_after_interrupt1 .and. .not. resume_after_interrupt2) then
            if (settings%max_function_evaluations /= 0 .and. &
                settings%function_evaluations >= settings%max_function_evaluations) then
                ind = -1
                return
            end if

            if (ind /= 6) then
                if (settings%reuse_first_derivative .and. settings%first_derivative_valid) then
                    ! w(:,1) still holds fcn(x, y) saved from the final stage of the
                    ! step that reached the previous xend.
                    settings%first_derivative_valid = .false.
                else
                    call fcn(EV, n, x, y, w(1, 1))
                    settings%function_evaluations = settings%function_evaluations + 1
                end if
            end if

            settings%computed_min_step_size = abs(settings%min_step_size)
            if (settings%computed_min_step_size == 0._dl) then
                select case (settings%error_control)
                case (RK45ErrorAbsolute)
                    settings%weighted_norm_y = maxval(abs(y(1:n)))
                case (RK45ErrorRelative)
                    settings%weighted_norm_y = 1._dl
                case (RK45ErrorRelativeFloor)
                    temp = maxval(abs(y(1:n)) / settings%error_floor)
                    settings%weighted_norm_y = min(temp, 1._dl)
                case default
                    temp = maxval(abs(y(1:n)))
                    settings%weighted_norm_y = min(temp, 1._dl)
                end select

                settings%computed_min_step_size = 10._dl * max( &
                    machine_tiny, machine_roundoff * max(settings%weighted_norm_y / tol, abs(x)) &
                    )
            end if

            settings%computed_scale = abs(settings%problem_scale)
            if (settings%computed_scale == 0._dl) settings%computed_scale = 1._dl

            if (settings%max_step_size /= 0._dl .and. settings%problem_scale /= 0._dl) then
                settings%computed_max_step_size = min(abs(settings%max_step_size), 2._dl / abs(settings%problem_scale))
            else if (settings%max_step_size /= 0._dl) then
                settings%computed_max_step_size = abs(settings%max_step_size)
            else if (settings%problem_scale /= 0._dl) then
                settings%computed_max_step_size = 2._dl / abs(settings%problem_scale)
            else
                settings%computed_max_step_size = default_max_step_size
            end if

            if (settings%computed_min_step_size > settings%computed_max_step_size) then
                ind = -2
                return
            end if

            if (ind <= 2) then
                settings%current_step_size = abs(settings%initial_step_size)
                if (settings%current_step_size == 0._dl) then
                    settings%current_step_size = settings%computed_max_step_size * tol**one_fifth
                end if
            else if (settings%successive_failures <= 1) then
                temp = 2._dl * settings%current_step_size
                if (tol < (2._dl / 0.9_dl)**5 * settings%estimated_error) then
                    temp = 0.9_dl * (tol / settings%estimated_error)**one_fifth * settings%current_step_size
                end if
                settings%current_step_size = max(temp, 0.5_dl * settings%current_step_size)
            else
                settings%current_step_size = 0.5_dl * settings%current_step_size
            end if

            settings%current_step_size = min(settings%current_step_size, settings%computed_max_step_size)
            settings%current_step_size = max(settings%current_step_size, settings%computed_min_step_size)

            if (settings%interrupt_before_trial_step) then
                ind = 4
                return
            end if
        end if

        if (.not. resume_after_interrupt2) then
            temp = abs(xend - x)
            if (settings%current_step_size < temp .and. &
                (temp > settings%xend_stretch * settings%current_step_size .or. &
                temp > settings%computed_max_step_size)) then
                settings%current_step_size = min(settings%current_step_size, 0.5_dl * temp)
                settings%trial_x = x + sign(settings%current_step_size, xend - x)
            else
                settings%current_step_size = temp
                settings%trial_x = xend
            end if

            settings%trial_step = settings%trial_x - x

            w(1:n, 9) = y(1:n) + settings%trial_step * dp_a21 * w(1:n, 1)
            call fcn(EV, n, x + settings%trial_step / 5._dl, w(1, 9), w(1, 2))

            w(1:n, 9) = y(1:n) + settings%trial_step * (dp_a31 * w(1:n, 1) + dp_a32 * w(1:n, 2))
            call fcn(EV, n, x + settings%trial_step * (3._dl / 10._dl), w(1, 9), w(1, 3))

            w(1:n, 9) = y(1:n) + settings%trial_step * (dp_a41 * w(1:n, 1) + dp_a42 * w(1:n, 2) + &
                dp_a43 * w(1:n, 3))
            call fcn(EV, n, x + settings%trial_step * (4._dl / 5._dl), w(1, 9), w(1, 4))

            w(1:n, 9) = y(1:n) + settings%trial_step * (dp_a51 * w(1:n, 1) + dp_a52 * w(1:n, 2) + &
                dp_a53 * w(1:n, 3) + dp_a54 * w(1:n, 4))
            call fcn(EV, n, x + settings%trial_step * (8._dl / 9._dl), w(1, 9), w(1, 5))

            w(1:n, 9) = y(1:n) + settings%trial_step * (dp_a61 * w(1:n, 1) + dp_a62 * w(1:n, 2) + &
                dp_a63 * w(1:n, 3) + dp_a64 * w(1:n, 4) + dp_a65 * w(1:n, 5))
            call fcn(EV, n, x + settings%trial_step, w(1, 9), w(1, 6))

            w(1:n, 9) = y(1:n) + settings%trial_step * (dp_b1 * w(1:n, 1) + dp_b3 * w(1:n, 3) + &
                dp_b4 * w(1:n, 4) + dp_b5 * w(1:n, 5) + dp_b6 * w(1:n, 6))
            call fcn(EV, n, x + settings%trial_step, w(1, 9), w(1, 7))

            settings%function_evaluations = settings%function_evaluations + 6

            w(1:n, 2) = dp_e1 * w(1:n, 1) + dp_e3 * w(1:n, 3) + dp_e4 * w(1:n, 4) + &
                dp_e5 * w(1:n, 5) + dp_e6 * w(1:n, 6) + dp_e7 * w(1:n, 7)

            select case (settings%error_control)
            case (RK45ErrorAbsolute)
                temp = maxval(abs(w(1:n, 2)))
            case (RK45ErrorRelative)
                temp = maxval(abs(w(1:n, 2) / y(1:n)))
            case (RK45ErrorRelativeFloor)
                temp = maxval(abs(w(1:n, 2)) / max(settings%error_floor, abs(y(1:n))))
            case default
                temp = maxval(abs(w(1:n, 2)) / max(1._dl, abs(y(1:n))))
            end select

            settings%estimated_error = temp * settings%current_step_size * settings%computed_scale

            ind = 5
            if (settings%estimated_error > tol) ind = 6

            if (settings%interrupt_after_trial_step) return
        end if

        if (ind == 5) then
            x = settings%trial_x
            y(1:n) = w(1:n, 9)
            settings%successful_steps = settings%successful_steps + 1
            settings%successive_failures = 0

            if (x == xend) then
                ind = 3
                settings%previous_xend = xend
                settings%xend_reached = .true.
                if (settings%reuse_first_derivative) then
                    ! FSAL: stage 7 was evaluated at (xend, y(xend)), so it can seed
                    ! the first stage of the next continuation call.
                    w(1:n, 1) = w(1:n, 7)
                    settings%first_derivative_valid = .true.
                end if
                return
            end if

            w(1:n, 1) = w(1:n, 7)
            ind = 6
        else
            settings%successive_failures = settings%successive_failures + 1

            if (settings%current_step_size <= settings%computed_min_step_size) then
                ind = -3
                return
            end if
        end if

        resume_after_interrupt1 = .false.
        resume_after_interrupt2 = .false.
    end do step_loop

    contains

    subroutine normalize_settings(settings)
    type(RungeKuttaDP45Settings), intent(inout) :: settings

    settings%error_control = abs(settings%error_control)
    settings%error_floor = abs(settings%error_floor)
    settings%min_step_size = abs(settings%min_step_size)
    settings%initial_step_size = abs(settings%initial_step_size)
    settings%problem_scale = abs(settings%problem_scale)
    settings%max_step_size = abs(settings%max_step_size)
    settings%max_function_evaluations = abs(settings%max_function_evaluations)
    settings%xend_stretch = abs(settings%xend_stretch)

    end subroutine normalize_settings

    subroutine abort_runge_kutta_dp45()
    write (*,*) 'Error in RungeKuttaDP45, x =', x, ' xend =', xend
    call GlobalError('RungeKuttaDP45 error', error_evolution)
    end subroutine abort_runge_kutta_dp45

    end subroutine RungeKuttaDP45
