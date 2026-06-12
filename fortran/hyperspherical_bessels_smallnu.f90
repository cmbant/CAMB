! Approximate open-universe hyperspherical Bessel functions for nu << L.
!
! Implements
!
!   u_l^nu(chi)  ~=  sqrt(sinh(pi nu)/(pi nu))
!                     * K_{i nu}( 2 Lambda_{l nu} exp(-chi) )
!
!   phi_l^nu(chi) = u_l^nu(chi)/sinh(chi)
!
! with
!
!   log Lambda_{l nu} = Im log Gamma(l+1+i nu)/nu,
!
! using a branch-continuous implementation.  The code evaluates the scaled
! Macdonald function directly, avoiding overflow in the large exp(pi nu/2)
! and small K_{i nu} factors.
!
! This is intended as a fast special branch for the open K=-1 case when
! nu/L is small.  This version is tuned for O(1e-4) peak-normalized use,
! not full double-precision special-function evaluation. The scaled
! K_{i nu} router now uses:
!
!   * the I_{i nu} power series while q=x^2/(4 nu) is modest;
!   * a leading uniform Airy turning-point approximation for nu >= 20
!     outside the safe series region;
!   * the old real-axis integral only for nu < 20.
!
! O(1e-4) validity:
!   l > 20 and nu <= max(1.2e-4 l^1.5, min(12*(l/1000)^3, 0.032*l)).
! For l >= 3000, a leading Liouville amplitude correction broadens the
! validated high-l gate to max(old_gate, 0.04*l, min(8*(l/1000)^2, 0.16*l)).
! Assume inputs all pre-validated

    module HypersphericalBesselSmallNu
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use MathUtils, only: airy_ai_fast
    implicit none
    private

    integer, parameter :: dp = real64

    real(dp), parameter :: pi_dp      = 3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: twopi_dp   = 6.283185307179586476925286766559005768394_dp
    real(dp), parameter :: euler_gamma = 0.577215664901532860606512090082402431043_dp
    real(dp), parameter :: log_two = 0.6931471805599453094172321214581765680755_dp
    real(dp), parameter :: log_tiny_safe = -700.0_dp
    real(dp), parameter :: nu_zero = 1.0e-8_dp

    ! Fast/robust routing constants for the scaled K_{i nu} evaluator.
    ! The power series is cheap and very accurate while q=x^2/(4 nu)
    ! stays modest.  Outside that range, for nu >= nu_airy_min, the
    ! leading uniform Airy approximation avoids the catastrophic cancellation
    ! of the I-series near and beyond the turning point x ~= nu.
    real(dp), parameter :: nu_airy_min = 20.0_dp
    real(dp), parameter :: series_q_max  = 25.0_dp
    real(dp), parameter :: series_term2_tol = 4.0e-32_dp
    integer, parameter :: smallnu_l_min = 20
    integer, parameter :: smallnu_amp_l_min = 3000

    public :: open_smallnu_loglambda
    public :: open_smallnu_action_u
    public :: open_smallnu_action_u_from_loglambda
    public :: open_smallnu_ok
    public :: open_smallnu_u

    contains

    function open_smallnu_action_u(l, nu, chi) result(u)
    integer,  intent(in) :: l
    real(dp), intent(in) :: nu, chi
    real(dp) :: u, loglambda

    loglambda = open_smallnu_loglambda(l, nu)
    u = open_smallnu_action_u_from_loglambda(l, nu, chi, loglambda)
    end function open_smallnu_action_u

    function open_smallnu_action_u_from_loglambda(l, nu, chi, loglambda) result(u)
    ! Action-matched reduced radial approximation.
    !
    ! It keeps the same scaled K_{i nu} comparison solution as
    ! open_smallnu_u, but replaces x=2 Lambda exp(-chi) by an argument
    ! x_*(chi) whose comparison action matches the exact open action, with
    ! the additive constant chosen so x_* -> 2 Lambda exp(-chi) at infinity.
    !
    ! No Liouville amplitude factor is applied.
    integer,  intent(in) :: l
    real(dp), intent(in) :: nu, chi, loglambda
    real(dp) :: u, x

    if (l < 1 .or. abs(nu) < nu_zero) then
        u = comparison_u_from_loglambda(nu, chi, loglambda)
        return
    end if

    x = action_argument_x(l, nu, chi, loglambda)
    if (x <= 0.0_dp) then
        ! Extremely far into the oscillatory tail with underflowed x.  The
        ! non-action version can still evaluate the leading small-x phase from
        ! log Lambda - chi, and the two coincide asymptotically.
        u = comparison_u_from_loglambda(nu, chi, loglambda)
    else
        u = scaled_ki_nu(abs(nu), x)
    end if
    end function open_smallnu_action_u_from_loglambda


    function open_smallnu_u(l, nu, chi, ok) result(u)
    ! Validated small-nu reduced radial approximation for the regular open
    ! hyperspherical Bessel.  The high-l branch applies only the leading
    ! Liouville amplitude A0 = sqrt(asinh(csch chi)/csch chi), which was
    ! calibrated separately from the action-matched x_* argument.
    integer, intent(in) :: l
    real(dp), intent(in) :: nu, chi
    logical, intent(out) :: ok
    real(dp) :: u, loglambda

    ok = open_smallnu_ok(l, nu)
    if (.not. ok) then
        u = 0.0_dp
        return
    end if

    loglambda = open_smallnu_loglambda(l, nu)
    u = open_smallnu_u_from_loglambda(l, nu, chi, loglambda)
    end function open_smallnu_u


    function open_smallnu_u_from_loglambda(l, nu, chi, loglambda) result(u)
    integer, intent(in) :: l
    real(dp), intent(in) :: nu, chi, loglambda
    real(dp) :: u, x, kval

    if (l >= smallnu_amp_l_min) then
        x = action_argument_x(l, nu, chi, loglambda)
        if (x <= 0.0_dp) then
            kval = comparison_u_from_loglambda(nu, chi, loglambda)
        else
            kval = scaled_ki_nu(abs(nu), x)
        end if
        u = open_smallnu_liouville_amp_a0(l, abs(nu), chi) * kval
    else
        u = open_smallnu_action_u_from_loglambda(l, nu, chi, loglambda)
    end if
    end function open_smallnu_u_from_loglambda


    elemental logical function open_smallnu_ok(l, nu) result(ok)
    integer, intent(in) :: l
    real(dp), intent(in) :: nu

    ok = l > smallnu_l_min .and. nu > 0.0_dp .and. nu <= open_smallnu_gate(l)
    end function open_smallnu_ok


    elemental real(dp) function open_smallnu_gate(l) result(nu_max)
    integer, intent(in) :: l
    real(dp) :: ell, old_gate

    ell = real(l, dp)
    old_gate = max(1.2e-4_dp * ell**1.5_dp, &
        min(12.0_dp * (ell / 1000.0_dp)**3, 0.032_dp * ell))
    if (l >= smallnu_amp_l_min) then
        nu_max = max(old_gate, 0.04_dp * ell, &
            min(8.0_dp * (ell / 1000.0_dp)**2, 0.16_dp * ell))
    else
        nu_max = old_gate
    end if
    end function open_smallnu_gate


    function action_argument_x(l, nu, chi, loglambda) result(x)
    ! Return the action-matched comparison argument x_*(chi).
    integer,  intent(in) :: l
    real(dp), intent(in) :: nu, chi, loglambda
    real(dp) :: x, b, A, sh, y, signed_action, cphase, target

    b = abs(nu)
    if (l < 1 .or. b < nu_zero) then
        if (loglambda - chi > 700.0_dp) then
            x = huge(1.0_dp)
        else if (loglambda - chi < -745.0_dp) then
            x = 0.0_dp
        else
            x = 2.0_dp*exp(loglambda - chi)
        end if
        return
    end if

    A = sqrt(real(l,dp)*(real(l,dp) + 1.0_dp))
    if (chi < 350.0_dp) then
        sh = sinh(max(chi, 1.0e-300_dp))
        y = A/sh
    else
        y = 2.0_dp*A*exp(-chi)
    end if

    ! Constant that enforces x_* ~ 2 Lambda exp(-chi) as chi -> infinity.
    cphase = b*(0.5_dp*log(A*A + b*b) - loglambda) - b + A*atan(b/A)

    if (y < b) then
        signed_action = true_action_osc(A, b, y)
    else
        signed_action = -true_action_forb(A, b, y)
    end if

    target = signed_action + cphase
    if (target >= 0.0_dp) then
        x = inv_comp_osc(b, target)
    else
        x = inv_comp_forb(b, -target)
    end if
    end function action_argument_x


    function true_action_osc(A, b, y) result(act)
    ! int_y^b sqrt(b^2-t^2)/(t sqrt(1+t^2/A^2)) dt, for 0 <= y <= b.
    real(dp), intent(in) :: A, b, y
    real(dp) :: act, yy, s, arg1, arg2, y_over_a

    yy = max(0.0_dp, min(y, b))
    if (yy <= 1.0e-4_dp*b) then
        ! Stable y -> 0 limit.  The next correction is O(y^2).
        yy = max(yy, tiny(1.0_dp))
        act = b*log(2.0_dp*b*A/(yy*sqrt(A*A+b*b))) - A*atan(b/A)
        return
    end if

    s = sqrt(max(0.0_dp, b*b - yy*yy))
    y_over_a = yy/A
    arg1 = s/(b*sqrt(1.0_dp + y_over_a*y_over_a))
    arg1 = min(max(arg1, 0.0_dp), 1.0_dp - 4.0_dp*epsilon(1.0_dp))
    arg2 = s/sqrt(A*A + b*b)
    arg2 = min(max(arg2, 0.0_dp), 1.0_dp)
    act = b*atanh(arg1) - A*asin(arg2)
    end function true_action_osc


    function true_action_forb(A, b, y) result(act)
    ! int_b^y sqrt(t^2-b^2)/(t sqrt(1+t^2/A^2)) dt, for y >= b.
    real(dp), intent(in) :: A, b, y
    real(dp) :: act, yy, s, arg1, arg2

    yy = max(y, b)
    s = sqrt(max(0.0_dp, yy*yy - b*b))
    arg1 = s/sqrt(A*A + yy*yy)
    arg1 = min(max(arg1, 0.0_dp), 1.0_dp - 4.0_dp*epsilon(1.0_dp))
    if (yy > 0.0_dp) then
        arg2 = A*s/(yy*sqrt(A*A + b*b))
    else
        arg2 = 0.0_dp
    end if
    arg2 = min(max(arg2, 0.0_dp), 1.0_dp)
    act = A*atanh(arg1) - b*asin(arg2)
    end function true_action_forb


    function inv_comp_osc(b, target) result(x)
    ! Invert int_x^b sqrt(b^2-t^2)/t dt = target without bisection.
    !
    ! Write s=sqrt(1-(x/b)^2).  Then
    !    target/b = atanh(s) - s,       0 <= s < 1,
    ! and x=b*sqrt(1-s^2).  Four Halley steps from asymptotic/series
    ! starts are overkill for the intended 1e-4 peak-normalized branch.
    real(dp), intent(in) :: b, target
    real(dp) :: x, T, s, p, p2, one, f, fp, fpp, den, ds, s_old, logx
    integer :: it

    if (target <= 0.0_dp) then
        x = b
        return
    end if

    T = target/b

    ! Deep oscillatory tail.  From target/b = log(2b/x) - 1 + O((x/b)^2).
    if (T >= 8.0_dp) then
        logx = log(2.0_dp*b) - (T + 1.0_dp)
        if (logx < -745.0_dp) then
            x = 0.0_dp
        else
            x = exp(logx)
        end if
        return
    end if

    if (T < 0.30_dp) then
        ! Reversion of atanh(s)-s = s^3/3 + s^5/5 + ...
        p = (3.0_dp*T)**(1.0_dp/3.0_dp)
        p2 = p*p
        s = p*(1.0_dp - 0.20_dp*p2 + (3.0_dp/175.0_dp)*p2*p2)
    else
        ! Good middle/large start; Halley rapidly removes the bias.
        s = tanh(T + 1.0_dp)
    end if
    s = min(max(s, 0.0_dp), 1.0_dp - 8.0_dp*epsilon(1.0_dp))

    do it = 1, 4
        s_old = s
        one = max(1.0_dp - s*s, tiny(1.0_dp))
        f   = atanh(s) - s - T
        fp  = s*s/one
        fpp = 2.0_dp*s/(one*one)
        den = 2.0_dp*fp*fp - f*fpp
        if (den == 0.0_dp) exit
        ds = 2.0_dp*f*fp/den
        s = s - ds
        if (s <= 0.0_dp) s = 0.5_dp*s_old
        if (s >= 1.0_dp) s = 0.5_dp*(1.0_dp + s_old)
    end do

    x = b*sqrt(max(0.0_dp, 1.0_dp - s*s))
    end function inv_comp_osc


    function inv_comp_forb(b, target) result(x)
    ! Invert int_b^x sqrt(t^2-b^2)/t dt = target without bisection.
    !
    ! Write s=sqrt((x/b)^2-1).  Then
    !    target/b = s - atan(s),        s >= 0,
    ! and x=b*sqrt(1+s^2).
    real(dp), intent(in) :: b, target
    real(dp) :: x, T, s, p, p2, a, one, f, fp, fpp, den, ds, s_old
    integer :: it

    if (target <= 0.0_dp) then
        x = b
        return
    end if

    T = target/b
    if (T < 0.30_dp) then
        ! Reversion of s-atan(s) = s^3/3 - s^5/5 + ...
        p = (3.0_dp*T)**(1.0_dp/3.0_dp)
        p2 = p*p
        s = p*(1.0_dp + 0.20_dp*p2 + (3.0_dp/175.0_dp)*p2*p2)
    else
        ! Large-s asymptotic: T = s - pi/2 + 1/s + O(s^-3).
        a = T + 0.5_dp*pi_dp
        s = a - 1.0_dp/a
    end if
    s = max(s, 0.0_dp)

    do it = 1, 3
        s_old = s
        one = 1.0_dp + s*s
        f   = s - atan(s) - T
        fp  = s*s/one
        fpp = 2.0_dp*s/(one*one)
        den = 2.0_dp*fp*fp - f*fpp
        if (den == 0.0_dp) exit
        ds = 2.0_dp*f*fp/den
        s = s - ds
        if (s < 0.0_dp) s = 0.5_dp*s_old
    end do

    x = b*sqrt(1.0_dp + s*s)
    end function inv_comp_forb

    function comparison_u_from_loglambda(nu, chi, loglambda) result(u)
    ! Reduced radial function, with log Lambda already precomputed.
    ! Centralized through the same scaled-K router as scaled_ki_nu().
    real(dp), intent(in) :: nu, chi, loglambda
    real(dp) :: u, b, logz2

    b = abs(nu)
    logz2 = loglambda - chi        ! x/2 = Lambda exp(-chi)
    u = scaled_ki_from_logz2(b, logz2)
    end function comparison_u_from_loglambda


    function open_smallnu_loglambda(l, nu) result(loglambda)
    ! Branch-continuous evaluation of
    !     log Lambda = Im log Gamma(l+1+i nu)/nu.
    !
    ! For nu/(l+1) small, use the convergent tail expansion
    !     H_l - gamma + sum_{m>=1} (-1)^(m+1) nu^(2m)/(2m+1)
    !                    * sum_{n=l+1}^infty n^{-(2m+1)}.
    ! Otherwise compute Im log Gamma(1+i nu) plus sum atan(nu/j).
    integer,  intent(in) :: l
    real(dp), intent(in) :: nu
    real(dp) :: loglambda
    real(dp) :: b, a, b2, corr, powb
    integer :: j

    b = abs(nu)
    if (b < nu_zero) then
        loglambda = harmonic_minus_euler(l)
        return
    end if

    a = real(l + 1, dp)

    ! This is the regime for which the approximation is intended.  It is
    ! faster and avoids branch issues for very large l.
    if (l >= 8 .and. b/a < 0.30_dp) then
        b2 = b*b
        powb = b2
        corr =  powb/3.0_dp * zeta_tail_em(3, a)
        powb = powb*b2
        corr = corr - powb/5.0_dp * zeta_tail_em(5, a)
        powb = powb*b2
        corr = corr + powb/7.0_dp * zeta_tail_em(7, a)
        powb = powb*b2
        corr = corr - powb/9.0_dp * zeta_tail_em(9, a)
        if (b/a >= 0.20_dp) then
            powb = powb*b2
            corr = corr + powb/11.0_dp * zeta_tail_em(11, a)
            powb = powb*b2
            corr = corr - powb/13.0_dp * zeta_tail_em(13, a)
        end if
        loglambda = harmonic_minus_euler(l) + corr
        return
    end if

    loglambda = arg_gamma_1_plus_i(b)
    do j = 1, l
        loglambda = loglambda + atan(b/real(j,dp))
    end do
    loglambda = loglambda/b
    end function open_smallnu_loglambda


    function scaled_ki_nu(nu, x) result(val)
    ! Standalone evaluator of
    !    sqrt(sinh(pi nu)/(pi nu))*K_{i nu}(x), x>0.
    !
    ! Accuracy target: fast O(1e-4) peak-normalized use, not full double
    ! precision everywhere.
    !   - I_{i nu} power series while q=x^2/(4 nu) <= series_q_max;
    !   - uniform Airy turning-point approximation otherwise.
    real(dp), intent(in) :: nu, x
    real(dp) :: val, b, logz2

    b = abs(nu)
    if (x <= 0.0_dp) then
        val = ieee_nan()
        return
    end if

    if (b < nu_zero) then
        val = k0_approx(x)
        return
    end if

    logz2 = log(0.5_dp*x)
    val = scaled_ki_from_xlog(b, x, logz2)
    end function scaled_ki_nu


    function scaled_ki_from_logz2(nu, logz2) result(val)
    ! Same evaluator as scaled_ki_nu(), but accepts log(x/2).  This avoids
    ! overflow/underflow in the open_smallnu_* wrappers.
    real(dp), intent(in) :: nu, logz2
    real(dp) :: val, b, x

    b = abs(nu)
    if (b < nu_zero) then
        val = k0_from_logz2(logz2)
        return
    end if

    ! x is so large that the scaled K is negligible for this approximation.
    if (logz2 > 120.0_dp) then
        val = 0.0_dp
        return
    end if

    ! Endpoint oscillatory limit.  Corrections are O(x^2/b), utterly
    ! negligible here; this also handles x underflow without forming x.
    if (logz2 < -20.0_dp) then
        val = scaled_ki_series(b, logz2, 0.0_dp, .true.)
        return
    end if

    x = 2.0_dp*exp(logz2)
    val = scaled_ki_from_xlog(b, x, logz2)
    end function scaled_ki_from_logz2


    function scaled_ki_from_xlog(nu, x, logz2) result(val)
    ! Core router.  nu must be positive and x>0.
    real(dp), intent(in) :: nu, x, logz2
    real(dp) :: val, b, q, x_series_max, x_asymp_min
    logical :: ok

    b = abs(nu)

    if (b < nu_airy_min) then
        ! Low nu: keep the original conservative routing.  The real-axis
        ! integral is only used here, where it is not multiplied by an
        ! astronomically large exp(pi*nu/2) factor.
        x_series_max = min(20.0_dp, max(3.0_dp, 0.85_dp*b + 1.0_dp))
        if (x <= x_series_max) then
            call scaled_ki_series_checked(b, logz2, x, .false., val, ok)
            if (ok) return
        end if

        x_asymp_min = 45.0_dp + 4.0_dp*b*b
        if (x > x_asymp_min) then
            val = scaled_ki_asymptotic(b, x)
        else
            val = scaled_ki_integral(b, x)
        end if
        return
    end if

    ! Large nu: the I-series is excellent and fast while q is modest, but
    ! it becomes a cancellation machine near and beyond the turning point.
    ! On the forbidden side, Airy is both faster and
    ! more accurate once we are a few units past the turning point.
    q = 0.25_dp*x*x/b

    ! If x is vastly larger than nu, the standard large-x K expansion is
    ! even cheaper than Airy and safely exponentially small.
    x_asymp_min = 45.0_dp + 4.0_dp*b*b
    if (x > x_asymp_min) then
        val = scaled_ki_asymptotic(b, x)
        return
    end if

    if (x > b + 5.0_dp) then
        val = scaled_ki_airy_leading(b, x)
        return
    end if

    if (q <= series_q_max) then
        call scaled_ki_series_checked(b, logz2, x, .false., val, ok)
        if (ok) return
    end if

    val = scaled_ki_airy_leading(b, x)
    end function scaled_ki_from_xlog


    function scaled_ki_series(nu, logz2, x, leading_only) result(val)
    ! Compatibility wrapper around the checked series evaluator.
    real(dp), intent(in) :: nu, logz2, x
    logical,  intent(in) :: leading_only
    real(dp) :: val
    logical :: ok

    call scaled_ki_series_checked(nu, logz2, x, leading_only, val, ok)
    end function scaled_ki_series


    subroutine scaled_ki_series_checked(nu, logz2, x, leading_only, val, ok)
    ! Scaled imaginary-order Macdonald function from the I_{i nu} series.
    !
    ! I_{i nu}(x) = exp(i nu log(x/2) - log Gamma(1+i nu))
    !                 * sum_k (x^2/4)^k/[k! (1+i nu)_k]
    !
    ! Multiplication by sqrt(sinh(pi nu)/(pi nu)) cancels the modulus of
    ! Gamma(1+i nu), giving the stable expression
    !
    ! scaled K = - Im[ exp(i phase) * series ]/nu,
    ! phase = nu log(x/2) - arg Gamma(1+i nu).
    real(dp), intent(in) :: nu, logz2, x
    logical,  intent(in) :: leading_only
    real(dp), intent(out) :: val
    logical,  intent(out) :: ok
    real(dp) :: b, y, phase, cph, sph
    real(dp) :: term_re, term_im, sum_re, sum_im
    real(dp) :: fac_re, fac_im, new_re, new_im, kk, den, term2, sum2
    integer :: k
    integer, parameter :: max_series_iter = 384

    b = abs(nu)
    ok = .true.
    if (b < nu_zero) then
        val = k0_from_logz2(logz2)
        return
    end if

    phase = b*logz2 - arg_gamma_1_plus_i(b)
    phase = modulo(phase + pi_dp, twopi_dp) - pi_dp
    cph = cos(phase)
    sph = sin(phase)

    if (leading_only) then
        val = -sph/b
        return
    end if

    y = 0.25_dp*x*x
    term_re = 1.0_dp
    term_im = 0.0_dp
    sum_re  = 1.0_dp
    sum_im  = 0.0_dp
    ok = .false.

    do k = 1, max_series_iter
        kk = real(k,dp)
        den = kk*kk + b*b
        fac_re = y/den
        fac_im = -y*b/(kk*den)
        new_re = term_re*fac_re - term_im*fac_im
        new_im = term_re*fac_im + term_im*fac_re
        term_re = new_re
        term_im = new_im
        sum_re = sum_re + term_re
        sum_im = sum_im + term_im

        term2 = term_re*term_re + term_im*term_im
        sum2 = sum_re*sum_re + sum_im*sum_im
        if (term2 <= series_term2_tol*max(1.0_dp, sum2)) then
            ok = .true.
            exit
        end if

        ! Avoid silently returning nonsense if an accidental call sends the
        ! series into the cancellation-prone region.
        if (.not. ieee_is_finite(sum_re) .or. .not. ieee_is_finite(sum_im)) exit
    end do

    val = -(cph*sum_im + sph*sum_re)/b
    end subroutine scaled_ki_series_checked


    function scaled_ki_airy_leading(nu, x) result(val)
    ! Leading uniform Airy approximation for
    !   sqrt(sinh(pi nu)/(pi nu))*K_{i nu}(x), nu > 20 large.
    real(dp), intent(in) :: nu, x
    real(dp) :: val
    real(dp) :: b, z, s, t, act, zeta_abs, phi, arg, ai
    real(dp) :: logb, b23, scale

    real(dp), parameter :: one_third  = 0.333333333333333333333333333333333333333_dp
    real(dp), parameter :: one_fifth  = 0.2_dp
    real(dp), parameter :: one_seventh = 0.142857142857142857142857142857142857143_dp
    real(dp), parameter :: one_ninth  = 0.111111111111111111111111111111111111111_dp
    real(dp), parameter :: one_eleventh = 0.090909090909090909090909090909090909091_dp
    real(dp), parameter :: two_third  = 0.666666666666666666666666666666666666667_dp
    real(dp), parameter :: five_sixth = 0.833333333333333333333333333333333333333_dp
    real(dp), parameter :: sqrt_pi_over_two = 1.253314137315500251207882642405522626504_dp
    real(dp), parameter :: phi_turn = 1.259921049894873164767210607278228350570_dp
    real(dp), parameter :: s_series_max = 0.10_dp

    b = abs(nu)
    if (b < nu_zero .or. x <= 0.0_dp) then
        val = ieee_nan()
        return
    end if

    logb = log(b)

    ! b^(2/3), used in the Airy argument.
    b23 = exp(two_third*logb)

    if (pi_dp*b >= 20.0_dp) then
        !log_sinh(pi*b) = pi*b - log(2) to double precision here
        scale = sqrt_pi_over_two * exp(-five_sixth*logb)
    else
        scale = exp(0.5_dp*(log_sinh(pi_dp*b) - log(pi_dp*b)) &
            + log(pi_dp) - 0.5_dp*pi_dp*b - logb/3.0_dp)
    end if

    z = x/b

    if (abs(z - 1.0_dp) < 1.0e-8_dp) then
        arg = 0.0_dp
        phi = phi_turn

    else if (z < 1.0_dp) then
        s = sqrt(max(0.0_dp, (1.0_dp - z)*(1.0_dp + z)))
        t = s*s

        if (s < s_series_max) then
            ! atanh(s)-s = s^3/3 + s^5/5 + s^7/7 + ...
            act = s*t*(one_third + t*(one_fifth + t*(one_seventh + &
                t*(one_ninth + t*one_eleventh))))
        else
            act = atanh(s) - s
        end if

        zeta_abs = exp(two_third*log(1.5_dp*act))
        phi = sqrt(sqrt(4.0_dp*zeta_abs/t))
        arg = -b23*zeta_abs

    else
        s = sqrt(max(0.0_dp, (z - 1.0_dp)*(z + 1.0_dp)))
        t = s*s

        if (s < s_series_max) then
            ! s-atan(s) = s^3/3 - s^5/5 + s^7/7 - ...
            act = s*t*(one_third - t*(one_fifth - t*(one_seventh - &
                t*(one_ninth - t*one_eleventh))))
        else
            act = s - atan(s)
        end if

        zeta_abs = exp(two_third*log(1.5_dp*act))
        phi = sqrt(sqrt(4.0_dp*zeta_abs/t))
        arg = b23*zeta_abs
    end if

    ! Ai(arg) is already below double underflow once
    ! (2/3)*arg^(3/2) > ~745, i.e. arg ~= 108.
    if (arg > 110.0_dp) then
        val = 0.0_dp
        return
    end if

    ai = airy_ai_fast(arg)
    val = scale*phi*ai
    end function scaled_ki_airy_leading


    function scaled_ki_integral(nu, x) result(val)
    ! Fallback quadrature for the scaled K_{i nu}:
    !   K_{i nu}(x) = int_0^infty exp(-x cosh t) cos(nu t) dt.
    !
    ! Uses composite 16-point Gauss-Legendre quadrature over [0,tmax].
    ! This branch is deliberately robust rather than maximally fast.
    real(dp), intent(in) :: nu, x
    real(dp) :: val, b, tmax, target, h, a, c, halfw, mid, integ
    real(dp) :: logscale, exponent
    integer :: n_panel, p, i

    real(dp), parameter :: gx(8) = [ &
        0.095012509837637440185319335424958063130_dp, &
        0.281603550779258913230460501460496106486_dp, &
        0.458016777657227386342419442983577573540_dp, &
        0.617876244402643748446671764048791018991_dp, &
        0.755404408355003033895101194847442268354_dp, &
        0.865631202387831743880467897712393132387_dp, &
        0.944575023073232576077988415534608345091_dp, &
        0.989400934991649932596154173450332627426_dp ]
    real(dp), parameter :: gw(8) = [ &
        0.189450610455068496285396723208283105146_dp, &
        0.182603415044923588866763667969219939384_dp, &
        0.169156519395002538189312079030359962211_dp, &
        0.149595988816576732081501730547478548970_dp, &
        0.124628971255533872052476282192016420144_dp, &
        0.095158511682492784809925107602246226355_dp, &
        0.062253523938647892862843836994377694274_dp, &
        0.027152459411754094851780572456018103512_dp ]

    b = abs(nu)
    if (b < nu_zero) then
        val = k0_approx(x)
        return
    end if

    ! If x is already very large, use the asymptotic branch.
    if (x > 45.0_dp + 4.0_dp*b*b) then
        val = scaled_ki_asymptotic(b, x)
        return
    end if

    target = 45.0_dp
    if (x >= target) then
        tmax = sqrt(max(2.0_dp*(target/x - 1.0_dp), 1.0e-12_dp))
    else
        tmax = acosh_safe(target/x)
    end if
    tmax = max(tmax, 0.5_dp)

    ! Resolve both the smooth exp(-x cosh t) envelope and cos(nu t).
    if (b > 0.0_dp) then
        h = min(0.25_dp, pi_dp/(10.0_dp*b))
    else
        h = 0.25_dp
    end if
    h = max(h, 0.0025_dp)
    n_panel = max(1, ceiling(tmax/h))
    h = tmax/real(n_panel,dp)

    integ = 0.0_dp
    do p = 0, n_panel-1
        a = real(p,dp)*h
        c = a + h
        mid = 0.5_dp*(a+c)
        halfw = 0.5_dp*h
        do i = 1, 8
            integ = integ + halfw*gw(i) * &
                ( exp(-x*cosh(mid + halfw*gx(i))) * cos(b*(mid + halfw*gx(i))) + &
                exp(-x*cosh(mid - halfw*gx(i))) * cos(b*(mid - halfw*gx(i))) )
        end do
    end do

    logscale = 0.5_dp*(log_sinh(pi_dp*b) - log(pi_dp*b))
    if (integ == 0.0_dp) then
        val = 0.0_dp
    else
        exponent = logscale + log(abs(integ))
        if (exponent < log_tiny_safe) then
            val = 0.0_dp
        else if (exponent > 700.0_dp) then
            ! Outside the intended parameter range for this simple fallback.
            val = sign(huge(1.0_dp), integ)
        else
            val = sign(exp(exponent), integ)
        end if
    end if
    end function scaled_ki_integral


    function scaled_ki_asymptotic(nu, x) result(val)
    ! Large-x asymptotic for scaled K_{i nu}(x).
    ! K_v(x) ~ sqrt(pi/(2x))*exp(-x) * sum_m a_m/x^m,
    ! with mu = 4 v^2 = -4 nu^2.
    real(dp), intent(in) :: nu, x
    real(dp) :: val, b, mu, term, sum, logscale, exponent
    integer :: m

    b = abs(nu)
    mu = -4.0_dp*b*b
    term = 1.0_dp
    sum = 1.0_dp
    do m = 1, 80
        term = term * (mu - real((2*m-1)*(2*m-1),dp)) / (real(m,dp)*8.0_dp*x)
        if (abs(term) > 0.5_dp*abs(sum) .and. m > 8) exit
        sum = sum + term
        if (abs(term) <= 2.0e-16_dp*abs(sum)) exit
    end do

    logscale = 0.5_dp*(log_sinh(pi_dp*b) - log(pi_dp*b))
    exponent = logscale + 0.5_dp*(log(pi_dp) - log(2.0_dp*x)) - x
    if (exponent < log_tiny_safe) then
        val = 0.0_dp
    else
        val = exp(exponent)*sum
    end if
    end function scaled_ki_asymptotic


    function k0_from_logz2(logz2) result(val)
    ! K_0(2 exp(logz2)), retaining accuracy when the argument underflows.
    real(dp), intent(in) :: logz2
    real(dp) :: val, x

    if (logz2 < -20.0_dp) then
        ! K0(x) = -log(x/2)-gamma + O(x^2 log x), x=2 exp(logz2)
        val = -logz2 - euler_gamma
    else if (logz2 > 120.0_dp) then
        val = 0.0_dp
    else
        x = 2.0_dp*exp(logz2)
        val = k0_approx(x)
    end if
    end function k0_from_logz2


    pure function open_smallnu_liouville_amp_a0(l, nu, chi) result(amp)
    ! Leading fast Liouville amplitude for the high-l small-nu branch:
    ! A0 = sqrt(asinh(eta)/eta), eta = csch(chi).  The higher-order
    ! exp((nu/sqrt(l(l+1)))^2 H) correction is deliberately omitted; it
    ! worsened the calibrated mid-l gate while giving negligible benefit.
    integer, intent(in) :: l
    real(dp), intent(in) :: nu, chi
    real(dp) :: amp
    real(dp), parameter :: chi_large = 40.0_dp
    real(dp) :: eta, p

    if (l <= 0 .or. nu <= 0.0_dp .or. chi <= 0.0_dp) then
        amp = 1.0_dp
        return
    end if

    if (chi > chi_large) then
        eta = 2.0_dp * exp(-chi)
    else
        eta = 1.0_dp / sinh(chi)
    end if

    if (eta < 1.e-3_dp) then
        amp = 1.0_dp
        return
    end if

    if (eta > 1.e100_dp) then
        p = log(eta) + log_two
    else
        p = log(eta + sqrt(1.0_dp + eta * eta))
    end if

    amp = sqrt(p / eta)
    end function open_smallnu_liouville_amp_a0


    function k0_approx(x) result(ans)
    ! Modified Bessel K0 approximation, adapted from the classic Cephes/NR
    ! piecewise rational-polynomial form.  Sufficient for the nu -> 0
    ! fallback of this approximation branch.
    real(dp), intent(in) :: x
    real(dp) :: ans, y

    if (x <= 0.0_dp) then
        ans = huge(1.0_dp)
    else if (x <= 2.0_dp) then
        y = 0.25_dp*x*x
        ans = -log(0.5_dp*x)*i0_approx(x) + &
            (-0.57721566_dp + y*(0.42278420_dp + y*(0.23069756_dp + &
            y*(0.03488590_dp + y*(0.00262698_dp + y*(0.00010750_dp + &
            y*0.00000740_dp))))))
    else
        y = 2.0_dp/x
        ans = exp(-x)/sqrt(x) * &
            (1.25331414_dp + y*(-0.07832358_dp + y*(0.02189568_dp + &
            y*(-0.01062446_dp + y*(0.00587872_dp + y*(-0.00251540_dp + &
            y*0.00053208_dp))))))
    end if
    end function k0_approx


    function i0_approx(x) result(ans)
    real(dp), intent(in) :: x
    real(dp) :: ans, ax, y

    ax = abs(x)
    if (ax < 3.75_dp) then
        y = (x/3.75_dp)**2
        ans = 1.0_dp + y*(3.5156229_dp + y*(3.0899424_dp + y*(1.2067492_dp + &
            y*(0.2659732_dp + y*(0.0360768_dp + y*0.0045813_dp)))))
    else
        y = 3.75_dp/ax
        ans = exp(ax)/sqrt(ax) * &
            (0.39894228_dp + y*(0.01328592_dp + y*(0.00225319_dp + &
            y*(-0.00157565_dp + y*(0.00916281_dp + y*(-0.02057706_dp + &
            y*(0.02635537_dp + y*(-0.01647633_dp + y*0.00392377_dp))))))))
    end if
    end function i0_approx


    function arg_gamma_1_plus_i(nu) result(arg)
    ! Imaginary part of log Gamma(1+i nu), evaluated with the same
    ! g=7, n=9 Lanczos approximation using real arithmetic specialized to z = 1+i nu.
    real(dp), intent(in) :: nu
    real(dp) :: arg, b, t_re, logabs_t, theta_t, x_re, x_im
    real(dp) :: a, den
    integer :: i
    real(dp), parameter :: g = 7.0_dp
    real(dp), parameter :: p(9) = [ &
        0.99999999999980993227684700473478_dp, &
        676.520368121885098567009190444019_dp, &
        -1259.13921672240287047156078755283_dp, &
        771.3234287776530788486528258894_dp, &
        -176.61502916214059906584551354_dp, &
        12.507343278686904814458936853_dp, &
        -0.13857109526572011689554707_dp, &
        0.000009984369578019570859563e0_dp, &
        0.000000150563273514931155834e0_dp ]

    b = abs(nu)
    t_re = g + 0.5_dp
    logabs_t = 0.5_dp*log(t_re*t_re + b*b)
    theta_t = atan2(b, t_re)

    x_re = p(1)
    x_im = 0.0_dp
    do i = 2, 9
        a = real(i-1, dp)
        den = a*a + b*b
        x_re = x_re + p(i)*a/den
        x_im = x_im - p(i)*b/den
    end do

    arg = 0.5_dp*theta_t + b*logabs_t - b + atan2(x_im, x_re)
    end function arg_gamma_1_plus_i


    function harmonic_minus_euler(l) result(hmg)
    ! H_l - gamma, with a direct sum for small l and an asymptotic expansion
    ! for large l.  For l=0 this is -gamma.
    integer, intent(in) :: l
    real(dp) :: hmg, n, invn, invn2
    integer :: j

    if (l <= 0) then
        hmg = -euler_gamma
    else if (l < 10000) then
        hmg = -euler_gamma
        do j = 1, l
            hmg = hmg + 1.0_dp/real(j,dp)
        end do
    else
        n = real(l,dp)
        invn = 1.0_dp/n
        invn2 = invn*invn
        ! H_l = log l + gamma + 1/(2l) - 1/(12l^2)
        !       + 1/(120l^4) - 1/(252l^6) + ...
        hmg = log(n) + 0.5_dp*invn - invn2/12.0_dp &
            + invn2*invn2/120.0_dp - invn2*invn2*invn2/252.0_dp
    end if
    end function harmonic_minus_euler


    function zeta_tail_em(p, a) result(s)
    ! Euler-Maclaurin estimate of sum_{n=a}^infty n^{-p}, p>1.
    ! Good for a >= roughly 8; only used in that regime above.
    integer,  intent(in) :: p
    real(dp), intent(in) :: a
    real(dp) :: s, pp

    pp = real(p,dp)
    s = a**(1.0_dp-pp)/(pp-1.0_dp) &
        + 0.5_dp*a**(-pp) &
        + (pp/12.0_dp)*a**(-pp-1.0_dp) &
        - (pp*(pp+1.0_dp)*(pp+2.0_dp)/720.0_dp)*a**(-pp-3.0_dp) &
        + (pp*(pp+1.0_dp)*(pp+2.0_dp)*(pp+3.0_dp)*(pp+4.0_dp)/30240.0_dp)*a**(-pp-5.0_dp)
    end function zeta_tail_em


    function log_sinh(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y

    if (x < 20.0_dp) then
        y = log(sinh(x))
    else
        y = x - log_two + log(1.0_dp - exp(-2.0_dp*x))
    end if
    end function log_sinh


    function acosh_safe(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y, xx

    xx = max(x, 1.0_dp)
    y = log(xx + sqrt(max(0.0_dp, xx*xx - 1.0_dp)))
    end function acosh_safe


    function ieee_nan() result(x)
    real(dp) :: x
    x = ieee_value(1.0_dp, ieee_quiet_nan)
    end function ieee_nan

    end module HypersphericalBesselSmallNu
