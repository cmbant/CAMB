    module HypersphericalBesselUtils
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64
    real(dp), parameter :: PI = acos(-1.0_dp)
    real(dp), parameter :: CACHE_EPS = 1.0e-12_dp

    public :: normalize_chi, turning_point, curved_radius, qintegral_exact

    contains

    elemental subroutine normalize_chi(l, K, nu, chi, achi, symm)
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, chi
    real(dp), intent(out) :: achi, symm
    integer :: inu

    achi = abs(chi)
    symm = 1.0_dp

    if (K /= 1) return

    inu = nint(nu)
    achi = modulo(achi, 2.0_dp * PI)
    if (achi > PI) then
        achi = 2.0_dp * PI - achi
        if (mod(l, 2) /= 0) symm = -symm
    end if
    if (achi > PI / 2.0_dp) then
        achi = PI - achi
        if (mod(inu - l - 1, 2) /= 0) symm = -symm
    end if
    end subroutine normalize_chi


    elemental real(dp) function turning_point(ell, nu, K)
    real(dp), intent(in) :: ell, nu
    integer, intent(in) :: K

    select case (K)
    case (-1)
        turning_point = asinh(ell / nu)
    case (0)
        turning_point = ell / nu
    case (1)
        turning_point = asin(ell / nu)
    case default
        turning_point = 0.0_dp
    end select
    end function turning_point


    elemental real(dp) function curved_radius(K, chi)
    integer, intent(in) :: K
    real(dp), intent(in) :: chi

    select case (K)
    case (-1)
        curved_radius = sinh(chi)
    case (0)
        curved_radius = chi
    case (1)
        curved_radius = sin(chi)
    case default
        curved_radius = chi
    end select
    end function curved_radius


    elemental real(dp) function qintegral_exact(sin_K, alpha, K) result(q)
    real(dp), intent(in) :: sin_K, alpha
    integer, intent(in) :: K

    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, half = 0.5_dp
    real(dp) :: x, x2, a2, ha, r1, r2, u, m

    x  = alpha * sin_K
    x2 = x * x
    a2 = alpha * alpha
    ha = half * alpha

    select case (K)

    case (0)
        if (x > one) then
            r1 = sqrt(x2 - one)
            q = self_minus_atan(r1)
        else
            r1 = sqrt(max(zero, one - x2))
            if (r1 < one) then
                q = atanh_minus_self(r1)
            else
                q = log(two / max(x, CACHE_EPS)) - one
            end if
        end if

    case (-1)
        if (x > one) then
            r1 = sqrt(x2 - one)
            r2 = sqrt(x2 + a2)
            q = alpha * asinh(r1 / sqrt(one + a2)) - atan2(alpha*r1, r2)
        else
            u = sqrt(max(zero, one - x2))
            r2 = sqrt(x2 + a2)
            r1 = u * r2
            q = ha * atan2(-two*r1, two*x2 + a2 - one) &
                + asinh(alpha*u / max(x*sqrt(one + a2), CACHE_EPS))
        end if

    case default
        if (x > one) then
            r1 = sqrt(x2 - one)
            r2 = sqrt(max(zero, a2 - x2))
            q = ha * atan2(two*r1*r2, a2 + one - two*x2) &
                - atan2(alpha*r1, r2)
        else
            m = a2 - one
            if (m > CACHE_EPS) then
                u = sqrt(max(zero, one - x2))
                q = asinh(alpha*u / max(x*sqrt(m), CACHE_EPS)) - alpha*asinh(u / sqrt(m))
            else
                q = -half*log(max(x2, CACHE_EPS))
            end if
        end if
    end select
    end function qintegral_exact


    elemental real(dp) function atanh_minus_self(y) result(res)
    real(dp), intent(in) :: y

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: small_cut = 1.0e-3_dp
    real(dp) :: y2

    if (abs(y) < small_cut) then
        y2 = y * y
        res = y * y2 * (one/3.0_dp + y2 * (one/5.0_dp + y2 * (one/7.0_dp)))
    else
        res = atanh(y) - y
    end if
    end function atanh_minus_self


    elemental real(dp) function self_minus_atan(y) result(res)
    real(dp), intent(in) :: y

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: small_cut = 1.0e-3_dp
    real(dp) :: y2

    if (abs(y) < small_cut) then
        y2 = y * y
        res = y * y2 * (one/3.0_dp + y2 * (-one/5.0_dp + y2 * (one/7.0_dp)))
    else
        res = y - atan(y)
    end if
    end function self_minus_atan

    end module HypersphericalBesselUtils


    module HypersphericalBesselAiry
    ! One-point second-order Airy/Olver approximation for hyperspherical Bessel functions.
    ! It solves the reduced equation u'' = lambda^2 q(chi) u + r(chi) u with
    ! u = S_K(chi) phi_l^nu(chi), normalized to the exact origin amplitude.
    ! The local second-order correction is calibrated to about 1e-4 peak-relative
    ! accuracy in its gate against phi_recurs.
    ! Assume L > 0 and other variables already checked for physical limits
    use iso_fortran_env, only: real64
    use HypersphericalBesselUtils, only: curved_radius, qintegral_exact, turning_point, normalize_chi
    use MathUtils, only: airy_fast
    implicit none
    private

    integer, parameter :: dp = real64
    real(dp), parameter :: PI = acos(-1.0_dp)
    real(dp), parameter :: LOG2PI = log(2.0_dp * PI)
    real(dp), parameter :: CACHE_EPS = 1.0e-12_dp

    ! Fast second-order one-point Airy/Olver patch.  The correction is
    ! evaluated from local samples on the single requested chi segment; no
    ! zeta->chi inversion, dense solve, or whole-window residual fit is used.
    ! The calibrated validity gate is set by the low-nu open and near-degenerate
    ! closed limits; the positive Airy tail is peak-negligible and returns zero.
    integer, parameter :: AIRY_SECOND_FAST_DEG = 4
    real(dp), parameter :: AIRY_SECOND_TAU(1:AIRY_SECOND_FAST_DEG) = [ &
        1.4644660940672625e-1_dp, &
        5.0000000000000000e-1_dp, &
        8.5355339059327373e-1_dp, &
        1.0000000000000000e0_dp ]
    integer, parameter :: AIRY_SECOND_L_MIN = 10
    real(dp), parameter :: AIRY_SECOND_TAIL_ZERO_X = 25.77_dp
    real(dp), parameter :: AIRY_SECOND_OPEN_NU_MIN = 7.0_dp
    real(dp), parameter :: AIRY_SECOND_CLOSED_GATE = 15.0_dp
    real(dp), parameter :: AIRY_PSI_ZETA_LOCAL = 2.0e-3_dp
    real(dp), parameter :: AIRY_SECOND_FIT_ZETA_MIN = 1.6e-2_dp

    public :: airy_u, airy_u_normalized
    public :: airy_ok

    contains

    pure real(dp) function log_origin_u0_fast(l, K, nu) result(logu0)
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu

    real(dp) :: logdf, logprod

    logdf = real(l + 1, dp) * log(2.0_dp) &
        + log_gamma(real(l, dp) + 1.5_dp) &
        - 0.5_dp * log(PI)

    select case (K)
    case (0)
        logprod = real(l, dp) * log(nu)
    case (1)
        logprod = 0.5_dp * ( log_gamma(nu + l + 1) - log_gamma(nu - l) - log(nu) )
    case (-1)
        logprod = 0.5_dp * log_prod_plus_one(l, nu)
    end select

    logu0 = logprod - logdf
    end function log_origin_u0_fast


    pure real(dp) function airy_q0(K, beta) result(q0)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: x, b2, invb

    b2 = beta * beta

    select case (K)
    case (1)
        if (beta <= 1.5_dp) then
            ! Algebraically equivalent to:
            ! log(2/sqrt(beta^2-1)) - beta*atanh(1/beta)
            ! but stable as beta -> 1+.
            x = beta - 1.0_dp

            if (x <= 0.0_dp) then
                q0 = 0.0_dp
            else
                q0 = log(2.0_dp) &
                    + 0.5_dp * ( x * log(x) &
                    - (beta + 1.0_dp) * log(beta + 1.0_dp) )
            end if
        else
            invb = 1.0_dp / beta
            q0 = log(2.0_dp) - log(beta) &
                - 0.5_dp * log(1.0_dp - invb * invb) &
                - beta * atanh(invb)
        end if

    case (-1)
        if (beta < 1.0_dp) then
            q0 = log(2.0_dp) &
                - 0.5_dp * log(1.0_dp + beta * beta) &
                - beta * atan2(1.0_dp, beta)
        else
            invb = 1.0_dp / beta
            q0 = log(2.0_dp) - log(beta) &
                - 0.5_dp * log(1.0_dp + invb * invb) &
                - beta * atan(invb)
        end if
    case default
        q0 = log(2.0_dp / beta) - 1.0_dp
    end select
    end function airy_q0


    pure real(dp) function airy_B0_turn(K, beta) result(b0)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: rk, b2, denom

    rk = real(K, dp)
    b2 = beta * beta

    denom = 280.0_dp * beta**(4.0_dp/3.0_dp) &
        * max(b2 - rk, tiny(1.0_dp))**(4.0_dp/3.0_dp)

    b0 = 2.0_dp**(1.0_dp/3.0_dp) &
        * (rk - 4.0_dp * b2) * (4.0_dp * rk - b2) / denom
    end function airy_B0_turn


    pure elemental subroutine airy_zeta_q(K, beta, chi, turn_chi, zeta, q)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, chi, turn_chi
    real(dp), intent(out) :: zeta, q

    real(dp) :: s, action, delta, aturn, qprime_t

    s = curved_radius(K, chi)

    if (abs(s) <= CACHE_EPS) then
        q = huge(1.0_dp)
        zeta = huge(1.0_dp)
        return
    end if

    q = 1.0_dp / (s*s) - beta*beta
    delta = chi - turn_chi

    if (abs(delta) <= 1.0e-7_dp * max(1.0_dp, turn_chi)) then
        aturn = airy_turn_a(K, beta)
        qprime_t = -aturn**3
        zeta = -aturn * delta
        q = qprime_t * delta
        return
    end if

    action = qintegral_exact(s, beta, K)

    if (chi <= turn_chi) then
        zeta = (1.5_dp * max(action, 0.0_dp))**(2.0_dp/3.0_dp)
    else
        zeta = -(1.5_dp * max(action, 0.0_dp))**(2.0_dp/3.0_dp)
    end if
    end subroutine airy_zeta_q


    pure real(dp) function airy_liouville_amp(K, beta, zeta, q) result(amp)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, zeta, q

    real(dp) :: aturn

    if (abs(zeta) <= 1.0e-10_dp .or. abs(q) <= 1.0e-10_dp * max(1.0_dp, beta*beta)) then
        aturn = airy_turn_a(K, beta)
        amp = aturn**(-0.5_dp)
    else
        amp = abs(zeta / q)**0.25_dp
    end if
    end function airy_liouville_amp


    pure real(dp) function airy_turn_a(K, beta) result(a)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: st, ct

    st = 1.0_dp / beta
    ct = sqrt(max(0.0_dp, 1.0_dp - real(K, dp) * st*st))
    a = (2.0_dp * ct * beta**3)**(1.0_dp/3.0_dp)
    end function airy_turn_a


    function airy_u(l, K, nu, chi, ok, log_norm_in) result(u)
    ! Fast second-order Airy/Olver one-point approximation to reduced u=S_K phi.
    ! It builds only the local psi data needed for
    ! the requested point, sampling directly in chi between the turn and chi.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, chi
    logical, intent(out), optional :: ok
    real(dp), intent(in), optional :: log_norm_in
    real(dp) :: u

    real(dp) :: achi, symm
    logical :: lok

    call normalize_chi(l, K, nu, chi, achi, symm)
    u = symm * airy_u_normalized(l, K, nu, achi, lok, log_norm_in)
    if (present(ok)) ok = lok
    end function airy_u


    function airy_u_normalized(l, K, nu, achi, ok, log_norm_in) result(u)
    ! Second-order approximation for already-normalized achi; no parity sign.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, achi
    logical, intent(out), optional :: ok
    real(dp), intent(in), optional :: log_norm_in
    real(dp) :: u

    real(dp) :: lambda, beta, turn_chi, zeta, q
    real(dp) :: xairy, amp, ai, aip, b0, a1, series, log_norm
    logical :: lok

    u = 0.0_dp
    if (present(ok)) ok = .false.

    if (.not. airy_second_base_ok(l, K, nu)) return
    if (achi <= CACHE_EPS) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu / lambda
    turn_chi = turning_point(lambda, nu, K)

    call airy_zeta_q(K, beta, achi, turn_chi, zeta, q)
    xairy = lambda**(2.0_dp/3.0_dp) * zeta

    lok = .true.
    if (present(ok)) ok = lok
    if (xairy > AIRY_SECOND_TAIL_ZERO_X) return

    amp = airy_liouville_amp(K, beta, zeta, q)

    call second_coeffs_onepoint_fast(K, beta, turn_chi, achi, zeta, b0, a1)

    call airy_fast(xairy, ai, aip)

    series = ai * (1.0_dp + a1 / (lambda * lambda)) &
        + b0 * aip / lambda**(4.0_dp/3.0_dp)

    if (present(log_norm_in)) then
        log_norm = log_norm_in
    else
        call compute_airy_second_norm_fast(l, K, nu, log_norm)
    end if

    if (log_norm <= log(tiny(1.0_dp))) then
        u = 0.0_dp
    else if (log_norm >= log(huge(1.0_dp))) then
        u = sign(huge(1.0_dp), amp * series)
    else
        u = exp(log_norm) * amp * series
    end if
    end function airy_u_normalized


    pure subroutine second_coeffs_onepoint_fast(K, beta, turn_chi, chi, zeta, b0, a1)
    ! Build B0(zeta) and A1(zeta) for this one point only.
    !
    ! This version avoids both expensive zeta->chi inversion and the small
    ! dense Vandermonde solve.  It samples psi directly along a same-side chi
    ! segment from the turning point, interpolates in scaled zeta using Newton
    ! divided differences, and then evaluates the Olver B0/A1 functionals
    ! analytically.  Very small |zeta| uses a slightly enlarged same-side fit
    ! segment to avoid cancellation in airy_psi_from_s near the turn.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, turn_chi, chi, zeta
    real(dp), intent(out) :: b0, a1

    integer :: i
    real(dp) :: psi0, dchi, tau, chis, zs, qs, zscale, aturn
    real(dp) :: w(0:AIRY_SECOND_FAST_DEG)
    real(dp) :: f(0:AIRY_SECOND_FAST_DEG)
    real(dp) :: c(0:AIRY_SECOND_FAST_DEG)

    psi0 = airy_B0_turn(K, beta)

    if (abs(zeta) <= 100.0_dp * tiny(1.0_dp)) then
        b0 = psi0
        a1 = 0.0_dp
        return
    end if

    zscale = zeta
    dchi = chi - turn_chi

    if (abs(zeta) < AIRY_SECOND_FIT_ZETA_MIN) then
        zscale = sign(AIRY_SECOND_FIT_ZETA_MIN, zeta)
        aturn = airy_turn_a(K, beta)
        dchi = -zscale / max(aturn, CACHE_EPS)

        if (turn_chi + dchi <= CACHE_EPS) then
            dchi = max(chi - turn_chi, -0.5_dp * turn_chi)
        else if (K == 1 .and. turn_chi + dchi >= PI/2.0_dp - 10.0_dp*CACHE_EPS) then
            dchi = max(chi - turn_chi, PI/2.0_dp - 10.0_dp*CACHE_EPS - turn_chi)
        end if
    end if

    w(0) = 0.0_dp
    f(0) = psi0

    do i = 1, AIRY_SECOND_FAST_DEG
        ! Fixed Lobatto-like chi-segment nodes; no per-call cosine.
        tau = AIRY_SECOND_TAU(i)
        chis = turn_chi + tau * dchi

        call airy_zeta_q(K, beta, chis, turn_chi, zs, qs)
        w(i) = zs / zscale

        if (abs(zs) <= AIRY_PSI_ZETA_LOCAL) then
            f(i) = psi0
        else
            f(i) = airy_psi_from_chi(K, beta, zs, chis)
        end if
    end do

    call interp_power_from_nodes(w, f, c)
    call eval_second_scaled_poly(c, zscale, zeta, b0, a1)
    end subroutine second_coeffs_onepoint_fast


    pure subroutine interp_power_from_nodes(x, y, c)
    ! Convert interpolation data (x_i,y_i), i=0..n, to power coefficients
    ! c such that p(x)=sum_i c(i)*x**i.  Uses Newton divided differences;
    ! for n<=4 this is cheaper and better conditioned than a dense solve.
    real(dp), intent(in) :: x(0:AIRY_SECOND_FAST_DEG), y(0:AIRY_SECOND_FAST_DEG)
    real(dp), intent(out) :: c(0:AIRY_SECOND_FAST_DEG)

    integer :: i, j, k
    real(dp) :: dd(0:AIRY_SECOND_FAST_DEG), basis(0:AIRY_SECOND_FAST_DEG)
    real(dp) :: newbasis(0:AIRY_SECOND_FAST_DEG), den

    dd = y
    do j = 1, AIRY_SECOND_FAST_DEG
        do i = AIRY_SECOND_FAST_DEG, j, -1
            den = x(i) - x(i-j)
            if (abs(den) <= 100.0_dp * tiny(1.0_dp)) then
                dd(i) = 0.0_dp
            else
                dd(i) = (dd(i) - dd(i-1)) / den
            end if
        end do
    end do

    c = 0.0_dp
    basis = 0.0_dp
    basis(0) = 1.0_dp

    do k = 0, AIRY_SECOND_FAST_DEG
        do i = 0, k
            c(i) = c(i) + dd(k) * basis(i)
        end do

        if (k < AIRY_SECOND_FAST_DEG) then
            newbasis = 0.0_dp
            newbasis(0) = -x(k) * basis(0)
            do i = 1, k+1
                newbasis(i) = basis(i-1)
                if (i <= k) newbasis(i) = newbasis(i) - x(k) * basis(i)
            end do
            basis = newbasis
        end if
    end do
    end subroutine interp_power_from_nodes


    pure subroutine eval_second_scaled_poly(c, zscale, zeta, b0, a1)
    ! psi(v) = sum_i c(i) * (v/zscale)**i.
    !
    ! Then
    !   B0(zeta) = sum_i c(i)*(zeta/zscale)**i/(2*i+1)
    ! and A1 follows from A1' = 0.5*(psi*B0 - B0''), A1(0)=0.
    ! This keeps the evaluation scaled, avoiding large powers of 1/zeta.
    real(dp), intent(in) :: c(0:AIRY_SECOND_FAST_DEG), zscale, zeta
    real(dp), intent(out) :: b0, a1

    integer :: i, j
    real(dp) :: rho, prod_int, d2_int
    real(dp) :: rp(0:2*AIRY_SECOND_FAST_DEG)

    if (abs(zeta) <= 100.0_dp * tiny(1.0_dp)) then
        b0 = c(0)
        a1 = 0.0_dp
        return
    end if

    rho = zeta / zscale
    rp(0) = 1.0_dp
    do i = 1, 2*AIRY_SECOND_FAST_DEG
        rp(i) = rp(i-1) * rho
    end do

    b0 = 0.0_dp
    do i = 0, AIRY_SECOND_FAST_DEG
        b0 = b0 + c(i) * rp(i) / real(2*i + 1, dp)
    end do

    prod_int = 0.0_dp
    do i = 0, AIRY_SECOND_FAST_DEG
        do j = 0, AIRY_SECOND_FAST_DEG
            prod_int = prod_int + c(i) * c(j) * rp(i+j) / &
                ( real(2*j + 1, dp) * real(i + j + 1, dp) )
        end do
    end do

    d2_int = 0.0_dp
    do i = 2, AIRY_SECOND_FAST_DEG
        d2_int = d2_int + c(i) * real(i, dp) * rp(i) / real(2*i + 1, dp)
    end do

    a1 = 0.5_dp * ( zeta * prod_int - d2_int / zeta )
    end subroutine eval_second_scaled_poly


    pure elemental logical function airy_ok(l, K, nu, achi) result(ok)
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, achi

    ok = .false.
    if (.not. airy_second_base_ok(l, K, nu)) return
    if (achi <= CACHE_EPS) return

    ok = .true.
    end function airy_ok


    pure elemental logical function airy_second_base_ok(l, K, nu) result(ok)
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu

    real(dp) :: lambda, beta

    ok = .false.
    if (l < AIRY_SECOND_L_MIN) return
    if (nu <= 0.0_dp) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu / lambda

    select case (K)
    case (0)
        ok = .true.
    case (-1)
        ok = nu >= AIRY_SECOND_OPEN_NU_MIN
    case (1)
        if (beta <= 1.0_dp) return
        ok = lambda * (beta*beta - 1.0_dp) >= AIRY_SECOND_CLOSED_GATE
    case default
        ok = .false.
    end select
    end function airy_second_base_ok


    pure subroutine compute_airy_second_norm_fast(l, K, nu, log_norm)
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu
    real(dp), intent(out) :: log_norm

    real(dp) :: lambda, beta, q0, logu0, log_r, lk, rk

    log_norm = -huge(1.0_dp)
    if (.not. airy_second_base_ok(l, K, nu)) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu / lambda
    rk = real(K, dp)

    logu0 = log_origin_u0_fast(l, K, nu)
    if (logu0 <= -0.5_dp * huge(1.0_dp)) return

    q0 = airy_q0(K, beta)

    if (K == 0) then
        lk = 1.0_dp / 12.0_dp
    else
        lk = 1.0_dp / 12.0_dp - rk / (24.0_dp * (beta*beta - rk))
    end if

    ! Second-order scalar normalization.  The -1/(225 lambda^2) term
    ! corresponds to the coefficient convention A1(0)=0 used below.
    log_r = lk / lambda &
        - 1.0_dp / (225.0_dp * lambda*lambda) &
        - 1.0_dp / (360.0_dp * lambda**3) &
        + 1.0_dp / (1260.0_dp * lambda**5)

    log_norm = log(2.0_dp * sqrt(PI)) &
        + log(lambda) / 6.0_dp &
        + logu0 &
        + lambda * q0 &
        + log_r
    end subroutine compute_airy_second_norm_fast



    pure elemental real(dp) function airy_psi_from_chi(K, beta, zeta, chi) result(psi)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, zeta, chi

    real(dp) :: s

    if (abs(zeta) <= AIRY_PSI_ZETA_LOCAL) then
        psi = airy_B0_turn(K, beta)
        return
    end if

    s = curved_radius(K, chi)
    psi = airy_psi_from_s(K, beta, s, zeta)
    end function airy_psi_from_chi


    pure elemental real(dp) function airy_psi_from_s(K, beta, s, zeta) result(psi)
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, s, zeta

    real(dp) :: c, s2, q, qp, qpp, r

    if (abs(zeta) <= AIRY_PSI_ZETA_LOCAL) then
        psi = airy_B0_turn(K, beta)
        return
    end if

    s2 = max(s*s, tiny(1.0_dp))
    c = sqrt(max(0.0_dp, 1.0_dp - real(K, dp) * s2))

    q = 1.0_dp / s2 - beta*beta
    if (abs(q) <= 1.0e-8_dp * max(1.0_dp, beta*beta)) then
        psi = airy_B0_turn(K, beta)
        return
    end if

    r = -0.25_dp / s2
    qp = -2.0_dp * c / (s2 * sqrt(s2))
    qpp = 6.0_dp / (s2*s2) - 4.0_dp * real(K, dp) / s2

    ! Origin-side cancellation was checked; affected points are deep in the Airy tail.
    psi = zeta * r / q &
        + zeta * qpp / (4.0_dp * q*q) &
        - 5.0_dp * zeta * qp*qp / (16.0_dp * q*q*q) &
        + 5.0_dp / (16.0_dp * zeta*zeta)
    end function airy_psi_from_s


    ! Fast one-shot approximation to log product_{j=1}^l (nu^2+j^2), accurate at O(1e-6)
    pure elemental function log_prod_plus_one(l, nu) result(lp)
    integer, intent(in) :: l
    real(dp), intent(in) :: nu
    real(dp) :: lp
    real(dp) :: a, x, t, logr2, corr, den, phase
    real(dp) :: x2, y2, r2, inv, q, q2
    real(dp) :: invr, invi, qr, qi, pr, pi_, tr
    real(dp) :: m, sx, sy, s2, earg, ee
    real(dp), parameter :: HUGE_GUARD = sqrt(huge(1.0_dp))

    if (l <= 3) then
        lp = exact_small_l(l, abs(nu))
        return
    end if

    a = abs(nu)
    x = real(l + 1, dp)

    if (a <= 0.0_dp) then
        inv = 1.0_dp / x
        q = inv * inv
        corr = inv * (1.0_dp / 6.0_dp + q * (-1.0_dp / 180.0_dp))
        lp = (2.0_dp * x - 1.0_dp) * log(x) - 2.0_dp * x + LOG2PI + corr
        return
    end if

    if (a < HUGE_GUARD) then
        x2 = x * x
        y2 = a * a
        r2 = x2 + y2
        logr2 = log(r2)
        invr = x / r2
        invi = -a / r2
        qr = invr * invr - invi * invi
        qi = 2.0_dp * invr * invi
    else
        m = a
        sx = x / m
        sy = 1.0_dp
        s2 = sx * sx + sy * sy
        logr2 = 2.0_dp * log(m) + log(s2)
        invr = sx / (m * s2)
        invi = -sy / (m * s2)
        qr = invr * invr - invi * invi
        qi = 2.0_dp * invr * invi
    end if

    pr = invr
    pi_ = invi
    corr = (1.0_dp / 6.0_dp) * pr
    tr = pr * qr - pi_ * qi
    corr = corr - (1.0_dp / 180.0_dp) * tr

    t = PI * a

    if (t < 0.5_dp) then
        y2 = t * t
        den = y2 * (-1.0_dp / 6.0_dp + y2 * (1.0_dp / 180.0_dp + y2 * ( &
            -1.0_dp / 2835.0_dp + y2 * (1.0_dp / 37800.0_dp))))

        q = a / x
        q2 = q * q
        phase = -2.0_dp * a * q * (1.0_dp + q2 * (-1.0_dp / 3.0_dp + q2 / 5.0_dp))

        lp = (x - 0.5_dp) * logr2 + phase - 2.0_dp * x + LOG2PI + corr - den
    else
        earg = 2.0_dp * t
        if (earg < 36.0_dp) then
            ee = log1p_stable(-exp(-earg))
        else
            ee = 0.0_dp
        end if

        if (a < x) then
            phase = PI * a - 2.0_dp * a * atan2(a, x)
        else
            phase = 2.0_dp * a * atan2(x, a)
        end if

        lp = (x - 0.5_dp) * logr2 - 2.0_dp * x - log(a) + phase + ee + corr
    end if
    end function log_prod_plus_one


    pure elemental function exact_small_l(l, a) result(s)
    integer, intent(in) :: l
    real(dp), intent(in) :: a
    real(dp) :: s

    s = log_y2_plus_j2(a, 1.0_dp)
    if (l >= 2) s = s + log_y2_plus_j2(a, 2.0_dp)
    if (l >= 3) s = s + log_y2_plus_j2(a, 3.0_dp)
    end function exact_small_l


    pure elemental function log_y2_plus_j2(a, j) result(v)
    real(dp), intent(in) :: a, j
    real(dp) :: v, q

    if (a <= 0.0_dp) then
        v = 2.0_dp * log(j)
    else if (a > j) then
        q = j / a
        v = 2.0_dp * log(a) + log1p_stable(q * q)
    else
        q = a / j
        v = 2.0_dp * log(j) + log1p_stable(q * q)
    end if
    end function log_y2_plus_j2


    elemental real(dp) function log1p_stable(y) result(res)
    real(dp), intent(in) :: y

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: tiny_cut = 1.0e-4_dp

    if (abs(y) < tiny_cut) then
        res = y * (one + y * (-0.5_dp + y * (one/3.0_dp - 0.25_dp*y)))
    else
        res = log(one + y)
    end if
    end function log1p_stable


    end module HypersphericalBesselAiry
