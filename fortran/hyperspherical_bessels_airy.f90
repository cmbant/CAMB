    module HypersphericalBesselUtils
    use iso_fortran_env, only: real64
    use constants, only: const_pi
    implicit none
    private

    integer, parameter :: dp = real64
    real(dp), parameter :: CACHE_EPS = 1.0e-12_dp

    public :: dp, CACHE_EPS
    public :: normalize_chi, turning_point, curved_radius, qintegral_exact

    contains

    elemental subroutine normalize_chi(l, K, nu, chi, achi, symm)
    ! Fold a nonnegative radial coordinate chi into the fundamental domain.
    ! Open/flat retain chi unchanged; closed space additionally uses the pi-
    ! and pi/2-reflection symmetries of phi_l^nu to reduce achi to [0, pi/2]
    ! and returns the corresponding parity sign symm.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, chi
    real(dp), intent(out) :: achi, symm
    integer :: inu

    achi = abs(chi)
    symm = 1.0_dp

    if (K /= 1) return

    inu = nint(nu)
    achi = modulo(achi, 2.0_dp*const_pi)
    if (achi > const_pi) then
        achi = 2.0_dp*const_pi - achi
        if (mod(l, 2) /= 0) symm = -symm
    end if
    if (achi > const_pi/2.0_dp) then
        achi = const_pi - achi
        if (mod(inu - l - 1, 2) /= 0) symm = -symm
    end if

    end subroutine normalize_chi

    elemental real(dp) function turning_point(ell, nu, K)
    ! Classical turning point chi_t, where S_K(chi_t) = ell/nu and the
    ! centrifugal term ell^2/S_K^2 balances nu^2.  Callers must supply
    ! ell < nu for K=1.
    real(dp), intent(in) :: ell, nu
    integer, intent(in) :: K

    select case (K)
    case (-1)
        turning_point = asinh(ell/nu)
    case (0)
        turning_point = ell/nu
    case (1)
        turning_point = asin(ell/nu)
    case default
        turning_point = 0.0_dp
    end select

    end function turning_point

    elemental real(dp) function curved_radius(K, chi)
    ! Comoving angular-diameter distance S_K(chi): sinh, chi, sin for K=-1,0,1.
    integer, intent(in) :: K
    real(dp), intent(in) :: chi

    select case (K)
    case (-1)
        curved_radius = sinh(chi)
    case (1)
        curved_radius = sin(chi)
    case default
        curved_radius = chi
    end select

    end function curved_radius

    elemental real(dp) function qintegral_exact(sin_K, alpha, K) result(q)
    ! Closed-form Liouville-Green action |int_chi^chi_t sqrt(|1/S_K^2 - alpha^2|) dchi'|
    ! between chi (with sin_K = S_K(chi)) and the turning point, for alpha = nu/ell.
    ! The x = alpha*S_K > 1 branch is the oscillatory side, x < 1 the evanescent side.
    real(dp), intent(in) :: sin_K, alpha
    integer, intent(in) :: K

    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, half = 0.5_dp
    real(dp) :: x, x2, a2, ha, r1, r2, u, m

    x = alpha*sin_K
    x2 = x*x
    a2 = alpha*alpha
    ha = half*alpha

    select case (K)

    case (0)
        if (x > one) then
            r1 = sqrt(x2 - one)
            q = self_minus_atan(r1)
        else
            r1 = sqrt(max(zero, one - x2))
            q = atanh_minus_self(r1, x)
        end if

    case (-1)
        if (x > one) then
            r1 = sqrt(x2 - one)
            r2 = sqrt(x2 + a2)
            q = alpha*asinh(r1/sqrt(one + a2)) - atan2(alpha*r1, r2)
        else
            u = sqrt(max(zero, one - x2))
            r2 = sqrt(x2 + a2)
            r1 = u*r2
            q = ha*atan2(-two*r1, two*x2 + a2 - one) &
                + asinh(alpha*u/max(x*sqrt(one + a2), CACHE_EPS))
        end if

    case default
        if (x > one) then
            r1 = sqrt(x2 - one)
            r2 = sqrt(max(zero, a2 - x2))
            q = ha*atan2(two*r1*r2, a2 + one - two*x2) &
                - atan2(alpha*r1, r2)
        else
            m = a2 - one
            if (m > CACHE_EPS) then
                u = sqrt(max(zero, one - x2))
                q = asinh(alpha*u/max(x*sqrt(m), CACHE_EPS)) - alpha*asinh(u/sqrt(m))
            else
                q = -half*log(max(x2, CACHE_EPS))
            end if
        end if
    end select

    end function qintegral_exact

    elemental real(dp) function atanh_minus_self(y, x) result(res)
    ! atanh(y) - y for y = sqrt(1 - x^2) with x > 0.  Since (1 + y)(1-y) = x^2,
    ! atanh(y) = log((1 + y)/x); using that identity instead of atanh avoids the
    ! cancellation in 1 - y, which otherwise destroys the result as x -> 0
    ! (about 0.4 absolute error at x = 1e-8).  Small y still needs the series.
    real(dp), intent(in) :: y, x

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: small_cut = 1.0e-3_dp
    real(dp) :: y2

    if (abs(y) < small_cut) then
        y2 = y*y
        res = y*y2*(one/3.0_dp + y2*(one/5.0_dp + y2*(one/7.0_dp)))
    else
        res = log((one + y)/max(x, tiny(one))) - y
    end if

    end function atanh_minus_self

    elemental real(dp) function self_minus_atan(y) result(res)
    ! y - atan(y), series-summed for small |y| where the two terms cancel.
    real(dp), intent(in) :: y

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: small_cut = 1.0e-3_dp
    real(dp) :: y2

    if (abs(y) < small_cut) then
        y2 = y*y
        res = y*y2*(one/3.0_dp + y2*(-one/5.0_dp + y2*(one/7.0_dp)))
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
    use constants, only: const_pi
    use HypersphericalBesselUtils, only: dp, CACHE_EPS, curved_radius, qintegral_exact
    use MathUtils, only: airy_fast
    implicit none
    private

    ! Analytic second-order one-point Airy/Olver coefficients, with a short
    ! local series to avoid cancellation at the turning point.
    ! Keep the existing calibrated gates and peak-negligible tail cutoff.
    integer, parameter :: AIRY_SECOND_L_MIN = 10
    real(dp), parameter :: AIRY_SECOND_TAIL_ZERO_X = 25.77_dp
    real(dp), parameter :: AIRY_SECOND_OPEN_NU_MIN = 7.0_dp
    real(dp), parameter :: AIRY_SECOND_CLOSED_GATE = 15.0_dp
    real(dp), parameter :: AIRY_SECOND_LOCAL_Y = 0.05_dp

    public :: airy_u_normalized
    public :: airy_ok

    contains

    pure real(dp) function airy_b0_prime_turn(kappa) result(bp)
    ! Derivative with respect to zeta at the turning point (kappa=K/beta^2).
    real(dp), intent(in) :: kappa

    bp = (1.0_dp + kappa)*(8.0_dp*kappa*kappa - 9.0_dp*kappa + 8.0_dp) &
        /(900.0_dp*(1.0_dp - kappa)**2)

    end function airy_b0_prime_turn

    pure subroutine airy_zeta_ratio(K, beta, s, zeta, ratio, ct)
    ! Action and Liouville amplitude ratio, without constructing chi_t.
    ! Near the turn use the analytic series in y, also used by the B0
    ! expansion.  Form zeta/q without dividing two small quantities.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, s
    real(dp), intent(out) :: zeta, ratio, ct
    real(dp) :: s2, b2, kappa, c2, delta, y, w, p, r, h2, h3, h4, h5, f, zscale

    s2 = s*s
    b2 = beta*beta
    kappa = real(K, dp)/b2
    c2 = max(0.0_dp, 1.0_dp - real(K, dp)*s2)
    delta = 1.0_dp - b2*s2
    ct = 0.0_dp
    if (abs(delta)*max(1.0_dp, abs(kappa)) <= AIRY_SECOND_LOCAL_Y*c2) then
        y = delta/c2
        w = kappa*y
        p = y + w
        r = y*w
        ! h_j = sum_{i = 0}^j y^(j-i) w^i; F = sum 3 h_j/(2j + 3).
        h2 = p*p - r
        h3 = p*h2 - r*p
        h4 = p*h3 - r*h2
        h5 = p*h4 - r*h3
        f = 1.0_dp + 0.6_dp*p + (3.0_dp/7.0_dp)*h2 + h3/3.0_dp &
            + (3.0_dp/11.0_dp)*h4 + (3.0_dp/13.0_dp)*h5
        ct = ((1.0_dp - kappa)/2.0_dp)**(2.0_dp/3.0_dp)
        zscale = ct*f**(2.0_dp/3.0_dp)
        zeta = y*zscale
        ratio = zscale*s2/c2
    else
        zeta = sign((1.5_dp*max(0.0_dp, qintegral_exact(s, beta, K)))**(2.0_dp/3.0_dp), delta)
        ratio = abs(zeta*s2/delta)
    end if

    end subroutine airy_zeta_ratio

    function airy_u_normalized(l, K, nu, achi, ok, log_norm_in, sin_k_in) result(u)
    ! Fast second-order Airy/Olver one-point approximation to the reduced
    ! u = S_K phi.  B0 and A1 are obtained directly from the analytic action.
    !
    ! achi must already have been folded into the fundamental domain by
    ! normalize_chi; the caller applies the resulting parity sign itself.
    ! ok is false where the calibrated gate does not accept the point.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, achi
    logical, intent(out), optional :: ok
    real(dp), intent(in), optional :: log_norm_in, sin_k_in
    real(dp) :: u

    real(dp) :: lambda, beta, zeta, ratio, sin_k, ct, lc, lc2
    real(dp) :: xairy, ai, aip, b0, a1, bp0, series, norm, inv_l2

    u = 0.0_dp
    if (present(ok)) ok = .false.

    if (.not. airy_ok(l, K, nu, achi)) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu/lambda
    if (present(sin_k_in)) then
        sin_k = sin_k_in
    else
        sin_k = curved_radius(K, achi)
    end if
    if (present(ok)) ok = .true.
    if (abs(sin_k) <= CACHE_EPS) return

    call airy_zeta_ratio(K, beta, sin_k, zeta, ratio, ct)
    lc = lambda**(1.0_dp/3.0_dp)
    lc2 = lc*lc
    xairy = lc2*zeta
    if (xairy > AIRY_SECOND_TAIL_ZERO_X) return

    call second_coeffs_onepoint_fast(K, beta, sin_k, zeta, ct, b0, a1, bp0)
    call airy_fast(xairy, ai, aip)

    inv_l2 = 1.0_dp/(lambda*lambda)
    series = ai*(1.0_dp + a1*inv_l2) + b0*aip*(lc2*inv_l2)

    if (present(log_norm_in)) then
        norm = exp(log_norm_in)
    else
        norm = airy_second_norm(K, nu, lambda, lc, bp0)
    end if
    u = norm*sqrt(sqrt(ratio))*series

    end function airy_u_normalized

    pure subroutine second_coeffs_onepoint_fast(K, beta, s, zeta, ct, b0, a1, bp0)
    ! B0+2*zeta*B0'=psi and 2*A1'=psi*B0-B0'' give
    ! A1 = (zeta*B0**2-B0'+B0'(0))/2.  The B0 integral has an
    ! elementary primitive for q=1/S_K**2-beta**2 and r=-1/(4*S_K**2).
    ! Use the real formula on both sides of the turn.  A quadratic
    ! expansion in y is sufficient for the peak-relative accuracy target
    ! near the turn, where the separate closed-form terms cancel.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, s, zeta, ct
    real(dp), intent(out) :: b0, a1, bp0

    real(dp) :: kappa, k2, s2, c2, delta, b2, scale, y, p0, p1, p2
    real(dp) :: dydz, bp, sig, t, zinv

    b2 = beta*beta
    kappa = real(K, dp)/b2
    k2 = kappa*kappa
    s2 = s*s
    c2 = max(0.0_dp, 1.0_dp - real(K, dp)*s2)
    delta = 1.0_dp - b2*s2
    bp0 = airy_b0_prime_turn(kappa)

    if (ct > 0.0_dp) then
        y = delta/c2
        scale = 1.0_dp/(8.0_dp*(1.0_dp - kappa)*sqrt(ct))
        p0 = (kappa - 4.0_dp)*(4.0_dp*kappa - 1.0_dp)/35.0_dp
        p1 = (1.0_dp + kappa)*(8.0_dp*k2 - 9.0_dp*kappa + 8.0_dp)/225.0_dp
        p2 = (11108.0_dp*(k2*k2 + 1.0_dp) - 2543.0_dp*kappa*(k2 + 1.0_dp) &
            + 823.0_dp*k2)/606375.0_dp

        b0 = scale*(p0 + y*(p1 + y*p2))
        ! dy/dzeta from the action identity; at the turn use its finite limit.
        if (abs(y)*max(1.0_dp, abs(kappa)) <= 1.0e-7_dp) then
            dydz = 1.0_dp/ct
        else
            dydz = 2.0_dp*(1.0_dp - y)*(1.0_dp - kappa*y) &
                *sqrt(abs(zeta/y))/(1.0_dp - kappa)
        end if
        bp = scale*(p1 + 2.0_dp*y*p2)*dydz
    else
        sig = sign(1.0_dp, zeta)
        delta = abs(delta)
        t = sqrt(c2/delta)
        zinv = 1.0_dp/zeta
        b0 = ((5.0_dp/3.0_dp)*t*t - sig*(1.0_dp + kappa))*t &
            /(8.0_dp*(1.0_dp - kappa)*sqrt(abs(zeta))) - (5.0_dp/48.0_dp)*zinv*zinv
        bp = -sig*(5.0_dp*t*t - sig*(1.0_dp + kappa))*b2*s2/(8.0_dp*delta*delta) &
            - 0.5_dp*b0*zinv + (5.0_dp/32.0_dp)*zinv**3
    end if

    a1 = 0.5_dp*(zeta*b0*b0 - bp + bp0)

    end subroutine second_coeffs_onepoint_fast

    pure elemental logical function airy_ok(l, K, nu, achi) result(ok)
    ! Public predicate: would airy_u_normalized return a usable value here?
    ! It is the sole entry guard of airy_u_normalized, so callers can test
    ! acceptance without evaluating.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu, achi

    ok = .false.
    if (.not. airy_second_base_ok(l, K, nu)) return
    if (achi <= CACHE_EPS) return

    ok = .true.

    end function airy_ok

    pure elemental logical function airy_second_base_ok(l, K, nu) result(ok)
    ! Calibrated chi-independent validity gate: l above AIRY_SECOND_L_MIN, and
    ! for K=-1 nu not too small, for K=1 the turning point safely inside the
    ! sphere (lambda*(beta^2-1) large enough that the Airy region is resolved).
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu

    real(dp) :: lambda, beta

    ok = .false.
    if (l < AIRY_SECOND_L_MIN) return
    if (nu <= 0.0_dp) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu/lambda

    select case (K)
    case (0)
        ok = .true.
    case (-1)
        ok = nu >= AIRY_SECOND_OPEN_NU_MIN
    case (1)
        if (beta <= 1.0_dp) return
        ok = lambda*(beta*beta - 1.0_dp) >= AIRY_SECOND_CLOSED_GATE
    case default
        ok = .false.
    end select

    end function airy_second_base_ok

    pure real(dp) function airy_second_norm(K, nu, lambda, lc, bp0) result(norm)
    ! Combine the origin gamma ratios, regularized action and Stirling
    ! corrections analytically.  For K=1 keep the cubic shifted-Stirling
    ! remainder (nu-lambda can be small); for K=-1 its effect is below
    ! 1e-6 even at the base gate l>=10, nu>=7.
    ! The log correction is small enough to linearize at the same tolerance.
    integer, intent(in) :: K
    real(dp), intent(in) :: nu, lambda, lc, bp0
    real(dp) :: corr, plus_inv, minus_inv

    corr = -0.5_dp*bp0/(lambda*lambda)
    if (K == 1) then
        plus_inv = 1.0_dp/(nu + lambda)
        minus_inv = 1.0_dp/(nu - lambda)
        corr = corr + (7.0_dp/5760.0_dp)*(plus_inv**3 - minus_inv**3)
    end if
    norm = sqrt(const_pi/nu)/lc*(1.0_dp + corr)

    end function airy_second_norm

    end module HypersphericalBesselAiry
