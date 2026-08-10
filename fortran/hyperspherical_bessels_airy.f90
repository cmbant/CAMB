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
    use HypersphericalBesselUtils, only: dp, CACHE_EPS, curved_radius, qintegral_exact, turning_point
    use MathUtils, only: airy_fast
    implicit none
    private

    real(dp), parameter :: LOG2PI = log(2.0_dp*const_pi)

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
        1.0000000000000000e0_dp]
    integer, parameter :: AIRY_SECOND_L_MIN = 10
    real(dp), parameter :: AIRY_SECOND_TAIL_ZERO_X = 25.77_dp
    real(dp), parameter :: AIRY_SECOND_OPEN_NU_MIN = 7.0_dp
    real(dp), parameter :: AIRY_SECOND_CLOSED_GATE = 15.0_dp
    real(dp), parameter :: AIRY_PSI_ZETA_LOCAL = 2.0e-3_dp
    real(dp), parameter :: AIRY_SECOND_FIT_ZETA_MIN = 1.6e-2_dp

    ! 1/(2i+1) and 1/(i+1), the fixed rational weights of the B0 and A1
    ! quadratures below, tabulated to keep divides out of the inner loops.
    real(dp), parameter :: INV_ODD(0:AIRY_SECOND_FAST_DEG) = &
        [1.0_dp, 1.0_dp/3.0_dp, 1.0_dp/5.0_dp, 1.0_dp/7.0_dp, 1.0_dp/9.0_dp]
    real(dp), parameter :: INV_NP1(0:2*AIRY_SECOND_FAST_DEG) = &
        [1.0_dp, 1.0_dp/2.0_dp, 1.0_dp/3.0_dp, 1.0_dp/4.0_dp, 1.0_dp/5.0_dp, &
        1.0_dp/6.0_dp, 1.0_dp/7.0_dp, 1.0_dp/8.0_dp, 1.0_dp/9.0_dp]

    public :: airy_u_normalized
    public :: airy_ok

    contains

    pure real(dp) function log_origin_u0_fast(l, K, nu) result(logu0)
    ! log of the exact small-chi coefficient u0 in u = S_K phi ~ u0 chi^(l+1),
    ! i.e. log[ sqrt(prod_{j=1}^l (nu^2 - K j^2)) / (2l+1)!! ], written with
    ! log_gamma so it stays finite for large l.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu

    real(dp) :: logdf, logprod

    logdf = real(l + 1, dp)*log(2.0_dp) &
        + log_gamma(real(l, dp) + 1.5_dp) &
        - 0.5_dp*log(const_pi)

    select case (K)
    case (0)
        logprod = real(l, dp)*log(nu)
    case (1)
        logprod = 0.5_dp*(log_gamma(nu + l + 1) - log_gamma(nu - l) - log(nu))
    case (-1)
        logprod = 0.5_dp*log_prod_plus_one(l, nu)
    case default
        logprod = -huge(1.0_dp)
    end select

    logu0 = logprod - logdf

    end function log_origin_u0_fast

    pure real(dp) function airy_q0(K, beta) result(q0)
    ! Regularized origin action q0 = lim_{chi->0} [ qintegral_exact(S_K,beta,K) + log S_K ],
    ! used to tie the Airy form to the exact origin amplitude in
    ! compute_airy_second_norm_fast.  Branches are algebraically equal and split
    ! only for stability at beta -> 1+ (K=1) and at large beta.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: x, invb

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
                    + 0.5_dp*(x*log(x) &
                    - (beta + 1.0_dp)*log(beta + 1.0_dp))
            end if
        else
            invb = 1.0_dp/beta
            q0 = log(2.0_dp) - log(beta) &
                - 0.5_dp*log(1.0_dp - invb*invb) &
                - beta*atanh(invb)
        end if

    case (-1)
        if (beta < 1.0_dp) then
            q0 = log(2.0_dp) &
                - 0.5_dp*log(1.0_dp + beta*beta) &
                - beta*atan2(1.0_dp, beta)
        else
            invb = 1.0_dp/beta
            q0 = log(2.0_dp) - log(beta) &
                - 0.5_dp*log(1.0_dp + invb*invb) &
                - beta*atan(invb)
        end if
    case default
        q0 = log(2.0_dp/beta) - 1.0_dp
    end select

    end function airy_q0

    pure real(dp) function airy_B0_turn(K, beta) result(b0)
    ! Value of Olver's psi at the turning point, psi(0) = B0(0), which is the
    ! finite limit of airy_psi_from_s as zeta -> 0.  Used both as the zeta=0
    ! interpolation node and as the near-turn value of psi.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: rk, b2, denom

    rk = real(K, dp)
    b2 = beta*beta

    denom = 280.0_dp*beta**(4.0_dp/3.0_dp) &
        *max(b2 - rk, tiny(1.0_dp))**(4.0_dp/3.0_dp)

    b0 = 2.0_dp**(1.0_dp/3.0_dp) &
        *(rk - 4.0_dp*b2)*(4.0_dp*rk - b2)/denom

    end function airy_B0_turn

    pure subroutine airy_zeta_q(K, beta, chi, turn_chi, zeta, q, sin_k)
    ! Liouville variable zeta (positive below the turning point, negative above)
    ! from zeta*(dzeta/dchi)^2 = q, together with q = 1/S_K^2 - beta^2.
    ! Within 1e-7 of the turning point both are replaced by their exact linear
    ! limits, zeta = -a*(chi-chi_t) and q = -a^3*(chi-chi_t) with a = airy_turn_a.
    ! Optionally returns S_K(chi) so callers need not recompute it.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, chi, turn_chi
    real(dp), intent(out) :: zeta, q
    real(dp), intent(out), optional :: sin_k

    real(dp) :: s, action, delta, aturn, qprime_t

    s = curved_radius(K, chi)
    if (present(sin_k)) sin_k = s

    if (abs(s) <= CACHE_EPS) then
        q = huge(1.0_dp)
        zeta = huge(1.0_dp)
        return
    end if

    q = 1.0_dp/(s*s) - beta*beta
    delta = chi - turn_chi

    if (abs(delta) <= 1.0e-7_dp*max(1.0_dp, turn_chi)) then
        aturn = airy_turn_a(K, beta)
        qprime_t = -aturn**3
        zeta = -aturn*delta
        q = qprime_t*delta
        return
    end if

    action = qintegral_exact(s, beta, K)

    if (chi <= turn_chi) then
        zeta = (1.5_dp*max(action, 0.0_dp))**(2.0_dp/3.0_dp)
    else
        zeta = -(1.5_dp*max(action, 0.0_dp))**(2.0_dp/3.0_dp)
    end if

    end subroutine airy_zeta_q

    pure real(dp) function airy_liouville_amp(K, beta, zeta, q) result(amp)
    ! Liouville-Green amplitude |zeta/q|^(1/4), with the removable 0/0 at the
    ! turning point replaced by its limit a^(-1/2).
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, zeta, q

    real(dp) :: aturn

    if (abs(zeta) <= 1.0e-10_dp .or. abs(q) <= 1.0e-10_dp*max(1.0_dp, beta*beta)) then
        aturn = airy_turn_a(K, beta)
        amp = aturn**(-0.5_dp)
    else
        amp = abs(zeta/q)**0.25_dp
    end if

    end function airy_liouville_amp

    pure real(dp) function airy_turn_a(K, beta) result(a)
    ! Turning-point scale a = (-dq/dchi)^(1/3) = (2 cos_K(chi_t) beta^3)^(1/3),
    ! the local slope linking zeta and chi near chi_t.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta

    real(dp) :: st, ct

    st = 1.0_dp/beta
    ct = sqrt(max(0.0_dp, 1.0_dp - real(K, dp)*st*st))
    a = (2.0_dp*ct*beta**3)**(1.0_dp/3.0_dp)

    end function airy_turn_a

    function airy_u_normalized(l, K, nu, achi, ok, log_norm_in) result(u)
    ! Fast second-order Airy/Olver one-point approximation to the reduced
    ! u = S_K phi.  Only the local psi data for the requested point is built,
    ! sampling directly in chi between the turning point and achi.
    !
    ! achi must already have been folded into the fundamental domain by
    ! normalize_chi; the caller applies the resulting parity sign itself.
    ! ok is false where the calibrated gate does not accept the point.
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

    if (.not. airy_ok(l, K, nu, achi)) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu/lambda
    turn_chi = turning_point(lambda, nu, K)

    call airy_zeta_q(K, beta, achi, turn_chi, zeta, q)
    xairy = lambda**(2.0_dp/3.0_dp)*zeta

    lok = .true.
    if (present(ok)) ok = lok
    if (xairy > AIRY_SECOND_TAIL_ZERO_X) return

    amp = airy_liouville_amp(K, beta, zeta, q)

    call second_coeffs_onepoint_fast(K, beta, turn_chi, achi, zeta, b0, a1)

    call airy_fast(xairy, ai, aip)

    series = ai*(1.0_dp + a1/(lambda*lambda)) &
        + b0*aip/lambda**(4.0_dp/3.0_dp)

    if (present(log_norm_in)) then
        log_norm = log_norm_in
    else
        call compute_airy_second_norm_fast(l, K, nu, log_norm)
    end if

    if (log_norm <= log(tiny(1.0_dp))) then
        u = 0.0_dp
    else if (log_norm >= log(huge(1.0_dp))) then
        u = sign(huge(1.0_dp), amp*series)
    else
        u = exp(log_norm)*amp*series
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
    logical :: rescaled
    real(dp) :: psi0, dchi, tau, chis, zs, qs, ss, zscale, aturn
    real(dp) :: w(0:AIRY_SECOND_FAST_DEG)
    real(dp) :: f(0:AIRY_SECOND_FAST_DEG)
    real(dp) :: c(0:AIRY_SECOND_FAST_DEG)

    psi0 = airy_B0_turn(K, beta)

    if (abs(zeta) <= 100.0_dp*tiny(1.0_dp)) then
        b0 = psi0
        a1 = 0.0_dp
        return
    end if

    zscale = zeta
    dchi = chi - turn_chi
    rescaled = abs(zeta) < AIRY_SECOND_FIT_ZETA_MIN

    if (rescaled) then
        zscale = sign(AIRY_SECOND_FIT_ZETA_MIN, zeta)
        aturn = airy_turn_a(K, beta)
        dchi = -zscale/max(aturn, CACHE_EPS)

        if (turn_chi + dchi <= CACHE_EPS) then
            dchi = max(chi - turn_chi, -0.5_dp*turn_chi)
        else if (K == 1 .and. turn_chi + dchi >= const_pi/2.0_dp - 10.0_dp*CACHE_EPS) then
            dchi = max(chi - turn_chi, const_pi/2.0_dp - 10.0_dp*CACHE_EPS - turn_chi)
        end if
    end if

    w(0) = 0.0_dp
    f(0) = psi0

    do i = 1, AIRY_SECOND_FAST_DEG
        ! Fixed Lobatto-like chi-segment nodes; no per-call cosine.
        tau = AIRY_SECOND_TAU(i)
        chis = turn_chi + tau*dchi

        call airy_zeta_q(K, beta, chis, turn_chi, zs, qs, ss)
        w(i) = zs/zscale
        ! airy_psi_from_s falls back to psi0 itself for small |zs|.
        f(i) = airy_psi_from_s(K, beta, ss, zs)
    end do

    call interp_power_from_nodes(w, f, c)
    call eval_second_scaled_poly(c, zscale, zeta, rescaled, b0, a1)

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
            den = x(i) - x(i - j)
            if (abs(den) <= 100.0_dp*tiny(1.0_dp)) then
                dd(i) = 0.0_dp
            else
                dd(i) = (dd(i) - dd(i - 1))/den
            end if
        end do
    end do

    c = 0.0_dp
    basis = 0.0_dp
    basis(0) = 1.0_dp

    do k = 0, AIRY_SECOND_FAST_DEG
        do i = 0, k
            c(i) = c(i) + dd(k)*basis(i)
        end do

        if (k < AIRY_SECOND_FAST_DEG) then
            newbasis = 0.0_dp
            newbasis(0) = -x(k)*basis(0)
            do i = 1, k + 1
                newbasis(i) = basis(i - 1)
                if (i <= k) newbasis(i) = newbasis(i) - x(k)*basis(i)
            end do
            basis = newbasis
        end if
    end do

    end subroutine interp_power_from_nodes

    pure subroutine eval_second_scaled_poly(c, zscale, zeta, rescaled, b0, a1)
    ! psi(v) = sum_i c(i) * (v/zscale)**i.
    !
    ! Then
    !   B0(zeta) = sum_i c(i)*(zeta/zscale)**i/(2*i + 1)
    ! and A1 follows from A1' = 0.5*(psi*B0 - B0''), A1(0)=0.
    ! This keeps the evaluation scaled, avoiding large powers of 1/zeta.
    !
    ! rescaled is false when zscale is just zeta, in which case every power of
    ! zeta/zscale is one and the power ladder is skipped; the caller knows this
    ! exactly, so it is passed in rather than tested for.
    real(dp), intent(in) :: c(0:AIRY_SECOND_FAST_DEG), zscale, zeta
    logical, intent(in) :: rescaled
    real(dp), intent(out) :: b0, a1

    integer :: i, j
    real(dp) :: rho, prod_int, d2_int, inner
    real(dp) :: rp(0:2*AIRY_SECOND_FAST_DEG), g(0:2*AIRY_SECOND_FAST_DEG)
    real(dp) :: bw(0:AIRY_SECOND_FAST_DEG), cb(0:AIRY_SECOND_FAST_DEG)

    if (abs(zeta) <= 100.0_dp*tiny(1.0_dp)) then
        b0 = c(0)
        a1 = 0.0_dp
        return
    end if

    cb = c*INV_ODD

    if (rescaled) then
        rho = zeta/zscale
        rp(0) = 1.0_dp
        do i = 1, 2*AIRY_SECOND_FAST_DEG
            rp(i) = rp(i - 1)*rho
        end do
        bw = cb*rp(0:AIRY_SECOND_FAST_DEG)
        g = rp*INV_NP1
    else
        bw = cb
        g = INV_NP1
    end if

    b0 = sum(bw)

    prod_int = 0.0_dp
    do i = 0, AIRY_SECOND_FAST_DEG
        inner = 0.0_dp
        do j = 0, AIRY_SECOND_FAST_DEG
            inner = inner + cb(j)*g(i + j)
        end do
        prod_int = prod_int + c(i)*inner
    end do

    d2_int = 0.0_dp
    do i = 2, AIRY_SECOND_FAST_DEG
        d2_int = d2_int + real(i, dp)*bw(i)
    end do

    a1 = 0.5_dp*(zeta*prod_int - d2_int/zeta)

    end subroutine eval_second_scaled_poly

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

    pure subroutine compute_airy_second_norm_fast(l, K, nu, log_norm)
    ! log of the single scalar that fixes the overall Airy normalization, by
    ! matching the chi -> 0 limit of amp*Ai(lambda^(2/3) zeta) to the exact
    ! origin behaviour u ~ u0 chi^(l+1).  Returned as a log because both
    ! factors overflow badly at large l.
    integer, intent(in) :: l, K
    real(dp), intent(in) :: nu
    real(dp), intent(out) :: log_norm

    real(dp) :: lambda, beta, q0, logu0, log_r, lk, rk

    log_norm = -huge(1.0_dp)
    if (.not. airy_second_base_ok(l, K, nu)) return

    lambda = real(l, dp) + 0.5_dp
    beta = nu/lambda
    rk = real(K, dp)

    logu0 = log_origin_u0_fast(l, K, nu)
    if (logu0 <= -0.5_dp*huge(1.0_dp)) return

    q0 = airy_q0(K, beta)

    if (K == 0) then
        lk = 1.0_dp/12.0_dp
    else
        lk = 1.0_dp/12.0_dp - rk/(24.0_dp*(beta*beta - rk))
    end if

    ! Second-order scalar normalization.  The -1/(225 lambda^2) term
    ! corresponds to the coefficient convention A1(0)=0 used below.
    log_r = lk/lambda &
        - 1.0_dp/(225.0_dp*lambda*lambda) &
        - 1.0_dp/(360.0_dp*lambda**3) &
        + 1.0_dp/(1260.0_dp*lambda**5)

    log_norm = log(2.0_dp*sqrt(const_pi)) &
        + log(lambda)/6.0_dp &
        + logu0 &
        + lambda*q0 &
        + log_r

    end subroutine compute_airy_second_norm_fast

    pure elemental real(dp) function airy_psi_from_s(K, beta, s, zeta) result(psi)
    ! Olver's error-control function
    !   psi = zeta r/q + zeta(4 q q'' - 5 q'^2)/(16 q^3) + 5/(16 zeta^2),
    ! for u'' = (lambda^2 q + r) u with q = 1/s^2 - beta^2 and r = -1/(4 s^2),
    ! primes being d/dchi.  Individually singular at the turning point, so
    ! near zeta = 0 the finite limit airy_B0_turn is returned instead.
    integer, intent(in) :: K
    real(dp), intent(in) :: beta, s, zeta

    real(dp) :: c, s2, q, qp, qpp, r

    if (abs(zeta) <= AIRY_PSI_ZETA_LOCAL) then
        psi = airy_B0_turn(K, beta)
        return
    end if

    s2 = max(s*s, tiny(1.0_dp))
    c = sqrt(max(0.0_dp, 1.0_dp - real(K, dp)*s2))

    q = 1.0_dp/s2 - beta*beta
    if (abs(q) <= 1.0e-8_dp*max(1.0_dp, beta*beta)) then
        psi = airy_B0_turn(K, beta)
        return
    end if

    r = -0.25_dp/s2
    qp = -2.0_dp*c/(s2*sqrt(s2))
    qpp = 6.0_dp/(s2*s2) - 4.0_dp*real(K, dp)/s2

    ! Origin-side cancellation was checked; affected points are deep in the Airy tail.
    psi = zeta*r/q &
        + zeta*qpp/(4.0_dp*q*q) &
        - 5.0_dp*zeta*qp*qp/(16.0_dp*q*q*q) &
        + 5.0_dp/(16.0_dp*zeta*zeta)

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
        inv = 1.0_dp/x
        q = inv*inv
        corr = inv*(1.0_dp/6.0_dp + q*(-1.0_dp/180.0_dp))
        lp = (2.0_dp*x - 1.0_dp)*log(x) - 2.0_dp*x + LOG2PI + corr
        return
    end if

    if (a < HUGE_GUARD) then
        x2 = x*x
        y2 = a*a
        r2 = x2 + y2
        logr2 = log(r2)
        invr = x/r2
        invi = -a/r2
        qr = invr*invr - invi*invi
        qi = 2.0_dp*invr*invi
    else
        m = a
        sx = x/m
        sy = 1.0_dp
        s2 = sx*sx + sy*sy
        logr2 = 2.0_dp*log(m) + log(s2)
        invr = sx/(m*s2)
        invi = -sy/(m*s2)
        qr = invr*invr - invi*invi
        qi = 2.0_dp*invr*invi
    end if

    pr = invr
    pi_ = invi
    corr = (1.0_dp/6.0_dp)*pr
    tr = pr*qr - pi_*qi
    corr = corr - (1.0_dp/180.0_dp)*tr

    t = const_pi*a

    if (t < 0.5_dp) then
        y2 = t*t
        den = y2*(-1.0_dp/6.0_dp + y2*(1.0_dp/180.0_dp + y2*( &
            -1.0_dp/2835.0_dp + y2*(1.0_dp/37800.0_dp))))

        q = a/x
        q2 = q*q
        phase = -2.0_dp*a*q*(1.0_dp + q2*(-1.0_dp/3.0_dp + q2/5.0_dp))

        lp = (x - 0.5_dp)*logr2 + phase - 2.0_dp*x + LOG2PI + corr - den
    else
        earg = 2.0_dp*t
        if (earg < 36.0_dp) then
            ee = log1p_stable(-exp(-earg))
        else
            ee = 0.0_dp
        end if

        ! pi*a - 2*a*atan2(a,x) written using atan2(a,x) + atan2(x,a) = pi/2,
        ! which has no cancellation when a >> x.
        phase = 2.0_dp*a*atan2(x, a)

        lp = (x - 0.5_dp)*logr2 - 2.0_dp*x - log(a) + phase + ee + corr
    end if

    end function log_prod_plus_one

    pure elemental function exact_small_l(l, a) result(s)
    ! Direct term-by-term sum of log(a^2+j^2) for the l <= 3 cases, where the
    ! Stirling form in log_prod_plus_one is not accurate enough.
    ! Assumes l >= 1: the empty l = 0 product would have to return 0, but every
    ! caller is gated on l >= AIRY_SECOND_L_MIN, so it is never requested.
    integer, intent(in) :: l
    real(dp), intent(in) :: a
    real(dp) :: s

    s = log_y2_plus_j2(a, 1.0_dp)
    if (l >= 2) s = s + log_y2_plus_j2(a, 2.0_dp)
    if (l >= 3) s = s + log_y2_plus_j2(a, 3.0_dp)

    end function exact_small_l

    pure elemental function log_y2_plus_j2(a, j) result(v)
    ! log(a^2 + j^2), factoring out the larger square so neither term overflows.
    real(dp), intent(in) :: a, j
    real(dp) :: v, q

    if (a <= 0.0_dp) then
        v = 2.0_dp*log(j)
    else if (a > j) then
        q = j/a
        v = 2.0_dp*log(a) + log1p_stable(q*q)
    else
        q = a/j
        v = 2.0_dp*log(j) + log1p_stable(q*q)
    end if

    end function log_y2_plus_j2

    elemental real(dp) function log1p_stable(y) result(res)
    ! log(1+y), series-summed for small |y| (Fortran has no intrinsic log1p).
    real(dp), intent(in) :: y

    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: tiny_cut = 1.0e-4_dp

    if (abs(y) < tiny_cut) then
        res = y*(one + y*(-0.5_dp + y*(one/3.0_dp - 0.25_dp*y)))
    else
        res = log(one + y)
    end if

    end function log1p_stable

    end module HypersphericalBesselAiry
