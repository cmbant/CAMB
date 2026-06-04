    module OlverHypersphericalBessel
    use Precision
    use MpiUtils
    use FlatBessels, only: bjl
    use SpherBessels, only: phi_recurs_stable
    implicit none
    private

    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: CACHE_EPS = 1.0e-12_dl
    ! Pointwise raw-Olver gate calibrated to keep peak-normalized errors below
    ! about 1e-4 on the open/closed validation grid.
    real(dl), parameter :: OLVER_GATE_ALPHA = 3._dl
    real(dl), parameter :: OLVER_GATE_OPEN_EPS = 2.6e-2_dl
    real(dl), parameter :: OLVER_GATE_CLOSED_EPS = 6.2e-3_dl
    real(dl), parameter :: SMALLCHI_GATE_ALPHA = 3._dl
    real(dl), parameter :: SMALLCHI_GATE_METRIC = 5.0e-2_dl

    public :: phi_olver, u_olver, phi_olver_smallchi
    public :: olver_cached_coordinate, compute_olver_z_amp_smallchi

    contains

    pure function olver_cached_coordinate(l, K, beta, chi) result(z)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: z, achi, symm

    call normalize_chi(l, K, beta, chi, achi, symm)
    call compute_olver_z_amp(l, K, beta, achi, z)
    end function olver_cached_coordinate


    function phi_olver(l, K, beta, chi) result(phi)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: phi

    phi = olver_value(l, K, beta, chi, reduced=.false., raw=.false.)
    end function phi_olver

    function u_olver(l, K, beta, chi) result(u)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: u

    u = olver_value(l, K, beta, chi, reduced=.true., raw=.false.)
    end function u_olver


    function olver_value(l, K, beta, chi, reduced, raw) result(val)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    logical, intent(in) :: reduced, raw
    real(dl) :: val
    real(dl) :: achi, symm, j_l

    call normalize_chi(l, K, beta, chi, achi, symm)

    if (K == 0) then
        call bjl(l, beta * achi, j_l)
        if (reduced) then
            val = achi * j_l
        else
            val = j_l
        end if
        return
    end if

    if (l <= 2) then
        val = symm * phi_recurs_stable(l, K, beta, achi)
        if (reduced) val = val * curved_radius(K, achi)
        return
    end if

    if (achi <= CACHE_EPS) then
        val = 0._dl
        return
    end if

    val = olver_reduced(l, K, beta, achi, symm, raw)
    if (.not. reduced) val = val / curved_radius(K, achi)
    end function olver_value


    function olver_reduced(l, K, beta, achi, symm, raw) result(u)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: achi, symm
    logical, intent(in) :: raw
    real(dl) :: u, beta
    real(dl) :: alpha_gate, metric, denom, z, amp, j_l

    alpha_gate = beta / real(l, dl)

    if (.not. raw) then
        if (use_smallchi_map(l, beta, achi, alpha_gate)) then
            u = olver_smallchi_reduced(l, K, beta, achi, symm)
            return
        end if

        ! Empirical stable-recursion fallback for the low-alpha corner of the
        ! full Olver map. With alpha ~= beta/l, the raw
        ! curved Olver remainder falls rapidly once alpha is moderately large
        ! (consistent with the expected curvature expansion order O(alpha^-4)).
        ! For alpha below the cutoff the peak-normalized error is governed more
        ! directly by the small endpoint parameter
        !
        !   open:   eps_peak = O(chi / (2 beta))
        !   closed: eps_peak = O(chi / (2 (beta-l))).
        !
        ! The constants below are grid-calibrated thresholds on these metrics.
        if (alpha_gate < OLVER_GATE_ALPHA) then
            if (K == 1) then
                metric = achi / (2._dl * (beta - real(l, dl)))
                if (metric > OLVER_GATE_CLOSED_EPS) then
                    u = symm * phi_recurs_stable(l, K, beta, achi) * curved_radius(K, achi)
                    return
                end if
            else if (K == -1) then
                metric = achi / (2._dl * max(beta, tiny(1._dl)))
                if (metric > OLVER_GATE_OPEN_EPS) then
                    u = symm * phi_recurs_stable(l, K, beta, achi) * curved_radius(K, achi)
                    return
                end if
            end if
        end if
    end if

    call compute_olver_z_amp(l, K, beta, achi, z, amp)
    call bjl(l, beta * z, j_l)
    u = symm * amp * z * j_l
    end function olver_reduced

    elemental logical function use_smallchi_map(l, beta, achi, alpha_gate) result(use_smallchi)
    integer, intent(in) :: l
    real(dl), intent(in) :: beta, achi, alpha_gate

    ! Pointwise version of the near-flat small-chi integration gate. The
    ! full integration uses 0.3 with chi_max; using the current achi as the
    ! local endpoint needs a smaller threshold to preserve the 1e-4 pointwise
    ! phi_olver envelope near alpha ~= 3. The gate uses beta/l; the map below
    ! still uses the more accurate sqrt(l(l+1)) curvature scale.
    use_smallchi = l > 0 .and. alpha_gate > SMALLCHI_GATE_ALPHA .and. &
        real(l, dl)**2 * achi**7 / beta < SMALLCHI_GATE_METRIC
    end function use_smallchi_map


    function olver_smallchi_reduced(l, K, beta, achi, symm) result(u)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, achi, symm
    real(dl) :: u
    real(dl) :: z, amp, j_l

    call compute_olver_z_amp_smallchi(l, K, beta, achi, z, amp)
    if (amp <= 0._dl) then
        u = 0._dl
        return
    end if

    call bjl(l, beta * z, j_l)
    u = symm * amp * z * j_l
    end function olver_smallchi_reduced


    function phi_olver_smallchi(l, K, beta, chi) result(phi)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: phi

    real(dl) :: achi, symm

    call validate_inputs(l, K, beta)
    call normalize_chi(l, K, beta, chi, achi, symm)

    if (achi <= CACHE_EPS) then
        phi = 0._dl
        return
    end if

    if (K == 0) then
        call bjl(l, beta * achi, phi)
        return
    end if

    phi = olver_smallchi_reduced(l, K, beta, achi, symm) / curved_radius(K, achi)
    end function phi_olver_smallchi


    pure subroutine compute_olver_z_amp(l, K, beta, achi, z, amp)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, achi
    real(dl), intent(out) :: z
    real(dl), intent(out), optional :: amp

    real(dl) :: ell, alpha, turn_z, turn_chi, sin_k, action, turn_scale

    if (K == 0) then
        z = achi
        if (present(amp)) amp = 1._dl
        return
    end if

    ell = sqrt(real(l, dl) * real(l + 1, dl))
    alpha = beta / ell
    turn_z = ell / beta
    turn_chi = turning_point(ell, beta, K)

    if (achi <= min(1.0e-6_dl, 1.0e-4_dl * max(turn_chi, 1._dl))) then
        z = achi
    else if (abs(achi - turn_chi) <= 1.0e-4_dl * max(turn_chi, turn_z)) then
        if (K == 1) then
            turn_scale = cos(turn_chi)
        else
            turn_scale = cosh(turn_chi)
        end if
        z = turn_z + turn_scale**(1._dl / 3._dl) * (achi - turn_chi)
        if (present(amp)) amp = turn_scale**(-1._dl / 6._dl)
        return
    else
        sin_k = curved_radius(K, achi)
        action = qintegral_exact(sin_k, alpha, K)
        z = invert_flat_action(action, turn_z, achi < turn_chi)
    end if

    if (present(amp)) then
        amp = analytic_amplitude(achi, z, K, alpha, turn_chi)
    end if
    end subroutine compute_olver_z_amp


    elemental subroutine compute_olver_z_amp_smallchi(l, K, beta, chi, z, amp)
    ! Small-chi curvature expansion for the Olver action map.
    !
    ! This solves the differentiated Liouville-Green map perturbatively for
    ! small curvature across the interval,
    !
    !   (dz/dchi)^2 * (alpha^2 - 1/z^2)
    !       = alpha^2 - 1/S_K(chi)^2,
    !
    ! through O((K/alpha^2)^3), writing z = chi * F(alpha*chi,K/alpha^2).
    ! For this local curvature expansion the exact centrifugal coefficient
    ! sqrt(l(l+1)) is more accurate than the Langer parameter at low ell, and
    ! is indistinguishable from it at high ell. The full cached Olver map uses
    ! the same coefficient in compute_olver_z_amp.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl), intent(out) :: z, amp

    real(dl) :: ell2, beta2, chi2
    real(dl) :: rk, h, h2, a, a2, ha
    real(dl) :: F, D

    real(dl), parameter :: c6     = 1._dl / 6._dl
    real(dl), parameter :: c360   = 1._dl / 360._dl
    real(dl), parameter :: c45360 = 1._dl / 45360._dl

    if (K == 0 .or. abs(chi) <= CACHE_EPS) then
        z = chi
        amp = 1._dl
        return
    end if

    rk = real(K, dl)

    ! ell^2 = l(l+1), avoiding sqrt(l(l+1)).
    ell2 = real(l * (l + 1), dl)

    beta2 = beta * beta
    chi2  = chi * chi

    ! h = K / alpha^2 = K * ell^2 / beta^2.
    h = rk * ell2 / beta2

    ! a = h * t^2 = K * chi^2.
    a  = rk * chi2
    a2 = a * a

    h2 = h * h
    ha = h * a

    F = 1._dl - h * ( &
        c6 &
        + c360   * (4._dl   * a  + 13._dl  * h) &
        + c45360 * (48._dl  * a2 + 148._dl * ha + 737._dl * h2) )

    D = 1._dl - h * ( &
        c6 &
        + c360   * (12._dl  * a  + 13._dl  * h) &
        + c45360 * (240._dl * a2 + 444._dl * ha + 737._dl * h2) )

    z = chi * F

    if (D > 0._dl) then
        amp = 1._dl / sqrt(D)
    else
        amp = 0._dl
    end if

    end subroutine compute_olver_z_amp_smallchi

    subroutine validate_inputs(l, K, beta)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta
    integer :: ibeta

    if (l < 0) call MpiStop("Bessel function index ell < 0")
    if (beta < 0._dl) call MpiStop("Wavenumber beta < 0")
    if ((abs(K) /= 1) .and. (K /= 0)) call MpiStop("K must be 1, 0 or -1")

    if (K == 1) then
        ibeta = nint(beta)
        if (ibeta < 3) call MpiStop("Wavenumber beta < 3 for K=1")
        if (ibeta <= l) call MpiStop("Wavenumber beta <= l")
    end if
    end subroutine validate_inputs


    elemental subroutine normalize_chi(l, K, beta, chi, achi, symm)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl), intent(out) :: achi, symm
    integer :: ibeta

    achi = abs(chi)
    symm = 1._dl

    if (K /= 1) return

    ibeta = nint(beta)
    achi = achi - 2._dl * PI * floor(achi / (2._dl * PI))
    if (achi > PI) then
        achi = 2._dl * PI - achi
        if (mod(l, 2) /= 0) symm = -symm
    end if
    if (achi > PI / 2._dl) then
        achi = PI - achi
        if (mod(ibeta - l - 1, 2) /= 0) symm = -symm
    end if
    end subroutine normalize_chi


    elemental real(dl) function turning_point(ell, beta, K)
    real(dl), intent(in) :: ell, beta
    integer, intent(in) :: K

    select case (K)
    case (-1)
        turning_point = asinh(ell / beta)
    case (0)
        turning_point = ell / beta
    case (1)
        turning_point = asin(ell / beta)
    case default
        turning_point = 0._dl
    end select
    end function turning_point


    elemental real(dl) function curved_radius(K, chi)
    integer, intent(in) :: K
    real(dl), intent(in) :: chi

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


    elemental real(dl) function analytic_amplitude(chi, z, K, alpha, turn_chi)
    real(dl), intent(in) :: chi, z, alpha, turn_chi
    integer, intent(in) :: K

    real(dl) :: sin_k, flat_term, curved_term, turn_scale

    if (K == 0) then
        analytic_amplitude = 1._dl
        return
    end if

    if (chi <= 1.0e-8_dl .or. z <= 1.0e-8_dl) then
        analytic_amplitude = 1._dl
        return
    end if

    if (abs(chi - turn_chi) <= 1.0e-10_dl * max(1._dl, turn_chi)) then
        if (K == 1) then
            turn_scale = cos(turn_chi)
        else
            turn_scale = cosh(turn_chi)
        end if
        analytic_amplitude = turn_scale**(-1._dl / 6._dl)
        return
    end if

    sin_k = curved_radius(K, chi)
    flat_term = alpha**2 - 1._dl / z**2
    curved_term = alpha**2 - 1._dl / sin_k**2

    if (abs(flat_term) + abs(curved_term) <= 100._dl * CACHE_EPS * max(1._dl, alpha**2)) then
        if (K == 1) then
            turn_scale = cos(turn_chi)
        else
            turn_scale = cosh(turn_chi)
        end if
        analytic_amplitude = turn_scale**(-1._dl / 6._dl)
        return
    end if

    analytic_amplitude = abs(flat_term / curved_term)**0.25_dl
    end function analytic_amplitude

    elemental real(dl) function invert_flat_action(action, z_turn, below_turn)
    real(dl), intent(in) :: action, z_turn
    logical, intent(in) :: below_turn

    real(dl) :: q, p, s, u
    real(dl) :: A, e
    real(dl) :: x, qplus, invqplus, invqplus2

    q = max(action, 0._dl)
    p = (3._dl * q)**(1._dl / 3._dl)

    ! The near-turning coefficients below are numerical polynomial fits for u(p^2),
    ! not Taylor coefficients. This avoids a root solve in the hot Olver map path.
    if (below_turn) then

        ! Evanescent branch:
        !   q = t - tanh(t),  u = sech(t)
        !
        ! Near the turning point use a polynomial fit in p^2.
        ! Farther out use the large-q asymptotic inversion for sech(t).
        if (p < 1.8_dl) then
            s = p * p

            u = (((((( &
                2.2687533176976695e-05_dl * s &
                - 1.7210470362266376e-04_dl) * s &
                - 3.1231653154203167e-04_dl) * s &
                + 2.1586618823186233e-04_dl) * s &
                + 7.5055874727872618e-02_dl) * s &
                - 5.0000720190143755e-01_dl) * s &
                + 1.0000000000000000_dl)

        else
            A = q + 1._dl
            x = exp(-A)
            e = x * x

            u = 2._dl * x * (1._dl + e * (1._dl + 3._dl * e))
        end if

    else

        ! oscillatory inverse:
        !
        !   q = tan(theta) - theta
        !   u = sec(theta)
        !
        ! Uses a polynomial fit in p^2 near the turning point and a high-order
        ! asymptotic expansion in q + pi/2 farther out.


        if (p < 2._dl) then
            s = p * p

            u = ((((((( &
                3.4828644185221886e-07_dl * s &
                - 8.3097833638549810e-06_dl) * s &
                + 8.8500265224195657e-05_dl) * s &
                - 4.8620948599887724e-04_dl) * s &
                - 3.5056578112088423e-04_dl) * s &
                + 7.4998103457285289e-02_dl) * s &
                + 5.0000018483290443e-01_dl) * s &
                + 1.0000000000000000_dl)

        else
            qplus = q + PI / 2._dl
            invqplus = 1._dl / qplus
            invqplus2 = invqplus * invqplus

            u = qplus - invqplus * ( &
                0.5_dl + invqplus2 * ( &
                7._dl / 24._dl + invqplus2 * ( &
                83._dl / 240._dl + invqplus2 * ( &
                6949._dl / 13440._dl + invqplus2 * ( &
                23399._dl / 26880._dl + invqplus2 * &
                266317._dl / 168960._dl)))))
        end if

    end if

    invert_flat_action = z_turn * u
    end function invert_flat_action

    elemental real(dl) function qintegral_exact(sin_K, alpha, K) result(q)
    real(dl), intent(in) :: sin_K, alpha
    integer, intent(in) :: K

    real(dl), parameter :: zero = 0._dl, one = 1._dl, two = 2._dl, half = 0.5_dl
    real(dl) :: x, x2, a2, ha, r1, r2, u, v, m

    x  = alpha * sin_K
    x2 = x * x
    a2 = alpha * alpha
    ha = half * alpha

    select case (K)

    case (0)
        if (x > one) then
            r1 = sqrt(x2 - one)
            q = r1 - acos(one / x)
        else
            r1 = sqrt(max(zero, one - x2))
            q = log((one + r1) / max(x, CACHE_EPS)) - r1
        end if

    case (-1)
        if (x > one) then
            r1 = sqrt(x2 - one)
            r2 = sqrt(x2 + a2)

            q = ha * log((two*x2 + a2 - one + two*r1*r2) / (one + a2)) &
                - atan2(alpha*r1, r2)
        else
            r1 = sqrt(max(zero, (one - x2) * (x2 + a2)))

            q = ha * atan2(-two*r1, two*x2 + a2 - one) &
                + half * log((two*alpha*r1 + two*a2 + x2*(one - a2)) / &
                max(x2*(one + a2), CACHE_EPS))
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
                v = sqrt(max(zero, a2 - x2))
                q = log(v + alpha*u) - alpha*log(v + u) &
                    - half*log(max(x2, CACHE_EPS)) &
                    + half*(alpha - one)*log(m)
            else
                q = -half*log(max(x2, CACHE_EPS))
            end if
        end if

    end select
    end function qintegral_exact

    end module OlverHypersphericalBessel
