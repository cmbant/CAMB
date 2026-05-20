    module OlverHypersphericalBessel
    use Precision
    use MpiUtils
    use SpherBessels, only: bjl, phi_recurs
    use splines, only: spline_def
    implicit none
    private

    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: CACHE_EPS = 1.0e-12_dl

    ! Universal flat-action splines, parameterized by p = (3 q)^(1/3),
    ! sampled uniformly in p so runtime lookup is direct-indexed (no bisection).
    ! Depend only on the sign of (z - z_turn), not on l, beta, or K.
    real(dl), allocatable, save :: univ_u_ev(:), univ_y2_ev(:)
    real(dl), allocatable, save :: univ_theta_os(:), univ_y2_theta_os(:)
    real(dl), save :: univ_dp_ev = 0._dl, univ_dp_os = 0._dl
    real(dl), save :: univ_p_max_ev = 0._dl, univ_p_max_os = 0._dl
    logical, save :: univ_cache_initialized = .false.

    public :: phi_olver_cached, phi_olver_smallchi, olver_cached_coordinate, compute_olver_z_amp_smallchi

    contains

    function olver_cached_coordinate(l, K, beta, chi) result(z)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: z, achi, symm

    call normalize_chi(l, K, beta, chi, achi, symm)
    call compute_olver_z_amp(l, K, beta, achi, z)
    end function olver_cached_coordinate


    function phi_olver_cached(l, K, beta, chi) result(phi)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: phi

    real(dl) :: achi, symm, z, amp, sin_k, j_l, u

    call validate_inputs(l, K, beta)
    call normalize_chi(l, K, beta, chi, achi, symm)

    if (l <= 1) then
        phi = phi_recurs(l, K, beta, achi)
        if (K == 1) phi = symm * phi
        return
    end if

    if (achi <= CACHE_EPS) then
        phi = 0._dl
        return
    end if

    if (K == 0) then
        call bjl(l, beta * achi, phi)
        return
    end if

    call compute_olver_z_amp(l, K, beta, achi, z, amp)
    sin_k = curved_radius(K, achi)

    call bjl(l, beta * z, j_l)
    u = amp * z * j_l
    phi = symm * u / sin_k
    end function phi_olver_cached


    function phi_olver_smallchi(l, K, beta, chi) result(phi)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: phi

    real(dl) :: achi, symm, z, amp, sin_k, j_l

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

    call compute_olver_z_amp_smallchi(l, K, beta, achi, z, amp)
    if (amp <= 0._dl) then
        phi = 0._dl
        return
    end if

    sin_k = curved_radius(K, achi)
    call bjl(l, beta * z, j_l)
    phi = symm * amp * z * j_l / sin_k
    end function phi_olver_smallchi


    subroutine compute_olver_z_amp(l, K, beta, achi, z, amp)
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

    if (l > 0) then
        ell = sqrt(real(l, dl) * (real(l, dl) + 1._dl))
    else
        ell = real(l, dl) + 0.5_dl
    end if
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
    else if (abs(achi - turn_chi) <= 1.0e-10_dl * max(1._dl, turn_chi)) then
        z = turn_z
    else
        sin_k = curved_radius(K, achi)
        action = qintegral_exact(sin_k, alpha, K)
        z = invert_flat_action(action, turn_z, achi < turn_chi)
    end if

    if (present(amp)) then
        amp = analytic_amplitude(achi, z, K, alpha, turn_chi)
    end if
    end subroutine compute_olver_z_amp


    subroutine compute_olver_z_amp_smallchi(l, K, beta, chi, z, amp)
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

    real(dl) :: ell, alpha, t, t2, t4, h, h2, h3, F, D

    if (K == 0) then
        z = chi
        amp = 1._dl
        return
    end if

    if (abs(chi) <= CACHE_EPS) then
        z = chi
        amp = 1._dl
        return
    end if

    if (l > 0) then
        ell = sqrt(real(l, dl) * (real(l, dl) + 1._dl))
    else
        ell = real(l, dl) + 0.5_dl
    end if
    alpha = beta / ell

    t = alpha * chi
    t2 = t * t
    t4 = t2 * t2

    h = real(K, dl) / (alpha * alpha)
    h2 = h * h
    h3 = h2 * h

    F = 1._dl &
        - h / 6._dl &
        - h2 * (4._dl * t2 + 13._dl) / 360._dl &
        - h3 * (48._dl * t4 + 148._dl * t2 + 737._dl) / 45360._dl

    D = 1._dl &
        - h / 6._dl &
        - h2 * (12._dl * t2 + 13._dl) / 360._dl &
        - h3 * (240._dl * t4 + 444._dl * t2 + 737._dl) / 45360._dl

    z = chi * F
    if (D <= 0._dl) then
        amp = 0._dl
    else
        amp = 1._dl / sqrt(D)
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


    subroutine normalize_chi(l, K, beta, chi, achi, symm)
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


    real(dl) function turning_point(ell, beta, K)
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


    real(dl) function curved_radius(K, chi)
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


    real(dl) function analytic_amplitude(chi, z, K, alpha, turn_chi)
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


    real(dl) function invert_flat_action(action, z_turn, below_turn)
    real(dl), intent(in) :: action, z_turn
    logical, intent(in) :: below_turn

    real(dl) :: q, p, u, theta

    if (.not. univ_cache_initialized) call init_universal_olver_cache()

    q = max(action, 0._dl)
    p = (3._dl * q)**(1._dl / 3._dl)

    if (below_turn) then
        if (p < univ_p_max_ev) then
            u = uniform_spline_eval(univ_u_ev, univ_y2_ev, univ_dp_ev, p)
        else
            u = 1._dl / cosh(q + 1._dl)
        end if
    else
        if (p < univ_p_max_os) then
            theta = uniform_spline_eval(univ_theta_os, univ_y2_theta_os, univ_dp_os, p)
            u = 1._dl / cos(theta)
        else
            u = q + PI / 2._dl
        end if
    end if

    invert_flat_action = z_turn * u
    end function invert_flat_action


    subroutine init_universal_olver_cache()
    integer, parameter :: n_ev = 1000, n_os = 2000
    real(dl), parameter :: t_max = 12._dl
    real(dl), parameter :: theta_max = PI / 2._dl - 1.0e-5_dl
    integer :: i
    real(dl) :: p, q
    real(dl) :: p_ev(n_ev), p_os(n_os)

    !$OMP CRITICAL (olver_univ_cache_init)
    if (.not. univ_cache_initialized) then
        univ_p_max_ev = (3._dl * (t_max - tanh(t_max)))**(1._dl / 3._dl)
        univ_p_max_os = (3._dl * (tan(theta_max) - theta_max))**(1._dl / 3._dl)
        univ_dp_ev = univ_p_max_ev / real(n_ev - 1, dl)
        univ_dp_os = univ_p_max_os / real(n_os - 1, dl)

        allocate(univ_u_ev(n_ev), univ_y2_ev(n_ev))
        p_ev(1) = 0._dl
        univ_u_ev(1) = 1._dl
        do i = 2, n_ev
            p = real(i - 1, dl) * univ_dp_ev
            q = p**3 / 3._dl
            p_ev(i) = p
            univ_u_ev(i) = 1._dl / cosh(solve_evanescent_t(q))
        end do
        call spline_def(p_ev, univ_u_ev, n_ev, univ_y2_ev)

        allocate(univ_theta_os(n_os), univ_y2_theta_os(n_os))
        p_os(1) = 0._dl
        univ_theta_os(1) = 0._dl
        do i = 2, n_os
            p = real(i - 1, dl) * univ_dp_os
            q = p**3 / 3._dl
            p_os(i) = p
            univ_theta_os(i) = solve_oscillatory_theta(q)
        end do
        call spline_def(p_os, univ_theta_os, n_os, univ_y2_theta_os)

        univ_cache_initialized = .true.
    end if
    !$OMP END CRITICAL (olver_univ_cache_init)
    end subroutine init_universal_olver_cache


    real(dl) function solve_evanescent_t(q) result(t)
    real(dl), intent(in) :: q
    real(dl) :: low, high, mid
    integer :: iter

    low = q
    high = q + 1._dl
    do iter = 1, 60
        mid = 0.5_dl * (low + high)
        if (mid - tanh(mid) > q) then
            high = mid
        else
            low = mid
        end if
    end do
    t = 0.5_dl * (low + high)
    end function solve_evanescent_t


    real(dl) function solve_oscillatory_theta(q) result(theta)
    real(dl), intent(in) :: q
    real(dl) :: low, high, mid
    integer :: iter

    low = 0._dl
    high = PI / 2._dl - 1.0e-15_dl
    do iter = 1, 70
        mid = 0.5_dl * (low + high)
        if (tan(mid) - mid > q) then
            high = mid
        else
            low = mid
        end if
    end do
    theta = 0.5_dl * (low + high)
    end function solve_oscillatory_theta


    real(dl) function uniform_spline_eval(y, y2, dx, x0) result(val)
    real(dl), intent(in) :: y(:), y2(:), dx, x0
    integer :: idx, n
    real(dl) :: a, b, t

    n = size(y)
    t = x0 / dx
    idx = min(max(int(t), 0), n - 2)
    b = t - real(idx, dl)
    a = 1._dl - b
    val = a * y(idx + 1) + b * y(idx + 2) + &
        ((a**3 - a) * y2(idx + 1) + (b**3 - b) * y2(idx + 2)) * dx**2 / 6._dl
    end function uniform_spline_eval


    real(dl) function qintegral_exact(sin_K, alpha, K)
    real(dl), intent(in) :: sin_K, alpha
    integer, intent(in) :: K

    real(dl) :: arg, root1, root2, dummyarg

    arg = alpha * sin_K

    if (K == 0) then
        if (arg > 1._dl) then
            qintegral_exact = sqrt(arg * arg - 1._dl) - acos(1._dl / arg)
        else
            root1 = sqrt(max(0._dl, 1._dl - arg * arg))
            qintegral_exact = log((1._dl + root1) / max(arg, CACHE_EPS)) - root1
        end if
        return
    end if

    if (K == -1) then
        if (arg > 1._dl) then
            root1 = sqrt(arg * arg - 1._dl)
            root2 = sqrt(arg * arg + alpha * alpha)
            qintegral_exact = alpha / 2._dl * log((2._dl * arg * arg + alpha * alpha - 1._dl + &
                2._dl * root1 * root2) / (1._dl + alpha * alpha)) - atan2(alpha * root1, root2)
        else
            root1 = sqrt(max(0._dl, (1._dl - arg * arg) * (arg * arg + alpha * alpha)))
            qintegral_exact = alpha / 2._dl * atan2(-2._dl * root1, 2._dl * arg * arg + alpha * alpha - 1._dl) + &
                0.5_dl * log((2._dl * alpha * root1 + 2._dl * alpha * alpha + arg * arg * (1._dl - alpha * alpha)) / &
                max(arg * arg * (1._dl + alpha * alpha), CACHE_EPS))
        end if
        return
    end if

    if (arg > 1._dl) then
        root1 = sqrt(arg * arg - 1._dl)
        root2 = sqrt(max(0._dl, alpha * alpha - arg * arg))
        qintegral_exact = alpha / 2._dl * atan2(2._dl * root1 * root2, alpha * alpha + 1._dl - 2._dl * arg * arg) - &
            atan2(root1 * alpha, root2)
    else
        root1 = sqrt(max(0._dl, (1._dl - arg * arg) * (alpha * alpha - arg * arg)))
        dummyarg = alpha * (1._dl - arg * arg) / max(root1, CACHE_EPS)
        dummyarg = min(max(dummyarg, -1._dl + 1.0e-12_dl), 1._dl - 1.0e-12_dl)
        qintegral_exact = 0.5_dl * log((1._dl + dummyarg) / (1._dl - dummyarg)) - alpha / 2._dl * &
            log((alpha * alpha - 2._dl * arg * arg + 1._dl + 2._dl * root1) / max(alpha * alpha - 1._dl, CACHE_EPS))
    end if
    end function qintegral_exact

    end module OlverHypersphericalBessel
