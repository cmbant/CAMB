    module SpherBessels
    !Accurate calculation of ultraspherical Bessel functions for non-flat universe
    !(and imports from FlatBessels for standard approximate spherical Bessel bjl).

    !Hyperspherical bessel routines assume arguments have been pre-validated to be valid.
    use Precision, only: dl
    use MathUtils, only: brentq
    use MpiUtils, only: MpiStop
    use FlatBessels, only: bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor, bjl, Bessels_Free
    implicit none
    private

    public bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor
    public phi_recurs, phi_derivative, phi_first_peak_chi, phi_first_peak_amplitude
    public bjl, Bessels_Free

    type, private :: phi_derivative_root_params
        integer :: l, K
        real(dl) :: nu
    end type phi_derivative_root_params

    contains

    function phi_derivative(l, K, nu, chi) result(dphi)
    ! Derivative d phi_l^nu(K, chi) / d chi from the adjacent-l recurrence.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    real(dl) :: dphi
    real(dl) :: cot_K, root_l

    if (l <= 0) then
        dphi = phi_derivative_l0(K, nu, chi)
        return
    end if

    cot_K = cot_curvature(K, chi)
    select case (K)
    case (-1)
        root_l = sqrt(nu**2 + real(l, dl)**2)
    case (1)
        root_l = sqrt(max(0._dl, nu**2 - real(l, dl)**2))
    case default
        root_l = nu
    end select

    dphi = root_l * phi_recurs(l - 1, K, nu, chi) - real(l + 1, dl) * cot_K * phi_recurs(l, K, nu, chi)
    end function phi_derivative


    function phi_first_peak_chi(l, K, nu, no_peak_found) result(chi_peak)
    ! First maximum at or after the classical turning point, matching the
    ! normalization convention used by the Python mathutils tests.
    !
    ! If the derivative is still positive when the finite search boundary is
    ! reached, no stationary peak has been found. In that case chi_peak is the
    ! search boundary and optional no_peak_found is returned true.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu
    logical, intent(out), optional :: no_peak_found
    real(dl) :: chi_peak

    type(phi_derivative_root_params) :: params
    integer :: iter, iflag
    real(dl) :: turn, upper, hi, delta, cot_t, fturn, fhi
    real(dl) :: froot
    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: ROOT_X = 1.0188_dl
    real(dl), parameter :: MIN_CHI = 1.0e-8_dl
    real(dl), parameter :: MAX_OPEN_CHI = 1.0e4_dl

    if (present(no_peak_found)) no_peak_found = .false.

    turn = hyperspherical_turning_point(l, K, nu)
    if (K == 1) then
        upper = PI / 2._dl
    else
        upper = MAX_OPEN_CHI
    end if
    turn = min(max(turn, MIN_CHI), upper)

    fturn = phi_derivative(l, K, nu, turn)
    if (fturn <= 0._dl .or. turn >= upper) then
        chi_peak = turn
        if (turn >= upper .and. present(no_peak_found)) no_peak_found = .true.
        return
    end if

    cot_t = cot_curvature(K, turn)
    delta = ROOT_X / max(2._dl * nu**2 * abs(cot_t), 1.0e-30_dl)**(1._dl / 3._dl)
    hi = min(turn + 2._dl * delta, upper)
    fhi = phi_derivative(l, K, nu, hi)
    do iter = 1, 100
        if (fhi <= 0._dl) exit
        if (hi >= upper) then
            chi_peak = upper
            if (present(no_peak_found)) no_peak_found = .true.
            return
        end if
        hi = min(turn + 2._dl * (hi - turn), upper)
        fhi = phi_derivative(l, K, nu, hi)
    end do

    if (fhi > 0._dl) then
        chi_peak = hi
        if (present(no_peak_found)) no_peak_found = .true.
    else
        params%l = l
        params%K = K
        params%nu = nu
        call brentq(params, phi_derivative_root, turn, hi, 1.0e-10_dl * max(1._dl, abs(hi)), &
            chi_peak, froot, iflag, fturn, fhi)
        if (iflag /= 0) call MpiStop("phi_first_peak_chi: brentq failed")
    end if
    end function phi_first_peak_chi


    function phi_first_peak_amplitude(l, K, nu, peak_chi, no_peak_found) result(peak)
    ! Absolute amplitude at the first maximum at or after the turning point.
    !
    ! If optional no_peak_found is true, the returned amplitude is evaluated at
    ! peak_chi/search-boundary rather than at a stationary point.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu
    real(dl), intent(out), optional :: peak_chi
    logical, intent(out), optional :: no_peak_found
    real(dl) :: peak
    real(dl) :: chi
    logical :: no_peak

    chi = phi_first_peak_chi(l, K, nu, no_peak)
    if (present(peak_chi)) peak_chi = chi
    if (present(no_peak_found)) no_peak_found = no_peak
    peak = abs(phi_recurs(l, K, nu, chi))
    end function phi_first_peak_amplitude


    pure function hyperspherical_turning_point(l, K, nu) result(turn)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu
    real(dl) :: turn
    real(dl) :: arg

    if (nu <= 0._dl) then
        turn = 0._dl
        return
    end if

    arg = sqrt(real(l, dl) * real(l + 1, dl)) / nu
    select case (K)
    case (-1)
        turn = asinh(arg)
    case (1)
        turn = asin(min(arg, 1._dl))
    case default
        turn = arg
    end select
    end function hyperspherical_turning_point


    pure function cot_curvature(K, chi) result(cot_K)
    integer, intent(in) :: K
    real(dl), intent(in) :: chi
    real(dl) :: cot_K

    select case (K)
    case (-1)
        cot_K = 1._dl / tanh(chi)
    case (1)
        cot_K = 1._dl / tan(chi)
    case default
        cot_K = 1._dl / chi
    end select
    end function cot_curvature


    pure function phi_derivative_l0(K, nu, chi) result(dphi)
    integer, intent(in) :: K
    real(dl), intent(in) :: nu, chi
    real(dl) :: dphi
    real(dl) :: sin_K, cos_K, cot_K, phi0

    if (chi == 0._dl .or. nu == 0._dl) then
        dphi = 0._dl
        return
    end if

    select case (K)
    case (0)
        sin_K = chi
        cos_K = 1._dl
    case (1)
        sin_K = sin(chi)
        cos_K = cos(chi)
    case default
        sin_K = sinh(chi)
        cos_K = cosh(chi)
    end select
    cot_K = cos_K / sin_K
    phi0 = sin(nu * chi) / (nu * sin_K)
    dphi = cos(nu * chi) / sin_K - phi0 * cot_K
    end function phi_derivative_l0


    function phi_derivative_root(obj, chi) result(dphi)
    class(*) :: obj
    real(dl) :: chi
    real(dl) :: dphi

    select type(params => obj)
    type is (phi_derivative_root_params)
        dphi = phi_derivative(params%l, params%K, params%nu, chi)
    class default
        call MpiStop("phi_derivative_root: unexpected parameter type")
    end select
    end function phi_derivative_root

    function phi_recurs(l, K, nu, chi) result(phi)
    ! Recursive evaluation of the regular hyperspherical Bessel function phi_l^nu(K,chi).
    ! Precondition: chi >= 0. For closed K=+1, non-negative chi is folded
    ! into [0,pi/2] using the closed-space parity relations.
    !
    ! The recurrence and exact l=0,1 seeds follow Abbott & Schaefer
    ! (1986, ApJ 308, 546). As in Tram (2017, arXiv:1311.0839) and
    ! Lesgourgues & Tram (2014, arXiv:1312.2697), upward recurrence is used
    ! only in the safe oscillatory region; elsewhere Miller backward
    ! recurrence is started from a stable top boundary condition.
    ! For K=+1 Miller starts use the finite endpoint where
    ! b_j=sqrt(nu^2-j^2) vanishes at j=nu when nu-l>64; closer to the endpoint
    ! they use the closed-space Gegenbauer representation, which is the CLASS
    ! stable-recursion cure for the finite closed-spectrum tail.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, chi
    real(dl) :: phi

    integer :: j, inu
    logical :: use_up, ok
    real(dl) :: nu_use, nu2, ell, chi_use, symm, amp_arg
    real(dl) :: sin_K, cot_K, cos_K, root_K, open_turning_ratio
    real(dl) :: phi0, phi1, phi_top
    real(dl) :: phi_minus, phi_zero, phi_plus
    real(dl) :: b_minus, b_zero
    real(dl) :: cf, bphi_plus
    real(dl) :: phi_cur, phi_lm1, phi0_down, phi1_down
    real(dl) :: scale
    real(dl), parameter :: BIG = 1.e100_dl, TINY = 1.e-280_dl
    real(dl), parameter :: UNDERFLOW_LOG = -744.4400719213812_dl
    real(dl), parameter :: OPEN_TURNING_TOL = 5.e-3_dl
    real(dl), parameter :: OPEN_LOW_BETA_RATIO = 2.0e-3_dl
    real(dl), parameter :: OPEN_LOW_BETA_MIN_TURNING_RATIO = 0.8_dl
    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    integer, parameter :: closed_endpoint_min_degree = 64

    if (K == 1) then
        inu = nint(nu)
        nu_use = real(inu, dl)
    else
        inu = huge(inu)
        nu_use = nu
    end if

    nu2 = nu_use**2
    ell = real(l, dl)

    chi_use = chi
    symm = 1._dl
    if (K == 1) then
        chi_use = modulo(chi_use, 2._dl * PI)
        if (chi_use > PI) then
            chi_use = 2._dl * PI - chi_use
            if (mod(l, 2) /= 0) symm = -symm
        end if
        if (chi_use > PI / 2._dl) then
            chi_use = PI - chi_use
            if (mod(inu - l - 1, 2) /= 0) symm = -symm
        end if
    end if

    if (chi_use == 0._dl) then
        if (l == 0) then
            phi = symm
        else
            phi = 0._dl
        end if
        return
    end if

    if (l > 0) then
        if (K == -1) then
            amp_arg = chi_use * (abs(nu_use) + ell)
        else
            amp_arg = chi_use * abs(nu_use)
        end if

        if (amp_arg == 0._dl) then
            phi = 0._dl
            return
        else if (amp_arg < 1._dl .and. ell * log(amp_arg) < UNDERFLOW_LOG) then
            phi = 0._dl
            return
        end if
    end if

    select case (K)
    case (0)
        sin_K = chi_use
        cot_K = 1._dl / chi_use
    case (1)
        sin_K = sin(chi_use)
        cos_K = cos(chi_use)
        cot_K = cos_K / sin_K
    case (-1)
        sin_K = sinh(chi_use)
        cot_K = 1._dl / tanh(chi_use)
    end select

    call phi01_exact(K, nu_use, chi_use, sin_K, cot_K, phi0, phi1)

    if (l == 0) then
        phi = symm * phi0
        return
    else if (l == 1) then
        phi = symm * phi1
        return
    else if (K == 0 .and. nu_use == 0._dl) then
        phi = 0._dl
        return
    end if

    if (K == 0) then
        root_K = nu
        use_up = (root_K > 0._dl) .and. (abs(cot_K) < root_K / max(1._dl, ell))

        if (use_up) then
            phi_minus = phi0
            phi_zero = phi1

            do j = 2, l
                phi_plus = ((2 * j - 1) * cot_K * phi_zero - nu * phi_minus) / nu

                phi_minus = phi_zero
                phi_zero = phi_plus
            end do

            phi = phi_zero
            return
        end if

        call phi_logderiv(l, K, nu_use, cot_K, cf, ok)
        if (.not. ok) call MpiStop("phi_recurs: failed to get log-derivative")

        phi_cur = 1._dl
        phi_top = 1._dl
        bphi_plus = ell * cot_K - cf
        call rescale_miller_state(phi_cur, bphi_plus, phi_top)
        phi1_down = 0._dl

        do j = l, 1, -1
            phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / nu

            if (j == 1) phi1_down = phi_cur

            bphi_plus = nu * phi_cur
            phi_cur = phi_lm1

            if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
                phi_cur = phi_cur / BIG
                bphi_plus = bphi_plus / BIG
                phi_top = phi_top / BIG
                if (j == 1) phi1_down = phi1_down / BIG
            end if
        end do

        phi0_down = phi_cur

        if (abs(phi0) >= abs(phi1)) then
            if (abs(phi0_down) > TINY) then
                scale = phi0 / phi0_down
            else if (abs(phi1_down) > TINY) then
                scale = phi1 / phi1_down
            else
                call MpiStop("phi_recurs: zero normalization")
            end if
        else
            if (abs(phi1_down) > TINY) then
                scale = phi1 / phi1_down
            else if (abs(phi0_down) > TINY) then
                scale = phi0 / phi0_down
            else
                call MpiStop("phi_recurs: zero normalization")
            end if
        end if

        phi = scale * phi_top
        return
    else if (K == -1) then
        ! For open space the oscillatory/upward-recursion condition
        ! coth(chi) < sqrt(nu^2+l^2)/l is equivalent to
        ! nu*sinh(chi)/l > 1.  Using this form avoids comparing two
        ! quantities that are both very close to one when nu << l.
        ! Just below the open-space turning boundary the continued
        ! fraction can converge slowly; upward recurrence remains
        ! well conditioned in this narrow boundary layer.
        open_turning_ratio = abs(nu_use) * sin_K / ell
        use_up = open_turning_ratio >= 1._dl - OPEN_TURNING_TOL
        if (.not. use_up .and. abs(nu_use) <= OPEN_LOW_BETA_RATIO * ell) then
            ! For very small beta/l, upward recurrence remains accurate
            ! farther below the formal turning boundary, avoiding slow or
            ! stalled continued fractions in open low-beta cases.
            use_up = open_turning_ratio >= OPEN_LOW_BETA_MIN_TURNING_RATIO
        end if
    else !closed
        root_K = sqrt(real(inu - l, dl) * real(inu + l, dl))
        use_up = (root_K > 0._dl) .and. (abs(cot_K) < root_K / max(1._dl, ell))
    end if

    if (use_up) then
        call phi_upward_recur(l, K, inu, nu2, cot_K, phi0, phi1, phi)
        phi = symm * phi
        return
    end if

    if (K == 1 .and. inu - l > closed_endpoint_min_degree) then
        call phi_closed_endpoint_down(l, inu, cot_K, phi0, phi1, phi, ok)
        if (ok) then
            phi = symm * phi
            return
        end if
    end if

    if (K == 1) then
        call phi_closed_gegenbauer_start(l, inu, sin_K, cos_K, phi_cur, bphi_plus, ok)
        if (.not. ok) call MpiStop("phi_recurs: failed to get closed Miller start")
        phi_top = phi_cur
    else
        call phi_logderiv(l, K, nu_use, cot_K, cf, ok)
        if (.not. ok) then
            if (K == -1 .and. abs(nu_use) <= OPEN_LOW_BETA_RATIO * ell) then
                call phi_upward_recur(l, K, inu, nu2, cot_K, phi0, phi1, phi)
                phi = symm * phi
                return
            end if
            call MpiStop("phi_recurs: failed to get log-derivative")
        end if

        phi_cur = 1._dl
        phi_top = 1._dl
        bphi_plus = ell * cot_K - cf
    end if

    call rescale_miller_state(phi_cur, bphi_plus, phi_top)
    phi1_down = 0._dl

    if (K == 1) then
        do j = l, 1, -1
            b_zero = sqrt(real(inu - j, dl) * real(inu + j, dl))

            if (b_zero <= 0._dl) call MpiStop("phi_recurs: zero recurrence coefficient")

            phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

            if (j == 1) phi1_down = phi_cur

            bphi_plus = b_zero * phi_cur
            phi_cur = phi_lm1

            if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
                phi_cur = phi_cur / BIG
                bphi_plus = bphi_plus / BIG
                phi_top = phi_top / BIG
                if (j == 1) phi1_down = phi1_down / BIG
            end if
        end do
    else
        do j = l, 1, -1
            b_zero = sqrt(nu2 + real(j, dl) * real(j, dl))

            phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

            if (j == 1) phi1_down = phi_cur

            bphi_plus = b_zero * phi_cur
            phi_cur = phi_lm1

            if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
                phi_cur = phi_cur / BIG
                bphi_plus = bphi_plus / BIG
                phi_top = phi_top / BIG
                if (j == 1) phi1_down = phi1_down / BIG
            end if
        end do
    end if

    phi0_down = phi_cur

    if (abs(phi0) >= abs(phi1)) then
        if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else
            call MpiStop("phi_recurs: zero normalization")
        end if
    else
        if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else
            call MpiStop("phi_recurs: zero normalization")
        end if
    end if

    phi = symm * scale * phi_top

    contains

    pure subroutine phi_upward_recur(l, K, inu, nu2, cot_K, phi0, phi1, phi)
    integer, intent(in) :: l, K, inu
    real(dl), intent(in) :: nu2, cot_K, phi0, phi1
    real(dl), intent(out) :: phi

    integer :: j
    real(dl) :: phi_minus, phi_zero, phi_plus, b_minus, b_zero

    phi_minus = phi0
    phi_zero = phi1

    if (K == 1) then
        b_minus = sqrt(real(inu - 1, dl) * real(inu + 1, dl))

        do j = 2, l
            b_zero = sqrt(real(inu - j, dl) * real(inu + j, dl))

            phi_plus = ((2 * j - 1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero

            phi_minus = phi_zero
            phi_zero = phi_plus
            b_minus = b_zero
        end do
    else
        b_minus = sqrt(nu2 + 1._dl)

        do j = 2, l
            b_zero = sqrt(nu2 + real(j, dl) * real(j, dl))

            phi_plus = ((2 * j - 1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero

            phi_minus = phi_zero
            phi_zero = phi_plus
            b_minus = b_zero
        end do
    end if

    phi = phi_zero
    end subroutine phi_upward_recur


    pure subroutine phi01_exact(K, nu, chi, sin_K, cot_K, phi0, phi1)
    ! Exact phi_0^nu and phi_1^nu seeds from Abbott & Schaefer
    ! (1986, ApJ 308, 546), with Taylor branches for the small-argument
    ! limits used to normalize Miller recurrence.
    integer, intent(in) :: K
    real(dl), intent(in) :: nu, chi, sin_K, cot_K
    real(dl), intent(out) :: phi0, phi1

    real(dl) :: nu2, kay, arg, arg2, arg4, chi2, chi_over_sin, chi_cot_m1, sinc, sinc_minus_cos, root1

    nu2 = nu**2
    kay = real(K, dl)
    arg = nu * chi
    arg2 = arg**2
    arg4 = arg2**2

    if (abs(arg) < 1.e-4_dl) then
        sinc = 1._dl - arg2 / 6._dl + arg4 / 120._dl
    else
        sinc = sin(arg) / arg
    end if

    if (K == 0) then
        phi0 = sinc

        if (abs(arg) <= 1.e-3_dl) then
            phi1 = arg * (1._dl / 3._dl - arg2 / 30._dl + arg4 / 840._dl)
        else
            phi1 = (sinc - cos(arg)) / arg
        end if

    else
        root1 = sqrt(max(0._dl, nu2 - kay))

        if (abs(chi) < 1.e-4_dl) then
            chi2 = chi**2
            if (abs(arg) < 1.e-4_dl) then
                phi0 = 1._dl - chi2 * (nu2 - kay) / 6._dl
                phi1 = chi * root1 / 3._dl * (1._dl - (3._dl * nu2 - 7._dl * kay) * chi2 / 30._dl)
            else
                chi_over_sin = 1._dl + kay * chi2 / 6._dl + 7._dl * chi2**2 / 360._dl
                chi_cot_m1 = -kay * chi2 / 3._dl - chi2**2 / 45._dl
                phi0 = sinc * chi_over_sin
                if (abs(arg) < 1.e-3_dl) then
                    sinc_minus_cos = arg2 / 3._dl - arg4 / 30._dl + arg4 * arg2 / 840._dl
                else
                    sinc_minus_cos = sinc - cos(arg)
                end if
                phi1 = (sinc_minus_cos + sinc * chi_cot_m1) * chi_over_sin / (root1 * chi)
            end if
        else
            chi_over_sin = chi / sin_K
            chi_cot_m1 = chi * cot_K - 1._dl
            phi0 = sinc * chi_over_sin

            if (abs(arg) < 1.e-3_dl) then
                sinc_minus_cos = arg2 / 3._dl - arg4 / 30._dl + arg4 * arg2 / 840._dl
            else
                sinc_minus_cos = sinc - cos(arg)
            end if

            if (root1 > 0._dl) then
                phi1 = (sinc_minus_cos + sinc * chi_cot_m1) * chi_over_sin / (root1 * chi)
            else
                phi1 = 0._dl
            end if
        end if
    end if

    end subroutine phi01_exact


    pure subroutine rescale_miller_state(phi_cur, bphi_plus, phi_top)
    real(dl), intent(inout) :: phi_cur, bphi_plus, phi_top

    do while (max(abs(phi_cur), abs(bphi_plus)) > BIG)
        phi_cur = phi_cur / BIG
        bphi_plus = bphi_plus / BIG
        phi_top = phi_top / BIG
    end do

    end subroutine rescale_miller_state


    pure subroutine phi_closed_endpoint_down(l, inu, cot_K, phi0, phi1, phi, ok)
    ! Closed-space Miller recurrence from the finite endpoint
    ! b_nu phi_nu=0, where b_j=sqrt(nu^2-j^2) for K=+1.
    ! This is algebraically equivalent to the closed finite-spectrum condition
    ! used by Tram (2017, arXiv:1311.0839) and Lesgourgues & Tram
    ! (2014, arXiv:1312.2697). It is used when the equivalent Gegenbauer
    ! start would require degree n=nu-l-1 >= 64.
    integer, intent(in) :: l, inu
    real(dl), intent(in) :: cot_K, phi0, phi1
    real(dl), intent(out) :: phi
    logical, intent(out) :: ok

    integer :: j
    real(dl) :: b_zero, bphi_plus
    real(dl) :: phi_cur, phi_lm1, phi_target, phi0_down, phi1_down, scale

    ok = .false.
    phi = 0._dl
    phi_cur = 1._dl
    bphi_plus = 0._dl
    phi1_down = 0._dl

    do j = inu - 1, l + 1, -1
        b_zero = sqrt(real(inu - j, dl) * real(inu + j, dl))
        if (b_zero <= 0._dl) return

        phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

        bphi_plus = b_zero * phi_cur
        phi_cur = phi_lm1

        if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
            phi_cur = phi_cur / BIG
            bphi_plus = bphi_plus / BIG
        end if
    end do

    phi_target = phi_cur

    do j = l, 1, -1
        b_zero = sqrt(real(inu - j, dl) * real(inu + j, dl))
        if (b_zero <= 0._dl) return

        phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

        if (j == 1) phi1_down = phi_cur

        bphi_plus = b_zero * phi_cur
        phi_cur = phi_lm1

        if (max(abs(phi_cur), abs(bphi_plus), abs(phi_target)) > BIG) then
            phi_cur = phi_cur / BIG
            bphi_plus = bphi_plus / BIG
            phi_target = phi_target / BIG
            if (j == 1) phi1_down = phi1_down / BIG
        end if
    end do

    phi0_down = phi_cur

    if (abs(phi0) >= abs(phi1)) then
        if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else
            return
        end if
    else
        if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else
            return
        end if
    end if

    phi = scale * phi_target
    ok = .true.

    end subroutine phi_closed_endpoint_down


    pure subroutine phi_logderiv(l, K, nu, cot_K, cf, ok)
    ! Continued-fraction logarithmic derivative used to start Miller recurrence
    ! for K=0,-1; this is the stable-recursion construction described by
    ! Tram (2017, arXiv:1311.0839) and used in the non-flat CLASS
    ! implementation of Lesgourgues & Tram (2014, arXiv:1312.2697).
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, cot_K
    real(dl), intent(out) :: cf
    logical, intent(out) :: ok

    integer :: iter, maxiter
    real(dl) :: nu2, aj, bj, fj, Cj, Dj, Delj, root_cur, root_next
    real(dl), parameter :: SMALL = 1.e-100_dl

    ok = .false.
    cf = 0._dl

    nu2 = nu**2
    maxiter = 1000000

    bj = real(l, dl) * cot_K
    fj = bj
    Cj = bj
    Dj = 0._dl

    if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

    root_cur = sqrt(max(0._dl, nu2 - real(K, dl) * real(l + 1, dl)**2))
    if (root_cur <= 0._dl) return

    do iter = 1, maxiter
        root_next = sqrt(max(0._dl, nu2 - real(K, dl) * real(l + iter + 1, dl)**2))
        if (root_next <= 0._dl) return

        aj = -root_cur / root_next
        if (iter == 1) aj = root_cur * aj

        bj = real(2 * (l + iter) + 1, dl) * cot_K / root_next

        Dj = bj + aj * Dj
        if (abs(Dj) < SMALL) Dj = sign(SMALL, Dj + SMALL)

        Cj = bj + aj / Cj
        if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

        Dj = 1._dl / Dj
        Delj = Cj * Dj
        fj = fj * Delj

        if (abs(Delj - 1._dl) < 10._dl * epsilon(1._dl)) then
            cf = fj
            ok = .true.
            return
        end if

        root_cur = root_next
    end do

    end subroutine phi_logderiv


    pure subroutine phi_closed_gegenbauer_start(l, inu, sin_K, cos_K, phi_l, bphi_plus, ok)
    ! Closed-space replacement for the open continued fraction. The finite
    ! Gegenbauer form gives the Miller start directly, following Tram
    ! (2017, arXiv:1311.0839) and Lesgourgues & Tram
    ! (2014, arXiv:1312.2697). The code divides by G=C_n^{l+1}(cos chi)
    ! only when that rescaling is numerically safe.
    integer, intent(in) :: l, inu
    real(dl), intent(in) :: sin_K, cos_K
    real(dl), intent(out) :: phi_l, bphi_plus
    logical, intent(out) :: ok

    integer :: n, alpha, k
    real(dl) :: x
    real(dl) :: Gm2, Gm1, Gk
    real(dl) :: dGm2, dGm1, dGk
    real(dl) :: G, dG
    real(dl), parameter :: SAFE_RATIO = 1.e-100_dl

    ok = .false.
    phi_l = 0._dl
    bphi_plus = 0._dl

    n = inu - l - 1
    if (n < 0) return

    alpha = l + 1
    x = cos_K

    if (n == 0) then
        G = 1._dl
        dG = 0._dl
    else
        Gm2 = 1._dl
        dGm2 = 0._dl

        Gm1 = 2._dl * real(alpha, dl) * x
        dGm1 = 2._dl * real(alpha, dl)

        if (n == 1) then
            G = Gm1
            dG = dGm1
        else
            do k = 2, n
                Gk = (2._dl * real(k + alpha - 1, dl) * x * Gm1 - &
                    real(k + 2 * alpha - 2, dl) * Gm2) / real(k, dl)

                dGk = (2._dl * real(k + alpha - 1, dl) * (Gm1 + x * dGm1) - &
                    real(k + 2 * alpha - 2, dl) * dGm2) / real(k, dl)

                if (max(abs(Gk), abs(dGk), abs(Gm1), abs(dGm1)) > BIG) then
                    Gm2 = Gm2 / BIG
                    Gm1 = Gm1 / BIG
                    Gk = Gk / BIG
                    dGm2 = dGm2 / BIG
                    dGm1 = dGm1 / BIG
                    dGk = dGk / BIG
                end if

                Gm2 = Gm1
                dGm2 = dGm1
                Gm1 = Gk
                dGm1 = dGk
            end do

            G = Gm1
            dG = dGm1
        end if
    end if

    bphi_plus = sin_K * dG
    if (abs(G) > SAFE_RATIO * max(1._dl, abs(bphi_plus))) then
        phi_l = 1._dl
        bphi_plus = bphi_plus / G
    else
        phi_l = G
    end if
    if (max(abs(phi_l), abs(bphi_plus)) <= TINY) return
    ok = .true.

    end subroutine phi_closed_gegenbauer_start

    end function phi_recurs

    end module SpherBessels
