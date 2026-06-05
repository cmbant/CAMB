    module SpherBessels
    !Calculation of ultraspherical Bessel functions for non-flat universe
    !(and imports from FlatBessels for standard spherical Bessel bjl).
    use Precision
    use results
    use RangeUtils
    use MpiUtils
    use splines
    use FlatBessels, only: bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor, bjl, Bessels_Free
    implicit none
    private

    public bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor
    public phi_recurs, bjl, Bessels_Free

    contains

    function phi_recurs(l, K, nu, chi) result(phi)
    ! Recursive evaluation of the regular hyperspherical Bessel function phi_l^nu(K,chi).
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
    real(dl) :: nu_use, nu2, kay, ell, chi_use, symm
    real(dl) :: sin_K, cot_K, cos_K, root_K
    real(dl) :: phi0, phi1, phi_top
    real(dl) :: phi_minus, phi_zero, phi_plus
    real(dl) :: b_minus, b_zero
    real(dl) :: cf, bphi_plus
    real(dl) :: phi_cur, phi_lm1, phi0_down, phi1_down
    real(dl) :: scale
    real(dl), parameter :: BIG = 1.d100, TINY = 1.d-280
    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    integer, parameter :: closed_endpoint_min_degree = 64

    if (K == 1) then
        inu = nint(nu)
        nu_use = dble(inu)
    else
        inu = huge(inu)
        nu_use = nu
    end if

    nu2 = nu_use**2
    kay = dble(K)
    ell = dble(l)

    chi_use = chi
    symm = 1._dl
    if (K == 1) then
        chi_use = abs(chi_use)
        chi_use = chi_use - 2._dl * PI * floor(chi_use / (2._dl * PI))
        if (chi_use > PI) then
            chi_use = 2._dl * PI - chi_use
            if (mod(l, 2) /= 0) symm = -symm
        end if
        if (chi_use > PI / 2._dl) then
            chi_use = PI - chi_use
            if (mod(inu - l - 1, 2) /= 0) symm = -symm
        end if
    end if

    if (abs(chi_use) < 1.d-14) then
        if (l == 0) then
            phi = symm
        else
            phi = 0._dl
        end if
        return
    end if

    select case (K)
    case (0)
        sin_K = chi_use
        cot_K = 1._dl / chi_use
    case (1)
        sin_K = sin(chi_use)
        cos_K = cos(chi_use)
        cot_K = 1._dl / tan(chi_use)
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
    end if

    if (K == 0) then
        root_K = nu
    else
        root_K = sqrt(max(0._dl, nu2 - kay * ell * ell))
    end if

    use_up = (root_K > 0._dl) .and. (abs(cot_K) < root_K / max(1._dl, ell))

    if (use_up) then
        phi_minus = phi0
        phi_zero = phi1

        if (K == 0) then
            b_minus = nu
        else
            b_minus = sqrt(max(0._dl, nu2 - kay))
        end if

        do j = 2, l
            if (K == 0) then
                b_zero = nu
            else
                b_zero = sqrt(max(0._dl, nu2 - kay * dble(j) * dble(j)))
            end if

            phi_plus = ((2 * j - 1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero

            phi_minus = phi_zero
            phi_zero = phi_plus
            b_minus = b_zero
        end do

        phi = symm * phi_zero
        return
    end if

    if (K == 1 .and. inu - l > closed_endpoint_min_degree) then
        call phi_closed_endpoint_down(l, inu, nu2, cot_K, phi0, phi1, phi, ok)
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
        call phi_logderiv(l, K, nu_use, sin_K, cot_K, cf, ok)
        if (.not. ok) call MpiStop("phi_recurs: failed to get log-derivative")

        phi_cur = 1._dl
        phi_top = 1._dl
        bphi_plus = ell * cot_K - cf
    end if

    phi1_down = 0._dl

    do j = l, 1, -1
        if (K == 0) then
            b_zero = nu
        else
            b_zero = sqrt(max(0._dl, nu2 - kay * dble(j) * dble(j)))
        end if

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

    pure subroutine phi01_exact(K, nu, chi, sin_K, cot_K, phi0, phi1)
    ! Exact phi_0^nu and phi_1^nu seeds from Abbott & Schaefer
    ! (1986, ApJ 308, 546), with Taylor branches for the small-argument
    ! limits used to normalize Miller recurrence.
    integer, intent(in) :: K
    real(dl), intent(in) :: nu, chi, sin_K, cot_K
    real(dl), intent(out) :: phi0, phi1

    real(dl) :: nu2, kay, arg, arg2, sinc, root1

    nu2 = nu**2
    kay = dble(K)
    arg = nu * chi
    arg2 = arg**2

    if (abs(arg) < 1.d-4) then
        sinc = 1._dl - arg2 / 6._dl + arg2 * arg2 / 120._dl
    else
        sinc = sin(arg) / arg
    end if

    if (K == 0) then
        phi0 = sinc

        if (abs(arg) < 1.d-4) then
            phi1 = arg / 3._dl * (1._dl - arg2 / 10._dl)
        else
            phi1 = (sinc - cos(arg)) / arg
        end if

    else
        root1 = sqrt(max(0._dl, nu2 - kay))

        if (abs(chi) < 1.d-4) then
            if (abs(arg) < 1.d-4) then
                phi0 = 1._dl - chi**2 * (nu2 - kay) / 6._dl
                phi1 = chi * root1 / 3._dl
            else
                phi0 = sinc
                phi1 = (sinc - cos(arg)) / (root1 * chi)
            end if
        else
            phi0 = sin(arg) / (nu * sin_K)
            phi1 = sin(arg) * cot_K / (nu * sin_K) - cos(arg) / sin_K
            phi1 = phi1 / root1
        end if
    end if

    end subroutine phi01_exact


    pure subroutine phi_closed_endpoint_down(l, inu, nu2, cot_K, phi0, phi1, phi, ok)
    ! Closed-space Miller recurrence from the finite endpoint
    ! b_nu phi_nu=0, where b_j=sqrt(nu^2-j^2) for K=+1.
    ! This is algebraically equivalent to the closed finite-spectrum condition
    ! used by Tram (2017, arXiv:1311.0839) and Lesgourgues & Tram
    ! (2014, arXiv:1312.2697). It is used when the equivalent Gegenbauer
    ! start would require degree n=nu-l-1 >= 64.
    integer, intent(in) :: l, inu
    real(dl), intent(in) :: nu2, cot_K, phi0, phi1
    real(dl), intent(out) :: phi
    logical, intent(out) :: ok

    integer :: j
    logical :: have_target
    real(dl) :: b_zero, bphi_plus
    real(dl) :: phi_cur, phi_lm1, phi_target, phi0_down, phi1_down, scale

    ok = .false.
    phi = 0._dl
    phi_cur = 1._dl
    bphi_plus = 0._dl
    phi_target = 0._dl
    phi1_down = 0._dl
    have_target = .false.

    do j = inu - 1, 1, -1
        if (j == l) then
            phi_target = phi_cur
            have_target = .true.
        end if

        b_zero = sqrt(max(0._dl, nu2 - dble(j) * dble(j)))
        if (b_zero <= 0._dl) return

        phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

        if (j == 1) phi1_down = phi_cur

        bphi_plus = b_zero * phi_cur
        phi_cur = phi_lm1

        if (max(abs(phi_cur), abs(bphi_plus), abs(phi_target)) > BIG) then
            phi_cur = phi_cur / BIG
            bphi_plus = bphi_plus / BIG
            if (have_target) phi_target = phi_target / BIG
            if (j == 1) phi1_down = phi1_down / BIG
        end if
    end do

    if (.not. have_target) return

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


    pure subroutine phi_logderiv(l, K, nu, sin_K, cot_K, cf, ok)
    ! Continued-fraction logarithmic derivative used to start Miller recurrence
    ! for K=0,-1; this is the stable-recursion construction described by
    ! Tram (2017, arXiv:1311.0839) and used in the non-flat CLASS
    ! implementation of Lesgourgues & Tram (2014, arXiv:1312.2697).
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu, sin_K, cot_K
    real(dl), intent(out) :: cf
    logical, intent(out) :: ok

    integer :: iter, maxiter
    real(dl) :: nu2, aj, bj, fj, Cj, Dj, Delj, sqrttmp
    real(dl), parameter :: SMALL = 1.d-100

    ok = .false.
    cf = 0._dl

    nu2 = nu**2
    maxiter = 1000000

    bj = dble(l) * cot_K
    fj = bj
    Cj = bj
    Dj = 0._dl

    if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

    do iter = 1, maxiter
        sqrttmp = sqrt(max(0._dl, nu2 - dble(K) * dble(l + iter + 1)**2))
        if (sqrttmp <= 0._dl) return

        aj = -sqrt(max(0._dl, nu2 - dble(K) * dble(l + iter)**2)) / sqrttmp
        if (iter == 1) aj = sqrt(max(0._dl, nu2 - dble(K) * dble(l + 1)**2)) * aj

        bj = dble(2 * (l + iter) + 1) * cot_K / sqrttmp

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
    real(dl), parameter :: SAFE_RATIO = 1.d-100

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

        Gm1 = 2._dl * dble(alpha) * x
        dGm1 = 2._dl * dble(alpha)

        if (n == 1) then
            G = Gm1
            dG = dGm1
        else
            do k = 2, n
                Gk = (2._dl * dble(k + alpha - 1) * x * Gm1 - dble(k + 2 * alpha - 2) * Gm2) / dble(k)

                dGk = (2._dl * dble(k + alpha - 1) * (Gm1 + x * dGm1) - dble(k + 2 * alpha - 2) * dGm2) / dble(k)

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
