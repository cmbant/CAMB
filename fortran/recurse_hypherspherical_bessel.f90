    function phi_recurs_stable(l, K, beta, chi) result(phi)
    use Precision
    use MpiUtils
    implicit none
    ! Stable version of phi_recurs(l,K,beta,chi).
!
! Uses ordinary upward recurrence only in the safe oscillatory region.
! In the evanescent/tail region it uses Miller backward recurrence,
! initialized by the logarithmic derivative phi'_l/phi_l.
!
! For K=+1 the logarithmic derivative is computed from the exact
! Gegenbauer representation, avoiding the unstable forced upward
! recurrence in the closed I-tail.
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, chi
    real(dl) :: phi

    integer :: j, ibeta
    logical :: use_up, ok
    real(dl) :: beta2, kay, ell, arg
    real(dl) :: sin_K, cot_K, root_K
    real(dl) :: phi0, phi1, phi_top
    real(dl) :: phi_minus, phi_zero, phi_plus
    real(dl) :: b_minus, b_zero
    real(dl) :: cf, bphi_plus
    real(dl) :: phi_cur, phi_lm1, phi0_down, phi1_down
    real(dl) :: scale
    real(dl), parameter :: BIG = 1.d100, TINY = 1.d-280

    if (l < 0) call MpiStop("Bessel function index ell < 0")
    if (beta < 0._dl) call MpiStop("Wavenumber beta < 0")
    if ((abs(K) /= 1) .and. (K /= 0)) call MpiStop("K must be 1, 0 or -1")

    beta2 = beta**2
    kay = dble(K)
    ell = dble(l)
    arg = beta*chi

    if (K == 1) then
        ibeta = nint(beta)
        if (ibeta < 3) call MpiStop("Wavenumber beta < 3 for K=1")
        if (ibeta <= l) call MpiStop("Wavenumber beta <= l")
    else
        ibeta = huge(ibeta)
    end if

    if (abs(chi) < 1.d-14) then
        if (l == 0) then
            phi = 1._dl
        else
            phi = 0._dl
        end if
        return
    end if

    select case (K)
    case (0)
        sin_K = chi
        cot_K = 1._dl/chi
    case (1)
        sin_K = sin(chi)
        cot_K = 1._dl/tan(chi)
    case (-1)
        sin_K = sinh(chi)
        cot_K = 1._dl/tanh(chi)
    end select

    call phi01_exact_stable(K, beta, chi, sin_K, cot_K, phi0, phi1)

    if (l == 0) then
        phi = phi0
        return
    else if (l == 1) then
        phi = phi1
        return
    end if

    if (K == 0) then
        root_K = beta
    else
        root_K = sqrt(max(0._dl, beta2 - kay*ell*ell))
    end if

    ! Same stability test as the old routine, but do NOT force K=+1 upward.
    use_up = (root_K > 0._dl) .and. (abs(cot_K) < root_K/max(1._dl, ell))

    if (use_up) then

        phi_minus = phi0
        phi_zero = phi1

        if (K == 0) then
            b_minus = beta
        else
            b_minus = sqrt(max(0._dl, beta2 - kay))
        end if

        do j = 2, l
            if (K == 0) then
                b_zero = beta
            else
                b_zero = sqrt(max(0._dl, beta2 - kay*dble(j)*dble(j)))
            end if

            phi_plus = ((2*j - 1)*cot_K*phi_zero - b_minus*phi_minus)/b_zero

            phi_minus = phi_zero
            phi_zero = phi_plus
            b_minus = b_zero
        end do

        phi = phi_zero
        return
    end if

    ! Tail branch: Miller backward recurrence from l to 0.
    !
    ! Start with arbitrary phi_l = 1 and use the logarithmic derivative
    !
    !   phi'_l/phi_l = l*cot_K - b_{l+1} phi_{l+1}/phi_l
    !
    ! so b_{l+1} phi_{l+1} = l*cot_K - cf.
    call phi_logderiv_stable(l, K, beta, sin_K, cot_K, cf, ok)
    if (.not. ok) call MpiStop("phi_recurs_stable: failed to get log-derivative")

    phi_cur = 1._dl
    phi_top = 1._dl
    bphi_plus = ell*cot_K - cf

    phi1_down = 0._dl

    do j = l, 1, -1
        if (K == 0) then
            b_zero = beta
        else
            b_zero = sqrt(max(0._dl, beta2 - kay*dble(j)*dble(j)))
        end if

        if (b_zero <= 0._dl) call MpiStop("phi_recurs_stable: zero recurrence coefficient")

        phi_lm1 = ((2*j + 1)*cot_K*phi_cur - bphi_plus)/b_zero

        if (j == 1) phi1_down = phi_cur

        bphi_plus = b_zero*phi_cur
        phi_cur = phi_lm1

        if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
            phi_cur = phi_cur/BIG
            bphi_plus = bphi_plus/BIG
            phi_top = phi_top/BIG
        end if
    end do

    phi0_down = phi_cur

    ! Normalize the arbitrary Miller solution.  Prefer the larger of phi0
    ! and phi1 to avoid loss near zeros.
    if (abs(phi0) >= abs(phi1)) then
        if (abs(phi0_down) > TINY) then
            scale = phi0/phi0_down
        else if (abs(phi1_down) > TINY) then
            scale = phi1/phi1_down
        else
            call MpiStop("phi_recurs_stable: zero normalization")
        end if
    else
        if (abs(phi1_down) > TINY) then
            scale = phi1/phi1_down
        else if (abs(phi0_down) > TINY) then
            scale = phi0/phi0_down
        else
            call MpiStop("phi_recurs_stable: zero normalization")
        end if
    end if

    phi = scale*phi_top

    contains

    subroutine phi01_exact_stable(K, beta, chi, sin_K, cot_K, phi0, phi1)
    integer, intent(in) :: K
    real(dl), intent(in) :: beta, chi, sin_K, cot_K
    real(dl), intent(out) :: phi0, phi1

    real(dl) :: beta2, kay, arg, arg2, sinc, root1

    beta2 = beta**2
    kay = dble(K)
    arg = beta*chi
    arg2 = arg**2

    if (abs(arg) < 1.d-4) then
        sinc = 1._dl - arg2/6._dl + arg2*arg2/120._dl
    else
        sinc = sin(arg)/arg
    end if

    if (K == 0) then
        phi0 = sinc

        if (abs(arg) < 1.d-4) then
            phi1 = arg/3._dl*(1._dl - arg2/10._dl)
        else
            phi1 = (sinc - cos(arg))/arg
        end if

    else
        root1 = sqrt(max(0._dl, beta2 - kay))

        if (abs(chi) < 1.d-4) then
            if (abs(arg) < 1.d-4) then
                phi0 = 1._dl - chi**2*(beta2 - kay)/6._dl
                phi1 = chi*root1/3._dl
            else
                phi0 = sinc
                phi1 = (sinc - cos(arg))/(root1*chi)
            end if
        else
            phi0 = sin(arg)/(beta*sin_K)
            phi1 = sin(arg)*cot_K/(beta*sin_K) - cos(arg)/sin_K
            phi1 = phi1/root1
        end if
    end if

    end subroutine phi01_exact_stable


    subroutine phi_logderiv_stable(l, K, beta, sin_K, cot_K, cf, ok)
    integer, intent(in) :: l, K
    real(dl), intent(in) :: beta, sin_K, cot_K
    real(dl), intent(out) :: cf
    logical, intent(out) :: ok

    integer :: iter, maxiter
    real(dl) :: beta2, aj, bj, fj, Cj, Dj, Delj, sqrttmp
    real(dl), parameter :: SMALL = 1.d-100

    ok = .false.
    cf = 0._dl

    if (K == 1) then
        call phi_logderiv_closed_gegenbauer(l, nint(beta), sin_K, cot_K, cf, ok)
        return
    end if

    beta2 = beta**2
    maxiter = 1000000

    bj = dble(l)*cot_K
    fj = bj
    Cj = bj
    Dj = 0._dl

    if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

    do iter = 1, maxiter
        sqrttmp = sqrt(max(0._dl, beta2 - dble(K)*dble(l + iter + 1)**2))
        if (sqrttmp <= 0._dl) return

        aj = -sqrt(max(0._dl, beta2 - dble(K)*dble(l + iter)**2))/sqrttmp
        if (iter == 1) aj = sqrt(max(0._dl, beta2 - dble(K)*dble(l + 1)**2))*aj

        bj = dble(2*(l + iter) + 1)*cot_K/sqrttmp

        Dj = bj + aj*Dj
        if (abs(Dj) < SMALL) Dj = sign(SMALL, Dj + SMALL)

        Cj = bj + aj/Cj
        if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

        Dj = 1._dl/Dj
        Delj = Cj*Dj
        fj = fj*Delj

        if (abs(Delj - 1._dl) < 10._dl*epsilon(1._dl)) then
            cf = fj
            ok = .true.
            return
        end if
    end do

    end subroutine phi_logderiv_stable


    subroutine phi_logderiv_closed_gegenbauer(l, ibeta, sin_K, cot_K, cf, ok)
    ! K=+1 exact logarithmic derivative:
    !
    !   Phi_l^beta(chi) proportional to sin(chi)^l
    !                              C_{beta-l-1}^{l+1}(cos chi)
    !
    ! so
    !
    !   phi'_l/phi_l = l*cot(chi) - sin(chi) C'(cos chi)/C(cos chi).
    integer, intent(in) :: l, ibeta
    real(dl), intent(in) :: sin_K, cot_K
    real(dl), intent(out) :: cf
    logical, intent(out) :: ok

    integer :: n, alpha, k
    real(dl) :: x
    real(dl) :: Gm2, Gm1, Gk
    real(dl) :: dGm2, dGm1, dGk
    real(dl) :: G, dG

    ok = .false.
    cf = 0._dl

    n = ibeta - l - 1
    if (n < 0) return

    alpha = l + 1
    x = sin_K*cot_K    ! cos(chi)

    if (n == 0) then
        G = 1._dl
        dG = 0._dl
    else
        Gm2 = 1._dl
        dGm2 = 0._dl

        Gm1 = 2._dl*dble(alpha)*x
        dGm1 = 2._dl*dble(alpha)

        if (n == 1) then
            G = Gm1
            dG = dGm1
        else
            do k = 2, n
                Gk = (2._dl*dble(k + alpha - 1)*x*Gm1 &
                    - dble(k + 2*alpha - 2)*Gm2)/dble(k)

                dGk = (2._dl*dble(k + alpha - 1)*(Gm1 + x*dGm1) &
                    - dble(k + 2*alpha - 2)*dGm2)/dble(k)

                if (max(abs(Gk), abs(dGk), abs(Gm1), abs(dGm1)) > BIG) then
                    Gm2 = Gm2/BIG
                    Gm1 = Gm1/BIG
                    Gk = Gk/BIG
                    dGm2 = dGm2/BIG
                    dGm1 = dGm1/BIG
                    dGk = dGk/BIG
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

    if (abs(G) <= TINY) return

    cf = dble(l)*cot_K - sin_K*dG/G
    ok = .true.

    end subroutine phi_logderiv_closed_gegenbauer

    end function phi_recurs_stable
