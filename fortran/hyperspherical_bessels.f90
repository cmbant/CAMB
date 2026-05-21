    module SpherBessels
    use Precision
    use results
    use RangeUtils
    use MpiUtils
    use splines
    use FlatBessels, only: bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor, bjl, Bessels_Free
    implicit none
    private

    public bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor
    public phi_recurs_stable, phi_langer, bjl, Bessels_Free

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    ! Calculation of ultraspherical Bessel functions.                      c
    ! Fortran version of the c program hyperjl.c by Arthur Kosowsky.       c
    ! WKB approx described in astro-ph/9805173                             c
    !                                                                      c
    ! Modifications by Anthony Challinor and Antony Lewis                  c
    ! Minor modifications to correct K=1 case outside [0,pi],              c
    ! the small chi approximations for lSamp%l=0 and lSamp%l=1, and                    c
    ! the quadratic approximation to Q(x) around Q(x)=0.                   c
    ! Bug fixed in downwards recursion for the hyperspherical solution     c
    !                                                                      c
    ! The stable recursion routine uses l-recursion to calculate           c
    ! the functions, which is accurate but relatively slow.                c
    !                                                                      c
    ! The routine phi_langer uses Langer's formula for a                   c
    ! uniform first-order asymptotic approximation in the open, closed     c
    ! and flat cases. This approximation is good at high L                 c
    !                                                                      c
    ! The routine qintegral calculates the closed-form answer              c
    ! to the eikonal integral used in the WKB approximation.               c
    !                                                                      c
    ! The routine airy_ai returns the Airy function Ai(x) of the argument  c
    ! passed. It employs a Pade-type approximation away from zero and      c
    ! a Taylor expansion around zero. Highly accurate.                     c
    !                                                                      c
    ! The routines polevl and p1evl are auxiliary polynomial               c
    ! evaluation routines used in the airy function calculation.           c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !Calculates y1,y2 (noramlized to a value near the turning point)
    !by integrating up the differential equation and normalizing to phi_recurs_stable
    !in the region in which phi_recurs_stable is stable
    !This allows closed functions to be computed where chi << turning point
    subroutine phi_small_closed_int(l,beta,chi,y1,y2)
    integer, intent(IN) :: l
    real(dl), intent(IN) :: beta, chi
    real(dl) y1,y2

    integer nsteps,i

    real(dl) ap1,nu2,dydchi1,dydchi2,yt1,yt2,dyt1,dyt2,dym1,dym2
    real(dl) x0, delchi,sh, h6,x
    real(dl) y1_x,y2_x, tmp,xh,hh

    nsteps = 200
    ap1 = l*(l+1)
    x0 = sqrt(ap1)/beta
    nu2 = beta**2

    if ((beta*chi)**2/l < 0.005) then
        !Series solution

        x = chi
        sh = sin(x)
        tmp=(ap1/sh**2 - nu2)
        y1=1e-20
        y2 = ((l+1)/x - (nu2-ap1/3)/(2*l+3)*x) * y1
    else

        x = max(1d-7,chi - 50._dl/l)
        y1=1e-20
        y2 = (l+1)*y1/x

        delchi = (chi-x)/nSteps
        h6=delchi/6
        hh=delchi/2
        sh = sin(x)
        tmp=(ap1/sh**2 - nu2)

        do i=1,nSteps
            ! One step in the ujl integration
            ! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2
            xh=x+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/sin(xh)**2 - nu2)


            dyt2=tmp*yt1       !deriv y2 at mid

            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid

            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+delchi*dym1 !y1 at end
            dym1=dyt1+dym1
            yt2=y2+delchi*dym2 !y2 at end
            dym2=dyt2+dym2

            x=x+delchi     !end point
            sh=sin(x)
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2*dym2)
            if (y1 > 1d10 .or. y2> 1d10) then

                y1=y1/1d10
                y2=y2/1d10

            end if
        end do

    end if

    y1_x = y1; y2_x = y2

    delchi = (x0 - chi)/nSteps
    h6=delchi/6
    hh=delchi/2

    do i=1,nSteps
        ! One step in the ujl integration
        ! fourth-order Runge-Kutta method to integrate equation for ujl

        dydchi1=y2         !deriv y1
        dydchi2=tmp*y1     !deriv y2
        xh=x+hh          !midpoint of step
        yt1=y1+hh*dydchi1  !y1 at midpoint
        yt2=y2+hh*dydchi2  !y2 at midpoint
        dyt1=yt2           !deriv y1 at mid
        tmp=(ap1/sin(xh)**2 - nu2)


        dyt2=tmp*yt1       !deriv y2 at mid

        yt1=y1+hh*dyt1     !y1 at mid
        yt2=y2+hh*dyt2     !y2 at mid

        dym1=yt2           !deriv y1 at mid
        dym2=tmp*yt1       !deriv y2 at mid
        yt1=y1+delchi*dym1 !y1 at end
        dym1=dyt1+dym1
        yt2=y2+delchi*dym2 !y2 at end
        dym2=dyt2+dym2

        x=x+delchi     !end point
        sh=sin(x)
        dyt1=yt2           !deriv y1 at end
        tmp=(ap1/sh**2 - nu2)
        dyt2=tmp*yt1       !deriv y2 at end
        y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
        y2=y2+h6*(dydchi2+dyt2+2*dym2)
        if (y1 > 1d10 .or. y2 > 1d10) then
            y1=y1/1d10
            y2=y2/1d10
            y1_x = y1_x/1d10
            y2_x = y2_x/1d10

        end if
    end do


    tmp = phi_recurs_stable(l,1,beta,x0)*sin(x0) / y1
    y1 = y1_x * tmp
    y2 = y2_x * tmp


    end subroutine phi_small_closed_int

    !***********************************************************************
    !                                                                      *
    ! Calculates Phi(l,beta,chi) using recursion on l.                     *
    ! See Abbot and Schaefer, ApJ 308, 546 (1986) for needed               *
    ! recursion relations and closed-form expressions for l=0,1.           *
    ! (Note: Their variable y is the same as chi here.)                    *
    !                                                                      *
    ! When the flag direction is negative, downwards recursion on l        *
    ! must be used because the upwards direction is unstable to roundoff   *
    ! errors. The downwards recursion begins with arbitrary values and     *
    ! continues downwards to l=1, where the desired l value is normalized  *
    ! using the closed form solution for l=1. (See, e.g., Numerical        *
    ! Recipes of Bessel functions for more detail)                         *
    !                                                                      *
    !***********************************************************************

    function phi_recurs_stable(l, K, beta, chi) result(phi)
    implicit none
    ! Stable hyperspherical recursion for Phi(l,K,beta,chi).
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
    arg = beta * chi

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
        cot_K = 1._dl / chi
    case (1)
        sin_K = sin(chi)
        cot_K = 1._dl / tan(chi)
    case (-1)
        sin_K = sinh(chi)
        cot_K = 1._dl / tanh(chi)
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
        root_K = sqrt(max(0._dl, beta2 - kay * ell * ell))
    end if

    use_up = (root_K > 0._dl) .and. (abs(cot_K) < root_K / max(1._dl, ell))

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
                b_zero = sqrt(max(0._dl, beta2 - kay * dble(j) * dble(j)))
            end if

            phi_plus = ((2 * j - 1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero

            phi_minus = phi_zero
            phi_zero = phi_plus
            b_minus = b_zero
        end do

        phi = phi_zero
        return
    end if

    call phi_logderiv_stable(l, K, beta, sin_K, cot_K, cf, ok)
    if (.not. ok) call MpiStop("phi_recurs_stable: failed to get log-derivative")

    phi_cur = 1._dl
    phi_top = 1._dl
    bphi_plus = ell * cot_K - cf

    phi1_down = 0._dl

    do j = l, 1, -1
        if (K == 0) then
            b_zero = beta
        else
            b_zero = sqrt(max(0._dl, beta2 - kay * dble(j) * dble(j)))
        end if

        if (b_zero <= 0._dl) call MpiStop("phi_recurs_stable: zero recurrence coefficient")

        phi_lm1 = ((2 * j + 1) * cot_K * phi_cur - bphi_plus) / b_zero

        if (j == 1) phi1_down = phi_cur

        bphi_plus = b_zero * phi_cur
        phi_cur = phi_lm1

        if (max(abs(phi_cur), abs(bphi_plus)) > BIG) then
            phi_cur = phi_cur / BIG
            bphi_plus = bphi_plus / BIG
            phi_top = phi_top / BIG
        end if
    end do

    phi0_down = phi_cur

    if (abs(phi0) >= abs(phi1)) then
        if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else
            call MpiStop("phi_recurs_stable: zero normalization")
        end if
    else
        if (abs(phi1_down) > TINY) then
            scale = phi1 / phi1_down
        else if (abs(phi0_down) > TINY) then
            scale = phi0 / phi0_down
        else
            call MpiStop("phi_recurs_stable: zero normalization")
        end if
    end if

    phi = scale * phi_top

    contains

    subroutine phi01_exact_stable(K, beta, chi, sin_K, cot_K, phi0, phi1)
    integer, intent(in) :: K
    real(dl), intent(in) :: beta, chi, sin_K, cot_K
    real(dl), intent(out) :: phi0, phi1

    real(dl) :: beta2, kay, arg, arg2, sinc, root1

    beta2 = beta**2
    kay = dble(K)
    arg = beta * chi
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
        root1 = sqrt(max(0._dl, beta2 - kay))

        if (abs(chi) < 1.d-4) then
            if (abs(arg) < 1.d-4) then
                phi0 = 1._dl - chi**2 * (beta2 - kay) / 6._dl
                phi1 = chi * root1 / 3._dl
            else
                phi0 = sinc
                phi1 = (sinc - cos(arg)) / (root1 * chi)
            end if
        else
            phi0 = sin(arg) / (beta * sin_K)
            phi1 = sin(arg) * cot_K / (beta * sin_K) - cos(arg) / sin_K
            phi1 = phi1 / root1
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

    bj = dble(l) * cot_K
    fj = bj
    Cj = bj
    Dj = 0._dl

    if (abs(Cj) < SMALL) Cj = sign(SMALL, Cj + SMALL)

    do iter = 1, maxiter
        sqrttmp = sqrt(max(0._dl, beta2 - dble(K) * dble(l + iter + 1)**2))
        if (sqrttmp <= 0._dl) return

        aj = -sqrt(max(0._dl, beta2 - dble(K) * dble(l + iter)**2)) / sqrttmp
        if (iter == 1) aj = sqrt(max(0._dl, beta2 - dble(K) * dble(l + 1)**2)) * aj

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

    end subroutine phi_logderiv_stable


    subroutine phi_logderiv_closed_gegenbauer(l, ibeta, sin_K, cot_K, cf, ok)
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
    x = sin_K * cot_K

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

    if (abs(G) <= TINY) return

    cf = dble(l) * cot_K - sin_K * dG / G
    ok = .true.

    end subroutine phi_logderiv_closed_gegenbauer

    end function phi_recurs_stable


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    ! Calculates Phi(l,beta,chi) using the Langer uniform approximation    c
    ! to the first-order WKB approximation.                                c
    ! See C.M. Bender and S.A. Orszag,  Mathematical Methods for           c
    ! Scientists and Engineers (McGraw-Hill, 1978; LC QA371.B43),          c
    ! chapter 10.                                                          c
    !                                                                      c
    ! Differential equation for needed function can be cast into the       c
    ! Schrodinger form      \epsilon^2 y'' = Q(x) y                        c
    ! where \epsilon^2 = 1/l(l+1) and Q(x) depends on the parameter        c
    ! alpha \equiv beta * epsilon.                                         c
    !                                                                      c
    ! In the K= +1 case, the function is                                   c
    ! determined by its value on the interval [0, pi/2] and the symmetry   c
    ! conditions Phi(chi + pi) = (-1)^{beta - l - 1} Phi(chi),             c
    !            Phi(pi - chi) = (-1)^{beta - l - 1} Phi(chi).             c
    ! This interval contains one turning point, so the Langer formula      c
    ! can be used.                                                         c
    ! Note that the second condition at chi = pi/2 gives an eigenvalue     c
    ! condition on beta, which  must corrected. For the lowest             c
    ! eigenvalue(s), the region between the turning points is not large    c
    ! enough for the asymptotic solution to be valid, so the functions     c
    ! have a small discontinuity or discontinuous derivative at pi/2;      c
    ! this behavior is corrected by employing a 4-term asymptotic          c
    ! series around the regular point chi=pi/2.                            c
    ! The exact eigenvalue condition requires that beta must be an         c
    ! integer >= 3 with beta > l. Note this implies alpha > 1.             c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



    function phi_langer(l,K,beta,chi)
    integer l,K,ibeta,kay
    real(dl) phi_langer
    real(dl) ell,symm, anu, alpha2
    real(dl) beta,chi,eikonal, wkb, arg, arg2, tmp
    real(dl) epsilon, alpha, chi0, x, a, b,achi

    real(dl) cot_K, sin_K
    real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
        PIOVER2=1.570796327d0

    ell=dble(l)
    achi=chi

    symm=1._dl
    !
    ! Test input values
    !
    if(l<0) call MpiStop("Bessel function index ell < 0")
    if(beta<0._dl) call MpiStop("Wavenumber beta < 0")
    if ((abs(K)/=1).and.(K/=0)) call MpiStop("K must be 1, 0 or -1")


    if(K == 1) then
        ibeta=nint(beta)
        if(ibeta<3) call MpiStop("Wavenumber beta < 3 for K=1")
        if(ibeta<=l) call MpiStop("Wavenumber beta <= l")
    endif

    kay=K


    ! For closed case, find equivalent chi in [0,pi/2]
    !
    if(K==1) then
        achi=achi-2._dl*Pi*int(achi/2._dl/PI)
        if(achi>PI) then
            achi=2._dl*PI-achi
            if(2*(l/2).eq.l) then
                symm=symm
            else
                symm=-symm
            endif
        endif
        if(achi>PI/2._dl) then
            achi=PI-achi
            if(2*((ibeta-l-1)/2).eq.(ibeta-l-1)) then
                symm=symm
            else
                symm=-symm
            endif
        endif
    endif

    ! Definitions
    if(K == 0) then
        sin_K = achi
    else
        if(K == -1) then
            sin_K = sinh(achi)
        else
            sin_K = sin(achi)
        end if
    endif

    ! Closed form solution for l=0
    !
    if(l == 0) then
        arg=beta*achi

        if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
                wkb=1._dl-achi*achi*(beta*beta-kay)/6._dl
            else
                wkb=sin(arg)/arg
            endif
        else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=1._dl-arg*arg/6._dl
        else
            wkb=sin(arg)/(beta*sin_K)
        endif
        phi_langer=symm*wkb
        return
    endif
    !
    ! Closed form solution for l=1
    !
    if(l==1) then
        arg=beta*achi

        if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
                wkb=achi*sqrt(beta*beta-kay)/3._dl
            else
                wkb=(sin(arg)/arg-cos(arg))/(sqrt(beta*beta-kay)*achi)
            endif
        else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=arg/3._dl
        else
            if(K/=0) then
                if(K==1) then
                    cot_K=1._dl/tan(achi)
                else
                    cot_K=1._dl/tanh(achi)
                endif
                wkb=sin(arg)*cot_K/(beta*sin_K)-cos(arg)/sin_K
                wkb=wkb/sqrt(beta*beta-kay)
            else
                wkb=(sin(arg)/arg-cos(arg))/arg
            endif
        end if
        phi_langer=symm*wkb
        return
    endif
    !
    ! Closed form solution for K=1 and beta = l+1 (lowest eigenfunction)
    !
    if((K==1).and.(ibeta == (l+1))) then
        wkb=(sin_K**ell)* &
            sqrt(sqrt(2._dl*PI/(2._dl*ell+1._dl))*ell/((ell+1._dl)*(2._dl*ell+1._dl)))
        wkb=wkb*(1+0.1875d0/ell-0.013671875/(ell*ell))
        phi_langer=symm*wkb
        return
    endif

    ! Very close to 0, return 0 (exponentially damped)
    !
    if(abs(achi)<1.d-8) then
        phi_langer=0._dl
        return
    endif


    ! For closed case, find corrected eigenvalue beta
    !
    if(K==1) then
        anu=dble(ibeta)-1._dl/(8._dl*ell)+1._dl/(16._dl*ell*ell)
    else
        anu=beta
    endif
    !
    ! Evaluate epsilon using asymptotic form for large l
    !
    if(l<20) then
        epsilon=1._dl/sqrt(ell*(ell+1._dl))
    else
        epsilon=1._dl/ell-0.5d0/(ell*ell)+0.375d0/(ell*ell*ell)
    endif

    alpha=epsilon*anu
    !
    ! Calculate the turning point where Q(x)=0.
    ! Function in question has only a single simple turning point.
    !
    if(K==-1) chi0=log((1._dl+sqrt(1._dl+alpha*alpha))/alpha)
    if(K==0) chi0=1._dl/alpha
    if(K==1) chi0=asin(1._dl/alpha)


    ! Very close to chi0, use usual wkb form to avoid dividing by zero
    !
    if(abs(achi-chi0)<1.d-5) then

        ! Calculate coefficients of linear and quadratic terms in Q(x) expansion
        ! in the neighborhood of the turning point
        ! Q(chi)=a*(chi0-chi)+b*(chi0-chi)**2
        alpha2=alpha*alpha

        if(K==-1) then
            a=2._dl*alpha2*sqrt(alpha2+1._dl)
            b=3._dl*alpha2**2+2._dl*alpha2
        endif
        if(K==0) then
            a=2._dl*alpha2*alpha
            b=3._dl*alpha2**2
        endif
        if(K==1) then
            a=2._dl*alpha2*sqrt(alpha2-1._dl)
            b=3._dl*alpha2**2-2._dl*alpha2
        endif

        ! Dependent variable x for which Q(x)=0 at x=0
        ! x>0 is the evanescent region
        !
        x=chi0-achi
        !
        ! Argument of Airy function
        !
        arg=(x+b*x*x/(5._dl*a))/(epsilon*epsilon/a)**(0.333333333d0)
        !
        ! Evaluate Airy function
        !
        wkb=airy_ai(arg)
        !
        ! Rest of functional dependence
        !
        wkb=wkb*(1._dl-b*x/(5._dl*a))/sin_K
        !  Normalization factor:

        wkb=wkb*symm*ROOTPI*((a*epsilon)**(-0.1666667d0))*sqrt(epsilon/anu)
        phi_langer=wkb

        return
    endif


    ! Langer approximation.
    !
    ! Transport factor:
    !
    tmp=sqrt(abs(1._dl/(sin_K*sin_K)-alpha*alpha))

    ! Eikonal factor
    !
    eikonal=qintegral(sin_K,alpha,K)



    arg=(1.5d0*eikonal/epsilon)**(1._dl/3._dl)

    arg2=arg*arg
    if(achi>chi0) arg2=-arg2

    ! Multiply factors together

    wkb=airy_ai(arg2)*symm*ROOTPI*sqrt(arg*epsilon/anu/tmp)/Sin_K
    phi_langer=wkb

    end function phi_langer

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                       c
    ! Evaluates the exact solution to  the integral giving the WKB          c
    ! eikonal solution,   \int^x sqrt(abs(Q(x))) dx                         c
    !                                                                       c
    ! In the open case, this integral costs 1 or 2 square roots, an atan    c
    ! and a log; its evaluation will be roughly as expensive as the rest    c
    ! of the Phi routine. An analytic fit cannot be substantially faster    c
    ! because the dependence on alpha of the y-intercept of the linear      c
    ! region of the integrand contains a log and an atan, so at best a fit  c
    ! can only save the computation of the square roots.                    c
    !                                                                       c
    ! The integrals are very bland functions of chi and alpha and could     c
    ! be precomputed and cached to save computation time; interpolation     c
    ! on a relatively small number of points should be very accurate        c
    !                                                                       c
    ! Note that for the closed case, the variable arg must be between 0     c
    ! and alpha; for arg > alpha, symmetry properties reduce the needed     c
    ! integral to this case.                                                c
    !                                                                       c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function qintegral(sin_K,alpha, K)
    implicit none
    real(dl) qintegral, sin_K
    integer K
    real(dl) alpha,exact,arg, root1, root2, dummyarg

    real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
        PIOVER2=1.570796327d0

    arg=alpha*sin_K

    if(K==0) then
        if(arg>1._dl) then
            exact=sqrt(arg*arg-1._dl)-acos(1._dl/arg)
            qintegral=exact
            return
        else
            root1=sqrt(1._dl-arg*arg)
            exact=log((1._dl+root1)/arg)-root1
            qintegral=exact
            return
        endif
    else if(K==-1) then
        if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(arg*arg+alpha*alpha)
            exact=alpha/2._dl*log((2._dl*arg*arg+alpha*alpha-1._dl+ &
                2._dl*root1*root2)/(1._dl+alpha*alpha))+atan(root2/ &
                (alpha*root1))-PIOVER2
            qintegral=exact
            return
        else
            root1=sqrt((1._dl-arg*arg)*(arg*arg+alpha*alpha))
            exact=alpha/2._dl*atan(-2._dl*root1/(2*arg*arg+alpha*alpha- &
                1._dl))+0.5d0*log((2._dl*alpha*root1+2._dl*alpha*alpha+ &
                arg*arg*(1._dl-alpha*alpha))/(arg*arg*(1._dl+ &
                alpha*alpha)))
            if(2._dl*arg*arg+alpha*alpha-1._dl<0._dl) then
                exact=exact-alpha*PIOVER2
            endif
            qintegral=exact
            return
        endif
    else
        if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(alpha*alpha-arg*arg)
            exact=alpha/2._dl*atan(-2._dl*root1*root2/ &
                (2._dl*arg*arg-alpha*alpha-1._dl))- &
                atan(-root2/(root1*alpha))-PIOVER2
            if(2._dl*arg*arg-alpha*alpha-1._dl>0._dl) then
                exact=exact+alpha*PIOVER2
            endif
        else
            root1=sqrt((1._dl-arg*arg)*(alpha*alpha-arg*arg))
            dummyarg=alpha*(1._dl-arg*arg)/root1
            exact=0.5d0*log((1._dl+dummyarg)/(1._dl-dummyarg))- &
                alpha/2._dl*log((alpha*alpha-2._dl*arg*arg+1._dl+2._dl*root1)/ &
                (alpha*alpha-1._dl))
        endif
        qintegral=exact
        return
    endif

    end function qintegral

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !       Airy function                                                  c
    !                                                                      c
    ! Modified from original routine by Stephen Moshier, available         c
    ! as part of the Cephes library at www.netlib.com                      c
    ! Modifications: eliminates calculation of Bi(x), Ai'(x), Bi'(x)       c
    ! and translation to Fortran                                           c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! DESCRIPTION:
    !
    ! Solution of the differential equation
    !
    !       y"(x) = xy.
    !
    ! The function returns the two independent solutions Ai, Bi
    ! and their first derivatives Ai'(x), Bi'(x).
    !
    ! Evaluation is by power series summation for small x,
    ! by rational minimax approximations for large x.
    !
    !
    !
    ! ACCURACY:
    ! Error criterion is absolute when function <= 1, relative
    ! when function > 1, except * denotes relative error criterion.
    ! For large negative x, the absolute error increases as x^1.5.
    ! For large positive x, the relative error increases as x^1.5.
    !
    ! Arithmetic  domain   function  # trials      peak         rms
    ! IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
    ! IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
    ! IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
    ! IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
    ! IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
    ! IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
    ! DEC         -10, 0     Ai         5000       1.7e-16     2.8e-17
    ! DEC           0, 10    Ai         5000       2.1e-15*    1.7e-16*
    ! DEC         -10, 0     Ai'        5000       4.7e-16     7.8e-17
    ! DEC           0, 10    Ai'       12000       1.8e-15*    1.5e-16*
    ! DEC         -10, 10    Bi        10000       5.5e-16     6.8e-17
    ! DEC         -10, 10    Bi'        7000       5.3e-16     8.7e-17
    !
    !
    ! Cephes Math Library Release 2.1:  January, 1989
    ! Copyright 1984, 1987, 1989 by Stephen lSamp%l. Moshier
    ! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
    !

    function airy_ai(x)
    implicit none
    real(dl) airy_ai
    real(dl) x,z, zz, t, f, g, uf, ug, zeta, theta
    real(dl) ak
    real(dl) AN(8),AD(8),AFN(9),AFD(9),AGN(11),AGD(10)
    real(dl), parameter :: AMAXAIRY=25.77d0,ACC=1.d-8,PI=3.1415926536d0
    real(dl), parameter :: c1=0.35502805388781723926d0, c2=0.258819403792806798405d0
    real(dl), parameter :: sqrt3=1.732050807568877293527d0,sqpii=5.64189583547756286948d-1


    AN(1)=3.46538101525629032477d-1
    AN(2)=1.20075952739645805542d1
    AN(3)=7.62796053615234516538d1
    AN(4)=1.68089224934630576269d2
    AN(5)=1.59756391350164413639d2
    AN(6)=7.05360906840444183113d1
    AN(7)=1.40264691163389668864d1
    AN(8)=9.99999999999999995305d-1

    AD(1)=5.67594532638770212846d-1
    AD(2)=1.47562562584847203173d1
    AD(3)=8.45138970141474626562d1
    AD(4)=1.77318088145400459522d2
    AD(5)=1.64234692871529701831d2
    AD(6)=7.14778400825575695274d1
    AD(7)=1.40959135607834029598d1
    AD(8)=1.00000000000000000470d0

    AFN(1)=-1.31696323418331795333d-1
    AFN(2)=-6.26456544431912369773d-1
    AFN(3)=-6.93158036036933542233d-1
    AFN(4)=-2.79779981545119124951d-1
    AFN(5)=-4.91900132609500318020d-2
    AFN(6)=-4.06265923594885404393d-3
    AFN(7)=-1.59276496239262096340d-4
    AFN(8)=-2.77649108155232920844d-6
    AFN(9)=-1.67787698489114633780d-8

    AFD(1)=1.33560420706553243746d1
    AFD(2)=3.26825032795224613948d1
    AFD(3)=2.67367040941499554804d1
    AFD(4)=9.18707402907259625840d0
    AFD(5)=1.47529146771666414581d0
    AFD(6)=1.15687173795188044134d-1
    AFD(7)=4.40291641615211203805d-3
    AFD(8)=7.54720348287414296618d-5
    AFD(9)=4.51850092970580378464d-7

    AGN(1)=1.97339932091685679179d-2
    AGN(2)=3.91103029615688277255d-1
    AGN(3)=1.06579897599595591108d0
    AGN(4)=9.39169229816650230044d-1
    AGN(5)=3.51465656105547619242d-1
    AGN(6)=6.33888919628925490927d-2
    AGN(7)=5.85804113048388458567d-3
    AGN(8)=2.82851600836737019778d-4
    AGN(9)=6.98793669997260967291d-6
    AGN(10)=8.11789239554389293311d-8
    AGN(11)=3.41551784765923618484d-10

    AGD(1)=9.30892908077441974853d0
    AGD(2)=1.98352928718312140417d1
    AGD(3)=1.55646628932864612953d1
    AGD(4)=5.47686069422975497931d0
    AGD(5)=9.54293611618961883998d-1
    AGD(6)=8.64580826352392193095d-2
    AGD(7)=4.12656523824222607191d-3
    AGD(8)=1.01259085116509135510d-4
    AGD(9)=1.17166733214413521882d-6
    AGD(10)=4.91834570062930015649d-9
    !
    ! Exponentially tiny for large enough argument
    !
    if(x>AMAXAIRY) then
        airy_ai=0._dl
        return
    endif
    !
    ! Pade fit for large negative arguments
    !
    if(x<-2.09d0) then
        t=sqrt(-x)
        zeta=-2._dl*x*t/3._dl
        t=sqrt(t)
        ak=sqpii/t
        z=1._dl/zeta
        zz=z*z
        uf=1._dl+zz*polevl(zz,AFN,8)/p1evl(zz,AFD,9)
        ug=z*polevl(zz,AGN,10)/p1evl(zz,AGD,10)
        theta=zeta+0.25d0*PI
        f=sin(theta)
        g=cos(theta)
        airy_ai=ak*(f*uf-g*ug)
        return
    endif
    !
    ! Pade fit for large positive arguments
    !
    if(x>=2.09) then
        t=sqrt(x)
        zeta=2._dl*x*t/3._dl
        g=exp(zeta)
        t=sqrt(t)
        ak=2._dl*t*g
        z=1._dl/zeta
        f=polevl(z,AN,7)/polevl(z,AD,7)
        airy_ai=sqpii*f/ak
        return
    endif
    !
    ! Taylor series for region around x=0
    !

    f=1._dl
    g=x
    t=1._dl
    uf=1._dl
    ug=x
    ak=1._dl
    z=x*x*x
    do while (t>ACC)
        uf=uf*z
        ak=ak+1._dl
        uf=uf/ak
        ug=ug*z
        ak=ak+1._dl
        ug=ug/ak
        uf=uf/ak
        f=f+uf
        ak=ak+1._dl
        ug=ug/ak
        g=g+ug
        t=abs(uf/f)
    end do


    airy_ai=c1*f-c2*g

    end function airy_ai

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  Evaluate polynomial                                          c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! DESCRIPTION:
    !
    ! Evaluates polynomial of degree N:
    !
    !                     2          N
    ! y  =  C  + C x + C x  +...+ C x
    !        0    1     2          N
    !
    ! Coefficients are stored in reverse order:
    !
    ! coef(1) = C  , ..., coef(N+1) = C  .
    !            N                     0
    !
    ! The function p1evl() assumes that C = 1.0 and is
    !                                    N
    ! omitted from the array.  Its calling arguments are
    ! otherwise the same as polevl().
    !
    !
    ! SPEED:
    !
    ! In the interest of speed, there are no checks for out
    ! of bounds arithmetic.  This routine is used by most of
    ! the functions in the library.  Depending on available
    ! equipment features, the user may wish to rewrite the
    ! program in microcode or assembly language.
    !
    ! Cephes Math Library Release 2.1:  December, 1988
    ! Copyright 1984, 1987, 1988 by Stephen lSamp%l. Moshier
    ! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
    !
    function polevl(x,coef,N)
    implicit none
    real(dl) polevl
    real(dl) x,ans
    real(dl) coef
    integer N,i

    dimension coef(N+1)

    ans=coef(1)
    do i=2,N+1
        ans=ans*x+coef(i)
    end do
    polevl=ans

    end function polevl

    !
    !
    ! Evaluate polynomial when coefficient of x  is 1.0.
    ! Otherwise same as polevl.
    !
    function p1evl(x,coef,N)
    implicit none
    real(dl) p1evl
    real(dl) x,coef,ans
    integer N,i
    dimension coef(N)

    ans=x+coef(1)
    do i=2,N
        ans=ans*x+coef(i)
    end do
    p1evl=ans

    end function p1evl



    end module SpherBessels !USpherBessels
