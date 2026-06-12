    module MathUtils
    use precision
    implicit none

    interface
    FUNCTION obj_function(obj, x)
    use precision
    class(*) :: obj
    real(dl) :: x, obj_function
    END FUNCTION  obj_function
    end interface

    contains

    function Integrate_Romberg(obj, fin, a, b, tol, maxit, minsteps, abs_tol)
    !  Rombint returns the integral from a to b of f(obj,x) using Romberg integration.
    !  The method converges provided that f is continuous in (a,b).
    !  f must be real(dl). The first argument is a class instance.
    !  tol indicates the desired relative accuracy in the integral.

    ! Modified by AL to specify max iterations and minimum number of steps
    ! (min steps useful to stop wrong results on periodic or sharp functions)
    use iso_c_binding
    use MiscUtils
    use config, only : global_error_flag, print_fortran_warnings
    class(*) :: obj
    real(dl), external :: fin !a class function
    procedure(obj_function), pointer :: f
    real(dl), intent(in) :: a,b,tol
    integer, intent(in), optional :: maxit,minsteps
    logical, intent(in), optional :: abs_tol
    integer max_it, min_steps
    real(dl) :: Integrate_Romberg
    integer, parameter :: MAXJ=5
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g(MAXJ+1), g0, g1, fourj
    logical abstol

    !convert the class function (un-type-checked) into correct type to call correctly for class argument
    call C_F_PROCPOINTER(c_funloc(fin), f)
    Integrate_Romberg = -1
    max_it = PresentDefault(25, maxit)
    min_steps = PresentDefault(0, minsteps)
    abstol = DefaultFalse(abs_tol)
    h=0.5d0*(b-a)
    gmax=h*(f(obj,a)+f(obj,b))
    if (global_error_flag /=0) return
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
    do
        i=i+1
        if (i > max_it.or.(i > 5.and.abs(error) < tol) .and. nint > min_steps) exit
        !  Calculate next trapezoidal rule approximation to integral.
        g0=0._dl
        do k=1,nint
            g0=g0+f(obj, a+(k+k-1)*h)
            if (global_error_flag /=0) return
        end do
        g0=0.5d0*g(1)+h*g0
        h=0.5d0*h
        nint=nint+nint
        jmax=min(i,MAXJ)
        fourj=1._dl
        do j=1,jmax
            !  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
        end do
        if (abstol) then
            error=abs(gmax-g0)
        else
            if (abs(g0).gt.tol) then
                error=1._dl-gmax/g0
            else
                error=gmax
            end if
        end if
        gmax=g0
        g(jmax+1)=g0
    end do

    Integrate_Romberg=g0
    if (i > max_it .and. abs(error) > tol .and. print_fortran_warnings)  then
        write(*,*) 'Warning: Integrate_Romberg failed to converge; '
        write (*,*)'integral, error, tol:', Integrate_Romberg,error, tol
    end if

    end function Integrate_Romberg


    subroutine brentq(obj,func,ax,bx,tol,xzero,fzero,iflag,fax,fbx)
    use iso_c_binding

    !>
    !  Find a zero of the function \( f(x) \) in the given interval
    !  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
    !  where \( \epsilon \) is the relative machine precision defined as
    !  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
    !
    !  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
    !
    !#References
    !  * R. P. Brent, "[An algorithm with guaranteed convergence for
    !    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
    !    The Computer Journal, Vol 14, No. 4., 1971.
    !  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
    !    Prentice-Hall, Inc., 1973.
    !
    !# See also
    !  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

    use iso_fortran_env, only: error_unit
    implicit none
    class(*) :: obj
    real(dl), external :: func !a class function f(obj,x)
    procedure(obj_function), pointer :: f

    real(dl),intent(in)              :: ax      !! left endpoint of initial interval
    real(dl),intent(in)              :: bx      !! right endpoint of initial interval
    real(dl),intent(in)              :: tol     !! desired length of the interval of uncertainty of the final result (>=0)
    real(dl),intent(out)             :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(dl),intent(out)             :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)              :: iflag   !! status flag (`-1`=error, `0`=root found)
    real(dl),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(dl),intent(in),optional     :: fbx     !! if `f(bx)` is already known, it can be input here
    real(dl), parameter :: one = 1._dl, zero = 0._dl, two =2._dl, three = 3._dl
    real(dl),parameter :: eps   = epsilon(one)  !! original code had d1mach(4)
    real(dl) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    !convert the class function (un-type-checked) into correct type to call correctly for class argument
    call C_F_PROCPOINTER(c_funloc(func), f)

    tol1 = eps+one

    a=ax
    b=bx

    if (present(fax)) then
        fa = fax
    else
        fa=f(obj,a)
    end if
    if (present(fbx)) then
        fb = fbx
    else
        fb=f(obj,b)
    end if

    !check trivial cases first:
    if (fa==zero) then

        iflag = 0
        xzero = a
        fzero = fa

    elseif (fb==zero) then

        iflag = 0
        xzero = b
        fzero = fb

    elseif (fa*(fb/abs(fb))<zero) then  ! check that f(ax) and f(bx) have different signs

        c=a
        fc=fa
        d=b-a
        e=d

        do

            if (abs(fc)<abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end if

            tol1=two*eps*abs(b)+0.5_dl*tol
            xm = 0.5_dl*(c-b)
            if ((abs(xm)<=tol1).or.(fb==zero)) exit

            ! see if a bisection is forced
            if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
                s=fb/fa
                if (a/=c) then
                    ! inverse quadratic interpolation
                    q=fa/fc
                    r=fb/fc
                    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q=(q-one)*(r-one)*(s-one)
                else
                    ! linear interpolation
                    p=two*xm*s
                    q=one-s
                end if
                if (p<=zero) then
                    p=-p
                else
                    q=-q
                end if
                s=e
                e=d
                if (((two*p)>=(three*xm*q-abs(tol1*q))) .or. &
                    (p>=abs(0.5_dl*s*q))) then
                    d=xm
                    e=d
                else
                    d=p/q
                end if
            else
                d=xm
                e=d
            end if

            a=b
            fa=fb
            if (abs(d)<=tol1) then
                if (xm<=zero) then
                    b=b-tol1
                else
                    b=b+tol1
                end if
            else
                b=b+d
            end if
            fb=f(obj,b)
            if ((fb*(fc/abs(fc)))>zero) then
                c=a
                fc=fa
                d=b-a
                e=d
            end if

        end do

        iflag = 0
        xzero = b
        fzero = fb

    else

        iflag = -1
        write(error_unit,'(A)')&
            'Error in zeroin: f(ax) and f(bx) do not have different signs.'

    end if

    end subroutine brentq

    function Newton_Raphson2(xxl,xxh,funcs, param, param2) result(xm)
    use Precision
    implicit none
    real(dl), intent(in) :: xxl     ! root bracket 1
    real(dl), intent(in) :: xxh     ! root bracket 2
    real(dl)  :: xl,xh, xm     ! root
    external funcs        ! subroutine for non-linear equation
    real(dl), intent(in) :: param, param2 !parameters for function
    integer  :: k                      ! iteration count
    real(dl) :: xn, f,f2,df, error
    real(dl), parameter :: half=0.5_dl
    integer, parameter :: ITERMAX=1000 ! max number of iteration
    real(dl), parameter :: tol=1.e-8_dl ! tolerance for error

    xl =xxl
    xh = xxh
    call funcs(f,df,xl, param, param2)  ! Set xm=f(xl)
    call funcs(f2,df,xh, param, param2)  ! Set xn=f(xh)
    if (f*f2 > 0._dl) then           ! check if function changes sign
        error stop 'Newton_Raphson: root is not bracketed'
    endif
    if (f > 0._dl) then               ! Rearrange so that f(xl)< 0.d0 < f(xh)
        xm = xl
        xl = xh
        xh = xm
    endif

    error = abs(xh-xl)                ! error is width of bracketing interval
    xm = half*(xl+xh)                 ! Initialize guess for root
    k = 0                             ! initialize iteration count
    do while (error > tol .and. k < ITERMAX) ! iterate
        k = k+1                         ! increment iteration count
        call funcs(f,df,xm, param, param2) ! calculate f(xm), df(xm)
        if (f > 0._dl) then              ! Update root bracketing
            xh = xm                       ! update high
        else
            xl = xm                       ! update low
        endif
        xn = xm - f/df                  ! Tentative newton-Raphson step
        if ( (xn-xl)*(xn-xh) > 0._dl ) then ! check if new root falls within bracket
            xm = half* (xh+xl)            ! if no use a Bisection step
            error = abs(xh-xl)            ! error is width of interval
        else
            error = abs(xn-xm)            ! if within bracket: error is change in root
            xm = xn                       ! update successful Newton-Raphson step
        endif
    enddo

    if (error > tol) then       ! Check if solution converged
        write(*,*) 'Newton_Raphson:solution did not converge, xn, funcs(xn),D(xn)'
        write(*,*) xn, f, error
    endif

    end function Newton_Raphson2


    subroutine Gauss_Legendre(x, w, n)
    ! Get n-point Gauss-Legendre nodes x and weights w on [-1, 1].

    use constants, only : dl, const_pi
    implicit none

    integer,  intent(in)  :: n
    real(dl), intent(out) :: x(n), w(n)

    real(dl), parameter :: tol = 1.0e-15_dl
    integer,  parameter :: max_iter = 50

    integer  :: i, j, k, m, iter
    real(dl) :: p1, p2, p3, pp, z, dz, wi

    if (n < 1) error stop "Gauss_Legendre: n must be positive"

    m = (n + 1) / 2

!$omp parallel do default(none) schedule(static) if(n >= 128) &
!$omp& private(i, j, k, iter, p1, p2, p3, pp, z, dz, wi) &
!$omp& shared(x, w, n, m)
    do i = 1, m

        z = cos(const_pi * (real(i, dl) - 0.25_dl) / &
            (real(n, dl) + 0.5_dl))

        do iter = 1, max_iter
            p1 = 1._dl
            p2 = 0._dl

            do j = 1, n
                p3 = p2
                p2 = p1
                p1 = ((2*j - 1) * z * p2 - (j - 1) * p3) / j
            end do

            pp = real(n, dl) * (z*p1 - p2) / (z*z - 1._dl)
            dz = p1 / pp
            z  = z - dz

            if (abs(dz) <= tol) exit
        end do

        if (iter > max_iter) error stop "Gauss_Legendre: Newton iteration failed"

        k  = n + 1 - i
        wi = 2._dl / ((1._dl - z*z) * pp*pp)

        x(i) = -z
        x(k) =  z
        w(i) = wi
        w(k) = wi

        if (i == k) x(i) = 0._dl
    end do
!$omp end parallel do

    end subroutine Gauss_Legendre

    subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input parameters.
    !
    ! Generates the set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3) for all allowed l1
    ! from l1min = max(abs(l2-l3),abs(m1)) to l1max = l2+l3, with m1 = -(m2+m3).
    ! The resulting 3j-coeffs are stored as thrcof(l1-l1min+1).
    !
    ! For numerical stability the recursion proceeds simultaneously forwards from
    ! l1min and backwards from l1max, with the two branches matched over three
    ! overlapping points where the forward recursion stops being stable.
    use MpiUtils, only : MpiStop
    implicit none
    integer, intent(in) :: l2in, l3in, m2in, m3in
    real(dl), intent(out) :: thrcof(*)
    integer, parameter :: i8 = selected_int_kind(18)
    integer(i8) :: l2, l3, m2, m3
    integer(i8) :: l1, m1, l1min, l1max, nfin, a1, a2, dv, dv0, m32

    real(dl) :: newfac, oldfac, sumfor, c1, c2, c1old, denom, x, sum1, sumuni
    real(dl) :: x1, x2, x3, y, y1, y2, y3, sum2, sumbac, ratio, cnorm, thresh
    integer :: i, idx, nlim, lstep, nstep2, n
    logical :: forward_done
    real(dl), parameter :: zero = 0._dl, one = 1._dl
    real(dl), parameter :: tiny_val = 1.0e-30_dl, srtiny = 1.0e-15_dl
    real(dl), parameter :: huge_val = 1.e30_dl, srhuge = 1.e15_dl

    l2 = l2in
    l3 = l3in
    m2 = m2in
    m3 = m3in
    newfac = 0
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    if (l2 < abs(m2) .or. l3 < m3) then
        call MpiStop('GetThreeJs: invalid input (l2 < |m2| or l3 < m3)')
        return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3), abs(m1))
    l1max = l2+l3

    if (l1min >= l1max) then
        if (l1min /= l1max) then
            call MpiStop('GetThreeJs: invalid input (no allowed l1)')
            return
        end if

        ! reached if l1 can take only one value, i.e. l1min = l1max
        thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
        return

    end if

    nfin = l1max-l1min+1
    m32 = m3-m2
    dv0 = (l3*(l3+1) - l2*(l2+1))*m1

    ! forward recursion from l1min
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1+1)*tiny_val
    lstep = 1
    forward_done = .false.

    forward: do
        lstep = lstep+1
        l1 = l1+1

        oldfac = newfac
        a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
        a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
        newfac = sqrt(a2*real(a1,dl))
        if (l1 == 1) then
            ! if l1 = 1, (l1-1) has to be factored out of dv, hence
            c1 = -(2*l1-1)*l1*m32/newfac
        else
            dv = dv0 + l1*(l1-1)*m32
            denom = (l1-1)*newfac

            if (lstep > 2) c1old = abs(c1)
            c1 = -(2*l1-1)*real(dv,dl)/denom
        end if

        if (lstep <= 2) then
            ! if l1 = l1min+1 the third term in the recursion eqn vanishes, hence
            x = srtiny*c1
            thrcof(2) = x
            sum1 = sum1 + tiny_val*(2*l1+1)*c1*c1
            if (lstep == nfin) then
                forward_done = .true.
                exit forward
            end if
            cycle forward
        end if

        c2 = -l1*oldfac/denom

        ! recursion to the next 3j-coeff x
        x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
        thrcof(lstep) = x
        sumfor = sum1
        sum1 = sum1 + (2*l1+1)*x*x
        if (lstep == nfin) exit forward

        ! rescale if the last 3j-coeff exceeds srhuge, to prevent overflow
        if (abs(x) >= srhuge) then
            do i = 1, lstep
                if (abs(thrcof(i)) < srtiny) thrcof(i) = zero
                thrcof(i) = thrcof(i)/srhuge
            end do
            sum1 = sum1/huge_val
            sumfor = sumfor/huge_val
            x = x/srhuge
        end if

        ! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
        ! 3j-values and so is numerically stable. Once an increase of abs(c1) is
        ! detected, the recursion direction is reversed.
        if (c1old <= abs(c1)) exit forward
    end do forward

    if (forward_done) then
        ! forward recursion already covered the whole range (nfin = 2)
        sumuni = sum1
    else
        ! keep three 3j-coeffs around the matching point for comparison with
        ! the backward recursion
        x1 = x
        x2 = thrcof(lstep-1)
        x3 = thrcof(lstep-2)
        nstep2 = nfin-lstep+3

        ! backward recursion from l1max taking nstep2 steps, so that
        ! forward and backward recursion overlap at the three points
        ! l1 = lmatch-1, lmatch, lmatch+1
        l1 = l1max
        thrcof(nfin) = srtiny
        sum2 = tiny_val*(2*l1+1)

        l1 = l1+2
        lstep = 1

        backward: do
            lstep = lstep+1
            l1 = l1-1

            oldfac = newfac
            a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
            a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
            newfac = sqrt(a1*real(a2,dl))

            dv = dv0 + l1*(l1-1)*m32
            denom = l1*newfac
            c1 = -(2*l1-1)*real(dv,dl)/denom

            if (lstep <= 2) then
                ! if l1 = l1max-1 the third term in the recursion vanishes
                y = srtiny*c1
                thrcof(nfin-1) = y
                sumbac = sum2
                sum2 = sum2 + tiny_val*(2*l1-3)*c1*c1
                cycle backward
            end if

            c2 = -(l1-1)*oldfac/denom

            ! recursion to the next 3j-coeff y
            y = c1*thrcof(nfin+2-lstep) + c2*thrcof(nfin+3-lstep)

            if (lstep == nstep2) exit backward

            thrcof(nfin+1-lstep) = y
            sumbac = sum2
            sum2 = sum2 + (2*l1-3)*y*y

            ! rescale if the last 3j-coeff exceeds srhuge, to prevent overflow
            if (abs(y) >= srhuge) then
                do i = 1, lstep
                    idx = nfin-i+1
                    if (abs(thrcof(idx)) < srtiny) thrcof(idx) = zero
                    thrcof(idx) = thrcof(idx)/srhuge
                end do
                sum2 = sum2/huge_val
                sumbac = sumbac/huge_val
            end if
        end do backward

        ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the
        ! corresponding backward recursion values y1, y2, y3, determining the
        ! ratio such that yi = ratio*xi (i=1,2,3) holds with minimal error
        y3 = y
        y2 = thrcof(nfin+2-lstep)
        y1 = thrcof(nfin+3-lstep)

        ratio = (x1*y1 + x2*y2 + x3*y3)/(x1*x1 + x2*x2 + x3*x3)
        nlim = nfin-nstep2+1

        if (abs(ratio) >= 1) then
            thrcof(1:nlim) = ratio*thrcof(1:nlim)
            sumuni = ratio*ratio*sumfor + sumbac
        else
            ratio = 1/ratio
            thrcof(nlim+1:nfin) = ratio*thrcof(nlim+1:nfin)
            sumuni = sumfor + ratio*ratio*sumbac
        end if
    end if

    ! normalise 3j-coeffs; sign convention for the last coeff fixes the overall phase
    cnorm = 1/sqrt(sumuni)
    if (sign(one,thrcof(nfin))*(-1)**abs(l2+m2-l3+m3) <= 0) cnorm = -cnorm

    if (abs(cnorm) >= one) then
        thrcof(1:nfin) = cnorm*thrcof(1:nfin)
    else
        thresh = tiny_val/abs(cnorm)
        do n = 1, int(nfin)
            if (abs(thrcof(n)) < thresh) thrcof(n) = zero
            thrcof(n) = cnorm*thrcof(n)
        end do
    end if

    end subroutine GetThreeJs


    function GetChiSquared(c_inv, Y, n) result(chi2)
    !get dot_product(matmul(C_inv,Y), Y) efficiently assuming c_inv symmetric
    integer, intent(in) :: n
    real(dl), intent(in) :: Y(n)
    real(dl), intent(in) :: c_inv(n,n)
    integer j
    real(dl) ztemp, chi2

    chi2 = 0
    if (n>=512) then
        !$OMP parallel do private(j,ztemp) reduction(+:chi2) schedule(static,16)
        do  j = 1, n
            ztemp= dot_product(Y(j+1:n), c_inv(j+1:n, j))
            chi2=chi2+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
        end do
    else
        do  j = 1, n
            ztemp= dot_product(Y(j+1:n), c_inv(j+1:n, j))
            chi2=chi2+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
        end do
    end if

    end function GetChiSquared

    subroutine integrate_3j(W,lmax_w, n, dopol, M, lmax)
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    !Get coupling matrix, eg for pesudo-CL
    integer, intent(in) :: lmax, lmax_w, n
    real(dl), intent(in) :: W(0:lmax_w,n)
    logical, intent(in) :: dopol
    real(dl), intent(out) :: M(0:lmax,0:lmax, n)
    integer l1, l2, lplus, lminus, thread_ix, ix
    real(dl), allocatable :: threejj0(:,:), threejj2(:,:)

    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()

    allocate(threejj0(0:2*lmax,thread_ix))
    if (dopol) then
        allocate(threejj2(0:2*lmax,thread_ix))
    end if

    !$OMP parallel do private(l1,l2,lminus,lplus,thread_ix,ix), schedule(dynamic)
    do l1 = 0, lmax
        thread_ix =1
        !$ thread_ix = OMP_GET_THREAD_NUM()+1
        do l2 = 0, l1
            lplus =  min(lmax_w,l1+l2)
            lminus = abs(l1-l2)

            call GetThreeJs(threejj0(lminus:,thread_ix),l1,l2,0,0)

            if (dopol .and. l1>=2 .and. l2>=2) then
                !note that lminus is correct, want max(abs(l1-l2),abs(m1)) where m1=0 here
                !(polarization coupling depends on lowest multipoles of the mask)
                call GetThreeJs(threejj2(lminus:,thread_ix),l1,l2,-2,2)
                M(l2,l1,2) = sum(W(lminus:lplus:2,2)*threejj0(lminus:lplus:2,thread_ix) &
                    *threejj2(lminus:lplus:2,thread_ix)) !TE
                M(l2,l1,3) = sum(W(lminus:lplus:2,3)*threejj2(lminus:lplus:2,thread_ix)**2) !EE
                M(l2,l1,4) = sum(W(lminus+1:lplus:2,3)*threejj2(lminus+1:lplus:2,thread_ix)**2) !EB
            end if
            if (n>1 .and. .not. dopol) then
                threejj0(lminus:lplus,thread_ix) = threejj0(lminus:lplus,thread_ix)**2
                do ix=1,n
                    M(l2,l1,ix) = sum(W(lminus:lplus,ix)* threejj0(lminus:lplus,thread_ix))
                end do
            else
                M(l2,l1,1) = sum(W(lminus:lplus,1)* threejj0(lminus:lplus,thread_ix)**2)
            end if
        end do
    end do

    do l1=0, lmax
        do l2 = l1+1,lmax
            M(l2,l1,:) = M(l1,l2,:)
        end do
    end do
    end subroutine integrate_3j

    end module MathUtils
