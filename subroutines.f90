    !General numerical routines and global accuracy. Includes modified dverk for CAMB.


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine splder(y,dy,n, g)
    use Precision
    !  Splder fits a cubic spline to y and returns the first derivatives at
    !  the grid points in dy.  Dy is equivalent to a 4th-order Pade
    !  difference formula for dy/di.
    implicit none
    integer, intent(in) :: n
    real(dl), intent(in) :: y(n),g(n)
    real(dl), intent(out) :: dy(n)
    integer :: n1, i
    real(dl), allocatable, dimension(:) :: f

    allocate(f(n))
    n1=n-1
    !  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
    f(1)=(-10._dl*y(1)+15._dl*y(2)-6._dl*y(3)+y(4))/6._dl
    f(n)=(10._dl*y(n)-15._dl*y(n1)+6._dl*y(n-2)-y(n-3))/6._dl
    !  Solve the tridiagonal system
    !  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
    !  with dy(1)=f(1), dy(n)=f(n).
    do i=2,n1
        f(i)=g(i)*(3._dl*(y(i+1)-y(i-1))-f(i-1))
    end do
    dy(n)=f(n)
    do i=n1,1,-1
        dy(i)=f(i)-g(i)*dy(i+1)
    end do
    deallocate(f)
    end subroutine splder
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine splini(g,n)
    use Precision
    !  Splini must be called before splder to initialize array g in common.
    implicit none
    integer, intent(in) :: n
    real(dl), intent(out):: g(n)
    integer :: i

    g(1)=0._dl
    do i=2,n
        g(i)=1/(4._dl-g(i-1))
    end do
    end subroutine splini


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rombint2(f,a,b,tol, maxit, minsteps)
    use precision
    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.

    ! Modified by AL to specify max iterations and minimum number of steps
    ! (min steps useful to stop wrong results on periodic or sharp functions)
    implicit none
    integer, parameter :: MAXITER=20,MAXJ=5
    dimension g(MAXJ+1)
    real(dl) f
    external f
    real(dl) :: rombint2
    real(dl), intent(in) :: a,b,tol
    integer, intent(in):: maxit,minsteps

    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g, g0, g1, fourj

    h=0.5d0*(b-a)
    gmax=h*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
    do
        i=i+1
        if (i > maxit.or.(i > 5.and.abs(error) < tol) .and. nint > minsteps) exit
        !  Calculate next trapezoidal rule approximation to integral.
        g0=0._dl
        do k=1,nint
            g0=g0+f(a+(k+k-1)*h)
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
        if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
        else
            error=gmax
        end if
        gmax=g0
        g(jmax+1)=g0
    end do

    rombint2=g0
    if (i > maxit .and. abs(error) > tol)  then
        write(*,*) 'Warning: Rombint2 failed to converge; '
        write (*,*)'integral, error, tol:', rombint2,error, tol
    end if

    end function rombint2

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rombint(f,a,b,tol)
    use Precision
    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.
    !
    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    dimension g(MAXJ+1)
    real(dl) f
    external f
    real(dl) :: rombint
    real(dl), intent(in) :: a,b,tol
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g, g0, g1, fourj
    !

    h=0.5d0*(b-a)
    gmax=h*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do 20 k=1,nint
        g0=g0+f(a+(k+k-1)*h)
20  continue
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4._dl*fourj
        g1=g0+(g0-g(j))/(fourj-1._dl)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1._dl-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  rombint=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: Rombint failed to converge; '
        write (*,*)'integral, error, tol:', rombint,error, tol
    end if

    end function rombint


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rombint_obj(obj,f,a,b,tol, maxit)
    use Precision
    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.
    !
    implicit none
    integer, intent(in), optional :: maxit
    integer :: MAXITER=20
    integer, parameter :: MAXJ=5
    dimension g(MAXJ+1)
    real obj !dummy
    real(dl) f
    external f
    real(dl) :: rombint_obj
    real(dl), intent(in) :: a,b,tol
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g, g0, g1, fourj
    !

    if (present(maxit)) then
        MaxIter = maxit
    end if
    h=0.5d0*(b-a)
    gmax=h*(f(obj,a)+f(obj,b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do 20 k=1,nint
        g0=g0+f(obj,a+(k+k-1)*h)
20  continue
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4._dl*fourj
        g1=g0+(g0-g(j))/(fourj-1._dl)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1._dl-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  rombint_obj=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: Rombint failed to converge; '
        write (*,*)'integral, error, tol:', rombint_obj,error, tol
    end if

    end function rombint_obj


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! calculates array of second derivatives used by cubic spline
    ! interpolation. y2 is array of second derivatives, yp1 and ypn are first
    ! derivatives at end points.

    !Thanks Martin Reinecke
    subroutine spline(x,y,n,d11,d1n,d2)
    use Precision
    integer, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n), d11, d1n
    real(dl), intent(out) :: d2(n)
    integer i
    real(dl) xp,qn,sig,un,xxdiv,u(n-1),d1l,d1r

    d1r= (y(2)-y(1))/(x(2)-x(1))
    if (d11>.99e30_dl) then
        d2(1)=0._dl
        u(1)=0._dl
    else
        d2(1)=-0.5_dl
        u(1)=(3._dl/(x(2)-x(1)))*(d1r-d11)
    endif

    do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1._dl/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1._dl/(sig*d2(i-1)+2._dl)

        d2(i)=(sig-1._dl)*xp

        u(i)=(6._dl*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
    end do
    d1l=d1r

    if (d1n>.99e30_dl) then
        qn=0._dl
        un=0._dl
    else
        qn=0.5_dl
        un=(3._dl/(x(n)-x(n-1)))*(d1n-d1l)
    endif

    d2(n)=(un-qn*u(n-1))/(qn*d2(n-1)+1._dl)
    do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
    end do
    end subroutine spline

    SUBROUTINE spline_deriv(x,y,y2,y1,n)
    !Get derivative y1 given array of x, y and y''
    use Precision
    implicit none
    INTEGER, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n), y2(n)
    real(dl), intent(out) :: y1(n)
    INTEGER i
    real(dl) dx

    do i=1, n-1

        dx = (x(i+1) - x(i))
        y1(i) = (y(i+1) - y(i))/dx - dx*(2*y2(i) + y2(i+1))/6
    end do
    dx = x(n) - x(n-1)
    y1(n) = (y(n) - y(n-1))/dx + dx* ( y2(i-1)  + 2*y2(i) )/6

    END SUBROUTINE spline_deriv

    subroutine spline_integrate(x,y,y2,yint,n)
    !Cumulative integral of cubic spline
    use Precision
    integer, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n), y2(n)
    real(dl), intent(out) :: yint(n)
    real(dl) dx
    integer i

    yint(1) = 0
    do i=2, n

        dx = (x(i) - x(i-1))
        yint(i) = yint(i-1) + dx*( (y(i)+y(i-1))/2 - dx**2/24*(y2(i)+y2(i-1)))

    end do

    end subroutine spline_integrate



    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  this is not the splint given in numerical recipes


    subroutine splint(y,z,n)
    use Precision
    !  Splint integrates a cubic spline, providing the output value
    !  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
    !  to y(i).
    !
    implicit none
    integer, intent(in) :: n
    real(dl), intent(in) :: y(n)
    real(dl), intent(out) :: z

    integer :: n1
    real(dl) :: dy1, dyn
    !
    n1=n-1
    !  Cubic fit to dy/di at boundaries.
    !       dy1=(-11._dl*y(1)+18._dl*y(2)-9._dl*y(3)+2._dl*y(4))/6._dl
    dy1=0._dl
    dyn=(11._dl*y(n)-18._dl*y(n1)+9._dl*y(n-2)-2._dl*y(n-3))/6._dl
    !
    z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12._dl
    z= z + sum(y(2:n1))
    end subroutine splint


    !This version is modified to pass an object parameter to the function on each call
    !Fortunately Fortran doesn't do type checking on functions, so we can pretend the
    !passed object parameter (EV) is any type we like. In reality it is just a pointer.
    subroutine dverk (EV,n, fcn, x, y, xend, tol, ind, c, nw, w)
    use Precision
    use AMLUtils
    integer n, ind, nw, k
    real(dl) x, y(n), xend, tol, c(*), w(nw,9), temp
    real EV !It isn't, but as long as it maintains it as a pointer we are OK
    !     class(*) EV !seems to be correct way to do this in Fortran 2003
    !
    !***********************************************************************
    !                                                                      *
    ! note added 11/14/85.                                                 *
    !                                                                      *
    ! if you discover any errors in this subroutine, please contact        *
    !                                                                      *
    !        kenneth r. jackson                                            *
    !        department of computer science                                *
    !        university of toronto                                         *
    !        toronto, ontario,                                             *
    !        canada   m5s 1a4                                              *
    !                                                                      *
    !        phone: 416-978-7075                                           *
    !                                                                      *
    !        electronic mail:                                              *
    !        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
    !        csnet:  krj@toronto                                           *
    !        arpa:   krj.toronto@csnet-relay                               *
    !        bitnet: krj%toronto@csnet-relay.arpa                          *
    !                                                                      *
    ! dverk is written in fortran 66.                                      *
    !                                                                      *
    ! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
    ! set for a  vax  in  double  precision.  they  should  be  reset,  as *
    ! described below, if this program is run on another machine.          *
    !                                                                      *
    ! the c array is declared in this subroutine to have one element only, *
    ! although  more  elements  are  referenced  in this subroutine.  this *
    ! causes some compilers to issue warning messages.  there is,  though, *
    ! no  error  provided  c is declared sufficiently large in the calling *
    ! program, as described below.                                         *
    !                                                                      *
    ! the following external statement  for  fcn  was  added  to  avoid  a *
    ! warning  message  from  the  unix  f77 compiler.  the original dverk *
    ! comments and code follow it.                                         *
    !                                                                      *
    !***********************************************************************
    !
    external fcn
    !
    !***********************************************************************
    !                                                                      *
    !     purpose - this is a runge-kutta  subroutine  based  on  verner's *
    ! fifth and sixth order pair of formulas for finding approximations to *
    ! the solution of  a  system  of  first  order  ordinary  differential *
    ! equations  with  initial  conditions. it attempts to keep the global *
    ! error proportional to  a  tolerance  specified  by  the  user.  (the *
    ! proportionality  depends  on the kind of error control that is used, *
    ! as well as the differential equation and the range of integration.)  *
    !                                                                      *
    !     various options are available to the user,  including  different *
    ! kinds  of  error control, restrictions on step sizes, and interrupts *
    ! which permit the user to examine the state of the  calculation  (and *
    ! perhaps make modifications) during intermediate stages.              *
    !                                                                      *
    !     the program is efficient for non-stiff systems.  however, a good *
    ! variable-order-adams  method  will probably be more efficient if the *
    ! function evaluations are very costly.  such a method would  also  be *
    ! more suitable if one wanted to obtain a large number of intermediate *
    ! solution values by interpolation, as might be the case  for  example *
    ! with graphical output.                                               *
    !                                                                      *
    !                                    hull-enright-jackson   1/10/76    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     use - the user must specify each of the following                *
    !                                                                      *
    !     n  number of equations                                           *
    !                                                                      *
    !   fcn  name of subroutine for evaluating functions - the  subroutine *
    !           itself must also be provided by the user - it should be of *
    !           the following form                                         *
    !              subroutine fcn(n, x, y, yprime)                         *
    !              integer n                                               *
    !              real(dl) x, y(n), yprime(n)                     *
    !                      *** etc ***                                     *
    !           and it should evaluate yprime, given n, x and y            *
    !                                                                      *
    !     x  independent variable - initial value supplied by user         *
    !                                                                      *
    !     y  dependent variable - initial values of components y(1), y(2), *
    !           ..., y(n) supplied by user                                 *
    !                                                                      *
    !  xend  value of x to which integration is to be carried out - it may *
    !           be less than the initial value of x                        *
    !                                                                      *
    !   tol  tolerance - the subroutine attempts to control a norm of  the *
    !           local  error  in  such  a  way  that  the  global error is *
    !           proportional to tol. in some problems there will be enough *
    !           damping  of  errors, as well as some cancellation, so that *
    !           the global error will be less than tol. alternatively, the *
    !           control   can   be  viewed  as  attempting  to  provide  a *
    !           calculated value of y at xend which is the exact  solution *
    !           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
    !           is proportional to tol.  (the norm  is  a  max  norm  with *
    !           weights  that  depend on the error control strategy chosen *
    !           by the user.  the default weight for the k-th component is *
    !           1/max(1,abs(y(k))),  which therefore provides a mixture of *
    !           absolute and relative error control.)                      *
    !                                                                      *
    !   ind  indicator - on initial entry ind must be set equal to  either *
    !           1  or  2. if the user does not wish to use any options, he *
    !           should set ind to 1 - all that remains for the user to  do *
    !           then  is  to  declare c and w, and to specify nw. the user *
    !           may also  select  various  options  on  initial  entry  by *
    !           setting ind = 2 and initializing the first 9 components of *
    !           c as described in the next section.  he may also  re-enter *
    !           the  subroutine  with ind = 3 as mentioned again below. in *
    !           any event, the subroutine returns with ind equal to        *
    !              3 after a normal return                                 *
    !              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
    !              -1, -2, or -3 after an error condition (see below)      *
    !                                                                      *
    !     c  communications vector - the dimension must be greater than or *
    !           equal to 24, unless option c(1) = 4 or 5 is used, in which *
    !           case the dimension must be greater than or equal to n+30   *
    !                                                                      *
    !    nw  first dimension of workspace w -  must  be  greater  than  or *
    !           equal to n                                                 *
    !                                                                      *
    !     w  workspace matrix - first dimension must be nw and second must *
    !           be greater than or equal to 9                              *
    !                                                                      *
    !     the subroutine  will  normally  return  with  ind  =  3,  having *
    ! replaced the initial values of x and y with, respectively, the value *
    ! of xend and an approximation to y at xend.  the  subroutine  can  be *
    ! called  repeatedly  with new values of xend without having to change *
    ! any other argument.  however, changes in tol, or any of the  options *
    ! described below, may also be made on such a re-entry if desired.     *
    !                                                                      *
    !     three error returns are also possible, in which  case  x  and  y *
    ! will be the most recently accepted values -                          *
    !     with ind = -3 the subroutine was unable  to  satisfy  the  error *
    !        requirement  with a particular step-size that is less than or *
    !        equal to hmin, which may mean that tol is too small           *
    !     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
    !        probably  means  that the requested tol (which is used in the *
    !        calculation of hmin) is too small                             *
    !     with ind = -1 the allowed maximum number of fcn evaluations  has *
    !        been  exceeded,  but  this  can only occur if option c(7), as *
    !        described in the next section, has been used                  *
    !                                                                      *
    !     there are several circumstances that will cause the calculations *
    ! to  be  terminated,  along with output of information that will help *
    ! the user determine the cause of  the  trouble.  these  circumstances *
    ! involve  entry with illegal or inconsistent values of the arguments, *
    ! such as attempting a normal  re-entry  without  first  changing  the *
    ! value of xend, or attempting to re-enter with ind less than zero.    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     options - if the subroutine is entered with ind = 1, the first 9 *
    ! components of the communications vector are initialized to zero, and *
    ! the subroutine uses only default values  for  each  option.  if  the *
    ! subroutine  is  entered  with ind = 2, the user must specify each of *
    ! these 9 components - normally he would first set them all  to  zero, *
    ! and  then  make  non-zero  those  that  correspond to the particular *
    ! options he wishes to select. in any event, options may be changed on *
    ! re-entry  to  the  subroutine  -  but if the user changes any of the *
    ! options, or tol, in the course of a calculation he should be careful *
    ! about  how  such changes affect the subroutine - it may be better to *
    ! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
    ! program  -  the information is available to the user, but should not *
    ! normally be changed by him.)                                         *
    !                                                                      *
    !  c(1)  error control indicator - the norm of the local error is  the *
    !           max  norm  of  the  weighted  error  estimate  vector, the *
    !           weights being determined according to the value of c(1) -  *
    !              if c(1)=1 the weights are 1 (absolute error control)    *
    !              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
    !                 control)                                             *
    !              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
    !                 (relative  error  control,  unless abs(y(k)) is less *
    !                 than the floor value, abs(c(2)) )                    *
    !              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
    !                 (here individual floor values are used)              *
    !              if c(1)=5 the weights are 1/abs(c(k+30))                *
    !              for all other values of c(1), including  c(1) = 0,  the *
    !                 default  values  of  the  weights  are  taken  to be *
    !                 1/max(1,abs(y(k))), as mentioned earlier             *
    !           (in the two cases c(1) = 4 or 5 the user must declare  the *
    !           dimension of c to be at least n+30 and must initialize the *
    !           components c(31), c(32), ..., c(n+30).)                    *
    !                                                                      *
    !  c(2)  floor value - used when the indicator c(1) has the value 3    *
    !                                                                      *
    !  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
    !           to be abs(c(3)) - otherwise it uses the default value      *
    !              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
    !           where dwarf is a very small positive  machine  number  and *
    !           rreb is the relative roundoff error bound                  *
    !                                                                      *
    !  c(4)  hstart specification - if not zero, the subroutine  will  use *
    !           an  initial  hmag equal to abs(c(4)), except of course for *
    !           the restrictions imposed by hmin and hmax  -  otherwise it *
    !           uses the default value of hmax*(tol)**(1/6)                *
    !                                                                      *
    !  c(5)  scale specification - this is intended to be a measure of the *
    !           scale of the problem - larger values of scale tend to make *
    !           the method more reliable, first  by  possibly  restricting *
    !           hmax  (as  described  below) and second, by tightening the *
    !           acceptance requirement - if c(5) is zero, a default  value *
    !           of  1  is  used.  for  linear  homogeneous  problems  with *
    !           constant coefficients, an appropriate value for scale is a *
    !           norm  of  the  associated  matrix.  for other problems, an *
    !           approximation to  an  average  value  of  a  norm  of  the *
    !           jacobian along the trajectory may be appropriate           *
    !                                                                      *
    !  c(6)  hmax specification - four cases are possible                  *
    !           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
    !              min(abs(c(6)),2/abs(c(5)))                              *
    !           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
    !           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
    !              2/abs(c(5))                                             *
    !           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
    !              of 2                                                    *
    !                                                                      *
    !  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
    !           error  return with ind = -1 will be caused when the number *
    !           of function evaluations exceeds abs(c(7))                  *
    !                                                                      *
    !  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  after  it  has  chosen  its *
    !           preliminary value of hmag, and just before choosing htrial *
    !           and  xtrial  in  preparation for taking a step (htrial may *
    !           differ from hmag in sign, and may  require  adjustment  if *
    !           xend  is  near) - the subroutine returns with ind = 4, and *
    !           will resume calculation at the point  of  interruption  if *
    !           re-entered with ind = 4                                    *
    !                                                                      *
    !  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  immediately  after  it  has *
    !           decided whether or not to accept the result  of  the  most *
    !           recent  trial step, with ind = 5 if it plans to accept, or *
    !           ind = 6 if it plans to reject -  y(*)  is  the  previously *
    !           accepted  result, while w(*,9) is the newly computed trial *
    !           value, and w(*,2) is the unweighted error estimate vector. *
    !           the  subroutine  will  resume calculations at the point of *
    !           interruption on re-entry with ind = 5 or 6. (the user  may *
    !           change ind in this case if he wishes, for example to force *
    !           acceptance of a step that would otherwise be rejected,  or *
    !           vice versa. he can also restart with ind = 1 or 2.)        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  summary of the components of the communications vector              *
    !                                                                      *
    !     prescribed at the option       determined by the program         *
    !           of the user                                                *
    !                                                                      *
    !                                    c(10) rreb(rel roundoff err bnd)  *
    !     c(1) error control indicator   c(11) dwarf (very small mach no)  *
    !     c(2) floor value               c(12) weighted norm y             *
    !     c(3) hmin specification        c(13) hmin                        *
    !     c(4) hstart specification      c(14) hmag                        *
    !     c(5) scale specification       c(15) scale                       *
    !     c(6) hmax specification        c(16) hmax                        *
    !     c(7) max no of fcn evals       c(17) xtrial                      *
    !     c(8) interrupt no 1            c(18) htrial                      *
    !     c(9) interrupt no 2            c(19) est                         *
    !                                    c(20) previous xend               *
    !                                    c(21) flag for xend               *
    !                                    c(22) no of successful steps      *
    !                                    c(23) no of successive failures   *
    !                                    c(24) no of fcn evals             *
    !                                                                      *
    !  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  an overview of the program                                          *
    !                                                                      *
    !     begin initialization, parameter checking, interrupt re-entries   *
    !  ......abort if ind out of range 1 to 6                              *
    !  .     cases - initial entry, normal re-entry, interrupt re-entries  *
    !  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
    !  v........abort if n.gt.nw or tol.le.0                               *
    !  .        if initial entry without options (ind .eq. 1)              *
    !  .           set c(1) to c(9) equal to zero                          *
    !  .        else initial entry with options (ind .eq. 2)               *
    !  .           make c(1) to c(9) non-negative                          *
    !  .           make floor values non-negative if they are to be used   *
    !  .        end if                                                     *
    !  .        initialize rreb, dwarf, prev xend, flag, counts            *
    !  .     case 2 - normal re-entry (ind .eq. 3)                         *
    !  .........abort if xend reached, and either x changed or xend not    *
    !  .        re-initialize flag                                         *
    !  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
    !  v        transfer control to the appropriate re-entry point.......  *
    !  .     end cases                                                  .  *
    !  .  end initialization, etc.                                      .  *
    !  .                                                                v  *
    !  .  loop through the following 4 stages, once for each trial step .  *
    !  .     stage 1 - prepare                                          .  *
    !***********error return (with ind=-1) if no of fcn evals too great .  *
    !  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
    !  .        calc hmin, scale, hmax                                  .  *
    !***********error return (with ind=-2) if hmin .gt. hmax            .  *
    !  .        calc preliminary hmag                                   .  *
    !***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
    !  .        calc hmag, xtrial and htrial                            .  *
    !  .     end stage 1                                                .  *
    !  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
    !  .     stage 3 - calc the error estimate                          .  *
    !  .     stage 4 - make decisions                                   .  *
    !  .        set ind=5 if step acceptable, else set ind=6            .  *
    !***********interrupt no 2 if requested....................re-entry.v  *
    !  .        if step accepted (ind .eq. 5)                              *
    !  .           update x, y from xtrial, ytrial                         *
    !  .           add 1 to no of successful steps                         *
    !  .           set no of successive failures to zero                   *
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend *
    !  .        else step not accepted (ind .eq. 6)                        *
    !  .           add 1 to no of successive failures                      *
    !**************error return (with ind=-3) if hmag .le. hmin            *
    !  .        end if                                                     *
    !  .     end stage 4                                                   *
    !  .  end loop                                                         *
    !  .                                                                   *
    !  begin abort action                                                  *
    !     output appropriate  message  about  stopping  the  calculations, *
    !        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
    !        previous xend,  no of  successful  steps,  no  of  successive *
    !        failures, no of fcn evals, and the components of y            *
    !     stop                                                             *
    !  end abort action                                                    *
    !                                                                      *
    !***********************************************************************
    !
    !     ******************************************************************
    !     * begin initialization, parameter checking, interrupt re-entries *
    !     ******************************************************************
    !
    !  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) go to 500
    !
    !        cases - initial entry, normal re-entry, interrupt re-entries
    !         go to (5, 5, 45, 1111, 2222, 2222), ind
    if (ind==3) goto 45
    if (ind==4) goto 1111
    if (ind==5 .or. ind==6) goto 2222

    !        case 1 - initial entry (ind .eq. 1 or 2)
    !  .........abort if n.gt.nw or tol.le.0
    if (n.gt.nw .or. tol.le.0._dl) go to 500
    if (ind.eq. 2) go to 15
    !              initial entry without options (ind .eq. 1)
    !              set c(1) to c(9) equal to 0
    do k = 1, 9
        c(k) = 0._dl
    end do
    go to 35
15  continue
    !              initial entry with options (ind .eq. 2)
    !              make c(1) to c(9) non-negative
    do k = 1, 9
        c(k) = dabs(c(k))
    end do
    !              make floor values non-negative if they are to be used
    if (c(1).ne.4._dl .and. c(1).ne.5._dl) go to 30
    do k = 1, n
        c(k+30) = dabs(c(k+30))
    end do
30  continue
35  continue
    !           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2._dl**(-56)
    c(11) = 1.d-35
    !           set previous xend initially to initial value of x
    c(20) = x
    do k = 21, 24
        c(k) = 0._dl
    end do
    go to 50
    !        case 2 - normal re-entry (ind .eq. 3)
    !  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0._dl .and. &
        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
    !           re-initialize flag
    c(21) = 0._dl
    go to 50
    !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
    !           transfer control to the appropriate re-entry point..........
    !           this has already been handled by the computed go to        .
    !        end cases                                                     v
50  continue
    !
    !     end initialization, etc.
    !
    !     ******************************************************************
    !     * loop through the following 4 stages, once for each trial  step *
    !     * until the occurrence of one of the following                   *
    !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
    !     *        stage 4                                                 *
    !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
    !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
    !     *        requested, in stage 1 or stage 4                        *
    !     ******************************************************************
    !
99999 continue
    !
    !        ***************************************************************
    !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
    !        * and some parameter  checking,  and  end  up  with  suitable *
    !        * values of hmag, xtrial and htrial in preparation for taking *
    !        * an integration step.                                        *
    !        ***************************************************************
    !
    !***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0._dl .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
    !
    !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call fcn(EV,n, x, y, w(1,1))
    c(24) = c(24) + 1._dl
105 continue
    !
    !           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0._dl) go to 165
    !              calculate default value of hmin
    !              first calculate weighted norm y - c(12) - as specified
    !              by the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 115
    !                 absolute error control - weights are 1
    do 110 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2._dl) go to 120
    !                 relative error control - weights are 1/dabs(y(k)) so
    !                 weighted norm y is 1
    c(12) = 1._dl
    go to 160
120 if (c(1) .ne. 3._dl) go to 130
    !                 weights are 1/max(c(2),abs(y(k)))
    do 125 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(2))
125 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
130 if (c(1) .ne. 4._dl) go to 140
    !                 weights are 1/max(c(k+30),abs(y(k)))
    do 135 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
135 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
140 if (c(1) .ne. 5._dl) go to 150
    !                 weights are 1/c(k+30)
    do 145 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
145 continue
    c(12) = temp
    go to 160
150 continue
    !                 default case - weights are 1/max(1,abs(y(k)))
    do 155 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
155 continue
    c(12) = dmin1(temp, 1._dl)
160 continue
    c(13) = 10._dl*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
    !
    !           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0._dl) c(15) = 1._dl
    !
    !           calculate hmax - consider 4 cases
    !           case 1 both hmax and scale prescribed
    if (c(6).ne.0._dl .and. c(5).ne.0._dl) &
        c(16) = dmin1(c(6), 2._dl/c(5))
    !           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0._dl .and. c(5).eq.0._dl) c(16) = c(6)
    !           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0._dl .and. c(5).ne.0._dl) c(16) = 2._dl/c(5)
    !           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0._dl .and. c(5).eq.0._dl) c(16) = 2._dl
    !
    !***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
    !
    !           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
    !           case 1 - initial entry - use prescribed value of hstart, if
    !              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0._dl) c(14) = c(16)*tol**(1._dl/6._dl)
    go to 185
175 if (c(23) .gt. 1._dl) go to 180
    !           case 2 - after a successful step, or at most  one  failure,
    !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
    !              overflow. then avoid reduction by more than half.
    temp = 2._dl*c(14)
    if (tol .lt. (2._dl/.9d0)**6*c(19)) &
        temp = .9d0*(tol/c(19))**(1._dl/6._dl)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
    !           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
    !
    !           check against hmax
    c(14) = dmin1(c(14), c(16))
    !
    !           check against hmin
    c(14) = dmax1(c(14), c(13))
    !
    !***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0._dl) go to 1111
    ind = 4
    return
    !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
    !
    !           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
    !              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
    !              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
    !
    !           calculate htrial
    c(18) = c(17) - x
    !
    !        end stage 1
    !
    !        ***************************************************************
    !        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
    !        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
    !        * stage 3. w(*,9) is temporary storage until finally it holds *
    !        * ytrial.                                                     *
    !        ***************************************************************
    !
    temp = c(18)/1398169080000._dl
    !
    do 200 k = 1, n
        w(k,9) = y(k) + temp*w(k,1)*233028180000._dl
200 continue
    call fcn(EV,n, x + c(18)/6._dl, w(1,9), w(1,2))
    !
    do 205 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*74569017600._dl &
            + w(k,2)*298276070400._dl  )
205 continue
    call fcn(EV,n, x + c(18)*(4._dl/15._dl), w(1,9), w(1,3))
    !
    do 210 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._dl &
            - w(k,2)*3728450880000._dl &
            + w(k,3)*3495422700000._dl )
210 continue
    call fcn(EV,n, x + c(18)*(2._dl/3._dl), w(1,9), w(1,4))
    !
    do 215 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._dl &
            + w(k,2)*12816549900000._dl &
            - w(k,3)*9284716546875._dl &
            + w(k,4)*1237962206250._dl )
215 continue
    call fcn(EV,n, x + c(18)*(5._dl/6._dl), w(1,9), w(1,5))
    !
    do 220 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._dl &
            - w(k,2)*11185352640000._dl &
            + w(k,3)*9172628850000._dl &
            - w(k,4)*427218330000._dl &
            + w(k,5)*482505408000._dl  )
220 continue
    call fcn(EV,n, x + c(18), w(1,9), w(1,6))
    !
    do 225 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*770204740536._dl &
            + w(k,2)*2311639545600._dl &
            - w(k,3)*1322092233000._dl &
            - w(k,4)*453006781920._dl &
            + w(k,5)*326875481856._dl  )
225 continue
    call fcn(EV,n, x + c(18)/15._dl, w(1,9), w(1,7))
    !
    do 230 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._dl &
            - w(k,2)*9754668000000._dl &
            + w(k,3)*7897110375000._dl &
            - w(k,4)*192082660000._dl &
            + w(k,5)*400298976000._dl &
            + w(k,7)*201586000000._dl  )
230 continue
    call fcn(EV,n, x + c(18), w(1,9), w(1,8))
    !
    !           calculate ytrial, the extrapolated approximation and store
    !              in w(*,9)
    do 235 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*104862681000._dl &
            + w(k,3)*545186250000._dl &
            + w(k,4)*446637345000._dl &
            + w(k,5)*188806464000._dl &
            + w(k,7)*15076875000._dl &
            + w(k,8)*97599465000._dl   )
235 continue
    !
    !           add 7 to the no of fcn evals
    c(24) = c(24) + 7._dl
    !
    !        end stage 2
    !
    !        ***************************************************************
    !        * stage 3 - calculate the error estimate est. first calculate *
    !        * the  unweighted  absolute  error  estimate vector (per unit *
    !        * step) for the unextrapolated approximation and store it  in *
    !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
    !        * specified by the error  control  indicator  c(1).  finally, *
    !        * modify  this result to produce est, the error estimate (per *
    !        * unit step) for the extrapolated approximation ytrial.       *
    !        ***************************************************************
    !
    !           calculate the unweighted absolute error estimate vector
    do 300 k = 1, n
        w(k,2) = (   w(k,1)*8738556750._dl &
            + w(k,3)*9735468750._dl &
            - w(k,4)*9709507500._dl &
            + w(k,5)*8582112000._dl &
            + w(k,6)*95329710000._dl &
            - w(k,7)*15076875000._dl &
            - w(k,8)*97599465000._dl)/1398169080000._dl
300 continue
    !
    !           calculate the weighted max norm of w(*,2) as specified by
    !           the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 310
    !              absolute error control
    do 305 k = 1, n
        temp = dmax1(temp,dabs(w(k,2)))
305 continue
    go to 360
310 if (c(1) .ne. 2._dl) go to 320
    !              relative error control
    do 315 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/y(k)))
315 continue
    go to 360
320 if (c(1) .ne. 3._dl) go to 330
    !              weights are 1/max(c(2),abs(y(k)))
    do 325 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(2), dabs(y(k))) )
325 continue
    go to 360
330 if (c(1) .ne. 4._dl) go to 340
    !              weights are 1/max(c(k+30),abs(y(k)))
    do 335 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(k+30), dabs(y(k))) )
335 continue
    go to 360
340 if (c(1) .ne. 5._dl) go to 350
    !              weights are 1/c(k+30)
    do 345 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
345 continue
    go to 360
350 continue
    !              default case - weights are 1/max(1,abs(y(k)))
    do 355 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(1._dl, dabs(y(k))) )
355 continue
360 continue
    !
    !           calculate est - (the weighted max norm of w(*,2))*hmag*scale
    !              - est is intended to be a measure of the error  per  unit
    !              step in ytrial
    c(19) = temp*c(14)*c(15)
    !
    !        end stage 3
    !
    !        ***************************************************************
    !        * stage 4 - make decisions.                                   *
    !        ***************************************************************
    !
    !           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
    !
    !***********interrupt no 2 if requested
    if (c(9) .eq. 0._dl) go to 2222
    return
    !           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
    !
    if (ind .eq. 6) go to 410
    !              step accepted (ind .eq. 5), so update x, y from xtrial,
    !                 ytrial, add 1 to the no of successful steps, and set
    !                 the no of successive failures to zero
    x = c(17)
    do 400 k = 1, n
        y(k) = w(k,9)
400 continue
    c(22) = c(22) + 1._dl
    c(23) = 0._dl
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1._dl
    return
405 continue
    go to 420
410 continue
    !              step not accepted (ind .eq. 6), so add 1 to the no of
    !                 successive failures
    c(23) = c(23) + 1._dl
    !**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
    !
    !        end stage 4
    !
    go to 99999
    !     end loop
    !
    !  begin abort action
500 continue
    !

    write (*,*) 'Error in dverk, x =',x, 'xend=', xend
    call MpiStop()
    !
    !  end abort action
    !
    end subroutine dverk

    function Newton_Raphson(xxl,xxh,funcs, param, param2) result(xm)
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

    end function Newton_Raphson


    subroutine Gauss_Legendre(x,w,n)
    !Get Gauss-Legendre points x and weights w, for n points
    use constants
    implicit none
    integer, intent(in) :: n
    real(dl), intent(out) :: x(n), w(n)
    real(dl), parameter :: eps=1.d-15
    integer i, j, m
    real(dl) p1, p2, p3, pp, z, z1

    m=(n+1)/2
    !$OMP PARALLEL DO DEFAULT(PRIVATE), SHARED(x,w,n,m)
    do i=1,m
        z=cos(const_pi*(i-0.25_dl)/(n+0.5_dl))
        z1 = 0._dl
        do while (abs(z-z1) > eps)
            p1=1._dl
            p2=0._dl
            do j=1,n
                p3=p2
                p2=p1
                p1=((2*j-1)*z*p2-(j-1)*p3)/j
            end do
            pp=n*(z*p1-p2)/(z**2-1._dl)
            z1=z
            z=z1-p1/pp
        end do
        x(i)=-z
        x(n+1-i)=z
        w(i)=2/((1._dl-z**2)*pp**2)
        w(n+1-i)=w(i)
    end do
    !$OMP END PARALLEL DO

    end subroutine Gauss_Legendre

