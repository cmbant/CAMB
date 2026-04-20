    !Low-level numerical routines for splines and dverk for differential equation integration.

    module splines
    use Precision
    use Interpolation
    implicit none
    contains

    subroutine spline_def(x,y,n,d2)
    !Low-level initialize spline arrays with default boundary conditions
    integer, intent(in) :: n
    real(sp_acc), intent(in) :: x(n), y(n)
    real(sp_acc), intent(out) :: d2(n)

    call spline(x,y,n,SPLINE_DANGLE,SPLINE_DANGLE,d2)

    end subroutine spline_def

    subroutine splder(y,dy,n, g)
    !  Splder fits a cubic spline to y and returns the first derivatives at
    !  the grid points in dy.  Dy is equivalent to a 4th-order Pade
    !  difference formula for dy/di.
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

    end subroutine splder
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine splini(g,n)
    !  Splini must be called before splder to initialize array g in common.
    integer, intent(in) :: n
    real(dl), intent(out):: g(n)
    integer :: i

    g(1)=0._dl
    do i=2,n
        g(i)=1/(4._dl-g(i-1))
    end do
    end subroutine splini

    SUBROUTINE spline_deriv(x,y,y2,y1,n)
    !Get derivative y1 given array of x, y and y''
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
    !  Splint integrates a cubic spline, providing the output value
    !  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
    !  to y(i).
    !
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

    end module splines


    !This version is modified to pass an object parameter to the function on each call
    !Fortunately Fortran doesn't do type checking on functions, so we can pretend the
    !passed object parameter (EV) is any type we like. In reality it is just a pointer.
    subroutine dverk (EV,n, fcn, x, y, xend, tol_in, ind, c, nw, w)
    use Precision
    use MpiUtils
    use Config, only : GlobalError, error_evolution
    implicit none
    integer, intent(in) :: n, nw
    integer, intent(inout) :: ind
    real(dl), intent(inout) :: x, y(n), c(*), w(nw,9)
    real(dl), intent(in) :: xend, tol_in
    real(dl) :: tol, temp
    real(dl), parameter :: one_fifth = 1._dl / 5._dl
    real(dl), parameter :: dp_a21 = 1._dl / 5._dl
    real(dl), parameter :: dp_a31 = 3._dl / 40._dl
    real(dl), parameter :: dp_a32 = 9._dl / 40._dl
    real(dl), parameter :: dp_a41 = 44._dl / 45._dl
    real(dl), parameter :: dp_a42 = -56._dl / 15._dl
    real(dl), parameter :: dp_a43 = 32._dl / 9._dl
    real(dl), parameter :: dp_a51 = 19372._dl / 6561._dl
    real(dl), parameter :: dp_a52 = -25360._dl / 2187._dl
    real(dl), parameter :: dp_a53 = 64448._dl / 6561._dl
    real(dl), parameter :: dp_a54 = -212._dl / 729._dl
    real(dl), parameter :: dp_a61 = 9017._dl / 3168._dl
    real(dl), parameter :: dp_a62 = -355._dl / 33._dl
    real(dl), parameter :: dp_a63 = 46732._dl / 5247._dl
    real(dl), parameter :: dp_a64 = 49._dl / 176._dl
    real(dl), parameter :: dp_a65 = -5103._dl / 18656._dl
    real(dl), parameter :: dp_b1 = 35._dl / 384._dl
    real(dl), parameter :: dp_b3 = 500._dl / 1113._dl
    real(dl), parameter :: dp_b4 = 125._dl / 192._dl
    real(dl), parameter :: dp_b5 = -2187._dl / 6784._dl
    real(dl), parameter :: dp_b6 = 11._dl / 84._dl
    real(dl), parameter :: dp_e1 = 71._dl / 57600._dl
    real(dl), parameter :: dp_e3 = -71._dl / 16695._dl
    real(dl), parameter :: dp_e4 = 71._dl / 1920._dl
    real(dl), parameter :: dp_e5 = -17253._dl / 339200._dl
    real(dl), parameter :: dp_e6 = 22._dl / 525._dl
    real(dl), parameter :: dp_e7 = -1._dl / 40._dl
    real(dl), parameter :: machine_roundoff = epsilon(1._dl)
    real(dl), parameter :: machine_tiny = tiny(1._dl)
    logical :: resume_after_interrupt1, resume_after_interrupt2
    real :: EV !it isn't, but as long as it maintains it as a pointer we are OK
    !
    !***********************************************************************
    !                                                                      *
    ! This entry point preserves the historical DVERK interface, interrupt *
    ! handling, and work-array layout used by CAMB. The integration        *
    ! formula itself is the Dormand-Prince 5(4) embedded pair, and this    *
    ! copy has been refactored to modern Fortran syntax for readability    *
    ! while keeping the original control flow semantics. A PI-controller   *
    ! experiment is left commented in the accepted-step branch for         *
    ! reference; it used w(1,8) as spare workspace for the previous        *
    ! accepted error estimate, and showed mixed timings with no robust     *
    ! default speedup.                                                     *
    !                                                                      *
    ! The machine constants stored in c(10) and c(11) are initialized from *
    ! epsilon(1._dl) and tiny(1._dl), rather than from machine-specific    *
    ! literals.                                                            *
    !                                                                      *
    ! the c array is declared in this subroutine to have one element only, *
    ! although  more  elements  are  referenced  in this subroutine.  this *
    ! causes some compilers to issue warning messages.  there is,  though, *
    ! no  error  provided  c is declared sufficiently large in the calling *
    ! program, as described below.                                         *
    !
    !***********************************************************************
    !
    external fcn
    tol = tol_in
    !
    !***********************************************************************
    !                                                                      *
    !     purpose - this is a runge-kutta  subroutine  based  on Dormand-  *
    ! Prince's fifth and fourth order pair of formulas for finding approx- *
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
    !      historical interface after hull-enright-jackson   1/10/76       *
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
    !           uses the default value of hmax*(tol)**(1/5)                *
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
    !                                    c(10) machine roundoff estimate   *
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
    !  .        initialize machine constants, prev xend, flag, counts      *
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
    !  v     stage 2 - calc ytrial (adding 6 to no of fcn evals)        .  *
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

    resume_after_interrupt1 = .false.
    resume_after_interrupt2 = .false.

    select case (ind)
    case (1, 2)
        if (n > nw .or. tol <= 0._dl) call abort_dverk()

        if (ind == 1) then
            c(1:9) = 0._dl
        else
            c(1:9) = abs(c(1:9))
            if (c(1) == 4._dl .or. c(1) == 5._dl) then
                c(31:n+30) = abs(c(31:n+30))
            end if
        end if

        ! Initialize machine constants, previous xend, flag, and counts.
        c(10)    = machine_roundoff
        c(11)    = machine_tiny
        c(20)    = x
        c(21:24) = 0._dl

    case (3)
        ! Abort if xend was previously reached and x has changed,
        ! or xend has an invalid value.
        if (c(21) /= 0._dl .and. (x /= c(20) .or. xend == c(20))) then
            call abort_dverk()
        end if
        c(21) = 0._dl

    case (4)
        resume_after_interrupt1 = .true.

    case (5, 6)
        resume_after_interrupt2 = .true.

    case default
        call abort_dverk()
    end select

    step_loop: do

        !------------------------------------------------------------------
        ! Stage 1: prepare step size etc.
        ! Skip this part when resuming from interrupt 1 or 2.
        !------------------------------------------------------------------
        if (.not. resume_after_interrupt1 .and. .not. resume_after_interrupt2) then

            ! Error return (ind = -1) if number of function evaluations too great.
            if (c(7) /= 0._dl .and. c(24) >= c(7)) then
                ind = -1
                return
            end if

            ! Calculate slope (adding 1 to number of function evals) if ind /= 6.
            if (ind /= 6) then
                call fcn(EV, n, x, y, w(1, 1))
                c(24) = c(24) + 1._dl
            end if

            ! Calculate hmin - use default unless prescribed.
            c(13) = c(3)
            if (c(3) == 0._dl) then
                if (c(1) == 1._dl) then
                    ! Absolute error control.
                    c(12) = maxval(abs(y(1:n)))

                else if (c(1) == 2._dl) then
                    ! Relative error control.
                    c(12) = 1._dl

                else if (c(1) == 3._dl) then
                    ! Weights are 1 / max(c(2), abs(y(k))).
                    temp  = maxval(abs(y(1:n)) / c(2))
                    c(12) = min(temp, 1._dl)

                else if (c(1) == 4._dl) then
                    ! Weights are 1 / max(c(k+30), abs(y(k))).
                    temp  = maxval(abs(y(1:n)) / c(31:n+30))
                    c(12) = min(temp, 1._dl)

                else if (c(1) == 5._dl) then
                    ! Weights are 1 / c(k+30).
                    c(12) = maxval(abs(y(1:n)) / c(31:n+30))

                else
                    ! Default: weights are 1 / max(1, abs(y(k))).
                    temp  = maxval(abs(y(1:n)))
                    c(12) = min(temp, 1._dl)
                end if

                c(13) = 10._dl * max(c(11), c(10) * max(c(12) / tol, abs(x)))
            end if

            ! Calculate scale - use default unless prescribed.
            c(15) = c(5)
            if (c(5) == 0._dl) c(15) = 1._dl

            ! Calculate hmax.
            if (c(6) /= 0._dl .and. c(5) /= 0._dl) then
                c(16) = min(c(6), 2._dl / c(5))
            else if (c(6) /= 0._dl) then
                c(16) = c(6)
            else if (c(5) /= 0._dl) then
                c(16) = 2._dl / c(5)
            else
                c(16) = 2._dl
            end if

            ! Error return (ind = -2) if hmin > hmax.
            if (c(13) > c(16)) then
                ind = -2
                return
            end if

            ! Calculate preliminary hmag.
            if (ind <= 2) then
                ! Initial entry.
                c(14) = c(4)
                if (c(4) == 0._dl) c(14) = c(16) * tol**one_fifth

            else if (c(23) <= 1._dl) then
                ! Successful step, or at most one failure: original
                ! proportional controller.
                ! PI-controller experiment left for reference only:
                ! mixed timings, no robust default speedup on tested CAMB workloads.
                ! If re-enabled, w(1,8) must be initialized on entry and
                ! updated after each accepted step with max(c(19), c(11)).
                ! if (c(19) <= c(11)) then
                !     temp = 2._dl
                ! else if (c(22) <= 1._dl .or. w(1, 8) <= c(11)) then
                !     temp = 0.9_dl * (tol / c(19))**0.20_dl
                ! else
                !     temp = 0.9_dl * (tol / c(19))**0.20_dl * (w(1, 8) / c(19))**0.04_dl
                ! end if
                ! temp = min(2._dl, max(0.5_dl, temp))
                ! c(14) = temp * c(14)
                temp = 2._dl * c(14)
                if (tol < (2._dl / 0.9_dl)**5 * c(19)) then
                    temp = 0.9_dl * (tol / c(19))**one_fifth * c(14)
                end if
                c(14) = max(temp, 0.5_dl * c(14))

            else
                ! Two or more successive failures.
                c(14) = 0.5_dl * c(14)
            end if

            c(14) = min(c(14), c(16))
            c(14) = max(c(14), c(13))

            ! Interrupt no. 1.
            if (c(8) /= 0._dl) then
                ind = 4
                return
            end if
        end if

        !------------------------------------------------------------------
        ! Resume point for ind = 4.
        ! Skip this block when resuming from interrupt 2.
        !------------------------------------------------------------------
        if (.not. resume_after_interrupt2) then

            ! Calculate hmag, xtrial.
            if (c(14) < abs(xend - x)) then
                c(14) = min(c(14), 0.5_dl * abs(xend - x))
                c(17) = x + sign(c(14), xend - x)
            else
                c(14) = abs(xend - x)
                c(17) = xend
            end if

            c(18) = c(17) - x

            !------------------------------------------------------------------
            ! Stage 2: Dormand-Prince 5(4) trial step.
            !------------------------------------------------------------------
            w(1:n, 9) = y(1:n) + c(18) * dp_a21 * w(1:n, 1)
            call fcn(EV, n, x + c(18) / 5._dl, w(1, 9), w(1, 2))

            w(1:n, 9) = y(1:n) + c(18) * (dp_a31 * w(1:n, 1) + dp_a32 * w(1:n, 2))
            call fcn(EV, n, x + c(18) * (3._dl / 10._dl), w(1, 9), w(1, 3))

            w(1:n, 9) = y(1:n) + c(18) * (dp_a41 * w(1:n, 1) + dp_a42 * w(1:n, 2) + &
                dp_a43 * w(1:n, 3))
            call fcn(EV, n, x + c(18) * (4._dl / 5._dl), w(1, 9), w(1, 4))

            w(1:n, 9) = y(1:n) + c(18) * (dp_a51 * w(1:n, 1) + dp_a52 * w(1:n, 2) + &
                dp_a53 * w(1:n, 3) + dp_a54 * w(1:n, 4))
            call fcn(EV, n, x + c(18) * (8._dl / 9._dl), w(1, 9), w(1, 5))

            w(1:n, 9) = y(1:n) + c(18) * (dp_a61 * w(1:n, 1) + dp_a62 * w(1:n, 2) + &
                dp_a63 * w(1:n, 3) + dp_a64 * w(1:n, 4) + &
                dp_a65 * w(1:n, 5))
            call fcn(EV, n, x + c(18), w(1, 9), w(1, 6))

            w(1:n, 9) = y(1:n) + c(18) * (dp_b1 * w(1:n, 1) + dp_b3 * w(1:n, 3) + &
                dp_b4 * w(1:n, 4) + dp_b5 * w(1:n, 5) + &
                dp_b6 * w(1:n, 6))
            call fcn(EV, n, x + c(18), w(1, 9), w(1, 7))

            c(24) = c(24) + 6._dl

            !------------------------------------------------------------------
            ! Stage 3: error estimate.
            !------------------------------------------------------------------
            w(1:n, 2) = dp_e1 * w(1:n, 1) + dp_e3 * w(1:n, 3) + dp_e4 * w(1:n, 4) + &
                dp_e5 * w(1:n, 5) + dp_e6 * w(1:n, 6) + dp_e7 * w(1:n, 7)

            if (c(1) == 1._dl) then
                temp = maxval(abs(w(1:n, 2)))

            else if (c(1) == 2._dl) then
                temp = maxval(abs(w(1:n, 2) / y(1:n)))

            else if (c(1) == 3._dl) then
                temp = maxval(abs(w(1:n, 2)) / max(c(2), abs(y(1:n))))

            else if (c(1) == 4._dl) then
                temp = maxval(abs(w(1:n, 2)) / max(c(31:n+30), abs(y(1:n))))

            else if (c(1) == 5._dl) then
                temp = maxval(abs(w(1:n, 2) / c(31:n+30)))

            else
                temp = maxval(abs(w(1:n, 2)) / max(1._dl, abs(y(1:n))))
            end if

            c(19) = temp * c(14) * c(15)

            !------------------------------------------------------------------
            ! Stage 4: accept/reject decision.
            !------------------------------------------------------------------
            ind = 5
            if (c(19) > tol) ind = 6

            ! Interrupt no. 2.
            if (c(9) /= 0._dl) return
        end if

        !------------------------------------------------------------------
        ! Resume point for ind = 5 or 6.
        !------------------------------------------------------------------
        if (ind == 5) then
            ! Step accepted.
            x        = c(17)
            y(1:n)   = w(1:n, 9)
            c(22)    = c(22) + 1._dl
            c(23)    = 0._dl

            if (x == xend) then
                ind   = 3
                c(20) = xend
                c(21) = 1._dl
                return
            end if

            w(1:n, 1) = w(1:n, 7)
            ind = 6

        else
            ! Step rejected.
            c(23) = c(23) + 1._dl

            if (c(14) <= c(13)) then
                ind = -3
                return
            end if
        end if

        resume_after_interrupt1 = .false.
        resume_after_interrupt2 = .false.

    end do step_loop

    contains

    subroutine abort_dverk()
    write (*,*) 'Error in dverk, x =', x, ' xend =', xend
    call GlobalError('DVERK error', error_evolution)
    end subroutine abort_dverk

    end subroutine dverk
