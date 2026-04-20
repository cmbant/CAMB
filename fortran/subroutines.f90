    !Low-level numerical routines for splines.

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
