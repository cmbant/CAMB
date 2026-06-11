program bjl_recurrence_validation
! Validates the fast BJL dispatcher (fortran/bessels.f90) against
! BJL_RECURRENCE (Miller/forward recurrence reference) for L=2..L_MAX.
!
! For each L this samples:
!   - a dense window around the turning point, x_turn-10*L^(1/3) to
!     x_turn+10*L^(1/3) (240 points for L<100, 60 points otherwise),
!   - a log-spaced tail out to X_FAR (~30000),
!   - pairs of points straddling each branch-selection gate in BJL.
!
! Errors are normalized by the peak |j_l| estimated from BJL_RECURRENCE near
! the turning point. The worst normalized error is reported and the program
! fails (error stop) if it exceeds TARGET_ERR.
!
! Build and run, e.g.:
!   gfortran -O3 -ffree-line-length-none -Ifortran/Releaselib \
!     fortran/tests/bjl_recurrence_validation.f90 fortran/Releaselib/camblib.so \
!     -Wl,-rpath,$PWD/fortran/Releaselib -o fortran/tests/bjl_recurrence_validation
!   fortran/tests/bjl_recurrence_validation
! Remove the generated executable afterwards.

use Precision
use FlatBessels, only: bjl, BJL_RECURRENCE
implicit none

integer, parameter :: L_MIN = 2
integer, parameter :: L_MAX = 5000
integer, parameter :: L_DENSE_MAX = 99
real(dl), parameter :: TARGET_ERR = 1.0e-5_dl
real(dl), parameter :: X_FAR = 30000._dl
real(dl), parameter :: GATE_REL_EPS = 1.0e-6_dl

! Gate constants mirrored from the BJL dispatcher in fortran/bessels.f90;
! keep these in sync if the gates there change.
integer,  parameter :: BJL_RECURRENCE_MAX_L  = 25
real(dl), parameter :: BJL_RECURRENCE_ETA_LOW  = -5.0_dl
real(dl), parameter :: BJL_RECURRENCE_ETA_HIGH =  5.6_dl
real(dl), parameter :: BJL_AIRY_ETA_LOW  = -2.4_dl
real(dl), parameter :: BJL_PEAK_ETA_LOW  = -0.65_dl
real(dl), parameter :: BJL_PEAK_ETA_HIGH =  0.65_dl
real(dl), parameter :: BJL_AIRY_ETA_HIGH =  3.85_dl
real(dl), parameter :: BJL_AIRY_U_LOW  = 0.26_dl
real(dl), parameter :: BJL_AIRY_U_HIGH = 0.42_dl

integer :: l, worst_l
integer :: n_points, n_bjl
real(dl) :: worst_err, worst_x
real(dl) :: t_start, t_end, sum_jl

worst_err = -1._dl
worst_l = -1
worst_x = 0._dl
n_points = 0

call cpu_time(t_start)

do l = L_MIN, L_MAX
    call check_l(l, worst_err, worst_x, worst_l, n_points)
end do

call cpu_time(t_end)

write(*,'(a,i0,a,i0)') 'L range checked: ', L_MIN, '..', L_MAX
write(*,'(a,i0)') 'total points checked: ', n_points
write(*,'(a,f10.3)') 'cpu seconds (bjl + bjl_recurrence): ', t_end - t_start
write(*,'(a,es12.4,a,i0,a,f14.4)') 'worst peak-normalized error: ', worst_err, &
    ' at L=', worst_l, ', x=', worst_x
write(*,'(a,es10.3)') 'target: ', TARGET_ERR

! Separate bjl-only timing pass over the same near+far x sampling, with no
! reference-recurrence calls.
sum_jl = 0._dl
n_bjl = 0
call cpu_time(t_start)
do l = L_MIN, L_MAX
    call time_l(l, sum_jl, n_bjl)
end do
call cpu_time(t_end)
write(*,'(a,es12.4,a,es12.4,a)') 'cpu seconds per bjl eval: ', &
    (t_end - t_start)/real(n_bjl, dl), ' (sum=', sum_jl, ', anti-DCE check)'

if (worst_err > TARGET_ERR) then
    write(*,'(a)') 'FAIL: worst error exceeds target'
    error stop 1
else
    write(*,'(a)') 'PASS'
end if

contains

subroutine check_l(l, worst_err, worst_x, worst_l, n_points)
integer, intent(in) :: l
real(dl), intent(inout) :: worst_err, worst_x
integer, intent(inout) :: worst_l, n_points

real(dl) :: nu, l3, nu23, x_turn, x_lo, x_hi_near, peak_abs
real(dl) :: x, frac
integer :: n_near, n_far, i

nu = real(l, dl) + 0.5_dl
l3 = nu**(1.0_dl/3.0_dl)
nu23 = l3*l3
x_turn = nu
x_lo = max(1.0e-3_dl, x_turn - 10.0_dl*l3)
x_hi_near = x_turn + 10.0_dl*l3

if (l <= L_DENSE_MAX) then
    n_near = 240
    n_far = 40
else
    n_near = 60
    n_far = 30
end if

! First positive peak amplitude of j_l(x), Apeak ~= j_l(x_peak) where
! d/dx j_l(x_peak) = 0. Max relative error <~ 0.7% for l = 1..5000 (better
! for large l), good enough for error normalization. (l = 0 is not reached
! here since L_MIN >= 2.)
peak_abs = 0.845843_dl*nu**(-5.0_dl/6.0_dl) &
    *(1.0_dl - 0.55424_dl/nu23 + 0.25865_dl/nu23**2)

! Dense window around the turning point.
do i = 0, n_near-1
    x = x_lo + (x_hi_near - x_lo)*real(i, dl)/real(n_near-1, dl)
    call check_point(l, x, peak_abs, worst_err, worst_x, worst_l, n_points)
end do

! Log-spaced tail out to X_FAR.
if (X_FAR > x_hi_near) then
    do i = 1, n_far
        frac = real(i, dl)/real(n_far, dl)
        x = x_hi_near*(X_FAR/x_hi_near)**frac
        call check_point(l, x, peak_abs, worst_err, worst_x, worst_l, n_points)
    end do
end if

call check_gates(l, nu, l3, nu23, peak_abs, worst_err, worst_x, worst_l, n_points)

end subroutine check_l


subroutine time_l(l, sum_jl, n_bjl)
! Same near+far x sampling as check_l, but timing bjl alone (no reference
! recurrence calls and no gate points).
integer, intent(in) :: l
real(dl), intent(inout) :: sum_jl
integer, intent(inout) :: n_bjl

real(dl) :: nu, l3, x_turn, x_lo, x_hi_near
real(dl) :: x, frac, jl
integer :: n_near, n_far, i

nu = real(l, dl) + 0.5_dl
l3 = nu**(1.0_dl/3.0_dl)
x_turn = nu
x_lo = max(1.0e-3_dl, x_turn - 10.0_dl*l3)
x_hi_near = x_turn + 10.0_dl*l3

if (l <= L_DENSE_MAX) then
    n_near = 240
    n_far = 40
else
    n_near = 60
    n_far = 30
end if

do i = 0, n_near-1
    x = x_lo + (x_hi_near - x_lo)*real(i, dl)/real(n_near-1, dl)
    call bjl(l, x, jl)
    sum_jl = sum_jl + jl
    n_bjl = n_bjl + 1
end do

if (X_FAR > x_hi_near) then
    do i = 1, n_far
        frac = real(i, dl)/real(n_far, dl)
        x = x_hi_near*(X_FAR/x_hi_near)**frac
        call bjl(l, x, jl)
        sum_jl = sum_jl + jl
        n_bjl = n_bjl + 1
    end do
end if

end subroutine time_l


subroutine check_point(l, x, peak_abs, worst_err, worst_x, worst_l, n_points)
integer, intent(in) :: l
real(dl), intent(in) :: x, peak_abs
real(dl), intent(inout) :: worst_err, worst_x
integer, intent(inout) :: worst_l, n_points
real(dl) :: jl_fast, jl_ref, err

if (x <= 0._dl) return

call bjl(l, x, jl_fast)
call BJL_RECURRENCE(l, x, jl_ref)

err = abs(jl_fast - jl_ref)/peak_abs
n_points = n_points + 1

if (err > worst_err) then
    worst_err = err
    worst_x = x
    worst_l = l
end if
end subroutine check_point


subroutine check_gates(l, nu, l3, nu23, peak_abs, worst_err, worst_x, worst_l, n_points)
integer, intent(in) :: l
real(dl), intent(in) :: nu, l3, nu23, peak_abs
real(dl), intent(inout) :: worst_err, worst_x
integer, intent(inout) :: worst_l, n_points
real(dl) :: eta_airy_low, eta_airy_high, g

if (l < 7) then
    ! Explicit low-L (2..6) small-x series/closed-form switch in BJL.
    select case (l)
    case (2)
        g = 3.0e-1_dl
    case (3)
        g = 4.0e-1_dl
    case (4)
        g = 6.0e-1_dl
    case default
        g = 1.0_dl
    end select
    call add_gate(l, g, peak_abs, worst_err, worst_x, worst_l, n_points)
    return
end if

! Deep pre-peak boundary: AX^2/L = 0.5
call add_gate(l, sqrt(0.5_dl*real(l, dl)), peak_abs, worst_err, worst_x, worst_l, n_points)

! Far-x oscillatory boundary: L^2/AX = 1.2
call add_gate(l, real(l, dl)**2/1.2_dl, peak_abs, worst_err, worst_x, worst_l, n_points)

! Postpeak-Debye fast-precheck boundaries.
call add_gate(l, 1.5_dl*nu, peak_abs, worst_err, worst_x, worst_l, n_points)
call add_gate(l, 2.5_dl*nu, peak_abs, worst_err, worst_x, worst_l, n_points)
if (l > BJL_RECURRENCE_MAX_L) then
    ! Internal simplified-Debye threshold for the far post-peak tail.
    call add_gate(l, 3.0_dl*nu, peak_abs, worst_err, worst_x, worst_l, n_points)
end if

! Recurrence-band eta gates (only reachable for l <= BJL_RECURRENCE_MAX_L).
call add_gate(l, nu + BJL_RECURRENCE_ETA_LOW*l3, peak_abs, worst_err, worst_x, worst_l, n_points)
call add_gate(l, nu + BJL_RECURRENCE_ETA_HIGH*l3, peak_abs, worst_err, worst_x, worst_l, n_points)

! Airy/peak eta gates, including the u-cap that binds for low l.
eta_airy_low = max(BJL_AIRY_ETA_LOW, -BJL_AIRY_U_LOW*nu23)
eta_airy_high = min(BJL_AIRY_ETA_HIGH, BJL_AIRY_U_HIGH*nu23)

call add_gate(l, nu + eta_airy_low*l3, peak_abs, worst_err, worst_x, worst_l, n_points)
call add_gate(l, nu + BJL_PEAK_ETA_LOW*l3, peak_abs, worst_err, worst_x, worst_l, n_points)
call add_gate(l, nu + BJL_PEAK_ETA_HIGH*l3, peak_abs, worst_err, worst_x, worst_l, n_points)
call add_gate(l, nu + eta_airy_high*l3, peak_abs, worst_err, worst_x, worst_l, n_points)

! BJL_RECURRENCE_MAX_L itself (recurrence band L cutoff): check the
! neighbouring L value too, so both sides of the L cutoff are covered.
if (l == BJL_RECURRENCE_MAX_L .or. l == BJL_RECURRENCE_MAX_L + 1) then
    call check_point(l, nu, peak_abs, worst_err, worst_x, worst_l, n_points)
end if

end subroutine check_gates


subroutine add_gate(l, x_gate, peak_abs, worst_err, worst_x, worst_l, n_points)
integer, intent(in) :: l
real(dl), intent(in) :: x_gate, peak_abs
real(dl), intent(inout) :: worst_err, worst_x
integer, intent(inout) :: worst_l, n_points
real(dl) :: dx

if (x_gate <= 0._dl .or. x_gate > X_FAR*1.5_dl) return

dx = max(1.0e-9_dl, x_gate*GATE_REL_EPS)
call check_point(l, x_gate - dx, peak_abs, worst_err, worst_x, worst_l, n_points)
call check_point(l, x_gate + dx, peak_abs, worst_err, worst_x, worst_l, n_points)
end subroutine add_gate

end program bjl_recurrence_validation
