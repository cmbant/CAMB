    program phi_olver_gate_validation
    use Precision
    use HypersphericalBesselOlver, only: phi_olver
    use SpherBessels, only: phi_recurs, phi_first_peak_amplitude
    use omp_lib
    implicit none

    integer, parameter :: TOP_N = 20
    integer, parameter :: MAX_RAW_POINTS = 12000
    integer, parameter :: NBINS = 10
    real(dl), parameter :: PI = 3.1415926535897932384626433832795_dl
    real(dl), parameter :: WARN_ERR = 1.0e-4_dl
    real(dl), parameter :: TARGET_ERR = 2.0e-4_dl
    real(dl), parameter :: SPLINE_TARGET_ERR = 1.0e-4_dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_FLOOR = 0.095_dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_JOIN = 0.12_dl
    real(dl), parameter :: OLVER_OPEN_L_JOIN = 500._dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_LOW_EXP = 0.70_dl
    real(dl), parameter :: OLVER_OPEN_ALPHA_HIGH_EXP = 0.14_dl
    real(dl), parameter :: OLVER_GATE_OPEN_EPS = 3.0e-3_dl
    real(dl), parameter :: OLVER_GATE_CLOSED_EPS = 7.0e-3_dl
    real(dl), parameter :: SMALLCHI_GATE_ALPHA = 2.5_dl
    integer, parameter :: SMALLCHI_GATE_LOW_ALPHA_L_MIN = 50
    real(dl), parameter :: SMALLCHI_GATE_METRIC = 5.0e-2_dl

    type worst_record
        character(len=8) :: geometry = ""
        integer :: l = 0
        integer :: k = 0
        integer :: index = 0
        real(dl) :: nu = 0._dl
        real(dl) :: alpha = 0._dl
        real(dl) :: chi = 0._dl
        real(dl) :: metric = 0._dl
        real(dl) :: err = -1._dl
        real(dl) :: expected = TARGET_ERR
        real(dl) :: peak = 0._dl
    end type worst_record

    logical :: quick
    integer :: total_cases, approx_warnings, approx_failures, spline_failures
    integer :: total_points, nthreads
    integer :: total_fast_points, total_fallback_points
    integer :: bin_total(2, NBINS), bin_fast(2, NBINS), bin_fallback(2, NBINS)
    real(dl) :: approx_cpu, recurs_cpu, spline_cpu, wall_start, wall_end
    type(worst_record) :: worst_approx(TOP_N), worst_spline(TOP_N)

    quick = command_argument_is("quick")
    nthreads = omp_get_max_threads()
    total_cases = 0
    total_points = 0
    total_fast_points = 0
    total_fallback_points = 0
    bin_total = 0
    bin_fast = 0
    bin_fallback = 0
    approx_failures = 0
    approx_warnings = 0
    spline_failures = 0
    approx_cpu = 0._dl
    recurs_cpu = 0._dl
    spline_cpu = 0._dl
    wall_start = omp_get_wtime()

    write(*, '(a)') "phi_olver gate validation against phi_recurs"
    write(*, '(a, l1)') "quick mode: ", quick
    write(*, '(a, i0)') "OpenMP max threads: ", nthreads
    write(*, '(a, es10.3)') "warning peak-normalized error: ", WARN_ERR
    write(*, '(a, es10.3)') "failure peak-normalized error: ", TARGET_ERR

    call run_open_cases()
    call run_closed_cases()

    wall_end = omp_get_wtime()

    write(*, '(/, a)') "Summary"
    write(*, '(a, i0)') "cases: ", total_cases
    write(*, '(a, i0)') "points: ", total_points
    write(*, '(a, i0)') "predicted fast-path points: ", total_fast_points
    write(*, '(a, i0)') "predicted recursive-fallback points: ", total_fallback_points
    write(*, '(a, i0)') "approx warnings above 1e-4: ", approx_warnings
    write(*, '(a, i0)') "approx failures above target: ", approx_failures
    write(*, '(a, i0)') "recursive spline warnings above target: ", spline_failures
    write(*, '(a, f12.4)') "phi_olver CPU seconds: ", approx_cpu
    write(*, '(a, f12.4)') "phi_recurs CPU seconds: ", recurs_cpu
    write(*, '(a, f12.4)') "spline-check CPU seconds: ", spline_cpu
    write(*, '(a, f12.4)') "wall seconds: ", wall_end - wall_start
    if (total_points > 0) then
        write(*, '(a, es12.4)') "phi_olver CPU seconds/eval: ", approx_cpu / real(total_points, dl)
        write(*, '(a, es12.4)') "phi_recurs CPU seconds/eval: ", recurs_cpu / real(total_points, dl)
    end if

    call print_fast_fraction_table()
    call print_records("Worst phi_olver peak-normalized errors", worst_approx)
    call print_records("Worst recursive 1-in-5 spline peak-normalized errors", worst_spline)
    call run_timing_benchmarks()

    if (approx_failures > 0 .or. spline_failures > 0) then
        error stop 1
    end if

    contains

    logical function command_argument_is(expected)
    character(len=*), intent(in) :: expected
    character(len=32) :: arg
    integer :: narg

    narg = command_argument_count()
    if (narg < 1) then
        command_argument_is = .false.
    else
        call get_command_argument(1, arg)
        command_argument_is = trim(arg) == expected
    end if
    end function command_argument_is


    subroutine run_open_cases()
    integer, allocatable :: ls(:)
    real(dl), allocatable :: alphas(:)
    integer :: il, ia, l
    real(dl) :: alpha, nu

    call make_l_grid(ls)
    call make_open_alpha_grid(alphas)

    do il = 1, size(ls)
        l = ls(il)
        do ia = 1, size(alphas)
            alpha = alphas(ia)
            nu = alpha * real(l, dl)
            if (nu <= 0._dl) cycle
            call run_case(l, -1, nu)
        end do
    end do
    end subroutine run_open_cases


    subroutine run_closed_cases()
    integer, allocatable :: ls(:)
    integer, allocatable :: deltas(:)
    integer :: il, id, l, delta
    real(dl) :: nu

    call make_l_grid(ls)
    call make_closed_delta_grid(deltas)

    do il = 1, size(ls)
        l = ls(il)
        do id = 1, size(deltas)
            delta = deltas(id)
            nu = real(l + delta, dl)
            if (nint(nu) < 3) cycle
            call run_case(l, 1, nu)
        end do
    end do
    end subroutine run_closed_cases


    subroutine make_l_grid(ls)
    integer, allocatable, intent(out) :: ls(:)
    integer :: n, l

    if (quick) then
        allocate(ls(18))
        ls = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
        return
    end if

    n = 0
    n = n + 250
    n = n + ((600 - 260) / 20 + 1)
    n = n + ((2000 - 700) / 100 + 1)
    n = n + 7
    allocate(ls(n))

    n = 0
    do l = 1, 250
        n = n + 1
        ls(n) = l
    end do
    do l = 260, 600, 20
        n = n + 1
        ls(n) = l
    end do
    do l = 700, 2000, 100
        n = n + 1
        ls(n) = l
    end do
    ls(n + 1:n + 7) = [2500, 3000, 4000, 5000, 6000, 8000, 10000]
    end subroutine make_l_grid


    subroutine make_open_alpha_grid(alphas)
    real(dl), allocatable, intent(out) :: alphas(:)

    if (quick) then
        allocate(alphas(8))
        alphas = [0.002_dl, 0.02_dl, 0.05_dl, 0.10_dl, 0.12_dl, 0.20_dl, 1.0_dl, 2.9_dl]
    else
        allocate(alphas(17))
        alphas = [0.002_dl, 0.005_dl, 0.01_dl, 0.02_dl, 0.05_dl, 0.08_dl, 0.10_dl, 0.12_dl, &
            0.15_dl, 0.20_dl, 0.30_dl, 0.50_dl, 0.75_dl, 1.0_dl, 1.5_dl, 2.5_dl, 2.9_dl]
    end if
    end subroutine make_open_alpha_grid


    subroutine make_closed_delta_grid(deltas)
    integer, allocatable, intent(out) :: deltas(:)

    if (quick) then
        allocate(deltas(9))
        deltas = [1, 2, 5, 20, 100, 175, 200, 500, 1500]
    else
        allocate(deltas(20))
        deltas = [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 80, 100, 150, 175, 200, 300, 500, &
            900, 1500, 3000]
    end if
    end subroutine make_closed_delta_grid


    subroutine run_case(l, k, nu)
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu
    real(dl), allocatable :: chi(:), approx(:), recurs(:)
    integer :: n, i, imax, ispline
    integer :: nfast, nfallback
    real(dl) :: t0, t1, peak, err, spline_err, spline_expected, metric
    type(worst_record) :: rec

    if (k == -1) then
        call make_open_chi_grid(l, nu, chi)
    else
        call make_closed_chi_grid(chi)
    end if
    n = size(chi)
    allocate(approx(n), recurs(n))

    call cpu_time(t0)
!$omp parallel do default(shared) private(i) schedule(static)
    do i = 1, n
        approx(i) = phi_olver(l, k, nu, chi(i))
    end do
!$omp end parallel do
    call cpu_time(t1)
    approx_cpu = approx_cpu + t1 - t0

    call cpu_time(t0)
!$omp parallel do default(shared) private(i) schedule(static)
    do i = 1, n
        recurs(i) = phi_recurs(l, k, nu, chi(i))
    end do
!$omp end parallel do
    call cpu_time(t1)
    recurs_cpu = recurs_cpu + t1 - t0

    peak = phi_first_peak_amplitude(l, k, nu)
    if (peak <= tiny(peak)) peak = tiny(peak)
    imax = maxloc(abs(approx - recurs), dim=1)
    err = abs(approx(imax) - recurs(imax)) / peak
    metric = case_metric(l, k, nu, chi(imax))

    rec = make_record(l, k, nu, chi(imax), metric, err, TARGET_ERR, peak, imax)
    call insert_record(worst_approx, rec)
    if (err > WARN_ERR) approx_warnings = approx_warnings + 1
    if (err > TARGET_ERR) approx_failures = approx_failures + 1

! l <= 2 uses exact seeds/fallback rather than the recurrence machinery this
! spline sentinel is meant to catch; the coarse every-fifth spline can be a
! poor surrogate for those low-nu seed shapes.
    if (l > 2) then
        call cpu_time(t0)
        call recursive_spline_error(l, k, nu, chi, recurs, peak, spline_err, spline_expected, ispline)
        call cpu_time(t1)
        spline_cpu = spline_cpu + t1 - t0
        metric = case_metric(l, k, nu, chi(ispline))
        rec = make_record(l, k, nu, chi(ispline), metric, spline_err, spline_expected, peak, ispline)
        call insert_record(worst_spline, rec)
        if (spline_err > spline_expected) spline_failures = spline_failures + 1
    end if

    total_cases = total_cases + 1
    total_points = total_points + n
    nfast = 0
    nfallback = 0
    do i = 1, n
        if (predicts_recursive_fallback(l, k, nu, chi(i))) then
            nfallback = nfallback + 1
        else
            nfast = nfast + 1
        end if
    end do
    total_fast_points = total_fast_points + nfast
    total_fallback_points = total_fallback_points + nfallback
    call accumulate_fast_fraction_bins(l, k, n, nfast, nfallback)
    if (mod(total_cases, 100) == 0) then
        write(*, '(a, i0, a, i0, a, es10.3, a, es10.3)') "completed cases: ", total_cases, &
            " points: ", total_points, " worst approx: ", worst_approx(1)%err, &
            " worst spline: ", worst_spline(1)%err
    end if

    deallocate(chi, approx, recurs)
    end subroutine run_case


    subroutine accumulate_fast_fraction_bins(l, k, n, nfast, nfallback)
    integer, intent(in) :: l, k, n, nfast, nfallback
    integer :: ig, ib

    if (k == -1) then
        ig = 1
    else if (k == 1) then
        ig = 2
    else
        return
    end if
    ib = l_bin(l)
    bin_total(ig, ib) = bin_total(ig, ib) + n
    bin_fast(ig, ib) = bin_fast(ig, ib) + nfast
    bin_fallback(ig, ib) = bin_fallback(ig, ib) + nfallback
    end subroutine accumulate_fast_fraction_bins


    integer function l_bin(l)
    integer, intent(in) :: l

    select case (l)
    case (1:10)
        l_bin = 1
    case (11:50)
        l_bin = 2
    case (51:100)
        l_bin = 3
    case (101:250)
        l_bin = 4
    case (251:500)
        l_bin = 5
    case (501:1000)
        l_bin = 6
    case (1001:2000)
        l_bin = 7
    case (2001:4000)
        l_bin = 8
    case (4001:6000)
        l_bin = 9
    case default
        l_bin = 10
    end select
    end function l_bin


    function l_bin_label(ib) result(label)
    integer, intent(in) :: ib
    character(len=12) :: label

    select case (ib)
    case (1)
        label = "1-10"
    case (2)
        label = "11-50"
    case (3)
        label = "51-100"
    case (4)
        label = "101-250"
    case (5)
        label = "251-500"
    case (6)
        label = "501-1000"
    case (7)
        label = "1001-2000"
    case (8)
        label = "2001-4000"
    case (9)
        label = "4001-6000"
    case default
        label = "6001-10000"
    end select
    end function l_bin_label


    subroutine print_fast_fraction_table()
    integer :: ig, ib
    real(dl) :: frac
    character(len=8) :: geometry

    write(*, '(/, a)') "Predicted phi_olver fast-path fraction by geometry and L bin"
    write(*, '(a)') "geometry Lbin total fast fallback fast_fraction"
    do ig = 1, 2
        geometry = merge("open    ", "closed  ", ig == 1)
        do ib = 1, NBINS
            if (bin_total(ig, ib) == 0) cycle
            frac = real(bin_fast(ig, ib), dl) / real(bin_total(ig, ib), dl)
            write(*, '(a8, 1x, a12, 1x, i10, 1x, i10, 1x, i10, 1x, f9.5)') &
                geometry, l_bin_label(ib), bin_total(ig, ib), bin_fast(ig, ib), bin_fallback(ig, ib), frac
        end do
    end do
    end subroutine print_fast_fraction_table


    logical function predicts_recursive_fallback(l, k, nu, chi)
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu, chi
    real(dl) :: alpha_gate, metric

    predicts_recursive_fallback = .false.
    if (k == 0) return
    if (l <= 2) then
        predicts_recursive_fallback = .true.
        return
    end if
    if (chi <= 1.0e-12_dl) return

    alpha_gate = nu / real(l, dl)
    if (use_smallchi_map_test(l, nu, chi, alpha_gate)) return

    if (k == 1) then
        metric = chi / (2._dl * (nu - real(l, dl)))
        predicts_recursive_fallback = metric > OLVER_GATE_CLOSED_EPS
    else if (k == -1) then
        if (alpha_gate < open_alpha_cut(l)) then
            metric = chi / (2._dl * max(nu, tiny(1._dl)))
            predicts_recursive_fallback = metric > OLVER_GATE_OPEN_EPS
        end if
    end if
    end function predicts_recursive_fallback


    real(dl) function open_alpha_cut(l)
    integer, intent(in) :: l
    real(dl) :: ell

    ell = real(l, dl)
    if (ell < OLVER_OPEN_L_JOIN) then
        open_alpha_cut = OLVER_OPEN_ALPHA_JOIN * (OLVER_OPEN_L_JOIN / ell)**OLVER_OPEN_ALPHA_LOW_EXP
    else
        open_alpha_cut = OLVER_OPEN_ALPHA_JOIN * (OLVER_OPEN_L_JOIN / ell)**OLVER_OPEN_ALPHA_HIGH_EXP
    end if
    open_alpha_cut = max(OLVER_OPEN_ALPHA_FLOOR, open_alpha_cut)
    end function open_alpha_cut


    logical function use_smallchi_map_test(l, nu, chi, alpha_gate)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, chi, alpha_gate

    use_smallchi_map_test = l > 0 .and. &
        (alpha_gate > SMALLCHI_GATE_ALPHA .or. &
        (l >= SMALLCHI_GATE_LOW_ALPHA_L_MIN .and. alpha_gate > 1._dl)) .and. &
        real(l, dl)**2 * chi**7 / nu < SMALLCHI_GATE_METRIC
    end function use_smallchi_map_test


    function make_record(l, k, nu, chi, metric, err, expected, peak, index) result(rec)
    integer, intent(in) :: l, k, index
    real(dl), intent(in) :: nu, chi, metric, err, expected, peak
    type(worst_record) :: rec

    if (k == -1) then
        rec%geometry = "open"
    else
        rec%geometry = "closed"
    end if
    rec%l = l
    rec%k = k
    rec%nu = nu
    rec%alpha = nu / real(l, dl)
    rec%chi = chi
    rec%metric = metric
    rec%err = err
    rec%expected = expected
    rec%peak = peak
    rec%index = index
    end function make_record


    real(dl) function case_metric(l, k, nu, chi)
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu, chi

    if (k == 1) then
        case_metric = chi / (2._dl * (nu - real(l, dl)))
    else
        case_metric = chi / (2._dl * max(nu, tiny(1._dl)))
    end if
    end function case_metric


    subroutine make_open_chi_grid(l, nu, chi)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu
    real(dl), allocatable, intent(out) :: chi(:)
    real(dl), allocatable :: raw(:)
    integer :: n
    real(dl) :: ell, turn, end_chi, width

    allocate(raw(MAX_RAW_POINTS))
    n = 0
    ell = sqrt(real(l, dl) * real(l + 1, dl))
    turn = asinh(ell / nu)
    end_chi = min(turn + 80._dl * PI / nu, 80._dl)
    width = max(0.08_dl, 0.02_dl * turn)

    call append_geom(raw, n, 1.0e-8_dl, min(0.5_dl, end_chi), 180)
    call append_lin(raw, n, 1.0e-7_dl, end_chi, merge(600, 1200, quick))
    call append_lin(raw, n, max(1.0e-8_dl, turn - width), min(end_chi, turn + width), &
        merge(1200, 3000, quick))
    if (end_chi > turn) then
        call append_lin(raw, n, turn, end_chi, merge(800, 1800, quick))
    end if
    call sort_unique(raw, n, chi)
    deallocate(raw)
    end subroutine make_open_chi_grid


    subroutine make_closed_chi_grid(chi)
    real(dl), allocatable, intent(out) :: chi(:)
    real(dl), allocatable :: raw(:), endpoint(:)
    integer :: n, i, nedge

    allocate(raw(MAX_RAW_POINTS))
    n = 0
    call append_geom(raw, n, 1.0e-8_dl, 0.5_dl, 180)
    call append_lin(raw, n, 1.0e-7_dl, PI / 2._dl, merge(800, 1800, quick))
    nedge = merge(1000, 2400, quick)
    allocate(endpoint(nedge))
    call geom_values(endpoint, 1.0e-12_dl, 0.25_dl)
    do i = 1, nedge
        call append_value(raw, n, PI / 2._dl - endpoint(i))
    end do
    deallocate(endpoint)
    call sort_unique(raw, n, chi)
    deallocate(raw)
    end subroutine make_closed_chi_grid


    subroutine append_lin(raw, n, a, b, m)
    real(dl), intent(inout) :: raw(:)
    integer, intent(inout) :: n
    real(dl), intent(in) :: a, b
    integer, intent(in) :: m
    integer :: i
    real(dl) :: x

    if (m <= 1) then
        call append_value(raw, n, a)
        return
    end if
    do i = 1, m
        x = a + (b - a) * real(i - 1, dl) / real(m - 1, dl)
        call append_value(raw, n, x)
    end do
    end subroutine append_lin


    subroutine append_geom(raw, n, a, b, m)
    real(dl), intent(inout) :: raw(:)
    integer, intent(inout) :: n
    real(dl), intent(in) :: a, b
    integer, intent(in) :: m
    real(dl), allocatable :: vals(:)
    integer :: i

    if (b <= a) then
        call append_value(raw, n, max(a, b))
        return
    end if
    allocate(vals(m))
    call geom_values(vals, a, b)
    do i = 1, m
        call append_value(raw, n, vals(i))
    end do
    deallocate(vals)
    end subroutine append_geom


    subroutine geom_values(vals, a, b)
    real(dl), intent(out) :: vals(:)
    real(dl), intent(in) :: a, b
    integer :: i, m
    real(dl) :: loga, logb

    m = size(vals)
    if (m == 1) then
        vals(1) = a
        return
    end if
    loga = log(a)
    logb = log(b)
    do i = 1, m
        vals(i) = exp(loga + (logb - loga) * real(i - 1, dl) / real(m - 1, dl))
    end do
    end subroutine geom_values


    subroutine append_value(raw, n, x)
    real(dl), intent(inout) :: raw(:)
    integer, intent(inout) :: n
    real(dl), intent(in) :: x

    if (x <= 0._dl) return
    if (n >= size(raw)) error stop "increase MAX_RAW_POINTS"
    n = n + 1
    raw(n) = x
    end subroutine append_value


    subroutine sort_unique(raw, n, values)
    real(dl), intent(inout) :: raw(:)
    integer, intent(in) :: n
    real(dl), allocatable, intent(out) :: values(:)
    integer :: i, m
    real(dl), allocatable :: work(:)

    allocate(work(n))
    work = raw(1:n)
    call quicksort(work, 1, n)
    m = 1
    do i = 2, n
        if (work(i) > work(m) * (1._dl + 1.0e-13_dl) + 1.0e-14_dl) then
            m = m + 1
            work(m) = work(i)
        end if
    end do
    allocate(values(m))
    values = work(1:m)
    deallocate(work)
    end subroutine sort_unique


    recursive subroutine quicksort(a, left, right)
    real(dl), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(dl) :: pivot, tmp

    if (left >= right) return
    i = left
    j = right
    pivot = a((left + right) / 2)
    do
        do while (a(i) < pivot)
            i = i + 1
        end do
        do while (pivot < a(j))
            j = j - 1
        end do
        if (i <= j) then
            tmp = a(i)
            a(i) = a(j)
            a(j) = tmp
            i = i + 1
            j = j - 1
        end if
        if (i > j) exit
    end do
    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
    end subroutine quicksort


    subroutine recursive_spline_error(l, k, nu, x, y, peak, max_err, expected_err, imax)
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu
    real(dl), intent(in) :: x(:), y(:), peak
    real(dl), intent(out) :: max_err, expected_err
    integer, intent(out) :: imax
    integer :: n, m, i, j
    real(dl), allocatable :: xk(:), yk(:), y2(:)
    real(dl) :: ys, err, max_phase, wave_number

    n = size(x)
    m = (n + 4) / 5
    if (mod(n - 1, 5) /= 0) m = m + 1
    allocate(xk(m), yk(m), y2(m))
    j = 0
    do i = 1, n, 5
        j = j + 1
        xk(j) = x(i)
        yk(j) = y(i)
    end do
    if (xk(j) /= x(n)) then
        j = j + 1
        xk(j) = x(n)
        yk(j) = y(n)
    end if
    m = j
    call natural_spline(xk(1:m), yk(1:m), y2(1:m))

    wave_number = spline_wave_number(l, k, nu)
    max_phase = 0._dl
    do i = 1, m - 1
        max_phase = max(max_phase, wave_number * (xk(i + 1) - xk(i)))
    end do
! Natural cubic interpolation has an O((k h)^4) envelope for resolved smooth
! oscillations. If the every-fifth knot grid is deliberately under-resolving a
! high-frequency closed mode, use this scale so the check flags unexpected
! spikes rather than ordinary interpolation undersampling.
    expected_err = max(SPLINE_TARGET_ERR, 2.0e-2_dl * max_phase**4)

    max_err = -1._dl
    imax = 1
    do i = 1, n
        ys = spline_value(xk(1:m), yk(1:m), y2(1:m), x(i))
        err = abs(ys - y(i)) / peak
        if (err > max_err) then
            max_err = err
            imax = i
        end if
    end do
    deallocate(xk, yk, y2)
    end subroutine recursive_spline_error


    real(dl) function spline_wave_number(l, k, nu)
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu

    if (k == 1) then
        spline_wave_number = max(1._dl, nu)
    else
        spline_wave_number = max(1._dl, sqrt(nu**2 + real(l, dl)**2))
    end if
    end function spline_wave_number


    subroutine natural_spline(x, y, y2)
    real(dl), intent(in) :: x(:), y(:)
    real(dl), intent(out) :: y2(:)
    integer :: n, i, k
    real(dl), allocatable :: u(:)
    real(dl) :: sig, p

    n = size(x)
    allocate(u(n))
    y2(1) = 0._dl
    u(1) = 0._dl
    do i = 2, n - 1
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
        p = sig * y2(i - 1) + 2._dl
        y2(i) = (sig - 1._dl) / p
        u(i) = (6._dl * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
            (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do
    y2(n) = 0._dl
    do k = n - 1, 1, -1
        y2(k) = y2(k) * y2(k + 1) + u(k)
    end do
    deallocate(u)
    end subroutine natural_spline


    real(dl) function spline_value(x, y, y2, xval)
    real(dl), intent(in) :: x(:), y(:), y2(:), xval
    integer :: klo, khi, k
    real(dl) :: h, a, b

    klo = 1
    khi = size(x)
    do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (x(k) > xval) then
            khi = k
        else
            klo = k
        end if
    end do
    h = x(khi) - x(klo)
    if (h <= 0._dl) error stop "bad spline knots"
    a = (x(khi) - xval) / h
    b = (xval - x(klo)) / h
    spline_value = a * y(klo) + b * y(khi) + ((a**3 - a) * y2(klo) + (b**3 - b) * y2(khi)) * h**2 / 6._dl
    end function spline_value


    subroutine insert_record(records, rec)
    type(worst_record), intent(inout) :: records(:)
    type(worst_record), intent(in) :: rec
    integer :: i, j

    do i = 1, size(records)
        if (rec%err > records(i)%err) then
            do j = size(records), i + 1, -1
                records(j) = records(j - 1)
            end do
            records(i) = rec
            exit
        end if
    end do
    end subroutine insert_record


    subroutine print_records(title, records)
    character(len=*), intent(in) :: title
    type(worst_record), intent(in) :: records(:)
    integer :: i

    write(*, '(/, a)') trim(title)
    write(*, '(a)') "rank geometry L K nu alpha chi metric err expected peak index"
    do i = 1, size(records)
        if (records(i)%err < 0._dl) exit
        write(*, '(i3, 1x, a8, 1x, i6, 1x, i2, 1x, es13.5, 1x, es10.3, 1x, es13.5, &
            1x, es10.3, 1x, es10.3, 1x, es10.3, 1x, es10.3, 1x, i8)') &
            i, records(i)%geometry, records(i)%l, records(i)%k, records(i)%nu, records(i)%alpha, &
            records(i)%chi, records(i)%metric, records(i)%err, records(i)%expected, records(i)%peak, records(i)%index
    end do
    end subroutine print_records


    subroutine run_timing_benchmarks()
    write(*, '(/, a)') "Timing benchmarks by predicted gate region"
    write(*, '(a)') "label region L K nu alpha selected repeats phi_olver_cpu/eval " // &
        "phi_recurs_cpu/eval phi_olver_wall/eval phi_recurs_wall/eval ratio_cpu ratio_wall"
    call timing_benchmark_case("open high-L accepted", 6000, -1, 0.12_dl * 6000._dl, want_fallback=.false.)
    call timing_benchmark_case("open high-L fallback", 6000, -1, 0.05_dl * 6000._dl, want_fallback=.true.)
    call timing_benchmark_case("closed high-L accepted", 6000, 1, 6200._dl, want_fallback=.false.)
    call timing_benchmark_case("closed high-L fallback", 6000, 1, 6100._dl, want_fallback=.true.)
    call timing_benchmark_case("closed mid-L accepted", 500, 1, 700._dl, want_fallback=.false.)
    call timing_benchmark_case("closed mid-L fallback", 500, 1, 600._dl, want_fallback=.true.)
    end subroutine run_timing_benchmarks


    subroutine timing_benchmark_case(label, l, k, nu, want_fallback)
    character(len=*), intent(in) :: label
    integer, intent(in) :: l, k
    real(dl), intent(in) :: nu
    logical, intent(in) :: want_fallback
    real(dl), allocatable :: chi(:)
    integer :: n, selected, repeats, repeat, i
    real(dl) :: t0_cpu, t1_cpu, t0_wall, t1_wall
    real(dl) :: approx_cpu_case, recurs_cpu_case, approx_wall_case, recurs_wall_case
    real(dl) :: approx_eval_cpu, recurs_eval_cpu, approx_eval_wall, recurs_eval_wall
    real(dl) :: sink
    character(len=8) :: region

    if (k == -1) then
        call make_open_chi_grid(l, nu, chi)
    else
        call make_closed_chi_grid(chi)
    end if
    n = size(chi)
    selected = 0
    do i = 1, n
        if (predicts_recursive_fallback(l, k, nu, chi(i)) .eqv. want_fallback) selected = selected + 1
    end do
    if (selected == 0) then
        deallocate(chi)
        return
    end if

    repeats = max(1, min(50, 200000 / selected))
    if (quick) repeats = max(1, min(10, 50000 / selected))
    region = merge("fallback", "fast    ", want_fallback)

    sink = 0._dl
    call cpu_time(t0_cpu)
    t0_wall = omp_get_wtime()
    do repeat = 1, repeats
        !$omp parallel do default(shared) private(i) reduction(+:sink) schedule(static)
        do i = 1, n
            if (predicts_recursive_fallback(l, k, nu, chi(i)) .eqv. want_fallback) then
                sink = sink + phi_olver(l, k, nu, chi(i))
            end if
        end do
        !$omp end parallel do
    end do
    t1_wall = omp_get_wtime()
    call cpu_time(t1_cpu)
    approx_cpu_case = t1_cpu - t0_cpu
    approx_wall_case = t1_wall - t0_wall

    call cpu_time(t0_cpu)
    t0_wall = omp_get_wtime()
    do repeat = 1, repeats
        !$omp parallel do default(shared) private(i) reduction(+:sink) schedule(static)
        do i = 1, n
            if (predicts_recursive_fallback(l, k, nu, chi(i)) .eqv. want_fallback) then
                sink = sink + phi_recurs(l, k, nu, chi(i))
            end if
        end do
        !$omp end parallel do
    end do
    t1_wall = omp_get_wtime()
    call cpu_time(t1_cpu)
    recurs_cpu_case = t1_cpu - t0_cpu
    recurs_wall_case = t1_wall - t0_wall

    approx_eval_cpu = approx_cpu_case / real(selected * repeats, dl)
    recurs_eval_cpu = recurs_cpu_case / real(selected * repeats, dl)
    approx_eval_wall = approx_wall_case / real(selected * repeats, dl)
    recurs_eval_wall = recurs_wall_case / real(selected * repeats, dl)

    write(*, '(a, 1x, a8, 1x, i6, 1x, i2, 1x, es12.4, 1x, es10.3, 1x, i8, 1x, i4, &
        1x, es12.4, 1x, es12.4, 1x, es12.4, 1x, es12.4, 1x, f9.3, 1x, f9.3)') &
        trim(label), region, l, k, nu, nu / real(l, dl), selected, repeats, &
        approx_eval_cpu, recurs_eval_cpu, approx_eval_wall, recurs_eval_wall, &
        recurs_eval_cpu / max(approx_eval_cpu, tiny(1._dl)), &
        recurs_eval_wall / max(approx_eval_wall, tiny(1._dl))
    if (abs(sink) > huge(1._dl)) write(*, *) sink

    deallocate(chi)
    end subroutine timing_benchmark_case

    end program phi_olver_gate_validation
