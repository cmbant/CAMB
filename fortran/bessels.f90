    !CAMB flat spherical Bessel function routines
    !June 2026: updated bjl, accurate to peak-normalized fraction <7e-6 at L>=28
    !           (max ~9e-6 at BJL_RECURRENCE_MAX_L+1, ~1.5e-6 at high L), pre-splining.
    !           Spline table accurate to 2e-4 in tail, better round peak.
    ! Tolerance note: first positive peak of spherical Bessel function j_l(x), l >= 1:
    ! nu = l + 1/2
    ! x_peak ~ nu + 0.80861652*nu^(1/3) - 0.23669965*nu^(-1/3) - 0.20430105*nu^(-1)
    ! A_peak ~ 0.845843*nu^(-5/6)*(1 - 0.55424*nu^(-2/3) + 0.25865*nu^(-4/3))
    ! where A_peak = j_l(x_peak); max relative error < 0.7%, better at high l.
    module FlatBessels
    use Precision
    use results
    use RangeUtils
    use MpiUtils
    use Interpolation, only: spline, SPLINE_DANGLE
    implicit none
    private

    ! Packed interval coefficients for the cubic Bessel interpolation in Horner form.
    ! The first index stores:
    !   1 -> base value y1
    !   2 -> Horner coefficient c1
    !   3 -> Horner coefficient c2
    !   4 -> Horner coefficient c3

    real(dl), dimension(:,:,:), allocatable ::  bessel_horner

    integer  num_xx, kmaxfile, max_ix
    Type(lSamples), save :: file_l
    ! parameter for working out where the flat Bessel functions are small
    ! Should increase for higher accuracy
    ! For x = l-delta below the turning point, j_l is suppressed roughly as
    ! exp[-(2*sqrt(2)/3)*d**(3/2)/sqrt(l)].  Requiring ~1e-4 of peak
    ! gives delta ~ 4.2*l**(1/3), with a small safety margin.
    real(dl), parameter :: bjl_pre_peak_start_factor  = 4.2_dl
    integer, parameter :: BJL_RECURRENCE_MAX_L = 25
    real(dl) file_acc, file_bessel_boost
    real(dl) bessel_xmaxfile

    type(TRanges), save:: BessRanges

    public bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor
    public bjl, Bessels_Free, airy_fast, BJL_RECURRENCE

    contains

    subroutine InitSpherBessels(lSamp, CP, max_bessels_l_index, max_bessels_etak)
    Type(lSamples) lSamp
    Type(CAMBParams) :: CP
    integer, intent(in) :: max_bessels_l_index
    real(dl), intent(in) :: max_bessels_etak

    integer :: new_kmaxfile, old_num_xx
    real(dl) :: requested_xmax
    logical :: same_lsamp_and_acc

    new_kmaxfile = int(max_bessels_etak) + 1
    if (do_bispectrum) new_kmaxfile = 2*new_kmaxfile
    requested_xmax = real(new_kmaxfile, dl)

! See if already loaded with enough and correct lSamp%l values, accuracy,
! and l coverage.
    same_lsamp_and_acc = .false.
    if (allocated(bessel_horner) .and. lSamp%nl <= file_l%nl) then
        if (all(file_l%l(1:lSamp%nl) == lSamp%l(1:lSamp%nl)) .and. &
            max_bessels_l_index <= max_ix .and. &
            abs(CP%Accuracy%BesselBoost - file_bessel_boost) < 1d-2 .and. &
            abs(CP%Accuracy%BesselBoost*CP%Accuracy%AccuracyBoost - file_acc) < 1d-2) then
            same_lsamp_and_acc = .true.
        end if
    end if

    if (same_lsamp_and_acc) then
        if (requested_xmax <= bessel_xmaxfile) return

        old_num_xx = num_xx
        call ExtendBessels(CP, old_num_xx, requested_xmax)

        if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Extended Bessels'
        return
    end if

! Haven't made them before, or l/accuracy requirements changed.
    kmaxfile = new_kmaxfile
    max_ix = min(max_bessels_l_index, lSamp%nl)

    call GenerateBessels(lSamp, CP, requested_xmax)

    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Calculated Bessels'

    end subroutine InitSpherBessels

    elemental subroutine bjl_deriv(l, x, jl, djl)
    integer, intent(in) :: l
    real(dl), intent(in) :: x, jl
    real(dl), intent(out) :: djl

    real(dl) :: jm1

    if (x == 0._dl) then
        if (l == 1) then
            djl = 1._dl/3._dl
        else
            djl = 0._dl
        end if
    else if (l == 0) then
        call bjl(1, x, jm1)
        djl = -jm1
    else
        call bjl(l - 1, x, jm1)
        djl = jm1 - real(l + 1, dl)*jl/x
    end if

    end subroutine bjl_deriv

    subroutine GenerateBessels(lSamp, CP, requested_xmax)
    Type(lSamples) lSamp
    Type(CAMBParams) :: CP
    real(dl), intent(in) :: requested_xmax

    integer j
    integer, parameter :: cut_max_l = 25
    real(dl), parameter :: cut(1:cut_max_l) = (/ &
        0.000000_dl, 0.000000_dl, 0.063316_dl, 0.208916_dl, &
        0.448187_dl, 0.769530_dl, 1.158338_dl, 1.601506_dl, 2.088363_dl, &
        2.610505_dl, 3.161369_dl, 3.735827_dl, 4.329841_dl, 4.940211_dl, &
        5.564374_dl, 6.200262_dl, 6.846189_dl, 7.500771_dl, 8.162861_dl, &
        8.831502_dl, 9.505890_dl,10.185346_dl,10.869292_dl,11.557230_dl, &
        12.248737_dl /)

    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'Generating flat Bessels...'

    file_l = lSamp

    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'x_max bessels', requested_xmax

    call BessRanges%Init()

    call BessRanges%Add_delta(0._dl, 1._dl, 0.01_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(1._dl, 5._dl, 0.1_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(5._dl, 25._dl, 0.2_dl/CP%Accuracy%BesselBoost)

    file_bessel_boost = CP%Accuracy%BesselBoost
    file_acc = file_bessel_boost*CP%Accuracy%AccuracyBoost

    call BessRanges%Add_delta(25._dl, 150._dl, 0.5_dl/file_acc)
    call BessRanges%Add_delta(150._dl, requested_xmax, 0.7_dl/file_acc) ! 2e-4 accuracy

    call BessRanges%GetArray(.false.)
    num_xx = BessRanges%npoints

    bessel_xmaxfile = BessRanges%points(num_xx)
    kmaxfile = ceiling(bessel_xmaxfile)

    if (allocated(bessel_horner)) then
        if (any(ubound(bessel_horner) < [4, num_xx - 1, max_ix])) deallocate(bessel_horner)
    end if
    if (.not. allocated(bessel_horner)) then
        allocate(bessel_horner(1:4, 1:num_xx-1, 1:max_ix))
    end if

!$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC)
    do j = 1, max_ix
        block
            real(dl) :: h2over6, y0, y1, d0, d1, xlim, d_end
            real(dl) :: knot_vals(num_xx), spline_y2(num_xx)
            integer :: min_ix, i

            xlim = max(lSamp%l(j) - bjl_pre_peak_start_factor*lSamp%l(j)**(1._dl/3._dl) - 1, &
                cut(min(cut_max_l, lSamp%l(j))))

            min_ix = max(1, BessRanges%IndexOf(xlim) - 1)

            knot_vals(1:max(min_ix-1, 1)) = 0
            do i = min_ix, num_xx
                call bjl(lSamp%l(j), BessRanges%points(i), knot_vals(i))
            end do

            call bjl_deriv(lSamp%l(j), BessRanges%points(num_xx), knot_vals(num_xx), d_end)
            call spline(BessRanges%points, knot_vals, num_xx, SPLINE_DANGLE, d_end, spline_y2)

            bessel_horner(:, 1:max(min_ix-1, 1), j) = 0

            do i = max(1, min_ix-1), num_xx-1
                y0 = knot_vals(i)
                y1 = knot_vals(i+1)
                d0 = spline_y2(i)
                d1 = spline_y2(i+1)

                h2over6 = (BessRanges%points(i+1) - BessRanges%points(i))**2/6

                bessel_horner(1, i, j) = y1
                bessel_horner(2, i, j) = y0 - y1 - h2over6*(d0 + 2*d1)
                bessel_horner(3, i, j) = 3*h2over6*d1
                bessel_horner(4, i, j) = y0 - y1 - &
                    bessel_horner(2, i, j) - bessel_horner(3, i, j)
            end do
        end block
    end do
!$OMP END PARALLEL DO

    end subroutine GenerateBessels

    subroutine ExtendBessels(CP, old_num_xx, requested_xmax)
! Extend bessel_horner to the new x range without invalidating the old
! Horner coefficients.
!
! The final [150, xmax] range is enlarged so that its endpoint is an integer
! multiple of the current final interval beyond 150.  Rebuilding BessRanges
! with this snapped endpoint gives the same final-range spacing as before,
! so the old grid remains a prefix of the new grid.
!
! We recompute an overlap window and overwrite the latter part of the old
! intervals in that window.  This avoids stitching a newly splined interval
! onto an old interval whose right endpoint was previously a spline boundary.
    Type(CAMBParams) :: CP
    integer, intent(in) :: old_num_xx
    real(dl), intent(in) :: requested_xmax

    integer :: j, ext_start_ix, n_ext
    integer :: old_final_n, new_final_n
    integer, parameter :: JUNCTION_N = 10

    real(dl) :: old_xmax, old_dx, snapped_xmax
    real(dl), parameter :: final_range_start = 150._dl

    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'Extending flat Bessels to x_max', requested_xmax

! Need a valid existing final interval to preserve.
    if (old_num_xx < 2) then
        call GenerateBessels(file_l, CP, requested_xmax)
        return
    end if

    old_xmax = BessRanges%points(old_num_xx)
    old_dx   = BessRanges%points(old_num_xx) - BessRanges%points(old_num_xx - 1)

! This extension logic assumes the existing table already reaches the final
! sparse range.  If not, regenerate safely.
    if (old_xmax <= final_range_start .or. old_dx <= 0._dl) then
        call GenerateBessels(file_l, CP, requested_xmax)
        return
    end if

    old_final_n = nint((old_xmax - final_range_start)/old_dx)
    new_final_n = ceiling((requested_xmax - final_range_start)/old_dx)

    if (new_final_n <= old_final_n) new_final_n = old_final_n + 1

    snapped_xmax = final_range_start + real(new_final_n, dl)*old_dx

    call BessRanges%Init()

    call BessRanges%Add_delta(0._dl,  1._dl,   0.01_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(1._dl,  5._dl,   0.1_dl /CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(5._dl,  25._dl,  0.2_dl /CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(25._dl, final_range_start, 0.5_dl/file_acc)
    call BessRanges%Add_delta(final_range_start, snapped_xmax, old_dx)

    call BessRanges%GetArray(.false.)
    num_xx = BessRanges%npoints

    bessel_xmaxfile = BessRanges%points(num_xx)
    kmaxfile = ceiling(bessel_xmaxfile)

    ext_start_ix = max(1, old_num_xx - JUNCTION_N)
    n_ext = num_xx - ext_start_ix + 1

! Reallocate bessel_horner, keeping all old intervals initially intact.
    block
        real(dl), allocatable :: old_horner(:,:,:)

        call move_alloc(bessel_horner, old_horner)

        allocate(bessel_horner(1:4, 1:num_xx-1, 1:max_ix))

        bessel_horner(:, 1:old_num_xx-1, 1:max_ix) = &
            old_horner(:, 1:old_num_xx-1, 1:max_ix)
    end block

!$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC)
    do j = 1, max_ix
        block
            real(dl) :: ext_x(n_ext), ext_y(n_ext), ext_y2(n_ext)
            real(dl) :: h2over6, y0, y1, d0, d1, d_end
            integer :: i, store_ix, first_new_i, overwrite_start

            do i = 1, n_ext
                ext_x(i) = BessRanges%points(ext_start_ix + i - 1)
            end do

            do i = 1, n_ext
                call bjl(file_l%l(j), ext_x(i), ext_y(i))
            end do

            call bjl_deriv(file_l%l(j), ext_x(n_ext), ext_y(n_ext), d_end)
            call spline(ext_x, ext_y, n_ext, SPLINE_DANGLE, d_end, ext_y2)

            ! Local index of the first strictly new interval.  Its global interval
            ! index is old_num_xx.
            first_new_i = old_num_xx - ext_start_ix + 1

            ! Move the stitch into the overlap, away from the old endpoint where
            ! the previous spline had a boundary condition.  This does not make the
            ! combined table exactly the same as a full re-spline, but it removes
            ! the large artificial kink at old_num_xx.
            overwrite_start = max(1, first_new_i - max(1, JUNCTION_N/2))

            do i = overwrite_start, n_ext - 1
                store_ix = ext_start_ix + i - 1

                y0 = ext_y(i)
                y1 = ext_y(i+1)
                d0 = ext_y2(i)
                d1 = ext_y2(i+1)

                h2over6 = (ext_x(i+1) - ext_x(i))**2/6

                bessel_horner(1, store_ix, j) = y1
                bessel_horner(2, store_ix, j) = y0 - y1 - h2over6*(d0 + 2*d1)
                bessel_horner(3, store_ix, j) = 3*h2over6*d1
                bessel_horner(4, store_ix, j) = y0 - y1 - &
                    bessel_horner(2, store_ix, j) - bessel_horner(3, store_ix, j)
            end do
        end block
    end do
!$OMP END PARALLEL DO

    end subroutine ExtendBessels

    subroutine Bessels_Free

    if (allocated(bessel_horner)) deallocate(bessel_horner)
    if (allocated(file_l%l)) deallocate(file_l%l)
    file_l%nl=0

    call BessRanges%Free()

    end subroutine Bessels_Free

! Optimized spherical Bessel wrapper.
! Strategy:
!   For low L use (v accurate and still fast) recursive result
!   Elsewhere use a two-term corrected uniform Airy asymptotic
!   in the transition bands, using fast approximations elsewhere
!   where they are accurate.
!
! accurate to peak-normalized fraction <7e-6 at L>=28 (worst at the
! Airy/Debye band edges), max ~9e-6 at BJL_RECURRENCE_MAX_L+1.

    ELEMENTAL SUBROUTINE BJL(L, X, JL)
    ! Optimized spherical Bessel j_l(x).
    !
    ! Branch order:
    !   low l           : explicit formulas/small-x series
    !   very small x    : zero for l >= 7
    !   deep pre-peak   : exponentially small ascending expansion
    !   far x           : large-x oscillatory expansion
    !   transition band : recurrence only for moderate l where needed;
    !                     otherwise Airy shoulders, peak polynomial,
    !                     or side asymptotics
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    REAL(dl), INTENT(IN) :: X
    REAL(dl), INTENT(OUT) :: JL

    REAL(dl), PARAMETER :: LN2=0.6931471805599453094_dl
    REAL(dl), PARAMETER :: ONEMLN2=0.30685281944005469058277_dl
    REAL(dl), PARAMETER :: PID2=1.5707963267948966192313217_dl
    REAL(dl), PARAMETER :: PID4=0.78539816339744830961566084582_dl
    REAL(dl), PARAMETER :: ROOTPI12 = 21.269446210866192327578_dl
    REAL(dl), PARAMETER :: GAMMA1 = 2.6789385347077476336556_dl
    REAL(dl), PARAMETER :: GAMMA2 = 1.3541179394264004169452_dl

    REAL(dl) :: AX, AX2, NU, NU2, BETA, BETA2, COSB
    REAL(dl) :: SX, SX2, COTB, COT3B, COT6B, SECB, SEC2B
    REAL(dl) :: EXPTERM, L3, ETA, ETA_AIRY_LOW, ETA_AIRY_HIGH

    ! Gates below are in eta = (x - nu)/nu^(1/3) units (so tau ~ -1.26*eta
    ! in the Airy band, nearly independent of l).
    ! Note: bjl_uniform_airy_fast uses a polynomial turning-point mapping fitted
    ! for the shoulder bands below; widening BJL_AIRY_ETA_LOW/HIGH beyond
    ! [-2.4, 3.85] (or the u caps beyond [-0.26, 0.42]) requires refitting its
    ! zeta/u and correction coefficients.
    real(dl), parameter :: BJL_RECURRENCE_ETA_LOW  = -5.0_dl
    real(dl), parameter :: BJL_RECURRENCE_ETA_HIGH =  5.6_dl
    real(dl), parameter :: BJL_AIRY_ETA_LOW  = -2.4_dl
    real(dl), parameter :: BJL_PEAK_ETA_LOW  = -0.65_dl
    real(dl), parameter :: BJL_PEAK_ETA_HIGH =  0.65_dl
    real(dl), parameter :: BJL_AIRY_ETA_HIGH =  3.85_dl
    ! Caps on u = x/nu - 1 so the Airy-branch zeta/u polynomial fit stays in
    ! range; they only bind for l <= 28 where the eta limits would exceed them.
    real(dl), parameter :: BJL_AIRY_U_LOW  = 0.26_dl
    real(dl), parameter :: BJL_AIRY_U_HIGH = 0.42_dl
    IF (L < 0) THEN
        ERROR STOP 'Can not evaluate Spherical Bessel Function with index l<0'
    END IF

    AX = ABS(X)
    AX2 = AX**2

    ! Low orders: explicit formulae, with small-x series to avoid cancellation.
    IF (L < 7) THEN
        SELECT CASE (L)
        CASE (0)
            IF (AX < 1.0E-1_dl) THEN
                JL = 1.0_dl - AX2/6.0_dl*(1.0_dl - AX2/20.0_dl)
            ELSE
                JL = SIN(AX)/AX
            END IF
        CASE (1)
            IF (AX < 2.0E-1_dl) THEN
                JL = AX/3.0_dl*(1.0_dl - AX2/10.0_dl*(1.0_dl - AX2/28.0_dl))
            ELSE
                JL = (SIN(AX)/AX - COS(AX))/AX
            END IF
        CASE (2)
            IF (AX < 3.0E-1_dl) THEN
                JL = AX2/15.0_dl*(1.0_dl - AX2/14.0_dl*(1.0_dl - AX2/36.0_dl))
            ELSE
                JL = (-3.0_dl*COS(AX)/AX - SIN(AX)*(1.0_dl - 3.0_dl/AX2))/AX
            END IF
        CASE (3)
            IF (AX < 4.0E-1_dl) THEN
                JL = AX*AX2/105.0_dl*(1.0_dl - AX2/18.0_dl*(1.0_dl - AX2/44.0_dl))
            ELSE
                JL = (COS(AX)*(1.0_dl - 15.0_dl/AX2) - SIN(AX)*(6.0_dl - 15.0_dl/AX2)/AX)/AX
            END IF
        CASE (4)
            IF (AX < 6.0E-1_dl) THEN
                JL = AX2**2/945.0_dl*(1.0_dl - AX2/22.0_dl*(1.0_dl - AX2/52.0_dl))
            ELSE
                JL = (SIN(AX)*(1.0_dl - (45.0_dl - 105.0_dl/AX2)/AX2) &
                    + COS(AX)*(10.0_dl - 105.0_dl/AX2)/AX)/AX
            END IF
        CASE (5)
            IF (AX < 1.0_dl) THEN
                JL = AX2**2*AX/10395.0_dl*(1.0_dl - AX2/26.0_dl*(1.0_dl - AX2/60.0_dl))
            ELSE
                JL = (SIN(AX)*(15.0_dl - (420.0_dl - 945.0_dl/AX2)/AX2)/AX &
                    - COS(AX)*(1.0_dl - (105.0_dl - 945.0_dl/AX2)/AX2))/AX
            END IF
        CASE DEFAULT
            IF (AX < 1.0_dl) THEN
                JL = AX2**3/135135.0_dl*(1.0_dl - AX2/30.0_dl*(1.0_dl - AX2/68.0_dl))
            ELSE
                JL = (SIN(AX)*(-1.0_dl + (210.0_dl - (4725.0_dl - 10395.0_dl/AX2)/AX2)/AX2) &
                    + COS(AX)*(-21.0_dl + (1260.0_dl - 10395.0_dl/AX2)/AX2)/AX)/AX
            END IF
        END SELECT

        IF (X < 0.0_dl .AND. MOD(L,2) /= 0) JL = -JL
        RETURN
    END IF

    NU = REAL(L, dl) + 0.5_dl
    NU2 = NU**2

    IF (AX < 1.0E-40_dl) THEN
        ! Very small x: j_l(x) is negligible here for l >= 7.
        JL = 0.0_dl
    ELSE IF ((AX2/REAL(L, dl)) < 5.0E-1_dl) THEN
        ! Deep pre-peak: exponentially small ascending expansion.
        JL = EXP(REAL(L,dl)*LOG(AX/NU) - LN2 + NU*ONEMLN2 &
            - (1.0_dl - (1.0_dl - 3.5_dl/NU2)/NU2/30.0_dl)/12.0_dl/NU) &
            /NU*(1.0_dl - AX2/(4.0_dl*NU + 4.0_dl)*(1.0_dl - AX2/(8.0_dl*NU + 16.0_dl) &
            *(1.0_dl - AX2/(12.0_dl*NU + 36.0_dl))))
    ELSE IF ((REAL(L, dl)**2/AX) < 1.2_dl) THEN
        ! Far past the peak: oscillatory large-x expansion.
        ! At the q = l^2/x = 1.2 gate the truncation error is < ~6e-6 of peak
        ! for l <= 15, < 3e-6 for l >= 25, falling rapidly with l and with x.
        BETA = AX - PID2*REAL(L + 1, dl)
        JL = (COS(BETA)*(1.0_dl - (NU2 - 0.25_dl)*(NU2 - 2.25_dl)/8.0_dl/AX2 &
            *(1.0_dl - (NU2 - 6.25_dl)*(NU2 - 12.25_dl)/48.0_dl/AX2)) &
            - SIN(BETA)*(NU2 - 0.25_dl)/2.0_dl/AX &
            *(1.0_dl - (NU2 - 2.25_dl)*(NU2 - 6.25_dl)/24.0_dl/AX2 &
            *(1.0_dl - (NU2 - 12.25_dl)*(NU2 - 20.25_dl)/80.0_dl/AX2)))/AX
    ELSE IF (AX > 1.5_dl*NU .AND. (L > BJL_RECURRENCE_MAX_L .OR. AX > 2.5_dl*NU)) THEN
        ! Well past the turning point every classification below lands in the
        ! post-peak Debye branch, so skip the slow nu**0.325 eta normalization.
        ! For l >= 26, ax > 1.5*nu gives u = ax/nu - 1 > 0.42, past the Airy
        ! gate; for 7 <= l <= 25, ax > 2.5*nu gives eta > 5.6, past the
        ! recurrence band.
        CALL BJL_POSTPEAK_DEBYE(NU, AX, JL)
    ELSE
        ! Turning-point neighbourhood: classify by normalized distance from peak.
        L3 = NU**(1.0_dl/3.0_dl)
        ETA = (AX - NU)/L3

        ! Moderate orders only need recurrence in the broad transition band.
        ! Deep pre-peak and far-x cases have already returned via cheaper branches.
        IF (L <= BJL_RECURRENCE_MAX_L .AND. &
            ETA >= BJL_RECURRENCE_ETA_LOW .AND. ETA <= BJL_RECURRENCE_ETA_HIGH) THEN
            CALL BJL_RECURRENCE(L, X, JL)
            RETURN
        END IF

        ETA_AIRY_LOW  = MAX(BJL_AIRY_ETA_LOW, -BJL_AIRY_U_LOW*NU/L3)
        ETA_AIRY_HIGH = MIN(BJL_AIRY_ETA_HIGH, BJL_AIRY_U_HIGH*NU/L3)

        IF ((ETA >= ETA_AIRY_LOW .AND. ETA < BJL_PEAK_ETA_LOW) .OR. &
            (ETA > BJL_PEAK_ETA_HIGH .AND. ETA <= ETA_AIRY_HIGH)) THEN
            ! Shoulder bands: corrected uniform Airy approximation.
            CALL BJL_UNIFORM_AIRY_FAST(L, X, L3*L3, JL)
            RETURN
        ELSE IF (ETA < ETA_AIRY_LOW) THEN
            ! Below the peak but outside the exponentially tiny region.
            COSB = NU/AX
            SX = SQRT(NU2 - AX2)
            COTB = NU/SX
            SECB = AX/NU
            BETA = LOG(COSB + SX/AX)
            COT3B = COTB**3
            COT6B = COT3B**2
            SEC2B = SECB**2
            EXPTERM = ((2.0_dl + 3.0_dl*SEC2B)*COT3B/24.0_dl &
                - ((4.0_dl + SEC2B)*SEC2B*COT6B/16.0_dl &
                + ((16.0_dl - (1512.0_dl + (3654.0_dl + 375.0_dl*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.0_dl &
                + (32.0_dl + (288.0_dl + (232.0_dl + 13.0_dl*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.0_dl/NU) &
                *COT6B/NU)/NU)/NU
            JL = SQRT(COTB*COSB)/(2.0_dl*NU)*EXP(-NU*BETA + NU/COTB - EXPTERM)
        ELSE IF (ETA > ETA_AIRY_HIGH) THEN
            ! Above the peak: oscillatory post-peak asymptotic form.
            CALL BJL_POSTPEAK_DEBYE(NU, AX, JL)
        ELSE
            ! Very close to the peak: polynomial transition expansion.
            BETA = AX - NU
            BETA2 = BETA**2
            SX = 6.0_dl/AX
            SX2 = SX**2
            SECB = SX**(1.0_dl/3.0_dl)
            SEC2B = SECB**2
            JL = (GAMMA1*SECB + BETA*GAMMA2*SEC2B &
                - (BETA2/18.0_dl - 1.0_dl/45.0_dl)*BETA*SX*SECB*GAMMA1 &
                - ((BETA2 - 1.0_dl)*BETA2/36.0_dl + 1.0_dl/420.0_dl)*SX*SEC2B*GAMMA2 &
                + (((BETA2/1620.0_dl - 7.0_dl/3240.0_dl)*BETA2 + 1.0_dl/648.0_dl)*BETA2 &
                - 1.0_dl/8100.0_dl)*SX2*SECB*GAMMA1 &
                + (((BETA2/4536.0_dl - 1.0_dl/810.0_dl)*BETA2 + 19.0_dl/11340.0_dl)*BETA2 &
                - 13.0_dl/28350.0_dl)*BETA*SX2*SEC2B*GAMMA2 &
                - ((((BETA2/349920.0_dl - 1.0_dl/29160.0_dl)*BETA2 + 71.0_dl/583200.0_dl)*BETA2 &
                - 121.0_dl/874800.0_dl)*BETA2 + 7939.0_dl/224532000.0_dl)*BETA*SX2*SX*SECB*GAMMA1) &
                *SQRT(SX)/ROOTPI12
        END IF
    END IF

    IF (X < 0.0_dl .AND. MOD(L,2) /= 0) JL = -JL

    END SUBROUTINE BJL


    elemental subroutine bjl_postpeak_debye(nu, ax, jl)
    ! Oscillatory Debye asymptotic above the turning point, valid above the
    ! Airy band (u = ax/nu - 1 > 0.42, or eta above the gates in BJL).
    ! Identical to the standard two-correction form but with
    !   sqrt(cotb*cosb)/nu = 1/sqrt(sx*ax),  nu/cotb = sx,
    ! exp(-expterm) expanded to second order (expterm < ~1e-3 over the
    ! accepted domain, so the expansion error is < ~1e-9), and
    ! nu*acos(cosb) = nu*(pi/2 - asin(cosb)) evaluated via a polynomial fit
    ! of asin(c)/c in c^2 on [0, 0.4975] (max |asin error| 3.7e-11) when in
    ! range, i.e. for ax > ~1.42*nu; acos is only called close to the
    ! turning point.
    implicit none
    real(dl), intent(in) :: nu, ax
    real(dl), intent(out) :: jl

    real(dl), parameter :: pid2 = 1.5707963267948966192313217_dl
    real(dl), parameter :: pid4 = 0.78539816339744830961566084582_dl
    real(dl), parameter :: a(0:10) = (/ &
        +1.00000000003854628e+00_dl, +1.66666648374737825e-01_dl, &
        +7.50014285476487963e-02_dl, +4.45995913792991625e-02_dl, &
        +3.10488804563286078e-02_dl, +1.64373788064744349e-02_dl, &
        +4.97790102349004449e-02_dl, -9.71772607183009063e-02_dl, &
        +2.46755256321489758e-01_dl, -2.76229262323224423e-01_dl, &
        +1.68861570948557360e-01_dl /)
    real(dl) :: nu2, sx, cotb, cot3b, cot6b, sec2b, trigarg, expterm
    real(dl) :: cb, s, p
    integer :: k

    nu2 = nu*nu
    sx = sqrt(ax*ax - nu2)

    cb = nu/ax
    s = cb*cb

    if (nu > real(BJL_RECURRENCE_MAX_L, dl) + 0.5_dl .and. ax >= 3.0_dl*nu) then
        ! In the far post-peak tail the first phase correction dominates the
        ! remaining Debye terms.  For L >= 26 and ax >= 3*nu this stays below
        ! 5e-6 peak-normalized error, improving rapidly with ax/nu.
        ! Here cb = nu/ax <= 1/3, so a short asin(cb)/cb series is accurate
        ! enough and faster than the broad-domain polynomial used below.
        p = 1.0_dl + s*(1.0_dl/6.0_dl + s*(3.0_dl/40.0_dl + s*(5.0_dl/112.0_dl &
            + s*(35.0_dl/1152.0_dl + s*(63.0_dl/2816.0_dl + s*(231.0_dl/13312.0_dl))))))
        trigarg = sx + nu*(cb*p - pid2)
        cotb = nu/sx
        cot3b = cotb*cotb*cotb
        sec2b = 1.0_dl/s
        trigarg = trigarg - pid4 - (2.0_dl + 3.0_dl*sec2b)*cot3b/(24.0_dl*nu)
        jl = cos(trigarg)/sqrt(sx*ax)
        return
    end if

    if (s < 0.4975_dl) then
        p = a(10)
        do k = 9, 0, -1
            p = p*s + a(k)
        end do
        trigarg = sx + nu*(cb*p - pid2)
    else
        trigarg = sx - nu*acos(cb)
    end if

    cotb = nu/sx
    sec2b = 1.0_dl/s
    cot3b = cotb*cotb*cotb
    cot6b = cot3b*cot3b

    trigarg = trigarg - pid4 &
        - ((2.0_dl + 3.0_dl*sec2b)*cot3b/24.0_dl &
        + (16.0_dl - (1512.0_dl + (3654.0_dl + 375.0_dl*sec2b)*sec2b)*sec2b) &
        *cot3b*cot6b/5760.0_dl/nu2)/nu
    expterm = ((4.0_dl + sec2b)*sec2b*cot6b/16.0_dl &
        - (32.0_dl + (288.0_dl + (232.0_dl + 13.0_dl*sec2b)*sec2b)*sec2b) &
        *sec2b*cot6b**2/128.0_dl/nu2)/nu2
    jl = (1.0_dl - expterm*(1.0_dl - 0.5_dl*expterm))*cos(trigarg)/sqrt(sx*ax)

    end subroutine bjl_postpeak_debye


    elemental subroutine bjl_uniform_airy_fast(l, x, nu23, jl)
    ! Two-term corrected Olver uniform Airy approximation:
    !
    !   j_l(x) ~= pref * [ Ai(tau)
    !       + eps  * (P1(tau) Ai(tau) + Q1(tau) Ai'(tau))
    !       + eps^2* (P2(tau) Ai(tau) + Q2(tau) Ai'(tau)) ]
    !
    ! nu23 must be (l + 0.5)**(2/3) (passed in to avoid a second power;
    ! the caller already has nu**(1/3)).
    !
    ! Intended domain: the BJL shoulder bands, eta = (x-nu)/nu^(1/3) in
    ! [-2.4,-0.65) U (0.65,3.85] (capped so u = x/nu - 1 is in [-0.26, 0.42])
    ! with l >= 26, i.e. tau = nu^(2/3) zeta in roughly [-4.85, 3.3].
    !
    ! The Olver mapping is evaluated through a single polynomial
    !   P(u) = zeta/u  (analytic through the turning point),
    ! with ratio = 4 zeta/(1-z^2) = -4 P(u)/(2+u), fitted on
    ! u in [-0.26, 0.42] with max |delta zeta| ~ 2e-10. This reproduces the
    ! exact log/acos mapping to ~9e-10 peak-normalized over the gate domain
    ! and is ~1.7x faster. Calls outside the fitted u range are inaccurate;
    ! widen the fit if the eta/u gates in BJL change.
    !
    ! The P1..Q2 corrections are a single weighted least-squares fit of
    ! exact j_l residuals over the gate domain for l >= 26 (so they are
    ! effective coefficients absorbing mapping and higher-order truncation
    ! effects, not the analytic Olver functions); the fit is accurate to
    ! ~3.5e-6 peak-normalized at l = 26, improving towards high l, on top
    ! of the ~1e-6 airy_fast floor (see below).

    implicit none
    integer, intent(in) :: l
    real(dl), intent(in) :: x, nu23
    real(dl), intent(out) :: jl

    real(dl), parameter :: pi = 3.141592653589793238462643383279502884197_dl
    real(dl) :: ax, nu, zeta, tau, eps, pref, ratio
    real(dl) :: u, pu
    real(dl) :: ai, aip
    real(dl) :: p1, q1, p2, q2
    integer :: k

    ! P(u) = zeta/u monomial fit on u in [-0.26, 0.42], max fit error 6.6e-10.
    real(dl), parameter :: pc(0:10) = (/ &
        -1.25992104996821208e+00_dl, +3.77976310917083169e-01_dl, &
        -2.30385514085488935e-01_dl, +1.65910344421374395e-01_dl, &
        -1.29319657869220062e-01_dl, +1.05647827110811374e-01_dl, &
        -8.89230826554308490e-02_dl, +7.73642310050780407e-02_dl, &
        -7.19131136272972149e-02_dl, +6.33363780675322147e-02_dl, &
        -3.21834319307111594e-02_dl /)

    real(dl), parameter :: c(24) = (/ &
        -3.68007293378092709e-02_dl, +1.52097125394687113e-02_dl, &
        +9.45343855707463099e-02_dl, +5.29457464968014282e-02_dl, &
        +1.01673865085949999e-02_dl, +6.32096478967918564e-04_dl, &
        -5.02972020102563366e-02_dl, +5.78573966116145630e-02_dl, &
        +7.95031046204105196e-02_dl, +2.60219848032685157e-02_dl, &
        +3.01851478297665045e-03_dl, +9.36429652248560787e-05_dl, &
        -4.78203638758779526e-01_dl, -8.70173349815381719e-02_dl, &
        +7.22863184145759674e-01_dl, +4.11559433361555371e-01_dl, &
        +7.15565959786843286e-02_dl, +3.57357294204752518e-03_dl, &
        -6.38867125689511428e-01_dl, +3.58715633433204339e-01_dl, &
        +6.31913102077036770e-01_dl, +1.96335539191201919e-01_dl, &
        +1.90209925594464375e-02_dl, +2.98832052141368223e-04_dl /)

    ax = abs(x)

    if (ax == 0.0_dl) then
        if (l == 0) then
            jl = 1.0_dl
        else
            jl = 0.0_dl
        end if
        return
    end if

    nu = real(l, dl) + 0.5_dl
    u = ax/nu - 1.0_dl

    pu = pc(10)
    do k = 9, 0, -1
        pu = pu*u + pc(k)
    end do
    zeta = u*pu
    ratio = -4.0_dl*pu/(2.0_dl + u)

    eps = 1.0_dl / nu23
    tau = nu23 * zeta

    call airy_fast(tau, ai, aip)

    pref = sqrt(pi/(2.0_dl*ax*nu23) * sqrt(ratio))

    ! Horner evaluation of correction polynomials.
    p1 = c(6)
    q1 = c(12)
    p2 = c(18)
    q2 = c(24)

    do k = 4, 0, -1
        p1 = p1*tau + c(1+k)
        q1 = q1*tau + c(7+k)
        p2 = p2*tau + c(13+k)
        q2 = q2*tau + c(19+k)
    end do

    jl = pref * (ai*(1.0_dl + eps*(p1 + eps*p2)) + aip*eps*(q1 + eps*q2))

    if (x < 0.0_dl .and. mod(l, 2) /= 0) jl = -jl

    end subroutine bjl_uniform_airy_fast



    elemental subroutine airy_fast(x, ai, aip)
    ! Fast real Airy Ai(x) and Ai'(x), optimized for < 1e-7 absolute error.
    !
    ! Branches:
    !   x < -6.5               simplified oscillatory asymptotic fallback
    !  -6.5 <= x < -4.4        inline Horner polynomial fit, negative tail
    !  -4.4 <= x < -2.09       inline Horner polynomial fit, negative middle interval
    !  -2.09 <= x <  0.0       inline Horner polynomial fit, central negative interval
    !   0.0 <= x <  2.09       inline Horner polynomial fit, central positive interval
    !   2.09 <= x <= 2.98      inline Horner polynomial fit, positive middle interval
    !   2.98 < x <= 6.5        inline Horner polynomial fit, positive tail
    !   6.5 < x <= 25.77       first-corrected decaying asymptotic fallback
    !   x > 25.77              returns zero for Ai and Ai'
    !
    !
    ! Hot-interval fits are ordinary Horner polynomials in y, where y maps
    ! each subinterval to [-1,1].  Coefficients were generated by Chebyshev-node
    ! interpolation of scipy.special.airy and then converted to power form to
    ! remove Clenshaw recurrence overhead.
    !
    ! Measured on 130001 points in [-6.5,6.5] against scipy.special.airy:
    !   full interval:          max |Ai error|  = 7.52e-08,
    !                           max |Ai' error| = 7.36e-08.
    !   negative tail:     max |Ai error|  = 7.52e-08,
    !                   max |Ai' error| = 6.32e-08.
    !   negative middle:   max |Ai error|  = 1.56e-08,
    !                   max |Ai' error| = 1.14e-08.
    !   central negative:  max |Ai error|  = 5.76e-08,
    !                   max |Ai' error| = 4.86e-09.
    !   central positive:  max |Ai error|  = 3.75e-08,
    !                   max |Ai' error| = 6.21e-08.
    !   positive middle:   max |Ai error|  = 3.14e-08,
    !                   max |Ai' error| = 7.36e-08.
    !   positive tail:     max |Ai error|  = 7.52e-08,
    !                   max |Ai' error| = 6.30e-08.
    !
    ! Relative Airy errors can be large near Airy zeros, so absolute error
    ! is the intended accuracy diagnostic for this routine.
    !
    ! Outside the hot interval the fallbacks are deliberately simplified for
    ! speed at the existing +/-6.5 cut points: the positive fallback keeps one
    ! correction term in the decaying asymptotic form, while the negative
    ! fallback uses first-degree polynomial corrections to the oscillatory
    ! asymptotic form.  They are not intended to be moved much closer to the
    ! origin without refitting.
    !
    ! Asymptotic-branch accuracy, measured against scipy.special.airy on
    ! [-50,-6.5) and (6.5,25.77], has max absolute errors about 2.3e-8/5.9e-8
    ! for Ai/Ai' on the negative side and 7.9e-10/2.4e-9 on the positive side.
    ! The positive decaying branch has worst relative errors about 2.9e-4 for
    ! Ai and 3.3e-4 for Ai' before the x > 25.77 zero cutoff; beyond the
    ! cutoff relative error is 100%, but absolute error is below 1e-38 at
    ! entry.  On the oscillatory negative side, relative error is meaningful
    ! only away from zeros; at sampled oscillation peaks it is below 6.4e-8.
    ! The positive-tail polynomial is endpoint-matched to the first-corrected
    ! asymptotic value at x = 6.5, avoiding a visible branch jump there.
    implicit none

    real(dl), intent(in)  :: x
    real(dl), intent(out) :: ai, aip

    real(dl), parameter :: xneg_tail = -6.5_dl
    real(dl), parameter :: xneg1 = -4.4_dl
    real(dl), parameter :: xneg2 = -2.09_dl
    real(dl), parameter :: xpos1 = 2.09_dl
    real(dl), parameter :: xpos2 = 2.98_dl
    real(dl), parameter :: xpos_tail = 6.5_dl
    real(dl), parameter :: amaxairy = 25.77_dl

    real(dl), parameter :: sqpii = 5.64189583547756286948e-1_dl
    real(dl), parameter :: half_sqpii = 2.82094791773878143474e-1_dl
    real(dl), parameter :: piq = 7.85398163397448309616e-1_dl
    real(dl), parameter :: twothird = 6.66666666666666666667e-1_dl
    real(dl), parameter :: threehalf = 1.5_dl

    real(dl) :: y, a0, p0
    real(dl) :: q, rt, qtr, zeta, z, zz, invq
    real(dl) :: theta, sn, cs, ak, h, t
    real(dl) :: rf, drf, rg, drg
    real(dl) :: uf, ug, duf_dz, dug_dz
    real(dl) :: dzdx, dakdx

    if (x >= 0.0_dl) then
        if (x < xpos1) then
            y = (x - 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 6.51840630617294678e-5_dl
            a0 = a0*y - 1.41820887227489902e-4_dl
            a0 = a0*y - 5.41899097043535300e-4_dl
            a0 = a0*y + 3.84764097766738741e-3_dl
            a0 = a0*y - 8.25368803111877630e-3_dl
            a0 = a0*y - 6.03495061405065480e-3_dl
            a0 = a0*y + 7.31872209255700207e-2_dl
            a0 = a0*y - 1.59974703563899756e-1_dl
            a0 = a0*y + 1.28267372044362643e-1_dl

            p0 = -7.67433498257096941e-5_dl
            p0 = p0*y + 4.94942129462562475e-4_dl
            p0 = p0*y - 8.15070380453042321e-4_dl
            p0 = p0*y - 3.10529357109378757e-3_dl
            p0 = p0*y + 1.83370369686281925e-2_dl
            p0 = p0*y - 3.15955876925410259e-2_dl
            p0 = p0*y - 1.73129913754885845e-2_dl
            p0 = p0*y + 1.40071486522548944e-1_dl
            p0 = p0*y - 1.53086150287526923e-1_dl

            ai  = a0
            aip = p0
            return
        else if (x <= xpos2) then
            y = (x - 2.53500000000000014e+0_dl) * 2.24719101123595477e+0_dl

            a0 = -1.19619719384582610e-6_dl
            a0 = a0*y + 7.30607372423890670e-5_dl
            a0 = a0*y - 7.09505585931523709e-4_dl
            a0 = a0*y + 3.72297934482035409e-3_dl
            a0 = a0*y - 1.10827654040494637e-2_dl
            a0 = a0*y + 1.48308424366850235e-2_dl

            p0 = -1.15848713346688742e-5_dl
            p0 = p0*y - 1.25032074830176718e-5_dl
            p0 = p0*y + 6.68294530383771890e-4_dl
            p0 = p0*y - 4.78388666680492399e-3_dl
            p0 = p0*y + 1.67303273172219359e-2_dl
            p0 = p0*y - 2.49050322117747754e-2_dl

            ai  = a0
            aip = p0
            return
        else if (x <= xpos_tail) then
            y = (x - 4.74000000000000021e+0_dl) * 5.68181818181818232e-1_dl

            a0 = 4.88795369996967941e-6_dl
            a0 = a0*y - 9.70309364545231360e-5_dl
            a0 = a0*y + 3.89113179462119549e-4_dl
            a0 = a0*y - 8.69120696960508647e-4_dl
            a0 = a0*y + 1.39931440824699808e-3_dl
            a0 = a0*y - 1.68525416034986312e-3_dl
            a0 = a0*y + 1.43003674608063172e-3_dl
            a0 = a0*y - 7.63904323725029025e-4_dl
            a0 = a0*y + 1.94752926593294709e-4_dl

            p0 = -2.77619561700030726e-5_dl
            p0 = p0*y + 5.01739366022296743e-5_dl
            p0 = p0*y + 7.70886695526835241e-5_dl
            p0 = p0*y - 4.82613277256875291e-4_dl
            p0 = p0*y + 1.29184781654615517e-3_dl
            p0 = p0*y - 2.41087049189967146e-3_dl
            p0 = p0*y + 3.18769569668455807e-3_dl
            p0 = p0*y - 2.88379934388740405e-3_dl
            p0 = p0*y + 1.62470517339803242e-3_dl
            p0 = p0*y - 4.33700533338538723e-4_dl

            ai  = a0
            aip = p0
            return
        else
            if (x > amaxairy) then
                ai = 0.0_dl
                aip = 0.0_dl
                return
            end if

            ! Positive fallback: x > 6.5 and x <= AMAXAIRY.
            ! Keep the first correction term in z = 1/zeta.  This is still
            ! much cheaper than the full rational fallback, reduces relative
            ! error near the cutoff, and matches the positive-tail polynomial
            ! endpoint at x = 6.5.
            rt = sqrt(x)

            zeta = twothird * x * rt
            qtr = sqrt(rt)
            z = 1.0_dl / zeta

            a0 = half_sqpii * exp(-zeta) / qtr
            ai = a0 * (1.0_dl - 6.94444444444444444e-2_dl*z)
            aip = -a0 * rt * (1.0_dl + 9.72222222222222222e-2_dl*z)
            return
        end if
    else
        if (x >= xneg2) then
            y = (x + 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 5.06580358435135449e-6_dl
            a0 = a0*y - 4.61082668202196131e-4_dl
            a0 = a0*y + 1.46688819390411451e-3_dl
            a0 = a0*y + 2.55013843603054671e-3_dl
            a0 = a0*y - 2.30836607936877845e-2_dl
            a0 = a0*y + 3.05136923162561725e-2_dl
            a0 = a0*y + 9.89646761942541003e-2_dl
            a0 = a0*y - 3.05531198486228894e-1_dl
            a0 = a0*y + 1.51357911442611820e-2_dl
            a0 = a0*y + 5.35467705092311230e-1_dl

            p0 = -5.15798967096810657e-5_dl
            p0 = p0*y + 2.36129299092830948e-4_dl
            p0 = p0*y + 1.49302570673894751e-4_dl
            p0 = p0*y - 4.00216684041114143e-3_dl
            p0 = p0*y + 9.75397014554868838e-3_dl
            p0 = p0*y + 1.49519953817788265e-2_dl
            p0 = p0*y - 1.10429690672384659e-1_dl
            p0 = p0*y + 1.16724977031110524e-1_dl
            p0 = p0*y + 2.84107726562980778e-1_dl
            p0 = p0*y - 5.84744087911635013e-1_dl
            p0 = p0*y + 1.44840201291639853e-2_dl

            ai  = a0
            aip = p0
            return
        else if (x >= xneg1) then
            y = (x + 3.24500000000000011e+0_dl) * 8.65800865800865571e-1_dl

            a0 = -9.54841120863266526e-5_dl
            a0 = a0*y - 2.40218737639588486e-4_dl
            a0 = a0*y + 2.18870768526106905e-3_dl
            a0 = a0*y - 7.26102098338769265e-4_dl
            a0 = a0*y - 2.16614327619476833e-2_dl
            a0 = a0*y + 4.14599225907353111e-2_dl
            a0 = a0*y + 9.39438729129915828e-2_dl
            a0 = a0*y - 3.26534621448200502e-1_dl
            a0 = a0*y - 1.11226136786098134e-1_dl
            a0 = a0*y + 9.06926447670876823e-1_dl
            a0 = a0*y + 5.02411976022500845e-3_dl
            a0 = a0*y - 4.19008537866886632e-1_dl

            p0 = 3.22834884762855353e-4_dl
            p0 = p0*y - 9.08270311970813633e-4_dl
            p0 = p0*y - 2.88893546837624927e-3_dl
            p0 = p0*y + 1.70524592992353889e-2_dl
            p0 = p0*y - 4.29883302410300362e-3_dl
            p0 = p0*y - 1.31279611846237987e-1_dl
            p0 = p0*y + 2.15091069681845409e-1_dl
            p0 = p0*y + 4.06682899587775981e-1_dl
            p0 = p0*y - 1.13081105767151313e+0_dl
            p0 = p0*y - 2.88899003685929689e-1_dl
            p0 = p0*y + 1.57043347576836512e+0_dl
            p0 = p0*y + 4.34988645555066042e-3_dl

            ai  = a0
            aip = p0
            return
        else if (x >= xneg_tail) then
            y = (x + 5.45000000000000018e+0_dl) * 9.52380952380952550e-1_dl

            a0 = -1.07510299195806721e-3_dl
            a0 = a0*y + 1.61081861596598392e-3_dl
            a0 = a0*y + 1.25898167350718959e-2_dl
            a0 = a0*y - 3.17287935514204492e-2_dl
            a0 = a0*y - 6.97215951205423234e-2_dl
            a0 = a0*y + 2.55547774507571435e-1_dl
            a0 = a0*y + 1.77851990438656710e-1_dl
            a0 = a0*y - 8.85682387350856803e-1_dl
            a0 = a0*y - 1.82574851147874584e-1_dl
            a0 = a0*y + 8.96114353477748704e-1_dl
            a0 = a0*y + 6.07711584563111612e-2_dl

            p0 = 7.27402496617779817e-4_dl
            p0 = p0*y + 1.71964876060855081e-4_dl
            p0 = p0*y - 1.19077060916135732e-2_dl
            p0 = p0*y + 1.34154904279652990e-2_dl
            p0 = p0*y + 9.72590787238585691e-2_dl
            p0 = p0*y - 2.11216276860003105e-1_dl
            p0 = p0*y - 3.98848523137400046e-1_dl
            p0 = p0*y + 1.21679578395402688e+0_dl
            p0 = p0*y + 6.77583856837294385e-1_dl
            p0 = p0*y - 2.53051030517908737e+0_dl
            p0 = p0*y - 3.47762941229583233e-1_dl
            p0 = p0*y + 8.53442055237177311e-1_dl

            ai  = a0
            aip = p0
            return
        else
            ! Negative fallback: x < -6.5.
            ! The correction functions are first-degree polynomial fits in
            ! zz = z**2 over the actual fallback domain.  This keeps the
            ! boundary errors matched to the hot interval while avoiding two
            ! high-degree rational evaluations and their derivatives.
            q = -x
            invq = 1.0_dl / q
            rt = sqrt(q)

            zeta = twothird * q * rt
            qtr = sqrt(rt)
            ak = sqpii / qtr

            z = 1.0_dl / zeta
            zz = z*z

            rf = -3.71305224980621948e-2_dl + zz*5.54256415602097374e-2_dl
            drf = 5.54256415602097374e-2_dl
            rg = 6.94432325868511585e-2_dl - zz*3.70901829479361872e-2_dl
            drg = -3.70901829479361872e-2_dl

            uf = 1.0_dl + zz*rf
            duf_dz = (2.0_dl*z) * (rf + zz*drf)

            ug = z*rg
            dug_dz = rg + (2.0_dl*zz)*drg

            theta = zeta + piq
            sn = sin(theta)
            cs = cos(theta)

            h = sn*uf - cs*ug
            ai = ak*h

            dzdx = threehalf*z*invq
            dakdx = 0.25_dl*ak*invq
            t = (-rt)*(cs*uf + sn*ug) + dzdx*(sn*duf_dz - cs*dug_dz)
            aip = dakdx*h + ak*t
            return
        end if
    end if


    end subroutine airy_fast



    ELEMENTAL SUBROUTINE BJL_RECURRENCE(L, X, JL)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: L
    REAL(dl), INTENT(IN) :: X
    REAL(dl), INTENT(OUT) :: JL

    INTEGER :: n, Nstart, margin
    REAL(dl) :: ax, j0, j1, jm1, jcur, jp1, scale

    ! Scalars for Miller downward recurrence
    REAL(dl) :: w_np1, w_n, w_nm1, w0, w1, wL
    LOGICAL :: have_wL

    REAL(dl), PARAMETER :: big = 1.0E200_dl
    REAL(dl), PARAMETER :: small = 1.0E-200_dl

    IF (L < 0) ERROR STOP 'Can not evaluate Spherical Bessel Function with index l<0'

    ax = ABS(X)

    IF (ax == 0.0E0_dl) THEN
        IF (L == 0) THEN
            JL = 1.0E0_dl
        ELSE
            JL = 0.0E0_dl
        END IF
        RETURN
    END IF

    IF (ax < 1.0E-4_dl) THEN
        j0 = 1.0E0_dl - ax**2/6.0E0_dl + ax**4/120.0E0_dl - ax**6/5040.0E0_dl
        j1 = ax/3.0E0_dl * (1.0E0_dl - ax**2/10.0E0_dl + ax**4/280.0E0_dl - ax**6/15120.0E0_dl)
    ELSE
        j0 = SIN(ax)/ax
        j1 = SIN(ax)/ax**2 - COS(ax)/ax
    END IF

    IF (L == 0) THEN
        JL = j0

    ELSE IF (L == 1) THEN
        JL = j1

    ELSE IF (ax > REAL(L, dl)) THEN
        jm1 = j0
        jcur = j1

        DO n = 1, L-1
            jp1 = (REAL(2*n+1, dl)/ax)*jcur - jm1
            jm1 = jcur
            jcur = jp1
        END DO

        JL = jcur

    ELSE
        margin = MAX(80, INT(12.0E0_dl*SQRT(REAL(L+1, dl))))
        Nstart = MAX(L + margin, INT(ax) + margin)

        ! Miller downward recurrence:
        w_np1 = 0.0E0_dl
        w_n   = 1.0E0_dl

        w0 = 0.0E0_dl
        w1 = 0.0E0_dl
        wL = 0.0E0_dl
        have_wL = .FALSE.

        DO n = Nstart, 1, -1
            w_nm1 = (REAL(2*n+1, dl)/ax)*w_n - w_np1

            IF (n-1 == L) THEN
                wL = w_nm1
                have_wL = .TRUE.
            END IF

            IF (ABS(w_nm1) > big) THEN
                w_nm1 = w_nm1 * small
                w_n   = w_n   * small
                w_np1 = w_np1 * small

                IF (have_wL) wL = wL * small
            END IF

            IF (n == 1) THEN
                w0 = w_nm1
                w1 = w_n
            END IF

            w_np1 = w_n
            w_n   = w_nm1
        END DO

        IF (ABS(w0) > ABS(w1)) THEN
            scale = j0 / w0
        ELSE
            scale = j1 / w1
        END IF

        JL = wL * scale
    END IF

    IF (X < 0.0E0_dl .AND. MOD(L, 2) /= 0) JL = -JL

    END SUBROUTINE BJL_RECURRENCE

    end module FlatBessels


    SUBROUTINE BJL_EXTERNAL(L,X,JL)
    use FlatBessels
    use Precision
    IMPLICIT NONE
    INTEGER L
    real(dl) X,JL

    call BJL(L,X,JL)

    END SUBROUTINE BJL_EXTERNAL
