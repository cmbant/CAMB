    !CAMB flat spherical Bessel function routines
    !May 2026: updated bjl, accurate to peak-normalized fraction 1e-5 at L>50, max ~4e-5 at BJL_RECURRENCE_MAX_L+1
    !          (pre-splining). Spline table accurate to 2e-4 in tail, better round peak.

    module FlatBessels
    use Precision
    use results
    use RangeUtils
    use MpiUtils
    use splines
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
    real(dl) file_acc

    type(TRanges), save:: BessRanges

    public bessel_horner, BessRanges, InitSpherBessels, bjl_pre_peak_start_factor
    public bjl, Bessels_Free

    contains

    subroutine InitSpherBessels(lSamp, CP,max_bessels_l_index, max_bessels_etak)
    Type(lSamples) lSamp
    Type(CAMBParams) :: CP
    integer, intent(in) :: max_bessels_l_index
    real(dl), intent(in) :: max_bessels_etak
    integer :: new_kmaxfile, old_num_xx
    logical :: same_lsamp_and_acc

    new_kmaxfile = int(max_bessels_etak)+1

    !See if already loaded with enough (and correct) lSamp%l values and k*eta values
    same_lsamp_and_acc = .false.
    if (allocated(bessel_horner) .and. lSamp%nl <= file_l%nl) then
        if (all(file_l%l(1:lSamp%nl)==lSamp%l(1:lSamp%nl)) .and. max_bessels_l_index <= max_ix .and. &
            abs(CP%Accuracy%BesselBoost*CP%Accuracy%AccuracyBoost - file_acc) < 1d-2) then
            same_lsamp_and_acc = .true.
        end if
    end if

    if (same_lsamp_and_acc) then
        if (new_kmaxfile <= kmaxfile) return
        ! Extend existing table to cover the larger x range
        old_num_xx = num_xx
        kmaxfile = new_kmaxfile
        call ExtendBessels(lSamp, CP, old_num_xx)
        if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Extended Bessels'
        return
    end if

    !Haven't made them before, so make them now
    kmaxfile = new_kmaxfile
    max_ix = min(max_bessels_l_index,lSamp%nl)

    call GenerateBessels(lSamp, CP)

    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Calculated Bessels'

    end subroutine InitSpherBessels

    subroutine GenerateBessels(lSamp, CP)
    Type(lSamples) lSamp
    Type(CAMBParams) :: CP
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
    if (do_bispectrum) kmaxfile = kmaxfile*2

    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'x_max bessels', kmaxfile

    call BessRanges%Init()

    call BessRanges%Add_delta(0._dl, 1._dl,0.01_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(1._dl, 5._dl,0.1_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(5._dl, 25._dl,0.2_dl/CP%Accuracy%BesselBoost)
    file_acc = CP%Accuracy%BesselBoost*CP%Accuracy%AccuracyBoost
    call BessRanges%Add_delta(25._dl, 150._dl,0.5_dl/file_acc)
    call BessRanges%Add_delta(150._dl, real(kmaxfile,dl),0.7_dl/file_acc) !2e-4 accuracy

    call BessRanges%GetArray(.false.)
    num_xx = BessRanges%npoints


    if (allocated(bessel_horner)) then
        if (any(ubound(bessel_horner) < [4, num_xx - 1, max_ix])) deallocate(bessel_horner)
    end if
    if (.not. allocated(bessel_horner)) then
        allocate(bessel_horner(1:4,1:num_xx-1,1:max_ix))
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC)
    do j=1,max_ix
        block
            real(dl) :: h2over6, y0, y1, d0, d1, xlim
            real(dl) :: knot_vals(num_xx), spline_y2(num_xx)
            integer :: min_ix, i

            xlim = max(lSamp%l(j) - bjl_pre_peak_start_factor*lSamp%l(j)**(1._dl/3._dl) - 1, &
                cut(min(cut_max_l, lSamp%l(j))))
            min_ix = max(1, BessRanges%IndexOf(xlim) - 1)
            knot_vals(1:max(min_ix-1,1)) = 0
            do  i=min_ix,num_xx
                call bjl(lSamp%l(j),BessRanges%points(i),knot_vals(i))
            end do

            call spline_def(BessRanges%points,knot_vals,num_xx,spline_y2)
            bessel_horner(:,1:max(min_ix-1,1),j) = 0
            do i=max(1,min_ix-1),num_xx-1
                y0 = knot_vals(i)
                y1 = knot_vals(i+1)
                d0 = spline_y2(i)
                d1 = spline_y2(i+1)
                h2over6 = (BessRanges%points(i+1) - BessRanges%points(i))**2/6
                bessel_horner(1,i,j) = y1
                bessel_horner(2,i,j) = y0 - y1 - h2over6*(d0 + 2*d1)
                bessel_horner(3,i,j) = 3*h2over6*d1
                bessel_horner(4,i,j) = y0 - y1 - bessel_horner(2,i,j) - bessel_horner(3,i,j)
            end do
        end block
    end do
    !$OMP END PARALLEL DO

    end subroutine GenerateBessels

    subroutine ExtendBessels(lSamp, CP, old_num_xx)
    ! Extend bessel_horner to the new (larger) kmaxfile without recomputing the
    ! already-stored intervals.  The x-grid points in [0, old_kmaxfile] are
    ! identical in old and new BessRanges (same step 0.7/file_acc in the sparse
    ! segment [150, kmaxfile]), so we only need to compute bjl for the junction
    ! overlap and the new tail.
    Type(lSamples) lSamp
    Type(CAMBParams) :: CP
    integer, intent(in) :: old_num_xx
    integer :: j, ext_start_ix, n_ext
    integer, parameter :: JUNCTION_N = 10  ! overlap points kept from old grid

    if (do_bispectrum) kmaxfile = kmaxfile*2

    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'Extending flat Bessels to x_max', kmaxfile

    call BessRanges%Init()
    call BessRanges%Add_delta(0._dl,  1._dl,   0.01_dl/CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(1._dl,  5._dl,   0.1_dl /CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(5._dl,  25._dl,  0.2_dl /CP%Accuracy%BesselBoost)
    call BessRanges%Add_delta(25._dl, 150._dl, 0.5_dl /file_acc)
    call BessRanges%Add_delta(150._dl, real(kmaxfile,dl), 0.7_dl/file_acc)
    call BessRanges%GetArray(.false.)
    num_xx = BessRanges%npoints

    ext_start_ix = old_num_xx - JUNCTION_N   ! first point of overlap region (1-based)
    n_ext = num_xx - ext_start_ix + 1        ! points: ext_start_ix .. num_xx

    ! Reallocate bessel_horner, keeping all old intervals intact
    block
        real(dl), allocatable :: old_horner(:,:,:)
        call move_alloc(bessel_horner, old_horner)
        allocate(bessel_horner(1:4,1:num_xx-1,1:max_ix))
        bessel_horner(:, 1:old_num_xx-1, :) = old_horner
    end block

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC)
    do j=1,max_ix
        block
            real(dl) :: ext_x(n_ext), ext_y(n_ext), ext_y2(n_ext)
            real(dl) :: h2over6, y0, y1, d0, d1
            integer :: i, store_ix

            do i=1,n_ext
                ext_x(i) = BessRanges%points(ext_start_ix + i - 1)
            end do
            do i=1,n_ext
                call bjl(lSamp%l(j), ext_x(i), ext_y(i))
            end do

            call spline_def(ext_x, ext_y, n_ext, ext_y2)

            ! Store Horner only for strictly new intervals (from old_num_xx onward);
            ! the JUNCTION_N overlap intervals serve only to anchor the spline.
            do i=JUNCTION_N+1,n_ext-1
                store_ix = ext_start_ix + i - 1   ! = old_num_xx when i = JUNCTION_N+1
                y0 = ext_y(i)
                y1 = ext_y(i+1)
                d0 = ext_y2(i)
                d1 = ext_y2(i+1)
                h2over6 = (ext_x(i+1) - ext_x(i))**2/6
                bessel_horner(1,store_ix,j) = y1
                bessel_horner(2,store_ix,j) = y0 - y1 - h2over6*(d0 + 2*d1)
                bessel_horner(3,store_ix,j) = 3*h2over6*d1
                bessel_horner(4,store_ix,j) = y0 - y1 - bessel_horner(2,store_ix,j) - bessel_horner(3,store_ix,j)
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
!   in the transition bands, using fast approximation elsewhere
!   where they are accurate.
!
! accurate to peak-normalized fraction 1e-5 at L>50, max ~4e-5 at BJL_RECURRENCE_MAX_L+1.

    SUBROUTINE BJL(L, X, JL)
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
    REAL(dl) :: TRIGARG, EXPTERM, L3, ETA

    real(dl), parameter :: BJL_RECURRENCE_ETA_LOW  = -5.0_dl
    real(dl), parameter :: BJL_RECURRENCE_ETA_HIGH =  5.6_dl
    real(dl), parameter :: BJL_AIRY_ETA_LOW  = -2.25_dl
    real(dl), parameter :: BJL_PEAK_ETA_LOW  = -0.65_dl
    real(dl), parameter :: BJL_PEAK_ETA_HIGH =  0.65_dl
    real(dl), parameter :: BJL_AIRY_ETA_HIGH =  3.7_dl
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
    ELSE IF ((REAL(L, dl)**2/AX) < 5.0E-1_dl) THEN
        ! Far past the peak: oscillatory large-x expansion.
        BETA = AX - PID2*REAL(L + 1, dl)
        JL = (COS(BETA)*(1.0_dl - (NU2 - 0.25_dl)*(NU2 - 2.25_dl)/8.0_dl/AX2 &
            *(1.0_dl - (NU2 - 6.25_dl)*(NU2 - 12.25_dl)/48.0_dl/AX2)) &
            - SIN(BETA)*(NU2 - 0.25_dl)/2.0_dl/AX &
            *(1.0_dl - (NU2 - 2.25_dl)*(NU2 - 6.25_dl)/24.0_dl/AX2 &
            *(1.0_dl - (NU2 - 12.25_dl)*(NU2 - 20.25_dl)/80.0_dl/AX2)))/AX
    ELSE
        ! Turning-point neighbourhood: classify by normalized distance from peak.
        L3 = NU**0.325_dl
        ETA = (AX - NU)/L3

        ! Moderate orders only need recurrence in the broad transition band.
        ! Deep pre-peak and far-x cases have already returned via cheaper branches.
        IF (L <= BJL_RECURRENCE_MAX_L .AND. &
            ETA >= BJL_RECURRENCE_ETA_LOW .AND. ETA <= BJL_RECURRENCE_ETA_HIGH) THEN
            CALL BJL_RECURRENCE(L, X, JL)
            RETURN
        END IF

        IF ((ETA >= BJL_AIRY_ETA_LOW .AND. ETA < BJL_PEAK_ETA_LOW) .OR. &
            (ETA > BJL_PEAK_ETA_HIGH .AND. ETA <= BJL_AIRY_ETA_HIGH)) THEN
            ! Shoulder bands: corrected uniform Airy approximation.
            CALL BJL_UNIFORM_AIRY_FAST(L, X, JL)
            RETURN
        ELSE IF (ETA < BJL_AIRY_ETA_LOW) THEN
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
        ELSE IF (ETA > BJL_AIRY_ETA_HIGH) THEN
            ! Above the peak: oscillatory post-peak asymptotic form.
            COSB = NU/AX
            SX = SQRT(AX2 - NU2)
            COTB = NU/SX
            SECB = AX/NU
            BETA = ACOS(COSB)
            COT3B = COTB**3
            COT6B = COT3B**2
            SEC2B = SECB**2
            TRIGARG = NU/COTB - NU*BETA - PID4 &
                - ((2.0_dl + 3.0_dl*SEC2B)*COT3B/24.0_dl &
                + (16.0_dl - (1512.0_dl + (3654.0_dl + 375.0_dl*SEC2B)*SEC2B)*SEC2B) &
                *COT3B*COT6B/5760.0_dl/NU2)/NU
            EXPTERM = ((4.0_dl + SEC2B)*SEC2B*COT6B/16.0_dl &
                - (32.0_dl + (288.0_dl + (232.0_dl + 13.0_dl*SEC2B)*SEC2B)*SEC2B) &
                *SEC2B*COT6B**2/128.0_dl/NU2)/NU2
            JL = SQRT(COTB*COSB)/NU*EXP(-EXPTERM)*COS(TRIGARG)
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




    subroutine bjl_uniform_airy_fast(l, x, jl)
    ! Two-term corrected Olver uniform Airy approximation:
    !
    !   j_l(x) ~= pref * [ Ai(tau)
    !       + eps  * (P1(tau) Ai(tau) + Q1(tau) Ai'(tau))
    !       + eps^2* (P2(tau) Ai(tau) + Q2(tau) Ai'(tau)) ]
    !
    ! Intended domain: turning-point region with
    !   tau = nu^(2/3) zeta in roughly [-6.5, 6.5].
    ! Dense-grid validation against scipy.special.spherical_jn, using
    ! l = 1,2,5,10,20,50,100,200,500,1000,2000,5000 and 2601 tau values:
    !
    !   full tau [-6.5,6.5], max absolute error:
    !     l>=10: 2.77e-5, l>=20: 4.75e-6,
    !     l>=50: 3.08e-7, l>=100: 2.78e-8, l>=200: 1.04e-8.
    !
    !   first positive maximum of j_l, found from exact derivative zero;
    !   maximum relative error over tested l>=lmin:
    !     l>=10: 2.28e-5, l>=20: 5.72e-6,
    !     l>=50: 6.76e-7, l>=100: 7.98e-8, l>=200: 6.35e-9.

    implicit none
    integer, intent(in) :: l
    real(dl), intent(in) :: x
    real(dl), intent(out) :: jl

    real(dl), parameter :: pi = 3.141592653589793238462643383279502884197_dl
    real(dl) :: ax, nu, nu13, nu23, z, zeta, tau, eps, eps2, pref, ratio
    real(dl) :: t, s, val, denom, u
    real(dl) :: ai, aip
    real(dl) :: p1, q1, p2, q2
    integer :: k

    real(dl), parameter :: c(24) = (/ &
        -1.57687029733452188e-04_dl, -2.86291841602745543e-04_dl, &
        -1.64967917049690195e-04_dl, -4.20144950724723239e-05_dl, &
        -4.91866087021318015e-06_dl, -2.16160522018663616e-07_dl, &
        -2.19528732498613553e-04_dl, -2.36964997600504554e-04_dl, &
        -9.00604510654162075e-05_dl, -1.59170174801830321e-05_dl, &
        -1.26798459743804658e-06_dl, -3.48372195594890281e-08_dl, &
        6.66207928274303168e-03_dl,  1.23061280051791379e-02_dl, &
        7.08775888397744549e-03_dl,  1.80375527531944258e-03_dl, &
        2.10966849276599796e-04_dl,  9.26086631302043278e-06_dl, &
        2.74389063894026689e-02_dl,  1.04229584215560795e-02_dl, &
        3.86950969205413983e-03_dl,  6.82819070773826560e-04_dl, &
        5.42682642617796167e-05_dl,  1.48342002616123936e-06_dl /)

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
    z = ax / nu

    ! Stable series through the turning point z=1.
    ! zeta = -2^(1/3)*(z-1) + O((z-1)^2)
    ! ratio = 4*zeta/(1-z^2), whose limit is 2^(4/3).
    if (abs(z - 1.0_dl) < 1.0e-3_dl) then
        u = z - 1.0_dl

        zeta = ((((-1.2931387086451008907e-1_dl * u &
            +1.6590960364964869484e-1_dl) * u &
            -2.3038556340934823584e-1_dl) * u &
            +3.7797631496846194943e-1_dl) * u &
            -1.2599210498948731648e0_dl) * u

        ratio = (((( 7.9171433706231762385e-1_dl * u &
            -1.0661731906665948914e0_dl) * u &
            +1.4687079667345950035e0_dl) * u &
            -2.0158736798317970636e0_dl) * u &
            +2.5198420997897463295e0_dl)

    else if (z < 1.0_dl) then
        t = sqrt(max(0.0_dl, 1.0_dl - z*z))
        val = log((1.0_dl + t)/z) - t
        zeta = (1.5_dl * val)**(2.0_dl/3.0_dl)
        denom = 1.0_dl - z*z
        ratio = 4.0_dl * zeta / denom

    else
        s = sqrt(max(0.0_dl, z*z - 1.0_dl))
        val = s - acos(1.0_dl/z)
        zeta = - (1.5_dl * val)**(2.0_dl/3.0_dl)
        denom = 1.0_dl - z*z
        ratio = 4.0_dl * zeta / denom
    end if

    nu13 = nu**(1.0_dl/3.0_dl)
    nu23 = nu13 * nu13
    eps = 1.0_dl / nu23
    eps2 = eps * eps

    tau = nu23 * zeta

    call airy_fast(tau, ai, aip)

    pref = sqrt(pi/(2.0_dl*ax)) * sqrt(sqrt(ratio)) / nu13

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

    jl = pref * (ai + eps*(p1*ai + q1*aip) + eps2*(p2*ai + q2*aip))

    if (x < 0.0_dl .and. mod(l, 2) /= 0) jl = -jl

    end subroutine bjl_uniform_airy_fast


    subroutine airy_fast(x, ai, aip)
    ! Fast real Airy Ai(x) and Ai'(x).
    !
    ! Branches:
    !   x < -4.4               Cephes-style oscillatory asymptotic fallback
    !  -4.4 <= x < -2.09       Chebyshev fit
    !  -2.09 <= x <  2.09      shortened Taylor/Horner series
    !   2.09 <= x <= 2.98      Chebyshev fit
    !   2.98 < x <= 25.77      Cephes-style decaying asymptotic fallback
    !   x > 25.77              returns zero for Ai and Ai'
    !
    ! Measured on 130001 points in [-6.5,6.5] against scipy.special.airy:
    !   full interval:        max |Ai error|  = 5.46e-7,
    !                         max |Ai' error| = 2.08e-6.
    !   central branch:       max |Ai error|  = 5.05e-8,
    !                         max |Ai' error| = 5.11e-7.
    !   negative Cheb branch: max |Ai error|  = 5.46e-7,
    !                         max |Ai' error| = 2.08e-6.
    !   positive Cheb branch: max |Ai error|  = 9.67e-8,
    !                         max |Ai' error| = 7.98e-7.
    !
    ! Relative Airy errors can be large near Airy zeros; Bessel accuracy is
    ! better characterized by the j_l first-maximum relative errors above.
    implicit none

    real(dl), intent(in)  :: x
    real(dl), intent(out) :: ai, aip

    real(dl), parameter :: amaxairy = 25.77_dl
    real(dl), parameter :: pi = 3.141592653589793238462643383279502884197_dl
    real(dl), parameter :: c1 = 0.35502805388781723926_dl
    real(dl), parameter :: c2 = 0.258819403792806798405_dl
    real(dl), parameter :: sqpii = 5.64189583547756286948e-1_dl

    real(dl), parameter :: an(8) = (/ &
        3.46538101525629032477e-1_dl, &
        1.20075952739645805542e1_dl, &
        7.62796053615234516538e1_dl, &
        1.68089224934630576269e2_dl, &
        1.59756391350164413639e2_dl, &
        7.05360906840444183113e1_dl, &
        1.40264691163389668864e1_dl, &
        9.99999999999999995305e-1_dl /)

    real(dl), parameter :: ad(8) = (/ &
        5.67594532638770212846e-1_dl, &
        1.47562562584847203173e1_dl, &
        8.45138970141474626562e1_dl, &
        1.77318088145400459522e2_dl, &
        1.64234692871529701831e2_dl, &
        7.14778400825575695274e1_dl, &
        1.40959135607834029598e1_dl, &
        1.00000000000000000470e0_dl /)

    real(dl), parameter :: afn(9) = (/ &
        -1.31696323418331795333e-1_dl, &
        -6.26456544431912369773e-1_dl, &
        -6.93158036036933542233e-1_dl, &
        -2.79779981545119124951e-1_dl, &
        -4.91900132609500318020e-2_dl, &
        -4.06265923594885404393e-3_dl, &
        -1.59276496239262096340e-4_dl, &
        -2.77649108155232920844e-6_dl, &
        -1.67787698489114633780e-8_dl /)

    real(dl), parameter :: afd(9) = (/ &
        1.33560420706553243746e1_dl, &
        3.26825032795224613948e1_dl, &
        2.67367040941499554804e1_dl, &
        9.18707402907259625840e0_dl, &
        1.47529146771666414581e0_dl, &
        1.15687173795188044134e-1_dl, &
        4.40291641615211203805e-3_dl, &
        7.54720348287414296618e-5_dl, &
        4.51850092970580378464e-7_dl /)

    real(dl), parameter :: agn(11) = (/ &
        1.97339932091685679179e-2_dl, &
        3.91103029615688277255e-1_dl, &
        1.06579897599595591108e0_dl, &
        9.39169229816650230044e-1_dl, &
        3.51465656105547619242e-1_dl, &
        6.33888919628925490927e-2_dl, &
        5.85804113048388458567e-3_dl, &
        2.82851600836737019778e-4_dl, &
        6.98793669997260967291e-6_dl, &
        8.11789239554389293311e-8_dl, &
        3.41551784765923618484e-10_dl /)

    real(dl), parameter :: agd(10) = (/ &
        9.30892908077441974853e0_dl, &
        1.98352928718312140417e1_dl, &
        1.55646628932864612953e1_dl, &
        5.47686069422975497931e0_dl, &
        9.54293611618961883998e-1_dl, &
        8.64580826352392193095e-2_dl, &
        4.12656523824222607191e-3_dl, &
        1.01259085116509135510e-4_dl, &
        1.17166733214413521882e-6_dl, &
        4.91834570062930015649e-9_dl /)

    real(dl) :: q, rt, qtr, zeta, z, zz
    real(dl) :: theta, sn, cs, ak, h
    real(dl) :: pn, pd, dpn, dpd
    real(dl) :: r, dr, rf, drf, rg, drg
    real(dl) :: uf, ug, duf_dz, dug_dz
    real(dl) :: dzdx, dthdx, dakdx
    real(dl) :: f, g, df, dg, z3

    if (x > amaxairy) then
        ai = 0.0_dl
        aip = 0.0_dl
        return
    end if

    if (x >= -4.4_dl .and. x < -2.09_dl) then
        call airy_neg_cheb_fast(x, ai, aip)
        return
    end if

    if (x >= 2.09_dl .and. x <= 2.98_dl) then
        call airy_pos_cheb_fast(x, ai, aip)
        return
    end if

    ! Negative fallback: x < -4.4.
    if (x < -4.4_dl) then
        q = -x
        rt = sqrt(q)

        zeta = 2.0_dl * q * rt / 3.0_dl
        qtr = sqrt(rt)
        ak = sqpii / qtr

        z = 1.0_dl / zeta
        zz = z * z

        call polevl_der_fast(zz, afn, 8, pn, dpn)
        call p1evl_der_fast(zz, afd, 9, pd, dpd)
        rf = pn / pd
        drf = (dpn*pd - pn*dpd) / (pd*pd)

        call polevl_der_fast(zz, agn, 10, pn, dpn)
        call p1evl_der_fast(zz, agd, 10, pd, dpd)
        rg = pn / pd
        drg = (dpn*pd - pn*dpd) / (pd*pd)

        uf = 1.0_dl + zz * rf
        duf_dz = 2.0_dl * z * (rf + zz*drf)

        ug = z * rg
        dug_dz = rg + 2.0_dl * zz * drg

        theta = zeta + 0.25_dl * pi
        sn = sin(theta)
        cs = cos(theta)

        h = sn*uf - cs*ug
        ai = ak * h

        dzdx = 1.5_dl * z / q
        dthdx = -rt
        dakdx = ak / (4.0_dl * q)

        aip = dakdx*h + ak * ( &
            dthdx*(cs*uf + sn*ug) &
            + sn*duf_dz*dzdx &
            - cs*dug_dz*dzdx )

        return
    end if

    ! Positive fallback: x > 2.98 and x <= AMAXAIRY.
    if (x > 2.98_dl) then
        rt = sqrt(x)

        zeta = 2.0_dl * x * rt / 3.0_dl
        qtr = sqrt(rt)
        z = 1.0_dl / zeta

        call polevl_der_fast(z, an, 7, pn, dpn)
        call polevl_der_fast(z, ad, 7, pd, dpd)

        r = pn / pd
        dr = (dpn*pd - pn*dpd) / (pd*pd)

        ai = 0.5_dl * sqpii * r * exp(-zeta) / qtr

        dzdx = -1.5_dl * z / x
        aip = ai * ((dr/r)*dzdx - 0.25_dl/x - rt)

        return
    end if

    ! Central branch: -2.09 <= x < 2.09.
    z3 = x*x*x

    f = 9.09662613850461810e-12_dl
    f = 2.78356759838241336e-09_dl + z3*f
    f = 5.84549195660306779e-07_dl + z3*f
    f = 7.71604938271604923e-05_dl + z3*f
    f = 5.55555555555555577e-03_dl + z3*f
    f = 1.66666666666666657e-01_dl + z3*f
    f = 1.0_dl + z3*f

    g = 1.72172984605299955e-12_dl
    g = 5.88831607350125831e-10_dl + z3*g
    g = 1.41319585764030204e-07_dl + z3*g
    g = 2.20458553791887140e-05_dl + z3*g
    g = 1.98412698412698402e-03_dl + z3*g
    g = 8.33333333333333287e-02_dl + z3*g
    g = x * (1.0_dl + z3*g)

    df = 1.63739270493083139e-10_dl
    df = 4.17535139757362004e-08_dl + z3*df
    df = 7.01459034792368135e-06_dl + z3*df
    df = 6.94444444444444471e-04_dl + z3*df
    df = 3.33333333333333329e-02_dl + z3*df
    df = 5.00000000000000000e-01_dl + z3*df
    df = x*x * df

    dg = 3.27128670750069899e-11_dl
    dg = 9.42130571760201330e-09_dl + z3*dg
    dg = 1.83715461493239276e-06_dl + z3*dg
    dg = 2.20458553791887113e-04_dl + z3*dg
    dg = 1.38888888888888881e-02_dl + z3*dg
    dg = 3.33333333333333315e-01_dl + z3*dg
    dg = 1.0_dl + z3*dg

    ai  = c1*f  - c2*g
    aip = c1*df - c2*dg

    CONTAINS


    subroutine airy_neg_cheb_fast(xv, aiv, aipv)
    implicit none
    real(dl), intent(in)  :: xv
    real(dl), intent(out) :: aiv, aipv

    integer, parameter :: nch = 10
    real(dl), parameter :: xmid  = -3.24500000000000011e+00_dl
    real(dl), parameter :: xhalf =  1.15500000000000025e+00_dl

    real(dl) :: y, twoy
    real(dl) :: a0, a1, a2
    real(dl) :: p0, p1, p2
    integer :: j

    real(dl), parameter :: cai(10) = (/ &
        -7.52972311373991954e-02_dl, -3.04924834824625880e-02_dl, &
        3.09214054941112593e-01_dl, -4.86933307406160927e-03_dl, &
        -3.32582282959081738e-02_dl,  3.79467424343539342e-03_dl, &
        1.22912822483253726e-03_dl, -2.66641673714729916e-04_dl, &
        -1.03644439619337352e-05_dl,  7.52394509644803996e-06_dl /)

    real(dl), parameter :: caip(10) = (/ &
        -2.41791445513994847e-02_dl,  8.53129932850260952e-01_dl, &
        4.44254809627120731e-03_dl, -2.17741253092986226e-01_dl, &
        2.97377848446421338e-02_dl,  1.26187697358148727e-02_dl, &
        -3.11653760934400706e-03_dl, -1.51393639068021042e-04_dl, &
        1.15482678107576387e-04_dl, -7.81692617584523665e-06_dl /)

    y = (xv - xmid) / xhalf
    twoy = 2.0_dl * y

    a1 = 0.0_dl
    a2 = 0.0_dl
    p1 = 0.0_dl
    p2 = 0.0_dl

    do j = nch, 2, -1
        a0 = twoy*a1 - a2 + cai(j)
        p0 = twoy*p1 - p2 + caip(j)

        a2 = a1
        a1 = a0

        p2 = p1
        p1 = p0
    end do

    aiv  = y*a1 - a2 + cai(1)
    aipv = y*p1 - p2 + caip(1)

    end subroutine airy_neg_cheb_fast


    subroutine airy_pos_cheb_fast(xv, aiv, aipv)
    implicit none
    real(dl), intent(in)  :: xv
    real(dl), intent(out) :: aiv, aipv

    integer, parameter :: nch = 5
    real(dl), parameter :: xmid  = 2.53500000000000014e+00_dl
    real(dl), parameter :: xhalf = 4.45000000000000062e-01_dl

    real(dl) :: y, twoy
    real(dl) :: a0, a1, a2
    real(dl) :: p0, p1, p2
    integer :: j

    real(dl), parameter :: cai(5) = (/ &
        1.67197298855612173e-02_dl, -1.16156422167520770e-02_dl, &
        1.89802004119299199e-03_dl, -1.77750208417504322e-04_dl, &
        9.13251766892205164e-06_dl /)

    real(dl), parameter :: caip(5) = (/ &
        -2.73016642479830610e-02_dl,  1.72243076704324746e-02_dl, &
        -2.39819493752986600e-03_dl,  1.63453367574147994e-04_dl, &
        -1.56291392306149089e-06_dl /)

    y = (xv - xmid) / xhalf
    twoy = 2.0_dl * y

    a1 = 0.0_dl
    a2 = 0.0_dl
    p1 = 0.0_dl
    p2 = 0.0_dl

    do j = nch, 2, -1
        a0 = twoy*a1 - a2 + cai(j)
        p0 = twoy*p1 - p2 + caip(j)

        a2 = a1
        a1 = a0

        p2 = p1
        p1 = p0
    end do

    aiv  = y*a1 - a2 + cai(1)
    aipv = y*p1 - p2 + caip(1)

    end subroutine airy_pos_cheb_fast


    subroutine polevl_der_fast(xp, coef, n, p, dp)
    implicit none
    integer, intent(in) :: n
    real(dl), intent(in) :: xp
    real(dl), intent(in) :: coef(n+1)
    real(dl), intent(out) :: p, dp

    integer :: i

    p = coef(1)
    dp = 0.0_dl

    do i = 2, n+1
        dp = dp*xp + p
        p = p*xp + coef(i)
    end do

    end subroutine polevl_der_fast


    subroutine p1evl_der_fast(xp, coef, n, p, dp)
    implicit none
    integer, intent(in) :: n
    real(dl), intent(in) :: xp
    real(dl), intent(in) :: coef(n)
    real(dl), intent(out) :: p, dp

    integer :: i

    p = 1.0_dl
    dp = 0.0_dl

    do i = 1, n
        dp = dp*xp + p
        p = p*xp + coef(i)
    end do

    end subroutine p1evl_der_fast

    end subroutine airy_fast





    SUBROUTINE BJL_RECURRENCE(L, X, JL)
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
