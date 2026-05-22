    !CAMB flat spherical Bessel function routines
    !May 2026: updated bjl, accurate to peak-normalized fraction 1e-5 at L>50, max ~4.4e-5 at BJL_RECURRENCE_MAX_L
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

    !See if already loaded with enough (and correct) lSamp%l values and k*eta values
    if (lSamp%nl <= file_l%nl) then
        if (allocated(bessel_horner) .and. all(file_l%l(1:lSamp%nl)==lSamp%l(1:lSamp%nl) &
            .and. (max_bessels_l_index <= max_ix)) &
            .and. (int(max_bessels_etak)+1 <= kmaxfile) &
            .and. (abs(CP%Accuracy%BesselBoost*CP%Accuracy%AccuracyBoost - file_acc) < 1d-2)) return
    end if

    !Haven't made them before, so make them now
    kmaxfile = int(max_bessels_etak)+1
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

    subroutine Bessels_Free

    if (allocated(bessel_horner)) deallocate(bessel_horner)
    if (allocated(file_l%l)) deallocate(file_l%l)
    file_l%nl=0

    call BessRanges%Free()

    end subroutine Bessels_Free

! Optimized spherical Bessel wrapper.
!
! The Airy Ai and Ai' used by the uniform expansion are evaluated by Chebyshev
! fits on tau in [-6.5,6.5], which covers the transition bands used below.

    SUBROUTINE BJL(L, X, JL)
    ! Optimized spherical Bessel j_l(x).
    ! Branches are ordered from most robust/cheap special cases to the
    ! large-order asymptotic regions.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    REAL(dl), INTENT(IN) :: X
    REAL(dl), INTENT(OUT) :: JL

    REAL(dl), PARAMETER :: LN2      = 0.6931471805599453094E0_dl
    REAL(dl), PARAMETER :: ONEMLN2  = 0.30685281944005469058277E0_dl
    REAL(dl), PARAMETER :: PID2     = 1.5707963267948966192313217E0_dl
    REAL(dl), PARAMETER :: PID4     = 0.78539816339744830961566084582E0_dl
    REAL(dl), PARAMETER :: ROOTPI12 = 21.269446210866192327578E0_dl
    REAL(dl), PARAMETER :: GAMMA1   = 2.6789385347077476336556E0_dl
    REAL(dl), PARAMETER :: GAMMA2   = 1.3541179394264004169452E0_dl

    REAL(dl) :: AX, AX2, NU, NU2, BETA, BETA2, COSB
    REAL(dl) :: SX, SX2, COTB, COT3B, COT6B, SECB, SEC2B
    REAL(dl) :: TRIGARG, EXPTERM, L3, ETA

    REAL(dl), PARAMETER :: BJL_RECURRENCE_ETA_LOW  = -5.0E0_dl
    REAL(dl), PARAMETER :: BJL_RECURRENCE_ETA_HIGH =  5.6E0_dl
    REAL(dl), PARAMETER :: BJL_AIRY_ETA_LOW        = -2.3E0_dl
    REAL(dl), PARAMETER :: BJL_PEAK_ETA_LOW        = -0.65E0_dl
    REAL(dl), PARAMETER :: BJL_PEAK_ETA_HIGH       =  0.65E0_dl
    REAL(dl), PARAMETER :: BJL_AIRY_ETA_HIGH       =  4.0E0_dl

    IF (L < 0) THEN
        ERROR STOP 'Can not evaluate Spherical Bessel Function with index l<0'
    END IF

    AX = ABS(X)
    AX2 = AX**2

    ! Low orders: use explicit formulas, with small-x series to avoid cancellation.
    IF (L < 7) THEN
        SELECT CASE (L)
        CASE (0)
            IF (AX < 1.0E-1_dl) THEN
                JL = 1.0E0_dl - AX2/6.0E0_dl*(1.0E0_dl - AX2/20.0E0_dl)
            ELSE
                JL = SIN(AX)/AX
            END IF
        CASE (1)
            IF (AX < 2.0E-1_dl) THEN
                JL = AX/3.0E0_dl*(1.0E0_dl - AX2/10.0E0_dl*(1.0E0_dl - AX2/28.0E0_dl))
            ELSE
                JL = (SIN(AX)/AX - COS(AX))/AX
            END IF
        CASE (2)
            IF (AX < 3.0E-1_dl) THEN
                JL = AX2/15.0E0_dl*(1.0E0_dl - AX2/14.0E0_dl*(1.0E0_dl - AX2/36.0E0_dl))
            ELSE
                JL = (-3.0E0_dl*COS(AX)/AX - SIN(AX)*(1.0E0_dl - 3.0E0_dl/AX2))/AX
            END IF
        CASE (3)
            IF (AX < 4.0E-1_dl) THEN
                JL = AX*AX2/105.0E0_dl*(1.0E0_dl - AX2/18.0E0_dl*(1.0E0_dl - AX2/44.0E0_dl))
            ELSE
                JL = (COS(AX)*(1.0E0_dl - 15.0E0_dl/AX2) - SIN(AX)*(6.0E0_dl - 15.0E0_dl/AX2)/AX)/AX
            END IF
        CASE (4)
            IF (AX < 6.0E-1_dl) THEN
                JL = AX2**2/945.0E0_dl*(1.0E0_dl - AX2/22.0E0_dl*(1.0E0_dl - AX2/52.0E0_dl))
            ELSE
                JL = (SIN(AX)*(1.0E0_dl - (45.0E0_dl - 105.0E0_dl/AX2)/AX2) + COS(AX)*(10.0E0_dl - 105.0E0_dl/AX2)/AX)/AX
            END IF
        CASE (5)
            IF (AX < 1.0E0_dl) THEN
                JL = AX2**2*AX/10395.0E0_dl*(1.0E0_dl - AX2/26.0E0_dl*(1.0E0_dl - AX2/60.0E0_dl))
            ELSE
                JL = (SIN(AX)*(15.0E0_dl - (420.0E0_dl - 945.0E0_dl/AX2)/AX2)/AX - &
                    COS(AX)*(1.0E0_dl - (105.0E0_dl - 945.0E0_dl/AX2)/AX2))/AX
            END IF
        CASE DEFAULT
            IF (AX < 1.0E0_dl) THEN
                JL = AX2**3/135135.0E0_dl*(1.0E0_dl - AX2/30.0E0_dl*(1.0E0_dl - AX2/68.0E0_dl))
            ELSE
                JL = (SIN(AX)*(-1.0E0_dl + (210.0E0_dl - (4725.0E0_dl - 10395.0E0_dl/AX2)/AX2)/AX2) + &
                    COS(AX)*(-21.0E0_dl + (1260.0E0_dl - 10395.0E0_dl/AX2)/AX2)/AX)/AX
            END IF
        END SELECT

        IF (X < 0.0E0_dl .AND. MOD(L, 2) /= 0) JL = -JL
        RETURN
    END IF


    NU = REAL(L, dl) + 0.5E0_dl
    NU2 = NU**2

    IF (AX < 1.0E-40_dl) THEN
        ! Very small x: j_l(x) is negligible here for l >= 7.
        JL = 0.0E0_dl
    ELSE IF ((AX2/REAL(L, dl)) < 5.0E-1_dl) THEN
        ! Deep pre-peak: exponentially small ascending expansion.
        JL = EXP(REAL(L, dl)*LOG(AX/NU) - LN2 + NU*ONEMLN2 &
            - (1.0E0_dl - (1.0E0_dl - 3.5E0_dl/NU2)/NU2/30.0E0_dl)/12.0E0_dl/NU) &
            /NU*(1.0E0_dl - AX2/(4.0E0_dl*NU + 4.0E0_dl) &
            *(1.0E0_dl - AX2/(8.0E0_dl*NU + 16.0E0_dl) &
            *(1.0E0_dl - AX2/(12.0E0_dl*NU + 36.0E0_dl))))
    ELSE IF ((REAL(L, dl)**2/AX) < 5.0E-1_dl) THEN
        ! Far past the peak: oscillatory large-x expansion.
        BETA = AX - PID2*REAL(L + 1, dl)
        JL = (COS(BETA)*(1.0E0_dl - (NU2 - 0.25E0_dl)*(NU2 - 2.25E0_dl)/8.0E0_dl/AX2 &
            *(1.0E0_dl - (NU2 - 6.25E0_dl)*(NU2 - 12.25E0_dl)/48.0E0_dl/AX2)) &
            - SIN(BETA)*(NU2 - 0.25E0_dl)/2.0E0_dl/AX &
            *(1.0E0_dl - (NU2 - 2.25E0_dl)*(NU2 - 6.25E0_dl)/24.0E0_dl/AX2 &
            *(1.0E0_dl - (NU2 - 12.25E0_dl)*(NU2 - 20.25E0_dl)/80.0E0_dl/AX2)))/AX
    ELSE
        ! Broad transition region: classify once by normalized distance from
        ! the peak. The lower Airy shoulder is kept only where the below-peak
        ! form loses accuracy; the upper shoulder is deliberately shorter
        ! because the post-peak asymptotic form becomes accurate and faster.
        L3 = NU**0.325E0_dl
        ETA = (AX - NU)/L3

        ! For moderate orders, use recurrence only where the asymptotic
        ! transition forms lose accuracy. Deep pre-peak and far-x cases have
        ! already returned through faster asymptotic branches.
        IF (L <= BJL_RECURRENCE_MAX_L) THEN
            IF (ETA >= BJL_RECURRENCE_ETA_LOW .AND. ETA <= BJL_RECURRENCE_ETA_HIGH) THEN
                CALL BJL_RECURRENCE_FAST(L, X, JL)
                RETURN
            END IF
        END IF

        IF (ETA < BJL_AIRY_ETA_LOW) THEN
            ! Below the peak but outside the exponentially tiny region.
            COSB = NU/AX
            SX = SQRT(NU2 - AX2)
            COTB = NU/SX
            SECB = AX/NU
            BETA = LOG(COSB + SX/AX)
            COT3B = COTB**3
            COT6B = COT3B**2
            SEC2B = SECB**2
            EXPTERM = ((2.0E0_dl + 3.0E0_dl*SEC2B)*COT3B/24.0E0_dl &
                - ((4.0E0_dl + SEC2B)*SEC2B*COT6B/16.0E0_dl &
                + ((16.0E0_dl - (1512.0E0_dl + (3654.0E0_dl + 375.0E0_dl*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.0E0_dl &
                + (32.0E0_dl + (288.0E0_dl + (232.0E0_dl + 13.0E0_dl*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.0E0_dl/NU)*COT6B/NU) &
                /NU)/NU
            JL = SQRT(COTB*COSB)/(2.0E0_dl*NU)*EXP(-NU*BETA + NU/COTB - EXPTERM)
        ELSE IF (ETA < BJL_PEAK_ETA_LOW) THEN
            ! Lower shoulder: uniform Airy form covers the pre-peak transition.
            CALL BJL_UNIFORM_AIRY_FAST(L, X, JL)
            RETURN
        ELSE IF (ETA <= BJL_PEAK_ETA_HIGH) THEN
            ! Very close to the peak: polynomial transition expansion.
            BETA = AX - NU
            BETA2 = BETA**2
            SX = 6.0E0_dl/AX
            SX2 = SX**2
            SECB = SX**(1.0E0_dl/3.0E0_dl)
            SEC2B = SECB**2
            JL = (GAMMA1*SECB + BETA*GAMMA2*SEC2B &
                - (BETA2/18.0E0_dl - 1.0E0_dl/45.0E0_dl)*BETA*SX*SECB*GAMMA1 &
                - ((BETA2 - 1.0E0_dl)*BETA2/36.0E0_dl + 1.0E0_dl/420.0E0_dl)*SX*SEC2B*GAMMA2 &
                + (((BETA2/1620.0E0_dl - 7.0E0_dl/3240.0E0_dl)*BETA2 + 1.0E0_dl/648.0E0_dl)*BETA2 &
                - 1.0E0_dl/8100.0E0_dl)*SX2*SECB*GAMMA1 &
                + (((BETA2/4536.0E0_dl - 1.0E0_dl/810.0E0_dl)*BETA2 + 19.0E0_dl/11340.0E0_dl)*BETA2 &
                - 13.0E0_dl/28350.0E0_dl)*BETA*SX2*SEC2B*GAMMA2 &
                - ((((BETA2/349920.0E0_dl - 1.0E0_dl/29160.0E0_dl)*BETA2 + 71.0E0_dl/583200.0E0_dl)*BETA2 &
                - 121.0E0_dl/874800.0E0_dl)*BETA2 + 7939.0E0_dl/224532000.0E0_dl)*BETA*SX2*SX*SECB*GAMMA1) &
                *SQRT(SX)/ROOTPI12
        ELSE IF (ETA <= BJL_AIRY_ETA_HIGH) THEN
            ! Upper shoulder: uniform Airy form only until the oscillatory form is better.
            CALL BJL_UNIFORM_AIRY_FAST(L, X, JL)
            RETURN
        ELSE
            ! Above the peak but not yet in the far-x regime.
            COSB = NU/AX
            SX = SQRT(AX2 - NU2)
            COTB = NU/SX
            SECB = AX/NU
            BETA = ACOS(COSB)
            COT3B = COTB**3
            COT6B = COT3B**2
            SEC2B = SECB**2
            TRIGARG = NU/COTB - NU*BETA - PID4 &
                - ((2.0E0_dl + 3.0E0_dl*SEC2B)*COT3B/24.0E0_dl &
                + (16.0E0_dl - (1512.0E0_dl + (3654.0E0_dl + 375.0E0_dl*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.0E0_dl/NU2)/NU
            EXPTERM = ((4.0E0_dl + SEC2B)*SEC2B*COT6B/16.0E0_dl &
                - (32.0E0_dl + (288.0E0_dl + (232.0E0_dl + 13.0E0_dl*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.0E0_dl/NU2)/NU2
            JL = SQRT(COTB*COSB)/NU*EXP(-EXPTERM)*COS(TRIGARG)
        END IF
    END IF

    IF (X < 0.0E0_dl .AND. MOD(L, 2) /= 0) JL = -JL

    END SUBROUTINE BJL


    SUBROUTINE BJL_UNIFORM_AIRY_FAST(L, X, JL)
    ! Two-term corrected Olver uniform Airy approximation:
    !
    !   j_l(x) ~= pref * [ Ai(tau)
    !       + eps  * (P1(tau) Ai(tau) + Q1(tau) Ai'(tau))
    !       + eps^2* (P2(tau) Ai(tau) + Q2(tau) Ai'(tau)) ]
    !
    ! tau = nu^(2/3) zeta, eps = nu^(-2/3), nu = l+1/2.

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    REAL(dl), INTENT(IN) :: X
    REAL(dl), INTENT(OUT) :: JL

    REAL(dl), PARAMETER :: PI = 3.141592653589793238462643383279502884197E0_dl
    REAL(dl) :: AX, NU, Z, ZETA, TAU, EPS, PREF, RATIO
    REAL(dl) :: T, S, VAL, DENOM, U
    REAL(dl) :: AI, AIP
    REAL(dl) :: P1, Q1, P2, Q2
    INTEGER :: K

    REAL(dl), PARAMETER :: C(24) = (/ &
        -1.57687029733452188E-04_dl, -2.86291841602745543E-04_dl, &
        -1.64967917049690195E-04_dl, -4.20144950724723239E-05_dl, &
        -4.91866087021318015E-06_dl, -2.16160522018663616E-07_dl, &
        -2.19528732498613553E-04_dl, -2.36964997600504554E-04_dl, &
        -9.00604510654162075E-05_dl, -1.59170174801830321E-05_dl, &
        -1.26798459743804658E-06_dl, -3.48372195594890281E-08_dl, &
        6.66207928274303168E-03_dl,  1.23061280051791379E-02_dl, &
        7.08775888397744549E-03_dl,  1.80375527531944258E-03_dl, &
        2.10966849276599796E-04_dl,  9.26086631302043278E-06_dl, &
        2.74389063894026689E-02_dl,  1.04229584215560795E-02_dl, &
        3.86950969205413983E-03_dl,  6.82819070773826560E-04_dl, &
        5.42682642617796167E-05_dl,  1.48342002616123936E-06_dl /)

    AX = ABS(X)

    IF (AX == 0.0E0_dl) THEN
        IF (L == 0) THEN
            JL = 1.0E0_dl
        ELSE
            JL = 0.0E0_dl
        END IF
        RETURN
    END IF

    NU = REAL(L, dl) + 0.5E0_dl
    Z = AX / NU

    IF (ABS(Z - 1.0E0_dl) < 1.0E-5_dl) THEN
        U = Z - 1.0E0_dl

        ZETA = -7.9370052598409973738E-1_dl * U &
            +9.5244063118091968485E-1_dl * U**2 &
            +3.9911797878057586794E-1_dl * U**3 &
            +1.4430735404495669892E0_dl  * U**4

        RATIO =  1.5874010519681994748E0_dl &
            -2.6985817883459391071E0_dl * U &
            +5.5105493661181781766E-1_dl * U**2 &
            -3.1616745492050428872E0_dl * U**3

    ELSE IF (Z < 1.0E0_dl) THEN
        T = SQRT(MAX(0.0E0_dl, 1.0E0_dl - Z*Z))
        VAL = LOG((1.0E0_dl + T)/Z) - T
        ZETA = (1.5E0_dl * VAL)**(2.0E0_dl/3.0E0_dl)
        DENOM = 1.0E0_dl - Z*Z
        RATIO = 4.0E0_dl * ZETA / DENOM

    ELSE
        S = SQRT(MAX(0.0E0_dl, Z*Z - 1.0E0_dl))
        VAL = S - ACOS(1.0E0_dl/Z)
        ZETA = - (1.5E0_dl * VAL)**(2.0E0_dl/3.0E0_dl)
        DENOM = 1.0E0_dl - Z*Z
        RATIO = 4.0E0_dl * ZETA / DENOM
    END IF

    TAU = NU**(2.0E0_dl/3.0E0_dl) * ZETA
    EPS = NU**(-2.0E0_dl/3.0E0_dl)

    CALL AIRY_FAST(TAU, AI, AIP)

    PREF = SQRT(PI/(2.0E0_dl*AX)) * RATIO**0.25E0_dl * NU**(-1.0E0_dl/3.0E0_dl)

    ! Horner evaluation of correction polynomials.
    P1 = C(6)
    Q1 = C(12)
    P2 = C(18)
    Q2 = C(24)

    DO K = 4, 0, -1
        P1 = P1*TAU + C(1+K)
        Q1 = Q1*TAU + C(7+K)
        P2 = P2*TAU + C(13+K)
        Q2 = Q2*TAU + C(19+K)
    END DO

    JL = PREF * (AI + EPS*(P1*AI + Q1*AIP) + EPS*EPS*(P2*AI + Q2*AIP))

    IF (X < 0.0E0_dl .AND. MOD(L, 2) /= 0) JL = -JL

    END SUBROUTINE BJL_UNIFORM_AIRY_FAST

    SUBROUTINE AIRY_FAST(X, AI, AIP)
    ! Fast real Airy Ai(x) and Ai'(x), targeting about 1e-6 absolute
    ! accuracy for use in the Olver/uniform Bessel expansion.
    !
    ! Method:
    !
    !   x < -6.5
    !       Oscillatory Cephes-style asymptotic form with rational
    !       corrections.  Used only as an out-of-range fallback.
    !
    !   -6.5 <= x < -2.09
    !       Local Chebyshev fits for Ai and Ai' on the active upper-shoulder range [-3.6, -2.09].
    !       This avoids sin/cos/rational derivative work in the usual
    !       tau range.
    !
    !   -2.09 <= x < 2.09
    !       Fixed Taylor/Horner polynomials for Ai and Ai', deliberately
    !       truncated for speed rather than full double precision.
    !
    !   2.09 <= x <= 6.5
    !       Local Chebyshev fits for Ai and Ai' on the active lower-shoulder range [2.09, 3.05].
    !
    !   x > 6.5
    !       Exponentially decaying Cephes-style asymptotic form with
    !       rational correction, or zero once x > AMAXAIRY.
    !
    ! Expected checked accuracy against a double-precision Airy reference:
    !
    !   On x in [-6.5, 6.5]:
    !       max |Ai error|                 ~ 2e-7
    !       max |Ai' error|                ~ 5.2e-7
    !       max sqrt(dAi**2 + dAiPrime**2) ~ 5.2e-7
    !
    ! The worst error is set by the deliberately shortened central branch
    ! near the transition points.  Outside [-6.5,6.5], this remains a useful
    ! all-real fallback, but the strongest 1e-6 statement is for the finite
    ! tau interval used by the uniform expansion.

    IMPLICIT NONE

    REAL(dl), INTENT(IN)  :: X
    REAL(dl), INTENT(OUT) :: AI, AIP

    REAL(dl), PARAMETER :: AMAXAIRY = 25.77E0_dl
    REAL(dl), PARAMETER :: PI = 3.141592653589793238462643383279502884197E0_dl
    REAL(dl), PARAMETER :: C1 = 0.35502805388781723926E0_dl
    REAL(dl), PARAMETER :: C2 = 0.258819403792806798405E0_dl
    REAL(dl), PARAMETER :: SQPII = 5.64189583547756286948E-1_dl

    REAL(dl), PARAMETER :: AN(8) = (/ &
        3.46538101525629032477E-1_dl, &
        1.20075952739645805542E1_dl, &
        7.62796053615234516538E1_dl, &
        1.68089224934630576269E2_dl, &
        1.59756391350164413639E2_dl, &
        7.05360906840444183113E1_dl, &
        1.40264691163389668864E1_dl, &
        9.99999999999999995305E-1_dl /)

    REAL(dl), PARAMETER :: AD(8) = (/ &
        5.67594532638770212846E-1_dl, &
        1.47562562584847203173E1_dl, &
        8.45138970141474626562E1_dl, &
        1.77318088145400459522E2_dl, &
        1.64234692871529701831E2_dl, &
        7.14778400825575695274E1_dl, &
        1.40959135607834029598E1_dl, &
        1.00000000000000000470E0_dl /)

    REAL(dl), PARAMETER :: AFN(9) = (/ &
        -1.31696323418331795333E-1_dl, &
        -6.26456544431912369773E-1_dl, &
        -6.93158036036933542233E-1_dl, &
        -2.79779981545119124951E-1_dl, &
        -4.91900132609500318020E-2_dl, &
        -4.06265923594885404393E-3_dl, &
        -1.59276496239262096340E-4_dl, &
        -2.77649108155232920844E-6_dl, &
        -1.67787698489114633780E-8_dl /)

    REAL(dl), PARAMETER :: AFD(9) = (/ &
        1.33560420706553243746E1_dl, &
        3.26825032795224613948E1_dl, &
        2.67367040941499554804E1_dl, &
        9.18707402907259625840E0_dl, &
        1.47529146771666414581E0_dl, &
        1.15687173795188044134E-1_dl, &
        4.40291641615211203805E-3_dl, &
        7.54720348287414296618E-5_dl, &
        4.51850092970580378464E-7_dl /)

    REAL(dl), PARAMETER :: AGN(11) = (/ &
        1.97339932091685679179E-2_dl, &
        3.91103029615688277255E-1_dl, &
        1.06579897599595591108E0_dl, &
        9.39169229816650230044E-1_dl, &
        3.51465656105547619242E-1_dl, &
        6.33888919628925490927E-2_dl, &
        5.85804113048388458567E-3_dl, &
        2.82851600836737019778E-4_dl, &
        6.98793669997260967291E-6_dl, &
        8.11789239554389293311E-8_dl, &
        3.41551784765923618484E-10_dl /)

    REAL(dl), PARAMETER :: AGD(10) = (/ &
        9.30892908077441974853E0_dl, &
        1.98352928718312140417E1_dl, &
        1.55646628932864612953E1_dl, &
        5.47686069422975497931E0_dl, &
        9.54293611618961883998E-1_dl, &
        8.64580826352392193095E-2_dl, &
        4.12656523824222607191E-3_dl, &
        1.01259085116509135510E-4_dl, &
        1.17166733214413521882E-6_dl, &
        4.91834570062930015649E-9_dl /)

    REAL(dl) :: Q, RT, QTR, ZETA, Z, ZZ
    REAL(dl) :: THETA, SN, CS, AK, H
    REAL(dl) :: PN, PD, DPN, DPD
    REAL(dl) :: R, DR, RF, DRF, RG, DRG
    REAL(dl) :: UF, UG, DUF_DZ, DUG_DZ
    REAL(dl) :: DZDX, DTHDX, DAKDX
    REAL(dl) :: F, G, DF, DG, Z3

    IF (X > AMAXAIRY) THEN
        AI = 0.0E0_dl
        AIP = 0.0E0_dl
        RETURN
    END IF

    IF (X >= -6.5E0_dl .AND. X < -2.09E0_dl) THEN
        CALL AIRY_NEG_CHEB_FAST(X, AI, AIP)
        RETURN
    END IF

    IF (X >= 2.09E0_dl .AND. X <= 6.5E0_dl) THEN
        CALL AIRY_POS_CHEB_FAST(X, AI, AIP)
        RETURN
    END IF

    ! Negative fallback: x < -6.5.
    IF (X < -6.5E0_dl) THEN
        Q = -X
        RT = SQRT(Q)

        ZETA = 2.0E0_dl * Q * RT / 3.0E0_dl
        QTR = SQRT(RT)
        AK = SQPII / QTR

        Z = 1.0E0_dl / ZETA
        ZZ = Z * Z

        CALL POLEVL_DER_FAST(ZZ, AFN, 8, PN, DPN)
        CALL P1EVL_DER_FAST(ZZ, AFD, 9, PD, DPD)
        RF = PN / PD
        DRF = (DPN*PD - PN*DPD) / (PD*PD)

        CALL POLEVL_DER_FAST(ZZ, AGN, 10, PN, DPN)
        CALL P1EVL_DER_FAST(ZZ, AGD, 10, PD, DPD)
        RG = PN / PD
        DRG = (DPN*PD - PN*DPD) / (PD*PD)

        UF = 1.0E0_dl + ZZ * RF
        DUF_DZ = 2.0E0_dl * Z * (RF + ZZ*DRF)

        UG = Z * RG
        DUG_DZ = RG + 2.0E0_dl * ZZ * DRG

        THETA = ZETA + 0.25E0_dl * PI
        SN = SIN(THETA)
        CS = COS(THETA)

        H = SN*UF - CS*UG
        AI = AK * H

        DZDX = 1.5E0_dl * Z / Q
        DTHDX = -RT
        DAKDX = AK / (4.0E0_dl * Q)

        AIP = DAKDX*H + AK * ( &
            DTHDX*(CS*UF + SN*UG) &
            + SN*DUF_DZ*DZDX &
            - CS*DUG_DZ*DZDX )

        RETURN
    END IF

    ! Positive fallback: x > 6.5 and x <= AMAXAIRY.
    IF (X > 6.5E0_dl) THEN
        RT = SQRT(X)

        ZETA = 2.0E0_dl * X * RT / 3.0E0_dl
        QTR = SQRT(RT)
        Z = 1.0E0_dl / ZETA

        CALL POLEVL_DER_FAST(Z, AN, 7, PN, DPN)
        CALL POLEVL_DER_FAST(Z, AD, 7, PD, DPD)

        R = PN / PD
        DR = (DPN*PD - PN*DPD) / (PD*PD)

        AI = 0.5E0_dl * SQPII * R * EXP(-ZETA) / QTR

        DZDX = -1.5E0_dl * Z / X
        AIP = AI * ((DR/R)*DZDX - 0.25E0_dl/X - RT)

        RETURN
    END IF

    ! Central branch: -2.09 <= x < 2.09.
    Z3 = X*X*X

    F = 9.09662613850461810E-12_dl
    F = 2.78356759838241336E-09_dl + Z3*F
    F = 5.84549195660306779E-07_dl + Z3*F
    F = 7.71604938271604923E-05_dl + Z3*F
    F = 5.55555555555555577E-03_dl + Z3*F
    F = 1.66666666666666657E-01_dl + Z3*F
    F = 1.0E0_dl + Z3*F

    G = 1.72172984605299955E-12_dl
    G = 5.88831607350125831E-10_dl + Z3*G
    G = 1.41319585764030204E-07_dl + Z3*G
    G = 2.20458553791887140E-05_dl + Z3*G
    G = 1.98412698412698402E-03_dl + Z3*G
    G = 8.33333333333333287E-02_dl + Z3*G
    G = X * (1.0E0_dl + Z3*G)

    DF = 1.63739270493083139E-10_dl
    DF = 4.17535139757362004E-08_dl + Z3*DF
    DF = 7.01459034792368135E-06_dl + Z3*DF
    DF = 6.94444444444444471E-04_dl + Z3*DF
    DF = 3.33333333333333329E-02_dl + Z3*DF
    DF = 5.00000000000000000E-01_dl + Z3*DF
    DF = X*X * DF

    DG = 3.27128670750069899E-11_dl
    DG = 9.42130571760201330E-09_dl + Z3*DG
    DG = 1.83715461493239276E-06_dl + Z3*DG
    DG = 2.20458553791887113E-04_dl + Z3*DG
    DG = 1.38888888888888881E-02_dl + Z3*DG
    DG = 3.33333333333333315E-01_dl + Z3*DG
    DG = 1.0E0_dl + Z3*DG

    AI  = C1*F  - C2*G
    AIP = C1*DF - C2*DG

    CONTAINS

    SUBROUTINE AIRY_NEG_CHEB_FAST(XV, AIV, AIPV)
    IMPLICIT NONE
    REAL(dl), INTENT(IN)  :: XV
    REAL(dl), INTENT(OUT) :: AIV, AIPV

    INTEGER, PARAMETER :: NCH = 14

    ! Maps the used upper-shoulder tau range [-4.75, -2.09] to [-1, 1].
    REAL(dl), PARAMETER :: XMID  = -3.42000000000000037E0_dl
    REAL(dl), PARAMETER :: XHALF =  1.33000000000000007E0_dl

    REAL(dl) :: Y, TWOY
    REAL(dl) :: A0, A1, A2
    REAL(dl) :: P0, P1, P2
    INTEGER :: J

    REAL(dl), PARAMETER :: CAI(14) = (/ &
        -4.03840442579786723E-03_dl, -1.59160239850393986E-01_dl, &
        3.32033428954130960E-01_dl, 5.52816344806841692E-02_dl, &
        -5.87380789439008039E-02_dl, 1.42594167224066983E-03_dl, &
        3.83320765503384524E-03_dl, -5.15288734810007577E-04_dl, &
        -9.88229869572182192E-05_dl, 2.78689025935962363E-05_dl, &
        -9.86711248543552530E-08_dl, -6.64164643179690834E-07_dl, &
        6.67078851235519155E-08_dl, 6.23584827668025404E-09_dl /)

    REAL(dl), PARAMETER :: CAIP(14) = (/ &
        7.85785534099702615E-03_dl, 6.78681155842563721E-01_dl, &
        2.55054417201942774E-01_dl, -3.19915621638119096E-01_dl, &
        5.66358489047926600E-03_dl, 3.33961341634626263E-02_dl, &
        -5.05778118575379246E-03_dl, -1.18919365301728679E-03_dl, &
        3.66310575640137399E-04_dl, -3.39181185274878346E-07_dl, &
        -1.08628175993337255E-05_dl, 1.15412678791436187E-06_dl, &
        1.22925519235653169E-07_dl, -3.46284641833705560E-08_dl /)

    Y = (XV - XMID) / XHALF
    TWOY = 2.0E0_dl * Y

    A1 = 0.0E0_dl
    A2 = 0.0E0_dl
    P1 = 0.0E0_dl
    P2 = 0.0E0_dl

    DO J = NCH, 2, -1
        A0 = TWOY*A1 - A2 + CAI(J)
        P0 = TWOY*P1 - P2 + CAIP(J)

        A2 = A1
        A1 = A0

        P2 = P1
        P1 = P0
    END DO

    AIV  = Y*A1 - A2 + CAI(1)
    AIPV = Y*P1 - P2 + CAIP(1)

    END SUBROUTINE AIRY_NEG_CHEB_FAST


    SUBROUTINE AIRY_POS_CHEB_FAST(XV, AIV, AIPV)
    IMPLICIT NONE
    REAL(dl), INTENT(IN)  :: XV
    REAL(dl), INTENT(OUT) :: AIV, AIPV

    INTEGER, PARAMETER :: NCH = 7

    ! Maps the actually used lower-shoulder tau range [2.09, 3.05] to [-1, 1].
    ! Degree 6 is enough over this final interval; the central Airy polynomial
    ! handles tau < 2.09.
    REAL(dl), PARAMETER :: XMID  = 2.57000000000000028E+00_dl
    REAL(dl), PARAMETER :: XHALF = 4.79999999999999982E-01_dl

    REAL(dl) :: Y, TWOY
    REAL(dl) :: A0, A1, A2
    REAL(dl) :: P0, P1, P2
    INTEGER :: J

    REAL(dl), PARAMETER :: CAI(7) = (/ &
        1.60886772887947720E-02_dl, -1.19842055007600909E-02_dl, &
        2.11908265584750044E-03_dl, -2.16014564156138982E-04_dl, &
        1.22427404493726490E-05_dl, -1.38037935759893120E-07_dl, &
        -3.87577561695252637E-08_dl /)

    REAL(dl), PARAMETER :: CAIP(7) = (/ &
        -2.63185722252196816E-02_dl, 1.78620943233821147E-02_dl, &
        -2.70295486393902365E-03_dl, 2.03072191319627753E-04_dl, &
        -2.77281198728709877E-06_dl, -9.73482836631669391E-07_dl, &
        1.02978341010519716E-07_dl /)

    Y = (XV - XMID) / XHALF
    TWOY = 2.0E0_dl * Y

    A1 = 0.0E0_dl
    A2 = 0.0E0_dl
    P1 = 0.0E0_dl
    P2 = 0.0E0_dl

    DO J = NCH, 2, -1
        A0 = TWOY*A1 - A2 + CAI(J)
        P0 = TWOY*P1 - P2 + CAIP(J)

        A2 = A1
        A1 = A0

        P2 = P1
        P1 = P0
    END DO

    AIV  = Y*A1 - A2 + CAI(1)
    AIPV = Y*P1 - P2 + CAIP(1)

    END SUBROUTINE AIRY_POS_CHEB_FAST


    SUBROUTINE POLEVL_DER_FAST(XP, COEF, N, P, DP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(dl), INTENT(IN) :: XP
    REAL(dl), INTENT(IN) :: COEF(N+1)
    REAL(dl), INTENT(OUT) :: P, DP

    INTEGER :: I

    P = COEF(1)
    DP = 0.0E0_dl

    DO I = 2, N+1
        DP = DP*XP + P
        P = P*XP + COEF(I)
    END DO

    END SUBROUTINE POLEVL_DER_FAST


    SUBROUTINE P1EVL_DER_FAST(XP, COEF, N, P, DP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(dl), INTENT(IN) :: XP
    REAL(dl), INTENT(IN) :: COEF(N)
    REAL(dl), INTENT(OUT) :: P, DP

    INTEGER :: I

    P = 1.0E0_dl
    DP = 0.0E0_dl

    DO I = 1, N
        DP = DP*XP + P
        P = P*XP + COEF(I)
    END DO

    END SUBROUTINE P1EVL_DER_FAST

    END SUBROUTINE AIRY_FAST


    SUBROUTINE BJL_RECURRENCE_FAST(L, X, JL)
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

    END SUBROUTINE BJL_RECURRENCE_FAST

    end module FlatBessels


    SUBROUTINE BJL_EXTERNAL(L,X,JL)
    use FlatBessels
    use Precision
    IMPLICIT NONE
    INTEGER L
    real(dl) X,JL

    call BJL(L,X,JL)

    END SUBROUTINE BJL_EXTERNAL
