    !CAMB spherical and hyperspherical Bessel function routines
    !This version May 2006 - minor changes to bjl (http://cosmocoffee.info/viewtopic.php?t=530)
    !Feb 2007: fixed for high l, uses Ranges
    !Feb 2009: minor fix for non-flat compiled with non-smart IF evaluation
    !Dec 2011: minor tweak to DoRecurs for smoother errors across flat for L~O(30)
    !May 2026: updated bjl; direct peak-normalized checks now give a ~4.4e-5 floor
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !Flat bessel function module

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
    real(dl) x, xlim
    integer i,j

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
    call BessRanges%Add_delta(150._dl, real(kmaxfile,dl),0.8_dl/file_acc)

    call BessRanges%GetArray(.false.)
    num_xx = BessRanges%npoints


    if (allocated(bessel_horner)) then
        if (any(ubound(bessel_horner) < [4, num_xx - 1, max_ix])) deallocate(bessel_horner)
    end if
    if (.not. allocated(bessel_horner)) then
        allocate(bessel_horner(1:4,1:num_xx-1,1:max_ix))
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(i, x, xlim)
    do j=1,max_ix
        block
            real(dl) :: h2over6, y0, y1, d0, d1
            real(dl) :: knot_vals(num_xx), spline_y2(num_xx)

            do  i=1,num_xx
                x=BessRanges%points(i)
                xlim = lSamp%l(j) - bjl_pre_peak_start_factor*lSamp%l(j)**(1._dl/3._dl)
                if (x > xlim) then
                    if ((lSamp%l(j)==3).and.(x <=0.2) .or. (lSamp%l(j) > 3).and.(x < 0.5) .or. &
                        (lSamp%l(j)>5).and.(x < 1.0)) then
                        knot_vals(i)=0
                    else
                        call bjl(lSamp%l(j),x,knot_vals(i))
                    end if
                else
                    knot_vals(i)=0
                end if
            end do

            call spline_def(BessRanges%points,knot_vals,num_xx,spline_y2)
            do i=1,num_xx-1
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
! Strategy:
!   For low L use (v accurate and still fast) recursive result
!   Elsewhere use a two-term corrected uniform Airy asymptotic
!   in the transition bands, using fast approximate bjl_approx
!   elsewhere where it is accurate.
!
! The Airy Ai and Ai' used by the uniform expansion are evaluated by Chebyshev
! fits on tau in [-6.5,6.5], which covers the transition bands used below.
! Max peak-normalized error is about 4.4e-5 in the current scan. Typical error < 2e-6.

    SUBROUTINE BJL(L, X, JL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    REAL(dl), INTENT(IN) :: X
    REAL(dl), INTENT(OUT) :: JL

    REAL(dl) :: AX, NU, L3, LOGELL
    REAL(dl) :: SMALL_BOUNDARY
    REAL(dl) :: PATCH_LOW, PATCH_HIGH, KEEP
    REAL(dl) :: PATCH_LOWER, KEEP_LOWER, KEEP_UPPER, PATCH_UPPER

    IF (L < 7) THEN
        CALL BJL_approx(L, X, JL) !acutually exact..
        RETURN
    ELSE IF (L <= 25) THEN
        CALL BJL_RECURRENCE_FAST(L, X, JL)
        RETURN
    END IF

    AX = ABS(X)
    NU = REAL(L, dl) + 0.5E0_dl
    L3 = NU**0.325E0_dl
    ! l-dependent switch scaling fitted to the original error envelope.
    LOGELL = LOG(REAL(L, dl))

    PATCH_LOW  = 3.25E0_dl + 0.195E0_dl * LOGELL
    PATCH_HIGH = 5.75E0_dl - 0.075E0_dl * LOGELL
    KEEP       = 1.02E0_dl - 0.065E0_dl * LOGELL

    PATCH_LOW  = MIN(MAX(PATCH_LOW,  4.00E0_dl), 5.00E0_dl)
    PATCH_HIGH = MIN(MAX(PATCH_HIGH, 5.10E0_dl), 5.60E0_dl)
    KEEP       = MIN(MAX(KEEP,       0.46E0_dl), 0.76E0_dl)

    PATCH_LOWER = NU - PATCH_LOW  * L3
    KEEP_LOWER  = NU - KEEP       * L3
    KEEP_UPPER  = NU + KEEP       * L3
    PATCH_UPPER = NU + PATCH_HIGH * L3

    IF ((AX >= PATCH_LOWER .AND. AX < KEEP_LOWER) .OR. &
        (AX > KEEP_UPPER .AND. AX <= PATCH_UPPER)) THEN
        CALL BJL_UNIFORM_AIRY_FAST(L, X, JL)
    ELSE
        CALL BJL_approx(L, X, JL)
    END IF

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
    !       Local Chebyshev fits for Ai and Ai' on [-6.5, -2.09].
    !       This avoids sin/cos/rational derivative work in the usual
    !       tau range.
    !
    !   -2.09 <= x < 2.09
    !       Fixed Taylor/Horner polynomials for Ai and Ai', deliberately
    !       truncated for speed rather than full double precision.
    !
    !   2.09 <= x <= 6.5
    !       Local Chebyshev fits for Ai and Ai' on [2.09, 6.5].
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
        CALL AIRY_NEG_CHEB15_FAST(X, AI, AIP)
        RETURN
    END IF

    IF (X >= 2.09E0_dl .AND. X <= 6.5E0_dl) THEN
        CALL AIRY_POS_CHEB10_FAST(X, AI, AIP)
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

    SUBROUTINE AIRY_NEG_CHEB15_FAST(XV, AIV, AIPV)
    IMPLICIT NONE
    REAL(dl), INTENT(IN)  :: XV
    REAL(dl), INTENT(OUT) :: AIV, AIPV

    INTEGER, PARAMETER :: NCH = 16

    ! Maps XV in [-6.5, -2.09] to Y in [-1, 1].
    REAL(dl), PARAMETER :: XMID  = -4.295E0_dl
    REAL(dl), PARAMETER :: XHALF =  2.205E0_dl

    REAL(dl) :: Y, TWOY
    REAL(dl) :: A0, A1, A2
    REAL(dl) :: P0, P1, P2
    INTEGER :: J

    REAL(dl), PARAMETER :: CAI(16) = (/ &
        -9.55550193128158892E-02_dl,  7.41177823586725154E-02_dl, &
        -7.87017962013092653E-02_dl,  2.44917931901850044E-01_dl, &
        1.61840769749720559E-01_dl, -1.42496606118389518E-01_dl, &
        -1.96787835450132759E-02_dl,  3.03118074085087612E-02_dl, &
        -2.53992015279218145E-03_dl, -2.92924621494435305E-03_dl, &
        7.10714331325374348E-04_dl,  1.11224040935810749E-04_dl, &
        -6.35725583267568239E-05_dl,  3.03034309182473455E-06_dl, &
        2.77692484870453487E-06_dl, -5.21087238783824334E-07_dl /)

    REAL(dl), PARAMETER :: CAIP(16) = (/ &
        1.28554846638544606E-01_dl,  3.24671075566074663E-01_dl, &
        1.89882679799835141E-01_dl,  4.67440773890445471E-01_dl, &
        -4.76560672314042921E-01_dl, -1.19736622026907338E-01_dl, &
        1.69682439333981833E-01_dl, -1.26412013737734301E-02_dl, &
        -2.27734807200423474E-02_dl,  5.78905823832411796E-03_dl, &
        1.13873327950330771E-03_dl, -6.57330254424061703E-04_dl, &
        2.90149572409112524E-05_dl,  3.46159586568294970E-05_dl, &
        -6.71697944266892732E-06_dl, -6.46579102281296447E-07_dl /)

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

    END SUBROUTINE AIRY_NEG_CHEB15_FAST


    SUBROUTINE AIRY_POS_CHEB10_FAST(XV, AIV, AIPV)
    IMPLICIT NONE
    REAL(dl), INTENT(IN)  :: XV
    REAL(dl), INTENT(OUT) :: AIV, AIPV

    INTEGER, PARAMETER :: NCH = 10

    ! Maps XV in [2.09, 6.5] to Y in [-1, 1].
    REAL(dl), PARAMETER :: XMID  = 4.295E0_dl
    REAL(dl), PARAMETER :: XHALF = 2.205E0_dl

    REAL(dl) :: Y, TWOY
    REAL(dl) :: A0, A1, A2
    REAL(dl) :: P0, P1, P2
    INTEGER :: J

    REAL(dl), PARAMETER :: CAI(10) = (/ &
        6.56835345885114271E-03_dl, -1.13188735137325670E-02_dl, &
        7.28898522549794800E-03_dl, -3.54092182656727842E-03_dl, &
        1.29632234600740061E-03_dl, -3.47294929101679922E-04_dl, &
        6.03481025954715256E-05_dl, -2.54640404330824147E-06_dl, &
        -2.29226532514133926E-06_dl,  8.50239953409589462E-07_dl /)

    REAL(dl), PARAMETER :: CAIP(10) = (/ &
        -1.07428727333418600E-02_dl,  1.82362021582016497E-02_dl, &
        -1.12191935268355983E-02_dl,  5.01355322437066573E-03_dl, &
        -1.58403210409663363E-03_dl,  3.10342928248546947E-04_dl, &
        -8.99822672257408064E-06_dl, -1.80849377841010130E-05_dl, &
        7.18604938632310802E-06_dl, -1.48852216585817192E-06_dl /)

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

    END SUBROUTINE AIRY_POS_CHEB10_FAST


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



    SUBROUTINE BJL_approx(L,X,JL)
    !!== MODIFIED SUBROUTINE FOR SPHERICAL BESSEL FUNCTIONS.                       ==!!
    !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!!
    !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
    !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!!
    !!== zqhuang@astro.utoronto.ca                                                 ==!!
    IMPLICIT NONE
    INTEGER L
    real(dl) X,JL
    real(dl) AX,AX2
    real(dl),PARAMETER::LN2=0.6931471805599453094D0
    real(dl),PARAMETER::ONEMLN2=0.30685281944005469058277D0
    real(dl),PARAMETER::PID2=1.5707963267948966192313217D0
    real(dl),PARAMETER::PID4=0.78539816339744830961566084582D0
    real(dl),parameter::ROOTPI12 = 21.269446210866192327578D0
    real(dl),parameter::GAMMA1 =   2.6789385347077476336556D0 !/* Gamma function of 1/3 */
    real(dl),parameter::GAMMA2 =   1.3541179394264004169452D0 !/* Gamma function of 2/3 */
    real(dl),PARAMETER::PI=3.141592653589793238463D0
    real(dl) NU,NU2,BETA,BETA2,COSB
    real(dl) sx,sx2
    real(dl) cotb,cot3b,cot6b,secb,sec2b
    real(dl) trigarg,expterm,L3

    IF(L.LT.0)THEN
        error stop 'Can not evaluate Spherical Bessel Function with index l<0'
    ENDIF
    AX=DABS(X)
    AX2=AX**2
    IF(L.LT.7)THEN
        IF(L.EQ.0)THEN
            IF(AX.LT.1.D-1)THEN
                JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
            ELSE
                JL=DSIN(AX)/AX
            ENDIF

        ELSEIF(L.EQ.1)THEN
            IF(AX.LT.2.D-1)THEN
                JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
            ELSE
                JL=(DSIN(AX)/AX-DCOS(AX))/AX
            ENDIF
        ELSEIF(L.EQ.2)THEN
            IF(AX.LT.3.D-1)THEN
                JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
            ELSE
                JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
            ENDIF
        ELSEIF(L.EQ.3)THEN
            IF(AX.LT.4.D-1)THEN
                JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
            ELSE
                JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
            ENDIF
        ELSEIF(L.EQ.4)THEN
            IF(AX.LT.6.D-1)THEN
                JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
            ELSE
                JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
            ENDIF
        ELSEIF(L.EQ.5)THEN
            IF(AX.LT.1.D0)THEN
                JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
            ELSE
                JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
            ENDIF
        ELSE
            IF(AX.LT.1.D0)THEN
                JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
            ELSE
                JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
                    DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
            ENDIF
        ENDIF
    ELSE
        NU=0.5D0+L
        NU2=NU**2
        IF(AX.LT.1.D-40)THEN
            JL=0.D0
        ELSEIF((AX2/L).LT.5.D-1)THEN
            JL=DEXP(L*DLOG(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
                /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
        ELSEIF((real(L,dl)**2/AX).LT.5.D-1)THEN
            BETA=AX-PID2*(L+1)
            JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
                -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
                (NU2-20.25)/80.D0/AX2)) )/AX
        ELSE
            L3=NU**0.325
            IF(AX .LT. NU-1.31*L3) then
                COSB=NU/AX
                SX = DSQRT(NU2-AX2)
                COTB=NU/SX
                SECB=AX/NU
                BETA=DLOG(COSB+SX/AX)
                COT3B=COTB**3
                COT6B=COT3B**2
                SEC2B=SECB**2
                EXPTERM=( (2.D0+3.D0*SEC2B)*COT3B/24.D0 &
                    - ( (4.D0+SEC2B)*SEC2B*COT6B/16.D0 &
                    + ((16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.D0 &
                    + (32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.D0/NU)*COT6B/NU) &
                    /NU)/NU
                JL=DSQRT(COTB*COSB)/(2.D0*NU)*DEXP(-NU*BETA+NU/COTB-EXPTERM)

                !          /**************** Region 2: x >> l ****************/

            ELSEIF (AX .GT. NU+1.48*L3) then
                COSB=NU/AX
                SX=DSQRT(AX2-NU2)
                COTB=NU/SX
                SECB=AX/NU
                BETA=DACOS(COSB)
                COT3B=COTB**3
                COT6B=COT3B**2
                SEC2B=SECB**2
                TRIGARG=NU/COTB-NU*BETA-PID4 &
                    -((2.0+3.0*SEC2B)*COT3B/24.D0  &
                    +(16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.D0/NU2)/NU
                EXPTERM=( (4.D0+sec2b)*sec2b*cot6b/16.D0 &
                    -(32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.D0/NU2)/NU2
                JL=DSQRT(COTB*COSB)/NU*DEXP(-EXPTERM)*DCOS(TRIGARG)

                !          /***************** Region 3: x near l ****************/

            ELSE
                BETA=AX-NU
                BETA2=BETA**2
                SX=6.D0/AX
                SX2=SX**2
                SECB=SX**0.3333333333333333d0
                SEC2B=SECB**2
                JL=( GAMMA1*SECB + BETA*GAMMA2*SEC2B &
                    -(BETA2/18.D0-1.D0/45.D0)*BETA*SX*SECB*GAMMA1 &
                    -((BETA2-1.D0)*BETA2/36.D0+1.D0/420.D0)*SX*SEC2B*GAMMA2   &
                    +(((BETA2/1620.D0-7.D0/3240.D0)*BETA2+1.D0/648.D0)*BETA2-1.D0/8100.D0)*SX2*SECB*GAMMA1 &
                    +(((BETA2/4536.D0-1.D0/810.D0)*BETA2+19.D0/11340.D0)*BETA2-13.D0/28350.D0)*BETA*SX2*SEC2B*GAMMA2 &
                    -((((BETA2/349920.D0-1.D0/29160.D0)*BETA2+71.D0/583200.D0)*BETA2-121.D0/874800.D0)* &
                    BETA2+7939.D0/224532000.D0)*BETA*SX2*SX*SECB*GAMMA1)*DSQRT(SX)/ROOTPI12
            ENDIF
        ENDIF
    ENDIF
    IF(X.LT.0.AND.MOD(L,2).NE.0)JL=-JL
    END SUBROUTINE BJL_approx


    end module FlatBessels



    SUBROUTINE BJL_EXTERNAL(L,X,JL)
    use FlatBessels
    use Precision
    !!== MODIFIED SUBROUTINE FOR SPHERICAL BESSEL FUNCTIONS.                       ==!!
    !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!!
    !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
    !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!!
    !!== zqhuang@astro.utoronto.ca                                                 ==!!
    IMPLICIT NONE
    INTEGER L
    real(dl) X,JL

    call BJL(L,X,JL)

    END SUBROUTINE BJL_EXTERNAL
