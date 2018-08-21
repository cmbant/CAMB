    !CAMB spherical and hyperspherical Bessel function routines
    !This version May 2006 - minor changes to bjl (http://cosmocoffee.info/viewtopic.php?t=530)
    !Feb 2007: fixed for high l, uses Ranges
    !Feb 2009: minor fix for non-flat compiled with non-smart IF evaluation
    !Dec 2011: minor tweak to DoRecurs for smoother errors across flat for L~O(30)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !Flat bessel function module

    module SpherBessels
    use Precision
    use ModelParams
    use Ranges
    implicit none
    private

    !     Bessel functions and their second derivatives for interpolation

    real(dl), dimension(:,:), allocatable ::  ajl,ajlpr, ddajlpr

    integer  num_xx, kmaxfile, file_numl,  file_l(lmax_arr)
    !      parameters for working out where the flat Bessel functions are small
    !      Both should increase for higher accuracy
    !        real(dl), parameter :: xlimmin=15._dl  , xlimfrac = 0.05_dl
    real(dl), parameter :: xlimmin=35._dl  , xlimfrac = 0.05_dl

    Type(Regions):: BessRanges

    public ajl, ajlpr, ddajlpr, BessRanges, InitSpherBessels, xlimmin, xlimfrac
    public USpherBesselWithDeriv, phi_recurs,phi_langer, bjl, Bessels_Free

    contains


    subroutine InitSpherBessels
    !     This subroutine reads the jl files from disk (or generates them if not on disk)
    use lvalues
    implicit none

    !See if already loaded with enough (and correct) lSamp%l values and k*eta values
    if (allocated(ajl) .and. (lSamp%l0 <= file_numl) .and. all(file_l(1:lSamp%l0)-lSamp%l(1:lSamp%l0)==0) &
        .and. (int(min(max_bessels_etak,CP%Max_eta_k))+1 <= kmaxfile)) return

    !Haven't made them before, so make them now
    call GenerateBessels

    if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Calculated Bessels'

    end subroutine InitSpherBessels

    subroutine GenerateBessels
    use lvalues
    real(dl) x
    real(dl) xlim
    integer i,j
    integer max_ix
    real(dl), parameter :: bessel_boost =1._dl


    if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'Generating flat Bessels...'


    file_numl= lSamp%l0
    file_l(1:lSamp%l0) = lSamp%l(1:lSamp%l0)
    kmaxfile = int(min(CP%Max_eta_k,max_bessels_etak))+1
    if (do_bispectrum) kmaxfile = kmaxfile*2


    call Ranges_Init(BessRanges)

    call Ranges_Add_delta(BessRanges,0._dl, 1._dl,0.01_dl/bessel_boost)
    call Ranges_Add_delta(BessRanges,1._dl, 5._dl,0.1_dl/bessel_boost)
    call Ranges_Add_delta(BessRanges,5._dl, 25._dl,0.2_dl/bessel_boost)
    call Ranges_Add_delta(BessRanges,25._dl, 150._dl,0.5_dl/bessel_boost/AccuracyBoost)
    call Ranges_Add_delta(BessRanges,150._dl, real(kmaxfile,dl),0.8_dl/bessel_boost/AccuracyBoost)

    call Ranges_GetArray(bessRanges, .false.)
    num_xx = BessRanges%npoints


    max_ix = min(max_bessels_l_index,lSamp%l0)

    if (allocated(ajl)) deallocate(ajl)
    if (allocated(ajlpr)) deallocate(ajlpr)
    if (allocated(ddajlpr)) deallocate(ddajlpr)
    Allocate(ajl(1:num_xx,1:max_ix))
    Allocate(ajlpr(1:num_xx,1:max_ix))
    Allocate(ddajlpr(1:num_xx,1:max_ix))

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC), PRIVATE(j,i,x,xlim)
    do j=1,max_ix

        do  i=1,num_xx
            x=BessRanges%points(i)
            xlim=xlimfrac*lSamp%l(j)
            xlim=max(xlim,xlimmin)
            xlim=lSamp%l(j)-xlim
            if (x > xlim) then
                if ((lSamp%l(j)==3).and.(x <=0.2) .or. (lSamp%l(j) > 3).and.(x < 0.5) .or. &
                    (lSamp%l(j)>5).and.(x < 1.0)) then
                    ajl(i,j)=0
                else
                    !if ( lSamp%l(j) > 40000) then
                    ! ajl(i,j) = phi_langer(lSamp%l(j),0,1._dl,x)
                    !else
                    call bjl(lSamp%l(j),x,ajl(i,j))
                    !end if
                end if
            else
                ajl(i,j)=0
            end if
        end do

        !     get the interpolation matrix for bessel functions
        call spline(BessRanges%points,ajl(1,j),num_xx,spl_large,spl_large,ajlpr(1,j))
        call spline(BessRanges%points,ajlpr(1,j),num_xx,spl_large,spl_large,ddajlpr(1,j))

    end do
    !$OMP END PARALLEL DO

    end subroutine GenerateBessels

    subroutine Bessels_Free

    if (allocated(ajl)) deallocate(ajl)
    if (allocated(ajlpr)) deallocate(ajlpr)
    if (allocated(ddajlpr)) deallocate(ddajlpr)
    call Ranges_Free(BessRanges)

    end  subroutine Bessels_Free


    SUBROUTINE BJL(L,X,JL)
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
    END SUBROUTINE BJL

    !    end module SpherBessels




    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    ! Calculation of ultraspherical Bessel functions.                      c
    ! Fortran version of the c program hyperjl.c by Arthur Kosowsky.       c
    ! WKB approx described in astro-ph/9805173                             c
    !                                                                      c
    ! Modifications by Anthony Challinor and Antony Lewis                  c
    ! Minor modifications to correct K=1 case outside [0,pi],              c
    ! the small chi approximations for lSamp%l=0 and lSamp%l=1, and                    c
    ! the quadratic approximation to Q(x) around Q(x)=0.                   c
    ! Bug fixed in downwards recursion (phi_recurs)                        c
    !                                                                      c
    ! The routine phi_recurs uses recursion relations to calculate         c
    ! the functions, which is accurate but relatively slow.                c
    !   ***NOT STABLE FOR K=1 or for all cases ***                         c
    !                                                                      c
    ! The routine phi_langer uses Langer's formula for a                   c
    ! uniform first-order asymptotic approximation in the open, closed     c
    ! and flat cases. This approximation is EXCELLENT for all lSamp%l >= 3.      c
    !                                                                      c
    ! The routine qintegral calculates the closed-form answer              c
    ! to the eikonal integral used in the WKB approximation.               c
    !                                                                      c
    ! The routine airy_ai returns the Airy function Ai(x) of the argument  c
    ! passed. It employs a Pade-type approximation away from zero and      c
    ! a Taylor expansion around zero. Highly accurate.                     c
    !                                                                      c
    ! The routines polevl and p1evl are auxiliary polynomial               c
    ! evaluation routines used in the airy function calculation.           c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



    ! module USpherBessels
    ! use Precision
    ! implicit none
    ! private


    !public USpherBesselWithDeriv, phi_recurs,phi_langer


    ! contains

    subroutine USpherBesselWithDeriv(closed,Chi,l,beta,y1,y2)
    !returns y1=ujl*sinhChi and y2=diff(y1,Chi)
    !aim for accuracy > 1% for all inputs
    real(dl) Chi,beta,y1,y2,sinhChi,cothChi
    real(dl) sin_K, cot_K
    integer l,K
    logical, intent(IN) :: closed
    logical DoRecurs

    if (closed) then
        sin_K = sin(Chi)
        cot_K= 1._dl/tan(Chi)
        K=1
    else
        sin_K=sinh(Chi)
        cot_K = 1._dl/tanh(Chi)
        K=-1
    end if

    sinhChi = sin_K
    cothChi = cot_K

    DoRecurs = ((l<=45*AccuracyBoost).OR.((.not.closed.or.(abs(Chi-pi/2)>0.2d0)).and.(beta*l<750) &
        .or.closed.and.(beta*l<4000)))

    !Deep in the tails the closed recursion relation is not stable
    !Added July 2003 to prevent problems with very nearly flat models
    if (DoRecurs .and. closed) then
        if  (Chi < asin(sqrt(l*(l+1._dl))/beta) - 2/beta) then
            if (phi_langer(l,K,beta,Chi) < 1e-7) then
                call phi_small_closed_int(l,beta,chi,y1,y2)
                return
            end if
        end if
    end if

    if (DoRecurs) then
        !use recursive evaluation where WKB is poor and recurs is fast anyway.
        y1=phi_recurs(l,K,beta,Chi)*sinhChi
        y2=y1*(l+1)*cothChi
        if (.not.closed.or.(l+1<nint(beta))) y2=y2 - &
            sqrt(beta**2-(K*(l+1)**2))*phi_recurs(l+1,K,beta,Chi)*sinhChi
        !of course we could get y2 much more quickly by modifying
        !phi_recurs to return l and l+1 for each beta,Chi...


    else !WKB approx
        y1=phi_langer(l,K,beta,Chi)*sinhChi
        y2=y1*(l+1)*cothChi
        if (.not.closed.or.(l+1<nint(beta))) y2=y2 - &
            sqrt(beta**2-(K*(l+1)**2))*phi_langer(l+1,K,beta,Chi)*sinhChi

    end if

    end subroutine USpherBesselWithDeriv



    !Calculates y1,y2 (noramlized to a value near the turning point)
    !by integrating up the differential equation and normalizing to phi_recurs
    !in the region in which phi_recurs is stable
    !This allows closed functions to be computed where chi << turning point
    subroutine phi_small_closed_int(l,beta,chi,y1,y2)
    integer, intent(IN) :: l
    real(dl), intent(IN) :: beta, chi
    real(dl) y1,y2

    integer nsteps,i

    real(dl) ap1,nu2,dydchi1,dydchi2,yt1,yt2,dyt1,dyt2,dym1,dym2
    real(dl) x0, delchi,sh, h6,x
    real(dl) y1_x,y2_x, tmp,xh,hh

    nsteps = 200
    ap1 = l*(l+1)
    x0 = sqrt(ap1)/beta
    nu2 = beta**2

    if ((beta*chi)**2/l < 0.005) then
        !Series solution

        x = chi
        sh = sin(x)
        tmp=(ap1/sh**2 - nu2)
        y1=1e-20
        y2 = ((l+1)/x - (nu2-ap1/3)/(2*l+3)*x) * y1
    else

        x = max(1d-7,chi - 50._dl/l)
        y1=1e-20
        y2 = (l+1)*y1/x

        delchi = (chi-x)/nSteps
        h6=delchi/6
        hh=delchi/2
        sh = sin(x)
        tmp=(ap1/sh**2 - nu2)

        do i=1,nSteps
            ! One step in the ujl integration
            ! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2
            xh=x+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/sin(xh)**2 - nu2)


            dyt2=tmp*yt1       !deriv y2 at mid

            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid

            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+delchi*dym1 !y1 at end
            dym1=dyt1+dym1
            yt2=y2+delchi*dym2 !y2 at end
            dym2=dyt2+dym2

            x=x+delchi     !end point
            sh=sin(x)
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2*dym2)
            if (y1 > 1d10 .or. y2> 1d10) then

                y1=y1/1d10
                y2=y2/1d10

            end if
        end do

    end if

    y1_x = y1; y2_x = y2

    delchi = (x0 - chi)/nSteps
    h6=delchi/6
    hh=delchi/2

    do i=1,nSteps
        ! One step in the ujl integration
        ! fourth-order Runge-Kutta method to integrate equation for ujl

        dydchi1=y2         !deriv y1
        dydchi2=tmp*y1     !deriv y2
        xh=x+hh          !midpoint of step
        yt1=y1+hh*dydchi1  !y1 at midpoint
        yt2=y2+hh*dydchi2  !y2 at midpoint
        dyt1=yt2           !deriv y1 at mid
        tmp=(ap1/sin(xh)**2 - nu2)


        dyt2=tmp*yt1       !deriv y2 at mid

        yt1=y1+hh*dyt1     !y1 at mid
        yt2=y2+hh*dyt2     !y2 at mid

        dym1=yt2           !deriv y1 at mid
        dym2=tmp*yt1       !deriv y2 at mid
        yt1=y1+delchi*dym1 !y1 at end
        dym1=dyt1+dym1
        yt2=y2+delchi*dym2 !y2 at end
        dym2=dyt2+dym2

        x=x+delchi     !end point
        sh=sin(x)
        dyt1=yt2           !deriv y1 at end
        tmp=(ap1/sh**2 - nu2)
        dyt2=tmp*yt1       !deriv y2 at end
        y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
        y2=y2+h6*(dydchi2+dyt2+2*dym2)
        if (y1 > 1d10 .or. y2 > 1d10) then
            y1=y1/1d10
            y2=y2/1d10
            y1_x = y1_x/1d10
            y2_x = y2_x/1d10

        end if
    end do


    tmp = phi_recurs(l,1,beta,x0)*sin(x0) / y1
    y1 = y1_x * tmp
    y2 = y2_x * tmp


    end subroutine phi_small_closed_int

    !***********************************************************************
    !                                                                      *
    ! Calculates Phi(l,beta,chi) using recursion on l.                     *
    ! See Abbot and Schaefer, ApJ 308, 546 (1986) for needed               *
    ! recursion relations and closed-form expressions for l=0,1.           *
    ! (Note: Their variable y is the same as chi here.)                    *
    !                                                                      *
    ! When the flag direction is negative, downwards recursion on l        *
    ! must be used because the upwards direction is unstable to roundoff   *
    ! errors. The downwards recursion begins with arbitrary values and     *
    ! continues downwards to l=1, where the desired l value is normalized  *
    ! using the closed form solution for l=1. (See, e.g., Numerical        *
    ! Recipes of Bessel functions for more detail)                         *
    !                                                                      *
    !***********************************************************************

    function phi_recurs(l, K, beta, chi)
    !doesn't like values which give exponentially small phi
    integer, intent(IN) :: l, K
    real(dl), intent(IN) :: beta, chi
    real(dl) phi_recurs
    integer j, direction, lstart,ibeta
    real(dl) ell, kay, arg, answer,beta2
    real(dl) root_K
    real(dl) phi0, phi1, phi_plus, phi_zero, phi_minus, b_zero, b_minus
    real(dl), parameter :: ACC=40._dl, BIG=1.d10
    real(dl) sin_K, cot_K

    ell=dble(l)


    ! Test input values

    if(l<0) then
        call MpiStop("Bessel function index ell < 0")
    endif
    if(beta<0._dl) then
        call MpiStop("Wavenumber beta < 0")
    endif
    if ((abs(K)/=1).and.(K/=0)) then
        call MpiStop("K must be 1, 0 or -1")
    end if

    if(K==1) then
        ibeta=nint(beta)
        if(ibeta<3) then
            call MpiStop("Wavenumber beta < 3 for K=1")
        endif
        if(ibeta<=l) then
            call MpiStop("Wavenumber beta <= l")
        endif
    endif

    if (chi<1/BIG) then
        phi_recurs=0
        return
    end if

    kay = dble(K)
    arg = beta * chi
    beta2 = beta**2

    if(K == 0) then
        cot_K = 1._dl/chi
        sin_K = chi
        root_K = beta
    else
        root_K = sqrt(beta2 -kay*ell*ell)

        if(K == -1) then
            cot_K = 1._dl/tanh(chi)
            sin_K = sinh(chi)
        else
            cot_K = 1._dl/tan(chi)
            sin_K = sin(chi)
        end if

    endif


    ! Closed form solution for l=0

    if (abs(chi) < 1.d-4) then
        if (abs(arg)<1.d-4) then
            phi0 = 1._dl-chi**2*(beta*beta-kay)/6._dl
        else
            phi0=sin(arg)/arg
        end if
    else
        phi0 = sin(arg) / (beta * sin_K)
    end if

    if (l==0) then
        phi_recurs=phi0
        return
    end if


    ! Closed form solution for l=1

    if((abs(chi) < 1.d-4).and.(K/=0)) then
        if(arg < 1.d-4) then
            phi1 = chi*sqrt(beta*beta-kay)/3._dl
            !beta2 * chi / (3._dl * sqrt(1._dl+ kay * beta2))
        else
            phi1 = (sin(arg)/arg-cos(arg))/(sqrt(beta*beta-kay)*chi)
            !(sin(arg)/arg - cos(arg))/arg
        end if
    elseif ((abs(arg) < 1.d-4).and.(K == 0)) then
        phi1 = arg / 3._dl
    else
        if (K /= 0 ) then
            phi1 = sin(arg) * cot_K / (beta * sin_K) - cos(arg) / sin_K
            phi1 = phi1/sqrt(beta2 - kay)
        else
            phi1 = (sin(arg)/arg - cos(arg))/arg
        end if
    end if
    if(l==1) then
        phi_recurs=phi1
        return
    end if
    ! Find recursion direction
    !  direction = +1 for upward recursion, -1 for downward

    if(abs(cot_K) < root_K / ell) then
        direction = 1
    else
        direction = -1
    end if

    ! For K=1, must do upwards recursion:
    ! NOT STABLE for all values of chi

    if(K==1) direction = 1

    ! Do upwards recursion on l


    if(direction == 1)then
        b_minus = sqrt(beta2 - kay)
        phi_minus = phi0
        phi_zero = phi1

        do j=2,l

            if(K == 0) then
                phi_plus = ((2*j-1) * cot_K * phi_zero - beta*phi_minus)/ beta
            else
                b_zero = sqrt(beta2 - (K*j*j))
                phi_plus = ((2*j-1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero
                b_minus = b_zero
            end if
            phi_minus = phi_zero
            phi_zero = phi_plus
        end do


        phi_recurs=phi_plus


        return

        ! Do downwards recursion on l

    else
        lstart = l + 2 * int(sqrt(ell*ACC))

        b_zero = sqrt(beta2 - dble(K*lstart*lstart))
        phi_plus = 0._dl
        phi_zero = 1._dl
        answer = 0._dl

        do j= lstart - 2,1,-1

            if(K == 0) then
                phi_minus = ((2*j + 3) * cot_K * phi_zero - beta * phi_plus) / beta
            else
                b_minus = sqrt(beta2 - (K*(j+1)**2))
                phi_minus = ((2*j + 3) * cot_K * phi_zero - b_zero * phi_plus) / b_minus
                b_zero = b_minus
            end if
            phi_plus = phi_zero
            phi_zero = phi_minus

            if(j == l) answer = phi_minus
            if((abs(phi_zero) > BIG).and.(j/=1)) then
                phi_plus = phi_plus/BIG
                phi_zero = phi_zero/BIG
                answer = answer/BIG
            end if
        end do

        ! Normalize answer to previously computed phi1

        answer = answer*phi1 / phi_minus
        phi_recurs=answer

    end if

    end function phi_recurs


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    ! Calculates Phi(l,beta,chi) using the Langer uniform approximation    c
    ! to the first-order WKB approximation.                                c
    ! See C.M. Bender and S.A. Orszag,  Mathematical Methods for           c
    ! Scientists and Engineers (McGraw-Hill, 1978; LC QA371.B43),          c
    ! chapter 10.                                                          c
    !                                                                      c
    ! Differential equation for needed function can be cast into the       c
    ! Schrodinger form      \epsilon^2 y'' = Q(x) y                        c
    ! where \epsilon^2 = 1/l(l+1) and Q(x) depends on the parameter        c
    ! alpha \equiv beta * epsilon.                                         c
    !                                                                      c
    ! In the K= +1 case, the function is                                   c
    ! determined by its value on the interval [0, pi/2] and the symmetry   c
    ! conditions Phi(chi + pi) = (-1)^{beta - l - 1} Phi(chi),             c
    !            Phi(pi - chi) = (-1)^{beta - l - 1} Phi(chi).             c
    ! This interval contains one turning point, so the Langer formula      c
    ! can be used.                                                         c
    ! Note that the second condition at chi = pi/2 gives an eigenvalue     c
    ! condition on beta, which  must corrected. For the lowest             c
    ! eigenvalue(s), the region between the turning points is not large    c
    ! enough for the asymptotic solution to be valid, so the functions     c
    ! have a small discontinuity or discontinuous derivative at pi/2;      c
    ! this behavior is corrected by employing a 4-term asymptotic          c
    ! series around the regular point chi=pi/2.                            c
    ! The exact eigenvalue condition requires that beta must be an         c
    ! integer >= 3 with beta > l. Note this implies alpha > 1.             c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



    function phi_langer(l,K,beta,chi)
    integer l,K,ibeta,kay
    real(dl) phi_langer
    real(dl) ell,symm, anu, alpha2
    real(dl) beta,chi,eikonal, wkb, arg, arg2, tmp
    real(dl) epsilon, alpha, chi0, x, a, b,achi

    real(dl) cot_K, sin_K
    real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
        PIOVER2=1.570796327d0

    ell=dble(l)
    achi=chi

    symm=1._dl
    !
    ! Test input values
    !
    if(l<0) call MpiStop("Bessel function index ell < 0")
    if(beta<0._dl) call MpiStop("Wavenumber beta < 0")
    if ((abs(K)/=1).and.(K/=0)) then
        call MpiStop("K must be 1, 0 or -1")
    end if


    if(K == 1) then
        ibeta=nint(beta)
        if(ibeta<3) call MpiStop("Wavenumber beta < 3 for K=1")
        if(ibeta<=l) call MpiStop("Wavenumber beta <= l")
    endif

    kay=K


    ! For closed case, find equivalent chi in [0,pi/2]
    !
    if(K==1) then
        achi=achi-2._dl*Pi*int(achi/2._dl/PI)
        if(achi>PI) then
            achi=2._dl*PI-achi
            if(2*(l/2).eq.l) then
                symm=symm
            else
                symm=-symm
            endif
        endif
        if(achi>PI/2._dl) then
            achi=PI-achi
            if(2*((ibeta-l-1)/2).eq.(ibeta-l-1)) then
                symm=symm
            else
                symm=-symm
            endif
        endif
    endif

    ! Definitions
    if(K == 0) then
        sin_K = achi
    else
        if(K == -1) then
            sin_K = sinh(achi)
        else
            sin_K = sin(achi)
        end if
    endif

    ! Closed form solution for l=0
    !
    if(l == 0) then
        arg=beta*achi

        if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
                wkb=1._dl-achi*achi*(beta*beta-kay)/6._dl
            else
                wkb=sin(arg)/arg
            endif
        else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=1._dl-arg*arg/6._dl
        else
            wkb=sin(arg)/(beta*sin_K)
        endif
        phi_langer=symm*wkb
        return
    endif
    !
    ! Closed form solution for l=1
    !
    if(l==1) then
        arg=beta*achi

        if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
                wkb=achi*sqrt(beta*beta-kay)/3._dl
            else
                wkb=(sin(arg)/arg-cos(arg))/(sqrt(beta*beta-kay)*achi)
            endif
        else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=arg/3._dl
        else
            if(K/=0) then
                if(K==1) then
                    cot_K=1._dl/tan(achi)
                else
                    cot_K=1._dl/tanh(achi)
                endif
                wkb=sin(arg)*cot_K/(beta*sin_K)-cos(arg)/sin_K
                wkb=wkb/sqrt(beta*beta-kay)
            else
                wkb=(sin(arg)/arg-cos(arg))/arg
            endif
        end if
        phi_langer=symm*wkb
        return
    endif
    !
    ! Closed form solution for K=1 and beta = l+1 (lowest eigenfunction)
    !
    if((K==1).and.(ibeta == (l+1))) then
        wkb=(sin_K**ell)* &
            sqrt(sqrt(2._dl*PI/(2._dl*ell+1._dl))*ell/((ell+1._dl)*(2._dl*ell+1._dl)))
        wkb=wkb*(1+0.1875d0/ell-0.013671875/(ell*ell))
        phi_langer=symm*wkb
        return
    endif

    ! Very close to 0, return 0 (exponentially damped)
    !
    if(abs(achi)<1.d-8) then
        phi_langer=0._dl
        return
    endif


    ! For closed case, find corrected eigenvalue beta
    !
    if(K==1) then
        anu=dble(ibeta)-1._dl/(8._dl*ell)+1._dl/(16._dl*ell*ell)
    else
        anu=beta
    endif
    !
    ! Evaluate epsilon using asymptotic form for large l
    !
    if(l<20) then
        epsilon=1._dl/sqrt(ell*(ell+1._dl))
    else
        epsilon=1._dl/ell-0.5d0/(ell*ell)+0.375d0/(ell*ell*ell)
    endif

    alpha=epsilon*anu
    !
    ! Calculate the turning point where Q(x)=0.
    ! Function in question has only a single simple turning point.
    !
    if(K==-1) chi0=log((1._dl+sqrt(1._dl+alpha*alpha))/alpha)
    if(K==0) chi0=1._dl/alpha
    if(K==1) chi0=asin(1._dl/alpha)


    ! Very close to chi0, use usual wkb form to avoid dividing by zero
    !
    if(abs(achi-chi0)<1.d-5) then

        ! Calculate coefficients of linear and quadratic terms in Q(x) expansion
        ! in the neighborhood of the turning point
        ! Q(chi)=a*(chi0-chi)+b*(chi0-chi)**2
        alpha2=alpha*alpha

        if(K==-1) then
            a=2._dl*alpha2*sqrt(alpha2+1._dl)
            b=3._dl*alpha2**2+2._dl*alpha2
        endif
        if(K==0) then
            a=2._dl*alpha2*alpha
            b=3._dl*alpha2**2
        endif
        if(K==1) then
            a=2._dl*alpha2*sqrt(alpha2-1._dl)
            b=3._dl*alpha2**2-2._dl*alpha2
        endif

        ! Dependent variable x for which Q(x)=0 at x=0
        ! x>0 is the evanescent region
        !
        x=chi0-achi
        !
        ! Argument of Airy function
        !
        arg=(x+b*x*x/(5._dl*a))/(epsilon*epsilon/a)**(0.333333333d0)
        !
        ! Evaluate Airy function
        !
        wkb=airy_ai(arg)
        !
        ! Rest of functional dependence
        !
        wkb=wkb*(1._dl-b*x/(5._dl*a))/sin_K
        !  Normalization factor:

        wkb=wkb*symm*ROOTPI*((a*epsilon)**(-0.1666667d0))*sqrt(epsilon/anu)
        phi_langer=wkb

        return
    endif


    ! Langer approximation.
    !
    ! Transport factor:
    !
    tmp=sqrt(abs(1._dl/(sin_K*sin_K)-alpha*alpha))

    ! Eikonal factor
    !
    eikonal=qintegral(sin_K,alpha,K)



    arg=(1.5d0*eikonal/epsilon)**(1._dl/3._dl)

    arg2=arg*arg
    if(achi>chi0) arg2=-arg2

    ! Multiply factors together

    wkb=airy_ai(arg2)*symm*ROOTPI*sqrt(arg*epsilon/anu/tmp)/Sin_K
    phi_langer=wkb

    end function phi_langer

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                       c
    ! Evaluates the exact solution to  the integral giving the WKB          c
    ! eikonal solution,   \int^x sqrt(abs(Q(x))) dx                         c
    !                                                                       c
    ! In the open case, this integral costs 1 or 2 square roots, an atan    c
    ! and a log; its evaluation will be roughly as expensive as the rest    c
    ! of the Phi routine. An analytic fit cannot be substantially faster    c
    ! because the dependence on alpha of the y-intercept of the linear      c
    ! region of the integrand contains a log and an atan, so at best a fit  c
    ! can only save the computation of the square roots.                    c
    !                                                                       c
    ! The integrals are very bland functions of chi and alpha and could     c
    ! be precomputed and cached to save computation time; interpolation     c
    ! on a relatively small number of points should be very accurate        c
    !                                                                       c
    ! Note that for the closed case, the variable arg must be between 0     c
    ! and alpha; for arg > alpha, symmetry properties reduce the needed     c
    ! integral to this case.                                                c
    !                                                                       c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function qintegral(sin_K,alpha, K)
    implicit none
    real(dl) qintegral, sin_K
    integer K
    real(dl) alpha,exact,arg, root1, root2, dummyarg

    real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
        PIOVER2=1.570796327d0

    arg=alpha*sin_K

    if(K==0) then
        if(arg>1._dl) then
            exact=sqrt(arg*arg-1._dl)-acos(1._dl/arg)
            qintegral=exact
            return
        else
            root1=sqrt(1._dl-arg*arg)
            exact=log((1._dl+root1)/arg)-root1
            qintegral=exact
            return
        endif
    else if(K==-1) then
        if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(arg*arg+alpha*alpha)
            exact=alpha/2._dl*log((2._dl*arg*arg+alpha*alpha-1._dl+ &
                2._dl*root1*root2)/(1._dl+alpha*alpha))+atan(root2/ &
                (alpha*root1))-PIOVER2
            qintegral=exact
            return
        else
            root1=sqrt((1._dl-arg*arg)*(arg*arg+alpha*alpha))
            exact=alpha/2._dl*atan(-2._dl*root1/(2*arg*arg+alpha*alpha- &
                1._dl))+0.5d0*log((2._dl*alpha*root1+2._dl*alpha*alpha+ &
                arg*arg*(1._dl-alpha*alpha))/(arg*arg*(1._dl+ &
                alpha*alpha)))
            if(2._dl*arg*arg+alpha*alpha-1._dl<0._dl) then
                exact=exact-alpha*PIOVER2
            endif
            qintegral=exact
            return
        endif
    else
        if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(alpha*alpha-arg*arg)
            exact=alpha/2._dl*atan(-2._dl*root1*root2/ &
                (2._dl*arg*arg-alpha*alpha-1._dl))- &
                atan(-root2/(root1*alpha))-PIOVER2
            if(2._dl*arg*arg-alpha*alpha-1._dl>0._dl) then
                exact=exact+alpha*PIOVER2
            endif
        else
            root1=sqrt((1._dl-arg*arg)*(alpha*alpha-arg*arg))
            dummyarg=alpha*(1._dl-arg*arg)/root1
            exact=0.5d0*log((1._dl+dummyarg)/(1._dl-dummyarg))- &
                alpha/2._dl*log((alpha*alpha-2._dl*arg*arg+1._dl+2._dl*root1)/ &
                (alpha*alpha-1._dl))
        endif
        qintegral=exact
        return
    endif

    end function qintegral

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !       Airy function                                                  c
    !                                                                      c
    ! Modified from original routine by Stephen Moshier, available         c
    ! as part of the Cephes library at www.netlib.com                      c
    ! Modifications: eliminates calculation of Bi(x), Ai'(x), Bi'(x)       c
    ! and translation to Fortran                                           c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! DESCRIPTION:
    !
    ! Solution of the differential equation
    !
    !       y"(x) = xy.
    !
    ! The function returns the two independent solutions Ai, Bi
    ! and their first derivatives Ai'(x), Bi'(x).
    !
    ! Evaluation is by power series summation for small x,
    ! by rational minimax approximations for large x.
    !
    !
    !
    ! ACCURACY:
    ! Error criterion is absolute when function <= 1, relative
    ! when function > 1, except * denotes relative error criterion.
    ! For large negative x, the absolute error increases as x^1.5.
    ! For large positive x, the relative error increases as x^1.5.
    !
    ! Arithmetic  domain   function  # trials      peak         rms
    ! IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
    ! IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
    ! IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
    ! IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
    ! IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
    ! IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
    ! DEC         -10, 0     Ai         5000       1.7e-16     2.8e-17
    ! DEC           0, 10    Ai         5000       2.1e-15*    1.7e-16*
    ! DEC         -10, 0     Ai'        5000       4.7e-16     7.8e-17
    ! DEC           0, 10    Ai'       12000       1.8e-15*    1.5e-16*
    ! DEC         -10, 10    Bi        10000       5.5e-16     6.8e-17
    ! DEC         -10, 10    Bi'        7000       5.3e-16     8.7e-17
    !
    !
    ! Cephes Math Library Release 2.1:  January, 1989
    ! Copyright 1984, 1987, 1989 by Stephen lSamp%l. Moshier
    ! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
    !

    function airy_ai(x)
    implicit none
    real(dl) airy_ai
    real(dl) x,z, zz, t, f, g, uf, ug, zeta, theta
    real(dl) ak
    real(dl) AN(8),AD(8),AFN(9),AFD(9),AGN(11),AGD(10)
    real(dl), parameter :: AMAXAIRY=25.77d0,ACC=1.d-8,PI=3.1415926536d0
    real(dl), parameter :: c1=0.35502805388781723926d0, c2=0.258819403792806798405d0
    real(dl), parameter :: sqrt3=1.732050807568877293527d0,sqpii=5.64189583547756286948d-1


    AN(1)=3.46538101525629032477d-1
    AN(2)=1.20075952739645805542d1
    AN(3)=7.62796053615234516538d1
    AN(4)=1.68089224934630576269d2
    AN(5)=1.59756391350164413639d2
    AN(6)=7.05360906840444183113d1
    AN(7)=1.40264691163389668864d1
    AN(8)=9.99999999999999995305d-1

    AD(1)=5.67594532638770212846d-1
    AD(2)=1.47562562584847203173d1
    AD(3)=8.45138970141474626562d1
    AD(4)=1.77318088145400459522d2
    AD(5)=1.64234692871529701831d2
    AD(6)=7.14778400825575695274d1
    AD(7)=1.40959135607834029598d1
    AD(8)=1.00000000000000000470d0

    AFN(1)=-1.31696323418331795333d-1
    AFN(2)=-6.26456544431912369773d-1
    AFN(3)=-6.93158036036933542233d-1
    AFN(4)=-2.79779981545119124951d-1
    AFN(5)=-4.91900132609500318020d-2
    AFN(6)=-4.06265923594885404393d-3
    AFN(7)=-1.59276496239262096340d-4
    AFN(8)=-2.77649108155232920844d-6
    AFN(9)=-1.67787698489114633780d-8

    AFD(1)=1.33560420706553243746d1
    AFD(2)=3.26825032795224613948d1
    AFD(3)=2.67367040941499554804d1
    AFD(4)=9.18707402907259625840d0
    AFD(5)=1.47529146771666414581d0
    AFD(6)=1.15687173795188044134d-1
    AFD(7)=4.40291641615211203805d-3
    AFD(8)=7.54720348287414296618d-5
    AFD(9)=4.51850092970580378464d-7

    AGN(1)=1.97339932091685679179d-2
    AGN(2)=3.91103029615688277255d-1
    AGN(3)=1.06579897599595591108d0
    AGN(4)=9.39169229816650230044d-1
    AGN(5)=3.51465656105547619242d-1
    AGN(6)=6.33888919628925490927d-2
    AGN(7)=5.85804113048388458567d-3
    AGN(8)=2.82851600836737019778d-4
    AGN(9)=6.98793669997260967291d-6
    AGN(10)=8.11789239554389293311d-8
    AGN(11)=3.41551784765923618484d-10

    AGD(1)=9.30892908077441974853d0
    AGD(2)=1.98352928718312140417d1
    AGD(3)=1.55646628932864612953d1
    AGD(4)=5.47686069422975497931d0
    AGD(5)=9.54293611618961883998d-1
    AGD(6)=8.64580826352392193095d-2
    AGD(7)=4.12656523824222607191d-3
    AGD(8)=1.01259085116509135510d-4
    AGD(9)=1.17166733214413521882d-6
    AGD(10)=4.91834570062930015649d-9
    !
    ! Exponentially tiny for large enough argument
    !
    if(x>AMAXAIRY) then
        airy_ai=0._dl
        return
    endif
    !
    ! Pade fit for large negative arguments
    !
    if(x<-2.09d0) then
        t=sqrt(-x)
        zeta=-2._dl*x*t/3._dl
        t=sqrt(t)
        ak=sqpii/t
        z=1._dl/zeta
        zz=z*z
        uf=1._dl+zz*polevl(zz,AFN,8)/p1evl(zz,AFD,9)
        ug=z*polevl(zz,AGN,10)/p1evl(zz,AGD,10)
        theta=zeta+0.25d0*PI
        f=sin(theta)
        g=cos(theta)
        airy_ai=ak*(f*uf-g*ug)
        return
    endif
    !
    ! Pade fit for large positive arguments
    !
    if(x>=2.09) then
        t=sqrt(x)
        zeta=2._dl*x*t/3._dl
        g=exp(zeta)
        t=sqrt(t)
        ak=2._dl*t*g
        z=1._dl/zeta
        f=polevl(z,AN,7)/polevl(z,AD,7)
        airy_ai=sqpii*f/ak
        return
    endif
    !
    ! Taylor series for region around x=0
    !

    f=1._dl
    g=x
    t=1._dl
    uf=1._dl
    ug=x
    ak=1._dl
    z=x*x*x
    do while (t>ACC)
        uf=uf*z
        ak=ak+1._dl
        uf=uf/ak
        ug=ug*z
        ak=ak+1._dl
        ug=ug/ak
        uf=uf/ak
        f=f+uf
        ak=ak+1._dl
        ug=ug/ak
        g=g+ug
        t=abs(uf/f)
    end do


    airy_ai=c1*f-c2*g

    end function airy_ai

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  Evaluate polynomial                                          c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! DESCRIPTION:
    !
    ! Evaluates polynomial of degree N:
    !
    !                     2          N
    ! y  =  C  + C x + C x  +...+ C x
    !        0    1     2          N
    !
    ! Coefficients are stored in reverse order:
    !
    ! coef(1) = C  , ..., coef(N+1) = C  .
    !            N                     0
    !
    ! The function p1evl() assumes that C = 1.0 and is
    !                                    N
    ! omitted from the array.  Its calling arguments are
    ! otherwise the same as polevl().
    !
    !
    ! SPEED:
    !
    ! In the interest of speed, there are no checks for out
    ! of bounds arithmetic.  This routine is used by most of
    ! the functions in the library.  Depending on available
    ! equipment features, the user may wish to rewrite the
    ! program in microcode or assembly language.
    !
    ! Cephes Math Library Release 2.1:  December, 1988
    ! Copyright 1984, 1987, 1988 by Stephen lSamp%l. Moshier
    ! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
    !
    function polevl(x,coef,N)
    implicit none
    real(dl) polevl
    real(dl) x,ans
    real(dl) coef
    integer N,i

    dimension coef(N+1)

    ans=coef(1)
    do i=2,N+1
        ans=ans*x+coef(i)
    end do
    polevl=ans

    end function polevl

    !
    !
    ! Evaluate polynomial when coefficient of x  is 1.0.
    ! Otherwise same as polevl.
    !
    function p1evl(x,coef,N)
    implicit none
    real(dl) p1evl
    real(dl) x,coef,ans
    integer N,i
    dimension coef(N)

    ans=x+coef(1)
    do i=2,N
        ans=ans*x+coef(i)
    end do
    p1evl=ans

    end function p1evl



    end module SpherBessels !USpherBessels



    SUBROUTINE BJL_EXTERNAL(L,X,JL)
    use SpherBessels
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

