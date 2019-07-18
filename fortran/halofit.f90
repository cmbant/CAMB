    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The `halofit' code models the nonlinear evolution of cold matter
    ! cosmological power spectra. The full details of the way in which
    ! this is done are presented in Smith et al. (2002), MNRAS, ?, ?.
    !
    ! The code `halofit' was written by R. E. Smith & J. A. Peacock.
    ! See http://www.astro.upenn.edu/~res,
    !
    ! Subsequent updates as below
    ! Only tested for basic models with power law initial power spectra
    ! References for variant versions are
    !   halofit_original: astro-ph/0207664
    !   halofit_peacock: http://www.roe.ac.uk/~jap/haloes/
    !   halofit_bird: arXiv: 1109.4416
    !   halofit_takahashi: arXiv: 1208.2701
    !   halofit_mead: arXiv:1505.07833,1602.02154
    !   halofit_casarini: arXiv:0810.0190, arXiv:1601.07230

    ! Adapted for F90 and CAMB, AL March 2005
    !!BR09 Oct 09: generalized expressions for om(z) and ol(z) to include w

    ! RT12 Oct: update some fitting parameters in the code to enhance
    !           the power spectrum at small scales (arXiv:1208.2701)

    !!JD 08/13: generalized expressions for om(z) and ol(z) to include
    !           w_0 and w_a
    ! SPB14 Feb: update the fitting parameters for neutrinos to work with RT12
    !           modifications
    ! AL Sept 14: added halofit_version parameter to change approximation used;
    !   separate halofit.f90 is no longer needed as equations.f90 defined fixed wa_ppf
    ! Jan 15: Suggested change from Simeon Bird to avoid issues with very large Omm and neutrinos
    !AM Mar 16: Added in HMcode
    !AM May 16: Fixed some small bugs and added better neutrino approximations
    !AL Jun16: put in partial openmp for HMcode (needs restructure to do properly)
    !AM Sep 16: Attempted fix of strange bug. No more modules with unallocated arrays as inputs
    !LC Oct 16: extended Halofit from w=const. models to w=w(a) with PKequal
    !AM May 17: Made the baryon feedback parameters more obvious in HMcode
    !AL Jul 17: fixed undefined z calling Tcb_Tcbnu_ratio
    !AM Jul 17: sped-up HMcode integration routines
    !AM May 18: Fixed bug in Dolag correction to c(M) power
    !AM Jul 19: Upgraded accuracy and bug fix for massive-neutrino models

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module NonLinear
    use results
    use DarkEnergyInterface
    use classes
    use Transfer
    implicit none
    private

    integer, parameter :: halofit_original = 1, halofit_bird=2, halofit_peacock=3, halofit_takahashi=4
    integer, parameter :: halofit_mead=5, halofit_halomodel=6, halofit_casarini=7, halofit_mead2015=8
    integer, parameter :: halofit_default = halofit_mead

    logical :: HM_verbose = .false.

    type, extends(TNonLinearModel) :: THalofit
        integer :: halofit_version = halofit_default
        !!TT - These are the baryon parameters of HMCode
        real(dl) :: HMCode_A_baryon=3.13
        real(dl) :: HMCode_eta_baryon=0.603
        !!AM - Added these types for HMcode
        integer, private :: imead !!AM - added these for HMcode, need to be visible to all subroutines and functions
        real(dl), private :: om_m,om_v,fnu,omm0, acur, w_hf, wa_hf
        real(dl), private :: om_c, om_b
    contains
    procedure :: ReadParams => THalofit_ReadParams
    procedure :: GetNonLinRatios => THalofit_GetNonLinRatios
    procedure :: halofit
    procedure :: HMcode
    procedure, nopass :: PythonClass => THalofit_PythonClass
    procedure, nopass :: SelfPointer => THalofit_SelfPointer
    procedure, private :: Delta_v
    procedure, private :: delta_c
    procedure, private :: eta
    procedure, private :: kstar
    procedure, private :: As
    procedure, private :: conc_bull
    procedure, private :: fdamp
    procedure, private :: p_1h
    procedure, private :: p_2h
    procedure, private :: alpha
    procedure, private :: halomod
    procedure, private :: halomod_init
    procedure, private :: write_parameters
    procedure, private :: zcoll_bull
    end type

    public THalofit, HM_verbose
    public halofit_default, halofit_original, halofit_bird, halofit_peacock, halofit_takahashi
    public halofit_mead, halofit_halomodel, halofit_casarini

    TYPE HM_cosmology
        !Contains only things that do not need to be recalculated with each new z
        REAL(dl) :: om_m, om_c, om_b, om_v, w, wa, f_nu, ns, h, Tcmb, Nnu
        REAL(dl), ALLOCATABLE :: log_r_sigma(:), log_sigma(:)
        REAL(dl), ALLOCATABLE :: growth(:), a_growth(:)
        REAL(dl), ALLOCATABLE :: log_k_plin(:), log_plin(:), log_plinc(:)
        INTEGER :: nk, ng, nsig
        !AM - Added feedback parameters below at fixed fiducial (DMONLY) values
        REAL(dl) :: A_baryon=3.13
        REAL(dl) :: eta_baryon=0.603
    END TYPE HM_cosmology

    TYPE HM_tables
        !Stuff that needs to be recalculated for each new z
        REAL(dl), ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
        REAL(dl) :: sigv, sigv100, knl, rnl, neff, sig8z, z, dc
        INTEGER :: n
    END TYPE HM_tables
    !!AM - End of my additions

    ! HMcode numerical parameters
    REAL(dl), PARAMETER :: mmin_HMcode=1e0 !Lower mass limit
    REAL(dl), PARAMETER :: mmax_HMcode=1e18 !Upper mass limit
    INTEGER, PARAMETER :: n_HMcode=256 !Number of entries in look-up HM_tables.

    ! HMcode linear P(k) numerical parameters
    INTEGER, PARAMETER :: imeth_pk=1 !AM Jul 19: Changed to imeth=1 from imeth=2
    REAL(dl), PARAMETER :: kmin_pk=1e-3
    REAL(dl), PARAMETER :: kmax_pk=1e2
    INTEGER, PARAMETER :: nk_pk=512 !AM Jul 19: Improved accuracy from 128 to 512 points

    ! HMcode numerical parameters for sigma(R)
    REAL(dl), PARAMETER :: rmin_sigma=1e-4
    REAL(dl), PARAMETER :: rmax_sigma=1e3
    INTEGER, PARAMETER :: n_sigma=64
    REAL(dl), PARAMETER :: acc_sigma=1e-4 !AM Jul 19: Upgraded accuracy from 1e-3 to 1e-4
    REAL(dl), PARAMETER :: alpha_sigma=3.
    INTEGER, PARAMETER :: iorder_sigma=3

    ! HMcode numerical parameters for sigmaV(R)
    REAL(dl), PARAMETER :: acc_sigmaV=1e-4
    REAL(dl), PARAMETER :: alpha_sigmaV=3.
    INTEGER, PARAMETER :: iorder_sigmaV=3

    ! HMcode numerical parameters for neff(R)
    REAL(dl), PARAMETER :: acc_neff=1e-4
    REAL(dl), PARAMETER :: alpha_neff=1.
    INTEGER, PARAMETER :: iorder_neff=3

    contains

    function THalofit_PythonClass()
    character(LEN=:), allocatable :: THalofit_PythonClass

    THalofit_PythonClass = 'Halofit'

    end function THalofit_PythonClass

    subroutine THalofit_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (THalofit), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine THalofit_SelfPointer

    subroutine THalofit_ReadParams(this,Ini)
    use IniObjects
    class(THalofit) :: this
    class(TIniFile), intent(in) :: Ini

    this%halofit_version = Ini%Read_Int('halofit_version', halofit_default)

    end subroutine THalofit_ReadParams

    subroutine THalofit_GetNonLinRatios(this,State,CAMB_Pk)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    !This implementation uses Halofit
    class(THalofit) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    integer itf
    real(dl) a,plin,pq,ph,pnl,rk
    real(dl) sig,rknl,rneff,rncur,d1,d2
    real(dl) diff,xlogr1,xlogr2,rmid, h2
    integer i

    select type (State)
    class is (CAMBdata)
        associate(Params => State%CP)

            IF(this%halofit_version==halofit_mead .OR. this%halofit_version==halofit_halomodel &
                .OR.  this%halofit_version==halofit_mead2015) THEN
                !AM - Call HMcode here
                CALL this%HMcode(State,CAMB_Pk)
            ELSE

                !!BR09 putting neutrinos into the matter as well, not sure if this is correct, but at least one will get a consisent omk.
                h2 = (Params%H0/100)**2
                this%omm0 = (Params%omch2+Params%ombh2+Params%omnuh2)/h2
                this%fnu = Params%omnuh2/h2/this%omm0

                CAMB_Pk%nonlin_ratio = 1

                do itf = 1, CAMB_Pk%num_z

                    call Params%DarkEnergy%Effective_w_wa(this%w_hf, this%wa_hf)
                    if (this%halofit_version == halofit_casarini) then
                        ! calculate equivalent w-constant models (w_hf,0) for w_lam+wa_ppf(1-a) models
                        ! [Casarini+ (2009,2016)].
                        call PKequal(State,CAMB_Pk%Redshifts(itf),this%w_hf,this%wa_hf,this%w_hf,this%wa_hf)
                    endif

                    ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and
                    ! curvature (rncur) of the power spectrum at the desired redshift, using method
                    ! described in Smith et al (2002).
                    a = 1/real(1+CAMB_Pk%Redshifts(itf),dl)
                    this%om_m = omega_m(a, this%omm0, State%omega_de, this%w_hf, this%wa_hf)
                    this%om_v = omega_v(a, this%omm0, State%omega_de, this%w_hf, this%wa_hf)
                    this%acur = a
                    xlogr1=-2.0
                    xlogr2=3.5
                    do
                        rmid=(xlogr2+xlogr1)/2.0
                        rmid=10**rmid
                        call wint(CAMB_Pk,itf,rmid,sig,d1,d2)
                        diff=sig-1.0
                        if (abs(diff).le.0.001) then
                            rknl=1./rmid
                            rneff=-3-d1
                            rncur=-d2
                            exit
                        elseif (diff.gt.0.001) then
                            xlogr1=log10(rmid)
                        elseif (diff.lt.-0.001) then
                            xlogr2=log10(rmid)
                        endif
                        if (xlogr2 < -1.9999) then
                            !is still linear, exit
                            goto 101
                        else if (xlogr1>3.4999) then
                            ! Totally crazy non-linear
                            global_error_flag=349
                            write(*,*) 'Error in halofit'
                            goto 101
                        end if
                    end do

                    ! now calculate power spectra for a logarithmic range of wavenumbers (rk)

                    do i=1, CAMB_PK%num_k
                        rk = exp(CAMB_Pk%log_kh(i))

                        if (rk > this%Min_kh_nonlinear) then

                            ! linear power spectrum !! Remeber => plin = k^3 * P(k) * constant
                            ! constant = 4*pi*V/(2*pi)^3

                            plin= MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3/(2*const_pi**2))

                            ! calculate nonlinear power according to halofit: pnl = pq + ph,
                            ! where pq represents the quasi-linear (halo-halo) power and
                            ! where ph is represents the self-correlation halo term.

                            call this%halofit(rk,rneff,rncur,rknl,plin,pnl,pq,ph)   ! halo fitting formula
                            CAMB_Pk%nonlin_ratio(i,itf) = sqrt(pnl/plin)

                        end if

                    enddo

101                 continue
                end do

            END IF
        end associate
    end select

    end subroutine THalofit_GetNonLinRatios

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subroutine halofit(this,rk,rn,rncur,rknl,plin,pnl,pq,ph)
    class(THalofit) :: this
    real(dl) gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
    real(dl) rk,rn,plin,pnl,pq,ph,plinaa
    real(dl) rknl,y,rncur
    real(dl) f1a,f2a,f3a,f1b,f2b,f3b,frac
    real(dl) extragam, peacock_fudge

    if (this%halofit_version ==halofit_original .or. this%halofit_version ==halofit_bird &
        .or. this%halofit_version == halofit_peacock) then
        ! halo model nonlinear fitting formula as described in
        ! Appendix C of Smith et al. (2002)
        !SPB11: Standard halofit underestimates the power on the smallest scales by a
        !factor of two. Add an extra correction from the simulations in Bird, Viel,
        !Haehnelt 2011 which partially accounts for this.
        if (this%halofit_version ==halofit_bird) then
            extragam = 0.3159 -0.0765*rn -0.8350*rncur
            gam=extragam+0.86485+0.2989*rn+0.1631*rncur
        else
            gam=0.86485+0.2989*rn+0.1631*rncur
        end if
        a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+ &
            0.1670756*rn*rn*rn*rn-0.620695*rncur
        a=10**a
        b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
        c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
        xmu=10**(-3.54419+0.19086*rn)
        xnu=10**(0.95897+1.2857*rn)
        alpha=1.38848+0.3701*rn-0.1452*rn*rn
        beta=0.8291+0.9854*rn+0.3400*rn**2+this%fnu*(-6.4868+1.4373*rn**2)
    elseif (this%halofit_version == halofit_takahashi .or. this%halofit_version == halofit_casarini) then
        !RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
        !than latest N-body simulations at small scales.
        !Update the following fitting parameters of gam,a,b,c,xmu,xnu,
        !alpha & beta from the simulations in Takahashi+ 2012.
        !The improved halofit accurately provide the power spectra for WMAP
        !cosmological models with constant w.
        !LC16 Jun: Casarini+ 2009,2016 extended constant w prediction for w(a).
        gam=0.1971-0.0843*rn+0.8460*rncur
        a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+ &
            0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*this%om_v*(1.+this%w_hf+this%wa_hf*(1-this%acur))
        a=10**a
        b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+ &
            0.2279*this%om_v*(1.+this%w_hf+this%wa_hf*(1-this%acur)))
        c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
        xmu=0.
        xnu=10**(5.2105+3.6902*rn)
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ &
            0.3980*rn**4-0.1682*rncur + this%fnu*(1.081 + 0.395*rn**2)
    else
        call MpiStop('Unknown halofit_version')
    end if

    if(abs(1-this%om_m).gt.0.01) then ! omega evolution
        f1a=this%om_m**(-0.0732)
        f2a=this%om_m**(-0.1423)
        f3a=this%om_m**(0.0725)
        f1b=this%om_m**(-0.0307)
        f2b=this%om_m**(-0.0585)
        f3b=this%om_m**(0.0743)
        frac=this%om_v/(1.-this%om_m)
        f1=frac*f1b + (1-frac)*f1a
        f2=frac*f2b + (1-frac)*f2a
        f3=frac*f3b + (1-frac)*f3a
    else
        f1=1.0
        f2=1.
        f3=1.
    endif

    y=(rk/rknl)


    ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
    ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+this%fnu*0.977)
    plinaa=plin*(1+this%fnu*47.48*rk**2/(1+1.5*rk**2))
    pq=plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)

    pnl=pq+ph

    if (this%halofit_version == halofit_peacock) then
        !From http://www.roe.ac.uk/~jap/haloes/
        !(P-P_linear) -> (P-P_linear) * (1+2y^2)/(1+y^2), where y = k/10 h Mpc^(-1).
        peacock_fudge = rk/10
        pnl = plin + (pnl-plin)*(1+2*peacock_fudge**2)/(1+peacock_fudge**2)
    end if

    end subroutine halofit


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The subroutine wint, finds the effective spectral quantities
    ! rknl, rneff & rncur. This it does by calculating the radius of
    ! the Gaussian filter at which the variance is unity = rknl.
    ! rneff is defined as the first derivative of the variance, calculated
    ! at the nonlinear wavenumber and similarly the rncur is the second
    ! derivative at the nonlinear wavenumber.

    subroutine wint(CAMB_Pk,itf,r,sig,d1,d2)
    integer, intent(in) :: itf
    type(MatterPowerData) :: CAMB_Pk
    real(dl) sum1,sum2,sum3,t,y,x,w1,w2,w3
    real(dl) x2,rk, fac,r, sig, d1,d2, anorm
    integer i,nint

    nint=3000
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    anorm = 1/(2*const_pi**2)
    do i=1,nint
        t=(i-0.5_dl)/nint
        y=-1.d0+1.d0/t
        rk=y
        d2=MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3*anorm)
        x=y*r
        x2=x*x
        w1=exp(-x2)
        w2=2*x2*w1
        w3=4*x2*(1-x2)*w1
        fac=d2/y/t/t
        sum1=sum1+w1*fac
        sum2=sum2+w2*fac
        sum3=sum3+w3*fac
    enddo
    sum1=sum1/nint
    sum2=sum2/nint
    sum3=sum3/nint
    sig=sqrt(sum1)
    d1=-sum2/sum1
    d2=-sum2*sum2/sum1/sum1 - sum3/sum1

    end subroutine wint

    !!JD 08/13 generalize to variable w

    function omega_m(aa,om_m0,om_v0,wval,waval)
    real(dl) omega_m,omega_t,om_m0,om_v0,aa,wval,waval,Qa2
    Qa2= aa**(-1.0-3.0*(wval+waval))*dexp(-3.0*(1-aa)*waval)
    omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*Qa2+om_m0/aa)
    omega_m=omega_t*om_m0/(om_m0+om_v0*aa*Qa2)
    end function omega_m

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! evolution of omega lambda with expansion factor

    function omega_v(aa,om_m0,om_v0,wval,waval)
    real(dl) aa,omega_v,om_m0,om_v0,omega_t,wval,waval,Qa2
    Qa2= aa**(-1.0-3.0*(wval+waval))*dexp(-3.0*(1-aa)*waval)
    omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*Qa2+om_m0/aa)
    omega_v=omega_t*om_v0*Qa2/(om_v0*Qa2+om_m0/aa)
    end function omega_v

    !!JD end generalize to variable w

    !!AM Below is for HMcode
    SUBROUTINE HMcode(this,State,CAMB_Pk)
    !!AM - A CAMB derived type that I need
    class(THalofit) :: this
    Class(CAMBdata) :: State
    TYPE(MatterPowerData) :: CAMB_Pk
    REAL(dl) :: z, k
    REAL(dl) :: p1h, p2h, pfull, plin
    INTEGER :: i, j, nk, nz
    TYPE(HM_cosmology) :: cosi
    TYPE(HM_tables) :: lut
    REAL(dl), PARAMETER :: pi=3.141592654

    !HMcode developed by Alexander Mead (alexander.j.mead@googlemail.com)
    !Please contact me if you have any questions whatsoever
    !If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
    !If you use the extensions (w(a) and massive neutrinos) then please cite: http://arxiv.org/abs/1602.02154
    !Also consider citing the source code at ASCL: http://ascl.net/1508.001

    !Use imead to switch between the standard and accurate halo-model calcuation
    !0 - Standard (this is just a vanilla halo model calculation with no accuracy tweaks)
    !1 - Accurate from Mead et al. (2015; arXiv 1505.07833)
    !2 - Accurate from Mead et al. (2015; arXiv 1505.07833; no calibration for massive neutrinos)
    IF(this%halofit_version==halofit_halomodel) this%imead=0
    IF(this%halofit_version==halofit_mead) this%imead=1
    IF(this%halofit_version==halofit_mead2015) this%imead=2

    HM_verbose = (FeedbackLevel>1)

    IF(HM_verbose) WRITE(*,*)
    IF(HM_verbose) WRITE(*,*) 'HMcode: Running HMcode'
    IF(HM_verbose) WRITE(*,*)

    !!AM - Translate from CAMB variables to my variables
    nz=CAMB_PK%num_z
    nk=CAMB_PK%num_k

    !!TT - Assign baryon parameters
    cosi%A_baryon = this%HMCode_A_baryon
    cosi%eta_baryon = this%HMCode_eta_baryon

    !!AM - Assign cosmological parameters for the halo model calculation
    CALL assign_HM_cosmology(State,cosi)

    !Fill growth function table (only needs to be done once)
    CALL fill_growtab(cosi)

    !Loop over redshifts
    DO j=1,nz

        !Initialise the specific HM_cosmology (fill sigma(R) and P_lin HM_tables)
        !Currently this needs to be done at each z (mainly because of scale-dependent growth with neutrinos)
        !For non-massive-neutrino models this could only be done once, which would speed things up a bit
        CALL initialise_HM_cosmology(j,cosi,CAMB_PK)

        !Sets the current redshift from the table
        z=CAMB_Pk%Redshifts(j)

        !Initiliasation for the halomodel calcualtion (needs to be done for each z)
        CALL this%halomod_init(z,lut,cosi)

        !Loop over k values and calculate P(k)
        !$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin, pfull,p1h,p2h)
        DO i=1,nk
            k=exp(CAMB_Pk%log_kh(i))
            plin=p_lin(k,z,0,cosi)
            CALL this%halomod(k,z,p1h,p2h,pfull,plin,lut,cosi)
            CAMB_Pk%nonlin_ratio(i,j)=sqrt(pfull/plin)
        END DO
        !$OMP END PARALLEL DO

    END DO

    END SUBROUTINE HMcode

    FUNCTION Delta_v(this,z,lut,cosm)
    class(THalofit) :: this
    !Function for the virialised overdensity
    REAL(dl) :: Delta_v
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0) THEN
        !Value that is normally used in halo model predictions
        Delta_v=200.
    ELSE IF(this%imead==1 .OR. this%imead==2) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        Delta_v=418.*(Omega_m_hm(z,cosm)**(-0.352))
        !Mead et al. (2016; arXiv 1602.02154) neutrino addition
        IF(this%imead==1) Delta_v=Delta_v*(1.+0.916*cosm%f_nu)
    END IF

    END FUNCTION Delta_v

    FUNCTION delta_c(this,z,lut,cosm)
    class(THalofit) :: this
    !Function for the linear collapse density
    REAL(dl) :: delta_c
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0) THEN
        delta_c=1.686
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        delta_c=1.59+0.0314*log(lut%sig8z)
        IF(this%imead==1) THEN
            delta_c=delta_c*(1.+0.262*cosm%f_nu) !Mead et al. (2016; arXiv 1602.02154) neutrino addition
            delta_c=delta_c*(1.+0.0123*log10(Omega_m_hm(z,cosm))) !Nakamura & Suto (1997) fitting formula for LCDM
        END IF
    END IF

    END FUNCTION delta_c

    FUNCTION eta(this,z,lut,cosm)
    class(THalofit) :: this
    !Function eta that puffs halo profiles
    REAL(dl) :: eta
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL(dl) :: eta0

    IF(this%imead==0) THEN
        eta=0.
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
        !eta=0.603-0.3*lut%sig8z
        !AM - made baryon feedback parameter obvious
        eta0=cosm%eta_baryon
        !eta0=1.03-0.11*cosm%A_baryon !This is the original one-parameter relation from 1505.07833
        !eta0=0.98-0.12*cosm%A_baryon !This is an updated one-parameter relation: Section 4.1.2 of 1707.06627
        eta=eta0-0.3*lut%sig8z
    END IF

    END FUNCTION eta

    FUNCTION kstar(this,z,lut,cosm)
    class(THalofit) :: this
    !Function k* that cuts off the 1-halo term at large scales
    REAL(dl) :: kstar
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0) THEN
        !Set to zero for the standard Poisson one-halo term
        kstar=0.
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !One-halo cut-off wavenumber
        !Mead et al. (2015; arXiv 1505.07833) value
        kstar=0.584*(lut%sigv)**(-1.)
    END IF

    END FUNCTION kstar

    FUNCTION As(this,z,lut,cosm)
    class(THalofit) :: this
    !Halo concentration pre-factor from Bullock et al. (2001) relation
    REAL(dl) :: As
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0) THEN
        !Set to 4 for the standard Bullock value
        As=4.
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
        !As=3.13
        !AM - added for easy modification of feedback parameter
        As=cosm%A_baryon
    END IF

    END FUNCTION As

    FUNCTION fdamp(this,z,lut,cosm)
    class(THalofit) :: this
    !Linear power damping function from Mead et al. (2015; arXiv 1505.07833)
    REAL(dl) ::fdamp
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    !Linear theory damping factor
    IF(this%imead==0) THEN
        !Set to 0 for the standard linear theory two halo term
        fdamp=0.
    ELSE IF(this%imead==1) THEN
        !Mead et al. (2016; arXiv 1602.02154) value
        fdamp=0.0095*lut%sigv100**1.37
    ELSE IF(this%imead==2) THEN
        !Mead et al. (2015) value
        fdamp=0.188*lut%sig8z**4.29
    END IF

    !Catches extreme values of fdamp
    IF(fdamp<1.e-3) fdamp=1.e-3
    IF(fdamp>0.99)  fdamp=0.99

    END FUNCTION fdamp

    FUNCTION alpha(this,z,lut,cosm)
    class(THalofit) :: this
    !Two- to one-halo transition smoothing from Mead et al. (2015; arXiv 1505.07833)
    REAL(dl) :: alpha
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(this%imead==0) THEN
        !Set to 1 for the standard halomodel sum of one- and two-halo terms
        alpha=1.
    ELSE IF(this%imead==1) THEN
        !This uses the top-hat defined neff (HALOFIT uses Gaussian filtered fields instead)
        !Mead et al. (2016; arXiv 1602.02154) value
        alpha=3.24*1.85**lut%neff
    ELSE IF(this%imead==2) THEN
        !Mead et al. (2015) value
        alpha=2.93*1.77**lut%neff
    END IF

    !Catches values of alpha that are crazy
    IF(alpha>2.)  alpha=2.
    IF(alpha<0.5) alpha=0.5

    END FUNCTION alpha

    FUNCTION r_nl(lut)
    !Calculates R_nl, defined by nu(R_nl)=1., nu=dc/sigma(R)
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL(dl) :: r_nl

    IF(lut%nu(1)>1.) THEN
        !This catches some very strange values that appear for odd cosmological models
        r_nl=lut%rr(1)
    ELSE
        r_nl=exp(find(log(1.d0),log(lut%nu),log(lut%rr),lut%n,3,3,2))
    END IF

    END FUNCTION r_nl

    SUBROUTINE halomod(this,k,z,p1h,p2h,pfull,plin,lut,cosm)
    class(THalofit) :: this
    !Calcuates 1-halo and 2-halo terms and combines them to form the full halomodel power
    REAL(dl), INTENT(OUT) :: p1h, p2h, pfull
    REAL(dl), INTENT(IN) :: plin, k, z
    REAL(dl) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut
    LOGICAL, PARAMETER :: sigma=.FALSE.

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN
        p1h=0.
        p2h=0.
    ELSE
        p1h=this%p_1h(k,z,lut,cosm)
        p2h=this%p_2h(k,z,plin,lut,cosm)
    END IF

    a=this%alpha(z,lut,cosm)
    pfull=(p2h**a+p1h**a)**(1./a)
    IF(sigma) THEN
        pfull=sigmac(k,z,cosm)
    END IF

    END SUBROUTINE halomod

    SUBROUTINE fill_table(min,max,arr,n)
    !Fills array 'arr' in equally spaced intervals
    INTEGER :: i
    REAL(dl), INTENT(IN) :: min, max
    REAL(dl), ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: n

    !Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
        arr(1)=min
    ELSE IF(n>1) THEN
        DO i=1,n
            arr(i)=min+(max-min)*float(i-1)/float(n-1)
        END DO
    END IF

    END SUBROUTINE fill_table

    ! SUBROUTINE fill_table8(min,max,arr,n)

    ! !Fills array 'arr' in equally spaced intervals
    ! IMPLICIT NONE
    ! INTEGER :: i
    ! REAL(dl), INTENT(IN) :: min, max
    ! real(dl), ALLOCATABLE :: arr(:)
    ! INTEGER, INTENT(IN) :: n

    ! !Allocate the array, and deallocate it if it is full
    ! IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ! ALLOCATE(arr(n))
    ! arr=0.

    ! IF(n==1) THEN
    !     arr(1)=min
    ! ELSE IF(n>1) THEN
    !     DO i=1,n
    !         arr(i)=min+(max-min)*float(i-1)/float(n-1)
    !     END DO
    ! END IF

    ! END SUBROUTINE fill_table8

    SUBROUTINE fill_plintab(iz,cosm,CAMB_PK)

    !Fills internal HMcode HM_tables for the linear power spectrum at z=0
    TYPE(MatterPowerData), INTENT(IN) :: CAMB_PK
    INTEGER, INTENT(IN) :: iz
    TYPE(HM_cosmology) :: cosm
    INTEGER :: i
    REAL(dl) :: z, g
    REAL(dl), ALLOCATABLE :: k(:), Pk(:), Pkc(:)
    INTEGER, PARAMETER :: imeth=imeth_pk
    REAL(dl), PARAMETER :: pi=3.141592654
    REAL(dl), PARAMETER :: kmin=kmin_pk
    REAL(dl), PARAMETER :: kmax=kmax_pk
    INTEGER :: nk=nk_pk

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Filling linear power HM_tables'

    !Fill arrays
    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
    IF(ALLOCATED(cosm%log_plin))   DEALLOCATE(cosm%log_plin)
    IF(ALLOCATED(cosm%log_plinc))  DEALLOCATE(cosm%log_plinc)

    IF(imeth==1) THEN

        !Fill k-table with the same k points as in the CAMB calculation
        !If a user has specified lots of points this could make the halo-model
        !calculation chug
        nk=CAMB_PK%num_k
        cosm%nk=nk
        ALLOCATE(cosm%log_k_plin(nk))
        cosm%log_k_plin=CAMB_Pk%log_kh

    ELSE IF(imeth==2) THEN

        !Fill a k-table with an equal-log-spaced k range
        !Note that the minimum should be such that the linear spectrum is accurately a power-law below this wavenumber
        cosm%nk=nk
        CALL fill_table(log(kmin),log(kmax),cosm%log_k_plin,nk)

    ELSE

      STOP 'FILL_PLINTAB: Error, imeth specified incorrectly'

    END IF

    ALLOCATE(k(nk))
    k=exp(cosm%log_k_plin)

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: k_min:', k(1)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: k_max:', k(nk)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: nk:', nk

    ALLOCATE(Pk(nk),Pkc(nk))

    !Find the redshift
    z=CAMB_Pk%Redshifts(iz)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: z of input:', z

    !Fill power table, both cold- and all-matter
    DO i=1,nk
        !Take the power from the current redshift choice
        Pk(i)=MatterPowerData_k(CAMB_PK,k(i),iz)*(k(i)**3/(2.*pi**2))
        Pkc(i)=Pk(i)*Tcb_Tcbnu_ratio(k(i),z,cosm)**2
    END DO

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Delta2_min:', Pk(1)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Delta2_max:', Pk(nk)

    !Calculate the growth factor at the redshift of interest
    g=grow(z,cosm)

    ALLOCATE(cosm%log_plin(nk),cosm%log_plinc(nk))

    !Grow the power to z=0
    cosm%log_plin=log(Pk/(g**2.))
    cosm%log_plinc=log(Pkc/(g**2.))

    DEALLOCATE(k,Pk,Pkc)

    !Check sigma_8 value
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: sigma_8:', sigma(8.d0,0.d0,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Done'
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE fill_plintab

    FUNCTION Tcb_Tcbnu_ratio(k,z,cosm)

    !Calculates the ratio of T(k) for cold vs. all matter
    !Uses approximations in Eisenstein & Hu (1999; arXiv 9710252)
    !Note that this assumes that there are exactly 3 species of neutrinos with
    !Nnu<=3 of these being massive, and with the mass split evenly between the number of massive species.

    REAL(dl) :: Tcb_Tcbnu_ratio
    REAL(dl), INTENT(IN) :: k, z
    REAL(dl) :: D, Dcb, Dcbnu, pcb, zeq, q, yfs
    REAL(dl) :: BigT
    TYPE(HM_cosmology) :: cosm

    IF(cosm%f_nu==0.) THEN

        Tcb_Tcbnu_ratio=1.

    ELSE

        !Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
        pcb=(5.-sqrt(1.+24.*(1.-cosm%f_nu)))/4.

        !Theta for temperature (BigT=T/2.7 K)
        BigT=cosm%Tcmb/2.7

        !The matter-radiation equality redshift
        zeq=(2.5e4)*cosm%om_m*(cosm%h**2.)*(BigT**(-4.))

        !The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m \approx 1)
        !For my purpose (just the ratio) seems to work better using the EdS growth function result, \propto a .
        !In any case, can't use grow at the moment because that is normalised by default.
        D=(1.+zeq)/(1.+z)

        !Wave number relative to the horizon scale at equality (equation 5)
        !Extra factor of h becauase all my k are in units of h/Mpc
        q=k*cosm%h*BigT**2./(cosm%om_m*cosm%h**2.)

        !Free streaming scale (equation 14)
        !Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
        !with Nnu of these being massve with the mass split evenly between Nnu species.
        yfs=17.2*cosm%f_nu*(1.+0.488*cosm%f_nu**(-7./6.))*(cosm%Nnu*q/cosm%f_nu)**2.

        !These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
        !Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
        !Equations (12) and (13)
        Dcb=(1.+(D/(1.+yfs))**0.7)**(pcb/0.7)
        Dcbnu=((1.-cosm%f_nu)**(0.7/pcb)+(D/(1.+yfs))**0.7)**(pcb/0.7)

        Tcb_Tcbnu_ratio=Dcb/Dcbnu

    END IF

    END FUNCTION Tcb_Tcbnu_ratio

    SUBROUTINE assign_HM_cosmology(State,cosm)
    class(CAMBdata) :: State
    !Assigns the internal HMcode cosmological parameters
    TYPE(HM_cosmology) :: cosm
    real(dl) h2

    associate(CP => State%CP)
        !Converts CAMB parameters to Meadfit parameters
        h2 = (CP%H0/100)**2
        cosm%om_m=(CP%omch2+CP%ombh2+CP%omnuh2)/h2
        cosm%om_c=CP%omch2/h2
        cosm%om_b=CP%ombh2/h2
        cosm%om_v=State%omega_de
        call CP%DarkEnergy%Effective_w_wa(cosm%w, cosm%wa)
        cosm%f_nu=CP%omnuh2/h2/cosm%om_m
        cosm%h=CP%H0/100
        cosm%Tcmb=CP%tcmb
        cosm%Nnu=CP%Num_Nu_massive

        !n_s is read in here. The non-linear CAMB module does not work if there is more than
        !one value in this array, so explicitly setting '1' here is fine.
        cosm%ns= CP%InitPower%Effective_ns()
    end associate

    !Write out cosmological parameters if necessary
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_m:', cosm%om_m
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_c:', cosm%om_c
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_b:', cosm%om_b
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_v:', cosm%om_v
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: w_0:', cosm%w
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: w_a:', cosm%wa
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: f_nu:', cosm%f_nu
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: n_s:', cosm%ns
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: h:', cosm%h
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: T_CMB [K]:', cosm%Tcmb
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: N_nu (massive):', cosm%Nnu
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: A_baryon:', cosm%A_baryon
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: eta_baryon:', cosm%eta_baryon
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE assign_HM_cosmology

    SUBROUTINE initialise_HM_cosmology(iz,cosm,CAMB_PK)

    !Sets up HM_tables of sigma, growth and linear power for the HM_cosmology
    TYPE(MatterPowerData), INTENT(IN) :: CAMB_PK
    TYPE(HM_cosmology) :: cosm
    INTEGER, INTENT(IN) :: iz

    !Fill linear power table and grows it to z=0
    CALL fill_plintab(iz,cosm,CAMB_PK)

    !Fill sigma(r) table
    CALL fill_sigtab(cosm)

    END SUBROUTINE initialise_HM_cosmology

    SUBROUTINE allocate_LUT(lut)

    !Allocates memory for the HMcode look-up HM_tables
    TYPE(HM_tables) :: lut
    INTEGER :: n

    n=lut%n
    ALLOCATE(lut%zc(n),lut%m(n),lut%c(n),lut%rv(n))
    ALLOCATE(lut%nu(n),lut%rr(n),lut%sigf(n),lut%sig(n))

    lut%zc=0.
    lut%m=0.
    lut%c=0.
    lut%rv=0.
    lut%nu=0.
    lut%rr=0.
    lut%sigf=0.
    lut%sig=0.

    END SUBROUTINE allocate_LUT

    SUBROUTINE deallocate_LUT(lut)

    !Deallocates HMcode look-up HM_tables
    TYPE(HM_tables) :: lut

    DEALLOCATE(lut%zc,lut%m,lut%c,lut%rv,lut%nu,lut%rr,lut%sigf,lut%sig)

    END SUBROUTINE deallocate_LUT

    SUBROUTINE halomod_init(this,z,lut,cosm)
    class(THalofit) :: this
    !Halo-model initialisation routine
    !Computes look-up HM_tables necessary for the halo model calculations
    REAL(dl), INTENT(IN) :: z
    INTEGER :: i
    REAL(dl) :: Dv, dc, f, m, nu, r, sig
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut
    REAL(dl), PARAMETER :: mmin=mmin_HMcode !Lower mass limit
    REAL(dl), PARAMETER :: mmax=mmax_HMcode !Upper mass limit
    INTEGER, PARAMETER :: n=n_HMcode !Number of entries in look-up HM_tables.

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: Filling look-up HM_tables'
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: HM_tables being filled at redshift:', z
    lut%z=z

    !Find value of sigma_v, sig8, etc.

    lut%sigv=sigmaV(0.d0,z,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
    lut%sigv100=sigmaV(100.d0,z,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
    lut%sig8z=sigma(8.d0,z,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sig8(z):', lut%sig8z

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    lut%n=n
    CALL allocate_lut(lut)

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M_min:', mmin
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M_max:', mmax

    dc=this%delta_c(z,lut,cosm)
    lut%dc=dc

    !$OMP PARALLEL DO default(shared), private(m,r,sig,nu)
    DO i=1,n

        m=exp(log(mmin)+log(mmax/mmin)*float(i-1)/float(n-1))
        r=radius_m(m,cosm)
        sig=sigmac(r,z,cosm)
        nu=dc/sig

        lut%m(i)=m
        lut%rr(i)=r
        lut%sig(i)=sig
        lut%nu(i)=nu

    END DO
    !$OMP END PARALLEL DO

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: m, r, nu, sig HM_tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
        lut%sigf(i)=sigmac(lut%rr(i)*f,z,cosm)
    END DO
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sigf HM_tables filled'

    !Fill virial radius table using real radius table
    Dv=this%Delta_v(z,lut,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: rv HM_tables filled'
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M min [Msun/h]:', lut%m(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M max [Msun/h]:', lut%m(lut%n)

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1./lut%rnl

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    !Calcuate the effective spectral index at the collapse scale
    lut%neff=neff(lut,cosm)

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: n_eff:', lut%neff

    !Get the concentration for all the haloes
    CALL this%conc_bull(z,lut,cosm)

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: c HM_tables filled'
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: c min [Msun/h]:', lut%c(lut%n)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: c max [Msun/h]:', lut%c(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: Done'
    IF(HM_verbose) WRITE(*,*)
    IF(HM_verbose) CALL this%write_parameters(z,lut,cosm)

    !Switch off verbose mode if doing multiple z
    HM_verbose= .false.

    END SUBROUTINE halomod_init

    SUBROUTINE write_parameters(this,z,lut,cosm)
    class(THalofit) :: this
    !This subroutine writes out the halomodel parameters at the current redshift
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(HM_verbose) WRITE(*,*) 'WRITE_PARAMETERS: at this redshift'
    IF(HM_verbose) WRITE(*,*) '=================================='
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'z:', z
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'Dv:', this%Delta_v(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'dc:', this%delta_c(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'eta:', this%eta(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'k*:', this%kstar(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'A:', this%As(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'fdamp:', this%fdamp(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'alpha:', this%alpha(z,lut,cosm)
    IF(HM_verbose) WRITE(*,*) '=================================='
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE write_parameters

    PURE FUNCTION radius_m(m,cosm)

    !Calculates the co-moving radius that encloses a mass 'm' in the homogeneous Universe
    REAL(dl) :: radius_m
    REAL(dl), INTENT(IN) :: m
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: pi=3.141592654

    radius_m=(3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

    END FUNCTION radius_m

    FUNCTION neff(lut,cosm)

    !Finds the effective spectral index at the collapse scale r_nl, where nu(r_nl)=1.
    REAL(dl) :: neff
    REAL(dl) :: ns
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut 
    REAL(dl), PARAMETER :: tmin=0.
    REAL(dl), PARAMETER :: tmax=1.
    INTEGER, PARAMETER :: itype=1 ! Cold matter here
    REAL(dl), PARAMETER :: acc=acc_neff
    INTEGER, PARAMETER :: iorder=iorder_neff

    !neff=-3.-2.*neff_integral(lut%rnl, lut%z, cosm)/lut%dc**2
    neff=-3.-2.*integrate(tmin,tmax,neff_integrand,lut%rnl,lut%z,itype,cosm,acc,iorder)/lut%dc**2

    !For some bizarre cosmological models r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    ns=cosm%ns
    IF(neff<ns-4.) neff=ns-4.
    IF(neff>ns)    neff=ns

    END FUNCTION neff

    SUBROUTINE conc_bull(this,z,lut,cosm)
    class(THalofit) :: this
    !Calculates the Bullock et al. (2001) concentration-mass relation
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology) :: cosm, cos_lcdm
    TYPE(HM_tables) :: lut
    REAL(dl) :: A, zf, ainf, zinf, g_lcdm, g_wcdm, pow
    INTEGER :: i

    !Amplitude of relation (4. in Bullock et al. 2001)
    A=this%As(z,lut,cosm)

    !Fill the collapse time look-up table
    CALL this%zcoll_bull(z,cosm,lut)

    !Fill the concentration look-up table
    DO i=1,lut%n
        zf=lut%zc(i)
        lut%c(i)=A*(1.+zf)/(1.+z)
    END DO

    !Dolag (2004) prescription for adding DE dependence

    !This is approximately z=infinity
    zinf=10.
    g_wcdm=grow(zinf,cosm)

    !Make a LCDM HM_cosmology
    !Only need to make sure model is flat with the same Omega_m and w=-1
    !This is *only* used for a calculation of the growth function
    cos_lcdm=cosm
    DEALLOCATE(cos_lcdm%growth)
    DEALLOCATE(cos_lcdm%a_growth)
    cos_lcdm%w=-1.
    cos_lcdm%wa=0.
    cos_lcdm%om_v=1.-cosm%om_m !Added this so that 'making a LCDM cosmology' works for curved models.

    ainf=1./(1.+zinf)

    !Needs to use grow_int explicitly here for LCDM model to avoid growth HM_tables
    g_lcdm=growint(ainf,cos_lcdm)

    !This is the Dolag et al. (2004) correction for halo concentrations
    IF(this%imead==0 .OR. this%imead==2) THEN
        ! Mead et al. (2015) used the Dolag (2004) correction
        pow=1.
    ELSE IF(this%imead==1) THEN
        ! Mead et al. (2016) changed the power to 1.5 to better accomodate more extreme dark-energy models
        pow=1.5
    END IF
    lut%c=lut%c*((g_wcdm/g_lcdm)**pow)

    END SUBROUTINE conc_bull

    FUNCTION growint(a,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL(dl) :: growint
    REAL(dl), INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: b
    INTEGER :: i, j
    INTEGER :: n
    REAL(dl) :: x, dx
    REAL(dl) :: f1, f2, fx
    REAL(dl) :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    real(dl), PARAMETER :: acc=1d-3
    INTEGER, PARAMETER :: iorder=3

    !Integration range for integration parameter
    !Note a -> 1
    b=1.d0

    IF(a==b) THEN

        !Fix the answer to zero if the integration limits are identical
        growint=exp(0.d0)

    ELSE

        !Reset the sum variable for the integration
        sum_2n=0.d0
        sum_n=0.d0
        sum_old=0.d0
        sum_new=0.d0

        DO j=1,jmax

            !Note, you need this to be 1+2**n for some integer n
            !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n=1+2**(j-1)

            !Calculate the dx interval for this value of 'n'
            dx=(b-a)/REAL(n-1)

            IF(j==1) THEN

                !The first go is just the trapezium of the end points
                f1=growint_integrand(a,cosm)
                f2=growint_integrand(b,cosm)
                sum_2n=0.5d0*(f1+f2)*dx
                sum_new=sum_2n

            ELSE

                !Loop over only new even points to add these to the integral
                DO i=2,n,2
                    x=a+(b-a)*REAL(i-1)/REAL(n-1)
                    fx=growint_integrand(x,cosm)
                    sum_2n=sum_2n+fx
                END DO

                !Now create the total using the old and new parts
                sum_2n=sum_n/2.d0+sum_2n*dx

                !Now calculate the new sum depending on the integration order
                IF(iorder==1) THEN
                    sum_new=sum_2n
                ELSE IF(iorder==3) THEN
                    sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
                ELSE
                    STOP 'GROWINT: Error, iorder specified incorrectly'
                END IF

            END IF

            IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
                !jmin avoids spurious early convergence
                growint=exp(sum_new)
                EXIT
            ELSE IF(j==jmax) THEN
                STOP 'GROWINT: Integration timed out'
            ELSE
                !Integral has not converged so store old sums and reset sum variables
                sum_old=sum_new
                sum_n=sum_2n
                sum_2n=0.d0
            END IF

        END DO

    END IF

    END FUNCTION growint

    FUNCTION growint_integrand(a,cosm)

    !Integrand for the approximate growth integral
    IMPLICIT NONE
    REAL(dl) :: growint_integrand
    REAL(dl), INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: gam

    IF(cosm%w<-1.) THEN
        gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
        gam=0.55+0.05*(1.+cosm%w)
    ELSE
        gam=0.55
    END IF

    !Note the minus sign here
    !AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    growint_integrand=-(Omega_cold_hm(-1.+1./a,cosm)**gam)/a

    END FUNCTION growint_integrand

    SUBROUTINE zcoll_bull(this,z,cosm,lut)
    class(THalofit) :: this
    !Calcuates the halo collapse redshift according to the Bullock (2001) prescription
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut
    REAL(dl) :: dc
    REAL(dl) :: af, zf, RHS, a, growz
    INTEGER :: i

    !This fills up the halo formation redshift table as per Bullock relations

    !Needs to interpolate g(z) which should be pretty linear for a<0.05
    !in 'g(a) vs a' space for all standard cosmologies

    dc=this%delta_c(z,lut,cosm)

    !Find the growth function at the current redshift
    a=1./(1.+z)
    growz=find(a,cosm%a_growth,cosm%growth,cosm%ng,3,3,2)

    !Do numerical inversion
    DO i=1,lut%n

        RHS=dc*growz/lut%sigf(i)

        IF(RHS>growz) THEN
            !This is the case of 'halo forms in the future'
            !in this case set formation redshift to current redshift
            zf=z
        ELSE
            af=find(RHS,cosm%growth,cosm%a_growth,cosm%ng,3,3,2)
            zf=-1.+1./af
        END IF

        lut%zc(i)=zf

    END DO

    END SUBROUTINE zcoll_bull

    FUNCTION mass_r(r,cosm)

    !Calcuates the average mass enclosed at co-moving radius r
    REAL(dl) :: mass_r, r
    TYPE(HM_cosmology) :: cosm
    REAL(dl), PARAMETER :: pi=3.141592654

    !Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*cosmic_density(cosm)*(r**3.)

    END FUNCTION mass_r

    PURE FUNCTION cosmic_density(cosm)

    !The z=0 cosmological matter density
    REAL(dl) :: cosmic_density
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !In M_sun per Mpc^3 with h factors included. The constant does this.
    cosmic_density=(2.775e11)*cosm%om_m

    END FUNCTION cosmic_density

    FUNCTION find_pk(k,itype,cosm)

    !Look-up and interpolation for P(k,z=0)
    REAL(dl) :: find_pk
    REAL(dl) :: kmax, ns
    REAL(dl), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: itype
    INTEGER :: n
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Set number of k points as well as min and max k values
    !Note that the min k value should be set to the same as the CAMB min k value
    n=SIZE(cosm%log_k_plin)
    kmax=exp(cosm%log_k_plin(n))

    !Spectral index used in the high-k extrapolation
    ns=cosm%ns

    IF(k>kmax) THEN
        !Do some interpolation here based on knowledge of things at high k
        IF(itype==0) THEN
            find_pk=exp(cosm%log_plin(n))*((log(k)/log(kmax))**2.)*((k/kmax)**(ns-1.))
        ELSE IF(itype==1) THEN
            find_pk=exp(cosm%log_plinc(n))*((log(k)/log(kmax))**2.)*((k/kmax)**(ns-1.))
        END IF
    ELSE
        !Otherwise use the standard find algorithm
        IF(itype==0) THEN
            find_pk=exp(find(log(k),cosm%log_k_plin,cosm%log_plin,cosm%nk,3,3,2))
        ELSE IF(itype==1) THEN
            find_pk=exp(find(log(k),cosm%log_k_plin,cosm%log_plinc,cosm%nk,3,3,2))
        END IF
    END IF

    END FUNCTION find_pk

    FUNCTION p_lin(k,z,itype,cosm)

    !Looks up the value for the linear power spectrum
    REAL(dl) :: p_lin
    REAL(dl), INTENT (IN) :: k, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised to z=0

    p_lin=(grow(z,cosm)**2.)*find_pk(k,itype,cosm)

    END FUNCTION p_lin

    FUNCTION p_2h(this,k,z,plin,lut,cosm)
    class(THalofit) :: this
    !Calculates the 2-halo term
    REAL(dl) :: p_2h
    REAL(dl), INTENT(IN) :: k, plin
    REAL(dl) :: sigv, frac
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Damping function
    frac=this%fdamp(z,lut,cosm)

    IF(this%imead==0 .OR. frac<1.e-3) THEN
        p_2h=plin
    ELSE
        sigv=lut%sigv
        p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2.)
    END IF

    !For some strange cosmologies frac>1. so this must be added to prevent p_2h<0.
    IF(p_2h<0.) p_2h=0.

    END FUNCTION p_2h

    FUNCTION p_1h(this,k,z,lut,cosm)
    class(THalofit) :: this
    !Calculates the 1-halo term
    REAL(dl) :: p_1h
    REAL(dl), INTENT(IN) :: k, z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: Dv, g, fac, et, ks, wk
    REAL(dl) :: integrand(lut%n)
    REAL(dl) :: sum
    INTEGER :: i
    REAL(dl), PARAMETER :: pi=3.141592654

    !Does the one-halo power integral

    !Only call eta once
    et=this%eta(z,lut,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
        g=gnu(lut%nu(i))
        wk=win(k*(lut%nu(i)**et),lut%rv(i),lut%c(i))
        integrand(i)=(lut%rv(i)**3.)*g*(wk**2.)
    END DO

    !Carries out the integration
    sum=inttab(lut%nu,integrand,lut%n,1)

    !Virial density
    Dv=this%Delta_v(z,lut,cosm)

    !These are just numerical factors from the 1-halo integral in terms of nu!
    p_1h=sum*2.*Dv*(k**3.)/(3.*pi)

    !Damping of the 1-halo term at very large scales
    ks=this%kstar(z,lut,cosm)

    !Prevents problems if k/ks is very large
    IF(ks==0.) THEN
        fac=0.
    ELSE IF((k/ks)**2.>7.) THEN
        fac=0.
    ELSE
        fac=exp(-((k/ks)**2.))
    END IF

    !Damping of the one-halo term at very large scales
    p_1h=p_1h*(1.-fac)

    END FUNCTION p_1h

    SUBROUTINE fill_sigtab(cosm)

    !Fills look-up HM_tables for sigma(R)
    INTEGER :: i
    TYPE(HM_cosmology) :: cosm
    REAL(dl), ALLOCATABLE :: r(:), sig(:)
    REAL(dl), PARAMETER :: rmin=rmin_sigma
    REAL(dl), PARAMETER :: rmax=rmax_sigma
    INTEGER, PARAMETER :: nsig=n_sigma

    !This fills up HM_tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigmac(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is approximately power-law below and above these values of R
    !This wouldn't be appropriate for models with a small-scale linear spectrum cut-off (e.g., WDM)

    !Allocate arrays
    IF(ALLOCATED(cosm%log_r_sigma)) DEALLOCATE(cosm%log_r_sigma)
    IF(ALLOCATED(cosm%log_sigma))   DEALLOCATE(cosm%log_sigma)

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    ALLOCATE(r(nsig),sig(nsig))

    IF(HM_verbose) WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
    IF(HM_verbose) WRITE(*,*) 'SIGTAB: R_min:', rmin
    IF(HM_verbose) WRITE(*,*) 'SIGTAB: R_max:', rmax
    IF(HM_verbose) WRITE(*,*) 'SIGTAB: Values:', nsig

    !$OMP PARALLEL DO default(shared)
    DO i=1,nsig

        !Equally spaced r in log
        r(i)=exp(log(rmin)+log(rmax/rmin)*float(i-1)/float(nsig-1))

        sig(i)=sigma(r(i),0.d0,1,cosm)

    END DO
    !$OMP END PARALLEL DO

    IF(HM_verbose) WRITE(*,*) 'SIGTAB: sigma_min:', sig(nsig)
    IF(HM_verbose) WRITE(*,*) 'SIGTAB: sigma_max:', sig(1)

    cosm%nsig=nsig
    ALLOCATE(cosm%log_r_sigma(nsig),cosm%log_sigma(nsig))
    cosm%log_r_sigma=log(r)
    cosm%log_sigma=log(sig)

    DEALLOCATE(r,sig)

    IF(HM_verbose) WRITE(*,*) 'SIGTAB: Done'
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE fill_sigtab

    FUNCTION sigmac(r,z,cosm)

    !Finds sigma_cold(R) from look-up table
    REAL(dl) :: sigmac
    REAL(dl), INTENT(IN) :: r, z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Assumes scale-independet growth for the cold matter
    !Uses the approximation sigma(R,z)=g(z)*sigma(R,z=0)

    sigmac=grow(z,cosm)*exp(find(log(r),cosm%log_r_sigma,cosm%log_sigma,cosm%nsig,3,3,2))

    END FUNCTION sigmac

    FUNCTION wk_tophat(x)

    !The normlaised Fourier Transform of a top-hat
    REAL(dl) :: wk_tophat
    REAL(dl), INTENT(IN) :: x
    REAL(dl), PARAMETER :: dx=1e-3 ! Taylor expansion for |x|<dx

    !Taylor expansion used for low x to avoid cancellation problems
    IF(abs(x)<dx) THEN
        wk_tophat=1.-(x**2.)/10.
    ELSE
        wk_tophat=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

    END FUNCTION wk_tophat

    FUNCTION wk_tophat_deriv(x)

    ! The derivative of a normlaised Fourier Transform of a spherical top-hat
    IMPLICIT NONE
    REAL(dl) :: wk_tophat_deriv
    REAL(dl), INTENT(IN) :: x
    REAL(dl), PARAMETER :: dx=1e-3 ! Taylor expansion for |x|<dx

    ! Taylor expansion used for low x to avoid cancelation problems
    IF (abs(x)<dx) THEN
        wk_tophat_deriv=-x/5.+x**3/70.
    ELSE
        wk_tophat_deriv=(3./x**4)*((x**2-3.)*sin(x)+3.*x*cos(x))
    END IF

    END FUNCTION wk_tophat_deriv

    FUNCTION inttab(x,y,n,iorder)

    !Integrates tables y(x)dx
    REAL(dl) :: inttab
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x(n), y(n)
    REAL(dl) :: a, b, c, d, h
    REAL(dl) :: q1, q2, q3, qi, qf
    REAL(dl) :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    real(dl) :: sum
    INTEGER :: i, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    sum=0.d0

    IF(iorder==1) THEN

        !Sums over all Trapezia (a+b)*h/2
        DO i=1,n-1
            a=y(i+1)
            b=y(i)
            h=x(i+1)-x(i)
            sum=sum+(a+b)*h/2.d0
        END DO

    ELSE IF(iorder==2) THEN

        DO i=1,n-2

            x1=x(i)
            x2=x(i+1)
            x3=x(i+2)

            y1=y(i)
            y2=y(i+1)
            y3=y(i+2)

            CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

            q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
            q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
            q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

            !Takes value for first and last sections but averages over sections where you
            !have two independent estimates of the area
            IF(n==3) THEN
                sum=sum+q3-q1
            ELSE IF(i==1) THEN
                sum=sum+(q2-q1)+(q3-q2)/2.d0
            ELSE IF(i==n-2) THEN
                sum=sum+(q2-q1)/2.d0+(q3-q2)
            ELSE
                sum=sum+(q3-q1)/2.
            END IF

        END DO

    ELSE IF(iorder==3) THEN

        DO i=1,n-1

            !First choose the integers used for defining cubics for each section
            !First and last are different because the section does not lie in the *middle* of a cubic

            IF(i==1) THEN

                i1=1
                i2=2
                i3=3
                i4=4

            ELSE IF(i==n-1) THEN

                i1=n-3
                i2=n-2
                i3=n-1
                i4=n

            ELSE

                i1=i-1
                i2=i
                i3=i+1
                i4=i+2

            END IF

            x1=x(i1)
            x2=x(i2)
            x3=x(i3)
            x4=x(i4)

            y1=y(i1)
            y2=y(i2)
            y3=y(i3)
            y4=y(i4)

            CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

            !These are the limits of the particular section of integral
            xi=x(i)
            xf=x(i+1)

            qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
            qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

            sum=sum+qf-qi

        END DO

    ELSE

        ERROR STOP 'INTTAB: Error, order not specified correctly'

    END IF

    inttab=REAL(sum)

    END FUNCTION inttab

    FUNCTION sigma(r,z,itype,cosm)

    !Gets sigma(R)
    IMPLICIT NONE
    REAL(dl) :: sigma
    REAL(dl), INTENT(IN) :: r, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: tmin=0.
    REAL(dl), PARAMETER :: tmax=1.
    REAL(dl), PARAMETER :: acc=acc_sigma !AM Jul 19: Upgraded accuracy from 1d-3 to 1d-4
    INTEGER, PARAMETER :: iorder=iorder_sigma

    !sigma=sqrt(sigint(r,z,itype,cosm,acc,iorder))
    sigma=sqrt(integrate(tmin,tmax,sigma_integrand,r,z,itype,cosm,acc,iorder))

    END FUNCTION sigma

    FUNCTION sigma_integrand(t,R,z,itype,cosm)

    !The integrand for the sigma(R) integrals
    IMPLICIT NONE
    REAL(dl) :: sigma_integrand
    REAL(dl), INTENT(IN) :: t, R, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat
    REAL(dl), PARAMETER :: alpha=alpha_sigma

    IF(t<=0. .OR. t>=1.) THEN
        !t=0 corresponds to k=infintiy when W(kR)=0.
        !t=1 corresponds to k=0. when P(k)=0.
        sigma_integrand=0.d0
    ELSE
        kR=(-1.+1./t)**alpha
        k=kR/R
        w_hat=wk_tophat(kR)
        sigma_integrand=p_lin(k,z,itype,cosm)*(w_hat**2)*alpha/(t*(1.-t))
    END IF

    END FUNCTION sigma_integrand

    FUNCTION integrate(a,b,f,y,z,itype,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL(dl) :: integrate
    REAL(dl), INTENT(IN) :: a
    REAL(dl), INTENT(IN) :: b
    REAL(dl), EXTERNAL :: f
    REAL(dl), INTENT(IN) :: y
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j
    INTEGER :: n
    REAL(dl) :: x, dx
    REAL(dl) :: f1, f2, fx
    real(dl) :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
    FUNCTION f(x, y, z, itype, cosm)
        USE precision
        IMPORT :: HM_cosmology    
        REAL(dl), INTENT(IN) :: x
        REAL(dl), INTENT(IN) :: y
        REAL(dl), INTENT(IN) :: z
        INTEGER, INTENT(IN) :: itype
        TYPE(HM_cosmology), INTENT(IN) :: cosm
    END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

        !Fix the answer to zero if the integration limits are identical
        integrate=0.d0

    ELSE

        !Reset the sum variable for the integration
        sum_2n=0.d0
        sum_n=0.d0
        sum_old=0.d0
        sum_new=0.d0

        DO j=1,jmax

            !Note, you need this to be 1+2**n for some integer n
            !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
            n=1+2**(j-1)

            !Calculate the dx interval for this value of 'n'
            dx=(b-a)/REAL(n-1)

            IF(j==1) THEN

                !The first go is just the trapezium of the end points
                f1=f(a,y,z,itype,cosm)
                f2=f(b,y,z,itype,cosm)
                sum_2n=0.5d0*(f1+f2)*dx
                sum_new=sum_2n

            ELSE

                !Loop over only new even points to add these to the integral
                DO i=2,n,2
                    x=a+(b-a)*real(i-1)/real(n-1)
                    fx=f(x,y,z,itype,cosm)
                    sum_2n=sum_2n+fx
                END DO

                !Now create the total using the old and new parts
                sum_2n=sum_n/2.d0+sum_2n*dx

                !Now calculate the new sum depending on the integration order
                IF(iorder==1) THEN
                    sum_new=sum_2n
                ELSE IF(iorder==3) THEN
                    sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
                ELSE
                    STOP 'INTEGRATE: Error, iorder specified incorrectly'
                END IF

            END IF

            IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
                !jmin avoids spurious early convergence
                integrate=REAL(sum_new)
                EXIT
            ELSE IF(j==jmax) THEN
                STOP 'INTEGRATE: Integration timed out'
            ELSE
                !Integral has not converged so store old sums and reset sum variables
                sum_old=sum_new
                sum_n=sum_2n
                sum_2n=0.d0
            END IF

        END DO

    END IF

    END FUNCTION integrate

    FUNCTION win(k,rv,c)

    !Selects the halo window function (k-space halo profile)
    REAL(dl) :: win
    REAL(dl), INTENT(IN) :: k, rv, c

    !Choose the NFW analytic form
    win=winnfw(k,rv,c)

    !Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    IF(win>1.) win=1.
    IF(win<0.) win=0.

    END FUNCTION win

    FUNCTION winnfw(k,rv,c)

    !The analytic Fourier Transform of the NFW profile
    REAL(dl) :: winnfw
    REAL(dl), INTENT(IN) :: k, rv, c
    REAL(dl) :: si1, si2, ci1, ci2, ks
    REAL(dl) :: p1, p2, p3

    !Define the scale wavenumber
    ks=k*rv/c

    !Sine and cosine integrals
    si1=Si(ks)
    si2=Si((1.+c)*ks)
    ci1=Ci(ks)
    ci2=Ci((1.+c)*ks)

    !These three parts sum to give the full W(k)
    p1=cos(ks)*(ci2-ci1)
    p2=sin(ks)*(si2-si1)
    p3=sin(ks*c)/(ks*(1.+c))

    !Create full W(k) and divide out mass factor
    winnfw=p1+p2-p3
    winnfw=winnfw/mass(c)

    END FUNCTION winnfw

    FUNCTION mass(c)

    !This calculates the (normalised) mass of a halo of concentration c
    !The 'normalised' mass is that divided by the prefactor r_s^3 4*pi rho_n
    !where rho_n is the profile normalisation [i.e, rho=rho_n/((r/r_s)*(1.+r/r_s)^2]
    REAL(dl) :: mass
    REAL(dl), INTENT(IN) :: c

    mass=log(1.+c)-c/(1.+c)

    END FUNCTION mass

    FUNCTION gnu(nu)

    !Select the mass function
    REAL(dl) :: gnu
    REAL(dl), INTENT(IN) :: nu

    !Sheth & Torman (1999)
    gnu=gst(nu)

    END FUNCTION gnu

    FUNCTION gst(nu)

    !Sheth & Tormen (1999) mass function!
    REAL(dl) :: gst
    REAL(dl), INTENT(IN) :: nu
    REAL(dl) :: p, a, bigA

    !Note I use nu=dc/sigma whereas ST (1999) use nu=(dc/sigma)^2
    !This accounts for the different pre-factor and slighly changed nu dependence
    !f(nu^2)d(nu^2)=2*nu*f(nu)dnu

    !Sheth & Tormen fitting numbers
    p=0.3
    a=0.707
    bigA=0.21616

    !Full mass function. Note this is normalised such that integral f(nu)dnu = 1
    gst=bigA*(1.+((a*nu*nu)**(-p)))*exp(-a*nu*nu/2.)

    END FUNCTION gst

    FUNCTION Hubble2(z,cosm)

    !This calculates the dimensionless squared hubble parameter at redshift z (H/H_0)^2!
    !Ignores contributions from radiation (not accurate at high z, but consistent with simulations)!
    REAL(dl) :: Hubble2
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: a

    a=1./(1.+z)

    Hubble2=cosm%om_m*(1.+z)**3+cosm%om_v*X_de(a,cosm)+(1.-cosm%om_m-cosm%om_v)*(1.+z)**2

    END FUNCTION Hubble2

    FUNCTION X_de(a,cosm)

    !The time evolution for dark energy: rho_de=rho_de,0 * X(a)
    !X(a)=1 for LCDM but changes for other models
    REAL(dl) :: X_de
    REAL(dl), INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

    END FUNCTION X_de

    FUNCTION w_de_hm(a,cosm)

    !The dark energy w(a) function
    REAL(dl) :: w_de_hm
    REAL(dl), INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    w_de_hm=cosm%w+(1.-a)*cosm%wa

    END FUNCTION w_de_hm

    FUNCTION Omega_m_hm(z,cosm)

    !This calculates omega_m variations with z!
    REAL(dl) :: Omega_m_hm
    REAL(dl), INTENT(IN) :: z
    REAL(dl) :: om_m
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    om_m=cosm%om_m
    Omega_m_hm=(om_m*(1.+z)**3.)/Hubble2(z,cosm)

    END FUNCTION Omega_m_hm

    FUNCTION Omega_cold_hm(z,cosm)

    !This calculates omega_cold variations with z!
    REAL(dl) :: Omega_cold_hm
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
  
    Omega_cold_hm=((cosm%om_c+cosm%om_b)*(1.+z)**3.)/Hubble2(z,cosm)
  
    END FUNCTION Omega_cold_hm

    FUNCTION grow(z,cosm)

    !Finds the scale-independent growth fuction at redshift z
    REAL(dl) :: grow
    REAL(dl), INTENT(IN) :: z
    REAL(dl) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(z==0.) THEN
        grow=1.
    ELSE
        a=1./(1.+z)
        grow=find(a,cosm%a_growth,cosm%growth,cosm%ng,3,3,2)
    END IF

    END FUNCTION grow

    FUNCTION sigmaV(R,z,itype,cosm)

    IMPLICIT NONE
    REAL(dl) :: sigmaV
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm   
    REAL(dl), PARAMETER :: tmin=0.
    REAL(dl), PARAMETER :: tmax=1.
    REAL(dl), PARAMETER :: acc=acc_sigmaV
    INTEGER, PARAMETER :: iorder=iorder_sigmaV

    sigmaV=sqrt(integrate(tmin,tmax,sigmaV_integrand,R,z,itype,cosm,acc,iorder)/3.)

    END FUNCTION sigmaV

    FUNCTION sigmaV_integrand(t,R,z,itype,cosm)

    !This is the integrand for the velocity dispersion integral
    IMPLICIT NONE
    REAL(dl) :: sigmaV_integrand
    REAL(dl), INTENT(IN) :: t
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat
    REAL(dl), PARAMETER :: alpha=alpha_sigmaV !Speeds up integral for large 'R'

    IF(t<=0. .OR. t>=1.) THEN
        sigmaV_integrand=0.
    ELSE
        IF(R==0.) THEN
            kR=0.
            k=(-1.+1./t)**alpha
        ELSE
            kR=(-1.+1./t)**alpha
            k=kR/R
        END IF
        w_hat=wk_tophat(kR)
        sigmaV_integrand=(p_lin(k,z,itype,cosm)/k**2)*(w_hat**2)*alpha/(t*(1.-t))
    END IF

    END FUNCTION sigmaV_integrand
  
    FUNCTION neff_integrand(t,R,z,itype,cosm)
    
    !This is the integrand for the velocity dispersion integral
    IMPLICIT NONE
    REAL(dl) :: neff_integrand
    REAL(dl), INTENT(IN) :: t
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat, w_hat_deriv
    REAL(dl), PARAMETER :: alpha=alpha_neff !Speeds up integral for large 'R'

    IF(t<=0. .OR. t>=1.) THEN
        neff_integrand=0.
    ELSE
        kR=(-1.+1./t)**alpha
        k=kR/R
        w_hat=wk_tophat(kR)
        w_hat_deriv=wk_tophat_deriv(kR)
        neff_integrand=p_lin(k,z,itype,cosm)*w_hat*w_hat_deriv*(alpha/t**2)*(-1.+1./t)**(alpha-1.)
    END IF
    
    END FUNCTION neff_integrand

    FUNCTION Si(x)

    !Calculates the 'sine integral' function Si(x)
    REAL(dl) :: Si
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: x2, y, f, g, si8
    REAL(dl), PARAMETER :: pi=3.1415926535897932384626433d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=4.) THEN

        x2=x*x

        si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

        Si=si8

    ELSE IF(ABS(x)>4.) THEN

        y=1.d0/(x*x)

        f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 +&
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10 &
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + &
            y*(1.00795182980368574617d13 + y*(4.94816688199951963482d12 +&
            y*(-4.94701168645415959931d11)))))))))))/ (x*(1. +y*(7.46437068161927678031d2 +&
            y*(1.97865247031583951450d5 +y*(2.41535670165126845144d7 + &
            y*(1.47478952192985464958d9 + y*(4.58595115847765779830d10 +&
            y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 + &
            y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))


        g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + &
            y*(3.12557570795778731d7 + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +y*(1.81004487464664575d13 +&
            y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))/&
            (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 + y*(3.26026661647090822d7 &
            + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 + y*(1.39866710696414565d12 &
            + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

        Si=pi/2.d0-f*cos(x)-g*sin(x)

    END IF

    END FUNCTION Si

    FUNCTION Ci(x)

    !Calculates the 'cosine integral' function Ci(x)
    REAL(dl) :: Ci
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: x2, y, f, g, ci8
    REAL(dl), PARAMETER :: em_const=0.577215664901532861d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=4.) THEN

        x2=x*x

        ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

        Ci=ci8

    ELSE IF(ABS(x)>4.) THEN

        y=1./(x*x)

        f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 + &
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10&
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + y*(1.00795182980368574617d13&
            + y*(4.94816688199951963482d12 +y*(-4.94701168645415959931d11)))))))))))/&
            (x*(1. +y*(7.46437068161927678031d2 +y*(1.97865247031583951450d5 +&
            y*(2.41535670165126845144d7 + y*(1.47478952192985464958d9 + &
            y*(4.58595115847765779830d10 +y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 &
            + y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))

        g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + y*(3.12557570795778731d7&
            + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +&
            y*(1.81004487464664575d13 + y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))&
            / (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 +&
            y*(3.26026661647090822d7 + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 &
            + y*(1.39866710696414565d12 + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

        Ci=f*sin(x)-g*cos(x)

    END IF

    END FUNCTION Ci

    FUNCTION find(x,xin,yin,n,iorder,ifind,imeth)

    !Given two arrays x and y this routine interpolates to find the y_i value at position x_i
    IMPLICIT NONE
    REAL(dl) :: find
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xin(n), yin(n)
    REAL(dl), ALLOCATABLE ::  xtab(:), ytab(:)
    REAL(dl) :: a, b, c, d
    REAL(dl) :: xarr(4)
    REAL(dl) :: yarr(4)
    INTEGER :: i
    INTEGER, INTENT(IN) :: iorder, ifind, imeth

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    !ifind = 1 => find x in xtab quickly assuming the table is linearly spaced
    !ifind = 2 => find x in xtab by crudely searching from x(1) to x(n)
    !ifind = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !imeth = 1 => Uses cubic polynomials for interpolation
    !imeth = 2 => Uses Lagrange polynomials for interpolation

    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
        !Reverse the arrays in this case
        CALL reverse(xtab,n)
        CALL reverse(ytab,n)
    END IF

    IF(x<xtab(1)) THEN

        !Do a linear interpolation beyond the table boundary
        IF(imeth==1) THEN
            CALL fit_line(a,b,xtab(1),ytab(1),xtab(2),ytab(2))
            find=a*x+b
        ELSE IF(imeth==2) THEN
            find=Lagrange_polynomial(x,1,xtab,ytab)
        ELSE
            STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(x>xtab(n)) THEN

        !Do a linear interpolation beyond the table boundary

        xarr(1:2) = xtab(n-1:n)
        yarr(1:2) = ytab(n-1:n)

        IF(imeth==1) THEN
            CALL fit_line(a,b,xarr(1),yarr(1),xarr(2),yarr(2))
            find=a*x+b
        ELSE IF(imeth==2) THEN
            find=Lagrange_polynomial(x,1,xarr,yarr)
        ELSE
            STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(iorder==1) THEN

        IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

        IF(x<=xtab(2)) THEN

            xarr(1:2) = xtab(1:2)
            yarr(1:2) = ytab(1:2)

        ELSE IF (x>=xtab(n-1)) THEN

            xarr(1:2) = xtab(n-1:n)
            yarr(1:2) = ytab(n-1:n)

        ELSE

            i=table_integer(x,xtab,n,ifind)

            xarr(1:2) = xtab(i:i+1)
            yarr(1:2) = ytab(i:i+1)

        END IF

        IF(imeth==1) THEN
            CALL fit_line(a,b,xarr(1),yarr(1),xarr(2),yarr(2))
            find=a*x+b
        ELSE IF(imeth==2) THEN
            find=Lagrange_polynomial(x,1,xarr,yarr)
        ELSE
            STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(iorder==2) THEN

        IF(n<3) STOP 'FIND: Not enough points in your table'

        IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

            IF(x<=xtab(2)) THEN

                xarr(1:3) = xtab(1:3)
                yarr(1:3) = ytab(1:3)

            ELSE IF (x>=xtab(n-1)) THEN

                xarr(1:3) = xtab(n-2:n)
                yarr(1:3) = ytab(n-2:n)

            END IF

            IF(imeth==1) THEN
                CALL fit_quadratic(a,b,c,xarr(1),yarr(1),xarr(2),yarr(2),xarr(3),yarr(3))
                find=a*(x**2)+b*x+c
            ELSE IF(imeth==2) THEN
                find=Lagrange_polynomial(x,2,xarr,yarr)
            ELSE
                STOP 'FIND: Error, method not specified correctly'
            END IF

        ELSE

            i=table_integer(x,xtab,n,ifind)

            xarr(1:4) = xtab(i-1:i+2)
            yarr(1:4) = ytab(i-1:i+2)

            IF(imeth==1) THEN
                !In this case take the average of two separate quadratic spline values
                CALL fit_quadratic(a,b,c,xarr(1),yarr(1),xarr(2),yarr(2),xarr(3),yarr(3))
                find=(a*x**2+b*x+c)/2.
                CALL fit_quadratic(a,b,c,xarr(2),yarr(2),xarr(3),yarr(3),xarr(4),yarr(4))
                find=find+(a*x**2+b*x+c)/2.
            ELSE IF(imeth==2) THEN
                !In this case take the average of two quadratic Lagrange polynomials
                find=(Lagrange_polynomial(x,2,xarr,yarr)+Lagrange_polynomial(x,2,xarr(2:),yarr(2:)))/2.
            ELSE
                STOP 'FIND: Error, method not specified correctly'
            END IF

        END IF

    ELSE IF(iorder==3) THEN

        IF(n<4) STOP 'FIND: Not enough points in your table'

        IF(x<=xtab(3)) THEN
            xarr(1:4) = xtab(1:4)
            yarr(1:4) = ytab(1:4)
        ELSE IF (x>=xtab(n-2)) THEN
            xarr(1:4) = xtab(n-3:n)
            yarr(1:4) = ytab(n-3:n)
        ELSE
            i=table_integer(x,xtab,n,ifind)

            xarr(1:4) = xtab(i-1:i+2)
            yarr(1:4) = ytab(i-1:i+2)
        END IF

        IF(imeth==1) THEN
            CALL fit_cubic(a,b,c,d,xarr(1),yarr(1),xarr(2),yarr(2),xarr(3),yarr(3),xarr(4),yarr(4))
            find=a*x**3+b*x**2+c*x+d
        ELSE IF(imeth==2) THEN
            find=Lagrange_polynomial(x,3,xarr,yarr)
        ELSE
            STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE

        STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

    END FUNCTION find

    FUNCTION table_integer(x,xtab,n,imeth)

    !Chooses between ways to find the integer location below some value in an array
    IMPLICIT NONE
    INTEGER :: table_integer
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER, INTENT(IN) :: imeth

    IF(imeth==1) THEN
        table_integer=linear_table_integer(x,xtab,n)
    ELSE IF(imeth==2) THEN
        table_integer=search_int(x,xtab,n)
    ELSE IF(imeth==3) THEN
        table_integer=int_split(x,xtab,n)
    ELSE
        STOP 'TABLE INTEGER: Method specified incorrectly'
    END IF

    END FUNCTION table_integer

    FUNCTION linear_table_integer(x,xtab,n)

    !Assuming the table is exactly linear this gives you the integer position
    IMPLICIT NONE
    INTEGER :: linear_table_integer
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    REAL(dl) :: x1, x2, xn
    REAL(dl), PARAMETER :: acc=1d-3 !Test for linear table

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    !n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

    END FUNCTION linear_table_integer

    FUNCTION search_int(x,xtab,n)

    !Does a stupid search through the table from beginning to end to find integer
    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER :: i

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
        IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

    END FUNCTION search_int

    FUNCTION int_split(x,xtab,n)

    !Finds the position of the value in the table by continually splitting it in half
    IMPLICIT NONE
    INTEGER :: int_split
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER :: i1, i2, imid

    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'

    i1=1
    i2=n

    DO

        imid=NINT((i1+i2)/2.)

        IF(x<xtab(imid)) THEN
            i2=imid
        ELSE
            i1=imid
        END IF

        IF(i2==i1+1) EXIT

    END DO

    int_split=i1

    END FUNCTION int_split

    SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)

    !Given xi, yi i=1,2 fits a line between these points
    IMPLICIT NONE
    REAL(dl), INTENT(OUT) :: a0, a1
    REAL(dl), INTENT(IN) :: x1, y1, x2, y2

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

    END SUBROUTINE fit_line

    SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    !Given xi, yi i=1,2,3 fits a quadratic between these points
    IMPLICIT NONE
    REAL(dl), INTENT(OUT) :: a0, a1, a2
    REAL(dl), INTENT(IN) :: x1, y1, x2, y2, x3, y3

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2)-a1*x1

    END SUBROUTINE fit_quadratic

    SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    !Given xi, yi i=1,2,3,4 fits a cubic between these points
    IMPLICIT NONE
    REAL(dl), INTENT(OUT) :: a, b, c, d
    REAL(dl), INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL(dl) :: f1, f2, f3

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2+x4*x1+x1**2)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3-b*x1**2-c*x1

    END SUBROUTINE fit_cubic

    FUNCTION Lagrange_polynomial(x,n,xv,yv)

    !Computes the result of the nth order Lagrange polynomial at point x, L(x)
    IMPLICIT NONE
    REAL(dl) :: Lagrange_polynomial
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xv(n+1), yv(n+1)
    REAL(dl) :: l(n+1)
    INTEGER :: i, j

    !Initialise variables, one for sum and one for multiplication
    Lagrange_polynomial=0.
    l=1.

    !Loops to find the polynomials, one is a sum and one is a multiple
    DO i=0,n
        DO j=0,n
            IF(i .NE. j) l(i+1)=l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
        END DO
        Lagrange_polynomial=Lagrange_polynomial+l(i+1)*yv(i+1)
    END DO

    END FUNCTION Lagrange_polynomial

    SUBROUTINE reverse(arry,n)

    !This reverses the contents of arry!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL(dl), ALLOCATABLE :: hold(:)

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
        arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

    END SUBROUTINE reverse

    SUBROUTINE fill_growtab(cosm)

    !Fills a table of values of the scale-independent growth function
    TYPE(HM_cosmology) :: cosm
    INTEGER :: i
    REAL(dl) :: a, norm
    REAL(dl), ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL(dl) :: ainit, amax, dinit, vinit
    REAL(dl) :: acc
    INTEGER, PARAMETER :: n=64

    !The calculation should start at a z when Om_m(z)=1., so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    ainit=0.001

    !Final should be a=1. unless considering models in the future
    amax=1.

    !These set the initial conditions to be the Om_m=1. growing mode
    !AM Jul 19: changed initial conditions to be appropriate for massive neutrino cosmologies
    dinit = ainit**(1.-3.*cosm%f_nu/5.)
    vinit = (1.-3.*cosm%f_nu/5.)*ainit**(-3.*cosm%f_nu/5.)

    !Overall accuracy for the ODE solver
    acc=1d-3

    IF(HM_verbose) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,0.d0,ainit,amax,dinit,vinit,acc,3,cosm)
    IF(HM_verbose) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    norm=find(1.d0,a_tab,d_tab,SIZE(a_tab),3,3,2)
    IF(HM_verbose) WRITE(*,*) 'GROWTH: Unnormalised g(a=1):', norm
    d_tab=d_tab/norm

    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_growth)) DEALLOCATE(cosm%a_growth)
    IF(ALLOCATED(cosm%growth)) DEALLOCATE(cosm%growth)

    cosm%ng=n
    ALLOCATE(cosm%a_growth(n),cosm%growth(n))
    DO i=1,n
        a=ainit+(amax-ainit)*float(i-1)/float(n-1)
        cosm%a_growth(i)=a
        cosm%growth(i)=find(a,a_tab,d_tab,SIZE(a_tab),3,3,2)
    END DO

    IF(HM_verbose) WRITE(*,*) 'GROWTH: Done'
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE fill_growtab

    SUBROUTINE ode_growth(x,v,t,kk,ti,tf,xi,vi,acc,imeth,cosm)

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    IMPLICIT NONE
    REAL(dl) :: xi, ti, tf, dt, acc, vi, x4, v4, t4, kk
    REAL(dl) :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL(dl), ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL(dl), ALLOCATABLE :: x(:), v(:), t(:)
    INTEGER :: i, j, k, n, np, ifail, kn, imeth
    TYPE(HM_cosmology) :: cosm
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)

    DO j=1,jmax

        !Set the number of points for the forward integration
        n=ninit*(2**(j-1))
        n=n+1

        !Allocate arrays
        ALLOCATE(x8(n),t8(n),v8(n))

        !Set the arrays to initialy be zeroes (is this neceseary?)
        x8=0.d0
        t8=0.d0
        v8=0.d0

        !Set the intial conditions at the intial time
        x8(1)=xi
        v8(1)=vi

        !Fill up a table for the time values
        CALL fill_table(ti,tf,t8,n)

        !Set the time interval
        dt=(tf-ti)/float(n-1)

        !Intially fix this to zero. It will change to 1 if method is a 'failure'
        ifail=0

        DO i=1,n-1

            x4=real(x8(i))
            v4=real(v8(i))
            t4=real(t8(i))

            IF(imeth==1) THEN

                !Crude method
                kx1=dt*fd(x4,v4,kk,t4,cosm)
                kv1=dt*fv(x4,v4,kk,t4,cosm)

                x8(i+1)=x8(i)+kx1
                v8(i+1)=v8(i)+kv1

            ELSE IF(imeth==2) THEN

                !Mid-point method
                !2017/06/18 - There was a bug in this part before. Luckily it was not used. Thanks Dipak Munshi.
                kx1=dt*fd(x4,v4,kk,t4,cosm)
                kv1=dt*fv(x4,v4,kk,t4,cosm)
                kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
                kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)

                x8(i+1)=x8(i)+kx2
                v8(i+1)=v8(i)+kv2

            ELSE IF(imeth==3) THEN

                !4th order Runge-Kutta method (fast!)
                kx1=dt*fd(x4,v4,kk,t4,cosm)
                kv1=dt*fv(x4,v4,kk,t4,cosm)
                kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
                kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
                kx3=dt*fd(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
                kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
                kx4=dt*fd(x4+kx3,v4+kv3,kk,t4+dt,cosm)
                kv4=dt*fv(x4+kx3,v4+kv3,kk,t4+dt,cosm)

                x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
                v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.

            END IF

        END DO

        IF(j==1) ifail=1

        IF(j .NE. 1) THEN

            np=1+(n-1)/2

            DO k=1,1+(n-1)/2

                kn=2*k-1

                IF(ifail==0) THEN

                    IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                    IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                    IF(ifail==1) THEN
                        DEALLOCATE(xh,th,vh)
                        EXIT
                    END IF

                END IF
            END DO

        END IF

        IF(ifail==0) THEN
            ALLOCATE(x(n),t(n),v(n))
            x=real(x8)
            v=real(v8)
            t=real(t8)
            EXIT
        END IF

        ALLOCATE(xh(n),th(n),vh(n))
        xh=x8
        vh=v8
        th=t8
        DEALLOCATE(x8,t8,v8)

    END DO

    END SUBROUTINE ode_growth

    FUNCTION fv(d,v,k,a,cosm)

    !v'=f(v) in ODE solver
    REAL(dl) :: fv
    REAL(dl), INTENT(IN) :: d, v, k, a
    REAL(dl) :: f1, f2, z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    z=-1.+(1./a)

    !AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    f1=3.*Omega_cold_hm(z,cosm)*d/(2.*(a**2.))
    f2=(2.+AH(z,cosm)/Hubble2(z,cosm))*(v/a)

    fv=f1-f2

    END FUNCTION fv

    FUNCTION fd(d,v,k,a,cosm)
    !d'=f(d) in ODE solver
    REAL(dl) :: fd
    REAL(dl), INTENT(IN) :: d, v, k, a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    fd=v

    END FUNCTION fd

    FUNCTION AH(z,cosm)

    !The Hubble acceleration function \ddot{a}/a
    REAL(dl) :: AH
    REAL(dl), INTENT(IN) :: z
    REAL(dl) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    a=1./(1.+z)
    AH=cosm%om_m*(a**(-3.))+cosm%om_v*(1.+3.*w_de_hm(a,cosm))*x_de(a,cosm)
    AH=-AH/2.

    END FUNCTION AH

    !!AM End HMcode

    subroutine PKequal(State,redshift,w_lam,wa_ppf,w_hf,wa_hf)
    !used by halofit_casarini: arXiv:0810.0190, arXiv:1601.07230
    implicit none
    Type(CAMBdata) :: State
    real(dl) :: redshift,w_lam,wa_ppf,w_hf,wa_hf
    real(dl) :: z_star,tau_star,dlsb,dlsb_eq,w_true,wa_true,error

    z_star=State%ThermoDerivedParams( derived_zstar )
    tau_star=State%TimeOfz(z_star)
    dlsb=State%TimeOfz(redshift)-tau_star
    w_true=w_lam
    wa_true=wa_ppf
    wa_ppf=0._dl
    do
        z_star=State%ThermoDerivedParams( derived_zstar )
        tau_star=State%TimeOfz(State%ThermoData%z_star)
        dlsb_eq=State%TimeOfz(redshift)-tau_star
        error=1.d0-dlsb_eq/dlsb
        if (abs(error).le.1e-7) exit
        w_lam=w_lam*(1+error)**10.d0
    enddo
    w_hf=w_lam
    wa_hf=0._dl
    w_lam=w_true
    wa_ppf=wa_true
    write(*,*)'at z = ',real(redshift),' equivalent w_const =', real(w_hf)

    end subroutine PKequal


    end module NonLinear
