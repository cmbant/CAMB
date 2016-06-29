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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module NonLinear
    use ModelParams
    use DarkEnergyInterface
    use transfer
    implicit none
    private

    real, parameter :: Min_kh_nonlinear = 0.005
    real(dl):: om_m,om_v,fnu,omm0, acur

    integer, parameter :: halofit_original = 1, halofit_bird=2, halofit_peacock=3, halofit_takahashi=4
    integer, parameter :: halofit_mead=5, halofit_halomodel=6
    integer, parameter :: halofit_default = halofit_takahashi
    integer :: halofit_version = halofit_default
    public Min_kh_nonlinear, NonLinear_GetNonLinRatios, NonLinear_ReadParams
    public halofit_version, halofit_default, halofit_original, halofit_bird, halofit_peacock, halofit_takahashi
    public halofit_mead, halofit_halomodel

    !!AM - Added these types for HMcode
    INTEGER :: imead, ihm !!AM - added these for HMcode, need to be visible to all subroutines and functions

    TYPE HM_cosmology
        !Contains only things that do not need to be recalculated with each new z
        REAL :: om_m, om_v, w, wa, f_nu, ns, h, Tcmb, Nnu
        REAL, ALLOCATABLE :: r_sigma(:), sigma(:)
        REAL, ALLOCATABLE :: growth(:), a_growth(:)
        REAL, ALLOCATABLE :: k_plin(:), plin(:), plinc(:)
    END TYPE HM_cosmology

    TYPE HM_tables
        !Stuff that needs to be recalculated for each new z
        REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
        REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
        INTEGER :: n
    END TYPE HM_tables
    !!AM - End of my additions

    contains

    subroutine NonLinear_ReadParams(Ini)
    use IniObjects
    Type(TIniFile) :: Ini

    halofit_version = Ini%Read_Int('halofit_version', halofit_default)

    end subroutine NonLinear_ReadParams

    subroutine NonLinear_GetNonLinRatios(CAMB_Pk)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    !This implementation uses Halofit
    type(MatterPowerData) :: CAMB_Pk
    integer itf
    real(dl) a,plin,pq,ph,pnl,rk
    real(dl) sig,rknl,rneff,rncur,d1,d2
    real(dl) diff,xlogr1,xlogr2,rmid
    integer i

    IF(halofit_version==halofit_mead .OR. halofit_version==halofit_halomodel) THEN

        !AM - Call HMcode here
        CALL HMcode(CAMB_Pk)

    ELSE

        !!BR09 putting neutrinos into the matter as well, not sure if this is correct, but at least one will get a consisent omk.
        omm0 = CP%omegac+CP%omegab+CP%omegan
        fnu = CP%omegan/omm0

        CAMB_Pk%nonlin_ratio = 1

        do itf = 1, CAMB_Pk%num_z


            ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and
            ! curvature (rncur) of the power spectrum at the desired redshift, using method
            ! described in Smith et al (2002).
            a = 1/real(1+CAMB_Pk%Redshifts(itf),dl)
            om_m = omega_m(a, omm0, CP%omegav, CP%DarkEnergy%w_lam, CP%DarkEnergy%wa_ppf)
            om_v = omega_v(a, omm0, CP%omegav, CP%DarkEnergy%w_lam, CP%DarkEnergy%wa_ppf)
            acur = a
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
                else if (xlogr2>3.4999) then
                    ! Totally crazy non-linear
                    global_error_flag=349
                    write(*,*) 'Error in halofit'
                    goto 101
                end if
            end do

            ! now calculate power spectra for a logarithmic range of wavenumbers (rk)

            do i=1, CAMB_PK%num_k
                rk = exp(CAMB_Pk%log_kh(i))

                if (rk > Min_kh_nonlinear) then

                    ! linear power spectrum !! Remeber => plin = k^3 * P(k) * constant
                    ! constant = 4*pi*V/(2*pi)^3

                    plin= MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3/(2*const_pi**2))

                    ! calculate nonlinear power according to halofit: pnl = pq + ph,
                    ! where pq represents the quasi-linear (halo-halo) power and
                    ! where ph is represents the self-correlation halo term.

                    call halofit(rk,rneff,rncur,rknl,plin,pnl,pq,ph)   ! halo fitting formula
                    CAMB_Pk%nonlin_ratio(i,itf) = sqrt(pnl/plin)

                end if

            enddo

101         continue
        end do

    END IF

    end subroutine NonLinear_GetNonLinRatios

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subroutine halofit(rk,rn,rncur,rknl,plin,pnl,pq,ph)

    real(dl) gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
    real(dl) rk,rn,plin,pnl,pq,ph,plinaa
    real(dl) rknl,y,rncur
    real(dl) f1a,f2a,f3a,f1b,f2b,f3b,frac
    real(dl) extragam, peacock_fudge

    if (halofit_version ==halofit_original .or. halofit_version ==halofit_bird &
        .or. halofit_version == halofit_peacock) then
    ! halo model nonlinear fitting formula as described in
    ! Appendix C of Smith et al. (2002)
    !SPB11: Standard halofit underestimates the power on the smallest scales by a
    !factor of two. Add an extra correction from the simulations in Bird, Viel,
    !Haehnelt 2011 which partially accounts for this.
    if (halofit_version ==halofit_bird) then
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
    beta=0.8291+0.9854*rn+0.3400*rn**2+fnu*(-6.4868+1.4373*rn**2)
    elseif (halofit_version == halofit_takahashi) then
        !RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
        !than latest N-body simulations at small scales.
        !Update the following fitting parameters of gam,a,b,c,xmu,xnu,
        !alpha & beta from the simulations in Takahashi+ 2012.
        !The improved halofit accurately provide the power spectra for WMAP
        !cosmological models with constant w.
        gam=0.1971-0.0843*rn+0.8460*rncur
        a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+ &
            0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*om_v* &
            (1.+CP%DarkEnergy%w_lam+CP%DarkEnergy%wa_ppf*(1-acur))
        a=10**a
        b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+ &
            0.2279*om_v*(1.+CP%DarkEnergy%w_lam+CP%DarkEnergy%wa_ppf*(1-acur)))
        c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
        xmu=0.
        xnu=10**(5.2105+3.6902*rn)
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ &
            0.3980*rn**4-0.1682*rncur + fnu*(1.081 + 0.395*rn**2)
    else
        call MpiStop('Unknown halofit_version')
    end if

    if(abs(1-om_m).gt.0.01) then ! omega evolution
        f1a=om_m**(-0.0732)
        f2a=om_m**(-0.1423)
        f3a=om_m**(0.0725)
        f1b=om_m**(-0.0307)
        f2b=om_m**(-0.0585)
        f3b=om_m**(0.0743)
        frac=om_v/(1.-om_m)
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
    ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fnu*0.977)
    plinaa=plin*(1+fnu*47.48*rk**2/(1+1.5*rk**2))
    pq=plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)

    pnl=pq+ph

    if (halofit_version == halofit_peacock) then
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
    SUBROUTINE HMcode(CAMB_Pk)
    !!AM - A CAMB derived type that I need
    TYPE(MatterPowerData) :: CAMB_Pk
    REAL :: z, k
    REAL :: p1h, p2h, pfull, plin
    INTEGER :: i, j, nk, nz
    REAL, PARAMETER :: pi=3.141592654
    TYPE(HM_cosmology) :: cosi
    TYPE(HM_tables) :: lut

    !HMcode developed by Alexander Mead (alexander.j.mead@googlemail.com)
    !Please contact me if you have any questions whatsoever
    !If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
    !If you use the extensions (w(a) and massive neutrinos) then please cite: http://arxiv.org/abs/1602.02154
    !Also consider citing the source code at ASCL: http://ascl.net/1508.001

    !Use imead to switch between the standard and accurate halo-model calcuation
    !0 - Standard (this is just a vanilla halo model calculation with no accuracy tweaks)
    !1 - Accurate from Mead et al. (2015; arXiv 1505.07833)
    IF(halofit_version==halofit_halomodel) imead=0
    IF(halofit_version==halofit_mead) imead=1

    !Use ihm to switch between verbose (diagnostic) and non-verbose mode
    !0 - Non-verbose
    !1 - Verbose
    ihm=0

    IF(ihm==1) WRITE(*,*)
    IF(ihm==1) WRITE(*,*) 'Running HMcode'
    IF(ihm==1) WRITE(*,*)

    !!AM - Translate from CAMB variables to my variables
    nz=CAMB_PK%num_z
    nk=CAMB_PK%num_k

    !!AM - Assign cosmological parameters for the halo model calculation
    CALL assign_HM_cosmology(cosi)

    !Fill growth function table (only needs to be done once)
    CALL fill_growtab(cosi)

    !Loop over redshifts
    DO j=1,nz

        !Initialise the specific HM_cosmology (fill sigma(R) and P_lin HM_tables)
        !Currently this needs to be done at each z (mainly because of scale-dependent growth with neutrinos)
        !For non neutrino models this could only be done once, which would speed things up a bit
        CALL initialise_HM_cosmology(j,cosi,CAMB_PK)

        !Sets the current redshift from the table
        z=CAMB_Pk%Redshifts(j)

        !Initiliasation for the halomodel calcualtion (needs to be done for each z)
        CALL halomod_init(z,lut,cosi)

        !Loop over k values and calculate P(k)
        !$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin, pfull,p1h,p2h)
        DO i=1,nk
            k=exp(CAMB_Pk%log_kh(i))
            plin=p_lin(k,z,0,cosi)
            CALL halomod(k,z,p1h,p2h,pfull,plin,lut,cosi)
            CAMB_Pk%nonlin_ratio(i,j)=sqrt(pfull/plin)
        END DO
        !$OMP END PARALLEL DO
    END DO

    END SUBROUTINE HMcode

    FUNCTION Delta_v(z,lut,cosm)

    !Function for the virialised overdensity
    REAL :: Delta_v
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(imead==0) THEN
        !Value that is normally used in halo model predictions
        Delta_v=200.
    ELSE IF(imead==1) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        Delta_v=418.*(omega_m_hm(z,cosm)**(-0.352))
        !Mead et al. (2016; arXiv 1602.02154) neutrino addition
        Delta_v=Delta_v*(1.+0.916*cosm%f_nu)
    END IF

    END FUNCTION Delta_v

    FUNCTION delta_c(z,lut,cosm)

    !Function for the linear collapse density
    REAL :: delta_c
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(imead==0) THEN
        delta_c=1.686
    ELSE IF(imead==1) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        delta_c=1.59+0.0314*log(lut%sig8z)
        !Mead et al. (2016; arXiv 1602.02154) neutrino addition
        delta_c=delta_c*(1.+0.262*cosm%f_nu)
    END IF

    !Nakamura & Suto (1997) fitting formula for LCDM models
    delta_c=delta_c*(1.+0.0123*log10(omega_m_hm(z,cosm)))

    END FUNCTION delta_c

    FUNCTION eta(z,lut,cosm)

    !Function eta that puffs halo profiles
    REAL :: eta
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(imead==0) THEN
        eta=0.
    ELSE IF(imead==1) THEN
        !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
        eta=0.603-0.3*lut%sig8z
    END IF

    END FUNCTION eta

    FUNCTION kstar(z,lut,cosm)

    !Function k* that cuts off the 1-halo term at large scales
    REAL :: kstar
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(imead==0) THEN
        !Set to zero for the standard Poisson one-halo term
        kstar=0.
    ELSE IF(imead==1) THEN
        !One-halo cut-off wavenumber
        !Mead et al. (2015; arXiv 1505.07833) value
        kstar=0.584*(lut%sigv)**(-1.)
    END IF

    END FUNCTION kstar

    FUNCTION As(z,lut,cosm)

    !Halo concentration pre-factor from Bullock et al. (2001) relation
    REAL :: As
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(imead==0) THEN
        !Set to 4 for the standard Bullock value
        As=4.
    ELSE IF(imead==1) THEN
        !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
        As=3.13
    END IF

    END FUNCTION As

    FUNCTION fdamp(z,lut,cosm)

    !Linear power damping function from Mead et al. (2015; arXiv 1505.07833)
    REAL ::fdamp
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    !Linear theory damping factor
    IF(imead==0) THEN
        !Set to 0 for the standard linear theory two halo term
        fdamp=0.
    ELSE IF(imead==1) THEN
        !Mead et al. (2016; arXiv 1602.02154) value
        fdamp=0.0095*lut%sigv100**1.37
    END IF

    !Catches extreme values of fdamp
    IF(fdamp<1.e-3) fdamp=1.e-3
    IF(fdamp>0.99)  fdamp=0.99

    END FUNCTION fdamp

    FUNCTION alpha(z,lut,cosm)

    !Two- to one-halo transition smoothing from Mead et al. (2015; arXiv 1505.07833)
    REAL :: alpha
    REAL, INTENT(IN) :: z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(imead==0) THEN
        !Set to 1 for the standard halomodel sum of one- and two-halo terms
        alpha=1.
    ELSE IF(imead==1) THEN
        !This uses the top-hat defined neff (HALOFIT uses Gaussian filtered fields instead)
        !Mead et al. (2016; arXiv 1602.02154) value
        alpha=3.24*1.85**lut%neff
    END IF

    !Catches values of alpha that are crazy
    IF(alpha>2.)  alpha=2.
    IF(alpha<0.5) alpha=0.5

    END FUNCTION alpha

    FUNCTION r_nl(lut)

    !Calculates R_nl, defined by nu(R_nl)=1., nu=dc/sigma(R)
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL :: r_nl

    IF(lut%nu(1)>1.) THEN
        !This catches some very strange values that appear for odd cosmological models
        r_nl=lut%rr(1)
    ELSE
        r_nl=exp(find(log(1.),log(lut%nu),log(lut%rr),3,3))
    END IF

    END FUNCTION r_nl

    SUBROUTINE halomod(k,z,p1h,p2h,pfull,plin,lut,cosm)

    !Calcuates 1-halo and 2-halo terms and combines them to form the full halomodel power
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    REAL :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN
        p1h=0.
        p2h=0.
    ELSE
        p1h=p_1h(k,z,lut,cosm)
        p2h=p_2h(k,z,plin,lut,cosm)
    END IF

    a=alpha(z,lut,cosm)
    pfull=(p2h**a+p1h**a)**(1./a)

    END SUBROUTINE halomod

    SUBROUTINE fill_table(min,max,arr,n,ilog)

    !Fills array 'arr' in equally spaced intervals
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL :: a, b
    REAL, ALLOCATABLE, INTENT(OUT) :: arr(:)
    INTEGER, INTENT(IN) :: ilog, n

    !ilog=0 does linear spacing
    !ilog=1 does log spacing

    !Allocate arrays
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))

    !This is probably unnecessary
    arr=0.

    !Decide on linear or log spacing
    IF(ilog==0) THEN
        a=min
        b=max
    ELSE IF(ilog==1) THEN
        a=log(min)
        b=log(max)
    END IF

    !Fill the array
    IF(n==1) THEN
        !This should not be necessary
        arr(1)=a
    ELSE IF(n>1) THEN
        DO i=1,n
            arr(i)=a+(b-a)*float(i-1)/float(n-1)
        END DO
    END IF

    !If your are filling in log
    IF(ilog==1) arr=exp(arr)

    END SUBROUTINE fill_table

    SUBROUTINE fill_plintab(iz,cosm,CAMB_PK)

    !Fills internal HMcode HM_tables for the linear power spectrum at z=0
    TYPE(MatterPowerData), INTENT(IN) :: CAMB_PK
    INTEGER, INTENT(IN) :: iz
    TYPE(HM_cosmology) :: cosm
    INTEGER :: i, nk
    INTEGER :: imeth
    REAL :: z, g, k, kmin, kmax
    REAL, PARAMETER :: pi=3.141592654

    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: Filling linear power HM_tables'

    !Fill arrays
    IF(ALLOCATED(cosm%k_plin)) DEALLOCATE(cosm%k_plin)
    IF(ALLOCATED(cosm%plin))   DEALLOCATE(cosm%plin)
    IF(ALLOCATED(cosm%plinc))  DEALLOCATE(cosm%plinc)

    imeth=2
    IF(imeth==1) THEN

        !Fill k-table with the same k points as in the CAMB calculation
        !If a user has specified lots of points this could make the halo-model
        !calculation chug
        nk=CAMB_PK%num_k
        ALLOCATE(cosm%k_plin(nk))
        DO i=1,nk
            cosm%k_plin(i)=exp(CAMB_Pk%log_kh(i))
        END DO

    ELSE IF(imeth==2) THEN

        !Fill a k-table with an equal-log-spaced k range
        !Note that the minimum should be such that the spectrum is accurately a power-law below this wavenumber
        kmin=1e-3
        kmax=1e2
        nk=128
        CALL fill_table(kmin,kmax,cosm%k_plin,nk,1)

    END IF

    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: k_min:', cosm%k_plin(1)
    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: k_max:', cosm%k_plin(nk)
    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: nk:', nk

    ALLOCATE(cosm%plin(nk),cosm%plinc(nk))

    !Fill power table
    DO i=1,nk
        !Take the power from the current redshift choice
        cosm%plin(i)=MatterPowerData_k(CAMB_PK,DBLE(cosm%k_plin(i)),iz)*(cosm%k_plin(i)**3/(2.*pi**2))
        cosm%plinc(i)=cosm%plin(i)*(Tcb_Tcbnu_ratio(cosm%k_plin(i),z,cosm))**2.
    END DO

    !Find the redshift
    z=CAMB_Pk%Redshifts(iz)
    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: z of input:', z

    !Calculate the growth factor at the redshift of interest
    g=grow(z,cosm)

    !Grow the power to z=0
    cosm%plin=cosm%plin/(g**2.)
    cosm%plinc=cosm%plinc/(g**2.)

    !Check sigma_8 value
    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: Sigma_8:', sigma(8.,0.,0,cosm)
    IF(ihm==1) WRITE(*,*) 'LINEAR POWER: Done'
    IF(ihm==1) WRITE(*,*)

    END SUBROUTINE fill_plintab

    FUNCTION Tcb_Tcbnu_ratio(k,z,cosm)

    !Calculates the ratio of T(k) for cold vs. all matter
    !Uses approximations in Eisenstein & Hu (1999; arXiv 9710252)
    !Note that this assumes that there are exactly 3 species of neutrinos with
    !Nnu<=3 of these being massive, and with the mass split evenly between the number of massive species.

    REAL :: Tcb_Tcbnu_ratio
    REAL, INTENT(IN) :: k, z
    REAL :: D, Dcb, Dcbnu, pcb, zeq, q, yfs
    REAL :: BigT
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
        !D=(1.+zeq)*grow(z,cosm)
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

    SUBROUTINE assign_HM_cosmology(cosm)

    !Assigns the internal HMcode cosmological parameters
    TYPE(HM_cosmology) :: cosm

    !Converts CAMB parameters to Meadfit parameters
    cosm%om_m=CP%omegac+CP%omegab+CP%omegan
    cosm%om_v=CP%omegav
    cosm%w=CP%DarkEnergy%w_lam
    cosm%wa=CP%DarkEnergy%wa_ppf
    cosm%f_nu=CP%omegan/cosm%om_m
    cosm%h=CP%H0/100.
    cosm%Tcmb=CP%tcmb
    cosm%Nnu=CP%Num_Nu_massive

    !n_s is read in here. The non-linear CAMB module does not work if there is more than
    !one value in this array, so explicitly setting '1' here is fine.
    cosm%ns=CP%InitPower%an(1)

    !Write out cosmological parameters if necessary
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: Om_m:', cosm%om_m
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: Om_v:', cosm%om_v
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: w_0:', cosm%w
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: w_a:', cosm%wa
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: f_nu:', cosm%f_nu
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: n_s:', cosm%ns
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: h:', cosm%h
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: T_cmb:', cosm%Tcmb
    IF(ihm==1) WRITE(*,*) 'HM_cosmology: N_nu (massive):', cosm%Nnu
    IF(ihm==1) WRITE(*,*)

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

    SUBROUTINE halomod_init(z,lut,cosm)

    !Halo-model initialisation routine
    !Computes look-up HM_tables necessary for the halo model calculations
    REAL, INTENT(IN) :: z
    INTEGER :: i, n
    REAL :: Dv, dc, f, m, mmin, mmax, nu, r, sig
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut

    IF(ihm==1) WRITE(*,*) 'HALOMOD: Filling look-up HM_tables'
    IF(ihm==1) WRITE(*,*) 'HALOMOD: HM_tables being filled at redshift:', z

    !Find value of sigma_v, sig8, etc.
    lut%sigv=sqrt(dispint(z,cosm))
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
    lut%sigv100=sigma_v(100.,z,cosm)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
    lut%sig8z=sigma(8.,z,0,cosm)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sig8(z):', lut%sig8z

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    !Number of entries in look-up HM_tables. Could be played with to improve speed or accuracy
    n=256

    lut%n=n
    CALL allocate_lut(lut)

    !Sets the mass range for halo model calculation
    !Default is 1e0 to 1e18, I cannot believe this would ever be insufficient
    mmin=1.e0
    mmax=1.e18

    IF(ihm==1) WRITE(*,*) 'HALOMOD: M_min:', mmin
    IF(ihm==1) WRITE(*,*) 'HALOMOD: M_max:', mmax

    dc=delta_c(z,lut,cosm)

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

    IF(ihm==1) WRITE(*,*) 'HALOMOD: m, r, nu, sig HM_tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
        lut%sigf(i)=sigmac(lut%rr(i)*f,z,cosm)
    END DO
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigf HM_tables filled'

    !Fill virial radius table using real radius table
    Dv=Delta_v(z,lut,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(ihm==1) WRITE(*,*) 'HALOMOD: rv HM_tables filled'
    IF(ihm==1) WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: M min [Msun/h]:', lut%m(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: M max [Msun/h]:', lut%m(lut%n)

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1./lut%rnl

    IF(ihm==1) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    IF(ihm==1) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    !Calcuate the effective spectral index at the collapse scale
    lut%neff=neff(lut,cosm)

    IF(ihm==1) WRITE(*,*) 'HALOMOD: n_eff:', lut%neff

    !Get the concentration for all the haloes
    CALL conc_bull(z,lut,cosm)

    IF(ihm==1) WRITE(*,*) 'HALOMOD: c HM_tables filled'
    IF(ihm==1) WRITE(*,*) 'HALOMOD: c min [Msun/h]:', lut%c(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: c max [Msun/h]:', lut%c(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: Done'
    IF(ihm==1) WRITE(*,*)
    IF(ihm==1) CALL write_parameters(z,lut,cosm)

    !Switch off verbose mode if doing multiple z
    ihm=0

    END SUBROUTINE halomod_init

    SUBROUTINE write_parameters(z,lut,cosm)

    !This subroutine writes out the halomodel parameters at the current redshift
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(ihm==1) WRITE(*,*) 'Parameters at your redshift'
    IF(ihm==1) WRITE(*,*) '==========================='
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'z:', z
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'Dv:', Delta_v(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'dc:', delta_c(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'eta:', eta(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'k*:', kstar(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'A:', As(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'fdamp:', fdamp(z,lut,cosm)
    IF(ihm==1) WRITE(*,fmt='(A10,F10.5)') 'alpha:', alpha(z,lut,cosm)
    IF(ihm==1) WRITE(*,*)

    END SUBROUTINE write_parameters

    PURE FUNCTION radius_m(m,cosm)

    !Calculates the co-moving radius that encloses a mass 'm' in the homogeneous Universe
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    radius_m=(3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

    END FUNCTION radius_m

    FUNCTION neff(lut,cosm)

    !Finds the effective spectral index at the collapse scale r_nl
    !Where nu(r_nl)=1.
    REAL :: neff
    REAL :: ns
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    !Numerical differentiation to find effective index at collapse
    neff=-3.-derivative_table(log(lut%rnl),log(lut%rr),log(lut%sig**2.),3,3)

    !For some bizarre cosmological models r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    ns=cosm%ns
    IF(neff<ns-4.) neff=ns-4.
    IF(neff>ns)    neff=ns

    END FUNCTION neff

    SUBROUTINE conc_bull(z,lut,cosm)

    !Calculates the Bullock et al. (2001) concentration-mass relation
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology) :: cosm, cos_lcdm
    TYPE(HM_tables) :: lut
    REAL :: A, zf, ainf, zinf, g_lcdm, g_wcdm, pow
    INTEGER :: i

    !Amplitude of relation (4. in Bullock et al. 2001)
    A=As(z,lut,cosm)

    !Fill the collapse time look-up table
    CALL zcoll_bull(z,cosm,lut)

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
    cos_lcdm=cosm
    DEALLOCATE(cos_lcdm%growth)
    DEALLOCATE(cos_lcdm%a_growth)
    cos_lcdm%w=-1.
    cos_lcdm%wa=0.

    ainf=1./(1.+zinf)

    !Needs to use grow_int explicitly here for LCDM model to avoid growth HM_tables
    g_lcdm=grow_int(ainf,0.001,cos_lcdm)

    pow=1.5
    !This is the Dolag et al. (2004) correction with a 1.5 power rather than 1
    lut%c=lut%c*((g_wcdm/g_lcdm)**pow)

    END SUBROUTINE conc_bull

    FUNCTION grow_int(b,acc,cosm)

    !Integrates the (very good) approximation for the growth rate (f=\Om^\gamma)
    !to get the growth function. This only works for wCDM modesls with constant w
    !but is only used for these models.
    INTEGER :: i, j, jmax
    REAL :: grow_int, a, b, acc, dx
    INTEGER :: nint
    REAL :: x, fac, func, gam
    REAL(dl) :: sum1, sum2
    TYPE(HM_cosmology) :: cosm

    sum1=0.d0
    sum2=0.d0

    jmax=20

    a=1.

    IF(a==b) THEN

        grow_int=1.

    ELSE

        DO j=1,jmax

            nint=10.*(2.**j)

            DO i=1,nint

                x=a+(b-a)*((float(i)-1)/(float(nint)-1))

                IF(i==1 .OR. i==nint) THEN
                    !multiple of 1 for beginning and end and multiple of 2 for middle points!
                    fac=1.
                ELSE
                    fac=2.
                END IF

                !Growth rate is Omega_m(z)^0.55 to a very good approximation
                !Some small deivations for wCDM are included
                !Not calibrated to w(a) models!
                IF(cosm%w<-1.) THEN
                    gam=0.55+0.02*(1.+cosm%w)
                ELSE IF(cosm%w>-1) THEN
                    gam=0.55+0.05*(1.+cosm%w)
                ELSE
                    gam=0.55
                END IF

                func=(omega_m_hm(-1.+1./x,cosm)**gam)/x

                sum2=sum2+fac*func

            END DO

            dx=((b-a)/(float(nint)-1.))
            sum2=sum2*dx/2.

            IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
                grow_int=exp(sum2)
                EXIT
            ELSE
                sum1=sum2
                sum2=0.d0
            END IF

        END DO

    END IF

    END FUNCTION grow_int

    SUBROUTINE zcoll_bull(z,cosm,lut)

    !Calcuates the halo collapse redshift according to the Bullock (2001) prescription
    REAL, INTENT(IN) :: z
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut
    REAL :: dc
    REAL :: af, zf, RHS, a, growz
    INTEGER :: i, j

    !This fills up the halo formation redshift table as per Bullock relations

    !Needs to interpolate g(z) which should be pretty linear for a<0.05
    !in 'g(a) vs a' space for all standard cosmologies

    dc=delta_c(z,lut,cosm)

    !Find the growth function at the current redshift
    a=1./(1.+z)
    growz=find(a,cosm%a_growth,cosm%growth,3,3)

    !Do numerical inversion
    DO i=1,lut%n

        RHS=dc*grow(z,cosm)/lut%sigf(i)

        IF(RHS>growz) THEN
            !This is the case of 'halo forms in the future'
            !in this case set formation redshift to current redshift
            zf=z
        ELSE
            af=find(RHS,cosm%growth,cosm%a_growth,3,3)
            zf=-1.+1./af
        END IF

        lut%zc(i)=zf

    END DO

    END SUBROUTINE zcoll_bull

    FUNCTION mass_r(r,cosm)

    !Calcuates the average mass enclosed at co-moving radius r
    REAL :: mass_r, r
    TYPE(HM_cosmology) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    !Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*cosmic_density(cosm)*(r**3.)

    END FUNCTION mass_r

    PURE FUNCTION cosmic_density(cosm)

    !The z=0 cosmological matter density
    REAL :: cosmic_density
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !In M_sun per Mpc^3 with h factors included. The constant does this.
    cosmic_density=(2.775e11)*cosm%om_m

    END FUNCTION cosmic_density

    FUNCTION find_pk(k,itype,cosm)

    !Look-up and interpolation for P(k,z=0)
    REAL :: find_pk
    REAL :: kmax, ns
    REAL, INTENT(IN) :: k
    INTEGER, INTENT(IN) :: itype
    INTEGER :: n
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Set number of k points as well as min and max k values
    !Note that the min k value should be set to the same as the CAMB min k value
    n=SIZE(cosm%k_plin)
    kmax=cosm%k_plin(n)

    !Spectral index used in the high-k extrapolation
    ns=cosm%ns

    IF(k>kmax) THEN
        !Do some interpolation here based on knowledge of things at high k
        IF(itype==0) THEN
            find_pk=cosm%plin(n)*((log(k)/log(kmax))**2.)*((k/kmax)**(ns-1.))
        ELSE IF(itype==1) THEN
            find_pk=cosm%plinc(n)*((log(k)/log(kmax))**2.)*((k/kmax)**(ns-1.))
        END IF
    ELSE
        !Otherwise use the standard find algorithm
        IF(itype==0) THEN
            find_pk=exp(find(log(k),log(cosm%k_plin),log(cosm%plin),3,3))
        ELSE IF(itype==1) THEN
            find_pk=exp(find(log(k),log(cosm%k_plin),log(cosm%plinc),3,3))
        END IF
    END IF

    !Old method, works fine for m_nu<0.5 eV
    !IF(itype==1) find_pk=find_pk/(1.-cosm%f_nu)**2.

    END FUNCTION find_pk

    FUNCTION p_lin(k,z,itype,cosm)

    !Looks up the value for the linear power spectrum
    REAL :: p_lin
    REAL, INTENT (IN) :: k, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised to z=0

    p_lin=(grow(z,cosm)**2.)*find_pk(k,itype,cosm)

    END FUNCTION p_lin

    FUNCTION p_2h(k,z,plin,lut,cosm)

    !Calculates the 2-halo term
    REAL :: p_2h
    REAL, INTENT(IN) :: k, plin
    REAL :: sigv, frac
    REAL, INTENT(IN) :: z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Damping function
    frac=fdamp(z,lut,cosm)

    IF(frac<1e-3) THEN
        p_2h=plin
    ELSE
        sigv=lut%sigv
        p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2.)
    END IF

    !For some strange cosmologies frac>1. so this must be added to prevent p_2h<0.
    IF(p_2h<0.) p_2h=0.

    END FUNCTION p_2h

    FUNCTION p_1h(k,z,lut,cosm)

    !Calculates the 1-halo term
    REAL :: p_1h
    REAL, INTENT(IN) :: k, z
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL :: Dv, g, fac, et, ks, wk
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: sum
    INTEGER :: i
    REAL, PARAMETER :: pi=3.141592654

    !Does the one-halo power integral

    !Allocates arrays for the integration HM_tables
    ALLOCATE(integrand(lut%n))
    integrand=0.

    !Only call eta once
    et=eta(z,lut,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
        g=gnu(lut%nu(i))
        wk=win(k*(lut%nu(i)**et),lut%rv(i),lut%c(i))
        integrand(i)=(lut%rv(i)**3.)*g*(wk**2.)
    END DO

    !Carries out the integration
    sum=inttab(lut%nu,REAL(integrand),1)

    !Deallocate arrays
    DEALLOCATE(integrand)

    !Virial density
    Dv=Delta_v(z,lut,cosm)

    !These are just numerical factors from the 1-halo integral in terms of nu!
    p_1h=sum*2.*Dv*(k**3.)/(3.*pi)

    !Damping of the 1-halo term at very large scales
    ks=kstar(z,lut,cosm)

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
    REAL :: rmin, rmax
    REAL :: r, sig
    INTEGER :: i, nsig
    TYPE(HM_cosmology) :: cosm

    !This fills up HM_tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigmac(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is approximately power-law below and above these values of R
    !This wouldn't be appropriate for models with a small-scale linear spectrum cut-off (e.g., WDM)

    !Allocate arrays
    IF(ALLOCATED(cosm%r_sigma)) DEALLOCATE(cosm%r_sigma)
    IF(ALLOCATED(cosm%sigma))   DEALLOCATE(cosm%sigma)

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    rmin=1e-4
    rmax=1e3
    nsig=64
    ALLOCATE(cosm%r_sigma(nsig),cosm%sigma(nsig))

    IF(ihm==1) WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
    IF(ihm==1) WRITE(*,*) 'SIGTAB: R_min:', rmin
    IF(ihm==1) WRITE(*,*) 'SIGTAB: R_max:', rmax
    IF(ihm==1) WRITE(*,*) 'SIGTAB: Values:', nsig

    !$OMP PARALLEL DO default(shared), private(sig, r)
    DO i=1,nsig

        !Equally spaced r in log
        r=exp(log(rmin)+log(rmax/rmin)*float(i-1)/float(nsig-1))
        sig=sigma(r,0.,1,cosm)

        cosm%r_sigma(i)=r
        cosm%sigma(i)=sig

    END DO
    !$OMP END PARALLEL DO

    IF(ihm==1) WRITE(*,*) 'SIGTAB: sigma_min:', cosm%sigma(nsig)
    IF(ihm==1) WRITE(*,*) 'SIGTAB: sigma_max:', cosm%sigma(1)
    IF(ihm==1) WRITE(*,*) 'SIGTAB: Done'
    IF(ihm==1) WRITE(*,*)

    END SUBROUTINE fill_sigtab

    FUNCTION sigma_v(R,z,cosm)

    !Calcuates the RMS in the displacement field at scale R
    REAL :: sigma_v
    REAL, INTENT(IN) :: z, R
    real(dl) :: sum
    REAL :: alpha
    REAL :: dtheta, k, theta, oldsum, acc
    REAL, PARAMETER :: pi=3.141592654
    INTEGER :: i, j, n, ninit, jmax
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    acc=0.001
    ninit=100
    jmax=30

    !Rapidising integral
    alpha=1.65

    DO j=1,jmax

        n=ninit*(2**(j-1))

        sum=0.d0
        dtheta=1./float(n)

        DO i=2,n-1

            !theta converts integral to 0->1 range
            !Values at the end points are 0 so removed for convenience
            theta=float(i-1)/float(n-1)
            k=(-1.+1./theta)/r**alpha
            sum=sum+p_lin(k,z,0,cosm)*(wk_tophat(k*r)**2.)/((k**2.)*theta*(1.-theta))

        END DO

        sum=sum*dtheta

        IF(j>1 .AND. ABS(-1.+sum/oldsum)<acc) THEN
            !Convert from sigma_v^2 and to 1D dispersion
            sigma_v=sqrt(sum/3.)
            EXIT
        ELSE
            oldsum=sum
        END IF

    END DO

    END FUNCTION sigma_v

    FUNCTION sigma(r,z,itype,cosm)

    !Chooses how best to do the sigma(R) integral depending on the R value
    REAL :: sigma
    REAL, INTENT(IN) :: r, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(r>=1.e-2) THEN
        sigma=sigint0(r,z,itype,cosm)
    ELSE IF(r<1.e-2) THEN
        sigma=sqrt(sigint1(r,z,itype,cosm)+sigint2(r,z,itype,cosm))
    END IF

    END FUNCTION sigma

    FUNCTION sigmac(r,z,cosm)

    !Finds sigma_cold(R) from look-up table
    REAL :: sigmac
    REAL, INTENT(IN) :: r, z
    REAL :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Uses the approximation sigma(R,z)=g(z)*sigma(R,z=0)

    !Approximate the effect of massive neutrinos using the small-scale limit of the ratio
    !between the perturbation of cold matter to all matter (assumes delta_nu=0)
    !This works remarkably well (error <1% in P(k) for standard LCDM for m_nu<0.3eV)

    sigmac=grow(z,cosm)*exp(find(log(r),log(cosm%r_sigma),log(cosm%sigma),3,3))!/(1.-cosm%f_nu)

    END FUNCTION sigmac

    FUNCTION wk_tophat(x)

    !The normlaised Fourier Transform of a top-hat
    Real :: wk_tophat, x

    !Taylor expansion used for low x to avoid cancellation problems
    IF(x<0.01) THEN
        wk_tophat=1.-(x**2.)/10.
    ELSE
        wk_tophat=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

    END FUNCTION wk_tophat

    FUNCTION inttab(x,y,iorder)

    !Routine to integrate HM_tables of data using the trapezium rule
    REAL :: inttab
    REAL, INTENT(IN) :: x(:), y(:)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    real(dl) :: sum
    INTEGER :: i, n, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    !Can either use linear, quadratic or cubic methods

    n=SIZE(x)

    IF(n .NE. SIZE(y)) ERROR STOP 'HM_tables must be of the same length'

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

    END IF

    inttab=sum

    END FUNCTION inttab

    FUNCTION sigma_integrand(t,R,f,z,itype,cosm)

    !The integrand for the sigma(R) intergal
    REAL :: sigma_integrand
    REAL, INTENT(IN) :: t, R, z
    REAL :: k, y, w_hat
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, INTENT(IN) :: itype

    INTERFACE
    REAL FUNCTION f(x)
    REAL, INTENT(IN) :: x
    END FUNCTION f
    END INTERFACE

    !Integrand to the sigma integral in terms of t. Defined by k=(1/t-1)/f(R) where f(R) is *any* function

    IF(t==0.) THEN
        !t=0 corresponds to k=infintiy when W(kR)=0.
        sigma_integrand=0.
    ELSE IF(t==1.) THEN
        !t=1 corresponds to k=0. when P(k)=0.
        sigma_integrand=0.
    ELSE
        !f(R) can be *any* function of R here to improve integration speed
        k=(-1.+1./t)/f(R)
        y=k*R
        w_hat=wk_tophat(y)
        sigma_integrand=p_lin(k,z,itype,cosm)*(w_hat**2.)/(t*(1.-t))
    END IF

    END FUNCTION sigma_integrand

    FUNCTION f_rapid(r)

    !Rapidizing function for the sigma(R) integrals
    REAL :: f_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error

    IF(r>1.e-2) THEN
        !alpha 0.3-0.5 works well
        alpha=0.5
    ELSE
        !alpha 0.7-0.9 works well
        alpha=0.8
    END IF

    f_rapid=r**alpha

    END FUNCTION f_rapid

    FUNCTION sigint0(r,z,itype,cosm)

    !One form of the sigma(R) integral
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint0, acc, dx
    INTEGER :: ninit, n
    REAL :: x
    real(dl) :: sum1, sum2
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, INTENT(IN) :: itype

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    DO j=1,jmax

        n=ninit*2**(j-1)

        !Avoids the end-points where the integrand is 0 anyway
        DO i=2,n-1

            !x is defined on the interval 0 -> 1
            x=float(i-1)/float(n-1)
            sum2=sum2+sigma_integrand(x,r,f_rapid,z,itype,cosm)

        END DO

        dx=1./float(n-1)
        sum2=sum2*dx
        sum2=sqrt(sum2)

        IF(j==1) THEN
            sum1=sum2
        ELSE IF(ABS(-1.+sum2/sum1)<acc) THEN
            sigint0=sum2
            EXIT
        ELSE IF(j==jmax) THEN
            WRITE(*,*)
            WRITE(*,*) 'SIGINT0: r:', r
            ERROR STOP 'SIGINT0: Integration timed out'
        ELSE
            sum1=sum2
            sum2=0.d0
        END IF

    END DO

    END FUNCTION sigint0

    FUNCTION sigint1(r,z,itype,cosm)

    !One form of the sigma(R) integral
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint1, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, k
    real(dl) :: sum1, sum2
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, INTENT(IN) :: itype

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    xmin=r/(r+r**.5)
    xmax=1.

    DO j=1,jmax

        n=ninit*2**(j-1)

        !Avoids the end-point where the integrand is 0 anyway
        DO i=1,n-1

            x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

            IF(i==1 .OR. i==n) THEN
                fac=0.5
            ELSE
                fac=1.
            END IF

            k=(-1.+1./x)/r**.5
            sum2=sum2+fac*p_lin(k,z,itype,cosm)*(wk_tophat(k*r)**2.)/(x*(1.-x))

        END DO

        dx=(xmax-xmin)/float(n-1)
        sum2=sum2*dx

        IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
            sigint1=sum2
            EXIT
        ELSE IF(j==jmax) THEN
            WRITE(*,*)
            WRITE(*,*) 'SIGINT1: r:', r
            ERROR STOP 'SIGINT1: Integration timed out'
        ELSE
            sum1=sum2
            sum2=0.d0
        END IF

    END DO

    END FUNCTION sigint1

    FUNCTION sigint2(r,z,itype,cosm)

    !One form of the sigma(R) integral
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint2, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, A
    real(dl) :: sum1, sum2
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, INTENT(IN) :: itype

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    !How far to go out in 1/r units for integral
    A=10.

    xmin=1./r
    xmax=A/r

    DO j=1,jmax

        n=ninit*2**(j-1)

        DO i=1,n

            x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

            IF(i==1 .OR. i==n) THEN
                fac=0.5
            ELSE
                fac=1.
            END IF

            !Integrate linearly in k for the rapidly oscillating part
            sum2=sum2+fac*p_lin(x,z,itype,cosm)*(wk_tophat(x*r)**2.)/x

        END DO

        dx=(xmax-xmin)/float(n-1)
        sum2=sum2*dx

        IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
            sigint2=sum2
            EXIT
        ELSE IF(j==jmax) THEN
            WRITE(*,*)
            WRITE(*,*) 'SIGINT2: r:', r
            ERROR STOP 'SIGINT2: Integration timed out'
        ELSE
            sum1=sum2
            sum2=0.d0
        END IF

    END DO

    END FUNCTION sigint2

    FUNCTION win(k,rv,c)

    !Selects the halo window function (k-space halo profile)
    REAL :: win, k, rv, c

    !Choose the NFW analytic form
    win=winnfw(k,rv,c)

    !Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    IF(win>1.) win=1.
    IF(win<0.) win=0.

    END FUNCTION win

    FUNCTION winnfw(k,rv,c)

    !The analytic Fourier Transform of the NFW profile
    REAL :: winnfw
    REAL, INTENT(IN) :: k, rv, c
    REAL :: si1, si2, ci1, ci2, ks
    REAL :: p1, p2, p3

    !Define the scale wavenumber
    ks=k*rv/c

    !Sine and cosine integrals
    si1=si(ks)
    si2=si((1.+c)*ks)
    ci1=ci(ks)
    ci2=ci((1.+c)*ks)

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
    REAL :: mass, c

    mass=log(1.+c)-c/(1.+c)

    END FUNCTION mass

    FUNCTION gnu(nu)

    !Select the mass function
    REAL :: gnu, nu

    !Sheth & Torman (1999)
    gnu=gst(nu)

    END FUNCTION gnu

    FUNCTION gst(nu)

    !Sheth & Tormen (1999) mass function!
    REAL :: nu, gst
    REAL :: p, a

    !Note I use nu=dc/sigma whereas ST (1999) use nu=(dc/sigma)^2
    !This accounts for the different pre-factor and slighly changed nu dependence
    !f(nu^2)d(nu^2)=2*nu*f(nu)dnu

    !Sheth & Tormen fitting numbers
    p=0.3
    a=0.707

    !Full mass function. Note this is normalised such that integral f(nu)dnu = 1
    gst=0.21616*(1.+((a*nu*nu)**(-p)))*exp(-a*nu*nu/2.)

    END FUNCTION gst

    FUNCTION Hubble2(z,cosm)

    !This calculates the dimensionless squared hubble parameter at redshift z (H/H_0)^2!
    !Ignores contributions from radiation (not accurate at high z, but consistent with simulations)!
    REAL :: Hubble2
    REAL, INTENT(IN) :: z
    REAL :: om_m, om_v, a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    om_m=cosm%om_m
    om_v=cosm%om_v

    a=1./(1.+z)

    Hubble2=(om_m*(1.+z)**3.)+om_v*x_de(a,cosm)+((1.-om_m-om_v)*(1.+z)**2.)

    END FUNCTION Hubble2

    FUNCTION x_de(a,cosm)

    !The time evolution for dark energy: rho_de=rho_de,0 * X(a)
    !X(a)=1 for LCDM but changes for other models
    REAL :: x_de
    REAL, INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    x_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

    END FUNCTION x_de

    FUNCTION wde(a,cosm)

    !The dark energy w(a) function
    IMPLICIT NONE
    REAL :: wde
    REAL, INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    wde=cosm%w+(1.-a)*cosm%wa

    END FUNCTION wde

    FUNCTION Omega_m_hm(z,cosm)

    !This calculates omega_m variations with z!
    REAL :: Omega_m_hm
    REAL, INTENT(IN) :: z
    REAL :: om_m, a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    om_m=cosm%om_m
    Omega_m_hm=(om_m*(1.+z)**3.)/Hubble2(z,cosm)

    END FUNCTION Omega_m_hm

    FUNCTION Omega_v_hm(z,cosm)

    !This calculates omega_v variations with z!
    REAL :: Omega_v_hm
    REAL, INTENT(IN) :: z
    REAL :: om_v, a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    om_v=cosm%om_v
    a=1./(1.+z)
    Omega_v_hm=om_v*x_de(a,cosm)/Hubble2(z,cosm)

    END FUNCTION Omega_v_hm

    FUNCTION grow(z,cosm)

    !Finds the scale-independent growth fuction at redshift z
    REAL :: grow
    REAL, INTENT(IN) :: z
    REAL :: a, acc
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(z==0.) THEN
        grow=1.
    ELSE
        a=1./(1.+z)
        grow=find(a,cosm%a_growth,cosm%growth,3,3)
    END IF

    END FUNCTION grow

    FUNCTION dispint(z,cosm)

    !Calculates the variance in the displacement field over all scales
    !This converges, unlike the same quantiy for the density field
    REAL :: dispint
    REAL, INTENT(IN) :: z
    real(dl) :: sum
    REAL :: dtheta, k, theta, oldsum, acc
    REAL, PARAMETER :: pi=3.141592654
    INTEGER :: i, j, n, ninit, jmax
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    acc=0.001
    ninit=100
    jmax=30

    DO j=1,jmax

        n=ninit*(2**(j-1))

        sum=0.d0
        dtheta=1./float(n)

        DO i=2,n-1

            !theta converts integral to 0->1 range
            !Values at the end points are 0 so removed for convenience
            theta=float(i-1)/float(n-1)
            k=(-1.+1./theta)
            sum=sum+((1.+k)**2.)*p_lin(k,z,0,cosm)/(k**3.)

        END DO

        sum=sum*dtheta

        IF(j>1 .AND. ABS(-1.+sum/oldsum)<acc) THEN
            !Division by 3 because we are interested in 1D, not 3D
            dispint=sum/3.
            EXIT
        ELSE
            oldsum=sum
        END IF

    END DO

    END FUNCTION dispint

    FUNCTION si(x)

    !Calculates the 'sine integral' function Si(x)
    REAL :: si, x
    real(dl) :: x2, y, f, g, si8
    real(dl), PARAMETER :: pi=3.1415926535897932384626433d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=4.) THEN

        x2=x*x

        si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

        si=si8

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

        si=pi/2.d0-f*cos(x)-g*sin(x)

    END IF

    END FUNCTION si

    FUNCTION ci(x)

    !Calculates the 'cosine integral' function Ci(x)
    REAL :: ci, x
    real(dl) :: x2, y, f, g, ci8
    real(dl), PARAMETER :: em_const=0.577215664901532861d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=4.) THEN

        x2=x*x

        ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

        ci=ci8

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

        ci=f*sin(x)-g*cos(x)

    END IF

    END FUNCTION ci

    FUNCTION derivative_table(x,xin,yin,iorder,imeth)

    !Calculates the derivative of a tabulated function (xin,yin) at point x
    REAL :: derivative_table
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !Finds the derivative f'(x) given HM_tables x, f(x)

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !imeth = 1 => find x in xtab by crudely searching
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    n=SIZE(xtab)

    maxorder=3
    maxmethod=3

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) ERROR STOP 'FIND: HM_tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
        !Reverse the arrays in this case
        CALL reverse(xtab)
        CALL reverse(ytab)
    END IF

    IF(iorder<1) ERROR STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) ERROR STOP 'FIND: find order not specified correctly'
    IF(imeth<1) ERROR STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) ERROR STOP 'FIND: Method of finding within a table not specified correctly'

    IF(iorder==1) THEN

        IF(n<2) ERROR STOP 'FIND: Not enough points in your table for linear interpolation'

        IF(x<=xtab(2)) THEN

            x2=xtab(2)
            x1=xtab(1)

            y2=ytab(2)
            y1=ytab(1)

        ELSE IF (x>=xtab(n-1)) THEN

            x2=xtab(n)
            x1=xtab(n-1)

            y2=ytab(n)
            y1=ytab(n-1)

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x2=xtab(i+1)
            x1=xtab(i)

            y2=ytab(i+1)
            y1=ytab(i)

        END IF

        CALL fit_line(a,b,x1,y1,x2,y2)
        derivative_table=a

    ELSE IF(iorder==2) THEN

        IF(n<3) ERROR STOP 'FIND_QUADRATIC: Not enough points in your table'

        IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

            IF(x<=xtab(2)) THEN

                x3=xtab(3)
                x2=xtab(2)
                x1=xtab(1)

                y3=ytab(3)
                y2=ytab(2)
                y1=ytab(1)

            ELSE IF (x>=xtab(n-1)) THEN

                x3=xtab(n)
                x2=xtab(n-1)
                x1=xtab(n-2)

                y3=ytab(n)
                y2=ytab(n-1)
                y1=ytab(n-2)

            END IF

            CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

            derivative_table=2.*a*x+b

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x1=xtab(i-1)
            x2=xtab(i)
            x3=xtab(i+1)
            x4=xtab(i+2)

            y1=ytab(i-1)
            y2=ytab(i)
            y3=ytab(i+1)
            y4=ytab(i+2)

            !In this case take the average of two separate quadratic spline values

            derivative_table=0.

            CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
            derivative_table=derivative_table+(2.*a*x+b)/2.

            CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
            derivative_table=derivative_table+(2.*a*x+b)/2.

        END IF

    ELSE IF(iorder==3) THEN

        IF(n<4) ERROR STOP 'FIND_CUBIC: Not enough points in your table'

        IF(x<=xtab(3)) THEN

            x4=xtab(4)
            x3=xtab(3)
            x2=xtab(2)
            x1=xtab(1)

            y4=ytab(4)
            y3=ytab(3)
            y2=ytab(2)
            y1=ytab(1)

        ELSE IF (x>=xtab(n-2)) THEN

            x4=xtab(n)
            x3=xtab(n-1)
            x2=xtab(n-2)
            x1=xtab(n-3)

            y4=ytab(n)
            y3=ytab(n-1)
            y2=ytab(n-2)
            y1=ytab(n-3)

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x1=xtab(i-1)
            x2=xtab(i)
            x3=xtab(i+1)
            x4=xtab(i+2)

            y1=ytab(i-1)
            y2=ytab(i)
            y3=ytab(i+1)
            y4=ytab(i+2)

        END IF

        CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
        derivative_table=3.*a*(x**2.)+2.*b*x+c

    END IF

    END FUNCTION derivative_table

    FUNCTION find(x,xin,yin,iorder,imeth)

    !Interpolation routine to get the value y(x) from HM_tables of x and y
    REAL :: find
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !Interpolation routine.

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    maxorder=3
    maxmethod=3

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) ERROR STOP 'FIND: HM_tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
        !Reverse the arrays in this case
        CALL reverse(xtab)
        CALL reverse(ytab)
    END IF

    IF(iorder<1) ERROR STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) ERROR STOP 'FIND: find order not specified correctly'
    IF(imeth<1) ERROR STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) ERROR STOP 'FIND: Method of finding within a table not specified correctly'

    IF(x<xtab(1)) THEN

        x1=xtab(1)
        x2=xtab(2)

        y1=ytab(1)
        y2=ytab(2)

        CALL fit_line(a,b,x1,y1,x2,y2)
        find=a*x+b

    ELSE IF(x>xtab(n)) THEN

        x1=xtab(n-1)
        x2=xtab(n)

        y1=ytab(n-1)
        y2=ytab(n)

        CALL fit_line(a,b,x1,y1,x2,y2)
        find=a*x+b

    ELSE IF(iorder==1) THEN

        IF(n<2) ERROR STOP 'FIND: Not enough points in your table for linear interpolation'

        IF(x<=xtab(2)) THEN

            x1=xtab(1)
            x2=xtab(2)

            y1=ytab(1)
            y2=ytab(2)

        ELSE IF (x>=xtab(n-1)) THEN

            x1=xtab(n-1)
            x2=xtab(n)

            y1=ytab(n-1)
            y2=ytab(n)

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x1=xtab(i)
            x2=xtab(i+1)

            y1=ytab(i)
            y2=ytab(i+1)

        END IF

        CALL fit_line(a,b,x1,y1,x2,y2)
        find=a*x+b

    ELSE IF(iorder==2) THEN

        IF(n<3) ERROR STOP 'FIND: Not enough points in your table'

        IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

            IF(x<=xtab(2)) THEN

                x1=xtab(1)
                x2=xtab(2)
                x3=xtab(3)

                y1=ytab(1)
                y2=ytab(2)
                y3=ytab(3)

            ELSE IF (x>=xtab(n-1)) THEN

                x1=xtab(n-2)
                x2=xtab(n-1)
                x3=xtab(n)

                y1=ytab(n-2)
                y2=ytab(n-1)
                y3=ytab(n)

            END IF

            CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

            find=a*(x**2.)+b*x+c

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x1=xtab(i-1)
            x2=xtab(i)
            x3=xtab(i+1)
            x4=xtab(i+2)

            y1=ytab(i-1)
            y2=ytab(i)
            y3=ytab(i+1)
            y4=ytab(i+2)

            !In this case take the average of two separate quadratic spline values

            find=0.

            CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
            find=find+(a*(x**2.)+b*x+c)/2.

            CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
            find=find+(a*(x**2.)+b*x+c)/2.

        END IF

    ELSE IF(iorder==3) THEN

        IF(n<4) ERROR STOP 'FIND: Not enough points in your table'

        IF(x<=xtab(3)) THEN

            x1=xtab(1)
            x2=xtab(2)
            x3=xtab(3)
            x4=xtab(4)

            y1=ytab(1)
            y2=ytab(2)
            y3=ytab(3)
            y4=ytab(4)

        ELSE IF (x>=xtab(n-2)) THEN

            x1=xtab(n-3)
            x2=xtab(n-2)
            x3=xtab(n-1)
            x4=xtab(n)

            y1=ytab(n-3)
            y2=ytab(n-2)
            y3=ytab(n-1)
            y4=ytab(n)

        ELSE

            IF(imeth==1) i=search_int(x,xtab)
            IF(imeth==2) i=linear_table_integer(x,xtab)
            IF(imeth==3) i=int_split(x,xtab)

            x1=xtab(i-1)
            x2=xtab(i)
            x3=xtab(i+1)
            x4=xtab(i+2)

            y1=ytab(i-1)
            y2=ytab(i)
            y3=ytab(i+1)
            y4=ytab(i+2)

        END IF

        CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
        find=a*x**3.+b*x**2.+c*x+d

    END IF

    END FUNCTION find

    FUNCTION linear_table_integer(x,xtab)

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 3
    !Assumes table is organised linearly (care for logs)
    INTEGER :: linear_table_integer
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: n
    REAL :: x1, x2, xn
    REAL :: acc

    n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    !Test for linear table
    acc=0.001

    IF(x1>xn) ERROR STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) ERROR STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

    END FUNCTION linear_table_integer

    FUNCTION search_int(x,xtab)

    !Finds the integer table position below x in the table xtab using a silly brute force method
    INTEGER :: search_int
    INTEGER :: i, n
    REAL, INTENT(IN) :: x, xtab(:)

    !Searches for the point in the table brute force.
    !This is usually a stupid thing to do

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) ERROR STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
        IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

    END FUNCTION search_int

    FUNCTION int_split(x,xtab)

    !Finds the integer table position by continually splitting it the table in half
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: i1, i2, imid, n
    INTEGER :: int_split

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) ERROR STOP 'INT_SPLIT: table in wrong order'

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
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

    END SUBROUTINE fit_line

    SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    !Given xi, yi i=1,2,3 fits a quadratic between these points
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

    END SUBROUTINE fit_quadratic

    SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    !Given xi, yi i=1,2,3,4 fits a cubic between these points
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

    END SUBROUTINE fit_cubic

    SUBROUTINE reverse(arry)

    !This reverses the contents (order) of the array 'arry'
    INTEGER :: n, i
    REAL, ALLOCATABLE :: hold(:)
    REAL :: arry(:)

    n=SIZE(arry)

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
    INTEGER :: i, n
    REAL :: a, norm
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit
    REAL :: acc

    !The calculation should start at a z when Om_m(z)=1., so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    ainit=0.001
    !Final should be a=1. unless considering models in the future
    amax=1.

    !These set the initial conditions to be the Om_m=1. growing mode
    dinit=ainit
    vinit=1.

    !Overall accuracy for the ODE solver
    acc=0.001

    IF(ihm==1) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,0.,ainit,amax,dinit,vinit,acc,3,cosm)
    IF(ihm==1) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    norm=find(1.,a_tab,d_tab,3,3)
    IF(ihm==1) WRITE(*,*) 'GROWTH: Unnormalised g(a=1):', norm
    d_tab=d_tab/norm

    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_growth)) DEALLOCATE(cosm%a_growth)
    IF(ALLOCATED(cosm%growth)) DEALLOCATE(cosm%growth)
    n=64
    ALLOCATE(cosm%a_growth(n),cosm%growth(n))
    DO i=1,n
        a=ainit+(amax-ainit)*float(i-1)/float(n-1)
        cosm%a_growth(i)=a
        cosm%growth(i)=find(a,a_tab,d_tab,3,3)
    END DO

    IF(ihm==1) WRITE(*,*) 'GROWTH: Done'
    IF(ihm==1) WRITE(*,*)

    END SUBROUTINE fill_growtab

    SUBROUTINE ode_growth(x,v,t,kk,ti,tf,xi,vi,acc,imeth,cosm)

    !ODE solver used for growth equations
    IMPLICIT NONE
    REAL :: xi, ti, tf, dt, acc, vi, x4, v4, t4, kk
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    real(dl), ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL, ALLOCATABLE :: x(:), v(:), t(:)
    INTEGER :: i, j, k, n, np, ifail, kn, imeth
    TYPE(HM_cosmology) :: cosm

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)

    DO j=1,30

        n=100*(2**(j-1))
        n=n+1

        ALLOCATE(x8(n),t8(n),v8(n))

        x8=0.d0
        t8=0.d0
        v8=0.d0

        dt=(tf-ti)/float(n-1)

        x8(1)=xi
        v8(1)=vi
        t8(1)=ti

        ifail=0

        DO i=1,n-1

            x4=real(x8(i))
            v4=real(v8(i))
            t4=real(t8(i))

            IF(imeth==1) THEN

                !Crude method!
                v8(i+1)=v8(i)+fv(x4,v4,kk,t4,cosm)*dt
                x8(i+1)=x8(i)+fd(x4,v4,kk,t4,cosm)*dt
                t8(i+1)=t8(i)+dt

            ELSE IF(imeth==2) THEN

                !Mid-point method!
                kx1=dt*fd(x4,v4,kk,t4,cosm)
                kv1=dt*fv(x4,v4,kk,t4,cosm)
                kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
                kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)

                v8(i+1)=v8(i)+kv2*dt
                x8(i+1)=x8(i)+kx2*dt
                t8(i+1)=t8(i)+dt

            ELSE IF(imeth==3) THEN

                !RK4 (Holy Christ, this is so fast compared to above methods)!
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
                t8(i+1)=t8(i)+dt

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
            x=x8
            v=v8
            t=t8
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
    REAL :: fv
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: f1, f2, z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    z=-1.+(1./a)

    f1=3.*omega_m_hm(z,cosm)*d/(2.*(a**2.))
    f2=(2.+AH(z,cosm)/hubble2(z,cosm))*(v/a)

    fv=f1-f2

    END FUNCTION fv

    FUNCTION fd(d,v,k,a,cosm)

    !d'=f(d) in ODE solver
    REAL :: fd
    REAL, INTENT(IN) :: d, v, k, a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    fd=v

    END FUNCTION fd

    FUNCTION AH(z,cosm)

    !The Hubble acceleration function \ddot{a}/a
    REAL :: AH
    REAL, INTENT(IN) :: z
    REAL :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    a=1./(1.+z)
    AH=cosm%om_m*(a**(-3.))+cosm%om_v*(1.+3.*wde(a,cosm))*x_de(a,cosm)
    AH=-AH/2.

    END FUNCTION AH

    !!AM End HMcode

    end module NonLinear


    !workaround for f90 circular-module reference
    subroutine NonLinear_GetRatios(CAMB_Pk)
    use Transfer
    use NonLinear
    type(MatterPowerData) :: CAMB_Pk

    call NonLinear_GetNonLinRatios(CAMB_Pk)

    end subroutine NonLinear_GetRatios



    subroutine NonLinear_GetRatios_all(CAMB_Pk)
    use Transfer
    use NonLinear
    type(MatterPowerData) :: CAMB_Pk

    call MpiStop('Halofit module doesn''t support non-linear velocities')

    end subroutine NonLinear_GetRatios_All

