    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The `halofit' code models the nonlinear evolution of cold matter
    ! cosmological power spectra. The full details of the way in which
    ! this is done are presented in Smith et al. (2003), MNRAS, 341, 4
    !
    ! The code `halofit' was written by R. E. Smith & J. A. Peacock.
    ! See http://www.astro.upenn.edu/~res
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
    !AL Jul 19: Speedups, use linear interpolation for pk; find index using fixed spacing; precompute growth(z)
    !AL Sep 19: Propagate errors rather than stop, decrease jmax for integration time out (prevent very slow error)
    !AM Sep 20: Added HMcode-2020 model
    !AM Jan 23: Fixed HMcode-2020 feedback low-k predictions
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module NonLinear
    use results
    use DarkEnergyInterface
    use classes
    use Transfer
    use constants
    use config
    implicit none
    private

    integer, parameter :: halofit_original=1, halofit_bird=2, halofit_peacock=3, halofit_takahashi=4
    integer, parameter :: halofit_casarini=7
    integer, parameter :: halofit_mead2016=5, halofit_halomodel=6, halofit_mead2015=8, halofit_mead2020=9
    integer, parameter :: halofit_mead2020_feedback=10
    integer, parameter :: halofit_mead=halofit_mead2016 ! AM Kept for backwards compatability
    integer, parameter :: halofit_default=halofit_mead2020

    logical :: HM_verbose = .false.

    type, extends(TNonLinearModel) :: THalofit
        integer :: halofit_version = halofit_default
        !!TT - These are the baryon parameters of HMCode
        real(dl) :: HMcode_A_baryon=3.13_dl
        real(dl) :: HMcode_eta_baryon=0.603_dl
        real(dl) :: HMcode_logT_AGN=7.8_dl
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
    public halofit_mead2016, halofit_mead2015, halofit_mead2020, halofit_halomodel, halofit_casarini
    public halofit_mead2020_feedback
    public halofit_mead ! AM for backwards compatability

    TYPE HM_cosmology
        !Contains only things that do not need to be recalculated with each new z
        REAL(dl) :: om_m, om_c, om_b, om_nu, om_v, w, wa, f_nu, ns, h, Tcmb, Nnu
        REAL(dl), ALLOCATABLE :: log_r_sigma(:), log_sigma(:)
        REAL(dl), ALLOCATABLE :: a_growth(:), growth(:), agrow(:)
        REAL(dl), ALLOCATABLE :: log_k_plin(:), log_plin(:), log_plinc(:)
        REAL(dl), ALLOCATABLE :: log_k_wiggle(:), pk_wiggle(:)
        real(dl) :: kmax
        real(dl) :: gnorm
        INTEGER :: nk, ng, nsig
        real(dl) :: grow2_z, this_z !cached value at redshift being calculated
        !AM - Added feedback parameters below at fixed fiducial (DMONLY) values
        REAL(dl) :: A_baryon=3.13
        REAL(dl) :: eta_baryon=0.603
        REAL(dl) :: logT_AGN=7.8
    END TYPE HM_cosmology

    TYPE HM_tables
        !Stuff that needs to be recalculated for each new z
        REAL(dl), ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
        REAL(dl) :: sigv, sigv100, knl, rnl, neff, sig8z, z, dc, sig8z_cold
        INTEGER :: n
    END TYPE HM_tables
    !!AM - End of my additions

    ! HMcode parameters
    REAL(dl), PARAMETER :: zc_Dolag=10._dl   ! Halo collapse redshift for Dolag
    REAL(dl), PARAMETER :: fdamp_min=1e-3_dl ! Minimum value of fdamp
    REAL(dl), PARAMETER :: fdamp_max=0.99_dl ! Maximum value of fdamp
    REAL(dl), PARAMETER :: alpha_min=0.5_dl  ! Minimum value of alpha transition
    REAL(dl), PARAMETER :: alpha_max=2._dl   ! Maximum value of alpha transition
    REAL(dl), PARAMETER :: ks_limit=7._dl    ! Limit for (k/ks)^2 in one-halo term
    REAL(dl), PARAMETER :: pi_HM=const_pi    ! Lovely pi

    ! HMcode linear P(k) numerical parameters
    ! AM: Jul 19: Updated nk_pk_interpolation from 128 to 512
    ! AM: Dec 20: Calculation time and accuracy are especially sensive to these parameters
    LOGICAL, PARAMETER :: rebin_pk=.TRUE.             ! Should the linear P(k) be rebinned?
    REAL(dl), PARAMETER :: kmin_pk_interpolation=1d-3 ! Minimum wavenumber if rebinning [h/Mpc]
    REAL(dl), PARAMETER :: kmax_pk_interpolation=1d2  ! Maximum wavenumber if rebinning [h/Mpc]
    INTEGER, PARAMETER :: nk_pk_interpolation=512     ! Number of points in k if rebining
    LOGICAL, PARAMETER :: plin_extrap=.FALSE.         ! Extrapolate at high-k via thoery or simple power law
    INTEGER, PARAMETER :: iorder_pk_interpolation=1   ! Polynomial order for P(k) interpolation
    INTEGER, PARAMETER :: ifind_pk_interpolation=1    ! Finding scheme for P(k) interpolation (if rebin_pk=True)
    INTEGER, PARAMETER :: imeth_pk_interpolation=1    ! Method for P(k) interpolation

    ! HMcode dewiggle numerical parameters
    REAL, PARAMETER :: kmin_wiggle=5e-3    ! Minimum wavenumber to calulate wiggle [Mpc/h]
    REAL, PARAMETER :: kmax_wiggle=5.      ! Maximum wavenumber to calulate wiggle [Mpc/h]
    INTEGER, PARAMETER :: nk_wiggle=512    ! Number of k points to store wiggle
    INTEGER, PARAMETER :: iorder_wiggle=3  ! Order for wiggle interpolation
    INTEGER, PARAMETER :: ifind_wiggle=3   ! 3 - Mid-point finding scheme for wiggle interpolation
    INTEGER, PARAMETER :: imeth_wiggle=2   ! 2- Lagrange polynomial interpolation
    REAL, PARAMETER :: wiggle_sigma=0.25   ! Smoothing width if using Gaussian smoothing
    REAL, PARAMETER :: knorm_nowiggle=0.03 ! Wavenumber at which to force linear and nowiggle to be identical [Mpc/h]

    ! Linear growth integral numerical parameters (LCDM only; only used in Dolag correction)
    ! AM: Jul 19: Updated acc_growint from 1e-3 to 1e-4
    ! AM: Sep 20: Changed cold_growth = .FALSE. to be in line with my code
    REAL(dl), PARAMETER :: acc_growth_integration=1e-4 ! Accuracy for growth function integral
    INTEGER, PARAMETER :: iorder_growth_integration=3  ! Polynomial order for growth integral
    LOGICAL, PARAMETER :: cold_growth=.FALSE.          ! Should growth be of cold or all matter?

    ! Linear growth factor tabulation and interpolation numerical parameters
    ! AM: TODO: Change finding scheme to assume linear spacing may save time
    REAL(dl), PARAMETER :: amin_growth_interpolation=1e-3 ! Minimum scale factor for growth interpolation
    REAL(dl), PARAMETER :: amax_growth_interpolation=1.   ! Maximum scale factor for growth interpolation
    INTEGER, PARAMETER :: n_growth_interpolation=64       ! Number of entries for growth look-up table
    INTEGER, PARAMETER :: iorder_growth_interpolation=3   ! Polynomial order for growth function interpolation
    INTEGER, PARAMETER :: ifind_growth_interpolation=3    ! Finding scheme for growth function interpolation
    INTEGER, PARAMETER :: imeth_growth_interpolation=2    ! Method for growth function interpolation

    ! Growth function ODE numerical parameters
    ! AM: Jul 19: Updated acc_growth_ODE from 1e-3 to 1e-4
    ! AM: Sep 20: Changed aini from 1e-3 to 1e-4
    REAL(dl), PARAMETER :: aini_growth_ODE=1e-4             ! Initial scale factor for growth ODE
    REAL(dl), PARAMETER :: afin_growth_ODE=1.               ! Final scale factor for growth ODE
    REAL(dl), PARAMETER :: acc_growth_ODE=1e-4              ! Accuracy for growth integral or ODE
    INTEGER, PARAMETER :: imeth_growth_ODE=3                ! Method for growth function ODE solving
    INTEGER, PARAMETER :: iorder_growth_ODE_interpolation=3 ! Polynomial order for growth function ODE interpolation
    INTEGER, PARAMETER :: ifind_growth_ODE_interpolation=3  ! Finding scheme for growth function ODE interpolation
    INTEGER, PARAMETER :: imeth_growth_ODE_interpolation=2  ! Method growth function ODE interpolation

    ! Linear growth function inversion numerical parameters (used for c(M) only)
    INTEGER, PARAMETER :: iorder_growth_inversion=3 ! Polynomial order for growth function ODE inversion
    INTEGER, PARAMETER :: ifind_growth_inversion=3  ! Finding scheme for growth function ODE inversion
    INTEGER, PARAMETER :: imeth_growth_inversion=2  ! Method growth function ODE inversion

    ! Accumulated growth parameters
    INTEGER, PARAMETER :: iorder_integration_agrow=3     ! Polynomial order for accumulated growth integration
    INTEGER, PARAMETER :: iorder_agrowth_interpolation=3 ! Polynomial order for accumlated growth interpolation
    INTEGER, PARAMETER :: ifind_agrowth_interpolation=3  ! Finding scheme for accumulated growth interpolation
    INTEGER, PARAMETER :: imeth_agrowth_interpolation=2  ! Method for accumulated growth interpolation

    ! HMcode numerical parameters for sigma(R) tabulation and interpolation
    REAL(dl), PARAMETER :: rmin_sigma_interpolation=1e-4 ! Minimum scale for sigma(R) look-up tables [Mpc/h]
    REAL(dl), PARAMETER :: rmax_sigma_interpolation=1e3  ! Maximum scale for sigma(R) look-up tables [Mpc/h]
    INTEGER, PARAMETER :: n_sigma_interpolation=64       ! Number of points in look-up tables
    INTEGER, PARAMETER :: iorder_sigma_interpolation=3   ! Polynomial order for sigma(R) interpolation
    INTEGER, PARAMETER :: ifind_sigma_interpolation=3    ! Finding scheme for sigma(R) interpolation
    INTEGER, PARAMETER :: imeth_sigma_interpolation=2    ! Method sigma(R) interpolation

    ! HMcode numerical parameters for sigma(R) integration (dominates run time as of Jul 2019)
    ! AM: Jul 19: Upgraded acc_sigma from 1e-3 to 3e-4
    ! AM: Sep 20: Upgraded acc_sigma from 3e-4 to 1e-4 to fix problems for some cosmologies
    REAL(dl), PARAMETER :: acc_sigma_integration=1e-4 ! Relative accuracy of numerical integration
    REAL(dl), PARAMETER :: alpha_sigma_integration=3. ! Exponent to speed up integration
    INTEGER, PARAMETER :: iorder_sigma_integration=3  ! Polynomail order for numerical integration

    ! HMcode numerical parameters for sigmaV(R) integration
    ! AM: Jul 19: Upgraded acc_sigmaV from 1e-3 to 1e-4
    REAL(dl), PARAMETER :: acc_sigmaV_integration=1e-4 ! Relative accuracy of numerical integration
    REAL(dl), PARAMETER :: alpha_sigmaV_integration=3. ! Exponent to speed up integration
    INTEGER, PARAMETER :: iorder_sigmaV_integration=3  ! Polynomial order for numerical integration

    ! HMcode numerical parameters for neff(R) integration
    REAL(dl), PARAMETER :: acc_neff_integration=1e-4 ! Relative accuracy of numerical integration
    REAL(dl), PARAMETER :: alpha_neff_integration=2. ! Exponent to speed up integration
    INTEGER, PARAMETER :: iorder_neff_integration=3  ! Polynomial order for numerical integration

    ! HMcode numerical parameters for cold transfer function approximation
    ! AM: Sep 20: Care here, before EdS_Tcold_growth=.TRUE.
    LOGICAL, PARAMETER :: EdS_Tcold_growth=.FALSE. ! Should the EdS growth function (incorrectly) be used?

    ! HMcode numerical parameters for one-halo term
    INTEGER, PARAMETER :: iorder_integration_1h=1 ! Should be linear order (i.e., trapezium rule)

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
    IF(this%halofit_version == halofit_mead2020_feedback) THEN
        this%HMcode_logT_AGN = Ini%Read_Double('HMcode_logT_AGN', 7.8_dl)
    END IF

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

    !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)

    select type (State)
    class is (CAMBdata)
        associate(Params => State%CP)

            IF(this%halofit_version==halofit_mead2016 .OR. &
                this%halofit_version==halofit_halomodel .OR. &
                this%halofit_version==halofit_mead2015 .OR. &
                this%halofit_version==halofit_mead2020 .OR. &
                this%halofit_version==halofit_mead2020_feedback) THEN
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
                            call GlobalError('Error in halofit (xlogr1>3.4999)', error_nonlinear)
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
    integer index_cache

    index_cache = 1
    nint=3000
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    anorm = 1/(2*const_pi**2)
    do i=1,nint
        t=(i-0.5_dl)/nint
        y=-1.d0+1.d0/t
        rk=y
        d2=MatterPowerData_k(CAMB_PK, rk, itf, index_cache)*(rk**3*anorm)
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
    REAL(dl), ALLOCATABLE :: p_den(:,:), p_num(:,:)
    INTEGER :: i, j, ii, nk, nz
    REAL :: t1, t2
    TYPE(HM_cosmology) :: cosi
    TYPE(HM_tables) :: lut
    REAL(dl), PARAMETER :: pi=pi_HM
    LOGICAL, PARAMETER :: timing_test = .FALSE.

    !HMcode developed by Alexander Mead (alexander.j.mead@googlemail.com)
    !Please contact me if you have any questions whatsoever
    !If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
    !If you use the extensions (w(a) and massive neutrinos) then please cite: http://arxiv.org/abs/1602.02154
    !Also consider citing the source code at ASCL: http://ascl.net/1508.001

    IF (timing_test) CALL CPU_TIME(t1)

    !Use imead to switch between the standard and accurate halo-model calcuation
    !0 - Standard (this is just a vanilla halo model calculation with no accuracy tweaks)
    !1 - Accurate from Mead et al. (2016; arXiv 1602.02154)
    !2 - Accurate from Mead et al. (2015; arXiv 1505.07833)
    !3 - Accurate from Mead et al. (2020; arXiv 2009.01858)
    !4 - Denominator for feedback reaction model from Mead et al. (2020; arXiv 2009.01858)
    !5 - Numerator for feedback reaction from Mead et al. (2020; arXiv 2009.01858)
    IF(this%halofit_version==halofit_halomodel) this%imead=0
    IF(this%halofit_version==halofit_mead2016) this%imead=1
    IF(this%halofit_version==halofit_mead2015) this%imead=2
    IF(this%halofit_version==halofit_mead2020) this%imead=3

    HM_verbose = (FeedbackLevel>1)

    IF(HM_verbose) WRITE(*,*)
    IF(HM_verbose) WRITE(*,*) 'HMcode: Running HMcode'
    IF(HM_verbose) WRITE(*,*)

    !!AM - Translate from CAMB variables to my variables
    nz=CAMB_PK%num_z
    nk=CAMB_PK%num_k
    IF(this%halofit_version==halofit_mead2020_feedback) THEN
        ALLOCATE(p_den(nk,nz), p_num(nk,nz))
    END IF

    !!AM - Assign cosmological parameters for the halo model calculation
    CALL assign_HM_cosmology(this,State,cosi)

    !Fill growth function table (only needs to be done once)
    CALL fill_growtab(cosi)

    !Loop over redshifts
    DO j=1,nz

        !Initialise the specific HM_cosmology (fill sigma(R) and P_lin HM_tables)
        !Currently this needs to be done at each z (mainly because of scale-dependent growth with neutrinos)
        !For non-massive-neutrino models this could only be done once, which would speed things up a bit
        CALL initialise_HM_cosmology(this,j,cosi,CAMB_PK)

        !Sets the current redshift from the table
        z=CAMB_Pk%Redshifts(j)

        IF(this%halofit_version==halofit_mead2020_feedback) THEN

            ! Loop over numerator, denominator and HMcode to make feedback response model
            DO ii = 1, 3

                IF(ii==1) this%imead=3 ! HMcode 2020
                IF(ii==2) this%imead=4 ! Denominator for response
                IF(ii==3) this%imead=5 ! Numerator for response

                !Initiliasation for the halomodel calculation (needs to be done for each z)
                CALL this%halomod_init(z,lut,cosi)
                if (global_error_flag/=0) return

                !Loop over k values and calculate P(k)
                !$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin,pfull,p1h,p2h)
                DO i=1,nk
                    k=exp(CAMB_Pk%log_kh(i))
                    plin=p_lin(k,z,0,cosi)
                    CALL this%halomod(k,p1h,p2h,pfull,plin,lut,cosi)
                    IF(this%imead==3) THEN
                        CAMB_Pk%nonlin_ratio(i,j)=sqrt(pfull/plin)
                    ELSE IF(this%imead==4) THEN
                        p_den(i,j)=pfull
                    ELSE IF(this%imead==5) THEN
                        p_num(i,j)=pfull
                    END IF
                END DO
                !$OMP END PARALLEL DO

            END DO

        ELSE

            !Initiliasation for the halomodel calculation (needs to be done for each z)
            CALL this%halomod_init(z,lut,cosi)
            if (global_error_flag/=0) return

            !Loop over k values and calculate P(k)
            !$OMP PARALLEL DO DEFAULT(SHARED), private(k,plin,pfull,p1h,p2h)
            DO i=1,nk
                k=exp(CAMB_Pk%log_kh(i))
                plin=p_lin(k,z,0,cosi)
                CALL this%halomod(k,p1h,p2h,pfull,plin,lut,cosi)
                CAMB_Pk%nonlin_ratio(i,j)=sqrt(pfull/plin)
            END DO
            !$OMP END PARALLEL DO

        END IF

    END DO

    ! Make the non-linear correction from the response for HMcode 2020
    IF(this%halofit_version==halofit_mead2020_feedback) THEN
        CAMB_Pk%nonlin_ratio=CAMB_Pk%nonlin_ratio*sqrt(p_num/p_den)
    END IF

    IF (timing_test) THEN
        CALL CPU_TIME(t2)
        WRITE(*, *) 'HMcode number of k:', nk
        WRITE(*, *) 'HMcode number of z:', nz
        WRITE(*, *) 'HMcode run time [s]:', t2-t1
        STOP 'HMcode timing test complete'
    END IF

    END SUBROUTINE HMcode

    FUNCTION Delta_v(this,z,cosm)
    class(THalofit) :: this
    !Function for the virialised overdensity
    REAL(dl) :: Delta_v
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    IF(this%imead==1 .OR. this%imead==2) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        Delta_v=418*Omega_m_hm(z,cosm)**(-0.352_dl)
        !Mead et al. (2016; arXiv 1602.02154) neutrino addition
        IF(this%imead==1) Delta_v=Delta_v*(1+0.916_dl*cosm%f_nu)
    ELSE IF(this%imead==0 .OR. this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
        Delta_v=Dv_Mead(z, cosm)
    END IF

    END FUNCTION Delta_v

    FUNCTION delta_c(this,z,lut,cosm)
    class(THalofit) :: this
    !Function for the linear collapse density
    REAL(dl) :: delta_c
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==1 .or. this%imead==2) THEN
        !Mead et al. (2015; arXiv 1505.07833) value
        delta_c=1.59+0.0314*log(lut%sig8z)
        IF(this%imead==1) THEN
            delta_c=delta_c*(1.+0.262*cosm%f_nu) !Mead et al. (2016; arXiv 1602.02154) neutrino addition
            delta_c=delta_c*(1.+0.0123*log10(Omega_m_hm(z,cosm))) !Nakamura & Suto (1997) fitting formula for LCDM
        END IF
    ELSE IF(this%imead==0 .OR. this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
        delta_c=dc_Mead(z, cosm)
    END IF

    END FUNCTION delta_c

    FUNCTION eta(this,lut,cosm)
    class(THalofit) :: this
    !Function eta that puffs halo profiles
    REAL(dl) :: eta
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL(dl) :: eta0

    IF(this%imead==0 .OR. this%imead==4 .OR. this%imead==5) THEN
        eta=0.
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
        !eta=0.603-0.3*lut%sig8z
        !AM - made baryon feedback parameter obvious
        eta0=cosm%eta_baryon
        !eta0=1.03-0.11*cosm%A_baryon !Original one-parameter relation from 1505.07833
        !eta0=0.98-0.12*cosm%A_baryon !Updated one-parameter relation: Section 4.1.2 of 1707.06627
        eta=eta0-0.3*lut%sig8z
    ELSE IF(this%imead==3) THEN
        eta=0.1281*lut%sig8z_cold**(-0.3644)
    END IF

    END FUNCTION eta

    FUNCTION kstar(this,lut)
    class(THalofit) :: this
    !Function k* that cuts off the 1-halo term at large scales
    REAL(dl) :: kstar
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0) THEN
        !Set to zero for the standard Poisson one-halo term
        kstar=0.
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !One-halo cut-off wavenumber
        !Mead et al. (2015; arXiv 1505.07833) value
        kstar=0.584*(lut%sigv)**(-1.)
    ELSE IF(this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
        kstar=0.05618*lut%sig8z_cold**(-1.013)
    END IF

    END FUNCTION kstar

    FUNCTION As(this,lut,cosm)
    class(THalofit) :: this
    !Halo concentration pre-factor from Bullock et al. (2001) relation
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: As
    REAL(dl) :: B0, Bz, theta

    IF(this%imead==0 .OR. this%imead==4) THEN
        !Set to 4 for the standard Bullock value
        As=4.
    ELSE IF(this%imead==5) THEN
        theta=cosm%logT_AGN-7.8
        B0=3.44-0.496*theta
        Bz=-0.0671-0.0371*theta
        As=B0*10**(lut%z*Bz)
    ELSE IF(this%imead==1 .or. this%imead==2) THEN
        !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
        !As=3.13
        !AM - added for easy modification of feedback parameter
        As=cosm%A_baryon
    ELSE IF(this%imead==3) THEN
        As=5.196
    END IF

    END FUNCTION As

    FUNCTION fdamp(this,lut)
    class(THalofit) :: this
    !Linear power damping function from Mead et al. (2015; arXiv 1505.07833)
    REAL(dl) ::fdamp
    TYPE(HM_tables), INTENT(IN) :: lut

    !Linear theory damping factor
    IF(this%imead==0 .OR. this%imead==4 .OR. this%imead==5) THEN
        !Set to 0 for the standard linear theory two halo term
        fdamp=0.
    ELSE IF(this%imead==1) THEN
        !Mead et al. (2016; arXiv 1602.02154) value
        fdamp=0.0095*lut%sigv100**1.37
    ELSE IF(this%imead==2) THEN
        !Mead et al. (2015) value
        fdamp=0.188*lut%sig8z**4.29
    ELSE IF(this%imead==3) THEN
        fdamp=0.2696*lut%sig8z_cold**0.9403
    END IF

    !Catches extreme values of fdamp
    IF(fdamp<fdamp_min) fdamp=fdamp_min
    IF(fdamp>fdamp_max) fdamp=fdamp_max

    END FUNCTION fdamp

    FUNCTION alpha(this,lut)
    class(THalofit) :: this
    !Two- to one-halo transition smoothing from Mead et al. (2015; arXiv 1505.07833)
    REAL(dl) :: alpha
    TYPE(HM_tables), INTENT(IN) :: lut

    IF(this%imead==0 .OR. this%imead==4 .OR. this%imead==5) THEN
        !Set to 1 for the standard halomodel sum of one- and two-halo terms
        alpha=1.
    ELSE IF(this%imead==1) THEN
        !This uses the top-hat defined neff (HALOFIT uses Gaussian filtered fields instead)
        !Mead et al. (2016; arXiv 1602.02154) value
        alpha=3.24*1.85**lut%neff
    ELSE IF(this%imead==2) THEN
        !Mead et al. (2015) value
        alpha=2.93*1.77**lut%neff
    ELSE IF (this%imead==3) THEN
        alpha=1.875*(1.603)**lut%neff
    END IF

    !Catches values of alpha that are crazy
    IF(alpha>alpha_max) alpha=alpha_max
    IF(alpha<alpha_min) alpha=alpha_min

    END FUNCTION alpha

    FUNCTION r_nl(lut)
    !Calculates R_nl, defined by nu(R_nl)=1., nu=dc/sigma(R)
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL(dl) :: r_nl
    INTEGER, PARAMETER :: iorder=3
    INTEGER, PARAMETER :: ifind=3
    INTEGER, PARAMETER :: imeth=2

    IF(lut%nu(1)>1.) THEN
        !This catches some very strange values that appear for odd cosmological models
        !This is a terrible fudge, but I cannot think of a better solution
        r_nl=lut%rr(1)
    ELSE
        r_nl=exp(find(log(1.d0),log(lut%nu),log(lut%rr),lut%n,iorder,ifind,imeth))
    END IF

    END FUNCTION r_nl

    SUBROUTINE halomod(this,k,p1h,p2h,pfull,plin,lut,cosm)
    class(THalofit) :: this
    !Calcuates 1-halo and 2-halo terms and combines them to form the full halomodel power
    REAL(dl), INTENT(OUT) :: p1h, p2h, pfull
    REAL(dl), INTENT(IN) :: plin, k
    REAL(dl) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN
        p1h=0.
        p2h=0.
    ELSE
        p1h=this%p_1h(k,lut,cosm)
        p2h=this%p_2h(k,plin,lut,cosm)
    END IF

    IF (this%imead==1 .OR. this%imead==2 .OR. this%imead==3) THEN
        a=this%alpha(lut)
        pfull=(p2h**a+p1h**a)**(1./a)
    ELSE
        pfull=p2h+p1h
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
    arr=0

    IF(n==1) THEN
        arr(1)=min
    ELSE IF(n>1) THEN
        DO i=1,n
            arr(i)=min+(max-min)*real(i-1,dl)/(n-1)
        END DO
    END IF

    END SUBROUTINE fill_table

    SUBROUTINE fill_plintab(iz,cosm,CAMB_PK)
    !Fills internal HMcode HM_tables for the linear power spectrum at z=0
    TYPE(MatterPowerData), INTENT(IN) :: CAMB_PK
    INTEGER, INTENT(IN) :: iz
    TYPE(HM_cosmology) :: cosm
    INTEGER :: i
    REAL(dl) :: z, g
    REAL(dl), ALLOCATABLE :: k(:), Pk(:), Pkc(:)
    REAL(dl), PARAMETER :: pi=pi_HM
    REAL(dl), PARAMETER :: kmin=kmin_pk_interpolation
    REAL(dl), PARAMETER :: kmax=kmax_pk_interpolation
    INTEGER :: nk=nk_pk_interpolation, index_cache

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Filling linear power HM_tables'

    !Fill arrays
    IF(ALLOCATED(cosm%log_k_plin)) DEALLOCATE(cosm%log_k_plin)
    IF(ALLOCATED(cosm%log_plin))   DEALLOCATE(cosm%log_plin)
    IF(ALLOCATED(cosm%log_plinc))  DEALLOCATE(cosm%log_plinc)

    IF(rebin_pk) THEN

        !Fill a k-table with an equal-log-spaced k range
        !Note that the minimum should be such that the linear spectrum is accurately a power-law below this wavenumber
        cosm%nk=nk
        CALL fill_table(log(kmin),log(kmax),cosm%log_k_plin,nk)

    ELSE
        if (ifind_pk_interpolation==1) error stop 'ifind_pk_interpolation=1 assumes rebin_pk'
        !Fill k-table with the same k points as in the CAMB calculation
        !If a user has specified lots of points this could make the halo-model
        !calculation chug
        nk=CAMB_PK%num_k
        cosm%nk=nk
        ALLOCATE(cosm%log_k_plin(nk))
        cosm%log_k_plin=CAMB_Pk%log_kh

    END IF

    ALLOCATE(k(nk))
    k=exp(cosm%log_k_plin)
    cosm%kmax = k(nk)

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: k_min:', k(1)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: k_max:', k(nk)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: nk:', nk

    ALLOCATE(Pk(nk),Pkc(nk))

    !Find the redshift
    z=CAMB_Pk%Redshifts(iz)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: z of input:', z
    index_cache = 1
    !Fill power table, both cold- and all-matter
    !$OMP PARALLEL DO DEFAULT(SHARED), FIRSTPRIVATE(index_cache)
    DO i=1,nk
        !Take the power from the current redshift choice
        Pk(i)=MatterPowerData_k(CAMB_PK,k(i),iz, index_cache)*(k(i)**3/(2*pi**2))
        Pkc(i)=Pk(i)*Tcb_Tcbnu_ratio(k(i),z,cosm)**2
    END DO

    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Delta2_min:', Pk(1)
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: Delta2_max:', Pk(nk)

    !Calculate the growth factor at the redshift of interest
    g=grow(z,cosm)
    cosm%grow2_z = g**2
    cosm%this_z = z
    ALLOCATE(cosm%log_plin(nk),cosm%log_plinc(nk))

    !Grow the power to z=0
    cosm%log_plin=log(Pk/(g**2))
    cosm%log_plinc=log(Pkc/(g**2))

    !Check sigma_8 value
    IF(HM_verbose) WRITE(*,*) 'LINEAR POWER: sigma_8:', sigma_integral(8.d0,0.d0,0,cosm)
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

    IF(cosm%f_nu==0._dl) THEN

        Tcb_Tcbnu_ratio=1.

    ELSE

        !Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
        pcb=(5.-sqrt(1+24*(1._dl-cosm%f_nu)))/4

        !Theta for temperature (BigT=T/2.7 K)
        BigT=cosm%Tcmb/2.7

        !The matter-radiation equality redshift
        zeq=(2.5e4)*cosm%om_m*(cosm%h**2)*(BigT**(-4.))

        !The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m \approx 1)
        IF (EdS_Tcold_growth) THEN
            D=(1.+zeq)/(1.+z) ! EdS solution
        ELSE
            D=(1.+zeq)*ungrow(z, cosm) ! General solution
        END IF

        !Wave number relative to the horizon scale at equality (equation 5)
        !Extra factor of h becauase all my k are in units of h/Mpc
        q=k*cosm%h*BigT**2/(cosm%om_m*cosm%h**2)

        !Free streaming scale (equation 14)
        !Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
        !with Nnu of these being massve with the mass split evenly between Nnu species.
        yfs=17.2*cosm%f_nu*(1.+0.488*cosm%f_nu**(-7./6.))*(cosm%Nnu*q/cosm%f_nu)**2

        !These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
        !Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
        !Equations (12) and (13)
        Dcb=(1.+(D/(1.+yfs))**0.7)**(pcb/0.7)
        Dcbnu=((1-cosm%f_nu)**(0.7/pcb)+(D/(1.+yfs))**0.7)**(pcb/0.7)

        Tcb_Tcbnu_ratio=Dcb/Dcbnu

    END IF

    END FUNCTION Tcb_Tcbnu_ratio

    SUBROUTINE assign_HM_cosmology(this,State,cosm)
    class(THalofit) :: this
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
        cosm%om_nu=CP%omnuh2/h2
        cosm%om_v=State%omega_de
        call CP%DarkEnergy%Effective_w_wa(cosm%w, cosm%wa)
        cosm%f_nu=cosm%om_nu/cosm%om_m
        cosm%h=CP%H0/100
        cosm%Tcmb=CP%tcmb
        cosm%Nnu=CP%Num_Nu_massive
        cosm%ns= CP%InitPower%Effective_ns()
    end associate

    ! Baryon feedback parameters
    IF(this%halofit_version==halofit_mead2015 .OR. this%halofit_version==halofit_mead2016)  THEN
        cosm%A_baryon = this%HMcode_A_baryon
        cosm%eta_baryon = this%HMcode_eta_baryon
    ELSE IF(this%halofit_version==halofit_mead2020_feedback) THEN
        cosm%logT_AGN = this%HMcode_logT_AGN
    END IF

    !Write out cosmological parameters if necessary
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_m:', cosm%om_m
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_c:', cosm%om_c
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_b:', cosm%om_b
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_nu:', cosm%om_nu
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: Om_v:', cosm%om_v
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: w_0:', cosm%w
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: w_a:', cosm%wa
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: f_nu:', cosm%f_nu
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: n_s:', cosm%ns
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: h:', cosm%h
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: T_CMB [K]:', cosm%Tcmb
    IF(HM_verbose) WRITE(*,*) 'HM_cosmology: N_nu (massive):', cosm%Nnu
    IF(HM_verbose .AND. (this%halofit_version==halofit_mead2015 .OR. this%halofit_version==halofit_mead2016)) THEN
        WRITE(*,*) 'HM_cosmology: A_baryon:', cosm%A_baryon
        WRITE(*,*) 'HM_cosmology: eta_baryon:', cosm%eta_baryon
    ELSE IF(HM_verbose .AND. this%halofit_version==halofit_mead2020_feedback) THEN
        WRITE(*,*) 'HM_cosmology: log10(T_AGN/K):', cosm%logT_AGN
    END IF
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE assign_HM_cosmology

    SUBROUTINE initialise_HM_cosmology(this,iz,cosm,CAMB_PK)
    class(THalofit) :: this
    !Sets up HM_tables of sigma, growth and linear power for the HM_cosmology
    TYPE(MatterPowerData), INTENT(IN) :: CAMB_PK
    TYPE(HM_cosmology) :: cosm
    INTEGER, INTENT(IN) :: iz

    !Fill linear power table and grows it to z=0
    CALL fill_plintab(iz,cosm,CAMB_PK)

    !Fill sigma(r) table
    CALL fill_sigtab(this,cosm)

    ! Extract BAO wiggle from P(k)
    ! AM: TODO: Maybe move this so that it is not done every z
    CALL init_wiggle(cosm)

    END SUBROUTINE initialise_HM_cosmology

    SUBROUTINE allocate_LUT(lut, n)
    !Allocates memory for the HMcode look-up HM_tables
    TYPE(HM_tables) :: lut
    INTEGER, INTENT(IN) :: n

    if (.not. allocated(lut%zc)) then
        lut%n =n
        ALLOCATE(lut%zc(n),lut%m(n),lut%c(n),lut%rv(n))
        ALLOCATE(lut%nu(n),lut%rr(n),lut%sigf(n),lut%sig(n))
    end if
    lut%zc=0
    lut%m=0
    lut%c=0
    lut%rv=0
    lut%nu=0
    lut%rr=0
    lut%sigf=0
    lut%sig=0

    END SUBROUTINE allocate_LUT

    SUBROUTINE halomod_init(this,z,lut,cosm)
    class(THalofit) :: this
    !Halo-model initialisation routine
    !Computes look-up HM_tables necessary for the halo model calculations
    REAL(dl), INTENT(IN) :: z
    INTEGER :: i,nm
    REAL(dl) :: Dv, dc, m, nu, r, sig, mmin, mmax
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut
    REAL(dl), PARAMETER :: f_Bullock=0.01_dl**(1/3._dl)

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: Filling look-up HM_tables'
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: HM_tables being filled at redshift:', z
    lut%z=z

    ! Mass range and number of points
    mmin=1e0
    mmax=1e18
    nm=256

    !Find value of sigma_v, sig8, etc.
    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    lut%sigv=sigmaV(0.d0,z,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
    !$OMP SECTION
    if (global_error_flag==0) lut%sigv100=sigmaV(100.d0,z,0,cosm)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
    !$OMP SECTION
    if (global_error_flag==0) then
        lut%sig8z=sigma_integral(8.d0,z,0,cosm)
        lut%sig8z_cold=sigma_integral(8.d0,z,1,cosm)
    end if
    IF(HM_verbose) THEN
        WRITE(*,*) 'HALOMOD: sig8(z):', lut%sig8z
        WRITE(*,*) 'HALOMOD: cold sig8(z):', lut%sig8z_cold
    END IF
    !$OMP END PARALLEL SECTIONS
    if (global_error_flag/=0) return

    CALL allocate_lut(lut, nm)

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M_min [log10(Msun/h)]:', log10(mmin)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: M_max [log10(Msun/h)]:', log10(mmax)

    dc=this%delta_c(z,lut,cosm)
    lut%dc=dc

    !$OMP PARALLEL DO default(shared), private(m,r,sig,nu)
    DO i=1,lut%n

        m=exp(log(mmin)+log(mmax/mmin)*real(i-1,dl)/(lut%n-1))
        r=radius_m(m,cosm)
        sig=sigma_lut(r,z,cosm)
        nu=dc/sig

        lut%m(i)=m
        lut%rr(i)=r
        lut%sig(i)=sig
        lut%nu(i)=nu
        lut%sigf(i)=sigma_lut(r*f_Bullock,z,cosm)

    END DO
    !$OMP END PARALLEL DO
    if (global_error_flag/=0) return

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: m, r, nu, sig, sigf HM_tables filled'

    !Fill virial radius table using real radius table
    Dv=this%Delta_v(z,cosm)
    lut%rv=lut%rr/(Dv**(1/3._dl))

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: rv HM_tables filled'
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sig min:', lut%sig(lut%n)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: sig max:', lut%sig(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1/lut%rnl

    IF(HM_verbose) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    IF(HM_verbose) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    !Calcuate the effective spectral index at the collapse scale
    lut%neff=neff(this,lut,cosm)

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
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'Dv:', this%Delta_v(z,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'dc:', this%delta_c(z,lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'eta:', this%eta(lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'k*:', this%kstar(lut)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'A:', this%As(lut,cosm)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'fdamp:', this%fdamp(lut)
    IF(HM_verbose) WRITE(*,fmt='(A10,F10.5)') 'alpha:', this%alpha(lut)
    IF(HM_verbose) WRITE(*,*) '=================================='
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE write_parameters

    PURE FUNCTION radius_m(m,cosm)
    !Calculates the co-moving radius that encloses a mass 'm' in the homogeneous Universe
    REAL(dl) :: radius_m
    REAL(dl), INTENT(IN) :: m
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: pi=pi_HM

    radius_m=(3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

    END FUNCTION radius_m

    FUNCTION neff(this,lut,cosm)
    class(THalofit) :: this
    !Finds the effective spectral index at the collapse scale r_nl, where nu(r_nl)=1.
    REAL(dl) :: neff
    REAL(dl) :: ns
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    TYPE(HM_tables), INTENT(IN) :: lut
    REAL(dl) :: sig
    INTEGER :: itype
    REAL(dl), PARAMETER :: tmin=0.
    REAL(dl), PARAMETER :: tmax=1.
    REAL(dl), PARAMETER :: acc=acc_neff_integration
    INTEGER, PARAMETER :: iorder=iorder_neff_integration

    ! Choose type of sigma(R) to tabulate depending on HMcode version
    IF (this%imead==1 .OR. this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
        itype=1 ! 1 - Cold matter
    ELSE
        itype=0 ! 0 - All matter
    END IF

    ! Choosings sig = delta_c should be equivalent to actually calculating it again, however
    ! Do the actual calculation to be consistent with HMx. Problems with weird cosmologies with
    ! low spectral indices such that no collapse has occured. R_nl very small
    !sig=lut%dc ! Take great care here. This should be the same as below, but won't be for strange models
    sig=sigma_integral(lut%rnl,lut%z,itype,cosm)
    neff=-3.-2.*integrate(tmin,tmax,neff_integrand,lut%rnl,lut%z,itype,cosm,acc,iorder)/sig**2

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
    TYPE(HM_cosmology) :: cosm, cosm_lcdm
    TYPE(HM_tables) :: lut
    REAL(dl) :: A, zf, pow
    REAL(dl) :: ginf_lcdm, ginf_wcdm, g_lcdm, g_wcdm
    INTEGER :: i
    REAL(dl), PARAMETER :: zinf=zc_Dolag

    !Amplitude of relation (4. in Bullock et al. 2001)
    A=this%As(lut,cosm)

    !Fill the collapse time look-up table
    CALL this%zcoll_bull(z,cosm,lut)

    !Fill the concentration look-up table
    DO i=1,lut%n
        zf=lut%zc(i)
        lut%c(i)=A*(1.+zf)/(1.+z)
    END DO

    IF(z<zinf) THEN

        !Dolag (2004) prescription for adding DE dependence

        !This is approximately z=infinity
        ginf_wcdm=grow(zinf,cosm)

        !Make a LCDM HM_cosmology
        !Need to make sure model is flat with the same Omega_m and w=-1
        !This is *only* used for a calculation of the growth function
        cosm_lcdm=cosm
        DEALLOCATE(cosm_lcdm%growth)
        DEALLOCATE(cosm_lcdm%a_growth)
        cosm_lcdm%w=-1.
        cosm_lcdm%wa=0.
        cosm_lcdm%om_v=1.-cosm%om_m !Enforce flatness

        !Needs to use grow_int explicitly here for LCDM model to avoid growth HM_tables
        ginf_lcdm=growint(zinf,cosm_lcdm)

        !This is the Dolag et al. (2004) correction for halo concentrations
        IF(this%imead==0 .OR. this%imead==2 .OR. this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
            ! Mead et al. (2015) used the Dolag (2004) correction
            pow=1.
        ELSE IF(this%imead==1) THEN
            ! Mead et al. (2016) changed the power to 1.5 to better accomodate more extreme dark-energy models
            pow=1.5
        END IF
        lut%c=lut%c*(ginf_wcdm/ginf_lcdm)**pow

        ! This is needed for the correction to make sense at high z
        IF(this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
            g_lcdm=growint(z,cosm_lcdm)
            g_wcdm=grow(z,cosm)
            lut%c=lut%c*(g_lcdm/g_wcdm)**pow
        END IF

    END IF

    END SUBROUTINE conc_bull

    FUNCTION growint(z,cosm)
    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    REAL(dl) :: growint
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL(dl) :: x, dx
    REAL(dl) :: f1, f2, fx
    REAL(dl) :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=20
    REAL(dl), PARAMETER :: acc=acc_growth_integration
    INTEGER, PARAMETER :: iorder=iorder_growth_integration

    a=1./(1.+z)

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
                    ERROR STOP 'GROWINT: Error, iorder specified incorrectly'
                END IF

            END IF

            IF(j>=jmin .and. sum_old /= 0) THEN
                if (ABS(-1.d0+sum_new/sum_old)<acc) THEN
                    !jmin avoids spurious early convergence
                    growint=exp(sum_new)
                    EXIT
                end if
            end if
            IF(j==jmax) THEN
                growint=exp(sum_new)
                call GlobalError('HMCode GROWINT, Integration timed out', error_nonlinear)
                return
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
    REAL(dl) :: growint_integrand
    REAL(dl), INTENT(IN) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: Om_m, z, gam

    IF(cosm%w<-1.) THEN
        gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
        gam=0.55+0.05*(1.+cosm%w)
    ELSE
        gam=0.55
    END IF

    !Note the minus sign here
    !AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    z = -1.+1./a
    IF (cold_growth) THEN
        Om_m = Omega_cold_hm(z,cosm)
    ELSE
        Om_m = Omega_m_hm(z,cosm)
    END IF
    growint_integrand=-(Om_m**gam)/a

    END FUNCTION growint_integrand

    SUBROUTINE zcoll_bull(this,z,cosm,lut)
    class(THalofit) :: this
    !Calcuates the halo collapse redshift according to the Bullock (2001) prescription
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology) :: cosm
    TYPE(HM_tables) :: lut
    REAL(dl) :: dc
    REAL(dl) :: af, zf, RHS, growz
    INTEGER :: i
    INTEGER, PARAMETER :: iorder=iorder_growth_inversion
    INTEGER, PARAMETER :: ifind=ifind_growth_inversion
    INTEGER, PARAMETER :: imeth=imeth_growth_inversion

    !This fills up the halo formation redshift table as per Bullock relations

    !Needs to interpolate g(z) which should be pretty linear for a<0.05
    !in 'g(a) vs a' space for all standard cosmologies

    dc=this%delta_c(z,lut,cosm)

    !Find the growth function at the current redshift
    growz=grow(z,cosm)

    !Do numerical inversion
    DO i=1,lut%n

        RHS=dc*growz/lut%sigf(i)

        IF(RHS>growz) THEN
            !This is the case of 'halo forms in the future'
            !in this case set formation redshift to current redshift
            zf=z
        ELSE
            af=find(RHS,cosm%growth,cosm%a_growth,cosm%ng,iorder,ifind,imeth)
            zf=-1.+1./af
        END IF

        lut%zc(i)=zf

    END DO

    END SUBROUTINE zcoll_bull

    FUNCTION mass_r(r,cosm)
    !Calcuates the average mass enclosed at co-moving radius r
    REAL(dl) :: mass_r, r
    TYPE(HM_cosmology) :: cosm
    REAL(dl), PARAMETER :: pi=pi_HM

    !Relation between mean cosmological mass and radius
    mass_r=(4*pi/3)*cosmic_density(cosm)*(r**3)

    END FUNCTION mass_r

    PURE FUNCTION cosmic_density(cosm)
    !The z=0 cosmological matter density
    REAL(dl) :: cosmic_density
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !In M_sun per Mpc^3 with h factors included. The constant does this.
    cosmic_density=(2.775d11)*cosm%om_m

    END FUNCTION cosmic_density

    FUNCTION baryonify_wk(wk, m, lut, cosm)

    ! Change halo window to account for feedback mass loss and star formation
    IMPLICIT NONE
    REAL(dl) :: baryonify_wk
    REAL(dl), INTENT(IN) :: wk
    REAL(dl), INTENT(IN) :: m
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: wkn, fs

    wkn = wk
    wkn = wkn*halo_mass_fraction(m, lut, cosm) ! Gas expulsion
    fs = f_star(lut, cosm)
    wkn = wkn+fs ! Star formation

    baryonify_wk = wkn

    END FUNCTION baryonify_wk

    FUNCTION halo_mass_fraction(m, lut, cosm)

    ! Simple baryon model where high-mass haloes have a mass fraction of 1 and low-mass haloes have Omega_c/Omega_m
    ! This also accounts for massive neutrinos since it is written in terms of Om_c and Om_b (rather than Om_m)
    ! TODO: Could just precompute this once for each M in the halomod init function
    IMPLICIT NONE
    REAL(dl) :: halo_mass_fraction
    REAL(dl), INTENT(IN) :: m
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: r, fc, fb, fs, mb, beta

    mb = m_baryon(lut, cosm) ! Feedback halo mass
    beta = 2.                ! Power-law index
    r = (m/mb)**beta         ! If m>>m0 then r becomes large, if m<<m0 then r=0
    fb = cosm%Om_b/cosm%Om_m ! Halo baryon fraction
    fc = cosm%Om_c/cosm%Om_m ! Halo CDM fraction
    fs = f_star(lut, cosm)   ! Halo star fraction
    halo_mass_fraction = fc+(fb-fs)*r/(1.+r) ! Remaining halo mass fraction

    END FUNCTION halo_mass_fraction

    FUNCTION m_baryon(lut, cosm)

    REAL(dl) :: m_baryon
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: mb0, mbz, theta

    theta=cosm%logT_AGN-7.8
    mb0=13.87+1.81*theta
    mb0=10**mb0
    mbz=-0.108+0.195*theta
    m_baryon=mb0*10**(mbz*lut%z)

    END FUNCTION m_baryon

    FUNCTION f_star(lut, cosm)

    REAL(dl) :: f_star
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: f0, fz, theta

    theta=cosm%logT_AGN-7.8
    f0=0.0201-0.003*theta
    fz=0.409+0.0224*theta
    f_star=f0*10**(fz*lut%z)

    IF(f_star > cosm%Om_b/cosm%Om_m) f_star = cosm%Om_b/cosm%Om_m

    END FUNCTION f_star

    FUNCTION find_pk(k,itype,cosm)
    !Look-up and interpolation for P(k,z=0)
    REAL(dl) :: find_pk
    REAL(dl) :: ns
    REAL(dl), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: itype
    INTEGER :: n
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, PARAMETER :: iorder=iorder_pk_interpolation
    INTEGER, PARAMETER :: ifind=ifind_pk_interpolation
    INTEGER, PARAMETER :: imeth=imeth_pk_interpolation

    !Set number of k points as well as min and max k values
    !Note that the min k value should be set to the same as the CAMB min k value
    n=SIZE(cosm%log_k_plin)

    IF(plin_extrap .AND. k>cosm%kmax) THEN
        !Do some extrapolation here based on knowledge of things at high k
        ns=cosm%ns !Spectral index used in the high-k extrapolation
        IF(itype==0) THEN
            find_pk=exp(cosm%log_plin(n))*((log(k)/cosm%log_k_plin(n))**2)*((k/cosm%kmax)**(ns-1))
        ELSE IF(itype==1) THEN
            find_pk=exp(cosm%log_plinc(n))*((log(k)/cosm%log_k_plin(n))**2)*((k/cosm%kmax)**(ns-1))
        END IF
    ELSE
        !Otherwise use the standard find algorithm
        IF(itype==0) THEN
            find_pk=exp(find(log(k),cosm%log_k_plin,cosm%log_plin,cosm%nk,iorder,ifind,imeth))
        ELSE IF(itype==1) THEN
            find_pk=exp(find(log(k),cosm%log_k_plin,cosm%log_plinc,cosm%nk,iorder,ifind,imeth))
        END IF
    END IF

    END FUNCTION find_pk

    FUNCTION p_lin(k,z,itype,cosm)
    !Looks up the value for the linear power spectrum
    REAL(dl) :: p_lin
    REAL(dl), INTENT (IN) :: k, z
    REAL(dl) growth2
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised to z=0
    if (z==0._dl) then
        growth2 = 1
    else if (z==cosm%this_z) then
        growth2 = cosm%grow2_z
    else
        growth2 = grow(z,cosm)**2 !never actually needed
    end if
    p_lin=growth2*find_pk(k,itype,cosm)

    END FUNCTION p_lin

    FUNCTION p_2h(this,k,plin,lut,cosm)
    class(THalofit) :: this
    !Calculates the 2-halo term
    REAL(dl) :: p_2h
    REAL(dl), INTENT(IN) :: k, plin
    REAL(dl) :: sigv, frac, kdamp, ndamp, x
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    !Damping function
    frac=this%fdamp(lut)

    IF(this%imead==0 .OR. this%imead==4 .OR. this%imead==5 .OR. frac<1.e-3) THEN
        p_2h=plin
    ELSE IF(this%imead==1 .OR. this%imead==2) THEN
        sigv=lut%sigv
        p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2.)
    ELSE IF(this%imead==3) THEN
        kdamp=0.05699*lut%sig8z_cold**(-1.089)
        ndamp=2.85
        x=(k/kdamp)**ndamp
        p_2h=p_dewiggle(k, lut%z, plin, lut%sigv, cosm)*(1.-frac*x/(1.+x))
    END IF

    !For some strange cosmologies frac>1. so this must be added to prevent p_2h<0.
    IF(p_2h<0.) p_2h=0.

    END FUNCTION p_2h

    FUNCTION p_1h(this,k,lut,cosm)
    class(THalofit) :: this
    !Calculates the 1-halo term
    REAL(dl) :: p_1h
    REAL(dl), INTENT(IN) :: k
    TYPE(HM_tables), INTENT(IN) :: lut
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: g, fac, et, ks, wk, x, m
    REAL(dl) :: integrand(lut%n)
    REAL(dl) :: sum
    INTEGER :: i
    REAL(dl), PARAMETER :: pi=pi_HM
    INTEGER, PARAMETER :: iorder=iorder_integration_1h

    !Only call eta once
    et=this%eta(lut,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
        m=lut%m(i)
        g=gnu(lut%nu(i))
        wk=win(k*lut%nu(i)**et,lut%rv(i),lut%c(i))
        IF(this%imead==5) wk=baryonify_wk(wk, m, lut, cosm)
        integrand(i)=g*(wk**2)*lut%m(i)
    END DO

    IF(this%imead==3 .OR. this%imead==4) integrand=integrand*(1.-cosm%f_nu)**2

    !Carries out the integration
    sum=inttab(lut%nu,integrand,1,lut%n,iorder)/cosmic_density(cosm)

    !Numerical factors to convert from P(k) to Delta^2(k)
    p_1h=sum*k**3/(2.*pi**2)

    IF(this%imead==1 .OR. this%imead==2) THEN
        !Damping of the 1-halo term at very large scales
        ks=this%kstar(lut)
        IF(ks==0.) THEN
            fac=0.
        ELSE IF((k/ks)**2.>ks_limit) THEN
            fac=0. !Prevents problems if k/ks is very large
        ELSE
            fac=exp(-((k/ks)**2.))
        END IF
        !Damping of the one-halo term at very large scales
        p_1h=p_1h*(1.-fac)
    ELSE IF(this%imead==3 .OR. this%imead==4 .OR. this%imead==5) THEN
        ks=this%kstar(lut)
        x=(k/ks)**4
        p_1h=p_1h*x/(1.+x)
    END IF

    END FUNCTION p_1h

    REAL FUNCTION p_dewiggle(k, z, p_linear, sigv, cosm)
    ! Call the dewiggled power spectrum, which is linear but with damped wiggles
    REAL(dl), INTENT(IN) :: k, z, p_linear, sigv
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: p_wiggle, f, logk
    INTEGER, PARAMETER :: iorder = iorder_wiggle
    INTEGER, PARAMETER :: ifind = ifind_wiggle
    INTEGER, PARAMETER :: imeth = imeth_wiggle

    logk = log(k)
    IF (logk < cosm%log_k_wiggle(1) .OR. logk > cosm%log_k_wiggle(nk_wiggle)) THEN
        p_wiggle = 0
    ELSE
        p_wiggle = find(logk, cosm%log_k_wiggle, cosm%pk_wiggle, nk_wiggle, iorder, ifind, imeth)
    END IF
    f = exp(-(k*sigv)**2)
    p_dewiggle = p_linear+(f-1)*p_wiggle*grow(z, cosm)**2

    END FUNCTION p_dewiggle

    SUBROUTINE init_wiggle(cosm)
    ! Isolate the power spectrum wiggle
    TYPE(HM_cosmology), INTENT(INOUT) :: cosm
    REAL(dl), ALLOCATABLE :: k(:), Pk(:)
    REAL(dl), ALLOCATABLE :: Pk_smooth(:), Pk_wiggle(:)
    INTEGER :: i
    REAL(dl), PARAMETER :: kmin = kmin_wiggle
    REAL(dl), PARAMETER :: kmax = kmax_wiggle
    INTEGER, PARAMETER  :: nk = nk_wiggle
    REAL(dl), PARAMETER :: z = 0. ! Only need to run this routine for z=0
    INTEGER, PARAMETER :: itype = 0 ! Matter spectrum

    ! Words
    IF (HM_verbose) WRITE(*, *) 'INIT_WIGGLE: Starting'

    ! Allocate arrays
    ALLOCATE (Pk(nk), Pk_wiggle(nk), Pk_smooth(nk))

    ! Allocate array for k
    CALL fill_table(log(kmin), log(kmax), k, nk)
    k=exp(k)

    ! Get the linear power spectrum in an array
    DO i = 1, nk
        Pk(i) = p_lin(k(i), z, itype, cosm)
    END DO

    IF (HM_verbose) THEN
        WRITE(*, *) 'INIT_WIGGLE: kmin [h/Mpc]:', k(1)
        WRITE(*, *) 'INIT_WIGGLE: kmax [h/Mpc]:', k(nk)
        WRITE(*, *) 'INIT_WIGGLE: nk:', nk
        WRITE(*, *) 'INIT_WIGGLE: Splitting into wiggle and broad-band'
    END IF

    CALL calculate_psmooth(k, z, Pk, Pk_smooth, cosm)

    IF (HM_verbose) WRITE(*, *) 'INIT_WIGGLE: Isolating wiggle'

    ! Isolate the wiggle
    Pk_wiggle = Pk-Pk_smooth

    IF (HM_verbose) WRITE(*, *) 'INIT_WIGGLE: Initialising interpolator'

    ! Fill look-up tables
    IF(ALLOCATED(cosm%log_k_wiggle)) DEALLOCATE(cosm%log_k_wiggle)
    IF(ALLOCATED(cosm%pk_wiggle)) DEALLOCATE(cosm%pk_wiggle)
    ALLOCATE(cosm%log_k_wiggle(nk), cosm%pk_wiggle(nk))
    cosm%log_k_wiggle = log(k)
    cosm%pk_wiggle = pk_wiggle

    IF (HM_verbose) THEN
        WRITE (*, *) 'INIT_WIGGLE: Done'
        WRITE (*, *)
    END IF

    END SUBROUTINE init_wiggle

    SUBROUTINE calculate_psmooth(k, z, Pk, Pk_smt, cosm)
    ! Calculate the normalised smoothed power spectrum at a range of k
    REAL(dl), INTENT(IN) :: k(:)
    REAL(dl), INTENT(IN) :: z
    REAL(dl), INTENT(IN) :: Pk(:)
    REAL(dl), ALLOCATABLE, INTENT(OUT) :: Pk_smt(:)
    REAL(dl), ALLOCATABLE :: Pk_nw(:)
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: sig = wiggle_sigma

    ! Reduce dynamic range
    CALL calculate_nowiggle(k, z, Pk, Pk_nw, cosm)
    Pk_smt = Pk/Pk_nw

    ! Smooth linear power
    CALL smooth_array_Gaussian(log(k), Pk_smt, sig)

    ! Return dynamic range
    Pk_smt = Pk_smt*Pk_nw

    END SUBROUTINE calculate_psmooth

    SUBROUTINE calculate_nowiggle(k, z, Pk, Pk_nw, cosm)
    ! Calculate the normalised no wiggle power spectrum at a range of k and a
    ! Comes from the Eisenstein & Hu approximation
    REAL(dl), INTENT(IN) :: k(:)
    REAL(dl), INTENT(IN) :: z
    REAL(dl), INTENT(IN) :: Pk(:)
    REAL(dl), ALLOCATABLE, INTENT(OUT) :: Pk_nw(:)
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER :: ik, nk
    REAL(dl) :: Pk_norm, Pk_nw_norm
    REAL(dl), PARAMETER :: knorm = knorm_nowiggle
    INTEGER, PARAMETER :: type = 0 ! Matter here
    INTEGER, PARAMETER :: iorder = 3
    INTEGER, PARAMETER :: ifind = 3
    INTEGER, PARAMETER :: imeth = 2

    ! Allocate arrays
    nk = size(k)
    IF (nk /= size(Pk)) STOP 'CALCULATE_NOWIGGLE: Error, Pk should be same size as k and a'
    ALLOCATE(Pk_nw(nk))

    ! Get the no-wiggle power spectrum
    DO ik = 1, nk
        Pk_nw(ik) = Pk_nowiggle(k(ik), cosm)
    END DO

    ! Calculate the no-wiggle power spectrum and force spectra to agree at the normalisation wavenumber
    Pk_norm = p_lin(knorm, z, type, cosm)
    Pk_nw_norm = find(knorm, k, Pk_nw, nk, iorder, ifind, imeth)
    Pk_nw = Pk_nw*Pk_norm/Pk_nw_norm

    END SUBROUTINE calculate_nowiggle

    REAL FUNCTION Pk_nowiggle(k, cosm)
    ! Calculates the un-normalised no-wiggle power spectrum
    ! Comes from the Eisenstein & Hu approximation
    REAL(dl), INTENT(IN) :: k
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    Pk_nowiggle = (k**(cosm%ns+3.))*Tk_nw(k, cosm)**2

    END FUNCTION Pk_nowiggle

    REAL FUNCTION Tk_nw(k, cosm)
    ! No-wiggle transfer function from Eisenstein & Hu: astro-ph:9709112
    REAL(dl), INTENT(IN) :: k ! Wavenumber [h/Mpc]
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: q, L, C, Gamma, wm, wb, s, h, alpha, rb
    REAL(dl), PARAMETER :: e = exp(1.)

    ! Useful parameters to make equations shorter
    wm = cosm%Om_m*cosm%h**2 ! Real matter density
    wb = cosm%Om_b*cosm%h**2 ! Real baryon density
    rb = cosm%Om_b/cosm%Om_m ! Baryon ratio
    h = cosm%h               ! Hubble factor

    ! These only needs to be calculated once
    s = 44.5*log(9.83/wm)/sqrt(1.+10.*wb**0.75)              ! Equation (26)
    alpha = 1.-0.328*log(431.*wm)*rb+0.38*log(22.3*wm)*rb**2 ! Equation (31)

    ! Functions of k
    Gamma = cosm%Om_m*cosm%h*(alpha+(1.-alpha)/(1.+(0.43*k*s*h)**4)) ! Equation (30)
    q = k*(cosm%Tcmb/2.7)**2/Gamma ! Equation (28)
    L = log(2.*e+1.8*q)             ! Equation (29)
    C = 14.2+731./(1.+62.5*q)       ! Equation (29)
    Tk_nw = L/(L+C*q**2)            ! Equation (29)

    END FUNCTION Tk_nw

    SUBROUTINE smooth_array_Gaussian(x, f, sigma)
    ! Smooth an array f(x) using a Gaussian kernel
    REAL(dl), INTENT(IN) :: x(:)    ! x coordinates
    REAL(dl), INTENT(INOUT) :: f(:) ! Array to smooth
    REAL(dl), INTENT(IN) :: sigma   ! Width of smoothing Gaussian
    INTEGER :: i, j, n
    REAL(dl), ALLOCATABLE :: ff(:)
    REAL(dl) :: weight, total
    REAL(dl), PARAMETER :: nsig = 3. ! Do not smooth if point lies within this number of sigma from edge

    IF (sigma  .NE. 0.) THEN

        n = size(x)
        IF (n /= size(f)) STOP 'GAUSSIAN_SMOOTH_ARRAY: Error, x and y should be the same size'

        ! Save the original input array
        ff = f

        ! Delete the original array
        f = 0.

        ! Apply Gaussian smoothing
        DO i = 1, n
            total = 0.
            IF (abs(x(i)-x(1)) < nsig*sigma .OR. abs(x(i)-x(n)) < nsig*sigma) THEN
                f(i) = ff(i)
            ELSE
                DO j = 1, n
                    weight = exp(-(x(i)-x(j))**2/(2.*sigma**2))
                    f(i) = f(i)+ff(j)*weight
                    total = total+weight
                END DO
                f(i) = f(i)/total
            END IF
        END DO

        DEALLOCATE(ff)

    END IF

    END SUBROUTINE smooth_array_Gaussian

    SUBROUTINE fill_sigtab(this,cosm)
    class(THalofit) :: this
    !This fills up HM_tables of r vs. sigma(r) across a range in r
    !It is used only in look-up for further calculations of sigmac(r) and not otherwise
    !and prevents a large number of calls to the sigint functions
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is approximately power-law below and above these values of R
    !This wouldn't be appropriate for models with a small-scale linear spectrum cut-off (e.g., WDM)
    INTEGER :: i
    TYPE(HM_cosmology) :: cosm
    REAL(dl), ALLOCATABLE :: r(:), sig(:)
    INTEGER :: itype
    REAL(dl), PARAMETER :: rmin=rmin_sigma_interpolation
    REAL(dl), PARAMETER :: rmax=rmax_sigma_interpolation
    INTEGER, PARAMETER :: nsig=n_sigma_interpolation

    ! Choose type of sigma(R) to tabulate depending on HMcode version
    IF(this%halofit_version == halofit_mead2015) THEN
        itype=0 ! 0 - All matter
    ELSE
        itype=1 ! 1 - Cold matter
    END IF

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

        sig(i)=sigma_integral(r(i),0.d0,itype,cosm)

    END DO
    !$OMP END PARALLEL DO

    IF(HM_verbose) WRITE(*,*) 'SIGTAB: sigma_min:', sig(nsig)
    IF(HM_verbose) WRITE(*,*) 'SIGTAB: sigma_max:', sig(1)

    cosm%nsig=nsig
    ALLOCATE(cosm%log_r_sigma(nsig),cosm%log_sigma(nsig))
    cosm%log_r_sigma=log(r)
    cosm%log_sigma=log(sig)

    IF(HM_verbose) WRITE(*,*) 'SIGTAB: Done'
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE fill_sigtab

    FUNCTION sigma_lut(r,z,cosm)
    !Finds sigma_cold(R) from look-up table
    REAL(dl) :: sigma_lut
    REAL(dl), INTENT(IN) :: r, z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, PARAMETER :: iorder=iorder_sigma_interpolation
    INTEGER, PARAMETER :: ifind=ifind_sigma_interpolation
    INTEGER, PARAMETER :: imeth=imeth_sigma_interpolation

    sigma_lut=grow(z,cosm)*exp(find(log(r),cosm%log_r_sigma,cosm%log_sigma,cosm%nsig,iorder,ifind,imeth))

    END FUNCTION sigma_lut

    FUNCTION wk_tophat(x)

    !The normlaised Fourier Transform of a top-hat
    REAL(dl) :: wk_tophat
    REAL(dl), INTENT(IN) :: x
    REAL(dl), PARAMETER :: dx=1e-3 ! Taylor expansion for |x|<dx

    !Taylor expansion used for low x to avoid cancellation problems
    IF(abs(x)<dx) THEN
        wk_tophat=1-(x**2)/10
    ELSE
        wk_tophat=3*(sin(x)-x*cos(x))/(x**3)
    END IF

    END FUNCTION wk_tophat

    FUNCTION wk_tophat_deriv(x)
    ! The derivative of a normlaised Fourier Transform of a spherical top-hat
    REAL(dl) :: wk_tophat_deriv
    REAL(dl), INTENT(IN) :: x
    REAL(dl), PARAMETER :: dx=1e-3 ! Taylor expansion for |x|<dx

    ! Taylor expansion used for low x to avoid cancelation problems
    IF (abs(x)<dx) THEN
        wk_tophat_deriv=-x/5+x**3/70
    ELSE
        wk_tophat_deriv=(3/x**4)*((x**2-3)*sin(x)+3*x*cos(x))
    END IF

    END FUNCTION wk_tophat_deriv

    FUNCTION inttab(x,y,n1,n2,iorder)
    !Integrates tables y(x)dx
    REAL(dl) :: inttab
    INTEGER, INTENT(IN) :: n1, n2
    REAL(dl), INTENT(IN) :: x(:), y(:)
    INTEGER, INTENT(IN) :: iorder
    REAL(dl) :: a, b, c, d, h
    REAL(dl) :: q1, q2, q3, qi, qf
    REAL(dl) :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    real(dl) :: sum
    INTEGER :: i, i1, i2, i3, i4, n

    n=size(x)

    sum=0.d0

    IF(n1==n2) THEN

        inttab=0.

    ELSE IF(iorder==1) THEN

        !Sums over all Trapezia (a+b)*h/2
        DO i=n1,n2-1
            a=y(i+1)
            b=y(i)
            h=x(i+1)-x(i)
            sum=sum+(a+b)*h/2.d0
        END DO

    ELSE IF(iorder==2) THEN

        DO i=n1,n2-2

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

        DO i=n1,n2-1

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

    FUNCTION sigma_integral(r,z,itype,cosm)
    !Gets sigma(R)
    REAL(dl) :: sigma_integral
    REAL(dl), INTENT(IN) :: r, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: tmin=0.
    REAL(dl), PARAMETER :: tmax=1.
    REAL(dl), PARAMETER :: acc=acc_sigma_integration
    INTEGER, PARAMETER :: iorder=iorder_sigma_integration

    sigma_integral=sqrt(integrate(tmin,tmax,sigma_integrand,r,z,itype,cosm,acc,iorder))

    END FUNCTION sigma_integral

    FUNCTION sigma_integrand(t,R,z,itype,cosm)
    !The integrand for the sigma(R) integrals
    REAL(dl) :: sigma_integrand
    REAL(dl), INTENT(IN) :: t, R, z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat
    REAL(dl), PARAMETER :: alpha=alpha_sigma_integration

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
    REAL(dl) :: integrate
    REAL(dl), INTENT(IN) :: a
    REAL(dl), INTENT(IN) :: b
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
    INTEGER, PARAMETER :: jmax=20

    INTERFACE
    FUNCTION f(x, y, z, itype, cosm)
    USE precision
    IMPORT :: HM_cosmology
    REAL(dl) :: f
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
            dx=(b-a)/REAL(n-1,dl)

            IF(j==1) THEN

                !The first go is just the trapezium of the end points
                f1=f(a,y,z,itype,cosm)
                f2=f(b,y,z,itype,cosm)
                sum_2n=0.5d0*(f1+f2)*dx
                sum_new=sum_2n

            ELSE

                !Loop over only new even points to add these to the integral
                DO i=2,n,2
                    x=a+(b-a)*real(i-1,dl)/real(n-1,dl)
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
                    ERROR STOP 'INTEGRATE: Error, iorder specified incorrectly'
                END IF

            END IF

            IF(j>=jmin .and. sum_old /= 0) THEN
                if (ABS(-1.d0+sum_new/sum_old)<acc) THEN
                    !jmin avoids spurious early convergence
                    integrate=sum_new
                    EXIT
                end if
            end if
            IF(j==jmax) THEN
                integrate=sum_new
                call GlobalError('HMCode INTEGRATE, Integration timed out', error_nonlinear)
                return
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
    IF(win>1._dl) win=1._dl
    IF(win<0._dl) win=0._dl

    END FUNCTION win

    FUNCTION winnfw(k,rv,c)
    !The analytic Fourier Transform of the NFW profile; note  W(k->0)=1
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

    mass=log(1+c)-c/(1+c)

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
    REAL(dl), PARAMETER :: p=0.3
    REAL(dl), PARAMETER :: a=0.707
    REAL(dl), PARAMETER :: bigA=0.21616

    !Note I use nu=dc/sigma whereas ST (1999) use nu=(dc/sigma)^2
    !This accounts for the different pre-factor and slighly changed nu dependence
    !f(nu^2)d(nu^2)=2*nu*f(nu)dnu

    !Full mass function. Note this is normalised such that integral f(nu)dnu = 1
    gst=bigA*(1+((a*nu*nu)**(-p)))*exp(-a*nu*nu/2)

    END FUNCTION gst

    FUNCTION Hubble2(z,cosm)
    !This calculates the dimensionless squared hubble parameter at redshift z (H/H_0)^2!
    !Ignores contributions from radiation (not accurate at high z, but consistent with simulations)!
    REAL(dl) :: Hubble2
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    Hubble2=cosm%om_m*(1+z)**3+cosm%om_v*X_de(z,cosm)+(1-cosm%om_m-cosm%om_v)*(1+z)**2

    END FUNCTION Hubble2

    FUNCTION AH(z,cosm)
    !The Hubble acceleration function \ddot{a}/a
    REAL(dl) :: AH
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    AH=cosm%om_m*(1+z)**3+cosm%om_v*(1.+3.*w_de_hm(z,cosm))*X_de(z,cosm)
    AH=-AH/2.

    END FUNCTION AH

    FUNCTION X_de(z,cosm)
    !The time evolution for dark energy: rho_de=rho_de,0 * X(a)
    !X(a)=1 for LCDM but changes for other models
    REAL(dl) :: X_de
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: a

    a=1./(1.+z)
    X_de=(a**(-3*(1+cosm%w+cosm%wa)))*exp(-3*cosm%wa*(1-a))

    END FUNCTION X_de

    FUNCTION w_de_hm(z,cosm)

    !The dark energy w(a) function
    REAL(dl) :: w_de_hm
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: a

    a=1./(1.+z)
    w_de_hm=cosm%w+(1-a)*cosm%wa

    END FUNCTION w_de_hm

    FUNCTION Omega_m_hm(z,cosm)
    !This calculates omega_m variations with z!
    REAL(dl) :: Omega_m_hm
    REAL(dl), INTENT(IN) :: z
    REAL(dl) :: om_m
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    om_m=cosm%om_m
    Omega_m_hm=(om_m*(1+z)**3)/Hubble2(z,cosm)

    END FUNCTION Omega_m_hm

    FUNCTION Omega_cold_hm(z,cosm)
    !This calculates omega_cold variations with z (no neutrinos)
    REAL(dl) :: Omega_cold_hm
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    Omega_cold_hm=((cosm%om_c+cosm%om_b)*(1+z)**3)/Hubble2(z,cosm)

    END FUNCTION Omega_cold_hm

    FUNCTION grow(z,cosm)
    !Finds the scale-independent growth fuction at redshift z
    REAL(dl) :: grow
    REAL(dl), INTENT(IN) :: z
    REAL(dl) :: a
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    INTEGER, PARAMETER :: iorder=iorder_growth_interpolation
    INTEGER, PARAMETER :: ifind=ifind_growth_interpolation
    INTEGER, PARAMETER :: imeth=imeth_growth_interpolation

    IF(z==0.) THEN
        grow=1.
    ELSE
        a=1./(1.+z)
        grow=find(a,cosm%a_growth,cosm%growth,cosm%ng,iorder,ifind,imeth)
    END IF

    END FUNCTION grow

    REAL(dl) RECURSIVE FUNCTION ungrow(z, cosm)

    ! Growth function normalised such that g(a) = a at early (matter-dominated) times
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    ungrow = cosm%gnorm*grow(z, cosm)

    END FUNCTION ungrow

    REAL(dl) RECURSIVE FUNCTION acc_growth(z, cosm)

    ! Accumulated growth function: int_0^a g(a)/a da
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: a
    INTEGER, PARAMETER :: iorder=iorder_agrowth_interpolation
    INTEGER, PARAMETER :: ifind=ifind_agrowth_interpolation
    INTEGER, PARAMETER :: imeth=imeth_agrowth_interpolation

    a=1./(1.+z)
    acc_growth = find(a, cosm%a_growth, cosm%agrow, cosm%ng, iorder, ifind, imeth)

    END FUNCTION acc_growth

    FUNCTION sigmaV(R,z,itype,cosm)
    REAL(dl) :: sigmaV
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl), PARAMETER :: tmin=0._dl
    REAL(dl), PARAMETER :: tmax=1._dl
    REAL(dl), PARAMETER :: acc=acc_sigmaV_integration
    INTEGER, PARAMETER :: iorder=iorder_sigmaV_integration

    sigmaV=sqrt(integrate(tmin,tmax,sigmaV_integrand,R,z,itype,cosm,acc,iorder)/3)

    END FUNCTION sigmaV

    FUNCTION sigmaV_integrand(t,R,z,itype,cosm)
    !This is the integrand for the velocity dispersion integral
    REAL(dl) :: sigmaV_integrand
    REAL(dl), INTENT(IN) :: t
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat
    REAL(dl), PARAMETER :: alpha=alpha_sigmaV_integration !Speeds up integral for large 'R'

    IF(t<=0._dl .OR. t>=1._dl) THEN
        sigmaV_integrand=0
    ELSE
        IF(R==0.) THEN
            kR=0
            k=(-1+1/t)**alpha
        ELSE
            kR=(-1+1/t)**alpha
            k=kR/R
        END IF
        w_hat=wk_tophat(kR)
        sigmaV_integrand=(p_lin(k,z,itype,cosm)/k**2)*(w_hat**2)*alpha/(t*(1-t))
    END IF

    END FUNCTION sigmaV_integrand

    FUNCTION neff_integrand(t,R,z,itype,cosm)
    !This is the integrand for the velocity dispersion integral
    REAL(dl) :: neff_integrand
    REAL(dl), INTENT(IN) :: t
    REAL(dl), INTENT(IN) :: R
    REAL(dl), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: itype
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: k, kR, w_hat, w_hat_deriv
    REAL(dl), PARAMETER :: alpha=alpha_neff_integration !Speeds up integral for large 'R'

    IF(t<=0. .OR. t>=1.) THEN
        neff_integrand=0.
    ELSE
        kR=(-1.+1./t)**alpha
        k=kR/R
        w_hat=wk_tophat(kR)
        w_hat_deriv=wk_tophat_deriv(kR)
        neff_integrand=p_lin(k,z,itype,cosm)*w_hat*w_hat_deriv*alpha*kR/(t*(1.-t))
    END IF

    END FUNCTION neff_integrand

    FUNCTION Si(x)

    !Calculates the 'sine integral' function Si(x)
    REAL(dl) :: Si
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: x2, y, f, g, si8
    REAL(dl), PARAMETER :: x0=4. ! Transition between two different approximations

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=x0) THEN

        x2=x*x

        si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

        Si=si8

    ELSE IF(ABS(x)>x0) THEN

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

        Si=pi_HM/2.d0-f*cos(x)-g*sin(x)

    END IF

    END FUNCTION Si

    FUNCTION Ci(x)

    !Calculates the 'cosine integral' function Ci(x)
    REAL(dl) :: Ci
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: x2, y, f, g, ci8
    REAL(dl), PARAMETER :: em_const=0.577215664901532861d0
    REAL(dl), PARAMETER :: x0=4. ! Transition between two different approximations

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IF(ABS(x)<=x0) THEN

        x2=x*x

        ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

        Ci=ci8

    ELSE IF(ABS(x)>x0) THEN

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

    recursive FUNCTION find(x,xtab,ytab,n,iorder,ifind,imeth) result(y)
    !Given two arrays x and y this routine interpolates to find the y_i value at position x_i
    REAL(dl) :: y
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(in) :: x
    REAL(dl), INTENT(IN) :: xtab(n), ytab(n)
    REAL(dl), ALLOCATABLE :: xtab_rev(:), ytab_rev(:)
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

    if (xtab(1)>xtab(n)) then
        WRITE(*,*) 'WARNING: HMCODE find arrays in order. Please report this'
        WRITE(*,*) 'x = ',x, 'n=',n, 'xtab(1)=',xtab(1), 'xtab(n)=',xtab(n)
        call GlobalError('Array order error in HMCode', error_nonlinear)
        CALL reverse(xtab,n, xtab_rev)
        CALL reverse(ytab,n, ytab_rev)
        y =  find(x, xtab_rev, ytab_rev, n, iorder, ifind, imeth)
        return
    end if

    IF(x<xtab(1)) THEN

        !Do a linear interpolation beyond the table boundary
        IF(imeth==1) THEN
            CALL fit_line(a,b,xtab(1),ytab(1),xtab(2),ytab(2))
            y=a*x+b
        ELSE IF(imeth==2) THEN
            y=Lagrange_polynomial(x,1,xtab,ytab)
        ELSE
            ERROR STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(x>xtab(n)) THEN

        !Do a linear interpolation beyond the table boundary

        xarr(1:2) = xtab(n-1:n)
        yarr(1:2) = ytab(n-1:n)

        IF(imeth==1) THEN
            CALL fit_line(a,b,xarr(1),yarr(1),xarr(2),yarr(2))
            y=a*x+b
        ELSE IF(imeth==2) THEN
            y=Lagrange_polynomial(x,1,xarr,yarr)
        ELSE
            ERROR STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(iorder==1) THEN

        IF(n<2) ERROR STOP 'FIND: Not enough points in your table for linear interpolation'

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
            y=a*x+b
        ELSE IF(imeth==2) THEN
            y=Lagrange_polynomial(x,1,xarr,yarr)
        ELSE
            ERROR STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE IF(iorder==2) THEN

        IF(n<3) ERROR STOP 'FIND: Not enough points in your table'

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
                y=a*(x**2)+b*x+c
            ELSE IF(imeth==2) THEN
                y=Lagrange_polynomial(x,2,xarr,yarr)
            ELSE
                ERROR STOP 'FIND: Error, method not specified correctly'
            END IF

        ELSE

            i=table_integer(x,xtab,n,ifind)

            xarr(1:4) = xtab(i-1:i+2)
            yarr(1:4) = ytab(i-1:i+2)

            IF(imeth==1) THEN
                !In this case take the average of two separate quadratic spline values
                CALL fit_quadratic(a,b,c,xarr(1),yarr(1),xarr(2),yarr(2),xarr(3),yarr(3))
                y=(a*x**2+b*x+c)/2.
                CALL fit_quadratic(a,b,c,xarr(2),yarr(2),xarr(3),yarr(3),xarr(4),yarr(4))
                y=y+(a*x**2+b*x+c)/2.
            ELSE IF(imeth==2) THEN
                !In this case take the average of two quadratic Lagrange polynomials
                y=(Lagrange_polynomial(x,2,xarr,yarr)+Lagrange_polynomial(x,2,xarr(2:),yarr(2:)))/2.
            ELSE
                ERROR STOP 'FIND: Error, method not specified correctly'
            END IF

        END IF

    ELSE IF(iorder==3) THEN

        IF(n<4) ERROR STOP 'FIND: Not enough points in your table'

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
            y=a*x**3+b*x**2+c*x+d
        ELSE IF(imeth==2) THEN
            y=Lagrange_polynomial(x,3,xarr,yarr)
        ELSE
            ERROR STOP 'FIND: Error, method not specified correctly'
        END IF

    ELSE

        ERROR STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

    END FUNCTION find

    FUNCTION table_integer(x,xtab,n,imeth)
    !Chooses between ways to find the integer location below some value in an array
    INTEGER :: table_integer
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER, INTENT(IN) :: imeth

    IF(x<xtab(1)) THEN
        table_integer=0
    ELSE IF(x>xtab(n)) THEN
        table_integer=n
    ELSE IF(imeth==1) THEN
        table_integer=linear_table_integer(x,xtab,n)
    ELSE IF(imeth==2) THEN
        table_integer=search_int(x,xtab,n)
    ELSE IF(imeth==3) THEN
        table_integer=int_split(x,xtab,n)
    ELSE
        ERROR STOP 'TABLE INTEGER: Method specified incorrectly'
    END IF

    END FUNCTION table_integer

    FUNCTION linear_table_integer(x,xtab,n)
    !Assuming the table is exactly linear this gives you the integer position
    INTEGER :: linear_table_integer
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    REAL(dl) :: x1, xn

    x1=xtab(1)
    xn=xtab(n)
    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

    END FUNCTION linear_table_integer

    FUNCTION search_int(x,xtab,n)
    !Does a stupid search through the table from beginning to end to find integer
    INTEGER :: search_int
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER :: i

    IF(xtab(1)>xtab(n)) ERROR STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
        IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

    END FUNCTION search_int

    FUNCTION int_split(x,xtab,n)
    !Finds the position of the value in the table by continually splitting it in half
    INTEGER :: int_split
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xtab(n)
    INTEGER :: i1, i2, imid

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
    REAL(dl), INTENT(OUT) :: a0, a1
    REAL(dl), INTENT(IN) :: x1, y1, x2, y2

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

    END SUBROUTINE fit_line

    SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)
    !Given xi, yi i=1,2,3 fits a quadratic between these points
    REAL(dl), INTENT(OUT) :: a0, a1, a2
    REAL(dl), INTENT(IN) :: x1, y1, x2, y2, x3, y3

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2)-a1*x1

    END SUBROUTINE fit_quadratic

    SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
    !Given xi, yi i=1,2,3,4 fits a cubic between these points
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
    REAL(dl) :: Lagrange_polynomial
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: x, xv(n+1), yv(n+1)
    REAL(dl) :: l(n+1)
    INTEGER :: i, j
    real(dl) dx(n+1)

    if (n==3) then
        !Hard coded cubic for speed
        dx = x- xv
        Lagrange_polynomial =  &
            + dx(2)*dx(3)*dx(4)/(xv(1)-xv(2))/(xv(1)-xv(3))/(xv(1)-xv(4)) * yv(1) &
            + dx(1)*dx(3)*dx(4)/(xv(2)-xv(1))/(xv(2)-xv(3))/(xv(2)-xv(4)) * yv(2) &
            + dx(1)*dx(2)*dx(4)/(xv(3)-xv(1))/(xv(3)-xv(2))/(xv(3)-xv(4)) * yv(3) &
            + dx(1)*dx(2)*dx(3)/(xv(4)-xv(1))/(xv(4)-xv(2))/(xv(4)-xv(3)) * yv(4)
        return
    end if

    !Initialise variables, one for sum and one for multiplication
    Lagrange_polynomial=0
    l=1

    !Loops to find the polynomials, one is a sum and one is a multiple
    DO i=0,n
        DO j=0,n
            IF(i .NE. j) l(i+1)=l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
        END DO
        Lagrange_polynomial=Lagrange_polynomial+l(i+1)*yv(i+1)
    END DO

    END FUNCTION Lagrange_polynomial

    SUBROUTINE reverse(arry,n, output)
    !This reverses the contents of arry!
    INTEGER, INTENT(IN) :: n
    REAL(dl), INTENT(IN) :: arry(n)
    INTEGER :: i
    REAL(dl), ALLOCATABLE, intent(OUT) :: output(:)

    ALLOCATE(output(n))
    DO i=1,n
        output(i)=arry(n-i+1)
    END DO

    END SUBROUTINE reverse

    SUBROUTINE fill_growtab(cosm)
    !Fills a table of values of the scale-independent growth function
    TYPE(HM_cosmology) :: cosm
    INTEGER :: i
    REAL(dl) :: a
    REAL(dl), ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL(dl) :: dinit, vinit, zinit, f
    REAL(dl), PARAMETER :: aini=aini_growth_ODE
    REAL(dl), PARAMETER :: afin=afin_growth_ODE
    REAL(dl), PARAMETER :: amin=amin_growth_interpolation
    REAL(dl), PARAMETER :: amax=amax_growth_interpolation
    !REAL(dl), PARAMETER :: ainit=ainit_growth_interpolation
    !REAL(dl), PARAMETER :: amax=amax_growth_interpolation
    INTEGER, PARAMETER :: n=n_growth_interpolation
    REAL(dl), PARAMETER :: acc_ODE=acc_growth_ODE
    INTEGER, PARAMETER :: imeth_ODE=imeth_growth_ODE
    INTEGER, PARAMETER :: iorder_int=iorder_growth_ODE_interpolation
    INTEGER, PARAMETER :: ifind_int=ifind_growth_ODE_interpolation
    INTEGER, PARAMETER :: imeth_int=imeth_growth_ODE_interpolation
    INTEGER, PARAMETER :: iorder_agrow=iorder_integration_agrow

    !These set the initial conditions to be the Om_m=1. growing mode
    !AM Jul 19: changed initial conditions to be appropriate for massive neutrino cosmologies
    !AM Sep 20: changed initial conditions to assume neutrinos cluster, but changed to take into account EDE
    zinit = -1.+1./aini
    f = 1.-Omega_m_hm(zinit, cosm)
    dinit = aini**(1.-3.*f/5.)
    vinit = (1.-3.*f/5.)*aini**(-3.*f/5.)

    IF(HM_verbose) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,aini,afin,dinit,vinit,acc_ODE,imeth_ODE,cosm)
    IF(HM_verbose) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    cosm%gnorm=find(1.d0,a_tab,d_tab,SIZE(a_tab),iorder_int,ifind_int,imeth_int)
    IF(HM_verbose) WRITE(*,*) 'GROWTH: Unnormalised g(a=1):', cosm%gnorm
    d_tab=d_tab/cosm%gnorm

    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_growth)) DEALLOCATE(cosm%a_growth)
    IF(ALLOCATED(cosm%growth)) DEALLOCATE(cosm%growth)

    cosm%ng=n
    ALLOCATE(cosm%a_growth(n),cosm%growth(n))
    DO i=1,n
        a=amin+(amax-amin)*(i-1)/real(n-1,dl)
        cosm%a_growth(i)=a
        cosm%growth(i)=find(a,a_tab,d_tab,SIZE(a_tab),iorder_int,ifind_int,imeth_int)
    END DO

    ! Table integration to calculate G(a)=int_0^a g(a')/a' da'
    IF(ALLOCATED(cosm%agrow)) DEALLOCATE(cosm%agrow)
    ALLOCATE(cosm%agrow(n))

    ! Set to zero, because there is an x=x+y thing later on
    cosm%agrow = 0.
    DO i = 1, n
        ! Do the integral using the arrays
        IF (i > 1) THEN
            cosm%agrow(i) = inttab(cosm%a_growth, cosm%gnorm*cosm%growth/cosm%a_growth, 1, i, iorder_agrow)
        END IF
        ! Add missing section; g(a=0)/0 = 1, so you just add on a rectangle of height g*a/a=g
        cosm%agrow(i) = cosm%agrow(i)+cosm%gnorm*cosm%growth(1)
    END DO

    IF(HM_verbose) WRITE(*,*) 'GROWTH: Accumulated G(a=1):', cosm%agrow(n)
    IF(HM_verbose) WRITE(*,*) 'GROWTH: Done'
    IF(HM_verbose) WRITE(*,*)

    END SUBROUTINE fill_growtab

    SUBROUTINE ode_growth(x,v,t,ti,tf,xi,vi,acc,imeth,cosm)
    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    REAL(dl) :: xi, ti, tf, dt, acc, vi, x4, v4, t4
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
        dt=(tf-ti)/real(n-1,dl)

        !Intially fix this to zero. It will change to 1 if method is a 'failure'
        ifail=0

        DO i=1,n-1

            x4=real(x8(i))
            v4=real(v8(i))
            t4=real(t8(i))

            IF(imeth==1) THEN

                !Crude method
                kx1=dt*fd(v4)
                kv1=dt*fv(x4,v4,t4,cosm)

                x8(i+1)=x8(i)+kx1
                v8(i+1)=v8(i)+kv1

            ELSE IF(imeth==2) THEN

                !Mid-point method
                !2017/06/18 - There was a bug in this part before. Luckily it was not used. Thanks Dipak Munshi.
                kx1=dt*fd(v4)
                kv1=dt*fv(x4,v4,t4,cosm)
                kx2=dt*fd(v4+kv1/2.)
                kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.,cosm)

                x8(i+1)=x8(i)+kx2
                v8(i+1)=v8(i)+kv2

            ELSE IF(imeth==3) THEN

                !4th order Runge-Kutta method (fast!)
                kx1=dt*fd(v4)
                kv1=dt*fv(x4,v4,t4,cosm)
                kx2=dt*fd(v4+kv1/2.)
                kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,t4+dt/2.,cosm)
                kx3=dt*fd(v4+kv2/2.)
                kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,t4+dt/2.,cosm)
                kx4=dt*fd(v4+kv3)
                kv4=dt*fv(x4+kx3,v4+kv3,t4+dt,cosm)

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

    FUNCTION fv(d,v,a,cosm)

    !v'=f(v) in ODE solver
    REAL(dl) :: fv
    REAL(dl), INTENT(IN) :: d, v, a
    REAL(dl) :: f1, f2, z, Om_m
    TYPE(HM_cosmology), INTENT(IN) :: cosm

    z=-1.+(1./a)

    !AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    IF (cold_growth) THEN
        Om_m = Omega_cold_hm(z,cosm)
    ELSE
        Om_m = Omega_m_hm(z,cosm)
    END IF
    f1=3.*Om_m*d/(2.*a**2.)
    f2=(2.+AH(z,cosm)/Hubble2(z,cosm))*(v/a)

    fv=f1-f2

    END FUNCTION fv

    FUNCTION fd(v)

    !d'=f(d) in ODE solver
    REAL(dl) :: fd
    REAL(dl), INTENT(IN) :: v

    fd=v

    END FUNCTION fd

    REAL(dl) FUNCTION dc_Mead(z, cosm)

    ! delta_c fitting function from Mead (2017; 1606.05345)
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: lg, bG, Om_m, a

    ! See Appendix A of Mead (2017) for naming convention
    REAL(dl), PARAMETER :: p10 = -0.0069
    REAL(dl), PARAMETER :: p11 = -0.0208
    REAL(dl), PARAMETER :: p12 = 0.0312
    REAL(dl), PARAMETER :: p13 = 0.0021
    INTEGER, PARAMETER :: a1 = 1
    REAL(dl), PARAMETER :: p20 = 0.0001
    REAL(dl), PARAMETER :: p21 = -0.0647
    REAL(dl), PARAMETER :: p22 = -0.0417
    REAL(dl), PARAMETER :: p23 = 0.0646
    INTEGER, PARAMETER :: a2 = 0
    REAL(dl), PARAMETER :: dc0 = (3./20.)*(12.*pi_HM)**(2./3.)

    IF (cold_growth) STOP 'DV_MEAD: Error, this will not work if you want cold growth'

    lg = ungrow(z, cosm)
    bG = acc_growth(z, cosm)
    Om_m = Omega_m_hm(z, cosm)
    a = 1./(1.+z)

    dc_Mead = 1.
    dc_Mead = dc_Mead+f_Mead(lg/a, bG/a, p10, p11, p12, p13)*log10(Om_m)**a1
    dc_Mead = dc_Mead+f_Mead(lg/a, bG/a, p20, p21, p22, p23)
    dc_Mead = dc_Mead*dc0*(1.-0.041*cosm%f_nu)

    END FUNCTION dc_Mead

    REAL(dl) FUNCTION Dv_Mead(z, cosm)

    ! Delta_v fitting function from Mead (2017; 1606.05345)
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: z
    TYPE(HM_cosmology), INTENT(IN) :: cosm
    REAL(dl) :: lg, bG, Om_m, a

    ! See Appendix A of Mead (2017) for naming convention
    REAL(dl), PARAMETER :: p30 = -0.79
    REAL(dl), PARAMETER :: p31 = -10.17
    REAL(dl), PARAMETER :: p32 = 2.51
    REAL(dl), PARAMETER :: p33 = 6.51
    INTEGER, PARAMETER :: a3 = 1
    REAL(dl), PARAMETER :: p40 = -1.89
    REAL(dl), PARAMETER :: p41 = 0.38
    REAL(dl), PARAMETER :: p42 = 18.8
    REAL(dl), PARAMETER :: p43 = -15.87
    INTEGER, PARAMETER :: a4 = 2
    REAL(dl), PARAMETER :: Dv0 = 18.*pi_HM**2

    IF (cold_growth) STOP 'DV_MEAD: Error, this will not work if you want cold growth'

    lg = ungrow(z, cosm)
    bG = acc_growth(z, cosm)
    Om_m = Omega_m_hm(z, cosm)
    a = 1./(1.+z)

    Dv_Mead = 1.
    Dv_Mead = Dv_Mead+f_Mead(lg/a, bG/a, p30, p31, p32, p33)*log10(Om_m)**a3
    Dv_Mead = Dv_Mead+f_Mead(lg/a, bG/a, p40, p41, p42, p43)*log10(Om_m)**a4
    Dv_Mead = Dv_Mead*Dv0*(1.+0.763*cosm%f_nu)

    END FUNCTION Dv_Mead

    REAL(dl) FUNCTION f_Mead(x, y, p0, p1, p2, p3)

    ! Equation A3 in Mead (2017)
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: x, y
    REAL(dl), INTENT(IN) :: p0, p1, p2, p3

    f_Mead = p0+p1*(1.-x)+p2*(1.-x)**2+p3*(1.-y)

    END FUNCTION f_Mead

    !!AM End HMcode

    subroutine PKequal(State,redshift,w_lam,wa_ppf,w_hf,wa_hf)
    !used by halofit_casarini: arXiv:0810.0190, arXiv:1601.07230
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
        if (abs(error) <= 1d-7) exit
        w_lam=w_lam*(1+error)**10.d0
    enddo
    w_hf=w_lam
    wa_hf=0._dl
    w_lam=w_true
    wa_ppf=wa_true
    write(*,*)'at z = ',real(redshift),' equivalent w_const =', real(w_hf)

    end subroutine PKequal


    end module NonLinear
