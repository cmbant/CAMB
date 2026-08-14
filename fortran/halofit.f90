    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    ! AM Mar 16: Added in HMcode
    ! AM May 16: Fixed some small bugs and added better neutrino approximations
    ! AL Jun16: put in partial openmp for HMcode (needs restructure to do properly)
    ! AM Sep 16: Attempted fix of strange bug. No more modules with unallocated arrays as inputs
    ! LC Oct 16: extended Halofit from w=const. models to w=w(a) with PKequal
    ! AM May 17: Made the baryon feedback parameters more obvious in HMcode
    ! AL Jul 17: fixed undefined z calling Tcb_Tcbnu_ratio
    ! AM Jul 17: sped-up HMcode integration routines
    ! AM May 18: Fixed bug in Dolag correction to c(M) power
    ! AM Jul 19: Upgraded accuracy and bug fix for massive-neutrino models
    ! AL Jul 19: Speedups, use linear interpolation for pk; find index using fixed spacing; precompute growth(z)
    ! AL Sep 19: Propagate errors rather than stop, decrease jmax for integration time out (prevent very slow error)
    ! AM Sep 20: Added HMcode-2020 model
    ! AM Jan 23: Fixed HMcode-2020 feedback low-k predictions
    ! AL Jul 26: Code optimizations and cleanups
    ! AL Jul 26: Fixed halofit_casarini; PKequal had been a no-op since w_lam/wa_ppf stopped
    !           being global variables, so it silently returned the takahashi result
    ! AL Aug 26: Only extract the HMcode BAO wiggle for the 2020 versions that use it, and only once
    !           (growth-scaled from z=0) for cosmologies with nearly scale-independent growth
    !           (see HMcode_wiggle_max_fnu, DarkEnergy%assume_scale_indep_lowz_growth); non-linear ratios
    !           for all matter power redshifts are now obtained in one call
    ! AL Aug 26: Parallelize over redshift rather than inside each redshift (no per-redshift state is
    !           kept on THalofit any more, so the loops are re-entrant); re-use the redshift set-up
    !           between the passes of the HMcode-2020 feedback response. Results are unchanged.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module NonLinear
    use results
    use DarkEnergyInterface
    use classes
    use Transfer
    use constants
    use config
    use MiscUtils, only: DefaultFalse
    implicit none
    private

    integer, parameter :: halofit_original = 1, halofit_bird = 2, halofit_peacock = 3, halofit_takahashi = 4
    integer, parameter :: halofit_casarini = 7
    integer, parameter :: halofit_mead2016 = 5, halofit_halomodel = 6, halofit_mead2015 = 8, halofit_mead2020 = 9
    integer, parameter :: halofit_mead2020_feedback = 10
    integer, parameter :: halofit_mead = halofit_mead2016 ! AM Kept for backwards compatibility
    integer, parameter :: halofit_default = halofit_mead2020

    logical :: HM_verbose = .false.

    ! Number of points used for the (Gaussian-filtered) sigma(R) integral in standard halofit
    integer, parameter :: nint_wint = 3000
    ! The integral is done in t=1/(1+k), sampled at the mid-points of nint_wint equal intervals
    ! in t; the wavenumbers are therefore fixed, and shared by wint_pk_table and wint
    integer :: i_wint ! only used to fill wint_k below
    real(dl), parameter :: wint_k(nint_wint) = [(nint_wint/(i_wint - 0.5_dl) - 1, i_wint=1, nint_wint)]
    ! Tolerance on sigma(R_nl)=1 defining the standard halofit non-linear scale. Older versions
    ! used 1e-3, which left R_nl (and hence P_nl at the 0.3% level) dependent on the solver's iterates
    real(dl), parameter :: tol_r_nl = 1e-4_dl

    type :: THalofit_zparams
        ! Per-redshift background quantities used by the standard halofit fitting formulae.
        ! Held in a local variable (not on THalofit) so that redshifts can be done in parallel.
        real(dl) :: om_m = 1._dl, om_v = 0._dl, fnu = 0._dl, acur = 1._dl
        real(dl) :: w_hf = -1._dl, wa_hf = 0._dl
    end type THalofit_zparams

    type :: TEquivalent_wCDM
        ! Workspace for PKequal (halofit_casarini); describes a trial model with constant dark
        ! energy equation of state w, but otherwise the same densities as State.
        class(CAMBdata), pointer :: State => null()
        real(dl) :: w = -1._dl ! equation of state of the trial model
        real(dl) :: a_star = 0._dl, a_z = 0._dl ! integration limits (last scattering, redshift of interest)
        real(dl) :: dlsb = 0._dl ! comoving distance between them in the actual model
        real(dl) :: tol = 1e-7_dl ! relative tolerance of the distance integrals
    end type TEquivalent_wCDM

    type, extends(TNonLinearModel) :: THalofit
        integer :: halofit_version = halofit_default
        !!TT - These are the baryon parameters of HMCode
        real(dl) :: HMcode_A_baryon = 3.13_dl
        real(dl) :: HMcode_eta_baryon = 0.603_dl
        real(dl) :: HMcode_logT_AGN = 7.8_dl
        ! Condition for extracting the HMcode-2020 BAO wiggle only once (see assign_HM_cosmology);
        ! relaxing it is faster, at the cost of assuming growth is more nearly scale-independent
        real(dl) :: HMcode_wiggle_max_fnu = 0.01_dl ! largest neutrino fraction for which the wiggle is re-used
        ! Note THalofit holds no per-redshift state: it is read-only while GetNonLinRatios runs, so
        ! that the loop over redshifts can be parallelized (HMcode's imead lives in HM_tables%imead,
        ! and the standard halofit background quantities in a local THalofit_zparams).
    contains
    procedure :: ReadParams => THalofit_ReadParams
    procedure :: GetNonLinRatios => THalofit_GetNonLinRatios
    procedure :: halofit
    procedure :: HMcode
    procedure, nopass :: PythonClass => THalofit_PythonClass
    procedure, nopass :: SelfPointer => THalofit_SelfPointer
    procedure, nopass, private :: Delta_v
    procedure, nopass, private :: delta_c
    procedure, nopass, private :: eta
    procedure, nopass, private :: kstar
    procedure, nopass, private :: As
    procedure, private :: conc_bull
    procedure, nopass, private :: fdamp
    procedure, nopass, private :: p_1h
    procedure, nopass, private :: p_2h
    procedure, nopass, private :: alpha
    procedure, private :: halomod
    procedure, private :: halomod_init
    procedure, private :: write_parameters
    procedure, private :: zcoll_bull
    end type

    public THalofit, HM_verbose
    public halofit_default, halofit_original, halofit_bird, halofit_peacock, halofit_takahashi
    public halofit_mead2016, halofit_mead2015, halofit_mead2020, halofit_halomodel, halofit_casarini
    public halofit_mead2020_feedback
    public halofit_mead ! AM for backwards compatibility

    type HM_cosmology
        ! Contains only things that do not need to be recalculated with each new z
        real(dl) :: om_m, om_c, om_b, om_nu, om_v, w, wa, f_nu, ns, h, Tcmb, Nnu
        real(dl), allocatable :: log_r_sigma(:), log_sigma(:)
        real(dl), allocatable :: a_growth(:), growth(:), agrow(:)
        real(dl), allocatable :: log_k_plin(:), log_plin(:), log_plinc(:)
        real(dl), allocatable :: log_k_wiggle(:), pk_wiggle(:)
        real(dl) :: kmax
        real(dl) :: gnorm
        integer :: nk, ng, nsig
        integer :: plin_iz = 0 ! redshift index currently stored in the linear power table (0 if none)
        real(dl) :: grow_z, this_z ! cached growth factor at redshift being calculated
        logical :: wiggle_once = .false. ! extract the BAO wiggle once, at z=0
        logical :: wiggle_each_z = .false. ! extract the BAO wiggle separately at each redshift
        ! AM - Added feedback parameters below at fixed fiducial (DMONLY) values
        real(dl) :: A_baryon = 3.13
        real(dl) :: eta_baryon = 0.603
        real(dl) :: logT_AGN = 7.8
    end type HM_cosmology

    type HM_tables
        ! Stuff that needs to be recalculated for each new z
        real(dl), allocatable :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
        real(dl), allocatable :: rs(:), nfw_norm(:), nfw_moment2(:) ! cached halo-profile quantities
        real(dl), allocatable :: p1h_quad_weight(:), nu_eta(:), baryon_mass_fraction(:)
        real(dl) :: sigv, sigv100, knl, rnl, neff, sig8z, z, dc, sig8z_cold
        real(dl) :: eta_hm, kstar_hm, alpha_hm, fdamp_hm, kdamp_hm, one_minus_fnu_sq, f_star_hm
        real(dl) :: dolag_inf = 1, dolag_z = 1 ! cached Dolag (2004) concentration corrections
        integer :: n
        integer :: imead = -1 ! HMcode variant these HM_tables were filled for (see HMcode)
    end type HM_tables
    !!AM - End of my additions

    ! HMcode parameters
    real(dl), parameter :: zc_Dolag = 10._dl ! Halo collapse redshift for Dolag
    real(dl), parameter :: fdamp_min = 1e-3_dl ! Minimum value of fdamp
    real(dl), parameter :: fdamp_max = 0.99_dl ! Maximum value of fdamp
    real(dl), parameter :: alpha_min = 0.5_dl ! Minimum value of alpha transition
    real(dl), parameter :: alpha_max = 2._dl ! Maximum value of alpha transition
    real(dl), parameter :: ks_limit = 7._dl ! Limit for (k/ks)^2 in one-halo term
    real(dl), parameter :: pi_HM = const_pi ! Lovely pi

    ! HMcode linear P(k) numerical parameters
    ! AM: Jul 19: Updated nk_pk_interpolation from 128 to 512
    ! AM: Dec 20: Calculation time and accuracy are especially sensive to these parameters
    real(dl), parameter :: kmin_pk_interpolation = 1d-3 ! Minimum wavenumber if rebinning [h/Mpc]
    real(dl), parameter :: kmax_pk_interpolation = 1d2 ! Maximum wavenumber if rebinning [h/Mpc]
    integer, parameter :: nk_pk_interpolation = 512 ! Number of points in k if rebining
    logical, parameter :: plin_extrap = .false. ! Extrapolate at high-k via thoery or simple power law

    ! HMcode dewiggle numerical parameters
    real(dl), parameter :: kmin_wiggle = 5e-3_dl ! Minimum wavenumber to calculate wiggle [Mpc/h]
    real(dl), parameter :: kmax_wiggle = 5._dl ! Maximum wavenumber to calculate wiggle [Mpc/h]
    integer, parameter :: nk_wiggle = 512 ! Number of k points to store wiggle
    real(dl), parameter :: wiggle_sigma = 0.25_dl ! Smoothing width if using Gaussian smoothing
    ! Wavenumber where linear and nowiggle, forced to be identical [Mpc/h]
    real(dl), parameter :: knorm_nowiggle = 0.03_dl

    ! Linear growth integral numerical parameters (LCDM only; only used in Dolag correction)
    ! AM: Jul 19: Updated acc_growint from 1e-3 to 1e-4
    ! AM: Sep 20: Changed cold_growth = .FALSE. to be in line with my code
    real(dl), parameter :: acc_growth_integration = 1e-4 ! Accuracy for growth function integral
    logical, parameter :: cold_growth = .false. ! Should growth be of cold or all matter?

    ! Linear growth factor tabulation and interpolation numerical parameters
    real(dl), parameter :: amin_growth_interpolation = 1e-3 ! Minimum scale factor for growth interpolation
    real(dl), parameter :: amax_growth_interpolation = 1. ! Maximum scale factor for growth interpolation
    integer, parameter :: n_growth_interpolation = 64 ! Number of entries for growth look-up table

    ! Growth function ODE numerical parameters
    ! AM: Sep 20: Changed aini from 1e-3 to 1e-4
    real(dl), parameter :: aini_growth_ODE = 1e-4 ! Initial scale factor for growth ODE
    integer, parameter :: nsub_growth_ODE = 8 ! RK4 steps per output point of the growth table

    ! HMcode numerical parameters for sigma(R) tabulation and interpolation
    real(dl), parameter :: rmin_sigma_interpolation = 1e-4 ! Minimum scale for sigma(R) look-up tables [Mpc/h]
    real(dl), parameter :: rmax_sigma_interpolation = 1e3 ! Maximum scale for sigma(R) look-up tables [Mpc/h]
    integer, parameter :: n_sigma_interpolation = 64 ! Number of points in look-up tables

    ! HMcode numerical parameters for sigma(R) integration (dominates run time as of Jul 2019)
    ! AM: Jul 19: Upgraded acc_sigma from 1e-3 to 3e-4
    ! AM: Sep 20: Upgraded acc_sigma from 3e-4 to 1e-4 to fix problems for some cosmologies
    real(dl), parameter :: acc_sigma_integration = 1e-4 ! Relative accuracy of numerical integration
    real(dl), parameter :: alpha_sigma_integration = 3. ! Exponent to speed up integration

    ! HMcode numerical parameters for sigmaV(R) integration
    ! AM: Jul 19: Upgraded acc_sigmaV from 1e-3 to 1e-4
    real(dl), parameter :: acc_sigmaV_integration = 1e-4 ! Relative accuracy of numerical integration
    real(dl), parameter :: alpha_sigmaV_integration = 3. ! Exponent to speed up integration

    ! HMcode numerical parameters for neff(R) integration
    real(dl), parameter :: acc_neff_integration = 1e-4 ! Relative accuracy of numerical integration
    real(dl), parameter :: alpha_neff_integration = 2. ! Exponent to speed up integration

    ! HMcode numerical parameters for cold transfer function approximation
    ! AM: Sep 20: Care here, before EdS_Tcold_growth=.TRUE.
    logical, parameter :: EdS_Tcold_growth = .false. ! Should the EdS growth function (incorrectly) be used?

    contains

    function THalofit_PythonClass()
    character(len=:), allocatable :: THalofit_PythonClass

    THalofit_PythonClass = 'Halofit'

    end function THalofit_PythonClass

    subroutine THalofit_SelfPointer(cptr, P)
    use iso_c_binding
    type(c_ptr) :: cptr
    type(THalofit), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine THalofit_SelfPointer

    subroutine THalofit_ReadParams(this, Ini)
    use IniObjects
    class(THalofit) :: this
    class(TIniFile), intent(in) :: Ini

    this%halofit_version = Ini%Read_Int('halofit_version', halofit_default)
    if (this%halofit_version == halofit_mead2015 .or. this%halofit_version == halofit_mead2016) then
        this%HMcode_A_baryon = Ini%Read_Double('HMcode_A_baryon', 3.13_dl)
        this%HMcode_eta_baryon = Ini%Read_Double('HMcode_eta_baryon', 0.603_dl)
    else if (this%halofit_version == halofit_mead2020_feedback) then
        this%HMcode_logT_AGN = Ini%Read_Double('HMcode_logT_AGN', 7.8_dl)
    end if
    this%HMcode_wiggle_max_fnu = Ini%Read_Double('HMcode_wiggle_max_fnu', this%HMcode_wiggle_max_fnu)

    end subroutine THalofit_ReadParams

    subroutine THalofit_GetNonLinRatios(this, State, CAMB_Pk)
    ! Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    ! for each redshift and wavenumber
    ! This implementation uses Halofit
    !$ use omp_lib, only: omp_set_num_threads
    class(THalofit) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    integer itf
    real(dl) a, plin, pq, ph, pnl, rk
    real(dl) sig, rknl, rneff, rncur, d1, d2
    real(dl) diff, xlogr1, xlogr2, xlogr, step, last_step, rmid, h2
    real(dl) w_eff, wa_eff, omm0, fnu
    real(dl), allocatable :: pk_fac(:)
    type(THalofit_zparams) :: hf
    integer i, index_cache
    logical found_nonlinear_scale

    !$ if (ThreadNum /= 0) call OMP_SET_NUM_THREADS(ThreadNum)

    select type (State)
    class is (CAMBdata)
        associate(Params => State%CP)

            if (this%halofit_version == halofit_mead2016 .or. &
                this%halofit_version == halofit_halomodel .or. &
                this%halofit_version == halofit_mead2015 .or. &
                this%halofit_version == halofit_mead2020 .or. &
                this%halofit_version == halofit_mead2020_feedback) then
                call this%HMcode(State, CAMB_Pk)
            else

                !!BR09 putting neutrinos into the matter as well, not sure if this is correct, but at least one
                !!will get a consistent omk.
                h2 = (Params%H0/100)**2
                omm0 = (Params%omch2 + Params%ombh2 + Params%omnuh2)/h2
                fnu = Params%omnuh2/h2/omm0

                CAMB_Pk%nonlin_ratio = 1

                call Params%DarkEnergy%Effective_w_wa(w_eff, wa_eff)

                ! Each redshift is independent; hf and pk_fac hold all the per-redshift state
                !$OMP PARALLEL DO IF(CAMB_Pk%num_z > 1), DEFAULT(SHARED), SCHEDULE(DYNAMIC) &
                !$OMP PRIVATE(itf, hf, pk_fac, a, i, rk, plin, pnl, pq, ph, index_cache), &
                !$OMP PRIVATE(sig, rknl, rneff, rncur, d1, d2, diff, xlogr1, xlogr2, xlogr, step, last_step), &
                !$OMP PRIVATE(rmid, found_nonlinear_scale)
                do itf = 1, CAMB_Pk%num_z
                    if (global_error_flag /= 0) cycle
                    if (.not. allocated(pk_fac)) allocate(pk_fac(nint_wint))

                    hf%fnu = fnu
                    hf%w_hf = w_eff
                    hf%wa_hf = wa_eff
                    if (this%halofit_version == halofit_casarini) then
                        ! calculate equivalent w-constant models (w_hf,0) for w_lam+wa_ppf(1-a) models
                        ! [Casarini+ (2009,2016)].
                        call PKequal(State, CAMB_Pk%redshifts(itf), w_eff, wa_eff, hf%w_hf, hf%wa_hf)
                        if (global_error_flag /= 0) cycle
                    end if

                    ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and
                    ! curvature (rncur) of the power spectrum at the desired redshift, using method
                    ! described in Smith et al (2002).
                    a = 1/real(1 + CAMB_Pk%redshifts(itf), dl)
                    call omegas_hf(a, omm0, State%Omega_de, hf%w_hf, hf%wa_hf, hf%om_m, hf%om_v)
                    hf%acur = a
                    call wint_pk_table(CAMB_Pk, itf, pk_fac)

                    ! Solve sigma(R)=1 for the Gaussian filter scale R. wint already returns
                    ! d1 = dln sigma^2/dln R, so Newton's method on log10(sigma^2)=0 needs no extra
                    ! derivative (a log-derivative is the same in any base). The step is safeguarded
                    ! by the bracket [xlogr1,xlogr2] in log10(R), in which sigma crosses one.
                    xlogr1 = -2.0
                    xlogr2 = 3.5
                    xlogr = (xlogr1 + xlogr2)/2
                    step = xlogr2 - xlogr1
                    found_nonlinear_scale = .true.
                    do
                        rmid = 10**xlogr
                        call wint(pk_fac, rmid, sig, d1, d2)
                        diff = sig - 1.0
                        ! Stop once sigma=1 to tolerance, or if the bracket cannot be narrowed further
                        if (abs(diff) <= tol_r_nl .or. xlogr2 - xlogr1 < 1e-10) exit
                        if (diff > 0) then
                            xlogr1 = xlogr ! sigma>1, so the non-linear scale is at larger R
                        else
                            xlogr2 = xlogr
                        end if
                        if (xlogr2 < -1.9999) then
                            ! is still linear, exit
                            found_nonlinear_scale = .false.
                            exit
                        else if (xlogr1 > 3.4999) then
                            ! Totally crazy non-linear
                            call GlobalError('Error in halofit (xlogr1>3.4999)', error_nonlinear)
                            found_nonlinear_scale = .false.
                            exit
                        end if
                        ! Newton step, bisecting instead if it leaves the bracket or does not at
                        ! least halve the step size (which guarantees the bracket collapses)
                        last_step = step
                        step = -2*log10(sig)/d1
                        if (xlogr + step <= xlogr1 .or. xlogr + step >= xlogr2 .or. &
                            abs(step) > abs(last_step)/2) step = (xlogr1 + xlogr2)/2 - xlogr
                        xlogr = xlogr + step
                    end do

                    if (.not. found_nonlinear_scale) cycle
                    rknl = 1./rmid
                    rneff = -3 - d1
                    rncur = -d2

                    ! now calculate power spectra for a logarithmic range of wavenumbers (rk)

                    index_cache = 1
                    do i = 1, CAMB_Pk%num_k
                        rk = exp(CAMB_Pk%log_kh(i))

                        if (rk > this%Min_kh_nonlinear) then

                            ! linear power spectrum !! Remember => plin = k^3 * P(k) * constant
                            ! constant = 4*pi*V/(2*pi)^3

                            plin = MatterPowerData_k(CAMB_Pk, rk, itf, index_cache)*(rk**3/(2*const_pi**2))

                            ! calculate nonlinear power according to halofit: pnl = pq + ph,
                            ! where pq represents the quasi-linear (halo-halo) power and
                            ! where ph is represents the self-correlation halo term.

                            call this%halofit(hf, rk, rneff, rncur, rknl, plin, pnl, pq, ph) ! halo fitting formula
                            CAMB_Pk%nonlin_ratio(i, itf) = sqrt(pnl/plin)

                        end if

                    end do
                end do
                !$OMP END PARALLEL DO

            end if
        end associate
    end select

    end subroutine THalofit_GetNonLinRatios

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine halofit(this, hf, rk, rn, rncur, rknl, plin, pnl, pq, ph)
    class(THalofit) :: this
    type(THalofit_zparams), intent(in) :: hf
    real(dl) gam, a, b, c, xmu, xnu, alpha, beta, f1, f2, f3
    real(dl) rk, rn, plin, pnl, pq, ph, plinaa
    real(dl) rknl, y, rncur
    real(dl) f1a, f2a, f3a, f1b, f2b, f3b, frac
    real(dl) extragam, peacock_fudge

    if (this%halofit_version == halofit_original .or. this%halofit_version == halofit_bird &
        .or. this%halofit_version == halofit_peacock) then
        ! halo model nonlinear fitting formula as described in
        ! Appendix C of Smith et al. (2002)
        ! SPB11: Standard halofit underestimates the power on the smallest scales by a
        ! factor of two. Add an extra correction from the simulations in Bird, Viel,
        ! Haehnelt 2011 which partially accounts for this.
        if (this%halofit_version == halofit_bird) then
            extragam = 0.3159 - 0.0765*rn - 0.8350*rncur
            gam = extragam + 0.86485 + 0.2989*rn + 0.1631*rncur
        else
            gam = 0.86485 + 0.2989*rn + 0.1631*rncur
        end if
        a = 1.4861 + 1.83693*rn + 1.67618*rn*rn + 0.7940*rn*rn*rn + &
            0.1670756*rn*rn*rn*rn - 0.620695*rncur
        a = 10**a
        b = 10**(0.9463 + 0.9466*rn + 0.3084*rn*rn - 0.940*rncur)
        c = 10**(-0.2807 + 0.6669*rn + 0.3214*rn*rn - 0.0793*rncur)
        xmu = 10**(-3.54419 + 0.19086*rn)
        xnu = 10**(0.95897 + 1.2857*rn)
        alpha = 1.38848 + 0.3701*rn - 0.1452*rn*rn
        beta = 0.8291 + 0.9854*rn + 0.3400*rn**2 + hf%fnu*(-6.4868 + 1.4373*rn**2)
    else if (this%halofit_version == halofit_takahashi .or. this%halofit_version == halofit_casarini) then
        ! RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
        ! than latest N-body simulations at small scales.
        ! Update the following fitting parameters of gam,a,b,c,xmu,xnu,
        ! alpha & beta from the simulations in Takahashi+ 2012.
        ! The improved halofit accurately provide the power spectra for WMAP
        ! cosmological models with constant w.
        ! LC16 Jun: Casarini+ 2009,2016 extended constant w prediction for w(a).
        gam = 0.1971 - 0.0843*rn + 0.8460*rncur
        a = 1.5222 + 2.8553*rn + 2.3706*rn*rn + 0.9903*rn*rn*rn + &
            0.2250*rn*rn*rn*rn - 0.6038*rncur + 0.1749*hf%om_v*(1. + hf%w_hf + hf%wa_hf*(1 - hf%acur))
        a = 10**a
        b = 10**(-0.5642 + 0.5864*rn + 0.5716*rn*rn - 1.5474*rncur + &
            0.2279*hf%om_v*(1. + hf%w_hf + hf%wa_hf*(1 - hf%acur)))
        c = 10**(0.3698 + 2.0404*rn + 0.8161*rn*rn + 0.5869*rncur)
        xmu = 0.
        xnu = 10**(5.2105 + 3.6902*rn)
        alpha = abs(6.0835 + 1.3373*rn - 0.1959*rn*rn - 5.5274*rncur)
        beta = 2.0379 - 0.7354*rn + 0.3157*rn**2 + 1.2490*rn**3 + &
            0.3980*rn**4 - 0.1682*rncur + hf%fnu*(1.081 + 0.395*rn**2)
    else
        call MpiStop('Unknown halofit_version')
    end if

    if (abs(1 - hf%om_m) > 0.01) then ! omega evolution
        f1a = hf%om_m**(-0.0732)
        f2a = hf%om_m**(-0.1423)
        f3a = hf%om_m**(0.0725)
        f1b = hf%om_m**(-0.0307)
        f2b = hf%om_m**(-0.0585)
        f3b = hf%om_m**(0.0743)
        frac = hf%om_v/(1. - hf%om_m)
        f1 = frac*f1b + (1 - frac)*f1a
        f2 = frac*f2b + (1 - frac)*f2a
        f3 = frac*f3b + (1 - frac)*f3a
    else
        f1 = 1.0
        f2 = 1.
        f3 = 1.
    end if

    y = (rk/rknl)

    ph = a*y**(f1*3)/(1 + b*y**(f2) + (f3*c*y)**(3 - gam))
    ph = ph/(1 + xmu*y**(-1) + xnu*y**(-2))*(1 + hf%fnu*0.977)
    plinaa = plin*(1 + hf%fnu*47.48*rk**2/(1 + 1.5*rk**2))
    pq = plin*(1 + plinaa)**beta/(1 + plinaa*alpha)*exp(-y/4.0 - y**2/8.0)

    pnl = pq + ph

    if (this%halofit_version == halofit_peacock) then
        ! From http://www.roe.ac.uk/~jap/haloes/
        ! (P-P_linear) -> (P-P_linear) * (1+2y^2)/(1+y^2), where y = k/10 h Mpc^(-1).
        peacock_fudge = rk/10
        pnl = plin + (pnl - plin)*(1 + 2*peacock_fudge**2)/(1 + peacock_fudge**2)
    end if

    end subroutine halofit

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The subroutine wint, finds the effective spectral quantities
    ! rknl, rneff & rncur. This it does by calculating the radius of
    ! the Gaussian filter at which the variance is unity = rknl.
    ! rneff is defined as the first derivative of the variance, calculated
    ! at the nonlinear wavenumber and similarly the rncur is the second
    ! derivative at the nonlinear wavenumber.

    subroutine wint_pk_table(CAMB_Pk, itf, fac)
    ! Tabulates the part of the wint integrand that does not depend on the filter scale r,
    ! so that the search for r_nl does not repeat the (expensive) power spectrum look-ups.
    ! With k=1/t-1 the measure dk/dt = (1+k)^2, so k^3 dk/(k dt^2) is k^2 (1+k)^2
    type(MatterPowerData) :: CAMB_Pk
    integer, intent(in) :: itf
    real(dl), intent(out) :: fac(nint_wint)
    real(dl) rk, anorm
    integer i, index_cache

    index_cache = 1
    anorm = 1/(2*const_pi**2)
    do i = 1, nint_wint
        rk = wint_k(i)
        fac(i) = MatterPowerData_k(CAMB_Pk, rk, itf, index_cache)*(anorm*(rk*(1 + rk))**2)
    end do

    end subroutine wint_pk_table

    subroutine wint(fac, r, sig, d1, d2)
    ! sigma(r) for a Gaussian filter, with d1 = dln sigma^2/dln r and d2 = d2ln sigma^2/dln r^2
    real(dl), intent(in) :: fac(nint_wint)
    real(dl), intent(in) :: r
    real(dl), intent(out) :: sig, d1, d2
    real(dl) sum1, sum2, sum3, x2, w1
    integer i

    sum1 = 0.d0
    sum2 = 0.d0
    sum3 = 0.d0
    do i = 1, nint_wint
        x2 = (wint_k(i)*r)**2
        w1 = exp(-x2)*fac(i)
        sum1 = sum1 + w1
        sum2 = sum2 + x2*w1
        sum3 = sum3 + x2*(1 - x2)*w1
    end do
    sig = sqrt(sum1/nint_wint)
    d1 = -2*sum2/sum1
    d2 = -d1*d1 - 4*sum3/sum1

    end subroutine wint

    !!JD 08/13 generalize to variable w

    subroutine omegas_hf(aa, om_m0, om_v0, wval, waval, om_m, om_v)
    ! Evolution of omega_m and omega_lambda with expansion factor, generalized to w0-wa
    real(dl), intent(in) :: aa, om_m0, om_v0, wval, waval
    real(dl), intent(out) :: om_m, om_v
    real(dl) omega_t, Qa2

    Qa2 = aa**(-1.0 - 3.0*(wval + waval))*exp(-3.0*(1 - aa)*waval)
    omega_t = 1.0 + (om_m0 + om_v0 - 1.0)/(1 - om_m0 - om_v0 + om_v0*Qa2 + om_m0/aa)
    om_m = omega_t*om_m0/(om_m0 + om_v0*aa*Qa2)
    om_v = omega_t*om_v0*Qa2/(om_v0*Qa2 + om_m0/aa)

    end subroutine omegas_hf

    !!JD end generalize to variable w

    logical function HM_par_inner()
    ! The loops over redshift are the outermost parallel regions, and give much better scaling than
    ! the small regions used to set up each redshift. Those inner regions are therefore only used
    ! when the redshift loop is not itself running in parallel (e.g. for a single redshift).
    !$ use omp_lib, only: omp_in_parallel

    HM_par_inner = .true.
    !$ HM_par_inner = .not. omp_in_parallel()

    end function HM_par_inner

    !!AM Below is for HMcode
    subroutine HMcode(this, State, CAMB_Pk)
    !!AM - A CAMB derived type that I need
    class(THalofit) :: this
    class(CAMBdata) :: State
    type(MatterPowerData) :: CAMB_Pk
    integer :: j, nz, iz_wiggle, npass, imead_base
    integer :: imead_pass(3)
    real(dl), allocatable :: kh(:)
    type(HM_cosmology) :: cosi, cosm
    type(HM_tables) :: lut

    ! HMcode developed by Alexander Mead (alexander.j.mead@googlemail.com)
    ! Please contact me if you have any questions whatsoever
    ! If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
    ! If you use the extensions (w(a) and massive neutrinos) then please cite: http://arxiv.org/abs/1602.02154
    ! Also consider citing the source code at ASCL: http://ascl.net/1508.001

    ! Use imead to switch between the standard and accurate halo-model calculation
    ! 0 - Standard (this is just a vanilla halo model calculation with no accuracy tweaks)
    ! 1 - Accurate from Mead et al. (2016; arXiv 1602.02154)
    ! 2 - Accurate from Mead et al. (2015; arXiv 1505.07833)
    ! 3 - Accurate from Mead et al. (2020; arXiv 2009.01858)
    ! 4 - Denominator for feedback reaction model from Mead et al. (2020; arXiv 2009.01858)
    ! 5 - Numerator for feedback reaction from Mead et al. (2020; arXiv 2009.01858)
    imead_base = -1
    if (this%halofit_version == halofit_halomodel) imead_base = 0
    if (this%halofit_version == halofit_mead2016) imead_base = 1
    if (this%halofit_version == halofit_mead2015) imead_base = 2
    if (this%halofit_version == halofit_mead2020) imead_base = 3

    ! The 2020 feedback model is a response: HMcode 2020, then the denominator and numerator
    if (this%halofit_version == halofit_mead2020_feedback) then
        npass = 3
        imead_pass = [3, 4, 5]
    else
        npass = 1
        imead_pass(1) = imead_base
    end if

    HM_verbose = (FeedbackLevel > 1)

    if (HM_verbose) write(*, *)
    if (HM_verbose) write(*, *) 'HMcode: Running HMcode'
    if (HM_verbose) write(*, *)

    !!AM - Translate from CAMB variables to my variables
    nz = CAMB_Pk%num_z

    !!AM - Assign cosmological parameters for the halo model calculation
    call assign_HM_cosmology(this, State, cosi)

    ! Requested wavenumbers are shared by every redshift and feedback pass.
    kh = exp(CAMB_Pk%log_kh)

    ! Fill growth function table (only needs to be done once)
    call fill_growtab(cosi)

    if (cosi%wiggle_once) then
        ! Extract the BAO wiggle from the available spectrum nearest z=0. The
        ! wiggle is stored at z=0 and scaled to each redshift in p_dewiggle.
        iz_wiggle = minloc(abs(CAMB_Pk%redshifts), dim=1)
        call fill_plintab(iz_wiggle, cosi, CAMB_Pk)
        if (global_error_flag /= 0) return
        call init_wiggle(cosi)
    end if

    ! Loop over redshifts. The redshifts are independent, and give much better scaling than the
    ! small parallel regions used to set each one up, so this is the outermost parallel region:
    ! each thread takes its own copy of the cosmology (intrinsic assignment deep-copies the
    ! interpolation tables) and its own look-up tables. The regions inside then run serially
    ! (see HM_par_inner). For a single redshift it is faster not to open the region at all.
    if (nz > 1 .and. .not. HM_verbose) then
        !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(cosm, lut, j)
        cosm = cosi
        !$OMP DO SCHEDULE(DYNAMIC)
        do j = 1, nz
            call HMcode_redshift(this, CAMB_Pk, j, npass, imead_pass, kh, cosm, lut)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    else
        cosm = cosi
        do j = 1, nz
            call HMcode_redshift(this, CAMB_Pk, j, npass, imead_pass, kh, cosm, lut)
            ! Only report the halo-model set-up for the first redshift
            HM_verbose = .false.
        end do
    end if

    end subroutine HMcode

    subroutine HMcode_redshift(this, CAMB_Pk, j, npass, imead_pass, kh, cosm, lut)
    ! Everything HMcode does for one redshift; cosm and lut are the caller's working space, so
    ! this can be called either serially or from a thread with its own private copies of them
    class(THalofit) :: this
    type(MatterPowerData) :: CAMB_Pk
    integer, intent(in) :: j, npass, imead_pass(:)
    real(dl), intent(in) :: kh(:)
    type(HM_cosmology) :: cosm
    type(HM_tables) :: lut
    real(dl), allocatable :: p_den(:), plin_k(:)
    real(dl) :: z, k, p1h, p2h, pfull, plin
    integer :: i, ii, imead

    if (global_error_flag /= 0) return

    ! Initialise the specific HM_cosmology (fill sigma(R) and P_lin HM_tables)
    ! Currently this needs to be done at each z (mainly because of scale-dependent growth with neutrinos)
    ! For non-massive-neutrino models this could only be done once, which would speed things up a bit
    call initialise_HM_cosmology(this, j, cosm, CAMB_Pk)
    if (global_error_flag /= 0) return

    ! Sets the current redshift from the table
    z = CAMB_Pk%redshifts(j)

    if (npass > 1) then
        allocate(p_den(CAMB_Pk%num_k), plin_k(CAMB_Pk%num_k))
        !$OMP PARALLEL DO IF(HM_par_inner()), DEFAULT(SHARED)
        do i = 1, CAMB_Pk%num_k
            plin_k(i) = p_lin(kh(i), z, 0, cosm)
        end do
        !$OMP END PARALLEL DO
    end if

    do ii = 1, npass
        imead = imead_pass(ii)

        ! Initialisation for the halomodel calculation (needs to be done for each z).
        ! imead = 3,4,5 share everything except the halo concentration amplitude, so only the
        ! first pass of the feedback response model does the expensive set-up.
        call this%halomod_init(imead, z, lut, cosm, reuse=ii > 1)
        if (global_error_flag /= 0) return

        ! Loop over k values and calculate P(k)
        !$OMP PARALLEL DO IF(HM_par_inner()), DEFAULT(SHARED), PRIVATE(k, plin, pfull, p1h, p2h)
        do i = 1, CAMB_Pk%num_k
            k = kh(i)
            if (npass > 1) then
                plin = plin_k(i)
            else
                plin = p_lin(k, z, 0, cosm)
            end if
            call this%halomod(k, p1h, p2h, pfull, plin, lut, cosm)
            if (npass == 1 .or. imead == 3) then
                CAMB_Pk%nonlin_ratio(i, j) = sqrt(pfull/plin)
            else if (imead == 4) then
                p_den(i) = pfull
            else
                ! Non-linear correction from the feedback response of HMcode 2020
                CAMB_Pk%nonlin_ratio(i, j) = CAMB_Pk%nonlin_ratio(i, j)*sqrt(pfull/p_den(i))
            end if
        end do
        !$OMP END PARALLEL DO
    end do

    end subroutine HMcode_redshift

    function Delta_v(z, lut, cosm)
    ! Function for the virialised overdensity
    real(dl) :: Delta_v
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut

    if (lut%imead == 1 .or. lut%imead == 2) then
        ! Mead et al. (2015; arXiv 1505.07833) value
        Delta_v = 418*Omega_m_hm(z, cosm)**(-0.352_dl)
        ! Mead et al. (2016; arXiv 1602.02154) neutrino addition
        if (lut%imead == 1) Delta_v = Delta_v*(1 + 0.916_dl*cosm%f_nu)
    else if (lut%imead == 0 .or. lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
        Delta_v = Dv_Mead(z, cosm)
    end if

    end function Delta_v

    function delta_c(z, lut, cosm)
    ! Function for the linear collapse density
    real(dl) :: delta_c
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut

    if (lut%imead == 1 .or. lut%imead == 2) then
        ! Mead et al. (2015; arXiv 1505.07833) value
        delta_c = 1.59 + 0.0314*log(lut%sig8z)
        if (lut%imead == 1) then
            delta_c = delta_c*(1. + 0.262*cosm%f_nu) ! Mead et al. (2016; arXiv 1602.02154) neutrino addition
            delta_c = delta_c*(1. + 0.0123*log10(Omega_m_hm(z, cosm))) ! Nakamura & Suto (1997) fitting formula for LCDM
        end if
    else if (lut%imead == 0 .or. lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
        delta_c = dc_Mead(z, cosm)
    end if

    end function delta_c

    function eta(lut, cosm)
    ! Function eta that puffs halo profiles
    real(dl) :: eta
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut
    real(dl) :: eta0

    if (lut%imead == 0 .or. lut%imead == 4 .or. lut%imead == 5) then
        eta = 0.
    else if (lut%imead == 1 .or. lut%imead == 2) then
        ! The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
        ! eta = 0.603-0.3*lut%sig8z
        ! AM - made baryon feedback parameter obvious
        eta0 = cosm%eta_baryon
        ! eta0 = 1.03-0.11*cosm%A_baryon !Original one-parameter relation from 1505.07833
        ! eta0 = 0.98-0.12*cosm%A_baryon !Updated one-parameter relation: Section 4.1.2 of 1707.06627
        eta = eta0 - 0.3*lut%sig8z
    else if (lut%imead == 3) then
        eta = 0.1281*lut%sig8z_cold**(-0.3644)
    end if

    end function eta

    function kstar(lut)
    ! Function k* that cuts off the 1-halo term at large scales
    real(dl) :: kstar
    type(HM_tables), intent(in) :: lut

    if (lut%imead == 0) then
        ! Set to zero for the standard Poisson one-halo term
        kstar = 0.
    else if (lut%imead == 1 .or. lut%imead == 2) then
        ! One-halo cut-off wavenumber
        ! Mead et al. (2015; arXiv 1505.07833) value
        kstar = 0.584*(lut%sigv)**(-1.)
    else if (lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
        kstar = 0.05618*lut%sig8z_cold**(-1.013)
    end if

    end function kstar

    function As(lut, cosm)
    ! Halo concentration pre-factor from Bullock et al. (2001) relation
    type(HM_tables), intent(in) :: lut
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: As
    real(dl) :: B0, Bz, theta

    if (lut%imead == 0 .or. lut%imead == 4) then
        ! Set to 4 for the standard Bullock value
        As = 4.
    else if (lut%imead == 5) then
        theta = cosm%logT_AGN - 7.8
        B0 = 3.44 - 0.496*theta
        Bz = -0.0671 - 0.0371*theta
        As = B0*10**(lut%z*Bz)
    else if (lut%imead == 1 .or. lut%imead == 2) then
        ! This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
        ! As = 3.13
        ! AM - added for easy modification of feedback parameter
        As = cosm%A_baryon
    else if (lut%imead == 3) then
        As = 5.196
    end if

    end function As

    function fdamp(lut)
    ! Linear power damping function from Mead et al. (2015; arXiv 1505.07833)
    real(dl) :: fdamp
    type(HM_tables), intent(in) :: lut

    ! Linear theory damping factor
    if (lut%imead == 0 .or. lut%imead == 4 .or. lut%imead == 5) then
        ! Set to 0 for the standard linear theory two halo term
        fdamp = 0.
    else
        if (lut%imead == 1) then
            ! Mead et al. (2016; arXiv 1602.02154) value
            fdamp = 0.0095*lut%sigv100**1.37
        else if (lut%imead == 2) then
            ! Mead et al. (2015) value
            fdamp = 0.188*lut%sig8z**4.29
        else if (lut%imead == 3) then
            fdamp = 0.2696*lut%sig8z_cold**0.9403
        end if

        ! Catches extreme values of fdamp (only for the models that actually damp)
        if (fdamp < fdamp_min) fdamp = fdamp_min
        if (fdamp > fdamp_max) fdamp = fdamp_max
    end if

    end function fdamp

    function alpha(lut)
    ! Two- to one-halo transition smoothing from Mead et al. (2015; arXiv 1505.07833)
    real(dl) :: alpha
    type(HM_tables), intent(in) :: lut

    if (lut%imead == 0 .or. lut%imead == 4 .or. lut%imead == 5) then
        ! Set to 1 for the standard halomodel sum of one- and two-halo terms
        alpha = 1.
    else if (lut%imead == 1) then
        ! This uses the top-hat defined neff (HALOFIT uses Gaussian filtered fields instead)
        ! Mead et al. (2016; arXiv 1602.02154) value
        alpha = 3.24*1.85**lut%neff
    else if (lut%imead == 2) then
        ! Mead et al. (2015) value
        alpha = 2.93*1.77**lut%neff
    else if (lut%imead == 3) then
        alpha = 1.875*(1.603)**lut%neff
    end if

    ! Catches values of alpha that are crazy
    if (alpha > alpha_max) alpha = alpha_max
    if (alpha < alpha_min) alpha = alpha_min

    end function alpha

    function r_nl(lut)
    ! Calculates R_nl, defined by nu(R_nl)=1., nu=dc/sigma(R)
    type(HM_tables), intent(in) :: lut
    real(dl) :: r_nl

    if (lut%nu(1) > 1.) then
        ! This catches some very strange values that appear for odd cosmological models
        ! This is a terrible fudge, but I cannot think of a better solution
        r_nl = lut%rr(1)
    else
        r_nl = exp(interp_cubic_sorted(log(1.d0), log(lut%nu), log(lut%rr), lut%n))
    end if

    end function r_nl

    subroutine halomod(this, k, p1h, p2h, pfull, plin, lut, cosm)
    class(THalofit) :: this
    ! Calculates 1-halo and 2-halo terms and combines them to form the full halomodel power
    real(dl), intent(out) :: p1h, p2h, pfull
    real(dl), intent(in) :: plin, k
    real(dl) :: a
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut

    ! Calls expressions for one- and two-halo terms and then combines
    ! to form the full power spectrum
    if (k == 0.) then
        p1h = 0.
        p2h = 0.
    else
        p1h = this%p_1h(k, lut, cosm)
        p2h = this%p_2h(k, plin, lut, cosm)
    end if

    if (lut%imead == 1 .or. lut%imead == 2 .or. lut%imead == 3) then
        a = lut%alpha_hm
        pfull = (p2h**a + p1h**a)**(1./a)
    else
        pfull = p2h + p1h
    end if

    end subroutine halomod

    subroutine ensure_table_size(arr, n)
    ! Allocate a work table once and only resize it if its required extent changes.
    real(dl), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: n

    if (allocated(arr)) then
        if (size(arr) == n) return
        deallocate(arr)
    end if
    allocate(arr(n))

    end subroutine ensure_table_size

    subroutine fill_table(min, max, arr, n)
    ! Fills array 'arr' in equally spaced intervals
    integer :: i
    real(dl), intent(in) :: min, max
    real(dl), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: n

    call ensure_table_size(arr, n)

    if (n == 1) then
        arr(1) = min
    else if (n > 1) then
        do i = 1, n
            arr(i) = min + (max - min)*real(i - 1, dl)/(n - 1)
        end do
    end if

    end subroutine fill_table

    subroutine fill_plintab(iz, cosm, CAMB_PK)
    ! Fills internal HMcode HM_tables for the linear power spectrum at z=0
    type(MatterPowerData), intent(in) :: CAMB_PK
    integer, intent(in) :: iz
    type(HM_cosmology) :: cosm
    integer :: i
    real(dl) :: z, g, g2, k, Pk, Pkc
    real(dl), parameter :: pi = pi_HM
    real(dl), parameter :: kmin = kmin_pk_interpolation
    real(dl), parameter :: kmax = kmax_pk_interpolation
    integer :: nk, index_cache

    ! The table for this redshift is already there (e.g. filled to extract the BAO wiggle)
    if (cosm%plin_iz == iz) return

    nk = nk_pk_interpolation

    if (HM_verbose) write(*, *) 'LINEAR POWER: Filling linear power HM_tables'

    ! Fill a k-table with an equal-log-spaced k range (find_pk assumes this spacing)
    ! Note that the minimum should be such that the linear spectrum is accurately a power-law below this wavenumber
    cosm%nk = nk
    call fill_table(log(kmin), log(kmax), cosm%log_k_plin, nk)
    call ensure_table_size(cosm%log_plin, nk)
    call ensure_table_size(cosm%log_plinc, nk)
    cosm%kmax = exp(cosm%log_k_plin(nk))

    if (HM_verbose) write(*, *) 'LINEAR POWER: k_min:', exp(cosm%log_k_plin(1))
    if (HM_verbose) write(*, *) 'LINEAR POWER: k_max:', cosm%kmax
    if (HM_verbose) write(*, *) 'LINEAR POWER: nk:', nk

    ! Find the redshift
    z = CAMB_PK%redshifts(iz)
    if (HM_verbose) write(*, *) 'LINEAR POWER: z of input:', z
    g = grow(z, cosm)
    g2 = g**2
    cosm%grow_z = g
    cosm%this_z = z
    index_cache = 1
    ! Fill power table, both cold- and all-matter
    !$OMP PARALLEL DO IF(HM_par_inner() .and. nk>256), DEFAULT(SHARED), PRIVATE(k, Pk, Pkc), &
    !$OMP FIRSTPRIVATE(index_cache), SCHEDULE(STATIC)
    do i = 1, nk
        ! Take the power from the current redshift choice
        k = exp(cosm%log_k_plin(i))
        Pk = MatterPowerData_k(CAMB_PK, k, iz, index_cache)*(k**3/(2*pi**2))
        ! Pk is zero if MatterPowerData_k flagged an error; skip the log and let the check below fire
        if (Pk > 0) then
            Pkc = Pk*Tcb_Tcbnu_ratio(k, z, cosm)**2
            cosm%log_plin(i) = log(Pk/g2)
            cosm%log_plinc(i) = log(Pkc/g2)
        end if
    end do
    if (global_error_flag /= 0) return

    if (HM_verbose) write(*, *) 'LINEAR POWER: Delta2_min:', exp(cosm%log_plin(1))*g2
    if (HM_verbose) write(*, *) 'LINEAR POWER: Delta2_max:', exp(cosm%log_plin(nk))*g2
    cosm%plin_iz = iz

    ! Check sigma_8 value
    if (HM_verbose) write(*, *) 'LINEAR POWER: sigma_8:', sigma_integral(8.d0, 0.d0, 0, cosm)
    if (HM_verbose) write(*, *) 'LINEAR POWER: Done'
    if (HM_verbose) write(*, *)

    end subroutine fill_plintab

    function Tcb_Tcbnu_ratio(k, z, cosm)
    ! Calculates the ratio of T(k) for cold vs. all matter
    ! Uses approximations in Eisenstein & Hu (1999; arXiv 9710252)
    ! Note that this assumes that there are exactly 3 species of neutrinos with
    ! Nnu<=3 of these being massive, and with the mass split evenly between the number of massive species.

    real(dl) :: Tcb_Tcbnu_ratio
    real(dl), intent(in) :: k, z
    real(dl) :: D, Dcb, Dcbnu, pcb, zeq, q, yfs
    real(dl) :: BigT
    type(HM_cosmology), intent(in) :: cosm

    if (cosm%f_nu == 0._dl) then

        Tcb_Tcbnu_ratio = 1.

    else

        ! Growth exponent under the assumption that neutrinos are completely unclustered (equation 11)
        pcb = (5. - sqrt(1 + 24*(1._dl - cosm%f_nu)))/4

        ! Theta for temperature (BigT=T/2.7 K)
        BigT = cosm%Tcmb/2.7

        ! The matter-radiation equality redshift
        zeq = (2.5e4)*cosm%om_m*(cosm%h**2)*(BigT**(-4))

        ! The growth function normalised such that D=(1.+z_eq)/(1+z) at early times (when Omega_m \approx 1)
        if (EdS_Tcold_growth) then
            D = (1. + zeq)/(1. + z) ! EdS solution
        else
            D = (1. + zeq)*ungrow(z, cosm) ! General solution
        end if

        ! Wave number relative to the horizon scale at equality (equation 5)
        ! Extra factor of h because all my k are in units of h/Mpc
        q = k*cosm%h*BigT**2/(cosm%om_m*cosm%h**2)

        ! Free streaming scale (equation 14)
        ! Note that Eisenstein & Hu (1999) only consider the case of 3 neutrinos
        ! with Nnu of these being massve with the mass split evenly between Nnu species.
        yfs = 17.2*cosm%f_nu*(1. + 0.488*cosm%f_nu**(-7./6.))*(cosm%Nnu*q/cosm%f_nu)**2

        ! These are (almost) the scale-dependent growth functions for each component in Eisenstein & Hu (1999)
        ! Some part is missing, but this cancels when they are divided by each other, which is all I need them for.
        ! Equations (12) and (13)
        Dcb = (1. + (D/(1. + yfs))**0.7)**(pcb/0.7)
        Dcbnu = ((1 - cosm%f_nu)**(0.7/pcb) + (D/(1. + yfs))**0.7)**(pcb/0.7)

        Tcb_Tcbnu_ratio = Dcb/Dcbnu

    end if

    end function Tcb_Tcbnu_ratio

    subroutine assign_HM_cosmology(this, State, cosm)
    class(THalofit) :: this
    class(CAMBdata) :: State
    ! Assigns the internal HMcode cosmological parameters
    type(HM_cosmology) :: cosm
    real(dl) h2
    logical use_wiggle

    associate(CP => State%CP)
        ! Converts CAMB parameters to Meadfit parameters
        h2 = (CP%H0/100)**2
        cosm%om_m = (CP%omch2 + CP%ombh2 + CP%omnuh2)/h2
        cosm%om_c = CP%omch2/h2
        cosm%om_b = CP%ombh2/h2
        cosm%om_nu = CP%omnuh2/h2
        cosm%om_v = State%Omega_de
        call CP%DarkEnergy%Effective_w_wa(cosm%w, cosm%wa)
        cosm%f_nu = cosm%om_nu/cosm%om_m
        cosm%h = CP%H0/100
        cosm%Tcmb = CP%TCMB
        cosm%Nnu = CP%Num_Nu_massive
        cosm%ns = CP%InitPower%Effective_ns()
        ! The de-wiggled spectrum is only used by the 2020 versions. The wiggle is stored at z=0 and
        ! scaled to other redshifts by the scale-independent growth factor, so it only needs extracting
        ! once when growth is nearly scale-independent (small neutrino fraction) and the dark energy
        ! model does not itself introduce scale-dependent growth. Otherwise re-extract it at each redshift.
        use_wiggle = this%halofit_version == halofit_mead2020 .or. this%halofit_version == halofit_mead2020_feedback
        cosm%wiggle_once = use_wiggle .and. cosm%f_nu < this%HMcode_wiggle_max_fnu .and. &
            CP%DarkEnergy%assume_scale_indep_lowz_growth()
        cosm%wiggle_each_z = use_wiggle .and. .not. cosm%wiggle_once
    end associate

    ! Baryon feedback parameters
    if (this%halofit_version == halofit_mead2015 .or. this%halofit_version == halofit_mead2016) then
        cosm%A_baryon = this%HMcode_A_baryon
        cosm%eta_baryon = this%HMcode_eta_baryon
    else if (this%halofit_version == halofit_mead2020_feedback) then
        cosm%logT_AGN = this%HMcode_logT_AGN
    end if

    ! Write out cosmological parameters if necessary
    if (HM_verbose) write(*, *) 'HM_cosmology: Om_m:', cosm%om_m
    if (HM_verbose) write(*, *) 'HM_cosmology: Om_c:', cosm%om_c
    if (HM_verbose) write(*, *) 'HM_cosmology: Om_b:', cosm%om_b
    if (HM_verbose) write(*, *) 'HM_cosmology: Om_nu:', cosm%om_nu
    if (HM_verbose) write(*, *) 'HM_cosmology: Om_v:', cosm%om_v
    if (HM_verbose) write(*, *) 'HM_cosmology: w_0:', cosm%w
    if (HM_verbose) write(*, *) 'HM_cosmology: w_a:', cosm%wa
    if (HM_verbose) write(*, *) 'HM_cosmology: f_nu:', cosm%f_nu
    if (HM_verbose) write(*, *) 'HM_cosmology: n_s:', cosm%ns
    if (HM_verbose) write(*, *) 'HM_cosmology: h:', cosm%h
    if (HM_verbose) write(*, *) 'HM_cosmology: T_CMB [K]:', cosm%Tcmb
    if (HM_verbose) write(*, *) 'HM_cosmology: N_nu (massive):', cosm%Nnu
    if (HM_verbose .and. (this%halofit_version == halofit_mead2015 .or. this%halofit_version == halofit_mead2016)) then
        write(*, *) 'HM_cosmology: A_baryon:', cosm%A_baryon
        write(*, *) 'HM_cosmology: eta_baryon:', cosm%eta_baryon
    else if (HM_verbose .and. this%halofit_version == halofit_mead2020_feedback) then
        write(*, *) 'HM_cosmology: log10(T_AGN/K):', cosm%logT_AGN
    end if
    if (HM_verbose) write(*, *)

    end subroutine assign_HM_cosmology

    subroutine initialise_HM_cosmology(this, iz, cosm, CAMB_PK)
    class(THalofit) :: this
    ! Sets up HM_tables of sigma, growth and linear power for the HM_cosmology
    type(MatterPowerData), intent(in) :: CAMB_PK
    type(HM_cosmology) :: cosm
    integer, intent(in) :: iz

    ! Fill linear power table and grows it to z=0
    call fill_plintab(iz, cosm, CAMB_PK)
    if (global_error_flag /= 0) return

    ! Fill sigma(r) table
    call fill_sigtab(this, cosm)

    ! Extract BAO wiggle from P(k); only needed at each z for models where the
    ! wiggle is not simply growth-scaled (otherwise done once in HMcode)
    if (cosm%wiggle_each_z) call init_wiggle(cosm)

    end subroutine initialise_HM_cosmology

    subroutine allocate_LUT(lut, n, new_mass_grid)
    ! Allocates memory for the HMcode look-up HM_tables. Nothing is initialised here: every
    ! array and cached scalar is fully overwritten before use, by halomod_tables/fill_conc (once
    ! per redshift) or by halomod_init immediately afterwards; the feedback-only fields
    ! (baryon_mass_fraction, f_star_hm) are only ever read by the imead variant that fills them.
    type(HM_tables) :: lut
    integer, intent(in) :: n
    logical, intent(out) :: new_mass_grid

    new_mass_grid = .false.
    if (allocated(lut%zc)) then
        if (size(lut%zc) /= n) then
            deallocate(lut%zc, lut%m, lut%c, lut%rv, lut%rs, lut%nfw_norm, lut%nfw_moment2)
            deallocate(lut%nu, lut%rr, lut%sigf, lut%sig)
            deallocate(lut%p1h_quad_weight, lut%nu_eta, lut%baryon_mass_fraction)
        end if
    end if
    if (.not. allocated(lut%zc)) then
        lut%n = n
        allocate(lut%zc(n), lut%m(n), lut%c(n), lut%rv(n), lut%rs(n), lut%nfw_norm(n), lut%nfw_moment2(n))
        allocate(lut%nu(n), lut%rr(n), lut%sigf(n), lut%sig(n))
        allocate(lut%p1h_quad_weight(n), lut%nu_eta(n), lut%baryon_mass_fraction(n))
        new_mass_grid = .true.
    end if

    end subroutine allocate_LUT

    subroutine halomod_init(this, imead, z, lut, cosm, reuse)
    class(THalofit) :: this
    ! Halo-model initialisation routine
    ! Computes look-up HM_tables necessary for the halo model calculations
    integer, intent(in) :: imead
    real(dl), intent(in) :: z
    ! reuse: the tables in lut were already filled at this z for a variant of HMcode that shares
    ! everything except the halo-model parameters set here (only used for the 2020 feedback response,
    ! where imead=3,4,5 give identical delta_c, Delta_v, n_eff, collapse redshifts and Dolag correction)
    logical, intent(in), optional :: reuse
    integer :: i
    real(dl) :: fb, fc, fs, mb, beta, ratio
    type(HM_cosmology) :: cosm
    type(HM_tables) :: lut

    lut%imead = imead
    lut%z = z

    if (.not. DefaultFalse(reuse)) then
        call halomod_tables(this, z, lut, cosm)
        if (global_error_flag /= 0) return
    end if

    ! Get the concentration for all the haloes (cheap, and the only part of the set-up that
    ! changes between the passes of the 2020 feedback response model)
    call fill_conc(this, z, lut, cosm)

    ! Cache redshift-only halo-model factors used in every k evaluation.
    lut%eta_hm = this%eta(lut, cosm)
    lut%kstar_hm = this%kstar(lut)
    lut%fdamp_hm = this%fdamp(lut)
    lut%alpha_hm = this%alpha(lut)
    if (lut%imead == 3) lut%kdamp_hm = 0.05699_dl*lut%sig8z_cold**(-1.089_dl)
    lut%one_minus_fnu_sq = (1._dl - cosm%f_nu)**2
    lut%nu_eta = lut%nu**lut%eta_hm
    if (lut%imead == 5) then
        mb = m_baryon(lut, cosm)
        beta = 2._dl
        fb = cosm%om_b/cosm%om_m
        fc = cosm%om_c/cosm%om_m
        fs = f_star(lut, cosm)
        lut%f_star_hm = fs
        do i = 1, lut%n
            ratio = (lut%m(i)/mb)**beta
            lut%baryon_mass_fraction(i) = fc + (fb - fs)*ratio/(1._dl + ratio)
        end do
    end if

    if (HM_verbose) write(*, *) 'HALOMOD: c HM_tables filled'
    if (HM_verbose) write(*, *) 'HALOMOD: c min [Msun/h]:', lut%c(lut%n)
    if (HM_verbose) write(*, *) 'HALOMOD: c max [Msun/h]:', lut%c(1)
    if (HM_verbose) write(*, *) 'HALOMOD: Done'
    if (HM_verbose) write(*, *)
    if (HM_verbose) call this%write_parameters(z, lut, cosm)

    end subroutine halomod_init

    subroutine halomod_tables(this, z, lut, cosm)
    class(THalofit) :: this
    ! Fills the look-up tables for the halo model at redshift z that do not depend on the
    ! halo concentration amplitude (see halomod_init)
    real(dl), intent(in) :: z
    type(HM_cosmology) :: cosm
    type(HM_tables) :: lut
    integer :: i, nm
    real(dl) :: Dv, dc, m, nu, r, sig, mmin, mmax
    real(dl), parameter :: f_Bullock = 0.01_dl**(1/3._dl)
    logical :: new_mass_grid

    if (HM_verbose) write(*, *) 'HALOMOD: Filling look-up HM_tables'
    if (HM_verbose) write(*, *) 'HALOMOD: HM_tables being filled at redshift:', z

    ! Mass range and number of points
    mmin = 1e0
    mmax = 1e18
    nm = 256

    ! Find value of sigma_v, sig8, etc.
    lut%sigv = sigmaV(0.d0, z, 0, cosm)
    if (HM_verbose) write(*, *) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
    if (global_error_flag == 0) lut%sigv100 = sigmaV(100.d0, z, 0, cosm)
    if (HM_verbose) write(*, *) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
    if (global_error_flag == 0) then
        lut%sig8z = sigma_integral(8.d0, z, 0, cosm)
        lut%sig8z_cold = sigma_integral(8.d0, z, 1, cosm)
    end if
    if (HM_verbose) then
        write(*, *) 'HALOMOD: sig8(z):', lut%sig8z
        write(*, *) 'HALOMOD: cold sig8(z):', lut%sig8z_cold
    end if
    if (global_error_flag /= 0) return

    call allocate_LUT(lut, nm, new_mass_grid)

    if (HM_verbose) write(*, *) 'HALOMOD: M_min [log10(Msun/h)]:', log10(mmin)
    if (HM_verbose) write(*, *) 'HALOMOD: M_max [log10(Msun/h)]:', log10(mmax)

    dc = this%delta_c(z, lut, cosm)
    lut%dc = dc

    !$OMP PARALLEL DO IF(HM_par_inner()), DEFAULT(SHARED), PRIVATE(m, r, sig, nu)
    do i = 1, lut%n

        if (new_mass_grid) then
            lut%m(i) = exp(log(mmin) + log(mmax/mmin)*real(i - 1, dl)/(lut%n - 1))
            lut%rr(i) = radius_m(lut%m(i), cosm)
        end if
        m = lut%m(i)
        r = lut%rr(i)
        sig = sigma_lut(r, z, cosm)
        nu = dc/sig

        lut%sig(i) = sig
        lut%nu(i) = nu
        lut%sigf(i) = sigma_lut(r*f_Bullock, z, cosm)

    end do
    !$OMP END PARALLEL DO
    if (global_error_flag /= 0) return

    ! Fold the trapezium-rule intervals into per-node weights once per redshift.
    lut%p1h_quad_weight(1) = 0.5_dl*gst(lut%nu(1))*lut%m(1)*(lut%nu(2) - lut%nu(1))
    do i = 2, lut%n - 1
        lut%p1h_quad_weight(i) = 0.5_dl*gst(lut%nu(i))*lut%m(i)*(lut%nu(i + 1) - lut%nu(i - 1))
    end do
    lut%p1h_quad_weight(lut%n) = &
        0.5_dl*gst(lut%nu(lut%n))*lut%m(lut%n)*(lut%nu(lut%n) - lut%nu(lut%n - 1))

    if (HM_verbose) write(*, *) 'HALOMOD: m, r, nu, sig, sigf HM_tables filled'

    ! Fill virial radius table using real radius table
    Dv = this%Delta_v(z, lut, cosm)
    lut%rv = lut%rr/(Dv**(1/3._dl))

    if (HM_verbose) write(*, *) 'HALOMOD: rv HM_tables filled'
    if (HM_verbose) write(*, *) 'HALOMOD: nu min:', lut%nu(1)
    if (HM_verbose) write(*, *) 'HALOMOD: nu max:', lut%nu(lut%n)
    if (HM_verbose) write(*, *) 'HALOMOD: sig min:', lut%sig(lut%n)
    if (HM_verbose) write(*, *) 'HALOMOD: sig max:', lut%sig(1)
    if (HM_verbose) write(*, *) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    if (HM_verbose) write(*, *) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)

    ! Find non-linear radius and scale
    lut%rnl = r_nl(lut)
    lut%knl = 1/lut%rnl

    if (HM_verbose) write(*, *) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    if (HM_verbose) write(*, *) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    ! Calculate the effective spectral index at the collapse scale
    lut%neff = neff(lut, cosm)

    if (HM_verbose) write(*, *) 'HALOMOD: n_eff:', lut%neff

    ! Get the halo collapse redshifts and the Dolag (2004) dark-energy correction
    call this%conc_bull(z, lut, cosm)

    end subroutine halomod_tables

    subroutine fill_conc(this, z, lut, cosm)
    class(THalofit) :: this
    ! Fills the halo concentration table from the collapse redshifts and the (already cached)
    ! Dolag (2004) dark-energy corrections. This is the only z set-up that depends on the
    ! halo-concentration amplitude A, which is what differs between the feedback response passes.
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables) :: lut
    real(dl) :: A, c
    integer :: i

    ! Amplitude of relation (4. in Bullock et al. 2001), including the Dolag corrections
    A = this%As(lut, cosm)*lut%dolag_inf*lut%dolag_z

    do i = 1, lut%n
        c = A*(1. + lut%zc(i))/(1. + z)
        lut%c(i) = c
        ! Scale radius and profile normalisation, cached since win() is called at every k
        lut%rs(i) = lut%rv(i)/c
        lut%nfw_norm(i) = 1/(log(1 + c) - c/(1 + c))
        lut%nfw_moment2(i) = 0.5_dl*c*c - 2._dl*c + 3._dl*log(1._dl + c) - c/(1._dl + c)
    end do

    end subroutine fill_conc

    subroutine write_parameters(this, z, lut, cosm)
    class(THalofit) :: this
    ! This subroutine writes out the halomodel parameters at the current redshift
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut

    if (HM_verbose) write(*, *) 'WRITE_PARAMETERS: at this redshift'
    if (HM_verbose) write(*, *) '=================================='
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'z:', z
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'Dv:', this%Delta_v(z, lut, cosm)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'dc:', this%delta_c(z, lut, cosm)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'eta:', this%eta(lut, cosm)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'k*:', this%kstar(lut)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'A:', this%As(lut, cosm)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'fdamp:', this%fdamp(lut)
    if (HM_verbose) write(*, fmt='(A10,F10.5)') 'alpha:', this%alpha(lut)
    if (HM_verbose) write(*, *) '=================================='
    if (HM_verbose) write(*, *)

    end subroutine write_parameters

    pure function radius_m(m, cosm)
    ! Calculates the co-moving radius that encloses a mass 'm' in the homogeneous Universe
    real(dl) :: radius_m
    real(dl), intent(in) :: m
    type(HM_cosmology), intent(in) :: cosm
    real(dl), parameter :: pi = pi_HM

    radius_m = (3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

    end function radius_m

    function neff(lut, cosm)
    ! Finds the effective spectral index at the collapse scale r_nl, where nu(r_nl)=1.
    real(dl) :: neff
    real(dl) :: ns
    type(HM_cosmology), intent(in) :: cosm
    type(HM_tables), intent(in) :: lut
    real(dl) :: sig
    integer :: itype
    real(dl), parameter :: tmin = 0.
    real(dl), parameter :: tmax = 1.
    real(dl), parameter :: acc = acc_neff_integration

    ! Choose type of sigma(R) to tabulate depending on HMcode version
    if (lut%imead == 1 .or. lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
        itype = 1 ! 1 - Cold matter
    else
        itype = 0 ! 0 - All matter
    end if

    ! Choosings sig = delta_c should be equivalent to actually calculating it again, however
    ! Do the actual calculation to be consistent with HMx. Problems with weird cosmologies with
    ! low spectral indices such that no collapse has occurred. R_nl very small
    ! sig = lut%dc ! Take great care here. This should be the same as below, but won't be for strange models
    sig = sigma_integral(lut%rnl, lut%z, itype, cosm)
    neff = -3. - 2.*integrate(tmin, tmax, neff_integrand, lut%rnl, lut%z, itype, cosm, acc)/sig**2

    ! For some bizarre cosmological models r_nl is very small, so almost no collapse has occurred
    ! In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    ns = cosm%ns
    if (neff < ns - 4.) neff = ns - 4.
    if (neff > ns) neff = ns

    end function neff

    subroutine conc_bull(this, z, lut, cosm)
    class(THalofit) :: this
    ! Sets up the Bullock et al. (2001) concentration-mass relation: fills the collapse redshift
    ! table and caches the Dolag (2004) dark-energy corrections. The concentrations themselves are
    ! then filled by fill_conc, which is all that has to be redone if only the amplitude A changes.
    real(dl), intent(in) :: z
    type(HM_cosmology) :: cosm, cosm_lcdm
    type(HM_tables) :: lut
    real(dl) :: pow
    real(dl) :: ginf_lcdm, ginf_wcdm, g_lcdm, g_wcdm
    real(dl), parameter :: zinf = zc_Dolag

    ! Fill the collapse time look-up table
    call this%zcoll_bull(z, cosm, lut)

    lut%dolag_inf = 1
    lut%dolag_z = 1

    if (z < zinf) then

        ! Dolag (2004) prescription for adding DE dependence

        ! This is approximately z=infinity
        ginf_wcdm = grow(zinf, cosm)

        ! Make a LCDM HM_cosmology
        ! Need to make sure model is flat with the same Omega_m and w=-1
        ! This is *only* used for a calculation of the growth function, which depends on the
        ! scalars below alone, so avoid copying the (large) interpolation tables in cosm
        cosm_lcdm%om_m = cosm%om_m
        cosm_lcdm%om_c = cosm%om_c
        cosm_lcdm%om_b = cosm%om_b
        cosm_lcdm%w = -1.
        cosm_lcdm%wa = 0.
        cosm_lcdm%om_v = 1. - cosm%om_m ! Enforce flatness

        ! Needs to use grow_int explicitly here for LCDM model to avoid growth HM_tables
        ginf_lcdm = growint(zinf, cosm_lcdm)

        ! This is the Dolag et al. (2004) correction for halo concentrations
        if (lut%imead == 0 .or. lut%imead == 2 .or. lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
            ! Mead et al. (2015) used the Dolag (2004) correction
            pow = 1.
        else if (lut%imead == 1) then
            ! Mead et al. (2016) changed the power to 1.5 to better accommodate more extreme dark-energy models
            pow = 1.5
        end if
        lut%dolag_inf = (ginf_wcdm/ginf_lcdm)**pow

        ! This is needed for the correction to make sense at high z
        if (lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
            g_lcdm = growint(z, cosm_lcdm)
            g_wcdm = grow(z, cosm)
            lut%dolag_z = (g_lcdm/g_wcdm)**pow
        end if

    end if

    end subroutine conc_bull

    function growint(z, cosm)
    ! Approximate growth function from integrating Omega_m(a)^gamma dln(a) up to a=1
    real(dl) :: growint
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: a
    real(dl), parameter :: b = 1._dl ! Integration range for integration parameter; note a -> 1
    real(dl), parameter :: acc = acc_growth_integration

    a = 1./(1. + z)
    growint = exp(integrate(a, b, growint_integrand, 0._dl, 0._dl, 0, cosm, acc))

    end function growint

    function growint_integrand(a, y, z, itype, cosm)
    ! Integrand for the approximate growth integral
    ! y, z and itype are unused, but are needed to match the integrate() interface
    real(dl) :: growint_integrand
    real(dl), intent(in) :: a
    real(dl), intent(in) :: y
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: Om_m, zz, gam

    if (cosm%w < -1.) then
        gam = 0.55 + 0.02*(1. + cosm%w)
    else if (cosm%w > -1) then
        gam = 0.55 + 0.05*(1. + cosm%w)
    else
        gam = 0.55
    end if

    ! Note the minus sign here
    ! AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    zz = -1. + 1./a
    if (cold_growth) then
        Om_m = Omega_cold_hm(zz, cosm)
    else
        Om_m = Omega_m_hm(zz, cosm)
    end if
    growint_integrand = -(Om_m**gam)/a

    end function growint_integrand

    subroutine zcoll_bull(this, z, cosm, lut)
    class(THalofit) :: this
    ! Calculates the halo collapse redshift according to the Bullock (2001) prescription
    real(dl), intent(in) :: z
    type(HM_cosmology) :: cosm
    type(HM_tables) :: lut
    real(dl) :: dc
    real(dl) :: af, zf, RHS, growz
    integer :: i

    ! This fills up the halo formation redshift table as per Bullock relations

    ! Needs to interpolate g(z) which should be pretty linear for a<0.05
    ! in 'g(a) vs a' space for all standard cosmologies

    dc = this%delta_c(z, lut, cosm)

    ! Find the growth function at the current redshift
    growz = grow(z, cosm)

    ! Do numerical inversion
    do i = 1, lut%n

        RHS = dc*growz/lut%sigf(i)

        if (RHS > growz) then
            ! This is the case of 'halo forms in the future'
            ! in this case set formation redshift to current redshift
            zf = z
        else
            ! Invert the growth table (increasing, but not equally spaced) for a
            af = interp_cubic_sorted(RHS, cosm%growth, cosm%a_growth, cosm%ng)
            zf = -1. + 1./af
        end if

        lut%zc(i) = zf

    end do

    end subroutine zcoll_bull

    pure function cosmic_density(cosm)
    ! The z=0 cosmological matter density
    real(dl) :: cosmic_density
    type(HM_cosmology), intent(in) :: cosm

    ! In M_sun per Mpc^3 with h factors included. The constant does this.
    cosmic_density = (2.775d11)*cosm%om_m

    end function cosmic_density

    function m_baryon(lut, cosm)

    real(dl) :: m_baryon
    type(HM_tables), intent(in) :: lut
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: mb0, mbz, theta

    theta = cosm%logT_AGN - 7.8
    mb0 = 13.87 + 1.81*theta
    mb0 = 10**mb0
    mbz = -0.108 + 0.195*theta
    m_baryon = mb0*10**(mbz*lut%z)

    end function m_baryon

    function f_star(lut, cosm)

    real(dl) :: f_star
    type(HM_tables), intent(in) :: lut
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: f0, fz, theta

    theta = cosm%logT_AGN - 7.8
    f0 = 0.0201 - 0.003*theta
    fz = 0.409 + 0.0224*theta
    f_star = f0*10**(fz*lut%z)

    if (f_star > cosm%om_b/cosm%om_m) f_star = cosm%om_b/cosm%om_m

    end function f_star

    function find_pk(k, itype, cosm)
    ! Look-up and interpolation for P(k,z=0); itype=0 for all matter, 1 for cold matter only
    real(dl) :: find_pk
    real(dl), intent(in) :: k
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    integer :: n
    real(dl) :: log_pk_max

    n = cosm%nk

    if (plin_extrap .and. k > cosm%kmax) then
        ! Do some extrapolation here based on knowledge of things at high k
        if (itype == 0) then
            log_pk_max = cosm%log_plin(n)
        else
            log_pk_max = cosm%log_plinc(n)
        end if
        find_pk = exp(log_pk_max)*((log(k)/cosm%log_k_plin(n))**2)*((k/cosm%kmax)**(cosm%ns - 1))
    else
        ! Otherwise interpolate the (equal-log-spaced) table
        if (itype == 0) then
            find_pk = exp(interp_linear_uniform(log(k), cosm%log_k_plin, cosm%log_plin, n))
        else
            find_pk = exp(interp_linear_uniform(log(k), cosm%log_k_plin, cosm%log_plinc, n))
        end if
    end if

    end function find_pk

    function cached_grow(z, cosm)
    ! The growth factor, using the value cached by fill_plintab where possible
    real(dl) :: cached_grow
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm

    if (z == 0._dl) then
        cached_grow = 1
    else if (z == cosm%this_z) then
        cached_grow = cosm%grow_z
    else
        cached_grow = grow(z, cosm) ! never actually needed
    end if

    end function cached_grow

    function p_lin(k, z, itype, cosm)
    ! Looks up the value for the linear power spectrum
    real(dl) :: p_lin
    real(dl), intent (in) :: k, z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm

    ! This gives the linear power spectrum for the model in question
    ! P(k) should have been previously normalised to z=0
    p_lin = cached_grow(z, cosm)**2*find_pk(k, itype, cosm)

    end function p_lin

    function p_2h(k, plin, lut, cosm)
    ! Calculates the 2-halo term
    real(dl) :: p_2h
    real(dl), intent(in) :: k, plin
    real(dl) :: sigv, frac, x
    real(dl), parameter :: ndamp = 2.85_dl
    type(HM_tables), intent(in) :: lut
    type(HM_cosmology), intent(in) :: cosm

    ! Damping function
    frac = lut%fdamp_hm

    ! frac<=fdamp_min means the damping has been clamped to its minimum, i.e. it is negligible
    if (lut%imead == 0 .or. lut%imead == 4 .or. lut%imead == 5 .or. frac <= fdamp_min) then
        p_2h = plin
    else if (lut%imead == 1 .or. lut%imead == 2) then
        sigv = lut%sigv
        p_2h = plin*(1. - frac*(tanh(k*sigv/sqrt(abs(frac))))**2)
    else if (lut%imead == 3) then
        x = (k/lut%kdamp_hm)**ndamp
        p_2h = p_dewiggle(k, lut%z, plin, lut%sigv, cosm)*(1. - frac*x/(1. + x))
    end if

    ! For some strange cosmologies frac>1. so this must be added to prevent p_2h<0.
    if (p_2h < 0.) p_2h = 0.

    end function p_2h

    function p_1h(k, lut, cosm)
    ! Calculates the 1-halo term
    real(dl) :: p_1h
    real(dl), intent(in) :: k
    type(HM_tables), intent(in) :: lut
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: fac, ks, x
    real(dl) :: sum, wk
    integer :: i
    real(dl), parameter :: pi = pi_HM

    ! Trapezium rule over nu; the redshift-only node weights were cached by halomod_tables.
    sum = 0
    do i = 1, lut%n
        wk = win(k*lut%nu_eta(i), lut%rs(i), lut%c(i), lut%nfw_norm(i), lut%nfw_moment2(i))
        if (lut%imead == 5) wk = wk*lut%baryon_mass_fraction(i) + lut%f_star_hm
        sum = sum + lut%p1h_quad_weight(i)*wk*wk
    end do
    sum = sum/cosmic_density(cosm)
    if (lut%imead == 3 .or. lut%imead == 4) sum = sum*lut%one_minus_fnu_sq

    ! Numerical factors to convert from P(k) to Delta^2(k)
    p_1h = sum*k**3/(2.*pi**2)

    if (lut%imead == 1 .or. lut%imead == 2) then
        ! Damping of the 1-halo term at very large scales
        ! Note kstar is only zero for imead==0, which does not reach here
        ks = lut%kstar_hm
        if ((k/ks)**2 > ks_limit) then
            fac = 0. ! Prevents problems if k/ks is very large
        else
            fac = exp(-((k/ks)**2))
        end if
        ! Damping of the one-halo term at very large scales
        p_1h = p_1h*(1. - fac)
    else if (lut%imead == 3 .or. lut%imead == 4 .or. lut%imead == 5) then
        ks = lut%kstar_hm
        x = (k/ks)**4
        p_1h = p_1h*x/(1. + x)
    end if

    end function p_1h

    real(dl) function p_dewiggle(k, z, p_linear, sigv, cosm)
    ! Call the dewiggled power spectrum, which is linear but with damped wiggles
    real(dl), intent(in) :: k, z, p_linear, sigv
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: p_wiggle, f, logk

    logk = log(k)
    if (logk < cosm%log_k_wiggle(1) .or. logk > cosm%log_k_wiggle(nk_wiggle)) then
        p_wiggle = 0
    else
        p_wiggle = interp_cubic_uniform(logk, cosm%log_k_wiggle, cosm%pk_wiggle, nk_wiggle)
    end if
    f = exp(-(k*sigv)**2)
    p_dewiggle = p_linear + (f - 1)*p_wiggle*cached_grow(z, cosm)**2

    end function p_dewiggle

    subroutine init_wiggle(cosm)
    ! Isolate the power spectrum wiggle
    type(HM_cosmology), intent(inout) :: cosm
    real(dl), allocatable :: k(:), Pk(:), Pk_smooth(:)
    integer :: i
    real(dl), parameter :: kmin = kmin_wiggle
    real(dl), parameter :: kmax = kmax_wiggle
    integer, parameter  :: nk = nk_wiggle
    real(dl), parameter :: z = 0. ! Only need to run this routine for z=0
    integer, parameter :: itype = 0 ! Matter spectrum

    ! Words
    if (HM_verbose) write(*, *) 'INIT_WIGGLE: Starting'

    ! Allocate arrays (Pk_smooth is allocated by calculate_psmooth)
    allocate(Pk(nk))

    ! Allocate array for k; fill log(k) directly (rather than as log of the k array) so that
    ! the look-up table is exactly equally spaced, as assumed by interp_cubic_uniform
    call fill_table(log(kmin), log(kmax), cosm%log_k_wiggle, nk)
    allocate(k(nk))
    k = exp(cosm%log_k_wiggle)

    ! Get the linear power spectrum in an array
    do i = 1, nk
        Pk(i) = p_lin(k(i), z, itype, cosm)
    end do

    if (HM_verbose) then
        write(*, *) 'INIT_WIGGLE: kmin [h/Mpc]:', k(1)
        write(*, *) 'INIT_WIGGLE: kmax [h/Mpc]:', k(nk)
        write(*, *) 'INIT_WIGGLE: nk:', nk
        write(*, *) 'INIT_WIGGLE: Splitting into wiggle and broad-band'
    end if

    call calculate_psmooth(k, z, Pk, Pk_smooth, cosm)

    if (HM_verbose) write(*, *) 'INIT_WIGGLE: Isolating wiggle'

    ! Isolate the wiggle in the look-up table (log_k_wiggle was filled above)
    call ensure_table_size(cosm%pk_wiggle, nk)
    cosm%pk_wiggle = Pk - Pk_smooth

    if (HM_verbose) then
        write(*, *) 'INIT_WIGGLE: Done'
        write(*, *)
    end if

    end subroutine init_wiggle

    subroutine calculate_psmooth(k, z, Pk, Pk_smt, cosm)
    ! Calculate the normalised smoothed power spectrum at a range of k
    real(dl), intent(in) :: k(:)
    real(dl), intent(in) :: z
    real(dl), intent(in) :: Pk(:)
    real(dl), allocatable, intent(out) :: Pk_smt(:)
    real(dl), allocatable :: Pk_nw(:)
    type(HM_cosmology), intent(in) :: cosm
    real(dl), parameter :: sig = wiggle_sigma

    ! Reduce dynamic range
    call calculate_nowiggle(k, z, Pk_nw, cosm)
    Pk_smt = Pk/Pk_nw

    ! Smooth linear power
    call smooth_array_Gaussian(log(k), Pk_smt, sig)

    ! Return dynamic range
    Pk_smt = Pk_smt*Pk_nw

    end subroutine calculate_psmooth

    subroutine calculate_nowiggle(k, z, Pk_nw, cosm)
    ! Calculate the normalised no wiggle power spectrum at a range of k and a
    ! Comes from the Eisenstein & Hu approximation
    real(dl), intent(in) :: k(:)
    real(dl), intent(in) :: z
    real(dl), allocatable, intent(out) :: Pk_nw(:)
    type(HM_cosmology), intent(in) :: cosm
    integer :: ik, nk
    real(dl) :: Pk_norm, Pk_nw_norm, s, alpha
    real(dl), parameter :: knorm = knorm_nowiggle
    integer, parameter :: type = 0 ! Matter here

    nk = size(k)
    allocate(Pk_nw(nk))

    ! Get the no-wiggle power spectrum (s and alpha do not depend on k)
    call Tk_nw_init(cosm, s, alpha)
    do ik = 1, nk
        Pk_nw(ik) = Pk_nowiggle(k(ik), s, alpha, cosm)
    end do

    ! Calculate the no-wiggle power spectrum and force spectra to agree at the normalisation wavenumber
    Pk_norm = p_lin(knorm, z, type, cosm)
    Pk_nw_norm = interp_cubic_sorted(knorm, k, Pk_nw, nk)
    Pk_nw = Pk_nw*Pk_norm/Pk_nw_norm

    end subroutine calculate_nowiggle

    real(dl) function Pk_nowiggle(k, s, alpha, cosm)
    ! Calculates the un-normalised no-wiggle power spectrum
    ! Comes from the Eisenstein & Hu approximation
    real(dl), intent(in) :: k
    real(dl), intent(in) :: s, alpha ! k-independent parameters from Tk_nw_init
    type(HM_cosmology), intent(in) :: cosm

    Pk_nowiggle = (k**(cosm%ns + 3.))*Tk_nw(k, s, alpha, cosm)**2

    end function Pk_nowiggle

    subroutine Tk_nw_init(cosm, s, alpha)
    ! The k-independent parameters of the no-wiggle transfer function; only needed once per cosmology
    type(HM_cosmology), intent(in) :: cosm
    real(dl), intent(out) :: s, alpha
    real(dl) :: wm, wb, rb

    wm = cosm%om_m*cosm%h**2 ! Real matter density
    wb = cosm%om_b*cosm%h**2 ! Real baryon density
    rb = cosm%om_b/cosm%om_m ! Baryon ratio

    s = 44.5*log(9.83/wm)/sqrt(1. + 10.*wb**0.75) ! Equation (26)
    alpha = 1. - 0.328*log(431.*wm)*rb + 0.38*log(22.3*wm)*rb**2 ! Equation (31)

    end subroutine Tk_nw_init

    real(dl) function Tk_nw(k, s, alpha, cosm)
    ! No-wiggle transfer function from Eisenstein & Hu: astro-ph:9709112
    real(dl), intent(in) :: k ! Wavenumber [h/Mpc]
    real(dl), intent(in) :: s, alpha ! k-independent parameters from Tk_nw_init
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: q, L, C, Gamma
    real(dl), parameter :: e = exp(1._dl)

    ! Functions of k
    Gamma = cosm%om_m*cosm%h*(alpha + (1. - alpha)/(1. + (0.43*k*s*cosm%h)**4)) ! Equation (30)
    q = k*(cosm%Tcmb/2.7)**2/Gamma ! Equation (28)
    L = log(2.*e + 1.8*q) ! Equation (29)
    C = 14.2 + 731./(1. + 62.5*q) ! Equation (29)
    Tk_nw = L/(L + C*q**2) ! Equation (29)

    end function Tk_nw

    subroutine smooth_array_Gaussian(x, f, sigma)
    ! Smooth an array f(x) using a Gaussian kernel. x is assumed equally spaced (guaranteed by
    ! the caller), so the kernel depends only on |i-j| and can be tabulated once, avoiding the
    ! n^2 exponentials of the general case
    real(dl), intent(in) :: x(:) ! x coordinates
    real(dl), intent(inout) :: f(:) ! Array to smooth
    real(dl), intent(in) :: sigma ! Width of smoothing Gaussian
    integer :: i, j, n
    real(dl), allocatable :: ff(:), kernel(:)
    real(dl) :: total, dx
    real(dl), parameter :: nsig = 3. ! Do not smooth if point lies within this number of sigma from edge

    if (sigma /= 0.) then

        n = size(x)
        if (n /= size(f)) error stop 'GAUSSIAN_SMOOTH_ARRAY: Error, x and y should be the same size'

        ! Save the original input array
        ff = f

        ! Delete the original array
        f = 0.

        ! Tabulate the kernel against index separation
        dx = (x(n) - x(1))/real(n - 1, dl)
        allocate(kernel(0:n - 1))
        do i = 0, n - 1
            kernel(i) = exp(-(i*dx)**2/(2.*sigma**2))
        end do

        ! Apply Gaussian smoothing
        do i = 1, n
            total = 0.
            if (abs(x(i) - x(1)) < nsig*sigma .or. abs(x(i) - x(n)) < nsig*sigma) then
                f(i) = ff(i)
            else
                do j = 1, n
                    f(i) = f(i) + ff(j)*kernel(abs(i - j))
                    total = total + kernel(abs(i - j))
                end do
                f(i) = f(i)/total
            end if
        end do

        deallocate(ff)

    end if

    end subroutine smooth_array_Gaussian

    subroutine fill_sigtab(this, cosm)
    class(THalofit) :: this
    ! This fills up HM_tables of r vs. sigma(r) across a range in r
    ! It is used only in look-up for further calculations of sigmac(r) and not otherwise
    ! and prevents a large number of calls to the sigint functions
    ! rmin and rmax need to be decided in advance and are chosen such that
    ! R vs. sigma(R) is approximately power-law below and above these values of R
    ! This wouldn't be appropriate for models with a small-scale linear spectrum cut-off (e.g., WDM)
    integer :: i
    type(HM_cosmology), intent(inout) :: cosm
    integer :: itype
    real(dl), parameter :: rmin = rmin_sigma_interpolation
    real(dl), parameter :: rmax = rmax_sigma_interpolation
    integer, parameter :: nsig = n_sigma_interpolation

    ! Choose type of sigma(R) to tabulate depending on HMcode version
    if (this%halofit_version == halofit_mead2015) then
        itype = 0 ! 0 - All matter
    else
        itype = 1 ! 1 - Cold matter
    end if

    ! These values of 'r' work fine for any power spectrum of cosmological importance
    ! Fill log(R) directly (rather than as log(r)) so that it is exactly equally spaced,
    ! as assumed by interp_cubic_uniform
    call fill_table(log(rmin), log(rmax), cosm%log_r_sigma, nsig)
    call ensure_table_size(cosm%log_sigma, nsig)

    if (HM_verbose) write(*, *) 'SIGTAB: Filling sigma interpolation table'
    if (HM_verbose) write(*, *) 'SIGTAB: R_min:', rmin
    if (HM_verbose) write(*, *) 'SIGTAB: R_max:', rmax
    if (HM_verbose) write(*, *) 'SIGTAB: Values:', nsig

    ! Cost per point varies a lot with R (adaptive integration), so balance dynamically
    !$OMP PARALLEL DO IF(HM_par_inner()), DEFAULT(SHARED), SCHEDULE(DYNAMIC)
    do i = 1, nsig
        cosm%log_sigma(i) = log(sigma_integral(exp(cosm%log_r_sigma(i)), 0.d0, itype, cosm))
    end do
    !$OMP END PARALLEL DO

    if (HM_verbose) write(*, *) 'SIGTAB: sigma_min:', exp(cosm%log_sigma(nsig))
    if (HM_verbose) write(*, *) 'SIGTAB: sigma_max:', exp(cosm%log_sigma(1))

    cosm%nsig = nsig

    if (HM_verbose) write(*, *) 'SIGTAB: Done'
    if (HM_verbose) write(*, *)

    end subroutine fill_sigtab

    function sigma_lut(r, z, cosm)
    ! Finds sigma_cold(R) from look-up table
    real(dl) :: sigma_lut
    real(dl), intent(in) :: r, z
    type(HM_cosmology), intent(in) :: cosm

    sigma_lut = cached_grow(z, cosm)*exp(interp_cubic_uniform(log(r), cosm%log_r_sigma, cosm%log_sigma, cosm%nsig))

    end function sigma_lut

    function wk_tophat(x)

    ! The normlaised Fourier Transform of a top-hat
    real(dl) :: wk_tophat
    real(dl), intent(in) :: x
    real(dl), parameter :: dx = 1e-3 ! Taylor expansion for |x|<dx

    ! Taylor expansion used for low x to avoid cancellation problems
    if (abs(x) < dx) then
        wk_tophat = 1 - (x**2)/10
    else
        wk_tophat = 3*(sin(x) - x*cos(x))/(x**3)
    end if

    end function wk_tophat

    function wk_tophat_deriv(x)
    ! The derivative of a normlaised Fourier Transform of a spherical top-hat
    real(dl) :: wk_tophat_deriv
    real(dl), intent(in) :: x
    real(dl), parameter :: dx = 1e-3 ! Taylor expansion for |x|<dx

    ! Taylor expansion used for low x to avoid cancellation problems
    if (abs(x) < dx) then
        wk_tophat_deriv = -x/5 + x**3/70
    else
        wk_tophat_deriv = (3/x**4)*((x**2 - 3)*sin(x) + 3*x*cos(x))
    end if

    end function wk_tophat_deriv

    function inttab(x, y, n1, n2)
    ! Integrates the table y(x)dx from x(n1) to x(n2) using a cubic through each section
    real(dl) :: inttab
    integer, intent(in) :: n1, n2
    real(dl), intent(in) :: x(:), y(:)
    real(dl) :: a, b, c, d
    real(dl) :: qi, qf
    real(dl) :: sum
    integer :: i, i1, n

    n = size(x)

    sum = 0.d0

    do i = n1, n2 - 1

        ! First choose the integers used for defining cubics for each section
        ! First and last are different because the section does not lie in the *middle* of a cubic
        i1 = min(max(i - 1, 1), n - 3)

        call fit_cubic(a, b, c, d, x(i1), y(i1), x(i1 + 1), y(i1 + 1), x(i1 + 2), y(i1 + 2), x(i1 + 3), y(i1 + 3))

        ! These are the limits of the particular section of integral
        qi = a*(x(i)**4)/4. + b*(x(i)**3)/3. + c*(x(i)**2)/2. + d*x(i)
        qf = a*(x(i + 1)**4)/4. + b*(x(i + 1)**3)/3. + c*(x(i + 1)**2)/2. + d*x(i + 1)

        sum = sum + qf - qi

    end do

    inttab = sum

    end function inttab

    function sigma_integral(r, z, itype, cosm)
    ! Gets sigma(R)
    real(dl) :: sigma_integral
    real(dl), intent(in) :: r, z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl), parameter :: tmin = 0.
    real(dl), parameter :: tmax = 1.
    real(dl), parameter :: acc = acc_sigma_integration

    sigma_integral = sqrt(integrate(tmin, tmax, sigma_integrand, r, z, itype, cosm, acc))

    end function sigma_integral

    function sigma_integrand(t, R, z, itype, cosm)
    ! The integrand for the sigma(R) integrals
    real(dl) :: sigma_integrand
    real(dl), intent(in) :: t, R, z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: k, kR, w_hat
    real(dl), parameter :: alpha = alpha_sigma_integration

    if (t <= 0. .or. t >= 1.) then
        ! t = 0 corresponds to k = infintiy when W(kR) = 0.
        ! t = 1 corresponds to k = 0. when P(k) = 0.
        sigma_integrand = 0.d0
    else
        kR = (-1. + 1./t)**alpha
        k = kR/R
        w_hat = wk_tophat(kR)
        sigma_integrand = p_lin(k, z, itype, cosm)*(w_hat**2)*alpha/(t*(1. - t))
    end if

    end function sigma_integrand

    function integrate(a, b, f, y, z, itype, cosm, acc)
    ! Integrates between a and b by Simpson's rule until the desired accuracy is reached
    ! Stores information to reduce function calls
    real(dl) :: integrate
    real(dl), intent(in) :: a
    real(dl), intent(in) :: b
    real(dl), intent(in) :: y
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl), intent(in) :: acc
    integer :: i, j
    integer :: n
    real(dl) :: x, dx
    real(dl) :: f1, f2, fx
    real(dl) :: sum_n, sum_2n, sum_new, sum_old
    integer, parameter :: jmin = 5
    integer, parameter :: jmax = 20

    interface
    function f(x, y, z, itype, cosm)
    use precision
    import :: HM_cosmology
    real(dl) :: f
    real(dl), intent(in) :: x
    real(dl), intent(in) :: y
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    end function f
    end interface

    if (a == b) then

        ! Fix the answer to zero if the integration limits are identical
        integrate = 0.d0

    else

        ! Reset the sum variable for the integration
        sum_2n = 0.d0
        sum_n = 0.d0
        sum_old = 0.d0
        sum_new = 0.d0

        do j = 1, jmax

            ! Note, you need this to be 1+2**n for some integer n
            ! j = 1 n = 2; j = 2 n = 3; j = 3 n = 5; j = 4 n = 9; ...'
            n = 1 + 2**(j - 1)

            ! Calculate the dx interval for this value of 'n'
            dx = (b - a)/real(n - 1, dl)

            if (j == 1) then

                ! The first go is just the trapezium of the end points
                f1 = f(a, y, z, itype, cosm)
                f2 = f(b, y, z, itype, cosm)
                sum_2n = 0.5d0*(f1 + f2)*dx
                sum_new = sum_2n

            else

                ! Loop over only new even points to add these to the integral
                do i = 2, n, 2
                    x = a + (b - a)*real(i - 1, dl)/real(n - 1, dl)
                    fx = f(x, y, z, itype, cosm)
                    sum_2n = sum_2n + fx
                end do

                ! Now create the total using the old and new parts
                sum_2n = sum_n/2.d0 + sum_2n*dx

                ! This is Simpson's rule, which cancels the leading trapezium error
                sum_new = (4.d0*sum_2n - sum_n)/3.d0

            end if

            if (j >= jmin .and. sum_old /= 0) then
                if (abs(-1.d0 + sum_new/sum_old) < acc) then
                    ! jmin avoids spurious early convergence
                    integrate = sum_new
                    exit
                end if
            end if
            if (j == jmax) then
                integrate = sum_new
                call GlobalError('HMCode INTEGRATE, Integration timed out', error_nonlinear)
                return
            else
                ! Integral has not converged so store old sums and reset sum variables
                sum_old = sum_new
                sum_n = sum_2n
                sum_2n = 0.d0
            end if

        end do

    end if

    end function integrate

    function win(k, rs, c, norm, moment2)
    ! The halo window function (k-space halo profile): the analytic Fourier transform of the
    ! NFW profile of scale radius rs and concentration c, normalised so that W(k->0)=1.
    ! norm is 1/(log(1+c)-c/(1+c)), the reciprocal of the halo mass in units of 4*pi*rho_n*rs^3,
    ! where rho_n is the profile normalisation [i.e. rho=rho_n/((r/rs)*(1+r/rs)^2]
    real(dl) :: win
    real(dl), intent(in) :: k, rs, c, norm, moment2
    real(dl) :: si1, si2, ci1, ci2, sin1, sin2, cos1, cos2, ks, ks2

    ks = k*rs
    ks2 = (1 + c)*ks

    ! The small-k expansion is the second radial moment of the NFW profile. Its
    ! parameter is k*r_vir=c*ks; below 0.2 the pointwise fourth-order remainder is
    ! bounded by (k*r_vir)^4/120 < 1.4e-5, safely below the HMcode accuracy target.
    if (c*ks < 0.2_dl) then
        win = 1._dl - ks*ks*moment2*norm/6._dl
        ! At large arguments SiCi already evaluates the sine and cosine needed below.
        ! Return those values and get sin(c*ks)=sin(ks2-ks) by angle subtraction, avoiding
        ! three otherwise redundant trigonometric calls in the common high-k branch.
    else if (ks2 > 4._dl) then
        sin1 = sin(ks); cos1 = cos(ks)
        sin2 = sin(ks2); cos2 = cos(ks2)
        call SiCi(ks, si1, ci1, sin1, cos1)
        call SiCi(ks2, si2, ci2, sin2, cos2)
        win = (cos1*(ci2 - ci1) + sin1*(si2 - si1) - (sin2*cos1 - cos2*sin1)/ks2)*norm
    else
        call SiCi(ks, si1, ci1)
        call SiCi(ks2, si2, ci2)
        win = (cos(ks)*(ci2 - ci1) + sin(ks)*(si2 - si1) - sin(c*ks)/ks2)*norm
    end if

    ! Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    if (win > 1._dl) win = 1._dl
    if (win < 0._dl) win = 0._dl

    end function win

    function gst(nu)
    ! Sheth & Tormen (1999) mass function!
    real(dl) :: gst
    real(dl), intent(in) :: nu
    real(dl), parameter :: p = 0.3
    real(dl), parameter :: a = 0.707
    real(dl), parameter :: bigA = 0.21616

    ! Note I use nu=dc/sigma whereas ST (1999) use nu=(dc/sigma)^2
    ! This accounts for the different pre-factor and slightly changed nu dependence
    ! f(nu^2)d(nu^2) = 2*nu*f(nu)dnu

    ! Full mass function. Note this is normalised such that integral f(nu)dnu = 1
    gst = bigA*(1 + ((a*nu*nu)**(-p)))*exp(-a*nu*nu/2)

    end function gst

    function Hubble2(z, cosm)
    ! This calculates the dimensionless squared hubble parameter at redshift z (H/H_0)^2!
    ! Ignores contributions from radiation (not accurate at high z, but consistent with simulations)!
    real(dl) :: Hubble2
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm

    Hubble2 = cosm%om_m*(1 + z)**3 + cosm%om_v*X_de(z, cosm) + (1 - cosm%om_m - cosm%om_v)*(1 + z)**2

    end function Hubble2

    function AH(z, cosm)
    ! The Hubble acceleration function \ddot{a}/a
    real(dl) :: AH
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm

    AH = cosm%om_m*(1 + z)**3 + cosm%om_v*(1. + 3.*w_de_hm(z, cosm))*X_de(z, cosm)
    AH = -AH/2.

    end function AH

    function X_de(z, cosm)
    ! The time evolution for dark energy: rho_de=rho_de,0 * X(a)
    ! X(a) = 1 for LCDM but changes for other models
    real(dl) :: X_de
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: a

    a = 1./(1. + z)
    X_de = (a**(-3*(1 + cosm%w + cosm%wa)))*exp(-3*cosm%wa*(1 - a))

    end function X_de

    function w_de_hm(z, cosm)

    ! The dark energy w(a) function
    real(dl) :: w_de_hm
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: a

    a = 1./(1. + z)
    w_de_hm = cosm%w + (1 - a)*cosm%wa

    end function w_de_hm

    function Omega_m_hm(z, cosm)
    ! This calculates omega_m variations with z!
    real(dl) :: Omega_m_hm
    real(dl), intent(in) :: z
    real(dl) :: om_m
    type(HM_cosmology), intent(in) :: cosm

    om_m = cosm%om_m
    Omega_m_hm = (om_m*(1 + z)**3)/Hubble2(z, cosm)

    end function Omega_m_hm

    function Omega_cold_hm(z, cosm)
    ! This calculates omega_cold variations with z (no neutrinos)
    real(dl) :: Omega_cold_hm
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm

    Omega_cold_hm = ((cosm%om_c + cosm%om_b)*(1 + z)**3)/Hubble2(z, cosm)

    end function Omega_cold_hm

    function grow(z, cosm)
    ! Finds the scale-independent growth function at redshift z
    real(dl) :: grow
    real(dl), intent(in) :: z
    real(dl) :: a
    type(HM_cosmology), intent(in) :: cosm

    if (z == 0.) then
        grow = 1.
    else
        a = 1./(1. + z)
        grow = interp_cubic_uniform(a, cosm%a_growth, cosm%growth, cosm%ng)
    end if

    end function grow

    real(dl) recursive function ungrow(z, cosm)

    ! Growth function normalised such that g(a) = a at early (matter-dominated) times
    implicit none
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm

    ungrow = cosm%gnorm*grow(z, cosm)

    end function ungrow

    real(dl) recursive function acc_growth(z, cosm)

    ! Accumulated growth function: int_0^a g(a)/a da
    implicit none
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: a

    a = 1./(1. + z)
    acc_growth = interp_cubic_uniform(a, cosm%a_growth, cosm%agrow, cosm%ng)

    end function acc_growth

    function sigmaV(R, z, itype, cosm)
    real(dl) :: sigmaV
    real(dl), intent(in) :: R
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl), parameter :: tmin = 0._dl
    real(dl), parameter :: tmax = 1._dl
    real(dl), parameter :: acc = acc_sigmaV_integration

    sigmaV = sqrt(integrate(tmin, tmax, sigmaV_integrand, R, z, itype, cosm, acc)/3)

    end function sigmaV

    function sigmaV_integrand(t, R, z, itype, cosm)
    ! This is the integrand for the velocity dispersion integral
    real(dl) :: sigmaV_integrand
    real(dl), intent(in) :: t
    real(dl), intent(in) :: R
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: k, kR, w_hat
    real(dl), parameter :: alpha = alpha_sigmaV_integration ! Speeds up integral for large 'R'

    if (t <= 0._dl .or. t >= 1._dl) then
        sigmaV_integrand = 0
    else
        if (R == 0.) then
            kR = 0
            k = (-1 + 1/t)**alpha
        else
            kR = (-1 + 1/t)**alpha
            k = kR/R
        end if
        w_hat = wk_tophat(kR)
        sigmaV_integrand = (p_lin(k, z, itype, cosm)/k**2)*(w_hat**2)*alpha/(t*(1 - t))
    end if

    end function sigmaV_integrand

    function neff_integrand(t, R, z, itype, cosm)
    ! This is the integrand for the velocity dispersion integral
    real(dl) :: neff_integrand
    real(dl), intent(in) :: t
    real(dl), intent(in) :: R
    real(dl), intent(in) :: z
    integer, intent(in) :: itype
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: k, kR, w_hat, w_hat_deriv
    real(dl), parameter :: alpha = alpha_neff_integration ! Speeds up integral for large 'R'

    if (t <= 0. .or. t >= 1.) then
        neff_integrand = 0.
    else
        kR = (-1. + 1./t)**alpha
        k = kR/R
        w_hat = wk_tophat(kR)
        w_hat_deriv = wk_tophat_deriv(kR)
        neff_integrand = p_lin(k, z, itype, cosm)*w_hat*w_hat_deriv*alpha*kR/(t*(1. - t))
    end if

    end function neff_integrand

    subroutine SiCi(x, Six, Cix, sinx, cosx)

    ! Calculates the 'sine integral' Si(x) and 'cosine integral' Ci(x) together
    ! The large-x expansions share the same auxiliary functions f and g, so are done once
    real(dl), intent(in) :: x
    real(dl), intent(out) :: Six, Cix
    real(dl), intent(in), optional :: sinx, cosx
    real(dl) :: x2, y, f, g, sin_x, cos_x
    real(dl), parameter :: em_const = 0.577215664901532861d0
    real(dl), parameter :: x0 = 4.d0 ! Transition between two different approximations

    ! Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    if (abs(x) <= x0) then

        x2 = x*x

        Six = x*(1.d0 + x2*(-4.54393409816329991d-2 + x2*(1.15457225751016682d-3 &
            + x2*(-1.41018536821330254d-5 + x2*(9.43280809438713025d-8 + x2*(-3.53201978997168357d-10 &
            + x2*(7.08240282274875911d-13 + x2*(-6.05338212010422477d-16))))))))/ &
            (1. + x2*(1.01162145739225565d-2 + x2*(4.99175116169755106d-5 + &
            x2*(1.55654986308745614d-7 + x2*(3.28067571055789734d-10 + x2*(4.5049097575386581d-13 &
            + x2*(3.21107051193712168d-16)))))))

        Cix = em_const + log(x) + x2*(-0.25d0 + x2*(7.51851524438898291d-3 + x2*(-1.27528342240267686d-4 &
            + x2*(1.05297363846239184d-6 + x2*(-4.68889508144848019d-9 + x2*(1.06480802891189243d-11 &
            + x2*(-9.93728488857585407d-15)))))))/(1. + x2*(1.1592605689110735d-2 + &
            x2*(6.72126800814254432d-5 + x2*(2.55533277086129636d-7 + x2*(6.97071295760958946d-10 + &
            x2*(1.38536352772778619d-12 + x2*(1.89106054713059759d-15 + x2*(1.39759616731376855d-18))))))))

    else

        y = 1.d0/(x*x)

        f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 + &
            y*(2.37750310125431834034d7 + y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10 &
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + &
            y*(1.00795182980368574617d13 + y*(4.94816688199951963482d12 + &
            y*(-4.94701168645415959931d11)))))))))))/(x*(1. + y*(7.46437068161927678031d2 + &
            y*(1.97865247031583951450d5 + y*(2.41535670165126845144d7 + &
            y*(1.47478952192985464958d9 + y*(4.58595115847765779830d10 + &
            y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 + &
            y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))

        g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + &
            y*(3.12557570795778731d7 + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 + &
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 + y*(1.81004487464664575d13 + &
            y*(6.43291613143049485d12 + y*(-1.36517137670871689d12)))))))))))/ &
            (1. + y*(8.19595201151451564d2 + y*(2.40036752835578777d5 + y*(3.26026661647090822d7 &
            + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 + y*(1.39866710696414565d12 &
            + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 + y*(3.99653257887490811d13))))))))))

        if (present(sinx) .and. present(cosx)) then
            sin_x = sinx
            cos_x = cosx
        else
            sin_x = sin(x)
            cos_x = cos(x)
        end if

        Six = pi_HM/2.d0 - f*cos_x - g*sin_x
        Cix = f*sin_x - g*cos_x

    end if

    end subroutine SiCi

    ! Interpolation in look-up tables that are in increasing order of x. Values off either end
    ! are extrapolated linearly from the two points at that end. The _uniform routines assume
    ! the table is exactly equally spaced in x, the _sorted ones only that it is increasing.
    ! Interpolate log(x) and/or log(y) where that gives better results.

    pure function lin_interp(x, x1, y1, x2, y2) result(y)
    ! Linear interpolation (or extrapolation) through two points
    real(dl) :: y
    real(dl), intent(in) :: x, x1, y1, x2, y2

    y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)

    end function lin_interp

    function interp_linear_uniform(x, xtab, ytab, n) result(y)
    ! Linear interpolation in an equally spaced table
    real(dl) :: y
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n), ytab(n)
    integer :: i

    ! linear_table_integer clamps to the end intervals, which extrapolates off the table
    i = linear_table_integer(x, xtab, n)
    y = lin_interp(x, xtab(i), ytab(i), xtab(i + 1), ytab(i + 1))

    end function interp_linear_uniform

    function interp_cubic_uniform(x, xtab, ytab, n) result(y)
    ! Cubic interpolation in an equally spaced table
    real(dl) :: y
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n), ytab(n)

    if (x < xtab(1) .or. x > xtab(n)) then
        y = interp_off_table(x, xtab, ytab, n)
    else
        y = cubic_Lagrange(x, xtab, ytab, linear_table_integer(x, xtab, n), n)
    end if

    end function interp_cubic_uniform

    function interp_cubic_sorted(x, xtab, ytab, n) result(y)
    ! Cubic interpolation in an increasing table, located by bisection
    real(dl) :: y
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n), ytab(n)

    if (x < xtab(1) .or. x > xtab(n)) then
        y = interp_off_table(x, xtab, ytab, n)
    else
        y = cubic_Lagrange(x, xtab, ytab, int_split(x, xtab, n), n)
    end if

    end function interp_cubic_sorted

    function interp_off_table(x, xtab, ytab, n) result(y)
    ! Linear extrapolation from whichever end of the table x lies beyond
    real(dl) :: y
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n), ytab(n)
    integer :: i

    if (x < xtab(1)) then
        i = 1
    else
        i = n - 1
    end if
    y = lin_interp(x, xtab(i), ytab(i), xtab(i + 1), ytab(i + 1))

    end function interp_off_table

    function cubic_Lagrange(x, xtab, ytab, i, n) result(y)
    ! Cubic Lagrange polynomial through the four points bracketing the interval starting at i,
    ! shifted at the ends of the table where the interval is not in the middle of the four
    real(dl) :: y
    integer, intent(in) :: i, n
    real(dl), intent(in) :: x, xtab(n), ytab(n)
    integer :: j
    real(dl) :: dx(4)

    j = min(max(i - 1, 1), n - 3)
    associate(xv => xtab(j:j + 3), yv => ytab(j:j + 3))
        dx = x - xv
        y = dx(2)*dx(3)*dx(4)/(xv(1) - xv(2))/(xv(1) - xv(3))/(xv(1) - xv(4))*yv(1) &
            + dx(1)*dx(3)*dx(4)/(xv(2) - xv(1))/(xv(2) - xv(3))/(xv(2) - xv(4))*yv(2) &
            + dx(1)*dx(2)*dx(4)/(xv(3) - xv(1))/(xv(3) - xv(2))/(xv(3) - xv(4))*yv(3) &
            + dx(1)*dx(2)*dx(3)/(xv(4) - xv(1))/(xv(4) - xv(2))/(xv(4) - xv(3))*yv(4)
    end associate

    end function cubic_Lagrange

    function linear_table_integer(x, xtab, n)
    ! Assuming the table is exactly linear this gives you the integer position
    ! Clamped because the tabulated points can differ from the exactly equally spaced
    ! values by a last bit, which could otherwise put x in a neighbouring cell
    integer :: linear_table_integer
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n)

    linear_table_integer = 1 + floor(real(n - 1, dl)*(x - xtab(1))/(xtab(n) - xtab(1)))
    linear_table_integer = max(1, min(n - 1, linear_table_integer))

    end function linear_table_integer

    function int_split(x, xtab, n)
    ! Finds the position of the value in the table by continually splitting it in half
    integer :: int_split
    integer, intent(in) :: n
    real(dl), intent(in) :: x, xtab(n)
    integer :: i1, i2, imid

    if (xtab(1) > xtab(n)) error stop 'INT_SPLIT: table in wrong order'

    i1 = 1
    i2 = n

    do

        imid = (i1 + i2)/2

        if (x < xtab(imid)) then
            i2 = imid
        else
            i1 = imid
        end if

        if (i2 == i1 + 1) exit

    end do

    int_split = i1

    end function int_split

    subroutine fit_cubic(a, b, c, d, x1, y1, x2, y2, x3, y3, x4, y4)
    ! Given xi, yi i=1,2,3,4 fits a cubic between these points
    real(dl), intent(out) :: a, b, c, d
    real(dl), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
    real(dl) :: f1, f2, f3

    f1 = (y4 - y1)/((x4 - x2)*(x4 - x1)*(x4 - x3))
    f2 = (y3 - y1)/((x3 - x2)*(x3 - x1)*(x4 - x3))
    f3 = (y2 - y1)/((x2 - x1)*(x4 - x3))*(1./(x4 - x2) - 1./(x3 - x2))

    a = f1 - f2 - f3

    f1 = (y3 - y1)/((x3 - x2)*(x3 - x1))
    f2 = (y2 - y1)/((x2 - x1)*(x3 - x2))
    f3 = a*(x3 + x2 + x1)

    b = f1 - f2 - f3

    f1 = (y4 - y1)/(x4 - x1)
    f2 = a*(x4**2 + x4*x1 + x1**2)
    f3 = b*(x4 + x1)

    c = f1 - f2 - f3

    d = y1 - a*x1**3 - b*x1**2 - c*x1

    end subroutine fit_cubic

    subroutine fill_growtab(cosm)
    ! Fills a table of values of the scale-independent growth function
    type(HM_cosmology) :: cosm
    integer :: i
    real(dl), parameter :: amin = amin_growth_interpolation
    real(dl), parameter :: amax = amax_growth_interpolation
    integer, parameter :: n = n_growth_interpolation
    real(dl) :: g_over_a(n)

    ! Equally spaced in a, as assumed by interp_cubic_uniform
    cosm%ng = n
    call fill_table(amin, amax, cosm%a_growth, n)
    if (amax /= 1) error stop 'FILL_GROWTAB: the growth table must end at a=1 to normalise there'
    call ensure_table_size(cosm%growth, n)

    if (HM_verbose) write(*, *) 'GROWTH: Solving growth equation'
    call solve_growth_ODE(cosm, cosm%a_growth, cosm%growth)
    if (HM_verbose) write(*, *) 'GROWTH: ODE done'

    ! Normalise so that g(z=0)=1
    cosm%gnorm = cosm%growth(n)
    if (HM_verbose) write(*, *) 'GROWTH: Unnormalised g(a=1):', cosm%gnorm
    cosm%growth = cosm%growth/cosm%gnorm

    ! Table integration to calculate G(a)=int_0^a g(a')/a' da'
    call ensure_table_size(cosm%agrow, n)

    ! Each interval contributes the same however far the integral goes, so accumulate them.
    ! The missing section below a_growth(1): g(a=0)/0 = 1, so just add a rectangle of height g*a/a=g
    g_over_a = cosm%gnorm*cosm%growth/cosm%a_growth
    cosm%agrow(1) = cosm%gnorm*cosm%growth(1)
    do i = 2, n
        cosm%agrow(i) = cosm%agrow(i - 1) + inttab(cosm%a_growth, g_over_a, i - 1, i)
    end do

    if (HM_verbose) write(*, *) 'GROWTH: Accumulated G(a=1):', cosm%agrow(n)
    if (HM_verbose) write(*, *) 'GROWTH: Done'
    if (HM_verbose) write(*, *)

    end subroutine fill_growtab

    subroutine solve_growth_ODE(cosm, a, d)
    ! Solves the linear growth ODE for d(a) at the increasing output points a(:), by RK4 with
    ! nsub_growth_ODE steps per output interval. In matter domination d is proportional to a,
    ! which RK4 integrates exactly, so a fixed step in a is accurate even at early times.
    type(HM_cosmology), intent(in) :: cosm
    real(dl), intent(in) :: a(:)
    real(dl), intent(out) :: d(:)
    real(dl) :: x, v, t, h, f
    integer :: i, j

    ! These set the initial conditions to be the Om_m=1. growing mode
    ! AM Jul 19: changed initial conditions to be appropriate for massive neutrino cosmologies
    ! AM Sep 20: changed initial conditions to assume neutrinos cluster, but changed to take into account EDE
    t = aini_growth_ODE
    f = 1. - Omega_m_hm(-1. + 1./t, cosm)
    x = t**(1. - 3.*f/5.)
    v = (1. - 3.*f/5.)*t**(-3.*f/5.)

    do i = 1, size(a)
        h = (a(i) - t)/nsub_growth_ODE
        do j = 1, nsub_growth_ODE
            call rk4_growth_step(x, v, t, h, cosm)
        end do
        t = a(i) ! the accumulated t can differ in the last bits
        d(i) = x
    end do

    end subroutine solve_growth_ODE

    subroutine rk4_growth_step(d, v, a, h, cosm)
    ! One 4th-order Runge-Kutta step of the growth ODE d''(a) = fv(d,d',a), with v = d'
    real(dl), intent(inout) :: d, v, a
    real(dl), intent(in) :: h
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: kd1, kd2, kd3, kd4, kv1, kv2, kv3, kv4

    kd1 = h*v
    kv1 = h*fv(d, v, a, cosm)
    kd2 = h*(v + kv1/2)
    kv2 = h*fv(d + kd1/2, v + kv1/2, a + h/2, cosm)
    kd3 = h*(v + kv2/2)
    kv3 = h*fv(d + kd2/2, v + kv2/2, a + h/2, cosm)
    kd4 = h*(v + kv3)
    kv4 = h*fv(d + kd3, v + kv3, a + h, cosm)

    d = d + (kd1 + 2*kd2 + 2*kd3 + kd4)/6
    v = v + (kv1 + 2*kv2 + 2*kv3 + kv4)/6
    a = a + h

    end subroutine rk4_growth_step

    function fv(d, v, a, cosm)

    ! d'' as a function of the growth d, its derivative v=d' and the scale factor
    real(dl) :: fv
    real(dl), intent(in) :: d, v, a
    real(dl) :: f1, f2, z, Om_m
    type(HM_cosmology), intent(in) :: cosm

    z = -1. + (1./a)

    ! AM Jul 19: changed Omega_m to Omega_cold for massive neutrino cosmologies
    if (cold_growth) then
        Om_m = Omega_cold_hm(z, cosm)
    else
        Om_m = Omega_m_hm(z, cosm)
    end if
    f1 = 3.*Om_m*d/(2.*a**2)
    f2 = (2. + AH(z, cosm)/Hubble2(z, cosm))*(v/a)

    fv = f1 - f2

    end function fv

    subroutine Mead_growth_terms(z, cosm, x, y, Om_m)

    ! Growth terms shared by the Mead (2017; 1606.05345) delta_c and Delta_v fitting functions
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl), intent(out) :: x, y, Om_m
    real(dl) :: a

    if (cold_growth) error stop 'MEAD_GROWTH_TERMS: Error, this will not work if you want cold growth'

    a = 1./(1. + z)
    x = ungrow(z, cosm)/a
    y = acc_growth(z, cosm)/a
    Om_m = Omega_m_hm(z, cosm)

    end subroutine Mead_growth_terms

    real(dl) function dc_Mead(z, cosm)

    ! delta_c fitting function from Mead (2017; 1606.05345)
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: x, y, Om_m

    ! See Appendix A of Mead (2017) for naming convention
    real(dl), parameter :: p10 = -0.0069
    real(dl), parameter :: p11 = -0.0208
    real(dl), parameter :: p12 = 0.0312
    real(dl), parameter :: p13 = 0.0021
    integer, parameter :: a1 = 1
    real(dl), parameter :: p20 = 0.0001
    real(dl), parameter :: p21 = -0.0647
    real(dl), parameter :: p22 = -0.0417
    real(dl), parameter :: p23 = 0.0646
    ! a2 = 0 in the paper, so log10(Om_m)**a2 = 1 and is omitted below
    real(dl), parameter :: dc0 = (3./20.)*(12.*pi_HM)**(2./3.)

    call Mead_growth_terms(z, cosm, x, y, Om_m)

    dc_Mead = 1.
    dc_Mead = dc_Mead + f_Mead(x, y, p10, p11, p12, p13)*log10(Om_m)**a1
    dc_Mead = dc_Mead + f_Mead(x, y, p20, p21, p22, p23)
    dc_Mead = dc_Mead*dc0*(1. - 0.041*cosm%f_nu)

    end function dc_Mead

    real(dl) function Dv_Mead(z, cosm)

    ! Delta_v fitting function from Mead (2017; 1606.05345)
    real(dl), intent(in) :: z
    type(HM_cosmology), intent(in) :: cosm
    real(dl) :: x, y, Om_m

    ! See Appendix A of Mead (2017) for naming convention
    real(dl), parameter :: p30 = -0.79
    real(dl), parameter :: p31 = -10.17
    real(dl), parameter :: p32 = 2.51
    real(dl), parameter :: p33 = 6.51
    integer, parameter :: a3 = 1
    real(dl), parameter :: p40 = -1.89
    real(dl), parameter :: p41 = 0.38
    real(dl), parameter :: p42 = 18.8
    real(dl), parameter :: p43 = -15.87
    integer, parameter :: a4 = 2
    real(dl), parameter :: Dv0 = 18.*pi_HM**2

    call Mead_growth_terms(z, cosm, x, y, Om_m)

    Dv_Mead = 1.
    Dv_Mead = Dv_Mead + f_Mead(x, y, p30, p31, p32, p33)*log10(Om_m)**a3
    Dv_Mead = Dv_Mead + f_Mead(x, y, p40, p41, p42, p43)*log10(Om_m)**a4
    Dv_Mead = Dv_Mead*Dv0*(1. + 0.763*cosm%f_nu)

    end function Dv_Mead

    real(dl) function f_Mead(x, y, p0, p1, p2, p3)

    ! Equation A3 in Mead (2017)
    implicit none
    real(dl), intent(in) :: x, y
    real(dl), intent(in) :: p0, p1, p2, p3

    f_Mead = p0 + p1*(1. - x) + p2*(1. - x)**2 + p3*(1. - y)

    end function f_Mead

    !!AM End HMcode

    subroutine PKequal(State, redshift, w_eff, wa_eff, w_hf, wa_hf)
    ! used by halofit_casarini: arXiv:0810.0190, arXiv:1601.07230
    ! Solve for the constant-w model (w_hf, wa_hf=0) that has the same comoving distance
    ! between redshift and last scattering as the actual dark energy model, keeping all the
    ! other densities (and the present dark energy density) fixed. The Takahashi halofit
    ! fit is then evaluated for that equivalent constant-w model.
    ! w_eff, wa_eff are the effective w0-wa of the actual model, used only to bracket the root.
    class(CAMBdata), target :: State
    real(dl), intent(in) :: redshift, w_eff, wa_eff
    real(dl), intent(out) :: w_hf, wa_hf
    ! Range of constant w searched; the equivalent w need not be accelerating, but for w>=1/3
    ! the dark energy would no longer be subdominant at last scattering
    real(dl), parameter :: w_search_min = -20._dl, w_search_max = 0._dl
    integer, parameter :: max_widen = 20
    type(TEquivalent_wCDM) :: trial
    real(dl) :: w1, w2, f1, f2, fzero
    integer :: iflag, i

    wa_hf = 0._dl
    w_hf = w_eff
    trial%State => State
    trial%a_star = 1/(1 + State%ThermoDerivedParams(derived_zstar))
    trial%a_z = 1/(1 + redshift)
    trial%tol = base_tol/1000/exp(State%CP%Accuracy%AccuracyBoost*State%CP%Accuracy%IntTolBoost - 1)
    trial%dlsb = State%DeltaTime(trial%a_star, trial%a_z, trial%tol)

    ! distance_error decreases monotonically with w (more dark energy early on means a larger
    ! H and so a smaller distance), so widen the bracket in whichever direction is short
    w1 = min(w_eff, w_eff + wa_eff) - 0.1_dl
    w2 = max(w_eff, w_eff + wa_eff) + 0.1_dl
    f1 = distance_error(trial, w1)
    f2 = distance_error(trial, w2)
    do i = 1, max_widen
        if (f1*f2 <= 0 .or. global_error_flag /= 0) exit
        if (w1 <= w_search_min .and. w2 >= w_search_max) exit
        if (f1 < 0) then
            w1 = max(w_search_min, w1 - max(0.5_dl, w2 - w1))
            f1 = distance_error(trial, w1)
        else
            w2 = min(w_search_max, w2 + max(0.5_dl, w2 - w1))
            f2 = distance_error(trial, w2)
        end if
    end do

    if (global_error_flag == 0) then
        if (f1*f2 > 0) then
            call GlobalError('halofit_casarini: no equivalent constant-w model in ' // &
                'searched range; try another halofit_version', error_unsupported_params)
        else
            call brentq(trial, distance_error, w1, w2, 1e-6_dl, w_hf, fzero, iflag, f1, f2)
            if (iflag /= 0) call GlobalError('halofit_casarini: equivalent constant-w solve failed', &
                error_unsupported_params)
        end if
    end if
    if (FeedbackLevel > 1) write(*, '(a,f8.3,a,f9.5)') &
        ' PKequal: at z = ', redshift, ' equivalent w_const = ', w_hf

    end subroutine PKequal

    function distance_error(obj, w) result(error)
    ! Fractional difference between the comoving distance from a_z to a_star in a model with
    ! constant dark energy equation of state w and the same distance in the actual model
    class(*) :: obj
    real(dl), intent(in) :: w
    real(dl) :: error

    select type (this => obj)
    type is (TEquivalent_wCDM)
        this%w = w
        error = Integrate_Romberg(this, dtauda_wconst, this%a_star, this%a_z, this%tol)/this%dlsb - 1
    class default
        error stop 'distance_error: expected TEquivalent_wCDM'
    end select

    end function distance_error

    function dtauda_wconst(obj, a) result(dtauda_w)
    ! d tau/d a for a model with constant dark energy equation of state, and all other
    ! densities as in the actual model. Agrees with dtauda to the last bit when w is the
    ! actual (constant) equation of state.
    class(*) :: obj
    real(dl), intent(in) :: a
    real(dl) :: dtauda_w, grhoa2

    select type (this => obj)
    type is (TEquivalent_wCDM)
        ! 8*pi*G*rho*a**4, using 8*pi*G*rho_de*a**4 = grhov*a**(1-3w)
        grhoa2 = this%State%grho_no_de(a) + this%State%grhov*a**(1 - 3*this%w)
        if (grhoa2 <= 0) then
            call GlobalError('halofit_casarini: trial model stops expanding before today', &
                error_unsupported_params)
            dtauda_w = 0
        else
            dtauda_w = sqrt(3/grhoa2)
        end if
    class default
        error stop 'dtauda_wconst: expected TEquivalent_wCDM'
    end select

    end function dtauda_wconst

    end module NonLinear
