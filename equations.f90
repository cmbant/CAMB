    ! Equations module for dark energy with constant equation of state parameter w
    ! allowing for perturbations based on a quintessence model
    ! by Antony Lewis (http://cosmologist.info/)

    ! Dec 2003, fixed (fatal) bug in tensor neutrino setup
    ! Changes to tight coupling approximation
    ! June 2004, fixed problem with large scale polarized tensors; support for vector modes
    ! Generate vector modes on their own. The power spectrum is taken from the scalar parameters.
    ! August 2004, fixed reionization term in lensing potential
    ! Nov 2004, change massive neutrino l_max to be consistent with massless if light
    ! Apr 2005, added DoLateRadTruncation option
    ! June 2006, added support for arbitary neutrino mass splittings
    ! Nov 2006, tweak to high_precision transfer function accuracy at lowish k
    ! June 2011, improved radiation approximations from arXiv: 1104.2933; Some 2nd order tight coupling terms
    !            merged fderivs and derivs so flat and non-flat use same equations; more precomputed arrays
    !            optimized neutrino sampling, and reorganised neutrino integration functions
    ! Feb 2013: fixed various issues with accuracy at larger neutrino masses
    ! Mar 2014: fixes for tensors with massive neutrinos


    !Return OmegaK - modify this if you add extra fluid components
    function GetOmegak()
    use precision
    use ModelParams
    real(dl)  GetOmegak

    GetOmegak = 1 - (CP%omegab+CP%omegac+CP%omegav+CP%omegan)

    end function GetOmegak


    subroutine Init_Backgrounds
    use ModelParams
    use DarkEnergyInterface
    !This is only called once per model, and is a good point to do any extra initialization.
    !It is called before first call to dtauda, but after
    !massive neutrinos are initialized and after GetOmegak

    call CP%DarkEnergy%Init_Background()

    end  subroutine Init_Backgrounds


    ! Background evolution
    function dtauda(a)
    use precision
    use ModelParams
    use MassiveNu
    use DarkEnergyInterface
    implicit none
    real(dl), intent(in) :: a
    real(dl) :: dtauda, rhonu, grhoa2, a2, grhov_t
    integer :: nu_i

    a2 = a ** 2
    call CP%DarkEnergy%BackgroundDensityAndPressure(a, grhov_t)

    !  8*pi*G*rho*a**4.
    grhoa2 = grhok * a2 + (grhoc + grhob) * a + grhog + grhornomass + &
        grhov_t * a2

    if (CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a * nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * grhormass(nu_i)
        end do
    end if

    dtauda = sqrt(3 / grhoa2)

    end function dtauda



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !Gauge-dependent perturbation equations

    module GaugeInterface
    use precision
    use ModelParams
    use MassiveNu
    use DarkEnergyInterface
    use Errors
    use Transfer
    implicit none
    public

    !Description of this file. Change if you make modifications.
    !Initialized in DarkEnergy%Init_Background
    character(LEN=:), allocatable :: Eqns_name

    logical, parameter :: plot_evolve = .false. !for outputing time evolution

    integer, parameter :: basic_num_eqns = 5

    logical :: DoTensorNeutrinos = .true.

    logical :: Evolve_baryon_cs = .false.
    !if true, evolves equation for Delta_{T_m} to get cs_2 = \delta p /\delta\rho for perfect gas

    logical :: Evolve_delta_xe = .false.

    logical :: Evolve_delta_Ts =.false. !Equilibrium result agree to sub-percent level

    logical :: DoLateRadTruncation = .true.
    !if true, use smooth approx to radition perturbations after decoupling on
    !small scales, saving evolution of irrelevant osciallatory multipole equations

    logical, parameter :: second_order_tightcoupling = .true.

    real(dl) :: Magnetic = 0._dl
    !Vector mode anisotropic stress in units of rho_gamma
    real(dl) :: vec_sig0 = 1._dl
    !Vector mode shear
    integer, parameter :: max_l_evolve = 256 !Maximum l we are ever likely to propagate
    !Note higher values increase size of Evolution vars, hence memory

    !Supported scalar initial condition flags
    integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
        initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0
    integer, parameter :: initial_nummodes =  initial_iso_neutrino_vel

    type EvolutionVars
        real(dl) q, q2
        real(dl) k_buf,k2_buf ! set in initial

        integer w_ix !Index of two quintessence equations
        integer Tg_ix !index of matter temerature perturbation
        integer reion_line_ix !index of matter temerature perturbation

        integer xe_ix !index of x_e perturbation
        integer Ts_ix !index of Delta_{T_s}

        integer r_ix !Index of the massless neutrino hierarchy
        integer g_ix !Index of the photon neutrino hierarchy

        integer q_ix !index into q_evolve array that gives the value q
        logical TransferOnly

        !       nvar  - number of scalar (tensor) equations for this k
        integer nvar,nvart, nvarv

        !Max_l for the various hierarchies
        integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
        integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt, MaxlNeededt
        logical EvolveTensorMassiveNu(max_nu)
        integer lmaxnrv, lmaxv, lmaxpolv
        integer lmaxline !21cm multipoles for getting reionization effect

        integer polind  !index into scalar array of polarization hierarchy

        !array indices for massive neutrino equations
        integer nu_ix(max_nu), nu_pert_ix
        integer nq(max_nu), lmaxnu_pert
        logical has_nu_relativistic

        !Initial values for massive neutrino v*3 variables calculated when switching
        !to non-relativistic approx
        real(dl) G11(max_nu),G30(max_nu)
        !True when using non-relativistic approximation
        logical MassiveNuApprox(max_nu)
        real(dl) MassiveNuApproxTime(max_nu)

        !True when truncating at l=2,3 when k*tau>>1 (see arXiv:1104.2933)
        logical high_ktau_neutrino_approx

        !Massive neutrino scheme being used at the moment
        integer NuMethod

        !True when using tight-coupling approximation (required for stability at early times)
        logical TightCoupling, TensTightCoupling
        real(dl) TightSwitchoffTime

        !Numer of scalar equations we are propagating
        integer ScalEqsToPropagate
        integer TensEqsToPropagate
        !beta > l for closed models
        integer FirstZerolForBeta
        !Tensor vars
        real(dl) aux_buf

        real(dl) pig, pigdot !For tight coupling
        real(dl) poltruncfac

        logical no_nu_multpoles, no_phot_multpoles
        integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
        logical nu_nonrelativistic(max_nu)

        real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
        real(dl) Kf(max_l_evolve)

        integer E_ix, B_ix !tensor polarization indices
        real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)

        logical :: saha !still high x_e
        logical :: evolve_TM !\delta T_g evolved separately

        real, pointer :: OutputTransfer(:) => null()
    end type EvolutionVars

    !precalculated arrays
    real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

    real(dl), parameter :: ep0=1.0d-2
    integer, parameter :: lmaxnu_high_ktau=4 !Jan2015, increased from 3 to fix mpk for mnu~6eV

    real(dl) epsw
    real(dl) nu_tau_notmassless(nqmax0+1,max_nu), nu_tau_nonrelativistic(max_nu),nu_tau_massive(max_nu)

    real(dl), private, external :: dtauda
    contains


    subroutine GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
    type(EvolutionVars) EV
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
    integer ind

    call dverk(EV,EV%ScalEqsToPropagate,derivs,tau,y,tauend,tol1,ind,c,EV%nvar,w)
    if (ind==-3) then
        call GlobalError('Dverk error -3: the subroutine was unable  to  satisfy  the  error ' &
            //'requirement  with a particular step-size that is less than or * ' &
            //'equal to hmin, which may mean that tol is too small' &
            //'--- but most likely you''ve messed up the y array indexing; ' &
            //'compiling with bounds checking may (or may not) help find the problem.',error_evolution)
    end if
    end subroutine GaugeInterface_ScalEv

    function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
        next_nq=1
    else
        q = int(nu_q(nq))
        if (q>=10) then
            next_nq = nqmax
        else
            next_nq = nq+1
        end if
    end if

    end function next_nu_nq

    recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    use RECDATA, only : CB1
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity, dopacity
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) noSwitch, smallTime
    !Sources
    real(dl) tau_switch_saha, Delta_TM, xe,a,tau_switch_evolve_TM

    noSwitch= CP%tau0+1
    smallTime =  min(tau, 1/EV%k_buf)/100

    tau_switch_ktau = noSwitch
    tau_switch_no_nu_multpoles= noSwitch
    tau_switch_no_phot_multpoles= noSwitch

    !Massive neutrino switches
    tau_switch_nu_massless = noSwitch
    tau_switch_nu_nonrel = noSwitch
    tau_switch_nu_massive= noSwitch

    !Sources
    tau_switch_saha=noSwitch
    if (Evolve_delta_xe .and. EV%saha)  tau_switch_saha = recombination_saha_tau
    tau_switch_evolve_TM=noSwitch
    if (Evolve_baryon_cs .and. .not. EV%Evolve_tm) tau_switch_evolve_TM = recombination_Tgas_tau

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
        tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%nq(nu_i) /= nqmax) then
                tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
            else if (.not. EV%nu_nonrelativistic(nu_i)) then
                tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
            else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
                tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
            end if
        end do
    end if

    if (DoLateRadTruncation) then
        if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
            tau_switch_no_nu_multpoles=max(15/EV%k_buf*AccuracyBoost,min(taurend,matter_verydom_tau))

        if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*AccuracyBoost)) &
            tau_switch_no_phot_multpoles =max(15/EV%k_buf,taurend)*AccuracyBoost
    end if

    next_switch = min(tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, tau_switch_nu_massive, &
        tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, tau_switch_nu_nonrel, noSwitch, &
        tau_switch_saha, tau_switch_evolve_TM)

    if (next_switch < tauend) then
        if (next_switch > tau+smallTime) then
            call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
            if (global_error_flag/=0) return
        end if

        EVout=EV

        if (next_switch == EV%TightSwitchoffTime) then
            !TightCoupling
            EVout%TightCoupling=.false.
            EVout%TightSwitchoffTime = noSwitch
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1
            !Set up variables with their tight coupling values
            y(EV%g_ix+2) = EV%pig
            call thermo(tau,cs2,opacity,dopacity)

            if (second_order_tightcoupling) then
                ! Francis-Yan Cyr-Racine November 2010

                y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity)*(1._dl+dopacity/opacity**2) + &
                    (3._dl/7._dl)*EV%pigdot*(EV%k_buf/opacity**2)*(-1._dl)

                y(EV%polind+2) = EV%pig/4 + EV%pigdot*(1._dl/opacity)*(-5._dl/8._dl- &
                    (25._dl/16._dl)*dopacity/opacity**2) + &
                    EV%pig*(EV%k_buf/opacity)**2*(-5._dl/56._dl)
                y(EV%polind+3) = (3._dl/7._dl)*(EV%k_buf/opacity)*y(EV%polind+2)*(1._dl + &
                    dopacity/opacity**2) + (3._dl/7._dl)*(EV%k_buf/opacity**2)*((EV%pigdot/4._dl)* &
                    (1._dl+(5._dl/2._dl)*dopacity/opacity**2))*(-1._dl)
            else
                y(EV%g_ix+3) = 3./7*y(EV%g_ix+2)*EV%k_buf/opacity
                y(EV%polind+2) = EV%pig/4
                y(EV%polind+3) =y(EV%g_ix+3)/4
            end if
        else if (next_switch==tau_switch_ktau) then
            !k tau >> 1, evolve massless neutrino effective fluid up to l=2
            EVout%high_ktau_neutrino_approx=.true.
            EV%nq(1:CP%Nu_mass_eigenstates) = nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch == tau_switch_nu_massless) then
            !Mass starts to become important, start evolving next momentum mode
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (EV%nq(nu_i) /= nqmax .and. &
                        next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                    EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                    call SetupScalarArrayIndices(EVout)
                    call CopyScalarVariableArray(y,yout, EV, EVout)
                    EV=EVout
                    y=yout
                    exit
                end if
            end do
        else if (next_switch == tau_switch_nu_nonrel) then
            !Neutrino becomes non-relativistic, don't need high L
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%nu_nonrelativistic(nu_i) .and.  next_switch==nu_tau_nonrelativistic(nu_i) ) then
                    EVout%nu_nonrelativistic(nu_i)=.true.
                    call SetupScalarArrayIndices(EVout)
                    call CopyScalarVariableArray(y,yout, EV, EVout)
                    EV=EVout
                    y=yout
                    exit
                end if
            end do
        else if (next_switch == tau_switch_nu_massive) then
            !Very non-relativistic neutrinos, switch to truncated velocity-weight hierarchy
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%MassiveNuApprox(nu_i) .and.  next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                    call SwitchToMassiveNuApprox(EV,y, nu_i)
                    exit
                end if
            end do
        else if (next_switch==tau_switch_no_nu_multpoles) then
            !Turn off neutrino hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_nu_multpoles=.true.
            EVOut%nq(1:CP%Nu_mass_eigenstates ) = nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_no_phot_multpoles) then
            !Turn off photon hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_phot_multpoles=.true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_saha) then
            !Sources
            ind=1
            EVout%saha = .false.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
            a=y(1)
            Delta_Tm = y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
            xe= Recombination_xe(a)
            y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(4) + (3./2+  CB1/(CP%TCMB/a))*Delta_TM)
        else if (next_switch==tau_switch_evolve_TM) then
            !Sources
            ind=1
            EVout%evolve_TM = .true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
            y(EV%Tg_ix) =y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
        end if

        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
        return
    end if

    call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)

    end subroutine GaugeInterface_EvolveScal

    subroutine GaugeInterface_EvolveTens(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    type(EvolutionVars) EV, EVOut
    real(dl) c(24),w(EV%nvart,9), y(EV%nvart),yout(EV%nvart), tol1, tau, tauend
    integer ind
    real(dl) opacity, cs2

    if (EV%TensTightCoupling .and. tauend > EV%TightSwitchoffTime) then
        if (EV%TightSwitchoffTime > tau) then
            call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,EV%TightSwitchoffTime,tol1,ind,c,EV%nvart,w)
        end if
        EVOut=EV
        EVOut%TensTightCoupling = .false.
        call SetupTensorArrayIndices(EVout)
        call CopyTensorVariableArray(y,yout,Ev, Evout)
        Ev = EvOut
        y=yout
        call thermo(tau,cs2,opacity)
        y(EV%g_ix+2)= 32._dl/45._dl*EV%k_buf/opacity*y(3)
        y(EV%E_ix+2) = y(EV%g_ix+2)/4
    end if

    call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,tauend,tol1,ind,c,EV%nvart,w)


    end subroutine GaugeInterface_EvolveTens

    function DeltaTimeMaxed(a1,a2, tol) result(t)
    real(dl) a1,a2,t
    real(dl), optional :: tol
    if (a1>1._dl) then
        t=0
    elseif (a2 > 1._dl) then
        t = DeltaTime(a1,1.01_dl, tol)
    else
        t= DeltaTime(a1,a2, tol)
    end if
    end function DeltaTimeMaxed

    subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/CP%tau0

    if (CP%WantScalars) then
        do j=2,max_l_evolve
            polfac(j)=real((j+3)*(j-1),dl)/(j+1)
        end do
    end if

    if (CP%WantVectors) then
        do j=2,max_l_evolve
            vecfac(j)=real((j+2),dl)/(j+1)
            vecfacpol(j)=real((j+3)*j,dl)*(j-1)*vecfac(j)/(j+1)**2
        end do
    end if

    do j=1,max_l_evolve
        denl(j)=1._dl/(2*j+1)
    end do

    do nu_i=1, CP%Nu_Mass_eigenstates
        nu_mass = max(0.1_dl,nu_masses(nu_i))
        a_mass =  1.e-1_dl/nu_mass/lAccuracyBoost
        !if (HighAccuracyDefault) a_mass=a_mass/4
        time=DeltaTime(0._dl,nu_q(1)*a_mass)
        nu_tau_notmassless(1, nu_i) = time
        do j=2,nqmax
            !times when each momentum mode becomes signficantly nonrelativistic
            time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
            nu_tau_notmassless(j, nu_i) = time
        end do

        a_nonrel =  2.5d0/nu_mass*AccuracyBoost !!!Feb13tweak
        nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
        a_massive =  17.d0/nu_mass*AccuracyBoost
        nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
    end do

    end subroutine GaugeInterface_Init


    subroutine SetupScalarArrayIndices(EV, max_num_eqns)
    !Set up array indices after the lmax have been decided
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    integer, intent(out), optional :: max_num_eqns
    integer neq, maxeq, nu_i

    neq=basic_num_eqns
    maxeq=neq
    if (.not. EV%no_phot_multpoles) then
        !Photon multipoles
        EV%g_ix=basic_num_eqns+1
        if (EV%TightCoupling) then
            neq=neq+2
        else
            neq = neq+ (EV%lmaxg+1)
            !Polarization multipoles
            EV%polind = neq -1 !polind+2 is L=2, for polarizationthe first calculated
            neq=neq + EV%lmaxgpol-1
        end if
    end if
    if (.not. EV%no_nu_multpoles) then
        !Massless neutrino multipoles
        EV%r_ix= neq+1
        if (EV%high_ktau_neutrino_approx) then
            neq=neq + 3
        else
            neq=neq + (EV%lmaxnr+1)
        end if
    end if
    maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

    !Dark energy
    if (.not. CP%DarkEnergy%is_cosmological_constant) then
        EV%w_ix = neq + 1
        neq = neq + CP%DarkEnergy%num_perturb_equations
        maxeq = maxeq + CP%DarkEnergy%num_perturb_equations
    else
        EV%w_ix = 0
    end if

    !Sources
    if (Evolve_delta_xe) then
        if (.not. EV%saha) then
            EV%xe_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
    end if

    if (Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            EV%Tg_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
        if (Do21cm .and. line_reionization) then
            EV%reion_line_ix = neq+1
            neq=neq+ EV%lmaxline+1 +  EV%lmaxline-1
            maxeq=maxeq+EV%lmaxline+1 +  EV%lmaxline-1
        end if
    end if

    if (Evolve_delta_Ts) then
        EV%Ts_ix = neq+1
        neq=neq+1
        maxeq=maxeq+1
    end if

    !Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
        EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=nqmax)
        if (EV%has_nu_relativistic) then
            EV%lmaxnu_pert=EV%lmaxnu
            EV%nu_pert_ix=neq+1
            neq = neq+ EV%lmaxnu_pert+1
            maxeq=maxeq+ EV%lmaxnu_pert+1
        else
            EV%lmaxnu_pert=0
        end if

        do nu_i=1, CP%Nu_Mass_eigenstates
            if (EV%high_ktau_neutrino_approx) then
                EV%lmaxnu_tau(nu_i) = int(lmaxnu_high_ktau *lAccuracyBoost)
            else
                EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i)*lAccuracyBoost),EV%lmaxnu),3)
                !!!Feb13tweak
                if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu_tau(nu_i),nint(4*lAccuracyBoost))
            end if
            if (nu_masses(nu_i) > 5000 .and. CP%Transfer%high_precision) EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i)*2 !megadamping
            EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

            EV%nu_ix(nu_i)=neq+1
            if (EV%MassiveNuApprox(nu_i)) then
                neq = neq+4
            else
                neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
            endif
            maxeq = maxeq + nqmax*(EV%lmaxnu+1)
        end do
    else
        EV%has_nu_relativistic = .false.
    end if

    EV%ScalEqsToPropagate = neq
    if (present(max_num_eqns)) then
        max_num_eqns=maxeq
    end if

    end subroutine SetupScalarArrayIndices

    subroutine CopyScalarVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvar)
    real(dl), intent(out) :: yout(EVout%nvar)
    integer lmax,i, nq
    integer nnueq,nu_i, ix_off, ix_off2, ind, ind2
    real(dl) q, pert_scale

    yout=0
    yout(1:basic_num_eqns) = y(1:basic_num_eqns)

    ! DarkEnergy
    if (CP%DarkEnergy%num_perturb_equations > 0) &
        yout(EVOut%w_ix:EVOut%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
            y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1)

    if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
        if (EV%TightCoupling .or. EVOut%TightCoupling) then
            lmax=1
        else
            lmax = min(EV%lmaxg,EVout%lmaxg)
        end if
        yout(EVout%g_ix:EVout%g_ix+lmax)=y(EV%g_ix:EV%g_ix+lmax)
        if (.not. EV%TightCoupling .and. .not. EVOut%TightCoupling) then
            lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
            yout(EVout%polind+2:EVout%polind+lmax)=y(EV%polind+2:EV%polind+lmax)
        end if
    end if

    if (.not. EV%no_nu_multpoles .and. .not. EVout%no_nu_multpoles) then
        if (EV%high_ktau_neutrino_approx .or. EVout%high_ktau_neutrino_approx) then
            lmax=2
        else
            lmax = min(EV%lmaxnr,EVout%lmaxnr)
        end if
        yout(EVout%r_ix:EVout%r_ix+lmax)=y(EV%r_ix:EV%r_ix+lmax)
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            ix_off=EV%nu_ix(nu_i)
            ix_off2=EVOut%nu_ix(nu_i)
            if (EV%MassiveNuApprox(nu_i) .and. EVout%MassiveNuApprox(nu_i)) then
                nnueq=4
                yout(ix_off2:ix_off2+nnueq-1)=y(ix_off:ix_off+nnueq-1)
            else if (.not. EV%MassiveNuApprox(nu_i) .and. .not. EVout%MassiveNuApprox(nu_i)) then
                lmax=min(EV%lmaxnu_tau(nu_i),EVOut%lmaxnu_tau(nu_i))
                nq = min(EV%nq(nu_i), EVOut%nq(nu_i))
                do i=1,nq
                    ind= ix_off + (i-1)*(EV%lmaxnu_tau(nu_i)+1)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(ind:ind+lmax)
                end do
                do i=nq+1, EVOut%nq(nu_i)
                    lmax = min(EVOut%lmaxnu_tau(nu_i), EV%lmaxnr)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax)

                    !Add leading correction for the mass
                    q=nu_q(i)
                    pert_scale=(nu_masses(nu_i)/q)**2/2
                    lmax = min(lmax,EV%lmaxnu_pert)
                    yout(ind2:ind2+lmax) = yout(ind2:ind2+lmax) &
                        + y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)*pert_scale
                end do
            end if
        end do

        if (EVOut%has_nu_relativistic .and. EV%has_nu_relativistic) then
            lmax = min(EVOut%lmaxnu_pert, EV%lmaxnu_pert)
            yout(EVout%nu_pert_ix:EVout%nu_pert_ix+lmax)=  y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)
        end if
    end if
    !Sources
    if (.not. EV%saha .and. .not. EVOut%saha) then
        yout(EVOut%xe_ix) =y(EV%xe_ix)
    end if
    if (Evolve_baryon_cs) then
        if (EV%Evolve_TM .and. EVout%Evolve_TM) yout(EVOut%Tg_ix) = y(EV%Tg_ix)
        if (Do21cm .and. line_reionization) then
            yout(EVOut%reion_line_ix:EVOut%reion_line_ix+EVout%lmaxline +  EVout%lmaxline-1) = &
                y(EV%reion_line_ix:EV%reion_line_ix+EV%lmaxline +  EV%lmaxline-1)
        end if
    end if
    if (Evolve_delta_Ts) then
        yout(EVOut%Ts_ix) = y(EV%Ts_ix)
    end if

    end subroutine CopyScalarVariableArray


    subroutine SetupTensorArrayIndices(EV, maxeq)
    type(EvolutionVars) EV
    integer nu_i, neq
    integer, optional, intent(out) :: maxeq
    neq=3
    EV%g_ix = neq-1 !EV%g_ix+2 is quadrupole
    if (.not. EV%TensTightCoupling) then
        EV%E_ix = EV%g_ix + (EV%lmaxt-1)
        EV%B_ix = EV%E_ix + (EV%lmaxpolt-1)
        neq = neq+ (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    if (present(maxeq)) then
        maxeq =3 + (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    EV%r_ix = neq -1
    if (DoTensorNeutrinos) then
        neq = neq + EV%lmaxnrt-1
        if (present(maxeq)) maxeq = maxeq+EV%lmaxnrt-1
        if (CP%Num_Nu_massive /= 0 ) then
            do nu_i=1, CP%nu_mass_eigenstates
                EV%EvolveTensorMassiveNu(nu_i) = nu_tau_nonrelativistic(nu_i) < 0.8*tau_maxvis*AccuracyBoost
                if (EV%EvolveTensorMassiveNu(nu_i)) then
                    EV%nu_ix(nu_i)=neq-1
                    neq = neq+ nqmax*(EV%lmaxnut-1)
                    if (present(maxeq)) maxeq = maxeq + nqmax*(EV%lmaxnut-1)
                end if
            end do
        end if
    end if

    EV%TensEqsToPropagate = neq

    end  subroutine SetupTensorArrayIndices

    subroutine CopyTensorVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvart)
    real(dl), intent(out) :: yout(EVout%nvart)
    integer lmaxpolt, lmaxt, nu_i, ind, ind2, i

    yout=0
    yout(1:3) = y(1:3)
    if (.not. EVOut%TensTightCoupling .and. .not.EV%TensTightCoupling) then
        lmaxt = min(EVOut%lmaxt,EV%lmaxt)
        yout(EVout%g_ix+2:EVout%g_ix+lmaxt)=y(EV%g_ix+2:EV%g_ix+lmaxt)
        lmaxpolt = min(EV%lmaxpolt, EVOut%lmaxpolt)
        yout(EVout%E_ix+2:EVout%E_ix+lmaxpolt)=y(EV%E_ix+2:EV%E_ix+lmaxpolt)
        yout(EVout%B_ix+2:EVout%B_ix+lmaxpolt)=y(EV%B_ix+2:EV%B_ix+lmaxpolt)
    end if
    if (DoTensorNeutrinos) then
        lmaxt=min(EV%lmaxnrt,EVOut%lmaxnrt)
        yout(EVout%r_ix+2:EVout%r_ix+lmaxt)=y(EV%r_ix+2:EV%r_ix+lmaxt)
        do nu_i =1, CP%nu_mass_eigenstates
            if (EV%EvolveTensorMassiveNu(nu_i)) then
                lmaxt=min(EV%lmaxnut,EVOut%lmaxnut)
                do i=1,nqmax
                    ind= EV%nu_ix(nu_i) + (i-1)*(EV%lmaxnut-1)
                    ind2=EVOut%nu_ix(nu_i)+ (i-1)*(EVOut%lmaxnut-1)
                    yout(ind2+2:ind2+lmaxt) = y(ind+2:ind+lmaxt)
                end do
            end if
        end do
    end if

    end subroutine CopyTensorVariableArray

    subroutine GetNumEqns(EV)
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    real(dl) scal, max_nu_mass
    integer nu_i,q_rel,j

    if (CP%Num_Nu_massive == 0) then
        EV%lmaxnu=0
        max_nu_mass=0
    else
        max_nu_mass = maxval(nu_masses(1:CP%Nu_mass_eigenstates))
        do nu_i = 1, CP%Nu_mass_eigenstates
            !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
            q_rel=0
            do j=1, nqmax
                !two different q's here EV%q ~k
                if (nu_q(j) > nu_masses(nu_i)*adotrad/EV%q) exit
                q_rel = q_rel + 1
            end do

            if (q_rel>= nqmax-2 .or. CP%WantTensors) then
                EV%nq(nu_i)=nqmax
            else
                EV%nq(nu_i)=q_rel
            end if
            !q_rel = nint(nu_masses(nu_i)*adotrad/EV%q) !two dffierent q's here EV%q ~k
            !EV%nq(nu_i)=max(0,min(nqmax0,q_rel)) !number of momentum modes to evolve intitially
            EV%nu_nonrelativistic(nu_i) = .false.
        end do

        EV%NuMethod = CP%MassiveNuMethod
        if (EV%NuMethod == Nu_Best) EV%NuMethod = Nu_Trunc
        !l_max for massive neutrinos
        if (CP%Transfer%high_precision) then
            EV%lmaxnu=nint(25*lAccuracyBoost)
        else
            EV%lmaxnu=max(3,nint(10*lAccuracyBoost))
            if (max_nu_mass>700) EV%lmaxnu=max(3,nint(15*lAccuracyBoost)) !Feb13 tweak
        endif
    end if

    if (CP%closed) then
        EV%FirstZerolForBeta = nint(EV%q*CP%r)
    else
        EV%FirstZerolForBeta=l0max !a large number
    end if

    EV%high_ktau_neutrino_approx = .false.
    if (CP%WantScalars) then
        EV%TightCoupling=.true.
        EV%no_phot_multpoles =.false.
        EV%no_nu_multpoles =.false.
        EV%MassiveNuApprox=.false.
        !Sources
        EV%saha = .true.
        EV%Evolve_TM = .false.

        if (HighAccuracyDefault .and. CP%AccuratePolarization) then
            EV%lmaxg  = max(nint(11*lAccuracyBoost),3)
        else
            EV%lmaxg  = max(nint(8*lAccuracyBoost),3)
        end if
        EV%lmaxnr = max(nint(14*lAccuracyBoost),3)
        if (max_nu_mass>700 .and. HighAccuracyDefault) EV%lmaxnr = max(nint(32*lAccuracyBoost),3) !Feb13 tweak

        EV%lmaxgpol = EV%lmaxg
        if (.not.CP%AccuratePolarization) EV%lmaxgpol=max(nint(4*lAccuracyBoost),3)

        if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            if (CP%AccuratePolarization) scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*lAccuracyBoost))
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*lAccuracyBoost))
            !Sources
            if (line_phot_quadrupole) then
                EV%lmaxg=EV%lmaxg*8
                EV%lmaxgpol=EV%lmaxgpol*4
            elseif (CP%AccurateReionization) then
                EV%lmaxg=EV%lmaxg*4
                EV%lmaxgpol=EV%lmaxgpol*2
            end if
        end if

        if (EV%TransferOnly) then
            EV%lmaxgpol = min(EV%lmaxgpol,nint(5*lAccuracyBoost))
            EV%lmaxg = min(EV%lmaxg,nint(6*lAccuracyBoost))
        end if
        if (CP%Transfer%high_precision .or. Do21cm) then
            if (HighAccuracyDefault) then
                EV%lmaxnr=max(nint(45*lAccuracyBoost),3)
            else
                EV%lmaxnr=max(nint(30*lAccuracyBoost),3)
            endif
            if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
                EV%lmaxg=max(EV%lmaxg,10)
            end if
        end if

        if (Do21cm .and. line_reionization) then
            EV%lmaxg =  EV%lmaxg*8
            EV%lmaxgpol = EV%lmaxgpol*3
        end if

        if (Do21cm .or.Evolve_delta_xe .or. Evolve_delta_Ts) Evolve_baryon_cs = .true.

        if (Do21cm .and. line_reionization) then
            EV%lmaxline  = EV%lmaxg
        end if

        if (CP%closed) then
            EV%lmaxnu=min(EV%lmaxnu, EV%FirstZerolForBeta-1)
            EV%lmaxnr=min(EV%lmaxnr, EV%FirstZerolForBeta-1)
            EV%lmaxg=min(EV%lmaxg, EV%FirstZerolForBeta-1)
            EV%lmaxgpol=min(EV%lmaxgpol, EV%FirstZerolForBeta-1)
        end if

        EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
        EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
        if (EV%MaxlNeeded > max_l_evolve) stop 'Need to increase max_l_evolve'
        call SetupScalarArrayIndices(EV,EV%nvar)
        if (CP%closed) EV%nvar=EV%nvar+1 !so can reference lmax+1 with zero coefficient
        EV%lmaxt=0
    else
        EV%nvar=0
    end if

    if (CP%WantTensors) then
        EV%TensTightCoupling = .true.
        EV%lmaxt=max(3,nint(8*lAccuracyBoost))
        EV%lmaxpolt = max(3,nint(4*lAccuracyBoost))
        ! if (EV%q < 1e-3) EV%lmaxpolt=EV%lmaxpolt+1
        if (DoTensorNeutrinos) then
            EV%lmaxnrt=nint(6*lAccuracyBoost)
            EV%lmaxnut=EV%lmaxnrt
        else
            EV%lmaxnut=0
            EV%lmaxnrt=0
        end if
        if (CP%closed) then
            EV%lmaxt=min(EV%FirstZerolForBeta-1,EV%lmaxt)
            EV%lmaxpolt=min(EV%FirstZerolForBeta-1,EV%lmaxpolt)
            EV%lmaxnrt=min(EV%FirstZerolForBeta-1,EV%lmaxnrt)
            EV%lmaxnut=min(EV%FirstZerolForBeta-1,EV%lmaxnut)
        end if
        EV%MaxlNeededt=max(EV%lmaxpolt,EV%lmaxt, EV%lmaxnrt, EV%lmaxnut)
        if (EV%MaxlNeededt > max_l_evolve) stop 'Need to increase max_l_evolve'
        call SetupTensorArrayIndices(EV, EV%nvart)
    else
        EV%nvart=0
    end if


    if (CP%WantVectors) then
        EV%lmaxv=max(10,nint(8*lAccuracyBoost))
        EV%lmaxpolv = max(5,nint(5*lAccuracyBoost))

        EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3

        EV%lmaxnrv=nint(30*lAccuracyBoost)

        EV%nvarv=EV%nvarv+EV%lmaxnrv
        if (CP%Num_Nu_massive /= 0 ) then
            stop 'massive neutrinos not supported for vector modes'
        end if
    else
        EV%nvarv=0
    end if

    end subroutine GetNumEqns

    !cccccccccccccccccccccccccccccccccc
    subroutine SwitchToMassiveNuApprox(EV,y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i

    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a=y(1)
    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

    !Integrate over q
    call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
    !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
    dpnu=dpnu/rhonu
    qnu=qnu/rhonu
    clxnu = clxnu/rhonu
    pinu=pinu/rhonu

    yout(EVout%nu_ix(nu_i))=clxnu
    yout(EVout%nu_ix(nu_i)+1)=dpnu
    yout(EVout%nu_ix(nu_i)+2)=qnu
    yout(EVout%nu_ix(nu_i)+3)=pinu

    call Nu_Intvsq(EV,y, a, nu_i, EVout%G11(nu_i),EVout%G30(nu_i))
    !Analytic solution for higher moments, proportional to a^{-3}
    EVout%G11(nu_i)=EVout%G11(nu_i)*a2*a/rhonu
    EVout%G30(nu_i)=EVout%G30(nu_i)*a2*a/rhonu

    EV=EVout
    y=yout

    end subroutine SwitchToMassiveNuApprox

    subroutine MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum,clxnu_all)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a
    real(dl), optional :: grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum,clxnu_all
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) adotoa, grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu, grhonu, dgrhonu

    grhonu=0
    dgrhonu=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            !dpnu = y(EV%iq0+1+off_ix)
            qnu=y(EV%nu_ix(nu_i)+2)
            pinu=y(EV%nu_ix(nu_i)+3)
            pinudot=yprime(EV%nu_ix(nu_i)+3)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !dpnu=dpnu/rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
            pinu=pinu/rhonu
            adotoa = 1/(a*dtauda(a))
            rhonudot = Nu_drho(a*nu_masses(nu_i),adotoa,rhonu)

            call Nu_pinudot(EV,y, yprime, a,adotoa, nu_i,pinudot)
            pinudot=pinudot/rhonu - rhonudot/rhonu*pinu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grhonu = grhonu  + grhonu_t
        if (present(gpres)) gpres= gpres + gpnu_t

        dgrhonu= dgrhonu + grhonu_t*clxnu
        if (present(dgq)) dgq  = dgq   + grhonu_t*qnu
        if (present(dgpi)) dgpi = dgpi  + grhonu_t*pinu
        if (present(gdpi_diff)) gdpi_diff = gdpi_diff + pinu*(3*gpnu_t-grhonu_t)
        if (present(pidot_sum)) pidot_sum = pidot_sum + grhonu_t*pinudot
    end do
    if (present(grho)) grho = grho  + grhonu
    if (present(dgrho)) dgrho= dgrho + dgrhonu
    if (present(clxnu_all)) clxnu_all = dgrhonu/grhonu

    end subroutine MassiveNuVarsOut

    subroutine Nu_Integrate_L012(EV,y,a,nu_i,drhonu,fnu,dpnu,pinu)
    type(EvolutionVars) EV
    !  Compute the perturbations of density and energy flux
    !  of one eigenstate of massive neutrinos, in units of the mean
    !  density of one eigenstate of massless neutrinos, by integrating over
    !  momentum.
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  drhonu,fnu
    real(dl), optional, intent(OUT) :: dpnu,pinu
    real(dl) tmp, am, aq,v, pert_scale
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.

    drhonu=0
    fnu=0
    if (present(dpnu)) then
        dpnu=0
        pinu=0
    end if
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    do iq=1,EV%nq(nu_i)
        aq=am/nu_q(iq)
        v=1._dl/sqrt(1._dl+aq*aq)
        drhonu=drhonu+ nu_int_kernel(iq)* y(ind)/v
        fnu=fnu+nu_int_kernel(iq)* y(ind+1)
        if (present(dpnu)) then
            dpnu=dpnu+  nu_int_kernel(iq)* y(ind)*v
            pinu=pinu+ nu_int_kernel(iq)*y(ind+2)*v
        end if
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix
    do iq=EV%nq(nu_i)+1,nqmax
        !Get the rest from perturbatively relativistic expansion
        aq=am/nu_q(iq)
        v=1._dl/sqrt(1._dl+aq*aq)
        pert_scale=(nu_masses(nu_i)/nu_q(iq))**2/2
        tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind)  )
        drhonu=drhonu+ tmp/v
        fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
        if (present(dpnu)) then
            dpnu=dpnu+ tmp*v
            pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
        end if
    end do

    if (present(dpnu)) then
        dpnu = dpnu/3
    end if

    end subroutine Nu_Integrate_L012

    subroutine Nu_pinudot(EV,y, ydot, a,adotoa, nu_i,pinudot)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a,adotoa, y(EV%nvar), ydot(EV%nvar)

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) pinudot
    real(dl) aq,q,v,aqdot,vdot
    real(dl) psi2,psi2dot
    real(dl) am, pert_scale
    integer iq,ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    pinudot=0._dl
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        pinudot=pinudot+nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,nqmax
        q=nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        pert_scale=(nu_masses(nu_i)/q)**2/2
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
        psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
        pinudot=pinudot+nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
    end do

    end subroutine Nu_pinudot

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function Nu_pi(EV, y, a, nu_i) result(pinu)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvart)
    real(dl) :: am
    real(dl) pinu,q,aq,v
    integer iq, ind

    if (EV%nq(nu_i)/=nqmax) stop 'Nu_pi: nq/=nqmax'
    pinu=0
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1, EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        pinu=pinu+nu_int_kernel(iq)*y(ind)*v
        ind =ind+EV%lmaxnut+1
    end do

    end function Nu_pi

    !cccccccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_Intvsq(EV,y, a, nu_i, G11,G30)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  G11,G30

    !  Compute the third order variables (in velocity dispersion)
    !by integrating over momentum.
    real(dl) aq,q,v, am
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=nqmax) stop 'Nu_Intvsq nq/=nqmax0'
    do iq=1, EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        G11=G11+nu_int_kernel(iq)*y(ind+1)*v**2
        if (EV%lmaxnu_tau(nu_i)>2) then
            G30=G30+nu_int_kernel(iq)*y(ind+3)*v**2
        end if
        ind = ind+EV%lmaxnu_tau(nu_i)+1
    end do

    end subroutine Nu_Intvsq


    subroutine MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq, wnu_arr)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), a, grho,gpres,dgrho,dgq
    real(dl), intent(out), optional :: wnu_arr(max_nu)
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    integer nu_i
    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            qnu=y(EV%nu_ix(nu_i)+2)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grho = grho  + grhonu_t
        gpres= gpres + gpnu_t
        dgrho= dgrho + grhonu_t*clxnu
        dgq  = dgq   + grhonu_t*qnu

        if (present(wnu_arr)) then
            wnu_arr(nu_i) =pnu/rhonu
        end if
    end do

    end subroutine MassiveNuVars


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    function Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, vterm )
    !Delta_Tspin - Delta_TCMB
    use constants
    !vterm = hdot/clh + k*n/3/clh
    real(dl), intent(in) :: a,Delta_source,Delta_TCMB,Delta_Tm,Tmat,Trad,xe, Delta_xe, vterm
    real(dl) :: Get21cm_source2
    real(dl) Rgamma,Rm
    real(dl) dC10, n_H,C10, C10_HH, C10_eH
    real(dl) kappa_HH,kappa_eH
    real(dl) tau_eps
    real(dl) dtauda, H
    external dtauda
    real(dl) TSpin
    n_H = NNow/a**3
    kappa_HH = kappa_HH_21cm(Tmat, .false.)
    kappa_eH = kappa_eH_21cm(Tmat, .false.)
    C10_HH = n_H*kappa_HH* (1- xe)
    C10_eH = n_H*kappa_eH*xe
    C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
    Rgamma = 1._dl/(C10+A10*Trad/T_21cm)
    Rm = 1._dl/(C10+A10*Tmat/T_21cm)

    !          TSpin=TsRecfast(a)
    !          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
    !                  n_H*kappa_pH_21cm(Tmat, .false.)*xe
    !          if (a>0.5) stop

    dC10 = (C10*Delta_source + &
        (C10_HH*kappa_HH_21cm(Tmat, .true.)+C10_eH*kappa_eH_21cm(Tmat, .true.)) * &
        Delta_Tm + (kappa_eH-kappa_HH)*xe*n_H*Delta_xe)



    Get21cm_source2 =  dC10*(Rgamma-Rm) +  C10*(Rm*Delta_tm - Delta_TCMB*Rgamma)

    TSpin=Recombination_Ts(a)
    H = (1/(a*dtauda(a)))
    tau_eps = a*line21_const*NNow/a**3/H/Tspin/1000

    Get21cm_source2 = Get21cm_source2 + &
        tau_eps/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) ) * &
        (Delta_source -vterm + dC10/C10 + 2*( - Rgamma*dC10 + Delta_TCMB*(C10*Rgamma-1)) &
        + Trad/(Tmat-Trad)*(Delta_tm-Delta_TCMB)   )

    end function Get21cm_source2



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    function Get21cm_dTs(a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
    !d Delta T_s / d eta dropping small \Delta_xe terms
    use constants
    real(dl), intent(in) :: a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe
    real(dl) :: Get21cm_dTs
    real(dl) n_H,C10, C10_HH, C10_eH, delta_C10
    real(dl) kappa_HH,kappa_eH, TSpin

    n_H = NNow/a**3
    kappa_HH = kappa_HH_21cm(Tmat, .false.)
    kappa_eH = kappa_eH_21cm(Tmat, .false.)
    C10_HH = n_H*kappa_HH* (1- xe)
    C10_eH = n_H*kappa_eH*xe
    C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
    TSpin=Recombination_ts(a)
    delta_C10 = C10*Delta_n + (C10_HH*kappa_HH_21cm(Tmat, .true.)+C10_eH*kappa_eH_21cm(Tmat, .true.))*Delta_Tm

    !          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
    !                  n_H*kappa_pH_21cm(Tmat, .false.)*xe

    Get21cm_dTs =  4*a*( TSpin/TMat*(Delta_Tm-Delta_ts)*C10 + (1-TSpin/TMat)*delta_C10 + &
        (Trad*Delta_TCMB - Tspin*Delta_Ts)*A10/T_21cm ) * MPC_in_sec

    end function Get21cm_dTs



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine output(EV,y, j,tau,sources)
    use ThermoData
    use lvalues
    use ModelData
    use constants, only : barssc0
    use RedshiftSpaceData
    use Recombination

    implicit none
    integer j
    type(EvolutionVars) EV
    real(dl), target :: y(EV%nvar),yprime(EV%nvar)
    real(dl), dimension(:),pointer :: ypol,ypolprime

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
    real(dl) qgdot,pigdot,pirdot,vbdot,dgrho
    real(dl) a,a2,dz,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) w_dark_energy_t

    real(dl) tau,x,divfac
    real(dl) dgpi_diff, pidot_sum
    real(dl), target :: pol(3),polprime(3)
    !dgpi_diff = sum (3*p_nu -rho_nu)*pi_nu

    real(dl) k,k2  ,adotoa, grho, gpres,etak,phi,dgpi
    real(dl) diff_rhopi, octg, octgprime
    real(dl) sources(CTransScal%NumSources)
    real(dl) ISW
    !Sources
    real(dl) phidot
    real(dl) Tmat,Trad, Tspin, Delta_source, Delta_source2
    real(dl) Delta_TCMB, Delta_tm, Delta_xe
    real(dl) polter_line, chi
    integer w_ix, lineoff,lineoffpol
    real(dl) cs2, xe,opacity, delta_p
    real(dl) s(0:10), t(0:10)
    real(dl) counts_radial_source, counts_velocity_source, counts_density_source, counts_ISW_source, &
        counts_redshift_source, counts_timedelay_source, counts_potential_source
    real(dl) help

    sources = 0
    yprime = 0

    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)

    if (EV%TightCoupling .or. EV%no_phot_multpoles) then
        pol=0
        polprime=0
        ypolprime => polprime
        ypol => pol
    else
        ypolprime => yprime(EV%polind+1:)
        ypol => y(EV%polind+1:)
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    a   =y(1)
    a2  =a*a
    etak=y(2)
    clxc=y(3)
    clxb=y(4)
    vb  =y(5)
    vbdot =yprime(5)

    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    call CP%DarkEnergy%BackgroundDensityAndPressure(a, grhov_t, w_dark_energy_t)

    !  8*pi*a*a*SUM[rho_i*clx_i] add radiation later
    dgrho=grhob_t*clxb+grhoc_t*clxc

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    if (.not. CP%DarkEnergy%is_cosmological_constant) &
        call CP%DarkEnergy%AddStressEnergy(dgrho, dgq, grhov_t, y, EV%w_ix, .true.)

    gpres = (grhog_t + grhor_t) / 3 + w_dark_energy_t*grhov_t
    grho = grhob_t + grhoc_t + grhor_t + grhog_t + grhov_t

    dgpi = 0
    dgpi_diff = 0
    pidot_sum = 0

    if (CP%Num_Nu_Massive /= 0) then
        call MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum)
    end if

    adotoa=sqrt((grho+grhok)/3)

    if (EV%no_nu_multpoles) then
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxr=-4*dz/k
        qr=-4._dl/3*z
        pir=0
        pirdot=0
    else
        clxr=y(EV%r_ix)
        qr  =y(EV%r_ix+1)
        pir =y(EV%r_ix+2)
        pirdot=yprime(EV%r_ix+2)
    end if

    if (EV%no_phot_multpoles) then
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxg=-4*dz/k -4/k*opac(j)*(vb+z)
        qg=-4._dl/3*z
        pig=0
        pigdot=0
        octg=0
        octgprime=0
        qgdot = -4*dz/3
    else
        if (EV%TightCoupling) then
            pig = EV%pig
            pigdot=EV%pigdot
            if (second_order_tightcoupling) then
                octg = (3._dl/7._dl)*pig*(EV%k_buf/opac(j))
                ypol(2) = EV%pig/4 + pigdot*(1._dl/opac(j))*(-5._dl/8._dl)
                ypol(3) = (3._dl/7._dl)*(EV%k_buf/opac(j))*ypol(2)
            else
                ypol(2) = EV%pig/4
                octg=0
            end if
            octgprime=0
        else
            pig =y(EV%g_ix+2)
            pigdot=yprime(EV%g_ix+2)
            octg=y(EV%g_ix+3)
            octgprime=yprime(EV%g_ix+3)
        end if
        clxg=y(EV%g_ix)
        qg  =y(EV%g_ix+1)
        qgdot =yprime(EV%g_ix+1)
    end if

    dgrho = dgrho + grhog_t*clxg+grhor_t*clxr
    dgq   = dgq   + grhog_t*qg+grhor_t*qr
    dgpi  = dgpi  + grhor_t*pir + grhog_t*pig


    !  Get sigma (shear) and z from the constraints
    !  have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)

    polter = 0.1_dl*pig+9._dl/15._dl*ypol(2)

    if (CP%flat) then
        x=k*(CP%tau0-tau)
        divfac=x*x
    else
        x=(CP%tau0-tau)/CP%r
        divfac=(CP%r*rofChi(x))**2*k2
    end if


    if (EV%TightCoupling) then
        if (second_order_tightcoupling) then
            pigdot = EV%pigdot
            ypolprime(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopac(j)/opac(j)**2))
        else
            pigdot = -dopac(j)/opac(j)*pig + 32._dl/45*k/opac(j)*(-2*adotoa*sigma  &
                +etak/EV%Kf(1)-  dgpi/k +vbdot )
            ypolprime(2)= pigdot/4
        end if
    end if

    pidot_sum =  pidot_sum + grhog_t*pigdot + grhor_t*pirdot
    diff_rhopi = pidot_sum - (4 * dgpi + dgpi_diff) * adotoa + &
        CP%DarkEnergy%diff_rhopi_Add_Term(grho, gpres, w_dark_energy_t, grhok, adotoa, &
            EV%kf(1), k, grhov_t, z, k2, yprime, y, EV%w_ix)

    !Maple's fortran output - see scal_eqs.map
    !2phi' term (\phi' + \psi' in Newtonian gauge)
    ISW = (4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
        -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*gpres+5.D0*grho)*sigma/k/3.D0 &
        -2.D0/k*adotoa/EV%Kf(1)*etak)*expmmu(j)

    !e.g. to get only late-time ISW
    !  if (1/a-1 < 30) ISW=0

    !The rest, note y(9)->octg, yprime(9)->octgprime (octopoles)
    sources(1)= ISW +  ((-9.D0/160.D0*pig-27.D0/80.D0*ypol(2))/k**2*opac(j)+ &
        (11.D0/10.D0*sigma- 3.D0/8.D0*EV%Kf(2)*ypol(3)+vb-9.D0/80.D0*EV%Kf(2)*octg+3.D0/40.D0*qg)/k- &
        (-180.D0*ypolprime(2)-30.D0*pigdot)/k**2/160.D0)*dvis(j) + &
        (-(9.D0*pigdot+ 54.D0*ypolprime(2))/k**2*opac(j)/160.D0+pig/16.D0+clxg/4.D0+3.D0/8.D0*ypol(2) + &
        (-21.D0/5.D0*adotoa*sigma-3.D0/8.D0*EV%Kf(2)*ypolprime(3) + &
        vbdot+3.D0/40.D0*qgdot- 9.D0/80.D0*EV%Kf(2)*octgprime)/k + &
        (-9.D0/160.D0*dopac(j)*pig-21.D0/10.D0*dgpi-27.D0/80.D0*dopac(j)*ypol(2))/k**2)*vis(j) + &
        (3.D0/16.D0*ddvis(j)*pig+9.D0/8.D0*ddvis(j)*ypol(2))/k**2+21.D0/10.D0/k/EV%Kf(1)*vis(j)*etak

    ! Doppler term
    !   sources(1)=  (sigma+vb)/k*dvis(j)+((-2.D0*adotoa*sigma+vbdot)/k-1.D0/k**2*dgpi)*vis(j) &
    !         +1.D0/k/EV%Kf(1)*vis(j)*etak

    !Equivalent full result
    !    t4 = 1.D0/adotoa
    !    t92 = k**2
    !    sources(1) = (4.D0/3.D0*EV%Kf(1)*expmmu(j)*sigma+2.D0/3.D0*(-sigma-t4*etak)*expmmu(j))*k+ &
    !        (3.D0/8.D0*ypol(2)+pig/16.D0+clxg/4.D0)*vis(j)
    !    sources(1) = sources(1)-t4*expmmu(j)*dgrho/3.D0+((11.D0/10.D0*sigma- &
    !         3.D0/8.D0*EV%Kf(2)*ypol(3)+vb+ 3.D0/40.D0*qg-9.D0/80.D0*EV%Kf(2)*y(9))*dvis(j)+(5.D0/3.D0*grho+ &
    !        gpres)*sigma*expmmu(j)+(-2.D0*adotoa*etak*expmmu(j)+21.D0/10.D0*etak*vis(j))/ &
    !        EV%Kf(1)+(vbdot-3.D0/8.D0*EV%Kf(2)*ypolprime(3)+3.D0/40.D0*qgdot-21.D0/ &
    !        5.D0*sigma*adotoa-9.D0/80.D0*EV%Kf(2)*yprime(9))*vis(j))/k+(((-9.D0/160.D0*pigdot- &
    !        27.D0/80.D0*ypolprime(2))*opac(j)-21.D0/10.D0*dgpi -27.D0/80.D0*dopac(j)*ypol(2) &
    !        -9.D0/160.D0*dopac(j)*pig)*vis(j) - diff_rhopi*expmmu(j)+((-27.D0/80.D0*ypol(2)-9.D0/ &
    !        160.D0*pig)*opac(j)+3.D0/16.D0*pigdot+9.D0/8.D0*ypolprime(2))*dvis(j)+9.D0/ &
    !        8.D0*ddvis(j)*ypol(2)+3.D0/16.D0*ddvis(j)*pig)/t92


    if (x > 0._dl) then
        !E polarization source
        sources(2)=vis(j)*polter*(15._dl/8._dl)/divfac
        !factor of four because no 1/16 later
    else
        sources(2)=0
    end if

    if (CTransScal%NumSources > 2) then
        !Get lensing sources
        phi = -(dgrho +3*dgq*adotoa/k)/(k2*EV%Kf(1)*2) - dgpi/k2/2
        !Can modify this here if you want to get power spectra for other tracer
        !CMB lensing sources
        if (tau > tau_maxvis .and. CP%tau0-tau > 0.1_dl) then
            !phi_lens = Phi - 1/2 kappa (a/k)^2 sum_i rho_i pi_i
            sources(3) = -2* phi *f_K(tau-tau_maxvis)/ &
                (f_K(CP%tau0-tau_maxvis)*f_K(CP%tau0-tau))
            !We include the lensing factor of two here

            !!!!!
            !!      if (lens21cm) then
            !        w_ix=1
            !      if (tau>Redshift_W(w_ix)%tau) then
            !          sources(3) = -2*phi*f_K(tau-Redshift_W(w_ix)%tau)/(f_K(CP%tau0-Redshift_W(w_ix)%tau)*f_K(CP%tau0-tau))
            !         else
            !          sources(3) = 0
            !        end if
            !       else
            !         sources(3) = 0
            !     end if
        else
            sources(3) = 0
        end if

        if (line_reionization) sources(2)=0

        if (tau>tau_start_redshiftwindows .or. plot_evolve) then
            !There are line of sight contributions...

            if (Do21cm) then
                Delta_TCMB = clxg/4
                Delta_source = clxb
                Trad = CP%TCMB/a

                xe = Recombination_xe(a)
                Tmat = Recombination_Tm(a)

                if (EV%Evolve_TM) then
                    Delta_tm = y(EV%Tg_ix)
                else
                    Delta_tm = Delta_TCMB + (1-Tmat/Trad)*4*Delta_TCMB
                end if
                if (Evolve_delta_xe .and. .not. EV%Saha) then
                    Delta_xe = y(EV%xe_ix)
                else
                    Delta_xe = 0
                end if
                Delta_source2 = Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, &
                    k*(z+vb)/adotoa/3)
            end if


            if (plot_evolve) then
                Tspin = Recombination_Ts(a)

                call thermo(tau,cs2,opacity)

                delta_p = barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe)*Tmat*(clxb + delta_tm)
                xe=Recombination_xe(a)
                stop 'write code in equations.f90::output'
                !       write(*,'(5e15.5)') 1/a-1, clxb + Trad/(Tspin-Trad)*delta_source2, delta_source2,  Trad/(Tspin-Trad), clxb

                write(*,'(9e15.5)') 1/a-1, y(EV%xe_ix), clxb, Recombination_xe(a), clxg, tau, delta_TM, delta_p, clxc
                return
            end if


            do w_ix = 1, num_redshiftwindows
                associate (W => Redshift_W(w_ix))

                    if (W%kind == window_lensing) then
                        sources(3+w_ix) =-2*phi*W%win_lens(j)
                    elseif (W%kind == window_counts) then
                        !assume zero velocity bias and relevant tracer is CDM perturbation
                        !neglect anisotropic stress in some places

                        !phidot neglecting anisotropic stress so phi=psi
                        phidot = ((4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
                            -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*gpres+5.D0*grho)*sigma/k/3.D0 &
                            -2.D0/k*adotoa/EV%Kf(1)*etak)) / 2

                        !Main density source
                        if (counts_density) then
                            counts_density_source= W%wing(j)*(clxc*W%bias + (W%comoving_density_ev(j) - 3*adotoa)*sigma/k)
                            !Newtonian gauge count density; bias assumed to be on synchronous gauge CDM density
                        else
                            counts_density_source= 0
                        endif


                        if (counts_redshift) then
                            !Main redshift distortion from kV_N/H j'' integrated by parts twice (V_N = sigma in synch gauge)
                            counts_redshift_source = ((4.D0*adotoa**2+gpres+grho/3.D0)/k*W%wing2(j)+ &
                                (-4.D0*W%dwing2(j)*adotoa+W%ddwing2(j))/k)*sigma+(-etak/adotoa*k/3.D0-dgrho/ &
                                adotoa/6.D0+(etak/adotoa*k/3.D0+dgrho/adotoa/6.D0+(dgq/2.D0-2.D0*etak*adotoa)/k) &
                                /EV%Kf(1))*W%wing2(j)+2.D0*W%dwing2(j)*etak/k/EV%Kf(1)
                            if (k>0.8e-2) then
                                !    write(*,'(8E15.5)') 1/a-1, k*sigma/adotoa, W%wing(j)*clxc, source_redshift, ISW, W%wing(j),W%wing2(j),W%wingtau(j)
                                !    if (1/a-1 < 0.01) stop
                            end if
                        else
                            counts_redshift_source= 0
                        end if

                        ! 2v j'/(H\chi) geometric term
                        if (CP%tau0-tau > 0.1_dl .and. counts_radial) then
                            chi =  CP%tau0-tau
                            counts_radial_source= (1-2.5*W%dlog10Ndm)*((-4.D0*W%wing2(j)/chi*adotoa &
                                -2.D0*(-W%dwing2(j)*chi-W%wing2(j))/chi**2)/ &
                                k*sigma+2.D0*W%wing2(j)*etak/chi/k/EV%Kf(1))
                        else
                            counts_radial_source = 0
                        end if

                        !        if (counts_evolve) then
                        !       !Just source evolution term if window is actual source distribution
                        !         counts_evolve_source =(2.D0*W%dwing2(j)*adotoa-W%ddwing2(j))/k*sigma-W%dwing2(j)*etak/k/EV%Kf(1)
                        !        else
                        !          counts_evolve_source = 0
                        !        end if

                        if (counts_timedelay) then
                            !time delay; WinV is int g/chi
                            counts_timedelay_source= 2*(1-2.5*W%dlog10Ndm)*W%WinV(j)*2*phi
                        else
                            counts_timedelay_source = 0
                        end if

                        if (counts_ISW) then
                            !WinF is int wingtau
                            counts_ISW_source = W%WinF(j)*2*phidot
                        else
                            counts_ISW_source = 0
                        end if

                        if (counts_potential) then
                            !approx phi = psi
                            counts_potential_source = ( phidot/adotoa + phi +(5*W%dlog10Ndm-2)*phi ) * W%wing(j) + phi * W%wingtau(j)
                        else
                            counts_potential_source = 0
                        end if

                        if (counts_velocity) then
                            counts_velocity_source =  (-2.D0*W%wingtau(j)*adotoa+W%dwingtau(j))/k*sigma+W%wingtau(j)*etak/k/EV%Kf(1) &
                                - counts_radial_source  !don't double count terms; counts_radial is part of counts_velocity with 1/H/chi
                        else
                            counts_velocity_source = 0
                        end if

                        sources(3+w_ix)=  counts_radial_source +  counts_density_source + counts_redshift_source &
                            + counts_timedelay_source + counts_potential_source &
                            + counts_ISW_source + counts_velocity_source

                        sources(3+w_ix)=sources(3+w_ix)/W%Fq

                        if (DoRedshiftLensing) &
                            sources(3+W%mag_index+num_redshiftwindows) = phi*W%win_lens(j)*(2-5*W%dlog10Ndm)
                    elseif (W%kind == window_21cm) then
                        if (line_basic) then
                            sources(3+w_ix)= expmmu(j)*(W%wing(j)*Delta_source + W%wing2(j)*Delta_source2 &
                                - W%Wingtau(j)*(clxb - (Delta_source2+clxg/4)))
                            !!    sources(3+w_ix)= expmmu(j)*W%wing(j)*phi
                        else
                            sources(3+w_ix)= 0
                        end if

                        if (line_distortions ) then
                            !With baryon velocity, dropping small terms
                            s(1) =  (sigma/adotoa/3.D0-etak/adotoa**2/3.D0)*W%wing(j)*expmmu(j)*k
                            s(2) =  -1.D0/adotoa**2*expmmu(j)*W%wing(j)*dgrho/6.D0+((((4.D0*sigma+ &
                                vb)*adotoa+(-grho*sigma/2.D0-vb*grho/3.D0)/adotoa+(sigma*grho**2/18.D0+ &
                                vb*grho**2/18.D0)/adotoa**3)*W%wing(j)-4.D0*W%dwing(j)*sigma+(W%ddwing(j)*sigma+ &
                                W%ddwing(j)*vb)/adotoa+(W%dwing(j)*sigma*grho/3.D0+W%dwing(j)*vb*grho/3.D0)/ &
                                adotoa**2-2.D0*W%dwing(j)*vb+((-2.D0*etak+etak*grho/adotoa**2/3.D0)*W%wing(j) &
                                + 2.D0*W%dwing(j)*etak/adotoa)/EV%Kf(1))*expmmu(j)+&
                                (-4.D0*vis(j)*sigma- 2.D0*vis(j)*vb+(dvis(j)*sigma+dvis(j)*vb)/adotoa+(vis(j)*grho*sigma/3.D0+ &
                                vis(j)*vb*grho/3.D0)/adotoa**2)*W%wing(j)+2.D0*vis(j)*etak/adotoa*W%wing(j)/ &
                                EV%Kf(1)+(2.D0*vis(j)*W%dwing(j)*sigma+2.D0*vis(j)*W%dwing(j)*vb)/adotoa)/k
                            t(0) =  s(1)+s(2)

                            sources(3+w_ix)= sources(3+w_ix) + t(0)
                        end if


                        if (line_extra) then
                            !All sources except below
                            if (line_basic .and. line_distortions) then
                                sources(3+w_ix) =  (-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)*expmmu(j)*k+ &
                                    (W%wing2(j)*Delta_source2+W%wing(j)*Delta_source+1.D0/adotoa*W%winV(j)*dgrho/3.D0)* &
                                    expmmu(j)+((-W%dwing(j)*vb+(-(3.D0*gpres+grho)*sigma/3.D0 &
                                    - 4.D0*adotoa**2*sigma)*W%winV(j)+4.D0*adotoa*W%dwinV(j)*sigma+(-sigma- &
                                    vb)*W%ddWinV(j)-vbdot*W%wing(j)-W%dwinV(j)*vbdot+(-2.D0*W%dwinV(j)*etak &
                                    + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*expmmu(j)-2.D0*vis(j)*sigma*W%dwinV(j)+ &
                                    (4.D0*vis(j)*sigma*adotoa-dvis(j)*sigma)*W%winV(j)-2.D0*vis(j)*W%winV(j)*etak/ &
                                    EV%Kf(1)-vis(j)*W%dwinV(j)*vb-vis(j)*W%wing(j)*vb)/k+((2.D0*W%dwinV(j)*dgpi+ &
                                    diff_rhopi*W%winV(j))*expmmu(j)+2.D0*vis(j)*W%winV(j)*dgpi)/k**2
                            else
                                s(1) =  ((-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)+(-sigma/adotoa/3.D0+ &
                                    etak/adotoa**2/3.D0)*W%wing(j))*expmmu(j)*k+(1.D0/adotoa*W%winV(j)*dgrho/3.D0 &
                                    + 1.D0/adotoa**2*W%wing(j)*dgrho/6.D0)*expmmu(j)
                                s(2) =  s(1)
                                s(6) =  ((-vb-sigma)*W%ddWinV(j)+(-4.D0*adotoa**2*sigma-&
                                    (18.D0*gpres+ 6.D0*grho)*sigma/18.D0)*W%winV(j)+((-4.D0*sigma-vb)*adotoa-vbdot+&
                                    (grho*sigma/ 2.D0+vb*grho/3.D0)/adotoa+(-grho**2*sigma/18.D0-vb*grho**2/18.D0)/ &
                                    adotoa**3)*W%wing(j)+W%dwing(j)*vb+(-W%ddwing(j)*sigma-W%ddwing(j)*vb)/adotoa &
                                    + 4.D0*W%dwinV(j)*sigma*adotoa+4.D0*W%dwing(j)*sigma+(-W%dwing(j)*grho*sigma/3.D0- &
                                    W%dwing(j)*vb*grho/3.D0)/adotoa**2-W%dwinV(j)*vbdot+((2.D0*etak-etak*grho/ &
                                    adotoa**2/3.D0)*W%wing(j)-2.D0*W%dwing(j)*etak/adotoa-2.D0*W%dwinV(j)*etak &
                                    + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*expmmu(j)-vis(j)*W%dwinV(j)*vb+ &
                                    (4.D0*vis(j)*sigma*adotoa-dvis(j)*sigma)*W%winV(j)
                                s(5) =  s(6)+(-2.D0*vis(j)*etak/adotoa*W%wing(j)-2.D0*vis(j)*W%winV(j)*etak)/ &
                                    EV%Kf(1)+(4.D0*vis(j)*sigma+(-vis(j)*grho*sigma/3.D0-vis(j)*vb*grho/3.D0)/ &
                                    adotoa**2+vis(j)*vb+(-dvis(j)*sigma-dvis(j)*vb)/adotoa)*W%wing(j)+ &
                                    (-2.D0*vis(j)*W%dwing(j)*sigma-2.D0*vis(j)*W%dwing(j)*vb)/adotoa &
                                    - 2.D0*vis(j)*W%dwinV(j)*sigma
                                s(6) =  1.D0/k
                                s(4) =  s(5)*s(6)
                                s(5) =  ((diff_rhopi*W%winV(j)+2.D0*W%dwinV(j)*dgpi)*expmmu(j) &
                                    + 2.D0*vis(j)*dgpi*W%winV(j))/k**2
                                s(3) =  s(4)+s(5)
                                t(0) =  s(2)+s(3)

                                sources(3+w_ix) =   sources(3+w_ix) + t(0)
                            end if
                        end if



                        if (line_reionization) then
                            if (num_redshiftwindows>1) stop 'reionization only for one window at the mo'
                            lineoff=EV%reion_line_ix
                            lineoffpol = lineoff+EV%lmaxline-1
                            polter_line = 0.1_dl*y(lineoff+2)+9._dl/15._dl*y(lineoffpol+2)

                            if (x > 0._dl) then
                                sources(2)=vis(j)*polter_line*(15._dl/2._dl)/divfac
                            else
                                sources(2)=0
                            end if

                            if (.not. use_mK) sources(2)= sources(2) /W%Fq

                            s(1) =  vis(j)*y(lineoff+2)/4.D0+vis(j)*y(lineoff)
                            s(2) =  s(1)
                            s(4) =  (-1.D0/EV%Kf(1)*vis(j)*W%winV(j)*etak/10.D0-vis(j)*sigma*W%dwinV(j)/10.D0 &
                                - 9.D0/20.D0*vis(j)*yprime(lineoff+2)-27.D0/100.D0*vis(j)*opac(j)*y(lineoff+1) &
                                - 9.D0/10.D0*dvis(j)*y(lineoff+3)-3.D0/20.D0*vis(j)*opac(j)*EV%Kf(2)*y(lineoffpol+3)+ &
                                vis(j)*W%dwinV(j)*vb+81.D0/200.D0*vis(j)*opac(j)*y(lineoff+3) &
                                +3.D0/5.D0*dvis(j)*y(lineoff+1)+3.D0/10.D0*vis(j)*yprime(lineoff+1)+ &
                                (vis(j)*adotoa*sigma/5.D0+(36.D0*vis(j)*opac(j)-80.D0*dvis(j))*sigma/400.D0+ &
                                dvis(j)*vb+vis(j)*vbdot)*W%winV(j))/k
                            s(5) =  (vis(j)*W%winV(j)*dgpi/10.D0+9.D0/20.D0*vis(j)*dopac(j)*y(lineoffpol+2) &
                                + 261.D0/400.D0*vis(j)*opac(j)**2.D0*y(lineoff+2)&
                                -117.D0/200.D0*vis(j)*opac(j)**2.D0*y(lineoffpol+2)+3.D0/4.D0*ddvis(j)*y(lineoff+2) &
                                - 27.D0/20.D0*dvis(j)*opac(j)*y(lineoff+2)+9.D0/10.D0*dvis(j)*opac(j)*y(lineoffpol+2)&
                                -27.D0/40.D0*vis(j)*dopac(j)*y(lineoff+2))/k**2
                            s(3) =  s(4)+s(5)
                            t(0) =  s(2)+s(3)

                            sources(3+w_ix)= sources(3+w_ix) + t(0)
                        end if

                        if (line_phot_quadrupole) then
                            s(1) =  (EV%kf(1)*W%wing2(j)*pig/2.D0+(-clxg/4.D0-5.D0/8.D0*pig)*W%wing2(j))*expmmu(j)
                            s(3) =  ((-1.D0/EV%kf(1)*W%wing2(j)*etak+(-sigma+9.D0/8.D0*EV%kf(2)*y(9) &
                                -3.D0/4.D0*qg)*W%dwing2(j)+(-opac(j)*vb+2.D0*adotoa*sigma+9.D0/8.D0*EV%kf(2)*yprime(9) &
                                + 3.D0/8.D0*opac(j)*EV%kf(2)*ypol(3)+3.D0/4.D0*opac(j)*qg)*W%wing2(j))*expmmu(j)+ &
                                (-3.D0/4.D0*vis(j)*qg-vis(j)*sigma+9.D0/8.D0*vis(j)*EV%kf(2)*y(9))*W%wing2(j))/k
                            s(4) =  (((27.D0/16.D0*opac(j)*pig-15.D0/8.D0*pigdot &
                                -9.D0/8.D0*opac(j)*ypol(2))*W%dwing2(j)+(27.D0/16.D0*dopac(j)*pig &
                                +9.D0/8.D0*opac(j)**2.D0*ypol(2)-9.D0/8.D0*opac(j)**2.D0*polter &
                                +27.D0/16.D0*opac(j)*pigdot+dgpi-9.D0/8.D0*dopac(j)*ypol(2))*W%wing2(j)&
                                -15.D0/8.D0*W%ddwing2(j)*pig)*expmmu(j)-15.D0/4.D0*vis(j)*W%dwing2(j)*pig+(- &
                                (-27.D0*vis(j)*opac(j)+30.D0*dvis(j))*pig/16.D0-9.D0/8.D0*vis(j)*opac(j)*ypol(2) &
                                - 15.D0/8.D0*vis(j)*pigdot)*W%wing2(j))/k**2
                            s(2) =  s(3)+s(4)
                            t(0) =  s(1)+s(2)

                            sources(3+w_ix)= sources(3+w_ix)+ t(0)
                        end if


                        if (line_phot_dipole) then
                            sources(3+w_ix)=sources(3+w_ix) + (EV%kf(1)*W%wing2(j)*pig/2.D0-W%wing2(j)*clxg/4.D0)*expmmu(j)+(((vbdot- &
                                opac(j)*vb+3.D0/4.D0*opac(j)*qg)*W%wing2(j)+(vb-3.D0/4.D0*qg)*W%dwing2(j))*expmmu(j)+&
                                (vis(j)*vb-3.D0/4.D0*vis(j)*qg)*W%wing2(j))/k
                        end if

                        if (.not. use_mK) sources(3+w_ix)= sources(3+w_ix) /W%Fq
                    end if
                end associate
            end do
        end if
    end if !num sources > 2

    end subroutine output


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine outputt(EV,yt,n,j,tau,dt,dte,dtb)
    !calculate the tensor sources for open and closed case
    use ThermoData

    implicit none
    integer j,n
    type(EvolutionVars) :: EV
    real(dl), target :: yt(n), ytprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot,polterddot,prefac
    real(dl) pig, pigdot, octg, aux, polter, shear, adotoa,a
    real(dl) sinhxr,cothxor
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Bprime,Eprime
    real(dl), target :: pol(3),polEprime(3), polBprime(3)

    call derivst(EV,EV%nvart,tau,yt,ytprime)

    k2=EV%k2_buf
    k=EV%k_buf
    aux=EV%aux_buf
    shear = yt(3)

    x=(CP%tau0-tau)/CP%r

    !  And the electric part of the Weyl.
    if (.not. EV%TensTightCoupling) then
        !  Use the full expression for pigdt
        pig=yt(EV%g_ix+2)
        pigdot=ytprime(EV%g_ix+2)
        E => yt(EV%E_ix+1:)
        Eprime=> ytprime(EV%E_ix+1:)
        Bprime => ytprime(EV%B_ix+1:)
        octg=ytprime(EV%g_ix+3)
    else
        !  Use the tight-coupling approximation
        a =yt(1)
        adotoa = 1/(a*dtauda(a))
        pigdot=32._dl/45._dl*k/opac(j)*(2._dl*adotoa*shear+ytprime(3))
        pig = 32._dl/45._dl*k/opac(j)*shear
        pol=0
        polEprime=0
        polBprime=0
        E=>pol
        EPrime=>polEPrime
        BPrime=>polBPrime
        E(2)=pig/4._dl
        EPrime(2)=pigdot/4
        octg=0
    endif

    sinhxr=rofChi(x)*CP%r

    if (EV%q*sinhxr > 1.e-8_dl) then
        prefac=sqrt(EV%q2*CP%r*CP%r-CP%Ksign)
        cothxor=cosfunc(x)/sinhxr

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
        polterddot = 9._dl/15._dl*(-dopac(j)*(E(2)-polter)-opac(j)*(  &
            Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - 5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
            +0.1_dl*(k*(-octg*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(3)) - &
            dopac(j)*(pig - polter) - opac(j)*(pigdot-polterdot))

        dt=(shear*expmmu(j) + (15._dl/8._dl)*polter*vis(j)/k)*CP%r/sinhxr**2/prefac

        dte=CP%r*15._dl/8._dl/k/prefac* &
            ((ddvis(j)*polter + 2._dl*dvis(j)*polterdot + vis(j)*polterddot)  &
            + 4._dl*cothxor*(dvis(j)*polter + vis(j)*polterdot) - &
            vis(j)*polter*(k2 -6*cothxor**2))

        dtb=15._dl/4._dl*EV%q*CP%r/k/prefac*(vis(j)*(2._dl*cothxor*polter + polterdot) + dvis(j)*polter)
    else
        dt=0._dl
        dte=0._dl
        dtb=0._dl
    end if

    end subroutine outputt

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine outputv(EV,yv,n,j,tau,dt,dte,dtb)
    !calculate the vector sources
    use ThermoData

    implicit none
    integer j,n
    type(EvolutionVars) :: EV
    real(dl), target :: yv(n), yvprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot
    real(dl) vb,qg, pig, polter, sigma
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Eprime

    call derivsv(EV,EV%nvarv,tau,yv,yvprime)

    k2=EV%k2_buf
    k=EV%k_buf
    sigma = yv(2)
    vb  = yv(3)
    qg  = yv(4)
    pig = yv(5)


    x=(CP%tau0-tau)*k

    if (x > 1.e-8_dl) then
        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)

        if (yv(1) < 1e-3) then
            dt = 1
        else
            dt =0
        end if
        dt= (4*(vb+sigma)*vis(j) + 15._dl/2/k*( vis(j)*polterdot + dvis(j)*polter) &
            + 4*(expmmu(j)*yvprime(2)) )/x

        dte= 15._dl/2*2*polter/x**2*vis(j) + 15._dl/2/k*(dvis(j)*polter + vis(j)*polterdot)/x

        dtb= -15._dl/2*polter/x*vis(j)
    else
        dt=0
        dte=0
        dtb=0
    end if

    end subroutine outputv


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initial(EV,y, tau)
    !  Initial conditions.
    use ThermoData
    implicit none

    type(EvolutionVars) EV
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
        Rc,Rb,Rv,Rg,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, iqg, rhomass,a_massive, ep
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
        i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxq=11,i_vq=12
    integer, parameter :: i_max = i_vq
    real(dl) initv(6,1:i_max), initvec(1:i_max)

    nullify(EV%OutputTransfer) !Should not be needed, but avoids issues in ifort 14

    if (CP%flat) then
        EV%k_buf=EV%q
        EV%k2_buf=EV%q2
        EV%Kf(1:EV%MaxlNeeded)=1._dl
    else
        EV%k2_buf=EV%q2-CP%curv
        EV%k_buf=sqrt(EV%k2_buf)

        do l=1,EV%MaxlNeeded
            EV%Kf(l)=1._dl-CP%curv*(l*(l+2))/EV%k2_buf
        end do
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    do j=1,EV%MaxlNeeded
        EV%denlk(j)=denl(j)*k*j
        EV%denlk2(j)=denl(j)*k*EV%Kf(j)*(j+1)
        EV%polfack(j)=polfac(j)*k*EV%Kf(j)*denl(j)
    end do

    !Get time to switch off tight coupling
    !The numbers here are a bit of guesswork
    !The high k increase saves time for very small loss of accuracy
    !The lower k ones are more delicate. Nead to avoid instabilities at same time
    !as ensuring tight coupling is accurate enough
    if (EV%k_buf > epsw) then
        if (EV%k_buf > epsw*5) then
            ep=ep0*5/AccuracyBoost
            if (HighAccuracyDefault) ep = ep*0.65
        else
            ep=ep0
        end if
    else
        ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))

    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+grhornomass

    om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))
    omtau=om*tau
    Rv=grhonu/(grhonu+grhog)

    Rg = 1-Rv
    Rc=CP%omegac/(CP%omegac+CP%omegab)
    Rb=1-Rc
    Rp15=4*Rv+15

    if (CP%Scalar_initial_condition > initial_nummodes) &
        stop 'Invalid initial condition for scalar modes'

    a=tau*adotrad*(1+omtau/4)
    a2=a*a

    initv=0

    !  Set adiabatic initial conditions

    chi=1  !Get transfer function for chi
    initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
    initv(1,i_clxr)= initv(1,i_clxg)
    initv(1,i_clxb)=0.75_dl*initv(1,i_clxg)
    initv(1,i_clxc)=initv(1,i_clxb)
    initv(1,i_qg)=initv(1,i_clxg)*x/9._dl
    initv(1,i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
    initv(1,i_vb)=0.75_dl*initv(1,i_qg)
    initv(1,i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
    initv(1,i_aj3r)=chi*4/21._dl/Rp15*x3
    initv(1,i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))

    if (CP%Scalar_initial_condition/= initial_adiabatic) then
        !CDM isocurvature

        initv(2,i_clxg)= Rc*omtau*(-2._dl/3 + omtau/4)
        initv(2,i_clxr)=initv(2,i_clxg)
        initv(2,i_clxb)=initv(2,i_clxg)*0.75_dl
        initv(2,i_clxc)=1+initv(2,i_clxb)
        initv(2,i_qg)=-Rc/9*omtau*x
        initv(2,i_qr)=initv(2,i_qg)
        initv(2,i_vb)=0.75_dl*initv(2,i_qg)
        initv(2,i_pir)=-Rc*omtau*x2/3/(2*Rv+15._dl)
        initv(2,i_eta)= Rc*omtau*(1._dl/3 - omtau/8)*EV%Kf(1)
        initv(2,i_aj3r)=0
        !Baryon isocurvature
        if (Rc==0) stop 'Isocurvature initial conditions assume non-zero dark matter'

        initv(3,:) = initv(2,:)*(Rb/Rc)
        initv(3,i_clxc) = initv(3,i_clxb)
        initv(3,i_clxb)= initv(3,i_clxb)+1

        !neutrino isocurvature density mode

        initv(4,i_clxg)=Rv/Rg*(-1 + x2/6)
        initv(4,i_clxr)=1-x2/6
        initv(4,i_clxc)=-omtau*x2/80*Rv*Rb/Rg
        initv(4,i_clxb)= Rv/Rg/8*x2
        iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
        initv(4,i_qg) =iqg
        initv(4,i_qr) = x/3
        initv(4,i_vb)=0.75_dl*iqg
        initv(4,i_pir)=x2/Rp15
        initv(4,i_eta)=EV%Kf(1)*Rv/Rp15/3*x2

        !neutrino isocurvature velocity mode

        initv(5,i_clxg)=Rv/Rg*x - 2*x*omtau/16*Rb*(2+Rg)/Rg**2
        initv(5,i_clxr)=-x -3*x*omtau*Rb/16/Rg
        initv(5,i_clxc)=-9*omtau*x/64*Rv*Rb/Rg
        initv(5,i_clxb)= 3*Rv/4/Rg*x - 9*omtau*x/64*Rb*(2+Rg)/Rg**2
        iqg = Rv/Rg*(-1 + 3*Rb/4/Rg*omtau+x2/6 +3*omtau**2/16*Rb/Rg**2*(Rg-3*Rb))
        initv(5,i_qg) =iqg
        initv(5,i_qr) = 1 - x2/6*(1+4*EV%Kf(1)/(4*Rv+5))
        initv(5,i_vb)=0.75_dl*iqg
        initv(5,i_pir)=2*x/(4*Rv+5)+omtau*x*6/Rp15/(4*Rv+5)
        initv(5,i_eta)=2*EV%Kf(1)*x*Rv/(4*Rv+5) + omtau*x*3*EV%Kf(1)*Rv/32*(Rb/Rg - 80/Rp15/(4*Rv+5))
        initv(5,i_aj3r) = 3._dl/7*x2/(4*Rv+5)

        !quintessence isocurvature mode
    end if

    if (CP%Scalar_initial_condition==initial_vector) then
        InitVec = 0
        do i=1,initial_nummodes
            InitVec = InitVec+ initv(i,:)*CP%InitialConditionVector(i)
        end do
    else
        InitVec = initv(CP%Scalar_initial_condition,:)
        if (CP%Scalar_initial_condition==initial_adiabatic) InitVec = -InitVec
        !So we start with chi=-1 as before
    end if

    y(1)=a
    y(2)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable

    !  CDM
    y(3)=InitVec(i_clxc)

    !  Baryons
    y(4)=InitVec(i_clxb)
    y(5)=InitVec(i_vb)

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    ! DarkEnergy: This initializes also i_vq, when num_perturb_equations is set
    !             to 2.
    if (CP%DarkEnergy%num_perturb_equations > 0) &
        y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
            InitVec(i_clxq:i_clxq + CP%DarkEnergy%num_perturb_equations - 1)

    !  Neutrinos
    y(EV%r_ix)=InitVec(i_clxr)
    y(EV%r_ix+1)=InitVec(i_qr)
    y(EV%r_ix+2)=InitVec(i_pir)

    if (EV%lmaxnr>2) then
        y(EV%r_ix+3)=InitVec(i_aj3r)
    endif

    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
        EV%MassiveNuApproxTime(nu_i) = Nu_tau_massive(nu_i)
        a_massive =  20000*k/nu_masses(nu_i)*AccuracyBoost*lAccuracyBoost
        if (a_massive >=0.99) then
            EV%MassiveNuApproxTime(nu_i)=CP%tau0+1
        else if (a_massive > 17.d0/nu_masses(nu_i)*AccuracyBoost) then
            EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),DeltaTime(0._dl,a_massive, 0.01_dl))
        end if
        ind = EV%nu_ix(nu_i)
        do  i=1,EV%nq(nu_i)
            y(ind:ind+2)=y(EV%r_ix:EV%r_ix+2)
            if (EV%lmaxnu_tau(nu_i)>2) y(ind+3)=InitVec(i_aj3r)
            ind = ind + EV%lmaxnu_tau(nu_i)+1
        end do
    end do

    end subroutine initial


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialt(EV,yt,tau)
    !  Initial conditions for tensors
    use ThermoData
    implicit none
    real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
    integer l
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yt(EV%nvart)
    real(dl) tens0, ep, tensfac

    if (CP%flat) then
        EV%aux_buf=1._dl
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
        EV%Kft(1:EV%MaxlNeededt)=1._dl !initialize for flat case
    else
        EV%k2_buf=EV%q2-3*CP%curv
        EV%k_buf=sqrt(EV%k2_buf)
        EV%aux_buf=sqrt(1._dl+3*CP%curv/EV%k2_buf)
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    do l=1,EV%MaxlNeededt
        if (.not. CP%flat) EV%Kft(l)=1._dl-CP%curv*((l+1)**2-3)/k2
        EV%denlkt(1,l)=k*denl(l)*l !term for L-1
        tensfac=real((l+3)*(l-1),dl)/(l+1)
        EV%denlkt(2,l)=k*denl(l)*tensfac*EV%Kft(l) !term for L+1
        EV%denlkt(3,l)=k*denl(l)*tensfac**2/(l+1)*EV%Kft(l) !term for polarization
        EV%denlkt(4,l)=k*4._dl/(l*(l+1))*EV%aux_buf !other for polarization
    end do

    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if

    !    finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep))
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))

    a=tau*adotrad
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+rhomass+grhornomass))

    if (DoTensorNeutrinos) then
        bigR = (rhomass+grhornomass)/(rhomass+grhornomass+grhog)
    else
        bigR = 0._dl
    end if

    x=k*tau

    yt(1)=a
    tens0 = 1

    yt(2)= tens0
    !commented things are for the compensated mode with magnetic fields; can be neglected
    !-15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))

    elec=-tens0*(1+2*CP%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1

    !shear
    yt(3)=-5._dl/2/(bigR+5)*x*elec
    !          + 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))

    yt(4:EV%nvart)=0._dl

    !  Neutrinos
    if (DoTensorNeutrinos) then
        pir=-2._dl/3._dl/(bigR+5)*x**2*elec
        !           + (bigR-1)/bigR*Magnetic*(1-15./14*x**2/(15+4*bigR))
        aj3r=  -2._dl/21._dl/(bigR+5)*x**3*elec !&
            !           + 3._dl/7*x*(bigR-1)/bigR*Magnetic
        yt(EV%r_ix+2)=pir
        yt(EV%r_ix+3)=aj3r
        !Should set up massive too, but small anyway..
    end if

    end subroutine initialt

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialv(EV,yv,tau)
    !  Initial conditions for vectors

    implicit none
    real(dl) bigR,Rc,tau,x,pir
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yv(EV%nvarv)

    if (CP%flat) then
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
    else
        stop 'Vectors not supported in non-flat models'
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+grhornomass))

    a=tau*adotrad*(1+omtau/4)

    x=k*tau

    bigR = (grhornomass)/(grhornomass+grhog)
    Rc=CP%omegac/(CP%omegac+CP%omegab)

    yv(1)=a


    yv(2)= vec_sig0*(1- 15._dl/2*omtau/(4*bigR+15)) + 45._dl/14*x*Magnetic*(BigR-1)/(4*BigR+15)
    !qg
    yv(4)= vec_sig0/3* (4*bigR + 5)/(1-BigR)*(1  -0.75_dl*omtau*(Rc-1)/(bigR-1)* &
        (1 - 0.25_dl*omtau*(3*Rc-2-bigR)/(BigR-1))) &
        -x/2*Magnetic
    yv(3)= 3._dl/4*yv(4)

    yv(5:EV%nvarv) = 0

    !        if (.false.) then
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = vec_sig0/6/bigR*x**2*(1+2*bigR*omtau/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2) = -2/3._dl*vec_sig0/bigR*x*(1 +3*omtau*bigR/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+3) = 1/4._dl*vec_sig0/bigR*(5+4*BigR)
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+4) =1/9.*x*vec_sig0*(5+4*bigR)/bigR
    !         yv(4) = 0
    !         yv(3)= 3._dl/4*yv(4)
    !          return
    !        end if

    !  Neutrinos
    !q_r
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = -1._dl/3*vec_sig0*(4*BigR+5)/bigR &
        + x**2*vec_sig0/6/BigR +0.5_dl*x*(1/bigR-1)*Magnetic
    !pi_r
    pir=-2._dl/3._dl*x*vec_sig0/BigR - (1/bigR-1)*Magnetic
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +1)=pir
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +2)=3._dl/7*x*Magnetic*(1-1/BigR)

    end subroutine initialv


    subroutine outtransf(EV, y,tau, Arr)
    !write out clxc, clxb, clxg, clxn
    implicit none
    type(EvolutionVars) EV
    real(dl), intent(in) :: tau
    real, target :: Arr(:)
    real(dl) y(EV%nvar),yprime(EV%nvar)

    yprime = 0
    EV%OutputTransfer =>  Arr
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputTransfer)

    Arr(Transfer_kh+1:Transfer_max) = Arr(Transfer_kh+1:Transfer_max)/EV%k2_buf

    end subroutine outtransf


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine derivs(EV,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the perturbations
    !  ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
    use constants, only : barssc0, Compton_CT, line21_const
    use ThermoData
    use MassiveNu
    use Recombination
    use RECDATA, only : CB1

    implicit none
    type(EvolutionVars) EV

    integer n,nu_i
    real(dl) ay(n),ayprime(n)
    real(dl) tau,w
    real(dl) k,k2

    !  Internal variables.

    real(dl) opacity
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
    real(dl) w_dark_energy_t !equation of state of dark energy
    real(dl) gpres_matter !pressure from massive neutrinos
    real(dl) gpres_noDE !Pressure with matter and radiation, no dark energy
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) E2, dopacity
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz !, ddz
    real(dl) dgpi,dgrho_matter,grho_matter, clxnu_all
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    real(dl) xe,Trad, Delta_TM, Tmat, Delta_TCMB
    real(dl) delta_p, wing_t, wing2_t,winv_t
    real(dl) Delta_source2, polter_line
    real(dl) Delta_xe, Tspin, tau_eps, tau_fac, Tb
    integer lineoff,lineoffpol


    k=EV%k_buf
    k2=EV%k2_buf

    a=ay(1)
    a2=a*a

    etak=ay(2)

    !  CDM variables
    clxc=ay(3)

    !  Baryon variables
    clxb=ay(4)
    vb=ay(5)

    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    gpres_matter = 0
    call CP%DarkEnergy%BackgroundDensityAndPressure(a, grhov_t, w_dark_energy_t)

    !  Get sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call thermo(tau,cs2,opacity,dopacity)
    else
        call thermo(tau,cs2,opacity)
    end if

    grho_matter=grhob_t+grhoc_t

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    if (CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,grho_matter,gpres_matter,dgrho_matter,dgq, wnu_arr)
    end if

    grho = grho_matter+grhor_t+grhog_t+grhov_t
    gpres_noDE = gpres_matter + (grhor_t + grhog_t)/3

    if (CP%flat) then
        adotoa=sqrt(grho/3)
        cothxor=1._dl/tau
    else
        adotoa=sqrt((grho+grhok)/3._dl)
        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
    end if

    dgrho = dgrho_matter
    if (.not. CP%DarkEnergy%is_cosmological_constant) &
        call CP%DarkEnergy%AddStressEnergy(dgrho, dgq, grhov_t, ay, EV%w_ix, .false.)

    if (EV%no_nu_multpoles) then
        !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
        !Approximate total density variables with just matter terms
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxr=-4*dz/k
        qr=-4._dl/3*z
        pir=0
    else
        !  Massless neutrinos
        clxr=ay(EV%r_ix)
        qr  =ay(EV%r_ix+1)
        pir =ay(EV%r_ix+2)
    endif

    pig=0
    if (EV%no_phot_multpoles) then
        if (.not. EV%no_nu_multpoles) then
            z=(0.5_dl*dgrho/k + etak)/adotoa
            dz= -adotoa*z - 0.5_dl*dgrho/k
            clxg=-4*dz/k-4/k*opacity*(vb+z)
            qg=-4._dl/3*z
        else
            clxg=clxr-4/k*opacity*(vb+z)
            qg=qr
        end if
    else
        !  Photons
        clxg=ay(EV%g_ix)
        qg=ay(EV%g_ix+1)
        if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if

    !  8*pi*a*a*SUM[rho_i*clx_i] - radiation terms
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=dgq + grhog_t*qg+grhor_t*qr

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    ayprime(1)=adotoa*a

    call CP%DarkEnergy%DerivsAddPreSigma(sigma, ayprime, dgq, dgrho, &
        grho, grhov_t, w_dark_energy_t, gpres_noDE, ay, EV%w_ix, etak, adotoa, k, k2, EV%kf(1))

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    if (CP%flat) then
        !eta*k equation
        sigma=(z+1.5_dl*dgq/k2)
        ayprime(2)=0.5_dl*dgq
    else
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
        ayprime(2)=0.5_dl*dgq + CP%curv*z
    end if

    call CP%DarkEnergy%PerturbationEvolve(ayprime, EV%w_ix, adotoa, k, z, ay)

    !  CDM equation of motion
    clxcdot=-k*z
    ayprime(3)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(4)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    !Sources
    if (Evolve_baryon_cs) then
        if (a > Do21cm_mina) then
            Tmat = Recombination_Tm(a)
        else
            Tmat = CP%TCMB/a
        end if
        if (EV%Evolve_TM) then
            Delta_TM = ay(EV%Tg_ix)
        else
            Delta_TM = clxg/4
        end if
        delta_p = barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*opacity*a2/akthom)*Tmat*(clxb + delta_tm)
    else
        Delta_TM = clxg/4
        delta_p = cs2*clxb
    end if


    if (Evolve_delta_xe) then
        if (EV%saha) then
            xe=Recombination_xe(a)
            Delta_xe = (1-xe)/(2-xe)*(-clxb + (3._dl/2+  CB1/Tmat)*Delta_TM)
        else
            Delta_xe = ay(EV%xe_ix)
        end if
    else
        Delta_xe = 0
    end if

    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
        !  ddota/a
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        adotdota=(adotoa*adotoa-gpres)/2

        pig = 32._dl/45/opacity*k*(sigma+vb)

        !  First-order approximation to baryon-photon splip
        slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

        if (second_order_tightcoupling) then
            ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming flat
            !AL: First order slip seems to be fine here to 2e-4

            !  8*pi*G*a*a*SUM[rho_i*sigma_i]
            dgs = grhog_t*pig+grhor_t*pir

            ! Define shear derivative to first order
            sigmadot = -2*adotoa*sigma-dgs/k+etak

            !Once know slip, recompute qgdot, pig, pigdot
            qgdot = k*(clxg/4._dl-pig/2._dl) +opacity*slip

            pig = 32._dl/45/opacity*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity*11._dl/6._dl/opacity**2)) &
                + (32._dl/45._dl/opacity**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

            pigdot = -(32._dl/45._dl)*(dopacity/opacity**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
                dopacity*11._dl/6._dl/opacity**2 ) &
                + (32._dl/45._dl/opacity)*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
                *(dopacity/opacity**2))

            EV%pigdot = pigdot
        end if

        !  Use tight-coupling approximation for vb
        !  zeroth order approximation to vbdot + the pig term
        vbdot=(-adotoa*vb+cs2*k*clxb  &
            +k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip

        EV%pig = pig
    else
        vbdot=-adotoa*vb+k*delta_p-photbar*opacity*(4._dl/3*vb-qg)
    end if

    ayprime(5)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
        ayprime(EV%g_ix+1)=qgdot

        !  Use explicit equations for photon moments if appropriate
        if (.not. EV%tightcoupling) then
            E2=ay(EV%polind+2)
            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
            ix= EV%g_ix+2
            if (EV%lmaxg>2) then
                pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity*(pig - polter) &
                    +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
                do  l=3,EV%lmaxg-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity*ay(ix)
                end do
                ix=ix+1
                !  Truncate the photon moment expansion
                ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity*ay(ix)
            else !closed case
                pigdot=EV%denlk(2)*qg-opacity*(pig - polter) +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
            endif
            !  Polarization
            !l=2
            ix=EV%polind+2
            if (EV%lmaxgpol>2) then
                ayprime(ix) = -opacity*(ay(ix) - polter) - k/3._dl*ay(ix+1)
                do l=3,EV%lmaxgpol-1
                    ix=ix+1
                    ayprime(ix)=-opacity*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                end do
                ix=ix+1
                !truncate
                ayprime(ix)=-opacity*ay(ix) + &
                    k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
            else !closed case
                ayprime(ix) = -opacity*(ay(ix) - polter)
            endif
        end if
    end if

    if (.not. EV%no_nu_multpoles) then
        !  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%r_ix)=clxrdot
        qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
        ayprime(EV%r_ix+1)=qrdot
        if (EV%high_ktau_neutrino_approx) then
            !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
            !Method from arXiv:1104.2933
            !                if (.not. EV%TightCoupling) then
            !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            !                 adotdota=(adotoa*adotoa-gpres)/2
            !                end if
            !                ddz=(2*adotoa**2 - adotdota)*z  &
            !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*(grhoc_t*clxc+grhob_t*clxb) ) &
            !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + grhoc_t*clxcdot + grhob_t*clxbdot )
            !                dz= -adotoa*z - 0.5_dl*dgrho/k
            !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
            pirdot= -3*pir*cothxor - clxrdot
            ayprime(EV%r_ix+2)=pirdot

            !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
            !                ayprime(EV%lmaxg+9)=pirdot
            !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
            !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
            !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
        else
            ix=EV%r_ix+2
            if (EV%lmaxnr>2) then
                pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
                do l=3,EV%lmaxnr-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
                end do
                !  Truncate the neutrino expansion
                ix=ix+1
                ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
            else
                pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
            end if
        end if
    end if ! no_nu_multpoles

    if (Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            Delta_TCMB = clxg/4
            xe = Recombination_xe(a)
            Trad = CP%TCMB/a

            !Matter temperature
            !Recfast_CT = (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc [a_R = radiation constant]
            ayprime(EV%Tg_ix) = -2*k*(z+vb)/3 - a*  Compton_CT * (Trad**4) * xe / (1._dl+xe+fHe) * &
                ((1- Trad/Tmat)*(Delta_TCMB*4 + Delta_xe/(1+xe/(1+fHe))) + Trad/Tmat*(Delta_Tm - Delta_TCMB)  )

            if (Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) =  Get21cm_dTs(a,clxb,ay(EV%Ts_ix),Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
            end if
        else
            if (Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) = -k*(4._dl/3._dl*z+qg)/4  !Assume follows Delta_TM which follows clxg
            end if
        end if
    end if

    if (Evolve_delta_xe .and. .not. EV%saha) then
        ayprime(EV%xe_ix) = dDeltaxe_dtau(a, Delta_xe,clxb, Delta_Tm, k*z/3,k*vb)
    end if

    if (Do21cm) then
        if (a > Do21cm_mina) then
            if (line_reionization) then
                lineoff = EV%reion_line_ix+1
                lineoffpol = lineoff+EV%lmaxline-1

                if (tau> tau_start_redshiftwindows) then
                    !Multipoles of 21cm

                    polter_line = ay(lineoff+2)/10+9._dl/15*ay(lineoffpol+2)

                    call interp_window(Redshift_W(1),tau,wing_t,wing2_t,winv_t)

                    delta_source2 = Get21cm_source2(a,clxb,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe,k*(z+vb)/adotoa/3)


                    !Drop some small terms since mulipoles only enter into reionzation anyway
                    !monopole
                    ayprime(lineoff) = -k*ay(lineoff+1) +  wing_t * clxb + wing2_t*delta_source2 + k*z/3*winV_t


                    !dipole
                    ayprime(lineoff+1)= EV%denlk(1)*ay(lineoff)-EV%denlk2(1)*ay(lineoff+2) - opacity*ay(lineoff+1) &
                        -wing2_t * ( qg/4 - vb/3)   ! vb/3*WinV_t)

                    !quadrupole
                    ayprime(lineoff+2)= EV%denlk(2)*ay(lineoff+1)-EV%denlk2(2)*ay(lineoff+3) &
                        +opacity*(polter_line -ay(lineoff+2) ) -   2._dl/15*k*sigma*winV_t &
                        - wing2_t * ay(6+2)/4

                    do  l=3,EV%lmaxline-1
                        ayprime(lineoff+l)=EV%denlk(l)*ay(lineoff+l-1)-EV%denlk2(l)*ay(lineoff+l+1)-opacity*ay(lineoff+l) &
                            - wing2_t * ay(6+l)/4
                    end do
                    !truncate
                    ayprime(lineoff+EV%lmaxline)=k*ay(lineoff+EV%lmaxline-1)-(EV%lmaxline+1)*cothxor*ay(lineoff+EV%lmaxline)  &
                        -opacity*ay(lineoff+EV%lmaxline) - wing2_t * ay(6+EV%lmaxline)/4

                    !  21cm Polarization
                    !l=2
                    ayprime(lineoffpol+2) = -opacity*(ay(lineoffpol+2) - polter_line) - k/3._dl*ay(lineoffpol+3)
                    !and the rest
                    do l=3,EV%lmaxline-1
                        ayprime(lineoffpol+l)=-opacity*ay(lineoffpol+l) + EV%denlk(l)*ay(lineoffpol+l-1) -&
                            EV%polfack(l)*ay(lineoffpol+l+1)
                    end do

                    !truncate
                    ayprime(lineoffpol+EV%lmaxline)=-opacity*ay(lineoffpol+EV%lmaxline) + &
                        k*EV%poltruncfac*ay(lineoffpol+EV%lmaxline-1)-(EV%lmaxline+3)*cothxor*ay(lineoffpol+EV%lmaxline)
                else
                    ayprime(lineoff:lineoffpol+EV%lmaxline)=0
                end if
            end if
        end if
    end if

    if (associated(EV%OutputTransfer)) then
        EV%OutputTransfer(Transfer_kh) = k/(CP%h0/100._dl)
        EV%OutputTransfer(Transfer_cdm) = clxc
        EV%OutputTransfer(Transfer_b) = clxb
        EV%OutputTransfer(Transfer_g) = clxg
        EV%OutputTransfer(Transfer_r) = clxr
        clxnu_all=0
        dgpi  = grhor_t*pir + grhog_t*pig
        if (CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, clxnu_all =clxnu_all, dgpi= dgpi)
        end if
        EV%OutputTransfer(Transfer_nu) = clxnu_all
        EV%OutputTransfer(Transfer_tot) =  dgrho_matter/grho_matter !includes neutrinos
        EV%OutputTransfer(Transfer_nonu) = (grhob_t*clxb+grhoc_t*clxc)/(grhob_t + grhoc_t)
        EV%OutputTransfer(Transfer_tot_de) =  dgrho/grho_matter
        !Transfer_Weyl is k^2Phi, where Phi is the Weyl potential
        EV%OutputTransfer(Transfer_Weyl) = -(dgrho +3*dgq*adotoa/k)/(EV%Kf(1)*2) - dgpi/2
        EV%OutputTransfer(Transfer_Newt_vel_cdm)=  -k*sigma/adotoa
        EV%OutputTransfer(Transfer_Newt_vel_baryon) = -k*(vb + sigma)/adotoa
        EV%OutputTransfer(Transfer_vel_baryon_cdm) = vb

        if (do21cm) then
            Tspin = Recombination_Ts(a)
            xe = Recombination_xe(a)

            tau_eps = a*line21_const*NNow/a**3/adotoa/Tspin/1000
            delta_source2 = Get21cm_source2(a,clxb,clxg/4,Delta_Tm,Delta_xe,Tmat,&
                CP%TCMB/a,xe,k*(z+vb)/adotoa/3)
            tau_fac = tau_eps/(exp(tau_eps)-1)
            EV%OutputTransfer(Transfer_monopole) = ( clxb + Trad/(Tspin-Trad)*delta_source2 ) /k2 &
                + (tau_fac-1)*(clxb - (delta_source2 + clxg/4)  ) / k2

            EV%OutputTransfer(Transfer_vnewt) = tau_fac*k*(vb+sigma)/adotoa/k2
            EV%OutputTransfer(Transfer_Tmat) =  delta_TM/k2
            if (use_mK) then
                Tb = (1-exp(-tau_eps))*a*(Tspin-Trad)*1000

                EV%OutputTransfer(Transfer_monopole) = EV%OutputTransfer(Transfer_monopole)*Tb
                EV%OutputTransfer(Transfer_vnewt) = EV%OutputTransfer(Transfer_vnewt)*Tb
                EV%OutputTransfer(Transfer_Tmat) = EV%OutputTransfer(Transfer_Tmat)*Tb
            end if
        end if
    end if

    !  Massive neutrino equations of motion.
    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
        if (EV%MassiveNuApprox(nu_i)) then
            !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
            !see astro-ph/0203507
            G11_t=EV%G11(nu_i)/a/a2
            G30_t=EV%G30(nu_i)/a/a2
            off_ix = EV%nu_ix(nu_i)
            w=wnu_arr(nu_i)
            ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
            ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
            ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
            ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
        else
            ind=EV%nu_ix(nu_i)
            do i=1,EV%nq(nu_i)
                q=nu_q(i)
                aq=a*nu_masses(nu_i)/q
                v=1._dl/sqrt(1._dl+aq*aq)

                ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
                ind=ind+1
                ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                ind=ind+1
                if (EV%lmaxnu_tau(nu_i)==2) then
                    ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                else
                    ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                        +k*8._dl/15._dl*sigma
                    do l=3,EV%lmaxnu_tau(nu_i)-1
                        ind=ind+1
                        ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                    end do
                    !  Truncate moment expansion.
                    ind = ind+1
                    ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                end if
                ind = ind+1
            end do
        end if
    end do

    if (EV%has_nu_relativistic) then
        ind=EV%nu_pert_ix
        ayprime(ind)=+k*a2*qr -k*ay(ind+1)
        ind2= EV%r_ix
        do l=1,EV%lmaxnu_pert-1
            ind=ind+1
            ind2=ind2+1
            ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
        end do
        ind=ind+1
        ind2=ind2+1
        ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
    end if

    end subroutine derivs



    subroutine derivsv(EV,n,tau,yv,yvprime)
    !  Evaluate the time derivatives of the vector perturbations, flat case
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l
    real(dl), target ::  yv(n),yvprime(n)
    real(dl) ep,tau,grho,rhopi,cs2,opacity,gpres
    logical finished_tightcoupling
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) sigma, qg,pig, qr, vb, rhoq, vbdot, photbar, pb43
    real(dl) k,k2,a,a2, adotdota
    real(dl) pir,adotoa
    real(dl) w_dark_energy_t

    k2=EV%k2_buf
    k=EV%k_buf

    !E and B start at l=2. Set up pointers accordingly to fill in y arrays
    E => yv(EV%lmaxv+3:)
    Eprime=> yvprime(EV%lmaxv+3:)
    B => E(EV%lmaxpolv:)
    Bprime => Eprime(EV%lmaxpolv:)
    neutprime => Bprime(EV%lmaxpolv+1:)
    neut => B(EV%lmaxpolv+1:)

    a=yv(1)

    sigma=yv(2)

    a2=a*a

    !  Get sound speed and opacity, and see if should use tight-coupling

    call thermo(tau,cs2,opacity)
    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if

    finished_tightcoupling = &
        ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4))


    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    call CP%DarkEnergy%BackgroundDensityAndPressure(a, grhov_t, w_dark_energy_t)

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
    gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_dark_energy_t

    adotoa=sqrt(grho/3._dl)
    adotdota=(adotoa*adotoa-gpres)/2

    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    yvprime(1)=adotoa*a

    vb = yv(3)
    qg = yv(4)
    qr = neut(1)

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    rhoq=grhob_t*vb+grhog_t*qg+grhor_t*qr
    !  sigma = 2*rhoq/k**2
    !for non-large k this expression for sigma is unstable at early times
    !so propagate sigma equation separately (near total cancellation in rhoq)
    ! print *,yv(2),2*rhoq/k**2

    if (finished_tightcoupling) then
        !  Use explicit equations:

        pig = yv(5)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        vbdot = -adotoa*vb-photbar*opacity*(4._dl/3*vb-qg) - 0.5_dl*k*photbar*Magnetic

        !  Equation for the photon heat flux stress

        yvprime(4)=-0.5_dl*k*pig + opacity*(4._dl/3*vb-qg)

        !  Equation for the photon anisotropic stress
        yvprime(5)=k*(2._dl/5*qg -8/15._dl*yv(6))+8._dl/15._dl*k*sigma  &
            -opacity*(pig - polter)
        ! And for the moments
        do  l=3,EV%lmaxv-1
            yvprime(l+3)=k*denl(l)*l*(yv(l+2)-   &
                vecfac(l)*yv(l+4))-opacity*yv(l+3)
        end do
        !  Truncate the hierarchy
        yvprime(EV%lmaxv+3)=k*EV%lmaxv/(EV%lmaxv-1._dl)*yv(EV%lmaxv+2)- &
            (EV%lmaxv+2._dl)*yv(EV%lmaxv+3)/tau-opacity*yv(EV%lmaxv+3)

        !E equations

        Eprime(2) = - opacity*(E(2) - polter) + k*(1/3._dl*B(2) - 8._dl/27._dl*E(3))
        do l=3,EV%lmaxpolv-1
            Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
                vecfacpol(l)*E(l+1)) + 2._dl/(l*(l+1))*B(l))
        end do
        !truncate
        Eprime(EV%lmaxpolv)=0._dl

        !B-bar equations

        do l=2,EV%lmaxpolv-1
            Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
                vecfacpol(l)*B(l+1)) - 2._dl/(l*(l+1))*E(l))
        end do
        !truncate
        Bprime(EV%lmaxpolv)=0._dl
    else
        !Tight coupling expansion results

        pig = 32._dl/45._dl*k/opacity*(vb + sigma)

        EV%pig = pig

        vbdot=(-adotoa*vb  -3._dl/8*pb43*k*Magnetic  -3._dl/8*k*pb43*pig &
            - pb43/(1+pb43)/opacity*(0.75_dl*k*adotoa*pb43**2/(pb43+1)*Magnetic + vb*&
            ( 2*pb43*adotoa**2/(1+pb43) + adotdota)) &
            )/(1+pb43)

        !  Equation for the photon heat flux
        ! Get drag from vbdot expression
        yvprime(4)=-0.5_dl*k*pig - &
            (vbdot+adotoa*vb)/photbar - 0.5_dl*k*Magnetic

        !  Set the derivatives to zero
        yvprime(5:n)=0._dl
        yv(5)=pig
        E(2)=  pig/4
    endif

    yvprime(3) = vbdot

    !  Neutrino equations:

    !  Massless neutrino anisotropic stress
    pir=neut(2)
    neutprime(1)= -0.5_dl*k*pir
    neutprime(2)=2._dl/5*k*qr -8._dl/15._dl*k*neut(3)+ 8._dl/15._dl*k*sigma
    !  And for the moments
    do  l=3,EV%lmaxnrv-1
        neutprime(l)=k*denl(l)*l*(neut(l-1)- vecfac(l)*neut(l+1))
    end do

    !  Truncate the hierarchy
    neutprime(EV%lmaxnrv)=k*EV%lmaxnrv/(EV%lmaxnrv-1._dl)*neut(EV%lmaxnrv-1)-  &
        (EV%lmaxnrv+2._dl)*neut(EV%lmaxnrv)/tau


    !  Get the propagation equation for the shear

    rhopi=grhog_t*pig+grhor_t*pir+ grhog_t*Magnetic

    yvprime(2)=-2*adotoa*sigma -rhopi/k

    end subroutine derivsv



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine derivst(EV,n,tau,ayt,aytprime)
    !  Evaluate the time derivatives of the tensor perturbations.
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l,i,ind, nu_i
    real(dl), target ::  ayt(n),aytprime(n)
    real(dl) tau,grho,rhopi,cs2,opacity,pirdt
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl) q,aq,v
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) Hchi,pinu, pig
    real(dl) k,k2,a,a2
    real(dl) pir, adotoa, rhonu, shear

    real(dl) cothxor

    k2=EV%k2_buf
    k= EV%k_buf

    a=ayt(1)

    Hchi=ayt(2)

    shear=ayt(3)

    a2=a*a

    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    call CP%DarkEnergy%BackgroundDensityAndPressure(a, grhov_t)

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t

    !Do massive neutrinos
    if (CP%Num_Nu_Massive >0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grho=grho+grhormass(nu_i)*rhonu/a2
        end do
    end if

    if (CP%flat) then
        cothxor=1._dl/tau
        adotoa=sqrt(grho/3._dl)
    else
        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
        adotoa=sqrt((grho+grhok)/3._dl)
    end if

    aytprime(1)=adotoa*a

    call thermo(tau,cs2,opacity)

    if (.not. EV%TensTightCoupling) then
        !  Don't use tight coupling approx - use explicit equations:
        !  Equation for the photon anisotropic stress


        !E and B start at l=2. Set up pointers accordingly to fill in ayt arrays
        E => ayt(EV%E_ix+1:)
        B => ayt(EV%B_ix+1:)
        Eprime=> aytprime(EV%E_ix+1:)
        Bprime => aytprime(EV%B_ix+1:)

        ind = EV%g_ix+2

        !  Photon anisotropic stress
        pig=ayt(ind)
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        if (EV%lmaxt > 2) then
            aytprime(ind)=-EV%denlkt(2,2)*ayt(ind+1)+k*8._dl/15._dl*shear  &
                -opacity*(pig - polter)

            do l=3, EV%lmaxt -1
                ind = ind+1
                aytprime(ind)=EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1)-opacity*ayt(ind)
            end do

            !Truncate the hierarchy
            ind=ind+1
            aytprime(ind)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(ind-1)- &
                (EV%lmaxt+3._dl)*cothxor*ayt(ind)-opacity*ayt(ind)

            !E and B-bar equations

            Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2) - &
                EV%denlkt(3,2)*E(3)

            do l=3, EV%lmaxpolt-1
                Eprime(l) =(EV%denlkt(1,L)*E(l-1)-EV%denlkt(3,L)*E(l+1) + EV%denlkt(4,L)*B(l)) &
                    -opacity*E(l)
            end do
            l= EV%lmaxpolt
            !truncate: difficult, but setting l+1 to zero seems to work OK
            Eprime(l) = (EV%denlkt(1,L)*E(l-1) + EV%denlkt(4,L)*B(l)) -opacity*E(l)

            Bprime(2) =-EV%denlkt(3,2)*B(3) - EV%denlkt(4,2)*E(2)  -opacity*B(2)
            do l=3, EV%lmaxpolt-1
                Bprime(l) =(EV%denlkt(1,L)*B(l-1) -EV%denlkt(3,L)*B(l+1) - EV%denlkt(4,L)*E(l)) &
                    -opacity*B(l)
            end do
            l=EV%lmaxpolt
            !truncate
            Bprime(l) =(EV%denlkt(1,L)*B(l-1) - EV%denlkt(4,L)*E(l))  -opacity*B(l)

        else !lmax=2

            aytprime(ind)=k*8._dl/15._dl*shear-opacity*(pig - polter)
            Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2)
            Bprime(2) = - EV%denlkt(4,2)*E(2)  -opacity*B(2)
        end if

    else  !Tight coupling
        pig = 32._dl/45._dl*k/opacity*shear
    endif

    rhopi=grhog_t*pig


    !  Neutrino equations:
    !  Anisotropic stress
    if (DoTensorNeutrinos) then
        neutprime => aytprime(EV%r_ix+1:)
        neut => ayt(EV%r_ix+1:)

        !  Massless neutrino anisotropic stress
        pir=neut(2)

        rhopi=rhopi+grhor_t*pir

        if (EV%lmaxnrt>2) then
            pirdt=-EV%denlkt(2,2)*neut(3) + 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
            !  And for the moments
            do  l=3, EV%lmaxnrt-1
                neutprime(l)= EV%denlkt(1,L)*neut(l-1) -EV%denlkt(2,L)*neut(l+1)
            end do

            !  Truncate the hierarchy
            neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
                (EV%lmaxnrt+3._dl)*cothxor*neut(EV%lmaxnrt)
        else
            pirdt= 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
        end if

        !  Massive neutrino equations of motion and contributions to anisotropic stress.
        if (CP%Num_Nu_massive > 0) then
            do nu_i=1,CP%Nu_mass_eigenstates
                if (.not. EV%EvolveTensorMassiveNu(nu_i)) then
                    rhopi=rhopi+ grhormass(nu_i)/a2*pir !- good approx, note no rhonu weighting
                else
                    ind=EV%nu_ix(nu_i)+2

                    pinu= Nu_pi(EV, ayt, a, nu_i)
                    rhopi=rhopi+ grhormass(nu_i)/a2*pinu

                    do i=1,nqmax
                        q=nu_q(i)
                        aq=a*nu_masses(nu_i)/q
                        v=1._dl/sqrt(1._dl+aq*aq)
                        if (EV%lmaxnut>2) then
                            aytprime(ind)=-v*EV%denlkt(2,2)*ayt(ind+1)+8._dl/15._dl*k*shear
                            do l=3,EV%lmaxnut-1
                                ind=ind+1
                                aytprime(ind)=v*(EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1))
                            end do
                            ind = ind+1
                            !Truncate moment expansion.
                            aytprime(ind)=k*v*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-1)-(EV%lmaxnut+3)*cothxor*ayt(ind)
                        else
                            aytprime(ind)=8._dl/15._dl*k*shear
                        end if
                        ind=ind+1
                    end do
                end if
            end do
        end if
    end if

    !  Get the propagation equation for the shear

    if (CP%flat) then
        aytprime(3)=-2*adotoa*shear+k*Hchi-rhopi/k
    else
        aytprime(3)=-2*adotoa*shear+k*Hchi*(1+2*CP%curv/k2)-rhopi/k
    endif

    aytprime(2)=-k*shear

    end subroutine derivst



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    end module GaugeInterface

    function isTmNeeded()
    use GaugeInterface
    logical :: isTmNeeded

    isTmNeeded = Evolve_baryon_cs .or. Evolve_delta_xe

    end function isTmNeeded
