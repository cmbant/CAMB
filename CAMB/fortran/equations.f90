    ! Equations module for background and ! To avoid circular module issues, some things are not part of module

    ! Background evolution, return d tau/ d a, where tau is the conformal time
    function dtauda(this,a)
    use results
    use DarkEnergyInterface
    implicit none
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) :: dtauda, grhoa2, grhov_t

    call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)

    !  8*pi*G*rho*a**4.
    grhoa2 = this%grho_no_de(a) +  grhov_t * a**2
    if (grhoa2 <= 0) then
        call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
        dtauda = 0
    else
        dtauda = sqrt(3 / grhoa2)
    end if

    end function dtauda

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !Gauge-dependent perturbation equations

    module GaugeInterface
    use precision
    use results
    use MassiveNu
    use DarkEnergyInterface
    use Transfer
    implicit none
    public

    !Description of this file. Change if you make modifications.
    character(LEN=*), parameter :: Eqns_name = 'cdm_gauge'

    logical, parameter :: plot_evolve = .false. !for outputing time evolution

    integer, parameter :: basic_num_eqns = 4
    integer, parameter :: ix_etak=1, ix_clxc=2, ix_clxb=3, ix_vb=4 !Scalar array indices for each quantity
    integer, parameter :: ixt_H = 1, ixt_shear = 2 !tensor indices

    logical :: DoTensorNeutrinos = .true.

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
        logical :: is_cosmological_constant

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

        real(dl) pig, pigdot
        real(dl) poltruncfac

        logical no_nu_multpoles, no_phot_multpoles
        integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
        logical nu_nonrelativistic(max_nu)

        real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
        real(dl) Kf(max_l_evolve)

        integer E_ix, B_ix !tensor polarizatisdon indices
        real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)

        logical :: saha !still high x_e
        logical :: evolve_TM !\delta T_g evolved separately
        logical :: evolve_baryon_cs !evolving sound speed rather than using background approx

        !Workaround for ifort, gives class pointer to avoid creating temps and huge slow down
        class(TThermoData), pointer :: ThermoData => null()

        real, pointer :: OutputTransfer(:) => null()
        real(dl), pointer :: OutputSources(:) => null()
        real(dl), pointer :: CustomSources(:) => null()
        integer :: OutputStep = 0

    end type EvolutionVars

    ABSTRACT INTERFACE
    SUBROUTINE TSource_func(sources, tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxr, clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis, Kf, f_K)
    use precision
    real(dl), intent(out) :: sources(:)
    real(dl), intent(in) :: tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc, clxr, clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3), &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis
    REAL(dl), intent(in) :: Kf(*)
    real(dl), external :: f_K
    END SUBROUTINE TSource_func
    END INTERFACE

    class(CAMBdata), pointer :: State !Current state.
    !(Due to ifort bug State needs to be class pointer to avoid making temp class pointers when calling functions)
    type(CAMBParams), pointer :: CP   !Current parameters (part of state)

    !precalculated arrays
    real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

    real(dl), parameter :: ep0=1.0d-2
    integer, parameter :: lmaxnu_high_ktau=4 !Jan2015, increased from 3 to fix mpk for mnu~6eV

    real(dl) epsw
    real(dl), allocatable :: nu_tau_notmassless(:,:)
    real(dl) nu_tau_nonrelativistic(max_nu), nu_tau_massive(max_nu)

    procedure(state_function), private :: dtauda
    contains

    subroutine SetActiveState(P)
    class(CAMBdata), target :: P

    State => P
    CP => P%CP

    end subroutine SetActiveState


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

    function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x

    f_K = State%curvature_radius*State%rofChi(x/State%curvature_radius)

    end function f_K

    function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
        next_nq=1
    else
        q = int(State%NuPerturbations%nu_q(nq))
        if (q>=10) then
            next_nq = State%NuPerturbations%nqmax
        else
            next_nq = nq+1
        end if
    end if

    end function next_nu_nq

    recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
    use Recombination, only : CB1
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity, dopacity
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) noSwitch, smallTime
    !Sources
    real(dl) tau_switch_saha, Delta_TM, xe,a,tau_switch_evolve_TM

    noSwitch= State%tau0+1
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
    if (CP%Evolve_delta_xe .and. EV%saha)  tau_switch_saha = EV%ThermoData%recombination_saha_tau
    tau_switch_evolve_TM=noSwitch
    if (EV%Evolve_baryon_cs .and. .not. EV%Evolve_tm) tau_switch_evolve_TM = EV%ThermoData%recombination_Tgas_tau

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
        tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
            else if (.not. EV%nu_nonrelativistic(nu_i)) then
                tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
            else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
                tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
            end if
        end do
    end if

    if (CP%DoLateRadTruncation) then
        if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
            tau_switch_no_nu_multpoles= &
            max(15/EV%k_buf*CP%Accuracy%AccuracyBoost,min(State%taurend,EV%ThermoData%matter_verydom_tau))

        if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*CP%Accuracy%AccuracyBoost)) &
            tau_switch_no_phot_multpoles =max(15/EV%k_buf,State%taurend)*CP%Accuracy%AccuracyBoost
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
            call EV%ThermoData%Values(tau,a, cs2,opacity,dopacity)

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
            EVout%nq(1:CP%Nu_mass_eigenstates) = State%NuPerturbations%nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch == tau_switch_nu_massless) then
            !Mass starts to become important, start evolving next momentum mode
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                    if (next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                        EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                        call SetupScalarArrayIndices(EVout)
                        call CopyScalarVariableArray(y,yout, EV, EVout)
                        EV=EVout
                        y=yout
                        exit
                    endif
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
            call EV%ThermoData%Values(tau,a, cs2,opacity)
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%MassiveNuApprox(nu_i) .and.  next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                    call SwitchToMassiveNuApprox(EV, a, y, nu_i)
                    exit
                end if
            end do
        else if (next_switch==tau_switch_no_nu_multpoles) then
            !Turn off neutrino hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_nu_multpoles=.true.
            EVOut%nq(1:CP%Nu_mass_eigenstates ) = State%NuPerturbations%nqmax
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
            call EV%ThermoData%Values(tau,a, cs2,opacity)
            y=yout
            EV=EVout
            Delta_Tm = y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
            xe= CP%Recomb%x_e(a)
            y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(ix_clxb) + (3./2+  CB1/(CP%TCMB/a))*Delta_TM)
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
    type(EvolutionVars) EV, EVOut
    real(dl) c(24),w(EV%nvart,9), y(EV%nvart),yout(EV%nvart), tol1, tau, tauend
    integer ind
    real(dl) opacity, cs2, a

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
        call EV%ThermoData%Values(tau,a,cs2,opacity)
        y(EV%g_ix+2)= 32._dl/45._dl*EV%k_buf/opacity*y(ixt_shear)
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
        t = State%DeltaTime(a1,1.01_dl, tol)
    else
        t = State%DeltaTime(a1,a2, tol)
    end if
    end function DeltaTimeMaxed

    subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/State%tau0

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

    if (CP%Nu_Mass_eigenstates>0) then
        associate(nqmax => State%NuPerturbations%nqmax, nu_q => State%NuPerturbations%nu_q)
            if (allocated(nu_tau_notmassless)) deallocate(nu_tau_notmassless)
            allocate(nu_tau_notmassless(nqmax,max_nu))
            do nu_i=1, CP%Nu_Mass_eigenstates
                nu_mass = max(0.1_dl,State%nu_masses(nu_i))
                a_mass =  1.e-1_dl/nu_mass/CP%Accuracy%lAccuracyBoost
                time=State%DeltaTime(0._dl,State%NuPerturbations%nu_q(1)*a_mass)
                nu_tau_notmassless(1, nu_i) = time
                do j=2,nqmax
                    !times when each momentum mode becomes signficantly nonrelativistic
                    time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
                    nu_tau_notmassless(j, nu_i) = time
                end do

                a_nonrel =  2.5d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
                a_massive =  17.d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
            end do
        end associate
    end if

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
    if (CP%Evolve_delta_xe) then
        if (.not. EV%saha) then
            EV%xe_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
    end if

    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            EV%Tg_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%reion_line_ix = neq+1
            neq=neq+ EV%lmaxline+1 +  EV%lmaxline-1
            maxeq=maxeq+EV%lmaxline+1 +  EV%lmaxline-1
        end if
    end if

    if (CP%Evolve_delta_Ts) then
        EV%Ts_ix = neq+1
        neq=neq+1
        maxeq=maxeq+1
    end if

    !Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
        EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=State%NuPerturbations%nqmax)
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
                EV%lmaxnu_tau(nu_i) = int(lmaxnu_high_ktau *CP%Accuracy%lAccuracyBoost)
                if (CP%Transfer%accurate_massive_neutrinos) EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i) *3
            else
                EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i) &
                    *CP%Accuracy%lAccuracyBoost),EV%lmaxnu),3)
                !!!Feb13tweak
                if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)= &
                    min(EV%lmaxnu_tau(nu_i),nint(4*CP%Accuracy%lAccuracyBoost))
            end if
            if (State%nu_masses(nu_i) > 5000 .and. CP%Transfer%high_precision) &
                EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i)*2 !megadamping
            EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

            EV%nu_ix(nu_i)=neq+1
            if (EV%MassiveNuApprox(nu_i)) then
                neq = neq+4
            else
                neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
            endif
            maxeq = maxeq + State%NuPerturbations%nqmax*(EV%lmaxnu+1)
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
                    q=State%NuPerturbations%nu_q(i)
                    pert_scale=(State%nu_masses(nu_i)/q)**2/2
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
    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM .and. EVout%Evolve_TM) yout(EVOut%Tg_ix) = y(EV%Tg_ix)
        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            yout(EVOut%reion_line_ix:EVOut%reion_line_ix+EVout%lmaxline +  EVout%lmaxline-1) = &
                y(EV%reion_line_ix:EV%reion_line_ix+EV%lmaxline +  EV%lmaxline-1)
        end if
    end if
    if (CP%Evolve_delta_Ts) then
        yout(EVOut%Ts_ix) = y(EV%Ts_ix)
    end if

    end subroutine CopyScalarVariableArray


    subroutine SetupTensorArrayIndices(EV, maxeq)
    type(EvolutionVars) EV
    integer nu_i, neq
    integer, optional, intent(out) :: maxeq
    neq=2
    EV%g_ix = neq-1 !EV%g_ix+2 is quadrupole
    if (.not. EV%TensTightCoupling) then
        EV%E_ix = EV%g_ix + (EV%lmaxt-1)
        EV%B_ix = EV%E_ix + (EV%lmaxpolt-1)
        neq = neq+ (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    if (present(maxeq)) then
        maxeq = 2 + (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    EV%r_ix = neq -1
    if (DoTensorNeutrinos) then
        neq = neq + EV%lmaxnrt-1
        if (present(maxeq)) maxeq = maxeq+EV%lmaxnrt-1
        if (CP%Num_Nu_massive /= 0 ) then
            do nu_i=1, CP%nu_mass_eigenstates
                EV%EvolveTensorMassiveNu(nu_i) = &
                    nu_tau_nonrelativistic(nu_i) < 0.8*State%tau_maxvis*CP%Accuracy%AccuracyBoost
                if (EV%EvolveTensorMassiveNu(nu_i)) then
                    EV%nu_ix(nu_i)=neq-1
                    neq = neq+ State%NuPerturbations%nqmax*(EV%lmaxnut-1)
                    if (present(maxeq)) maxeq = maxeq + State%NuPerturbations%nqmax*(EV%lmaxnut-1)
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
    yout(1:2) = y(1:2)
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
                do i=1,State%NuPerturbations%nqmax
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
        max_nu_mass = maxval(State%nu_masses(1:CP%Nu_mass_eigenstates))
        do nu_i = 1, CP%Nu_mass_eigenstates
            !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
            q_rel=0
            do j=1, State%NuPerturbations%nqmax
                !two different q's here EV%q ~k
                if (State%NuPerturbations%nu_q(j) > State%nu_masses(nu_i)*State%adotrad/EV%q) exit
                q_rel = q_rel + 1
            end do

            if (q_rel>= State%NuPerturbations%nqmax-2 .or. CP%WantTensors) then
                EV%nq(nu_i)=State%NuPerturbations%nqmax
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
            EV%lmaxnu=nint(25*CP%Accuracy%lAccuracyBoost)
        else
            EV%lmaxnu=max(3,nint(10*CP%Accuracy%lAccuracyBoost))
            if (max_nu_mass>700) EV%lmaxnu=max(3,nint(15*CP%Accuracy%lAccuracyBoost)) !Feb13 tweak
        endif
    end if

    if (State%closed) then
        EV%FirstZerolForBeta = nint(EV%q*State%curvature_radius)
    else
        EV%FirstZerolForBeta= 100000 !a large number
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

        if (CP%Accuracy%AccuratePolarization) then
            EV%lmaxg  = max(nint(11*CP%Accuracy%lAccuracyBoost),3)
        else
            EV%lmaxg  = max(nint(8*CP%Accuracy%lAccuracyBoost),3)
        end if
        EV%lmaxnr = max(nint(14*CP%Accuracy%lAccuracyBoost),3)
        if (max_nu_mass>700) EV%lmaxnr = max(nint(32*CP%Accuracy%lAccuracyBoost),3) !Feb13 tweak

        EV%lmaxgpol = EV%lmaxg
        if (.not.CP%Accuracy%AccuratePolarization) EV%lmaxgpol=max(nint(4*CP%Accuracy%lAccuracyBoost),3)

        if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            if (CP%Accuracy%AccuratePolarization) scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*CP%Accuracy%lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*CP%Accuracy%lAccuracyBoost))
            if (EV%lmaxnr < EV%lmaxnu) then
                ! Nov 2020 change following Pavel Motloch report
                EV%lmaxnr = EV%lmaxnu
                !EV%lmaxnu = min(EV%lmaxnu, EV%lmaxnr) ! may be better but have not tested and makes small result changes
            endif
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*CP%Accuracy%lAccuracyBoost))
            !Sources
            if (CP%SourceTerms%line_phot_quadrupole) then
                EV%lmaxg=EV%lmaxg*8
                EV%lmaxgpol=EV%lmaxgpol*4
            elseif (CP%Accuracy%AccurateReionization) then
                EV%lmaxg=EV%lmaxg*4
                EV%lmaxgpol=EV%lmaxgpol*2
            end if
        end if

        if (EV%TransferOnly) then
            EV%lmaxgpol = min(EV%lmaxgpol,nint(5*CP%Accuracy%lAccuracyBoost))
            EV%lmaxg = min(EV%lmaxg,nint(6*CP%Accuracy%lAccuracyBoost))
        end if
        if (CP%Transfer%high_precision .or. CP%Do21cm) then
            EV%lmaxnr=max(nint(45*CP%Accuracy%lAccuracyBoost),3)
            if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
                EV%lmaxg=max(EV%lmaxg,10)
            end if
        end if

        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%lmaxg =  EV%lmaxg*8
            EV%lmaxgpol = EV%lmaxgpol*3
        end if

        EV%Evolve_baryon_cs = CP%Do21cm .or.CP%Evolve_delta_xe .or. CP%Evolve_delta_Ts

        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%lmaxline  = EV%lmaxg
        end if

        if (State%closed) then
            EV%lmaxnu=min(EV%lmaxnu, EV%FirstZerolForBeta-1)
            EV%lmaxnr=min(EV%lmaxnr, EV%FirstZerolForBeta-1)
            EV%lmaxg=min(EV%lmaxg, EV%FirstZerolForBeta-1)
            EV%lmaxgpol=min(EV%lmaxgpol, EV%FirstZerolForBeta-1)
        end if

        EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
        EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
        if (EV%MaxlNeeded > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupScalarArrayIndices(EV,EV%nvar)
        if (State%closed) EV%nvar=EV%nvar+1 !so can reference lmax+1 with zero coefficient
        EV%lmaxt=0
    else
        EV%nvar=0
    end if

    if (CP%WantTensors) then
        EV%TensTightCoupling = .true.
        EV%lmaxt=max(3,nint(8*CP%Accuracy%lAccuracyBoost))
        EV%lmaxpolt = max(3,nint(4*CP%Accuracy%lAccuracyBoost))
        ! if (EV%q < 1e-3) EV%lmaxpolt=EV%lmaxpolt+1
        if (DoTensorNeutrinos) then
            EV%lmaxnrt=nint(6*CP%Accuracy%lAccuracyBoost)
            EV%lmaxnut=EV%lmaxnrt
        else
            EV%lmaxnut=0
            EV%lmaxnrt=0
        end if
        if (State%closed) then
            EV%lmaxt=min(EV%FirstZerolForBeta-1,EV%lmaxt)
            EV%lmaxpolt=min(EV%FirstZerolForBeta-1,EV%lmaxpolt)
            EV%lmaxnrt=min(EV%FirstZerolForBeta-1,EV%lmaxnrt)
            EV%lmaxnut=min(EV%FirstZerolForBeta-1,EV%lmaxnut)
        end if
        EV%MaxlNeededt=max(EV%lmaxpolt,EV%lmaxt, EV%lmaxnrt, EV%lmaxnut)
        if (EV%MaxlNeededt > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupTensorArrayIndices(EV, EV%nvart)
    else
        EV%nvart=0
    end if


    if (CP%WantVectors) then
        EV%lmaxv=max(10,nint(8*CP%Accuracy%lAccuracyBoost))
        EV%lmaxpolv = max(5,nint(5*CP%Accuracy%lAccuracyBoost))

        EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3

        EV%lmaxnrv=nint(30*CP%Accuracy%lAccuracyBoost)

        EV%nvarv=EV%nvarv+EV%lmaxnrv
        if (CP%Num_Nu_massive /= 0 ) then
            call MpiStop('massive neutrinos not supported for vector modes')
        end if
    else
        EV%nvarv=0
    end if

    end subroutine GetNumEqns

    subroutine SwitchToMassiveNuApprox(EV,a, y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i
    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

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

    subroutine MassiveNuVarsOut(EV,y,yprime,a,adotoa,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a, adotoa
    real(dl), optional :: grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu, grhonu, dgrhonu

    grhonu=0
    dgrhonu=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

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
            rhonudot = ThermalNuBack%drho(a*State%nu_masses(nu_i),adotoa)

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
        if (present(dgpi_diff)) dgpi_diff = dgpi_diff + pinu*(3*gpnu_t-grhonu_t)
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
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    associate(nu_q=>State%NuPerturbations%nu_q, nu_int_kernel=>State%NuPerturbations%nu_int_kernel)
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
        do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
            !Get the rest from perturbatively relativistic expansion
            aq=am/nu_q(iq)
            v=1._dl/sqrt(1._dl+aq*aq)
            pert_scale=(State%nu_masses(nu_i)/nu_q(iq))**2/2
            tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind))
            drhonu=drhonu+ tmp/v
            fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
            if (present(dpnu)) then
                dpnu=dpnu+ tmp*v
                pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
            end if
        end do
    end associate
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
    am=a*State%nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        pert_scale=(State%nu_masses(nu_i)/q)**2/2
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
        psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
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

    if (EV%nq(nu_i)/=State%NuPerturbations%nqmax) call MpiStop('Nu_pi: nq/=nqmax')
    pinu=0
    ind=EV%nu_ix(nu_i)+2
    am=a*State%nu_masses(nu_i)
    do iq=1, EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        pinu=pinu+State%NuPerturbations%nu_int_kernel(iq)*y(ind)*v
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
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=State%NuPerturbations%nqmax) call MpiStop('Nu_Intvsq nq/=nqmax0')
    do iq=1, EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        G11=G11+State%NuPerturbations%nu_int_kernel(iq)*y(ind+1)*v**2
        if (EV%lmaxnu_tau(nu_i)>2) then
            G30=G30+State%NuPerturbations%nu_int_kernel(iq)*y(ind+3)*v**2
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
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

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


    function Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, vterm )
    !Delta_Tspin - Delta_TCMB
    use constants
    use DarkAge21cm
    !vterm = hdot/clh + k*n/3/clh
    real(dl), intent(in) :: a,Delta_source,Delta_TCMB,Delta_Tm,Tmat,Trad,xe, Delta_xe, vterm
    real(dl) :: Get21cm_source2
    real(dl) Rgamma,Rm
    real(dl) dC10, n_H,C10, C10_HH, C10_eH
    real(dl) kappa_HH,kappa_eH
    real(dl) tau_eps
    real(dl) H
    real(dl) TSpin
    n_H = State%NNow/a**3
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

    dC10 = (C10*Delta_source + (C10_HH*kappa_HH_21cm(Tmat, .true.) &
        + C10_eH*kappa_eH_21cm(Tmat, .true.)) * &
        Delta_Tm + (kappa_eH-kappa_HH)*xe*n_H*Delta_xe)

    Get21cm_source2 =  dC10*(Rgamma-Rm) +  C10*(Rm*Delta_tm - Delta_TCMB*Rgamma)

    TSpin=CP%Recomb%T_s(a)
    H = (1/(a*dtauda(State,a)))
    tau_eps = a*line21_const*State%NNow/a**3/H/Tspin/1000

    Get21cm_source2 = Get21cm_source2 + &
        tau_eps/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) ) * &
        (Delta_source -vterm + dC10/C10 + 2*( - Rgamma*dC10 + Delta_TCMB*(C10*Rgamma-1)) &
        + Trad/(Tmat-Trad)*(Delta_tm-Delta_TCMB)   )

    end function Get21cm_source2


    function Get21cm_dTs(a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
    !d Delta T_s / d eta dropping small \Delta_xe terms
    use constants
    use DarkAge21cm
    real(dl), intent(in) :: a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe
    real(dl) :: Get21cm_dTs
    real(dl) n_H,C10, C10_HH, C10_eH, delta_C10
    real(dl) kappa_HH,kappa_eH, TSpin

    n_H = State%NNow/a**3
    kappa_HH = kappa_HH_21cm(Tmat, .false.)
    kappa_eH = kappa_eH_21cm(Tmat, .false.)
    C10_HH = n_H*kappa_HH* (1- xe)
    C10_eH = n_H*kappa_eH*xe
    C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
    TSpin=CP%Recomb%T_s(a)
    delta_C10 = C10*Delta_n + (C10_HH*kappa_HH_21cm(Tmat, .true.) &
        +C10_eH*kappa_eH_21cm(Tmat, .true.))*Delta_Tm

    !          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
    !                  n_H*kappa_pH_21cm(Tmat, .false.)*xe

    Get21cm_dTs =  4*a*( TSpin/TMat*(Delta_Tm-Delta_ts)*C10 + (1-TSpin/TMat)*delta_C10 + &
        (Trad*Delta_TCMB - Tspin*Delta_Ts)*A10/T_21cm ) * MPC_in_sec

    end function Get21cm_dTs


    subroutine output_window_sources(EV, sources, y, yprime, &
        tau, a, adotoa, grho, gpres, &
        k, etak, z, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxnu, Delta_TM, Delta_xe,  &
        dgq, qg,  vb, qgdot, vbdot, &
        dgpi, pig, pigdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau)
    !Line of sight sources for number counts, lensing and 21cm redshift windows
    type(EvolutionVars) EV
    real(dl) y(EV%nvar), yprime(EV%nvar)
    real(dL), intent(out) :: sources(:)
    real(dL), intent(in) :: tau, a, adotoa, grho, gpres, &
        k,etak, z, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxnu,  &
        dgq, qg, vb, qgdot, vbdot, &
        dgpi, pig, pigdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3), &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau
    real(dl), intent(in) :: Delta_TM, Delta_xe
    real(dl) s(0:10), t(0:10)
    real(dl) counts_radial_source, counts_velocity_source, counts_density_source, counts_ISW_source, &
        counts_redshift_source, counts_timedelay_source, counts_potential_source
    integer w_ix, lineoff,lineoffpol
    real(dl) Delta_TCMB
    integer j
    real(dl) Tmat,Trad, Delta_source, Delta_source2
    real(dl) xe, chi, polter_line

    j = EV%OutputStep
    if (CP%SourceTerms%line_reionization) sources(2)=0

    if (tau <= EV%ThermoData%tau_start_redshiftwindows) return

    !There are line of sight contributions...
    if (CP%Do21cm) then
        Delta_TCMB = clxg/4
        Delta_source = clxb
        Trad = CP%TCMB/a

        xe = CP%Recomb%x_e(a)
        Tmat = CP%Recomb%T_m(a)

        Delta_source2 = Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, &
            k*(z+vb)/adotoa/3)
    end if

    do w_ix = 1, State%num_redshiftwindows
        associate (W => State%Redshift_W(w_ix))

            if (W%kind == window_lensing) then
                sources(3+w_ix) =-2*phi*W%win_lens(j)
            elseif (W%kind == window_counts) then
                !assume zero velocity bias and relevant tracer is CDM perturbation
                !neglect anisotropic stress in some places

                !Main density source
                if (CP%SourceTerms%counts_density) then
                    counts_density_source= W%wing(j)*(clxc*W%Window%GetBias(k,a) + (W%comoving_density_ev(j) - 3*adotoa)*sigma/k)
                    !Newtonian gauge count density; bias assumed to be on synchronous gauge CDM density
                else
                    counts_density_source= 0
                endif

                if (CP%SourceTerms%counts_redshift) then
                    !Main redshift distortion from kV_N/H j'' integrated by parts twice (V_N = sigma in synch gauge)
                    counts_redshift_source = ((4.D0*adotoa**2+gpres+grho/3.D0)/k*W%wing2(j)+ &
                        (-4.D0*W%dwing2(j)*adotoa+W%ddwing2(j))/k)*sigma+(-etak/adotoa*k/3.D0-dgrho/ &
                        adotoa/6.D0+(etak/adotoa*k/3.D0+dgrho/adotoa/6.D0+(dgq/2.D0-2.D0*etak*adotoa)/k) &
                        /EV%Kf(1))*W%wing2(j)+2.D0*W%dwing2(j)*etak/k/EV%Kf(1)
                else
                    counts_redshift_source= 0
                end if

                ! 2v j'/(H\chi) geometric term
                if (State%tau0-tau > 0.1_dl .and. CP%SourceTerms%counts_radial) then
                    chi =State%tau0-tau
                    counts_radial_source= (1-2.5*W%Window%dlog10Ndm)*((-4.D0*W%wing2(j)/chi*adotoa &
                        -2.D0*(-W%dwing2(j)*chi-W%wing2(j))/chi**2)/ &
                        k*sigma+2.D0*W%wing2(j)*etak/chi/k/EV%Kf(1))
                else
                    counts_radial_source = 0
                end if

                if (CP%SourceTerms%counts_timedelay) then
                    !time delay; WinV is int g/chi
                    counts_timedelay_source= 2*(1-2.5*W%Window%dlog10Ndm)*W%WinV(j)*2*phi
                else
                    counts_timedelay_source = 0
                end if

                if (CP%SourceTerms%counts_ISW) then
                    !WinF is int wingtau
                    counts_ISW_source = W%WinF(j)*2*phidot
                else
                    counts_ISW_source = 0
                end if

                if (CP%SourceTerms%counts_potential) then
                    !approx phi = psi
                    counts_potential_source = ( phidot/adotoa + phi +(5*W%Window%dlog10Ndm-2)*phi ) * W%wing(j) &
                        + phi * W%wingtau(j)
                else
                    counts_potential_source = 0
                end if

                if (CP%SourceTerms%counts_velocity) then
                    counts_velocity_source =  (-2.D0*W%wingtau(j)*adotoa+W%dwingtau(j))/k*sigma &
                        +W%wingtau(j)*etak/k/EV%Kf(1) &
                        - counts_radial_source  !don't double count terms; counts_radial is part of counts_velocity with 1/H/chi
                else
                    counts_velocity_source = 0
                end if

                sources(3+w_ix)=  counts_radial_source +  counts_density_source + counts_redshift_source &
                    + counts_timedelay_source + counts_potential_source &
                    + counts_ISW_source + counts_velocity_source

                sources(3+w_ix)=sources(3+w_ix)/W%Fq

                if (CP%SourceTerms%counts_lensing) &
                    sources(3+W%mag_index+State%num_redshiftwindows) = phi*W%win_lens(j)*(2-5*W%Window%dlog10Ndm)
            elseif (W%kind == window_21cm) then
                if (CP%SourceTerms%line_basic) then
                    sources(3+w_ix)= exptau*(W%wing(j)*Delta_source + W%wing2(j)*Delta_source2 &
                        - W%Wingtau(j)*(clxb - (Delta_source2+clxg/4)))
                    !!    sources(3+w_ix)= exptau*W%wing(j)*phi
                else
                    sources(3+w_ix)= 0
                end if

                if (CP%SourceTerms%line_distortions ) then
                    !With baryon velocity, dropping small terms
                    s(1) =  (sigma/adotoa/3.D0-etak/adotoa**2/3.D0)*W%wing(j)*exptau*k
                    s(2) =  -1.D0/adotoa**2*exptau*W%wing(j)*dgrho/6.D0+((((4.D0*sigma+ &
                        vb)*adotoa+(-grho*sigma/2.D0-vb*grho/3.D0)/adotoa+(sigma*grho**2/18.D0+ &
                        vb*grho**2/18.D0)/adotoa**3)*W%wing(j)-4.D0*W%dwing(j)*sigma+(W%ddwing(j)*sigma+ &
                        W%ddwing(j)*vb)/adotoa+(W%dwing(j)*sigma*grho/3.D0+W%dwing(j)*vb*grho/3.D0)/ &
                        adotoa**2-2.D0*W%dwing(j)*vb+((-2.D0*etak+etak*grho/adotoa**2/3.D0)*W%wing(j) &
                        + 2.D0*W%dwing(j)*etak/adotoa)/EV%Kf(1))*exptau+ &
                        (-4.D0*visibility*sigma- 2.D0*visibility*vb+ &
                        (dvisibility*sigma+dvisibility*vb)/adotoa+(visibility*grho*sigma/3.D0+ &
                        visibility*vb*grho/3.D0)/adotoa**2)*W%wing(j)+2.D0*visibility*etak/adotoa*W%wing(j)/ &
                        EV%Kf(1)+(2.D0*visibility*W%dwing(j)*sigma+2.D0*visibility*W%dwing(j)*vb)/adotoa)/k
                    t(0) =  s(1)+s(2)

                    sources(3+w_ix)= sources(3+w_ix) + t(0)
                end if


                if (CP%SourceTerms%line_extra) then
                    !All sources except below
                    if (CP%SourceTerms%line_basic .and. CP%SourceTerms%line_distortions) then
                        sources(3+w_ix) =  (-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)*exptau*k+ &
                            (W%wing2(j)*Delta_source2+W%wing(j)*Delta_source+1.D0/adotoa*W%winV(j)*dgrho/3.D0)* &
                            exptau+((-W%dwing(j)*vb+(-(3.D0*gpres+grho)*sigma/3.D0 &
                            - 4.D0*adotoa**2*sigma)*W%winV(j)+4.D0*adotoa*W%dwinV(j)*sigma+(-sigma- &
                            vb)*W%ddWinV(j)-vbdot*W%wing(j)-W%dwinV(j)*vbdot+(-2.D0*W%dwinV(j)*etak &
                            + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*exptau-2.D0*visibility*sigma*W%dwinV(j)+ &
                            (4.D0*visibility*sigma*adotoa-dvisibility*sigma)*W%winV(j)-2.D0*visibility*W%winV(j)*etak/ &
                            EV%Kf(1)-visibility*W%dwinV(j)*vb-visibility*W%wing(j)*vb)/k+((2.D0*W%dwinV(j)*dgpi+ &
                            diff_rhopi*W%winV(j))*exptau+2.D0*visibility*W%winV(j)*dgpi)/k**2
                    else
                        s(1) =  ((-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)+(-sigma/adotoa/3.D0+ &
                            etak/adotoa**2/3.D0)*W%wing(j))*exptau*k+(1.D0/adotoa*W%winV(j)*dgrho/3.D0 &
                            + 1.D0/adotoa**2*W%wing(j)*dgrho/6.D0)*exptau
                        s(2) =  s(1)
                        s(6) =  ((-vb-sigma)*W%ddWinV(j)+(-4.D0*adotoa**2*sigma-&
                            (18.D0*gpres+ 6.D0*grho)*sigma/18.D0)*W%winV(j)+((-4.D0*sigma-vb)*adotoa-vbdot+&
                            (grho*sigma/ 2.D0+vb*grho/3.D0)/adotoa+(-grho**2*sigma/18.D0-vb*grho**2/18.D0)/ &
                            adotoa**3)*W%wing(j)+W%dwing(j)*vb+(-W%ddwing(j)*sigma-W%ddwing(j)*vb)/adotoa &
                            + 4.D0*W%dwinV(j)*sigma*adotoa+4.D0*W%dwing(j)*sigma+(-W%dwing(j)*grho*sigma/3.D0- &
                            W%dwing(j)*vb*grho/3.D0)/adotoa**2-W%dwinV(j)*vbdot+((2.D0*etak-etak*grho/ &
                            adotoa**2/3.D0)*W%wing(j)-2.D0*W%dwing(j)*etak/adotoa-2.D0*W%dwinV(j)*etak &
                            + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*exptau-visibility*W%dwinV(j)*vb+ &
                            (4.D0*visibility*sigma*adotoa-dvisibility*sigma)*W%winV(j)
                        s(5) =  s(6)+(-2.D0*visibility*etak/adotoa*W%wing(j)-2.D0*visibility*W%winV(j)*etak)/ &
                            EV%Kf(1)+(4.D0*visibility*sigma+(-visibility*grho*sigma/3.D0-visibility*vb*grho/3.D0)/ &
                            adotoa**2+visibility*vb+(-dvisibility*sigma-dvisibility*vb)/adotoa)*W%wing(j)+ &
                            (-2.D0*visibility*W%dwing(j)*sigma-2.D0*visibility*W%dwing(j)*vb)/adotoa &
                            - 2.D0*visibility*W%dwinV(j)*sigma
                        s(6) =  1.D0/k
                        s(4) =  s(5)*s(6)
                        s(5) =  ((diff_rhopi*W%winV(j)+2.D0*W%dwinV(j)*dgpi)*exptau &
                            + 2.D0*visibility*dgpi*W%winV(j))/k**2
                        s(3) =  s(4)+s(5)
                        t(0) =  s(2)+s(3)

                        sources(3+w_ix) =   sources(3+w_ix) + t(0)
                    end if
                end if


                if (CP%SourceTerms%line_reionization) then
                    if (State%num_redshiftwindows>1) stop 'reionization only for one window at the mo'
                    lineoff=EV%reion_line_ix
                    lineoffpol = lineoff+EV%lmaxline-1

                    if (tau  < State%tau0) then
                        polter_line = 0.1_dl*y(lineoff+2)+9._dl/15._dl*y(lineoffpol+2)
                        sources(2)=visibility*polter_line*(15._dl/2._dl)/(f_K(State%tau0-tau)*k)**2
                    else
                        sources(2)=0
                    end if

                    if (.not. CP%SourceTerms%use_21cm_mK) sources(2)= sources(2) /W%Fq

                    s(1) =  visibility*y(lineoff+2)/4.D0+visibility*y(lineoff)
                    s(2) =  s(1)
                    s(4) =  (-1.D0/EV%Kf(1)*visibility*W%winV(j)*etak/10.D0-visibility*sigma*W%dwinV(j)/10.D0 &
                        - 9.D0/20.D0*visibility*yprime(lineoff+2)-27.D0/100.D0*visibility*opacity*y(lineoff+1) &
                        - 9.D0/10.D0*dvisibility*y(lineoff+3)-3.D0/20.D0*visibility*opacity*EV%Kf(2)*y(lineoffpol+3)+ &
                        visibility*W%dwinV(j)*vb+81.D0/200.D0*visibility*opacity*y(lineoff+3) &
                        +3.D0/5.D0*dvisibility*y(lineoff+1)+3.D0/10.D0*visibility*yprime(lineoff+1)+ &
                        (visibility*adotoa*sigma/5.D0+(36.D0*visibility*opacity-80.D0*dvisibility)*sigma/400.D0+ &
                        dvisibility*vb+visibility*vbdot)*W%winV(j))/k
                    s(5) =  (visibility*W%winV(j)*dgpi/10.D0+9.D0/20.D0*visibility*dopacity*y(lineoffpol+2) &
                        + 261.D0/400.D0*visibility*opacity**2.D0*y(lineoff+2)&
                        -117.D0/200.D0*visibility*opacity**2.D0*y(lineoffpol+2)+3.D0/4.D0*ddvisibility*y(lineoff+2) &
                        - 27.D0/20.D0*dvisibility*opacity*y(lineoff+2)+9.D0/10.D0*dvisibility*opacity*y(lineoffpol+2)&
                        -27.D0/40.D0*visibility*dopacity*y(lineoff+2))/k**2
                    s(3) =  s(4)+s(5)
                    t(0) =  s(2)+s(3)

                    sources(3+w_ix)= sources(3+w_ix) + t(0)
                end if

                if (CP%SourceTerms%line_phot_quadrupole) then
                    s(1) =  (EV%kf(1)*W%wing2(j)*pig/2.D0+(-clxg/4.D0-5.D0/8.D0*pig)*W%wing2(j))*exptau
                    s(3) =  ((-1.D0/EV%kf(1)*W%wing2(j)*etak+(-sigma+9.D0/8.D0*EV%kf(2)*y(9) &
                        -3.D0/4.D0*qg)*W%dwing2(j)+(-opacity*vb+2.D0*adotoa*sigma+9.D0/8.D0*EV%kf(2)*yprime(9) &
                        + 3.D0/8.D0*opacity*EV%kf(2)*E(3)+3.D0/4.D0*opacity*qg)*W%wing2(j))*exptau+ &
                        (-3.D0/4.D0*visibility*qg-visibility*sigma+9.D0/8.D0*visibility*EV%kf(2)*y(9))*W%wing2(j))/k
                    s(4) =  (((27.D0/16.D0*opacity*pig-15.D0/8.D0*pigdot &
                        -9.D0/8.D0*opacity*E(2))*W%dwing2(j)+(27.D0/16.D0*dopacity*pig &
                        +9.D0/8.D0*opacity**2.D0*E(2)-9.D0/8.D0*opacity**2.D0*polter &
                        +27.D0/16.D0*opacity*pigdot+dgpi-9.D0/8.D0*dopacity*E(2))*W%wing2(j)&
                        -15.D0/8.D0*W%ddwing2(j)*pig)*exptau-15.D0/4.D0*visibility*W%dwing2(j)*pig+(- &
                        (-27.D0*visibility*opacity+30.D0*dvisibility)*pig/16.D0-9.D0/8.D0*visibility*opacity*E(2) &
                        - 15.D0/8.D0*visibility*pigdot)*W%wing2(j))/k**2
                    s(2) =  s(3)+s(4)
                    t(0) =  s(1)+s(2)

                    sources(3+w_ix)= sources(3+w_ix)+ t(0)
                end if


                if (CP%SourceTerms%line_phot_dipole) then
                    sources(3+w_ix)=sources(3+w_ix) + (EV%kf(1)*W%wing2(j)*pig/2.D0-W%wing2(j)*clxg/4.D0)*exptau &
                        +(((vbdot- opacity*vb+3.D0/4.D0*opacity*qg)*&
                        W%wing2(j)+(vb-3.D0/4.D0*qg)*W%dwing2(j))*exptau+&
                        (visibility*vb-3.D0/4.D0*visibility*qg)*W%wing2(j))/k
                end if

                if (.not. CP%SourceTerms%use_21cm_mK) sources(3+w_ix)= sources(3+w_ix) /W%Fq
            end if
        end associate
    end do
    end subroutine output_window_sources

    subroutine output(EV, y, j, tau,sources, num_custom_sources)
    type(EvolutionVars) EV
    real(dl) y(EV%nvar), yprime(EV%nvar)
    integer, intent(in) :: j
    real(dl) tau
    real(dl), target :: sources(:)
    integer, intent(in) :: num_custom_sources

    yprime = 0
    EV%OutputSources => Sources
    EV%OutputStep = j
    if (num_custom_sources>0) &
        EV%CustomSources => sources(State%CLdata%CTransScal%NumSources - num_custom_sources+1:)
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputSources, EV%CustomSources)

    end subroutine output

    subroutine outputt(EV,yt,n,tau,dt,dte,dtb)
    !calculate the tensor sources for open and closed case
    implicit none
    integer n
    type(EvolutionVars) :: EV
    real(dl), target :: yt(n), ytprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot,polterddot,prefac
    real(dl) pig, pigdot, octg, aux, polter, shear, adotoa,a
    real(dl) sinhxr,cothxor
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Bprime,Eprime
    real(dl), target :: pol(3),polEprime(3), polBprime(3)
    real(dl) opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow

    call derivst(EV,EV%nvart,tau,yt,ytprime)

    k2=EV%k2_buf
    k=EV%k_buf
    aux=EV%aux_buf
    shear = yt(ixt_shear)

    x=(State%tau0-tau)/State%curvature_radius
    call EV%ThermoData%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow)

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
        adotoa = 1/(a*dtauda(State,a))
        pigdot=32._dl/45._dl*k/opacity*(2._dl*adotoa*shear+ytprime(ixt_shear))
        pig = 32._dl/45._dl*k/opacity*shear
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

    sinhxr=State%rofChi(x)*State%curvature_radius

    if (EV%q*sinhxr > 1.e-8_dl) then
        prefac=sqrt(EV%q2*State%curvature_radius*State%curvature_radius-State%Ksign)
        cothxor=State%cosfunc(x)/sinhxr

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
        polterddot = 9._dl/15._dl*(-dopacity*(E(2)-polter)-opacity*(  &
            Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - 5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
            +0.1_dl*(k*(-octg*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(ixt_shear)) - &
            dopacity*(pig - polter) - opacity*(pigdot-polterdot))

        dt=(shear*exptau + (15._dl/8._dl)*polter*visibility/k)*State%curvature_radius/sinhxr**2/prefac

        dte=State%curvature_radius*15._dl/8._dl/k/prefac* &
            ((ddvisibility*polter + 2._dl*dvisibility*polterdot + visibility*polterddot)  &
            + 4._dl*cothxor*(dvisibility*polter + visibility*polterdot) - &
            visibility*polter*(k2 -6*cothxor**2))

        dtb=15._dl/4._dl*EV%q*State%curvature_radius/k/prefac*(visibility*(2._dl*cothxor*polter + polterdot) + dvisibility*polter)
    else
        dt=0._dl
        dte=0._dl
        dtb=0._dl
    end if

    end subroutine outputt


    subroutine outputv(EV,yv,n,tau,dt,dte,dtb)
    !calculate the vector sources
    implicit none
    integer n
    type(EvolutionVars) :: EV
    real(dl), target :: yv(n), yvprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot
    real(dl) vb,qg, pig, polter, sigma
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Eprime
    real(dl) a, opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow


    call derivsv(EV,EV%nvarv,tau,yv,yvprime)

    k2=EV%k2_buf
    k=EV%k_buf
    sigma = yv(2)
    vb  = yv(3)
    qg  = yv(4)
    pig = yv(5)


    x=(State%tau0-tau)*k

    if (x > 1.e-8_dl) then
        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)

        call EV%ThermoData%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
            visibility, dvisibility, ddvisibility, exptau, lenswindow)

        if (yv(1) < 1e-3) then
            dt = 1
        else
            dt =0
        end if
        dt= (4*(vb+sigma)*visibility + 15._dl/2/k*( visibility*polterdot + dvisibility*polter) &
            + 4*(exptau*yvprime(2)) )/x

        dte= 15._dl/2*2*polter/x**2*visibility + 15._dl/2/k*(dvisibility*polter + visibility*polterdot)/x

        dtb= -15._dl/2*polter/x*visibility
    else
        dt=0
        dte=0
        dtb=0
    end if

    end subroutine outputv

    subroutine initial(EV,y, tau)
    !  Scalar initial conditions.
    implicit none

    type(EvolutionVars) EV
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
        Rc,Rb,Rv,Rg,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, iqg, rhomass,a_massive, ep
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
        i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxde=11,i_vde=12
    integer, parameter :: i_max = i_vde
    real(dl) initv(6,1:i_max), initvec(1:i_max)

    nullify(EV%OutputTransfer) !Should not be needed, but avoids issues in ifort 14
    nullify(EV%OutputSources)
    nullify(EV%CustomSources)

    EV%is_cosmological_constant = State%CP%DarkEnergy%is_cosmological_constant

    if (State%flat) then
        EV%k_buf=EV%q
        EV%k2_buf=EV%q2
        EV%Kf(1:EV%MaxlNeeded)=1._dl
    else
        EV%k2_buf=EV%q2-State%curv
        EV%k_buf=sqrt(EV%k2_buf)

        do l=1,EV%MaxlNeeded
            EV%Kf(l)=1._dl-State%curv*(l*(l+2))/EV%k2_buf
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
            ep=ep0*5/CP%Accuracy%AccuracyBoost*0.65
        else
            ep=ep0
        end if
    else
        ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(EV%ThermoData%tight_tau, EV%ThermoData%OpacityToTime(EV%k_buf/ep))

    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(State%grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+State%grhornomass

    om = (State%grhob+State%grhoc)/sqrt(3*(State%grhog+grhonu))
    omtau=om*tau
    Rv=grhonu/(grhonu+State%grhog)

    Rg = 1-Rv
    Rc=CP%omch2/(CP%omch2+CP%ombh2)
    Rb=1-Rc
    Rp15=4*Rv+15

    if (CP%Scalar_initial_condition > initial_nummodes) &
        call MpiStop('Invalid initial condition for scalar modes')

    a=tau*State%adotrad*(1+omtau/4)
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
        if (Rc==0) call MpiStop('Isocurvature initial conditions assume non-zero dark matter')

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

    y(ix_etak)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable

    !  CDM
    y(ix_clxc)=InitVec(i_clxc)

    !  Baryons
    y(ix_clxb)=InitVec(i_clxb)
    y(ix_vb)=InitVec(i_vb)

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    ! DarkEnergy: This initializes also i_vq, when num_perturb_equations is set
    !             to 2.
    if (CP%DarkEnergy%num_perturb_equations > 0) then
        call CP%DarkEnergy%PerturbationInitial(InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1), &
            a, tau,  k)
        y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
            InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1)
    end if

    if (CP%Evolve_delta_Ts) then
        y(EV%Ts_ix) = y(EV%g_ix)/4
    end if

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
        a_massive =  20000*k/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost*CP%Accuracy%lAccuracyBoost
        if (a_massive >=0.99) then
            EV%MassiveNuApproxTime(nu_i)=State%tau0+1
        else if (a_massive > 17.d0/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost) then
            EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),State%DeltaTime(0._dl,a_massive, 0.01_dl))
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
    implicit none
    real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
    integer l
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yt(EV%nvart)
    real(dl) tens0, ep, tensfac

    if (State%flat) then
        EV%aux_buf=1._dl
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
        EV%Kft(1:EV%MaxlNeededt)=1._dl !initialize for flat case
    else
        EV%k2_buf=EV%q2-3*State%curv
        EV%k_buf=sqrt(EV%k2_buf)
        EV%aux_buf=sqrt(1._dl+3*State%curv/EV%k2_buf)
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    do l=1,EV%MaxlNeededt
        if (.not. State%flat) EV%Kft(l)=1._dl-State%curv*((l+1)**2-3)/k2
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
    EV%TightSwitchoffTime = min(EV%ThermoData%tight_tau,EV%ThermoData%OpacityToTime(EV%k_buf/ep))

    rhomass =  sum(State%grhormass(1:CP%Nu_mass_eigenstates))
    omtau = tau*(State%grhob+State%grhoc)/sqrt(3*(State%grhog+rhomass+State%grhornomass))
    a=tau*State%adotrad*(1+omtau/4)

    if (DoTensorNeutrinos) then
        bigR = (rhomass+State%grhornomass)/(rhomass+State%grhornomass+State%grhog)
    else
        bigR = 0._dl
    end if

    x=k*tau

    tens0 = 1

    yt(ixt_H)= tens0
    !commented things are for the compensated mode with magnetic fields; can be neglected
    !-15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))

    elec=-tens0*(1+2*State%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1

    !shear
    yt(ixt_shear)=-5._dl/2/(bigR+5)*x*elec
    !          + 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))

    yt(ixt_shear+1:EV%nvart)=0._dl

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

    if (State%flat) then
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
    else
        call MpiStop('Vectors not supported in non-flat models')
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    omtau = tau*(State%grhob+State%grhoc)/sqrt(3*(State%grhog+State%grhornomass))

    a=tau*State%adotrad*(1+omtau/4)

    x=k*tau

    bigR = (State%grhornomass)/(State%grhornomass+State%grhog)
    Rc=CP%omch2/(CP%omch2+CP%ombh2)

    yv(1)=a !Could eliminate this, but rarely used anyway

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
    use Transfer
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

    subroutine derivs(EV,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the scalar perturbations
    use constants, only : barssc0, Compton_CT, line21_const
    use MassiveNu
    use Recombination
    implicit none
    type(EvolutionVars) EV
    integer n,nu_i
    real(dl) ay(n),ayprime(n)
    real(dl) tau, w
    real(dl) k,k2
    real(dl) opacity
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t,sigma,polter
    real(dl) w_dark_energy_t !equation of state of dark energy
    real(dl) gpres_noDE !Pressure with matter and radiation, no dark energy
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) E2, dopacity
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz
    real(dl) dgpi,dgrho_matter,grho_matter, clxnu, gpres_nu
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    real(dl) xe,Trad, Delta_TM, Tmat, Delta_TCMB
    real(dl) delta_p_b, wing_t, wing2_t,winv_t
    real(dl) Delta_source2, polter_line
    real(dl) Delta_xe, Tspin, tau_eps, tau_fac, Tb
    integer lineoff,lineoffpol
    !Variables for source calculation
    real(dl) diff_rhopi, pidot_sum, dgpi_diff, phi
    real(dl) E(2:3), Edot(2:3)
    real(dl) phidot, polterdot, polterddot, octg, octgdot
    real(dl) ddopacity, visibility, dvisibility, ddvisibility, exptau, lenswindow
    real(dl) ISW, quadrupole_source, doppler, monopole_source, tau0, ang_dist
    real(dl) dgrho_de, dgq_de, cs2_de

    k=EV%k_buf
    k2=EV%k2_buf

    !  Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call EV%ThermoData%Values(tau,a,cs2,opacity,dopacity)
    else
        call EV%ThermoData%Values(tau,a,cs2,opacity)
    end if
    a2=a*a

    etak=ay(ix_etak)

    !  CDM variables
    clxc=ay(ix_clxc)

    !  Baryon variables
    clxb=ay(ix_clxb)
    vb=ay(ix_vb)
    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2

    if (EV%is_cosmological_constant) then
        grhov_t = State%grhov * a2
        w_dark_energy_t = -1_dl
    else
        call State%CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)
    end if

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    gpres_nu=0
    grhonu_t=0

    if (State%CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,grhonu_t,gpres_nu,dgrho_matter,dgq, wnu_arr)
    end if

    grho_matter=grhonu_t+grhob_t+grhoc_t
    grho = grho_matter+grhor_t+grhog_t+grhov_t
    gpres_noDE = gpres_nu + (grhor_t + grhog_t)/3

    if (State%flat) then
        adotoa=sqrt(grho/3)
        cothxor=1._dl/tau
    else
        adotoa=sqrt((grho+State%grhok)/3._dl)
        cothxor=1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
    end if

    dgrho = dgrho_matter

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

    if (.not. EV%is_cosmological_constant) then
        call State%CP%DarkEnergy%PerturbedStressEnergy(dgrho_de, dgq_de, &
            a, dgq, dgrho, grho, grhov_t, w_dark_energy_t, gpres_noDE, etak, &
            adotoa, k, EV%Kf(1), ay, ayprime, EV%w_ix)
        dgrho = dgrho + dgrho_de
        dgq = dgq + dgq_de
    end if

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    if (State%flat) then
        !eta*k equation
        sigma=(z+1.5_dl*dgq/k2)
        ayprime(ix_etak)=0.5_dl*dgq
    else
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
        ayprime(ix_etak)=0.5_dl*dgq + State%curv*z
    end if

    if (.not. EV%is_cosmological_constant) &
        call State%CP%DarkEnergy%PerturbationEvolve(ayprime, w_dark_energy_t, &
        EV%w_ix, a, adotoa, k, z, ay)

    !  CDM equation of motion
    clxcdot=-k*z
    ayprime(ix_clxc)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(ix_clxb)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    !Sources
    if (EV%Evolve_baryon_cs) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            Tmat = State%CP%Recomb%T_m(a)
        else
            Tmat = State%CP%TCMB/a
        end if
        if (EV%Evolve_TM) then
            Delta_TM = ay(EV%Tg_ix)
        else
            Delta_TM = clxg/4
        end if
        delta_p_b = barssc0*(1._dl-0.75d0*State%CP%yhe+(1._dl-State%CP%yhe)*opacity*a2/State%akthom)*Tmat*(clxb + delta_tm)
    else
        Delta_TM = clxg/4
        delta_p_b = cs2*clxb
    end if

    if (State%CP%Evolve_delta_xe) then
        if (EV%saha) then
            xe=State%CP%Recomb%x_e(a)
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
        vbdot=(-adotoa*vb+cs2*k*clxb + k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip
        EV%pig = pig

    else
        vbdot=-adotoa*vb+k*delta_p_b-photbar*opacity*(4._dl/3*vb-qg)
    end if

    ayprime(ix_vb)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+delta_p_b*k)/pb43 &
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

    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            Delta_TCMB = clxg/4
            xe = State%CP%Recomb%x_e(a)
            Trad = State%CP%TCMB/a

            !Matter temperature
            !Recfast_CT = (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc [a_R = radiation constant]
            ayprime(EV%Tg_ix) = -2*k*(z+vb)/3 - a*  Compton_CT * (Trad**4) * xe / (1._dl+xe+State%fHe) * &
                ((1- Trad/Tmat)*(Delta_TCMB*4 + Delta_xe/(1+xe/(1+State%fHe))) &
                + Trad/Tmat*(Delta_Tm - Delta_TCMB)  )

            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) =  Get21cm_dTs(a,clxb,ay(EV%Ts_ix),Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
            end if
        else
            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) = -k*(4._dl/3._dl*z+qg)/4  !Assume follows Delta_TM which follows clxg
            end if
        end if
    end if

    if (State%CP%Evolve_delta_xe .and. .not. EV%saha) then
        ayprime(EV%xe_ix) = &
            State%CP%Recomb%dDeltaxe_dtau(a, Delta_xe,clxb, Delta_Tm, k*z/3,k*vb, adotoa)
    end if

    if (State%CP%Do21cm) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            if (State%CP%SourceTerms%line_reionization) then
                lineoff = EV%reion_line_ix+1
                lineoffpol = lineoff+EV%lmaxline-1

                if (tau> EV%ThermoData%tau_start_redshiftwindows) then
                    !Multipoles of 21cm

                    polter_line = ay(lineoff+2)/10+9._dl/15*ay(lineoffpol+2)

                    call interp_window(State%TimeSteps,State%Redshift_W(1),tau,wing_t,wing2_t,winv_t)

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
                        - wing2_t * ay(EV%g_ix+2)/4

                    do  l=3,EV%lmaxline-1
                        ayprime(lineoff+l)=EV%denlk(l)*ay(lineoff+l-1)-EV%denlk2(l)*ay(lineoff+l+1)-opacity*ay(lineoff+l) &
                            - wing2_t * ay(EV%g_ix+l)/4
                    end do
                    !truncate
                    ayprime(lineoff+EV%lmaxline)=k*ay(lineoff+EV%lmaxline-1)-(EV%lmaxline+1)*cothxor*ay(lineoff+EV%lmaxline)  &
                        -opacity*ay(lineoff+EV%lmaxline) - wing2_t * ay(EV%g_ix+EV%lmaxline)/4

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



    !  Massive neutrino equations of motion.
    if (State%CP%Num_Nu_massive >0) then
        !DIR$ LOOP COUNT MIN(1), AVG(1)
        do nu_i = 1, State%CP%Nu_mass_eigenstates
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
                !DIR$ LOOP COUNT MIN(3), AVG(3)
                do i=1,EV%nq(nu_i)
                    q=State%NuPerturbations%nu_q(i)
                    aq=a*State%nu_masses(nu_i)/q
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
    end if

    if (associated(EV%OutputTransfer) .or. associated(EV%OutputSources)) then
        if (EV%TightCoupling .or. EV%no_phot_multpoles) then
            E=0
            Edot=0
        else
            E = ay(EV%polind+2:EV%polind+3)
            Edot = ayprime(EV%polind+2:EV%polind+3)
        end if
        if (EV%no_nu_multpoles) then
            pirdot=0
            qrdot = -4*dz/3
        end if
        if (EV%no_phot_multpoles) then
            pigdot=0
            octg=0
            octgdot=0
            qgdot = -4*dz/3
        else
            if (EV%TightCoupling) then
                if (second_order_tightcoupling) then
                    octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity)
                    E(2) = pig/4 + pigdot*(1._dl/opacity)*(-5._dl/8._dl)
                    E(3) = (3._dl/7._dl)*(EV%k_buf/opacity)*E(2)
                    Edot(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity/opacity**2))
                else
                    pigdot = -dopacity/opacity*pig + 32._dl/45*k/opacity*(-2*adotoa*sigma  &
                        +etak/EV%Kf(1)-  dgpi/k +vbdot )
                    Edot(2) = pigdot/4
                    E(2) = pig/4
                    octg=0
                end if
                octgdot=0
            else
                octg=ay(EV%g_ix+3)
                octgdot=ayprime(EV%g_ix+3)
            end if
        end if
        if (EV%is_cosmological_constant) then
            dgrho_de=0
            dgq_de=0
        end if

        dgpi  = grhor_t*pir + grhog_t*pig
        dgpi_diff = 0  !sum (3*p_nu -rho_nu)*pi_nu
        pidot_sum = grhog_t*pigdot + grhor_t*pirdot
        clxnu =0
        if (State%CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, adotoa, dgpi=dgpi, clxnu_all=clxnu, &
                dgpi_diff=dgpi_diff, pidot_sum=pidot_sum)
        end if
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff)*adotoa + &
            State%CP%DarkEnergy%diff_rhopi_Add_Term(dgrho_de, dgq_de, grho, &
            gpres, w_dark_energy_t, State%grhok, adotoa, &
            EV%kf(1), k, grhov_t, z, k2, ayprime, ay, EV%w_ix)
        phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)

        if (associated(EV%OutputTransfer)) then
            EV%OutputTransfer(Transfer_kh) = k/(State%CP%h0/100._dl)
            EV%OutputTransfer(Transfer_cdm) = clxc
            EV%OutputTransfer(Transfer_b) = clxb
            EV%OutputTransfer(Transfer_g) = clxg
            EV%OutputTransfer(Transfer_r) = clxr
            EV%OutputTransfer(Transfer_nu) = clxnu
            EV%OutputTransfer(Transfer_tot) =  dgrho_matter/grho_matter !includes neutrinos
            EV%OutputTransfer(Transfer_nonu) = (grhob_t*clxb+grhoc_t*clxc)/(grhob_t + grhoc_t)
            EV%OutputTransfer(Transfer_tot_de) =  dgrho/grho_matter
            !Transfer_Weyl is k^2Phi, where Phi is the Weyl potential
            EV%OutputTransfer(Transfer_Weyl) = k2*phi
            EV%OutputTransfer(Transfer_Newt_vel_cdm)=  -k*sigma/adotoa
            EV%OutputTransfer(Transfer_Newt_vel_baryon) = -k*(vb + sigma)/adotoa
            EV%OutputTransfer(Transfer_vel_baryon_cdm) = vb
            if (State%CP%do21cm) then
                Tspin = State%CP%Recomb%T_s(a)
                xe = State%CP%Recomb%x_e(a)

                tau_eps = a*line21_const*State%NNow/a**3/adotoa/Tspin/1000
                delta_source2 = Get21cm_source2(a,clxb,clxg/4,Delta_Tm,Delta_xe,Tmat,&
                    State%CP%TCMB/a,xe,k*(z+vb)/adotoa/3)
                tau_fac = tau_eps/(exp(tau_eps)-1)
                EV%OutputTransfer(Transfer_monopole) = ( clxb + Trad/(Tspin-Trad)*delta_source2 )  &
                    + (tau_fac-1)*(clxb - (delta_source2 + clxg/4)  )

                EV%OutputTransfer(Transfer_vnewt) = tau_fac*k*(vb+sigma)/adotoa
                EV%OutputTransfer(Transfer_Tmat) =  delta_TM
                if (State%CP%SourceTerms%use_21cm_mK) then
                    Tb = (1-exp(-tau_eps))*a*(Tspin-Trad)*1000

                    EV%OutputTransfer(Transfer_monopole) = EV%OutputTransfer(Transfer_monopole)*Tb
                    EV%OutputTransfer(Transfer_vnewt) = EV%OutputTransfer(Transfer_vnewt)*Tb
                    EV%OutputTransfer(Transfer_Tmat) = EV%OutputTransfer(Transfer_Tmat)*Tb
                end if
            end if
        end if
        if (associated(EV%OutputSources)) then

            EV%OutputSources = 0
            call EV%ThermoData%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
                visibility, dvisibility, ddvisibility, exptau, lenswindow)

            tau0 = State%tau0
            phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                diff_rhopi+ k*sigma*(gpres + grho))/k2
            !time derivative of shear
            sigmadot = -adotoa*sigma - 1.0d0/2.0d0*dgpi/k + k*phi
            !quadrupole source derivatives; polter = pi_g/10 + 3/5 E_2
            polter = pig/10+9._dl/15*E(2)
            polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &
                + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity
            !Temperature source terms, after integrating by parts in conformal time

            !2phi' term (\phi' + \psi' in Newtonian gauge), phi is the Weyl potential
            ISW = 2*phidot*exptau
            monopole_source =  (-etak/(k*EV%Kf(1)) + 2*phi + clxg/4)*visibility
            doppler = ((sigma + vb)*dvisibility + (sigmadot + vbdot)*visibility)/k
            quadrupole_source = (5.0d0/8.0d0)*(3*polter*ddvisibility + 6*polterdot*dvisibility &
                + (k**2*polter + 3*polterddot)*visibility)/k**2

            EV%OutputSources(1) = ISW + doppler + monopole_source + quadrupole_source
            ang_dist = f_K(tau0-tau)
            if (tau < tau0) then
                !E polarization source
                EV%OutputSources(2)=visibility*polter*(15._dl/8._dl)/(ang_dist**2*k2)
                !factor of four because no 1/16 later
            end if

            if (size(EV%OutputSources) > 2) then
                !Get lensing sources
                if (tau>State%tau_maxvis .and. tau0-tau > 0.1_dl) then
                    EV%OutputSources(3) = -2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
                    !We include the lensing factor of two here
                end if
            end if
            if (State%num_redshiftwindows > 0) then
                call output_window_sources(EV, EV%OutputSources, ay, ayprime, &
                    tau, a, adotoa, grho, gpres, &
                    k, etak, z, ayprime(ix_etak), phi, phidot, sigma, sigmadot, &
                    dgrho, clxg,clxb,clxc,clxnu, Delta_TM, Delta_xe, &
                    dgq, qg, vb, qgdot, vbdot, &
                    dgpi, pig, pigdot, diff_rhopi, &
                    polter, polterdot, polterddot, octg, octgdot, E, Edot, &
                    opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau)
            end if
            if (associated(EV%CustomSources)) then
                select type(DE=>State%CP%DarkEnergy)
                class is (TDarkEnergyEqnOfState)
                    cs2_de = DE%cs2_lam
                class default
                    cs2_de=1
                end select
                block
                    procedure(TSource_func), pointer :: custom_sources_func

                    call c_f_procpointer(CP%CustomSources%c_source_func,custom_sources_func)

                    call custom_sources_func(EV%CustomSources, tau, a, adotoa, grho, gpres,w_dark_energy_t, cs2_de, &
                        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
                        k, etak, ayprime(ix_etak), phi, phidot, sigma, sigmadot, &
                        dgrho, clxg,clxb,clxc,clxr,clxnu, dgrho_de/grhov_t, delta_p_b, &
                        dgq, qg, qr, dgq_de/grhov_t, vb, qgdot, qrdot, vbdot, &
                        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
                        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
                        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
                        tau0, State%tau_maxvis, EV%Kf,f_K)
                end block
            end if
        end if
    end if

    end subroutine derivs



    subroutine derivsv(EV,n,tau,yv,yvprime)
    !  Evaluate the time derivatives of the vector perturbations, flat case
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


    !  Get sound speed and opacity, and see if should use tight-coupling

    call EV%ThermoData%Values(tau,a, cs2,opacity)
    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if
    a2=a*a

    finished_tightcoupling = &
        ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4))


    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2
    call CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)

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


    subroutine derivst(EV,n,tau,ayt,aytprime)
    !  Evaluate the time derivatives of the tensor perturbations.
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l,i,ind, nu_i
    real(dl), target ::  ayt(n),aytprime(n)
    real(dl) tau, rhopi,opacity,pirdt
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl) q,aq,v
    real(dl) Hchi,pinu, pig, polter
    real(dl) k,k2,a,a2,grhog_t, grhor_t
    real(dl) pir, adot, adotoa, rhonu, shear
    real(dl) cothxor

    k2=EV%k2_buf
    k= EV%k_buf

    call EV%ThermoData%Expansion_Values(tau,a, adot,opacity)

    Hchi=ayt(ixt_H)

    shear=ayt(ixt_shear)

    a2=a*a
    adotoa = adot/a

    if (State%flat) then
        cothxor=1._dl/tau
    else
        cothxor=1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
    end if

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
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2

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
        if (State%CP%Num_Nu_massive > 0) then
            do nu_i=1,State%CP%Nu_mass_eigenstates
                if (.not. EV%EvolveTensorMassiveNu(nu_i)) then
                    rhopi=rhopi+ State%grhormass(nu_i)/a2*pir !- good approx, note no rhonu weighting
                else
                    ind=EV%nu_ix(nu_i)+2

                    pinu= Nu_pi(EV, ayt, a, nu_i)
                    rhopi=rhopi+ State%grhormass(nu_i)/a2*pinu

                    do i=1,State%NuPerturbations%nqmax
                        q=State%NuPerturbations%nu_q(i)
                        aq=a*State%nu_masses(nu_i)/q
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

    if (State%flat) then
        aytprime(ixt_shear)=-2*adotoa*shear+k*Hchi-rhopi/k
    else
        aytprime(ixt_shear)=-2*adotoa*shear+k*Hchi*(1+2*State%curv/k2)-rhopi/k
    endif

    aytprime(ixt_H)=-k*shear

    end subroutine derivst

    end module GaugeInterface
