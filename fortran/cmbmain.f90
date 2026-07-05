    !     This this is the main CAMB program module.
    !
    !     Code for Anisotropies in the Microwave Background
    !     by Antony lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.

    !     Note that though the code is internally parallelised, it is not thread-safe
    !     so you cannot generate more than one model at the same time in different threads.
    !
    !     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
    !     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
    !     Original CMBFAST copyright and disclaimer:
    !
    !     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
    !     the Massachusetts Institute of Technology.  All rights reserved.
    !
    !     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
    !     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPlIED.
    !     By way of example, but not limitation,
    !     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
    !     MERCHANTABIlITY OR FITNESS FOR ANY PARTICUlAR PURPOSE OR THAT
    !     THE USE OF THE lICENSED SOFTWARE OR DOCUMENTATION WIll NOT INFRINGE
    !     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !     portions of this software are based on the COSMICS package of
    !     E. Bertschinger.  See the lICENSE file of the COSMICS distribution
    !     for restrictions on the modification and distribution of this software.

    module CAMBmain

    !     This code evolves the linearized perturbation equations of general relativity,
    !     the Boltzmann equations and the fluid equations for perturbations
    !     of a Friedmann-Robertson-Walker universe with a supplied system of gauge-dependent equation
    !     in a modules called GaugeInterface.  The sources for the line of sight integral are
    !     computed at sampled times during the evolution for various of wavenumbers. The sources
    !     are then interpolated to a denser wavenumber sampling for computing the line of
    !     sight integrals of the form Integral d(conformal time) source_k * bessel_k_l.
    !     For flat models the bessel functions are interpolated from a pre-computed table, for
    !     non-flat models the hyperspherical Bessel functions are computed by integrating their
    !     differential equation. Both phases ('Evolution' and 'Integration') can do separate
    !     wavenumbers in parallel.

    !     The time variable is conformal  time dtau=dt/a(t) and the spatial dependence is Fourier transformed
    !     with q=sqrt(k**2 + (|m|+1)K), comoving distances are x=r/a(t), with a(t)=1 today.
    !     The units of both length and time are Mpc.

    !    Many elements are part of derived types (to make thread safe or to allow non-sequential code use
    !    CP = CAMB parameters
    !    EV = Time evolution variables
    !    IV = Source integration variables
    !    CT = Cl transfer data
    !    MT = matter transfer data


    use precision
    use results
    use Interpolation, only : cubic_spline_second_derivs
    use GaugeInterface
    use HypersphericalBesselOlver, only: u_olver, phi_olver
    use SpherBessels
    use MassiveNu
    use InitialPower
    use SourceWindows
    use Recombination, only: TRecfast
    use RangeUtils
    use constants, only: const_fourpi, const_pi
    use DarkEnergyInterface
    use MathUtils
    use RungeKuttaDP45Module, only : RungeKuttaDP45Settings
    implicit none
    private

    logical :: WantLateTime = .false. !if lensing or redshift windows

    logical ExactClosedSum  !do all nu values in sum for Cls for Omega_k>0.1

    !Variables for integrating the sources with the bessel functions for each wavenumber
    type IntegrationVars
        integer q_ix
        real(dl) q, dq    !q value we are doing and delta q
        !Contribution to C_l integral from this k
        real(dl), dimension(:,:), pointer :: Source_q
        !Interpolated sources for this k
        integer SourceSteps !number of steps up to where source is zero
    end type IntegrationVars

    real(dl), dimension(:,:), allocatable :: iCl_scalar, iCl_vector, iCl_tensor
    ! Cls at the l values we actually compute,  iCl_xxx(l_index, Cl_type, initial_power_index)

    real(dl), dimension(:,:,:), allocatable :: iCl_Array
    !Full covariance at each L (alternative more general arrangement to above)

    real(dl), parameter :: qmin0=0.1_dl

    real(dl), parameter :: NEARFLAT_TOL = 1e-3_dl !determines whether near-flat approximations used

    real(dl) :: dtaurec_q

    !     qmax - CP%Max_eta_k/State%tau0, qmin = qmin0/State%tau0 for flat case
    real(dl) qmin, qmax

    real(dl) max_etak_tensor , max_etak_vector, max_etak_scalar
    !     Will only be calculated if k*tau < max_etak_xx

    integer maximum_l !Max value of l to compute
    real(dl) :: maximum_qeta = 3000._dl

    integer :: l_smooth_sample = 3000 !assume transfer functions effectively small for q*DMt0>2*l_smooth_sample

    integer :: max_bessels_l_index  = 1000000
    real(dl) :: max_bessels_etak = 1000000*2

    !Evaluate u_olver directly at each source time step, instead of Numerov stepping in chi,
    !once a single time step would need at least this many Numerov substeps (one u_olver
    !evaluation costs about the same as 13 substeps, and avoids accumulated phase drift).
    integer, parameter :: direct_olver_min_substeps = 10

    Type(TTimeSources) , pointer :: ThisSources => null()
    Type(TTimeSources), allocatable, target :: TempSources

    real(dl), dimension(:,:,:), allocatable :: ddScaledSrc !temporary source second derivative array
    real(dl), dimension(:,:,:), pointer ::  ScaledSrc => null() !linear source with optional non-linear scaling

    integer n_source_points ! number of CL source wavenumbers (for use when calculted remaining non-CL transfers)

    procedure(obj_function), private :: dtauda

    public cmbmain, TimeSourcesToCl, ClTransferToCl, InitVars, GetTauStart, CAMBmain_Free !InitVars for BAO hack

    contains


    subroutine cmbmain
    integer q_ix
    type(EvolutionVars) EV
    Type(TTimer) :: Timer
    real(dl) starttime
    Type(ClTransferData), pointer :: ThisCT

    WantLateTime =  CP%DoLensing .or. State%num_redshiftwindows > 0 .or. CP%CustomSources%num_custom_sources>0

    if (CP%WantCls) then
        if (CP%WantTensors .and. CP%WantScalars) call MpiStop('CMBMAIN cannot generate tensors and scalars')
        !Use CAMB_GetResults instead

        if (CP%WantTensors) then
            maximum_l = CP%Max_l_tensor
            maximum_qeta = CP%Max_eta_k_tensor
        else
            maximum_l = CP%Max_l
            maximum_qeta = CP%Max_eta_k
        end if
    end if

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%Start(starttime)

    call InitVars(State) !Most of single thread time spent here (in InitRECFAST)
    if (global_error_flag/=0) return

    if (DebugMsgs .and. Feedbacklevel > 0) then
        call Timer%WriteTime('Timing for InitVars')
        if (.not. State%flat) call WriteFormat('r = %f, scale = %f',State%curvature_radius, State%scale)
    end if

    if (.not. State%HasScalarTimeSources .and. (.not. State%OnlyTransfer .or. &
        CP%NonLinear==NonLinear_Lens .or. CP%NonLinear==NonLinear_both)) &
        call CP%InitPower%Init(CP)
    if (global_error_flag/=0) return

    !Calculation of the CMB and other sources.
    if (CP%WantCls) then
        if (CP%WantScalars) then
            !Allow keeping time source for scalars, so e.g. different initial power spectra
            ! and non-linear corrections can be applied later
            if (.not. allocated(State%ScalarTimeSources)) allocate(State%ScalarTimeSources)
            ThisSources => State%ScalarTimeSources
        else
            if (.not. allocated(TempSources)) allocate(TempSources)
            ThisSources => TempSources
        end if
        call SetkValuesForSources
    end if

    if (CP%WantTransfer) call InitTransfer

    !***note that !$ is the prefix for conditional multi-processor compilation***
    !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)

    if (CP%WantCls) then
        if (DebugMsgs .and. Feedbacklevel > 0) call WriteFormat('Set %d source k values', &
            ThisSources%Evolve_q%npoints)

        if (CP%WantScalars) then
            ThisCT => State%ClData%CTransScal
        else if (CP%WantVectors) then
            ThisCT => State%ClData%CTransVec
        else
            ThisCT => State%ClData%CTransTens
        end if

        call GetSourceMem

        ThisCT%NumSources = ThisSources%SourceNum
        call ThisCT%ls%Init(State,CP%Min_l, maximum_l)

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(DYNAMIC), PRIVATE(EV, q_ix)
        do q_ix= ThisSources%Evolve_q%npoints,1,-1
            if (global_error_flag==0) call DoSourcek(EV,q_ix)
        end do
        !$OMP END PARALLEL DO

        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for source calculation')

    endif !WantCls

    ! If transfer functions are requested, set remaining k values and output
    if (CP%WantTransfer .and. global_error_flag==0) then
        call TransferOut
        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for transfer k values')
    end if

    if (CP%WantTransfer .and. .not. State%OnlyTransfer .and. global_error_flag==0) &
        call Transfer_Get_sigmas(State, State%MT)

    !     if CMB calculations are requested, calculate the Cl by
    !     integrating the sources over time and over k.

    if (CP%WantCls .and. (.not. CP%WantScalars .or. .not. State%HasScalarTimeSources)) then
        call TimeSourcesToCl(ThisCT)

        if (CP%WantScalars) then
            deallocate(State%ScalarTimeSources)
        else
            deallocate(TempSources)
        end if
        nullify(ThisSources)
    end if
    if (DebugMsgs .and. Feedbacklevel > 0) then
        call Timer%WriteTime('Timing for whole of cmbmain', starttime)
    end if

    end subroutine cmbmain

    subroutine TimeSourcesToCl(ThisCT)
    Type(ClTransferData) :: ThisCT
    integer q_ix
    real(dl) :: spline_bessels_etak
    Type(TTimer) :: Timer

    if (CP%WantScalars) ThisSources => State%ScalarTimeSources

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%Start()

    if (CP%WantScalars .and. WantLateTime &
        .and. (CP%NonLinear==NonLinear_Lens .or. CP%NonLinear==NonLinear_both) .and. global_error_flag==0) then
        call MakeNonlinearSources
        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for NonLinear sources')
    else
        ScaledSrc => ThisSources%LinearSrc
    end if

    if (global_error_flag==0) then
        call InitSourceInterpolation

        ExactClosedSum = State%curv > 5e-9_dl .or. State%curv_flat_DM_ratio < 0.93_dl

        max_bessels_l_index = ThisCT%ls%nl
        max_bessels_etak  = maximum_qeta

        if (CP%WantScalars) call GetLimberTransfers(ThisCT)
        ThisCT%max_index_nonlimber = max_bessels_l_index

        if (State%flat .or. State%scale >= 1._dl - near_flat_scale_tol) then
            spline_bessels_etak = max_bessels_etak
            if (.not. State%flat) then
                spline_bessels_etak = ShiftedNuBesselTableMaxEtak(ThisCT%ls%l(max_bessels_l_index), max_bessels_etak)
            end if
            call InitSpherBessels(ThisCT%ls, CP, max_bessels_l_index, spline_bessels_etak)
        end if
        !This is only slow if not called before with same (or higher) Max_l, Max_eta_k
        !Preferably stick to Max_l being a multiple of 50

        call SetkValuesForInt(ThisCT)

        if (DebugMsgs .and. Feedbacklevel > 0) call WriteFormat('Set %d integration k values',ThisCT%q%npoints)

        !Begin k-loop and integrate Sources*Bessels over time
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,4)
        do q_ix=1,ThisCT%q%npoints
            call SourceToTransfers(ThisCT, q_ix)
        end do !q loop
        !$OMP END PARALLEL DO

        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for Integration')
    end if

    if (allocated(ddScaledSrc)) deallocate(ddScaledSrc)
    if (associated(ScaledSrc) .and. .not. associated(ScaledSrc,ThisSources%LinearSrc)) then
        deallocate(ScaledSrc)
        nullify(ScaledSrc)
    end if

    !Final calculations for CMB output unless want the Cl transfer functions only.
    if (.not. State%OnlyTransfer .and. global_error_flag==0) &
        call ClTransferToCl(State)

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for final CL output')

    end subroutine TimeSourcesToCl

    subroutine ClTransferToCl(ActiveState)
    class(CAMBdata) :: ActiveState

    call SetActiveState(ActiveState)
    if (State%CP%WantScalars .and. State%CP%WantCls .and. global_error_flag==0) then
        allocate(iCl_Scalar(State%CLdata%CTransScal%ls%nl,C_Temp:State%Scalar_C_last), source=0._dl)
        if (State%CP%want_cl_2D_array) then
            allocate(iCl_Array(State%CLdata%CTransScal%ls%nl, &
                State%CLdata%CTransScal%NumSources,State%CLdata%CTransScal%NumSources))
            iCl_Array = 0
        end if

        call CalcLimberScalCls(State%CLdata%CTransScal)
        call CalcScalCls(State%CLdata%CTransScal)
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcScalCls'
    end if

    if (State%CP%WantVectors .and. global_error_flag==0) then
        allocate(iCl_vector(State%CLdata%CTransVec%ls%nl,C_Temp:CT_Cross), source=0._dl)
        call CalcVecCls(State%CLdata%CTransVec,GetInitPowerArrayVec)
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcVecCls'
    end if

    if (CP%WantTensors .and. global_error_flag==0) then
        allocate(iCl_Tensor(State%CLdata%CTransTens%ls%nl,CT_Temp:CT_Cross), source=0._dl)
        call CalcTensCls(State%CLdata%CTransTens,GetInitPowerArrayTens)
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcTensCls'
    end if

    if (global_error_flag==0) then
        call State%CLdata%InitCls(State)
        !     Calculating Cls for every l.
        call InterpolateCls()
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'InterplolateCls'
    end if

    if (CP%WantScalars .and. allocated(iCl_Scalar)) deallocate(iCl_scalar)
    if (CP%WantScalars .and. allocated(iCl_Array)) deallocate(iCl_Array)
    if (CP%WantVectors .and. allocated(iCl_Vector)) deallocate(iCl_vector)
    if (CP%WantTensors .and. allocated(iCl_Tensor)) deallocate(iCl_tensor)

    if (global_error_flag/=0) return

    if (CP%OutputNormalization >=2) call State%CLData%NormalizeClsAtl(CP,CP%OutputNormalization)
    !Normalize to C_l=1 at l=OutputNormalization

    end subroutine ClTransferToCl

    subroutine CalcLimberScalCls(CTrans)
    Type(ClTransferData), target :: CTrans
    integer ell, i, s_ix
    real(dl) Cl, reall,fac
    integer s_ix2,j,n
    integer winmin
    Type(LimberRec), pointer :: LimbRec, LimbRec2

    if (CTrans%NumSources <3) return
    if (CP%SourceTerms%limber_phi_lmin>0) then
        winmin = 0
    else
        winmin = 1
    end if

    do i =winmin, State%num_redshiftwindows
        s_ix = 3+i
        if (CTrans%limber_l_min(s_ix) /=0) then
            do j= i, State%num_redshiftwindows
                s_ix2 = 3+j
                if (CTrans%limber_l_min(s_ix2) /=0) then
                    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(Cl,ell,reall,fac,n, LimbRec, LimbRec2)
                    do ell = max(CTrans%limber_l_min(s_ix), CTrans%limber_l_min(s_ix2)), Ctrans%ls%nl
                        !Don't use associate to avoid ifort omp bug
                        LimbRec => CTrans%Limber_windows(s_ix,ell)
                        LimbRec2 => CTrans%Limber_windows(s_ix2,ell)
                        Cl = 0

                        do n = max(LimbRec%n1,LimbRec2%n1), min(LimbRec%n2,LimbRec2%n2)
                            !Actually integral over chi; source has sqrt( chi dchi)
                            !Same n corresponds to same k since ell fixed here
                            Cl = Cl + LimbRec%Source(n)*LimbRec2%Source(n) * CP%InitPower%ScalarPower(LimbRec%k(n))
                        end do

                        reall = real(CTrans%ls%l(ell),dl)
                        fac = (2 * const_pi ** 2)/const_fourpi/(reall+0.5_dl)**3 !fourpi because multipled by fourpi later
                        if (j >= 1) then
                            if (State%Redshift_w(j)%kind == window_lensing) &
                                fac = fac / 2 * reall * (reall + 1)
                        end if
                        if (i >= 1) then
                            if (State%Redshift_w(i)%kind == window_lensing) &
                                fac = fac / 2 * reall * (reall + 1)
                        end if
                        Cl = Cl*fac

                        if(j==0 .and. i==0) iCl_scalar(ell,C_Phi) = Cl
                        if (State%CP%want_cl_2D_array) then
                            iCl_Array(ell,s_ix,s_ix2) = Cl
                            if (i/=j) iCl_Array(ell,s_ix2,s_ix) = Cl
                        end if
                    end do
                    !$OMP END PARALLEL DO
                end if
            end do
        end if
    end do

    end subroutine CalcLimberScalCls

    subroutine GetLimberTransfers(ThisCT)
    Type(ClTransferData), target :: ThisCT
    integer ell, ell_needed
    integer i, s_ix, s_ix_lens
    type(TRedWin), pointer :: W
    integer n1,n2,n, ell_limb
    real(dl) int,k, chi
    integer klo, khi
    real(dl) a0,b0,ho2o6,a03,b03,ho,reall
    Type(LimberRec), pointer :: LimbRec

    call Init_Limber(ThisCT,State)

    if (.not. CP%SourceTerms%limber_windows .or. State%num_redshiftwindows==0 &
        .and. CP%SourceTerms%limber_phi_lmin==0) return

    ! Only if there are extra source windows. Note only flat models supported for sources.

    if (ThisCT%ls%l(ThisCT%ls%nl) > 5000) then
        max_bessels_l_index = ThisCT%ls%indexOf(5000)
    else
        max_bessels_l_index = ThisCT%ls%nl
    end if

    if (CP%Want_CMB) then
        max_bessels_etak= min(ThisCT%ls%l(ThisCT%ls%nl),3000)*2.5_dl*CP%Accuracy%AccuracyBoost
    else
        max_bessels_etak = 5000
    end if

    do i = 0, State%num_redshiftwindows
        s_ix = 3+i

        if (i==0) then
            ell_limb = CP%SourceTerms%limber_phi_lmin
        else
            W => State%Redshift_w(i)
            ell_limb = Win_limber_ell(W, CP,ThisCT%ls%l(ThisCT%ls%nl))
        end if

        ell_needed = ThisCT%ls%l(ThisCT%ls%nl)
        do ell = 1, ThisCT%ls%nl
            if (ThisCT%ls%l(ell) >= ell_limb) then
                ThisCT%limber_l_min(s_ix) =  ell
                ell_needed = ThisCT%ls%l(ell)
                max_bessels_l_index = max(max_bessels_l_index,ThisCT%limber_l_min(s_ix)-1)
                if (FeedbackLevel > 1 .or. DebugMsgs) &
                    call WriteFormat('Limber switch %d: %d', i,  ell_needed)
                exit
            end if
        end do

        if (i==0) then
            if (CP%Want_CMB_lensing) &
                max_bessels_etak = max(max_bessels_etak, min(CP%Max_eta_k, &
                ell_needed * 25._dl * CP%Accuracy%AccuracyBoost))
        else
            max_bessels_etak = max(max_bessels_etak, WindowKmaxForL(W,CP,ell_needed)*State%tau0)
        end if

        if (ThisCT%limber_l_min(s_ix)/=0) then
            s_ix_lens = 0
            if (i == 0) then
                n1 = State%TimeSteps%IndexOf(State%tau_maxvis)
                n2 = State%TimeSteps%npoints - 1
            else
                n1 = State%TimeSteps%IndexOf(W%tau_start)
                if (W%kind == window_lensing .or. W%kind == window_counts &
                    .and. CP%SourceTerms%counts_lensing) then
                    n2 = State%TimeSteps%npoints - 1
                else
                    n2 = min(State%TimeSteps%npoints - 1, State%TimeSteps%IndexOf(W%tau_end))
                end if
                if (W%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                    s_ix_lens = 3 + W%mag_index + State%num_redshiftwindows
                end if
            end if

            do ell = ThisCT%limber_l_min(s_ix), ThisCT%ls%nl
                LimbRec => ThisCT%Limber_windows(s_ix,ell)
                LimbRec%n1 = n1
                LimbRec%n2 = n2
                reall = ThisCT%ls%l(ell)

                allocate(LimbRec%k(n1:n2))
                allocate(LimbRec%Source(n1:n2))

                int = 0
                do n = n1,n2
                    chi = (State%tau0-State%TimeSteps%points(n))
                    k = (reall + 0.5_dl) / chi
                    LimbRec%k(n) = k
                    if (k<=qmax) then
                        klo = ThisSources%Evolve_q%IndexOf(k)
                        khi=klo+1
                        ho=ThisSources%Evolve_q%points(khi)-ThisSources%Evolve_q%points(klo)
                        a0=(ThisSources%Evolve_q%points(khi)-k)/ho
                        b0=(k-ThisSources%Evolve_q%points(klo))/ho
                        ho2o6 = ho**2/6
                        a03=(a0**3-a0)
                        b03=(b0**3-b0)

                        LimbRec%Source(n)= sqrt(chi*State%TimeSteps%dpoints(n))* (a0*ScaledSrc(klo,s_ix,n)+&
                            b0*ScaledSrc(khi,s_ix,n)+(a03 *ddScaledSrc(klo,s_ix,n)+ b03*ddScaledSrc(khi,s_ix,n)) *ho2o6)
                        if (s_ix_lens>0) then
                            LimbRec%Source(n) = LimbRec%Source(n) + reall * (reall + 1) * &
                                sqrt(chi * State%TimeSteps%dpoints(n)) * (a0 * ScaledSrc(klo, s_ix_lens, n) + &
                                b0 * ScaledSrc(khi, s_ix_lens, n) + (a03 * ddScaledSrc(klo, s_ix_lens, n) + &
                                b03 * ddScaledSrc(khi, s_ix_lens, n)) * ho2o6)
                        end if
                    else
                        LimbRec%Source(n)=0
                    end if
                end do
            end do
        else
            max_bessels_l_index = ThisCT%ls%nl
        end if
    end do

    end subroutine GetLimberTransfers

    subroutine SourceToTransfers(ThisCT, q_ix)
    type(ClTransferData), target :: ThisCT
    integer q_ix
    type(IntegrationVars) :: IV

    allocate(IV%Source_q(State%TimeSteps%npoints,ThisSources%SourceNum))

    call IntegrationVars_init(IV)

    IV%q_ix = q_ix
    IV%q =ThisCT%q%points(q_ix)
    IV%dq= ThisCT%q%dpoints(q_ix)

    call InterpolateSources(IV)

    call DoSourceIntegration(IV, ThisCT)

    deallocate(IV%Source_q)

    end subroutine SourceToTransfers


    subroutine InitTransfer
    integer nu,lastnu, ntodo, nq, q_ix, first_i
    real(dl) dlog_lowk1,dlog_lowk, d_osc,dlog_osc, dlog_highk, boost
    real(dl) amin,q_switch_lowk,q_switch_lowk1,q_switch_osc,q_switch_highk
    real(dl), dimension(:), allocatable :: q_transfer
    Type(MatterTransferData), pointer :: MT
    real(dl) k_per_logint

    MT => State%MT
    if (CP%Transfer%k_per_logint==0 .or. State%needs_good_pk_sampling) then
        !Optimized spacing
        !Large log spacing on superhorizon scales
        !Linear spacing for horizon scales and first few baryon oscillations
        !Log spacing for last few oscillations
        !large log spacing for small scales
        !All have at least k_per_logint steps per log k

        boost = CP%Accuracy%AccuracyBoost * CP%Accuracy%TransferkBoost
        k_per_logint =CP%Transfer%k_per_logint
        if (CP%Transfer%high_precision) boost = boost*1.5

        q_switch_lowk1 = 0.7/State%taurst
        dlog_lowk1=max(2*boost, k_per_logint)

        q_switch_lowk = 8/State%taurst
        dlog_lowk=max(8*boost*2.5, k_per_logint)

        q_switch_osc = min(CP%Transfer%kmax,30/State%taurst)
        d_osc= 200*boost*1.8

        if (CP%Transfer%k_per_logint>0) then
            ! only keep linear steps where smaller than needed for k_per_logint
            q_switch_lowk =  min(q_switch_osc, max(q_switch_lowk,  &
                1/d_osc/(exp(1./k_per_logint)-1)))
        end if

        dlog_osc = max(17*boost, k_per_logint)
        q_switch_highk = min(CP%Transfer%kmax,90/State%taurst)

        !Then up to kmax
        dlog_highk = max(3*boost,k_per_logint)

        amin = 5e-5_dl

        nq=int((log(CP%Transfer%kmax/amin))*max(d_osc, k_per_logint))+1
        allocate(q_transfer(nq))

        nq=int((log(q_switch_lowk1/amin))*dlog_lowk1)+1
        do q_ix=1, nq
            q_transfer(q_ix) = amin*exp((q_ix-1)/dlog_lowk1)
        end do
        MT%num_q_trans = nq

        nq=int(log( q_switch_lowk/q_transfer(MT%num_q_trans))*dlog_lowk) +1
        do q_ix=1, nq
            q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_lowk)
        end do
        MT%num_q_trans = MT%num_q_trans + nq

        if (q_switch_osc> q_switch_lowk) then
            nq=int((q_switch_osc-q_transfer(MT%num_q_trans))*d_osc)+1
            do q_ix=1, nq
                q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)+ q_ix/d_osc
            end do
            MT%num_q_trans = MT%num_q_trans + nq
        end if

        if (CP%Transfer%kmax > q_transfer(MT%num_q_trans)) then
            nq=int(log( q_switch_highk/q_transfer(MT%num_q_trans))*dlog_osc) +1
            do q_ix=1, nq
                q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_osc)
            end do
            MT%num_q_trans = MT%num_q_trans + nq
        end if

        if (CP%Transfer%kmax > q_transfer(MT%num_q_trans)) then
            nq=int(log(CP%Transfer%kmax/q_transfer(MT%num_q_trans))*dlog_highk)+1
            dlog_highk = nq/log(CP%Transfer%kmax*1.05/q_transfer(MT%num_q_trans)) !Don't go more than 5% past kmax
            do q_ix=1, nq
                q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_highk)
            end do
            MT%num_q_trans = MT%num_q_trans + nq
        end if

    else
        !Fixed spacing
        MT%num_q_trans=int((log(CP%Transfer%kmax)-log(qmin))*CP%Transfer%k_per_logint)+1
        allocate(q_transfer(MT%num_q_trans))
        do q_ix=1, MT%num_q_trans
            q_transfer(q_ix) = qmin*exp(real(q_ix,dl)/CP%Transfer%k_per_logint)
        end do
    end if

    if (State%closed) then
        lastnu=0
        ntodo = 0
        do q_ix=1,MT%num_q_trans
            nu =nint(State%curvature_radius*q_transfer(q_ix))
            if (.not. ((nu<3).or.(nu<=lastnu))) then
                ntodo=ntodo+1
                q_transfer(ntodo)= nu/State%curvature_radius
                lastnu=nu
            end if
        end do
        MT%num_q_trans = ntodo
    end if

    if (CP%WantCls) then
        ntodo = MT%num_q_trans
        first_i = ntodo+1
        n_source_points = ThisSources%Evolve_q%npoints

        do q_ix = 1,ntodo
            if (q_transfer(q_ix) > ThisSources%Evolve_q%highest+1d-4) then
                !Feb13 fix for case with closed universe where qmax is not neccessarily right quantized value
                !   if (q_transfer(q_ix) > qmax) then
                first_i=q_ix
                exit
            end if
        end do

        if (first_i > ntodo) then
            MT%num_q_trans = ThisSources%Evolve_q%npoints
        else
            MT%num_q_trans = ThisSources%Evolve_q%npoints + (ntodo - first_i+1)
        end if
        call Transfer_Allocate(MT,State)

        MT%q_trans(1:n_source_points) = ThisSources%Evolve_q%points(1:n_source_points)
        if (MT%num_q_trans > n_source_points) then
            MT%q_trans(n_source_points+1:MT%num_q_trans) = q_transfer(first_i:ntodo)
        end if
    else
        n_source_points = 0
        call Transfer_Allocate(MT,State)
        MT%q_trans = q_transfer(1:MT%num_q_trans)
    end if

    deallocate(q_transfer)

    end subroutine InitTransfer

    function GetTauStart(q)
    real(dl), intent(IN) :: q
    real(dl) taustart, GetTauStart

    !     Begin when wave is far outside horizon.
    !     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
    !     of relativistic neutrinos.
    if (State%flat) then
        taustart=0.001_dl/q
    else
        taustart=0.001_dl/sqrt(q**2-State%curv)
    end if

    !     Make sure to start early in the radiation era.
    taustart=min(taustart,0.1_dl)

    !     Start when massive neutrinos are strongly relativistic.
    if (CP%Num_nu_massive>0 .and. any(State%nu_masses(1:CP%Nu_mass_eigenstates)/=0)) then
        taustart=min(taustart,1.d-3/maxval(State%nu_masses(1:CP%Nu_mass_eigenstates))/State%adotrad)
    end if

    GetTauStart=taustart
    end function GetTauStart

    subroutine DoSourcek(EV,q_ix)
    integer q_ix
    real(dl) taustart
    type(EvolutionVars) EV

    EV%q=ThisSources%Evolve_q%points(q_ix)

    EV%q2=EV%q**2

    EV%q_ix = q_ix
    EV%TransferOnly=.false.

    EV%ThermoData => State%ThermoData

    taustart = GetTauStart(EV%q)

    call GetNumEqns(EV)

    if (CP%WantScalars .and. global_error_flag==0) call CalcScalarSources(EV,taustart)
    if (CP%WantVectors .and. global_error_flag==0) call CalcVectorSources(EV,taustart)
    if (CP%WantTensors .and. global_error_flag==0) call CalcTensorSources(EV,taustart)

    end subroutine DoSourcek

    subroutine GetSourceMem
    integer :: err

    if (CP%WantScalars) then
        if (WantLateTime) then
            ThisSources%SourceNum=3
            State%Scalar_C_last = C_PhiE
            ThisSources%NonCustomSourceNum=ThisSources%SourceNum + State%num_redshiftwindows + &
                State%num_extra_redshiftwindows
            ThisSources%SourceNum = ThisSources%NonCustomSourceNum + CP%CustomSources%num_custom_sources
        else
            ThisSources%SourceNum=2
            ThisSources%NonCustomSourceNum = ThisSources%SourceNum
            State%Scalar_C_last = C_Cross
        end if
    else
        ThisSources%SourceNum=3
        ThisSources%NonCustomSourceNum = ThisSources%SourceNum
    end if

    if (allocated(ThisSources%LinearSrc)) &
        deallocate(ThisSources%LinearSrc)
    allocate(ThisSources%LinearSrc(ThisSources%Evolve_q%npoints,&
        ThisSources%SourceNum,State%TimeSteps%npoints), source=0._dl, stat=err)
    if (err/=0) call GlobalError('Sources requires too much memory to allocate', &
        error_unsupported_params)

    end subroutine GetSourceMem



    !  initial variables, number of steps, etc.
    subroutine InitVars(ActiveState)
    type(CAMBdata) :: ActiveState
    real(dl) taumin, maxq, initAccuracyBoost
    integer itf

    call SetActiveState(ActiveState)

    initAccuracyBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%TimeStepBoost

    ! Maximum and minimum k-values.
    if (State%flat) then
        qmax=maximum_qeta/State%tau0
        qmin=qmin0/State%tau0/initAccuracyBoost
    else
        qmax=maximum_qeta/State%DMt0
        qmin=qmin0/State%DMt0/initAccuracyBoost
    end if
    !     Timesteps during recombination (tentative, the actual
    !     timestep is the minimum between this value and taurst/40,
    !     where taurst is the time when recombination starts - see Thermo_Init

    dtaurec_q=4/qmax/initAccuracyBoost
    !AL:Changed Dec 2003, dtaurec feeds back into the non-flat integration via the step size
    State%dtaurec = dtaurec_q
    !dtau rec may be changed by Thermo_Init

    max_etak_tensor = initAccuracyBoost*maximum_qeta /10
    max_etak_scalar = initAccuracyBoost*max(1700._dl,maximum_qeta) /20
    if (maximum_qeta <3500 .and. CP%Accuracy%AccuracyBoost < 2) max_etak_scalar = max_etak_scalar * 1.5
    !tweak to get large scales right
    max_etak_vector = max_etak_scalar

    if (CP%WantCls) then
        maxq = qmax
    else
        maxq=CP%Transfer%kmax
    end if
    if (CP%WantTransfer .and. CP%Transfer%k_per_logint==0) then
        maxq = max(maxq, CP%Transfer%kmax*1.05)
    elseif (CP%WantTransfer .and. CP%Transfer%k_per_logint/=0) then
        maxq = max(maxq, CP%Transfer%kmax*exp(1._dl/CP%Transfer%k_per_logint))
    end if

    taumin=GetTauStart(maxq)

    !     Initialize baryon temperature and ionization fractions vs. time.
    !     This subroutine also fixes the timesteps where the sources are
    !     saved in order to do the integration. So TimeSteps is set here.
    !These routines in ThermoData (modules.f90)
    call State%ThermoData%Init(State,taumin)
    if (global_error_flag/=0) return

    if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'ThermoData.Init'

    !Do any array initialization for propagation equations
    call GaugeInterface_Init

    if (Feedbacklevel > 0)  &
        write(*,'("tau_recomb/Mpc       = ",f7.2,"  tau_now/Mpc = ",f8.1)') State%tau_maxvis,State%tau0

    do itf=1, State%num_redshiftwindows
        associate (Win => State%Redshift_w(itf))
            if (Feedbacklevel > 0) &
                write(*,'("z = ", f7.2,"  tau = ",f8.1,"  chi = ",f8.1, " tau_min =",f7.1, " tau_max =",f7.1)') &
                Win%Redshift, Win%tau,Win%chi0, Win%tau_start,Win%tau_end
        end associate
    end do

    ! Calculating the times for the outputs of the transfer functions.
    if (CP%WantTransfer .or. .not. allocated(State%Transfer_times)) then
        if (allocated(State%Transfer_times)) deallocate(State%Transfer_times)
        !allocate even if no transfer to prevent debugging access errors
        allocate(State%Transfer_Times(0:State%num_transfer_redshifts+1))
        if (CP%WantTransfer) &
            call State%TimeOfzArr(State%Transfer_times(1:),  State%Transfer_redshifts, State%num_transfer_redshifts)
    endif

    end subroutine InitVars

    subroutine SetkValuesForSources
    implicit none
    real(dl) dlnk0, dkn1, dkn2, q_switch, q_cmb, dksmooth
    real(dl) qmax_log
    real(dl) SourceAccuracyBoost, transferLogDensity
    logical highPrecisionTransferSources
    !     set k values for which the sources for the anisotropy and
    !     polarization will be calculated. For low values of k we
    !     use a logarithmic spacing. closed case dealt with by SetClosedkValues

    SourceAccuracyBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%SourcekAccuracyBoost
    highPrecisionTransferSources = AccuracyTarget > 0 .and. CP%Want_CMB .and. CP%WantScalars .and. CP%WantTransfer &
        .and. CP%Transfer%high_precision
    if (CP%WantScalars .and. CP%Reion%Reionization .and. CP%Accuracy%AccuratePolarization) then
        dlnk0=2._dl/10/SourceAccuracyBoost
        !Need this to get accurate low l polarization
    else
        dlnk0=5._dl/10/SourceAccuracyBoost
        if (State%closed) dlnk0=dlnk0/2
    end if

    if (CP%Accuracy%AccurateReionization) dlnk0 = dlnk0/2
    !Low-l reionization polarization benefits from a modest cap on the low-k logarithmic source grid.
    if (AccuracyTarget > 0 .and. CP%Want_CMB .and. CP%WantScalars .and. CP%Reion%Reionization .and. &
        CP%Accuracy%AccuratePolarization .and. CP%Accuracy%AccurateReionization) dlnk0 = min(dlnk0, 0.075_dl)
    if (highPrecisionTransferSources) then
        !High-precision transfer output reuses this CMB source q grid up to q_cmb.
        !Use 20/log interval for automatic transfer sampling, or the requested denser transfer grid.
        transferLogDensity = max(20._dl, real(CP%Transfer%k_per_logint, dl))
        dlnk0 = min(dlnk0, 1._dl/transferLogDensity)
    end if

    dkn1=0.6_dl/State%taurst/SourceAccuracyBoost
    dkn2=0.9_dl/State%taurst/SourceAccuracyBoost/1.2

    if (CP%WantTensors .or. CP%WantVectors) then
        dkn1=dkn1  *0.8_dl
        dlnk0=dlnk0/2 !*0.3_dl
        dkn2=dkn2*0.85_dl
    end if

    q_switch = 2*6.3/State%taurst
    !If the low-k logarithmic grid reaches the horizon-entry switch, skip the intermediate linear interval.
    qmax_log = min(dkn1/dlnk0, q_switch)
    !Want linear spacing for wavenumbers which come inside horizon
    !Could use sound horizon, but for tensors that is not relevant

    q_cmb = 2*l_smooth_sample/State%DMt0*SourceAccuracyBoost  !assume everything is smooth at l > l_smooth_sample
    if (CP%Want_CMB .and. maximum_l > 5000 .and. CP%Accuracy%AccuratePolarization) q_cmb = q_cmb*1.4
    q_cmb = max(q_switch*2, q_cmb)
    !prevent EE going wild in tail
    dksmooth = q_cmb/2/(SourceAccuracyBoost)**2
    if (CP%Want_CMB) dksmooth = dksmooth/6

    associate(Evolve_q => ThisSources%Evolve_q)
        call Evolve_q%Init()
        call Evolve_q%Add_delta(qmin, qmax_log, dlnk0, IsLog = .true.)
        if (qmax > qmax_log .and. q_switch > qmax_log) &
            call Evolve_q%Add_delta(qmax_log, min(qmax,q_switch), dkn1)
        if (qmax > q_switch) then
            call Evolve_q%Add_delta(q_switch, min(q_cmb,qmax), dkn2)
            if (qmax > q_cmb) then
                dksmooth = log(1 + dksmooth/q_cmb)
                call Evolve_q%Add_delta(q_cmb, qmax, dksmooth, IsLog = .true.)
            end if
        end if

        call Evolve_q%GetArray(.false.)

        if (State%closed) call SetClosedkValuesFromArr(Evolve_q, .false.)
    end associate

    end subroutine SetkValuesForSources


    subroutine SetClosedkValuesFromArr(R, forInt)
    type(TRanges), intent(inout) :: R
    integer i,nu,lastnu,nmax
    !nu = 3,4,5... in State%closed case, so set nearest integers from arr array
    logical, intent(in) :: forInt
    integer ix
    real(dl) dnu
    integer, allocatable :: nu_array(:)

    if (forInt .and. nint(R%points(1)*State%curvature_radius)<=3) then
        !quantization is important
        call R%Getdpoints(half_ends = .false.)
        R%dpoints = max(1,int(R%dpoints*State%curvature_radius+0.02))
        lastnu=2
        ix=1
        dnu =R%dpoints(ix)
        nmax=0
        lastnu=2
        allocate(nu_array(R%npoints*2))
        do
            do while (R%dpoints(ix)==dnu .and. ix <R%npoints)
                ix=ix+1
            end do
            do nu=lastnu+1,nint(R%points(ix)*State%curvature_radius), nint(dnu)
                nmax=nmax+1
                nu_array(nmax)= nu
            end do
            lastnu=nu_array(nmax) + nint(dnu)-1
            if (ix==R%npoints) exit
            dnu = R%dpoints(ix)
        end do
        if (nint(R%points(R%npoints)*State%curvature_radius) > nu_array(nmax)) then
            nmax=nmax+1
            nu_array(nmax) = nint(R%points(R%npoints)*State%curvature_radius)
        end if
        deallocate(R%points)
        allocate(R%points(nmax))
        R%points = nu_array(1:nmax)/State%curvature_radius
        deallocate(nu_array)
    else
        lastnu=3
        nmax=1

        do i=2,R%npoints
            nu=nint(R%points(i)*State%curvature_radius)
            if (nu > lastnu) then
                nmax=nmax+1
                lastnu=nu
                R%points(nmax)=nu/State%curvature_radius
            end if
        end do
        R%points(1)=3/State%curvature_radius
    end if

    R%Lowest = R%points(1)
    R%Highest = R%points(nmax)
    R%npoints=nmax

    end subroutine SetClosedkValuesFromArr

    subroutine CalcScalarSources(EV,taustart)
    use FileUtils
    use MpiUtils
    implicit none
    type(EvolutionVars) EV
    real(dl) tau,tol1,tauend, taustart
    integer j,ind,itf
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) w(EV%nvar,9), y(EV%nvar), sources(ThisSources%SourceNum)

    w=0
    y=0
    call initial(EV,y, taustart)
    if (global_error_flag/=0) return

    tau=taustart
    ind=1

    !     Begin timestep loop.
    itf=1
    tol1=base_tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)

    if (CP%WantTransfer) then
        if  (CP%Transfer%high_precision) tol1=tol1/100
        do while (itf <= State%num_transfer_redshifts .and. State%TimeSteps%points(2) > State%Transfer_Times(itf))
            !Just in case someone wants to get the transfer outputs well before recombination
            call GaugeInterface_EvolveScal(EV, tau, y, State%Transfer_Times(itf), tol1, ind, rk_settings, w)
            if (global_error_flag/=0) return
            call outtransf(EV,y, tau, State%MT%TransferData(:,EV%q_ix,itf))
            itf = itf+1
        end do
    end if

    do j=2,State%TimeSteps%npoints
        tauend=State%TimeSteps%points(j)

        if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > State%taurend) &
            .and. .not. WantLateTime .and. (.not.CP%WantTransfer.or.tau > State%Transfer_Times(State%num_transfer_redshifts))) then
            ThisSources%LinearSrc(EV%q_ix,:,j)=0
        else
            !Integrate over time, calulate end point derivs and calc output
            call GaugeInterface_EvolveScal(EV, tau, y, tauend, tol1, ind, rk_settings, w)
            if (global_error_flag/=0) return

            call output(EV,y,j, tau,sources, CP%CustomSources%num_custom_sources)
            ThisSources%LinearSrc(EV%q_ix,:,j)=sources

            !     Calculation of transfer functions.
101         if (CP%WantTransfer.and.itf <= State%num_transfer_redshifts) then
                if (j < State%TimeSteps%npoints) then
                    if (tauend < State%Transfer_Times(itf) .and. State%TimeSteps%points(j+1)  > State%Transfer_Times(itf)) then
                        call GaugeInterface_EvolveScal(EV, tau, y, State%Transfer_Times(itf), tol1, ind, &
                            rk_settings, w)
                        if (global_error_flag/=0) return
                    endif
                end if
                !     output transfer functions for this k-value.

                if (abs(tau-State%Transfer_Times(itf)) < 1.e-5_dl .or. j==State%TimeSteps%npoints) then
                    call outtransf(EV,y, tau, State%MT%TransferData(:,EV%q_ix,itf))
                    itf=itf+1
                    if (j < State%TimeSteps%npoints) then
                        if (itf <= State%num_transfer_redshifts.and. &
                            State%TimeSteps%points(j+1) > State%Transfer_Times(itf)) goto 101
                    else
                        if (abs(tau-State%Transfer_Times(itf-1)) > 5.e-5_dl .and. print_fortran_warnings) then
                            write(*,*) 'WARNING: mismatch in integrated times (CAMB: CalcScalarSources)'
                        end if
                    end if
                endif
            end if
        end if
    end do !time step loop

    end subroutine CalcScalarSources


    subroutine CalcTensorSources(EV,taustart)
    implicit none
    type(EvolutionVars) EV
    real(dl) tau,tol,tauend, taustart
    integer j,ind
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) wt(EV%nvart,9), yt(EV%nvart)

    call initialt(EV,yt, taustart)

    tau=taustart
    ind=1
    tol=base_tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)

    !     Begin timestep loop.
    do j=2,State%TimeSteps%npoints
        tauend=State%TimeSteps%points(j)
        if (EV%q*tauend > max_etak_tensor) then
            ThisSources%LinearSrc(EV%q_ix,:,j) = 0
        else
            call GaugeInterface_EvolveTens(EV, tau, yt, tauend, tol, ind, rk_settings, wt)

            call outputt(EV,yt,EV%nvart,tau,ThisSources%LinearSrc(EV%q_ix,CT_Temp,j),&
                ThisSources%LinearSrc(EV%q_ix,CT_E,j),ThisSources%LinearSrc(EV%q_ix,CT_B,j))
        end if
    end do

    end subroutine CalcTensorSources


    subroutine CalcVectorSources(EV,taustart)
    implicit none
    type(EvolutionVars) EV
    real(dl) tau,tol,tauend, taustart
    integer j,ind
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) wt(EV%nvarv,9), yv(EV%nvarv)

    call initialv(EV,yv, taustart)

    tau=taustart
    ind=1
    tol=base_tol/100/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)


    !     Begin timestep loop.
    do j=2,State%TimeSteps%npoints
        tauend=State%TimeSteps%points(j)

        if ( EV%q*tauend > max_etak_vector) then
            ThisSources%LinearSrc(EV%q_ix,:,j) = 0
        else
            rk_settings%reuse_first_derivative = .true.
            call RungeKuttaDP45(EV, EV%nvarv, derivsv, tau, yv, tauend, tol, ind, rk_settings, EV%nvarv, wt)

            call outputv(EV,yv,EV%nvarv,tau,ThisSources%LinearSrc(EV%q_ix,CT_Temp,j),&
                ThisSources%LinearSrc(EV%q_ix,CT_E,j),ThisSources%LinearSrc(EV%q_ix,CT_B,j))
        end if
    end do

    end subroutine CalcVectorSources


    subroutine TransferOut
    !Output transfer functions for k larger than used for C_l computation
    implicit none
    integer q_ix
    real(dl) tau
    type(EvolutionVars) EV


    if (DebugMsgs .and. Feedbacklevel > 0) &
        call WriteFormat('Transfer k values: %d',State%MT%num_q_trans-n_source_points)

    !     loop over wavenumbers.
    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(DYNAMIC), PRIVATE(EV, tau, q_ix)
    do q_ix=State%MT%num_q_trans, n_source_points+1, -1
        EV%TransferOnly=.true. !in case we want to do something to speed it up

        EV%q= State%MT%q_trans(q_ix)

        EV%q2=EV%q**2
        EV%q_ix = q_ix
        EV%ThermoData => State%ThermoData

        tau = GetTauStart(EV%q)

        call GetNumEqns(EV)

        call GetTransfer(EV, tau)
    end do
    !$OMP END PARALLEL DO

    end subroutine TransferOut

    subroutine GetTransfer(EV,tau)
    type(EvolutionVars) EV
    real(dl) tau
    integer ind, i
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) w(EV%nvar,9), y(EV%nvar)
    real(dl) atol

    atol=base_tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    if (CP%Transfer%high_precision) atol=atol/100

    ind=1
    call initial(EV,y, tau)
    if (global_error_flag/=0) return

    do i=1,State%num_transfer_redshifts
        call GaugeInterface_EvolveScal(EV, tau, y, State%Transfer_Times(i), atol, ind, rk_settings, w)
        if (global_error_flag/=0) return
        call outtransf(EV,y,tau,State%MT%TransferData(:,EV%q_ix,i))
    end do

    end subroutine GetTransfer


    subroutine MakeNonlinearSources
    !Scale lensing and any other source terms by non-linear scaling at each redshift and wavenumber
    use NonLinear
    integer i,ik,first_step
    real (dl) tau
    real(dl) scaling(State%num_transfer_redshifts), ddScaling(State%num_transfer_redshifts)
    real(dl) ho,a0,b0, ascale
    integer tf_lo, tf_hi

    if (allocated(State%CAMB_Pk)) deallocate(State%CAMB_PK)
    allocate(State%CAMB_PK)
    call Transfer_GetMatterPowerData(State, State%MT, State%CAMB_PK)

    call CP%NonLinearModel%GetNonLinRatios(State, State%CAMB_PK)
    first_step=1
    do while(State%TimeSteps%points(first_step) < State%Transfer_Times(1))
        first_step = first_step + 1
    end do
    allocate(ScaledSrc, source = ThisSources%LinearSrc)

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
    !$OMP & PRIVATE(i, scaling, ddScaling, tf_lo, tf_hi, tau, ho, a0, b0, ascale)
    do ik=1, ThisSources%Evolve_q%npoints
        if (CP%Do21cm) then
            ScaledSrc(ik, 4:, :) = ScaledSrc(ik, 4:, :) * State%CAMB_Pk%nonlin_ratio(ik,1)
        elseif (ThisSources%Evolve_q%points(ik)/(CP%H0/100) >  CP%NonLinearModel%Min_kh_nonlinear) then
            !Interpolate non-linear scaling in conformal time
            !Do not use an associate for scaling. It does not work.
            scaling = State%CAMB_Pk%nonlin_ratio(ik,1:State%num_transfer_redshifts)
            if (all(abs(scaling-1) < 5e-4)) cycle
            call cubic_spline_second_derivs(State%Transfer_Times, scaling, State%num_transfer_redshifts,ddScaling)

            tf_lo=1
            tf_hi=tf_lo+1

            do i= first_step, State%TimeSteps%npoints-1
                tau = State%TimeSteps%points(i)

                do while (tau > State%Transfer_Times(tf_hi))
                    tf_lo = tf_lo + 1
                    tf_hi = tf_hi + 1
                end do

                ho=State%Transfer_Times(tf_hi)-State%Transfer_Times(tf_lo)
                a0=(State%Transfer_Times(tf_hi)-tau)/ho
                b0=1-a0

                ascale = a0*scaling(tf_lo)+ b0*scaling(tf_hi)+&
                    ((a0**3-a0)* ddscaling(tf_lo) &
                    +(b0**3-b0)*ddscaling(tf_hi))*ho**2/6

                ScaledSrc(ik,3:,i) = ScaledSrc(ik,3:,i) * ascale
            end  do
        end if
    end do
    !$OMP END PARALLEL DO

    end subroutine MakeNonlinearSources


    subroutine InitSourceInterpolation
    integer i,j
    !     get the interpolation matrix for the sources to interpolate them
    !     for other k-values

    if (allocated(ddScaledSrc)) deallocate(ddScaledSrc)
    allocate(ddScaledSrc(ThisSources%Evolve_q%npoints,ThisSources%SourceNum,State%TimeSteps%npoints))
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(i,j)
    do  i=1,State%TimeSteps%npoints
        do j=1, ThisSources%SourceNum
            call cubic_spline_second_derivs(ThisSources%Evolve_q%points,ScaledSrc(:,j,i), &
                ThisSources%Evolve_q%npoints, ddScaledSrc(:,j,i))
        end do
    end do
    !$OMP END PARALLEL DO

    end subroutine InitSourceInterpolation


    subroutine SetkValuesForInt(ThisCT)
    Type(ClTransferData) :: ThisCT
    integer no
    real(dl) dk,dk0,dlnk1, dk2, max_k_dk, k_max_log, k_max_0
    integer lognum
    real(dl)  qmax_int,IntSampleBoost,LowQIntBoost, LowQCurvBoost


    qmax_int = min(qmax,max_bessels_etak/State%tau0)

    IntSampleBoost=CP%Accuracy%AccuracyBoost*CP%Accuracy%IntkAccuracyBoost
    if (do_bispectrum) then
        IntSampleBoost = IntSampleBoost * 2
        if (hard_bispectrum) IntSampleBoost = IntSampleBoost * 2
    end if

    !     Fixing the # of k for the integration.

    call ThisCT%q%Init()

    if (State%closed.and.ExactClosedSum) then
        call ThisCT%q%Add(3/State%curvature_radius, nint(qmax_int*State%curvature_radius)/State%curvature_radius, &
            nint(qmax_int*State%curvature_radius)-3)
        call Init_ClTransfer(ThisCT)
        call ThisCT%q%Getdpoints(half_ends = .false.)
    else
        !Accurate-reionization polarization needs denser CMB integration sampling at low q (reionization bump
        !is at low l). Only the log block (lognum) and the first linear block step (dk0) are scaled here.
        LowQIntBoost = IntSampleBoost
        if (AccuracyTarget > 0 .and. CP%Want_CMB .and. CP%WantScalars .and. CP%Accuracy%AccuratePolarization &
            .and. CP%Accuracy%AccurateReionization) then
            LowQIntBoost = LowQIntBoost*1.8_dl
            if (State%closed) then
                LowQCurvBoost = min(1._dl, abs(State%curv_flat_DM_ratio - 1._dl)/0.01_dl)
                LowQIntBoost = LowQIntBoost*(1._dl + 0.4_dl*LowQCurvBoost)
            end if
        end if
        !Split up into logarithmically spaced intervals from qmin up to k=lognum*dk0
        !then no-lognum*dk0 linearly spaced at dk0 up to no*dk0
        !then at dk up to max_k_dk, then dk2 up to qmax_int
        lognum=nint(10*LowQIntBoost)
        dlnk1=1._dl/lognum
        no=nint(600*IntSampleBoost)
        dk0=1.8_dl/State%DMt0/LowQIntBoost
        dk=3._dl/State%DMt0/IntSampleBoost/1.6

        k_max_log = lognum*dk0
        k_max_0  = no*dk0

        if (do_bispectrum) k_max_0 = max(10.d0,k_max_0)

        dk2 = 0.04/IntSampleBoost  !very small scales
        if (State%num_redshiftwindows>0) dk2 = dk  !very small scales

        call ThisCT%q%Add_delta(qmin, k_max_log, dlnk1, IsLog = .true.)
        call ThisCT%q%Add_delta(k_max_log, min(qmax_int,k_max_0), dk0)

        if (qmax_int > k_max_0) then
            max_k_dk = max(3000, 2*maximum_l)/State%tau0

            call ThisCT%q%Add_delta(k_max_0, min(qmax_int, max_k_dk), dk)
            if (qmax_int > max_k_dk) then
                !This allows inclusion of high k modes for computing BB lensed spectrum accurately
                !without taking ages to compute.
                call ThisCT%q%Add_delta(max_k_dk, qmax_int, dk2)
            end if
        end if

        call Init_ClTransfer(ThisCT)

        if (State%closed) then
            call SetClosedkValuesFromArr(ThisCT%q,.true.)
            ! half end at top for consistency with flat if tail not converged
            call ThisCT%q%Getdpoints(half_ends = .true.)
            ThisCT%q%dpoints(1) = 1/State%curvature_radius
            deallocate(ThisCT%Delta_p_l_k) !Re-do this from Init_ClTransfer because number of points changed
            allocate(ThisCT%Delta_p_l_k(ThisCT%NumSources, &
                min(ThisCT%max_index_nonlimber,ThisCT%ls%nl), ThisCT%q%npoints))
            ThisCT%Delta_p_l_k = 0
        end if

    end if !ExactClosedSum


    end subroutine setkValuesForInt

    subroutine InterpolateSources(IV)
    implicit none

    type(IntegrationVars), intent(inout) :: IV

    integer :: i, j
    integer :: klo, khi, lo, hi, mid
    integer :: ntime, nsrc, npoints
    integer :: step

    real(dl) :: q, tau, qtau, xf
    real(dl) :: ho, ho2o6
    real(dl) :: a0, b0, a03, b03
    real(dl) :: max_etak
    real(dl), parameter :: eps_xf = 1.e-8_dl

    logical :: ignore_etak

    type(TRanges), pointer :: Evolve_q

    Evolve_q => ThisSources%Evolve_q

    q       = IV%q
    ntime   = State%TimeSteps%npoints
    npoints = Evolve_q%npoints
    nsrc    = size(IV%Source_q, 2)

    ! Preserve original default behaviour: step starts at 2.
    step = 2
    IV%SourceSteps = step

    ! Zero all rows that this routine is responsible for.
    if (ntime >= 2) then
        IV%Source_q(2:ntime, :) = 0._dl
    else
        return
    end if

    ! Determine active source type.
    ! This preserves the original overwrite order:
    ! tensors first, then vectors, then scalars.
    ignore_etak = .false.
    max_etak = 0._dl

    if (CP%WantScalars) then
        if (DebugEvolution .or. WantLateTime) then
            ignore_etak = .true.
        else
            max_etak = max_etak_scalar
        end if
    else if (CP%WantTensors) then
        max_etak = max_etak_tensor
    else if (CP%WantVectors) then
        max_etak = max_etak_vector
    end if

    ! ------------------------------------------------------------------
    ! Thread-safe binary search for interpolation interval.
    ! Assumes Evolve_q%points is monotonically increasing.
    ! Clamp outside range to the nearest interval, matching typical
    ! interpolation-table behaviour.
    ! ------------------------------------------------------------------
    if (q <= Evolve_q%points(1)) then
        klo = 1
    else if (q >= Evolve_q%points(npoints)) then
        klo = npoints - 1
    else
        lo = 1
        hi = npoints

        do while (hi - lo > 1)
            mid = (lo + hi) / 2

            if (q >= Evolve_q%points(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do

        klo = lo
    end if

    khi = klo + 1

    ho = Evolve_q%points(khi) - Evolve_q%points(klo)

    if (ho == 0._dl) then
        error stop 'InterpolateSources: duplicate Evolve_q points give zero interval'
    end if

    a0 = (Evolve_q%points(khi) - q) / ho
    b0 = (q - Evolve_q%points(klo)) / ho

    ho2o6 = ho * ho / 6._dl

    a03 = a0 * a0 * a0 - a0
    b03 = b0 * b0 * b0 - b0

    ! ------------------------------------------------------------------
    ! Interpolate source as a function of time for current q.
    !
    ! Because IV%Source_q(2:ntime,:) was already zeroed, only valid
    ! source rows need to be filled.
    !
    ! If State%TimeSteps%points is monotonically increasing, once the
    ! condition fails it will keep failing, so we can exit early.
    ! ------------------------------------------------------------------
    do i = 2, ntime
        tau  = State%TimeSteps%points(i)
        qtau = q * tau
        xf   = q * (State%tau0 - tau)

        if ((ignore_etak .or. qtau < max_etak) .and. xf > eps_xf) then
            step = i

            do j = 1, nsrc
                IV%Source_q(i, j) = a0 * ScaledSrc(klo, j, i) + &
                    b0 * ScaledSrc(khi, j, i) + &
                    (a03 * ddScaledSrc(klo, j, i) + &
                    b03 * ddScaledSrc(khi, j, i)) * ho2o6
            end do
        else
            exit
        end if
    end do

    IV%SourceSteps = step

    end subroutine InterpolateSources


    subroutine IntegrationVars_Init(IV)
    type(IntegrationVars), intent(INOUT) :: IV

    IV%Source_q(1,:)=0
    IV%Source_q(State%TimeSteps%npoints,:) = 0
    IV%Source_q(State%TimeSteps%npoints-1,:) = 0

    end subroutine IntegrationVars_Init


    subroutine DoSourceIntegration(IV, ThisCT) !for particular wave number q
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    integer j,ll,llmax,n
    real(dl) nu, chi_full_max, chi_here, sh
    real(dl) :: near_flat_boost
    real(dl) :: sixpibynu
    real(dl) :: near_flat_time_weight(IV%SourceSteps)
    logical :: use_near_flat_gate, use_shifted_nu_near_flat
    logical :: have_near_flat_time_weight

    nu=IV%q*State%curvature_radius
    sixpibynu  = 6._dl*const_pi/nu
    chi_full_max = (State%tau0 - State%TimeSteps%points(2)) / State%curvature_radius
    near_flat_boost = NearFlatIntegrationBoost()
    have_near_flat_time_weight = .false.

    if (State%closed) then
        if (nu<20 .or. State%tau0/State%curvature_radius+sixpibynu > const_pi/2) then
            llmax=nint(nu)-1
        else
            llmax=nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
            llmax=min(llmax,nint(nu)-1)  !nu >= l+1
        end if
    else
        llmax=nint(IV%q*State%DMt0)
        if (llmax<15) then
            llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
        else
            llmax=nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
        end if
    end if

    if (State%flat) then
        call DoFlatIntegration(IV,ThisCT, llmax)
    else
        do j=1,ThisCT%ls%nl
            ll=ThisCT%ls%l(j)
            if (ll>llmax) exit
            if (enable_olver_source_integration .and. CP%WantScalars .and. ThisSources%SourceNum <= 3) then
                call DoCurvedOlverIntegration(IV, ThisCT, j, ll, nu)
                cycle
            end if
            use_shifted_nu_near_flat = .false.
            use_near_flat_gate = CP%WantScalars .and. ThisSources%SourceNum <= 3 .and. &
                enable_do_near_flat_integration .and. &
                UseNearFlatSmallChiApprox(ll, nu, chi_full_max, near_flat_boost)
            if (use_near_flat_gate) then
                use_shifted_nu_near_flat = enable_shifted_nu_scalar_approx .and. &
                    UseShiftedNuNearFlatIntegration(ll, nu, chi_full_max, near_flat_boost)
                use_near_flat_gate = use_shifted_nu_near_flat .or. enable_near_flat_smallchi_integration
            end if

            if (use_near_flat_gate) then
                if (.not. have_near_flat_time_weight) then
                    do n = 1, IV%SourceSteps
                        chi_here = (State%tau0 - State%TimeSteps%points(n)) / State%curvature_radius
                        if (abs(chi_here) < 1.e-8_dl) then
                            near_flat_time_weight(n) = State%TimeSteps%dpoints(n)
                        else
                            sh = State%rofChi(chi_here)
                            near_flat_time_weight(n) = State%TimeSteps%dpoints(n) * chi_here / sh
                        end if
                    end do
                    have_near_flat_time_weight = .true.
                end if

                if (use_shifted_nu_near_flat) then
                    call DoNearFlatIntegration(IV, ThisCT, j, ll, nu, near_flat_time_weight)
                else
                    call DoNearFlatSmallChiIntegration(IV, ThisCT, j, ll, nu, near_flat_time_weight)
                end if
            else
                call IntegrateSourcesBessels(IV,ThisCT,j,ll,nu)
            end if
        end do !j loop
    end if

    end subroutine DoSourceIntegration

    real(dl) function NearFlatIntegrationBoost() result(boost_scale)

    boost_scale = max(CP%Accuracy%NonFlatIntAccuracyBoost * CP%Accuracy%AccuracyBoost, 1._dl)

    end function NearFlatIntegrationBoost

    logical function UseShiftedNuNearFlatIntegration(l, nu, chi_max, boost_scale) result(use_shifted_nu)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, chi_max
    real(dl), intent(in) :: boost_scale
    integer ell2
    real(dl) :: ell, alpha, tmax
    real(dl) :: arg_err, amp_err, err_shift

    ell2 = l * (l + 1)
    ell = sqrt(real(ell2, dl))
    alpha = nu / ell

    tmax = alpha * chi_max

    ! Residual after absorbing the constant Olver-map correction into F0.
    arg_err = ell * tmax**3 / (90._dl * alpha**4)
    amp_err = tmax**2 / (180._dl * alpha**4)
    err_shift = 0.5_dl * arg_err + amp_err

    use_shifted_nu =  (err_shift < NEARFLAT_TOL / boost_scale)

    end function UseShiftedNuNearFlatIntegration


    logical function UseNearFlatSmallChiApprox(l, nu, chi_max, boost_scale) result(use_smallchi)
    integer, intent(in) :: l
    real(dl), intent(in) :: nu, chi_max
    real(dl), intent(in) :: boost_scale

    real(dl) :: alpha_gate, smallchi_metric

    if (State%scale < 1._dl - near_flat_scale_tol) then
        use_smallchi = .false.
        return
    end if

    alpha_gate = nu / real(l, dl)
    smallchi_metric = real(l, dl)**2 * chi_max**7 / nu
    use_smallchi = (alpha_gate > 2.5_dl .or. (l >= 50 .and. alpha_gate > 1._dl)) .and. &
        smallchi_metric < 0.1_dl / boost_scale

    end function UseNearFlatSmallChiApprox

    function UseLimber(l)
    !Calculate lensing potential power using Limber rather than j_l integration
    !even when sources calculated as part of temperature calculation
    !(Limber better on small scales unless step sizes made much smaller)
    !This affects speed, esp. of non-flat case
    logical :: UseLimber
    integer l

    !note increasing non-limber is not neccessarily more accurate unless AccuracyBoost much higher
    !use **0.5 to at least give some sensitivity to Limber effects
    !Could be lower but care with phi-T correlation at lower L
    if (CP%SourceTerms%limber_windows) then
        UseLimber = l >= CP%SourceTerms%limber_phi_lmin
    else
        UseLimber = l > 400 * (CP%Accuracy%AccuracyBoost * CP%Accuracy%LimberBoost)** 0.5
    end if

    end function UseLimber

    !flat source integration
    subroutine DoFlatIntegration(IV, ThisCT, llmax)
    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    integer llmax
    integer j, win_ix
    logical DoInt
    real(dl) xlim,xlmax1
    real(dl) tmin, tmax
    real(dl) a2, base, coeff1, coeff2, coeff3, dx, J_l, x_hi, left_weight(IV%SourceSteps)
    real(dl) xf, sums(ThisSources%SourceNum)
    real(dl) qmax_int
    integer bes_ix,n, bes_index(IV%SourceSteps)
    integer custom_source_off, s_ix
    integer nwin
    real(dl) :: BessIntBoost, BessIntBoostL, source_window_tail_lmax
    logical :: has_narrow_source_window

    BessIntBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%BessIntBoost
    custom_source_off = State%num_redshiftwindows + State%num_extra_redshiftwindows + 4
    has_narrow_source_window = .false.
    source_window_tail_lmax = 0._dl
    if (State%num_redshiftwindows > 0 .and. CP%WantScalars) then
        source_window_tail_lmax = 20._dl
        do win_ix = 1, State%num_redshiftwindows
            associate(Win => State%Redshift_w(win_ix))
                if (Win%kind == window_21cm) then
                    has_narrow_source_window = .true.
                    source_window_tail_lmax = max(source_window_tail_lmax, 160._dl)
                else if (Win%kind /= window_lensing .and. Win%chi0 > 0._dl .and. Win%sigma_tau > 0._dl) then
                    source_window_tail_lmax = max(source_window_tail_lmax, &
                        min(160._dl, max(20._dl, 4._dl*Win%chi0/Win%sigma_tau)))
                    if (Win%sigma_tau < 0.02_dl*Win%chi0) has_narrow_source_window = .true.
                end if
            end associate
        end do
    end if

    !     Find the position in the xx table for the x correponding to each
    !     timestep

    associate(BesselPoints => BessRanges%points, TimePoints => State%TimeSteps%points)
        !Precompute xf into left_weight as temporary storage, then batch-index and fix weights
        do j=1,IV%SourceSteps
            left_weight(j) = abs(IV%q*(State%tau0-TimePoints(j)))
        end do
        call BessRanges%IndexOfOrdered(left_weight, IV%SourceSteps, bes_index)
        do j=1,IV%SourceSteps
            bes_ix = bes_index(j)
            x_hi = BesselPoints(bes_ix+1)
            dx = x_hi - BesselPoints(bes_ix)
            left_weight(j) = (x_hi - left_weight(j)) / dx
        end do
    end associate

    do j=1,max_bessels_l_index
        if (ThisCT%ls%l(j) > llmax) return
        BessIntBoostL = BessIntBoost
        ! Very narrow redshift windows need a minimum low-L Bessel tail, xlmax1 ~ L*BessIntBoost.
        if (has_narrow_source_window) &
            BessIntBoostL = max(BessIntBoostL, 20._dl/ThisCT%ls%l(j))
        ! Cut where the hyperspherical Bessel before peak, approximated by
        ! j_l(x_eff) with x_eff=nu_eff*chi, is <~1e-4 of its peak.
        xlim = max(0._dl, ThisCT%ls%l(j) - bjl_pre_peak_start_factor*ThisCT%ls%l(j)**(1._dl/3._dl))
        if (full_bessel_integration .or. do_bispectrum) then
            tmin = State%TimeSteps%points(2)
        else
            xlmax1=80*ThisCT%ls%l(j)*BessIntBoostL
            if (source_window_tail_lmax > 0._dl) &
                xlmax1 = xlmax1 * max(1._dl, min(8._dl, source_window_tail_lmax/ThisCT%ls%l(j)))
            tmin=State%tau0-xlmax1/IV%q
            tmin=max(State%TimeSteps%points(2),tmin)
        end if
        tmax=State%tau0-xlim/IV%q
        tmin=max(State%TimeSteps%points(2),tmin)
        if (.not. CP%Want_CMB .and. .not. CP%Want_CMB_lensing) &
            tmin = max(tmin, State%ThermoData%tau_start_redshiftwindows)
        if (tmax < State%TimeSteps%points(2)) exit

        sums = 0
        !As long as we sample the source well enough, it is sufficient to
        !interpolate the Bessel functions only

        if (ThisSources%SourceNum==2) then
            !This is the innermost loop, so we separate the no lensing scalar case to optimize it
            do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                a2=left_weight(n)
                bes_ix=bes_index(n)
                base = bessel_horner(1,bes_ix,j)
                coeff1 = bessel_horner(2,bes_ix,j)
                coeff2 = bessel_horner(3,bes_ix,j)
                coeff3 = bessel_horner(4,bes_ix,j)

                J_l = base + a2*(coeff1 + a2*(coeff2 + a2*coeff3)) !cubic spline in Horner form

                J_l = J_l*State%TimeSteps%dpoints(n)
                sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                sums(2) = sums(2) + IV%Source_q(n,2)*J_l
            end do
        else
            qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoostL/State%tau0*1.2
            DoInt = .not. CP%WantScalars .or. IV%q < qmax_int
            !Do integral if any useful contribution to the CMB, or large scale effects

            if (DoInt) then
                if (CP%CustomSources%num_custom_sources==0 .and. State%num_redshiftwindows==0) then
                    do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                        !Full Bessel integration
                        a2=left_weight(n)
                        bes_ix=bes_index(n)
                        base = bessel_horner(1,bes_ix,j)
                        coeff1 = bessel_horner(2,bes_ix,j)
                        coeff2 = bessel_horner(3,bes_ix,j)
                        coeff3 = bessel_horner(4,bes_ix,j)

                        J_l = base + a2*(coeff1 + a2*(coeff2 + a2*coeff3)) !cubic spline in Horner form
                        J_l = J_l*State%TimeSteps%dpoints(n)

                        !The unwrapped form is faster
                        sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                        sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                        sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                    end do
                else
                    if (State%num_redshiftwindows>0) then
                        nwin = State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows)
                    else
                        nwin = State%TimeSteps%npoints+1
                    end if
                    if (CP%CustomSources%num_custom_sources==0) then
                        do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                            !Full Bessel integration
                            a2=left_weight(n)
                            bes_ix=bes_index(n)
                            base = bessel_horner(1,bes_ix,j)
                            coeff1 = bessel_horner(2,bes_ix,j)
                            coeff2 = bessel_horner(3,bes_ix,j)
                            coeff3 = bessel_horner(4,bes_ix,j)

                            J_l = base + a2*(coeff1 + a2*(coeff2 + a2*coeff3)) !cubic spline in Horner form
                            J_l = J_l*State%TimeSteps%dpoints(n)

                            !The unwrapped form is faster
                            sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                            sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                            sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                            if (n >= nwin) then
                                do s_ix = 4, ThisSources%SourceNum
                                    sums(s_ix) = sums(s_ix) + IV%Source_q(n,s_ix)*J_l
                                end do
                            end if
                        end do
                    else
                        do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                            !Full Bessel integration
                            a2=left_weight(n)
                            bes_ix=bes_index(n)
                            base = bessel_horner(1,bes_ix,j)
                            coeff1 = bessel_horner(2,bes_ix,j)
                            coeff2 = bessel_horner(3,bes_ix,j)
                            coeff3 = bessel_horner(4,bes_ix,j)

                            J_l = base + a2*(coeff1 + a2*(coeff2 + a2*coeff3)) !cubic spline in Horner form
                            J_l = J_l*State%TimeSteps%dpoints(n)

                            !The unwrapped form is faster
                            sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                            sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                            sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                            sums(custom_source_off) = sums(custom_source_off) +  IV%Source_q(n,custom_source_off)*J_l
                            if (n >= nwin) then
                                do s_ix = 4, ThisSources%NonCustomSourceNum
                                    sums(s_ix) = sums(s_ix) + IV%Source_q(n,s_ix)*J_l
                                end do
                            end if
                            do s_ix = custom_source_off+1, custom_source_off+CP%CustomSources%num_custom_sources -1
                                sums(s_ix) = sums(s_ix)  + IV%Source_q(n,s_ix)*J_l
                            end do
                        end do
                    end if
                end if
            end if
            if (.not. DoInt .or. UseLimber(ThisCT%ls%l(j)) .and. CP%WantScalars) then
                !Limber approximation for small scale lensing (better than poor version of above integral)
                xf = State%tau0-(ThisCT%ls%l(j)+0.5_dl)/IV%q
                if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
                    n=State%TimeSteps%IndexOf(xf)
                    xf= (xf-State%TimeSteps%points(n))/(State%TimeSteps%points(n+1)-State%TimeSteps%points(n))
                    sums(3) = (IV%Source_q(n,3)*(1-xf) + xf*IV%Source_q(n+1,3))*&
                        sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/IV%q
                else
                    sums(3)=0
                end if
            end if
            if (.not. DoInt .and. ThisSources%NonCustomSourceNum>3) then
                if (any(ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum)==0 .or. &
                    ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum) > j)) then
                    !When CMB does not need integral but other sources do
                    do n= State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows), &
                        min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))
                        !Full Bessel integration
                        a2 = left_weight(n)
                        bes_ix = bes_index(n)
                        base = bessel_horner(1,bes_ix,j)
                        coeff1 = bessel_horner(2,bes_ix,j)
                        coeff2 = bessel_horner(3,bes_ix,j)
                        coeff3 = bessel_horner(4,bes_ix,j)

                        J_l = base + a2*(coeff1 + a2*(coeff2 + a2*coeff3)) !cubic spline in Horner form
                        J_l = J_l * State%TimeSteps%dpoints(n)

                        sums(4) = sums(4) + IV%Source_q(n, 4) * J_l
                        do s_ix = 5, ThisSources%NonCustomSourceNum
                            sums(s_ix) = sums(s_ix) + IV%Source_q(n, s_ix) * J_l
                        end do
                    end do
                end if
            end if
        end if

        ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums
    end do

    end subroutine DoFlatIntegration


    subroutine DoNearFlatIntegration(IV, ThisCT, j, l, nu, near_flat_time_weight)
    use SpherBessels, only: BessRanges, bessel_horner
    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    integer, intent(in) :: j, l
    real(dl), intent(in) :: nu
    real(dl), intent(in) :: near_flat_time_weight(IV%SourceSteps)

    integer :: n, bes_ix, bes_index(IV%SourceSteps)
    logical :: DoInt
    real(dl) :: xlim, xlmax1, tmin, tmax
    real(dl) :: a2, base, coeff1, coeff2, coeff3, dx, J_l, weight, x_hi, xf
    real(dl) :: nu_eff, nu_eff_tau, nu_amp, chi_here, sh, chi_over_sh
    real(dl) :: qmax_int, sums(ThisSources%SourceNum), left_weight(IV%SourceSteps)
    real(dl) :: BessIntBoost
    real(dl) :: h, h2, F0

    BessIntBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%BessIntBoost

    ! nu_eff2 = nu**2 - State%Ksign * real(l * (l + 1), dl) / 3._dl
    ! if (nu_eff2 <= 0._dl) return
    ! nu_eff = sqrt(nu_eff2)
    ! nu_amp = sqrt(nu_eff / nu)

    h  = State%Ksign / (nu**2 / (l * (l + 1)))
    h2 = h*h
    F0 = 1._dl - h/6._dl - h2*(13._dl/360._dl + 737._dl/45360._dl*h)
    nu_eff = nu * F0
    nu_amp = sqrt(F0)

    nu_eff_tau = nu_eff / State%curvature_radius
    ! Cut where the near-flat hyperspherical Bessel before peak, approximated by
    ! j_l(x_eff) with x_eff=nu_eff*chi, is <~1e-4 of its peak.
    xlim = real(l, dl) - bjl_pre_peak_start_factor * real(l, dl)**(1._dl/3._dl)
    xlim = max(0._dl, xlim)
    tmax = State%tau0 - xlim / nu_eff_tau
    if (tmax < State%TimeSteps%points(2)) return

    if (full_bessel_integration .or. do_bispectrum) then
        tmin = State%TimeSteps%points(2)
    else
        xlmax1 = 80._dl * l * BessIntBoost
        tmin = State%tau0 - xlmax1 / nu_eff_tau
        tmin = max(State%TimeSteps%points(2), tmin)
    end if
    tmin = max(State%TimeSteps%points(2), tmin)

    sums = 0._dl
    DoInt = .true.
    if (ThisSources%SourceNum /= 2) then
        qmax_int = max(850, l) * 3 * BessIntBoost / State%tau0 * 1.2_dl
        DoInt = nu_eff_tau < qmax_int
    end if

    if (DoInt) then
        ! Avoid the per-l spline lookup setup entirely when this mode only contributes via Limber.
        associate(BesselPoints => BessRanges%points, TimePoints => State%TimeSteps%points)
            do n=1,IV%SourceSteps
                left_weight(n) = abs(nu_eff_tau * (State%tau0 - TimePoints(n)))
            end do
            call BessRanges%IndexOfOrdered(left_weight, IV%SourceSteps, bes_index)
            do n=1,IV%SourceSteps
                bes_ix = bes_index(n)
                x_hi = BesselPoints(bes_ix + 1)
                dx = x_hi - BesselPoints(bes_ix)
                left_weight(n) = (x_hi - left_weight(n)) / dx
            end do
        end associate
    end if

    if (ThisSources%SourceNum == 2) then
        do n = State%TimeSteps%IndexOf(tmin), min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))
            a2 = left_weight(n)
            bes_ix = bes_index(n)
            base = bessel_horner(1, bes_ix, j)
            coeff1 = bessel_horner(2, bes_ix, j)
            coeff2 = bessel_horner(3, bes_ix, j)
            coeff3 = bessel_horner(4, bes_ix, j)
            J_l = base + a2 * (coeff1 + a2 * (coeff2 + a2 * coeff3))

            weight = near_flat_time_weight(n) * nu_amp * J_l
            sums(1) = sums(1) + IV%Source_q(n,1) * weight
            sums(2) = sums(2) + IV%Source_q(n,2) * weight
        end do
    else
        if (DoInt) then
            do n = State%TimeSteps%IndexOf(tmin), min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))
                a2 = left_weight(n)
                bes_ix = bes_index(n)
                base = bessel_horner(1, bes_ix, j)
                coeff1 = bessel_horner(2, bes_ix, j)
                coeff2 = bessel_horner(3, bes_ix, j)
                coeff3 = bessel_horner(4, bes_ix, j)
                J_l = base + a2 * (coeff1 + a2 * (coeff2 + a2 * coeff3))

                weight = near_flat_time_weight(n) * nu_amp * J_l
                sums(1) = sums(1) + IV%Source_q(n,1) * weight
                sums(2) = sums(2) + IV%Source_q(n,2) * weight
                sums(3) = sums(3) + IV%Source_q(n,3) * weight
            end do
        end if
        if (.not. DoInt .or. UseLimber(l)) then
            xf = State%tau0 - State%invsinfunc((l + 0.5_dl) / nu) * State%curvature_radius
            if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
                n = State%TimeSteps%IndexOf(xf)
                xf = (xf - State%TimeSteps%points(n)) / (State%TimeSteps%points(n + 1) - State%TimeSteps%points(n))
                sums(3) = (IV%Source_q(n,3) * (1 - xf) + xf * IV%Source_q(n + 1,3)) * &
                    sqrt(const_pi / 2 / (l + 0.5_dl) / sqrt(1._dl - State%Ksign * real(l**2, dl) / nu**2)) / IV%q
            else
                sums(3) = 0._dl
            end if
        end if
    end if

    ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums

    end subroutine DoNearFlatIntegration


    subroutine DoNearFlatSmallChiIntegration(IV, ThisCT, j, l, nu, near_flat_time_weight)
    use SpherBessels, only: BessRanges, bessel_horner
    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    integer, intent(in) :: j, l
    real(dl), intent(in) :: nu
    real(dl), intent(in) :: near_flat_time_weight(IV%SourceSteps)

    integer :: n, n_first, n_last, bes_ix, bes_index(IV%SourceSteps)
    logical :: DoInt
    real(dl) :: xlim, xlmax1, qmax_int, BessIntBoost
    real(dl) :: a2, base, coeff1, coeff2, coeff3, dx, J_l, x_hi, xf
    real(dl) :: alpha2, chi, chi2, chi4, D, D0, D2, D4, ell2, F, F0, F2, F4, h, h2, h3
    real(dl) :: x_vals(IV%SourceSteps), weights(IV%SourceSteps), left_weight(IV%SourceSteps)
    real(dl) :: sums(ThisSources%SourceNum), weight

    BessIntBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%BessIntBoost

    ! phi small where chi <~  chi_t- \frac{5\eta^{1/3}}{\nu\left(1-K\eta^2/\nu^2\right)^{1/6}}
    ! but here x is the j_l argument, so looks like flat formula.
    xlim = max(0._dl, real(l, dl) - bjl_pre_peak_start_factor * real(l, dl)**(1._dl / 3._dl))
    if (full_bessel_integration .or. do_bispectrum) then
        xlmax1 = huge(1._dl)
    else
        xlmax1 = 80._dl * l * BessIntBoost
    end if

    ! The small-chi map has z = chi*F and Jacobian amplitude D^(-1/2),
    ! where F and D are even quartic polynomials in chi for fixed l and nu.
    ! Precompute the curvature-dependent coefficients once per source integral.
    ell2 = real(l * (l + 1), dl)
    alpha2 = nu**2 / ell2
    h = State%Ksign / alpha2
    h2 = h * h
    h3 = h2 * h

    F0 = 1._dl - h / 6._dl - 13._dl * h2 / 360._dl - 737._dl * h3 / 45360._dl
    F2 = -4._dl * h2 * alpha2 / 360._dl - 148._dl * h3 * alpha2 / 45360._dl
    F4 = -48._dl * h3 * alpha2 * alpha2 / 45360._dl

    D0 = F0
    D2 = -12._dl * h2 * alpha2 / 360._dl - 444._dl * h3 * alpha2 / 45360._dl
    D4 = -240._dl * h3 * alpha2 * alpha2 / 45360._dl

    n_first = 0
    n_last = 0
    associate(TimePoints => State%TimeSteps%points)
        do n = 1, IV%SourceSteps
            chi = (State%tau0 - TimePoints(n)) / State%curvature_radius
            if (abs(chi) < 1.e-8_dl) then
                x_vals(n) = nu * chi
                weights(n) = near_flat_time_weight(n)
            else
                chi2 = chi * chi
                chi4 = chi2 * chi2
                F = F0 + F2 * chi2 + F4 * chi4
                D = D0 + D2 * chi2 + D4 * chi4

                x_vals(n) = nu * chi * F
                weights(n) = near_flat_time_weight(n) * F / sqrt(max(D, 1.e-30_dl))
            end if

            if (x_vals(n) <= xlmax1) then
                if (n_first == 0) n_first = n
                if (x_vals(n) >= xlim) then
                    n_last = n
                else if (n_last > 0) then
                    exit
                end if
            end if
        end do
    end associate

    if (n_first == 0 .or. n_last == 0) return

    sums = 0._dl
    DoInt = .true.
    if (ThisSources%SourceNum /= 2) then
        qmax_int = max(850, l) * 3 * BessIntBoost / State%tau0 * 1.2_dl
        DoInt = nu / State%curvature_radius < qmax_int
    end if

    if (DoInt) then
        call BessRanges%IndexOfOrdered(x_vals(n_first:n_last), n_last - n_first + 1, bes_index(n_first:n_last))
        associate(BesselPoints => BessRanges%points)
            do n = n_first, n_last
                bes_ix = min(bes_index(n), BessRanges%npoints - 1)
                x_hi = BesselPoints(bes_ix + 1)
                dx = x_hi - BesselPoints(bes_ix)
                left_weight(n) = (x_hi - x_vals(n)) / dx
            end do
        end associate
    end if

    if (ThisSources%SourceNum == 2) then
        do n = n_first, n_last
            a2 = left_weight(n)
            bes_ix = min(bes_index(n), BessRanges%npoints - 1)
            base = bessel_horner(1, bes_ix, j)
            coeff1 = bessel_horner(2, bes_ix, j)
            coeff2 = bessel_horner(3, bes_ix, j)
            coeff3 = bessel_horner(4, bes_ix, j)
            J_l = base + a2 * (coeff1 + a2 * (coeff2 + a2 * coeff3))
            weight = weights(n) * J_l

            sums(1) = sums(1) + IV%Source_q(n,1) * weight
            sums(2) = sums(2) + IV%Source_q(n,2) * weight
        end do
    else
        if (DoInt) then
            do n = n_first, n_last
                a2 = left_weight(n)
                bes_ix = min(bes_index(n), BessRanges%npoints - 1)
                base = bessel_horner(1, bes_ix, j)
                coeff1 = bessel_horner(2, bes_ix, j)
                coeff2 = bessel_horner(3, bes_ix, j)
                coeff3 = bessel_horner(4, bes_ix, j)
                J_l = base + a2 * (coeff1 + a2 * (coeff2 + a2 * coeff3))
                weight = weights(n) * J_l

                sums(1) = sums(1) + IV%Source_q(n,1) * weight
                sums(2) = sums(2) + IV%Source_q(n,2) * weight
                sums(3) = sums(3) + IV%Source_q(n,3) * weight
            end do
        end if
        if (.not. DoInt .or. UseLimber(l)) then
            xf = State%tau0 - State%invsinfunc((l + 0.5_dl) / nu) * State%curvature_radius
            if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
                n = State%TimeSteps%IndexOf(xf)
                xf = (xf - State%TimeSteps%points(n)) / (State%TimeSteps%points(n + 1) - State%TimeSteps%points(n))
                sums(3) = (IV%Source_q(n,3) * (1 - xf) + xf * IV%Source_q(n + 1,3)) * &
                    sqrt(const_pi / 2 / (l + 0.5_dl) / sqrt(1._dl - State%Ksign * real(l**2, dl) / nu**2)) / IV%q
            else
                sums(3) = 0._dl
            end if
        end if
    end if

    ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums

    end subroutine DoNearFlatSmallChiIntegration


    subroutine DoCurvedOlverIntegration(IV, ThisCT, j, l, nu)
    ! Not used, just for accuracy checking
    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    integer, intent(in) :: j, l
    real(dl), intent(in) :: nu

    integer :: K, n, n_first, n_last
    logical :: DoInt
    real(dl) :: BessIntBoost, chi, chi_start, sh, phi, xlmax1, tmin, tmax, xf
    real(dl) :: qmax_int, sums(ThisSources%SourceNum), phi_vals(IV%SourceSteps)

    BessIntBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%BessIntBoost

    chi_start = HypersphericalPhiOnsetChi(l, nu, merge(1, -1, State%closed))
    tmax = State%tau0 - chi_start * State%curvature_radius
    if (tmax < State%TimeSteps%points(2)) return

    if (full_bessel_integration .or. do_bispectrum) then
        tmin = State%TimeSteps%points(2)
    else
        xlmax1 = 80._dl * l * BessIntBoost
        tmin = State%tau0 - xlmax1 / IV%q
        tmin = max(State%TimeSteps%points(2), tmin)
    end if

    sums = 0._dl
    qmax_int = max(850, l) * 3 * BessIntBoost / State%DMt0 * 1.2_dl
    DoInt = ThisSources%SourceNum /= 3 .or. IV%q < qmax_int

    if (DoInt) then
        n_first = State%TimeSteps%IndexOf(tmin)
        n_last = min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))

        if (n_first <= n_last) then
            K = merge(1, -1, State%closed)
            do n = n_first, n_last
                chi = (State%tau0 - State%TimeSteps%points(n)) / State%curvature_radius
                sh = State%rofChi(chi)
                if (abs(sh) > 1.e-30_dl) then
                    phi = u_olver(l, K, nu, chi) / sh
                else
                    phi = 0._dl
                end if
                phi_vals(n) = phi * State%TimeSteps%dpoints(n)
            end do

            if (ThisSources%SourceNum == 2) then
                do n = n_first, n_last
                    sums(1) = sums(1) + IV%Source_q(n,1) * phi_vals(n)
                    sums(2) = sums(2) + IV%Source_q(n,2) * phi_vals(n)
                end do
            else
                do n = n_first, n_last
                    sums(1) = sums(1) + IV%Source_q(n,1) * phi_vals(n)
                    sums(2) = sums(2) + IV%Source_q(n,2) * phi_vals(n)
                    sums(3) = sums(3) + IV%Source_q(n,3) * phi_vals(n)
                end do
            end if
        end if
    end if

    if (ThisSources%SourceNum == 3 .and. (.not. DoInt .or. UseLimber(l))) then
        xf = State%tau0 - State%invsinfunc((l + 0.5_dl) / nu) * State%curvature_radius
        if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
            n = State%TimeSteps%IndexOf(xf)
            xf = (xf - State%TimeSteps%points(n)) / (State%TimeSteps%points(n + 1) - State%TimeSteps%points(n))
            sums(3) = (IV%Source_q(n,3) * (1 - xf) + xf * IV%Source_q(n + 1,3)) * &
                sqrt(const_pi / 2 / (l + 0.5_dl) / sqrt(1._dl - State%Ksign * real(l**2, dl) / nu**2)) / IV%q
        else
            sums(3) = 0._dl
        end if
    end if

    ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums

    end subroutine DoCurvedOlverIntegration


    real(dl) function HypersphericalPhiOnsetChi(l, nu, K) result(chi_start)
    implicit none
    integer, intent(in) :: l, K
    real(dl), intent(in) :: nu

    real(dl) :: c_start, denom, lambda, ratio, turn_chi

    lambda = real(l, dl) + 0.5_dl
    if (nu <= 0._dl) then
        chi_start = 0._dl
        return
    end if

    ratio = nu / lambda
    c_start = 4.2_dl
    if (K == -1) then
        if (ratio >= 0.2_dl) then
            c_start = 4.21_dl
        else if (ratio >= 0.05_dl) then
            c_start = 4.30_dl
        else
            c_start = min(5.15_dl, 4.21_dl + 0.60_dl * log10(0.05_dl / max(ratio, tiny(1._dl))))
        end if
    end if

    if (1._dl - c_start * lambda**(-2._dl / 3._dl) <= 0._dl) then
        chi_start = 0._dl
        return
    end if

    turn_chi = State%invsinfunc(lambda / nu)
    denom = (1._dl - real(K, dl) * (lambda / nu)**2)**(1._dl / 6._dl)
    chi_start = max(0._dl, turn_chi - c_start * lambda**(1._dl / 3._dl) / (nu * denom))

    end function HypersphericalPhiOnsetChi


    !non-flat source integration

    subroutine IntegrateSourcesBessels(IV,ThisCT,j,l,nu)
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT
    logical DoInt
    integer l,j, nstart,nDissipative,ntop,nbot,nrange,nnow,nmin_osc,nbot_eff,K
    real(dl) nu,ChiDissipative,ChiStart,tDissipative,y1,y1dis
    real(dl) xf,x,chi, miny1,xlmax1,tmin
    real(dl) sums(ThisSources%SourceNum),out_arr(ThisSources%SourceNum), qmax_int
    real(dl) BessIntBoost

    BessIntBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%BessIntBoost

    !Calculate chi where for smaller chi it is dissipative
    x=sqrt(real(l*(l+1),dl))/nu

    ChiDissipative=State%invsinfunc(x)

    ChiStart=ChiDissipative
    !Move down a bit to get smaller value (better accuracy integrating up from small values)
    if (nu<300) ChiStart = max(ChiDissipative-1._dl/nu,1d-6)   !max(ChiDissipative-1._dl/nu,1d-6)

    !Then get nearest source point with lower Chi...
    tDissipative=State%tau0 - State%curvature_radius*ChiStart
    if (tDissipative<State%TimeSteps%points(1)) then
        nDissipative=2
    else
        nDissipative = State%TimeSteps%IndexOf(tDissipative)+1
    endif
    nDissipative=min(nDissipative,State%TimeSteps%npoints-1)

    tDissipative = State%TimeSteps%points(nDissipative)

    ChiStart =  max(1d-8,(State%tau0-tDissipative)/State%curvature_radius)

    !Get values at ChiStart
    K = merge(1, -1, State%closed)
    y1dis = u_olver(l, K, nu, ChiStart)

    nstart=nDissipative
    chi=ChiStart

    if (CP%WantScalars) then !Do Scalars
        if (ThisSources%SourceNum > 3) call MpiStop('Non-flat not implemented for extra sources')
        !Integrate chi down in dissipative region
        ! cuts off when phi gets small
        ! Fix for phi discontinuity in near-flat models - adaptive formula provides superior performance
        ! Adaptive formula for Omega_k discontinuity fix, see https://github.com/cmbant/CAMB/pull/185
        miny1= 0.5d-4/(l+max(0,30-l))/BessIntBoost
        sums=0
        qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/State%DMt0*1.2
        DoInt =  ThisSources%SourceNum/=3 .or. IV%q < qmax_int
        if (DoInt) then
            if ((nstart < min(State%TimeSteps%npoints-1,IV%SourceSteps)).and.(y1dis > miny1)) then
                y1=y1dis
                nnow=nstart
                do nrange = 1,State%TimeSteps%Count
                    if (nrange == State%TimeSteps%count) then
                        ntop = State%TimeSteps%npoints -1
                    else
                        ntop = State%TimeSteps%R(nrange+1)%start_index
                    end if
                    if (nnow < ntop) then
                        call DoRangeInt(IV,chi,ChiDissipative,nnow,ntop,State%TimeSteps%R(nrange)%delta, &
                            nu,j,l,y1,out_arr)
                        sums  = sums + out_arr
                        nnow = ntop
                        if (chi==0) exit !small enough to cut off
                    end if
                end do
            end if !integrate down chi

            !Integrate chi up in oscillatory region, but optionally cut off very high-x tail
            if (nstart > 2) then
                y1=y1dis
                chi=ChiStart
                nnow=nstart

                if (full_bessel_integration .or. do_bispectrum) then
                    nmin_osc = 2
                else
                    xlmax1 = 80._dl * l * BessIntBoost
                    if (State%num_redshiftwindows > 0 .and. CP%WantScalars) then
                        xlmax1 = xlmax1 * 8._dl
                    end if

                    tmin = State%tau0 - xlmax1 / IV%q
                    tmin = max(State%TimeSteps%points(2), tmin)

                    if (tmin >= State%TimeSteps%points(nstart)) then
                        nmin_osc = nstart
                    else
                        nmin_osc = max(2, State%TimeSteps%IndexOf(tmin))
                    end if
                end if

                if (nstart > nmin_osc) then
                    do nrange = State%TimeSteps%Count,1,-1
                        if (nnow <= nmin_osc) exit

                        nbot = State%TimeSteps%R(nrange)%start_index
                        nbot_eff = max(nbot, nmin_osc)
                        if (nnow > nbot_eff) then
                            call DoRangeInt(IV,chi,ChiDissipative,nnow,nbot_eff,State%TimeSteps%R(nrange)%delta, &
                                nu,j,l,y1,out_arr)
                            sums=sums+out_arr
                            if (chi==0) exit !small for remaining region
                            nnow = nbot_eff
                        end if

                        if (nnow == nmin_osc) exit
                    end do
                end if
            end if
        end if !DoInt
        if (ThisSources%SourceNum==3 .and. (.not. DoInt .or. UseLimber(l))) then
            !Limber approximation for small scale lensing (better than poor version of above integral)
            xf = State%tau0-State%invsinfunc((l+0.5_dl)/nu)*State%curvature_radius
            if (xf < State%TimeSteps%Highest .and. xf > State%TimeSteps%Lowest) then
                nbot=State%TimeSteps%IndexOf(xf)
                xf= (xf-State%TimeSteps%points(nbot))/(State%TimeSteps%points(nbot+1)-State%TimeSteps%points(nbot))
                sums(3) = (IV%Source_q(nbot,3)*(1-xf) + xf*IV%Source_q(nbot+1,3))*&
                    sqrt(const_pi/2/(l+0.5_dl)/sqrt(1-State%Ksign*real(l**2, dl)/nu**2))/IV%q
            else
                sums(3) = 0
            end if
        end if

        ThisCT%Delta_p_l_k(:,j,IV%q_ix)=ThisCT%Delta_p_l_k(:,j,IV%q_ix)+sums

    end if !Do Scalars

    if ((CP%WantTensors)) then !Do Tensors
        chi=ChiStart

        !Integrate chi down in dissipative region
        !DoRangeInt cuts off when phi gets small
        miny1= 1.d-6/l/BessIntBoost
        if ((nstart < State%TimeSteps%npoints-1).and.(y1dis>miny1)) then
            y1=y1dis
            nnow=nstart
            do nrange = 1,State%TimeSteps%Count
                if (nrange == State%TimeSteps%count) then
                    ntop = State%TimeSteps%npoints -1
                else
                    ntop = State%TimeSteps%R(nrange+1)%start_index
                end if
                if (nnow < ntop) then
                    call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,ntop,State%TimeSteps%R(nrange)%delta, &
                        nu,l,y1,out_arr)

                    ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + out_arr

                    nnow = ntop
                    if (chi==0) exit
                end if
            end do
        end if


        !Integrate chi up in oscillatory region
        if (nstart > 2) then
            y1=y1dis
            chi=ChiStart

            nnow=nstart
            do nrange = State%TimeSteps%Count,1,-1
                nbot = State%TimeSteps%R(nrange)%start_index
                if (nnow >  nbot) then
                    call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,nbot,State%TimeSteps%R(nrange)%delta, &
                        nu,l,y1,out_arr)
                    ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + out_arr

                    nnow = nbot
                    if (chi==0) exit !small for remaining region
                end if
            end do
        end if

    end if !Do Tensors

    end subroutine IntegrateSourcesBessels


    real(dl) function ShiftedNuBesselTableMaxEtak(lmax, base_etak) result(bessel_etak)
    integer, intent(in) :: lmax
    real(dl), intent(in) :: base_etak
    real(dl), parameter :: safety_margin = 10._dl
    real(dl) :: ell, chi_max, scale_x, curve_x

    ell = real(lmax, dl)

    ! Take the larger of the model-independent conservative bound and the
    ! actual comoving chi range for this specific run.
    chi_max = max(sqrt(max(0._dl, 6._dl * near_flat_scale_tol / (1._dl - near_flat_scale_tol))), &
        State%tau0 / State%curvature_radius)

    scale_x = base_etak / (1._dl - near_flat_scale_tol)
    curve_x = chi_max * sqrt(ell * (ell + 1._dl) / 3._dl)

    bessel_etak = hypot(scale_x, curve_x) + safety_margin
    end function ShiftedNuBesselTableMaxEtak

    subroutine FillShiftedNuPhiVals(j, l, nu, chi_start, dchisource, sgn, nIntSteps, phi_vals, y1, chi)
    use SpherBessels, only: BessRanges, bessel_horner

    integer, intent(in) :: j, l, nIntSteps
    real(dl), intent(in) :: nu, chi_start, dchisource, sgn
    real(dl), intent(out) :: phi_vals(nIntSteps + 1)
    real(dl), intent(inout) :: y1, chi

    integer :: step_ix, ix, nVals, max_ix
    integer :: bes_ix(nIntSteps + 1)

    real(dl) :: nu_eff, nu_amp
    real(dl) :: chi_here, chi_step
    real(dl) :: chi_over_sh, sh
    real(dl) :: x, x_step, x_hi, inv_dx, w
    real(dl) :: jl
    real(dl) :: h, h2, F0, ell
    real(dl) :: x_vals(nIntSteps + 1)

    nVals = nIntSteps + 1
    chi_step = dchisource * sgn

    ell = real(l, dl) * real(l + 1, dl)

    h  = real(State%Ksign, dl) * ell / (nu * nu)
    h2 = h * h

    F0 = 1._dl - h / 6._dl - h2 * (13._dl / 360._dl + 737._dl * h / 45360._dl)

    nu_eff = nu * F0
    nu_amp = sqrt(F0)

    x = nu_eff * chi_start
    x_step = nu_eff * chi_step

    do step_ix = 1, nVals
        x_vals(step_ix) = x
        x = x + x_step
    end do

    call BessRanges%IndexOfOrdered(x_vals, nVals, bes_ix)

    max_ix = BessRanges%npoints - 1
    chi_here = chi_start

    do step_ix = 1, nVals
        ix = min(bes_ix(step_ix), max_ix)

        x_hi = BessRanges%points(ix + 1)
        inv_dx = 1._dl / (x_hi - BessRanges%points(ix))

        w = (x_hi - x_vals(step_ix)) * inv_dx

        jl = bessel_horner(1, ix, j) + w * ( &
            bessel_horner(2, ix, j) + w * ( &
            bessel_horner(3, ix, j) + w *   &
            bessel_horner(4, ix, j)))

        if (abs(chi_here) < 1.e-8_dl) then
            chi_over_sh = 1._dl
        else
            sh = State%rofChi(chi_here)
            chi_over_sh = chi_here / sh
        end if

        phi_vals(step_ix) = nu_amp * chi_over_sh * jl

        chi_here = chi_here + chi_step
    end do

    chi = chi_start + real(nIntSteps, dl) * chi_step
    y1 = 0._dl

    end subroutine FillShiftedNuPhiVals


    subroutine FillSmallChiPhiVals(j, l, nu, chi_start, dchisource, sgn, nIntSteps, phi_vals, y1, chi)
    use SpherBessels, only: BessRanges, bessel_horner

    integer, intent(in) :: j, l, nIntSteps
    real(dl), intent(in) :: nu, chi_start, dchisource, sgn
    real(dl), intent(out) :: phi_vals(nIntSteps + 1)
    real(dl), intent(inout) :: y1, chi

    integer :: step_ix, ix, npts_m1
    integer :: bes_ix(nIntSteps + 1)

    real(dl) :: alpha2, base, chi2, chi_here, chi_over_sh
    real(dl) :: coeff1, coeff2, coeff3, D, D0, D2, D4
    real(dl) :: dx, ell2, F, F0, F2, F4
    real(dl) :: h, h2, h3, jl, sh, w, x_hi
    real(dl) :: dchi
    real(dl) :: x_vals(nIntSteps + 1)

    real(dl), parameter :: small_chi = 1.e-8_dl
    real(dl), parameter :: min_D = 1.e-30_dl

    ell2 = real(l, dl) * real(l + 1, dl)
    alpha2 = (nu * nu) / ell2

    h  = State%Ksign / alpha2
    h2 = h * h
    h3 = h2 * h

    F0 = 1._dl - h / 6._dl - 13._dl * h2 / 360._dl - 737._dl * h3 / 45360._dl
    F2 = -4._dl  * h2 * alpha2 / 360._dl - 148._dl * h3 * alpha2 / 45360._dl
    F4 = -48._dl * h3 * alpha2 * alpha2 / 45360._dl

    D0 = F0
    D2 = -12._dl  * h2 * alpha2 / 360._dl - 444._dl  * h3 * alpha2 / 45360._dl
    D4 = -240._dl * h3 * alpha2 * alpha2 / 45360._dl

    dchi = dchisource * sgn
    npts_m1 = BessRanges%npoints - 1

    ! Build spline x positions first so IndexOfOrdered can process the ordered batch.
    do step_ix = 1, nIntSteps + 1
        chi_here = chi_start + real(step_ix - 1, dl) * dchi

        if (abs(chi_here) < small_chi) then
            x_vals(step_ix) = nu * chi_here
        else
            chi2 = chi_here * chi_here
            F = F0 + chi2 * (F2 + F4 * chi2)
            x_vals(step_ix) = nu * chi_here * F
        end if
    end do

    call BessRanges%IndexOfOrdered(x_vals, nIntSteps + 1, bes_ix)

    do step_ix = 1, nIntSteps + 1
        chi_here = chi_start + real(step_ix - 1, dl) * dchi

        if (abs(chi_here) < small_chi) then
            F = 1._dl
            D = 1._dl
            chi_over_sh = 1._dl
        else
            chi2 = chi_here * chi_here

            F = F0 + chi2 * (F2 + F4 * chi2)
            D = D0 + chi2 * (D2 + D4 * chi2)

            sh = State%rofChi(chi_here)
            chi_over_sh = chi_here / sh
        end if

        ix = min(bes_ix(step_ix), npts_m1)

        x_hi = BessRanges%points(ix + 1)
        dx = x_hi - BessRanges%points(ix)
        w = (x_hi - x_vals(step_ix)) / dx

        base   = bessel_horner(1, ix, j)
        coeff1 = bessel_horner(2, ix, j)
        coeff2 = bessel_horner(3, ix, j)
        coeff3 = bessel_horner(4, ix, j)

        jl = base + w * (coeff1 + w * (coeff2 + w * coeff3))

        phi_vals(step_ix) = chi_over_sh * F * jl / sqrt(max(D, min_D))
    end do

    chi = chi_start + real(nIntSteps, dl) * dchi
    y1 = 0._dl

    end subroutine FillSmallChiPhiVals



    subroutine DoRangeInt(IV, chi, chiDisp, nstart, nend, dtau, nu, j, l, y1, out)
    ! Non-flat version.
    !
    ! Returns chi at the end of the integral, where the integral stops, not
    ! necessarily at nend.
    !
    ! Integrates source*phi from nstart to nend. phi is obtained from Numerov
    ! stepping of the hyperspherical Bessel ODE, bootstrapped by Olver
    ! evaluations, or by direct Olver evaluation when too many Numerov substeps
    ! would be needed. dtau is the spacing of the equally spaced time steps.

    use precision
    implicit none

    type(IntegrationVars), intent(in) :: IV

    integer, intent(in) :: nstart, nend, j, l
    real(dl), intent(in) :: chiDisp, dtau, nu
    real(dl), intent(inout) :: chi, y1
    real(dl), intent(out) :: out(ThisSources%SourceNum)

    integer :: i, isgn, step_ix, src_ix
    integer :: start_n, nIntSteps, nSubSteps, last_ix
    integer :: numerov_step_count, rebootstrap_interval
    integer :: K

    real(dl) :: IntAccuracyBoost, near_flat_boost, numerov_step_boost
    real(dl) :: num1, num2, scalel
    real(dl) :: dchisource, dchimax, delchi, sgndelchi, sgn
    real(dl) :: nu2, ap1, minphi
    real(dl) :: chi_end, chi_max_range, chiDispTop

    real(dl) :: sh, ch, sh_new, ch_new
    real(dl) :: step_c, step_s, step_s_for_ch
    real(dl) :: h2_12
    real(dl) :: phi
    real(dl) :: q_now, q_prev, q_next
    real(dl) :: y_prev, y_next, y2_est

    real(dl) :: alpha_rebootstrap

    real(dl), parameter :: MINPHI1 = 0.5e-4_dl
    integer, parameter :: numerov_rebootstrap_steps = 200
    integer, parameter :: numerov_low_l_rebootstrap_steps = 25
    integer, parameter :: numerov_very_low_l_rebootstrap_steps = 25

    real(dl) :: phi_vals(abs(nstart - nend) + 1)

    logical :: need_y2
    logical :: use_near_flat_gate, use_shifted_nu_near_flat
    logical :: do_rebootstrap

    out = 0._dl
    last_ix = 0

    IntAccuracyBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%NonFlatIntAccuracyBoost
    near_flat_boost = max(IntAccuracyBoost, 1._dl)

    if (nend == nstart) return

    if (nstart > IV%SourceSteps .and. nend > IV%SourceSteps) then
        y1 = 0._dl
        chi = (State%tau0 - State%TimeSteps%points(nend)) / State%curvature_radius
        return
    end if

    dchisource = dtau / State%curvature_radius

    num1 = 1._dl / nu
    scalel = real(l, dl) / State%scale
    K = merge(1, -1, State%closed)

    if (scalel >= 2400._dl) then
        num2 = 2.5_dl * num1
    else if (scalel < 50._dl) then
        num2 = 0.8_dl * num1
    else
        num2 = 1.5_dl * num1

        if (AccuracyTarget > 0) then
            num2 = 1.2_dl * num2
            if (scalel >= 700._dl .and. scalel < 1200._dl) num2 = 1.3_dl * num2
        end if
    end if

    if (num2 * IntAccuracyBoost < dchisource .and. &
        (.not. WantLateTime .or. UseLimber(l))) then
        y1 = 0._dl
        chi = (State%tau0 - State%TimeSteps%points(nend)) / State%curvature_radius
        return
    end if

    ! Mid-L modes need both enough source ranges retained and enough Numerov substeps.
    numerov_step_boost = 1._dl

    if (scalel < 700._dl) then
        numerov_step_boost = 1.25_dl

    else if (State%closed) then
        if (scalel < 1100._dl) then
            ! Closed models already converge with the shorter mid-L taper.
            numerov_step_boost = 1.25_dl - 0.25_dl * (scalel - 700._dl) / 400._dl
        end if

    else if (scalel < 1200._dl) then
        ! Open models need denser mid-L non-flat source integration extending
        ! through the TT-sensitive l~800-1000 band.
        numerov_step_boost = 1.25_dl - 0.25_dl * (scalel - 700._dl) / 500._dl
    end if

    start_n = nstart
    if (nstart > IV%SourceSteps .and. nend < IV%SourceSteps) then
        chi = (State%tau0 - State%TimeSteps%points(IV%SourceSteps)) / State%curvature_radius
        start_n = IV%SourceSteps
    end if

    chiDispTop = 1.0e20_dl
    if (State%closed) chiDispTop = const_pi - chiDisp

    isgn = sign(1, start_n - nend)   ! Direction of chi integration.
    sgn = real(isgn, dl)

    nIntSteps = isgn * (start_n - nend)
    if (nIntSteps <= 0) return

    minphi = MINPHI1 / real(l, dl) / IntAccuracyBoost
    need_y2 = isgn < 0 .or. State%closed

    nu2 = nu * nu
    ap1 = real(l, dl) * real(l + 1, dl)

    if (scalel < 1100._dl) then
        dchimax = 0.3_dl * num1
    else if (scalel < 1400._dl) then
        dchimax = 0.25_dl * num1 * 1.5_dl
    else
        dchimax = 0.35_dl * num1 * 1.5_dl
    end if

    dchimax = dchimax / (IntAccuracyBoost * numerov_step_boost)

    chi_end = chi + real(nIntSteps, dl) * dchisource * sgn
    chi_max_range = max(abs(chi), abs(chi_end))

    use_shifted_nu_near_flat = .false.
    use_near_flat_gate = UseNearFlatSmallChiApprox(l, nu, chi_max_range, near_flat_boost)

    if (use_near_flat_gate) then
        ! In the local chi << 1 regime, fill source-grid values directly from
        ! the near-flat approximation instead of stepping the full ODE.
        use_shifted_nu_near_flat = enable_shifted_nu_scalar_approx .and. &
            UseShiftedNuNearFlatIntegration(l, nu, chi_max_range, near_flat_boost)

        use_near_flat_gate = use_shifted_nu_near_flat .or. &
            enable_near_flat_smallchi_integration
    end if

    if (use_near_flat_gate) then
        if (use_shifted_nu_near_flat) then
            call FillShiftedNuPhiVals(j, l, nu, chi, dchisource, sgn, nIntSteps, &
                phi_vals(1:nIntSteps + 1), y1, chi)
        else
            call FillSmallChiPhiVals(j, l, nu, chi, dchisource, sgn, nIntSteps, &
                phi_vals(1:nIntSteps + 1), y1, chi)
        end if

        last_ix = nIntSteps + 1

    else
        if (State%closed) then
            sh = sin(chi)
            ch = cos(chi)
        else
            sh = sinh(chi)
            ch = cosh(chi)
        end if

        if (y1 == 0._dl) y1 = u_olver(l, K, nu, chi)

        phi = y1 / sh
        phi_vals(1) = phi

        nSubSteps = max(1, ceiling(dchisource / dchimax))

        if (nSubSteps >= direct_olver_min_substeps) then
            call FillDirectOlverPhiVals(l, K, nu, dchisource, sgn, isgn, nIntSteps, &
                chiDispTop, minphi, chi, y1, phi_vals, last_ix)

        else
            delchi = dchisource / real(nSubSteps, dl)
            sgndelchi = sgn * delchi
            h2_12 = sgndelchi * sgndelchi / 12._dl

            ! Step sh,ch by angle addition to avoid sin/sinh calls in the hot loop.
            if (State%closed) then
                step_c = cos(sgndelchi)
                step_s = sin(sgndelchi)
                step_s_for_ch = -step_s
            else
                step_c = cosh(sgndelchi)
                step_s = sinh(sgndelchi)
                step_s_for_ch = step_s
            end if

            q_now = ap1 / sh**2 - nu2

            rebootstrap_interval = numerov_rebootstrap_steps

            ! Low-L EE is sensitive to small Numerov phase drift when nu/l is
            ! moderately large. There are too few source-weighted oscillations
            ! for the integration to average it down, so re-sync more often.
            alpha_rebootstrap = nu / max(1._dl, real(l, dl))

            if (scalel < 25._dl .or. &
                (nu < 1000._dl .and. alpha_rebootstrap > 2._dl)) then

                rebootstrap_interval = numerov_low_l_rebootstrap_steps

                if (AccuracyTarget > 0 .and. scalel < 25._dl) then
                    rebootstrap_interval = numerov_very_low_l_rebootstrap_steps
                end if
            end if

            rebootstrap_interval = max(1, nint(rebootstrap_interval / IntAccuracyBoost))
            numerov_step_count = 0

            do i = 1, nIntSteps
                do step_ix = 1, nSubSteps
                    numerov_step_count = numerov_step_count + 1

                    do_rebootstrap = i == 1 .and. step_ix == 1
                    do_rebootstrap = do_rebootstrap .or. &
                        mod(numerov_step_count, rebootstrap_interval) == 0

                    if (do_rebootstrap) then
                        if (i == 1 .and. step_ix == 1) then
                            y_prev = y1
                        else
                            y_prev = u_olver(l, K, nu, chi)
                        end if
                        q_prev = q_now
                    end if

                    chi = chi + sgndelchi

                    sh_new = sh * step_c + ch * step_s
                    ch_new = ch * step_c + sh * step_s_for_ch
                    sh = sh_new
                    ch = ch_new

                    q_next = ap1 / sh**2 - nu2

                    if (do_rebootstrap) then
                        y1 = u_olver(l, K, nu, chi)
                        q_now = q_next
                    else
                        y_next = (2._dl * (1._dl + 5._dl * h2_12 * q_now) * y1 - &
                            (1._dl - h2_12 * q_prev) * y_prev) / &
                            (1._dl - h2_12 * q_next)

                        y_prev = y1
                        q_prev = q_now
                        y1 = y_next
                        q_now = q_next
                    end if

                    phi = y1 / sh

                    if (need_y2) then
                        ! Two-point backward-difference proxy for dy/dchi.
                        ! Only the sign matters for the cut-off.
                        y2_est = (y1 - y_prev) / sgndelchi

                        if ((isgn < 0 .and. y1 * y2_est < 0._dl) .or. &
                            (chi > chiDispTop .and. &
                            (chi > 3.14_dl .or. y1 * y2_est > 0._dl))) then
                            chi = 0._dl
                            exit
                        end if
                    end if

                    if (((isgn < 0) .or. chi > chiDispTop) .and. abs(phi) < minphi) then
                        chi = 0._dl
                        exit
                    end if
                end do

                if (chi == 0._dl) exit
                phi_vals(i + 1) = phi
            end do

            last_ix = i
        end if
    end if

    if (last_ix > 0) then
        phi_vals(1) = 0.5_dl * phi_vals(1)
        if (last_ix > 1) phi_vals(last_ix) = 0.5_dl * phi_vals(last_ix)

        src_ix = start_n

        if (ThisSources%SourceNum == 3) then
            do step_ix = 1, last_ix
                phi = phi_vals(step_ix)

                out(1) = out(1) + phi * IV%Source_q(src_ix, 1)
                out(2) = out(2) + phi * IV%Source_q(src_ix, 2)
                out(3) = out(3) + phi * IV%Source_q(src_ix, 3)

                src_ix = src_ix - isgn
            end do
        else
            do step_ix = 1, last_ix
                out = out + phi_vals(step_ix) * IV%Source_q(src_ix, :)
                src_ix = src_ix - isgn
            end do
        end if
    end if

    out = out * dchisource * State%curvature_radius

    end subroutine DoRangeInt

    subroutine FillDirectOlverPhiVals(l, K, nu, dchisource, sgn, isgn, nIntSteps, &
        chiDispTop, minphi, chi, y1, phi_vals, last_ix)
    ! Fill phi at each source time step by direct Olver evaluation, rather than Numerov
    ! stepping the hyperspherical Bessel ODE on a finer chi grid. The caller has already
    ! set phi_vals(1) from y1 at the incoming chi; chi=0 on exit signals the cut-off.
    ! The Numerov loop's derivative-sign test is not needed here (it guarded against
    ! Numerov instability in the dissipative tail); the smallness cut alone suffices.
    integer, intent(in) :: l, K, isgn, nIntSteps
    real(dl), intent(in) :: nu, dchisource, sgn, chiDispTop, minphi
    real(dl), intent(inout) :: chi
    real(dl), intent(out) :: y1
    real(dl), intent(out) :: phi_vals(:)
    integer, intent(out) :: last_ix
    integer i
    real(dl) phi

    do i = 1, nIntSteps
        chi = chi + dchisource*sgn
        phi = phi_olver(l, K, nu, chi)
        if (((isgn<0).or.(chi>chiDispTop)).and.(abs(phi) < minphi)) then
            chi=0._dl
            exit !break when getting exponentially small in dissipative region
        end if
        phi_vals(i + 1) = phi
    end do
    last_ix = i

    if (chi /= 0._dl) then
        !u value at the final chi for Numerov continuation into the next range
        if (State%closed) then
            y1 = phi*sin(chi)
        else
            y1 = phi*sinh(chi)
        end if
    else
        y1 = 0._dl
    end if

    end subroutine FillDirectOlverPhiVals

    subroutine DoRangeIntTensor(IV,chi,chiDisp,nstart,nend,dtau,nu,l,y1,out)
    ! Tensor analogue of DoRangeInt: integrates source*phi for steps nstart to nend.
    ! phi is obtained from Numerov stepping of the hyperspherical Bessel ODE,
    ! bootstrapped by two Olver evaluations (at chi and chi+delchi).
    ! dtau is the spacing of the timesteps (they must be equally spaced).
    ! (The direct-Olver high-substep path in DoRangeInt is not needed here: the denser
    ! tensor time sampling and lower tensor k range keep nSubSteps well below the gate.)

    type(IntegrationVars) :: IV
    integer l,nIntSteps,nstart,nend,isgn,i,step_ix,nSubSteps,src_ix,last_ix,K
    real(dl) nu,dtau,num1,num2,Deltachi,delchi,scalel
    real(dl) nu2,chi,chiDisp,chiDispTop

    real(dl) q_now,y1,ap1,sh,phi
    real(dl) dchimax,dchisource,sgn,sgndelchi,minphi
    real(dl) y_prev, y_next, y2_est, q_prev, q_next
    real(dl), parameter:: MINPHI1 = 1.D-6  !cut-off point for small phi l=1
    real(dl) out(ThisSources%SourceNum), phi_vals(abs(nstart - nend) + 1)
    real(dl) IntAccuracyBoost
    real(dl) h2_12, ch, step_c, step_s, step_s_for_ch, sh_new, ch_new
    logical need_y2


    if (nend==nstart) then
        out=0
        return
    end if

    IntAccuracyBoost=CP%Accuracy%AccuracyBoost*CP%Accuracy%NonFlatIntAccuracyBoost
    minphi=MINPHI1/l/IntAccuracyBoost
    isgn=sign(1,nstart-nend)!direction of chi integration
    !higher n, later time, smaller chi

    if (State%closed) then
        !Need to cut off when phi gets exponentially small as it approaches Pi
        chiDispTop = const_pi - chiDisp
    else
        chiDispTop = 1d20
    end if

    num1=1._dl/nu
    dchisource=dtau/State%curvature_radius

    scalel=l/State%scale
    if (scalel>=2000) then
        num2=num1*4
    else if (scalel>=1000) then
        num2=num1*2.5
    else if (scalel< 75) then
        num2=num1*0.1_dl
    else if (scalel<180) then
        num2=num1*0.3_dl
    else if (scalel < 600) then
        num2=num1*0.8_dl
    else
        num2=num1
    end if

    if ((isgn==1).and.(num2*IntAccuracyBoost < dchisource)) then  !Oscillating fast
        out = 0
        y1=0._dl !So we know to calculate starting y1 if there is next range
        chi=(State%tau0-State%TimeSteps%points(nend))/State%curvature_radius
        return
    end if

    sgn=isgn
    need_y2 = isgn < 0 .or. State%closed  !y2_est is only consulted by the dissipative-region cut-off

    nu2=nu*nu
    ap1=l*(l+1)
    K = merge(1, -1, State%closed)

    if (State%closed) then
        sh = sin(chi); ch = cos(chi)
    else
        sh = sinh(chi); ch = cosh(chi)
    end if

    if (scalel < 120) then
        dchimax=0.6_dl*num1
    else if (scalel < 1400) then
        dchimax=0.25_dl*num1
    else
        dchimax=0.35_dl*num1
    end if

    dchimax=dchimax/IntAccuracyBoost

    if (y1 == 0._dl) y1 = u_olver(l, K, nu, chi)
    phi=y1/sh
    out = 0
    phi_vals(1) = phi

    nIntSteps=isgn*(nstart-nend)
    if (nIntSteps <= 0) return

    nSubSteps = max(1, int(dchisource/dchimax + 0.99_dl))
    delchi = dchisource/nSubSteps


    sgndelchi=delchi*sgn
    h2_12 = sgndelchi*sgndelchi / 12._dl
    ! Angle-addition step factors: advance sh,ch by sgndelchi per substep
    ! without calling sin/sinh in the hot loop.
    if (State%closed) then
        step_c = cos(sgndelchi)
        step_s = sin(sgndelchi)
        step_s_for_ch = -step_s
    else
        step_c = cosh(sgndelchi)
        step_s = sinh(sgndelchi)
        step_s_for_ch = step_s
    end if
    q_now=(ap1/sh**2 - nu2)

    do i=1,nIntSteps
        do step_ix = 1, nSubSteps
            if (i == 1 .and. step_ix == 1) then
                !Bootstrap Numerov with a second Olver evaluation at chi+sgndelchi
                y_prev = y1
                q_prev = q_now
                chi = chi + sgndelchi
                sh_new = sh*step_c + ch*step_s
                ch_new = ch*step_c + sh*step_s_for_ch
                sh = sh_new; ch = ch_new
                q_now = ap1 / sh**2 - nu2
                y1 = u_olver(l, K, nu, chi)
            else
                chi = chi + sgndelchi
                sh_new = sh*step_c + ch*step_s
                ch_new = ch*step_c + sh*step_s_for_ch
                sh = sh_new; ch = ch_new
                q_next = ap1 / sh**2 - nu2
                y_next = (2._dl * (1._dl + 5._dl * h2_12 * q_now) * y1 - &
                    (1._dl - h2_12 * q_prev) * y_prev) / &
                    (1._dl - h2_12 * q_next)
                y_prev = y1
                q_prev = q_now
                y1 = y_next
                q_now = q_next
            end if

            phi=y1/sh
            if (need_y2) then
                !2-point backward-difference proxy for dy/dchi; only the sign matters for the cut-off
                y2_est = (y1 - y_prev) / sgndelchi
                if ((isgn<0 .and. y1*y2_est<0._dl) .or. &
                    (chi>chiDispTop .and. (chi>3.14 .or. y1*y2_est>0))) then
                    chi=0._dl
                    exit   !exit because phi now small
                end if
            end if

            if (((isgn<0).or.(chi>chiDispTop)).and.(abs(phi) < minphi)) then
                chi=0._dl
                exit  !break when getting exponentially small in dissipative region
            end if
        end do

        if (chi==0._dl) exit
        phi_vals(i + 1) = phi
    end do

    last_ix = i
    if (last_ix > 0) then
        phi_vals(1) = 0.5_dl * phi_vals(1)
        if (last_ix > 1) phi_vals(last_ix) = 0.5_dl * phi_vals(last_ix)
        src_ix  = nstart
        do step_ix = 1, last_ix
            phi = phi_vals(step_ix)
            out(1) = out(1) + phi * IV%Source_q(src_ix,1)
            out(2) = out(2) + phi * IV%Source_q(src_ix,2)
            out(3) = out(3) + phi * IV%Source_q(src_ix,3)
            src_ix = src_ix - isgn
        end do
    end if

    out = out*dchisource*State%curvature_radius

    end subroutine DoRangeIntTensor

    subroutine GetInitPowerArrayVec(pows,ks, numks)
    integer, intent(in) :: numks
    real(dl) pows(numks), ks(numks)
    integer i

    do i = 1, numks
        !!change to vec...
        pows(i) =  CP%InitPower%ScalarPower(ks(i))
        if (global_error_flag/=0) exit
    end do

    end subroutine GetInitPowerArrayVec


    subroutine GetInitPowerArrayTens(pows,ks, numks)
    integer, intent(in) :: numks
    real(dl) pows(numks), ks(numks)
    integer i

    do i = 1, numks
        pows(i) =  CP%InitPower%TensorPower(ks(i))
        if (global_error_flag/=0) exit
    end do

    end subroutine GetInitPowerArrayTens


    subroutine CalcScalCls(CTrans)
    use Bispectrum
    implicit none
    Type(ClTransferData) :: CTrans
    integer j, q_ix, w_ix, w_ix2
    real(dl) apowers
    real(dl) dlnk, ell, ctnorm, dbletmp, Delta1, Delta2
    real(dl), allocatable :: ks(:), dlnks(:), pows(:)
    real(dl) fac(3 + State%num_redshiftwindows + State%CP%CustomSources%num_custom_sources)
    integer nscal, i

    allocate(ks(CTrans%q%npoints),dlnks(CTrans%q%npoints), pows(CTrans%q%npoints))
    do q_ix = 1, CTrans%q%npoints
        if (State%flat) then
            ks(q_ix) = CTrans%q%points(q_ix)
            dlnks(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
        else
            ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - State%curv)
            dlnks(q_ix) = CTrans%q%dpoints(q_ix)*CTrans%q%points(q_ix)/ks(q_ix)**2
        end if
        pows(q_ix) = CP%InitPower%ScalarPower(ks(q_ix))
        if (global_error_flag/=0) return
    end do

#ifndef __INTEL_COMPILER
    !Can't use OpenMP here on ifort. Does not terminate.
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,4), &
    !$OMP PRIVATE(ell,q_ix,dlnk,apowers,ctnorm,dbletmp,Delta1,Delta2,w_ix,w_ix2,fac, nscal, i)
#endif
    do j=1,CTrans%ls%nl
        !Integrate dk/k Delta_l_q**2 * Power(k)
        ell = real(CTrans%ls%l(j),dl)
        if (j<= CTrans%max_index_nonlimber) then
            do q_ix = 1, CTrans%q%npoints
                if (.not.(State%closed.and.nint(CTrans%q%points(q_ix)*State%curvature_radius)<=CTrans%ls%l(j))) then
                    !cut off at nu = l + 1
                    dlnk = dlnks(q_ix)
                    apowers = pows(q_ix)

                    iCl_scalar(j,C_Temp:C_E) = iCl_scalar(j,C_Temp:C_E) +  &
                        apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)**2*dlnk
                    iCl_scalar(j,C_Cross) = iCl_scalar(j,C_Cross) + &
                        apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk

                    if (CTrans%NumSources>2 .and. State%CP%want_cl_2D_array) then

                        do w_ix=1,3 + State%num_redshiftwindows
                            Delta1= CTrans%Delta_p_l_k(w_ix,j,q_ix)
                            if (w_ix>3) then
                                associate (Win => State%Redshift_w(w_ix - 3))
                                    if (Win%kind == window_lensing) &
                                        Delta1 = Delta1 / 2 * ell * (ell + 1)
                                    if (Win%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                                        !want delta f/f - 2kappa;
                                        ! grad^2 = -l(l+1);
                                        Delta1 = Delta1 + ell * (ell + 1) * &
                                            CTrans%Delta_p_l_k(3 + Win%mag_index + &
                                            State%num_redshiftwindows, j, q_ix)
                                    end if
                                end associate
                            end if
                            do w_ix2=w_ix,3 + State%num_redshiftwindows
                                if (w_ix2>= 3.and. w_ix>=3) then
                                    !Skip if the auto or cross-correlation is included in direct Limber result
                                    !Otherwise we need to include the sources e.g. to get counts-Temperature correct
                                    if (CTrans%limber_l_min(w_ix2)/= 0 .and. j>=CTrans%limber_l_min(w_ix2) &
                                        .and. CTrans%limber_l_min(w_ix)/= 0 .and. j>=CTrans%limber_l_min(w_ix)) cycle

                                end if
                                Delta2 = CTrans%Delta_p_l_k(w_ix2, j, q_ix)
                                if (w_ix2 > 3) then
                                    associate (Win => State%Redshift_w(w_ix2 - 3))
                                        if (Win%kind == window_lensing) &
                                            Delta2 = Delta2 / 2 * ell * (ell + 1)
                                        if (Win%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                                            !want delta f/f - 2kappa;
                                            ! grad^2 = -l(l+1);
                                            Delta2 = Delta2 + ell * (ell + 1) * &
                                                CTrans%Delta_p_l_k(3 + Win%mag_index + &
                                                State%num_redshiftwindows, j, q_ix)
                                        end if
                                    end associate
                                end if
                                iCl_Array(j,w_ix,w_ix2) = iCl_Array(j,w_ix,w_ix2)+Delta1*Delta2*apowers*dlnk
                            end do
                        end do
                        if (CP%CustomSources%num_custom_sources >0) then
                            do w_ix=1,3 + State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                                if (w_ix > 3 + State%num_redshiftwindows) then
                                    Delta1= CTrans%Delta_p_l_k(w_ix+State%num_extra_redshiftwindows,j,q_ix)
                                else
                                    Delta1= CTrans%Delta_p_l_k(w_ix,j,q_ix)
                                end if
                                do w_ix2=max(w_ix,3 + State%num_redshiftwindows +1), &
                                    3 + State%num_redshiftwindows +CP%CustomSources%num_custom_sources
                                    Delta2=  CTrans%Delta_p_l_k(w_ix2+State%num_extra_redshiftwindows,j,q_ix)
                                    iCl_Array(j,w_ix,w_ix2) = iCl_Array(j,w_ix,w_ix2) &
                                        +Delta1*Delta2*apowers*dlnk
                                end do
                            end do
                        end if
                    end if

                    if (CTrans%NumSources>2 ) then
                        if (CP%SourceTerms%limber_phi_lmin==0 .or.  &
                            CTrans%limber_l_min(3)== 0 .or. j<CTrans%limber_l_min(3)) then
                            iCl_scalar(j,C_Phi) = iCl_scalar(j,C_Phi) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)**2*dlnk
                            iCl_scalar(j,C_PhiTemp) = iCl_scalar(j,C_PhiTemp) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
                            iCl_scalar(j,C_PhiE) = iCl_scalar(j,C_PhiE) +  &
                                apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
                        end if
                    end if
                end if
            end do

        end if !limber (j<= max_bessels_l_index)

        !Output l(l+1)C_l/OutputDenominator
        ctnorm=(ell*ell-1)*(ell+2)*ell
        dbletmp=(ell*(ell+1))/OutputDenominator*const_fourpi
        if (State%CP%want_cl_2D_array) then
            fac=1
            fac(2) = sqrt(ctnorm)
            if (CTrans%NumSources > 2) then
                fac(3) = sqrt(ell*(ell+1)*CP%ALens) !Changed Dec18 for consistency
                do w_ix=3 + State%num_redshiftwindows+1,3 + State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                    nscal= CP%CustomSources%custom_source_ell_scales(w_ix - State%num_redshiftwindows -3)
                    do i=1, nscal
                        fac(w_ix) = fac(w_ix)*(ell+i)*(ell-i+1)
                    end do
                    fac(w_ix) = sqrt(fac(w_ix))
                end do
            end if

            do w_ix=1, CTrans%NumSources - State%num_extra_redshiftwindows
                do w_ix2=w_ix,CTrans%NumSources - State%num_extra_redshiftwindows
                    iCl_Array(j,w_ix,w_ix2) =iCl_Array(j,w_ix,w_ix2) &
                        *fac(w_ix)*fac(w_ix2)*dbletmp
                    iCl_Array(j,w_ix2,w_ix) = iCl_Array(j,w_ix,w_ix2)
                end do
            end do
        end if

        iCl_scalar(j,C_Temp)  =  iCl_scalar(j,C_Temp)*dbletmp
        iCl_scalar(j,C_E) =  iCl_scalar(j,C_E)*dbletmp*ctnorm
        iCl_scalar(j,C_Cross) =  iCl_scalar(j,C_Cross)*dbletmp*sqrt(ctnorm)
        if (CTrans%NumSources>2) then
            iCl_scalar(j,C_Phi) = CP%ALens*iCl_scalar(j,C_Phi)*const_fourpi*ell**4
            !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
            !We put pix extra factors of l here to improve interpolation in l
            iCl_scalar(j,C_PhiTemp) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiTemp)*const_fourpi*ell**3
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
            iCl_scalar(j,C_PhiE) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiE)*const_fourpi*ell**3*sqrt(ctnorm)
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi E}
        end if
    end do
#ifndef __INTEL_COMPILER
    !$OMP END PARALLEL DO
#endif

    end subroutine CalcScalCls

    subroutine CalcScalCls2(CTrans)
    !Calculate C_ll' for non-isotropic models
    !Run with l_sample_boost=50 to get every l
    !not used in normal CAMB
    use FileUtils
    implicit none
    Type(ClTransferData) :: CTrans
    integer j,j2
    real(dl) apowers, pows(CTrans%q%npoints)
    integer q_ix
    real(dl)  ks(CTrans%q%npoints),dlnks(CTrans%q%npoints),dlnk
    real(dl) ctnorm,dbletmp
    real(dl), allocatable :: iCl_Scalar2(:,:,:)
    type(TTextFile) :: F

    allocate(iCl_Scalar2(CTranS%ls%nl,CTrans%ls%nl,C_Temp:State%Scalar_C_last))
    iCl_scalar2 = 0

    do q_ix = 1, CTrans%q%npoints
        if (State%flat) then
            ks(q_ix) = CTrans%q%points(q_ix)
            dlnks(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
        else
            ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - State%curv)
            dlnks(q_ix) = CTrans%q%dpoints(q_ix)*CTrans%q%points(q_ix)/ks(q_ix)**2
        end if

        pows(q_ix) =  CP%InitPower%ScalarPower(ks(q_ix))
        if (global_error_flag/=0) return
    end do

    do j=1,CTrans%ls%nl
        do j2=1,CTrans%ls%nl
            !Integrate dk/k Delta_l_q**2 * Power(k)
            do q_ix = 1, CTrans%q%npoints
                if (.not.(State%closed.and.nint(CTrans%q%points(q_ix)*State%curvature_radius)<= CTrans%ls%l(j))) then
                    !cut off at nu = l + 1
                    dlnk = dlnks(q_ix)
                    apowers = pows(q_ix)

                    iCl_scalar2(j,j2,C_Temp:C_E) = iCl_scalar2(j,j2,C_Temp:C_E) +  &
                        apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)*CTrans%Delta_p_l_k(1:2,j2,q_ix)*dlnk
                    iCl_scalar2(j,j2,C_Cross) = iCl_scalar2(j,j2,C_Cross) + &
                        apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j2,q_ix)*dlnk
                end if
            end do

            !Output l(l+1)C_l/OutputDenominator

            !ctnorm = (CTrans%ls%l+2)!/(CTrans%ls%l-2)! - beware of int overflow
            ctnorm=(CTrans%ls%l(j)*CTrans%ls%l(j)-1)*real((CTrans%ls%l(j)+2)*CTrans%ls%l(j),dl)
            ctnorm=sqrt(ctnorm*(CTrans%ls%l(j2)*CTrans%ls%l(j2)-1)*real((CTrans%ls%l(j2)+2)*CTrans%ls%l(j2),dl))

            dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*const_fourpi
            dbletmp=sqrt(dbletmp*(CTrans%ls%l(j2)*(CTrans%ls%l(j2)+1))/OutputDenominator*const_fourpi)

            iCl_scalar2(j,j2,C_Temp)  =  iCl_scalar2(j,j2,C_Temp)*dbletmp
            iCl_scalar2(j,j2,C_E)     =  iCl_scalar2(j,j2,C_E)*dbletmp*ctnorm
            iCl_scalar2(j,j2,C_Cross) =  iCl_scalar2(j,j2,C_Cross)*dbletmp*sqrt(ctnorm)
        end do
    end do

    call F%CreateFile('cl2.dat')
    do j=1,CTrans%ls%nl
        do j2=1,CTrans%ls%nl
            call F%write(CTrans%ls%l(j),CTrans%ls%l(j2),iCl_scalar2(j,j2,1)*7.4311e12)
        end do
    end do
    call F%close()
    call F%CreateFile('cl1l2.dat')
    do j=1,999
        call F%write(iCl_scalar2(j,1:999,1)*7.4311e12)
    end do
    call F%close()
    stop

    end subroutine CalcScalCls2


    subroutine CalcTensCls(CTrans, GetInitPowers)
    implicit none
    Type(ClTransferData) :: CTrans
    external GetInitPowers
    integer j, q_ix
    real(dl) nu
    real(dl) apowert,  measure
    real(dl) ctnorm,dbletmp
    real(dl) pows(CTrans%q%npoints)
    real(dl)  ks(CTrans%q%npoints),measures(CTrans%q%npoints)

    !For tensors we want Integral dnu/nu (nu^2-3)/(nu^2-1) Delta_l_k^2 P(k) for State%closed

    do q_ix = 1, CTrans%q%npoints
        if (State%flat) then
            ks(q_ix) = CTrans%q%points(q_ix)
            measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
        else
            nu = CTrans%q%points(q_ix)*State%curvature_radius
            ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - 3*State%curv)
            measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)*(nu**2-3*State%Ksign)/(nu**2-State%Ksign)
        end if
    end do

    call GetInitPowers(pows,ks,CTrans%q%npoints)

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC,4) &
    !$OMP & PRIVATE(q_ix, measure, apowert, ctnorm, dbletmp)
    do j=1,CTrans%ls%nl
        do q_ix = 1, CTrans%q%npoints
            if (.not.(State%closed.and. nint(CTrans%q%points(q_ix)*State%curvature_radius)<=CTrans%ls%l(j))) then
                !cut off at nu = l+1
                apowert = pows(q_ix)
                measure = measures(q_ix)

                iCl_tensor(j,CT_Temp:CT_B) = iCl_tensor(j,CT_Temp:CT_B) + &
                    apowert*CTrans%Delta_p_l_k(CT_Temp:CT_B,j,q_ix)**2*measure

                iCl_tensor(j,CT_cross ) = iCl_tensor(j,CT_cross) &
                    +apowert*CTrans%Delta_p_l_k(CT_Temp,j,q_ix)*CTrans%Delta_p_l_k(CT_E,j,q_ix)*measure
            end if
        end do

        ctnorm=(CTrans%ls%l(j)*CTrans%ls%l(j)-1)*real((CTrans%ls%l(j)+2)*CTrans%ls%l(j),dl)
        dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*const_pi/4
        iCl_tensor(j, CT_Temp) = iCl_tensor(j, CT_Temp)*dbletmp*ctnorm
        if (CTrans%ls%l(j)==1) dbletmp=0
        iCl_tensor(j, CT_E:CT_B) = iCl_tensor(j, CT_E:CT_B)*dbletmp
        iCl_tensor(j, CT_Cross)  = iCl_tensor(j, CT_Cross)*dbletmp*sqrt(ctnorm)
    end do
    !$OMP END PARALLEL DO

    end subroutine CalcTensCls


    subroutine CalcVecCls(CTrans, GetInitPowers)
    implicit none
    Type(ClTransferData) :: CTrans
    external GetInitPowers
    integer j, q_ix
    real(dl) power,  measure
    real(dl) ctnorm,lfac,dbletmp
    real(dl) pows(CTrans%q%npoints)
    real(dl)  ks(CTrans%q%npoints),measures(CTrans%q%npoints)

    do q_ix = 1, CTrans%q%npoints
        ks(q_ix) = CTrans%q%points(q_ix)
        measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
    end do

    call GetInitPowers(pows,ks,CTrans%q%npoints)

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC,4), &
    !$OMP & PRIVATE(j,q_ix,measure,power,ctnorm,dbletmp,lfac)
    do j=1,CTrans%ls%nl
        do q_ix = 1, CTrans%q%npoints
            if (.not.(State%closed.and. nint(CTrans%q%points(q_ix)*State%curvature_radius)<=CTrans%ls%l(j))) then
                !cut off at nu = l+1
                power = pows(q_ix)
                measure = measures(q_ix)

                iCl_vector(j,CT_Temp:CT_B) = iCl_vector(j,CT_Temp:CT_B) + &
                    power*CTrans%Delta_p_l_k(CT_Temp:CT_B,j,q_ix)**2*measure

                iCl_vector(j,CT_cross) = iCl_vector(j,CT_cross) &
                    +power*CTrans%Delta_p_l_k(CT_Temp,j,q_ix)*CTrans%Delta_p_l_k(CT_E,j,q_ix)*measure
            end if
        end do

        ctnorm=CTrans%ls%l(j)*(CTrans%ls%l(j)+1)
        dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*const_pi/8
        iCl_vector(j, CT_Temp)   = iCl_vector(j, CT_Temp)*dbletmp*ctnorm
        lfac = (CTrans%ls%l(j) + 2)*(CTrans%ls%l(j) - 1)
        iCl_vector(j, CT_E:CT_B) = iCl_vector(j, CT_E:CT_B)*dbletmp*lfac
        iCl_vector(j, CT_Cross)  = iCl_vector(j, CT_Cross)*dbletmp*sqrt(lfac*ctnorm)
    end do
    !$OMP END PARALLEL DO

    end subroutine CalcVecCls


    subroutine InterpolateCls()
    implicit none
    integer i,j
    integer, dimension(2,2), parameter :: ind = reshape( (/ 1,3,3,2 /), shape(ind))
    !use verbose form above for gfortran consistency  [[1,3],[3,2]]

    !Note using log interpolation is worse

    if (CP%WantScalars) then
        associate(lSet=>State%CLdata%CTransScal%ls)
            do i = C_Temp, State%Scalar_C_last
                call lSet%InterpolateClArrTemplated(iCl_scalar(1,i),State%CLData%Cl_scalar(lSet%lmin, i), &
                    lSet%nl,i)
            end do

            if (State%CLdata%CTransScal%NumSources>2 .and. State%CP%want_cl_2D_array) then
                do i=1,3+State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                    do j=i,3+State%num_redshiftwindows + CP%CustomSources%num_custom_sources
                        if (i<3 .and. j<3) then
                            State%CLData%Cl_scalar_array(:,i,j) = State%CLData%Cl_scalar(:, ind(i,j))
                        else
                            call lSet%InterpolateClArr(iCl_array(1,i,j), &
                                State%CLData%Cl_scalar_array(lSet%lmin, i,j))
                        end if
                        if (i/=j) State%CLData%Cl_scalar_array(:,j,i) = State%CLData%Cl_scalar_array(:,i,j)
                    end do
                end do
            end if
        end associate
    end if

    if (CP%WantVectors) then
        associate(ls=>State%CLdata%CTransVec%ls)
            do i = C_Temp, CT_cross
                call ls%InterpolateClArr(iCl_vector(1,i),State%CLData%Cl_vector(ls%lmin, i))
            end do
        end associate
    end if

    if (CP%WantTensors) then
        associate(ls=>State%CLdata%CTransTens%ls)
            do i = CT_Temp, CT_Cross
                call ls%InterpolateClArr(iCl_tensor(1,i),State%CLData%Cl_tensor(ls%lmin, i))
            end do
        end associate
    end if

    end subroutine InterpolateCls

    subroutine CAMBMain_Free()
    ! Not needed, only used to clean up for mem leak testing
    if (allocated(iCl_scalar)) deallocate(iCl_scalar)
    if (allocated(iCl_vector)) deallocate(iCl_vector)
    if (allocated(iCl_tensor)) deallocate(iCl_tensor)
    if (allocated(iCl_Array)) deallocate(iCl_Array)
    if (allocated(TempSources)) deallocate(TempSources)
    end subroutine CAMBMain_Free


    end module CAMBmain
