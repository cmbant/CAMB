    !Interface module for CAMB. Call CAMB_GetResults to do the work.

    module CAMB
    use Precision
    use MathUtils, only: brentq
    use results
    use GaugeInterface
    use InitialPower
    use Reionization
    use Recombination, only : TRecFast
    use lensing
    use DarkEnergyFluid
    implicit none

    type TCAMBThetaH0Solver
        type(CAMBparams) :: Params
        type(CAMBdata) :: Data
        real(dl) :: theta = 0._dl
        real(dl) :: zstar = 0._dl
        logical :: cosmomc_approx = .false.
    end type TCAMBThetaH0Solver

    contains

    subroutine CAMB_TransfersToPowers(CData)
    !From Delta_p_l_k or time transfers to CMB powers and transfers to P(k)
    use CAMBmain
    use lensing
    type (CAMBdata) :: CData
    logical :: want_tensors, want_vectors

    call SetActiveState(CData)
    CData%OnlyTransfer = .false.
    call CData%CP%InitPower%Init(CData%CP)
    if (global_error_flag/=0) return
    if (allocated(Cdata%CAMB_Pk)) deallocate(Cdata%CAMB_PK)

    if (CData%CP%WantCls) then
        if (allocated(CData%ScalarTimeSources) .and. CData%CP%WantScalars) then
            want_tensors = CData%CP%WantTensors
            want_vectors = CData%CP%WantVectors
            Cdata%OnlyTransfer = .true. !prevent ClTransferToCl
            Cdata%CP%WantTensors = .false.
            CData%CP%WantVectors = .false.
            call TimeSourcesToCl(CData%ClData%CTransScal)
            Cdata%CP%WantTensors = want_tensors
            CData%CP%WantVectors = want_vectors
            Cdata%OnlyTransfer = .false.
        end if
        call ClTransferToCl(CData)
        if (CData%CP%DoLensing .and. global_error_flag==0) call lens_Cls(Cdata)
        if (global_error_flag/=0) return
    end if

    if (CData%CP%WantTransfer) call Transfer_Get_sigmas(Cdata, Cdata%MT)

    end subroutine CAMB_TransfersToPowers

    !Call this routine with a set of parameters to generate the results you want.
    subroutine CAMB_GetResults(OutData, Params, error, onlytransfer, onlytimesources, &
        BispectrumConfig, BispectrumResult, bispectrum_output_root)
    use CAMBmain
    use lensing
    use Bispectrum
    type(CAMBdata)  :: OutData
    type(CAMBparams) :: Params
    integer, optional :: error !Zero if OK
    logical, optional :: onlytransfer, onlytimesources
    Type(TBispectrumParams), optional, intent(in) :: BispectrumConfig
    Type(TBispectrumResult), optional :: BispectrumResult
    character(LEN=*), optional, intent(in) :: bispectrum_output_root
    type(CAMBparams) P
    logical :: call_again, calc_bispectrum, old_do_bispectrum

    global_error_flag = 0
    old_do_bispectrum = do_bispectrum
    if (present(BispectrumConfig)) then
        calc_bispectrum = BispectrumConfig%do_lensing_bispectrum .or. BispectrumConfig%do_primordial_bispectrum
        do_bispectrum = calc_bispectrum
    else
        calc_bispectrum = do_bispectrum
    end if
    call_again = .false.
    call OutData%Free()
    call SetActiveState(OutData)
    OutData%HasScalarTimeSources= DefaultFalse(onlytimesources)
    OutData%OnlyTransfer = DefaultFalse(onlytransfer) .or. OutData%HasScalarTimeSources

    !Vector and tensors first, so at end time steps in OutData are for scalars
    if (Params%WantCls .and. Params%WantTensors) then
        P=Params
        P%WantTransfer = .false.
        P%Transfer%high_precision = .false.
        P%WantScalars = .false.
        P%WantVectors = .false.
        call OutData%SetParams(P, call_again=call_again)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            do_bispectrum = old_do_bispectrum
            return
        end if
        call_again = .true.
    end if

    if (Params%WantCls .and. Params%WantVectors) then
        P=Params
        P%WantTransfer = .false.
        P%Transfer%high_precision = .false.
        P%WantScalars = .false.
        P%WantTensors = .false.
        call OutData%SetParams(P, call_again=call_again)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            do_bispectrum = old_do_bispectrum
            return
        end if
        call_again = .true.
    end if

    if (Params%WantCls .and. Params%WantScalars) then
        P = Params
        if (calc_bispectrum) P%Accuracy%lSampleBoost = max(P%Accuracy%lSampleBoost, 50._dl)
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        P%WantTensors = .false.
        P%WantVectors = .false.
        if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. &
            (P%DoLensing .or. OutData%num_redshiftwindows > 0)) then
            P%WantTransfer  = .true.
        end if
        call OutData%SetParams(P)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            do_bispectrum = old_do_bispectrum
            return
        end if
        call_again = .true.
    end if


    if (Params%WantTransfer .and. .not. (Params%WantCls .and. Params%WantScalars)) then
        P=Params
        P%WantCls = .false.
        P%WantScalars = .false.
        P%WantTensors = .false.
        P%WantVectors = .false.
        call OutData%SetParams(P, call_again=call_again)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            do_bispectrum = old_do_bispectrum
            return
        end if
    end if
    OutData%CP%WantCls = Params%WantCls
    OutData%CP%WantScalars = Params%WantScalars
    OutData%CP%WantTensors = Params%WantTensors
    OutData%CP%WantVectors = Params%WantVectors
    OutData%CP%WantTransfer = Params%WantTransfer
    OutData%CP%Accuracy = Params%Accuracy
    if (calc_bispectrum) OutData%CP%Accuracy%lSampleBoost = max(OutData%CP%Accuracy%lSampleBoost, 50._dl)
    OutData%CP%Reion%Reionization = Params%Reion%Reionization
    OutData%CP%Transfer%high_precision = Params%Transfer%high_precision
    OutData%CP%WantDerivedParameters = Params%WantDerivedParameters

    if (.not. OutData%OnlyTransfer .and. Params%WantCls .and. Params%WantScalars) then
        if (Params%DoLensing .and. global_error_flag==0) then
            call lens_Cls(OutData)
        end if

        if (calc_bispectrum .and. global_error_flag==0) then
            if (present(BispectrumConfig)) then
                if (present(bispectrum_output_root)) then
                    call GetBispectrum(OutData, OutData%CLData%CTransScal, &
                        BispectrumConfig, BispectrumResult, bispectrum_output_root)
                else
                    call GetBispectrum(OutData, OutData%CLData%CTransScal, BispectrumConfig, BispectrumResult)
                end if
            else
                call GetBispectrum(OutData,OutData%CLData%CTransScal)
            end if
        end if
    end if
    if (global_error_flag/=0 .and. present(error)) error =global_error_flag
    do_bispectrum = old_do_bispectrum

    end subroutine CAMB_GetResults


    !Return real (NOT double precision) arrays of the computed CMB  Cls
    !Output is l(l+1)C_l/2pi
    !If GC_Conventions = .false. use E-B conventions (as the rest of CAMB does)
    !Used by WriteFits only
    subroutine CAMB_GetCls(OutData, Cls, lmax,GC_conventions)
    Type(CAMBdata) :: OutData
    integer, intent(IN) :: lmax
    logical, intent(IN) :: GC_conventions
    real, intent(OUT) :: Cls(2:lmax,1:4)
    integer l

    Cls = 0
    do l=2, lmax
        if (OutData%CP%WantScalars .and. l<= OutData%CP%Max_l) then
            if (OutData%CP%DoLensing) then
                if (l<=OutData%CLData%lmax_lensed) &
                    Cls(l,1:4) = OutData%CLData%Cl_lensed(l, CT_Temp:CT_Cross)
            else
                Cls(l,1:2) = OutData%CLData%Cl_scalar(l, C_Temp:C_E)
                Cls(l,4) = OutData%CLData%Cl_scalar(l, C_Cross)
            endif
        end if
        if (OutData%CP%WantTensors .and. l <= OutData%CP%Max_l_tensor) then
            Cls(l,1:4) = Cls(l,1:4) + OutData%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
        end if
    end do
    if (GC_conventions) then
        Cls(:,2:3) = Cls(:,2:3)/2
        Cls(:,4)   = Cls(:,4)/sqrt(2.0)
    end if

    end subroutine CAMB_GetCls

    function CAMB_GetAge(P)
    !Return age in Julian gigayears, returns -1 on error
    type(CAMBparams), intent(in) :: P
    real(dl) CAMB_GetAge
    integer error
    Type(CAMBdata) :: OutData

    call  OutData%SetParams(P, error, .false., .false., .true.)

    if (error/=0) then
        CAMB_GetAge = -1
    else
        CAMB_GetAge = OutData%DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
    end if

    end function CAMB_GetAge

    subroutine CAMB_SetDefParams(P)
    use NonLinear
    type(CAMBparams), intent(inout), target :: P
    type(CAMBparams) :: emptyP

    P= emptyP !Set default values set in type definitions
    P%Nu_mass_numbers=0
    P%Nu_mass_degeneracies=0
    P%Nu_mass_fractions=0
    P%lens_output_margin = 200

    allocate(THalofit::P%NonLinearModel)
    allocate(TDarkEnergyFluid::P%DarkEnergy)
    allocate(TInitialPowerLaw::P%InitPower)
    allocate(TRecfast::P%Recomb)
    allocate(TTanhReionization::P%Reion)

    end subroutine CAMB_SetDefParams

    real(dl) function CAMB_ThetaH0Difference(obj, H0)
    class(*) :: obj
    real(dl), intent(in) :: H0
    integer :: error

    select type (solver => obj)
    type is (TCAMBThetaH0Solver)
        solver%Params%H0 = H0
        call solver%Data%SetParams(solver%Params, error=error, background_only=.true.)
        if (error /= 0) then
            CAMB_ThetaH0Difference = huge(1._dl)
            return
        end if
        if (solver%cosmomc_approx) then
            CAMB_ThetaH0Difference = solver%Data%CosmomcTheta() - solver%theta
        else
            CAMB_ThetaH0Difference = solver%Data%sound_horizon(solver%zstar) / &
                (solver%Data%AngularDiameterDistance(solver%zstar) * (1 + solver%zstar)) - solver%theta
        end if
    class default
        error stop 'CAMB_ThetaH0Difference: unexpected solver state'
    end select

    end function CAMB_ThetaH0Difference

    logical function CAMB_SetH0ForTheta(P, theta, ErrMsg, cosmomc_approx, theta_H0_range, est_H0, iteration_threshold)
    type(CAMBparams), intent(inout) :: P
    real(dl), intent(in) :: theta
    character(LEN=*), intent(out) :: ErrMsg
    logical, intent(in), optional :: cosmomc_approx
    real(dl), intent(in), optional :: theta_H0_range(2)
    real(dl), intent(in), optional :: est_H0, iteration_threshold
    type(TCAMBThetaH0Solver) :: solver
    real(dl) :: H0_range(2), initial_H0, H0_iteration_threshold
    real(dl) :: xzero, fzero
    integer :: error, iflag

    CAMB_SetH0ForTheta = .false.
    ErrMsg = ''

    H0_range = [10._dl, 100._dl]
    if (present(theta_H0_range)) H0_range = theta_H0_range
    initial_H0 = 67._dl
    if (present(est_H0)) initial_H0 = est_H0
    H0_iteration_threshold = 8._dl
    if (present(iteration_threshold)) H0_iteration_threshold = iteration_threshold

    if (.not. (theta > 0.001_dl .and. theta < 0.1_dl)) then
        ErrMsg = 'theta looks wrong (parameter is just theta, not 100*theta)'
        return
    end if

    solver%Params = P
    solver%theta = theta
    solver%cosmomc_approx = PresentDefault(.false., cosmomc_approx)
    if (.not. solver%cosmomc_approx) solver%Params%H0 = initial_H0

    call solver%Data%SetParams(solver%Params, error=error, background_only=.true.)
    if (error /= 0) then
        ErrMsg = trim(global_error_message)
        return
    end if

    if (.not. solver%cosmomc_approx) solver%zstar = solver%Data%get_zstar()

    call brentq(solver, CAMB_ThetaH0Difference, H0_range(1), H0_range(2), 5e-5_dl, xzero, fzero, iflag)
    if (iflag /= 0) then
        ErrMsg = 'No solution for H0 inside of theta_H0_range'
        return
    end if

    solver%Params%H0 = xzero
    if (.not. solver%cosmomc_approx .and. abs(xzero - initial_H0) > H0_iteration_threshold) then
        call solver%Data%SetParams(solver%Params, error=error, background_only=.true.)
        if (error /= 0) then
            ErrMsg = trim(global_error_message)
            return
        end if
        solver%zstar = solver%Data%get_zstar()
        call brentq(solver, CAMB_ThetaH0Difference, H0_range(1), H0_range(2), 5e-5_dl, xzero, fzero, iflag)
        if (iflag /= 0) then
            ErrMsg = 'No solution for H0 inside of theta_H0_range'
            return
        end if
    end if

    P%H0 = xzero
    CAMB_SetH0ForTheta = .true.

    end function CAMB_SetH0ForTheta

    logical function CAMB_ReadParamFile(P, InputFile, InpLen)
    Type(CAMBParams) :: P
    integer, intent(in) :: InpLen
    character(LEN=InpLen), intent(in) :: InputFile
    character(LEN=len(global_error_message)) :: ErrMsg
    Type(TIniFile) :: Ini
    logical bad

    call Ini%Open(InputFile, bad, .false.)
    ErrMsg = ''
    CAMB_ReadParamFile = CAMB_ReadParams(P, Ini, ErrMsg)
    call Ini%Close()
    if (ErrMsg/='') call GlobalError(ErrMsg,error_ini)

    end function CAMB_ReadParamFile

    logical function CAMB_ReadParams(P, Ini, ErrMsg)
    use NonLinear
    use SPkNonLinear, only: TSPkNonLinear
    use DarkEnergyFluid
    use DarkEnergyPPF
    use Quintessence
    use results
#ifdef COSMOREC
    use CosmoRec
#endif
#ifdef HYREC
    use HyRec
#endif
    class(TIniFile) :: Ini
    Type(CAMBParams) :: P
    integer num_redshiftwindows
    logical PK_WantTransfer
    integer i, status, num_theta_inputs
    real(dl) nmassive, theta
    character(LEN=*), intent(out) :: ErrMsg
    character(LEN=:), allocatable :: NumStr, S, DarkEneryModel, RecombinationModel
    logical :: DoCounts, has_hubble, has_thetastar, has_cosmomc_theta

    ErrMsg = ''
    global_error_flag = 0
    global_error_message = ''

    CAMB_ReadParams = .false.
    call CAMB_SetDefParams(P)

    P%WantScalars = Ini%Read_Logical('get_scalar_cls')
    P%WantVectors = Ini%Read_Logical('get_vector_cls', .false.)
    P%WantTensors = Ini%Read_Logical('get_tensor_cls', .false.)

    P%Want_CMB =  Ini%Read_Logical('want_CMB',.true.)
    P%Want_CMB_lensing =  P%Want_CMB .or. Ini%Read_Logical('want_CMB_lensing',.true.)

    if (P%WantScalars) then
        num_redshiftwindows = Ini%Read_Int('num_redshiftwindows',0)
    else
        num_redshiftwindows = 0
    end if
    call Ini%Read('limber_windows', P%SourceTerms%limber_windows)
    if (P%SourceTerms%limber_windows) call Ini%Read('limber_phiphi', P%SourceTerms%limber_phi_lmin)
    if (num_redshiftwindows > 0) then
        allocate(P%SourceWindows(num_redshiftwindows))
        P%SourceTerms%counts_lensing = Ini%Read_Logical('DoRedshiftLensing', .false.)
        call ReadAccuracyReal(P%Accuracy%KmaxBoost, 'KmaxBoost', 'Kmax_Boost')
        if (ErrMsg /= '') return
    end if
    P%Do21cm = Ini%Read_Logical('Do21cm', .false.)
    DoCounts = .false.
    do i=1, num_redshiftwindows
        allocate(TGaussianSourceWindow::P%SourceWindows(i)%Window)
        select type (RedWin=>P%SourceWindows(i)%Window)
        class is (TGaussianSourceWindow)
            RedWin%Redshift = Ini%Read_Double_Array('redshift', i)
            S = Ini%Read_String_Array('redshift_kind', i)
            if (S == '21cm') then
                RedWin%source_type = window_21cm
            elseif (S == 'counts') then
                RedWin%source_type = window_counts
            elseif (S == 'lensing') then
                RedWin%source_type = window_lensing
            else
                ErrMsg = 'Error: unknown type of window '//trim(S)
                return
            end if
            if (RedWin%source_type /= window_21cm) then
                RedWin%sigma = Ini%Read_Double_Array('redshift_sigma', i)
            else
                P%Do21cm = .true.
                RedWin%sigma = Ini%Read_Double_Array('redshift_sigma_Mhz', i)
                if (RedWin%sigma < 0.003 .and. print_fortran_warnings) then
                    write(*,*) 'WARNING:Window very narrow.'
                    write(*,*) ' --> use transfer functions and transfer_21cm_cl =T ?'
                end if
                !with 21cm widths are in Mhz, make dimensionless scale factor
                RedWin%sigma = RedWin%sigma / (f_21cm / 1e6)
                if (FeedbackLevel>0) write(*,*) i,'delta_z = ',  RedWin%sigma * (1 + RedWin%RedShift) ** 2
            end if
            if (RedWin%source_type == window_counts) then
                DoCounts = .true.
                RedWin%bias = Ini%Read_Double_Array('redshift_bias', i)
                RedWin%dlog10Ndm = Ini%Read_Double_Array('redshift_dlog10Ndm', i ,0.d0)
            end if
        class default
            call MpiStop('Probable compiler bug')
        end select
    end do


    if (P%Do21cm) then
        call Ini%Read('line_basic',P%SourceTerms%line_basic)
        call Ini%Read('line_distortions',P%SourceTerms%line_distortions)
        call Ini%Read('line_extra',P%SourceTerms%line_extra)
        call Ini%Read('line_phot_dipole',P%SourceTerms%line_phot_dipole)
        call Ini%Read('line_phot_quadrupole',P%SourceTerms%line_phot_quadrupole)
        call Ini%Read('line_reionization',P%SourceTerms%line_reionization)

        call Ini%Read('use_mK',P%SourceTerms%use_21cm_mK)
        if (DebugMsgs) then
            write (*,*) 'Doing 21cm'
            write (*,*) 'dipole = ',P%SourceTerms%line_phot_dipole, ' quadrupole =', P%SourceTerms%line_phot_quadrupole
        end if
    else
        P%SourceTerms%line_extra = .false.
    end if

    if (DoCounts) then
        call Ini%Read('counts_density', P%SourceTerms%counts_density)
        call Ini%Read('counts_redshift', P%SourceTerms%counts_redshift)
        call Ini%Read('counts_radial', P%SourceTerms%counts_radial)
        call Ini%Read('counts_evolve', P%SourceTerms%counts_evolve)
        call Ini%Read('counts_timedelay', P%SourceTerms%counts_timedelay)
        call Ini%Read('counts_ISW', P%SourceTerms%counts_ISW)
        call Ini%Read('counts_potential', P%SourceTerms%counts_potential)
        call Ini%Read('counts_velocity', P%SourceTerms%counts_velocity)
    end if

    P%OutputNormalization=outNone

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    PK_WantTransfer = Ini%Read_Logical('get_transfer')

    call ReadAccuracyReal(P%Accuracy%AccuracyBoost, 'AccuracyBoost', 'accuracy_boost')
    call ReadAccuracyReal(P%Accuracy%lAccuracyBoost, 'lAccuracyBoost', 'l_accuracy_boost')
    call ReadAccuracyReal(P%Accuracy%TimeStepBoost, 'TimeStepBoost')
    call ReadAccuracyReal(P%Accuracy%BackgroundTimeStepBoost, 'BackgroundTimeStepBoost')
    call ReadAccuracyReal(P%Accuracy%TimeSwitchBoost, 'TimeSwitchBoost')
    call ReadAccuracyReal(P%Accuracy%IntTolBoost, 'IntTolBoost')
    call ReadAccuracyReal(P%Accuracy%SourcekAccuracyBoost, 'SourcekAccuracyBoost')
    call ReadAccuracyReal(P%Accuracy%IntkAccuracyBoost, 'IntkAccuracyBoost')
    call ReadAccuracyReal(P%Accuracy%TransferkBoost, 'TransferkBoost')
    call ReadAccuracyReal(P%Accuracy%NonFlatIntAccuracyBoost, 'NonFlatIntAccuracyBoost')
    call ReadAccuracyReal(P%Accuracy%BessIntBoost, 'BessIntBoost')
    call ReadAccuracyReal(P%Accuracy%LensingBoost, 'LensingBoost')
    call ReadAccuracyReal(P%Accuracy%NonlinSourceBoost, 'NonlinSourceBoost')
    call ReadAccuracyReal(P%Accuracy%BesselBoost, 'BesselBoost')
    call ReadAccuracyReal(P%Accuracy%LimberBoost, 'LimberBoost')
    call ReadAccuracyReal(P%Accuracy%SourceLimberBoost, 'SourceLimberBoost')
    call ReadAccuracyReal(P%Accuracy%neutrino_q_boost, 'neutrino_q_boost')
    if (ErrMsg /= '') return

    P%NonLinear = Ini%Read_Int('do_nonlinear', NonLinear_none)

    P%Evolve_baryon_cs = Ini%Read_Logical('evolve_baryon_cs', .false.)
    P%Evolve_delta_xe = Ini%Read_Logical('evolve_delta_xe', .false.)
    P%Evolve_delta_Ts = Ini%Read_Logical('evolve_delta_ts', .false.)

    P%DoLensing = .false.
    P%Min_l = Ini%Read_int('l_min',2)
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini%Read_Int('l_max_scalar')
            P%Max_eta_k = Ini%Read_Double('k_eta_max_scalar', P%Max_l*2.5_dl)
            P%lens_output_margin = Ini%Read_Int('lens_output_margin', P%lens_output_margin)
            if (P%WantScalars) then
                P%DoLensing = Ini%Read_Logical('do_lensing', .false.)
                if (P%DoLensing) lensing_method = Ini%Read_Int('lensing_method', 1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) then
                    ErrMsg = 'Must generate vector modes on their own'
                    return
                end if
                i = Ini%Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    ErrMsg = 'vector_mode must be 0 (regular) or 1 (magnetic)'
                    return
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini%Read_Int('l_max_tensor')
            P%Max_eta_k_tensor = Ini%Read_Double('k_eta_max_tensor', Max(500._dl, P%Max_l_tensor * 2._dl))
        end if
    endif

    !  Read initial parameters.
    DarkEneryModel = UpperCase(Ini%Read_String_Default('dark_energy_model', 'fluid'))
    if (allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
    if (DarkEneryModel == 'FLUID') then
        allocate (TDarkEnergyFluid::P%DarkEnergy)
    else if (DarkEneryModel == 'PPF') then
        allocate (TDarkEnergyPPF::P%DarkEnergy)
    else if (DarkEneryModel == 'AXIONEFFECTIVEFLUID') then
        allocate (TAxionEffectiveFluid::P%DarkEnergy)
    else if (DarkEneryModel == 'EARLYQUINTESSENCE') then
        allocate (TEarlyQuintessence::P%DarkEnergy)
    else
        ErrMsg = 'Unknown dark energy model: '//trim(DarkEneryModel)
        return
    end if
    call P%DarkEnergy%ReadParams(Ini)

    RecombinationModel = UpperCase(Ini%Read_String_Default('recombination_model', 'Recfast'))
    if (RecombinationModel == 'COSMOREC') then
#ifdef COSMOREC
        deallocate(P%Recomb)
        allocate(TCosmoRec::P%Recomb)
#else
        ErrMsg = 'Compile with CosmoRec to use recombination_model=CosmoRec'
        return
#endif
    else if (RecombinationModel == 'HYREC') then
#ifdef HYREC
        deallocate(P%Recomb)
        allocate(THyRec::P%Recomb)
#else
        ErrMsg = 'Compile with HyRec to use recombination_model=HyRec'
        return
#endif
    else if (RecombinationModel /= 'RECFAST') then
        ErrMsg =  'Unknown recombination_model: '//trim(RecombinationModel)
        return
    end if

    call P%Recomb%ReadParams(Ini)
    if (global_error_flag /= 0) then
        ErrMsg = trim(global_error_message)
        return
    end if

    has_hubble = Ini%HasKey('hubble')
    has_thetastar = Ini%HasKey('thetastar')
    has_cosmomc_theta = Ini%HasKey('cosmomc_theta')
    num_theta_inputs = merge(1, 0, has_hubble) + merge(1, 0, has_thetastar) + merge(1, 0, has_cosmomc_theta)
    if (num_theta_inputs > 1) then
        ErrMsg = 'Can only set one of hubble, thetastar, cosmomc_theta'
        return
    else if (num_theta_inputs == 0) then
        ErrMsg = 'Must set one of hubble, thetastar, cosmomc_theta'
        return
    else if (has_hubble) then
        P%h0 = Ini%Read_Double('hubble')
    else if (has_thetastar) then
        theta = Ini%Read_Double('thetastar')
    else
        theta = Ini%Read_Double('cosmomc_theta')
    end if

    if (Ini%Read_Logical('use_physical', .true.)) then
        P%ombh2 = Ini%Read_Double('ombh2')
        P%omch2 = Ini%Read_Double('omch2')
        P%omnuh2 = Ini%Read_Double('omnuh2')
        P%omk = Ini%Read_Double('omk')
    else
        ErrMsg = 'use_physical = F no longer supported. Use ombh2, omch2, omnuh2, omk'
        return
    end if

    P%tcmb = Ini%Read_Double('temp_cmb', COBE_CMBTemp)
    P%yhe = Ini%Read_Double('helium_fraction', 0.24_dl)
    P%Num_Nu_massless = Ini%Read_Double('massless_neutrinos')

    P%Nu_mass_eigenstates = Ini%Read_Int('nu_mass_eigenstates', 1)
    if (P%Nu_mass_eigenstates > max_nu) then
        ErrMsg = 'too many mass eigenstates'
        return
    end if

    numstr = Ini%Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) then
        ErrMsg =  'massive_neutrinos should now be integer (or integer array)'
        return
    end if
    read(numstr,*, iostat=status) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    if (status/=0) then
        ErrMsg = 'Must give num_massive number of integer physical neutrinos for each eigenstate'
        return
    end if
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini%Read_Logical('share_delta_neff', .true.)
        numstr = Ini%Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='' .and. print_fortran_warnings) &
                write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') then
                ErrMsg = 'must give degeneracies for each eigenstate if share_delta_neff=F'
                return
            end if
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini%Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) then
                ErrMsg =  'must give nu_mass_fractions for the eigenstates'
                return
            end if
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if

    if (has_thetastar .or. has_cosmomc_theta) then
        if (.not. CAMB_SetH0ForTheta(P, theta, ErrMsg, cosmomc_approx=has_cosmomc_theta)) return
    end if

    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) .or. PK_WantTransfer) then
        P%Transfer%high_precision = Ini%Read_Logical('transfer_high_precision', .false.)
    else
        P%transfer%high_precision = .false.
    endif
    if (PK_WantTransfer) then
        P%Transfer%accurate_massive_neutrinos = Ini%Read_Logical('accurate_massive_neutrino_transfers',.false.)
    else
        P%Transfer%accurate_massive_neutrinos = .false.
    end if
    if (P%NonLinear/=NonLinear_none) then
        S = UpperCase(Ini%Read_String_Default('nonlinear_model', 'Halofit'))
        if (allocated(P%NonLinearModel)) deallocate(P%NonLinearModel)
        if (S == 'SPKNONLINEAR') then
            allocate(TSPkNonLinear::P%NonLinearModel)
        else if (S == 'HALOFIT') then
            allocate(THalofit::P%NonLinearModel)
        else
            ErrMsg = 'Unknown non-linear model: '//trim(S)
            return
        end if
        call P%NonLinearModel%ReadParams(Ini)
    end if

    if (PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax = Ini%Read_Double('transfer_kmax')*(P%h0 / 100._dl)
        P%transfer%k_per_logint = Ini%Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts = Ini%Read_Int('transfer_num_redshifts')

        if (P%Do21cm) P%transfer_21cm_cl = Ini%Read_Logical('transfer_21cm_cl',.false.)
        if (P%transfer_21cm_cl .and. P%transfer%kmax > 800 .and. print_fortran_warnings) then
            !Actually line widths are important at significantly larger scales too
            write (*,*) 'WARNING: kmax very large. '
            write(*,*) ' -- Neglected line width effects will dominate'
        end if

        call Ini%Read('transfer_interp_matterpower', transfer_interp_matterpower)
        call Ini%Read('transfer_power_var', transfer_power_var)
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) then
            ErrMsg = 'Too many redshifts, increase max_transfer_redshifts'
            return
        end if
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini%Read_Double_Array('transfer_redshift', i, 0._dl)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    call Ini%Read('Alens', P%Alens)

    call P%Reion%ReadParams(Ini)
    call P%InitPower%ReadParams(Ini)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini%Read_Int('initial_condition', initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            allocate(P%InitialCOnditionVector(initial_nummodes))
            numstr = Ini%Read_String('initial_vector', .true.)
            read (numstr,*) P%InitialConditionVector
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) P%use_cl_spline_template = .false.
    end if
    if (P%Scalar_initial_condition== initial_adiabatic) &
        call Ini%Read('use_cl_spline_template', P%use_cl_spline_template)

    P%WantDerivedParameters = Ini%Read_Logical('derived_parameters', .true.)

    !optional parameters controlling the computation

    call ReadAccuracyLogical(P%Accuracy%AccuratePolarization, 'AccuratePolarization', 'accurate_polarization')
    call ReadAccuracyLogical(P%Accuracy%AccurateReionization, 'AccurateReionization', 'accurate_reionization')
    call ReadAccuracyLogical(P%Accuracy%AccurateBB, 'AccurateBB', 'accurate_BB')
    if (ErrMsg /= '') return
    if (print_fortran_warnings .and. P%Accuracy%AccurateBB .and. P%WantCls .and. (P%Max_l < 3500 .or. &
        (P%NonLinear/=NonLinear_lens .and. P%NonLinear/=NonLinear_both) .or. P%Max_eta_k < 18000)) &
        write(*,*) 'WARNING: for accurate lensing BB you need high l_max_scalar, k_eta_max_scalar and non-linear lensing'

    !Mess here to fix typo with backwards compatibility
    if (Ini%HasKey('do_late_rad_trunction')) then
        P%DoLateRadTruncation = Ini%Read_Logical('do_late_rad_trunction', .true.)
        if (Ini%HasKey('do_late_rad_truncation')) error stop 'check do_late_rad_xxxx'
    else
        P%DoLateRadTruncation = Ini%Read_Logical('do_late_rad_truncation', .true.)
    end if
    P%MassiveNuMethod = Ini%Read_Int('massive_nu_approx', Nu_best)

    call ReadAccuracyReal(P%Accuracy%lSampleBoost, 'lSampleBoost', 'l_sample_boost')
    call Ini%Read('min_l_logl_sampling', P%min_l_logl_sampling)
    if (ErrMsg /= '') return

    CAMB_ReadParams = .true.

    contains

    subroutine ReadAccuracyReal(Value, Name, IniName)
    real(dl), intent(inout) :: Value
    character(LEN=*), intent(in) :: Name
    character(LEN=*), optional, intent(in) :: IniName

    if (present(IniName)) then
        if (Ini%HasKey(Name) .and. Ini%HasKey(IniName)) then
            ErrMsg = 'Cannot set both '//trim(Name)//' and '//trim(IniName)
            return
        end if
    end if

    if (Ini%HasKey(Name)) then
        call Ini%Read(Name, Value)
    else if (present(IniName)) then
        call Ini%Read(IniName, Value)
    end if

    end subroutine ReadAccuracyReal

    subroutine ReadAccuracyLogical(Value, Name, IniName)
    logical, intent(inout) :: Value
    character(LEN=*), intent(in) :: Name
    character(LEN=*), optional, intent(in) :: IniName

    if (present(IniName)) then
        if (Ini%HasKey(Name) .and. Ini%HasKey(IniName)) then
            ErrMsg = 'Cannot set both '//trim(Name)//' and '//trim(IniName)
            return
        end if
    end if

    if (Ini%HasKey(Name)) then
        call Ini%Read(Name, Value)
    else if (present(IniName)) then
        call Ini%Read(IniName, Value)
    end if

    end subroutine ReadAccuracyLogical

    end function CAMB_ReadParams

    logical function CAMB_RunFromIni(Ini, InputFile, ErrMsg)
    use IniObjects
    use Lensing
    use constants
    use Bispectrum
    use CAMBmain
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: InputFile
    character(LEN=*), intent(inout) :: ErrMsg
    type(CAMBparams) P
    character(len=:), allocatable :: outroot, VectorFileName, &
        ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName, ScalarCovFileName, &
        version_check, ArrayKey
    integer :: i
    character(len=Ini_max_string_len), allocatable :: TransferFileNames(:), &
        MatterPowerFileNames(:), TransferClFileNames(:)
    real(dl) :: output_factor
    logical PK_WantTransfer
    Type(CAMBdata) :: ActiveState

    call SetActiveState(ActiveState)
    CAMB_RunFromIni = .false.

    if (.not. CAMB_ReadParams(P, Ini, ErrMsg)) return

    outroot = Ini%Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    PK_WantTransfer = Ini%Read_Logical('get_transfer')
    if (PK_WantTransfer)  then
        call Ini%Read('transfer_interp_matterpower', transfer_interp_matterpower)
        call Ini%Read('transfer_power_var', transfer_power_var)
        allocate (TransferFileNames(P%Transfer%PK_num_redshifts))
        allocate (MatterPowerFileNames(P%Transfer%PK_num_redshifts))
        allocate (TransferClFileNames(P%Transfer%PK_num_redshifts))
        do i=1, P%transfer%PK_num_redshifts
            ArrayKey = Ini%Key_To_Arraykey('transfer_filename', i)
            if (i == 1) then
                TransferFileNames(i) = Ini%Read_String_Default(ArrayKey, 'transfer_out.dat')
            else
                TransferFileNames(i) = Ini%Read_String_Default(ArrayKey, trim(numcat('transfer_',i))//'.dat')
            end if

            ArrayKey = Ini%Key_To_Arraykey('transfer_matterpower', i)
            if (i == 1) then
                MatterPowerFilenames(i) = Ini%Read_String_Default(ArrayKey, 'matterpower.dat', .true.)
            else
                MatterPowerFilenames(i) = Ini%Read_String_Default(ArrayKey, trim(numcat('matterpower_',i))//'.dat', .true.)
            end if

            TransferFileNames(i) = outroot // TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') MatterPowerFilenames(i) = outroot // MatterPowerFilenames(i)

            if (P%Do21cm) then
                TransferClFileNames(i) = Ini%Read_String_Array('transfer_cl_filename',i)
                if (TransferClFileNames(i) == '') &
                    TransferClFileNames(i) =  trim(numcat('sharp_cl_',i))//'.dat'
            else
                TransferClFileNames(i) = ''
            end if

            if (TransferClFileNames(i)/= '') &
                TransferClFileNames(i) = outroot // TransferClFileNames(i)
        end do
    end if

    call Bispectrum_ReadParams(BispectrumParams, Ini, outroot)

    output_factor = Ini%Read_Double('CMB_outputscale', 1.d0)

    if (P%WantScalars) then
        ScalarFileName = Ini%Read_String_Default('scalar_output_file', 'scalCls.dat', .true.)
        if (ScalarFileName /= '') ScalarFileName = outroot // ScalarFileName
        LensedFileName = Ini%Read_String_Default('lensed_output_file', 'lensedCls.dat')
        LensedFileName = outroot // LensedFileName
        LensPotentialFileName = Ini%Read_String_Default('lens_potential_output_file', 'lenspotentialCls.dat')
        LensPotentialFileName = outroot // LensPotentialFileName
        ScalarCovFileName = Ini%Read_String_Default('scalar_covariance_output_file', &
            'scalarCovCls.dat', .false.)
        if (ScalarCovFileName /= '') then
            P%want_cl_2D_array = .true.
            ScalarCovFileName = outroot // ScalarCovFileName
        end if
    else
        ScalarFileName = ''
        LensedFileName = ''
        LensPotentialFileName = ''
        ScalarCovFileName = ''
    end if
    if (P%WantTensors) then
        TensorFileName = Ini%Read_String_Default('tensor_output_file', 'tensCls.dat')
        TensorFileName = outroot // TensorFileName
        if (P%WantScalars) then
            TotalFileName = Ini%Read_String_Default('total_output_file', 'totCls.dat', .true.)
            if (TotalFileName /= '') TotalFileName = outroot // TotalFileName
            LensedTotFileName = Ini%Read_String_Default('lensed_total_output_file', 'lensedtotCls.dat')
            LensedTotFileName = outroot // LensedTotFileName
        else
            TotalFileName = ''
            LensedTotFileName = ''
        end if
    else
        TensorFileName = ''
        TotalFileName = ''
        LensedTotFileName = ''
    end if
    if (P%WantVectors) then
        VectorFileName = Ini%Read_String_Default('vector_output_file', 'vecCls.dat')
        VectorFileName = outroot // VectorFileName
    else
        VectorFileName = ''
    end if

    version_check = Ini%Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call Ini%ReadValues%Add('version_check', version)
    else if (version_check /= version .and. print_fortran_warnings) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if

    call Ini%Read('output_file_headers', output_file_headers)

    if (do_bispectrum) then
        P%Accuracy%lSampleBoost   = 50
    end if
    if (outroot /= '') then
        if (InputFile /= trim(outroot) // 'params.ini') then
            call Ini%SaveReadValues(trim(outroot) // 'params.ini')
        else
            write(*,*) 'Output _params.ini not created as would overwrite input'
        end if
    end if

    if (.not. P%Validate()) then
        ErrMsg = 'Invalid parameter value'
        return
    end if

    if (global_error_flag==0) call CAMB_GetResults(State,P)
    if (global_error_flag/=0) then
        ErrMsg =  trim(global_error_message)
        return
    endif

    if (PK_WantTransfer) then
        call Transfer_SaveToFiles(State%MT,State, TransferFileNames)
        call Transfer_SaveMatterPower(State%MT,State,MatterPowerFileNames)
        call Transfer_output_sig8(State%MT, State)
        if (P%do21cm .and. P%transfer_21cm_cl) call Transfer_Get21cmCls(State%MT, State,TransferClFileNames)
    end if

    if (P%WantCls) then
        call State%CLData%output_cl_files(State,ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
            LensedFileName, LensedTotFilename, output_factor)

        call State%CLData%output_lens_pot_files(State%CP,LensPotentialFileName, output_factor)

        if (P%WantVectors) then
            call State%CLData%output_veccl_files(State%CP,VectorFileName, output_factor)
        end if

    end if

    CAMB_RunFromIni = .true.

    end function CAMB_RunFromIni

    logical function CAMB_RunIniFile(InputFile, InpLen)
    integer, intent(in) :: InpLen
    character(LEN=InpLen), intent(in) :: InputFile
    character(LEN=len(global_error_message)) :: ErrMsg
    Type(TIniFile) :: Ini
    logical bad

    !Same as CAMB_CommandLineRun but does not read variables that change global state
    !Intended for use from python

    call Ini%Open(InputFile, bad, .false.)
    Ini%Fail_on_not_found = .false.
    ErrMsg = ''
    CAMB_RunIniFile = CAMB_RunFromIni(Ini, InputFile, ErrMsg)
    call Ini%Close()
    if (ErrMsg/='') call GlobalError(ErrMsg,error_ini)

    end function CAMB_RunIniFile

    subroutine CAMB_CommandLineValidate(InputFile)
    !Error stop if any problem
    character(LEN=*), intent(in) :: InputFile
    Type(TIniFile) :: Ini
    logical bad
    Type(CAMBParams) :: P
    character(LEN=1024) :: ErrMsg

    call Ini%Open(InputFile, bad, .false.)
    if (bad) then
        error stop 'File does not exist'
    end if
    if (.not. CAMB_ReadParams(P, Ini, ErrMsg)) then
        write(*,*) trim(ErrMsg)
        error stop
    end if

    call Ini%Close()

    end subroutine CAMB_CommandLineValidate

    subroutine CAMB_CommandLineRun(InputFile)
    character(LEN=*), intent(in) :: InputFile
    Type(TIniFile) :: Ini
    character(LEN=1024) :: ErrMsg
    logical bad

    call Ini%Open(InputFile, bad, .false.)
    if (bad) then
        write(*,*) 'File not found: '//trim(InputFile)
        error stop
    end if

    highL_unlensed_cl_template = Ini%Read_String_Default( &
        'highL_unlensed_cl_template', highL_unlensed_cl_template)
    call Ini%Read('number_of_threads', ThreadNum)
    call Ini%Read('AccuracyTarget', AccuracyTarget)
    call Ini%Read('DebugParam', DebugParam)
    if (Ini%HasKey('enable_olver_source_integration')) &
        call Ini%Read('enable_olver_source_integration', enable_olver_source_integration)
    call Ini%Read('feedback_level', FeedbackLevel)
    call Ini%Read('print_fortran_warnings', print_fortran_warnings)
    if (Ini%HasKey('DebugMsgs')) call Ini%Read('DebugMsgs', DebugMsgs)

    Ini%Fail_on_not_found = .false.
    if (.not. CAMB_RunFromIni(Ini, InputFile, ErrMsg)) then
        write(*,*) trim(ErrMsg)
        error stop 'Invalid parameter'
    end if
    call Ini%Close()

    end subroutine CAMB_CommandLineRun

    subroutine CAMB_GetVersion(ver)
    character(LEN=*) :: ver
    ver = version
    end subroutine CAMB_GetVersion

    subroutine CAMB_FreeGlobalMemory()
    use SpherBessels
    use CAMBmain
    call Bessels_Free()
    !flat sky lensing also allocates bessels, but not normally used.
    call CAMBMain_Free()
    if (allocated(highL_CL_template)) deallocate(highL_CL_template)
    if (allocated(nu_tau_notmassless)) deallocate(nu_tau_notmassless)

    end subroutine CAMB_FreeGlobalMemory

    end module CAMB
