    !Interface module for CAMB. Call CAMB_GetResults to do the work.

    module CAMB
    use Precision
    use results
    use GaugeInterface
    use InitialPower
    use Reionization
    use Recombination
    use lensing
    use DarkEnergyFluid
    implicit none
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
        if (State%CP%DoLensing .and. global_error_flag==0) call lens_Cls(Cdata)
        if (global_error_flag/=0) return
    end if

    if (CData%CP%WantTransfer) call Transfer_Get_sigmas(Cdata, Cdata%MT)

    end subroutine CAMB_TransfersToPowers

    !Call this routine with a set of parameters to generate the results you want.
    subroutine CAMB_GetResults(OutData, Params, error, onlytransfer, onlytimesources)
    use CAMBmain
    use lensing
    use Bispectrum
    type(CAMBdata)  :: OutData
    type(CAMBparams) :: Params
    integer, optional :: error !Zero if OK
    logical, optional :: onlytransfer, onlytimesources
    type(CAMBparams) P
    logical :: call_again

    global_error_flag = 0
    call_again = .false.
    call OutData%Free()
    call SetActiveState(OutData)
    OutData%HasScalarTimeSources= DefaultFalse(onlytimesources)
    OutData%OnlyTransfer = DefaultFalse(onlytransfer) .or. OutData%HasScalarTimeSources

    !Vector and tensors first, so at end time steps in state are for scalars
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
            return
        end if
        call_again = .true.
    end if

    if (Params%WantCls .and. Params%WantScalars) then
        P = Params
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        P%WantTensors = .false.
        P%WantVectors = .false.
        if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. &
            (P%DoLensing .or. State%num_redshiftwindows > 0)) then
            P%WantTransfer  = .true.
        end if
        call OutData%SetParams(P)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
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
            return
        end if
    end if
    OutData%CP%WantCls = Params%WantCls
    OutData%CP%WantScalars = Params%WantScalars
    OutData%CP%WantTensors = Params%WantTensors
    OutData%CP%WantVectors = Params%WantVectors
    OutData%CP%WantTransfer = Params%WantTransfer
    OutData%CP%Accuracy = Params%Accuracy
    OutData%CP%Reion%Reionization = Params%Reion%Reionization
    OutData%CP%Transfer%high_precision = Params%Transfer%high_precision
    OutData%CP%WantDerivedParameters = Params%WantDerivedParameters

    if (.not. OutData%OnlyTransfer .and. Params%WantCls .and. Params%WantScalars) then
        if (Params%DoLensing .and. global_error_flag==0) then
            call lens_Cls(OutData)
        end if

        if (do_bispectrum .and. global_error_flag==0) &
            call GetBispectrum(OutData,OutData%CLData%CTransScal)
    end if
    if (global_error_flag/=0 .and. present(error)) error =global_error_flag

    end subroutine CAMB_GetResults


    !Return real (NOT double precision) arrays of the computed CMB  Cls
    !Output is l(l+1)C_l/2pi
    !If GC_Conventions = .false. use E-B conventions (as the rest of CAMB does)
    !Used by WriteFits only
    subroutine CAMB_GetCls(State, Cls, lmax,GC_conventions)
    Type(CAMBdata) :: State
    integer, intent(IN) :: lmax
    logical, intent(IN) :: GC_conventions
    real, intent(OUT) :: Cls(2:lmax,1:4)
    integer l

    Cls = 0
    do l=2, lmax
        if (State%CP%WantScalars .and. l<= State%CP%Max_l) then
            if (State%CP%DoLensing) then
                if (l<=State%CLData%lmax_lensed) &
                    Cls(l,1:4) = State%CLData%Cl_lensed(l, CT_Temp:CT_Cross)
            else
                Cls(l,1:2) = State%CLData%Cl_scalar(l, C_Temp:C_E)
                Cls(l,4) = State%CLData%Cl_scalar(l, C_Cross)
            endif
        end if
        if (State%CP%WantTensors .and. l <= State%CP%Max_l_tensor) then
            Cls(l,1:4) = Cls(l,1:4) + State%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
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
    Type(CAMBdata) :: State

    call  State%SetParams(P, error, .false., .false., .true.)

    if (error/=0) then
        CAMB_GetAge = -1
    else
        CAMB_GetAge = State%DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
    end if

    end function CAMB_GetAge

    subroutine CAMB_SetDefParams(P)
    use NonLinear
    use Recombination
    type(CAMBparams), intent(inout), target :: P
    type(CAMBparams) :: emptyP

    P= emptyP !Set default values set in type definitions
    P%Nu_mass_numbers=0
    P%Nu_mass_degeneracies=0
    P%Nu_mass_fractions=0

    allocate(THalofit::P%NonLinearModel)
    allocate(TDarkEnergyFluid::P%DarkEnergy)
    allocate(TInitialPowerLaw::P%InitPower)
    allocate(TRecfast::P%Recomb)
    allocate(TTanhReionization::P%Reion)

    end subroutine CAMB_SetDefParams

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
    integer i, status
    real(dl) nmassive
    character(LEN=*), intent(inout) :: ErrMsg
    character(LEN=:), allocatable :: NumStr, S, DarkEneryModel, RecombinationModel
    logical :: DoCounts

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
        call Ini%Read('Kmax_Boost', P%Accuracy%KmaxBoost)
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
                if (RedWin%sigma < 0.003) then
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

    call Ini%Read('accuracy_boost', P%Accuracy%AccuracyBoost)
    call Ini%Read('l_accuracy_boost', P%Accuracy%lAccuracyBoost)

    P%NonLinear = Ini%Read_Int('do_nonlinear', NonLinear_none)

    P%Evolve_baryon_cs = Ini%Read_Logical('evolve_baryon_cs', .false.)
    P%Evolve_delta_xe = Ini%Read_Logical('evolve_delta_xe', .false.)
    P%Evolve_delta_Ts = Ini%Read_Logical('evolve_delta_ts', .false.)

    P%DoLensing = .false.
    P%Min_l = Ini%Read_int('l_min',2)
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini%Read_Int('l_max_scalar')
            P%Max_eta_k = Ini%Read_Double('k_eta_max_scalar', P%Max_l*2._dl)
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

    P%h0 = Ini%Read_Double('hubble')

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
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
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
    if (P%NonLinear/=NonLinear_none) call P%NonLinearModel%ReadParams(Ini)

    if (PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax = Ini%Read_Double('transfer_kmax')*(P%h0 / 100._dl)
        P%transfer%k_per_logint = Ini%Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts = Ini%Read_Int('transfer_num_redshifts')

        if (P%Do21cm) P%transfer_21cm_cl = Ini%Read_Logical('transfer_21cm_cl',.false.)
        if (P%transfer_21cm_cl .and. P%transfer%kmax > 800) then
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

    P%Accuracy%AccuratePolarization = Ini%Read_Logical('accurate_polarization', .true.)
    P%Accuracy%AccurateReionization = Ini%Read_Logical('accurate_reionization', .true.)
    P%Accuracy%AccurateBB = Ini%Read_Logical('accurate_BB', .false.)
    if (P%Accuracy%AccurateBB .and. P%WantCls .and. (P%Max_l < 3500 .or. &
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

    call Ini%Read('l_sample_boost', P%Accuracy%lSampleBoost)

    CAMB_ReadParams = .true.

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
        version_check
    integer :: i
    character(len=Ini_max_string_len), allocatable :: TransferFileNames(:), &
        MatterPowerFileNames(:), TransferClFileNames(:)
    real(dl) :: output_factor
#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif
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
            transferFileNames(i)     = Ini%Read_String_Array('transfer_filename', i)
            MatterPowerFilenames(i)  = Ini%Read_String_Array('transfer_matterpower', i)
            if (TransferFileNames(i) == '') then
                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
            end if
            if (MatterPowerFilenames(i) == '') then
                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
            end if
            if (TransferFileNames(i)/= '') &
                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') &
                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)

            if (P%Do21cm) then
                TransferClFileNames(i) = Ini%Read_String_Array('transfer_cl_filename',i)
                if (TransferClFileNames(i) == '') &
                    TransferClFileNames(i) =  trim(numcat('sharp_cl_',i))//'.dat'
            else
                TransferClFileNames(i) = ''
            end if

            if (TransferClFileNames(i)/= '') &
                TransferClFileNames(i) = trim(outroot)//TransferClFileNames(i)
        end do
    end if

    call Bispectrum_ReadParams(BispectrumParams, Ini, outroot)

    output_factor = Ini%Read_Double('CMB_outputscale', 1.d0)

    if (P%WantScalars) then
        ScalarFileName = trim(outroot) // Ini%Read_String('scalar_output_file')
        LensedFileName = trim(outroot) // Ini%Read_String('lensed_output_file')
        LensPotentialFileName = Ini%Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName = Ini%Read_String_Default('scalar_covariance_output_file', &
            'scalarCovCls.dat', .false.)
        if (ScalarCovFileName /= '') then
            P%want_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot, ScalarCovFileName)
        end if
    else
        ScalarFileName = ''
        LensedFileName = ''
        LensPotentialFileName = ''
        ScalarCovFileName = ''
    end if
    if (P%WantTensors) then
        TensorFileName = trim(outroot) // Ini%Read_String('tensor_output_file')
        if (P%WantScalars) then
            TotalFileName = trim(outroot) // Ini%Read_String('total_output_file')
            LensedTotFileName = Ini%Read_String('lensed_total_output_file')
            if (LensedTotFileName /= '') LensedTotFileName = trim(outroot) // trim(LensedTotFileName)
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
        VectorFileName = trim(outroot) // Ini%Read_String('vector_output_file')
    else
        VectorFileName = ''
    end if

#ifdef WRITE_FITS
    if (P%WantCls) then
        FITSfilename = trim(outroot) // Ini%Read_String('FITS_filename', .true.)
        if (FITSfilename /= '') then
            inquire(file=FITSfilename, exist=bad)
            if (bad) then
                open(unit=18, file=FITSfilename, status='old')
                close(18, status='delete')
            end if
        end if
    end if
#endif

    version_check = Ini%Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call Ini%ReadValues%Add('version_check', version)
    else if (version_check /= version) then
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

#ifdef WRITE_FITS
        if (FITSfilename /= '') call WriteFitsCls(State, FITSfilename, CP%Max_l)
#endif
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
    call Ini%Read('DebugParam', DebugParam)
    call Ini%Read('feedback_level', FeedbackLevel)
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

    end module CAMB
