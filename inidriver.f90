    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     See readme.html for documentation. This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use IniObjects
    use CAMB
    use Lensing
    use constants
    use Bispectrum
    use CAMBmain
    use NonLinear
    use DarkEnergyFluid
    use DarkEnergyPPF
    use CAMBsettings
    implicit none
    type(CAMBparams) P

    character(len=:), allocatable :: numstr, outroot, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName, ScalarCovFileName, &
        version_check, DarkEneryModel, S

    integer :: i
    ! max_transfer_redshifts
    character(len=Ini_max_string_len), allocatable :: TransferFileNames(:), &
        MatterPowerFileNames(:), TransferClFileNames(:)
    real(dl) :: output_factor, nmassive

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif
    type(TIniFile) :: Ini
    logical bad
    integer status
    logical :: DoCounts = .false.
    Type(CAMBstate) :: ActiveState
    logical PK_WantTransfer

    call SetActiveState(ActiveState)
    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call Ini%Open(InputFile, bad, .false.)
    if (bad) error stop 'Error opening parameter file'

    Ini%Fail_on_not_found = .false.

    outroot = Ini%Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini%Read_String_Default('highL_unlensed_cl_template', highL_unlensed_cl_template)

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
    if (P%SourceTerms%limber_windows) call Ini%Read('limber_phiphi', limber_phiphi)
    if (num_redshiftwindows > 0) then
        P%SourceTerms%do_counts_lensing = Ini%Read_Logical('DoRedshiftLensing', .false.)
        call Ini%Read('Kmax_Boost', P%Accuracy%Kmax_Boost)
    end if
    P%Do21cm = Ini%Read_Logical('Do21cm', .false.)
    num_extra_redshiftwindows = 0
    do i=1, num_redshiftwindows
        associate (RedWin => Redshift_w(i))
            call InitRedshiftWindow(RedWin)
            RedWin%Redshift = Ini%Read_Double_Array('redshift', i)
            S = Ini%Read_String_Array('redshift_kind', i)
            if (S == '21cm') then
                RedWin%kind = window_21cm
            elseif (S == 'counts') then
                RedWin%kind = window_counts
            elseif (S == 'lensing') then
                RedWin%kind = window_lensing
            else
                write (*,*) i, 'Error: unknown type of window '//trim(S)
                error stop
            end if
            RedWin%a = 1 / (1 + RedWin%Redshift)
            if (RedWin%kind /= window_21cm) then
                RedWin%sigma = Ini%Read_Double_Array('redshift_sigma', i)
                RedWin%sigma_z = RedWin%sigma
            else
                P%Do21cm = .true.
                RedWin%sigma = Ini%Read_Double_Array('redshift_sigma_Mhz', i)
                if (RedWin%sigma < 0.003) then
                    write(*,*) 'WARNING:Window very narrow.'
                    write(*,*) ' --> use transfer functions and transfer_21cm_cl =T ?'
                end if
                !with 21cm widths are in Mhz, make dimensionless scale factor
                RedWin%sigma = RedWin%sigma / (f_21cm / 1e6)
                RedWin%sigma_z = RedWin%sigma * (1 + RedWin%RedShift) ** 2
                write(*,*) i,'delta_z = ', RedWin%sigma_z
            end if
            if (RedWin%kind == window_counts) then
                DoCounts = .true.
                RedWin%bias = Ini%Read_Double_Array('redshift_bias', i)
                RedWin%dlog10Ndm = Ini%Read_Double_Array('redshift_dlog10Ndm', i ,0.d0)
                if (P%SourceTerms%do_counts_lensing) then
                    num_extra_redshiftwindows = num_extra_redshiftwindows + 1
                    RedWin%mag_index = num_extra_redshiftwindows
                end if
            end if
        end associate
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
    output_factor = Ini%Read_Double('CMB_outputscale', 1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    PK_WantTransfer = Ini%Read_Logical('get_transfer')

    call Ini%Read('accuracy_boost', P%Accuracy%AccuracyBoost)
    call Ini%Read('l_accuracy_boost', P%Accuracy%lAccuracyBoost)

    P%NonLinear = Ini%Read_Int('do_nonlinear', NonLinear_none)

    P%Evolve_baryon_cs = Ini%Read_Logical('evolve_baryon_cs', .false.)
    P%Evolve_delta_xe = Ini%Read_Logical('evolve_delta_xe', .false.)
    P%Evolve_delta_Ts = Ini%Read_Logical('evolve_delta_ts', .false.)

    P%DoLensing = .false.
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini%Read_Int('l_max_scalar')
            P%Max_eta_k = Ini%Read_Double('k_eta_max_scalar', P%Max_l*2._dl)
            if (P%WantScalars) then
                P%DoLensing = Ini%Read_Logical('do_lensing', .false.)
                if (P%DoLensing) lensing_method = Ini%Read_Int('lensing_method', 1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) error stop 'Must generate vector modes on their own'
                i = Ini%Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    error stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini%Read_Int('l_max_tensor')
            P%Max_eta_k_tensor = Ini%Read_Double('k_eta_max_tensor', Max(500._dl, P%Max_l_tensor * 2._dl))
        end if
    endif

    !  Read initial parameters.
    DarkEneryModel = Ini%Read_String_Default('dark_energy_model', 'fluid')
    if (allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
    if (DarkEneryModel == 'fluid') then
        allocate (TDarkEnergyFluid::P%DarkEnergy)
    else if (DarkEneryModel == 'PPF') then
        allocate (TDarkEnergyPPF::P%DarkEnergy)
    else if (DarkEneryModel == 'AxionEffectiveFluid') then
        allocate (TAxionEffectiveFluid::P%DarkEnergy)
    else
        error stop 'Unknown dark energy model'
    end if
    call P%DarkEnergy%ReadParams(Ini)

    P%h0 = Ini%Read_Double('hubble')

    if (Ini%Read_Logical('use_physical', .false.)) then
        P%omegab = Ini%Read_Double('ombh2') / (P%H0 / 100) ** 2
        P%omegac = Ini%Read_Double('omch2') / (P%H0 / 100) ** 2
        P%omegan = Ini%Read_Double('omnuh2') / (P%H0 / 100) ** 2
        P%omegav = 1- Ini%Read_Double('omk') - P%omegab - P%omegac - P%omegan
    else
        P%omegab = Ini%Read_Double('omega_baryon')
        P%omegac = Ini%Read_Double('omega_cdm')
        P%omegav = Ini%Read_Double('omega_lambda')
        P%omegan = Ini%Read_Double('omega_neutrino')
    end if

    P%tcmb = Ini%Read_Double('temp_cmb', COBE_CMBTemp)
    P%yhe = Ini%Read_Double('helium_fraction', 0.24_dl)
    P%Num_Nu_massless = Ini%Read_Double('massless_neutrinos')

    P%Nu_mass_eigenstates = Ini%Read_Int('nu_mass_eigenstates', 1)
    if (P%Nu_mass_eigenstates > max_nu) error stop 'too many mass eigenstates'

    numstr = Ini%Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) error stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, iostat=status) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    if (status/=0) error stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'

    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini%Read_Logical('share_delta_neff', .true.)
        numstr = Ini%Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') error stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini%Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) error stop 'must give nu_mass_fractions for the eigenstates'
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
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        allocate (TransferFileNames(P%Transfer%PK_num_redshifts))
        allocate (MatterPowerFileNames(P%Transfer%PK_num_redshifts))
        allocate (TransferClFileNames(P%Transfer%PK_num_redshifts))
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini%Read_Double_Array('transfer_redshift', i, 0._dl)
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
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if
    
    Ini%Fail_on_not_found = .false.

    call Ini%Read('DebugParam', DebugParam)
    call Ini%Read('Alens', P%Alens)

    call P%Reion%ReadParams(Ini)
    call P%InitPower%ReadParams(Ini, P%WantTensors)
    call P%Recomb%ReadParams(Ini)
    if (Ini%HasKey('recombination')) then
        i = Ini%Read_Int('recombination', 1)
        if (i/=1) error stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, Ini, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini%Read_Int('initial_condition', initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            allocate(P%InitialCOnditionVector(initial_nummodes))
            numstr = Ini%Read_String('initial_vector', .true.)
            read (numstr,*) P%InitialConditionVector
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
    end if

    if (P%WantScalars) then
        ScalarFileName = trim(outroot) // Ini%Read_String('scalar_output_file')
        LensedFileName = trim(outroot) // Ini%Read_String('lensed_output_file')
        LensPotentialFileName = Ini%Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName = Ini%Read_String_Default('scalar_covariance_output_file', &
            'scalCovCls.dat', .false.)
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

    Ini%Fail_on_not_found = .false.

    !optional parameters controlling the computation

    P%Accuracy%AccuratePolarization = Ini%Read_Logical('accurate_polarization', .true.)
    P%Accuracy%AccurateReionization = Ini%Read_Logical('accurate_reionization', .true.)
    P%Accuracy%AccurateBB = Ini%Read_Logical('accurate_BB', .false.)
    if (P%Accuracy%AccurateBB .and. P%WantCls .and. (P%Max_l < 3500 .or. &
        (P%NonLinear/=NonLinear_lens .and. P%NonLinear/=NonLinear_both) .or. P%Max_eta_k < 18000)) &
        write(*,*) 'WARNING: for accurate lensing BB you need high l_max_scalar, k_eta_max_scalar and non-linear lensing'

    P%DerivedParameters = Ini%Read_Logical('derived_parameters', .true.)

    version_check = Ini%Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call Ini%ReadValues%Add('version_check', version)
    else if (version_check /= version) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if
    !Mess here to fix typo with backwards compatibility
    if (Ini%HasKey('do_late_rad_trunction')) then
        P%DoLateRadTruncation = Ini%Read_Logical('do_late_rad_trunction', .true.)
        if (Ini%HasKey('do_late_rad_truncation')) error stop 'check do_late_rad_xxxx'
    else
        P%DoLateRadTruncation = Ini%Read_Logical('do_late_rad_truncation', .true.)
    end if

    call Ini%Read('feedback_level', FeedbackLevel)

    call Ini%Read('output_file_headers', output_file_headers)

    P%MassiveNuMethod = Ini%Read_Int('massive_nu_approx', Nu_best)

    call Ini%Read('number_of_threads', ThreadNum)
    call Ini%Read('use_spline_template', use_spline_template)

    if (do_bispectrum) then
        P%Accuracy%lSampleBoost   = 50
    else
        call Ini%Read('l_sample_boost', P%Accuracy%lSampleBoost)
    end if
    if (outroot /= '') then
        if (InputFile /= trim(outroot) // 'params.ini') then
            call Ini%SaveReadValues(trim(outroot) // 'params.ini')
        else
            write(*,*) 'Output _params.ini not created as would overwrite input'
        end if
    end if

    call Ini%Close()

    if (.not. P%Validate()) error stop 'Stopped due to parameter error'

#ifdef RUNIDLE
    call SetIdle
#endif

    if (global_error_flag==0) call CAMB_GetResults(State,P)
    if (global_error_flag/=0) then
        write (*,*) 'Error result '//trim(global_error_message)
        error stop
    endif

    if (PK_WantTransfer) then
        call Transfer_SaveToFiles(State%MT,TransferFileNames)
        call Transfer_SaveMatterPower(State%MT,MatterPowerFileNames)
        call Transfer_output_sig8(State%MT)
        if (P%do21cm .and. P%transfer_21cm_cl) call Transfer_Get21cmCls(State%MT, TransferClFileNames)
    end if

    if (P%WantCls) then
        call State%CLData%output_cl_files(ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
            LensedFileName, LensedTotFilename, output_factor)

        call State%CLData%output_lens_pot_files(LensPotentialFileName, output_factor)

        if (P%WantVectors) then
            call State%CLData%output_veccl_files(VectorFileName, output_factor)
        end if

#ifdef WRITE_FITS
        if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
    end if

    if (allocated(MatterPowerFileNames)) deallocate (MatterPowerFileNames, &
        TransferFileNames, TransferClFileNames)

    end program driver


#ifdef RUNIDLE
    !If in Windows and want to run with low priorty so can multitask
    subroutine SetIdle
    USE DFWIN
    Integer dwPriority
    Integer CheckPriority

    dwPriority = 64 ! idle priority
    CheckPriority = SetPriorityClass(GetCurrentProcess(), dwPriority)

    end subroutine SetIdle
#endif
