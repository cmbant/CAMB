    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     See readme.html for documentation. This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use IniObjects
    use CAMB
    use Lensing
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use NonLinear
#ifdef NAGF95
    use F90_UNIX
#endif
    implicit none

    type(CAMBparams) P

    character(len=:), allocatable :: numstr, outroot, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName, ScalarCovFileName, version_check
    integer :: i
    ! max_transfer_redshifts
    character(len=Ini_max_string_len), allocatable :: TransferFileNames(:), MatterPowerFileNames(:)
    real(dl) :: output_factor, nmassive

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif
    type(TIniFile) :: Ini
    logical bad

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') stop 'No parameter input file'

    call Ini%Open(InputFile, bad, .false.)
    if (bad) stop 'Error opening parameter file'

    Ini%Fail_on_not_found = .false.

    outroot = Ini%Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini%Read_String_Default('highL_unlensed_cl_template', highL_unlensed_cl_template)

    call CAMB_SetDefParams(P)

    P%WantScalars = Ini%Read_Logical('get_scalar_cls')
    P%WantVectors = Ini%Read_Logical('get_vector_cls', .false.)
    P%WantTensors = Ini%Read_Logical('get_tensor_cls', .false.)

    P%OutputNormalization=outNone
    output_factor = Ini%Read_Double('CMB_outputscale', 1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    P%PK_WantTransfer = Ini%Read_Logical('get_transfer')

    call Ini%Read('accuracy_boost', AccuracyBoost)
    call Ini%Read('l_accuracy_boost', lAccuracyBoost)
    call Ini%Read('high_accuracy_default', HighAccuracyDefault)

    P%NonLinear = Ini%Read_Int('do_nonlinear', NonLinear_none)

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
                if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
                i = Ini%Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini%Read_Int('l_max_tensor')
            P%Max_eta_k_tensor = Ini%Read_Double('k_eta_max_tensor', Max(500._dl, P%Max_l_tensor * 2._dl))
        end if
    endif

    !  Read initial parameters.

	allocate (TDarkEnergy::DarkEnergy)
    call DarkEnergy%ReadParams(Ini)
    call DarkEnergy%Init_Background()

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
    if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

    numstr = Ini%Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive - nint(nmassive))>1e-6) stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini%Read_Logical('share_delta_neff', .true.)
        numstr = Ini%Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini%Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if

    !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
    !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
    !in the P%WantTransfer loop.
    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) &
        .or. P%PK_WantTransfer) then
    P%Transfer%high_precision = Ini%Read_Logical('transfer_high_precision', .false.)
    else
        P%transfer%high_precision = .false.
    endif
    if (P%NonLinear/=NonLinear_none) call NonLinear_ReadParams(Ini)

    if (P%PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax = Ini%Read_Double('transfer_kmax')
        P%transfer%k_per_logint = Ini%Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts = Ini%Read_Int('transfer_num_redshifts')

        call Ini%Read('transfer_interp_matterpower ', transfer_interp_matterpower)
        call Ini%Read('transfer_power_var', transfer_power_var)
        !        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        allocate (TransferFileNames(P%Transfer%PK_num_redshifts))
        allocate (MatterPowerFileNames(P%Transfer%PK_num_redshifts))
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
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
        P%WantTransfer  = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if

    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !JD 08/13 end changes

    P%transfer%kmax = P%transfer%kmax*(P%h0 / 100._dl)

    Ini%Fail_on_not_found = .false.

    call Ini%Read('DebugParam', DebugParam)
    call Ini%Read('Alens', Alens)

    call Reionization_ReadParams(P%Reion, Ini)
    call InitialPower_ReadParams(P%InitPower, Ini, P%WantTensors)
    call Recombination_ReadParams(P%Recomb, Ini)
    if (Ini%HasKey('recombination')) then
        i = Ini%Read_Int('recombination', 1)
        if (i/=1) stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, Ini, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini%Read_Int('initial_condition', initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            P%InitialConditionVector=0
            numstr = Ini%Read_String('initial_vector', .true.)
            read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
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
            has_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot, ScalarCovFileName)
        end if
    end if
    if (P%WantTensors) then
        TensorFileName = trim(outroot) // Ini%Read_String('tensor_output_file')
        if (P%WantScalars) then
            TotalFileName = trim(outroot) // Ini%Read_String('total_output_file')
            LensedTotFileName = Ini%Read_String('lensed_total_output_file')
            if (LensedTotFileName /= '') LensedTotFileName = trim(outroot) // trim(LensedTotFileName)
        end if
    end if
    if (P%WantVectors) then
        VectorFileName = trim(outroot) // Ini%Read_String('vector_output_file')
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

    P%AccuratePolarization = Ini%Read_Logical('accurate_polarization', .true.)
    P%AccurateReionization = Ini%Read_Logical('accurate_reionization', .false.)
    P%AccurateBB = Ini%Read_Logical('accurate_BB', .false.)
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
        DoLateRadTruncation = Ini%Read_Logical('do_late_rad_trunction', .true.)
        if (Ini%HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
    else
        DoLateRadTruncation = Ini%Read_Logical('do_late_rad_truncation', .true.)
    end if

    if (HighAccuracyDefault) then
        DoTensorNeutrinos = .true.
    else
        call Ini%Read('do_tensor_neutrinos', DoTensorNeutrinos )
    end if
    call Ini%Read('feedback_level', FeedbackLevel)

    P%MassiveNuMethod = Ini%Read_Int('massive_nu_approx', Nu_best)

    call Ini%Read('number_of_threads', ThreadNum)
    call Ini%Read('use_spline_template', use_spline_template)

    if (do_bispectrum) then
        lSampleBoost   = 50
    else
        call Ini%Read('l_sample_boost', lSampleBoost)
    end if
    if (outroot /= '') then
        if (InputFile /= trim(outroot) // 'params.ini') then
            call Ini%SaveReadValues(trim(outroot) // 'params.ini')
        else
            write(*,*) 'Output _params.ini not created as would overwrite input'
        end if
    end if

    call Ini%Close()

    if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

#ifdef RUNIDLE
    call SetIdle
#endif

    if (global_error_flag==0) call CAMB_GetResults(P)
    if (global_error_flag/=0) then
        write (*,*) 'Error result '//trim(global_error_message)
        stop
    endif

    if (P%PK_WantTransfer) then
        call Transfer_SaveToFiles(MT,TransferFileNames)
        call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
        call Transfer_output_sig8(MT)
    end if

    if (P%WantCls) then
        call output_cl_files(ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
            LensedFileName, LensedTotFilename, output_factor)

        call output_lens_pot_files(LensPotentialFileName, output_factor)

        if (P%WantVectors) then
            call output_veccl_files(VectorFileName, output_factor)
        end if

#ifdef WRITE_FITS
        if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
    end if

    if (allocated(MatterPowerFileNames)) deallocate (MatterPowerFileNames, TransferFileNames)
    call CAMB_cleanup
    stop

100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'
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
