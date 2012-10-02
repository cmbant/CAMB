!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output. 

    program driver
        use IniFile
        use CAMB
        use LambdaGeneral
        use Lensing
        use AMLUtils
        use Transfer
        use constants
        use Bispectrum
#ifdef NAGF95 
        use F90_UNIX
#endif
        implicit none
      
        Type(CAMBparams) P
        
        character(LEN=Ini_max_string_len) numstr, VectorFileName, &
            InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
            LensedTotFileName, LensPotentialFileName
        integer i
        character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
               MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
        real(dl) output_factor, Age, nmassive

#ifdef WRITE_FITS
       character(LEN=Ini_max_string_len) FITSfilename
#endif

        logical bad

        InputFile = GetParam(1)
        if (InputFile == '') stop 'No parameter input file'

        call Ini_Open(InputFile, 1, bad, .false.)
        if (bad) stop 'Error opening parameter file'

        Ini_fail_on_not_found = .false.
    
        outroot = Ini_Read_String('output_root')
        if (outroot /= '') outroot = trim(outroot) // '_'
        
        call CAMB_SetDefParams(P)

        P%WantScalars = Ini_Read_Logical('get_scalar_cls')
        P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
        P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)
        
        P%OutputNormalization=outNone
        if (Ini_Read_Logical('COBE_normalize',.false.))  P%OutputNormalization=outCOBE
        output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

        P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

        P%WantTransfer=Ini_Read_Logical('get_transfer')
        
        P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)
   
        P%DoLensing = .false.
        if (P%WantCls) then
          if (P%WantScalars  .or. P%WantVectors) then
           P%Max_l = Ini_Read_Int('l_max_scalar')
           P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
           if (P%WantScalars) then
             P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
             if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
           end if
           if (P%WantVectors) then
            if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
            i = Ini_Read_Int('vector_mode')
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
           P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
           P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
          end if
        endif

                
!  Read initial parameters.
       
       call DarkEnergy_ReadParams(DefIni)

       P%h0     = Ini_Read_Double('hubble')
 
       if (Ini_Read_Logical('use_physical',.false.)) then 

        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
        P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
  
       else
       
        P%omegab = Ini_Read_Double('omega_baryon')
        P%omegac = Ini_Read_Double('omega_cdm')
        P%omegav = Ini_Read_Double('omega_lambda')
        P%omegan = Ini_Read_Double('omega_neutrino')

       end if

       P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)
       P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)
       P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')
       nmassive = Ini_Read_Double('massive_neutrinos')
       !Store fractional numbers in the massless total
       P%Num_Nu_massive   = int(nmassive+1e-6)
       P%Num_Nu_massless  = P%Num_Nu_massless + nmassive-P%Num_Nu_massive
   
       P%nu_mass_splittings = .true.
       P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
       if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'
       numstr = Ini_Read_String('nu_mass_degeneracies')
       if (numstr=='') then
         P%Nu_mass_degeneracies(1)= P%Num_nu_massive
       else
        read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
       end if
       numstr = Ini_read_String('nu_mass_fractions')
       if (numstr=='') then
        P%Nu_mass_fractions(1)=1  
        if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
       else
        read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
       end if

       if (P%NonLinear==NonLinear_lens .and. P%DoLensing) then
          if (P%WantTransfer) &
             write (*,*) 'overriding transfer settings to get non-linear lensing'
          P%WantTransfer  = .true.
          call Transfer_SetForNonlinearLensing(P%Transfer)
          P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
       
       else if (P%WantTransfer)  then
        P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%num_redshifts =  Ini_Read_Int('transfer_num_redshifts')
        
        transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
        if (P%transfer%num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        do i=1, P%transfer%num_redshifts
             P%transfer%redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
             transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
             MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
             
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


        P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)
                
       else
         P%transfer%high_precision = .false.
       endif
  
        Ini_fail_on_not_found = .false. 
  
        call Reionization_ReadParams(P%Reion, DefIni)
        call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors) 
        call Recombination_ReadParams(P%Recomb, DefIni)
        if (Ini_HasKey('recombination')) then
         i = Ini_Read_Int('recombination',1)
         if (i/=1) stop 'recombination option deprecated'
        end if
        
        call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)
        
        if (P%WantScalars .or. P%WantTransfer) then
            P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
            if (P%Scalar_initial_condition == initial_vector) then
                P%InitialConditionVector=0
              numstr = Ini_Read_String('initial_vector',.true.)
              read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
            end if
            if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
        end if
        
       if (P%WantScalars) then
          ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
          LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
          LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
          if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        end if
        if (P%WantTensors) then
          TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
         if (P%WantScalars)  then
          TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
          LensedTotFileName = Ini_Read_String('lensed_total_output_file')
          if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
         end if
        end if
        if (P%WantVectors) then
          VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
        end if
         
#ifdef WRITE_FITS
        if (P%WantCls) then
        FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
        if (FITSfilename /='') then
        inquire(file=FITSfilename, exist=bad)
        if (bad) then
         open(unit=18,file=FITSfilename,status='old')
         close(18,status='delete')
        end if
       end if
        end if
#endif        
       

       Ini_fail_on_not_found = .false. 

!optional parameters controlling the computation

       P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
       P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
       P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
       P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

       version_check = Ini_Read_String('version_check')
       if (version_check == '') then
          !tag the output used parameters .ini file with the version of CAMB being used now
           call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
       else if (version_check /= version) then
           write(*,*) 'WARNING: version_check does not match this CAMB version'
        end if
       !Mess here to fix typo with backwards compatibility
       if (Ini_HasKey('do_late_rad_trunction')) then
         DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
         if (Ini_HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
       else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
       end if
       DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
       FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)
       
       P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

       ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
       AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
       lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
       HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)
       use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)
       if (HighAccuracyDefault) then
         P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
       end if
       DoTensorNeutrinos = DoTensorNeutrinos .or. HighAccuracyDefault
       if (do_bispectrum) then
        lSampleBoost   = 50
       else
        lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
       end if
       if (outroot /= '') then
         if (InputFile /= trim(outroot) //'params.ini') then   
          call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
         else
          write(*,*) 'Output _params.ini not created as would overwrite input'    
         end if
       end if

       call Ini_Close

       if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

#ifdef RUNIDLE
       call SetIdle
#endif 

       if (FeedbackLevel > 0) then
         Age = CAMB_GetAge(P) 
         write (*,'("Age of universe/GYr  = ",f7.3)') Age  
       end if 

       if (global_error_flag==0) call CAMB_GetResults(P)
       if (global_error_flag/=0) then
        write (*,*) 'Error result '//trim(global_error_message)
        stop
       endif
    
        if (P%WantTransfer .and. .not. (P%NonLinear==NonLinear_lens .and. P%DoLensing)) then
         call Transfer_SaveToFiles(MT,TransferFileNames)
         call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
         if ((P%OutputNormalization /= outCOBE) .or. .not. P%WantCls)  call Transfer_output_sig8(MT)
        end if

        if (P%WantCls) then
  
         if (P%OutputNormalization == outCOBE) then

            if (P%WantTransfer) call Transfer_output_Sig8AndNorm(MT)
           
          end if

         call output_cl_files(ScalarFileName, TensorFileName, TotalFileName, &
              LensedFileName, LensedTotFilename, output_factor)
              
         call output_lens_pot_files(LensPotentialFileName, output_factor)

         if (P%WantVectors) then
           call output_veccl_files(VectorFileName, output_factor)
         end if


#ifdef WRITE_FITS
         if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
        end if

        call CAMB_cleanup             
  
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
