!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output. 

    program driver
        use IniFile
        use CAMB
        use LambdaGeneral
        use Lensing
        use NonLinear
        use AMLutils
        use constants
        implicit none
      
        Type(CAMBparams) P
        
        character(LEN=Ini_max_string_len) numstr, S, VectorFileName, &
            InputFile, ScalarFileName, ScalarCovFileName,TensorFileName, TotalFileName, LensedFileName
        integer i
        character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
               MatterPowerFileNames(max_transfer_redshifts), outroot, &
               TransferClFileNames(max_transfer_redshifts)
        

        real(dl) output_factor, Age
        Type (TRedWin), pointer :: RedWin

#ifdef WRITE_FITS
       character(LEN=Ini_max_string_len) FITSfilename
#endif

        logical bad
        logical :: DoCounts = .false.


        InputFile = ''

        if (GetParamCount() /= 0)  InputFile = GetParam(1)
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
        
        P%Want_CMB =  Ini_Read_Logical('want_CMB',.true.)

        if (P%WantScalars) then
         num_redshiftwindows = Ini_Read_Int('num_redshiftwindows',0)
        else
         num_redshiftwindows = 0
        end if
        if (num_redshiftwindows>0) then
         DoRedshiftLensing = Ini_Read_Logical('DoRedshiftLensing',.false.)
         limber_windows = Ini_Read_Logical('limber_windows',.false.)
        end if
        Do21cm = Ini_Read_Logical('Do21cm', .false.)
        num_extra_redshiftwindows = 0
        do i=1, num_redshiftwindows
             RedWin => Redshift_w(i)
             call InitRedshiftWindow(RedWin)
             write (numstr,*) i 
             numstr=adjustl(numstr)
             RedWin%Redshift = Ini_Read_Double('redshift('//trim(numstr)//')')
             S = Ini_Read_String('redshift_kind('//trim(numstr)//')')
             if (S=='21cm') then
              RedWin%kind = window_21cm
            elseif (S=='counts') then
              RedWin%kind = window_counts
            elseif (S=='lensing') then
              RedWin%kind = window_lensing
            else
              write (*,*) i, 'Error: unknown type of window '//trim(S)
              stop
            end if
            RedWin%a = 1/(1+RedWin%Redshift)
            if (RedWin%kind /= window_21cm) then  
              RedWin%sigma = Ini_Read_Double('redshift_sigma('//trim(numstr)//')')
              RedWin%sigma_z = RedWin%sigma
            else
              Do21cm = .true.
              RedWin%sigma = Ini_Read_Double('redshift_sigma_Mhz('//trim(numstr)//')')
              if (RedWin%sigma < 0.003) then
               write(*,*) 'WARNING:Window very narrow.'
               write(*,*) ' --> use transfer functions and transfer_21cm_cl =T ?'
              end if
              !with 21cm widths are in Mhz, make dimensionless scale factor
              RedWin%sigma = RedWin%sigma/(f_21cm/1e6)
              RedWin%sigma_z = RedWin%sigma*(1+RedWin%RedShift)**2
              write(*,*) i,'delta_z = ', RedWin%sigma_z
            end if
            if (RedWin%kind == window_counts) then
              DoCounts = .true.
              RedWin%bias = Ini_Read_Double('redshift_bias('//trim(numstr)//')')
              RedWin%dlog10Ndm = Ini_Read_Double('redshift_dlog10Ndm('//trim(numstr)//')',0.d0)
              if (DoRedshiftLensing) then 
               num_extra_redshiftwindows=num_extra_redshiftwindows+1
               RedWin%mag_index = num_extra_redshiftwindows
              end if
             end if
        end do

 
        if (Do21cm) then
          line_basic = Ini_Read_Logical('line_basic')
          line_distortions = Ini_read_Logical('line_distortions')
          line_extra = Ini_Read_Logical('line_extra')

          line_phot_dipole = Ini_read_Logical('line_phot_dipole')
          line_phot_quadrupole = Ini_Read_Logical('line_phot_quadrupole')
          line_reionization = Ini_Read_Logical('line_reionization')

          use_mK = Ini_read_Logical('use_mK')
          if (DebugMsgs) then
           write (*,*) 'Doing 21cm'
           write (*,*) 'dipole = ',line_phot_dipole, ' quadrupole =', line_phot_quadrupole  
          end if    
        else
         line_extra = .false.
        end if

        if (DoCounts) then
          counts_density = Ini_read_Logical('counts_density')
          counts_redshift = Ini_read_Logical('counts_redshift')
          counts_radial = Ini_read_Logical('counts_radial')
          counts_evolve = Ini_read_Logical('counts_evolve')
          counts_timedelay = Ini_read_Logical('counts_timedelay')
          counts_ISW = Ini_read_Logical('counts_ISW')
          counts_potential = Ini_read_Logical('counts_potential')
          counts_velocity = Ini_read_Logical('counts_velocity')
          
        end if
          
        P%OutputNormalization=outNone
        if (Ini_Read_Logical('COBE_normalize',.false.))  P%OutputNormalization=outCOBE
        output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

        P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

        P%WantTransfer=Ini_Read_Logical('get_transfer')
        
        P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

        evolve_delta_xe = Ini_read_Logical('evolve_delta_xe', .false.)
   
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
       
       w_lam = Ini_Read_Double('w', -1.d0)   
       cs2_lam = Ini_Read_Double('cs2_lam',1.d0)

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
       P%Num_Nu_massive   = Ini_Read_Double('massive_neutrinos')
   
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

       if (P%NonLinear==NonLinear_lens .and. (P%DoLensing .or. num_redshiftwindows>0)) then
          if (P%WantTransfer) &
             write (*,*) 'overriding transfer settings to get non-linear lensing'
          P%WantTransfer  = .true.
          call Transfer_SetForNonlinearLensing(P%Transfer, P%Max_eta_k)
          P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
       
       else if (P%WantTransfer)  then
        P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%num_redshifts =  Ini_Read_Int('transfer_num_redshifts')
        if (Do21cm) transfer_21cm_cl = Ini_Read_Logical('transfer_21cm_cl',.false.)
        if (transfer_21cm_cl .and. P%transfer%kmax > 800) then
        !Actually line widths are important at significantly larger scales too
         write (*,*) 'WARNING: kmax very large. '
         write(*,*) ' -- Neglected line width effects will dominate'
        end if
        if (P%transfer%num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        do i=1, P%transfer%num_redshifts
             write (numstr,*) i 
             numstr=adjustl(numstr)
             P%transfer%redshifts(i)  = Ini_Read_Double('transfer_redshift('//trim(numstr)//')',0._dl)
             TransferFileNames(i)     = Ini_Read_String('transfer_filename('//trim(numstr)//')')
             MatterPowerFilenames(i)  = Ini_Read_String('transfer_matterpower('//trim(numstr)//')')
             if (Do21cm) then
               TransferClFileNames(i)     = Ini_Read_String('transfer_cl_filename('//trim(numstr)//')')
               if (TransferClFileNames(i) == '') then
                  TransferClFileNames(i) =  trim(numcat('sharp_cl_',i))//'.dat' 
               end if
             end if
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
            
             if (TransferClFileNames(i)/= '') &
                   TransferClFileNames(i) = trim(outroot)//TransferClFileNames(i)
          
        end do
        P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)
                
       else
         P%transfer%high_precision = .false.
       endif
  
     
      
          call Reionization_ReadParams(P%Reion, DefIni)
          call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors) 
          call Recombination_ReadParams(P%Recomb, DefIni)


           i = Ini_Read_Int('recombination',1)
           if (i>1) then
             stop 'recombination option deprecated'
           end if
      
            if (P%WantScalars .or. P%WantTransfer) then
            P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
            if (P%Scalar_initial_condition == initial_vector) then
                P%InitialConditionVector=0
              numstr = Ini_Read_String('initial_vector',.true.)
              read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
            end if
        end if

        
       if (P%WantScalars) then
          ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
          LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
          ScalarCovFileName = trim(outroot)//Ini_Read_String('scalar_covariance_output_file')
        end if
        if (P%WantTensors) then
          TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
         if (P%WantScalars) TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
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
        
       !Mess here to fix typo with backwards compatibility
       if (Ini_Read_String('do_late_rad_trunction') /= '') then
         DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
         if (Ini_Read_String('do_late_rad_truncation')/='') stop 'check do_late_rad_xxxx'
       else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
       end if
       DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',.false.)
       FeedbackLevel = Ini_Read_Int('feedback_level',0)
       
       P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

       ThreadNum      = Ini_Read_Int('number_of_threads',0)
       AccuracyBoost  = Ini_Read_Double('accuracy_boost',1.d0)
       lAccuracyBoost = Ini_Read_Double('l_accuracy_boost',1.d0)
       lSampleBoost   = Ini_Read_Double('l_sample_boost',1.d0)

       if (outroot /= '') then
         call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
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

       call CAMB_GetResults(P)
    
        if (P%WantTransfer .and. .not. (P%NonLinear==NonLinear_lens &
            .and. (P%DoLensing .or. num_redshiftwindows>0))) then

         call Transfer_SaveToFiles(MT,TransferFileNames)
         call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
!        call Transfer_SaveMatterPower(MT,MatterPowerFileNames, .true.)
         
         if (do21cm .and. transfer_21cm_cl) call Transfer_Get21cmCls(MT,TransferClFileNames) 

         if ((P%OutputNormalization /= outCOBE) .or. .not. P%WantCls)  call Transfer_output_sig8(MT)
        end if

        if (P%WantCls) then
  
         if (P%OutputNormalization == outCOBE) then

            if (P%WantTransfer) call Transfer_output_Sig8AndNorm(MT)
           
          end if

         call output_cl_files(ScalarFileName, ScalarCovFileName,TensorFileName, TotalFileName, &
              LensedFileName, output_factor)

         if (P%WantVectors) then
           call output_veccl_files(VectorFileName, output_factor)
         end if


#ifdef WRITE_FITS
         if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
        end if

             
  
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
