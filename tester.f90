!Simple program to get the scalar and tensor Cls, and print them out (and the sum)

      program tester
       use CAMB
        implicit none
        integer l
        real(dl) ratio

        type(CAMBparams)  P !defined in ModelParams in modules.f90

        
        call CAMB_SetDefParams(P)

        P%omegab  = .045
        P%omegac  = 0.355
        P%omegav  = 0.6
        P%omegan  = 0.0
        P%H0      = 65
   
        P%InitPower%nn     = 1 !number of initial power spectra
        P%InitPower%an(1)  = 1 !scalar spectral index
        P%InitPower%ant(1) = 0 !tensor spectra index
        P%InitPower%rat(1) = 1 !ratio of initial amplitudes
         !actually we don't use this here since we generate the Cls separately
         !so set to 1, and then put in the ratio after calculating the Cls

 
        P%OutputNormalization = outNone
      
        !Generate scalars first so that start with maximum Max_l that is used
        P%WantScalars = .true.
        P%WantTensors = .true.
      
        P%Max_l=1500
        P%Max_eta_k=3000
        P%Max_l_tensor=200
        P%Max_eta_k_tensor=500

        P%AccuratePolarization = .false. !We are only interested in the temperature here

        call CAMB_GetResults(P) 

        ratio =0.1

        do l=2,P%Max_l
           !print out scalar and tensor temperature, then sum
           if (l <= P%Max_l_tensor) then
              !The indices of the Cl_xxx arrays are l, initial power spectrum index, Cl type
             write(*,'(1I5,3E15.5)') l, Cl_scalar(l,1, C_Temp), ratio*Cl_tensor(l,1,CT_Temp), &
                   ratio*Cl_tensor(l,1,C_Temp)+Cl_scalar(l,1,C_Temp) 
           else
             write(*,'(1I5,3E15.5)') l,  Cl_scalar(l,1, C_Temp),0._dl,Cl_scalar(l,1,C_Temp) 
           end if
        end do   


        end program Tester


