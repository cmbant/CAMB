     !Simple test program to print out sigma_8 as a function of the CDM density
     program GetSigma8
        use CAMB
        implicit none
        integer i
    
        type(CAMBparams)  P !defined in ModelParams in modules.f90

        call CAMB_SetDefParams(P)

        P%WantTransfer= .true.

        P%WantCls = .false.

        P%omegab  = .045
        P%omegac  = 0.155
        P%omegav  = 0.8
        P%omegan  = 0.0
        P%H0      = 65
       
        P%InitPower%ScalarPowerAmp = 2e-9
        P%InitPower%nn     = 1 !number of initial power spectra
        P%InitPower%an(1)  = 1 !scalar spectral index
        P%InitPower%ant(1) = 0 !Not used here
        P%InitPower%rat(1) = 1 !ditto

        !these settings seem good enough for sigma8 to a percent or so
        P%Transfer%high_precision=.false.
        P%Transfer%kmax=0.5
        P%Transfer%k_per_logint=3
        P%Transfer%num_redshifts=1
        P%Transfer%redshifts(1)=0

        do i=1,10
         P%Omegav=P%Omegav-0.05
         P%Omegac=P%Omegac+0.05
        call CAMB_GetResults(P) 

        !Results are in the Transfer module in modules.f90
             
        write (*,*) 'Omc = ',real(P%Omegac),'OmLam=',real(P%Omegav) &
           , 'sigma_8 = ', real(MT%sigma_8(1,1))
        end do

        end program GetSigma8


