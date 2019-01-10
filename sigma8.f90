    !Simple test program to print out sigma_8 as a function of the CDM density
    !However usually it is much easier to do this kind of thing via the python wrapper.
    program GetSigma8
    use CAMB
    implicit none
    integer i

    type(CAMBparams) P !defined in model.f90
    type(CAMBdata) results !contains computed quantities and functions for results (results.f90)

    !Need to set default classes for initial power spectrum, recombination etc.
    call CAMB_SetDefParams(P)

    !Set one massive neutrino with mass~ 0.06eV
    call P%SetNeutrinoHierarchy(0.00064_dl, 0._dl, 3.046_dl, neutrino_hierarchy_normal)

    P%WantTransfer= .true.

    P%WantCls = .false.

    P%ombh2  = .0222_dl
    P%omk  = 0._dl
    P%H0      = 67._dl
    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
        InitPower%As = 2.1e-9
        InitPower%ns  = 1
    end select

    !these settings seem good enough for sigma8 to a percent or so
    P%Transfer%high_precision=.true.
    P%Transfer%kmax=0.5
    P%Transfer%k_per_logint=3
    P%Transfer%PK_num_redshifts=1
    P%Transfer%PK_redshifts(1)=0

    do i=1,10
        P%omch2 = 0.05_dl + i*0.01_dl

        call CAMB_GetResults(results, P)

        write (*,*) 'Omegac h^2 = ',real(P%omch2), 'sigma_8 = ', real(results%MT%sigma_8(1))
    end do

    end program GetSigma8


