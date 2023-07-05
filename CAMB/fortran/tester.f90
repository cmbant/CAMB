    !Simple demo program to get the scalar and tensor temperature Cls, and print them out (and the sum)
    !Usually it's easier to do things like this via the python wrapper

    program tester
    use CAMB
    implicit none
    integer l
    real(dl) ratio

    type(CAMBparams) P !defined in model.f90
    type(CAMBdata) results !contains computed quantities and functions for results (results.f90)

    !Need to set default classes for initial power spectrum, recombination etc.
    call CAMB_SetDefParams(P)

    !Set one massive neutrino with mass~ 0.06eV
    call P%SetNeutrinoHierarchy(0.00064_dl, 0._dl, 3.046_dl, neutrino_hierarchy_normal)


    P%ombh2  = .0222_dl
    P%omk  = 0._dl
    P%H0      = 67._dl
    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
        InitPower%As = 2.1e-9
        InitPower%ns  = 1
        InitPower%r = 1
        !we don't use r here since we generate the Cls separately
        !so set to 1, and then put in the ratio after calculating the Cls
    end select

    P%WantScalars = .true.
    P%WantTensors = .true.

    P%Max_l=2500
    P%Max_eta_k=6000
    P%Max_l_tensor=200
    P%Max_eta_k_tensor=500

    call CAMB_GetResults(results, P)

    ratio =0.1
    associate(CL_scalar => results%CLData%CL_scalar, &
        CL_tensor => results%CLData%CL_tensor)
        do l=2,P%Max_l
            !print out scalar and tensor temperature, then sum. Units are dimensionless here.
            if (l <= P%Max_l_tensor) then
                !The indices of the Cl_xxx arrays are L, Cl type
                write(*,'(1I5,3E15.5)') l, Cl_scalar(l, C_Temp), ratio*Cl_tensor(l,CT_Temp), &
                    ratio*Cl_tensor(l,C_Temp)+Cl_scalar(l,C_Temp)
            else
                write(*,'(1I5,3E15.5)') l,  Cl_scalar(l,C_Temp),0._dl,Cl_scalar(l,C_Temp)
            end if
        end do
    end associate

    end program Tester


