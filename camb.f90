    !Interface module for CAMB. Call CAMB_GetResults to do the work.

    module CAMB
    use Precision
    use CambSettings
    use GaugeInterface
    use InitialPower
    use Reionization
    use Recombination
    use lensing
    use DarkEnergyFluid
    implicit none
    contains

    subroutine CAMB_GetTransfers(Params, OutData, error, onlytransfer)
    type(CAMBparams) :: Params
    type(CAMBstate) :: OutData
    integer :: error !Zero if OK
    logical onlytransfer

    call CAMB_FreeCAMBstate(OutData)
    call CAMB_GetResults(OutData, Params, error,onlytransfer)

    end subroutine CAMB_GetTransfers

    subroutine CAMB_FreeCAMBstate(Dat)
    type (CAMBstate) :: Dat

    call Free_ClTransfer(Dat%CLdata%CTransScal)
    call Free_ClTransfer(Dat%ClData%CTransVec)
    call Free_ClTransfer(Dat%ClData%CTransTens)
    call Transfer_Free(Dat%MT)

    end subroutine CAMB_FreeCAMBstate


    subroutine CAMB_TransfersToPowers(CData)
    use CAMBmain
    use lensing
    type (CAMBstate) :: CData

    call SetActiveState(CData)
    CData%OnlyTransfer = .false.
    call CData%CP%InitPower%Init(CData%CP, Cdata%omegak)
    if (global_error_flag/=0) return

    if (CData%CP%WantCls) then
        call ClTransferToCl(CData)
        if (CP%DoLensing .and. global_error_flag==0) call lens_Cls(Cdata)
        if (global_error_flag/=0) return
    end if

    if (CData%CP%WantTransfer) call Transfer_Get_sigmas(Cdata, Cdata%MT)

    end subroutine CAMB_TransfersToPowers


    !Call this routine with a set of parameters to generate the results you want.
    subroutine CAMB_GetResults(OutData, Params, error, onlytransfer)
    use CAMBmain
    use lensing
    use Bispectrum
    use Errors
    type(CAMBstate)  :: OutData
    type(CAMBparams) :: Params
    integer, optional :: error !Zero if OK
    logical, optional :: onlytransfer
    type(CAMBparams) P
    logical :: call_again

    global_error_flag = 0
    call_again = .false.
    call SetActiveState(OutData)
    OutData%OnlyTransfer = DefaultFalse(onlytransfer)
    if (Params%WantCls .and. Params%WantScalars) then
        P = Params
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        P%WantTensors = .false.
        P%WantVectors = .false.
        if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. (P%DoLensing .or. num_redshiftwindows > 0)) then
            P%WantTransfer  = .true.
        end if
        call OutData%CAMBParams_Set(P)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return
        end if
        call_again = .true.
    end if

    if (Params%WantCls .and. Params%WantTensors) then
        P=Params
        P%WantTransfer = .false.
        P%Transfer%high_precision = .false.
        P%WantScalars = .false.
        P%WantVectors = .false.
        call OutData%CAMBParams_Set(P, call_again=call_again)
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
        call OutData%CAMBParams_Set(P, call_again=call_again)
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
        call OutData%CAMBParams_Set(P, call_again=call_again)
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
    OutData%CP%DerivedParameters = Params%DerivedParameters

    if (.not. OutData%OnlyTransfer .and. Params%WantCls .and. Params%WantScalars) then
        if (Params%DoLensing .and. global_error_flag==0) then
            call lens_Cls(OutData)
        end if

        if (do_bispectrum .and. global_error_flag==0) &
            call GetBispectrum(OutData%CLData%CTransScal)
    end if

    end subroutine CAMB_GetResults


    !Return real (NOT double precision) arrays of the computed CMB  Cls
    !Output is l(l+1)C_l/2pi
    !If GC_Conventions = .false. use E-B conventions (as the rest of CAMB does)
    subroutine CAMB_GetCls(State, Cls, lmax,GC_conventions)
    Type(CAMBstate) :: State
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
    !Return age in gigayears, returns -1 on error
    type(CAMBparams), intent(in) :: P
    real(dl) CAMB_GetAge
    integer error
    Type(CAMBstate) :: State

    call  State%CAMBParams_Set(P, error, .false.)

    if (error/=0) then
        CAMB_GetAge = -1
    else
        CAMB_GetAge = State%DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
    end if

    end function CAMB_GetAge


    function CAMB_GetZreFromTau(P, tau)
    type(CAMBparams) :: P
    real(dl) tau
    real(dl) CAMB_GetZreFromTau
    integer error
    type(CAMBstate) :: State

    P%Reion%use_optical_depth = .true.
    P%Reion%optical_depth = tau
    call State%CAMBParams_Set(P,error)
    if (error/=0)  then
        CAMB_GetZreFromTau = -1
    else
        CAMB_GetZreFromTau = State%CP%Reion%redshift
    end if

    end function CAMB_GetZreFromTau


    subroutine CAMB_SetDefParams(P)
    use NonLinear
    type(CAMBparams), intent(inout), target :: P
    type(CAMBparams) :: emptyP

    P= emptyP !Set default values set in type definitions
    P%Nu_mass_numbers=0
    P%Nu_mass_degeneracies=0
    P%Nu_mass_fractions=0

    allocate(THalofit::P%NonLinearModel)
    allocate(TDarkEnergyFluid::P%DarkEnergy)
    allocate(TInitialPowerLaw::P%InitPower)

    end subroutine CAMB_SetDefParams

   
    end module CAMB
