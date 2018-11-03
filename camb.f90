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

    subroutine CAMB_GetTransfers(Params, OutData, error)
    type(CAMBparams) :: Params
    type(CAMBstate) :: OutData
    integer :: error !Zero if OK

    call Transfer_Free(OutData%MT)

    call Free_ClTransfer(OutData%CLData%CTransScal)
    call Free_ClTransfer(OutData%CLData%CTransVec)
    call Free_ClTransfer(OutData%CLData%CTransTens)

    call CAMB_GetResults(OutData, Params, error)

    end subroutine CAMB_GetTransfers

    subroutine CAMB_InitCAMBstate(Dat)
    type (CAMBstate) :: Dat

    nullify(Dat%CLData%CTransScal%Delta_p_l_k)
    nullify(Dat%CLData%CTransVec%Delta_p_l_k)
    nullify(Dat%CLData%CTransTens%Delta_p_l_k)
    nullify(Dat%MT%sigma_8,Dat%MT%TransferData,Dat%MT%q_trans)

    end subroutine CAMB_InitCAMBstate


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
    call CData%CP%InitPower%Init(Cdata%curv)
    if (global_error_flag/=0) return

    if (CData%CP%WantCls) then
        call ClTransferToCl(CData)
        if (CP%DoLensing .and. global_error_flag==0) call lens_Cls(Cdata)
        if (global_error_flag/=0) return
    end if

    if (CData%CP%WantTransfer) call Transfer_Get_sigmas(Cdata%MT, Cdata%CP)

    end subroutine CAMB_TransfersToPowers


    !Call this routine with a set of parameters to generate the results you want.
    subroutine CAMB_GetResults(OutData, Params, error)
    use CAMBmain
    use lensing
    use Bispectrum
    use Errors
    type(CAMBstate)  :: OutData
    type(CAMBparams) :: Params
    integer, optional :: error !Zero if OK
    type(CAMBparams) P
    logical :: call_again

    global_error_flag = 0
    call_again = .false.

    call SetActiveState(OutData)
    if (Params%WantCls .and. Params%WantScalars) then
        P = Params
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        P%WantTensors = .false.
        P%WantVectors = .false.
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

    if (.not. CP%OnlyTransfers) then
        if (CP%DoLensing .and. global_error_flag==0) then
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
        CAMB_GetAge = DeltaPhysicalTimeGyr(State,0.0_dl,1.0_dl)
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
    type(CAMBparams), intent(inout) :: P
    type(CAMBparams) :: emptyP

    P= emptyP !Set default values set in type definitions
    P%Nu_mass_numbers=0

    allocate(THalofit::P%NonLinearModel)
    allocate(TDarkEnergyFluid::P%DarkEnergy)
    allocate(TInitialPowerLaw::P%InitPower)

    !Set up transfer just enough to get sigma_8 OK
    P%Transfer%num_redshifts = 1
    P%Transfer%redshifts=0
    !JD 08/13 CAMB Fix for for nonlinear lensing of CMB + MPK compatibility
    P%Transfer%PK_num_redshifts=1
    P%Transfer%PK_redshifts=0
    P%Transfer%NLL_num_redshifts=0 !AL 11/13, def to zero
    P%Transfer%NLL_redshifts=0
    !End JD

    end subroutine CAMB_SetDefParams

    subroutine CAMB_SetNeutrinoHierarchy(P, omnuh2, omnuh2_sterile, nnu, neutrino_hierarchy, num_massive_neutrinos)
    use constants
    use MathUtils
    type(CAMBparams), intent(inout) :: P
    real(dl), intent(in) :: omnuh2, omnuh2_sterile, nnu
    integer, intent(in) :: neutrino_hierarchy
    integer, intent(in), optional :: num_massive_neutrinos  !for degenerate hierarchy
    integer, parameter :: neutrino_hierarchy_normal = 1, neutrino_hierarchy_inverted = 2, neutrino_hierarchy_degenerate = 3
    real(dl) normal_frac, m3, neff_massive_standard, mnu, m1

    if (omnuh2==0) return
    P%Nu_mass_eigenstates=0
    if ( omnuh2 > omnuh2_sterile) then
        normal_frac =  (omnuh2-omnuh2_sterile)/omnuh2
        if (neutrino_hierarchy == neutrino_hierarchy_degenerate) then
            neff_massive_standard = num_massive_neutrinos*default_nnu/3
            P%Num_Nu_Massive = num_massive_neutrinos
            P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
            if (nnu > neff_massive_standard) then
                P%Num_Nu_Massless = nnu - neff_massive_standard
            else
                P%Num_Nu_Massless = 0
                neff_massive_standard=nnu
            end if
            P%Nu_mass_numbers(P%Nu_mass_eigenstates) = num_massive_neutrinos
            P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = neff_massive_standard
            P%Nu_mass_fractions(P%Nu_mass_eigenstates) = normal_frac
        else
            !Use normal or inverted hierarchy, approximated as two eigenstates in physical regime, 1 at minimum an below
            mnu = (omnuh2 - omnuh2_sterile)*neutrino_mass_fac / (default_nnu / 3) ** 0.75_dl
            if (neutrino_hierarchy == neutrino_hierarchy_normal) then
                if (mnu > mnu_min_normal + 1e-4_dl) then
                    !Two eigenstate approximation.
                    m1=Newton_Raphson(0._dl, mnu, sum_mnu_for_m1, mnu, 1._dl)
                    P%Num_Nu_Massive = 3
                else
                    !One eigenstate
                    P%Num_Nu_Massive = 1
                end if
            else if (neutrino_hierarchy == neutrino_hierarchy_inverted) then
                if (mnu > sqrt(delta_mnu31)+sqrt(delta_mnu31+delta_mnu21) + 1e-4_dl ) then
                    !Valid case, two eigenstates
                    m1=Newton_Raphson(sqrt(delta_mnu31), mnu, sum_mnu_for_m1, mnu, -1._dl)
                    P%Num_Nu_Massive = 3
                else
                    !Unphysical low mass case: take one (2-degenerate) eigenstate
                    P%Num_Nu_Massive = 2
                end if
            else
                error stop 'Unknown neutrino_hierarchy setting'
            end if
            neff_massive_standard = P%Num_Nu_Massive *default_nnu/3
            if (nnu > neff_massive_standard) then
                P%Num_Nu_Massless = nnu - neff_massive_standard
            else
                P%Num_Nu_Massless = 0
                neff_massive_standard=nnu
            end if
            if (P%Num_Nu_Massive==3) then
                !two with mass m1, one with m3
                P%Nu_mass_eigenstates = 2
                P%Nu_mass_degeneracies(1) = neff_massive_standard*2/3._dl
                P%Nu_mass_degeneracies(2) = neff_massive_standard*1/3._dl
                m3 = mnu - 2*m1
                P%Nu_mass_fractions(1) = 2*m1/mnu*normal_frac
                P%Nu_mass_fractions(2) = m3/mnu*normal_frac
                P%Nu_mass_numbers(1) = 2
                P%Nu_mass_numbers(2) = 1
            else
                P%Nu_mass_degeneracies(1) = neff_massive_standard
                P%Nu_mass_numbers(1) = P%Num_Nu_Massive
                P%Nu_mass_eigenstates = 1
                P%Nu_mass_fractions(1) = normal_frac
            end if
        end if
    else
        neff_massive_standard=0
    end if
    if (omnuh2_sterile>0) then
        if (nnu<default_nnu) call MpiStop('nnu < 3.046 with massive sterile')
        P%Num_Nu_Massless = default_nnu - neff_massive_standard
        P%Num_Nu_Massive=P%Num_Nu_Massive+1
        P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
        P%Nu_mass_numbers(P%Nu_mass_eigenstates) = 1
        P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = max(1d-6,nnu - default_nnu)
        P%Nu_mass_fractions(P%Nu_mass_eigenstates) = omnuh2_sterile/omnuh2
    end if
    end subroutine CAMB_SetNeutrinoHierarchy


    !Stop with error is not good
    function CAMB_ValidateParams(P) result(OK)
    type(CAMBparams), intent(in) :: P
    logical OK

    OK = .true.
    if (.not. P%WantTransfer .and. .not. P%WantCls) then
        OK = .false.
        write(*,*) 'There is nothing to do! Do transfer functions or Cls.'
    end if

    if (P%h0 < 20._dl.or.P%h0 > 100._dl) then
        OK = .false.
        write(*,*) '  Warning: H0 has units of km/s/Mpc. You have:', P%h0
    end if
    if (P%tcmb < 2.7d0.or.P%tcmb > 2.8d0) then
        write(*,*) '  Warning: Tcmb has units of K.  Your have:', P%tcmb
    end if

    if (P%yhe < 0.2d0.or.P%yhe > 0.8d0) then
        OK = .false.
        write(*,*) &
            '  Warning: YHe is the Helium fraction of baryons.', &
            '  Your have:', P%yhe
    end if
    if (P%Num_Nu_massive < 0) then
        OK = .false.
        write(*,*) &
            'Warning: Num_Nu_massive is strange:',P%Num_Nu_massive
    end if
    if (P%Num_Nu_massless < 0) then
        OK = .false.
        write(*,*) &
            'Warning: Num_nu_massless is strange:', P%Num_Nu_massless
    end if
    if (P%Num_Nu_massive < 1 .and. P%omegan > 0.0) then
        OK = .false.
        write(*,*) &
            'Warning: You have omega_neutrino > 0, but no massive species'
    end if


    if (P%omegab<0.001 .or. P%omegac<0 .or. P%omegab>1 .or. P%omegac>3) then
        OK = .false.
        write(*,*) 'Your matter densities are strange'
    end if

    if (P%WantScalars .and. P%Max_eta_k < P%Max_l .or.  &
        P%WantTensors .and. P%Max_eta_k_tensor < P%Max_l_tensor) then
        OK = .false.
        write(*,*) 'You need Max_eta_k larger than Max_l to get good results'
    end if

    call P%Reion%Validate(OK)
    call P%Recomb%Validate(OK)

    if (P%WantTransfer) then
        if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
            OK = .false.
            write(*,*) 'Maximum ',  max_transfer_redshifts, &
                'redshifts. You have: ', P%transfer%num_redshifts
        end if
        if (P%transfer%kmax < 0.01 .or. P%transfer%kmax > 50000 .or. &
            P%transfer%k_per_logint>0 .and.  P%transfer%k_per_logint <1) then
            !            OK = .false.
            write(*,*) 'Strange transfer function settings.'
        end if
        if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
            OK = .false.
            write(*,*) 'Maximum ',  max_transfer_redshifts, &
                'redshifts. You have: ', P%transfer%num_redshifts
        end if
    end if

    end function CAMB_ValidateParams

    end module CAMB
