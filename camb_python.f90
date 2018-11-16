
    module handles
    use CAMB
    use Precision
    use CambSettings
    use iso_c_binding
    use DarkEnergyFluid
    use DarkEnergyPPF
    implicit none

    Type c_MatterTransferData
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        type(c_ptr) :: q_trans, sigma_8, sigma2_vdelta_8, TransferData
        integer :: sigma_8_size
        integer :: sigma2_vdelta_8_size
        integer :: TransferData_size(3)
    end Type c_MatterTransferData


    Type c_ClTransferData
        integer :: NumSources
        integer :: q_size
        type(c_ptr) :: q
        integer delta_size(3)
        type(c_ptr) Delta_p_l_k
        integer l_size
        type(c_ptr) ls
        !skip limber for now...
    end Type c_ClTransferData

    contains

    !SPECIAL BRIDGE ROUTINES FOR PYTHON

    subroutine CAMBdata_new(handle)
    type(c_ptr), intent(inout) :: handle
    Type (CAMBstate), pointer :: pSource,pCAMBstate

    if (associated(pSource)) then
        allocate(pCAMBstate, source=pSource)
    else
        allocate(pCAMBstate)
        call CAMB_InitCAMBstate(pCAMBstate)
        call CAMB_setDefParams(pCAMBstate%CP)
    end if
    handle = c_loc(pCAMBstate)

    end subroutine CAMBdata_new

    subroutine CAMBdata_free(cptr)
    type(c_ptr)  :: cptr
    Type (CAMBstate), pointer :: pCAMBstate

    call c_f_pointer(cptr, pCAMBstate)
    call CAMB_FreeCAMBstate(pCAMBstate)
    deallocate(pCAMBstate)

    end subroutine CAMBdata_free

    subroutine CAMBdata_setParams(State, Params)
    Type (CAMBstate), target :: State
    type(CAMBparams) :: Params

    State%CP = Params

    end subroutine CAMBdata_setParams

    subroutine CAMBdata_getParams(State, handle)
    Type (CAMBstate), target :: State
    type(c_ptr), intent(out)  ::  handle

    handle = c_loc(State%CP)

    end subroutine CAMBdata_getParams

    subroutine CAMBdata_getDerived(State, derived)
    Type (CAMBstate), target :: State
    real(dl) :: derived(nthermo_derived)

    derived = State%ThermoDerivedParams

    end subroutine CAMBdata_getDerived

    integer function CAMBdata_getlmax_lensed(State)
    Type (CAMBstate), target :: State
    CAMBdata_getlmax_lensed = State%CLdata%lmax_lensed
    end function CAMBdata_getlmax_lensed

    function CAMBdata_GetTransfers(State, Params, onlytransfer) result(error)
    Type (CAMBstate):: State
    type(CAMBparams) :: Params, P
    logical(kind=c_bool)  :: onlytransfer
    integer :: error

    call SetActiveState(State)
    P = Params
    if (P%DoLensing .and. (P%NonLinear == NonLinear_lens .or. P%NonLinear == NonLinear_both)) then
        P%WantTransfer = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if
    call Transfer_SortAndIndexRedshifts(P%Transfer)
    error = 0

    P%OnlyTransfers = onlytransfer
    call CAMB_GetTransfers(P, State, error)

    end function CAMBdata_GetTransfers

    subroutine CAMBdata_SetParamsForBackground(State, P)
    Type (CAMBstate):: State
    type(CAMBparams) :: P

    global_error_flag = 0
    call State%CAMBParams_Set(P)
    end subroutine CAMBdata_SetParamsForBackground

    function CAMBdata_CalcBackgroundTheory(State, P) result(error)
    use cambmain, only: initvars
    Type (CAMBstate):: State
    type(CAMBparams) :: P
    integer error

    global_error_flag = 0
    call State%CAMBParams_Set(P)
    if (global_error_flag==0) call InitVars(State) !calculate thermal history, e.g. z_drag etc.
    error=global_error_flag

    end function CAMBdata_CalcBackgroundTheory


    subroutine CAMBdata_MatterTransferData(State, cData)
    Type(CAMBstate), target :: State
    Type(c_MatterTransferData) :: cData

    cData%num_q_trans = State%MT%num_q_trans
    cData%q_trans = c_loc(State%MT%q_trans)
    cData%sigma_8 = c_loc(State%MT%sigma_8)
    cData%sigma2_vdelta_8 = c_loc(State%MT%sigma2_vdelta_8)
    cData%TransferData = c_loc(State%MT%TransferData)
    cData%q_trans = c_loc(State%MT%q_trans)
    cData%sigma_8_size = size(State%MT%sigma_8)
    cData%sigma2_vdelta_8_size = size(State%MT%sigma2_vdelta_8)
    cData%TransferData_size = shape(State%MT%TransferData)

    end subroutine CAMBdata_MatterTransferData

    subroutine CAMBdata_ClTransferData(State, cData, i)
    Type(CAMBstate), target :: State
    Type(c_ClTransferData) :: cData
    integer, intent(in) :: i

    if (i==0) then
        call Convert_ClTransferData(State%CLdata%CTransScal, cData)
    else if (i==1) then
        call Convert_ClTransferData(State%CLdata%CTransVec, cData)
    else if (i==2) then
        call Convert_ClTransferData(State%CLdata%CTransTens, cData)
    else
        error stop 'Unknown ClTransferData index'
    end if

    end subroutine CAMBdata_ClTransferData

    subroutine Convert_ClTransferData(CTrans, cData)
    Type(ClTransferData), target :: CTrans
    Type(c_ClTransferData) :: cData

    cData%NumSources = CTrans%NumSources
    if (allocated(CTrans%q%points)) then
        cData%q_size = size(CTrans%q%points)
        cData%q = c_loc(CTrans%q%points)
    else
        cData%q_size = 0
    end if
    if (associated(CTrans%Delta_p_l_k)) then
        cData%delta_size = shape(CTrans%Delta_p_l_k)
        cData%delta_p_l_k = c_loc(CTrans%Delta_p_l_k)
    else
        cData%delta_size = 0
    end if
    cData%l_size = CTrans%ls%nl
    cData%ls = c_loc(CTrans%ls%l)

    end subroutine Convert_ClTransferData


    subroutine CAMBdata_GetLinearMatterPower(State, PK, var1, var2, hubble_units)
    Type(CAMBstate) :: State
    real(dl) :: PK(State%MT%num_q_trans,State%CP%Transfer%PK_num_redshifts)
    integer, intent(in) :: var1, var2
    logical :: hubble_units

    call Transfer_GetUnsplinedPower(State%MT, State%CP, PK, var1, var2, hubble_units)

    end subroutine CAMBdata_GetLinearMatterPower

    subroutine CAMBdata_GetNonLinearMatterPower(State, PK, var1, var2, hubble_units)
    Type(CAMBstate) :: State
    real(dl) :: PK(State%MT%num_q_trans,State%CP%Transfer%PK_num_redshifts)
    integer, intent(in) :: var1, var2
    logical :: hubble_units

    call Transfer_GetUnsplinedNonlinearPower(State%MT, State%CP,PK, var1, var2, hubble_units)

    end subroutine CAMBdata_GetNonLinearMatterPower


    subroutine CAMBdata_GetMatterPower(State, outpower, minkh, dlnkh, npoints, var1, var2)
    Type(CAMBstate) :: State
    real(dl), intent(out) :: outpower(npoints,State%CP%Transfer%PK_num_redshifts)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in) :: npoints, var1, var2
    integer i

    do i=1,State%CP%Transfer%PK_num_redshifts
        call Transfer_GetMatterPowerD(State%MT, State%CP, outpower(:,i), &
            State%CP%Transfer%PK_num_redshifts-i+1, minkh, dlnkh, npoints, var1, var2)
    end do

    end subroutine CAMBdata_GetMatterPower

    real(dl) function CAMBdata_get_tau_maxvis(State)
    Type(CAMBstate) :: State
    CAMBdata_get_tau_maxvis = State%tau_maxvis
    end function CAMBdata_get_tau_maxvis

    subroutine CAMBParams_setinitialpower(Params, cptr, cls)
    type(CAMBparams) :: Params
    type(c_ptr)  :: cptr
    Type (TInitialPowerLaw), pointer :: pInitialPowerLaw
    Type (TSplinedInitialPower), pointer :: pInitialSplined
    integer, intent(in) :: cls

    if (allocated(Params%InitPower)) deallocate(Params%InitPower)
    if (cls==0) then
        call c_f_pointer(cptr, pInitialPowerLaw)
        allocate(Params%InitPower, source = pInitialPowerLaw)
    elseif (cls==1) then
        call c_f_pointer(cptr, pInitialSplined)
        allocate(Params%InitPower, source = pInitialSplined)
    else
        call MpiStop('Unknown initial power')
    end if

    end subroutine CAMBParams_setinitialpower


    subroutine CAMB_SetTotCls(State,lmax, tot_scalar_Cls)
    type(CAMBstate) State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: tot_scalar_cls(4, 0:lmax)
    integer l

    tot_scalar_cls = 0
    do l=lmin, lmax
        if (State%CP%WantScalars .and. l<= State%CP%Max_l) then
            if (State%CP%DoLensing) then
                if (l<=State%CLData%lmax_lensed) &
                    tot_scalar_cls(1:4,l) = State%CLData%Cl_lensed(l, CT_Temp:CT_Cross)
            else
                tot_scalar_cls(1:2,l) = State%CLData%Cl_scalar(l,C_Temp:C_E)
                tot_scalar_cls(4,l) = State%CLData%Cl_scalar(l, C_Cross)
            endif
        end if
        if (CP%WantTensors .and. l <= CP%Max_l_tensor) then
            tot_scalar_cls(1:4,l) = tot_scalar_cls(1:4,l) &
                + State%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetTotCls

    subroutine CAMB_SetUnlensedCls(State,lmax, unlensed_cls)
    Type(CAMBstate) :: State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: unlensed_cls(4,0:lmax)
    integer l

    unlensed_cls = 0
    do l=lmin, lmax
        if (State%CP%WantScalars .and. l<= CP%Max_l) then
            unlensed_cls(1:2,l) = State%CLData%Cl_scalar(l, C_Temp:C_E)
            unlensed_cls(4,l) = State%CLData%Cl_scalar(l, C_Cross)
        end if
        if (State%CP%WantTensors &
            .and. l <= State%CP%Max_l_tensor) then
            unlensed_cls(1:4,l) = unlensed_cls(1:4,l) &
                + State%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetUnlensedCls

    subroutine CAMB_SetLensPotentialCls(State,lmax, cls)
    use constants
    Type(CAMBstate) :: State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: cls(3, 0:lmax) !phi-phi, phi-T, phi-E
    integer l

    cls = 0
    if (State%CP%WantScalars .and. State%CP%DoLensing) then
        do l=lmin, min(lmax,State%CP%Max_l)
            cls(1,l) = State%CLData%Cl_scalar(l,C_Phi) * (real(l+1)/l)**2/const_twopi
            cls(2:3,l) = State%CLData%Cl_scalar(l,C_PhiTemp:C_PhiE) &
                * ((real(l+1)/l)**1.5/const_twopi)
        end do
    end if

    end subroutine CAMB_SetLensPotentialCls

    subroutine CAMB_SetUnlensedScalCls(State,lmax, scalar_Cls)
    Type(CAMBState) :: State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: scalar_Cls(4, 0:lmax)
    integer lmx

    scalar_Cls = 0
    if (State%CP%WantScalars) then
        lmx = min(State%CP%Max_l, lmax)
        scalar_Cls(1:2,lmin:lmx) = &
            transpose(State%CLData%Cl_Scalar(lmin:lmx, C_Temp:C_E))
        scalar_Cls(4,lmin:lmx) = State%CLData%Cl_Scalar(lmin:lmx, C_Cross)
    end if

    end subroutine CAMB_SetUnlensedScalCls

    subroutine CAMB_SetlensedScalCls(State,lmax, lensed_Cls)
    type(CAMBState) State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: lensed_Cls(4, 0:lmax)
    integer lmx

    lensed_Cls = 0
    if (State%CP%WantScalars .and. State%CP%DoLensing) then
        lmx = min(lmax,State%CLData%lmax_lensed)
        lensed_Cls(1:4,lmin:lmx) = &
            transpose(State%CLData%Cl_lensed(lmin:lmx, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetlensedScalCls

    subroutine CAMB_SetTensorCls(State,lmax, tensor_Cls)
    Type(CAMBstate) :: State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: tensor_Cls(4, 0:lmax)
    integer lmx

    tensor_Cls = 0
    if (State%CP%WantTensors) then
        lmx = min(lmax,State%CP%Max_l_tensor)
        tensor_Cls(1:4,lmin:lmx) = &
            transpose(State%CLData%Cl_Tensor(lmin:lmx, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetTensorCls


    subroutine CAMB_SetUnlensedScalarArray(State,lmax, ScalarArray, n)
    Type(CAMBstate) :: State
    integer, intent(IN) :: lmax, n
    real(dl), intent(OUT) :: ScalarArray(n, n, 0:lmax)
    integer l

    ScalarArray = 0
    if (State%CP%WantScalars) then
        do l=lmin, min(lmax,State%CP%Max_l)
            ScalarArray(1:n,1:n,l) = State%CLData%Cl_scalar_array(l, 1:n,1:n)
        end do
    end if

    end subroutine CAMB_SetUnlensedScalarArray

    subroutine CAMBParams_SetBackgroundOutputs_z(P,redshifts,n)
    type(CAMBparams):: P
    integer, intent(in) :: n
    real(dl), intent(in) :: redshifts(n)

    if (allocated(P%z_outputs)) &
        deallocate(P%z_outputs)
    if (n>0) then
        allocate(P%z_outputs, source=redshifts)
    end if

    end subroutine CAMBParams_SetBackgroundOutputs_z

    function CAMB_GetNumBackgroundOutputs(State)
    Type(CAMBstate) :: State
    integer CAMB_GetNumBackgroundOutputs

    if (.not. allocated(State%CP%z_outputs)) then
        CAMB_GetNumBackgroundOutputs = 0
    else
        CAMB_GetNumBackgroundOutputs = size(State%CP%z_outputs)
    end if

    end function CAMB_GetNumBackgroundOutputs

    subroutine CAMB_GetBackgroundOutputs(State,outputs, n)
    use constants
    Type(CAMBstate) :: State
    integer, intent(in) :: n
    real(dl), intent(out) :: outputs(4,n)
    integer i

    if (allocated(State%CP%z_outputs)) then
        do i=1, size(State%CP%z_outputs)
            outputs(1,i) = State%BackgroundOutputs%rs_by_D_v(i)
            outputs(2,i) = State%BackgroundOutputs%H(i)*c/1e3_dl
            outputs(3,i) = State%BackgroundOutputs%DA(i)
            outputs(4,i) = (1+State%CP%z_outputs(i))* &
                State%BackgroundOutputs%DA(i) * State%BackgroundOutputs%H(i) !F_AP parameter
        end do
    end if

    end subroutine CAMB_GetBackgroundOutputs


    subroutine set_cls_template(cls_template)
    character(len=*), intent(in) :: cls_template

    highL_unlensed_cl_template = trim(cls_template)

    end subroutine set_cls_template

    function CAMB_PrimordialPower(Params, k, powers, n,  i) result(err)
    use constants
    type(CAMBparams) :: Params
    integer, intent(in) :: i,n
    real(dl), intent(in) :: k(n)
    real(dl), intent(out) :: powers(n)
    real(dl) curv
    integer err,ix
    real(dl), external :: GetOmegak

    global_error_flag = 0
    curv =-GetOmegak(Params)/((c/1000)/Params%h0)**2
    call Params%InitPower%Init(curv)
    if (global_error_flag==0) then
        do ix =1, n
            if (i==0) then
                powers(ix) = Params%InitPower%ScalarPower(k(ix))
            elseif (i==2) then
                powers(ix) = Params%InitPower%TensorPower(k(ix))
            else
                error stop 'Unknown power type index'
            end if
            if (global_error_flag /= 0) exit
        end do
    end if
    err= global_error_flag

    end function CAMB_PrimordialPower

    subroutine GetOutputEvolutionFork(State, EV, times, outputs, nsources,ncustomsources)
    use CAMBmain
    implicit none
    type(CAMBState) :: State
    type(EvolutionVars) EV
    real(dl), intent(in) :: times(:)
    real(dl), intent(out) :: outputs(:,:,:)
    integer, intent(in) :: nsources, ncustomsources
    real(dl) tau,tol1,tauend, taustart
    integer j,ind
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar)
    real(dl) yprime(EV%nvar), ddelta, delta, adotoa,dtauda, growth, a
    real(dl), target :: sources(nsources), custom_sources(ncustomsources)
    external dtauda
    real, target :: Arr(Transfer_max)

    w=0
    y=0
    taustart = GetTauStart(min(500._dl,EV%q))
    call initial(EV,y, taustart)

    tau=taustart
    ind=1
    tol1=tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    do j=1,size(times)
        tauend = times(j)
        if (tauend<taustart) cycle

        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
        yprime = 0
        EV%OutputTransfer =>  Arr
        EV%OutputSources => sources
        EV%OutputStep = 0
        if (ncustomsources>0) EV%CustomSources => custom_sources
        call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
        nullify(EV%OutputTransfer, EV%OutputSources, EV%CustomSources)

        a = y(1)
        outputs(1:Transfer_Max, j, EV%q_ix) = Arr
        outputs(Transfer_Max+1, j, EV%q_ix) = a
        outputs(Transfer_Max+2, j, EV%q_ix) = y(2) !etak
        adotoa = 1/(y(1)*dtauda(y(1)))
        ddelta= (yprime(3)*State%grhoc+yprime(4)*State%grhob)/(State%grhob+State%grhoc)
        delta=(State%grhoc*y(3)+State%grhob*y(4))/(State%grhob+State%grhoc)
        growth= ddelta/delta/adotoa
        outputs(Transfer_Max+3, j, EV%q_ix) = adotoa !hubble
        outputs(Transfer_Max+4, j, EV%q_ix) = growth
        if (.not. EV%no_phot_multpoles) then
            outputs(Transfer_Max+5, j, EV%q_ix) = y(EV%g_ix+1) !v_g
            if (EV%TightCoupling) then
                outputs(Transfer_Max+6, j, EV%q_ix) = EV%pig
                outputs(Transfer_Max+7, j, EV%q_ix) = EV%pig/4 !just first order result
            else
                outputs(Transfer_Max+6, j, EV%q_ix) = y(EV%g_ix+2) !pi_g
                outputs(Transfer_Max+7, j, EV%q_ix) = y(EV%polind+2) !E_2
            end if
        end if
        if (.not. EV%no_nu_multpoles) then
            outputs(Transfer_Max+8, j, EV%q_ix) = y(EV%r_ix+1) !v_r
        end if
        outputs(Transfer_max + 9:Transfer_max + 9 + nsources-1, j, EV%q_ix) = sources
        if (ncustomsources > 0) then
            outputs(Transfer_max + 9+nsources: &
                Transfer_max + 9 + nsources + ncustomsources-1, j, EV%q_ix) = custom_sources
        end if

        if (global_error_flag/=0) return
    end do
    end subroutine GetOutputEvolutionFork

    function CAMB_TimeEvolution(this,nq, q, ntimes, times, noutputs, outputs, &
        ncustomsources,c_source_func) result(err)
    use GaugeInterface
    Type(CAMBstate) :: this
    integer, intent(in) :: nq, ntimes, noutputs, ncustomsources
    real(dl), intent(in) :: q(nq), times(ntimes)
    real(dl), intent(out) :: outputs(noutputs, ntimes, nq)
    TYPE(C_FUNPTR), INTENT(IN) :: c_source_func
    integer err
    integer q_ix, i
    Type(EvolutionVars) :: Ev
    procedure(TSource_func), pointer :: old_sources

    call SetActiveState(this)
    if (ncustomsources > 0) then
        ! Convert C to Fortran procedure pointer.
        old_sources => custom_sources_func
        CALL C_F_PROCPOINTER (c_source_func, custom_sources_func)
    end if

    global_error_flag = 0
    outputs = 0
    !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC), PRIVATE(EV, q_ix)
    do q_ix= 1, nq
        if (global_error_flag==0) then
            EV%q_ix = q_ix
            EV%q = q(q_ix)
            EV%TransferOnly=.false.
            EV%q2=EV%q**2
            call GetNumEqns(EV)
            call GetOutputEvolutionFork(State,EV, times, outputs, 3, ncustomsources)
        end if
    end do
    !$OMP END PARALLEL DO
    if (ncustomsources>0) custom_sources_func => old_sources
    err = global_error_flag
    end function CAMB_TimeEvolution

    function GetAllocatableSize() result(sz)
    use iso_c_binding
    Type dummyAllocatable
        class(c_ClTransferData), allocatable :: P
    end Type dummyAllocatable
    Type(dummyAllocatable) :: T
    integer sz

    sz = storage_size(T)

    end function GetAllocatableSize

    subroutine CAMBparams_SetDarkEnergy(P, i, handle)
    Type(CAMBparams), target :: P
    integer, intent(in) :: i
    type(c_ptr), intent(out)  ::  handle

    if (allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
    if (i==0) then
        allocate(TDarkEnergyFluid::P%DarkEnergy)
    else if (i==1) then
        allocate(TDarkEnergyPPF::P%DarkEnergy)
    else
        error stop 'Unknown dark energy model'
    end if

    !can't reference polymorphic type, but can reference first data entry (which is same thing here)
    handle = c_loc(P%DarkEnergy%w_lam)

    end subroutine CAMBparams_SetDarkEnergy

    subroutine CAMBparams_GetAllocatables(P, i,de_handle, nonlin_handle, j, power_handle)
    Type(CAMBparams), target :: P
    integer, intent(out) :: i,j
    type(c_ptr), intent(out)  ::  de_handle, nonlin_handle, power_handle

    if (allocated(P%DarkEnergy)) then
        select type (point => P%DarkEnergy)
        class is (TDarkEnergyFluid)
            i =0
        class is (TDarkEnergyPPF)
            i =1
            class default
            i = -1
        end select
        de_handle = c_loc(P%DarkEnergy%w_lam)
    else
        de_handle = c_null_ptr
    end if
    if (allocated(P%NonLinearModel)) then
        nonlin_handle = c_loc(P%NonLinearModel%Min_kh_nonlinear)
    else
        nonlin_handle = c_null_ptr
    end if
    if (allocated(P%InitPower)) then
        select type (PK => P%InitPower)
        class is (TInitialPowerLaw)
            j=0
        class is (TSplinedInitialPower)
            j=1
            class default
            j=-1
        end select
        power_handle = c_loc(P%InitPower%curv)
    else
        power_handle = c_null_ptr
    end if

    end subroutine CAMBparams_GetAllocatables

    subroutine CAMBparams_SetDarkEnergyTable(DE, a, w, n)
    Type(TDarkEnergyBase) :: DE
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), w(n)

    call DE%SetwTable(a,w)

    end subroutine CAMBparams_SetDarkEnergyTable

    subroutine CAMBparams_DarkEnergyStressEnergy(P, a, grhov_t, w, n)
    Type(CAMBparams) :: P
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n)
    real(dl), intent(out) :: grhov_t(n), w(n)
    real(dl) grhov
    integer i

    call P%DarkEnergy%Init(P%omegav)
    do i=1, n
        call P%DarkEnergy%BackgroundDensityAndPressure(1._dl, a(i), grhov_t(i), w(i))
    end do
    grhov_t = grhov_t/a**2

    end subroutine CAMBparams_DarkEnergyStressEnergy

    subroutine CAMBparams_SetPKTable(P, n, nt, k, PK, PKt, power_handle)
    use interpolation
    integer, intent(in) :: n, nt
    real(dl), intent(in) :: k(max(n,nt)), PK(n), PKt(nt)
    Type(CAMBparams), target :: P
    type(c_ptr), intent(out)  :: power_handle

    if (allocated(P%InitPower)) deallocate(P%InitPower)
    allocate(TSplinedInitialPower::P%InitPower)

    select type (InitPower => P%InitPower)
    class is (TSplinedInitialPower)
        call InitPower%SetScalarTable(n,k, PK)
        call InitPower%SetTensorTable(nt,k, PKt)
    end select
    power_handle = c_loc(P%InitPower%curv)

    end subroutine CAMBparams_SetPKTable

    subroutine CAMBparams_new(handle)
    type(c_ptr), intent(inout) :: handle
    Type (CAMBparams), pointer :: pSource, pData

    call c_f_pointer(handle, pSource)
    if (associated(pSource)) then
        allocate(pData, source=pSource)
    else
        allocate(pData)
        call CAMB_SetDefParams(pData)
    end if
    handle = c_loc(pData)

    end subroutine CAMBparams_new

    subroutine CAMBParams_Free(cptr)
    type(c_ptr)  :: cptr
    Type(CAMBparams), pointer :: P

    call c_f_pointer(cptr, P)
    deallocate(P)

    end subroutine CAMBParams_Free

    subroutine CAMB_SetCustomSourcesFunc(ncustomsources, c_source_func, ell_scales)
    use GaugeInterface
    integer, intent(in) :: ncustomsources
    integer, intent(in) :: ell_scales(ncustomsources)
    TYPE(C_FUNPTR), INTENT(IN) :: c_source_func

    num_custom_sources = ncustomsources
    if (allocated(custom_source_ell_scales)) deallocate(custom_source_ell_scales)
    if (ncustomsources > 0) then
        ! Convert C to Fortran procedure pointer.
        CALL C_F_PROCPOINTER (c_source_func, custom_sources_func)
        allocate(custom_source_ell_scales(num_custom_sources))
        custom_source_ell_scales=ell_scales
    else
        nullify(custom_sources_func)
    end if

    end subroutine CAMB_SetCustomSourcesFunc

    subroutine InitialPowerLaw_new(handle)
    type(c_ptr), intent(inout) :: handle
    Type (TInitialPowerLaw), pointer :: pSource,pInitialPowerLaw

    call c_f_pointer(handle, pSource)
    if (associated(pSource)) then
        allocate(pInitialPowerLaw, source=pSource)
    else
        allocate(pInitialPowerLaw)
    end if
    handle = c_loc(pInitialPowerLaw)

    end subroutine InitialPowerLaw_new

    subroutine InitialPowerLaw_free(cptr)
    type(c_ptr)  :: cptr
    Type (TInitialPowerLaw), pointer :: pInitialPowerLaw

    call c_f_pointer(cptr, pInitialPowerLaw)
    deallocate(pInitialPowerLaw)

    end subroutine InitialPowerLaw_free

    subroutine SplinedInitialPower_new(handle)
    type(c_ptr), intent(inout) :: handle
    Type (TSplinedInitialPower), pointer :: pSource,pSplinedInitialPower

    if (associated(pSource)) then
        allocate(pSplinedInitialPower, source=pSource)
    else
        allocate(pSplinedInitialPower)
    end if
    handle = c_loc(pSplinedInitialPower)

    end subroutine SplinedInitialPower_new

    subroutine SplinedInitialPower_free(cptr)
    type(c_ptr)  :: cptr
    Type (TSplinedInitialPower), pointer :: pSplinedInitialPower

    call c_f_pointer(cptr, pSplinedInitialPower)
    deallocate(pSplinedInitialPower)

    end subroutine SplinedInitialPower_free

    subroutine SplinedInitialPower_setscalartable(this, n, k, PK)
    Type(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    call this%SetScalarTable(n,k,Pk)

    end subroutine SplinedInitialPower_setscalartable

    subroutine SplinedInitialPower_settensortable(this, n, k, PK)
    Type(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    call this%SetTensorTable(n,k,Pk)

    end subroutine SplinedInitialPower_settensortable

    subroutine SplinedInitialPower_setscalarlogregular(this, kmin, kmax, n, PK)
    Type(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    call this%SetScalarLogRegular(kmin, kmax, n, PK)

    end subroutine SplinedInitialPower_setscalarlogregular

    subroutine SplinedInitialPower_settensorlogregular(this, kmin, kmax, n, PK)
    Type(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    call this%SetTensorLogRegular(kmin, kmax, n, PK)

    end subroutine SplinedInitialPower_settensorlogregular


    logical function SplinedInitialPower_HasTensors(this)
    Type(TSplinedInitialPower) :: this

    SplinedInitialPower_HasTensors = allocated(this%Ptensor)

    end function SplinedInitialPower_HasTensors

    subroutine GetBackgroundThermalEvolution(this,ntimes, times, outputs)
    Type(CAMBstate) :: this
    integer, intent(in) :: ntimes
    real(dl), intent(in) :: times(ntimes)
    real(dl) :: outputs(5, ntimes)
    real(dl), allocatable :: spline_data(:), ddxe(:), ddTb(:)
    real(dl) :: d, tau, cs2b, opacity, vis, Tbaryon
    integer i, ix

    associate(T=>this%ThermoData)
        allocate(spline_data(T%nthermo), ddxe(T%nthermo), ddTb(T%nthermo))
        call splini(spline_data,T%nthermo)
        call splder(T%xe,ddxe,T%nthermo,spline_data)
        call splder(T%Tb,ddTb,T%nthermo,spline_data)

        outputs = 0
        do ix = 1, ntimes
            tau = times(ix)
            if (tau < T%tauminn) cycle
            d=log(tau/T%tauminn)/T%dlntau+1._dl
            i=int(d)
            d=d-i
            call T%Values(tau,cs2b, opacity)

            if (i < T%nthermo) then
                outputs(1,ix)=T%xe(i)+d*(ddxe(i)+d*(3._dl*(T%xe(i+1)-T%xe(i)) &
                    -2._dl*ddxe(i)-ddxe(i+1)+d*(ddxe(i)+ddxe(i+1) &
                    +2._dl*(T%xe(i)-T%xe(i+1)))))
                vis=T%emmu(i)+d*(T%demmu(i)+d*(3._dl*(T%emmu(i+1)-T%emmu(i)) &
                    -2._dl*T%demmu(i)-T%demmu(i+1)+d*(T%demmu(i)+T%demmu(i+1) &
                    +2._dl*(T%emmu(i)-T%emmu(i+1)))))
                Tbaryon = T%tb(i)+d*(ddtb(i)+d*(3._dl*(T%tb(i+1)-T%tb(i)) &
                    -2._dl*ddtb(i)-ddtb(i+1)+d*(ddtb(i)+ddtb(i+1) &
                    +2._dl*(T%tb(i)-T%tb(i+1)))))
            else
                outputs(1,ix)=T%xe(T%nthermo)
                vis = T%emmu(T%nthermo)
                Tbaryon = T%Tb(T%nthermo)
            end if

            outputs(2, ix) = opacity
            outputs(3, ix) = opacity*vis
            outputs(4, ix) = cs2b
            outputs(5, ix) = Tbaryon
        end do
    end associate

    end subroutine GetBackgroundThermalEvolution

    subroutine CAMBdata_GetBackgroundDensities(this, n, a_arr, densities)
    ! return array of 8*pi*G*rho*a**4 for each species
    Type(CAMBstate) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a_arr(n)
    real(dl) :: grhov_t, rhonu, grhonu, a
    real(dl), intent(out) :: densities(8,n)
    integer nu_i,i

    do i=1, n
        a = a_arr(i)
        call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)
        grhonu = 0

        if (this%CP%Num_Nu_massive /= 0) then
            !Get massive neutrino density relative to massless
            do nu_i = 1, this%CP%nu_mass_eigenstates
                call ThermalNuBackground%rho(a * this%nu_masses(nu_i), rhonu)
                grhonu = grhonu + rhonu * this%grhormass(nu_i)
            end do
        end if

        densities(2,i) = this%grhok * a**2
        densities(3,i) = this%grhoc * a
        densities(4,i) = this%grhob * a
        densities(5,i) = this%grhog
        densities(6,i) = this%grhornomass
        densities(7,i) = grhonu
        densities(8,i) = grhov_t*a**2
        densities(1,i) = sum(densities(2:8,i))
    end do

    end subroutine CAMBdata_GetBackgroundDensities

    end module handles
