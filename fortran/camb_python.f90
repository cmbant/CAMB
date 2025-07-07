
    module handles
    use CAMB
    use Precision
    use results
    use iso_c_binding
    use DarkEnergyFluid
    use DarkEnergyPPF
    use ObjectLists
    use classes
    use Interpolation
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

    Type RealAllocatableArray
        double precision, allocatable :: R(:)
    end Type

    Type IntAllocatableArray
        integer, allocatable :: R(:)
    end Type

    Type PythonClassAllocatableArray
        Type(PythonClassAllocatable), allocatable :: R(:)
    end Type PythonClassAllocatableArray

    abstract interface
    subroutine TSelfPointer(cptr, P)
    use iso_c_binding
    import TPythonInterfacedClass
    Type(c_ptr) :: cptr
    class (TPythonInterfacedClass), pointer :: P
    end subroutine TSelfPointer
    end interface

    integer, private, save, target :: dummy

    contains


    subroutine GetAllocatableSize(sz, sz_array, sz_object_array)
    Type(PythonClassAllocatable) :: T
    Type(RealAllocatableArray) :: T2
    Type(PythonClassAllocatableArray) :: T3
    integer, intent(out) :: sz, sz_array, sz_object_array

    sz = storage_size(T) / 8
    sz_array = storage_size(T2) / 8
    sz_object_array = storage_size(T3) / 8

    end subroutine GetAllocatableSize

    function better_c_loc(R)
    !Longwinded way to get pointer to R (needed in gfortran, can do directly in ifort)
    class(TPythonInterfacedClass), target :: R
    Type(c_ptr) better_c_loc
    type(PythonClassPtr), target :: P
    Type Dum
        integer :: Self
    end type Dum
    Type DummyClassPtr
        class(Dum), pointer :: Ref
    end type
    type(DummyClassPtr), pointer :: mock

    P%Ref => R
    call c_f_pointer(c_loc(P), mock)
    better_c_loc = c_loc(mock%Ref%Self)

    end function better_c_loc

    function get_allocatable_1D_array(array, ptr) result(sz)
    Type(RealAllocatableArray), target :: array
    Type(c_Ptr), intent(out) :: ptr
    integer sz

    if (allocated(array%R)) then
        ptr = c_loc(array%R)
        sz = size(array%R)
    else
        sz=0
    end if

    end function get_allocatable_1D_array

    subroutine set_allocatable_1D_array(array, vals, sz)
    Type(RealAllocatableArray), target :: array
    integer, intent(in) :: sz
    real(dl) :: vals(sz)

    if (allocated(array%R)) deallocate(array%R)
    if (sz>0) then
        allocate(array%R, source = vals)
    end if

    end subroutine set_allocatable_1D_array

    function get_allocatable_1D_array_int(array, ptr) result(sz)
    Type(IntAllocatableArray), target :: array
    Type(c_Ptr), intent(out) :: ptr
    integer sz

    if (allocated(array%R)) then
        ptr = c_loc(array%R)
        sz = size(array%R)
    else
        sz=0
    end if

    end function get_allocatable_1D_array_int

    subroutine set_allocatable_1D_array_int(array, vals, sz)
    Type(IntAllocatableArray), target :: array
    integer, intent(in) :: sz
    integer :: vals(sz)

    if (allocated(array%R)) deallocate(array%R)
    if (sz>0) then
        allocate(array%R, source = vals)
    end if

    end subroutine set_allocatable_1D_array_int

    function get_allocatable_object_1D_array(array, ptr) result(sz)
    Type(PythonClassAllocatableArray), target :: array
    Type(c_Ptr), intent(out) :: ptr
    integer sz

    if (allocated(array%R)) then
        ptr = c_loc(array%R)
        sz = size(array%R)
    else
        sz=0
    end if

    end function get_allocatable_object_1D_array


    subroutine set_allocatable_object_1D_array(array, vals, sz)
    Type(PythonClassAllocatableArray), target :: array
    integer, intent(in) :: sz
    Type(PythonClassAllocatable) :: vals(sz)

    if (allocated(array%R)) deallocate(array%R)
    if (sz>0) then
        allocate(array%R(sz))
        array%R = vals
    end if

    end subroutine set_allocatable_object_1D_array

    function get_effective_null()
    type(c_ptr) get_effective_null

    !Can use actual null pointer in gfortran. ifort requires non-null (but doesn't have to be valid memory)
    get_effective_null = c_loc(dummy)

    end function get_effective_null

    subroutine F2003Class_get_id(SelfPtr, pSource)
    TYPE(C_FUNPTR), INTENT(IN) :: SelfPtr
    procedure(TSelfPointer), pointer :: SelfFunc
    class(TPythonInterfacedClass), pointer :: pSource

    CALL C_F_PROCPOINTER (SelfPtr, SelfFunc)
    call SelfFunc(get_effective_null(), pSource)

    end subroutine F2003Class_get_id

    subroutine F2003Class_GetAllocatable(classobject, id, handle)
    Type(PythonClassAllocatable), target :: classobject
    type(c_ptr), intent(out)  :: handle
    class (TPythonInterfacedClass), pointer :: id

    if (allocated(classobject%P)) then
        call classobject%P%SelfPointer(get_effective_null(), id)
        ! handle =  c_loc(classobject%P) ! only works in ifort
        handle =  better_c_loc(classobject%P)
    else
        handle = c_null_ptr
    end if

    end subroutine F2003Class_GetAllocatable

    subroutine F2003Class_SetAllocatable(f_allocatable, source)
    Type(PythonClassAllocatable), target :: f_allocatable
    class(TPythonInterfacedClass), optional :: source

    if (allocated(f_allocatable%P)) deallocate(f_allocatable%P)
    if (present(Source)) then
        allocate(f_allocatable%P,source = source)
    end if

    end subroutine F2003Class_SetAllocatable

    subroutine F2003Class_new(handle, selfPtr)
    type(c_ptr), intent(inout) :: handle
    TYPE(C_FUNPTR), INTENT(IN) :: SelfPtr
    procedure(TSelfPointer), pointer :: SelfFunc
    class(TPythonInterfacedClass), pointer :: pSource
    class(TPythonInterfacedClass), pointer :: p

    CALL C_F_PROCPOINTER (SelfPtr, SelfFunc)
    if (C_ASSOCIATED(handle)) then
        call SelfFunc(handle, pSource)
        allocate(p, source=pSource)
    else
        call SelfFunc(get_effective_null(), pSource)
        allocate(p, mold = pSource)
    end if
    handle = better_c_loc(p)

    end subroutine F2003Class_new

    subroutine F2003Class_free(cptr,SelfPtr)
    type(c_ptr)  :: cptr
    TYPE(C_FUNPTR), INTENT(IN) :: SelfPtr
    procedure(TSelfPointer), pointer :: SelfFunc
    class(TPythonInterfacedClass), pointer :: pSource

    CALL C_F_PROCPOINTER (SelfPtr, SelfFunc)
    call SelfFunc(cptr, pSource)
    if (.not. associated(pSource)) error stop 'Null in F2003Class_free'
    deallocate(pSource)

    end subroutine F2003Class_free

    function CAMBdata_GetTransfers(Data, Params, onlytransfer, onlytimesources) result(error)
    Type (CAMBdata):: Data
    type(CAMBparams) :: Params
    logical :: onlytransfer, onlytimesources
    integer :: error

    error = 0
    call CAMB_GetResults(Data, Params, error, onlytransfer, onlytimesources)

    end function CAMBdata_GetTransfers

    function CAMBdata_CalcBackgroundTheory(Data, P) result(error)
    use cambmain, only: initvars
    Type (CAMBdata):: Data
    type(CAMBparams) :: P
    integer error

    global_error_flag = 0
    call Data%SetParams(P)
    if (global_error_flag==0) call InitVars(Data) !calculate thermal history, e.g. z_drag etc.
    error=global_error_flag

    end function CAMBdata_CalcBackgroundTheory

    subroutine CAMBdata_GetSigma8(Data, s8, i)
    Type(CAMBdata) :: Data
    integer i
    real(dl) s8(Data%CP%Transfer%PK_num_redshifts)

    if (i==0) then
        s8= Data%MT%sigma_8
    elseif (i==1 .and. allocated(Data%MT%sigma2_vdelta_8)) then
        s8 = Data%MT%sigma2_vdelta_8/Data%MT%sigma_8
    else
        s8 = 0
    end if

    end subroutine CAMBdata_GetSigma8

    subroutine CAMBdata_GetSigmaRArray(Data, sigma, R, nR, z_ix, nz, var1, var2)
    Type(CAMBdata) :: Data
    integer, intent(in) :: nR, nz
    real(dl) :: sigma(nR,nz)
    real(dl) :: R(nR)
    integer var1, var2, z_ix(nz)
    integer i, ix

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(ix)
    do i=1, nz
        ix = z_ix(i)
        if (ix==-1) ix = Data%PK_redshifts_index(Data%CP%Transfer%PK_num_redshifts)
        call Transfer_GetSigmaRArray(Data, Data%MT, R, sigma(:,i), ix, var1, var2)
    end do

    end subroutine CAMBdata_GetSigmaRArray


    subroutine CAMBdata_GetMatterTransferks(Data, nk, ks)
    Type(CAMBdata) :: Data
    integer nk
    real(dl) ks(nk)

    if (nk/=0) then
        if (.not. allocated(Data%MT%TransferData)) then
            global_error_flag = error_allocation
            global_error_message = 'TransferData not allocated: make sure transfer functions computed'
            return
        end if
        ks =  Data%MT%TransferData(Transfer_kh,:,1) * (Data%CP%H0/100)
    end if
    nk = Data%MT%num_q_trans

    end subroutine CAMBdata_GetMatterTransferks

    subroutine CAMBdata_MatterTransferData(Data, cData)
    Type(CAMBdata), target :: Data
    Type(c_MatterTransferData) :: cData

    if (allocated(Data%MT%sigma2_vdelta_8)) then
        cData%sigma2_vdelta_8_size = size(Data%MT%sigma2_vdelta_8)
        cData%sigma2_vdelta_8 = c_loc(Data%MT%sigma2_vdelta_8)
    else
        cData%sigma2_vdelta_8_size =0
    end if
    cData%sigma_8_size = size(Data%MT%sigma_8)
    cData%sigma_8 = c_loc(Data%MT%sigma_8)
    cData%TransferData_size = shape(Data%MT%TransferData)
    cData%num_q_trans = Data%MT%num_q_trans
    cData%q_trans = c_loc(Data%MT%q_trans)
    cData%TransferData = c_loc(Data%MT%TransferData)
    cData%q_trans = c_loc(Data%MT%q_trans)

    end subroutine CAMBdata_MatterTransferData

    subroutine CAMBdata_ClTransferData(Data, cData, i)
    Type(CAMBdata), target :: Data
    Type(c_ClTransferData) :: cData
    integer, intent(in) :: i

    if (i==0) then
        call Convert_ClTransferData(Data%CLdata%CTransScal, cData)
    else if (i==1) then
        call Convert_ClTransferData(Data%CLdata%CTransVec, cData)
    else if (i==2) then
        call Convert_ClTransferData(Data%CLdata%CTransTens, cData)
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
    if (allocated(CTrans%Delta_p_l_k)) then
        cData%delta_size = shape(CTrans%Delta_p_l_k)
        cData%delta_p_l_k = c_loc(CTrans%Delta_p_l_k)
    else
        cData%delta_size = 0
    end if
    cData%l_size = CTrans%ls%nl
    cData%ls = c_loc(CTrans%ls%l)

    end subroutine Convert_ClTransferData


    subroutine CAMBdata_GetLinearMatterPower(Data, PK, var1, var2, hubble_units)
    Type(CAMBdata) :: Data
    real(dl) :: PK(Data%MT%num_q_trans,Data%CP%Transfer%PK_num_redshifts)
    integer, intent(in) :: var1, var2
    logical :: hubble_units

    call Transfer_GetUnsplinedPower(Data, Data%MT, PK, var1, var2, hubble_units)

    end subroutine CAMBdata_GetLinearMatterPower

    subroutine CAMBdata_GetNonLinearMatterPower(Data, PK, var1, var2, hubble_units)
    Type(CAMBdata) :: Data
    real(dl) :: PK(Data%MT%num_q_trans,Data%CP%Transfer%PK_num_redshifts)
    integer, intent(in) :: var1, var2
    logical :: hubble_units

    call Transfer_GetUnsplinedNonlinearPower(Data, Data%MT, PK, var1, var2, hubble_units)

    end subroutine CAMBdata_GetNonLinearMatterPower


    subroutine CAMBdata_GetMatterPower(Data, outpower, minkh, dlnkh, npoints, var1, var2)
    Type(CAMBdata) :: Data
    integer, intent(in) :: npoints, var1, var2
    real(dl), intent(out) :: outpower(npoints,Data%CP%Transfer%PK_num_redshifts)
    real(dl), intent(in) :: minkh, dlnkh
    integer i

    do i=1,Data%CP%Transfer%PK_num_redshifts
        call Transfer_GetMatterPowerD(Data, Data%MT, outpower(:,i), &
            Data%CP%Transfer%PK_num_redshifts-i+1, minkh, dlnkh, npoints, var1, var2)
    end do

    end subroutine CAMBdata_GetMatterPower

    subroutine CAMB_SetTotCls(Data,lmax, tot_scalar_Cls)
    type(CAMBdata) Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: tot_scalar_cls(4, 0:lmax)
    integer l

    tot_scalar_cls = 0
    do l=Data%CP%Min_l, lmax
        if (Data%CP%WantScalars .and. l<= Data%CP%Max_l) then
            if (Data%CP%DoLensing) then
                if (l<=Data%CLData%lmax_lensed) &
                    tot_scalar_cls(1:4,l) = Data%CLData%Cl_lensed(l, CT_Temp:CT_Cross)
            else
                tot_scalar_cls(1:2,l) = Data%CLData%Cl_scalar(l,C_Temp:C_E)
                tot_scalar_cls(4,l) = Data%CLData%Cl_scalar(l, C_Cross)
            endif
        end if
        if (Data%CP%WantTensors .and. l <= Data%CP%Max_l_tensor) then
            tot_scalar_cls(1:4,l) = tot_scalar_cls(1:4,l) &
                + Data%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetTotCls

    subroutine CAMB_SetUnlensedCls(Data,lmax, unlensed_cls)
    Type(CAMBdata) :: Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: unlensed_cls(4,0:lmax)
    integer l

    unlensed_cls = 0
    do l=Data%CP%Min_l, lmax
        if (Data%CP%WantScalars .and. l<= Data%CP%Max_l) then
            unlensed_cls(1:2,l) = Data%CLData%Cl_scalar(l, C_Temp:C_E)
            unlensed_cls(4,l) = Data%CLData%Cl_scalar(l, C_Cross)
        end if
        if (Data%CP%WantTensors &
            .and. l <= Data%CP%Max_l_tensor) then
            unlensed_cls(1:4,l) = unlensed_cls(1:4,l) &
                + Data%CLData%Cl_tensor(l, CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetUnlensedCls

    subroutine CAMB_SetLensPotentialCls(Data,lmax, cls)
    use constants
    Type(CAMBdata) :: Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: cls(3, 0:lmax) !phi-phi, phi-T, phi-E
    integer l

    cls = 0
    if (Data%CP%WantScalars .and. Data%CP%DoLensing) then
        do l=Data%CP%Min_l, min(lmax,Data%CP%Max_l)
            cls(1,l) = Data%CLData%Cl_scalar(l,C_Phi) * (real(l+1)/l)**2/const_twopi
            cls(2:3,l) = Data%CLData%Cl_scalar(l,C_PhiTemp:C_PhiE) &
                * ((real(l+1)/l)**1.5/const_twopi)
        end do
    end if

    end subroutine CAMB_SetLensPotentialCls

    subroutine CAMB_SetUnlensedScalCls(Data,lmax, scalar_Cls)
    Type(CAMBdata) :: Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: scalar_Cls(4, 0:lmax)
    integer lmx

    scalar_Cls = 0
    if (Data%CP%WantScalars) then
        lmx = min(Data%CP%Max_l, lmax)
        scalar_Cls(1:2,Data%CP%Min_l:lmx) = &
            transpose(Data%CLData%Cl_Scalar(Data%CP%Min_l:lmx, C_Temp:C_E))
        scalar_Cls(4,Data%CP%Min_l:lmx) = Data%CLData%Cl_Scalar(Data%CP%Min_l:lmx, C_Cross)
    end if

    end subroutine CAMB_SetUnlensedScalCls

    subroutine CAMB_SetlensedScalCls(Data,lmax, lensed_Cls)
    type(CAMBdata) Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: lensed_Cls(4, 0:lmax)
    integer lmx

    lensed_Cls = 0
    if (Data%CP%WantScalars .and. Data%CP%DoLensing) then
        lmx = min(lmax,Data%CLData%lmax_lensed)
        lensed_Cls(1:4,Data%CP%Min_l:lmx) = &
            transpose(Data%CLData%Cl_lensed(Data%CP%Min_l:lmx, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetlensedScalCls

    subroutine CAMB_SetTensorCls(Data,lmax, tensor_Cls)
    Type(CAMBdata) :: Data
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: tensor_Cls(4, 0:lmax)
    integer lmx

    tensor_Cls = 0
    if (Data%CP%WantTensors) then
        lmx = min(lmax,Data%CP%Max_l_tensor)
        tensor_Cls(1:4,Data%CP%Min_l:lmx) = &
            transpose(Data%CLData%Cl_Tensor(Data%CP%Min_l:lmx, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetTensorCls
    
    subroutine CAMB_SetVectorCls(State,lmax, vector_Cls)
    Type(CAMBdata) :: State
    integer, intent(IN) :: lmax
    real(dl), intent(OUT) :: vector_Cls(4, 0:lmax)
    integer lmx

    vector_Cls = 0
    if (State%CP%WantVectors) then
        lmx = min(lmax,State%CP%Max_l)
        vector_Cls(1:4,State%CP%Min_l:lmx) = &
            transpose(State%CLData%Cl_Vector(State%CP%Min_l:lmx, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetVectorCls


    subroutine CAMB_SetUnlensedScalarArray(Data,lmax, ScalarArray, n)
    Type(CAMBdata) :: Data
    integer, intent(IN) :: lmax, n
    real(dl), intent(OUT) :: ScalarArray(n, n, 0:lmax)
    integer l

    ScalarArray = 0
    if (Data%CP%WantScalars) then
        do l=Data%CP%Min_l, min(lmax,Data%CP%Max_l)
            ScalarArray(1:n,1:n,l) = Data%CLData%Cl_scalar_array(l, 1:n,1:n)
        end do
    end if

    end subroutine CAMB_SetUnlensedScalarArray

    subroutine CAMB_GetBackgroundOutputs(Data,outputs, n)
    use constants
    Type(CAMBdata) :: Data
    integer, intent(in) :: n
    real(dl), intent(out) :: outputs(4,n)
    integer i

    if (allocated(Data%CP%z_outputs)) then
        do i=1, size(Data%CP%z_outputs)
            outputs(1,i) = Data%BackgroundOutputs%rs_by_D_v(i)
            outputs(2,i) = Data%BackgroundOutputs%H(i)*c/1e3_dl
            outputs(3,i) = Data%BackgroundOutputs%DA(i)
            outputs(4,i) = (1+Data%CP%z_outputs(i))* &
                Data%BackgroundOutputs%DA(i) * Data%BackgroundOutputs%H(i) !F_AP parameter
        end do
    end if

    end subroutine CAMB_GetBackgroundOutputs


    subroutine set_cls_template(cls_template)
    character(len=*), intent(in) :: cls_template

    if (allocated(highL_CL_template)) deallocate(highL_CL_template)
    highL_unlensed_cl_template = trim(cls_template)
    call CheckLoadedHighLTemplate

    end subroutine set_cls_template

    subroutine GetOutputEvolutionFork(Data, EV, times, outputs, nsources,ncustomsources)
    use CAMBmain
    type(CAMBdata) :: Data
    type(EvolutionVars) EV
    real(dl), intent(in) :: times(:)
    real(dl), intent(out) :: outputs(:,:,:)
    integer, intent(in) :: nsources, ncustomsources
    real(dl) tau,tol1,tauend, taustart
    integer j,ind
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), cs2, opacity
    real(dl) yprime(EV%nvar), ddelta, delta, adotoa,growth, a
    real(dl), target :: sources(nsources), custom_sources(ncustomsources)
    real, target :: Arr(Transfer_max)
    procedure(obj_function) :: dtauda

    w=0
    y=0
    taustart = GetTauStart(min(500._dl,EV%q))
    call initial(EV,y, taustart)

    tau=taustart
    ind=1
    tol1=base_tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
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
        call Data%ThermoData%Values(tau,a, cs2,opacity)
        outputs(1:Transfer_Max, j, EV%q_ix) = Arr
        outputs(Transfer_Max+1, j, EV%q_ix) = a
        outputs(Transfer_Max+2, j, EV%q_ix) = y(ix_etak) !etak
        adotoa = 1/(a*dtauda(Data,a))
        ddelta= (yprime(ix_clxc)*Data%grhoc+yprime(ix_clxb)*Data%grhob)/(Data%grhob+Data%grhoc)
        delta=(Data%grhoc*y(ix_clxc)+Data%grhob*y(ix_clxb))/(Data%grhob+Data%grhoc)
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

    function CAMB_TimeEvolution(this, nq, q, ntimes, times, noutputs, outputs, &
        ncustomsources,c_source_func) result(err)
    use GaugeInterface
    use CAMBmain
    Type(CAMBdata),target :: this
    integer, intent(in) :: nq, ntimes, noutputs, ncustomsources
    real(dl), intent(in) :: q(nq), times(ntimes)
    real(dl), intent(out) :: outputs(noutputs, ntimes, nq)
    TYPE(C_FUNPTR), INTENT(IN) :: c_source_func
    integer err, q_ix
    real(dl) taustart
    Type(EvolutionVars) :: Ev
    Type(TCUstomSourceParams) :: Old

    call SetActiveState(this)
    if (ncustomsources > 0) then
        ! Convert C to Fortran procedure pointer.
        Old = State%CP%CustomSources
        State%CP%CustomSources%c_source_func = c_source_func
        State%CP%CustomSources%num_custom_sources = ncustomsources
    end if

    global_error_flag = 0
    outputs = 0
    taustart = min(times(1),GetTauStart(maxval(q)))
    if (.not. this%ThermoData%HasTHermoData .or. taustart < this%ThermoData%tauminn) call this%ThermoData%Init(this,taustart)
    !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC), PRIVATE(EV, q_ix)
    do q_ix= 1, nq
        if (global_error_flag==0) then
            EV%q_ix = q_ix
            EV%q = q(q_ix)
            EV%TransferOnly=.false.
            EV%q2=EV%q**2
            EV%ThermoData => this%ThermoData
            call GetNumEqns(EV)
            call GetOutputEvolutionFork(State,EV, times, outputs, 3, ncustomsources)
        end if
    end do
    !$OMP END PARALLEL DO
    if (ncustomsources>0) State%CP%CustomSources = Old
    err = global_error_flag

    end function CAMB_TimeEvolution


    subroutine GetBackgroundThermalEvolution(this, ntimes, times, outputs)
    use splines
    Type(CAMBdata) :: this
    integer, intent(in) :: ntimes
    real(dl), intent(in) :: times(ntimes)
    real(dl) :: outputs(9, ntimes)
    real(dl), allocatable :: spline_data(:), ddxe(:), ddTb(:)
    real(dl) :: a, d, tau, cs2b, opacity, Tbaryon, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow
    integer i, ix

    if (.not. this%ThermoData%HasTHermoData) call this%ThermoData%Init(this,min(1d-3,max(1d-5,minval(times))))

    associate(T=>this%ThermoData)
        allocate(spline_data(T%nthermo), ddxe(T%nthermo), ddTb(T%nthermo))
        call splini(spline_data,T%nthermo)
        call splder(T%xe,ddxe,T%nthermo,spline_data)
        call splder(T%Tb,ddTb,T%nthermo,spline_data)

        outputs = 0
        do ix = 1, ntimes
            tau = times(ix)
            if (tau < T%tauminn*1.01) cycle
            d=log(tau/T%tauminn)/T%dlntau+1._dl
            i=int(d)
            d=d-i
            call T%Values(tau,a,cs2b, opacity)
            call T%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
                visibility, dvisibility, ddvisibility, exptau, lenswindow)

            if (i < T%nthermo) then
                outputs(1,ix)=T%xe(i)+d*(ddxe(i)+d*(3._dl*(T%xe(i+1)-T%xe(i)) &
                    -2._dl*ddxe(i)-ddxe(i+1)+d*(ddxe(i)+ddxe(i+1) &
                    +2._dl*(T%xe(i)-T%xe(i+1)))))
                Tbaryon = T%tb(i)+d*(ddtb(i)+d*(3._dl*(T%tb(i+1)-T%tb(i)) &
                    -2._dl*ddtb(i)-ddtb(i+1)+d*(ddtb(i)+ddtb(i+1) &
                    +2._dl*(T%tb(i)-T%tb(i+1)))))
            else
                outputs(1,ix)=T%xe(T%nthermo)
                Tbaryon = T%Tb(T%nthermo)
            end if

            outputs(2, ix) = opacity
            outputs(3, ix) = visibility
            outputs(4, ix) = cs2b
            outputs(5, ix) = Tbaryon
            outputs(6, ix) = dopacity
            outputs(7, ix) = ddopacity
            outputs(8, ix) = dvisibility
            outputs(9, ix) = ddvisibility
        end do
    end associate

    end subroutine GetBackgroundThermalEvolution

    end module handles
