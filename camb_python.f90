    module handles
    use CAMB
    use Precision
    use ModelParams
    use Transfer
    use iso_c_binding
    implicit none

    Type c_MatterTransferData
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        type(c_ptr) :: q_trans, sigma_8, sigma2_vdelta_8, TransferData
        integer :: sigma_8_size(2)
        integer :: sigma2_vdelta_8_size(2)
        integer :: TransferData_size(3)
    end Type c_MatterTransferData

    contains

    !SPECIAL BRIDGE ROUTINES FOR PYTHON

    subroutine CAMBdata_new(handle)
    type(c_ptr), intent(out) :: handle
    Type (CAMBdata), pointer :: pCAMBdata

    allocate(pCAMBdata)
    call CAMB_InitCAMBdata(pCAMBdata)
    call CAMB_setDefParams(pCAMBdata%Params)
    handle = c_loc(pCAMBdata)

    end subroutine CAMBdata_new

    subroutine CAMBdata_free(cptr)
    type(c_ptr)  :: cptr
    Type (CAMBdata), pointer :: pCAMBdata

    call c_f_pointer(cptr, pCAMBdata)
    call CAMB_FreeCAMBdata(pCAMBdata)
    deallocate(pCAMBdata)

    end subroutine CAMBdata_free

    subroutine CAMBdata_setParams(data, Params)
    Type (CAMBdata), target :: data
    type(CAMBparams) :: Params

    data%Params = Params

    end subroutine CAMBdata_setParams

    subroutine CAMBdata_getParams(data, handle) 
    Type (CAMBdata), target :: data
    type(c_ptr), intent(out)  ::  handle

    handle = c_loc(data%Params)

    end subroutine CAMBdata_getParams

    function CAMBdata_GetTransfers(data, Params, onlytransfer) result(error)
    Type (CAMBdata):: data
    type(CAMBparams) :: Params, P
    logical(kind=c_bool)  :: onlytransfer
    integer :: error

    P = Params
    if (P%DoLensing .and. (P%NonLinear == NonLinear_lens .or. P%NonLinear == NonLinear_both)) then
        P%WantTransfer = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if
    call Transfer_SortAndIndexRedshifts(P%Transfer)
    error = 0

    P%OnlyTransfers = onlytransfer
    call CAMB_GetTransfers(P, data, error)

    end function CAMBdata_GetTransfers

    subroutine CAMBdata_SetParamsForBackground(data, P)
    Type (CAMBdata):: data
    type(CAMBparams) :: P

    global_error_flag = 0
    data%Params = P
    call CAMBParams_Set(data%Params)
    end subroutine CAMBdata_SetParamsForBackground

    function CAMBdata_CalcBackgroundTheory(data, P) result(error)
    use cambmain, only: initvars
    Type (CAMBdata):: data
    type(CAMBparams) :: P
    integer error

    global_error_flag = 0
    data%Params = P
    call CAMBParams_Set(data%Params)
    call InitVars !calculate thermal history, e.g. z_drag etc.
    error=global_error_flag

    end function CAMBdata_CalcBackgroundTheory


    subroutine CAMBdata_MatterTransferData(data, cData)
    Type(CAMBdata) :: data
    Type(c_MatterTransferData) :: cData

    cData%num_q_trans = data%MTrans%num_q_trans
    cData%q_trans = c_loc(data%MTrans%q_trans)
    cData%sigma_8 = c_loc(data%MTrans%sigma_8)
    cData%sigma2_vdelta_8 = c_loc(data%MTrans%sigma2_vdelta_8)
    cData%TransferData = c_loc(data%MTrans%TransferData)
    cData%q_trans = c_loc(data%MTrans%q_trans)
    cData%sigma_8_size(1) = size(data%MTrans%sigma_8,1)
    cData%sigma_8_size(2) = size(data%MTrans%sigma_8,2)
    cData%sigma2_vdelta_8_size(1) = size(data%MTrans%sigma2_vdelta_8,1)
    cData%sigma2_vdelta_8_size(2) = size(data%MTrans%sigma2_vdelta_8,2)
    cData%TransferData_size(1) = size(data%MTrans%TransferData,1)
    cData%TransferData_size(2) = size(data%MTrans%TransferData,2)
    cData%TransferData_size(3) = size(data%MTrans%TransferData,3)

    end subroutine CAMBdata_MatterTransferData

    subroutine CAMBdata_GetLinearMatterPower(data, PK, var1, var2, hubble_units)
    Type(CAMBdata) :: data
    real(dl) :: PK(data%MTrans%num_q_trans,data%Params%Transfer%PK_num_redshifts)
    integer, intent(in) :: var1, var2
    logical :: hubble_units

    call Transfer_GetUnsplinedPower(data%MTrans, PK, var1, var2, hubble_units)

    end subroutine CAMBdata_GetLinearMatterPower


    subroutine CAMBdata_GetMatterPower(data, outpower, minkh, dlnkh, npoints, var1, var2)
    Type(CAMBdata) :: data
    real(dl), intent(out) :: outpower(npoints,data%Params%Transfer%PK_num_redshifts)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in) :: npoints, var1, var2
    integer i

    do i=1,data%Params%Transfer%PK_num_redshifts
        call Transfer_GetMatterPowerD(data%MTrans,outpower(:,i), data%Params%Transfer%PK_num_redshifts-i+1, &
            & 1, minkh, dlnkh, npoints, var1, var2)
    end do

    end subroutine CAMBdata_GetMatterPower


    subroutine CAMB_setinitialpower(Params, P)
    type(CAMBparams) :: Params
    type(InitialPowerParams) :: P

    Params%InitPower = P

    end subroutine CAMB_setinitialpower


    subroutine CAMB_SetTotCls(lmax, tot_scalar_Cls, in)
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: tot_scalar_cls(4, 0:lmax)
    integer l

    tot_scalar_cls = 0
    do l=lmin, lmax
        if (CP%WantScalars .and. l<= CP%Max_l) then
            if (CP%DoLensing) then
                if (l<=lmax_lensed) tot_scalar_cls(1:4,l) = Cl_lensed(l, in, CT_Temp:CT_Cross)
            else
                tot_scalar_cls(1:2,l) = Cl_scalar(l, in,  C_Temp:C_E)
                tot_scalar_cls(4,l) = Cl_scalar(l, in,  C_Cross)
            endif
        end if
        if (CP%WantTensors .and. l <= CP%Max_l_tensor) then
            tot_scalar_cls(1:4,l) = tot_scalar_cls(1:4,l) + Cl_tensor(l, in,  CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetTotCls

    subroutine CAMB_SetUnlensedCls(lmax, unlensed_cls, in)
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: unlensed_cls(4,0:lmax)
    integer l

    unlensed_cls = 0
    do l=lmin, lmax
        if (CP%WantScalars .and. l<= CP%Max_l) then
            unlensed_cls(1:2,l) = Cl_scalar(l, in,  C_Temp:C_E)
            unlensed_cls(4,l) = Cl_scalar(l, in,  C_Cross)
        end if
        if (CP%WantTensors .and. l <= CP%Max_l_tensor) then
            unlensed_cls(1:4,l) = unlensed_cls(1:4,l) + Cl_tensor(l, in,  CT_Temp:CT_Cross)
        end if
    end do

    end subroutine CAMB_SetUnlensedCls

    subroutine CAMB_SetLensPotentialCls(lmax, cls, in)
    use constants
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: cls(3, 0:lmax) !phi-phi, phi-T, phi-E
    integer l

    cls = 0
    if (CP%WantScalars .and. CP%DoLensing) then
        do l=lmin, min(lmax,CP%Max_l)
            cls(1,l) = Cl_scalar(l,in,C_Phi) * (real(l+1)/l)**2/const_twopi
            cls(2:3,l) = Cl_scalar(l,in,C_PhiTemp:C_PhiE) * ((real(l+1)/l)**1.5/const_twopi)
        end do
    end if

    end subroutine CAMB_SetLensPotentialCls

    subroutine CAMB_SetUnlensedScalCls(lmax, scalar_Cls, in)
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: scalar_Cls(4, 0:lmax)
    integer lmx

    scalar_Cls = 0
    if (CP%WantScalars) then
        lmx = min(CP%Max_l, lmax)
        scalar_Cls(1:2,lmin:lmx) = transpose(Cl_Scalar(lmin:lmx, in,C_Temp:C_E))
        scalar_Cls(4,lmin:lmx) = Cl_Scalar(lmin:lmx, in,C_Cross)
    end if

    end subroutine CAMB_SetUnlensedScalCls

    subroutine CAMB_SetlensedScalCls(lmax, lensed_Cls, in)
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: lensed_Cls(4, 0:lmax)
    integer lmx

    lensed_Cls = 0
    if (CP%WantScalars .and. CP%DoLensing) then
        lmx = min(lmax,lmax_lensed)
        lensed_Cls(1:4,lmin:lmx) = transpose(Cl_lensed(lmin:lmx, in,CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetlensedScalCls

    subroutine CAMB_SetTensorCls(lmax, tensor_Cls, in)
    integer, intent(IN) :: lmax, in
    real(dl), intent(OUT) :: tensor_Cls(4, 0:lmax)
    integer lmx

    tensor_Cls = 0
    if (CP%WantTensors) then
        lmx = min(lmax,CP%Max_l_tensor)
        tensor_Cls(1:3,lmin:lmx) = transpose(Cl_Tensor(lmin:lmx, in, CT_Temp:CT_Cross))
    end if

    end subroutine CAMB_SetTensorCls


    subroutine CAMB_SetUnlensedScalarArray(lmax, ScalarArray, in, n)
    integer, intent(IN) :: lmax, in, n
    real(dl), intent(OUT) :: ScalarArray(n, n, 0:lmax)
    integer l

    ScalarArray = 0
    if (CP%WantScalars) then
        do l=lmin, min(lmax,CP%Max_l)
            ScalarArray(1:n,1:n,l) = Cl_scalar_array(l, in, 1:n,1:n)
        end do
    end if

    end subroutine CAMB_SetUnlensedScalarArray

    subroutine CAMB_SetBackgroundOutputs_z(redshifts,n)
    integer, intent(in) :: n
    real(dl), intent(in) :: redshifts(n)

    if (associated(BackgroundOutputs%z_outputs)) deallocate(BackgroundOutputs%z_outputs)
    if (n>0) then
        allocate(BackgroundOutputs%z_outputs(n))
        BackgroundOutputs%z_outputs = redshifts
    else
        nullify(BackgroundOutputs%z_outputs)
    end if

    end subroutine CAMB_SetBackgroundOutputs_z

    function CAMB_GetNumBackgroundOutputs()
    integer CAMB_GetNumBackgroundOutputs

    if (.not. associated(BackgroundOutputs%z_outputs)) then
        CAMB_GetNumBackgroundOutputs = 0
    else
        CAMB_GetNumBackgroundOutputs = size(BackgroundOutputs%z_outputs)
    end if

    end function CAMB_GetNumBackgroundOutputs

    subroutine CAMB_GetBackgroundOutputs(outputs, n)
    use constants
    integer, intent(in) :: n
    real(dl), intent(out) :: outputs(4,n)
    integer i

    if (associated(BackgroundOutputs%z_outputs)) then
        do i=1, size(BackgroundOutputs%z_outputs)
            outputs(1,i) = BackgroundOutputs%rs_by_D_v(i)
            outputs(2,i) = BackgroundOutputs%H(i)*c/1e3_dl
            outputs(3,i) = BackgroundOutputs%DA(i)
            outputs(4,i) = (1+BackgroundOutputs%z_outputs(i))* &
                BackgroundOutputs%DA(i) * BackgroundOutputs%H(i) !F_AP parameter
        end do
    end if

    end subroutine CAMB_GetBackgroundOutputs


    subroutine set_cls_template(cls_template)
    character(len=*), intent(in) :: cls_template

    highL_unlensed_cl_template = trim(cls_template)

    end subroutine set_cls_template


    ! END BRIDGE FOR PYTHON


    end module handles
