module ActiveSources
    use precision
    use interpolation, only: TInterpGrid2D
    use classes
    implicit none

    private
    public :: TActiveSources, TActiveSources_SetCorrelatorTable, TActiveSources_SetActiveEigenmode

    type, extends(TPythonInterfacedClass) :: TActiveSources

        ! For Correlator Table
        real(dl), allocatable    :: k_grid(:)
        real(dl), allocatable    :: tau_grid(:)
        integer                  :: nk = 0
        integer                  :: ntau = 0

        ! Raw Eigenfunction table
        real(dl), allocatable    :: eigenfunctions_table(:,:,:,:)

        ! Arrays of TInterpGrid2D objects for each eigenfunction type
        type(TInterpGrid2D), allocatable :: ef_interp_00(:)
        type(TInterpGrid2D), allocatable :: ef_interp_S(:)
        type(TInterpGrid2D), allocatable :: ef_interp_V(:)
        type(TInterpGrid2D), allocatable :: ef_interp_T(:)

        ! Eigenvalues
        real(dl), allocatable    :: eigenvalues_S(:,:)
        real(dl), allocatable    :: eigenvalues_00(:,:)

        real(dl), allocatable    :: eigenvalues_V(:,:)
        real(dl), allocatable    :: eigenvalues_T(:,:)

        ! For eigenvalue interpolation
        type(TCubicSpline), allocatable :: lambda_interp_S(:)
        type(TCubicSpline), allocatable :: lambda_interp_00(:)

        type(TCubicSpline), allocatable :: lambda_interp_V(:)
        type(TCubicSpline), allocatable :: lambda_interp_T(:)
        logical                         :: eigenvalue_interpolators_set = .false.

        integer                  :: nmodes = 0
        integer                  :: ntypes = 0
        real(dl)                 :: string_mu = 0.0_dl
        real(dl)                 :: weighting = 0.0_dl

        logical                  :: interp_objects_are_set = .false. ! Flag for interpolators

        integer                  :: active_mode_idx = 0 ! 0 means no active source, 1 to N for modes

        !Storage for derivatives
        real(dl), allocatable :: eigenfunc_derivs_table(:,:,:,:) ! Raw table (nk, ntypes, nmodes, ntau)
        type(TInterpGrid2D), allocatable :: ef_deriv_interp_00(:)     ! Interpolators for 00 component derivatives
        type(TInterpGrid2D), allocatable :: ef_deriv_interp_S(:)      ! Interpolators for S component derivatives
        ! We might add V and T derivatives later if needed
        logical :: deriv_interp_objects_are_set = .false. ! Flag for derivative interpolators

    contains
        procedure, nopass :: PythonClass => TActiveSources_PythonClass
        procedure, nopass :: SelfPointer => TActiveSources_SelfPointer
        procedure :: SetCorrelatorTable => TActiveSources_SetCorrelatorTable
        procedure :: SetActiveEigenmode => TActiveSources_SetActiveEigenmode
    end type TActiveSources

contains

    subroutine TActiveSources_SetActiveEigenmode(this, mode_to_set)
        class(TActiveSources), intent(inout) :: this
        integer, intent(in)           :: mode_to_set

        if (mode_to_set >= 0 .and. (mode_to_set == 0 .or. mode_to_set <= this%nmodes) ) then
            this%active_mode_idx = mode_to_set
            if (mode_to_set > 0) then
                print *, "Fortran TActiveSources: Active UETC eigenmode set to: ", this%active_mode_idx
            else
                print *, "Fortran TActiveSources: UETC sources turned OFF (active_mode_idx = 0)."
            endif
        else
            print *, "Fortran TActiveSources Warning: Invalid active_mode_idx requested: ", mode_to_set, &
                     " Max modes available: ", this%nmodes, ". Setting to 0 (OFF)."
            this%active_mode_idx = 0
        endif
    end subroutine TActiveSources_SetActiveEigenmode

    subroutine TActiveSources_SetCorrelatorTable(this, &
        k_grid_in, nk_in, &
        tau_grid_in, ntau_in, &
        eigenfuncs_flat_in, num_eigen_types_in, nmodes_in, &
        evals_S_flat_in, evals_00_flat_in, evals_V_flat_in, evals_T_flat_in, &
        eigenfunc_derivs_logkt_flat_in, &
        mu_in, weighting_param_in)

        class(TActiveSources), intent(inout) :: this
        integer, intent(in)                :: nk_in, ntau_in, num_eigen_types_in, nmodes_in
        real(dl), intent(in)               :: k_grid_in(nk_in), tau_grid_in(ntau_in)
        real(dl), intent(in)               :: eigenfuncs_flat_in(*)
        real(dl), intent(in)               :: eigenfunc_derivs_logkt_flat_in(*)
        real(dl), intent(in)               :: evals_S_flat_in(*),evals_00_flat_in(*), evals_V_flat_in(*), evals_T_flat_in(*)
        real(dl), intent(in)               :: mu_in, weighting_param_in

        integer :: nk_local, ntau_local, nmodes_local, num_eigen_types_local
        integer :: i, j, l, m, mode_idx
        real(dl), allocatable :: slice_2d_for_interp(:,:)
        real(dl), allocatable :: slice_2d_deriv_for_interp(:,:)

        call DeallocateCorrelatorTable(this)
        this%interp_objects_are_set = .false.
        this%eigenvalue_interpolators_set = .false.
        this%deriv_interp_objects_are_set = .false.

        nk_local = nk_in; ntau_local = ntau_in
        nmodes_local = nmodes_in; num_eigen_types_local = num_eigen_types_in

        this%nk = nk_local; this%ntau = ntau_local
        this%ntypes = num_eigen_types_local; this%nmodes = nmodes_local
        this%string_mu = mu_in; this%weighting = weighting_param_in

        print *, "Fortran TActiveSources_SetCorrelatorTable: Received nk=", nk_local, ", ntau=", ntau_local, &
                 ", ntypes=", num_eigen_types_local, ", nmodes=", nmodes_local

        if (nk_local < 2 .or. ntau_local < 2 .or. nmodes_local <= 0 .or. num_eigen_types_local /= 4) then
            print *, "Fortran SetCorrelatorTable: Insufficient dimensions. Aborting setup."
            call EnsureZeroSizeAllocations(this)
            return
        endif

        allocate(this%k_grid(nk_local)); this%k_grid = k_grid_in(1:nk_local)
        allocate(this%tau_grid(ntau_local)); this%tau_grid = tau_grid_in(1:ntau_local)

        ! Store raw eigenfunctions
        allocate(this%eigenfunctions_table(nk_local, num_eigen_types_local, nmodes_local, ntau_local))
        do i = 1, nk_local; do j = 1, num_eigen_types_local; do l = 1, nmodes_local; do m = 1, ntau_local
            this%eigenfunctions_table(i,j,l,m) = eigenfuncs_flat_in( &
                ((( (i-1)*num_eigen_types_local + (j-1) )*nmodes_local + (l-1) )*ntau_local + (m-1) ) + 1 )
        end do; end do; end do; end do

        ! Store raw eigenfunction derivatives
        allocate(this%eigenfunc_derivs_table(nk_local, num_eigen_types_local, nmodes_local, ntau_local))
        do i = 1, nk_local; do j = 1, num_eigen_types_local; do l = 1, nmodes_local; do m = 1, ntau_local
            this%eigenfunc_derivs_table(i,j,l,m) = eigenfunc_derivs_logkt_flat_in( &
                ((( (i-1)*num_eigen_types_local + (j-1) )*nmodes_local + (l-1) )*ntau_local + (m-1) ) + 1 )
        end do; end do; end do; end do

        ! Store raw eigenvalues
        allocate(this%eigenvalues_S(nk_local, nmodes_local))
        allocate(this%eigenvalues_00(nk_local, nmodes_local))

        allocate(this%eigenvalues_V(nk_local, nmodes_local))
        allocate(this%eigenvalues_T(nk_local, nmodes_local))

        ! Reshape flat arrays to original structure
        this%eigenvalues_S = TRANSPOSE(RESHAPE(evals_S_flat_in(1:nk_local*nmodes_local), [nmodes_local, nk_local]))
        this%eigenvalues_00 = TRANSPOSE(RESHAPE(evals_00_flat_in(1:nk_local*nmodes_local), [nmodes_local, nk_local]))
        this%eigenvalues_V = TRANSPOSE(RESHAPE(evals_V_flat_in(1:nk_local*nmodes_local), [nmodes_local, nk_local]))
        this%eigenvalues_T = TRANSPOSE(RESHAPE(evals_T_flat_in(1:nk_local*nmodes_local), [nmodes_local, nk_local]))
        print *, "Fortran SetCorrelatorTable: Raw data tables (functions, derivatives, eigenvalues) stored."

        ! Initialize eigenfunction interpolators
        allocate(this%ef_interp_00(nmodes_local)); allocate(this%ef_interp_S(nmodes_local))
        allocate(this%ef_interp_V(nmodes_local)); allocate(this%ef_interp_T(nmodes_local))
        allocate(slice_2d_for_interp(nk_local, ntau_local))

        do mode_idx = 1, nmodes_local
            do i = 1, nk_local; slice_2d_for_interp(i, :) = this%eigenfunctions_table(i, 1, mode_idx, :); end do
            call this%ef_interp_00(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_for_interp)

            do i = 1, nk_local; slice_2d_for_interp(i, :) = this%eigenfunctions_table(i, 2, mode_idx, :); end do
            call this%ef_interp_S(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_for_interp)

            do i = 1, nk_local; slice_2d_for_interp(i, :) = this%eigenfunctions_table(i, 3, mode_idx, :); end do
            call this%ef_interp_V(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_for_interp)

            do i = 1, nk_local; slice_2d_for_interp(i, :) = this%eigenfunctions_table(i, 4, mode_idx, :); end do
            call this%ef_interp_T(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_for_interp)
        end do
        deallocate(slice_2d_for_interp)
        this%interp_objects_are_set = .true.
        print *, "Fortran SetCorrelatorTable: Initialized TInterpGrid2D for eigenfunctions."

        ! Initialize derivative interpolators
        allocate(this%ef_deriv_interp_00(nmodes_local))
        allocate(this%ef_deriv_interp_S(nmodes_local))
        allocate(slice_2d_deriv_for_interp(nk_local, ntau_local))

        do mode_idx = 1, nmodes_local
            ! Type 1 (index 1 in Fortran for raw table) is '00' component
            do i = 1, nk_local; slice_2d_deriv_for_interp(i, :) = this%eigenfunc_derivs_table(i, 1, mode_idx, :); end do
            call this%ef_deriv_interp_00(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_deriv_for_interp)

            ! Type 2 (index 2 in Fortran for raw table) is 'S' component
            do i = 1, nk_local; slice_2d_deriv_for_interp(i, :) = this%eigenfunc_derivs_table(i, 2, mode_idx, :); end do
            call this%ef_deriv_interp_S(mode_idx)%Init(this%k_grid, this%tau_grid, slice_2d_deriv_for_interp)
        end do
        deallocate(slice_2d_deriv_for_interp)
        this%deriv_interp_objects_are_set = .true.
        print *, "Fortran SetCorrelatorTable: Initialized TInterpGrid2D for eigenfunction derivatives (00, S)."

        ! Initialize eigenvalue interpolators (TCubicSpline)
        allocate(this%lambda_interp_S(nmodes_local))
        allocate(this%lambda_interp_00(nmodes_local))

        allocate(this%lambda_interp_V(nmodes_local))
        allocate(this%lambda_interp_T(nmodes_local))
        do mode_idx = 1, nmodes_local
            call this%lambda_interp_S(mode_idx)%Init(Xarr=this%k_grid, values=this%eigenvalues_S(:, mode_idx), n=nk_local)
            call this%lambda_interp_00(mode_idx)%Init(Xarr=this%k_grid, values=this%eigenvalues_00(:, mode_idx), n=nk_local)
            call this%lambda_interp_V(mode_idx)%Init(Xarr=this%k_grid, values=this%eigenvalues_V(:, mode_idx), n=nk_local)
            call this%lambda_interp_T(mode_idx)%Init(Xarr=this%k_grid, values=this%eigenvalues_T(:, mode_idx), n=nk_local)
        end do
        this%eigenvalue_interpolators_set = .true.
        print *, "Fortran SetCorrelatorTable: Initialized TCubicSpline for eigenvalues."

    contains
                subroutine DeallocateCorrelatorTable(obj_internal)
            class(TActiveSources), intent(inout) :: obj_internal
            integer :: mode_idx_loop

            if (allocated(obj_internal%k_grid)) deallocate(obj_internal%k_grid)
            if (allocated(obj_internal%tau_grid)) deallocate(obj_internal%tau_grid)

            if (allocated(obj_internal%eigenfunctions_table)) deallocate(obj_internal%eigenfunctions_table)
            if (allocated(obj_internal%eigenfunc_derivs_table)) deallocate(obj_internal%eigenfunc_derivs_table)

            if (allocated(obj_internal%eigenvalues_S)) deallocate(obj_internal%eigenvalues_S)
            if (allocated(obj_internal%eigenvalues_00)) deallocate(obj_internal%eigenvalues_00)
            if (allocated(obj_internal%eigenvalues_V)) deallocate(obj_internal%eigenvalues_V)
            if (allocated(obj_internal%eigenvalues_T)) deallocate(obj_internal%eigenvalues_T)

            obj_internal%nk = 0; obj_internal%ntau = 0;
            obj_internal%ntypes = 0; obj_internal%nmodes = 0

            ! Deallocation for eigenfunction TInterpGrid2D arrays
            if (allocated(obj_internal%ef_interp_00)) then
                do mode_idx_loop=1,size(obj_internal%ef_interp_00); call obj_internal%ef_interp_00(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_interp_00)
            endif
            if (allocated(obj_internal%ef_interp_S))  then
                do mode_idx_loop=1,size(obj_internal%ef_interp_S); call obj_internal%ef_interp_S(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_interp_S)
            endif
            if (allocated(obj_internal%ef_interp_V))  then
                do mode_idx_loop=1,size(obj_internal%ef_interp_V); call obj_internal%ef_interp_V(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_interp_V)
            endif
            if (allocated(obj_internal%ef_interp_T))  then
                do mode_idx_loop = 1, size(obj_internal%ef_interp_T); call obj_internal%ef_interp_T(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_interp_T)
            endif
            obj_internal%interp_objects_are_set = .false.

            ! Deallocate eigenfunction derivative TInterpGrid2D arrays
            if (allocated(obj_internal%ef_deriv_interp_00)) then
                do mode_idx_loop=1,size(obj_internal%ef_deriv_interp_00); call obj_internal%ef_deriv_interp_00(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_deriv_interp_00)
            endif
            if (allocated(obj_internal%ef_deriv_interp_S)) then
                do mode_idx_loop=1,size(obj_internal%ef_deriv_interp_S); call obj_internal%ef_deriv_interp_S(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%ef_deriv_interp_S)
            endif
            obj_internal%deriv_interp_objects_are_set = .false.

            ! Deallocate eigenvalue TCubicSpline interpolators
            if (allocated(obj_internal%lambda_interp_S)) then
                do mode_idx_loop=1,size(obj_internal%lambda_interp_S); call obj_internal%lambda_interp_S(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%lambda_interp_S)
            endif
            if (allocated(obj_internal%lambda_interp_V)) then
                do mode_idx_loop=1,size(obj_internal%lambda_interp_V); call obj_internal%lambda_interp_V(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%lambda_interp_V)
            endif
            if (allocated(obj_internal%lambda_interp_T)) then
                do mode_idx_loop=1,size(obj_internal%lambda_interp_T); call obj_internal%lambda_interp_T(mode_idx_loop)%Clear(); enddo
                deallocate(obj_internal%lambda_interp_T)
            endif
            obj_internal%eigenvalue_interpolators_set = .false.

        end subroutine DeallocateCorrelatorTable

        subroutine EnsureZeroSizeAllocations(obj_internal)
             class(TActiveSources), intent(inout) :: obj_internal
             if (.not. allocated(obj_internal%k_grid)) allocate(obj_internal%k_grid(0))
             if (.not. allocated(obj_internal%tau_grid)) allocate(obj_internal%tau_grid(0))

             if (.not. allocated(obj_internal%eigenfunctions_table)) allocate(obj_internal%eigenfunctions_table(0,0,0,0))
             ! For derivatives
             if (.not. allocated(obj_internal%eigenfunc_derivs_table)) allocate(obj_internal%eigenfunc_derivs_table(0,0,0,0))

             if (.not. allocated(obj_internal%eigenvalues_S)) allocate(obj_internal%eigenvalues_S(0,0))
             if (.not. allocated(obj_internal%eigenvalues_V)) allocate(obj_internal%eigenvalues_V(0,0))
             if (.not. allocated(obj_internal%eigenvalues_T)) allocate(obj_internal%eigenvalues_T(0,0))

             if (.not. allocated(obj_internal%ef_interp_00)) allocate(obj_internal%ef_interp_00(0))
             if (.not. allocated(obj_internal%ef_interp_S)) allocate(obj_internal%ef_interp_S(0))
             if (.not. allocated(obj_internal%ef_interp_V)) allocate(obj_internal%ef_interp_V(0))
             if (.not. allocated(obj_internal%ef_interp_T)) allocate(obj_internal%ef_interp_T(0))

             ! !JR NEW: For derivative interpolators
             if (.not. allocated(obj_internal%ef_deriv_interp_00)) allocate(obj_internal%ef_deriv_interp_00(0))
             if (.not. allocated(obj_internal%ef_deriv_interp_S)) allocate(obj_internal%ef_deriv_interp_S(0))

             ! For eigenvalue interpolators
             if (.not. allocated(obj_internal%lambda_interp_S)) allocate(obj_internal%lambda_interp_S(0))
             if (.not. allocated(obj_internal%lambda_interp_V)) allocate(obj_internal%lambda_interp_V(0))
             if (.not. allocated(obj_internal%lambda_interp_T)) allocate(obj_internal%lambda_interp_T(0))
                 end subroutine EnsureZeroSizeAllocations
    end subroutine TActiveSources_SetCorrelatorTable

    function TActiveSources_PythonClass() result(pyClassName)
        character(LEN=:), allocatable :: pyClassName; pyClassName = 'ActiveSources'
    end function TActiveSources_PythonClass
    subroutine TActiveSources_SelfPointer(cptr, P)
        use iso_c_binding; Type(C_PTR) :: cptr
        class(TPythonInterfacedClass), pointer :: P; Type(TActiveSources), pointer :: P_Specific_TActiveSources
        call c_f_pointer(cptr, P_Specific_TActiveSources); P => P_Specific_TActiveSources
    end subroutine TActiveSources_SelfPointer
end module ActiveSources
