    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HYREC-2 (July 2020)
    ! HYREC-2 is available at https://github.com/nanoomlee/HYREC-2 
    ! Download HYREC-2 and place it in the parent directory of CAMB
    ! If the model defined in history.h of HYREC-2 is not the default one, SWIFT,
    ! then Nz should be changed (See line 21 and 22)
    !---------------------------------------------------------------------------------------------------

    module HyRec
    use precision
    use constants
    use classes
    use MathUtils
    use results
    use config
    use MpiUtils, only : MpiStop
    implicit none
    private
    
	real(dl), parameter ::  zinitial = 8e3_dl !highest redshift
    real(dl), parameter ::  zfinal=0._dl
    integer, parameter :: Nz=2248            !For SWIFT model of HYREC-2
    !integer,  parameter :: Nz=105859          !For the rest of models
    
    Type RecombinationData
        real(dl), private :: xhyrec(Nz), tmhyrec(Nz)
        real(dl), private :: Tnow
    end Type RecombinationData

    type, extends(TRecombinationModel) :: THyRec
        Type(RecombinationData), allocatable :: Calc
    contains
    procedure :: Init => THyRec_init
    procedure :: x_e => THyRec_xe
    procedure :: T_m => THyRec_tm !baryon temperature
    procedure :: xe_Tm => THyRec_xe_Tm
    procedure, nopass :: SelfPointer => THyRec_SelfPointer
    end type THyRec

    class(CAMBdata), pointer :: CurrentState

    public THyRec
    contains


    function THyRec_tm(this,a)
    class(THyRec) :: this
    real(dl), intent(in) :: a
    real(dl) THyRec_tm,hyrec_tm
	real(dl) z
    external hyrec_tm

    z=1/a-1
    associate( Calc => this%Calc)
        if (z >= 8000) then
            THyRec_tm=Calc%Tnow/a
        else
            if (z <= 0) then
                THyRec_tm=Calc%tmhyrec(nz)
            else
                THyRec_tm=hyrec_tm(a,Calc%tmhyrec)
            endif
        endif
    end associate

    end function THyRec_tm


    function THyRec_xe(this,a)
    class(THyRec) :: this
    real(dl), intent(in) :: a
    real(dl) THyRec_xe,hyrec_xe
	real(dl) z
    external hyrec_xe
    
    z=1/a-1
    associate( Calc => this%Calc)
        if (z >= zinitial) then
            THyRec_xe=Calc%xhyrec(1)
        else
            if (z <= zfinal) then
                THyRec_xe=Calc%xhyrec(nz)
            else
                THyRec_xe=hyrec_xe(a,Calc%xhyrec)
            endif
        endif
    end associate

    end function THyRec_xe

    subroutine THyRec_xe_Tm(this,a, xe, Tm)
    class(THyRec) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: xe, Tm
    real(dl) hyrec_xe, hyrec_tm
    real(dl) z
    external hyrec_xe, hyrec_tm

    z=1/a-1
    associate(Calc => this%Calc)
        if (z >= zinitial) then
            xe=Calc%xhyrec(1)
            Tm=Calc%Tnow/a
        else
            if (z <= zfinal) then
                xe=Calc%xhyrec(nz)
                Tm=Calc%tmhyrec(nz)
            else
                xe=hyrec_xe(a,Calc%xhyrec)
                Tm=hyrec_tm(a,Calc%tmhyrec)
            endif
        endif
    
    end associate
    end subroutine THyRec_xe_Tm


    real(dl) function THyrec_dtauda(a) BIND(C, NAME='exported_dtauda')
    real(dl), intent(in) :: a
    procedure(obj_function) :: dtauda

    THyrec_dtauda = dtauda(CurrentState,a)
    end function THyrec_dtauda

    subroutine THyRec_init(this,State, WantTSpin)
    class(THyRec), target :: this
    class(TCAMBdata), target :: State
    Type(RecombinationData), pointer :: Calc
	logical, intent(in), optional :: WantTSpin
    real(dl) OmegaB, OmegaC, OmegaN, h2
    external rec_build_history_camb

    if (DefaultFalse(WantTSpin)) call MpiStop('HyRec does not support 21cm')

    if (.not. allocated(this%Calc)) allocate(this%Calc)
    Calc => this%Calc

    select type(State)
    class is (CAMBdata)
        CurrentState => State

        if (State%CP%Evolve_delta_xe) &
            call MpiStop('HyRec currently does not support evolving Delta x_e')

        h2 = (State%CP%H0/100)**2
        OmegaB = State%CP%ombh2/h2
        OmegaC = State%CP%omch2/h2
        
        Calc%Tnow=State%CP%tcmb

        call rec_build_history_camb(OmegaC, OmegaB, &
            State%CP%H0, State%CP%tcmb, State%CP%Yhe, State%CP%N_eff(), &
            this%Calc%xhyrec, this%Calc%tmhyrec, Nz)
    end select
    end subroutine THyRec_init

    subroutine THyRec_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (THyRec), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine THyRec_SelfPointer

    end module HyRec

