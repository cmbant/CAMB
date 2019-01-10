    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HyRec
    ! Note you will need to rename dtauda_ in history.c to exported_dtauda.
    ! To use with the python wrapper add -fPIC to the HYREC CCFLAGS (for gcc)
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

    type, extends(TRecombinationModel) :: THyRec
    contains
    procedure :: Init => THyRec_init
    procedure :: x_e => THyRec_xe
    procedure :: T_m => THyRec_tm !baryon temperature
    procedure, nopass :: SelfPointer => THyRec_SelfPointer
    end type THyRec

    class(CAMBdata), pointer :: CurrentState

    public THyRec
    contains


    function Thyrec_tm(this,a)
    class(THyRec) :: this
    real(dl), intent(in) :: a
    real(dl) Thyrec_tm,hyrec_tm
    external hyrec_tm

    Thyrec_tm =  hyrec_tm(a)

    end function Thyrec_tm

    function THyRec_xe(this,a)
    class(THyRec) :: this
    real(dl), intent(in) :: a
    real(dl) THyRec_xe,hyrec_xe
    external hyrec_xe

    THyRec_xe = hyrec_xe(a);

    end function THyRec_xe

    real(dl) function THyrec_dtauda(a) BIND(C, NAME='exported_dtauda')
    real(dl), intent(in) :: a
    procedure(obj_function) :: dtauda

    THyrec_dtauda = dtauda(CurrentState,a)
    end function THyrec_dtauda

    subroutine THyRec_init(this,State, WantTSpin)
    class(THyRec), target :: this
    class(TCAMBdata), target :: State
    logical, intent(in), optional :: WantTSpin
    real(dl) OmegaB, OmegaC, OmegaN, h2
    external rec_build_history_camb

    if (DefaultFalse(WantTSpin)) call MpiStop('HyRec does not support 21cm')

    select type(State)
    class is (CAMBdata)
        CurrentState => State

        if (State%CP%Evolve_delta_xe) &
            call MpiStop('HyRec currently does not support evolving Delta x_e')

        h2 = (State%CP%H0/100)**2
        OmegaB = State%CP%ombh2/h2
        OmegaN = State%CP%omnuh2/h2
        OmegaC = State%CP%omch2/h2

        call rec_build_history_camb(OmegaC, OmegaB, OmegaN, State%Omega_de, &
            State%CP%H0, State%CP%tcmb, State%CP%Yhe, State%CP%N_eff())

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

