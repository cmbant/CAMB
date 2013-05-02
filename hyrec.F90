    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HyRec
    ! Author: Antony Lewis
    !---------------------------------------------------------------------------------------------------

    module Recombination
    use constants
    use AMLUtils
    implicit none
    private

    type RecombinationParams

    end type RecombinationParams

    character(LEN=*), parameter :: Recombination_Name = 'HyRec'

    public RecombinationParams, Recombination_xe, Recombination_tm, Recombination_init,   &
    Recombination_ReadParams, Recombination_SetDefParams, &
    Recombination_Validate, Recombination_Name


    contains

    subroutine Recombination_ReadParams(R, Ini)
    use IniFile
    Type(RecombinationParams) :: R
    Type(TIniFile) :: Ini


    end subroutine Recombination_ReadParams


    subroutine Recombination_SetDefParams(R)
    type (RecombinationParams) ::R


    end subroutine Recombination_SetDefParams


    subroutine Recombination_Validate(R, OK)
    Type(RecombinationParams), intent(in) :: R
    logical, intent(inout) :: OK


    end subroutine Recombination_Validate



    function Recombination_tm(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_tm, hyrec_tm
    external hyrec_tm

    Recombination_tm =  hyrec_tm(a);

    end function Recombination_tm


    function Recombination_xe(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_xe,hyrec_xe
    external hyrec_xe

    Recombination_xe = hyrec_xe(a);

    end function Recombination_xe


    subroutine Recombination_init(Recomb, OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu)
    use AMLUtils
    implicit none
    Type (RecombinationParams), intent(in) :: Recomb
    real(dl), intent(in) :: OmegaC, OmegaB, OmegaN, OmegaV, h0inp, tcmb, yp, num_nu
    external rec_build_history_camb

    call rec_build_history_camb(OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu)

    end subroutine Recombination_init




    end module Recombination

