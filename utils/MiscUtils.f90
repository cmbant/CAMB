    module MiscUtils
    implicit none

    INTERFACE PresentDefault
    module procedure PresentDefault_S, PresentDefault_L, PresentDefault_I, PresentDefault_R, PresentDefault_D
    END INTERFACE PresentDefault

    INTERFACE IfThenElse
    module procedure IfThenElse_S, IfThenElse_L, IfThenElse_I, IfThenElse_R, IfThenElse_D
    END INTERFACE IfThenElse

    INTERFACE IsFloat
    module procedure IsFloat0, IsFloat1, IsFloat2
    END INTERFACE IsFloat

    contains

    function DefaultFalse(S) result(Sout)
    logical, intent(in), optional :: S
    logical :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = .false.
    end if

    end function DefaultFalse

    function DefaultTrue(S) result(Sout)
    logical, intent(in), optional :: S
    logical :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = .true.
    end if

    end function DefaultTrue


    function PresentDefault_S(default, S) result(Sout)
    character(LEN=*), intent(in), target :: default
    character(LEN=*), intent(in), target, optional :: S
    character(LEN=:), pointer :: Sout

    if (present(S)) then
        SOut => S
    else
        SOut => default
    end if
    end function PresentDefault_S

    function PresentDefault_L(default, S) result(Sout)
    logical, intent(in) :: default
    logical, intent(in), optional :: S
    logical :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_L

    function PresentDefault_I(default, S) result(Sout)
    integer, intent(in) :: default
    integer, intent(in), optional :: S
    integer :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_I

    function PresentDefault_R(default, S) result(Sout)
    real, intent(in) :: default
    real, intent(in), optional :: S
    real:: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_R

    function PresentDefault_D(default, S) result(Sout)
    double precision, intent(in) :: default
    double precision, intent(in), optional :: S
    double precision:: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_D


    function IfThenElse_S(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    character(LEN=:), pointer :: IfThenElse
    character(LEN=*), target :: either, or

    if (flag) then
        IfThenElse => either
    else
        IfThenElse => or
    end if

    end function

    function IfThenElse_L(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    logical :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    function IfThenElse_I(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    integer :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    function IfThenElse_R(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    real :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    function IfThenElse_D(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    double precision :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    logical function isFloat0(R)
    class(*), intent(in) :: R

    select type(R)
    type is (real)
        isFloat0 = .true.
    type is (double precision)
        isFloat0 = .true.
        class default
        isFloat0 = .false.
    end select
    end function isFloat0

    logical function isFloat1(R)
    class(*), intent(in) :: R(:)

    isFloat1 = isFloat0(R(LBOUND(R,1)))

    end function isFloat1

    logical function isFloat2(R)
    class(*), intent(in) :: R(:,:)

    isFloat2 = isFloat0(R(LBOUND(R,1),LBOUND(R,2)))

    end function isFloat2

    function IndexOf(aval,arr, n)
        integer, intent(in) :: n, arr(n), aval
        integer :: IndexOf, i

        do i=1,n
            if (arr(i)==aval) then
                IndexOf= i
                return
            end if
        end do
        IndexOf = 0

    end function IndexOf
end module MiscUtils
