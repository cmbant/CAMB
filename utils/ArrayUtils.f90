    module ArrayUtils
    implicit none

    contains

    function IndexOf(aval,arr, n)
    integer, intent(in) :: n, arr(n), aval
    integer IndexOf, i

    do i=1,n
        if (arr(i)==aval) then
            IndexOf= i
            return
        end if
    end do
    IndexOf = 0

    end function IndexOf

    function MaxIndex(arr, n)
    integer, intent(in) :: n
    real, intent(in) :: arr(n)
    integer locs(1:1), MaxIndex

    locs = maxloc(arr(1:n))
    MaxIndex = locs(1)

    end function MaxIndex


    function MinIndex(arr, n)
    integer, intent(in) :: n
    real, intent(in) :: arr(n)
    integer locs(1:1), MinIndex

    locs = minloc(arr(1:n))
    MinIndex = locs(1)

    end function MinIndex

    end module ArrayUtils
