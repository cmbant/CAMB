    module StringUtils
    use MiscUtils
    implicit none


    INTERFACE RealToStr
    module procedure SingleToStr, DoubleToStr
    END INTERFACE RealToStr

    contains

    function IsWhiteSpace(C)
    character, intent(in) :: C
    logical IsWhiteSpace

    IsWhiteSpace = (C==' ') .or. (C==char(9))

    end function IsWhiteSpace

    function GetParamCount()
    integer GetParamCount

    GetParamCount = command_argument_count()

    end function GetParamCount

    function GetParam(i)
    character(LEN=:), allocatable :: GetParam
    integer, intent(in) :: i
    character(LEN=:), allocatable :: tmp
    integer l

    if (GetParamCount() < i) then
        GetParam = ''
    else
        call get_command_argument(i,length=l)
        allocate(character(l)::tmp)
        call get_command_argument(i,value=tmp)
        GetParam = trim(tmp)
    end if

    end function GetParam

    function GetEnvironmentVariable(name, is_present) result(value)
    character(LEN=*), intent(in) :: name
    character(LEN=:), allocatable :: value
    logical, intent(out), optional :: is_present
    integer L, status

    call get_environment_variable(name, length=L, status=status)
    if (present(is_present)) is_present = status==0
    if (status==0) then
        allocate(character(L)::value)
        call get_environment_variable(name, value, status=status)
    end if
    if (status/=0) value=''

    end function GetEnvironmentVariable


    function StringStarts(S, substring, index) result(OK)
    character(LEN=*), intent(in) :: S, substring
    integer, intent(in), optional :: index
    logical OK
    integer start

    if (present(index)) then
        start = index
    else
        start =1
    end if

    OK = S(start:min(len(S),start+len_trim(substring)-1))==substring

    end function StringStarts

    function StringTrimmed(S, trimmed) result(newS)
    character(LEN=*), intent(in) :: S
    logical, intent(in), optional :: trimmed
    character(LEN=:), allocatable :: newS

    if (DefaultFalse(trimmed)) then
        newS = trim(S)
    else
        newS = S
    end if

    end function StringTrimmed

    subroutine StringReplace(FindS, RepS, S)
    character(LEN=*), intent(in) :: FindS, RepS
    character(LEN=:), allocatable, intent(inout) :: S
    integer i

    i = index(S,FindS)
    if (i>0) then
        S = S(1:i-1)//trim(RepS)//S(i+len_trim(FindS):len_trim(S))
    end if

    end subroutine StringReplace

    function StringEscape(S, C, escape) result(newS)
    character(LEN=*), intent(in) :: S
    character, intent(in) :: C
    character, intent(in), optional :: escape
    character(LEN=:), allocatable :: newS, esc
    integer i
    character, parameter :: backslash = char(92)

    if (present(escape)) then
        esc= escape
    else
        esc = backslash
    end if
    newS = ''
    do i=1, len_trim(S)
        if (S(i:i)==C) then
            newS = newS //esc// C
        else
            newS = newS //S(i:i)
        end if
    end do

    end function StringEscape

    function Join(separator, S, S1,S2,S3,S4,S5,S6, trimmed) result(newS)
    character(LEN=*), intent(in) :: Separator, S, S1
    character(LEN=*), optional :: S2,S3,S4,S5,S6
    character(LEN=:), allocatable :: newS
    logical, intent(in), optional :: trimmed

    newS = StringTrimmed(S,trimmed)//Separator//StringTrimmed(S1,trimmed)
    if (present(S2)) newS = newS //Separator //StringTrimmed(S2,trimmed)
    if (present(S3)) newS = newS //Separator //StringTrimmed(S3,trimmed)
    if (present(S4)) newS = newS //Separator //StringTrimmed(S4,trimmed)
    if (present(S5)) newS = newS //Separator //StringTrimmed(S5,trimmed)
    if (present(S6)) newS = newS //Separator //StringTrimmed(S6,trimmed)

    end function Join

    function numcat(S, num)
    character(LEN=*) S
    character(LEN=:), allocatable :: numcat
    integer num

    numcat = concat(S,num)
    end function numcat

    function IntToStr(I, minlen)
    integer , intent(in) :: I
    character(LEN=:), allocatable :: IntToStr
    integer, intent(in), optional :: minlen
    integer n
    character (LEN=128) :: form, tmp

    if (present(minlen)) then
        n = minlen
        if (I<0) n=n+1
        form = concat('(I',n,'.',minlen,')')
        write (tmp,form) i
        IntToStr = trim(tmp)
    else
        write (tmp,*) i
        IntToStr = trim(adjustl(tmp))
    end if

    end function IntToStr

    function StrToInt(S)
    integer :: StrToInt
    character(LEN=*), intent(in) :: S

    read(S,*) StrToInt

    end function StrToInt

    subroutine StringAppend(S,X)
    character(LEN=:), allocatable, intent(inout) :: S
    class(*), intent(in) :: X

    if (.not. allocated(S)) S=''
    select type (X)
    type is (character(LEN=*))
        S = S // trim(X)
    type is (integer)
        S = S // IntToStr(X)
    type is (real)
        S = S // RealToStr(X)
    type is (double precision)
        S=S //RealToStr(X)
        class default
        stop 'StringAppend: Unknown type'
    end select
    end subroutine

    function concat(S1,S2,S3,S4,S5,S6,S7,S8) result(outstr)
    character(LEN=*), intent(in) :: S1
    class(*), intent(in) :: S2
    class(*), intent(in), optional :: S3, S4, S5, S6,S7,S8
    character(LEN=:), allocatable :: outstr

    outstr=S1
    call StringAppend(outstr,S2)
    if (present(S3)) then
        call StringAppend(outstr,S3)
        if (present(S4)) then
            call StringAppend(outstr,S4)
            if (present(S5)) then
                call StringAppend(outstr,S5)
                if (present(S6)) then
                    call StringAppend(outstr,S6)
                    if (present(S7)) then
                        call StringAppend(outstr,S7)
                        if (present(S8)) then
                            call StringAppend(outstr,S8)
                        end if
                    end if
                end if
            end if
        end if
    end if

    end function concat


    function DoubleToStr(R, figs)
    double precision, intent(in) :: R
    integer, intent(in), optional :: figs
    character(LEN=:), allocatable :: DoubleToStr

    DoubleToStr = SingleToStr(real(R),figs)

    end function DoubleToStr

    function SingleToStr(R, figs)
    real, intent(in) :: R
    integer, intent(in), optional :: figs
    character(LEN=:), allocatable :: SingleToStr
    character(LEN=30) tmp

    if (abs(R)>=0.001 .or. R==0.) then
        write (tmp,'(f12.6)') R

        tmp = adjustl(tmp)
        if (present(figs)) then
            SingleToStr = tmp(1:figs)
        else
            if (abs(R)>10000) then
                write(tmp,*) R
                SingleToStr =  trim(adjustl(tmp))
            else
                SingleToStr = tmp(1:6)
            end if
        end if
    else
        if (present(figs)) then
            write (tmp,trim(numcat('(E',figs))//'.2)') R
        else
            write (tmp,'(G9.2)') R
        end if
        SingleToStr = trim(adjustl(tmp))
    end if

    end function SingleToStr

    function SubNextFormat(S, X) result(OK)
    character(LEN=:), allocatable :: S
    class(*) X
    logical OK
    integer ix, P, n
    character c
    character(LEN=:), allocatable :: form, fform, rep

    P=1
    do
        ix=scan(S(P:),'%')
        OK = ix/=0 .and. ix < len(S)
        if (.not. OK) return
        c = S(ix+P:ix+P)
        if (c=='%') then
            P=P+Ix+1
        else
            exit
        end if
    end do
    form = ''
    do while( verify(c,'0123456789') == 0)
        form = form // c
        P=P+1
        c= S(ix+P:ix+P)
    end do
    select type (X)
    type is (integer)
        if (len(form)>0) then
            n= StrToInt(form)
            fform = 'I'//IntToStr(n)
            if (form(1:1)=='0') fform=fform//'.'//IntToStr(n)
            allocate(character(n)::rep)
            write(rep,'('//fform//')') X
        else
            rep = IntToStr(X)
        end if
        if (c=='d' .or. c=='u') then
            call StringReplace('%'//form//c, rep, S)
        else
            write(*,*) 'Wrong format for type: '//trim(S)
            stop
        end if
    type is (Character(LEN=*))
        if (c/='s') then
            write(*,*) 'Wrong format for type: '//trim(S)
            stop
        end if
        call StringReplace('%s', X, S)
    type is (double precision)
        if (c/='f') then
            write(*,*) 'Wrong format for type: '//trim(S)
            stop
        end if
        call StringReplace('%f', RealToStr(X),S)
    type is (real)
        if (c/='f') then
            write(*,*) 'Wrong format for type: '//trim(S)
            stop
        end if
        call StringReplace('%f', RealToStr(X),S)
        class default
        stop 'Unsupported format type'
    end select

    end function SubNextFormat

    function FormatString(formatst, i1,i2,i3,i4,i5,i6,i7,i8, allow_unused) result(S)
    character(LEN=*), intent(in) :: formatst
    class(*), intent(in),optional :: i1,i2,i3,i4,i5,i6,i7,i8
    logical, optional, intent(in) :: allow_unused
    character(LEN=:), allocatable :: S
    logical OK
    !Note that this routine is incomplete and very simple (so buggy in complex cases)
    !(should not substitute on the previously substituted string, etc, etc..)
    !Can do things like FormatString('case %d, ans = %03d%%',i,percent)
    S = formatst
    OK = .true.
    if (present(i1)) OK = SubNextFormat(S, i1)
    if (OK .and. present(i2)) OK = SubNextFormat(S, i2)
    if (OK .and. present(i3)) OK = SubNextFormat(S, i3)
    if (OK .and. present(i4)) OK = SubNextFormat(S, i4)
    if (OK .and. present(i5)) OK = SubNextFormat(S, i5)
    if (OK .and. present(i6)) OK = SubNextFormat(S, i6)
    if (OK .and. present(i7)) OK = SubNextFormat(S, i7)
    if (OK .and. present(i8)) OK = SubNextFormat(S, i8)
    if (.not. OK .and. .not. DefaultFalse(allow_unused)) stop 'FormatString: Wrong number or kind of formats in string'
    call StringReplace('%%', '%', S)

    end function FormatString


    end module StringUtils
