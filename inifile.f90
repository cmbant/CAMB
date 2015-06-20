    !Module to read in name/value pairs from a file, with each line of the form line 'name = value'
    !Should correctly interpret FITS headers
    !Antony Lewis (http://cosmologist.info/). Released to the public domain.
    !Apr 11, added support for INCLUDE(file); check for duplicate keys
    !Jun 12, added support for DEFAULT(file), provides default values for keys not defined in the read file

    module IniFile
    implicit none
    public

    integer, parameter :: Ini_max_name_len = 128

    integer, parameter :: Ini_max_string_len = 1024
    logical :: Ini_fail_on_not_found = .false.

    logical :: Ini_Echo_Read = .false.

    logical :: Ini_AllowDuplicateKeys = .false.

    type TNameValue
        !no known way to make character string pointers..
        character(Ini_max_name_len)  :: Name
        character(Ini_max_string_len):: Value
    end type TNameValue

    type TNameValue_pointer
        Type(TNameValue), pointer :: P
    end type TNameValue_pointer

    Type TNameValueList
        integer Count
        integer Delta
        integer Capacity
        logical ignoreDuplicates
        type(TNameValue_pointer), dimension(:), pointer :: Items
    end Type TNameValueList

    Type TIniFile
        logical SlashComments
        Type (TNameValueList) :: L, ReadValues
    end Type TIniFile

    Type(TIniFile) :: DefIni

    contains

    subroutine TNameValueList_Init(L, ignoreDuplicates)
    Type (TNameValueList) :: L
    logical, intent(in), optional :: ignoreDuplicates

    L%Count = 0
    L%Capacity = 0
    L%Delta = 128
    L%ignoreDuplicates = .false.
    if (present(ignoreDuplicates)) L%ignoreDuplicates=ignoreDuplicates
    nullify(L%Items)

    end subroutine TNameValueList_Init

    subroutine TNameValueList_Clear(L)
    Type (TNameValueList) :: L
    integer i, status

    do i=L%count,1,-1
        deallocate (L%Items(i)%P, stat = status)
    end do
    deallocate (L%Items, stat = status)
    call TNameValueList_Init(L)

    end subroutine TNameValueList_Clear

    subroutine TNameValueList_ValueOf(L, AName, AValue)
    Type (TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    CHARACTER(LEN=*), intent(out) :: AValue
    integer i

    do i=1, L%Count
        if (L%Items(i)%P%Name == AName) then
            AValue = L%Items(i)%P%Value
            return
        end if
    end do
    AValue = ''

    end subroutine TNameValueList_ValueOf


    function TNameValueList_IndexOf(L, AName) result (AValue)
    Type (TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    integer :: AValue
    integer i

    do i=1, L%Count
        if (L%Items(i)%P%Name == AName) then
            AValue = i
            return
        end if
    end do
    AValue = -1

    end function TNameValueList_IndexOf

    function TNameValueList_HasKey(L, AName) result (AValue)
    Type (TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    logical :: AValue

    AValue = TNameValueList_IndexOf(L,AName) /= -1

    end function TNameValueList_HasKey

    subroutine TNameValueList_Add(L, AName, AValue, only_if_undefined)
    Type (TNameValueList) :: L
    character(LEN=*), intent(in) :: AName, AValue
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .true.
    end if
    if ((.not. Ini_AllowDuplicateKeys .or. isDefault) .and. TNameValueList_HasKey(L,AName)) then
        if (L%ignoreDuplicates .or. isDefault) return
        write (*,*) 'IniFile,TNameValueList_Add: duplicate key name in file: '//trim(AName)
        stop
    end if
    if (L%Count == L%Capacity) call TNameValueList_SetCapacity(L, L%Capacity + L%Delta)
    L%Count = L%Count + 1
    allocate(L%Items(L%Count)%P)
    L%Items(L%Count)%P%Name = AName
    L%Items(L%Count)%P%Value = AValue

    end subroutine TNameValueList_Add

    subroutine TNameValueList_SetCapacity(L, C)
    Type (TNameValueList) :: L
    integer C
    type(TNameValue_pointer), dimension(:), pointer :: TmpItems

    if (L%Count > 0) then
        if (C < L%Count) stop 'TNameValueList_SetCapacity: smaller than Count'
        allocate(TmpItems(L%Count))
        TmpItems = L%Items(1:L%Count)
        deallocate(L%Items)
        allocate(L%Items(C))
        L%Items(1:L%Count) = TmpItems
        deallocate(TmpItems)
    else
        allocate(L%Items(C))
    end if
    L%Capacity = C

    end subroutine TNameValueList_SetCapacity

    subroutine TNameValueList_Delete(L, i)
    Type (TNameValueList) :: L
    integer, intent(in) :: i

    deallocate(L%Items(i)%P)
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1

    end subroutine TNameValueList_Delete

    subroutine Ini_NameValue_Add(Ini,AInLine,only_if_undefined)
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: AInLine
    integer EqPos, slashpos, lastpos
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault
    character (LEN=len(AInLine)) :: AName, S, InLine

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .false.
    end if
    InLine=trim(adjustl(AInLine))
    EqPos = scan(InLine,'=')
    if (EqPos/=0 .and. InLine(1:1)/='#' .and. InLine(1:7) /= 'COMMENT' ) then

        AName = trim(InLine(1:EqPos-1))

        S = adjustl(InLine(EqPos+1:))
        if (Ini%SlashComments) then
            slashpos=scan(S,'/')
            if (slashpos /= 0) then
                S  = S(1:slashpos-1)
            end if
        end if
        lastpos=len_trim(S)
        if (lastpos>1) then
            if (S(1:1)=='''' .and. S(lastpos:lastpos)=='''') then
                S = S(2:lastpos-1)
            end if
        end if
        call TNameValueList_Add(Ini%L, AName, S,only_if_undefined = isDefault )

    end if

    end subroutine Ini_NameValue_Add

    subroutine Ini_Open(filename, unit_id,  error, slash_comments)
    character (LEN=*), intent(IN) :: filename
    integer, intent(IN) :: unit_id
    logical, optional, intent(OUT) :: error
    logical, optional, intent(IN) :: slash_comments
    logical aerror

    call TNameValueList_Init(DefIni%L)
    call TNameValueList_Init(DefIni%ReadValues, .true.)

    if (present(slash_comments)) then
        call Ini_Open_File(DefIni,filename,unit_id,aerror,slash_comments)
    else
        call Ini_Open_File(DefIni,filename,unit_id,aerror)
    end if

    if (present(error)) then
        error = aerror
    else
        if (aerror) then
            write (*,*) 'Ini_Open: Error opening file ' // trim(filename)
            stop
        end if
    end if

    end subroutine Ini_Open

    function Ini_ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=Ini_max_string_len) Ini_ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='/') then
            Ini_ExtractFilePath = aname(1:i)
            return
        end if
    end do
    Ini_ExtractFilePath = ''

    end function Ini_ExtractFilePath

    recursive subroutine Ini_Open_File(Ini, filename, unit_id,  &
        error, slash_comments, append,only_if_undefined)
    Type(TIniFile) :: Ini

    character (LEN=*), intent(IN) :: filename
    integer, intent(IN) :: unit_id
    logical, intent(OUT) :: error
    logical, optional, intent(IN) :: slash_comments
    logical, optional, intent(in) :: append, only_if_undefined
    character (LEN=Ini_max_string_len) :: InLine, IncludeFile
    integer lastpos, i
    Type (TNameValueList) IncudeFiles, DefaultValueFiles
    logical doappend, FileExists, isDefault

    if (present(append)) then
        doappend=append
    else
        doappend=.false.
    end if

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .false.
    end if

    if (.not. doappend) then
        call TNameValueList_Init(Ini%L)
        call TNameValueList_Init(Ini%ReadValues, .true.)
    end if

    call TNameValueList_Init(IncudeFiles)
    call TNameValueList_Init(DefaultValueFiles)

    if (present(slash_comments)) then
        Ini%SlashComments = slash_comments
    else
        Ini%SlashComments = .false.
    end if

    open(unit=unit_id,file=filename,form='formatted',status='old', err=500)

    do
        read (unit_id,'(a)',end=400) InLine
        if (InLine == 'END') exit;
        if (InLine(1:8) == 'INCLUDE(') then
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call TNameValueList_Add(IncudeFiles, trim(adjustl(InLine(9:lastpos-1))),'')
            else
                stop 'Ini_Open_File: error in INCLUDE line'
            end if
        elseif (InLine(1:8) == 'DEFAULT(') then
            !Settings to read in as defaults, overridden by subsequent re-definitions
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call TNameValueList_Add(DefaultValueFiles, trim(adjustl(InLine(9:lastpos-1))),'')
            else
                stop 'Ini_Open_File: error in DEFAULT line'
            end if
        elseif (InLine /= '') then
            call Ini_NameValue_Add(Ini,InLine, only_if_undefined=isDefault)
        end if
    end do

400 close(unit_id)
    error=.false.

    do i=1, IncudeFiles%Count
        if (error) exit
        IncludeFile=IncudeFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(Ini_ExtractFilePath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                write(*,*) 'Ini_Open_File: INCLUDE file not found: '//trim(IncudeFiles%Items(i)%P%Name)
                stop
            end if
        end if
        call Ini_Open_File(Ini, IncludeFile, unit_id,  &
            error, slash_comments, append=.true.,only_if_undefined=isDefault)
    end do
    do i=1, DefaultValueFiles%Count
        if (error) exit
        IncludeFile=DefaultValueFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(Ini_ExtractFilePath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                write(*,*) 'Ini_Open_File: DEFAULT file not found:' //trim(DefaultValueFiles%Items(i)%P%Name)
                stop
            end if
        end if
        call Ini_Open_File(Ini, IncludeFile, unit_id,  &
            error, slash_comments, append=.true., only_if_undefined=.true.)
    end do

    call TNameValueList_Clear(IncudeFiles)
    call TNameValueList_Clear(DefaultValueFiles)

    return

500 error=.true.
    call TNameValueList_Clear(IncudeFiles)
    call TNameValueList_Clear(DefaultValueFiles)

    end subroutine Ini_Open_File

    subroutine Ini_Open_Fromlines(Ini, Lines, NumLines, slash_comments)
    Type(TIniFile) :: Ini

    integer, intent(IN) :: NumLines
    character (LEN=*), dimension(NumLines), intent(IN) :: Lines
    logical, intent(IN) :: slash_comments
    integer i

    call TNameValueList_Init(Ini%L)
    call TNameValueList_Init(Ini%ReadValues, .true.)

    Ini%SlashComments = slash_comments

    do i=1,NumLines
        call Ini_NameValue_Add(Ini,Lines(i))
    end do

    end  subroutine Ini_Open_Fromlines

    subroutine Ini_Close

    call Ini_close_File(DefIni)

    end subroutine Ini_Close


    subroutine Ini_Close_File(Ini)
    Type(TIniFile) :: Ini

    call TNameValueList_Clear(Ini%L)
    call TNameValueList_Clear(Ini%ReadValues)

    end  subroutine Ini_Close_File



    function Ini_Read_String(Key, NotFoundFail) result(AValue)
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=Ini_max_string_len) :: AValue

    if (present(NotFoundFail)) then
        AValue = Ini_Read_String_File(DefIni, Key, NotFoundFail)
    else
        AValue = Ini_Read_String_File(DefIni, Key)
    end if

    end function Ini_Read_String

    function Ini_Read_String_Default(Key, Default,AllowBlank) result(AValue)
    character (LEN=*), intent(IN) :: Key, Default
    character(LEN=Ini_max_string_len) :: AValue
    logical, intent(in), optional :: AllowBlank

    if (present(AllowBlank)) then
        AValue = Ini_Read_String_Default_File(DefIni, Key, Default, AllowBlank)
    else
        AValue = Ini_Read_String_Default_File(DefIni, Key, Default)
    end if

    end function Ini_Read_String_Default


    function Ini_Read_String_File(Ini, Key, NotFoundFail) result(AValue)
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=Ini_max_string_len) :: AValue

    call TNameValueList_ValueOf(Ini%L, Key, AValue)

    if (AValue/='') then

        call  TNameValueList_Add(Ini%ReadValues, Key, AValue)
        if (Ini_Echo_Read) write (*,*) trim(Key)//' = ',trim(AValue)
        return

    end if
    if (present(NotFoundFail)) then
        if (NotFoundFail) then
            write(*,*) 'key not found : '//trim(Key)
            stop
        end if
    else if (Ini_fail_on_not_found) then
        write(*,*) 'key not found : '//trim(Key)
        stop
    end if

    end function Ini_Read_String_File

    function Ini_Read_String_Default_File(Ini, Key, Default, AllowBlank) result(AValue)
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key, Default
    character(LEN=Ini_max_string_len) :: AValue
    logical, intent(in), optional :: AllowBlank

    if (Ini_HasKey_File(Ini,Key)) then
        AValue = Ini_Read_String_file(Ini, Key, .false.)
        if (present(AllowBlank)) then
            if (AllowBlank) return
        end if
        if (AValue=='') AValue = Default
    else
        AValue = Default
    end if

    end function Ini_Read_String_Default_File

    function Ini_HasKey(Key) result(AValue)
    character (LEN=*), intent(IN) :: Key
    logical AValue

    AValue = Ini_HasKey_File(DefIni, Key)

    end function Ini_HasKey

    function Ini_HasKey_File(Ini, Key) result(AValue)
    type(TIniFile), intent(in) :: Ini
    character (LEN=*), intent(IN) :: Key
    logical AValue

    Avalue = TNameValueList_HasKey(Ini%L, Key)

    end function Ini_HasKey_File

    function Ini_Key_To_Arraykey(Key, index)  result(AValue)
    character (LEN=*), intent(IN) :: Key
    integer, intent(in) :: index
    character(LEN=Ini_max_string_len) :: AValue

    character(LEN=32) :: numstr
    write (numstr,*) index
    numstr=adjustl(numstr)
    AValue = trim(Key) // '(' // trim(numStr) // ')'

    end function Ini_Key_To_Arraykey

    function Ini_Read_String_Array(Key, index, NotFoundFail) result(AValue)
    character (LEN=*), intent(IN) :: Key
    integer, intent(in) :: index
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=Ini_max_string_len) :: AValue

    if (present(NotFoundFail)) then
        AValue = Ini_Read_String_Array_File(DefIni, Key, index, NotFoundFail)
    else
        AValue = Ini_Read_String_Array_File(DefIni, Key, index)
    end if

    end function Ini_Read_String_Array

    function Ini_Read_String_Array_File(Ini, Key, index, NotFoundFail) result(AValue)
    Type(TIniFile) :: Ini
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=Ini_max_string_len) :: AValue
    character(LEN=Ini_max_string_len) :: ArrayKey

    ArrayKey = Ini_Key_To_Arraykey(Key,index)
    if (present(NotFoundFail)) then
        AValue = Ini_Read_String_File(Ini, ArrayKey, NotFoundFail)
    else
        AValue = Ini_Read_String_File(Ini, ArrayKey)
    end if

    end function Ini_Read_String_Array_File

    function Ini_Read_Int_Array(Key, index, Default)
    integer, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    integer Ini_Read_Int_Array

    if (present(Default)) then
        Ini_Read_Int_Array = Ini_Read_Int_Array_File(DefIni, Key, index, Default)
    else
        Ini_Read_Int_Array = Ini_Read_Int_Array_File(DefIni, Key, index)
    end if

    end function Ini_Read_Int_Array

    function Ini_Read_Int_Array_File(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    Type(TIniFile) :: Ini
    integer Ini_Read_Int_Array_File
    integer, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: ArrrayKey
    ArrrayKey = Ini_Key_To_Arraykey(Key,index)
    if (present(Default)) then
        Ini_Read_Int_Array_File = Ini_Read_Int_File(Ini, ArrrayKey, Default)
    else
        Ini_Read_Int_Array_File = Ini_Read_Int_File(Ini, ArrrayKey)
    end if
    end function Ini_Read_Int_Array_File


    function Ini_Read_Int(Key, Default)
    integer, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    integer Ini_Read_Int

    if (present(Default)) then
        Ini_Read_Int = Ini_Read_Int_File(DefIni, Key, Default)
    else
        Ini_Read_Int = Ini_Read_Int_File(DefIni, Key)
    end if
    end function Ini_Read_Int

    function Ini_Read_Int_File(Ini, Key, Default)
    Type(TIniFile) :: Ini
    integer Ini_Read_Int_File
    integer, optional, intent(IN) :: Default
    character  (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini_Read_String_File(Ini, Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Int_File = Default
        write (S,*) Default
        call  TNameValueList_Add(Ini%ReadValues, Key, S)
    else
        if (verify(trim(S),'-+0123456789') /= 0) goto 10
        read (S,*, err = 10) Ini_Read_Int_File
    end if
    return
10  write (*,*) 'error reading integer for key: '//Key
    stop

    end function Ini_Read_Int_File

    function Ini_Read_Double(Key, Default)
    double precision, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    double precision Ini_Read_Double

    if (present(Default)) then
        Ini_Read_Double = Ini_Read_Double_File(DefIni, Key, Default)
    else
        Ini_Read_Double = Ini_Read_Double_File(DefIni, Key)
    end if

    end function Ini_Read_Double

    function Ini_Read_Double_File(Ini,Key, Default)
    Type(TIniFile) :: Ini
    double precision Ini_Read_Double_File
    double precision, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini_Read_String_File(Ini,Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Double_File = Default
        write (S,*) Default

        call  TNameValueList_Add(Ini%ReadValues, Key, S)

    else
        read (S,*, err=10) Ini_Read_Double_File
    end if

    return

10  write (*,*) 'error reading double for key: '//Key
    stop

    end function Ini_Read_Double_File



    function Ini_Read_Double_Array(Key, index, Default)
    double precision, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    double precision Ini_Read_Double_Array

    if (present(Default)) then
        Ini_Read_Double_Array = Ini_Read_Double_Array_File(DefIni, Key, index, Default)
    else
        Ini_Read_Double_Array = Ini_Read_Double_Array_File(DefIni, Key, index)
    end if

    end function Ini_Read_Double_Array


    function Ini_Read_Double_Array_File(Ini,Key, index, Default)

    !Reads Key(1), Key(2), etc.

    Type(TIniFile) :: Ini

    double precision Ini_Read_Double_Array_File
    double precision, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) ::  ArrrayKey

    ArrrayKey = Ini_Key_To_Arraykey(Key,index)
    if (present(Default)) then

        Ini_Read_Double_Array_File = Ini_Read_Double_File(Ini, ArrrayKey, Default)
    else
        Ini_Read_Double_Array_File = Ini_Read_Double_File(Ini, ArrrayKey)
    end if
    end function Ini_Read_Double_Array_File

    function Ini_Read_Real(Key, Default)
    real, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    real Ini_Read_Real

    if (present(Default)) then
        Ini_Read_Real = Ini_Read_Real_File(DefIni, Key, Default)
    else
        Ini_Read_Real = Ini_Read_Real_File(DefIni, Key)
    end if

    end function Ini_Read_Real

    function Ini_Read_Real_File(Ini,Key, Default)
    Type(TIniFile) :: Ini
    real Ini_Read_Real_File
    real, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini_Read_String_File(Ini,Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Real_File = Default
        write (S,*) Default
        call  TNameValueList_Add(Ini%ReadValues, Key, S)

    else
        read (S,*, err=10) Ini_Read_Real_File
    end if

    return

10  write (*,*) 'error reading double for key: '//Key
    stop

    end function Ini_Read_Real_File



    function Ini_Read_Real_Array(Key, index, Default)
    real, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    real Ini_Read_Real_Array

    if (present(Default)) then
        Ini_Read_Real_Array = Ini_Read_Real_Array_File(DefIni, Key, index, Default)
    else
        Ini_Read_Real_Array = Ini_Read_Real_Array_File(DefIni, Key, index)
    end if
    end function Ini_Read_Real_Array

    function Ini_Read_Real_Array_File(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    Type(TIniFile) :: Ini
    real Ini_Read_Real_Array_File
    real, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: ArrrayKey

    ArrrayKey = Ini_Key_To_Arraykey(Key,index)
    if (present(Default)) then
        Ini_Read_Real_Array_File = Ini_Read_Real_File(Ini, ArrrayKey, Default)
    else
        Ini_Read_Real_Array_File = Ini_Read_Real_File(Ini, ArrrayKey)
    end if
    end function Ini_Read_Real_Array_File

    function Ini_Read_Logical(Key, Default)
    Logical, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    logical Ini_Read_Logical

    if (present(Default)) then
        Ini_Read_Logical = Ini_Read_Logical_File(DefIni, Key, Default)
    else
        Ini_Read_Logical = Ini_Read_Logical_File(DefIni, Key)
    end if
    end function Ini_Read_Logical

    function Ini_Read_Logical_File(Ini, Key, Default)
    Type(TIniFile) :: Ini

    logical Ini_Read_Logical_File
    logical, optional, intent(IN) :: Default
    character  (LEN=*), intent(IN) :: Key

    character(LEN=Ini_max_string_len) :: S

    S = Ini_Read_String_File(Ini,Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Logical_File = Default
        write (S,*) Default

        call  TNameValueList_Add(Ini%ReadValues, Key, S)

    else

        if (verify(trim(S),'10TF') /= 0) goto 10
        read (S,*, err = 10) Ini_Read_Logical_File
    end if

    return

10  write (*,*) 'error reading logical for key: '//Key
    stop
    end function Ini_Read_Logical_File



    subroutine Ini_SaveReadValues(afile,unit_id)
    character(LEN=*)  :: afile
    integer, intent(in) :: unit_id

    call Ini_SaveReadValues_File(DefIni, afile, unit_id)

    end subroutine Ini_SaveReadValues



    subroutine Ini_SaveReadValues_File(Ini, afile, unit_id)
    Type(TIniFile) :: Ini
    character(LEN=*), intent(in) :: afile
    integer, intent(in) :: unit_id
    integer i

    open(unit=unit_id,file=afile,form='formatted',status='replace', err=500)

    do i=1, Ini%ReadValues%Count

        write (unit_id,'(a)') trim(Ini%ReadValues%Items(i)%P%Name) // ' = ' &
            //trim(Ini%ReadValues%Items(i)%P%Value)

    end do

    close(unit_id)
    return

500 write(*,*) 'Ini_SaveReadValues_File: Error creating '//trim(afile)

    end subroutine Ini_SaveReadValues_File

    end module IniFile

