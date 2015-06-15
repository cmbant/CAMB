    module FileUtils
    use MpiUtils
    use MiscUtils
    use StringUtils
    use, intrinsic :: ISO_FORTRAN_ENV, only : INT64
    implicit none
    !Utils using F2008 features
    private

    integer, parameter :: File_size_int = KIND(INT64)

    Type :: TFileStream
        integer :: unit = 0
        logical :: ReadOnly = .true.
        character(LEN=:), allocatable :: mode
        character(LEN=:), allocatable :: access !sequential (fortran or text) or stream (binary, with file position)
        character(LEN=:), allocatable :: FileName
    contains
    procedure :: SetDefaultModes
    procedure :: Open => TFileStream_Open
    procedure :: OpenFile => TFileStream_OpenFile
    procedure :: Flush => TFileStream_Flush
    procedure :: CreateFile
    procedure :: CreateOpenFile
    procedure :: Position => TFileStream_Position
    procedure :: Size => TFileStream_Size
    procedure :: Close => TFileStream_Close
    procedure :: Rewind => TFileStream_Rewind
    procedure :: Error
    procedure :: Opened
    procedure :: WriteTrim
    procedure :: CheckOpen
    procedure, private :: ReadStringFunc
    procedure, private :: ReadStringSub
    procedure, private :: ReadItemSub
    procedure, private :: ReadItemFunc
    procedure, private :: ReadItems
    procedure, private :: ReadArray
    procedure, private :: ReadArray2
    procedure, private :: ReadArray2Func
    procedure, private :: ReadArrayFunc
    procedure, private :: WriteItemSub
    procedure, private :: WriteArray
    procedure, private :: WriteOneAndArray
    procedure, private :: WriteArray2
    procedure, private :: WriteItems
    procedure, private  :: WriteSizedArray1
    procedure, private  :: WriteSizedArray2
    procedure, private  :: ReadSizedArray_R
    procedure, private  :: ReadSizedArray_D
    procedure, private  :: ReadSizedArray_I
    procedure, private  :: ReadSizedArray2_I
    procedure, private  :: ReadSizedArray2_D
    procedure, private  :: ReadSizedArray2_R
    generic :: Write => WriteItems, WriteArray, WriteArray2, WriteOneAndArray
    generic :: Read => ReadItems, ReadArray, ReadArray2
    generic :: ReadItem => ReadItemFunc, ReadArrayFunc, ReadArray2Func
    generic :: ReadString => ReadStringSub
    generic :: ReadStringItem => ReadStringFunc
    generic :: WriteSizedArray => WriteSizedArray1,WriteSizedArray2
    generic :: ReadSizedArray => ReadSizedArray_R,ReadSizedArray_D,ReadSizedArray_I, &
        & ReadSizedArray2_I, ReadSizedArray2_R, ReadSizedArray2_D
    final :: TFileStream_Free
    end type

    Type, extends(TFileStream) :: TBinaryFile
    end type

    Type, extends(TFileStream) :: TSequentialFile
    contains
    procedure :: SetDefaultModes => TSequentialFile_SetDefaultModes
    end type TSequentialFile

    Type, extends(TSequentialFile) :: TTextFile
        character(LEN=20) :: RealFormat = '(*(E17.7))'
        character(LEN=20) :: IntegerFormat = '(*(I10))'
        logical :: AdvanceDefault = .true.
    contains
    procedure :: SetDefaultModes => TTextFile_SetDefaultModes
    procedure :: ReadLine
    procedure :: ReadLineSkipEmptyAndComments
    procedure :: ReadNextContentLine !as above, but opens and closes file
    procedure :: NewLine
    procedure :: SkipLines
    procedure :: Lines
    procedure :: Columns
    procedure :: WriteLeftAligned
    procedure :: WriteItemTxt
    procedure :: WriteArrayTxt
    procedure :: WriteInLineItems
    procedure :: WriteInLineTrim
    procedure :: WriteFormat
    procedure :: WriteTrim => WriteTrimTxt
    procedure, private :: ReadStringSub => ReadStringTxt
    procedure, private :: ReadItemTxt
    procedure, private :: ReadArrayTxt
    procedure, private :: WriteInLineItem
    procedure, private :: WriteInLineArray
    procedure, private :: DefaultAdvance
    procedure, private :: WriteItem => TTextFile_WriteItem
    procedure, private :: WriteArray => TTextFile_WriteArray
    procedure, private :: WriteOneAndArray => TTextFile_WriteOneAndArray
    procedure, private :: WriteArray2 => WriteArray2Txt
    procedure, private :: ReadItemSub => ReadItemTxt
    procedure, private :: ReadArray => ReadArrayTxt
    procedure, private :: ReadArray2 => ReadArray2Txt
    procedure, private :: WriteItems => WriteItemsTxt
    generic :: WriteInLine => WriteInLineItem, WriteInLineArray
    end type

    !Functions on filenames and text
    !File instance below acts as a namespace
    Type TFile
    contains
    procedure, nopass :: TxtNumberColumns
    procedure, nopass :: TxtColumns
    procedure, nopass :: TxtFileColumns
    procedure, nopass :: TxtFileLines
    procedure, nopass :: LastTopComment
    procedure, nopass :: TopCommentLine
    procedure, nopass :: LastLine => LastFileLine
    procedure, nopass :: Size => FileSize
    procedure, nopass :: Exists => FileExists
    procedure, nopass :: ExtractName => ExtractFileName
    procedure, nopass :: ExtractPath => ExtractFilePath
    procedure, nopass :: Join => File_Join
    procedure, nopass :: ChangeExt => ChangeFileExt
    procedure, nopass :: CheckTrailingSlash
    procedure, nopass :: IsFullPath
    procedure, nopass :: ExtractExt => ExtractFileExt
    procedure, nopass :: Delete => DeleteFile
    procedure, nopass :: ReadTextMatrix
    procedure, nopass :: ReadTextVector
    procedure, nopass :: WriteTextMatrix
    procedure, nopass :: WriteTextVector
    procedure, nopass :: LoadTxt_1D
    procedure, nopass :: LoadTxt_2D
    procedure, nopass :: OpenTextFile
    procedure, nopass :: CreateTextFile
    procedure, nopass :: CharIsSlash
    generic  :: LoadTxt => LoadTxt_2D, LoadTxt_1D
    generic  :: SaveTxt => WriteTextMatrix, WriteTextVector
    end type

    type(TFile), public, save :: File

    public TFileStream, TBinaryFile, TTextFile, TFile, File_size_int
    contains

    function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists

    inquire(file=aname, exist = FileExists)

    end function FileExists

    function FileSize(name) result(fsize)
    integer(file_size_int) fsize
    character(LEN=*), intent(in)::name

    inquire(file=name, size=fsize)

    end function FileSize

    function TxtNumberColumns(InLine) result(n)
    character(LEN=*) :: InLine
    integer n,i
    logical isNum

    n=0
    isNum=.false.
    do i=1, len_trim(InLIne)
        if (verify(InLine(i:i),'-+eE.0123456789') == 0) then
            if (.not. IsNum) n=n+1
            IsNum=.true.
        else
            IsNum=.false.
        end if
    end do

    end function TxtNumberColumns

    function TxtColumns(InLine) result(n)
    character(LEN=*) :: InLine
    integer n,i
    logical isNum

    n=0
    isNum=.false.
    do i=1, len_trim(InLine)
        if (InLine(i:i) > char(32)) then
            if (.not. IsNum) n=n+1
            IsNum=.true.
        else
            IsNum=.false.
        end if
    end do

    end function TxtColumns


    subroutine SetDefaultModes(this)
    class(TFileStream) :: this

    if (.not. allocated(this%mode)) this%mode = 'unformatted'
    if (.not. allocated(this%access)) this%access = 'stream'

    end subroutine SetDefaultModes

    subroutine TSequentialFile_SetDefaultModes(this)
    class(TSequentialFile) :: this

    if (.not. allocated(this%access)) this%access= 'sequential'
    call this%TFileStream%SetDefaultModes

    end subroutine TSequentialFile_SetDefaultModes

    subroutine TTextFile_SetDefaultModes(this)
    class(TTextFile) :: this

    if (.not. allocated(this%mode)) this%mode = 'formatted'
    call this%TSequentialFile%SetDefaultModes

    end subroutine TTextFile_SetDefaultModes


    subroutine CheckOpen(this, forWrite)
    class(TFileStream) :: this
    logical, intent(in), optional :: forWrite

    if (this%unit/=0) return
    if (DefaultFalse(forWrite) .and. this%ReadOnly) then
        call this%Error('File not open for write')
    else
        call this%Error('File not opened')
    end if

    end subroutine

    function Opened(this)
    class(TFileStream) :: this
    logical Opened
    Opened = this%unit /=0
    end function Opened

    subroutine Error(this, msg, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: msg
    character(LEN=*), intent(IN), optional :: errormsg
    character(LEN=:), allocatable :: Filename

    if (.not. allocated(this%FileName)) then
        FileName = '(no filename set)'
    else
        FileName = this%FileName
    end if
    if (present(errormsg)) then
        call MpiStop(trim(errormsg)//' : '//FileName )
    else
        call MpiStop(trim(msg)//' : '// FileName )
    end if

    end subroutine

    subroutine TFileStream_Close(this)
    class(TFileStream), intent(inout) :: this
    if (this%unit/=0) close(this%unit)
    this%unit=0
    end subroutine TFileStream_Close

#ifdef __GFORTRAN__
    impure elemental subroutine TFileStream_Free(this)
#else
    subroutine TFileStream_Free(this)
#endif
    Type(TFileStream), intent(inout) :: this

    call this%Close()
    end subroutine TFileStream_Free

    subroutine TFileStream_Flush(this)
    class(TFileStream) :: this

    call this%CheckOpen
    flush(this%unit)

    end subroutine TFileStream_Flush

    subroutine TFileStream_Rewind(this)
    class(TFileStream) :: this

    if (this%Opened()) rewind(this%unit)

    end subroutine TFileStream_Rewind

    function TFileStream_Position(this)
    class(TFileStream) :: this
    integer(file_size_int) TFileStream_Position

    call this%CheckOpen
    if (this%access /= 'stream') call this%Error('Position requires access=stream')
    inquire(this%unit, pos=TFileStream_Position)

    end function TFileStream_Position


    function TFileStream_Size(this)
    class(TFileStream) :: this
    integer(file_size_int) TFileStream_Size

    if (this%Opened()) then
        inquire(this%unit, size=TFileStream_Size)
    else if (allocated(this%FileName)) then
        TFileStream_Size = File%Size(this%FileName)
    else
        call this%Error('File not defined for size')
    end if

    end function TFileStream_Size

    subroutine TFileStream_Open(this, aname, errormsg, status)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg
    integer, intent(out), optional :: status

    call this%OpenFile(aname,errormsg=errormsg, status=status)

    end subroutine TFileStream_Open

    subroutine TFileStream_OpenFile(this, aname, mode, errormsg, forwrite, append, status)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: mode
    character(LEN=*), intent(IN), optional :: errormsg
    integer, intent(out), optional :: status
    logical, intent(in), optional :: forwrite, append
    character(LEN=:), allocatable :: amode, state, action, pos
    integer out_status

    call this%Close()
    call this%SetDefaultModes()

    this%FileName = trim(aname)

    amode = PresentDefault(this%mode, mode)
    if (DefaultFalse(forwrite)) then
        state = 'replace'
        action = 'readwrite'
        this%ReadOnly = .false.
    else
        state = 'old'
        action = 'read'
        this%ReadOnly = .true.
    end if
    if (DefaultFalse(append) .and. FileExists(aname)) then
        pos = 'append'
        state = 'old'
        action = 'readwrite'
        this%ReadOnly = .false.
    else
        pos='asis'
    end if

    open(file=aname,form=amode,status=state, action=action, newunit=this%unit, &
        & iostat=out_status, position =pos,  access=this%access)
    if (present(status)) then
        status=out_status
        if (out_status/=0) this%unit = 0
    else
        if (out_status/=0) then
            if (state == 'replace') then
                call this%Error('Error creating file', errormsg)
            else
                call this%Error('File not found', errormsg)
            end if
            this%unit = 0
        end if
    end if
    end subroutine TFileStream_OpenFile

    subroutine CreateFile(this,aname, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, errormsg=errormsg, forwrite =.true.)

    end subroutine CreateFile

    subroutine CreateOpenFile(this, aname, append, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    logical, optional, intent(in) :: append
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, forwrite =.true., append=append, errormsg=errormsg)

    end subroutine CreateOpenFile


    subroutine ReadStringSub(this, S, OK)
    class(TFileStream) :: this
    character(LEN=:), allocatable :: S
    logical, optional :: OK
    logical isOK
    integer i, status

    call this%CheckOpen
    read(this%unit, iostat=status) i
    isOK = status==0
    if (isOK) then
        if (allocated(S)) deallocate(S)
        allocate(character(LEN=i)::S)
        read(this%unit, iostat=status) S
        isOK = status==0
    end if
    if (present(OK)) OK = isOK
    end subroutine ReadStringSub

    function ReadStringFunc(this) result(S)
    class(TFileStream) :: this
    character(LEN=:), allocatable :: S
    call this%ReadString(S)
    end function ReadStringFunc

    subroutine ReadItemSub(this, R, OK)
    class(TFileStream) :: this
    class(*), intent(out) :: R
    logical, optional :: OK
    logical res
    integer status

    call this%CheckOpen
    select type(R)
    type is (real)
        read(this%unit, iostat=status) R
    type is (double precision)
        read(this%unit, iostat=status) R
    type is (integer)
        read(this%unit, iostat=status) R
    type is (logical)
        read(this%unit, iostat=status) R
    type is (character(LEN=*))
        read(this%unit, iostat=status) R
        class default
        call this%Error('Unknown type to read')
    end select

    res = status==0
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
        & call this%Error('Error reading item')
    if (present(OK)) OK = res
    end subroutine ReadItemSub


    function ReadItemFunc(this, R) result(res)
    class(TFileStream) :: this
    class(*), intent(out) :: R
    logical :: res

    call this%ReadItemSub(R, res)
    end function ReadItemFunc

    subroutine ReadArray(this, R, n, OK)
    class(TFileStream) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical, optional :: OK
    integer status

    call this%CheckOpen
    select type(R)
    type is (real)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (double precision)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (integer)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (logical)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
        class default
        call this%Error('Unknown type to read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
        & call this%Error('Error reading item')
    if (present(OK)) OK = status==0
    end subroutine ReadArray


    function ReadArrayFunc(this, R, n) result(res)
    class(TFileStream) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical :: res

    call this%ReadArray(R,n,res)

    end function ReadArrayFunc


    function ReadArray2Func(this, R) result(res)
    class(TFileStream) :: this
    class(*) :: R(:,:)
    logical :: res
    call this%ReadArray2(R,res)
    end function ReadArray2Func


    subroutine ReadArray2(this, R, OK)
    class(TFileStream) :: this
    class(*) :: R(:,:)
    logical, optional :: OK
    integer status

    call this%CheckOpen
    select type(R)
    type is (real)
        read(this%unit, iostat=status) R
    type is (double precision)
        read(this%unit, iostat=status) R
    type is (integer)
        read(this%unit, iostat=status) R
    type is (logical)
        read(this%unit, iostat=status) R
        class default
        call this%Error('Unknown type to read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
        & call this%Error('Error reading item')
    if (present(OK)) OK = status==0
    end subroutine ReadArray2

    subroutine ReadSizedArray_R(this, R)
    class(TFileStream) :: this
    Real, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_R


    subroutine ReadSizedArray_D(this, R)
    class(TFileStream) :: this
    double precision, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_D

    subroutine ReadSizedArray_I(this, R)
    class(TFileStream) :: this
    integer, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_I

    subroutine ReadSizedArray2_R(this, R)
    class(TFileStream) :: this
    real, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_R

    subroutine ReadSizedArray2_D(this, R)
    class(TFileStream) :: this
    double precision, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_D

    subroutine ReadSizedArray2_I(this, R)
    class(TFileStream) :: this
    integer, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_I

    subroutine WriteItemSub(this, R)
    class(TFileStream) :: this
    class(*), intent(in) :: R

    call this%CheckOpen(.true.)
    select type(R)
    type is (real)
        Write(this%unit) R
    type is (double precision)
        Write(this%unit) R
    type is (integer)
        Write(this%unit) R
    type is (logical)
        Write(this%unit) R
    type is (character(LEN=*))
        Write(this%unit) len(R)
        Write(this%unit) R
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteItemSub

    subroutine WriteTrim(this, S)
    class(TFileStream) :: this
    character(LEN=*), intent(in) :: S
    call this%WriteItemSub(trim(S))
    end subroutine WriteTrim

    subroutine WriteSizedArray1(this, R, n)
    class(TFileStream) :: this
    class(*), intent(in) :: R(1:)
    integer, intent(in), optional :: n
    integer sz

    sz=PresentDefault(size(R),n)
    call this%Write(sz)
    call this%WriteArray(R,n)

    end subroutine WriteSizedArray1

    subroutine WriteSizedArray2(this, R)
    class(TFileStream) :: this
    class(*), intent(in) :: R(:,:)

    call this%Write(size(R,dim=1))
    call this%Write(size(R,dim=2))
    call this%WriteArray2(R)

    end subroutine WriteSizedArray2

    subroutine WriteOneAndArray(this, I, R)
    class(TFileStream) :: this
    class(*), intent(in) :: I
    class(*), intent(in) :: R(:)

    call this%WriteItemSub(I)
    call this%WriteArray(R)

    end subroutine WriteOneAndArray

    subroutine WriteArray(this, R, n)
    class(TFileStream) :: this
    class(*), intent(in) :: R(1:)
    integer, intent(in), optional :: n
    integer sz

    sz=PresentDefault(size(R),n)
    call this%CheckOpen(.true.)
    select type(R)
    type is (real)
        Write(this%unit) R(1:sz)
    type is (double precision)
        Write(this%unit) R(1:sz)
    type is (integer)
        Write(this%unit) R(1:sz)
    type is (logical)
        Write(this%unit) R(1:sz)
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteArray

    subroutine WriteArray2(this, R)
    class(TFileStream) :: this
    class(*), intent(in) :: R(:,:)

    call this%CheckOpen(.true.)
    select type(R)
    type is (real)
        Write(this%unit) R
    type is (double precision)
        Write(this%unit) R
    type is (integer)
        Write(this%unit) R
    type is (logical)
        Write(this%unit) R
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteArray2

    subroutine WriteItems(this, S1, S2,S3,S4,S5,S6)
    class(TFileStream) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteItemSub(S1)
    if (present(S2)) call this%WriteItemSub(S2)
    if (present(S3)) call this%WriteItemSub(S3)
    if (present(S4)) call this%WriteItemSub(S4)
    if (present(S5)) call this%WriteItemSub(S5)
    if (present(S6)) call this%WriteItemSub(S6)

    end subroutine WriteItems

    subroutine ReadItems(this, S1, S2,S3,S4,S5,S6,OK)
    class(TFileStream) :: this
    class(*) S1
    class(*), optional :: S2,S3,S4,S5,S6
    logical, optional :: OK

    call this%ReadItemSub(S1,OK)
    if (present(S2) .and. DefaultTrue(OK)) call this%ReadItemSub(S2,OK)
    if (present(S3) .and. DefaultTrue(OK)) call this%ReadItemSub(S3,OK)
    if (present(S4) .and. DefaultTrue(OK)) call this%ReadItemSub(S4,OK)
    if (present(S5) .and. DefaultTrue(OK)) call this%ReadItemSub(S5,OK)
    if (present(S6) .and. DefaultTrue(OK)) call this%ReadItemSub(S6,OK)

    end subroutine ReadItems

    !Text unformatted files

    function ReadLine(this, InLine, trimmed) result(OK)
    class(TTextFile) :: this
    character(LEN=:), allocatable, optional :: InLine
    logical, intent(in), optional :: trimmed
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size

    call this%CheckOpen
    OK = .false.
    set = .true.
    do
        read (this%unit,'(a)',advance='NO',iostat=status, size=size) InS
        OK = .not. IS_IOSTAT_END(status)
        if (.not. OK) return
        if (present(InLine)) then
            if (set) then
                InLine = InS(1:size)
                set=.false.
            else
                InLine = InLine // InS(1:size)
            end if
        end if
        if (IS_IOSTAT_EOR(status)) exit
    end do
    if (present(trimmed) .and. present(InLine)) then
        if (trimmed) InLine = trim(adjustl(InLine))
    end if

    end function ReadLine

    function ReadLineSkipEmptyAndComments(this, InLine, comment) result(OK)
    class(TTextFile) :: this
    logical :: OK
    character(LEN=:), allocatable :: InLine
    character(LEN=:), allocatable, optional, intent(inout) :: comment

    if (present(comment)) then
        if (.not. allocated(comment)) comment=''
    end if
    do
        OK = this%ReadLine(InLine, trimmed=.true.)
        if (.not. OK) return
        if (InLine=='') cycle
        if (InLine(1:1)/='#') then
            return
        else
            if (present(comment)) comment = trim(InLine(2:))
        end if
    end do

    end function ReadLineSkipEmptyAndComments

    function ReadNextContentLine(this, filename, InLine) result(OK)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: filename
    character(LEN=:), intent(out), allocatable :: InLine
    logical :: OK

    if (.not. this%Opened()) call this%Open(filename)
    OK = this%ReadLineSkipEmptyAndComments(InLine)
    if (.not. OK) call this%Close()

    end function ReadNextContentLine

    function SkipLines(this, n) result(OK)
    class(TTextFile) :: this
    integer, intent(in) :: n
    logical OK
    integer ix

    do ix = 1, n
        if (.not. this%ReadLine()) then
            OK = .false.
            return
        end if
    end do
    OK = .true.

    end function SkipLines

    function Columns(this) result(n)
    class(TTextFile) :: this
    integer n
    character(LEN=:), allocatable :: InLine

    if (this%ReadLineSkipEmptyAndComments(InLine)) then
        n = File%TxtNumberColumns(InLine)
    else
        n=0
    end if
    call this%Rewind()

    end function Columns

    function Lines(this, nocomments) result(n)
    class(TTextFile) :: this
    logical, intent(in), optional :: nocomments
    integer n
    character(LEN=:), allocatable :: InLine

    n=0
    if (DefaultTrue(nocomments)) then
        do while (this%ReadLineSkipEmptyAndComments(InLine))
            n = n+1
        end do
    else
        do while (this%ReadLine())
            n = n+1
        end do
    end if
    call this%Rewind()

    end function Lines

    function DefaultAdvance(this,advance)
    class(TTextFile) :: this
    logical, intent(in), optional :: advance
    character(3) :: DefaultAdvance

    if (PresentDefault(this%AdvanceDefault,advance)) then
        DefaultAdvance='YES'
    else
        DefaultAdvance='NO'
    end if

    end function

    subroutine NewLine(this)
    class(TTextFile) :: this
    call this%WriteItemTxt('')
    end subroutine

    subroutine WriteLeftAligned(this, Form, str)
    class(TTextFile) :: this
    character(LEN=*) str, Form
    Character(LEN=max(len(str),128)) tmp

    call this%CheckOpen(.true.)
    tmp = str
    write(this%unit,form, advance='NO') tmp

    end subroutine WriteLeftAligned


    subroutine WriteInLineItems(this, S1, S2,S3,S4,S5,S6)
    class(TTextFile) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteInLine(S1)
    if (present(S2)) call this%WriteInLine(S2)
    if (present(S3)) call this%WriteInLine(S3)
    if (present(S4)) call this%WriteInLine(S4)
    if (present(S5)) call this%WriteInLine(S5)
    if (present(S6)) call this%WriteInLine(S6)

    end subroutine WriteInLineItems

    subroutine WriteItemsTxt(this, S1, S2,S3,S4,S5,S6)
    class(TTextFile) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteInLineItems(S1, S2,S3,S4,S5,S6)
    call this%NewLine()

    end subroutine WriteItemsTxt

    subroutine WriteInLineItem(this,str,form)
    class(TTextFile) :: this
    class(*), intent(in) :: str
    character(LEN=*), intent(in), optional :: form
    call this%WriteItemTxt(str,form,.false.)
    end subroutine WriteInLineItem

    subroutine WriteInLineArray(this,str,form,n)
    class(TTextFile) :: this
    class(*), intent(in) :: str(:)
    character(LEN=*), intent(in), optional :: form
    integer, intent(in), optional :: n

    call this%WriteArrayTxt(str,form,.false.,number=n)
    end subroutine WriteInLineArray

    subroutine TTextFile_WriteItem(this,R)
    class(TTextFile) :: this
    class(*), intent(in) :: R
    call this%WriteItemTxt(R,advance=.true.)
    end subroutine TTextFile_WriteItem

    subroutine TTextFile_WriteArray(this,R,n)
    class(TTextFile) :: this
    class(*), intent(in) :: R(:)
    integer, intent(in), optional :: n
    call this%WriteArrayTxt(R,number=n,advance=.true.)
    end subroutine TTextFile_WriteArray


    subroutine TTextFile_WriteOneAndArray(this, I, R)
    class(TTextFile) :: this
    class(*), intent(in) :: I
    class(*), intent(in) :: R(:)

    call this%WriteInlineItem(I)
    call this%WriteArrayTxt(R)

    end subroutine TTextFile_WriteOneAndArray

    subroutine WriteItemTxt(this, str, form, advance)
    class(TTextFile) :: this
    class(*), intent(in) :: str
    character(LEN=*), intent(in), optional :: form
    logical, intent(in), optional :: advance
    character(LEN=3) :: Ad

    call this%CheckOpen(.true.)
    Ad = this%DefaultAdvance(advance)
    select type(str)
    type is (character(LEN=*))
        if (Ad=='YES') then
            write(this%unit,PresentDefault('(a)',form)) trim(str)
        else
            write(this%unit,PresentDefault('(a)',form), advance=Ad) str
        end if
    type is (real)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str
    type is (double precision)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str
    type is (integer)
        write(this%unit,PresentDefault(this%IntegerFormat,form), advance=Ad) str
    type is (logical)
        write(this%unit,PresentDefault('(L2)',form), advance=Ad) str
        class default
        call this%Error('unknown type to write')
    end select

    end subroutine WriteItemTxt

    subroutine WriteArrayTxt(this, str, form, advance, number)
    class(TTextFile) :: this
    class(*), intent(in) :: str(:)
    character(LEN=*), intent(in), optional :: form
    logical, intent(in), optional :: advance
    integer, intent(in), optional :: number
    integer n
    character(LEN=3) :: Ad

    Ad = this%DefaultAdvance(advance)
    n = PresentDefault(size(str),number)
    call this%CheckOpen(.true.)
    select type(str)
    type is (character(LEN=*))
        write(this%unit,PresentDefault('(a)',form), advance=Ad) str(1:n)
    type is (real)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str(1:n)
    type is (double precision)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str(1:n)
    type is (integer)
        write(this%unit,PresentDefault(this%IntegerFormat,form), advance=Ad) str(1:n)
    type is (logical)
        write(this%unit,PresentDefault('(*(L2))',form), advance=Ad) str(1:n)
        class default
        call this%Error('unknown type to write')
    end select

    end subroutine WriteArrayTxt

    subroutine WriteArray2Txt(this, R)
    class(TTextFile) :: this
    class(*), intent(in) :: R(:,:)
    integer i

    do i=1, size(R,1)
        call this%WriteArrayTxt(R(i,:))
    end do

    end subroutine WriteArray2Txt


    subroutine WriteTrimTxt(this, S)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: S

    call this%WriteItemTxt(S, advance=.true.)

    end subroutine WriteTrimTxt

    subroutine WriteInLineTrim(this, string)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: string

    call this%WriteItemTxt(trim(string), advance=.false.)
    end subroutine WriteInLineTrim


    subroutine ReadArray2Txt(this, R, OK)
    class(TTextFile) :: this
    class(*) ::  R(:,:)
    logical, optional :: OK
    integer i

    do i=1, size(R,1)
        call this%ReadArrayTxt(R(i,:), OK=OK)
        if (present(OK)) then
            if (.not. OK) return
        end if
    end do

    end subroutine ReadArray2Txt


    subroutine ReadItemTxt(this, R, OK)
    class(TTextFile) :: this
    class(*), intent(out) :: R
    logical, optional :: OK
    integer status

    call this%CheckOpen
    select type(R)
    type is (character(LEN=*))
        Read(this%unit,'(a)', iostat=status) R
    type is (real)
        Read(this%unit,*, iostat=status) R
    type is (double precision)
        Read(this%unit,*, iostat=status) R
    type is (integer)
        Read(this%unit,*, iostat=status) R
        class default
        call this%Error('unknown type to Read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
        & call this%Error('Error reading item')
    if (present(OK)) OK = status==0

    end subroutine ReadItemTxt

    subroutine ReadArrayTxt(this, R, n, OK)
    class(TTextFile) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical, optional :: OK
    integer status

    call this%CheckOpen
    select type(R)
    type is (real)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
    type is (double precision)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
    type is (integer)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
        class default
        call this%Error('unknown type to Read')
    end select

    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
        & call this%Error('Error reading item')
    if (present(OK)) OK = status==0

    end subroutine ReadArrayTxt

    subroutine ReadStringTxt(this, S, OK)
    class(TTextFile) :: this
    character(LEN=:), allocatable :: S
    logical, optional :: OK
    logical isOK

    isOK = this%ReadLine(S)
    if (present(OK)) OK = isOK

    end subroutine  ReadStringTxt

    subroutine WriteFormat(this, formatst, i1,i2,i3,i4,i5,i6,i7,i8)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: formatst
    class(*), intent(in) :: i1
    class(*), intent(in),optional :: i2,i3,i4,i5,i6,i7,i8

    call this%CheckOpen(.true.)
    write(this%unit,'(a)') FormatString(formatst,i1,i2,i3,i4,i5,i6,i7,i8)

    end subroutine WriteFormat

    !Misc functions

    function TopCommentLine(aname) result(res)
    !Get top comment line in file, including #
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: res
    Type(TTextFile) :: F

    call F%Open(aname)
    res=''
    do while (res == '')
        if (.not. F%ReadLine(res)) exit
    end do
    If (res(1:1)/='#') then
        res = ''
    end if
    call F%Close()

    end function TopCommentLine


    function LastTopComment(aname) result(res)
    !Get content of last commented line at the top of file (e.g. column header), without #
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: res, InLine
    Type(TTextFile) :: F

    call F%Open(aname)
    res=''
    do while (F%ReadLine(InLine))
        if (trim(InLine)=='') cycle
        if (InLine(1:1)=='#') then
            res = trim(adjustl(InLine(2:)))
        else
            exit
        end if
    end do
    call F%Close()

    end function LastTopComment


    function TxtFileColumns(aname) result(n)
    character(LEN=*), intent(IN) :: aname
    integer n
    Type(TTextFile) :: F

    call F%Open(aname)
    n = F%Columns()
    call F%Close()

    end function TxtFileColumns

    function TxtFileLines(aname) result(n)
    character(LEN=*), intent(IN) :: aname
    integer n
    Type(TTextFile) :: F

    call F%Open(aname)
    n = F%Lines()
    call F%Close()

    end function TxtFileLines


    function LastFileLine(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: LastFileLine
    Type(TTextFile) :: F

    LastFileLine = ''
    call F%Open(aname)
    do while (F%ReadLine(LastFileLine))
    end do
    call F%Close()
    end function LastFileLine


    function CharIsSlash(C)
    character, intent(in) :: C
    logical CharIsSlash
    character, parameter :: win_slash = char(92)

    CharIsSlash = C == win_slash .or. C == '/'

    end function CharIsSlash

    function IsFullPath(aname)
    character(LEN=*), intent(IN) :: aname
    logical IsFullPath

    IsFullPath = aname/=''
    if (.not. IsFullPath) return
    IsFullPath = CharIsSlash(aname(1:1))

    end function IsFullpath

    function ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFilePath
    integer i

    do i = len_trim(aname), 1, -1
        if (CharIsSlash(aname(i:i))) then
            ExtractFilePath = aname(1:i)
            return
        end if
    end do
    ExtractFilePath = ''

    end function ExtractFilePath

    function ExtractFileExt(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFileExt
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (CharIsSlash(aname(i:i))) then
            ExtractFileExt = ''
            return
        else if (aname(i:i)=='.') then
            ExtractFileExt= aname(i:len)
            return
        end if
    end do
    ExtractFileExt = ''

    end function ExtractFileExt


    function ExtractFileName(aname, no_ext, all_ext)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFileName
    logical, intent(in), optional :: no_ext, all_ext
    integer alen, i

    alen = len_trim(aname)
    do i = alen, 1, -1
        if (CharIsSlash(aname(i:i))) then
            ExtractFileName = aname(i+1:alen)
            exit
        end if
    end do
    if (.not. allocated(ExtractFileName)) ExtractFileName = trim(aname)
    if (DefaultFalse(no_ext)) then
        do i = len(ExtractFileName), 1, -1
            if (ExtractFileName(i:i)=='.') then
                ExtractFileName = ExtractFileName(1:i-1)
                if (.not. DefaultFalse(all_ext)) exit
            end if
        end do
    end if

    end function ExtractFileName

    function ChangeFileExt(aname,ext)
    character(LEN=*), intent(IN) :: aname,ext
    character(LEN=:), allocatable :: ChangeFileExt
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='.') then
            ChangeFileExt = aname(1:i) // trim(ext)
            return
        end if
    end do
    ChangeFileExt = trim(aname) // '.' // trim(ext)

    end function ChangeFileExt

    function CheckTrailingSlash(aname)
    character(LEN=*), intent(in) :: aname
    character(LEN=:), allocatable :: CheckTrailingSlash
    integer alen

    alen = len_trim(aname)
    if (CharIsSlash(aname(alen:alen))) then
        CheckTrailingSlash = aname(1:alen)
    else
        CheckTrailingSlash = trim(aname)//'/'
    end if

    end function CheckTrailingSlash


    function File_Join(path, aname)
    character(LEN=*), intent(in) :: path, aname
    character(LEN=:), allocatable :: File_Join

    File_Join = CheckTrailingSlash(path)//trim(aname)

    end function File_Join


    subroutine DeleteFile(aname)
    character(LEN=*), intent(IN) :: aname
    integer file_id, status

    if (FileExists(aname)) then
        open(newunit = file_id, file = aname, iostat=status)
        if (status/=0) return
        close(unit = file_id, status = 'DELETE')
    end if

    end subroutine DeleteFile


    subroutine ReadTextVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    class(*), intent(out) :: vec(n)
    integer j
    Type(TTextFile) :: F

    call F%Open(aname)
    do j=1,n
        call F%Read(vec(j))
    end do
    call F%Close()

    end subroutine ReadTextVector

    subroutine WriteTextVector(aname, vec, n, fmt)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in), optional :: n
    class(*), intent(in) :: vec(:)
    character(LEN=*), intent(in), optional :: fmt
    integer j
    Type(TTextFile) :: F

    call F%CreateFile(aname)
    if (present(fmt)) then
        if (isFLoat(vec)) then
            F%RealFormat = fmt
        else
            F%IntegerFormat = fmt
        end if
    end if
    do j=1, PresentDefault(size(vec),n)
        call F%Write(vec(j))
    end do
    call F%Close()

    end subroutine WriteTextVector

    subroutine WriteTextMatrix(aname, mat, m, n, fmt)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in), optional :: m, n
    character(LEN=*), intent(in), optional :: fmt
    class(*), intent(in) :: mat(:,:)
    Type(TTextFile) :: F

    call F%CreateFile(aname)
    if (present(fmt)) then
        if (isFLoat(mat)) then
            F%RealFormat = fmt
        else
            F%IntegerFormat = fmt
        end if
    end if
    if (present(m) .or. present(n)) then
        call F%Write(mat(1:PresentDefault(size(mat,1),m),1:PresentDefault(size(mat,2),n)))
    else
        call F%Write(mat)
    end if
    call F%Close()

    end subroutine WriteTextMatrix


    subroutine ReadTextMatrix(aname, mat, inm,inn)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in), optional :: inm,inn
    real(kind(1.d0)), intent(out) :: mat(:,:)
    character(LEN=:), allocatable :: InLine
    integer j,k, status, n,m
    Type(TTextFile) :: F

    m = PresentDefault(size(mat,dim=1),inm)
    n = PresentDefault(size(mat,dim=2),inn)
    call F%Open(aname)

    do j=1,m
        status = 1
        if (.not. F%ReadLineSkipEmptyAndComments(InLine)) exit
        read (InLine,*, iostat=status) mat(j,1:n)
        if (status/=0) exit
    end do
    if (status/=0) then
        call F%Rewind()  !Try other possible format
        do j=1,m
            do k=1,n
                status = 1
                if (.not. F%ReadLineSkipEmptyAndComments(InLine)) exit
                read (InLine,*, iostat=status) mat(j,k)
                if (status/=0) call F%Error( 'matrix file is the wrong size')
            end do
        end do
    end if

    if (F%ReadLineSkipEmptyAndComments(InLine)) call F%Error( 'matrix file is the wrong size (too big)')

    call F%Close()

    end subroutine ReadTextMatrix

    subroutine LoadTxt_2D(aname, mat, m, n, comment)
    character(LEN=*), intent(IN) :: aname
    real(kind(1.d0)), allocatable :: mat(:,:)
    integer mm, nn, j
    Type(TTextFile) :: F
    character(LEN=:), allocatable :: InLine
    integer, optional, intent(out) :: m, n
    character(LEN=:), allocatable, optional, intent(out) :: comment
    integer status

    call F%Open(aname)
    nn = F%Columns()
    mm = F%Lines()
    allocate(mat(mm,nn))
    j=1
    do while (F%ReadLineSkipEmptyAndComments(InLine, comment = comment))
        read (InLine,*, iostat=status) mat(j,1:nn)
        if (status/=0) call F%Error( 'LoadTxt: error reading line:' //trim(InLine))
        j = j+1
    end do
    call F%Close()
    if (present(m)) m = mm
    if (present(n)) n = nn

    end subroutine LoadTxt_2D

    subroutine LoadTxt_1D(aname, vec, n, comment)
    character(LEN=*), intent(IN) :: aname
    real(kind(1.d0)), allocatable :: vec(:)
    integer nn, j
    Type(TTextFile) :: F
    character(LEN=:), allocatable :: InLine
    integer, optional, intent(out) :: n
    character(LEN=:), allocatable, optional, intent(out) :: comment
    integer status

    call F%Open(aname)
    nn = F%Lines()
    allocate(vec(nn))
    j=1
    do while (F%ReadLineSkipEmptyAndComments(InLine, comment = comment))
        read (InLine,*, iostat=status) vec(j)
        if (status/=0) call F%Error( 'LoadTxt: error reading line:' //trim(InLine))
        j = j+1
    end do
    call F%Close()
    if (present(n)) n = nn

    end subroutine LoadTxt_1D

    function CreateTextFile(fname)
    character(LEN=*), intent(in) :: fname
    Type(TTextFile) :: CreateTextFile

    call CreateTextFile%CreateFile(fname)

    end function CreateTextFile

    function OpenTextFile(fname)
    character(LEN=*), intent(in) :: fname
    Type(TTextFile) :: OpenTextFile

    call OpenTextFile%OpenFile(fname)

    end function OpenTextFile


    end module FileUtils
