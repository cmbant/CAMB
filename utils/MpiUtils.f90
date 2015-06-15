    module MpiUtils
    implicit none

#ifdef MPI
    include "mpif.h"
#endif

    integer, parameter :: TTimer_dp = Kind(1.d0)

    Type TTimer
        real(TTimer_dp) start_time
    contains
    procedure :: Start => TTimer_Start
    procedure :: Time => TTimer_Time
    procedure :: WriteTime => TTimer_WriteTime
    end type TTimer

    contains

    function GetMpiRank()
    integer GetMpiRank
#ifdef MPI
    integer ierror
    call mpi_comm_rank(mpi_comm_world,GetMPIrank,ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MPI fail')
#else
    GetMpiRank=0
#endif

    end function GetMpiRank

    function IsMainMPI()
    logical IsMainMPI

    IsMainMPI =  GetMpiRank() == 0

    end function IsMainMPI

    subroutine MpiStop(Msg)
    character(LEN=*), intent(in), optional :: Msg
    integer i
#ifdef MPI
    integer ierror, MpiRank
#endif

    if (present(Msg)) write(*,*) trim(Msg)

#ifdef MPI
    call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)
    write (*,*) 'MpiStop: ', MpiRank
    call MPI_ABORT(MPI_COMM_WORLD,i)
#endif
    i=1     !put breakpoint on this line to debug
    stop

    end subroutine MpiStop

    subroutine MpiStat(MpiID, MpiSize)
    implicit none
    integer MpiID,MpiSize
#ifdef MPI
    integer ierror
    call mpi_comm_rank(mpi_comm_world,MpiID,ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MpiStat: MPI rank')
    call mpi_comm_size(mpi_comm_world,MpiSize,ierror)
#else
    MpiID=0
    MpiSize=1
#endif
    end subroutine MpiStat

    subroutine MpiQuietWait
    !Set MPI thread to sleep, e.g. so can run openmp on cpu instead
#ifdef MPI
    integer flag, ierr, STATUS(MPI_STATUS_SIZE)
    integer i, MpiId, MpiSize

    call MpiStat(MpiID, MpiSize)
    if (MpiID/=0) then
        do
            call MPI_IPROBE(0,0,MPI_COMM_WORLD,flag, MPI_STATUS_IGNORE,ierr)
            if (flag/=0) then
                call MPI_RECV(i,1,MPI_INTEGER, 0,0,MPI_COMM_WORLD,status,ierr)
                exit
            end if
            call sleep(1)
        end do
    end if
#endif
    end subroutine

    subroutine MpiWakeQuietWait
#ifdef MPI
    integer j, MpiId, MpiSize, ierr,r

    call MpiStat(MpiID, MpiSize)
    if (MpiID==0) then
        do j=1, MpiSize-1
            call MPI_ISSEND(MpiId,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,r,ierr)
        end do
    end if
#endif
    end subroutine MpiWakeQuietWait

    subroutine MpiShareString(S, from)
    character(LEN=:), allocatable :: S
    integer from
#ifdef MPI
    integer inlen, rank, ierror

    rank = GetMpiRank()

    if (rank==from) inlen=len(S)

    CALL MPI_Bcast(inlen, 1, MPI_INTEGER, from, MPI_COMM_WORLD, ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MpiShareString: fail')

    if (rank /= from ) allocate(character(inlen)::S)
    CALL MPI_Bcast(S, LEN(S), MPI_CHARACTER, from, MPI_COMM_WORLD, ierror)
#endif
    end subroutine MpiShareString


    function TimerTime()
    real(TTimer_dp) time
    real(TTimer_dp) :: TimerTime
#ifdef MPI
    TimerTime = MPI_WTime()
#else
    call cpu_time(time)
    TimerTime=  time
#endif
    end function TimerTime

    subroutine TTimer_Start(this)
    class(TTimer) :: this
    this%start_time = TimerTime()
    end subroutine TTimer_Start

    real(TTimer_dp) function TTimer_Time(this)
    class(TTimer) :: this
    TTimer_Time =  TimerTime() - this%start_time
    end function TTimer_Time

    subroutine TTimer_WriteTime(this,Msg, start)
    class(TTimer) :: this
    character(LEN=*), intent(in), optional :: Msg
    real(TTimer_dp), optional :: start
    real(TTimer_dp) T, DeltaT
    character(LEN=:), allocatable :: tmp

    if (present(start)) then
        T=start
    else
        T=this%start_time
    end if

    DeltaT = TimerTime() - T
    if (present(Msg)) then
        tmp = trim(Msg)//': '
        if (DeltaT > 0.00002 .and. DeltaT < 1000 .and. len_trim(tmp)<24) then
            write (*,'(a25,f10.5)') tmp, DeltaT
        else
            write (*,*) trim(Msg)//': ', DeltaT
        end if
    end if
    if (.not. present(start)) this%start_time = TimerTime()

    end subroutine TTimer_WriteTime


    end module MpiUtils
