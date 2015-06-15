    !A collection of ranges, consisting of sections of minimum step size
    !Useful for getting piecewise uniform point samplings for numerical integration
    !(with possible mix of log and linear spacing)
    !Antony Lewis, http://cosmologist.info/

    module RangeUtils
    use MiscUtils
    implicit none
    private

    type, public :: TRange
        integer :: start_index
        integer :: steps
        logical :: IsLog
        double precision :: Low, High
        double precision :: delta
        double precision :: delta_max, delta_min !for log spacing, the non-log max and min step size
    contains

    end type TRange

    type, public :: TRanges
        integer :: count = 0
        integer :: npoints = 0
        double precision :: Lowest, Highest
        type(TRange), allocatable :: R(:)
        logical :: has_dpoints = .false.
        double precision, dimension(:), allocatable :: points, dpoints
        !dpoints is (points(i+1)-points(i-1))/2
        double precision :: RangeTol = 0.1d0
        !fraction of bin width we are prepared for merged bin widths to increase by
    contains
    procedure :: Init => TRanges_Free
    procedure :: Free => TRanges_Free
    procedure :: IndexOf => TRanges_IndexOf
    procedure :: GetArray => TRanges_GetArray
    procedure :: Getdpoints => TRanges_Getdpoints
    procedure :: Add_delta => TRanges_Add_delta
    procedure :: Add => TRanges_Add
    procedure :: Write => TRanges_Write
    end type TRanges

    contains

#ifdef __GFORTRAN__
    impure elemental subroutine TRanges_Free(this)
#else
    subroutine TRanges_Free(this)
#endif
    class(TRanges), intent(inout) :: this

    if (allocated(this%R)) deallocate(this%R)
    if (allocated(this%points)) deallocate(this%points)
    if (allocated(this%dpoints)) deallocate(this%dpoints)
    this%count = 0
    this%npoints = 0
    this%has_dpoints = .false.
    end subroutine TRanges_Free


    function TRanges_IndexOf(this, tau) result(pointstep)
    class(TRanges), intent(in) :: this
    double precision, intent(in) :: tau
    integer :: pointstep, i

    pointstep=0
    do i=1, this%count
        associate(AReg => this%R(i))
            if (tau < AReg%High .and. tau >= AReg%Low) then
                if (AReg%IsLog) then
                    pointstep = AReg%start_index + int(log(tau / AReg%Low) / AReg%delta)
                else
                    pointstep = AReg%start_index + int((tau - AReg%Low) / AReg%delta)
                end if
                return
            end if
        end associate
    end do

    if (tau >= this%Highest) then
        pointstep = this%npoints
    else
        stop 'TRanges_IndexOf: value out of range'
    end if
    end function TRanges_IndexOf


    subroutine TRanges_GetArray(this, want_dpoints)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: want_dpoints
    integer :: i,j,ix

    this%has_dpoints = DefaultTrue(want_dpoints)

    ! Dealloc/Realloc only, when not enough space present.
    if (allocated(this%points) .and. size(this%points) < this%npoints) &
        deallocate(this%points)
    if (.not. allocated(this%points)) allocate(this%points(this%npoints))

    ix=0
    do i=1, this%count
        associate (AReg => this%R(i))
            do j = 0, AReg%steps-1
                ix=ix+1
                if (AReg%IsLog) then
                    this%points(ix) = AReg%Low*exp(j*AReg%delta)
                else
                    this%points(ix) = AReg%Low + AReg%delta*j
                end if
            end do
        end associate
    end do
    ix =ix+1
    this%points(ix) = this%Highest
    if (ix /= this%npoints) stop 'TRanges_GetArray: ERROR'

    if (this%has_dpoints) call this%Getdpoints()
    end subroutine TRanges_GetArray


    subroutine TRanges_Getdpoints(this, half_ends)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: half_ends
    integer :: i
    logical :: halfs

    halfs = DefaultTrue(half_ends)

    ! Dealloc/Realloc only, when not enough space present.
    if (allocated(this%dpoints) .and. size(this%dpoints) < this%npoints) &
        deallocate(this%dpoints)
    if (.not. allocated(this%dpoints)) allocate(this%dpoints(this%npoints))

    do i=2, this%npoints-1
        this%dpoints(i) = (this%points(i+1) - this%points(i-1))/2
    end do
    if (halfs) then
        this%dpoints(1) = (this%points(2) - this%points(1))/2
        this%dpoints(this%npoints) = (this%points(this%npoints) - this%points(this%npoints-1))/2
    else
        this%dpoints(1) = (this%points(2) - this%points(1))
        this%dpoints(this%npoints) = (this%points(this%npoints) - this%points(this%npoints-1))
    end if
    end subroutine TRanges_Getdpoints


    subroutine TRanges_Add_delta(this, t_start, t_end, t_approx_delta, IsLog)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: IsLog
    double precision, intent(in) :: t_start, t_end, t_approx_delta
    integer :: n
    logical :: WantLog

    WantLog = DefaultFalse(IsLog)

    if (t_end <= t_start) &
        stop 'TRanges_Add_delta: end must be larger than start'
    if (t_approx_delta <=0) stop 'TRanges_Add_delta: delta must be > 0'

    if (WantLog) then
        n  = max(1,int(log(t_end/t_start)/t_approx_delta + 1.d0 - this%RangeTol))
    else
        n  = max(1,int((t_end-t_start)/t_approx_delta + 1.d0 - this%RangeTol))
    end if
    call this%Add(t_start, t_end, n, WantLog)
    end subroutine TRanges_Add_delta


    subroutine TRanges_Add(this, t_start, t_end, nstep, IsLog)
    class(TRanges), intent(inout) :: this
    logical, intent(in), optional :: IsLog
    double precision, intent(in) :: t_start, t_end
    integer, intent(in) :: nstep
    type(TRange), allocatable, target :: NewRanges(:)
    double precision, allocatable :: EndPoints(:), RequestDelta(:)
    integer :: ixin, nreg, ix, i,j, nsteps
    double precision :: min_request, max_request, min_log_step, max_log_step
    double precision :: delta, diff, max_delta
    logical :: WantLog

    WantLog = DefaultFalse(IsLog)

    if (WantLog) then
        delta = log(t_end / t_start) / nstep
    else
        delta = (t_end - t_start) / nstep
    end if

    if (t_end <= t_start) stop 'TRanges_Add: end must be larger than start'
    if (nstep <= 0) stop 'TRanges_Add: nstep must be > 0'

    nreg = this%count + 1
    allocate(NewRanges(nreg))
    if (allocated(this%R) .and. this%count > 0) then
#if defined(__IBMCPP__) || defined(__xlC__) || defined(__xlc__)
        !           avoid IBM compiler bug, from Angel de Vicente
        !           detection of IBM compiler by preprocessors symbols taken from boost
        do i= 1, this%count
            NewRanges(i) = this%R(i)
        end do
#else
        NewRanges(1:this%count) = this%R(1:this%count)
#endif
    end if

    allocate(EndPoints(0:nreg * 2))
    associate (AReg => NewRanges(nreg))
        AReg%Low   = t_start
        AReg%High  = t_end
        AReg%delta = delta
        AReg%steps = nstep
        AReg%IsLog = WantLog
    end associate

    !Get end point in order
    ix = 0
    do i=1, nreg
        associate (AReg => NewRanges(i))
            if (ix==0) then
                EndPoints(1) = AReg%Low
                EndPoints(2) = AReg%High
                ix = 2
            else
                ixin = ix
                do j=1,ixin
                    if (AReg%Low < EndPoints(j)) then
                        EndPoints(j+1:ix+1) = EndPoints(j:ix)
                        EndPoints(j) = AReg%Low
                        ix=ix+1
                        exit
                    end if
                end do
                if (ixin == ix) then
                    ix = ix+1
                    EndPoints(ix) = AReg%Low
                    ix = ix+1
                    EndPoints(ix) = AReg%High
                else
                    ixin = ix
                    do j=1,ixin
                        if (AReg%High < EndPoints(j)) then
                            EndPoints(j+1:ix+1) = EndPoints(j:ix)
                            EndPoints(j) = AReg%High
                            ix=ix+1
                            exit
                        end if
                    end do
                    if (ixin == ix) then
                        ix = ix+1
                        EndPoints(ix) = AReg%High
                    end if
                end if
            end if
        end associate
    end do

    !remove duplicate points
    ixin = ix
    ix = 1
    do i=2, ixin
        if (EndPoints(i) /= EndPoints(ix)) then
            ix=ix+1
            EndPoints(ix) = EndPoints(i)
        end if
    end do


    !ix is the number of end points
    this%Lowest = EndPoints(1)
    this%Highest = EndPoints(ix)
    this%count = 0

    max_delta = this%Highest - this%Lowest

    if (.not. allocated(this%R) .or. size(this%R) < ix + 1) then
        if (allocated (this%R)) deallocate (this%R)
        allocate (this%R(ix + 1))
    end if

    allocate(RequestDelta(ix))

    do i=1, ix - 1
        associate (AReg => this%R(i))
            AReg%Low = EndPoints(i)
            AReg%High = EndPoints(i+1)

            !               max_delta = EndPoints(i+1) - EndPoints(i)
            delta = max_delta
            AReg%IsLog = .false.

            do j=1, nreg
                if (AReg%Low >= NewRanges(j)%Low .and. AReg%Low < NewRanges(j)%High) then
                    if (NewRanges(j)%IsLog) then
                        if (AReg%IsLog) then
                            delta = min(delta,NewRanges(j)%delta)
                        else
                            min_log_step = AReg%Low*(exp(NewRanges(j)%delta)-1)
                            if (min_log_step < delta) then
                                max_log_step = AReg%High*(1-exp(-NewRanges(j)%delta))
                                if  (delta < max_log_step) then
                                    delta = min_log_step
                                else
                                    AReg%IsLog = .true.
                                    delta = NewRanges(j)%delta
                                end if
                            end if
                        end if
                    else !New Range is not log
                        if (AReg%IsLog) then
                            max_log_step = AReg%High*(1-exp(-delta))
                            if (NewRanges(j)%delta < max_log_step) then
                                min_log_step = AReg%Low*(exp(delta)-1)
                                if (min_log_step <  NewRanges(j)%delta) then
                                    AReg%IsLog = .false.
                                    delta =  min_log_step
                                else
                                    delta = - log(1- NewRanges(j)%delta/AReg%High)
                                end if
                            end if
                        else
                            delta = min(delta, NewRanges(j)%delta)
                        end if
                    end if
                end if
            end do

            if (AReg%IsLog) then
                Diff = log(AReg%High/AReg%Low)
            else
                Diff = AReg%High - AReg%Low
            endif
            if (delta >= Diff) then
                AReg%delta = Diff
                AReg%steps = 1
            else
                AReg%steps  = max(1,int(Diff/delta + 1.d0 - this%RangeTol))
                AReg%delta = Diff / AReg%steps
            end if

            this%count = this%count + 1
            RequestDelta(this%count) = delta

            if (AReg%IsLog) then
                if (AReg%steps ==1) then
                    AReg%Delta_min = AReg%High - AReg%Low
                    AReg%Delta_max = AReg%Delta_min
                else
                    AReg%Delta_min = AReg%Low*(exp(AReg%delta)-1)
                    AReg%Delta_max = AReg%High*(1-exp(-AReg%delta))
                end if
            else
                AReg%Delta_max = AReg%delta
                AReg%Delta_min = AReg%delta
            end if
        end associate
    end do


    !Get rid of tiny TRanges
    ix = this%count
    do i=ix, 1, -1
        associate (AReg => this%R(i))
            if (AReg%steps ==1) then
                Diff = AReg%High - AReg%Low
                if (AReg%IsLog) then
                    min_request = AReg%Low*(exp(RequestDelta(i))-1)
                    max_request = AReg%High*(1-exp(-RequestDelta(i)))
                else
                    min_request = RequestDelta(i)
                    max_request = min_request
                end if
                if (i/= this%count) then  !from i/= ix Mar08
                    associate (LastReg => this%R(i+1))
                        if (RequestDelta(i) >= AReg%delta .and. Diff <= LastReg%Delta_min &
                            .and. LastReg%Delta_min <= max_request) then

                        LastReg%Low = AReg%Low
                        if (Diff > LastReg%Delta_min*this%RangeTol) then
                            LastReg%steps =  LastReg%steps + 1
                        end if
                        if (LastReg%IsLog) then
                            LastReg%delta = log(LastReg%High/LastReg%Low) / LastReg%steps
                        else
                            LastReg%delta = (LastReg%High -LastReg%Low) / LastReg%steps
                        end if
                        this%R(i:this%Count-1) = this%R(i+1:this%Count)
                        this%Count = this%Count -1
                        cycle
                        end if
                    end associate
                end if
                if (i/=1) then
                    associate (LastReg => this%R(i-1))
                        if (RequestDelta(i) >= AReg%delta .and. Diff <= LastReg%Delta_max &
                            .and. LastReg%Delta_max <= min_request) then
                        LastReg%High = AReg%High
                        !AlMat08 LastReg%Low = AReg%Low
                        if (Diff > LastReg%Delta_max*this%RangeTol) then
                            LastReg%steps =  LastReg%steps + 1
                        end if
                        if (LastReg%IsLog) then
                            LastReg%delta = log(LastReg%High/LastReg%Low) / LastReg%steps
                        else
                            LastReg%delta = (LastReg%High -LastReg%Low) / LastReg%steps
                        end if
                        this%R(i:this%Count-1) = this%R(i+1:this%Count)
                        this%Count = this%Count -1
                        end if
                    end associate
                end if
            end if
        end associate
    end do


    !Set up start indices and get total number of steps
    nsteps = 1
    do i = 1, this%Count
        associate (AReg => this%R(i))
            AReg%Start_index = nsteps
            nsteps = nsteps + AReg%steps
            if (AReg%IsLog) then
                if (AReg%steps ==1) then
                    AReg%Delta_min = AReg%High - AReg%Low
                    AReg%Delta_max = AReg%Delta_min
                else
                    AReg%Delta_min = AReg%Low*(exp(AReg%delta)-1)
                    AReg%Delta_max = AReg%High*(1-exp(-AReg%delta))
                end if
            else
                AReg%Delta_max = AReg%delta
                AReg%Delta_min = AReg%delta
            end if
        end associate
    end do

    this%npoints = nsteps

    deallocate(NewRanges, EndPoints, RequestDelta)
    end subroutine TRanges_Add


    subroutine TRanges_Write(this)
    class(TRanges), intent(in), target :: this
    integer :: i

    do i=1,this%count
        associate (AReg => this%R(i))
            if (AReg%IsLog) then
                Write (*,'("Range ",I3,":", 3E14.4," log")') i, AReg%Low, AReg%High, AReg%delta
            else
                Write (*,'("Range ",I3,":", 3E14.4," linear")') i, AReg%Low, AReg%High, AReg%delta
            end if
        end associate
    end do
    end subroutine TRanges_Write


    end module RangeUtils
