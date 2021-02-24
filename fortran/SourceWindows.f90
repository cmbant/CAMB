    module SourceWindows
    use precision
    use Classes
    use MpiUtils
    implicit none

    integer, parameter :: window_21cm = 1, window_counts = 2, window_lensing = 3

    Type, extends(TPythonInterfacedClass) :: TSourceWindow
        integer :: source_type = window_counts
        real(dl) :: bias = 1._dl
        real(dl) :: dlog10Ndm = 0._dl
    contains
    procedure :: count_obs_window_z
    procedure :: Window_f_a
    procedure :: counts_background_z
    procedure :: GetScales
    procedure :: GetBias
    end Type TSourceWindow

    Type, extends(TSourceWindow) :: TGaussianSourceWindow
        real(dl) :: redshift
        real(dl) :: sigma !for 21cm, width in a, otherwise width in z
    contains
    procedure, nopass :: SelfPointer => TGaussianSourceWindow_SelfPointer
    procedure :: count_obs_window_z => TGaussianSourceWindow_count_obs_window_z
    procedure :: Window_f_a => TGaussianSourceWindow_Window_f_a
    procedure :: GetScales => TGaussianSourceWindow_GetScales
    end Type TGaussianSourceWindow

    Type, extends(TSourceWindow) :: TSplinedSourceWindow
        Type(TCubicSpline), allocatable :: Window
        Type(TCubicSpline), allocatable :: Bias_z
        real(dl) :: maxwin
    contains
    procedure, nopass :: SelfPointer => TSplinedSourceWindow_SelfPointer
    procedure :: count_obs_window_z => TSplinedSourceWindow_count_obs_window_z
    procedure :: GetScales => TSplinedSourceWindow_GetScales
    procedure :: SetTable => TSplinedSourceWindow_SetTable
    procedure :: GetBias => TSplinedSourceWindow_GetBias
    end Type TSplinedSourceWindow

    Type TSourceWindowHolder
        class(TSourceWindow), allocatable :: Window
    end Type TSourceWindowHolder

    Type SourceTermParams
        logical :: limber_windows = .true.
        integer :: limber_phi_lmin = 100
        !for l>limber_phi use limber approx for lensing potential when sourceTerms%LimberWindows=True;
        !Limber is also used if LimberWindows=False but via the time integrals instead.
        logical :: counts_density = .true.
        logical :: counts_redshift = .true.
        logical :: counts_lensing = .false.
        logical :: counts_velocity  = .true.
        logical :: counts_radial = .false. !does not include time delay; subset of counts_velocity, just 1/(chi*H) term
        logical :: counts_timedelay = .true. !time delay terms * 1/(H*chi)
        logical :: counts_ISW = .true.
        logical :: counts_potential = .true. !terms in potentials at source
        logical :: counts_evolve = .false.
        logical :: line_phot_dipole = .false.
        logical :: line_phot_quadrupole= .false.
        logical :: line_basic = .true.
        logical :: line_distortions = .true.
        logical :: line_extra = .false.
        logical :: line_reionization = .false.
        logical :: use_21cm_mK = .true.
    end type SourceTermParams

    Type TRedWin !internal type
        class(TSourceWindow), pointer :: Window => null()
        integer kind
        real(dl) Redshift
        real(dl) tau, tau_start, tau_end, tau_peakstart, tau_peakend
        real(dl) sigma_tau !approx width in conformal time (set by code)
        real(dl) chi0, chimin
        integer :: mag_index =0 !The index into the extra sources used for adding magnification to counts
        real(dl), dimension(:), allocatable :: winF, wing,wing2,wingtau,dwing,dwing2,dwingtau,ddwing,ddwing2,ddwingtau,&
            winV,dwinV,ddwinV, win_lens, comoving_density_ev
        real(dl) Fq, optical_depth_21
        logical has_lensing_window
    end Type TRedWin

    contains

    function count_obs_window_z(this, z, winamp)
    !distribution function W(z) for the observed sources, used for lensing and number count spectrum
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    !note this is the total count distribution observed, not a fractional selection function on an underlying distribution
    class(TSourceWindow) :: this
    real(dl), intent(in) :: z
    real(dl) count_obs_window_z, winamp
    count_obs_window_z =0
    end function count_obs_window_z


    function counts_background_z(this, z)
    !if counts_evolve = T this function is used to get n(z) for the source population
    !(not the same as the number actually observed)
    class(TSourceWindow) :: this
    real(dl), intent(in) :: z
    real(dl) counts_background_z, winamp

    counts_background_z = this%count_obs_window_z(z, winamp)
    !This is the special case where you observe all of the sources

    end function counts_background_z

    subroutine GetScales(this, zpeak, sigma_z, zpeakstart, zpeakend)
    class(TSourceWindow) :: this
    real(dl), intent(out) :: zpeak, sigma_z, zpeakstart, zpeakend

    zpeak=0
    sigma_z=0
    zpeakstart=0
    zpeakend=0
    call MpiStop('Must define GetScales function')

    end subroutine GetScales

    real(dl) function GetBias(this,k,a)
    class(TSourceWindow) :: this
    real(dl), intent(in) :: k,a
    GetBias = this%Bias !Simplest scale-independent and time independent model
    end function

    function Window_f_a(this, a, winamp)
    !distribution function as function of scale factor a
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    class(TSourceWindow) :: this
    real(dl), intent(in) :: a
    real(dl) Window_f_a, winamp

    if (this%source_type == window_21cm) then
        call MpiStop('Must define Window_f_a function for 21cm')
    else
        Window_f_a = this%count_obs_window_z(1/a-1,winamp)/a**2
    end if

    end function Window_f_a

    subroutine TGaussianSourceWindow_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TGaussianSourceWindow), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TGaussianSourceWindow_SelfPointer

    function TGaussianSourceWindow_count_obs_window_z(this, z, winamp)
    !distribution function W(z) for the observed sources, used for lensing and number count spectrum
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    !note this is the total count distribution observed, not a fractional selection function on an underlying distribution
    class(TGaussianSourceWindow) :: this
    real(dl), intent(in) :: z
    real(dl) TGaussianSourceWindow_count_obs_window_z, dz,winamp
    real(dl), parameter :: root2pi = 2.506628274_dl

    dz = z-this%Redshift
    winamp =  exp(-(dz/this%sigma)**2/2)
    TGaussianSourceWindow_count_obs_window_z =winamp/this%sigma/root2pi

    end function TGaussianSourceWindow_count_obs_window_z

    function TGaussianSourceWindow_Window_f_a(this, a, winamp)
    !distribution function as function of scale factor a
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    class(TGaussianSourceWindow) :: this
    real(dl), intent(in) :: a
    real(dl) TGaussianSourceWindow_Window_f_a, winamp
    real(dl), parameter :: root2pi = 2.506628274_dl

    if (this%source_type == window_21cm) then
        !Take W_T(S) = W_f(S)/S to be Gaussain, for frequency integration of T_b
        winamp =  exp(-((a-(1/(this%redshift+1)))/this%sigma)**2/2)
        TGaussianSourceWindow_Window_f_a = a*winamp/this%sigma/root2pi
    else
        TGaussianSourceWindow_Window_f_a = this%count_obs_window_z(1/a-1,winamp)/a**2
    end if

    end function TGaussianSourceWindow_Window_f_a


    subroutine TGaussianSourceWindow_GetScales(this, zpeak, sigma_z,  zpeakstart, zpeakend)
    class(TGaussianSourceWindow) :: this
    real(dl), intent(out) :: zpeak, sigma_z, zpeakstart, zpeakend

    if (this%source_type == Window_21cm) then
        sigma_z = this%sigma* (1 + this%RedShift) **2
    else
        sigma_z = this%sigma
    end if
    zpeak = this%redshift
    zpeakstart = this%redshift + sigma_z*3
    zpeakend = max(0._dl,this%redshift - sigma_z*3)

    end subroutine TGaussianSourceWindow_GetScales


    subroutine TSplinedSourceWindow_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TSplinedSourceWindow), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSplinedSourceWindow_SelfPointer

    real(dl) function TSplinedSourceWindow_GetBias(this,k,a)
    class(TSplinedSourceWindow) :: this
    real(dl), intent(in) :: k,a
    real(dl) z
    if (allocated(this%bias_z))  then
        z = 1/a-1
        if (z > this%Window%X(this%Window%n) .or. z < this%Window%X(1)) then
            TSplinedSourceWindow_GetBias = 0
        else
            TSplinedSourceWindow_GetBias = this%bias_z%value(z)
        end if
    else
        TSplinedSourceWindow_GetBias = this%Bias !Simplest scale-independent and time independent model
    end if
    end function


    subroutine  TSplinedSourceWindow_SetTable(this, n, z, W, bias_z)
    class(TSplinedSourceWindow) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: z(n), W(n)
    real(dl), intent(in), optional :: bias_z(n)

    if (allocated(this%Window)) deallocate(this%Window)
    if (n>0) then
        allocate(this%Window)
        call this%Window%Init(z,W)
        this%maxwin = maxval(this%Window%F)
    end if
    if (present(bias_z)) then
        if (allocated(this%Bias_z)) deallocate(this%Bias_z)
        if (n>0) then
            allocate(this%Bias_z)
            call this%Bias_z%Init(z,bias_z)
        end if
    end if

    end subroutine TSplinedSourceWindow_SetTable

    function TSplinedSourceWindow_count_obs_window_z(this, z, winamp)
    !distribution function W(z) for the observed sources, used for lensing and number count spectrum
    !Winamp is amplitude normalized to 1 so the code can tell when the window is very small
    !note this is the total count distribution observed, not a fractional selection function on an underlying distribution
    class(TSplinedSourceWindow) :: this
    real(dl), intent(in) :: z
    real(dl) TSplinedSourceWindow_count_obs_window_z, dz,winamp

    if (z > this%Window%X(this%Window%n) .or. z < this%Window%X(1)) then
        TSplinedSourceWindow_count_obs_window_z =0
    else
        TSplinedSourceWindow_count_obs_window_z = this%Window%Value(z)
    end if
    winamp = TSplinedSourceWindow_count_obs_window_z/this%maxwin

    end function TSplinedSourceWindow_count_obs_window_z

    subroutine TSplinedSourceWindow_GetScales(this, zpeak, sigma_z, zpeakstart, zpeakend)
    class(TSplinedSourceWindow) :: this
    real(dl), intent(out) :: zpeak, sigma_z, zpeakstart, zpeakend
    integer ix, i, j
    real(dl) z1, zstart, zend, targ

    associate(z => this%Window%X, W => this%Window%F, n=>this%Window%n)
        zstart = z(n)
        zend = z(1)
        zpeak = z(maxloc(W, dim=1))
        zpeakstart = zstart
        do i=n,1,-1
            if (W(i) > this%maxwin/15) then
                zpeakstart = z(i)
                exit
            end if
        end do
        zpeakend = zend
        do i=1,n
            if (W(i) > this%maxwin/15) then
                if (zpeakstart > z(i)) zpeakend = z(i)
                exit
            end if
        end do
        z1 = 0._dl
        !sigma_z sets the scale of structure in the window function
        !Used for determining when Limber should be valid and required time step size
        !Try minimum of some simple heuristics
        sigma_z = (zpeakstart - zpeakend)/3
        do i = 1, n
            if (W(i) > this%maxwin/2 .and. z1==0._dl) then
                z1 = z(i)
            else if (W(i) < this%maxwin/2 .and. z1/= 0._dl) then
                sigma_z = min(sigma_z, z(i)-z1)
                z1 = 0._dl
            end if
        end do
        do i = 1, n
            if (W(i) < this%maxwin*0.45) then
                targ = w(i) + this%maxwin/2
                do j=i+1, n
                    if (W(j) >= targ) then
                        sigma_z = min(sigma_z, (z(j)-z(i))/0.85)
                        exit
                    end if
                end do
                do j=i-1, 1,-1
                    if (W(j) >= targ) then
                        sigma_z = min(sigma_z, (z(i)-z(j))/0.85)
                        exit
                    end if
                end do
            end if
        end do
    end associate

    end subroutine TSplinedSourceWindow_GetScales

    end module SourceWindows