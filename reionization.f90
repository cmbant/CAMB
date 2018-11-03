
    module Reionization
    use Precision
    use MiscUtils
    use classes
    implicit none

    !This module puts smooth tanh reionization of specified mid-point (z_{re}) and width
    !The tanh function is in the variable (1+z)**Rionization_zexp

    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth

    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    character(LEN=*), parameter :: Reionization_Name = 'CAMB_reionization'

    real(dl), parameter :: Reionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other

    real(dl) :: Reionization_AccuracyBoost = 1._dl
    real(dl) :: Rionization_zexp = 1.5_dl

    logical :: include_helium_fullreion = .true.

    type, extends(TCambComponent) :: ReionizationParams
        logical    :: Reionization = .true.
        logical    :: use_optical_depth = .false.
        real(dl)   :: redshift = 10._dl
        real(dl)   :: delta_redshift = 0.5_dl
        real(dl)   :: fraction = Reionization_DefFraction
        real(dl)   :: optical_depth = 0._dl
        !Parameters for the second reionization of Helium
        real(dl)   :: helium_redshift  = 3.5_dl
        real(dl)   :: helium_delta_redshift  = 0.5
        real(dl)   :: helium_redshiftstart  = 5._dl
    contains
    procedure :: ReadParams => Reionization_ReadParams
    procedure :: Validate => Reionization_Validate
    end type ReionizationParams

    type ReionizationHistory
        !These two are used by main code to bound region where xe changing
        real(dl) :: tau_start, tau_complete
        !This is set from main code
        real(dl) :: akthom, fHe

        !The rest are internal to this module.
        real(dl) :: WindowVarMid, WindowVarDelta
        class(ReionizationParams), pointer :: Params
    contains 
    procedure :: Init => Reionization_Init
    procedure :: x_e => Reionization_xe
    procedure :: timesteps => Reionization_timesteps
    procedure, private :: SetParamsForZre => Reionization_SetParamsForZre
    procedure, private :: SetFromOptDepth => Reionization_SetFromOptDepth
    procedure, private :: zreFromOptDepth => Reionization_zreFromOptDepth
    procedure :: GetOptDepth => Reionization_GetOptDepth
    end type ReionizationHistory

    real(dl), parameter :: Reionization_maxz = 50._dl
    real(dl), private, parameter :: Reionization_tol = 1d-5

    real(dl), private, external :: dtauda, rombint
    contains


    function Reionization_xe(this,a, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(ReionizationHistory) :: this
    real(dl), intent(in) :: a
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) Reionization_xe
    real(dl) tgh, xod
    real(dl) xstart

    xstart = PresentDefault( 0._dl, xe_recomb)

    xod = (this%WindowVarMid - 1._dl/a**Rionization_zexp)/this%WindowVarDelta
    if (xod > 100) then
        tgh=1.d0
    else
        tgh=tanh(xod)
    end if
    Reionization_xe =(this%Params%fraction-xstart)*(tgh+1._dl)/2._dl+xstart

    if (include_helium_fullreion .and. a > (1/(1+ this%Params%helium_redshiftstart))) then

        !Effect of Helium becoming fully ionized is small so details not important
        xod = (1+this%Params%helium_redshift - 1._dl/a)/this%Params%helium_delta_redshift
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if

        Reionization_xe =  Reionization_xe + this%fHe*(tgh+1._dl)/2._dl

    end if

    end function Reionization_xe

    function Reionization_timesteps(this)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(ReionizationHistory) :: this
    integer Reionization_timesteps

    Reionization_timesteps = 50

    end  function Reionization_timesteps

    subroutine Reionization_ReadParams(this, Ini)
    use IniObjects
    class(ReionizationParams) :: this
    class(TIniFile), intent(in) :: Ini

    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

        this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')

        if (this%use_optical_depth) then
            this%optical_depth = Ini%Read_Double('re_optical_depth')
        else
            this%redshift = Ini%Read_Double('re_redshift')
        end if

        this%delta_redshift = Ini%Read_Double('re_delta_redshift', 0.5_dl) !default similar to CMBFAST original
        this%fraction = Ini%Read_Double('re_ionization_frac',Reionization_DefFraction)

        this%helium_redshift  = Ini%Read_Double('re_helium_redshift', 3.5_dl)
        this%helium_delta_redshift  = Ini%Read_Double('re_helium_delta_redshift', 0.5_dl)
        this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 3*this%helium_delta_redshift)

    end if

    end subroutine Reionization_ReadParams

    subroutine Reionization_SetParamsForZre(this, Reion)
    class(ReionizationHistory) :: this
    class(ReionizationParams) :: Reion

    this%WindowVarMid = (1._dl+Reion%redshift)**Rionization_zexp
    this%WindowVarDelta = &
        Rionization_zexp*(1._dl+Reion%redshift)**(Rionization_zexp-1._dl)*Reion%delta_redshift

    end subroutine Reionization_SetParamsForZre

    subroutine Reionization_Init(this, Reion, Yhe, akthom, tau0, FeedbackLevel)
    use constants
    use errors
    class(ReionizationHistory) :: this
    class(ReionizationParams), target :: Reion
    real(dl), intent(in) :: akthom, tau0, Yhe
    integer, intent(in) :: FeedbackLevel
    real(dl) astart

    this%Params => Reion
    this%akthom = akthom
    this%fHe =  YHe/(mass_ratio_He_H*(1.d0-YHe))

    this%tau_start=tau0
    this%tau_complete=tau0

    if (Reion%Reionization) then

        if (Reion%optical_depth /= 0._dl .and. .not. Reion%use_optical_depth) &
            write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'

        if (Reion%use_optical_depth.and.Reion%optical_depth<0.001 &
            .or. .not.Reion%use_optical_depth .and. Reion%Redshift<0.001) then
            Reion%Reionization = .false.
        end if

    end if

    if (Reion%Reionization) then

        if (Reion%fraction==Reionization_DefFraction) &
            Reion%fraction = 1._dl + this%fHe  !H + singly ionized He

        if (Reion%use_optical_depth) then
            call this%SetFromOptDepth(Reion)
            if (global_error_flag/=0) return
            if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') Reion%redshift
        end if

        call this%SetParamsForZre(Reion)

        !this is a check, agrees very well in default parameterization
        if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') this%GetOptDepth()

        !Get relevant times
        astart=1.d0/(1.d0+Reion%redshift + Reion%delta_redshift*8)
        this%tau_start = max(0.05_dl, rombint(dtauda,0._dl,astart,1d-3))
        !Time when a very small reionization fraction (assuming tanh fitting)

        this%tau_complete = min(tau0, &
            this%tau_start+ rombint(dtauda,astart,1.d0/(1.d0+max(0.d0,Reion%redshift-Reion%delta_redshift*8)),1d-3))

    end if

    end subroutine Reionization_Init

    subroutine Reionization_Validate(Reion, OK)
    class(ReionizationParams),intent(in) :: Reion
    logical, intent(inout) :: OK

    if (Reion%Reionization) then
        if (Reion%use_optical_depth) then
            if (Reion%optical_depth<0 .or. Reion%optical_depth > 0.9  .or. &
                include_helium_fullreion .and. Reion%optical_depth<0.01) then
                OK = .false.
                write(*,*) 'Optical depth is strange. You have:', Reion%optical_depth
            end if
        else
            if (Reion%redshift < 0 .or. Reion%Redshift +Reion%delta_redshift*3 > Reionization_maxz .or. &
                include_helium_fullreion .and. Reion%redshift < Reion%helium_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',Reion%Redshift
            end if
        end if
        if (Reion%fraction/= Reionization_DefFraction .and. (Reion%fraction < 0 .or. Reion%fraction > 1.5)) then
            OK = .false.
            write(*,*) 'Reionization fraction strange. You have: ',Reion%fraction
        end if
        if (Reion%delta_redshift > 3 .or. Reion%delta_redshift<0.1 ) then
            !Very narrow windows likely to cause problems in interpolation etc.
            !Very broad likely to conflict with quasar data at z=6
            OK = .false.
            write(*,*) 'Reionization delta_redshift is strange. You have: ',Reion%delta_redshift
        end if
    end if

    end subroutine Reionization_Validate


    function Reionization_doptdepth_dz(this,z)
    class(ReionizationHistory) :: this
    real(dl) :: Reionization_doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) a

    a = 1._dl/(1._dl+z)

    Reionization_doptdepth_dz = this%x_e(a)*this%akthom*dtauda(a)

    end function Reionization_doptdepth_dz

    function Reionization_GetOptDepth(this)
    use MathUtils
    class(ReionizationHistory) :: this
    real(dl) Reionization_GetOptDepth

    Reionization_GetOptDepth = Integrate_Romberg_classfunc(this, Reionization_doptdepth_dz,0.d0,Reionization_maxz,&
        Reionization_tol, 20, nint(Reionization_maxz/this%Params%delta_redshift*5))

    end function Reionization_GetOptDepth

    subroutine Reionization_zreFromOptDepth(this, Reion)
    !General routine to find zre parameter given optical depth
    use Errors
    class(ReionizationHistory) :: this
    class(ReionizationParams), intent(in) :: Reion
    real(dl) try_b, try_t
    real(dl) tau
    integer i

    try_b = 0
    try_t = Reionization_maxz
    i=0
    do
        i=i+1
        this%Params%redshift = (try_t + try_b)/2
        call this%SetParamsForZre(this%Params)
        tau = this%GetOptDepth()

        if (tau > Reion%optical_depth) then
            try_t = this%Params%redshift
        else
            try_b = this%Params%redshift
        end if
        if (abs(try_b - try_t) < 2e-3/Reionization_AccuracyBoost) exit
        if (i>100) call GlobalError('Reionization_zreFromOptDepth: failed to converge',error_reionization)
    end do

    if (abs(tau - Reion%optical_depth) > 0.002 .and. global_error_flag==0) then
        write (*,*) 'Reionization_zreFromOptDepth: Did not converge to optical depth'
        write (*,*) 'tau =',tau, 'optical_depth = ', Reion%optical_depth
        write (*,*) try_t, try_b
        write (*,*) '(If running a chain, have you put a constraint on tau?)'
        call GlobalError('Reionization did not converge to optical depth',error_reionization)
    end if

    end subroutine Reionization_zreFromOptDepth


    subroutine Reionization_SetFromOptDepth(this, Reion)
    class(ReionizationHistory) :: this
    Type(ReionizationParams) :: Reion
    ! This subroutine calculates the redshift of reionizatio
    ! This implementation is approximate but quite accurate and fast

    Reion%redshift = 0

    if (Reion%Reionization .and. Reion%optical_depth /= 0) then
        !Do binary search to find zre from z
        !This is general method
        call this%zreFromOptDepth(Reion)
    else
        Reion%Reionization = .false.
    end if

    end subroutine Reionization_SetFromOptDepth

    end module Reionization

