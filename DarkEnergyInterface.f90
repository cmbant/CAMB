    module DarkEnergyInterface
    use precision
    use interpolation

    implicit none

    private

    type TDarkEnergyBase
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa = 0._dl !may not be used, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed, though may not be used
        logical :: use_tabulated_w = .false.  !Use interpolated table; note this is quite slow.
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
        logical :: is_cosmological_constant = .true.
        integer :: num_perturb_equations = 0
        !Interpolations if use_tabulated_w=.true.
        Type(TCubicSpline) :: equation_of_state, logdensity
    contains
    procedure :: ReadParams
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: PerturbationEvolve
    procedure :: SetwTable => TDarkEnergyBase_SetwTable
    procedure :: PrintFeedback
    procedure :: w_de
    procedure :: grho_de
    end type TDarkEnergyBase

    public TDarkEnergyBase
    contains

    function w_de(this, a)
    class(TDarkEnergyBase) :: this
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        w_de= this%w_lam+ this%wa*(1._dl-a)
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            w_de= this%equation_of_state%F(1)
        elseif(al >= this%equation_of_state%Xmax_interp) then
            w_de= this%equation_of_state%F(this%equation_of_state%n)
        else
            w_de = this%equation_of_state%Value(al)
        endif
    endif

    end function w_de  ! equation of state of the PPF DE

    function grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyBase) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        grho_de = a ** (1._dl - 3. * this%w_lam - 3. * this%wa)
        if (this%wa/=0) grho_de=grho_de*exp(-3. * this%wa * (1._dl - a))
    else
        if(a == 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            if (a>=1) then
                fint= 1
            else
                al = dlog(a)
                if(al <= this%logdensity%X(1)) then
                    ! assume here w=w_de(a_min)
                    fint = exp(this%logdensity%F(1) + (1. - 3. * this%equation_of_state%F(1))*(al - this%logdensity%X(1)))
                else
                    fint = exp(this%logdensity%Value(al))
                endif
            end if
            grho_de = fint
        endif
    endif

    end function grho_de

    subroutine PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyBase) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa

    end subroutine PrintFeedback

    subroutine ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyBase), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)

    this%use_tabulated_w = Ini%Read_Logical('use_tabulated_w', .false.)
    if(.not. this%use_tabulated_w)then
        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa = Ini%Read_Double('wa', 0.d0)
    else
        call File%LoadTxt(Ini%Read_String('wafile'), table)
        call this%SetwTable(table(:,1),table(:,2))
    endif

    end subroutine ReadParams


    subroutine Init(this, omegav)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(in) :: omegav
    
    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
        &  abs(this%w_lam + 1._dl) < 1.e-6_dl .and. this%wa==0._dl

    end subroutine Init


    subroutine TDarkEnergyBase_SetwTable(this, a, w)
    class(TDarkEnergyBase) :: this
    real(dl), intent(in) :: a(:), w(:)
    real(dl), allocatable :: integral(:)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'w table must end at a=1'

    this%use_tabulated_w = .true.
    call this%equation_of_state%Init(log(a), w)

    allocate(integral(this%equation_of_state%n))
    ! log (rho) =  -3 int dlna (1+w)
    call spline_integrate(this%equation_of_state%X,this%equation_of_state%F, &
        this%equation_of_state%ddF,integral,this%equation_of_state%n)
    integral  = -3*( (this%equation_of_state%X-this%equation_of_state%X(1)) + integral) + 4*this%equation_of_state%X
    integral = integral - integral(this%equation_of_state%n) !log(a^4 rho_de)) normalized to 0 at a=1
    call this%logdensity%Init(this%equation_of_state%X, integral)
    !Set w and wa to values today (e.g. as too simple first guess for approx fittings etc).
    this%w_lam = w(size(a))
    this%wa = -this%equation_of_state%Derivative(0._dl)

    end subroutine TDarkEnergyBase_SetwTable

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    else
        ! Ensure a valid result
        if (a > 1e-10) then
            grhov_t = grhov * this%grho_de(a) / (a * a)
        else
            grhov_t = 0._dl
        end if
        if (present(w)) w = this%w_de(a)
    end if

    end subroutine BackgroundDensityAndPressure


    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe=0
    dgqe=0

    end subroutine PerturbedStressEnergy


    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, y(:), w
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    end module DarkEnergyInterface
