    module DarkEnergyPPF
    use precision
    use ModelParams
    use RandUtils
    use DarkEnergyInterface
    use interpolation
    implicit none
    
    private

    integer, private, parameter :: nde = 2000
    real(dl), private, parameter :: amin = 1.d-9

    type, extends(TDarkEnergyBase) :: TDarkEnergyPPF
        logical :: use_tabulated_w = .false.
        real(dl) :: c_Gamma_ppf = 0.4_dl
        Type(TCubicSpline) :: equation_of_state
        real(dl), private :: rde(nde),ade(nde),ddrde(nde)

    contains
    procedure :: ReadParams => TDarkEnergyPPF_ReadParams
    procedure :: Init => TDarkEnergyPPF_Init
    procedure :: BackgroundDensityAndPressure => TDarkEnergyPPF_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TDarkEnergyPPF_PerturbedStressEnergy
    procedure :: diff_rhopi_Add_Term => TDarkEnergyPPF_diff_rhopi_Add_Term
    procedure :: w_de
    procedure :: grho_de
    procedure, private :: interpolrde
    procedure, private :: setcgammappf
    final :: TDarkEnergyPPF_Finalize
    end type TDarkEnergyPPF

    public TDarkEnergyPPF
    contains

    subroutine TDarkEnergyPPF_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyPPF), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini
    character(len=:), allocatable :: wafile
    real(dl), allocatable :: table(:,:)
    integer i

    if (Ini%HasKey('usew0wa')) then
        error stop 'input variables changed from usew0wa: now use_tabulated_w or w, wa'
    end if

    this%use_tabulated_w = Ini%Read_Logical('use_tabulated_w', .false.)
    if(.not. this%use_tabulated_w)then
        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa_ppf = Ini%Read_Double('wa', 0.d0)
        if (Rand_Feedback >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
            &   this%w_lam, this%wa_ppf
    else
        wafile = Ini%Read_String('loga_w_file')
        call this%equation_of_state%InitFromFile(wafile) 
        call this%interpolrde
    endif
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    if (this%cs2_lam /= 1.d0) error stop 'cs2_lam not supported by PPF model'
    call this%setcgammappf

    call this%Init()

    end subroutine TDarkEnergyPPF_ReadParams
    
    subroutine TDarkEnergyPPF_Init(this)
    class(TDarkEnergyPPF), intent(inout) :: this

    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
        &  abs(this%w_lam + 1._dl) < 1.e-6_dl .and. this%wa_ppf==0._dl

    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 1
    end if

    end subroutine TDarkEnergyPPF_Init


    subroutine TDarkEnergyPPF_Finalize(this)
    type(TDarkEnergyPPF), intent(inout) :: this

    end subroutine


    function w_de(this, a)
    type(TDarkEnergyPPF) :: this
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        w_de= this%w_lam+ this%wa_ppf*(1._dl-a)
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            w_de= this%equation_of_state%X(1)                  
        elseif(al >= this%equation_of_state%Xmax_interp) then
            w_de= this%equation_of_state%X(this%equation_of_state%n)        
        else
            w_de = this%equation_of_state%Value(al)      
        endif
    endif
    
    end function w_de  ! equation of state of the PPF DE


    function drdlna_de(this, a)
    type(TDarkEnergyPPF) :: this
    real(dl) :: drdlna_de
    real(dl), intent(IN) :: a

    drdlna_de = 3._dl * (1._dl + this%w_de(log(a)))

    end function drdlna_de


    subroutine interpolrde(this)
    class(TDarkEnergyPPF) :: this
    real(dl), parameter :: rlo=1.d30, rhi=1.d30
    real(dl) :: atol, almin, al, rombint, fint
    integer :: i
    external rombint

    atol = 1.d-5
    almin = dlog(amin)
    do i = 1, nde
        al = almin - almin / (nde - 1) * (i - 1) !interpolate between amin and today
        fint = rombint(drdlna_de, al, 0._dl, atol) + 4._dl * al
        this%ade(i) = al
        this%rde(i) = dexp(fint) !rho_de*a^4 normalize to its value at today
    enddo
    call spline(this%ade, this%rde, nde, rlo, rhi, this%ddrde)

    end subroutine interpolrde


    function grho_de(this, a)  !8 pi G a^4 rho_de
    class(TDarkEnergyPPF) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        grho_de = grhov * a ** (1._dl - 3. * this%w_lam - 3. * this%wa_ppf) * &
            exp(-3. * this%wa_ppf * (1._dl - a))
    else
        if(a .eq. 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            al = dlog(a)
            if(al .lt. this%ade(1))then
                !if a<amin, assume here w=w_de(amin)
                fint = this%rde(1) * (a / amin) ** (1. - 3. * this%w_de(amin))
            else        !if amin is small enough, this extrapolation will be unnecessary.
                call cubicsplint(this%ade, this%rde, this%ddrde, nde, al, fint)
            endif
            grho_de = grhov*fint
        endif
    endif

    end function grho_de


    !-------------------------------------------------------------------
    SUBROUTINE cubicsplint(xa,ya,y2a,n,x,y)
    INTEGER n
    real(dl)x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    real(dl)a,b,h
    klo=1
    khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
            khi=k
        else
            klo=k
        endif
        goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) error stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+&
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
    END SUBROUTINE cubicsplint
    !--------------------------------------------------------------------


    subroutine setcgammappf(this)
    class(TDarkEnergyPPF) :: this

    this%c_Gamma_ppf = 0.4_dl * sqrt(this%cs2_lam)

    end subroutine setcgammappf


    subroutine TDarkEnergyPPF_BackgroundDensityAndPressure(this, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    else
        ! Ensure a valid result
        if (a /= 0._dl) then
            grhov_t = this%grho_de(a) / (a * a)
        else
            grhov_t = 0._dl
        end if
        if (present(w)) w = this%w_de(a)
    end if

    end subroutine TDarkEnergyPPF_BackgroundDensityAndPressure



    function TDarkEnergyPPF_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    !Get derivative of anisotropic stress
    class(TDarkEnergyPPF), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh

    if (this%is_cosmological_constant) then
        ppiedot = 0
    else
        hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
        ppiedot = 3._dl * dgrhoe + dgqe * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
        ppiedot = ppiedot * adotoa / Kf1
    end if

    end function TDarkEnergyPPF_diff_rhopi_Add_Term


    subroutine TDarkEnergyPPF_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
    real(dl) :: vT, grhoT, k2
    
    k2=k**2
    !ppf
    grhoT = grho - grhov_t
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = this%c_Gamma_ppf * k / adotoa

    if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 ?????????
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe = S_Gamma - Gammadot / adotoa - Gamma
    dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
    dgrhoe = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe
    
    end subroutine TDarkEnergyPPF_PerturbedStressEnergy



    end module DarkEnergyPPF
