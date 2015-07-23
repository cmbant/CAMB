    module DarkEnergyPPFModule
    use precision
    use ModelParams
    use RandUtils
    use DarkEnergyInterface
    implicit none

    integer, parameter :: nwmax = 5000
    integer, private, parameter :: nde = 2000
    real(dl), private, parameter :: amin = 1.d-9

    type, extends(TDarkEnergyBase) :: TDarkEnergyPPF
        ! w_lam is now w0
        real(dl) :: cs2_lam = 1_dl
        !comoving sound speed. Always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)

        logical :: use_tabulated_w = .false.
        real(dl) :: c_Gamma_ppf = 0.4_dl
        integer :: nw_ppf
        real(dl) :: w_ppf(nwmax), a_ppf(nwmax)
        real(dl), private :: ddw_ppf(nwmax)
        real(dl), private :: rde(nde),ade(nde),ddrde(nde)

        !PPF parameters for dervis
        real(dl), private, dimension(:), allocatable :: dgrho_e_ppf, dgq_e_ppf
    contains
    procedure :: ReadParams => TDarkEnergyPPF_ReadParams
    procedure :: Init => TDarkEnergyPPF_Init
    procedure :: Init_Background => TDarkEnergyPPF_Init_Background
    procedure :: BackgroundDensityAndPressure => TDarkEnergyPPF_BackgroundDensityAndPressure
    procedure :: AddStressEnergy => TDarkEnergyPPF_AddStressEnergy
    procedure :: diff_rhopi_Add_Term => TDarkEnergyPPF_diff_rhopi_Add_Term
    procedure :: DerivsAddPreSigma => TDarkEnergyPPF_DerivsAddPreSigma
    procedure, private :: setddwa
    procedure, private :: interpolrde
    procedure, private :: grho_de
    procedure, private :: setcgammappf
    final :: TDarkEnergyPPF_Finalize
    end type TDarkEnergyPPF

    private w_de, cubicsplint
    contains

    subroutine TDarkEnergyPPF_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyPPF), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini
    character(len=:), allocatable :: wafile
    integer i

    if (Ini%HasKey('usew0wa')) then
        stop 'input variables changed from usew0wa: now use_tabulated_w or w, wa'
    end if

    this%use_tabulated_w = Ini%Read_Logical('use_tabulated_w', .false.)
    if(.not. this%use_tabulated_w)then
        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa_ppf = Ini%Read_Double('wa', 0.d0)
        if (Rand_Feedback >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa_ppf
    else
        wafile = Ini%Read_String('wafile')
        open(unit=10, file=wafile, status='old')
        this%nw_ppf=0
        do i=1, nwmax + 1
            read(10, *, end=100) this%a_ppf(i), this%w_ppf(i)
            this%a_ppf(i) = dlog(this%a_ppf(i))
            this%nw_ppf = this%nw_ppf + 1
        enddo
        write(*,'("Note: ", a, " has more than ", I8, " data points")') &
        &   trim(wafile), nwmax
        write(*,*)'Increase nwmax in LambdaGeneral'
        stop
100     close(10)
        write(*,'("read in ", I8, " (a, w) data points from ", a)') &
        &   this%nw_ppf, trim(wafile)
        call this%setddwa
        call this%interpolrde
    endif
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    call this%setcgammappf

    call this%Init()

    end subroutine TDarkEnergyPPF_ReadParams


    subroutine TDarkEnergyPPF_Init(this)
    class(TDarkEnergyPPF), intent(inout) :: this
    integer :: numThreads

    numThreads = ThreadNum
    call GetNumThreads(numThreads)
    numThreads = numThreads - 1
    if (.not. allocated(this%dgq_e_ppf)) &
        allocate(this%dgq_e_ppf(0:numThreads), this%dgrho_e_ppf(0:numThreads))

    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
    &   this%w_lam==-1_dl .and. this%wa_ppf==0._dl

    ! Set both cases to be always on the safe side.
    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 1
    end if

    end subroutine TDarkEnergyPPF_Init


    subroutine TDarkEnergyPPF_Finalize(this)
    type(TDarkEnergyPPF), intent(inout) :: this

    if (allocated(this%dgq_e_ppf)) deallocate(&
        this%dgq_e_ppf, this%dgrho_e_ppf)

    end subroutine


    subroutine setddwa(this)
    class(TDarkEnergyPPF) :: this
    real(dl), parameter :: wlo = 1.d30, whi = 1.d30

    call spline(this%a_ppf, this%w_ppf, this%nw_ppf, wlo, whi, this%ddw_ppf) !a_ppf is lna here

    end subroutine setddwa


    function w_de(curr, a)
    type(TDarkEnergyPPF), intent(in) :: curr
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    if(.not. curr%use_tabulated_w) then
        w_de= curr%w_lam+ curr%wa_ppf*(1._dl-a)
    else
        al=dlog(a)
        if(al .lt. curr%a_ppf(1)) then
            w_de= curr%w_ppf(1)                   !if a < minimum a from wa.dat
        elseif(al .gt. curr%a_ppf(curr%nw_ppf)) then
            w_de= curr%w_ppf(curr%nw_ppf)         !if a > maximus a from wa.dat
        else
            call cubicsplint(curr%a_ppf, curr%w_ppf, curr%ddw_ppf, &
            &   curr%nw_ppf, al, w_de)
        endif
    endif
    end function w_de  ! equation of state of the PPF DE


    function drdlna_de(curr, al)
    type(TDarkEnergyPPF), intent(in) :: curr
    real(dl) :: drdlna_de, a
    real(dl), intent(IN) :: al

    a = dexp(al)
    drdlna_de = 3._dl * (1._dl + w_de(curr, a))

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
                fint = this%rde(1) * (a / amin) ** (1. - 3. * w_de(this, amin))
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
    if (h.eq.0.) stop 'bad xa input in splint'
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


    subroutine TDarkEnergyPPF_Init_Background(this)
    use GaugeInterface
    class(TDarkEnergyPPF) :: this

    ! Set the name to export for equations.
    Eqns_name = 'equations_ppf-Jan15'

    end subroutine TDarkEnergyPPF_Init_Background

    subroutine TDarkEnergyPPF_BackgroundDensityAndPressure(this, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t
    real(dl), intent(out), optional :: w

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    else
        ! Ensure a valid result
        if (a /= 0._dl) then
            grhov_t = grho_de(this, a) / (a * a)
        else
            grhov_t = 0._dl
        end if
        if (present(w)) w = w_de(this, a)
    end if

    end subroutine TDarkEnergyPPF_BackgroundDensityAndPressure



    subroutine TDarkEnergyPPF_AddStressEnergy(this, dgrho, dgq, &
        grhov_t, y, w_ix, output)
    class(TDarkEnergyPPF), intent(inout) :: this
        real(dl), intent(inout) :: dgrho, dgq
    real(dl), intent(in) :: grhov_t, y(:)
    integer, intent(in) :: w_ix
    logical, intent(in) :: output
    integer :: tID

    if (output) then
        !If not in output, this is calculated in DerivsAddPreSigma
        tID = GetThreadID()
        dgrho = dgrho + this%dgrho_e_ppf(tID)
        dgq = dgq + this%dgq_e_ppf(tID)
    end if

    end subroutine TDarkEnergyPPF_AddStressEnergy


    function TDarkEnergyPPF_diff_rhopi_Add_Term(this, grho, gpres, w,  grhok, adotoa, &
        EV_KfAtOne, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    !Get derivative of anisotropic stress
    class(TDarkEnergyPPF), intent(in) :: this
    real(dl), intent(in) :: grho, gpres, w, grhok, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), EV_KfAtOne
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh
    integer :: tID

    if (this%is_cosmological_constant) then
        ppiedot = 0
    else
        tID = GetThreadID()
        hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
        ppiedot = 3._dl * this%dgrho_e_ppf(tID) + this%dgq_e_ppf(tID) * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * EV_KfAtOne * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
        ppiedot = ppiedot * adotoa / EV_KfAtOne
    end if

    end function TDarkEnergyPPF_diff_rhopi_Add_Term


    subroutine TDarkEnergyPPF_DerivsAddPreSigma(this, sigma, &
        ayprime, dgq, dgrho, &
        grho, grhov_t, w, gpres_noDE, ay, w_ix, etak, adotoa, k, k2, EV_kf1)
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(inout) :: sigma, ayprime(:), dgq, dgrho
    real(dl), intent(in) :: grho, grhov_t, w, gpres_noDE, ay(:), etak
    real(dl), intent(in) :: adotoa, k, k2, EV_kf1
    integer, intent(in) :: w_ix
    real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, dgqe, dgrhoe
    real(dl) :: vT, grhoT
    integer :: tID

    if (.not. this%is_cosmological_constant) then
        tID = GetThreadID()
        !ppf
        grhoT = grho - grhov_t
        vT = dgq / (grhoT + gpres_noDE)
        Gamma = ay(w_ix)

        !sigma for ppf
        sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / EV_kf1 - &
            k * Gamma
        sigma = sigma / adotoa

        S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
        ckH = this%c_Gamma_ppf * k / adotoa
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
        ayprime(w_ix) = Gammadot

        if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 ?????????
            Gamma = 0
            Gammadot = 0.d0
            ayprime(w_ix) = Gammadot
        endif

        Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / EV_kf1
        dgqe = S_Gamma - Gammadot / adotoa - Gamma
        dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
        dgrhoe = -2 * k2 * EV_kf1 * Gamma - 3 / k * adotoa * dgqe
        dgrho = dgrho + dgrhoe
        dgq = dgq + dgqe

        this%dgrho_e_ppf(tID) = dgrhoe
        this%dgq_e_ppf(tID) = dgqe
    end if

    end subroutine TDarkEnergyPPF_DerivsAddPreSigma


    end module DarkEnergyPPFModule
