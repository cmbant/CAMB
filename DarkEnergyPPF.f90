    module DarkEnergyPPF
    use precision
    use ModelParams
    use RandUtils
    use DarkEnergyInterface
    implicit none

    type, extends(TDarkEnergy) :: TDarkEnergyPPF
        ! w_lam is now w0
        !cs2_lam now is ce^2

        logical :: use_tabulated_w = .false.
        real(dl) :: c_Gamma_ppf = 0.4_dl
        integer, parameter :: nwmax = 5000
        integer, private, parameter :: nde = 2000
        integer :: nw_ppf
        real(dl) w_ppf(nwmax), a_ppf(nwmax)
        real(dl), private :: ddw_ppf(nwmax)
        real(dl), private :: rde(nde),ade(nde),ddrde(nde)
        real(dl), private, parameter :: amin = 1.d-9
        logical :: is_cosmological_constant
    contains
    procedure ReadParams => TDarkEnergyPPF_ReadParams
    procedure Init_Background => TDarkEnergyPPF_Init_Background
    procedure dtaudaFctn => TDarkEnergyPPF_dtauda
    procedure, private :: setddwa
    procedure, private :: interpolrde
    procedure, private :: grho_de
    procedure, private :: setcgammappf
    end type TDarkEnergyPPF

    private w_de, cubicsplint
    contains

    type(TDarkEnergyPPF), private, pointer :: curr

    subroutine TDarkEnergyPPF_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyPPF), intent(inout) :: this
    Type(TIniFile) :: Ini
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
        do i=1, this%nwmax+1
            read(10, *, end=100) this%a_ppf(i), this%w_ppf(i)
            this%a_ppf(i) = dlog(this%a_ppf(i))
            this%nw_ppf = this%nw_ppf + 1
        enddo
        write(*,'("Note: ", a, " has more than ", I8, " data points")') &
        &   trim(wafile), this%nwmax
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

    end subroutine DarkEnergy_ReadParams


    subroutine setddwa(this)
    class(TDarkEnergyPPF) :: this
    real(dl), parameter :: wlo=1.d30, whi=1.d30

    call spline(this%a_ppf, this%w_ppf, this%nw_ppf, wlo, whi, this%ddw_ppf) !a_ppf is lna here

    end subroutine setddwa


    function w_de(a)
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    if(.not. use_tabulated_w) then
        w_de= curr%w_lam+ curr%wa_ppf*(1._dl-a)
    else
        al=dlog(a)
        if(al .lt. curr%a_ppf(1)) then
            w_de= curr%w_ppf(1)                   !if a < minimum a from wa.dat
        elseif(al .gt. curr%a_ppf(curr%nw_ppf)) then
            w_de= curr%w_ppf(curr%nw_ppf)         !if a > maximus a from wa.dat
        else
            call cubicsplint(cur%a_ppf, curr%w_ppf, curr%ddw_ppf, &
            &   curr%nw_ppf, al, w_de)
        endif
    endif
    end function w_de  ! equation of state of the PPF DE


    function drdlna_de(al)
    real(dl) :: drdlna_de, a
    real(dl), intent(IN) :: al

    a=dexp(al)
    drdlna_de=3._dl*(1._dl+w_de(a))

    end function drdlna_de


    subroutine interpolrde(this)
    class(TDarkEnergyPPF), target :: this
    real(dl), parameter :: rlo=1.d30, rhi=1.d30
    real(dl) :: atol, almin, al, rombint, fint
    integer :: i
    external rombint
    curr => this
    atol= 1.d-5
    almin= dlog(this%amin)
    do i=1, this%nde
        al= almin-almin/(this%nde-1)*(i-1)    !interpolate between amin and today
        fint=rombint(drdlna_de, al, 0._dl, atol)+4._dl*al
        ade(i)=al
        rde(i)=dexp(fint) !rho_de*a^4 normalize to its value at today
    enddo
    call spline(ade,rde,nde,rlo,rhi,ddrde)
    nullify(curr)
    end subroutine interpolrde


    function grho_de(this, a)  !8 pi G a^4 rho_de
    class(TDarkEnergyPPF), target :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        grho_de= grhov * a**(1._dl-3.* this%w_lam-3.* this%wa_ppf)* &
        &   exp(-3.* this%wa_ppf*(1._dl-a))
    else
        if(a.eq.0.d0)then
            grho_de=0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            al=dlog(a)
            curr => this
            if(al.lt.this%ade(1))then
                fint=this%rde(1)*(a/this%amin)**(1.-3.*w_de(this%amin))    !if a<amin, assume here w=w_de(amin)
            else              !if amin is small enough, this extrapolation will be unnecessary.
                call cubicsplint(this%ade,this%rde,this%ddrde,this%nde,al,fint)
            endif
            grho_de=grhov*fint
            nullify(curr)
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

    this%c_Gamma_ppf=0.4d0*sqrt(this%cs2_lam)

    end subroutine setcgammappf


    subroutine TDarkEnergyPPF_Init_Background(this)
    class(TDarkEnergyPPF) :: this
    !This is only called once per model, and is a good point to do any extra initialization.
    !It is called before first call to dtauda, but after
    !massive neutrinos are initialized and after GetOmegak
    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
    &   this%w_lam==-1_dl .and. this%wa_ppf==0._dl
    end  subroutine TDarkEnergyPPF_Init_Background


    !Background evolution
    function TDarkEnergyPPF_dtauda(this, a) result(dtauda)
    !get d tau / d a
    use precision
    use ModelParams
    use MassiveNu
    implicit none
    class (TDarkEnergyPPF), intent(in) :: this
    real(dl) dtauda
    real(dl), intent(IN) :: a
    real(dl) rhonu,grhoa2, a2
    integer nu_i

    a2=a**2

    !  8*pi*G*rho*a**4.
    grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
    if (this%is_cosmological_constant) then
        grhoa2=grhoa2+grhov*a2**2
    else
        grhoa2=grhoa2+ grho_de(a)
    end if

    if (CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grhoa2=grhoa2+rhonu*grhormass(nu_i)
        end do
    end if

    dtauda=sqrt(3/grhoa2)

    end function TDarkEnergyPPF_dtauda

    end module DarkEnergyPPF
