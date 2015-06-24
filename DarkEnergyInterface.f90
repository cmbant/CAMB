    module DarkEnergyInterface
    use precision
    implicit none

    type :: TDarkEnergy
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (assumed constant)
        real(dl) :: cs2_lam = 1_dl
        !comoving sound speed. Always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)

        real(dl) :: wa_ppf = 0._dl !Not used here, just for compatibility with e.g. halofit

        logical :: w_perturb = .true.
        !If you are tempted to set this = .false. read
        ! http://cosmocoffee.info/viewtopic.php?t=811
        ! http://cosmocoffee.info/viewtopic.php?t=512
    contains
    procedure :: ReadParams => TDarkEnergy_ReadParams
    procedure :: Init_Background => TDarkEnergy_Init_Background
    procedure :: dtaudaFctn => TDarkEnergy_dtauda
    end type TDarkEnergy

    class(TDarkEnergy), allocatable :: DarkEnergy

    contains

    subroutine TDarkEnergy_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergy), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    this%w_lam = Ini%Read_Double('w', -1.d0)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergy_ReadParams

    subroutine TDarkEnergy_Init_Background(this)
    class(TDarkEnergy) :: this
    !This is only called once per model, and is a good point to do any extra initialization.
    !It is called before first call to dtauda, but after
    !massive neutrinos are initialized and after GetOmegak
    end  subroutine TDarkEnergy_Init_Background

    !Background evolution
    function TDarkEnergy_dtauda(this, a) result(dtauda)
    !get d tau / d a
    use precision
    use ModelParams
    use MassiveNu
    implicit none
    class (TDarkEnergy), intent(in) :: this
    real(dl) dtauda
    real(dl), intent(IN) :: a
    real(dl) rhonu,grhoa2, a2
    integer nu_i

    a2=a**2

    !  8*pi*G*rho*a**4.
    grhoa2 = grhok *a2+(grhoc+grhob)*a+grhog+grhornomass
    if (this%w_lam == -1._dl) then
        grhoa2=grhoa2+grhov*a2**2
    else
        grhoa2=grhoa2+grhov*a**(1-3*this%w_lam)
    end if
    if (CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grhoa2=grhoa2+rhonu*grhormass(nu_i)
        end do
    end if

    dtauda=sqrt(3/grhoa2)

    end function TDarkEnergy_dtauda

    end module DarkEnergyInterface


	! Wrapper to get dtauda of DarkEnergy
    function dtauda(a)
	use DarkEnergyInterface
	implicit none
	real(dl), intent(in) :: a
	real(dl) :: dtauda

	dtauda = DarkEnergy%dtaudaFctn(a)

	end function dtauda
