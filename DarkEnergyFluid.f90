
    module DarkEnergyFluidModule
    use DarkEnergyInterface
    implicit none

    type, extends(TDarkEnergyBase) :: TDarkEnergyFluid
        real(dl) :: cs2_lam = 1_dl
        !comoving sound speed. Always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure :: Init =>TDarkEnergyFluid_Init
    procedure :: Init_Background => TDarkEnergyFluid_Init_Background
    procedure :: BackgroundDensity => TDarkEnergyFluid_BackgroundDensity
    procedure :: AddStressEnergy => TDarkEnergyFluid_AddStressEnergy
    procedure :: PrepDerivs => TDarkEnergyFluid_PrepDerivs
    procedure :: BackgroundDensityAndPressure => TDarkEnergyFluid_BackgroundDensityAndPressure
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    procedure :: DerivsAdd2Gpres => TDarkEnergyFluid_DerivsAdd2Gpres
    final :: TDarkEnergyFluid_Finalize
    end type TDarkEnergyFluid

    contains


    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    this%w_lam = Ini%Read_Double('w', -1.d0)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    call this%Init()
    end subroutine TDarkEnergyFluid_ReadParams


    subroutine TDarkEnergyFluid_Init(this)
    use ModelParams
    class(TDarkEnergyFluid), intent(inout) :: this

    ! is_cosmological_constant = this%w_lam == -1.dl
    ! to prevent issues with fp comparision:
    this%is_cosmological_constant = abs(this%w_lam + 1._dl) < 1.e-6_dl

    ! Set in both cases to be on the safe side.
    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyFluid_Init

    subroutine TDarkEnergyFluid_Init_Background(this)
    use GaugeInterface
    class(TDarkEnergyFluid) :: this

    ! Set the name to export for equations.
    Eqns_name = 'gauge_inv'

    end  subroutine TDarkEnergyFluid_Init_Background


    !Background evolution
    function TDarkEnergyFluid_BackgroundDensity(this, a) result(grhoa2)
    !get additional terms for d tau / d a
    use precision
    use ModelParams
    implicit none
    class (TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(IN) :: a
    real(dl) :: grhoa2

    if (this%is_cosmological_constant) then
        grhoa2 = grhov * a ** 4
    else
        grhoa2 = grhov * a ** (1 - 3 * this%w_lam)
    end if

    end function TDarkEnergyFluid_BackgroundDensity


    function TDarkEnergyFluid_AddStressEnergy(this, gpres, dgq, dgrho, &
        a, grhov) result(grhov_t)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl) :: grhov_t
    real(dl), intent(inout) :: gpres, dgq, dgrho
    real(dl), intent(in) :: a, grhov

    grhov_t = grhov * a ** (-1 - 3 * this%w_lam)
    gpres = gpres + grhov_t * this%w_lam

    end function TDarkEnergyFluid_AddStressEnergy


    subroutine TDarkEnergyFluid_PrepDerivs(this, dgrho, dgq, &
        grhov_t, y, w_ix)
    use ModelParams
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(inout) :: dgrho, dgq
    real(dl), intent(in) :: grhov_t, y(:)
    integer, intent(in) :: w_ix

    if (.not. this%is_cosmological_constant) then
        dgrho = dgrho + y(w_ix) * grhov_t
        dgq = dgq + y(w_ix + 1) * grhov_t * (1 + this%w_lam)
    end if

    end subroutine TDarkEnergyFluid_PrepDerivs


    subroutine TDarkEnergyFluid_BackgroundDensityAndPressure(this, grhov_t, &
        a, grhov, grhor_t, grhog_t, gpres)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(inout) :: grhov_t
    real(dl), intent(inout), optional :: gpres
    real(dl), intent(in) :: a, grhov, grhor_t, grhog_t

! From BackgroundDensity
!    if (this%is_cosmological_constant) then
!        grhoa2 = grhov * a ** 4
!    else
!        grhoa2 = grhov * a ** (1 - 3 * this%w_lam)
!    end if

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
    else
        grhov_t = grhov * a ** (-1 - 3 * this%w_lam)
    end if
    if (present(gpres)) gpres = 0

    end subroutine TDarkEnergyFluid_BackgroundDensityAndPressure


    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w_ix, &
        adotoa, k, z, y)
    use ModelParams
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: adotoa, k, z, y(:)
    integer, intent(in) :: w_ix
    integer :: threadID

     if (.not. this%is_cosmological_constant) then
        threadID = GetThreadID()
        ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - this%w_lam) * &
			(y(w_ix) + 3 * adotoa * (1 + this%w_lam) * y(w_ix + 1) / k) - &
            (1 + this%w_lam) * k * y(w_ix + 1) - (1 + this%w_lam) * k * z

        ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
			k * this%cs2_lam * y(w_ix) / (1 + this%w_lam)
    end if

    end subroutine TDarkEnergyFluid_PerturbationEvolve


    function TDarkEnergyFluid_DerivsAdd2Gpres(this, grhog_t, grhor_t, grhov_t) result(gpres)
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(in) :: grhog_t, grhor_t, grhov_t
    real(dl) :: gpres

    gpres = (grhog_t + grhor_t) / 3 + grhov_t * this%w_lam

    end function TDarkEnergyFluid_DerivsAdd2Gpres


    subroutine TDarkEnergyFluid_Finalize(this)
    type(TDarkEnergyFluid), intent(inout) :: this

    end subroutine TDarkEnergyFluid_Finalize

    end module DarkEnergyFluidModule
