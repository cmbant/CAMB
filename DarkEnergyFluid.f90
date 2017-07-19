
    module DarkEnergyFluid
    use DarkEnergyInterface
    implicit none

    type, extends(TDarkEnergyBase) :: TDarkEnergyFluid
        !comoving sound speed is always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure :: Init =>TDarkEnergyFluid_Init
    procedure :: BackgroundDensityAndPressure => TDarkEnergyFluid_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TDarkEnergyFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    final :: TDarkEnergyFluid_Finalize
    end type TDarkEnergyFluid

    contains


    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    this%w_lam = Ini%Read_Double('w', -1.d0)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    if (Ini%Read_Double('wa', 0.d0) /= 0) error stop 'wa not supported by fluid model'
    call this%Init()

    end subroutine TDarkEnergyFluid_ReadParams


    subroutine TDarkEnergyFluid_Init(this)
    use ModelParams
    class(TDarkEnergyFluid), intent(inout) :: this

    this%is_cosmological_constant = abs(this%w_lam + 1._dl) < 1.e-6_dl

    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyFluid_Init


    subroutine TDarkEnergyFluid_BackgroundDensityAndPressure(this, a, grhov_t, w)
    use ModelParams
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w


    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
    else
        grhov_t = grhov * a ** (-1 - 3 * this%w_lam)
    end if
    if (present(w)) w = this%w_lam

    end subroutine TDarkEnergyFluid_BackgroundDensityAndPressure

    subroutine TDarkEnergyFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t * (1 + this%w_lam)

    end subroutine TDarkEnergyFluid_PerturbedStressEnergy


    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w_ix, &
        adotoa, k, z, y)
    use ModelParams
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: adotoa, k, z, y(:)
    integer, intent(in) :: w_ix

    ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - this%w_lam) * &
        (y(w_ix) + 3 * adotoa * (1 + this%w_lam) * y(w_ix + 1) / k) - &
        (1 + this%w_lam) * k * y(w_ix + 1) - (1 + this%w_lam) * k * z

    ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
        k * this%cs2_lam * y(w_ix) / (1 + this%w_lam)

    end subroutine TDarkEnergyFluid_PerturbationEvolve



    subroutine TDarkEnergyFluid_Finalize(this)
    type(TDarkEnergyFluid), intent(inout) :: this

    end subroutine TDarkEnergyFluid_Finalize

    end module DarkEnergyFluid
