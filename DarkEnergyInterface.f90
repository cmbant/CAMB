    module DarkEnergyInterface
    use precision
    implicit none

    type, abstract :: TDarkEnergyBase
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (assumed constant)

        real(dl) :: wa_ppf = 0._dl !Not used here, just for compatibility with e.g. halofit
        integer :: num_perturb_equations = 0
        logical :: is_cosmological_constant
    contains
    procedure :: ReadParams
    procedure :: Init
    procedure :: Init_Background
    procedure :: BackgroundDensity
    procedure :: AddStressEnergy
    procedure :: PrepDerivs
    procedure :: diff_rhopi_Add_Term
    procedure :: InitializeYfromVec
    procedure :: BackgroundDensityAndPressure
    procedure :: DerivsAddPreSigma
    procedure :: PerturbationEvolve
    procedure :: DerivsAdd2Gpres
    end type TDarkEnergyBase

    integer, external :: GetThreadID
    contains

    subroutine ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyBase), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini
    end subroutine ReadParams


    subroutine Init(this)
    class(TDarkEnergyBase), intent(inout) :: this
    end subroutine Init


    subroutine Init_Background(this)
    class(TDarkEnergyBase) :: this
    end subroutine Init_Background


    function BackgroundDensity(this, a)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: BackgroundDensity

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    BackgroundDensity = 0._dl

    end function BackgroundDensity


    function AddStressEnergy(this, gpres, dgq, dgrho, &
        a, grhov) result(grhov_t)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: gpres, dgq, dgrho
    real(dl) :: grhov_t
    real(dl), intent(in) :: a, grhov
    end function AddStressEnergy


    subroutine PrepDerivs(this, dgrho, dgq, &
        grhov_t, y, w_ix)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: dgrho, dgq
    real(dl), intent(in) :: grhov_t, y(:)
    integer, intent(in) :: w_ix
    end subroutine PrepDerivs


    function diff_rhopi_Add_Term(this, grho, gpres, grhok, adotoa, &
        EV_KfAtOne, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(in) :: grho, gpres, grhok, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), EV_KfAtOne
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term


    subroutine InitializeYfromVec(this, y, EV_w_ix, InitVecAti_clxq, InitVecAti_vq)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    integer, intent(in) :: EV_w_ix
    real(dl), intent(in) :: InitVecAti_clxq, InitVecAti_vq
    end subroutine InitializeYfromVec


    subroutine BackgroundDensityAndPressure(this, grhov_t, &
        a, grhov, grhor_t, grhog_t, gpres)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: grhov_t
    real(dl), intent(inout), optional :: gpres
    real(dl), intent(in) :: a, grhov, grhor_t, grhog_t
    end subroutine BackgroundDensityAndPressure


    subroutine DerivsAddPreSigma(this, sigma, ayprime, dgq, dgrho, &
        grho, grhov_t,  gpres, ay, w_ix, etak, adotoa, k, k2, EV_kf1)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: sigma, ayprime(:), dgq, dgrho
    real(dl), intent(in) :: grho, grhov_t, gpres, ay(:), etak
    real(dl), intent(in) :: adotoa, k, k2, EV_kf1
    integer, intent(in) :: w_ix
    end subroutine DerivsAddPreSigma


    subroutine PerturbationEvolve(this, ayprime, w_ix, &
        adotoa, k, z, y)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: adotoa, k, z, y(:)
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve


    function DerivsAdd2Gpres(this, grhog_t, grhor_t, grhov_t) result(gpres)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(in) :: grhog_t, grhor_t, grhov_t
    real(dl) :: gpres

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    gpres = 0._dl

    end function DerivsAdd2Gpres

    end module DarkEnergyInterface
