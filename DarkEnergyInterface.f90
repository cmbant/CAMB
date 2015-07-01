    module DarkEnergyInterface
    use precision
    implicit none

    type, abstract :: TDarkEnergyBase
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)

        real(dl) :: wa_ppf = 0._dl !Not used here, just for compatibility with e.g. halofit
        integer :: num_perturb_equations = 0
        logical :: is_cosmological_constant
    contains
    procedure :: ReadParams
    procedure :: Init
    procedure :: Init_Background
    procedure :: BackgroundDensityAndPressure
    procedure :: AddStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: InitializeYfromVec
    procedure :: DerivsAddPreSigma
    procedure :: PerturbationEvolve
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

    subroutine BackgroundDensityAndPressure(this, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t
    real(dl), intent(out), optional :: w
    
    grhov_t  =0
    
    end subroutine BackgroundDensityAndPressure


    subroutine AddStressEnergy(this, dgrho, dgq, grhov_t, y, w_ix, output) 
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: dgrho, dgq
    real(dl), intent(in) :: grhov_t, y(:)
    integer, intent(in) :: w_ix
    logical, intent(in) :: output
    end subroutine AddStressEnergy


    function diff_rhopi_Add_Term(this, grho, gpres, w, grhok, adotoa, &
        EV_KfAtOne, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(in) :: grho, gpres, grhok, w, adotoa, &
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


    subroutine DerivsAddPreSigma(this, sigma, ayprime, dgq, dgrho, &
        grho, grhov_t, w, gpres_noDE, ay, w_ix, etak, adotoa, k, k2, EV_kf1)
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(inout) :: sigma, ayprime(:), dgq, dgrho
    real(dl), intent(in) :: grho, grhov_t, w, gpres_noDE, ay(:), etak
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

    end module DarkEnergyInterface
