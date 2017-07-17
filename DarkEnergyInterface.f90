    module DarkEnergyInterface
    use precision
    implicit none

    type, abstract :: TDarkEnergyBase
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa_ppf = 0._dl !Not used here, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed
        integer :: num_perturb_equations = 0
        logical :: is_cosmological_constant = .true.
    contains
    procedure :: ReadParams
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: InitializeYfromVec
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


    subroutine BackgroundDensityAndPressure(this, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyBase), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    grhov_t  =0
    if (present(w)) w=-1

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


    subroutine InitializeYfromVec(this, y, EV_w_ix, InitVecAti_clxq, InitVecAti_vq)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    integer, intent(in) :: EV_w_ix
    real(dl), intent(in) :: InitVecAti_clxq, InitVecAti_vq
    end subroutine InitializeYfromVec


    subroutine PerturbationEvolve(this, ayprime, w_ix, &
        adotoa, k, z, y)
    class(TDarkEnergyBase), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: adotoa, k, z, y(:)
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    end module DarkEnergyInterface
