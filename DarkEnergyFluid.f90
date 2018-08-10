
    module DarkEnergyFluid
    use DarkEnergyInterface
    use constants
    implicit none

    type, extends(TDarkEnergyBase) :: TDarkEnergyFluid
        !comoving sound speed is always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure :: Init =>TDarkEnergyFluid_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    end type TDarkEnergyFluid

    contains


    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyBase%ReadParams(Ini)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergyFluid_ReadParams


    subroutine TDarkEnergyFluid_Init(this)
    use ModelParams
    class(TDarkEnergyFluid), intent(inout) :: this

    call this%TDarkEnergyBase%Init()

    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        if (this%use_tabulated_w) then
            if (any(this%equation_of_state%F<-1)) &
                error stop 'Fluid dark energy model does not allow w crossing -1'
        elseif (this%wa/=0 .and. &
            ((1+this%w_lam < -1.e-6_dl) .or. 1+this%w_lam + this%wa < -1.e-6_dl)) then
            error stop 'Fluid dark energy model does not allow w crossing -1'
        end if
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyFluid_Init


    subroutine TDarkEnergyFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t * (1 + w)

    end subroutine TDarkEnergyFluid_PerturbedStressEnergy


    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    use ModelParams
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, loga

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - w) *  (y(w_ix) + (1 + w) * Hv3_over_k) &
        -  (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
    if (this%use_tabulated_w) then
        !account for derivatives of w
        loga = log(a)
        if (loga > this%equation_of_state%Xmin_interp .and. loga < this%equation_of_state%Xmax_interp) then
            ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)* Hv3_over_k
        end if
    elseif (this%wa/=0) then
        ayprime(w_ix) = ayprime(w_ix) + Hv3_over_k*this%wa*adotoa*a
    end if
    !velocity
    if (abs(w+1) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
            k * this%cs2_lam * y(w_ix) / (1 + w)
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TDarkEnergyFluid_PerturbationEvolve


    end module DarkEnergyFluid
