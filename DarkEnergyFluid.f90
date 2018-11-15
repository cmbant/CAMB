
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

    !Example implementation of fluid model using specific analytic form
    !(c_s^2=1 approximate effective axion fluid model from arXiv:1806.10608)
    type, extends(TDarkEnergyFluid) :: TAxionEffectiveFluid
        real(dl) :: w_n = 1._dl !Effective equation of state when oscillating
        real(dl) :: Om = 0._dl !Omega of the early DE component today (assumed to be negligible compared to omega_lambda)
        real(dl) :: a_c  !transition scale factor
        real(dl), private :: pow, omL, acpow !cached internally
    contains
    procedure :: ReadParams =>  TAxionEffectiveFluid_ReadParams
    procedure :: Init => TAxionEffectiveFluid_Init
    procedure :: w_de => TAxionEffectiveFluid_w_de
    procedure :: grho_de => TAxionEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TAxionEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionEffectiveFluid_PerturbationEvolve
    end type TAxionEffectiveFluid

    contains

    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyBase%ReadParams(Ini)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergyFluid_ReadParams


    subroutine TDarkEnergyFluid_Init(this, omegav)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(in) :: omegav

    call this%TDarkEnergyBase%Init(omegav)

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



    subroutine TAxionEffectiveFluid_ReadParams(this, Ini)
    use IniObjects
    class(TAxionEffectiveFluid), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyBase%ReadParams(Ini)
    call Ini%Read('cs2_lam', this%cs2_lam)
    this%w_n  = Ini%Read_Double('AxionEffectiveFluid_w_n')
    this%om  = Ini%Read_Double('AxionEffectiveFluid_om')
    this%a_c  = Ini%Read_Double('AxionEffectiveFluid_a_c')

    end subroutine TAxionEffectiveFluid_ReadParams

    subroutine TAxionEffectiveFluid_Init(this, omegav)
    class(TAxionEffectiveFluid), intent(inout) :: this
    real(dl), intent(in) :: omegav

    this%use_tabulated_w = .false.
    this%is_cosmological_constant = this%om==0
    this%pow = 3*(1+this%w_n)
    this%omL = omegav
    this%acpow = this%a_c**this%pow
    this%num_perturb_equations = 2

    end subroutine TAxionEffectiveFluid_Init


    function TAxionEffectiveFluid_w_de(this, a)
    class(TAxionEffectiveFluid) :: this
    real(dl) :: TAxionEffectiveFluid_w_de, al
    real(dl), intent(IN) :: a
    real(dl) :: rho, apow, acpow

    apow = a**this%pow
    acpow = this%acpow
    rho = this%omL+ this%om*(1+acpow)/(apow+acpow)
    TAxionEffectiveFluid_w_de = this%om*(1+acpow)/(apow+acpow)**2*(1+this%w_n)*apow/rho - 1

    end function TAxionEffectiveFluid_w_de

    function TAxionEffectiveFluid_grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TAxionEffectiveFluid) :: this
    real(dl) :: TAxionEffectiveFluid_grho_de, apow
    real(dl), intent(IN) :: a

    if(a == 0.d0)then
        TAxionEffectiveFluid_grho_de = 0.d0
    else
        apow = a**this%pow
        TAxionEffectiveFluid_grho_de = (this%omL*(apow+this%acpow)+this%om*(1+this%acpow))*a**4 &
            /((apow+this%acpow)*(this%omL+this%om))
    endif

    end function TAxionEffectiveFluid_grho_de

    subroutine TAxionEffectiveFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TAxionEffectiveFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv, apow, acpow

    apow = a**this%pow
    acpow = this%acpow
    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    ! dw/dlog a/(1+w)
    deriv  = (acpow**2*(this%om+this%omL)+this%om*acpow-apow**2*this%omL)*this%pow &
        /((apow+acpow)*(this%omL*(apow+acpow)+this%om*(1+acpow)))
    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam - deriv) * y(w_ix + 1) + &
        k * this%cs2_lam * y(w_ix)

    end subroutine TAxionEffectiveFluid_PerturbationEvolve


    subroutine TAxionEffectiveFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TAxionEffectiveFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TAxionEffectiveFluid_PerturbedStressEnergy

    end module DarkEnergyFluid
