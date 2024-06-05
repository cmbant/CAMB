    module DarkEnergyPPF
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyPPF
        real(dl) :: c_Gamma_ppf = 0.4_dl
    contains
    procedure :: ReadParams => TDarkEnergyPPF_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyPPF_PythonClass
    procedure :: Init => TDarkEnergyPPF_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyPPF_PerturbedStressEnergy
    procedure :: diff_rhopi_Add_Term => TDarkEnergyPPF_diff_rhopi_Add_Term
    procedure, nopass :: SelfPointer => TDarkEnergyPPF_SelfPointer
    procedure, private :: setcgammappf
    end type TDarkEnergyPPF

    public TDarkEnergyPPF
    contains

    subroutine TDarkEnergyPPF_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyPPF) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    if (this%cs2_lam /= 1._dl) error stop 'cs2_lam not supported by PPF model'
    call this%setcgammappf

    end subroutine TDarkEnergyPPF_ReadParams

    function TDarkEnergyPPF_PythonClass()
    character(LEN=:), allocatable :: TDarkEnergyPPF_PythonClass

    TDarkEnergyPPF_PythonClass = 'DarkEnergyPPF'
    end function TDarkEnergyPPF_PythonClass


    subroutine TDarkEnergyPPF_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyPPF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyPPF_SelfPointer

    subroutine TDarkEnergyPPF_Init(this, State)
    use classes
    use config
    class(TDarkEnergyPPF), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    call this%TDarkEnergyEqnOfState%Init(State)
    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 1
    end if
    if (this%cs2_lam /= 1._dl) &
        call GlobalError('DarkEnergyPPF does not support varying sound speed',error_unsupported_params)

    end subroutine TDarkEnergyPPF_Init

    subroutine setcgammappf(this)
    class(TDarkEnergyPPF) :: this

    this%c_Gamma_ppf = 0.4_dl * sqrt(this%cs2_lam)

    end subroutine setcgammappf


    function TDarkEnergyPPF_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    !Get derivative of anisotropic stress
    class(TDarkEnergyPPF), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh

    if (this%is_cosmological_constant .or. this%no_perturbations) then
        ppiedot = 0
    else
        hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
        ppiedot = 3._dl * dgrhoe + dgqe * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
        ppiedot = ppiedot * adotoa / Kf1
    end if

    end function TDarkEnergyPPF_diff_rhopi_Add_Term


    subroutine TDarkEnergyPPF_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
    real(dl) :: vT, grhoT, k2

    if (this%no_perturbations) then
        dgrhoe=0
        dgqe=0
        return
    end if

    k2=k**2
    !ppf
    grhoT = grho - grhov_t
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = this%c_Gamma_ppf * k / adotoa

    if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 ?????????
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe = S_Gamma - Gammadot / adotoa - Gamma
    dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
    dgrhoe = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe

    end subroutine TDarkEnergyPPF_PerturbedStressEnergy



    end module DarkEnergyPPF
