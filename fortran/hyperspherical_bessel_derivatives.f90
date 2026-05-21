    module HypersphericalBesselDerivatives
    use Precision
    use results
    use SpherBessels, only: phi_langer
    use OlverHypersphericalBessel, only: u_olver
    implicit none
    private

    public :: USpherBesselWithDeriv

    contains

    subroutine USpherBesselWithDeriv(closed,CP,Chi,l,beta,y1,y2)
    !returns y1=ujl*sinhChi and y2=diff(y1,Chi)
    !Olver should be accurate to 1e-4 (falls back to recursive if needed)
    real(dl) Chi,beta,y1,y2
    real(dl) sin_K, cot_K
    integer l,K
    logical, intent(IN) :: closed
    Type(CAMBParams), intent(in) :: CP

    if (closed) then
        cot_K = 1._dl / tan(Chi)
        K = 1
    else
        cot_K = 1._dl / tanh(Chi)
        K = -1
    end if

    ! if ((l > 60 * (1 + AccuracyTarget) * CP%Accuracy%AccuracyBoost) .and. &
    !     ((.not. closed) .or. (abs(Chi - const_pi / 2) > 0.2d0))) then
    !     y1 = phi_langer(l, K, beta, Chi) * sinhChi
    !     y2 = y1 * (l + 1) * cothChi
    !     if (.not. closed .or. (l + 1 < nint(beta))) y2 = y2 - &
    !         sqrt(beta**2 - (K * (l + 1)**2)) * phi_langer(l + 1, K, beta, Chi) * sinhChi
    ! else
    y1 = u_olver(l, K, beta, Chi)
    y2 = y1 * (l + 1) * cot_K
    if (.not. closed .or. (l + 1 < nint(beta))) y2 = y2 - &
        sqrt(beta**2 - (K * (l + 1)**2)) * u_olver(l + 1, K, beta, Chi)
    !end if

    end subroutine USpherBesselWithDeriv

    end module HypersphericalBesselDerivatives
