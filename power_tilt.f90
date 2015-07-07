    !This module provides the initial power spectra, parameterized as an expansion in ln k
    !
    ! ln P_s = ln A_s + (n_s -1)*ln(k/k_0_scalar) + n_{run}/2 * ln(k/k_0_scalar)^2 + n_{runrun}/6 * ln(k/k_0_scalar)^3
    !
    ! so if n_{run} = 0, n_{runrun}=0
    !
    ! P_s = A_s (k/k_0_scalar)^(n_s-1)
    !
    !for the scalar spectrum, when n_s=an(in) is the in'th spectral index. k_0_scalar
    !is a pivot scale, fixed here to 0.05/Mpc (change it below as desired or via .ini file).
    !
    !The tensor spectrum has three different supported parameterizations giving
    !
    ! ln P_t = ln A_t + n_t*ln(k/k_0_tensor) + n_{t,run}/2 * ln(k/k_0_tensor)^2
    !
    ! tensor_parameterization==tensor_param_indeptilt (=1) (default, same as CAMB pre-April 2014)
    !
    ! A_t = r A_s
    !
    ! tensor_parameterization==tensor_param_rpivot (=2)
    !
    ! A_t = r P_s(k_0_tensor)
    !
    ! tensor_parameterization==tensor_param_AT (=3)
    !
    ! A_t =  tensor_amp
    !
    !The absolute normalization of the Cls is unimportant here, but the relative ratio
    !of the tensor and scalar Cls generated with this module will be correct for general models
    !
    !December 2003 - changed default tensor pivot to 0.05 (consistent with CMBFAST 4.5)
    !April 2014 added different tensor parameterizations, running of running and running of tensors

    module InitialPower
    use Precision
    use AMLutils
    implicit none

    private

    character(LEN=*), parameter :: Power_Name = 'power_tilt'

    integer, parameter :: nnmax= 5
    !Maximum possible number of different power spectra to use

    integer, parameter, public :: tensor_param_indeptilt=1,  tensor_param_rpivot = 2, tensor_param_AT = 3

    Type InitialPowerParams
        integer :: tensor_parameterization = tensor_param_indeptilt
        integer nn  !Must have this variable
        !The actual number of power spectra to use

        !For the default implementation return power spectra based on spectral indices
        real(dl) an(nnmax) !scalar spectral indices
        real(dl) n_run(nnmax) !running of spectral index
        real(dl) n_runrun(nnmax) !running of spectral index
        real(dl) ant(nnmax) !tensor spectral indices
        real(dl) nt_run(nnmax) !tensor spectral index running
        real(dl) rat(nnmax) !ratio of scalar to tensor initial power spectrum amplitudes
        real(dl) k_0_scalar, k_0_tensor !pivot scales
        real(dl) ScalarPowerAmp(nnmax)
        real(dl) TensorPowerAmp(nnmax) !A_T at k_0_tensor if tensor_parameterization==tensor_param_AT
    end Type InitialPowerParams

    real(dl) curv  !Curvature contant, set in InitializePowers

    Type(InitialPowerParams), save :: P

    !Make things visible as neccessary...

    public InitialPowerParams, InitialPower_ReadParams, InitializePowers, ScalarPower, TensorPower
    public nnmax,Power_Descript, Power_Name, SetDefPowerParams
    contains


    subroutine SetDefPowerParams(AP)
    Type (InitialPowerParams) :: AP

    AP%nn     = 1 !number of initial power spectra
    AP%an     = 1 !scalar spectral index
    AP%n_run   = 0 !running of scalar spectral index
    AP%n_runrun   = 0 !running of running of scalar spectral index
    AP%ant    = 0 !tensor spectra index
    AP%nt_run   = 0 !running of tensor spectral index
    AP%rat    = 1
    AP%k_0_scalar = 0.05
    AP%k_0_tensor = 0.05
    AP%ScalarPowerAmp = 1
    AP%TensorPowerAmp = 1
    AP%tensor_parameterization = tensor_param_indeptilt

    end subroutine SetDefPowerParams

    subroutine InitializePowers(AParamSet,acurv)
    Type (InitialPowerParams) :: AParamSet
    !Called before computing final Cls in cmbmain.f90
    !Could read spectra from disk here, do other processing, etc.

    real(dl) acurv

    if (AParamSet%nn > nnmax) then
        write (*,*) 'To use ',AParamSet%nn,'power spectra you need to increase'
        write (*,*) 'nnmax in power_tilt.f90, currently ',nnmax
    end if
    P = AParamSet

    curv=acurv

    !Write implementation specific code here...

    end subroutine InitializePowers


    function ScalarPower(k,ix)

    !"ix" gives the index of the power to return for this k
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + curv)/|curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.

    real(dl) ScalarPower,k, lnrat
    integer ix

    lnrat = log(k/P%k_0_scalar)
    ScalarPower=P%ScalarPowerAmp(ix)*exp(lnrat*( P%an(ix)-1 + lnrat*(P%n_run(ix)/2 + P%n_runrun(ix)/6*lnrat)))

    !         ScalarPower = ScalarPower * (1 + 0.1*cos( lnrat*30 ) )

    end function ScalarPower


    function TensorPower(k,ix)

    !TensorPower= const for scale invariant spectrum
    !The normalization is defined so that
    ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
    !for a closed model
    ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
    !for an open model
    !"ix" gives the index of the power spectrum to return
    !Here nu^2 = (k^2 + 3*curv)/|curv|

    real(dl) TensorPower,k
    real(dl), parameter :: PiByTwo=3.14159265d0/2._dl
    integer ix
    real(dl) lnrat, k_dep

    lnrat = log(k/P%k_0_tensor)
    k_dep = exp(lnrat*(P%ant(ix) + P%nt_run(ix)/2*lnrat))
    if (P%tensor_parameterization==tensor_param_indeptilt) then
        TensorPower = P%rat(ix)*P%ScalarPowerAmp(ix)*k_dep
    else if (P%tensor_parameterization==tensor_param_rpivot) then
        TensorPower = P%rat(ix)*ScalarPower(P%k_0_tensor,ix) * k_dep
    else if (P%tensor_parameterization==tensor_param_AT) then
        TensorPower = P%TensorPowerAmp(ix) * k_dep
    end if
    if (curv < 0) TensorPower=TensorPower*tanh(PiByTwo*sqrt(-k**2/curv-3))

    end function TensorPower

    !Get parameters describing parameterisation (for FITS file)
    !Does not support running extensions
    function Power_Descript(in, Scal, Tens, Keys, Vals)
    character(LEN=8), intent(out) :: Keys(*)
    real(dl), intent(out) :: Vals(*)
    integer, intent(IN) :: in
    logical, intent(IN) :: Scal, Tens
    integer num, Power_Descript
    num=0
    if (Scal) then
        num=num+1
        Keys(num) = 'n_s'
        Vals(num) = P%an(in)
        num=num+1
        Keys(num) = 'n_run'
        Vals(num) = P%n_run(in)
        num=num+1
        Keys(num) = 's_pivot'
        Vals(num) = P%k_0_scalar
    end if
    if (Tens) then
        num=num+1
        Keys(num) = 'n_t'
        Vals(num) = P%ant(in)
        num=num+1
        Keys(num) = 't_pivot'
        Vals(num) = P%k_0_tensor
        if (Scal) then
            num=num+1
            Keys(num) = 'p_ratio'
            Vals(num) = P%rat(in)
        end if
    end if
    Power_Descript = num

    end  function Power_Descript

    subroutine InitialPower_ReadParams(InitPower, Ini, WantTensors)
    use IniFile
    Type(InitialPowerParams) :: InitPower
    Type(TIniFile) :: Ini
    logical, intent(in) :: WantTensors
    integer i

    InitPower%k_0_scalar = Ini_Read_Double_File(Ini,'pivot_scalar',InitPower%k_0_scalar)
    InitPower%k_0_tensor = Ini_Read_Double_File(Ini,'pivot_tensor',InitPower%k_0_tensor)
    InitPower%nn = Ini_Read_Int_File(Ini,'initial_power_num',1)
    if (InitPower%nn>nnmax) call MpiStop('Too many initial power spectra - increase nnmax in InitialPower')
    if (WantTensors) then
        InitPower%tensor_parameterization =  Ini_Read_Int_File(Ini, 'tensor_parameterization',tensor_param_indeptilt)
        if (InitPower%tensor_parameterization < tensor_param_indeptilt .or. &
            & InitPower%tensor_parameterization > tensor_param_AT) &
            & call MpiStop('InitialPower: unknown tensor_parameterization')
    end if
    InitPower%rat(:) = 1
    do i=1, InitPower%nn
        InitPower%an(i) = Ini_Read_Double_Array_File(Ini,'scalar_spectral_index', i)
        InitPower%n_run(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrun',i,0._dl)
        InitPower%n_runrun(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrunrun',i,0._dl)

        if (WantTensors) then
            InitPower%ant(i) = Ini_Read_Double_Array_File(Ini,'tensor_spectral_index',i)
            InitPower%nt_run(i) = Ini_Read_Double_Array_File(Ini,'tensor_nrun',i,0._dl)
            if (InitPower%tensor_parameterization == tensor_param_AT) then
                InitPower%TensorPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'tensor_amp',i)
            else
                InitPower%rat(i) = Ini_Read_Double_Array_File(Ini,'initial_ratio',i)
            end if
        end if

        InitPower%ScalarPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'scalar_amp',i,1._dl)
        !Always need this as may want to set tensor amplitude even if scalars not computed
    end do

    end  subroutine InitialPower_ReadParams


    end module InitialPower
