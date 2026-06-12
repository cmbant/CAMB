    module SPkNonLinear
! Wrapper non-linear model: apply SP(k) suppression on top of a base non-linear model.
    use precision
    use results
    use NonLinear, only: THalofit, halofit_mead, halofit_mead2015, halofit_mead2016, halofit_mead2020_feedback
    use config
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    real(dl), parameter :: SPk_calibrated_z_min = 0.0_dl
    real(dl), parameter :: SPk_calibrated_z_max = 3.0_dl
    real(dl), parameter :: SPk_calibrated_k_min = 1e-12_dl
    real(dl), parameter :: SPk_calibrated_k_max = 12.0_dl
    real(dl), parameter :: SPk_log_calibrated_k_min = log(SPk_calibrated_k_min)
    real(dl), parameter :: SPk_log_calibrated_k_max = log(SPk_calibrated_k_max)
    real(dl), parameter :: SPk_inv_log_10 = 1.0_dl / log(10.0_dl)
    real(dl), parameter :: SPk_min_suppression = 1e-6_dl

    integer, parameter :: SPk_rel_power_law = 1
    integer, parameter :: SPk_rel_cosmo_power_law = 2
    integer, parameter :: SPk_rel_double_power_law = 3

    real(dl), parameter :: SPk_limit_z_nodes(6) = [0.0_dl, 0.125_dl, 0.5_dl, 1.0_dl, 2.0_dl, 3.0_dl]

    real(dl), parameter :: SPk_limit_min_x0_200(6) = [63.59373179416563_dl, 59.88726319810792_dl, &
        56.365373020207954_dl, 39.64033211739476_dl, 91.48777680660496_dl, 45.013496639467114_dl]
    real(dl), parameter :: SPk_limit_min_x1_200(6) = [-9.731727022847117_dl, -9.176876134517682_dl, &
        -8.677101127391419_dl, -6.141984569473165_dl, -14.545324008239655_dl, -7.155837194116757_dl]
    real(dl), parameter :: SPk_limit_min_x2_200(6) = [0.36717360571115487_dl, 0.34646698848026913_dl, &
        0.32901138538950075_dl, 0.23315608004243354_dl, 0.5737424941339003_dl, 0.280547072910215_dl]

    real(dl), parameter :: SPk_limit_max_x0_200(6) = [18.37309024305502_dl, 19.39185050748328_dl, &
        23.581135709114342_dl, 9.921243320018773_dl, 5.886510333331434_dl, 3.2505595228962156_dl]
    real(dl), parameter :: SPk_limit_max_x1_200(6) = [-2.759876785639794_dl, -2.9057029039529163_dl, &
        -3.5325526537664476_dl, -1.4680671748329108_dl, -0.8608020248618151_dl, -0.4637637299766807_dl]
    real(dl), parameter :: SPk_limit_max_x2_200(6) = [0.10310665850449997_dl, 0.10833932967483784_dl, &
        0.13187267761804813_dl, 0.05397858264481507_dl, 0.031183646792386614_dl, 0.016278390350696954_dl]

    real(dl), parameter :: SPk_limit_min_x0_500(6) = [75.5677443552107_dl, 79.78651701865594_dl, &
        12.453646158973672_dl, 119.15086918977593_dl, -44.84238282598863_dl, -7.6276638395020955_dl]
    real(dl), parameter :: SPk_limit_min_x1_500(6) = [-11.520848743752344_dl, -12.198643828164334_dl, &
        -2.0331966873703737_dl, -18.404036176901812_dl, 7.0469849096252615_dl, 1.391204616961_dl]
    real(dl), parameter :: SPk_limit_min_x2_500(6) = [0.4340110697019119_dl, 0.4611079754178873_dl, &
        0.07774265268071975_dl, 0.705717707613519_dl, -0.2811444520707123_dl, -0.06643766025515274_dl]

    real(dl), parameter :: SPk_limit_max_x0_500(6) = [29.7098168708338_dl, 34.60582857431024_dl, &
        43.06158434616325_dl, 20.02440685404008_dl, 10.119926258411336_dl, 7.046905277066627_dl]
    real(dl), parameter :: SPk_limit_max_x1_500(6) = [-4.468289444788325_dl, -5.194226654546699_dl, &
        -6.461887652137118_dl, -2.9536733023513286_dl, -1.4994660666136326_dl, -1.020435592470318_dl]
    real(dl), parameter :: SPk_limit_max_x2_500(6) = [0.16735646291136727_dl, 0.1942560839764333_dl, &
        0.2418483292735806_dl, 0.10842596379693702_dl, 0.055362809275934344_dl, 0.0367046507001731_dl]

    type, extends(TNonLinearModel) :: TSPkNonLinear
        class(TNonLinearModel), allocatable :: BaseModel
        logical :: SPk_feedback = .false.
        integer :: SPk_SO = 200
        integer :: SPk_relation_kind = SPk_rel_power_law
        real(dl) :: SPk_fb_a = 1.0_dl
        real(dl) :: SPk_fb_pow = 0.0_dl
        real(dl) :: SPk_fb_pivot = 1.0_dl
        real(dl) :: SPk_alpha = 0.0_dl
        real(dl) :: SPk_beta = 0.0_dl
        real(dl) :: SPk_gamma = 0.0_dl
        real(dl) :: SPk_epsilon = 0.0_dl
        real(dl) :: SPk_m_pivot = 1.0_dl
    contains
    procedure :: Init => TSPkNonLinear_Init
    procedure :: ReadParams => TSPkNonLinear_ReadParams
    procedure :: GetNonLinRatios => TSPkNonLinear_GetNonLinRatios
    procedure, nopass :: PythonClass => TSPkNonLinear_PythonClass
    procedure, nopass :: SelfPointer => TSPkNonLinear_SelfPointer
    procedure, private :: ValidateConfig => TSPkNonLinear_ValidateConfig
    end type TSPkNonLinear

    public TSPkNonLinear

    contains

    pure function SPk_Poly2(x, c0, c1, c2) result(y)
    real(dl), intent(in) :: x, c0, c1, c2
    real(dl) :: y

    y = c2 * x * x + c1 * x + c0

    end function SPk_Poly2

    pure subroutine SPk_GetParams(SO, z, spk_a_out, spk_b_out, spk_g_out, lambda_a, lambda_b, mu_a, mu_b, mu_c, &
        nu_a, nu_b, nu_c)
    integer, intent(in) :: SO
    real(dl), intent(in) :: z
    real(dl), intent(out) :: spk_a_out, spk_b_out, spk_g_out, lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c
    real(dl) :: x

    x = 1.0_dl + z

    if (SO == 500) then
        spk_a_out = SPk_Poly2(x, 14.783423122120318_dl, -0.999062404857228_dl, 0.12062854541689262_dl)
        spk_b_out = SPk_Poly2(x, 14.620528368613265_dl, -0.9136466201011957_dl, 0.10835389086945699_dl)
        spk_g_out = SPk_Poly2(x, 0.9671320682693298_dl, -0.03185388045484575_dl, 0.02650236152450093_dl)
        lambda_a = SPk_Poly2(x, 0.019349810078190303_dl, -0.007410668383424459_dl, 0.0008334762393555539_dl)
        lambda_b = SPk_Poly2(x, 2.9566773924238143_dl, 0.6205340408676114_dl, -0.001928273640110775_dl)
        mu_a = SPk_Poly2(x, 0.715853343781141_dl, -0.19276613600825665_dl, 0.04948240117059147_dl)
        mu_b = SPk_Poly2(x, 3.385355123440431_dl, 0.9658906605139421_dl, -0.06825861100375574_dl)
        mu_c = SPk_Poly2(x, 4.457257708010122_dl, -2.191853871334233_dl, 0.45457701107254733_dl)
        nu_a = SPk_Poly2(x, 478.86477329610375_dl, 429.88795783439946_dl, 249.25655627821902_dl)
        nu_b = SPk_Poly2(x, -11.227459319819815_dl, -0.5581080204509223_dl, 0.4489962047114509_dl)
        nu_c = SPk_Poly2(x, 3.499449440557995_dl, -0.08488559389068073_dl, -0.0923847866118189_dl)
    else
        spk_a_out = SPk_Poly2(x, 15.24311120000861_dl, -1.2436699435560352_dl, 0.14837558774401766_dl)
        spk_b_out = SPk_Poly2(x, 14.969187892657688_dl, -1.0993025612653198_dl, 0.12905587245129102_dl)
        spk_g_out = SPk_Poly2(x, 0.8000441576980428_dl, -0.01715621131893159_dl, 0.06131887249968379_dl)
        lambda_a = SPk_Poly2(x, 0.02178116280689233_dl, -0.0077564325654746955_dl, 0.0007915576054589781_dl)
        lambda_b = SPk_Poly2(x, 3.0878286643613437_dl, 0.4529677646796634_dl, 0.001552571083240605_dl)
        mu_a = SPk_Poly2(x, 0.6930259177449359_dl, -0.16913553700233935_dl, 0.04263185199898842_dl)
        mu_b = SPk_Poly2(x, 3.161914061444856_dl, 0.8616834297321924_dl, 0.011346427353554053_dl)
        mu_c = SPk_Poly2(x, 5.532188503256583_dl, -3.0864672185252537_dl, 0.5083422518560442_dl)
        nu_a = SPk_Poly2(x, 413.00988701513904_dl, 311.63957063032285_dl, 37.89105940901369_dl)
        nu_b = SPk_Poly2(x, -11.243859405779181_dl, -0.34421412616421965_dl, 0.3343548325485801_dl)
        nu_c = SPk_Poly2(x, 3.476463891168505_dl, -0.018333059687988575_dl, -0.08276237963970698_dl)
    end if

    end subroutine SPk_GetParams

    pure function SPk_AkimaInterp(z, x_nodes, y_nodes) result(y)
    real(dl), intent(in) :: z
    real(dl), intent(in) :: x_nodes(:), y_nodes(:)
    real(dl) :: y
    integer :: i, n
    real(dl) :: h, s, w1, w2
    ! Fixed sizes for 6 z-nodes: delta(n-1), mext(n+3), t(n).
    real(dl) :: delta(5), mext(9), t(6)
    real(dl) :: h00, h10, h01, h11

    n = size(x_nodes)
    if (z <= x_nodes(1)) then
        y = y_nodes(1)
        return
    end if
    if (z >= x_nodes(n)) then
        y = y_nodes(n)
        return
    end if

    if (n < 2) then
        y = y_nodes(1)
        return
    end if

    do i = 1, n - 1
        delta(i) = (y_nodes(i + 1) - y_nodes(i)) / (x_nodes(i + 1) - x_nodes(i))
    end do

    ! Extend slopes at both ends using Akima endpoint construction.
    mext(3:n + 1) = delta
    mext(2) = 2.0_dl * mext(3) - mext(4)
    mext(1) = 2.0_dl * mext(2) - mext(3)
    mext(n + 2) = 2.0_dl * mext(n + 1) - mext(n)
    mext(n + 3) = 2.0_dl * mext(n + 2) - mext(n + 1)

    do i = 1, n
        w1 = abs(mext(i + 3) - mext(i + 2))
        w2 = abs(mext(i + 1) - mext(i))
        if (w1 + w2 > 0.0_dl) then
            t(i) = (w1 * mext(i + 1) + w2 * mext(i + 2)) / (w1 + w2)
        else
            t(i) = 0.5_dl * (mext(i + 1) + mext(i + 2))
        end if
    end do

    do i = 1, n - 1
        if (z < x_nodes(i + 1)) then
            h = x_nodes(i + 1) - x_nodes(i)
            s = (z - x_nodes(i)) / h
            h00 = (1.0_dl + 2.0_dl * s) * (1.0_dl - s) * (1.0_dl - s)
            h10 = s * (1.0_dl - s) * (1.0_dl - s)
            h01 = s * s * (3.0_dl - 2.0_dl * s)
            h11 = s * s * (s - 1.0_dl)
            y = h00 * y_nodes(i) + h10 * h * t(i) + h01 * y_nodes(i + 1) + h11 * h * t(i + 1)
            return
        end if
    end do

    y = y_nodes(n)

    end function SPk_AkimaInterp

    pure subroutine SPk_ComputeFbFromParams(log_kh, spk_a, spk_b, spk_g, relation_kind, fb_a, fb_pow, fb_pivot, &
        rel_alpha, rel_beta, rel_gamma, rel_epsilon, rel_m_pivot, e_ratio, fb, m_opt)
    integer, intent(in) :: relation_kind
    real(dl), intent(in) :: log_kh, spk_a, spk_b, spk_g, fb_a, fb_pow, fb_pivot
    real(dl), intent(in) :: rel_alpha, rel_beta, rel_gamma, rel_epsilon, rel_m_pivot
    real(dl), intent(in) :: e_ratio
    real(dl), intent(out) :: fb, m_opt
    real(dl) :: best_mass, e_ratio_eff

    best_mass = spk_a - (spk_a - spk_b) * exp(spk_g * log_kh)
    m_opt = 10.0_dl ** best_mass
    e_ratio_eff = max(e_ratio, SPk_calibrated_k_min)

    select case (relation_kind)
    case (SPk_rel_cosmo_power_law)
        fb = (exp(rel_alpha) / 100.0_dl) * ((m_opt / 1.0e14_dl) ** (rel_beta - 1.0_dl)) * &
            (e_ratio_eff ** rel_gamma)
    case (SPk_rel_double_power_law)
        fb = 0.5_dl * rel_epsilon * (((m_opt / rel_m_pivot) ** rel_alpha) + &
            ((m_opt / rel_m_pivot) ** rel_beta)) * (e_ratio_eff ** rel_gamma)
    case default
        fb = fb_a * ((m_opt / fb_pivot) ** fb_pow)
    end select

    end subroutine SPk_ComputeFbFromParams

    pure subroutine SPk_GetFbLimitCoefficients(SO, z, min_c0, min_c1, min_c2, max_c0, max_c1, max_c2)
    integer, intent(in) :: SO
    real(dl), intent(in) :: z
    real(dl), intent(out) :: min_c0, min_c1, min_c2, max_c0, max_c1, max_c2

    if (SO == 500) then
        min_c0 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x0_500)
        min_c1 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x1_500)
        min_c2 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x2_500)
        max_c0 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x0_500)
        max_c1 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x1_500)
        max_c2 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x2_500)
    else
        min_c0 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x0_200)
        min_c1 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x1_200)
        min_c2 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_min_x2_200)
        max_c0 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x0_200)
        max_c1 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x1_200)
        max_c2 = SPk_AkimaInterp(z, SPk_limit_z_nodes, SPk_limit_max_x2_200)
    end if

    end subroutine SPk_GetFbLimitCoefficients

    pure subroutine SPk_GetFbLimitsFromCoefficients(m_halo, min_c0, min_c1, min_c2, max_c0, max_c1, max_c2, &
        min_fb, max_fb)
    real(dl), intent(in) :: m_halo, min_c0, min_c1, min_c2, max_c0, max_c1, max_c2
    real(dl), intent(out) :: min_fb, max_fb
    real(dl) :: logm

    logm = log10(m_halo)
    min_fb = 0.8_dl * 10.0_dl ** (min_c0 + min_c1 * logm + min_c2 * logm * logm)
    max_fb = 1.2_dl * 10.0_dl ** (max_c0 + max_c1 * logm + max_c2 * logm * logm)

    end subroutine SPk_GetFbLimitsFromCoefficients

    pure function SPk_Suppression(log10_kh, fb, lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c) result(sup)
    ! Compute the SP(k) suppression factor for a single (k, z) point.
    ! Shape function uses fitting parameters from Salcido et al. (2023, MNRAS 523, 2247).
    real(dl), intent(in) :: log10_kh, fb
    real(dl), intent(in) :: lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c
    real(dl) :: sup
    real(dl) :: x, x0, x1, x2

    x = log10_kh
    x0 = 1.0_dl + lambda_a * exp(lambda_b * x)
    x1 = mu_a + ((1.0_dl - mu_a) / (1.0_dl + exp(mu_b * x + mu_c)))
    x2 = nu_a * exp(-0.5_dl * ((x - nu_b) / nu_c) ** 2)

    sup = x0 - (x0 - x1) * exp(-x2 * fb)
    sup = max(sup, SPk_min_suppression)

    end function SPk_Suppression

    function TSPkNonLinear_PythonClass()
    character(LEN=:), allocatable :: TSPkNonLinear_PythonClass

    TSPkNonLinear_PythonClass = 'SPkNonLinear'

    end function TSPkNonLinear_PythonClass

    subroutine TSPkNonLinear_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TSPkNonLinear), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSPkNonLinear_SelfPointer

    subroutine TSPkNonLinear_Init(this, State)
    class(TSPkNonLinear) :: this
    class(TCAMBdata), target :: State

    if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)
    this%Min_kh_nonlinear = this%BaseModel%Min_kh_nonlinear
    call this%BaseModel%Init(State)

    end subroutine TSPkNonLinear_Init

    subroutine TSPkNonLinear_ReadParams(this, Ini)
    use IniObjects
    class(TSPkNonLinear) :: this
    class(TIniFile), intent(in) :: Ini

    if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)
    call this%BaseModel%ReadParams(Ini)

    this%SPk_feedback = Ini%Read_Logical('SPk_feedback', .false.)
    this%SPk_SO = Ini%Read_Int('SPk_SO', 200)
    this%SPk_relation_kind = Ini%Read_Int('SPk_relation_kind', SPk_rel_power_law)
    this%SPk_fb_a = Ini%Read_Double('SPk_fb_a', 1.0_dl)
    this%SPk_fb_pow = Ini%Read_Double('SPk_fb_pow', 0.0_dl)
    this%SPk_fb_pivot = Ini%Read_Double('SPk_fb_pivot', 1.0_dl)
    this%SPk_alpha = Ini%Read_Double('SPk_alpha', 0.0_dl)
    this%SPk_beta = Ini%Read_Double('SPk_beta', 0.0_dl)
    this%SPk_gamma = Ini%Read_Double('SPk_gamma', 0.0_dl)
    this%SPk_epsilon = Ini%Read_Double('SPk_epsilon', 0.0_dl)
    this%SPk_m_pivot = Ini%Read_Double('SPk_m_pivot', 1.0_dl)

    end subroutine TSPkNonLinear_ReadParams

    subroutine TSPkNonLinear_ValidateConfig(this)
    class(TSPkNonLinear), intent(in) :: this

    if (this%SPk_SO /= 200 .and. this%SPk_SO /= 500) then
        call MpiStop('SP(k): SPk_SO must be 200 or 500')
    end if
    if (this%SPk_relation_kind < SPk_rel_power_law .or. this%SPk_relation_kind > SPk_rel_double_power_law) then
        call MpiStop('SP(k): SPk_relation_kind must be 1 (power_law), 2 (cosmo_power_law), or 3 (double_power_law)')
    end if

    if (this%SPk_relation_kind == SPk_rel_power_law) then
        if (this%SPk_fb_pivot <= 0.0_dl) call MpiStop('SP(k): SPk_fb_pivot must be > 0 for power_law relation')
    end if
    if (this%SPk_relation_kind == SPk_rel_double_power_law) then
        if (this%SPk_m_pivot <= 0.0_dl) call MpiStop('SP(k): SPk_m_pivot must be > 0 for double_power_law relation')
    end if

    select type (base => this%BaseModel)
    type is (THalofit)
        if (this%SPk_feedback .and. base%halofit_version == halofit_mead2020_feedback) then
            call MpiStop('SP(k) is not compatible with halofit_version=mead2020_feedback')
        end if
        if (this%SPk_feedback .and. &
            (base%halofit_version == halofit_mead .or. base%halofit_version == halofit_mead2015 &
            .or. base%halofit_version == halofit_mead2016) .and. &
            (abs(base%HMcode_A_baryon - 3.13_dl) > 1e-12_dl .or. abs(base%HMcode_eta_baryon - 0.603_dl) > 1e-12_dl)) &
            then
            call MpiStop('SP(k) cannot be combined with HMCode_A_baryon/' &
                // 'HMCode_eta_baryon baryonic corrections in HMCode 2015/2016')
        end if
    class default
        continue
    end select

    end subroutine TSPkNonLinear_ValidateConfig

    subroutine TSPkNonLinear_GetNonLinRatios(this, State, CAMB_Pk)
    class(TSPkNonLinear) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    integer :: itf, i
    real(dl) :: rk, log_rk, log_rk_eval, spk_sup, spk_href, spk_eratio
    real(dl) :: spk_fb, spk_m_opt, spk_fb_min, spk_fb_max
    real(dl) :: spk_a, spk_b, spk_g
    real(dl) :: lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c
    real(dl) :: min_c0, min_c1, min_c2, max_c0, max_c1, max_c2
    logical, save :: warned_spk_z_outside = .false.
    logical, save :: warned_spk_k_clamped = .false.
    logical, save :: warned_spk_fb_outside = .false.

    if (.not. allocated(this%BaseModel)) allocate(THalofit::this%BaseModel)

    call this%ValidateConfig()
    call this%BaseModel%GetNonLinRatios(State, CAMB_Pk)

    if (.not. this%SPk_feedback) return

    select type (State)
    class is (CAMBdata)
        if (this%SPk_relation_kind == SPk_rel_cosmo_power_law &
            .or. this%SPk_relation_kind == SPk_rel_double_power_law) then
            spk_href = State%Hofz(0.3_dl)
        else
            spk_href = 1.0_dl
        end if
    class default
        call MpiStop('SP(k): unsupported state type for Hofz evaluation')
    end select

    do itf = 1, CAMB_Pk%num_z
        if (CAMB_Pk%redshifts(itf) < SPk_calibrated_z_min .or. CAMB_Pk%redshifts(itf) > SPk_calibrated_z_max) then
            if (FeedbackLevel > 0 .and. .not. warned_spk_z_outside) then
                write(*,'(A,F8.3,A,F6.2,A,F6.2,A)') 'WARNING: SP(k) skipped outside calibrated redshift range. z=', &
                    CAMB_Pk%redshifts(itf), ', &
                    calibrated range=[', SPk_calibrated_z_min, ',', SPk_calibrated_z_max, '].'
                warned_spk_z_outside = .true.
            end if
            cycle
        end if
        select type (State)
        class is (CAMBdata)
            if (this%SPk_relation_kind == SPk_rel_cosmo_power_law &
                .or. this%SPk_relation_kind == SPk_rel_double_power_law) then
                spk_eratio = State%Hofz(CAMB_Pk%redshifts(itf)) / spk_href
            else
                spk_eratio = 1.0_dl
            end if
        class default
            call MpiStop('SP(k): unsupported state type for Hofz evaluation')
        end select
        call SPk_GetParams(this%SPk_SO, CAMB_Pk%redshifts(itf), spk_a, spk_b, spk_g, lambda_a, lambda_b, &
            mu_a, mu_b, mu_c, nu_a, nu_b, nu_c)
        call SPk_GetFbLimitCoefficients(this%SPk_SO, CAMB_Pk%redshifts(itf), min_c0, min_c1, min_c2, &
            max_c0, max_c1, max_c2)
        do i = 1, CAMB_Pk%num_k
            log_rk = CAMB_Pk%log_kh(i)
            if (log_rk < SPk_log_calibrated_k_min) cycle
            if (log_rk > SPk_log_calibrated_k_max) then
                if (FeedbackLevel > 0 .and. .not. warned_spk_k_clamped) then
                    rk = exp(log_rk)
                    write(*,'(A,F8.3,A,F6.2,A)') &
                        'WARNING: SP(k) input k exceeds calibrated range; clamping to k_max=', &
                        rk, ' -> ', SPk_calibrated_k_max, ' h/Mpc.'
                    warned_spk_k_clamped = .true.
                end if
            end if
            log_rk_eval = min(log_rk, SPk_log_calibrated_k_max)
            call SPk_ComputeFbFromParams(log_rk_eval, spk_a, spk_b, spk_g, this%SPk_relation_kind, &
                this%SPk_fb_a, this%SPk_fb_pow, this%SPk_fb_pivot, this%SPk_alpha, this%SPk_beta, &
                this%SPk_gamma, this%SPk_epsilon, this%SPk_m_pivot, spk_eratio, spk_fb, spk_m_opt)
            call SPk_GetFbLimitsFromCoefficients(spk_m_opt, min_c0, min_c1, min_c2, max_c0, max_c1, max_c2, &
                spk_fb_min, spk_fb_max)
            if (spk_fb < spk_fb_min .or. spk_fb > spk_fb_max) then
                if (FeedbackLevel > 0 .and. .not. warned_spk_fb_outside) then
                    write(*,'(A)') 'WARNING: SP(k) baryon fraction outside calibrated fitting limits; ' &
                        // 'setting out-of-range points to NaN.'
                    warned_spk_fb_outside = .true.
                end if
                CAMB_Pk%nonlin_ratio(i, itf) = ieee_value(1.0_dl, ieee_quiet_nan)
                cycle
            end if
            spk_sup = SPk_Suppression(log_rk_eval * SPk_inv_log_10, spk_fb, &
                lambda_a, lambda_b, mu_a, mu_b, mu_c, nu_a, nu_b, nu_c)
            CAMB_Pk%nonlin_ratio(i, itf) = CAMB_Pk%nonlin_ratio(i, itf) * sqrt(spk_sup)
        end do
    end do

    end subroutine TSPkNonLinear_GetNonLinRatios

    end module SPkNonLinear
