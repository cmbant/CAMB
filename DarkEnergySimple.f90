
    module DarkEnergyModule
    use DarkEnergyInterface
    implicit none

    type, extends(TDarkEnergyBase) :: TDarkEnergy
        real(dl) :: cs2_lam = 1_dl
        !comoving sound speed. Always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)

        logical :: w_perturb = .true.
        !If you are tempted to set this = .false. read
        ! http://cosmocoffee.info/viewtopic.php?t=811
        ! http://cosmocoffee.info/viewtopic.php?t=512

        ! for output and derivs
        real(dl), private, dimension(:), allocatable :: clxq, vq
    contains
    procedure :: ReadParams => TDarkEnergy_ReadParams
    procedure :: Init =>TDarkEnergy_Init
    procedure :: Init_Background => TDarkEnergy_Init_Background
    procedure :: dtauda_Add_Term => TDarkEnergy_dtauda_Add_Term
    procedure :: SetupIndices => TDarkEnergy_SetupIndices
    procedure :: PrepYout => TDarkEnergy_PrepYout
    procedure :: OutputPreMassiveNu => TDarkEnergy_OutputPreMassiveNu
    procedure :: PrepDerivs => TDarkEnergy_PrepDerivs
    procedure :: InitializeYfromVec => TDarkEnergy_InitializeYfromVec
    procedure :: DerivsPrep => TDarkEnergy_DerivsPrep
    procedure :: DerivsAddPostSigma => TDarkEnergy_DerivsAddPostSigma
    procedure :: DerivsAdd2Gpres => TDarkEnergy_DerivsAdd2Gpres
    final :: TDarkEnergy_Finalize
    end type TDarkEnergy

    contains


    subroutine TDarkEnergy_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergy), intent(inout) :: this
    type(TIniFile), intent(in) :: Ini

    this%w_lam = Ini%Read_Double('w', -1.d0)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)

    end subroutine TDarkEnergy_ReadParams


    subroutine TDarkEnergy_Init(this)
    use ModelParams
    class(TDarkEnergy), intent(inout) :: this
    integer :: numThreads

    numThreads = ThreadNum
    call GetNumThreads(numThreads)
    if (.not. allocated(this%clxq)) allocate(this%clxq(0:numThreads), &
        this%vq(0:numThreads))

    end subroutine TDarkEnergy_Init

    subroutine TDarkEnergy_Init_Background(this)
    use GaugeInterface
    class(TDarkEnergy) :: this

    !This is only called once per model, and is a good point to do any extra initialization.
    !It is called before first call to dtauda, but after
    !massive neutrinos are initialized and after GetOmegak

    ! Set the name to export for equations.
    Eqns_name = 'gauge_inv'

    end  subroutine TDarkEnergy_Init_Background


    !Background evolution
    function TDarkEnergy_dtauda_Add_Term(this, a) result(grhoa2)
    !get additional terms for d tau / d a
    use precision
    use ModelParams
    use MassiveNu
    implicit none
    class (TDarkEnergy), intent(in) :: this
    real(dl), intent(IN) :: a
    real(dl) :: rhonu, grhoa2
    integer :: nu_i

    grhoa2 = 0._dl
    if (this%w_lam == -1._dl) then
        grhoa2 = grhov * a ** 4
    else
        grhoa2 = grhov * a ** (1 - 3 * this%w_lam)
    end if
    if (CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a * nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * grhormass(nu_i)
        end do
    end if

    end function TDarkEnergy_dtauda_Add_Term


    subroutine TDarkEnergy_SetupIndices(this, w_ix, neq, maxeq)
    class(TDarkEnergy), intent(in) :: this
    integer, intent(inout) :: w_ix, neq, maxeq

    if (this%w_lam /= -1 .and. this%w_Perturb) then
        w_ix = neq+1
        neq = neq + 2
        maxeq = maxeq + 2
    else
        w_ix = 0
    end if

    end subroutine TDarkEnergy_SetupIndices


    subroutine TDarkEnergy_PrepYout(this, w_ix_out, w_ix, yout, y)
    class(TDarkEnergy), intent(in) :: this
    integer, intent(in) :: w_ix_out, w_ix
    real(dl), intent(inout) :: yout(:)
    real(dl), intent(in) :: y(:)

    if (this%w_lam /= -1 .and. this%w_Perturb) then
        yout(w_ix_out) = y(w_ix)
        yout(w_ix_out+1) = y(w_ix+1)
    end if

    end subroutine TDarkEnergy_PrepYout


    subroutine TDarkEnergy_OutputPreMassiveNu(this, grhov_t, grho, gpres, dgq, dgrho, &
        a, grhov, grhob_t, grhoc_t, grhor_t, grhog_t)
    class(TDarkEnergy), intent(inout) :: this
    real(dl), intent(inout) :: grhov_t, grho, gpres, dgq, dgrho
    real(dl), intent(in) :: a, grhov, grhob_t, grhoc_t, grhor_t, grhog_t

    grhov_t = grhov * a ** (-1 - 3 * this%w_lam)
    grho = grhob_t + grhoc_t + grhor_t + grhog_t + grhov_t
    gpres = (grhog_t + grhor_t) / 3 + grhov_t * this%w_lam

    end subroutine TDarkEnergy_OutputPreMassiveNu


    subroutine TDarkEnergy_PrepDerivs(this, dgrho, dgq, &
        grhov_t, y, w_ix)
    use ModelParams
    class(TDarkEnergy), intent(inout) :: this
    real(dl), intent(inout) :: dgrho, dgq
    real(dl), intent(in) :: grhov_t, y(:)
    integer, intent(in) :: w_ix
    integer :: threadID

    if (this%w_lam /= -1 .and. this%w_Perturb) then
        threadID = GetThreadID()
        this%clxq(threadID) = y(w_ix)
        this%vq(threadID) = y(w_ix + 1)
        dgrho = dgrho + this%clxq(threadID) * grhov_t
        dgq = dgq + this%vq(threadID) * grhov_t * (1 + this%w_lam)
    end if

    end subroutine TDarkEnergy_PrepDerivs


    subroutine TDarkEnergy_InitializeYfromVec(this, y, EV_w_ix, InitVecAti_clxq, &
        InitVecAti_vq)
    class(TDarkEnergy), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    integer, intent(in) :: EV_w_ix
    real(dl), intent(in) :: InitVecAti_clxq, InitVecAti_vq

    if (this%w_lam /= -1 .and. this%w_Perturb) then
        y(EV_w_ix) = InitVecATi_clxq
        y(EV_w_ix + 1) = InitVecAti_vq
    end if

    end subroutine TDarkEnergy_InitializeYfromVec


    subroutine TDarkEnergy_DerivsPrep(this, grhov_t, &
        a, grhov, grhor_t, grhog_t, gpres)
    class(TDarkEnergy), intent(inout) :: this
    real(dl), intent(inout) :: grhov_t
    real(dl), intent(inout), optional :: gpres
    real(dl), intent(in) :: a, grhov, grhor_t, grhog_t

    if (this%w_lam == -1._dl) then
        grhov_t = grhov * a * a
    else
        grhov_t = grhov * a ** (-1 - 3 * this%w_lam)
    end if
    if (present(gpres)) gpres = 0

    end subroutine TDarkEnergy_DerivsPrep


    subroutine TDarkEnergy_DerivsAddPostSigma(this, ayprime, w_ix, &
        adotoa, k, z)
    use ModelParams
    class(TDarkEnergy), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: adotoa, k, z
    integer, intent(in) :: w_ix
    integer :: threadID

     if (this%w_lam /= -1 .and. this%w_Perturb) then
        threadID = GetThreadID()
        ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - this%w_lam) * &
			(this%clxq(threadID) + 3 * adotoa * (1 + this%w_lam) * this%vq(threadID) / k) - &
            (1 + this%w_lam) * k * this%vq(threadID) - (1 + this%w_lam) * k * z

        ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * this%vq(threadID) + &
			k * this%cs2_lam * this%clxq(threadID) / (1 + this%w_lam)
    end if

    end subroutine TDarkEnergy_DerivsAddPostSigma


    function TDarkEnergy_DerivsAdd2Gpres(this, grhog_t, grhor_t, grhov_t) result(gpres)
    class(TDarkEnergy), intent(in) :: this
    real(dl), intent(in) :: grhog_t, grhor_t, grhov_t
    real(dl) :: gpres

    gpres = (grhog_t + grhor_t) / 3 + grhov_t * this%w_lam

    end function TDarkEnergy_DerivsAdd2Gpres


    subroutine TDarkEnergy_Finalize(this)
    type(TDarkEnergy), intent(inout) :: this

    if (allocated(this%clxq)) deallocate(this%clxq, this%vq)
    end subroutine TDarkEnergy_Finalize

    end module DarkEnergyModule
