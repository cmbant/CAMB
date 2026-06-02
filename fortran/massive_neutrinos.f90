    module MassiveNu
    use constants
    implicit none
    private

    real(dl), parameter  :: fermi_dirac_const  = 7._dl/120*const_pi**4 ! 5.68219698_dl
    !fermi_dirac_const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7._dl/const_pi**2   !0.072372274_dl

    ! Smallest/largest a*m_nu values using the rho/P direct fits rather than series expansions.
    real(dl), parameter :: nu_rhop_am_min = 0.42_dl
    real(dl), parameter :: nu_rhop_am_max = 70._dl

    ! Cut values used by the D = rho - 3P fit in ThermalNuBackground_drho.
    real(dl), parameter :: am_min = 0.3_dl
    real(dl), parameter :: am_max = 70._dl

    ! Direct smooth fit for thermal neutrino rho/P over nu_rhop_am_min < am < nu_rhop_am_max.
    ! z maps [nu_rhop_am_min,nu_rhop_am_max] to [-1,1].
    real(dl), parameter :: nu_rhop_fit_c = sqrt(nu_rhop_am_min*nu_rhop_am_max)
    real(dl), parameter :: nu_rhop_fit_inv_zmax = &
        (nu_rhop_am_max + nu_rhop_fit_c)/(nu_rhop_am_max - nu_rhop_fit_c)

    ! D-fit mapping over am_min < am < am_max.
    real(dl), parameter :: nu_fit_c = sqrt(am_min*am_max)
    real(dl), parameter :: nu_fit_inv_zmax = (am_max + nu_fit_c)/(am_max - nu_fit_c)

    ! Large-am scaling coefficients:
    ! rho ~ nu_fit_rho_scale * am
    ! P   ~ (1/nu_fit_p_denom_scale) / am
    real(dl), parameter :: nu_fit_rho_scale = 3._dl*zeta3/(2._dl*fermi_dirac_const)
    real(dl), parameter :: nu_fit_p_denom_scale = fermi_dirac_const/((900._dl/120._dl)*zeta5)

    ! Power coefficients in z, increasing order. rho is fitted for absolute error;
    ! P is fitted for relative error. The rho/P cut points are separate from
    ! the D-fit cut points so drho keeps the original fit below.
    ! Max relative errors over 1e-4 <= am <= 1e4 are about 2.4e-5 in rho
    ! and 4.7e-5 in P against direct quadrature.
    real(dl), parameter :: nu_fit_rho_c(0:5) = (/ &
        7.499710425926954249e-01_dl,  1.188066888581748443e-01_dl, &
        1.545077526625578401e-01_dl, -8.378996019212282820e-02_dl, &
        2.123075360029347963e-02_dl, -2.544964203886362908e-03_dl /)

    real(dl), parameter :: nu_fit_p_c(0:6) = (/ &
        9.283417043601330798e-01_dl,  3.156134695375645838e-01_dl, &
        -2.435843189264263742e-01_dl, -2.509793530357270000e-02_dl, &
        4.067958090718486880e-02_dl,  2.378742143669706002e-03_dl, &
        -1.965863329520473827e-03_dl /)

    ! Low-am polynomial replacing the small-am logarithmic expansion in rho/P.
    ! rho = 1 + const2*am**2 + am**4*(c0 + c1*am**2)
    ! P   = (1 - const2*am**2)/3 + am**4*(c0 + c1*am**2)
    ! Coefficients are chosen to match the mid fit at nu_rhop_am_min.
    real(dl), parameter :: nu_low_rho_c(0:1) = (/ &
        -1.872929199948838302e-02_dl,  1.831728681837861694e-02_dl /)
    real(dl), parameter :: nu_low_p_c(0:1) = (/ &
        1.667044572156534468e-02_dl, -2.277471498716537868e-02_dl /)

    ! ~1e-4 fit for D = rho - 3P over am_min < am < am_max.
    ! D(am) = am**2/(1 + 2*am) * poly(D_c,z)
    real(dl), parameter :: nu_fit_D_c(0:5) = (/ &
        5.800497880734221123e-01_dl,  1.801974470841860576e-01_dl, &
        -1.602709049602861757e-01_dl,  3.167230341650149189e-02_dl, &
        1.056591590995091534e-02_dl, -3.877662402660121861e-03_dl /)

    Type TNuPerturbations
        !Sample for massive neutrino momentum
        !Default settings appear to be OK for P_k accuate at 1e-3 level
        integer nqmax !actual number of q modes evolves
        real(dl), allocatable ::  nu_q(:), nu_int_kernel(:)
    contains
    procedure :: init => TNuPerturbations_init
    end type TNuPerturbations

    Type TThermalNuBackground
        ! Thermal massive-neutrino background. rho/P are direct fit evaluations;
        ! target_rho is only temporary state for the mass-inversion root solve.
        real(dl), private :: target_rho
    contains
    procedure :: rho_P => ThermalNuBackground_rho_P
    procedure :: rho => ThermalNuBackground_rho
    procedure :: drho => ThermalNuBackground_drho
    procedure :: find_nu_mass_for_rho => ThermalNuBackground_find_nu_mass_for_rho
    end type TThermalNuBackground

    Type(TThermalNuBackground) :: ThermalNuBackground

    public fermi_dirac_const,  sum_mnu_for_m1, neutrino_mass_fac, TNuPerturbations, &
        ThermalNuBackground
    contains

    pure subroutine sum_mnu_for_m1(summnu,dsummnu, m1, targ, sgn)
    use constants
    real(dl), intent(in) :: m1, targ, sgn
    real(dl), intent(out) :: summnu, dsummnu
    real(dl) :: m2,m3

    m2 = sqrt(m1**2 + delta_mnu21)
    m3 = sqrt(m1**2 + sgn*delta_mnu31)
    summnu = m1 + m2 + m3 - targ
    dsummnu = m1/m2+m1/m3 + 1

    end subroutine sum_mnu_for_m1

    subroutine TNuPerturbations_init(this,Accuracy)
    !Set up which momenta to integrate the neutrino perturbations, depending on accuracy
    !Using three optimized momenta works very well in most cases
    class(TNuPerturbations) :: this
    real(dl), intent(in) :: Accuracy
    real(dl) :: dq,dlfdlq, q
    integer i

    if (Accuracy > 3) then
        this%nqmax = nint(Accuracy*10)
    else if (Accuracy > 2) then
        this%nqmax = 5
    else if (Accuracy > 1) then
        this%nqmax = 4
    else
        this%nqmax = 3
    end if
    !note this may well be worse than the 5 optimized points

    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes and https://camb.info/maple/nu_integration_kernels.py
    if (allocated(this%nu_q)) deallocate(this%nu_q, this%nu_int_kernel)
    allocate(this%nu_q(this%nqmax))
    allocate(this%nu_int_kernel(this%nqmax))

    if (this%nqmax==3) then
        !Accurate at 2e-4 level
        this%nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
        this%nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
    else if (this%nqmax==4) then
        !Free-node least-squares fit for n=-4,-2..2 and v(am/q), 1/v(am/q)
        !Original rule kept here for reference:
        !this%nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
        !this%nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
        this%nu_q(1:4) = (/0.5802007037903776_dl, 2.2150938570691223_dl, 4.948032138986023_dl, 9.65253759848097_dl/)
        this%nu_int_kernel(1:4) = (/0.0082119845039711_dl, 1.1143258498419168_dl, &
            3.6819104154615907_dl, 0.8777790167504481_dl/)
    else if (this%nqmax==5) then
        !Exact for n=-4,-2..2 with remaining freedom fit to v(am/q), 1/v(am/q)
        !Original rule kept here for reference:
        !this%nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
        !this%nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)
        this%nu_q(1:5) = (/0.4620995950854295_dl, 1.7331898360630928_dl, 3.7956972681313816_dl, &
            7.2113928588584990_dl, 13.2665914595911080_dl/)
        this%nu_int_kernel(1:5) = (/0.0026946402277849193_dl, 0.46041071394952310_dl, 2.9207114780286405_dl, &
            2.1821643017186352_dl, 0.11621584305889110_dl/)
    else
        dq = (12 + this%nqmax/5)/real(this%nqmax)
        do i=1,this%nqmax
            q=(i-0.5d0)*dq
            this%nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            this%nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
    end if
    this%nu_int_kernel=this%nu_int_kernel/fermi_dirac_const

    end subroutine TNuPerturbations_init

    pure subroutine ThermalNuBackground_rho_P(this,am,rhonu,pnu)
    class(TThermalNuBackground), intent(in) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    real(dl) am2, am4, z, rfit, pfit
    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos. Use series solutions or
    !  direct smooth fits.

    if (am <= nu_rhop_am_min) then
        am2 = am**2
        am4 = am2**2
        rhonu = 1._dl + const2*am2 + am4*(nu_low_rho_c(0) + am2*nu_low_rho_c(1))
        pnu = (1._dl - const2*am2)/3._dl + am4*(nu_low_p_c(0) + am2*nu_low_p_c(1))
        return
    else if (am >= nu_rhop_am_max) then
        !Simple series solution (expanded in 1/(a*m))
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        pnu = 900._dl/120._dl/fermi_dirac_const*(zeta5-63._dl/4*Zeta7/am**2)/am
        return
    end if

    z = ((am - nu_rhop_fit_c)/(am + nu_rhop_fit_c))*nu_rhop_fit_inv_zmax

    rfit = (((((nu_fit_rho_c(5)*z + nu_fit_rho_c(4))*z + nu_fit_rho_c(3))*z &
        + nu_fit_rho_c(2))*z + nu_fit_rho_c(1))*z + nu_fit_rho_c(0))

    pfit = ((((((nu_fit_p_c(6)*z + nu_fit_p_c(5))*z + nu_fit_p_c(4))*z &
        + nu_fit_p_c(3))*z + nu_fit_p_c(2))*z + nu_fit_p_c(1))*z &
        + nu_fit_p_c(0))

    rhonu = (1._dl + nu_fit_rho_scale*am)*rfit
    pnu = pfit/(1._dl + nu_fit_p_denom_scale*am)

    end subroutine ThermalNuBackground_rho_P

    pure subroutine ThermalNuBackground_rho(this,am,rhonu)
    class(TThermalNuBackground), intent(in) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    real(dl) am2, am4, z, rfit

    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos. Use series solutions or
    !  direct smooth fits.

    if (am <= nu_rhop_am_min) then
        am2 = am**2
        am4 = am2**2
        rhonu = 1._dl + const2*am2 + am4*(nu_low_rho_c(0) + am2*nu_low_rho_c(1))
        return
    else if (am >= nu_rhop_am_max) then
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        return
    end if

    z = ((am - nu_rhop_fit_c)/(am + nu_rhop_fit_c))*nu_rhop_fit_inv_zmax

    rfit = (((((nu_fit_rho_c(5)*z + nu_fit_rho_c(4))*z + nu_fit_rho_c(3))*z &
        + nu_fit_rho_c(2))*z + nu_fit_rho_c(1))*z + nu_fit_rho_c(0))

    rhonu = (1._dl + nu_fit_rho_scale*am)*rfit

    end subroutine ThermalNuBackground_rho

    pure function rho_err(this, nu_mass)
    class(TThermalNuBackground), intent(in) :: this
    real(dl), intent(in) :: nu_mass
    real(dl) rho_err, rhonu

    call this%rho(nu_mass, rhonu)
    rho_err = rhonu - this%target_rho

    end function rho_err

    function ThermalNuBackground_find_nu_mass_for_rho(this,rho) result(nu_mass)
    !  Get eigenstate mass given input density (rho is neutrino density in units of one massless)
    !  nu_mass=m_n*c**2/(k_B*T_nu0).
    !  Get number density n of neutrinos from
    !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
    !  then m = Omega_nu/N_nu rho_crit /n if non-relativistic
    use MathUtils
    use config
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: rho
    real(dl) nu_mass, rhonu, rhonu1, delta
    real(dl) fzero
    integer iflag

    if (rho <= 1.001_dl) then
        !energy density all accounted for by massless result
        nu_mass=0
    else
        !Get mass assuming fully non-relativistic
        nu_mass=fermi_dirac_const/(1.5d0*zeta3)*rho

        if (nu_mass>4) then
            !  perturbative correction for velocity when nearly non-relativistic
            !  Error due to velocity < 1e-5 for mnu~0.06 but can easily correct (assuming non-relativistic today)
            !  Note that python does not propagate mnu to omnuh2 consistently to the same accuracy; but this makes
            !  fortran more internally consistent between input and computed Omega_nu h^2

            !Make perturbative correction for the tiny error due to the neutrino velocity
            call this%rho(nu_mass, rhonu)
            call this%rho(nu_mass*0.9, rhonu1)
            delta = rhonu - rho
            nu_mass = nu_mass*(1 + delta/((rhonu1 - rhonu)/0.1) )
        else
            !Directly solve to avoid issues with perturbative result when no longer very relativistic
            this%target_rho = rho
            call brentq(this,rho_err,0._dl,nu_mass,0.01_dl,nu_mass,fzero,iflag)
            if (iflag/=0) call GlobalError('find_nu_mass_for_rho failed to find neutrino mass')
        end if
    end if

    end function ThermalNuBackground_find_nu_mass_for_rho


    pure function ThermalNuBackground_drho(this,am,adotoa) result (rhonudot)
    !  Compute the time derivative of the mean density in massive neutrinos
    class(TThermalNuBackground), intent(in) :: this
    real(dl), intent(in) :: adotoa
    real(dl) rhonudot
    real(dl) am2, z, Dfit
    real(dl), intent(IN) :: am

    if (am <= am_min) then
        !rhonudot = 2*const2*am**2*adotoa
        am2 = am**2
        rhonudot = am2 * (2 * const2 + am2 * (.4399706676d-1 * log(am) &
            - .2970400378d-2 - .29331377855d-1 * am)) * adotoa
    else if (am >= am_max) then
        rhonudot = 3/(2*fermi_dirac_const)*(zeta3*am +( -(15*zeta5)/2 + 2835._dl/16*zeta7/am**2)/am)*adotoa
    else
        z = ((am - nu_fit_c)/(am + nu_fit_c))*nu_fit_inv_zmax

        Dfit = (((((nu_fit_D_c(5)*z + nu_fit_D_c(4))*z + nu_fit_D_c(3))*z &
            + nu_fit_D_c(2))*z + nu_fit_D_c(1))*z + nu_fit_D_c(0))

        rhonudot = am**2/(1._dl + 2._dl*am)*Dfit*adotoa
    end if

    end function ThermalNuBackground_drho

    end module MassiveNu
