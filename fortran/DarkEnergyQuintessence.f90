    ! Equations module allowing for fairly general quintessence models
    !
    ! by Antony Lewis (http://cosmologist.info/)

    !!FIX March 2005: corrected update to next treatment of tight coupling
    !!Fix Oct 2011: a^2 factor in ayprime(EV%w_ix+1) [thanks Raphael Flauger]
    ! Oct 2013: update for latest CAMB, thanks Nelson Lima, Martina Schwind
    ! May 2020: updated for CAMB 1.x+

    ! Notes at http://antonylewis.com/notes/CAMB.pdf

    !This module is not well tested, use at your own risk!

    !Need to specify Vofphi function, and also initial_phi
    !You may also need to change other things to get it to work with different types of quintessence model

    !It works backwards, in that it assumes Omega_de is Omega_Q today, then does a binary search on the
    !initial conditions to find what is required to give that Omega_Q today after evolution.

    module Quintessence
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    private

    real(dl), parameter :: Tpl= sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time


    !This is an example, it's not supposed to be a rigorous model!  (not very well tested)
    type, extends(TDarkEnergyModel) :: TQuintessence
        real(dl) :: n = 3._dl
        real(dl) :: f =0.05 ! sqrt(8*pi*G)*f
        real(dl) :: m = 5d-54 !m in reduced Planck mass units
        real(dl) :: theta_i = 3.1_dl !initial value of phi/f
        real(dl) :: frac_lambda0 = 1._dl !fraction of dark energy density that is cosmological constant today
        integer :: npoints = 5000 !baseline number of log a steps; will be increased if needed when there are oscillations
        integer :: min_steps_per_osc = 10
        integer :: npoints_linear, npoints_log
        real(dl) :: dloga, da, astart, log_astart, max_a_log
        real(dl), dimension(:), allocatable :: aVals, phi_a, phidot_a, fde, ddphi_a, ddphidot_a
        class(CAMBdata), pointer, private :: State
    contains
    procedure :: ReadParams =>  TQuintessence_ReadParams
    procedure, nopass :: PythonClass => TQuintessence_PythonClass
    procedure, nopass :: SelfPointer => TQuintessence_SelfPointer
    procedure :: Init => TQuintessence_Init
    procedure :: PerturbedStressEnergy => TQuintessence_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TQuintessence_PerturbationEvolve
    procedure :: BackgroundDensityAndPressure => TQuintessence_BackgroundDensityAndPressure
    procedure :: EvolveBackground
    procedure :: EvolveBackgroundLog
    procedure :: Vofphi
    !    procedure, private :: GetOmegaFromInitial
    procedure, private :: ValsAta
    procedure, private :: phidot_start => TQuintessence_phidot_start
    end type TQuintessence

    procedure(TClassDverk) :: dverk

    public TQuintessence
    contains


    function TQuintessence_PythonClass()
    character(LEN=:), allocatable :: TQuintessence_PythonClass

    TQuintessence_PythonClass = 'Quintessence'

    end function TQuintessence_PythonClass

    subroutine TQuintessence_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TQuintessence), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TQuintessence_SelfPointer

    subroutine TQuintessence_ReadParams(this, Ini)
    use IniObjects
    class(TQuintessence) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyModel%ReadParams(Ini)

    end subroutine TQuintessence_ReadParams


    function Vofphi(this, phi, deriv)
    !The input variable phi is sqrt(8*Pi*G)*psi
    !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
    !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
    class(TQuintessence) :: this
    real(dl) phi,Vofphi
    integer deriv
    real(dl) theta, costheta

    ! Assume f = sqrt(kappa*f_theory) = f_theory/M_pl
    ! m = m_theory/M_Pl
    theta = phi/this%f
    if (deriv==0) then
        Vofphi = this%m**2*this%f**2*(1 - cos(theta))**this%n
    else if (deriv ==1) then
        Vofphi = this%m**2*this%f*this%n*(1 - cos(theta))**(this%n-1)*sin(theta)
    else if (deriv ==2) then
        costheta = cos(theta)
        Vofphi = this%m**2*this%n*(1 - costheta)**(this%n-1)*(this%n*(1+costheta) -1)
    end if
    Vofphi = Vofphi * MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2


    !!dimensionless normalization kappa*V*T_pl^2
    !!         norm= 1
    !norm= 1d-122
    !
    !norm = norm * (Mpc/c)**2 /Tpl**2 !convert to units of 1/Mpc^2
    !
    !if (deriv==0) then
    !    Vofphi= norm*this%m*exp(-this%sigma_model*phi)
    !else if (deriv ==1) then
    !    Vofphi=-norm*this%m*sigma_model*exp(-this%sigma_model*phi)
    !else if (deriv ==2) then
    !    Vofphi=norm*this%m*sigma_model**2*exp(-this%sigma_model*phi)
    !else
    !    stop 'Invalid deriv in Vofphi'
    !end if

    end function Vofphi


    subroutine TQuintessence_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TQuintessence), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) V, a2, grhov_lambda, phi, phidot

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    elseif (a >= this%astart) then
        a2 = a**2
        call this%ValsAta(a,phi,phidot)
        V = this%Vofphi(phi,0)
        grhov_t = phidot**2/2 + a2*V
        if (this%frac_lambda0>0) then
            grhov_lambda = this%frac_lambda0*grhov*a2
            grhov_t = grhov_t +  grhov_lambda
            if (present(w)) then
                w = ((phidot**2/2 - a2*V) - grhov_lambda)/grhov_t
            end if
        elseif (present(w)) then
            w = (phidot**2/2 - a2*V)/grhov_t
        end if
    else
        grhov_t=0
        if (present(w)) w = -1
    end if

    end subroutine TQuintessence_BackgroundDensityAndPressure

    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
    ! Evolve the background equation in terms of loga.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TQuintessence) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) loga, a, a2, tot
    real(dl) phi, tmp, phidot, grhode, adotoa
    integer nu_i

    a = exp(loga)
    a2=a**2
    phi = y(1)
    phidot = y(2)/a2

    grhode=a2*(0.5d0*phidot**2 + a2*(this%Vofphi(phi,0) + this%frac_lambda0*this%State%grhov))
    tot = this%state%grho_no_de(a) + grhode

    adotoa=sqrt(tot/3.0d0)/a
    yprime(1)=phidot/adotoa !d phi /d ln a
    yprime(2)= -a2**2*this%Vofphi(phi,1)/adotoa

    end subroutine EvolveBackgroundLog

    subroutine EvolveBackground(this,num,a,y,yprime)
    ! Evolve the background equation in terms of a.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TQuintessence) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) a, a2, tot
    real(dl) phi, grhode, phidot, rhonu, adot
    integer nu_i

    a2=a**2
    phi = y(1)
    phidot = y(2)/a2

    grhode=a2*(0.5d0*phidot**2 + a2*(this%Vofphi(phi,0) + this%frac_lambda0*this%State%grhov))
    tot = this%state%grho_no_de(a) + grhode

    adot=sqrt(tot/3.0d0)
    yprime(1)=phidot/adot !d phi /d a
    yprime(2)= -a2**2*this%Vofphi(phi,1)/adot

    end subroutine EvolveBackground

    !real(dl) function GetOmegaFromInitial(this, astart,phi,phidot,atol)
    !!Get omega_de today given particular conditions phi and phidot at a=astart
    !class(TQuintessence) :: this
    !real(dl), intent(IN) :: astart, phi,phidot, atol
    !integer, parameter ::  NumEqs=2
    !real(dl) c(24),w(NumEqs,9), y(NumEqs), ast
    !integer ind, i
    !
    !ast=astart
    !ind=1
    !y(1)=phi
    !y(2)=phidot*astart**2
    !call dverk(this,NumEqs,EvolveBackground,ast,y,1._dl,atol,ind,c,NumEqs,w)
    !call EvolveBackground(this,NumEqs,1._dl,y,w(:,1))
    !
    !GetOmegaFromInitial=(0.5d0*y(2)**2 + Vofphi(y(1),0))/this%State%grhocrit !(3*adot**2)
    !
    !end function GetOmegaFromInitial


    real(dl) function TQuintessence_phidot_start(this,phi)
    class(TQuintessence) :: this
    real(dl) :: phi

    TQuintessence_phidot_start = 0

    end function TQuintessence_phidot_start

    subroutine TQuintessence_Init(this, State)
    use FileUtils
    class(TQuintessence), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) aend, afrom
    integer, parameter ::  NumEqs=2
    real(dl) c(24),w(NumEqs,9), y(NumEqs),atol
    integer ind, i, ix
    real(dl), parameter :: splZero = 0._dl
    real(dl) rhofrac, lastsign, da_osc, last_a
    real(dl) initial_phi, initial_phidot, a2
    Type(TTextFile) Fout
    real(dl), dimension(:), allocatable :: aVals, phi_a, phidot_a, fde
    integer npoints, tot_points

    !Make interpolation table, etc,
    !At this point massive neutrinos have been initialized
    !so nu1 can be used to get their density and pressure at scale factor a
    !Other built-in components have density and pressure scaling trivially with a

    select type(State)
    class is (CAMBdata)
        this%State => State
    end select

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2

    initial_phi = this%theta_i*this%f
    this%astart=1d-7
    this%log_astart = log(this%astart)

    this%dloga = (-this%log_astart)/(this%npoints-1)

    !use log spacing in a up to max_a_log, then linear. Switch where step matches
    this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)
    npoints = (log(this%max_a_log)-this%log_astart)/this%dloga + 1

    if (allocated(this%phi_a)) then
        deallocate(this%phi_a,this%phidot_a)
        deallocate(this%ddphi_a,this%ddphidot_a, this%aVals)
    end if
    allocate(phi_a(npoints),phidot_a(npoints), aVals(npoints), fde(npoints))

    !initial_phi  = 10  !  0.3*grhom/m**3
    !initial_phi2 = 100 !   6*grhom/m**3
    !
    !!           initial_phi  = 65 !  0.3*grhom/m**3
    !!           initial_phi2 = 65 !   6*grhom/m**3
    !
    !astart=1d-9
    !
    !!See if initial conditions are giving correct omega_de now
    !atol=1d-8
    !initial_phidot =  astart*this%phidot_start(this%initial_phi)
    !om1= this%GetOmegaFromInitial(astart,initial_phi,initial_phidot, atol)
    !
    !print*, State%omega_de, 'first trial:', om1
    !if (abs(om1-State%omega_de > this%omega_tol)) then
    !    !if not, do binary search in the interval
    !    OK=.false.
    !    initial_phidot = astart*this%phidot_start(initial_phi2)
    !    om2= this%GetOmegaFromInitial(astart,initial_phi2,initial_phidot, atol)
    !    if (om1 > State%omega_de .or. om2 < State%omega_de) then
    !        write (*,*) 'initial phi values must bracket required value.  '
    !        write (*,*) 'om1, om2 = ', real(om1), real(om2)
    !        stop
    !    end if
    !    do iter=1,100
    !        deltaphi=initial_phi2-initial_phi
    !        phi =initial_phi + deltaphi/2
    !        initial_phidot =  astart*Quint_phidot_start(phi)
    !        om = this%GetOmegaFromInitial(astart,phi,initial_phidot,atol)
    !        if (om < State%omega_de) then
    !            om1=om
    !            initial_phi=phi
    !        else
    !            om2=om
    !            initial_phi2=phi
    !        end if
    !        if (om2-om1 < 1d-3) then
    !            OK=.true.
    !            initial_phi = (initial_phi2+initial_phi)/2
    !            if (FeedbackLevel > 0) write(*,*) 'phi_initial = ',initial_phi
    !            exit
    !        end if
    !
    !    end do !iterations
    !    if (.not. OK) stop 'Search for good intial conditions did not converge' !this shouldn't happen
    !
    !end if !Find initial

    atol=1d-6
    y(1)=initial_phi
    initial_phidot =  this%astart*this%phidot_start(initial_phi)
    y(2)= initial_phidot*this%astart**2

    phi_a(1)=y(1)
    phidot_a(1)=y(2)/this%astart**2
    aVals(1)=this%astart
    da_osc = 1
    last_a = this%astart

    ind=1
    afrom=this%log_astart
    call Fout%CreateOpenFile('z:\test.txt')
    do i=1, npoints-1
        aend = this%log_astart + this%dloga*i
        ix = i+1
        aVals(ix)=exp(aend)
        a2 = aVals(ix)**2
        call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,atol,ind,c,NumEqs,w)
        call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
        phi_a(ix)=y(1)
        phidot_a(ix)=y(2)/a2
        if (i==1) then
            lastsign = y(2)
        elseif (y(2)*lastsign < 0) then
            !derivative has changed sign. Use to probe any oscillation scale:
            da_osc = min(da_osc, exp(aend) - last_a)
            last_a = exp(aend)
            lastsign= y(2)
        end if

        !Define fde is ratio of energy density to total assuming the neutrinos fully relativistic
        !adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates)))/3)

        fde(ix) = 1/((this%state%grho_no_de(aVals(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
            /(a2*(0.5d0* phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
        call Fout%Write(exp(aend), phi_a(ix), phidot_a(ix), fde(ix))
        if (aVals(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) then
            !Step size getting too big to sample oscillations well
            exit
        end if
    end do

    ! Do remaining steps with linear spacing in a, trying to be small enough
    this%npoints_log = ix
    this%max_a_log = aVals(ix)
    this%da = min(this%max_a_log *(exp(this%dloga)-1), &
        da_osc/this%min_steps_per_osc, (1- this%max_a_log)/(this%npoints-this%npoints_log))
    this%npoints_linear = int((1- this%max_a_log)/ this%da)+1
    this%da = (1- this%max_a_log)/this%npoints_linear

    tot_points = this%npoints_log+this%npoints_linear
    allocate(this%phi_a(tot_points),this%phidot_a(tot_points))
    allocate(this%ddphi_a(tot_points),this%ddphidot_a(tot_points))
    allocate(this%aVals(tot_points), this%fde(tot_points))
    this%aVals(1:ix) = aVals(1:ix)
    this%phi_a(1:ix) = phi_a(1:ix)
    this%phidot_a(1:ix) = phidot_a(1:ix)
    this%aVals(1:ix) = aVals(1:ix)
    this%fde(1:ix) = fde(1:ix)

    ind=1
    afrom = this%max_a_log
    do i=1, this%npoints_linear
        ix = this%npoints_log + i
        aend = this%max_a_log + this%da*i
        a2 =aend**2
        this%aVals(ix)=aend
        call dverk(this,NumEqs,EvolveBackground,afrom,y,aend,atol,ind,c,NumEqs,w)
        call EvolveBackground(this,NumEqs,aend,y,w(:,1))
        this%phi_a(ix)=y(1)
        this%phidot_a(ix)=y(2)/a2

        this%fde(ix) = 1/((this%state%grho_no_de(aend) +  this%frac_lambda0*this%State%grhov*a2**2) &
            /(a2*(0.5d0* this%phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
        call Fout%Write(aend, this%phi_a(ix), this%phidot_a(ix), this%fde(ix))
    end do

    call Fout%Close()

    call spline(this%aVals,this%phi_a,tot_points,splZero,splZero,this%ddphi_a)
    call spline(this%aVals,this%phidot_a,tot_points,splZero,splZero,this%ddphidot_a)

    end subroutine TQuintessence_Init

    subroutine ValsAta(this,a,aphi,aphidot)
    class(TQuintessence) :: this
    !Do interpolation for background phi and phidot at a
    real(dl) a, aphi, aphidot
    real(dl) a0,b0,ho2o6,delta,da
    integer ix

    if (a >= 0.9999999d0) then
        aphi= this%phi_a(this%npoints_linear+this%npoints_log)
        aphidot= this%phidot_a(this%npoints_linear+this%npoints_log)
        return
    elseif (a > this%max_a_log) then
        delta= a-this%max_a_log
        ix = this%npoints_log + int(delta/this%da)
    else
        delta= log(a)-this%log_astart
        ix = int(delta/this%dloga)+1
    end if
    da = this%aVals(ix+1) - this%aVals(ix)
    a0 = (this%aVals(ix+1) - a)/da
    b0 = 1 - a0
    ho2o6 = da**2/6._dl
    aphi=b0*this%phi_a(ix+1) + a0*(this%phi_a(ix)-b0*((a0+1)*this%ddphi_a(ix)+(2-a0)*this%ddphi_a(ix+1))*ho2o6)
    aphidot=b0*this%phidot_a(ix+1) + a0*(this%phidot_a(ix)-b0*((a0+1)*this%ddphidot_a(ix)+(2-a0)*this%ddphidot_a(ix+1))*ho2o6)

    end subroutine ValsAta


    subroutine TQuintessence_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TQuintessence), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) phi, phidot, clxq, vq

    call this%ValsAta(a,phi,phidot)
    clxq=ay(w_ix)
    vq=ay(w_ix+1)
    dgrhoe= phidot*vq +clxq*a**2*this%Vofphi(phi,1)
    dgqe= k*phidot*clxq

    end subroutine TQuintessence_PerturbedStressEnergy


    subroutine TQuintessence_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TQuintessence), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) clxq, vq, phi, phidot

    call this%ValsAta(a,phi,phidot) !wasting time calling this again..
    clxq=y(w_ix)
    vq=y(w_ix+1)
    ayprime(w_ix)= vq
    ayprime(w_ix+1) = - 2*adotoa*vq - k*z*phidot - k**2*clxq - a**2*clxq*this%Vofphi(phi,2)

    end subroutine TQuintessence_PerturbationEvolve

    end module Quintessence
