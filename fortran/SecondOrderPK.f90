
    module SecondOrderPK
    !Third-order Newtonian perturbation theory results for the non-linear correction
    !Only intended for use at very high redshift (z>10) where corrections are perturbative
    !May help to demonstrate how more general non-linear model classes could be implemented
    !See Appendix F of astro-ph/0702600 for equations and references
    !
    !The code is not intended for production use and may be quite slow
    use results
    use transfer
    use MathUtils
    implicit none
    private

    type, extends(TNonLinearModel) :: TSecondOrderPK
        type(MatterPowerData), pointer :: CAMB_Pk
        real(dl) :: min_store,min_k,max_k, k, PK, r, epsilon
        integer itf
        integer term
    contains
    procedure :: Init => TSecondOrderPK_Init
    procedure :: GetNonLinRatios => TSecondOrderPK_GetNonLinRatios
    procedure :: GetNonLinRatios_All => TSecondOrderPK_GetNonLinRatios_All
    procedure :: GetRatios
    procedure, private :: Integrand
    procedure, private :: Integrand_x
    procedure, private :: Integrand_Log
    procedure, private :: Integrand_series
    procedure, nopass :: SelfPointer => TSecondOrderPK_SelfPointer
    end type TSecondOrderPK

    integer, parameter :: term_dd=1,term_vv=2,term_dv=3

    public TSecondOrderPK
    contains

    subroutine TSecondOrderPK_Init(this, State)
    class(TSecondOrderPK) :: this
    class(TCAMBdata), target :: State

    this%Min_kh_nonlinear = 0.001

    end subroutine TSecondOrderPK_Init

    subroutine TSecondOrderPK_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TSecondOrderPK), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSecondOrderPK_SelfPointer

    function Integrand_x(this, x) result(int)
    class(TSecondOrderPK) :: this
    real(dl), intent(in) :: x
    real(dl) k, int

    associate(this_r => this%r, this_k => this%k)
        k = sqrt((1+this_r*(this_r-2*x)))*This_k
        if (k > this%min_store .and. k< this%max_k) then

            if (this%term==term_dd) then
                int = ((3*this_r + x*(7-10*this_r*x))/(1+this_r*(this_r-2*x)))**2
            elseif (this%term==term_vv) then
                int = ((this_r - x*(7-6*this_r*x))/(1+this_r*(this_r-2*x)))**2
            elseif (this%term==term_dv) then
                int = (this_r - x*(7-6*this_r*x))*(-3*this_r - x*(7-10*this_r*x))/(1+this_r*(this_r-2*x))**2

            end if
            int = int * MatterPowerData_k(this%CAMB_PK, k, this%itf)
        else
            int = 0
        end if
    end associate
    end function Integrand_x

    function Integrand_Log(this,p) result(int)
    !p = log r
    class(TSecondOrderPK) :: this
    real(dl), intent(in) ::p
    real(dl) r,r2,int, int22
    real(dl) :: xtol = 1.e-4_dl
    integer i

    r = exp(p)
    int = this%Integrand(r)*r
    end function Integrand_log

    function Integrand(this,r) result(int)
    class(TSecondOrderPK) :: this
    real(dl), intent(in) ::r
    real(dl) r2,int, int22
    real(dl) :: xtol = 1.e-4_dl
    integer i

    this%r = r
    r2=r**2
    if (this%term==term_dd) then

        Int = (12._dl/r2-158._dl+100._dl*r2-42._dl*r2**2)
        if (abs(r-1._dl) > 1e-6) then
            Int=Int  +3._dl/r2/r*(r2-1)**3*(7*r2 +2)*log((1._dl+r)/abs(1._dl-r))
        end if
        Int=Int*this%pk/252._dl

    elseif (this%term==term_vv) then

        Int = (12._dl/r2-82._dl+4._dl*r2-6._dl*r2**2)
        if (abs(r-1._dl) > 1e-6) then
            Int=Int  +3._dl/r2/r*(r2-1)**3*(r2 +2)*log((1._dl+r)/abs(1._dl-r))
        end if
        Int=Int*this%pk/84._dl

    elseif (this%term==term_dv) then

        Int = (24._dl/r2-202._dl+56._dl*r2-30._dl*r2**2)
        if (abs(r-1._dl) > 1e-6) then
            Int=Int  +3._dl/r2/r*(r2-1)**3*(5*r2 +4)*log((1._dl+r)/abs(1._dl-r))
        end if
        Int=Int*this%pk/252._dl

    end if

    if (r<this%epsilon) then
        int22=2*Integrate_Romberg(this,Integrand_x,-1._dl,1._dl, xtol)/98._dl
    else if (r >= 1-this%epsilon .and. r<= 1+this%epsilon) then
        int22=Integrate_Romberg(this,Integrand_x,-1._dl,(1._dl+r2-this%epsilon**2)/(2._dl*r), xtol)/98._dl
    else
        int22=Integrate_Romberg(this,Integrand_x,-1._dl,1._dl, xtol)/98._dl
    end if

    Int = Int+ int22
    Int=  Int * this%k**3 * MatterPowerData_k(this%CAMB_PK, r*this%k, this%itf)/(2._dl*const_pi)**2
    !put in k^3 here to keep answer sensible size
    end function Integrand


    function Integrand_series(this, p) result(int)
    class(TSecondOrderPK) :: this
    !For low r
    real(dl), intent(in) ::p
    real(dl) :: int, r
    integer i

    r = exp(p)
    Int=  r* r**2* this%pk* this%k**3 * MatterPowerData_k(this%CAMB_PK, r*this%k, this%itf)/(2._dl*const_pi)**2
    !put in k^3 here to keep answer sensible size
    !extra r because change to dr = r dp
    end function Integrand_series

    subroutine TSecondOrderPK_GetNonLinRatios(this, State, CAMB_Pk)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    class(TSecondOrderPK) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk

    call this%GetRatios(CAMB_PK, .false.)

    end subroutine TSecondOrderPK_GetNonLinRatios

    subroutine TSecondOrderPK_GetNonLinRatios_all(this,State, CAMB_Pk)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    class(TSecondOrderPK) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk

    call this%GetRatios(CAMB_PK, .true.)

    end subroutine TSecondOrderPK_GetNonLinRatios_all

    subroutine GetRatios(this,CAMB_Pk, DoVel)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    use splines
    class(TSecondOrderPK) :: this
    type(MatterPowerData), target :: CAMB_Pk
    logical, intent(in):: DoVel
    integer i, it,j
    real(dl) tmp,pnl, max_store,t1,t2,t3
    real(dl) :: rtol = 1e-5
    real(dl), parameter :: r_series = 0.003_dl
    real(dl), allocatable, dimension(:) :: dPdLogK
    real(dl) sc
    integer term

    this%CAMB_PK => CAMB_Pk

    CAMB_Pk%nonlin_ratio = 1

    if (doVel) then
        allocate(CAMB_Pk%nonlin_ratio_vv(CAMB_Pk%num_k,CAMB_Pk%num_z))
        allocate(CAMB_Pk%nonlin_ratio_vd(CAMB_Pk%num_k,CAMB_Pk%num_z))
        CAMB_Pk%nonlin_ratio_vv = 1
        CAMB_Pk%nonlin_ratio_vd = 1
    end if

    do term = term_dd, term_dv
        this%term = term

        max_store =  exp(CAMB_Pk%log_kh(CAMB_PK%num_k))
        this%min_store =  exp(CAMB_Pk%log_kh(1))

        this%min_k = this%min_store
        this%max_k = max_store

        do it = 1, CAMB_Pk%num_z
            this%itf = it

            allocate(dPdLogK(CAMB_PK%num_k))
            !Get first derivative needed for series expansion at low r
            call spline_deriv(CAMB_Pk%log_kh, CAMB_Pk%matpower(:,it), CAMB_Pk%ddmat(:,it), dPdLogK,CAMB_PK%num_k)

            do i=1, CAMB_PK%num_k

                this%k = exp(CAMB_Pk%log_kh(i))

                if (this%k > this%Min_kh_nonlinear) then

                    this%min_k = max(this%min_store,0.4 * this%k)
                    this%max_k = min(max_store,300* this%k)

                    this%epsilon = this%min_k/this%k
                    if (this%epsilon >=0.5) stop 'epsilon >=0.5'


                    this%pk =MatterPowerData_k(CAMB_PK, this%k, this%itf)

                    if (this%min_k > this%min_store) then

                        if (this%min_store/this%k < r_series) then
                            !Series result
                            if (term==term_vv) then
                                pnl = 94./245*Integrate_romberg(this,Integrand_series,log(this%min_store/this%k), &
                                    log(r_series),rtol*this%pk,abs_tol=.true.)
                                pnl = pnl*(1 - 217._dl/141._dl*dPdLogK(i) &
                                    + 49._dl/94._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                            else if (term==term_dv) then
                                pnl = 2558./2205*Integrate_romberg(this,Integrand_series, &
                                    log(this%min_store/this%k),log(r_series),rtol*this%pk,abs_tol=.true.)
                                pnl = pnl*(1 - 819._dl/1279._dl*dPdLogK(i)  &
                                    + 441._dl/2558._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                            else if (term==term_dd) then
                                pnl = 5038./2205*Integrate_romberg(this,Integrand_series,&
                                    log(this%min_store/this%k),log(r_series),rtol*this%pk,abs_tol=.true.)
                                pnl = pnl*(1 - 987._dl/2519._dl*dPdLogK(i) &
                                    + 441._dl/5038._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                            end if
                            !plus integral with log spacing
                            pnl = pnl+Integrate_romberg(this,Integrand_Log,log(r_series), &
                                log(this%epsilon),rtol*this%pk,20,abs_tol=.true.)
                        else
                            pnl = Integrate_romberg(this,Integrand_Log,log(this%min_store/this%k), &
                                log(this%epsilon),rtol*this%pk,20,abs_tol=.true.)
                        end if

                    else
                        pnl = 0
                    end if


                    t1 = Integrate_romberg(this,Integrand,this%epsilon,1-this%epsilon,rtol*this%pk,abs_tol=.true.)

                    if (1+this%epsilon*2<this%max_k/this%k) then
                        t2 = Integrate_romberg(this,Integrand,1+this%epsilon,1+this%epsilon*2,&
                            rtol*this%pk,abs_tol=.true.) !function falls quite rapidly
                        t2 =t2+ Integrate_romberg(this,Integrand,1+this%epsilon*2, &
                            this%max_k/this%k,rtol*this%pk,abs_tol=.true.)
                    else
                        t2 = Integrate_romberg(this,Integrand,1+this%epsilon,1+this%epsilon*2,rtol*this%pk,abs_tol=.true.)
                    end if
                    t3 = Integrate_romberg(this,Integrand,1-this%epsilon,1+this%epsilon,rtol*this%pk,abs_tol=.true.)
                    ! sc = this_K**3/(2*pi**2)
                    ! write (*,'(5e15.5)') this_k,this_PK*sc,pnl*sc, (this_PK+pnl)*sc, pnl/this_PK

                    pnl = pnl+t1+t2+t3

                    ! sc = this_K**3/(2*pi**2)
                    ! write (*,'(5e15.5)') this_k,this_PK*sc,pnl*sc, (this_PK+pnl)*sc, pnl/this_PK


                    if (term==term_dd) then
                        CAMB_Pk%nonlin_ratio(i,this%itf) = sqrt(abs(1+pnl/this%pk))
                    elseif (term==term_vv) then
                        CAMB_Pk%nonlin_ratio_vv(i,this%itf) = sqrt(abs(1+pnl/this%pk))
                    else
                        CAMB_Pk%nonlin_ratio_vd(i,this%itf) = sqrt(abs(1+pnl/this%pk))
                    end if
                end if

            enddo

            deallocate(dPdLogk)

        end do !redshifts

        if (.not. doVel) exit

    end do !terms

    end subroutine GetRatios

    end module SecondOrderPK

