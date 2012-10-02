    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using CosmoRec 
    ! Author: Richard Shaw (CITA)
    !---------------------------------------------------------------------------------------------------
    ! 08.06.2012: added possibility to communicate Hubble (Jens Chluba)
    ! 12.06.2012: AL, changed interface to pass nnu directly; fixed spline extrapolation

    module Recombination
    use constants
    use AMLUtils
    implicit none
    private

    type RecombinationParams

        integer :: runmode
        real(dl) :: fdm  ! Dark matter annihilation efficiency

        ! Internal accuracy of CosmoRec (0 - normal, 3 - most accurate
        ! other values defined in CosmoRec.cpp source file)
        real(dl) :: accuracy 

    end type RecombinationParams

    character(LEN=*), parameter :: Recombination_Name = 'CosmoRec'

    logical :: first_run = .true.
    integer, parameter :: Nz = 10000
    real(dl) :: zmax = 1d4
    real(dl) :: zmin = 0._dl
    real(dl), dimension(Nz) :: zrec, arec, Hz
    real(dl), dimension(Nz) :: xrec, tmrec, x2rec, tm2rec

    public RecombinationParams, Recombination_xe, Recombination_tm, Recombination_init,   &
    Recombination_ReadParams, Recombination_SetDefParams, &
    Recombination_Validate, Recombination_Name


    contains

    subroutine Recombination_ReadParams(R, Ini)
    use IniFile
    Type(RecombinationParams) :: R
    Type(TIniFile) :: Ini

    R%runmode = Ini_Read_Int_File(Ini, 'cosmorec_runmode', 0)
    R%accuracy = Ini_Read_Double_File(Ini, 'cosmorec_accuracy', 0.0D0)
    R%fdm = Ini_Read_Double_File(Ini, 'cosmorec_fdm', 0.0D0)

    end subroutine Recombination_ReadParams


    subroutine Recombination_SetDefParams(R)
    type (RecombinationParams) ::R

    R%runmode = 0
    R%fdm = 0.0
    R%accuracy = 0

    end subroutine Recombination_SetDefParams


    subroutine Recombination_Validate(R, OK)
    Type(RecombinationParams), intent(in) :: R
    logical, intent(inout) :: OK

    if(R%runmode < 0 .or. R%runmode > 3) then
        write(*,*) "Invalid runmode for CosmoRec,"
        OK = .false.
    end if

    if(R%runmode < 2 .and. R%fdm > 1d-23) then
        write(*,*) "Dark matter annihilation rate too high. Will crash CosmoRec."
        OK = .false.
    end if

    if(R%accuracy < 0.0 .or. R%accuracy > 3.0) then
        write(*,*) "CosmoRec accuracy mode undefined."
        OK = .false.
    end if

    end subroutine Recombination_Validate



    function Recombination_tm(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_tm

    Recombination_tm = spline_val(a, arec, tmrec, tm2rec, Nz)

    end function Recombination_tm


    function Recombination_xe(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_xe

    Recombination_xe = spline_val(a, arec, xrec, x2rec, Nz)

    end function Recombination_xe




    subroutine Recombination_init(Recomb, OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu)
    !Would love to pass structure as arguments, but F90 would give circular reference...
    !hence mess passing parameters explcitly and non-generally

    use AMLUtils
    implicit none
    Type (RecombinationParams), intent(in) :: Recomb
    real(dl), intent(in) :: OmegaC, OmegaB, OmegaN, OmegaV, h0inp, tcmb, yp, num_nu

    real(dl) OmegaK
    integer :: i, label
    real(dl), dimension(5) :: runpars

    real(dl) dtauda
    external dtauda

    ! Calculate the curvature
    OmegaK = 1._dl - OmegaC - OmegaB - OmegaV
#ifdef MPI
    label = GetMpiRank()
#else
    label = 0
#endif
    ! Some feedback

    if (Feedback >1) then
        print *, "" ;
        print *, "==== CosmoRec parameters ====" ;

        print *, "Runmode: ", Recomb%runmode
        print "(a,f10.5)", " Omega_c: ", OmegaC
        print "(a,f10.5)", " Omega_b: ", OmegaB
        print "(a,f10.5)", " Omega_k: ", OmegaK
        print "(a,f10.5)", " Num_nu : ", num_nu
        print "(a,f10.5)", " Hubble : ", h0inp
        print "(a,f10.5)", " T_cmb  : ", tcmb
        print "(a,f10.5)", " Y_He   : ", yp
        print "(a,f10.5)", " f_dm   : ", Recomb%fdm
    end if

    ! Set runtime parameters
    runpars = 0._dl
    runpars(1) = Recomb%fdm ! Set dark matter annihilation efficiency
    runpars(2) = Recomb%accuracy

    ! Set redshifts to calculate at.
    do i=1,Nz
        zrec(i) = zmax - i*((zmax - zmin) / Nz)
        arec(i) = 1d0 / (1.0D0 + zrec(i))
        Hz(i) = 1/dtauda(1/(1._dl+zrec(i)))*(1._dl+zrec(i))**2/MPC_in_sec  
    end do

    ! internal Hubble function of CosmoRec is used
    !call CosmoRec_calc_cpp(Recomb%runmode, runpars, &
    !     OmegaC, OmegaB, OmegaK, num_nu, h0inp, tcmb, yp, &
    !     zrec, xrec, tmrec, Nz, label)

    ! version which uses camb Hubble function
    call CosmoRec_calc_h_cpp(Recomb%runmode, runpars, &
    OmegaC, OmegaB, OmegaK, num_nu, h0inp, tcmb, yp, &
    zrec, Hz, Nz, zrec, xrec, tmrec, Nz, label)

    call spline_double(arec, xrec, Nz, x2rec)
    call spline_double(arec, tmrec, Nz, tm2rec)

    ! print some output
    !open(unit=267,file="CosmoRec.out.Xe.II.dat")
    !do i=1,Nz
    !   write (267,*) zrec(i), xrec(i), tmrec(i)
    !end do
    !close(267)

    end subroutine Recombination_init



    ! General routine for cubic spline interpolation (see NR)
    real(dl) function spline_val(x, xv, yv, y2, n)

    real(dl), intent(in) :: x
    real(dl), intent(in) :: xv(n), yv(n), y2(n)
    integer, intent(in) :: n

    integer :: kh,kl,kn
    real(dl) :: h,a,b,c,d

    ! Extrapolate if value is above or below interval
    if(x < xv(1)) then
        spline_val = yv(1)
    else if(x > xv(n)) then
        spline_val = yv(n)
    else
        ! Bisection to find correct interval
        kh = n
        kl = 1
        do while(kh - kl > 1)
            kn = (kh + kl) / 2
            if(xv(kn) > x) then
                kh = kn
            else
                kl = kn
            end if
        end do

        ! Set up constants (a la NR)
        h = xv(kh) - xv(kl)

        a = (xv(kh) - x) / h
        b = (x - xv(kl)) / h
        c = (a**3 - a)* h**2 / 6
        d = (b**3 - b)* h**2 / 6

        spline_val = (a*yv(kl) + b*yv(kh) + c*y2(kl) + d*y2(kh))

    end if
    end function spline_val

    end module Recombination

