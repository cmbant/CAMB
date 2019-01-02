    module DarkAge21cm
    use constants
    !Functions for collision rates.
    implicit none
    contains

    function kappa_HH_21cm(T, deriv)
    !Polynomial fit to Hydrogen-Hydrogen collision rate as function of Tmatter, from astro-ph/0608032
    !if deriv return d log kappa / d log T
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    !        real(dl), dimension(8), parameter :: fit = &
    !         (/ 0.00120402_dl, -0.0322247_dl,0.339581_dl, -1.75094_dl,4.3528_dl,-4.03562_dl, 1.26899_dl, -29.6113_dl /)
    integer, parameter :: n_table = 27
    integer, dimension(n_table), parameter :: Temps = &
        (/ 1, 2, 4, 6,8,10,15,20,25,30,40,50,60,70,80,90,100,200,300,500,700,1000,2000,3000,5000,7000,10000/)
    real, dimension(n_table), parameter :: rates = &
        (/ 1.38e-13, 1.43e-13,2.71e-13, 6.60e-13,1.47e-12,2.88e-12,9.10e-12,1.78e-11,2.73e-11,&
        3.67e-11,5.38e-11,6.86e-11,8.14e-11,9.25e-11, &
        1.02e-10,1.11e-10,1.19e-10,1.75e-10,2.09e-10,2.56e-10,2.91e-10,3.31e-10,4.27e-10,&
        4.97e-10,6.03e-10,6.87e-10,7.87e-10/)

    real(dl) kappa_HH_21cm, logT, logRate
    real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
    integer xlo, xhi
    real(dl) :: a0, b0, ho

    if (.not. allocated(logRates)) then

        allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
        logRates = log(real(rates,dl)*0.01**3)
        logTemps = log(real(Temps,dl))
        call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
    end if

    if (T<=Temps(1)) then
        if (deriv) then
            kappa_HH_21cm = 0
        else
            kappa_HH_21cm = rates(1)*0.01**3
        end if
        return
    elseif (T >=Temps(n_table)) then
        if (deriv) then
            kappa_HH_21cm = 0
        else
            kappa_HH_21cm = rates(n_table)*0.01**3
        end if
        return
    end if

    logT = log(T)
    xlo=0
    do xhi=2, n_table
        if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
        end  if
    end do
    xhi = xlo+1

    ho=logTemps(xhi)-logTemps(xlo)
    a0=(logTemps(xhi)-logT)/ho
    b0=1-a0

    if (deriv) then
        kappa_HH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
            ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
        !          kappa_HH_21cm = derivpolevl(logT,fit,7)
    else
        logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
        kappa_HH_21cm = exp(logRate)
        !          kappa_HH_21cm = exp(polevl(logT,fit,7))*0.01**3

    end if

    end function kappa_HH_21cm


    function kappa_eH_21cm(T, deriv)
    !Polynomail fit to electron-Hydrogen collision rate as function of Tmatter; from astro-ph/0608032
    !if deriv return d log kappa / d log T
    ! from astro-ph/0608032
    !    1 2.39e-10
    !    2 3.37e-10
    !    5 5.3e-10
    !    10 7.46e-10
    !    20 1.05e-9
    !    50 1.63e-9
    !    100 2.26e-9
    !    200 3.11e-9
    !    500 4.59e-9
    !    1000 5.92e-9
    !    2000 7.15e-9
    !    5000 8.17e-9
    !    10000 8.37e-9
    !    15000 8.29e-9
    !    20000 8.11e-9
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    real(dl), dimension(6), parameter :: fit = &
        (/5.86236d-005,  -0.00171375_dl, 0.0137303_dl, -0.0435277_dl, 0.540905_dl,-22.1596_dl /)
    real(dl) kappa_eH_21cm, logT

    logT = log(T)
    if (deriv) then
        kappa_eH_21cm = derivpolevl(logT,fit,5)
    else
        kappa_eH_21cm = exp(polevl(logT,fit,5))*0.01**3
    end if

    end function kappa_eH_21cm

    function kappa_pH_21cm(T, deriv) ! from astro-ph/0702487
    !Not actually used
    !Polynomail fit to proton-Hydrogen collision rate as function of Tmatter
    !if deriv return d log kappa / d log T
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    integer, parameter :: n_table = 17
    integer, dimension(n_table), parameter :: Temps = &
        (/ 1, 2, 5, 10,20,50,100,200,500,1000,2000,3000,5000,7000,10000,15000,20000/)
    real, dimension(n_table), parameter :: rates = &
        (/ 0.4028, 0.4517,0.4301,0.3699,0.3172,0.3047, 0.3379, 0.4043, 0.5471, 0.7051, 0.9167, 1.070, &
        1.301, 1.48,1.695,1.975,2.201/)

    real(dl) kappa_pH_21cm, logT, logRate
    real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
    integer xlo, xhi
    real(dl) :: a0, b0, ho
    real(dl):: factor = 0.01**3*1e-9

    if (.not. allocated(logRates)) then

        allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
        logRates = log(real(rates,dl)*factor)
        logTemps = log(real(Temps,dl))
        call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
    end if

    if (T<=Temps(1)) then
        if (deriv) then
            kappa_pH_21cm = 0
        else
            kappa_pH_21cm = rates(1)*factor
        end if
        return
    elseif (T >=Temps(n_table)) then
        if (deriv) then
            kappa_pH_21cm = 0
        else
            kappa_pH_21cm = rates(n_table)*factor
        end if
        return
    end if

    logT = log(T)
    xlo=0
    do xhi=2, n_table
        if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
        end  if
    end do
    xhi = xlo+1

    ho=logTemps(xhi)-logTemps(xlo)
    a0=(logTemps(xhi)-logT)/ho
    b0=1-a0

    if (deriv) then
        kappa_pH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
            ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
    else
        logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
        kappa_pH_21cm = exp(logRate)
    end if

    end function kappa_pH_21cm

    function polevl(x,coef,N)
    integer N
    real(dl) polevl
    real(dl) x,ans
    real(dl) coef(N+1)
    integer i

    ans=coef(1)
    do i=2,N+1
        ans=ans*x+coef(i)
    end do
    polevl=ans

    end function polevl


    function derivpolevl(x,coef,N)
    integer N
    real(dl) derivpolevl
    real(dl) x,ans
    real(dl) coef(N+1)
    integer i

    ans=coef(1)*N
    do i=2,N
        ans=ans*x+coef(i)*(N-i+1)
    end do
    derivpolevl=ans

    end function derivpolevl

    end module DarkAge21cm

