
      module NonLinear
      use ModelParams
      use transfer
      implicit none 
      private
       
       real, parameter :: Min_kh_nonlinear = 0.001
       real(dl) :: min_store,min_k,max_k, this_k, this_PK, this_r, epsilon
       integer itf
       integer, parameter :: term_dd=1,term_vv=2,term_dv=3
       integer term

       type(MatterPowerData) , pointer :: CAMB_Pk

      public Min_kh_nonlinear,NonLinear_GetNonLinRatios
      
      contains

      function Integrand_x(x) result(int)
       real(dl), intent(in) :: x
       real(dl) k, int

        k = sqrt((1+this_r*(this_r-2*x)))*This_k
        if (k > min_store .and. k< max_k) then
  
         if (term==term_dd) then
          int = ((3*this_r + x*(7-10*this_r*x))/(1+this_r*(this_r-2*x)))**2
         elseif (term==term_vv) then
          int = ((this_r - x*(7-6*this_r*x))/(1+this_r*(this_r-2*x)))**2
         elseif (term==term_dv) then
          int = (this_r - x*(7-6*this_r*x))*(-3*this_r - x*(7-10*this_r*x))/(1+this_r*(this_r-2*x))**2
         
         end if 
          int = int * MatterPowerData_k(CAMB_PK, k, itf)  
        else
         int = 0
        end if
      end function Integrand_x


      
      function Integrand_Log(p) result(int)
     !p = log r
       real(dl), intent(in) ::p
       real(dl) r,r2,int, int22
       real(dl) :: xtol = 1.e-4_dl
       integer i

        r = exp(p)
        int = Integrand(r)*r
      end function Integrand_log

      function Integrand(r) result(int)
       real(dl), intent(in) ::r
       real(dl) r2,int, int22
       real(dl) :: xtol = 1.e-4_dl
       integer i


        this_r = r    
        r2=r**2
        if (term==term_dd) then

         Int = (12._dl/r2-158._dl+100._dl*r2-42._dl*r2**2)
         if (abs(r-1._dl) > 1e-6) then
          Int=Int  +3._dl/r2/r*(r2-1)**3*(7*r2 +2)*log((1._dl+r)/abs(1._dl-r))        
         end if
         Int=Int*This_PK/252._dl

        elseif (term==term_vv) then

         Int = (12._dl/r2-82._dl+4._dl*r2-6._dl*r2**2)
         if (abs(r-1._dl) > 1e-6) then
          Int=Int  +3._dl/r2/r*(r2-1)**3*(r2 +2)*log((1._dl+r)/abs(1._dl-r))        
         end if
         Int=Int*This_PK/84._dl

        elseif (term==term_dv) then

         Int = (24._dl/r2-202._dl+56._dl*r2-30._dl*r2**2)
         if (abs(r-1._dl) > 1e-6) then
          Int=Int  +3._dl/r2/r*(r2-1)**3*(5*r2 +4)*log((1._dl+r)/abs(1._dl-r))        
         end if
         Int=Int*This_PK/252._dl

        end if

         if (r<epsilon) then
          int22=2*rombint2(Integrand_x,-1._dl,1._dl, xtol)/98._dl
         else if (r >= 1-epsilon .and. r<= 1+epsilon) then 
          int22=rombint2(Integrand_x,-1._dl,(1._dl+r2-epsilon**2)/(2._dl*r), xtol)/98._dl
         else
          int22=rombint2(Integrand_x,-1._dl,1._dl, xtol)/98._dl
         end if

        Int = Int+ int22 
        Int=  Int * this_k**3 * MatterPowerData_k(CAMB_PK, r*this_k, itf)/(2._dl*pi)**2 
           !put in k^3 here to keep answer sensible size      
      end function Integrand


     function Integrand_series(p) result(int)
      !For low r
       real(dl), intent(in) ::p
       real(dl) :: int, r
       integer i

        r = exp(p)
        Int=  r* r**2* This_PK* this_k**3 * MatterPowerData_k(CAMB_PK, r*this_k, itf)/(2._dl*pi)**2 
           !put in k^3 here to keep answer sensible size      
           !extra r because change to dr = r dp
      end function Integrand_series


     subroutine NonLinear_GetNonLinRatios(inCAMB_Pk, velocities)
     !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
     !for each redshift and wavenumber
     !This implementation uses Halofit
      type(MatterPowerData), target :: inCAMB_Pk
      integer i, it,j
      logical, intent(in), optional :: velocities
      logical doVel
      real(dl) tmp,pnl, max_store,t1,t2,t3
      real(dl) :: rtol = 1e-5
      real(dl), parameter :: r_series = 0.003_dl
      real(dl), allocatable, dimension(:) :: dPdLogK
      real(dl) sc

       CAMB_PK => inCAMB_Pk

       term = term_dd

       CAMB_Pk%nonlin_ratio = 1
       doVel = .false.
       if (present(velocities)) then
        doVel = velocities
       end if

       if (doVel) then
         allocate(CAMB_Pk%nonlin_ratio_vv(CAMB_Pk%num_k,CAMB_Pk%num_z))  
         allocate(CAMB_Pk%nonlin_ratio_vd(CAMB_Pk%num_k,CAMB_Pk%num_z))  
         CAMB_Pk%nonlin_ratio_vv = 1
         CAMB_Pk%nonlin_ratio_vd = 1

       end if   


       do term = term_dd, term_dv 

       max_store =  exp(CAMB_Pk%log_kh(CAMB_PK%num_k))
       min_store =  exp(CAMB_Pk%log_kh(1))

       min_k = min_store
       max_k = max_store

       do it = 1, CAMB_Pk%num_z
          itf = it

          allocate(dPdLogK(CAMB_PK%num_k))
         !Get first derivative needed for series expansion at low r
          call spline_deriv(CAMB_Pk%log_kh, CAMB_Pk%matpower(1,it), CAMB_Pk%ddmat(1,it), dPdLogK,CAMB_PK%num_k)

          do i=1, CAMB_PK%num_k
             
             this_k = exp(CAMB_Pk%log_kh(i))
         

             if (this_k > Min_kh_nonlinear) then

                 min_k = max(min_store,0.4 * this_K)
                 max_k = min(max_store,300* this_K)
             
                 epsilon = min_k/This_k
                 if (epsilon >=0.5) stop 'epsilon >=0.5'

 
                 this_PK =MatterPowerData_k(CAMB_PK, this_k, itf) 

                 if (min_k > min_store) then
  
                   if (min_store/this_k < r_series) then
                    !Series result
                    if (term==term_vv) then
                     pnl = 94./245*Rombint_abs(Integrand_series,log(min_store/this_k),log(r_series),rtol*this_PK)
                     pnl = pnl*(1 - 217._dl/141._dl*dPdLogK(i) + 49._dl/94._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                    else if (term==term_dv) then
                     pnl = 2558./2205*Rombint_abs(Integrand_series,log(min_store/this_k),log(r_series),rtol*this_PK)
                     pnl = pnl*(1 - 819._dl/1279._dl*dPdLogK(i) + 441._dl/2558._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                    else if (term==term_dd) then
                     pnl = 5038./2205*Rombint_abs(Integrand_series,log(min_store/this_k),log(r_series),rtol*this_PK)
                     pnl = pnl*(1 - 987._dl/2519._dl*dPdLogK(i) + 441._dl/5038._dl*( CAMB_Pk%ddmat(i,it) + dPdLogK(i)**2))
                    end if
                    !plus integral with log spacing
                    pnl = pnl+Rombint_abs(Integrand_Log,log(r_series),log(epsilon),rtol*this_PK,20)
                   else
                    pnl = Rombint_abs(Integrand_Log,log(min_store/this_k),log(epsilon),rtol*this_PK,20)
                   end if 

                 else 
                   pnl = 0
                 end if
                 

               t1 = Rombint_abs(Integrand,epsilon,1-epsilon,rtol*this_PK)
     ! t1 = Rombint_abs(Integrand_Log,log(epsilon),log(1-epsilon),rtol*this_PK)

               if (1+epsilon*2<max_k/this_k) then
                t2 = Rombint_abs(Integrand,1+epsilon,1+epsilon*2,rtol*this_PK) !function falls quite rapidly
                t2 =t2+ Rombint_abs(Integrand,1+epsilon*2,max_k/this_k,rtol*this_PK)
               else
                t2 = Rombint_abs(Integrand,1+epsilon,1+epsilon*2,rtol*this_PK)
               end if
               t3 = Rombint_abs(Integrand,1-epsilon,1+epsilon,rtol*this_PK)
              ! sc = this_K**3/(2*pi**2)
              ! write (*,'(5e15.5)') this_k,this_PK*sc,pnl*sc, (this_PK+pnl)*sc, pnl/this_PK

                pnl = pnl+t1+t2+t3
              
              
              ! sc = this_K**3/(2*pi**2)
              ! write (*,'(5e15.5)') this_k,this_PK*sc,pnl*sc, (this_PK+pnl)*sc, pnl/this_PK


               if (term==term_dd) then
                CAMB_Pk%nonlin_ratio(i,itf) = sqrt(abs(pnl/this_PK))
               elseif (term==term_vv) then
                CAMB_Pk%nonlin_ratio_vv(i,itf) = sqrt(abs(pnl/this_PK))
               else
                CAMB_Pk%nonlin_ratio_vd(i,itf) = sqrt(abs(pnl/this_PK))
               end if 
             end if

          enddo

          deallocate(dPdLogk)

      end do !redshifts

      if (.not. doVel) exit

      end do !terms
            
      end subroutine NonLinear_GetNonLinRatios
       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        recursive function rombint2(f,a,b,tol, maxit) 
!Use absolute error
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, intent(in), optional :: maxit
        integer:: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint2
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!
        if (present(maxit)) then
          MAXITER=maxit
        end if

        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
        do
          i=i+1
          if (i > MAXITER.or.(i > 5.and.abs(error) < tol)) exit
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
          do k=1,nint
            g0=g0+f(a+(k+k-1)*h)
          end do
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
          do j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
          end do
          if (abs(g0) > tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        end do
        rombint2=g0
        if (i > MAXITER.and.abs(error) > tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint2,error, tol
        end if
        
        end function rombint2

     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        recursive function rombint_abs(f,a,b,tol, maxit) 
!Use absolute error
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, intent(in), optional :: maxit
        integer:: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint_abs
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!
        if (present(maxit)) then
          MAXITER=maxit
        end if

        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
        do
          i=i+1
          if (i > MAXITER.or.(i > 5.and.abs(error) < tol)) exit
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
          do k=1,nint
            g0=g0+f(a+(k+k-1)*h)
          end do
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
          do j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
          end do
            error=abs(gmax-g0)
          gmax=g0
          g(jmax+1)=g0
        end do
        rombint_abs=g0
        if (i > MAXITER.and.abs(error) > tol)  then
          write(*,*) 'Warning: Rombint_abs failed to converge; '
          write (*,*)'integral, error, tol:', rombint_abs,error, tol
        end if
        
        end function rombint_abs 

end module NonLinear


!workaround for f90 circular-module reference
     subroutine NonLinear_GetRatios(CAMB_Pk)
      use Transfer
      use NonLinear
      type(MatterPowerData) :: CAMB_Pk

      call NonLinear_GetNonLinRatios(CAMB_Pk)      

     end subroutine NonLinear_GetRatios
     

     subroutine NonLinear_GetRatios_all(CAMB_Pk)
      use Transfer
      use NonLinear
      type(MatterPowerData) :: CAMB_Pk

      call NonLinear_GetNonLinRatios(CAMB_Pk, .true.)      

     end subroutine NonLinear_GetRatios_All
     

  subroutine NonLinear_test
    use Random
    use Transfer
    type(MatterPowerData) :: Pk


      call MatterPowerData_Load(Pk,'c:\work\camb_dist\WMAP3_matterpower.dat')
 
      call NonLinear_GetRatios(Pk)

  end  subroutine NonLinear_test
