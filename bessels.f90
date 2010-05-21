!CAMB spherical and hyperspherical Bessel function routines
!This version June 2004 - bessels calculated on the fly rather than cached to disk

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Flat bessel function module

        module SpherBessels
        use Precision
        use ModelParams
        implicit none
        private
        
  !     Bessel functions and their second derivatives for interpolation
    
        real(dl), dimension(:,:), allocatable ::  ajl,ajlpr

  !     xx stores the value of x, the argument of the Bessel function,
  !     at each point the Bessel functions are calculated   
        real(dl), dimension(:), allocatable :: xx
        integer  num_xx, kmaxfile, file_numl,  file_l(lmax_arr)

!     parameters for working out where the flat Bessel functions are small
        real(dl), parameter :: xlimmin=15._dl, xlimfrac = 0.05_dl

        public ajl, ajlpr, xx, InitSpherBessels, GetBesselIndex, xlimmin, xlimfrac
       contains

         
      subroutine InitSpherBessels
!     This subroutine reads the jl files from disk (or generates them if not on disk)
      use lvalues
      implicit none
    
      !See if already loaded with enough (and correct) lSamp%l values and k*eta values
      if (allocated(ajl) .and. (lSamp%l0 <= file_numl) .and. all(file_l(1:lSamp%l0)-lSamp%l(1:lSamp%l0)==0) &
                    .and. (int(CP%Max_eta_k)+1 <= kmaxfile)) return

      !Haven't made them before, so make them now
      call GenerateBessels
    
      if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Calculated Bessels'

      end subroutine InitSpherBessels


       function GetBesselIndex(xf)
         integer GetBesselIndex
         real(dl), intent(IN) :: xf

            if (xf <= 25) then
                  if (xf <= 5) then
                     if (xf<=1) then
                      GetBesselIndex=int(100*xf)+1
                     else
                      GetBesselIndex=int(10*(xf-1))+101
                     end if
                  else
                     GetBesselIndex=int((xf-5._dl)*5)+141
                  end if
               else
                  GetBesselIndex=int((xf-25))+241
               end if

       end function GetBesselIndex


     subroutine GenerateBessels
       use lvalues
       real(dl) x
       real(dl) xlim
       integer i,j
 
        if (DebugMsgs .and. FeedbackLevel > 0) write (*,*) 'Generating flat Bessels...'
      
        file_numl= lSamp%l0 
        file_l(1:lSamp%l0) = lSamp%l(1:lSamp%l0)
        kmaxfile = int(CP%Max_eta_k)+1
     
        if (allocated(xx)) deallocate(xx)
        num_xx = (kmaxfile-25) +241 
        allocate(xx(num_xx))

        do i=1,num_xx
         if (i <= 241) then
            if (i <= 141) then
             if (i<= 101) then
              !0 to 1 - need good accuracy to get tensors right
               xx(i) = (i-1)/100._dl
             else
              !1.1 to 5
               xx(i)=(i-101)/10._dl + 1
             end if
            else
             ! 5.2 to 25
               xx(i)=(i-141)/5._dl+5
            end if
         else
            xx(i)=(i-241)+25
         end if
        end do

       if (allocated(ajl)) deallocate(ajl)
       if (allocated(ajlpr)) deallocate(ajlpr)
       Allocate(ajl(1:num_xx,1:lSamp%l0))
       Allocate(ajlpr(1:num_xx,1:lSamp%l0))

       do j=1,lSamp%l0
       
         do  i=1,num_xx
            x=xx(i)
            xlim=xlimfrac*lSamp%l(j)
            xlim=max(xlim,xlimmin)
            xlim=lSamp%l(j)-xlim
            if (x > xlim) then
               if ((lSamp%l(j)==3).and.(x <=0.2) .or. (lSamp%l(j) > 3).and.(x < 0.6) .or. &
                            (lSamp%l(j)>=5).and.(x < 1.0)) then
                   ajl(i,j)=0
                 else
                   call bjl(lSamp%l(j),x,ajl(i,j))
               end if
            else
               ajl(i,j)=0
            end if
         end do
      end do

!     get the interpolation matrix for bessel functions

      do j=1,lSamp%l0
         call spline(xx,ajl(1,j),num_xx,spl_large,spl_large,ajlpr(1,j))
      end do

     end subroutine GenerateBessels

      subroutine bjl(l,x,jl)
         use Precision
!  Calculates the spherical bessel function j_l(x)

!  and optionally its derivative for real x and integer l>=0.

!  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from

!  G.N.Watson, A Treatise on the Theory of Bessel Functions,
!  2nd Edition (Cambridge University Press, 1944).
!  Higher terms in expansion for x near l given by
!  Airey in Phil. Mag. 31, 520 (1916).

!  This approximation is accurate to near 0.1% at the boundaries
!  between the asymptotic regions; well away from the boundaries
!  the accuracy is better than 10^{-5}. The derivative accuracy
!  is somewhat worse than the function accuracy but still better
!  than 1%.

!  Point *jlp initially to a negative value to forego calculating
!  the derivative; point it to a positive value to do the derivative
!  also (Note: give it a definite value before the calculation
!  so it's not pointing at junk.) The derivative calculation requires

!  only arithmetic operations, plus evaluation of one sin() for the
!  x>>l region.


!  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu
!  This fortran version only computes j_l(x)
        implicit none
        integer l
        real(dl) nu, nu2,ax,ax2,beta,beta2,beta4,beta6
        real(dl) sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1
        real(dl) cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b
        real(dl) trigarg,trigcos,expterm,prefactor,llimit,ulimit,fl
        real(dl) x,jl

!        PI = 3.1415926536
        real(dl), parameter :: ROOTPI = 1.772453851
        real(dl), parameter :: GAMMA1 = 2.6789385347 ! Gamma function of 1/3
        real(dl), parameter :: GAMMA2 = 1.3541179394 ! Gamma function of 2/3


        ax = abs(x)
        fl = l


        beta = fl**0.325
        llimit=1.31*beta   !/* limits of asymptotic regions; fitted */
        ulimit=1.48*beta

         nu= fl + 0.5

         nu2=nu*nu

        if (l .lt. 0) then
                print*, 'Bessel function index < 0\n'
                stop
        end if

!          /************* Use closed form for l<6 **********/

        if (l .lt. 6) then

        sx=sin(ax)
        cx=cos(ax)
        ax2=ax*ax

            if(l .eq. 0) then
            if(ax .gt. 0.001) then
                 jl=sx/ax
                 else

                 jl=1._dl-ax2/6._dl
                 end if    !   /* small x trap */
            endif


            if(l .eq. 1) then

            if(ax .gt. 0.001) then
                 jl=(sx/ax -cx)/ax
                 else

                 jl=ax/3._dl
                 end if
            endif

            if(l .eq. 2) then
              if(ax .gt. 0.001) then
                  jl=(-3._dl*cx/ax &
                          -sx*(1._dl-3._dl/ax2))/ax
              else

                  jl=ax2/15._dl
                  end if
            endif

            if(l .eq. 3) then
          if(ax .gt. 0.001) then
                jl=(cx*(1._dl-15._dl/ax2) &
                          -sx*(6._dl-15._dl/ax2)/ax)/ax
           else

                jl=ax*ax2/105._dl
                endif
            endif

            if(l .eq. 4) then
          if(ax .gt. 0.001) then
        jl=(sx*(1._dl-45._dl/(ax*ax)+105._dl &
             /(ax*ax*ax*ax)) +cx*(10._dl-105._dl/(ax*ax))/ax)/ax
          else

                jl=ax2*ax2/945._dl
                end if
            endif

             if(l .eq. 5) then

          if(ax .gt. 0.001) then
        jl=(sx*(15._dl-420._dl/(ax*ax)+945._dl &
          /(ax*ax*ax*ax))/ax -cx*(1.0-105._dl/(ax*ax)+945._dl &
                                          /(ax*ax*ax*ax)))/ax
           else

                jl=ax2*ax2*ax/10395._dl
                endif
             endif


!          /********************** x=0 **********************/

        else if (ax .lt. 1.d-30) then
        jl=0.0

!          /*************** Region 1: x << l ****************/

        else if (ax .le. fl+0.5-llimit) then


!       beta=acosh(nu/ax)
        if (nu/ax .lt. 1._dl) print*, 'trouble with acosh'
        beta = dlog(nu/ax + sqrt((nu/ax)**2 - 1._dl) )
                !(4.6.21)
        cotb=nu/sqrt(nu*nu-ax*ax)      ! /* cotb=coth(beta) */
        cot3b=cotb*cotb*cotb
        cot6b=cot3b*cot3b
        secb=ax/nu
        sec2b=secb*secb
        sec4b=sec2b*sec2b
        sec6b=sec4b*sec2b
        sum1=2.0+3.0*sec2b
        expterm=sum1*cot3b/(24.0*nu)
        sum2=4.0+sec2b
        expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2)
        sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
        expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2)
        sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
        expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
        expterm=exp(-nu*beta+nu/cotb-expterm)
        prefactor=sqrt(cotb/secb)/(2.0*nu)
        jl=prefactor*expterm

!          /**************** Region 2: x >> l ****************/


        else if (ax .ge. fl+0.5+ulimit) then


        beta=acos(nu/ax)
        cotb=nu/sqrt(ax*ax-nu*nu)      !/* cotb=cot(beta) */
        cot3b=cotb*cotb*cotb
        cot6b=cot3b*cot3b
        secb=ax/nu
        sec2b=secb*secb
        sec4b=sec2b*sec2b
        sec6b=sec4b*sec2b
        trigarg=nu/cotb - nu*beta - PI/4.0
        sum1=2.0+3.0*sec2b
        trigarg = trigarg - sum1*cot3b/(24.0*nu)
        sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
        trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2)
        trigcos=cos(trigarg)
        sum2=4.0+sec2b
        expterm=sum2*sec2b*cot6b/(16.0*nu2)
        sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
        expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
        expterm=exp(-expterm)
        prefactor=sqrt(cotb/secb)/nu
        jl=prefactor*expterm*trigcos

!          /***************** Region 3: x near l ****************/

        else



        beta=ax-nu

        beta2=beta*beta
        beta4=beta2*beta2
        beta6=beta2*beta4
        sx=6.0/ax
        sx2=sx*sx
        cx=sqrt(sx)

        secb=sx**0.333333333

        sec2b=secb*secb

        deriv1=GAMMA1*secb
        deriv1= deriv1+ beta*GAMMA2*sec2b
        sum1=(beta2/6.0-1.0/15.0)*beta
        deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0
        sum2=beta4/24.0-beta2/24.0+1.0/280.0
        deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0
        sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0
        deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0
        sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0- &
        13.0/31500.0)*beta
        deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0
        sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0 &
                     -121.0*beta2/907200.0 + 7939.0/232848000.0)*beta
        deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0

        jl=deriv1*cx/(12.0*ROOTPI)

        end if

        end subroutine bjl
     

     end module SpherBessels




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Calculation of ultraspherical Bessel functions.                      c
! Fortran version of the c program hyperjl.c by Arthur Kosowsky.       c
! WKB approx described in astro-ph/9805173                             c
!                                                                      c
! Modifications by Anthony Challinor and Antony Lewis                  c
! Minor modifications to correct K=1 case outside [0,pi],              c
! the small chi approximations for lSamp%l=0 and lSamp%l=1, and                    c
! the quadratic approximation to Q(x) around Q(x)=0.                   c
! Bug fixed in downwards recursion (phi_recurs)                        c
!                                                                      c
! The routine phi_recurs uses recursion relations to calculate         c
! the functions, which is accurate but relatively slow.                c  
!   ***NOT STABLE FOR K=1 or for all cases ***                         c
!                                                                      c
! The routine phi_langer uses Langer's formula for a                   c
! uniform first-order asymptotic approximation in the open, closed     c
! and flat cases. This approximation is EXCELLENT for all lSamp%l >= 3.      c
!                                                                      c
! The routine qintegral calculates the closed-form answer              c
! to the eikonal integral used in the WKB approximation.               c
!                                                                      c
! The routine airy_ai returns the Airy function Ai(x) of the argument  c
! passed. It employs a Pade-type approximation away from zero and      c
! a Taylor expansion around zero. Highly accurate.                     c
!                                                                      c
! The routines polevl and p1evl are auxiliary polynomial               c
! evaluation routines used in the airy function calculation.           c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



 module USpherBessels
 use Precision
 implicit none
 private 
 
 
 public USpherBesselWithDeriv, phi_recurs,phi_langer
 
 
 contains

        subroutine USpherBesselWithDeriv(closed,Chi,l,beta,y1,y2)
          !returns y1=ujl*sinhChi and y2=diff(y1,Chi)
          !aim for accuracy > 1% for all inputs
        real(dl) Chi,beta,y1,y2,sinhChi,cothChi
        real(dl) sin_K, cot_K
        integer l,K
        logical, intent(IN) :: closed
        logical DoRecurs
        
        if (closed) then
            sin_K = sin(Chi) 
            cot_K= 1._dl/tan(Chi)
            K=1
         else 
            sin_K=sinh(Chi)
            cot_K = 1._dl/tanh(Chi)
            K=-1
         end if
      
        sinhChi = sin_K
        cothChi = cot_K
   
       DoRecurs = ((l<=30).OR.((.not.closed.or.(abs(Chi-pi/2)>0.2d0)).and.(beta*l<750) &
            .or.closed.and.(beta*l<4000)))

!Deep in the tails the closed recursion relation is not stable
!Added July 2003 to prevent problems with very nearly flat models
       if (DoRecurs .and. closed .and. Chi < asin(sqrt(l*(l+1._dl))/beta) - 2/beta) then
           if (phi_langer(l,K,beta,Chi) < 1e-7) then
             call phi_small_closed_int(l,beta,chi,y1,y2)
             return
           end if 
     end if        

     if (DoRecurs) then 
             !use recursive evaluation where WKB is poor and recurs is fast anyway.
         y1=phi_recurs(l,K,beta,Chi)*sinhChi
         y2=y1*(l+1)*cothChi
         if (.not.closed.or.(l+1<nint(beta))) y2=y2 - &
              sqrt(beta**2-(K*(l+1)**2))*phi_recurs(l+1,K,beta,Chi)*sinhChi
          !of course we could get y2 much more quickly by modifying
          !phi_recurs to return l and l+1 for each beta,Chi...


         else !WKB approx
         y1=phi_langer(l,K,beta,Chi)*sinhChi
         y2=y1*(l+1)*cothChi   
         if (.not.closed.or.(l+1<nint(beta))) y2=y2 - &
              sqrt(beta**2-(K*(l+1)**2))*phi_langer(l+1,K,beta,Chi)*sinhChi  
         
         end if
        
        end subroutine USpherBesselWithDeriv



!Calculates y1,y2 (noramlized to a value near the turning point)  
!by integrating up the differential equation and normalizing to phi_recurs
!in the region in which phi_recurs is stable
!This allows closed functions to be computed where chi << turning point
   subroutine phi_small_closed_int(l,beta,chi,y1,y2)
   integer, intent(IN) :: l
   real(dl), intent(IN) :: beta, chi
   real(dl) y1,y2
   
   integer nsteps,i
  
   real(dl) ap1,nu2,dydchi1,dydchi2,yt1,yt2,dyt1,dyt2,dym1,dym2
   real(dl) x0, delchi,sh, h6,x
   real(dl) y1_x,y2_x, tmp,xh,hh
 
   nsteps = 200
   ap1 = l*(l+1)
   x0 = sqrt(ap1)/beta
   nu2 = beta**2

   if ((beta*chi)**2/l < 0.005) then
 !Series solution
    
     x = chi
     sh = sin(x)
     tmp=(ap1/sh**2 - nu2)
     y1=1e-20  
     y2 = ((l+1)/x - (nu2-ap1/3)/(2*l+3)*x) * y1
   else

   x = max(1d-7,chi - 50._dl/l)
   y1=1e-20  
   y2 = (l+1)*y1/x
   
   delchi = (chi-x)/nSteps
   h6=delchi/6
   hh=delchi/2
   sh = sin(x)
   tmp=(ap1/sh**2 - nu2) 
   
          do i=1,nSteps          
! One step in the ujl integration
! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2     
            xh=x+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/sin(xh)**2 - nu2) 
            
            
            dyt2=tmp*yt1       !deriv y2 at mid
       
            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid
           
            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+delchi*dym1 !y1 at end
            dym1=dyt1+dym1     
            yt2=y2+delchi*dym2 !y2 at end
            dym2=dyt2+dym2
            
            x=x+delchi     !end point
            sh=sin(x)    
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2*dym2)       
           if (y1 > 1d10 .or. y2> 1d10) then
         
              y1=y1/1d10
              y2=y2/1d10
        
           end if
         end do

      end if     
   
       y1_x = y1; y2_x = y2 

       delchi = (x0 - chi)/nSteps
       h6=delchi/6
       hh=delchi/2

          do i=1,nSteps          
! One step in the ujl integration
! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2     
            xh=x+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/sin(xh)**2 - nu2) 
            
            
            dyt2=tmp*yt1       !deriv y2 at mid
       
            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid
           
            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+delchi*dym1 !y1 at end
            dym1=dyt1+dym1     
            yt2=y2+delchi*dym2 !y2 at end
            dym2=dyt2+dym2
            
            x=x+delchi     !end point
            sh=sin(x)    
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2*dym2)       
            if (y1 > 1d10 .or. y2 > 1d10) then
              y1=y1/1d10
              y2=y2/1d10
              y1_x = y1_x/1d10
              y2_x = y2_x/1d10
        
            end if
         end do


         tmp = phi_recurs(l,1,beta,x0)*sin(x0) / y1
         y1 = y1_x * tmp
         y2 = y2_x * tmp


   end subroutine phi_small_closed_int
 
!***********************************************************************
!                                                                      *
! Calculates Phi(l,beta,chi) using recursion on l.                     *
! See Abbot and Schaefer, ApJ 308, 546 (1986) for needed               *
! recursion relations and closed-form expressions for l=0,1.           *
! (Note: Their variable y is the same as chi here.)                    *
!                                                                      *
! When the flag direction is negative, downwards recursion on l        *
! must be used because the upwards direction is unstable to roundoff   *
! errors. The downwards recursion begins with arbitrary values and     *
! continues downwards to l=1, where the desired l value is normalized  *
! using the closed form solution for l=1. (See, e.g., Numerical        *
! Recipes of Bessel functions for more detail)                         *
!                                                                      *
!***********************************************************************

  function phi_recurs(l, K, beta, chi)
  !doesn't like values which give exponentially small phi
  integer, intent(IN) :: l, K
  real(dl), intent(IN) :: beta, chi
  real(dl) phi_recurs
  integer j, direction, lstart,ibeta
  real(dl) ell, kay, arg, answer,beta2
  real(dl) root_K
  real(dl) phi0, phi1, phi_plus, phi_zero, phi_minus, b_zero, b_minus 
  real(dl), parameter :: ACC=40._dl, BIG=1.d10
  real(dl) sin_K, cot_K

  ell=dble(l)
 
 
  ! Test input values

  if(l<0) then
     write(*,*) "Bessel function index ell < 0"
     stop
  endif
  if(beta<0._dl) then
     write(*,*) "Wavenumber beta < 0"
     stop
  endif
  if ((abs(K)/=1).and.(K/=0)) then
     write(*,*) "K must be 1, 0 or -1"
     stop
  end if
  
  if(K==1) then    
     ibeta=nint(beta)
     if(ibeta<3) then
        write(*,*) "Wavenumber beta < 3 for K=1"
        stop
     endif
     if(ibeta<=l) then
        write(*,*) "Wavenumber beta <= l"
        stop
     endif
  endif
 
  if (chi<1/BIG) then
     phi_recurs=0
     return
  end if
 
  kay = dble(K)
  arg = beta * chi
  beta2 = beta**2

  if(K == 0) then
     cot_K = 1._dl/chi
     sin_K = chi
     root_K = beta
  else 
   root_K = sqrt(beta2 -kay*ell*ell)
  
    if(K == -1) then
    cot_K = 1._dl/tanh(chi)
    sin_K = sinh(chi)   
    else
    cot_K = 1._dl/tan(chi)
    sin_K = sin(chi)   
    end if
  
  endif
 
  
  ! Closed form solution for l=0
 
  if (abs(chi) < 1.d-4) then
    if (abs(arg)<1.d-4) then
     phi0 = 1._dl-chi**2*(beta*beta-kay)/6._dl
    else
     phi0=sin(arg)/arg
    end if
  else     
     phi0 = sin(arg) / (beta * sin_K)
  end if

  if (l==0) then
       phi_recurs=phi0
       return 
   end if


     ! Closed form solution for l=1

     if((abs(chi) < 1.d-4).and.(K/=0)) then
        if(arg < 1.d-4) then
           phi1 = chi*sqrt(beta*beta-kay)/3._dl
           !beta2 * chi / (3._dl * sqrt(1._dl+ kay * beta2))
        else 
           phi1 = (sin(arg)/arg-cos(arg))/(sqrt(beta*beta-kay)*chi)
             !(sin(arg)/arg - cos(arg))/arg
        end if
     elseif ((abs(arg) < 1.d-4).and.(K == 0)) then
        phi1 = arg / 3._dl
     else 
        if (K /= 0 ) then
           phi1 = sin(arg) * cot_K / (beta * sin_K) - cos(arg) / sin_K
           phi1 = phi1/sqrt(beta2 - kay)
        else
           phi1 = (sin(arg)/arg - cos(arg))/arg
        end if
     end if
     if(l==1) then
        phi_recurs=phi1
        return
     end if
     ! Find recursion direction
     !  direction = +1 for upward recursion, -1 for downward

     if(abs(cot_K) < root_K / ell) then 
        direction = 1
     else 
        direction = -1
     end if
    
     ! For K=1, must do upwards recursion:
     ! NOT STABLE for all values of chi

     if(K==1) direction = 1

     ! Do upwards recursion on l
    
     
     if(direction == 1)then
        b_minus = sqrt(beta2 - kay)
        phi_minus = phi0
        phi_zero = phi1
      
        do j=2,l 
           
           if(K == 0) then
              phi_plus = ((2*j-1) * cot_K * phi_zero - beta*phi_minus)/ beta
           else 
              b_zero = sqrt(beta2 - (K*j*j))  
              phi_plus = ((2*j-1) * cot_K * phi_zero - b_minus * phi_minus) / b_zero
              b_minus = b_zero
           end if
           phi_minus = phi_zero
           phi_zero = phi_plus
        end do
      
      
        phi_recurs=phi_plus
        
     
        return

        ! Do downwards recursion on l

     else 
        lstart = l + 2 * int(sqrt(ell*ACC))
      
        b_zero = sqrt(beta2 - dble(K*lstart*lstart))
        phi_plus = 0._dl
        phi_zero = 1._dl
        answer = 0._dl
       
        do j= lstart - 2,1,-1
           
           if(K == 0) then
              phi_minus = ((2*j + 3) * cot_K * phi_zero - beta * phi_plus) / beta
           else
              b_minus = sqrt(beta2 - (K*(j+1)**2))  
              phi_minus = ((2*j + 3) * cot_K * phi_zero - b_zero * phi_plus) / b_minus
              b_zero = b_minus
           end if
           phi_plus = phi_zero
           phi_zero = phi_minus

           if(j == l) answer = phi_minus
           if((abs(phi_zero) > BIG).and.(j/=1)) then
              phi_plus = phi_plus/BIG
              phi_zero = phi_zero/BIG            
              answer = answer/BIG
           end if
        end do

        ! Normalize answer to previously computed phi1

        answer = answer*phi1 / phi_minus
        phi_recurs=answer    
       
     end if
 
   end function phi_recurs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Calculates Phi(l,beta,chi) using the Langer uniform approximation    c
! to the first-order WKB approximation.                                c
! See C.M. Bender and S.A. Orszag,  Mathematical Methods for           c
! Scientists and Engineers (McGraw-Hill, 1978; LC QA371.B43),          c
! chapter 10.                                                          c
!                                                                      c
! Differential equation for needed function can be cast into the       c
! Schrodinger form      \epsilon^2 y'' = Q(x) y                        c
! where \epsilon^2 = 1/l(l+1) and Q(x) depends on the parameter        c
! alpha \equiv beta * epsilon.                                         c
!                                                                      c
! In the K= +1 case, the function is                                   c
! determined by its value on the interval [0, pi/2] and the symmetry   c
! conditions Phi(chi + pi) = (-1)^{beta - l - 1} Phi(chi),             c
!            Phi(pi - chi) = (-1)^{beta - l - 1} Phi(chi).             c
! This interval contains one turning point, so the Langer formula      c
! can be used.                                                         c
! Note that the second condition at chi = pi/2 gives an eigenvalue     c
! condition on beta, which  must corrected. For the lowest             c
! eigenvalue(s), the region between the turning points is not large    c
! enough for the asymptotic solution to be valid, so the functions     c
! have a small discontinuity or discontinuous derivative at pi/2;      c
! this behavior is corrected by employing a 4-term asymptotic          c
! series around the regular point chi=pi/2.                            c
! The exact eigenvalue condition requires that beta must be an         c
! integer >= 3 with beta > l. Note this implies alpha > 1.             c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      


      function phi_langer(l,K,beta,chi)
      integer l,K,ibeta,kay
      real(dl) phi_langer
      real(dl) ell,symm, anu, alpha2
      real(dl) beta,chi,eikonal, wkb, arg, arg2, tmp
      real(dl) epsilon, alpha, chi0, x, a, b,achi 
    
      real(dl) cot_K, sin_K
      real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
                           PIOVER2=1.570796327d0

      ell=dble(l)
      achi=chi
   
      symm=1._dl
!
! Test input values
!
      if(l<0) then
         write(*,*) "Bessel function index ell < 0"
         stop
      endif
      if(beta<0._dl) then
         write(*,*) "Wavenumber beta < 0"
         stop
      endif
      if ((abs(K)/=1).and.(K/=0)) then
        write(*,*) "K must be 1, 0 or -1"
        stop
      end if

      
      if(K == 1) then        
         ibeta=nint(beta)
         if(ibeta<3) then
            write(*,*) "Wavenumber beta < 3 for K=1"
            stop
         endif
         if(ibeta<=l) then
            write(*,*) "Wavenumber beta <= l"
            stop
         endif
      endif
   
      kay=K


! For closed case, find equivalent chi in [0,pi/2]
!
      if(K==1) then
         achi=achi-2._dl*Pi*int(achi/2._dl/PI)
         if(achi>PI) then
            achi=2._dl*PI-achi
            if(2*(l/2).eq.l) then
               symm=symm
            else
               symm=-symm
            endif
         endif
         if(achi>PI/2._dl) then
            achi=PI-achi
            if(2*((ibeta-l-1)/2).eq.(ibeta-l-1)) then
               symm=symm
            else
               symm=-symm
            endif
         endif
      endif

! Definitions
      if(K == 0) then
          sin_K = achi         
      else 
            if(K == -1) then
               sin_K = sinh(achi)   
             else
               sin_K = sin(achi)   
            end if           
      endif 

! Closed form solution for l=0
!
      if(l == 0) then
         arg=beta*achi
       
         if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
               wkb=1._dl-achi*achi*(beta*beta-kay)/6._dl
            else
               wkb=sin(arg)/arg
            endif
         else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=1._dl-arg*arg/6._dl
         else
            wkb=sin(arg)/(beta*sin_K)
         endif
         phi_langer=symm*wkb
         return
      endif
!
! Closed form solution for l=1
!
      if(l==1) then
         arg=beta*achi
      
         if((abs(achi)<1.d-4).and.(K/=0)) then
            if(abs(arg)<1.d-4) then
               wkb=achi*sqrt(beta*beta-kay)/3._dl
            else
               wkb=(sin(arg)/arg-cos(arg))/(sqrt(beta*beta-kay)*achi)
            endif
         else if((abs(arg)<1.d-4).and.(K==0)) then
            wkb=arg/3._dl
         else 
          if(K/=0) then
            if(K==1) then
               cot_K=1._dl/tan(achi)
            else
               cot_K=1._dl/tanh(achi)
            endif
            wkb=sin(arg)*cot_K/(beta*sin_K)-cos(arg)/sin_K
            wkb=wkb/sqrt(beta*beta-kay)
         else
            wkb=(sin(arg)/arg-cos(arg))/arg
         endif
         end if
         phi_langer=symm*wkb
         return
      endif
!
! Closed form solution for K=1 and beta = l+1 (lowest eigenfunction)
!
      if((K==1).and.(ibeta == (l+1))) then
         wkb=(sin_K**ell)* &
              sqrt(sqrt(2._dl*PI/(2._dl*ell+1._dl))*ell/((ell+1._dl)*(2._dl*ell+1._dl)))
         wkb=wkb*(1+0.1875d0/ell-0.013671875/(ell*ell))
         phi_langer=symm*wkb
         return
      endif

! Very close to 0, return 0 (exponentially damped)
!
      if(abs(achi)<1.d-8) then
         phi_langer=0._dl
         return
      endif


! For closed case, find corrected eigenvalue beta
!
      if(K==1) then
         anu=dble(ibeta)-1._dl/(8._dl*ell)+1._dl/(16._dl*ell*ell)
      else
         anu=beta
      endif
!
! Evaluate epsilon using asymptotic form for large l
!
      if(l<20) then
         epsilon=1._dl/sqrt(ell*(ell+1._dl))
      else
         epsilon=1._dl/ell-0.5d0/(ell*ell)+0.375d0/(ell*ell*ell)
      endif
     
      alpha=epsilon*anu
!
! Calculate the turning point where Q(x)=0.
! Function in question has only a single simple turning point.
!
      if(K==-1) chi0=log((1._dl+sqrt(1._dl+alpha*alpha))/alpha)
      if(K==0) chi0=1._dl/alpha
      if(K==1) chi0=asin(1._dl/alpha)
 

! Very close to chi0, use usual wkb form to avoid dividing by zero
!
      if(abs(achi-chi0)<1.d-5) then

! Calculate coefficients of linear and quadratic terms in Q(x) expansion
! in the neighborhood of the turning point
! Q(chi)=a*(chi0-chi)+b*(chi0-chi)**2
         alpha2=alpha*alpha
       
         if(K==-1) then
            a=2._dl*alpha2*sqrt(alpha2+1._dl)
            b=3._dl*alpha2**2+2._dl*alpha2
         endif
         if(K==0) then
            a=2._dl*alpha2*alpha
            b=3._dl*alpha2**2
         endif
         if(K==1) then
            a=2._dl*alpha2*sqrt(alpha2-1._dl)
            b=3._dl*alpha2**2-2._dl*alpha2
         endif
       
! Dependent variable x for which Q(x)=0 at x=0
! x>0 is the evanescent region
!
         x=chi0-achi
!
! Argument of Airy function
!
         arg=(x+b*x*x/(5._dl*a))/(epsilon*epsilon/a)**(0.333333333d0)
!
! Evaluate Airy function
!
         wkb=airy_ai(arg)
!
! Rest of functional dependence
!
         wkb=wkb*(1._dl-b*x/(5._dl*a))/sin_K 
!  Normalization factor:

         wkb=wkb*symm*ROOTPI*((a*epsilon)**(-0.1666667d0))*sqrt(epsilon/anu)
         phi_langer=wkb  
        
         return
      endif


! Langer approximation.
!
! Transport factor:
!
      tmp=sqrt(abs(1._dl/(sin_K*sin_K)-alpha*alpha))
     
! Eikonal factor
!
      eikonal=qintegral(sin_K,alpha,K)

    
     
      arg=(1.5d0*eikonal/epsilon)**(1._dl/3._dl)
   
      arg2=arg*arg
      if(achi>chi0) arg2=-arg2
      
! Multiply factors together

      wkb=airy_ai(arg2)*symm*ROOTPI*sqrt(arg*epsilon/anu/tmp)/Sin_K
      phi_langer=wkb

      end function phi_langer

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                       c
! Evaluates the exact solution to  the integral giving the WKB          c
! eikonal solution,   \int^x sqrt(abs(Q(x))) dx                         c
!                                                                       c
! In the open case, this integral costs 1 or 2 square roots, an atan    c
! and a log; its evaluation will be roughly as expensive as the rest    c
! of the Phi routine. An analytic fit cannot be substantially faster    c
! because the dependence on alpha of the y-intercept of the linear      c
! region of the integrand contains a log and an atan, so at best a fit  c
! can only save the computation of the square roots.                    c
!                                                                       c
! The integrals are very bland functions of chi and alpha and could     c
! be precomputed and cached to save computation time; interpolation     c
! on a relatively small number of points should be very accurate        c
!                                                                       c
! Note that for the closed case, the variable arg must be between 0     c
! and alpha; for arg > alpha, symmetry properties reduce the needed     c
! integral to this case.                                                c
!                                                                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function qintegral(sin_K,alpha, K)
      implicit none
      real(dl) qintegral, sin_K
      integer K
      real(dl) alpha,exact,arg, root1, root2, dummyarg
     
      real(dl), parameter :: PI=3.1415926536d0,ROOTPI=1.772453851d0,ROOT2PI=2.506628275d0, &
                           PIOVER2=1.570796327d0

      arg=alpha*sin_K

      if(K==0) then
         if(arg>1._dl) then
            exact=sqrt(arg*arg-1._dl)-acos(1._dl/arg)
            qintegral=exact
            return
         else
            root1=sqrt(1._dl-arg*arg)
            exact=log((1._dl+root1)/arg)-root1
            qintegral=exact
            return
         endif
      else if(K==-1) then
         if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(arg*arg+alpha*alpha)
            exact=alpha/2._dl*log((2._dl*arg*arg+alpha*alpha-1._dl+ &
                  2._dl*root1*root2)/(1._dl+alpha*alpha))+atan(root2/ &
                  (alpha*root1))-PIOVER2
            qintegral=exact
            return
         else
            root1=sqrt((1._dl-arg*arg)*(arg*arg+alpha*alpha))
            exact=alpha/2._dl*atan(-2._dl*root1/(2*arg*arg+alpha*alpha- &
                  1._dl))+0.5d0*log((2._dl*alpha*root1+2._dl*alpha*alpha+ &
                  arg*arg*(1._dl-alpha*alpha))/(arg*arg*(1._dl+ &
                  alpha*alpha)))
            if(2._dl*arg*arg+alpha*alpha-1._dl<0._dl) then
               exact=exact-alpha*PIOVER2
            endif
            qintegral=exact
            return
         endif
      else
         if(arg>1._dl) then
            root1=sqrt(arg*arg-1._dl)
            root2=sqrt(alpha*alpha-arg*arg)
            exact=alpha/2._dl*atan(-2._dl*root1*root2/ &
               (2._dl*arg*arg-alpha*alpha-1._dl))- &
               atan(-root2/(root1*alpha))-PIOVER2
            if(2._dl*arg*arg-alpha*alpha-1._dl>0._dl) then
               exact=exact+alpha*PIOVER2
            endif
         else
            root1=sqrt((1._dl-arg*arg)*(alpha*alpha-arg*arg))
            dummyarg=alpha*(1._dl-arg*arg)/root1
            exact=0.5d0*log((1._dl+dummyarg)/(1._dl-dummyarg))- &
             alpha/2._dl*log((alpha*alpha-2._dl*arg*arg+1._dl+2._dl*root1)/ &
             (alpha*alpha-1._dl))
         endif
         qintegral=exact
         return
      endif
   
      end function qintegral 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!       Airy function                                                  c
!                                                                      c
! Modified from original routine by Stephen Moshier, available         c
! as part of the Cephes library at www.netlib.com                      c
! Modifications: eliminates calculation of Bi(x), Ai'(x), Bi'(x)       c
! and translation to Fortran                                           c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! DESCRIPTION:
!
! Solution of the differential equation
!
!       y"(x) = xy.
!
! The function returns the two independent solutions Ai, Bi
! and their first derivatives Ai'(x), Bi'(x).
!
! Evaluation is by power series summation for small x,
! by rational minimax approximations for large x.
!
!
!
! ACCURACY:
! Error criterion is absolute when function <= 1, relative
! when function > 1, except * denotes relative error criterion.
! For large negative x, the absolute error increases as x^1.5.
! For large positive x, the relative error increases as x^1.5.
!
! Arithmetic  domain   function  # trials      peak         rms
! IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
! IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
! IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
! IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
! IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
! IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
! DEC         -10, 0     Ai         5000       1.7e-16     2.8e-17
! DEC           0, 10    Ai         5000       2.1e-15*    1.7e-16*
! DEC         -10, 0     Ai'        5000       4.7e-16     7.8e-17
! DEC           0, 10    Ai'       12000       1.8e-15*    1.5e-16*
! DEC         -10, 10    Bi        10000       5.5e-16     6.8e-17
! DEC         -10, 10    Bi'        7000       5.3e-16     8.7e-17
!
!
! Cephes Math Library Release 2.1:  January, 1989
! Copyright 1984, 1987, 1989 by Stephen lSamp%l. Moshier
! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
!

      function airy_ai(x)
      implicit none
      real(dl) airy_ai
      real(dl) x,z, zz, t, f, g, uf, ug, zeta, theta
      real(dl) ak
      real(dl) AN(8),AD(8),AFN(9),AFD(9),AGN(11),AGD(10)
      real(dl), parameter :: AMAXAIRY=25.77d0,ACC=1.d-8,PI=3.1415926536d0
      real(dl), parameter :: c1=0.35502805388781723926d0, c2=0.258819403792806798405d0
      real(dl), parameter :: sqrt3=1.732050807568877293527d0,sqpii=5.64189583547756286948d-1
     
      
      AN(1)=3.46538101525629032477d-1
      AN(2)=1.20075952739645805542d1
      AN(3)=7.62796053615234516538d1
      AN(4)=1.68089224934630576269d2
      AN(5)=1.59756391350164413639d2
      AN(6)=7.05360906840444183113d1
      AN(7)=1.40264691163389668864d1
      AN(8)=9.99999999999999995305d-1

      AD(1)=5.67594532638770212846d-1
      AD(2)=1.47562562584847203173d1
      AD(3)=8.45138970141474626562d1
      AD(4)=1.77318088145400459522d2
      AD(5)=1.64234692871529701831d2
      AD(6)=7.14778400825575695274d1
      AD(7)=1.40959135607834029598d1
      AD(8)=1.00000000000000000470d0

      AFN(1)=-1.31696323418331795333d-1
      AFN(2)=-6.26456544431912369773d-1
      AFN(3)=-6.93158036036933542233d-1
      AFN(4)=-2.79779981545119124951d-1
      AFN(5)=-4.91900132609500318020d-2
      AFN(6)=-4.06265923594885404393d-3
      AFN(7)=-1.59276496239262096340d-4
      AFN(8)=-2.77649108155232920844d-6
      AFN(9)=-1.67787698489114633780d-8

      AFD(1)=1.33560420706553243746d1
      AFD(2)=3.26825032795224613948d1
      AFD(3)=2.67367040941499554804d1
      AFD(4)=9.18707402907259625840d0
      AFD(5)=1.47529146771666414581d0
      AFD(6)=1.15687173795188044134d-1
      AFD(7)=4.40291641615211203805d-3
      AFD(8)=7.54720348287414296618d-5
      AFD(9)=4.51850092970580378464d-7

      AGN(1)=1.97339932091685679179d-2
      AGN(2)=3.91103029615688277255d-1
      AGN(3)=1.06579897599595591108d0
      AGN(4)=9.39169229816650230044d-1
      AGN(5)=3.51465656105547619242d-1
      AGN(6)=6.33888919628925490927d-2
      AGN(7)=5.85804113048388458567d-3
      AGN(8)=2.82851600836737019778d-4
      AGN(9)=6.98793669997260967291d-6
      AGN(10)=8.11789239554389293311d-8
      AGN(11)=3.41551784765923618484d-10

      AGD(1)=9.30892908077441974853d0
      AGD(2)=1.98352928718312140417d1
      AGD(3)=1.55646628932864612953d1
      AGD(4)=5.47686069422975497931d0
      AGD(5)=9.54293611618961883998d-1
      AGD(6)=8.64580826352392193095d-2
      AGD(7)=4.12656523824222607191d-3
      AGD(8)=1.01259085116509135510d-4
      AGD(9)=1.17166733214413521882d-6
      AGD(10)=4.91834570062930015649d-9
!
! Exponentially tiny for large enough argument
!
      if(x>AMAXAIRY) then
         airy_ai=0._dl
         return
      endif
!
! Pade fit for large negative arguments
!
      if(x<-2.09d0) then
         t=sqrt(-x)
         zeta=-2._dl*x*t/3._dl
         t=sqrt(t)
         ak=sqpii/t
         z=1._dl/zeta
         zz=z*z
         uf=1._dl+zz*polevl(zz,AFN,8)/p1evl(zz,AFD,9)
         ug=z*polevl(zz,AGN,10)/p1evl(zz,AGD,10)
         theta=zeta+0.25d0*PI
         f=sin(theta)
         g=cos(theta)
         airy_ai=ak*(f*uf-g*ug)
         return
      endif
!
! Pade fit for large positive arguments
!
      if(x>=2.09) then
         t=sqrt(x)
         zeta=2._dl*x*t/3._dl
         g=exp(zeta)
         t=sqrt(t)
         ak=2._dl*t*g
         z=1._dl/zeta
         f=polevl(z,AN,7)/polevl(z,AD,7)
         airy_ai=sqpii*f/ak
         return
      endif
!
! Taylor series for region around x=0
!

      f=1._dl
      g=x
      t=1._dl
      uf=1._dl
      ug=x
      ak=1._dl
      z=x*x*x
     do while (t>ACC) 
         uf=uf*z
         ak=ak+1._dl
         uf=uf/ak
         ug=ug*z
         ak=ak+1._dl
         ug=ug/ak
         uf=uf/ak
         f=f+uf
         ak=ak+1._dl
         ug=ug/ak
         g=g+ug
         t=abs(uf/f)
      end do
     

      airy_ai=c1*f-c2*g
      
      end function airy_ai

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Evaluate polynomial                                          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! DESCRIPTION:
!
! Evaluates polynomial of degree N:
!
!                     2          N
! y  =  C  + C x + C x  +...+ C x
!        0    1     2          N
!
! Coefficients are stored in reverse order:
!
! coef(1) = C  , ..., coef(N+1) = C  .
!            N                     0
!
! The function p1evl() assumes that C = 1.0 and is
!                                    N
! omitted from the array.  Its calling arguments are
! otherwise the same as polevl().
!
!
! SPEED:
!
! In the interest of speed, there are no checks for out
! of bounds arithmetic.  This routine is used by most of
! the functions in the library.  Depending on available
! equipment features, the user may wish to rewrite the
! program in microcode or assembly language.
!
! Cephes Math Library Release 2.1:  December, 1988
! Copyright 1984, 1987, 1988 by Stephen lSamp%l. Moshier
! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
!
      function polevl(x,coef,N)
      implicit none
      real(dl) polevl
      real(dl) x,ans
      real(dl) coef
      integer N,i

      dimension coef(N+1)

      ans=coef(1)
      do i=2,N+1
         ans=ans*x+coef(i)
      end do
      polevl=ans
      
      end function polevl

!
!                                      N
! Evaluate polynomial when coefficient of x  is 1.0.
! Otherwise same as polevl.
!
      function p1evl(x,coef,N)
      implicit none
      real(dl) p1evl
      real(dl) x,coef,ans
      integer N,i
      dimension coef(N)

      ans=x+coef(1)
      do i=2,N
         ans=ans*x+coef(i)
      end do
      p1evl=ans
     
      end function p1evl



  end module USpherBessels



