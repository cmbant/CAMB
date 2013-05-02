!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The `halofit' code models the nonlinear evolution of cold matter 
! cosmological power spectra. The full details of the way in which 
! this is done are presented in Smith et al. (2002), MNRAS, ?, ?. 
!
! The code `halofit' was written by R. E. Smith & J. A. Peacock. 
! See http://www.astro.upenn.edu/~res, 
! Last edited 8/5/2002.

! Only tested for plain LCDM models with power law initial power spectra

! Adapted for F90 and CAMB, AL March 2005
!!BR09 Oct 09: generalized expressions for om(z) and ol(z) to include w

! RT12 Oct: update some fitting parameters in the code to enhance
!           the power spectrum at small scales (arXiv:1208.2701)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module NonLinear
      use ModelParams
      use transfer
      use LambdaGeneral
      implicit none 
      private
       
       real, parameter :: Min_kh_nonlinear = 0.005
       real(dl):: om_m,om_v,fnu,omm0

      public Min_kh_nonlinear,NonLinear_GetNonLinRatios
      
      contains

     subroutine NonLinear_GetNonLinRatios(CAMB_Pk)
     !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
     !for each redshift and wavenumber
     !This implementation uses Halofit
      type(MatterPowerData) :: CAMB_Pk
      integer itf
      real(dl) a,plin,pq,ph,pnl,rk
      real(dl) sig,rknl,rneff,rncur,d1,d2
      real(dl) diff,xlogr1,xlogr2,rmid
      integer i

       !!BR09 putting neutrinos into the matter as well, not sure if this is correct, but at least one will get a consisent omk.
       omm0 = CP%omegac+CP%omegab+CP%omegan
       fnu = CP%omegan/omm0

       CAMB_Pk%nonlin_ratio = 1

       do itf = 1, CAMB_Pk%num_z


! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and 
! curvature (rncur) of the power spectrum at the desired redshift, using method 
! described in Smith et al (2002).
       a = 1/real(1+CAMB_Pk%Redshifts(itf),dl)
       om_m = omega_m(a, omm0, CP%omegav, w_lam)  
       om_v = omega_v(a, omm0, CP%omegav, w_lam)

      xlogr1=-2.0
      xlogr2=3.5
      do
          rmid=(xlogr2+xlogr1)/2.0
          rmid=10**rmid
          call wint(CAMB_Pk,itf,rmid,sig,d1,d2)
          diff=sig-1.0
          if (abs(diff).le.0.001) then
             rknl=1./rmid
             rneff=-3-d1
             rncur=-d2                  
             exit
          elseif (diff.gt.0.001) then
             xlogr1=log10(rmid)
          elseif (diff.lt.-0.001) then
             xlogr2=log10(rmid)
          endif
          if (xlogr2 < -1.9999) then
               !is still linear, exit
               goto 101
         end if
      end do

! now calculate power spectra for a logarithmic range of wavenumbers (rk)

      do i=1, CAMB_PK%num_k
         rk = exp(CAMB_Pk%log_kh(i))

         if (rk > Min_kh_nonlinear) then

    ! linear power spectrum !! Remeber => plin = k^3 * P(k) * constant
    ! constant = 4*pi*V/(2*pi)^3 

             plin= MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3/(2*pi**2)) 

    ! calculate nonlinear power according to halofit: pnl = pq + ph,
    ! where pq represents the quasi-linear (halo-halo) power and 
    ! where ph is represents the self-correlation halo term. 
 
             call halofit(rk,rneff,rncur,rknl,plin,pnl,pq,ph)   ! halo fitting formula 
             CAMB_Pk%nonlin_ratio(i,itf) = sqrt(pnl/plin)

         end if

      enddo

101   continue
      end do
            
      end subroutine NonLinear_GetNonLinRatios
       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! halo model nonlinear fitting formula as described in 
! Appendix C of Smith et al. (2002)

      subroutine halofit(rk,rn,rncur,rknl,plin,pnl,pq,ph)
      implicit none

      real(dl) gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
      real(dl) rk,rn,plin,pnl,pq,ph,plinaa
      real(dl) rknl,y,rncur
      real(dl) f1a,f2a,f3a,f1b,f2b,f3b,frac

!RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
!than latest N-body simulations at small scales.
!Update the following fitting parameters of gam,a,b,c,xmu,xnu,
!alpha & beta from the simulations in Takahashi+ 2012.
!The improved halofit accurately provide the power spectra for WMAP
!cosmological models with constant w.
      gam=0.1971-0.0843*rn+0.8460*rncur
      a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+ &
           0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*om_v*(1.+w_lam)
      a=10**a      
      b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+ &
           0.2279*om_v*(1.+w_lam))
      c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
      xmu=0.
      xnu=10**(5.2105+3.6902*rn)
      alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
      beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ &
           0.3980*rn**4-0.1682*rncur

      if(abs(1-om_m).gt.0.01) then ! omega evolution 
         f1a=om_m**(-0.0732)
         f2a=om_m**(-0.1423)
         f3a=om_m**(0.0725)
         f1b=om_m**(-0.0307)
         f2b=om_m**(-0.0585)
         f3b=om_m**(0.0743)       
         frac=om_v/(1.-om_m) 
         f1=frac*f1b + (1-frac)*f1a
         f2=frac*f2b + (1-frac)*f2a
         f3=frac*f3b + (1-frac)*f3a
      else         
         f1=1.0
         f2=1.
         f3=1.
      endif

      y=(rk/rknl)

      ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
      ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fnu*(2.080-12.39*(omm0-0.3))/(1+1.201e-03*y**3))
      plinaa=plin*(1+fnu*26.29*rk**2/(1+1.5*rk**2))
      pq=plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)

      pnl=pq+ph

      end subroutine halofit       


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! The subroutine wint, finds the effective spectral quantities
! rknl, rneff & rncur. This it does by calculating the radius of 
! the Gaussian filter at which the variance is unity = rknl.
! rneff is defined as the first derivative of the variance, calculated 
! at the nonlinear wavenumber and similarly the rncur is the second
! derivative at the nonlinear wavenumber. 

      subroutine wint(CAMB_Pk,itf,r,sig,d1,d2)
      implicit none
      integer, intent(in) :: itf
      type(MatterPowerData) :: CAMB_Pk
      real(dl) sum1,sum2,sum3,t,y,x,w1,w2,w3
      real(dl) x2,rk, fac,r, sig, d1,d2, anorm
      integer i,nint

      nint=3000
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      anorm = 1/(2*pi**2)
      do i=1,nint
         t=(i-0.5_dl)/nint
         y=-1.d0+1.d0/t
         rk=y
         d2=MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3*anorm) 
         x=y*r
         x2=x*x
         w1=exp(-x2)
         w2=2*x2*w1
         w3=4*x2*(1-x2)*w1
         fac=d2/y/t/t
         sum1=sum1+w1*fac
         sum2=sum2+w2*fac
         sum3=sum3+w3*fac
      enddo
      sum1=sum1/nint
      sum2=sum2/nint
      sum3=sum3/nint
      sig=sqrt(sum1)
      d1=-sum2/sum1
      d2=-sum2*sum2/sum1/sum1 - sum3/sum1
      
      end subroutine wint
      
!!BR09 generalize to constant w

      function omega_m(aa,om_m0,om_v0,wval)
      implicit none
      real(dl) omega_m,omega_t,om_m0,om_v0,aa,wval
      omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*((aa)**(-1.0-3.0*wval))+om_m0/aa)
      omega_m=omega_t*om_m0/(om_m0+om_v0*((aa)**(-3.0*wval)))
      end function omega_m

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! evolution of omega lambda with expansion factor

      function omega_v(aa,om_m0,om_v0,wval)      
      implicit none
      real(dl) aa,omega_v,om_m0,om_v0,omega_t,wval
      omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*((aa)**(-1.0-3.0*wval))+om_m0/aa)
      omega_v=omega_t*om_v0*((aa)**(-3.0-3.0*wval))/(om_v0*((aa)**(-3.0-3.0*wval))+om_m0/aa/aa/aa)
      end function omega_v

!!BR09 end generalize to constant w

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

       stop 'Halofit module doesn''t support non-linear velocities'

     end subroutine NonLinear_GetRatios_All

