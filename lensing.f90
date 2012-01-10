!Lensing the C_l using the deflection angle from the computed lensing potential
!power spectrum. 
!lensing_method=1: using an accurate curved-sky correlation function method
!lensing_method=2: using the flat-sky lower order result of astro-ph/9505109 
!                  and astro-ph/9803150 as in CMBFAST 
!lensing_method=3: using inaccurate full sky harmonic method of astro-ph/0001303

!The flat sky result is accurate to about 0.1% in TT, and 0.4% in EE and is
!about a factor of two faster than lensing_method=1.
!lensing_method=3 is only present for comparison and is not recommended in any regime

!Set accurate_BB=T if you want BB accurately by integrating the full angular range
!otherwise it saves a large amount of time by only integrating the small scales
!accute_BB only does *not* include any non-linear corrections or ensure you have
!chosen sufficiently high l_max and k_max, so does not neccessarily give an accurate
!result

!Uses the un-lensed Cls and the computed lensing potential power spectrum.
!Usual values of k_max are fine for all but the lensed BB Cls
!To get the lensed BB accurate around l=1000 you need to go to l_max >2000, and
!higher for higher l. Since this probes small scales in the lensing power spectrum you
!also need to go to higher k_max - for concordance models something like 
!k_eta_max_scalar=10000. At l>1000 you can expect to need higher k_max, and for 
!non-linear evolution to cause a significant error.

!Correlation function routines by AL+AC Nov 2004 with flat-sky borrowings from CMBFAST
!Curved sky results use the method of astro-ph/xxx.

!Full sky harmonic lensing routines by Gayoung Chon and AC.
!Ref: astro-ph/0001303 by W. Hu. 
!For better derivations see also astro-ph/0301064 and astro-ph/0301031
!Adapted for CAMB and optimized by AL.
!Uses f90 version of "J1-RECURSION OF 3J-COEFFICIENTS" by K. Schulten and R.G. Gordon 
!obtainable from the CPC program library (www.cpc.cs.qub.ac.uk).

!March 2006: fixed problem with l_max when generating with tensors (thanks Chad Fendt)

module lensing
use Precision
use ModelParams
use AmlUtils
implicit none
 integer, parameter :: lensing_method_curv_corr=1,lensing_method_flat_corr=2, &
                       lensing_method_harmonic=3
 
 integer :: lensing_method = lensing_method_curv_corr

private

 logical  :: lensing_includes_tensors = .false.

!flat method stores
 real(dl), parameter :: dbessel = 0.05_dl 
 real(dl), dimension(:), allocatable :: Bess0, ddBess0
 real(dl), dimension(:), allocatable :: Bess2, ddBess2
 real(dl), dimension(:), allocatable :: Bess4, ddBess4
 real(dl), dimension(:), allocatable :: Bess6, ddBess6

 integer, parameter :: lensed_convolution_margin = 100
   !Number of L less than L max at which the lensed power spectrum is calculated 

!Harmonic method stores
 integer :: lmax_donelnfa = 0
 real(dl), dimension(:), allocatable  :: lnfa

public lens_Cls, lensing_includes_tensors, lensing_method, lensing_method_flat_corr,&
      lensing_method_curv_corr,lensing_method_harmonic, BessI, bessj0
contains


subroutine lens_Cls
 use lvalues

 !Must set l again in case computed tessors (thanks to Chad)
 call initlval(lSamp,CP%Max_l)
 if (lensing_method == lensing_method_curv_corr) then
    call CorrFuncFullSky()
  elseif (lensing_method == lensing_method_flat_corr) then 
    call CorrFuncFlatSky()
  elseif (lensing_method == lensing_method_harmonic) then 
    call BadHarmonic
  else
    stop 'Unknown lensing method'
 end if
end subroutine lens_Cls


subroutine CorrFuncFullSky

  integer :: lmax_extrap 
  
  lmax_extrap = CP%Max_l - lensed_convolution_margin + 450  
  if (HighAccuracyDefault) lmax_extrap=lmax_extrap+300
  lmax_extrap = min(lmax_extrap_highl,lmax_extrap)
  call CorrFuncFullSkyImpl(max(lmax_extrap,CP%max_l))

end subroutine CorrFuncFullSky


subroutine CorrFuncFullSkyImpl(lmax)
 !Accurate curved sky correlation function method
 !Uses non-perturbative isotropic term with 2nd order expansion in C_{gl,2}
 !Neglects C_{gl}(theta) terms (very good approx)
  use ModelParams
  use ModelData
  use lvalues
  implicit none
  integer, intent(in) :: lmax
  integer l, i, in
  integer :: npoints 
  real(dl) corr(4), Cg2, sigmasq, theta
  real(dl) dtheta
  real(dl) llp1,fac, fac1,fac2,fac3, rootllp1, rootfac1, rootfac2, rootfac3
  integer max_lensed_ix
  real(dl) P(lmax),dP(lmax)
  real(dl) sinth,halfsinth, x, T2,T4
  real(dl) roots(-1:lmax+4), lfacs(lmax), lfacs2(lmax), lrootfacs(lmax)
  real(dl) d_11(lmax),d_m11(lmax)
  real(dl) d_22(lmax),d_2m2(lmax),d_20(lmax)
  real(dl) Cphil3(lmin:lmax), CTT(lmin:lmax), CTE(lmin:lmax),CEE(lmin:lmax)
  real(dl) ls(lmax)
  real(dl) xl(lmax),ddcontribs(lmax,4),corrcontribs(lmax,4)
  real(dl), allocatable, dimension(:,:,:) :: lens_contrib(:,:,:)
  integer thread_ix
  real(dl) pmm, pmmp1
  real(dl) d4m4,d11,dm11,d2m2,d22,d20,d23,d2m3,d33,d3m3,d04,d1m2,d12,d13,d1m3,d2m4
  real(dl) sinfac, Cg2sq
  real(dl) X000,X022,X220,X121,X132,X242
  real(dl) dX000,dX022
  real(sp) timeprev
  integer  interp_fac
  integer j,jmax
  integer llo, lhi
  real(dl) a0,b0,ho, sc
  logical :: short_integral_range
  integer, parameter :: slow_highL = 5000 !Lmax at which to do full range to prevent ringing etc

  logical, parameter :: approx = .false.

!$ integer  OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
!$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

    if (lensing_includes_tensors) stop 'Haven''t implemented tensor lensing'

    max_lensed_ix = lSamp%l0-1
!    do while(lSamp%l(max_lensed_ix) > CP%Max_l -250)
    do while(lSamp%l(max_lensed_ix) > CP%Max_l - lensed_convolution_margin) 
      max_lensed_ix = max_lensed_ix -1
    end do
    lmax_lensed = lSamp%l(max_lensed_ix)
    if (allocated(Cl_lensed)) deallocate(Cl_lensed)
    allocate(Cl_lensed(lmin:lmax_lensed,CP%InitPower%nn,1:4))
    
    Cl_Lensed = 0
   
    npoints = CP%Max_l  * 2    
    short_integral_range = .not. CP%AccurateBB .and. CP%Max_l<=slow_highL
    if (.not. short_integral_range ) npoints = npoints * 2 

    dtheta = pi / npoints
    if (short_integral_range) then
      npoints = int(npoints /32 *min(32._dl,AccuracyBoost)) 
      !OK for TT, EE, TE but inaccurate for low l BB
      !this induces high frequency ringing on very small scales
    end if

    if (DebugMsgs) timeprev=GetTestTime()

    if (.not. short_integral_range) then
     !There is an odd serious problem with interpolating if you do a large
     !angular range.
      
     interp_fac =1
    else
     interp_fac = max(1,nint(10/AccuracyBoost))
    end if

    jmax = 0
    do l=lmin,lmax
       if (l<=15 .or. mod(l-15,interp_fac)==interp_fac/2) then
         jmax =jmax+1
         ls(jmax)=l
         xl(jmax)=l
       end if
       lfacs(l) = real(l*(l+1),dl)
       lfacs2(l) = real((l+2)*(l-1),dl)
       lrootfacs(l) = sqrt(lfacs(l)*lfacs2(l))
    end do

    roots(-1)=0 !just so dipole doesn't screw up
    do l=0,lmax+4
     roots(l) = sqrt(real(l,dl))
    end do


    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()  
    allocate(lens_contrib(4,lmax_lensed,thread_ix))

    do in = 1, CP%InitPower%nn

    do l=lmin,CP%Max_l
     ! (2*l+1)l(l+1)/4pi C_phi_phi: Cl_scalar(l,1,C_Phi) is l^4 C_phi_phi
       Cphil3(l) = Cl_scalar(l,in,C_Phi)*(2*l+1)*(l+1)/real(l,dl)**3/(4*pi) 
       fac = (2*l+1)/(4*pi) * 2*pi/(l*(l+1))
       CTT(l) =  Cl_scalar(l,in,C_Temp)*fac
       CEE(l) =  Cl_scalar(l,in,C_E)*fac
       CTE(l) =  Cl_scalar(l,in,C_Cross)*fac
    end do
    if (Cphil3(10) > 1e-7) then
     write (*,*) 'You need to normalize realistically to use lensing.'
     write (*,*) 'see http://cosmocoffee.info/viewtopic.php?t=94'
     stop
    end if
    if (lmax > CP%Max_l) then
     l=CP%Max_l
     sc = (2*l+1)/(4*pi) * 2*pi/(l*(l+1))     
     fac2=CTT(CP%Max_l)/(sc*highL_CL_template(CP%Max_l, C_Temp))
     fac=Cphil3(CP%Max_l)/(sc*highL_CL_template(CP%Max_l, C_Phi))  
     do l=CP%Max_l+1, lmax
       !Fill in tail from template
       sc = (2*l+1)/(4*pi) * 2*pi/(l*(l+1))  
       Cphil3(l) = highL_CL_template(l, C_Phi)*fac*sc
       
       CTT(l) =  highL_CL_template(l, C_Temp)*fac2*sc
       CEE(l) =  highL_CL_template(l, C_E)*fac2 *sc
       CTE(l) =  highL_CL_template(l, C_Cross)*fac2*sc 
      if (Cphil3(CP%Max_l+1) > 1e-7) then
       write (*,*) 'You need to normalize the high-L template so it is dimensionless'
       stop
      end if
     end do
   end if
  lens_contrib=0

  !uncomment second line for PGF90 workaround
  !$OMP PARALLEL DO DEFAULT(PRIVATE),  &
  !OMP PRIVATE(P,dP,d11,dm11,d22,d2m2,d20,corrcontribs,ddcontribs),& 
  !$OMP SHARED(lfacs,lfacs2,lrootfacs,Cphil3,CTT,CTE,CEE,lens_contrib, lmax), &
  !$OMP SHARED(dtheta,CP,lmax_lensed,roots, npoints,interp_fac,jmax,ls,xl,short_integral_range) 
      do i=1,npoints-1

      theta = i * dtheta 
      x = cos(theta)
      sinth = sin(theta)
      halfsinth = sinth/2

      pmm=1
      pmmp1=x

      Cg2=0
      sigmasq=0
      if (lmin==1) then
        d_11(1) = cos(theta/2)**2
        d_m11(1) = sin(theta/2)**2
        sigmasq = sigmasq  +  (1-d_11(1))*Cphil3(lmin) 
        Cg2 = Cg2  + d_m11(1)*Cphil3(lmin)
        P(1) = x
        d_22(1)=0
        d_2m2(1)=0
        d_20(1)=0
      end if
      do l=2,lmax

        P(l)= ((2*l-1)* x *pmmp1 - (l-1)*Pmm)/ l
        dP(l) = l*(pmmp1-x*P(l))/sinth**2
        Pmm=pmmp1
        pmmp1=P(l)
        llp1 = lfacs(l)
  
        fac1 = (1-x)
        fac2 = (1+x)
        fac = fac1/fac2

        d_11(l) =  fac1*dP(l)/llp1 + P(l)
        d_m11(l) = fac2*dP(l)/llp1 - P(l)

        sigmasq = sigmasq  +  (1-d_11(l))*Cphil3(l) 
        Cg2 = Cg2  + d_m11(l)*Cphil3(l)
        
        d_22(l) = ( ((4*x-8)/fac2 + llp1)*P(l) &
            + 4*fac*( fac2 + (x - 2)/llp1)*dP(l) )/ lfacs2(l)
                  
        d_2m2(l) = ( (llp1- (4*x+8)/fac1) *P(l) &
            +4/fac*( -fac1 + (x+2)/llp1) *dP(l) )/lfacs2(l)              

        d_20(l) = (2*x*dP(l) - llp1*P(l) ) / lrootfacs(l)

      end do
    
       do j=1,jmax
        l =ls(j)
 
        fac1 = (1-x)
        fac2 = (1+x)
        llp1 = lfacs(l)
  
        rootllp1 = roots(l)*roots(l+1) 
        rootfac1 = roots(l+2)*roots(l-1)
        rootfac2 = roots(l+3)*roots(l-2)

        llp1=lfacs(l)
        dm11=d_m11(l)
        d11=d_11(l)
        if (l<2) then
         d2m2=0
         d22=0
         d20=0
         d1m2 = 0
         d12 =  0     
        else
         d2m2=d_2m2(l)
         d22=d_22(l)
         d20=d_20(l)
         d1m2 = sinth/rootfac1*(dP(l) -2/fac1*dm11)
         d12 =  sinth/rootfac1*(dP(l) -2/fac2*d11)
        end if
        if (l<3) then
         d1m3=0
         d2m3=0
         d3m3=0
         d13 =0 
         d23 =0
         d33 =0 
        else
         sinfac=4/sinth
         d1m3 = (-(x+0.5_dl)*d1m2*sinfac - lfacs2(l)*dm11/rootfac1 )/rootfac2
         d2m3 = (-fac2*d2m2*sinfac - rootfac1*d1m2)/rootfac2
         d3m3 = (-(x+1.5_dl)*d2m3*sinfac - rootfac1*d1m3)/rootfac2
         d13  =  ((x-0.5_dl)*d12*sinfac - lfacs2(l)*d11/rootfac1 ) /rootfac2
         d23  = (-fac1*d22*sinfac + rootfac1*d12 ) / rootfac2
         d33  = (-(x-1.5_dl)*d23*sinfac - rootfac1*d13)/rootfac2
        end if 
        if (l<4) then
         d04=0
         d2m4=0
         d4m4=0
         rootfac3=0
        else
         rootfac3=roots(l-3)*roots(l+4)
         d04=( (-llp1 + (18*x**2 + 6)/sinth**2 )*d20  -&
             6*x*lfacs2(l)*dP(l)/lrootfacs(l) ) / (rootfac2*rootfac3)
         d2m4= (-(6*x+4)*d2m3/sinth - rootfac2*d2m2 ) / rootfac3
         d4m4 = (-7/5._dl*(llp1-6)*d2m2 + &
                12/5._dl*( -llp1+(9*x+26)/fac1)*d3m3 ) / (llp1-12)
        end if

       !Non perturbative isotropic integrals
       !these are approx, but extremely good approximations
         X000 = exp(-llp1*sigmasq/4)
         if (approx) then

         X022 = X000  
         X220 = rootllp1**2/4*X000
         X121 = -0.5_dl*rootllp1*X000
         X132 = -0.5_dl*rootllp1*X000
         X242 = 0.25_dl*rootllp1**2*X022 
         
         dX000 = -llp1/4*X000
         dX022 = -llp1/4*X022
        

         else
         X022 = X000*(1+sigmasq)   !exp(-(llp1-4)*sigmasq/4)
         X220 = lrootfacs(l)/4*X000
         X121 = -0.5_dl*rootfac1*X000
         X132 = -0.5_dl*rootfac2*X000
         X242 = 0.25_dl*rootfac2*rootfac3*X022 
         
         dX000 = -llp1/4*X000
         dX022 = (1-llp1/4)*X022
         end if  
!second order
         !TT
         fac1 = dX000**2
         fac3 = X220**2
         Cg2sq = Cg2**2

!Here we drop terms in Cgt which are down by powers of l
!Approx good to 1e-4 level
         fac = ( (X000**2-1) + Cg2sq*fac1)*P(l)+ Cg2sq*fac3*d2m2 &
                    + 8/llp1* fac1*Cg2*dm11 
     
         corrcontribs(j,1)=  CTT(l) * fac 

         fac2=(Cg2*dX022)**2+(X022**2-1)
!Q+U
         fac = 2*Cg2*X121*X132*d13 + fac2*d22 +Cg2sq*X242*X220*d04 

         corrcontribs(j,2)= CEE(l) * fac 

!Q-U 
         fac = ( fac3*P(l) + X242**2*d4m4)*Cg2sq/2 &
              + Cg2*(X121**2*dm11+ X132**2*d3m3) + fac2*d2m2 

         corrcontribs(j,3)= CEE(l) * fac 

!TE
        fac = (X000*X022-1)*d20+ &
          2*dX000*Cg2*(X121*d11 + X132*d1m3)/rootllp1 &
             + Cg2sq*(X220/2*d2m4*X242 +( fac3/2 + dX022*dX000)*d20) 

        corrcontribs(j,4)= CTE(l) * fac 

      end do

do j=1,4
  corr(j) = sum(corrcontribs(1:14,j))+interp_fac*sum(corrcontribs(15:jmax,j))
end do

if (short_integral_range .and. i>npoints-20) &
        corr=corr*exp(-(i-npoints+20)**2/150.0) !taper the end to help prevent ringing

!Interpolate contributions
!Increasing interp_fac and using this seems to be slower than above
if (.false.) then
      if (abs(sum(corrcontribs(1:jmax,1)))>1e-11) print *,i,sum(corrcontribs(1:jmax,1))
      do j=1,4
       call spline(xl,corrcontribs(1,j),jmax,1d30,1d30,ddcontribs(1,j))
      end do 
      corr=0
      llo=1
      do l=lmin,lmax
           if ((l > ls(llo+1)).and.(llo < jmax)) then
              llo=llo+1
           end if
           lhi=llo+1
           ho=ls(lhi)-ls(llo)
           a0=(ls(lhi)-l)/ho
           b0=(l-ls(llo))/ho
           fac1 = ho**2/6
           fac2 = (b0**3-b0)*fac1
           fac1 = (a0**3-a0)*fac1 
  
           corr(1) = Corr(1)+ a0*corrcontribs(llo,1)+ b0*corrcontribs(lhi,1)+ &
            fac1* ddcontribs(llo,1) +fac2*ddcontribs(lhi,1)
           corr(2) = Corr(2)+ a0*corrcontribs(llo,2)+ b0*corrcontribs(lhi,2)+ &
            fac1* ddcontribs(llo,2) +fac2*ddcontribs(lhi,2)
           corr(3) = Corr(3)+ a0*corrcontribs(llo,3)+ b0*corrcontribs(lhi,3)+ &
            fac1* ddcontribs(llo,3) +fac2*ddcontribs(lhi,3)
           corr(4) = Corr(4)+ a0*corrcontribs(llo,4)+ b0*corrcontribs(lhi,4)+ &
            fac1* ddcontribs(llo,4) +fac2*ddcontribs(lhi,4)
         
      end do
end if 
      
 !$   thread_ix = OMP_GET_THREAD_NUM()+1

      do l=lmin, lmax_lensed
       !theta factors were put in earlier (already in corr)


       lens_contrib(C_Temp, l, thread_ix)= lens_contrib(C_Temp,l, thread_ix) + &
                                          corr(1)*P(l)*sinth 

       T2 = corr(2)* d_22(l)
       T4 = corr(3)* d_2m2(l)


       lens_contrib(CT_E, l, thread_ix)= lens_contrib(CT_E,l, thread_ix) + &
                                          (T2+T4)*halfsinth 
       lens_contrib(CT_B, l, thread_ix)= lens_contrib(CT_B,l, thread_ix) + &
                                          (T2-T4)*halfsinth 
 
       lens_contrib(CT_Cross, l, thread_ix)= lens_contrib(CT_Cross,l, thread_ix) + &
                                          corr(4)*d_20(l)*sinth 
 
      end do

     end do
  !$OMP END PARALLEL DO
     
      do l=lmin, lmax_lensed
         !sign from d(cos theta) = -sin theta dtheta
       fac = l*(l+1)/OutputDenominator*dtheta *2*pi
       Cl_lensed(l,in,CT_Temp) = sum(lens_contrib(CT_Temp,l,:))*fac &
                 + Cl_scalar(l,in,C_Temp) 
       Cl_lensed(l,in,CT_E) = sum(lens_contrib(CT_E,l,:))*fac &
                 + Cl_scalar(l,in,C_E) 
       Cl_lensed(l,in,CT_B) = sum(lens_contrib(CT_B,l,:))*fac
       Cl_lensed(l,in,CT_Cross) = sum(lens_contrib(CT_Cross,l,:))*fac &
                 + Cl_scalar(l,in,C_Cross) 

      end do

      end do !loop over different initial power spectra
     deallocate(lens_contrib)

     if (DebugMsgs) write(*,*) GetTestTime()-timeprev, 'Time for corr lensing'

end subroutine CorrFuncFullSkyImpl



subroutine CorrFuncFlatSky
 !Do flat sky approx partially non-perturbative lensing, lensing_method=2
   use ModelParams
  use ModelData
  use lvalues
  integer l, i
  integer :: npoints 
  real(dl) Cgl2,  sigmasq, theta
  real(dl) dtheta
  real(dl) dbessfac, fac, fac1,fac2,  C2term, expsig, corr(4)
  real(sp) timeprev
  real(dl) Bessel0(lmin:CP%Max_l),Bessel2(lmin:CP%Max_l)
  real(dl) Bessel4(lmin:CP%Max_l),Bessel6(lmin:CP%Max_l)
  real(dl) Cphil3(lmin:CP%Max_l), CTT(lmin:CP%Max_l), CTE(lmin:CP%Max_l),CEE(lmin:CP%Max_l)
  integer max_lensed_ix
  integer b_lo
  integer in
  real(dl) T2,T4,a0, b0
  real(dl) lfacs(CP%Max_l)
  real(dl), allocatable, dimension(:,:,:) :: lens_contrib(:,:,:)
  integer thread_ix
!$ integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
!$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

    if (lensing_includes_tensors) stop 'Haven''t implemented tensor lensing'

    max_lensed_ix = lSamp%l0-1
    do while(lSamp%l(max_lensed_ix) > CP%Max_l -250)
      max_lensed_ix = max_lensed_ix -1
    end do
    lmax_lensed = lSamp%l(max_lensed_ix)
    if (allocated(Cl_lensed)) deallocate(Cl_lensed)
    allocate(Cl_lensed(lmin:lmax_lensed,CP%InitPower%nn,1:4))
    
    Cl_Lensed = 0
   
    npoints = CP%Max_l  * 2   
    if (CP%AccurateBB) npoints = npoints * 2

    dtheta = pi / npoints
    if (.not. CP%AccurateBB) then
     npoints = int(npoints /32 *min(32._dl,AccuracyBoost)) 
      !OK for TT, EE, TE but inaccurate for low l BB
      !this induces high frequency ringing on very small scales
    end if

    call GetBessels(npoints*dtheta*CP%Max_l)

    if (DebugMsgs) timeprev=GetTestTime()

    dbessfac = dbessel**2/6

    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()  
    allocate(lens_contrib(4,lmax_lensed,thread_ix))

    do in = 1, CP%InitPower%nn

    do l=lmin,CP%Max_l
     ! l^3 C_phi_phi/2/pi: Cl_scalar(l,1,C_Phi) is l^4 C_phi_phi
       Cphil3(l) = Cl_scalar(l,in,C_Phi)/l /(2*pi)
       fac = l/(2*pi)*2*pi/(l*(l+1))
       CTT(l) =  Cl_scalar(l,in,C_Temp)*fac
       CEE(l) =  Cl_scalar(l,in,C_E)*fac
       CTE(l) =  Cl_scalar(l,in,C_Cross)*fac
       lfacs(l) = l**2*0.5_dl
    end do

    if (Cphil3(10) > 1e-7) then
     write (*,*) 'You need to normalize realistically to use lensing.'
     write (*,*) 'see http://cosmocoffee.info/viewtopic.php?t=94'
     stop
    end if

  lens_contrib=0

  !$OMP PARALLEL DO DEFAULT(SHARED),  &
  !$OMP PRIVATE(theta, sigmasq,cgl2,b_lo,a0,b0,fac,fac1,fac2), &
  !$OMP PRIVATE(Bessel0,Bessel2,Bessel4,Bessel6), &
  !$OMP PRIVATE(corr,expsig,C2term,T2,T4,i,l, thread_ix)     

    do i=1,npoints-1

      theta = i * dtheta 
      sigmasq =0
      Cgl2=0
      fac = theta /dbessel
     
      do l=lmin,CP%Max_l

!Interpolate the Bessel functions, and compute sigma^2 and C_{gl,2} 
        b0 = l*fac
        b_lo = int(b0) +1 
        a0=  b_lo - b0                
        b0=  1._dl - a0 
        fac1 = a0*b0*dbessfac
        fac2 = fac1*(a0-2)
        fac1 = fac1*(b0-2)

        Bessel0(l) = a0*Bess0(b_lo)+ b0*Bess0(b_lo+1) +fac1*ddBess0(b_lo) &
                       +fac2*ddBess0(b_lo+1)
        sigmasq = sigmasq + (1-Bessel0(l))*Cphil3(l) 


        Bessel2(l) = a0*Bess2(b_lo)+ b0*Bess2(b_lo+1) +fac1*ddBess2(b_lo) &
                      +fac2*ddBess2(b_lo+1)
        Cgl2 =  Cgl2 + Bessel2(l)*Cphil3(l)

        Bessel4(l) = a0*Bess4(b_lo)+ b0*Bess4(b_lo+1) +fac1*ddBess4(b_lo) &
                      +fac2*ddBess4(b_lo+1)
        Bessel6(l) = a0*Bess6(b_lo)+ b0*Bess6(b_lo+1) +fac1*ddBess6(b_lo) &
                      +fac2*ddBess6(b_lo+1)

      end do

!Get difference between lensed and unlensed correlation function
     corr = 0
      do l=lmin,CP%Max_l
!For 2nd order perturbative result use 
!         expsig = 1 -sigmasq*l**2/2._dl
!         C2term = l**2*Cgl2/2._dl
          fac = sigmasq*lfacs(l)
          expsig = exp(-fac) 
          C2term = Cgl2*lfacs(l)
!Put theta factor later  in here
          fac1 = expsig*theta
          fac2 = C2term*fac1
          fac1 = fac1 - theta  !we want expsig-1 to get lensing difference

          fac = fac1*Bessel0(l) + fac2*Bessel2(l) 

          !TT
          corr(1) = corr(1) + CTT(l) * fac                              

          !Q + U
          corr(2) = corr(2) + CEE(l) * fac                              
          fac2 = fac2*0.5_dl
          !Q-U
          corr(3) = corr(3) + CEE(l) * &
              (fac1*Bessel4(l) + fac2*(Bessel2(l)+Bessel6(l)))                               
          !Cross
          corr(4) = corr(4) + CTE(l) * &
              (fac1*Bessel2(l) + fac2*(Bessel0(l)+Bessel4(l)))                               
 

      end do

      
 !$   thread_ix = OMP_GET_THREAD_NUM()+1

      do l=lmin, lmax_lensed
       !theta factors were put in earlier (already in corr)
       lens_contrib(C_Temp, l, thread_ix)= lens_contrib(C_Temp,l, thread_ix) + &
                                          corr(1)*Bessel0(l) 
       T2 = corr(2)*Bessel0(l)
       T4 = corr(3)*Bessel4(l)
       lens_contrib(CT_E,l,thread_ix)  = lens_contrib(CT_E,l, thread_ix) + T2+T4
       lens_contrib(CT_B,l,thread_ix)  = lens_contrib(CT_B,l, thread_ix) + T2-T4
       lens_contrib(CT_Cross,l, thread_ix) = lens_contrib(CT_Cross,l, thread_ix) + &
                                              corr(4)*Bessel2(l)
      end do

     end do
  !$OMP END PARALLEL DO
     
      do l=lmin, lmax_lensed
       fac = l*(l+1)* 2*pi/OutputDenominator*dtheta
       Cl_lensed(l,in,CT_Temp) = sum(lens_contrib(CT_Temp,l,:))*fac &
                 + Cl_scalar(l,in,CT_Temp) 
       Cl_lensed(l,in,CT_Cross) = sum(lens_contrib(CT_Cross,l,:))*fac &
                 +Cl_scalar(l,in,C_Cross)
       fac = fac /2 !(factor of 1/2 should have been in T2+/-T4 above           
       Cl_lensed(l,in,CT_E) = sum(lens_contrib(CT_E,l,:))*fac &
                 + Cl_scalar(l,in,CT_E) 
       Cl_lensed(l,in,CT_B) = sum(lens_contrib(CT_B,l,:))*fac
      end do

      end do !loop over different initial power spectra
     deallocate(lens_contrib)

     if (DebugMsgs) write(*,*) GetTestTime()-timeprev, 'Time for corr lensing'

end subroutine CorrFuncFlatSky

subroutine BadHarmonic
  use ModelParams
  use ModelData
  use lvalues
  use InitialPower
  integer maxl, i, in, almin, max_lensed_ix, maxl_phi
  real(dl) , dimension (:,:,:), allocatable :: bare_cls
  real(dl) pp(CP%InitPower%nn,CP%Max_l)
  real(dl) asum(CP%InitPower%nn), RR(CP%InitPower%nn), roots(CP%Max_l)
  real(dl) asum_TE(CP%InitPower%nn), asum_EE(CP%InitPower%nn), asum_BB(CP%InitPower%nn)
  integer l1,l2,al,j, j1, k, hk, llp_1, llp_al, g1
  real(dl)  F, fct
  real(dl) g2l,g2l1, norm
  real(dl) a3j(CP%Max_l*2+1), tF, expF
  logical DoPol
  real(dl) iContribs(lSamp%l0,CP%InitPower%nn, 1:4), intcontrib(lmin:lSamp%l(lSamp%l0))
  real(dl) , dimension (:,:,:), allocatable :: iCl_lensed
  integer max_j_contribs

  real(sp) timeprev
    
!Otherwise use second order perturbative harmonic method

  if (DebugMsgs) timeprev=GetTestTime()

  DoPol = CP%AccuratePolarization

  maxl = CP%Max_l
 
  if (allocated(Cl_lensed)) deallocate(Cl_lensed)


  allocate(bare_cls(CP%InitPower%nn,maxl,1:4))
  
  RR = 0
  do j=lmin,maxl
     norm = OutputDenominator/(j*(j+1))
     if (lensing_includes_tensors .and. CP%WantTensors .and. j<= CP%Max_l_tensor) then !Use total Cls
      bare_cls(:,j,CT_Temp:CT_E) = (Cl_scalar(j,:,C_Temp:C_E) + &
           Cl_tensor(j,:,CT_Temp:CT_E))*norm
      bare_cls(:,j,CT_B) = Cl_tensor(j,:,CT_B)*norm
      bare_cls(:,j,CT_Cross) =  (Cl_scalar(j,:,C_Cross) + &
          Cl_tensor(j,:,CT_Cross))*norm
     else
      bare_cls(:,j,CT_Temp:CT_E) = Cl_scalar(j,:,C_Temp:C_E)*norm
      bare_cls(:,j,CT_B) = 0
      bare_cls(:,j,CT_Cross) =  Cl_scalar(j,:,C_Cross)*norm
     end if
     pp(:,j) = Cl_scalar(j,:,C_Phi)/real(j**2,dl)**2
     RR = RR + j*(j+1)*real(2*j+1,dl)*pp(:,j)
     roots(j) = sqrt(real(2*j+1,dl))
  end do

  RR = RR/2/fourpi
  if (RR(1) > 1e-5) then
     write (*,*) 'You need to normalize realistically to use lensing.'
     write (*,*) 'see http://cosmocoffee.info/viewtopic.php?t=94'
   stop
  end if
  if (maxl > lmax_donelnfa) then 
   !Get ln factorials
   if (allocated(lnfa)) deallocate(lnfa)
   allocate(lnfa(0:maxl*3+1))
   lmax_donelnfa = maxl 
   lnfa(0) = 0
   do i=1,CP%Max_l*3+1 
     lnfa(i)=lnfa(i-1) + log(real(i,dl))
   end do
  end if
   
  max_lensed_ix = lSamp%l0-1
  do while(lSamp%l(max_lensed_ix) > maxl -250)
     max_lensed_ix = max_lensed_ix -1
  end do
  lmax_lensed = lSamp%l(max_lensed_ix)

  allocate(iCl_lensed(max_lensed_ix, CP%InitPower%nn, 1:4))

  max_j_contribs = lSamp%l0-1
  if (.not. DoPol) then
           maxl_phi = min(maxl,nint(max(600,(maxl*2)/5)*scale*AccuracyBoost))
           do while (lSamp%l(max_j_contribs) > maxl_phi)
              max_j_contribs=max_j_contribs-1
           end do
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC), SHARED(max_j_contribs) &
  !$OMP PRIVATE(al,g1,llp_al,llp_1,g2l,asum,l1,g2l1,l2,k,hk,F,fct,almin), &
  !$OMP PRIVATE(asum_EE,asum_BB,asum_TE,expF,tF, a3j, iContribs,in,intcontrib)
  do j=max_lensed_ix,1,-1  
     !Only compute lensed spectra at lSamp%l(j). Start with slow ones.
     
     al=lSamp%l(j)

     llp_al = al*(al+1)
     g2l=sqrt((2*al+1)/fourpi)

     asum = 0
     asum_EE = 0
     asum_BB = 0
     asum_TE = 0

   
     do j1 = 1, max_j_contribs
        !  Contributions to C_al are a smooth function of l_1 - so interpolate
        l1=lSamp%l(j1)

        llp_1 = l1*(l1+1)
        g2l1=roots(l1)

        almin = max(abs(al-l1),2)

        if (DoPol) then
          call GetThreeJs(a3j(almin),l1,al,0,2)
          do l2= almin, min(maxl,al+l1)
              g1 = llp_1+l2*(l2+1)-llp_al
              if (g1 == 0 ) cycle

              k=al+l1+l2
              fct=g1*g2l*g2l1*roots(l2)/2
              tF = fct*a3j(l2)

              if (mod(k,2)==0) then

                 hk = k/2
                 F = lnfa(hk)-lnfa(hk-al)-lnfa(hk-l1)-lnfa(hk-l2)+(lnfa(k-2*al)+lnfa(k-2*l1)&
                  & +lnfa(k-2*l2)-lnfa(k+1))/2
                 
                 expF = exp(F)

                 asum=asum + bare_cls(:,l2,C_Temp)*(expF*fct)**2
        
                 asum_EE = asum_EE + bare_cls(:,l2,CT_E)*tF**2
                 asum_BB = asum_BB + bare_cls(:,l2,CT_B)*tF**2
                 if (mod(hk,2)/=0) tF=-tF
                 asum_TE = asum_TE + bare_cls(:,l2,CT_Cross)*expF*fct*tF

              else
                
                 asum_BB = asum_BB + bare_cls(:,l2,CT_E)*tF**2
                 asum_EE = asum_EE +bare_cls(:,l2,CT_B)*tF**2 

              end if
              
          end do

        else !No polarization
          do l2= almin +mod(al+l1+almin,2),min(maxl,al+l1), 2 
             !Only do lSamp%l's where al + l1 + l2 is even
             
              g1 = llp_1+l2*(l2+1)-llp_al
            
              if (g1 == 0 ) cycle  !Contribution is zero

              k=al+l1+l2
              hk=k/2
           
              fct=g1*g2l*g2l1*roots(l2)/2
              expF = exp(2*(lnfa(hk)-lnfa(hk-al)-lnfa(hk-l1)-lnfa(hk-l2))+lnfa(k-2*al)+lnfa(k-2*l1)&
                  & +lnfa(k-2*l2)-lnfa(k+1))
              asum=asum + bare_cls(:,l2,CT_Temp)*expF *fct**2
          
          end do
          end if !No polarization


             iContribs(j1,:,CT_Temp) = asum*pp(:,l1)
             if (DoPol) then
                iContribs(j1,:,CT_E) = asum_EE*pp(:,l1)
                iContribs(j1,:,CT_B) = asum_BB*pp(:,l1)
                iContribs(j1,:,CT_Cross) = asum_TE*pp(:,l1)
             end if
             asum = 0
             asum_EE = 0
             asum_BB = 0
             asum_TE = 0

        
       end do

          
       !Interpolate contributions to sum and add up
         do in=1, CP%InitPower%nn
          
            call InterpolateClArr(lSamp,iContribs(1,in,CT_Temp),intcontrib,max_j_contribs)
            asum(in) = sum(intcontrib(lmin:lSamp%l(max_j_contribs)))
            if (DoPol) then
               call InterpolateClArr(lSamp,iContribs(1,in,CT_E),intcontrib,max_j_contribs)
               asum_EE(in) = sum(intcontrib(lmin:lSamp%l(max_j_contribs)))
               call InterpolateClArr(lSamp,iContribs(1,in,CT_B),intcontrib,max_j_contribs)
               asum_BB(in) = sum(intcontrib(lmin:lSamp%l(max_j_contribs)))
               call InterpolateClArr(lSamp,iContribs(1,in,CT_Cross),intcontrib,max_j_contribs)
               asum_TE(in) = sum(intcontrib(lmin:lSamp%l(max_j_contribs)))
            end if
         end do

     iCl_lensed(j,:,CT_Temp) =  ((1-al*(al+1)*RR)*bare_cls(:,al,CT_Temp)  & !Linear part
              + asum/(2*al+1))*llp_al/OutputDenominator !add quadratic part and *l(l+1)/2pi
     if (DoPol) then
        iCl_lensed(j,:,CT_E) = ((1-(al**2+al-4)*RR)*bare_cls(:,al,CT_E)  & 
              + asum_EE/(2*al+1))*llp_al/OutputDenominator
        iCl_lensed(j,:,CT_B) = ((1-(al**2+al-4)*RR)*bare_cls(:,al,CT_B)  & 
              + asum_BB/(2*al+1))*llp_al/OutputDenominator
        iCl_lensed(j,:,CT_Cross) =  ((1-(al**2+al-2)*RR)*bare_cls(:,al,CT_Cross) &
                + asum_TE/(2*al+1))*llp_al/OutputDenominator

     else
        iCl_lensed(j,:,CT_E:CT_Cross) = bare_cls(:,al,CT_E:CT_Cross)
     end if

  end do
  !$OMP END PARALLEL DO

  deallocate(bare_cls)

  allocate(Cl_lensed(lmin:lmax_lensed,CP%InitPower%nn,1:4))

  !Interpolate to get final spectrum
  do in=1, CP%InitPower%nn
     do j = CT_Temp, CT_Cross
      call InterpolateClArr(lSamp,iCl_lensed(1,in,j),Cl_lensed(lmin, in, j),max_lensed_ix)
     end do
  end do

  deallocate(iCl_lensed)

  if (DebugMsgs) then
        if (FeedbackLevel>0) write(*,*) GetTestTime()-timeprev,' Timing for lensing'
   end if

end subroutine BadHarmonic

      subroutine GetBessels(MaxArg)
       real(dl), intent(in):: MaxArg
       integer i
       real(dl), allocatable, dimension(:) :: x
       integer max_bes_ix
       integer, save :: last_max = 0

       max_bes_ix = nint(MaxArg / dbessel) + 3
       if (max_bes_ix > last_max) then
           last_max = max_bes_ix
           if (allocated(Bess0)) then
             deallocate(Bess0,ddBess0)
             deallocate(Bess2,ddBess2)
             deallocate(Bess4,ddBess4)
             deallocate(Bess6,ddBess6)
           end if
           allocate(Bess0(max_bes_ix),ddBess0(max_bes_ix))
           allocate(Bess2(max_bes_ix),ddBess2(max_bes_ix))
           allocate(Bess4(max_bes_ix),ddBess4(max_bes_ix))
           allocate(Bess6(max_bes_ix),ddBess6(max_bes_ix))

           allocate(x(max_bes_ix))
           Bess0(1)=1
           Bess2(1)=0; Bess4(1)=0; Bess6(1)=0
           x(1)=0
           do i=2, max_bes_ix
             x(i) = (i-1)*dbessel
             Bess0(i) = Bessj0(x(i)) 
             Bess2(i) = Bessj(2,x(i)) 
             Bess4(i) = Bessj(4,x(i)) 
             Bess6(i) = Bessj(6,x(i)) 
           end do 
           call spline(x,Bess0,max_bes_ix,spl_large,spl_large,ddBess0)
           call spline(x,Bess2,max_bes_ix,spl_large,spl_large,ddBess2)
           call spline(x,Bess4,max_bes_ix,spl_large,spl_large,ddBess4)
           call spline(x,Bess6,max_bes_ix,spl_large,spl_large,ddBess6)

           deallocate(x)
       end if

      end subroutine GetBessels



      FUNCTION bessj0(x)
      real(dl) bessj0,x
      real(dl) ax,xx,z
      real(dl) p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
        s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4, &
       s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
      -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
      .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
       651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2, &
       s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0, &
       59272.64853d0,267.8532712d0,1.d0/

      if(abs(x).lt.8.d0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
          (s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
           p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
!  (C) Copr. 1986-92 Numerical Recipes Software

      END FUNCTION bessj0



      FUNCTION BESSJ1(X)
      real(dl), intent(in) :: x
      real(dl) bessj1,ax,z,xx
      real(dl) Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5, &
         R6,S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0,&
         -2972611.439D0,15704.48260D0,-30.16036606D0/, &
         S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0, &
         18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5, & 
         -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3, &     
         .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))) &
           /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8.0d0/AX
        Y=Z**2
        XX=AX-2.356194491d0
        BESSJ1=SQRT(.636619772d0/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y &
           *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5))))) &
           *SIGN(1._dl,x)
      ENDIF

      END FUNCTION BESSJ1


      FUNCTION BESSJ(N,X)
      real(dl) bessj
      real(dl), intent(in) :: x
      integer, intent(in) :: n
      integer, parameter :: IACC = 40
      real(dl), parameter :: BIGNO=1.d10,BIGNI=1.d-10
      integer jsum,j,m
      real(dl) bj,bjm, bjp, tox, sum
      
      IF(N.LT.2)STOP 'bad argument N in BESSJ'

      TOX=2/X
      IF(X.GT.FLOAT(N))THEN
        BJM=BESSJ0(X)
        BJ=BESSJ1(X)
        DO J=1,N-1
          BJP=J*TOX*BJ-BJM
          BJM=BJ
          BJ=BJP
        END DO
        BESSJ=BJ
      ELSE
        M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
        BESSJ=0.0d0
        JSUM=0
        SUM=0._dl
        BJP=0._dl
        BJ=1._dl
        DO J=M,1,-1
          BJM=J*TOX*BJ-BJP
          BJP=BJ
          BJ=BJM
          IF(ABS(BJ).GT.BIGNO)THEN
            BJ=BJ*BIGNI
            BJP=BJP*BIGNI
            BESSJ=BESSJ*BIGNI
            SUM=SUM*BIGNI
          ENDIF
          IF(JSUM.NE.0)SUM=SUM+BJ
          JSUM=1-JSUM
          IF(J.EQ.N)BESSJ=BJP
        end do
        SUM=2.*SUM-BJ
        BESSJ=BESSJ/SUM
      ENDIF
      END FUNCTION BESSJ

! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      double precision X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END FUNCTION BESSI0
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      double precision X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END FUNCTION BESSI1


      FUNCTION BESSI(N,X)
      !from http://perso.orange.fr/jean-pierre.moreau/Fortran/tbessi_f90.txt
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
      integer, intent(in) :: N
      integer, PARAMETER :: IACC = 40
      integer m,j
      double precision, parameter ::  BIGNO = 1.D10, BIGNI = 1.D-10
      double precision X,BESSI,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO J = M,1,-1
      BIM = BIP+ J*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
      END DO
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END FUNCTION BESSI


end module lensing

