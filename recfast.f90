!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C Integrator for Cosmic Recombination of Hydrogen and Helium,
!C developed by Douglas Scott (dscott@astro.ubc.ca)
!C based on calculations in the paper Seager, Sasselov & Scott
!C (ApJ, 523, L1, 1999).
!C
!C Permission to use, copy, modify and distribute without fee or royalty at
!C any tier, this software and its documentation, for any purpose and without
!C fee or royalty is hereby granted, provided that you agree to comply with
!C the following copyright notice and statements, including the disclaimer,
!C and that the same appear on ALL copies of the software and documentation,
!C including modifications that you make for internal use or for distribution:
!C
!C Copyright 1999 by University of British Columbia.  All rights reserved.
!C
!C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
!C REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  
!C BY WAY OF EXAMPLE, BUT NOT LIMITATION,
!c U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF 
!C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
!C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE 
!C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   
!C
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!CN     Name:        RECFAST
!CV     Version: 1.0
!C 
!CP     Purpose:  Calculate ionised fraction as a function of redshift.
!CP            Solves for H and He simultaneously, and includes
!CP            "fudge factor" for low z effect.
!C
!CD     Description: Solves for ionisation history since recombination
!CD     using the equations in Seager, Sasselov & Scott (ApJ, 1999).
!CD     The Cosmological model can be flat or open.
!CD     The matter temperature is also followed.
!CD     The values for \alpha_B for H are from Hummer (1994).
!CD     The singlet HeI coefficient is a fit from the full code.
!CD     Care is taken to use the most accurate constants.
!C            
!CA     Arguments:
!CA     Name, Description
!CA     real(dl) throughout
!CA
!CA     z is redshift - W is sqrt(1+z), like conformal time
!CA     x is total ionised fraction, relative to H
!CA     x_H is ionized fraction of H - y(1) in R-K routine
!CA     x_He is ionized fraction of He - y(2) in R-K routine
!CA     Tmat is matter temperature - y(3) in R-K routine
!CA     f's are the derivatives of the Y's
!CA     alphaB is case B recombination rate
!CA     alpHe is the singlet only HeII recombination rate
!CA     a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     a_VF is Verner and Ferland type fitting parameter for Helium
!CA     b_VF is Verner and Ferland type fitting parameter for Helium
!CA     T_0 is Verner and Ferland type fitting parameter for Helium
!CA     T_1 is Verner and Ferland type fitting parameter for Helium
!CA     Tnow is the observed CMB temperature today
!CA     Yp is the primordial helium abundace
!CA     fHe is He/H number ratio = Yp/4(1-Yp)
!CA     Trad and Tmat are radiation and matter temperatures
!CA     OmegaB is Omega in baryons today
!CA     H is Hubble constant in units of 100 km/s/Mpc
!CA     HO is Hubble constant in SI units
!CA     bigH is 100 km/s/Mpc in SI units
!CA     G is grvitational constant
!CA     n is number density of hydrogen
!CA     Nnow is number density today
!CA     x0 is initial ionized fraction
!CA     x_H0 is initial ionized fraction of Hydrogen
!CA     x_He0 is initial ionized fraction of Helium
!CA     rhs is dummy for calculating x0
!CA     zinitial and zfinal are starting and ending redshifts
!CA     zeq is the redshift of matter-radiation equality
!CA     zstart and zend are for each pass to the integrator
!CA     w0 and w1 are conformal-time-like initial and final zi and zf's
!CA     Lw0 and Lw1 are logs of w0 and w1
!CA     hw is the interval in W
!CA     C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
!CA     m_e,m_H: electron mass and mass of H atom in SI
!CA     sigma: Thomson cross-section
!CA     a: radiation constant for u=aT^4
!CA     Lambda: 2s-1s two photon rate for Hydrogen
!CA     Lambda_He: 2s-1s two photon rate for Helium
!CA     DeltaB: energy of first excited state from continuum = 3.4eV
!CA     DeltaB_He: energy of first excited state from cont. for He = 3.4eV
!CA     L_H_ion: level for H ionization in m^-1
!CA     L_H_alpha: level for H Ly alpha in m^-1
!CA     L_He1_ion: level for HeI ionization
!CA     L_He2_ion: level for HeII ionization
!CA     L_He_2s: level for HeI 2s
!CA     L_He_2p: level for HeI 2p
!CA     Lalpha: Ly alpha wavelength in SI
!CA     Lalpha_He: Helium I 2p-1s wavelength in SI
!CA     mu_H,mu_T: mass per H atom and mass per particle
!CA     H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
!CA     CDB=DeltaB/k_B                     Constants derived from B1,B2,R
!CA     CDB_He=DeltaB_He/k_B  n=2-infinity for He in Kelvin
!CA     CB1=CDB*4.         Lalpha and sigma_Th, calculated
!CA     CB1_He1: CB1 for HeI ionization potential
!CA     CB1_He2: CB1 for HeII ionization potential
!CA     CR=2*Pi*(m_e/h_P)*(k_B/h_P)  once and passed in a common block
!CA     CK=Lalpha**3/(8.*Pi)
!CA     CK_He=Lalpha_He**3/(8.*Pi)
!CA     CL=C*h_P/(k_B*Lalpha)
!CA     CL_He=C*h_P/(k_B*Lalpha_He)
!CA     CT=(8./3.)*(sigma/(m_e*C))*a
!CA     Bfact=exp((E_2p-E_2s)/kT)    Extra Boltzmann factor
!CA     tol: tolerance for the integrator
!CA     cw(24),w(3,9): work space for DVERK
!CA     Ndim: number of d.e.'s to solve (integer)
!CA     Nz: number of output redshitf (integer)
!CA     I: loop index (integer)
!CA     ind,nw: work-space for DVERK (integer)
!C
!CG     Global data (common blocks) referenced:
!CG     /Cfund/C,k_B,h_P,m_e,m_H,sigma,a,Pi
!CG     /data/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
!CG          fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He
!C
!CF     File & device access:
!CF     Unit /I,IO,O  /Name (if known)
!C
!CM     Modules called:
!CM     DVERK (numerical integrator)
!CM     GET_INIT (initial values for ionization fractions)
!CM     ION (ionization and Temp derivatices)
!C
!CC     Comments:
!CC     none
!C
!CH     History:
!CH     CREATED            (simplest version) 19th March 1989
!CH     RECREATED    11th January 1995
!CH               includes variable Cosmology
!CH               uses DVERK integrator
!CH               initial conditions are Saha
!CH     TESTED              a bunch, well, OK, not really
!CH     MODIFIED     January 1995 (include Hummer's 1994 alpha table)
!CH               January 1995 (include new value for 2s-1s rate)
!CH               January 1995 (expand comments)
!CH               March 1995 (add Saha for Helium)
!CH               August 1997 (add HeII alpha table)
!CH               July 1998 (include OmegaT correction and H fudge factor)
!CH               Nov 1998 (change Trad to Tmat in Rup)
!CH               Jan 1999 (tidied up for public consumption)
!CH               Sept 1999 (switch to formula for alpha's, fix glitch)
!CH                  Sept 1999 modified to CMBFAST by US & MZ          
!CH                     Nov 1999 modified for F90 and CAMB (AML)
!CH                     Aug 2000 modified to prevent overflow erorr in He_Boltz (AML)
!CH                     Feb 2001 corrected fix of Aug 2000 (AML)
!CH                     Oct 2001 fixed error in hubble parameter, now uses global function (AML)
!                       March 2003 fixed bugs reported by savita gahlaut
!                       March 2005 added option for corrections from astro-ph/0501672.
!                                  thanks to V.K.Dubrovich, S.I.Grachev
!                       June 2006 defined RECFAST_fudge as free parameter (AML)
!!      ===============================================================

       module RECDATA
        use Precision
        implicit none
         
        real(dl) C,k_B,h_P,m_e,m_H,sigma,a,G
        real(dl) Lambda,DeltaB,DeltaB_He,Lalpha,mu_H,mu_T,H_frac
        real(dl) Lambda_He,Lalpha_He,Bfact,CK_He,CL_He
        real(dl) L_H_ion,L_H_alpha,L_He1_ion,L_He2_ion,L_He_2s,L_He_2p
        real(dl) CB1,CDB,CR,CK,CL,CT,fHe,CB1_He1,CB1_He2,CDB_He,fu
        real(dl), parameter :: bigH=100.0D3/(1.0D6*3.0856775807D16) !Ho in s-1
        real(dl) Tnow,HO,Nnow

!       --- Data
        data    C,k_B,h_P /2.99792458D8,1.380658D-23,6.6260755D-34/
        data    m_e,m_H     /9.1093897D-31,1.673725D-27/ !av. H atom
        data    sigma,a     /6.6524616D-29,7.565914D-16/ !a=4/c*sigma_B
        data    G   /6.67259D-11/
!       Fundamental constants in SI units

        data    Lambda            /8.2245809d0/
        data    Lambda_He   /51.3d0/   !new value from Dalgarno
        data    L_H_ion          /1.096787737D7/      !level for H ion. (in m^-1)
        data    L_H_alpha   /8.225916453D6/ !averaged over 2 levels
        data    L_He1_ion   /1.98310772D7/     !from Drake (1993)
        data    L_He2_ion   /4.389088863D7/    !from JPhysChemRefData (1987)
        data    L_He_2s          /1.66277434D7/       !from Drake (1993)
        data    L_He_2p          /1.71134891D7/       !from Drake (1993)
!       2 photon rates and atomic levels in SI units
       end module RECDATA


        module RECFAST
        use Precision
        implicit none
        private
        integer, parameter :: Nz0=10000
        real(dl) zrec(Nz0),xrec(Nz0),dxrec(Nz0)
        integer Nz
        real(dl) :: RECFAST_fudge = 1.14

        logical :: use_Dubrovich = .false. !use astro-ph/0501672 corrections
 
        public xeRECFAST, InitRECFAST, use_Dubrovich, RECFAST_fudge
       contains

        function xeRECFAST(a)
        real(dl) a,z,az,bz,xeRECFAST
        integer ilo,ihi
        
        z=1/a-1
        if (z.ge.zrec(1)) then
          xeRECFAST=xrec(1)
        else
         if (z.le.zrec(nz)) then
          xeRECFAST=xrec(nz)
         else
          ilo=nz-z
          ihi=ilo+1
          az=z-int(z)
          bz=1._dl-az        
          xeRECFAST=az*xrec(ilo)+bz*xrec(ihi)+ &
           ((az**3-az)*dxrec(ilo)+(bz**3-bz)*dxrec(ihi))/6._dl
         endif
        endif

        end function xeRECFAST

        subroutine InitRECFAST(OmegaB,h0inp,tcmb,yp)
        use RECDATA
        implicit none
        real(dl), save :: last_OmB =0, Last_YHe=0, Last_H0=0, Last_dtauda=0, last_fudge 

        real(dl) Trad,Tmat,d0hi,d0lo
        integer Ndim,I

        real(dl) OmegaB,H
        real(dl) z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0,h0inp
        real(dl) zinitial,zfinal
        real(dl) zstart,zend,w0,w1,Lw0,Lw1,hw,tcmb
        real(dl) cw(24),w(3,9)
        real(dl) y(3)
        real(dl) yp
        real dum
   
        integer ind,nw

!       --- Parameter statements
        real(dl), parameter :: tol=1.D-5                !Tolerance for R-K

        real(dl) dtauda
        external dtauda

!       ===============================================================

        if (Last_OmB==OmegaB .and. Last_H0 == h0inp .and. yp == Last_YHe .and. & 
             dtauda(0.2352375823_dl) == Last_dtauda .and. last_fudge == RECFAST_fudge) return
           !This takes up most of the single thread time, so cache if at all possible
           !For example if called with different reionization, or tensor rather than scalar
        
        Last_dtauda =  dtauda(0.2352375823_dl) !Just get it at a random scale factor
        Last_OmB = OmegaB
        Last_H0 = h0inp
        Last_YHe=yp
        last_fudge = RECFAST_FUDGE


!       dimensions for integrator
        Ndim = 3

!       write(*,*)'recfast version 1.0'
!       write(*,*)'Using Hummer''s case B recombination rates for H'
!       write(*,*)' with fudge factor = 1.14'
!       write(*,*)'and tabulated HeII singlet recombination rates'
!       write(*,*)

        Tnow=tcmb
!       These are easy to inquire as input, but let's use simple values
        zinitial = 1.d4
        z = zinitial
        zfinal=0._dl
!       will output every 1 in z, but this is easily changed also

  
!       convert the Hubble constant units
        H = H0inp/100._dl
        HO = H*bigH


!       sort out the helium abundance parameters
        mu_H = 1._dl/(1._dl-Yp)         !Mass per H atom
        mu_T = 4._dl/(4._dl-3._dl*Yp)            !Mass per atom
        fHe = Yp/(4._dl*(1._dl-Yp))              !n_He_tot / n_H_tot

        Nnow = 3._dl*HO*HO*OmegaB/(8._dl*Pi*G*mu_H*m_H)
        n = Nnow * (1._dl+z)**3
        !fnu = (21._dl/8._dl)*(4._dl/11._dl)**(4._dl/3._dl)
        !z_eq = 3._dl*(HO*C)**2/(8._dl*Pi*G*a*(1._dl+fnu)*Tnow**4) !This is wrong

      
!       Set up some constants so they don't have to be calculated later
        Lalpha = 1._dl/L_H_alpha
        Lalpha_He = 1._dl/L_He_2p
        DeltaB = h_P*C*(L_H_ion-L_H_alpha)
        CDB = DeltaB/k_B
        DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)   !2s, not 2p
        CDB_He = DeltaB_He/k_B
        CB1 = h_P*C*L_H_ion/k_B
        CB1_He1 = h_P*C*L_He1_ion/k_B   !ionization for HeI
        CB1_He2 = h_P*C*L_He2_ion/k_B   !ionization for HeII
        CR = 2._dl*Pi*(m_e/h_P)*(k_B/h_P)
        CK = Lalpha**3/(8._dl*Pi)
        CK_He = Lalpha_He**3/(8._dl*Pi)
        CL = C*h_P/(k_B*Lalpha)
        CL_He = C*h_P/(k_B/L_He_2s)     !comes from det.bal. of 2s-1s
        CT = (8._dl/3._dl)*(sigma/(m_e*C))*a
        Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

!       Matter departs from radiation when t(Th) > H_frac * t(H)
!       choose some safely small number
        H_frac = 1.D-3

!       Fudge factor to approximate for low z out of equilibrium effect
        fu=RECFAST_fudge

!       Set initial matter temperature
        y(3) = Tnow*(1._dl+z)            !Initial rad. & mat. temperature
        Tmat = y(3)

        call get_init(z,x_H0,x_He0,x0)
    
        y(1) = x_H0
        y(2) = x_He0

!       OK that's the initial conditions, now start writing output file

        w0=1._dl/ sqrt(1._dl + zinitial) !like a conformal time
        w1=1._dl/ sqrt(1._dl + zfinal)
        Lw0 = log(w0)
        Lw1 = log(w1)
        Nz=Nz0
        hW=(Lw1-Lw0)/real(Nz,dl)  !interval in log of conf time

!       Set up work-space stuff for DVERK
        ind  = 1
        nw   = 3
        do i = 1,24
          cw(i) = 0._dl
        end do

        do i = 1,Nz
!       calculate the start and end redshift for the interval at each z
!       or just at each z
          zstart = zinitial  + real(i-1,dl)*(zfinal-zinitial)/real(Nz,dl)
          zend   = zinitial  + real(i,dl)*(zfinal-zinitial)/real(Nz,dl)

! Use Saha to get x_e, using the equation for x_e for ionized helium
! and for neutral helium.
! Everything ionized above z=8000.  First ionization over by z=5000.
! Assume He all singly ionized down to z=3500, then use He Saha until
! He is 99% singly ionized, and *then* switch to joint H/He recombination.

          z = zend
        
          if (zend > 8000._dl) then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = 1._dl+2._dl*fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)

          else if(z > 5000._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He2/(Tnow*(1._dl+z)) ) / Nnow
            rhs = rhs*1._dl            !ratio of g's is 1 for He++ <-> He+
            x0 = 0.5d0 * ( sqrt( (rhs-1._dl-fHe)**2 &
                + 4._dl*(1._dl+2._dl*fHe)*rhs) - (rhs-1._dl-fHe) )
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)

          else if(z > 3500._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = x_H0 + fHe*x_He0
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)

          else if(y(2) > 0.99)then

            x_H0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He1/(Tnow*(1._dl+z)) ) / Nnow
            rhs = rhs*4._dl            !ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5d0 * ( sqrt( (rhs-1._dl)**2 &
                + 4._dl*(1._dl+fHe)*rhs )- (rhs-1._dl))
            x0 = x_He0
            x_He0 = (x0 - 1._dl)/fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)

          else if (y(1) > 0.99d0) then

            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1/(Tnow*(1._dl+z)) ) / Nnow
            x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )

            call DVERK(dum,nw,ION,zstart,y,zend,tol,ind,cw,nw,w)
            y(1) = x_H0
            x0 = y(1) + fHe*y(2)

          else
            
            call DVERK(dum,nw,ION,zstart,y,zend,tol,ind,cw,nw,w)
          
            x0 = y(1) + fHe*y(2)
          
          end if
          

          Trad = Tnow * (1._dl+zend)
          Tmat = y(3)
          x_H = y(1)
          x_He = y(2)
          x = x0

          zrec(i)=zend
          xrec(i)=x
          
        !  write(*,'(4E15.5)') zend, x, y(1), y(2)
     
        end do
      ! stop
        d0hi=1.0d40
        d0lo=1.0d40
        call spline(zrec,xrec,nz,d0lo,d0hi,dxrec)
        
        end subroutine InitRECFAST

!       ===============================================================
        subroutine GET_INIT(z,x_H0,x_He0,x0)

!       Set up the initial conditions so it will work for general,
!       but not pathological choices of zstart
!       Initial ionization fraction using Saha for relevant species
        use RECDATA
        implicit none
  
        
        real(dl) z,x0,rhs,x_H0,x_He0
  

        if(z > 8000._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = 1._dl+2._dl*fHe

        else if(z > 3500._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He2/(Tnow*(1._dl+z)) ) / Nnow
        rhs = rhs*1._dl    !ratio of g's is 1 for He++ <-> He+
        x0 = 0.5d0 * ( sqrt( (rhs-1._dl-fHe)**2 &
                + 4._dl*(1._dl+2._dl*fHe)*rhs) - (rhs-1._dl-fHe) )

        else if(z > 2000._dl)then

        x_H0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He1/(Tnow*(1._dl+z)) ) / Nnow
        rhs = rhs*4._dl    !ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5d0  * ( sqrt( (rhs-1._dl)**2 + 4._dl*(1._dl+fHe)*rhs )- (rhs-1._dl))
            x0 = x_He0
            x_He0 = (x0 - 1._dl)/fHe

        else

            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1/(Tnow*(1._dl+z)) ) / Nnow
            x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )
            x_He0 = 0._dl
            x0 = x_H0

        end if

        
        end subroutine GET_INIT



        subroutine ION(Dum,Ndim,z,Y,f)
        use RECDATA
        implicit none

        integer Ndim

        real(dl) z,x,n,n_He,Trad,Tmat,x_H,x_He, Hz
        real(dl) y(Ndim),f(Ndim)
        real(dl) Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
        real(dl) timeTh,timeH
        real(dl) a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB

        real(dl) dtauda, Dum
        external dtauda
!Stas        
        real(dl) SUMS,SUMP,DL2,DLE,XX,D2S,D2P
        real(dl) C2p1P,C2p3P,C1P3P,cc3P1P,cccP,A2p3P,hck,cc3P,tau,g3P
        integer l,l2,ll
!Stas


!       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen     
        a_PPB = 4.309
        b_PPB =- 0.6166
        c_PPB = 0.6703
        d_PPB = 0.5300
!       the Verner and Ferland type fitting parameters for Helium
        T_0 = 3._dl

   !    a_VF = 10._dl**(-11.7718d0)
   !    b_VF = 1.51930d0
   !    T_1 = 10._dl**(4.50550d0)
   !Bug reported by savita gahlaut
        a_VF = 10.d0**(-16.744d0)
        b_VF = 0.711d0
        T_1 = 10.d0**(5.114d0)
       
        x_H = y(1)
        x_He = y(2)
        x = x_H + fHe * x_He
        Tmat = y(3)

        n = Nnow * (1._dl+z)**3
        n_He = fHe * Nnow * (1._dl+z)**3
        Trad = Tnow * (1._dl+z)

        Hz = 1/dtauda(1/(1._dl+z))*(1._dl+z)**2/1.02928d14
        ! Now use universal background function, fixing factor of Omega_total
        !  Hz = HO * sqrt((1._dl+z)**4/(1+z_eq)*OmegaT + OmegaT*(1._dl+z)**3 &
        !       + OmegaK*(1._dl+z)**2 + OmegaL)
     

!       Get the radiative rates using PPQ fit, identical to Hummer's table
        
        Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
                /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)
        Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)
      
!       calculate He using a fit to a Verner & Ferland type formula
        sq_0 = sqrt(Tmat/T_0)
        sq_1 = sqrt(Tmat/T_1)
!        Rdown_He = a_VF/(sq_0*(1._dl+sq_1)**(1._dl-b_VF))
!        Rdown_He = rdown_He/(1._dl+sq_1)**(1._dl+b_VF)
 !         !Bug reported by savita gahlaut
        Rdown_He = a_VF/(sq_0 * (1.d0+sq_0)**(1.d0-b_VF)* &
                    (1.d0 + sq_1)**(1.d0 + b_VF))


        
        Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*exp(-CDB_He/Tmat)
        Rup_He = 4._dl*Rup_He    !statistical weights factor for HeI   
   
 
        K = CK/Hz              !Peebles coefficient K=lambda_a^3/8piH
        K_He = CK_He/Hz  !Peebles coefficient for Helium

!       Estimates of Thomson scattering time and Hubble time
        timeTh=(1._dl/(CT*Trad**4))*(1._dl+x+fHe)/x       !Thomson time
        timeH=2./(3.*HO*(1._dl+z)**1.5)      !Hubble time

!       calculate the derivatives
!       turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
!       (clunky, but seems to work)
        if (x_H > 0.99) then   !use Saha rate for Hydrogen
                f(1) = 0.
        else if (x_H > 0.98) then
                f(1) = (x*x_H*n*Rdown - Rup*(1.-x_H)*exp(-CL/Tmat)) &
                /(Hz*(1._dl+z))
        else

        if (use_Dubrovich) then
            SUMS=0D0
            SUMP=0D0
            DO L=3,40
            DL2=1D0/(L*L)
            DLE=DEXP(-39450D0*(1D0-4D0*DL2)/(Tnow*(1D0+z)))
             XX=(L-1D0)/(L+1D0)
             XX=XX*XX
             XX=XX**L
             XX=XX*(11D0*L*L-41D0)/L
             SUMS=SUMS+XX*DLE
             SUMP=SUMP+DLE*(1D0-DL2)**3
             ENDDO
             D2S=Lambda*(1D0+54D0*SUMS)
             D2P=K/(1D0+SUMP*64D0/27D0)
         else
           D2S=Lambda
           D2P=K
         end if

        f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat)) &
                *(1.d0 + D2P*D2S*n*(1.d0-x_H))) &
                /(Hz*(1.d0+z)*(1.d0/fu+D2P*D2S*n*(1.d0-x)/fu &
                +D2P*Rup*n*(1.d0-x)))

        end if
   
!       turn off the He once it is small
        if (x_He < 1.e-15) then 
                f(2)=0.
        else

        if (use_Dubrovich) then
            SUMS=0D0
             DO L=6,40
              LL=L+L
              L2=L*L
              DLE=exp(-1.4388*(32033.33D0-109678.77D0/L2)/(Tnow*(1D0+z)))
              XX=(L-1D0)/(L+1D0)
              XX=(XX**LL)*(11D0*L2-41D0)/L
              SUMS=SUMS+XX*DLE
             ENDDO
            D2S=Lambda_He+12D0*1045D0*SUMS
            C2p3P=1.69087D7
            C2p1P=1.71135D7
            C1P3P=C2p1P-C2p3P
            cc3P1P=C2p3P/C2p1P 
            cccP=cc3P1P**3
            A2p3P=233D0
            g3P=3D0
            hck=h_P*C/k_B 
            cc3P=cccP*exp(hck*C1P3P/Trad)
            tau=A2p3P*n_He*(1D0-x_He)*g3P/8D0/PI/(C2p3P**3)/Hz
            D2P=K_He/(1D0+cc3P*(1D0-exp(-tau)))
        else
            D2S=Lambda_He
            D2P=K_He
        end if
!  Stas

        He_Boltz = exp(min(Bfact/Tmat,680._dl))  !Changed by AML Aug/00
           
        f(2) = ((x*x_He*n*Rdown_He &
            - Rup_He*(1-x_He)*exp(-CL_He/Tmat)) &
                *(1 + D2P*D2S*n_He*(1.d0-x_He)*He_Boltz)) &
                /(Hz*(1+z) &
                * (1 + D2P*(D2S+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))

        end if
       
!       follow the matter temperature once it has a chance of diverging
        if (timeTh < H_frac*timeH) then
                f(3)=Tmat/(1._dl+z)      !Tmat follows Trad
        else
                f(3)= CT * (Trad**4) * x / (1._dl+x+fHe) &
                        * (Tmat-Trad) / (Hz*(1._dl+z)) + 2._dl*Tmat/(1._dl+z)
        end if
      

        end subroutine ION


!  ===============================================================
        

        end module RECFAST

