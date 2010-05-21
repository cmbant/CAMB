!     Modules used by cmbmain and other routines.

!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
!     See readme.html for documentation. This version June 2008.
!
!     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
!     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
!     Original CMBFAST copyright and disclaimer:
!
!     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
!     the Massachusetts Institute of Technology.  All rights reserved.
!
!     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
!     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
!     By way of example, but not limitation,
!     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
!     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
!     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
!     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
!
!     portions of this software are based on the COSMICS package of
!     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
!     for restrictions on the modification and distribution of this software.

        module constants
         use precision
        
         real(dl), parameter :: Mpc = 3.085678e22_dl, G=6.6742e-11_dl, kappa=2*fourpi*G, &
                   sigma_thomson = 6.6524616e-29_dl, c = 2.99792458e8_dl, m_p = 1.672623e-27_dl, &
                   sigma_boltz = 5.67051e-8_dl, k_B = 1.380658e-23_dl, m_H = 1.673575e-27_dl

        real(dl), parameter:: barssc0= k_B / m_p / c**2

        end module constants

        module ModelParams
        use precision
        use InitialPower
        use Ranges
        use Reionization
        
        implicit none    
        public

        character(LEN=*), parameter :: version = 'Jun_08'
        
        integer :: FeedbackLevel = 0 !if >0 print out useful information about the model

        logical, parameter :: DebugMsgs=.false. !Set to true to view progress and timing

        logical, parameter :: DebugEvolution = .false. !Set to true to do all the evolution for all k

        integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
         !For CAMBparams%MassiveNuMethod
         !Nu_int: always integrate distribution function
         !Nu_trunc: switch to expansion in velocity once non-relativistic
         !Nu_approx: approximate scheme - good for CMB, but not formally correct and no good for matter power
         !Nu_best: automatically use mixture which is fastest and most accurate

        integer, parameter :: max_Nu = 5 !Maximum number of neutrino species    
        integer, parameter :: max_transfer_redshifts = 128
        integer, parameter :: fileio_unit = 13 !Any number not used elsewhere will do       
        integer, parameter :: outCOBE=0, outNone=1
    
        integer :: max_bessels_l_index  = 1000000
        real(dl) :: max_bessels_etak = 1000000*2


        real(dl), parameter ::  OutputDenominator =twopi
       !When using outNone the output is l(l+1)Cl/OutputDenominator


        Type(Regions) :: TimeSteps


        type TransferParams
            logical     ::  high_precision
            integer     ::  num_redshifts
            real(dl)    ::  kmax         !these are acutally q values, but same as k for CP%flat
            integer     ::  k_per_logint ! ..
            real(dl)    ::  redshifts(max_transfer_redshifts)         
        end type TransferParams

!other variables, options, derived variables, etc.

         integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2

! Main parameters type
        type CAMBparams
    
         logical   :: WantCls, WantTransfer
         logical   :: WantScalars, WantTensors, WantVectors
         logical   :: DoLensing
         integer   :: NonLinear

         integer   :: Max_l, Max_l_tensor
         real(dl)  :: Max_eta_k, Max_eta_k_tensor
          ! _tensor settings only used in initialization, 
          !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

         real(dl)  :: omegab, omegac, omegav, omegan
         !Omega baryon, CDM, Lambda and massive neutrino
         real(dl)  :: H0,TCMB,yhe,Num_Nu_massless,Num_Nu_massive

         logical :: Nu_mass_splittings
         integer   :: Nu_mass_eigenstates  !1 for degenerate masses
         real(dl)  :: Nu_mass_degeneracies(max_nu)
         real(dl)  :: Nu_mass_fractions(max_nu)
             !The ratios of the masses

         integer   :: Scalar_initial_condition 
         !must be one of the initial_xxx values defined in GaugeInterface
         
         integer   :: OutputNormalization  
         !outNone, outCOBE, or C_OutputNormalization=1 if > 1

         logical   :: AccuratePolarization
           !Do you care about the accuracy of the polarization Cls?
  
         logical   :: AccurateBB
           !Do you care about BB accuracy (e.g. in lensing)

!Reionization settings - used if Reion%Reionization=.true.
         logical   :: AccurateReionization
           !Do you care about pecent level accuracy on EE signal from reionization?

         integer   :: MassiveNuMethod
        
         type(InitialPowerParams) :: InitPower  !see power_tilt.f90 - you can change this
         type(ReionizationParams) :: Reion
         type(TransferParams)     :: Transfer 

         real(dl) ::  InitialConditionVector(1:10) !Allow up to 10 for future extensions
          !ignored unless Scalar_initial_condition == initial_vector

         logical OnlyTransfers !Don't use initial power spectrum data, instead get Delta_q_l array
           !If trye, sigma_8 is not calculated either

!Derived parameters, not set initially
         type(ReionizationHistory) :: ReionHist
         logical flat,closed,open
         real(dl) omegak
         real(dl) curv,r, Ksign !CP%r = 1/sqrt(|CP%curv|), CP%Ksign = 1,0 or -1
         real(dl) tau0,chi0 !time today and rofChi(CP%tau0/CP%r) 
    
         
         end type CAMBparams

        type(CAMBparams) CP  !Global collection of parameters


       real(dl) scale !relative to CP%flat. e.g. for scaling lSamp%l sampling.

       logical ::call_again = .false.
          !if being called again with same parameters to get different thing

 
!     grhom =kappa*a^2*rho_m0
!     grhornomass=grhor*number of massless neutrino species
!     taurst,taurend - time at start/end of recombination
!     dtaurec - dtau during recombination
!     adotrad - a(tau) in radiation era

        real(dl) grhom,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
        real(dl) taurst,dtaurec,taurend, tau_maxvis,adotrad

!Neutrinos
        real(dl) grhormass(max_nu)
     
!     nu_masses=m_nu*c**2/(k_B*T_nu0)      
       real(dl) :: nu_masses(max_nu) 
        
       real(dl) akthom !sigma_T * (number density of protons now)

      integer :: ThreadNum = 0 
       !If zero assigned automatically, obviously only used if parallelised
   
!Parameters for checking/changing overall accuracy
!1._dl corresponds to target 1% scalar accuracy for non-extreme models

      real(dl) :: lSampleBoost=1._dl
          !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation
          
      real(dl) :: AccuracyBoost =1._dl  


          !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
          !Can also be used to improve speed significantly if less accuracy is required.              
          !or improving accuracy for extreme models. 
          !Note this does not increase lSamp%l sampling or massive neutrino q-sampling

      real(sp) :: lAccuracyBoost=1. 
          !Boost number of multipoles integrated in Boltzman heirarchy

      integer, parameter :: lmin = 2
          !must be either 1 or 2       

      real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code

      real(dl),parameter :: tol=1.0d-4 !Base tolerance for integrations

!     used as parameter for spline - tells it to use 'natural' end values
      real(dl), parameter :: spl_large=1.e40_dl

      integer, parameter:: l0max=4000

!     lmax is max possible number of l's evaluated
      integer, parameter :: lmax_arr = 100+l0max/7
 
        contains
      

         subroutine CAMBParams_Set(P, error, DoReion)
           use constants
           type(CAMBparams), intent(in) :: P
           real(dl) GetOmegak
           integer, optional :: error !Zero if OK
           logical, optional :: DoReion
           logical WantReion
           integer nu_i
           external GetOmegak
           real(dl), save :: last_tau0
           !Constants in SI units

            if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
              write (*,*) 'Cannot generate tensor C_l and transfer without scalar C_l'
              if (present(error)) then
                error = 1
                return
              else
                stop
              end if
           end if
           
           if (present(DoReion)) then
            WantReion = DoReion
           else
            WantReion = .true.
           end if
        
           CP=P
          
           CP%Max_eta_k = max(CP%Max_eta_k,CP%Max_eta_k_tensor)
           
           if (CP%WantTransfer) then
              CP%WantScalars=.true.
              if (.not. CP%WantCls) then
                 CP%AccuratePolarization = .false.
                 CP%Reion%Reionization = .false.
              end if
           else
              CP%transfer%num_redshifts=0
           end if

           if (CP%Omegan == 0 .and. CP%Num_Nu_Massive /=0) then
              CP%Num_Nu_Massless = CP%Num_Nu_Massless + CP%Num_Nu_Massive
              CP%Num_Nu_Massive  = 0
           end if

           if (CP%Num_nu_massive > 0) then
               if (.not. CP%Nu_mass_splittings) then
                 !Default totally degenerate masses
                 CP%Nu_mass_eigenstates = 1
                 CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive 
                 CP%Nu_mass_fractions(1) = 1
               else
                 if (CP%Nu_mass_degeneracies(1)==0) CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive 
                 if (abs(sum(CP%Nu_mass_fractions(1:CP%Nu_mass_eigenstates))-1) > 1e-4) &
                   stop 'Nu_mass_fractions do not add up to 1'

                 if (abs(sum(CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates))-CP%Num_nu_massive) >1e-4 ) &
                    stop 'nu_mass_eigenstates do not add up to num_nu_massive'
                 if (CP%Nu_mass_eigenstates==0) stop 'Have Num_nu_massive>0 but no nu_mass_eigenstates'

               end if
           else
            CP%Nu_mass_eigenstates = 0
           end if
           
           if ((CP%WantTransfer).and. CP%MassiveNuMethod==Nu_approx) then
              CP%MassiveNuMethod = Nu_trunc
           end if

           CP%omegak = GetOmegak()
          
           CP%flat = (abs(CP%omegak) <= OmegaKFlat)
           CP%closed = CP%omegak < -OmegaKFlat
        
           CP%open = .not.CP%flat.and..not.CP%closed
           if (CP%flat) then
              CP%curv=0
              CP%Ksign=0
              CP%r=1._dl !so we can use tau/CP%r, etc, where CP%r's cancel
           else   
           CP%curv=-CP%omegak/((c/1000)/CP%h0)**2
           CP%Ksign =sign(1._dl,CP%curv)
           CP%r=1._dl/sqrt(abs(CP%curv))
           end if
!  grho gives the contribution to the expansion rate from: (g) photons,
!  (r) one flavor of relativistic neutrino (2 degrees of freedom),
!  (m) nonrelativistic matter (for Omega=1).  grho is actually
!  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
!  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
!  (Used only to set the initial conformal time.)

           !H0 is in km/s/Mpc

           grhom = 3*CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
        
          !grhom=3.3379d-11*h0*h0 
           grhog = kappa/c**2*4*sigma_boltz/c**3*CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
          ! grhog=1.4952d-13*tcmb**4
           grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
          !grhor=3.3957d-14*tcmb**4
           grhornomass=grhor*CP%Num_Nu_massless
           grhormass=0
           do nu_i = 1, CP%Nu_mass_eigenstates
            grhormass(nu_i)=grhor*CP%Nu_mass_degeneracies(nu_i)
           end do
           grhoc=grhom*CP%omegac
           grhob=grhom*CP%omegab
           grhov=grhom*CP%omegav
           grhok=grhom*CP%omegak
!  adotrad gives the relation a(tau) in the radiation era:
           adotrad = sqrt((grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates)))/3)
       
           akthom = sigma_thomson*CP%omegab*(1-CP%yhe)*grhom*c**2/kappa/m_H/Mpc
              !sigma_T * (number density of protons now)
    
           if (CP%omegan==0) then
              CP%Num_nu_massless = CP%Num_nu_massless + CP%Num_nu_massive
              CP%Num_nu_massive = 0
           end if

           if (.not.call_again) then
      
            call init_massive_nu(CP%omegan /=0)
            call init_background
            CP%tau0=TimeOfz(0._dl)
            last_tau0=CP%tau0
            if (WantReion) call Reionization_Init(CP%Reion,CP%ReionHist, CP%YHe, akthom, CP%tau0, FeedbackLevel)
           else
              CP%tau0=last_tau0
           end if         
           
           if ( CP%NonLinear==NonLinear_Lens) then
             CP%Transfer%kmax = max(CP%Transfer%kmax, CP%Max_eta_k/CP%tau0) 
             if (FeedbackLevel > 0 .and. CP%Transfer%kmax== CP%Max_eta_k/CP%tau0) &
                  write (*,*) 'max_eta_k changed to ', CP%Max_eta_k
           end if


           if (CP%closed .and. CP%tau0/CP%r >3.14) then
             if (present(error)) then
              error = 2
              return
             else
              stop 'chi >= pi in closed model not supported'
             end if
           end if
    
           if (present(error)) then
              error = 0
           else if (FeedbackLevel > 0 .and. .not. call_again) then
              write(*,'("Om_b h^2             = ",f9.6)') CP%omegab*(CP%H0/100)**2
              write(*,'("Om_c h^2             = ",f9.6)') CP%omegac*(CP%H0/100)**2
              write(*,'("Om_nu h^2            = ",f9.6)') CP%omegan*(CP%H0/100)**2
              write(*,'("Om_Lambda            = ",f9.6)') CP%omegav
              write(*,'("Om_K                 = ",f9.6)') CP%omegak
              write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-CP%omegak-CP%omegav
              write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta()
              if (CP%Num_Nu_Massive > 0) then
                do nu_i=1, CP%Nu_mass_eigenstates 
                 write(*,'(f5.2, " nu, m_nu*c^2/k_B/T_nu0   = ",f8.2," (m_nu = ",f6.3," eV)")') &
                     CP%nu_mass_degeneracies(nu_i), nu_masses(nu_i),1.68e-4*nu_masses(nu_i)
                end do
              end if
           end if
           CP%chi0=rofChi(CP%tau0/CP%r)
           scale= CP%chi0*CP%r/CP%tau0  !e.g. changel sampling depending on approx peak spacing      
           
         end subroutine CAMBParams_Set

    
         function GetTestTime()
           real(sp) GetTestTime
           real(sp) atime

!           GetTestTime = etime(tarray)
         !Can replace this if etime gives problems
         !Or just comment out - only used if DebugMsgs = .true.
           call cpu_time(atime)
           GetTestTime = atime
            
         end function GetTestTime

        
         function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
         real(dl) Chi,rofChi

         if (CP%closed) then
            rofChi=sin(chi)
         else if (CP%open) then
            rofChi=sinh(chi)
         else
            rofChi=chi
         endif
         end function rofChi  
         

         function cosfunc (Chi)
         real(dl) Chi,cosfunc

         if (CP%closed) then
            cosfunc= cos(chi)
         else if (CP%open) then
            cosfunc=cosh(chi)
         else
            cosfunc = 1._dl
         endif
         end function cosfunc  

         function tanfunc(Chi)
         real(dl) Chi,tanfunc
         if (CP%closed) then
            tanfunc=tan(Chi)
         else if (CP%open) then
            tanfunc=tanh(Chi)
         else
            tanfunc=Chi
         end if

         end  function tanfunc

         function invsinfunc(x)
         real(dl) invsinfunc,x

         if (CP%closed) then
          invsinfunc=asin(x)
          else if (CP%open) then
          invsinfunc=log((x+sqrt(1._dl+x**2)))  
          else
          invsinfunc = x
         endif
         end function invsinfunc    

        function f_K(x)
         real(dl) :: f_K
         real(dl), intent(in) :: x
         f_K = CP%r*rofChi(x/CP%r)
          
        end function f_K


        function DeltaTime(a1,a2)
        implicit none
        real(dl) DeltaTime, atol
        real(dl), intent(IN) :: a1,a2
        real(dl) dtauda, rombint !diff of tau w.CP%r.t a and integration
        external dtauda, rombint

        atol = tol/1000/exp(AccuracyBoost-1)
        DeltaTime=rombint(dtauda,a1,a2,atol)
      
        end function DeltaTime

        function TimeOfz(z)
        implicit none
        real(dl) TimeOfz
        real(dl), intent(IN) :: z
        
        TimeOfz=DeltaTime(0._dl,1._dl/(z+1._dl))
        end function TimeOfz

        function AngularDiameterDistance(z)
          real(dl) AngularDiameterDistance
          real(dl), intent(in) :: z

          AngularDiameterDistance = CP%r/(1+z)*rofchi(DeltaTime(1/(1+z),1._dl)/CP%r)

        end function AngularDiameterDistance

       function dsound_da(a)
          implicit none
          real(dl) dsound_da,dtauda,a,R,cs
          external dtauda

           R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
           cs=1.0d0/sqrt(3*(1+R))
           dsound_da=dtauda(a)*cs
        
       end function dsound_da


       function CosmomcTheta()
         real(dl) zstar, astar, atol, rs, DA
         real(dl) CosmomcTheta
         real(dl) ombh2, omdmh2
         real(dl) rombint
         external rombint

         ombh2 = CP%omegab*(CP%h0/100.0d0)**2
         omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2


    !!From Hu & Sugiyama
           zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
            (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
               (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))
     
           astar = 1/(1+zstar)
           atol = 1e-6
           rs = rombint(dsound_da,1d-8,astar,atol)
           DA = AngularDiameterDistance(zstar)/astar
           CosmomcTheta = rs/DA
    !       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

      end function CosmomcTheta

   end module ModelParams



!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module lvalues
        use precision
        use ModelParams
        implicit none
        public

        Type lSamples
            integer l0
           ! integer, dimension(:), pointer :: lSamp%l 
            integer l(lmax_arr)
        end Type lSamples

        Type(lSamples) :: lSamp

       contains


        subroutine initlval(lSet,max_l)

! This subroutines initializes lSet%l arrays. Other values will be interpolated.
  
        implicit none
        type(lSamples) :: lSet
         
        integer, intent(IN) :: max_l
        integer lind, lvar, step,top,bot,ls(lmax_arr)
        real(dl) AScale
    
        Ascale=scale/lSampleBoost       
      
        lind=0
        do lvar=lmin, 10
           lind=lind+1
           ls(lind)=lvar 
        end do

        if (CP%AccurateReionization) then
             do lvar=11, 37,2
               lind=lind+1
               ls(lind)=lvar 
             end do       

            step = max(nint(5*Ascale),2)           
            bot=40
            top=bot + step*10
        else

            if (lSampleBoost >1) then
             do lvar=11, 15
               lind=lind+1
               ls(lind)=lvar 
             end do           
            else
             lind=lind+1
             ls(lind)=12
             lind=lind+1
             ls(lind)=15
            end if
            step = max(nint(10*Ascale),3)           
            bot=15+max(step/2,2)
            top=bot + step*7
        end if

        do lvar=bot, top, step
           lind=lind+1
           ls(lind)=lvar          
        end do

        step=max(nint(20*Ascale),4)
        bot=ls(lind)+step
        top=bot+step*2

        do lvar = bot,top,step 
          lind=lind+1
          ls(lind)=lvar
        end do

        if (ls(lind)>=max_l) then
           do lvar=lind,1,-1
            if (ls(lvar)<=max_l) exit  
           end do
           lind=lvar
           if (ls(lind)<max_l) then
              lind=lind+1
              ls(lind)=max_l
           end if
        else

        step=max(nint(25*Ascale),4)
!Get EE right around l=200 by putting extra point at 175
        bot=ls(lind)+step
        top=bot+step

        do lvar = bot,top,step 
          lind=lind+1
          ls(lind)=lvar
        end do


        if (ls(lind)>=max_l) then
           do lvar=lind,1,-1
            if (ls(lvar)<=max_l) exit  
           end do
           lind=lvar
           if (ls(lind)<max_l) then
              lind=lind+1
              ls(lind)=max_l
           end if
        else

        step=max(nint(50*Ascale),7)
        bot=ls(lind)+step
        top=min(5000,max_l)

         do lvar = bot,top,step
          lind=lind+1
          ls(lind)=lvar
         end do

         if (max_l > 5000) then
             !Should be pretty smooth or tiny out here   
             step=max(nint(400*Ascale),50)
             lvar = ls(lind)
            
             do
              lvar = lvar + step
              if (lvar > max_l) exit
              lind=lind+1
              ls(lind)=lvar
              step = nint(step*1.5) !log spacing
             end do

         end if

         if (ls(lind) /=max_l) then          
           lind=lind+1
           ls(lind)=max_l
         end if
        if (.not. CP%flat) ls(lind-1)=int(max_l+ls(lind-2))/2
        !Not in CP%flat case so interpolation table is the same when using lower l_max
        end if
        end if
        lSet%l0=lind
        lSet%l(1:lind) = ls(1:lind)
        
      end subroutine initlval

      subroutine InterpolateClArr(lSet,iCl, all_Cl, max_ind)
      type (lSamples), intent(in) :: lSet        
      real(dl), intent(in) :: iCl(*)
      real(dl), intent(out):: all_Cl(lmin:*)
      integer, intent(in) :: max_ind
      integer il,llo,lhi, xi
      real(dl) ddCl(lSet%l0)
      real(dl) xl(lSet%l0)

      real(dl) a0,b0,ho
      real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

      if (max_ind > lSet%l0) stop 'Wrong max_ind in InterpolateClArr'

      xl = real(lSet%l(1:lSet%l0),dl)
      call spline(xl,iCl(1),max_ind,cllo,clhi,ddCl(1))
     
            llo=1
            do il=lmin,lSet%l(max_ind)
               xi=il
               if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
                  llo=llo+1
               end if
               lhi=llo+1
               ho=lSet%l(lhi)-lSet%l(llo)
               a0=(lSet%l(lhi)-xi)/ho
               b0=(xi-lSet%l(llo))/ho
      
               all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
                       +(b0**3-b0)*ddCl(lhi))*ho**2/6
              
            end do

      end subroutine InterpolateClArr

 
    

!ccccccccccccccccccccccccccc

        end module lvalues
        


!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module ModelData
        use precision
        use ModelParams
        use InitialPower
        use lValues
        use Ranges
        use AMlUtils
        implicit none
        public

         Type ClTransferData
      !Cl transfer function variables
       !values of q for integration over q to get C_ls
          Type (lSamples) :: ls ! scalar and tensor l that are computed
          integer :: NumSources 
          !Changes -scalars:  2 for just CMB, 3 for lensing
          !- tensors: T and E and phi (for lensing), and T, E, B respectively
        
          Type (Regions) :: q

          real(dl), dimension(:,:,:), pointer :: Delta_p_l_k
      
         end Type ClTransferData

         Type(ClTransferData), target :: CTransScal, CTransTens, CTransVec

        !Computed output power spectra data
                         
        integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4, C_PhiTemp = 5
        integer :: C_last = C_PhiTemp
        integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4
  

        real(dl), dimension (:,:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , intial_power_index, Cl_type)
        !where Cl_type is one of the above constants

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this 
        real(dl) , dimension (:,:,:), allocatable :: Cl_lensed
          !Cl_lensed(l, power_index, Cl_type) are the interpolated Cls
    
        real(dl), dimension (:), allocatable ::  COBElikelihoods,COBE_scales
        !Set by COBEnormalize if using outCOBE
        contains


        subroutine Init_ClTransfer(CTrans)
        !Need to set the Ranges array q before calling this
          Type(ClTransferData) :: CTrans
          integer st

          deallocate(CTrans%Delta_p_l_k, STAT = st)
          call Ranges_getArray(CTrans%q, .true.)

          allocate(CTrans%Delta_p_l_k(CTrans%NumSources,min(max_bessels_l_index,CTrans%ls%l0), CTrans%q%npoints))
          CTrans%Delta_p_l_k = 0

  
         end subroutine Init_ClTransfer


        subroutine Free_ClTransfer(CTrans)
          Type(ClTransferData) :: CTrans
          integer st

           deallocate(CTrans%Delta_p_l_k, STAT = st)
           nullify(CTrans%Delta_p_l_k)
           call Ranges_Free(CTrans%q)

        end subroutine Free_ClTransfer



        subroutine Init_Cls
      
        if (CP%WantScalars) then
         if (allocated(Cl_scalar)) deallocate(Cl_scalar)
         allocate(Cl_scalar(lmin:CP%Max_l, CP%InitPower%nn, C_Temp:C_last))
         Cl_scalar = 0
        end if

        if (CP%WantVectors) then
         if (allocated(Cl_vector)) deallocate(Cl_vector)
         allocate(Cl_vector(lmin:CP%Max_l, CP%InitPower%nn, CT_Temp:CT_Cross))
         Cl_vector = 0
        end if


        if (CP%WantTensors) then
          if (allocated(Cl_tensor)) deallocate(Cl_tensor)
          allocate(Cl_tensor(lmin:CP%Max_l_tensor, CP%InitPower%nn, CT_Temp:CT_Cross))
          Cl_tensor = 0
        end if

        end subroutine Init_Cls
       
        subroutine output_cl_files(ScalFile,TensFile, TotFile, LensFile, LensTotFile, factor)
        implicit none
        integer in,il
        character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile
        real(dl), intent(in), optional :: factor
        real(dl) fact


        if (present(factor)) then
          fact = factor
        else
          fact =1
        end if

         if (CP%WantScalars .and. ScalFile /= '') then
 
           open(unit=fileio_unit,file=ScalFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,min(10000,CP%Max_l)
               write(fileio_unit,trim(numcat('(1I6,',C_last))//'E15.5)')il ,fact*Cl_scalar(il,in,C_Temp:C_last)
             end do
             do il=10100,CP%Max_l, 100
               write(fileio_unit,trim(numcat('(1E15.5,',C_last))//'E15.5)') real(il) ,fact*Cl_scalar(il,in,C_Temp:C_last)
             end do
            end do
            close(fileio_unit)
         end if
  
       if (CP%WantTensors .and. TensFile /= '') then
           open(unit=fileio_unit,file=TensFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l_tensor
               write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_tensor(il, in, CT_Temp:CT_Cross)
             end do
            end do
           close(fileio_unit)
        end if
 
        if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
           open(unit=fileio_unit,file=TotFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l_tensor

                write(fileio_unit,'(1I6,4E15.5)')il, fact*(Cl_scalar(il, in, C_Temp:C_E)+ Cl_tensor(il,in, C_Temp:C_E)), &
                   fact*Cl_tensor(il,in, CT_B), fact*(Cl_scalar(il, in, C_Cross) + Cl_tensor(il, in, CT_Cross))
             end do
             do il=CP%Max_l_tensor+1,CP%Max_l
                  write(fileio_unit,'(1I6,4E15.5)')il ,fact*Cl_scalar(il,in,C_Temp:C_E), 0._dl, fact*Cl_scalar(il,in,C_Cross)
             end do
           end do
           close(fileio_unit)
        end if
 
        if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
           open(unit=fileio_unit,file=LensFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin, lmax_lensed
               write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
             end do
            end do
           close(fileio_unit)      
        end if


       if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
           open(unit=fileio_unit,file=LensTotFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,min(CP%Max_l_tensor,lmax_lensed)
                write(fileio_unit,'(1I6,4E15.5)')il, fact*(Cl_lensed(il, in, CT_Temp:CT_Cross)+ Cl_tensor(il,in, CT_Temp:CT_Cross))
             end do
             do il=min(CP%Max_l_tensor,lmax_lensed)+1,lmax_lensed
                write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
             end do
           end do
     
        end if



        end subroutine output_cl_files

        subroutine output_veccl_files(VecFile, factor)
        implicit none
        integer in,il
        character(LEN=*) VecFile
        real(dl), intent(in), optional :: factor
        real(dl) fact


        if (present(factor)) then
          fact = factor
        else
          fact =1
        end if

  
       if (CP%WantVectors .and. VecFile /= '') then
           open(unit=fileio_unit,file=VecFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l
               write(fileio_unit,'(1I5,4E15.5)')il, fact*Cl_vector(il, in, CT_Temp:CT_Cross)
             end do
            end do

           close(fileio_unit)
        end if
 
        end subroutine output_veccl_files


        subroutine output_COBElikelihood
          integer in
          do in=1, CP%InitPower%nn
             write(*,*)'COBE Likelihood relative to CP%flat=',COBElikelihoods(in)
          end do
        end  subroutine output_COBElikelihood

        
      subroutine NormalizeClsAtL(lnorm)
        implicit none
        integer, intent(IN) :: lnorm
        integer in
        real(dl) Norm

         do in=1,CP%InitPower%nn
             
             if (CP%WantScalars) then
                Norm=1/Cl_scalar(lnorm,in, C_Temp)
                Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) = Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) * Norm
             end if

             if (CP%WantTensors) then
                  if (.not.CP%WantScalars) Norm = 1/Cl_tensor(lnorm,in, C_Temp)
                  !Otherwise Norm already set correctly
                  Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) =  &
                    Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) * Norm
             end if
       end do

      end  subroutine NormalizeClsAtL

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine COBEnormalize
        use precision
        use ModelParams
        

        integer in
        real(dl) xlog10
        real(dl) c10, d1,d2,d3,d4,d5,d6,d7, xlogl, COBE_scale
        real(dl) x1, x2,x3,x4,x5,x6,x7,sy,s,sx,sxy,sxx,delt,d1pr,d1ppr
        real(dl) Ctot(lmin:20)

    
           if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
           if (allocated(COBE_scales)) deallocate(COBE_scales)
           allocate(COBElikelihoods(CP%InitPower%nn))
           allocate(COBE_scales(CP%InitPower%nn))
    

          
        xlog10=log(10._dl)
  

! COBE normalization
! fit the spectrum to a quadratic around C_10 with equal weights in logl

        do in=1,CP%InitPower%nn

           if (CP%WantTensors) then
              Ctot =  Cl_tensor(lmin:20, in, C_Temp)
           else
              Ctot = 0
           end if
           if (CP%WantScalars) then
              Ctot=Ctot + Cl_scalar(lmin:20, in, C_Temp)
     
           end if
           c10=Ctot(10)
      
           d1=(Ctot(3))/c10-1._dl
           d2=(Ctot(4))/c10-1._dl
           d3=(Ctot(6))/c10-1._dl
           d4=(Ctot(8))/c10-1._dl
           d5=(Ctot(12))/c10-1._dl
           d6=(Ctot(15))/c10-1._dl
           d7=(Ctot(20))/c10-1._dl

     
           x1=log(3._dl)/xlog10-1._dl
           x2=log(4._dl)/xlog10-1._dl
           x3=log(6._dl)/xlog10-1._dl
           x4=log(8._dl)/xlog10-1._dl
           x5=log(12._dl)/xlog10-1._dl
           x6=log(15._dl)/xlog10-1._dl
           x7=log(20._dl)/xlog10-1._dl
           sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7
           s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7
           sx=x1**3+x2**3+x3**3+x4**3+x5**3+x6**3+x7**3
           sxy=x1**2*d1+x2**2*d2+x3**2*d3+x4**2*d4+ &
              x5**2*d5+x6**2*d6+x7**2*d7
           sxx=x1**4+x2**4+x3**4+x4**4+x5**4+x6**4+x7**4
           delt=s*sxx-sx*sx
           d1pr=(sxx*sy-sx*sxy)/delt
           d1ppr=2._dl*(s*sxy-sx*sy)/delt

! Bunn and White fitting formula
           c10=(0.64575d0+0.02282d0*d1pr+0.01391d0*d1pr*d1pr &
           -0.01819d0*d1ppr-0.00646d0*d1pr*d1ppr &
           +0.00103d0*d1ppr*d1ppr)/c10
! logl
           xlogl=-0.01669d0+1.19895d0*d1pr-0.83527d0*d1pr*d1pr &
                 -0.43541d0*d1ppr-0.03421d0*d1pr*d1ppr &
                 +0.01049d0*d1ppr*d1ppr
          ! write(*,*)'COBE Likelihood relative to CP%flat=',exp(xlogl)
           COBElikelihoods(in) = exp(xlogl)

! density power spectrum normalization;

           COBE_scale=c10/OutputDenominator*1.1d-9
           COBE_scales(in)=COBE_scale

!!$!delta^2 = k^4*(tf)^2*ScalarPower(k,in)*COBE_scale where (tf) is output in the transfer function file
!!$!delta^2 = 4*pi*k^3 P(k)


! C_l normalization; output l(l+1)C_l/twopi
           c10=c10*2.2d-9/fourpi

           if (CP%WantScalars) Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last) = &
                        Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last)*c10
           if (CP%WantTensors) Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) = &
                                    Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross)*c10
           
          end do !in
         end subroutine COBEnormalize
        
         subroutine ModelData_Free

             call Free_ClTransfer(CTransScal)
             call Free_ClTransfer(CTransVec)
             call Free_ClTransfer(CTransTens)
             if (allocated(Cl_vector)) deallocate(Cl_vector)
             if (allocated(Cl_tensor)) deallocate(Cl_tensor)
             if (allocated(Cl_scalar)) deallocate(Cl_scalar)
             if (allocated(Cl_lensed)) deallocate(Cl_lensed)
             if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
             if (allocated(COBE_scales)) deallocate(COBE_scales) 
 
         end subroutine ModelData_Free

        end module ModelData


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    module MassiveNu
      use precision
      use ModelParams
      implicit none
        private 
        
          real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
             !const = int q^3 F(q) dq = 7/120*pi^4
          real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
          real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
          real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
          real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

          integer, parameter  :: nrhopn=2000  
          real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
            !smallest a*m_nu to integrate distribution function rather than using series
          real(dl), parameter :: am_max = 600._dl 
            !max a*m_nu to integrate
          
          real(dl),parameter  :: am_minp=am_min*1.1
          real(dl), parameter :: am_maxp=am_max*0.9
   
          real(dl) dlnam

          real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1,qdn
          
          real(dl), parameter :: dq=1._dl  
          !Sample for massive neutrino momentum; increase nqmax0 appropriately
          !These settings appear to be OK for P_k accuate at 1e-3 level
          integer, parameter :: nqmax0=15 !number of q to sample for each l

          real(dl) dlfdlq(nqmax0) !pre-computed array
          integer nqmax
 
       public Nu_Init,Nu_background,Nu_Integrate,Nu_Integrate01,Nu_Intvsq, &
              Nu_Shear,Nu_derivs, Nu_rho, nqmax, dlfdlq, dq, nqmax0
       contains
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine Nu_init
      
!  Initialize interpolation tables for massive neutrinos.
!  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
 
         integer i
         real(dl) q, am, rhonu,pnu
         real(dl) spline_data(nrhopn)

     
!  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
!  Get number density n of neutrinos from
!  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
!  then m = Omega_nu/N_nu rho_crit /n
!  Error due to velocity < 1e-5
        
        do i=1, CP%Nu_mass_eigenstates 
         nu_masses(i)=const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
               /CP%Nu_mass_degeneracies(i)
        end do

        if (allocated(r1)) return
         

        allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn),qdn(nqmax0))

         do i=1,nqmax0
            q=(i-0.5d0)*dq
            dlfdlq(i)=-q/(1._dl+exp(-q))
            qdn(i)=dq*q**3/(exp(q)+1._dl)
         end do

        dlnam=-(log(am_min/am_max))/(nrhopn-1)
 

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
        !$OMP & PRIVATE(am, rhonu,pnu) 
        do i=1,nrhopn
          am=am_min*exp((i-1)*dlnam)
          call nuRhoPres(am,rhonu,pnu)
          r1(i)=log(rhonu)
          p1(i)=log(pnu)
        end do
        !$OMP END PARALLEL DO


        call splini(spline_data,nrhopn)
        call splder(r1,dr1,nrhopn,spline_data)
        call splder(p1,dp1,nrhopn,spline_data)
        call splder(dr1,ddr1,nrhopn,spline_data)       

       
        end subroutine Nu_init


      

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nuRhoPres(am,rhonu,pnu)
!  Compute the density and pressure of one eigenstate of massive neutrinos,
!  in units of the mean density of one flavor of massless neutrinos.

        real(dl),  parameter :: qmax=30._dl
        integer, parameter :: nq=1000
        real(dl) dum1(nq+1),dum2(nq+1)
        real(dl), intent(in) :: am
        real(dl), intent(out) ::  rhonu,pnu
        integer i
        real(dl) q,aq,v,aqdn,adq
       

!  q is the comoving momentum in units of k_B*T_nu0/c.
!  Integrate up to qmax and then use asymptotic expansion for remainder.
        adq=qmax/nq
        dum1(1)=0._dl
        dum2(1)=0._dl
        do  i=1,nq
          q=i*adq
          aq=am/q
          v=1._dl/sqrt(1._dl+aq*aq)
          aqdn=adq*q*q*q/(exp(q)+1._dl)
          dum1(i+1)=aqdn/v
          dum2(i+1)=aqdn*v
        end do
        call splint(dum1,rhonu,nq+1)
        call splint(dum2,pnu,nq+1)
!  Apply asymptotic corrrection for q>qmax and normalize by relativistic
!  energy density.
        rhonu=(rhonu+dum1(nq+1)/adq)/const
        pnu=(pnu+dum2(nq+1)/adq)/const/3._dl
       
        end subroutine nuRhoPres

!cccccccccccccccccccccccccccccccccccccccccc
       subroutine Nu_background(am,rhonu,pnu)
        use precision
        use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu, pnu

!  Compute massive neutrino density and pressure in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i
      
        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2  
          pnu=(2-rhonu)/3._dl
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          pnu = 900._dl/120._dl/const*(zeta5-63._dl/4*Zeta7/am**2)/am
          return
        end if

        
        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i
       
!  Cubic spline interpolation.
          rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
          pnu=p1(i)+d*(dp1(i)+d*(3._dl*(p1(i+1)-p1(i))-2._dl*dp1(i) &
               -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2._dl*(p1(i)-p1(i+1)))))
          rhonu=exp(rhonu)
          pnu=exp(pnu)

        end subroutine Nu_background

!cccccccccccccccccccccccccccccccccccccccccc
       subroutine Nu_rho(am,rhonu)
        use precision
        use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu

!  Compute massive neutrino density in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i
      
        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2  
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          return
        end if
        
        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i
       
!  Cubic spline interpolation.
        rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
        rhonu=exp(rhonu)
       end subroutine Nu_rho


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Nu_Integrate01(am,drhonu,fnu,psi0,psi1)
        use precision
        use ModelParams
    
!  Compute the perturbations of density and energy flux
!  of one eigenstate of massive neutrinos, in units of the mean
!  density of one eigenstate of massless neutrinos, by integrating over
!  momentum.
        real(dl), intent(IN)  :: am,psi0(nqmax0),psi1(nqmax0) 
        real(dl), intent(OUT) ::  drhonu,fnu
        real(dl), parameter   :: qmax=nqmax0-0.5d0
    

        real(dl) g0(4),g1(nqmax0+1),g3(nqmax0+1)
        real(dl) aq,v,q
        integer iq
   
!  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0
        g3(1)=0
        do iq=2,(nqmax0+1)
            q=(iq-1.5d0)*dq
            aq=am/q
            v=1._dl/sqrt(1._dl+aq*aq)          
            g1(iq)=qdn(iq-1)*psi0(iq-1)/v           
            g3(iq)=qdn(iq-1)*psi1(iq-1)           
        end do
        call splint(g1,g0(1),nqmax0+1)
        call splint(g3,g0(3),nqmax0+1)
  
        drhonu=(g0(1)+g1(nqmax0+1)*2._dl/qmax)/const
        fnu=(g0(3)+g3(nqmax0+1)*2._dl/qmax)/const
     
        end subroutine Nu_Integrate01


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Nu_Integrate(am,drhonu,fnu,dpnu,shearnu,psi0,psi1,psi2)
        use precision
        use ModelParams
    
!  Compute the perturbations of density, energy flux, pressure, and
!  shear stress of one eigenstate of massive neutrinos, in units of the mean
!  density of one eigenstate of massless neutrinos, by integrating over
!  momentum.
        real(dl), intent(IN)  :: am,psi0(nqmax0),psi1(nqmax0),psi2(nqmax0) 
        real(dl), intent(OUT) ::  drhonu,fnu,dpnu,shearnu
        real(dl), parameter   :: qmax=nqmax0-0.5d0
    

        real(dl) g0(4),g1(nqmax0+1),g2(nqmax0+1)
        real(dl) g3(nqmax0+1),g4(nqmax0+1)
        real(dl) aq,v,q
        integer iq
   

!  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0._dl
        g2(1)=0._dl
        g3(1)=0._dl
        g4(1)=0._dl
        do iq=2,(nqmax0+1)
            q=(iq-1.5d0)*dq
            aq=am/q
            v=1._dl/sqrt(1._dl+aq*aq)          
            g1(iq)=qdn(iq-1)*psi0(iq-1)/v
            g2(iq)=qdn(iq-1)*psi0(iq-1)*v
            g3(iq)=qdn(iq-1)*psi1(iq-1)
            g4(iq)=qdn(iq-1)*psi2(iq-1)*v         
        end do
        call splint(g1,g0(1),nqmax0+1)
        call splint(g2,g0(2),nqmax0+1)
        call splint(g3,g0(3),nqmax0+1)
        call splint(g4,g0(4),nqmax0+1)
  
        drhonu=(g0(1)+g1(nqmax0+1)*2._dl/qmax)/const
        fnu=(g0(3)+g3(nqmax0+1)*2._dl/qmax)/const
        dpnu=(g0(2)+g2(nqmax0+1)*2._dl/qmax)/const/3._dl
        shearnu=(g0(4)+g4(nqmax0+1)*2._dl/qmax)/const*2._dl/3._dl

        end subroutine Nu_Integrate

!cccccccccccccccccccccccccccccccccccccccccccccc
   subroutine Nu_Intvsq(am,G11,G30,psi1,psi3)
        use precision
        use ModelParams

!  Compute the third order variables (in velocity dispersion) 
!by integrating over momentum.

        real(dl), intent(IN) :: am
        real(dl), parameter :: qmax=nqmax0-0.5d0
        real(dl)  psi1(nqmax0),psi3(nqmax0)

        real(dl) g0(4),g1(nqmax0+1),g2(nqmax0+1)
            
        real(dl) G11,G30
        real(dl) aq,q,v
        integer iq

!  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0._dl
        g2(1)=0._dl
   
        do iq=2,(nqmax0+1)
            q=(iq-1.5d0)*dq
            aq=am/q
            v=1._dl/sqrt(1._dl+aq*aq)          
            g1(iq)=qdn(iq-1)*psi1(iq-1)*v**2
            g2(iq)=qdn(iq-1)*psi3(iq-1)*v**2
         
        end do
        call splint(g1,g0(1),nqmax0+1)
        call splint(g2,g0(2),nqmax0+1)
             
        G11=(g0(1)+g1(nqmax0+1)*2._dl/qmax)/const
        G30=(g0(2)+g2(nqmax0+1)*2._dl/qmax)/const
       
        end subroutine Nu_Intvsq

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Nu_Shear(am,shearnu,psi2)
        use precision
        use ModelParams

!  Compute the perturbations of 
!  shear stress of one eigenstate of massive neutrinos, in units of the mean
!  density of one eigenstate of massless neutrinos, by integrating over
!  momentum.
        real(dl), intent(IN) :: am 
     
        real(dl), parameter :: qmax=nqmax0-0.5d0
        real(dl) psi2(nqmax0)

        real(dl) g0(4),g4(nqmax0+1)
        real(dl) shearnu,q,aq,v
        integer iq

        if (nqmax==0) then
          shearnu=0._dl
          return
        end if
!
!  q is the comoving momentum in units of k_B*T_nu0/c.
        g4(1)=0._dl      
        do  iq=2,(nqmax0+1)
            q=(iq-1.5d0)*dq
            aq=am/q
            v=1._dl/sqrt(1._dl+aq*aq)                     
            g4(iq)=qdn(iq-1)*psi2(iq-1)*v         
        end do
        call splint(g4,g0(4),nqmax0+1)       
        shearnu=(g0(4)+g4(nqmax0+1)*2._dl/qmax)/const*2._dl/3._dl

        end subroutine Nu_Shear
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine Nu_derivs(am,adotoa,rhonu,rhonudot,shearnudot,psi2,psi2dot)
        use precision
        use ModelParams

!  Compute the time derivative of the mean density in massive neutrinos
!  and the shear perturbation.
!
      
        real(dl), parameter :: qmax=nqmax0-0.5d0
        real(dl) psi2(nqmax0),psi2dot(nqmax0)

        real(dl) g1(nqmax0+1)
        real(dl) adotoa,rhonu,rhonudot,shearnudot
        real(dl) aq,q,v,d,aqdot,vdot,g0
        real(dl), intent(IN) :: am
        integer iq,i
     

!  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0._dl
        do iq=2,(nqmax0+1)
            q=(iq-1.5d0)*dq
            aq=am/q
            aqdot=aq*adotoa
            v=1._dl/sqrt(1._dl+aq*aq)
            vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
            g1(iq)=qdn(iq-1)*(psi2dot(iq-1)*v+psi2(iq-1)*vdot)
        end do
        call splint(g1,g0,nqmax0+1)
     
        shearnudot=(g0+g1(nqmax0+1)*2._dl/qmax)/const*2._dl/3._dl

        if (am< am_minp) then

           rhonudot = 2*const2*am**2*adotoa

        else if (am>am_maxp) then

           rhonudot = 3/(2*const)*(zeta3*am - (15*zeta5)/2/am)*adotoa

        else
           
           d=log(am/am_min)/dlnam+1._dl
           i=int(d)
           d=d-i
           !  Cubic spline interpolation for rhonudot.
           rhonudot=dr1(i)+d*(ddr1(i)+d*(3._dl*(dr1(i+1)-dr1(i)) &
                -2._dl*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1) &
                +2._dl*(dr1(i)-dr1(i+1)))))
     
           rhonudot=rhonu*adotoa*rhonudot/dlnam
        end if
      
        end subroutine Nu_derivs

      end module MassiveNu

! wrapper function to avoid cirular module references
      subroutine init_massive_nu(has_massive_nu)
        use MassiveNu
        use ModelParams
        implicit none
        logical, intent(IN) :: has_massive_nu

        if (has_massive_nu) then
             nqmax=nqmax0 
             call Nu_Init  
        else
             nu_masses = 0
             nqmax= 0
        end if
      end subroutine init_massive_nu


!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module Transfer
        use ModelData
        implicit none
        public
        integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
                              Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
                              Transfer_tot=7
       
        integer, parameter :: Transfer_max = Transfer_tot

        logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
         !set to false to output calculated values for later interpolation

        integer :: transfer_power_var = Transfer_tot 
         !What to use to calulcate the output matter power spectrum and sigma_8
         !Transfer_tot uses total matter perturbation

        Type MatterTransferData
         !Computed data
         integer   ::  num_q_trans   !    number of steps in k for transfer calculation
         real(dl), dimension (:), pointer :: q_trans
         real(dl), dimension (:,:), pointer ::  sigma_8
         real, dimension(:,:,:), pointer :: TransferData
         !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
        end Type MatterTransferData

        Type MatterPowerData
         !everything is a function of k/h
          integer   ::  num_k, num_z          
          real(dl), dimension(:), pointer :: log_kh, redshifts
          !matpower is log(P_k)
          real(dl), dimension(:,:), pointer :: matpower, ddmat
          !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
          !function of k and redshift NonLinearScaling(k_index,z_index)         
          real(dl), dimension(:,:), pointer :: nonlin_ratio
        end Type MatterPowerData

        Type (MatterTransferData) :: MT        

      contains

        subroutine Transfer_GetMatterPowerData(MTrans, PK_data, in, itf_only)
         !Does *NOT* include non-linear corrections
          !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PK_data
          integer, intent(in) :: in
          integer, intent(in), optional :: itf_only
          real(dl) h, kh, k
          integer ik
          integer nz,itf, itf_start, itf_end
          
          if (present(itf_only)) then
              itf_start=itf_only
              itf_end = itf_only
              nz = 1
          else
              itf_start=1
              nz= size(MTrans%TransferData,3)
              itf_end = nz
          end if
          PK_data%num_k = MTrans%num_q_trans
          PK_Data%num_z = nz

          allocate(PK_data%matpower(PK_data%num_k,nz))
          allocate(PK_data%ddmat(PK_data%num_k,nz))
          allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
          allocate(PK_data%log_kh(PK_data%num_k))
          allocate(PK_data%redshifts(nz))
          PK_data%redshifts = CP%Transfer%Redshifts(itf_start:itf_end)
       
          h = CP%H0/100

          do ik=1,MTrans%num_q_trans
                 kh = MTrans%TransferData(Transfer_kh,ik,1)
                 k = kh*h
                 PK_data%log_kh(ik) = log(kh)
                 do itf = 1, nz
                   PK_data%matpower(ik,itf) = &
                    log(MTrans%TransferData(transfer_power_var,ik,itf_start+itf-1)**2*k & 
                                   *pi*twopi*h**3*ScalarPower(k,in))
                 end do
          end do
     
          call MatterPowerdata_getsplines(PK_data)

        end subroutine Transfer_GetMatterPowerData

        subroutine MatterPowerData_Load(PK_data,fname)
          !Loads in kh, P_k from file for one redshiftr and one initial power spectrum
          !Not redshift is not stored in file, so not set correctly
          !Also note that output _matterpower file is already interpolated, so re-interpolating is probs not a good idea

          !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          use AmlUtils
          character(LEN=*) :: fname
          Type(MatterPowerData) :: PK_data
          real(dl)kh, Pk
          integer ik
          integer nz
          

          nz = 1
          call openTxtFile(fname, fileio_unit)
         
          PK_data%num_k = FileLines(fileio_unit)
          PK_Data%num_z = 1

          allocate(PK_data%matpower(PK_data%num_k,nz))
          allocate(PK_data%ddmat(PK_data%num_k,nz))
          allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
          allocate(PK_data%log_kh(PK_data%num_k))
       
          allocate(PK_data%redshifts(nz))
          PK_data%redshifts = 0

          do ik=1,PK_data%num_k
              read (fileio_unit,*) kh, Pk
              PK_data%matpower(ik,1) = log(Pk) 
              PK_data%log_kh(ik) = log(kh)
          end do
     
          call MatterPowerdata_getsplines(PK_data)

        end subroutine MatterPowerData_Load


        subroutine MatterPowerdata_getsplines(PK_data)
          Type(MatterPowerData) :: PK_data
          integer i
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

          do i = 1,PK_Data%num_z
          
           call spline(PK_data%log_kh,PK_data%matpower(1,i),PK_data%num_k,&
                               cllo,clhi,PK_data%ddmat(1,i))
          end do

        end subroutine MatterPowerdata_getsplines
        
        subroutine MatterPowerdata_MakeNonlinear(PK_data)
          Type(MatterPowerData) :: PK_data

          call NonLinear_GetRatios(PK_data)
          PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
          call MatterPowerdata_getsplines(PK_data)

        end subroutine MatterPowerdata_MakeNonlinear

        subroutine MatterPowerdata_Free(PK_data)
          Type(MatterPowerData) :: PK_data
          integer i

          deallocate(PK_data%log_kh,stat=i)
          deallocate(PK_data%matpower,stat=i)
          deallocate(PK_data%ddmat,stat=i)
          deallocate(PK_data%nonlin_ratio,stat=i)
          deallocate(PK_data%redshifts)

        end subroutine MatterPowerdata_Free

        function MatterPowerData_k(PK,  kh, itf) result(outpower)
         !Get matter power spectrum at particular k/h by interpolation
          Type(MatterPowerData) :: PK
          integer, intent(in) :: itf
          real (dl), intent(in) :: kh
          real(dl) :: logk
          integer llo,lhi
          real(dl) outpower, dp
          real(dl) ho,a0,b0
          integer, save :: i_last = 1          
          
           logk = log(kh)
           if (logk < PK%log_kh(1)) then
              dp = (PK%matpower(2,itf) -  PK%matpower(1,itf)) / &
                 ( PK%log_kh(2)-PK%log_kh(1) )
              outpower = PK%matpower(1,itf) + dp*(logk - PK%log_kh(1))
           else if (logk > PK%log_kh(PK%num_k)) then
            !Do dodgy linear extrapolation on assumption accuracy of result won't matter
           
             dp = (PK%matpower(PK%num_k,itf) -  PK%matpower(PK%num_k-1,itf)) / &
                 ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
             outpower = PK%matpower(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
           else 

            llo=min(i_last,PK%num_k)
            do while (PK%log_kh(llo) > logk)
               llo=llo-1
            end do
            do while (PK%log_kh(llo+1)< logk)
               llo=llo+1
            end do
            i_last =llo  
            lhi=llo+1
            ho=PK%log_kh(lhi)-PK%log_kh(llo) 
            a0=(PK%log_kh(lhi)-logk)/ho
            b0=1-a0
              
            outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
                  ((a0**3-a0)* PK%ddmat(llo,itf) &
                       +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6
              
          end if

          outpower = exp(max(-30._dl,outpower))

        end function MatterPowerData_k


        subroutine Transfer_GetMatterPower(MTrans,outpower, itf, in, minkh, dlnkh, npoints)
          !Allows for non-smooth priordial spectra
          !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
          !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
          !in units of (h Mpc^{-1})^3.   
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PK
        
          integer, intent(in) :: itf, in, npoints
          real, intent(out) :: outpower(npoints)
          real, intent(in) :: minkh, dlnkh
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
          integer ik, llo,il,lhi,lastix
          real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
          real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
          

          if (npoints < 2) stop 'Need at least 2 points in Transfer_GetMatterPower'

!         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
!            stop 'Transfer_GetMatterPower: kh out of computed region'
!          end if
          if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
                .and. FeedbackLevel > 0 ) &
                    write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'

          
          if (CP%NonLinear/=NonLinear_None) then
           call Transfer_GetMatterPowerData(MTrans, PK, in, itf)
           call NonLinear_GetRatios(PK)
          end if
           
          h = CP%H0/100
          logmink = log(minkh)
          do ik=1,MTrans%num_q_trans
             kh = MTrans%TransferData(Transfer_kh,ik,itf)
             k = kh*h
             kvals(ik) = log(kh)
             atransfer=MTrans%TransferData(transfer_power_var,ik,itf)
             if (CP%NonLinear/=NonLinear_None) &
                 atransfer = atransfer* PK%nonlin_ratio(ik,1) !only one element, this itf
             matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
                 !Put in power spectrum later: transfer functions should be smooth, initial power may not be                
          end do
             
          call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)

            llo=1
            lastix = npoints + 1
            do il=1, npoints
               xi=logmink + dlnkh*(il-1)
               if (xi < kvals(1)) then
                 outpower(il)=-30.
                 cycle
               end if
               do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
                  llo=llo+1
                  if (llo >= MTrans%num_q_trans) exit
               end do
               if (llo == MTrans%num_q_trans) then
                   lastix = il
                   exit
               end if
               lhi=llo+1
               ho=kvals(lhi)-kvals(llo) 
               a0=(kvals(lhi)-xi)/ho
               b0=(xi-kvals(llo))/ho
              
               outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
                       +(b0**3-b0)*ddmat(lhi))*ho**2/6
              
            end do

            do while (lastix <= npoints)
               !Do linear extrapolation in the log
               !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
               outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
               lastix = lastix+1
            end do

            outpower = exp(max(-30.,outpower))

            do il = 1, npoints
               k = exp(logmink + dlnkh*(il-1))*h
               outpower(il) = outpower(il) * ScalarPower(k,in) 
            end do

          if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_Free(PK)

        end subroutine Transfer_GetMatterPower
     
        subroutine Transfer_Get_sigma8(MTrans, sigr8)
          use MassiveNu
          Type(MatterTransferData) :: MTrans
          integer ik, itf, in
          real(dl) kh, k, h, x, win, delta
          real(dl) lnk, dlnk, lnko
          real(dl) dsig8, dsig8o, sig8, sig8o, powers
          real(dl), intent(IN) :: sigr8
         
          !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
          !of radius sigr8 h^{-1} Mpc
          
         H=CP%h0/100._dl
         do in = 1, CP%InitPower%nn
          do itf=1,CP%Transfer%num_redshifts
            lnko=0
            dsig8o=0
            sig8=0
            sig8o=0
          do ik=1, MTrans%num_q_trans
               kh = MTrans%TransferData(Transfer_kh,ik,itf)
               if (kh==0) cycle
               k = kh*H
               
               delta = k**2*MTrans%TransferData(transfer_power_var,ik,itf)
               !if (CP%NonLinear/=NonLinear_None) delta= delta* MTrans%NonLinearScaling(ik,itf)
               !sigma_8 defined "as though it were linear"

               x= kh *sigr8
               win =3*(sin(x)-x*cos(x))/x**3
               lnk=log(k)
               if (ik==1) then
                  dlnk=0.5_dl 
                 !Approx for 2._dl/(CP%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
                 !Contribution should be very small in any case 
               else
                  dlnk=lnk-lnko
               end if
               powers = ScalarPower(k,in)
               dsig8=(win*delta)**2*powers
               sig8=sig8+(dsig8+dsig8o)*dlnk/2
               dsig8o=dsig8
               lnko=lnk


          end do
    
          MTrans%sigma_8(itf,in) = sqrt(sig8)
          end do
         end do

        end subroutine Transfer_Get_sigma8

        subroutine Transfer_output_Sig8(MTrans)
           Type(MatterTransferData), intent(in) :: MTrans
           
           integer in, j
       
           do in=1, CP%InitPower%nn
            if (CP%InitPower%nn>1)  write(*,*) 'Power spectrum : ', in
            do j = 1, CP%Transfer%num_redshifts
               write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), ' sigma8 (all matter)=', real(MTrans%sigma_8(j,in))
            end do
           end do

         end subroutine Transfer_output_Sig8


        subroutine Transfer_output_Sig8AndNorm(MTrans)
           Type(MatterTransferData), intent(in) :: MTrans
           integer in, j

           do in=1, CP%InitPower%nn
             write(*,*) 'Power spectrum ',in, ' COBE_scale = ',real(COBE_scales(in))
            do j = 1, CP%Transfer%num_redshifts
               write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), ' sigma8(all matter) = ', &
                    real(MTrans%sigma_8(j,in)*sqrt(COBE_scales(in)))
            end do
           end do
                
         end subroutine Transfer_output_Sig8AndNorm


        subroutine Transfer_Allocate(MTrans)
         Type(MatterTransferData) :: MTrans
         integer st

          deallocate(MTrans%q_trans, STAT = st)
          deallocate(MTrans%TransferData, STAT = st)
          deallocate(MTrans%sigma_8, STAT = st)
          allocate(MTrans%q_trans(MTrans%num_q_trans))            
          allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans,CP%Transfer%num_redshifts))  
          allocate(MTrans%sigma_8(CP%Transfer%num_redshifts, CP%InitPower%nn))
           
        end  subroutine Transfer_Allocate

        subroutine Transfer_Free(MTrans)
          Type(MatterTransferData):: MTrans
          integer st

          deallocate(MTrans%q_trans, STAT = st)
          deallocate(MTrans%TransferData, STAT = st)
          deallocate(MTrans%sigma_8, STAT = st)
          nullify(MTrans%q_trans)
          nullify(MTrans%TransferData)
          nullify(MTrans%sigma_8)
          
        end subroutine Transfer_Free

       subroutine Transfer_SetForNonlinearLensing(P)
          Type(TransferParams) :: P
          integer i

          P%kmax = 5*AccuracyBoost
          P%k_per_logint  = 0
          P%num_redshifts =  nint(10*AccuracyBoost)
          if (P%num_redshifts > max_transfer_redshifts) &
                stop 'Transfer_SetForNonlinearLensing: Too many redshifts'
          do i=1,P%num_redshifts
           P%redshifts(i) = real(P%num_redshifts-i)/(P%num_redshifts/10)
          end do

       end subroutine Transfer_SetForNonlinearLensing



        subroutine Transfer_SaveToFiles(MTrans,FileNames)
          use IniFile
          Type(MatterTransferData), intent(in) :: MTrans
          integer i,ik
          character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)

          do i=1, CP%Transfer%num_redshifts
            if (FileNames(i) /= '') then
            open(unit=fileio_unit,file=FileNames(i),form='formatted',status='replace')
             do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                 write(fileio_unit,'(7E14.6)') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
             end do
            close(fileio_unit)
            end if
          end do

          
        end subroutine Transfer_SaveToFiles

        subroutine Transfer_SaveMatterPower(MTrans, FileNames)
          use IniFile
          !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
          Type(MatterTransferData), intent(in) :: MTrans
          character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
          integer itf,in,i
          integer points
          real, dimension(:,:), allocatable :: outpower
          character(LEN=80) fmt
          real minkh,dlnkh
          Type(MatterPowerData) :: PK_data


          write (fmt,*) CP%InitPower%nn+1
          fmt = '('//trim(adjustl(fmt))//'E15.5)'
          do itf=1, CP%Transfer%num_redshifts
            if (FileNames(itf) /= '') then

            
             if (.not. transfer_interp_matterpower ) then
             
             points = MTrans%num_q_trans
             allocate(outpower(points,CP%InitPower%nn))
       
                 do in = 1, CP%InitPower%nn

                   call Transfer_GetMatterPowerData(MTrans, PK_data, in, itf)

                  if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_MakeNonlinear(PK_Data)

                   outpower(:,in) = exp(PK_data%matpower(:,1))
                   call MatterPowerdata_Free(PK_Data)
                 end do

                 open(unit=fileio_unit,file=FileNames(itf),form='formatted',status='replace')
                 do i=1,points
                  write (fileio_unit, fmt) MTrans%TransferData(Transfer_kh,i,1),outpower(i,1:CP%InitPower%nn)
                 end do
                 close(fileio_unit)

             else


             minkh = 1e-4
             dlnkh = 0.02
             points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
!             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
             allocate(outpower(points,CP%InitPower%nn))
             do in = 1, CP%InitPower%nn
              call Transfer_GetMatterPower(MTrans,outpower(1,in), itf, in, minkh,dlnkh, points)
              if (CP%OutputNormalization == outCOBE) then
                 if (allocated(COBE_scales)) then
                  outpower(:,in) = outpower(:,in)*COBE_scales(in)
                 else
                  if (FeedbackLevel>0) write (*,*) 'Cannot COBE normalize - no Cls generated'
                 end if
             end if
             end do
     
             open(unit=fileio_unit,file=FileNames(itf),form='formatted',status='replace')
             do i=1,points
              write (fileio_unit, fmt) minkh*exp((i-1)*dlnkh),outpower(i,1:CP%InitPower%nn)
             end do
             close(fileio_unit)
             
             end if

             deallocate(outpower) 
             
            end if
          end do

        end subroutine Transfer_SaveMatterPower

        end module Transfer


!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module ThermoData
        use ModelData
        implicit none
        private
        integer,parameter :: nthermo=10000
        
        real(dl) tb(nthermo),cs2(nthermo),xe(nthermo)
        real(dl) dcs2(nthermo)
        real(dl) dotmu(nthermo), ddotmu(nthermo)
        real(dl) sdotmu(nthermo),emmu(nthermo)
        real(dl) demmu(nthermo)
        real(dl) dddotmu(nthermo),ddddotmu(nthermo)
        real(dl) tauminn,dlntau,Maxtau
        real(dl), dimension(:), allocatable :: vis,dvis,ddvis,expmmu,dopac, opac
    
        real(dl) :: tight_tau, actual_opt_depth
         !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
         
        public thermo,inithermo,vis,opac,expmmu,dvis,dopac,ddvis, tight_tau,&
               Thermo_OpacityToTime,matter_verydom_tau, ThermoData_Free
       contains

        subroutine thermo(tau,cs2b,opacity, dopacity)
        !Compute unperturbed sound speed squared,
        !and ionization fraction by interpolating pre-computed tables.
        !If requested also get time derivative of opacity
        implicit none
        real(dl) tau,cs2b,opacity
        real(dl), intent(out), optional :: dopacity

        integer i
        real(dl) d
        
        d=log(tau/tauminn)/dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
          cs2b=cs2(1)+(d+i-1)*dcs2(1)
          opacity=dotmu(1)+(d-1)*ddotmu(1)
!!!
           stop 'thermo out of bounds'
        else if (i >= nthermo) then
          cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
          opacity=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
          if (present(dopacity)) then
             dopacity = 0
             stop 'thermo: shouldn''t happen'
           end if
        else
        !Cubic spline interpolation.
          cs2b=cs2(i)+d*(dcs2(i)+d*(3*(cs2(i+1)-cs2(i))  &
              -2*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)  &
              +2*(cs2(i)-cs2(i+1)))))
          opacity=dotmu(i)+d*(ddotmu(i)+d*(3*(dotmu(i+1)-dotmu(i)) &
              -2*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
              +2*(dotmu(i)-dotmu(i+1)))))

         if (present(dopacity)) then
!!!
          dopacity=(ddotmu(i)+d*(dddotmu(i)+d*(3*(ddotmu(i+1)  &
              -ddotmu(i))-2*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
              +dddotmu(i+1)+2*(ddotmu(i)-ddotmu(i+1))))))/(tau*dlntau)

         end if
        end if
        end subroutine thermo
        


       function Thermo_OpacityToTime(opacity)
         real(dl), intent(in) :: opacity
         integer j
         real(dl) Thermo_OpacityToTime
         !Do this the bad slow way for now..
          !The answer is approximate
         j =1
         do while(dotmu(j)> opacity)
            j=j+1
         end do
          
         Thermo_OpacityToTime = exp((j-1)*dlntau)*tauminn

       end function Thermo_OpacityToTime

     subroutine inithermo(taumin,taumax)
!  Compute and save unperturbed baryon temperature and ionization fraction
!  as a function of time.  With nthermo=10000, xe(tau) has a relative 
! accuracy (numerical integration precision) better than 1.e-5.
        use constants
        use precision
        use ModelParams
        use RECFAST
        use MassiveNu
        real(dl) taumin,taumax
   
   
        real(dl) tau01,adot0,a0,a02,x1,x2,barssc,dtau
        real(dl) xe0,tau,a,a2
        real(dl) adot,tg0,ahalf,adothalf,fe,thomc,thomc0,etc,a2t
        real(dl) dtbdla,vfi,cf1,maxvis, vis
        integer ncount,i,j1,j2,iv,ns
        real(dl) spline_data(nthermo)
        real(dl) last_dotmu
        real(dl) dtauda  !diff of tau w.CP%r.t a and integration
        external dtauda
 
        real(dl) a_verydom

         call InitRECFAST(CP%omegab,CP%h0,CP%tcmb,CP%yhe)
          !almost all the time spent here

        Maxtau=taumax
        tight_tau = 0
        actual_opt_depth = 0
        ncount=0
        thomc0=5.0577d-8*CP%tcmb**4
        tauminn=0.05d0*taumin
        dlntau=log(CP%tau0/tauminn)/(nthermo-1)
        last_dotmu = 0

        matter_verydom_tau = 0
        a_verydom = AccuracyBoost*10*(grhog+grhornomass)/(grhoc+grhob)
        if (CP%Num_Nu_massive /= 0) a_verydom=a_verydom*1.5 

!  Initial conditions: assume radiation-dominated universe.
        tau01=tauminn
        adot0=adotrad
        a0=adotrad*tauminn
        a02=a0*a0
!  Assume that any entropy generation occurs before tauminn.
!  This gives wrong temperature before pair annihilation, but
!  the error is harmless.
        tb(1)=CP%tcmb/a0
        xe0=1._dl
        x1=0._dl
        x2=1._dl
        xe(1)=xe0+0.25d0*CP%yhe/(1._dl-CP%yhe)*(x1+2*x2)
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(1))
        cs2(1)=4._dl/3._dl*barssc*tb(1)
        dotmu(1)=xe(1)*akthom/a02
        sdotmu(1)=0
  
          do i=2,nthermo
          tau=tauminn*exp((i-1)*dlntau)
          dtau=tau-tau01
!  Integrate Friedmann equation using inverse trapezoidal rule.
      
          a=a0+adot0*dtau
          a2=a*a

          adot=1/dtauda(a)

          if (matter_verydom_tau ==0 .and. a > a_verydom) then
             matter_verydom_tau = tau  
          end if
          
          a=a0+2._dl*dtau/(1._dl/adot0+1._dl/adot)         
!  Baryon temperature evolution: adiabatic except for Thomson cooling.
!  Use  quadrature solution.
! This is redundant as also calculated in REFCAST, but agrees well before reionization
          tg0=CP%tcmb/a0
          ahalf=0.5d0*(a0+a)
          adothalf=0.5d0*(adot0+adot)
!  fe=number of free electrons divided by total number of free baryon
!  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
!  more accuracy is required (unlikely) then this can be iterated with
!  the solution of the ionization equation.
          fe=(1._dl-CP%yhe)*xe(i-1)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i-1))
          thomc=thomc0*fe/adothalf/ahalf**3
          etc=exp(-thomc*(a-a0))
          a2t=a0*a0*(tb(i-1)-tg0)*etc-CP%tcmb/thomc*(1._dl-etc)
          tb(i)=CP%tcmb/a+a2t/(a*a)
       
! If there is re-ionization, smoothly increase xe to the 
! requested value.
          if (CP%Reion%Reionization .and. tau > CP%ReionHist%tau_start) then
             if(ncount == 0) then
                ncount=i-1
             end if   

            xe(i) = Reionization_xe(a, tau, xe(ncount))
            !print *,1/a-1,xe(i)
            if (CP%AccurateReionization .and. FeedbackLevel > 0) then                         
                dotmu(i)=(xeRECFAST(a) - xe(i))*akthom/a2
                if (last_dotmu /=0) then
                 actual_opt_depth = actual_opt_depth - 2._dl*dtau/(1._dl/dotmu(i)+1._dl/last_dotmu)
                end if
                last_dotmu = dotmu(i) 
            end if
           
          else
            xe(i)=xeRECFAST(a)
          end if 
       
       
!  Baryon sound speed squared (over c**2).
          dtbdla=-2._dl*tb(i)-thomc*adothalf/adot*(a*tb(i)-CP%tcmb)
          barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i))
          cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3._dl)
! Calculation of the visibility function
          dotmu(i)=xe(i)*akthom/a2

          if (tight_tau==0 .and. 1/(tau*dotmu(i)) > 0.005) tight_tau = tau !0.005
           !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)
           
          if (tau < 0.001) then
             sdotmu(i)=0
          else
             sdotmu(i)=sdotmu(i-1)+2._dl*dtau/(1._dl/dotmu(i)+1._dl/dotmu(i-1))
          end if

          a0=a
          tau01=tau
          adot0=adot
          end do !i
                 
          if (CP%Reion%Reionization .and. (xe(nthermo) < 0.999d0)) then
             write(*,*)'Warning: xe at redshift zero is < 1'
             write(*,*) 'Check input parameters an Reionization_xe'
             write(*,*) 'function in the Reionization module'
          end if
   
        do j1=1,nthermo
           if (sdotmu(j1) - sdotmu(nthermo)< -69) then
           emmu(j1)=1.d-30
           else
           emmu(j1)=exp(sdotmu(j1)-sdotmu(nthermo))
           if (.not. CP%AccurateReionization .and. &
               actual_opt_depth==0 .and. xe(j1) < 1e-3) then
              actual_opt_depth = -sdotmu(j1)+sdotmu(nthermo) 
           end if
          end if
        end do  

        if (CP%AccurateReionization .and. FeedbackLevel > 0) then                         
         write(*,'("Reion opt depth      = ",f7.4)') actual_opt_depth
        end if

        iv=0
        vfi=0._dl
! Getting the starting and finishing times for decoupling a+nd time of maximum visibility
        if (ncount == 0) then
           cf1=1._dl
           ns=nthermo
           else
              cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
              ns=ncount
           end if
         maxvis = 0
         do j1=1,ns
           vis = emmu(j1)*dotmu(j1)
           tau = tauminn*exp((j1-1)*dlntau)
           vfi=vfi+vis*cf1*dlntau*tau
           if ((iv == 0).and.(vfi > 1.0d-7/AccuracyBoost)) then  
              taurst=9._dl/10._dl*tau
              iv=1
           elseif (iv == 1) then 
               if (vis > maxvis) then
                maxvis=vis
                tau_maxvis = tau
               end if
               if (vfi > 0.995) then 
                taurend=tau
                iv=2
                exit
               end if
           end if
         end do

           if (iv /= 2) then
             stop 'inithermo: failed to find end of recombination'
!              taurend=1.5d0*(tauminn*exp((ncount-1)*dlntau))
           end if
          

! Calculating the timesteps during recombination.
           
    
           if (CP%WantTensors) then
              dtaurec=min(dtaurec,taurst/160)/AccuracyBoost 
           else
              dtaurec=min(dtaurec,taurst/40)/AccuracyBoost 
           end if
 
           if (CP%Reion%Reionization) taurend=min(taurend,CP%ReionHist%tau_start)

         if (DebugMsgs) then
           write (*,*) 'taurst, taurend = ', taurst, taurend
         end if
  
        matter_verydom_tau = max(matter_verydom_tau,taurend)

        call splini(spline_data,nthermo)
        call splder(cs2,dcs2,nthermo,spline_data)
        call splder(dotmu,ddotmu,nthermo,spline_data)
        call splder(ddotmu,dddotmu,nthermo,spline_data)  
        call splder(dddotmu,ddddotmu,nthermo,spline_data)
        call splder(emmu,demmu,nthermo,spline_data)
   
        call SetTimeSteps

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) 
        do j2=1,TimeSteps%npoints
             call DoThermoSpline(j2,TimeSteps%points(j2))
        end do 
         !$OMP END PARALLEL DO 

        end subroutine inithermo        


        subroutine SetTimeSteps
        real(dl) dtau0
        integer nri0, nstep

         call Ranges_Init(TimeSteps)

         call Ranges_Add_delta(TimeSteps, taurst, taurend, dtaurec)

        ! Calculating the timesteps after recombination
           if (CP%WantTensors) then
              dtau0=max(taurst/40,Maxtau/2000._dl/AccuracyBoost)
           else       
              dtau0=Maxtau/500._dl/AccuracyBoost
             !Don't need this since adding in Limber on small scales
              !  if (CP%DoLensing) dtau0=dtau0/2 
              !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
           end if
    
         call Ranges_Add_delta(TimeSteps,taurend, CP%tau0, dtau0)

         if (CP%Reion%Reionization) then
           
              nri0=int(Reionization_timesteps(CP%ReionHist)*AccuracyBoost) 
                !Steps while reionization going from zero to maximum
              call Ranges_Add(TimeSteps,CP%ReionHist%tau_start,CP%ReionHist%tau_complete,nri0) 

         end if

!Create arrays out of the region information.
        call Ranges_GetArray(TimeSteps)
        nstep = TimeSteps%npoints

        if (allocated(vis)) then
           deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
        end if
        allocate(vis(nstep),dvis(nstep),ddvis(nstep),expmmu(nstep),dopac(nstep),opac(nstep))

        if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Set ',nstep, ' time steps'
    
        end subroutine SetTimeSteps


        subroutine ThermoData_Free
         if (allocated(vis)) then
           deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
         end if
         call Ranges_Free(TimeSteps)

        end subroutine ThermoData_Free

!cccccccccccccc
        subroutine DoThermoSpline(j2,tau)
        integer j2,i
        real(dl) d,ddopac,tau
        
!     Cubic-spline interpolation.
           d=log(tau/tauminn)/dlntau+1._dl
           i=int(d)
      
           d=d-i
           if (i < nthermo) then
          opac(j2)=dotmu(i)+d*(ddotmu(i)+d*(3._dl*(dotmu(i+1)-dotmu(i)) &
              -2._dl*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
              +2._dl*(dotmu(i)-dotmu(i+1)))))
          dopac(j2)=(ddotmu(i)+d*(dddotmu(i)+d*(3._dl*(ddotmu(i+1)  &
              -ddotmu(i))-2._dl*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
              +dddotmu(i+1)+2._dl*(ddotmu(i)-ddotmu(i+1))))))/(tau &
              *dlntau)
          ddopac=(dddotmu(i)+d*(ddddotmu(i)+d*(3._dl*(dddotmu(i+1) &
              -dddotmu(i))-2._dl*ddddotmu(i)-ddddotmu(i+1)  &
              +d*(ddddotmu(i)+ddddotmu(i+1)+2._dl*(dddotmu(i) &
              -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2)) &
              /(tau*dlntau)**2
          expmmu(j2)=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
              -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
              +2._dl*(emmu(i)-emmu(i+1)))))
 
          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3*opac(j2)*dopac(j2)+ddopac)
          else
          opac(j2)=dotmu(nthermo)
          dopac(j2)=ddotmu(nthermo)
          ddopac=dddotmu(nthermo)
          expmmu(j2)=emmu(nthermo)
          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3._dl*opac(j2)*dopac(j2)+ddopac)

          end if
        end subroutine DoThermoSpline

      
 
      end module ThermoData
