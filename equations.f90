! Equations module for dark energy with constant equation of state parameter w
! allowing for perturbations based on a quintessence model
! by Antony Lewis (http://cosmologist.info/)
! This version November 2006.

! Dec 2003, fixed (fatal) bug in tensor neutrino setup
! Changes to tight coupling approximatione
! June 2004, fixed problem with large scale polarized tensors; support for vector modes
! Generate vector modes on their own. The power spectrum is taken from the scalar parameters.
! August 2004, fixed reionization term in lensing potential
! Nov 2004, change massive neutrino l_max to be consistent with massless if light
! Apr 2005, added DoLateRadTruncation option
! June 2006, added support for arbitary neutrino mass splittings
! Nov 2006, tweak to high_precision transfer function accuracy at lowish k

       module LambdaGeneral
         use precision
         implicit none
          
         real(dl)  :: w_lam = -1 !p/rho for the dark energy (assumed constant) 
         real(dl) :: cs2_lam = 1_dl 
          !comoving sound speed. Always exactly 1 for quintessence 
          !(otherwise assumed constant, though this is almost certainly unrealistic)

         logical :: w_perturb = .true.

       end module LambdaGeneral


    
!Return OmegaK - modify this if you add extra fluid components
        function GetOmegak()
        use precision
        use ModelParams
        real(dl)  GetOmegak
         GetOmegak = 1 - (CP%omegab+CP%omegac+CP%omegav+CP%omegan) 
          
        end function GetOmegak
  
  
       subroutine init_background
         !This is only called once per model, and is a good point to do any extra initialization.
         !It is called before first call to dtauda, but after
         !massive neutrinos are initialized and after GetOmegak
       end  subroutine init_background


!Background evolution
        function dtauda(a)
         !get d tau / d a
        use precision
        use ModelParams
        use MassiveNu
        use LambdaGeneral
        implicit none
        real(dl) dtauda
        real(dl), intent(IN) :: a
        real(dl) rhonu,grhoa2, a2
        integer nu_i

        a2=a**2

!  8*pi*G*rho*a**4.
        grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
         if (w_lam == -1._dl) then
           grhoa2=grhoa2+grhov*a2**2
         else
           grhoa2=grhoa2+grhov*a**(1-3*w_lam)
         end if
        if (CP%Num_Nu_massive /= 0) then
!Get massive neutrino density relative to massless
           do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grhoa2=grhoa2+rhonu*grhormass(nu_i)
           end do
        end if

        dtauda=sqrt(3/grhoa2)
     
        end function dtauda
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!Gauge-dependent perturbation equations

        module GaugeInterface
        use precision
        use ModelParams
        use MassiveNu
        use LambdaGeneral
        implicit none
        public

        !Description of this file. Change if you make modifications.
        character(LEN=*), parameter :: Eqns_name = 'gauge_inv'

        logical :: DoTensorNeutrinos = .false.
        
        logical :: DoLateRadTruncation = .true.
            !if true, use approx to radition perturbations after matter domination on
            !small scales, saving evolution of irrelevant osciallatory multipole equations

        real(dl) :: Magnetic = 0._dl
            !Vector mode anisotropic stress in units of rho_gamma
        real(dl) :: vec_sig0 = 1._dl
            !Vector mode shear      
        integer, parameter :: max_l_evolve = 1024 !Maximum l we are ever likely to propagate

        !Supported scalar initial condition flags
         integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
         initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0
         integer, parameter :: initial_nummodes =  initial_iso_neutrino_vel

        type EvolutionVars
            real(dl) q, q2
            real(dl) k_buf,k2_buf ! set in initial

            integer w_ix !Index of two quintessence equations

            integer q_ix !index into q_evolve array that gives the value q
            logical TransferOnly

    !       nvar  - number of scalar (tensor) equations for this k       
            integer nvar,nvart, nvarv

           !Max_l for the various hierarchies
            integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
            integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt
            integer lmaxnrv, lmaxv, lmaxpolv

            integer polind  !index into scalar array of polarization hierarchy

    !array indices for massive neutrino equations
            integer iq0,iq1,iq2

    !array index for tensor massive neutrino equations
            integer iqt
 
    !Initial values for massive neutrino v*3 variables calculated when switching 
    !to non-relativistic approx
            real(dl) G11(max_nu),G30(max_nu)
            real(dl) w_nu !equation of state parameter for massive neutrinos
    !True when using non-relativistic approximation
            logical MassiveNuApprox

    !Massive neutrino scheme being used at the moment        
            integer NuMethod

    !Tru when using tight-coupling approximation (required for stability at early times)

            logical TightCoupling
  
    !Numer of scalar equations we are propagating
            integer ScalEqsToPropagate
    !beta > l for closed models 
            integer FirstZerolForBeta
    !Tensor vars
            real(dl) tenspigdot, aux_buf

            real(dl) pig !For tight coupling

            real(dl) poltruncfac
        
            logical no_rad_multpoles 

    !Buffers for non-flat vars
            real(dl) Kf(max_l_evolve),Kft(max_l_evolve)      

        end type EvolutionVars

!precalculated arrays
        real(dl) polfac(max_l_evolve),tensfac(max_l_evolve),tensfacpol(max_l_evolve), &
               denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve) 
        
       real(dl), parameter :: ep0=1.0d-2 

       real(dl) epsw
    
       integer debug

       contains


        subroutine GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
         type(EvolutionVars) EV
         real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
         integer ind

          if (CP%flat) then
            call dverk(EV,EV%ScalEqsToPropagate,fderivs,tau,y,tauend,tol1,ind,c,EV%nvar,w)
          else
            call dverk(EV,EV%ScalEqsToPropagate,derivs,tau,y,tauend,tol1,ind,c,EV%nvar,w)
          end if
        end subroutine GaugeInterface_ScalEv

         subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
         use ThermoData
         type(EvolutionVars) EV
         real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
         integer ind
         real(dl) ep, tau_switch, tau_check
         real(dl) cs2, opacity

         !Evolve equations from tau to tauend, performing switch-off of
         !tight coupling if necessary.
         !In principle the tight coupling evolution routine could be different
         !which could probably be used to improve the speed a bit
         
         !It's possible due to instabilities that changing later is actually more
         !accurate, so cannot expect accuracy to vary monotonically with the switch

         if (EV%TightCoupling) then
    
             tau_switch = tight_tau !when 1/(opacity*tau) = 0.01ish
    
            !The numbers here are a bit of guesswork
            !The high k increase saves time for very small loss of accuracy
            !The lower k ones are more delicate. Nead to avoid instabilities at same time
            !as ensuring tight coupling is accurate enough
             if (EV%k_buf > epsw) then
               if (EV%k_buf > epsw*5) then
                ep=ep0*5/AccuracyBoost
               else
                ep=ep0
               end if
             else
               ep=ep0 
             end if

            !Check the k/opacity criterion
             tau_check = min(tauend, tau_switch)   
             call thermo(tau_check,cs2,opacity)
             if (EV%k_buf/opacity > ep) then
                !so need to switch in this time interval
                 tau_switch = Thermo_OpacityToTime(EV%k_buf/ep) 
             end if

             if (tauend > tau_switch) then 
           
                 if (tau_switch > tau) then
                  call GaugeInterface_ScalEv(EV, y, tau,tau_switch,tol1,ind,c,w)
                 end if
            
               !Set up variables with their tight coupling values
                 y(8) = EV%pig
                 y(9) = 3./7*y(8)*EV%k_buf/opacity
                 y(EV%polind+2) = EV%pig/4   
                 y(EV%polind+3) =y(9)/4      

                 EV%TightCoupling = .false.

            end if
         end if

         !Turn off radiation hierarchies at late time where slow and 
         !not needed. Must be matter domainted and well inside the horizon
         !Not tested for non-adiabatic         
         tau_switch=max(15/EV%k_buf,matter_verydom_tau)
         if (.not. EV%no_rad_multpoles .and. tauend > tau_switch .and. DoLateRadTruncation &
           .and. (.not.CP%WantCls .or. EV%k_buf>0.02*AccuracyBoost) &
            .and. (CP%Scalar_initial_condition==initial_adiabatic .or. &
               CP%Scalar_initial_condition==initial_vector .and. &
                 all(CP%InitialConditionVector(2:initial_nummodes) ==0))) then

                 if (tau_switch > tau) then
                  call GaugeInterface_ScalEv(EV, y, tau,tau_switch,tol1,ind,c,w)
                 end if
                 EV%no_rad_multpoles = .true.
                 y(6:EV%polind+EV%lmaxgpol)=0
                if (CP%Num_Nu_massive == 0) then
                 EV%ScalEqsToPropagate=5
                 if (w_lam /= -1 .and. w_Perturb) then
                  !actually DE perturbations probably irrelvant and could set to zero too
                    y(6)=y(EV%w_ix)
                    y(7)=y(EV%w_ix+1)
                    EV%w_ix = 6
                    EV%ScalEqsToPropagate=7                    
                 end if
                end if
                  
         end if       
         
         call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
        
        end subroutine GaugeInterface_EvolveScal

 
        subroutine GaugeInterface_Init
          !Precompute various arrays and other things independent of wavenumber
          integer j

          epsw = 100/CP%tau0
         
          if (CP%WantScalars) then
            do j=2,max_l_evolve
              polfac(j)=real((j+3)*(j-1),dl)/(j+1)
            end do               
          end if
      
          if (CP%WantVectors) then
            do j=2,max_l_evolve
             vecfac(j)=real((j+2),dl)/(j+1)
             vecfacpol(j)=real((j+3)*j,dl)*(j-1)*vecfac(j)/(j+1)**2
           end do
   
          end if

          if (CP%WantTensors) then
           do j=2,max_l_evolve
           tensfac(j)=real((j+3)*(j-1),dl)/(j+1)
           tensfacpol(j)=tensfac(j)**2/(j+1)
          end do
          end if

         do j=1,max_l_evolve
           denl(j)=1._dl/(2*j+1)
         end do     
       
        end subroutine GaugeInterface_Init


        subroutine GetNumEqns(EV)
          use MassiveNu
          !Set the numer of equations in each hierarchy, and get total number of equations for this k
          type(EvolutionVars) EV
          real(dl) scal
        
         if (CP%Num_Nu_massive == 0) then
            EV%lmaxnu=0
         else 
             EV%NuMethod = CP%MassiveNuMethod
             if (EV%NuMethod == Nu_Best) then
              if (all(nu_masses(1:CP%Nu_mass_eigenstates) < 1800._dl/AccuracyBoost) .and. &
                  .not. CP%WantTransfer .and. .not. CP%DoLensing) then
                !If light then approx is very good for CMB
                EV%NuMethod = Nu_approx
               else
                EV%NuMethod = Nu_Trunc
               end if
             end if
            !l_max for massive neutrinos
            !if relativistic at recombination set to massless lmaxnr below
             EV%lmaxnu=max(3,nint(6*lAccuracyBoost))   
            if (CP%Transfer%high_precision) EV%lmaxnu=nint(25*lAccuracyBoost)
         end if

        if (CP%WantScalars) then
         EV%lmaxg  = max(nint(8*lAccuracyBoost),3)  
         EV%lmaxnr = max(nint(14*lAccuracyBoost),3)  
            !need quite high to avoid ~ 0.5% systematics at l>1000
         EV%lmaxgpol = EV%lmaxg  
         if (.not.CP%AccuratePolarization) EV%lmaxgpol=max(nint(4*lAccuracyBoost),3)
     
         if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            if (CP%AccuratePolarization) scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*lAccuracyBoost))
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*lAccuracyBoost)) 
            if (CP%AccurateReionization) then
             EV%lmaxg=EV%lmaxg*4
             EV%lmaxgpol=EV%lmaxgpol*2
            end if
         end if                  
         if (CP%Transfer%high_precision) then
           EV%lmaxnr=max(nint(25*lAccuracyBoost),3) 
           EV%lmaxg=max(EV%lmaxg,nint(min(8,nint(sqrt(scal) *600 * EV%q))*lAccuracyBoost)) 
         end if
         EV%nvar=5+ (EV%lmaxg+1) + EV%lmaxgpol-1 +(EV%lmaxnr+1) 
         if (w_lam /= -1 .and. w_Perturb) then
            EV%w_ix = EV%nvar+1
            EV%nvar=EV%nvar+2 
         else
            EV%w_ix=0
         end if

         if (CP%Num_Nu_massive /= 0) then
           if (any(nu_masses(1:CP%Nu_mass_eigenstates) < 7000)) then              
               EV%lmaxnu=EV%lmaxnr 
           end if
 
           if (EV%NuMethod == Nu_approx) then
              EV%iq0=EV%nvar + 1 
              EV%nvar= EV%nvar+(EV%lmaxnu+1)*CP%Nu_mass_eigenstates
           else
            EV%iq0=EV%nvar + 1 
            EV%iq1=EV%iq0+nqmax
            EV%iq2=EV%iq1+nqmax
            EV%nvar=EV%nvar+ nqmax*(EV%lmaxnu+1)*CP%Nu_mass_eigenstates
           endif
         end if
        
         EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
         if (EV%MaxlNeeded > max_l_evolve) stop 'Need to increase max_l_evolve'
        
         EV%poltruncfac=real(EV%lmaxgpol,dl)/(EV%lmaxgpol-2)
         EV%polind = 6+EV%lmaxnr+EV%lmaxg
         EV%lmaxt=0
         
        else
          EV%nvar=0
        end if
    
       if (CP%WantTensors) then
          EV%lmaxt=max(3,nint(8*lAccuracyBoost))
          EV%lmaxpolt = max(3,nint(5*lAccuracyBoost)) 
          if (EV%q < 0.05) then
            if (.not. CP%AccuratePolarization) then
             scal  = 1
             EV%lmaxt=max(3,nint(min(8,nint(scal *150 * EV%q))*lAccuracyBoost))
             EV%lmaxpolt=max(3,nint(min(5,nint( scal*150 * EV%q))*lAccuracyBoost))
            end if
          end if      
          EV%nvart=(EV%lmaxt-1)+(EV%lmaxpolt-1)*2+3
          if (DoTensorNeutrinos) then
           
            EV%lmaxnrt=nint(6*lAccuracyBoost)
            EV%lmaxnut=EV%lmaxnrt
            EV%nvart=EV%nvart+EV%lmaxnrt-1
             if (CP%Num_Nu_massive /= 0 ) then
                 EV%iqt=EV%nvart + 1                  
                 EV%nvart=EV%nvart+ nqmax*(EV%lmaxnut-1)*CP%Nu_mass_eigenstates         
              end if 
          end if
        else
          EV%nvart=0
        end if       


       if (CP%WantVectors) then
          EV%lmaxv=max(10,nint(8*lAccuracyBoost))
          EV%lmaxpolv = max(5,nint(5*lAccuracyBoost)) 
          
           EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3
          
           EV%lmaxnrv=nint(30*lAccuracyBoost)
  
           EV%nvarv=EV%nvarv+EV%lmaxnrv
           if (CP%Num_Nu_massive /= 0 ) then
                  stop 'massive neutrinos not supported for vector modes'
           end if 
        else
          EV%nvarv=0
        end if       

        end subroutine GetNumEqns

!cccccccccccccccccccccccccccccccccc
        subroutine SwitchToMassiveNuApprox(EV,y)
!When the neutrinos are no longer highly relativistic we use a truncated
!energy-integrated hierarchy going up to third order in velocity dispersion
           type(EvolutionVars) EV

        real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
        real(dl) shearnu, qnu
        real(dl) y(EV%nvar)
        integer nu_i, off_ix, qoff_ix
   
        a=y(1)
        a2=a*a

        do nu_i = 1, CP%Nu_mass_eigenstates
          

           off_ix = (nu_i-1)*4
           qoff_ix = (nu_i-1)*nqmax*(EV%lmaxnu+1)

           !Get density and pressure as ratio to massles by interpolation from table
           call Nu_background(a*nu_masses(nu_i),rhonu,pnu)
           !Integrate over q

           call Nu_Integrate(a*nu_masses(nu_i),clxnu,qnu,dpnu,shearnu, &
                 y(EV%iq0+qoff_ix),y(EV%iq1+qoff_ix),y(EV%iq2+qoff_ix))
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !Could save time by only calculating shearnu in output
           dpnu=dpnu/rhonu
           qnu=qnu/rhonu
           clxnu = clxnu/rhonu
           pinu=1.5_dl*shearnu/rhonu                     
    
 
           EV%MassiveNuApprox=.true.
        
           y(EV%iq0+off_ix)=clxnu
           y(EV%iq0+off_ix+1)=dpnu
           y(EV%iq0+off_ix+2)=qnu
           y(EV%iq0+off_ix+3)=pinu

           call Nu_Intvsq(a*nu_masses(nu_i),EV%G11(nu_i),EV%G30(nu_i),y(EV%iq1+qoff_ix),y(EV%iq0+qoff_ix+3*nqmax))
!Analytic solution for higher moments, proportional to a^{-3}
           EV%G11(nu_i)=EV%G11(nu_i)*a2*a/rhonu  
           EV%G30(nu_i)=EV%G30(nu_i)*a2*a/rhonu  

         end do
        
          EV%ScalEqsToPropagate=(EV%iq0-1)+4*CP%Nu_mass_eigenstates

        end  subroutine SwitchToMassiveNuApprox

       subroutine MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum)
        implicit none
        type(EvolutionVars) EV
        real(dl) :: y(EV%nvar), yprime(EV%nvar),a, grho,gpres,dgrho,dgq,dgpi, gdpi_diff,pidot_sum
          !grho = a^2 kappa rho
          !gpres = a^2 kappa p
          !dgrho = a^2 kappa \delta\rho
          !dgp =  a^2 kappa \delta p
          !dgq = a^2 kappa q (heat flux)
          !dgpi = a^2 kappa pi (anisotropic stress) 
          !dgpi_diff = a^2 kappa (3*p -rho)*pi

        integer nu_i, off_ix
        real(dl) pinudot,grhormass_t, rhonu, pnu, shearnu, rhonudot
        real(dl) shearnudot, adotoa, grhonu_t,gpnu_t
        real(dl) clxnu, qnu, pinu, dpnu
        real(dl) dtauda

         EV%w_nu = 0
         do nu_i = 1, CP%Nu_mass_eigenstates

           grhormass_t=grhormass(nu_i)/a**2

           !Get density and pressure as ratio to massless by interpolation from table
           call Nu_background(a*nu_masses(nu_i),rhonu,pnu)
           
           if (EV%NuMethod == Nu_approx) then

             off_ix = EV%iq0+(nu_i-1)*(EV%lmaxnu+1)
             clxnu=y(off_ix)
             qnu=y(off_ix+1)
             pinu=y(off_ix+2)
             pinudot=yprime(off_ix+2)
    
           else
           

           if (EV%MassiveNuApprox) then    
              off_ix = (nu_i-1)*4
              clxnu=y(EV%iq0+off_ix)
              !dpnu = y(EV%iq0+1+off_ix)
              qnu=y(EV%iq0+off_ix+2)
              pinu=y(EV%iq0+off_ix+3)
              pinudot=yprime(EV%iq0+off_ix+3)
        
           else
            off_ix = (nu_i-1)*nqmax*(EV%lmaxnu+1) 
            !Integrate over q

            call Nu_Integrate(a*nu_masses(nu_i),clxnu,qnu,dpnu,shearnu, &
                  y(EV%iq0+off_ix),y(EV%iq1+off_ix),y(EV%iq2+off_ix))
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !dpnu=dpnu/rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
            pinu=1.5_dl*shearnu/rhonu       
            adotoa = 1/(a*dtauda(a))
            call Nu_derivs(a*nu_masses(nu_i),adotoa,rhonu,rhonudot,shearnudot, &
                                y(EV%iq2+off_ix),yprime(EV%iq2+off_ix))
            pinudot=1.5_dl*shearnudot/rhonu - rhonudot/rhonu*pinu    
           endif
          end if

          grhonu_t=grhormass_t*rhonu
          gpnu_t=grhormass_t*pnu
          
          grho = grho  + grhonu_t
          gpres= gpres + gpnu_t

          EV%w_nu = max(EV%w_nu,pnu/rhonu)

          dgrho= dgrho + grhonu_t*clxnu
          dgq  = dgq   + grhonu_t*qnu
          dgpi = dgpi  + grhonu_t*pinu
          gdpi_diff = gdpi_diff + pinu*(3*gpnu_t-grhonu_t)
          pidot_sum = pidot_sum + grhonu_t*pinudot

        end do

     end subroutine MassiveNuVarsOut

           
     subroutine MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq, wnu_arr)
        implicit none
        type(EvolutionVars) EV
        real(dl) :: y(EV%nvar), a, grho,gpres,dgrho,dgq
        real(dl), intent(out), optional :: wnu_arr(max_nu)
          !grho = a^2 kappa rho
          !gpres = a^2 kappa p
          !dgrho = a^2 kappa \delta\rho
          !dgp =  a^2 kappa \delta p
          !dgq = a^2 kappa q (heat flux)
        integer nu_i, off_ix
        real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

         do nu_i = 1, CP%Nu_mass_eigenstates

          grhormass_t=grhormass(nu_i)/a**2

          !Get density and pressure as ratio to massless by interpolation from table
          call Nu_background(a*nu_masses(nu_i),rhonu,pnu)
 
          if (EV%NuMethod == Nu_approx) then
             off_ix = EV%iq0+(nu_i-1)*(EV%lmaxnu+1)
             clxnu=y(off_ix)
             qnu=y(off_ix+1)
          else

           if (EV%MassiveNuApprox) then    
              off_ix = (nu_i-1)*4
              clxnu=y(EV%iq0+off_ix)
              qnu=y(EV%iq0+off_ix+2)
           else
            off_ix = (nu_i-1)*nqmax*(EV%lmaxnu+1) 
            !Integrate over q
            call Nu_Integrate01(a*nu_masses(nu_i),clxnu,qnu,y(EV%iq0+off_ix),y(EV%iq1+off_ix))
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
           endif

          end if

          grhonu_t=grhormass_t*rhonu
          gpnu_t=grhormass_t*pnu
          
          grho = grho  + grhonu_t
          gpres= gpres + gpnu_t
          dgrho= dgrho + grhonu_t*clxnu
          dgq  = dgq   + grhonu_t*qnu

          if (present(wnu_arr)) then
           wnu_arr(nu_i) =pnu/rhonu
          end if

        end do

     end subroutine MassiveNuVars


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine output(EV,y, n,j,tau,sources)
        use ThermoData
        use lvalues
        use ModelData 
        implicit none
        integer n,j
        type(EvolutionVars) EV
        real(dl), target :: y(EV%nvar),yprime(EV%nvar)
        real(dl), dimension(:),pointer :: ypol,ypolprime
        
        real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
        real(dl) qgdot,pigdot,pirdot,vbdot,dgrho
        real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir

        real(dl) tau,x,divfac
        real(dl) dgpi_diff, pidot_sum
          !dgpi_diff = sum (3*p_nu -rho_nu)*pi_nu

        real(dl) k,k2  ,adotoa, grho, gpres,etak,phi,dgpi
        real(dl) clxq, vq, diff_rhopi
        real(dl) sources(CTransScal%NumSources)
        real(dl) ISW

        !real(dl) t4,t92
 
        yprime = 0
        if (CP%flat) then
            call fderivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
        else
            call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)        
        end if
        
        ypolprime => yprime(EV%polind+1:)
        ypol => y(EV%polind+1:)

        k=EV%k_buf
        k2=EV%k2_buf
     
        a   =y(1)
        a2  =a*a
        etak=y(2)
        clxc=y(3)
        clxb=y(4)
        vb  =y(5)
        vbdot =yprime(5)

!  Compute expansion rate from: grho 8*pi*rho*a**2

        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*a**(-1-3*w_lam)
    
        if (EV%no_rad_multpoles) then
            clxg=2*(grhoc_t*clxc+grhob_t*clxb)/3/k**2
            clxr=clxg
            qg= clxg*k/sqrt((grhoc_t+grhob_t)/3)*(2/3._dl)
            qr=qg
            pig=0
            pir=0
            pigdot=0
            pirdot=0
            qgdot=k/3*clxg
        else
            if (EV%TightCoupling) then
             y(8) = EV%pig
             ypol(2) = EV%pig/4
            end if 
            clxg=y(6)
            qg  =y(7)
            pig =y(8)
            clxr=y(7+EV%lmaxg)
            qr  =y(8+EV%lmaxg)
            pir =y(9+EV%lmaxg)
            pigdot=yprime(8)
            pirdot=yprime(EV%lmaxg+9)
            qgdot =yprime(7)
 
        end if

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres=(grhog_t+grhor_t)/3+grhov_t*w_lam

!  8*pi*a*a*SUM[rho_i*clx_i]
        dgrho=grhob_t*clxb+grhoc_t*clxc + grhog_t*clxg+grhor_t*clxr

!  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        dgq=grhob_t*vb+grhog_t*qg+grhor_t*qr

        dgpi = grhor_t*pir + grhog_t*pig

        dgpi_diff = 0
        pidot_sum = 0

        if (CP%Num_Nu_Massive /= 0) then
         call MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum)
        end if

        adotoa=sqrt((grho+grhok)/3)

        if (w_lam /= -1 .and. w_Perturb) then
          
           clxq=y(EV%w_ix)
           vq=y(EV%w_ix+1) 
           dgrho=dgrho + clxq*grhov_t
           dgq = dgq + vq*grhov_t*(1+w_lam)
         
        end if

!  Get sigma (shear) and z from the constraints
!  have to get z from eta for numerical stability       
        z=(0.5_dl*dgrho/k + etak)/adotoa 
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)

        polter = 0.1_dl*pig+9._dl/15._dl*ypol(2)

        if (CP%flat) then
        x=k*(CP%tau0-tau)
        divfac=x*x    
        else   
        x=(CP%tau0-tau)/CP%r
        divfac=(CP%r*rofChi(x))**2*k2 
        end if


        if (EV%TightCoupling) then
          pigdot = -dopac(j)/opac(j)*pig + 32._dl/45*k/opac(j)*(-2*adotoa*sigma  &
                 +etak/EV%Kf(1)-  dgpi/k +vbdot )
          ypolprime(2)= pigdot/4
        end if

        pidot_sum =  pidot_sum + grhog_t*pigdot + grhor_t*pirdot
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff )*adotoa

!Maple's fortran output - see scal_eqs.map
!2phi' term (\phi' + \psi' in Newtonian gauge)
        ISW = (4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
              -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*gpres+5.D0*grho)*sigma/k/3.D0 &
              -2.D0/k*adotoa/EV%Kf(1)*etak)*expmmu(j)

!The rest
        sources(1)= ISW +  ((-9.D0/160.D0*pig-27.D0/80.D0*ypol(2))/k**2*opac(j)+(11.D0/10.D0*sigma- &
    3.D0/8.D0*EV%Kf(2)*ypol(3)+vb-9.D0/80.D0*EV%Kf(2)*y(9)+3.D0/40.D0*qg)/k-(- &
    180.D0*ypolprime(2)-30.D0*pigdot)/k**2/160.D0)*dvis(j)+(-(9.D0*pigdot+ &
    54.D0*ypolprime(2))/k**2*opac(j)/160.D0+pig/16.D0+clxg/4.D0+3.D0/8.D0*ypol(2)+(- &
    21.D0/5.D0*adotoa*sigma-3.D0/8.D0*EV%Kf(2)*ypolprime(3)+vbdot+3.D0/40.D0*qgdot- &
    9.D0/80.D0*EV%Kf(2)*yprime(9))/k+(-9.D0/160.D0*dopac(j)*pig-21.D0/10.D0*dgpi-27.D0/ &
    80.D0*dopac(j)*ypol(2))/k**2)*vis(j)+(3.D0/16.D0*ddvis(j)*pig+9.D0/ &
    8.D0*ddvis(j)*ypol(2))/k**2+21.D0/10.D0/k/EV%Kf(1)*vis(j)*etak   


!Equivalent result
!    t4 = 1.D0/adotoa
!    t92 = k**2
!   sources(1) = (4.D0/3.D0*EV%Kf(1)*expmmu(j)*sigma+2.D0/3.D0*(-sigma-t4*etak)*expmmu(j))*k+ &
!       (3.D0/8.D0*ypol(2)+pig/16.D0+clxg/4.D0)*vis(j)
!    sources(1) = sources(1)-t4*expmmu(j)*dgrho/3.D0+((11.D0/10.D0*sigma- &
!         3.D0/8.D0*EV%Kf(2)*ypol(3)+vb+ 3.D0/40.D0*qg-9.D0/80.D0*EV%Kf(2)*y(9))*dvis(j)+(5.D0/3.D0*grho+ &
!        gpres)*sigma*expmmu(j)+(-2.D0*adotoa*etak*expmmu(j)+21.D0/10.D0*etak*vis(j))/ &
!        EV%Kf(1)+(vbdot-3.D0/8.D0*EV%Kf(2)*ypolprime(3)+3.D0/40.D0*qgdot-21.D0/ &
!        5.D0*sigma*adotoa-9.D0/80.D0*EV%Kf(2)*yprime(9))*vis(j))/k+(((-9.D0/160.D0*pigdot- &
!        27.D0/80.D0*ypolprime(2))*opac(j)-21.D0/10.D0*dgpi -27.D0/80.D0*dopac(j)*ypol(2) &
!        -9.D0/160.D0*dopac(j)*pig)*vis(j) - diff_rhopi*expmmu(j)+((-27.D0/80.D0*ypol(2)-9.D0/ &
!        160.D0*pig)*opac(j)+3.D0/16.D0*pigdot+9.D0/8.D0*ypolprime(2))*dvis(j)+9.D0/ &
!        8.D0*ddvis(j)*ypol(2)+3.D0/16.D0*ddvis(j)*pig)/t92


      if (x > 0._dl) then
         !E polarization source
           sources(2)=vis(j)*polter*(15._dl/8._dl)/divfac 
               !factor of four because no 1/16 later
        else
           sources(2)=0
        end if
 
      if (CTransScal%NumSources > 2) then
         !Get lensing sources
         !Can modify this here if you want to get power spectra for other tracer
       if (tau>taurend .and. CP%tau0-tau > 0.1_dl) then
         
         !phi_lens = Phi - 1/2 kappa (a/k)^2 sum_i rho_i pi_i
         !Neglect pi contributions because not accurate at late time anyway
         phi = -(dgrho +3*dgq*adotoa/k)/(k2*EV%Kf(1)*2) 
            ! - (grhor_t*pir + grhog_t*pig+ pinu*gpnu_t)/k2
         
         sources(3) = -2*phi*f_K(tau-tau_maxvis)/(f_K(CP%tau0-tau_maxvis)*f_K(CP%tau0-tau))
 
!         sources(3) = -2*phi*(tau-tau_maxvis)/((CP%tau0-tau_maxvis)*(CP%tau0-tau))
          !We include the lensing factor of two here
       else
         sources(3) = 0
       end if
      end if
     end subroutine output


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine outputt(EV,yt,n,j,tau,dt,dte,dtb)
!calculate the tensor sources for open and closed case
        use ThermoData

        implicit none
        integer j,n
        type(EvolutionVars) :: EV
        real(dl), target :: yt(n), ytprime(n)
        real(dl) tau,dt,dte,dtb,x,polterdot,polterddot,prefac
        real(dl) pig, pigdot, aux, polter, shear
        real(dl) sinhxr,cothxor
        real(dl) k,k2
        real(dl), dimension(:),pointer :: E,Bprime,Eprime

        if (CP%flat) then
           call fderivst(EV,EV%nvart,tau,yt,ytprime)
        else
           call  derivst(EV,EV%nvart,tau,yt,ytprime)
        end if
      
        k2=EV%k2_buf
        k=EV%k_buf 
        pigdot=EV%tenspigdot
        pig=yt(4)
        aux=EV%aux_buf  
        shear = yt(3)

        x=(CP%tau0-tau)/CP%r
       
        sinhxr=rofChi(x)*CP%r
    
        if (EV%q*sinhxr > 1.e-8_dl) then  

        prefac=sqrt(EV%q2*CP%r*CP%r-CP%Ksign)
        cothxor=cosfunc(x)/sinhxr
        E => yt(EV%lmaxt+3-1:)   
        Eprime=> ytprime(EV%lmaxt+3-1:) 
      
        Bprime => Eprime(EV%lmaxpolt:)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
        polterddot = 9._dl/15._dl*(-dopac(j)*(E(2)-polter)-opac(j)*(  &
                   Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - &
                   5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
                   +0.1_dl*(k*(-ytprime(5)*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(3)) - &
                   dopac(j)*(pig - polter) - opac(j)*(pigdot-polterdot))

        dt=(shear*expmmu(j) + (15._dl/8._dl)*polter*vis(j)/k)*CP%r/sinhxr**2/prefac 
   
        dte=CP%r*15._dl/8._dl/k/prefac* &
            ((ddvis(j)*polter + 2._dl*dvis(j)*polterdot + vis(j)*polterddot)  &
              + 4._dl*cothxor*(dvis(j)*polter + vis(j)*polterdot) - &
                   vis(j)*polter*(k2 -6*cothxor**2))
      
        dtb=15._dl/4._dl*EV%q*CP%r/k/prefac*(vis(j)*(2._dl*cothxor*polter + polterdot) + dvis(j)*polter)
   
        else
        dt=0._dl
        dte=0._dl
        dtb=0._dl
        end if

        end subroutine outputt

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine outputv(EV,yv,n,j,tau,dt,dte,dtb)
!calculate the vector sources 
        use ThermoData

        implicit none
        integer j,n
        type(EvolutionVars) :: EV
        real(dl), target :: yv(n), yvprime(n)
        real(dl) tau,dt,dte,dtb,x,polterdot
        real(dl) vb,qg, pig, polter, sigma
        real(dl) k,k2
        real(dl), dimension(:),pointer :: E,Eprime

        call fderivsv(EV,EV%nvarv,tau,yv,yvprime)
      
        k2=EV%k2_buf
        k=EV%k_buf 
        sigma = yv(2)
        vb  = yv(3)
        qg  = yv(4)
        pig = yv(5)

  
        x=(CP%tau0-tau)*k
              
        if (x > 1.e-8_dl) then  

        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:) 
      
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)
 
        if (yv(1) < 1e-3) then
         dt = 1
        else
         dt =0
        end if   
        dt= (4*(vb+sigma)*vis(j) + 15._dl/2/k*( vis(j)*polterdot + dvis(j)*polter) &
               + 4*(expmmu(j)*yvprime(2)) )/x 
   
        dte= 15._dl/2*2*polter/x**2*vis(j) + 15._dl/2/k*(dvis(j)*polter + vis(j)*polterdot)/x 
      
        dtb= -15._dl/2*polter/x*vis(j)

        else
         dt=0
         dte=0
         dtb=0
        end if

        end subroutine outputv


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initial(EV,y, tau)
!  Initial conditions.
        implicit none
        type(EvolutionVars) EV
        real(dl) y(EV%nvar)
        real(dl) Rp15,tau,x,x2,x3,om,omtau, &
             Rc,Rb,Rv,Rg,grhonu,chi
        real(dl) k,k2 
        real(dl) a,a2, iqg, rhomass
        integer l,i, nu_i, off_ix
        integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
                i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxq=11,i_vq=12
        integer, parameter :: i_max = i_vq
        real(dl) initv(6,1:i_max), initvec(1:i_max)

           
        EV%ScalEqsToPropagate=EV%nvar
        if (CP%flat) then
         EV%k_buf=EV%q
         EV%k2_buf=EV%q2     
         EV%Kf(1:EV%MaxlNeeded)=1._dl !initialize for CP%flat case 
        else
         EV%k2_buf=EV%q2-CP%curv
         EV%k_buf=sqrt(EV%k2_buf)
                    
         do l=1,EV%MaxlNeeded
           EV%Kf(l)=1._dl-CP%curv*(l*(l+2))/EV%k2_buf           
         end do
        end if

        k=EV%k_buf
        k2=EV%k2_buf
   
        if (CP%closed) then
           EV%FirstZerolForBeta = nint(EV%q*CP%r) 
        else 
           EV%FirstZerolForBeta=l0max !a large number
        end if

   
!  k*tau, (k*tau)**2, (k*tau)**3
        x=k*tau
        x2=x*x
        x3=x2*x
        rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates)) 
        grhonu=rhomass+grhornomass
                                                        
        om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))       
        omtau=om*tau
        Rv=grhonu/(grhonu+grhog)
        Rg = 1-Rv
        Rc=CP%omegac/(CP%omegac+CP%omegab)
        Rb=1-Rc
        Rp15=4*Rv+15

        if (CP%Scalar_initial_condition > initial_nummodes) &
          stop 'Invalid initial condition for scalar modes'

        a=tau*adotrad*(1+omtau/4)
        a2=a*a
        
        initv=0

!  Set adiabatic initial conditions

        chi=1  !Get transfer function for chi
        initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
        initv(1,i_clxr)= initv(1,i_clxg)
        initv(1,i_clxb)=0.75_dl*initv(1,i_clxg)
        initv(1,i_clxc)=initv(1,i_clxb)
        initv(1,i_qg)=initv(1,i_clxg)*x/9._dl
        initv(1,i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
        initv(1,i_vb)=0.75_dl*initv(1,i_qg)
        initv(1,i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
        initv(1,i_aj3r)=chi*4/21._dl/Rp15*x3
        initv(1,i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))
      
        if (CP%Scalar_initial_condition/= initial_adiabatic) then
!CDM isocurvature   
       
         initv(2,i_clxg)= Rc*omtau*(-2._dl/3 + omtau/4)
         initv(2,i_clxr)=initv(2,i_clxg)
         initv(2,i_clxb)=initv(2,i_clxg)*0.75_dl
         initv(2,i_clxc)=1+initv(2,i_clxb)
         initv(2,i_qg)=-Rc/9*omtau*x
         initv(2,i_qr)=initv(2,i_qg)
         initv(2,i_vb)=0.75_dl*initv(2,i_qg)
         initv(2,i_pir)=-Rc*omtau*x2/3/(2*Rv+15._dl)
         initv(2,i_eta)= Rc*omtau*(1._dl/3 - omtau/8)*EV%Kf(1)
         initv(2,i_aj3r)=0
!Baryon isocurvature
         if (Rc==0) stop 'Isocurvature initial conditions assume non-zero dark matter'

         initv(3,:) = initv(2,:)*(Rb/Rc)
         initv(3,i_clxc) = initv(3,i_clxb)
         initv(3,i_clxb)= initv(3,i_clxb)+1
      
!neutrino isocurvature density mode
       
         initv(4,i_clxg)=Rv/Rg*(-1 + x2/6)
         initv(4,i_clxr)=1-x2/6
         initv(4,i_clxc)=-omtau*x2/80*Rv*Rb/Rg
         initv(4,i_clxb)= Rv/Rg/8*x2
         iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
         initv(4,i_qg) =iqg
         initv(4,i_qr) = x/3
         initv(4,i_vb)=0.75_dl*iqg
         initv(4,i_pir)=x2/Rp15
         initv(4,i_eta)=EV%Kf(1)*Rv/Rp15/3*x2
     
!neutrino isocurvature velocity mode

         initv(5,i_clxg)=Rv/Rg*x - 2*x*omtau/16*Rb*(2+Rg)/Rg**2
         initv(5,i_clxr)=-x -3*x*omtau*Rb/16/Rg
         initv(5,i_clxc)=-9*omtau*x/64*Rv*Rb/Rg
         initv(5,i_clxb)= 3*Rv/4/Rg*x - 9*omtau*x/64*Rb*(2+Rg)/Rg**2
         iqg = Rv/Rg*(-1 + 3*Rb/4/Rg*omtau+x2/6 +3*omtau**2/16*Rb/Rg**2*(Rg-3*Rb))
         initv(5,i_qg) =iqg
         initv(5,i_qr) = 1 - x2/6*(1+4*EV%Kf(1)/(4*Rv+5))
         initv(5,i_vb)=0.75_dl*iqg
         initv(5,i_pir)=2*x/(4*Rv+5)+omtau*x*6/Rp15/(4*Rv+5)
         initv(5,i_eta)=2*EV%Kf(1)*x*Rv/(4*Rv+5) + omtau*x*3*EV%Kf(1)*Rv/32*(Rb/Rg - 80/Rp15/(4*Rv+5))
         initv(5,i_aj3r) = 3._dl/7*x2/(4*Rv+5)

!quintessence isocurvature mode


         end if

       if (CP%Scalar_initial_condition==initial_vector) then
          InitVec = 0
          do i=1,initial_nummodes
          InitVec = InitVec+ initv(i,:)*CP%InitialConditionVector(i)
          end do
       else
          InitVec = initv(CP%Scalar_initial_condition,:)
          if (CP%Scalar_initial_condition==initial_adiabatic) InitVec = -InitVec
            !So we start with chi=-1 as before
       end if

        y(1)=a
        y(2)= -InitVec(i_eta)*k/2
        !get eta_s*k, where eta_s is synchronous gauge variable

!  CDM
        y(3)=InitVec(i_clxc)
  
!  Baryons
        y(4)=InitVec(i_clxb)
        y(5)=InitVec(i_vb)

!  Photons
        y(6)=InitVec(i_clxg)
        y(7)=InitVec(i_qg) 
              
        y(8:EV%nvar)=0._dl
 
        if (w_lam /= -1 .and. w_Perturb) then
         y(EV%w_ix) = InitVec(i_clxq)
         y(EV%w_ix+1) = InitVec(i_vq)
        end if

!  Neutrinos
        y(7+EV%lmaxg)=InitVec(i_clxr)
        y(8+EV%lmaxg)=InitVec(i_qr)
        y(9+EV%lmaxg)=InitVec(i_pir)

        if (EV%FirstZerolForBeta==3) then
         y(10+EV%lmaxg)=0._dl
        else
         y(10+EV%lmaxg)=InitVec(i_aj3r)
        endif


        EV%TightCoupling=.true.
        EV%no_rad_multpoles =.false.
        EV%MassiveNuApprox=.false.
        if (CP%Num_Nu_massive == 0) return 

        do nu_i = 1, CP%Nu_mass_eigenstates       
       
        if (EV%NuMethod == Nu_approx) then
        
         !Values are just same as massless neutrino ones since highly relativistic
         off_ix = EV%iq0+(nu_i-1)*(EV%lmaxnu+1)
         y(off_ix)=InitVec(i_clxr)
         y(off_ix+1)=InitVec(i_qr)
         y(off_ix+2)=InitVec(i_pir)
         if (EV%FirstZerolForBeta/=3) then
           y(off_ix+3)=InitVec(i_aj3r)
         endif
        
        else
        
          off_ix =  (nu_i-1)*nqmax*(EV%lmaxnu+1) 
          do  i=1,nqmax
           y(EV%iq0+off_ix+i-1)=-0.25_dl*dlfdlq(i)*InitVec(i_clxr)
           y(EV%iq1+off_ix+i-1)=-0.25_dl*dlfdlq(i)*InitVec(i_qr)
           y(EV%iq2+off_ix+i-1)=-0.25_dl*dlfdlq(i)*InitVec(i_pir)
          end do

        end if
        end do       

        end subroutine initial


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialt(EV,yt,tau)
!  Initial conditions for tensors

        implicit none
        real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
        integer l
         type(EvolutionVars) EV
        real(dl) k,k2 ,a, omtau
        real(dl) yt(EV%nvart)
        real(dl) tens0
        
        if (CP%flat) then
         EV%aux_buf=1._dl
         EV%k2_buf=EV%q2
         EV%k_buf=EV%q       
         EV%Kft(1:EV%lmaxt)=1._dl !initialize for flat case
        else
      
        EV%k2_buf=EV%q2-3*CP%curv
        EV%k_buf=sqrt(EV%k2_buf) 
        EV%aux_buf=sqrt(1._dl+3*CP%curv/EV%k2_buf)  
        do l=1,EV%lmaxt
           EV%Kft(l)=1._dl-CP%curv*((l+1)**2-3)/EV%k2_buf
        end do
        endif
        
        k=EV%k_buf
        k2=EV%k2_buf


        if (CP%closed) then
           EV%FirstZerolForBeta = nint(EV%q*CP%r) 
        else 
           EV%FirstZerolForBeta=l0max !a large number
        end if
         
        a=tau*adotrad
        rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates)) 
        omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+rhomass+grhornomass))       
            
        if (DoTensorNeutrinos) then
         bigR = (rhomass+grhornomass)/(rhomass+grhornomass+grhog)
        else
         bigR = 0._dl
        end if

        x=k*tau
 
        yt(1)=a
        tens0 = 1

        yt(2)= tens0 
!commented things are for the compensated mode with magnetic fields; can be neglected
          ! -15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))
    
        elec=-tens0*(1+2*CP%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1
       
        !shear
        yt(3)=-5._dl/2/(bigR+5)*x*elec 
          !+ 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))
        
        yt(4:EV%nvart)=0._dl
     
!  Neutrinos 
        if (DoTensorNeutrinos) then
         pir=-2._dl/3._dl/(bigR+5)*x**2*elec 
           !+ (bigR-1)/bigR*Magnetic*(1-15./14*x**2/(15+4*bigR))
         aj3r=  -2._dl/21._dl/(bigR+5)*x**3*elec 
           !+ 3._dl/7*x*(bigR-1)/bigR*Magnetic
         yt((EV%lmaxt-1)+(EV%lmaxpolt-1)*2+3+1)=pir
         yt((EV%lmaxt-1)+(EV%lmaxpolt-1)*2+3+2)=aj3r
        end if
    
        end subroutine initialt

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialv(EV,yv,tau)
!  Initial conditions for vectors

        implicit none
        real(dl) bigR,Rc,tau,x,pir
        type(EvolutionVars) EV
        real(dl) k,k2 ,a, omtau
        real(dl) yv(EV%nvarv)
        
        if (CP%flat) then
         EV%k2_buf=EV%q2
         EV%k_buf=EV%q       
        else
         stop 'Vectors not supported in non-flat models'
        endif
        
        k=EV%k_buf
        k2=EV%k2_buf

        omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+grhornomass))       

        a=tau*adotrad*(1+omtau/4)
        
        x=k*tau
       
        bigR = (grhornomass)/(grhornomass+grhog)
        Rc=CP%omegac/(CP%omegac+CP%omegab)

        yv(1)=a

        
        yv(2)= vec_sig0*(1- 15._dl/2*omtau/(4*bigR+15)) + 45._dl/14*x*Magnetic*(BigR-1)/(4*BigR+15)
        !qg
        yv(4)= vec_sig0/3* (4*bigR + 5)/(1-BigR)*(1  -0.75_dl*omtau*(Rc-1)/(bigR-1)* &
               (1 - 0.25_dl*omtau*(3*Rc-2-bigR)/(BigR-1))) &
                 -x/2*Magnetic
        yv(3)= 3._dl/4*yv(4)
        
        yv(5:EV%nvarv) = 0

!        if (.false.) then
!         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = vec_sig0/6/bigR*x**2*(1+2*bigR*omtau/(4*bigR+15))
!         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2) = -2/3._dl*vec_sig0/bigR*x*(1 +3*omtau*bigR/(4*bigR+15))
!         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+3) = 1/4._dl*vec_sig0/bigR*(5+4*BigR) 
!         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+4) =1/9.*x*vec_sig0*(5+4*bigR)/bigR
!         yv(4) = 0
!         yv(3)= 3._dl/4*yv(4)
!          return 
!        end if

!  Neutrinos
!q_r
         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = -1._dl/3*vec_sig0*(4*BigR+5)/bigR &
             + x**2*vec_sig0/6/BigR +0.5_dl*x*(1/bigR-1)*Magnetic 
!pi_r
         pir=-2._dl/3._dl*x*vec_sig0/BigR - (1/bigR-1)*Magnetic
         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +1)=pir
         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +2)=3._dl/7*x*Magnetic*(1-1/BigR)
  
        end subroutine initialv


      subroutine outtransf(EV, y, Arr)
 !write out clxc, clxb, clxg, clxn
        use Transfer
        implicit none
        type(EvolutionVars) EV
   
        real(dl) clxc, clxb, clxg, clxr, k,k2
        real(dl) grho,gpres,dgrho,dgq,a
        real Arr(Transfer_max)
        real(dl) y(EV%nvar)

        a    = y(1)
        clxc = y(3)
        clxb = y(4)
        clxg = y(6)
        clxr = y(7+EV%lmaxg)
        k    = EV%k_buf
        k2   = EV%k2_buf
 
        Arr(Transfer_kh) = k/(CP%h0/100._dl)
        Arr(Transfer_cdm) = clxc/k2
        Arr(Transfer_b) = clxb/k2
        Arr(Transfer_g) = clxg/k2
        Arr(Transfer_r) = clxr/k2
  
        dgrho = 0 
        grho =  0

        if (CP%Num_Nu_Massive > 0) then
          call MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq)
           Arr(Transfer_nu) = dgrho/grho/k2
        else
           Arr(Transfer_nu) = 0
        end if

!!!If we want DE perturbations to get \delta\rho/\rho_m    
     !  dgrho=dgrho+y(EV%w_ix)*grhov*a**(-1-3*w_lam)
     !   Arr(Transfer_r) = y(EV%w_ix)/k2

!        dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a 
!        grho =  grho+(grhoc+grhob)/a + grhov*a**(-1-3*w_lam)

  
        dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a 
        grho =  grho+(grhoc+grhob)/a
        
        Arr(Transfer_tot) = dgrho/grho/k2 
             
    
     end subroutine outtransf

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fderivs(EV,n,tau,ay,ayprime)
!  Evaluate the time derivatives of the perturbations, flat case
!  ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
        use ThermoData
        use MassiveNu
        implicit none      
        type(EvolutionVars) EV

        integer n,nu_i,ix_off
        real(dl) ay(n),ayprime(n)
        real(dl) tau,w
        real(dl) k,k2 

!  Internal variables.

        real(dl) opacity
        real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
                      clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
        real(dl) q,aq,v,akv(nqmax0)
        real(dl) G11_t,G30_t, wnu_arr(max_nu)

        real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
        real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
        real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
        real(dl) f
        real(dl) clxq, vq,  E2, dopacity
        integer l,i,ind, off_ix
                    

        k=EV%k_buf
        k2=EV%k2_buf     
  
        a=ay(1)
        a2=a*a

        etak=ay(2)

!  CDM variables
        clxc=ay(3)

!  Baryon variables
        clxb=ay(4)
        vb=ay(5)

!  Compute expansion rate from: grho 8*pi*rho*a**2

        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        if (w_lam==-1._dl) then
         grhov_t=grhov*a2
         else
         grhov_t=grhov*a**(-1-3*w_lam)
        end if


       if (EV%no_rad_multpoles) then
        clxg=2*(grhoc_t*clxc+grhob_t*clxb)/3/k**2
        clxr=clxg
        qg= clxg*k/sqrt((grhoc_t+grhob_t)/3)*(2/3._dl)
        qr=qg
        pig=0
        E2=0
        pir=0
       else
!  Photons
        clxg=ay(6)
        qg=ay(7)
        pig=ay(8)
        E2=ay(EV%polind+2)

!  Massless neutrinos
        clxr=ay(7+EV%lmaxg)
        qr  =ay(8+EV%lmaxg)
        pir =ay(9+EV%lmaxg)
       end if

!  Get sound speed and ionisation fraction.
        if (EV%TightCoupling) then
          call thermo(tau,cs2,opacity,dopacity)
         else
          call thermo(tau,cs2,opacity)        
        end if

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t

!  Photon mass density over baryon mass density
        photbar=grhog_t/grhob_t
        pb43=4._dl/3*photbar

!  8*pi*a*a*SUM[rho_i*clx_i]
        dgrho=grhob_t*clxb+grhoc_t*clxc + grhog_t*clxg+grhor_t*clxr 
       
!  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        dgq=grhob_t*vb+grhog_t*qg+grhor_t*qr

        gpres=0

        if (CP%Num_Nu_Massive > 0) then
           call MassiveNuVars(EV,ay,a,grho,gpres,dgrho,dgq, wnu_arr)
        end if

        
        adotoa=sqrt(grho/3)
        ayprime(1)=adotoa*a
       
        if (w_lam /= -1 .and. w_Perturb) then
           clxq=ay(EV%w_ix) 
           vq=ay(EV%w_ix+1) 
           dgrho=dgrho + clxq*grhov_t
           dgq = dgq + vq*grhov_t*(1+w_lam)
       end if

!  Get sigma (shear) and z from the constraints
! have to get z from eta for numerical stability
        z=(0.5_dl*dgrho/k + etak)/adotoa 
        sigma=z+1.5_dl*dgq/k2
      
        if (w_lam /= -1 .and. w_Perturb) then

           ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxq+3*adotoa*(1+w_lam)*vq/k) &
               -(1+w_lam)*k*vq -(1+w_lam)*k*z

           ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_lam)

        end if

      
 !eta*k equation
        ayprime(2)=dgq/2

!  CDM equation of motion
        clxcdot=-k*z
        ayprime(3)=clxcdot

!  Baryon equation of motion.
        clxbdot=-k*(z+vb)
        ayprime(4)=clxbdot
!  Photon equation of motion
        clxgdot=-k*(4._dl/3._dl*z+qg)

! Small k: potential problem with stability, using full equations earlier is NOT more accurate in general
! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

!  Use explicit equation for vb if appropriate

         if (EV%TightCoupling) then
   
            pig = 32._dl/45/opacity*k*(sigma+vb)
            E2 = pig/4
            EV%pig = pig
   
    !  Use tight-coupling approximation for vb
    !  zeroth order approximation to vbdot + the pig term
            vbdot=(-adotoa*vb+cs2*k*clxb  &
                 +k/4*pb43*(clxg-2*pig))/(1+pb43)

           !  ddota/a
             gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
             adotdota=(adotoa*adotoa-gpres)/2
!            delta = -1/(4*opacity*(1+pb43))*(k*clxg -4*k*cs2*clxb+4*adotoa*vb)
!            delta = (vb-3._dl/4*qg)

    !  First-order approximation to baryon-photon splip

             slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
             +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

! This approx using n_e \propto 1/S^3 is not good enough in some cases
!            slip=2*pb43/(1+pb43)*adotoa*(vb-3._dl/4*qg) &
!               +1._dl/opacity*(-adotdota*vb-k/2*adotoa*clxg  &
!               +k*(cs2*clxbdot-clxgdot/4))/(1+pb43)

         vbdot=vbdot+pb43/(1+pb43)*slip

        else
            vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity*(4._dl/3*vb-qg)
        end if

        ayprime(5)=vbdot
    
     if (EV%no_rad_multpoles) then

        if (CP%Num_Nu_massive/=0) then
          ayprime(6:EV%polind+EV%lmaxgpol)=0._dl
        end if

     else

 !  Photon equations of motion
        ayprime(6)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +k*(clxg-2*pig)/3
        ayprime(7)=qgdot
        
!  Use explicit equations for photon moments if appropriate       
        if (EV%tightcoupling) then

         !  Use tight-coupling approximation where moments are zero for l>1
          ayprime(8:EV%lmaxg+6)=0._dl
          ayprime(EV%polind+2:EV%polind+EV%lmaxgpol)=0._dl

        else

            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)

            pigdot=0.4_dl*k*qg-0.6_dl*k*ay(9)-opacity*(pig - polter) &
                   +8._dl/15._dl*k*sigma
            ayprime(8)=pigdot
            do  l=3,EV%lmaxg-1
               ayprime(l+6)=k*denl(l)*(l*ay(l+5)-(l+1)*ay(l+7))-opacity*ay(l+6)
            end do
          !  Truncate the photon moment expansion
            ayprime(EV%lmaxg+6)=k*ay(EV%lmaxg+5)-(EV%lmaxg+1)/tau*ay(EV%lmaxg+6)  &
                           -opacity*ay(EV%lmaxg+6)

!  Polarization
            !l=2 
            ayprime(EV%polind+2) = -opacity*(ay(EV%polind+2) - polter) - k/3._dl*ay(EV%polind+3)
            !and the rest
            do l=3,EV%lmaxgpol-1
               ayprime(EV%polind+l)=-opacity*ay(EV%polind+l) + k*denl(l)*(l*ay(EV%polind+l-1) -&
                              polfac(l)*ay(EV%polind+l+1))
            end do  
    
            !truncate
            ayprime(EV%polind+EV%lmaxgpol)=-opacity*ay(EV%polind+EV%lmaxgpol) + &
               k*EV%poltruncfac*ay(EV%polind+EV%lmaxgpol-1)-(EV%lmaxgpol+3)*ay(EV%polind+EV%lmaxgpol)/tau
      
        end if
    
!  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%lmaxg+7)=clxrdot
        qrdot=k*(clxr-2._dl*pir)/3._dl
        ayprime(EV%lmaxg+8)=qrdot
        pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
        ayprime(EV%lmaxg+9)=pirdot
        do l=3,EV%lmaxnr-1
           ayprime(l+EV%lmaxg+7)=k*denl(l)*(l*ay(l+EV%lmaxg+6) -(l+1)*ay(l+EV%lmaxg+8))
        end do   
!  Truncate the neutrino expansion
        ayprime(EV%lmaxnr+EV%lmaxg+7)=k*ay(EV%lmaxnr+EV%lmaxg+6)- &
                                (EV%lmaxnr+1)/tau*ay(EV%lmaxnr+EV%lmaxg+7)
     
        end if ! no_rad_multpoles

!  Massive neutrino equations of motion.
         if (CP%Num_Nu_massive == 0) return
          
          do nu_i = 1, CP%Nu_mass_eigenstates

         if (EV%NuMethod == Nu_approx) then

            !crude approximation scheme
            off_ix = EV%iq0+(nu_i-1)*(EV%lmaxnu+1)

            w=wnu_arr(nu_i)
            f = (5._dl/3)**(a*nu_masses(nu_i)/(a*nu_masses(nu_i)+200)) &
                *(3*w)**((2+a*nu_masses(nu_i))/(4+a*nu_masses(nu_i))) 

            !clxnudot
            ayprime(off_ix)=-(f-3*w)*adotoa*ay(off_ix)-k*((1+w)*z+ay(off_ix+1))
            !qnudot
            ayprime(off_ix+1)=-adotoa*(1-3*w)*ay(off_ix+1)+k/3._dl*(f*ay(off_ix)-2._dl*ay(off_ix+2))
            !pinudot

            ayprime(off_ix+2)=-adotoa*(-f +2-3*w)*ay(off_ix+2) +  &
                 k*(2._dl/5*f*ay(off_ix+1)-0.6_dl*ay(off_ix+3)+ 2._dl/5*w*(5-3*w)*sigma)
 
            do l=3,EV%lmaxnu-1
               ayprime(off_ix+l)=-adotoa*((1-l)*f+l-3*w)*ay(off_ix+l) + &
                    k*denl(l)*(f*l*ay(off_ix+l-1) -(l+1)*ay(off_ix+l+1))
            end do
            !  Truncate the neutrino expansion
            ayprime(off_ix+EV%lmaxnu)=k*ay(off_ix+EV%lmaxnu-1)- &
                 (EV%lmaxnu+1)/tau*ay(off_ix+EV%lmaxnu)
                        
                     
         else !Not approx scheme


          if (EV%MassiveNuApprox) then
             !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
             !see astro-ph/0203507
             G11_t=EV%G11(nu_i)/a/a2 
             G30_t=EV%G30(nu_i)/a/a2
             off_ix = (nu_i-1)*4 + EV%iq0
             w=wnu_arr(nu_i)
             ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
             ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
             ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*ay(off_ix+3)-ay(off_ix+1))
             ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*G30_t-2*G11_t)
     
          else
          
              do i=1,nqmax
                q=(i-0.5_dl)*dq
                aq=a*nu_masses(nu_i)/q
                v=1._dl/sqrt(1._dl+aq*aq)
                akv(i)=k*v        
              end do
        
              ix_off=nqmax*(EV%lmaxnu+1)*(nu_i-1)
    !  l = 0, 1, 2,EV%lmaxnu.
              do i=1,nqmax
               ind=EV%iq0+i-1 + ix_off
               ayprime(ind)=-akv(i)*ay(ind+nqmax)+z*k*dlfdlq(i)/3
               ind=EV%iq1+i-1 + ix_off
               ayprime(ind)=akv(i)*(ay(ind-nqmax)-2*ay(ind+nqmax))/3
               ind=EV%iq2+i-1 + ix_off
               ayprime(ind)=akv(i)*(2*ay(ind-nqmax)-3*ay(ind+nqmax))/5 &
                     -k*2._dl/15._dl*sigma*dlfdlq(i)
               ind=EV%iq0+i-1+EV%lmaxnu*nqmax + ix_off
    !  Truncate moment expansion.
               ayprime(ind)=akv(i)*ay(ind-nqmax)-(EV%lmaxnu+1)/tau*ay(ind)
              end do
              do l=3,EV%lmaxnu-1
                ind=EV%iq0-1+l*nqmax + ix_off
                do i=1,nqmax
                 ind=ind+1
                 ayprime(ind)=akv(i)*denl(l)*(l*ay(ind-nqmax)-(l+1)*ay(ind+nqmax))
                end do
              end do

          end if
          end if

          end do

        end subroutine fderivs

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fderivst(EV,n,tau,ayt,aytprime)
!  Evaluate the time derivatives of the tensor perturbations, CP%flat case
        use ThermoData
        use MassiveNu
        implicit none      
        type(EvolutionVars) EV
        integer n,l,i,ind, nu_i, off_ix
        real(dl), target ::  ayt(n),aytprime(n)
        real(dl) ep,tau,grho,rhopi,cs2,opacity,pirdt
        logical finished_tightcoupling
        real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
        real(dl) q,aq,v,akv(nqmax0)
        real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
        real(dl) Hchi,shear, pig
        real(dl) k,k2,a,a2
        real(dl) pir,adotoa, rhonu, shearnu

         k2=EV%k2_buf
         k=EV%k_buf       

        !E and B start at l=2. Set up pointers accordingly to fill in y arrays
        E => ayt(EV%lmaxt+3-1:)
        Eprime=> aytprime(EV%lmaxt+3-1:) 
        B => E(EV%lmaxpolt:)
        Bprime => Eprime(EV%lmaxpolt:)
     
        a=ayt(1)        

!Hchi = metric perturbation variable, conserved on large scales
        Hchi=ayt(2)


!shear (Hchi' = - k*shear)
        shear=ayt(3)

!  Photon anisotropic stress
        pig=ayt(4)

        a2=a*a

!  Get sound speed and opacity, and see if should use tight-coupling
             
        call thermo(tau,cs2,opacity)
        if (k > 0.06_dl*epsw) then
           ep=ep0
        else
           ep=0.2_dl*ep0
        end if
   
        finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4)) 


! Compute expansion rate from: grho=8*pi*rho*a**2
! Also calculate gpres: 8*pi*p*a**2
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*a**(-1-3*w_lam)

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
     
!Do massive neutrinos
        if (CP%Num_Nu_Massive >0) then
         do nu_i=1,CP%Nu_mass_eigenstates
           call Nu_rho(a*nu_masses(nu_i),rhonu)
           grho=grho+grhormass(nu_i)*rhonu/a2
         end do
        end if


        adotoa=sqrt(grho/3._dl)
       
        aytprime(1)=adotoa*a
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        if (finished_tightcoupling) then
!  Use explicit equations:
!  Equation for the photon anisotropic stress
        aytprime(4)=k*(-1._dl/3._dl*ayt(5)+8._dl/15._dl*shear)  &
                  -opacity*(pig - polter)
! And for the moments            
        do  l=3,EV%lmaxt-1
           aytprime(l+2)=k*denl(l)*(l*ayt(l+1)-   &
                  tensfac(l)*ayt(l+3))-opacity*ayt(l+2)
        end do
!  Truncate the hierarchy
        
        aytprime(EV%lmaxt+2)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(EV%lmaxt+1)- &
                       (EV%lmaxt+3._dl)*ayt(EV%lmaxt+2)/tau-opacity*ayt(EV%lmaxt+2)
     
!E equations

        Eprime(2) = - opacity*(E(2) - polter) + k*(4._dl/6._dl*B(2) - &
                        5._dl/27._dl*E(3))
        do l=3,EV%lmaxpolt-1
        Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
                        tensfacpol(l)*E(l+1)) + 4._dl/(l*(l+1))*B(l))
        end do
!truncate
        Eprime(EV%lmaxpolt)=0._dl
        
        
!B-bar equations
        
        do l=2,EV%lmaxpolt-1
        Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
                        tensfacpol(l)*B(l+1)) - 4._dl/(l*(l+1))*E(l))
        end do
!truncate
        Bprime(EV%lmaxpolt)=0._dl

       else
        pig = 32._dl/45._dl*k/opacity*shear 
        ayt(4) = pig   
        E(2)=pig/4
!  Set the derivatives to zero
        aytprime(4:n)=0._dl
        
        endif
     

        rhopi=grhog_t*pig 

!  Neutrino equations: 
!  Anisotropic stress  
        if (DoTensorNeutrinos) then
        
        neutprime => Bprime(EV%lmaxpolt:)
        neut => B(EV%lmaxpolt:)
               
!  Massless neutrino anisotropic stress
        pir=neut(2)
        
        rhopi=rhopi+grhor_t*pir

        pirdt=-1._dl/3._dl*k*neut(3)+ 8._dl/15._dl*k*shear
        neutprime(2)=pirdt
!  And for the moments
        do  l=3,EV%lmaxnrt-1
           neutprime(l)=k*denl(l)*(l*neut(l-1)- tensfac(l)*neut(l+1))
        end do
        !end if
!  Truncate the hierarchy
        neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
                       (EV%lmaxnrt+3._dl)*neut(EV%lmaxnrt)/tau
    

        !  Massive neutrino equations of motion.
         if (CP%Num_Nu_massive > 0) then
          
          do nu_i=1,CP%Nu_mass_eigenstates

              off_ix = (nu_i-1)*nqmax*(EV%lmaxnut-1)
              do i=1,nqmax
                q=(i-0.5_dl)*dq
                aq=a*nu_masses(nu_i)/q
                v=1._dl/sqrt(1._dl+aq*aq)
                akv(i)=k*v
              end do
              do i=1,nqmax
               ind=EV%iqt+i-1+off_ix                          
               aytprime(ind)=-1._dl/3._dl*akv(i)*ayt(ind+nqmax)- 2._dl/15._dl*k*shear*dlfdlq(i)
               ind=EV%iqt+i-1+(EV%lmaxnut-2)*nqmax +off_ix
    !  Truncate moment expansion.
               aytprime(ind)=akv(i)*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-nqmax)-(EV%lmaxnut+3)/tau*ayt(ind)
              end do
              do l=3,EV%lmaxnut-1
               ind=EV%iqt-1+(l-2)*nqmax +off_ix
               do i=1,nqmax
                ind=ind+1
                aytprime(ind)=akv(i)*denl(l)*(l*ayt(ind-nqmax)-tensfac(l)*ayt(ind+nqmax))     
               end do
              end do

              call Nu_Shear(a*nu_masses(nu_i),shearnu,ayt(EV%iqt+off_ix))                 
              rhopi=rhopi+ grhormass(nu_i)/a2*1.5_dl*shearnu  
 
          end do
                      
         end if
       end if
      
!  Get the propagation equation for the shear
                         
        aytprime(3)=-2*adotoa*shear+k*Hchi-rhopi/k   

        if (finished_tightcoupling) then
!  Use the full expression for pigdt
           EV%tenspigdot=aytprime(4)
        else
!  Use the tight-coupling approximation
           EV%tenspigdot= 32._dl/45._dl*k/opacity*(2._dl*adotoa*shear+aytprime(3))
         endif

        aytprime(2)=-k*shear


end subroutine fderivst
      

        subroutine fderivsv(EV,n,tau,yv,yvprime)
!  Evaluate the time derivatives of the vector perturbations, flat case
        use ThermoData
        use MassiveNu
        implicit none      
        type(EvolutionVars) EV
        integer n,l
        real(dl), target ::  yv(n),yvprime(n)
        real(dl) ep,tau,grho,rhopi,cs2,opacity,gpres
        logical finished_tightcoupling
        real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
        real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
        real(dl) sigma, qg,pig, qr, vb, rhoq, vbdot, photbar, pb43
        real(dl) k,k2,a,a2, adotdota
        real(dl) pir,adotoa  
   
         k2=EV%k2_buf
         k=EV%k_buf       

        !E and B start at l=2. Set up pointers accordingly to fill in y arrays
        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:) 
        B => E(EV%lmaxpolv:)
        Bprime => Eprime(EV%lmaxpolv:)
        neutprime => Bprime(EV%lmaxpolv+1:)
        neut => B(EV%lmaxpolv+1:)
            
        a=yv(1)        

        sigma=yv(2)

        a2=a*a

!  Get sound speed and opacity, and see if should use tight-coupling
             
        call thermo(tau,cs2,opacity)
        if (k > 0.06_dl*epsw) then
           ep=ep0
        else
           ep=0.2_dl*ep0
        end if
   
        finished_tightcoupling = &
         ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4)) 


! Compute expansion rate from: grho=8*pi*rho*a**2
! Also calculate gpres: 8*pi*p*a**2
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*a**(-1-3*w_lam)

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_lam 

        adotoa=sqrt(grho/3._dl)
        adotdota=(adotoa*adotoa-gpres)/2

        photbar=grhog_t/grhob_t
        pb43=4._dl/3*photbar
       
        yvprime(1)=adotoa*a
     
        vb = yv(3)
        qg = yv(4)         
        qr = neut(1)

!  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        rhoq=grhob_t*vb+grhog_t*qg+grhor_t*qr
     !  sigma = 2*rhoq/k**2
        !for non-large k this expression for sigma is unstable at early times
        !so propagate sigma equation separately (near total cancellation in rhoq)
        ! print *,yv(2),2*rhoq/k**2

        if (finished_tightcoupling) then
!  Use explicit equations:

        pig = yv(5) 

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        vbdot = -adotoa*vb-photbar*opacity*(4._dl/3*vb-qg) - 0.5_dl*k*photbar*Magnetic

!  Equation for the photon heat flux stress
  
         yvprime(4)=-0.5_dl*k*pig + opacity*(4._dl/3*vb-qg) 
      
!  Equation for the photon anisotropic stress
        yvprime(5)=k*(2._dl/5*qg -8/15._dl*yv(6))+8._dl/15._dl*k*sigma  &
                  -opacity*(pig - polter)
! And for the moments            
        do  l=3,EV%lmaxv-1
           yvprime(l+3)=k*denl(l)*l*(yv(l+2)-   &
                  vecfac(l)*yv(l+4))-opacity*yv(l+3)
        end do
!  Truncate the hierarchy
        yvprime(EV%lmaxv+3)=k*EV%lmaxv/(EV%lmaxv-1._dl)*yv(EV%lmaxv+2)- &
                       (EV%lmaxv+2._dl)*yv(EV%lmaxv+3)/tau-opacity*yv(EV%lmaxv+3)
     
!E equations
   
        Eprime(2) = - opacity*(E(2) - polter) + k*(1/3._dl*B(2) - &
                        8._dl/27._dl*E(3))
        do l=3,EV%lmaxpolv-1
        Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
                        vecfacpol(l)*E(l+1)) + 2._dl/(l*(l+1))*B(l))
        end do
!truncate
        Eprime(EV%lmaxpolv)=0._dl
                
!B-bar equations
        
        do l=2,EV%lmaxpolv-1
        Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
                        vecfacpol(l)*B(l+1)) - 2._dl/(l*(l+1))*E(l))
        end do
!truncate
        Bprime(EV%lmaxpolv)=0._dl

       else

!Tight coupling expansion results

        pig = 32._dl/45._dl*k/opacity*(vb + sigma)

        EV%pig = pig

        vbdot=(-adotoa*vb  -3._dl/8*pb43*k*Magnetic  -3._dl/8*k*pb43*pig &
                - pb43/(1+pb43)/opacity*(0.75_dl*k*adotoa*pb43**2/(pb43+1)*Magnetic + vb*&
              ( 2*pb43*adotoa**2/(1+pb43) + adotdota)) &
                  )/(1+pb43) 

!  Equation for the photon heat flux
! Get drag from vbdot expression 
        yvprime(4)=-0.5_dl*k*pig - &
           (vbdot+adotoa*vb)/photbar - 0.5_dl*k*Magnetic

!  Set the derivatives to zero
        yvprime(5:n)=0._dl
        yv(5)=pig
        E(2)=  pig/4 
        
        endif
 
        yvprime(3) = vbdot

!  Neutrino equations: 
               
!  Massless neutrino anisotropic stress
        pir=neut(2)
        neutprime(1)= -0.5_dl*k*pir
        neutprime(2)=2._dl/5*k*qr -8._dl/15._dl*k*neut(3)+ 8._dl/15._dl*k*sigma
!  And for the moments
        do  l=3,EV%lmaxnrv-1
           neutprime(l)=k*denl(l)*l*(neut(l-1)- vecfac(l)*neut(l+1))
        end do
        
!  Truncate the hierarchy
        neutprime(EV%lmaxnrv)=k*EV%lmaxnrv/(EV%lmaxnrv-1._dl)*neut(EV%lmaxnrv-1)-  &
                       (EV%lmaxnrv+2._dl)*neut(EV%lmaxnrv)/tau

    
!  Get the propagation equation for the shear
            
        rhopi=grhog_t*pig+grhor_t*pir+ grhog_t*Magnetic
                      
        yvprime(2)=-2*adotoa*sigma -rhopi/k

       end subroutine fderivsv


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine derivs(EV,n,tau,ay,ayprime)
!  Evaluate the time derivatives of the perturbations.
!ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
        use ThermoData
        use MassiveNu
        implicit none      
        type(EvolutionVars) EV

        integer n, nu_i, ix_off
        real(dl) ay(n),ayprime(n)
        real(dl) tau,w
        real(dl) k,k2 

!  Internal variables.

        real(dl) opacity
        real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
                      clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
        real(dl) q,aq,v,akv(nqmax0)
        real(dl) G11_t,G30_t, wnu_arr(max_nu)

        real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,sigma,polter
        real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
        real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
        real(dl) f
        real(dl) clxq, vq,  E2, dopacity
        integer l,i,ind, off_ix
  
        real(dl) cothxor

   
        k2=EV%k2_buf
        k=EV%k_buf

        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r      
      
        a=ay(1)
        a2=a*a

        etak=ay(2)

!  CDM variables
        clxc=ay(3)

!  Baryon variables
        clxb=ay(4)
        vb=ay(5)

!  Compute expansion rate from: grho 8*pi*rho*a**2

        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*a**(-1-3*w_lam)
   
        if (EV%no_rad_multpoles) then
         clxg=2*(grhoc_t*clxc+grhob_t*clxb)/3/k**2
         clxr=clxg
         qg= clxg*k/sqrt((grhoc_t+grhob_t)/3)*(2/3._dl)
         qr=qg
         pig=0
         E2=0
         pir=0
       else
!  Photons
        clxg=ay(6)
        qg=ay(7)
        pig=ay(8)
        E2=ay(EV%polind+2)
!  Massless neutrinos
        clxr=ay(7+EV%lmaxg)
        qr=ay(8+EV%lmaxg)
        pir=ay(9+EV%lmaxg)
        end if

!  Get sound speed and ionisation fraction.
        if (EV%TightCoupling) then
          call thermo(tau,cs2,opacity,dopacity)
         else
          call thermo(tau,cs2,opacity)        
        end if
        
        
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t

!  Photon mass density over baryon mass density
        photbar=grhog_t/grhob_t
        pb43=4._dl/3._dl*photbar        
     
!  8*pi*a*a*SUM[rho_i*clx_i]
        dgrho=grhob_t*clxb+grhoc_t*clxc + grhog_t*clxg+grhor_t*clxr

!  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        dgq=grhob_t*vb+grhog_t*qg+grhor_t*qr

        gpres=0
        
        if (CP%Num_Nu_Massive > 0) then
           call MassiveNuVars(EV,ay,a,grho,gpres,dgrho,dgq, wnu_arr)
        end if

       
        adotoa=sqrt((grho+grhok)/3._dl)
        ayprime(1)=adotoa*a


        if (EV%FirstZerolForBeta <= EV%MaxlNeeded) ayprime(8:EV%ScalEqsToPropagate)=0 !bit lazy this

        if (w_lam /= -1 .and. w_Perturb) then
           clxq=ay(EV%w_ix)
           vq=ay(EV%w_ix+1) 
           dgrho=dgrho + clxq*grhov_t
           dgq = dgq + vq*grhov_t*(1+w_lam)
        end if

!  Get sigma (shear) and z from the constraints
!  have to get z from eta for numerical stability       
        z=(0.5_dl*dgrho/k + etak)/adotoa 
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
        
        if (w_lam /= -1.and. w_Perturb ) then
           ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxq+3*adotoa*(1+w_lam)*vq/k) &
               -(1+w_lam)*k*vq -(1+w_lam)*k*z
           ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_lam)
        end if
 
!  eta*k equation
        ayprime(2)=0.5_dl*dgq + CP%curv*z

!  CDM equation of motion
        clxcdot=-k*z
        ayprime(3)=clxcdot
!  Baryon equation of motion.
        clxbdot=-k*(z+vb)
        ayprime(4)=clxbdot
!  Photon equation of motion
        clxgdot=-k*(4._dl/3._dl*z+qg)

!  Use explicit equation for vb if appropriate
        
        if (EV%TightCoupling) then
  
            pig = 32._dl/45/opacity*k*(sigma+vb)
            E2 = pig/4
            EV%pig = pig

    !  Use tight-coupling approximation for vb
    !  zeroth order (in t_c) approximation to vbdot
            vbdot=(-adotoa*vb+cs2*k*clxb  &
                 +k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1._dl+pb43) 

    !  ddota/a
            gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            adotdota=0.5_dl*(adotoa*adotoa-gpres)

            slip = - (2*adotoa/(1+pb43) + dopacity/opacity)*(vb-3._dl/4*qg) &
                 +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

    !  First-order approximation to vbdot
            vbdot=vbdot+pb43/(1._dl+pb43)*slip

        else

         vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity*(4._dl/3._dl*vb-qg)
          
        end if
        ayprime(5)=vbdot

       if (EV%no_rad_multpoles) then

         if (CP%Num_Nu_massive/=0) then
           ayprime(6:EV%polind+EV%lmaxgpol)=0._dl
         end if

       else

!  Photon equations of motion
        ayprime(6)=clxgdot
        qgdot=4._dl/3._dl*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +k*(clxg-2._dl*pig*EV%Kf(1))/3._dl
        ayprime(7)=qgdot

!  Use explicit equations for photon moments if appropriate
        
        if (EV%tightcoupling) then

    !  Use tight-coupling approximation where moments are zero for l>1
            pigdot=0        
            ayprime(8:EV%lmaxg+6)=0
       
        else

            polter = 0.1_dl*pig+9._dl/15._dl*E2 !2/15*(3/4 pig + 9/2 E2)

            pigdot=0.4_dl*k*qg-0.6_dl*k*EV%Kf(2)*ay(9)-opacity*(pig - polter) &
                   +8._dl/15._dl*k*sigma
            ayprime(8)=pigdot
            do  l=3,min(EV%FirstZerolForBeta,EV%lmaxg)-1
               ayprime(l+6)=k*denl(l)*(l*ay(l+5)-(l+1)*EV%Kf(l)*ay(l+7))-opacity*ay(l+6)
            end do
    !  Truncate the photon moment expansion
            if (EV%lmaxg/=EV%FirstZerolForBeta) then
          
              ayprime(EV%lmaxg+6)=k*ay(EV%lmaxg+5)-(EV%lmaxg+1)*cothxor*ay(EV%lmaxg+6)  &
                           -opacity*ay(EV%lmaxg+6)
            end if
        
        end if

!  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%lmaxg+7)=clxrdot
        qrdot=k*(clxr-2._dl*pir*EV%Kf(1))/3._dl
        ayprime(EV%lmaxg+8)=qrdot
        pirdot=k*(0.4_dl*qr-0.6_dl*EV%Kf(2)*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
        ayprime(EV%lmaxg+9)=pirdot
        do l=3,min(EV%FirstZerolForBeta,EV%lmaxnr)-1
           ayprime(l+EV%lmaxg+7)=k*denl(l)*(l*ay(l+EV%lmaxg+6) -(l+1)*EV%Kf(l)*ay(l+EV%lmaxg+8))
        end do   
!  Truncate the neutrino expansion  
        if (EV%lmaxnr/=EV%FirstZerolForBeta) then
        
         ayprime(EV%lmaxnr+EV%lmaxg+7)=k*ay(EV%lmaxnr+EV%lmaxg+6)- &
                                (EV%lmaxnr+1)*cothxor*ay(EV%lmaxnr+EV%lmaxg+7)
        end if

!  Polarization
        if (EV%TightCoupling) then     

             ayprime(EV%polind+2:EV%polind+EV%lmaxgpol)=0

        else
            !l=2 (defined to be zero for l<2)
             ayprime(EV%polind+2) = -opacity*(ay(EV%polind+2) - polter) - k/3._dl*EV%Kf(2)*ay(EV%polind+3)
            !and the rest
            do l=3,min(EV%FirstZerolForBeta,EV%lmaxgpol)-1
               ayprime(EV%polind+l)=-opacity*ay(EV%polind+l) + k*denl(l)*(l*ay(EV%polind+l-1) -&
                              polfac(l)*EV%Kf(l)*ay(EV%polind+l+1))
            end do  
        
            !truncate
            if (EV%lmaxgpol/=EV%FirstZerolForBeta) then   
         
             ayprime(EV%polind+EV%lmaxgpol)=-opacity*ay(EV%polind+EV%lmaxgpol) +  &
              k*EV%poltruncfac*ay(EV%polind+EV%lmaxgpol-1)-(EV%lmaxgpol+3)*cothxor*ay(EV%polind+EV%lmaxgpol)
            end if
        
        end if

        end if !no_rad_multpoles

!  Massive neutrino equations of motion.
        if (CP%Num_Nu_massive == 0) return

         do nu_i = 1, CP%Nu_mass_eigenstates

         if (EV%NuMethod == Nu_approx) then

            !crude approximation scheme
            off_ix = EV%iq0+(nu_i-1)*(EV%lmaxnu+1)

            w=wnu_arr(nu_i)
            f = (5._dl/3)**(a*nu_masses(nu_i)/(a*nu_masses(nu_i)+200)) &
                *(3*w)**((2+a*nu_masses(nu_i))/(4+a*nu_masses(nu_i))) 

            !clxnudot
            ayprime(off_ix)=-(f-3*w)*adotoa*ay(off_ix)-k*((1+w)*z+ay(off_ix+1))
            !qnudot
            ayprime(off_ix+1)=-adotoa*(1-3*w)*ay(off_ix+1)+k/3._dl*(f*ay(off_ix)-2._dl*EV%Kf(1)*ay(off_ix+2))
            !pinudot

            ayprime(off_ix+2)=-adotoa*(-f +2-3*w)*ay(off_ix+2) +  &
                 k*(2._dl/5*f*ay(off_ix+1)-0.6_dl*EV%Kf(2)*ay(off_ix+3)+ 2._dl/5*w*(5-3*w)*sigma)
 
            do l=3,EV%lmaxnu-1
               ayprime(off_ix+l)=-adotoa*((1-l)*f+l-3*w)*ay(off_ix+l) + &
                    k*denl(l)*(f*l*ay(off_ix+l-1) -EV%Kf(l)*(l+1)*ay(off_ix+l+1))
            end do
            !  Truncate the neutrino expansion
            ayprime(off_ix+EV%lmaxnu)=k*ay(off_ix+EV%lmaxnu-1)- &
                 (EV%lmaxnu+1)*cothxor*ay(off_ix+EV%lmaxnu)
         
         else !Not approx scheme


          if (EV%MassiveNuApprox) then
             !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
             G11_t=EV%G11(nu_i)/a2/a
             G30_t=EV%G30(nu_i)/a2/a
             off_ix = (nu_i-1)*4 + EV%iq0
             w=wnu_arr(nu_i)
             ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
             ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
             ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
             ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
     
          else
          
          
          do i=1,nqmax
            q=(i-0.5_dl)*dq
            aq=a*nu_masses(nu_i)/q
            v=1._dl/sqrt(1._dl+aq*aq)
            akv(i)=k*v
          end do

          ix_off=nqmax*(EV%lmaxnu+1)*(nu_i-1)

!  l = 0, 1, 2,EV%lmaxnu.
          do i=1,nqmax             
           ind=EV%iq0+i-1+ ix_off
           ayprime(ind)=-akv(i)*ay(ind+nqmax)+z*k*dlfdlq(i)/3
           ind=EV%iq1+i-1+ ix_off
           ayprime(ind)=akv(i)*(ay(ind-nqmax)-2*EV%Kf(1)*ay(ind+nqmax))/3
           ind=EV%iq2+i-1+ ix_off
           ayprime(ind)=akv(i)*(2*ay(ind-nqmax)-3*EV%Kf(2)*ay(ind+nqmax))/5 &
                 -k*2._dl/15._dl*sigma*dlfdlq(i)
           ind=EV%iq0+i-1+EV%lmaxnu*nqmax+ ix_off
!  Truncate moment expansion.
           if (EV%lmaxnu<EV%FirstZerolForBeta) then
             ayprime(ind)=akv(i)*ay(ind-nqmax)-(EV%lmaxnu+1)*cothxor*ay(ind)
        
           end if         

          end do
          do l=3,min(EV%FirstZerolForBeta,EV%lmaxnu)-1
           ind=EV%iq0-1+l*nqmax + ix_off
            do i=1,nqmax
             ind=ind+1
             ayprime(ind)=akv(i)*denl(l)*(l*ay(ind-nqmax)-(l+1)*EV%Kf(l)*ay(ind+nqmax))
            end do
          end do  

          end if
          end if
          end do

        end subroutine derivs

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
        subroutine derivst(EV,n,tau,ayt,aytprime)
!  Evaluate the time derivatives of the tensor perturbations.
        use ThermoData
        use MassiveNu
        implicit none      
        type(EvolutionVars) EV
        integer n,l,i,ind, nu_i, off_ix
        real(dl), target ::  ayt(n),aytprime(n)
        real(dl) ep,tau,grho,rhopi,cs2,opacity,pirdt
        logical finished_tightcoupling
        real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
        real(dl) q,aq,v,akv(nqmax0)
        real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
        real(dl) Hchi,shear, pig
        real(dl) k,k2,a,a2
        real(dl) pir, adotoa, rhonu, shearnu

        real(dl) aux, cothxor

      
        k2=EV%k2_buf
        k= EV%k_buf   
        aux=EV%aux_buf 


        !E and B start at l=2. Set up pointers accordingly to fill in ayt arrays
        E => ayt(EV%lmaxt+3-1:)
        Eprime=> aytprime(EV%lmaxt+3-1:) 
        B => E(EV%lmaxpolt:)
        Bprime => Eprime(EV%lmaxpolt:)
      
        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r  
       
        a=ayt(1)        

!  Electric part of Weyl tensor and shear tensor
!        elec=ayt(2)

        Hchi=ayt(2)  

        shear=ayt(3)
!  Photon anisotropic stress
        pig=ayt(4)

        a2=a*a

!  Get sound speed and opacity, and see if should use tight-coupling
             
        call thermo(tau,cs2,opacity)
        if (k > 0.06_dl*epsw) then
           ep=ep0
        else
           ep=0.2_dl*ep0
        end if

        finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep)) 
        

! Compute expansion rate from: grho=8*pi*rho*a**2
! Also calculate gpres: 8*pi*p*a**2
        grhob_t=grhob/a
        grhoc_t=grhoc/a
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        grhov_t=grhov*a**(-1-3*w_lam)

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
     
!Do massive neutrinos
        if (CP%Num_Nu_Massive >0) then
         do nu_i=1,CP%Nu_mass_eigenstates
           call Nu_rho(a*nu_masses(nu_i),rhonu)
           grho=grho+grhormass(nu_i)*rhonu/a2
         end do
        end if

        adotoa=sqrt((grho+grhok)/3._dl)

        aytprime(1)=adotoa*a
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        if (finished_tightcoupling) then
!  Don't use tight coupling approx - use explicit equations:
!  Equation for the photon anisotropic stress
        aytprime(4)=k*(-1._dl/3._dl*EV%Kft(2)*ayt(5)+8._dl/15._dl*shear)  &
                  -opacity*(pig - polter)


        do l=3,min(EV%FirstZerolForBeta,EV%lmaxt)-1
         aytprime(l+2)=k*denl(l)*(l*ayt(l+1)-   &
                  tensfac(l)*EV%Kft(l)*ayt(l+3))-opacity*ayt(l+2)
 
        end do

        !Truncate the hierarchy 
        if (EV%lmaxt/=EV%FirstZerolForBeta) then
         aytprime(EV%lmaxt+2)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(EV%lmaxt+1)- &
                       (EV%lmaxt+3._dl)*cothxor*ayt(EV%lmaxt+2)-opacity*ayt(EV%lmaxt+2)
        end if

!E and B-bar equations
     
        Eprime(2) = - opacity*(E(2) - polter) + k*(4._dl/6._dl*aux*B(2) - &
                        5._dl/27._dl*EV%Kft(2)*E(3))
        
        do l=3,min(EV%FirstZerolForBeta,EV%lmaxpolt)-1
           Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
                        tensfacpol(l)*EV%Kft(l)*E(l+1)) + 4._dl/(l*(l+1))*aux*B(l))
   
       end do

!truncate: difficult, but zetting to zero seems to work OK
        Eprime(EV%lmaxpolt)=0._dl
                    
        do l=2,min(EV%FirstZerolForBeta,EV%lmaxpolt)-1
        Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
                   tensfacpol(l)*EV%Kft(l)*B(l+1)) - 4._dl/(l*(l+1))*aux*E(l))
        end do

!truncate
        Bprime(EV%lmaxpolt)=0._dl
       
        else  !Tight coupling
        ayt(4)=32._dl/45._dl*k/opacity*shear  
        E(2)=ayt(4)/4._dl
        aytprime(4:n)=0._dl
        
        endif
      
      rhopi=grhog_t*pig 


!  Neutrino equations: 
!  Anisotropic stress
        if (DoTensorNeutrinos) then
        
        neutprime => Bprime(EV%lmaxpolt:)
        neut => B(EV%lmaxpolt:)
        
       
!  Massless neutrino anisotropic stress
        pir=neut(2)

        rhopi=rhopi+grhor_t*pir

        pirdt=k*(-1._dl/3._dl*EV%Kft(2)*neut(3)+ 8._dl/15._dl*shear)
        neutprime(2)=pirdt
!  And for the moments
        do  l=3,min(EV%FirstZerolForBeta,EV%lmaxnrt)-1
           neutprime(l)=k*denl(l)*(l*neut(l-1)- tensfac(l)*EV%Kft(l)*neut(l+1))
        end do
        
!  Truncate the hierarchy
        if (EV%lmaxnrt/=EV%FirstZerolForBeta) then
        neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
                       (EV%lmaxnrt+3._dl)*cothxor*neut(EV%lmaxnrt)
        endif


         !  Massive neutrino equations of motion.
         if (CP%Num_Nu_massive > 0) then
          
          do nu_i=1,CP%Nu_mass_eigenstates

              off_ix = (nu_i-1)*nqmax*(EV%lmaxnut-1)

              do i=1,nqmax
                q=(i-0.5_dl)*dq
                aq=a*nu_masses(nu_i)/q
                v=1._dl/sqrt(1._dl+aq*aq)
                akv(i)=k*v
              end do
              do i=1,nqmax
               ind=EV%iqt+i-1+off_ix                          
               aytprime(ind)=-1._dl/3._dl*akv(i)*EV%Kft(2)*ayt(ind+nqmax)- 2._dl/15._dl*k*shear*dlfdlq(i)
               ind=EV%iqt+i-1+(EV%lmaxnut-2)*nqmax+off_ix
    !  Truncate moment expansion.
               if (EV%lmaxnut/=EV%FirstZerolForBeta) then
               aytprime(ind)=akv(i)*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-nqmax)-(EV%lmaxnut+3)*cothxor*ayt(ind)
               end if         
              end do
              do l=3,min(EV%FirstZerolForBeta,EV%lmaxnut)-1
               ind=EV%iqt-1+(l-2)*nqmax +off_ix
                do i=1,nqmax
                ind=ind+1
                aytprime(ind)=akv(i)*denl(l)*(l*ayt(ind-nqmax)-tensfac(l)*EV%Kft(l)*ayt(ind+nqmax))     
               end do
              end do

               call Nu_Shear(a*nu_masses(nu_i),shearnu,ayt(EV%iqt+off_ix))                 
               rhopi=rhopi+ grhormass(nu_i)/a2*1.5_dl*shearnu  
          
          end do
                      
         end if
        end if
             
!  Get the propagation equation for the shear
        
        aytprime(3)=-2*adotoa*shear+k*Hchi*(1+2*CP%curv/k2)-rhopi/k   


!  And the electric part of the Weyl.
        if (finished_tightcoupling) then
!  Use the full expression for pigdt
           EV%tenspigdot=aytprime(4)
        else
!  Use the tight-coupling approximation
           EV%tenspigdot=32._dl/45._dl*k/opacity*(2._dl*adotoa*shear+aytprime(3))
        endif

        aytprime(2)=-k*shear

        end subroutine derivst



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        end module GaugeInterface
