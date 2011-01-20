! Equations module for dark energy with constant equation of state parameter w
! allowing for perturbations based on a quintessence model
! by Antony Lewis (http://cosmologist.info/)

! Feb 2007 changes for 21cm and other power spectra
! July 2007 added perturbed recombination, self-absorption, changes to transfer function output

       module LambdaGeneral
         use precision
         implicit none
          
         real(dl)  :: w_lam = -1 !p/rho for the dark energy (assumed constant) 
         real(dl) :: cs2_lam = 1_dl 
          !comoving sound speed. Always exactly 1 for quintessence 
          !(otherwise assumed constant, though this is almost certainly unrealistic)

         logical :: w_perturb = .true.
          !If you are tempted to set this = .false. read
          ! http://cosmocoffee.info/viewtopic.php?t=811
          ! http://cosmocoffee.info/viewtopic.php?t=512

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
        use RedshiftSpaceData
        use Recombination
        implicit none
        public

   
        !Description of this file. Change if you make modifications.
        character(LEN=*), parameter :: Eqns_name = 'gauge_inv'

        logical, parameter :: plot_evolve = .false. !for outputing time evolution
        real(dl), parameter :: plot_evolve_k = 1_dl
         !note may need to set k_eta_max_scalar  large enough if k for evole is large

        logical :: DoTensorNeutrinos = .false.

        logical :: Evolve_baryon_cs = .true.
         !if true, evolves equation for Delta_{T_m} to get cs_2 = \delta p /\delta\rho for perfect gas
        
        logical :: Evolve_delta_xe =.true. 

        logical :: Evolve_delta_Ts =.false. !Equilibrium result agree to sub-percent level

        logical :: DoLateRadTruncation = .true.
            !if true, use approx to radition perturbations after matter domination on
            !small scales, saving evolution of irrelevant osciallatory multipole equations

        real(dl) :: Magnetic = 0._dl
            !Vector mode anisotropic stress in units of rho_gamma
        real(dl) :: vec_sig0 = 1._dl
            !Vector mode shear      
        integer, parameter :: max_l_evolve = 2024 !Maximum l we are ever likely to propagate

        !Supported scalar initial condition flags
         integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
         initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0
         integer, parameter :: initial_nummodes =  initial_iso_neutrino_vel

        type EvolutionVars
            real(dl) q, q2
            real(dl) k_buf,k2_buf ! set in initial

            integer w_ix !Index of two quintessence equations
            integer line_ix !index of matter temerature perturbation
            integer xe_ix !index of x_e perturbation
            integer Ts_ix !index of Delta_{T_s}

            integer q_ix !index into q_evolve array that gives the value q
            logical TransferOnly

    !       nvar  - number of scalar (tensor) equations for this k       
            integer nvar,nvart, nvarv

           !Max_l for the various hierarchies
            integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
            integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt
            integer lmaxnrv, lmaxv, lmaxpolv
            integer lmaxline !21cm multipoles for getting reionization effect

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

            !Newtonian gauge v_baryon
            real(dl) vb_Newt

            logical :: saha !still high x_e
            
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
         use RECDATA, only : CB1        
         use ThermoData
         type(EvolutionVars) EV
         real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
         integer ind
         real(dl) a, ep, tau_switch, tau_check, Delta_TM
         real(dl) cs2, opacity, xe

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

             if (Evolve_delta_xe .and. tau_switch > recombination_saha_tau) then
                tau_switch = recombination_saha_tau
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

         if (Evolve_delta_xe .and. tau<= recombination_saha_tau .and. tauend > recombination_saha_tau ) then

             if (recombination_saha_tau > tau) call GaugeInterface_ScalEv(EV, y, tau,recombination_saha_tau,tol1,ind,c,w)
             EV%saha = .false.        
             a=y(1)
             Delta_Tm = y(6)/4 ! assume delta_TM = delta_T_gamma
             xe= Recombination_xe(a)
             y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(4) + (3./2+  CB1/(CP%TCMB/a))*Delta_TM)  

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
                if (Evolve_baryon_cs) then
                   y(EV%ScalEqsToPropagate+1) = y(EV%line_ix)
                   EV%line_ix = EV%ScalEqsToPropagate+1
                   EV%ScalEqsToPropagate = EV%ScalEqsToPropagate+1
                end if
                if (Evolve_delta_xe) then
                   y(EV%ScalEqsToPropagate+1) = y(EV%xe_ix)
                   EV%xe_ix = EV%ScalEqsToPropagate+1
                   EV%ScalEqsToPropagate = EV%ScalEqsToPropagate+1
                end if

                if (Evolve_delta_Ts) then
                   y(EV%ScalEqsToPropagate+1) = y(EV%Ts_ix)
                   EV%Ts_ix = EV%ScalEqsToPropagate+1
                   EV%ScalEqsToPropagate = EV%ScalEqsToPropagate+1
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
                  .not. CP%WantTransfer .and. .not. CP%DoLensing .and. num_redshiftwindows==0) then
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
            if (line_phot_quadrupole) then

             EV%lmaxg=EV%lmaxg*8
             EV%lmaxgpol=EV%lmaxgpol*4
         
            elseif (CP%AccurateReionization) then
             EV%lmaxg=EV%lmaxg*4
             EV%lmaxgpol=EV%lmaxgpol*2
            
            end if
         end if       

         if (CP%Transfer%high_precision .or. Do21cm)  EV%lmaxnr=max(nint(25*lAccuracyBoost),3)

         if (Do21cm .and. line_reionization) then
              EV%lmaxg =  EV%lmaxg*8
              EV%lmaxgpol = EV%lmaxgpol*3
         end if
    
         EV%nvar=5+ (EV%lmaxg+1) + EV%lmaxgpol-1 +(EV%lmaxnr+1) 
         if (w_lam /= -1 .and. w_Perturb) then
            EV%w_ix = EV%nvar+1
            EV%nvar=EV%nvar+2
         else
            EV%w_ix=0
         end if

         if (Do21cm .or.Evolve_delta_xe .or. Evolve_delta_Ts) Evolve_baryon_cs = .true.


         if (Evolve_baryon_cs) then
           EV%line_ix = EV%nvar + 1 ! \Delta T_matter
           EV%nvar=EV%nvar+1
         end if

         if (Evolve_delta_xe) then
           EV%xe_ix = EV%nvar + 1 ! \Delta T_matter
           EV%nvar=EV%nvar+1
         end if

         if (Evolve_delta_Ts) then
           EV%Ts_ix = EV%nvar + 1 ! \Delta T_s
           EV%nvar=EV%nvar+1
         end if

         if (Do21cm .and. line_reionization) then           
            EV%lmaxline  = EV%lmaxg
            EV%nvar=EV%nvar+ EV%lmaxline+1 +  EV%lmaxline-1
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
          
          if (Do21cm .and. line_reionization) then
            EV%lmaxt = EV%lmaxt * 8
          end if  
               
          EV%nvart=(EV%lmaxt-1)+(EV%lmaxpolt-1)*2+3

         if (Do21cm) then
            if (line_reionization) then           
             EV%line_ix = EV%nvart + 1 
             EV%lmaxline  = EV%lmaxt
             EV%nvart=EV%nvart+ EV%lmaxline-1 +  (EV%lmaxline-1)*2
           end if
         end if

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
        real(dl) grhormass_t, rhonu, pnu, shearnu
        real(dl) clxnu, qnu, grhonu_t, gpnu_t

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

      function Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, vterm )
      !Delta_Tspin - Delta_TCMB
          use constants
      !vterm = hdot/clh + k*n/3/clh
          real(dl), intent(in) :: a,Delta_source,Delta_TCMB,Delta_Tm,Tmat,Trad,xe, Delta_xe, vterm
          real(dl) :: Get21cm_source2
          real(dl) Rgamma,Rm
          real(dl) dC10, n_H,C10, C10_HH, C10_eH
          real(dl) kappa_HH,kappa_eH
          real(dl) tau_eps
          real(dl) dtauda, H
          external dtauda
          real(dl) TSpin
          n_H = NNow/a**3
          kappa_HH = kappa_HH_21cm(Tmat, .false.)
          kappa_eH = kappa_eH_21cm(Tmat, .false.)
          C10_HH = n_H*kappa_HH* (1- xe)
          C10_eH = n_H*kappa_eH*xe
          C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
          Rgamma = 1._dl/(C10+A10*Trad/T_21cm)
          Rm = 1._dl/(C10+A10*Tmat/T_21cm)

!          TSpin=TsRecfast(a)
!          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
!                  n_H*kappa_pH_21cm(Tmat, .false.)*xe
!          if (a>0.5) stop
  
          dC10 =(C10*Delta_source + &
           (C10_HH*kappa_HH_21cm(Tmat, .true.)+C10_eH*kappa_eH_21cm(Tmat, .true.))*Delta_Tm &
             + (kappa_eH-kappa_HH)*xe*n_H*Delta_xe   ) 

    

          Get21cm_source2 =  dC10*(Rgamma-Rm) +  C10*(Rm*Delta_tm - Delta_TCMB*Rgamma)  

          TSpin=Recombination_Ts(a)
          H = (1/(a*dtauda(a)))
          tau_eps = a*line21_const*NNow/a**3/H/Tspin/1000

          Get21cm_source2 = Get21cm_source2 + &
             tau_eps/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) ) * &
              (Delta_source -vterm + dC10/C10 + 2*( - Rgamma*dC10 + Delta_TCMB*(C10*Rgamma-1)) &
                + Trad/(Tmat-Trad)*(Delta_tm-Delta_TCMB)   ) 

       end function Get21cm_source2



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function Get21cm_dTs(a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
       !d Delta T_s / d eta dropping small \Delta_xe terms
          use constants
          real(dl), intent(in) :: a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe
          real(dl) :: Get21cm_dTs
          real(dl) n_H,C10, C10_HH, C10_eH, delta_C10
          real(dl) kappa_HH,kappa_eH, TSpin

          n_H = NNow/a**3
          kappa_HH = kappa_HH_21cm(Tmat, .false.)
          kappa_eH = kappa_eH_21cm(Tmat, .false.)
          C10_HH = n_H*kappa_HH* (1- xe)
          C10_eH = n_H*kappa_eH*xe
          C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
          TSpin=Recombination_ts(a)
          delta_C10 = C10*Delta_n + (C10_HH*kappa_HH_21cm(Tmat, .true.)+C10_eH*kappa_eH_21cm(Tmat, .true.))*Delta_Tm

!          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
!                  n_H*kappa_pH_21cm(Tmat, .false.)*xe
  
          Get21cm_dTs =  4*a*( TSpin/TMat*(Delta_Tm-Delta_ts)*C10 + (1-TSpin/TMat)*delta_C10 + &
                    (Trad*Delta_TCMB - Tspin*Delta_Ts)*A10/T_21cm ) * MPC_in_sec

       end function Get21cm_dTs



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine output(EV,y, n,j,tau,sources)
        use constants, only : barssc0
        use ThermoData
        use RedshiftSpaceData
        use Recombination
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

        real(dl) k,k2  ,adotoa, grho, gpres,etak,phi,phidot,dgpi
        real(dl) clxq, vq, diff_rhopi
        real(dl) sources(CTransScal%NumSources)
        real(dl) t4,t92
        real(dl) Tmat,Trad, Tspin, Delta_source, Delta_source2
        real(dl) Delta_TCMB, Delta_tm, Delta_xe
        real(dl) polter_line, chi
        integer w_ix, lineoff,lineoffpol
        Type (TRedWin), pointer :: W
       
        real(dl) cs2, xe,opacity, delta_p
        real(dl) s(0:10), t(0:10), ISW
        real(dl) counts_radial_source, counts_velocity_source, counts_density_source, counts_ISW_source, &
                 counts_redshift_source, counts_timedelay_source, counts_potential_source
        sources = 0

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


 if (.not. plot_evolve) then
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
     
     
!doppler term
   !   sources(1)=  (sigma+vb)/k*dvis(j)+((-2.D0*adotoa*sigma+vbdot)/k-1.D0/k**2*dgpi)*vis(j) &
   !         +1.D0/k/EV%Kf(1)*vis(j)*etak
       if (1/a-1 > 20) then !remove numercal problems with velocity one
   !        sources(1)= 0
       end if 

     ! sources(1)= ISW  

      if (x > 0._dl) then
         !E polarization source
           sources(2)=vis(j)*polter*(15._dl/8._dl)/divfac 
               !factor of four because no 1/16 later
        else
           sources(2)=0
        end if

end if !!!!!
 
      if (CTransScal%NumSources > 2) then


         !phi_lens = Phi - 1/2 kappa (a/k)^2 sum_i rho_i pi_i
         !Neglect pi contributions because not accurate at late time anyway
         phi = -(dgrho +3*dgq*adotoa/k)/(k2*EV%Kf(1)*2)
         
       !CMB lensing sources
       if (tau>taurend .and. CP%tau0-tau > 0.1_dl) then
         
         sources(3) = -2*phi*f_K(tau-tau_maxvis)/(f_K(CP%tau0-tau_maxvis)*f_K(CP%tau0-tau))
        !We include the lensing factor of two here

!!!!!
!!      if (lens21cm) then
!        w_ix=1
!      if (tau>Redshift_W(w_ix)%tau) then
!          sources(3) = -2*phi*f_K(tau-Redshift_W(w_ix)%tau)/(f_K(CP%tau0-Redshift_W(w_ix)%tau)*f_K(CP%tau0-tau))
!         else
!          sources(3) = 0   
!        end if
!       else
!         sources(3) = 0
  !     end if
       else 
         sources(3)=0
       end if

   if (line_reionization) sources(2)=0
       
     if (tau>tau_start_redshiftwindows .or. plot_evolve) then
      !There are line of sight contributions...

     if (Do21cm) then

         Delta_TCMB = clxg/4
         Delta_source = clxb  
         Trad = CP%TCMB/a
  
         xe = Recombination_xe(a)   
         Tmat = Recombination_Tm(a)
     
         if (a < Do21cm_mina) then
          Delta_tm = Delta_TCMB + (1-Tmat/Trad)*4*Delta_TCMB
          y(EV%line_ix) = Delta_tm
         else
          Delta_tm = y(EV%line_ix)
         end if
         if (Evolve_delta_xe) then
          Delta_xe = y(EV%xe_ix)
         else
          Delta_xe = 0
         end if
         delta_source2 = Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, &
                                         k*(z+vb)/adotoa/3)

     end if


if (plot_evolve) then
       Tspin = Recombination_Ts(a)

       call thermo(tau,cs2,opacity)        
      
       delta_p = barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe)*Tmat*(clxb + delta_tm)
       xe=Recombination_xe(a)
       stop 'write code in equations.f90::output'
!       write(*,'(5e15.5)') 1/a-1, clxb + Trad/(Tspin-Trad)*delta_source2, delta_source2,  Trad/(Tspin-Trad), clxb

       write(*,'(9e15.5)') 1/a-1, y(EV%xe_ix), clxb, Recombination_xe(a), clxg, tau, delta_TM, delta_p, clxc
       return
end if


    do w_ix = 1, num_redshiftwindows

      W => Redshift_W(w_ix)

      if (W%kind == window_lensing) then

       sources(3+w_ix) =-2*phi*W%win_lens(j)

      elseif (W%kind == window_counts) then
     
        !assume zero velocity bias and relevant tracer is CDM perturbation  
        !neglect anisotropic stress in some places
       
        !phidot neglecting anisotropic stress so phi=psi
         phidot = ((4.D0/3.D0*k*EV%Kf(1)*sigma+(-2.D0/3.D0*sigma-2.D0/3.D0*etak/adotoa)*k &
              -diff_rhopi/k**2-1.D0/adotoa*dgrho/3.D0+(3.D0*gpres+5.D0*grho)*sigma/k/3.D0 &
              -2.D0/k*adotoa/EV%Kf(1)*etak)) / 2
        
        !Main density source
        if (counts_density) then
         counts_density_source= W%wing(j)*(clxc*W%bias + (W%comoving_density_ev(j) - 3*adotoa)*sigma/k) 
          !Newtonian gauge count density; bias assumed to be on synchronous gauge CDM density
        else
         counts_density_source= 0
        endif
        
              
        if (counts_redshift) then    
      !Main redshift distortion from kV_N/H j'' integrated by parts twice (V_N = sigma in synch gauge)
        counts_redshift_source = ((4.D0*adotoa**2+gpres+grho/3.D0)/k*W%wing2(j)+(- &
         4.D0*W%dwing2(j)*adotoa+W%ddwing2(j))/k)*sigma+(-etak/adotoa*k/3.D0-dgrho/ &
         adotoa/6.D0+(etak/adotoa*k/3.D0+dgrho/adotoa/6.D0+(dgq/2.D0-2.D0*etak*adotoa)/k) &
            /EV%Kf(1))*W%wing2(j)+2.D0*W%dwing2(j)*etak/k/EV%Kf(1)         
         if (k>0.8e-2) then
       !    write(*,'(8E15.5)') 1/a-1, k*sigma/adotoa, W%wing(j)*clxc, source_redshift, ISW, W%wing(j),W%wing2(j),W%wingtau(j)
       !    if (1/a-1 < 0.01) stop     
         end if       
        else
          counts_redshift_source= 0   
        end if 
          
       ! 2v j'/(H\chi) geometric term 
        if (CP%tau0-tau > 0.1_dl .and. counts_radial) then
       
             chi =  CP%tau0-tau
             counts_radial_source= (1-2.5*W%dlog10Ndm)*((-4.D0*W%wing2(j)/chi*adotoa-2.D0*(-W%dwing2(j)*chi-W%wing2(j))/chi**2)/ &
                              k*sigma+2.D0*W%wing2(j)*etak/chi/k/EV%Kf(1))
     
        else
            counts_radial_source = 0
        end if
        
!        if (counts_evolve) then
!       !Just source evolution term if window is actual source distribution    
!         counts_evolve_source =(2.D0*W%dwing2(j)*adotoa-W%ddwing2(j))/k*sigma-W%dwing2(j)*etak/k/EV%Kf(1)
!        else
!          counts_evolve_source = 0
!        end if
       
        if (counts_timedelay) then
         !time delay; WinV is int g/chi
           counts_timedelay_source= 2*(1-2.5*W%dlog10Ndm)*W%WinV(j)*2*phi 
        else 
           counts_timedelay_source = 0
        end if
        
        if (counts_ISW) then
          !WinF is int wingtau
           counts_ISW_source = W%WinF(j)*2*phidot
        else
           counts_ISW_source = 0        
        end if
        
        if (counts_potential) then
           !approx phi = psi
            counts_potential_source = ( phidot/adotoa + phi +(5*W%dlog10Ndm-2)*phi ) * W%wing(j) + phi * W%wingtau(j)
        else
            counts_potential_source = 0 
        end if       

        if (counts_velocity) then
          counts_velocity_source =  (-2.D0*W%wingtau(j)*adotoa+W%dwingtau(j))/k*sigma+W%wingtau(j)*etak/k/EV%Kf(1) &
           - counts_radial_source  !don't double count terms; counts_radial is part of counts_velocity with 1/H/chi          
        else
          counts_velocity_source = 0
        end if
  
        sources(3+w_ix)=  counts_radial_source +  counts_density_source + counts_redshift_source &
            + counts_timedelay_source + counts_potential_source &
            + counts_ISW_source + counts_velocity_source
                
        sources(3+w_ix)=sources(3+w_ix)/W%Fq
 
       if (DoRedshiftLensing) &
          sources(3+W%mag_index+num_redshiftwindows) = phi*W%win_lens(j)*(2-5*W%dlog10Ndm) 
  
      elseif (W%kind == window_21cm) then

       if (line_basic) then
         sources(3+w_ix)= expmmu(j)*(W%wing(j)*Delta_source + W%wing2(j)*Delta_source2 &
               - W%Wingtau(j)*(clxb - (Delta_source2+clxg/4))) 
     !!    sources(3+w_ix)= expmmu(j)*W%wing(j)*phi
       else
         sources(3+w_ix)= 0
       end if

     if (line_distortions ) then

!With baryon velocity, dropping small terms
s(1) =  (sigma/adotoa/3.D0-etak/adotoa**2/3.D0)*W%wing(j)*expmmu(j)*k
s(2) =  -1.D0/adotoa**2*expmmu(j)*W%wing(j)*dgrho/6.D0+((((4.D0*sigma+ &
    vb)*adotoa+(-grho*sigma/2.D0-vb*grho/3.D0)/adotoa+(sigma*grho**2/18.D0+ &
    vb*grho**2/18.D0)/adotoa**3)*W%wing(j)-4.D0*W%dwing(j)*sigma+(W%ddwing(j)*sigma+ &
    W%ddwing(j)*vb)/adotoa+(W%dwing(j)*sigma*grho/3.D0+W%dwing(j)*vb*grho/3.D0)/ &
    adotoa**2-2.D0*W%dwing(j)*vb+((-2.D0*etak+etak*grho/adotoa**2/3.D0)*W%wing(j)+ &
    2.D0*W%dwing(j)*etak/adotoa)/EV%Kf(1))*expmmu(j)+(-4.D0*vis(j)*sigma- &
    2.D0*vis(j)*vb+(dvis(j)*sigma+dvis(j)*vb)/adotoa+(vis(j)*grho*sigma/3.D0+ &
    vis(j)*vb*grho/3.D0)/adotoa**2)*W%wing(j)+2.D0*vis(j)*etak/adotoa*W%wing(j)/ &
    EV%Kf(1)+(2.D0*vis(j)*W%dwing(j)*sigma+2.D0*vis(j)*W%dwing(j)*vb)/adotoa)/k
t(0) =  s(1)+s(2)

sources(3+w_ix)= sources(3+w_ix) + t(0)

     end if


    if (line_extra) then


!All sources except below
   if (line_basic .and. line_distortions) then
  sources(3+w_ix) =  (-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)*expmmu(j)*k+ &
    (W%wing2(j)*Delta_source2+W%wing(j)*Delta_source+1.D0/adotoa*W%winV(j)*dgrho/ &
    3.D0)*expmmu(j)+((-W%dwing(j)*vb+(-(3.D0*gpres+grho)*sigma/3.D0- &
    4.D0*adotoa**2*sigma)*W%winV(j)+4.D0*adotoa*W%dwinV(j)*sigma+(-sigma- &
    vb)*W%ddWinV(j)-vbdot*W%wing(j)-W%dwinV(j)*vbdot+(-2.D0*W%dwinV(j)*etak+ &
    2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*expmmu(j)-2.D0*vis(j)*sigma*W%dwinV(j)+ &
    (4.D0*vis(j)*sigma*adotoa-dvis(j)*sigma)*W%winV(j)-2.D0*vis(j)*W%winV(j)*etak/ &
    EV%Kf(1)-vis(j)*W%dwinV(j)*vb-vis(j)*W%wing(j)*vb)/k+((2.D0*W%dwinV(j)*dgpi+ &
    diff_rhopi*W%winV(j))*expmmu(j)+2.D0*vis(j)*W%winV(j)*dgpi)/k**2

  else

s(1) =  ((-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)+(-sigma/adotoa/3.D0+ &
    etak/adotoa**2/3.D0)*W%wing(j))*expmmu(j)*k+(1.D0/adotoa*W%winV(j)*dgrho/3.D0+ &
    1.D0/adotoa**2*W%wing(j)*dgrho/6.D0)*expmmu(j)
s(2) =  s(1)
s(6) =  ((-vb-sigma)*W%ddWinV(j)+(-4.D0*adotoa**2*sigma-(18.D0*gpres+ &
    6.D0*grho)*sigma/18.D0)*W%winV(j)+((-4.D0*sigma-vb)*adotoa-vbdot+(grho*sigma/ &
    2.D0+vb*grho/3.D0)/adotoa+(-grho**2*sigma/18.D0-vb*grho**2/18.D0)/ &
    adotoa**3)*W%wing(j)+W%dwing(j)*vb+(-W%ddwing(j)*sigma-W%ddwing(j)*vb)/adotoa+ &
    4.D0*W%dwinV(j)*sigma*adotoa+4.D0*W%dwing(j)*sigma+(-W%dwing(j)*grho*sigma/3.D0- &
    W%dwing(j)*vb*grho/3.D0)/adotoa**2-W%dwinV(j)*vbdot+((2.D0*etak-etak*grho/ &
    adotoa**2/3.D0)*W%wing(j)-2.D0*W%dwing(j)*etak/adotoa-2.D0*W%dwinV(j)*etak+ &
    2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*expmmu(j)-vis(j)*W%dwinV(j)*vb+ &
    (4.D0*vis(j)*sigma*adotoa-dvis(j)*sigma)*W%winV(j)
s(5) =  s(6)+(-2.D0*vis(j)*etak/adotoa*W%wing(j)-2.D0*vis(j)*W%winV(j)*etak)/ &
    EV%Kf(1)+(4.D0*vis(j)*sigma+(-vis(j)*grho*sigma/3.D0-vis(j)*vb*grho/3.D0)/ &
    adotoa**2+vis(j)*vb+(-dvis(j)*sigma-dvis(j)*vb)/adotoa)*W%wing(j)+(- &
    2.D0*vis(j)*W%dwing(j)*sigma-2.D0*vis(j)*W%dwing(j)*vb)/adotoa- &
    2.D0*vis(j)*W%dwinV(j)*sigma
s(6) =  1.D0/k
s(4) =  s(5)*s(6)
s(5) =  ((diff_rhopi*W%winV(j)+2.D0*W%dwinV(j)*dgpi)*expmmu(j)+ &
    2.D0*vis(j)*dgpi*W%winV(j))/k**2
s(3) =  s(4)+s(5)
t(0) =  s(2)+s(3)

      sources(3+w_ix) =   sources(3+w_ix) + t(0)
   
  end if

    end if



     if (line_reionization) then
        if (num_redshiftwindows>1) stop 'reionization only for one window at the mo'
        lineoff=EV%line_ix+1
        lineoffpol = lineoff+EV%lmaxline-1
        polter_line = 0.1_dl*y(lineoff+2)+9._dl/15._dl*y(lineoffpol+2)

       if (x > 0._dl) then
        sources(2)=vis(j)*polter_line*(15._dl/2._dl)/divfac 
       else
        sources(2)=0
       end if        

     if (.not. use_mK) sources(2)= sources(2) /W%Fq 


s(1) =  vis(j)*y(lineoff+2)/4.D0+vis(j)*y(lineoff)
s(2) =  s(1)
s(4) =  (-1.D0/EV%Kf(1)*vis(j)*W%winV(j)*etak/10.D0-vis(j)*sigma*W%dwinV(j)/10.D0- &
    9.D0/20.D0*vis(j)*yprime(lineoff+2)-27.D0/100.D0*vis(j)*opac(j)*y(lineoff+1)- &
    9.D0/10.D0*dvis(j)*y(lineoff+3)-3.D0/20.D0*vis(j)*opac(j)*EV%Kf(2)*y(lineoffpol+3)+ &
    vis(j)*W%dwinV(j)*vb+81.D0/200.D0*vis(j)*opac(j)*y(lineoff+3)+3.D0/ &
    5.D0*dvis(j)*y(lineoff+1)+3.D0/10.D0*vis(j)*yprime(lineoff+1)+ &
    (vis(j)*adotoa*sigma/5.D0+(36.D0*vis(j)*opac(j)-80.D0*dvis(j))*sigma/400.D0+ &
    dvis(j)*vb+vis(j)*vbdot)*W%winV(j))/k
s(5) =  (vis(j)*W%winV(j)*dgpi/10.D0+9.D0/20.D0*vis(j)*dopac(j)*y(lineoffpol+2)+ &
    261.D0/400.D0*vis(j)*opac(j)**2.D0*y(lineoff+2)-117.D0/ &
    200.D0*vis(j)*opac(j)**2.D0*y(lineoffpol+2)+3.D0/4.D0*ddvis(j)*y(lineoff+2)- &
    27.D0/20.D0*dvis(j)*opac(j)*y(lineoff+2)+9.D0/ &
    10.D0*dvis(j)*opac(j)*y(lineoffpol+2)-27.D0/40.D0*vis(j)*dopac(j)*y(lineoff+2))/ &
    k**2
s(3) =  s(4)+s(5)
t(0) =  s(2)+s(3)

        sources(3+w_ix)= sources(3+w_ix) + t(0)
           
    
     end if

     if (line_phot_quadrupole) then
     
         s(1) =  (EV%kf(1)*W%wing2(j)*pig/2.D0+(-clxg/4.D0-5.D0/ &
        8.D0*pig)*W%wing2(j))*expmmu(j)
    s(3) =  ((-1.D0/EV%kf(1)*W%wing2(j)*etak+(-sigma+9.D0/8.D0*EV%kf(2)*y(9)-3.D0/ &
        4.D0*qg)*W%dwing2(j)+(-opac(j)*vb+2.D0*adotoa*sigma+9.D0/8.D0*EV%kf(2)*yprime(9)+ &
        3.D0/8.D0*opac(j)*EV%kf(2)*ypol(3)+3.D0/4.D0*opac(j)*qg)*W%wing2(j))*expmmu(j)+(- &
        3.D0/4.D0*vis(j)*qg-vis(j)*sigma+9.D0/8.D0*vis(j)*EV%kf(2)*y(9))*W%wing2(j))/k
    s(4) =  (((27.D0/16.D0*opac(j)*pig-15.D0/8.D0*pigdot-9.D0/ &
        8.D0*opac(j)*ypol(2))*W%dwing2(j)+(27.D0/16.D0*dopac(j)*pig+9.D0/ &
        8.D0*opac(j)**2.D0*ypol(2)-9.D0/8.D0*opac(j)**2.D0*polter+27.D0/ &
        16.D0*opac(j)*pigdot+dgpi-9.D0/8.D0*dopac(j)*ypol(2))*W%wing2(j)-15.D0/ &
        8.D0*W%ddwing2(j)*pig)*expmmu(j)-15.D0/4.D0*vis(j)*W%dwing2(j)*pig+(-(- &
        27.D0*vis(j)*opac(j)+30.D0*dvis(j))*pig/16.D0-9.D0/8.D0*vis(j)*opac(j)*ypol(2)- &
        15.D0/8.D0*vis(j)*pigdot)*W%wing2(j))/k**2
    s(2) =  s(3)+s(4)
    t(0) =  s(1)+s(2)

    sources(3+w_ix)= sources(3+w_ix)+ t(0)
 
     end if


     if (line_phot_dipole) then

     sources(3+w_ix)=sources(3+w_ix) + (EV%kf(1)*W%wing2(j)*pig/2.D0-W%wing2(j)*clxg/4.D0)*expmmu(j)+(((vbdot- &
    opac(j)*vb+3.D0/4.D0*opac(j)*qg)*W%wing2(j)+(vb-3.D0/ &
    4.D0*qg)*W%dwing2(j))*expmmu(j)+(vis(j)*vb-3.D0/4.D0*vis(j)*qg)*W%wing2(j))/k
 
     end if

     if (.not. use_mK) sources(3+w_ix)= sources(3+w_ix) /W%Fq 

    end if
 
    end do

    end if

     end if !num sources > 2

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
        integer lineoff, lineoff_E, lineoff_B


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

        if (Do21cm .and. line_reionization) then

          if (tau < tau_start_redshiftwindows) then

            dt=0._dl
            dte=0._dl
            dtb=0._dl
            return

          else 
                lineoff=EV%line_ix  -2
                lineoff_E = EV%line_ix + (EV%lmaxline-1) - 2
                lineoff_B = EV%line_ix + (EV%lmaxline-1)*2 - 2
              
                E => yt(lineoff_E+2-1:)   
                Eprime=> ytprime(lineoff_E+2-1:) 
                Bprime => Eprime(lineoff_B+2-1:)
                pig = yt(lineoff+2)
                polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
                pigdot = ytprime(lineoff+2)

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
   
          end if

        else

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

        end if
   
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

        EV%saha = .true.

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


       if (Evolve_baryon_cs) then
         y(EV%line_ix) = y(6)/4
       end if

       if (Evolve_delta_xe) then
         y(EV%xe_ix) = 0
       end if


       if (Evolve_delta_Ts) then
         y(EV%Ts_ix) = y(EV%line_ix) 
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
        use constants
        implicit none
        type(EvolutionVars) EV
   
        real(dl) clxc, clxb, clxg, clxr, H, k,k2, Tb
        real(dl) grho,gpres,dgrho,dgq,a, tau_fac
        real Arr(Transfer_max)
        real(dl) y(EV%nvar)
        real(dl) Delta_TM,Tmat, Trad,xe, Tspin, delta_source2, Delta_xe,vb
        real(dl) dtauda
        real(dl) tau_eps,z
        external dtauda

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
!       dgrho=dgrho+y(EV%w_ix)*grhov*a**(-1-3*w_lam)

         dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a 
         grho =  grho+(grhoc+grhob)/a
        
         Arr(Transfer_tot) = dgrho/grho/k2 


        if (do21cm) then

         Tmat = Recombination_Tm(a)
         Tspin = Recombination_Ts(a)
         xe = Recombination_xe(a)
         Trad = CP%TCMB/a
         Delta_TM = y(EV%line_ix)
         vb = y(5)
         H = (1/(a*dtauda(a)))

         tau_eps = a*line21_const*NNow/a**3/H/Tspin/1000

         if (Evolve_delta_xe) then
          Delta_xe = y(EV%xe_ix)
         else
          Delta_xe = 0
         end if

         z=(0.5_dl*dgrho/k + y(2))/H 
         delta_source2 = Get21cm_source2(a,clxb,clxg/4,Delta_Tm,Delta_xe,Tmat,Trad,xe,k*(z+vb)/H/3) 
         tau_fac = tau_eps/(exp(tau_eps)-1)
         Arr(Transfer_monopole) = ( clxb + Trad/(Tspin-Trad)*delta_source2 ) /k2 &
              + (tau_fac-1)*(clxb - (delta_source2 + clxg/4)  ) / k2

         Arr(Transfer_vnewt) = tau_fac*k*(EV%vb_newt)/H /k2
         Arr(Transfer_Tmat) =  delta_TM/k2
         if (use_mK) then
           Tb = (1-exp(-tau_eps))*a*(Tspin-Trad)*1000
           
           Arr(Transfer_monopole) = Arr(Transfer_monopole)*Tb
           Arr(Transfer_vnewt) = Arr(Transfer_vnewt)*Tb
           Arr(Transfer_Tmat) = Arr(Transfer_Tmat)*Tb
         end if

  !      Arr(Transfer_nu) = k*(EV%vb_newt)/(1/(a*dtauda(a))) /k2  !k*v/H
  !       Arr(Transfer_tot) = (clxb + Trad/(Tspin-Trad)*delta_source2 ) /k2

        end if        

    
     end subroutine outtransf

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fderivs(EV,n,tau,ay,ayprime)
!  Evaluate the time derivatives of the perturbations, flat case
!  ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
        use constants, only : barssc0, Compton_CT
        use ThermoData
        use MassiveNu
        use Recombination
        use RECDATA, only : CB1        
         
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
        real(dl) rhonu,pnu,clxnu,qnu,f
        real(dl) clxq, vq,  E2, dopacity
        integer l,i,ind, nu_ix, off_ix
        real(dl) xe,Trad, Delta_Ts, Delta_TM, Tmat, Delta_TCMB
        real(dl) delta_p, wing_t, wing2_t,winv_t
        real(dl) Delta_source2, polter_line
        real(dl) Delta_xe
        integer lineoff,lineoffpol
       
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
           
        EV%vb_Newt = vb+sigma
        if (w_lam /= -1 .and. w_Perturb) then
           ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxq+3*adotoa*(1+w_lam)*vq/k) &
               -(1+w_lam)*k*vq -(1+w_lam)*k*z

           ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_lam)

        end if

   
!  Delta_Newt = clxb - 3*adotoa*sigma/k 


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

    if (Evolve_baryon_cs) then
         if (a > Do21cm_mina) then
           Tmat = Recombination_Tm(a)
         else
           Tmat = CP%TCMB/a
         end if
         Delta_TM = ay(EV%line_ix)
         delta_p = barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*opacity*a2/akthom)*Tmat*(clxb + delta_tm)

    else 
        Delta_TM = clxg/4
        delta_p = cs2*clxb
    end if


    if (Evolve_delta_xe) then
      if (EV%saha) then
       xe=Recombination_xe(a)
       Delta_xe = (1-xe)/(2-xe)*(-clxb + (3._dl/2+  CB1/Tmat)*Delta_TM)
      else
       Delta_xe = ay(EV%xe_ix)
      end if
    else
      Delta_xe = 0
    end if

! Small k: potential problem with stability, using full equations earlier is NOT more accurate in general
! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

!  Use explicit equation for vb if appropriate

         if (EV%TightCoupling) then
   
            pig = 32._dl/45/opacity*k*(sigma+vb)
            E2 = pig/4
            EV%pig = pig
   
    !  Use tight-coupling approximation for vb
    !  zeroth order approximation to vbdot + the pig term
            vbdot=(-adotoa*vb+delta_p*k  &
                 +k/4*pb43*(clxg-2*pig))/(1+pb43)

           !  ddota/a
             gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
             adotdota=(adotoa*adotoa-gpres)/2
!            delta = -1/(4*opacity*(1+pb43))*(k*clxg -4*k*cs2*clxb+4*adotoa*vb)
!            delta = (vb-3._dl/4*qg)

    !  First-order approximation to baryon-photon splip

!!!Note really need to generalise for cs2 term
             slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
             +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))


! This approx using n_e \propto 1/S^3 is not good enough in some cases
!            slip=2*pb43/(1+pb43)*adotoa*(vb-3._dl/4*qg) &
!               +1._dl/opacity*(-adotdota*vb-k/2*adotoa*clxg  &
!               +k*(cs2*clxbdot-clxgdot/4))/(1+pb43)

         vbdot=vbdot+pb43/(1+pb43)*slip

        else
            vbdot=-adotoa*vb+k*delta_p-photbar*opacity*(4._dl/3*vb-qg)
        end if

        ayprime(5)=vbdot
    
     if (EV%no_rad_multpoles) then

        if (CP%Num_Nu_massive/=0) then
          ayprime(6:EV%polind+EV%lmaxgpol)=0._dl
        end if

     else

 !  Photon equations of motion
        ayprime(6)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+k*delta_p)/pb43 &
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


       if (Evolve_baryon_cs) then
         if (a > Do21cm_mina) then

          Delta_TCMB = clxg/4
          xe = Recombination_xe(a)
          Trad = CP%TCMB/a
        
!Matter temperature
 !Recfast_CT = (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc [a_R = radiation constant]
           ayprime(EV%line_ix) = -2*k*(z+vb)/3 - a*  Compton_CT * (Trad**4) * xe / (1._dl+xe+fHe) * &
           (  (1- Trad/Tmat)*(Delta_TCMB*4 + Delta_xe/(1+xe/(1+fHe))) + Trad/Tmat*(Delta_Tm - Delta_TCMB)  ) 
 
           if (Evolve_delta_Ts) then
             ayprime(EV%Ts_ix) =  Get21cm_dTs(a,clxb,ay(EV%Ts_ix),Delta_TCMB,Delta_Tm,Tmat,Trad,xe ) 
           end if
         else
           ayprime(EV%line_ix) = -k*(4._dl/3._dl*z+qg)/4  !Assume follows clxg
           if (Evolve_delta_Ts) then
              ayprime(EV%Ts_ix) = ayprime(EV%line_ix)
           end if
         end if
       end if

       if (Evolve_delta_xe) then
            if (EV%saha) then
              ayprime(EV%xe_ix) = 0
            else
              ayprime(EV%xe_ix) = dDeltaxe_dtau(a, Delta_xe,clxb, Delta_Tm, k*z/3,k*vb)  
            end if
       end if


       if (Do21cm) then

         if (a > Do21cm_mina) then

           if (line_reionization) then

            lineoff = EV%line_ix+1
            lineoffpol = lineoff+EV%lmaxline-1

            if (tau> tau_start_redshiftwindows) then

              !Multipoles of 21cm
        
             polter_line = ay(lineoff+2)/10+9._dl/15*ay(lineoffpol+2) 

              call interp_window(Redshift_W(1),tau,wing_t,wing2_t,winv_t)
       
              delta_source2 = Get21cm_source2(a,clxb,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe,k*(z+vb)/adotoa/3)


           !Drop some small terms since mulipoles only enter into reionzation anyway
              !monopole
              ayprime(lineoff) = -k*ay(lineoff+1) +  wing_t * clxb + wing2_t*delta_source2 + k*z/3*winV_t

 
              !dipole
              ayprime(lineoff+1)= k*denl(1)*(ay(lineoff)-2*ay(lineoff+2)) - opacity*ay(lineoff+1) & 
                        -wing2_t * ( qg/4 - vb/3)   ! vb/3*WinV_t)
       
             !quadrupole
              ayprime(lineoff+2)= k*denl(2)*(2*ay(lineoff+1)-3*ay(lineoff+3)) &
                        +opacity*(polter_line -ay(lineoff+2) ) -   2._dl/15*k*sigma*winV_t &
                           - wing2_t * ay(6+2)/4
     
              do  l=3,EV%lmaxline-1
                ayprime(lineoff+l)=k*denl(l)*(l*ay(lineoff+l-1)-(l+1)*ay(lineoff+l+1))-opacity*ay(lineoff+l) &
                           - wing2_t * ay(6+l)/4
       
              end do
               !truncate
               ayprime(lineoff+EV%lmaxline)=k*ay(lineoff+EV%lmaxline-1)-(EV%lmaxline+1)/tau*ay(lineoff+EV%lmaxline)  &
                           -opacity*ay(lineoff+EV%lmaxline) - wing2_t * ay(6+EV%lmaxline)/4

!  21cm Polarization
            !l=2 
            ayprime(lineoffpol+2) = -opacity*(ay(lineoffpol+2) - polter_line) - k/3._dl*ay(lineoffpol+3)
            !and the rest
            do l=3,EV%lmaxline-1
               ayprime(lineoffpol+l)=-opacity*ay(lineoffpol+l) + k*denl(l)*(l*ay(lineoffpol+l-1) -&
                              polfac(l)*ay(lineoffpol+l+1))
            end do  
    
            !truncate
            ayprime(lineoffpol+EV%lmaxline)=-opacity*ay(lineoffpol+EV%lmaxline) + &
               k*EV%poltruncfac*ay(lineoffpol+EV%lmaxline-1)-(EV%lmaxline+3)*ay(lineoffpol+EV%lmaxline)/tau


            else
                ayprime(lineoff:lineoffpol+EV%lmaxline)=0                          
            end if

           end if

         end if
 
       end if


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
        real(dl) pir,pnu,adotoa, rhonu, shearnu

        real(dl) polter_line, wing_t,wing2_t,winv_t
        integer lineoff, lineoff_E, lineoff_B


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
     


        if (Do21cm .and. line_reionization) then
           if (tau> tau_start_redshiftwindows) then

              call interp_window(Redshift_W(1),tau,wing_t,wing2_t,winv_t)

             lineoff=EV%line_ix-2
             lineoff_E = EV%line_ix+ EV%lmaxline -1 -2
             lineoff_B = EV%line_ix+(EV%lmaxline -1)*2 -2

             polter_line = 0.1_dl*ayt(lineoff+2) + 9._dl/15._dl*ayt(lineoff_E+2)

             aytprime(lineoff+2)=k*(-1._dl/3._dl*ayt(lineoff+3)+8._dl/15._dl*shear*winv_t)  &
                         -opacity*(ayt(lineoff+2) - polter_line) - wing2_t*ayt(2+2)

            ! And for the moments            
                    do  l=3,EV%lmaxline-1
                       aytprime(lineoff+l)=k*denl(l)*(l*ayt(lineoff+1)-   &
                              tensfac(l)*ayt(lineoff+3))-opacity*ayt(lineoff+2) - wing2_t*ayt(l+2)
                    end do
            !  Truncate the hierarchy
              aytprime(lineoff+EV%lmaxline)=k*EV%lmaxline/(EV%lmaxline-2)*ayt(lineoff+EV%lmaxline-1)- &
                                   (EV%lmaxline+3._dl)*ayt(lineoff+EV%lmaxline)/tau-opacity*ayt(lineoff+EV%lmaxline) &
                                    - wing2_t*ayt(2+EV%lmaxline) 
     
            !E equations

                    aytprime(lineoff_E+2) = - opacity*(ayt(lineoff_E+2) - polter_line) + &
                      k*(4._dl/6._dl*ayt(lineoff_B+2) -  5._dl/27._dl*ayt(lineoff_E+3))
                    do l=3,EV%lmaxline-1
                    aytprime(lineoff_E+l) =-opacity*ayt(lineoff_E+l) + k*(denl(l)*(l*ayt(lineoff_E+l-1) - &
                                    tensfacpol(l)*ayt(lineoff_E+l+1)) + 4._dl/(l*(l+1))*ayt(lineoff_B+l))
                    end do
            !truncate
                    aytprime(lineoff_E+EV%lmaxline)=0._dl
        
        
            !B-bar equations
        
                    do l=2,EV%lmaxline-1
                    aytprime(lineoff_B+l) =-opacity*ayt(lineoff_B+l) + k*(denl(l)*(l*ayt(lineoff_B+l-1) - &
                                    tensfacpol(l)*ayt(lineoff_B+l+1)) - 4._dl/(l*(l+1))*ayt(lineoff_E+l))
                    end do
            !truncate
                    aytprime(lineoff_B+EV%lmaxline)=0._dl
   

           else
                    aytprime(EV%line_ix:EV%line_ix + (EV%lmaxline-1)*3-1)=0
           end if
        end if

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
        real(dl) sigma, qg,pig, qr, vb, qnu, rhoq, vbdot, photbar, pb43
        real(dl) k,k2,a,a2, adotdota
        real(dl) pir,adotoa,  pnu
   
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
        real(dl) rhonu,pnu,clxnu,qnu,f
        real(dl) clxq, vq,  E2, dopacity
        integer l,i,ind, nu_ix, off_ix
  
        real(dl) cothxor


        stop 'need to add 21cm cs^2 and delta_m for non-flat'

   
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
        real(dl) pir,pnu,adotoa, rhonu, shearnu
        real(dl) aux, cothxor

      
               stop 'not done non-flat'

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
