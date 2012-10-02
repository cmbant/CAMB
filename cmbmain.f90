!     This this is the main CAMB program module.
!
!     Code for Anisotropies in the Microwave Background
!     by Antony lewis (http://cosmologist.info) and Anthony Challinor
!     See readme.html for documentation. 

!     Note that though the code is internally parallelised, it is not thread-safe
!     so you cannot generate more than one model at the same time in different threads.
!
!     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
!     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
!     Original CMBFAST copyright and disclaimer:
!
!     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
!     the Massachusetts Institute of Technology.  All rights reserved.
!
!     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
!     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPlIED.
!     By way of example, but not limitation,
!     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
!     MERCHANTABIlITY OR FITNESS FOR ANY PARTICUlAR PURPOSE OR THAT
!     THE USE OF THE lICENSED SOFTWARE OR DOCUMENTATION WIll NOT INFRINGE
!     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
!
!     portions of this software are based on the COSMICS package of
!     E. Bertschinger.  See the lICENSE file of the COSMICS distribution
!     for restrictions on the modification and distribution of this software.

  module CAMBmain

!     This code evolves the linearized perturbation equations of general relativity,
!     the Boltzmann equations and the fluid equations for perturbations
!     of a Friedmann-Robertson-Walker universe with a supplied system of gauge-dependent equation
!     in a modules called GaugeInterface.  The sources for the line of sight integral are
!     computed at sampled times during the evolution for various of wavenumbers. The sources
!     are then interpolated to a denser wavenumber sampling for computing the line of
!     sight integrals of the form Integral d(conformal time) source_k * bessel_k_l. 
!     For CP%flat models the bessel functions are interpolated from a pre-computed table, for
!     non-CP%flat models the hyperspherical Bessel functions are computed by integrating their
!     differential equation. Both phases ('Evolution' and 'Integration') can do separate 
!     wavenumbers in parallel.

!     The time variable is conformal  time dtau=dt/a(t) and the spatial dependence is Fourier transformed 
!     with q=sqrt(k**2 + (|m|+1)K), comoving distances are x=CP%r/a(t), with a(t)=1 today.  
!     The units of both length and time are Mpc.

!    Many elements are part of derived types (to make thread safe or to allow non-sequential code use
!    CP = CAMB parameters
!    EV = Time evolution variables
!    IV = Source integration variables
!    CT = Cl transfer data
!    MT = matter transfer data

! Modules are defined in modules.f90, except GaugeInterface which gives the background and gauge-dependent 
! perturbation equations, and InitialPower which provides the initial power spectrum.

      use precision
      use ModelParams
      use ModelData
      use GaugeInterface
      use Transfer
      use SpherBessels
      use lvalues
      use MassiveNu
      use InitialPower
      use Errors

      implicit none

      private

      logical ExactClosedSum  !do all nu values in sum for Cls for Omega_k>0.1
  

      !Variables for integrating the sources with the bessel functions for each wavenumber
       type IntegrationVars
          integer q_ix
          real(dl) q, dq    !q value we are doing and delta q
  !        real(dl), dimension(:,:), pointer :: Delta_l_q
            !Contribution to C_l integral from this k 
          real(dl), dimension(:,:), pointer :: Source_q, ddSource_q
            !Interpolated sources for this k
    
          integer SourceSteps !number of steps up to where source is zero

      end type IntegrationVars
      
      integer SourceNum
      !SourceNum is total number sources (2 or 3 for scalars, 3 for tensors).

      real(dl) tautf(0:max_transfer_redshifts)  !Time of Trasfer%redshifts
        

      real(dl), dimension(:,:,:), allocatable :: Src, ddSrc !Sources and second derivs
        ! indices  Src( k_index, source_index, time_step_index )    
  
      real(dl), dimension(:,:,:), allocatable :: iCl_scalar, iCl_vector,iCl_tensor
       ! Cls at the l values we actually compute,  iCl_xxx(l_index, Cl_type, initial_power_index)
      
      Type(Regions) :: Evolve_q

      real(dl),parameter :: qmin0=0.1_dl

      real(dl) :: dtaurec_q
    
!     qmax - CP%Max_eta_k/CP%tau0, qmin = qmin0/CP%tau0 for flat case

      real(dl) qmin, qmax 

      real(dl) max_etak_tensor , max_etak_vector, max_etak_scalar
 !     Will only be calculated if k*tau < max_etak_xx

      integer maximum_l !Max value of l to compute
      real(dl) :: maximum_qeta = 3000._dl

      real(dl) :: fixq = 0._dl !Debug output of one q
      
      real(dl) :: ALens = 1._dl
      
      Type(ClTransferData), pointer :: ThisCT
                    
      public cmbmain, ALens, ClTransferToCl, InitVars !InitVars for BAO hack

contains  

     
      subroutine cmbmain
      integer q_ix 
      type(EvolutionVars) EV
  
!     Timing variables for testing purposes. Used if DebugMsgs=.true. in ModelParams
      real(sp) actual,timeprev,starttime

      if (CP%WantCls) then
        
         if (CP%WantTensors .and. CP%WantScalars) stop 'CMBMAIN cannot generate tensors and scalars'
         !Use CAMB_GetResults instead

         if (CP%WantTensors) then
            maximum_l = CP%Max_l_tensor 
            maximum_qeta = CP%Max_eta_k_tensor
         else
            maximum_l = CP%Max_l
            maximum_qeta = CP%Max_eta_k
         end if

   
         call initlval(lSamp, maximum_l)
         
         if (CP%flat)  call InitSpherBessels
         !This is only slow if not called before with same (or higher) Max_l, Max_eta_k
         !Preferably stick to Max_l being a multiple of 50
      end if 


      if (DebugMsgs .and. Feedbacklevel > 0) then
         actual=GetTestTime()
         starttime=actual !times don't include reading the Bessel file
       end if
    
      call InitVars !Most of single thread time spent here (in InitRECFAST)
      if (global_error_flag/=0) return

      if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*) actual-timeprev,' Timing for InitVars'
         write (*,*) 'r = ',real(CP%r),' scale = ',real(scale), 'age = ', real(CP%tau0)  
      end if 

       if (.not. CP%OnlyTransfers)  call InitializePowers(CP%InitPower,CP%curv)
       if (global_error_flag/=0) return

!     Calculation of the CMB sources.


      if (CP%WantCls) call SetkValuesForSources 

      if (CP%WantTransfer) call InitTransfer
 
!      ***note that !$ is the prefix for conditional multi-processor compilation***
      !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)
        
   
      if (CP%WantCls) then

         if (DebugMsgs .and. Feedbacklevel > 0) write(*,*) 'Set ',Evolve_q%npoints,' source k values'
         
         call GetSourceMem

         if (CP%WantScalars) then
             ThisCT => CTransScal
         else if (CP%WantVectors) then
             ThisCT => CTransVec
         else
             ThisCT => CTransTens
         end if

         ThisCT%NumSources = SourceNum
         ThisCT%ls = lSamp

         !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
         !$OMP & PRIVATE(EV, q_ix)
         do q_ix= 1,Evolve_q%npoints
             if (global_error_flag==0) call DoSourcek(EV,q_ix)
         end do
         !$OMP END PARAllEl DO
        
         if (DebugMsgs .and. Feedbacklevel > 0) then
            timeprev=actual
            actual=GetTestTime()
            write(*,*) actual-timeprev,' Timing for source calculation'
         end if

      endif !WantCls

   
!     If transfer functions are requested, set remaining k values and output
      if (CP%WantTransfer .and. global_error_flag==0) then
        call TransferOut
         if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*) actual-timeprev,' Timing for transfer k values'
         end if  
      end if

       if (CP%WantTransfer .and. CP%WantCls .and. CP%DoLensing &
            .and. CP%NonLinear==NonLinear_Lens .and. global_error_flag==0) then
          
          call NonLinearLensing
          if (DebugMsgs .and. Feedbacklevel > 0) then
             timeprev=actual
             actual=GetTestTime()
             write(*,*) actual-timeprev,' Timing for NonLinear'
          end if

       end if

       if (CP%WantTransfer .and. .not. CP%OnlyTransfers .and. global_error_flag==0) &
          call Transfer_Get_sigma8(MT,8._dl) 
           !Can call with other arguments if need different size
 
!     if CMB calculations are requested, calculate the Cl by
!     integrating the sources over time and over k.


      if (CP%WantCls) then
         
         if (global_error_flag==0) then
          
         call InitSourceInterpolation   
       
         ExactClosedSum = CP%curv > 5e-9_dl .or. scale < 0.93_dl
     
         call SetkValuesForInt

         if (DebugMsgs .and. Feedbacklevel > 0) write(*,*) 'Set ',ThisCT%q%npoints,' integration k values'
   
      !Begin k-loop and integrate Sources*Bessels over time

      
      !$OMP PARAllEl DO DEFAUlT(SHARED),SHARED(TimeSteps), SCHEDUlE(STATIC,4) 
          do q_ix=1,ThisCT%q%npoints
            call SourceToTransfers(q_ix)
         end do !q loop
       !$OMP END PARAllEl DO 
  
        if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*)actual-timeprev,' Timing For Integration'
        end if
 
        end if
        
        call FreeSourceMem

        !Final calculations for CMB output unless want the Cl transfer functions only.

        if (.not. CP%OnlyTransfers .and. global_error_flag==0) &
          call ClTransferToCl(CTransScal,CTransTens, CTransVec)

      end if

 
      if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual = GetTestTime()
         write(*,*) actual - timeprev,' Timing for final output'
         write(*,*) actual -starttime,' Timing for whole of cmbmain'
      end if

      end subroutine cmbmain

     subroutine ClTransferToCl(CTransS,CTransT, CTransV)
        Type(ClTransferData) :: CTransS,CTransT, CTransV

       if (CP%WantScalars .and. global_error_flag==0) then
           lSamp = CTransS%ls
           allocate(iCl_Scalar(CTransS%ls%l0,C_Temp:C_last,CP%InitPower%nn))
           iCl_scalar = 0
           
           call CalcScalCls(CTransS)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcScalCls'
       end if    

       if (CP%WantVectors .and. global_error_flag==0) then
           allocate(iCl_vector(CTransV%ls%l0,C_Temp:CT_Cross,CP%InitPower%nn))
           iCl_vector = 0
           call CalcVecCls(CTransV,GetInitPowerArrayVec)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcVecCls'
       end if    


       if (CP%WantTensors .and. global_error_flag==0) then
           allocate(iCl_Tensor(CTransT%ls%l0,CT_Temp:CT_Cross,CP%InitPower%nn))
           iCl_tensor = 0
           call CalcTensCls(CTransT,GetInitPowerArrayTens)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcTensCls'
       end if

       if (global_error_flag==0) then
        call Init_Cls
 
  !     Calculating Cls for every l.
        call InterpolateCls(CTransS,CTransT, CTransV)
   
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'InterplolateCls'

       end if

       if (CP%WantScalars .and. allocated(iCl_Scalar)) deallocate(iCl_scalar)
       if (CP%WantVectors .and. allocated(iCl_Vector)) deallocate(iCl_vector)
       if (CP%WantTensors .and. allocated(iCl_Tensor)) deallocate(iCl_tensor)
       
       if (global_error_flag/=0) return
 
       if (CP%OutputNormalization >=2) call NormalizeClsAtl(CP%OutputNormalization)
       !Normalize to C_l=1 at l=OutputNormalization 
     
     end subroutine ClTransferToCl
    

     subroutine SourceToTransfers(q_ix)
      integer q_ix
      type(IntegrationVars) :: IV
    
          allocate(IV%Source_q(TimeSteps%npoints,SourceNum))
          if (.not.CP%flat) allocate(IV%ddSource_q(TimeSteps%npoints,SourceNum))

            call IntegrationVars_init(IV)

            IV%q_ix = q_ix
            IV%q =ThisCT%q%points(q_ix)
            IV%dq= ThisCT%q%dpoints(q_ix)

            call InterpolateSources(IV)
           
            call DoSourceIntegration(IV)

          if (.not.CP%flat) deallocate(IV%ddSource_q) 
          deallocate(IV%Source_q)
         
     end subroutine SourceToTransfers
    

      subroutine InitTransfer
        integer nu,lastnu, ntodo, nq, q_ix, first_i
        real(dl) dlog_lowk1,dlog_lowk, d_osc,dlog_osc, dlog_highk, boost
        real(dl) amin,q_switch_lowk,q_switch_lowk1,q_switch_osc,q_switch_highk
        real(dl), dimension(:), allocatable :: q_transfer

       if (CP%Transfer%k_per_logint==0) then
        !Optimized spacing
        !Large log spacing on superhorizon scales
        !Linear spacing for horizon scales and first few baryon osciallations
        !Log spacing for last few oscillations
        !large log spacing for small scales

           boost = AccuracyBoost 
           if (CP%Transfer%high_precision) boost = boost*1.5

           q_switch_lowk1 = 0.7/taurst
           dlog_lowk1=2*boost

           q_switch_lowk = 8/taurst
           dlog_lowk=8*boost
           if (HighAccuracyDefault) dlog_lowk = dlog_lowk*2.5

           q_switch_osc = min(CP%Transfer%kmax,30/taurst)
           d_osc= 200*boost
           if (HighAccuracyDefault) d_osc = d_osc*1.8
    
           q_switch_highk = min(CP%Transfer%kmax,60/taurst)
           dlog_osc = 17*boost 
           if (HighAccuracyDefault) q_switch_highk = min(CP%Transfer%kmax,90/taurst)

           !Then up to kmax
           dlog_highk = 3*boost
            
           amin = 5e-5_dl

           nq=int((log(CP%Transfer%kmax/amin))*d_osc)+1 
           allocate(q_transfer(nq))
     
           nq=int((log(q_switch_lowk1/amin))*dlog_lowk1)+1 
           do q_ix=1, nq
            q_transfer(q_ix) = amin*exp((q_ix-1)/dlog_lowk1)
           end do
           MT%num_q_trans = nq

           nq=int(log( q_switch_lowk/q_transfer(MT%num_q_trans))*dlog_lowk) +1
           do q_ix=1, nq
            q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_lowk)
           end do
           MT%num_q_trans = MT%num_q_trans + nq

           nq=int((q_switch_osc-q_transfer(MT%num_q_trans))*d_osc)+1 
           do q_ix=1, nq
            q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)+ q_ix/d_osc
           end do
           MT%num_q_trans = MT%num_q_trans + nq

           if (CP%Transfer%kmax > q_transfer(MT%num_q_trans)) then
            nq=int(log( q_switch_highk/q_transfer(MT%num_q_trans))*dlog_osc) +1
            do q_ix=1, nq
             q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_osc)
            end do
            MT%num_q_trans = MT%num_q_trans + nq
           end if

           if (CP%Transfer%kmax > q_transfer(MT%num_q_trans)) then
            nq=int(log(CP%Transfer%kmax/q_transfer(MT%num_q_trans))*dlog_highk)+1 
            do q_ix=1, nq
             q_transfer(MT%num_q_trans+q_ix) = q_transfer(MT%num_q_trans)*exp(q_ix/dlog_highk)
            end do
            MT%num_q_trans = MT%num_q_trans + nq
           end if
   
        else
         !Fixed spacing
          MT%num_q_trans=int((log(CP%Transfer%kmax)-log(qmin))*CP%Transfer%k_per_logint)+1
          allocate(q_transfer(MT%num_q_trans))
          do q_ix=1, MT%num_q_trans
            q_transfer(q_ix) = qmin*exp(real(q_ix)/CP%Transfer%k_per_logint)
          end do
        end if

         if (CP%closed) then
              lastnu=0
              ntodo = 0
               do q_ix=1,MT%num_q_trans
                nu =nint(CP%r*q_transfer(q_ix))
                if (.not. ((nu<3).or.(nu<=lastnu))) then
                   ntodo=ntodo+1
                   q_transfer(ntodo)= nu/CP%r
                   lastnu=nu
                end if
               end do
              MT%num_q_trans = ntodo
         end if

         if (CP%WantCls) then
               ntodo = MT%num_q_trans
               first_i = ntodo+1
               do q_ix = 1,ntodo
                if (q_transfer(q_ix) > qmax) then
                      first_i=q_ix 
                      exit
                end if
               end do            
           
               if (first_i > ntodo) then
               MT%num_q_trans = Evolve_q%npoints 
              else
               MT%num_q_trans = Evolve_q%npoints + (ntodo - first_i+1) 
              end if
              call Transfer_Allocate(MT)

              MT%q_trans(1:Evolve_q%npoints) = Evolve_q%points(1:Evolve_q%npoints)
              if (MT%num_q_trans > Evolve_q%npoints) then
               MT%q_trans(Evolve_q%npoints+1:MT%num_q_trans) = q_transfer(first_i:ntodo)
          end if

         else
             Evolve_q%npoints = 0
             call Transfer_Allocate(MT)
             MT%q_trans = q_transfer(1:MT%num_q_trans)
         end if
          
         deallocate(q_transfer)
  
      end  subroutine InitTransfer

      function GetTauStart(q)
        real(dl), intent(IN) :: q
        real(dl) taustart, GetTauStart

!     Begin when wave is far outside horizon.
!     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
!     of relativistic neutrinos.
            if (CP%flat) then
              taustart=0.001_dl/q
            else
              taustart=0.001_dl/sqrt(q**2-CP%curv)
            end if

            if (fixq/=0._dl) then
              taustart = 0.001_dl/fixq  
            end if
            
!     Make sure to start early in the radiation era.
           taustart=min(taustart,0.1_dl)

!     Start when massive neutrinos are strongly relativistic.
            if (CP%Num_nu_massive>0) then
               taustart=min(taustart,1.d-3/maxval(nu_masses(1:CP%Nu_mass_eigenstates))/adotrad)
            end if

            GetTauStart=taustart
      end function GetTauStart


      subroutine DoSourcek(EV,q_ix)
        integer q_ix
        real(dl) taustart
        type(EvolutionVars) EV

            EV%q=Evolve_q%points(q_ix) 
 
            if (fixq/=0._dl) then
                EV%q=fixq !for testing
            end if
            EV%q2=EV%q**2

            EV%q_ix = q_ix
            EV%TransferOnly=.false.
      
            taustart = GetTauStart(EV%q)
   
            call GetNumEqns(EV)

            if (CP%WantScalars .and. global_error_flag==0) call CalcScalarSources(EV,taustart)
            if (CP%WantVectors .and. global_error_flag==0) call CalcVectorSources(EV,taustart)
            if (CP%WantTensors .and. global_error_flag==0) call CalcTensorSources(EV,taustart)

      end subroutine DoSourcek

       subroutine GetSourceMem
     
        if (CP%WantScalars) then
           if (CP%Dolensing) then
            SourceNum=3
            C_last = C_PhiE
         else
            SourceNum=2
            C_last = C_Cross
           end if
        else
           SourceNum=3 
        end if
       
        allocate(Src(Evolve_q%npoints,SourceNum,TimeSteps%npoints))
        Src=0
        allocate(ddSrc(Evolve_q%npoints,SourceNum,TimeSteps%npoints))
    
       end subroutine GetSourceMem


       subroutine FreeSourceMem
         
        if (allocated(Src))deallocate(Src, ddSrc)
        call Ranges_Free(Evolve_q)

       end subroutine FreeSourceMem


!  initial variables, number of steps, etc.
      subroutine InitVars
      use ThermoData
      use precision
      use ModelParams
      
      implicit none
      real(dl) taumin, maxq, initAccuracyBoost
      integer itf
      
      initAccuracyBoost = AccuracyBoost 

 ! Maximum and minimum k-values.      
      if (CP%flat) then
      qmax=maximum_qeta/CP%tau0
      qmin=qmin0/CP%tau0/initAccuracyBoost 
      else              
        qmax=maximum_qeta/CP%r/CP%chi0
        qmin=qmin0/CP%r/CP%chi0/initAccuracyBoost
      end if
!     Timesteps during recombination (tentative, the actual
!     timestep is the minimum between this value and taurst/40,
!     where taurst is the time when recombination starts - see inithermo

      dtaurec_q=4/qmax/initAccuracyBoost 
      if (.not. CP%flat) dtaurec_q=dtaurec_q/6
      !AL:Changed Dec 2003, dtaurec feeds back into the non-flat integration via the step size
      dtaurec = dtaurec_q
      !dtau rec may be changed by inithermo

      max_etak_tensor = initAccuracyBoost*maximum_qeta /10  
      max_etak_scalar = initAccuracyBoost*max(1700._dl,maximum_qeta) /20 
      if (maximum_qeta <3500 .and. AccuracyBoost < 2) max_etak_scalar = max_etak_scalar * 1.5
        !tweak to get large scales right
      max_etak_vector = max_etak_scalar
      
      if (CP%WantCls) then     
         maxq = qmax
         if (CP%WantTransfer) maxq=max(qmax,CP%Transfer%kmax)
      else
         maxq=CP%Transfer%kmax
      end if


      taumin=GetTauStart(maxq)
   
!     Initialize baryon temperature and ionization fractions vs. time.
!     This subroutine also fixes the timesteps where the sources are
!     saved in order to do the integration. So TimeSteps is set here.
      !These routines in ThermoData (modules.f90)
      call inithermo(taumin,CP%tau0)
      if (global_error_flag/=0) return
   
      if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'inithermo'

!Do any array initialization for propagation equations
      call GaugeInterface_Init
 
      if (Feedbacklevel > 0)  &
           write(*,'("tau_recomb/Mpc       = ",f7.2,"  tau_now/Mpc = ",f8.1)') tau_maxvis,CP%tau0

!     Calculating the times for the outputs of the transfer functions.
!
      if (CP%WantTransfer) then
         do itf=1,CP%Transfer%num_redshifts
            tautf(itf)=min(TimeOfz(CP%Transfer%redshifts(itf)),CP%tau0)
            if (itf>1) then
             if (tautf(itf) <= tautf(itf-1)) then
               stop 'Transfer redshifts not set or out of order'
             end if
            end if
         end do
      endif
 
      end subroutine InitVars

      subroutine SetkValuesForSources
      implicit none
      real(dl) dlnk0, dkn1, dkn2, q_switch
      real(dl) qmax_log
      real(dl) SourceAccuracyBoost
!     set k values for which the sources for the anisotropy and
!     polarization will be calculated. For low values of k we
!     use a logarithmic spacing. closed case dealt with by SetClosedkValues

         SourceAccuracyBoost = AccuracyBoost  
         if (CP%WantScalars .and. CP%Reion%Reionization .and. CP%AccuratePolarization) then
            dlnk0=2._dl/10/SourceAccuracyBoost 
            !Need this to get accurate low l polarization
         else
            dlnk0=5._dl/10/SourceAccuracyBoost
            if (CP%closed) dlnk0=dlnk0/2
         end if

         if (CP%AccurateReionization) dlnk0 = dlnk0/2

         dkn1=0.6_dl/taurst/SourceAccuracyBoost   
         dkn2=0.9_dl/taurst/SourceAccuracyBoost 
         if (HighAccuracyDefault) dkn2=dkn2/1.2
         if (CP%WantTensors .or. CP%WantVectors) then
              dkn1=dkn1  *0.8_dl
              dlnk0=dlnk0/2 !*0.3_dl
              dkn2=dkn2*0.85_dl
          end if

         qmax_log = dkn1/dlnk0
         q_switch = 2*6.3/taurst 
           !Want linear spacing for wavenumbers which come inside horizon
           !Could use sound horizon, but for tensors that is not relevant

         call Ranges_Init(Evolve_q)
         call Ranges_Add_delta(Evolve_q, qmin, qmax_log, dlnk0, IsLog = .true.)
         call Ranges_Add_delta(Evolve_q, qmax_log, min(qmax,q_switch), dkn1)
         if (qmax > q_switch) then
           call Ranges_Add_delta(Evolve_q, q_switch, qmax, dkn2)
         end if

         call Ranges_GetArray(Evolve_q, .false.)

         if (CP%closed) &
           call SetClosedkValuesFromArr(Evolve_q, .false.) 
 
      end subroutine SetkValuesForSources


      subroutine SetClosedkValuesFromArr(R, forInt)
      Type(Regions) :: R
      integer i,nu,lastnu,nmax
       !nu = 3,4,5... in CP%closed case, so set nearest integers from arr array
      logical, intent(in) :: forInt
      integer ix
      real(dl) dnu
      integer, allocatable :: nu_array(:)
             
       if (forInt .and. nint(R%points(1)*CP%r)<=3) then
        
         !quantization is important       
       call Ranges_Getdpoints(R,half_ends = .false.)
       R%dpoints = max(1,int(R%dpoints*CP%r+0.02))
       lastnu=2
       ix=1
       dnu =R%dpoints(ix)
       nmax=0
       lastnu=2   
       allocate(nu_array(R%npoints*2))    
       do
        
       do while (R%dpoints(ix)==dnu .and. ix <R%npoints) 
        ix=ix+1
       end do
       do nu=lastnu+1,nint(R%points(ix)*CP%r), nint(dnu)
          nmax=nmax+1 
          nu_array(nmax)= nu
       end do
       lastnu=nu_array(nmax) + nint(dnu)-1
       if (ix==R%npoints) exit
       dnu = R%dpoints(ix)
       end do
       if (nint(R%points(R%npoints)*CP%r) > nu_array(nmax)) then
          nmax=nmax+1
          nu_array(nmax) = nint(R%points(R%npoints)*CP%r)
       end if 
       deallocate(R%points)
       allocate(R%points(nmax))
       R%points = nu_array(1:nmax)/CP%r
       deallocate(nu_array)
       
      else
                
       lastnu=3
       nmax=1
      
       do i=2,R%npoints
         nu=nint(R%points(i)*CP%r)
         if (nu > lastnu) then
          nmax=nmax+1 
          lastnu=nu         
          R%points(nmax)=nu/CP%r 
         end if
       
       end do  
       R%points(1)=3/CP%r
       
       end if
       
       R%Lowest = R%points(1) 
       R%Highest = R%points(nmax)
       R%npoints=nmax

      end subroutine SetClosedkValuesFromArr



      subroutine CalcScalarSources(EV,taustart)
      use Transfer
      implicit none
      type(EvolutionVars) EV
      real(dl) tau,tol1,tauend, taustart
      integer j,ind,itf
      real(dl) c(24),w(EV%nvar,9), y(EV%nvar), sources(SourceNum)
      
      real(dl) yprime(EV%nvar), ddelta, delta, adotoa,dtauda, growth
      external dtauda

        if (fixq/=0._dl) then
            !evolution output
            EV%q=fixq
            EV%q2=EV%q**2
        endif

         w=0
         y=0
         call initial(EV,y, taustart)
         if (global_error_flag/=0) return
       
         tau=taustart
         ind=1

!!Example code for plotting out variable evolution
       if (fixq/=0._dl) then
        tol1=tol/exp(AccuracyBoost-1)
        call CreateTxtFile('evolve_q005.txt',1)
         do j=1,1000       
          tauend = taustart+(j-1)*(CP%tau0-taustart)/1000 
          call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
          yprime = 0
          call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)      
          adotoa = 1/(y(1)*dtauda(y(1)))
          ddelta= (yprime(3)*grhoc+yprime(4)*grhob)/(grhob+grhoc)
          delta=(grhoc*y(3)+grhob*y(4))/(grhob+grhoc)
          growth= ddelta/delta/adotoa
          write (1,'(7E15.5)') tau, delta, growth, y(3), y(4), y(EV%g_ix), y(1)
         end do
         close(1)
         stop
      end if

!     Begin timestep loop.

           itf=1
           tol1=tol/exp(AccuracyBoost-1) 
           if (CP%WantTransfer .and. CP%Transfer%high_precision) tol1=tol1/100

           do j=2,TimeSteps%npoints               
             tauend=TimeSteps%points(j)  

             if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > taurend) &
                  .and. .not. CP%Dolensing .and. &
                  (.not.CP%WantTransfer.or.tau > tautf(CP%Transfer%num_redshifts))) then
      
              Src(EV%q_ix,1:SourceNum,j)=0
 
             else
             
             !Integrate over time, calulate end point derivs and calc output
             call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
             if (global_error_flag/=0) return
 
             call output(EV,y,j,tau,sources)
             Src(EV%q_ix,1:SourceNum,j)=sources
            
!     Calculation of transfer functions.
101          if (CP%WantTransfer.and.itf <= CP%Transfer%num_redshifts) then
                if (j < TimeSteps%npoints) then
                  if (tauend < tautf(itf) .and.TimeSteps%points(j+1)  > tautf(itf)) then
                          
                     call GaugeInterface_EvolveScal(EV,tau,y,tautf(itf),tol1,ind,c,w)
                     if (global_error_flag/=0) return

                  endif
                end if  
!     output transfer functions for this k-value.
                      
                  if (abs(tau-tautf(itf)) < 1.e-5_dl) then
                           call outtransf(EV,y, MT%TransferData(:,EV%q_ix,itf))
                         
                           itf=itf+1
                           if (j < TimeSteps%npoints) then
                            if (itf <= CP%Transfer%num_redshifts.and. &
                                TimeSteps%points(j+1) > tautf(itf)) goto 101
                           end if    
                  endif

                  end if
             end if

            end do !time step loop

      end subroutine


      subroutine CalcTensorSources(EV,taustart)

      implicit none
      type(EvolutionVars) EV
      real(dl) tau,tol1,tauend, taustart
      integer j,ind
      real(dl) c(24),wt(EV%nvart,9), yt(EV%nvart)


           call initialt(EV,yt, taustart)
     
           tau=taustart
           ind=1 
           tol1=tol/exp(AccuracyBoost-1)

!     Begin timestep loop.    
           do j=2,TimeSteps%npoints
                  tauend=TimeSteps%points(j)        
                  if (EV%q*tauend > max_etak_tensor) then
                     Src(EV%q_ix,1:SourceNum,j) = 0
                   else
      
                      call GaugeInterface_EvolveTens(EV,tau,yt,tauend,tol1,ind,c,wt)
             
                      call outputt(EV,yt,EV%nvart,j,tau,Src(EV%q_ix,CT_Temp,j),&
                                  Src(EV%q_ix,CT_E,j),Src(EV%q_ix,CT_B,j))
           
                   end if
              end do
    
      end subroutine CalcTensorSources


      subroutine CalcVectorSources(EV,taustart)

      implicit none
      type(EvolutionVars) EV
      real(dl) tau,tol1,tauend, taustart
      integer j,ind
      real(dl) c(24),wt(EV%nvarv,9), yv(EV%nvarv)!,yvprime(EV%nvarv)
        
           
!EV%q=0.2
!EV%q2=EV%q**2

           call initialv(EV,yv, taustart)
     
           tau=taustart
           ind=1 
           tol1=tol*0.01/exp(AccuracyBoost-1)

!!Example code for plotting out variable evolution
!if (.false.) then
!        do j=1,6000      
!          tauend = taustart * exp(j/6000._dl*log(CP%tau0/taustart))
!         call dverk(EV,EV%nvarv,fderivsv,tau,yv,tauend,tol1,ind,c,EV%nvarv,wt) !tauend
!          call fderivsv(EV,EV%nvarv,tau,yv,yvprime)
!
!          write (*,'(7E15.5)') yv(1), yv(2), yv(3),yv(4), &
!                   yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1), &
!                yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2),yv(5)                            
!         end do
!       stop
!nd if

!     Begin timestep loop.            
               do j=2,TimeSteps%npoints
                  tauend=TimeSteps%points(j)

                  if ( EV%q*tauend > max_etak_vector) then
                     Src(EV%q_ix,1:SourceNum,j) = 0
                   else
      
                      call dverk(EV,EV%nvarv,derivsv,tau,yv,tauend,tol1,ind,c,EV%nvarv,wt) !tauend
                 
                      call outputv(EV,yv,EV%nvarv,j,tau,Src(EV%q_ix,CT_Temp,j),&
                                  Src(EV%q_ix,CT_E,j),Src(EV%q_ix,CT_B,j))
           
                   end if
              end do

      end subroutine CalcVectorSources


     subroutine TransferOut
    !Output transfer functions for k larger than used for C_l computation
      implicit none
      integer q_ix
      real(dl) tau
      type(EvolutionVars) EV
    

       if (DebugMsgs .and. Feedbacklevel > 0) & 
         write(*,*) MT%num_q_trans-Evolve_q%npoints, 'transfer k values'

      !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
      !$OMP & PRIVATE(EV, tau, q_ix) 

!     loop over wavenumbers.
         do q_ix=Evolve_q%npoints+1,MT%num_q_trans
      
            EV%TransferOnly=.true. !in case we want to do something to speed it up
          
            EV%q= MT%q_trans(q_ix)

            EV%q2=EV%q**2
            EV%q_ix = q_ix
          
            tau = GetTauStart(EV%q)

            call GetNumEqns(EV)

            call GetTransfer(EV, tau)

         end do
     !$OMP END PARAllEl DO 

     end subroutine TransferOut
      
     subroutine GetTransfer(EV,tau)
       type(EvolutionVars) EV
       real(dl) tau
       integer ind, i
       real(dl) c(24),w(EV%nvar,9), y(EV%nvar)
       real(dl) atol
       
       atol=tol/exp(AccuracyBoost-1)
       if (CP%Transfer%high_precision) atol=atol/10000

       ind=1
       call initial(EV,y, tau) 
       if (global_error_flag/=0) return 
         
       do i=1,CP%Transfer%num_redshifts
          call GaugeInterface_EvolveScal(EV,tau,y,tautf(i),atol,ind,c,w)
          if (global_error_flag/=0) return
          call outtransf(EV,y,MT%TransferData(:,EV%q_ix,i))
       end do
 
     end subroutine GetTransfer


     subroutine NonLinearLensing
        !Scale lensing source terms by non-linear scaling at each redshift and wavenumber
        use NonLinear
        integer i,ik,first_step
        real (dl) tau
        real(dl) scaling(CP%Transfer%num_redshifts), ddScaling(CP%Transfer%num_redshifts)
        real(dl) ho,a0,b0, ascale
        integer tf_lo, tf_hi
        type(MatterPowerData) :: CAMB_Pk

        call Transfer_GetMatterPowerData(MT, CAMB_PK, 1)

        call NonLinear_GetNonLinRatios(CAMB_PK)
 
         if (CP%InitPower%nn > 1) stop 'Non-linear lensing only does one initial power'

        first_step=1
        do while(TimeSteps%points(first_step) < tautf(1))
           first_step = first_step + 1 
        end do
 
        do ik=1, Evolve_q%npoints
          if (Evolve_q%points(ik)/(CP%H0/100) >  Min_kh_nonlinear) then

            !Interpolate non-linear scaling in conformal time
            do i = 1, CP%Transfer%num_redshifts
                scaling(i) = CAMB_Pk%nonlin_ratio(ik,i)
            end do
            if (all(abs(scaling-1) < 5e-4)) cycle
            call spline(tautf(1),scaling(1),CP%Transfer%num_redshifts,&
                                 spl_large,spl_large,ddScaling(1))
       
            tf_lo=1
            tf_hi=tf_lo+1

            do i=first_step,TimeSteps%npoints-1
  
              tau = TimeSteps%points(i)
            
              do while (tau > tautf(tf_hi))
                  tf_lo = tf_lo + 1
                  tf_hi = tf_hi + 1
              end do

              ho=tautf(tf_hi)-tautf(tf_lo) 
              a0=(tautf(tf_hi)-tau)/ho
              b0=1-a0
              
              ascale = a0*scaling(tf_lo)+ b0*scaling(tf_hi)+&
                  ((a0**3-a0)* ddscaling(tf_lo) &
                       +(b0**3-b0)*ddscaling(tf_hi))*ho**2/6
                        
              Src(ik,3,i) = Src(ik,3,i) * ascale
            end  do

          end if
       end do
       
       call MatterPowerdata_Free(CAMB_pk)
 
      end subroutine NonLinearLensing


      subroutine InitSourceInterpolation
      integer i,j
!     get the interpolation matrix for the sources to interpolate them
!     for other k-values
      !$OMP PARAllEl DO DEFAUlT(SHARED), SCHEDUlE(STATIC), PRIVATE(i,j) , SHARED(Evolve_q)
        do  i=1,TimeSteps%npoints
           do j=1, SourceNum
               call spline(Evolve_q%points,Src(1,j,i),Evolve_q%npoints,spl_large,spl_large,ddSrc(1,j,i))
          end do
        end do
       !$OMP END PARAllEl DO 
       end subroutine InitSourceInterpolation


     subroutine SetkValuesForInt
      implicit none
       
      integer no
      real(dl) dk,dk0,dlnk1, dk2, max_k_dk, k_max_log, k_max_0
      integer lognum
      real(dl)  qmax_int,IntSampleBoost
      

       qmax_int = min(qmax,max_bessels_etak/CP%tau0)

       IntSampleBoost=AccuracyBoost 
       if (do_bispectrum) then
        IntSampleBoost = IntSampleBoost * 2  
        if (hard_bispectrum) IntSampleBoost = IntSampleBoost * 2  
       end if
       
!     Fixing the # of k for the integration. 
      
       call Ranges_Init(ThisCT%q)

      if (CP%closed.and.ExactClosedSum) then
       
        call Ranges_Add(ThisCT%q,3/CP%r, nint(qmax_int*CP%r)/CP%r, nint(qmax_int*CP%r)-3) !fix jun08
        call Init_ClTransfer(ThisCT)
        call Ranges_Getdpoints(ThisCT%q,half_ends = .false.) !Jun08
      else

      !Split up into logarithmically spaced intervals from qmin up to k=lognum*dk0
      !then no-lognum*dk0 linearly spaced at dk0 up to no*dk0
      !then at dk up to max_k_dk, then dk2 up to qmax_int
         lognum=nint(10*IntSampleBoost)
         dlnk1=1._dl/lognum  
         no=nint(600*IntSampleBoost)
         dk0=1.8_dl/CP%r/CP%chi0/IntSampleBoost   
         dk=3._dl/CP%r/CP%chi0/IntSampleBoost 
         
         if (HighAccuracyDefault) dk=dk/1.6 

         k_max_log = lognum*dk0
         k_max_0  = no*dk0
         
         if (do_bispectrum) k_max_0 = max(10.d0,k_max_0) 

         dk2 = 0.04/IntSampleBoost  !very small scales  

         call Ranges_Add_delta(ThisCT%q, qmin, k_max_log, dlnk1, IsLog = .true.)
         call Ranges_Add_delta(ThisCT%q, k_max_log, min(qmax_int,k_max_0), dk0)     
       
         if (qmax_int > k_max_0) then

         max_k_dk = max(3000, 2*maximum_l)/CP%tau0
       
          call Ranges_Add_delta(ThisCT%q, k_max_0, min(qmax_int, max_k_dk), dk)
          if (qmax_int > max_k_dk) then
           !This allows inclusion of high k modes for computing BB lensed spectrum accurately
           !without taking ages to compute.
            call Ranges_Add_delta(ThisCT%q, max_k_dk, qmax_int, dk2)
          end if    

         end if           

         call Init_ClTransfer(ThisCT)

       if (CP%closed) then
        call SetClosedkValuesFromArr(ThisCT%q,.true.)
        call Ranges_Getdpoints(ThisCT%q,half_ends = .false.)
        ThisCT%q%dpoints(1) = 1/CP%r
        !!!  
        deallocate(ThisCT%Delta_p_l_k) !Re-do this from Init_ClTransfer because number of points changed
        allocate(ThisCT%Delta_p_l_k(ThisCT%NumSources,min(max_bessels_l_index,ThisCT%ls%l0), ThisCT%q%npoints))
        ThisCT%Delta_p_l_k = 0     
       end if

       end if !ExactClosedSum

         
      end subroutine setkValuesForInt

      subroutine InterpolateSources(IV)
      implicit none
      integer i,khi,klo, step
      real(dl) xf,b0,ho,a0,ho2o6,a03,b03
      type(IntegrationVars) IV


!     finding position of k in table Evolve_q to do the interpolation.

!Can't use the following in closed case because regions are not set up (only points)  
!           klo = min(Evolve_q%npoints-1,Ranges_IndexOf(Evolve_q, IV%q))
           !This is a bit inefficient, but thread safe
            klo=1
            do while ((IV%q > Evolve_q%points(klo+1)).and.(klo < (Evolve_q%npoints-1))) 
               klo=klo+1
            end do

            khi=klo+1


            ho=Evolve_q%points(khi)-Evolve_q%points(klo)
            a0=(Evolve_q%points(khi)-IV%q)/ho
            b0=(IV%q-Evolve_q%points(klo))/ho           
            ho2o6 = ho**2/6
            a03=(a0**3-a0)
            b03=(b0**3-b0)
            IV%SourceSteps = 0

!     Interpolating the source as a function of time for the present
!     wavelength.
             step=2
               do i=2, TimeSteps%npoints
                  xf=IV%q*(CP%tau0-TimeSteps%points(i))
                  if (CP%WantTensors) then
                   if (IV%q*TimeSteps%points(i) < max_etak_tensor.and. xf > 1.e-8_dl) then
                    step=i
                     IV%Source_q(i,1:SourceNum) =a0*Src(klo,1:SourceNum,i)+&
                          b0*Src(khi,1:SourceNum,i)+(a03 *ddSrc(klo,1:SourceNum,i)+ &
                          b03*ddSrc(khi,1:SourceNum,i)) *ho2o6
                    else
                     IV%Source_q(i,1:SourceNum) = 0
                   end if
                  end if
                  if (CP%WantVectors) then
                   if (IV%q*TimeSteps%points(i) < max_etak_vector.and. xf > 1.e-8_dl) then
                     step=i
                     IV%Source_q(i,1:SourceNum) =a0*Src(klo,1:SourceNum,i)+&
                          b0*Src(khi,1:SourceNum,i)+(a03 *ddSrc(klo,1:SourceNum,i)+ &
                          b03*ddSrc(khi,1:SourceNum,i)) *ho2o6
                    else
                     IV%Source_q(i,1:SourceNum) = 0
                   end if
                  end if

                  if (CP%WantScalars) then
                     if ((DebugEvolution .or. CP%Dolensing .or. IV%q*TimeSteps%points(i) < max_etak_scalar) &
                          .and. xf > 1.e-8_dl) then
                        step=i
                        IV%Source_q(i,1:SourceNum)=a0*Src(klo,1:SourceNum,i)+ &
                         b0*Src(khi,1:SourceNum,i) + (a03*ddSrc(klo,1:SourceNum,i) &
                         +b03*ddSrc(khi,1:SourceNum,i))*ho2o6
  
                     else
                       IV%Source_q(i,1:SourceNum) = 0 
                     end if
                  end if
               end do
               IV%SourceSteps = step


          if (.not.CP%flat) then
             do i=1, SourceNum
                call spline(TimeSteps%points,IV%Source_q(1,i),TimeSteps%npoints,&
                     spl_large,spl_large,IV%ddSource_q(1,i))
             end do
          end if
           
           IV%SourceSteps = IV%SourceSteps*1
           !This is a fix for a compiler bug on Seaborg 

      end subroutine

      
      subroutine IntegrationVars_Init(IV)
      type(IntegrationVars), intent(INOUT) :: IV
        
        IV%Source_q(1,1:SourceNum)=0
        IV%Source_q(TimeSteps%npoints,1:SourceNum) = 0
        IV%Source_q(TimeSteps%npoints-1,1:SourceNum) = 0
 
      end  subroutine IntegrationVars_Init


!cccccccccccccccccccccccccccccccccccc


      subroutine DoSourceIntegration(IV) !for particular wave number q
      integer j,ll,llmax
      real(dl) nu
      type(IntegrationVars) IV
       
         
         nu=IV%q*CP%r
        
         if (CP%closed) then
          
          if (nu<20 .or. CP%tau0/CP%r+6*pi/nu > pi/2) then
           llmax=nint(nu)-1
          else
           llmax=nint(nu*rofChi(CP%tau0/CP%r + 6*pi/nu))
           llmax=min(llmax,nint(nu)-1)  !nu >= l+1
          end if       

         else
          llmax=nint(nu*CP%chi0)
          
          if (llmax<15) then
           llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
          else
           llmax=nint(nu*rofChi(CP%tau0/CP%r + 6*pi/nu))
          end if
          
      
         end if 

         if (CP%flat) then
            call DoFlatIntegration(IV,llmax)
         else
           do j=1,lSamp%l0
            ll=lSamp%l(j)
            if (ll>llmax) exit  
         
            call IntegrateSourcesBessels(IV,j,ll,nu)      
           end do !j loop
         end if

       

      end subroutine DoSourceIntegration     

      function UseLimber(l,k)
       !Calculate lensing potential power using Limber rather than j_l integration
       !even when sources calculated as part of temperature calculation
       !(Limber better on small scales unless step sizes made much smaller)
       !This affects speed, esp. of non-flat case
        logical :: UseLimber
        integer l
        real(dl) :: k
 
        !note increasing non-limber to l>700 is not neccessarily more accurate unless AccruacyBoost much higher
        !use **0.2 to at least give some sensitivity to Limber effects
        if (CP%AccurateBB .or. CP%flat) then
         UseLimber = l > 700*AccuracyBoost**0.2 .and. k > 0.05
        else
         !This is accurate at percent level only (good enough here)
         UseLimber = l > 300*min(AccuracyBoost,2.4_dl) .or. k>0.05
        end if

      end function UseLimber

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!flat source integration
        subroutine DoFlatIntegration(IV, llmax)
        implicit none
        type(IntegrationVars) IV
        integer llmax
        integer j
        logical DoInt
        real(dl) xlim,xlmax1 
        real(dl) tmin, tmax
        real(dl) a2, J_l, aa(IV%SourceSteps), fac(IV%SourceSteps)
        real(dl) xf, sums(SourceNum)
        real(dl) qmax_int
        integer bes_ix,n, bes_index(IV%SourceSteps)
     
!     Find the position in the xx table for the x correponding to each
!     timestep

         do j=1,IV%SourceSteps !Precompute arrays for this k
            xf=abs(IV%q*(CP%tau0-TimeSteps%points(j)))
            bes_index(j)=Ranges_indexOf(BessRanges,xf)
          !Precomputed values for the interpolation
            bes_ix= bes_index(j)
            fac(j)=BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
            aa(j)=(BessRanges%points(bes_ix+1)-xf)/fac(j)
            fac(j)=fac(j)**2*aa(j)/6
         end do

         do j=1,lSamp%l0 
             if (lSamp%l(j) > llmax) return
             xlim=xlimfrac*lSamp%l(j)
             xlim=max(xlim,xlimmin)
             xlim=lSamp%l(j)-xlim
             if (full_bessel_integration .or. do_bispectrum) then
                 tmin = TimeSteps%points(2)
             else
                 xlmax1=80*lSamp%l(j)*AccuracyBoost 
                 tmin=CP%tau0-xlmax1/IV%q
                 tmin=max(TimeSteps%points(2),tmin)                 
             end if 
             tmax=CP%tau0-xlim/IV%q
             tmax=min(CP%tau0,tmax)
        
             if (tmax < TimeSteps%points(2)) exit
             sums(1:SourceNum) = 0
 
            !As long as we sample the source well enough, it is sufficient to
            !interpolate the Bessel functions only

             if (SourceNum==2) then
              !This is the innermost loop, so we separate the no lensing scalar case to optimize it
                 do n= Ranges_IndexOf(TimeSteps,tmin),min(IV%SourceSteps,Ranges_IndexOf(TimeSteps,tmax))

                 a2=aa(n)
                 bes_ix=bes_index(n) 

                 J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                        *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline
 
                 J_l = J_l*TimeSteps%dpoints(n)
                 sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                 sums(2) = sums(2) + IV%Source_q(n,2)*J_l

                end do
              else 
                 qmax_int= max(850,lSamp%l(j))*3*AccuracyBoost/CP%tau0
                 if (HighAccuracyDefault) qmax_int=qmax_int*1.2
                 DoInt = .not. CP%WantScalars .or. IV%q < qmax_int 
                 if (DoInt) then
                  do n= Ranges_IndexOf(TimeSteps,tmin),min(IV%SourceSteps,Ranges_IndexOf(TimeSteps,tmax))
                  !Full Bessel integration
                     a2=aa(n)
                     bes_ix=bes_index(n) 

                     J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                            *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline
                     J_l = J_l*TimeSteps%dpoints(n)
              
                    !The unwrapped form is faster
     
                     sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                     sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                     sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                     
                  end do
                 end if
                 if (.not. DoInt .or. UseLimber(lsamp%l(j),IV%q) .and. CP%WantScalars) then
                  !Limber approximation for small scale lensing (better than poor version of above integral)
                  xf = CP%tau0-lSamp%l(j)/IV%q
                  if (xf < TimeSteps%Highest .and. xf > TimeSteps%Lowest) then
                   n=Ranges_IndexOf(TimeSteps,xf)
                   xf= (xf-TimeSteps%points(n))/(TimeSteps%points(n+1)-TimeSteps%points(n))                  
                   sums(3) = (IV%Source_q(n,3)*(1-xf) + xf*IV%Source_q(n+1,3))*sqrt(pi/2/lSamp%l(j))/IV%q 
                  else
                   sums(3)=0
                  end if
                 end if

              end if
  
              ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) = ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) + sums(1:SourceNum)
 !             IV%Delta_l_q(1:SourceNum,j) = IV%Delta_l_q(1:SourceNum,j) + sums(1:SourceNum)
         
          end do

      
        end subroutine DoFlatIntegration
   

    
!non-flat source integration

      subroutine IntegrateSourcesBessels(IV,j,l,nu)  
      use SpherBessels
      type(IntegrationVars) IV
      logical DoInt   
      integer l,j, nstart,nDissipative,ntop,nbot,nrange,nnow
      real(dl) nu,ChiDissipative,ChiStart,tDissipative,y1,y2,y1dis,y2dis     
      real(dl) xf,x,chi, miny1
      real(dl) sums(SourceNum),out_arr(SourceNum), qmax_int   
    
      !Calculate chi where for smaller chi it is dissipative
      x=sqrt(real(l*(l+1),dl))/nu
    
      ChiDissipative=invsinfunc(x) 
    
      ChiStart=ChiDissipative
      !Move down a bit to get smaller value (better accuracy integrating up from small values)
      if (nu<300) ChiStart = max(ChiDissipative-1._dl/nu,1d-6)   !max(ChiDissipative-1._dl/nu,1d-6)
 
        !Then get nearest source point with lower Chi...
      tDissipative=CP%tau0 - CP%r*ChiStart
      if (tDissipative<TimeSteps%points(1)) then
         nDissipative=2
       else 
         nDissipative = Ranges_IndexOf(TimeSteps,tDissipative)+1
      endif 
      nDissipative=min(nDissipative,TimeSteps%npoints-1)

      tDissipative = TimeSteps%points(nDissipative)

      ChiStart =  max(1d-8,(CP%tau0-tDissipative)/CP%r)
       
      !Get values at ChiStart
   
      call USpherBesselWithDeriv(CP%closed,ChiStart,l,nu,y1dis,y2dis)  

     nstart=nDissipative         
     chi=ChiStart

     if ((CP%WantScalars)) then !Do Scalars

      !Integrate chi down in dissipative region
      ! cuts off when ujl gets small
         miny1= 0.5d-4/l/AccuracyBoost
         sums=0
         qmax_int= max(850,lSamp%l(j))*3*AccuracyBoost/(CP%chi0*CP%r)
         if (HighAccuracyDefault) qmax_int=qmax_int*1.2 
         DoInt =  SourceNum/=3 .or. IV%q < qmax_int
         if (DoInt) then
         if ((nstart < min(TimeSteps%npoints-1,IV%SourceSteps)).and.(y1dis > miny1)) then
     
            y1=y1dis
            y2=y2dis
            nnow=nstart
            do nrange = 1,TimeSteps%Count
               if (nrange == TimeSteps%count) then
                ntop = TimeSteps%npoints -1
               else
                ntop = TimeSteps%R(nrange+1)%start_index
               end if
               if (nnow < ntop) then
                  call DoRangeInt(IV,chi,ChiDissipative,nnow,ntop,TimeSteps%R(nrange)%delta, &
                              nu,l,y1,y2,out_arr)        
                 sums  = sums + out_arr
                 nnow = ntop
                 if (chi==0) exit !small enough to cut off
               end if
             end do
           
           end if !integrate down chi
             
         !Integrate chi up in oscillatory region
          if (nstart > 2) then
           y1=y1dis
           y2=y2dis
           chi=ChiStart
            nnow=nstart
            do nrange = TimeSteps%Count,1,-1 
               nbot = TimeSteps%R(nrange)%start_index          
               if (nnow >  nbot) then
                  call DoRangeInt(IV,chi,ChiDissipative,nnow,nbot,TimeSteps%R(nrange)%delta, &
                              nu,l,y1,y2,out_arr)
                 sums=sums+out_arr
         
                 if (chi==0) exit !small for remaining region
                 nnow = nbot
               end if
               
            end do

           end if

           end if !DoInt
         if (SourceNum==3 .and. (.not. DoInt .or. UseLimber(l,IV%q))) then
            !Limber approximation for small scale lensing (better than poor version of above integral)
             xf = CP%tau0-invsinfunc(l/nu)*CP%r

!Feb09 fix screw up introduced Feb 2008 version
!       if (xf > TimeSteps%Lowest .and. xf > TimeSteps%Highest) then
             if (xf < TimeSteps%Highest .and. xf > TimeSteps%Lowest) then
             nbot=Ranges_IndexOf(TimeSteps,xf)
             xf= (xf-TimeSteps%points(nbot))/(TimeSteps%points(nbot+1)-TimeSteps%points(nbot))                  
             sums(3) = (IV%Source_q(nbot,3)*(1-xf) + xf*IV%Source_q(nbot+1,3))*&
                           sqrt(pi/2/l/sqrt(1-CP%Ksign*real(l**2)/nu**2))/IV%q
             else
              sums(3) = 0
             end if                
         end if

         ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix)=ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix)+sums
           
      end if !Do Scalars
           
     if ((CP%WantTensors)) then !Do Tensors
         chi=ChiStart
        
      !Integrate chi down in dissipative region
      !DoRangeInt cuts off when ujl gets small
         miny1= 1.d-6/l/AccuracyBoost
         if ((nstart < TimeSteps%npoints-1).and.(y1dis>miny1)) then
            y1=y1dis
            y2=y2dis
            nnow=nstart
            do nrange = 1,TimeSteps%Count
               if (nrange == TimeSteps%count) then
                ntop = TimeSteps%npoints -1
               else
                ntop = TimeSteps%R(nrange+1)%start_index
               end if
               if (nnow < ntop) then
                 call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,ntop,TimeSteps%R(nrange)%delta, &
                              nu,l,y1,y2,out_arr)

                 ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) = ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) + out_arr
               
                 nnow = ntop
                 if (chi==0) exit
               end if
            end do
                   
         end if 
        
                 
!Integrate chi up in oscillatory region
          if (nstart > 2) then
           y1=y1dis
           y2=y2dis
           chi=ChiStart
         
           nnow=nstart
            do nrange = TimeSteps%Count,1,-1 
               nbot = TimeSteps%R(nrange)%start_index          
               if (nnow >  nbot) then
                 call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,nbot,TimeSteps%R(nrange)%delta, &
                              nu,l,y1,y2,out_arr)
                 ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) = ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) + out_arr
                
                 nnow = nbot
                 if (chi==0) exit !small for remaining region
               end if               
            end do

          end if
       
        end if !Do Tensors
     
      end subroutine IntegrateSourcesBessels

   

 subroutine DoRangeInt(IV,chi,chiDisp,nstart,nend,dtau,nu,l,y1,y2,out)
 !Non-flat version

!returns chi at end of integral (where integral stops, not neccessarily end)
! This subroutine integrates the source*ujl for steps nstart to nend
! It calculates ujl by integrating a second order
! differential equation from initial values.
! dtau is the spacing of the timesteps (they must be equally spaced)

      use precision
      use ModelParams
      type(IntegrationVars) IV
      integer l,nIntSteps,nstart,nend,nlowest,isgn,i,is,Startn
      real(dl) nu,dtau,num1,num2,Deltachi,aux1,aux2
      real(dl) a,b,tmpa,tmpb,hh,h6,xh,delchi,taui
      real(dl) nu2,chi,chiDisp,dydchi1,dydchi2,yt1,yt2,dyt1,dyt2,dym1,dym2
  
      real(dl) tmp,dtau2o6,y1,y2,ap1,sh,ujl,chiDispTop
      real(dl) dchimax,dchisource,sgn,sgndelchi,minujl
      real(dl), parameter:: MINUJl1 = 0.5d-4  !cut-off point for small ujl l=1
      logical Interpolate
      real(dl) scalel
      real(dl) IntAccuracyBoost
      real(dl) sources(SourceNum), out(SourceNum)
   
      IntAccuracyBoost=AccuracyBoost 

! atau0 is the array with the time where the sources are stored.
      if (nend==nstart) then  
            out = 0
            return
         end if
      
      dchisource=dtau/CP%r
     
      num1=1._dl/nu     

      scalel=l/scale
      if (scalel>=2400) then
         num2=num1*2.5
      else if (scalel< 50) then
         num2=num1*0.8_dl
      else 
        num2=num1*1.5_dl 
      end if    
      !Dec 2003, since decrease dtaurec, can make this smaller
      if (dtau==dtaurec_q) then
       num2=num2/4
      end if

      if (HighAccuracyDefault .and. scalel<1500 .and. scalel > 150) &
          IntAccuracyBoost=IntAccuracyBoost*(1+(2000-scalel)*0.6/2000 ) 

      if (num2*IntAccuracyBoost < dchisource .and. (.not. CP%DoLensing .or. UseLimber(l,IV%q)) & 
!       if ((num2*IntAccuracyBoost < dchisource ) & !Oscillating fast 
        .or. (nstart>IV%SourceSteps.and.nend>IV%SourceSteps)) then  
         out = 0
         y1=0._dl !So we know to calculate starting y1,y2 if there is next range
         y2=0._dl
         chi=(CP%tau0-TimeSteps%points(nend))/CP%r
         return
        end if
     
      Startn=nstart
      if (nstart>IV%SourceSteps .and. nend < IV%SourceSteps) then
         chi=(CP%tau0-TimeSteps%points(IV%SourceSteps))/CP%r
         Startn=IV%SourceSteps
         call USpherBesselWithDeriv(CP%closed,chi,l,nu,y1,y2)
      else if ((y2==0._dl).and.(y1==0._dl)) then 
         call USpherBesselWithDeriv(CP%closed,chi,l,nu,y1,y2)
      end if
    
      if (CP%closed) then
        !Need to cut off when ujl gets exponentially small as it approaches Pi
         chiDispTop = pi - chiDisp
      else
         chiDispTop = 1d20
      end if

      minujl=MINUJl1/l/IntAccuracyBoost 
      isgn=sign(1,Startn-nend)!direction of chi integration 
        !higher n, later time, smaller chi
    
      sgn= isgn

      nlowest=min(Startn,nend)
      aux1=1._dl*CP%r/dtau  !used to calculate nearest timestep quickly
      aux2=(CP%tau0-TimeSteps%points(nlowest))/dtau + nlowest
          
      nu2=nu*nu
      ap1=l*(l+1)
      sh=rofChi(chi)
           
      if (scalel < 1100) then
         dchimax= 0.3*num1 
      else if (scalel < 1400) then
         dchimax=0.25_dl*num1 *1.5
      else
         dchimax=0.35_dl*num1 *1.5
      end if

      dchimax=dchimax/IntAccuracyBoost 
    
      ujl=y1/sh
      sources = IV%Source_q(Startn,1:SourceNum)

      out = 0.5_dl*ujl*sources

      Interpolate = dchisource > dchimax
      if (Interpolate) then !split up smaller than source step size
         delchi=dchimax
         Deltachi=sgn*(TimeSteps%points(Startn)-TimeSteps%points(nend))/CP%r
         nIntSteps=int(Deltachi/delchi+0.99_dl)        
         delchi=Deltachi/nIntSteps 
         dtau2o6=(CP%r*delchi)**2/6._dl
       else !step size is that of source
         delchi=dchisource
         nIntSteps=isgn*(Startn-nend)      
       end if
        
         sgndelchi=delchi*sgn
         tmp=(ap1/sh**2 - nu2) 
         hh=0.5_dl*sgndelchi  
         h6=sgndelchi/6._dl

      
          do i=1,nIntSteps          
! One step in the ujl integration
! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2     
            xh=chi+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/rofChi(xh)**2 - nu2) 
            
            
            dyt2=tmp*yt1       !deriv y2 at mid
       
            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid
           
            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+sgndelchi*dym1 !y1 at end
            dym1=dyt1+dym1     
            yt2=y2+sgndelchi*dym2 !y2 at end
            dym2=dyt2+dym2
            
            chi=chi+sgndelchi     !end point
            sh=rofChi(chi)    
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2*dym2)       

            ujl=y1/sh
            if ((isgn<0).and.(y1*y2<0._dl).or.((chi>chiDispTop).and.((chi>3.14).or.(y1*y2>0)))) then
                chi=0._dl 
                exit   !If this happens we are small, so stop integration
            end if

     
            if (Interpolate) then
! Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
             b=taui-is
    
            if (b > 0.998) then 
               !may save time, and prevents numerical error leading to access violation of IV%Source_q(0)
             sources = IV%Source_q(is+1,1:SourceNum)
             else

             a=1._dl-b          
             tmpa=(a**3-a)
             tmpb=(b**3-b)
             sources=a*IV%Source_q(is,1:SourceNum)+b*IV%Source_q(is+1,1:SourceNum)+ &
                  (tmpa*IV%ddSource_q(is,1:SourceNum)+ &
                   tmpb*IV%ddSource_q(is+1,1:SourceNum))*dtau2o6
             end if

            else
             sources = IV%Source_q(Startn - i*isgn,1:SourceNum)
      
            end if

            out = out + ujl*sources
                
            if (((isgn<0).or.(chi>chiDispTop)).and.(abs(ujl) < minujl)) then
           
             chi=0
            exit !break when getting  exponentially small in dissipative region

            end if
            
         end do

         out = (out - sources*ujl/2)*delchi*CP%r
         
         end subroutine DoRangeInt

 


      subroutine DoRangeIntTensor(IV,chi,chiDisp,nstart,nend,dtau,nu,l,y1,y2,out)
! It calculates ujl by integrating a second order
! differential equation from initial values for calculating ujl.
! nstart and nend are the starting and finishing values of the
! integration.
! dtau is the spacing of the timesteps (they must be equally spaced)

      use precision
      use ModelParams
      type(IntegrationVars), target :: IV
      integer l,nIntSteps,nstart,nend,nlowest,isgn,i,is
      real(dl) nu,dtau,num1,num2,Deltachi,aux1,aux2
      real(dl) a,b,tmpa,tmpb,hh,h6,xh,delchi,taui,scalel
      real(dl) nu2,chi,chiDisp,chiDispTop
      real(dl) dydchi1,dydchi2,yt1,yt2,dyt1,dyt2,dym1,dym2
  
      real(dl) tmp,dtau2o6,y1,y2,ap1,sh,ujl
      real(dl) dchimax,dchisource,sgn,sgndelchi,minujl
      real(dl), parameter:: MINUJl1 = 1.D-6  !cut-off point for smal ujl l=1
      logical Interpolate
      real(dl) out(SourceNum), source(SourceNum)
      real(dl), dimension(:,:), pointer :: sourcep, ddsourcep

      sourcep => IV%Source_q(:,1:)
      ddsourcep => IV%ddSource_q(:,1:)
      
     
      if (nend==nstart) then  
            out=0
            return
         end if    
      minujl=MINUJL1*AccuracyBoost/l
      isgn=sign(1,nstart-nend)!direction of chi integration 
        !higher n, later time, smaller chi

      if (CP%closed) then
        !Need to cut off when ujl gets exponentially small as it approaches Pi
         chiDispTop = pi - chiDisp
      else
         chiDispTop = 1d20
      end if
      
     
      num1=1._dl/nu
      dchisource=dtau/CP%r

      scalel=l/scale
      if (scalel>=2000) then
         num2=num1*4
      else if (scalel>=1000) then
         num2=num1*2.5 
      else if (scalel< 75) then
         num2=num1*0.1_dl
      else if (scalel<180) then
         num2=num1*0.3_dl
      else if (scalel < 600) then
         num2=num1*0.8_dl
      else
         num2=num1
      end if
  
      if ((isgn==1).and.(num2*AccuracyBoost < dchisource)) then  !Oscillating fast
         out = 0
         y1=0._dl !!So we know to calculate starting y1,y2 if there is next range
         y2=0._dl
         chi=(CP%tau0-TimeSteps%points(nend))/CP%r
         return
        end if
      if ((y2==0._dl).and.(y1==0._dl)) call USpherBesselWithDeriv(CP%closed,chi,l,nu,y1,y2)
   
      sgn=isgn

      nlowest=min(nstart,nend)
      aux1=1._dl*CP%r/dtau  !used to calculate nearest timestep quickly
      aux2=(CP%tau0-TimeSteps%points(nlowest))/dtau + nlowest
          
     
      nu2=nu*nu
      ap1=l*(l+1)

      sh=rofChi(chi)
      
      if (scalel < 120) then
         dchimax=0.6_dl*num1
      else if (scalel < 1400) then
         dchimax=0.25_dl*num1
      else
         dchimax=0.35_dl*num1 
      end if

      dchimax=dchimax/AccuracyBoost
     
      ujl=y1/sh
      out = ujl * sourcep(nstart,1:SourceNum)/2
   
      Interpolate = dchisource > dchimax
      if (Interpolate) then !split up smaller than source step size
         delchi=dchimax
         Deltachi=sgn*(TimeSteps%points(nstart)-TimeSteps%points(nend))/CP%r
         nIntSteps=int(Deltachi/delchi+0.99_dl)
         delchi=Deltachi/nIntSteps 
         dtau2o6=(CP%r*delchi)**2/6._dl
       else !step size is that of source
         delchi=dchisource
         nIntSteps=isgn*(nstart-nend)      
       end if
    
       
         sgndelchi=delchi*sgn
         tmp=(ap1/sh**2 - nu2) 
         hh=0.5_dl*sgndelchi  
         h6=sgndelchi/6._dl
    
                
         do i=1,nIntSteps
            

! One step in the ujl integration
! fourth-order Runge-Kutta method to integrate equation for ujl

            dydchi1=y2         !deriv y1
            dydchi2=tmp*y1     !deriv y2     
            xh=chi+hh          !midpoint of step
            yt1=y1+hh*dydchi1  !y1 at midpoint
            yt2=y2+hh*dydchi2  !y2 at midpoint
            dyt1=yt2           !deriv y1 at mid
            tmp=(ap1/rofChi(xh)**2 - nu2) 
          
            
            dyt2=tmp*yt1       !deriv y2 at mid        
            yt1=y1+hh*dyt1     !y1 at mid
            yt2=y2+hh*dyt2     !y2 at mid
           
            dym1=yt2           !deriv y1 at mid
            dym2=tmp*yt1       !deriv y2 at mid
            yt1=y1+sgndelchi*dym1 !y1 at end
            dym1=dyt1+dym1     
            yt2=y2+sgndelchi*dym2 !y2 at end
            dym2=dyt2+dym2
         
            chi=chi+sgndelchi     !end point
            sh=rofChi(chi)    
            dyt1=yt2           !deriv y1 at end
            tmp=(ap1/sh**2 - nu2)
            dyt2=tmp*yt1       !deriv y2 at end
            y1=y1+h6*(dydchi1+dyt1+2._dl*dym1) !add up
            y2=y2+h6*(dydchi2+dyt2+2._dl*dym2)

            ujl=y1/sh 
            if ((isgn<0).and.(y1*y2<0._dl).or.((chi>chiDispTop).and.((chi>3.14).or.(y1*y2>0)))) then
                chi=0._dl
                exit   !exit because ujl now small
                end if
            
            if (Interpolate) then
! Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
            b=taui-is
            if (b > 0.995) then 
               !may save time, and prevents numerical error leading to access violation of zero index
             is=is+1
             source = sourcep(is,1:SourceNum)
            
             else
             a=1._dl-b            
             tmpa=(a**3-a)
             tmpb=(b**3-b)
             source = a*sourcep(is,1:SourceNum)+b*sourcep(is+1,1:SourceNum)+ &
                   (tmpa*ddsourcep(is,1:SourceNum) +  tmpb*ddsourcep(is+1,1:SourceNum))*dtau2o6

             end if           

            else
             source = sourcep(nstart - i*isgn,1:SourceNum)
          
            end if
            out = out + source * ujl
                  
            if (((isgn<0).or.(chi>chiDispTop)).and.(abs(ujl) < minujl)) then
            chi=0
            exit  !break when getting  exponentially small in dissipative region
            end if
         end do

         out = (out - source * ujl /2)*delchi*CP%r
       
         end subroutine DoRangeIntTensor

        subroutine GetInitPowerArrayVec(pows,ks, numks,pix)
         integer, intent(in) :: numks, pix
         real(dl) pows(numks), ks(numks)
         integer i
      
         do i = 1, numks
         !!change to vec...
            pows(i) =  ScalarPower(ks(i) ,pix)
            if (global_error_flag/=0) exit
         end do

        end subroutine GetInitPowerArrayVec


        subroutine GetInitPowerArrayTens(pows,ks, numks,pix)
         integer, intent(in) :: numks, pix
         real(dl) pows(numks), ks(numks)
         integer i
      
         do i = 1, numks
            pows(i) =  TensorPower(ks(i) ,pix)
            if (global_error_flag/=0) exit
         end do

        end subroutine GetInitPowerArrayTens


        subroutine CalcScalCls(CTrans)
        use Bispectrum
        implicit none
        Type(ClTransferData) :: CTrans
        integer pix,j
        real(dl) apowers, pows(CTrans%q%npoints)
        integer q_ix
        real(dl)  ks(CTrans%q%npoints),dlnks(CTrans%q%npoints),dlnk
        real(dl) ctnorm,dbletmp


         do pix=1,CP%InitPower%nn

          do q_ix = 1, CTrans%q%npoints 

             if (CP%flat) then
                     ks(q_ix) = CTrans%q%points(q_ix)
                     dlnks(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
             else
                     ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - CP%curv)
                     dlnks(q_ix) = CTrans%q%dpoints(q_ix)*CTrans%q%points(q_ix)/ks(q_ix)**2
             end if

             pows(q_ix) =  ScalarPower(ks(q_ix) ,pix)
             if (global_error_flag/=0) return
      
          end do


        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,dlnk,apowers,ctnorm,dbletmp)

         do j=1,CTrans%ls%l0

        !Integrate dk/k Delta_l_q**2 * Power(k)
          do q_ix = 1, CTrans%q%npoints 

             if (.not.(CP%closed.and.nint(CTrans%q%points(q_ix)*CP%r)<=CTrans%ls%l(j))) then 
               !cut off at nu = l + 1
             dlnk = dlnks(q_ix)
             apowers = pows(q_ix)

             iCl_scalar(j,C_Temp:C_E,pix) = iCl_scalar(j,C_Temp:C_E,pix) +  &
                          apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)**2*dlnk
             iCl_scalar(j,C_Cross,pix) = iCl_scalar(j,C_Cross,pix) + &
                          apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
             if (CTrans%NumSources>2) then
                        iCl_scalar(j,C_Phi,pix) = iCl_scalar(j,C_Phi,pix) +  &
                                                       apowers*CTrans%Delta_p_l_k(3,j,q_ix)**2*dlnk
                        iCl_scalar(j,C_PhiTemp,pix) = iCl_scalar(j,C_PhiTemp,pix) +  &
                                          apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
                        iCl_scalar(j,C_PhiE,pix) = iCl_scalar(j,C_PhiE,pix) +  &
                                          apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
             end if

             end if
             
           end do

!Output l(l+1)C_l/OutputDenominator

           !ctnorm = (CTrans%ls%l+2)!/(CTrans%ls%l-2)! - beware of int overflow
            ctnorm=(CTrans%ls%l(j)*CTrans%ls%l(j)-1)*real((CTrans%ls%l(j)+2)*CTrans%ls%l(j),dl)
            dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*fourpi  
                 
            iCl_scalar(j,C_Temp,pix)  =  iCl_scalar(j,C_Temp,pix)*dbletmp
            iCl_scalar(j,C_E,pix)     =  iCl_scalar(j,C_E,pix)*dbletmp*ctnorm
            iCl_scalar(j,C_Cross,pix) =  iCl_scalar(j,C_Cross,pix)*dbletmp*sqrt(ctnorm)
            if (CTrans%NumSources>2) then
                     iCl_scalar(j,C_Phi,pix)   =  &
                          ALens*iCl_scalar(j,C_Phi,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)**2    
                     !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
                     !We put pix extra factors of l here to improve interpolation in CTrans%ls%l
                     iCl_scalar(j,C_PhiTemp,pix)   =  &
                          sqrt(ALens)*  iCl_scalar(j,C_PhiTemp,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)*CTrans%ls%l(j)
                      !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
                     iCl_scalar(j,C_PhiE,pix)   =  &
                          sqrt(ALens)*  iCl_scalar(j,C_PhiE,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)*CTrans%ls%l(j)*sqrt(ctnorm)
                      !Cross-correlation is CTrans%ls%l^3 C_l^{\phi E}
             end if

           end do
         !$OMP END PARAllEl DO

          end do

        end subroutine CalcScalCls

        subroutine CalcScalCls2(CTrans)
        !Calculate C_ll' for non-isotropic models
        !Run with l_sample_boost=50 to get every l
        !not used in normal CAMB
        implicit none
        Type(ClTransferData) :: CTrans
        integer j,j2,in
        real(dl) apowers, pows(CTrans%q%npoints)
        integer q_ix
        real(dl)  ks(CTrans%q%npoints),dlnks(CTrans%q%npoints),dlnk
        real(dl) ctnorm,dbletmp
        real(dl), allocatable :: iCl_Scalar2(:,:,:,:) 

         allocate(iCl_Scalar2(CTranS%ls%l0,CTrans%ls%l0,C_Temp:C_last,CP%InitPower%nn))
         iCl_scalar2 = 0
    
         do in=1,CP%InitPower%nn
          do q_ix = 1, CTrans%q%npoints 

             if (CP%flat) then
                     ks(q_ix) = CTrans%q%points(q_ix)
                     dlnks(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
             else
                     ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - CP%curv)
                     dlnks(q_ix) = CTrans%q%dpoints(q_ix)*CTrans%q%points(q_ix)/ks(q_ix)**2
             end if

             pows(q_ix) =  ScalarPower(ks(q_ix) ,in)
             if (global_error_flag/=0) return

          end do

         do j=1,CTrans%ls%l0
          do j2=1,CTrans%ls%l0

        !Integrate dk/k Delta_l_q**2 * Power(k)

          do q_ix = 1, CTrans%q%npoints 

             if (.not.(CP%closed.and.nint(CTrans%q%points(q_ix)*CP%r)<= CTrans%ls%l(j))) then 
               !cut off at nu = l + 1
             dlnk = dlnks(q_ix)
             apowers = pows(q_ix)

             iCl_scalar2(j,j2,C_Temp:C_E,in) = iCl_scalar2(j,j2,C_Temp:C_E,in) +  &
                          apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)*CTrans%Delta_p_l_k(1:2,j2,q_ix)*dlnk
             iCl_scalar2(j,j2,C_Cross,in) = iCl_scalar2(j,j2,C_Cross,in) + &
                          apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j2,q_ix)*dlnk
    
             end if

           end do

!Output l(l+1)C_l/OutputDenominator

           !ctnorm = (CTrans%ls%l+2)!/(CTrans%ls%l-2)! - beware of int overflow
            ctnorm=(CTrans%ls%l(j)*CTrans%ls%l(j)-1)*real((CTrans%ls%l(j)+2)*CTrans%ls%l(j),dl)
            ctnorm=sqrt(ctnorm*(CTrans%ls%l(j2)*CTrans%ls%l(j2)-1)*real((CTrans%ls%l(j2)+2)*CTrans%ls%l(j2),dl))
            
            dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*fourpi  
            dbletmp=sqrt(dbletmp*(CTrans%ls%l(j2)*(CTrans%ls%l(j2)+1))/OutputDenominator*fourpi  )
                 
            iCl_scalar2(j,j2,C_Temp,in)  =  iCl_scalar2(j,j2,C_Temp,in)*dbletmp
            iCl_scalar2(j,j2,C_E,in)     =  iCl_scalar2(j,j2,C_E,in)*dbletmp*ctnorm
            iCl_scalar2(j,j2,C_Cross,in) =  iCl_scalar2(j,j2,C_Cross,in)*dbletmp*sqrt(ctnorm)
    
           end do

          end do

        end do
        
        call CreateTxtFile('z:\cl2.dat',1)
        do j=1,CTrans%ls%l0
        do j2=1,CTrans%ls%l0
         write (1,*) CTrans%ls%l(j),CTrans%ls%l(j2),iCl_scalar2(j,j2,1,1)*7.4311e12
        end do
        end do  
        close(1)
        call CreateTxtFile('cl1l2.dat',1)
        do j=1,999
         write (1,'(999E15.5)') iCl_scalar2(j,1:999,1,1)*7.4311e12
        end do          
        stop

        end subroutine CalcScalCls2

     
        subroutine CalcTensCls(CTrans, GetInitPowers)
        implicit none
        Type(ClTransferData) :: CTrans
        external GetInitPowers       
        integer in,j, q_ix
        real(dl) nu
        real(dl) apowert,  measure
        real(dl) ctnorm,dbletmp
        real(dl) pows(CTrans%q%npoints)
        real(dl)  ks(CTrans%q%npoints),measures(CTrans%q%npoints)

        !For tensors we want Integral dnu/nu (nu^2-3)/(nu^2-1) Delta_l_k^2 P(k) for CP%closed

        do in=1,CP%InitPower%nn
        
         do q_ix = 1, CTrans%q%npoints 

             if (CP%flat) then
                     ks(q_ix) = CTrans%q%points(q_ix)
                     measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
             else
                     nu = CTrans%q%points(q_ix)*CP%r
                     ks(q_ix) = sqrt(CTrans%q%points(q_ix)**2 - 3*CP%curv)
                     measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)*(nu**2-3*CP%Ksign)/(nu**2-CP%Ksign)
             end if

          end do

         call GetInitPowers(pows,ks,CTrans%q%npoints,in)

        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,measure,apowert,ctnorm,dbletmp)
         do j=1,CTrans%ls%l0

          do q_ix = 1, CTrans%q%npoints 

                 if (.not.(CP%closed.and. nint(CTrans%q%points(q_ix)*CP%r)<=CTrans%ls%l(j))) then
                     !cut off at nu = l+1

                 apowert = pows(q_ix)
                 measure = measures(q_ix)
               
                 iCl_tensor(j,CT_Temp:CT_B,in) = iCl_tensor(j,CT_Temp:CT_B,in) + &
                      apowert*CTrans%Delta_p_l_k(CT_Temp:CT_B,j,q_ix)**2*measure
                 
                 iCl_tensor(j,CT_cross, in ) = iCl_tensor(j,CT_cross, in ) &
                      +apowert*CTrans%Delta_p_l_k(CT_Temp,j,q_ix)*CTrans%Delta_p_l_k(CT_E,j,q_ix)*measure
                 end if 
           end do

            ctnorm=(CTrans%ls%l(j)*CTrans%ls%l(j)-1)*real((CTrans%ls%l(j)+2)*CTrans%ls%l(j),dl)
            dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*pi/4
            iCl_tensor(j, CT_Temp, in)   = iCl_tensor(j, CT_Temp, in)*dbletmp*ctnorm
            if (CTrans%ls%l(j)==1) dbletmp=0
            iCl_tensor(j, CT_E:CT_B, in) = iCl_tensor(j, CT_E:CT_B, in)*dbletmp
            iCl_tensor(j, CT_Cross, in)  = iCl_tensor(j, CT_Cross, in)*dbletmp*sqrt(ctnorm)

         end do

         end do
    
        end subroutine CalcTensCls


        subroutine CalcVecCls(CTrans, GetInitPowers)
        implicit none
        Type(ClTransferData) :: CTrans
        external GetInitPowers       
        integer in,j, q_ix
        real(dl) power,  measure
        real(dl) ctnorm,lfac,dbletmp
        real(dl) pows(CTrans%q%npoints)
        real(dl)  ks(CTrans%q%npoints),measures(CTrans%q%npoints)

        do in=1,CP%InitPower%nn
        
         do q_ix = 1, CTrans%q%npoints 

               ks(q_ix) = CTrans%q%points(q_ix)
               measures(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)

          end do

         call GetInitPowers(pows,ks,CTrans%q%npoints,in)

        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,measure,power,ctnorm,dbletmp,lfac)
         do j=1,CTrans%ls%l0

          do q_ix = 1, CTrans%q%npoints 

             if (.not.(CP%closed.and. nint(CTrans%q%points(q_ix)*CP%r)<=CTrans%ls%l(j))) then
                     !cut off at nu = l+1

                 power = pows(q_ix)
                 measure = measures(q_ix)
               
                 iCl_vector(j,CT_Temp:CT_B,in) = iCl_vector(j,CT_Temp:CT_B,in) + &
                      power*CTrans%Delta_p_l_k(CT_Temp:CT_B,j,q_ix)**2*measure
                 
                 iCl_vector(j,CT_cross, in ) = iCl_vector(j,CT_cross, in ) &
                      +power*CTrans%Delta_p_l_k(CT_Temp,j,q_ix)*CTrans%Delta_p_l_k(CT_E,j,q_ix)*measure
             end if 
          end do

            ctnorm=CTrans%ls%l(j)*(CTrans%ls%l(j)+1)
            dbletmp=(CTrans%ls%l(j)*(CTrans%ls%l(j)+1))/OutputDenominator*pi/8
            iCl_vector(j, CT_Temp, in)   = iCl_vector(j, CT_Temp, in)*dbletmp*ctnorm
            lfac = (CTrans%ls%l(j) + 2)*(CTrans%ls%l(j) - 1)
            iCl_vector(j, CT_E:CT_B, in) = iCl_vector(j, CT_E:CT_B, in)*dbletmp*lfac
            iCl_vector(j, CT_Cross, in)  = iCl_vector(j, CT_Cross, in)*dbletmp*sqrt(lfac*ctnorm)

         end do

         end do
    
        end subroutine CalcVecCls


    subroutine InterpolateCls(CTransS,CTransT,CTransV)
      implicit none
       Type(ClTransferData) :: CTransS, CTransT, CTransV
      integer in,i
  
!Note using log interpolation is worse]

      !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,in), SHARED(CTransS,CTransT),IF(CP%InitPower%nn > 1)
      do in=1,CP%InitPower%nn
         if (CP%WantScalars) then
           do i = C_Temp, C_last
              call InterpolateClArrTemplated(CTransS%ls,iCl_scalar(1,i,in),Cl_scalar(lmin, in, i), &
                 CTransS%ls%l0,i)
           end do
         end if
      
         if (CP%WantVectors) then
           do i = C_Temp, CT_cross
            call InterpolateClArr(CTransV%ls,iCl_vector(1,i,in),Cl_vector(lmin, in, i),CTransV%ls%l0)
           end do
         end if
         
         if (CP%WantTensors) then
           do i = CT_Temp, CT_Cross
             call InterpolateClArr(CTransT%ls,iCl_tensor(1,i,in),Cl_tensor(lmin, in, i), &
               CTransT%ls%l0)
           end do
         end if
      end do
      !$OMP END PARALLEL DO
    end subroutine InterpolateCls


end module CAMBmain

    
     

    
