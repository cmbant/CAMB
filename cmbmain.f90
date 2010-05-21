!     This this is the main CAMB program module.
!
!     Code for Anisotropies in the Microwave Background
!     by Antony lewis (http://cosmologist.info) and Anthony Challinor
!     See readme.html for documentation. 
!     This version August 2006

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
!     of a Friedmann-Robertson-Walker universe with a supplied system of gauge-dependent
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
      use TimeSteps
      use Transfer
      use SpherBessels
      use lvalues
      use MassiveNu
      use InitialPower

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
      
      ! values of q to evolve the propagation equations to compute the sources
      integer num_q_evolve
      real(dl), dimension(:), allocatable :: q_evolve !array of q=nu/r (=k if K=0) values for sources

      real(dl),parameter :: qmin0=0.1_dl

      real(dl) :: dtaurec_q
    
!     qmax - CP%Max_eta_k/CP%tau0, qmin = qmin0/CP%tau0 for CP%flat case

      real(dl) qmin, qmax 

      real(dl) max_etak_tensor , max_etak_vector, max_etak_scalar
 !     Will only be calculated if k*tau < max_etak_xx

      integer maximum_l !Max value of l to compute
      real(dl) :: maximum_qeta = 3000._dl

      Type(ClTransferData), pointer :: ThisCT
                    
      public cmbmain, ClTransferToCl 

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

      if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*) actual-timeprev,' Timing for InitVars'
         write (*,*) 'r = ',real(CP%r),' scale = ',real(scale), 'age = ', real(CP%tau0)  
      end if 

       if (.not. CP%OnlyTransfers)  call InitializePowers(CP%InitPower,CP%curv)

!     Calculation of the CMB sources.


      if (CP%WantCls) call SetkValuesForSources 

      if (CP%WantTransfer) call InitTransfer
 

!      ***note that !$ is the prefix for conditional multi-processor compilation***
      !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)
        
    
      if (CP%WantCls) then

         if (DebugMsgs .and. Feedbacklevel > 0) write(*,*) 'Set ',num_q_evolve,' source k values'
         
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
         do q_ix= 1,num_q_evolve
             call DoSourcek(EV,q_ix)
         end do
         !$OMP END PARAllEl DO
        
         if (DebugMsgs .and. Feedbacklevel > 0) then
            timeprev=actual
            actual=GetTestTime()
            write(*,*) actual-timeprev,' Timing for source calculation'
         end if

      endif !WantCls
   
!     If transfer functions are requested, set remaining k values and output
      if (CP%WantTransfer) then
        call TransferOut
         if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*) actual-timeprev,' Timing for transfer k values'
         end if  
      end if

       if (CP%WantTransfer .and. CP%WantCls .and. CP%DoLensing &
            .and. CP%NonLinear==NonLinear_Lens) then
          
          call NonLinearLensing
          if (DebugMsgs .and. Feedbacklevel > 0) then
             timeprev=actual
             actual=GetTestTime()
             write(*,*) actual-timeprev,' Timing for NonLinear'
          end if

       end if

       if (CP%WantTransfer .and. .not. CP%OnlyTransfers) &
          call Transfer_Get_sigma8(MT,8._dl) 
           !Can call with other arguments if need different size
 
!     if CMB calculations are requested, calculate the Cl by
!     integrating the sources over time and over k.

      if (CP%WantCls) then

         call InitSourceInterpolation   
       
         ExactClosedSum = CP%curv > 5e-9_dl .or. scale < 0.93_dl


         call SetkValuesForInt

         if (DebugMsgs .and. Feedbacklevel > 0) write(*,*) 'Set ',ThisCT%num_q_int,' integration k values'

    
      !Begin k-loop and integrate Sources*Bessels over time
      
      !$OMP PARAllEl DO DEFAUlT(SHARED),SHARED(nstep), SCHEDUlE(STATIC,4) 
          do q_ix=1,ThisCT%num_q_int
            call SourceToTransfers(q_ix)
         end do !q loop
  
       !$OMP END PARAllEl DO 
  
        if (DebugMsgs .and. Feedbacklevel > 0) then
         timeprev=actual
         actual=GetTestTime()
         write(*,*)actual-timeprev,' Timing For Integration'
        end if
  

        call FreeSourceMem

        !Final calculations for CMB output unless want the Cl transfer functions only.

        if (.not. CP%OnlyTransfers) call ClTransferToCl(CTransScal,CTransTens, CTransVec)

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

       if (CP%WantScalars) then
           lSamp = CTransS%ls
           allocate(iCl_Scalar(CTransS%ls%l0,C_Temp:C_last,CP%InitPower%nn))
           iCl_scalar = 0
           call CalcScalCls(CTransS,GetInitPowerArrayScal)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcScalCls'
       end if    

       if (CP%WantVectors) then
           allocate(iCl_vector(CTransV%ls%l0,C_Temp:CT_Cross,CP%InitPower%nn))
           iCl_vector = 0
           call CalcVecCls(CTransV,GetInitPowerArrayVec)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcVecCls'
       end if    


       if (CP%WantTensors) then
           allocate(iCl_Tensor(CTransT%ls%l0,CT_Temp:CT_Cross,CP%InitPower%nn))
           iCl_tensor = 0
           call CalcTensCls(CTransT,GetInitPowerArrayTens)
           if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcTensCls'
       end if

       call Init_Cls
 !     Calculating Cls for every l.
       call InterpolateCls(CTransS,CTransT, CTransV)
   
       if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'InterplolateCls'

       if (CP%WantScalars) deallocate(iCl_scalar)
       if (CP%WantVectors) deallocate(iCl_vector)
       if (CP%WantTensors) deallocate(iCl_tensor)
 
       if (CP%OutputNormalization >=2) call NormalizeClsAtl(CP%OutputNormalization)
       !Normalize to C_l=1 at l=OutputNormalization 
  
     
     end subroutine ClTransferToCl
    

     subroutine SourceToTransfers(q_ix)
      integer q_ix
      type(IntegrationVars) :: IV
    
          allocate(IV%Source_q(nstep,SourceNum))
          if (.not.CP%flat) allocate(IV%ddSource_q(nstep,SourceNum))

            call IntegrationVars_init(IV)

            IV%q_ix = q_ix
            IV%q =ThisCT%q_int(q_ix)
            IV%dq= ThisCT%dq_int(q_ix)

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
        !Log spacing for last few osciallations
        !large log spacing for small scales

           boost = AccuracyBoost
           if (CP%Transfer%high_precision) boost = boost*1.5

           q_switch_lowk1 = 0.7/taurst
           dlog_lowk1=2*boost

           q_switch_lowk = 8/taurst
           dlog_lowk=8*boost

           q_switch_osc = min(CP%Transfer%kmax,30/taurst)
           d_osc= 200*boost

           q_switch_highk = min(CP%Transfer%kmax,60/taurst)
           dlog_osc = 17*boost

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
               MT%num_q_trans = num_q_evolve 
              else
               MT%num_q_trans = num_q_evolve + (ntodo - first_i+1) 
              end if
              call Transfer_Allocate(MT)

              MT%q_trans(1:num_q_evolve) = q_evolve(1:num_q_evolve)
              if (MT%num_q_trans > num_q_evolve) then
               MT%q_trans(num_q_evolve+1:MT%num_q_trans) = q_transfer(first_i:ntodo)
          end if

         else
             num_q_evolve = 0
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

            EV%q=q_evolve(q_ix) 
 
            EV%q2=EV%q**2

            EV%q_ix = q_ix
            EV%TransferOnly=.false.
      
            taustart = GetTauStart(EV%q)
   
            call GetNumEqns(EV)

            if (CP%WantScalars) call CalcScalarSources(EV,taustart)
            if (CP%WantVectors) call CalcVectorSources(EV,taustart)
            if (CP%WantTensors) call CalcTensorSources(EV,taustart)

      end subroutine DoSourcek

       subroutine GetSourceMem
     
        if (CP%WantScalars) then
           if (CP%Dolensing) then
            SourceNum=3
            C_last = C_PhiTemp
         else
            SourceNum=2
            C_last = C_Cross
           end if
        else
           SourceNum=3 
        end if
       
        allocate(Src(num_q_evolve,SourceNum,nstep))
        Src=0
        allocate(ddSrc(num_q_evolve,SourceNum,nstep))
    
       end subroutine GetSourceMem


       subroutine FreeSourceMem
         
        deallocate(Src, ddSrc)
        deallocate(q_evolve)

       end subroutine FreeSourceMem


!  initial variables, number of steps, etc.
      subroutine InitVars
      use ThermoData
      use precision
      use ModelParams
      
      implicit none
      real(dl) taumin, maxq
      integer itf

 ! Maximum and minimum k-values.      
      if (CP%flat) then
      qmax=maximum_qeta/CP%tau0
      qmin=qmin0/CP%tau0/AccuracyBoost
      else              
        qmax=maximum_qeta/CP%r/CP%chi0
        qmin=qmin0/CP%r/CP%chi0/AccuracyBoost
      end if

!     Timesteps during recombination (tentative, the actual
!     timestep is the minimum between this value and taurst/40,
!     where taurst is the time when recombination starts - see inithermo

      dtaurec_q=4/qmax/AccuracyBoost  
      if (.not. CP%flat) dtaurec_q=dtaurec_q/6
      !AL:Changed Dec 2003, dtaurec feeds back into the non-flat integration via the step size
      dtaurec = dtaurec_q
      !dtau rec may be changed by inithermo

      max_etak_tensor = AccuracyBoost*maximum_qeta /10  
      max_etak_scalar = AccuracyBoost*max(1700._dl,maximum_qeta) /20 
      if (maximum_qeta <3500 .and. AccuracyBoost < 2) max_etak_scalar = max_etak_scalar * 1.5
        !tweak to get large scales right
      max_etak_vector = max_etak_scalar
      call Reionization_Init
 
      if (CP%WantCls) then     
         maxq = qmax
         if (CP%WantTransfer) maxq=max(qmax,CP%Transfer%kmax)
      else
         maxq=CP%Transfer%kmax
      end if

      taumin=GetTauStart(maxq)
   
!     Initialize baryon temperature and ionization fractions vs. time.
!     This subroutine also fixes the timesteps where the sources are
!     saved in order to do the integration. So nstep is set here.
      !These routines in ThermoData (modules.f90)
  
      call inithermo(taumin,CP%tau0)
    
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
      integer nk1,nk2, q_ix

         num_q_evolve=0

!     set k values for which the sources for the anisotropy and
!     polarization will be calculated. For low values of k we
!     use a logarithmic spacing. closed case dealt with by SetClosedkValues
         if (CP%WantScalars .and. CP%Reionization .and. CP%AccuratePolarization) then
            dlnk0=2._dl/10/AccuracyBoost
            !Need this to get accurate low l polarization
         else
            dlnk0=5._dl/10/AccuracyBoost
            if (CP%closed) dlnk0=dlnk0/2
         end if

         if (CP%AccurateReionization) dlnk0 = dlnk0/2

         dkn1=0.6_dl/taurst/AccuracyBoost 
         dkn2=1.1_dl/taurst/AccuracyBoost 
          if (CP%WantTensors .or. CP%WantVectors) then
              dkn1=dkn1  *0.8_dl
              dlnk0=dlnk0/2 !*0.3_dl
              dkn2=dkn2*0.7_dl
          end if
         nk1=int(log(dkn1/qmin/dlnk0)/dlnk0)+1

         q_switch = 2*6.3/taurst
           !Want linear spacing for wavenumbers which come inside horizon
           !Could use sound horizon, but for tensors that is not relevant
         if (qmax > q_switch) then
            nk2=int((q_switch-qmin*exp((nk1-1)*dlnk0)) /dkn1)+nk1+1
            num_q_evolve=int((qmax-q_switch)/dkn2)+nk2+1
         else
             num_q_evolve=int((qmax-qmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1
            nk2=num_q_evolve+1
         end if

         allocate(q_evolve(num_q_evolve))

         do q_ix=1,num_q_evolve
            if (q_ix <= nk1) then
               q_evolve(q_ix)=qmin*exp((q_ix-1)*dlnk0)
            else
               if (q_ix > nk2) then
                  q_evolve(q_ix)=q_evolve(nk2)+(q_ix-nk2)*dkn2
               else
                  q_evolve(q_ix)=q_evolve(nk1)+(q_ix-nk1)*dkn1
               end if
            end if

         end do

         if (CP%closed) call SetClosedkValuesFromArr(q_evolve,num_q_evolve)

      end subroutine SetkValuesForSources


      subroutine SetClosedkValuesFromArr(arr,size)
      real(dl) arr(*)
      integer i,nu,lastnu,size,nmax
       !nu = 3,4,5... in CP%closed case, so set nearest integers from arr array
                
       lastnu=3
       nmax=1
      
       do i=2,size
         nu=nint(arr(i)*CP%r)
         if (nu > lastnu) then
          nmax=nmax+1
          lastnu=nu         
          arr(nmax)=nu/CP%r 
         end if
       
       end do  
       arr(1)=3/CP%r
     
        size=nmax

      end subroutine SetClosedkValuesFromArr



      subroutine CalcScalarSources(EV,taustart)
      use Transfer
      implicit none
      type(EvolutionVars) EV
      real(dl) tau,tol1,tauend, taustart
      integer j,ind,itf
      real(dl) c(24),w(EV%nvar,9), y(EV%nvar), sources(SourceNum)

!!
! EV%q=0.3
! EV%q2=EV%q**2

         w=0
         y=0
         call initial(EV,y, taustart)

       
         tau=taustart
         ind=1

!!Example code for plotting out variable evolution
!if (.false.) then
   !        tol1=tol/exp(AccuracyBoost-1)
  !      call GaugeInterface_EvolveScal(EV,tau,y,CP%tau0,tol1,ind,c,w)
       
     !!   do j=1,6000      
      !    tauend = taustart * exp(j/6000._dl*log(CP%tau0/taustart))
      !    call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
         
       !   write (*,'(8E15.5)') y(1:7), grhoc*y(1)/(grhog)                            
       !  end do
 !     stop
!end if


!     Begin timestep loop.
!
!     d contains the sources for the anisotropy and dp
!     for the polarization. t means tensor.

           itf=1
           tol1=tol/exp(AccuracyBoost-1)
           if (CP%WantTransfer .and. CP%Transfer%high_precision) tol1=tol1/100

           do j=2,nstep               
             tauend=atau0(j)  
             if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > taurend) &
                  .and. .not. CP%Dolensing .and. &
                  (.not.CP%WantTransfer.or.tau > tautf(CP%Transfer%num_redshifts))) then
      
              Src(EV%q_ix,1:SourceNum,j)=0
 
             else
             
             !Integrate over time, calulate end point derivs and calc output
             call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
 
             call output(EV,y, EV%ScalEqsToPropagate,j,tau,sources)
             Src(EV%q_ix,1:SourceNum,j)=sources
 
             if (CP%Num_Nu_Massive > 0 .and.(EV%NuMethod==Nu_trunc).and..not.EV%MassiveNuApprox.and. &
                  .not.CP%Transfer%high_precision.and. &
                 ((EV%q<0.1_dl .and.EV%w_nu < 0.015/AccuracyBoost/lAccuracyBoost).or.&
                  (EV%w_nu < 0.008/AccuracyBoost/lAccuracyBoost))) then 
     
                !Neutrinos no longer highly relativistic so make approximation               
                if (.not.CP%DoLensing.and..not. CP%WantTransfer &
                       .or. EV%w_nu < 1e-8/EV%q**2/AccuracyBoost/lAccuracyBoost)  &
                           call SwitchToMassiveNuApprox(EV, y)
              end if
             
!     Calculation of transfer functions.
101          if (CP%WantTransfer.and.itf <= CP%Transfer%num_redshifts) then
                  if (j < nstep.and.tauend < tautf(itf) .and.atau0(j+1) > tautf(itf)) then
                          
                     call GaugeInterface_EvolveScal(EV,tau,y,tautf(itf),tol1,ind,c,w)

                  endif
!     output transfer functions for this k-value.
                      
                  if (abs(tau-tautf(itf)) < 1.e-5_dl) then
                           call outtransf(EV,y, MT%TransferData(:,EV%q_ix,itf))
                         
                           itf=itf+1
                           if (j < nstep.and.itf <= CP%Transfer%num_redshifts.and. &
                                atau0(j+1) > tautf(itf)) goto 101
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
               do j=2,nstep
                  tauend=atau0(j)
                  if (EV%q*tauend > max_etak_tensor) then
                     Src(EV%q_ix,1:SourceNum,j) = 0
                   else
      
                     if (CP%flat) then
                      call dverk(EV,EV%nvart,fderivst,tau,yt,tauend,tol1,ind,c,EV%nvart,wt) !tauend
                     else 
                      call dverk(EV,EV%nvart, derivst,tau,yt,tauend,tol1,ind,c,EV%nvart,wt)
                     end if
 
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
               do j=2,nstep
                  tauend=atau0(j)

                  if ( EV%q*tauend > max_etak_vector) then
                     Src(EV%q_ix,1:SourceNum,j) = 0
                   else
      
                      call dverk(EV,EV%nvarv,fderivsv,tau,yv,tauend,tol1,ind,c,EV%nvarv,wt) !tauend
                 
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
         write(*,*) MT%num_q_trans-num_q_evolve, 'transfer k values'

      !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
      !$OMP & PRIVATE(EV, tau, q_ix) 

!     loop over wavenumbers.
         do q_ix=num_q_evolve+1,MT%num_q_trans
      
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
         
       do i=1,CP%Transfer%num_redshifts
          call GaugeInterface_EvolveScal(EV,tau,y,tautf(i),atol,ind,c,w)
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
        do while(atau0(first_step) < tautf(1))
           first_step = first_step + 1 
        end do

        do ik=1, num_q_evolve
         if (q_evolve(ik)/(CP%H0/100) >  Min_kh_nonlinear) then

            !Interpolate non-linear scaling in conformal time
            do i = 1, CP%Transfer%num_redshifts
                scaling(i) = CAMB_Pk%nonlin_ratio(ik,i)
            end do
            if (all(abs(scaling-1) < 5e-4)) cycle
            call spline(tautf,scaling(1),CP%Transfer%num_redshifts,&
                                 spl_large,spl_large,ddScaling(1))
       
            tf_lo=1
            tf_hi=tf_lo+1

            do i=first_step,nstep-1
  
              tau = atau0(i)
              
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
      !$OMP PARAllEl DO DEFAUlT(SHARED), SCHEDUlE(STATIC), PRIVATE(i,j) , SHARED(num_q_evolve)
        do  i=1,nstep
          do j=1, SourceNum
               call spline(q_evolve,Src(1,j,i),num_q_evolve,spl_large,spl_large,ddSrc(1,j,i))
          end do
        end do
       !$OMP END PARAllEl DO 
       end subroutine InitSourceInterpolation

      subroutine SetkValuesForInt
      implicit none
       
      integer k,no,no1,no2
      real(dl) dk,dk0,dlnk1, dk2, max_k_dk
      integer lognum
      real(dl)  IntSampleBoost

      IntSampleBoost=AccuracyBoost

!     Fixing the # of k for the integration. 
      
      if (CP%closed.and.ExactClosedSum) then
       
        ThisCT%num_q_int = nint(qmax*CP%r)-2
        call Init_ClTransfer(ThisCT)
 
        do k=1,ThisCT%num_q_int
               ThisCT%q_int(k) = (k+2)/CP%r
               ThisCT%dq_int(k) = 1/CP%r  
        end do 

       else
      !Split up into logarithmically spaced intervals from qmin up to k=lognum*dk0
      !then no-lognum*dk0 linearly spaced at dk0 up to no*dk0
      !then at dk up to qmax
         if (CP%closed) then
            lognum=nint(20*IntSampleBoost)
            dlnk1=1._dl/lognum 
            no=nint(1300*IntSampleBoost)
            dk=1.2/CP%r/CP%chi0/IntSampleBoost
            dk0=0.4_dl/CP%r/CP%chi0/IntSampleBoost
         else 
            lognum=nint(10*IntSampleBoost)
            dlnk1=1._dl/lognum  
            no=nint(600*IntSampleBoost)
            dk0=1.8_dl/CP%r/CP%chi0/IntSampleBoost
            dk=3._dl/CP%r/CP%chi0/IntSampleBoost
         end if

         no1=int(log(lognum*dk0/qmin)/dlnk1)+1

         dk2 = 0.04/IntSampleBoost

         if (qmax > (no*dk0)) then
            max_k_dk = max(3000, 2*maximum_l)/CP%tau0
            if (qmax > max_k_dk) then
               no2 = int((max_k_dk-no*dk0)/dk)+no
               ThisCT%num_q_int=int((qmax-max_k_dk)/dk2)+no2
            else
               ThisCT%num_q_int=int((qmax-no*dk0)/dk)+no
               no2 = ThisCT%num_q_int+1
            end if
         else
            no=int((qmax-lognum*dk0)/dk0)+no1
            ThisCT%num_q_int=no
            no2 = ThisCT%num_q_int+1
         end if

         call Init_ClTransfer(ThisCT)

         do  k=1,ThisCT%num_q_int 
           if (k <=no2) then
            if (k <= no) then
               if (k <= no1) then
                  ThisCT%q_int(k)=lognum*dk0*exp(-(no1-k)*dlnk1)
                  if (k ==no1) then
                   ThisCT%dq_int(no1)=lognum*dk0*(1._dl-exp(-dlnk1/2)) + dk0/2
                  else
                   ThisCT%dq_int(k)=ThisCT%q_int(k)*2*sinh(dlnk1/2)
                  end if
               else
                  ThisCT%q_int(k)=ThisCT%q_int(no1)+(k-no1)*dk0
                  ThisCT%dq_int(k)=dk0
               end if
            else
               ThisCT%q_int(k)=ThisCT%q_int(no)+(k-no)*dk
               ThisCT%dq_int(k)=dk
            end if
           else
              !This allows inclusion of high k modes for computing BB lensed spectrum accurately
              !without taking ages to compute.
               ThisCT%q_int(k)=ThisCT%q_int(no2)+(k-no2)*dk2
               ThisCT%dq_int(k)=dk2
           end if
         end do
        
        ThisCT%dq_int(no)=(dk+dk0)/2
        if (no2<ThisCT%num_q_int) ThisCT%dq_int(no2) = (dk+dk2)/2

 
       if (CP%closed) then
        call SetClosedkValuesFromArr(ThisCT%q_int,ThisCT%num_q_int)
        do k=2,ThisCT%num_q_int-1
           ThisCT%dq_int(k)=(ThisCT%q_int(k+1)-ThisCT%q_int(k-1))/2
        end do
        ThisCT%dq_int(1)=1/CP%r
        ThisCT%dq_int(ThisCT%num_q_int)=ThisCT%dq_int(ThisCT%num_q_int-1)        
       end if

       end if !ExactClosedSum

         
      end subroutine setkValuesForInt

      subroutine InterpolateSources(IV)
      implicit none
      integer i,khi,klo, step
      real(dl) xf,b0,ho,a0,ho2o6,a03,b03
      type(IntegrationVars) IV


!     finding position of k in table q_evolve to do the interpolation.

            klo=1
            !This is a bit inefficient, but thread safe
            do while ((IV%q > q_evolve(klo+1)).and.(klo < (num_q_evolve-1))) 
               klo=klo+1
            end do
            khi=klo+1

            ho=q_evolve(khi)-q_evolve(klo)
            a0=(q_evolve(khi)-IV%q)/ho
            b0=(IV%q-q_evolve(klo))/ho           
            ho2o6 = ho**2/6
            a03=(a0**3-a0)
            b03=(b0**3-b0)
            IV%SourceSteps = 0

!     Interpolating the source as a function of time for the present
!     wavelength.
             step=2
               do i=2, nstep
                  xf=IV%q*(CP%tau0-atau0(i))
                  if (CP%WantTensors) then
                   if (IV%q*atau0(i) < max_etak_tensor.and. xf > 1.e-8_dl) then
                     step=i
                     IV%Source_q(i,1:SourceNum) =a0*Src(klo,1:SourceNum,i)+&
                          b0*Src(khi,1:SourceNum,i)+(a03 *ddSrc(klo,1:SourceNum,i)+ &
                          b03*ddSrc(khi,1:SourceNum,i)) *ho2o6
                    else
                     IV%Source_q(i,1:SourceNum) = 0
                   end if
                  end if
                  if (CP%WantVectors) then
                   if (IV%q*atau0(i) < max_etak_vector.and. xf > 1.e-8_dl) then
                     step=i
                     IV%Source_q(i,1:SourceNum) =a0*Src(klo,1:SourceNum,i)+&
                          b0*Src(khi,1:SourceNum,i)+(a03 *ddSrc(klo,1:SourceNum,i)+ &
                          b03*ddSrc(khi,1:SourceNum,i)) *ho2o6
                    else
                     IV%Source_q(i,1:SourceNum) = 0
                   end if
                  end if

                  if (CP%WantScalars) then
                     if ((CP%Dolensing .or. IV%q*atau0(i) < max_etak_scalar) .and. xf > 1.e-8_dl) then
                        step=i
                        IV%Source_q(i,1:SourceNum)=a0*Src(klo,1:SourceNum,i)+ &
                         b0*Src(khi,1:SourceNum,i) + (a03*ddSrc(klo,1:SourceNum,i) &
                         +b03*ddSrc(khi,1:SourceNum,i))*ho2o6
  
                     else
                       IV%Source_q(i,1:SourceNum) = 0 
                     end if
                  end if
               end do
               step=max(step,n1)
               IV%SourceSteps = step
  

          if (.not.CP%flat) then
             do i=1, SourceNum
                call spline(atau0,IV%Source_q(1,i),nstep,spl_large,spl_large,IV%ddSource_q(1,i))
             end do
          end if
           
           IV%SourceSteps = IV%SourceSteps*1
           !This is a fix for a compiler bug on Seaborg 

      end subroutine

      
      subroutine IntegrationVars_Init(IV)
      type(IntegrationVars), intent(INOUT) :: IV
        
        IV%Source_q(1,1:SourceNum)=0
        IV%Source_q(nstep,1:SourceNum) = 0
        IV%Source_q(nstep-1,1:SourceNum) = 0
 
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
           llmax=min(llmax,nint(nu)-1)  !nu >= lSamp%l+1
          end if       

         else
          llmax=nint(nu*CP%chi0)
          
          if (llmax<15) then
           llmax=15
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
 
        if (CP%AccurateBB .or. CP%flat) then
         UseLimber = l > 700*AccuracyBoost .and. k > 0.05
        else
         !This is accurate at percent level only (good enough here)
         UseLimber = l > 300*AccuracyBoost .or. k>0.05
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
        integer bes_ix,n, bes_index(IV%SourceSteps)
     
!     Find the position in the xx table for the x correponding to each
!     timestep

         do j=1,IV%SourceSteps !Precompute arrays for this k
            xf=abs(IV%q*(CP%tau0-atau0(j)))
            bes_index(j)=GetBesselIndex(xf)
          !Precomputed values for the interpolation
            bes_ix= bes_index(j)
            fac(j)=xx(bes_ix+1)-xx(bes_ix)
            aa(j)=(xx(bes_ix+1)-xf)/fac(j)
            fac(j)=fac(j)**2*aa(j)/6
         end do

         do j=1,lSamp%l0 
             if (lSamp%l(j) > llmax) return
             xlim=xlimfrac*lSamp%l(j)
             xlim=max(xlim,xlimmin)
             xlim=lSamp%l(j)-xlim
             xlmax1=80*lSamp%l(j)
             tmin=CP%tau0-xlmax1/IV%q
             tmax=CP%tau0-xlim/IV%q
             tmax=min(CP%tau0,tmax)
             tmin=max(atau0(2),tmin)
               
             if (tmax < atau0(2)) exit
             sums(1:SourceNum) = 0
 
            !As long as we sample the source well enough, it is sufficient to
            !interpolate the Bessel functions only

             if (SourceNum==2) then
              !This is the innermost loop, so we separate the no lensing scalar case to optimize it
                do n= TimeToTimeStep(tmin),min(IV%SourceSteps,TimeToTimeStep(tmax))

                 a2=aa(n)
                 bes_ix=bes_index(n) 

                 J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                        *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline
                 J_l = J_l*dtau2(n)
                 sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                 sums(2) = sums(2) + IV%Source_q(n,2)*J_l

                end do
              else 
                 DoInt = .not. CP%WantScalars .or. IV%q < max(850,lSamp%l(j))*3*AccuracyBoost/CP%tau0  
                 if (DoInt) then
                  do n= TimeToTimeStep(tmin),min(IV%SourceSteps,TimeToTimeStep(tmax))
                   !Full Bessel integration
                     a2=aa(n)
                     bes_ix=bes_index(n) 

                     J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                            *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline
                     J_l = J_l*dtau2(n)
              
                    !The unwrapped form is faster
     
                     sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                     sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                     sums(3) = sums(3) + IV%Source_q(n,3)*J_l
    
                  end do
                 end if
                 if (.not. DoInt .or. UseLimber(lsamp%l(j),IV%q) .and. CP%WantScalars) then
                  !Limber approximation for small scale lensing (better than poor version of above integral)
                  xf = CP%tau0-lSamp%l(j)/IV%q
                  n=TimeToTimeStep(xf)
                  xf= (xf-atau0(n))/(atau0(n+1)-atau0(n))                  
                  sums(3) = (IV%Source_q(n,3)*(1-xf) + xf*IV%Source_q(n+1,3))*sqrt(pi/2/lSamp%l(j))/IV%q 
                 end if

              end if
  
              ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) = ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) + sums(1:SourceNum)
 !             IV%Delta_l_q(1:SourceNum,j) = IV%Delta_l_q(1:SourceNum,j) + sums(1:SourceNum)
         
          end do

      
        end subroutine DoFlatIntegration
   

    
!non-flat source integration

      subroutine IntegrateSourcesBessels(IV,j,l,nu)  
      use USpherBessels
      type(IntegrationVars) IV
      logical DoInt   
      integer l,j, nstart,nDissipative,ntop,nbot,nrange,nnow
      real(dl) nu,ChiDissipative,ChiStart,tDissipative,y1,y2,y1dis,y2dis     
      real(dl) xf,x,chi, miny1
      real(dl) sums(SourceNum),out_arr(SourceNum)     
    
      !Calculate chi where for smaller chi it is dissipative
      x=sqrt(real(l*(l+1),dl))/nu
    
      ChiDissipative=invsinfunc(x) 
    
      ChiStart=ChiDissipative
      !Move down a bit to get smaller value (better accuracy integrating up from small values)
      if (nu<300) ChiStart = max(ChiDissipative-1._dl/nu,1d-6)   !max(ChiDissipative-1._dl/nu,1d-6)
 
        !Then get nearest source point with lower Chi...
      tDissipative=CP%tau0 - CP%r*ChiStart
      if (tDissipative<atau0(1)) then
         nDissipative=2
       else 
         nDissipative = TimeToTimeStep(tDissipative)+1
      endif 
      nDissipative=min(nDissipative,nstep-1)

      tDissipative = atau0(nDissipative)

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

         DoInt =  SourceNum/=3 .or. IV%q < max(850,l)*3*AccuracyBoost/(CP%chi0*CP%r) 
         if (DoInt) then
         if ((nstart < min(nstep-1,IV%SourceSteps)).and.(y1dis > miny1)) then
     
            y1=y1dis
            y2=y2dis
            nnow=nstart
            do nrange = 1,nr
               ntop = min(nreg(nrange+1),nstep-1)
               if (nnow < ntop) then
                  call DoRangeInt(IV,chi,ChiDissipative,nnow,ntop,dtaureg(nrange), &
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
            do nrange = nr,1,-1 
               nbot = nreg(nrange)         
               if (nnow >  nbot) then
                  call DoRangeInt(IV,chi,ChiDissipative,nnow,nbot,dtaureg(nrange), &
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
             nbot=TimeToTimeStep(xf)
             xf= (xf-atau0(nbot))/(atau0(nbot+1)-atau0(nbot))                  
             sums(3) = (IV%Source_q(nbot,3)*(1-xf) + xf*IV%Source_q(nbot+1,3))*&
                           sqrt(pi/2/l/sqrt(1-CP%Ksign*real(l**2)/nu**2))/IV%q 
         end if

         ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix)=ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix)+sums
           
      end if !Do Scalars
           
     if ((CP%WantTensors)) then !Do Tensors
         chi=ChiStart
        
      !Integrate chi down in dissipative region
      !DoRangeInt cuts off when ujl gets small
         miny1= 1.d-6/l/AccuracyBoost
         if ((nstart < nstep-1).and.(y1dis>miny1)) then
            y1=y1dis
            y2=y2dis
            nnow=nstart
            do nrange = 1,nr
               ntop = min(nreg(nrange+1),nstep-1)
               if (nnow < ntop) then
                 call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,ntop,dtaureg(nrange), &
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
           do nrange = nr,1,-1 
               nbot = nreg(nrange)         
               if (nnow >  nbot) then
                 call DoRangeIntTensor(IV,chi,ChiDissipative,nnow,nbot,dtaureg(nrange), &
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
      use USpherBessels
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
      else if (scalel>=1100) then
         num2=num1*1.5
      else if (scalel< 50) then
         num2=num1*0.8_dl
      else if (scalel < 170) then 
        num2=num1*1.5_dl
      else if (scalel < 1100) then 
         num2=num1*1.2_dl
      else
         num2=num1
      end if    
      !Dec 2003, since decrease dtaurec, can make this smaller
      if (dtau==dtaurec_q) then
       num2=num2/4
      end if
      if (num2*IntAccuracyBoost < dchisource .and. (.not. CP%DoLensing .or. UseLimber(l,IV%q)) & 
!       if ((num2*IntAccuracyBoost < dchisource ) & !Oscillating fast 
        .or. (nstart>IV%SourceSteps.and.nend>IV%SourceSteps)) then  
         out = 0
         y1=0._dl !So we know to calculate starting y1,y2 if there is next range
         y2=0._dl
         chi=(CP%tau0-atau0(nend))/CP%r
         return
        end if
     
      Startn=nstart
      if (nstart>IV%SourceSteps .and. nend < IV%SourceSteps) then
         chi=(CP%tau0-atau0(IV%SourceSteps))/CP%r
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
      aux2=(CP%tau0-atau0(nlowest))/dtau + nlowest
          
      nu2=nu*nu
      ap1=l*(l+1)
      sh=rofChi(chi)
           
      if (scalel < 120) then
         dchimax=0.4_dl*num1
      else if (scalel < 1400) then
         dchimax=0.25_dl*num1
      else
         dchimax=0.35_dl*num1 
      end if

      dchimax=dchimax/IntAccuracyBoost
    
      ujl=y1/sh
      sources = IV%Source_q(Startn,1:SourceNum)

      out = 0.5_dl*ujl*sources

      Interpolate = dchisource > dchimax
      if (Interpolate) then !split up smaller than source step size
         delchi=dchimax
         Deltachi=sgn*(atau0(Startn)-atau0(nend))/CP%r
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
      use USpherBessels
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
      
     
! atau0 is the array with the time where the sources are stored.
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
         chi=(CP%tau0-atau0(nend))/CP%r
         return
        end if
      if ((y2==0._dl).and.(y1==0._dl)) call USpherBesselWithDeriv(CP%closed,chi,l,nu,y1,y2)
   
      sgn=isgn

      nlowest=min(nstart,nend)
      aux1=1._dl*CP%r/dtau  !used to calculate nearest timestep quickly
      aux2=(CP%tau0-atau0(nlowest))/dtau + nlowest
     
     
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
         Deltachi=sgn*(atau0(nstart)-atau0(nend))/CP%r
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


        subroutine GetInitPowerArrayScal(pows,ks, numks,pix)
         integer, intent(in) :: numks, pix
         real(dl) pows(numks), ks(numks)
         integer i
      
         do i = 1, numks
            pows(i) =  ScalarPower(ks(i) ,pix)
         end do

        end subroutine GetInitPowerArrayScal

        subroutine GetInitPowerArrayVec(pows,ks, numks,pix)
         integer, intent(in) :: numks, pix
         real(dl) pows(numks), ks(numks)
         integer i
      
         do i = 1, numks
         !!change to vec...
            pows(i) =  ScalarPower(ks(i) ,pix)
         end do

        end subroutine GetInitPowerArrayVec


        subroutine GetInitPowerArrayTens(pows,ks, numks,pix)
         integer, intent(in) :: numks, pix
         real(dl) pows(numks), ks(numks)
         integer i
      
         do i = 1, numks
            pows(i) =  TensorPower(ks(i) ,pix)
         end do

        end subroutine GetInitPowerArrayTens


        subroutine CalcScalCls(CTrans, GetInitPowers)
        implicit none
        Type(ClTransferData) :: CTrans
        external GetInitPowers
        integer pix,j
        real(dl) apowers, pows(CTrans%num_q_int)
        integer q_ix
        real(dl)  ks(CTrans%num_q_int),dlnks(CTrans%num_q_int),dlnk
        real(dl) ctnorm,dbletmp

         do pix=1,CP%InitPower%nn

          do q_ix = 1, CTrans%num_q_int 

             if (CP%flat) then
                     ks(q_ix) = CTrans%q_int(q_ix)
                     dlnks(q_ix) = CTrans%dq_int(q_ix)/CTrans%q_int(q_ix)
             else
                     ks(q_ix) = sqrt(CTrans%q_int(q_ix)**2 - CP%curv)
                     dlnks(q_ix) = CTrans%dq_int(q_ix)*CTrans%q_int(q_ix)/ks(q_ix)**2
             end if
      
          end do

         call GetInitPowers(pows,ks,CTrans%num_q_int,pix)


        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,dlnk,apowers,ctnorm,dbletmp)

         do j=1,CTrans%ls%l0

        !Integrate dk/k Delta_l_q**2 * Power(k)

          do q_ix = 1, CTrans%num_q_int 

             if (.not.(CP%closed.and.nint(CTrans%q_int(q_ix)*CP%r)<=CTrans%ls%l(j))) then 
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
                            iCl_scalar(j,C_Phi,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)**2    
                     !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
                     !We put pix extra factors of %l here to improve interpolation in CTrans%ls%l
                     iCl_scalar(j,C_PhiTemp,pix)   =  &
                            iCl_scalar(j,C_PhiTemp,pix)*fourpi*real(CTrans%ls%l(j)**2,dl)*CTrans%ls%l(j)
                      !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
             end if

           end do
         !$OMP END PARAllEl DO

          end do

        end subroutine CalcScalCls

     
        subroutine CalcTensCls(CTrans, GetInitPowers)
        implicit none
        Type(ClTransferData) :: CTrans
        external GetInitPowers       
        integer in,j, q_ix
        real(dl) nu
        real(dl) apowert,  measure
        real(dl) ctnorm,dbletmp
        real(dl) pows(CTrans%num_q_int)
        real(dl)  ks(CTrans%num_q_int),measures(CTrans%num_q_int)

        !For tensors we want Integral dnu/nu (nu^2-3)/(nu^2-1) Delta_l_k^2 P(k) for CP%closed

        do in=1,CP%InitPower%nn
        
         do q_ix = 1, CTrans%num_q_int 

             if (CP%flat) then
                     ks(q_ix) = CTrans%q_int(q_ix)
                     measures(q_ix) = CTrans%dq_int(q_ix)/CTrans%q_int(q_ix)
             else
                     nu = CTrans%q_int(q_ix)*CP%r
                     ks(q_ix) = sqrt(CTrans%q_int(q_ix)**2 - 3*CP%curv)
                     measures(q_ix) = CTrans%dq_int(q_ix)/CTrans%q_int(q_ix)*(nu**2-3*CP%Ksign)/(nu**2-CP%Ksign)
             end if

          end do

         call GetInitPowers(pows,ks,CTrans%num_q_int,in)

        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,measure,apowert,ctnorm,dbletmp)
         do j=1,CTrans%ls%l0

          do q_ix = 1, CTrans%num_q_int 

                 if (.not.(CP%closed.and. nint(CTrans%q_int(q_ix)*CP%r)<=CTrans%ls%l(j))) then
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
        real(dl) pows(CTrans%num_q_int)
        real(dl)  ks(CTrans%num_q_int),measures(CTrans%num_q_int)

        do in=1,CP%InitPower%nn
        
         do q_ix = 1, CTrans%num_q_int 

               ks(q_ix) = CTrans%q_int(q_ix)
               measures(q_ix) = CTrans%dq_int(q_ix)/CTrans%q_int(q_ix)

          end do

         call GetInitPowers(pows,ks,CTrans%num_q_int,in)

        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
        !$OMP & PRIVATE(j,q_ix,measure,power,ctnorm,dbletmp)
         do j=1,CTrans%ls%l0

          do q_ix = 1, CTrans%num_q_int 

             if (.not.(CP%closed.and. nint(CTrans%q_int(q_ix)*CP%r)<=CTrans%ls%l(j))) then
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
  
      !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,in), SHARED(CTransS,CTransT),IF(CP%InitPower%nn > 1)
      do in=1,CP%InitPower%nn
          if (CP%WantScalars) then
               do i = C_Temp, C_last
                call InterpolateClArr(CTransS%ls,iCl_scalar(1,i,in),Cl_scalar(lmin, in, i),CTransS%ls%l0)
               end do
             end if
          
              if (CP%WantVectors) then
               do i = C_Temp, CT_cross
                call InterpolateClArr(CTransV%ls,iCl_vector(1,i,in),Cl_vector(lmin, in, i),CTransV%ls%l0)
               end do
             end if
             
             if (CP%WantTensors) then
               do i = CT_Temp, CT_Cross
                call InterpolateClArr(CTransT%ls,iCl_tensor(1,i,in),Cl_tensor(lmin, in, i),CTransT%ls%l0)
               end do
             end if
      end do
      !$OMP END PARALLEL DO
      end subroutine InterpolateCls


end module CAMBmain

    
     

    
