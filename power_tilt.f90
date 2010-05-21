!This module provides the initial power spectra, parameterized as an expansion in ln k
!
! ln P = ln A_s + (n_s -1)*ln(k/k_0) + n_{run}/2 * ln(k/k_0)^2 
!
! so if n_run = 0
!
! P = A_s (k/k_0_scalar)^(n_s-1)
!
!for the scalar spectrum, when an(in) is the in'th spectral index. k_0_scalar
!is a pivot scale, fixed here to 0.05/Mpc (change it below as desired).
!
!This module uses the same inputs an(in), ant(in) and rat(in) as CMBFAST, however here
!rat(in) is used to set the ratio of the initial power spectra, so here
!
!** rat(in) is not the Cl quadrupole ratio ***
!
!in general models the quadrupole ratio depends in a complicated way on the ratio of the initial
!power spectra

!The absolute normalization of the Cls is unimportant here, but the relative ratio
!of the tensor and scalar Cls generated with this module will be correct for general models


!The OutputNormalization parameter controls the final output
!Absolute Cls can be obtained by setting OuputNormalization=outNone, otherwise the overall normalization
!of the power spectra doesn't matter

!This version December 2003 - changed default tensor pivot to 0.05 (consistent with CMBFAST 4.5)

     module InitialPower   
     use Precision
     implicit none   

     private
    
      character(LEN=*), parameter :: Power_Name = 'power_tilt'

      integer, parameter :: nnmax= 5 
       !Maximum possible number of different power spectra to use

      Type InitialPowerParams

       integer nn  !Must have this variable
      !The actual number of power spectra to use
  
      !For the default implementation return power spectra based on spectral indices
       real(dl) an(nnmax) !scalar spectral indices
       real(dl) n_run(nnmax) !running of spectral index 
       real(dl) ant(nnmax) !tensor spectral indices
       real(dl) rat(nnmax) !ratio of scalar to tensor initial power spectrum amplitudes
       real(dl) k_0_scalar, k_0_tensor
       real(dl) ScalarPowerAmp(nnmax)
 
      end Type InitialPowerParams

      real(dl) curv  !Curvature contant, set in InitializePowers     
      
      Type(InitialPowerParams) :: P
  
!Make things visible as neccessary...
 
      public InitialPowerParams, InitialPower_ReadParams, InitializePowers, ScalarPower, TensorPower
      public nnmax,Power_Descript, Power_Name, SetDefPowerParams
!      public 
    contains
       
  
       subroutine SetDefPowerParams(AP)
        Type (InitialPowerParams) :: AP

         AP%nn     = 1 !number of initial power spectra
         AP%an     = 1 !scalar spectral index
         AP%n_run   = 0 !running of scalar spectral index
         AP%ant    = 0 !tensor spectra index
         AP%rat    = 1
         AP%k_0_scalar = 0.05
         AP%k_0_tensor = 0.05
         AP%ScalarPowerAmp = 1

       end subroutine SetDefPowerParams

       subroutine InitializePowers(AParamSet,acurv)
         Type (InitialPowerParams) :: AParamSet
         !Called before computing final Cls in cmbmain.f90
         !Could read spectra from disk here, do other processing, etc.

        real(dl) acurv

        if (AParamSet%nn > nnmax) then
           write (*,*) 'To use ',AParamSet%nn,'power spectra you need to increase'
           write (*,*) 'nnmax in power_tilt.f90, currently ',nnmax
        end if
        P = AParamSet

        curv=acurv         

!Write implementation specific code here...

       end subroutine InitializePowers
       

      function ScalarPower(k,in)

       !"in" gives the index of the power to return for this k
       !ScalarPower = const for scale invariant spectrum
       !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci 
       !scalar on co-moving hypersurfaces receives power
       ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k) 
       !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
       !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
       !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
       !Near the end of inflation chi is equal to 3/2 Psi.
       !Here nu^2 = (k^2 + curv)/|curv| 

       !This power spectrum is also used for isocurvature modes where 
       !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
       !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.

        real(dl) ScalarPower,k, lnrat
        integer in

          lnrat = log(k/P%k_0_scalar)
          ScalarPower=P%ScalarPowerAmp(in)*exp((P%an(in)-1)*lnrat + P%n_run(in)/2*lnrat**2)   
     
!         ScalarPower = ScalarPower * (1 + 0.1*cos( lnrat*30 ) )
 
      end function ScalarPower

      
      function TensorPower(k,in)
      
       !TensorPower= const for scale invariant spectrum
       !The normalization is defined so that
       ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
       !for a closed model
       ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
       !for an open model
       !"in" gives the index of the power spectrum to return 
       !Here nu^2 = (k^2 + 3*curv)/|curv| 


        real(dl) TensorPower,k   
        real(dl), parameter :: PiByTwo=3.14159265d0/2._dl
   
        integer in

        TensorPower=P%rat(in)*P%ScalarPowerAmp(in)*exp(P%ant(in)*log(k/P%k_0_tensor))
        if (curv < 0) TensorPower=TensorPower*tanh(PiByTwo*sqrt(-k**2/curv-3)) 

       
      end function TensorPower

      !Get parameters describing parameterisation (for FITS file)
     function Power_Descript(in, Scal, Tens, Keys, Vals)
         character(LEN=8), intent(out) :: Keys(*)
         real(dl), intent(out) :: Vals(*)
         integer, intent(IN) :: in
         logical, intent(IN) :: Scal, Tens
         integer num, Power_Descript
         num=0
         if (Scal) then
         num=num+1
         Keys(num) = 'n_s'
         Vals(num) = P%an(in)
         num=num+1
         Keys(num) = 'n_run'
         Vals(num) = P%n_run(in)
         num=num+1
         Keys(num) = 's_pivot'
         Vals(num) = P%k_0_scalar
         end if
         if (Tens) then
         num=num+1
         Keys(num) = 'n_t'
         Vals(num) = P%ant(in)
         num=num+1
         Keys(num) = 't_pivot'
         Vals(num) = P%k_0_tensor
         if (Scal) then
           num=num+1
           Keys(num) = 'p_ratio'
           Vals(num) = P%rat(in)
         end if
         end if
         Power_Descript = num

       end  function Power_Descript
        
       subroutine InitialPower_ReadParams(InitPower, Ini, WantTensors)
          use IniFile
          Type(InitialPowerParams) :: InitPower
          Type(TIniFile) :: Ini
          logical, intent(in) :: WantTensors
          integer i
          
           InitPower%k_0_scalar = Ini_Read_Double_File(Ini,'pivot_scalar',InitPower%k_0_scalar)
           InitPower%k_0_tensor = Ini_Read_Double_File(Ini,'pivot_tensor',InitPower%k_0_tensor) 
           InitPower%nn = Ini_Read_Int('initial_power_num')
           if (InitPower%nn>nnmax) stop 'Too many initial power spectra - increase nnmax in InitialPower'
           InitPower%rat(:) = 1
           do i=1, InitPower%nn

              InitPower%an(i) = Ini_Read_Double_Array_File(Ini,'scalar_spectral_index', i)
              InitPower%n_run(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrun',i,0._dl) 
              
              if (WantTensors) then
                 InitPower%ant(i) = Ini_Read_Double_Array_File(Ini,'tensor_spectral_index',i)
                 InitPower%rat(i) = Ini_Read_Double_Array_File(Ini,'initial_ratio',i)
              end if              

              InitPower%ScalarPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'scalar_amp',i,1._dl) 
              !Always need this as may want to set tensor amplitude even if scalars not computed
           end do
          
       end  subroutine InitialPower_ReadParams 


     end module InitialPower
