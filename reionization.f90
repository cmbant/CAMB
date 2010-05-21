
module Reionization
 use Precision
 use AMLutils
 implicit none  

!This module puts smooth tanh reionization of specified mid-point (z_{re}) and width
!The tanh function is in the variable (1+z)**Rionization_zexp

!Rionization_zexp=1.5 has the property that for the same z_{re} 
!the optical depth agrees with infinitely sharp model for matter domination 
!So tau and zre can be mapped into each other easily (for any symmetric window)
!However for generality the module maps tau into z_{re} using a binary search
!so could be easily modified for other monatonic parameterizations.

!AL March 2008
!AL July 2008 - added trap for setting optical depth without use_optical_depth

!See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

       character(LEN=*), parameter :: Reionization_Name = 'CAMB_reionization'

       real(dl), parameter :: Reionization_DefFraction = -1._dl 
        !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other
       
       real(dl) :: Reionization_AccuracyBoost = 1._dl
       real(dl) :: Rionization_zexp = 1.5_dl

       logical :: include_helium_fullreion = .true.
       real(dl) :: helium_fullreion_redshift  = 3.5_dl
       real(dl) :: helium_fullreion_deltaredshift  = 0.5
       real(dl) :: helium_fullreion_redshiftstart  = 5._dl
       
       
       type ReionizationParams
             logical    :: Reionization
             logical    :: use_optical_depth
             real(dl)   :: redshift, delta_redshift, fraction
             real(dl)   :: optical_depth
        end type ReionizationParams
  
        type ReionizationHistory
!These two are used by main code to bound region where xe changing 
          real(dl) :: tau_start, tau_complete
!This is set from main code          
          real(dl) :: akthom, fHe

!The rest are internal to this module.
          real(dl) :: WindowVarMid, WindowVarDelta

        end type ReionizationHistory

      real(dl), parameter :: Reionization_maxz = 40._dl
      real(dl), private, parameter :: Reionization_tol = 1d-5

      real(dl), private, external :: dtauda, rombint,rombint2

    Type(ReionizationParams), private, pointer ::  ThisReion
    Type(ReionizationHistory), private, pointer :: ThisReionHist

contains

 
 function Reionization_xe(a, tau, xe_recomb)
 !a and time tau and redundant, both provided for convenience
 !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
 !xe should map smoothly onto xe_recomb
  real(dl), intent(in) :: a
  real(dl), intent(in), optional :: tau, xe_recomb
  real(dl) Reionization_xe
  real(dl) tgh, xod
  real(dl) xstart
  
        if (present(xe_recomb)) then
        xstart = xe_recomb
        else
        xstart = 0._dl
        end if
        
        xod = (ThisReionHist%WindowVarMid - 1._dl/a**Rionization_zexp)/ThisReionHist%WindowVarDelta
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if
        Reionization_xe =(ThisReion%fraction-xstart)*(tgh+1._dl)/2._dl+xstart
     
        if (include_helium_fullreion .and. a > (1/(1+ helium_fullreion_redshiftstart))) then
        
         !Effect of Helium becoming fully ionized at z <~ 3.5 is very small so details not important
          xod = (1+helium_fullreion_redshift - 1._dl/a)/helium_fullreion_deltaredshift
          if (xod > 100) then
             tgh=1.d0
          else
             tgh=tanh(xod)
          end if
      
        Reionization_xe =  Reionization_xe + ThisReionHist%fHe*(tgh+1._dl)/2._dl
  
       end if
      
 end function Reionization_xe  
 
 function Reionization_timesteps(ReionHist)
 !minimum number of time steps to use between tau_start and tau_complete
 !Scaled by AccuracyBoost later
 !steps may be set smaller than this anyway
  Type(ReionizationHistory) :: ReionHist
  integer Reionization_timesteps
  
  Reionization_timesteps = 50 
  
 end  function Reionization_timesteps
 
 subroutine Reionization_ReadParams(Reion, Ini)
  use IniFile
  Type(ReionizationParams) :: Reion
  Type(TIniFile) :: Ini

   Reion%Reionization = Ini_Read_Logical_File(Ini,'reionization')
   if (Reion%Reionization) then
   
    Reion%use_optical_depth = Ini_Read_Logical_File(Ini,'re_use_optical_depth') 
  
    if (Reion%use_optical_depth) then
              Reion%optical_depth = Ini_Read_Double_File(Ini,'re_optical_depth')
    else 
              Reion%redshift = Ini_Read_Double_File(Ini,'re_redshift')
    end if 

    Reion%delta_redshift = Ini_Read_Double_File(Ini,'re_delta_redshift', 0.5_dl) !default similar to CMBFAST original
    Reion%fraction = Ini_Read_Double_File(Ini,'re_ionization_frac',Reionization_DefFraction)
  
  end if

 end subroutine Reionization_ReadParams 

 subroutine Reionization_SetParamsForZre(Reion,ReionHist)
  Type(ReionizationParams), target :: Reion
  Type(ReionizationHistory), target :: ReionHist 
     
       ReionHist%WindowVarMid = (1._dl+Reion%redshift)**Rionization_zexp
       ReionHist%WindowVarDelta = &
         Rionization_zexp*(1._dl+Reion%redshift)**(Rionization_zexp-1._dl)*Reion%delta_redshift
 
 end subroutine Reionization_SetParamsForZre

 subroutine Reionization_Init(Reion, ReionHist, Yhe, akthom, tau0, FeedbackLevel)
  use constants
  Type(ReionizationParams), target :: Reion
  Type(ReionizationHistory), target :: ReionHist
  real(dl), intent(in) :: akthom, tau0, Yhe 
  integer, intent(in) :: FeedbackLevel
  real(dl) astart

     ReionHist%akthom = akthom  
     ReionHist%fHe =  YHe/(mass_ratio_He_H*(1.d0-YHe))
 
     ReionHist%tau_start=tau0
     ReionHist%tau_complete=tau0
      
     ThisReion => Reion
     ThisReionHist => ReionHist

     if (Reion%Reionization) then
 
            if (Reion%optical_depth /= 0._dl .and. .not. Reion%use_optical_depth) &
             write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'
    

           if (Reion%use_optical_depth.and.Reion%optical_depth<0.001 &
                .or. .not.Reion%use_optical_depth .and. Reion%Redshift<0.001) then
               Reion%Reionization = .false.
           end if
          
      end if   
     
     if (Reion%Reionization) then
           
        if (Reion%fraction==Reionization_DefFraction) &
                 Reion%fraction = 1._dl + ReionHist%fHe  !H + singly ionized He
          
       if (Reion%use_optical_depth) then
        call Reionization_SetFromOptDepth(Reion,ReionHist)
        if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') Reion%redshift
       end if

       call Reionization_SetParamsForZre(ThisReion,ThisReionHist)
       
      !this is a check, agrees very well in default parameterization
       if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') &
            Reionization_GetOptDepth(Reion, ReionHist) 

      !Get relevant times       
       astart=1.d0/(1.d0+Reion%redshift + Reion%delta_redshift*8)
       ReionHist%tau_start = max(0.05_dl, rombint(dtauda,0._dl,astart,1d-3))
          !Time when a very small reionization fraction (assuming tanh fitting)

       ReionHist%tau_complete = min(tau0, &
          ReionHist%tau_start+ rombint(dtauda,astart,1.d0/(1.d0+max(0.d0,Reion%redshift-Reion%delta_redshift*8)),1d-3))

    end if   
       
 end subroutine Reionization_Init
 
 
 subroutine Reionization_SetDefParams(Reion)
  Type(ReionizationParams) :: Reion
 
       Reion%Reionization = .true.
       Reion%use_optical_depth = .false.
       Reion%optical_depth = 0._dl
       Reion%redshift = 10
       Reion%fraction = Reionization_DefFraction
       Reion%delta_redshift = 0.5_dl

 end subroutine Reionization_SetDefParams

 subroutine Reionization_Validate(Reion, OK)
  Type(ReionizationParams),intent(in) :: Reion
  logical, intent(inout) :: OK
 
      if (Reion%Reionization) then
        if (Reion%use_optical_depth) then
            if (Reion%optical_depth<0 .or. Reion%optical_depth > 0.9  .or. &
               include_helium_fullreion .and. Reion%optical_depth<0.01) then
             OK = .false.
             write(*,*) 'Optical depth is strange. You have:', Reion%optical_depth 
            end if
        else
            if (Reion%redshift < 0 .or. Reion%Redshift +Reion%delta_redshift*3 > Reionization_maxz .or. &
              include_helium_fullreion .and. Reion%redshift < helium_fullreion_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',Reion%Redshift
            end if
        end if
        if (Reion%fraction/= Reionization_DefFraction .and. (Reion%fraction < 0 .or. Reion%fraction > 1.5)) then
                OK = .false.
                write(*,*) 'Reionization fraction strange. You have: ',Reion%fraction
        end if
        if (Reion%delta_redshift > 3 .or. Reion%delta_redshift<0.1 ) then
        !Very narrow windows likely to cause problems in interpolation etc.
        !Very broad likely to conflic with quasar data at z=6
                OK = .false.
                write(*,*) 'Reionization delta_redshift is strange. You have: ',Reion%delta_redshift
        end if


      end if
             
  end  subroutine Reionization_Validate


 function Reionization_doptdepth_dz(z)
   real(dl) :: Reionization_doptdepth_dz
   real(dl), intent(in) :: z
   real(dl) a
   
   a = 1._dl/(1._dl+z)
   
   Reionization_doptdepth_dz = Reionization_xe(a)*ThisReionHist%akthom*dtauda(a)

 end function Reionization_doptdepth_dz

function Reionization_GetOptDepth(Reion, ReionHist) 
 Type(ReionizationParams), target :: Reion
 Type(ReionizationHistory), target :: ReionHist
 real(dl) Reionization_GetOptDepth      
  
  ThisReion => Reion
  ThisReionHist => ReionHist
  Reionization_GetOptDepth = rombint2(Reionization_doptdepth_dz,0.d0,Reionization_maxz,&
         Reionization_tol, 20, nint(Reionization_maxz/Reion%delta_redshift*5))

end function Reionization_GetOptDepth

 subroutine Reionization_zreFromOptDepth(Reion, ReionHist)
 !General routine to find zre parameter given optical depth
 !Not used for Rionization_zexp = 1.5
  Type(ReionizationParams) :: Reion
  Type(ReionizationHistory) :: ReionHist
  real(dl) try_b, try_t
  real(dl) tau
  integer i
  
  try_b = 0
  try_t = Reionization_maxz
  i=0
  do 
       i=i+1  
       Reion%redshift = (try_t + try_b)/2
       call Reionization_SetParamsForZre(Reion,ReionHist)
       tau = Reionization_GetOptDepth(Reion, ReionHist)
       
       if (tau > Reion%optical_depth) then
                  try_t = Reion%redshift
          else
                  try_b = Reion%redshift
       end if
       if (abs(try_b - try_t) < 2e-3/Reionization_AccuracyBoost) exit
       if (i>100) call mpiStop('Reionization_zreFromOptDepth: failed to converge')        
  end do
  
  
   if (abs(tau - Reion%optical_depth) > 0.002) then
    write (*,*) 'Reionization_zreFromOptDepth: Did not converge to optical depth'
    write (*,*) 'tau =',tau, 'optical_depth = ', Reion%optical_depth
    write (*,*) try_t, try_b
    call mpiStop()
  end if
    
 end subroutine Reionization_zreFromOptDepth 


 
 subroutine Reionization_SetFromOptDepth(Reion, ReionHist)
  Type(ReionizationParams) :: Reion
  Type(ReionizationHistory) :: ReionHist
   
! This subroutine calculates the redshift of reionization

! This implementation is approximate but quite accurate and fast

      real(dl) dz, optd
      real(dl) z, tmp, tmpHe
      integer na
      
      Reion%redshift = 0

      if (Reion%Reionization .and. Reion%optical_depth /= 0) then

           !Do binary search to find zre from z
           !This is general method
           call Reionization_zreFromOptDepth(Reion, ReionHist)

        if (.false.) then
          !Use equivalence with sharp for special case
            optd=0
            na=1
            dz=1._dl/2000/Reionization_AccuracyBoost
            tmp = dz*Reion%fraction*ThisReionHist%akthom
            tmpHe = dz*(Reion%fraction+ReionHist%fHe)*ThisReionHist%akthom
            z=0
            do while (optd < Reion%optical_depth)
                z=na*dz
                if (include_helium_fullreion .and. z < helium_fullreion_redshift) then
                optd=optd+ tmpHe*dtauda(1._dl/(1._dl+z))
                else
                optd=optd+tmp*dtauda(1._dl/(1._dl+z))
                end if
                na=na+1
            end do
         end if
      else
         Reion%Reionization = .false.
      end if
 
 end  subroutine Reionization_SetFromOptDepth 



end module Reionization

