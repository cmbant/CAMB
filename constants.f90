

   module Precision
      implicit none

      integer, parameter :: dl = KIND(1.d0)
      integer, parameter :: sp = KIND(1.0)
 
   end module Precision



   module constants
       use precision
       implicit none

      
       real(dl), parameter :: const_pi = 3.1415926535897932384626433832795_dl
       real(dl), parameter :: const_twopi=2._dl*const_pi, const_fourpi=4._dl*const_pi
       real(dl), parameter :: const_sqrt6=2.4494897427831780981972840747059_dl

       real(dl), parameter :: c = 2.99792458e8_dl
       real(dl), parameter :: h_P = 6.62606896e-34_dl
        
       real(dl), parameter :: G=6.67428e-11_dl
       real(dl), parameter :: sigma_thomson = 6.6524616e-29_dl
       real(dl), parameter :: sigma_boltz = 5.6704e-8_dl  
       real(dl), parameter :: k_B = 1.3806504e-23_dl 


       real(dl), parameter :: m_p = 1.672621637e-27_dl  ! 1.672623e-27_dl
       real(dl), parameter :: m_H = 1.673575e-27_dl !av. H atom
       real(dl), parameter :: m_e = 9.10938215e-31_dl
       real(dl), parameter :: mass_ratio_He_H = 3.9715_dl
 

       real(dl), parameter :: Gyr=3.1556926e16_dl
       real(dl), parameter :: Mpc = 3.085678e22_dl !seem to be different definitions of this?
       real(dl), parameter :: MPC_in_sec = Mpc/c ! Mpc/c = 1.029272d14 in SI units      

       real(dl), parameter :: barssc0= k_B / m_p / c**2
       real(dl), parameter :: kappa=8._dl*const_pi*G
       real(dl), parameter :: a_rad = 8._dl*const_pi**5*k_B**4/15/c**3/h_p**3  
                !7.565914e-16_dl !radiation constant for u=aT^4


       real(dl), parameter :: Compton_CT = MPC_in_sec*(8.d0/3.d0)*(sigma_thomson/(m_e*c))*a_rad
        !Compton_CT is CT in Mpc units, (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc
        !Used to get evolution of matter temperature 
        
       !For 21cm
       real(dl), parameter :: f_21cm = 1420.40575e6_dl, l_21cm= c/f_21cm, T_21cm = h_P*f_21cm/k_B
       real(dl), parameter :: A10 = 2.869e-15, B10 = l_21cm**3/2/h_P/c*A10

       real(dl), parameter :: line21_const = 3*l_21cm**2*C*h_P/32/const_pi/k_B*A10 * Mpc_in_sec * 1000
        !1000 to get in MiliKelvin


   end module constants



