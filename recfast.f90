    !Recombination module for CAMB, using RECFAST

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !C Integrator for Cosmic Recombination of Hydrogen and Helium,
    !C developed by Douglas Scott (dscott@astro.ubc.ca)
    !C based on calculations in the paper Seager, Sasselov & Scott
    !C (ApJ, 523, L1, 1999).
    !and "fudge" updates in Wong, Moss & Scott (2008).
    !C
    !C Permission to use, copy, modify and distribute without fee or royalty at
    !C any tier, this software and its documentation, for any purpose and without
    !C fee or royalty is hereby granted, provided that you agree to comply with
    !C the following copyright notice and statements, including the disclaimer,
    !C and that the same appear on ALL copies of the software and documentation,
    !C including modifications that you make for internal use or for distribution:
    !C
    !C Copyright 1999-2010 by University of British Columbia.  All rights reserved.
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
    !CV     Version: 1.5.2
    !C
    !CP     Purpose:  Calculate ionised fraction as a function of redshift.
    !CP            Solves for H and He simultaneously, and includes
    !CP           H "fudge factor" for low z effect, as well as
    !CP           HeI fudge factor.
    !C
    !CD     Description: Solves for ionisation history since recombination
    !CD     using the equations in Seager, Sasselov & Scott (ApJ, 1999).
    !CD     The Cosmological model can be flat or open.
    !CD	 The matter temperature is also followed, with an update from
    !CD	 Scott & Scott (2009).
    !CD	 The values for \alpha_B for H are from Hummer (1994).
    !CD	 The singlet HeI coefficient is a fit from the full code.
    !CD	 Additional He "fudge factors" are as described in Wong, Moss
    !CD	 and Scott (2008).
    !CD	 Extra fitting function included (in optical depth) to account
    !CD	 for extra H physics described in Rubino-Martin et al. (2010).
    !CD	 Care is taken to use the most accurate constants.
    !C
    !CA     Arguments:
    !CA     Name, Description
    !CA     real(dl) throughout
    !CA
    !CA     z is redshift - W is sqrt(1+z), like conformal time
    !CA     x is total ionised fraction, relative to H
    !CA     x_H is ionized fraction of H - y(1) in R-K routine
    !CA     x_He is ionized fraction of He - y(2) in R-K routine
    !CA       (note that x_He=n_He+/n_He here and not n_He+/n_H)
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
    !CA	    epsilon is the approximate difference (=Trad-Tmat) at high z
    !CA     OmegaB is Omega in baryons today
    !CA     H is Hubble constant in units of 100 km/s/Mpc
    !CA     HO is Hubble constant in SI units
    !CA     bigH is 100 km/s/Mpc in SI units
    !CA	    Hz is the value of H at the specific z (in ION)
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
    !CA     C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
    !CA     m_e,m_H: electron mass and mass of H atom in SI
    !CA     not4: ratio of 4He atomic mass to 1H atomic mass
    !CA     sigma: Thomson cross-section
    !CA     a_rad: radiation constant for u=aT^4
    !CA     Lambda: 2s-1s two photon rate for Hydrogen
    !CA     Lambda_He: 2s-1s two photon rate for Helium
    !CA     DeltaB: energy of first excited state from continuum = 3.4eV
    !CA     DeltaB_He: energy of first excited state from cont. for He = 3.4eV
    !CA     L_H_ion: level for H ionization in m^-1
    !CA     L_H_alpha: level for H Ly alpha in m^-1
    !CA     L_He1_ion: level for HeI ionization
    !CA     L_He2_ion: level for HeII ionization
    !CA     L_He_2s: level for HeI 2s
    !CA     L_He_2p: level for HeI 2p (21P1-11S0) in m^-1
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
    !CA b_He= "fudge factor" for HeI, to approximate higher z behaviour
    !CA Heswitch=integer for modifying HeI recombination
    !CA Parameters and quantities to describe the extra triplet states
    !CA  and also the continuum opacity of H, with a fitting function
    !CA  suggested by KIV, astro-ph/0703438
    !CA a_trip: used to fit HeI triplet recombination rate
    !CA b_trip: used to fit HeI triplet recombination rate
    !CA L_He_2Pt: level for 23P012-11S0 in m^-1
    !CA L_He_2St: level for 23S1-11S0 in m^-1
    !CA L_He2St_ion: level for 23S1-continuum in m^-1
    !CA A2P_s: Einstein A coefficient for He 21P1-11S0
    !CA A2P_t: Einstein A coefficient for He 23P1-11S0
    !CA sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
    !CA sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
    !CA CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
    !CA CfHe_t: triplet statistical correction
    !CA	Hswitch is an boolean for modifying the H recombination
    !CA	AGauss1 is the amplitude of the 1st Gaussian for the H fudging
    !CA	AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
    !CA	zGauss1 is the ln(1+z) central value of the 1st Gaussian
    !CA	zGauss2 is the ln(1+z) central value of the 2nd Gaussian
    !CA	wGauss1 is the width of the 1st Gaussian
    !CA	wGauss2 is the width of the 2nd Gaussian


    !CA     tol: tolerance for the integrator
    !CA     cw(24),w(3,9): work space for DVERK
    !CA     Ndim: number of d.e.'s to solve (integer)
    !CA     Nz: number of output redshitf (integer)
    !CA     I: loop index (integer)
    !CA     ind,nw: work-space for DVERK (integer)
    !C
    !CF     File & device access:
    !CF     Unit /I,IO,O  /Name (if known)
    !C
    !CM     Modules called:
    !CM     DVERK (numerical integrator)
    !CM     GET_INIT (initial values for ionization fractions)
    !CM     ION (ionization and Temp derivatives)
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
    !                       October 2006 (included new value for G)
    !                       October 2006 (improved m_He/m_H to be "not4")
    !                       October 2006 (fixed error, x for x_H in part of f(1))
    !CH              January 2008 (improved HeI recombination effects,
    !CH                       including HeI rec. fudge factor)
    !                Feb 2008   Recfast 1.4 changes above added (AML)
    !                           removed Dubrovich option (wrong anyway)
    !CH   			 Sept 2008 (added extra term to make transition, smoother for Tmat evolution)
    !                Sept 2008 Recfast 1.4.2 changes above added (AML)
    !                          General recombination module structure, fix to make He x_e smooth also in recfast (AML)
    !CH		 Jan 2010 (added fitting function to modify K
    !CH			 	 to match x_e(z) for new H physics)
    !AL             June 2012 updated fudge parameters to match HyRec and CosmoRec (AML)
    !AL             Sept 2012 changes now in public recfast, version number changed to match Recfast 1.5.2.

    !!      ===============================================================

    module RECDATA
    use constants
    implicit none


    real(dl) Lambda,DeltaB,DeltaB_He,Lalpha,mu_H,mu_T,H_frac
    real(dl) Lambda_He,Lalpha_He,Bfact,CK_He,CL_He
    real(dl) L_H_ion,L_H_alpha,L_He1_ion,L_He2_ion,L_He_2s,L_He_2p
    real(dl) CB1,CDB,CR,CK,CL,CT,CB1_He1,CB1_He2,CDB_He,fu
    real(dl) A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt
    real(dl)  L_He_2Pt,L_He_2St,L_He2St_ion


    real(dl), parameter :: bigH=100.0D3/Mpc !Ho in s-1
    real(dl), parameter :: sigma = sigma_thomson
    real(dl), parameter :: not4  = mass_ratio_He_H    !mass He/H atom

    real(dl) Tnow,HO
    integer :: n_eq = 3

    !The following only used for approximations where small effect
    real(dl) OmegaK, OmegaT, z_eq


    !Fundamental constants in SI units
    !      ("not4" pointed out by Gary Steigman)

    data    Lambda      /8.2245809d0/
    data    Lambda_He   /51.3d0/    !new value from Dalgarno
    data    L_H_ion     /1.096787737D7/ !level for H ion. (in m^-1)
    data    L_H_alpha   /8.225916453D6/ !averaged over 2 levels
    data    L_He1_ion   /1.98310772D7/  !from Drake (1993)
    data    L_He2_ion   /4.389088863D7/ !from JPhysChemRefData (1987)
    data    L_He_2s     /1.66277434D7/  !from Drake (1993)
    data    L_He_2p     /1.71134891D7/  !from Drake (1993)
    !   2 photon rates and atomic levels in SI units

    data    A2P_s       /1.798287D9/    !Morton, Wu & Drake (2006)
    data    A2P_t       /177.58D0/      !Lach & Pachuski (2001)
    data    L_He_2Pt    /1.690871466D7/ !Drake & Morton (2007)
    data    L_He_2St    /1.5985597526D7/ !Drake & Morton (2007)
    data    L_He2St_ion /3.8454693845D6/ !Drake & Morton (2007)
    data    sigma_He_2Ps    /1.436289D-22/  !Hummer & Storey (1998)
    data    sigma_He_2Pt    /1.484872D-22/  !Hummer & Storey (1998)
    !    Atomic data for HeI


    end module RECDATA


    module Recombination
    use constants
    use AMLUtils
    implicit none
    private

    real(dl), parameter ::  zinitial = 1e4_dl !highest redshift
    real(dl), parameter ::  zfinal=0._dl
    integer,  parameter :: Nz=10000
    real(dl), parameter :: delta_z = (zinitial-zfinal)/Nz

    integer, parameter ::  RECFAST_Heswitch_default = 6
    real(dl), parameter :: RECFAST_fudge_He_default = 0.86_dl !Helium fudge parameter
    logical, parameter  :: RECFAST_Hswitch_default = .true. !include H corrections (v1.5, 2010)
    real(dl), parameter :: RECFAST_fudge_default = 1.14_dl !1.14_dl
    real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 + 0.02d0
    !fudge parameter if RECFAST_Hswitch

    real(dl) :: AGauss1 =      -0.14D0  !Amplitude of 1st Gaussian
    real(dl) :: AGauss2 =       0.079D0 ! 0.05D0  !Amplitude of 2nd Gaussian
    real(dl) :: zGauss1 =       7.28D0  !ln(1+z) of 1st Gaussian
    real(dl) :: zGauss2=        6.73D0  !ln(1+z) of 2nd Gaussian
    real(dl) :: wGauss1=        0.18D0  !Width of 1st Gaussian
    real(dl) :: wGauss2=        0.33D0  !Width of 2nd Gaussian
    !Gaussian fits for extra H physics (fit by Adam Moss , modified by Antony Lewis)

    type RecombinationParams

        real(dl) :: RECFAST_fudge
        real(dl) :: RECFAST_fudge_He
        integer  :: RECFAST_Heswitch
        logical  :: RECFAST_Hswitch
        !0) no change from old Recfast'
        !1) full expression for escape probability for singlet'
        !'   1P-1S transition'
        !2) also including effect of contiuum opacity of H on HeI'
        !'   singlet (based in fitting formula suggested by'
        !'   Kholupenko, Ivanchik & Varshalovich, 2007)'
        !3) only including recombination through the triplets'
        !4) including 3 and the effect of the contiuum '
        !'   (although this is probably negligible)'
        !5) including only 1, 2 and 3'
        !6) including all of 1 to 4'

    end  type RecombinationParams

    character(LEN=*), parameter :: Recombination_Name = 'Recfast_1.5.2'

    real(dl) zrec(Nz),xrec(Nz),dxrec(Nz), Tsrec(Nz) ,dTsrec(Nz), tmrec(Nz),dtmrec(Nz)

    real(dl), parameter :: Do21cm_mina = 1/(1+900.) !at which to start evolving Delta_TM
    logical, parameter :: evolve_Ts = .false. !local equilibrium is very accurate
    real(dl), parameter :: Do21cm_minev = 1/(1+400.) !at which to evolve T_s


    real(dl), parameter :: B01 = 3*B10
    real(dl) :: NNow, fHe


    logical :: Do21cm = .false.
    logical :: doTmatTspin = .false.

    real(dl) :: recombination_saha_z !Redshift at which saha OK
    real(dl) :: recombination_saha_tau !set externally


    public RecombinationParams, Recombination_xe, Recombination_tm,Recombination_ts ,Recombination_init,   &
        Recombination_ReadParams, Recombination_SetDefParams, Recombination_Validate, Recombination_Name, &
        kappa_HH_21cm,kappa_eH_21cm,kappa_pH_21cm, &
        Do21cm, doTmatTspin, Do21cm_mina, dDeltaxe_dtau, &
        recombination_saha_tau, recombination_saha_z

    contains



    subroutine Recombination_ReadParams(R, Ini)
    use IniFile
    Type(RecombinationParams) :: R
    Type(TIniFile) :: Ini


    R%RECFAST_fudge_He = Ini_Read_Double_File(Ini,'RECFAST_fudge_He',RECFAST_fudge_He_default)
    R%RECFAST_Heswitch = Ini_Read_Int_File(Ini, 'RECFAST_Heswitch',RECFAST_Heswitch_default)
    R%RECFAST_Hswitch = Ini_Read_Logical_File(Ini, 'RECFAST_Hswitch',RECFAST_Hswitch_default)
    R%RECFAST_fudge = Ini_Read_Double_File(Ini,'RECFAST_fudge',RECFAST_fudge_default)
    AGauss1 = Ini_REad_Double_File(Ini,'AGauss1',AGauss1)
    AGauss2 = Ini_REad_Double_File(Ini,'AGauss2',AGauss2)
    zGauss1 = Ini_REad_Double_File(Ini,'zGauss1',zGauss1)
    zGauss2 = Ini_REad_Double_File(Ini,'zGauss2',zGauss2)
    wGauss1 = Ini_REad_Double_File(Ini,'wGauss1',wGauss1)
    wGauss2 = Ini_REad_Double_File(Ini,'wGauss2',wGauss2)
    if (R%RECFAST_Hswitch) then
        R%RECFAST_fudge = R%RECFAST_fudge - (RECFAST_fudge_default - RECFAST_fudge_default2)
    end if
    end subroutine Recombination_ReadParams

    subroutine Recombination_SetDefParams(R)
    type (RecombinationParams) ::R


    R%RECFAST_fudge = RECFAST_fudge_default
    R%RECFAST_fudge_He = RECFAST_fudge_He_default !Helium fudge parameter
    R%RECFAST_Heswitch = RECFAST_Heswitch_default
    R%RECFAST_Hswitch =  RECFAST_Hswitch_default
    if (R%RECFAST_Hswitch) then
        R%RECFAST_fudge = RECFAST_fudge_default2
    end if

    end subroutine Recombination_SetDefParams


    subroutine Recombination_Validate(R, OK)
    Type(RecombinationParams),intent(in) :: R
    logical, intent(inout) :: OK

    if (R%RECFAST_Heswitch<0 .or. R%RECFAST_Heswitch > 6) then
        OK = .false.
        write(*,*) 'RECFAST_Heswitch unknown'
    end if

    end subroutine Recombination_Validate


    function Recombination_tm(a)
    use RECDATA, only : Tnow
    real(dl) zst,a,z,az,bz,Recombination_tm
    integer ilo,ihi

    if (.not. doTmatTspin) call MpiStop('RECFAST: Recombination_tm not stored')
    z=1/a-1
    if (z >= zrec(1)) then
        Recombination_tm=Tnow/a
    else
        if (z <=zrec(nz)) then
            Recombination_tm=Tmrec(nz)
        else
            zst=(zinitial-z)/delta_z
            ihi= int(zst)
            ilo = ihi+1
            az=zst - int(zst)
            bz=1-az
            Recombination_tm=az*Tmrec(ilo)+bz*Tmrec(ihi)+ &
                ((az**3-az)*dTmrec(ilo)+(bz**3-bz)*dTmrec(ihi))/6._dl
        endif
    endif

    end function Recombination_tm


    function Recombination_ts(a)
    !zrec(1) is zinitial-delta_z
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,Recombination_ts
    integer ilo,ihi

    z=1/a-1
    if (z.ge.zrec(1)) then
        Recombination_ts=tsrec(1)
    else
        if (z.le.zrec(nz)) then
            Recombination_ts=tsrec(nz)
        else
            zst=(zinitial-z)/delta_z
            ihi= int(zst)
            ilo = ihi+1
            az=zst - int(zst)
            bz=1-az

            Recombination_ts=az*tsrec(ilo)+bz*tsrec(ihi)+ &
                ((az**3-az)*dtsrec(ilo)+(bz**3-bz)*dtsrec(ihi))/6._dl
        endif
    endif

    end function Recombination_ts


    function Recombination_xe(a)
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,Recombination_xe
    integer ilo,ihi

    z=1/a-1
    if (z.ge.zrec(1)) then
        Recombination_xe=xrec(1)
    else
        if (z.le.zrec(nz)) then
            Recombination_xe=xrec(nz)
        else
            zst=(zinitial-z)/delta_z
            ihi= int(zst)
            ilo = ihi+1
            az=zst - int(zst)
            bz=1-az
            Recombination_xe=az*xrec(ilo)+bz*xrec(ihi)+ &
                ((az**3-az)*dxrec(ilo)+(bz**3-bz)*dxrec(ihi))/6._dl
        endif
    endif

    end function Recombination_xe



    subroutine Recombination_init(Recomb, OmegaC, OmegaB, Omegan, Omegav, h0inp,tcmb,yp, nnu)
    !Would love to pass structure as arguments, but F90 would give circular reference...
    !hence mess passing parameters explcitly and non-generally
    !Note recfast only uses OmegaB, h0inp, tcmb and yp - others used only for Tmat approximation where effect small
    !nnu currently not used here
    use RECDATA
    use AMLUtils
    implicit none
    Type (RecombinationParams) :: Recomb

    real(dl), save :: last_OmB =0, Last_YHe=0, Last_H0=0, Last_dtauda=0, last_fudge, last_fudgeHe

    real(dl) Trad,Tmat,Tspin,d0hi,d0lo
    integer I

    real(dl), intent(in) :: OmegaB,OmegaC, Omegan, Omegav, h0inp, yp
    real(dl), intent(in), optional :: nnu
    real(dl) z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0,H
    real(dl) zstart,zend,tcmb
    real(dl) cw(24)
    real(dl), dimension(:,:), allocatable :: w
    real(dl) y(4)
    real(dl) C10, tau_21Ts
    real(dl) fnu
    integer ind,nw

    !       --- Parameter statements
    real(dl), parameter :: tol=1.D-5                !Tolerance for R-K

    real(dl) dtauda
    external dtauda, dverk

    !       ===============================================================

    if (Last_OmB==OmegaB .and. Last_H0 == h0inp .and. yp == Last_YHe .and. &
        dtauda(0.2352375823_dl) == Last_dtauda .and. last_fudge == Recomb%RECFAST_fudge &
        .and. last_fudgeHe==Recomb%RECFAST_fudge_He) return
    !This takes up most of the single thread time, so cache if at all possible
    !For example if called with different reionization, or tensor rather than scalar

    Last_dtauda =  dtauda(0.2352375823_dl) !Just get it at a random scale factor
    Last_OmB = OmegaB
    Last_H0 = h0inp
    Last_YHe=yp
    last_fudge = Recomb%RECFAST_FUDGE
    last_fudgeHe = Recomb%RECFAST_FUDGE_He

    if (Do21cm) doTmatTspin = .true.


    !       write(*,*)'recfast version 1.0'
    !       write(*,*)'Using Hummer''s case B recombination rates for H'
    !       write(*,*)' with fudge factor = 1.14'
    !       write(*,*)'and tabulated HeII singlet recombination rates'
    !       write(*,*)

    n_eq = 3
    if (Evolve_Ts) n_eq=4
    allocate(w(n_eq,9))

    recombination_saha_z=0.d0

    Tnow=tcmb
    !       These are easy to inquire as input, but let's use simple values
    z = zinitial
    !       will output every 1 in z, but this is easily changed also

    !Not general, but only for approx
    OmegaT=OmegaC+OmegaB            !total dark matter + baryons
    OmegaK=1.d0-OmegaT-OmegaV       !curvature


    !       convert the Hubble constant units
    H = H0inp/100._dl
    HO = H*bigH


    !       sort out the helium abundance parameters
    mu_H = 1.d0/(1.d0-Yp)           !Mass per H atom
    mu_T = not4/(not4-(not4-1.d0)*Yp)   !Mass per atom
    fHe = Yp/(not4*(1.d0-Yp))       !n_He_tot / n_H_tot


    Nnow = 3._dl*HO*HO*OmegaB/(8._dl*Pi*G*mu_H*m_H)

    n = Nnow * (1._dl+z)**3
    fnu = (21.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
    !	(this is explictly for 3 massless neutrinos - change if N_nu.ne.3; but only used for approximation so not critical)
    z_eq = (3.d0*(HO*C)**2/(8.d0*Pi*G*a_rad*(1.d0+fnu)*Tnow**4))*(OmegaB+OmegaC)
    z_eq = z_eq - 1.d0


    !       Set up some constants so they don't have to be calculated later
    Lalpha = 1.d0/L_H_alpha
    Lalpha_He = 1.d0/L_He_2p
    DeltaB = h_P*C*(L_H_ion-L_H_alpha)
    CDB = DeltaB/k_B
    DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)   !2s, not 2p
    CDB_He = DeltaB_He/k_B
    CB1 = h_P*C*L_H_ion/k_B
    CB1_He1 = h_P*C*L_He1_ion/k_B   !ionization for HeI
    CB1_He2 = h_P*C*L_He2_ion/k_B   !ionization for HeII
    CR = 2.d0*Pi*(m_e/h_P)*(k_B/h_P)
    CK = Lalpha**3/(8.d0*Pi)
    CK_He = Lalpha_He**3/(8.d0*Pi)
    CL = C*h_P/(k_B*Lalpha)
    CL_He = C*h_P/(k_B/L_He_2s) !comes from det.bal. of 2s-1s
    CT = Compton_CT / MPC_in_sec

    Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B


    !       Matter departs from radiation when t(Th) > H_frac * t(H)
    !       choose some safely small number
    H_frac = 1D-3

    !       Fudge factor to approximate for low z out of equilibrium effect
    fu=Recomb%RECFAST_fudge

    !       Set initial matter temperature
    y(3) = Tnow*(1._dl+z)            !Initial rad. & mat. temperature
    Tmat = y(3)
    y(4) = Tmat
    Tspin = Tmat

    call get_init(z,x_H0,x_He0,x0)

    y(1) = x_H0
    y(2) = x_He0

    !       OK that's the initial conditions, now start writing output file


    !       Set up work-space stuff for DVERK
    ind  = 1
    nw   = n_eq
    do i = 1,24
        cw(i) = 0._dl
    end do

    do i = 1,Nz
        !       calculate the start and end redshift for the interval at each z
        !       or just at each z
        zstart = zinitial  - real(i-1,dl)*delta_z
        zend   = zinitial  - real(i,dl)*delta_z

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
            y(4) = y(3)

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
            y(4) = y(3)

        else if(z > 3500._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = x_H0 + fHe*x_He0
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)
            y(4) = y(3)

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
            y(4) = y(3)

        else if (y(1) > 0.99d0) then

            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1/(Tnow*(1._dl+z)) ) / Nnow
            x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )

            call DVERK(Recomb,3,ION,zstart,y,zend,tol,ind,cw,nw,w)
            y(1) = x_H0
            x0 = y(1) + fHe*y(2)
            y(4)=y(3)
        else

            call DVERK(Recomb,nw,ION,zstart,y,zend,tol,ind,cw,nw,w)

            x0 = y(1) + fHe*y(2)

        end if

        Trad = Tnow * (1._dl+zend)
        Tmat = y(3)
        x_H = y(1)
        x_He = y(2)
        x = x0

        zrec(i)=zend
        xrec(i)=x


        if (doTmatTspin) then
            if (Evolve_Ts .and. zend< 1/Do21cm_minev-1 ) then
                Tspin = y(4)
            else
                C10 = Nnow * (1._dl+zend)**3*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) + kappa_eH_21cm(Tmat,.false.)*x)
                tau_21Ts = line21_const*NNow*(1+zend)*dtauda(1/(1+zend))/1000

                Tspin = Trad*( C10/Trad + A10/T_21cm)/(C10/Tmat + A10/T_21cm) + &
                    tau_21Ts/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) )

                y(4) = Tspin
            end if

            tsrec(i) = Tspin
            tmrec(i) = Tmat

        end if

        !          write (*,'(5E15.5)') zend, Trad, Tmat, Tspin, x

    end do

    d0hi=1.0d40
    d0lo=1.0d40
    call spline(zrec,xrec,nz,d0lo,d0hi,dxrec)
    if (doTmatTspin) then
        call spline(zrec,tsrec,nz,d0lo,d0hi,dtsrec)
        call spline(zrec,tmrec,nz,d0lo,d0hi,dtmrec)
    end if
    deallocate(w)

    end subroutine Recombination_init

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



    subroutine ION(Recomb,Ndim,z,Y,f)
    use RECDATA
    implicit none

    integer Ndim
    Type (RecombinationParams) :: Recomb

    real(dl) z,x,n,n_He,Trad,Tmat,Tspin,x_H,x_He, Hz
    real(dl) y(Ndim),f(Ndim)
    real(dl) Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
    real(dl) timeTh,timeH
    real(dl) a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB
    real(dl) tauHe_s,pHe_s
    real(dl) a_trip,b_trip,Rdown_trip,Rup_trip
    real(dl) Doppler,gamma_2Ps,pb,qb,AHcon
    real(dl) tauHe_t,pHe_t,CL_PSt,CfHe_t,gamma_2Pt
    real(dl) epsilon
    integer Heflag
    real(dl) dtauda
    real(dl) C10, dHdz
    external dtauda

    !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    a_PPB = 4.309d0
    b_PPB = -0.6166d0
    c_PPB = 0.6703d0
    d_PPB = 0.5300d0
    !       the Verner and Ferland type fitting parameters for Helium
    !       fixed to match those in the SSS papers, and now correct
    a_VF = 10.d0**(-16.744d0)
    b_VF = 0.711d0
    T_0 = 10.d0**(0.477121d0)   !3K
    T_1 = 10.d0**(5.114d0)
    !      fitting parameters for HeI triplets
    !      (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)

    a_trip = 10.d0**(-16.306d0)
    b_trip = 0.761D0


    x_H = y(1)
    x_He = y(2)
    x = x_H + fHe * x_He
    Tmat = y(3)
    !        Tspin = y(4)

    n = Nnow * (1._dl+z)**3
    n_He = fHe * Nnow * (1._dl+z)**3
    Trad = Tnow * (1._dl+z)

    Hz = 1/dtauda(1/(1._dl+z))*(1._dl+z)**2/MPC_in_sec


    !       Get the radiative rates using PPQ fit, identical to Hummer's table

    Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
        /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)
    Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)

    !       calculate He using a fit to a Verner & Ferland type formula
    sq_0 = sqrt(Tmat/T_0)
    sq_1 = sqrt(Tmat/T_1)
    !       typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
    Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
    Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*exp(-CDB_He/Tmat)
    Rup_He = 4.d0*Rup_He    !statistical weights factor for HeI
    !       Avoid overflow (pointed out by Jacques Roland)
    if((Bfact/Tmat) > 680.d0)then
        He_Boltz = exp(680.d0)
    else
        He_Boltz = exp(Bfact/Tmat)
    end if
    !	now deal with H and its fudges
    if (.not. Recomb%RECFAST_Hswitch) then
        K = CK/Hz !Peebles coefficient K=lambda_a^3/8piH
    else
        !c	fit a double Gaussian correction function
        K = CK/Hz*(1.0d0 &
            +AGauss1*exp(-((log(1.0d0+z)-zGauss1)/wGauss1)**2.d0) &
            +AGauss2*exp(-((log(1.0d0+z)-zGauss2)/wGauss2)**2.d0))
    end if


    !  add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip/(sq_0*(1.d0+sq_0)**(1.0-b_trip))
    Rdown_trip = Rdown_trip/((1.d0+sq_1)**(1.d0+b_trip))
    Rup_trip = Rdown_trip*dexp(-h_P*C*L_He2St_ion/(k_B*Tmat))
    Rup_trip = Rup_trip*((CR*Tmat)**(1.5d0))*(4.d0/3.d0)
    !   last factor here is the statistical weight

    !       try to avoid "NaN" when x_He gets too small
    if ((x_He.lt.5.d-9) .or. (x_He.gt.0.98d0)) then
        Heflag = 0
    else
        Heflag = Recomb%RECFAST_Heswitch
    end if
    if (Heflag.eq.0)then        !use Peebles coeff. for He
        K_He = CK_He/Hz
    else    !for Heflag>0       !use Sobolev escape probability
        tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0-x_He)/Hz
        pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
        K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
        !      if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.99999d0) then
        if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.9999999d0) then
            !AL changed July 08 to get smoother Helium

            !   use fitting formula for continuum opacity of H
            !   first get the Doppler width parameter
            Doppler = 2.D0*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2p*dsqrt(Doppler)
            gamma_2Ps = 3.d0*A2P_s*fHe*(1.d0-x_He)*C*C &
                /(dsqrt(Pi)*sigma_He_2Ps*8.d0*Pi*Doppler*(1.d0-x_H)) &
                /((C*L_He_2p)**2.d0)
            pb = 0.36d0  !value from KIV (2007)
            qb = Recomb%RECFAST_fudge_He
            !   calculate AHcon, the value of A*p_(con,H) for H continuum opacity
            AHcon = A2P_s/(1.d0+pb*(gamma_2Ps**qb))
            K_He=1.d0/((A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
        end if
        if (Heflag.ge.3) then     !include triplet effects
            tauHe_t = A2P_t*n_He*(1.d0-x_He)*3.d0
            tauHe_t = tauHe_t /(8.d0*Pi*Hz*L_He_2Pt**(3.d0))
            pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
            CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
            if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then !Recfast 1.4.2 (?)
                !        if ((Heflag.eq.3) .or. (Heflag.eq.5) .or. x_H >= 0.9999999d0) then    !no H cont. effect
                CfHe_t = A2P_t*pHe_t*dexp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
            else                  !include H cont. effect
                Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
                Doppler = C*L_He_2Pt*dsqrt(Doppler)
                gamma_2Pt = 3.d0*A2P_t*fHe*(1.d0-x_He)*C*C &
                    /(dsqrt(Pi)*sigma_He_2Pt*8.d0*Pi*Doppler*(1.d0-x_H)) &
                    /((C*L_He_2Pt)**2.d0)
                !   use the fitting parameters from KIV (2007) in this case
                pb = 0.66d0
                qb = 0.9d0
                AHcon = A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
                CfHe_t = (A2P_t*pHe_t+AHcon)*dexp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
            end if
        end if
    end if


    !       Estimates of Thomson scattering time and Hubble time
    timeTh=(1._dl/(CT*Trad**4))*(1._dl+x+fHe)/x       !Thomson time
    timeH=2./(3.*HO*(1._dl+z)**1.5)      !Hubble time

    !       calculate the derivatives
    !       turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    !       (clunky, but seems to work)
    if (x_H > 0.99) then   !don't change at all
        f(1) = 0._dl
        !!        else if (x_H > 0.98_dl) then
    else if (x_H.gt.0.985d0) then     !use Saha rate for Hydrogen
        f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat)) /(Hz*(1.d0+z))
        recombination_saha_z = z
        !AL: following commented as not used
        !   for interest, calculate the correction factor compared to Saha
        !   (without the fudge)
        !       factor=(1.d0 + K*Lambda*n*(1.d0-x_H))
        !       /(Hz*(1.d0+z)*(1.d0+K*Lambda*n*(1.d0-x)
        !       +K*Rup*n*(1.d0-x)))
    else !use full rate for H

        f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat)) &
            *(1.d0 + K*Lambda*n*(1.d0-x_H))) &
            /(Hz*(1.d0+z)*(1.d0/fu+K*Lambda*n*(1.d0-x_H)/fu &
            +K*Rup*n*(1.d0-x_H)))

    end if

    !       turn off the He once it is small
    if (x_He < 1.e-15) then
        f(2)=0.d0
    else

        f(2) = ((x*x_He*n*Rdown_He &
            - Rup_He*(1-x_He)*exp(-CL_He/Tmat)) &
            *(1 + K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz)) &
            /(Hz*(1+z) &
            * (1 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))

        !   Modification to HeI recombination including channel via triplets
        if (Heflag.ge.3) then
            f(2) = f(2)+ (x*x_He*n*Rdown_trip &
                - (1.d0-x_He)*3.d0*Rup_trip*dexp(-h_P*C*L_He_2st/(k_B*Tmat))) &
                *CfHe_t/(Hz*(1.d0+z))
        end if

    end if

    if (timeTh < H_frac*timeH) then
        !                f(3)=Tmat/(1._dl+z)      !Tmat follows Trad
        !	additional term to smooth transition to Tmat evolution,
        !	(suggested by Adam Moss)
        dHdz = (HO**2/2.d0/Hz)*(4.d0*(1.d0+z)**3/(1.d0+z_eq)*OmegaT &
            + 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z) )

        epsilon = Hz*(1.d0+x+fHe)/(CT*Trad**3*x)
        f(3) = Tnow &
            + epsilon*((1.d0+fHe)/(1.d0+fHe+x))*((f(1)+fHe*f(2))/x) &
            - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z)

    else
        f(3)= CT * (Trad**4) * x / (1._dl+x+fHe) &
            * (Tmat-Trad) / (Hz*(1._dl+z)) + 2._dl*Tmat/(1._dl+z)
    end if

    ! print *, z, f(3)*(1+z)/Tmat

    if (Do21cm .and. evolve_Ts) then

        !       follow the matter temperature once it has a chance of diverging
        if (timeTh < H_frac*timeH) then
            f(4) = Tnow !spin follows Trad and Tmat
        else
            if (z< 1/Do21cm_minev-1) then

                Tspin = y(4)
                C10 = n*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) + kappa_eH_21cm(Tmat,.false.)*x)

                f(4) = 4*Tspin/Hz/(1+z)*( (Tspin/Tmat-1._dl)*C10 + Trad/T_21cm*(Tspin/Trad-1._dl)*A10) - f(1)*Tspin/(1-x_H)
            else
                f(4)=f(3)
            end if
        end if

    end if

    end subroutine ION



    function dDeltaxe_dtau(a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb)
    !d x_e/d tau assuming Helium all neutral and temperature perturbations negligible
    !it is not accurate for x_e of order 1
    use RECDATA
    implicit none
    real(dl) dDeltaxe_dtau
    real(dl), intent(in):: a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb
    real(dl) Delta_Tg
    real(dl) xedot,z,x,n,n_He,Trad,Tmat,x_H,Hz, C_r, dlnC_r
    real(dl) Rup,Rdown,K
    real(dl) a_PPB,b_PPB,c_PPB,d_PPB
    real(dl) delta_alpha, delta_beta, delta_K, clh
    real(dl) dtauda
    external dtauda


    Delta_tg =Delta_Tm
    x_H = min(1._dl,Recombination_xe(a))

    !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    a_PPB = 4.309d0
    b_PPB = -0.6166d0
    c_PPB = 0.6703d0
    d_PPB = 0.5300d0

    z=1/a-1

    x = x_H

    n = Nnow /a**3
    n_He = fHe * n
    Trad = Tnow /a
    clh = 1/dtauda(a)/a !conformal time
    Hz = clh/a/MPC_in_sec !normal time in seconds

    Tmat = Recombination_tm(a)

    !       Get the radiative rates using PPQ fit, identical to Hummer's table

    Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
        /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)   !alpha
    Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)

    K = CK/Hz              !Peebles coefficient K=lambda_a^3/8piH


    Rdown = Rdown*fu
    Rup = Rup*fu
    C_r =  a*(1.d0 + K*Lambda*n*(1.d0-x_H)) /( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) )*MPC_in_sec

    xedot = -(x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat))*C_r

    delta_alpha = (b_PPB + c_PPB*(Tmat/1d4)**d_PPB*(b_PPB-d_PPB))/(1+c_PPB*(Tmat/1d4)**d_PPB)*Delta_Tg
    delta_beta = delta_alpha + (3./2 + CDB/Tmat)*delta_Tg !(Rup = beta)
    delta_K = - hdot/clh - kvb/clh/3


    dlnC_r = -Rup*K*n*( (Delta_nH+Delta_K + Delta_beta*(1+K*Lambda*n*(1-x_H)))*(1-x_H) - x_H*Delta_xe) &
        / ( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) ) /(1.d0 + K*Lambda*n*(1.d0-x_H))

    dDeltaxe_dtau= xedot/x_H*(dlnC_r +Delta_alpha - Delta_xe) &
        - C_r*( (2*Delta_xe + Delta_nH)*x_H*n*Rdown + (Delta_xe - (3./2+ CB1/Tmat)*(1/x_H-1)*Delta_Tg)*Rup*exp(-CL/Tmat))


    !Approximate form valid at late times
    !        dDeltaxe_dtau= xedot/x_H*(Delta_alpha + Delta_xe + Delta_nH)


    end function dDeltaxe_dtau

    !  ===============================================================


    function polevl(x,coef,N)
    implicit none
    integer N
    real(dl) polevl
    real(dl) x,ans
    real(dl) coef(N+1)

    integer i

    ans=coef(1)
    do i=2,N+1
        ans=ans*x+coef(i)
    end do
    polevl=ans

    end function polevl


    function derivpolevl(x,coef,N)
    implicit none
    integer N
    real(dl) derivpolevl
    real(dl) x,ans
    real(dl) coef(N+1)
    integer i

    ans=coef(1)*N
    do i=2,N
        ans=ans*x+coef(i)*(N-i+1)
    end do
    derivpolevl=ans

    end function derivpolevl


    function kappa_HH_21cm(T, deriv)
    !Polynomail fit to Hydrogen-Hydrogen collision rate as function of Tmatter, from astro-ph/0608032
    !if deriv return d log kappa / d log T
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    !        real(dl), dimension(8), parameter :: fit = &
    !         (/ 0.00120402_dl, -0.0322247_dl,0.339581_dl, -1.75094_dl,4.3528_dl,-4.03562_dl, 1.26899_dl, -29.6113_dl /)
    integer, parameter :: n_table = 27
    integer, dimension(n_table), parameter :: Temps = &
        (/ 1, 2, 4, 6,8,10,15,20,25,30,40,50,60,70,80,90,100,200,300,500,700,1000,2000,3000,5000,7000,10000/)
    real, dimension(n_table), parameter :: rates = &
        (/ 1.38e-13, 1.43e-13,2.71e-13, 6.60e-13,1.47e-12,2.88e-12,9.10e-12,1.78e-11,2.73e-11,&
        3.67e-11,5.38e-11,6.86e-11,8.14e-11,9.25e-11, &
        1.02e-10,1.11e-10,1.19e-10,1.75e-10,2.09e-10,2.56e-10,2.91e-10,3.31e-10,4.27e-10,&
        4.97e-10,6.03e-10,6.87e-10,7.87e-10/)

    real(dl) kappa_HH_21cm, logT, logRate
    real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
    integer xlo, xhi
    real(dl) :: a0, b0, ho

    if (.not. allocated(logRates)) then

        allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
        logRates = log(real(rates,dl)*0.01**3)
        logTemps = log(real(Temps,dl))
        call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
    end if

    if (T<=Temps(1)) then
        if (deriv) then
            kappa_HH_21cm = 0
        else
            kappa_HH_21cm = rates(1)*0.01**3
        end if
        return
    elseif (T >=Temps(n_table)) then
        if (deriv) then
            kappa_HH_21cm = 0
        else
            kappa_HH_21cm = rates(n_table)*0.01**3
        end if
        return
    end if

    logT = log(T)
    xlo=0
    do xhi=2, n_table
        if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
        end  if
    end do
    xhi = xlo+1

    ho=logTemps(xhi)-logTemps(xlo)
    a0=(logTemps(xhi)-logT)/ho
    b0=1-a0

    if (deriv) then
        kappa_HH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
            ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
        !          kappa_HH_21cm = derivpolevl(logT,fit,7)
    else
        logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
        kappa_HH_21cm = exp(logRate)
        !          kappa_HH_21cm = exp(polevl(logT,fit,7))*0.01**3

    end if

    end function kappa_HH_21cm


    function kappa_eH_21cm(T, deriv)
    !Polynomail fit to electron-Hydrogen collision rate as function of Tmatter; from astro-ph/0608032
    !if deriv return d log kappa / d log T
    ! from astro-ph/0608032
    !    1 2.39e-10
    !    2 3.37e-10
    !    5 5.3e-10
    !    10 7.46e-10
    !    20 1.05e-9
    !    50 1.63e-9
    !    100 2.26e-9
    !    200 3.11e-9
    !    500 4.59e-9
    !    1000 5.92e-9
    !    2000 7.15e-9
    !    5000 8.17e-9
    !    10000 8.37e-9
    !    15000 8.29e-9
    !    20000 8.11e-9
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    real(dl), dimension(6), parameter :: fit = &
        (/5.86236d-005,  -0.00171375_dl, 0.0137303_dl, -0.0435277_dl, 0.540905_dl,-22.1596_dl /)

    real(dl) kappa_eH_21cm, logT

    logT = log(T)
    if (deriv) then
        kappa_eH_21cm = derivpolevl(logT,fit,5)
    else
        kappa_eH_21cm = exp(polevl(logT,fit,5))*0.01**3
    end if

    end function kappa_eH_21cm




    function kappa_pH_21cm(T, deriv) ! from astro-ph/0702487
    !Not actually used
    !Polynomail fit to proton-Hydrogen collision rate as function of Tmatter
    !if deriv return d log kappa / d log T
    real(dl), intent(in) :: T
    logical, intent(in) :: deriv
    integer, parameter :: n_table = 17
    integer, dimension(n_table), parameter :: Temps = &
        (/ 1, 2, 5, 10,20,50,100,200,500,1000,2000,3000,5000,7000,10000,15000,20000/)
    real, dimension(n_table), parameter :: rates = &
        (/ 0.4028, 0.4517,0.4301,0.3699,0.3172,0.3047, 0.3379, 0.4043, 0.5471, 0.7051, 0.9167, 1.070, &
        1.301, 1.48,1.695,1.975,2.201/)

    real(dl) kappa_pH_21cm, logT, logRate
    real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
    integer xlo, xhi
    real(dl) :: a0, b0, ho
    real(dl):: factor = 0.01**3*1e-9

    if (.not. allocated(logRates)) then

        allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
        logRates = log(real(rates,dl)*factor)
        logTemps = log(real(Temps,dl))
        call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
    end if

    if (T<=Temps(1)) then
        if (deriv) then
            kappa_pH_21cm = 0
        else
            kappa_pH_21cm = rates(1)*factor
        end if
        return
    elseif (T >=Temps(n_table)) then
        if (deriv) then
            kappa_pH_21cm = 0
        else
            kappa_pH_21cm = rates(n_table)*factor
        end if
        return
    end if

    logT = log(T)
    xlo=0
    do xhi=2, n_table
        if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
        end  if
    end do
    xhi = xlo+1

    ho=logTemps(xhi)-logTemps(xlo)
    a0=(logTemps(xhi)-logT)/ho
    b0=1-a0

    if (deriv) then
        kappa_pH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
            ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
    else
        logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
        kappa_pH_21cm = exp(logRate)
    end if

    end function kappa_pH_21cm


    end module Recombination

