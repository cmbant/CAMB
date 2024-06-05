    ! Modules used by cmbmain and other routines.

    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.
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

    module results
    use constants, only : const_pi, const_twopi
    use MiscUtils
    use RangeUtils
    use StringUtils
    use MathUtils
    use config
    use model
    use splines
    implicit none
    public

    Type TBackgroundOutputs
        real(dl), allocatable :: H(:), DA(:), rs_by_D_v(:)
    end Type TBackgroundOutputs

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13
    integer, parameter :: nthermo_derived = 13

    Type lSamples
        integer :: nl = 0
        integer :: lmin = 2
        integer, allocatable :: l(:)
        logical :: use_spline_template = .true.
    contains
    procedure :: Init => lSamples_init
    procedure :: IndexOf => lSamples_indexOf
    procedure :: InterpolateClArr
    procedure :: InterpolateClArrTemplated
    end Type lSamples

    logical, parameter :: dowinlens = .false. !not used, test getting CMB lensing using visibility
    integer, parameter :: thermal_history_def_timesteps = 20000

    Type TThermoData
        logical :: HasThermoData = .false. !Has it been computed yet for current parameters?
        !Background thermal history, interpolated from precomputed tables
        integer :: nthermo !Number of table steps
        !baryon temperature, sound speed, ionization fractions, and opacity
        real(dl), dimension(:), allocatable :: tb, cs2, xe, dotmu
        ! e^(-tau) and derivatives
        real(dl), dimension(:), allocatable :: emmu, dcs2,demmu, ddotmu, dddotmu, ddddotmu
        real(dl), dimension(:), allocatable :: ScaleFactor, dScaleFactor, adot, dadot
        real(dl), dimension(:), allocatable :: winlens, dwinlens
        real(dl) tauminn,dlntau
        real(dl) :: tight_tau, actual_opt_depth
        !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
        real(dl) :: recombination_saha_tau
        !sound horizon and recombination redshifts
        real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.
        real(dl), dimension(:), allocatable :: step_redshift, rhos_fac, drhos_fac
        real(dl) :: tau_start_redshiftwindows,tau_end_redshiftwindows
        logical :: has_lensing_windows = .false.
        real(dl) recombination_Tgas_tau
        Type(TCubicSpline) :: ScaleFactorAtTime
        !Mapping between redshift and time
        real(dl), private,dimension(:), allocatable :: redshift_time, dredshift_time
        real(dl), private, dimension(:), allocatable :: arhos_fac, darhos_fac, ddarhos_fac
    contains
    procedure :: Init => Thermo_Init
    procedure :: OpacityToTime => Thermo_OpacityToTime
    procedure :: values => Thermo_values
    procedure :: expansion_values => Thermo_expansion_values
    procedure :: IonizationFunctionsAtTime
    procedure, private :: DoWindowSpline
    procedure, private :: SetTimeSteps
    procedure, private :: SetTimeStepWindows
    end type TThermoData

    !Sources
    Type CalWins
        real(dl), allocatable :: awin_lens(:), dawin_lens(:)
    end Type CalWins

    Type LimberRec
        integer n1,n2 !corresponding time step array indices
        real(dl), dimension(:), allocatable :: k
        real(dl), dimension(:), allocatable :: Source
    end Type LimberRec

    Type ClTransferData
        !Cl transfer function variables
        !values of q for integration over q to get C_ls
        Type (lSamples) :: ls ! l that are computed
        integer :: NumSources
        !Changes -scalars:  2 for just CMB, 3 for lensing
        !- tensors: T and E and phi (for lensing), and T, E, B respectively

        type (TRanges) :: q
        real(dl), dimension(:,:,:), allocatable :: Delta_p_l_k

        !The L index of the lowest L to use for Limber
        integer, dimension(:), allocatable :: Limber_l_min
        !For each l, the set of k in each limber window
        !indices LimberWindow(SourceNum,l)
        Type(LimberRec), dimension(:,:), allocatable :: Limber_windows

        !The maximum L needed for non-Limber
        integer max_index_nonlimber

    end Type ClTransferData

    type TCLdata
        Type(ClTransferData) :: CTransScal, CTransTens, CTransVec

        real(dl), dimension (:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , Cl_type)
        !where Cl_type is one of the above constants

        real(dl), dimension (:,:,:), allocatable :: Cl_Scalar_Array
        !Indices are Cl_xxx( l , field1,field2)
        !where ordering of fields is T, E, \psi (CMB lensing potential), window_1, window_2...

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this
        real(dl) , dimension (:,:), allocatable :: Cl_lensed
        !Cl_lensed(l, Cl_type) are the interpolated Cls
    contains
    procedure :: InitCls => TCLdata_InitCls
    procedure :: output_cl_files => TCLdata_output_cl_files
    procedure :: output_lens_pot_files => TCLdata_output_lens_pot_files
    procedure :: NormalizeClsAtL => TCLdata_NormalizeClsAtL
    procedure :: output_veccl_files => TCLdata_output_veccl_files
    end type TCLdata

    Type TTimeSources
        ! values of q to evolve the propagation equations to compute the sources
        type(TRanges) :: Evolve_q
        real(dl), dimension(:,:,:), allocatable :: LinearSrc !Sources and second derivs
        !LinearSrc indices  ( k_index, source_index, time_step_index )
        integer SourceNum, NonCustomSourceNum
        !SourceNum is total number sources (2 or 3 for scalars, 3 for tensors).
    end type TTimeSources

    type, extends(TCAMBdata) :: CAMBdata

        type(CAMBparams) :: CP

        real(dl) ThermoDerivedParams(nthermo_derived)

        logical flat,closed

        !     grhocrit =kappa*a^2*rho_crit(0)
        !     grhornomass=grhor*number of massless neutrino species
        !     taurst,taurend - time at start/end of recombination
        !     dtaurec - dtau during recombination
        !     adotrad - a(tau) in radiation era
        real(dl) grhocrit,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
        real(dl) taurst,dtaurec,taurend,tau_maxvis,adotrad

        real(dl) Omega_de
        real(dl) curv, curvature_radius, Ksign !curvature_radius = 1/sqrt(|curv|), Ksign = 1,0 or -1
        real(dl) tau0,chi0 !time today and rofChi(tau0/curvature_radius)
        real(dl) scale !relative to flat. e.g. for scaling lSamp%l sampling.

        real(dl) akthom !sigma_T * (number density of protons now)
        real(dl) fHe !n_He_tot / n_H_tot
        real(dl) Nnow
        real(dl) z_eq !defined assuming all neutrinos massless
        !Neutrinos
        real(dl) grhormass(max_nu)
        !     nu_masses=m_nu*c**2/(k_B*T_nu0)
        real(dl) nu_masses(max_nu)
        integer ::  num_transfer_redshifts = 1
        real(dl), allocatable  ::  transfer_redshifts(:)
        integer  ::  PK_redshifts_index(max_transfer_redshifts)

        logical :: OnlyTransfer = .false. !C_L/PK not computed; initial power spectrum data, instead get Delta_q_l array
        !If true, sigma_8 is not calculated either]]
        logical :: HasScalarTimeSources = .false. !No power spectra, only time transfer functions

        logical :: get_growth_sigma8 = .true.
        !gets sigma_vdelta, like sigma8 but using velocity-density cross power,
        !in late LCDM f*sigma8 = sigma_vdelta^2/sigma8

        logical :: needs_good_pk_sampling = .false.

        logical ::call_again = .false.
        !if being called again with same parameters to get different thing

        real(dl) :: reion_tau_start, reion_tau_complete
        integer :: reion_n_steps

        Type(TNuPerturbations) :: NuPerturbations

        Type(TBackgroundOutputs) :: BackgroundOutputs

        !Time steps for sampling sources
        Type(TRanges) :: TimeSteps
        !Background interpolation tables for thermal history etc.
        Type(TThermoData) :: ThermoData

        real(dl), allocatable :: transfer_times(:)

        !Matter transfer data
        Type (MatterTransferData):: MT

        !Matter power spectrum for default variable (used for non-linear corrections)
        Type(MatterPowerData), allocatable :: CAMB_PK

        Type(TClData) :: CLdata

        integer :: num_redshiftwindows = 0
        integer :: num_extra_redshiftwindows = 0
        Type(TRedWin), allocatable :: Redshift_W(:)
        real(dl), dimension(:), allocatable :: optical_depths_for21cm

        Type(TTimeSources), allocatable :: ScalarTimeSources
        integer :: Scalar_C_last = C_PhiE


    contains
    procedure :: DeltaTime => CAMBdata_DeltaTime
    procedure :: DeltaTimeArr => CAMBdata_DeltaTimeArr
    procedure :: TimeOfz => CAMBdata_TimeOfz
    procedure :: TimeOfzArr => CAMBdata_TimeOfzArr
    procedure :: DeltaPhysicalTimeGyr => CAMBdata_DeltaPhysicalTimeGyr
    procedure :: DeltaPhysicalTimeGyrArr => CAMBdata_DeltaPhysicalTimeGyrArr
    procedure :: AngularDiameterDistance => CAMBdata_AngularDiameterDistance
    procedure :: AngularDiameterDistanceArr => CAMBdata_AngularDiameterDistanceArr
    procedure :: AngularDiameterDistance2 => CAMBdata_AngularDiameterDistance2
    procedure :: AngularDiameterDistance2Arr => CAMBdata_AngularDiameterDistance2Arr
    procedure :: LuminosityDistance => CAMBdata_LuminosityDistance
    procedure :: ComovingRadialDistance => CAMBdata_ComovingRadialDistance
    procedure :: ComovingRadialDistanceArr => CAMBdata_ComovingRadialDistanceArr
    procedure :: GetBackgroundDensities => CAMBdata_GetBackgroundDensities
    procedure :: Hofz => CAMBdata_Hofz
    procedure :: HofzArr => CAMBdata_HofzArr
    procedure :: sound_horizon => CAMBdata_sound_horizon
    procedure :: sound_horizon_zArr => CAMBdata_sound_horizon_zArr
    procedure :: RedshiftAtTimeArr => CAMBdata_RedshiftAtTimeArr
    procedure :: BAO_D_v => CAMBdata_BAO_D_v
    procedure :: CosmomcTheta => CAMBdata_CosmomcTheta
    procedure :: get_lmax_lensed => CAMBdata_get_lmax_lensed
    procedure :: get_zstar => CAMBdata_get_zstar
    procedure :: DarkEnergyStressEnergy => CAMBdata_DarkEnergyStressEnergy
    procedure :: SetParams => CAMBdata_SetParams
    procedure :: Free => CAMBdata_Free
    procedure :: grho_no_de
    procedure :: GetReionizationOptDepth
    procedure :: rofChi
    procedure :: cosfunc
    procedure :: tanfunc
    procedure :: invsinfunc
    procedure :: GetComputedPKRedshifts
    procedure :: binary_search
    procedure, nopass :: PythonClass => CAMBdata_PythonClass
    procedure, nopass :: SelfPointer => CAMBdata_SelfPointer
    end type CAMBdata

    interface
    FUNCTION state_function(obj, a)
    use precision
    import
    class(CAMBdata) :: obj
    real(dl), intent(in) :: a
    real(dl) :: state_function
    END FUNCTION  state_function
    end interface

    procedure(obj_function), private :: dtauda

    contains

    function CAMBdata_PythonClass()
    character(LEN=:), allocatable :: CAMBdata_PythonClass
    CAMBdata_PythonClass = 'CAMBdata'
    end function CAMBdata_PythonClass

    subroutine CAMBdata_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (CAMBdata), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine CAMBdata_SelfPointer

    subroutine CAMBdata_SetParams(this, P, error, DoReion, call_again, background_only)
    !Initialize background variables; does not yet calculate thermal history
    use constants
    class(CAMBdata), target :: this
    type(CAMBparams), intent(in) :: P
    real(dl) fractional_number, conv
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical, optional :: call_again, background_only
    logical WantReion, calling_again
    integer nu_i,actual_massless
    real(dl) nu_massless_degeneracy, neff_i, eta_k, h2
    real(dl) zpeak, sigma_z, zpeakstart, zpeakend
    Type(TRedWin), pointer :: Win
    logical back_only
    !Constants in SI units

    global_error_flag = 0

    if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
        call GlobalError( 'Cannot generate tensor C_l and transfer without scalar C_l',error_unsupported_params)
    end if

    if (.not. allocated(P%DarkEnergy)) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if

    if (present(error)) error = global_error_flag
    if (global_error_flag/=0) return

    WantReion = DefaultTrue(DoReion)
    calling_again= DefaultFalse(call_again)
    back_only = DefaultFalse(background_only)

    if (calling_again) then
        this%CP%WantDerivedParameters = .false.
        this%CP%Accuracy = P%Accuracy
        this%CP%Transfer%high_precision = P%Transfer%high_precision
        this%CP%Reion%Reionization = P%Reion%Reionization
        this%CP%WantTransfer =P%WantTransfer
        this%CP%WantScalars =P%WantScalars
        this%CP%WantTensors =P%WantTensors
        this%CP%WantVectors =P%WantVectors
        this%CP%WantCls = P%WantCls
    else
        this%CP=P
        this%CP%Max_eta_k = max(this%CP%Max_eta_k,this%CP%Max_eta_k_tensor)
    end if

    if (P%WantTransfer .and. .not. back_only) then
        this%CP%WantScalars=.true.
        if (.not. P%WantCls) then
            this%CP%Accuracy%AccuratePolarization = .false.
            !Sources
            this%CP%Reion%Reionization = this%CP%transfer_21cm_cl
        end if
        call this%GetComputedPKRedshifts(this%CP)
    end if
    if (this%CP%WantTransfer.and. this%CP%MassiveNuMethod==Nu_approx) then
        this%CP%MassiveNuMethod = Nu_trunc
    end if

    if (.not. calling_again) then
        this%ThermoData%HasThermoData = .false.
        if (this%CP%Num_Nu_Massive /= sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))) then
            if (sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))/=0) &
                call GlobalError('Num_Nu_Massive is not sum of Nu_mass_numbers', error_unsupported_params)
        end if
10      if (this%CP%Omnuh2 < 1.e-7_dl) this%CP%Omnuh2 = 0
        if (this%CP%Omnuh2==0 .and. this%CP%Num_Nu_Massive /=0) then
            if (this%CP%share_delta_neff) then
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + this%CP%Num_Nu_Massive
            else
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            end if
            this%CP%Num_Nu_Massive  = 0
            this%CP%Nu_mass_numbers = 0
        end if

        nu_massless_degeneracy = this%CP%Num_Nu_massless !N_eff for massless neutrinos
        if (this%CP%Num_nu_massive > 0) then
            if (this%CP%Nu_mass_eigenstates==0) &
                call GlobalError('Have Num_nu_massive>0 but no nu_mass_eigenstates', error_unsupported_params)
            if (this%CP%Nu_mass_eigenstates==1 .and. this%CP%Nu_mass_numbers(1)==0) &
                this%CP%Nu_mass_numbers(1) = this%CP%Num_Nu_Massive
            if (all(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)==0)) this%CP%Nu_mass_numbers=1 !just assume one for all
            if (this%CP%share_delta_neff) then
                !default case of equal heating of all neutrinos
                fractional_number = this%CP%Num_Nu_massless + this%CP%Num_Nu_massive
                actual_massless = int(this%CP%Num_Nu_massless + 1e-6_dl)
                neff_i = fractional_number/(actual_massless + this%CP%Num_Nu_massive)
                nu_massless_degeneracy = neff_i*actual_massless
                this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates) = &
                    this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)*neff_i
            end if
            if (abs(sum(this%CP%Nu_mass_fractions(1:this%CP%Nu_mass_eigenstates))-1) > 1e-4) &
                call GlobalError('Nu_mass_fractions do not add up to 1', error_unsupported_params)
        else
            this%CP%Nu_mass_eigenstates = 0
        end if

        if (global_error_flag/=0) then
            if (present(error)) error = global_error_flag
            return
        end if


        call This%ThermoData%ScaleFactorAtTime%Clear()

        this%flat = (abs(this%CP%omk) <= OmegaKFlat)
        this%closed = this%CP%omk < -OmegaKFlat

        if (this%flat) then
            this%curv=0
            this%Ksign=0
            this%curvature_radius=1._dl !so we can use tau/curvature_radius, etc, where r's cancel
        else
            this%curv=-this%CP%omk/((c/1000)/this%CP%h0)**2
            this%Ksign =sign(1._dl,this%curv)
            this%curvature_radius=1._dl/sqrt(abs(this%curv))
        end if
        !  grho gives the contribution to the expansion rate from: (g) photons,
        !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
        !  (m) nonrelativistic matter (for Omega=1).  grho is actually
        !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
        !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
        !  (Used only to set the initial conformal time.)

        !H0 is in km/s/Mpc

        this%grhocrit = 3*this%CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

        this%grhog = kappa/c**2*4*sigma_boltz/c**3*this%CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
        ! grhog=1.4952d-13*tcmb**4
        this%grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*this%grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        !grhor=3.3957d-14*tcmb**4

        !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
        !for massive Nu_mass_degeneracies parameters account for heating from grhor

        this%grhornomass=this%grhor*nu_massless_degeneracy
        this%grhormass=0
        do nu_i = 1, this%CP%Nu_mass_eigenstates
            this%grhormass(nu_i)=this%grhor*this%CP%Nu_mass_degeneracies(nu_i)
        end do
        h2 = (this%CP%H0/100)**2
        this%grhoc=this%grhocrit*this%CP%omch2/h2
        this%grhob=this%grhocrit*this%CP%ombh2/h2
        this%grhok=this%grhocrit*this%CP%omk
        this%Omega_de = 1 -(this%CP%omch2 + this%CP%ombh2 + this%CP%omnuh2)/h2 - this%CP%omk  &
            - (this%grhornomass + this%grhog)/this%grhocrit
        this%grhov=this%grhocrit*this%Omega_de

        !  adotrad gives da/dtau in the asymptotic radiation-dominated era:
        this%adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates)))/3)

        this%Nnow = this%CP%ombh2/h2*(1-this%CP%yhe)*this%grhocrit*c**2/kappa/m_H/Mpc**2

        this%akthom = sigma_thomson*this%Nnow*Mpc
        !sigma_T * (number density of protons now)

        this%fHe = this%CP%YHe/(mass_ratio_He_H*(1.d0-this%CP%YHe))  !n_He_tot / n_H_tot

        this%z_eq = (this%grhob+this%grhoc)/(this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates))) -1

        if (this%CP%omnuh2/=0) then
            !Initialize things for massive neutrinos
            call ThermalNuBackground%Init()
            call this%NuPerturbations%Init(P%Accuracy%AccuracyBoost*P%Accuracy%neutrino_q_boost)
            !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0)
            do nu_i=1, this%CP%Nu_mass_eigenstates
                this%nu_masses(nu_i)= ThermalNuBackground%find_nu_mass_for_rho(this%CP%omnuh2/h2*this%CP%Nu_mass_fractions(nu_i)&
                    *this%grhocrit/this%grhormass(nu_i))
            end do
            if (all(this%nu_masses(1:this%CP%Nu_mass_eigenstates)==0)) then
                !All density accounted for by massless, so just use massless
                this%CP%Omnuh2 = 0
                goto 10
            end if
            !Just prevent divide by zero
            this%nu_masses(1:this%CP%Nu_mass_eigenstates) = max(this%nu_masses(1:this%CP%Nu_mass_eigenstates),1e-3_dl)
        else
            this%nu_masses = 0
        end if
        call this%CP%DarkEnergy%Init(this)
        if (global_error_flag==0) this%tau0=this%TimeOfz(0._dl)
        if (global_error_flag==0) then
            this%chi0=this%rofChi(this%tau0/this%curvature_radius)
            this%scale= this%chi0*this%curvature_radius/this%tau0  !e.g. change l sampling depending on approx peak spacing
            if (this%closed .and. this%tau0/this%curvature_radius >3.14) then
                call GlobalError('chi >= pi in closed model not supported',error_unsupported_params)
            end if
            if (WantReion) call this%CP%Reion%Init(this)
            if (this%CP%NonLinear/=NonLinear_None .and. .not. back_only) &
                call this%CP%NonLinearModel%Init(this)
        end if
    end if
    if (allocated(this%CP%SourceWindows) .and. .not. back_only) then
        if (.not. this%CP%WantScalars) then
            this%num_redshiftwindows=0
        else
            this%num_redshiftwindows = size(this%CP%SourceWindows)
        end if
    else
        this%num_redshiftwindows = 0
        this%CP%SourceTerms%limber_windows = .false.
    endif

    if (this%CP%WantScalars .and. this%CP%WantCls .and. this%num_redshiftwindows>0) then
        eta_k = this%CP%Max_eta_k
        if (allocated(this%Redshift_W)) deallocate(this%Redshift_W)
        allocate(this%Redshift_W(this%num_redshiftwindows))
        this%num_extra_redshiftwindows = 0
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(zpeak, sigma_z, zpeakstart, zpeakend, nu_i, Win)
        do nu_i = 1, this%num_redshiftwindows
            Win => this%Redshift_w(nu_i)
            Win%Window => this%CP%SourceWindows(nu_i)%Window
            Win%kind = Win%Window%source_type
            call Win%Window%GetScales(zpeak, sigma_z, zpeakstart, zpeakend)
            if (FeedbackLevel > 1) then
                write(*,*) FormatString('Window scales: %d peak: %f, sigma: %f, start:%f, end %f', &
                    nu_i, zpeak, sigma_z, zpeakstart, zpeakend)
            end if
            Win%Redshift = zpeak
            Win%tau = this%TimeOfz(zpeak, tol=1e-4_dl)
            Win%sigma_tau = sigma_z*dtauda(this,1/(1+zpeak))/(1+zpeak)**2
            Win%tau_peakstart=this%TimeOfZ(zpeakstart, tol=1e-4_dl)
            Win%tau_peakend = this%TimeOfZ(max(0._dl,zpeakend), tol=1e-4_dl)
            Win%chi0 = this%tau0-Win%tau
            Win%chimin = min(Win%chi0,this%tau0 - this%TimeOfz(max(0.05_dl,zpeakend), tol=1e-4_dl))
            !$OMP CRITICAL
            this%CP%Max_eta_k = max(this%CP%Max_eta_k, this%tau0*WindowKmaxForL(Win,this%CP,this%CP%max_l))
            if (Win%Window%source_type==window_21cm) this%CP%Do21cm = .true.
            if (Win%Window%source_type==window_counts .and. P%SourceTerms%counts_lensing) then
                this%num_extra_redshiftwindows = this%num_extra_redshiftwindows + 1
                Win%mag_index = this%num_extra_redshiftwindows
            end if
            !$OMP END CRITICAL
        end do
        if (eta_k /= this%CP%Max_eta_k .and. FeedbackLevel>0) &
            write (*,*) 'source max_eta_k: ', this%CP%Max_eta_k,'kmax = ', this%CP%Max_eta_k/this%tau0
    end if

    if ((this%CP%NonLinear==NonLinear_Lens .or. this%CP%NonLinear==NonLinear_both) .and. &
        this%CP%Max_eta_k/this%tau0 > this%CP%Transfer%kmax) then
        this%CP%Transfer%kmax =this%CP%Max_eta_k/this%tau0
        if (FeedbackLevel > 0) write (*,*) 'kmax changed to ', this%CP%Transfer%kmax
    end if

    if (global_error_flag/=0) then
        if (present(error)) error = global_error_flag
        return
    end if

    if (present(error)) then
        error = 0
    else if (FeedbackLevel > 0 .and. .not. calling_again) then
        write(*,'("Om_b h^2             = ",f9.6)') P%ombh2
        write(*,'("Om_c h^2             = ",f9.6)') P%omch2
        write(*,'("Om_nu h^2            = ",f9.6)') P%omnuh2
        write(*,'("Om_darkenergy        = ",f9.6)') this%Omega_de
        write(*,'("Om_K                 = ",f9.6)') P%omk
        write(*,'("Om_m (inc Om_u)      = ",f9.6)') (P%ombh2+P%omch2+P%omnuh2)/h2
        write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*this%CosmomcTheta()
        if (this%CP%Num_Nu_Massive > 0) then
            write(*,'("N_eff (total)        = ",f9.6)') nu_massless_degeneracy + &
                sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            do nu_i=1, this%CP%Nu_mass_eigenstates
                conv = k_B*(8*this%grhor/this%grhog/7)**0.25*this%CP%tcmb/eV * &
                    (this%CP%nu_mass_degeneracies(nu_i)/this%CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                write(*,'(I2, " nu, g=",f7.4," m_nu*c^2/k_B/T_nu0= ",f9.2," (m_nu= ",f6.3," eV)")') &
                    this%CP%nu_mass_numbers(nu_i), this%CP%nu_mass_degeneracies(nu_i), &
                    this%nu_masses(nu_i),conv*this%nu_masses(nu_i)
            end do
        end if
    end if

    end subroutine CAMBdata_SetParams

    subroutine CAMBdata_Free(this)
    class(CAMBdata) :: this

    call Free_ClTransfer(this%CLdata%CTransScal)
    call Free_ClTransfer(this%ClData%CTransVec)
    call Free_ClTransfer(this%ClData%CTransTens)
    call this%MT%Free()
    if (allocated(this%CAMB_Pk)) deallocate(this%CAMB_PK)

    end subroutine CAMBdata_Free

    function CAMBdata_DeltaTime(this, a1,a2, in_tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol

    atol = PresentDefault(tol/1000/exp(this%CP%Accuracy%AccuracyBoost*this%CP%Accuracy%IntTolBoost-1), in_tol)
    CAMBdata_DeltaTime = Integrate_Romberg(this, dtauda,a1,a2,atol)

    end function CAMBdata_DeltaTime

    subroutine CAMBdata_DeltaTimeArr(this, arr, a1, a2, n, tol)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: a1(n), a2(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%DeltaTime(a1(i), a2(i), tol)
    end do

    end subroutine CAMBdata_DeltaTimeArr

    function CAMBdata_TimeOfz(this, z, tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_TimeOfz
    real(dl), intent(in), optional :: tol
    real(dl), intent(IN) :: z

    CAMBdata_TimeOfz= this%DeltaTime(0._dl,1._dl/(z+1._dl), tol)
    end function CAMBdata_TimeOfz

    subroutine CAMBdata_TimeOfzArr(this, arr, z, n, tol)
    !z array must be monotonically *decreasing* so times increasing
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        if (i==1) then
            arr(i) = this%DeltaTime(0._dl, 1/(1+z(1)), tol)
        else
            if (z(i) < z(i-1)) then
                arr(i) = this%DeltaTime(1/(1+z(i-1)),1/(1+z(i)),tol)
            elseif (z(i) < z(i-1) + 1e-6_dl) then
                arr(i)=0
            else
                error stop 'CAMBdata_TimeOfzArr redshifts must be decreasing'
            end if
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine CAMBdata_TimeOfzArr

    function CAMBdata_DeltaPhysicalTimeGyr(this, a1,a2, in_tol)
    use constants
    class(CAMBdata) :: this
    real(dl), intent(in) :: a1, a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) CAMBdata_DeltaPhysicalTimeGyr, atol

    atol = PresentDefault(1d-4/exp(this%CP%Accuracy%AccuracyBoost-1), in_tol)
    CAMBdata_DeltaPhysicalTimeGyr = Integrate_Romberg(this, dtda,a1,a2,atol)*Mpc/c/Gyr
    end function CAMBdata_DeltaPhysicalTimeGyr

    subroutine CAMBdata_DeltaPhysicalTimeGyrArr(this, arr, a1, a2, n, tol)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: a1(n), a2(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%DeltaPhysicalTimeGyr(a1(i), a2(i), tol)
    end do

    end subroutine CAMBdata_DeltaPhysicalTimeGyrArr


    function CAMBdata_AngularDiameterDistance(this,z)
    class(CAMBdata) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc
    real(dl) CAMBdata_AngularDiameterDistance
    real(dl), intent(in) :: z

    CAMBdata_AngularDiameterDistance = this%curvature_radius/(1+z)* &
        this%rofchi(this%ComovingRadialDistance(z) /this%curvature_radius)

    end function CAMBdata_AngularDiameterDistance

    subroutine CAMBdata_AngularDiameterDistanceArr(this, arr, z, n)
    class(CAMBdata) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc for array of z
    !z array must be monotonically increasing
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    call this%ComovingRadialDistanceArr(arr, z, n, 1e-4_dl)
    if (this%flat) then
        arr = arr/(1+z)
    else
        do i=1, n
            arr(i) =  this%curvature_radius/(1+z(i))*this%rofchi(arr(i)/this%curvature_radius)
        end do
    end if

    end subroutine CAMBdata_AngularDiameterDistanceArr


    function CAMBdata_AngularDiameterDistance2(this,z1, z2)
    ! z1 < z2, otherwise returns zero
    !From http://www.slac.stanford.edu/~amantz/work/fgas14/#cosmomc
    class(CAMBdata) :: this
    real(dl) CAMBdata_AngularDiameterDistance2
    real(dl), intent(in) :: z1, z2

    if (z2 < z1 + 1e-4) then
        CAMBdata_AngularDiameterDistance2=0
    else
        CAMBdata_AngularDiameterDistance2 = this%curvature_radius/(1+z2)* &
            this%rofchi( this%DeltaTime(1/(1+z2),1/(1+z1))/this%curvature_radius)
    end if

    end function CAMBdata_AngularDiameterDistance2

    subroutine CAMBdata_AngularDiameterDistance2Arr(this, arr, z1, z2, n)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z1(n), z2(n)
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%AngularDiameterDistance2(z1(i),z2(i))
    end do
    !$OMP END PARALLEL DO

    end subroutine CAMBdata_AngularDiameterDistance2Arr


    function CAMBdata_LuminosityDistance(this,z)
    class(CAMBdata) :: this
    real(dl) CAMBdata_LuminosityDistance
    real(dl), intent(in) :: z

    CAMBdata_LuminosityDistance = this%AngularDiameterDistance(z)*(1+z)**2

    end function CAMBdata_LuminosityDistance

    function CAMBdata_ComovingRadialDistance(this, z)
    class(CAMBdata) :: this
    real(dl) CAMBdata_ComovingRadialDistance
    real(dl), intent(in) :: z

    CAMBdata_ComovingRadialDistance = this%DeltaTime(1/(1+z),1._dl)

    end function CAMBdata_ComovingRadialDistance

    subroutine CAMBdata_ComovingRadialDistanceArr(this, arr, z, n, tol)
    !z array must be monotonically increasing
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    real(dl), intent(in) :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        if (i==1) then
            if (z(i) < 1e-6_dl) then
                arr(i) = 0
            else
                arr(i) = this%DeltaTime(1/(1+z(i)),1._dl, tol)
            end if
        else
            if (z(i) < z(i-1)) error stop 'ComovingRadialDistanceArr redshifts out of order'
            arr(i) = this%DeltaTime(1/(1+z(i)),1/(1+z(i-1)),tol)
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine CAMBdata_ComovingRadialDistanceArr

    function CAMBdata_Hofz(this,z)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    class(CAMBdata) :: this
    real(dl) CAMBdata_Hofz, a
    real(dl), intent(in) :: z

    a = 1/(1+z)
    CAMBdata_Hofz = 1/(a**2*dtauda(this,a))

    end function CAMBdata_Hofz

    subroutine CAMBdata_HofzArr(this, arr, z, n)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i
    real(dl) :: a

    do i=1, n
        a = 1/(1+z(i))
        arr(i) = 1/(a**2*dtauda(this,a))
    end do

    end subroutine CAMBdata_HofzArr

    real(dl) function CAMBdata_sound_horizon(this, z)
    class(CAMBdata) :: this
    real(dl), intent(in) :: z

    CAMBdata_sound_horizon = Integrate_Romberg(this,dsound_da_exact,1d-9,1/(z+1),1e-6_dl)

    end function CAMBdata_sound_horizon

    subroutine CAMBdata_sound_horizon_zArr(this,arr, z,n)
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), IF(n>4)
    do i=1,n
        arr(i) = this%sound_horizon(z(i))
    end do

    end subroutine CAMBdata_sound_horizon_zArr

    subroutine CAMBdata_RedshiftAtTimeArr(this, arr, tau, n)
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: tau(n)
    integer i
    real(dl) om

    if (this%ThermoData%ScaleFactorAtTime%n==0) &
        call GlobalError('RedshiftAtTimeArr: background history not calculated', error_unsupported_params)
    if (global_error_flag/=0) return
    !$OMP PARALLEL DO DEFAULT(SHARED), private(om, i)
    do i=1, n
        if (tau(i) < this%ThermoData%tauminn*1.1) then
            om = (this%grhob+this%grhoc)/&
                sqrt(3*(this%grhog+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates))+this%grhornomass))
            arr(i) = 1/(this%adotrad*tau(i)*(1+om*tau(i)/4))-1
        else
            arr(i) = 1/this%ThermoData%ScaleFactorAtTime%Value(tau(i))-1
        end if
    end do

    end subroutine CAMBdata_RedshiftAtTimeArr

    real(dl) function BAO_D_v_from_DA_H(z, DA, Hz)
    real(dl), intent(in) :: z, DA, Hz
    real(dl) ADD

    ADD = DA*(1.d0+z)
    BAO_D_v_from_DA_H = ((ADD)**2.d0*z/Hz)**(1.d0/3.d0)

    end function BAO_D_v_from_DA_H

    real(dl) function CAMBdata_BAO_D_v(this,z)
    class(CAMBdata) :: this
    real(dl), intent(IN) :: z

    CAMBdata_BAO_D_v = BAO_D_v_from_DA_H(z,this%AngularDiameterDistance(z), this%Hofz(z))

    end function CAMBdata_BAO_D_v

    function CAMBdata_get_zstar(this)
    class(CAMBdata) :: this
    real(dl) CAMBdata_get_zstar
    real(dl) z_scale

    call this%CP%Recomb%Init(this)
    z_scale =  COBE_CMBTemp/this%CP%TCMB
    CAMBdata_get_zstar=this%binary_search(noreion_optdepth, 1.d0, 700.d0*z_scale, &
        2000.d0*z_scale, 1d-3*z_scale,100.d0*z_scale,5000.d0*z_scale)

    end function CAMBdata_get_zstar

    function CAMBdata_CosmomcTheta(this)
    class(CAMBdata) :: this
    real(dl) zstar, astar, atol, rs, DA
    real(dl) CAMBdata_CosmomcTheta
    real(dl) ombh2, omdmh2

    ombh2 = this%CP%ombh2
    omdmh2 = (this%CP%omch2+this%CP%omnuh2)

    !!From Hu & Sugiyama
    zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
        (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
        (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

    astar = 1/(1+zstar)
    atol = 1e-6
    rs = Integrate_Romberg(this,dsound_da_approx,1d-8,astar,atol)
    DA = this%AngularDiameterDistance(zstar)/astar
    CAMBdata_CosmomcTheta = rs/DA

    end function CAMBdata_CosmomcTheta


    subroutine CAMBdata_GetBackgroundDensities(this, n, a_arr, densities)
    ! return array of 8*pi*G*rho*a**4 for each species
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a_arr(n)
    real(dl) :: grhov_t, rhonu, grhonu, a
    real(dl), intent(out) :: densities(8,n)
    integer nu_i,i

    do i=1, n
        a = a_arr(i)
        call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)
        grhonu = 0

        if (this%CP%Num_Nu_massive /= 0) then
            !Get massive neutrino density relative to massless
            do nu_i = 1, this%CP%nu_mass_eigenstates
                call ThermalNuBackground%rho(a * this%nu_masses(nu_i), rhonu)
                grhonu = grhonu + rhonu * this%grhormass(nu_i)
            end do
        end if

        densities(2,i) = this%grhok * a**2
        densities(3,i) = this%grhoc * a
        densities(4,i) = this%grhob * a
        densities(5,i) = this%grhog
        densities(6,i) = this%grhornomass
        densities(7,i) = grhonu
        densities(8,i) = grhov_t*a**2
        densities(1,i) = sum(densities(2:8,i))
    end do

    end subroutine CAMBdata_GetBackgroundDensities

    integer function CAMBdata_get_lmax_lensed(this)
    class(CAMBdata) :: this
    CAMBdata_get_lmax_lensed = this%CLdata%lmax_lensed
    end function CAMBdata_get_lmax_lensed

    !JD 08/13 New function for nonlinear lensing of CMB + MPK compatibility
    !Build master redshift array from array of desired Nonlinear lensing (NLL)
    !redshifts and an array of desired Power spectrum (PK) redshifts.
    !At the same time fill arrays for NLL and PK that indicate indices
    !of their desired redshifts in the master redshift array.
    !Finally define number of redshifts in master array. This is usually given by:
    !P%num_redshifts = P%PK_num_redshifts + NLL_num_redshifts - 1.  The -1 comes
    !from the fact that z=0 is in both arrays (when non-linear is on)
    subroutine GetComputedPKRedshifts(this, Params,eta_k_max)
    use MpiUtils, only : MpiStop
    class(CAMBdata) :: this
    Type(CAMBParams) :: Params
    integer i, iPK, iNLL
    real(dl), parameter :: tol = 1.d-5
    real(dl) maxRedshift, NL_Boost
    integer   ::  NLL_num_redshifts
    real(dl), allocatable    ::  NLL_redshifts(:), redshifts(:)
    !Sources, but unused currently
    real(dl), intent(in), optional :: eta_k_max

    NLL_num_redshifts = 0
    this%needs_good_pk_sampling = .false.
    associate(P => Params%Transfer)
        if ((Params%NonLinear==NonLinear_lens .or. Params%NonLinear==NonLinear_both) .and. &
            (Params%DoLensing .or. this%num_redshiftwindows > 0)) then
            ! Want non-linear lensing or other sources
            this%needs_good_pk_sampling = .false.
            NL_Boost = Params%Accuracy%AccuracyBoost*Params%Accuracy%NonlinSourceBoost
            if (Params%Do21cm) then
                !Sources
                if (maxval(this%Redshift_w(1:this%num_redshiftwindows)%Redshift) &
                    /= minval(this%Redshift_w(1:this%num_redshiftwindows)%Redshift))  &
                    stop 'Non-linear 21cm currently only for narrow window at one redshift'
                if (.not. present(eta_k_max)) stop 'bad call to GetComputedPKRedshifts'
                P%kmax = eta_k_max/10000.
                NLL_num_redshifts =  1
                allocate(NLL_redshifts(NLL_num_redshifts+1))
                NLL_redshifts(1) = this%Redshift_w(1)%Redshift
            else
                P%kmax = max(P%kmax,5*NL_Boost)
                maxRedshift = 10
                NLL_num_redshifts =  nint(10*5*NL_Boost)
                if (NL_Boost>=2.5) then
                    !only notionally more accuracy, more stable for RS
                    maxRedshift =15
                end if
                allocate(NLL_redshifts(NLL_num_redshifts+1)) !+1 to stop access issues below
                do i=1,NLL_num_redshifts
                    NLL_redshifts(i) = real(NLL_num_redshifts-i)/(NLL_num_redshifts/maxRedshift)
                end do
            end if
        end if
        if (allocated(this%transfer_redshifts)) deallocate(this%transfer_redshifts)
        if (NLL_num_redshifts==0) then
            this%num_transfer_redshifts=P%PK_num_redshifts
            allocate(this%transfer_redshifts(this%num_transfer_redshifts))
            this%transfer_redshifts = P%PK_redshifts(:this%num_transfer_redshifts)
            this%PK_redshifts_index(:this%num_transfer_redshifts) = (/ (i, i=1, this%num_transfer_redshifts ) /)
        else
            i=0
            iPK=1
            iNLL=1
            allocate(redshifts(NLL_num_redshifts+P%PK_num_redshifts))
            do while (iPk<=P%PK_num_redshifts .or. iNLL<=NLL_num_redshifts)
                !JD write the next line like this to account for roundoff issues with ==. Preference given to PK_Redshift
                i=i+1
                if(iNLL>NLL_num_redshifts .or. P%PK_redshifts(iPK)>NLL_redshifts(iNLL)+tol) then
                    redshifts(i)=P%PK_redshifts(iPK)
                    this%PK_redshifts_index(iPK)=i
                    iPK=iPK+1
                else if(iPK>P%PK_num_redshifts .or. NLL_redshifts(iNLL)>P%PK_redshifts(iPK)+tol) then
                    redshifts(i)=NLL_redshifts(iNLL)
                    iNLL=iNLL+1
                else
                    redshifts(i)=P%PK_redshifts(iPK)
                    this%PK_redshifts_index(iPK)=i
                    iPK=iPK+1
                    iNLL=iNLL+1
                end if
            end do
            this%num_transfer_redshifts=i
            allocate(this%transfer_redshifts(this%num_transfer_redshifts))
            this%transfer_redshifts = redshifts(:this%num_transfer_redshifts)
        end if
    end associate

    end subroutine GetComputedPKRedshifts

    subroutine CAMBdata_DarkEnergyStressEnergy(this, a, grhov_t, w, n)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n)
    real(dl), intent(out) :: grhov_t(n), w(n)
    integer i

    do i=1, n
        call this%CP%DarkEnergy%BackgroundDensityAndPressure(1._dl, a(i), grhov_t(i), w(i))
    end do
    grhov_t = grhov_t/a**2

    end subroutine CAMBdata_DarkEnergyStressEnergy

    function rofChi(this,Chi) !sinh(chi) for open, sin(chi) for closed.
    class(CAMBdata) :: this
    real(dl) Chi,rofChi

    if (this%flat) then
        rofChi=chi
    else if (this%closed) then
        rofChi=sin(chi)
    else
        rofChi=sinh(chi)
    endif
    end function rofChi


    function cosfunc (this,Chi)
    class(CAMBdata) :: this
    real(dl) Chi,cosfunc

    if (this%flat) then
        cosfunc = 1._dl
    else if (this%closed) then
        cosfunc= cos(chi)
    else
        cosfunc=cosh(chi)
    endif
    end function cosfunc

    function tanfunc(this,Chi)
    class(CAMBdata) :: this
    real(dl) Chi,tanfunc
    if (this%flat) then
        tanfunc=Chi
    else if (this%closed) then
        tanfunc=tan(Chi)
    else
        tanfunc=tanh(Chi)
    end if

    end function tanfunc

    function invsinfunc(this,x)
    class(CAMBdata) :: this
    real(dl) invsinfunc,x

    if (this%flat) then
        invsinfunc = x
    else if (this%closed) then
        invsinfunc=asin(x)
    else
        invsinfunc=log((x+sqrt(1._dl+x**2)))
    endif
    end function invsinfunc

    function dsound_da_exact(this,a)
    class(CAMBdata) :: this
    real(dl) dsound_da_exact,a,R,cs

    R = 3*this%grhob*a / (4*this%grhog)
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_exact=dtauda(this,a)*cs

    end function dsound_da_exact

    function dsound_da_approx(this,a)
    !approximate form used e.g. by CosmoMC for theta
    class(CAMBdata) :: this
    real(dl) dsound_da_approx,a,R,cs

    R=3.0d4*a*this%CP%ombh2
    !          R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and used for previous consistency
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_approx=dtauda(this,a)*cs

    end function dsound_da_approx

    function dtda(this,a)
    class(CAMBdata) :: this
    real(dl) dtda,a

    dtda= dtauda(this,a)*a
    end function

    function ddamping_da(this, a)
    class(CAMBdata) :: this
    real(dl) :: ddamping_da
    real(dl), intent(in) :: a
    real(dl) :: R

    R=this%ThermoData%r_drag0*a
    !ignoring reionisation, not relevant for distance measures
    ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(this,a)*a**2/(this%CP%Recomb%x_e(a)*this%akthom)

    end function ddamping_da

    function noreion_doptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: noreion_doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)

    !ignoring reionisation, not relevant for distance measures
    noreion_doptdepth_dz = this%CP%Recomb%x_e(a)*this%akthom*dtauda(this,a)

    end function noreion_doptdepth_dz

    function noreion_optdepth(this,z)
    class(CAMBdata) :: this
    real(dl) noreion_optdepth
    real(dl),intent(in) :: z

    noreion_optdepth = Integrate_Romberg(this, noreion_doptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function noreion_optdepth

    function ddragoptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: ddragoptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)
    ddragoptdepth_dz = noreion_doptdepth_dz(this,z)/this%ThermoData%r_drag0/a

    end function ddragoptdepth_dz

    function dragoptdepth(this,z)
    class(CAMBdata) :: this
    real(dl) dragoptdepth
    real(dl),intent(in) :: z

    dragoptdepth =  Integrate_Romberg(this,ddragoptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function dragoptdepth

    function reion_doptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: reion_doptdepth_dz
    real(dl), intent(in) :: z

    reion_doptdepth_dz = this%CP%Reion%x_e(z)*this%akthom*dtauda(this,1._dl/(1._dl+z))

    end function reion_doptdepth_dz

    function grho_no_de(this, a) result(grhoa2)
    !  Return 8*pi*G*rho_no_de*a**4 where rho_no_de includes everything except dark energy.
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) grhoa2, rhonu
    integer nu_i

    grhoa2 = this%grhok * a**2 + (this%grhoc + this%grhob) * a + this%grhog + this%grhornomass

    if (this%CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, this%CP%nu_mass_eigenstates
            call ThermalNuBack%rho(a * this%nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * this%grhormass(nu_i)
        end do
    end if

    end function grho_no_de

    function GetReionizationOptDepth(this)
    class(CAMBdata) :: this
    real(dl) GetReionizationOptDepth
    integer n
    real(dl) zstart, zend

    call this%CP%Reion%get_timesteps(n, zstart, zend)
    GetReionizationOptDepth = Integrate_Romberg(this, reion_doptdepth_dz,0.d0,zstart,&
        1d-5/this%CP%Accuracy%AccuracyBoost)

    end function GetReionizationOptDepth

    real(dl) function binary_search(this,func, goal, x1, x2, tol, widex1, widex2)
    !This is about twice as inefficient as Brent
    class(CAMBdata) :: this
    procedure(state_function) :: func
    real(dl), intent(in) :: goal,x1,x2,tol
    real(dl), intent(in), optional :: widex1, widex2 !Wider range in case of failure
    real(dl) try_t, try_b, avg, D_try, last_bot, last_top, diff
    integer count
    logical wide

    try_b = x1
    try_t = x2
    diff = tol*2
    count = 0
    wide = .false.
    do while (diff > tol)
        if (count>100) then
            if (.not. wide .and. present(widex1)) then
                count=0
                wide=.true.
                try_b=widex1
                try_t=widex2
            else
                call GlobalError(FormatString('binary_search (e.g for optical depth) did not converge: ' //&
                    'Base range %f-%f.',x1,x2),error_reionization)
                binary_search = 0
                return
            end if
        end if
        avg = (try_b+try_t)/2
        D_try = func(this,avg)
        count = count+1
        if (D_try < goal) then
            try_b = avg
            last_bot = D_try
        else
            try_t = avg
            last_top = D_try
        end if
        diff = abs(D_try - goal)
    end do
    if (try_b==x1) last_bot = func(this,x1)
    if (try_t==x2) last_top = func(this,x2)
    binary_search =  (try_t*(goal-last_bot) + try_b*(last_top-goal))/(last_top-last_bot)

    end function binary_search

    !Sources
    function WindowKmaxForL(W,CP, ell) result(res)
    class(CAMBparams), intent(in) :: CP
    Type(TRedWin), intent(in) :: W
    real(dl) res
    integer, intent(in)::  ell

    if (W%kind == window_lensing) then
        res = CP%Accuracy%AccuracyBoost*18*ell/W%chi0
    else
        !On large scales power can be aliased from smaller, so make sure k goes up until at least the turnover
        !in the matter power spectrum
        res = CP%Accuracy%AccuracyBoost*max(0.05_dl,2.5*ell/W%chimin)
    end if

    res = res* CP%Accuracy%KmaxBoost
    end function WindowKmaxForL


    function lSamples_indexOf(lSet,l)
    class(lSamples) :: lSet
    integer, intent(in) :: l
    integer lSamples_indexOf, i

    do i=2,lSet%nl
        if (l < lSet%l(i)) then
            lSamples_indexOf = i-1
            return
        end if
    end do
    lSamples_indexOf = lSet%nl

    end function  lSamples_indexOf

    subroutine lSamples_init(this, State, lmin, max_l)
    ! This subroutines initializes lSet%l arrays. Other values will be interpolated.
    class(lSamples) :: this
    class(CAMBdata), target :: State
    integer, intent(IN) :: lmin,max_l
    integer lind, lvar, step, top, bot, lmin_log
    integer, allocatable :: ls(:)
    real(dl) AScale

    allocate(ls(max_l))
    if (allocated(this%l)) deallocate(this%l)
    this%lmin = lmin
    this%use_spline_template = State%CP%use_cl_spline_template
    lmin_log = State%CP%min_l_logl_sampling
    associate(Accuracy => State%CP%Accuracy)
        Ascale=State%scale/Accuracy%lSampleBoost

        if (Accuracy%lSampleBoost >=50) then
            !just do all of them
            lind=0
            do lvar=lmin, max_l
                lind=lind+1
                ls(lind)=lvar
            end do
            this%nl=lind
            allocate(this%l(lind))
            this%l = ls(1:lind)
            return
        end if

        lind=0
        do lvar=lmin, 10
            lind=lind+1
            ls(lind)=lvar
        end do

        if (Accuracy%AccurateReionization) then
            do lvar=11, 14
                lind=lind+1
                ls(lind)=lvar
            end do
            if (Accuracy%lSampleBoost > 1) then
                do lvar=15, 37, 1
                    lind=lind+1
                    ls(lind)=lvar
                end do
            else
                do lvar=15, 37, 2
                    lind=lind+1
                    ls(lind)=lvar
                end do
            end if

            step = max(nint(5*Ascale),2)
            bot=40
            top=bot + step*10
        else
            if (Accuracy%lSampleBoost >1) then
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

        if (State%CP%Log_lvalues) then
            !Useful for generating smooth things like 21cm to high l
            step=max(nint(20*Ascale),4)
            do
                lvar = lvar + step
                if (lvar > max_l) exit
                lind=lind+1
                ls(lind)=lvar
                step = nint(step*1.2) !log spacing
            end do
        else
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
                    if (.not. this%use_spline_template) then
                        step=max(nint(42*Ascale),7)
                    else
                        step=max(nint(50*Ascale),7)
                    end if
                    bot=ls(lind)+step
                    top=min(lmin_log,max_l)

                    do lvar = bot,top,step
                        lind=lind+1
                        ls(lind)=lvar
                    end do

                    if (max_l > lmin_log) then
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
                        if (ls(lind) < max_l - 100) then
                            !Try to keep lensed spectra up to specified lmax
                            lind=lind+1
                            ls(lind)=max_l - lensed_convolution_margin
                        else if (ls(lind) - ls(lind-1) > lensed_convolution_margin) then
                            ls(lind)=max_l - lensed_convolution_margin
                        end if
                    end if
                end if !log_lvalues
                if (ls(lind) /=max_l) then
                    lind=lind+1
                    ls(lind)=max_l
                end if
                if (.not. State%flat .and. max_l<=lmin_log) ls(lind-1)=int(max_l+ls(lind-2))/2
                !Not in flat case so interpolation table is the same when using lower l_max
            end if
        end if
    end associate
    this%nl=lind
    allocate(this%l, source=ls(1:lind))

    end subroutine lSamples_init

    subroutine InterpolateClArr(lSet,iCl, all_Cl, max_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(1:*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in), optional :: max_index
    integer il,llo,lhi, xi
    real(dl) ddCl(lSet%nl)
    real(dl) xl(lSet%nl)
    real(dl) a0,b0,ho
    integer max_ind

    max_ind = PresentDefault(lSet%nl, max_index)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArr')

    xl = real(lSet%l(1:lSet%nl),dl)
    call spline_def(xl,iCL,max_ind,ddCl)

    llo=1
    do il=lSet%lmin,lSet%l(max_ind)
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

    subroutine InterpolateClArrTemplated(lSet,iCl, all_Cl, max_ind, template_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in) :: max_ind
    integer, intent(in), optional :: template_index
    integer maxdelta, il
    real(dl) DeltaCL(lSet%nl)
    real(dl), allocatable :: tmpall(:)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArrTemplated')
    if (lSet%use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
            !interpolate only the difference between the C_l and an accurately interpolated template.
            !Using unlensed for template, seems to be good enough
            maxdelta=max_ind
            do while (lSet%l(maxdelta) > lmax_extrap_highl)
                maxdelta=maxdelta-1
            end do
            DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

            call lSet%InterpolateClArr(DeltaCl, all_Cl, maxdelta)

            do il=lSet%lmin,lSet%l(maxdelta)
                all_Cl(il) = all_Cl(il) +  highL_CL_template(il,template_index)
            end do

            if (maxdelta < max_ind) then
                !directly interpolate high L where no t  emplate (doesn't effect lensing spectrum much anyway)
                allocate(tmpall(lSet%lmin:lSet%l(max_ind)))
                call InterpolateClArr(lSet,iCl, tmpall, max_ind)
                !overlap to reduce interpolation artefacts
                all_cl(lSet%l(maxdelta-2):lSet%l(max_ind) ) = tmpall(lSet%l(maxdelta-2):lSet%l(max_ind))
                deallocate(tmpall)
            end if
            return
        end if
    end if

    call InterpolateClArr(lSet,iCl, all_Cl, max_ind)

    end subroutine InterpolateClArrTemplated


    subroutine Thermo_values(this,tau, a, cs2b, opacity, dopacity)
    !Compute unperturbed sound speed squared,
    !and ionization fraction by interpolating pre-computed tables.
    !If requested also get time derivative of opacity
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out) :: a, cs2b, opacity
    real(dl), intent(out), optional :: dopacity
    integer i
    real(dl) d

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
        write(*,*) 'tau, taumin = ', tau, this%tauminn
        call MpiStop('thermo out of bounds')
    else if (i >= this%nthermo) then
        cs2b=this%cs2(this%nthermo)
        opacity=this%dotmu(this%nthermo)
        a=1
        if (present(dopacity)) then
            dopacity = this%ddotmu(this%nthermo)/(tau*this%dlntau)
        end if
    else
        cs2b=this%cs2(i)+d*(this%dcs2(i)+d*(3*(this%cs2(i+1)-this%cs2(i))  &
            -2*this%dcs2(i)-this%dcs2(i+1)+d*(this%dcs2(i)+this%dcs2(i+1)  &
            +2*(this%cs2(i)-this%cs2(i+1)))))
        opacity=this%dotmu(i)+d*(this%ddotmu(i)+d*(3*(this%dotmu(i+1)-this%dotmu(i)) &
            -2*this%ddotmu(i)-this%ddotmu(i+1)+d*(this%ddotmu(i)+this%ddotmu(i+1) &
            +2*(this%dotmu(i)-this%dotmu(i+1)))))
        a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
            -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
            +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
        if (present(dopacity)) then
            dopacity=(this%ddotmu(i)+d*(this%dddotmu(i)+d*(3*(this%ddotmu(i+1)  &
                -this%ddotmu(i))-2*this%dddotmu(i)-this%dddotmu(i+1)+d*(this%dddotmu(i) &
                +this%dddotmu(i+1)+2*(this%ddotmu(i)-this%ddotmu(i+1))))))/(tau*this%dlntau)
        end if
    end if
    end subroutine Thermo_values

    subroutine Thermo_expansion_values(this, tau, a, adot, opacity)
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out) :: a, adot, opacity
    integer i
    real(dl) d

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
        call MpiStop('thermo out of bounds')
    else if (i >= this%nthermo) then
        opacity=this%dotmu(this%nthermo)
        a=1
        adot=this%adot(this%nthermo)
    else
        a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
            -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
            +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau

        adot = (this%adot(i)+d*(this%dadot(i)+d*(3*(this%adot(i+1)-this%adot(i)) &
            -2*this%dadot(i)-this%dadot(i+1)+d*(this%dadot(i)+this%dadot(i+1) &
            +2*(this%adot(i)-this%adot(i+1))))))

        opacity=this%dotmu(i)+d*(this%ddotmu(i)+d*(3*(this%dotmu(i+1)-this%dotmu(i)) &
            -2*this%ddotmu(i)-this%ddotmu(i+1)+d*(this%ddotmu(i)+this%ddotmu(i+1) &
            +2*(this%dotmu(i)-this%dotmu(i+1)))))
    end if

    end subroutine Thermo_expansion_values

    function Thermo_OpacityToTime(this,opacity)
    class(TThermoData) :: this
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(this%dotmu(j)> opacity)
        j=j+1
    end do

    Thermo_OpacityToTime = exp((j-1)*this%dlntau)*this%tauminn

    end function Thermo_OpacityToTime

    subroutine Thermo_Init(this, State,taumin)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use StringUtils
    class(TThermoData) :: this
    class(CAMBdata), target :: State
    real(dl), intent(in) :: taumin
    integer nthermo
    real(dl) tau01,a0,barssc,dtau
    real(dl) tau,a,a2
    real(dl) adot,fe,thomc0
    real(dl) dtbdla,vfi,cf1,maxvis, vis, z_scale
    integer ncount,i,j1,iv,ns
    real(dl), allocatable :: spline_data(:)
    real(dl) last_dotmu, om
    real(dl) a_verydom
    real(dl) awin_lens1p,awin_lens2p,dwing_lens, rs, DA
    real(dl) a_eq, rs_eq, tau_eq, rstar
    integer noutput
    Type(CalWins), dimension(:), allocatable, target :: RW
    real(dl) awin_lens1(State%num_redshiftwindows),awin_lens2(State%num_redshiftwindows)
    real(dl) Tspin, Trad, rho_fac, window, tau_eps
    integer transfer_ix(State%CP%Transfer%PK_num_redshifts)
    integer RW_i, j2
    real(dl) Tb21cm, winamp, z, background_boost
    character(len=:), allocatable :: outstr
    real(dl), allocatable ::  taus(:)
    real(dl), allocatable :: xe_a(:), sdotmu(:), opts(:)
    real(dl), allocatable :: scale_factors(:), times(:), dt(:)
    Type(TCubicSpline) :: dotmuSp
    integer ninverse, nlin
    real(dl) dlna, zstar_min, zstar_max
    real(dl) reion_z_start, reion_z_complete
    Type(CAMBParams), pointer :: CP

    CP => State%CP

    !Allocate memory outside parallel region to keep ifort happy
    background_boost = CP%Accuracy%BackgroundTimeStepBoost*CP%Accuracy%AccuracyBoost
    if (background_boost > 20) then
        write(*,*) 'Warning: very small time steps can give less accurate spline derivatives'
        write(*,*) 'e.g. around reionization if not matched very smoothly'
    end if
    !Higher q starts earlier; scale by log(taumin) so actual step size is not made worse by increasing k_max
    nthermo = nint(thermal_history_def_timesteps*log(1.4e4/taumin)/log(1.4e4/2e-4)*background_boost)
    this%tauminn=0.95d0*taumin
    this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)

    do RW_i = 1, State%num_redshiftwindows
        !Make sure steps small enough for any features in source window functions
        associate (Win => State%Redshift_w(RW_i))
            if ((Win%kind /= window_21cm .or. .not. CP%transfer_21cm_cl) .and. &
                Win%sigma_tau/5/background_boost < Win%tau*(exp(this%dlntau)-1)) then
                this%dlntau = log(Win%sigma_tau/5/background_boost/Win%tau+1)
                nthermo = nint(log(State%tau0/this%tauminn)/this%dlntau) + 1
                this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)
            end if
        end associate
    end do
    this%nthermo = nthermo
    allocate(spline_data(nthermo), sdotmu(nthermo))

    if (allocated(this%tb) .and. this%nthermo/=size(this%tb)) then
        deallocate(this%scaleFactor, this%cs2, this%dcs2, this%ddotmu)
        deallocate(this%dscaleFactor, this%adot, this%dadot)
        deallocate(this%tb, this%xe, this%emmu, this%dotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    endif
    if (.not. allocated(this%tb)) then
        allocate(this%scaleFactor(nthermo), this%cs2(nthermo), this%dcs2(nthermo), this%ddotmu(nthermo))
        allocate(this%dscaleFactor(nthermo), this%adot(nthermo), this%dadot(nthermo))
        allocate(this%tb(nthermo), this%xe(nthermo), this%emmu(nthermo),this%dotmu(nthermo))
        allocate(this%demmu(nthermo), this%dddotmu(nthermo), this%ddddotmu(nthermo))
        if (dowinlens) allocate(this%winlens(nthermo), this%dwinlens(nthermo))
    end if

    if (State%num_redshiftwindows >0) then
        allocate(this%redshift_time(nthermo),this%dredshift_time(nthermo))
        allocate(this%arhos_fac(nthermo), this%darhos_fac(nthermo), this%ddarhos_fac(nthermo))
        allocate(RW(State%num_redshiftwindows))
    end if


    do RW_i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_w(RW_i))
            RedWin%tau_start = 0
            RedWin%tau_end = State%tau0
            if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                allocate(RW(RW_i)%awin_lens(nthermo))
                allocate(RW(RW_i)%dawin_lens(nthermo))
            end if
        end associate
    end do
    om = (State%grhob+State%grhoc)/&
        sqrt(3*(State%grhog+sum(State%grhormass(1:CP%Nu_mass_eigenstates))+State%grhornomass))
    a0=this%tauminn*State%adotrad*(1+om*this%tauminn/4)
    ninverse = nint(background_boost*log(1/a0)/log(1/2d-10)*4000)
    if (.not. CP%DarkEnergy%is_cosmological_constant) ninverse = ninverse*2

    nlin = ninverse/2
    allocate(scale_factors(ninverse+nlin))
    allocate(times(ninverse+nlin))
    allocate(dt(ninverse+nlin))
    allocate(taus(nthermo), xe_a(nthermo))

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    call CP%Recomb%Init(State,WantTSpin=CP%Do21cm)    !almost all the time spent here

    if (CP%Evolve_delta_xe) this%recombination_saha_tau  = State%TimeOfZ(CP%Recomb%get_saha_z(), tol=1e-4_dl)
    if (CP%Evolve_baryon_cs .or. CP%Evolve_delta_xe .or. CP%Evolve_delta_Ts .or. CP%Do21cm) &
        this%recombination_Tgas_tau = State%TimeOfz(1/CP%Recomb%min_a_evolve_Tm-1, tol=1e-4_dl)

    !$OMP SECTION
    !Do other stuff while recombination calculating
    awin_lens1=0
    awin_lens2=0
    transfer_ix =0

    call splini(spline_data,nthermo)

    this%tight_tau = 0
    this%actual_opt_depth = 0
    ncount=0
    this%z_drag=0.d0
    thomc0= Compton_CT * CP%tcmb**4
    this%r_drag0 = 3.d0/4.d0*State%grhob/State%grhog
    last_dotmu = 0

    this%matter_verydom_tau = 0
    a_verydom = CP%Accuracy%AccuracyBoost*5*(State%grhog+State%grhornomass)/(State%grhoc+State%grhob)
    if (CP%Reion%Reionization) then
        call CP%Reion%get_timesteps(State%reion_n_steps, reion_z_start, reion_z_complete)
        State%reion_tau_start = max(0.05_dl, State%TimeOfZ(reion_z_start, 1d-3))
        !Time when a very small reionization fraction (assuming tanh fitting)
        State%reion_tau_complete = min(State%tau0, &
            State%reion_tau_start+ State%DeltaTime(1/(1+reion_z_start),1/(1.d0+reion_z_complete),1d-3))
    else
        State%reion_tau_start = State%tau0
        State%reion_tau_complete = State%tau0
    end  if
    !  Initial conditions: assume radiation-dominated universe.
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.

    !Get scale factor as function of time by inverting tau(a)
    dlna = log(0.2_dl/a0)/(ninverse-1)
    do i=2, ninverse-1
        scale_factors(1+i) = a0*exp((i-1)*dlna)
    end do
    scale_factors(1) = a0
    scale_factors(2) = a0*exp(dlna/3)
    da = 0.8_dl/(nlin-2)
    do i=1, nlin-2
        scale_factors(ninverse+i) = 0.2_dl + (i-1)*da
    end do
    scale_factors(ninverse+nlin-1) = 0.9_dl + 0.1_dl*scale_factors(ninverse+nlin-2)
    scale_factors(ninverse+nlin) = 1
    do i=1, ninverse+nlin
        dt(i) = dtauda(State,scale_factors(i))
    end do
    call this%ScaleFactorAtTime%Init(scale_factors, dt)
    call this%ScaleFactorATTime%IntegralArray(times(2), first_index=2)
    times(1) = this%tauminn
    times(2:) = times(2:) + 2*(sqrt(1 + om*scale_factors(2)/ State%adotrad) -1)/om
    times(ninverse+nlin) = State%tau0
    call This%ScaleFactorAtTime%Init(times, scale_factors)
    taus(1) = this%tauminn
    do i=2,nthermo-1
        taus(i) = this%tauminn*exp((i-1)*this%dlntau)
    end do
    taus(nthermo) = State%tau0
    call this%ScaleFactorAtTime%Array(taus(2:), this%scaleFactor(2:))
    this%scaleFactor(1) = a0
    this%scaleFactor(nthermo) = 1
    this%adot(1) = 1/dtauda(State,a0)

    tau01=this%tauminn
    do i=2,nthermo
        !Get recombination-independent parts of background now as function of conformal time tau
        !Now at the log spaced time steps
        tau=taus(i)
        dtau = tau-tau01
        a = this%scaleFactor(i)
        adot = 1/dtauda(State,a)
        this%adot(i) = adot
        if (this%matter_verydom_tau ==0 .and. a > a_verydom) then
            this%matter_verydom_tau = tau
        end if
        z= 1._dl/a-1._dl
        if (State%num_redshiftwindows>0) then
            this%redshift_time(i) = z
            do RW_i = 1, State%num_redshiftwindows
                associate (Win => RW(RW_i), RedWin => State%Redshift_w(RW_i))
                    if (a > 1d-4) then
                        window = RedWin%Window%Window_f_a(a, winamp)

                        if  (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts  &
                            .and. CP%SourceTerms%counts_lensing) then
                            if (State%tau0 - tau > 2) then
                                dwing_lens =  adot * window *dtau
                                awin_lens1(RW_i) = awin_lens1(RW_i) + dwing_lens
                                awin_lens2(RW_i) = awin_lens2(RW_i) + dwing_lens/(State%tau0-tau)
                                Win%awin_lens(i) = awin_lens1(RW_i)/(State%tau0-tau) - awin_lens2(RW_i)
                            else
                                Win%awin_lens(i) = 0
                            end if
                        end if

                        if (RedWin%tau_start ==0 .and. winamp > 1e-8) then
                            RedWin%tau_start = tau01
                        else if (RedWin%tau_start /=0 .and. RedWin%tau_end==State%tau0 .and. winamp < 1e-8) then
                            RedWin%tau_end = min(State%tau0,tau + dtau)
                            if (DebugMsgs) call WriteFormat('Window %d: tau1 = %f, tau2 = %f',&
                                RW_i, RedWin%tau_start, RedWin%tau_end)
                        end if
                    else
                        if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts &
                            .and. CP%SourceTerms%counts_lensing) then
                            Win%awin_lens(i)=0
                        end if
                    end if
                end associate
            end do
        end if
        if (CP%WantTransfer .and.  CP%do21cm .and. CP%transfer_21cm_cl) then
            do RW_i = 1, CP%Transfer%PK_num_redshifts
                if (z< CP%Transfer%PK_redshifts(RW_i) .and. transfer_ix(RW_i)==0) then
                    transfer_ix(RW_i) = i
                end if
            end do
        end if
        tau01 =tau
    end do
    do RW_i = 1, State%num_redshiftwindows
        associate(Win => RW(RW_i))
            if (State%Redshift_w(RW_i)%kind == window_lensing .or. &
                State%Redshift_w(RW_i)%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                this%has_lensing_windows = .true.
                State%Redshift_w(RW_i)%has_lensing_window = .true.
                if (FeedbackLevel>0)  write(*,'(I1," Int W              = ",f9.6)') RW_i, awin_lens1(RW_i)
                Win%awin_lens=Win%awin_lens/awin_lens1(RW_i)
            else
                State%Redshift_w(RW_i)%has_lensing_window = .false.
            end if
        end associate
    end do
    !$OMP END PARALLEL SECTIONS

    if (global_error_flag/=0) return

    call CP%Recomb%xe_tm(a0,this%xe(1), this%tb(1))
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(1))
    this%cs2(1)=4._dl/3._dl*barssc*this%tb(1)
    this%dotmu(1)=this%xe(1)*State%akthom/a0**2


    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,16)
    do i=2,nthermo
        call CP%Recomb%xe_tm(this%scaleFactor(i), xe_a(i), this%tb(i))
    end do

    do i=2,nthermo
        tau =taus(i)
        a = this%scaleFactor(i)
        a2=a*a
        adot=this%adot(i)

        if (State%num_redshiftwindows>0) then
            if (a > 1d-4) then
                if (CP%Do21cm ) then
                    Tspin = CP%Recomb%T_s(a)
                    Trad = CP%TCMB/a
                    rho_fac = line21_const*State%NNow/a**3
                    tau_eps = a*line21_const*State%NNow/a**3/(adot/a)/Tspin/1000
                    this%arhos_fac(i) = (1-exp(-tau_eps))/tau_eps*a*rho_fac*(1 - Trad/Tspin)/(adot/a)
                    !         arhos_fac(i) = a*rho_fac*(1 - Trad/Tspin)/(adot/a)
                else
                    rho_fac = State%grhoc/(a2*a)
                    this%arhos_fac(i) = a*rho_fac/(adot/a)
                end if
            end if
        end if

        ! If there is re-ionization, smoothly increase xe to the
        ! requested value.
        if (CP%Reion%Reionization .and. tau > State%reion_tau_start) then
            if(ncount == 0) then
                ncount=i-1
            end if
            this%xe(i) = CP%Reion%x_e(1/a-1, tau, this%xe(ncount))
            if (CP%Accuracy%AccurateReionization .and. CP%WantDerivedParameters) then
                this%dotmu(i)=(xe_a(i) - this%xe(i))*State%akthom/a2

                if (last_dotmu /=0) then
                    this%actual_opt_depth = this%actual_opt_depth - 2._dl*(tau-taus(i-1))/(1._dl/this%dotmu(i)+1._dl/last_dotmu)
                end if
                last_dotmu = this%dotmu(i)
            end if
        else
            this%xe(i)=xe_a(i)
        end if

        !  approximate Baryon sound speed squared (over c**2).
        fe=(1._dl-CP%yhe)*this%xe(i)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        dtbdla=-2._dl*this%tb(i)
        if (a*this%tb(i)-CP%tcmb < -1e-8) then
            dtbdla= dtbdla -thomc0*fe/adot*(a*this%tb(i)-CP%tcmb)/a**3
        end if
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        this%cs2(i)=barssc*this%tb(i)*(1-dtbdla/this%tb(i)/3._dl)

        ! Calculation of the visibility function
        this%dotmu(i)=this%xe(i)*State%akthom/a2

        if (this%tight_tau==0 .and. 1/(tau*this%dotmu(i)) > 0.005) this%tight_tau = tau !0.005
        !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)
    end do

    if (CP%Reion%Reionization .and. (this%xe(nthermo) < 0.999d0)) then
        write(*,*)'Warning: xe at redshift zero is < 1'
        write(*,*) 'Check input parameters an Reionization_xe'
        write(*,*) 'function in the Reionization module'
    end if

    !Integrate for optical depth
    call dotmuSp%Init(taus(nthermo:1:-1), this%dotmu(nthermo:1:-1))
    allocate(opts(nthermo))
    call dotmuSp%IntegralArray(opts)
    sdotmu = opts(nthermo:1:-1)
    do j1=1,nthermo
        if (sdotmu(j1)< -69) then
            this%emmu(j1)=1.d-30
        else
            this%emmu(j1)=exp(sdotmu(j1))
            if (CP%Reion%Reionization .and. .not. CP%Accuracy%AccurateReionization .and. &
                this%actual_opt_depth==0 .and. this%xe(j1) < 1e-3) then
                this%actual_opt_depth = -sdotmu(j1)
            end if
        end if
    end do
    z_scale =  COBE_CMBTemp/CP%TCMB
    zstar_min = 700._dl * z_scale
    zstar_max = 2000._dl * z_scale
    if ((.not. CP%Reion%Reionization .or. CP%Accuracy%AccurateReionization) .and. CP%WantDerivedParameters) then
        do j1=nint(log(100/this%tauminn)/this%dlntau),nthermo
            if (-sdotmu(j1) - this%actual_opt_depth < 1) then
                !Bracket z_star
                zstar_min = 1/this%scaleFactor(j1+1)-1
                zstar_max = 1/this%scaleFactor(j1-2)-1
                exit
            end if
        end do
    end if

    if (CP%WantTransfer .and.  CP%do21cm .and. CP%transfer_21cm_cl) then
        if (allocated(State%optical_depths_for21cm)) deallocate(State%optical_depths_for21cm)
        allocate(State%optical_depths_for21cm(CP%Transfer%PK_num_redshifts))
        do RW_i = 1, CP%Transfer%PK_num_redshifts
            if (CP%Transfer%PK_Redshifts(RW_i) < 1e-3) then
                State%optical_depths_for21cm(RW_i) = 0 !zero may not be set correctly in transfer_ix
            else
                State%optical_depths_for21cm(RW_i) =  -sdotmu(transfer_ix(RW_i))
            end if
        end do
    end if

    if (CP%Reion%Reionization .and. CP%Accuracy%AccurateReionization &
        .and. FeedbackLevel > 0 .and. CP%WantDerivedParameters) then
        write(*,'("Reion opt depth      = ",f7.4)') this%actual_opt_depth
    end if

    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time of maximum visibility
    if (ncount == 0) then
        cf1=1._dl
        ns=nthermo
    else
        cf1=exp(-sdotmu(ncount))
        ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
        vis = this%emmu(j1)*this%dotmu(j1)
        tau = taus(j1)
        vfi=vfi+vis*cf1*this%dlntau*tau
        if ((iv == 0).and.(vfi > 1.0d-7/CP%Accuracy%AccuracyBoost)) then
            State%taurst=9._dl/10._dl*tau
            iv=1
        elseif (iv == 1) then
            if (vis > maxvis) then
                maxvis=vis
                State%tau_maxvis = tau
            end if
            if (vfi > 0.995) then
                State%taurend=tau
                iv=2
                exit
            end if
        end if
    end do

    if (iv /= 2) then
        call GlobalError('ThemoData Init: failed to find end of recombination',error_reionization)
        return
    end if

    if (dowinlens) then
        vfi=0
        awin_lens1p=0
        awin_lens2p=0
        this%winlens=0
        do j1=1,nthermo-1
            vis = this%emmu(j1)*this%dotmu(j1)
            tau = this%tauminn* taus(j1)
            vfi=vfi+vis*cf1*this%dlntau*tau
            if (vfi < 0.995) then
                dwing_lens =  vis*cf1*this%dlntau*tau / 0.995

                awin_lens1p = awin_lens1p + dwing_lens
                awin_lens2p = awin_lens2p + dwing_lens/(State%tau0-tau)
            end if
            this%winlens(j1)= awin_lens1p/(State%tau0-tau) - awin_lens2p
        end do
    end if

    ! Calculating the timesteps during recombination.

    if (CP%WantTensors) then
        State%dtaurec=min(State%dtaurec,State%taurst/160)/CP%Accuracy%AccuracyBoost
    else
        State%dtaurec=min(State%dtaurec,State%taurst/40)/CP%Accuracy%AccuracyBoost
        if (do_bispectrum .and. hard_bispectrum) State%dtaurec = State%dtaurec / 4
    end if

    if (CP%Reion%Reionization) State%taurend=min(State%taurend,State%reion_tau_start)

    if (DebugMsgs) then
        write (*,*) 'taurst, taurend = ', State%taurst, State%taurend
    end if

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    call splder(this%dotmu,this%ddotmu,nthermo,spline_data)
    call splder(this%ddotmu,this%dddotmu,nthermo,spline_data)
    call splder(this%dddotmu,this%ddddotmu,nthermo,spline_data)
    if (CP%want_zstar .or. CP%WantDerivedParameters) &
        this%z_star = State%binary_search(noreion_optdepth, 1.d0, zstar_min, zstar_max, &
        & 1d-3/background_boost, 100._dl*z_scale, 4000._dl*z_scale)
    !$OMP SECTION
    call splder(this%cs2,this%dcs2,nthermo,spline_data)
    call splder(this%emmu,this%demmu,nthermo,spline_data)
    call splder(this%adot,this%dadot,nthermo,spline_data)
    if (dowinlens) call splder(this%winlens,this%dwinlens,nthermo,spline_data)
    if (CP%want_zdrag .or. CP%WantDerivedParameters) &
        this%z_drag = State%binary_search(dragoptdepth, 1.d0, 800*z_scale, &
        & max(zstar_max*1.1_dl,1200._dl*z_scale), 2d-3/background_boost, 100.d0*z_scale, 4000._dl*z_scale)
    !$OMP SECTION
    this%ScaleFactor(:) = this%scaleFactor/taus !a/tau
    this%dScaleFactor(:) = (this%adot - this%ScaleFactor)*this%dlntau !derivative of a/tau
    if (State%num_redshiftwindows >0) then
        call splder(this%redshift_time,this%dredshift_time,nthermo,spline_data)
        call splder(this%arhos_fac,this%darhos_fac,nthermo,spline_data)
        call splder(this%darhos_fac,this%ddarhos_fac,nthermo,spline_data)
        do j2 = 1, State%num_redshiftwindows
            if (State%Redshift_w(j2)%has_lensing_window) then
                call splder(RW(j2)%awin_lens,RW(j2)%dawin_lens,nthermo,spline_data)
            end if
        end do
    end if
    call this%SetTimeSteps(State,State%TimeSteps)
    !$OMP END PARALLEL SECTIONS

    if (State%num_redshiftwindows>0) then
        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
        do j2=1,State%TimeSteps%npoints
            call this%DoWindowSpline(State,j2,State%TimeSteps%points(j2), RW)
        end do
        !$OMP END PARALLEL DO
        call this%SetTimeStepWindows(State,State%TimeSteps)
    end if

    if (CP%WantDerivedParameters) then
        associate(ThermoDerivedParams => State%ThermoDerivedParams)
            !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
            !$OMP SECTION
            ThermoDerivedParams( derived_Age ) = State%DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
            rstar =State%sound_horizon(this%z_star)
            ThermoDerivedParams( derived_rstar ) = rstar
            DA = State%AngularDiameterDistance(this%z_star)/(1/(this%z_star+1))
            ThermoDerivedParams( derived_zdrag ) = this%z_drag
            !$OMP SECTION
            rs =State%sound_horizon(this%z_drag)
            ThermoDerivedParams( derived_rdrag ) = rs
            ThermoDerivedParams( derived_kD ) =  &
                sqrt(1.d0/(Integrate_Romberg(State,ddamping_da, 1d-8, 1/(this%z_star+1), 1d-6)/6))
            !$OMP SECTION
            ThermoDerivedParams( derived_zEQ ) = State%z_eq
            a_eq = 1/(1+State%z_eq)
            ThermoDerivedParams( derived_kEQ ) = 1/(a_eq*dtauda(State,a_eq))
            rs_eq = State%sound_horizon(State%z_eq)
            tau_eq = State%timeOfz(State%z_eq)
            !$OMP SECTION
            !$OMP END PARALLEL SECTIONS

            ThermoDerivedParams( derived_zstar ) = this%z_star
            ThermoDerivedParams( derived_thetastar ) = 100*rstar/DA
            ThermoDerivedParams( derived_DAstar ) = DA/1000
            ThermoDerivedParams( derived_thetaEQ ) = 100*tau_eq/DA
            ThermoDerivedParams( derived_theta_rs_EQ ) = 100*rs_EQ/DA
            ThermoDerivedParams( derived_thetaD ) =  100*const_pi/ThermoDerivedParams( derived_kD )/DA

            if (allocated(CP%z_outputs)) then
                if (allocated(State%BackgroundOutputs%H)) &
                    deallocate(State%BackgroundOutputs%H, State%BackgroundOutputs%DA, State%BackgroundOutputs%rs_by_D_v)
                noutput = size(CP%z_outputs)
                allocate(State%BackgroundOutputs%H(noutput), State%BackgroundOutputs%DA(noutput), &
                    State%BackgroundOutputs%rs_by_D_v(noutput))
                !$OMP PARALLEL DO DEFAULT(shared)
                do i=1,noutput
                    State%BackgroundOutputs%H(i) = State%HofZ(CP%z_outputs(i))
                    State%BackgroundOutputs%DA(i) = State%AngularDiameterDistance(CP%z_outputs(i))
                    State%BackgroundOutputs%rs_by_D_v(i) = rs/BAO_D_v_from_DA_H(CP%z_outputs(i), &
                        State%BackgroundOutputs%DA(i),State%BackgroundOutputs%H(i))
                end do
            end if

            if (FeedbackLevel > 0) then
                write(*,'("Age of universe/GYr  = ",f7.3)') ThermoDerivedParams( derived_Age )
                write(*,'("zstar                = ",f8.2)') ThermoDerivedParams( derived_zstar )
                write(*,'("r_s(zstar)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rstar )
                write(*,'("100*theta            = ",f9.6)') ThermoDerivedParams( derived_thetastar )
                write(*,'("DA(zstar)/Gpc        = ",f9.5)') ThermoDerivedParams( derived_DAstar )

                write(*,'("zdrag                = ",f8.2)') ThermoDerivedParams( derived_zdrag )
                write(*,'("r_s(zdrag)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rdrag )

                write(*,'("k_D(zstar) Mpc       = ",f7.4)') ThermoDerivedParams( derived_kD )
                write(*,'("100*theta_D          = ",f9.6)') ThermoDerivedParams( derived_thetaD )

                write(*,'("z_EQ (if v_nu=1)     = ",f8.2)') ThermoDerivedParams( derived_zEQ )
                write(*,'("k_EQ Mpc (if v_nu=1) = ",f9.6)') ThermoDerivedParams( derived_kEQ )
                write(*,'("100*theta_EQ         = ",f9.6)') ThermoDerivedParams( derived_thetaEQ )
                write(*,'("100*theta_rs_EQ      = ",f9.6)') ThermoDerivedParams( derived_theta_rs_EQ )
            end if
        end associate
    end if

    !Sources
    if(State%num_redshiftwindows>0) then
        deallocate(RW,this%arhos_fac, this%darhos_fac, this%ddarhos_fac)
        deallocate(this%redshift_time, this%dredshift_time)
        do RW_i = 1, State%num_redshiftwindows
            associate (RedWin => State%Redshift_W(RW_i))
                if (RedWin%kind == window_21cm) then
                    outstr = 'z= '//trim(RealToStr(real(RedWin%Redshift),4))//': T_b = '//trim(RealToStr(real(RedWin%Fq),6))// &
                        'mK; tau21 = '//trim(RealToStr(real(RedWin%optical_depth_21),5))
                    write (*,*) RW_i,trim(outstr)
                end if
            end associate
        end do
    end if

    if (CP%Do21cm .and. CP%transfer_21cm_cl) then
        do RW_i=1,CP%Transfer%PK_num_redshifts
            a=1._dl/(1+CP%Transfer%PK_redshifts(RW_i))
            Tspin = CP%Recomb%T_s(a)
            Trad = CP%TCMB/a
            adot = 1/dtauda(State,a)
            tau_eps = a**2*line21_const*State%NNow/a**3/adot/Tspin/1000
            Tb21cm = 1000*(1-exp(-tau_eps))*a*(Tspin-Trad)
            if (FeedbackLevel>0) then
                outstr = 'z= '//trim(RealToStr(real(CP%Transfer%PK_redshifts(RW_i))))// &
                    ': tau_21cm = '//trim(RealToStr(real(tau_eps),5))//'; T_b = '//trim(RealToStr(real(Tb21cm),6))//'mK'
                write (*,*) trim(outstr)
            end if
        end do
    end if

    this%HasThermoData = .true.
    end subroutine Thermo_Init


    subroutine SetTimeSteps(this,State,TimeSteps)
    !Set time steps to use for sampling the source functions for the CMB power spectra
    class(TThermoData) :: this
    Type(TRanges) :: TimeSteps
    class(CAMBdata) State
    real(dl) dtau0
    integer nri0, nstep
    !Sources
    integer ix,i,nwindow, L_limb
    real(dl) keff, win_end, TimeSampleBoost, delta

    TimeSampleBoost = State%CP%Accuracy%AccuracyBoost*State%CP%Accuracy%TimeStepBoost
    call TimeSteps%Init()

    call TimeSteps%Add_delta(State%taurst, State%taurend, State%dtaurec)

    ! Calculating the timesteps after recombination
    if (State%CP%WantTensors) then
        dtau0=max(State%taurst/40,State%tau0/2000._dl/TimeSampleBoost)
    else
        dtau0=State%tau0/500._dl/TimeSampleBoost
        if (do_bispectrum) dtau0 = dtau0/3
        !Don't need this since adding in Limber on small scales
        !  if (CP%DoLensing) dtau0=dtau0/2
        !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
    end if

    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)

    !Sources

    this%tau_start_redshiftwindows = State%tau0
    this%tau_end_redshiftwindows = 0
    if (State%CP%WantScalars .or. State%CP%SourceTerms%line_reionization) then
        do ix=1, State%num_redshiftwindows
            associate (Win => State%Redshift_W(ix))

                Win%tau_start = max(Win%tau_start,state%taurst)

                if (Win%kind /= window_lensing) then
                    !Have to be careful to integrate dwinV as the window tails off
                    this%tau_end_redshiftwindows = max(Win%tau_end,this%tau_end_redshiftwindows)
                    nwindow = nint(150*TimeSampleBoost)
                    win_end = Win%tau_end
                else !lensing
                    nwindow = nint(TimeSampleBoost*Win%chi0/100)
                    win_end = State%tau0
                end if

                if (Win%kind == window_21cm .and. (State%CP%SourceTerms%line_phot_dipole .or. &
                    State%CP%SourceTerms%line_phot_quadrupole)) nwindow = nwindow *3

                L_limb = Win_limber_ell(Win,State%CP,State%CP%max_l)
                keff = WindowKmaxForL(Win,State%CP,L_limb)

                !Keep sampling in x better than Nyquist
                nwindow = max(nwindow, nint(TimeSampleBoost *(win_end- Win%tau_start)* keff/3))
                if (Win%kind /= window_lensing .and. Win%sigma_tau < (win_end - Win%tau_start)/nwindow*8) then
                    nwindow = nint(TimeSampleBoost*(win_end - Win%tau_start)/Win%sigma_tau*8)
                end if
                if (Feedbacklevel > 1 .or. DebugMsgs) call WriteFormat('nwindow %d: %d', ix, nwindow)
                delta = (win_end-Win%tau_start)/nwindow
                Win%tau_start = Win%tau_start - delta*3
                win_end = min(State%tau0, win_end+delta*2)

                this%tau_start_redshiftwindows = min(Win%tau_start,this%tau_start_redshiftwindows)

                call TimeSteps%Add(Win%tau_start, win_end, nwindow)
                !This should cover whole range where not tiny

                if (Win%kind /= window_lensing  .and. &
                    Win%tau_peakend-Win%tau_peakstart < nint(60*TimeSampleBoost) * delta) then
                    call TimeSteps%Add(Win%tau_peakstart,Win%tau_peakend, nint(60*TimeSampleBoost))
                    !This should be over peak
                end if
            end associate
        end do
    end if


    if (State%CP%Reion%Reionization) then
        nri0=int(State%reion_n_steps*State%CP%Accuracy%AccuracyBoost)
        !Steps while reionization going from zero to maximum
        call TimeSteps%Add(State%reion_tau_start,State%reion_tau_complete,nri0)
    end if

    !Sources
    if (.not. State%CP%Want_CMB .and. State%CP%WantCls) then
        if (State%num_redshiftwindows==0) then
            call GlobalError('Want_CMB=false and WantCls=true, but no redshift windows either', &
                error_unsupported_params)
        else
            call TimeSteps%Add_delta(this%tau_start_redshiftwindows, State%tau0, dtau0)
        endif
    end if

    if (global_error_flag/=0) then
        return
    end if



    !Create arrays out of the region information.
    call TimeSteps%GetArray()
    nstep = TimeSteps%npoints

    if (State%num_redshiftwindows > 0) then
        if (allocated(this%step_redshift)) deallocate(this%step_redshift, this%rhos_fac, this%drhos_fac)
        allocate(this%step_redshift(nstep), this%rhos_fac(nstep), this%drhos_fac(nstep))
        do i=1,State%num_redshiftwindows
            associate (Win => State%Redshift_W(i))
                allocate(Win%winF(nstep),Win%wing(nstep),Win%dwing(nstep),Win%ddwing(nstep), &
                    Win%winV(nstep),Win%dwinV(nstep),Win%ddwinV(nstep))
                allocate(Win%win_lens(nstep),Win%wing2(nstep),Win%dwing2(nstep),Win%ddwing2(nstep))
                allocate(Win%wingtau(nstep),Win%dwingtau(nstep),Win%ddwingtau(nstep))
                if (Win%kind == window_counts) then
                    allocate(Win%comoving_density_ev(nstep))
                end if
            end associate
        end do
    end if

    if (DebugMsgs .and. FeedbackLevel > 0) call WriteFormat('Set %d time steps', nstep)

    end subroutine SetTimeSteps

    subroutine SetTimeStepWindows(this,State,TimeSteps)
    use constants
    class(TThermoData) :: this
    class(CAMBdata) :: State
    Type(TRanges) :: TimeSteps
    integer i, j, jstart, ix
    real(dl) tau, a, a2
    real(dl) Tspin, Trad, rho_fac, tau_eps
    real(dl) window, winamp
    real(dl) z,rhos, adot, exp_fac
    real(dl) tmp(TimeSteps%npoints), tmp2(TimeSteps%npoints), hubble_tmp(TimeSteps%npoints)
    real(dl), allocatable , dimension(:,:) :: int_tmp, back_count_tmp
    integer ninterp

    ! Prevent false positive warnings for uninitialized
    Tspin = 0._dl
    Trad = 0._dl
    tau_eps = 0._dl
    exp_fac = 0._dl

    jstart = TimeSteps%IndexOf(this%tau_start_redshiftwindows)
    ninterp = TimeSteps%npoints - jstart + 1

    do i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_W(i))
            RedWin%wing=0
            RedWin%winV=0
            RedWin%winF=0
            RedWin%wing2=0
            RedWin%dwing=0
            RedWin%dwinV=0
            RedWin%dwing2=0
            RedWin%ddwing=0
            RedWin%ddwinV=0
            RedWin%ddwing2=0
            RedWin%wingtau=0
            RedWin%dwingtau=0
            RedWin%ddwingtau=0
            RedWin%Fq = 0
            if (RedWin%kind == window_counts) then
                RedWin%comoving_density_ev  = 0
            end if
        end associate
    end do

    allocate(int_tmp(jstart:TimeSteps%npoints,State%num_redshiftwindows))
    int_tmp = 0
    allocate(back_count_tmp(jstart:TimeSteps%npoints,State%num_redshiftwindows))
    back_count_tmp = 0

    do j=jstart, TimeSteps%npoints
        tau = TimeSteps%points(j)
        z = this%step_redshift(j)
        a = 1._dl/(1._dl+z)
        a2=a**2
        adot=1._dl/dtauda(State,a)


        if (State%CP%Do21cm) then
            Tspin = State%CP%Recomb%T_s(a)
            Trad = State%CP%TCMB/a
            rho_fac = line21_const*State%NNow/a**3 !neglect ionization fraction
            tau_eps = a*rho_fac/(adot/a)/Tspin/1000
            exp_fac =  (1-exp(-tau_eps))/tau_eps
        else
            rho_fac = State%grhoc/a**3
        end if

        hubble_tmp(j) = adot/a

        do i = 1, State%num_redshiftwindows
            associate (RedWin => State%Redshift_W(i))
                if (tau < RedWin%tau_start) cycle

                window = RedWin%Window%Window_f_a(a, winamp)

                if (RedWin%kind == window_21cm) then
                    rhos = rho_fac*(1 - Trad/Tspin)

                    !Want to integrate this...
                    int_tmp(j,i) = this%drhos_fac(j)*a*window

                    RedWin%WinV(j) = -exp(-tau_eps)*a2*rhos*window/(adot/a)

                    RedWin%wing(j) = exp_fac*a2*rhos*window

                    !The window that doesn't go to zero at T_s = T_gamma
                    RedWin%wing2(j) = exp_fac*a2*rho_fac*Trad/Tspin*window

                    !Window for tau_s for the self-absoption term
                    RedWin%wingtau(j) =  RedWin%wing(j)*(1 - exp(-tau_eps)/exp_fac)
                elseif (RedWin%kind == window_counts) then

                    !window is n(a) where n is TOTAL not fractional number
                    !delta = int wing(eta) delta(eta) deta
                    RedWin%wing(j) = adot *window

                    !Window with 1/H in
                    RedWin%wing2(j) = RedWin%wing(j)/(adot/a)

                    !winv is g/chi for the ISW and time delay terms
                    RedWin%WinV(j) = 0
                    if (tau < State%tau0 -0.1) then
                        int_tmp(j,i) = RedWin%wing(j)/(State%tau0 - tau)
                    else
                        int_tmp(j,i)=0
                    end if

                    if (State%CP%SourceTerms%counts_evolve) then
                        back_count_tmp(j,i) =  RedWin%Window%counts_background_z(1/a-1)/a
                        if (tau < State%tau0 -0.1) then
                            RedWin%comoving_density_ev(j) = back_count_tmp(j,i)*(adot/a)/(State%tau0 - tau)**2
                        else
                            RedWin%comoving_density_ev(j) = 0
                        end if
                    end if
                end if
            end associate
        end do
    end do

    do i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_W(i))

            ! int (a*rho_s/H)' a W_f(a) d\eta, or for counts int g/chi deta
            call spline_def(TimeSteps%points(jstart:),int_tmp(jstart:,i),ninterp,tmp)
            call spline_integrate(TimeSteps%points(jstart:),int_tmp(jstart:,i),tmp, tmp2(jstart:),ninterp)
            RedWin%WinV(jstart:TimeSteps%npoints) =  &
                RedWin%WinV(jstart:TimeSteps%npoints) + tmp2(jstart:TimeSteps%npoints)

            call spline_def(TimeSteps%points(jstart:),RedWin%WinV(jstart:),ninterp,RedWin%ddWinV(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%WinV(jstart:),RedWin%ddWinV(jstart:), RedWin%dWinV(jstart:), ninterp)

            call spline_def(TimeSteps%points(jstart:),RedWin%Wing(jstart:),ninterp,RedWin%ddWing(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%Wing(jstart:),RedWin%ddWing(jstart:), RedWin%dWing(jstart:), ninterp)

            call spline_def(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),ninterp,RedWin%ddWing2(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),RedWin%ddWing2(jstart:), &
                RedWin%dWing2(jstart:), ninterp)

            call spline_integrate(TimeSteps%points(jstart:),RedWin%Wing(jstart:), &
                RedWin%ddWing(jstart:), RedWin%WinF(jstart:),ninterp)
            RedWin%Fq = RedWin%WinF(TimeSteps%npoints)

            if (RedWin%kind == window_21cm) then
                call spline_integrate(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),&
                    RedWin%ddWing2(jstart:), tmp(jstart:),ninterp)
                RedWin%optical_depth_21 = tmp(TimeSteps%npoints) / (State%CP%TCMB*1000)
                !WinF not used.. replaced below

                call spline_def(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),ninterp,RedWin%ddWingtau(jstart:))
                call spline_deriv(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),RedWin%ddWingtau(jstart:), &
                    RedWin%dWingtau(jstart:), ninterp)
            elseif (RedWin%kind == window_counts) then

                if (State%CP%SourceTerms%counts_evolve) then
                    call spline_def(TimeSteps%points(jstart:),back_count_tmp(jstart:,i),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),back_count_tmp(jstart:,i),tmp,tmp2(jstart:),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%Wingtau(ix) = 0
                        else
                            !evo bias is computed with total derivative
                            RedWin%Wingtau(ix) =  -tmp2(ix) * RedWin%Wing(ix) / (back_count_tmp(ix,i)*hubble_tmp(ix)) &
                                !+ 5*RedWin%dlog10Ndm * ( RedWin%Wing(ix)- int_tmp(ix,i)/hubble_tmp(ix))
                                !The correction from total to partial derivative takes 1/adot(tau0-tau) cancels
                                + 10*RedWin%Window%dlog10Ndm * RedWin%Wing(ix)
                        end if
                    end do

                    !comoving_density_ev is d log(a^3 n_s)/d eta * window
                    call spline_def(TimeSteps%points(jstart:),RedWin%comoving_density_ev(jstart:),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),RedWin%comoving_density_ev(jstart:),tmp,tmp2(jstart:),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%comoving_density_ev(ix) = 0
                        elseif (RedWin%comoving_density_ev(ix)/=0._dl) then
                            !correction needs to be introduced from total derivative to partial derivative
                            RedWin%comoving_density_ev(ix) =   tmp2(ix) / RedWin%comoving_density_ev(ix) &
                                -5*RedWin%Window%dlog10Ndm * ( hubble_tmp(ix) + int_tmp(ix,i)/RedWin%Wing(ix))
                        end if
                    end do
                else
                    RedWin%comoving_density_ev=0
                    call spline_def(TimeSteps%points(jstart:),hubble_tmp(jstart:),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),hubble_tmp(jstart:),tmp, tmp2(jstart:), ninterp)

                    !assume d( a^3 n_s) of background population is zero, so remaining terms are
                    !wingtau =  g*(2/H\chi + Hdot/H^2)  when s=0; int_tmp = window/chi
                    RedWin%Wingtau(jstart:TimeSteps%npoints) = &
                        2*(1-2.5*RedWin%Window%dlog10Ndm)*int_tmp(jstart:TimeSteps%npoints,i)/&
                        hubble_tmp(jstart:TimeSteps%npoints)&
                        + 5*RedWin%Window%dlog10Ndm*RedWin%Wing(jstart:TimeSteps%npoints) &
                        + tmp2(jstart:TimeSteps%npoints)/hubble_tmp(jstart:TimeSteps%npoints)**2 &
                        *RedWin%Wing(jstart:TimeSteps%npoints)
                endif

                call spline_def(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),ninterp, &
                    RedWin%ddWingtau(jstart:))
                call spline_deriv(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),RedWin%ddWingtau(jstart:), &
                    RedWin%dWingtau(jstart:), ninterp)

                !WinF is int[ g*(...)]
                call spline_integrate(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),&
                    RedWin%ddWingtau(jstart:), RedWin%WinF(jstart:),ninterp)
            end if
        end associate
    end do

    end subroutine SetTimeStepWindows


    subroutine interp_window(TimeSteps,RedWin,tau,wing_t, wing2_t, winv_t)
    !for evolving sources for reionization we neglect wingtau self-absorption
    Type(TRanges) :: TimeSteps
    Type(TRedWin)  :: RedWin
    integer i
    real(dl) :: tau, wing_t, wing2_t,winv_t
    real(dl) a0,b0,ho

    i = TimeSteps%IndexOf(tau)
    if (i< TimeSteps%npoints) then
        ho=TimeSteps%points(i+1)-TimeSteps%points(i)
        a0=(TimeSteps%points(i+1)-tau)/ho
        b0=1-a0
        wing_t = a0*RedWin%wing(i)+ b0*RedWin%wing(i+1)+((a0**3-a0)* RedWin%ddwing(i) &
            +(b0**3-b0)*RedWin%ddwing(i+1))*ho**2/6
        wing2_t = a0*RedWin%wing2(i)+ b0*RedWin%wing2(i+1)+((a0**3-a0)* RedWin%ddwing2(i) &
            +(b0**3-b0)*RedWin%ddwing2(i+1))*ho**2/6
        winv_t = a0*RedWin%winv(i)+ b0*RedWin%winv(i+1)+((a0**3-a0)* RedWin%ddwinv(i) &
            +(b0**3-b0)*RedWin%ddwinv(i+1))*ho**2/6
    else
        wing_t = 0
        wing2_t = 0
        winv_t = 0
    end if
    end subroutine interp_window

    subroutine DoWindowSpline(this,State,j2,tau, RW)
    class(TThermoData) :: this
    class(CAMBdata) :: State
    integer j2, i, RW_i
    real(dl) d, tau
    Type(CalWins) :: RW(:)

    !     Cubic-spline interpolation.
    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)

    d=d-i
    if (i < this%nthermo) then

        this%step_redshift(j2) = this%redshift_time(i)+d*(this%dredshift_time(i)+ &
            d*(3._dl*(this%redshift_time(i+1)-this%redshift_time(i)) &
            -2._dl*this%dredshift_time(i)-this%dredshift_time(i+1)+d*(this%dredshift_time(i)+this%dredshift_time(i+1) &
            +2._dl*(this%redshift_time(i)-this%redshift_time(i+1)))))

        this%rhos_fac(j2) = this%arhos_fac(i)+d*(this%darhos_fac(i)+d*(3._dl*(this%arhos_fac(i+1)-this%arhos_fac(i)) &
            -2._dl*this%darhos_fac(i)-this%darhos_fac(i+1)+d*(this%darhos_fac(i)+this%darhos_fac(i+1) &
            +2._dl*(this%arhos_fac(i)-this%arhos_fac(i+1)))))
        this%drhos_fac(j2) = (this%darhos_fac(i)+d*(this%ddarhos_fac(i)+d*(3._dl*(this%darhos_fac(i+1)  &
            -this%darhos_fac(i))-2._dl*this%ddarhos_fac(i)-this%ddarhos_fac(i+1)+d*(this%ddarhos_fac(i) &
            +this%ddarhos_fac(i+1)+2._dl*(this%darhos_fac(i)-this%darhos_fac(i+1))))))/(tau &
            *this%dlntau)

        do RW_i=1, State%num_redshiftwindows
            if (State%Redshift_w(RW_i)%has_lensing_window) then
                associate(W => State%Redshift_W(RW_i), C=> RW(RW_i))

                    W%win_lens(j2) = C%awin_lens(i)+d*(C%dawin_lens(i)+d*(3._dl*(C%awin_lens(i+1)-C%awin_lens(i)) &
                        -2._dl*C%dawin_lens(i)-C%dawin_lens(i+1)+d*(C%dawin_lens(i)+C%dawin_lens(i+1) &
                        +2._dl*(C%awin_lens(i)-C%awin_lens(i+1)))))
                end associate
            end if
        end do

    else
        this%step_redshift(j2) = 0
        this%rhos_fac(j2)=0
        this%drhos_fac(j2)=0
        do RW_i=1, State%num_redshiftwindows
            associate (W => State%Redshift_W(RW_i))
                W%win_lens(j2)=0
            end associate
        end do
    end if

    end subroutine DoWindowSpline

    subroutine IonizationFunctionsAtTime(this,tau, a, opac, dopac, ddopac, &
        vis, dvis, ddvis, expmmu, lenswin)
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out):: a, opac, dopac, ddopac, vis, dvis, ddvis, expmmu, lenswin
    real(dl) d, cs2
    integer i

    call this%Values(tau,a,cs2,opac,dopac)

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i

    if (i < this%nthermo) then
        ddopac=(this%dddotmu(i)+d*(this%ddddotmu(i)+d*(3._dl*(this%dddotmu(i+1) &
            -this%dddotmu(i))-2._dl*this%ddddotmu(i)-this%ddddotmu(i+1)  &
            +d*(this%ddddotmu(i)+this%ddddotmu(i+1)+2._dl*(this%dddotmu(i) &
            -this%dddotmu(i+1)))))-(this%dlntau**2)*tau*dopac) &
            /(tau*this%dlntau)**2
        expmmu=this%emmu(i)+d*(this%demmu(i)+d*(3._dl*(this%emmu(i+1)-this%emmu(i)) &
            -2._dl*this%demmu(i)-this%demmu(i+1)+d*(this%demmu(i)+this%demmu(i+1) &
            +2._dl*(this%emmu(i)-this%emmu(i+1)))))

        if (dowinlens) then
            lenswin=this%winlens(i)+d*(this%dwinlens(i)+d*(3._dl*(this%winlens(i+1)-this%winlens(i)) &
                -2._dl*this%dwinlens(i)-this%dwinlens(i+1)+d*(this%dwinlens(i)+this%dwinlens(i+1) &
                +2._dl*(this%winlens(i)-this%winlens(i+1)))))
        end if
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3*opac*dopac+ddopac)
    else
        ddopac=this%dddotmu(this%nthermo)
        expmmu=this%emmu(this%nthermo)
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3._dl*opac*dopac+ddopac)
    end if

    end subroutine IonizationFunctionsAtTime

    subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call CTrans%q%getArray(.true.)

    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,&
        min(CTrans%max_index_nonlimber,CTrans%ls%nl), CTrans%q%npoints), STAT = st)
    if (st /= 0) call MpiStop('Init_ClTransfer: Error allocating memory for transfer functions')
    CTrans%Delta_p_l_k = 0

    end subroutine Init_ClTransfer

    subroutine Init_Limber(CTrans,State)
    Type(ClTransferData) :: CTrans
    class(CAMBdata) :: State

    call Free_Limber(Ctrans)
    allocate(CTrans%Limber_l_min(CTrans%NumSources))
    CTrans%Limber_l_min = 0
    if (State%num_redshiftwindows>0 .or. State%CP%SourceTerms%limber_phi_lmin>0) then
        allocate(CTrans%Limber_windows(CTrans%NumSources,CTrans%ls%nl))
    end if

    end subroutine Init_Limber

    subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Delta_p_l_k)) deallocate(CTrans%Delta_p_l_k)
    call CTrans%q%Free()
    call Free_Limber(CTrans)

    end subroutine Free_ClTransfer

    subroutine Free_Limber(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Limber_l_min)) deallocate(CTrans%Limber_l_min)
    if (allocated(CTrans%Limber_windows)) deallocate(CTrans%Limber_windows)

    end subroutine Free_Limber


    function Win_Limber_ell(W,CP,lmax) result(ell_limb)
    Type(TRedWin) :: W
    Type(CAMBParams) :: CP
    integer, intent(in) :: lmax
    integer ell_limb
    real(dl) LimBoost

    if (CP%SourceTerms%limber_windows) then
        LimBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%SourceLimberBoost
        !Turn on limber when k is a scale smaller than window width
        if (W%kind==window_lensing) then
            ell_limb = max(CP%SourceTerms%limber_phi_lmin,nint(50*LimBoost))
        else
            ell_limb = max(CP%SourceTerms%limber_phi_lmin, nint(LimBoost*6*W%chi0/W%sigma_tau))
        end if
    else
        ell_limb = lmax
    end if
    end function Win_Limber_ell


    subroutine TCLdata_InitCls(this, State)
    class(TCLData) :: this
    class(CAMBdata) :: State

    associate(CP=>State%CP)
        call CheckLoadedHighLTemplate
        if (CP%WantScalars) then
            if (allocated(this%Cl_scalar)) deallocate(this%Cl_scalar)
            allocate(this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:State%Scalar_C_last), source=0._dl)
            if (CP%want_cl_2D_array) then
                if (allocated(this%Cl_scalar_array)) deallocate(this%Cl_scalar_array)
                allocate(this%Cl_scalar_Array(CP%Min_l:CP%Max_l, &
                    3+State%num_redshiftwindows+CP%CustomSources%num_custom_sources, &
                    3+State%num_redshiftwindows+CP%CustomSources%num_custom_sources))
                this%Cl_scalar_array = 0
            end if
        end if

        if (CP%WantVectors) then
            if (allocated(this%Cl_vector)) deallocate(this%Cl_vector)
            allocate(this%Cl_vector(CP%Min_l:CP%Max_l, CT_Temp:CT_Cross), source=0._dl)
        end if

        if (CP%WantTensors) then
            if (allocated(this%Cl_tensor)) deallocate(this%Cl_tensor)
            allocate(this%Cl_tensor(CP%Max_l_tensor, CT_Temp:CT_Cross), source=0._dl)
        end if
    end associate

    end subroutine TCLdata_InitCls

    function open_file_header(filename, Col1, Columns, n) result(unit)
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: col1
    character(LEN=name_tag_len), intent(in) :: Columns(:)
    integer, intent(in), optional :: n
    integer :: unit, nn

    nn = PresentDefault(6, n)

    open(newunit=unit,file=filename,form='formatted',status='replace')
    if (output_file_headers) then
        write(unit,'("#",1A'//Trim(IntToStr(nn-1))//'," ",*(A15))') Col1, Columns
    end if

    end function open_file_header


    function scalar_fieldname(i)
    integer, intent(in) :: i
    character(LEN=5) :: scalar_fieldname
    character(LEN=3), parameter :: scalar_fieldnames = 'TEP'

    if (i<=3) then
        scalar_fieldname = scalar_fieldnames(i:i)
    else
        scalar_fieldname = 'W'//trim(IntToStr(i-3))
    end if

    end function scalar_fieldname

    subroutine TCLdata_output_cl_files(this, State, ScalFile,ScalCovFile,TensFile, TotFile, LensFile, LensTotFile, factor)
    class(TCLData) :: this
    class(CAMBdata), target :: State
    character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile,ScalCovfile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact
    integer :: last_C, il, i, j, unit
    real(dl), allocatable :: outarr(:,:)
    character(LEN=name_tag_len) :: cov_names((3+State%num_redshiftwindows)**2)
    Type(CAMBParams), pointer :: CP
    integer lmin

    CP=> State%CP
    lmin= CP%Min_l

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. ScalFile /= '') then
        last_C=min(C_PhiTemp,State%Scalar_C_last)
        unit = open_file_header(ScalFile, 'L', C_name_tags(:last_C))
        do il=lmin,min(10000,CP%Max_l)
            write(unit,trim(numcat('(1I6,',last_C))//'E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        do il=10100,CP%Max_l, 100
            write(unit,trim(numcat('(1E15.6,',last_C))//'E15.6)') real(il),&
                fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%want_cl_2D_array .and. ScalCovFile /= '' .and. this%CTransScal%NumSources>2) then
        allocate(outarr(1:3+State%num_redshiftwindows,1:3+State%num_redshiftwindows))

        do i=1, 3+State%num_redshiftwindows
            do j=1, 3+State%num_redshiftwindows
                cov_names(j + (i-1)*(3+State%num_redshiftwindows)) = trim(scalar_fieldname(i))//'x'//trim(scalar_fieldname(j))
            end do
        end do
        unit = open_file_header(ScalCovFile, 'L', cov_names)

        do il=lmin,min(10000,CP%Max_l)
            outarr=this%Cl_scalar_array(il,1:3+State%num_redshiftwindows,1:3+State%num_redshiftwindows)
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            write(unit,trim(numcat('(1I6,',(3+State%num_redshiftwindows)**2))//'E15.6)') il, real(outarr)
        end do
        do il=10100,CP%Max_l, 100
            outarr=this%Cl_scalar_array(il,1:3+State%num_redshiftwindows,1:3+State%num_redshiftwindows)
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            write(unit,trim(numcat('(1E15.5,',(3+State%num_redshiftwindows)**2))//'E15.6)') real(il), real(outarr)
        end do
        close(unit)
        deallocate(outarr)
    end if

    if (CP%WantTensors .and. TensFile /= '') then
        unit = open_file_header(TensFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_tensor(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if

    if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
        unit = open_file_header(TotFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_scalar(il, C_Temp:C_E)+ this%Cl_tensor(il,C_Temp:C_E)), &
                fact*this%Cl_tensor(il, CT_B), fact*(this%Cl_scalar(il, C_Cross) + this%Cl_tensor(il, CT_Cross))
        end do
        do il=CP%Max_l_tensor+1,CP%Max_l
            write(unit,'(1I6,4E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:C_E), 0._dl, fact*this%Cl_scalar(il,C_Cross)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
        unit = open_file_header(LensFile, 'L', CT_name_tags)
        do il=lmin, this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if


    if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
        unit = open_file_header(LensTotFile, 'L', CT_name_tags)
        do il=lmin,min(CP%Max_l_tensor,this%lmax_lensed)
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_lensed(il,CT_Temp:CT_Cross)+ &
                this%Cl_tensor(il, CT_Temp:CT_Cross))
        end do
        do il=min(CP%Max_l_tensor,this%lmax_lensed)+1,this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_cl_files

    subroutine TCLdata_output_lens_pot_files(this,CP, LensPotFile, factor)
    !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed
    !This input supported by LensPix from 2010
    class(TCLdata) :: this
    Type(CAMBParams), intent(in) :: CP
    integer :: il, unit
    real(dl), intent(in), optional :: factor
    real(dl) fact, scale, BB, TT, TE, EE
    character(LEN=*), intent(in) :: LensPotFile
    !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi and [l(l+1)]^(3/2) C_phi_T/2pi
    !This is the format used by Planck_like but original LensPix uses scalar_output_file.

    !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi
    ! - for historical reasons)

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then
        unit = open_file_header(LensPotFile, 'L', lens_pot_name_tags)
        do il=CP%Min_l, min(10000,CP%Max_l)
            TT = this%Cl_scalar(il, C_Temp)
            EE = this%Cl_scalar(il, C_E)
            TE = this%Cl_scalar(il, C_Cross)
            if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                TT= TT+this%Cl_tensor(il, CT_Temp)
                EE= EE+this%Cl_tensor(il, CT_E)
                TE= TE+this%Cl_tensor(il, CT_Cross)
                BB= this%Cl_tensor(il, CT_B)
            else
                BB=0
            end if
            scale = (real(il+1)/il)**2/OutputDenominator !Factor to go from old l^4 factor to new

            write(unit,'(1I6,7E15.6)') il , fact*TT, fact*EE, fact*BB, fact*TE, scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        do il=10100,CP%Max_l, 100
            scale = (real(il+1)/il)**2/OutputDenominator
            write(unit,'(1E15.5,7E15.6)') real(il), fact*this%Cl_scalar(il,C_Temp:C_E),0.,&
                fact*this%Cl_scalar(il,C_Cross), scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_lens_pot_files


    subroutine TCLdata_output_veccl_files(this, CP,VecFile, factor)
    class(TCLData) :: this
    integer :: il, unit
    Type(CAMBParams), intent(in) :: CP
    character(LEN=*), intent(in) :: VecFile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact

    fact = PresentDefault(1._dl, factor)

    if (CP%WantVectors .and. VecFile /= '') then
        unit = open_file_header(VecFile, 'L', CT_name_tags)
        do il=CP%Min_l,CP%Max_l
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_vector(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if

    end subroutine TCLdata_output_veccl_files

    subroutine TCLdata_NormalizeClsAtL(this,CP,lnorm)
    class(TCLData) :: this
    Type(CAMBParams), intent(in) :: CP
    integer, intent(IN) :: lnorm
    real(dl) Norm

    if (CP%WantScalars) then
        Norm=1/this%Cl_scalar(lnorm, C_Temp)
        this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:C_Cross) = this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:C_Cross) * Norm
    end if

    if (CP%WantTensors) then
        if (.not.CP%WantScalars) Norm = 1/this%Cl_tensor(lnorm, C_Temp)
        !Otherwise Norm already set correctly
        this%Cl_tensor(CP%Min_l:CP%Max_l_tensor, CT_Temp:CT_Cross) =  &
            this%Cl_tensor(CP%Min_l:CP%Max_l_tensor, CT_Temp:CT_Cross) * Norm
    end if

    end subroutine TCLdata_NormalizeClsAtL

    end module results


    module Transfer
    !Module for matter transfer function/matter power spectrum
    use results
    implicit none
    public
    integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
        Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
        Transfer_tot=7, Transfer_nonu=8, Transfer_tot_de=9,  &
        ! total perturbations with and without neutrinos, with neutrinos+dark energy in the numerator
        Transfer_Weyl = 10, & ! the Weyl potential, for lensing and ISW
        Transfer_Newt_vel_cdm=11, Transfer_Newt_vel_baryon=12,   & ! -k v_Newtonian/H
        Transfer_vel_baryon_cdm = 13 !relative velocity of baryons and CDM
    !Sources
    !Alternatively for 21cm
    integer, parameter :: Transfer_monopole=4, Transfer_vnewt=5, Transfer_Tmat = 6

    integer, parameter :: Transfer_max = Transfer_vel_baryon_cdm
    character(LEN=name_tag_len) :: Transfer_name_tags(Transfer_max-1) = &
        ['CDM     ', 'baryon  ', 'photon  ', 'nu      ', 'mass_nu ', 'total   ', &
        'no_nu   ', 'total_de', 'Weyl    ', 'v_CDM   ', 'v_b     ', 'v_b-v_c ']
    character(LEN=name_tag_len) :: Transfer21cm_name_tags(Transfer_max-1) = &
        ['CDM      ', 'baryon   ', 'photon   ', 'monopole ', 'v_newt   ', 'delta_T_g', &
        'no_nu    ', 'total_de ', 'Weyl     ', 'v_CDM    ', 'v_b      ', 'v_b-v_c  ']

    logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
    !set to false to output calculated values for later interpolation

    integer :: transfer_power_var = Transfer_tot
    !What to use to calulcate the output matter power spectrum and sigma_8
    !Transfer_tot uses total matter perturbation

    !Sources
    Type Cl21cmVars
        Type(MatterPowerData), pointer :: PK
        integer l, itf
        logical logs
        real(dl) chi
    end Type Cl21cmVars

    interface Transfer_GetMatterPower
    module procedure Transfer_GetMatterPowerD,Transfer_GetMatterPowerS
    end interface

    contains

    subroutine Transfer_GetUnsplinedPower(State, M, PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k)
    Type(MatterTransferData) :: M
    Type(CAMBdata) :: State
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    real(dl) :: h, k
    integer :: nz, nk, zix, ik, s1, s2
    logical :: hnorm

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)
    hnorm = DefaultTrue (hubble_units)

    nk=M%num_q_trans
    nz=State%CP%Transfer%PK_num_redshifts
    if (nk/= size(PK,1) .or. nz/=size(PK,2)) call MpiStop('Trasfer_GetUnsplinedPower wrong size')
    h = State%CP%H0/100

    if (s1==transfer_power_var .and. s2==transfer_power_var .and. allocated(State%CAMB_Pk)) then
        !Already computed
        do zix=1,nz
            PK(:,zix) = exp(State%CAMB_pk%matpower(:, State%PK_redshifts_index(nz-zix+1)))
        end do
        if (.not. hnorm) PK=  PK / h**3
    else

        do ik=1,nk
            k = M%TransferData(Transfer_kh,ik,1)*h
            do zix=1,nz
                PK(ik,zix) = M%TransferData(s1,ik,State%PK_redshifts_index(nz-zix+1))*&
                    M%TransferData(s2,ik,State%PK_redshifts_index(nz-zix+1))*k*&
                    const_pi*const_twopi*State%CP%InitPower%ScalarPower(k)
            end do
        end do

        if (hnorm) PK=  PK * h**3
    end if

    end subroutine Transfer_GetUnsplinedPower

    subroutine Transfer_GetNonLinRatio_index(State,M, ratio, itf)
    Type(MatterTransferData), intent(in) :: M
    Type(CAMBdata) :: State
    real(dl), allocatable, intent(out) :: ratio(:)
    integer, intent(in) :: itf
    Type(MatterPowerData) :: PKdata

    if (allocated(State%CAMB_PK)) then
        allocate(ratio, source = State%CAMB_PK%nonlin_ratio(:,itf))
    else
        call Transfer_GetMatterPowerData(State, M, PKdata,itf)
        call State%CP%NonLinearModel%GetNonLinRatios(State,PKdata)
        allocate(ratio, source = PKdata%nonlin_ratio(:,1))
    end if
    end subroutine Transfer_GetNonLinRatio_index


    subroutine Transfer_GetUnsplinedNonlinearPower(State,M, PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k) after re-scaling for non-linear evolution (if turned on)
    Type(MatterTransferData), intent(in) :: M
    Type(CAMBdata) :: State
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    integer zix
    real(dl), allocatable :: ratio(:)

    if (.not. allocated(State%CAMB_Pk) .and. State%CP%Transfer%PK_num_redshifts == State%num_transfer_redshifts &
        .and. .not. State%OnlyTransfer) then
        allocate(State%CAMB_PK)
        call Transfer_GetMatterPowerData(State, State%MT, State%CAMB_PK)
        call State%CP%NonLinearModel%GetNonLinRatios(State, State%CAMB_PK)
    end if

    call Transfer_GetUnsplinedPower(State,M,PK,var1,var2, hubble_units)
    do zix=1, State%CP%Transfer%PK_num_redshifts
        call Transfer_GetNonLinRatio_index(State, M, ratio,State%PK_redshifts_index(State%CP%Transfer%PK_num_redshifts-zix+1))
        PK(:,zix) =  PK(:,zix) *ratio**2
    end do

    end subroutine Transfer_GetUnsplinedNonlinearPower

    subroutine Transfer_GetMatterPowerData(State, MTrans, PK_data, itf_only, var1, var2)
    !Does *NOT* include non-linear corrections
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here the definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !spectrum is generated to beyond the CMB k_max
    class(CAMBdata) :: State
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK_data
    integer, intent(in), optional :: itf_only
    integer, intent(in), optional :: var1, var2
    real(dl) :: h, kh, k, power
    integer :: ik, nz, itf, itf_start, itf_end, s1, s2

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

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
    PK_data%redshifts = State%Transfer_Redshifts(itf_start:itf_end)

    h = State%CP%H0/100

    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        k = kh*h
        PK_data%log_kh(ik) = log(kh)
        power = State%CP%InitPower%ScalarPower(k)
        if (global_error_flag/=0) then
            call MatterPowerdata_Free(PK_data)
            return
        end if
        do itf = 1, nz
            PK_data%matpower(ik,itf) = &
                log(MTrans%TransferData(s1,ik,itf_start+itf-1)*&
                MTrans%TransferData(s2,ik,itf_start+itf-1)*k &
                *const_pi*const_twopi*h**3*power)
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
    use FileUtils
    character(LEN=*) :: fname
    type(MatterPowerData) :: PK_data
    type(TTextFile) :: F
    real(dl)kh, Pk
    integer ik
    integer nz

    nz = 1
    call F%Open(fname)

    PK_data%num_k = F%Lines()
    PK_Data%num_z = 1

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))

    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = 0

    do ik=1,PK_data%num_k
        read (F%unit, *) kh, Pk
        PK_data%matpower(ik,1) = log(Pk)
        PK_data%log_kh(ik) = log(kh)
    end do

    call MatterPowerdata_getsplines(PK_data)
    call F%Close()

    end subroutine MatterPowerData_Load


    subroutine MatterPowerdata_getsplines(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    do i = 1,PK_Data%num_z
        call spline_def(PK_data%log_kh,PK_data%matpower(:,i),PK_data%num_k,&
            PK_data%ddmat(:,i))
    end do

    end subroutine MatterPowerdata_getsplines

    !Sources
    subroutine MatterPowerdata_getsplines21cm(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    do i = 1,PK_Data%num_z
        call spline_def(PK_data%log_k,PK_data%matpower(:,i),PK_data%num_k,&
            PK_data%ddmat(:,i))
        call spline_def(PK_data%log_k,PK_data%vvpower(:,i),PK_data%num_k,&
            PK_data%ddvvpower(:,i))
        call spline_def(PK_data%log_k,PK_data%vdpower(:,i),PK_data%num_k,&
            PK_data%ddvdpower(:,i))
    end do

    end subroutine MatterPowerdata_getsplines21cm

    subroutine MatterPowerdata_Free(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    !this shouldn't be needed when releasing the object.
    deallocate(PK_data%log_kh,stat=i)
    deallocate(PK_data%matpower,stat=i)
    deallocate(PK_data%ddmat,stat=i)
    deallocate(PK_data%nonlin_ratio,stat=i)
    deallocate(PK_data%redshifts,stat=i)
    !Sources
    deallocate(PK_data%log_k,stat=i)
    deallocate(PK_data%nonlin_ratio_vv,stat=i)
    deallocate(PK_data%nonlin_ratio_vd,stat=i)
    deallocate(PK_data%vvpower,stat=i)
    deallocate(PK_data%ddvvpower,stat=i)
    deallocate(PK_data%vdpower,stat=i)
    deallocate(PK_data%ddvdpower,stat=i)

    end subroutine MatterPowerdata_Free

    function MatterPowerData_k(PK, kh, itf, index_cache) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: kh
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho,a0,b0
    integer, optional :: index_cache
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
        llo=min(PresentDefault(i_last, index_cache),PK%num_k)
        do while (PK%log_kh(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_kh(llo+1) < logk)
            llo=llo+1
        end do
        if (present(index_cache)) then
            index_cache = llo
        else
            i_last =llo
        end if
        lhi=llo+1
        ho=PK%log_kh(lhi)-PK%log_kh(llo)
        a0=(PK%log_kh(lhi)-logk)/ho
        b0=1-a0

        outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            ((a0**3-a0)* PK%ddmat(llo,itf) &
            +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6
    end if

    outpower = exp(outpower)

    end function MatterPowerData_k

    !Sources
    subroutine MatterPower21cm_k(PK, k, itf, monopole, vv, vd)
    !Get monopole and velocity power at particular k by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: k
    real(dl), intent(out) :: monopole, vv, vd
    real(dl) :: logk
    integer llo,lhi
    real(dl) ho,a0,b0,f1,f2,f3
    integer, save :: i_last = 1

    logk = log(k)
    if (logk < PK%log_k(1)) then
        monopole = 0
        vv=0
        return
    end if
    if (logk > PK%log_k(PK%num_k)) then
        monopole=0
        vv=0
        return
        !stop 'MatterPower21cm_k: out of bounds'
    else
        llo=min(i_last,PK%num_k)
        do while (PK%log_k(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_k(llo+1)< logk)
            llo=llo+1
        end do
        i_last =llo
        lhi=llo+1
        ho=PK%log_k(lhi)-PK%log_k(llo)
        a0=(PK%log_k(lhi)-logk)/ho
        b0=1-a0
        f1= (a0**3-a0)
        f2= (b0**3-b0)
        f3= ho**2/6


        monopole = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            (f1* PK%ddmat(llo,itf) &
            +f2*PK%ddmat(lhi,itf))*f3
        vv = a0*PK%vvpower(llo,itf)+ b0*PK%vvpower(lhi,itf)+&
            (f1* PK%ddvvpower(llo,itf) &
            +f2*PK%ddvvpower(lhi,itf))*f3

        vd = a0*PK%vdpower(llo,itf)+ b0*PK%vdpower(lhi,itf)+&
            (f1* PK%ddvdpower(llo,itf) &
            +f2*PK%ddvdpower(lhi,itf))*f3
    end if

    monopole = exp(max(-30._dl,monopole))
    vv = exp(max(-30._dl,vv))
    vd = exp(max(-30._dl,vd))

    end subroutine MatterPower21cm_k


    subroutine Transfer_GetMatterPowerS(State, MTrans, outpower, itf, minkh, dlnkh, npoints, var1, var2)
    class(CAMBdata) :: state
    Type(MatterTransferData), intent(in) :: MTrans
    integer, intent(in) :: itf, npoints
    integer, intent(in), optional :: var1, var2
    real, intent(out) :: outpower(*)
    real, intent(in) :: minkh, dlnkh
    real(dl) :: outpowerd(npoints)
    real(dl):: minkhd, dlnkhd

    minkhd = minkh; dlnkhd = dlnkh
    call Transfer_GetMatterPowerD(State, MTrans, outpowerd, itf, minkhd, dlnkhd, npoints,var1, var2)
    outpower(1:npoints) = outpowerd(1:npoints)

    end subroutine Transfer_GetMatterPowerS

    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    !Changed input variable from itf to itf_PK because we are looking for the itf_PK'th
    !redshift in the PK_redshifts array.  The position of this redshift in the master redshift
    !array, itf, is given by itf = CP%Transfer%Pk_redshifts_index(itf_PK)
    !Also changed (CP%NonLinear/=NonLinear_None) to
    !CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
    subroutine Transfer_GetMatterPowerD(State, MTrans, outpower, itf_PK, minkh, dlnkh, npoints, var1, var2)
    !Allows for non-smooth priordial spectra
    !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
    !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
    !in units of (h Mpc^{-1})^3.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    class(CAMBdata) :: state
    Type(MatterTransferData), intent(in) :: MTrans

    integer, intent(in) :: itf_PK, npoints
    real(dl), intent(out) :: outpower(npoints)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in), optional :: var1, var2

    integer ik, llo,il,lhi,lastix
    real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
    real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
    integer itf
    integer :: s1,s2, sign
    logical log_interp
    real(dl), allocatable :: ratio(:)

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

    itf = State%PK_redshifts_index(itf_PK)

    if (npoints < 2) call MpiStop('Need at least 2 points in Transfer_GetMatterPower')

    if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
        .and. FeedbackLevel > 0 ) &
        write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'


    if (State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) then
        call Transfer_GetNonLinRatio_index(State, Mtrans, ratio, itf)
    end if

    h = State%CP%H0/100
    logmink = log(minkh)
    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,itf)
        k = kh*h
        kvals(ik) = log(kh)
        atransfer=MTrans%TransferData(s1,ik,itf)*MTrans%TransferData(s2,ik,itf)
        if (State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) &
            atransfer = atransfer* ratio(ik)**2 !only one element, this itf
        matpower(ik) = atransfer*k*const_pi*const_twopi*h**3
        !Put in power spectrum later: transfer functions should be smooth, initial power may not be
    end do
    sign = 1
    log_interp = .true.
    if (any(matpower <= 0)) then
        if (all(matpower < 0)) then
            sign = -1
        else
            log_interp = .false.
        endif
    endif
    if (log_interp) matpower = log(sign*matpower)
    call spline_def(kvals,matpower,MTrans%num_q_trans, ddmat)

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

    if (log_interp) then
        outpower = sign*exp(max(-30.d0,outpower))
    end if
    associate(InitPower => State%CP%InitPower)
        do il = 1, npoints
            k = exp(logmink + dlnkh*(il-1))*h
            outpower(il) = outpower(il) * InitPower%ScalarPower(k)
            if (global_error_flag /= 0) exit
        end do
    end associate

    end subroutine Transfer_GetMatterPowerD

    subroutine Transfer_Get_SigmaR(State, MTrans, R, outvals, var1, var2, root)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc, for all requested redshifts
    !set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R
    integer, intent(in), optional :: var1, var2
    logical, intent(in), optional :: root !if true, give sigma8, otherwise sigma8^2
    real(dl), intent(out) :: outvals(:)
    real(dl) :: kh, k, h, x, win, lnk, dlnk, lnko, powers
    real(dl), dimension(State%CP%Transfer%PK_num_redshifts) :: dsig8, dsig8o, sig8, sig8o
    integer :: s1, s2, ik

    s1 = PresentDefault(transfer_power_var, var1)
    s2 = PresentDefault(transfer_power_var, var2)
    H=State%CP%h0/100._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    do ik=1, MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        if (kh==0) cycle
        k = kh*H

        dsig8 = MTrans%TransferData(s1,ik, State%PK_redshifts_index(1:State%CP%Transfer%PK_num_redshifts))
        if (s1==s2) then
            dsig8 = dsig8**2
        else
            dsig8 = dsig8*MTrans%TransferData(s2,ik, State%PK_redshifts_index(1:State%CP%Transfer%PK_num_redshifts))
        end if
        x= kh *R
        if (x < 1e-2_dl) then
            win = 1._dl - x**2/10
        else
            win = 3*(sin(x)-x*cos(x))/x**3
        end if
        lnk=log(k)
        if (ik==1) then
            dlnk=0.5_dl
            !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
            !Contribution should be very small in any case
        else
            dlnk=lnk-lnko
        end if
        powers = State%CP%InitPower%ScalarPower(k)
        dsig8=(win*k**2)**2*powers*dsig8
        sig8=sig8+(dsig8+dsig8o)*dlnk/2
        dsig8o=dsig8
        lnko=lnk
    end do

    if (present(root)) then
        if (root) sig8 =sqrt(sig8)
    else
        sig8 =sqrt(sig8)
    end if
    outvals(1:State%CP%Transfer%PK_num_redshifts) = sig8

    end subroutine Transfer_Get_SigmaR

    subroutine Transfer_GetSigmaRArray(State, MTrans, R, sigmaR, redshift_ix, var1, var2)
    !Get array of SigmaR at (by default) redshift zero, for all values of R (in h^{-1}Mpc units)
    class(CAMBdata), target :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R(:)
    real(dl), intent(out) :: SigmaR(:)
    integer, intent(in), optional :: redshift_ix, var1, var2
    integer red_ix, ik, subk
    real(dl) kh, k, h, dkh, k_step
    real(dl) lnk, dlnk, lnko, minR, maxR
    real(dl), dimension(size(R)) ::  x, win, dsig8, dsig8o, sig8, sig8o
    type(MatterPowerData), target:: PKspline
    type(MatterPowerData), pointer :: PK
    integer PK_ix
    integer :: nsub
    integer index_cache, nextra

    index_cache = 1
    nsub = 5 !interpolation steps
    h=State%CP%h0/100._dl
    minR = minval(R)/ h
    maxR = maxval(R)
    red_ix = PresentDefault(State%PK_redshifts_index(State%CP%Transfer%PK_num_redshifts), redshift_ix)

    if (allocated(State%CAMB_PK) .and. var1==transfer_power_var .and. var2==transfer_power_var) then
        PK => State%CAMB_PK
        PK_ix = red_ix
    else
        call Transfer_GetMatterPowerData(State, MTrans, PKspline, red_ix, var1, var2)
        PK => PKspline
        PK_ix = 1
    end if
    dkh = 0._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    if (MTrans%TransferData(Transfer_kh,1,1)==0) call MpiStop('Transfer_GetSigmaRArray kh zero')

    !Steps to extrapolate beyond kmax for tail [could do analytically]
    nextra = (4 /minR - State%CP%Transfer%kmax)/h/ &
        (MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,1)- MTrans%TransferData(Transfer_kh,MTrans%num_q_trans-1,1))
    do ik=1, MTrans%num_q_trans + max(2, nextra)
        if (ik < MTrans%num_q_trans) then
            k_step= MTrans%TransferData(Transfer_kh,ik+1,1)-MTrans%TransferData(Transfer_kh,ik,1)
            kh = MTrans%TransferData(Transfer_kh,ik,1)
            nsub = max(1, nint(k_step*min(maxR,4/(kh*h))/0.5))
        else
            !after last step just extrapolate with previous size
            nsub = max(1, nint(k_step*min(maxR,4/(kh*h))))
        end if
        nsub = nint(nsub*State%CP%Accuracy%AccuracyBoost)
        dkh = k_step/nsub
        do subk = 1, nsub
            k = kh*h
            lnk=log(k)

            x= kh *R
            where (x < 1e-2_dl)
                win = 1._dl - x**2/5
            elsewhere
                win =(3*(sin(x)-x*cos(x))/x**3)**2
            end where
            if (ik==1 .and. subk==1) then
                dlnk=0.5_dl
                !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
                !Contribution should be very small in any case
            else
                dlnk=lnk-lnko
            end if
            dsig8=win*(MatterPowerData_k(PK, kh,PK_ix,index_cache)*k**3)
            sig8=sig8+(dsig8+dsig8o)*dlnk/2
            dsig8o=dsig8
            lnko=lnk
            kh = kh + dkh
        end do
    end do

    SigmaR=sqrt(sig8/(const_pi*const_twopi*h**3))

    end subroutine Transfer_GetSigmaRArray

    subroutine Transfer_Get_sigma8(State, MTrans, R, var1, var2)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc
    ! set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var1, var2
    real(dl) :: radius

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)

    call Transfer_Get_SigmaR(State, MTrans, radius, MTrans%sigma_8, var1,var2)

    end subroutine Transfer_Get_sigma8

    subroutine Transfer_Get_sigmas(State, MTrans, R, var_delta, var_v)
    !Get sigma8 and sigma_{delta v} (for growth, like f sigma8 in LCDM)
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var_delta, var_v
    real(dl) :: radius
    integer :: s1, s2

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)
    s1 = PresentDefault (transfer_power_var, var_delta)
    s2 = PresentDefault (Transfer_Newt_vel_cdm, var_v)

    call Transfer_Get_SigmaR(State, MTrans, radius, MTrans%sigma_8, s1,s1)
    if (State%get_growth_sigma8) call Transfer_Get_SigmaR(State, MTrans, radius, &
        MTrans%sigma2_vdelta_8(:), s1, s2, root=.false.)

    end subroutine Transfer_Get_sigmas

    subroutine Transfer_output_Sig8(MTrans, State)
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata), intent(in) :: State
    integer j_PK

    do j_PK=1, State%CP%Transfer%PK_num_redshifts
        write(*,'("at z =",f7.3," sigma8 (all matter) = ",f7.4)') &
            State%CP%Transfer%PK_redshifts(j_PK), MTrans%sigma_8(j_PK)
    end do
    if (State%get_growth_sigma8) then
        do j_PK=1, State%CP%Transfer%PK_num_redshifts
            write(*,'("at z =",f7.3," sigma8^2_vd/sigma8  = ",f7.4)') &
                State%CP%Transfer%PK_redshifts(j_PK), MTrans%sigma2_vdelta_8(j_PK)/MTrans%sigma_8(j_PK)
        end do
    end if

    end subroutine Transfer_output_Sig8


    subroutine Transfer_Allocate(MTrans, State)
    Type(MatterTransferData) :: MTrans
    class(CAMBdata) :: State

    call MTrans%Free()
    allocate(MTrans%q_trans(MTrans%num_q_trans))
    allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans,State%num_transfer_redshifts))
    allocate(MTrans%sigma_8(State%CP%Transfer%PK_num_redshifts))
    if (State%get_growth_sigma8) allocate(MTrans%sigma2_vdelta_8(State%CP%Transfer%PK_num_redshifts))

    end subroutine Transfer_Allocate

    subroutine Transfer_SaveToFiles(MTrans,State,FileNames)
    use constants
    Type(MatterTransferData), intent(in) :: MTrans
    class(CAMBdata) :: State
    integer i,ik
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer i_PK
    integer unit

    do i_PK=1, State%CP%Transfer%PK_num_redshifts
        if (FileNames(i_PK) /= '') then
            i = State%PK_redshifts_index(i_PK)
            if (State%CP%do21cm) then
                unit = open_file_header(FileNames(i_PK), 'k/h', Transfer21cm_name_tags, 14)
            else
                unit = open_file_header(FileNames(i_PK), 'k/h', transfer_name_tags, 14)
            end if
            do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                    write(unit,'(*(E15.6))') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
            end do
            close(unit)
        end if
    end do

    end subroutine Transfer_SaveToFiles

    subroutine Transfer_SaveMatterPower(MTrans, State,FileNames, all21cm)
    use constants
    !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata) :: State
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    character(LEN=name_tag_len) :: columns(3)
    integer itf, i, unit
    integer points
    real, dimension(:,:), allocatable :: outpower
    real minkh,dlnkh
    Type(MatterPowerData) :: PK_data
    integer ncol
    logical, intent(in), optional :: all21cm
    logical all21
    !JD 08/13 Changes in here to PK arrays and variables
    integer itf_PK

    all21 = DefaultFalse(all21cm)
    if (all21) then
        ncol = 3
    else
        ncol = 1
    end if

    do itf=1, State%CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = State%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
                allocate(outpower(points,ncol))

                !Sources
                if (all21) then
                    call Transfer_Get21cmPowerData(MTrans, State, PK_data, itf_PK)
                else
                    call Transfer_GetMatterPowerData(State, MTrans, PK_data, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) then
                        call State%CP%NonLinearModel%GetNonLinRatios(State, PK_data)
                        PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
                        call MatterPowerdata_getsplines(PK_data)
                    end if
                end if

                outpower(:,1) = exp(PK_data%matpower(:,1))
                !Sources
                if (all21) then
                    outpower(:,3) = exp(PK_data%vvpower(:,1))
                    outpower(:,2) = exp(PK_data%vdpower(:,1))

                    outpower(:,1) = outpower(:,1)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,2) = outpower(:,2)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,3) = outpower(:,3)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                end if

                call MatterPowerdata_Free(PK_Data)
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData(Transfer_kh,i,1),outpower(i,:)
                end do
                close(unit)
            else
                if (all21) stop 'Transfer_SaveMatterPower: if output all assume not interpolated'
                minkh = 1e-4
                dlnkh = 0.02
                points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
                !             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
                allocate(outpower(points,1))
                call Transfer_GetMatterPowerS(State, MTrans, outpower(1,1), itf, minkh,dlnkh, points)

                columns(1) = 'P'
                unit = open_file_header(FileNames(itf), 'k/h', columns(:1), 15)

                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1)
                end do
                close(unit)
            end if

            deallocate(outpower)
        end if
    end do

    end subroutine Transfer_SaveMatterPower


    subroutine Transfer_Get21cmPowerData(MTrans, State, PK_data, z_ix)
    !In terms of k, not k/h, and k^3 P_k /2pi rather than P_k
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata) :: State
    Type(MatterPowerData) :: PK_data, PK_cdm
    real(dl) h, k, pow
    integer ik
    integer z_ix,nz

    nz = 1
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%vvpower(PK_data%num_k,nz))
    allocate(PK_data%ddvvpower(PK_data%num_k,nz))
    allocate(PK_data%vdpower(PK_data%num_k,nz))
    allocate(PK_data%ddvdpower(PK_data%num_k,nz))
    allocate(PK_data%log_k(PK_data%num_k))
    allocate(PK_data%redshifts(nz))

    PK_data%redshifts = State%Transfer_Redshifts(z_ix)

    h = State%CP%H0/100

    if (State%CP%NonLinear/=NonLinear_None .and. State%CP%NonLinear/=NonLinear_Lens) then
        if (z_ix>1) stop 'not tested more than one redshift with Nonlinear 21cm'
        call Transfer_GetMatterPowerData(State, MTrans, PK_cdm, z_ix)
        call State%CP%NonLinearModel%GetNonLinRatios_All(State,PK_cdm)
    end if

    do ik=1,MTrans%num_q_trans
        k = MTrans%TransferData(Transfer_kh,ik,z_ix)*h
        PK_data%log_k(ik) = log(k)
        pow = State%CP%InitPower%ScalarPower(k)*1d10

        PK_data%matpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2)**2 * pow)
        PK_data%vvpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)**2 * pow)
        PK_data%vdpower(ik,1) = &
            log( abs((MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)*&
            (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2))* pow)

        if (State%CP%NonLinear/=NonLinear_None) then
            PK_data%matpower(ik,1) = PK_data%matpower(ik,1) + 2*log(PK_cdm%nonlin_ratio(ik,z_ix))
            PK_data%vvpower(ik,1) = PK_data%vvpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vv(ik,z_ix))
            PK_data%vdpower(ik,1) = PK_data%vdpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vd(ik,z_ix))
        end if

    end do

    if (State%CP%NonLinear/=NonLinear_None)  call MatterPowerdata_Free(PK_cdm)


    call MatterPowerdata_getsplines21cm(PK_data)

    end subroutine Transfer_Get21cmPowerData

    function Get21cmCl_l(Vars,kin)
    !Direct integration with j^2/r, etc.
    class(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,k, jlm1
    real(dl) Get21cmCl_l
    real(dl) monopole, vv , vd
    external BJL_EXTERNAL
    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK, k, Vars%itf, monopole, vv, vd)
    call bjl_external(Vars%l, x, jl)
    call bjl_external(Vars%l-1, x, jlm1)
    ddjl = -( 2/x*jlm1-(Vars%l+2)*real(Vars%l+1,dl)/x**2*jl + jl)

    Get21cmCl_l = jl**2*monopole + ddjl**2*vv - 2._dl *ddjl*jl*vd
    if (.not. Vars%logs)  Get21cmCl_l =  Get21cmCl_l / k

    end function Get21cmCl_l


    function Get21cmCl_l_avg(Vars,kin)
    !Asymptotic results where we take <cos^2>=1/2 assuming smooth power spectrum
    class(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,cross,k
    real(dl) Get21cmCl_l_avg
    real(dl) monopole, vv , vd,lphalf
    external BJL_EXTERNAL

    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK, k, Vars%itf, monopole, vv, vd)
    lphalf=Vars%l+0.5_dl

    jl = 1/(2*x**2) /sqrt(1-(lphalf/x)**2)

    !  ddjl = (4/x**4+1)/(2*x**2)
    !
    ddjl = (x**4-2*x**2*lphalf**2+lphalf**4)/(x**4*sqrt(x**2-lphalf**2)*x)/2

    !    cross = (2-x**2)/(2*x**4)

    cross = (-x**2+lphalf**2)/(x**2*sqrt(x**2-lphalf**2)*x)/2

    Get21cmCl_l_avg = jl*monopole + ddjl*vv - 2._dl *cross*vd
    if (.not. Vars%logs)  Get21cmCl_l_avg =  Get21cmCl_l_avg / k

    !       Get21cmCl_l_avg=Get21cmCl_l_avg
    end function Get21cmCl_l_avg


    subroutine Transfer_Get21cmCls(MTrans, State,FileNames)
    use constants
    !Get 21cm C_l from sharp shell, using only monopole source and redshift distortions
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata), target :: State
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer itf,ik, itf_PK
    integer points
    character(LEN=name_tag_len), dimension(3), parameter :: Transfer_21cm_name_tags = &
        ['CL  ','P   ','P_vv']
    Type(MatterPowerData), target ::PK_data
    real(dl)  tol,atol, chi, Cl
    integer l, lastl, unit
    real(dl) k_min, k_max,k, avg_fac
    Type(Cl21cmVars) vars
    Type(CAMBParams), pointer :: CP

    CP=>State%CP

    tol = 1e-5/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    do itf_PK=1, CP%Transfer%PK_num_redshifts
        itf = State%PK_redshifts_index(itf_PK)
        if (FileNames(itf) /= '') then
            chi = State%tau0-State%TimeOfz(CP%Transfer%PK_redshifts(itf_PK))

            points = MTrans%num_q_trans

            lastl=0

            call Transfer_Get21cmPowerData(MTrans, State, PK_data, itf)

            unit = open_file_header(FileNames(itf_PK), 'L', Transfer_21cm_name_tags, 8)

            do ik =1, points-1
                k =exp(PK_data%log_k(ik))
                l=nint(k*chi)
                !This is not an approximation, we are just chosing to sample l at values around (k samples)*chi

                if (l>1 .and. l/= lastl) then
                    lastl=l
                    Vars%l=l
                    Vars%chi = chi
                    Vars%PK => PK_data
                    Vars%itf = 1
                    Cl=0
                    atol = tol
                    avg_fac = 200
                    k_min = max(exp(PK_data%log_k(1)), k*(1-20*CP%Accuracy%AccuracyBoost/chi))
                    k_max = CP%Accuracy%AccuracyBoost*max(k*(1+avg_fac/chi), k*(1._dl+real(l,dl)**(-0.666_dl)))

                    if (k_max*chi < l+10) k_max = (l+10)/chi

                    Vars%logs = .false.
                    if (k_max < exp(PK_data%log_k(points))) then
                        !Integrate bessels properly
                        Cl = Integrate_Romberg(Vars,Get21cmCl_l,k_min,k_max, atol, 25)

                        Vars%logs = .true.
                        k_min = log(k_max)


                        if (l>2e6) then
                            !In baryon damping
                            !                     Vars%logs = .false.
                            !                     atol = tol/10
                            !                     k_min = max(exp(PK_data%log_k(1)), k*(1-10/chi) )
                            !                     k_max = k*(1+100/chi)

                            k_max = min(log(5*k), PK_data%log_k(points))

                        elseif (l>1e4) then
                            Vars%logs = .false.
                            k_min = k_max

                            k_max = min(k*35*CP%Accuracy%AccuracyBoost, exp(PK_data%log_k(points)))
                        else
                            !In white noise regime
                            k_max = min(log(max(0.3_dl,k)*18*CP%Accuracy%AccuracyBoost), PK_data%log_k(points))
                        end if

                        Cl = Cl+Integrate_Romberg(Vars,Get21cmCl_l_avg,k_min,k_max, atol, 25)
                    else
                        k_max = exp(PK_data%log_k(points))
                        Cl = Integrate_Romberg(Vars,Get21cmCl_l,k_min,k_max, atol, 25)
                    end if


                    Cl=exp(-2*State%optical_depths_for21cm(itf_PK))*const_fourpi*Cl* &
                        real(l,dl)*(l+1)/const_twopi/1d10

                    write (unit, '(1I8,3E15.5)') l, Cl, exp(PK_data%matpower(ik,1)/1d10), exp(PK_data%vvpower(ik,1)/1d10)
                end if
            end do

            close(unit)

            call MatterPowerdata_Free(PK_Data)
        end if
    end do

    end subroutine Transfer_Get21cmCls


    end module Transfer
