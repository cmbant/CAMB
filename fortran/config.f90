    module config
    use precision
    use constants, only: const_twopi
    implicit none

    character(LEN=*), parameter :: version = '1.5.5'

    integer :: FeedbackLevel = 0 !if >0 print out useful information about the model

    logical :: output_file_headers = .true.

    logical :: DebugMsgs =.false. !Set to true to view progress and timing

    logical, parameter :: DebugEvolution = .false. !Set to true to do all the evolution for all k

    real(dl) :: DebugParam = 0._dl !not used but read in, useful for parameter-dependent tests

    integer :: ThreadNum = 0
    !If zero assigned automatically, obviously only used if parallelised

    logical ::  do_bispectrum  = .false.
    logical, parameter :: hard_bispectrum = .false. ! e.g. warm inflation where delicate cancellations

    logical, parameter :: full_bessel_integration = .false. !(go into the tails when calculating the sources)

    real(dl), parameter ::  OutputDenominator = const_twopi
    !When using outNone the output is l(l+1)Cl/OutputDenominator

    !Computed output power spectra data

    integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4, C_PhiTemp = 5, C_PhiE=6
    integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4
    integer, parameter :: name_tag_len = 12
    character(LEN=name_tag_len), dimension(C_PhiE), parameter :: C_name_tags = ['TT','EE','TE','PP','TP','EP']
    character(LEN=name_tag_len), dimension(CT_Cross), parameter :: CT_name_tags = ['TT','EE','BB','TE']
    character(LEN=name_tag_len), dimension(7), parameter :: lens_pot_name_tags = ['TT','EE','BB','TE','PP','TP','EP']

    real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code

    real(dl), parameter :: tol=1.0d-4 !Base tolerance for perturbation integrations

    character(LEN=1024) :: highL_unlensed_cl_template = 'HighLExtrapTemplate_lenspotentialCls.dat'
    !fiducial high-accuracy high-L C_L used for making small cosmology-independent numerical corrections
    !to lensing and C_L interpolation. Ideally close to models of interest, but dependence is weak.
    integer, parameter :: lmax_extrap_highl = 8000
    real(dl), allocatable :: highL_CL_template(:,:)

    integer, parameter :: lensed_convolution_margin = 100
    !Number of L less than L max at which the lensed power spectrum is calculated

    integer :: global_error_flag=0

    character(LEN=1024) :: global_error_message = ''

    integer, parameter :: error_reionization=1
    integer, parameter :: error_recombination=2
    integer, parameter :: error_inital_power=3
    integer, parameter :: error_evolution=4
    integer, parameter :: error_unsupported_params=5
    integer, parameter :: error_darkenergy=6
    integer, parameter :: error_ini=7
    integer, parameter :: error_nonlinear=8

    contains

    subroutine GlobalError(message, id)
    use MiscUtils, only : PresentDefault
    character(LEN=*), intent(IN), optional :: message
    integer, intent(in), optional :: id

    global_error_message = PresentDefault('', message)
    if (present(id)) then
        if (id==0) error stop 'Error id must be non-zero'
        global_error_flag=id
    else
        global_error_flag=-1
    end if

    end subroutine GlobalError

    subroutine CheckLoadedHighLTemplate
    use FileUtils
    use MpiUtils
    integer :: L
    real(dl) :: array(7)
    type(TTextFile) :: F
    character(LEN=:), allocatable :: InLine

    if (.not. allocated(highL_CL_template)) then
        allocate(highL_CL_template(1:lmax_extrap_highl, C_Temp:C_Phi))
        highL_CL_template(1,:)=0

        do while (F%ReadNextContentLine(highL_unlensed_cl_template, InLine))
            read(InLine, *) L, array
            if (L>lmax_extrap_highl) exit
            highL_CL_template(L, C_Temp:C_E) =array(1:2)
            highL_CL_template(L, C_Cross) =array(4)
            highL_CL_template(L, C_Phi) =array(5)
        end do
        if (L< lmax_extrap_highl) &
            call MpiStop('CheckLoadedHighLTemplate: template file does not go up to lmax_extrap_highl')
    end if

    end subroutine CheckLoadedHighLTemplate

    end module config