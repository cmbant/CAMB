    module MathUtils
    use precision
    implicit none

    integer, parameter, private :: MATHUTILS_OMP_VECTOR_THRESHOLD = 256

    interface
    FUNCTION obj_function(obj, x)
    use precision
    class(*) :: obj
    real(dl) :: x, obj_function
    END FUNCTION  obj_function
    end interface

    contains

    elemental function airy_ai_fast(x)
    ! Fast real Airy Ai(x), optimized for < 5e-8 absolute error.
    !
    ! Ai-only version of airy_fast.  Uses the same branch cuts, Ai polynomial
    ! coefficients, and simplified asymptotic fallbacks.
    !
    ! Maximum absolute Ai error is about 4.9e-08 on [-6.5,6.5].  On [-50,-6.5) and
    ! (6.5,25.77], the corresponding asymptotic-branch max absolute Ai errors
    ! are about 2.3e-08 and 7.9e-10.  Relative errors can be large near zeros;
    ! absolute error is the intended accuracy diagnostic.
    implicit none

    real(dl), intent(in) :: x
    real(dl) :: airy_ai_fast

    real(dl), parameter :: xneg_tail = -6.5_dl
    real(dl), parameter :: xneg1 = -4.4_dl
    real(dl), parameter :: xneg2 = -2.09_dl
    real(dl), parameter :: xpos1 = 2.09_dl
    real(dl), parameter :: xpos2 = 2.98_dl
    real(dl), parameter :: xpos_tail = 6.5_dl
    real(dl), parameter :: amaxairy = 25.77_dl

    real(dl), parameter :: sqpii = 5.64189583547756286948e-1_dl
    real(dl), parameter :: half_sqpii = 2.82094791773878143474e-1_dl
    real(dl), parameter :: piq = 7.85398163397448309616e-1_dl
    real(dl), parameter :: twothird = 6.66666666666666666667e-1_dl

    real(dl) :: y, a0
    real(dl) :: q, rt, qtr, zeta, z, zz
    real(dl) :: theta, sn, cs, ak
    real(dl) :: rf, rg, uf, ug

    if (x >= 0.0_dl) then
        if (x < xpos1) then
            y = (x - 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 6.44822251131735552e-5_dl
            a0 = a0*y - 1.40391640977827854e-4_dl
            a0 = a0*y - 5.40241034087955122e-4_dl
            a0 = a0*y + 3.84494287280219777e-3_dl
            a0 = a0*y - 8.25505826923749429e-3_dl
            a0 = a0*y - 6.03341756314770741e-3_dl
            a0 = a0*y + 7.31876777283342883e-2_dl
            a0 = a0*y - 1.59974952689210193e-1_dl
            a0 = a0*y + 1.28267325475824201e-1_dl

            airy_ai_fast = a0
            return
        else if (x <= xpos2) then
            y = (x - 2.53500000000000014e+0_dl) * 2.24719101123595477e+0_dl

            a0 = -1.15727770414950082e-6_dl
            a0 = a0*y + 7.33423129841212084e-5_dl
            a0 = a0*y - 7.09543991974192004e-4_dl
            a0 = a0*y + 3.72265039571619601e-3_dl
            a0 = a0*y - 1.10827615737743822e-2_dl
            a0 = a0*y + 1.48309115523028411e-2_dl

            airy_ai_fast = a0
            return
        else if (x <= xpos_tail) then
            y = (x - 4.74000000000000021e+0_dl) * 5.68181818181818232e-1_dl

            a0 = 4.78267577513511941e-6_dl
            a0 = a0*y - 9.97689952930635397e-5_dl
            a0 = a0*y + 3.89342672538322220e-4_dl
            a0 = a0*y - 8.64063569753577743e-4_dl
            a0 = a0*y + 1.39928174226917976e-3_dl
            a0 = a0*y - 1.68796847845387589e-3_dl
            a0 = a0*y + 1.42992565990869981e-3_dl
            a0 = a0*y - 7.63522116978156422e-4_dl
            a0 = a0*y + 1.94785506580427345e-4_dl

            airy_ai_fast = a0
            return
        else
            if (x > amaxairy) then
                airy_ai_fast = 0.0_dl
                return
            end if

            rt = sqrt(x)
            zeta = twothird * x * rt
            qtr = sqrt(rt)
            z = 1.0_dl / zeta

            a0 = half_sqpii * exp(-zeta) / qtr
            airy_ai_fast = a0 * (1.0_dl - 6.94444444444444444e-2_dl*z)
            return
        end if
    else
        if (x >= xneg2) then
            y = (x + 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 2.63179257460660381e-6_dl
            a0 = a0*y - 4.61163744426384145e-4_dl
            a0 = a0*y + 1.47236778923881778e-3_dl
            a0 = a0*y + 2.55030111258419071e-3_dl
            a0 = a0*y - 2.30877738024996872e-2_dl
            a0 = a0*y + 3.05135902402533278e-2_dl
            a0 = a0*y + 9.89658201868661991e-2_dl
            a0 = a0*y - 3.05531177969466272e-1_dl
            a0 = a0*y + 1.51357051117096963e-2_dl
            a0 = a0*y + 5.35467704429718250e-1_dl

            airy_ai_fast = a0
            return
        else if (x >= xneg1) then
            y = (x + 3.24500000000000011e+0_dl) * 8.65800865800865571e-1_dl

            a0 = -9.32008077263815366e-5_dl
            a0 = a0*y - 2.62795924501147912e-4_dl
            a0 = a0*y + 2.18154044382433019e-3_dl
            a0 = a0*y - 6.67953041092565650e-4_dl
            a0 = a0*y - 2.16530263983360609e-2_dl
            a0 = a0*y + 4.14069894714043721e-2_dl
            a0 = a0*y + 9.39393868134793808e-2_dl
            a0 = a0*y - 3.26514640075362406e-1_dl
            a0 = a0*y - 1.11225097399635667e-1_dl
            a0 = a0*y + 9.06923730177615073e-1_dl
            a0 = a0*y + 5.02404417784273020e-3_dl
            a0 = a0*y - 4.19008474446736834e-1_dl

            airy_ai_fast = a0
            return
        else if (x >= xneg_tail) then
            y = (x + 5.45000000000000018e+0_dl) * 9.52380952380952550e-1_dl

            a0 = -1.05869308595918040e-3_dl
            a0 = a0*y + 1.60649285553263622e-3_dl
            a0 = a0*y + 1.25512673666490319e-2_dl
            a0 = a0*y - 3.17192417854017875e-2_dl
            a0 = a0*y - 6.96908798779071192e-2_dl
            a0 = a0*y + 2.55540735382093598e-1_dl
            a0 = a0*y + 1.77842695165789766e-1_dl
            a0 = a0*y - 8.85680442913919341e-1_dl
            a0 = a0*y - 1.82574102657576925e-1_dl
            a0 = a0*y + 8.96114200938392225e-1_dl
            a0 = a0*y + 6.07711779442129535e-2_dl

            airy_ai_fast = a0
            return
        else
            q = -x
            rt = sqrt(q)

            zeta = twothird * q * rt
            qtr = sqrt(rt)
            ak = sqpii / qtr

            z = 1.0_dl / zeta
            zz = z*z

            rf = -3.71305224980621948e-2_dl + zz*5.54256415602097374e-2_dl
            rg = 6.94432325868511585e-2_dl - zz*3.70901829479361872e-2_dl

            uf = 1.0_dl + zz*rf
            ug = z*rg

            theta = zeta + piq
            sn = sin(theta)
            cs = cos(theta)

            airy_ai_fast = ak*(sn*uf - cs*ug)
            return
        end if
    end if

    end function airy_ai_fast

    elemental subroutine airy_fast(x, ai, aip)
    ! Fast real Airy Ai(x) and Ai'(x)
    ! Optimized for < 5e-8 absolute error on Ai and < 7e-8 on Ai'.
    !
    !
    ! Relative Airy errors can be large near Airy zeros, so absolute error
    ! is the intended accuracy diagnostic for this routine.
    !
    ! Branch joins are endpoint-matched for smoothness.  Measured Ai jumps
    ! at the hot-interval joins are below 1.0e-12, and Ai' jumps are below
    ! 2.0e-14 in double precision.
    !
    ! Outside the hot interval the fallbacks are deliberately simplified for
    ! speed at the existing +/-6.5 cut points: the positive fallback keeps one
    ! correction term in the decaying asymptotic form, while the negative
    ! fallback uses first-degree polynomial corrections to the oscillatory
    ! asymptotic form.  They are not intended to be moved much closer to the
    ! origin without refitting.

    ! Asymptotic-branch accuracy, measured against the high-accuracy Ai
    ! routine on [-50,-6.5) and (6.5,25.77], has max absolute Ai errors about
    ! 2.3e-8 on the negative side and 7.9e-10 on the positive side.  Checked
    ! against scipy.special.airy, the corresponding max absolute Ai' errors
    ! are about 5.9e-8 and 2.4e-9.  The positive decaying branch has worst
    ! relative errors about 2.9e-4 for Ai and 3.3e-4 for Ai' before the
    ! x > 25.77 zero cutoff; beyond the cutoff relative error is 100%, but
    ! absolute error is below 1e-38 at entry.
    implicit none

    real(dl), intent(in)  :: x
    real(dl), intent(out) :: ai, aip

    real(dl), parameter :: xneg_tail = -6.5_dl
    real(dl), parameter :: xneg1 = -4.4_dl
    real(dl), parameter :: xneg2 = -2.09_dl
    real(dl), parameter :: xpos1 = 2.09_dl
    real(dl), parameter :: xpos2 = 2.98_dl
    real(dl), parameter :: xpos_tail = 6.5_dl
    real(dl), parameter :: amaxairy = 25.77_dl

    real(dl), parameter :: sqpii = 5.64189583547756286948e-1_dl
    real(dl), parameter :: half_sqpii = 2.82094791773878143474e-1_dl
    real(dl), parameter :: piq = 7.85398163397448309616e-1_dl
    real(dl), parameter :: twothird = 6.66666666666666666667e-1_dl
    real(dl), parameter :: threehalf = 1.5_dl

    real(dl) :: y, a0, p0
    real(dl) :: q, rt, qtr, zeta, z, zz, invq
    real(dl) :: theta, sn, cs, ak, h, t
    real(dl) :: rf, drf, rg, drg
    real(dl) :: uf, ug, duf_dz, dug_dz
    real(dl) :: dzdx, dakdx

    if (x >= 0.0_dl) then
        if (x < xpos1) then
            y = (x - 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 6.44822251131735552e-5_dl
            a0 = a0*y - 1.40391640977827854e-4_dl
            a0 = a0*y - 5.40241034087955122e-4_dl
            a0 = a0*y + 3.84494287280219777e-3_dl
            a0 = a0*y - 8.25505826923749429e-3_dl
            a0 = a0*y - 6.03341756314770741e-3_dl
            a0 = a0*y + 7.31876777283342883e-2_dl
            a0 = a0*y - 1.59974952689210193e-1_dl
            a0 = a0*y + 1.28267325475824201e-1_dl

            p0 = -7.60118673475642553e-5_dl
            p0 = p0*y + 4.95893375869159003e-4_dl
            p0 = p0*y - 8.17450074227650341e-4_dl
            p0 = p0*y - 3.10763492600551767e-3_dl
            p0 = p0*y + 1.83393323137599658e-2_dl
            p0 = p0*y - 3.15939011570565628e-2_dl
            p0 = p0*y - 1.73137514312311155e-2_dl
            p0 = p0*y + 1.40071151898272322e-1_dl
            p0 = p0*y - 1.53086081995339490e-1_dl

            ai  = a0
            aip = p0
            return
        else if (x <= xpos2) then
            y = (x - 2.53500000000000014e+0_dl) * 2.24719101123595477e+0_dl

            a0 = -1.15727770414950082e-6_dl
            a0 = a0*y + 7.33423129841212084e-5_dl
            a0 = a0*y - 7.09543991974192004e-4_dl
            a0 = a0*y + 3.72265039571619601e-3_dl
            a0 = a0*y - 1.10827615737743822e-2_dl
            a0 = a0*y + 1.48309115523028411e-2_dl

            p0 = -1.16702644245305844e-5_dl
            p0 = p0*y - 1.25042649360228874e-5_dl
            p0 = p0*y + 6.68401309180268409e-4_dl
            p0 = p0*y - 4.78388560690440074e-3_dl
            p0 = p0*y + 1.67303006042858714e-2_dl
            p0 = p0*y - 2.49050323424267611e-2_dl

            ai  = a0
            aip = p0
            return
        else if (x <= xpos_tail) then
            y = (x - 4.74000000000000021e+0_dl) * 5.68181818181818232e-1_dl

            a0 = 4.78267577513511941e-6_dl
            a0 = a0*y - 9.97689952930635397e-5_dl
            a0 = a0*y + 3.89342672538322220e-4_dl
            a0 = a0*y - 8.64063569753577743e-4_dl
            a0 = a0*y + 1.39928174226917976e-3_dl
            a0 = a0*y - 1.68796847845387589e-3_dl
            a0 = a0*y + 1.42992565990869981e-3_dl
            a0 = a0*y - 7.63522116978156422e-4_dl
            a0 = a0*y + 1.94785506580427345e-4_dl

            p0 = -2.76202942968613047e-5_dl
            p0 = p0*y + 4.33736346824519322e-5_dl
            p0 = p0*y + 7.67982862148008753e-5_dl
            p0 = p0*y - 4.69061800312796341e-4_dl
            p0 = p0*y + 1.29204706931427800e-3_dl
            p0 = p0*y - 2.41891622233840692e-3_dl
            p0 = p0*y + 3.18764633720774451e-3_dl
            p0 = p0*y - 2.88254497277541655e-3_dl
            p0 = p0*y + 1.62470672928840122e-3_dl
            p0 = p0*y - 4.33663076753041949e-4_dl

            ai  = a0
            aip = p0
            return
        else
            if (x > amaxairy) then
                ai = 0.0_dl
                aip = 0.0_dl
                return
            end if

            rt = sqrt(x)
            zeta = twothird * x * rt
            qtr = sqrt(rt)
            z = 1.0_dl / zeta

            a0 = half_sqpii * exp(-zeta) / qtr
            ai = a0 * (1.0_dl - 6.94444444444444444e-2_dl*z)
            aip = -a0 * rt * (1.0_dl + 9.72222222222222222e-2_dl*z)
            return
        end if
    else
        if (x >= xneg2) then
            y = (x + 1.04499999999999993e+0_dl) * 9.56937799043062309e-1_dl

            a0 = 2.63179257460660381e-6_dl
            a0 = a0*y - 4.61163744426384145e-4_dl
            a0 = a0*y + 1.47236778923881778e-3_dl
            a0 = a0*y + 2.55030111258419071e-3_dl
            a0 = a0*y - 2.30877738024996872e-2_dl
            a0 = a0*y + 3.05135902402533278e-2_dl
            a0 = a0*y + 9.89658201868661991e-2_dl
            a0 = a0*y - 3.05531177969466272e-1_dl
            a0 = a0*y + 1.51357051117096963e-2_dl
            a0 = a0*y + 5.35467704429718250e-1_dl

            p0 = -5.80598484961455786e-5_dl
            p0 = p0*y + 2.28734744824058274e-4_dl
            p0 = p0*y + 1.62690783405808570e-4_dl
            p0 = p0*y - 3.98697578092560394e-3_dl
            p0 = p0*y + 9.74573721615645122e-3_dl
            p0 = p0*y + 1.49424780652776101e-2_dl
            p0 = p0*y - 1.10428766186618330e-1_dl
            p0 = p0*y + 1.16726560906630195e-1_dl
            p0 = p0*y + 2.84108196353359932e-1_dl
            p0 = p0*y - 5.84744021609954867e-1_dl
            p0 = p0*y + 1.44839531108773359e-2_dl

            ai  = a0
            aip = p0
            return
        else if (x >= xneg1) then
            y = (x + 3.24500000000000011e+0_dl) * 8.65800865800865571e-1_dl

            a0 = -9.32008077263815366e-5_dl
            a0 = a0*y - 2.62795924501147912e-4_dl
            a0 = a0*y + 2.18154044382433019e-3_dl
            a0 = a0*y - 6.67953041092565650e-4_dl
            a0 = a0*y - 2.16530263983360609e-2_dl
            a0 = a0*y + 4.14069894714043721e-2_dl
            a0 = a0*y + 9.39393868134793808e-2_dl
            a0 = a0*y - 3.26514640075362406e-1_dl
            a0 = a0*y - 1.11225097399635667e-1_dl
            a0 = a0*y + 9.06923730177615073e-1_dl
            a0 = a0*y + 5.02404417784273020e-3_dl
            a0 = a0*y - 4.19008474446736834e-1_dl

            p0 = 3.29269486272559364e-4_dl
            p0 = p0*y - 9.16766639018285617e-4_dl
            p0 = p0*y - 2.90896820080877852e-3_dl
            p0 = p0*y + 1.70709844169529452e-2_dl
            p0 = p0*y - 4.27581516784487827e-3_dl
            p0 = p0*y - 1.31291783500191994e-1_dl
            p0 = p0*y + 2.15079574277514368e-1_dl
            p0 = p0*y + 4.06684848837665203e-1_dl
            p0 = p0*y - 1.13080899916065758e+0_dl
            p0 = p0*y - 2.88898694907671671e-1_dl
            p0 = p0*y + 1.57043350541346061e+0_dl
            p0 = p0*y + 4.34982024717498091e-3_dl

            ai  = a0
            aip = p0
            return
        else if (x >= xneg_tail) then
            y = (x + 5.45000000000000018e+0_dl) * 9.52380952380952550e-1_dl

            a0 = -1.05869308595918040e-3_dl
            a0 = a0*y + 1.60649285553263622e-3_dl
            a0 = a0*y + 1.25512673666490319e-2_dl
            a0 = a0*y - 3.17192417854017875e-2_dl
            a0 = a0*y - 6.96908798779071192e-2_dl
            a0 = a0*y + 2.55540735382093598e-1_dl
            a0 = a0*y + 1.77842695165789766e-1_dl
            a0 = a0*y - 8.85680442913919341e-1_dl
            a0 = a0*y - 1.82574102657576925e-1_dl
            a0 = a0*y + 8.96114200938392225e-1_dl
            a0 = a0*y + 6.07711779442129535e-2_dl

            p0 = 7.35022683916569100e-4_dl
            p0 = p0*y + 1.63175788291658046e-4_dl
            p0 = p0*y - 1.19285056616511183e-2_dl
            p0 = p0*y + 1.34376266863632989e-2_dl
            p0 = p0*y + 9.72788044161575216e-2_dl
            p0 = p0*y - 2.11235612423675390e-1_dl
            p0 = p0*y - 3.98856047156416094e-1_dl
            p0 = p0*y + 1.21680249694542830e+0_dl
            p0 = p0*y + 6.77584868086557712e-1_dl
            p0 = p0*y - 2.53051108822859394e+0_dl
            p0 = p0*y - 3.47762968362381963e-1_dl
            p0 = p0*y + 8.53442069032961381e-1_dl

            ai  = a0
            aip = p0
            return
        else
            q = -x
            invq = 1.0_dl / q
            rt = sqrt(q)

            zeta = twothird * q * rt
            qtr = sqrt(rt)
            ak = sqpii / qtr

            z = 1.0_dl / zeta
            zz = z*z

            rf = -3.71305224980621948e-2_dl + zz*5.54256415602097374e-2_dl
            drf = 5.54256415602097374e-2_dl
            rg = 6.94432325868511585e-2_dl - zz*3.70901829479361872e-2_dl
            drg = -3.70901829479361872e-2_dl

            uf = 1.0_dl + zz*rf
            duf_dz = (2.0_dl*z) * (rf + zz*drf)

            ug = z*rg
            dug_dz = rg + (2.0_dl*zz)*drg

            theta = zeta + piq
            sn = sin(theta)
            cs = cos(theta)

            h = sn*uf - cs*ug
            ai = ak*h

            dzdx = threehalf*z*invq
            dakdx = 0.25_dl*ak*invq
            t = (-rt)*(cs*uf + sn*ug) + dzdx*(sn*duf_dz - cs*dug_dz)
            aip = dakdx*h + ak*t
            return
        end if
    end if

    end subroutine airy_fast


    subroutine AiryAiFastArray(ai, x, n)
    implicit none

    integer, intent(in) :: n
    real(dl), intent(out) :: ai(n)
    real(dl), intent(in) :: x(n)
    integer :: i

    if (n >= MATHUTILS_OMP_VECTOR_THRESHOLD) then
        !$OMP parallel do default(shared) private(i) schedule(static)
        do i = 1, n
            ai(i) = airy_ai_fast(x(i))
        end do
        !$OMP end parallel do
    else
        ai = airy_ai_fast(x)
    end if

    end subroutine AiryAiFastArray


    subroutine AiryFastArray(ai, aip, x, n)
    implicit none

    integer, intent(in) :: n
    real(dl), intent(out) :: ai(n), aip(n)
    real(dl), intent(in) :: x(n)
    integer :: i

    if (n >= MATHUTILS_OMP_VECTOR_THRESHOLD) then
        !$OMP parallel do default(shared) private(i) schedule(static)
        do i = 1, n
            call airy_fast(x(i), ai(i), aip(i))
        end do
        !$OMP end parallel do
    else
        call airy_fast(x, ai, aip)
    end if

    end subroutine AiryFastArray



    pure elemental function airy_ai_reference(x) result(ai)
    ! High-precision real Airy Ai(x) reference, adapted from the Cephes Airy routine.
    ! This keeps only Ai(x); Bi and derivative branches from the original routine are omitted.
    implicit none

    real(dl), intent(in) :: x
    real(dl) :: ai

    real(dl), parameter :: amaxairy = 25.77_dl
    real(dl), parameter :: acc = 1.0e-14_dl
    real(dl), parameter :: c1 = 0.35502805388781723926_dl
    real(dl), parameter :: c2 = 0.258819403792806798405_dl
    real(dl), parameter :: sqpii = 5.64189583547756286948e-1_dl
    real(dl), parameter :: pi = 3.141592653589793238462643383279502884_dl
    real(dl), parameter :: an(8) = [ &
        3.46538101525629032477e-1_dl, 1.20075952739645805542e1_dl, &
        7.62796053615234516538e1_dl, 1.68089224934630576269e2_dl, &
        1.59756391350164413639e2_dl, 7.05360906840444183113e1_dl, &
        1.40264691163389668864e1_dl, 9.99999999999999995305e-1_dl]
    real(dl), parameter :: ad(8) = [ &
        5.67594532638770212846e-1_dl, 1.47562562584847203173e1_dl, &
        8.45138970141474626562e1_dl, 1.77318088145400459522e2_dl, &
        1.64234692871529701831e2_dl, 7.14778400825575695274e1_dl, &
        1.40959135607834029598e1_dl, 1.00000000000000000470_dl]
    real(dl), parameter :: afn(9) = [ &
        -1.31696323418331795333e-1_dl, -6.26456544431912369773e-1_dl, &
        -6.93158036036933542233e-1_dl, -2.79779981545119124951e-1_dl, &
        -4.91900132609500318020e-2_dl, -4.06265923594885404393e-3_dl, &
        -1.59276496239262096340e-4_dl, -2.77649108155232920844e-6_dl, &
        -1.67787698489114633780e-8_dl]
    real(dl), parameter :: afd(9) = [ &
        1.33560420706553243746e1_dl, 3.26825032795224613948e1_dl, &
        2.67367040941499554804e1_dl, 9.18707402907259625840_dl, &
        1.47529146771666414581_dl, 1.15687173795188044134e-1_dl, &
        4.40291641615211203805e-3_dl, 7.54720348287414296618e-5_dl, &
        4.51850092970580378464e-7_dl]
    real(dl), parameter :: agn(11) = [ &
        1.97339932091685679179e-2_dl, 3.91103029615688277255e-1_dl, &
        1.06579897599595591108_dl, 9.39169229816650230044e-1_dl, &
        3.51465656105547619242e-1_dl, 6.33888919628925490927e-2_dl, &
        5.85804113048388458567e-3_dl, 2.82851600836737019778e-4_dl, &
        6.98793669997260967291e-6_dl, 8.11789239554389293311e-8_dl, &
        3.41551784765923618484e-10_dl]
    real(dl), parameter :: agd(10) = [ &
        9.30892908077441974853_dl, 1.98352928718312140417e1_dl, &
        1.55646628932864612953e1_dl, 5.47686069422975497931_dl, &
        9.54293611618961883998e-1_dl, 8.64580826352392193095e-2_dl, &
        4.12656523824222607191e-3_dl, 1.01259085116509135510e-4_dl, &
        1.17166733214413521882e-6_dl, 4.91834570062930015649e-9_dl]

    real(dl) :: z, zz, t, f, g, uf, ug, zeta, theta, ak

    if (x > amaxairy) then
        ai = 0.0_dl
        return
    end if

    if (x < -2.09_dl) then
        t = sqrt(-x)
        zeta = -2.0_dl*x*t/3.0_dl
        t = sqrt(t)
        ak = sqpii/t
        z = 1.0_dl/zeta
        zz = z*z
        uf = 1.0_dl + zz*airy_polevl(zz, afn)/airy_p1evl(zz, afd)
        ug = z*airy_polevl(zz, agn)/airy_p1evl(zz, agd)
        theta = zeta + 0.25_dl*pi
        ai = ak*(sin(theta)*uf - cos(theta)*ug)
        return
    end if

    if (x >= 2.09_dl) then
        t = sqrt(x)
        zeta = 2.0_dl*x*t/3.0_dl
        g = exp(zeta)
        t = sqrt(t)
        ak = 2.0_dl*t*g
        z = 1.0_dl/zeta
        f = airy_polevl(z, an)/airy_polevl(z, ad)
        ai = sqpii*f/ak
        return
    end if

    f = 1.0_dl
    g = x
    t = 1.0_dl
    uf = 1.0_dl
    ug = x
    ak = 1.0_dl
    z = x*x*x
    do while (t > acc)
        uf = uf*z
        ak = ak + 1.0_dl
        uf = uf/ak
        ug = ug*z
        ak = ak + 1.0_dl
        ug = ug/ak
        uf = uf/ak
        f = f + uf
        ak = ak + 1.0_dl
        ug = ug/ak
        g = g + ug
        t = abs(uf/f)
    end do

    ai = c1*f - c2*g
    contains

    pure function airy_polevl(x, coef) result(value)
    implicit none

    real(dl), intent(in) :: x
    real(dl), intent(in) :: coef(:)
    real(dl) :: value
    integer :: i

    value = coef(1)
    do i = 2, size(coef)
        value = value*x + coef(i)
    end do

    end function airy_polevl


    pure function airy_p1evl(x, coef) result(value)
    implicit none

    real(dl), intent(in) :: x
    real(dl), intent(in) :: coef(:)
    real(dl) :: value
    integer :: i

    value = 1.0_dl
    do i = 1, size(coef)
        value = value*x + coef(i)
    end do

    end function airy_p1evl


    end function airy_ai_reference



    function Integrate_Romberg(obj, fin, a, b, tol, maxit, minsteps, abs_tol)
    !  Rombint returns the integral from a to b of f(obj,x) using Romberg integration.
    !  The method converges provided that f is continuous in (a,b).
    !  f must be real(dl). The first argument is a class instance.
    !  tol indicates the desired relative accuracy in the integral.

    ! Modified by AL to specify max iterations and minimum number of steps
    ! (min steps useful to stop wrong results on periodic or sharp functions)
    use iso_c_binding
    use MiscUtils
    use config, only : global_error_flag, print_fortran_warnings
    class(*) :: obj
    real(dl), external :: fin !a class function
    procedure(obj_function), pointer :: f
    real(dl), intent(in) :: a,b,tol
    integer, intent(in), optional :: maxit,minsteps
    logical, intent(in), optional :: abs_tol
    integer max_it, min_steps
    real(dl) :: Integrate_Romberg
    integer, parameter :: MAXJ=5
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g(MAXJ+1), g0, g1, fourj
    logical abstol

    !convert the class function (un-type-checked) into correct type to call correctly for class argument
    call C_F_PROCPOINTER(c_funloc(fin), f)
    Integrate_Romberg = -1
    max_it = PresentDefault(25, maxit)
    min_steps = PresentDefault(0, minsteps)
    abstol = DefaultFalse(abs_tol)
    h=0.5d0*(b-a)
    gmax=h*(f(obj,a)+f(obj,b))
    if (global_error_flag /=0) return
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
    do
        i=i+1
        if (i > max_it.or.(i > 5.and.abs(error) < tol) .and. nint > min_steps) exit
        !  Calculate next trapezoidal rule approximation to integral.
        g0=0._dl
        do k=1,nint
            g0=g0+f(obj, a+(k+k-1)*h)
            if (global_error_flag /=0) return
        end do
        g0=0.5d0*g(1)+h*g0
        h=0.5d0*h
        nint=nint+nint
        jmax=min(i,MAXJ)
        fourj=1._dl
        do j=1,jmax
            !  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
        end do
        if (abstol) then
            error=abs(gmax-g0)
        else
            if (abs(g0).gt.tol) then
                error=1._dl-gmax/g0
            else
                error=gmax
            end if
        end if
        gmax=g0
        g(jmax+1)=g0
    end do

    Integrate_Romberg=g0
    if (i > max_it .and. abs(error) > tol .and. print_fortran_warnings)  then
        write(*,*) 'Warning: Integrate_Romberg failed to converge; '
        write (*,*)'integral, error, tol:', Integrate_Romberg,error, tol
    end if

    end function Integrate_Romberg


    subroutine brentq(obj,func,ax,bx,tol,xzero,fzero,iflag,fax,fbx)
    use iso_c_binding

    !>
    !  Find a zero of the function \( f(x) \) in the given interval
    !  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
    !  where \( \epsilon \) is the relative machine precision defined as
    !  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
    !
    !  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
    !
    !#References
    !  * R. P. Brent, "[An algorithm with guaranteed convergence for
    !    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
    !    The Computer Journal, Vol 14, No. 4., 1971.
    !  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
    !    Prentice-Hall, Inc., 1973.
    !
    !# See also
    !  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

    use iso_fortran_env, only: error_unit
    implicit none
    class(*) :: obj
    real(dl), external :: func !a class function f(obj,x)
    procedure(obj_function), pointer :: f

    real(dl),intent(in)              :: ax      !! left endpoint of initial interval
    real(dl),intent(in)              :: bx      !! right endpoint of initial interval
    real(dl),intent(in)              :: tol     !! desired length of the interval of uncertainty of the final result (>=0)
    real(dl),intent(out)             :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(dl),intent(out)             :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)              :: iflag   !! status flag (`-1`=error, `0`=root found)
    real(dl),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(dl),intent(in),optional     :: fbx     !! if `f(bx)` is already known, it can be input here
    real(dl), parameter :: one = 1._dl, zero = 0._dl, two =2._dl, three = 3._dl
    real(dl),parameter :: eps   = epsilon(one)  !! original code had d1mach(4)
    real(dl) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    !convert the class function (un-type-checked) into correct type to call correctly for class argument
    call C_F_PROCPOINTER(c_funloc(func), f)

    tol1 = eps+one

    a=ax
    b=bx

    if (present(fax)) then
        fa = fax
    else
        fa=f(obj,a)
    end if
    if (present(fbx)) then
        fb = fbx
    else
        fb=f(obj,b)
    end if

    !check trivial cases first:
    if (fa==zero) then

        iflag = 0
        xzero = a
        fzero = fa

    elseif (fb==zero) then

        iflag = 0
        xzero = b
        fzero = fb

    elseif (fa*(fb/abs(fb))<zero) then  ! check that f(ax) and f(bx) have different signs

        c=a
        fc=fa
        d=b-a
        e=d

        do

            if (abs(fc)<abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end if

            tol1=two*eps*abs(b)+0.5_dl*tol
            xm = 0.5_dl*(c-b)
            if ((abs(xm)<=tol1).or.(fb==zero)) exit

            ! see if a bisection is forced
            if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
                s=fb/fa
                if (a/=c) then
                    ! inverse quadratic interpolation
                    q=fa/fc
                    r=fb/fc
                    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q=(q-one)*(r-one)*(s-one)
                else
                    ! linear interpolation
                    p=two*xm*s
                    q=one-s
                end if
                if (p<=zero) then
                    p=-p
                else
                    q=-q
                end if
                s=e
                e=d
                if (((two*p)>=(three*xm*q-abs(tol1*q))) .or. &
                    (p>=abs(0.5_dl*s*q))) then
                    d=xm
                    e=d
                else
                    d=p/q
                end if
            else
                d=xm
                e=d
            end if

            a=b
            fa=fb
            if (abs(d)<=tol1) then
                if (xm<=zero) then
                    b=b-tol1
                else
                    b=b+tol1
                end if
            else
                b=b+d
            end if
            fb=f(obj,b)
            if ((fb*(fc/abs(fc)))>zero) then
                c=a
                fc=fa
                d=b-a
                e=d
            end if

        end do

        iflag = 0
        xzero = b
        fzero = fb

    else

        iflag = -1
        write(error_unit,'(A)')&
            'Error in zeroin: f(ax) and f(bx) do not have different signs.'

    end if

    end subroutine brentq

    function Newton_Raphson2(xxl,xxh,funcs, param, param2) result(xm)
    use Precision
    implicit none
    real(dl), intent(in) :: xxl     ! root bracket 1
    real(dl), intent(in) :: xxh     ! root bracket 2
    real(dl)  :: xl,xh, xm     ! root
    external funcs        ! subroutine for non-linear equation
    real(dl), intent(in) :: param, param2 !parameters for function
    integer  :: k                      ! iteration count
    real(dl) :: xn, f,f2,df, error
    real(dl), parameter :: half=0.5_dl
    integer, parameter :: ITERMAX=1000 ! max number of iteration
    real(dl), parameter :: tol=1.e-8_dl ! tolerance for error

    xl =xxl
    xh = xxh
    call funcs(f,df,xl, param, param2)  ! Set xm=f(xl)
    call funcs(f2,df,xh, param, param2)  ! Set xn=f(xh)
    if (f*f2 > 0._dl) then           ! check if function changes sign
        error stop 'Newton_Raphson: root is not bracketed'
    endif
    if (f > 0._dl) then               ! Rearrange so that f(xl)< 0.d0 < f(xh)
        xm = xl
        xl = xh
        xh = xm
    endif

    error = abs(xh-xl)                ! error is width of bracketing interval
    xm = half*(xl+xh)                 ! Initialize guess for root
    k = 0                             ! initialize iteration count
    do while (error > tol .and. k < ITERMAX) ! iterate
        k = k+1                         ! increment iteration count
        call funcs(f,df,xm, param, param2) ! calculate f(xm), df(xm)
        if (f > 0._dl) then              ! Update root bracketing
            xh = xm                       ! update high
        else
            xl = xm                       ! update low
        endif
        xn = xm - f/df                  ! Tentative newton-Raphson step
        if ( (xn-xl)*(xn-xh) > 0._dl ) then ! check if new root falls within bracket
            xm = half* (xh+xl)            ! if no use a Bisection step
            error = abs(xh-xl)            ! error is width of interval
        else
            error = abs(xn-xm)            ! if within bracket: error is change in root
            xm = xn                       ! update successful Newton-Raphson step
        endif
    enddo

    if (error > tol) then       ! Check if solution converged
        write(*,*) 'Newton_Raphson:solution did not converge, xn, funcs(xn),D(xn)'
        write(*,*) xn, f, error
    endif

    end function Newton_Raphson2


    subroutine Gauss_Legendre(x, w, n)
    ! Get n-point Gauss-Legendre nodes x and weights w on [-1, 1].

    use constants, only : dl, const_pi
    implicit none

    integer,  intent(in)  :: n
    real(dl), intent(out) :: x(n), w(n)

    real(dl), parameter :: tol = 1.0e-15_dl
    integer,  parameter :: max_iter = 50

    integer  :: i, j, k, m, iter
    real(dl) :: p1, p2, p3, pp, z, dz, wi

    if (n < 1) error stop "Gauss_Legendre: n must be positive"

    m = (n + 1) / 2

!$omp parallel do default(none) schedule(static) if(n >= 128) &
!$omp& private(i, j, k, iter, p1, p2, p3, pp, z, dz, wi) &
!$omp& shared(x, w, n, m)
    do i = 1, m

        z = cos(const_pi * (real(i, dl) - 0.25_dl) / &
            (real(n, dl) + 0.5_dl))

        do iter = 1, max_iter
            p1 = 1._dl
            p2 = 0._dl

            do j = 1, n
                p3 = p2
                p2 = p1
                p1 = ((2*j - 1) * z * p2 - (j - 1) * p3) / j
            end do

            pp = real(n, dl) * (z*p1 - p2) / (z*z - 1._dl)
            dz = p1 / pp
            z  = z - dz

            if (abs(dz) <= tol) exit
        end do

        if (iter > max_iter) error stop "Gauss_Legendre: Newton iteration failed"

        k  = n + 1 - i
        wi = 2._dl / ((1._dl - z*z) * pp*pp)

        x(i) = -z
        x(k) =  z
        w(i) = wi
        w(k) = wi

        if (i == k) x(i) = 0._dl
    end do
!$omp end parallel do

    end subroutine Gauss_Legendre

    subroutine Legendre_Table(x, P, dP, lmax, npoints)
    ! Legendre polynomials P_l(x_i) and derivatives dP_l(x_i)/dx for all 0 <= l <= lmax
    ! at each of the npoints x values, requiring |x| < 1.
    ! P and dP are (0:lmax, npoints) arrays (C-ordered (npoints, lmax+1) from python).

    use constants, only : dl
    implicit none

    integer, intent(in) :: lmax, npoints
    real(dl), intent(in) :: x(npoints)
    real(dl), intent(out) :: P(0:lmax, npoints), dP(0:lmax, npoints)
    integer :: i, l
    real(dl) :: xi, pmm, pmmp1, pl, rsin2

    !$omp parallel do default(shared) schedule(static) if(npoints >= 8) &
    !$omp& private(i, l, xi, pmm, pmmp1, pl, rsin2)
    do i = 1, npoints
        xi = x(i)
        P(0, i) = 1._dl
        dP(0, i) = 0._dl
        if (lmax < 1) cycle
        P(1, i) = xi
        dP(1, i) = 1._dl
        rsin2 = 1._dl / (1._dl - xi*xi)
        pmm = 1._dl
        pmmp1 = xi
        do l = 2, lmax
            pl = ((2*l - 1)*xi*pmmp1 - (l - 1)*pmm) / l
            P(l, i) = pl
            dP(l, i) = l*(pmmp1 - xi*pl)*rsin2
            pmm = pmmp1
            pmmp1 = pl
        end do
    end do
    !$omp end parallel do

    end subroutine Legendre_Table

    subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input parameters.
    !
    ! Generates the set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3) for all allowed l1
    ! from l1min = max(abs(l2-l3),abs(m1)) to l1max = l2+l3, with m1 = -(m2+m3).
    ! The resulting 3j-coeffs are stored as thrcof(l1-l1min+1).
    !
    ! For numerical stability the recursion proceeds simultaneously forwards from
    ! l1min and backwards from l1max, with the two branches matched over three
    ! overlapping points where the forward recursion stops being stable.
    use MpiUtils, only : MpiStop
    implicit none
    integer, intent(in) :: l2in, l3in, m2in, m3in
    real(dl), intent(out) :: thrcof(*)
    integer, parameter :: i8 = selected_int_kind(18)
    integer(i8) :: l2, l3, m2, m3
    integer(i8) :: l1, m1, l1min, l1max, nfin, a1, a2, dv, dv0, m32

    real(dl) :: newfac, oldfac, sumfor, c1, c2, c1old, denom, x, sum1, sumuni
    real(dl) :: x1, x2, x3, y, y1, y2, y3, sum2, sumbac, ratio, cnorm, thresh
    integer :: i, idx, nlim, lstep, nstep2, n
    logical :: forward_done
    real(dl), parameter :: zero = 0._dl, one = 1._dl
    real(dl), parameter :: tiny_val = 1.0e-30_dl, srtiny = 1.0e-15_dl
    real(dl), parameter :: huge_val = 1.e30_dl, srhuge = 1.e15_dl

    l2 = l2in
    l3 = l3in
    m2 = m2in
    m3 = m3in
    newfac = 0
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    if (l2 < abs(m2) .or. l3 < m3) then
        call MpiStop('GetThreeJs: invalid input (l2 < |m2| or l3 < m3)')
        return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3), abs(m1))
    l1max = l2+l3

    if (l1min >= l1max) then
        if (l1min /= l1max) then
            call MpiStop('GetThreeJs: invalid input (no allowed l1)')
            return
        end if

        ! reached if l1 can take only one value, i.e. l1min = l1max
        thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
        return

    end if

    nfin = l1max-l1min+1
    m32 = m3-m2
    dv0 = (l3*(l3+1) - l2*(l2+1))*m1

    ! forward recursion from l1min
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1+1)*tiny_val
    lstep = 1
    forward_done = .false.

    forward: do
        lstep = lstep+1
        l1 = l1+1

        oldfac = newfac
        a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
        a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
        newfac = sqrt(a2*real(a1,dl))
        if (l1 == 1) then
            ! if l1 = 1, (l1-1) has to be factored out of dv, hence
            c1 = -(2*l1-1)*l1*m32/newfac
        else
            dv = dv0 + l1*(l1-1)*m32
            denom = (l1-1)*newfac

            if (lstep > 2) c1old = abs(c1)
            c1 = -(2*l1-1)*real(dv,dl)/denom
        end if

        if (lstep <= 2) then
            ! if l1 = l1min+1 the third term in the recursion eqn vanishes, hence
            x = srtiny*c1
            thrcof(2) = x
            sum1 = sum1 + tiny_val*(2*l1+1)*c1*c1
            if (lstep == nfin) then
                forward_done = .true.
                exit forward
            end if
            cycle forward
        end if

        c2 = -l1*oldfac/denom

        ! recursion to the next 3j-coeff x
        x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
        thrcof(lstep) = x
        sumfor = sum1
        sum1 = sum1 + (2*l1+1)*x*x
        if (lstep == nfin) exit forward

        ! rescale if the last 3j-coeff exceeds srhuge, to prevent overflow
        if (abs(x) >= srhuge) then
            do i = 1, lstep
                if (abs(thrcof(i)) < srtiny) thrcof(i) = zero
                thrcof(i) = thrcof(i)/srhuge
            end do
            sum1 = sum1/huge_val
            sumfor = sumfor/huge_val
            x = x/srhuge
        end if

        ! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
        ! 3j-values and so is numerically stable. Once an increase of abs(c1) is
        ! detected, the recursion direction is reversed.
        if (c1old <= abs(c1)) exit forward
    end do forward

    if (forward_done) then
        ! forward recursion already covered the whole range (nfin = 2)
        sumuni = sum1
    else
        ! keep three 3j-coeffs around the matching point for comparison with
        ! the backward recursion
        x1 = x
        x2 = thrcof(lstep-1)
        x3 = thrcof(lstep-2)
        nstep2 = nfin-lstep+3

        ! backward recursion from l1max taking nstep2 steps, so that
        ! forward and backward recursion overlap at the three points
        ! l1 = lmatch-1, lmatch, lmatch+1
        l1 = l1max
        thrcof(nfin) = srtiny
        sum2 = tiny_val*(2*l1+1)

        l1 = l1+2
        lstep = 1

        backward: do
            lstep = lstep+1
            l1 = l1-1

            oldfac = newfac
            a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
            a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
            newfac = sqrt(a1*real(a2,dl))

            dv = dv0 + l1*(l1-1)*m32
            denom = l1*newfac
            c1 = -(2*l1-1)*real(dv,dl)/denom

            if (lstep <= 2) then
                ! if l1 = l1max-1 the third term in the recursion vanishes
                y = srtiny*c1
                thrcof(nfin-1) = y
                sumbac = sum2
                sum2 = sum2 + tiny_val*(2*l1-3)*c1*c1
                cycle backward
            end if

            c2 = -(l1-1)*oldfac/denom

            ! recursion to the next 3j-coeff y
            y = c1*thrcof(nfin+2-lstep) + c2*thrcof(nfin+3-lstep)

            if (lstep == nstep2) exit backward

            thrcof(nfin+1-lstep) = y
            sumbac = sum2
            sum2 = sum2 + (2*l1-3)*y*y

            ! rescale if the last 3j-coeff exceeds srhuge, to prevent overflow
            if (abs(y) >= srhuge) then
                do i = 1, lstep
                    idx = nfin-i+1
                    if (abs(thrcof(idx)) < srtiny) thrcof(idx) = zero
                    thrcof(idx) = thrcof(idx)/srhuge
                end do
                sum2 = sum2/huge_val
                sumbac = sumbac/huge_val
            end if
        end do backward

        ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the
        ! corresponding backward recursion values y1, y2, y3, determining the
        ! ratio such that yi = ratio*xi (i=1,2,3) holds with minimal error
        y3 = y
        y2 = thrcof(nfin+2-lstep)
        y1 = thrcof(nfin+3-lstep)

        ratio = (x1*y1 + x2*y2 + x3*y3)/(x1*x1 + x2*x2 + x3*x3)
        nlim = nfin-nstep2+1

        if (abs(ratio) >= 1) then
            thrcof(1:nlim) = ratio*thrcof(1:nlim)
            sumuni = ratio*ratio*sumfor + sumbac
        else
            ratio = 1/ratio
            thrcof(nlim+1:nfin) = ratio*thrcof(nlim+1:nfin)
            sumuni = sumfor + ratio*ratio*sumbac
        end if
    end if

    ! normalise 3j-coeffs; sign convention for the last coeff fixes the overall phase
    cnorm = 1/sqrt(sumuni)
    if (sign(one,thrcof(nfin))*(-1)**abs(l2+m2-l3+m3) <= 0) cnorm = -cnorm

    if (abs(cnorm) >= one) then
        thrcof(1:nfin) = cnorm*thrcof(1:nfin)
    else
        thresh = tiny_val/abs(cnorm)
        do n = 1, int(nfin)
            if (abs(thrcof(n)) < thresh) thrcof(n) = zero
            thrcof(n) = cnorm*thrcof(n)
        end do
    end if

    end subroutine GetThreeJs


    function GetChiSquared(c_inv, Y, n) result(chi2)
    !get dot_product(matmul(C_inv,Y), Y) efficiently assuming c_inv symmetric
    integer, intent(in) :: n
    real(dl), intent(in) :: Y(n)
    real(dl), intent(in) :: c_inv(n,n)
    integer j
    real(dl) ztemp, chi2

    chi2 = 0
    if (n>=512) then
        !$OMP parallel do private(j,ztemp) reduction(+:chi2) schedule(static,16)
        do  j = 1, n
            ztemp= dot_product(Y(j+1:n), c_inv(j+1:n, j))
            chi2=chi2+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
        end do
    else
        do  j = 1, n
            ztemp= dot_product(Y(j+1:n), c_inv(j+1:n, j))
            chi2=chi2+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
        end do
    end if

    end function GetChiSquared

    subroutine integrate_3j(W,lmax_w, n, dopol, M, lmax)
    !$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    !Get coupling matrix, eg for pesudo-CL
    integer, intent(in) :: lmax, lmax_w, n
    real(dl), intent(in) :: W(0:lmax_w,n)
    logical, intent(in) :: dopol
    real(dl), intent(out) :: M(0:lmax,0:lmax, n)
    integer l1, l2, lplus, lminus, thread_ix, ix
    real(dl), allocatable :: threejj0(:,:), threejj2(:,:)

    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()

    allocate(threejj0(0:2*lmax,thread_ix))
    if (dopol) then
        allocate(threejj2(0:2*lmax,thread_ix))
    end if

    !$OMP parallel do private(l1,l2,lminus,lplus,thread_ix,ix), schedule(dynamic)
    do l1 = 0, lmax
        thread_ix =1
        !$ thread_ix = OMP_GET_THREAD_NUM()+1
        do l2 = 0, l1
            lplus =  min(lmax_w,l1+l2)
            lminus = abs(l1-l2)

            call GetThreeJs(threejj0(lminus:,thread_ix),l1,l2,0,0)

            if (dopol .and. l1>=2 .and. l2>=2) then
                !note that lminus is correct, want max(abs(l1-l2),abs(m1)) where m1=0 here
                !(polarization coupling depends on lowest multipoles of the mask)
                call GetThreeJs(threejj2(lminus:,thread_ix),l1,l2,-2,2)
                M(l2,l1,2) = sum(W(lminus:lplus:2,2)*threejj0(lminus:lplus:2,thread_ix) &
                    *threejj2(lminus:lplus:2,thread_ix)) !TE
                M(l2,l1,3) = sum(W(lminus:lplus:2,3)*threejj2(lminus:lplus:2,thread_ix)**2) !EE
                M(l2,l1,4) = sum(W(lminus+1:lplus:2,3)*threejj2(lminus+1:lplus:2,thread_ix)**2) !EB
            end if
            if (n>1 .and. .not. dopol) then
                threejj0(lminus:lplus,thread_ix) = threejj0(lminus:lplus,thread_ix)**2
                do ix=1,n
                    M(l2,l1,ix) = sum(W(lminus:lplus,ix)* threejj0(lminus:lplus,thread_ix))
                end do
            else
                M(l2,l1,1) = sum(W(lminus:lplus,1)* threejj0(lminus:lplus,thread_ix)**2)
            end if
        end do
    end do

    do l1=0, lmax
        do l2 = l1+1,lmax
            M(l2,l1,:) = M(l1,l2,:)
        end do
    end do
    end subroutine integrate_3j

    end module MathUtils
