    !First version AL October 2010

    !Calculates local fnl and CMB lensing bispectra
    !CMB lensing bispectra are calculated in the approximation in which the first-order
    !result is used replacing the unlensed with the lensing power spectra
    !This is non-perturbatively correct to about 1%
    !Note the lensing bispectrum only includes the linear potentials, no Rees-Sciama or SZ
    !See Lewis, Challinor & Hanson 2010 for details

    !Note that the primordial local bispectra are unlensed (see arXiv: 0905.4732)

    !Compile with LAPACK and -DFISHER if you want to get the fisher matrix outputs
    !This is disabled by default in order not to require LAPACK

    !Fisher results are with and without the cosmic variance from low L_1

    module Bispectrum
    use ModelParams
    use ModelData
    use InitialPower
    use SpherBessels
    use IniFile
    implicit none

    integer, parameter :: max_bispectrum_deltas = 5, max_bispectrum_fields=3

    Type TBispectrumParams
        logical do_lensing_bispectrum
        logical do_primordial_bispectrum
        integer nfields
        integer Slice_Base_L, ndelta, deltas(max_bispectrum_deltas)
        logical do_parity_odd
        logical DoFisher
        logical export_alpha_beta
        real(dl) FisherNoise, FisherNoisePol, FisherNoiseFwhmArcmin
        character(LEN=Ini_max_string_len)  FullOutputFile
        logical SparseFullOutput
    end Type

    !global parameter for now, only intend for this module to be used interactively for the moment
    Type(TBispectrumParams)  :: BispectrumParams

    Type TBispectrum
        real(dl), pointer :: b(:,:)
    end Type
    Type TCov
        real(dl), pointer :: C(:,:)
    end type TCov

    Type TCov2
        real(dl) :: C(2,2)
    end type TCov2

    Type TCov3
        real(dl) :: C(3,3)
    end type TCov3

    real(dl), allocatable :: dJl(:,:), dddJl(:,:)
    real(dl), parameter :: InternalScale = 1d10
    character(LEN=1024) :: output_root =''
    integer, parameter :: shape_local = 1, shape_warm = 2, shape_warm2 = 3
    integer, parameter :: shape = shape_local

    real(dl), allocatable :: TransferPolFac(:)      !sqrt((l+2)!/(l-2)!)
    contains

    subroutine InitBesselDerivs(CTrans)
    ! j_l' array for interpolation if needed; not needed for local fnl
    Type(ClTransferData) :: CTrans
    integer i,l1,j
    real(dl) Jm, Jp

    if (allocated(dJl)) then
        deallocate(dJL, dddJl)
    end if
    allocate(dJl(BessRanges%npoints,CTrans%ls%l0),dddJl(BessRanges%npoints,CTrans%ls%l0))

    do i=1, CTrans%ls%l0
        !Spline agrees well
        !  call spline_deriv(BessRanges%points,ajl(1,i),ajlpr(1,i),dJl(1,i),BessRanges%npoints)
        !  call spline(BessRanges%points,dJl(1,i),BessRanges%npoints,spl_large,spl_large,dddJl(1,i))

        l1 = CTrans%ls%l(i)
        do j=1, BessRanges%npoints
            call BJL(l1-1,BessRanges%points(j),Jm)
            call BJL(l1+1,BessRanges%points(j),Jp)
            dJl(j,i)=  ( l1*Jm - (l1+1)*Jp)/(2*l1+1)
        end do
        call spline(BessRanges%points,dJl(1,i),BessRanges%npoints,spl_large,spl_large,dddJl(1,i))

    end do

    end subroutine InitBesselDerivs

    subroutine NonGauss_l_r_localOpt(CTrans, ind, indP, res, resP, nfields, r)
    !functions of the form int dk k^2 k^i j_l(kr) Delta_l(k) [P]
    !ind and indP are arrays of required k^i powers
    !res and resP are the results without and with the power spectrum P in the integrand
    Type(ClTransferData) :: CTrans
    integer, intent(in) :: ind(:), indP(:)
    integer :: nfields
    real(dl) res(CTrans%ls%l0,size(ind),nfields), resP(CTrans%ls%l0,size(indP),nfields)
    real(dl), intent(in) :: r
    integer q_ix, j, bes_ix
    integer n, nP, ellmax
    real(dl) xf , J_l, fac, a2, k, dlnk, term, P, kpow, kpowP

    n = size(ind)
    nP =size(indP)
    res=0
    resP = 0
    do q_ix = 1, CTrans%q%npoints
        k = CTrans%q%points(q_ix)
        xf = k*r
        bes_ix=Ranges_indexOf(BessRanges,xf)
        fac=BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
        a2=(BessRanges%points(bes_ix+1)-xf)/fac
        fac=fac**2*a2/6
        dlnk = CTrans%q%dpoints(q_ix) /k
        P = ScalarPower(k, 1)*InternalScale  !!only first index for now

        ellmax = max(xf/(1-xlimfrac), xf + xlimmin) * AccuracyBoost
        kpow =  k**(ind(1)+3)
        kpowP = k**indP(1) * P
        do j=1,CTrans%ls%l0
            if (CTrans%ls%l(j) <= ellmax) then
                J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                    *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac) !cubic spline
                term = CTrans%Delta_p_l_k(1,j,q_ix)*J_l*dlnk
                res(j,1,1) = res(j,1,1) + term * kpow
                resP(j,1,1) = resP(j,1,1) + term * kpowP
                if (nfields>1) then
                    !E pol
                    term = CTrans%Delta_p_l_k(2,j,q_ix)*J_l*dlnk * TransferPolFac(CTrans%ls%l(j))
                    res(j,1,2) = res(j,1,2) + term * kpow
                    resP(j,1,2) = resP(j,1,2) + term * kpowP
                    if (nfields>2) then
                        !lensing potential
                        term = CTrans%Delta_p_l_k(3,j,q_ix)*J_l*dlnk
                        res(j,1,3) = res(j,1,3) + term * kpow
                        resP(j,1,3) = resP(j,1,3) + term * kpowP
                    end if
                end if

            end if
        end do
    end do
    resP = resP * fourpi
    res = res * 2/pi

    end subroutine NonGauss_l_r_localOpt

    subroutine NonGauss_l_r(CTrans, ind, indP,res, resP,nfields, r)
    !functions of the form int dk k^2 k^i j_l(kr) Delta_l(k) [P]
    !ind and indP are arrays of required k^i powers
    !res and resP are the results without and with the power spectrum P in the integrand
    !Output of P scaled by 1d10 (so bispectrum by 1d20)
    Type(ClTransferData) :: CTrans
    integer:: nfields
    integer, intent(in) :: ind(:), indP(:)
    real(dl) res(CTrans%ls%l0,size(ind),nfields), resP(CTrans%ls%l0,size(indP),nfields)
    real(dl), intent(in) :: r
    integer q_ix, j, bes_ix, i
    integer n, nP, ellmax
    real(dl) xf , J_l, fac, a2, k, dlnk, term, P, kpow(size(ind)), kpow2(size(indP))

    if (shape == shape_local) then
        call NonGauss_l_r_localOpt(CTrans, ind, indP,res, resP, nfields, r)
        return
    end if

    n = size(ind)
    nP =size(indP)
    res=0
    resP = 0
    do q_ix = 1, CTrans%q%npoints
        k = CTrans%q%points(q_ix)
        xf = k*r
        bes_ix=Ranges_indexOf(BessRanges,xf)
        fac=BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
        a2=(BessRanges%points(bes_ix+1)-xf)/fac
        fac=fac**2*a2/6
        dlnk = CTrans%q%dpoints(q_ix) /k
        P = ScalarPower(k, 1)*InternalScale  !!only first index for now

        ellmax = max(xf/(1-xlimfrac), xf + xlimmin) * AccuracyBoost
        do i=1,n
            kpow(i)=k**(ind(i)+3)
        end do
        do i=1,nP
            kpow2(i)=k**indP(i) * P
        end do

        do j=1,CTrans%ls%l0
            if (CTrans%ls%l(j) <= ellmax) then

                J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                    *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac) !cubic spline
                !call BJL(CTrans%ls%l(j), xf, J_l)
                term = CTrans%Delta_p_l_k(1,j,q_ix)*J_l*dlnk
                do i=1,n
                    res(j,i,1) = res(j,i,1) + term *kpow(i)
                end do
                do i=1,nP
                    resP(j,i,1) = resP(j,i,1) + term * kpow2(i)
                end do
                !     if (CTrans%ls%l(j)==8) write (1,'(9D20.7)') &
                !      k, xf, real(term * k**3/dlnk), real(term * k**indP(1) * P), &
                !      real(res(j,1)),real(resP(j,1)), J_l, real(term), real(CTrans%Delta_p_l_k(1,j,q_ix))
                if (nfields>1) then
                    !E pol
                    term = CTrans%Delta_p_l_k(2,j,q_ix)*J_l*dlnk* TransferPolFac(CTrans%ls%l(j))
                    do i=1,n
                        res(j,i,2) = res(j,i,2) + term *kpow(i)
                    end do
                    do i=1,nP
                        resP(j,i,2) = resP(j,i,2) + term * kpow2(i)
                    end do
                    if (nfields>2) then
                        !lensing potential
                        term = CTrans%Delta_p_l_k(3,j,q_ix)*J_l*dlnk
                        do i=1,n
                            res(j,i,3) = res(j,i,3) + term *kpow(i)
                        end do
                        do i=1,nP
                            resP(j,i,3) = resP(j,i,3) + term * kpow2(i)
                        end do
                    end if
                end if

            end if
        end do
    end do
    resP = resP * fourpi
    res = res * 2/pi

    end subroutine NonGauss_l_r


    subroutine GetBispectrum(CTrans)
    !Note: may need high maxetak to make sure oscillatory k integrals cancel correctly
    !for accurate alpha(r), beta(r), e.g. 8000; not so important for bispectrum
    !increase accuracy_boost
    use lensing
    use lvalues
    use constants
    use Ranges
#ifdef FISHER
    use MatrixUtils
#endif          
    integer, parameter :: max_bispectra = 2  !fnl, lensing
    Type(ClTransferData) :: CTrans
    Type(Regions) :: TimeStepsNongauss
    integer, allocatable ::  ind(:), indP(:), indPd(:)
    real(dl), allocatable :: res(:,:,:), resP(:,:,:), resPd(:,:)
    real(dl), allocatable :: res_l(:,:,:), resP_l(:,:,:), resPd_l(:,:)
    real(dl) r, term
    Type(TBispectrum), target,allocatable :: Bispectra(:,:,:,:), OddBispectra (:)
    !TTT, TTE, etc; last index is bispectrum kind, default 1=fnl, 2=lensing
    !OddBispectra are parity odd terms like TEB (if do_parity_odd requested for lensing)
    Type(TBispectrum), pointer :: Bispectrum
    !For use in Fisher approximations
    real(dl) test(lmin:CTrans%ls%l(CTrans%ls%l0))
    integer i, j, l1,l2,l3, il1, n,np, npd
    integer min_l, max_l, lmax
    real(dl) tmp, tmp1, tmp2, tmp3
    real(dl) a3j(0:CTrans%ls%l(CTrans%ls%l0)*2+1)
    real(dl) a3j2(0:CTrans%ls%l(CTrans%ls%l0)*2+1,4,2)
    real(dl) CLForLensingIn(4,lmin:CTrans%ls%l(CTrans%ls%l0)),CPhi(3,lmin:CTrans%ls%l(CTrans%ls%l0))
    Type(lSamples) :: SampleL
    real starttime
    real(dl) Bscale
    integer field, field1,field2,field3, bi_ix,bix
    Type(TCov2), allocatable :: CForLensing(:)
    integer nfields, nbispectra, bispectrum_type
    integer :: fnl_bispectrum_ix = 1
    integer :: lens_bispectrum_ix = 2
    character(LEN=256) ::  file_tag = ''
    integer idelta, fileid
    character(LEN=26) :: BispectrumNames(max_bispectra)
    integer :: parities(3), oddix
    integer, parameter :: lmax_lensing_corrT = 300
    !assume C^{T\psi} zero above this for CMB lensing; also neglect lensing contributions to variance
    integer, parameter :: lmax_lensing_corrE = 40 !assume C^{E\psi} zero above this for CMB lensing

    integer, parameter :: first_order_unlensed = 1, first_order_lensed = 2
    integer, parameter :: lens_bispectrum_approx = first_order_lensed
#ifdef FISHER
    Type(TBispectrum), pointer :: Bispectrum2
    real(dl) Cl(4,lmin:CTrans%ls%l(CTrans%ls%l0))
    real(dl)  a3j_00(0:CTrans%ls%l(CTrans%ls%l0)*2+1)
    integer lstart
    real(dl) Noise, NoiseP, bias
    real(dl), allocatable:: fish_contribs(:,:,:)
    real(dl) fish_contribs_sig(lmin:CTrans%ls%l(CTrans%ls%l0))
    real(dl), allocatable :: ifish_contribs(:,:,:,:,:), Fisher(:,:), tmpFisher(:,:),OptimalFisher(:,:)
    real(dl), allocatable :: tmpBigFisher(:,:), Fisher_L1(:,:,:),tmpProjFisher(:,:)
    real(dl), allocatable :: fish_l1(:,:,:,:), fish_L_ij(:,:), fish_L_noise(:,:)
    Type(TBispectrum), target,allocatable :: SqueezedLensingKernel(:,:)
    real(dl) sigma2, xlc, tmpf(3)
    integer  f1,f2,f3, minl2,bigi,bigj, bispectrum_type2, lmaxcuti
    integer sz,corrsize
    Type(TCov), allocatable :: InvC(:)
    integer ix1,ix2
    real(dl) tmpArr(lmin:CTrans%ls%l(CTrans%ls%l0))
#endif

    parities(1)=1  !T
    parities(2)=1  !E
    parities(3)=-1 !B


    if (BispectrumParams%do_primordial_bispectrum) then
        fnl_bispectrum_ix = 1
        nbispectra=1
        BispectrumNames(fnl_bispectrum_ix)='fnl'
    else
        fnl_bispectrum_ix = 0
        nbispectra=0
    end if
    if (BispectrumParams%do_lensing_bispectrum) then
        lens_bispectrum_ix = fnl_bispectrum_ix+1
        nbispectra=nbispectra+1
        BispectrumNames(lens_bispectrum_ix)='lensing'
    end if
    if (nbispectra>max_bispectra) call MpiStop('check max_bispectra')

    if (CP%InitPower%nn>1) call MpiStop('Bispectrum: multiple initial power spectra not supported')

    nfields=BispectrumParams%nfields

    if (lSampleBoost <50) call MpiStop('Bispectrum assumes lSampleBoost=50 (all L sampled)')

    if (lens_bispectrum_approx == first_order_unlensed) file_tag='_unlens'

    lmax = CTrans%ls%l(CTrans%ls%l0)
    if (CP%DoLensing) lmax = lmax_lensed
    SampleL%l0=0
    l1=1
    do
        if (l1<=lmax_lensing_corrE) then
            l1 = l1+1
        else if (l1<120) then
            l1 =l1+nint(7/AccuracyBoost)
        else
            l1 =l1+nint(50/AccuracyBoost)
        end if
        if (l1>lmax) then
            l1 =lmax
        end if
        if (BispectrumParams%Slice_Base_L>0 .and. SampleL%l0>0) then
            !Make sure requested slice base is actually calculated
            if ( BispectrumParams%Slice_Base_L <l1 .and. BispectrumParams%Slice_Base_L>SampleL%l(SampleL%l0)) then
                SampleL%l0= SampleL%l0 + 1
                SampleL%l(SampleL%l0) = BispectrumParams%Slice_Base_L
            end if
        end if
        SampleL%l0= SampleL%l0 + 1
        ! print *,l1
        SampleL%l(SampleL%l0) = l1
        if (l1 == lmax) exit
    end do

    allocate(Bispectra(nfields,nfields,nfields,nbispectra))
    do field1=1,nfields
        do field2=1,nfields
            do field3=1,nfields
                !Only store l2,l3 that are non-zero, array size is approx
                do bispectrum_type=1,nbispectra
                    allocate(Bispectra(field1,field2,field3,bispectrum_type)%b((lmax*(lmax+1))/4,SampleL%l0))
                    Bispectra(field1,field2,field3,bispectrum_type)%b=0
                end do
            end do
        end do
    end do

    if (BispectrumParams%do_lensing_bispectrum) then

        if (.not. CP%DoLensing) call MpiStop('Must turn on lensing to get lensing bispectra')
        print *,'Getting lensing reduced bispectra'

        allocate(CForLensing(lmax))

        CPhi=0
        do i=lmin,lmax
            CPhi(1,i) = Cl_scalar(i,1,C_Phi)/real(i,dl)**4 * InternalScale
            !set correlations to zero where very small to avoid numerical issues
            if (i<=lmax_lensing_corrT) then
                CPhi(2,i) = Cl_scalar(i,1,C_PhiTemp) /real(i,dl)**3 * InternalScale
            end if
            if (i<=lmax_lensing_corrE) then
                CPhi(3,i) = Cl_scalar(i,1,C_PhiE) /real(i,dl)**3 * InternalScale
            end if
            tmp = i*(i+1)/(2*pi)
            CLForLensingIn(:,i) = CL_lensed(i,1,CT_Temp:CT_Cross) * InternalScale/tmp
            ! CForLensing(i)%C=0
            CForLensing(i)%C(1,1)=CLForLensingIn(1,i)
            CForLensing(i)%C(1,2)=CLForLensingIn(4,i)
            CForLensing(i)%C(2,1)=CLForLensingIn(4,i)
            CForLensing(i)%C(2,2)=CLForLensingIn(2,i)
            ! CForLensing(i)%C(3,3)=CL_lensed(i,1,CT_B) * InternalScale/tmp
        end do

#ifdef FISHER             
        allocate(SqueezedLensingKernel(nfields,nfields))
        do field2=1,nfields
            do field3=1,nfields
                allocate(SqueezedLensingKernel(field2,field3)%b((lmax*(lmax+1))/4,SampleL%l0))
                SqueezedLensingKernel(field2,field3)%b=0
            end do
        end do
#endif
        if (DebugMsgs) starttime=GetTestTime()


        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDULE(STATIC,3) &
        !$OMP PRIVATE(il1,l1,l2,l3,max_l,min_l,bix,bi_ix, tmp1,tmp2,tmp3), &
        !$OMP PRIVATE(field1,field2,field3, Bispectrum, a3j,a3j2)

        do il1= 1, SampleL%l0
            l1 = SampleL%l(il1)
            if (l1 > lmax_lensing_corrT) cycle !no exit in parallel loop
            tmp1=l1*(l1+1)
            bi_ix=0
            do l2= max(lmin,l1), lmax
                tmp2=l2*(l2+1)
                min_l = max(abs(l1-l2),l2)
                if (mod(l1+l2+min_l,2)/=0) then
                    min_l = min_l+1
                end if
                max_l = min(lmax,l1+l2)
                bix=bi_ix
                a3j2(:,:,1)=0.5d0
                if (nfields>1) then
                    call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
                    call GetThreeJs(a3j2(max(2,abs(l2-l1)),1,2),l1,l2,2,0)
                    call GetThreeJs(a3j2(max(2,abs(l2-l1)),2,2),l1,l2,0,2)
                    call GetThreeJs(a3j2(max(0,abs(l2-l1)),3,2),l1,l2,2,-2)
                    do l3=min_l,max_l ,2
                        a3j2(l3,:,2) = a3j2(l3,:,2)/a3j(l3)*0.5d0
                    end do
                end if
                do field1=1,nfields
                    do field2=1,nfields
                        do field3=1,nfields
                            Bispectrum=> Bispectra(field1,field2,field3, lens_bispectrum_ix)
                            bi_ix=bix
                            do l3=min_l,max_l ,2
                                bi_ix=bi_ix+1
                                tmp3=l3*(l3+1)
                                !bispectrum is the reduced bispectrum
#ifdef FISHER             
                                if (field1==1) then
                                    SqueezedLensingKernel(field2,field3)%b(bi_ix,il1)=  &
                                        (-tmp2+tmp3+tmp1)*(a3j2(l3,2,field2)*CForLensing(l3)%C(field2,field3)) + &
                                        (-tmp3+tmp1+tmp2)*(a3j2(l3,2,field3)*CForLensing(l2)%C(field3,field2))
                                end if
#endif            

                                Bispectrum%b(bi_ix,il1)=  &
                                    (-tmp1+tmp2+tmp3) *  &
                                    (a3j2(l3,1,field1)*CPhi(1+field2,l2)*CForLensing(l3)%C(field1,field3) + &
                                    a3j2(l3,3,field1)*CPhi(1+field3,l3)*CForLensing(l2)%C(field1,field2) ) + &
                                    (-tmp2+tmp3+tmp1)* &
                                    (a3j2(l3,3,field2)*CPhi(1+field3,l3)*CForLensing(l1)%C(field2,field1) + &
                                    a3j2(l3,2,field2)*CPhi(1+field1,l1)*CForLensing(l3)%C(field2,field3) ) + &
                                    (-tmp3+tmp1+tmp2)* &
                                    (a3j2(l3,2,field3)*CPhi(1+field1,l1)*CForLensing(l2)%C(field3,field2) + &
                                    a3j2(l3,1,field3)*CPhi(1+field2,l2)*CForLensing(l1)%C(field3,field1) )

                            end do

                        end do
                    end do
                end do

            end do
        end do
        !$OMP END PARAllEl DO
        if (DebugMsgs) print *,'Time for lensing:', GetTestTime()-starttime

        if (nfields==1) BispectrumParams%do_parity_odd=.false.

        if (BispectrumParams%do_parity_odd) then

            allocate(OddBispectra(12))
            oddix=0
            do field1=1,3
                do field2=1,3
                    do field3=1,3
                        if (parities(field1)+parities(field2)+parities(field3)/=1) cycle
                        oddix=oddix+1
                        !Only store l2,l3 that are non-zero, array size is approx
                        allocate(OddBispectra(oddix)%b((lmax*(lmax+1))/4,SampleL%l0))
                        OddBispectra(oddix)%b=0
                    end do
                end do
            end do

            !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDULE(STATIC,3) &
            !$OMP PRIVATE(il1,l1,l2,l3,max_l,min_l,bix,bi_ix, tmp1,tmp2,tmp3), &
            !$OMP PRIVATE(field1,field2,field3, Bispectrum, a3j2, oddix)

            do il1= 1, SampleL%l0
                l1 = SampleL%l(il1)
                if (l1 > lmax_lensing_corrT) cycle !no exit in parallel loop
                tmp1=l1*(l1+1)
                bi_ix=0
                do l2= max(lmin,l1), lmax
                    tmp2=l2*(l2+1)
                    min_l = max(abs(l1-l2),l2)
                    if (mod(l1+l2+min_l,2)/=1) then
                        min_l = min_l+1
                    end if
                    max_l = min(lmax,l1+l2)
                    bix=bi_ix

                    a3j2(:,:,1)=0.5d0
                    call GetThreeJs(a3j2(max(2,abs(l2-l1)),1,2),l1,l2,2,0)
                    call GetThreeJs(a3j2(max(2,abs(l2-l1)),2,2),l1,l2,0,2)
                    call GetThreeJs(a3j2(max(0,abs(l2-l1)),3,2),l1,l2,2,-2)
                    do l3=min_l,max_l ,2
                        a3j2(l3,:,2) = a3j2(l3,:,2)*0.5d0*sqrt(real((2*L1+1)*(2*L2+1),dl)*(2*L3+1)/(3.1415926535*4))
                    end do


                    oddix=0
                    do field1=1,3
                        do field2=1,3
                            do field3=1,3
                                !Only calculate terms with one B
                                if (parities(field1)+parities(field2)+parities(field3)/=1) cycle
                                oddix=oddix+1
                                Bispectrum=> OddBispectra(oddix)
                                bi_ix=bix
                                do l3=min_l,max_l ,2
                                    bi_ix=bi_ix+1

                                    tmp3=l3*(l3+1)
                                    !bispectrum the non-reduced bispectrum without the i
                                    if (parities(field1)==-1) then
                                        Bispectrum%b(bi_ix,il1)=  &
                                            (-tmp1+tmp2+tmp3) *  &
                                            (a3j2(l3,1,min(2,field1))*CPhi(1+field2,l2)*CForLensing(l3)%C(2,field3)  &
                                            -a3j2(l3,3,min(2,field1))*CPhi(1+field3,l3)*CForLensing(l2)%C(2,field2) )
                                    elseif (parities(field2)==-1) then
                                        Bispectrum%b(bi_ix,il1)=  &
                                            (-tmp2+tmp3+tmp1)* &
                                            (-a3j2(l3,3,min(2,field2))*CPhi(1+field3,l3)*CForLensing(l1)%C(2,field1)  &
                                            -a3j2(l3,2,min(2,field2))*CPhi(1+field1,l1)*CForLensing(l3)%C(2,field3) )
                                    else if (parities(field3)==-1) then
                                        Bispectrum%b(bi_ix,il1)=  &
                                            (-tmp3+tmp1+tmp2)* &
                                            (-a3j2(l3,2,min(2,field3))*CPhi(1+field1,l1)*CForLensing(l2)%C(2,field2) + &
                                            a3j2(l3,1,min(2,field3))*CPhi(1+field2,l2)*CForLensing(l1)%C(2,field1) )
                                    end if

                                end do

                            end do
                        end do
                    end do

                end do
            end do
            !$OMP END PARAllEl DO
        end if

    end if

    if (BispectrumParams%do_primordial_bispectrum) then

        print *,'getting reduced local fnl bispectra'

        allocate(TransferPolFac(lmax))
        do i=2,lmax
            TransferPolFac(i) =sqrt( real((i+1)*i,dl)*(i+2)*(i-1))
        end do

        if (shape /= shape_local) call MpiStop('Non-local shapes not working')

        if (shape == shape_local) then
            n=1
            np=1
            npd=0 !derivatives of function
        else if (shape == shape_warm) then
            n=2
            np=3
            npd=0
        else if (shape == shape_warm2) then
            n=1
            np=2
            npd=2
        else
            call MpiStop('unknown shape')
        end if

        allocate(ind(n))
        allocate(indP(np))

        if (npd>0) then
            call InitBesselDerivs(CTrans)
            allocate(indPd(npd))
        end if

        if (shape==shape_warm) then
            !Separable form is very unstable and unworkable probably
            ind(1) = 0
            ind(2) = 2
            indP(1) = 0
            indP(2) = 2
            indP(3) = -2
        else if (shape==shape_warm2) then
            ind(1) = 0
            indP(1) = 0
            indP(2) = -2
            indPd(1) = 0
            indPd(2) = -2
        else
            ind(1) = 0
            indP(1) = 0
        end if

        test=0
        call Ranges_Nullify(TimeStepsNongauss)
        call Ranges_Assign(TimeStepsNongauss,TimeSteps)
        call Ranges_Add_delta(TimeStepsNongauss, -taurst*10*AccuracyBoost, taurst, dtaurec)
        call Ranges_getArray(TimeStepsNongauss, .true.)

        !$        if (BispectrumParams%export_alpha_beta) call OMP_SET_NUM_THREADS(1)
        if (BispectrumParams%export_alpha_beta) then
            !Note that all the points outside recombination are not really needed
            !And these are for curvature perturbation, so do not include 3/5 factor
            call CreateTxtFile(trim(output_root)//'_alpha.txt',100)
            call CreateTxtFile(trim(output_root)//'_beta.txt',101)
            call CreateTxtFile(trim(output_root)//'_alpha_beta_r.txt',102)
        end if

        if (DebugMsgs) starttime=GetTestTime()

        !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(STATIC,3) &
        !$OMP PRIVATE(i,r,res,resP,resPd,res_l,resP_l,resPd_l,term,j), &
        !$OMP PRIVATE(il1,l1,l2,l3,min_l,max_l,tmp,tmp1,tmp2,Bispectrum), &
        !$OMP PRIVATE(bi_ix,bix,field1,field2,field3,field)

        do i= TimeStepsNongauss%npoints-1, 2,-1
            r=(CP%tau0-TimeStepsNongauss%points(i))

            allocate(res(CTrans%ls%l0,n,nfields))
            allocate(resP(CTrans%ls%l0,np,nfields))

            allocate(res_l(1:CTrans%ls%l(CTrans%ls%l0),n,nfields))
            allocate(resP_l(1:CTrans%ls%l(CTrans%ls%l0),np,nfields))
            if (npd>0) then
                allocate(resPd(CTrans%ls%l0,npd))
                allocate(resPd_l(1:CTrans%ls%l(CTrans%ls%l0),npd))
            end if

            call NonGauss_l_r(CTrans, ind, indP,res, resP, nfields, r)
            if (npd>0) call NonGauss_deriv_l_r(CTrans, indPd,resPd, r, dJl,dddJl)

            do field=1,nfields
                do j=1,n
                    call InterpolateClArr(CTransScal%ls,res(1,j,field),res_l(lmin,j,field),CTransScal%ls%l0)
                end do
                do j=1,np
                    call InterpolateClArr(CTransScal%ls,resP(1,j,field),resP_l(lmin,j,field),CTransScal%ls%l0)
                end do
            end do
            deallocate(res,resP)

            if (BispectrumParams%export_alpha_beta) then
                write(100,concat('(',lmax-lmin+1 ,'E15.5)')) res_l(lmin:lmax,1,1)
                write(101,concat('(',lmax-lmin+1 ,'E15.5)')) resP_l(lmin:lmax,1,1)
                write(102,'(1E15.5)') r
            end if

            if (npd>0) then
                do j=1,npd
                    call InterpolateClArr(CTransScal%ls,resPd(1,j),resPd_l(lmin,j),CTransScal%ls%l0)
                end do
                deallocate(resPd)
            end if

            term = r**2 * TimeStepsNongauss%dpoints(i) * (3./5)


            !Restrict to l1<=l2<=l3
            do il1= 1, SampleL%l0
                l1 = SampleL%l(il1)
                bi_ix=0
                do l2= max(lmin,l1), lmax
                    min_l = max(abs(l1-l2),l2)
                    if (mod(l1+l2+min_l,2)/=0) then
                        min_l = min_l+1
                    end if
                    max_l = min(lmax,l1+l2)
                    do field1=1,nfields
                        do field2=1,nfields
                            tmp1 = 2*term*(res_l(l1,1,field1)*resP_l(l2,1,field2) + &
                                res_l(l2,1,field2)*resP_l(l1,1,field1))
                            tmp2 = 2*term*(resP_l(l1,1,field1)*resP_l(l2,1,field2))
                            do field3=1,nfields
                                Bispectrum => Bispectra(field1,field2,field3,fnl_bispectrum_ix)
                                bix=bi_ix
                                do l3=min_l,max_l ,2
                                    bix=bix+1
                                    Bispectrum%b(bix,il1) = Bispectrum%b(bix,il1) + &
                                        (tmp1*resP_l(l3,1,field3) +   tmp2*res_l(l3,1,field3))
                                end do
                            end do
                        end do
                    end do
                    bi_ix=bix

                end do !l2
            end do !il1

            deallocate(res_l,resP_l)
            if (npd>0) deallocate(resPd_l)
        end do !TimeStepsNongauss
        !$OMP END PARAllEl DO
        if (BispectrumParams%export_alpha_beta) then
            close(100)
            close(101)
            close(102)
        end if
        deallocate(TransferPolFac)
        call Ranges_Free(TimeStepsNongauss)

        if (DebugMsgs) print *,'Time for fnl bispectrum:', GetTestTime()-starttime

    end if !DoPrimordial

    if (BispectrumParams%Slice_Base_L>0 .or. BispectrumParams%FullOutputFile/='') then
        !write out slice in (muK)^3 units
        Bscale=(COBE_CMBTemp*1d6)**3/InternalScale**2;
        do bispectrum_type=1,nbispectra
            if (BispectrumParams%Slice_Base_L>0) then
                do idelta=1,BispectrumParams%ndelta
                    if (mod(BispectrumParams%Slice_Base_L + BispectrumParams%deltas(idelta),2)==1 &
                        .and. bispectrum_type/=lens_bispectrum_ix) cycle
                    call CreateTxtFile(concat(trim(output_root)//'bispectrum_'//trim(BispectrumNames(bispectrum_type))//'_base_', &
                        BispectrumParams%Slice_Base_L,'_delta_',BispectrumParams%deltas(idelta),trim(file_tag)//'.dat'),&
                        nbispectra +BispectrumParams%ndelta*(bispectrum_type-1)+idelta)
                end do
            end if
            if (BispectrumParams%FullOutputFile/='') then
                call CreateTxtFile(concat(output_root,BispectrumParams%FullOutputFile, &
                    '_', BispectrumNames(bispectrum_type), file_tag, '.dat'),bispectrum_type)
            end if
        end do
        do il1= 1, SampleL%l0
            l1 = SampleL%l(il1)
            bi_ix=0
            do l2= max(lmin,l1), lmax
                min_l = max(abs(l1-l2),l2)
                if (mod(l1+l2+min_l,2)/=0) then
                    min_l = min_l+1
                end if
                max_l = min(lmax,l1+l2)
                do l3=min_l, max_l ,2
                    bi_ix=bi_ix+1
                    if (l1==BispectrumParams%Slice_Base_L &
                        .and. any(l3-l2==BispectrumParams%deltas(1:BispectrumParams%ndelta))) then
                        !Particular slice
                        idelta=IndexOf(l3-l2,BispectrumParams%deltas,BispectrumParams%ndelta)
                        do bispectrum_type=1,nbispectra
                            fileid=nbispectra +BispectrumParams%ndelta*(bispectrum_type-1)+idelta
                            write (fileid,'(1I5)', advance='NO') L2
                            do field1=1,nfields
                                do field2=1,nfields
                                    do field3=1,nfields
                                        write(fileid,'(1E15.5)', advance='NO') &
                                            Bispectra(field1,field2,field3,bispectrum_type)%b(bi_ix,il1)*Bscale
                                    end do
                                end do
                            end do
                            write (fileid,'(a)') ''
                        end do
                    end if !slice
                    if (BispectrumParams%FullOutputFile/='') then
                        if (BispectrumParams%SparseFullOutput .and. .not. any( SampleL%l(1:SampleL%l0)==L2) .or. &
                            l1 > 30 .and. mod(l3-min_l,10)/=0 .and. l3 /= max_l) cycle

                        do bispectrum_type=1,nbispectra
                            if (bispectrum_type==lens_bispectrum_ix .and. L1 > lmax_lensing_corrT) cycle
                            write(bispectrum_type,'(3I5)', advance='NO') L1, L2, L3
                            do field1=1,nfields
                                do field2=1,nfields
                                    do field3=1,nfields
                                        write(bispectrum_type,'(1E14.5)', advance='NO') &
                                            Bispectra(field1,field2,field3,bispectrum_type)%b(bi_ix,il1)*Bscale
                                    end do
                                end do
                            end do
                            write (bispectrum_type,'(a)') ''
                        end do
                    end if

                end do
            end do
        end do
        if (BispectrumParams%do_parity_odd) then
            do il1= 1, SampleL%l0
                l1 = SampleL%l(il1)
                bi_ix=0
                do l2= max(lmin,l1), lmax
                    min_l = max(abs(l1-l2),l2)
                    if (mod(l1+l2+min_l,2)/=1) then
                        min_l = min_l+1
                    end if
                    max_l = min(lmax,l1+l2)
                    do l3=min_l, max_l ,2
                        bi_ix=bi_ix+1
                        if (l1==BispectrumParams%Slice_Base_L &
                            .and. any(l3-l2==BispectrumParams%deltas(1:BispectrumParams%ndelta))) then
                            !Particular slice
                            idelta=IndexOf(l3-l2,BispectrumParams%deltas,BispectrumParams%ndelta)
                            do bispectrum_type=1,nbispectra
                                if (bispectrum_type/=lens_bispectrum_ix) cycle
                                fileid=nbispectra +BispectrumParams%ndelta*(bispectrum_type-1)+idelta
                                write (fileid,'(1I5)', advance='NO') L2
                                oddix=0
                                do field1=1,3
                                    do field2=1,3
                                        do field3=1,3
                                            if (parities(field1)+parities(field2)+parities(field3)/=1) cycle
                                            oddix=oddix+1
                                            write(fileid,'(1E15.5)', advance='NO') OddBispectra(oddix)%b(bi_ix,il1)*Bscale
                                        end do
                                    end do
                                end do
                                write (fileid,'(a)') ''
                            end do
                        end if !slice

                    end do
                end do
            end do
        end if

        do bispectrum_type=1,nbispectra
            if (BispectrumParams%Slice_Base_L>0) then
                do idelta=1,BispectrumParams%ndelta
                    if (mod(BispectrumParams%Slice_Base_L + BispectrumParams%deltas(idelta),2)==1 &
                        .and. bispectrum_type/=lens_bispectrum_ix) cycle
                    close(nbispectra +BispectrumParams%ndelta*(bispectrum_type-1)+idelta)
                end do
            end if
            if (BispectrumParams%FullOutputFile/='') close(bispectrum_type)
        end do

    end if

#ifdef FISHER
    if (BispectrumParams%DoFisher) then
        !Get stuff for Fisher etc.

        print *,'Getting Fisher for lmax = ', lmax

        Noise = BispectrumParams%FisherNoise/ (COBE_CMBTemp*1e6)**2  !Planckish, dimensionless units
        NoiseP = BispectrumParams%FisherNoisePol/ (COBE_CMBTemp*1e6)**2

        do i=lmin,lmax
            if (CP%DoLensing) then
                cl(:,i) = CL_lensed(i,1,CT_Temp:CT_Cross)
            else
                cl(1,i) = CL_Scalar(i,1,C_Temp)
                cl(2,i) = CL_Scalar(i,1,C_E)
                cl(4,i) = CL_Scalar(i,1,C_Cross)
                cl(3,i) = 0
            end if
            if (CP%WantTensors .and. i<= CP%Max_l_tensor .and. i>=2) then
                cl(:,i) = cl(:,i) + Cl_tensor(i,1,CT_Temp:CT_Cross)
            end if
        end do
        if (.false.) then
            call OpenTxtFile('CAMBdefault_lensedCls.dat',3)
            do i=lmin,lmax
                !Assume T,E,B,X ordering
                read(3,*) j, cl(1:4,i)
                if (j<lmin) read(3,*) j, cl(1:4,i)
                cl(:,i)=cl(:,i)/(COBE_CMBTemp*1e6)**2
            end do
            close(3)
        end if

        if (Noise >0) then
            file_tag = concat(file_tag,'_noise')
        end if
        xlc= 180*sqrt(8.*log(2.))/3.14159
        sigma2 = (BispectrumParams%FisherNoiseFwhmArcmin/60/xlc)**2
        allocate(InvC(lmax))
        do l1= lmin, lmax
            tmp = l1*(l1+1)/(2*pi)
            Cl(1,l1) = Cl(1,l1)/tmp + Noise*exp(l1*(l1+1)*sigma2)
            Cl(2:3,l1) = Cl(2:3,l1)/tmp + NoiseP*exp(l1*(l1+1)*sigma2)
            Cl(4,l1) = Cl(4,l1)/tmp
            allocate(InvC(l1)%C(nfields,nfields))
            if (nfields > 2) call MpiStop('Not implemented nfields>2 in detail')
            if (nfields==1) then
                InvC(l1)%C(1,1)=(2*l1+1)/cl(1,l1)/InternalScale
            else
                InvC(l1)%C(1,1)=cl(2,l1)
                InvC(l1)%C(1,2)=-cl(4,l1)
                InvC(l1)%C(2,1)=-cl(4,l1)
                InvC(l1)%C(2,2)=cl(1,l1)
                InvC(l1)%C= InvC(l1)%C * (2*l1+1)/(cl(1,l1)*cl(2,l1)-cl(4,l1)**2)/InternalScale
            end if
        end do

        if (debugMsgs) starttime=GetTestTime()
        allocate(ifish_contribs(SampleL%l0,nbispectra,nbispectra,nfields,nfields) )
        !This loop is just in case want to plot out lmax dependence
        do lmaxcuti=SampleL%l0, SampleL%l0
            !    call CreateTxtFile('TE-Planck-LensFish.txt',20)
            !   do lmaxcuti=1, SampleL%l0
            !   if (SampleL%l(lmaxcuti) < 425) cycle

            lmax= SampleL%l(lmaxcuti)

            ifish_contribs=0
            lstart = 2 !lmin
            !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDULE(STATIC,3) &
            !$OMP PRIVATE(il1,l1,l2,l3,fish_l1,bi_ix,min_l,max_l,a3j_00,a3j), &
            !$OMP PRIVATE(Bispectrum,Bispectrum2,minl2,bix,tmp,tmp1,tmp2,tmpf), &
            !$OMP PRIVATE(field1,field2,field3,f1,f2,f3,bispectrum_type,bispectrum_type2)

            do il1= 1,  lmaxcuti !!!SampleL%l0
                allocate(fish_l1(nbispectra,nbispectra,nfields,nfields)) !last indices are field1,f1
                l1 = SampleL%l(il1)
                if (l1< lstart) cycle
                fish_l1=0
                bi_ix=0
                do l2 = l1,lmax
                    if (l2< lstart) cycle
                    min_l = max(lstart,max(abs(l1-l2),l2))
                    if (mod(l1+l2+min_l,2)/=0) then
                        min_l = min_l+1
                    end if
                    max_l = min(lmax,l1+l2)
                    call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
                    do l3=min_l,max_l ,2
                        a3j_00(l3)=a3j(l3)**2
                    end do

                    tmp1= 1.d0/(4*pi)  !(2l+1) factors included in InvC
                    minl2=min_l
                    bix=bi_ix
                    do field1=1,nfields
                        do f1=1,nfields
                            tmpf(1)= InvC(l1)%C(field1,f1)*tmp1
                            do field2=1,nfields
                                do f2=1,nfields
                                    tmpf(2)= InvC(l2)%C(field2,f2)*tmpf(1)
                                    do field3=1,nfields
                                        do bispectrum_type=1,nbispectra
                                            if (bispectrum_type==lens_bispectrum_ix) then
                                                Bispectrum=>SqueezedLensingKernel(field2,field3)
                                            else
                                                Bispectrum=>bispectra(field1,field2,field3,bispectrum_type)
                                            end if
                                            do f3=1,nfields

                                                do bispectrum_type2=bispectrum_type,nbispectra
                                                    if (bispectrum_type2==lens_bispectrum_ix) then
                                                        Bispectrum2=>SqueezedLensingKernel(f2,f3)
                                                    else
                                                        Bispectrum2=>Bispectra(f1,f2,f3,bispectrum_type2)
                                                    end if


                                                    min_l=minl2
                                                    bi_ix=bix
                                                    if (min_l==l2) then
                                                        !Symmetry factors
                                                        bi_ix=bi_ix+1
                                                        l3=l2
                                                        if (l2==l1) then
                                                            !l1=l2=l3
                                                            tmp = Bispectrum%b(bi_ix,il1)*tmpf(2)*Bispectrum2%b(bi_ix,il1) &
                                                                *InvC(l3)%C(field3,f3)*a3j_00(l3)/6
                                                        else
                                                            !l3=l2 (l3=l1<>l2 can't happen because l1<=l2<=l3)
                                                            tmp = Bispectrum%b(bi_ix,il1)*tmpf(2)*Bispectrum2%b(bi_ix,il1) &
                                                                * InvC(l3)%C(field3,f3)*a3j_00(l3)/2
                                                        end if
                                                        min_l = min_l+2
                                                    else
                                                        tmp=0
                                                    end if
                                                    tmp2=0
                                                    do l3=min_l,max_l ,2
                                                        bi_ix=bi_ix+1
                                                        tmp2 = tmp2 + Bispectrum%b(bi_ix,il1)*Bispectrum2%b(bi_ix,il1) &
                                                            * InvC(l3)%C(field3,f3)*a3j_00(l3)
                                                    end do
                                                    if (l2==l1) then
                                                        tmp2=tmp2*tmpf(2)/2
                                                    else
                                                        tmp2=tmp2*tmpf(2)
                                                    end if
                                                    fish_l1(bispectrum_type,bispectrum_type2,field1,f1)= &
                                                        fish_l1(bispectrum_type,bispectrum_type2,field1,f1)+(tmp+tmp2)

                                                end do !bispectrum_type2

                                            end do
                                        end do !bispectrum_type
                                    end do !field3
                                end do !f2
                            end do !field2
                        end do !f1
                    end do  !field1

                end do !l2

                ifish_contribs(il1,:,:,:,:)=fish_l1(1:nbispectra,1:nbispectra,:,:) /InternalScale
                deallocate(fish_L1)

            end do
            !$OMP END PARAllEl DO
            if (DebugMsgs) print *,'Time for Fisher:', GetTestTime()-starttime


            allocate(fish_contribs(lmin:CTrans%ls%l(CTrans%ls%l0),nfields,nfields))
            allocate(Fisher(nbispectra,nbispectra))
            allocate(tmpFisher(nbispectra,nbispectra))
            allocate(Fisher_L1(lmin:CTrans%ls%l(CTrans%ls%l0),nbispectra*nfields,nbispectra*nfields))

            do bispectrum_type=1,nbispectra
                do bispectrum_type2=bispectrum_type,nbispectra

                    fish_contribs=0
                    do field1=1,nfields
                        do f1=1,nfields
                            call InterpolateClArr(SampleL,ifish_contribs(1,bispectrum_type,bispectrum_type2,field1,f1), &
                                fish_contribs(lmin,field1,f1),lmaxcuti)  !SampleL%l0)
                        end do
                    end do
                    Fisher(bispectrum_type,bispectrum_type2)=0
                    do i=lmin, CTrans%ls%l(CTrans%ls%l0)
                        do field1=1,nfields
                            do f1=1,nfields
                                Fisher_L1(i,(bispectrum_type-1)*nfields+field1,(bispectrum_type2-1)*nfields+f1) = fish_contribs(i,field1,f1)
                                Fisher_L1(i,(bispectrum_type2-1)*nfields+f1,(bispectrum_type-1)*nfields+field1) = fish_contribs(i,field1,f1)
                                tmp=fish_contribs(i,field1,f1)
                                if (bispectrum_type==lens_bispectrum_ix) then
                                    tmp = tmp * CPhi(1+field1,i)
                                end if
                                if (bispectrum_type2==lens_bispectrum_ix) then
                                    tmp = tmp * CPhi(1+f1,i)
                                end if
                                Fisher(bispectrum_type,bispectrum_type2)=Fisher(bispectrum_type,bispectrum_type2)+ tmp
                            end do
                        end do
                    end do
                    Fisher(bispectrum_type2,bispectrum_type) = Fisher(bispectrum_type,bispectrum_type2)

                    print *,'Zero-signal Fisher ',trim(BispectrumNames(bispectrum_type))//'-'//trim(BispectrumNames(bispectrum_type2)), &
                        ':', Fisher(bispectrum_type2,bispectrum_type)

                    !!!! contribution of lensing to the fnl variance for temperature:
                    !           if (bispectrum_type == fnl_bispectrum_ix .and. bispectrum_type2 == fnl_bispectrum_ix ) then
                    !              fish_contribs_sig=0
                    !              tmpArr=0
                    !              call InterpolateClArr(SampleL,ifish_contribs(1,1,2,1,1),tmpArr(lmin),lmaxcuti)
                    !               do i=lmin, lmax_lensing_corrT
                    !               if (CPhi(1+1,i)/=0) then
                    !                 fish_contribs_sig(i) = &
                    !                  tmpArr(i)**2*(1+ CPhi(1,i)*CForLensing(i)%C(1,1)/(CPhi(1+1,i)*CPhi(1+1,i)))/(2*i+1)
                    !               end if
                    !              end do
                    !              print *,'signal contribution to fnl variance', sum(fish_contribs_sig)/Fisher(bispectrum_type,bispectrum_type2)**2
                    !           end if
                    !!! same with polarization
                    !  if (bispectrum_type == fnl_bispectrum_ix .and. bispectrum_type2 == lens_bispectrum_ix ) then
                    !            fish_contribs_sig=0
                    !            do i=lmin, lmax_lensing_corrT
                    !             corrsize = count(CPhi(2:1+nfields,i)/=0)
                    !             allocate(fish_L_ij(corrsize,corrsize))
                    !             allocate(fish_L_noise(corrsize,corrsize))
                    !             fish_L_ij=0
                    !             ix1=0
                    !             do field1=1,nfields
                    !              ix1=ix1+1
                    !              ix2=0
                    !              do f1=1,nfields
                    !                 ix2=ix2+1
                    !                 fish_L_noise(ix1,ix2)= fish_contribs(i,field1, f1)
                    !                 fish_L_ij(ix1,ix2) = (CPhi(1+f1,i)*CPhi(1+field1,i) + &
                    !                     CPhi(1,i)*CForLensing(i)%C(field1,f1) )/(2*i+1)
                    !              end do
                    !             end do
                    !             fish_L_ij=matmul(matmul((fish_L_noise(:,:)),fish_L_ij),transpose(fish_L_noise(:,:)))
                    !             fish_contribs_sig(i) = sum(fish_L_ij)!Matrix_trace(fish_L_ij)
                    !             deallocate(fish_L_ij)
                    !             deallocate(fish_L_noise)
                    !            end do
                    !            print *,'fnl signal factor', sum(fish_contribs_sig)
                    !       end if


                    if (bispectrum_type == lens_bispectrum_ix .and. bispectrum_type2 == lens_bispectrum_ix ) then
                        print *,'doing signal part of the lensing variance'
                        fish_contribs_sig=0
                        do i=lmin, lmax_lensing_corrT
                            corrsize = count(CPhi(2:1+nfields,i)/=0)
                            allocate(fish_L_ij(corrsize,corrsize))
                            allocate(fish_L_noise(corrsize,corrsize))
                            fish_L_ij=0
                            ix1=0
                            do field1=1,nfields
                                if (CPhi(1+field1,i)/=0) then
                                    ix1=ix1+1
                                    ix2=0
                                    do f1=1,nfields
                                        if (CPhi(1+f1,i)/=0) then
                                            ix2=ix2+1
                                            fish_L_noise(ix1,ix2)= fish_contribs(i,field1, f1)*CPhi(1+f1,i)*CPhi(1+field1,i)
                                            fish_L_ij(ix1,ix2) =  (1 + &
                                                CPhi(1,i)*CForLensing(i)%C(field1,f1)/(CPhi(1+f1,i)*CPhi(1+field1,i)))/(2*i+1)

                                        end if
                                    end do
                                end if
                            end do
                            call Matrix_Inverse(fish_L_noise)
                            fish_L_ij=fish_L_ij+fish_L_noise
                            call Matrix_Inverse(fish_L_ij)
                            fish_contribs_sig(i) = sum(fish_L_ij)
                            deallocate(fish_L_ij)
                            deallocate(fish_L_noise)
                        end do
                        print *,'Lensing Fisher including lensing variance ', sum(fish_contribs_sig)

                    end if

                end do
            end do
            deallocate(fish_contribs)

            print *, 'Results assuming zero fiducial bispectra'
            do bispectrum_type=1,nbispectra
                print *,trim(IntToStr(bispectrum_type))//'-'//trim(BispectrumNames(bispectrum_type)), &
                    ': 1/sqrt(F_ii) = ',1/sqrt(Fisher(bispectrum_type,bispectrum_type))
            end do

            do bispectrum_type=1,nbispectra
                do bispectrum_type2=bispectrum_type+1,nbispectra
                    bias = Fisher(bispectrum_type2,bispectrum_type)/Fisher(bispectrum_type,bispectrum_type)
                    print *,'Bias of ',trim(BispectrumNames(bispectrum_type2)),' on ', &
                        trim(BispectrumNames(bispectrum_type)),':', bias
                end do
                !           if (bispectrum_type==1) write(1,concat('(1I10,',nbispectra,'E15.5)')) lmax, Fisher(1,1), &
                !                                  Fisher(2:nbispectra,1)/Fisher(1,1)
            end do

            tmpFisher=Fisher
            do bispectrum_type=1,nbispectra
                tmp = sqrt(tmpFisher(bispectrum_type,bispectrum_type))
                tmpFisher(bispectrum_type,:)=tmpFisher(bispectrum_type,:)/tmp
                tmpFisher(:,bispectrum_type)=tmpFisher(:,bispectrum_type)/tmp
            end do
            if (nbispectra>1) then
                print *,'Zero-signal Bispectrum correlation matrix:'
                do bispectrum_type=1,nbispectra
                    print *,tmpFisher(:,bispectrum_type)
                end do

                if (nbispectra > 1) then
                    allocate(OptimalFisher(nbispectra,nbispectra))
                    allocate(tmpBigFisher(nbispectra*nfields,nbispectra*nfields))
                    !              allocate(diag(nbispectra*nfields))
                    OptimalFisher=0
                    do i=lmin, lmax
                        tmpBigFisher=Fisher_L1(i,:,:)
                        if (lens_bispectrum_ix/=0 .and. i<=lmax_lensing_corrT) then

                            call Matrix_Inverse(tmpBigFisher)
                            do field1 = 1,nfields
                                do field2 =1, nfields
                                    tmpBigFisher((lens_bispectrum_ix-1)*nfields+field1, &
                                        (lens_bispectrum_ix-1)*nfields+field2) = &
                                        tmpBigFisher((lens_bispectrum_ix-1)*nfields+field1, &
                                        (lens_bispectrum_ix-1)*nfields+field2) + &
                                        (CPhi(1+field1,i)*CPhi(1+field2,i) + CPhi(1,i)*CForLensing(i)%C(field1,field2))/(2*i+1)
                                end do
                            end do
                            call Matrix_Inverse(tmpBigFisher)

                        end if

                        do field1=1,nfields
                            do field2=1,nfields
                                do bispectrum_type=1,nbispectra
                                    do bispectrum_type2=1,nbispectra
                                        if (bispectrum_type==lens_bispectrum_ix) then
                                            tmp=CPhi(1+field1,i)
                                        else
                                            tmp=1
                                        end if
                                        if (bispectrum_type2==lens_bispectrum_ix) then
                                            tmp=tmp*CPhi(1+field2,i)
                                        end if
                                        OptimalFisher(bispectrum_type,bispectrum_type2)=OptimalFisher(bispectrum_type,bispectrum_type2)+ &
                                            tmp*tmpBigFisher((bispectrum_type-1)*nfields+field1,(bispectrum_type2-1)*nfields+field2)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    deallocate(tmpBigFisher)
                end if

                do bispectrum_type=1,nbispectra
                    print *,'Optimal Inc. lensing:', trim(IntToStr(bispectrum_type))//'-'//trim(BispectrumNames(bispectrum_type)), &
                        ': 1/sqrt(F_ii) = ',1/sqrt(OptimalFisher(bispectrum_type,bispectrum_type))
                end do
                tmpFisher=OptimalFisher
                call Matrix_Inverse(tmpFIsher)
                do bispectrum_type=1,nbispectra
                    print *,'Optimal Inc. lensing:', trim(IntToStr(bispectrum_type))//'-'//trim(BispectrumNames(bispectrum_type)), &
                        ': Cov_ii = ', sqrt(tmpFIsher(bispectrum_type,bispectrum_type))
                end do
                tmpFisher=OptimalFisher
                do bispectrum_type=1,nbispectra
                    tmp = sqrt(tmpFisher(bispectrum_type,bispectrum_type))
                    tmpFisher(bispectrum_type,:)=tmpFisher(bispectrum_type,:)/tmp
                    tmpFisher(:,bispectrum_type)=tmpFisher(:,bispectrum_type)/tmp
                end do
                print *,'Optimal Bispectrum correlation matrix with lensing variance:'
                do bispectrum_type=1,nbispectra
                    print *,tmpFisher(:,bispectrum_type)
                end do
            end if

            deallocate(Fisher, tmpFisher)

        end do

        deallocate(ifish_contribs)
        deallocate(InvC)

    end if !DoFIsher
#else
    if (BispectrumParams%DoFisher) call MpiStop('compile with FISHER defined')
#endif                    

    !Tidy up a bit
    do field1=1,nfields
        do field2=1,nfields
            do field3=1,nfields
                do bispectrum_type=1,nbispectra
                    deallocate(Bispectra(field1,field2,field3,bispectrum_type)%b)
                end do
            end do
        end do
    end do
    deallocate(Bispectra)
    if (allocated(CForLensing)) deallocate(CForLensing)

    if (allocated(OddBispectra)) then
        oddix=0
        do field1=1,3
            do field2=1,3
                do field3=1,3
                    if (parities(field1)+parities(field2)+parities(field3)/=1) cycle
                    oddix=oddix+1
                    deallocate(OddBispectra(oddix)%b)
                end do
            end do
        end do
        deallocate(OddBispectra)
    end if

    end subroutine GetBispectrum


    !not needed for local NG
    subroutine NonGauss_deriv_l_r(CTrans, indP,resP, r, dJl, dddJl)
    !As above, but integral against derivative of bessel function to get derivative of function
    Type(ClTransferData) :: CTrans
    real(dl), intent(in) :: dJl(BessRanges%npoints,CTrans%ls%l0), dddJl(BessRanges%npoints,CTrans%ls%l0)
    integer, intent(in) :: indP(:)
    real(dl) resP(CTrans%ls%l0,size(indP))
    real(dl), intent(in) :: r
    integer q_ix, j, bes_ix, i
    integer nP, ellmax
    real(dl) xf , dJ_l, fac, a2,  k, dlnk, term, P

    nP =size(indP)
    resP = 0
    do q_ix = 1, CTrans%q%npoints

        k = CTrans%q%points(q_ix)
        xf = k*r  !kr
        bes_ix=Ranges_indexOf(BessRanges,xf)
        fac=BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
        a2=(BessRanges%points(bes_ix+1)-xf)/fac
        fac=fac**2*a2/6
        dlnk = CTrans%q%dpoints(q_ix) /k
        P = ScalarPower(k, 1)  !!only first index for now
        ellmax = max(xf/(1-xlimfrac), xf + xlimmin) * AccuracyBoost

        do j=1,CTrans%ls%l0
            if (CTrans%ls%l(j) <= ellmax) then
                dJ_l=a2*djl(bes_ix,j)+(1-a2)*(djl(bes_ix+1,j) - ((a2+1) &
                    *dddjl(bes_ix,j)+(2-a2)*dddjl(bes_ix+1,j))* fac) !cubic spline

                term = CTrans%Delta_p_l_k(1,j,q_ix)*dJ_l*dlnk *k
                do i=1,nP
                    resP(j,i) = resP(j,i) + term * k**indP(i) * P
                end do

            end if
        end do
    end do
    resP = resP * fourpi

    end subroutine NonGauss_deriv_l_r

    subroutine Bispectrum_SetDefParams(B)
    Type(TBispectrumParams) :: B

    B%nfields=2
    B%Slice_Base_L=0
    B%ndelta=0
    B%DoFisher = .true.
    B%FisherNoise = 0 !2d-4 !Planckish
    B%FisherNoisePol = 4*B%FisherNoise
    B%FisherNoiseFwhmArcmin = 7
    B%FullOutputFile='' !quite large
    B%SparseFullOutput = .false.
    B%do_lensing_bispectrum = .true.
    B%do_primordial_bispectrum = .true.
    B%do_parity_odd = .false.
    B%export_alpha_beta = .false.

    end subroutine Bispectrum_SetDefParams


    subroutine Bispectrum_ReadParams(B, Ini, outroot)
    use IniFile
    Type(TBispectrumParams) :: B
    character(LEN=*), intent(in) :: outroot
    Type(TIniFile) :: Ini
    integer i

    call Bispectrum_SetDefParams(B)

    B%do_lensing_bispectrum = Ini_Read_Logical_File(Ini,'do_lensing_bispectrum',.false.)
    B%do_primordial_bispectrum = Ini_Read_Logical_File(Ini,'do_primordial_bispectrum',.false.)

    do_bispectrum= B%do_lensing_bispectrum .or. B%do_primordial_bispectrum

    if (do_bispectrum) then

        output_root = outroot

        B%nfields = Ini_Read_Int_File(Ini,'bispectrum_nfields',B%nfields)
        if (B%nfields /= 2 .and. B%nfields/=1) call MpiStop('Bispectrum: nfields=1 for T only or 2 for polarization')
        B%do_parity_odd = Ini_Read_Logical_File(Ini,'do_parity_odd',.false.)
        if (B%do_parity_odd .and. (.not.  B%do_lensing_bispectrum .or. B%nfields==1)) then
            B%do_parity_odd = .false.
            write(*,*) 'Ignoring do_parity_odd since do_lensing_bispectrum=F or no polarization'
        end if
        B%Slice_Base_L = Ini_Read_Int_File(Ini,'bispectrum_slice_base_L',B%Slice_Base_L)
        if (B%Slice_Base_L>0) then
            B%ndelta = Ini_Read_Int_File(Ini,'bispectrum_ndelta',B%ndelta)
            if (B%ndelta > max_bispectrum_deltas) call MpiStop('Bispectrum : increase max_bispectrum_deltas')
            do i=1, B%ndelta
                B%deltas(i) = Ini_Read_Int_Array_File(Ini,'bispectrum_delta', i)
            end do
            if (.not. B%do_parity_odd .and. B%Slice_Base_L>0 .and. &
                any(mod(B%Slice_Base_L + B%deltas(1:B%ndelta),2) /= 0)) &
                call MpiStop('Slice is zero for even parity with L1+L2+L3 odd, i.e. Base+DeltaL3 odd')

        end if
        B%DoFisher = Ini_Read_Logical_File(Ini,'bispectrum_do_fisher',B%DoFisher)
        if (B%DoFisher) then
            B%FisherNoise = Ini_Read_Double_File(Ini,'bispectrum_fisher_noise',B%FisherNoise)
            B%FisherNoisePol = Ini_Read_Double_File(Ini,'bispectrum_fisher_noise_pol',B%FisherNoisePol )
            B%FisherNoiseFwhmArcmin = Ini_Read_Double_File(Ini,'bispectrum_fisher_fwhm_arcmin',B%FisherNoiseFwhmArcmin)
        end if
        B%FullOutputFile = Ini_Read_String_File(Ini,'bispectrum_full_output_file')
        if (B%FullOutputFile /='') then
            B%SparseFullOutput = Ini_Read_Logical_file(Ini,'bispectrum_full_output_sparse',B%SparseFullOutput)
        end if

        if (B%do_primordial_bispectrum) then
            B%export_alpha_beta = Ini_Read_Logical_file(Ini,'bispectrum_export_alpha_beta', B%export_alpha_beta)
        end if

    end if

    end subroutine Bispectrum_ReadParams

    end module Bispectrum