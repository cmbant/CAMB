    ! ####### constrained and unconstrained minimization ######
    ! F90 translated version of plato.asu.edu/ftp/other_software/bobyqa.zip and newuoa
    ! Powell 2009 Method http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf

    ! AL Sept 2012/2020: translated to F90, changed to input "funkk" class function argument
    ! to minimize rather than calfun subroutine as original

    ! Typical Usage (typically Npt= 2*dimension+1), to get best fit in params(:) starting from guess

    ! Type(TBOBYQA) :: Minimize !bounded minimization

    ! if (Minimize%BOBYQA(this_class, class_function, dimension, Npt, params,param_min, &
    !           param_max, step_radius, tolerance, DebugLevel, max_eval)) ...
    !

    ! Type(TNEWUOA) :: Minimize !unbounded minimization

    ! if (Minimize%NEWUOA(this_class, class_function, dimension, Npt, param_guess, step_radius, &
    !           tolerance, DebugLevel, max_eval)) ...
    !
    ! The class_function is of the form
    !
    ! real(dp) function func(this, x)
    ! class(TMyType) :: this
    ! real(dp), intent(in) :: x(:)
    ! func = ...f(x)...
    ! end function func

    module Powell
    use precision, only: dl
    implicit none
    private

    integer, parameter :: dp = dl
    integer, parameter :: Powell_CO_prec = dp

    ! the type of function minimized, a class procedure taking in an array of values
    interface
    function obj_vec_function(obj, x)
    import :: dp
    class(*), intent(inout) :: obj
    real(dp), intent(in) :: x(:)
    real(dp) :: obj_vec_function
    end function obj_vec_function
    end interface

    type TMinimizer
        real(dp) :: Last_bestfit = huge(1._dp)
        procedure(obj_vec_function), pointer, nopass :: funkk => null()
        class(*), pointer :: obj => null()
    end type

    type, extends(TMinimizer) :: TBOBYQA
        real(dp) :: FVAL_Converge_difference = 0._dp
    contains
    procedure :: BOBYQA
    end type

    type, extends(TMinimizer) :: TNEWUOA
    contains
    procedure :: NEWUOA
    procedure, private :: NEWUOB
    end type

    public :: TBOBYQA, TNEWUOA, Powell_CO_prec

    contains

    logical function BOBYQA(this, obj, funk, n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun)
    ! Main function for bounded minimization
    class(TBOBYQA) :: this
    class(*), target :: obj
    procedure(obj_vec_function) :: funk ! funk(obj, x), unwrapping obj with select type
    integer, intent(in) :: n, npt, maxfun, iprint
    real(dp), intent(in) :: rhobeg, rhoend
    real(dp), intent(inout) :: x(*)
    real(dp), intent(in) :: xl(*), xu(*)
    integer :: ibmat, id, ifv, igo, ihq, ipq, isl, isu
    integer :: ivl, iw, ixa, ixb, ixn, ixo, ixp, izmat, j, jsl, jsu, ndim, np
    real(dp) :: temp
    real(dp), allocatable :: w(:)
    !   This subroutine seeks the least value of a function of many variables,
    !   by applying a trust region method that forms quadratic models by
    !   interpolation. There is usually some freedom in the interpolation
    !   conditions, which is taken up by minimizing the Frobenius norm of
    !   the change to the second derivative of the model, beginning with the
    !   zero matrix. The values of the variables are constrained by upper and
    !   lower bounds. The arguments of the subroutine are as follows.
    !
    !   N must be set to the number of variables and must be at least two.
    !   NPT is the number of interpolation conditions. Its value must be in
    !     the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
    !     recommended.
    !   Initial values of the variables must be set in X(1),X(2),...,X(N). They
    !     will be changed to the values that give the least calculated F.
    !   For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
    !     bounds, respectively, on X(I). The construction of quadratic models
    !     requires XL(I) to be strictly less than XU(I) for each I. Further,
    !     the contribution to a model from changes to the I-th variable is
    !     damaged severely by rounding errors if XU(I)-XL(I) is too small.
    !   RHOBEG and RHOEND must be set to the initial and final values of a trust
    !     region radius, so both must be positive with RHOEND no greater than
    !     RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
    !     expected change to a variable, while RHOEND should indicate the
    !     accuracy that is required in the final values of the variables. An
    !     error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
    !     is less than 2*RHOBEG.
    !   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
    !     amount of printing. Specifically, there is no output if IPRINT=0 and
    !     there is output only at the return if IPRINT=1. Otherwise, each new
    !     value of RHO is printed, with the best vector of variables so far and
    !     the corresponding value of the objective function. Further, each new
    !     value of F with its variables are output if IPRINT=3.
    !   MAXFUN must be set to an upper bound on the number of calls of CALFUN.
    !   The array W will be used for working space. Its length must be at least
    !     (NPT+5)*(NPT+N)+3*N*(N+5)/2.
    !
    !   funkk(N,X) has to be provided by the user. It must set
    !   F to the value of the objective function for the current values of the
    !   variables X(1),X(2),...,X(N), which are generated automatically in a
    !   way that satisfies the bounds given in XL and XU.
    !
    !   Return if the value of NPT is unacceptable.

    this%funkk => funk
    this%obj => obj

    BOBYQA = .false.
    this%Last_bestfit = huge(1._dp)

    np = n + 1
    if (npt < n + 2 .or. npt > ((n + 2)*np)/2) then
        write(*, *) "Return from BOBYQA because NPT is not in the required interval"
        return
    end if

    allocate(w((npt + 5)*(npt + n) + 3*n*(n + 5)/2))

    !  Partition the working space array, so that different parts of it can
    !  be treated separately during the calculation of BOBYQB. The partition
    !  requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
    !  space that is taken by the last array in the argument list of BOBYQB.
    !
    ndim = npt + n
    ixb = 1
    ixp = ixb + n
    ifv = ixp + n*npt
    ixo = ifv + npt
    igo = ixo + n
    ihq = igo + n
    ipq = ihq + (n*np)/2
    ibmat = ipq + npt
    izmat = ibmat + ndim*n
    isl = izmat + npt*(npt - np)
    isu = isl + n
    ixn = isu + n
    ixa = ixn + n
    id = ixa + n
    ivl = id + n
    iw = ivl + ndim
    !
    !  Return if there is insufficient space between the bounds. Modify the
    !  initial X if necessary in order to avoid conflicts between the bounds
    !  and the construction of the first quadratic model. The lower and upper
    !  bounds on moves from the updated X are set now, in the ISL and ISU
    !  partitions of W, in order to provide useful and exact information about
    !  components of X that become within distance RHOBEG from their bounds.
    !
    do j = 1, n
        temp = xu(j) - xl(j)
        if (temp < rhobeg + rhobeg) then
            write(*, *) "Return from BOBYQA because one of the differences XU(I)-XL(I) is less than 2*RHOBEG."
            return
        end if
        jsl = isl + j - 1
        jsu = jsl + n
        w(jsl) = xl(j) - x(j)
        w(jsu) = xu(j) - x(j)
        if (w(jsl) >= -rhobeg) then
            if (w(jsl) >= 0._dp) then
                x(j) = xl(j)
                w(jsl) = 0._dp
                w(jsu) = temp
            else
                x(j) = xl(j) + rhobeg
                w(jsl) = -rhobeg
                w(jsu) = max(xu(j) - x(j), rhobeg)
            end if
        else if (w(jsu) <= rhobeg) then
            if (w(jsu) <= 0._dp) then
                x(j) = xu(j)
                w(jsl) = -temp
                w(jsu) = 0._dp
            else
                x(j) = xu(j) - rhobeg
                w(jsl) = min(xl(j) - x(j), -rhobeg)
                w(jsu) = rhobeg
            end if
        end if
    end do
    !
    !     Make the call of BOBYQB.
    !
    BOBYQA = BOBYQB(this, n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w(ixb), &
        w(ixp), w(ifv), w(ixo), w(igo), w(ihq), w(ipq), w(ibmat), w(izmat), &
        ndim, w(isl), w(isu), w(ixn), w(ixa), w(id), w(ivl), w(iw))

    end function BOBYQA

    function BOBYQB (this, N, NPT, X, XL, XU, Rhobeg, Rhoend, Iprint, &
        Maxfun, XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, &
        SL, SU, XNEW, XALT, D, VLAG, W)
    class(TBOBYQA) this
    integer NDIM, NPT
    real(dp) X(*), XL(*), XU(*), XBASE(*), XPT(NPT, *), FVAL(*), &
        XOPT(*), GOPT(*), HQ(*), PQ(*), BMAT(NDIM, *), ZMAT(NPT, *), &
        SL(*), SU(*), XNEW(*), XALT(*), D(*), VLAG(*), W(*)
    logical BOBYQB
    real(dp) adelt, alpha, bdtest, bdtol, beta, biglsq, &
        bsum, cauchy, crvmin, curv, delsq, delta, den, &
        denom, densav, diff, diffa, diffb
    real(dp) diffc, dist, distsq, dnorm, dsq, dx, errbig, f, &
        fopt, fracsq, frhosq, fsave, gisq, &
        gqsq, half, hdiag, one
    real(dp) pqold, ratio, rho, Rhobeg, Rhoend, scaden, &
        sum, suma, sumb, sumpq, sumw, sumz, temp, ten, &
        tenth, two
    real(dp) vquad, xoptsq, zero
    integer i, ih, ip, Iprint, itest, j, jj, jp, k, kbase, &
        knew, kopt, ksav, Maxfun, N, nf, nfsav, nh, np
    integer nptm, nresc, ntrits

    !  The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
    !    are identical to the corresponding arguments in SUBROUTINE BOBYQA.
    !  XBASE holds a shift of origin that should reduce the contributions
    !    from rounding errors to values of the model and Lagrange functions.
    !  XPT is a two-dimensional array that holds the coordinates of the
    !    interpolation points relative to XBASE.
    !  FVAL holds the values of F at the interpolation points.
    !  XOPT is set to the displacement from XBASE of the trust region centre.
    !  GOPT holds the gradient of the quadratic model at XBASE+XOPT.
    !  HQ holds the explicit second derivatives of the quadratic model.
    !  PQ contains the parameters of the implicit second derivatives of the
    !    quadratic model.
    !  BMAT holds the last N columns of H.
    !  ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
    !    this factorization being ZMAT times ZMAT^T, which provides both the
    !    correct rank and positive semi-definiteness.
    !  NDIM is the first dimension of BMAT and has the value NPT+N.
    !  SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
    !    All the components of every XOPT are going to satisfy the bounds
    !    SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
    !    XOPT is on a constraint boundary.
    !  XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
    !    vector of variables for the next call of CALFUN. XNEW also satisfies
    !    the SL and SU constraints in the way that has just been mentioned.
    !  XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
    !    in order to increase the denominator in the updating of UPDATE.
    !  D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
    !  VLAG contains the values of the Lagrange functions at a new point X.
    !    They are part of a product that requires VLAG to be of length NDIM.
    !  W is a one-dimensional array that is used for working space. Its length
    !    must be at least 3*NDIM = 3*(NPT+N).
    !
    !     Set some constants.
    !
    half = 0.5d0
    one = 1.0d0
    ten = 10.0d0
    tenth = 0.1d0
    two = 2.0d0
    zero = 0.0d0
    np = N + 1
    nptm = NPT - np
    nh = (N*np)/2
    BOBYQB = .true.
    !
    !  The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    !  BMAT and ZMAT for the first iteration, with the corresponding values of
    !  of NF and KOPT, which are the number of calls of CALFUN so far and the
    !  index of the interpolation point at the trust region centre. Then the
    !  initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
    !  less than NPT. GOPT will be updated if KOPT is different from KBASE.
    !
    call PRELIM (this, N, NPT, X, XL, XU, Rhobeg, Iprint, Maxfun, XBASE, XPT, &
        FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, nf, kopt)
    xoptsq = zero
    do i = 1, N
        XOPT(i) = XPT(kopt, i)
        xoptsq = xoptsq + XOPT(i)**2
    end do
    fsave = FVAL(1)
    if (nf < NPT) then
        if (Iprint > 0) write(*, *) "Return from BOBYQA because funkk has been called MAXFUN times."
        goto 720
    end if
    kbase = 1
    !
    !     Complete the settings that are required for the iterative procedure
    !
    rho = Rhobeg
    delta = rho
    nresc = nf
    ntrits = 0
    diffa = zero
    diffb = zero
    itest = 0
    nfsav = nf
    !
    !     Update GOPT if necessary before the first iteration and after each
    !     call of RESCUE that makes a call of CALFUN.
    !
20  if (kopt /= kbase) then
        ih = 0
        do j = 1, N
            do i = 1, j
                ih = ih + 1
                if (i < j) GOPT(j) = GOPT(j) + HQ(ih)*XOPT(i)
                GOPT(i) = GOPT(i) + HQ(ih)*XOPT(j)
            end do
        end do
        if (nf > NPT) then
            do k = 1, NPT
                temp = zero
                do j = 1, N
                    temp = temp + XPT(k, j)*XOPT(j)
                end do
                temp = PQ(k)*temp
                do i = 1, N
                    GOPT(i) = GOPT(i) + temp*XPT(k, i)
                end do
            end do
        end if
    end if
    !
    !  Generate the next point in the trust region that provides a small value
    !  of the quadratic model subject to the constraints on the variables.
    !  The integer NTRITS is set to the number "trust region" iterations that
    !  have occurred since the last "alternative" iteration. If the length
    !  of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
    !  label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
    !
60  call TRSBOX (N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL, SU, delta, XNEW, D, &
        W, W(np), W(np + N), W(np + 2*N), W(np + 3*N), dsq, crvmin)
    dnorm = min(delta, sqrt(dsq))
    if (dnorm < half*rho) then
        ntrits = -1
        distsq = (ten*rho)**2
        if (nf <= nfsav + 2) goto 650
        !
        !  The following choice between labels 650 and 680 depends on whether or
        !  not our work with the current RHO seems to be complete. Either RHO is
        !  decreased or termination occurs if the errors in the quadratic model at
        !  the last three interpolation points compare favourably with predictions
        !  of likely improvements to the model within distance HALF*RHO of XOPT.
        !
        errbig = max(diffa, diffb, diffc)
        frhosq = 0.125d0*rho*rho
        if (crvmin > zero .and. errbig > frhosq*crvmin) &
            goto 650
        bdtol = errbig/rho
        do j = 1, N
            bdtest = bdtol
            if (XNEW(j) == SL(j)) bdtest = W(j)
            if (XNEW(j) == SU(j)) bdtest = -W(j)
            if (bdtest < bdtol) then
                curv = HQ((j + j*j)/2)
                do k = 1, NPT
                    curv = curv + PQ(k)*XPT(k, j)**2
                end do
                bdtest = bdtest + half*curv*rho
                if (bdtest < bdtol) goto 650
            end if
        end do
        goto 680
    end if
    ntrits = ntrits + 1
    !
    !  Severe cancellation is likely to occur if XOPT is too far from XBASE.
    !  If the following test holds, then XBASE is shifted so that XOPT becomes
    !  zero. The appropriate changes are made to BMAT and to the second
    !  derivatives of the current model, beginning with the changes to BMAT
    !  that do not depend on ZMAT. VLAG is used temporarily for working space.
    !
90  if (dsq <= 1.0d-3*xoptsq) then
        fracsq = 0.25d0*xoptsq
        sumpq = zero
        do k = 1, NPT
            sumpq = sumpq + PQ(k)
            sum = -half*xoptsq
            do i = 1, N
                sum = sum + XPT(k, i)*XOPT(i)
            end do
            W(NPT + k) = sum
            temp = fracsq - half*sum
            do i = 1, N
                W(i) = BMAT(k, i)
                VLAG(i) = sum*XPT(k, i) + temp*XOPT(i)
                ip = NPT + i
                do j = 1, i
                    BMAT(ip, j) = BMAT(ip, j) + W(i)*VLAG(j) + VLAG(i)*W(j)
                end do
            end do
        end do
        !
        !     Then the revisions of BMAT that depend on ZMAT are calculated.
        !
        do jj = 1, nptm
            sumz = zero
            sumw = zero
            do k = 1, NPT
                sumz = sumz + ZMAT(k, jj)
                VLAG(k) = W(NPT + k)*ZMAT(k, jj)
                sumw = sumw + VLAG(k)
            end do
            do j = 1, N
                sum = (fracsq*sumz - half*sumw)*XOPT(j)
                do k = 1, NPT
                    sum = sum + VLAG(k)*XPT(k, j)
                end do
                W(j) = sum
                do k = 1, NPT
                    BMAT(k, j) = BMAT(k, j) + sum*ZMAT(k, jj)
                end do
            end do
            do i = 1, N
                ip = i + NPT
                temp = W(i)
                do j = 1, i
                    BMAT(ip, j) = BMAT(ip, j) + temp*W(j)
                end do
            end do
        end do
        !
        !  The following instructions complete the shift, including the changes
        !  to the second derivative parameters of the quadratic model.
        !
        ih = 0
        do j = 1, N
            W(j) = -half*sumpq*XOPT(j)
            do k = 1, NPT
                W(j) = W(j) + PQ(k)*XPT(k, j)
                XPT(k, j) = XPT(k, j) - XOPT(j)
            end do
            do i = 1, j
                ih = ih + 1
                HQ(ih) = HQ(ih) + W(i)*XOPT(j) + XOPT(i)*W(j)
                BMAT(NPT + i, j) = BMAT(NPT + j, i)
            end do
        end do
        do i = 1, N
            XBASE(i) = XBASE(i) + XOPT(i)
            XNEW(i) = XNEW(i) - XOPT(i)
            SL(i) = SL(i) - XOPT(i)
            SU(i) = SU(i) - XOPT(i)
            XOPT(i) = zero
        end do
        xoptsq = zero
    end if
    if (ntrits == 0) goto 210
    goto 230
    !
    !  XBASE is also moved to XOPT by a call of RESCUE. This calculation is
    !  more expensive than the previous shift, because new matrices BMAT and
    !  ZMAT are generated from scratch, which may include the replacement of
    !  interpolation points whose positions seem to be causing near linear
    !  dependence in the interpolation conditions. Therefore RESCUE is called
    !  only if rounding errors have reduced by at least a factor of two the
    !  denominator of the formula for updating the H matrix. It provides a
    !  useful safeguard, but is not invoked in most applications of BOBYQA.
    !
190 nfsav = nf
    kbase = kopt
    call RESCUE (this, N, NPT, XL, XU, Iprint, Maxfun, XBASE, XPT, FVAL, &
        XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, nf, delta, kopt, &
        VLAG, W, W(N + np), W(NDIM + np))
    !
    !  XOPT is updated now in case the branch below to label 720 is taken.
    !  Any updating of GOPT occurs after the branch below to label 20, which
    !  leads to a trust region iteration as does the branch to label 60.
    !
    xoptsq = zero
    if (kopt /= kbase) then
        do i = 1, N
            XOPT(i) = XPT(kopt, i)
            xoptsq = xoptsq + XOPT(i)**2
        end do
    end if
    if (nf < 0) then
        nf = Maxfun
        if (Iprint > 0) write(*, *) "Return from BOBYQA because funkk has been called MAXFUN times."
        goto 720
    end if
    nresc = nf
    if (nfsav < nf) then
        nfsav = nf
        goto 20
    end if
    if (ntrits > 0) goto 60
    !
    !  Pick two alternative vectors of variables, relative to XBASE, that
    !  are suitable as new positions of the KNEW-th interpolation point.
    !  Firstly, XNEW is set to the point on a line through XOPT and another
    !  interpolation point that minimizes the predicted value of the next
    !  denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
    !  and SU bounds. Secondly, XALT is set to the best feasible point on
    !  a constrained version of the Cauchy step of the KNEW-th Lagrange
    !  function, the corresponding value of the square of this function
    !  being returned in CAUCHY. The choice between these alternatives is
    !  going to be made when the denominator is calculated.
    !
210 call ALTMOV (N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL, SU, kopt, &
        knew, adelt, XNEW, XALT, alpha, cauchy, W, W(np), W(NDIM + 1))
    do i = 1, N
        D(i) = XNEW(i) - XOPT(i)
    end do
    !
    !  Calculate VLAG and BETA for the current choice of D. The scalar
    !  product of D with XPT(K,.) is going to be held in W(NPT+K) for
    !  use when VQUAD is calculated.
    !
230 do k = 1, NPT
        suma = zero
        sumb = zero
        sum = zero
        do j = 1, N
            suma = suma + XPT(k, j)*D(j)
            sumb = sumb + XPT(k, j)*XOPT(j)
            sum = sum + BMAT(k, j)*D(j)
        end do
        W(k) = suma*(half*suma + sumb)
        VLAG(k) = sum
        W(NPT + k) = suma
    end do
    beta = zero
    do jj = 1, nptm
        sum = zero
        do k = 1, NPT
            sum = sum + ZMAT(k, jj)*W(k)
        end do
        beta = beta - sum*sum
        do k = 1, NPT
            VLAG(k) = VLAG(k) + sum*ZMAT(k, jj)
        end do
    end do
    dsq = zero
    bsum = zero
    dx = zero
    do j = 1, N
        dsq = dsq + D(j)**2
        sum = zero
        do k = 1, NPT
            sum = sum + W(k)*BMAT(k, j)
        end do
        bsum = bsum + sum*D(j)
        jp = NPT + j
        do i = 1, N
            sum = sum + BMAT(jp, i)*D(i)
        end do
        VLAG(jp) = sum
        bsum = bsum + sum*D(j)
        dx = dx + D(j)*XOPT(j)
    end do
    beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta - bsum
    VLAG(kopt) = VLAG(kopt) + one
    !
    !  If NTRITS is zero, the denominator may be increased by replacing
    !  the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
    !  rounding errors have damaged the chosen denominator.
    !
    if (ntrits == 0) then
        denom = VLAG(knew)**2 + alpha*beta
        if (denom < cauchy .and. cauchy > zero) then
            do i = 1, N
                XNEW(i) = XALT(i)
                D(i) = XNEW(i) - XOPT(i)
            end do
            cauchy = zero
            goto 230
        end if
        if (denom <= half*VLAG(knew)**2) then
            if (nf > nresc) goto 190
            if (Iprint > 0) write(*, *) "Return from BOBYQA because of much cancellation in a denominator."
            goto 720
        end if
        !
        !  Alternatively, if NTRITS is positive, then set KNEW to the index of
        !  the next interpolation point to be deleted to make room for a trust
        !  region step. Again RESCUE may be called if rounding errors have damaged
        !  the chosen denominator, which is the reason for attempting to select
        !  KNEW before calculating the next value of the objective function.
        !
    else
        delsq = delta*delta
        scaden = zero
        biglsq = zero
        knew = 0
        do k = 1, NPT
            if (k == kopt) continue
            hdiag = zero
            do jj = 1, nptm
                hdiag = hdiag + ZMAT(k, jj)**2
            end do
            den = beta*hdiag + VLAG(k)**2
            distsq = zero
            do j = 1, N
                distsq = distsq + (XPT(k, j) - XOPT(j))**2
            end do
            temp = max(one, (distsq/delsq)**2)
            if (temp*den > scaden) then
                scaden = temp*den
                knew = k
                denom = den
            end if
            biglsq = max(biglsq, temp*VLAG(k)**2)
        end do
        if (scaden <= half*biglsq) then
            if (nf > nresc) goto 190
            if (Iprint > 0) write(*, *) "Return from BOBYQA because of much cancellation in a denominator."
            goto 720
        end if
    end if
    !
    !  Put the variables for the next calculation of the objective function
    !    in XNEW, with any adjustments for the bounds.
    !
    !
    !  Calculate the value of the objective function at XBASE+XNEW, unless
    !    the limit on the number of calculations of F has been reached.
    !
360 do i = 1, N
        X(i) = min(max(XL(i), XBASE(i) + XNEW(i)), XU(i))
        if (XNEW(i) == SL(i)) X(i) = XL(i)
        if (XNEW(i) == SU(i)) X(i) = XU(i)
    end do
    if (nf >= Maxfun) then
        if (Iprint > 0) write(*, *) "Return from BOBYQA because funkk has been called MAXFUN times."
        BOBYQB = .false.
        goto 720
    end if
    nf = nf + 1
    f = this%funkk(this%obj, X(1:N))
    if (Iprint == 3) then
        write(*, "(/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))") &
            nf, f, (X(i), i=1, N)
    end if
    if (ntrits == -1) then
        fsave = f
        goto 720
    end if
    !
    !  Use the quadratic model to predict the change in F due to the step D,
    !    and set DIFF to the error of this prediction.
    !
    fopt = FVAL(kopt)
    vquad = zero
    ih = 0
    do j = 1, N
        vquad = vquad + D(j)*GOPT(j)
        do i = 1, j
            ih = ih + 1
            temp = D(i)*D(j)
            if (i == j) temp = half*temp
            vquad = vquad + HQ(ih)*temp
        end do
    end do
    do k = 1, NPT
        vquad = vquad + half*PQ(k)*W(NPT + k)**2
    end do
    diff = f - fopt - vquad
    diffc = diffb
    diffb = diffa
    diffa = abs(diff)
    if (dnorm > rho) nfsav = nf
    !
    !     Pick the next value of DELTA after a trust region step.
    !
    if (ntrits > 0) then
        if (vquad >= zero) then
            if (Iprint > 0) write(*, *) "Return from BOBYQA because a trust region step has failed to reduce Q."
            BOBYQB = .false.
            goto 720
        end if
        ratio = (f - fopt)/vquad
        if (ratio <= tenth) then
            delta = min(half*delta, dnorm)
        else if (ratio <= 0.7d0) then
            delta = max(half*delta, dnorm)
        else
            delta = max(half*delta, dnorm + dnorm)
        end if
        if (delta <= 1.5d0*rho) delta = rho
        !
        !     Recalculate KNEW and DENOM if the new F is less than FOPT.
        !
        if (f < fopt) then
            ksav = knew
            densav = denom
            delsq = delta*delta
            scaden = zero
            biglsq = zero
            knew = 0
            do k = 1, NPT
                hdiag = zero
                do jj = 1, nptm
                    hdiag = hdiag + ZMAT(k, jj)**2
                end do
                den = beta*hdiag + VLAG(k)**2
                distsq = zero
                do j = 1, N
                    distsq = distsq + (XPT(k, j) - XNEW(j))**2
                end do
                temp = max(one, (distsq/delsq)**2)
                if (temp*den > scaden) then
                    scaden = temp*den
                    knew = k
                    denom = den
                end if
                biglsq = max(biglsq, temp*VLAG(k)**2)
            end do
            if (scaden <= half*biglsq) then
                knew = ksav
                denom = densav
            end if
        end if
    end if
    !
    !  Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
    !  moved. Also update the second derivative terms of the model.
    !
    call UPDATE (N, NPT, BMAT, ZMAT, NDIM, VLAG, beta, denom, knew, W)
    ih = 0
    pqold = PQ(knew)
    PQ(knew) = zero
    do i = 1, N
        temp = pqold*XPT(knew, i)
        do j = 1, i
            ih = ih + 1
            HQ(ih) = HQ(ih) + temp*XPT(knew, j)
        end do
    end do
    do jj = 1, nptm
        temp = diff*ZMAT(knew, jj)
        do k = 1, NPT
            PQ(k) = PQ(k) + temp*ZMAT(k, jj)
        end do
    end do
    !
    !  Include the new interpolation point, and make the changes to GOPT at
    !  the old XOPT that are caused by the updating of the quadratic model.
    !
    FVAL(knew) = f
    do i = 1, N
        XPT(knew, i) = XNEW(i)
        W(i) = BMAT(knew, i)
    end do
    do k = 1, NPT
        suma = zero
        do jj = 1, nptm
            suma = suma + ZMAT(knew, jj)*ZMAT(k, jj)
        end do
        sumb = zero
        do j = 1, N
            sumb = sumb + XPT(k, j)*XOPT(j)
        end do
        temp = suma*sumb
        do i = 1, N
            W(i) = W(i) + temp*XPT(k, i)
        end do
    end do
    do i = 1, N
        GOPT(i) = GOPT(i) + diff*W(i)
    end do
    !
    !     Update XOPT, GOPT and KOPT if the new calculated F is less than FO
    !
    if (f < fopt) then
        kopt = knew
        xoptsq = zero
        ih = 0
        do j = 1, N
            XOPT(j) = XNEW(j)
            xoptsq = xoptsq + XOPT(j)**2
            do i = 1, j
                ih = ih + 1
                if (i < j) GOPT(j) = GOPT(j) + HQ(ih)*D(i)
                GOPT(i) = GOPT(i) + HQ(ih)*D(j)
            end do
        end do
        do k = 1, NPT
            temp = zero
            do j = 1, N
                temp = temp + XPT(k, j)*D(j)
            end do
            temp = PQ(k)*temp
            do i = 1, N
                GOPT(i) = GOPT(i) + temp*XPT(k, i)
            end do
        end do
    end if
    !
    !  Calculate the parameters of the least Frobenius norm interpolant to
    !  the current data, the gradient of this interpolant at XOPT being put
    !  into VLAG(NPT+I), I=1,2,...,N.
    !
    if (ntrits > 0) then
        do k = 1, NPT
            VLAG(k) = FVAL(k) - FVAL(kopt)
            W(k) = zero
        end do
        do j = 1, nptm
            sum = zero
            do k = 1, NPT
                sum = sum + ZMAT(k, j)*VLAG(k)
            end do
            do k = 1, NPT
                W(k) = W(k) + sum*ZMAT(k, j)
            end do
        end do
        do k = 1, NPT
            sum = zero
            do j = 1, N
                sum = sum + XPT(k, j)*XOPT(j)
            end do
            W(k + NPT) = W(k)
            W(k) = sum*W(k)
        end do
        gqsq = zero
        gisq = zero
        do i = 1, N
            sum = zero
            do k = 1, NPT
                sum = sum + BMAT(k, i)*VLAG(k) + XPT(k, i)*W(k)
            end do
            if (XOPT(i) == SL(i)) then
                gqsq = gqsq + min(zero, GOPT(i))**2
                gisq = gisq + min(zero, sum)**2
            else if (XOPT(i) == SU(i)) then
                gqsq = gqsq + max(zero, GOPT(i))**2
                gisq = gisq + max(zero, sum)**2
            else
                gqsq = gqsq + GOPT(i)**2
                gisq = gisq + sum*sum
            end if
            VLAG(NPT + i) = sum
        end do
        !
        !  Test whether to replace the new quadratic model by the least Frobenius
        !  norm interpolant, making the replacement if the test is satisfied.
        !
        itest = itest + 1
        if (gqsq < ten*gisq) itest = 0
        if (itest >= 3) then
            do i = 1, max(NPT, nh)
                if (i <= N) GOPT(i) = VLAG(NPT + i)
                if (i <= NPT) PQ(i) = W(NPT + i)
                if (i <= nh) HQ(i) = zero
                itest = 0
            end do
        end if
    end if
    !
    !  If a trust region step has provided a sufficient decrease in F, then
    !  branch for another trust region calculation. The case NTRITS=0 occurs
    !  when the new interpolation point was reached by an alternative step.
    !
    if (ntrits == 0) goto 60
    if (f <= fopt + tenth*vquad) goto 60
    !
    !     Alternatively, find out if the interpolation points are close enough
    !       to the best point so far.
    !
    distsq = max((two*delta)**2, (ten*rho)**2)
650 knew = 0
    do k = 1, NPT
        sum = zero
        do j = 1, N
            sum = sum + (XPT(k, j) - XOPT(j))**2
        end do
        if (sum > distsq) then
            knew = k
            distsq = sum
        end if
    end do
    !
    !  If KNEW is positive, then ALTMOV finds alternative new positions for
    !  the KNEW-th interpolation point within distance ADELT of XOPT. It is
    !  reached via label 90. Otherwise, there is a branch to label 60 for
    !  another trust region iteration, unless the calculations with the
    !  current RHO are complete.
    !
    if (knew > 0) then
        dist = sqrt(distsq)
        if (ntrits == -1) then
            delta = min(tenth*delta, half*dist)
            if (delta <= 1.5d0*rho) delta = rho
        end if
        ntrits = 0
        adelt = max(min(tenth*dist, delta), rho)
        dsq = adelt*adelt
        goto 90
    end if
    if (ntrits == -1) goto 680
    if (ratio > zero) goto 60
    if (max(delta, dnorm) > rho) goto 60
    !
    !  The calculations with the current value of RHO are complete. Pick the
    !    next values of RHO and DELTA.
    !
680 if (rho > Rhoend .and. (this%FVAL_Converge_difference == 0._dp .or. &
        abs(FVAL(kopt) - this%Last_bestfit) > this%FVAL_Converge_difference)) then
        delta = half*rho
        ratio = rho/Rhoend
        if (ratio <= 16.0d0) then
            rho = Rhoend
        else if (ratio <= 250.0d0) then
            rho = sqrt(ratio)*Rhoend
        else
            rho = tenth*rho
        end if
        delta = max(delta, rho)
        this%Last_bestfit = FVAL(kopt)
        if (Iprint >= 2) then
            if (Iprint >= 3) print "(5X)"
            write(*, "(/4X,'New RHO =',1PD11.4,5X,'Number of function values =',I6)") rho, nf
            write(*, "(4X,'Least value of F =',1PD23.15,9X,'The corresponding X is:'/(2X,5D15.6))") &
                FVAL(kopt), (XBASE(i) + XOPT(i), i=1, N)
        end if
        ntrits = 0
        nfsav = nf
        goto 60
    end if
    !
    !     Return from the calculation, after another Newton-Raphson step, if
    !       it is too short to have been tried before.
    !
    if (ntrits == -1) goto 360
720 if (FVAL(kopt) <= fsave) then
        do i = 1, N
            X(i) = min(max(XL(i), XBASE(i) + XOPT(i)), XU(i))
            if (XOPT(i) == SL(i)) X(i) = XL(i)
            if (XOPT(i) == SU(i)) X(i) = XU(i)
        end do
        f = FVAL(kopt)
    end if
    if (Iprint >= 1) then
        if (rho > Rhoend) write(*, *) '    FVAL_Converge_difference reached'
        write(*, "(/4X,'At the return from BOBYQA',5X,'Number of function values =',I6)") nf
        write(*, "(4X,'Least value of F =',1PD23.15,9X,'The corresponding X is:'/(2X,5D15.6))") &
            f, (X(i), i=1, N)
    end if

    end function BOBYQB

    subroutine ALTMOV (N, Npt, XPT, XOPT, BMAT, ZMAT, Ndim, SL, SU, Kopt, &
        Knew, Adelt, XNEW, XALT, Alpha, Cauchy, GLAG, HCOL, W)
    real(dp) Adelt, Alpha, bigstp, Cauchy, const, csave, &
        curv, dderiv, diff, distsq, ggfree, gw, ha, &
        half, one, predsq, presav
    real(dp) scale, slbd, step, stpsav, subd, sumin, &
        temp, tempa, tempb, tempd, vlag, wfixsq, wsqsav, zero
    integer i, ibdsav, iflag, ilbd, isbd, iubd, j, k, Knew, &
        Kopt, ksav, N, Ndim, Npt
    real(dp) XPT(Npt, *), XOPT(*), BMAT(Ndim, *), ZMAT(Npt, *), SL(*), &
        SU(*), XNEW(*), XALT(*), GLAG(*), HCOL(*), W(*)
    !
    !  The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
    !    the same meanings as the corresponding arguments of BOBYQB.
    !  KOPT is the index of the optimal interpolation point.
    !  KNEW is the index of the interpolation point that is going to be moved.
    !  ADELT is the current trust region bound.
    !  XNEW will be set to a suitable new position for the interpolation point
    !    XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
    !    bounds and it should provide a large denominator in the next call of
    !    UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
    !    straight lines through XOPT and another interpolation point.
    !  XALT also provides a large value of the modulus of the KNEW-th Lagrange
    !    function subject to the constraints that have been mentioned, its main
    !    difference from XNEW being that XALT-XOPT is a constrained version of
    !    the Cauchy step within the trust region. An exception is that XALT is
    !    not calculated if all components of GLAG (see below) are zero.
    !  ALPHA will be set to the KNEW-th diagonal element of the H matrix.
    !  CAUCHY will be set to the square of the KNEW-th Lagrange function at
    !    the step XALT-XOPT from XOPT for the vector XALT that is returned,
    !    except that CAUCHY is set to zero if XALT is not calculated.
    !  GLAG is a working space vector of length N for the gradient of the
    !    KNEW-th Lagrange function at XOPT.
    !  HCOL is a working space vector of length NPT for the second derivative
    !    coefficients of the KNEW-th Lagrange function.
    !  W is a working space vector of length 2N that is going to hold the
    !    constrained Cauchy step from XOPT of the Lagrange function, followed
    !    by the downhill version of XALT when the uphill step is calculated.
    !
    !  Set the first NPT components of W to the leading elements of the
    !  KNEW-th column of the H matrix.
    !
    half = 0.5d0
    one = 1.0d0
    zero = 0.0d0
    const = one + sqrt(2.0d0)
    do k = 1, Npt
        HCOL(k) = zero
    end do
    do j = 1, Npt - N - 1
        temp = ZMAT(Knew, j)
        do k = 1, Npt
            HCOL(k) = HCOL(k) + temp*ZMAT(k, j)
        end do
    end do
    Alpha = HCOL(Knew)
    ha = half*Alpha
    !
    !     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
    !
    do i = 1, N
        GLAG(i) = BMAT(Knew, i)
    end do
    do k = 1, Npt
        temp = zero
        do j = 1, N
            temp = temp + XPT(k, j)*XOPT(j)
        end do
        temp = HCOL(k)*temp
        do i = 1, N
            GLAG(i) = GLAG(i) + temp*XPT(k, i)
        end do
    end do
    !
    !  Search for a large denominator along the straight lines through XOPT
    !  and another interpolation point. SLBD and SUBD will be lower and upper
    !  bounds on the step along each of these lines in turn. PREDSQ will be
    !  set to the square of the predicted denominator for each line. PRESAV
    !  will be set to the largest admissible value of PREDSQ that occurs.
    presav = zero
    do k = 1, Npt
        if (k == Kopt) continue
        dderiv = zero
        distsq = zero
        do i = 1, N
            temp = XPT(k, i) - XOPT(i)
            dderiv = dderiv + GLAG(i)*temp
            distsq = distsq + temp*temp
        end do
        subd = Adelt/sqrt(distsq)
        slbd = -subd
        ilbd = 0
        iubd = 0
        sumin = min(one, subd)
        !
        !     Revise SLBD and SUBD if necessary because of the bounds in SL and
        !
        do i = 1, N
            temp = XPT(k, i) - XOPT(i)
            if (temp > zero) then
                if (slbd*temp < SL(i) - XOPT(i)) then
                    slbd = (SL(i) - XOPT(i))/temp
                    ilbd = -i
                end if
                if (subd*temp > SU(i) - XOPT(i)) then
                    subd = max(sumin, (SU(i) - XOPT(i))/temp)
                    iubd = i
                end if
            else if (temp < zero) then
                if (slbd*temp > SU(i) - XOPT(i)) then
                    slbd = (SU(i) - XOPT(i))/temp
                    ilbd = i
                end if
                if (subd*temp < SL(i) - XOPT(i)) then
                    subd = max(sumin, (SL(i) - XOPT(i))/temp)
                    iubd = -i
                end if
            end if
        end do
        !
        !  Seek a large modulus of the KNEW-th Lagrange function when the index
        !  of the other interpolation point on the line through XOPT is KNEW.
        !
        if (k == Knew) then
            diff = dderiv - one
            step = slbd
            vlag = slbd*(dderiv - slbd*diff)
            isbd = ilbd
            temp = subd*(dderiv - subd*diff)
            if (abs(temp) > abs(vlag)) then
                step = subd
                vlag = temp
                isbd = iubd
            end if
            tempd = half*dderiv
            tempa = tempd - diff*slbd
            tempb = tempd - diff*subd
            if (tempa*tempb < zero) then
                temp = tempd*tempd/diff
                if (abs(temp) > abs(vlag)) then
                    step = tempd/diff
                    vlag = temp
                    isbd = 0
                end if
            end if
            !
            !     Search along each of the other lines through XOPT and another point
            !
        else
            step = slbd
            vlag = slbd*(one - slbd)
            isbd = ilbd
            temp = subd*(one - subd)
            if (abs(temp) > abs(vlag)) then
                step = subd
                vlag = temp
                isbd = iubd
            end if
            if (subd > half) then
                if (abs(vlag) < 0.25d0) then
                    step = half
                    vlag = 0.25d0
                    isbd = 0
                end if
            end if
            vlag = vlag*dderiv
        end if
        !
        !     Calculate PREDSQ for the current line search and maintain PRESAV.
        !
        temp = step*(one - step)*distsq
        predsq = vlag*vlag*(vlag*vlag + ha*temp*temp)
        if (predsq > presav) then
            presav = predsq
            ksav = k
            stpsav = step
            ibdsav = isbd
        end if
    end do
    !
    !     Construct XNEW in a way that satisfies the bound constraints exact
    !
    do i = 1, N
        temp = XOPT(i) + stpsav*(XPT(ksav, i) - XOPT(i))
        XNEW(i) = max(SL(i), min(SU(i), temp))
    end do
    if (ibdsav < 0) XNEW(-ibdsav) = SL(-ibdsav)
    if (ibdsav > 0) XNEW(ibdsav) = SU(ibdsav)
    !
    !  Prepare for the iterative method that assembles the constrained Cauchy
    !  step in W. The sum of squares of the fixed components of W is formed in
    !  WFIXSQ, and the free components of W are set to BIGSTP.
    !
    bigstp = Adelt + Adelt
    iflag = 0
100 wfixsq = zero
    ggfree = zero
    do i = 1, N
        W(i) = zero
        tempa = min(XOPT(i) - SL(i), GLAG(i))
        tempb = max(XOPT(i) - SU(i), GLAG(i))
        if (tempa > zero .or. tempb < zero) then
            W(i) = bigstp
            ggfree = ggfree + GLAG(i)**2
        end if
    end do
    if (ggfree == zero) then
        Cauchy = zero
        return
    end if
    !
    !     Investigate whether more components of W can be fixed.
    !
120 temp = Adelt*Adelt - wfixsq
    if (temp > zero) then
        wsqsav = wfixsq
        step = sqrt(temp/ggfree)
        ggfree = zero
        do i = 1, N
            if (W(i) == bigstp) then
                temp = XOPT(i) - step*GLAG(i)
                if (temp <= SL(i)) then
                    W(i) = SL(i) - XOPT(i)
                    wfixsq = wfixsq + W(i)**2
                else if (temp >= SU(i)) then
                    W(i) = SU(i) - XOPT(i)
                    wfixsq = wfixsq + W(i)**2
                else
                    ggfree = ggfree + GLAG(i)**2
                end if
            end if
        end do
        if (wfixsq > wsqsav .and. ggfree > zero) goto 120
    end if
    !
    !     Set the remaining free components of W and all components of XALT,
    !     except that W may be scaled later.
    !
    gw = zero
    do i = 1, N
        if (W(i) == bigstp) then
            W(i) = -step*GLAG(i)
            XALT(i) = max(SL(i), min(SU(i), XOPT(i) + W(i)))
        else if (W(i) == zero) then
            XALT(i) = XOPT(i)
        else if (GLAG(i) > zero) then
            XALT(i) = SL(i)
        else
            XALT(i) = SU(i)
        end if
        gw = gw + GLAG(i)*W(i)
    end do
    !
    !  Set CURV to the curvature of the KNEW-th Lagrange function along W.
    !  Scale W by a factor less than one if that can reduce the modulus of
    !  the Lagrange function at XOPT+W. Set CAUCHY to the final value of
    !  the square of this function.
    !
    curv = zero
    do k = 1, Npt
        temp = zero
        do j = 1, N
            temp = temp + XPT(k, j)*W(j)
        end do
        curv = curv + HCOL(k)*temp*temp
    end do
    if (iflag == 1) curv = -curv
    if (curv > -gw .and. curv < -const*gw) then
        scale = -gw/curv
        do i = 1, N
            temp = XOPT(i) + scale*W(i)
            XALT(i) = max(SL(i), min(SU(i), temp))
        end do
        Cauchy = (half*gw*scale)**2
    else
        Cauchy = (gw + half*curv)**2
    end if
    !
    !  If IFLAG is zero, then XALT is calculated as before after reversing
    !  the sign of GLAG. Thus two XALT vectors become available. The one that
    !  is chosen is the one that gives the larger value of CAUCHY.
    !
    if (iflag == 0) then
        do i = 1, N
            GLAG(i) = -GLAG(i)
            W(N + i) = XALT(i)
        end do
        csave = Cauchy
        iflag = 1
        goto 100
    end if
    if (csave > Cauchy) then
        do i = 1, N
            XALT(i) = W(N + i)
        end do
        Cauchy = csave
    end if

    end subroutine ALTMOV

    subroutine PRELIM (this, N, Npt, X, XL, XU, Rhobeg, Iprint, Maxfun, XBASE, &
        XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, Ndim, SL, SU, Nf, Kopt)
    class(TBOBYQA) this
    real(dp) diff, f, fbeg, half, one, &
        recip, Rhobeg, rhosq, stepa, stepb, temp, two, zero
    integer i, ih, Iprint, ipt, itemp, j, jpt, k, Kopt, &
        Maxfun, N, Ndim, Nf, nfm, nfx, np, Npt
    real(dp) X(*), XL(*), XU(*), XBASE(*), XPT(Npt, *), FVAL(*), GOPT(*), &
        HQ(*), PQ(*), BMAT(Ndim, *), ZMAT(Npt, *), SL(*), SU(*)

    !
    !  The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
    !    same as the corresponding arguments in SUBROUTINE BOBYQA.
    !  The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
    !    are the same as the corresponding arguments in BOBYQB, the elements
    !    of SL and SU being set in BOBYQA.
    !  GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
    !    it is set by PRELIM to the gradient of the quadratic model at XBASE.
    !    If XOPT is nonzero, BOBYQB will change it to its usual value later.
    !  NF is maintained as the number of calls of CALFUN so far.
    !  KOPT will be such that the least calculated value of F so far is at
    !    the point XPT(KOPT,.)+XBASE in the space of the variables.
    !
    !  SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    !  BMAT and ZMAT for the first iteration, and it maintains the values of
    !  NF and KOPT. The vector X is also changed by PRELIM.
    !
    !  Set some constants.
    !
    half = 0.5d0
    one = 1.0d0
    two = 2.0d0
    zero = 0.0d0
    rhosq = Rhobeg*Rhobeg
    recip = one/rhosq
    np = N + 1
    !
    !     Set XBASE to the initial vector of variables, and set the initial
    !     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    !
    XBASE(1:N) = X(1:N)
    XPT(1:Npt, 1:N) = zero
    BMAT(1:Ndim, 1:N) = zero
    HQ(1:(N*np)/2) = zero
    PQ(1:Npt) = zero
    ZMAT(1:Npt, 1:Npt - np) = zero
    !
    !  Begin the initialization procedure. NF becomes one more than the number
    !  of function values so far. The coordinates of the displacement of the
    !  next initial interpolation point from XBASE are set in XPT(NF+1,.).
    !
    Nf = 0
50  nfm = Nf
    nfx = Nf - N
    Nf = Nf + 1
    if (nfm <= 2*N) then
        if (nfm >= 1 .and. nfm <= N) then
            stepa = Rhobeg
            if (SU(nfm) == zero) stepa = -stepa
            XPT(Nf, nfm) = stepa
        else if (nfm > N) then
            stepa = XPT(Nf - N, nfx)
            stepb = -Rhobeg
            if (SL(nfx) == zero) stepb = min(two*Rhobeg, SU(nfx))
            if (SU(nfx) == zero) stepb = max(-two*Rhobeg, SL(nfx))
            XPT(Nf, nfx) = stepb
        end if
    else
        itemp = (nfm - np)/N
        jpt = nfm - itemp*N - N
        ipt = jpt + itemp
        if (ipt > N) then
            itemp = jpt
            jpt = ipt - N
            ipt = itemp
        end if
        XPT(Nf, ipt) = XPT(ipt + 1, ipt)
        XPT(Nf, jpt) = XPT(jpt + 1, jpt)
    end if
    !
    !     Calculate the next value of F. The least function value so far and
    !     its index are required.
    !
    do j = 1, N
        X(j) = min(max(XL(j), XBASE(j) + XPT(Nf, j)), XU(j))
        if (XPT(Nf, j) == SL(j)) X(j) = XL(j)
        if (XPT(Nf, j) == SU(j)) X(j) = XU(j)
    end do
    f = this%funkk(this%obj, X(1:N))
    if (Iprint == 3) then
        write(*, "(/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))") &
            Nf, f, (X(i), i=1, N)
    end if
    FVAL(Nf) = f
    if (Nf == 1) then
        fbeg = f
        Kopt = 1
    else if (f < FVAL(Kopt)) then
        Kopt = Nf
    end if
    !
    !  Set the nonzero initial elements of BMAT and the quadratic model in the
    !  cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
    !  of the NF-th and (NF-N)-th interpolation points may be switched, in
    !  order that the function value at the first of them contributes to the
    !  off-diagonal second derivative terms of the initial quadratic model.
    !
    if (Nf <= 2*N + 1) then
        if (Nf >= 2 .and. Nf <= N + 1) then
            GOPT(nfm) = (f - fbeg)/stepa
            if (Npt < Nf + N) then
                BMAT(1, nfm) = -one/stepa
                BMAT(Nf, nfm) = one/stepa
                BMAT(Npt + nfm, nfm) = -half*rhosq
            end if
        else if (Nf >= N + 2) then
            ih = (nfx*(nfx + 1))/2
            temp = (f - fbeg)/stepb
            diff = stepb - stepa
            HQ(ih) = two*(temp - GOPT(nfx))/diff
            GOPT(nfx) = (GOPT(nfx)*stepb - temp*stepa)/diff
            if (stepa*stepb < zero) then
                if (f < FVAL(Nf - N)) then
                    FVAL(Nf) = FVAL(Nf - N)
                    FVAL(Nf - N) = f
                    if (Kopt == Nf) Kopt = Nf - N
                    XPT(Nf - N, nfx) = stepb
                    XPT(Nf, nfx) = stepa
                end if
            end if
            BMAT(1, nfx) = -(stepa + stepb)/(stepa*stepb)
            BMAT(Nf, nfx) = -half/XPT(Nf - N, nfx)
            BMAT(Nf - N, nfx) = -BMAT(1, nfx) - BMAT(Nf, nfx)
            ZMAT(1, nfx) = sqrt(two)/(stepa*stepb)
            ZMAT(Nf, nfx) = sqrt(half)/rhosq
            ZMAT(Nf - N, nfx) = -ZMAT(1, nfx) - ZMAT(Nf, nfx)
        end if
        !
        !     Set the off-diagonal second derivatives of the Lagrange functions
        !     the initial quadratic model.
        !
    else
        ih = (ipt*(ipt - 1))/2 + jpt
        ZMAT(1, nfx) = recip
        ZMAT(Nf, nfx) = recip
        ZMAT(ipt + 1, nfx) = -recip
        ZMAT(jpt + 1, nfx) = -recip
        temp = XPT(Nf, ipt)*XPT(Nf, jpt)
        HQ(ih) = (fbeg - FVAL(ipt + 1) - FVAL(jpt + 1) + f)/temp
    end if
    if (Nf < Npt .and. Nf < Maxfun) goto 50

    end subroutine PRELIM

    subroutine RESCUE (this, N, Npt, XL, XU, Iprint, Maxfun, XBASE, XPT, &
        FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, Ndim, SL, SU, Nf, Delta, &
        Kopt, VLAG, PTSAUX, PTSID, W)
    class(TBOBYQA) this

    real(dp) beta, bsum, Delta, den, denom, diff, distsq, &
        dsqmin, f, fbase, half, hdiag, one
    real(dp) sfrac, sum, sumpq, temp, vlmxsq, &
        vquad, winc, xp, xq, zero
    integer i, ih, ihp, ihq, ip, Iprint, iq, iw, j, jp, &
        jpn, k, knew, kold, Kopt, kpt, Maxfun, N, Ndim, &
        Nf
    integer np, Npt, nptm, nrem
    real(dp) XL(*), XU(*), XBASE(*), XPT(Npt, *), FVAL(*), XOPT(*), &
        GOPT(*), HQ(*), PQ(*), BMAT(Ndim, *), ZMAT(Npt, *), SL(*), SU(*), &
        VLAG(*), PTSAUX(2, *), PTSID(*), W(*)
    !  The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
    !    GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
    !    the corresponding arguments of BOBYQB on the entry to RESCUE.
    !  NF is maintained as the number of calls of CALFUN so far, except that
    !    NF is set to -1 if the value of MAXFUN prevents further progress.
    !  KOPT is maintained so that FVAL(KOPT) is the least calculated function
    !    value. Its correct value must be given on entry. It is updated if a
    !    new least function value is found, but the corresponding changes to
    !    XOPT and GOPT have to be made later by the calling program.
    !  DELTA is the current trust region radius.
    !  VLAG is a working space vector that will be used for the values of the
    !    provisional Lagrange functions at each of the interpolation points.
    !    They are part of a product that requires VLAG to be of length NDIM.
    !  PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
    !    PTSAUX(2,J) specify the two positions of provisional interpolation
    !    points when a nonzero step is taken along e_J (the J-th coordinate
    !    direction) through XBASE+XOPT, as specified below. Usually these
    !    steps have length DELTA, but other lengths are chosen if necessary
    !    in order to satisfy the given bounds on the variables.
    !  PTSID is also a working space array. It has NPT components that denote
    !    provisional new positions of the original interpolation points, in
    !    case changes are needed to restore the linear independence of the
    !    interpolation conditions. The K-th point is a candidate for change
    !    if and only if PTSID(K) is nonzero. In this case let p and q be the
    !    integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
    !    and q are both positive, the step from XBASE+XOPT to the new K-th
    !    interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
    !    the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
    !    p = 0, respectively.
    !  The first NDIM+NPT elements of the array W are used for working space.
    !  The final elements of BMAT and ZMAT are set in a well-conditioned way
    !    to the values that are appropriate for the new interpolation points.
    !  The elements of GOPT, HQ and PQ are also revised to the values that are
    !    appropriate to the final quadratic model.
    !
    !     Set some constants.
    !
    half = 0.5d0
    one = 1.0d0
    zero = 0.0d0
    np = N + 1
    sfrac = half/np
    nptm = Npt - np
    !
    !  Shift the interpolation points so that XOPT becomes the origin, and set
    !  the elements of ZMAT to zero. The value of SUMPQ is required in the
    !  updating of HQ below. The squares of the distances from XOPT to the
    !  other interpolation points are set at the end of W. Increments of WINC
    !  may be added later to these squares to balance the consideration of
    !  the choice of point that is going to become current.
    !
    sumpq = zero
    winc = zero
    do k = 1, Npt
        distsq = zero
        do j = 1, N
            XPT(k, j) = XPT(k, j) - XOPT(j)
            distsq = distsq + XPT(k, j)**2
        end do
        sumpq = sumpq + PQ(k)
        W(Ndim + k) = distsq
        winc = max(winc, distsq)
        do j = 1, nptm
            ZMAT(k, j) = zero
        end do
    end do
    !
    !     Update HQ so that HQ and PQ define the second derivatives of the model
    !     after XBASE has been shifted to the trust region centre.
    !
    ih = 0
    do j = 1, N
        W(j) = half*sumpq*XOPT(j)
        do k = 1, Npt
            W(j) = W(j) + PQ(k)*XPT(k, j)
        end do
        do i = 1, j
            ih = ih + 1
            HQ(ih) = HQ(ih) + W(i)*XOPT(j) + W(j)*XOPT(i)
        end do
    end do
    !
    !     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
    !     also set the elements of PTSAUX.
    !
    do j = 1, N
        XBASE(j) = XBASE(j) + XOPT(j)
        SL(j) = SL(j) - XOPT(j)
        SU(j) = SU(j) - XOPT(j)
        XOPT(j) = zero
        PTSAUX(1, j) = min(Delta, SU(j))
        PTSAUX(2, j) = max(-Delta, SL(j))
        if (PTSAUX(1, j) + PTSAUX(2, j) < zero) then
            temp = PTSAUX(1, j)
            PTSAUX(1, j) = PTSAUX(2, j)
            PTSAUX(2, j) = temp
        end if
        if (abs(PTSAUX(2, j)) < half*abs(PTSAUX(1, j))) then
            PTSAUX(2, j) = half*PTSAUX(1, j)
        end if
        do i = 1, Ndim
            BMAT(i, j) = zero
        end do
    end do
    fbase = FVAL(Kopt)
    !
    !     Set the identifiers of the artificial interpolation points that ard
    !     along a coordinate direction from XOPT, and set the corresponding
    !     nonzero elements of BMAT and ZMAT.
    !
    PTSID(1) = sfrac
    do j = 1, N
        jp = j + 1
        jpn = jp + N
        PTSID(jp) = j + sfrac
        if (jpn <= Npt) then
            PTSID(jpn) = real(j, dp)/np + sfrac
            temp = one/(PTSAUX(1, j) - PTSAUX(2, j))
            BMAT(jp, j) = -temp + one/PTSAUX(1, j)
            BMAT(jpn, j) = temp + one/PTSAUX(2, j)
            BMAT(1, j) = -BMAT(jp, j) - BMAT(jpn, j)
            ZMAT(1, j) = sqrt(2.0d0)/abs(PTSAUX(1, j)*PTSAUX(2, j))
            ZMAT(jp, j) = ZMAT(1, j)*PTSAUX(2, j)*temp
            ZMAT(jpn, j) = -ZMAT(1, j)*PTSAUX(1, j)*temp
        else
            BMAT(1, j) = -one/PTSAUX(1, j)
            BMAT(jp, j) = one/PTSAUX(1, j)
            BMAT(j + Npt, j) = -half*PTSAUX(1, j)**2
        end if
    end do
    !
    !     Set any remaining identifiers with their nonzero elements of ZMAT.
    !
    if (Npt >= N + np) then
        do k = 2*np, Npt
            iw = (real(k - np, dp) - half)/N
            ip = k - np - iw*N
            iq = ip + iw
            if (iq > N) iq = iq - N
            PTSID(k) = ip + real(iq, dp)/np + sfrac
            temp = one/(PTSAUX(1, ip)*PTSAUX(1, iq))
            ZMAT(1, k - np) = temp
            ZMAT(ip + 1, k - np) = -temp
            ZMAT(iq + 1, k - np) = -temp
            ZMAT(k, k - np) = temp
        end do
    end if
    nrem = Npt
    kold = 1
    knew = Kopt
    !
    !     Reorder the provisional points in the way that exchanges PTSID(KOL
    !     with PTSID(KNEW).
    !
80  do j = 1, N
        temp = BMAT(kold, j)
        BMAT(kold, j) = BMAT(knew, j)
        BMAT(knew, j) = temp
    end do
    do j = 1, nptm
        temp = ZMAT(kold, j)
        ZMAT(kold, j) = ZMAT(knew, j)
        ZMAT(knew, j) = temp
    end do
    PTSID(kold) = PTSID(knew)
    PTSID(knew) = zero
    W(Ndim + knew) = zero
    nrem = nrem - 1
    if (knew /= Kopt) then
        temp = VLAG(kold)
        VLAG(kold) = VLAG(knew)
        VLAG(knew) = temp
        !
        !  Update the BMAT and ZMAT matrices so that the status of the KNEW-th
        !  interpolation point can be changed from provisional to original. The
        !  return occurs if all the original points are reinstated.
        !  The nonnegative values of W(NDIM+K) are required in the search below.
        !
        call UPDATE (N, Npt, BMAT, ZMAT, Ndim, VLAG, beta, denom, knew, W)
        if (nrem == 0) return
        do k = 1, Npt
            W(Ndim + k) = abs(W(Ndim + k))
        end do
    end if
    !
    !  Pick the index KNEW of an original interpolation point that has not
    !  yet replaced one of the provisional interpolation points, giving
    !  attention to the closeness to XOPT and to previous tries with KNEW.
    !
120 dsqmin = zero
    do k = 1, Npt
        if (W(Ndim + k) > zero) then
            if (dsqmin == zero .or. W(Ndim + k) < dsqmin) then
                knew = k
                dsqmin = W(Ndim + k)
            end if
        end if
    end do
    if (dsqmin == zero) goto 260
    !
    !     Form the W-vector of the chosen original interpolation point.
    !
    do j = 1, N
        W(Npt + j) = XPT(knew, j)
    end do
    do k = 1, Npt
        sum = zero
        if (k == Kopt) then
            continue
        else if (PTSID(k) == zero) then
            do j = 1, N
                sum = sum + W(Npt + j)*XPT(k, j)
            end do
        else
            ip = PTSID(k)
            if (ip > 0) sum = W(Npt + ip)*PTSAUX(1, ip)
            iq = np*PTSID(k) - ip*np
            if (iq > 0) then
                iw = 1
                if (ip == 0) iw = 2
                sum = sum + W(Npt + iq)*PTSAUX(iw, iq)
            end if
        end if
        W(k) = half*sum*sum
    end do
    !
    !     Calculate VLAG and BETA for the required updating of the H matrix
    !     XPT(KNEW,.) is reinstated in the set of interpolation points.
    !
    do k = 1, Npt
        sum = zero
        do j = 1, N
            sum = sum + BMAT(k, j)*W(Npt + j)
        end do
        VLAG(k) = sum
    end do
    beta = zero
    do j = 1, nptm
        sum = zero
        do k = 1, Npt
            sum = sum + ZMAT(k, j)*W(k)
        end do
        beta = beta - sum*sum
        do k = 1, Npt
            VLAG(k) = VLAG(k) + sum*ZMAT(k, j)
        end do
    end do
    bsum = zero
    distsq = zero
    do j = 1, N
        sum = zero
        do k = 1, Npt
            sum = sum + BMAT(k, j)*W(k)
        end do
        jp = j + Npt
        bsum = bsum + sum*W(jp)
        do ip = Npt + 1, Ndim
            sum = sum + BMAT(ip, j)*W(ip)
        end do
        bsum = bsum + sum*W(jp)
        VLAG(jp) = sum
        distsq = distsq + XPT(knew, j)**2
    end do
    beta = half*distsq*distsq + beta - bsum
    VLAG(Kopt) = VLAG(Kopt) + one
    !
    !  KOLD is set to the index of the provisional interpolation point that is
    !  going to be deleted to make way for the KNEW-th original interpolation
    !  point. The choice of KOLD is governed by the avoidance of a small value
    !  of the denominator in the updating calculation of UPDATE.
    !
    denom = zero
    vlmxsq = zero
    do k = 1, Npt
        if (PTSID(k) /= zero) then
            hdiag = zero
            do j = 1, nptm
                hdiag = hdiag + ZMAT(k, j)**2
            end do
            den = beta*hdiag + VLAG(k)**2
            if (den > denom) then
                kold = k
                denom = den
            end if
        end if
        vlmxsq = max(vlmxsq, VLAG(k)**2)
    end do
    if (denom <= 1.0d-2*vlmxsq) then
        W(Ndim + knew) = -W(Ndim + knew) - winc
        goto 120
    end if
    goto 80
    !
    !  When label 260 is reached, all the final positions of the interpolation
    !  points have been chosen although any changes have not been included yet
    !  in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
    !  from the shift of XBASE, the updating of the quadratic model remains to
    !  be done. The following cycle through the new interpolation points begins
    !  by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
    !  except that a RETURN occurs if MAXFUN prohibits another value of F.
    !
260 do kpt = 1, Npt
        if (PTSID(kpt) == zero) continue
        if (Nf >= Maxfun) then
            Nf = -1
            return
        end if
        ih = 0
        do j = 1, N
            W(j) = XPT(kpt, j)
            XPT(kpt, j) = zero
            temp = PQ(kpt)*W(j)
            do i = 1, j
                ih = ih + 1
                HQ(ih) = HQ(ih) + temp*W(i)
            end do
        end do
        PQ(kpt) = zero
        ip = PTSID(kpt)
        iq = np*PTSID(kpt) - ip*np
        if (ip > 0) then
            xp = PTSAUX(1, ip)
            XPT(kpt, ip) = xp
        end if
        if (iq > 0) then
            xq = PTSAUX(1, iq)
            if (ip == 0) xq = PTSAUX(2, iq)
            XPT(kpt, iq) = xq
        end if
        !
        !     Set VQUAD to the value of the current model at the new point.
        !
        vquad = fbase
        if (ip > 0) then
            ihp = (ip + ip*ip)/2
            vquad = vquad + xp*(GOPT(ip) + half*xp*HQ(ihp))
        end if
        if (iq > 0) then
            ihq = (iq + iq*iq)/2
            vquad = vquad + xq*(GOPT(iq) + half*xq*HQ(ihq))
            if (ip > 0) then
                iw = max(ihp, ihq) - abs(ip - iq)
                vquad = vquad + xp*xq*HQ(iw)
            end if
        end if
        do k = 1, Npt
            temp = zero
            if (ip > 0) temp = temp + xp*XPT(k, ip)
            if (iq > 0) temp = temp + xq*XPT(k, iq)
            vquad = vquad + half*PQ(k)*temp*temp
        end do
        !
        !  Calculate F at the new interpolation point, and set DIFF to the factor
        !  that is going to multiply the KPT-th Lagrange function when the model
        !  is updated to provide interpolation to the new function value.
        !
        do i = 1, N
            W(i) = min(max(XL(i), XBASE(i) + XPT(kpt, i)), XU(i))
            if (XPT(kpt, i) == SL(i)) W(i) = XL(i)
            if (XPT(kpt, i) == SU(i)) W(i) = XU(i)
        end do
        Nf = Nf + 1
        f = this%funkk(this%obj, W(1:N))
        if (Iprint == 3) then
            write(*, "(/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))") &
                Nf, f, (W(i), i=1, N)
        end if
        FVAL(kpt) = f
        if (f < FVAL(Kopt)) Kopt = kpt
        diff = f - vquad
        !
        !     Update the quadratic model. The RETURN from the subroutine occurs when
        !     all the new interpolation points are included in the model.
        !
        do i = 1, N
            GOPT(i) = GOPT(i) + diff*BMAT(kpt, i)
        end do
        do k = 1, Npt
            sum = zero
            do j = 1, nptm
                sum = sum + ZMAT(k, j)*ZMAT(kpt, j)
            end do
            temp = diff*sum
            if (PTSID(k) == zero) then
                PQ(k) = PQ(k) + temp
            else
                ip = PTSID(k)
                iq = np*PTSID(k) - ip*np
                ihq = (iq*iq + iq)/2
                if (ip == 0) then
                    HQ(ihq) = HQ(ihq) + temp*PTSAUX(2, iq)**2
                else
                    ihp = (ip*ip + ip)/2
                    HQ(ihp) = HQ(ihp) + temp*PTSAUX(1, ip)**2
                    if (iq > 0) then
                        HQ(ihq) = HQ(ihq) + temp*PTSAUX(1, iq)**2
                        iw = max(ihp, ihq) - abs(iq - ip)
                        HQ(iw) = HQ(iw) + temp*PTSAUX(1, ip)*PTSAUX(1, iq)
                    end if
                end if
            end if
        end do
        PTSID(kpt) = zero
    end do

    end subroutine RESCUE

    subroutine TRSBOX (N, Npt, XPT, XOPT, GOPT, HQ, PQ, SL, SU, Delta, &
        XNEW, D, GNEW, XBDI, S, HS, HRED, Dsq, Crvmin)
    real(dp) angbd, angt, beta, blen, Crvmin, cth, delsq, &
        Delta, dhd, dhs, dredg, dredsq, ds, Dsq, ggsav, &
        gredsq, half
    real(dp) one, onemin, qred, ratio, &
        rdnext, rdprev, redmax, rednew, redsav, resid, &
        sdec, shs, sqstp, sredg
    real(dp) ssq, stepsq, sth, stplen, temp, tempa, tempb, &
        xsav, xsum, zero
    integer i, iact, ih, isav, itcsav, iterc, itermax, iu, j, &
        k, N, nact, Npt
    real(dp) XPT(Npt, *), XOPT(*), GOPT(*), HQ(*), PQ(*), SL(*), SU(*), &
        XNEW(*), D(*), GNEW(*), XBDI(*), S(*), HS(*), HRED(*)
    !
    !  The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
    !    meanings as the corresponding arguments of BOBYQB.
    !  DELTA is the trust region radius for the present calculation, which
    !    seeks a small value of the quadratic model within distance DELTA of
    !    XOPT subject to the bounds on the variables.
    !  XNEW will be set to a new vector of variables that is approximately
    !    the one that minimizes the quadratic model within the trust region
    !    subject to the SL and SU constraints on the variables. It satisfies
    !    as equations the bounds that become active during the calculation.
    !  D is the calculated trial step from XOPT, generated iteratively from an
    !    initial value of zero. Thus XNEW is XOPT+D after the final iteration.
    !  GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
    !    when D is updated.
    !  XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
    !    set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
    !    I-th variable has become fixed at a bound, the bound being SL(I) or
    !    SU(I) in the case XBDI(I) = -1.0 or XBDI(I) = 1.0, respectively. This
    !    information is accumulated during the construction of XNEW.
    !  The arrays S, HS and HRED are also used for working space. They hold the
    !    current search direction, and the changes in the gradient of Q along S
    !    and the reduced D, respectively, where the reduced D is the same as D,
    !    except that the components of the fixed variables are zero.
    !  DSQ will be set to the square of the length of XNEW-XOPT.
    !  CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
    !    it is set to the least curvature of H that occurs in the conjugate
    !    gradient searches that are not restricted by any constraints. The
    !    value CRVMIN=-1.0D0 is set, however, if all of these searches are
    !    constrained.
    !
    !  A version of the truncated conjugate gradient is applied. If a line
    !  search is restricted by a constraint, then the procedure is restarted,
    !  the values of the variables that are at their bounds being fixed. If
    !  the trust region boundary is reached, then further changes may be made
    !  to D, each one being in the two dimensional space that is spanned
    !  by the current D and the gradient of Q at XOPT+D, staying on the trust
    !  region boundary. Termination occurs when the reduction in Q seems to
    !  be close to the greatest reduction that can be achieved.
    !
    !     Set some constants.
    !
    half = 0.5d0
    one = 1.0d0
    onemin = -1.0d0
    zero = 0.0d0
    !
    !  The sign of GOPT(I) gives the sign of the change to the I-th variable
    !  that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
    !  or not to fix the I-th variable at one of its bounds initially, with
    !  NACT being set to the number of fixed variables. D and GNEW are also
    !  set for the first iteration. DELSQ is the upper bound on the sum of
    !  squares of the free variables. QRED is the reduction in Q so far.
    !
    iterc = 0
    nact = 0
    sqstp = zero
    do i = 1, N
        XBDI(i) = zero
        if (XOPT(i) <= SL(i)) then
            if (GOPT(i) >= zero) XBDI(i) = onemin
        else if (XOPT(i) >= SU(i)) then
            if (GOPT(i) <= zero) XBDI(i) = one
        end if
        if (XBDI(i) /= zero) nact = nact + 1
        D(i) = zero
        GNEW(i) = GOPT(i)
    end do
    delsq = Delta*Delta
    qred = zero
    Crvmin = onemin
    !
    !  Set the next search direction of the conjugate gradient method. It is
    !  the steepest descent direction initially and when the iterations are
    !  restarted because a variable has just been fixed by a bound, and of
    !  course the components of the fixed variables are zero. ITERMAX is an
    !  upper bound on the indices of the conjugate gradient iterations.
    !
20  beta = zero
30  stepsq = zero
    do i = 1, N
        if (XBDI(i) /= zero) then
            S(i) = zero
        else if (beta == zero) then
            S(i) = -GNEW(i)
        else
            S(i) = beta*S(i) - GNEW(i)
        end if
        stepsq = stepsq + S(i)**2
    end do
    if (stepsq == zero) goto 190
    if (beta == zero) then
        gredsq = stepsq
        itermax = iterc + N - nact
    end if
    if (gredsq*delsq <= 1.0d-4*qred*qred) goto 190
    !
    !  Multiply the search direction by the second derivative matrix of Q and
    !  calculate some scalars for the choice of steplength. Then set BLEN to
    !  the length of the the step to the trust region boundary and STPLEN to
    !  the steplength, ignoring the simple bounds.
    !
    goto 210
50  resid = delsq
    ds = zero
    shs = zero
    do i = 1, N
        if (XBDI(i) == zero) then
            resid = resid - D(i)**2
            ds = ds + S(i)*D(i)
            shs = shs + S(i)*HS(i)
        end if
    end do
    if (resid <= zero) goto 90
    temp = sqrt(stepsq*resid + ds*ds)
    if (ds < zero) then
        blen = (temp - ds)/stepsq
    else
        blen = resid/(temp + ds)
    end if
    stplen = blen
    if (shs > zero) then
        stplen = min(blen, gredsq/shs)
    end if

    !
    !     Reduce STPLEN if necessary in order to preserve the simple bounds,
    !     letting IACT be the index of the new constrained variable.
    !
    iact = 0
    do i = 1, N
        if (S(i) /= zero) then
            xsum = XOPT(i) + D(i)
            if (S(i) > zero) then
                temp = (SU(i) - xsum)/S(i)
            else
                temp = (SL(i) - xsum)/S(i)
            end if
            if (temp < stplen) then
                stplen = temp
                iact = i
            end if
        end if
    end do
    !
    !     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in
    !
    sdec = zero
    if (stplen > zero) then
        iterc = iterc + 1
        temp = shs/stepsq
        if (iact == 0 .and. temp > zero) then
            Crvmin = min(Crvmin, temp)
            if (Crvmin == onemin) Crvmin = temp
        end if
        ggsav = gredsq
        gredsq = zero
        do i = 1, N
            GNEW(i) = GNEW(i) + stplen*HS(i)
            if (XBDI(i) == zero) gredsq = gredsq + GNEW(i)**2
            D(i) = D(i) + stplen*S(i)
        end do
        sdec = max(stplen*(ggsav - half*stplen*shs), zero)
        qred = qred + sdec
    end if
    !
    !     Restart the conjugate gradient method if it has hit a new bound.
    !
    if (iact > 0) then
        nact = nact + 1
        XBDI(iact) = one
        if (S(iact) < zero) XBDI(iact) = onemin
        delsq = delsq - D(iact)**2
        if (delsq <= zero) goto 90
        goto 20
    end if
    !
    !     If STPLEN is less than BLEN, then either apply another conjugate
    !     gradient iteration or RETURN.
    !
    if (stplen < blen) then
        if (iterc == itermax) goto 190
        if (sdec <= 0.01d0*qred) goto 190
        beta = gredsq/ggsav
        goto 30
    end if
90  Crvmin = zero
    !
    !  Prepare for the alternative iteration by calculating some scalars
    !  and by multiplying the reduced D by the second derivative matrix of
    !  Q, where S holds the reduced D in the call of GGMULT.
    !
100 if (nact >= N - 1) goto 190
    dredsq = zero
    dredg = zero
    gredsq = zero
    do i = 1, N
        if (XBDI(i) == zero) then
            dredsq = dredsq + D(i)**2
            dredg = dredg + D(i)*GNEW(i)
            gredsq = gredsq + GNEW(i)**2
            S(i) = D(i)
        else
            S(i) = zero
        end if
    end do
    itcsav = iterc
    goto 210
    !
    !     Let the search direction S be a linear combination of the reduced
    !     and the reduced G that is orthogonal to the reduced D.
    !
120 iterc = iterc + 1
    temp = gredsq*dredsq - dredg*dredg
    if (temp <= 1.0d-4*qred*qred) goto 190
    temp = sqrt(temp)
    do i = 1, N
        if (XBDI(i) == zero) then
            S(i) = (dredg*D(i) - dredsq*GNEW(i))/temp
        else
            S(i) = zero
        end if
    end do
    sredg = -temp
    !
    !  By considering the simple bounds on the variables, calculate an upper
    !  bound on the tangent of half the angle of the alternative iteration,
    !  namely ANGBD, except that, if already a free variable has reached a
    !  bound, there is a branch back to label 100 after fixing that variable.
    !
    angbd = one
    iact = 0
    do i = 1, N
        if (XBDI(i) == zero) then
            tempa = XOPT(i) + D(i) - SL(i)
            tempb = SU(i) - XOPT(i) - D(i)
            if (tempa <= zero) then
                nact = nact + 1
                XBDI(i) = onemin
                goto 100
            else if (tempb <= zero) then
                nact = nact + 1
                XBDI(i) = one
                goto 100
            end if
            ratio = one
            ssq = D(i)**2 + S(i)**2
            temp = ssq - (XOPT(i) - SL(i))**2
            if (temp > zero) then
                temp = sqrt(temp) - S(i)
                if (angbd*temp > tempa) then
                    angbd = tempa/temp
                    iact = i
                    xsav = onemin
                end if
            end if
            temp = ssq - (SU(i) - XOPT(i))**2
            if (temp > zero) then
                temp = sqrt(temp) + S(i)
                if (angbd*temp > tempb) then
                    angbd = tempb/temp
                    iact = i
                    xsav = one
                end if
            end if
        end if
    end do
    !
    !     Calculate HHD and some curvatures for the alternative iteration.
    !
    goto 210
150 shs = zero
    dhs = zero
    dhd = zero
    do i = 1, N
        if (XBDI(i) == zero) then
            shs = shs + S(i)*HS(i)
            dhs = dhs + D(i)*HS(i)
            dhd = dhd + D(i)*HRED(i)
        end if
    end do
    !
    !  Seek the greatest reduction in Q for a range of equally spaced values
    !  of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
    !  the alternative iteration.
    !
    redmax = zero
    isav = 0
    redsav = zero
    iu = 17.0d0*angbd + 3.1d0
    do i = 1, iu
        angt = angbd*real(i, dp)/iu
        sth = (angt + angt)/(one + angt*angt)
        temp = shs + angt*(angt*dhd - dhs - dhs)
        rednew = sth*(angt*dredg - sredg - half*sth*temp)
        if (rednew > redmax) then
            redmax = rednew
            isav = i
            rdprev = redsav
        else if (i == isav + 1) then
            rdnext = rednew
        end if
        redsav = rednew
    end do
    !
    !  Return if the reduction is zero. Otherwise, set the sine and cosine
    !  of the angle of the alternative iteration, and calculate SDEC.
    !
    if (isav == 0) goto 190
    if (isav < iu) then
        temp = (rdnext - rdprev)/(redmax + redmax - rdprev - rdnext)
        angt = angbd*(isav + half*temp)/iu
    end if
    cth = (one - angt*angt)/(one + angt*angt)
    sth = (angt + angt)/(one + angt*angt)
    temp = shs + angt*(angt*dhd - dhs - dhs)
    sdec = sth*(angt*dredg - sredg - half*sth*temp)
    if (sdec <= zero) goto 190
    !
    !  Update GNEW, D and HRED. If the angle of the alternative iteration
    !  is restricted by a bound on a free variable, that variable is fixed
    !  at the bound.
    !
    dredg = zero
    gredsq = zero
    do i = 1, N
        GNEW(i) = GNEW(i) + (cth - one)*HRED(i) + sth*HS(i)
        if (XBDI(i) == zero) then
            D(i) = cth*D(i) + sth*S(i)
            dredg = dredg + D(i)*GNEW(i)
            gredsq = gredsq + GNEW(i)**2
        end if
        HRED(i) = cth*HRED(i) + sth*HS(i)
    end do
    qred = qred + sdec
    if (iact > 0 .and. isav == iu) then
        nact = nact + 1
        XBDI(iact) = xsav
        goto 100
    end if
    !  If SDEC is sufficiently small, then RETURN after setting XNEW to
    !  XOPT+D, giving careful attention to the bounds.
    !
    if (sdec > 0.01d0*qred) goto 120
190 Dsq = zero
    do i = 1, N
        XNEW(i) = max(min(XOPT(i) + D(i), SU(i)), SL(i))
        if (XBDI(i) == onemin) XNEW(i) = SL(i)
        if (XBDI(i) == one) XNEW(i) = SU(i)
        D(i) = XNEW(i) - XOPT(i)
        Dsq = Dsq + D(i)**2
    end do
    return

    !  The following instructions multiply the current S-vector by the second
    !  derivative matrix of the quadratic model, putting the product in HS.
    !  They are reached from three different parts of the software above and
    !  they can be regarded as an external subroutine.
    !
210 ih = 0
    do j = 1, N
        HS(j) = zero
        do i = 1, j
            ih = ih + 1
            if (i < j) HS(j) = HS(j) + HQ(ih)*S(i)
            HS(i) = HS(i) + HQ(ih)*S(j)
        end do
    end do
    do k = 1, Npt
        if (PQ(k) /= zero) then
            temp = zero
            do j = 1, N
                temp = temp + XPT(k, j)*S(j)
            end do
            temp = temp*PQ(k)
            do i = 1, N
                HS(i) = HS(i) + temp*XPT(k, i)
            end do
        end if
    end do
    if (Crvmin /= zero) goto 50
    if (iterc > itcsav) goto 150
    do i = 1, N
        HRED(i) = HS(i)
    end do
    goto 120

    end subroutine TRSBOX

    subroutine UPDATE (N, Npt, BMAT, ZMAT, Ndim, VLAG, Beta, Denom, Knew, W)
    real(dp) alpha, Beta, Denom, one, tau, temp, tempa, &
        tempb, zero, ztest
    integer i, j, jl, jp, k, Knew, N, Ndim, Npt, nptm
    real(dp) BMAT(Ndim, *), ZMAT(Npt, *), VLAG(*), W(*)
    !
    !  The arrays BMAT and ZMAT are updated, as required by the new position
    !  of the interpolation point that has the index KNEW. The vector VLAG has
    !  N+NPT components, set on entry to the first NPT and last N components
    !  of the product Hw in equation (4.11) of the Powell (2006) paper on
    !  NEWUOA. Further, BETA is set on entry to the value of the parameter
    !  with that name, and DENOM is set to the denominator of the updating
    !  formula. Elements of ZMAT may be treated as zero if their moduli are
    !  at most ZTEST. The first NDIM elements of W are used for working space.
    !
    !     Set some constants.
    !
    one = 1.0d0
    zero = 0.0d0
    nptm = Npt - N - 1
    ztest = maxval(abs(ZMAT(1:Npt, 1:nptm)))
    ztest = 1.0d-20*ztest
    !
    !     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    !
    jl = 1
    do j = 2, nptm
        if (abs(ZMAT(Knew, j)) > ztest) then
            temp = sqrt(ZMAT(Knew, 1)**2 + ZMAT(Knew, j)**2)
            tempa = ZMAT(Knew, 1)/temp
            tempb = ZMAT(Knew, j)/temp
            do i = 1, Npt
                temp = tempa*ZMAT(i, 1) + tempb*ZMAT(i, j)
                ZMAT(i, j) = tempa*ZMAT(i, j) - tempb*ZMAT(i, 1)
                ZMAT(i, 1) = temp
            end do
        end if
        ZMAT(Knew, j) = zero
    end do
    !
    !     Put the first NPT components of the KNEW-th column of HLAG into W,
    !     and calculate the parameters of the updating formula.
    !
    W(1:Npt) = ZMAT(Knew, 1)*ZMAT(1:Npt, 1)
    alpha = W(Knew)
    tau = VLAG(Knew)
    VLAG(Knew) = VLAG(Knew) - one
    !
    !     Complete the updating of ZMAT.
    !
    temp = sqrt(Denom)
    tempb = ZMAT(Knew, 1)/temp
    tempa = tau/temp
    do i = 1, Npt
        ZMAT(i, 1) = tempa*ZMAT(i, 1) - tempb*VLAG(i)
    end do
    !
    !     Finally, update the matrix BMAT.
    !
    do j = 1, N
        jp = Npt + j
        W(jp) = BMAT(Knew, j)
        tempa = (alpha*VLAG(jp) - tau*W(jp))/Denom
        tempb = (-Beta*W(jp) - tau*VLAG(jp))/Denom
        do i = 1, jp
            BMAT(i, j) = BMAT(i, j) + tempa*VLAG(i) + tempb*W(i)
            if (i > Npt) BMAT(jp, i - Npt) = BMAT(i, j)
        end do
    end do

    end subroutine UPDATE

    ! *****************************************************************************************
    ! >
    !  This subroutine seeks the least value of a function of many variables,
    !  by a trust region method that forms quadratic models by interpolation.
    !  There can be some freedom in the interpolation conditions, which is
    !  taken up by minimizing the Frobenius norm of the change to the second
    !  derivative of the quadratic model, beginning with a zero matrix.

    logical function newuoa (this, obj, funk, n, npt, x, rhobeg, rhoend, iprint, maxfun)
    class(TNEWUOA) this
    class(*), target :: obj
    procedure(obj_vec_function) :: funk ! funk(obj, x), unwrapping obj with select type

    integer, intent(in) :: n       !! the number of variables. must be at least 2.
    integer, intent(in) :: npt     !! The number of interpolation conditions.
    !! Its value must be in the interval `[N+2,(N+1)(N+2)/2]`.
    real(dp), dimension(*), intent(inout) :: x       !! Initial values must be set in X(1),X(2),...,X(N). They
    !! will be changed to the values that give the least calculated F.
    real(dp), intent(in) :: rhobeg  !! RHOBEG and RHOEND must be set to the initial and final values
    !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
    !! RHOBEG should be about one tenth of the greatest expected change to a
    !! variable, and RHOEND should indicate the accuracy that is required in
    !! the final values of the variables.
    real(dp), intent(in) :: rhoend  !! RHOBEG and RHOEND must be set to the initial and final values
    !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
    !! RHOBEG should be about one tenth of the greatest expected change to a
    !! variable, and RHOEND should indicate the accuracy that is required in
    !! the final values of the variables.
    integer, intent(in)  :: iprint  !! IPRINT should be set to 0, 1, 2 or 3, which controls the
    !! amount of printing. Specifically, there is no output if IPRINT=0 and
    !! there is output only at the return if IPRINT=1. Otherwise, each new
    !! value of RHO is printed, with the best vector of variables so far and
    !! the corresponding value of the objective function. Further, each new
    !! value of F with its variables are output if IPRINT=3.
    integer, intent(in)  :: maxfun  !! an upper bound on the number of calls of CALFUN.
    !! for the variables `X(1),X(2),...,X(N)`.

    real(dp), dimension(:), allocatable :: w
    integer :: np, nptm, ndim, ixb, ixo, ixn, ixp, ifv, igq, ihq, ipq, ibmat, izmat, id, ivl, iw

    ! Partition the working space array, so that different parts of it can be
    ! treated separately by the subroutine that performs the main calculation.

    this%funkk => funk
    this%obj => obj
    this%Last_bestfit = huge(1._dp)

    np = n + 1
    nptm = npt - np
    if (npt < n + 2 .or. npt > ((n + 2)*np)/2) then
        write(*, *) 'Return from NEWUOA because NPT is not in the required interval'
        newuoa = .false.
        return
    end if

    ! The array W will be used for working space
    allocate(w((npt + 13)*(npt + n) + 3*n*(n + 3)/2))

    ndim = npt + n
    ixb = 1
    ixo = ixb + n
    ixn = ixo + n
    ixp = ixn + n
    ifv = ixp + n*npt
    igq = ifv + npt
    ihq = igq + n
    ipq = ihq + (n*np)/2
    ibmat = ipq + npt
    izmat = ibmat + ndim*n
    id = izmat + npt*nptm
    ivl = id + n
    iw = ivl + ndim

    ! The above settings provide a partition of W for subroutine NEWUOB.
    ! The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
    ! W plus the space that is needed by the last array of NEWUOB.

    newuoa = this%NEWUOB (n, npt, x, rhobeg, rhoend, iprint, maxfun, w(ixb), w(ixo), w(ixn), &
        w(ixp), w(ifv), w(igq), w(ihq), w(ipq), w(ibmat), w(izmat), ndim, &
        w(id), w(ivl), w(iw))

    deallocate(w)

    end function newuoa
    ! *****************************************************************************************

    logical function newuob (this, N, Npt, x, Rhobeg, Rhoend, Iprint, Maxfun, xbase, xopt, &
        xnew, xpt, fval, gq, hq, pq, bmat, zmat, Ndim, d, vlag, w)
    class(TNEWUOA) this

    real(dp) alpha, beta, bsum, crvmin, delta, detrat, &
        diff, diffa, diffb, diffc, distsq, dnorm, dsq, &
        dstep, dx, f, fbeg, fopt
    real(dp) fsave, gisq, gqsq, half, hdiag, one, ratio, recip, reciq, &
        rho, Rhobeg, Rhoend, rhosq, sum, suma, sumb
    real(dp) sumz, temp, tempq, tenth, vquad, xjpt, xoptsq, zero, xipt
    integer i, idz, ih, ip, Iprint, ipt, itemp, itest, &
        j, jp, jpt, k, knew, kopt, ksave, ktemp, Maxfun, &
        N, Ndim
    integer nf, nfm, nfmm, nfsav, nftest, nh, np, Npt, nptm
    real(dp) x (*), xbase (*), xopt (*), xnew (*), xpt (Npt, *), fval (*), gq (*), hq(*), &
        pq (*), bmat (Ndim, *), zmat (Npt, *), d (*), vlag (*), w (*)
    !
    !     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
    !       to the corresponding arguments in SUBROUTINE NEWUOA.
    !     XBASE will hold a shift of origin that should reduce the contributions
    !       from rounding errors to values of the model and Lagrange functions.
    !     XOPT will be set to the displacement from XBASE of the vector of
    !       variables that provides the least calculated F so far.
    !     XNEW will be set to the displacement from XBASE of the vector of
    !       variables for the current calculation of F.
    !     XPT will contain the interpolation point coordinates relative to XBASE.
    !     FVAL will hold the values of F at the interpolation points.
    !     GQ will hold the gradient of the quadratic model at XBASE.
    !     HQ will hold the explicit second derivatives of the quadratic model.
    !     PQ will contain the parameters of the implicit second derivatives of
    !       the quadratic model.
    !     BMAT will hold the last N columns of H.
    !     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
    !       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
    !       the elements of DZ are plus or minus one, as specified by IDZ.
    !     NDIM is the first dimension of BMAT and has the value NPT+N.
    !     D is reserved for trial steps from XOPT.
    !     VLAG will contain the values of the Lagrange functions at a new point X.
    !       They are part of a product that requires VLAG to be of length NDIM.
    !     The array W will be used for working space. Its length must be at least
    !       10*NDIM = 10*(NPT+N).
    !
    !     Set some constants.
    !
    newuob = .true.
    half = 0.5_dp
    one = 1.0_dp
    tenth = 0.1_dp
    zero = 0.0_dp
    np = N + 1
    nh = (N*np)/2
    nptm = Npt - np
    nftest = max (Maxfun, 1)
    !
    !     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    !
    xbase(1:N) = x(1:N)
    xpt(1:Npt, 1:N) = zero
    bmat(1:Ndim, 1:N) = zero
    hq(1:nh) = zero
    pq(1:Npt) = zero
    zmat(1:Npt, 1:nptm) = zero
    !
    !     Begin the initialization procedure. NF becomes one more than the number
    !     of function values so far. The coordinates of the displacement of the
    !     next initial interpolation point from XBASE are set in XPT(NF,.).
    !
    rhosq = Rhobeg*Rhobeg
    recip = one/rhosq
    reciq = sqrt (half)/rhosq
    nf = 0
50  nfm = nf
    nfmm = nf - N
    nf = nf + 1
    if (nfm <= 2*N) then
        if (nfm >= 1 .and. nfm <= N) then
            xpt (nf, nfm) = Rhobeg
        else if (nfm > N) then
            xpt (nf, nfmm) = -Rhobeg
        end if
    else
        itemp = (nfmm - 1)/N
        jpt = nfm - itemp*N - N
        ipt = jpt + itemp
        if (ipt > N) then
            itemp = jpt
            jpt = ipt - N
            ipt = itemp
        end if
        xipt = Rhobeg
        if (fval(ipt + np) < fval(ipt + 1)) xipt = -xipt
        xjpt = Rhobeg
        if (fval(jpt + np) < fval(jpt + 1)) xjpt = -xjpt
        xpt (nf, ipt) = xipt
        xpt (nf, jpt) = xjpt
    end if
    !
    !     Calculate the next value of F, label 70 being reached immediately
    !     after this calculation. The least function value so far and its index
    !     are required.
    !
    do j = 1, N
        x (j) = xpt (nf, j) + xbase (j)
    end do
    goto 310
70  fval (nf) = f
    if (nf == 1) then
        fbeg = f
        fopt = f
        kopt = 1
    else if (f < fopt) then
        fopt = f
        kopt = nf
    end if
    !
    !     Set the nonzero initial elements of BMAT and the quadratic model in
    !     the cases when NF is at most 2*N+1.
    !
    if (nfm <= 2*N) then
        if (nfm >= 1 .and. nfm <= N) then
            gq (nfm) = (f - fbeg)/Rhobeg
            if (Npt < nf + N) then
                bmat (1, nfm) = -one/Rhobeg
                bmat (nf, nfm) = one/Rhobeg
                bmat (Npt + nfm, nfm) = -half*rhosq
            end if
        else if (nfm > N) then
            bmat (nf - N, nfmm) = half/Rhobeg
            bmat (nf, nfmm) = -half/Rhobeg
            zmat (1, nfmm) = -reciq - reciq
            zmat (nf - N, nfmm) = reciq
            zmat (nf, nfmm) = reciq
            ih = (nfmm*(nfmm + 1))/2
            temp = (fbeg - f)/Rhobeg
            hq (ih) = (gq(nfmm) - temp)/Rhobeg
            gq (nfmm) = half*(gq(nfmm) + temp)
        end if
        !
        !     Set the off-diagonal second derivatives of the Lagrange functions and
        !     the initial quadratic model.
        !
    else
        ih = (ipt*(ipt - 1))/2 + jpt
        if (xipt < zero) ipt = ipt + N
        if (xjpt < zero) jpt = jpt + N
        zmat (1, nfmm) = recip
        zmat (nf, nfmm) = recip
        zmat (ipt + 1, nfmm) = -recip
        zmat (jpt + 1, nfmm) = -recip
        hq (ih) = (fbeg - fval(ipt + 1) - fval(jpt + 1) + f)/(xipt*xjpt)
    end if
    if (nf < Npt) goto 50
    !
    !     Begin the iterative procedure, because the initial model is complete.
    !
    rho = Rhobeg
    delta = rho
    idz = 1
    diffa = zero
    diffb = zero
    itest = 0
    xoptsq = zero
    do i = 1, N
        xopt (i) = xpt (kopt, i)
        xoptsq = xoptsq + xopt (i)**2
    end do
90  nfsav = nf
    !
    !     Generate the next trust region step and test its length. Set KNEW
    !     to -1 if the purpose of the next F will be to improve the model.
    !
100 knew = 0
    call trsapp (N, Npt, xopt, xpt, gq, hq, pq, delta, d, w, w(np), w(np + N), &
        w(np + 2*N), crvmin)
    dsq = zero
    do i = 1, N
        dsq = dsq + d (i)**2
    end do
    dnorm = min (delta, sqrt(dsq))
    if (dnorm < half*rho) then
        knew = -1
        delta = tenth*delta
        ratio = -1.0_dp
        if (delta <= 1.5_dp*rho) delta = rho
        if (nf <= nfsav + 2) goto 460
        temp = 0.125_dp*crvmin*rho*rho
        if (temp <= max(diffa, diffb, diffc)) goto 460
        goto 490
    end if
    !
    !     Shift XBASE if XOPT may be too far from XBASE. First make the changes
    !     to BMAT that do not depend on ZMAT.
    !
120 if (dsq <= 1.0e-3_dp*xoptsq) then
        tempq = 0.25_dp*xoptsq
        do k = 1, Npt
            sum = zero
            do i = 1, N
                sum = sum + xpt (k, i)*xopt (i)
            end do
            temp = pq (k)*sum
            sum = sum - half*xoptsq
            w (Npt + k) = sum
            do i = 1, N
                gq (i) = gq (i) + temp*xpt (k, i)
                xpt (k, i) = xpt (k, i) - half*xopt (i)
                vlag (i) = bmat (k, i)
                w (i) = sum*xpt (k, i) + tempq*xopt (i)
                ip = Npt + i
                do j = 1, i
                    bmat (ip, j) = bmat (ip, j) + vlag (i)*w (j) + w (i)*vlag (j)
                end do
            end do
        end do
        !
        !     Then the revisions of BMAT that depend on ZMAT are calculated.
        !
        do k = 1, nptm
            sumz = zero
            do i = 1, Npt
                sumz = sumz + zmat (i, k)
                w (i) = w (Npt + i)*zmat (i, k)
            end do
            do j = 1, N
                sum = tempq*sumz*xopt (j)
                do i = 1, Npt
                    sum = sum + w (i)*xpt (i, j)
                end do
                vlag (j) = sum
                if (k < idz) sum = -sum
                do i = 1, Npt
                    bmat (i, j) = bmat (i, j) + sum*zmat (i, k)
                end do
            end do
            do i = 1, N
                ip = i + Npt
                temp = vlag (i)
                if (k < idz) temp = -temp
                do j = 1, i
                    bmat (ip, j) = bmat (ip, j) + temp*vlag (j)
                end do
            end do
        end do
        !
        !     The following instructions complete the shift of XBASE, including
        !     the changes to the parameters of the quadratic model.
        !
        ih = 0
        do j = 1, N
            w (j) = zero
            do k = 1, Npt
                w (j) = w (j) + pq (k)*xpt (k, j)
                xpt (k, j) = xpt (k, j) - half*xopt (j)
            end do
            do i = 1, j
                ih = ih + 1
                if (i < j) gq (j) = gq (j) + hq (ih)*xopt (i)
                gq (i) = gq (i) + hq (ih)*xopt (j)
                hq (ih) = hq (ih) + w (i)*xopt (j) + xopt (i)*w (j)
                bmat (Npt + i, j) = bmat (Npt + j, i)
            end do
        end do
        do j = 1, N
            xbase (j) = xbase (j) + xopt (j)
            xopt (j) = zero
        end do
        xoptsq = zero
    end if
    !
    !     Pick the model step if KNEW is positive. A different choice of D
    !     may be made later, if the choice of D by BIGLAG causes substantial
    !     cancellation in DENOM.
    !
    if (knew > 0) then
        call biglag (N, Npt, xopt, xpt, bmat, zmat, idz, Ndim, knew, dstep, d, alpha, &
            vlag, vlag(Npt + 1), w, w(np), w(np + N))
    end if
    !
    !     Calculate VLAG and BETA for the current choice of D. The first NPT
    !     components of W_check will be held in W.
    !
    do k = 1, Npt
        suma = zero
        sumb = zero
        sum = zero
        do j = 1, N
            suma = suma + xpt (k, j)*d (j)
            sumb = sumb + xpt (k, j)*xopt (j)
            sum = sum + bmat (k, j)*d (j)
        end do
        w (k) = suma*(half*suma + sumb)
        vlag (k) = sum
    end do
    beta = zero
    do k = 1, nptm
        sum = zero
        do i = 1, Npt
            sum = sum + zmat (i, k)*w (i)
        end do
        if (k < idz) then
            beta = beta + sum*sum
            sum = -sum
        else
            beta = beta - sum*sum
        end if
        do i = 1, Npt
            vlag (i) = vlag (i) + sum*zmat (i, k)
        end do
    end do
    bsum = zero
    dx = zero
    do j = 1, N
        sum = zero
        do i = 1, Npt
            sum = sum + w (i)*bmat (i, j)
        end do
        bsum = bsum + sum*d (j)
        jp = Npt + j
        do k = 1, N
            sum = sum + bmat (jp, k)*d (k)
        end do
        vlag (jp) = sum
        bsum = bsum + sum*d (j)
        dx = dx + d (j)*xopt (j)
    end do
    beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta - bsum
    vlag (kopt) = vlag (kopt) + one
    !
    !     If KNEW is positive and if the cancellation in DENOM is unacceptable,
    !     then BIGDEN calculates an alternative model step, XNEW being used for
    !     working space.
    !
    if (knew > 0) then
        temp = one + alpha*beta/vlag (knew)**2
        if (abs(temp) <= 0.8_dp) then
            call bigden (N, Npt, xopt, xpt, bmat, zmat, idz, Ndim, kopt, knew, d, w, &
                vlag, beta, xnew, w(Ndim + 1), w(6*Ndim + 1))
        end if
    end if
    !
    !     Calculate the next value of the objective function.
    !
290 do i = 1, N
        xnew (i) = xopt (i) + d (i)
        x (i) = xbase (i) + xnew (i)
    end do
    nf = nf + 1
310 if (nf > nftest) then
        nf = nf - 1
        if (Iprint > 0) write(*, *) "Return from NEWUOA because CALFUN has been called MAXFUN times."
        newuob = .false.
        goto 530
    end if
    f = this%funkk(this%obj, x(1:N))
    if (Iprint == 3) then
        write(*, "(/4X,'Function number',I6,'    F =',1PD18.10,'    The corresponding X is:'/(2X,5D15.6))") &
            nf, f, (x(i), i=1, N)
    end if
    if (nf <= Npt) goto 70
    if (knew == -1) goto 530
    !
    !     Use the quadratic model to predict the change in F due to the step D,
    !     and set DIFF to the error of this prediction.
    !
    vquad = zero
    ih = 0
    do j = 1, N
        vquad = vquad + d (j)*gq (j)
        do i = 1, j
            ih = ih + 1
            temp = d (i)*xnew (j) + d (j)*xopt (i)
            if (i == j) temp = half*temp
            vquad = vquad + temp*hq (ih)
        end do
    end do
    do k = 1, Npt
        vquad = vquad + pq (k)*w (k)
    end do
    diff = f - fopt - vquad
    diffc = diffb
    diffb = diffa
    diffa = abs (diff)
    if (dnorm > rho) nfsav = nf
    !
    !     Update FOPT and XOPT if the new F is the least value of the objective
    !     function so far. The branch when KNEW is positive occurs if D is not
    !     a trust region step.
    !
    fsave = fopt
    if (f < fopt) then
        fopt = f
        xoptsq = zero
        do i = 1, N
            xopt (i) = xnew (i)
            xoptsq = xoptsq + xopt (i)**2
        end do
    end if
    ksave = knew
    if (knew > 0) goto 410
    !
    !     Pick the next value of DELTA after a trust region step.
    !
    if (vquad >= zero) then
        if (Iprint > 0) write(*, *) "Return from NEWUOA because a trust region step has failed to reduce Q."
        newuob = .false.
        goto 530
    end if
    ratio = (f - fsave)/vquad
    if (ratio <= tenth) then
        delta = half*dnorm
    else if (ratio <= 0.7_dp) then
        delta = max (half*delta, dnorm)
    else
        delta = max (half*delta, dnorm + dnorm)
    end if
    if (delta <= 1.5_dp*rho) delta = rho
    !
    !     Set KNEW to the index of the next interpolation point to be deleted.
    !
    rhosq = max (tenth*delta, rho)**2
    ktemp = 0
    detrat = zero
    if (f >= fsave) then
        ktemp = kopt
        detrat = one
    end if
    do k = 1, Npt
        hdiag = zero
        do j = 1, nptm
            temp = one
            if (j < idz) temp = -one
            hdiag = hdiag + temp*zmat (k, j)**2
        end do
        temp = abs (beta*hdiag + vlag(k)**2)
        distsq = zero
        do j = 1, N
            distsq = distsq + (xpt(k, j) - xopt(j))**2
        end do
        if (distsq > rhosq) temp = temp*(distsq/rhosq)**3
        if (temp > detrat .and. k /= ktemp) then
            detrat = temp
            knew = k
        end if
    end do
    if (knew == 0) goto 460
    !
    !     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
    !     can be moved. Begin the updating of the quadratic model, starting
    !     with the explicit second derivative term.
    !
410 call update_newuoa (N, Npt, bmat, zmat, idz, Ndim, vlag, beta, knew, w)
    fval (knew) = f
    ih = 0
    do i = 1, N
        temp = pq (knew)*xpt (knew, i)
        do j = 1, i
            ih = ih + 1
            hq (ih) = hq (ih) + temp*xpt (knew, j)
        end do
    end do
    pq (knew) = zero
    !
    !     Update the other second derivative parameters, and then the gradient
    !     vector of the model. Also include the new interpolation point.
    !
    do j = 1, nptm
        temp = diff*zmat (knew, j)
        if (j < idz) temp = -temp
        do k = 1, Npt
            pq (k) = pq (k) + temp*zmat (k, j)
        end do
    end do
    gqsq = zero
    do i = 1, N
        gq (i) = gq (i) + diff*bmat (knew, i)
        gqsq = gqsq + gq (i)**2
        xpt (knew, i) = xnew (i)
    end do
    !
    !     If a trust region step makes a small change to the objective function,
    !     then calculate the gradient of the least Frobenius norm interpolant at
    !     XBASE, and store it in W, using VLAG for a vector of right hand sides.
    !
    if (ksave == 0 .and. delta == rho) then
        if (abs(ratio) > 1.0e-2_dp) then
            itest = 0
        else
            do k = 1, Npt
                vlag (k) = fval (k) - fval (kopt)
            end do
            gisq = zero
            do i = 1, N
                sum = zero
                do k = 1, Npt
                    sum = sum + bmat (k, i)*vlag (k)
                end do
                gisq = gisq + sum*sum
                w (i) = sum
            end do
            !
            !     Test whether to replace the new quadratic model by the least Frobenius
            !     norm interpolant, making the replacement if the test is satisfied.
            !
            itest = itest + 1
            if (gqsq < 100.0_dp*gisq) itest = 0
            if (itest >= 3) then
                do i = 1, N
                    gq (i) = w (i)
                end do
                do ih = 1, nh
                    hq (ih) = zero
                end do
                do j = 1, nptm
                    w (j) = zero
                    do k = 1, Npt
                        w (j) = w (j) + vlag (k)*zmat (k, j)
                    end do
                    if (j < idz) w (j) = -w (j)
                end do
                do k = 1, Npt
                    pq (k) = zero
                    do j = 1, nptm
                        pq (k) = pq (k) + zmat (k, j)*w (j)
                    end do
                end do
                itest = 0
            end if
        end if
    end if
    if (f < fsave) kopt = knew
    !
    !     If a trust region step has provided a sufficient decrease in F, then
    !     branch for another trust region calculation. The case KSAVE>0 occurs
    !     when the new function value was calculated by a model step.
    !
    if (f <= fsave + tenth*vquad) goto 100
    if (ksave > 0) goto 100
    !
    !     Alternatively, find out if the interpolation points are close enough
    !     to the best point so far.
    !
    knew = 0
460 distsq = 4.0_dp*delta*delta
    do k = 1, Npt
        sum = zero
        do j = 1, N
            sum = sum + (xpt(k, j) - xopt(j))**2
        end do
        if (sum > distsq) then
            knew = k
            distsq = sum
        end if
    end do
    !
    !     If KNEW is positive, then set DSTEP, and branch back for the next
    !     iteration, which will generate a "model step".
    !
    if (knew > 0) then
        dstep = max (min(tenth*sqrt(distsq), half*delta), rho)
        dsq = dstep*dstep
        goto 120
    end if
    if (ratio > zero) goto 100
    if (max(delta, dnorm) > rho) goto 100
    !
    !     The calculations with the current value of RHO are complete. Pick the
    !     next values of RHO and DELTA.
    !
490 if (rho > Rhoend) then
        delta = half*rho
        ratio = rho/Rhoend
        if (ratio <= 16.0_dp) then
            rho = Rhoend
        else if (ratio <= 250.0_dp) then
            rho = sqrt (ratio)*Rhoend
        else
            rho = tenth*rho
        end if
        delta = max (delta, rho)
        this%Last_bestfit = fopt
        if (Iprint >= 2) then
            if (Iprint >= 3) print "(5X)"
            write(*, "(/4X,'New RHO =',1PD11.4,5X,'Number of function values =',I6)") rho, nf
            write(*, "(4X,'Least value of F =',1PD23.15,9X,'The corresponding X is:'/(2X,5D15.6))") &
                fopt, (xbase(i) + xopt(i), i=1, N)
        end if
        goto 90
    end if
    !
    !     Return from the calculation, after another Newton-Raphson step, if
    !     it is too short to have been tried before.
    !
    if (knew == -1) goto 290
530 if (fopt <= f) then
        do i = 1, N
            x (i) = xbase (i) + xopt (i)
        end do
        f = fopt
    end if
    if (Iprint >= 1) then
        write(*, "(/4X,'At the return from NEWUOA',5X,'Number of function values =',I6)") nf
        write(*, "(4X,'Least value of F =',1PD23.15,9X,'The corresponding X is:'/(2X,5D15.6))") &
            f, (x(i), i=1, N)
    end if

    end function newuob

    subroutine bigden (N, Npt, xopt, xpt, bmat, zmat, Idz, Ndim, Kopt, Knew, d, w, vlag, &
        Beta, s, wvec, prod)
    real(dp) alpha, angle, Beta, dd, denmax, denold, densav, diff, ds, dstemp, &
        dtest, half, one, quart
    real(dp) ss, ssden, sstemp, step, sum, sumold, tau, &
        temp, tempa, tempb, tempc, two, twopi, xoptd, xopts
    real(dp) xoptsq, zero
    integer i, Idz, ip, isave, iterc, iu, j, jc, k, Knew, &
        Kopt, ksav, N, Ndim, Npt, nptm, nw
    real(dp) xopt (*), xpt (Npt, *), bmat (Ndim, *), zmat (Npt, *), d (*), w (*), vlag &
        (*), s (*), wvec (Ndim, *), prod (Ndim, *)
    real(dp) den (9), denex (9), par (9)
    !
    !     N is the number of variables.
    !     NPT is the number of interpolation equations.
    !     XOPT is the best interpolation point so far.
    !     XPT contains the coordinates of the current interpolation points.
    !     BMAT provides the last N columns of H.
    !     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
    !     NDIM is the first dimension of BMAT and has the value NPT+N.
    !     KOPT is the index of the optimal interpolation point.
    !     KNEW is the index of the interpolation point that is going to be moved.
    !     D will be set to the step from XOPT to the new point, and on entry it
    !       should be the D that was calculated by the last call of BIGLAG. The
    !       length of the initial D provides a trust region bound on the final D.
    !     W will be set to Wcheck for the final choice of D.
    !     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
    !     BETA will be set to the value that will occur in the updating formula
    !       when the KNEW-th interpolation point is moved to its new position.
    !     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
    !       for working space.
    !
    !     D is calculated in a way that should provide a denominator with a large
    !     modulus in the updating formula when the KNEW-th interpolation point is
    !     shifted to the new position XOPT+D.
    !
    !     Set some constants.
    !
    half = 0.5_dp
    one = 1.0_dp
    quart = 0.25_dp
    two = 2.0_dp
    zero = 0.0_dp
    twopi = 8.0_dp*atan (one)
    nptm = Npt - N - 1
    !
    !     Store the first NPT elements of the KNEW-th column of H in W(N+1)
    !     to W(N+NPT).
    !
    do k = 1, Npt
        w (N + k) = zero
    end do
    do j = 1, nptm
        temp = zmat (Knew, j)
        if (j < Idz) temp = -temp
        do k = 1, Npt
            w (N + k) = w (N + k) + temp*zmat (k, j)
        end do
    end do
    alpha = w (N + Knew)
    !
    !     The initial search direction D is taken from the last call of BIGLAG,
    !     and the initial S is set below, usually to the direction from X_OPT
    !     to X_KNEW, but a different direction to an interpolation point may
    !     be chosen, in order to prevent S from being nearly parallel to D.
    !
    dd = zero
    ds = zero
    ss = zero
    xoptsq = zero
    do i = 1, N
        dd = dd + d (i)**2
        s (i) = xpt (Knew, i) - xopt (i)
        ds = ds + d (i)*s (i)
        ss = ss + s (i)**2
        xoptsq = xoptsq + xopt (i)**2
    end do
    if (ds*ds > 0.99_dp*dd*ss) then
        ksav = Knew
        dtest = ds*ds/ss
        do k = 1, Npt
            if (k /= Kopt) then
                dstemp = zero
                sstemp = zero
                do i = 1, N
                    diff = xpt (k, i) - xopt (i)
                    dstemp = dstemp + d (i)*diff
                    sstemp = sstemp + diff*diff
                end do
                if (dstemp*dstemp/sstemp < dtest) then
                    ksav = k
                    dtest = dstemp*dstemp/sstemp
                    ds = dstemp
                    ss = sstemp
                end if
            end if
        end do
        do i = 1, N
            s (i) = xpt (ksav, i) - xopt (i)
        end do
    end if
    ssden = dd*ss - ds*ds
    iterc = 0
    densav = zero
    !
    !     Begin the iteration by overwriting S with a vector that has the
    !     required length and direction.
    !
70  iterc = iterc + 1
    temp = one/sqrt (ssden)
    xoptd = zero
    xopts = zero
    do i = 1, N
        s (i) = temp*(dd*s(i) - ds*d(i))
        xoptd = xoptd + xopt (i)*d (i)
        xopts = xopts + xopt (i)*s (i)
    end do
    !
    !     Set the coefficients of the first two terms of BETA.
    !
    tempa = half*xoptd*xoptd
    tempb = half*xopts*xopts
    den (1) = dd*(xoptsq + half*dd) + tempa + tempb
    den (2) = two*xoptd*dd
    den (3) = two*xopts*dd
    den (4) = tempa - tempb
    den (5) = xoptd*xopts
    do i = 6, 9
        den (i) = zero
    end do
    !
    !     Put the coefficients of Wcheck in WVEC.
    !
    do k = 1, Npt
        tempa = zero
        tempb = zero
        tempc = zero
        do i = 1, N
            tempa = tempa + xpt (k, i)*d (i)
            tempb = tempb + xpt (k, i)*s (i)
            tempc = tempc + xpt (k, i)*xopt (i)
        end do
        wvec (k, 1) = quart*(tempa*tempa + tempb*tempb)
        wvec (k, 2) = tempa*tempc
        wvec (k, 3) = tempb*tempc
        wvec (k, 4) = quart*(tempa*tempa - tempb*tempb)
        wvec (k, 5) = half*tempa*tempb
    end do
    do i = 1, N
        ip = i + Npt
        wvec (ip, 1) = zero
        wvec (ip, 2) = d (i)
        wvec (ip, 3) = s (i)
        wvec (ip, 4) = zero
        wvec (ip, 5) = zero
    end do
    !
    !     Put the coefficients of THETA*Wcheck in PROD.
    !
    do jc = 1, 5
        nw = Npt
        if (jc == 2 .or. jc == 3) nw = Ndim
        do k = 1, Npt
            prod (k, jc) = zero
        end do
        do j = 1, nptm
            sum = zero
            do k = 1, Npt
                sum = sum + zmat (k, j)*wvec (k, jc)
            end do
            if (j < Idz) sum = -sum
            do k = 1, Npt
                prod (k, jc) = prod (k, jc) + sum*zmat (k, j)
            end do
        end do
        if (nw == Ndim) then
            do k = 1, Npt
                sum = zero
                do j = 1, N
                    sum = sum + bmat (k, j)*wvec (Npt + j, jc)
                end do
                prod (k, jc) = prod (k, jc) + sum
            end do
        end if
        do j = 1, N
            sum = zero
            do i = 1, nw
                sum = sum + bmat (i, j)*wvec (i, jc)
            end do
            prod (Npt + j, jc) = sum
        end do
    end do
    !
    !     Include in DEN the part of BETA that depends on THETA.
    !
    do k = 1, Ndim
        sum = zero
        do i = 1, 5
            par (i) = half*prod (k, i)*wvec (k, i)
            sum = sum + par (i)
        end do
        den (1) = den (1) - par (1) - sum
        tempa = prod (k, 1)*wvec (k, 2) + prod (k, 2)*wvec (k, 1)
        tempb = prod (k, 2)*wvec (k, 4) + prod (k, 4)*wvec (k, 2)
        tempc = prod (k, 3)*wvec (k, 5) + prod (k, 5)*wvec (k, 3)
        den (2) = den (2) - tempa - half*(tempb + tempc)
        den (6) = den (6) - half*(tempb - tempc)
        tempa = prod (k, 1)*wvec (k, 3) + prod (k, 3)*wvec (k, 1)
        tempb = prod (k, 2)*wvec (k, 5) + prod (k, 5)*wvec (k, 2)
        tempc = prod (k, 3)*wvec (k, 4) + prod (k, 4)*wvec (k, 3)
        den (3) = den (3) - tempa - half*(tempb - tempc)
        den (7) = den (7) - half*(tempb + tempc)
        tempa = prod (k, 1)*wvec (k, 4) + prod (k, 4)*wvec (k, 1)
        den (4) = den (4) - tempa - par (2) + par (3)
        tempa = prod (k, 1)*wvec (k, 5) + prod (k, 5)*wvec (k, 1)
        tempb = prod (k, 2)*wvec (k, 3) + prod (k, 3)*wvec (k, 2)
        den (5) = den (5) - tempa - half*tempb
        den (8) = den (8) - par (4) + par (5)
        tempa = prod (k, 4)*wvec (k, 5) + prod (k, 5)*wvec (k, 4)
        den (9) = den (9) - half*tempa
    end do
    !
    !     Extend DEN so that it holds all the coefficients of DENOM.
    !
    sum = zero
    do i = 1, 5
        par (i) = half*prod (Knew, i)**2
        sum = sum + par (i)
    end do
    denex (1) = alpha*den (1) + par (1) + sum
    tempa = two*prod (Knew, 1)*prod (Knew, 2)
    tempb = prod (Knew, 2)*prod (Knew, 4)
    tempc = prod (Knew, 3)*prod (Knew, 5)
    denex (2) = alpha*den (2) + tempa + tempb + tempc
    denex (6) = alpha*den (6) + tempb - tempc
    tempa = two*prod (Knew, 1)*prod (Knew, 3)
    tempb = prod (Knew, 2)*prod (Knew, 5)
    tempc = prod (Knew, 3)*prod (Knew, 4)
    denex (3) = alpha*den (3) + tempa + tempb - tempc
    denex (7) = alpha*den (7) + tempb + tempc
    tempa = two*prod (Knew, 1)*prod (Knew, 4)
    denex (4) = alpha*den (4) + tempa + par (2) - par (3)
    tempa = two*prod (Knew, 1)*prod (Knew, 5)
    denex (5) = alpha*den (5) + tempa + prod (Knew, 2)*prod (Knew, 3)
    denex (8) = alpha*den (8) + par (4) - par (5)
    denex (9) = alpha*den (9) + prod (Knew, 4)*prod (Knew, 5)
    !
    !     Seek the value of the angle that maximizes the modulus of DENOM.
    !
    sum = denex (1) + denex (2) + denex (4) + denex (6) + denex (8)
    denold = sum
    denmax = sum
    isave = 0
    iu = 49
    temp = twopi/real (iu + 1, dp)
    par (1) = one
    do i = 1, iu
        angle = real (i, dp)*temp
        par (2) = cos (angle)
        par (3) = sin (angle)
        do j = 4, 8, 2
            par (j) = par (2)*par (j - 2) - par (3)*par (j - 1)
            par (j + 1) = par (2)*par (j - 1) + par (3)*par (j - 2)
        end do
        sumold = sum
        sum = zero
        do j = 1, 9
            sum = sum + denex (j)*par (j)
        end do
        if (abs(sum) > abs(denmax)) then
            denmax = sum
            isave = i
            tempa = sumold
        else if (i == isave + 1) then
            tempb = sum
        end if
    end do
    if (isave == 0) tempa = sum
    if (isave == iu) tempb = denold
    step = zero
    if (tempa /= tempb) then
        tempa = tempa - denmax
        tempb = tempb - denmax
        step = half*(tempa - tempb)/(tempa + tempb)
    end if
    angle = temp*(real(isave, dp) + step)
    !
    !     Calculate the new parameters of the denominator, the new VLAG vector
    !     and the new D. Then test for convergence.
    !
    par (2) = cos (angle)
    par (3) = sin (angle)
    do j = 4, 8, 2
        par (j) = par (2)*par (j - 2) - par (3)*par (j - 1)
        par (j + 1) = par (2)*par (j - 1) + par (3)*par (j - 2)
    end do
    Beta = zero
    denmax = zero
    do j = 1, 9
        Beta = Beta + den (j)*par (j)
        denmax = denmax + denex (j)*par (j)
    end do
    do k = 1, Ndim
        vlag (k) = zero
        do j = 1, 5
            vlag (k) = vlag (k) + prod (k, j)*par (j)
        end do
    end do
    tau = vlag (Knew)
    dd = zero
    tempa = zero
    tempb = zero
    do i = 1, N
        d (i) = par (2)*d (i) + par (3)*s (i)
        w (i) = xopt (i) + d (i)
        dd = dd + d (i)**2
        tempa = tempa + d (i)*w (i)
        tempb = tempb + w (i)*w (i)
    end do
    if (iterc >= N) goto 340
    if (iterc > 1) densav = max (densav, denold)
    if (abs(denmax) <= 1.1_dp*abs(densav)) goto 340
    densav = denmax
    !
    !     Set S to half the gradient of the denominator with respect to D.
    !     Then branch for the next iteration.
    !
    do i = 1, N
        temp = tempa*xopt (i) + tempb*d (i) - vlag (Npt + i)
        s (i) = tau*bmat (Knew, i) + alpha*temp
    end do
    do k = 1, Npt
        sum = zero
        do j = 1, N
            sum = sum + xpt (k, j)*w (j)
        end do
        temp = (tau*w(N + k) - alpha*vlag(k))*sum
        do i = 1, N
            s (i) = s (i) + temp*xpt (k, i)
        end do
    end do
    ss = zero
    ds = zero
    do i = 1, N
        ss = ss + s (i)**2
        ds = ds + d (i)*s (i)
    end do
    ssden = dd*ss - ds*ds
    if (ssden >= 1.0e-8_dp*dd*ss) goto 70
    !
    !     Set the vector W before the RETURN from the subroutine.
    !
340 do k = 1, Ndim
        w (k) = zero
        do j = 1, 5
            w (k) = w (k) + wvec (k, j)*par (j)
        end do
    end do
    vlag (Kopt) = vlag (Kopt) + one

    end subroutine bigden

    subroutine biglag (N, Npt, xopt, xpt, bmat, zmat, Idz, Ndim, Knew, Delta, d, Alpha, &
        hcol, gc, gd, s, w)

    real(dp) Alpha, angle, cf1, cf2, cf3, cf4, cf5, cth, &
        dd, delsq, Delta, denom, dhd, gg, half
    integer i, Idz, isave, iterc, iu, j, k, Knew, N, Ndim, &
        Npt, nptm
    real(dp) one, scale, sp, ss, step, sth, sum, tau, &
        taubeg, taumax, tauold, temp, tempa, tempb, twopi, zero
    real(dp) xopt (*), xpt (Npt, *), bmat (Ndim, *), zmat (Npt, *), d (*), hcol (*), &
        gc(*), gd (*), s (*), w (*)
    !
    !     N is the number of variables.
    !     NPT is the number of interpolation equations.
    !     XOPT is the best interpolation point so far.
    !     XPT contains the coordinates of the current interpolation points.
    !     BMAT provides the last N columns of H.
    !     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
    !     NDIM is the first dimension of BMAT and has the value NPT+N.
    !     KNEW is the index of the interpolation point that is going to be moved.
    !     DELTA is the current trust region bound.
    !     D will be set to the step from XOPT to the new point.
    !     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
    !     HCOL, GC, GD, S and W will be used for working space.
    !
    !     The step D is calculated in a way that attempts to maximize the modulus
    !     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
    !     the KNEW-th Lagrange function.
    !
    !     Set some constants.
    !
    half = 0.5_dp
    one = 1.0_dp
    zero = 0.0_dp
    twopi = 8.0_dp*atan (one)
    delsq = Delta*Delta
    nptm = Npt - N - 1
    !
    !     Set the first NPT components of HCOL to the leading elements of the
    !     KNEW-th column of H.
    !
    iterc = 0
    do k = 1, Npt
        hcol (k) = zero
    end do
    do j = 1, nptm
        temp = zmat (Knew, j)
        if (j < Idz) temp = -temp
        do k = 1, Npt
            hcol (k) = hcol (k) + temp*zmat (k, j)
        end do
    end do
    Alpha = hcol (Knew)
    !
    !     Set the unscaled initial direction D. Form the gradient of LFUNC at
    !     XOPT, and multiply D by the second derivative matrix of LFUNC.
    !
    dd = zero
    do i = 1, N
        d (i) = xpt (Knew, i) - xopt (i)
        gc (i) = bmat (Knew, i)
        gd (i) = zero
        dd = dd + d (i)**2
    end do
    do k = 1, Npt
        temp = zero
        sum = zero
        do j = 1, N
            temp = temp + xpt (k, j)*xopt (j)
            sum = sum + xpt (k, j)*d (j)
        end do
        temp = hcol (k)*temp
        sum = hcol (k)*sum
        do i = 1, N
            gc (i) = gc (i) + temp*xpt (k, i)
            gd (i) = gd (i) + sum*xpt (k, i)
        end do
    end do
    !
    !     Scale D and GD, with a sign change if required. Set S to another
    !     vector in the initial two dimensional subspace.
    !
    gg = zero
    sp = zero
    dhd = zero
    do i = 1, N
        gg = gg + gc (i)**2
        sp = sp + d (i)*gc (i)
        dhd = dhd + d (i)*gd (i)
    end do
    scale = Delta/sqrt (dd)
    if (sp*dhd < zero) scale = -scale
    temp = zero
    if (sp*sp > 0.99_dp*dd*gg) temp = one
    tau = scale*(abs(sp) + half*scale*abs(dhd))
    if (gg*delsq < 0.01_dp*tau*tau) temp = one
    do i = 1, N
        d (i) = scale*d (i)
        gd (i) = scale*gd (i)
        s (i) = gc (i) + temp*gd (i)
    end do
    !
    !     Begin the iteration by overwriting S with a vector that has the
    !     required length and direction, except that termination occurs if
    !     the given D and S are nearly parallel.
    !
80  iterc = iterc + 1
    dd = zero
    sp = zero
    ss = zero
    do i = 1, N
        dd = dd + d (i)**2
        sp = sp + d (i)*s (i)
        ss = ss + s (i)**2
    end do
    temp = dd*ss - sp*sp
    if (temp <= 1.0e-8_dp*dd*ss) goto 160
    denom = sqrt (temp)
    do i = 1, N
        s (i) = (dd*s(i) - sp*d(i))/denom
        w (i) = zero
    end do
    !
    !     Calculate the coefficients of the objective function on the circle,
    !     beginning with the multiplication of S by the second derivative matrix.
    !
    do k = 1, Npt
        sum = zero
        do j = 1, N
            sum = sum + xpt (k, j)*s (j)
        end do
        sum = hcol (k)*sum
        do i = 1, N
            w (i) = w (i) + sum*xpt (k, i)
        end do
    end do
    cf1 = zero
    cf2 = zero
    cf3 = zero
    cf4 = zero
    cf5 = zero
    do i = 1, N
        cf1 = cf1 + s (i)*w (i)
        cf2 = cf2 + d (i)*gc (i)
        cf3 = cf3 + s (i)*gc (i)
        cf4 = cf4 + d (i)*gd (i)
        cf5 = cf5 + s (i)*gd (i)
    end do
    cf1 = half*cf1
    cf4 = half*cf4 - cf1
    !
    !     Seek the value of the angle that maximizes the modulus of TAU.
    !
    taubeg = cf1 + cf2 + cf4
    taumax = taubeg
    tauold = taubeg
    isave = 0
    iu = 49
    temp = twopi/real (iu + 1, dp)
    do i = 1, iu
        angle = real (i, dp)*temp
        cth = cos (angle)
        sth = sin (angle)
        tau = cf1 + (cf2 + cf4*cth)*cth + (cf3 + cf5*cth)*sth
        if (abs(tau) > abs(taumax)) then
            taumax = tau
            isave = i
            tempa = tauold
        else if (i == isave + 1) then
            tempb = tau
        end if
        tauold = tau
    end do
    if (isave == 0) tempa = tau
    if (isave == iu) tempb = taubeg
    step = zero
    if (tempa /= tempb) then
        tempa = tempa - taumax
        tempb = tempb - taumax
        step = half*(tempa - tempb)/(tempa + tempb)
    end if
    angle = temp*(real(isave, dp) + step)
    !
    !     Calculate the new D and GD. Then test for convergence.
    !
    cth = cos (angle)
    sth = sin (angle)
    tau = cf1 + (cf2 + cf4*cth)*cth + (cf3 + cf5*cth)*sth
    do i = 1, N
        d (i) = cth*d (i) + sth*s (i)
        gd (i) = cth*gd (i) + sth*w (i)
        s (i) = gc (i) + gd (i)
    end do
    if (abs(tau) <= 1.1_dp*abs(taubeg)) goto 160
    if (iterc < N) goto 80
160 return

    end subroutine biglag

    subroutine trsapp (N, Npt, xopt, xpt, gq, hq, pq, Delta, step, d, g, hd, hs, Crvmin)
    real(dp) alpha, angle, angtest, bstep, cf, Crvmin, cth, &
        dd, delsq, Delta, dg, dhd, dhs, ds, gg, ggbeg, &
        ggsav
    real(dp) half, qadd, qbeg, qmin, qnew, qred, qsav, ratio, reduc, &
        sg, sgk, shs, ss, sth, temp, tempa, tempb, twopi, zero
    integer i, ih, isave, iterc, itermax, itersw, iu, j, k, &
        N, Npt
    real(dp) xopt (*), xpt (Npt, *), gq (*), hq (*), pq (*), step (*), d (*), g (*), &
        hd (*), hs (*)
    !
    !     N is the number of variables of a quadratic objective function, Q say.
    !     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
    !       in order to define the current quadratic model Q.
    !     DELTA is the trust region radius, and has to be positive.
    !     STEP will be set to the calculated trial step.
    !     The arrays D, G, HD and HS will be used for working space.
    !     CRVMIN will be set to the least curvature of H along the conjugate
    !       directions that occur, except that it is set to zero if STEP goes
    !       all the way to the trust region boundary.
    !
    !     The calculation of STEP begins with the truncated conjugate gradient
    !     method. If the boundary of the trust region is reached, then further
    !     changes to STEP may be made, each one being in the 2D space spanned
    !     by the current STEP and the corresponding gradient of Q. Thus STEP
    !     should provide a substantial reduction to Q within the trust region.
    !
    !     Initialization, which includes setting HD to H times XOPT.
    !
    half = 0.5_dp
    zero = 0.0_dp
    twopi = 8.0_dp*atan (1.0_dp)
    delsq = Delta*Delta
    iterc = 0
    itermax = N
    itersw = itermax
    do i = 1, N
        d (i) = xopt (i)
    end do
    goto 170
    !
    !     Prepare for the first line search.
    !
20  qred = zero
    dd = zero
    do i = 1, N
        step (i) = zero
        hs (i) = zero
        g (i) = gq (i) + hd (i)
        d (i) = -g (i)
        dd = dd + d (i)**2
    end do
    Crvmin = zero
    if (dd == zero) goto 160
    ds = zero
    ss = zero
    gg = dd
    ggbeg = gg
    !
    !     Calculate the step to the trust region boundary and the product HD.
    !
40  iterc = iterc + 1
    temp = delsq - ss
    bstep = temp/(ds + sqrt(ds*ds + dd*temp))
    goto 170
50  dhd = zero
    do j = 1, N
        dhd = dhd + d (j)*hd (j)
    end do
    !
    !     Update CRVMIN and set the step-length ALPHA.
    !
    alpha = bstep
    if (dhd > zero) then
        temp = dhd/dd
        if (iterc == 1) Crvmin = temp
        Crvmin = min (Crvmin, temp)
        alpha = min (alpha, gg/dhd)
    end if
    qadd = alpha*(gg - half*alpha*dhd)
    qred = qred + qadd
    !
    !     Update STEP and HS.
    !
    ggsav = gg
    gg = zero
    do i = 1, N
        step (i) = step (i) + alpha*d (i)
        hs (i) = hs (i) + alpha*hd (i)
        gg = gg + (g(i) + hs(i))**2
    end do
    !
    !     Begin another conjugate direction iteration if required.
    !
    if (alpha < bstep) then
        if (qadd <= 0.01_dp*qred) goto 160
        if (gg <= 1.0e-4_dp*ggbeg) goto 160
        if (iterc == itermax) goto 160
        temp = gg/ggsav
        dd = zero
        ds = zero
        ss = zero
        do i = 1, N
            d (i) = temp*d (i) - g (i) - hs (i)
            dd = dd + d (i)**2
            ds = ds + d (i)*step (i)
            ss = ss + step (i)**2
        end do
        if (ds <= zero) goto 160
        if (ss < delsq) goto 40
    end if
    Crvmin = zero
    itersw = iterc
    !
    !     Test whether an alternative iteration is required.
    !
90  if (gg <= 1.0e-4_dp*ggbeg) goto 160
    sg = zero
    shs = zero
    do i = 1, N
        sg = sg + step (i)*g (i)
        shs = shs + step (i)*hs (i)
    end do
    sgk = sg + shs
    angtest = sgk/sqrt (gg*delsq)
    if (angtest <= -0.99_dp) goto 160
    !
    !     Begin the alternative iteration by calculating D and HD and some
    !     scalar products.
    !
    iterc = iterc + 1
    temp = sqrt (delsq*gg - sgk*sgk)
    tempa = delsq/temp
    tempb = sgk/temp
    do i = 1, N
        d (i) = tempa*(g(i) + hs(i)) - tempb*step (i)
    end do
    goto 170
120 dg = zero
    dhd = zero
    dhs = zero
    do i = 1, N
        dg = dg + d (i)*g (i)
        dhd = dhd + hd (i)*d (i)
        dhs = dhs + hd (i)*step (i)
    end do
    !
    !     Seek the value of the angle that minimizes Q.
    !
    cf = half*(shs - dhd)
    qbeg = sg + cf
    qsav = qbeg
    qmin = qbeg
    isave = 0
    iu = 49
    temp = twopi/real (iu + 1, dp)
    do i = 1, iu
        angle = real (i, dp)*temp
        cth = cos (angle)
        sth = sin (angle)
        qnew = (sg + cf*cth)*cth + (dg + dhs*cth)*sth
        if (qnew < qmin) then
            qmin = qnew
            isave = i
            tempa = qsav
        else if (i == isave + 1) then
            tempb = qnew
        end if
        qsav = qnew
    end do
    if (isave == zero) tempa = qnew
    if (isave == iu) tempb = qbeg
    angle = zero
    if (tempa /= tempb) then
        tempa = tempa - qmin
        tempb = tempb - qmin
        angle = half*(tempa - tempb)/(tempa + tempb)
    end if
    angle = temp*(real(isave, dp) + angle)
    !
    !     Calculate the new STEP and HS. Then test for convergence.
    !
    cth = cos (angle)
    sth = sin (angle)
    reduc = qbeg - (sg + cf*cth)*cth - (dg + dhs*cth)*sth
    gg = zero
    do i = 1, N
        step (i) = cth*step (i) + sth*d (i)
        hs (i) = cth*hs (i) + sth*hd (i)
        gg = gg + (g(i) + hs(i))**2
    end do
    qred = qred + reduc
    ratio = reduc/qred
    if (iterc < itermax .and. ratio > 0.01_dp) goto 90
160 return
    !
    !     The following instructions act as a subroutine for setting the vector
    !     HD to the vector D multiplied by the second derivative matrix of Q.
    !     They are called from three different places, which are distinguished
    !     by the value of ITERC.
    !
170 do i = 1, N
        hd (i) = zero
    end do
    do k = 1, Npt
        temp = zero
        do j = 1, N
            temp = temp + xpt (k, j)*d (j)
        end do
        temp = temp*pq (k)
        do i = 1, N
            hd (i) = hd (i) + temp*xpt (k, i)
        end do
    end do
    ih = 0
    do j = 1, N
        do i = 1, j
            ih = ih + 1
            if (i < j) hd (j) = hd (j) + hq (ih)*d (i)
            hd (i) = hd (i) + hq (ih)*d (j)
        end do
    end do
    if (iterc == 0) goto 20
    if (iterc <= itersw) goto 50
    goto 120

    end subroutine trsapp

    subroutine update_newuoa (N, Npt, bmat, zmat, Idz, Ndim, vlag, Beta, Knew, w)
    real(dp) alpha, Beta, denom, one, scala, scalb, tau, &
        tausq, temp, tempa, tempb, zero
    integer i, Idz, iflag, j, ja, jb, jl, jp, Knew, N, &
        Ndim, Npt, nptm
    real(dp) bmat (Ndim, *), zmat (Npt, *), vlag (*), w (*)
    !
    !     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
    !     interpolation point that has index KNEW. On entry, VLAG contains the
    !     components of the vector Theta*Wcheck+e_b of the updating formula
    !     (6.11), and BETA holds the value of the parameter that has this name.
    !     The vector W is used for working space.
    !
    !     Set some constants.
    !
    one = 1.0_dp
    zero = 0.0_dp
    nptm = Npt - N - 1
    !
    !     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    !
    jl = 1
    do j = 2, nptm
        if (j == Idz) then
            jl = Idz
        else if (zmat(Knew, j) /= zero) then
            temp = sqrt (zmat(Knew, jl)**2 + zmat(Knew, j)**2)
            tempa = zmat (Knew, jl)/temp
            tempb = zmat (Knew, j)/temp
            do i = 1, Npt
                temp = tempa*zmat (i, jl) + tempb*zmat (i, j)
                zmat (i, j) = tempa*zmat (i, j) - tempb*zmat (i, jl)
                zmat (i, jl) = temp
            end do
            zmat (Knew, j) = zero
        end if
    end do
    !
    !     Put the first NPT components of the KNEW-th column of HLAG into W,
    !     and calculate the parameters of the updating formula.
    !
    tempa = zmat (Knew, 1)
    if (Idz >= 2) tempa = -tempa
    if (jl > 1) tempb = zmat (Knew, jl)
    do i = 1, Npt
        w (i) = tempa*zmat (i, 1)
        if (jl > 1) w (i) = w (i) + tempb*zmat (i, jl)
    end do
    alpha = w (Knew)
    tau = vlag (Knew)
    tausq = tau*tau
    denom = alpha*Beta + tausq
    vlag (Knew) = vlag (Knew) - one
    !
    !     Complete the updating of ZMAT when there is only one nonzero element
    !     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
    !     then the first column of ZMAT will be exchanged with another one later.
    !
    iflag = 0
    if (jl == 1) then
        temp = sqrt (abs(denom))
        tempb = tempa/temp
        tempa = tau/temp
        do i = 1, Npt
            zmat (i, 1) = tempa*zmat (i, 1) - tempb*vlag (i)
        end do
        if (Idz == 1 .and. temp < zero) Idz = 2
        if (Idz >= 2 .and. temp >= zero) iflag = 1
    else
        !
        !     Complete the updating of ZMAT in the alternative case.
        !
        ja = 1
        if (Beta >= zero) ja = jl
        jb = jl + 1 - ja
        temp = zmat (Knew, jb)/denom
        tempa = temp*Beta
        tempb = temp*tau
        temp = zmat (Knew, ja)
        scala = one/sqrt (abs(Beta)*temp*temp + tausq)
        scalb = scala*sqrt (abs(denom))
        do i = 1, Npt
            zmat (i, ja) = scala*(tau*zmat(i, ja) - temp*vlag(i))
            zmat (i, jb) = scalb*(zmat(i, jb) - tempa*w(i) - tempb*vlag(i))
        end do
        if (denom <= zero) then
            if (Beta < zero) Idz = Idz + 1
            if (Beta >= zero) iflag = 1
        end if
    end if
    !
    !     IDZ is reduced in the following case, and usually the first column
    !     of ZMAT is exchanged with a later one.
    !
    if (iflag == 1) then
        Idz = Idz - 1
        do i = 1, Npt
            temp = zmat (i, 1)
            zmat (i, 1) = zmat (i, Idz)
            zmat (i, Idz) = temp
        end do
    end if
    !
    !     Finally, update the matrix BMAT.
    !
    do j = 1, N
        jp = Npt + j
        w (jp) = bmat (Knew, j)
        tempa = (alpha*vlag(jp) - tau*w(jp))/denom
        tempb = (-Beta*w(jp) - tau*vlag(jp))/denom
        do i = 1, jp
            bmat (i, j) = bmat (i, j) + tempa*vlag (i) + tempb*w (i)
            if (i > Npt) bmat (jp, i - Npt) = bmat (i, j)
        end do
    end do

    end subroutine update_newuoa

    end module Powell
