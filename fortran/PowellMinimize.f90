    !####### constrained and unconstrained minimization ######
    !F90 translated version of plato.asu.edu/ftp/other_software/bobyqa.zip and newuoa
    !Powell 2009 Method http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf

    !AL Sept 2012/2020: translated to F90, changed to input "funkk" class function argument
    !to minimize rather than calfun subroutine as original

    !Typical Usage (typically Npt= 2*dimension+1), to get best fit in params(:) starting from guess

    ! Type(TBOBYQA) :: Minimize !bounded minimization

    ! if (Minimize%BOBYQA(this_class, class_function, dimension, Npt, params,param_min, &
    !           param_max, step_radius, tolerance, DebugLevel, max_eval)) ...
    !

    !Type(TNEWUOA) :: Minimize !unbounded minimization

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

    MODULE Powell
    implicit none
    private

    INTEGER, PARAMETER  :: dp = KIND(1.d0) !SELECTED_REAL_KIND(12, 60)
    INTEGER, PARAMETER  :: func_dp = dp  !SELECTED_REAL_KIND(12, 60)
    INTEGER, PARAMETER :: Powell_CO_prec = dp

    ! the type of function minimized, a class procedure taking in an array of values
    interface
    FUNCTION obj_vec_function(obj, x)
    use precision
    class(*) :: obj
    real(dl) :: obj_vec_function
    real(dl) :: x(:)
    END FUNCTION  obj_vec_function
    end interface

    type TMinimizer
        REAL(dp) :: Last_bestfit
        procedure(obj_vec_function), pointer, nopass :: funkk
        class(*), pointer :: obj
    end type

    type, extends(TMinimizer) :: TBOBYQA
        REAL(dp) :: FVAL_Converge_difference = 0._dp
    contains
    procedure :: BOBYQA
    end type

    type, extends(TMinimizer) :: TNEWUOA
    contains
    procedure :: NEWUOA
    procedure, private :: newuob
    end type

    PUBLIC  :: TBOBYQA , TNEWUOA, Powell_CO_prec
    CONTAINS


    function BOBYQA (this, obj, funk, N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN)
    !Main function for bounded minimization
    use iso_c_binding
    class(TBOBYQA) this
    class(*), target :: obj
    real(dp), external :: funk !a class function f(obj,x) where obj is any class instance
    integer, intent(in)::n, npt, maxfun, iPrint
    real(dp), intent(in) :: RHOBEG, RHOEND
    logical BOBYQA
    INTEGER ibmat , id , ifv , igo , ihq , ipq , isl , isu , &
        & ivl , iw , ixa , ixb , ixn , ixo , ixp , izmat , j , jsl ,&
        & jsu, ndim , np
    REAL(dp) temp ,  zero
    REAL(dp) X(*),XL(*),XU(*)
    real(dp), allocatable :: W(:)
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

    call C_F_PROCPOINTER(c_funloc(funk), this%funkk)
    this%obj => obj

    BOBYQA = .false.
    this%Last_bestfit = 1d30

    NP=N+1
    IF (NPT  <  N+2 .OR. NPT  >  ((N+2)*NP)/2) THEN
        PRINT 10
10      FORMAT (/4X,'Return from BOBYQA because NPT is not in',&
            &      ' the required interval')
        return
    END IF

    allocate(W((NPT+5)*(NPT+N)+3*N*(N+5)/2))

    !  Partition the working space array, so that different parts of it can
    !  be treated separately during the calculation of BOBYQB. The partition
    !  requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
    !  space that is taken by the last array in the argument list of BOBYQB.
    !
    NDIM=NPT+N
    IXB=1
    IXP=IXB+N
    IFV=IXP+N*NPT
    IXO=IFV+NPT
    IGO=IXO+N
    IHQ=IGO+N
    IPQ=IHQ+(N*NP)/2
    IBMAT=IPQ+NPT
    IZMAT=IBMAT+NDIM*N
    ISL=IZMAT+NPT*(NPT-NP)
    ISU=ISL+N
    IXN=ISU+N
    IXA=IXN+N
    ID=IXA+N
    IVL=ID+N
    IW=IVL+NDIM
    !
    !  Return if there is insufficient space between the bounds. Modify the
    !  initial X if necessary in order to avoid conflicts between the bounds
    !  and the construction of the first quadratic model. The lower and upper
    !  bounds on moves from the updated X are set now, in the ISL and ISU
    !  partitions of W, in order to provide useful and exact information about
    !  components of X that become within distance RHOBEG from their bounds.
    !
    ZERO=0.0D0
    DO J=1,N
        TEMP=XU(J)-XL(J)
        IF (TEMP  <  RHOBEG+RHOBEG) THEN
            PRINT 20
20          FORMAT (/4X,'Return from BOBYQA because one of the',&
                &      ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
            return
        END IF
        JSL=ISL+J-1
        JSU=JSL+N
        W(JSL)=XL(J)-X(J)
        W(JSU)=XU(J)-X(J)
        IF (W(JSL)  >=  -RHOBEG) THEN
            IF (W(JSL)  >=  ZERO) THEN
                X(J)=XL(J)
                W(JSL)=ZERO
                W(JSU)=TEMP
            ELSE
                X(J)=XL(J)+RHOBEG
                W(JSL)=-RHOBEG
                W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
            END IF
        ELSE IF (W(JSU)  <=  RHOBEG) THEN
            IF (W(JSU)  <=  ZERO) THEN
                X(J)=XU(J)
                W(JSL)=-TEMP
                W(JSU)=ZERO
            ELSE
                X(J)=XU(J)-RHOBEG
                W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
                W(JSU)=RHOBEG
            END IF
        END IF
    END DO
    !
    !     Make the call of BOBYQB.
    !
    BOBYQA = BOBYQB (this,N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),&
        &  W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),&
        &  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW))

    end function BOBYQA


    function BOBYQB (this, N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,&
        &  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,&
        &  SL,SU,XNEW,XALT,D,VLAG,W)
    class(TBOBYQA) this
    integer NDIM, NPT
    REAL(dp) X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),&
        &  XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),&
        &  SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
    logical BOBYQB
    REAL(dp) adelt , alpha , bdtest , bdtol , beta , biglsq ,    &
        & bsum , cauchy , crvmin , curv , delsq , delta , den ,  &
        & denom , densav , diff , diffa , diffb
    REAL(dp) diffc , dist , distsq , dnorm , dsq , dx , errbig , f , &
        & fopt , fracsq , frhosq , fsave , gisq , &
        & gqsq , half , hdiag , one
    REAL(dp) pqold , ratio , rho , Rhobeg , Rhoend , scaden ,  &
        & sum , suma , sumb , sumpq , sumw , sumz , temp , ten ,&
        & tenth , two
    REAL(dp) vquad , xoptsq , zero
    INTEGER i , ih , ip , Iprint , itest , j , jj , jp , k , kbase ,  &
        & knew , kopt , ksav , Maxfun , N , nf , nfsav , nh , np
    INTEGER nptm , nresc , ntrits


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
    HALF=0.5D0
    ONE=1.0D0
    TEN=10.0D0
    TENTH=0.1D0
    TWO=2.0D0
    ZERO=0.0D0
    NP=N+1
    NPTM=NPT-NP
    NH=(N*NP)/2
    BOBYQB = .true.
    !
    !  The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    !  BMAT and ZMAT for the first iteration, with the corresponding values of
    !  of NF and KOPT, which are the number of calls of CALFUN so far and the
    !  index of the interpolation point at the trust region centre. Then the
    !  initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
    !  less than NPT. GOPT will be updated if KOPT is different from KBASE.
    !
    CALL PRELIM (this,N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,&
        &  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
    XOPTSQ=ZERO
    DO I=1,N
        XOPT(I)=XPT(KOPT,I)
        XOPTSQ=XOPTSQ+XOPT(I)**2
    END DO
    FSAVE=FVAL(1)
    IF (NF  <  NPT) THEN
        IF (IPRINT  >  0) PRINT 390
        GOTO 720
    END IF
    KBASE=1
    !
    !     Complete the settings that are required for the iterative procedure
    !
    RHO=RHOBEG
    DELTA=RHO
    NRESC=NF
    NTRITS=0
    DIFFA=ZERO
    DIFFB=ZERO
    ITEST=0
    NFSAV=NF
    !
    !     Update GOPT if necessary before the first iteration and after each
    !     call of RESCUE that makes a call of CALFUN.
    !
20  IF (KOPT .NE. KBASE) THEN
        IH=0
        DO J=1,N
            DO I=1,J
                IH=IH+1
                IF (I  <  J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
                GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
            END DO
        END DO
        IF (NF  >  NPT) THEN
            DO K=1,NPT
                TEMP=ZERO
                DO J=1,N
                    TEMP=TEMP+XPT(K,J)*XOPT(J)
                END DO
                TEMP=PQ(K)*TEMP
                DO I=1,N
                    GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
                END DO
            END DO
        END IF
    END IF
    !
    !  Generate the next point in the trust region that provides a small value
    !  of the quadratic model subject to the constraints on the variables.
    !  The integer NTRITS is set to the number "trust region" iterations that
    !  have occurred since the last "alternative" iteration. If the length
    !  of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
    !  label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
    !
60  CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,&
        &  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
    DNORM=DMIN1(DELTA,DSQRT(DSQ))
    IF (DNORM  <  HALF*RHO) THEN
        NTRITS=-1
        DISTSQ=(TEN*RHO)**2
        IF (NF  <=  NFSAV+2) GOTO 650
        !
        !  The following choice between labels 650 and 680 depends on whether or
        !  not our work with the current RHO seems to be complete. Either RHO is
        !  decreased or termination occurs if the errors in the quadratic model at
        !  the last three interpolation points compare favourably with predictions
        !  of likely improvements to the model within distance HALF*RHO of XOPT.
        !
        ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
        FRHOSQ=0.125D0*RHO*RHO
        IF (CRVMIN  >  ZERO .AND. ERRBIG  >  FRHOSQ*CRVMIN)&
            &       GOTO 650
        BDTOL=ERRBIG/RHO
        DO J=1,N
            BDTEST=BDTOL
            IF (XNEW(J) .EQ. SL(J)) BDTEST=W(J)
            IF (XNEW(J) .EQ. SU(J)) BDTEST=-W(J)
            IF (BDTEST  <  BDTOL) THEN
                CURV=HQ((J+J*J)/2)
                DO K=1,NPT
                    CURV=CURV+PQ(K)*XPT(K,J)**2
                END DO
                BDTEST=BDTEST+HALF*CURV*RHO
                IF (BDTEST  <  BDTOL) GOTO 650
            END IF
        END DO
        GOTO 680
    END IF
    NTRITS=NTRITS+1
    !
    !  Severe cancellation is likely to occur if XOPT is too far from XBASE.
    !  If the following test holds, then XBASE is shifted so that XOPT becomes
    !  zero. The appropriate changes are made to BMAT and to the second
    !  derivatives of the current model, beginning with the changes to BMAT
    !  that do not depend on ZMAT. VLAG is used temporarily for working space.
    !
90  IF (DSQ  <=  1.0D-3*XOPTSQ) THEN
        FRACSQ=0.25D0*XOPTSQ
        SUMPQ=ZERO
        DO K=1,NPT
            SUMPQ=SUMPQ+PQ(K)
            SUM=-HALF*XOPTSQ
            DO I=1,N
                SUM=SUM+XPT(K,I)*XOPT(I)
            END DO
            W(NPT+K)=SUM
            TEMP=FRACSQ-HALF*SUM
            DO I=1,N
                W(I)=BMAT(K,I)
                VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
                IP=NPT+I
                DO J=1,I
                    BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
                END DO
            END DO
        END DO
        !
        !     Then the revisions of BMAT that depend on ZMAT are calculated.
        !
        DO JJ=1,NPTM
            SUMZ=ZERO
            SUMW=ZERO
            DO K=1,NPT
                SUMZ=SUMZ+ZMAT(K,JJ)
                VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
                SUMW=SUMW+VLAG(K)
            END DO
            DO J=1,N
                SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
                DO K=1,NPT
                    SUM=SUM+VLAG(K)*XPT(K,J)
                END DO
                W(J)=SUM
                DO K=1,NPT
                    BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
                END DO
            END DO
            DO I=1,N
                IP=I+NPT
                TEMP=W(I)
                DO J=1,I
                    BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
                END DO
            END DO
        END DO
        !
        !  The following instructions complete the shift, including the changes
        !  to the second derivative parameters of the quadratic model.
        !
        IH=0
        DO J=1,N
            W(J)=-HALF*SUMPQ*XOPT(J)
            DO K=1,NPT
                W(J)=W(J)+PQ(K)*XPT(K,J)
                XPT(K,J)=XPT(K,J)-XOPT(J)
            END DO
            DO I=1,J
                IH=IH+1
                HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
                BMAT(NPT+I,J)=BMAT(NPT+J,I)
            END DO
        END DO
        DO I=1,N
            XBASE(I)=XBASE(I)+XOPT(I)
            XNEW(I)=XNEW(I)-XOPT(I)
            SL(I)=SL(I)-XOPT(I)
            SU(I)=SU(I)-XOPT(I)
            XOPT(I)=ZERO
        END DO
        XOPTSQ=ZERO
    END IF
    IF (NTRITS .EQ. 0) GOTO 210
    GOTO 230
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
190 NFSAV=NF
    KBASE=KOPT
    CALL RESCUE (this, N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,&
        &  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,&
        &  VLAG,W,W(N+NP),W(NDIM+NP))
    !
    !  XOPT is updated now in case the branch below to label 720 is taken.
    !  Any updating of GOPT occurs after the branch below to label 20, which
    !  leads to a trust region iteration as does the branch to label 60.
    !
    XOPTSQ=ZERO
    IF (KOPT .NE. KBASE) THEN
        DO I=1,N
            XOPT(I)=XPT(KOPT,I)
            XOPTSQ=XOPTSQ+XOPT(I)**2
        END DO
    END IF
    IF (NF  <  0) THEN
        NF=MAXFUN
        IF (IPRINT  >  0) PRINT 390
        GOTO 720
    END IF
    NRESC=NF
    IF (NFSAV  <  NF) THEN
        NFSAV=NF
        GOTO 20
    END IF
    IF (NTRITS  >  0) GOTO 60
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
210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,&
        &  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
    DO I=1,N
        D(I)=XNEW(I)-XOPT(I)
    END DO
    !
    !  Calculate VLAG and BETA for the current choice of D. The scalar
    !  product of D with XPT(K,.) is going to be held in W(NPT+K) for
    !  use when VQUAD is calculated.
    !
230 DO K=1,NPT
        SUMA=ZERO
        SUMB=ZERO
        SUM=ZERO
        DO J=1,N
            SUMA=SUMA+XPT(K,J)*D(J)
            SUMB=SUMB+XPT(K,J)*XOPT(J)
            SUM=SUM+BMAT(K,J)*D(J)
        END DO
        W(K)=SUMA*(HALF*SUMA+SUMB)
        VLAG(K)=SUM
        W(NPT+K)=SUMA
    END DO
    BETA=ZERO
    DO JJ=1,NPTM
        SUM=ZERO
        DO K=1,NPT
            SUM=SUM+ZMAT(K,JJ)*W(K)
        END DO
        BETA=BETA-SUM*SUM
        DO K=1,NPT
            VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
        END DO
    END DO
    DSQ=ZERO
    BSUM=ZERO
    DX=ZERO
    DO J=1,N
        DSQ=DSQ+D(J)**2
        SUM=ZERO
        DO K=1,NPT
            SUM=SUM+W(K)*BMAT(K,J)
        END DO
        BSUM=BSUM+SUM*D(J)
        JP=NPT+J
        DO I=1,N
            SUM=SUM+BMAT(JP,I)*D(I)
        END DO
        VLAG(JP)=SUM
        BSUM=BSUM+SUM*D(J)
        DX=DX+D(J)*XOPT(J)
    END DO
    BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
    VLAG(KOPT)=VLAG(KOPT)+ONE
    !
    !  If NTRITS is zero, the denominator may be increased by replacing
    !  the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
    !  rounding errors have damaged the chosen denominator.
    !
    IF (NTRITS .EQ. 0) THEN
        DENOM=VLAG(KNEW)**2+ALPHA*BETA
        IF (DENOM  <  CAUCHY .AND. CAUCHY  >  ZERO) THEN
            DO I=1,N
                XNEW(I)=XALT(I)
                D(I)=XNEW(I)-XOPT(I)
            END DO
            CAUCHY=ZERO
            GO TO 230
        END IF
        IF (DENOM  <=  HALF*VLAG(KNEW)**2) THEN
            IF (NF  >  NRESC) GOTO 190
            IF (IPRINT  >  0) PRINT 320
320         FORMAT (/5X,'Return from BOBYQA because of much',&
                &          ' cancellation in a denominator.')
            GOTO 720
        END IF
        !
        !  Alternatively, if NTRITS is positive, then set KNEW to the index of
        !  the next interpolation point to be deleted to make room for a trust
        !  region step. Again RESCUE may be called if rounding errors have damaged
        !  the chosen denominator, which is the reason for attempting to select
        !  KNEW before calculating the next value of the objective function.
        !
    ELSE
        DELSQ=DELTA*DELTA
        SCADEN=ZERO
        BIGLSQ=ZERO
        KNEW=0
        DO K=1,NPT
            IF (K .EQ. KOPT) continue
            HDIAG=ZERO
            DO JJ=1,NPTM
                HDIAG=HDIAG+ZMAT(K,JJ)**2
            END DO
            DEN=BETA*HDIAG+VLAG(K)**2
            DISTSQ=ZERO
            DO J=1,N
                DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
            END DO
            TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
            IF (TEMP*DEN  >  SCADEN) THEN
                SCADEN=TEMP*DEN
                KNEW=K
                DENOM=DEN
            END IF
            BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
        END DO
        IF (SCADEN  <=  HALF*BIGLSQ) THEN
            IF (NF  >  NRESC) GOTO 190
            IF (IPRINT  >  0) PRINT 320
            GOTO 720
        END IF
    END IF
    !
    !  Put the variables for the next calculation of the objective function
    !    in XNEW, with any adjustments for the bounds.
    !
    !
    !  Calculate the value of the objective function at XBASE+XNEW, unless
    !    the limit on the number of calculations of F has been reached.
    !
360 DO I=1,N
        X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
        IF (XNEW(I) .EQ. SL(I)) X(I)=XL(I)
        IF (XNEW(I) .EQ. SU(I)) X(I)=XU(I)
    END DO
    IF (NF  >=  MAXFUN) THEN
        IF (IPRINT  >  0) PRINT 390
390     FORMAT (/4X,'Return from BOBYQA because funkk has been',&
            &      ' called MAXFUN times.')
        BOBYQB = .false.
        GOTO 720
    END IF
    NF=NF+1
    F = this%funkk(this%obj, X(1:N))
    IF (IPRINT .EQ. 3) THEN
        PRINT 400, NF,F,(X(I),I=1,N)
400     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
            &       '    The corresponding X is:'/(2X,5D15.6))
    END IF
    IF (NTRITS .EQ. -1) THEN
        FSAVE=F
        GOTO 720
    END IF
    !
    !  Use the quadratic model to predict the change in F due to the step D,
    !    and set DIFF to the error of this prediction.
    !
    FOPT=FVAL(KOPT)
    VQUAD=ZERO
    IH=0
    DO J=1,N
        VQUAD=VQUAD+D(J)*GOPT(J)
        DO I=1,J
            IH=IH+1
            TEMP=D(I)*D(J)
            IF (I .EQ. J) TEMP=HALF*TEMP
            VQUAD=VQUAD+HQ(IH)*TEMP
        END DO
    END DO
    DO K=1,NPT
        VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
    END DO
    DIFF=F-FOPT-VQUAD
    DIFFC=DIFFB
    DIFFB=DIFFA
    DIFFA=DABS(DIFF)
    IF (DNORM  >  RHO) NFSAV=NF
    !
    !     Pick the next value of DELTA after a trust region step.
    !
    IF (NTRITS  >  0) THEN
        IF (VQUAD  >=  ZERO) THEN
            IF (IPRINT  >  0) PRINT 430
430         FORMAT (/4X,'Return from BOBYQA because a trust',&
                &          ' region step has failed to reduce Q.')
            BOBYQB = .false.
            GOTO 720
        END IF
        RATIO=(F-FOPT)/VQUAD
        IF (RATIO  <=  TENTH) THEN
            DELTA=DMIN1(HALF*DELTA,DNORM)
        ELSEIF (RATIO  <=  0.7D0) THEN
            DELTA=DMAX1(HALF*DELTA,DNORM)
        ELSE
            DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
        END IF
        IF (DELTA  <=  1.5D0*RHO) DELTA=RHO
        !
        !     Recalculate KNEW and DENOM if the new F is less than FOPT.
        !
        IF (F  <  FOPT) THEN
            KSAV=KNEW
            DENSAV=DENOM
            DELSQ=DELTA*DELTA
            SCADEN=ZERO
            BIGLSQ=ZERO
            KNEW=0
            DO K=1,NPT
                HDIAG=ZERO
                DO JJ=1,NPTM
                    HDIAG=HDIAG+ZMAT(K,JJ)**2
                END DO
                DEN=BETA*HDIAG+VLAG(K)**2
                DISTSQ=ZERO
                DO J=1,N
                    DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
                END DO
                TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
                IF (TEMP*DEN  >  SCADEN) THEN
                    SCADEN=TEMP*DEN
                    KNEW=K
                    DENOM=DEN
                END IF
                BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
            END DO
            IF (SCADEN  <=  HALF*BIGLSQ) THEN
                KNEW=KSAV
                DENOM=DENSAV
            END IF
        END IF
    END IF
    !
    !  Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
    !  moved. Also update the second derivative terms of the model.
    !
    CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
    IH=0
    PQOLD=PQ(KNEW)
    PQ(KNEW)=ZERO
    DO I=1,N
        TEMP=PQOLD*XPT(KNEW,I)
        DO J=1,I
            IH=IH+1
            HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
        END DO
    END DO
    DO  JJ=1,NPTM
        TEMP=DIFF*ZMAT(KNEW,JJ)
        DO  K=1,NPT
            PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
        END DO
    END DO
    !
    !  Include the new interpolation point, and make the changes to GOPT at
    !  the old XOPT that are caused by the updating of the quadratic model.
    !
    FVAL(KNEW)=F
    DO  I=1,N
        XPT(KNEW,I)=XNEW(I)
        W(I)=BMAT(KNEW,I)
    END DO
    DO K=1,NPT
        SUMA=ZERO
        DO JJ=1,NPTM
            SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
        END DO
        SUMB=ZERO
        DO J=1,N
            SUMB=SUMB+XPT(K,J)*XOPT(J)
        END DO
        TEMP=SUMA*SUMB
        DO I=1,N
            W(I)=W(I)+TEMP*XPT(K,I)
        END DO
    END DO
    DO I=1,N
        GOPT(I)=GOPT(I)+DIFF*W(I)
    END DO
    !
    !     Update XOPT, GOPT and KOPT if the new calculated F is less than FO
    !
    IF (F  <  FOPT) THEN
        KOPT=KNEW
        XOPTSQ=ZERO
        IH=0
        DO J=1,N
            XOPT(J)=XNEW(J)
            XOPTSQ=XOPTSQ+XOPT(J)**2
            DO I=1,J
                IH=IH+1
                IF (I  <  J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
                GOPT(I)=GOPT(I)+HQ(IH)*D(J)
            END DO
        END DO
        DO K=1,NPT
            TEMP=ZERO
            DO J=1,N
                TEMP=TEMP+XPT(K,J)*D(J)
            END DO
            TEMP=PQ(K)*TEMP
            DO I=1,N
                GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
            END DO
        END DO
    END IF
    !
    !  Calculate the parameters of the least Frobenius norm interpolant to
    !  the current data, the gradient of this interpolant at XOPT being put
    !  into VLAG(NPT+I), I=1,2,...,N.
    !
    IF (NTRITS  >  0) THEN
        DO K=1,NPT
            VLAG(K)=FVAL(K)-FVAL(KOPT)
            W(K)=ZERO
        END DO
        DO J=1,NPTM
            SUM=ZERO
            DO K=1,NPT
                SUM=SUM+ZMAT(K,J)*VLAG(K)
            END DO
            DO K=1,NPT
                W(K)=W(K)+SUM*ZMAT(K,J)
            END DO
        END DO
        DO  K=1,NPT
            SUM=ZERO
            DO J=1,N
                SUM=SUM+XPT(K,J)*XOPT(J)
            END DO
            W(K+NPT)=W(K)
            W(K)=SUM*W(K)
        END DO
        GQSQ=ZERO
        GISQ=ZERO
        DO I=1,N
            SUM=ZERO
            DO K=1,NPT
                SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
            END DO
            IF (XOPT(I) .EQ. SL(I)) THEN
                GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
                GISQ=GISQ+DMIN1(ZERO,SUM)**2
            ELSE IF (XOPT(I) .EQ. SU(I)) THEN
                GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
                GISQ=GISQ+DMAX1(ZERO,SUM)**2
            ELSE
                GQSQ=GQSQ+GOPT(I)**2
                GISQ=GISQ+SUM*SUM
            END IF
            VLAG(NPT+I)=SUM
        END DO
        !
        !  Test whether to replace the new quadratic model by the least Frobenius
        !  norm interpolant, making the replacement if the test is satisfied.
        !
        ITEST=ITEST+1
        IF (GQSQ  <  TEN*GISQ) ITEST=0
        IF (ITEST  >=  3) THEN
            DO I=1,MAX0(NPT,NH)
                IF (I  <=  N) GOPT(I)=VLAG(NPT+I)
                IF (I  <=  NPT) PQ(I)=W(NPT+I)
                IF (I  <=  NH) HQ(I)=ZERO
                ITEST=0
            END DO
        END IF
    END IF
    !
    !  If a trust region step has provided a sufficient decrease in F, then
    !  branch for another trust region calculation. The case NTRITS=0 occurs
    !  when the new interpolation point was reached by an alternative step.
    !
    IF (NTRITS .EQ. 0) GOTO 60
    IF (F  <=  FOPT+TENTH*VQUAD) GOTO 60
    !
    !     Alternatively, find out if the interpolation points are close enough
    !       to the best point so far.
    !
    DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
650 KNEW=0
    DO K=1,NPT
        SUM=ZERO
        DO J=1,N
            SUM=SUM+(XPT(K,J)-XOPT(J))**2
        END DO
        IF (SUM  >  DISTSQ) THEN
            KNEW=K
            DISTSQ=SUM
        END IF
    END DO
    !
    !  If KNEW is positive, then ALTMOV finds alternative new positions for
    !  the KNEW-th interpolation point within distance ADELT of XOPT. It is
    !  reached via label 90. Otherwise, there is a branch to label 60 for
    !  another trust region iteration, unless the calculations with the
    !  current RHO are complete.
    !
    IF (KNEW  >  0) THEN
        DIST=DSQRT(DISTSQ)
        IF (NTRITS .EQ. -1) THEN
            DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
            IF (DELTA  <=  1.5D0*RHO) DELTA=RHO
        END IF
        NTRITS=0
        ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
        DSQ=ADELT*ADELT
        GOTO 90
    END IF
    IF (NTRITS .EQ. -1) GOTO 680
    IF (RATIO  >  ZERO) GOTO 60
    IF (DMAX1(DELTA,DNORM)  >  RHO) GOTO 60
    !
    !  The calculations with the current value of RHO are complete. Pick the
    !    next values of RHO and DELTA.
    !
680 IF (RHO  >  RHOEND .and. (this%FVAL_Converge_difference==0._dp .or. &
        abs(FVAL(KOPT)-this%Last_bestfit) > this%FVAL_Converge_difference)) THEN
        DELTA=HALF*RHO
        RATIO=RHO/RHOEND
        IF (RATIO  <=  16.0D0) THEN
            RHO=RHOEND
        ELSE IF (RATIO  <=  250.0D0) THEN
            RHO=DSQRT(RATIO)*RHOEND
        ELSE
            RHO=TENTH*RHO
        END IF
        DELTA=DMAX1(DELTA,RHO)
        this%Last_bestfit = FVAL(KOPT)
        IF (IPRINT  >=  2) THEN
            IF (IPRINT  >=  3) PRINT 690
690         FORMAT (5X)
            PRINT 700, RHO,NF
700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',&
                &          ' function values =',I6)
            PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
710         FORMAT (4X,'Least value of F =',1PD23.15,9X,&
                &          'The corresponding X is:'/(2X,5D15.6))
        END IF
        NTRITS=0
        NFSAV=NF
        GOTO 60
    END IF
    !
    !     Return from the calculation, after another Newton-Raphson step, if
    !       it is too short to have been tried before.
    !
    IF (NTRITS .EQ. -1) GOTO 360
720 IF (FVAL(KOPT)  <=  FSAVE) THEN
        DO I=1,N
            X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
            IF (XOPT(I) .EQ. SL(I)) X(I)=XL(I)
            IF (XOPT(I) .EQ. SU(I)) X(I)=XU(I)
        END DO
        F=FVAL(KOPT)
    END IF
    IF (IPRINT  >=  1) THEN
        if (RHO  >  RHOEND) write(*,*) '    FVAL_Converge_difference reached'
        PRINT 740, NF
740     FORMAT (/4X,'At the return from BOBYQA',5X,&
            &      'Number of function values =',I6)
        PRINT 710, F,(X(I),I=1,N)
    END IF
    end function BOBYQB


    SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,&
        &  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
    REAL(dp) Adelt , Alpha , bigstp , Cauchy , const , csave ,   &
        & curv , dderiv , diff , distsq , ggfree , gw , ha ,  &
        & half , one , predsq , presav
    REAL(dp) scale , slbd , step , stpsav , subd , sumin ,    &
        & temp , tempa , tempb , tempd , vlag , wfixsq , wsqsav, zero
    INTEGER i , ibdsav , iflag , ilbd , isbd , iubd , j , k , Knew ,  &
        & Kopt , ksav , N , Ndim , Npt
    REAL(dp) XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),&
        &  SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
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
    HALF=0.5D0
    ONE=1.0D0
    ZERO=0.0D0
    CONST=ONE+DSQRT(2.0D0)
    DO K=1,NPT
        HCOL(K)=ZERO
    END DO
    DO J=1,NPT-N-1
        TEMP=ZMAT(KNEW,J)
        DO K=1,NPT
            HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
        END DO
    END DO
    ALPHA=HCOL(KNEW)
    HA=HALF*ALPHA
    !
    !     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
    !
    DO I=1,N
        GLAG(I)=BMAT(KNEW,I)
    END DO
    DO  K=1,NPT
        TEMP=ZERO
        DO J=1,N
            TEMP=TEMP+XPT(K,J)*XOPT(J)
        END DO
        TEMP=HCOL(K)*TEMP
        DO I=1,N
            GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
        END DO
    END DO
    !
    !  Search for a large denominator along the straight lines through XOPT
    !  and another interpolation point. SLBD and SUBD will be lower and upper
    !  bounds on the step along each of these lines in turn. PREDSQ will be
    !  set to the square of the predicted denominator for each line. PRESAV
    !  will be set to the largest admissible value of PREDSQ that occurs.
    PRESAV=ZERO
    DO K=1,NPT
        IF (K .EQ. KOPT) continue
        DDERIV=ZERO
        DISTSQ=ZERO
        DO I=1,N
            TEMP=XPT(K,I)-XOPT(I)
            DDERIV=DDERIV+GLAG(I)*TEMP
            DISTSQ=DISTSQ+TEMP*TEMP
        END DO
        SUBD=ADELT/DSQRT(DISTSQ)
        SLBD=-SUBD
        ILBD=0
        IUBD=0
        SUMIN=DMIN1(ONE,SUBD)
        !
        !     Revise SLBD and SUBD if necessary because of the bounds in SL and
        !
        DO I=1,N
            TEMP=XPT(K,I)-XOPT(I)
            IF (TEMP  >  ZERO) THEN
                IF (SLBD*TEMP  <  SL(I)-XOPT(I)) THEN
                    SLBD=(SL(I)-XOPT(I))/TEMP
                    ILBD=-I
                END IF
                IF (SUBD*TEMP  >  SU(I)-XOPT(I)) THEN
                    SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
                    IUBD=I
                END IF
            ELSE IF (TEMP  <  ZERO) THEN
                IF (SLBD*TEMP  >  SU(I)-XOPT(I)) THEN
                    SLBD=(SU(I)-XOPT(I))/TEMP
                    ILBD=I
                END IF
                IF (SUBD*TEMP  <  SL(I)-XOPT(I)) THEN
                    SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
                    IUBD=-I
                END IF
            END IF
        END DO
        !
        !  Seek a large modulus of the KNEW-th Lagrange function when the index
        !  of the other interpolation point on the line through XOPT is KNEW.
        !
        IF (K .EQ. KNEW) THEN
            DIFF=DDERIV-ONE
            STEP=SLBD
            VLAG=SLBD*(DDERIV-SLBD*DIFF)
            ISBD=ILBD
            TEMP=SUBD*(DDERIV-SUBD*DIFF)
            IF (DABS(TEMP)  >  DABS(VLAG)) THEN
                STEP=SUBD
                VLAG=TEMP
                ISBD=IUBD
            END IF
            TEMPD=HALF*DDERIV
            TEMPA=TEMPD-DIFF*SLBD
            TEMPB=TEMPD-DIFF*SUBD
            IF (TEMPA*TEMPB  <  ZERO) THEN
                TEMP=TEMPD*TEMPD/DIFF
                IF (DABS(TEMP)  >  DABS(VLAG)) THEN
                    STEP=TEMPD/DIFF
                    VLAG=TEMP
                    ISBD=0
                END IF
            END IF
            !
            !     Search along each of the other lines through XOPT and another point
            !
        ELSE
            STEP=SLBD
            VLAG=SLBD*(ONE-SLBD)
            ISBD=ILBD
            TEMP=SUBD*(ONE-SUBD)
            IF (DABS(TEMP)  >  DABS(VLAG)) THEN
                STEP=SUBD
                VLAG=TEMP
                ISBD=IUBD
            END IF
            IF (SUBD  >  HALF) THEN
                IF (DABS(VLAG)  <  0.25D0) THEN
                    STEP=HALF
                    VLAG=0.25D0
                    ISBD=0
                END IF
            END IF
            VLAG=VLAG*DDERIV
        END IF
        !
        !     Calculate PREDSQ for the current line search and maintain PRESAV.
        !
        TEMP=STEP*(ONE-STEP)*DISTSQ
        PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
        IF (PREDSQ  >  PRESAV) THEN
            PRESAV=PREDSQ
            KSAV=K
            STPSAV=STEP
            IBDSAV=ISBD
        END IF
    END DO
    !
    !     Construct XNEW in a way that satisfies the bound constraints exact
    !
    DO I=1,N
        TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
        XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
    END DO
    IF (IBDSAV  <  0) XNEW(-IBDSAV)=SL(-IBDSAV)
    IF (IBDSAV  >  0) XNEW(IBDSAV)=SU(IBDSAV)
    !
    !  Prepare for the iterative method that assembles the constrained Cauchy
    !  step in W. The sum of squares of the fixed components of W is formed in
    !  WFIXSQ, and the free components of W are set to BIGSTP.
    !
    BIGSTP=ADELT+ADELT
    IFLAG=0
100 WFIXSQ=ZERO
    GGFREE=ZERO
    DO I=1,N
        W(I)=ZERO
        TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
        TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
        IF (TEMPA  >  ZERO .OR. TEMPB  <  ZERO) THEN
            W(I)=BIGSTP
            GGFREE=GGFREE+GLAG(I)**2
        END IF
    END DO
    IF (GGFREE .EQ. ZERO) THEN
        CAUCHY=ZERO
        return
    END IF
    !
    !     Investigate whether more components of W can be fixed.
    !
120 TEMP=ADELT*ADELT-WFIXSQ
    IF (TEMP  >  ZERO) THEN
        WSQSAV=WFIXSQ
        STEP=DSQRT(TEMP/GGFREE)
        GGFREE=ZERO
        DO I=1,N
            IF (W(I) .EQ. BIGSTP) THEN
                TEMP=XOPT(I)-STEP*GLAG(I)
                IF (TEMP  <=  SL(I)) THEN
                    W(I)=SL(I)-XOPT(I)
                    WFIXSQ=WFIXSQ+W(I)**2
                ELSE IF (TEMP  >=  SU(I)) THEN
                    W(I)=SU(I)-XOPT(I)
                    WFIXSQ=WFIXSQ+W(I)**2
                ELSE
                    GGFREE=GGFREE+GLAG(I)**2
                END IF
            END IF
        END DO
        IF (WFIXSQ  >  WSQSAV .AND. GGFREE  >  ZERO) GOTO 120
    END IF
    !
    !     Set the remaining free components of W and all components of XALT,
    !     except that W may be scaled later.
    !
    GW=ZERO
    DO I=1,N
        IF (W(I) .EQ. BIGSTP) THEN
            W(I)=-STEP*GLAG(I)
            XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
        ELSE IF (W(I) .EQ. ZERO) THEN
            XALT(I)=XOPT(I)
        ELSE IF (GLAG(I)  >  ZERO) THEN
            XALT(I)=SL(I)
        ELSE
            XALT(I)=SU(I)
        END IF
        GW=GW+GLAG(I)*W(I)
    END DO
    !
    !  Set CURV to the curvature of the KNEW-th Lagrange function along W.
    !  Scale W by a factor less than one if that can reduce the modulus of
    !  the Lagrange function at XOPT+W. Set CAUCHY to the final value of
    !  the square of this function.
    !
    CURV=ZERO
    DO K=1,NPT
        TEMP=ZERO
        DO  J=1,N
            TEMP=TEMP+XPT(K,J)*W(J)
        END DO
        CURV=CURV+HCOL(K)*TEMP*TEMP
    END DO
    IF (IFLAG .EQ. 1) CURV=-CURV
    IF (CURV  >  -GW .AND. CURV  <  -CONST*GW) THEN
        SCALE=-GW/CURV
        DO I=1,N
            TEMP=XOPT(I)+SCALE*W(I)
            XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
        END DO
        CAUCHY=(HALF*GW*SCALE)**2
    ELSE
        CAUCHY=(GW+HALF*CURV)**2
    END IF
    !
    !  If IFLAG is zero, then XALT is calculated as before after reversing
    !  the sign of GLAG. Thus two XALT vectors become available. The one that
    !  is chosen is the one that gives the larger value of CAUCHY.
    !
    IF (IFLAG .EQ. 0) THEN
        DO I=1,N
            GLAG(I)=-GLAG(I)
            W(N+I)=XALT(I)
        END DO
        CSAVE=CAUCHY
        IFLAG=1
        GOTO 100
    END IF
    IF (CSAVE  >  CAUCHY) THEN
        DO I=1,N
            XALT(I)=W(N+I)
        END DO
        CAUCHY=CSAVE
    END IF

    END SUBROUTINE ALTMOV


    SUBROUTINE PRELIM (this,N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,&
        &  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
    class(TBOBYQA) this
    REAL(dp) diff , f , fbeg , half , one ,   &
        & recip , Rhobeg , rhosq , stepa , stepb, temp , two,  zero
    INTEGER i , ih , Iprint , ipt , itemp , j , jpt , k , Kopt ,      &
        & Maxfun , N , Ndim , Nf , nfm , nfx , np , Npt
    REAL(DP) X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),&
        &  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)

    !
    !  The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
    !    same as the corresponding arguments in SUBROUTINE BOBYQA.
    !  The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
    !    are the same as the corresponding arguments in BOBYQB, the elements
    !    of SL and SU being set in BOBYQA.
    !  GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
    !    it is set by PRELIM to the gradient of the quadratic model at XBASE.
    !    If XOPT is nonzero, BOBYQB will change it to its usual value later.
    !  NF is maintaned as the number of calls of CALFUN so far.
    !  KOPT will be such that the least calculated value of F so far is at
    !    the point XPT(KOPT,.)+XBASE in the space of the variables.
    !
    !  SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    !  BMAT and ZMAT for the first iteration, and it maintains the values of
    !  NF and KOPT. The vector X is also changed by PRELIM.
    !
    !  Set some constants.
    !
    HALF=0.5D0
    ONE=1.0D0
    TWO=2.0D0
    ZERO=0.0D0
    RHOSQ=RHOBEG*RHOBEG
    RECIP=ONE/RHOSQ
    NP=N+1
    !
    !     Set XBASE to the initial vector of variables, and set the initial
    !     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    !
    DO J=1,N
        XBASE(J)=X(J)
        DO K=1,NPT
            XPT(K,J)=ZERO
        END DO
        DO I=1,NDIM
            BMAT(I,J)=ZERO
        END DO
    END DO
    DO IH=1,(N*NP)/2
        HQ(IH)=ZERO
    END DO
    DO K=1,NPT
        PQ(K)=ZERO
        DO J=1,NPT-NP
            ZMAT(K,J)=ZERO
        END DO
    END DO
    !
    !  Begin the initialization procedure. NF becomes one more than the number
    !  of function values so far. The coordinates of the displacement of the
    !  next initial interpolation point from XBASE are set in XPT(NF+1,.).
    !
    NF=0
50  NFM=NF
    NFX=NF-N
    NF=NF+1
    IF (NFM  <=  2*N) THEN
        IF (NFM  >=  1 .AND. NFM  <=  N) THEN
            STEPA=RHOBEG
            IF (SU(NFM) .EQ. ZERO) STEPA=-STEPA
            XPT(NF,NFM)=STEPA
        ELSE IF (NFM  >  N) THEN
            STEPA=XPT(NF-N,NFX)
            STEPB=-RHOBEG
            IF (SL(NFX) .EQ. ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
            IF (SU(NFX) .EQ. ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
            XPT(NF,NFX)=STEPB
        END IF
    ELSE
        ITEMP=(NFM-NP)/N
        JPT=NFM-ITEMP*N-N
        IPT=JPT+ITEMP
        IF (IPT  >  N) THEN
            ITEMP=JPT
            JPT=IPT-N
            IPT=ITEMP
        END IF
        XPT(NF,IPT)=XPT(IPT+1,IPT)
        XPT(NF,JPT)=XPT(JPT+1,JPT)
    END IF
    !
    !     Calculate the next value of F. The least function value so far and
    !     its index are required.
    !
    DO J=1,N
        X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
        IF (XPT(NF,J) .EQ. SL(J)) X(J)=XL(J)
        IF (XPT(NF,J) .EQ. SU(J)) X(J)=XU(J)
    END DO
    F = this%funkk(this%obj,X(1:N))
    IF (IPRINT .EQ. 3) THEN
        PRINT 70, NF,F,(X(I),I=1,N)
70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
            &       '    The corresponding X is:'/(2X,5D15.6))
    END IF
    FVAL(NF)=F
    IF (NF .EQ. 1) THEN
        FBEG=F
        KOPT=1
    ELSE IF (F  <  FVAL(KOPT)) THEN
        KOPT=NF
    END IF
    !
    !  Set the nonzero initial elements of BMAT and the quadratic model in the
    !  cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
    !  of the NF-th and (NF-N)-th interpolation points may be switched, in
    !  order that the function value at the first of them contributes to the
    !  off-diagonal second derivative terms of the initial quadratic model.
    !
    IF (NF  <=  2*N+1) THEN
        IF (NF  >=  2 .AND. NF  <=  N+1) THEN
            GOPT(NFM)=(F-FBEG)/STEPA
            IF (NPT  <  NF+N) THEN
                BMAT(1,NFM)=-ONE/STEPA
                BMAT(NF,NFM)=ONE/STEPA
                BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
            END IF
        ELSE IF (NF  >=  N+2) THEN
            IH=(NFX*(NFX+1))/2
            TEMP=(F-FBEG)/STEPB
            DIFF=STEPB-STEPA
            HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
            GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
            IF (STEPA*STEPB  <  ZERO) THEN
                IF (F  <  FVAL(NF-N)) THEN
                    FVAL(NF)=FVAL(NF-N)
                    FVAL(NF-N)=F
                    IF (KOPT .EQ. NF) KOPT=NF-N
                    XPT(NF-N,NFX)=STEPB
                    XPT(NF,NFX)=STEPA
                END IF
            END IF
            BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
            BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
            BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
            ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
            ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
            ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
        END IF
        !
        !     Set the off-diagonal second derivatives of the Lagrange functions
        !     the initial quadratic model.
        !
    ELSE
        IH=(IPT*(IPT-1))/2+JPT
        ZMAT(1,NFX)=RECIP
        ZMAT(NF,NFX)=RECIP
        ZMAT(IPT+1,NFX)=-RECIP
        ZMAT(JPT+1,NFX)=-RECIP
        TEMP=XPT(NF,IPT)*XPT(NF,JPT)
        HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
    END IF
    IF (NF  <  NPT .AND. NF  <  MAXFUN) GOTO 50

    END SUBROUTINE PRELIM


    SUBROUTINE RESCUE (this,N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,&
        &  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,&
        &  KOPT,VLAG,PTSAUX,PTSID,W)
    class(TBOBYQA) this

    REAL(dp) beta , bsum , Delta , den , denom , diff , distsq , &
        & dsqmin , f , fbase , half , hdiag , one
    REAL(dp) sfrac , sum , sumpq , temp , vlmxsq ,     &
        & vquad , winc , xp , xq , zero
    INTEGER i , ih , ihp , ihq , ip , Iprint , iq , iw , j , jp ,     &
        & jpn , k , knew , kold , Kopt , kpt , Maxfun , N , Ndim ,  &
        & Nf
    INTEGER np , Npt , nptm , nrem
    REAL(dp) XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),&
        &  GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),&
        &  VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
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
    !    p=0, respectively.
    !  The first NDIM+NPT elements of the array W are used for working space.
    !  The final elements of BMAT and ZMAT are set in a well-conditioned way
    !    to the values that are appropriate for the new interpolation points.
    !  The elements of GOPT, HQ and PQ are also revised to the values that are
    !    appropriate to the final quadratic model.
    !
    !     Set some constants.
    !
    HALF=0.5D0
    ONE=1.0D0
    ZERO=0.0D0
    NP=N+1
    SFRAC=HALF/NP
    NPTM=NPT-NP
    !
    !  Shift the interpolation points so that XOPT becomes the origin, and set
    !  the elements of ZMAT to zero. The value of SUMPQ is required in the
    !  updating of HQ below. The squares of the distances from XOPT to the
    !  other interpolation points are set at the end of W. Increments of WINC
    !  may be added later to these squares to balance the consideration of
    !  the choice of point that is going to become current.
    !
    SUMPQ=ZERO
    WINC=ZERO
    DO K=1,NPT
        DISTSQ=ZERO
        DO J=1,N
            XPT(K,J)=XPT(K,J)-XOPT(J)
            DISTSQ=DISTSQ+XPT(K,J)**2
        END DO
        SUMPQ=SUMPQ+PQ(K)
        W(NDIM+K)=DISTSQ
        WINC=DMAX1(WINC,DISTSQ)
        DO J=1,NPTM
            ZMAT(K,J)=ZERO
        END DO
    END DO
    !
    !     Update HQ so that HQ and PQ define the second derivatives of the model
    !     after XBASE has been shifted to the trust region centre.
    !
    IH=0
    DO J=1,N
        W(J)=HALF*SUMPQ*XOPT(J)
        DO K=1,NPT
            W(J)=W(J)+PQ(K)*XPT(K,J)
        END DO
        DO I=1,J
            IH=IH+1
            HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
        END DO
    END DO
    !
    !     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
    !     also set the elements of PTSAUX.
    !
    DO J=1,N
        XBASE(J)=XBASE(J)+XOPT(J)
        SL(J)=SL(J)-XOPT(J)
        SU(J)=SU(J)-XOPT(J)
        XOPT(J)=ZERO
        PTSAUX(1,J)=DMIN1(DELTA,SU(J))
        PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
        IF (PTSAUX(1,J)+PTSAUX(2,J)  <  ZERO) THEN
            TEMP=PTSAUX(1,J)
            PTSAUX(1,J)=PTSAUX(2,J)
            PTSAUX(2,J)=TEMP
        END IF
        IF (DABS(PTSAUX(2,J))  <  HALF*DABS(PTSAUX(1,J))) THEN
            PTSAUX(2,J)=HALF*PTSAUX(1,J)
        END IF
        DO I=1,NDIM
            BMAT(I,J)=ZERO
        END DO
    END DO
    FBASE=FVAL(KOPT)
    !
    !     Set the identifiers of the artificial interpolation points that ard
    !     along a coordinate direction from XOPT, and set the corresponding
    !     nonzero elements of BMAT and ZMAT.
    !
    PTSID(1)=SFRAC
    DO J=1,N
        JP=J+1
        JPN=JP+N
        PTSID(JP)=J+SFRAC
        IF (JPN  <=  NPT) THEN
            PTSID(JPN)=real(J,dp)/NP+SFRAC
            TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
            BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
            BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
            BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
            ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
            ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
            ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
        ELSE
            BMAT(1,J)=-ONE/PTSAUX(1,J)
            BMAT(JP,J)=ONE/PTSAUX(1,J)
            BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
        END IF
    END DO
    !
    !     Set any remaining identifiers with their nonzero elements of ZMAT.
    !
    IF (NPT  >=  N+NP) THEN
        DO K=2*NP,NPT
            IW=(real(K-NP,dp)-HALF)/N
            IP=K-NP-IW*N
            IQ=IP+IW
            IF (IQ  >  N) IQ=IQ-N
            PTSID(K)=IP+real(IQ,dp)/NP+SFRAC
            TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
            ZMAT(1,K-NP)=TEMP
            ZMAT(IP+1,K-NP)=-TEMP
            ZMAT(IQ+1,K-NP)=-TEMP
            ZMAT(K,K-NP)=TEMP
        END DO
    END IF
    NREM=NPT
    KOLD=1
    KNEW=KOPT
    !
    !     Reorder the provisional points in the way that exchanges PTSID(KOL
    !     with PTSID(KNEW).
    !
80  DO J=1,N
        TEMP=BMAT(KOLD,J)
        BMAT(KOLD,J)=BMAT(KNEW,J)
        BMAT(KNEW,J)=TEMP
    END DO
    DO J=1,NPTM
        TEMP=ZMAT(KOLD,J)
        ZMAT(KOLD,J)=ZMAT(KNEW,J)
        ZMAT(KNEW,J)=TEMP
    END DO
    PTSID(KOLD)=PTSID(KNEW)
    PTSID(KNEW)=ZERO
    W(NDIM+KNEW)=ZERO
    NREM=NREM-1
    IF (KNEW .NE. KOPT) THEN
        TEMP=VLAG(KOLD)
        VLAG(KOLD)=VLAG(KNEW)
        VLAG(KNEW)=TEMP
        !
        !  Update the BMAT and ZMAT matrices so that the status of the KNEW-th
        !  interpolation point can be changed from provisional to original. The
        !  return occurs if all the original points are reinstated.
        !  The nonnegative values of W(NDIM+K) are required in the search below.
        !
        CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
        IF (NREM .EQ. 0) return
        DO K=1,NPT
            W(NDIM+K)=DABS(W(NDIM+K))
        END DO
    END IF
    !
    !  Pick the index KNEW of an original interpolation point that has not
    !  yet replaced one of the provisional interpolation points, giving
    !  attention to the closeness to XOPT and to previous tries with KNEW.
    !
120 DSQMIN=ZERO
    DO  K=1,NPT
        IF (W(NDIM+K)  >  ZERO) THEN
            IF (DSQMIN .EQ. ZERO .OR. W(NDIM+K)  <  DSQMIN) THEN
                KNEW=K
                DSQMIN=W(NDIM+K)
            END IF
        END IF
    END DO
    IF (DSQMIN .EQ. ZERO) GOTO 260
    !
    !     Form the W-vector of the chosen original interpolation point.
    !
    DO J=1,N
        W(NPT+J)=XPT(KNEW,J)
    END DO
    DO K=1,NPT
        SUM=ZERO
        IF (K .EQ. KOPT) THEN
            CONTINUE
        ELSE IF (PTSID(K) .EQ. ZERO) THEN
            DO  J=1,N
                SUM=SUM+W(NPT+J)*XPT(K,J)
            END DO
        ELSE
            IP=PTSID(K)
            IF (IP  >  0) SUM=W(NPT+IP)*PTSAUX(1,IP)
            IQ=NP*PTSID(K)-IP*NP
            IF (IQ  >  0) THEN
                IW=1
                IF (IP .EQ. 0) IW=2
                SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
            END IF
        END IF
        W(K)=HALF*SUM*SUM
    END DO
    !
    !     Calculate VLAG and BETA for the required updating of the H matrix
    !     XPT(KNEW,.) is reinstated in the set of interpolation points.
    !
    DO K=1,NPT
        SUM=ZERO
        DO  J=1,N
            SUM=SUM+BMAT(K,J)*W(NPT+J)
        END DO
        VLAG(K)=SUM
    END DO
    BETA=ZERO
    DO J=1,NPTM
        SUM=ZERO
        DO K=1,NPT
            SUM=SUM+ZMAT(K,J)*W(K)
        END DO
        BETA=BETA-SUM*SUM
        DO K=1,NPT
            VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
        END DO
    END DO
    BSUM=ZERO
    DISTSQ=ZERO
    DO J=1,N
        SUM=ZERO
        DO K=1,NPT
            SUM=SUM+BMAT(K,J)*W(K)
        END DO
        JP=J+NPT
        BSUM=BSUM+SUM*W(JP)
        DO IP=NPT+1,NDIM
            SUM=SUM+BMAT(IP,J)*W(IP)
        END DO
        BSUM=BSUM+SUM*W(JP)
        VLAG(JP)=SUM
        DISTSQ=DISTSQ+XPT(KNEW,J)**2
    END DO
    BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
    VLAG(KOPT)=VLAG(KOPT)+ONE
    !
    !  KOLD is set to the index of the provisional interpolation point that is
    !  going to be deleted to make way for the KNEW-th original interpolation
    !  point. The choice of KOLD is governed by the avoidance of a small value
    !  of the denominator in the updating calculation of UPDATE.
    !
    DENOM=ZERO
    VLMXSQ=ZERO
    DO K=1,NPT
        IF (PTSID(K) .NE. ZERO) THEN
            HDIAG=ZERO
            DO J=1,NPTM
                HDIAG=HDIAG+ZMAT(K,J)**2
            END DO
            DEN=BETA*HDIAG+VLAG(K)**2
            IF (DEN  >  DENOM) THEN
                KOLD=K
                DENOM=DEN
            END IF
        END IF
        VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
    END DO
    IF (DENOM  <=  1.0D-2*VLMXSQ) THEN
        W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
        GOTO 120
    END IF
    GOTO 80
    !
    !  When label 260 is reached, all the final positions of the interpolation
    !  points have been chosen although any changes have not been included yet
    !  in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
    !  from the shift of XBASE, the updating of the quadratic model remains to
    !  be done. The following cycle through the new interpolation points begins
    !  by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
    !  except that a RETURN occurs if MAXFUN prohibits another value of F.
    !
260 DO KPT=1,NPT
        IF (PTSID(KPT) .EQ. ZERO) continue
        IF (NF  >=  MAXFUN) THEN
            NF=-1
            return
        END IF
        IH=0
        DO J=1,N
            W(J)=XPT(KPT,J)
            XPT(KPT,J)=ZERO
            TEMP=PQ(KPT)*W(J)
            DO I=1,J
                IH=IH+1
                HQ(IH)=HQ(IH)+TEMP*W(I)
            END DO
        END DO
        PQ(KPT)=ZERO
        IP=PTSID(KPT)
        IQ=NP*PTSID(KPT)-IP*NP
        IF (IP  >  0) THEN
            XP=PTSAUX(1,IP)
            XPT(KPT,IP)=XP
        END IF
        IF (IQ  >  0) THEN
            XQ=PTSAUX(1,IQ)
            IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
            XPT(KPT,IQ)=XQ
        END IF
        !
        !     Set VQUAD to the value of the current model at the new point.
        !
        VQUAD=FBASE
        IF (IP  >  0) THEN
            IHP=(IP+IP*IP)/2
            VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
        END IF
        IF (IQ  >  0) THEN
            IHQ=(IQ+IQ*IQ)/2
            VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
            IF (IP  >  0) THEN
                IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
                VQUAD=VQUAD+XP*XQ*HQ(IW)
            END IF
        END IF
        DO K=1,NPT
            TEMP=ZERO
            IF (IP  >  0) TEMP=TEMP+XP*XPT(K,IP)
            IF (IQ  >  0) TEMP=TEMP+XQ*XPT(K,IQ)
            VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
        END DO
        !
        !  Calculate F at the new interpolation point, and set DIFF to the factor
        !  that is going to multiply the KPT-th Lagrange function when the model
        !  is updated to provide interpolation to the new function value.
        !
        DO I=1,N
            W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
            IF (XPT(KPT,I) .EQ. SL(I)) W(I)=XL(I)
            IF (XPT(KPT,I) .EQ. SU(I)) W(I)=XU(I)
        END DO
        NF=NF+1
        F = this%funkk(this%obj, W(1:N))
        IF (IPRINT .EQ. 3) THEN
            PRINT 300, NF,F,(W(I),I=1,N)
300         FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
                &      '    The corresponding X is:'/(2X,5D15.6))
        END IF
        FVAL(KPT)=F
        IF (F  <  FVAL(KOPT)) KOPT=KPT
        DIFF=F-VQUAD
        !
        !     Update the quadratic model. The RETURN from the subroutine occurs when
        !     all the new interpolation points are included in the model.
        !
        DO I=1,N
            GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
        END DO
        DO K=1,NPT
            SUM=ZERO
            DO J=1,NPTM
                SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
            END DO
            TEMP=DIFF*SUM
            IF (PTSID(K) .EQ. ZERO) THEN
                PQ(K)=PQ(K)+TEMP
            ELSE
                IP=PTSID(K)
                IQ=NP*PTSID(K)-IP*NP
                IHQ=(IQ*IQ+IQ)/2
                IF (IP .EQ. 0) THEN
                    HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
                ELSE
                    IHP=(IP*IP+IP)/2
                    HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
                    IF (IQ  >  0) THEN
                        HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                        IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                        HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
                    END IF
                END IF
            END IF
        END DO
        PTSID(KPT)=ZERO
    END DO
    END SUBROUTINE RESCUE


    SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,&
        &  XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
    REAL(dp) angbd , angt , beta , blen , Crvmin , cth , delsq ,    &
        & Delta , dhd , dhs , dredg , dredsq , ds , Dsq , ggsav ,    &
        & gredsq , half
    REAL(dp) one , onemin , qred , ratio ,        &
        & rdnext , rdprev , redmax , rednew , redsav , resid ,  &
        & sdec , shs , sqstp , sredg
    REAL(dp) ssq , stepsq , sth , stplen , temp , tempa , tempb ,  &
        & xsav , xsum , zero
    INTEGER i , iact , ih , isav , itcsav , iterc , itermax , iu , j ,&
        & k , N , nact , Npt
    REAL(DP) XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),&
        &  XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
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
    !    SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
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
    HALF=0.5D0
    ONE=1.0D0
    ONEMIN=-1.0D0
    ZERO=0.0D0
    !
    !  The sign of GOPT(I) gives the sign of the change to the I-th variable
    !  that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
    !  or not to fix the I-th variable at one of its bounds initially, with
    !  NACT being set to the number of fixed variables. D and GNEW are also
    !  set for the first iteration. DELSQ is the upper bound on the sum of
    !  squares of the free variables. QRED is the reduction in Q so far.
    !
    ITERC=0
    NACT=0
    SQSTP=ZERO
    DO I=1,N
        XBDI(I)=ZERO
        IF (XOPT(I)  <=  SL(I)) THEN
            IF (GOPT(I)  >=  ZERO) XBDI(I)=ONEMIN
        ELSE IF (XOPT(I)  >=  SU(I)) THEN
            IF (GOPT(I)  <=  ZERO) XBDI(I)=ONE
        END IF
        IF (XBDI(I) .NE. ZERO) NACT=NACT+1
        D(I)=ZERO
        GNEW(I)=GOPT(I)
    END DO
    DELSQ=DELTA*DELTA
    QRED=ZERO
    CRVMIN=ONEMIN
    !
    !  Set the next search direction of the conjugate gradient method. It is
    !  the steepest descent direction initially and when the iterations are
    !  restarted because a variable has just been fixed by a bound, and of
    !  course the components of the fixed variables are zero. ITERMAX is an
    !  upper bound on the indices of the conjugate gradient iterations.
    !
20  BETA=ZERO
30  STEPSQ=ZERO
    DO I=1,N
        IF (XBDI(I) .NE. ZERO) THEN
            S(I)=ZERO
        ELSE IF (BETA .EQ. ZERO) THEN
            S(I)=-GNEW(I)
        ELSE
            S(I)=BETA*S(I)-GNEW(I)
        END IF
        STEPSQ=STEPSQ+S(I)**2
    END DO
    IF (STEPSQ .EQ. ZERO) GOTO 190
    IF (BETA .EQ. ZERO) THEN
        GREDSQ=STEPSQ
        ITERMAX=ITERC+N-NACT
    END IF
    IF (GREDSQ*DELSQ  <=  1.0D-4*QRED*QRED) GO TO 190
    !
    !  Multiply the search direction by the second derivative matrix of Q and
    !  calculate some scalars for the choice of steplength. Then set BLEN to
    !  the length of the the step to the trust region boundary and STPLEN to
    !  the steplength, ignoring the simple bounds.
    !
    GOTO 210
50  RESID=DELSQ
    DS=ZERO
    SHS=ZERO
    DO 60 I=1,N
        IF (XBDI(I) .EQ. ZERO) THEN
            RESID=RESID-D(I)**2
            DS=DS+S(I)*D(I)
            SHS=SHS+S(I)*HS(I)
        END IF
60  CONTINUE
    IF (RESID  <=  ZERO) GOTO 90
    TEMP=DSQRT(STEPSQ*RESID+DS*DS)
    IF (DS  <  ZERO) THEN
        BLEN=(TEMP-DS)/STEPSQ
    ELSE
        BLEN=RESID/(TEMP+DS)
    END IF
    STPLEN=BLEN
    IF (SHS  >  ZERO) THEN
        STPLEN=DMIN1(BLEN,GREDSQ/SHS)
    END IF

    !
    !     Reduce STPLEN if necessary in order to preserve the simple bounds,
    !     letting IACT be the index of the new constrained variable.
    !
    IACT=0
    DO 70 I=1,N
        IF (S(I) .NE. ZERO) THEN
            XSUM=XOPT(I)+D(I)
            IF (S(I)  >  ZERO) THEN
                TEMP=(SU(I)-XSUM)/S(I)
            ELSE
                TEMP=(SL(I)-XSUM)/S(I)
            END IF
            IF (TEMP  <  STPLEN) THEN
                STPLEN=TEMP
                IACT=I
            END IF
        END IF
70  CONTINUE
    !
    !     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in
    !
    SDEC=ZERO
    IF (STPLEN  >  ZERO) THEN
        ITERC=ITERC+1
        TEMP=SHS/STEPSQ
        IF (IACT .EQ. 0 .AND. TEMP  >  ZERO) THEN
            CRVMIN=DMIN1(CRVMIN,TEMP)
            IF (CRVMIN .EQ. ONEMIN) CRVMIN=TEMP
        END IF
        GGSAV=GREDSQ
        GREDSQ=ZERO
        DO I=1,N
            GNEW(I)=GNEW(I)+STPLEN*HS(I)
            IF (XBDI(I) .EQ. ZERO) GREDSQ=GREDSQ+GNEW(I)**2
            D(I)=D(I)+STPLEN*S(I)
        END DO
        SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
        QRED=QRED+SDEC
    END IF
    !
    !     Restart the conjugate gradient method if it has hit a new bound.
    !
    IF (IACT  >  0) THEN
        NACT=NACT+1
        XBDI(IACT)=ONE
        IF (S(IACT)  <  ZERO) XBDI(IACT)=ONEMIN
        DELSQ=DELSQ-D(IACT)**2
        IF (DELSQ  <=  ZERO) GOTO 90
        GOTO 20
    END IF
    !
    !     If STPLEN is less than BLEN, then either apply another conjugate
    !     gradient iteration or RETURN.
    !
    IF (STPLEN  <  BLEN) THEN
        IF (ITERC .EQ. ITERMAX) GOTO 190
        IF (SDEC  <=  0.01D0*QRED) GOTO 190
        BETA=GREDSQ/GGSAV
        GOTO 30
    END IF
90  CRVMIN=ZERO
    !
    !  Prepare for the alternative iteration by calculating some scalars
    !  and by multiplying the reduced D by the second derivative matrix of
    !  Q, where S holds the reduced D in the call of GGMULT.
    !
100 IF (NACT  >=  N-1) GOTO 190
    DREDSQ=ZERO
    DREDG=ZERO
    GREDSQ=ZERO
    DO I=1,N
        IF (XBDI(I) .EQ. ZERO) THEN
            DREDSQ=DREDSQ+D(I)**2
            DREDG=DREDG+D(I)*GNEW(I)
            GREDSQ=GREDSQ+GNEW(I)**2
            S(I)=D(I)
        ELSE
            S(I)=ZERO
        END IF
    END DO
    ITCSAV=ITERC
    GOTO 210
    !
    !     Let the search direction S be a linear combination of the reduced
    !     and the reduced G that is orthogonal to the reduced D.
    !
120 ITERC=ITERC+1
    TEMP=GREDSQ*DREDSQ-DREDG*DREDG
    IF (TEMP  <=  1.0D-4*QRED*QRED) GOTO 190
    TEMP=DSQRT(TEMP)
    DO I=1,N
        IF (XBDI(I) .EQ. ZERO) THEN
            S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
        ELSE
            S(I)=ZERO
        END IF
    END DO
    SREDG=-TEMP
    !
    !  By considering the simple bounds on the variables, calculate an upper
    !  bound on the tangent of half the angle of the alternative iteration,
    !  namely ANGBD, except that, if already a free variable has reached a
    !  bound, there is a branch back to label 100 after fixing that variable.
    !
    ANGBD=ONE
    IACT=0
    DO I=1,N
        IF (XBDI(I) .EQ. ZERO) THEN
            TEMPA=XOPT(I)+D(I)-SL(I)
            TEMPB=SU(I)-XOPT(I)-D(I)
            IF (TEMPA  <=  ZERO) THEN
                NACT=NACT+1
                XBDI(I)=ONEMIN
                GOTO 100
            ELSE IF (TEMPB  <=  ZERO) THEN
                NACT=NACT+1
                XBDI(I)=ONE
                GOTO 100
            END IF
            RATIO=ONE
            SSQ=D(I)**2+S(I)**2
            TEMP=SSQ-(XOPT(I)-SL(I))**2
            IF (TEMP  >  ZERO) THEN
                TEMP=DSQRT(TEMP)-S(I)
                IF (ANGBD*TEMP  >  TEMPA) THEN
                    ANGBD=TEMPA/TEMP
                    IACT=I
                    XSAV=ONEMIN
                END IF
            END IF
            TEMP=SSQ-(SU(I)-XOPT(I))**2
            IF (TEMP  >  ZERO) THEN
                TEMP=DSQRT(TEMP)+S(I)
                IF (ANGBD*TEMP  >  TEMPB) THEN
                    ANGBD=TEMPB/TEMP
                    IACT=I
                    XSAV=ONE
                END IF
            END IF
        END IF
    END DO
    !
    !     Calculate HHD and some curvatures for the alternative iteration.
    !
    GOTO 210
150 SHS=ZERO
    DHS=ZERO
    DHD=ZERO
    DO I=1,N
        IF (XBDI(I) .EQ. ZERO) THEN
            SHS=SHS+S(I)*HS(I)
            DHS=DHS+D(I)*HS(I)
            DHD=DHD+D(I)*HRED(I)
        END IF
    END DO
    !
    !  Seek the greatest reduction in Q for a range of equally spaced values
    !  of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
    !  the alternative iteration.
    !
    REDMAX=ZERO
    ISAV=0
    REDSAV=ZERO
    IU=17.0D0*ANGBD+3.1D0
    DO I=1,IU
        ANGT=ANGBD*real(I,dp)/IU
        STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
        TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
        REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
        IF (REDNEW  >  REDMAX) THEN
            REDMAX=REDNEW
            ISAV=I
            RDPREV=REDSAV
        ELSE IF (I .EQ. ISAV+1) THEN
            RDNEXT=REDNEW
        END IF
        REDSAV=REDNEW
    END DO
    !
    !  Return if the reduction is zero. Otherwise, set the sine and cosine
    !  of the angle of the alternative iteration, and calculate SDEC.
    !
    IF (ISAV .EQ. 0) GOTO 190
    IF (ISAV  <  IU) THEN
        TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
        ANGT=ANGBD*(ISAV+HALF*TEMP)/IU
    END IF
    CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
    STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
    TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
    SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
    IF (SDEC  <=  ZERO) GOTO 190
    !
    !  Update GNEW, D and HRED. If the angle of the alternative iteration
    !  is restricted by a bound on a free variable, that variable is fixed
    !  at the bound.
    !
    DREDG=ZERO
    GREDSQ=ZERO
    DO I=1,N
        GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
        IF (XBDI(I) .EQ. ZERO) THEN
            D(I)=CTH*D(I)+STH*S(I)
            DREDG=DREDG+D(I)*GNEW(I)
            GREDSQ=GREDSQ+GNEW(I)**2
        END IF
        HRED(I)=CTH*HRED(I)+STH*HS(I)
    END DO
    QRED=QRED+SDEC
    IF (IACT  >  0 .AND. ISAV .EQ. IU) THEN
        NACT=NACT+1
        XBDI(IACT)=XSAV
        GOTO 100
    END IF
    !  If SDEC is sufficiently small, then RETURN after setting XNEW to
    !  XOPT+D, giving careful attention to the bounds.
    !
    IF (SDEC  >  0.01D0*QRED) GOTO 120
190 DSQ=ZERO
    DO I=1,N
        XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
        IF (XBDI(I) .EQ. ONEMIN) XNEW(I)=SL(I)
        IF (XBDI(I) .EQ. ONE) XNEW(I)=SU(I)
        D(I)=XNEW(I)-XOPT(I)
        DSQ=DSQ+D(I)**2
    END DO
    RETURN

    !  The following instructions multiply the current S-vector by the second
    !  derivative matrix of the quadratic model, putting the product in HS.
    !  They are reached from three different parts of the software above and
    !  they can be regarded as an external subroutine.
    !
210 IH=0
    DO J=1,N
        HS(J)=ZERO
        DO I=1,J
            IH=IH+1
            IF (I  <  J) HS(J)=HS(J)+HQ(IH)*S(I)
            HS(I)=HS(I)+HQ(IH)*S(J)
        END DO
    END DO
    DO K=1,NPT
        IF (PQ(K) .NE. ZERO) THEN
            TEMP=ZERO
            DO J=1,N
                TEMP=TEMP+XPT(K,J)*S(J)
            END DO
            TEMP=TEMP*PQ(K)
            DO I=1,N
                HS(I)=HS(I)+TEMP*XPT(K,I)
            END DO
        END IF
    END DO
    IF (CRVMIN .NE. ZERO) GOTO 50
    IF (ITERC  >  ITCSAV) GOTO 150
    DO I=1,N
        HRED(I)=HS(I)
    END DO
    GOTO 120
    END SUBROUTINE TRSBOX


    SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
    REAL(dp) alpha , Beta ,  Denom , one , tau , temp , tempa ,   &
        & tempb , zero , ztest
    INTEGER i , j , jl , jp , k , Knew , N , Ndim , Npt , nptm
    REAL(dp) BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
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
    ONE=1.0D0
    ZERO=0.0D0
    NPTM=NPT-N-1
    ZTEST=ZERO
    DO K=1,NPT
        DO J=1,NPTM
            ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
        END DO
    END DO
    ZTEST=1.0D-20*ZTEST
    !
    !     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    !
    JL=1
    DO J=2,NPTM
        IF (DABS(ZMAT(KNEW,J))  >  ZTEST) THEN
            TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
            TEMPA=ZMAT(KNEW,1)/TEMP
            TEMPB=ZMAT(KNEW,J)/TEMP
            DO I=1,NPT
                TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
                ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
                ZMAT(I,1)=TEMP
            END DO
        END IF
        ZMAT(KNEW,J)=ZERO
    END DO
    !
    !     Put the first NPT components of the KNEW-th column of HLAG into W,
    !     and calculate the parameters of the updating formula.
    !
    DO I=1,NPT
        W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
    END DO
    ALPHA=W(KNEW)
    TAU=VLAG(KNEW)
    VLAG(KNEW)=VLAG(KNEW)-ONE
    !
    !     Complete the updating of ZMAT.
    !
    TEMP=DSQRT(DENOM)
    TEMPB=ZMAT(KNEW,1)/TEMP
    TEMPA=TAU/TEMP
    DO I=1,NPT
        ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
    END DO
    !
    !     Finally, update the matrix BMAT.
    !
    DO J=1,N
        JP=NPT+J
        W(JP)=BMAT(KNEW,J)
        TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
        TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
        DO I=1,JP
            BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
            IF (I  >  NPT) BMAT(JP,I-NPT)=BMAT(I,J)
        END DO
    END DO
    END SUBROUTINE UPDATE



    !*****************************************************************************************
    !>
    !  This subroutine seeks the least value of a function of many variables,
    !  by a trust region method that forms quadratic models by interpolation.
    !  There can be some freedom in the interpolation conditions, which is
    !  taken up by minimizing the Frobenius norm of the change to the second
    !  derivative of the quadratic model, beginning with a zero matrix.

    logical function newuoa (this, obj, funk, n, npt, x, rhobeg, rhoend, iprint, maxfun)
    use iso_c_binding
    class(TNEWUOA) this
    class(*), target :: obj
    real(dp), external :: funk !a class function f(obj,x) where obj is any class instance

    integer,intent(in)                  :: n       !! the number of variables. must be at least 2.
    integer,intent(in)                  :: npt     !! The number of interpolation conditions.
    !! Its value must be in the interval `[N+2,(N+1)(N+2)/2]`.
    real(dp),dimension(*),intent(inout) :: x       !! Initial values of the variables must be set in X(1),X(2),...,X(N). They
    !! will be changed to the values that give the least calculated F.
    real(dp),intent(in)                 :: rhobeg  !! RHOBEG and RHOEND must be set to the initial and final values of a trust
    !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
    !! RHOBEG should be about one tenth of the greatest expected change to a
    !! variable, and RHOEND should indicate the accuracy that is required in
    !! the final values of the variables.
    real(dp),intent(in)                 :: rhoend  !! RHOBEG and RHOEND must be set to the initial and final values of a trust
    !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
    !! RHOBEG should be about one tenth of the greatest expected change to a
    !! variable, and RHOEND should indicate the accuracy that is required in
    !! the final values of the variables.
    integer,intent(in)                  :: iprint  !! The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
    !! amount of printing. Specifically, there is no output if IPRINT=0 and
    !! there is output only at the return if IPRINT=1. Otherwise, each new
    !! value of RHO is printed, with the best vector of variables so far and
    !! the corresponding value of the objective function. Further, each new
    !! value of F with its variables are output if IPRINT=3.
    integer,intent(in)                  :: maxfun  !! an upper bound on the number of calls of CALFUN.
    !! for the variables `X(1),X(2),...,X(N)`.

    real(dp),dimension(:),allocatable :: w
    integer :: np,nptm,ndim,ixb,ixo,ixn,ixp,ifv,igq,ihq,ipq,ibmat,izmat,id,ivl,iw

    ! Partition the working space array, so that different parts of it can be
    ! treated separately by the subroutine that performs the main calculation.

    call C_F_PROCPOINTER(c_funloc(funk), this%funkk)
    this%obj => obj
    this%Last_bestfit = 1d30

    np = n + 1
    nptm = npt - np
    if (npt < n+2 .or. npt > ((n+2)*np)/2) then
        write(*,*) 'Return from NEWUOA because NPT is not in the required interval'
        newuoa = .false.
        return
    end if

    ! The array W will be used for working space
    allocate(w((NPT+13)*(NPT+N)+3*N*(N+3)/2))

    ndim = npt + n
    ixb = 1
    ixo = ixb + n
    ixn = ixo + n
    ixp = ixn + n
    ifv = ixp + n * npt
    igq = ifv + npt
    ihq = igq + n
    ipq = ihq + (n*np) / 2
    ibmat = ipq + npt
    izmat = ibmat + ndim * n
    id = izmat + npt * nptm
    ivl = id + n
    iw = ivl + ndim

    ! The above settings provide a partition of W for subroutine NEWUOB.
    ! The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
    ! W plus the space that is needed by the last array of NEWUOB.

    newuoa = this%newuob (n, npt, x, rhobeg, rhoend, iprint, maxfun, w(ixb), w(ixo), w(ixn), &
        w(ixp), w(ifv), w(igq), w(ihq), w(ipq), w(ibmat), w(izmat), ndim, &
        w(id), w(ivl), w(iw))


    deallocate(w)

    end function newuoa
    !*****************************************************************************************

    logical function newuob (this, n, npt, x, rhobeg, rhoend, iprint, maxfun, xbase, xopt, &
        xnew, xpt, fval, gq, hq, pq, bmat, zmat, ndim, d, vlag, w)
    class(TNEWUOA) this

    real(dp) alpha , beta , bsum , crvmin , delta , detrat , &
        & diff , diffa , diffb , diffc , distsq , dnorm , dsq ,      &
        & dstep , dx , f , fbeg , fopt
    real(dp) fsave , gisq , gqsq , half , hdiag , one , ratio , recip , reciq , &
        rho , Rhobeg , Rhoend , rhosq , sum , suma , sumb
    real(dp) sumz , temp , tempq , tenth , vquad , xjpt , xoptsq , zero, xipt
    INTEGER i , idz , ih , ip , Iprint , ipt , itemp , itest , &
        & j , jp , jpt , k , knew , kopt , ksave , ktemp , Maxfun , &
        & N , Ndim
    INTEGER nf , nfm , nfmm , nfsav , nftest , nh , np , Npt , nptm
    real(dp) x (*), xbase (*), xopt (*), xnew (*), xpt (npt,*), fval (*), gq (*), hq(*), &
        pq (*), bmat (ndim,*), zmat (npt,*), d (*), vlag (*), w (*)
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
    np = n + 1
    nh = (n*np) / 2
    nptm = npt - np
    nftest = max (maxfun, 1)
    !
    !     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    !
    do j = 1, n
        xbase (j) = x (j)
        do k = 1, npt
            xpt (k, j) = zero
        end do
        do i = 1, ndim
            bmat (i, j) = zero
        end do
    end do
    do ih = 1, nh
        hq (ih) = zero
    end do
    do k = 1, npt
        pq (k) = zero
        do j = 1, nptm
            zmat (k, j) = zero
        end do
    end do
    !
    !     Begin the initialization procedure. NF becomes one more than the number
    !     of function values so far. The coordinates of the displacement of the
    !     next initial interpolation point from XBASE are set in XPT(NF,.).
    !
    rhosq = rhobeg * rhobeg
    recip = one / rhosq
    reciq = sqrt (half) / rhosq
    nf = 0
50  nfm = nf
    nfmm = nf - n
    nf = nf + 1
    if (nfm <= 2*n) then
        if (nfm >= 1 .and. nfm <= n) then
            xpt (nf, nfm) = rhobeg
        else if (nfm > n) then
            xpt (nf, nfmm) = - rhobeg
        end if
    else
        itemp = (nfmm-1) / n
        jpt = nfm - itemp * n - n
        ipt = jpt + itemp
        if (ipt > n) then
            itemp = jpt
            jpt = ipt - n
            ipt = itemp
        end if
        xipt = rhobeg
        if (fval(ipt+np) < fval(ipt+1)) xipt = - xipt
        xjpt = rhobeg
        if (fval(jpt+np) < fval(jpt+1)) xjpt = - xjpt
        xpt (nf, ipt) = xipt
        xpt (nf, jpt) = xjpt
    end if
    !
    !     Calculate the next value of F, label 70 being reached immediately
    !     after this calculation. The least function value so far and its index
    !     are required.
    !
    do j = 1, n
        x (j) = xpt (nf, j) + xbase (j)
    end do
    go to 310
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
    if (nfm <= 2*n) then
        if (nfm >= 1 .and. nfm <= n) then
            gq (nfm) = (f-fbeg) / rhobeg
            if (npt < nf+n) then
                bmat (1, nfm) = - one / rhobeg
                bmat (nf, nfm) = one / rhobeg
                bmat (npt+nfm, nfm) = - half * rhosq
            end if
        else if (nfm > n) then
            bmat (nf-n, nfmm) = half / rhobeg
            bmat (nf, nfmm) = - half / rhobeg
            zmat (1, nfmm) = - reciq - reciq
            zmat (nf-n, nfmm) = reciq
            zmat (nf, nfmm) = reciq
            ih = (nfmm*(nfmm+1)) / 2
            temp = (fbeg-f) / rhobeg
            hq (ih) = (gq(nfmm)-temp) / rhobeg
            gq (nfmm) = half * (gq(nfmm)+temp)
        end if
        !
        !     Set the off-diagonal second derivatives of the Lagrange functions and
        !     the initial quadratic model.
        !
    else
        ih = (ipt*(ipt-1)) / 2 + jpt
        if (xipt < zero) ipt = ipt + n
        if (xjpt < zero) jpt = jpt + n
        zmat (1, nfmm) = recip
        zmat (nf, nfmm) = recip
        zmat (ipt+1, nfmm) = - recip
        zmat (jpt+1, nfmm) = - recip
        hq (ih) = (fbeg-fval(ipt+1)-fval(jpt+1)+f) / (xipt*xjpt)
    end if
    if (nf < npt) go to 50
    !
    !     Begin the iterative procedure, because the initial model is complete.
    !
    rho = rhobeg
    delta = rho
    idz = 1
    diffa = zero
    diffb = zero
    itest = 0
    xoptsq = zero
    do i = 1, n
        xopt (i) = xpt (kopt, i)
        xoptsq = xoptsq + xopt (i) ** 2
    end do
90  nfsav = nf
    !
    !     Generate the next trust region step and test its length. Set KNEW
    !     to -1 if the purpose of the next F will be to improve the model.
    !
100 knew = 0
    call trsapp (n, npt, xopt, xpt, gq, hq, pq, delta, d, w, w(np), w(np+n), &
        w(np+2*n), crvmin)
    dsq = zero
    do i = 1, n
        dsq = dsq + d (i) ** 2
    end do
    dnorm = min (delta, sqrt(dsq))
    if (dnorm < half*rho) then
        knew = - 1
        delta = tenth * delta
        ratio = - 1.0_dp
        if (delta <= 1.5_dp*rho) delta = rho
        if (nf <= nfsav+2) go to 460
        temp = 0.125_dp * crvmin * rho * rho
        if (temp <= max(diffa, diffb, diffc)) go to 460
        go to 490
    end if
    !
    !     Shift XBASE if XOPT may be too far from XBASE. First make the changes
    !     to BMAT that do not depend on ZMAT.
    !
120 if (dsq <= 1.0e-3_dp*xoptsq) then
        tempq = 0.25_dp * xoptsq
        do k = 1, npt
            sum = zero
            do i = 1, n
                sum = sum + xpt (k, i) * xopt (i)
            end do
            temp = pq (k) * sum
            sum = sum - half * xoptsq
            w (npt+k) = sum
            do i = 1, n
                gq (i) = gq (i) + temp * xpt (k, i)
                xpt (k, i) = xpt (k, i) - half * xopt (i)
                vlag (i) = bmat (k, i)
                w (i) = sum * xpt (k, i) + tempq * xopt (i)
                ip = npt + i
                do j = 1, i
                    bmat (ip, j) = bmat (ip, j) + vlag (i) * w (j) + w (i) * vlag (j)
                end do
            end do
        end do
        !
        !     Then the revisions of BMAT that depend on ZMAT are calculated.
        !
        do k = 1, nptm
            sumz = zero
            do i = 1, npt
                sumz = sumz + zmat (i, k)
                w (i) = w (npt+i) * zmat (i, k)
            end do
            do j = 1, n
                sum = tempq * sumz * xopt (j)
                do i = 1, npt
                    sum = sum + w (i) * xpt (i, j)
                end do
                vlag (j) = sum
                if (k < idz) sum = - sum
                do i = 1, npt
                    bmat (i, j) = bmat (i, j) + sum * zmat (i, k)
                end do
            end do
            do i = 1, n
                ip = i + npt
                temp = vlag (i)
                if (k < idz) temp = - temp
                do j = 1, i
                    bmat (ip, j) = bmat (ip, j) + temp * vlag (j)
                end do
            end do
        end do
        !
        !     The following instructions complete the shift of XBASE, including
        !     the changes to the parameters of the quadratic model.
        !
        ih = 0
        do j = 1, n
            w (j) = zero
            do k = 1, npt
                w (j) = w (j) + pq (k) * xpt (k, j)
                xpt (k, j) = xpt (k, j) - half * xopt (j)
            end do
            do i = 1, j
                ih = ih + 1
                if (i < j) gq (j) = gq (j) + hq (ih) * xopt (i)
                gq (i) = gq (i) + hq (ih) * xopt (j)
                hq (ih) = hq (ih) + w (i) * xopt (j) + xopt (i) * w (j)
                bmat (npt+i, j) = bmat (npt+j, i)
            end do
        end do
        do j = 1, n
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
        call biglag (n, npt, xopt, xpt, bmat, zmat, idz, ndim, knew, dstep, d, alpha, &
            vlag, vlag(npt+1), w, w(np), w(np+n))
    end if
    !
    !     Calculate VLAG and BETA for the current choice of D. The first NPT
    !     components of W_check will be held in W.
    !
    do k = 1, npt
        suma = zero
        sumb = zero
        sum = zero
        do j = 1, n
            suma = suma + xpt (k, j) * d (j)
            sumb = sumb + xpt (k, j) * xopt (j)
            sum = sum + bmat (k, j) * d (j)
        end do
        w (k) = suma * (half*suma+sumb)
        vlag (k) = sum
    end do
    beta = zero
    do k = 1, nptm
        sum = zero
        do i = 1, npt
            sum = sum + zmat (i, k) * w (i)
        end do
        if (k < idz) then
            beta = beta + sum * sum
            sum = - sum
        else
            beta = beta - sum * sum
        end if
        do i = 1, npt
            vlag (i) = vlag (i) + sum * zmat (i, k)
        end do
    end do
    bsum = zero
    dx = zero
    do j = 1, n
        sum = zero
        do i = 1, npt
            sum = sum + w (i) * bmat (i, j)
        end do
        bsum = bsum + sum * d (j)
        jp = npt + j
        do k = 1, n
            sum = sum + bmat (jp, k) * d (k)
        end do
        vlag (jp) = sum
        bsum = bsum + sum * d (j)
        dx = dx + d (j) * xopt (j)
    end do
    beta = dx * dx + dsq * (xoptsq+dx+dx+half*dsq) + beta - bsum
    vlag (kopt) = vlag (kopt) + one
    !
    !     If KNEW is positive and if the cancellation in DENOM is unacceptable,
    !     then BIGDEN calculates an alternative model step, XNEW being used for
    !     working space.
    !
    if (knew > 0) then
        temp = one + alpha * beta / vlag (knew) ** 2
        if (abs(temp) <= 0.8_dp) then
            call bigden (n, npt, xopt, xpt, bmat, zmat, idz, ndim, kopt, knew, d, w, &
                vlag, beta, xnew, w(ndim+1), w(6*ndim+1))
        end if
    end if
    !
    !     Calculate the next value of the objective function.
    !
290 do i = 1, n
        xnew (i) = xopt (i) + d (i)
        x (i) = xbase (i) + xnew (i)
    end do
    nf = nf + 1
310 if (nf > nftest) then
        nf = nf - 1
        if (iprint > 0) print 320
320     format (/ 4 x, 'Return from NEWUOA because CALFUN has been',&
            ' called MAXFUN times.')
        newuob = .false.
        go to 530
    end if
    f = this%funkk(this%obj, X(1:N))
    if (iprint == 3) then
        print 330, nf, f, (x(i), i=1, n)
330     format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
            '    The corresponding X is:' / (2 x, 5d15.6))
    end if
    if (nf <= npt) go to 70
    if (knew ==-1) go to 530
    !
    !     Use the quadratic model to predict the change in F due to the step D,
    !     and set DIFF to the error of this prediction.
    !
    vquad = zero
    ih = 0
    do j = 1, n
        vquad = vquad + d (j) * gq (j)
        do i = 1, j
            ih = ih + 1
            temp = d (i) * xnew (j) + d (j) * xopt (i)
            if (i == j) temp = half * temp
            vquad = vquad + temp * hq (ih)
        end do
    end do
    do k = 1, npt
        vquad = vquad + pq (k) * w (k)
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
        do i = 1, n
            xopt (i) = xnew (i)
            xoptsq = xoptsq + xopt (i) ** 2
        end do
    end if
    ksave = knew
    if (knew > 0) go to 410
    !
    !     Pick the next value of DELTA after a trust region step.
    !
    if (vquad >= zero) then
        if (iprint > 0) print 370
370     format (/ 4 x, 'Return from NEWUOA because a trust',&
            ' region step has failed to reduce Q.')
        newuob = .false.
        go to 530
    end if
    ratio = (f-fsave) / vquad
    if (ratio <= tenth) then
        delta = half * dnorm
    else if (ratio <= 0.7_dp) then
        delta = max (half*delta, dnorm)
    else
        delta = max (half*delta, dnorm+dnorm)
    end if
    if (delta <= 1.5_dp*rho) delta = rho
    !
    !     Set KNEW to the index of the next interpolation point to be deleted.
    !
    rhosq = max (tenth*delta, rho) ** 2
    ktemp = 0
    detrat = zero
    if (f >= fsave) then
        ktemp = kopt
        detrat = one
    end if
    do k = 1, npt
        hdiag = zero
        do j = 1, nptm
            temp = one
            if (j < idz) temp = - one
            hdiag = hdiag + temp * zmat (k, j) ** 2
        end do
        temp = abs (beta*hdiag+vlag(k)**2)
        distsq = zero
        do j = 1, n
            distsq = distsq + (xpt(k, j)-xopt(j)) ** 2
        end do
        if (distsq > rhosq) temp = temp * (distsq/rhosq) ** 3
        if (temp > detrat .and. k /= ktemp) then
            detrat = temp
            knew = k
        end if
    end do
    if (knew == 0) go to 460
    !
    !     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
    !     can be moved. Begin the updating of the quadratic model, starting
    !     with the explicit second derivative term.
    !
410 call update_newuoa (n, npt, bmat, zmat, idz, ndim, vlag, beta, knew, w)
    fval (knew) = f
    ih = 0
    do i = 1, n
        temp = pq (knew) * xpt (knew, i)
        do j = 1, i
            ih = ih + 1
            hq (ih) = hq (ih) + temp * xpt (knew, j)
        end do
    end do
    pq (knew) = zero
    !
    !     Update the other second derivative parameters, and then the gradient
    !     vector of the model. Also include the new interpolation point.
    !
    do j = 1, nptm
        temp = diff * zmat (knew, j)
        if (j < idz) temp = - temp
        do k = 1, npt
            pq (k) = pq (k) + temp * zmat (k, j)
        end do
    end do
    gqsq = zero
    do i = 1, n
        gq (i) = gq (i) + diff * bmat (knew, i)
        gqsq = gqsq + gq (i) ** 2
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
            do k = 1, npt
                vlag (k) = fval (k) - fval (kopt)
            end do
            gisq = zero
            do i = 1, n
                sum = zero
                do k = 1, npt
                    sum = sum + bmat (k, i) * vlag (k)
                end do
                gisq = gisq + sum * sum
                w (i) = sum
            end do
            !
            !     Test whether to replace the new quadratic model by the least Frobenius
            !     norm interpolant, making the replacement if the test is satisfied.
            !
            itest = itest + 1
            if (gqsq < 100.0_dp*gisq) itest = 0
            if (itest >= 3) then
                do i = 1, n
                    gq (i) = w (i)
                end do
                do ih = 1, nh
                    hq (ih) = zero
                end do
                do j = 1, nptm
                    w (j) = zero
                    do k = 1, npt
                        w (j) = w (j) + vlag (k) * zmat (k, j)
                    end do
                    if (j < idz) w (j) = - w (j)
                end do
                do k = 1, npt
                    pq (k) = zero
                    do j = 1, nptm
                        pq (k) = pq (k) + zmat (k, j) * w (j)
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
    if (f <= fsave+tenth*vquad) go to 100
    if (ksave > 0) go to 100
    !
    !     Alternatively, find out if the interpolation points are close enough
    !     to the best point so far.
    !
    knew = 0
460 distsq = 4.0_dp * delta * delta
    do k = 1, npt
        sum = zero
        do j = 1, n
            sum = sum + (xpt(k, j)-xopt(j)) ** 2
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
        dsq = dstep * dstep
        go to 120
    end if
    if (ratio > zero) go to 100
    if (max(delta, dnorm) > rho) go to 100
    !
    !     The calculations with the current value of RHO are complete. Pick the
    !     next values of RHO and DELTA.
    !
490 if (rho > rhoend) then
        delta = half * rho
        ratio = rho / rhoend
        if (ratio <= 16.0_dp) then
            rho = rhoend
        else if (ratio <= 250.0_dp) then
            rho = sqrt (ratio) * rhoend
        else
            rho = tenth * rho
        end if
        delta = max (delta, rho)
        this%Last_bestfit = fopt
        if (iprint >= 2) then
            if (iprint >= 3) print 500
500         format (5 x)
            print 510, rho, nf
510         format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
                ' function values =', i6)
            print 520, fopt, (xbase(i)+xopt(i), i=1, n)
520         format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
                'The corresponding X is:'/(2 x, 5d15.6))
        end if
        go to 90
    end if
    !
    !     Return from the calculation, after another Newton-Raphson step, if
    !     it is too short to have been tried before.
    !
    if (knew ==-1) go to 290
530 if (fopt <= f) then
        do i = 1, n
            x (i) = xbase (i) + xopt (i)
        end do
        f = fopt
    end if
    if (iprint >= 1) then
        print 550, nf
550     format (/ 4 x, 'At the return from NEWUOA', 5 x,&
            'Number of function values =', i6)
        print 520, f, (x(i), i=1, n)
    end if

    end function newuob

    subroutine bigden (n, npt, xopt, xpt, bmat, zmat, idz, ndim, kopt, knew, d, w, vlag, &
        beta, s, wvec, prod)
    real(dp) alpha , angle , Beta , dd , denmax , denold , densav , diff , ds , dstemp , &
        dtest , half , one , quart
    real(dp) ss , ssden , sstemp , step , sum , sumold , tau ,      &
        & temp , tempa , tempb , tempc , two , twopi , xoptd , xopts
    real(dp) xoptsq , zero
    INTEGER i , Idz , ip , isave , iterc , iu , j , jc , k , Knew ,   &
        & Kopt , ksav , N , Ndim , Npt , nptm , nw
    real(dp) xopt (*), xpt (npt,*), bmat (ndim,*), zmat (npt,*), d (*), w (*), vlag &
        (*), s (*), wvec (ndim,*), prod (ndim,*)
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
    twopi = 8.0_dp * atan (one)
    nptm = npt - n - 1
    !
    !     Store the first NPT elements of the KNEW-th column of H in W(N+1)
    !     to W(N+NPT).
    !
    do k = 1, npt
        w (n+k) = zero
    end do
    do j = 1, nptm
        temp = zmat (knew, j)
        if (j < idz) temp = - temp
        do k = 1, npt
            w (n+k) = w (n+k) + temp * zmat (k, j)
        end do
    end do
    alpha = w (n+knew)
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
    do i = 1, n
        dd = dd + d (i) ** 2
        s (i) = xpt (knew, i) - xopt (i)
        ds = ds + d (i) * s (i)
        ss = ss + s (i) ** 2
        xoptsq = xoptsq + xopt (i) ** 2
    end do
    if (ds*ds > 0.99_dp*dd*ss) then
        ksav = knew
        dtest = ds * ds / ss
        do k = 1, npt
            if (k /= kopt) then
                dstemp = zero
                sstemp = zero
                do i = 1, n
                    diff = xpt (k, i) - xopt (i)
                    dstemp = dstemp + d (i) * diff
                    sstemp = sstemp + diff * diff
                end do
                if (dstemp*dstemp/sstemp < dtest) then
                    ksav = k
                    dtest = dstemp * dstemp / sstemp
                    ds = dstemp
                    ss = sstemp
                end if
            end if
        end do
        do i = 1, n
            s (i) = xpt (ksav, i) - xopt (i)
        end do
    end if
    ssden = dd * ss - ds * ds
    iterc = 0
    densav = zero
    !
    !     Begin the iteration by overwriting S with a vector that has the
    !     required length and direction.
    !
70  iterc = iterc + 1
    temp = one / sqrt (ssden)
    xoptd = zero
    xopts = zero
    do i = 1, n
        s (i) = temp * (dd*s(i)-ds*d(i))
        xoptd = xoptd + xopt (i) * d (i)
        xopts = xopts + xopt (i) * s (i)
    end do
    !
    !     Set the coefficients of the first two terms of BETA.
    !
    tempa = half * xoptd * xoptd
    tempb = half * xopts * xopts
    den (1) = dd * (xoptsq+half*dd) + tempa + tempb
    den (2) = two * xoptd * dd
    den (3) = two * xopts * dd
    den (4) = tempa - tempb
    den (5) = xoptd * xopts
    do i = 6, 9
        den (i) = zero
    end do
    !
    !     Put the coefficients of Wcheck in WVEC.
    !
    do k = 1, npt
        tempa = zero
        tempb = zero
        tempc = zero
        do i = 1, n
            tempa = tempa + xpt (k, i) * d (i)
            tempb = tempb + xpt (k, i) * s (i)
            tempc = tempc + xpt (k, i) * xopt (i)
        end do
        wvec (k, 1) = quart * (tempa*tempa+tempb*tempb)
        wvec (k, 2) = tempa * tempc
        wvec (k, 3) = tempb * tempc
        wvec (k, 4) = quart * (tempa*tempa-tempb*tempb)
        wvec (k, 5) = half * tempa * tempb
    end do
    do i = 1, n
        ip = i + npt
        wvec (ip, 1) = zero
        wvec (ip, 2) = d (i)
        wvec (ip, 3) = s (i)
        wvec (ip, 4) = zero
        wvec (ip, 5) = zero
    end do
    !
    !     Put the coefficents of THETA*Wcheck in PROD.
    !
    do jc = 1, 5
        nw = npt
        if (jc == 2 .or. jc == 3) nw = ndim
        do k = 1, npt
            prod (k, jc) = zero
        end do
        do j = 1, nptm
            sum = zero
            do k = 1, npt
                sum = sum + zmat (k, j) * wvec (k, jc)
            end do
            if (j < idz) sum = - sum
            do k = 1, npt
                prod (k, jc) = prod (k, jc) + sum * zmat (k, j)
            end do
        end do
        if (nw == ndim) then
            do k = 1, npt
                sum = zero
                do j = 1, n
                    sum = sum + bmat (k, j) * wvec (npt+j, jc)
                end do
                prod (k, jc) = prod (k, jc) + sum
            end do
        end if
        do j = 1, n
            sum = zero
            do i = 1, nw
                sum = sum + bmat (i, j) * wvec (i, jc)
            end do
            prod (npt+j, jc) = sum
        end do
    end do
    !
    !     Include in DEN the part of BETA that depends on THETA.
    !
    do k = 1, ndim
        sum = zero
        do i = 1, 5
            par (i) = half * prod (k, i) * wvec (k, i)
            sum = sum + par (i)
        end do
        den (1) = den (1) - par (1) - sum
        tempa = prod (k, 1) * wvec (k, 2) + prod (k, 2) * wvec (k, 1)
        tempb = prod (k, 2) * wvec (k, 4) + prod (k, 4) * wvec (k, 2)
        tempc = prod (k, 3) * wvec (k, 5) + prod (k, 5) * wvec (k, 3)
        den (2) = den (2) - tempa - half * (tempb+tempc)
        den (6) = den (6) - half * (tempb-tempc)
        tempa = prod (k, 1) * wvec (k, 3) + prod (k, 3) * wvec (k, 1)
        tempb = prod (k, 2) * wvec (k, 5) + prod (k, 5) * wvec (k, 2)
        tempc = prod (k, 3) * wvec (k, 4) + prod (k, 4) * wvec (k, 3)
        den (3) = den (3) - tempa - half * (tempb-tempc)
        den (7) = den (7) - half * (tempb+tempc)
        tempa = prod (k, 1) * wvec (k, 4) + prod (k, 4) * wvec (k, 1)
        den (4) = den (4) - tempa - par (2) + par (3)
        tempa = prod (k, 1) * wvec (k, 5) + prod (k, 5) * wvec (k, 1)
        tempb = prod (k, 2) * wvec (k, 3) + prod (k, 3) * wvec (k, 2)
        den (5) = den (5) - tempa - half * tempb
        den (8) = den (8) - par (4) + par (5)
        tempa = prod (k, 4) * wvec (k, 5) + prod (k, 5) * wvec (k, 4)
        den (9) = den (9) - half * tempa
    end do
    !
    !     Extend DEN so that it holds all the coefficients of DENOM.
    !
    sum = zero
    do i = 1, 5
        par (i) = half * prod (knew, i) ** 2
        sum = sum + par (i)
    end do
    denex (1) = alpha * den (1) + par (1) + sum
    tempa = two * prod (knew, 1) * prod (knew, 2)
    tempb = prod (knew, 2) * prod (knew, 4)
    tempc = prod (knew, 3) * prod (knew, 5)
    denex (2) = alpha * den (2) + tempa + tempb + tempc
    denex (6) = alpha * den (6) + tempb - tempc
    tempa = two * prod (knew, 1) * prod (knew, 3)
    tempb = prod (knew, 2) * prod (knew, 5)
    tempc = prod (knew, 3) * prod (knew, 4)
    denex (3) = alpha * den (3) + tempa + tempb - tempc
    denex (7) = alpha * den (7) + tempb + tempc
    tempa = two * prod (knew, 1) * prod (knew, 4)
    denex (4) = alpha * den (4) + tempa + par (2) - par (3)
    tempa = two * prod (knew, 1) * prod (knew, 5)
    denex (5) = alpha * den (5) + tempa + prod (knew, 2) * prod (knew, 3)
    denex (8) = alpha * den (8) + par (4) - par (5)
    denex (9) = alpha * den (9) + prod (knew, 4) * prod (knew, 5)
    !
    !     Seek the value of the angle that maximizes the modulus of DENOM.
    !
    sum = denex (1) + denex (2) + denex (4) + denex (6) + denex (8)
    denold = sum
    denmax = sum
    isave = 0
    iu = 49
    temp = twopi / real (iu+1, dp)
    par (1) = one
    do i = 1, iu
        angle = real (i, dp) * temp
        par (2) = cos (angle)
        par (3) = sin (angle)
        do j = 4, 8, 2
            par (j) = par (2) * par (j-2) - par (3) * par (j-1)
            par (j+1) = par (2) * par (j-1) + par (3) * par (j-2)
        end do
        sumold = sum
        sum = zero
        do j = 1, 9
            sum = sum + denex (j) * par (j)
        end do
        if (abs(sum) > abs(denmax)) then
            denmax = sum
            isave = i
            tempa = sumold
        else if (i == isave+1) then
            tempb = sum
        end if
    end do
    if (isave == 0) tempa = sum
    if (isave == iu) tempb = denold
    step = zero
    if (tempa /= tempb) then
        tempa = tempa - denmax
        tempb = tempb - denmax
        step = half * (tempa-tempb) / (tempa+tempb)
    end if
    angle = temp * (real(isave, dp)+step)
    !
    !     Calculate the new parameters of the denominator, the new VLAG vector
    !     and the new D. Then test for convergence.
    !
    par (2) = cos (angle)
    par (3) = sin (angle)
    do j = 4, 8, 2
        par (j) = par (2) * par (j-2) - par (3) * par (j-1)
        par (j+1) = par (2) * par (j-1) + par (3) * par (j-2)
    end do
    beta = zero
    denmax = zero
    do j = 1, 9
        beta = beta + den (j) * par (j)
        denmax = denmax + denex (j) * par (j)
    end do
    do k = 1, ndim
        vlag (k) = zero
        do j = 1, 5
            vlag (k) = vlag (k) + prod (k, j) * par (j)
        end do
    end do
    tau = vlag (knew)
    dd = zero
    tempa = zero
    tempb = zero
    do i = 1, n
        d (i) = par (2) * d (i) + par (3) * s (i)
        w (i) = xopt (i) + d (i)
        dd = dd + d (i) ** 2
        tempa = tempa + d (i) * w (i)
        tempb = tempb + w (i) * w (i)
    end do
    if (iterc >= n) go to 340
    if (iterc > 1) densav = max (densav, denold)
    if (abs(denmax) <= 1.1_dp*abs(densav)) go to 340
    densav = denmax
    !
    !     Set S to half the gradient of the denominator with respect to D.
    !     Then branch for the next iteration.
    !
    do i = 1, n
        temp = tempa * xopt (i) + tempb * d (i) - vlag (npt+i)
        s (i) = tau * bmat (knew, i) + alpha * temp
    end do
    do k = 1, npt
        sum = zero
        do j = 1, n
            sum = sum + xpt (k, j) * w (j)
        end do
        temp = (tau*w(n+k)-alpha*vlag(k)) * sum
        do i = 1, n
            s (i) = s (i) + temp * xpt (k, i)
        end do
    end do
    ss = zero
    ds = zero
    do i = 1, n
        ss = ss + s (i) ** 2
        ds = ds + d (i) * s (i)
    end do
    ssden = dd * ss - ds * ds
    if (ssden >= 1.0e-8_dp*dd*ss) go to 70
    !
    !     Set the vector W before the RETURN from the subroutine.
    !
340 do k = 1, ndim
        w (k) = zero
        do j = 1, 5
            w (k) = w (k) + wvec (k, j) * par (j)
        end do
    end do
    vlag (kopt) = vlag (kopt) + one

    end subroutine bigden

    subroutine biglag (n, npt, xopt, xpt, bmat, zmat, idz, ndim, knew, delta, d, alpha, &
        hcol, gc, gd, s, w)

    real(dp) Alpha , angle , cf1 , cf2 , cf3 , cf4 , cf5 , cth , &
        & dd , delsq , Delta , denom , dhd , gg , half
    INTEGER i , Idz , isave , iterc , iu , j , k , Knew , N , Ndim ,  &
        & Npt , nptm
    real(dp) one , scale , sp , ss , step , sth , sum , tau ,       &
        & taubeg , taumax , tauold , temp , tempa , tempb , twopi ,  zero
    real(dp) xopt (*), xpt (npt,*), bmat (ndim,*), zmat (npt,*), d (*), hcol (*), &
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
    twopi = 8.0_dp * atan (one)
    delsq = delta * delta
    nptm = npt - n - 1
    !
    !     Set the first NPT components of HCOL to the leading elements of the
    !     KNEW-th column of H.
    !
    iterc = 0
    do k = 1, npt
        hcol (k) = zero
    end do
    do j = 1, nptm
        temp = zmat (knew, j)
        if (j < idz) temp = - temp
        do k = 1, npt
            hcol (k) = hcol (k) + temp * zmat (k, j)
        end do
    end do
    alpha = hcol (knew)
    !
    !     Set the unscaled initial direction D. Form the gradient of LFUNC at
    !     XOPT, and multiply D by the second derivative matrix of LFUNC.
    !
    dd = zero
    do i = 1, n
        d (i) = xpt (knew, i) - xopt (i)
        gc (i) = bmat (knew, i)
        gd (i) = zero
        dd = dd + d (i) ** 2
    end do
    do k = 1, npt
        temp = zero
        sum = zero
        do j = 1, n
            temp = temp + xpt (k, j) * xopt (j)
            sum = sum + xpt (k, j) * d (j)
        end do
        temp = hcol (k) * temp
        sum = hcol (k) * sum
        do i = 1, n
            gc (i) = gc (i) + temp * xpt (k, i)
            gd (i) = gd (i) + sum * xpt (k, i)
        end do
    end do
    !
    !     Scale D and GD, with a sign change if required. Set S to another
    !     vector in the initial two dimensional subspace.
    !
    gg = zero
    sp = zero
    dhd = zero
    do i = 1, n
        gg = gg + gc (i) ** 2
        sp = sp + d (i) * gc (i)
        dhd = dhd + d (i) * gd (i)
    end do
    scale = delta / sqrt (dd)
    if (sp*dhd < zero) scale = - scale
    temp = zero
    if (sp*sp > 0.99_dp*dd*gg) temp = one
    tau = scale * (abs(sp)+half*scale*abs(dhd))
    if (gg*delsq < 0.01_dp*tau*tau) temp = one
    do i = 1, n
        d (i) = scale * d (i)
        gd (i) = scale * gd (i)
        s (i) = gc (i) + temp * gd (i)
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
    do i = 1, n
        dd = dd + d (i) ** 2
        sp = sp + d (i) * s (i)
        ss = ss + s (i) ** 2
    end do
    temp = dd * ss - sp * sp
    if (temp <= 1.0e-8_dp*dd*ss) go to 160
    denom = sqrt (temp)
    do i = 1, n
        s (i) = (dd*s(i)-sp*d(i)) / denom
        w (i) = zero
    end do
    !
    !     Calculate the coefficients of the objective function on the circle,
    !     beginning with the multiplication of S by the second derivative matrix.
    !
    do k = 1, npt
        sum = zero
        do j = 1, n
            sum = sum + xpt (k, j) * s (j)
        end do
        sum = hcol (k) * sum
        do i = 1, n
            w (i) = w (i) + sum * xpt (k, i)
        end do
    end do
    cf1 = zero
    cf2 = zero
    cf3 = zero
    cf4 = zero
    cf5 = zero
    do i = 1, n
        cf1 = cf1 + s (i) * w (i)
        cf2 = cf2 + d (i) * gc (i)
        cf3 = cf3 + s (i) * gc (i)
        cf4 = cf4 + d (i) * gd (i)
        cf5 = cf5 + s (i) * gd (i)
    end do
    cf1 = half * cf1
    cf4 = half * cf4 - cf1
    !
    !     Seek the value of the angle that maximizes the modulus of TAU.
    !
    taubeg = cf1 + cf2 + cf4
    taumax = taubeg
    tauold = taubeg
    isave = 0
    iu = 49
    temp = twopi / real (iu+1, dp)
    do i = 1, iu
        angle = real (i, dp) * temp
        cth = cos (angle)
        sth = sin (angle)
        tau = cf1 + (cf2+cf4*cth) * cth + (cf3+cf5*cth) * sth
        if (abs(tau) > abs(taumax)) then
            taumax = tau
            isave = i
            tempa = tauold
        else if (i == isave+1) then
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
        step = half * (tempa-tempb) / (tempa+tempb)
    end if
    angle = temp * (real(isave, dp)+step)
    !
    !     Calculate the new D and GD. Then test for convergence.
    !
    cth = cos (angle)
    sth = sin (angle)
    tau = cf1 + (cf2+cf4*cth) * cth + (cf3+cf5*cth) * sth
    do i = 1, n
        d (i) = cth * d (i) + sth * s (i)
        gd (i) = cth * gd (i) + sth * w (i)
        s (i) = gc (i) + gd (i)
    end do
    if (abs(tau) <= 1.1_dp*abs(taubeg)) go to 160
    if (iterc < n) go to 80
160 return
    end subroutine biglag

    subroutine trsapp (n, npt, xopt, xpt, gq, hq, pq, delta, step, d, g, hd, hs, crvmin)
    real(dp) alpha , angle , angtest , bstep , cf , Crvmin , cth , &
        & dd , delsq , Delta , dg , dhd , dhs , ds , gg , ggbeg ,&
        & ggsav
    real(dp) half , qadd , qbeg , qmin , qnew , qred , qsav , ratio , reduc , &
        sg , sgk , shs , ss , sth , temp, tempa , tempb , twopi , zero
    INTEGER i , ih , isave , iterc , itermax , itersw , iu , j , k ,  &
        & N , Npt
    real(dp) xopt (*), xpt (npt,*), gq (*), hq (*), pq (*), step (*), d (*), g (*), &
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
    twopi = 8.0_dp * atan (1.0_dp)
    delsq = delta * delta
    iterc = 0
    itermax = n
    itersw = itermax
    do i = 1, n
        d (i) = xopt (i)
    end do
    go to 170
    !
    !     Prepare for the first line search.
    !
20  qred = zero
    dd = zero
    do i = 1, n
        step (i) = zero
        hs (i) = zero
        g (i) = gq (i) + hd (i)
        d (i) = - g (i)
        dd = dd + d (i) ** 2
    end do
    crvmin = zero
    if (dd == zero) go to 160
    ds = zero
    ss = zero
    gg = dd
    ggbeg = gg
    !
    !     Calculate the step to the trust region boundary and the product HD.
    !
40  iterc = iterc + 1
    temp = delsq - ss
    bstep = temp / (ds+sqrt(ds*ds+dd*temp))
    go to 170
50  dhd = zero
    do j = 1, n
        dhd = dhd + d (j) * hd (j)
    end do
    !
    !     Update CRVMIN and set the step-length ALPHA.
    !
    alpha = bstep
    if (dhd > zero) then
        temp = dhd / dd
        if (iterc == 1) crvmin = temp
        crvmin = min (crvmin, temp)
        alpha = min (alpha, gg/dhd)
    end if
    qadd = alpha * (gg-half*alpha*dhd)
    qred = qred + qadd
    !
    !     Update STEP and HS.
    !
    ggsav = gg
    gg = zero
    do i = 1, n
        step (i) = step (i) + alpha * d (i)
        hs (i) = hs (i) + alpha * hd (i)
        gg = gg + (g(i)+hs(i)) ** 2
    end do
    !
    !     Begin another conjugate direction iteration if required.
    !
    if (alpha < bstep) then
        if (qadd <= 0.01_dp*qred) go to 160
        if (gg <= 1.0e-4_dp*ggbeg) go to 160
        if (iterc == itermax) go to 160
        temp = gg / ggsav
        dd = zero
        ds = zero
        ss = zero
        do i = 1, n
            d (i) = temp * d (i) - g (i) - hs (i)
            dd = dd + d (i) ** 2
            ds = ds + d (i) * step (i)
            ss = ss + step (i) ** 2
        end do
        if (ds <= zero) go to 160
        if (ss < delsq) go to 40
    end if
    crvmin = zero
    itersw = iterc
    !
    !     Test whether an alternative iteration is required.
    !
90  if (gg <= 1.0e-4_dp*ggbeg) go to 160
    sg = zero
    shs = zero
    do i = 1, n
        sg = sg + step (i) * g (i)
        shs = shs + step (i) * hs (i)
    end do
    sgk = sg + shs
    angtest = sgk / sqrt (gg*delsq)
    if (angtest <=-0.99_dp) go to 160
    !
    !     Begin the alternative iteration by calculating D and HD and some
    !     scalar products.
    !
    iterc = iterc + 1
    temp = sqrt (delsq*gg-sgk*sgk)
    tempa = delsq / temp
    tempb = sgk / temp
    do i = 1, n
        d (i) = tempa * (g(i)+hs(i)) - tempb * step (i)
    end do
    go to 170
120 dg = zero
    dhd = zero
    dhs = zero
    do i = 1, n
        dg = dg + d (i) * g (i)
        dhd = dhd + hd (i) * d (i)
        dhs = dhs + hd (i) * step (i)
    end do
    !
    !     Seek the value of the angle that minimizes Q.
    !
    cf = half * (shs-dhd)
    qbeg = sg + cf
    qsav = qbeg
    qmin = qbeg
    isave = 0
    iu = 49
    temp = twopi / real (iu+1, dp)
    do i = 1, iu
        angle = real (i, dp) * temp
        cth = cos (angle)
        sth = sin (angle)
        qnew = (sg+cf*cth) * cth + (dg+dhs*cth) * sth
        if (qnew < qmin) then
            qmin = qnew
            isave = i
            tempa = qsav
        else if (i == isave+1) then
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
        angle = half * (tempa-tempb) / (tempa+tempb)
    end if
    angle = temp * (real(isave, dp)+angle)
    !
    !     Calculate the new STEP and HS. Then test for convergence.
    !
    cth = cos (angle)
    sth = sin (angle)
    reduc = qbeg - (sg+cf*cth) * cth - (dg+dhs*cth) * sth
    gg = zero
    do i = 1, n
        step (i) = cth * step (i) + sth * d (i)
        hs (i) = cth * hs (i) + sth * hd (i)
        gg = gg + (g(i)+hs(i)) ** 2
    end do
    qred = qred + reduc
    ratio = reduc / qred
    if (iterc < itermax .and. ratio > 0.01_dp) go to 90
160 return
    !
    !     The following instructions act as a subroutine for setting the vector
    !     HD to the vector D multiplied by the second derivative matrix of Q.
    !     They are called from three different places, which are distinguished
    !     by the value of ITERC.
    !
170 do i = 1, n
        hd (i) = zero
    end do
    do k = 1, npt
        temp = zero
        do j = 1, n
            temp = temp + xpt (k, j) * d (j)
        end do
        temp = temp * pq (k)
        do i = 1, n
            hd (i) = hd (i) + temp * xpt (k, i)
        end do
    end do
    ih = 0
    do j = 1, n
        do i = 1, j
            ih = ih + 1
            if (i < j) hd (j) = hd (j) + hq (ih) * d (i)
            hd (i) = hd (i) + hq (ih) * d (j)
        end do
    end do
    if (iterc == 0) go to 20
    if (iterc <= itersw) go to 50
    go to 120
    end subroutine trsapp

    subroutine update_newuoa (n, npt, bmat, zmat, idz, ndim, vlag, beta, knew, w)
    real(dp) alpha , Beta , denom , one , scala , scalb , tau ,  &
        & tausq , temp , tempa , tempb , zero
    INTEGER i , Idz , iflag , j , ja , jb , jl , jp , Knew , N ,      &
        & Ndim , Npt , nptm
    real(dp) bmat (ndim,*), zmat (npt,*), vlag (*), w (*)
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
    nptm = npt - n - 1
    !
    !     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    !
    jl = 1
    do j = 2, nptm
        if (j == idz) then
            jl = idz
        else if (zmat(knew, j) /= zero) then
            temp = sqrt (zmat(knew, jl)**2+zmat(knew, j)**2)
            tempa = zmat (knew, jl) / temp
            tempb = zmat (knew, j) / temp
            do i = 1, npt
                temp = tempa * zmat (i, jl) + tempb * zmat (i, j)
                zmat (i, j) = tempa * zmat (i, j) - tempb * zmat (i, jl)
                zmat (i, jl) = temp
            end do
            zmat (knew, j) = zero
        end if
    end do
    !
    !     Put the first NPT components of the KNEW-th column of HLAG into W,
    !     and calculate the parameters of the updating formula.
    !
    tempa = zmat (knew, 1)
    if (idz >= 2) tempa = - tempa
    if (jl > 1) tempb = zmat (knew, jl)
    do i = 1, npt
        w (i) = tempa * zmat (i, 1)
        if (jl > 1) w (i) = w (i) + tempb * zmat (i, jl)
    end do
    alpha = w (knew)
    tau = vlag (knew)
    tausq = tau * tau
    denom = alpha * beta + tausq
    vlag (knew) = vlag (knew) - one
    !
    !     Complete the updating of ZMAT when there is only one nonzero element
    !     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
    !     then the first column of ZMAT will be exchanged with another one later.
    !
    iflag = 0
    if (jl == 1) then
        temp = sqrt (abs(denom))
        tempb = tempa / temp
        tempa = tau / temp
        do i = 1, npt
            zmat (i, 1) = tempa * zmat (i, 1) - tempb * vlag (i)
        end do
        if (idz == 1 .and. temp < zero) idz = 2
        if (idz >= 2 .and. temp >= zero) iflag = 1
    else
        !
        !     Complete the updating of ZMAT in the alternative case.
        !
        ja = 1
        if (beta >= zero) ja = jl
        jb = jl + 1 - ja
        temp = zmat (knew, jb) / denom
        tempa = temp * beta
        tempb = temp * tau
        temp = zmat (knew, ja)
        scala = one / sqrt (abs(beta)*temp*temp+tausq)
        scalb = scala * sqrt (abs(denom))
        do i = 1, npt
            zmat (i, ja) = scala * (tau*zmat(i, ja)-temp*vlag(i))
            zmat (i, jb) = scalb * (zmat(i, jb)-tempa*w(i)-tempb*vlag(i))
        end do
        if (denom <= zero) then
            if (beta < zero) idz = idz + 1
            if (beta >= zero) iflag = 1
        end if
    end if
    !
    !     IDZ is reduced in the following case, and usually the first column
    !     of ZMAT is exchanged with a later one.
    !
    if (iflag == 1) then
        idz = idz - 1
        do i = 1, npt
            temp = zmat (i, 1)
            zmat (i, 1) = zmat (i, idz)
            zmat (i, idz) = temp
        end do
    end if
    !
    !     Finally, update the matrix BMAT.
    !
    do j = 1, n
        jp = npt + j
        w (jp) = bmat (knew, j)
        tempa = (alpha*vlag(jp)-tau*w(jp)) / denom
        tempb = (-beta*w(jp)-tau*vlag(jp)) / denom
        do i = 1, jp
            bmat (i, j) = bmat (i, j) + tempa * vlag (i) + tempb * w (i)
            if (i > npt) bmat (jp, i-npt) = bmat (i, j)
        end do
    end do

    end subroutine update_newuoa

    END MODULE Powell
