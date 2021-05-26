
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS MODULE CONTAINS THE FUNCTIONS AND SUBROUTINES REQUIRED TO SOLVE THE FIRM'S
! PROBLEM, COMPUTE THE INVARIANT DISTRIBUTION AND COMPUTE THE MOMENTS USED
! IN CALIBRATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE FIRMSPROBLEM
USE params
USE omp_lib
IMPLICIT NONE
REAL(rp) , ALLOCATABLE :: V0(:,:)
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine guides the program in solving the firm's problem.
SUBROUTINE SOLVEPROBLEM( )
  IMPLICIT NONE
  CALL ALLOCATESTATE( )
  CALL SOLVEFIRMSPROBLEM( )
  CALL SOLVEFIRMSDIST( )
  CALL AGGREGATE( )
  CALL SETMOMENTS( )
  RETURN
END SUBROUTINE SOLVEPROBLEM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine uses a Newton-based method to find the equilibrium wage rate.
SUBROUTINE FIND_EQUIL_WRATE(SHW)
  USE toolkit , ONLY : LMMIN
  IMPLICIT NONE
  INTEGER  , INTENT(IN) :: SHW
  REAL(rp)              :: WAGE0(1),WAGE1(1),FOCS(1)
  INTEGER               :: ITERWRATE,ICODE
  WAGE0(1) = one
  CALL LMMIN(FOCWRATE,WAGE1,FOCS,ITERWRATE,ICODE,WAGE0,IPRINT=MIN(2,SHW),SHCK=DBLE(0.01))
  wrate = WAGE1(1)
  CALL SOLVEPROBLEM( )
  RETURN
END SUBROUTINE FIND_EQUIL_WRATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This function returns the residual in the households' first-order condition.
! The subroutine FIND_EQUIL_wrate find the minimum squared-error of this function.
FUNCTION FOCWRATE(WAGE) RESULT(FOC)
  IMPLICIT NONE
  REAL(rp)               :: WAGE(:),CONS
  REAL(rp) , ALLOCATABLE :: FOC(:)
  ALLOCATE(FOC(SIZE(WAGE)))
  wrate  = WAGE(1) ; CALL SOLVEPROBLEM( )
  CONS   = TOTY - TOTLI
  FOC(1) = ( wrate - theta*((CONS**psic)*(TOTL**psi)))*cien
  RETURN
END FUNCTION FOCWRATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine iterates over the value functions until a fixed point is found.
SUBROUTINE SOLVEFIRMSPROBLEM( )
  IMPLICIT NONE
  REAL(rp) :: test
  INTEGER  :: in,id,iter
  IF (ALLOCATED(V0)) DEALLOCATE(V0) ; ALLOCATE(V0(Nd,Nn))
  V0 = zero ; SKIPPOINT = 0 ; test = 99.99 ; iter = 0
  DO WHILE ( test.gt.DBLE(0.0001) .AND. iter.lt.300 .AND. SKIPPOINT.eq.0)
    iter = iter + 1
    !$OMP PARALLEL NUM_THREADS(threads) DEFAULT(SHARED)
    !$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
    DO id = 1,Nd ; DO in = 1,Nn
      CALL FindOptimalN(id,in)
    END DO ; END DO
    !$OMP END DO
    !$OMP END PARALLEL
    test = 0.000000
    DO id = 1,Nd ; DO in = 1,Nn
      IF ( abs(V(id,in)-V0(id,in)).gt.test ) test = abs(V(id,in) - V0(id,in))
    END DO ; END DO
    V0 = V*DBLE(0.8) + V0*DBLE(0.2)
  END DO
  IF (test.GT.DBLE(0.0001)) PRINT * , '  VALUE FUNCTION DID NOT CONVERGE'
  RETURN
END SUBROUTINE SOLVEFIRMSPROBLEM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given a vector of state variables, this subroutines finds the optimal employment
! for the next period using a Golden Search algorithm.
! It assumes monotonicity of the value fucntion
SUBROUTINE FINDOPTIMALN(id,in)
  IMPLICIT NONE
  INTEGER  , INTENT(IN) :: id,in
  REAL(rp) , PARAMETER  :: alpha = 0.61803399
  REAL(rp)              :: f0,f1,f2,n0,n1,n2,n3,fopt,nopt
  INTEGER               :: iter

  ! If adjustment cost = 0, take the optimal n'
  IF (ABS(fc).LT.tol) THEN
    n0 = MIN(ngrid(Nn),EXP(dgrid(id))*((gamma*bigA/wrate)**(one/(one-gamma))))
    CALL ValueFunc(f0,n0,id,in,1)
    GOTO 1
  END IF

  ! Value of hiring = firing = 0 --> f0
  CALL ValueFunc(f0,ngrid(in),id,in,1)

  ! Optimal labor conditional on hiring
  ! If a point >f0 is found, return (monotonicity)
  IF (in.lt.Nn) THEN
    n0 = ngrid(in) ; n3 = ngrid(Nn)
    n1 = alpha*n0 + (one-alpha)*n3 ; CALL ValueFunc(f1,n1,id,in,0)
    n2 = alpha*n3 + (one-alpha)*n1 ; CALL ValueFunc(f2,n2,id,in,0) ; iter = 0
    DO WHILE ( abs(n0-n3).gt.tol*(abs(n2)+abs(n1)) .and. iter.lt.900 )
      iter = iter + 1
      IF (f2.gt.f1) THEN ; n0 = n1 ; n1 = n2 ; f1 = f2
        n2 = alpha*n1 + (one-alpha)*n3 ; CALL ValueFunc(f2,n2,id,in,0)
      ELSE ; n3 = n2 ; n2 = n1 ; f2 = f1
        n1 = alpha*n2 + (one-alpha)*n0 ; CALL ValueFunc(f1,n1,id,in,0)
      END IF
    END DO
    IF (f1.gt.f2) THEN
      IF (f1.gt.f0) THEN ; CALL ValueFunc(f0,n1,id,in,1) ; GOTO 1 ; END IF;
    ELSE IF (f1.le.f2) THEN
      IF (f2.gt.f0) THEN ; CALL ValueFunc(f0,n2,id,in,1) ; GOTO 1 ; END IF;
    END IF
  END IF

  ! Optimal labor conditional on firing
  ! If a point >f0 is found, return (monotonicity)
  IF (in.gt.1) THEN
    n0 = ngrid(1) ; n3 = ngrid(in)
    n1 = alpha*n0 + (one-alpha)*n3 ; CALL ValueFunc(f1,n1,id,in,0)
    n2 = alpha*n3 + (one-alpha)*n1 ; CALL ValueFunc(f2,n2,id,in,0) ; iter = 0
    DO WHILE ( abs(n0-n3).gt.tol*(abs(n2)+abs(n1)) .and. iter.lt.900 )
      iter = iter + 1
      IF (f2.gt.f1) THEN ; n0 = n1 ; n1 = n2 ; f1 = f2
        n2 = alpha*n1 + (one-alpha)*n3 ; CALL ValueFunc(f2,n2,id,in,0)
      ELSE ; n3 = n2 ; n2 = n1 ; f2 = f1
        n1 = alpha*n2 + (one-alpha)*n0 ; CALL ValueFunc(f1,n1,id,in,0)
      END IF
    END DO
    IF (f1.gt.f2) THEN
      IF (f1.gt.f0) THEN ; CALL ValueFunc(f0,n1,id,in,1) ; GOTO 1 ; END IF;
    ELSE IF (f1.le.f2) THEN
      IF (f2.gt.f0) THEN ; CALL ValueFunc(f0,n2,id,in,1) ; GOTO 1 ; END IF;
    END IF
  END IF

  1 RETURN
END SUBROUTINE FINDOPTIMALN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Given a vector of state vairables and a choice for next period's empoloyment, this
! subroutine solves the innovation problem and computes the value of the firm.
SUBROUTINE VALUEFUNC(val,npn,id,in,sav)
  USE toolkit , ONLY : INTERPOLATION,INTERPOLATE
  IMPLICIT NONE
  REAL(rp), INTENT(OUT) :: val
  REAL(rp), INTENT(IN)  :: npn
  INTEGER , INTENT(IN)  :: id,in,sav
  REAL(rp)              :: prof,w_np,cpi,nh,nf,rev,Oval,ckap,crho
  REAL(rp)              :: VN(Nd),VN0(Nd),Tdd1(Nd),mV,VIn,VNIn,rho,Eval,exx
  INTEGER               :: i_np,i

  ! Interpolate the value of next period's employment (linear interpolation)
  CALL INTERPOLATION(i_np,w_np,npn,ngrid)
  DO i=1,Nd
    VN(i) = V0(i,i_np)*w_np + V0(i,i_np-1)*(one-w_np)
  END DO

  ! Labor and profits
  nh   = MAX(zero,npn-ngrid(in))
  nf   = MAX(zero,ngrid(in)-npn)
  rev  = REVENUES(dgrid(id),npn)
  prof = rev - wrate*npn - wrate*fc*ABS(nf)

  ! Continuation value if firm does not innovate
  VNIn = SUM(Td0(id,:)*VN(:))
  Eval = -ABS(wrate*fc*npn)

  ! If innovation choices can be adjusted
  IF (EGROWTH.eq.0) THEN

    ! Intensive margin
    ckap    = kappa_0*EXP(kappa_1*dgrid(id))                        ! Cost of innovation:
    VN0(:)  = VN(:)/ckap+LOGMV(Td0(id,:)) ; mv = MAXVAL(VN0(:))     !
    VN0(:)  = EXPMV(VN0(:)-mv)                                      !
    Tdd1(:) = VN0(:)/SUM(VN0(:))                                    ! Chosen distribution from bechmark
    cpi     = SUM(Tdd1(:)*VN(:)) - ckap*mv - ckap*LOG(SUM(VN0(:)))  ! Cost of intensive margin
    VIn     = SUM(Tdd1(:)*VN(:)) - cpi                              ! Value of innovating

    ! Extensive margin
    rho  = rho_0/(rho_0+(one-rho_0)*EXPM(-(VIn-VNIn)/ckap))
    crho = - ckap*LOGM(rho_0+(one-rho_0)*EXPM(-(VIn-VNIn)/ckap)) - (one-rho)*(VIn-VNIn)
    Oval = VNIn + rho*(VIn-VNIn) - crho

  ! If innovation choices CANNOT be adjusted (for the experiment)
  ELSEIF (EGROWTH.eq.1) THEN

    ! Intensive margin
    Tdd1(:) = pdd0(id,in,:)                 ! Chosen distribution from bechmark
    cpi     = pcpi0(id,in)                  ! Cost of chosen distribution from bechmark
    VIn     = SUM(Tdd1(:)*VN(:)) - cpi      ! Value of innovating

    ! Extensive margin
    rho  = prho0(id,in)                   ! Probability of innovation from bechmark
    crho = pcrho0(id,in)                  ! Cost of probability of innovation from bechmark
    Oval = VNIn + rho*(VIn-VNIn) - crho   ! Value of innovation stage

  ! Dynamics of productivity ~ AR(1)
  ELSEIF (EGROWTH.eq.2) THEN

    cpi     = pcpi0(id,in)
    crho    = pcrho0(id,in)
    rho     = one
    Tdd1(:) = Tdf(id,:)
    Oval    = SUM(Tdf(id,:)*VN(:)) - prho0(id,in)*cpi - crho

  END IF

  ! Value function
  val = prof + (one-delta0)*beta*Oval + delta0*beta*Eval

  ! If the value function is NaN, the algorithm is interrumpted
  ! This is used in the calibration to skip invalid sets of parameters
  IF (ISNAN(val)) SKIPPOINT = 1

  ! Save choices and value function in the corresponding matrices
  IF (sav.eq.1) THEN
    V(id,in)     = val                                      ! Value
    pnp(id,in)   = npn                                      ! Next period's employment
    pin(id,in)   = i_np                                     ! Closest point in the grid to nn
    ppn(id,in)   = w_np                                     ! Relative position of npn in the grid
    pnh(id,in)   = nh                                       ! Hirings
    pnf(id,in)   = nf                                       ! Firings
    phr(id,in)   = nh/ngrid(in)                             ! Hiring rate
    pfr(id,in)   = nf/ngrid(in)                             ! Firing rate
    phs(id,in)   = zero ; IF (nh.gt.zero) phs(id,in) = one  ! Indicator of whether hiring
    pfs(id,in)   = zero ; IF (nf.gt.zero) pfs(id,in) = one  ! Indicator of whether firing
    pre(id,in)   = rev                                      ! Revenues
    ppi(id,in)   = prof                                     ! Profits
    prho(id,in)  = rho                                      ! Prob. of innovation
    pcrho(id,in) = crho                                     ! Cost of extensive margin
    pdd(id,in,:) = Tdd1(:)                                  ! Chosen distribution
    pcpi(id,in)  = cpi                                      ! Cost of intensive margin
    px(id,in)    = (one-delta0)*(crho + rho*cpi)            ! Cost of innovation
    ppd(id,in,:) = rho*Tdd1(:) + (one-rho)*Td0(id,:)        ! Efective distribution of next period's productivity
    pV(id,in,:)  = VN(:)/MAXVAL(VN(:))
  END IF
  RETURN
END SUBROUTINE VALUEFUNC
FUNCTION LOGM(x) RESULT(y)
  IMPLICIT NONE
  REAL(rp) :: x,y
  y = LOG(x)
  IF (X.LT.tol) y = -huge(x)
  RETURN
END FUNCTION LOGM
FUNCTION LOGMV(x) RESULT(y)
  IMPLICIT NONE
  REAL(rp) :: x(:),y(SIZE(x))
  INTEGER  :: i
  DO i=1,SIZE(x)
    y(i) = LOGM(x(i))
  END DO
  RETURN
END FUNCTION LOGMV
FUNCTION EXPM(x) RESULT(y)
  IMPLICIT NONE
  REAL(rp) :: x,y
  y = EXP(x)
  IF (ISNAN(y) .AND. x.LT.zero) y = zero
  IF (ISNAN(y) .AND. x.GT.zero) y = DBLE(99999.9999)
  RETURN
END FUNCTION EXPM
FUNCTION EXPMV(x) RESULT(y)
  IMPLICIT NONE
  REAL(rp) :: x(:),y(SIZE(x))
  INTEGER  :: i
  DO i=1,SIZE(x)
    y(i) = EXP(x(i))
  END DO
  RETURN
END FUNCTION EXPMV
FUNCTION REVENUES(d,n) RESULT(y)
  IMPLICIT NONE
  REAL(rp) :: d,n,y
  y = bigA*exp(d*(one-gamma))*(n**gamma)
  RETURN
END FUNCTION REVENUES

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subtoroutine finds the stationary distribution of firms.
SUBROUTINE SOLVEFIRMSDIST( )
  IMPLICIT NONE
  INTEGER  :: id,in,iter,idd,inn
  REAL(rp) :: test,Dist0(Nd,Nn),DistE(Nd,Nn),TMat(Nd,Nn)
  iter  = 0 ; test = mil
  DistE = zero ; DistE(:,1) = T0(:) ; Dist = zero ; Dist0 = DistE
  DO WHILE (test.GT.tol .and. iter.lt.5000) ; iter = iter + 1
    DO id = 1,Nd ; DO in = 1,Nn
      Dist(:,pin(id,in)  ) = Dist(:,pin(id,in)  ) + Dist0(id,in)*ppd(id,in,:)*ppn(id,in)
      Dist(:,pin(id,in)-1) = Dist(:,pin(id,in)-1) + Dist0(id,in)*ppd(id,in,:)*(one-ppn(id,in))
    END DO ; END DO
    Dist = (one-delta0)*Dist + delta0*DistE
    test = zero ; DO id = 1,Nd ; DO in = 1,Nn
      IF ( ABS(Dist(id,in)-Dist0(id,in)).GT.test) test = ABS(Dist(id,in)-Dist0(id,in))
    END DO ; END DO
    Dist0 = Dist ; Dist = zero
    IF (iter.lt.200) test = 99.999
  END DO
  Dist = Dist0/SUM(Dist0(:,:))
  RETURN
END SUBROUTINE SOLVEFIRMSDIST

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subtoroutine computes some aggregates (TFP, employment, output, ...).
SUBROUTINE AGGREGATE( )
  IMPLICIT NONE
  INTEGER  :: id,in
  REAL(rp) :: aux1,aux2

  ! Aggregate variables
  TOTL = zero ; TOTY = zero ; TOTF = zero ; TOTH = zero ; TOTLI = zero
  DO id = 1,Nd ; DO in = 1,Nn
    TOTL  = TOTL  + Dist(id,in)*pnp(id,in)                       ! Aggregate labor
    TOTY  = TOTY  + Dist(id,in)*REVENUES(dgrid(id),pnp(id,in))   ! Aggregate output
    TOTH  = TOTH  + Dist(id,in)*pnh(id,in)                       ! Aggregate hirings
    TOTF  = TOTF  + Dist(id,in)*pnf(id,in)                       ! Aggregate firings
    TOTLI = TOTLI + Dist(id,in)*px(id,in)                        ! Aggregate cost of innovation
  END DO ; END DO
  TOTFe = TOTL*delta0                                            ! Aggregate firings due to exit
  TOTPI = TOTY - wrate*TOTL - wrate*fc*( TOTF + TOTFe )          ! Aggregate profits

  ! Aggregate productivity
  TFP = zero ; aux2 = zero
  DO id = 1,Nd ; DO in = 1,Nn
    aux2 = aux2 + Dist(id,in)*(pnp(id,in)**gamma)
  END DO ; END DO
  TFP = (TOTY/aux2)**(one/(one-gamma))

  ! Average firm productivity and firm productivity growth
  TOTD = zero ; TOTDg = zero
  DO id = 1,Nd ; DO in = 1,Nn
    TOTD  = TOTD  + Dist(id,in)*EXP(dgrid(id))
    TOTDg = TOTDg + Dist(id,in)*( SUM(ppd(id,in,:)*EXP(dgrid(:)))/EXP(dgrid(id)) - one )
  END DO ; END DO
  TOTD  = TOTD/SUM(Dist(:,:))
  TOTDg = TOTDg/SUM(Dist(:,:))

  ! Share of hiring and firing firms
  TOTHF = zero ; TOTFF = zero
  DO id = 1,Nd ; DO in = 1,Nn
    TOTHF = TOTHF + Dist(id,in)*phs(id,in)
    TOTFF = TOTFF + Dist(id,in)*pfs(id,in)
  END DO ; END DO
  TOTHF = TOTHF/SUM(Dist(:,:))
  TOTFF = TOTFF/SUM(Dist(:,:))

  ! Average probability of innovation and average cost of extensive and intensive margin
  TOTR = zero ; TOTCRHO = zero ; TOTCPI = zero
  DO id = 1,Nd ; DO in = 1,Nn
    TOTR    = TOTR    + Dist(id,in)*prho(id,in)
    TOTCRHO = TOTCRHO + Dist(id,in)*pcrho(id,in)
    TOTCPI  = TOTCPI  + Dist(id,in)*pcpi(id,in)*prho(id,in)*(one-delta0)
  END DO ; END DO
  TOTR    = TOTR/SUM(Dist(:,:))
  TOTCRHO = TOTCRHO/SUM(Dist(:,:))
  TOTCPI  = TOTCPI/SUM(Dist(:,:))

  ! Transition matrix for productivity
  TRD(:,:) = zero
  DO id = 1,Nd
    TRD(id,:) = zero
    DO in=1,Nn
      TRD(id,:) = TRD(id,:) + ppd(id,in,:)*Dist(id,in)
    END DO
    TRD(id,:) = TRD(id,:)/SUM(Dist(id,:))
  END DO

  RETURN
END SUBROUTINE AGGREGATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subtoroutine computes the model-based moments used for calibration, and construct
! the vector of differences.
SUBROUTINE SETMOMENTS( )
  USE toolkit , ONLY : OLS1VAR
  IMPLICIT NONE
  INTEGER  :: i,id,in
  REAL(rp) :: aux1,aux2,totalfirms

  ! Read data moments from file
  OPEN(unit=14,file=TRIM(ADJUSTL(path))//"txtfiles/moments.txt",action= "read")
    READ(14,*) d_ncat(1)
    READ(14,*) d_ncat(2)
    READ(14,*) d_ncat(3)
    READ(14,*) d_ncat(4)
    READ(14,*) d_ncat(5)
    READ(14,*) d_ncat(6)
    READ(14,*) d_asf
    READ(14,*) d_arf
    READ(14,*) d_ash
    READ(14,*) d_arh
    READ(14,*) d_asize
    READ(14,*) d_nage0
    READ(14,*) d_stdn
    READ(14,*) d_stdn0
    READ(14,*) d_prof
  CLOSE(14)

  ! Averge firm size (population and entrants)
  m_asize = TOTL/SUM(Dist(:,:))
  m_nage0 = SUM(T0(:)*pnp(:,1))/SUM(T0(:))

  ! Size distribution
  m_ncat(:) = zero
  DO in=1,Nn ; i = NCAT(in)
    m_ncat(i) = m_ncat(i) + SUM(Dist(:,in))/SUM(Dist(:,:))
  END DO

  ! Share of firing and hiring firms, and firing and hiring rates
  m_arh = zero ; m_arf = zero
  m_ash = zero ; m_asf = zero
  DO id=1,Nd ; DO in=2,Nn
    m_arh = m_arh + Dist(id,in)*phr(id,in)
    m_arf = m_arf + Dist(id,in)*pfr(id,in)
    m_ash = m_ash + Dist(id,in)*phs(id,in)
    m_asf = m_asf + Dist(id,in)*pfs(id,in)
  END DO ; END DO
  m_arh = m_arh/(SUM(Dist)-SUM(Dist(:,1)))
  m_arf = m_arf/(SUM(Dist)-SUM(Dist(:,1)))
  m_ash = m_ash/(SUM(Dist)-SUM(Dist(:,1)))
  m_asf = m_asf/(SUM(Dist)-SUM(Dist(:,1)))
  IF (m_ash.GE.TOL) m_arh = m_arh/m_ash
  IF (m_asf.GE.TOL) m_arf = m_arf/m_asf
  IF (m_ash.LT.TOL) m_arh = zero
  IF (m_asf.LT.TOL) m_arf = zero

  ! Coef. of Variation of employment (population and entrants)
  aux1 = zero
  aux2 = zero
  DO id=1,Nd
    DO in=1,Nn
      aux1 = aux1 + Dist(id,in)*( pnp(id,in) - m_asize )**two
    END DO
    aux2 = aux2 + T0(id)*( pnp(id,1) - m_nage0 )**two
  END DO
  m_stdn  = ((aux1/SUM(Dist(:,:)))**half)/m_asize
  m_stdn0 = (aux2**half)/m_nage0

  m_prof  = (TOTPI - TOTLI)/TOTY

  ! Moment name                 Model                  Data                   Weight
  MNOM(:) = '               ' ; MOMSM(:) = zero      ; MOMSD(:) = zero      ; WMAT(:) = one        ; i = 1
  MNOM(i) = '  Size Ent     ' ; MOMSM(i) = m_nage0   ; MOMSD(i) = d_nage0   ; WMAT(i) = DBLE(5.00) ; i = i + 1
  MNOM(i) = '  CV Size      ' ; MOMSM(i) = m_stdn    ; MOMSD(i) = d_stdn    ; WMAT(i) = DBLE(5.00) ; i = i + 1
  MNOM(i) = '  CV Entrants  ' ; MOMSM(i) = m_stdn0   ; MOMSD(i) = d_stdn0   ; WMAT(i) = DBLE(5.00) ; i = i + 1
  MNOM(i) = '  Share firing ' ; MOMSM(i) = m_asf     ; MOMSD(i) = d_asf     ; WMAT(i) = DBLE(5.01) ; i = i + 1
  MNOM(i) = '  Share hiring ' ; MOMSM(i) = m_ash     ; MOMSD(i) = d_ash     ; WMAT(i) = DBLE(5.01) ; i = i + 1
  MNOM(i) = '  Rate firing  ' ; MOMSM(i) = m_arf     ; MOMSD(i) = d_arf     ; WMAT(i) = DBLE(5.01) ; i = i + 1
  MNOM(i) = '  Rate hiring  ' ; MOMSM(i) = m_arh     ; MOMSD(i) = d_arh     ; WMAT(i) = DBLE(5.01) ; i = i + 1
  MNOM(i) = '  % firms 0-5  ' ; MOMSM(i) = m_ncat(1) ; MOMSD(i) = d_ncat(1) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  % firms 6-10 ' ; MOMSM(i) = m_ncat(2) ; MOMSD(i) = d_ncat(2) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  % firms 11-15' ; MOMSM(i) = m_ncat(3) ; MOMSD(i) = d_ncat(3) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  % firms 16-20' ; MOMSM(i) = m_ncat(4) ; MOMSD(i) = d_ncat(4) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  % firms 21-25' ; MOMSM(i) = m_ncat(5) ; MOMSD(i) = d_ncat(5) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  % firms >25  ' ; MOMSM(i) = m_ncat(6) ; MOMSD(i) = d_ncat(6) ; WMAT(i) = DBLE(1.01) ; i = i + 1
  MNOM(i) = '  Profit share ' ; MOMSM(i) = m_prof    ; MOMSD(i) = d_prof    ; WMAT(i) = DBLE(0.01) ; i = i + 1

  ! Define targets as % deviation in moments
  DO i=1,SIZE(MOMS)
    MOMS(i) = MOMSM(i)/MOMSD(i) - one
  END DO

  ! If value functions is NaN at some point, return very high error
  IF (SKIPPOINT.eq.1) THEN
    DO i=1,SIZE(MOMS)
      MOMS(i) = -cien
    END DO
  END IF

  RETURN
END SUBROUTINE SETMOMENTS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE FIRMSPROBLEM
