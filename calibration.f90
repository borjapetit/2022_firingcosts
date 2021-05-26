
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS MODULE CONTAINS THE FUNCTIONS AND SUBROUTINES REQUIRED TO CALIBRATE THE PARAMETERS
! OF THE MODEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE calibration
USE params
USE firmsproblem , ONLY : SOLVEPROBLEM,FIND_EQUIL_wrate
IMPLICIT NONE
REAL(rp) :: MAXERRORCALIB
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine normalize the parameter values to lie within a specific range
SUBROUTINE SETPARAMS(VECPARS,INDICATOR)
  USE toolkit , ONLY : NORMALIZE
  IMPLICIT NONE
  INTEGER  , INTENT(IN)    :: INDICATOR
  REAL(rp) , INTENT(INOUT) :: VECPARS(:)
  INTEGER                  :: i ; i = 1
  CALL NORMALIZE ( VECPARS(i) , fc       , DBLE(0.800) , DBLE(0.000) , INDICATOR ) ; i = i + 1 ! 1
  CALL NORMALIZE ( VECPARS(i) , sigma_z  , DBLE(0.500) , DBLE(0.010) , INDICATOR ) ; i = i + 1 ! 2
  CALL NORMALIZE ( VECPARS(i) , mu_z     , DBLE(0.150) , DBLE(0.010) , INDICATOR ) ; i = i + 1 ! 3
  CALL NORMALIZE ( VECPARS(i) , bigA     , DBLE(4.000) , DBLE(1.000) , INDICATOR ) ; i = i + 1 ! 4
  CALL NORMALIZE ( VECPARS(i) , sigma0_z , DBLE(2.000) , DBLE(0.000) , INDICATOR ) ; i = i + 1 ! 5
  CALL NORMALIZE ( VECPARS(i) , rho_0    , DBLE(1.000) , DBLE(0.000) , INDICATOR ) ; i = i + 1 ! 6
  CALL NORMALIZE ( VECPARS(i) , kappa_0  , DBLE(70.00) , DBLE(0.000) , INDICATOR ) ; i = i + 1 ! 7
  CALL NORMALIZE ( VECPARS(i) , kappa_1  , DBLE(3.000) , DBLE(-1.00) , INDICATOR ) ; i = i + 1 ! 8
  RETURN
END SUBROUTINE SETPARAMS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine controls the calibration process
SUBROUTINE CALIBRATE( )
  USE toolkit , ONLY : SIMPLEX,LMMIN,TIMING
  IMPLICIT NONE
  INTEGER  , PARAMETER       :: NUMP   = 8    , SHOW0  = 2
  REAL(rp) , PARAMETER       :: SHOCK0 = 0.05 , DAMP0  = 0.10
  REAL(rp) , DIMENSION(NUMP) :: PARS1,PARS0,PARSS,PARSSS,SHOCKS
  REAL(rp)                   :: TERROR,TERROR0,TERROR1
  REAL(rp)                   :: RES(SIZE(MOMS))
  REAL(rp)                   :: TTIME,MAXTIME
  INTEGER                    :: ECODE,NITER,IITER,MAXITER

  WRITE(*,FMT="(A)",ADVANCE="YES") '  ---------------------------------------------------------------'
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
  WRITE(*,FMT="(A)",ADVANCE="YES") '  CALIBRATION                                                    '
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
  WRITE(*,FMT="(A)",ADVANCE="YES") '  Choose calibratrion algorithm                                  '
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 1 ] Levenberg-Marquardt                                    '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 2 ] Simplex                                                '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 3 ] Random search                                          '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 4 ] Levenberg-Marquardt, shocking initial values           '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 5 ] Simplex, shocking initial values                       '
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
  WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Your choice: '; READ (*,*) CALIBMOD
  IF (CALIBMOD.le.2) THEN
    WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Max iterations: '; READ (*,*) MAXITER
  ELSEIF (CALIBMOD.eq.3) THEN
    WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Max time (hours): '; READ (*,*) MAXTIME
  ELSEIF (CALIBMOD.gt.3) THEN
    WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Max iterations (each param): '; READ (*,*) MAXITER
  END IF
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
  WRITE(*,FMT="(A)",ADVANCE="YES") '  ---------------------------------------------------------------'
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '

  ! Start timing
  timing0 = TIMING(1)

  bigA = (DBLE(4.10)**(one-gamma))/gamma

  ! Set initial values
  CALL SETPARAMS(PARS0,0)

  ! Initialize max. error
  MAXERRORCALIB = 999999999.999

  ! Levenberg-Marquardt algorithm
  IF(INT(CALIBMOD).EQ.1) THEN
    CALL LMMIN(FUNC_CALIB,PARS1,RES,NITER,ECODE,PARS0,SHCK=SHOCK0,DAMP=DAMP0,IPRINT=SHOW0)

  ! SIMPLEX algorithm
  ELSE IF (INT(CALIBMOD).EQ.2) THEN
    CALL SIMPLEX(SUM_FUNC_CALIB,PARS1,TERROR,NITER,ECODE,PARS0,IPRINT=SHOW0)

  ! Random search
  ELSE IF ( INT(CALIBMOD) .EQ.  3 ) THEN
    TTIME   = TIMING(3)
    TERROR0 = SUM_FUNC_CALIB(PARS0) ! Initial error
    TERROR1 = TERROR0               ! Best error so far
    TERROR  = TERROR0               ! Error current iteration
    PRINT ('(I6,5(F12.6))') , NITER , TIMING(3)-TTIME , TERROR1 , TERROR , TERROR0
    DO WHILE ( TIMING(3)-TTIME .LT. MAXTIME ) ; CALL RANDOM_NUMBER(SHOCKS)
      PARS1(:) = PARS0(:) + PARS0(:)*SHOCK0*(SHOCKS(:)-half)
      TERROR1  = SUM_FUNC_CALIB(PARS1)
      PRINT ('(I6,5(F12.6))') , NITER , TIMING(3)-TTIME , TERROR1 , TERROR , TERROR0
      IF (TERROR1.LE.TERROR) THEN
        TERROR = TERROR1
        PARS0  = PARS1
      END IF
    END DO
    PARS1 = PARS0

  ! Levenberg-Marquardt or Simplex , shocking initial values
  ELSE IF ( INT(CALIBMOD) .GT. 3 ) THEN ! Simplex with large shocks
    IF (INT(CALIBMOD).EQ.4) CALL LMMIN(FUNC_CALIB,PARSS,RES,NITER,ECODE,PARS0,SHCK=SHOCK0,DAMP=DAMP0,IPRINT=SHOW0)
    IF (INT(CALIBMOD).EQ.5) CALL SIMPLEX(SUM_FUNC_CALIB,PARSS,TERROR,NITER,ECODE,PARS0,IPRINT=SHOW0)
      PARSSS = PARSS
    1 TERROR = SUM_FUNC_CALIB(PARSS)
      PRINT * , 0 , TERROR
    DO IITER=1,NUMP
      PARS0 = PARSS ; PARS0(IITER) = PARSS(IITER)*DBLE(1.50)
      IF (INT(CALIBMOD).EQ.4) CALL LMMIN(FUNC_CALIB,PARS1,RES,NITER,ECODE,PARS0,SHCK=SHOCK0,DAMP=DAMP0,IPRINT=SHOW0)
      IF (INT(CALIBMOD).EQ.5) CALL SIMPLEX(SUM_FUNC_CALIB,PARS1,TERROR0,NITER,ECODE,PARS0,IPRINT=SHOW0)
      TERROR0 = SUM_FUNC_CALIB(PARS1)
      PRINT * , IITER , TERROR0 , TERROR
      IF (TERROR0.LT.TERROR) THEN
        PARSSS = PARS1
        TERROR = TERROR0
      END IF
      PARS0 = PARSS ; PARS0(IITER) = PARSS(IITER)*DBLE(0.50)
      IF (INT(CALIBMOD).EQ.4) CALL LMMIN(FUNC_CALIB,PARS1,RES,NITER,ECODE,PARS0,SHCK=SHOCK0,DAMP=DAMP0,IPRINT=SHOW0)
      IF (INT(CALIBMOD).EQ.5) CALL SIMPLEX(SUM_FUNC_CALIB,PARS1,TERROR0,NITER,ECODE,PARS0,IPRINT=SHOW0)
      TERROR0 = SUM_FUNC_CALIB(PARS1)
      PRINT * , -IITER , TERROR0 , TERROR
      IF (TERROR0.LT.TERROR) THEN
        PARSSS = PARS1
        TERROR = TERROR0
      END IF
    END DO
    IF (SUM(ABS(PARSSS-PARSS)).GT.TOL) THEN
      PRINT * , '  ' ; PRINT * , '  ' ; PRINT * , '  Re-initialize ' ; PRINT * , '  '
      PARSS = PARSSS ; GOTO 1
    END IF
    PARS1 = PARSSS
  END IF

  ! Fix resulting parameter values
  CALL SETPARAMS(PARS1,1)

  ! Solve the model with new parameters, assuming wrate=1
  CALL SOLVEPROBLEM( )
  theta = FINDTHETA(TOTL,TOTY-TOTLI,one)

  RETURN
END SUBROUTINE CALIBRATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This function returns the vector differences between the model- and data-generated
! moments. It is used when applying the Levenberg-Marquardt method
FUNCTION FUNC_CALIB(VECPARS) RESULT(ER)
  IMPLICIT NONE
  REAL(rp)               :: VECPARS(:),SUMER
  REAL(rp) , ALLOCATABLE :: ER(:)
  ALLOCATE(ER(SIZE(MOMS)))
  CALL SETPARAMS(VECPARS,1)
  CALL SOLVEPROBLEM( )
  ER(:) = MOMS(:)*WMAT(:)
  SUMER = SUM(ER(:)*ER(:))
  IF ( SUMER.LT.MAXERRORCALIB ) THEN
    theta = FINDTHETA(TOTL,TOTY-TOTLI,one)
    CALL WRITE_RESULTS(1,"calibration.txt")
    MAXERRORCALIB = SUMER
  END IF
  RETURN
END FUNCTION FUNC_CALIB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This function returns the sum of squared differences between the model- and
! data-generated moments. It is used when applying the Simplex algorithm
FUNCTION SUM_FUNC_CALIB(VECPARS) RESULT(SUMER)
  IMPLICIT NONE
  REAL(rp) :: VECPARS(:),SUMER,ER(SIZE(MOMS))
  ER = FUNC_CALIB(VECPARS) ; SUMER = SUM(ER(:)*ER(:))
  RETURN
END FUNCTION SUM_FUNC_CALIB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE calibration
