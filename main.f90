
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! cd /Users/borjapetit/Dropbox/projects/2018_firingcost/code/
! gfortran -fopenmp -O3 -ffixed-line-length-132 -J $(pwd)/compiledfiles toolkit.f90 params.f90 firmsproblem.f90 calibration.f90 experiment.f90 main.f90 -o model

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS FILE CONTAINS THE PROGRAM THAT DRIVES THE SOLUTION OF THE MODEL.
! THE COMPILATION COMMAND IS:
! gfortran -fopenmp -O3 -ffixed-line-length-132 -J $(pwd)/compiledfiles toolkit.f90 params.f90 firmsproblem.f90 calibration.f90 experiment.f90 main.f90 -o model
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM main

USE params
USE firmsproblem , ONLY : SOLVEPROBLEM,FIND_EQUIL_wrate
USE experiment   , ONLY : RUNEXPERIMENT,SENSITIVITY
USE calibration  , ONLY : CALIBRATE
USE toolkit      , ONLY : TIMING

IMPLICIT NONE

REAL(rp) :: aux1, aux2

! Set parameter values from "params.txt" file
CALL SETPARAMETERS( )

! Print header
CALL PRINTHEADER( )

! ******************************************************************************
! SOLVING THE MODEL

1 WRITE(*,FMT="(A)",ADVANCE="YES") '  What do you want to do?                                   '
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                            '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 0 ] Solve the model for given parameters (params.txt) '
  WRITE(*,FMT="(A)",ADVANCE="YES") '          -> take value of theta from params.txt & find wage'
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 1 ] Solve the model for given parameters (params.txt) '
  WRITE(*,FMT="(A)",ADVANCE="YES") '          -> assume wage rate = 1, and compute the theta    '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 2 ] Run experiment                                    '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 3 ] Calibrate the parameters                          '
  WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 4 ] Run sensitivity analysis                          '
  WRITE(*,FMT="(A)",ADVANCE="YES") '                                                            '
  WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Your choice: '; READ (*,*) MODELMOD
  WRITE(*,FMT="(A)",ADVANCE="YES") '     '
  IF ( MODELMOD.GT.4 .OR. MODELMOD.LT.0 ) THEN
    WRITE(*,FMT="(A)",ADVANCE="YES") '  ERROR!! Choose a number between 0 and 4. Try again... '
    WRITE(*,FMT="(A)",ADVANCE="YES") '     '
    GOTO 1
  END IF

! Start timing
timing0 = TIMING(1)

! Solve the model for given parameters, finding the equilibrium wage
IF (INT(MODELMOD).EQ.0) THEN
  CALL FIND_EQUIL_WRATE(2)
  CALL WRITE_RESULTS(0)
  CALL WRITE_RESULTS(1,"solution0.txt")
END IF

! Solve the model for given parameters, assuming wage = 1
IF (INT(MODELMOD).EQ.1) THEN
  EGROWTH = 0
  wrate = one
  CALL SOLVEPROBLEM( )
  theta = FINDTHETA(TOTL,TOTY-TOTLI,wrate)
  CALL WRITE_RESULTS(0)
  CALL WRITE_RESULTS(1,"solution1.txt")
END IF

! Run experiment
IF (INT(MODELMOD).EQ.2) THEN
  CALL RUNEXPERIMENT( )
END IF

! Calibrate the parameters
IF (INT(MODELMOD).EQ.3) THEN
  wrate = one
  CALL CALIBRATE( )
  CALL WRITE_RESULTS(0)
END IF

! Run sensitivity analysis
IF (INT(MODELMOD).EQ.4) THEN
  CALL SENSITIVITY( )
END IF

! If solving the model for given parameters, you may generate some graphs
IF ( INT(MODELMOD).EQ.0 .OR. INT(MODELMOD).EQ.1) THEN
2 WRITE(*,FMT="(A)",ADVANCE="YES") '  Do you want to generate graphs? [1] Yes, [0] No '
  WRITE(*,FMT="(A)",ADVANCE="YES") '     '
  WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Your choice: '; READ (*,*) PLOTMOD
  WRITE(*,FMT="(A)",ADVANCE="YES") '     '
  IF ( PLOTMOD.NE.0 .AND. PLOTMOD.NE.1 ) THEN
    WRITE(*,FMT="(A)",ADVANCE="YES") '  ERROR!! Choose 0 or 1. Try again... '
    WRITE(*,FMT="(A)",ADVANCE="YES") '     '
    GOTO 2
  END IF
  IF (INT(PLOTMOD).EQ.1) THEN
    CALL PRINTGRAPHS( )
  END IF
END IF


! ******************************************************************************
! PRINT TIME ELAPSED

timing1 = TIMING(1)-timing0
WRITE(0,'(A1)') '  '
IF (timing1.LT.DBLE(60)) THEN
  WRITE(0,'(A15,F10.4,A10)') '  TIME ELAPSED:',timing1,' seconds  '
ELSE IF(timing1.GE.DBLE(60) .AND. timing1.LT.DBLE(60*60)) THEN
  WRITE(0,'(A15,F10.4,A10)') '  TIME ELAPSED:',timing1/DBLE(60),' minutes  '
ELSE IF (timing1.GE.DBLE(60*60)) THEN
  WRITE(0,'(A15,F10.4,A10)') '  TIME ELAPSED:',timing1/DBLE(60*60) ,' hours    '
END IF
WRITE(0,'(A1)') '  '

PAUSE 'Press any key to exit'

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END PROGRAM main
