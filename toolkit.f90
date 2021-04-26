
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TOOLKIT.F90, A TOOLKIT FOR FORTRAN90 PROGRAMMING
! BORJA PETIT, © 2019
!
! PLEASE, LOOK AT THE FILE "TOOLKIT_DOCUMENTATION.PDF" FOR A COMPLETE LIST OF
! FUNCTIONS AND SUBROUTINES INCLUDED AND DETAILS ON THEIR SINTAX
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE toolkit
IMPLICIT NONE

INTERFACE INTERPOLATE
  MODULE PROCEDURE INTERPOLATE1D,INTERPOLATE2D,INTERPOLATE3D,INTERPOLATE4D
END INTERFACE

CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MISCAELLANEOUS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE NORMALIZE(Y,X,XMAX,XMIN,S)
  IMPLICIT NONE
  DOUBLE PRECISION , INTENT(IN)    :: XMAX,XMIN
  DOUBLE PRECISION , INTENT(INOUT) :: Y,X
  INTEGER , OPTIONAL               :: S
  INTEGER                          :: SS
  IF (XMAX.LT.XMIN) THEN
    WRITE(*,*) ' ERRROR in NORMALIZE: XMAX>XMIN',XMAX,XMIN ; READ(*,*) ; RETURN
  END IF
  SS = 0
  IF (PRESENT(S) .AND. S.EQ.1) SS = 1
  IF (SS.EQ.1) THEN
    X = (EXP(Y)/(DBLE(1.00)+EXP(Y)))*(XMAX-XMIN) + XMIN
    IF (ISNAN(X) .and. Y.GT.DBLE(100.0)) X = XMAX
    RETURN
  ELSE IF (SS.EQ.0) THEN
    IF (X.GT.XMAX) THEN
      Y = DBLE(100.000)
    ELSE IF (X.LT.XMIN) THEN
      Y = - DBLE(100.000)
    ELSE
      Y = LOG((X-XMIN)/(XMAX-XMIN)/max(0.0001,DBLE(1.00)-(X-XMIN)/(XMAX-XMIN)))
    END IF
    RETURN
  END IF
END SUBROUTINE NORMALIZE

! ------------------------------------------------------------------------------

FUNCTION GRID(MAXV,MINV,N,S) RESULT(V)
  IMPLICIT NONE
  INTEGER                     :: I,N
  DOUBLE PRECISION , OPTIONAL :: S
  DOUBLE PRECISION            :: MAXV,MINV,GRID0(N),V(N),SS
  V(:)=0.0000
  SS = 1.0; IF (PRESENT(S)) SS = S
  IF (SS.Le.DBLE(0.0)) THEN
    WRITE(*,*) ' ERRROR in GRID: Spacing parameter is nonpositive'
    READ(*,*)
    RETURN
  END IF
  DO I=1,N
    GRID0(I)=DBLE(I-1)/DBLE(N-1)
  END DO
  IF (SS.gt.0.00) THEN
    DO I = 1,N
        V(I) = (GRID0(I)**SS)*(MAXV-MINV) + MINV
    END DO
  END IF
  RETURN
END FUNCTION GRID

! ------------------------------------------------------------------------------

SUBROUTINE INTERPOLATION(POS,WTH,XNOW,XGRID)
  IMPLICIT NONE
  INTEGER                       :: J,N
  DOUBLE PRECISION, INTENT(IN)  :: XNOW,XGRID(:)
  INTEGER,          INTENT(OUT) :: POS
  DOUBLE PRECISION, INTENT(OUT) :: WTH
  N = SIZE(XGRID)
  IF (XNOW.le.XGRID(1)) THEN
    POS = 2
    WTH = DBLE(0.0000)
  ELSE IF (XNOW.ge.XGRID(N)) THEN
    POS = N
    WTH = 1.00
  ELSE IF (XNOW.gt.XGRID(1) .and. XNOW.lt.XGRID(N)) THEN
    J = 2
    DO WHILE (XNOW.gt.XGRID(J))
      J = J + 1
    END DO
    POS = J
    WTH = (XNOW-XGRID(J-1))/(XGRID(J)-XGRID(J-1))
  END IF
  RETURN
END SUBROUTINE INTERPOLATION

! ------------------------------------------------------------------------------

FUNCTION INTERPOLATE1D(X1,Y1,M) RESULT(XI)
  IMPLICIT NONE
  INTEGER          :: POS
  DOUBLE PRECISION :: Y1(:),M(:),X1,XI,WTH
  IF (SIZE(Y1).NE.SIZE(M)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE: Size of X and Y do not coincide'
    READ(*,*)
    RETURN
  END IF
  CALL INTERPOLATION(POS,WTH,X1,Y1)
  XI = M(POS)*WTH + M(POS-1)*(DBLE(1.0)-WTH)
  RETURN
END FUNCTION INTERPOLATE1D

! ------------------------------------------------------------------------------

FUNCTION INTERPOLATE2D(X1,X2,Y1,Y2,M) RESULT(XI)
  IMPLICIT NONE
  INTEGER          :: POS1,POS2
  DOUBLE PRECISION :: WTH1,WTH2
  DOUBLE PRECISION :: X1,X2,XI
  DOUBLE PRECISION :: Y1(:),Y2(:),M(:,:)
  IF (SIZE(M,1).ne.SIZE(Y1)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 1st dim' , SIZE(M,1) , SIZE(Y1) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,2).ne.SIZE(Y2)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 2nd dim' , SIZE(M,2) , SIZE(Y2) ; READ(*,*) ; RETURN
  END IF
  CALL INTERPOLATION(POS1,WTH1,X1,Y1)
  CALL INTERPOLATION(POS2,WTH2,X2,Y2)
  XI = M( POS1   , POS2   )*WTH1*WTH2                   + &
       M( POS1   , POS2-1 )*WTH1*(DBLE(1.00)-WTH2)      + &
       M( POS1-1 , POS2   )*(DBLE(1.00)-WTH1)*WTH2      + &
       M( POS1-1 , POS2-1 )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)
  RETURN
END FUNCTION INTERPOLATE2D
FUNCTION INTERPOLATE3D(X1,X2,X3,Y1,Y2,Y3,M) RESULT(XI)
  IMPLICIT NONE
  INTEGER          :: POS1,POS2,POS3
  DOUBLE PRECISION :: Y1(:),Y2(:),Y3(:)
  DOUBLE PRECISION :: X1,X2,X3,XI
  DOUBLE PRECISION :: WTH1,WTH2,WTH3
  DOUBLE PRECISION :: M(:,:,:)
  IF (SIZE(M,1).ne.SIZE(Y1)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE3D: incorrect 1st dim' , SIZE(M,1) , SIZE(Y1) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,2).ne.SIZE(Y2)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE3D: incorrect 2nd dim' , SIZE(M,2) , SIZE(Y2) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,3).ne.SIZE(Y3)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE3D: incorrect 3rd dim' , SIZE(M,3) , SIZE(Y3) ; READ(*,*) ; RETURN
  END IF
  CALL INTERPOLATION(POS1,WTH1,X1,Y1)
  CALL INTERPOLATION(POS2,WTH2,X2,Y2)
  CALL INTERPOLATION(POS3,WTH3,X3,Y3)
  XI = M( POS1   , POS2   , POS3   )*WTH1*WTH2*WTH3 + &
       M( POS1   , POS2-1 , POS3   )*WTH1*(DBLE(1.00)-WTH2)*WTH3 + &
       M( POS1-1 , POS2   , POS3   )*(DBLE(1.00)-WTH1)*WTH2*WTH3 + &
       M( POS1-1 , POS2-1 , POS3   )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*WTH3 + &
       M( POS1   , POS2   , POS3-1 )*WTH1*WTH2*(DBLE(1.00)-WTH3) + &
       M( POS1   , POS2-1 , POS3-1 )*WTH1*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3) + &
       M( POS1-1 , POS2   , POS3-1 )*(DBLE(1.00)-WTH1)*WTH2*(DBLE(1.00)-WTH3) + &
       M( POS1-1 , POS2-1 , POS3-1 )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3)
  RETURN
END FUNCTION INTERPOLATE3D
FUNCTION INTERPOLATE4D(X1,X2,X3,X4,Y1,Y2,Y3,Y4,M) RESULT(XI)
  IMPLICIT NONE
  INTEGER          :: POS1,POS2,POS3,POS4
  DOUBLE PRECISION :: Y1(:),Y2(:),Y3(:),Y4(:)
  DOUBLE PRECISION :: X1,X2,X3,X4,XI
  DOUBLE PRECISION :: WTH1,WTH2,WTH3,WTH4
  DOUBLE PRECISION :: M(:,:,:,:)
  IF (SIZE(M,1).ne.SIZE(Y1)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 1st dim' , SIZE(M,1) , SIZE(Y1) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,2).ne.SIZE(Y2)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 2nd dim' , SIZE(M,2) , SIZE(Y2) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,3).ne.SIZE(Y3)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 3rd dim' , SIZE(M,3) , SIZE(Y3) ; READ(*,*) ; RETURN
  END IF
  IF (SIZE(M,4).ne.SIZE(Y4)) THEN
    WRITE(*,*) ' ERROR in INTERPOLATE4D: incorrect 4th dim' , SIZE(M,4) , SIZE(Y4) ; READ(*,*) ; RETURN
  END IF
  CALL INTERPOLATION(POS1,WTH1,X1,Y1)
  CALL INTERPOLATION(POS2,WTH2,X2,Y2)
  CALL INTERPOLATION(POS3,WTH3,X3,Y3)
  CALL INTERPOLATION(POS4,WTH4,X4,Y4)
  XI = M( POS1   , POS2   , POS3   , POS4   )*WTH1*WTH2*WTH3*WTH4 + &
       M( POS1   , POS2-1 , POS3   , POS4   )*WTH1*(DBLE(1.00)-WTH2)*WTH3*WTH4 + &
       M( POS1-1 , POS2   , POS3   , POS4   )*(DBLE(1.00)-WTH1)*WTH2*WTH3*WTH4 + &
       M( POS1-1 , POS2-1 , POS3   , POS4   )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*WTH3*WTH4 + &
       M( POS1   , POS2   , POS3-1 , POS4   )*WTH1*WTH2*(DBLE(1.00)-WTH3)*WTH4 + &
       M( POS1   , POS2-1 , POS3-1 , POS4   )*WTH1*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3)*WTH4 + &
       M( POS1-1 , POS2   , POS3-1 , POS4   )*(DBLE(1.00)-WTH1)*WTH2*(DBLE(1.00)-WTH3)*WTH4 + &
       M( POS1-1 , POS2-1 , POS3-1 , POS4   )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3)*WTH4 + &
       M( POS1   , POS2   , POS3   , POS4-1 )*WTH1*WTH2*WTH3*(DBLE(1.00)-WTH4) + &
       M( POS1   , POS2-1 , POS3   , POS4-1 )*WTH1*(DBLE(1.00)-WTH2)*WTH3*(DBLE(1.00)-WTH4) + &
       M( POS1-1 , POS2   , POS3   , POS4-1 )*(DBLE(1.00)-WTH1)*WTH2*WTH3*(DBLE(1.00)-WTH4) + &
       M( POS1-1 , POS2-1 , POS3   , POS4-1 )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*WTH3*(DBLE(1.00)-WTH4) + &
       M( POS1   , POS2   , POS3-1 , POS4-1 )*WTH1*WTH2*(DBLE(1.00)-WTH3)*(DBLE(1.00)-WTH4) + &
       M( POS1   , POS2-1 , POS3-1 , POS4-1 )*WTH1*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3)*(DBLE(1.00)-WTH4) + &
       M( POS1-1 , POS2   , POS3-1 , POS4-1 )*(DBLE(1.00)-WTH1)*WTH2*(DBLE(1.00)-WTH3)*(DBLE(1.00)-WTH4) + &
       M( POS1-1 , POS2-1 , POS3-1 , POS4-1 )*(DBLE(1.00)-WTH1)*(DBLE(1.00)-WTH2)*(DBLE(1.00)-WTH3)*(DBLE(1.00)-WTH4)
  RETURN
END FUNCTION INTERPOLATE4D

! ------------------------------------------------------------------------------

FUNCTION TIMING(MODE) RESULT(TIME)
  IMPLICIT NONE
  INTEGER , OPTIONAL :: MODE
  DOUBLE PRECISION   :: TIME,V2(8)
  INTEGER            :: V1(8),MOD
  CALL DATE_AND_TIME(VALUES=V1)
  MOD   = 1 ; IF (PRESENT(MODE)) MOD = MODE
  V2(1) = DBLE(0.000)
  V2(4) = DBLE(0.000)
  V2(8) = DBLE(0.000)
  V2(2) = DBLE(V1(2)*30*24*60*60)
  V2(3) = DBLE(V1(3)   *24*60*60)
  V2(5) = DBLE(V1(5)*      60*60)
  V2(6) = DBLE(V1(6)*         60)
  V2(7) = DBLE(V1(7)            )
  IF (MOD.EQ.1) TIME  = SUM(V2)               ! Seconds
  IF (MOD.EQ.2) TIME  = SUM(V2)/DBLE(60)      ! Minutes
  IF (MOD.EQ.3) TIME  = SUM(V2)/DBLE(60*60)   ! Hours
  RETURN
END FUNCTION TIMING

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! STATISTICS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE STATS(VAR,WEIGHT,MEANVAR,SDVAR)
  IMPLICIT NONE
  DOUBLE PRECISION , INTENT(IN)           :: VAR(:)
  DOUBLE PRECISION , INTENT(IN), OPTIONAL :: WEIGHT(:)
  DOUBLE PRECISION , INTENT(OUT)          :: MEANVAR,SDVAR
  DOUBLE PRECISION                        :: AUX1,AUX2,WEIG(SIZE(VAR))
  INTEGER                                 :: i
  WEIG(:) = DBLE(1.00000)
  MEANVAR = DBLE(0.00)
  SDVAR   = DBLE(0.00)
  IF (PRESENT(WEIGHT)) THEN
    IF (SIZE(VAR).NE.SIZE(WEIGHT)) THEN
      PRINT *, 'Error in stats!! var and weight have different size'
      RETURN
    END IF
    IF (SUM(WEIGHT).LT.DBLE(0.000000001)) THEN
      PRINT *, 'Error in stats!! weights have different size'
      RETURN
    END IF
    WEIG(:) = WEIGHT(:)
  END IF
  MEANVAR = SUM(VAR(:)*WEIG(:))/SUM(WEIG)
  AUX1    = DBLE(0.00)
  DO i=1,SIZE(VAR)
    AUX1 = AUX1 + WEIGHT(i)*((VAR(i)-MEANVAR)**DBLE(2.00))
  END DO
  AUX2  = AUX1/(SUM(WEIGHT)*DBLE(SIZE(WEIGHT)-1)/DBLE(SIZE(WEIGHT)))
  SDVAR = SQRT(AUX2)
  RETURN
END SUBROUTINE STATS

! ------------------------------------------------------------------------------

SUBROUTINE OLS1VAR(CONSTANT,SLOPE,YVEC,XVEC,WITHC,WVEC)
  IMPLICIT NONE
  DOUBLE PRECISION , INTENT(OUT)              :: CONSTANT,SLOPE
  INTEGER          , INTENT(IN)               :: WITHC
  DOUBLE PRECISION , INTENT(IN)               :: YVEC(:),XVEC(:)
  DOUBLE PRECISION , INTENT(INOUT) , OPTIONAL :: WVEC(:)
  DOUBLE PRECISION                            :: MY,MX,AUX1,AUX2,WEIGHT(SIZE(YVEC))
  INTEGER                                     :: J,NOBS
  SLOPE    = DBLE(0.00000)
  CONSTANT = DBLE(0.00000)
  AUX1     = DBLE(0.00000)
  AUX2     = DBLE(0.00000)
  MY       = DBLE(0.00000)
  MX       = DBLE(0.00000)
  NOBS     = SIZE(XVEC)
  IF (SIZE(XVEC).NE.SIZE(YVEC)) THEN
    WRITE(*,*) ' ERRROR in OLS1VAR: YVEC and/or XVEC are of incorrect size'
    READ(*,*)
    RETURN
  END IF
  IF (PRESENT(WVEC)) THEN
    WEIGHT(:) = WVEC(:)
  ELSE
    WEIGHT(:) = DBLE(1.00000)
  END IF
  IF (WITHC.eq.1) MY = SUM(YVEC(:)*WEIGHT(:))/SUM(WEIGHT(:))
  IF (WITHC.eq.1) MX = SUM(XVEC(:)*WEIGHT(:))/SUM(WEIGHT(:))
  DO J=1,NOBS
    AUX1 = AUX1 + WEIGHT(J)*(YVEC(j)-MY)*(XVEC(j)-MX)
    AUX2 = AUX2 + WEIGHT(J)*(XVEC(j)-MX)*(XVEC(j)-MX)
  END DO
  SLOPE    = AUX1/AUX2
  CONSTANT = MY - SLOPE*MX
  RETURN
END SUBROUTINE OLS1VAR

! ------------------------------------------------------------------------------

FUNCTION RANDOMNORMAL(MU,STD,SEED) RESULT(SHOCK)
  IMPLICIT NONE
  DOUBLE PRECISION , PARAMETER :: S  = 0.449871 , T  = -0.386595 , A = 0.19600
  DOUBLE PRECISION , PARAMETER :: R1 = 0.27597  , R2 = 0.27846   , B = 0.25472
  DOUBLE PRECISION :: U,V,Q,MU,STD,SHOCK
  INTEGER , OPTIONAL :: SEED(:)
  IF (PRESENT(SEED)) THEN
    CALL RANDOM_SEED(PUT=SEED)
  END IF
  DO
    CALL RANDOM_NUMBER(U)
    CALL RANDOM_NUMBER(V)
    V = 1.7156 * (V - DBLE(0.5))
    Q = (U-S)**2 + (ABS(V)-T)*(A*(ABS(V)-T) - B*(U-S))
    IF (Q < R1) EXIT
    IF (Q > R2) CYCLE
    IF (V**2 < -DBLE(4.0)*LOG(U)*U**DBLE(2.00)) EXIT
  END DO
  SHOCK = STD*V/U + MU
  RETURN
END FUNCTION RANDOMNORMAL

! ------------------------------------------------------------------------------

SUBROUTINE TAUCHEN(XVEC,RHO,MU,SIGMA,N,PMAT)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: N
  DOUBLE PRECISION, INTENT(IN)  :: RHO,MU,SIGMA,XVEC(N)
  DOUBLE PRECISION, INTENT(OUT) :: PMAT(N,N)
  DOUBLE PRECISION              :: STEP
  INTEGER                       :: I,J
  DO I=1,N
    STEP      = ( XVEC(2)-XVEC(1) )/2
    PMAT(I,1) = CDFN( (XVEC(1)+STEP-RHO*XVEC(I)-MU) / SIGMA )
    DO J=2,N-1
      STEP      = ( XVEC(J)-XVEC(J-1) )/2.0
      PMAT(I,J) = CDFN( (XVEC(J)+STEP-RHO*XVEC(I)-MU) / SIGMA ) - CDFN( (XVEC(J)-STEP-RHO*XVEC(I)-MU) / SIGMA )
    END DO
    STEP      = ( XVEC(N)-XVEC(N-1) )/2.0
    PMAT(I,N) = 1.0 - CDFN( (XVEC(N)-STEP-RHO*XVEC(I)-MU) / SIGMA )
  END DO
END SUBROUTINE TAUCHEN

! ------------------------------------------------------------------------------

FUNCTION CDFN(X) RESULT(F)
  IMPLICIT NONE
  DOUBLE PRECISION :: X,F
  DOUBLE PRECISION :: XABS, XSQ
  DOUBLE PRECISION , PARAMETER :: a0  = 0.500000000000d0 , a1  = 0.398942280444d0
  DOUBLE PRECISION , PARAMETER :: a2  = 0.399903438504d0 , a3  = 5.758854804580d0
  DOUBLE PRECISION , PARAMETER :: a4  = 29.82135578080d0 , a5  = 2.624331216790d0
  DOUBLE PRECISION , PARAMETER :: a6  = 48.69599306920d0 , a7  = 5.928857244380d0
  DOUBLE PRECISION , PARAMETER :: b0  = 0.398942280385d0 , b1  = 3.8052d-8
  DOUBLE PRECISION , PARAMETER :: b2  = 1.000006153020d0 , b3  = 3.98064794d-4
  DOUBLE PRECISION , PARAMETER :: b4  = 1.986153813640d0 , b5  = 0.151679116635d0
  DOUBLE PRECISION , PARAMETER :: b6  = 5.293303249260d0 , b7  = 4.8385912808d0
  DOUBLE PRECISION , PARAMETER :: b8  = 15.15089724510d0 , b9  = 0.742380924027d0
  DOUBLE PRECISION , PARAMETER :: b10 = 30.78993303400d0 , b11 = 3.99019417011d0
  XABS = abs(X)
  XSQ  = a0*X**2
  IF (XABS <= 1.28d0)then
      F = a0-XABS*(a1-a2*XSQ/(XSQ+a3-a4/(XSQ+a5+a6/(XSQ+a7))))
  ELSE IF (XABS <= 12.7d0)then
      F = b0*exp(-XSQ)/(XABS-b1+b2/(XABS+b3+b4/(XABS-b5+b6/(XABS+b7-b8/ &
              (XABS+b9+b10/(XABS+b11))))))
  ELSE
      F = 0d0
  END IF
  IF (X > 0d0) F = 1d0-F
  RETURN
END FUNCTION CDFN

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! LINEAR ALGEBRA
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION DIAG(MAT) RESULT(VEC)
  IMPLICIT NONE
  INTEGER          :: I
  DOUBLE PRECISION :: MAT(:,:),VEC(SIZE(MAT,DIM=1))
  DO I=1,SIZE(MAT,DIM=1)
    VEC(I) = MAT(I,I)
  END DO
  RETURN
END FUNCTION DIAG

! ------------------------------------------------------------------------------

FUNCTION MULMV(MAT1,VEC1) RESULT(MAT)
  IMPLICIT NONE
  INTEGER          :: I
  DOUBLE PRECISION :: MAT1(:,:),VEC1(:)
  DOUBLE PRECISION :: MAT(SIZE(MAT1,1))
  IF (SIZE(MAT1,2).NE.SIZE(VEC1)) THEN
    WRITE(*,*) ' ERROR in MULMV: Sizes do not coincide'
    READ(*,*)
    RETURN
  END IF
  DO I=1,SIZE(MAT1,1)
    MAT(I) = SUM(MAT1(I,:)*VEC1(:))
  END DO
  RETURN
END FUNCTION MULMV

! ------------------------------------------------------------------------------

FUNCTION MULMM(MAT1,MAT2) RESULT(MAT)
  IMPLICIT NONE
  INTEGER          :: I,J
  DOUBLE PRECISION :: MAT1(:,:),MAT2(:,:)
  DOUBLE PRECISION :: MAT(SIZE(MAT1,1),SIZE(MAT2,2))
  IF (SIZE(MAT1,2).NE.SIZE(MAT2,1)) THEN
    WRITE(*,*) ' ERROR in MULMM: Sizes do not coincide'
    READ(*,*)
    RETURN
  END IF
  DO I=1,SIZE(MAT1,1)
    DO J=1,SIZE(MAT2,2)
      MAT(I,J) = SUM(MAT1(I,:)*MAT2(:,J))
    END DO
  END DO
  RETURN
END FUNCTION MULMM

! ------------------------------------------------------------------------------

FUNCTION TRANSMAT(MAT) RESULT(MATT)
  IMPLICIT NONE
  INTEGER          :: I,J
  DOUBLE PRECISION :: MAT(:,:)
  DOUBLE PRECISION :: MATT(SIZE(MAT,2),SIZE(MAT,1))
  DO I=1,SIZE(MAT,1)
    DO J=1,SIZE(MAT,2)
      MATT(J,I) = MAT(I,J)
    END DO
  END DO
  RETURN
END FUNCTION TRANSMAT

! ------------------------------------------------------------------------------

FUNCTION INVERSE(M) RESULT(IM)
  IMPLICIT NONE
  INTEGER          :: N,I,J,K
  DOUBLE PRECISION :: M(:,:),COEFF
  DOUBLE PRECISION :: IM(SIZE(M,1),SIZE(M,1)),B(SIZE(M,1)),D(SIZE(M,1)),X(SIZE(M,1))
  DOUBLE PRECISION :: L(SIZE(M,1),SIZE(M,1)),U(SIZE(M,1),SIZE(M,1)),MB(SIZE(M,1),SIZE(M,1))
  IF (SIZE(M,1).NE.SIZE(M,2)) THEN
    WRITE(*,*) ' ERROR in INVERSE: Matrix is not square'
    READ(*,*)
    RETURN
  END IF
  N  = SIZE(M,1)
  L  = 0.0
  U  = 0.0
  B  = 0.0
  MB = M
  DO K=1,N-1
    DO I=K+1,N
      COEFF  = M(I,K)/M(K,K)
      L(I,K) = COEFF
      DO J=K+1,N
        M(I,J) = M(I,J)-COEFF*M(K,J)
      END DO
    END DO
  END DO
  DO I=1,N
    L(I,I) = 1.0
  END DO
  DO J=1,N
    DO I=1,J
      U(I,J) = M(I,J)
    END DO
  END DO
  DO K=1,N
    B(K) = 1.0
    D(1) = B(1)
    DO I=2,N
      D(I) = B(I)
      DO J=1,I-1
        D(I) = D(I) - L(I,J)*D(J)
      END DO
    END DO
    X(N)=D(N)/U(N,N)
    DO I = N-1,1,-1
      X(I) = D(I)
      DO J=N,I+1,-1
        X(I)=X(I)-U(I,J)*X(J)
      END DO
      X(I) = X(I)/U(I,I)
      END DO
     DO I=1,N
       IM(I,K) = X(I)
     END DO
    B(K) = 0.0
  END DO
  M = MB
  RETURN
END FUNCTION INVERSE

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! OPTIMIZATION
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE GOLDEN(FUNC,X,Y,XMAX,XMIN,ITERMAX,TOL)
  IMPLICIT NONE

  DOUBLE PRECISION , INTENT(IN)            :: XMAX,XMIN
  DOUBLE PRECISION , INTENT(OUT)           :: X,Y
  DOUBLE PRECISION , INTENT(IN) , OPTIONAL :: TOL
  INTEGER          , INTENT(IN) , OPTIONAL :: ITERMAX

  DOUBLE PRECISION , PARAMETER             :: ALPHA=0.61803399
  DOUBLE PRECISION , PARAMETER             :: BETA=DBLE(1.0)-ALPHA
  DOUBLE PRECISION                         :: X0,X1,X2,X3,F1,F2,TOLER
  INTEGER                                  :: ITER,MAXITER

  INTERFACE
    FUNCTION FUNC(X) RESULT(F)
    DOUBLE PRECISION :: X,F
    END FUNCTION FUNC
  END INTERFACE

  TOLER   = 1.0D-8 ; IF (PRESENT(TOL)    ) TOLER   = TOL
  MAXITER = 500    ; IF (PRESENT(ITERMAX)) MAXITER = ITERMAX

  ITER = 0

  X0 = XMIN
  X3 = XMAX
  X1 = ALPHA*X0 + BETA*X3
  X2 = ALPHA*X3 + BETA*X1

  F1 = FUNC(X1)
  F2 = FUNC(X2)

  DO WHILE (ABS(X0-X3).GT.TOLER*(ABS(X2)+ABS(X1)) .AND. ITER.LT.MAXITER)
    iter = iter + 1
    IF (F2.gt.F1) THEN
      X0 = X1
      X1 = X2
      F1 = F2
      X2 = ALPHA*X1 + BETA*X3
      F2 = FUNC(X2)
    ELSE
      X3 = X2
      X2 = X1
      F2 = F1
      X1 = ALPHA*X2 + BETA*X0
      F1 = FUNC(X1)
    END IF
  END DO

  IF (F1.gt.F2) THEN
    Y = F1
    X = X1
  ELSE
    Y = F2
    X = X2
  END IF

  RETURN
END SUBROUTINE GOLDEN

! ------------------------------------------------------------------------------

SUBROUTINE GRIDSEARCH(FUNC,X,Y,IM,XMAX,XMIN,N)
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: N
  DOUBLE PRECISION , INTENT(IN)  :: XMIN,XMAX
  DOUBLE PRECISION , INTENT(OUT) :: X,Y
  INTEGER          , INTENT(OUT) :: IM
  INTEGER                        :: I,M(1)
  DOUBLE PRECISION               :: XGRID(N),F(N)
  INTERFACE
    FUNCTION FUNC(X) RESULT(Y)
      DOUBLE PRECISION :: X,Y
    END FUNCTION FUNC
  END INTERFACE
  XGRID = GRID(XMAX,XMIN,N,DBLE(1.00))
  DO I=1,N
    F(I) = FUNC(XGRID(I))
  END DO
  M  = MAXLOC(F)
  IM = M(1)
  Y  = F(IM)
  X  = XGRID(IM)
  RETURN
END SUBROUTINE GRIDSEARCH

! ------------------------------------------------------------------------------

SUBROUTINE SIMPLEX(FUNC,X,Y,IY,IND,X0,ITERMAX,TOL,IPRINT)
  IMPLICIT NONE

  INTEGER          , INTENT(OUT)           :: IY,IND
  DOUBLE PRECISION , INTENT(OUT)           :: X(:),Y
  DOUBLE PRECISION , INTENT(IN)            :: X0(:)
  DOUBLE PRECISION , INTENT(IN) , OPTIONAL :: TOL
  INTEGER          , INTENT(IN) , OPTIONAL :: ITERMAX
  INTEGER          , INTENT(IN) , OPTIONAL :: IPRINT

  INTEGER                      :: N,I,J,ILOW,IHIGH,IHIGH2,MAXITER
  DOUBLE PRECISION             :: XP(SIZE(X,1),SIZE(X,1)+1),YP(SIZE(X,1)+1),XR(SIZE(X,1))
  DOUBLE PRECISION             :: XM(SIZE(X,1)),XE(SIZE(X,1)),XC(SIZE(X,1))
  DOUBLE PRECISION             :: Y0,YR,YE,YC,CENT,YHIGH,TOLER
  DOUBLE PRECISION , PARAMETER :: ALPHA=1.0,BETA=0.5,GAMMA=2.0

  INTERFACE
    FUNCTION FUNC(X) RESULT(RESID)
    DOUBLE PRECISION :: X(:),RESID
    END FUNCTION FUNC
  END INTERFACE

  XP = 0.0D-10
  XM = 0.0D-10
  XE = 0.0D-10
  XC = 0.0D-10
  XR = 0.0D-10
  YP = 0.0D-10
  Y0 = 0.0D-10
  YR = 0.0D-10
  YE = 0.0D-10
  YC = 0.0D-10

  TOLER   = 1.0D-8 ; IF (PRESENT(TOL)    ) TOLER   = TOL
  MAXITER = 500    ; IF (PRESENT(ITERMAX)) MAXITER = ITERMAX

  IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) THEN
    WRITE(*,*) '  '
    WRITE(*,*) ' STARTING SIMPLEX ALGORTIHM'
    WRITE(*,*) '  '
  END IF

  IND = 0
  IY  = 0
  N   = SIZE(X,1)

  CALL EVALUATE(Y0,X0)

  XP(:,1) = X0
  YP(1)   = Y0

  IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) PRINT 93 , Y0

  IF (Y0.LT.TOLER) THEN
    XP(:,1) = X0
    YP(1)   = Y0
    ILOW    = 1
    GOTO 4
  END IF
  IF (IY.GE.MAXITER) GOTO 2

  IF (PRESENT(IPRINT) .AND. IPRINT.GT.1) THEN
    WRITE(*,*) '  '
    WRITE(*,*) ' Evaluating initial Simplex'
    WRITE(*,*) '  '
  END IF

  XP = DBLE(0.0)
  DO J=1,N+1
    XP(:,J) = X0(:)
  END DO

  DO I=2,N+1

    IF (X0(i-1).LT.DBLE(0.00)) XP(I-1,I) = X0(I-1) - DBLE(0.20)
    IF (X0(i-1).GE.DBLE(0.00)) XP(I-1,I) = X0(I-1) + DBLE(0.20)

    CALL EVALUATE(YP(I),XP(:,I))

    IF (PRESENT(IPRINT) .AND. IPRINT.GT.1 .and. N.LT.5) PRINT 97 , I-1, YP(I) , XP(:,I)
    IF (PRESENT(IPRINT) .AND. IPRINT.GT.1 .and. N.GE.5) PRINT 97 , I-1, YP(I)

    IF (YP(I).LT.TOLER) THEN
      YP(ILOW)   = YP(I)
      XP(:,ILOW) = XP(:,I)
      GOTO 4
    END IF
    IF (IY.GE.MAXITER) GOTO 2

  END DO

  IF (PRESENT(IPRINT) .AND. IPRINT.GT.1) THEN
    WRITE(*,*) '  '
    WRITE(*,*) ' Beginning iterations'
    WRITE(*,*) '  '
  END IF

  1 ILOW = 1
    DO I=1,N+1
      IF (YP(I).LT.YP(ILOW)) ILOW = I
    END DO
    IHIGH = 1
    DO I=1,N+1
      IF (YP(I).GT.YP(IHIGH)) IHIGH = I
    END DO
    IF (PRESENT(IPRINT) .AND. IPRINT.GT.1) PRINT 97 , IY, YP(ILOW)

    IHIGH2 = ILOW
    DO I=1,N+1
      IF (YP(I).GT.YP(IHIGH2) .AND. YP(I).LT.YP(IHIGH)) IHIGH2 = IHIGH
    END DO
    YHIGH = YP(IHIGH)
    IF (MAXVAL(ABS(XP(:,IHIGH)-XP(:,ILOW))).LT.TOLER) GOTO 3

    DO J=1,N
      CENT = 0.0
      DO I=1,N+1
        IF (I.NE.IHIGH) THEN
          CENT = CENT + XP(J,I)
        END IF
      END DO
      XM(J) = CENT/DBLE(N)
    END DO
    DO J=1,N
      XR(J) = (DBLE(1.0)+ALPHA)*XM(J) - ALPHA*XP(J,IHIGH)
    END DO

    CALL EVALUATE(YR,XR)

    IF (YR.LT.TOLER) THEN
      YP(ILOW)   = YR
      XP(:,ILOW) = XR(:)
      GOTO 4
    END IF
    IF (IY.GE.MAXITER ) GOTO 2

    IF (YR.LT.YP(ILOW)) THEN

      DO J=1,N
        XE(J) = (DBLE(1.0)-GAMMA)*XM(J) - GAMMA*XP(J,IHIGH)
      END DO

      CALL EVALUATE(YE,XE)

      IF (YE.LT.TOLER) THEN
        YP(ILOW)   = YE
        XP(:,ILOW) = XE(:)
        GOTO 4
      END IF
      IF (IY.GE.MAXITER) GOTO 2

      IF (YE.LT.YR) THEN
        XP(:,IHIGH) = XE(:)
        YP(IHIGH)   = YE
        YHIGH       = YE
      ELSE
        XP(:,IHIGH) = XR(:)
        YP(IHIGH)   = YR
        YHIGH       = YR
      END IF

    ELSEIF (YR.GT.YP(IHIGH2)) THEN

      IF (YR.LT.YP(IHIGH)) THEN
        XP(:,IHIGH) = XR(:)
        YP(IHIGH)   = YR
        YHIGH       = YR
      END IF

      DO J=1,N
        XC(J) = (DBLE(1.0)-BETA)*XM(J) - BETA*XP(J,IHIGH)
      END DO

      CALL EVALUATE(YC,XC)

      IF (YC.LT.TOLER) THEN
        YP(ILOW)   = YC
        XP(:,ILOW) = XC
        GOTO 4
      END IF
      IF (IY.GE.MAXITER ) GOTO 2

      IF (YC.LT.YP(IHIGH)) THEN

        XP(:,IHIGH) = XC(:)
        YP(IHIGH)   = YC
        YHIGH       = YC

      ELSE

        DO I=1,N+1
          IF (I.NE.ILOW) THEN
            DO J=1,N
              XM(J)   = DBLE(0.5)*( XP(J,I) + XP(J,ILOW) )
              XP(J,I) = XM(J)
            END DO
            CALL EVALUATE(YP(I),XP(:,I))
            IF (YP(I).LT.TOLER) THEN
              YP(ILOW)   = YP(I)
              XP(:,ILOW) = XP(:,I)
              GOTO 6
            END IF
            IF (IY.GE.MAXITER) GOTO 2
          END IF
        END DO

      END IF

    ELSE

      XP(:,IHIGH) = XR(:)
      YP(IHIGH)   = YR
      YHIGH       = YR

    END IF

    GOTO 1

  2 IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) '  '
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) ' SIMPLEX FINISHED: Max number of iterations'
    IND = 2
    GOTO 5

  3 IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) '  '
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) ' SIMPLEX FINISHED: System has converged, Simplex is very small'
    IND = 1
    GOTO 5

  4 IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) '  '
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,*) ' SIMPLEX FINISHED: Solution found'
    IND = 1
    GOTO 6

  5 IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*, *)  '  '
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,94) IY

    ILOW = 1
    DO I=1,N+1
      IF (YP(I).LT.YP(ILOW)) ILOW = I
    END DO

  6 IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*, *) XP(:,ILOW)
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,93) YP(ILOW)
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*,92) YP(ILOW)/Y0 - 1.0
    IF (PRESENT(IPRINT) .AND. IPRINT.GE.1) WRITE (*, *) '  '

    X = XP(:,ILOW)
    Y = YP(ILOW)

  RETURN
  97 FORMAT ('  New point    ',I4,'  ',5(F10.4))
  94 FORMAT ('  Iterations = ',I4)
  93 FORMAT ('  Error      = ',F10.4)
  92 FORMAT ('  Reduction  = ',F10.4)
  CONTAINS
  SUBROUTINE EVALUATE(YEVAL,XEVAL)
    IMPLICIT NONE
    DOUBLE PRECISION , INTENT(IN)  :: XEVAL(:)
    DOUBLE PRECISION , INTENT(OUT) :: YEVAL
    YEVAL = FUNC(XEVAL)
    IF (ISNAN(YEVAL)) YEVAL = DBLE(1.10)*YHIGH
    IY = IY + 1
  END SUBROUTINE
END SUBROUTINE SIMPLEX

! ------------------------------------------------------------------------------

SUBROUTINE LMMIN(FUNC,X,Y,IY,IND,X0,ITERMAX,DAMP,TOL,SHCK,USEBRO,IPRINT)
  IMPLICIT NONE

  DOUBLE PRECISION , INTENT(OUT)           :: X(:),Y(:)
  INTEGER          , INTENT(OUT)           :: IY,IND

  DOUBLE PRECISION , INTENT(IN)            :: X0(:)
  DOUBLE PRECISION , INTENT(IN) , OPTIONAL :: SHCK
  DOUBLE PRECISION , INTENT(IN) , OPTIONAL :: DAMP
  DOUBLE PRECISION , INTENT(IN) , OPTIONAL :: TOL
  INTEGER          , INTENT(IN) , OPTIONAL :: ITERMAX
  INTEGER          , INTENT(IN) , OPTIONAL :: IPRINT
  INTEGER          , INTENT(IN) , OPTIONAL :: USEBRO

  INTEGER                      ::  I,K,Q,MAXITER,N,M,IPRIN
  DOUBLE PRECISION             ::  J(SIZE(Y,1),SIZE(X,1))
  DOUBLE PRECISION             :: IJ(SIZE(Y,1),SIZE(X,1))
  DOUBLE PRECISION             :: JJ(SIZE(X,1),SIZE(X,1))
  DOUBLE PRECISION             :: JA(SIZE(X,1),SIZE(X,1))
  DOUBLE PRECISION             :: JT(SIZE(X,1))
  DOUBLE PRECISION             :: Y0(SIZE(Y,1))
  DOUBLE PRECISION             :: YB(SIZE(Y,1)),XB(SIZE(X,1))
  DOUBLE PRECISION             :: YA(SIZE(Y,1)),XA(SIZE(X,1))
  DOUBLE PRECISION             :: Y1(SIZE(Y,1)),X1(SIZE(X,1))
  DOUBLE PRECISION             :: DAMP0,TOLER,SHOCK
  DOUBLE PRECISION , PARAMETER :: NU = DBLE(2.0)

  INTERFACE
    FUNCTION FUNC(X) RESULT(RESID)
      DOUBLE PRECISION :: X(:)
      DOUBLE PRECISION,ALLOCATABLE :: RESID(:)
    END FUNCTION FUNC
  END INTERFACE

  DAMP0   = 1.0D-0 ; IF (PRESENT(DAMP)   ) DAMP0   = DAMP
  TOLER   = 1.0D-8 ; IF (PRESENT(TOL)    ) TOLER   = TOL
  SHOCK   = 5.0D-2 ; IF (PRESENT(SHCK)   ) SHOCK   = SHCK
  MAXITER = 500    ; IF (PRESENT(ITERMAX)) MAXITER = ITERMAX
  IPRIN   = 0      ; IF (PRESENT(IPRINT) ) IPRIN   = IPRINT

  IF (IPRIN.GE.1) THEN
    WRITE(*,*) '                          '
    WRITE(*,*) ' STARTING LMMIN ALGORTIHM '
    WRITE(*,*) '                          '
  END IF

  N = SIZE(X,1)
  M = SIZE(Y,1)

  ! Initial point
  Y0 = FUNC(X0)
  IY = 1

  YB = Y0
  XB = X0

  IF ( IPRIN.GE.1 ) WRITE(*,93) SUM(Y0(:)*Y0(:))
  IF ( IPRIN.EQ.1 ) WRITE(*, *) '  '

  ! Initialize iteration
  1 Q  = 0

    ! Check convergence
    IF (SUM(YB(:)*YB(:)).LT.TOLER) THEN
      GOTO 4
    END IF

    ! Check number of iterations
    IF (IY.GE.MAXITER) GOTO 3

    ! Compute Numerical Jacobian
    IF (PRESENT(USEBRO) .AND. USEBRO.EQ.1 .AND. IY.GT.1) THEN
      CALL BROYDEN(JJ,J,X1,XA,Y1,YA)
      J = JJ
    ELSE
      IF ( IPRIN.EQ.2 ) THEN
        J = JACOBIAN(FUNC,XB,YB,MIN(SHOCK,SHOCK*MAX(0.300,SUM(YB(:)*YB(:)))),1)
      ELSE
        J = JACOBIAN(FUNC,XB,YB,MIN(SHOCK,SHOCK*MAX(0.300,SUM(YB(:)*YB(:)))),0)
      END IF
    END IF

    ! Start iteration from best point
    XA = XB
    YA = YB

    ! Compute new point
    DO I=1,N
      DO K=1,N
        JJ(I,K) = SUM(J(:,I)*J(:,K))
      END DO
      JT(I) = SUM(J(:,I)*YA(:))
    END DO

    ! Check is JAC is close to 0
    IF (MAXVAL(ABS(JT)).LT.TOLER) GOTO 5

  2 JA = JJ
    DO I=1,N
      JA(I,I) = (DBLE(1.0)+DAMP0)*JJ(I,I)
    END DO
    IJ = INVERSE(JA)
    DO I=1,N
      X1(I) = XA(I) - SUM(IJ(I,:)*JT(:))
    END DO

    ! Check is step in X is close to 0
    IF (SUM(ABS(X1-XA)).LE.TOLER) GOTO 6

    ! Evaluate new point
    Y1 = FUNC(X1)
    IY = IY + 1

    ! Print iteration
    IF (IPRIN.GT.1) PRINT 97 , IY , SUM(Y1(:)*Y1(:)) , SUM(YB(:)*YB(:)) , DAMP0

    ! Check convergence
    IF (SUM(Y1(:)*Y1(:)).LT.TOLER) THEN
      YB = Y1
      XB = X1
      GOTO 4
    END IF

    ! Check number of iterations
    IF (IY.GE.MAXITER) GOTO 3

    ! Q = 0: initial point
    ! Q = 1: previous point was good
    ! Q = 2: previous point was bad

    IF (SUM(Y1(:)*Y1(:)).LT.SUM(YB(:)*YB(:))) THEN

      ! Update the best point
      YB = Y1 ; XB = X1

      ! If previous point was bad, start the process again
      IF (Q.EQ.2) GOTO 1

      ! If prevoous point was good but dumping factor is too small
      IF (Q.EQ.1 .AND. DAMP0.LT.0.0001) GOTO 1

      ! If previous point was good, try again with new damping factor (until finding a bad point)
      DAMP0 = DAMP0/NU ; Q = 1 ; GOTO 2

    ELSEIF (SUM(Y1(:)*Y1(:)).GE.SUM(YB(:)*YB(:))) THEN

      ! Increase damping factor
      DAMP0 = DAMP0*NU

      ! If previous point was good, initialize the process again
      IF (Q.EQ.1) GOTO 1

      IF (DAMP0.GT.DBLE(5000.0)) GOTO 6
      ! If previous point was bad, try with new damping factor (until finding a good point)
      Q = 2 ; GOTO 2

    END IF

 3 IND = 3 ; GOTO 7
 4 IND = 0 ; GOTO 7
 5 IND = 1 ; GOTO 7
 6 IND = 2 ; GOTO 7

 7 IF (IPRIN.GE.1) THEN

     WRITE(*,*) '  '

     IF (IND.EQ.3) WRITE(*,*) ' LMMIN FINISHED: Max. number of iterations reached'
     IF (IND.EQ.0) WRITE(*,*) ' LMMIN FINISHED: function is close to 0 '
     IF (IND.EQ.1) WRITE(*,*) ' LMMIN FINISHED: Jacobian is close to 0'
     IF (IND.EQ.2) WRITE(*,*) ' LMMIN FINISHED: step in X close to 0'

     WRITE(*, *) '  '
     WRITE(*,94) IY
     WRITE(*,93) SUM(YB(:)*YB(:))
     WRITE(*,92) SUM(YB(:)*YB(:))/SUM(Y0(:)*Y0(:)) - 1.0
     WRITE(*, *) '  '

   END IF

    Y = YB
    X = XB

  RETURN
  97 FORMAT ('  New point | Iteration = ',I4,' | New point = ',F10.4,' | Best so far = ',F10.4,' | Damping = ',F10.4)
  94 FORMAT ('  Iterations = ',I4)
  93 FORMAT ('  Error      = ',F10.4)
  92 FORMAT ('  Reduction  = ',F10.4)
END SUBROUTINE LMMIN

! ------------------------------------------------------------------------------

FUNCTION JACOBIAN(FUNC,X,Y,SHCK,IPRINT) RESULT(JAC)
  IMPLICIT NONE
  DOUBLE PRECISION              :: X(:),Y(:),SHCK,SHCK2
  INTEGER,OPTIONAL              :: IPRINT
  DOUBLE PRECISION, ALLOCATABLE :: JAC(:,:),YP(:,:),XP(:,:)
  INTEGER                       :: N,M,I,J,IJ
  INTERFACE
    FUNCTION FUNC(X) RESULT(RESID)
    DOUBLE PRECISION :: X(:)
    DOUBLE PRECISION, ALLOCATABLE :: RESID(:)
    END FUNCTION FUNC
  END INTERFACE
  N = SIZE(X) ; M = SIZE(Y) ; ALLOCATE(JAC(M,N),YP(M,N),XP(N,N))
  IF (PRESENT(IPRINT) .AND. IPRINT.EQ.1) PRINT *, '  '
  DO I=1,N ; IJ = 0
    XP(:,I) = X
    SHCK2   = DBLE(1.000)
  7 IF (X(I).GT.DBLE(0.00)) THEN
      IF (ABS(X(I)).GE.DBLE(1.00)) XP(I,I) = X(I) + SHCK2*SHCK*ABS(X(I))
      IF (ABS(X(I)).LT.DBLE(1.00)) XP(I,I) = X(I) + SHCK2*SHCK
    ELSE IF (X(I).LE.DBLE(0.00)) THEN
      IF (ABS(X(I)).GE.DBLE(1.00)) XP(I,I) = X(I) - SHCK2*SHCK*ABS(X(I))
      IF (ABS(X(I)).LT.DBLE(1.00)) XP(I,I) = X(I) - SHCK2*SHCK
    END IF
    YP(:,I) = FUNC(XP(:,I))
    IF (PRESENT(IPRINT) .AND. IPRINT.EQ.1) WRITE (*,92) I,SUM(YP(:,I)*YP(:,I))
    IF (SUM(YP(:,I)*YP(:,I)).GT.DBLE(50.0)*SUM(Y(:)*Y(:)) .AND. SUM(Y(:)*Y(:)).GT.DBLE(0.001) .AND. IJ.LT.10) THEN
      SHCK2 = DBLE(0.5)*SHCK2
      IJ    = IJ + 1
      GOTO 7
    END IF
    IF (MAXVAL(ABS(YP(:,I)-Y(:))).LE.DBLE(0.0001) .AND. IJ.LT.10) THEN
      SHCK2 = -DBLE(2.0)*SHCK2
      IJ    = IJ + 1
      GOTO 7
    END IF
  END DO
  DO J=1,M
    DO I=1,N
      JAC(J,I) = ( YP(J,I) - Y(J) )/( XP(I,I) - X(I) )
    END DO
  END DO
  IF (PRESENT(IPRINT) .AND. IPRINT.EQ.1) PRINT *, '  '
  92 FORMAT ('   Jacobian, Param =',I4,F10.4)
  RETURN
END FUNCTION JACOBIAN

! ------------------------------------------------------------------------------

SUBROUTINE BROYDEN(J1,J0,X1,X0,F1,F0,S)
  IMPLICIT NONE

  DOUBLE PRECISION , INTENT(IN)            :: X1(:),F1(:),X0(:),F0(:),J0(:,:)
  INTEGER          , INTENT(IN) , OPTIONAL :: S
  DOUBLE PRECISION , INTENT(OUT)           :: J1(:,:)

  DOUBLE PRECISION                         :: JDF(SIZE(F0,1)),DF(SIZE(F0,1))
  DOUBLE PRECISION                         :: DX(SIZE(X1,1)),DX2(SIZE(X1,1))
  INTEGER                                  :: I,K,N,M

  N = SIZE(X1)
  M = SIZE(F0)

  DX(:) = X1(:)-X0(:)
  DF(:) = F1(:)-F0(:)

  IF (PRESENT(S) .AND. S.eq.0 .and. N.eq.M) THEN

    DO I=1,M
     JDF(i) = SUM(J0(I,:)*DF(:))
    END DO
    DO I=1,N
      DO K=1,N
        J1(I,1:N) = J0(I,1:N) +  (DX(I) - JDF(I))*SUM(DX(:)*J0(:,K))/SUM(DX(:)*JDF(:))
      END DO
    END DO

  ELSE

    DX2(:) = DX(:)*DX(:)
    DO I=1,M
      DO K=1,N
        J1(I,K) = J0(I,K) + DX(K)*( DF(I)-SUM(J0(I,:)*DX(:)))/SUM(DX(:)*DX(:))
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE BROYDEN

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE toolkit
