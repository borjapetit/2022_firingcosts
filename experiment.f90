
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS MODULE CONTAINS TWO IMPORTANT SUBROUTINES: ONE CONTROLLING THE EXPERIMENTS,
! AND ANOTHER ONE CONTROLLING THE SENSITIVITTY ANALYSIS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE experiment
USE params
USE firmsproblem , ONLY : SOLVEPROBLEM,FIND_EQUIL_wrate
IMPLICIT NONE
REAL (rp) :: a_d0
CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine controls the experiment
SUBROUTINE RUNEXPERIMENT( )

  USE toolkit , ONLY : GRID,TIMING
  IMPLICIT NONE
  INTEGER , PARAMETER    :: NUMMOMS = 14
  INTEGER                :: ifc,j,i,NUMEXP,in,id,ii
  REAL(rp)               :: fc0,wage0,aux1,aux2,GN0(Nn,2),GN1(Nn,2)
  REAL(rp)               :: Di1(Nd,Nn),prho1(Nd,Nn),GD0(Nd),GD1(Nd),VD0(Nd),VD1(Nd),PP1(Nd,Nd),PP0(Nd,Nd),R_w_1,R_w0_1
  REAL(rp) , ALLOCATABLE :: R_PE(:,:),R_GE(:,:),R_GE0(:,:),R_fc(:),R_w(:),R_w0(:)
  REAL(rp) , ALLOCATABLE :: R_PE_1(:),R_GE_1(:),R_GE0_1(:),Di0(:,:)

    WRITE(*,FMT="(A)",ADVANCE="YES") '  ---------------------------------------------------------------'
    WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
    WRITE(*,FMT="(A)",ADVANCE="YES") '  EXPERIMENT                                                     '
    WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
    WRITE(*,FMT="(A)",ADVANCE="YES") '  Which experiment do you want to run                            '
    WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
    WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 1 ] Varying firing costs: GE & PE                          '
    WRITE(*,FMT="(A)",ADVANCE="YES") '    [ 2 ] Varying firing costs: with & w/o innovation            '
    WRITE(*,FMT="(A)",ADVANCE="YES") '    [ n ] n points in (0,...,fc_0,...,2*fc_0), min 5             '
  1 WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
    WRITE(*,FMT="(A)",ADVANCE="NO" ) '  Your choice: '; READ (*,*) NUMEXP
    WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '
    IF ( NUMEXP.LT.1 .AND. NUMEXP.GT.2 .AND. NUMEXP.LT.5 ) THEN
      WRITE(*,FMT="(A)") '  ERROR!! Chose 1, 2 or a number larger than 5. Try again... '
      GOTO 1
    END IF
    WRITE(*,FMT="(A)",ADVANCE="YES") '  ---------------------------------------------------------------'
    WRITE(*,FMT="(A)",ADVANCE="YES") '                                                                 '

  ! Firing cost euqal to calibrated value
  PRINT ('(A27)') , '  Preparing the solution...'
  PRINT ('(A37)') , '  '
  wrate = one
  CALL SOLVEPROBLEM( )
  theta = FINDTHETA(TOTL,TOTY-TOTLI,wrate)

  ! Start experiment
  fc0     = fc
  wage0   = one
  EGROWTH = 0
  timing0 = TIMING(1)

  IF (NUMEXP.EQ.1 .OR. NUMEXP.EQ.2) THEN

    ! Allocate matrices to store results
    IF (ALLOCATED(R_PE)) DEALLOCATE(R_PE) ; ALLOCATE(R_PE(NUMMOMS,4)) ; R_PE(:,:) = zero    ! General equilibrium results
    IF (ALLOCATED(R_GE)) DEALLOCATE(R_GE) ; ALLOCATE(R_GE(NUMMOMS,4)) ; R_GE(:,:) = zero    ! Partial equilibrium or exogenous innovation results
    IF (ALLOCATED(R_w))  DEALLOCATE(R_w)  ; ALLOCATE(R_w(4))          ; R_w(:)    = zero    ! Equilibrium qages
    IF (ALLOCATED(R_w0)) DEALLOCATE(R_w0) ; ALLOCATE(R_w0(4))         ; R_w0(:)   = zero    ! Equilibrium wages if exogenous innovation

    IF (NUMEXP.EQ.1) PRINT ('(A27)') , '  General equilibrium      '
    IF (NUMEXP.EQ.2) PRINT ('(A27)') , '  With innovation          '

    ! No firing costs
    PRINT ('(A37)') , '   - Solving economy with fc = 0.0   ' ; wrate = one ; fc = zero
    CALL FIND_EQUIL_WRATE(0)
    CALL WRITE_RESULTS(1,'experiment_0')
    CALL FILLRESULTS(R_GE(:,1)) ; R_w(1) = wrate
    CALL FILLINNOVATION(GD0,VD0,PP0)
    CALL FILLINNOVATIONbis(GN0)

    ! Allocate matrices to store innovation choices from the frictionless economy if experiment 2
    IF (NUMEXP.eq.2) THEN
      IF (ALLOCATED(pdd0))   DEALLOCATE(pdd0)   ; ALLOCATE(pdd0(Nd,Nn,Nd)) ; pdd0(:,:,:) = pdd(:,:,:)
      IF (ALLOCATED(pcpi0))  DEALLOCATE(pcpi0)  ; ALLOCATE(pcpi0(Nd,Nn))   ; pcpi0(:,:)  = pcpi(:,:)
      IF (ALLOCATED(prho0))  DEALLOCATE(prho0)  ; ALLOCATE(prho0(Nd,Nn))   ; prho0(:,:)  = prho(:,:)
      IF (ALLOCATED(pcrho0)) DEALLOCATE(pcrho0) ; ALLOCATE(pcrho0(Nd,Nn))  ; pcrho0(:,:) = pcrho(:,:)
      IF (ALLOCATED(Di0))    DEALLOCATE(Di0)    ; ALLOCATE(Di0(Nd,Nn))     ; Di0(:,:)    = Dist(:,:)
    END IF

    ! Solving the general equilibrium economy with fc = 0.2, fc = 0.4 and fc = 1.0
    DO ii = 1,3 ; wrate = one

      IF (ii.eq.1) fc = fc0
      IF (ii.eq.2) fc = fc0*two
      IF (ii.eq.3) fc = one

      PRINT ('(A31),F7.4') , '   - Solving economy with fc = ' , fc

      ! Find equilibrium wage
      CALL FIND_EQUIL_WRATE(0)

      ! Print results to txt file
      IF (ii.eq.1) CALL WRITE_RESULTS(1,'experiment_fc')
      IF (ii.eq.2) CALL WRITE_RESULTS(1,'experiment_2fc')
      IF (ii.eq.3) CALL WRITE_RESULTS(1,'experiment_1')

      ! Fill results andd print innovation choices
      CALL FILLRESULTS(R_GE(:,ii+1)) ; R_w(ii+1) = wrate
      CALL FILLINNOVATION(GD1,VD1,PP1)
      CALL PRINTINNOVATION(GD1,VD1,PP1,GD0,VD0,PP0,ii)

    END DO

    PRINT ('(A37)') , '     '

    ! Call GNUplot to generate the graphs
    CALL SYSTEM('gnuplot -p figures/codes/experiment_innovation.gnu')

    ! Solve counterfactual economies in partial equilibrium
    IF (NUMEXP.EQ.1) THEN ; PRINT ('(A27)') , '  Partial equilibrium      '
      PRINT ('(A37)') , '   - Solving economy with fc = fc_0  '
      wrate = R_w(1) ; fc = fc0 ; CALL SOLVEPROBLEM( ) ; CALL FILLRESULTS(R_PE(:,2)) ; R_w0(2) = wrate
      PRINT ('(A37)') , '   - Solving economy with fc = 2*fc_0'
      wrate = R_w(1) ; fc = fc0*two ; CALL SOLVEPROBLEM( ) ; CALL FILLRESULTS(R_PE(:,3)) ; R_w0(3) = wrate
      PRINT ('(A37)') , '   - Solving economy with fc = 1.0   '
      wrate = R_w(1) ; fc = one ; CALL SOLVEPROBLEM( ) ; CALL FILLRESULTS(R_PE(:,4)) ; R_w0(4) = wrate

    ! Solve counterfactual economies in general equilibrium, fixing innovation choices from frictionless economy (EGROWTH=-1)
    ELSEIF (NUMEXP.EQ.2) THEN ; PRINT ('(A27)') , '  Without innovation       ' ; EGROWTH = -1
      PRINT ('(A37)') , '   - Solving economy with fc = fc_0  '
      wrate = one ; fc = fc0 ; CALL FIND_EQUIL_WRATE(0) ; CALL FILLRESULTS(R_PE(:,2)) ; R_w0(2) = wrate
      PRINT ('(A37)') , '   - Solving economy with fc = 2*fc_0'
      wrate = one ; fc = fc0*two ; CALL FIND_EQUIL_WRATE(0) ; CALL FILLRESULTS(R_PE(:,3)) ; R_w0(3) = wrate
      PRINT ('(A37)') , '   - Solving economy with fc = 1.0   '
      wrate = one ; fc = one ; CALL FIND_EQUIL_WRATE(0) ; CALL FILLRESULTS(R_PE(:,4)) ; R_w0(4) = wrate
    END IF

    ! Print table with results
    IF (NUMEXP.EQ.1) OPEN(unit=25,file="results/experiment_1.txt",action="write")
    IF (NUMEXP.EQ.2) OPEN(unit=25,file="results/experiment_2.txt",action="write")
      WRITE(25,34) '                    '
      WRITE(25,36) ('-',i=1,80)
      WRITE(25,'(A28,A3,A24,A3,A24)') '  ','  |','         Endogenous     ','  |','        Exogenous       '
      WRITE(25,35) '                    ' , 'fc=0.0','fc=0.2','fc=0.4','fc=1.0','fc=0.2','fc=0.4','fc=1.0'
      WRITE(25,36) ('-',i=1,80)
      IF (NUMEXP.EQ.1) WRITE(25,34) '  Wage rate         ' , R_w(1:4)             ; j = 1
      IF (NUMEXP.EQ.2) WRITE(25,34) '  Wage rate         ' , R_w(1:4) , R_w0(2:4) ; j = 1
      WRITE(25,34) '  Agg. Productivity ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Firm Productivity ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Productivity grow ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Innovation Exp.   ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Output            ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Employment        ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Output per worker ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Profits           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Consumption       ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Hirings           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Firings           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Share hiring      ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Share firing      ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Var MPL           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,36) ('-',i=1,80)
      DO i=4,1,-1
        R_w(i)  = cien*(R_w(i)/R_w(1)  - one)
        R_w0(i) = cien*(R_w0(i)/R_w(1) - one)
        DO j=1,NUMMOMS
          R_GE(j,i) = cien*(R_GE(j,i)/R_GE(j,1) - one)
          R_PE(j,i) = cien*(R_PE(j,i)/R_GE(j,1) - one)
        END DO
      END DO
      IF (NUMEXP.EQ.1) WRITE(25,34) '  Wage rate         ' , R_w(1:4)             ; j = 1
      IF (NUMEXP.EQ.2) WRITE(25,34) '  Wage rate         ' , R_w(1:4) , R_w0(2:4) ; j = 1
      WRITE(25,34) '  Agg. Productivity ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Firm Productivity ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Productivity grow ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Innovation Exp.   ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Output            ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Employment        ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Output per worker ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Profits           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Consumption       ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Firings           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Hirings           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Share hiring      ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Share firing      ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,34) '  Var MPL           ' , R_GE(j,1:4) , R_PE(j,2:4) ; j = j + 1
      WRITE(25,36) ('-',i=1,80)
      WRITE(25,34) '                    '
    CLOSE(25)

  END IF


  ! Experiment 3: varying firing costs from 0 to 0.40
  IF (NUMEXP.GE.5) THEN

    ! Make sure number of points is odd (to include bechmark value in the grid)
    IF (NUMEXP.gt.2) NUMEXP = NUMEXP+1-MOD(NUMEXP,2)

    ! Allocate matrices and vectors to store restuls
    IF (ALLOCATED(R_PE))  DEALLOCATE(R_PE)  ; ALLOCATE(R_PE(NUMMOMS,NUMEXP))  ; R_PE(:,:)  = zero
    IF (ALLOCATED(R_GE))  DEALLOCATE(R_GE)  ; ALLOCATE(R_GE(NUMMOMS,NUMEXP))  ; R_GE(:,:)  = zero
    IF (ALLOCATED(R_GE0)) DEALLOCATE(R_GE0) ; ALLOCATE(R_GE0(NUMMOMS,NUMEXP)) ; R_GE0(:,:) = zero
    IF (ALLOCATED(R_w))   DEALLOCATE(R_w)   ; ALLOCATE(R_w(NUMEXP))           ; R_w(:)     = zero
    IF (ALLOCATED(R_w0))  DEALLOCATE(R_w0)  ; ALLOCATE(R_w0(NUMEXP))          ; R_w0(:)    = zero
    IF (ALLOCATED(R_fc))  DEALLOCATE(R_fc)  ; ALLOCATE(R_fc(NUMEXP))          ; R_fc(:)    = zero

    ! Define grid for firing cost parameter
    R_fc = GRID( two*fc0 , zero , NUMEXP )

    ! Allocate firm choices
    IF (ALLOCATED(pdd0))   DEALLOCATE(pdd0)   ; ALLOCATE(pdd0(Nd,Nn,Nd)) ; pdd0   = DBLE(0.00)
    IF (ALLOCATED(pcpi0))  DEALLOCATE(pcpi0)  ; ALLOCATE(pcpi0(Nd,Nn))   ; pcpi0  = DBLE(0.00)
    IF (ALLOCATED(prho0))  DEALLOCATE(prho0)  ; ALLOCATE(prho0(Nd,Nn))   ; prho0  = DBLE(0.00)
    IF (ALLOCATED(pcrho0)) DEALLOCATE(pcrho0) ; ALLOCATE(pcrho0(Nd,Nn))  ; pcrho0 = DBLE(0.00)
    IF (ALLOCATED(Di0))    DEALLOCATE(Di0)    ; ALLOCATE(Di0(Nd,Nn))     ; Di0    = DBLE(0.00)

    ! Open text files to store results (input for generating graphs)
    OPEN(unit=16,file="results/experiment_3_ge.txt",action="write")
    OPEN(unit=17,file="results/experiment_3_pe.txt",action="write")
    OPEN(unit=18,file="results/experiment_3_ge_ni.txt",action="write")
    OPEN(unit=19,file="results/experiment_3_ge_raw.txt",action="write")
    OPEN(unit=20,file="results/experiment_3_pe_raw.txt",action="write")
    OPEN(unit=21,file="results/experiment_3_ge_ni_raw.txt",action="write")

    ! Print header of text files
    WRITE(16,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(17,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(18,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(19,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(20,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(21,'(40(A10))') 'fc','wrate','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'

    ! Solve counterfactual economies
    DO ifc=1,NUMEXP ; fc = R_fc(ifc)

      PRINT ('(A31,F4.2,A3,I2,A4,I2)') , '   - Solving economy with fc = ',fc,' | ',ifc,' of ',NUMEXP

      ! Solve GE with innovation
      EGROWTH =  0 ; wrate = one ; CALL FIND_EQUIL_WRATE(0) ; CALL FILLRESULTS(R_GE(:,ifc)) ; R_w(ifc) = wrate

      ! If frictionless economy, save innovation choices for counterfactuals
      IF (fc.lt.tol) THEN
        IF (ALLOCATED(pdd0))   DEALLOCATE(pdd0)   ; ALLOCATE(pdd0(Nd,Nn,Nd)) ; pdd0(:,:,:) = pdd(:,:,:)
        IF (ALLOCATED(pcpi0))  DEALLOCATE(pcpi0)  ; ALLOCATE(pcpi0(Nd,Nn))   ; pcpi0(:,:)  = pcpi(:,:)
        IF (ALLOCATED(prho0))  DEALLOCATE(prho0)  ; ALLOCATE(prho0(Nd,Nn))   ; prho0(:,:)  = prho(:,:)
        IF (ALLOCATED(pcrho0)) DEALLOCATE(pcrho0) ; ALLOCATE(pcrho0(Nd,Nn))  ; pcrho0(:,:) = pcrho(:,:)
        IF (ALLOCATED(Di0))    DEALLOCATE(Di0)    ; ALLOCATE(Di0(Nd,Nn))     ; Di0(:,:)    = Dist(:,:)
      END IF

      ! Solve GE without innovation (innovation choices fixed from frictionless economy)
      EGROWTH = -1 ; wrate = one ; CALL FIND_EQUIL_WRATE(0) ; CALL FILLRESULTS(R_GE0(:,ifc)) ; R_w0(ifc) = wrate

      ! Solve PE (with innovation)
      EGROWTH =  0 ; wrate = R_w(1) ; CALL SOLVEPROBLEM( ) ; CALL FILLRESULTS(R_PE(:,ifc))

      ! Print results to text files
      WRITE(16,'(40(F10.4))') R_fc(ifc) , cien*(R_w(ifc)-R_w(1))  , (cien*R_GE( j,ifc)/R_GE(j,1)-cien,j=1,SIZE(R_GE(:,1)))
      WRITE(17,'(40(F10.4))') R_fc(ifc) , zero                    , (cien*R_PE( j,ifc)/R_GE(j,1)-cien,j=1,SIZE(R_GE(:,1)))
      WRITE(18,'(40(F10.4))') R_fc(ifc) , cien*(R_w0(ifc)-R_w(1)) , (cien*R_GE0(j,ifc)/R_GE(j,1)-cien,j=1,SIZE(R_GE(:,1)))
      WRITE(19,'(40(F10.4))') R_fc(ifc) , R_w(ifc)                , (R_GE( j,ifc),j=1,SIZE(R_GE(:,1)))
      WRITE(20,'(40(F10.4))') R_fc(ifc) , R_w(ifc)                , (R_PE( j,ifc),j=1,SIZE(R_GE(:,1)))
      WRITE(21,'(40(F10.4))') R_fc(ifc) , R_w0(ifc)               , (R_GE0(j,ifc),j=1,SIZE(R_GE(:,1)))

    END DO

    ! Close text files
    CLOSE (16) ; CLOSE (17) ; CLOSE (18) ; CLOSE (19) ; CLOSE (20) ; CLOSE (21)

    ! Call GNUplot to generate the graphs
    CALL SYSTEM('gnuplot -p figures/codes/experiment_change.gnu')
  END IF

  34 FORMAT (A20,F8.2,'  |',3(F8.2),'  |',3(F8.2))
  35 FORMAT (A20,A8  ,'  |',3(A8  ),'  |',3(A8  ))
  36 FORMAT ('  ',100(A1))

  RETURN
  CONTAINS
    SUBROUTINE FILLINNOVATION(GD,VD,PP)
      IMPLICIT NONE
      REAL(rp) , INTENT(OUT) :: GD(Nd),VD(Nd),PP(Nd,Nd)
      REAL(rp) :: aux
      GD(:) = zero
      DO id = 1,Nd
       PP(id,:) = zero
       DO in=1,Nn
         PP(id,:) = PP(id,:) + pdd(id,in,:)*Dist(id,in)
       END DO
       PP(id,:) = PP(id,:)/SUM(Dist(id,:))
       GD(id)   = SUM(PP(id,:)*EXP(dgrid(:)))/EXP(dgrid(id)) - one
      END DO
      DO id = 1,Nd
        aux = zero
        DO in=1,Nd
          aux = aux + PP(id,in)*( EXP(dgrid(in))/EXP(dgrid(id)) - one - GD(id) )**two
        END DO
        VD(id) = aux**half
      END DO
      RETURN
    END SUBROUTINE FILLINNOVATION
    SUBROUTINE PRINTINNOVATION(GDa,Vda,PPa,GDb,VDb,PPb,IND)
      IMPLICIT NONE
      REAL(rp) , INTENT(IN) :: GDa(Nd),Vda(Nd),PPa(Nd,Nd),GDb(Nd),VDb(Nd),PPb(Nd,Nd)
      INTEGER  , INTENT(IN) :: IND
      IF (IND.eq.1) OPEN(unit=15, file="results/experiment_innovation_d_fc0.txt", action= "write")
      IF (IND.eq.2) OPEN(unit=15, file="results/experiment_innovation_d_2fc0.txt", action= "write")
      IF (IND.eq.3) OPEN(unit=15, file="results/experiment_innovation_d_one.txt", action= "write")
       WRITE(15,'(100(A20))') 'dc0','fc1','d','ddd','ed0','ed1','vvv','vv0','vvv1'
       DO id=1,Nd
         WRITE(15,'(100(F20.15))') zero , fc , dgrid(id)   , &
           GDa(id)-GDb(id)     , GDb(id)    , GDa(id)      , &
           Vda(id)-VDb(id)     , VDb(id)    , Vda(id)
       END DO
      CLOSE(15)
      RETURN
    END SUBROUTINE PRINTINNOVATION
END SUBROUTINE RUNEXPERIMENT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine fills the vector of results with the relevant aggregates
SUBROUTINE FILLRESULTS(RVEC)
  IMPLICIT NONE
  REAL(rp) , INTENT(INOUT) :: RVEC(:)
  RVEC(1)  = TFP
  RVEC(2)  = TOTD
  RVEC(3)  = TOTDg
  RVEC(4)  = TOTLI
  RVEC(5)  = TOTY
  RVEC(6)  = TOTL
  RVEC(7)  = TOTY/TOTL
  RVEC(8)  = TOTPI
  RVEC(9)  = TOTY - TOTLI
  RVEC(10) = cien*TOTF/TOTL
  RVEC(11) = cien*TOTH/TOTL
  RVEC(12) = cien*TOTHF
  RVEC(13) = cien*TOTFF
  RVEC(14) = NOPT
  RETURN
END SUBROUTINE FILLRESULTS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine runs the sensitivity analysis included in the appendix of the paper
SUBROUTINE SENSITIVITY( )
  IMPLICIT NONE
  INTEGER  :: j
  REAL(rp) :: R_0B(14),R_1B(14),R_0(14,7,2),R_1(14,7,2),fc0,SHOCK,aux1,FCBASE

  SHOCK  = DBLE(0.10)
  FCBASE = fc

  WRITE(0,"(A)") '  ---------------------------------------------------------------'
  WRITE(0,"(A)") '                                                                 '
  WRITE(0,"(A)") '  SENSITIVITY ANALYSIS                                           '
  WRITE(0,"(A)") '                                                                 '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'

  CALL FILLRESULTSVEC(R_0B,R_1B,FCBASE)

  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = bigA
  bigA = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,1,1),R_1(:,1,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ bigA    ' , cien*(R_1(:,1,1)/R_0(:,1,1) - one)
  bigA = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,1,2),R_1(:,1,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- bigA    ' , cien*(R_1(:,1,2)/R_0(:,1,2) - one)
  bigA = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = sigma0_z
  sigma0_z = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,2,1),R_1(:,2,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ sigma0_z' , cien*(R_1(:,2,1)/R_0(:,2,1) - one)
  sigma0_z = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,2,2),R_1(:,2,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- sigma0_z' , cien*(R_1(:,2,2)/R_0(:,2,2) - one)
  sigma0_z = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = mu_z
  mu_z = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,3,1),R_1(:,3,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ mu_z    ' , cien*(R_1(:,3,1)/R_0(:,3,1) - one)
  mu_z = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,3,2),R_1(:,3,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- mu_z    ' , cien*(R_1(:,3,2)/R_0(:,3,2) - one)
  mu_z = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = sigma_z
  sigma_z = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,4,1),R_1(:,4,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ sigma_z ' , cien*(R_1(:,4,1)/R_0(:,4,1) - one)
  sigma_z = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,4,2),R_1(:,4,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- sigma_z ' , cien*(R_1(:,4,2)/R_0(:,4,2) - one)
  sigma_z = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = kappa_0
  kappa_0 = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,5,1),R_1(:,5,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ kappa_0 ' , cien*(R_1(:,5,1)/R_0(:,5,1) - one)
  kappa_0 = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,5,2),R_1(:,5,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- kappa_0 ' , cien*(R_1(:,5,2)/R_0(:,5,2) - one)
  kappa_0 = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = kappa_1
  kappa_1 = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,6,1),R_1(:,6,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ kappa_1 ' , cien*(R_1(:,6,1)/R_0(:,6,1) - one)
  kappa_1 = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,6,2),R_1(:,6,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- kappa_1 ' , cien*(R_1(:,6,2)/R_0(:,6,2) - one)
  kappa_1 = aux1

  WRITE(0,'(A10,40(F10.2))') '   '
  WRITE(0,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
  WRITE(0,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
  aux1 = rho_0
  rho_0 = aux1*(one+SHOCK) ; CALL FILLRESULTSVEC(R_0(:,7,1),R_1(:,7,1),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '+ rho_0   ' , cien*(R_1(:,7,1)/R_0(:,7,1) - one)
  rho_0 = aux1*(one-SHOCK) ; CALL FILLRESULTSVEC(R_0(:,7,2),R_1(:,7,2),FCBASE)
  WRITE(0,'(A10,40(F10.2))') '- rho_0   ' , cien*(R_1(:,7,2)/R_0(:,7,2) - one)
  rho_0 = aux1

  OPEN(unit=22,file="results/sensitivity.txt",action="write")
    WRITE(22,'(40(A10))') '                      '
    WRITE(22,'(   A65 )') '  SENSITIVITY ANALYSIS'
    WRITE(22,'(40(A10))') '                      '
    WRITE(22,'(A65)'    ) ' % fall relative to frinctionless economy. Positive shock '
    WRITE(22,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(22,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
    WRITE(22,'(A10,40(F10.2))') 'bigA      ' , cien*(R_1(:,1,1)/R_0(:,1,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'sigma0_z  ' , cien*(R_1(:,2,1)/R_0(:,2,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'mu_z      ' , cien*(R_1(:,3,1)/R_0(:,3,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'sigma_z   ' , cien*(R_1(:,4,1)/R_0(:,4,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'kappa_0   ' , cien*(R_1(:,5,1)/R_0(:,5,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'kappa_1   ' , cien*(R_1(:,6,1)/R_0(:,6,1) - one)
    WRITE(22,'(A10,40(F10.2))') 'rho_0     ' , cien*(R_1(:,7,1)/R_0(:,7,1) - one)
    WRITE(22,'(40(A10))') '                      '
    WRITE(22,'(A65)'    ) ' % fall relative to frinctionless economy. Negative shock '
    WRITE(22,'(40(A10))') '  ','tfp','d','dg','x','y','n','yn','pi','c','f','h','sh','sf','mpl'
    WRITE(22,'(A10,40(F10.2))') 'Baseline  ' , cien*(R_1B(:)/R_0B(:) - one)
    WRITE(22,'(A10,40(F10.2))') 'bigA      ' , cien*(R_1(:,1,2)/R_0(:,1,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'sigma0_z  ' , cien*(R_1(:,2,2)/R_0(:,2,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'mu_z      ' , cien*(R_1(:,3,2)/R_0(:,3,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'sigma_z   ' , cien*(R_1(:,4,2)/R_0(:,4,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'kappa_0   ' , cien*(R_1(:,5,2)/R_0(:,5,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'kappa_1   ' , cien*(R_1(:,6,2)/R_0(:,6,2) - one)
    WRITE(22,'(A10,40(F10.2))') 'rho_0     ' , cien*(R_1(:,7,2)/R_0(:,7,2) - one)
  CLOSE(22)

  RETURN
  CONTAINS
    SUBROUTINE FILLRESULTSVEC(VEC0,VEC1,FCBASE)
      IMPLICIT NONE
      REAL(rp) , INTENT(OUT) :: VEC0(:),VEC1(:)
      REAL(rp)               :: FCBASE
      INTEGER                :: i
      wrate = one
      fc = FCBASE
      CALL SOLVEPROBLEM( )
      theta = FINDTHETA(TOTL,TOTY-TOTLI,wrate)
      CALL FILLRESULTS(VEC1)
      wrate = one
      fc = zero
      CALL FIND_EQUIL_wrate(0)
      CALL FILLRESULTS(VEC0)
      RETURN
    END SUBROUTINE FILLRESULTSVEC
END SUBROUTINE SENSITIVITY

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE experiment

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
