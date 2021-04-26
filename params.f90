
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS MODULE CONTAINS ALL THE GLOBAR VARIABLES, AS WELL AS A NUMBER OF SIMPLE FUNCTIONS
! AND SUBROUTINE TO PRINT RESULTS AND GENERATE GRAPHS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE PARAMS
IMPLICIT NONE

CHARACTER(LEN=100) , PARAMETER :: path = "/Users/borjapetit/Dropbox/projects/2018_firingcost/code/"

! COMPUTATION PARAMETERS *******************************************************
INTEGER  , PARAMETER :: rp      = kind(1.0d0)     ! Precision for real numbers
REAL(rp) , PARAMETER :: tol     = 10.0D-7         ! Tolerance level for iterations
INTEGER  , PARAMETER :: maxiter = 900             ! max number of iterations
REAL(rp) , PARAMETER :: zero    = 0.0d10          ! 0 in REAL(rp)
REAL(rp) , PARAMETER :: one     = DBLE(1.00)      ! 1 in REAL(rp)
REAL(rp) , PARAMETER :: cien    = DBLE(100.0)     ! 100 in REAL(rp)
REAL(rp) , PARAMETER :: mil     = DBLE(1000.0)    ! 100 in REAL(rp)
REAL(rp) , PARAMETER :: half    = DBLE(0.5000)    ! 100 in REAL(rp)
REAL(rp) , PARAMETER :: two     = DBLE(2.0000)    ! 100 in REAL(rp)

! SOLUTION PARAMETERS **********************************************************
INTEGER(4)        :: threads        ! Number of threads for parallelization
INTEGER           :: Nn             ! Dimension: state of "np"
INTEGER           :: Nd             ! Dimension: state of "d"
INTEGER           :: MODELMOD       ! Solution options
INTEGER           :: CALIBMOD       ! Calibration options
INTEGER           :: PLOTMOD
INTEGER           :: EGROWTH
CHARACTER(LEN=14) :: fecha          ! Data and time of execution
CHARACTER(LEN=2)  :: day0
CHARACTER(LEN=2)  :: month0
CHARACTER(LEN=2)  :: year0
CHARACTER(LEN=2)  :: hour0
CHARACTER(LEN=2)  :: minute0
REAL(rp)          :: timing0
REAL(rp)          :: timing1

! MODEL PARAMETERS *************************************************************
REAL(rp) :: beta            ! discount rate
REAL(rp) :: irate           ! interest rate
REAL(rp) :: delta0          ! exogenous probability of exit
REAL(rp) :: fc              ! firing cost
REAL(rp) :: gamma           ! decreasing returns to scale parameter
REAL(rp) :: bigA            ! aggregate productivity term
REAL(rp) :: sigma_z         ! standard deviaiton of productivity shocks, benchmarks distribution
REAL(rp) :: rho_0           ! bechanmark probability of innovation
REAL(rp) :: mu_z            ! productivity trend, bechmark distribution
REAL(rp) :: sigma0_z        ! standard eviation of initial draw of productiviy
REAL(rp) :: kappa_0         ! innovation cost, level
REAL(rp) :: kappa_1         ! innovation cost, shape
REAL(rp) :: theta           ! disutility of labor supply
REAL(rp) :: wrate           ! wage rate
REAL(rp) :: psi             ! elasticity of labor supply
REAL(rp) :: psic            ! intertemporal elasticity of substitution

! STATE VARIABLES **************************************************************
REAL(rp) , ALLOCATABLE :: dgrid(:)  ! Grid for productivity
REAL(rp) , ALLOCATABLE :: ngrid(:)  ! Grid for permanent labor force
REAL(rp) , ALLOCATABLE :: Td0(:,:)  ! Transition matrix for productivity
REAL(rp) , ALLOCATABLE :: Td1(:,:)  ! Distribution of initial productivities
REAL(rp) , ALLOCATABLE :: T0(:)     ! Distribution of initial productivities

! FIRMS PROBLEM ****************************************************************
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: V
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pnp
INTEGER  , ALLOCATABLE , DIMENSION(:,:)   :: pin
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: ppn
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pnh
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pnf
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: phr
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pfr
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: phs
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pfs
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pre
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: ppi
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: prho
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: prho0
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pcrho
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pcrho0
REAL(rp) , ALLOCATABLE , DIMENSION(:,:,:) :: pdd
REAL(rp) , ALLOCATABLE , DIMENSION(:,:,:) :: pdd0
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pcpi
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pcpi0
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: px
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: pe
REAL(rp) , ALLOCATABLE , DIMENSION(:,:,:) :: ppd
REAL(rp) , ALLOCATABLE , DIMENSION(:,:)   :: Dist

! AGGREGATES *******************************************************************
REAL(rp) :: TFPR
REAL(rp) :: TFPL
REAL(rp) :: TFP
REAL(rp) :: TOTL
REAL(rp) :: TOTLI
REAL(rp) :: TOTPI
REAL(rp) :: TOTH
REAL(rp) :: TOTF
REAL(rp) :: TOTFe
REAL(rp) :: TOTY
REAL(rp) :: TOTK
REAL(rp) :: TOTHF
REAL(rp) :: TOTFF
REAL(rp) :: TOTD
REAL(rp) :: TOTDg
REAL(rp) :: TOTDg0
REAL(rp) :: TOTR
REAL(rp) :: TOTCPI
REAL(rp) :: TOTCRHO
REAL(rp) :: MMPL
REAL(rp) :: VMPL

! CALIBRATION ******************************************************************
INTEGER , PARAMETER :: NOMS = 14
REAL(rp)            :: MOMS(NOMS)            ! Vector with deviations
REAL(rp)            :: MOMSD(NOMS)           ! Vector with moments from the data
REAL(rp)            :: MOMSM(NOMS)           ! Vector with moments from the model
REAL(rp)            :: WMAT(NOMS)            ! Vector with moments weights
CHARACTER(LEN=15)   :: MNOM(NOMS)            ! Vector with moments names
REAL(rp)            :: d_asize   , m_asize   ! Average firm size (data and model)
REAL(rp)            :: d_ncat(6) , m_ncat(6) ! Share of firms in each size category (data and model)
REAL(rp)            :: d_ash     , m_ash     ! Share of hiring firms (data and model)
REAL(rp)            :: d_arh     , m_arh     ! Hiring rate among hiring firms (data and model)
REAL(rp)            :: d_asf     , m_asf     ! Share of firing firms (data and model)
REAL(rp)            :: d_arf     , m_arf     ! Firing rate among firing firms (data and model)
REAL(rp)            :: d_stdn    , m_stdn    ! CV of employment (data and model)
REAL(rp)            :: d_stdn0   , m_stdn0   ! CV of employment of entrants (data and model)
REAL(rp)            :: d_nage0   , m_nage0   ! Avergae firm size of entrants (data and model)
REAL(rp)            :: d_prof    , m_prof    ! Profits share


INTEGER  :: SKIPPOINT

! ******************************************************************************
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SETPARAMETERS( )

  USE omp_lib , ONLY : OMP_GET_MAX_THREADS
  IMPLICIT NONE
  LOGICAL :: ex1,ex2,ex3
  INTEGER :: i,dat(8)

  CALL chdir(TRIM(ADJUSTL(path)))

  Nd = 50
  Nn = 60

  MODELMOD = 0
  CALIBMOD = 0
  EGROWTH  = 0

  psi  = zero ! inverse of the labro supply elasticity. Set to 0 for log utility
  psic = one  ! intertemporal elasticity of substitution. Set to 1 for log utility

  OPEN(unit=1, file="txtfiles/params.txt", action= "read")
    READ (1,*) irate       ! Interest rate
    READ (1,*) delta0      ! Exogenous probability of exit
    READ (1,*) fc          ! Firing costs
    READ (1,*) gamma       ! Degree of returns to scale
    READ (1,*) bigA        ! Aggregate productivity term
    READ (1,*) sigma0_z    ! Std of initial producitvity draw
    READ (1,*) mu_z        ! Autocorrelation of productivity
    READ (1,*) sigma_z     ! Std of productivity shocks
    READ (1,*) rho_0       ! Benchmark probability of innovation
    READ (1,*) kappa_0     ! Innovation cost, size
    READ (1,*) kappa_1     ! Innovation cost, size
    READ (1,*) theta       ! (Dis)utility of labor, PWs
  CLOSE (1)

  beta    = 1.00/DBLE(1.00+irate)
  threads = OMP_GET_MAX_THREADS()

  CALL DATE_AND_TIME(VALUES=dat)

                    WRITE ( year0   , '(I2)'    )     dat(1)-2000
  IF (dat(2).LT.10) WRITE ( month0  , '(A1,I1)' ) '0',dat(2)
  IF (dat(2).GE.10) WRITE ( month0  , '(I2)'    )     dat(2)
  IF (dat(3).LT.10) WRITE ( day0    , '(A1,I1)' ) '0',dat(3)
  IF (dat(3).GE.10) WRITE ( day0    , '(I2)'    )     dat(3)
  IF (dat(5).LT.10) WRITE ( hour0   , '(A1,I1)' ) '0',dat(5)
  IF (dat(5).GE.10) WRITE ( hour0   , '(I2)'    )     dat(5)
  IF (dat(6).LT.10) WRITE ( minute0 , '(A1,I1)' ) '0',dat(6)
  IF (dat(6).GE.10) WRITE ( minute0 , '(I2)'    )     dat(6)

  WRITE (fecha,'(A2,A2,A2,A1,A2,A2)') year0,month0,day0,'_',hour0,minute0

  CALL ALLOCATESTATE( )

  RETURN
END SUBROUTINE SETPARAMETERS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION FINDTHETA(LAB,CONS,WAGE) RESULT(THETANEW)
  IMPLICIT NONE
  REAL(rp) :: LAB,CONS,WAGE,THETANEW
  THETANEW = WAGE/((CONS**psic)*(LAB**psi))
  RETURN
END FUNCTION FINDTHETA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ALLOCATESTATE( )

  USE toolkit , ONLY : GRID , TAUCHEN
  IMPLICIT NONE
  REAL(rp) :: P0(Nd,Nd)
  INTEGER  :: i

  IF (ALLOCATED(ngrid))  DEALLOCATE(ngrid)  ; ALLOCATE(ngrid(Nn))      ; ngrid  = DBLE(0.00)
  IF (ALLOCATED(dgrid))  DEALLOCATE(dgrid)  ; ALLOCATE(dgrid(Nd))      ; dgrid  = DBLE(0.00)
  IF (ALLOCATED(Td0))    DEALLOCATE(Td0)    ; ALLOCATE(Td0(Nd,Nd))     ; Td0    = DBLE(0.00)
  IF (ALLOCATED(Td1))    DEALLOCATE(Td1)    ; ALLOCATE(Td1(Nd,Nd))     ; Td1    = DBLE(0.00)
  IF (ALLOCATED(T0))     DEALLOCATE(T0)     ; ALLOCATE(T0(Nd))         ; T0     = DBLE(0.00)
  IF (ALLOCATED(V))      DEALLOCATE(V)      ; ALLOCATE(V(Nd,Nn))       ; V      = DBLE(0.00)
  IF (ALLOCATED(pnp))    DEALLOCATE(pnp)    ; ALLOCATE(pnp(Nd,Nn))     ; pnp    = DBLE(0.00)
  IF (ALLOCATED(ppn))    DEALLOCATE(ppn)    ; ALLOCATE(ppn(Nd,Nn))     ; ppn    = DBLE(0.00)
  IF (ALLOCATED(pin))    DEALLOCATE(pin)    ; ALLOCATE(pin(Nd,Nn))     ; pin    = 1
  IF (ALLOCATED(pnh))    DEALLOCATE(pnh)    ; ALLOCATE(pnh(Nd,Nn))     ; pnh    = DBLE(0.00)
  IF (ALLOCATED(pnf))    DEALLOCATE(pnf)    ; ALLOCATE(pnf(Nd,Nn))     ; pnf    = DBLE(0.00)
  IF (ALLOCATED(phr))    DEALLOCATE(phr)    ; ALLOCATE(phr(Nd,Nn))     ; phr    = DBLE(0.00)
  IF (ALLOCATED(pfr))    DEALLOCATE(pfr)    ; ALLOCATE(pfr(Nd,Nn))     ; pfr    = DBLE(0.00)
  IF (ALLOCATED(phs))    DEALLOCATE(phs)    ; ALLOCATE(phs(Nd,Nn))     ; phs    = DBLE(0.00)
  IF (ALLOCATED(pfs))    DEALLOCATE(pfs)    ; ALLOCATE(pfs(Nd,Nn))     ; pfs    = DBLE(0.00)
  IF (ALLOCATED(ppi))    DEALLOCATE(ppi)    ; ALLOCATE(ppi(Nd,Nn))     ; ppi    = DBLE(0.00)
  IF (ALLOCATED(pre))    DEALLOCATE(pre)    ; ALLOCATE(pre(Nd,Nn))     ; pre    = DBLE(0.00)
  IF (ALLOCATED(prho))   DEALLOCATE(prho)   ; ALLOCATE(prho(Nd,Nn))    ; prho   = DBLE(0.00)
  IF (ALLOCATED(pcrho))  DEALLOCATE(pcrho)  ; ALLOCATE(pcrho(Nd,Nn))   ; pcrho  = DBLE(0.00)
  IF (ALLOCATED(pdd))    DEALLOCATE(pdd)    ; ALLOCATE(pdd(Nd,Nn,Nd))  ; pdd    = DBLE(0.00)
  IF (ALLOCATED(pcpi))   DEALLOCATE(pcpi)   ; ALLOCATE(pcpi(Nd,Nn))    ; pcpi   = DBLE(0.00)
  IF (ALLOCATED(px))     DEALLOCATE(px)     ; ALLOCATE(px(Nd,Nn))      ; px     = DBLE(0.00)
  IF (ALLOCATED(pe))     DEALLOCATE(pe)     ; ALLOCATE(pe(Nd,Nn))      ; pe     = DBLE(0.00)
  IF (ALLOCATED(ppd))    DEALLOCATE(ppd)    ; ALLOCATE(ppd(Nd,Nn,Nd))  ; ppd    = DBLE(0.00)
  IF (ALLOCATED(Dist))   DEALLOCATE(Dist)   ; ALLOCATE(Dist(Nd,Nn))    ; Dist   = DBLE(0.00)

  ngrid = GRID( DBLE(Nn)  ,  DBLE(0.0) , Nn+1 , one )
  dgrid = GRID( LOG(ngrid(Nn)*DBLE(1.20)/((bigA*gamma)**(one/(one-gamma)))) , -DBLE(3.0) , Nd , one )
  dgrid = GRID( DBLE(6.0) , -DBLE(3.0) , Nd , one )

  CALL TAUCHEN( dgrid , one  , -mu_z-half*sigma_z*sigma_z , sigma_z  , Nd , Td0 )
  CALL TAUCHEN( dgrid , one  ,      -half*sigma_z*sigma_z , sigma_z  , Nd , Td1 )
  CALL TAUCHEN( dgrid , zero ,    -half*sigma0_z*sigma0_z , sigma0_z , Nd , P0  ) ; T0 = P0(1,:)

  RETURN
END SUBROUTINE ALLOCATESTATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION NCAT(size) RESULT(sizecat)
  IMPLICIT NONE
  INTEGER :: size,sizecat
  IF ( ngrid(size).le.DBLE( 5)  .and. ngrid(size).ge.DBLE( 0) ) sizecat = 1
  IF ( ngrid(size).le.DBLE(10)  .and. ngrid(size).ge.DBLE( 6) ) sizecat = 2
  IF ( ngrid(size).le.DBLE(15)  .and. ngrid(size).ge.DBLE(11) ) sizecat = 3
  IF ( ngrid(size).le.DBLE(20)  .and. ngrid(size).ge.DBLE(16) ) sizecat = 4
  IF ( ngrid(size).le.DBLE(25)  .and. ngrid(size).ge.DBLE(21) ) sizecat = 5
  IF ( ngrid(size).le.ngrid(Nn) .and. ngrid(size).ge.DBLE(26) ) sizecat = 6
  RETURN
END FUNCTION NCAT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE PRINTHEADER(IND0)
  IMPLICIT NONE
  INTEGER , INTENT(IN) , OPTIONAL :: IND0
  INTEGER :: IND
  IND = 0 ; IF (PRESENT(IND0) .AND. IND0.eq.1) IND  = 1
  IF (IND.eq.0) CALL SYSTEM('clear')
  WRITE (IND,13) '  '
  WRITE (IND,13) '  '
  WRITE (IND,13) '  ***************************************************************'
  WRITE (IND,13) '  AGGREGATE EFFECTS OF FIRING COSTS WITH END. PRODUCTIVITY GROWTH'
  WRITE (IND,13) '  Borja Petit, CEMFI.                                            '
  WRITE (IND,13) '  2019                                                           '
  WRITE (IND,13) '  ***************************************************************'
  WRITE (IND,13) ' '
  WRITE (IND,13) '  Date:     ',day0,'/',month0,'/',year0
  WRITE (IND,12) '  Threads:  ',threads
  WRITE (IND,13) ' '
  12 FORMAT (A,I3)
  13 FORMAT (A,A2,A1,A2,A1,A2,A2,A2,A1,A2)
  RETURN
END SUBROUTINE PRINTHEADER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WRITE_RESULTS(IND0,NAME0)
  IMPLICIT NONE
  INTEGER                                  :: j,i,IND
  INTEGER          , INTENT(IN) , OPTIONAL :: IND0
  CHARACTER(LEN=*) , INTENT(IN) , OPTIONAL :: NAME0
  REAL(rp) :: x1,x2
  CHARACTER(LEN=100) :: NAME

  IND = 0            ; IF (PRESENT(IND0) .AND. IND0.eq.1) IND  = 2
  NAME = "calib.txt" ; IF (PRESENT(NAME0)               ) NAME = NAME0(:)
  IF (IND.eq.2) OPEN(unit=2,file='txtfiles/'//TRIM(ADJUSTL(fecha))//'_'//TRIM(ADJUSTL(NAME)), action= "write")
  IF (IND.eq.0) CALL PRINTHEADER(IND)

    WRITE (IND, 89) MODELMOD    ,' ! Solution option                            '
    WRITE (IND, 89) CALIBMOD    ,' ! Calibration option                         '
    WRITE (IND, 91) '                                                           '
    WRITE (IND, 90) irate       ,' ! Interest rate                              '
    WRITE (IND, 90) delta0      ,' ! Exit probability                           '
    WRITE (IND, 90) fc          ,' ! Firing cost                                '
    WRITE (IND, 90) gamma       ,' ! Returns to scale                           '
    WRITE (IND, 90) bigA        ,' ! Aggregate productivity                     '
    WRITE (IND, 90) sigma0_z    ,' ! Std of initial productivity draw           '
    WRITE (IND, 90) mu_z        ,' ! Productivity trend                         '
    WRITE (IND, 90) sigma_z     ,' ! Std of productivity shocks                 '
    WRITE (IND, 90) rho_0       ,' ! Default probability of innovation          '
    WRITE (IND, 90) kappa_0     ,' ! Innovation cost, size                      '
    WRITE (IND, 90) kappa_1     ,' ! Innovation cost, shape                     '
    WRITE (IND, 90) theta       ,' ! (Dis)utility of labor                      '
    WRITE (IND, 99) ' '
    WRITE (IND, 99) ' '

    WRITE (IND, 99) ' ',' ',('*',j=1,32)
    WRITE (IND, 91) '  AGGREGATE FIGURES '
    WRITE (IND, 99) ' ',' ',('*',j=1,32)
    WRITE (IND, 91) '  Wage        ' , wrate
    WRITE (IND, 91) '  Ave. D      ' , TOTD
    WRITE (IND, 91) '  D growth    ' , TOTDg
    WRITE (IND, 99) ' ',' ',('-',j=1,32)
    WRITE (IND, 91) '  Output      ' , TOTY
    WRITE (IND, 91) '  Employment  ' , TOTL
    WRITE (IND, 91) '  Y/L         ' , TOTY/TOTL , TOTY/(TOTL**gamma)
    WRITE (IND, 91) '  Consumption ' , TOTY - TOTLI
    WRITE (IND, 91) '  TFP         ' , TFP
    WRITE (IND, 91) '  Var MPL     ' , cien*VMPL
    WRITE (IND, 99) ' ',' ',('-',j=1,32)
    WRITE (IND, 91) '  Innovation  ' , TOTLI        , cien*TOTLI/TOTY
    WRITE (IND, 91) '  Intensive   ' , TOTCPI       , cien*TOTCPI/TOTLI
    WRITE (IND, 91) '  Lambda      ' , TOTR*cien    , cien*(TOTR - rho_0)
    WRITE (IND, 99) ' ',' ',('-',j=1,32)
    WRITE (IND, 91) '  Hiring      ' , cien*TOTH/TOTL , cien*TOTHF
    WRITE (IND, 91) '  Firing      ' , cien*TOTF/TOTL , cien*TOTFF
    WRITE (IND, 99) ' ',' ',('*',j=1,32)
    WRITE (IND, 99) ' '
    WRITE (IND, 99) ' '

    WRITE (IND, 99) ' ',' ',('*',j=1,49)
    WRITE (IND, 95) '  MOMENT       ','MODEL','DATA','DESV','MOM'
    WRITE (IND, 99) ' ',' ',('*',j=1,49)
    DO i=1,SIZE(MOMS)
      WRITE (IND, 94) MNOM(i) , MOMSM(I) , MOMSD(I) , MOMS(I) , MOMS(I)*WMAT(I)
    END Do
    WRITE (IND, 99) ' ',' ',('-',j=1,49)
    WRITE (IND, 94) '  TOTAL ERROR  ', SUM(MOMS(:)*WMAT(:)*WMAT(:)*MOMS(:))
    WRITE (IND, 99) ' ',' ',('*',j=1,49)
    WRITE (IND, 99) ' '
    WRITE (IND, 99) ' '

  IF (IND.eq.2) CLOSE(2)
  89 FORMAT (I12,A50)
  90 FORMAT (F12.8,A50)
  91 FORMAT (A,3(F10.4))
  94 FORMAT (A,9(F9.2))
  95 FORMAT (A,9(A9))
  99 FORMAT (130(A1))

  RETURN
END SUBROUTINE WRITE_RESULTS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE PRINTGRAPHS( )
  USE toolkit , ONLY : STATS , OLS1VAR
  IMPLICIT NONE
  INTEGER  :: id,in,iid,id25,id75,j
  REAL(rp) :: aux1,aux2,aux3,aux4,PP(Nd,Nd),TTD(Nd),TTN(Nn),DNMAT(Nd,Nn),DN(Nn),auxv(Nd)
  REAL(rp) :: TTD1(Nd),TTD2(Nd),TTD3(Nd)
  REAL(rp) :: TTN1(Nn),TTN2(Nn),TTN3(Nn)

  ! ----------------------------------------------------------------------------
  ! Print sample of firms to file
  !OPEN(unit=4, file ="results/baseline_sample.txt", action="write")
  !  WRITE(4,'(50(A10,A1))') 'd' , ',','n', ',' , 'w' , ',', 'np'
  !  DO id = 1,Nd ; DO in = 1,Nn
  !    WRITE(4,'(100(F10.4,A1))') dgrid(id), ',' , ngrid(in), ',' , mil*Dist(id,in), ',' , pnp(id,in)
  !  END DO ; END DO
  !CLOSE(4)

  ! ----------------------------------------------------------------------------
  ! Print sample of firms to file
  OPEN(unit=4, file ="results/baseline_sample.txt", action="write")
    WRITE(4,'(50(A10,A1))') 'd' , ',' , 'dp' , ',' , 'n' , ',' , 'w' , ',' , 'np'
    DO id = 1,Nd ; DO iid=1,Nd ; DO in = 1,Nn
      WRITE(4,'(100(F10.4,A1))') dgrid(id), ',' , dgrid(iid) , ',', ngrid(in), ',' , &
      mil*ppd(id,in,iid)*Dist(id,in) , ',' , pnp(id,in)
    END DO ; END DO ; END DO
  CLOSE(4)

  ! ----------------------------------------------------------------------------
  ! Distribution of firms by productivity, and by size
  DO id = 1,Nd
    TTD(id) = zero
    DO in = 1,Nn
      TTD(id) = TTD(id) + (SUM(ppd(id,in,:)*dgrid(:)) - dgrid(id))*Dist(id,in)
    END DO
    IF (SUM(Dist(id,:)).gt.tol) TTD(id) = TTD(id)/SUM(Dist(id,:))
  END DO
  DO in = 1,Nn
    TTN(in) = zero
    DO id = 1,Nn
      TTN(in) = TTN(in) + (SUM(ppd(id,in,:)*dgrid(:)) - dgrid(id))*Dist(id,in)
    END DO
    IF (SUM(Dist(:,in)).gt.tol) TTN(in) = TTN(in)/SUM(Dist(:,in))
  END DO

  OPEN(unit=4, file ="results/baseline_dist_d.txt", action="write")
    WRITE(4,'(50(A10))') 'd' , 'Dd' , 'g', 'D0'
    DO id = 1,Nd
      WRITE(4,'(100(F10.4))') dgrid(id) , SUM(Dist(id,:)) , TTD(id) , T0(id) , SUM(Dist(id,:)*pnp(id,:))
    END DO
  CLOSE(4)

  OPEN(unit=5, file ="results/baseline_dist_n.txt", action="write")
    WRITE(5,'(50(A10))') 'n', 'Dn' , 'g'
    DO in = 1,Nn
      WRITE(5,'(100(F10.4))') ngrid(in) , SUM(Dist(:,in)) , TTN(in) , SUM(Dist(:,in)*EXP(dgrid(:)))
    END DO
  CLOSE(5)
  CALL SYSTEM('gnuplot -p figures/codes/baseline_dist_d_n.gnu')

  ! ----------------------------------------------------------------------------
  ! Innovation choices: large vs small firms
  DO id = 1,Nd
    PP(id,:) = zero
    aux1     = zero
    DO in=1,Nn
      PP(id,:) = PP(id,:) + pdd(id,in,:)*Dist(id,in)
    END DO
    PP(id,:) = PP(id,:)/SUM(Dist(id,:))
  END DO
  OPEN(unit=7, file="results/baseline_innovation_large_vs_small.txt", action= "write")
    WRITE(7,'(100(A20))') 'dlow','td0low','td1low','dhigh','td0high','td1high'
    aux1 = zero ; id = 0
    DO WHILE (aux1.lt.0.25)
      id  = id + 1
      aux1 = aux1 + SUM(Dist(id,:))
    END DO
    id25 = id - 1
    aux1 = zero ; id = 0
    DO WHILE (aux1.lt.0.75)
      id  = id + 1
      aux1 = aux1 + SUM(Dist(id,:))
    END DO
    id75 = id - 1
    DO id=1,Nd
      WRITE(7,'(100(F20.15))') dgrid(id)-dgrid(id25) , Td0(id25,id) , PP(id25,id) , &
                               dgrid(id)-dgrid(id75) , Td0(id75,id) , PP(id75,id)
    END DO
  CLOSE(7)
  CALL SYSTEM('gnuplot -p figures/codes/baseline_innovation_large_vs_small.gnu')

  ! ----------------------------------------------------------------------------
  ! Growhth rates by productivity
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(ppd(id,in,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO id=1,Nd
    TTD(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(pdd(id,in,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO id=1,Nd
    TTD1(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(Td0(id,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO id=1,Nd
    TTD2(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  OPEN(unit=8, file ="results/baseline_growth_d_mean.txt", action="write")
    WRITE(8,'(23(A10))') 'd','dist','agd','igd','ngd'
    DO id=1,Nd
      IF (.NOT.ISNAN(TTD(id))) WRITE(8,'(23(F10.4))')  dgrid(id),SUM(Dist(id,:)),TTD(id),TTD1(id),TTD2(id)
    END DO
  CLOSE(8)

  ! ----------------------------------------------------------------------------
  ! Growhth rates by size
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(ppd(id,in,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO in=1,Nn
    TTN(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(pdd(id,in,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO in=1,Nn
    TTN1(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    DNMAT(id,in) = SUM(Td0(id,:)*dgrid(:)) - dgrid(id)
  END DO ; END DO
  DO in=1,Nn
    TTN2(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  OPEN(unit=9, file ="results/baseline_growth_n_mean.txt", action="write")
    WRITE(9,'(23(A10))') 'n','dist','agd','igd','ngd'
    DO in=1,Nn
      IF (.NOT.ISNAN(TTN(in))) WRITE(9,'(23(F10.4))')  ngrid(in),SUM(Dist(:,in)),TTN(in),TTN1(in),TTN2(in)
    END DO
  CLOSE(9)

  ! ----------------------------------------------------------------------------
  ! Std of growhth rates by productivity
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),ppd(id,in,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO id=1,Nd
    TTD(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),pdd(id,in,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO id=1,Nd
    TTD1(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),Td0(id,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO id=1,Nd
    TTD2(id) = SUM(Dist(id,:)*DNMAT(id,:))/SUM(Dist(id,:))
  END DO
  OPEN(unit=10, file ="results/baseline_growth_d_sd.txt", action="write")
    WRITE(10,'(23(A10))') 'd','dist','agd','igd','ngd'
    DO id=1,Nd
      IF (.NOT.ISNAN(TTD(id))) WRITE(10,'(23(F10.4))')  dgrid(id),SUM(Dist(id,:)),TTD(id),TTD1(id),TTD2(id)
    END DO
  CLOSE(10)

  ! ----------------------------------------------------------------------------
  ! Std of growhth rates by size
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),ppd(id,in,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO in=1,Nn
    TTN(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),pdd(id,in,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO in=1,Nn
    TTN1(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  DO in=1,Nn ; DO id=1,Nd
    CALL STATS(dgrid(:)-dgrid(id),Td0(id,:),aux1,DNMAT(id,in))
  END DO ; END DO
  DO in=1,Nn
    TTN2(in) = SUM(Dist(:,in)*DNMAT(:,in))/SUM(Dist(:,in))
  END DO
  OPEN(unit=11, file ="results/baseline_growth_n_sd.txt", action="write")
    WRITE(11,'(23(A10))') 'n','dist','agd','igd','ngd'
    DO in=1,Nn
      IF (.NOT.ISNAN(TTN(in))) WRITE(11,'(23(F10.4))')  ngrid(in),SUM(Dist(:,in)),TTN(in),TTN1(in),TTN2(in)
    END DO
  CLOSE(11)

  ! ----------------------------------------------------------------------------
  ! Prob. of innovation by productivity
  OPEN(unit=12, file ="results/baseline_growth_d_rho.txt", action="write")
    WRITE(12,'(23(A10))') 'd','dist','rho','rho0'
    DO id=1,Nd
      IF (.NOT.ISNAN(TTD(id))) WRITE(12,'(23(F10.4))') &
      dgrid(id),SUM(Dist(id,:)),SUM(Dist(id,:)*prho(id,:))/SUM(Dist(id,:)),rho_0
    END DO
  CLOSE(12)

  ! ----------------------------------------------------------------------------
  ! Prob. of innovation by size
  OPEN(unit=13, file ="results/baseline_growth_n_rho.txt", action="write")
    WRITE(13,'(23(A10))') 'n','dist','rho','rho0'
    DO in=1,Nn
      IF (.NOT.ISNAN(TTN(in))) WRITE(13,'(23(F10.4))') &
      ngrid(in),SUM(Dist(:,in)),SUM(Dist(:,in)*prho(:,in))/SUM(Dist(:,in)),rho_0
    END DO
  CLOSE(13)

  ! ----------------------------------------------------------------------------
  ! Generate graphs
  CALL SYSTEM('gnuplot -p figures/codes/baseline_growth_d.gnu')
  CALL SYSTEM('gnuplot -p figures/codes/baseline_growth_n.gnu')

  PP(:,:) = zero
  DO id = 1,Nd
    PP(id,:) = zero
    DO in=1,Nn
      PP(id,:) = PP(id,:) + pdd(id,in,:)*Dist(id,in)
    END DO
    PP(id,:) = PP(id,:)/SUM(Dist(id,:))
  END DO

  ! ----------------------------------------------------------------------------
  ! Distribution of productivity
  OPEN(unit=13, file ="results/baseline_transitionmatrix.txt", action="write")
    DO id=1,Nd
      WRITE(13,'(100(F10.4))') PP(id,:)
    END DO
  CLOSE(13)

  ! ----------------------------------------------------------------------------
  ! Distribution of productivity
  OPEN(unit=13, file ="results/baseline_expenses_size.txt", action="write")
  WRITE(13,'(100(A10))') 'Lab','Exp','Rev','Share'
    DO in=1,Nn
      IF (SUM(Dist(:,in)).gt.tol) WRITE(13,'(100(F10.4))') ngrid(in) , &
      SUM(Dist(:,in)*px(:,in))/SUM(Dist(:,in)) , &
      SUM(Dist(:,in)*pre(:,in))/SUM(Dist(:,in)) , &
      cien*SUM(Dist(:,in)*px(:,in))/SUM(Dist(:,in)*pre(:,in))
      IF (SUM(Dist(:,in)).le.tol) WRITE(13,'(100(F10.4))') ngrid(in),zero,zero,zero
    END DO
  CLOSE(13)

  RETURN
END SUBROUTINE PRINTGRAPHS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE PARAMS
