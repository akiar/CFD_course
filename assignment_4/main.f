*
*        file main.f
************************************************************************************
*
*     Course:  MME 9710 - Advanced CFD
*
*     Mainline for computations of one-dimensional diffusion 
*     problems.
*
************************************************************************************
*
*==============================
*     DECLARATION OF VARIABLES
*==============================
*
      PARAMETER(ID=120)
      REAL YNE(ID),ZUE(ID),XP(ID),XE(ID)
      REAL DIEP(ID),DISW(ID),DISE(ID),DSXY(ID),DSXZ(ID)
      REAL AREP(ID),ARO(ID),VOLP(ID)
      REAL WORK1(ID),WORK2(ID)
      REAL DE(ID),QT(ID),RT(ID)
      REAL ATW(ID),ATE(ID),ATP(ID),BT(ID)
      REAL T(ID),TOLD(ID),DEOLD(ID)
*
      INTEGER I,KNTIN,IE1
      INTEGER IDATI,IRSI,IDATO,IRSO,ITERMO
      INTEGER IB,IE,IDTYP,KNTTM
      INTEGER LVLGEO,LVLCOF,KNTOUT,KNTNL
      REAL DI,RHO,COND,CP,EMIS,VISC,T0,DTIME,CRIT
      REAL HCONV,TINF
*
      INTEGER ADVSCM                      ! Advection scheme
      REAL RSD(ID),AVRSD,RESIDUALS(ID)    ! Residual variables
      REAL U(ID),U0,UHE(ID),UHE0          ! New initialization variables
      REAL ACUE(ID),ACUW(ID),BC(ID)       ! for flow field
      REAL DCCE(ID)                       ! deferred correction array
      REAL ME(ID),ALFAE(ID)               ! Advection properties
*
      REAL DPDX(ID),P(ID),POLD(ID),P0
      REAL QU(ID),RU(ID)
      REAL AUP(2,2,ID),AUW(2,2,ID),AUE(2,2,ID),BU(2,ID)
      REAL DHUE(ID),UOLD(ID)
      REAL RESIDC,RESIDP,CD
      INTEGER PRBLMT
*
*==============================
*     INPUT AND INITIALIZATION
*==============================
*
*     READ INPUT PARAMETERS
*
      PRINT *, "*******************************************************
     C          *******************************************************"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "--------------------NEW SOLUTION-----------------------"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "*"
      PRINT *, "*******************************************************
     C          *******************************************************"   
      CALL FILDEF(IDATI,IRSI,IDATO,IRSO,ITERMO)
      CALL INPUT(IB,IE,IDTYP,DI,
     C           RHO,COND,CP,VISC,EMIS,
     C           T0,DTIME,KNTTM,KNTNL,CRIT,
     C           LVLGEO,LVLCOF,HCONV,TINF,IDATI,
     C           OMEG,
     C           U0,UHE0,ADVSCM,
     C           P0,PRBLMT)
*
      CALL ECHO(IB,IE,IDTYP,DI,
     C          RHO,COND,CP,VISC,EMIS,
     C          T0,DTIME,KNTTM,KNTNL,CRIT,
     C          LVLGEO,LVLCOF,HCONV,TINF,IDATO,IDATI,
     C          OMEG,
     C          U,U0,UHE,UHE0)
*
*     DEFINE INTEGRATION POINTS AT C.V. CORNERS
*
      CALL MAKGRD(XE,YNE,ZUE, DI,IDTYP,
     C            IB,IE,ID,ITERMO,IDATO)
*
*     COMPUTE GEOMETRY
*
      CALL GRDGEO(XP,DIEP,DISE,DISW,AREP,ARO,VOLP,
     C     XE,YNE,ZUE,IDTYP,DSXY,DSXZ,IB,IE,ID,IDATO,ITERMO,LVLGEO)
*
*     CHECK INTEGRATION POINTS AND GEOMETRY (set LVLGEO=1 to check)
*
      IF( LVLGEO .GE. 1) THEN
         CALL OUT1D(XE,  ' XE     ',IDATO,IB-1,IE+1,1,ID)
         CALL OUT1D(YNE, ' YNE    ',IDATO,IB-1,IE+1,1,ID)
         CALL OUT1D(ZUE, ' ZUE    ',IDATO,IB-1,IE+1,1,ID)
         CALL OUT1D(XP,  ' XP     ',IDATO,IB-1,IE+1,1,ID)
         CALL OUT1D(DIEP,' DIEP   ',IDATO,IB-1,IE  ,1,ID)
         CALL OUT1D(DISE,' DISE   ',IDATO,IB-1,IE  ,1,ID)
         CALL OUT1D(AREP,' AREP   ',IDATO,IB-1,IE+1,1,ID)
         CALL OUT1D(ARO, ' ARO    ',IDATO,IB  ,IE  ,1,ID)
         CALL OUT1D(VOLP,' VOLP   ',IDATO,IB-1,IE+1,1,ID)
      END IF
*
*     INITIALIZE T AND DE fields
*
      CALL INITAL(T,U,UHE,P, T0,U0,UHE0,P0,IRSI,IB,IE,ID)
*     
      KNTOUT = 0
      IF (PROBLMT == 1) THEN
        CALL DIFPHI(DE, COND/CP,AREP,DIEP,IB,IE,ID)                ! CHANGE VISC TO COND/CP
        CALL OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T,XP,KNTOUT)
        CALL OUT1D(T  ,'T(init)',IDATO,IB-1,IE+1,1,ID)
      ELSE
        CALL DIFPHI(DE, VISC,AREP,DIEP,IB,IE,ID)
        CALL OUTPYM(ID,IB,IE,DE,ATW,ATP,BT,P,U,XP,KNTOUT)
        CALL OUT1D(P  ,'P(init)',IDATO,IB-1,IE+1,1,ID)
      ENDIF
      PRINT *, "INITAL"
*
*==============================
*     SOLVER
*==============================
*
      DO 2000 KNTOUT=1,KNTTM                                  ! TIME LOOP
*      
         DO 100 I=IB-1,IE+1                                   ! SAVE OLD SOLUTIONS
            TOLD(I) = T(I)
            DEOLD(I) = DE(I)
            POLD(I) = P(I)
            UOLD(I) = U(I)
 100     CONTINUE
*
         DO 1700 M=1,KNTNL                                    ! NON LINEAR LOOP
           PRINT *, "==========================="
           PRINT *, "==========================="
*
*=================================================
*     --   MODULARIZE FOR PROBLEM TYPE
*            1: Temperature only
*            2: Momentum only
*            3: Combined temperature and momentum
*=================================================
*
*********************************************
*     --   1 OR 3: TEMPERATURE OR COMBINED
*
           IF (PRBLMT==1 .OR. PRBLMT==3) THEN
             PRINT *, "TEMP OR COMBINED"
*
*********************************************
*
*     ---    COMPUTE DIFFUSION COEFFICIENTS
*
             CALL NULL(BT, IB,IE,ID)                            ! NULL COEFFICIENTS
             CALL DIFPHI(DE, COND/CP,AREP,DIEP,IB,IE,ID)  
             PRINT *, "DIFPHI"
*
*     ---    COMPUTE ADVECTION CORRECTIONS
*
             CALL HOCONV(DCCE, ADVSCM,IB,IE,ID,T,RHO,  
     C                   ME,ALFAE,XP,XE)
             PRINT *, "HOCONV"
*
*     ---    COMPUTE SOURCE TERMS & ACTIVE COEFFICIENTS
*
             CALL SRCT(QT,RT, T,VOLP,ARO,HCONV/CP,TINF,IB,IE,ID,
     C                 EMIS,
     C                 OMEG,DEOLD,TOLD,
     C                 DCCE)
             PRINT *, "SRCT"
*
             CALL COEFF(ATP,ATW,ATE,BT,
     C                  ME,DE,QT,RT,
     C                  TOLD,VOLP,ALFAE,RHO,DTIME,
     C                  IB,IE,ID,OMEG)  
             PRINT *, "COEFF"
*
*     ---    SET BOUNDARY CONDITIONS
*
             CALL BNDCT(ATP,ATW,ATE,BT, DE,AREP,IB,IE,ID,HCONV,TINF)
             PRINT *, "BNDCT"
*
*     ---    CHECK COMPUTED, ACTIVE COEFFICIENTS (LVLCOF=1 to check)
*
             IF (LVLCOF .GE. 1) THEN
               CALL OUT1D(DE,   ' DE     ',IDATO,IB-1, IE  ,1,ID)
               CALL OUT1D(ATW , ' ATW    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUT1D(ATP , ' ATP    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUT1D(ATE , ' ATE    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUT1D(BT  , ' BT     ',IDATO,IB-1, IE+1,1,ID)
             ENDIF
*
*     ---    CALCULATE RESIDUALS
*
             IF (M > 1) THEN
               CALL RESID(AVRSD,RSD, T,ATP,ATW,ATE,BT,IB,IE,ID)       ! TEMPERATURE RESIDUAL
               RESIDUALS(M) = AVRSD                                   ! Save all residual calculations
               PRINT *, "AVERAGE RESIDUAL ",M," = ",AVRSD             ! Print avg. residual
*
*     ----     CHECK CONVERGENCE
*
               IF (ABS(AVRSD) < CRIT) THEN
                 PRINT *, "CONVERGENCE"
                 EXIT
*
               END IF
             END IF
*
*     ---    COMPUTE TEMPERATURE SOLUTION USING DIRECT SOLVER 
*
             CALL TDMA(T, ATP,ATE,ATW,BT,IB-1,IE+1,WORK1,WORK2,ID) 
             PRINT *, "TDMA"
*
*     ---    PRINT SOLUTION
*
             CALL OUT1D(T  , ' T      ',IDATO,IB-1, IE+1,1,ID)
*
*************************************************
*     --   2 OR 3: MOMENTUM OR COMBINED
*
           ELSEIF (PRBLMT==2 .OR. PRBLMT==3) THEN
             PRINT *, "MOMENTUM OR COMBINED"
*
*************************************************
*
             PRINT *, "RADIUS"
             DO 11 I=IB-1,IE+1
               PRINT *, I,YNE(I)
 11	         CONTINUE
             CALL NULLVB(BU, IB,IE,1,ID)                             ! TD: N=2 CORRECT? NULL VECTOR
             CALL NULLVB(BU, IB,IE,2,ID)                             ! TD: N=2 CORRECT? NULL VECTOR
*
*     --   COMPUTE PRESSURE GRADIENTS, UHAT & MASS FLUX
*
             CALL COEFCN(ACUW,ACUE,BC, RHO,AREP,IB,IE,ID)
             CALL MASFLX(ME, UHE,ACUW,ACUE,BC,IB,IE,ID)
             PRINT *, "COEFCN & MASFLX"
             CALL WEIGHT(ALFAE, ME,DE,IB,IE,ID)
             PRINT *, "WEIGHT"
*
*     ---    COMPUTE PRESSURE GRADIENT
*
             CALL GRADP(DPDX, P,DIEP,IB,IE,ID)
             PRINT *, "GRADP"
*
*     ---    COMPUTE DIFFUSION COEFFICIENTS
*
             CALL DIFPHI(DE, VISC,AREP,DIEP,IB,IE,ID)
             PRINT *, "DIFPHI"
*
*     ---    COMPUTE ADVECTION CORRECTIONS
*
             CALL HOCONV(DCCE, ADVSCM,IB,IE,ID,U,RHO,
     C                   ME,ALFAE,XP,XE)  
             PRINT *, "HOCONV"
*
*     ---    COMPUTE SOURCE TERMS & ACTIVE COEFFICIENTS
*
             CALL SRCU(QU,RU, U,UHE,RHO,VISC,XE,ARO,AREP,IB,IE,ID,
     C                 YNE,ZNE,DCCE)
             PRINT *, "SRCU"
*
             CALL SRCUP(AUP,AUW,AUE,BU, DIEP,VOLP,IB,IE,ID)
             PRINT *, "SRCUP"
*
             CALL COEFFM(AUP,AUW,AUE,BU,
     C                   ME,DE,QU,RU,UOLD,VOLP,ALFAE,RHO,DTIME,
     C                   IB,IE,ID)
             PRINT *, "COEFFM"
*
*     ---    COMPUTE DHAT & UHAT
*
             CALL DHAT(DHUE, AUP,VOLP,RHO,IB,IE,ID)
             PRINT *, "DHAT"
*
             CALL UHAT(UHE, U,P,DHUE,RHO,DPDX,DIEP,IB,IE,ID)
             PRINT *, "UHAT"
*
*     ---    COMPUTE TERMS IN CONSERVATION OF MASS EQUATION
*
             CALL ADCONT(AUP,AUW,AUE,BU,ACUW,ACUE,BC,DHUE,DPDX,
     C                   RHO,DIEP,DISE,IB,IE,ID)
             PRINT *, "ADCONT"
*
*     ---    SET BOUNDARY CONDITIONS
*
             CALL BNDCU(AUP,AUW,AUE,BU, IB,IE,ID)
             PRINT *, "BNDCU"
*
             CALL BNDCP(AUP,AUW,AUE,BU, P,XP,IB,IE,ID,DE)   !DIEP TO XP
             PRINT *, "BNDCP"
*
*     ---    CHECK COMPUTED, ACTIVE COEFFICIENTS (LVLCOF=1 to check)
*
             IF (LVLCOF .GE. 1) THEN
               CALL OUTCON(AUW , ' AUW    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUTCON(AUP , ' AUP    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUTCON(AUE , ' AUE    ',IDATO,IB-1, IE+1,1,ID)
               CALL OUTCON(BU  , ' BU     ',IDATO,IB-1, IE+1,1,ID)
             ENDIF
*
*     ---    CALCULATE RESIDUALS
*
             IF (M > 1) THEN
               CALL RESIDM(AVRSD,RSD, P,U,1,AUP,AUW,AUE,BU,IB,IE,ID)
               RESIDC = AVRSD
               PRINT *, "AVERAGE CONSERVATION RES. ",M," = ",RESIDC           ! Print avg. residual
               CALL RESIDM(AVRSD,RSD, P,U,2,AUP,AUW,AUE,BU,IB,IE,ID)
               RESIDP = AVERSD
               PRINT *, "AVERAGE MOMENTUM RES. ",M," = ",RESIDP               ! Print avg. residual
*
*     ----     CHECK CONVERGENCE
*
               IF (RESIDP < CRIT .AND. RESIDC < CRIT) THEN
                 PRINT *, "CONVERGENCE"
                 EXIT
*
               END IF
             END IF
*
*     ---    COMPUTE MOMENTUM SOLUTION USING DIRECT SOLVER 
*
*             DO 3000 I=IB,IE
*              PRINT *, I,":",AUP(2,1,I),AUW(2,1,I),AUE(2,1,I),BU(2,I)
*              PRINT *, I,":",BU(1,I),BU(2,I)
* 3000        CONTINUE
             CALL SOLUP(AUP,AUW,AUE,BU,P,U,IB,IE,ID)                      ! Updates P & U fields
             PRINT *, "SOLUP"
*             DO 10 I=IB-1,IE+1
*               PRINT *, "I: ",I," P: ",P(I)," U: ",U(I)
* 10          CONTINUE
*
*     ---    PRINT SOLUTION
*
             CALL OUT1D(P  , ' P      ',IDATO,IB-1, IE+1,1,ID)            ! TD: CHECK IF P IS CORRECT
             CALL OUT1D(U  , ' U      ',IDATO,IB-1, IE+1,1,ID)            ! TD: CHECK IF P IS CORRECT
*
*********************************************
*     --   OTHER: INVALID PROBLEM TYPE
*
           ELSE
             PRINT *, "ENTER VALID PROBLEM TYPE"
*
*********************************************
*
           ENDIF
 1700    CONTINUE
*
*=======================
*     FINAL FORMATTING
*=======================
*
*     -  PRINT TO PYTHON FILE AND SAVE RESULT TO UNFORMATTED FILE
*
         IF (PRBLMT==1) THEN
           CALL OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T,XP,KNTOUT,"T")
         ELSE
           CALL OUTPYM(ID,IB,IE,DE,ATW,ATP,BT,P,U,XP,KNTOUT,"P") ! TD: CHECK IF P IS CORRECT OR U
         ENDIF
         CD = (P(IB-1) - P(IE+1))/0.5/RHO/(U(IB-1)**2)
         PRINT *, "CD: ",CD
*
*     -  SAVE TO UNFORMATTED FILE
*
         CALL SAVE(T,U,UHE,P,IRSO,IB,IE,ID)
*
*     -  PRINT TIMESTEP ITERATION
*
         PRINT *, "TIMESTEP: ", KNTOUT, " FINISHED"
         PRINT *, "-----------------"
         PRINT *, "-----------------"
         PRINT *, "-----------------"
*
*********************************************
*
 2000 CONTINUE
      STOP
      END
*