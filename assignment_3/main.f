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
*============================
*  Declaration of Variables
*============================
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
      REAL RSD(ID),AVRSD,RESIDUALS(ID)    !Residual variables
      REAL TEND(ID),ERROR(ID),SUMERR
      REAL U(ID),U0,UHE(ID),UHE0,P(ID),P0 !New initialization variables
      REAL ACUE(ID),ACUW(ID),BC(ID)       ! for flow field
      REAL DCCE(ID)                       ! deferred correction array
      REAL ME(ID),ALFAE(ID)               ! Advection properties
*
*============================
*  Initialization and input
*============================
*
*--Read input parameters
*
      PRINT *, "*********************"
      CALL FILDEF(IDATI,IRSI,IDATO,IRSO,ITERMO)
      CALL INPUT(IB,IE,IDTYP,DI,
     C           RHO,COND,CP,VISC,EMIS,
     C           T0,DTIME,KNTTM,KNTNL,CRIT,
     C           LVLGEO,LVLCOF,HCONV,TINF,IDATI,
     C           OMEG,
     C           U0,UHE0,ADVSCM)
      CALL ECHO(IB,IE,IDTYP,DI,
     C          RHO,COND,CP,VISC,EMIS,
     C          T0,DTIME,KNTTM,KNTNL,CRIT,
     C          LVLGEO,LVLCOF,HCONV,TINF,IDATO,IDATI,
     C          OMEG,
     C          U,U0,UHE,UHE0)
*
*--Define the integration points at the C.V. corners.
*
      CALL MAKGRD(XE,YNE,ZUE, DI,IDTYP,
     C            IB,IE,ID,ITERMO,IDATO)
*
*--Compute geometry parameters.
*
      CALL GRDGEO(XP,DIEP,DISE,DISW,AREP,ARO,VOLP,
     C     XE,YNE,ZUE,IDTYP,DSXY,DSXZ,IB,IE,ID,IDATO,ITERMO,LVLGEO)
*
*--Check integration points and geometry (set LVLGEO=1 to check)
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
*--Initialize T AND DE fields
*
      CALL INITAL(T,U,UHE,P, T0,U0,UHE0,P0,IRSI,IB,IE,ID)
      PRINT *, "INITAL"
*
*--Initialize DE field
*
      CALL DIFPHI(DE, COND/CP,AREP,DIEP,IB,IE,ID)
*
      KNTOUT = 0
      CALL OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T,XP,KNTOUT)
*
*--Print initialized field
*
      CALL OUT1D(T  ,'T(init)',IDATO,IB-1,IE+1,1,ID)
*
*============================
*  Solver
*============================
*
*--Begin time loop
*
      DO 2000 KNTOUT=1,KNTTM
*
*  --Set TOLD and DEOLD
*
         DO 100 I=IB-1,IE+1
            TOLD(I) = T(I)
            DEOLD(I) = DE(I)
 100     CONTINUE
*
*  --Begin non-linear loop
*
         DO 1700 M=1,KNTNL
*
*     --Compute Diffusion coefficients
*
           CALL NULL(BT, IB,IE,ID)
           CALL DIFPHI(DE, COND/CP,AREP,DIEP,IB,IE,ID)
           PRINT *, "DIFPHI"
*
*     --Compute mass flux through east face
*
           CALL COEFCN(ACUW,ACUE,BC, RHO,AREP,IB,IE,ID)
           CALL MASFLX(ME, UHE,ACUW,ACUE,BC,IB,IE,ID)
           PRINT *, "COEFCN & MASFLX"
*
*     --Compute convective weighting factors on east face
*
           CALL WEIGHT(ALFAE, ME,DE,IB,IE,ID)
           PRINT *, "WEIGHT"
*
*     --Compute advection correction terms
*
           CALL HOCONV(DCCE, ADVSCM,IB,IE,ID,DE,T,RHO,  
     C                 ME,ALFAE,XP,XE)
           PRINT *, "HOCONV"
*
*     --Compute Sources and active coefficients
*
           CALL SRCT(QT,RT, T,VOLP,ARO,HCONV/CP,TINF,IB,IE,ID,
     C               EMIS,
     C               OMEG,DEOLD,TOLD,
     C               DCCE)
           PRINT *, "SRCT"
*
           CALL COEFF(ATP,ATW,ATE,BT,
     C                DE,QT,RT,VOLP,RHO,CP,
     C                IB,IE,ID,
     C                OMEG,DTIME,TOLD,
     C                ME,ALFAE,HCONV,ARO,TINF)
           PRINT *, "COEFF"
*
*     --Set boundary conditions
*
           CALL BNDCT(ATP,ATW,ATE,BT, DE,AREP,IB,IE,ID,HCONV,TINF)
           PRINT *, "BNDCT"
*
*     --Check computed, active coefficients (LVLCOF=1 to check)
*
           IF (LVLCOF .GE. 1) THEN
             CALL OUT1D(DE,   ' DE     ',IDATO,IB-1, IE  ,1,ID)
             CALL OUT1D(ATW , ' ATW    ',IDATO,IB-1, IE+1,1,ID)
             CALL OUT1D(ATP , ' ATP    ',IDATO,IB-1, IE+1,1,ID)
             CALL OUT1D(ATE , ' ATE    ',IDATO,IB-1, IE+1,1,ID)
             CALL OUT1D(BT  , ' BT     ',IDATO,IB-1, IE+1,1,ID)
           ENDIF
*
*     --Calculate Residuals and check convergence
*
           IF (M > 1) THEN      ! Only calculate residuals if first loop is passed
*
             CALL RESID(RSD,AVRSD, T,ATP,ATW,ATE,BT,IB,IE,ID) ! use old temperature
             RESIDUALS(M) = ABS(AVRSD)    ! Save all residual calculations
             PRINT *, "AVERAGE RESIDUAL ",M," = ",ABS(AVRSD) ! Print avg. residual
*
          ! Check if convergence has been met, exit loop if it has
*
             IF (ABS(AVRSD) < CRIT) THEN
               PRINT *, "Convergence Reached"
               EXIT
             END IF
           END IF
*
*     --Compute solution using direct solver if convergence is not met 
*
           CALL TDMA(T, ATP,ATE,ATW,BT,IB-1,IE+1,WORK1,WORK2,ID) 
           PRINT *, "TDMA"
*
*     --Print final solution
*
           CALL OUT1D(T  , ' T      ',IDATO,IB-1, IE+1,1,ID)
*
*     --Continue through Non-linearity loop
*
 1700     CONTINUE
*
*  --Print to python file for plotting
*
         CALL OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T,XP,KNTOUT)
*
*  --Save result to unformatted output file
*
         CALL SAVE(T,U,UHE,P,IRSO,IB,IE,ID)
         PRINT *, "Timestep: ", KNTOUT, " Finished"
         PRINT *, "-----------------"
*
*--Continue through Time-Step loop
*
 2000 CONTINUE
*
      STOP
      END
*