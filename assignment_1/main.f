*
*        file main.f
*********************************************************************
*
*     Course:  MME 9710 - Advanced CFD
*
*     Mainline for computations of one-dimensional diffusion 
*     problems.
*
********************************************************************
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
      REAL T(ID),TOLD(ID)
*
      INTEGER I,KNTIN,IE1
      INTEGER IDATI,IRSI,IDATO,IRSO,ITERMO
      INTEGER IB,IE,IDTYP,KNTTM
      INTEGER LVLGEO,LVLCOF,KNTOUT,KNTNL
      REAL DI,RHO,COND,CP,EMIS,VISC,T0,DTIME,CRIT
      REAL HCONV,TINF
      REAL CRITERIA,HEATFLUX,TEMPGRAD
      REAL RSD(ID),AVRSD,RESIDUALS(ID)
      INTEGER CENTERNODE,INTER
*
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
     C          RHO,COND,CP,VISC,EMIS,
     C          T0,DTIME,KNTTM,KNTNL,CRIT,
     C          LVLGEO,LVLCOF,HCONV,TINF,IDATI)
      CALL ECHO(IB,IE,IDTYP,DI,
     C          RHO,COND,CP,VISC,EMIS,
     C          T0,DTIME,KNTTM,KNTNL,CRIT,
     C          LVLGEO,LVLCOF,HCONV,TINF,IDATO)
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
*--Initialize T field
*
      CALL INITAL(T, T0,IRSI,IB,IE,ID)
*     
*--Print initialized field
*
      CALL OUT1D(T  ,'T(init)',IDATO,IB-1,IE+1,1,ID)
*      
*--Begin loop for non-linearities
*
      DO 10 M=1,KNTNL
*
*--Compute active coefficients for T
*
        CALL NULL(BT, IB,IE,ID)
        CALL DIFPHI(DE, COND,AREP,DIEP,IB,IE,ID)
        PRINT *, "Diffusion Coefficients"
*
        CALL SRCT(QT,RT, T,VOLP,ARO,HCONV,TINF,IB,IE,ID,EMIS)
        PRINT *, "Source terms"
*
        CALL COEFF(ATP,ATW,ATE,BT,
     C             DE,QT,RT,VOLP,RHO,CP,
     C             IB,IE,ID)
        PRINT *, "Active Coefficients"
*
*--Set boundary conditions
*
        CALL BNDCT(ATP,ATW,ATE,BT, DE,AREP,IB,IE,ID,HCONV,TINF)
        PRINT *, "Boundary Conditions"
*
*--Check computed, active coefficients (LVLCOF=1 to check)
*
        IF (LVLCOF .GE. 1) THEN
          CALL OUT1D(DE,   ' DE     ',IDATO,IB-1, IE  ,1,ID)
          CALL OUT1D(ATW , ' ATW    ',IDATO,IB-1, IE+1,1,ID)
          CALL OUT1D(ATP , ' ATP    ',IDATO,IB-1, IE+1,1,ID)
          CALL OUT1D(ATE , ' ATE    ',IDATO,IB-1, IE+1,1,ID)
          CALL OUT1D(BT  , ' BT     ',IDATO,IB-1, IE+1,1,ID)
        ENDIF
*
*--Calculate Residuals and check convergence
*
        IF (M > 1) THEN      ! Only calculate residuals if first loop is passed
*
          CALL RESID(RSD,AVRSD, T,ATP,ATW,ATE,BT,IB,IE,ID) ! use old temperature
          RESIDUALS(M) = AVRSD    ! Save all residual calculations
          PRINT *, "AVERAGE RESIDUAL ",M," = ",AVRSD
*
          ! Check if convergence has been met, exit if it has
*
          IF (ABS(AVRSD) < CRIT) THEN
            PRINT *, "Convergence Reached"
            EXIT
          END IF
        END IF
*
*--Compute solution using direct solver 
*
        CALL TDMA(T, ATP,ATE,ATW,BT,IB-1,IE+1,WORK1,WORK2,ID) 
        PRINT *, "Solved"
*
*--Print final solution
*
        CALL OUT1D(T  , ' T      ',IDATO,IB-1, IE+1,1,ID)
*
*--Continue through Non-linearity loop
*
 10	  CONTINUE
*
*--Print to python file for plotting
*
      CALL OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T,XP)
*
*     QUESTION 1 & 4: Calculate Heat flux at base of fin
*
      HEATFLUX = -DE(1)*(T(2) - T(1))
      PRINT *, "HEATFLUX = ", HEATFLUX,"W/m^2"
*        
*     QUESTION 2: Calculate temperature gradient across first half of fin
*
*      CENTERNODE = (IE+IB)/2
*      PRINT *, T(IE+2-(CENTERNODE+1)/2), ((CENTERNODE+1)/2)-1
*      
*      TEMPGRAD=(T(CENTERNODE)-T((CENTERNODE+1)/2))/
*     C         ((((CENTERNODE+1)/2)-1)*DIEP(CENTERNODE))  ! difference between CENTERNODE-1 and center node divided by
*      PRINT *, "TEMPGRAD FIRST = ", TEMPGRAD," K/m"       ! CV length
*      TEMPGRAD=(T(CENTERNODE)-T(IE+2-(CENTERNODE+1)/2))/
*     C         ((((CENTERNODE+1)/2)-1)*DIEP(CENTERNODE))  ! difference between CENTERNODE-1 and center node divided by
*      PRINT *, "TEMPGRAD SECOND = ", TEMPGRAD," K/m"       ! CV length
*
*--Save result to unformatted output file
*
      CALL SAVE(T,IRSO,IB,IE,ID)
*
      STOP
      END
