*
*         file main.f
************************************************************************
*
*  Program MAIN solves 2-dimensional laminar flows on 2-dimensional 
*  cartesian grids with equal-sized control-volumes.
*
*  The field equations are discretized using the finite-volume 
*  formulation of Patankar, but using collocated variable arrangement.
*  The default advection scheme is UDS.
*
***********************************************************************
*
*============================
*  Declaration of Variables
*============================
*
      PARAMETER(ID=502)
      PARAMETER(JD=502)
      PARAMETER(N=3)
      PARAMETER(GEE=9.80665)
      PARAMETER(THETA=00.0)
      PARAMETER(PI=3.14159)
*  
*--Arrays for geometry description
*
      REAL XNE(ID),YNE(JD)
      REAL XP(ID),YP(JD),XE(ID),YN(JD)
      REAL DIEP(ID),DJNP(JD),DISE(ID),DISN(JD)
      REAL AREP(JD),ARNP(ID),VOLP(ID,JD)
      REAL X_IB,X_IE,Y_JB,Y_JE
      INTEGER IBW,IEW,JBW,JEW
*
*--Arrays for active coefficients and source terms
*
      REAL DE(ID,JD),DN(ID,JD),ME(ID,JD),MN(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD),QT(ID,JD),RT(ID,JD)
      REAL QU(ID,JD),RU(ID,JD),QV(ID,JD),RV(ID,JD),BC(ID,JD)
*
      REAL ACUW(ID,JD),ACUE(ID,JD),ACVS(ID,JD),ACVN(ID,JD)
*
      REAL ATW(ID,JD),ATE(ID,JD),ATS(ID,JD)
      REAL ATN(ID,JD),ATP(ID,JD),BT(ID,JD)
      REAL AUW(N,N,ID,JD),AUE(N,N,ID,JD),AUS(N,N,ID,JD)
      REAL AUN(N,N,ID,JD),AUP(N,N,ID,JD),BU(N,ID,JD)
      REAL DCCE(ID,JD),DCCN(ID,JD)
*
*--Arrays of active variables
*
      REAL U(ID,JD),V(ID,JD),P(ID,JD),T(ID,JD)
      REAL DPDX(ID,JD),DPDY(ID,JD)
      REAL UHE(ID,JD),VHN(ID,JD),TOLD(ID,JD)
      REAL UOLD(ID,JD),VOLD(ID,JD)
      REAL U0,V0,P0,T0
*
*--Arrays of coefficients for coupling mass-momentum equations
*
      REAL DHUE(ID,JD),DHVN(ID,JD)
*
*--Arrays for reporting convergence of equations
*
      REAL NORM(4),RSD(ID,JD),ARSD(4),MAXRSD
      REAL RSDMAX(4),RSDOLD(4),CRIT
*
*--Declaration of all other variables
*
      INTEGER IB,IE1,IE,JB,JE1,JE,I,J
      INTEGER IDATI,IRSI,IDATO,IRSO,ITERMO
      INTEGER IRSTRT,PFTIME,KNTTM,KNTUVP
      INTEGER LVLGEO,LVLCOF,LVLMGD,ADVSCM
      INTEGER KNTOUT,KNTIN,PF,I_MAX,J_MAX
*
      REAL GRDXI,GRDX,GRDYI,GRDY
      REAL IRATX1,IRATX2,IRATY1,IRATY2
      REAL RHO,COND,CP,VISC,BETA,GEX,GEY
      REAL ITIME,TTIME,DTIME
      REAL WORK3(N,N),WORK4(N)
      CHARACTER*8 CHTIME
*
      REAL BOTTOMHT,TOPHT
      REAL SUMTOP,SUMBOTTOM,NUSSELT,RAYLEIGH
*
*============================
*  Initialization and input
*============================
*
      GEX= GEE*SIN(THETA*PI/180.0)
      GEY= GEE*COS(THETA*PI/180.0)
*
*--Read input parameters
*
      CALL FILDEF(IDATI,IRSI,IDATO,IRSO,ITERMO)
      CALL INPUT(IB,IE,JB,JE,
     C           GRDX,GRDY,
     C           RHO,COND,CP,VISC,BETA,
     C           T0,U0,V0,P0,IRSTRT,ITIME,
     C           DTIME,PFTIME,KNTTM,KNTUVP,CRIT,
     C           LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM,IDATI)
      CALL ECHO(IB,IE,JB,JE,
     C          GRDX,GRDY,
     C          RHO,COND,CP,VISC,BETA,
     C          T0,U0,V0,P0,IRSTRT,ITIME,
     C          DTIME,PFTIME,KNTTM,KNTUVP,CRIT,
     C          LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM,IDATO)
*
*--Create grid
*
      CALL CARGRD(XNE,YNE, GRDXI,GRDX,GRDYI,GRDY,
     C            IRATX1,IRATX2,IRATY1,IRATY2,
     C            IBW,IEW,JBW,JEW,X_IB,X_IE,Y_JB,Y_JE,
     C            XP,YP,XE,YN,DIEP,DJNP,DISE,DISN,
     C            AREP,ARNP,VOLP,
     C            IDATO,LVLGEO,IB,IE1,IE,JB,JE1,JE,ID,JD)
*
*--Initialize field variables
*
      CALL INITAL(T,P,U,V,UHE,VHN,
     C            T0,P0,U0,V0,IRSTRT,IRSI,IB,IE,JB,JE,ID,JD)
*
*--Print intialized field variables to output file
*
      WRITE(IDATO,6000)
 6000 FORMAT('       Field Variables   ',/)
      CALL OUT2D(T  ,'T(init) ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      CALL OUT2D(U  ,'U(init) ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      CALL OUT2D(V  ,'V(init) ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      CALL OUT2D(P  ,'P(init) ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
*
*--Compute pressure gradient field from initial conditions
*
      CALL GRADP(DPDX,DPDY, P,DIEP,DJNP,IB,IE,JB,JE,ID,JD)
*
*--Initialize residuals (for RATE calculations)
*
      DO 22 I=1,4
        ARSD(I)=   0.0
	    RSDOLD(I)= 1.0E+10
	    NORM(I)=   1.0E+10
  22  CONTINUE
*
*--Compute coefficients in conservation of mass equation (constant
*  for fixed-grid, incompressible flow)
*
      CALL COEFCN(ACUW,ACUE,ACVS,ACVN,BC,
     C            RHO,AREP,ARNP,IB,IE,JB,JE,ID,JD)
      IF (LVLCOF .GE. 4) THEN
	    CALL OUT2D(ACUW,'ACUW    ',IDATO,IB,IE,1,JB,JE,1,ID,JD)
	    CALL OUT2D(ACUE,'ACUE    ',IDATO,IB,IE,1,JB,JE,1,ID,JD)
	    CALL OUT2D(ACVS,'ACVS    ',IDATO,IB,IE,1,JB,JE,1,ID,JD)
	    CALL OUT2D(ACVN,'ACVN    ',IDATO,IB,IE,1,JB,JE,1,ID,JD)
      ENDIF
*
*=================================
*  Begin outer, time loop
*=================================
*
      PF=0
      TTIME= ITIME
*
      DO 2000 KNTOUT=1,KNTTM
*
      PF=PF+1
      TTIME= TTIME+DTIME
*
*--Transfer variables to old values
*
      DO 10 J=JB-1,JE+1
      DO 10 I=IB-1,IE+1
        TOLD(I,J)=  T(I,J)
  	    UOLD(I,J)=  U(I,J)
	    VOLD(I,J)=  V(I,J)
 10   CONTINUE
*
*----------------------------------
*  Begin inner, linearization loop
*----------------------------------
*     
      DO 1700 KNTIN=1,KNTUVP
*
      WRITE(IDATO,7220) KNTOUT,KNTIN,TTIME
      WRITE(ITERMO,7220) KNTOUT,KNTIN,TTIME
*
*--Compute active coefficients for T
*
      IF(ISOTHM.GT.0) THEN
        CALL NULLM(BT, IB,IE,JB,JE,ID,JD)
        CALL DIFPHI(DE,DN, COND/CP,AREP,ARNP,DIEP,DJNP,
     C              IB,IE,JB,JE,ID,JD)
        CALL MASFLX(ME,MN, UHE,VHN,ACUW,ACUE,ACVS,ACVN,
     C              IB,IE,JB,JE,ID,JD)
        CALL WEIGHT(ALFAE,ALFAN, ME,MN,DE,DN,
     C              IB,IE,JB,JE,ID,JD)
        CALL HOCONV(DCCE,DCCN, ALFAE,ALFAN,ME,MN,DE,DN,
     C          XP,XNE,YP,YNE,DISE,DIEP,DISN,DJNP,RHO,VOLP,
     C          U,UHE,V,VHN,T,
     C          ADVSCM,1.0,IB,IE,JB,JE,ID,JD,XE,YN)
        CALL SRCT(QT,RT, DCCE,DCCN,IB,IE,JB,JE,ID,JD)
        CALL COEFF(ATP,ATW,ATE,ATS,ATN,BT,
     C             ME,MN,DE,DN,QT,RT,TOLD,VOLP,
     C             ALFAE,ALFAN,DIEP,DJNP,
     C             RHO,DTIME,
     C             IB,IE,JB,JE,ID,JD)
*
*--Set boundary conditions in T equation
*
        CALL BNDCT(ATP,ATW,ATE,ATS,ATN,BT,
     C             IB,IE,JB,JE,ID,JD)
*
*--Print coefficients for energy equation (LVLCOF>0 to print)
*
        CALL TCOF(DE,DN,ALFAE,ALFAN,QT,RT,
     C            ATP,ATW,ATE,ATS,ATN,BT,
     C            LVLCOF,IDATO,IB,IE,JB,JE,ID,JD)
      ENDIF
*
*--Compute active coefficients for U & V momentum (nulling required
*  for adcont.f)
*
        CALL NULLMM(BU, IB,IE,JB,JE,N,ID,JD)
        CALL NULLMN(AUP,IB,IE,JB,JE,N,ID,JD)
        CALL NULLMN(AUW,IB,IE,JB,JE,N,ID,JD)
        CALL NULLMN(AUE,IB,IE,JB,JE,N,ID,JD)
        CALL NULLMN(AUS,IB,IE,JB,JE,N,ID,JD)
        CALL NULLMN(AUN,IB,IE,JB,JE,N,ID,JD)
*
        CALL DIFPHI(DE,DN, VISC,AREP,ARNP,DIEP,DJNP,
     C              IB,IE,JB,JE,ID,JD)
        CALL MASFLX(ME,MN, UHE,VHN,ACUW,ACUE,ACVS,ACVN,
     C              IB,IE,JB,JE,ID,JD)
        CALL WEIGHT(ALFAE,ALFAN, ME,MN,DE,DN,
     C              IB,IE,JB,JE,ID,JD)
        CALL HOCONV(DCCE,DCCN, ALFAE,ALFAN,ME,MN,DE,DN,
     C          XP,XNE,YP,YNE,DISE,DIEP,DISN,DJNP,RHO,VOLP,
     C          U,UHE,V,VHN,U,
     C          ADVSCM,1.0,IB,IE,JB,JE,ID,JD,XE,YN)
        CALL SRCU(QU,RU, DCCE,DCCN,T,RHO,VOLP,GEX,BETA,
     C            IB,IE,JB,JE,ID,JD)
        CALL HOCONV(DCCE,DCCN, ALFAE,ALFAN,ME,MN,DE,DN,
     C          XP,XNE,YP,YNE,DISE,DIEP,DISN,DJNP,RHO,VOLP,
     C          U,UHE,V,VHN,V,
     C          ADVSCM,1.0,IB,IE,JB,JE,ID,JD,XE,YN)
        CALL SRCV(QV,RV, DCCE,DCCN,T,RHO,VOLP,GEY,BETA,
     C            IB,IE,JB,JE,ID,JD)
        CALL COEFFM(AUP,AUW,AUE,AUS,AUN,BU,
     C              ME,MN,DE,DN,QU,RU,UOLD,VOLP,
     C              ALFAE,ALFAN,DIEP,DJNP,
     C              RHO,DTIME,
     C              2,IB,IE,JB,JE,N,ID,JD)
        CALL COEFFM(AUP,AUW,AUE,AUS,AUN,BU,
     C              ME,MN,DE,DN,QV,RV,VOLD,VOLP,
     C              ALFAE,ALFAN,DIEP,DJNP,
     C              RHO,DTIME,
     C              3,IB,IE,JB,JE,N,ID,JD)
*
*--Add the pressure terms to the momentum equations
*
        CALL SRCUVP(AUP,AUW,AUE,AUS,AUN,BU,
     C             VOLP,DIEP,DISE,DJNP,DISN,
     C             IB,IE,JB,JE,N,ID,JD)
*
*--Set the boundary conditions in the momentum equations
*
        CALL BNDCU(AUP,AUW,AUE,AUS,AUN,BU,
     C             IB,IE,JB,JE,N,ID,JD)
        CALL BNDCV(AUP,AUW,AUE,AUS,AUN,BU,
     C             IB,IE,JB,JE,N,ID,JD)
*
*--Compute terms required to couple mass-momentum
*
        CALL DHAT(DHUE,DHVN, AUP,DIEP,DISE,DJNP,DISN,
     C            VOLP,RHO,IB,IE,JB,JE,N,ID,JD)
*
*--Insert the mass coefficients into the coefficient blocks
*
        CALL ADCONT(AUP,AUW,AUE,AUS,AUN,BU,
     C              ACUW,ACUE,ACVS,ACVN,DHUE,DHVN,
     C              DPDX,DPDY,RHO,
     C              DIEP,DISE,DJNP,DISN,
     C              IB,IE,JB,JE,N,ID,JD)
*
*--Set the Boundary conditions in the P (mass) equation
*
        CALL BNDCP(AUP,AUW,AUE,AUS,AUN,BU,
     C             P,DIEP,DJNP,IB,IE,JB,JE,N,ID,JD)
*
*--Print coefficients for mass-momentum equations (LVLCOF>0 to print)
*
        CALL UVPCOF(DE,DN,ME,MN,ALFAE,ALFAN,QU,RU,QV,RV,
     C              DHUE,DHVN,
     C              DPDX,DPDY,AUP,AUW,AUE,AUS,AUN,BU,
     C              LVLCOF,IDATO,IB,IE,JB,JE,N,ID,JD)
*
*--Compute and report normalized residuals for energy, mass and momentum
*
        CALL UVPTRS(T,P,U,V,
     C              NORM,RSD,ARSD,RSDMAX,RSDOLD,I_MAX,J_MAX,
     C              AUP,AUW,AUE,AUS,AUN,BU,
     C              ATP,ATW,ATE,ATS,ATN,BT,
     C              IDATO,ITERMO,ISOTHM,IB,IE,JB,JE,N,ID,JD)
*
*--Check residuals of all field variables against convergence criterion
*
        IF(KNTOUT.GT.1.OR.KNTIN.GT.1) THEN
          MAXRSD= ABS(ARSD(1)/NORM(1))
  	      DO 1710 I=2,4
	        IF(ABS(RSDMAX(I)/NORM(I)).GT.MAXRSD) 
     C      MAXRSD=ABS(RSDMAX(I)/NORM(I))
 1710     CONTINUE
*
	      IF(MAXRSD.LT.CRIT) THEN
            IF(KNTUVP.LT.10) THEN
              WRITE(ITERMO,7216) KNTOUT
              GOTO 1840
            ELSE
              WRITE(ITERMO,7216) KNTOUT
              GOTO 1770
            ENDIF
	      ELSE IF(MAXRSD.GT.1.E+04) THEN
            WRITE(ITERMO,7215) KNTOUT
	        GOTO 9999
	      ENDIF
        ENDIF
*
*--Store current residual levels
*
	    IF(ISOTHM.GT.0) RSDOLD(1)=ARSD(1)/NORM(1)
        DO 1720 I=2,4
          RSDOLD(I)=ARSD(I)/NORM(I)
 1720   CONTINUE
*
*--Solve for T using a direct solver
*
        IF(ISOTHM.GT.0) THEN
           CALL SOLT(T,
     C             ATP,ATW,ATE,ATS,ATN,BT,
     C             IB,IE,JB,JE,ID,JD)
	    ENDIF
*
*--Solve P,U,V as a coupled set using a direct solver
*
        CALL SOLUVP(U,V,P,
     C        AUP,AUW,AUE,AUS,AUN,BU,WORK3,WORK4,
     C        IB,IE,JB,JE,N,ID,JD)
*
*--Update pressure gradient field
*
        CALL GRADP(DPDX,DPDY, P,DIEP,DJNP,IB,IE,JB,JE,ID,JD)
*
*--Update the face velocities for the mass fluxes using the 
*  new solution
*
        CALL UHAT(UHE, U,P,
     C            DHUE,DPDX,RHO,
     C            DIEP,DISE,IB,IE,JB,JE,ID,JD)
        CALL VHAT(VHN, V,P,
     C            DHVN,DPDY,RHO,
     C            DJNP,DISN,IB,IE,JB,JE,ID,JD)
*
*----------------------------------
*  End inner, linearization loop
*----------------------------------
*
 1700 CONTINUE
*
*     COMPUTE HEAT DIFFERENCES
      NUSSELT = 0.0
      RAYLEIGH = 0.0
*
*     QUESTION 3
*
      DO 8000 J=JB,JE
       BOTTOMHT = ABS(DE(IB-1,J)*(T(IB,J)-T(IB-1,J)))
       TOPHT = DE(IE,J)*(T(IE,J)-T(IE+1,J))
       SUMBOTTOM = SUMBOTTOM+ABS(BOTTOMHT)
       SUMTOP = SUMTOP +ABS(TOPHT)
 8000 CONTINUE
*      NUSSELT = ABS(SUMBOTTOM*GRDY/GRDX/1/(T(IB,JB)-T(IB,JB-1))/COND) ! Q2
      NUSSELT = ABS(SUMBOTTOM/GRDY/1/(T(IB-1,JB)-T(IB,JB-1))*GRDX/COND)
      RAYLEIGH = ABS(GEY*BETA*(GRDX*GRDX*GRDX)*
     C            (T(IB-1,JB)-T(IB,JB-1))*CP/COND/VISC*RHO*RHO)
      PRINT *, "NUSSELT: ",NUSSELT
      PRINT *, "RAYLEIGH: ",RAYLEIGH
      PRINT *, "LEFT HEAT: ",SUMBOTTOM
      PRINT *, "RIGHT HEAT",SUMTOP
      IF (ABS(SUMBOTTOM)==ABS(SUMTOP)) THEN
       EXIT
      ENDIF
*
*  Print solution for present time step
*
 1770 IF((PF/PFTIME).EQ.1) THEN
      WRITE(IDATO,6010)
 6010 FORMAT('       Field Variables   ',/)
       CALL OUT2D(T  ,'T       ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
       CALL OUT2D(U  ,'U       ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
       CALL OUT2D(V  ,'V       ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
       CALL OUT2D(P  ,'P       ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      PF=0
      ENDIF
*
      IF(KNTOUT.GE.KNTTM.AND.KNTTM.GT.1) THEN
       WRITE(IDATO,7213) KNTIN-1,KNTOUT
       WRITE(ITERMO,7213) KNTIN-1,KNTOUT
       GOTO 1840
      ENDIF
*
*=================================
*  End outer, time loop
*=================================
*
 2000 CONTINUE
*
*--Write all field variables to an unformatted output file, rso.bin,
*  to to tec.dat, which can be read directly by TECplot, and to 
*  mtl.dat, which can be used to postprocess using MATlab
*
 1840 CALL SAVE(T,P,U,V,UHE,VHN,XP,YP,IRSO,IB,IE,JB,JE,ID,JD)
*
*-----------------------------------------------------------------------
*  Format statements
*
 7213 FORMAT(/,'ITERATION LIMIT REACHED (',I3,
     C         ') WITHOUT CONVERGENCE AT TIME STEP ',I4,/)
 7214 FORMAT(//,'Time Step: ',I5,/)
 7215 FORMAT(/,'JOB STOPPED: MAXRSD EXCEEDS 1.0E+04 AT TIMESTEP ',I4,/)
 7216 FORMAT(/,'JOB CONVERGED AT TIMESTEP ',I4,/)
 7217 FORMAT(1PE13.5,2X,1PE13.5)
 7220 FORMAT(
     C '=======================================================',/
     C,'TIME STEP= 'I4'('I2')     SIMULATION TIME = '1PE13.5    ,/
     C,'-------------------------------------------------------',/
     C,'|Equation  | Rate  | Ave RSD | Max RSD | Max Location |',/
     C,'+----------+-------+---------+---------+--------------+')
 7224 FORMAT(1PE11.4,3(1X,1PE11.4))
 7225 FORMAT(1PE11.4,4(1X,1PE11.4))
 7226 FORMAT(1PE11.4,5(1X,1PE11.4))
*
*--End of main
*
 9999 STOP
      END

