*      This file contains two subroutines: INITAL and SAVE
*
*************************************************************************
*
      SUBROUTINE INITAL(T,P,U,V,UHE,VHN,
     C             T0,P0,U0,V0,IRSTRT,IRSI,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to set initial fields by either fixing values at 
*     0 (IRSTRT=0) or by reading in values from a restart file
*     (IRSTRT=1).
*
*     U0   initial guess for uniform velocity distribution; input
*     V0   initial guess for uniform velocity distribution; input
*     P0   initial guess for uniform pressure distribution; input
*     UU0  initial guess for uniform RS distribution; input
*     EPS0 initial guess for uniform EPS distribution; input
*
*     INTEGER IRSTRT restart parameter =0 uniform =1 readin; input
*     INTEGER IRSI   file number for restart data file; input
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*
*************************************************************************
*
      REAL U(ID,JD),V(ID,JD),P(ID,JD),T(ID,JD)
      REAL UHE(ID,JD),VHN(ID,JD)
      REAL U0,V0,P0,T0
      INTEGER IRSTRT,IRSI,IB,IE,JB,JE
      INTEGER I,J,ID,JD
      INTEGER IBM1,IEP1,JBM1,JEP1
*     
      IBM1=IB-1
      IEP1=IE+1
      JBM1=JB-1
      JEP1=JE+1
*
      IF (IRSTRT .EQ. 1) GO TO 100
*
      DO 40 J=JBM1,JEP1
        DO 30 I=IBM1,IEP1
          T(I,J)=  T0
          P(I,J)=  P0
          U(I,J)=  U0
          V(I,J)=  V0
          UHE(I,J)=U0
          VHN(I,J)=V0
 30     CONTINUE
 40   CONTINUE
      GOTO 200
*
 100  CONTINUE
      READ(IRSI) ((T(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      READ(IRSI) ((P(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      READ(IRSI) ((U(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      READ(IRSI) ((V(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      READ(IRSI) ((UHE(I,J),I=IBM1,IE),J=JB,JE)
      READ(IRSI) ((VHN(I,J),I=IB,IE),J=JBM1,JE)
*
 200  CONTINUE
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE SAVE(T,P,U,V,UHE,VHN,XP,YP,
     C                IRSO,IB,IE,JB,JE,ID,JD)
*
*     Routine to save solution arrays in unformatted form to restarting
*     code.
*
*     INTEGER IRSO  unit number of file for saved data; input
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*
*************************************************************************
*
      REAL T(ID,JD),U(ID,JD),V(ID,JD),P(ID,JD)
      REAL UHE(ID,JD),VHN(ID,JD),XP(ID),YP(JD)
      INTEGER IRSO,IB,IE,JB,JE,ID,JD
      INTEGER IBM1,IEP1,JBM1,JEP1,I,J
*
      OPEN(UNIT=31,FILE='mtl.dat')
      OPEN(UNIT=20,FILE='tec.dat')
*     
      IBM1=IB-1
      IEP1=IE+1
      JBM1=JB-1
      JEP1=JE+1
*     
*  Write unformatted output for restart (if IRSTRT=1)
*
      WRITE(IRSO) ((T(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      WRITE(IRSO) ((P(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      WRITE(IRSO) ((U(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      WRITE(IRSO) ((V(I,J),I=IBM1,IEP1),J=JBM1,JEP1)
      WRITE(IRSO) ((UHE(I,J),I=IBM1,IE),J=JB,JE)
      WRITE(IRSO) ((VHN(I,J),I=IB,IE),J=JBM1,JE)
*
*  Write output file for MATlab postprocessing
*
      WRITE(31,3000)
      DO 1 J=JB-1,JE+1
        DO 2 I=IB-1,IE+1
          WRITE(31,3001) 
     C          I,J,XP(I),YP(J),P(I,J),U(I,J),V(I,J),T(I,J)
 2      CONTINUE
 1    CONTINUE
*
 3000 FORMAT(4X,'I    ',1X,'J    ',4X,
     C          'X',13X,'Y',13X,'P',13X,'U',13X,'V',13X,'T')
 3001 FORMAT(I5,1X,I5,6(1X,1PE13.5))
*
*  Write output in TECplot format for postprocessing
*
      WRITE(20,320)
      WRITE(20,322) 0,0,1,(JE-JB+3),(IE-IB+3)
      DO I=IB-1,IE+1
        DO J=JB-1,JE+1
          WRITE(20,325) XP(I),YP(J),P(I,J),U(I,J),V(I,J),T(I,J)
        END DO
      END DO
 320  FORMAT('VARIABLES="XP","YP","P","U","V","T"')
 322  FORMAT('ZONE T="Zone ',I1,I1,I1,'" I=',I6,' J=',I6,' F=POINT')
 325  FORMAT(6(1X,1PE10.3))
*
      RETURN
      END
