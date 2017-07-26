*
*     This file contains 2 subroutines;  RESID and RESIDM
*
***********************************************************************
*
      SUBROUTINE RESID(RSD,AVRSD,RSDMAX,I_MAX,J_MAX, 
     C                 PHI,AP,AW,AE,AS,AN,B,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the residual at each interior c.v and
*     the average of the absolute residuals over all interior c.v.  
*     Used for independently solved variables.
*
*     RSD(ID,JD) residual array for each interior c.v; output
*     AVRSD average residual for all interior c.v.; output
*
*     PHI(ID,JD) updated estimate of phi field; input
*     AP(ID,JD) active coefficient for P node; input
*     AW(ID,JD) active coefficient for W node; input
*     AE(ID,JD) active coefficient for E node; input
*     AS(ID,JD) active coefficient for S node; input
*     AN(ID,JD) active coefficient for N node; input
*     B(ID,JD) accumulated fixed source term; input
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*
***********************************************************************
*
      REAL RSD(ID,JD),PHI(ID,JD),AP(ID,JD),AW(ID,JD)
      REAL AE(ID,JD),AS(ID,JD),AN(ID,JD),B(ID,JD)
      REAL AVRSD,RSDMAX
      INTEGER IB,IE,JB,JE,I,J,NPTS,I_MAX,J_MAX,ID,JD
*
      AVRSD=0.0
      RSDMAX= 0.0
      NPTS=(IE-IB+1)*(JE-JB+1)
*
      DO 20 J=JB,JE
        DO 10 I=IB,IE
          RSD(I,J)=AW(I,J)*PHI(I-1,J)+AE(I,J)*PHI(I+1,J)
     C            +AS(I,J)*PHI(I,J-1)+AN(I,J)*PHI(I,J+1)
     C            +B(I,J)-AP(I,J)*PHI(I,J)
          IF(ABS(RSD(I,J)).GT.RSDMAX) THEN
            RSDMAX=ABS(RSD(I,J))
            I_MAX= I
            J_MAX= J
          ENDIF
          AVRSD=AVRSD+ABS(RSD(I,J))
 10     CONTINUE
 20   CONTINUE
      AVRSD=AVRSD/NPTS
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE RESIDM(RSD,AVRSD,RSDMAX,I_MAX,J_MAX,
     C                  P,U,V,AP,AW,AE,AS,AN,B,
     C                  IV,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to calculate the residual at each interior c.v and
*     the average of the absolute residuals over all interior c.v.
*     Used for simulataneous solution variables.
*
*     RSD(ID,JD) residual array for each interior c.v; output
*     AVRSD average residual for all interior c.v.; output
*
*     AP(N,N,ID,JD) active coefficient for P node; input
*     AW(N,N,ID,JD) active coefficient for W node; input
*     AE(N,N,ID,JD) active coefficient for E node; input
*     AS(N,N,ID,JD) active coefficient for S node; input
*     AN(N,N,ID,JD) active coefficient for N node; input
*     B(N,ID,JD)    accumulated fixed source term; input
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER N      number of simultaneously solved variables; input
*     INTEGER ID,JD  array dimensions; input
*
***********************************************************************
*
      REAL RSD(ID,JD),AP(N,N,ID,JD),AW(N,N,ID,JD)
      REAL AE(N,N,ID,JD),AS(N,N,ID,JD),AN(N,N,ID,JD)
      REAL B(N,ID,JD),AVRSD,RSDMAX
      REAL P(ID,JD),U(ID,JD),V(ID,JD)
      INTEGER IV,IB,IE,JB,JE,N,ID,JD
      INTEGER I,J,NPTS,I_MAX,J_MAX
*
      AVRSD=0.0
      RSDMAX= 0.0
      NPTS=(IE-IB+1)*(JE-JB+1)
*
      DO 20 J=JB,JE
        DO 10 I=IB,IE
          RSD(I,J)=AW(IV,1,I,J)*P(I-1,J)+AE(IV,1,I,J)*P(I+1,J)
     C            +AS(IV,1,I,J)*P(I,J-1)+AN(IV,1,I,J)*P(I,J+1)
     C            -AP(IV,1,I,J)*P(I,J)
     C            +AW(IV,2,I,J)*U(I-1,J)+AE(IV,2,I,J)*U(I+1,J)
     C            +AS(IV,2,I,J)*U(I,J-1)+AN(IV,2,I,J)*U(I,J+1)
     C            -AP(IV,2,I,J)*U(I,J)
     C            +AW(IV,3,I,J)*V(I-1,J)+AE(IV,3,I,J)*V(I+1,J)
     C            +AS(IV,3,I,J)*V(I,J-1)+AN(IV,3,I,J)*V(I,J+1)
     C            -AP(IV,3,I,J)*V(I,J)
     C            +B(IV,I,J)
          IF(ABS(RSD(I,J)).GT.RSDMAX) THEN
            RSDMAX= ABS(RSD(I,J))
            I_MAX= I
            J_MAX= J
          ENDIF
          AVRSD=AVRSD+ABS(RSD(I,J))
 10     CONTINUE
 20   CONTINUE
      AVRSD=AVRSD/NPTS
*
      RETURN
      END

