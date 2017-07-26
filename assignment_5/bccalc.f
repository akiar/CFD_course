*
*     This file contains 2 subroutines: BCCLCI and BCCLCU
*
**********************************************************************
      SUBROUTINE BCCLCI(PHI, AP,AW,AE,AS,AN,B,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to update the phi field boundary values.  Used for 
*     independent variables.
*
*     PHI(ID,JD) updated estimate of phi field; output
*
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
*     NOTE:  this routine must be altered to give the correct
*            conditions for the corners of the grid.  Although
*            unimportant for the calculations, the output plots
*            may appear wrong if the corners are improperly set.
*
**********************************************************************
*
      REAL AP(ID,JD),AW(ID,JD),AE(ID,JD),AS(ID,JD)
      REAL AN(ID,JD),B(ID,JD),PHI(ID,JD)
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*     
*  West and East edges of domain
*
      DO 10 J=JB,JE
*
      I=IB-1
      PHI(I,J)=(AE(I,J)*PHI(I+1,J)+B(I,J))/AP(I,J)
*
      I=IE+1
      PHI(I,J)=(AW(I,J)*PHI(I-1,J)+B(I,J))/AP(I,J)
 10   CONTINUE
*
*  South and North edges of domain
*     
      DO 20 I=IB,IE
*
      J=JB-1
      PHI(I,J)=(AN(I,J)*PHI(I,J+1)+B(I,J))/AP(I,J)
*
      J=JE+1
      PHI(I,J)=(AS(I,J)*PHI(I,J-1)+B(I,J))/AP(I,J)
 20   CONTINUE
*
*  Lower and upper corners of domain
*
      PHI(IB-1,JB-1)= PHI(IB,JB-1)
      PHI(IE+1,JB-1)= PHI(IE,JB-1)
      PHI(IB-1,JE+1)= PHI(IB,JE+1)
      PHI(IE+1,JE+1)= PHI(IE,JE+1)
*
      RETURN
      END
*
**********************************************************************
      SUBROUTINE BCCLCU(PHI, AP,AW,AE,AS,AN,B,
     C                  IV,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to update the phi field boundary values.  Used for
*     N simulataneously solved variables.
*
*     PHI(ID,JD) updated estimate of phi field; output
*
*     AP(N,N,ID,JD) active coefficient for P node; input
*     AW(N,N,ID,JD) active coefficient for W node; input
*     AE(N,N,ID,JD) active coefficient for E node; input
*     AS(N,N,ID,JD) active coefficient for S node; input
*     AN(N,N,ID,JD) active coefficient for N node; input
*     B(N,ID,JD) accumulated fixed source term; input
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*     INTEGER N      number of simultaneous variables; input
*
**********************************************************************
*
      REAL AP(N,N,ID,JD),AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),B(N,ID,JD)
      REAL PHI(ID,JD)
      INTEGER IV,IB,IE,JB,JE,I,J,N,ID,JD
*     
*  West and East edges of domain
*
      DO 10 J=JB,JE
*
      I=IB-1
      PHI(I,J)= (AE(IV,IV,I,J)*PHI(I+1,J)
     C         +B(IV,I,J))/AP(IV,IV,I,J)
*
      I=IE+1
      PHI(I,J)= (AW(IV,IV,I,J)*PHI(I-1,J)
     C         +B(IV,I,J))/AP(IV,IV,I,J)
 10   CONTINUE
*
*  South and North edges of domain
*     
      DO 20 I=IB,IE
*
      J=JB-1
      PHI(I,J)= (AN(IV,IV,I,J)*PHI(I,J+1)
     C         +B(IV,I,J))/AP(IV,IV,I,J)
*
      J=JE+1
      PHI(I,J)= (AS(IV,IV,I,J)*PHI(I,J-1)
     C         +B(IV,I,J))/AP(IV,IV,I,J)
 20   CONTINUE
*
*  Lower and upper corners of domain
*
      PHI(IB-1,JB-1)= PHI(IB,JB-1)
      PHI(IE+1,JB-1)= PHI(IE,JB-1)
      PHI(IB-1,JE+1)= PHI(IB,JE+1)
      PHI(IE+1,JE+1)= PHI(IE,JE+1)
*
      RETURN
      END
