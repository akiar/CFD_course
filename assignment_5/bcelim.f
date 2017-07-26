*
*     This file contains 2 main subroutines: BCELIM and BCELEP.  The
*     additional routines are used in BCELIM for matrix algebra.
*
**********************************************************************
      SUBROUTINE BCELIM(AP,AW,AE,AS,AN,BP,
     C                  A,B,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to absorb the boundary conditions into the equations
*     at the first interior nodes.  Used for N simultaneous variables.
*
*     AP(N,N,ID,JD) active coefficient for P node; output
*     AW(N,N,ID,JD) active coefficient for W node; output
*     AE(N,N,ID,JD) active coefficient for E node; output
*     AS(N,N,ID,JD) active coefficient for S node; output
*     AN(N,N,ID,JD) active coefficient for N node; output
*     B(N,ID,JD)    accumulated fixed source term; output
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*     INTEGER N      number of simultaneous variables; input
*
*     NOTES: 1.) corner control-volumes are visited twice since
*                they are adjacent to two edges.
*
*            2.) for solutions of simultaneous variables, b.c.'s are
*                absorbed using matrix algebra.
*
*            3.) the indices in the active coefficients are:
*
*                  A(row, column, i-position, j-position)
*
*                  B(row, i-position, j-position)
*
***********************************************************************
*
      REAL AP(N,N,ID,JD),AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),BP(N,ID,JD)
      REAL A(N,N),B(N)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,N,ID,JD
*     
*  West and East edges of domain
*
      DO 10 J=JB,JE
         I=IB
         CALL MMULT(A,  AW,AE,I,J,I-1,J,N,ID,JD)
         CALL MVMULT(B, AW,BP,I,J,I-1,J,N,ID,JD)
         CALL MSUBT(AP, AP,A,I,J,N,ID,JD)
         CALL VPLUS(BP, BP,B,I,J,N,ID,JD)
         CALL ZEROM(AW, I,J,N,ID,JD)
*
         I=IE
         CALL MMULT(A,  AE,AW,I,J,I+1,J,N,ID,JD)
         CALL MVMULT(B, AE,BP,I,J,I+1,J,N,ID,JD)
         CALL MSUBT(AP, AP,A,I,J,N,ID,JD)
         CALL VPLUS(BP, BP,B,I,J,N,ID,JD)
         CALL ZEROM(AE, I,J,N,ID,JD)
 10   CONTINUE
*
*  South and North edges of domain
*
      DO 20 I=IB,IE
         J=JB
         CALL MMULT(A,  AS,AN,I,J,I,J-1,N,ID,JD)
         CALL MVMULT(B, AS,BP,I,J,I,J-1,N,ID,JD)
         CALL MSUBT(AP, AP,A,I,J,N,ID,JD)
         CALL VPLUS(BP, BP,B,I,J,N,ID,JD)
         CALL ZEROM(AS, I,J,N,ID,JD)
*
         J=JE
         CALL MMULT(A,  AN,AS,I,J,I,J+1,N,ID,JD)
         CALL MVMULT(B, AN,BP,I,J,I,J+1,N,ID,JD)
         CALL MSUBT(AP, AP,A,I,J,N,ID,JD)
         CALL VPLUS(BP, BP,B,I,J,N,ID,JD)
         CALL ZEROM(AN, I,J,N,ID,JD)
 20   CONTINUE
*
      RETURN
      END
*
**********************************************************************
      SUBROUTINE BCELEP(AP,AW,AE,AS,AN,B,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to absorb the boundary conditions into the equations
*     at the first interior nodes.  Used for independent variables.
*
*     AP(ID,JD) active coefficient for P node; output
*     AW(ID,JD) active coefficient for W node; output
*     AE(ID,JD) active coefficient for E node; output
*     AS(ID,JD) active coefficient for S node; output
*     AN(ID,JD) active coefficient for N node; output
*     B(ID,JD) accumulated fixed source term; output
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER ID,JD  array dimensions; input
*
*     NOTE:  corner control-volumes are visited twice since
*            they are adjacent to two edges.
*
***********************************************************************
      REAL AP(ID,JD),AW(ID,JD),AE(ID,JD),AS(ID,JD)
      REAL AN(ID,JD),B(ID,JD)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,ID,JD
*
*  West and East edges of domain
*
      DO 10 J=JB,JE
         I=IB
         AP(I,J)=AP(I,J)-AW(I,J)*AE(I-1,J)/AP(I-1,J)
         B(I,J)=B(I,J)+AW(I,J)*B(I-1,J)/AP(I-1,J)
         AW(I,J)=0.0
*
         I=IE
         AP(I,J)=AP(I,J)-AE(I,J)*AW(I+1,J)/AP(I+1,J)
         B(I,J)=B(I,J)+AE(I,J)*B(I+1,J)/AP(I+1,J)
         AE(I,J)=0.0
 10   CONTINUE
*
*  South and North edges of domain
*
      DO 20 I=IB,IE
         J=JB
         AP(I,J)=AP(I,J)-AS(I,J)*AN(I,J-1)/AP(I,J-1)
         B(I,J)=B(I,J)+AS(I,J)*B(I,J-1)/AP(I,J-1)
         AS(I,J)=0.0
*
         J=JE
         AP(I,J)=AP(I,J)-AN(I,J)*AS(I,J+1)/AP(I,J+1)
         B(I,J)=B(I,J)+AN(I,J)*B(I,J+1)/AP(I,J+1)
         AN(I,J)=0.0
 20   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE VPLUS(B3, B1,B2,I,J,N,ID,JD)
*
*     Routine to add two vectors, B1 and B2. The result is stored in
*     vector B3 on output.  B2 has dimension B2(N).
*
************************************************************************
      REAL B3(N,ID,JD),B1(N,ID,JD),B2(N)
      INTEGER I,J,K,N,ID,JD
*
      DO 10 K=1,N
        B3(K,I,J)=B1(K,I,J)+B2(K)
 10   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MSUBT(A3, A1,A2,I,J,N,ID,JD)
*
*     Routine to subtract two matrices, A1 and A2. The result (A1-A2) 
*     is stored in matrix A3 on output.  A2 has dimension A3(N,N).
*
************************************************************************
      REAL A3(N,N,ID,JD),A1(N,N,ID,JD),A2(N,N)
      INTEGER I,J,K,L,N,ID,JD
*
      DO 20 K=1,N
        DO 10 L=1,N
          A3(K,L,I,J)= A1(K,L,I,J)-A2(K,L)
 10     CONTINUE
 20   CONTINUE
* 
      RETURN
      END
*
************************************************************************
      SUBROUTINE MMULT(A3, A1,A2,I1,J1,I2,J2,N,ID,JD)
*
*     Routine to multiply two matrices, A1 and A2. The result (A1*A2)
*     is stored in A3 on output.  A3 has dimension A3(N,N).
*
*************************************************************************
      REAL A3(N,N),A1(N,N,ID,JD),A2(N,N,ID,JD)
      INTEGER I1,J1,I2,J2,K,L,M,N,ID,JD
*
      DO 30 K=1,N
        DO 20 L=1,N
          A3(K,L)= 0.0
          DO 10 M=1,N
            A3(K,L)= A1(K,M,I1,J1)*A2(M,L,I2,J2)+A3(K,L)
 10       CONTINUE
 20     CONTINUE
 30   CONTINUE
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE MVMULT(B2, A1,B1,I1,J1,I2,J2,N,ID,JD)
*
*     Routine to multiply a vector, B1, by a matrix, A1. The result
*     is stored in the vector B2 on output.  B2 has dimension B2(N).
*
***********************************************************************
      REAL B2(N),A1(N,N,ID,JD),B1(N,ID,JD)
      INTEGER I1,J1,I2,J2,K,L,M,N,ID,JD
*
      DO 20 K=1,N
        B2(K)= 0.0
        DO 10 M=1,N
          B2(K)= A1(K,M,I1,J1)*B1(M,I2,J2)+B2(K)
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE ZEROM(A1, I,J,N,ID,JD)
*
*     Routine to zero a matrix A1.
*
***********************************************************************
      REAL A1(N,N,ID,JD)
      INTEGER I,J,K,L,N,ID,JD
*
      DO 20 K=1,N
        DO 10 L=1,N
          A1(K,L,I,J)= 0.0
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END

