*     
*     This file contains 2 subroutines: ASSGNI and ASSGNS
*
********************************************************************
      SUBROUTINE ASSGNI(PHI, SOLN,IB,IE,JB,JE,ID,JD,M)
*
*     Subroutine to assign the nodal values of PHI(ID,JD) from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      REAL PHI(ID,JD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,I,J,E,M,ID,JD
*
      E=0
      DO 1 J=JB,JE
        DO 2 I=IB,IE
          E=E+1
          PHI(I,J)= REAL(SOLN(E))
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE ASSGNS(P,U,V, SOLN,IB,IE,JB,JE,N,ID,JD,M)
*
*     Subroutine to assign the nodal values of U,V,P from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      REAL P(ID,JD),U(ID,JD),V(ID,JD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,N,M,ID,JD
      INTEGER I,J,E
*
      E=1
      DO 1 J=JB,JE
        DO 2 I=IB,IE
          P(I,J)= REAL(SOLN(E  ))
          U(I,J)= REAL(SOLN(E+1))
          V(I,J)= REAL(SOLN(E+2))
          E=E+N
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE IGUESS(SOLN, PHI,IB,IE,JB,JE,ID,JD,M)
*
*     Subroutine to insert the nodal values of PHI(ID,JD) into the
*     solution vector SOLN(M) as an initial guess to solver.
*
********************************************************************
*
      REAL PHI(ID,JD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,I,J,E,M,ID,JD
*
      E=0
      DO 1 J=JB,JE
        DO 2 I=IB,IE
          E=E+1
          SOLN(E)= DBLE(PHI(I,J))
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE SGUESS(SOLN, P,U,V,IB,IE,JB,JE,N,ID,JD,M)
*
*     Subroutine to assign the nodal values of U,V,P from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      REAL P(ID,JD),U(ID,JD),V(ID,JD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,N,M,ID,JD
      INTEGER I,J,E
*
      E=1
      DO 1 J=JB,JE
        DO 2 I=IB,IE
          SOLN(E  )= DBLE(P(I,J))
          SOLN(E+1)= DBLE(U(I,J))
          SOLN(E+2)= DBLE(V(I,J))
          E=E+N
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END