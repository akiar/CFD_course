*
*        file solt.f 
************************************************************************
*
      SUBROUTINE SOLT(T,
     C        ATP,ATW,ATE,ATS,ATN,BT,
     C        IB,IE,JB,JE,ID,JD)
*
*  Routine to organize the solution of T. Routine uses BiConjugate
*  Gradient solver DBCG with either Diagonal or Incomplete LU
*  scaling.  ISOLVE is set in this routine since it is not something
*  that is changed often.
*
*  ISOLVE=1 uses Incomplete LU BiConjugate Gradient solver.
*  ISOLVE=2 uses Diagonally-scaled BiConjugate Gradient solver.
*
************************************************************************
*
*  Declaration of current solution and active coefficients 
*
      REAL T(ID,JD)
      REAL ATW(ID,JD),ATE(ID,JD),ATS(ID,JD)
      REAL ATN(ID,JD),ATP(ID,JD),BT(ID,JD)
      INTEGER IB,IE,JB,JE,ID,JD
*  
*  Declaration of storage arrays for solver; DLAP Triad format is used;
*  MR is the maximum number of equations, and MG should be 5*MR
*
      PARAMETER(MR=250000)
      PARAMETER(MG=5*MR)
      DOUBLE PRECISION X(MR),B(MR)
      DOUBLE PRECISION A(MG)
      INTEGER IA(MG),JA(MG)
      INTEGER NEQ,NELT,ISYM,ITOL,ITMAX,ISOLVE
      INTEGER ITER,IERR,IUNIT,LENW,LENIW,IWORK(MG+4*MR)
      DOUBLE PRECISION ERR,TOL,RWORK(MG+8*MR)
*
*  Select solver by setting ISOLVE:
*
      ISOLVE= 1
*
*  Absorb boundary conditions to reduce the size of the global 
*  coefficient matrix
*
      CALL BCELEP(ATP,ATW,ATE,ATS,ATN,BT,IB,IE,JB,JE,ID,JD)
*
*  Calculate order of matrix for solver routines
*
      NEQ = (IE-IB+1)*(JE-JB+1)
*
*  Convert the coeff matrix to DLAP Triad format storing in GLOBM(NELT),
*  IA(NELT),JA(NELT), and store the B vector in RHS(NEQ)
*
      CALL GLOBMI(A,B,IA,JA,
     C            ATP,ATW,ATE,ATS,ATN,BT,
     C            IB,IE,JB,JE,ID,JD,
     C            MG,MR,NELT,NEQ)
*
*  Provide initial guess for DSDBCG solver. Current solution is passed 
*  as an initial guess in X(); solution to system is then also 
*  returned in X().
*
      CALL IGUESS(X, T,IB,IE,JB,JE,ID,JD,MR)
*
*  Call the appropriate solver routine
*
      IF(ISOLVE.EQ.1) THEN
*
*  Solve system using Incomplete LU BiConjugate Gradient Sparse Ax=b 
*  solver DSLUBC. Routine to solve a linear system  Ax = b  using the
*  BiConjugate Gradient  method  with  Incomplete LU decomposition 
*  preconditioning.
*
      ISYM= 0
      ITOL= 1
      TOL=  1.0D-08
      ITMAX=500
      IUNIT=73
      LENW= NELT+NEQ+8*NEQ
      LENIW=NELT+NEQ+4*NEQ+12
*
      CALL DSLUBC(NEQ,B,X,NELT,IA,JA,A,ISYM,ITOL,TOL,
     C            ITMAX,ITER,ERR,IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
      ELSE
*
*  Solve system using Diagonally Scaled BiConjugate Gradient Sparse 
*  Ax=b solver.  Solves a linear system  Ax = b  using the BiConjugate 
*  Gradient method with diagonal scaling.
*
      ISYM= 0
      ITOL= 1
      TOL=  1.0D-08
      ITMAX=500
      IUNIT=73
      LENW= 8*NEQ+2
      LENIW=10+2
*
      CALL DSDBCG(NEQ,B,X,NELT,IA,JA,A,ISYM,ITOL,TOL,
     C            ITMAX,ITER,ERR,IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
      ENDIF
*
*  Insert new solution into T array
*
      CALL ASSGNI(T, X,IB,IE,JB,JE,ID,JD,MR)
*
*  Update the boundary values
*
      CALL BCCLCI(T, ATP,ATW,ATE,ATS,ATN,BT,
     C            IB,IE,JB,JE,ID,JD)
*
      RETURN
      END

