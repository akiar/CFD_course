*
*       file soluvp.f
************************************************************************
*
      SUBROUTINE SOLUVP(U,V,P,
     C        AUP,AUW,AUE,AUS,AUN,BU,WORK3,WORK4,
     C        IB,IE,JB,JE,N,ID,JD)
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
      REAL U(ID,JD),V(ID,JD),P(ID,JD)
      REAL AUW(N,N,ID,JD),AUE(N,N,ID,JD),AUS(N,N,ID,JD)
      REAL AUN(N,N,ID,JD),AUP(N,N,ID,JD),BU(N,ID,JD)
      REAL WORK3(N,N),WORK4(N)
      INTEGER IB,IE,JB,JE,ID,JD,N
      INTEGER IDATO,ITERMO
*  
*  Declaration of storage arrays for solver; DLAP Triad format is used;
*  MR is the maximum number of equations, and MG should be 5*MR
*
      PARAMETER(MR=1000000)
      PARAMETER(MG=5*MR)
      DOUBLE PRECISION X(MR),B(MR)
      DOUBLE PRECISION A(MG)
      INTEGER IA(MG),JA(MG)
      INTEGER NEQ,NELT,ISYM,ITOL,ITMAX
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
      CALL BCELIM(AUP,AUW,AUE,AUS,AUN,BU,
     C            WORK3,WORK4,IB,IE,JB,JE,N,ID,JD)
*
*  Calculate order of matrix for solver routines
*
      NEQ = N*(IE-IB+1)*(JE-JB+1)
*
*  Convert the coeff matrix to DLAP Triad format storing in GLOBM(NELT),
*  IA(NELT),JA(NELT), and store the B vector in RHS(NEQ)
*
      CALL GLOBMS(A,B,IA,JA,
     C            AUP,AUW,AUE,AUS,AUN,BU,
     C            IB,IE,JB,JE,N,ID,JD,
     C            MG,MR,NELT,NEQ)
*
*  Provide initial guess for DSDBCG solver. Current solution is passed 
*  as an initial guess in X(); solution to system is then also 
*  returned in X().
*
      CALL SGUESS(X, P,U,V,IB,IE,JB,JE,N,ID,JD,MR)
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
      TOL=  1.0D-06
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
*--Insert the new solution into the P(),U() and V() arrays
*
      CALL ASSGNS(P,U,V, X,IB,IE,JB,JE,N,ID,JD,M)
*
*--Update the boundary values
*
      CALL BCCLCU(P, AUP,AUW,AUE,AUS,AUN,BU,
     C            1,IB,IE,JB,JE,N,ID,JD)
      CALL BCCLCU(U, AUP,AUW,AUE,AUS,AUN,BU,
     C            2,IB,IE,JB,JE,N,ID,JD)
      CALL BCCLCU(V, AUP,AUW,AUE,AUS,AUN,BU,
     C            3,IB,IE,JB,JE,N,ID,JD)
*
      RETURN
      END

