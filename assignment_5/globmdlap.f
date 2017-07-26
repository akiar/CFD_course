*
*     file globmdlap.f
*********************************************************************
*
      SUBROUTINE GLOBMI(GLOBM,RHS,IA,JA,
     C                  ATP,ATW,ATE,ATS,ATN,BT,
     C                  IB,IE,JB,JE,ID,JD,
     C                  MG,MR,NELT,N)
*
*     Subroutine to construct the coefficient array in the structure 
*     required by DLAP.
*
*     GLOBM(NELT)    global array in IA(),JA() format; output
*     RHS(N)         right-hand-side vector; output
*     NELT           number of non-zero elements in GLOBM(); output
*     N              order of global coefficient matrix; input
*     
*     AT#(ID,JD)     active coefficient block arrays
*     IB,IE          starting and ending indices in 'x' direction
*     JB,JE          starting and ending indices in 'y' direction
*     ID,JD          array dimensions of block coefficients
*
*     MG,MR          maximum array lengths of GLOBM() and RHS()     
*
*     NOTE:  The structure of the matrix is set up such that the
*            diagonal of each row is the first element in the packed
*            global matrix.
*
*       Here is an example of the SLAP Triad storage format for a
*       5x5 Matrix.  The entries may appear in any order.
*
*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
*                              1  2  3  4  5  6  7  8  9 10 11
*       |11 12  0  0 15| GLOBM: 11 12 15 22 21 33 35 44 55 51 53
*       |21 22  0  0  0|    IA:  1  1  1  2  2  3  3  4  5  5  5
*       | 0  0 33  0 35|    JA:  1  2  3  2  1  3  5  4  5  1  3
*       | 0  0  0 44  0|
*       |51  0 53  0 55|
*
*                            Here; N=5, NELT=11
*
*********************************************************************
*
      DOUBLE PRECISION GLOBM(MG),RHS(MR)
      REAL ATP(ID,JD),ATW(ID,JD),ATE(ID,JD)
      REAL ATS(ID,JD),ATN(ID,JD),BT(ID,JD)
      INTEGER IA(MG),JA(MG)
      INTEGER IB,IE,JB,JE,E,ID,JD
      INTEGER MG,MR,NELT,N
      INTEGER W,I,J
*
*  Compute the half-bandwidth of the block matrix
*
      W=IE-IB+1
*
*-------------------------------------------------------
*  [1] Form the IA(), JA() arrays
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      ROW=1
      E=1
*
*  Bottom row of control-volumes
*
      IA(E)=ROW
      JA(E)=ROW
      E=E+1
      IF(IE.GT.IB) THEN
        IA(E)=ROW
	    JA(E)=ROW+1
	    E=E+1
      ENDIF
      IF(JE.GT.JB) THEN
        IA(E)=ROW
	    JA(E)=ROW+W
	    E=E+1
      ENDIF
*
      DO 5 I=IB+1,IE-1
        ROW=ROW+1
	    IA(E  )=ROW
	    JA(E  )=ROW
        IA(E+1)=ROW
	    JA(E+1)=ROW-1
        IA(E+2)=ROW
	    JA(E+2)=ROW+1
	    E=E+3
        IF(JE.GT.JB) THEN
          IA(E)=ROW
	      JA(E)=ROW+W
	      E=E+1
        ENDIF
  5   CONTINUE
*
      IF(IE.GT.IB) THEN
        I=IE
        ROW=ROW+1
	    IA(E  )=ROW
	    JA(E  )=ROW
        IA(E+1)=ROW
	    JA(E+1)=ROW-1
	    E=E+2
        IF(JE.GT.JB) THEN
          IA(E)=ROW
	      JA(E)=ROW+W
	      E=E+1
        ENDIF
      ENDIF
*
*  Interior control-volumes
*
      DO 10 J=JB+1,JE-1
        I=IB
	    ROW=ROW+1
	    IA(E  )=ROW
	    JA(E  )=ROW
        IA(E+1)=ROW
	    JA(E+1)=ROW-W
	    E=E+2
        IF(IE.GT.IB) THEN
          IA(E)=ROW
	      JA(E)=ROW+1
	      E=E+1
        ENDIF
        IA(E)=ROW
	    JA(E)=ROW+W
	    E=E+1
*
        DO 15 I=IB+1,IE-1
	      ROW=ROW+1
	      IA(E  )=ROW
	      JA(E  )=ROW
          IA(E+1)=ROW
	      JA(E+1)=ROW-W
          IA(E+2)=ROW
	      JA(E+2)=ROW-1
          IA(E+3)=ROW
	      JA(E+3)=ROW+1
          IA(E+4)=ROW
	      JA(E+4)=ROW+W
	      E=E+5
 15     CONTINUE
*
        IF(IE.GT.IB) THEN
          I=IE
	      ROW=ROW+1
	      IA(E  )=ROW
	      JA(E  )=ROW
          IA(E+1)=ROW
	      JA(E+1)=ROW-W
          IA(E+2)=ROW
	      JA(E+2)=ROW-1
          IA(E+3)=ROW
	      JA(E+3)=ROW+W
	      E=E+4
        ENDIF
 10   CONTINUE
*
*  Top row of control-volumes
*
      IF(JE.GT.JB) THEN
        J=JE
        I=IB
	    ROW=ROW+1
	    IA(E  )=ROW
	    JA(E  )=ROW
        IA(E+1)=ROW
	    JA(E+1)=ROW-W
	    E=E+2
        IF(IE.GT.IB) THEN
          IA(E)=ROW
	      JA(E)=ROW+1
	      E=E+1
        ENDIF
*
        DO 20 I=IB+1,IE-1
	      ROW=ROW+1
	      IA(E  )=ROW
	      JA(E  )=ROW
          IA(E+1)=ROW
	      JA(E+1)=ROW-W
          IA(E+2)=ROW
	      JA(E+2)=ROW-1
          IA(E+3)=ROW
	      JA(E+3)=ROW+1
	      E=E+4
 20     CONTINUE
*
        IF(IE.GT.IB) THEN
          I=IE
	      ROW=ROW+1
	      IA(E  )=ROW
	      JA(E  )=ROW
          IA(E+1)=ROW
	      JA(E+1)=ROW-W
          IA(E+2)=ROW
	      JA(E+2)=ROW-1
	      E=E+3
        ENDIF
      ENDIF
*
*  Store the number of non-zero entries in GLOBM()
*
      NELT=E-1
*
*-------------------------------------------------------
*  [3] Fill in the non-zero elements in matrix GLOBM()
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      ROW=1
      E=1
*
*  Bottom row of control-volumes
*
      GLOBM(E)= DBLE( ATP(I,J))
      E=E+1
      IF(IE.GT.IB) THEN
        GLOBM(E)= DBLE(-ATE(I,J))
	    E=E+1
      ENDIF
      IF(JE.GT.JB) THEN
        GLOBM(E)= DBLE(-ATN(I,J))
	    E=E+1
      ENDIF
      RHS(ROW)= DBLE(BT(I,J))
*
      DO 50 I=IB+1,IE-1
        ROW=ROW+1
        GLOBM(E  )= DBLE( ATP(I,J))
        GLOBM(E+1)= DBLE(-ATW(I,J))
        GLOBM(E+2)= DBLE(-ATE(I,J))
	    E=E+3
        IF(JE.GT.JB) THEN
          GLOBM(E)= DBLE(-ATN(I,J))
	      E=E+1
        ENDIF
        RHS(ROW)= DBLE(BT(I,J))
 50   CONTINUE
*
      IF(IE.GT.IB) THEN
        I=IE
        ROW=ROW+1
        GLOBM(E  )= DBLE( ATP(I,J))
        GLOBM(E+1)= DBLE(-ATW(I,J))
	    E=E+2
        IF(JE.GT.JB) THEN
          GLOBM(E)= DBLE(-ATN(I,J))
	      E=E+1
        ENDIF
        RHS(ROW)= DBLE(BT(I,J))
      ENDIF
*
*  [2] Interior control-volumes
*
      DO 55 J=JB+1,JE-1
        I=IB
	    ROW=ROW+1
        GLOBM(E  )= DBLE( ATP(I,J))
        GLOBM(E+1)= DBLE(-ATS(I,J))
	    E=E+2
        IF(IE.GT.IB) THEN
          GLOBM(E)= DBLE(-ATE(I,J))
	      E=E+1
        ENDIF
        GLOBM(E)= DBLE(-ATN(I,J))
	    E=E+1
        RHS(ROW)= DBLE(BT(I,J))
*
        DO 60 I=IB+1,IE-1
	      ROW=ROW+1
          GLOBM(E  )= DBLE( ATP(I,J))
          GLOBM(E+1)= DBLE(-ATS(I,J))
          GLOBM(E+2)= DBLE(-ATW(I,J))
          GLOBM(E+3)= DBLE(-ATE(I,J))
          GLOBM(E+4)= DBLE(-ATN(I,J))
	      E=E+5
          RHS(ROW)= DBLE(BT(I,J))
 60     CONTINUE
*
        IF(IE.GT.IB) THEN
          I=IE
	      ROW=ROW+1
          GLOBM(E  )= DBLE( ATP(I,J))
          GLOBM(E+1)= DBLE(-ATS(I,J))
          GLOBM(E+2)= DBLE(-ATW(I,J))
          GLOBM(E+3)= DBLE(-ATN(I,J))
	      E=E+4
          RHS(ROW)= DBLE(BT(I,J))
        ENDIF
 55   CONTINUE
*
*  Top row of control-volumes
*
      IF(JE.GT.JB) THEN
        J=JE
        I=IB
	    ROW=ROW+1
        GLOBM(E  )= DBLE( ATP(I,J))
        GLOBM(E+1)= DBLE(-ATS(I,J))
	    E=E+2
        IF(IE.GT.IB) THEN
          GLOBM(E)= DBLE(-ATE(I,J))
	      E=E+1
        ENDIF
        RHS(ROW)= DBLE(BT(I,J))
*
        DO 65 I=IB+1,IE-1
	      ROW=ROW+1
          GLOBM(E  )= DBLE( ATP(I,J))
          GLOBM(E+1)= DBLE(-ATS(I,J))
          GLOBM(E+2)= DBLE(-ATW(I,J))
          GLOBM(E+3)= DBLE(-ATE(I,J))
	      E=E+4
          RHS(ROW)= DBLE(BT(I,J))
 65     CONTINUE
*
        IF(IE.GT.IB) THEN
          I=IE
	      ROW=ROW+1
          GLOBM(E  )= DBLE( ATP(I,J))
          GLOBM(E+1)= DBLE(-ATS(I,J))
          GLOBM(E+2)= DBLE(-ATW(I,J))
	      E=E+3
          RHS(ROW)= DBLE(BT(I,J))
        ENDIF
      ENDIF
*
      RETURN
      END
*
*********************************************************************
*
      SUBROUTINE GLOBMS(GLOBM,RHS,IA,JA,
     C                  AUP,AUW,AUE,AUS,AUN,BU,
     C                  IB,IE,JB,JE,N,ID,JD,
     C                  MG,MR,NELT,NEQ)
*
*     Subroutine to construct the coefficient array in the structure 
*     required by DLAP.
*
*     GLOBM(NELT)    global array in IA(),JA() format; output
*     RHS(NEQ)       right-hand-side vector; output
*     NELT           number of non-zero elements in GLOBM(); output
*     NEQ            order of global coefficient matrix; input
*     
*     AU#(N,N,ID,JD) NxN active coefficient block arrays
*     IB,IE          starting and ending indices in 'x' direction
*     JB,JE          starting and ending indices in 'y' direction
*     ID,JD          array dimensions of block coefficients
*
*     MG,MR          maximum array lengths of GLOBM() and RHS()     
*
*     NOTE:  The structure of the matrix is set up such that the
*            diagonal of each row is the first element in the packed
*            global matrix.
*
*       Here is an example of the SLAP Triad storage format for a
*       5x5 Matrix.  The entries may appear in any order.
*
*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
*                              1  2  3  4  5  6  7  8  9 10 11
*       |11 12  0  0 15| GLOBM: 11 12 15 22 21 33 35 44 55 51 53
*       |21 22  0  0  0|    IA:  1  1  1  2  2  3  3  4  5  5  5
*       | 0  0 33  0 35|    JA:  1  2  3  2  1  3  5  4  5  1  3
*       | 0  0  0 44  0|
*       |51  0 53  0 55|
*
*                            Here; NEQ=5, NELT=11
*
*********************************************************************
*
      DOUBLE PRECISION GLOBM(MG),RHS(MR)
      REAL AUP(N,N,ID,JD),AUW(N,N,ID,JD),AUE(N,N,ID,JD)
      REAL AUS(N,N,ID,JD),AUN(N,N,ID,JD),BU(N,ID,JD)
      INTEGER IA(MG),JA(MG)
      INTEGER IB,IE,JB,JE,E,N,ID,JD
      INTEGER MG,MR,NELT,NEQ
      INTEGER W,I,J
*
*  Compute the half-bandwidth
*
      W=N*(IE-IB+1)
*
*  Set the values of the non-zero matrix elements
*
*  [1] Bottom row of control-volumes
*
      J=JB
      I=IB
      E=1
      F=1
      DO 23 K=1,N
       DO 24 L=1,N
         IA(F)= E+K-1
         JA(F)= E+L-1
         GLOBM(F)= DBLE( AUP(K,L,I,J))
         F=F+1
         IF(IE.GT.IB) THEN
          IA(F)= E+K-1
          JA(F)= E+N+L-1
          GLOBM(F)= DBLE(-AUE(K,L,I,J))
          F=F+1
         ENDIF
         IF(JE.GT.JB) THEN
          IA(F)= E+K-1
          JA(F)= E+W+L-1
          GLOBM(F)= DBLE(-AUN(K,L,I,J))
          F=F+1
         ENDIF
  24   CONTINUE
       RHS(E+K-1)= DBLE( BU(K,I,J))
  23  CONTINUE
*
      DO 25 I=IB+1,IE-1
        E=E+N
        DO 26 K=1,N
         DO 27 L=1,N
           IA(F)= E+K-1
           JA(F)= E-N+L-1
           GLOBM(F)= DBLE(-AUW(K,L,I,J))
           F=F+1
           IA(F)= E+K-1
           JA(F)= E+L-1
           GLOBM(F)= DBLE( AUP(K,L,I,J))
           F=F+1
           IA(F)= E+K-1
           JA(F)= E+N+L-1
           GLOBM(F)= DBLE(-AUE(K,L,I,J))
           F=F+1
           IF(JE.GT.JB) THEN
            IA(F)= E+K-1
            JA(F)= E+W+L-1
            GLOBM(F)= DBLE(-AUN(K,L,I,J))
            F=F+1
           ENDIF
  27     CONTINUE
         RHS(E+K-1)= DBLE( BU(K,I,J))
  26    CONTINUE
  25  CONTINUE
*
      IF(IE.GT.IB) THEN
        I=IE
        E=E+N
        DO 28 K=1,N
         DO 29 L=1,N
           IA(F)= E+K-1
           JA(F)= E-N+L-1
           GLOBM(F)= DBLE(-AUW(K,L,I,J))
           F=F+1
           IA(F)= E+K-1
           JA(F)= E+L-1
           GLOBM(F)= DBLE( AUP(K,L,I,J))
           F=F+1
           IF(JE.GT.JB) THEN
            IA(F)= E+K-1
            JA(F)= E+W+L-1
            GLOBM(F)= DBLE(-AUN(K,L,I,J))
            F=F+1
           ENDIF
  29     CONTINUE
         RHS(E+K-1)= DBLE( BU(K,I,J))
  28    CONTINUE
      ENDIF
*
*  [2] Interior control-volumes
*
      DO 30 J=JB+1,JE-1
        I=IB
        E=E+N
        DO 31 K=1,N
         DO 32 L=1,N
           IA(F)= E+K-1
           JA(F)= E-W+L-1
           GLOBM(F)= DBLE(-AUS(K,L,I,J))
           F=F+1
           IA(F)= E+K-1
           JA(F)= E+L-1
           GLOBM(F)= DBLE( AUP(K,L,I,J))
           F=F+1
           IF(IE.GT.IB) THEN
            IA(F)= E+K-1
            JA(F)= E+N+L-1
            GLOBM(F)= DBLE(-AUE(K,L,I,J))
            F=F+1
            IA(F)= E+K-1
            JA(F)= E+W+L-1
            GLOBM(F)= DBLE(-AUN(K,L,I,J))
            F=F+1
           ENDIF
  32     CONTINUE
         RHS(E+K-1)= DBLE( BU(K,I,J))
  31    CONTINUE
        DO 33 I=IB+1,IE-1
          E=E+N
          DO 34 K=1,N
           DO 35 L=1,N
             IA(F)= E+K-1
             JA(F)= E-W+L-1
             GLOBM(F)= DBLE(-AUS(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E-N+L-1
             GLOBM(F)= DBLE(-AUW(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+L-1
             GLOBM(F)= DBLE( AUP(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+N+L-1
             GLOBM(F)= DBLE(-AUE(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+W+L-1
             GLOBM(F)= DBLE(-AUN(K,L,I,J))
             F=F+1
  35       CONTINUE
           RHS(E+K-1)= DBLE( BU(K,I,J))
  34      CONTINUE
  33    CONTINUE
        IF(IE.GT.IB) THEN
          I=IE
          E=E+N
          DO 36 K=1,N
           DO 37 L=1,N
             IA(F)= E+K-1
             JA(F)= E-W+L-1
             GLOBM(F)= DBLE(-AUS(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E-N+L-1
             GLOBM(F)= DBLE(-AUW(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+L-1
             GLOBM(F)= DBLE( AUP(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+W+L-1
             GLOBM(F)= DBLE(-AUN(K,L,I,J))
             F=F+1
  37       CONTINUE
           RHS(E+K-1)= DBLE( BU(K,I,J))
  36      CONTINUE
        ENDIF
  30  CONTINUE
*
*  [3] Top row of control-volumes
*
      IF(JE.GT.JB) THEN
        J=JE
        I=IB
        E=E+N
        DO 38 K=1,N
         DO 39 L=1,N
           IA(F)= E+K-1
           JA(F)= E-W+L-1
           GLOBM(F)= DBLE(-AUS(K,L,I,J))
           F=F+1
           IA(F)= E+K-1
           JA(F)= E+L-1
           GLOBM(F)= DBLE( AUP(K,L,I,J))
           F=F+1
           IF(IE.GT.IB) THEN
            IA(F)= E+K-1
            JA(F)= E+N+L-1
            GLOBM(F)= DBLE(-AUE(K,L,I,J))
            F=F+1
           ENDIF
  39     CONTINUE
         RHS(E+K-1)= DBLE( BU(K,I,J))
  38    CONTINUE
        DO 40 I=IB+1,IE-1
          E=E+N
          DO 41 K=1,N
           DO 42 L=1,N
             IA(F)= E+K-1
             JA(F)= E-W+L-1
             GLOBM(F)= DBLE(-AUS(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E-N+L-1
             GLOBM(F)= DBLE(-AUW(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+L-1
             GLOBM(F)= DBLE( AUP(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+N+L-1
             GLOBM(F)= DBLE(-AUE(K,L,I,J))
             F=F+1
  42       CONTINUE
           RHS(E+K-1)= DBLE( BU(K,I,J))
  41      CONTINUE
  40    CONTINUE
        IF(IE.GT.IB) THEN
          I=IE
          E=E+N
          DO 43 K=1,N
           DO 44 L=1,N
             IA(F)= E+K-1
             JA(F)= E-W+L-1
             GLOBM(F)= DBLE(-AUS(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E-N+L-1
             GLOBM(F)= DBLE(-AUW(K,L,I,J))
             F=F+1
             IA(F)= E+K-1
             JA(F)= E+L-1
             GLOBM(F)= DBLE( AUP(K,L,I,J))
             F=F+1
  44       CONTINUE
           RHS(E+K-1)= DBLE( BU(K,I,J))
  43      CONTINUE
        ENDIF
      ENDIF
*
      NELT=F-1
*
      RETURN
      END
