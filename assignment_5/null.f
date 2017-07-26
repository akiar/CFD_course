*
*      This file contains 2 subroutines: NULLM and NULLV
*
************************************************************************
       SUBROUTINE NULLM(A, ISTRT,IFIN,JSTRT,JFIN,ID,JD)
*
*      Subroutine to null a matrix.
*
*      A(ID,JD)    matrix
*      ISTRT,IFIN  beginning and ending indices in I direction
*      JSTRT,JFIN  beginning and ending indices in J direction
*      ID,JD       array dimensions
*
************************************************************************
*
      REAL A(ID,JD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,I,J,ID,JD
*
      DO 1 I=ISTRT,IFIN
        DO 2 J=JSTRT,JFIN
          A(I,J)= 0.0
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE NULLMM(A, ISTRT,IFIN,JSTRT,JFIN,N,ID,JD)
*
*     Subroutine to null a block matrix.
*
*     A(N,ID,JD)  matrix
*     ISTRT,IFIN  beginning and ending indices in I direction
*     JSTRT,JFIN  beginning and ending indices in J direction
*     ID,JD       array dimensions
*
************************************************************************
*
      REAL A(N,ID,JD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,I,J,K,N,ID,JD
*
      DO 1 K=1,N
        DO 2 I=ISTRT,IFIN
          DO 3 J=JSTRT,JFIN
            A(K,I,J)= 0.0
  3       CONTINUE
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END
*
************************************************************************
       SUBROUTINE NULLV(A, ISTRT,IFIN,ID)
*
*      Subroutine to null a vector.
*
*      A(ID)     vector
*      ISRT,IFIN beginning and ending indices
*      ID        array dimension
*
************************************************************************
*
      REAL A(ID)
      INTEGER ISTRT,IFIN,I,ID
*
      DO 1 I=ISTRT,IFIN
        A(I)= 0.0
  1   CONTINUE
*
      RETURN
      END
*
************************************************************************
       SUBROUTINE NULLMN(A, ISTRT,IFIN,JSTRT,JFIN,N,ID,JD)
*
*      Subroutine to null a block matrix.
*
*      A(N,N,ID,JD)    matrix
*      ISTRT,IFIN  beginning and ending indices in I direction
*      JSTRT,JFIN  beginning and ending indices in J direction
*      ID,JD       array dimensions
*
************************************************************************
*
      REAL A(N,N,ID,JD)
      INTEGER ISTRT,IFIN,JSTRT,JFIN,N,ID,JD
      INTEGER I,J,K,L
*
      DO 1 K=1,N
        DO 2 L=1,N
          DO 3 I=ISTRT,IFIN
            DO 4 J=JSTRT,JFIN
              A(K,L,I,J)= 0.0
  4         CONTINUE
  3       CONTINUE
  2     CONTINUE
  1   CONTINUE
*
      RETURN
      END

