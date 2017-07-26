*
*       file meana.f
*****************************************************************************
*
      SUBROUTINE MEANA(MEAN, PHI,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to find the magnitude of the average of PHI.
*
*****************************************************************************
*
      REAL MEAN,PHI(ID,JD),SUM
      INTEGER IB,IE,JB,JE,NPTS,I,J,ID,JD
*
      NPTS= (IE-IB+1)*(JE-JB+1)
      SUM= 0.0
*
      DO 1 J=JB,JE
        DO 2 I=IB,IE
          SUM= SUM+ABS(PHI(I,J))
  2     CONTINUE
  1   CONTINUE
      MEAN= SUM/NPTS
      IF(MEAN.LT.1.E-20) MEAN= 1.0
*
      RETURN
      END
