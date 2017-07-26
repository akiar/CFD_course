*
*     This file contains 2 subroutines: WRT2D and WRT2DV
************************************************************************
*
      SUBROUTINE WRT2D(AP,AW,AE,AS,AN,B,IDATO,
     C                 IV1,IV2,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to write active coefficients selected from the
*     blocks AP...B using IV.
*
************************************************************************
*
      PARAMETER(M=200)
      REAL AP(N,N,ID,JD),AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),B(N,ID,JD)
      REAL ATP(M,M),ATW(M,M),ATE(M,M)
      REAL ATS(M,M),ATN(M,M),BB(M,M)
      INTEGER IDATO,IV1,IV2,IB,IE,JB,JE,I,J,N,ID,JD
*
      DO 1 J=JB-1,JE+1
        DO 2 I=IB-1,IE+1
          ATP(I,J)= AP(IV1,IV2,I,J)
          ATW(I,J)= AW(IV1,IV2,I,J)
          ATE(I,J)= AE(IV1,IV2,I,J)
          ATS(I,J)= AS(IV1,IV2,I,J)
          ATN(I,J)= AN(IV1,IV2,I,J)
          BB(I,J)=  B(IV1,I,J)
  2     CONTINUE
  1   CONTINUE
*
      CALL OUT2D(ATP,'AP      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
      CALL OUT2D(ATW,'AW      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
      CALL OUT2D(ATE,'AE      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
      CALL OUT2D(ATS,'AS      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
      CALL OUT2D(ATN,'AN      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
      CALL OUT2D(BB ,'B       ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE WRT2DV(UU, IJ1,IJ2,TIT,IDATO,IB,IE,JB,JE,IJ,ID,JD)
*
*     Subroutine to write a particular Reynolds-stress
*     UU(IJ1,IJ2,ID,JD).
*
***********************************************************************
*
      PARAMETER(M=200)
      REAL UU(IJ,IJ,ID,JD),A(M,M)
      CHARACTER*8 TIT
      INTEGER IJ1,IJ2,IDATO,IB,IE,JB,JE,IJ,ID,JD
      INTEGER I,J
*
      DO 1 I=IB-1,IE+1
        DO 2 J=JB-1,JE+1
          A(I,J)= UU(IJ1,IJ2,I,J)
  2     CONTINUE
  1   CONTINUE
*
      CALL OUT2D(A,TIT,IDATO,IB-1,IE+1,1,JB-1,JE+1,1,M,M)
*
      RETURN
      END

