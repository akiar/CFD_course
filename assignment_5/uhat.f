*
*                    file uhat.f
**********************************************************************
*
      SUBROUTINE UHAT(UHE, U,P,
     C                DHUE,DPDX,RHO,
     C                DIEP,DISE,IB,IE,JB,JE,ID,JD)
*     
*     Subroutine to calculate Uhe at the C.V. faces.
*
**********************************************************************
*
      REAL UHE(ID,JD),U(ID,JD),P(ID,JD)
      REAL DHUE(ID,JD),DPDX(ID,JD)
      REAL RHO,DIEP(ID),DISE(ID)
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*
      DO 1 J=JB,JE
       UHE(IB-1,J)= U(IB-1,J)
       DO 2 I=IB,IE-1
	     UHE(I,J)= (U(I,J)+U(I+1,J))/2.0
     C             -DHUE(I,J)*( (P(I+1,J)-P(I,J))/DIEP(I)
     C                         -(DPDX(I,J)+DPDX(I+1,J))/2.0 )
  2    CONTINUE
       UHE(IE,J)= U(IE+1,J)
  1   CONTINUE
*
      RETURN
      END

