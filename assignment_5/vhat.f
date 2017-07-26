*
*                    file vhat.f
**********************************************************************
*
      SUBROUTINE VHAT(VHN, V,P,
     C                DHVN,DPDY,RHO,
     C                DJNP,DISN,IB,IE,JB,JE,ID,JD)
*     
*     Subroutine to calculate Uhe at the C.V. faces.
*
**********************************************************************
*
      REAL VHN(ID,JD),V(ID,JD),P(ID,JD)
      REAL DHVN(ID,JD),DPDY(ID,JD)
      REAL RHO,DJNP(JD),DISN(JD)
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*
      DO 1 I=IB,IE
       VHN(I,JB-1)= V(I,JB-1)
       DO 2 J=JB,JE-1
	     VHN(I,J)= (V(I,J)+V(I,J+1))/2.0
     C              -DHVN(I,J)*( (P(I,J+1)-P(I,J))/DJNP(J)
     C                          -(DPDY(I,J)+DPDY(I,J+1))/2.0 )
  2    CONTINUE
       VHN(I,JE)= V(I,JE+1)
  1   CONTINUE
*
      RETURN
      END

