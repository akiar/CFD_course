*             file dhat.f
*     
************************************************************************
*
      SUBROUTINE DHAT(DHUE,DHVN, AUP,DIEP,DISE,DJNP,DISN,
     C                VOLP,RHO,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to calculate the diffusion coefficients for  
*     the special momentum equation required to compute the 
*     advecting velocity, UHE.
*
*     DHUE(ID,JD)    diffusion coefficient for east face; output
*     DHVN(ID,JD)    diffusion coefficient for north face; output
*
*     IB,IE        first and last interior indices in i; input
*     JB,JE        first and last interior indices in j; input
*     ID,JD        array dimensions; input     
*
***********************************************************************
*
      REAL DHUE(ID,JD),DHVN(ID,JD),AUP(N,N,ID,JD)
      REAL DIEP(ID),DISE(ID),DJNP(JD),DISN(JD)
      REAL VOLP(ID,JD),RHO
      INTEGER IB,IE,JB,JE,I,J,N,ID,JD
*
*  compute diffusion coefficients for i direction
*
      DO 10 J=JB,JE
       DO 20 I=IB,IE-1
	    DHUE(I,J)= ( VOLP(I  ,J)+VOLP(I+1,J) )/
     C             ( AUP(2,2,I  ,J)+AUP(2,2,I+1,J)+1.E-20 )
  20   CONTINUE
  10  CONTINUE
*
*  compute diffusion coefficients for j direction
*
      DO 30 I=IB,IE
       DO 40 J=JB,JE-1
	    DHVN(I,J)= ( VOLP(I,J  )+VOLP(I,J+1) )/
     C             ( AUP(3,3,I,J  )+AUP(3,3,I,J+1)+1.E-20 )
  40   CONTINUE
  30  CONTINUE
*
      RETURN
      END

