*             file dhat.f
*     
************************************************************************
*
      SUBROUTINE DHAT(DHUE, AUP,VOLP,RHO,IB,IE,ID)
*
*     Subroutine to calculate the diffusion coefficients for  
*     the special momentum equation required to compute UHAT.
*
*     DHUE(ID)  diffusion coefficient for east face; output
*
*     AUP(ID)   active coefficient on node P; input
*     IB,IE     first and last interior indices in i; input
*     ID        array dimensions; input     
*
***********************************************************************
*
      REAL DHUE(ID)
      REAL AUP(2,2,ID),VOLP(ID),RHO
      REAL AE,VE
      INTEGER IB,IE,ID,I
*
*      PRINT *, "DHUE AE VE"
      DO 10 I=IB,IE-1
        AE = 0.5*(AUP(2,2,I)+AUP(2,2,I+1))
        VE = 0.5*(VOLP(I)+VOLP(I+1))
        DHUE(I) = VE/AE
*        PRINT *, I,DHUE(I), AE, VE, VOLP(I), AUP(2,2,I)
 10	  CONTINUE
*
      RETURN
      END