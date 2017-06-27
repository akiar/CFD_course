*             file gradp.f
*     
************************************************************************
*
      SUBROUTINE GRADP(DPDX, P,DIEP,IB,IE,ID)
*
*     Subroutine to calculate the pressure gradient at each nodal
*     location using the most current pressure field.
*
*     DPDX(ID)  derivative of pressure at nodal location; output
*
*     P(ID)     current pressure field; input
*     DIEP(IE)  distance from P to E nodes; input
*     IB,IE     first and last interior indices in i; input
*     ID        array dimensions; input     
*
************************************************************************
*
      REAL DPDX(ID),DIEP(ID),P(ID)
      INTEGER IB,IE,ID
*
*      PRINT *, "DPDX"
      
      DO 10 I=IB-1,IE
        DPDX(I) = (P(I+1)-P(I))/DIEP(I)
*        PRINT *, I,DPDX(I)
 10	  CONTINUE
*
      RETURN
      END
