*            file uhat.f
*******************************************************************
*
      SUBROUTINE UHAT(UHE, U,P,DHUE,RHO,DPDX,DIEP,IB,IE,ID)
*     
*     Subroutine to calculate velocity to use in advecting term.
*
*     UHE(ID)    advecting velocity; output
*
*     U(ID)      current velocity at nodal points; input
*     P(ID)      current pressure at nodal points; input
*     DHUE(ID)   coefficients for UHE(ID); input
*     DIEP(ID)   x distance between nodes P and E; input
*     RHO        fluid density
*     DTIME      time step; input
*     DPDX(ID)   nodal values of pressure gradient
*
*     IB,IE      first and last interior indices; input
*     ID         array size; input
*
*******************************************************************
*
      REAL UHE(ID)
      REAL U(ID),P(ID),DHUE(ID),DPDX(ID)
      REAL RHO,DIEP(ID)
      INTEGER IB,IE,ID,I
*
*      PRINT *, "UHE"
      DO 10 I=IB,IE-1
        UHE(I) = (U(I)+U(I+1))/2 - DHUE(I)*((P(I+1)-P(I))/DIEP(I)
     C                                    - 0.5*(DPDX(I)+DPDX(I+1)))
*        PRINT *, I,UHE(I)
 10	  CONTINUE     
*
      RETURN
      END

