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
      REAL ...
      INTEGER ...
*
      RETURN
      END

