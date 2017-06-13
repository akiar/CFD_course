*                  file bndcp.f
********************************************************************
*
      SUBROUTINE BNDCP(AUP,AUW,AUE,BU, P,DIEP,IB,IE,ID)
*
*     Subroutine to put the boundary condition information for P
*     at each boundary node into equation coefficients.
*
*     AUP(2,2,ID)   active coefficient for P node; output
*     AUW(2,2,ID)   active coefficient for W node; output
*     AUE(2,2,ID)   active coefficient for E node; output
*     BU(2,ID)    accumulated fixed source term; output
*
*     IB,IE     first and last interior indices in i; input
*     ID        array dimensions; input
*
*
********************************************************************
*
      REAL ...
      INTEGER ...
*
      RETURN
      END

