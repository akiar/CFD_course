*                  file bndcu.f
********************************************************************
*
      SUBROUTINE BNDCU(AUP,AUW,AUE,BU, IB,IE,ID)
*
*     Subroutine to put the boundary condition information for U
*     at each boundary node into equation coefficients.
*
*     AUP(2,2,ID)   active coefficient for P node; output
*     AUW(2,2,ID)   active coefficient for W node; output
*     AUE(2,2,ID)   active coefficient for E node; output
*     BU(2,ID)    accumulated fixed source term; output
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER ID     array dimensions; input
*
*
********************************************************************
*
      REAL ...
      INTEGER ...
*
      RETURN
      END
