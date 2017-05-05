*
*        file bndct.f
***********************************************************************
*
      SUBROUTINE BNDCT(ATP,ATW,ATE,BT, DE,AREP,IB,IE,ID)
*
*     Subroutine to put the boundary condition information for T
*     at each boundary node into equation coefficients.
*
*     REAL ATP(ID) active coefficient for P node; output
*     REAL ATW(ID) active coefficient for W node; output
*     REAL ATE(ID) active coefficient for E node; output
*     REAL BT(ID)  accumulated fixed source term; output
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER ID     array dimensions; input
*
***********************************************************************
*
*     Variable Declaration
*
      REAL ATP(ID),ATW(ID),ATE(ID),BT(ID)
      INTEGER IB,IE,ID
*
***********************************************************************
*
*     Boundary Condition Calculations
*     Set coefficients for Dirichlet, Neuman, and Robin conditions
*      
*     Beginning Node
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Tspec
*     Neumann:    AW = 0, AE = 1, AP = 1, BP = Qspec/DE(IB-1)
*     Robin:      AW = 0, AE = DE(IB-1), AP = HCONV*ARE(IB-1) + DE(IB-1),
*                         BP = HCONV*ARE(IB-1)*TINFC
*                         
      ATW(IB-1) = 0
      ATE(IB-1) = 0
      ATP(IB-1) = 1
      BT(IB-1) = 100 + 273.15
*
*     End Node
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Tspec
*     Neumann:    AW = 1, AE = 0, AP = 1, BP = -Qspec/DE(IE)
*     Robin:      AW = 0, AE = DE(IE),  AP = HCONV*ARE(IE) + DE(IE),
*                         BP = HCONV*ARE(IB-1)*TINFC
*
      ATW(IE+1) = 0
      ATE(IE+1) = 0
      ATP(IE+1) = 1
      BT(IE+1) = 0 + 273.15
*     
      RETURN
      END
