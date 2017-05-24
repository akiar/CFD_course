*
*        file bndct.f
************************************************************************************
*
      SUBROUTINE BNDCT(ATP,ATW,ATE,BT, DE,AREP,IB,IE,ID,HCONV,TINF)
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
************************************************************************************
*
*     Variable Declaration
*
      REAL DE(ID),ATP(ID),ATW(ID),ATE(ID),BT(ID),HCONV,AREP(ID)
      INTEGER IB,IE,ID
*
************************************************************************************
*
*     Boundary Condition Calculations
*     Set coefficients for Dirichlet, Neuman, and Robin conditions
*
*     Beginning Node - 1
*      Node that has zero volume and only contains boundary information
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Tspec
*     Neumann:    AW = 0, AE = 1, AP = 1, BP = Qspec/DE(IB-1)
*     Robin:      AW = 0, AE = DE(IB-1), AP = HCONV*AREP(IB-1) + DE(IB-1),
*                         BP = HCONV*AREP(IB-1)*TINFC
*
      ATW(IB-1) = 0
      ATE(IB-1) = DE(IB-1) ! QUESTION 3 CONDITION
      ATP(IB-1) = 1 * AREP(IB-1) + DE(IB-1) ! QUESTION 3 CONDITION
      BT(IB-1) = 0 ! HCONV * AREP(IB-1) * 323.15 ! QUESTION 3 CONDITION
*
*     End Node + 1
*      Node that has zero volume and only contains boundary information
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Tspec
*     Neumann:    AW = 1, AE = 0, AP = 1, BP = -Qspec/DE(IE)
*     Robin:      AW = 0, AE = DE(IE),  AP = HCONV*AREP(IE) + DE(IE),
*                         BP = HCONV*AREP(IE)*TINFC
*
      ATW(IE+1) = 0 ! QUESTION 3 CONDITION
      ATE(IE+1) = DE(IE)
      ATP(IE+1) = 1 * AREP(IE) + DE(IE) ! QUESTION 3 CONDITION
      BT(IE+1) = 0 ! HCONV * AREP(IE) * 313.15 ! QUESTION 3 CONDITION
*     
      RETURN
      END
