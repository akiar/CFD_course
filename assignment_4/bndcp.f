*                  file bndcp.f
********************************************************************
*
      SUBROUTINE BNDCP(AUP,AUW,AUE,BU, P,XP,IB,IE,ID,DE)  !DIEP TO XP
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
      REAL AUP(2,2,ID),AUW(2,2,ID),AUE(2,2,ID),BU(2,ID)
      REAL P(ID),XP(ID),DE(ID)
      INTEGER IB,IE,ID
      PRINT *, "AUW(1,1) AUP(1,1) AUE(1,1) BU(1)"
*
*     BEGINNING NODE
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Tspec
*     Neumann:    AW = 0, AE = 1, AP = 1, BP = Qspec/DE(IB-1)
*     Robin:      AW = 0, AE = DE(IB-1), AP = HCONV*AREP(IB-1) + DE(IB-1),
*                         BP = HCONV*AREP(IB-1)*TINFC
*
      AUW(1,1,IB-1) = 0
      AUP(1,1,IB-1) = 1
      AUE(1,1,IB-1) = 1
      BU(1,IB-1)= 0 !(P(IB)+(XP(IB-1)-XP(IB))/                   !Q3: 0  
*     C            (XP(IB+1)-XP(IB))*(P(IB+1)-P(IB)))-P(IB)          !Q3: 0
      PRINT *,IB-1,AUW(1,1,IB-1),AUP(1,1,IB-1),AUE(1,1,IB-1), BU(1,IB-1)
      AUW(1,2,IB-1) = 0
      AUP(1,2,IB-1) = 0
      AUE(1,2,IB-1) = 0
*
*     END NODE
*
*     Dirichlet:  AW = 0, AE = 0, AP = 1, BP = Pspec
*     Neumann:    AW = 1, AE = 0, AP = 1, BP = -Pspec/DE(IE)
*     Robin:      AW = DE(IE), AE = 0,  AP = HCONV*AREP(IE) + DE(IE),
*                         BP = HCONV*AREP(IE)*TINFC
*
      AUW(1,1,IE+1) = 0
      AUP(1,1,IE+1) = 1
      AUE(1,1,IE+1) = 0
      BU(1,IE+1) = 0
*      
      PRINT *,IE+1,AUW(1,1,IE+1),AUP(1,1,IE+1),AUE(1,1,IE+1), BU(1,IE+1)
*
      AUW(1,2,IE+1) = 0
      AUP(1,2,IE+1) = 0
      AUE(1,2,IE+1) = 0
      RETURN
      END