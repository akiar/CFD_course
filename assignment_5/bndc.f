*
*     This file contains 4 subroutines: BNDCT, BNDCP, BNDCU and BNDCV
*
************************************************************************
*
      SUBROUTINE BNDCT(ATP,ATW,ATE,ATS,ATN,BT,
     C                 IB,IE,JB,JE,ID,JD)
*
*     Subroutine to put the boundary condition information for each
*     P at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
* 
*            2) A boundary condition must be set for pressure in 
*            all flow domains.  In enclosures, a position inside the
*            domain must be selected.
*
************************************************************************
*
      REAL ATP(ID,JD),ATW(ID,JD),ATE(ID,JD)
      REAL ATS(ID,JD),ATN(ID,JD),BT(ID,JD)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,IBM1,IEP1,JBM1,JEP1
*
      IBM1=IB-1
      JBM1=JB-1
      IEP1=IE+1
      JEP1=JE+1
*
      DO 10 J=JB,JE
*
*     West face boundary conditions
*
        I=IBM1
        ATW(I,J)=0.0
        ATE(I,J)=0.0  !Q2 Q1: 1 Q3: 0
        ATS(I,J)=0.0
        ATN(I,J)=0.0
        ATP(I,J)=1.0
        BT(I,J)= 313.15   !Q1 Q2: 0 Q3: 313.15
*
*     East face boundary conditions
*
        I=IEP1
        ATW(I,J)=0.0  !Q2 Q1: 1 Q3: 0
        ATE(I,J)=0.0
        ATS(I,J)=0.0
        ATN(I,J)=0.0
        ATP(I,J)=1.0
        BT(I,J)= 313.15   !Q2 Q1: 0 Q3: 313.15
 10   CONTINUE
*
      DO 20 I=IB,IE
*
*     South face boundary conditions
*
        J=JBM1
        ATW(I,J)=0.0
        ATE(I,J)=0.0
        ATS(I,J)=0.0
        ATN(I,J)=0.0
        ATP(I,J)=1.0
        BT(I,J)= 293.15    !Q2 360 Q3 293.15
*
*     North face boundary conditions
*
        J=JEP1
        ATW(I,J)=0.0
        ATE(I,J)=0.0
        ATS(I,J)=1.0  !Q1 Q3 1 Q2 0
        ATN(I,J)=0.0
        ATP(I,J)=1.0
        BT(I,J)= 0.0  !Q1 Q3 0 Q2 300
*
 20   CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCP(AUP,AUW,AUE,AUS,AUN,BU,
     C                 P,DIEP,DJNP,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to put the boundary condition information for each
*     P at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
* 
*            2) A boundary condition must be set for pressure in 
*            all flow domains.  In enclosures, a position inside the
*            domain must be selected.
*
************************************************************************
*
      REAL AUP(N,N,ID,JD),AUW(N,N,ID,JD),AUE(N,N,ID,JD)
      REAL AUS(N,N,ID,JD),AUN(N,N,ID,JD),BU(N,ID,JD),P(ID,JD)
      REAL DIEP(ID),DJNP(JD)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,K,IBM1,IEP1,JBM1,JEP1
*
      K=1
      IBM1=IB-1
      JBM1=JB-1
      IEP1=IE+1
      JEP1=JE+1
*
      DO 10 J=JB,JE
*
*     West face boundary conditions
*
        I=IBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=1.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)=-DIEP(IB-1)*(P(IB+1,J)-P(IB,J))/DIEP(IB)
*
*     East face boundary conditions
*
        I=IEP1
        AUW(K,K,I,J)=1.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= DIEP(IE)*(P(IE,J)-P(IE-1,J))/DIEP(IE-1)
 10   CONTINUE
*
      DO 20 I=IB,IE
*
*     South face boundary conditions
*
        J=JBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=1.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)=-DJNP(JB-1)*(P(I,JB+1)-P(I,JB))/DJNP(JB) 
*
*     North face boundary conditions
*
        J=JEP1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0  !Q1 Q2 1 Q3 0 
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0 !DJNP(JE)*(P(I,JE)-P(I,JE-1))/DJNP(JE-1)
*
 20   CONTINUE
*
*  Pressure level for some point in the domain (be careful!)
*     ONLY FOR ENCLOSURES
*
*      I=INT(IEP1/2)
*      J=INT(JEP1/2)
*      DO 30 L=1,3
*        AUW(K,L,I,J)=0.0
*        AUE(K,L,I,J)=0.0
*        AUS(K,L,I,J)=0.0
*        AUN(K,L,I,J)=0.0
*        AUP(K,L,I,J)=0.0
*        BU(K,I,J)= 0.0
* 30   CONTINUE
*      AUP(K,K,I,J)= 1.0
*      BU(K,I,J)= 0.0
*     
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCU(AUP,AUW,AUE,AUS,AUN,BU,
     C                 IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to put the boundary condition information for U
*     at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
*
************************************************************************
*
      REAL AUP(N,N,ID,JD),AUW(N,N,ID,JD),AUE(N,N,ID,JD)
      REAL AUS(N,N,ID,JD),AUN(N,N,ID,JD),BU(N,ID,JD)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,K,IBM1,IEP1,JBM1,JEP1
*
      K=2
      IBM1=IB-1
      JBM1=JB-1
      IEP1=IE+1
      JEP1=JE+1
*
      DO 10 J=JB,JE
*
*     West face boundary conditions
*
        I=IBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0
*
*     East face boundary conditions
*
        I=IEP1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)=0.0
 10   CONTINUE
*
      DO 20 I=IB,IE
*
*     South face boundary conditions
*
        J=JBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=1.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0 
*
*     North face boundary conditions
*
        J=JEP1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=1.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0
*
 20   CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCV(AUP,AUW,AUE,AUS,AUN,BU,
     C                 IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to put the boundary condition information for V
*     at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
*
************************************************************************
*
      REAL AUP(N,N,ID,JD),AUW(N,N,ID,JD),AUE(N,N,ID,JD)
      REAL AUS(N,N,ID,JD),AUN(N,N,ID,JD),BU(N,ID,JD)
      INTEGER IB,IE,JB,JE
      INTEGER I,J,K,IBM1,IEP1,JBM1,JEP1
*
      K=3
      IBM1=IB-1
      JBM1=JB-1
      IEP1=IE+1
      JEP1=JE+1
*
      DO 10 J=JB,JE
*
*     West face boundary conditions
*
        I=IBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0
*
*     East face boundary conditions
*
        I=IEP1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0
 10   CONTINUE
*
      DO 20 I=IB,IE
*
*     South face boundary conditions
*
        J=JBM1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=0.0
        AUN(K,K,I,J)=1.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)= 0.0 
*
*     North face boundary conditions
*
        J=JEP1
        AUW(K,K,I,J)=0.0
        AUE(K,K,I,J)=0.0
        AUS(K,K,I,J)=1.0
        AUN(K,K,I,J)=0.0
        AUP(K,K,I,J)=1.0
        BU(K,I,J)=0.0
 20   CONTINUE
*
      RETURN
      END
