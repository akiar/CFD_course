*
*        file coeff.f
***********************************************************************
*
      SUBROUTINE COEFF(AP,AW,AE,B,
     C                 DE,Q,R,VOLP,RHO,CP,
     C                 IB,IE,ID)
*
*     Subroutine to calculate the transport coefficients for the
*     variable PHI on a co-located non-orthogonal grid.
*
*     REAL AP(ID) active coefficient for P node; output
*     REAL ASUM(ID) sum of active neighbour coefficients; output
*     REAL AW(ID) active coefficient for W node; output
*     REAL AE(ID) active coefficient for E node; output
*     REAL B(ID)  accumulated fixed source term; output
*
*     REAL DE(ID) diffusion coefficient for east face; input
*     REAL Q(ID) fixed source term; input
*     REAL R(ID) linearized source term; input
*     REAL PHI(ID) estimate of phi at old time; input
*     REAL VOLP(ID) c.v. volume; input
*     REAL RHO  fluid density; input
*     REAL CP specific heat of fluid; input
*     REAL DTIME  time step for implicit scheme; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER ID array dimensions; input     
*
**********************************************************************
*
*     Declaration of variables
*
      REAL AP(ID),ASUM(ID),AW(ID),AE(ID),B(ID),   ! coefficient declaration
     C     DE(ID),Q(ID),R(ID),                    ! source term declaration
     C     PHI(ID),VOLP(ID),DTIME,                ! time variable declaration
     C     RHO,CP                                 ! properties declaration
      INTEGER IB,IE,ID,   ! first, last interior indices, length of array
     C        I           ! loop integer
*
**********************************************************************
*
*     Transport coefficient calculations
*     Assume eqation of the form:
*         AP*TP = AW*TW + AE*TE + bP
*         TW = T(I-1)
*
      DO 10 I = IB,IE
          AW(I) = DE(I-1)
          AE(I) = DE(I)
          ASUM(I) = AW(I) + AE(I)
          AP(I) = ASUM(I) + R(I)
          B(I) = Q(I)
          
   10 CONTINUE
      RETURN
      END
