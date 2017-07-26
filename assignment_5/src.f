*
*     This file contains 4 subroutines: SRCT, SRCU, SRCV and SRCDC
*
************************************************************************
      SUBROUTINE SRCT(QT,RT, DCCE,DCCN,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the net source of T in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*T 
*
*     Q(ID,JD) fixed source coefficient; output
*     R(ID,JD) linearized source coefficient; output
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input
*
************************************************************************
*
      REAL QT(ID,JD),RT(ID,JD)
      REAL DCCE(ID,JD),DCCN(ID,JD)
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*     
      DO 20 J=JB,JE
      DO 10 I=IB,IE
        QT(I,J)= DCCE(I-1,J) - DCCE(I,J) + DCCN(I,J-1) - DCCN(I,J)
        RT(I,J)= 0.0
 10   CONTINUE
 20   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCU(QU,RU, DCCE,DCCN,T,RHO,VOLP,GEX,BETA,
     C                IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the net source of U in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*U 
*
*     Q(ID,JD) fixed source coefficient; output
*     R(ID,JD) linearized source coefficient; output
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input
*
************************************************************************
*
      PARAMETER(TREF= 293.15) !Q3: 293.15
      REAL QU(ID,JD),RU(ID,JD)
      REAL DCCE(ID,JD),DCCN(ID,JD)
      REAL T(ID,JD),VOLP(ID,JD),RHO,GEX,BETA
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*     
      DO 20 J=JB,JE
      DO 10 I=IB,IE
        QU(I,J)= DCCE(I-1,J) - DCCE(I,J) + DCCN(I,J-1) - DCCN(I,J)
     C           -RHO*BETA*GEX*(T(I,J)-TREF)*VOLP(I,J) 
        RU(I,J)= 0.0
 10   CONTINUE
 20   CONTINUE
      RETURN
      END
*
************************************************************************
      SUBROUTINE SRCV(QV,RV, DCCE,DCCN,T,RHO,VOLP,GEY,BETA,
     C                IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the net source of V in each interior
*     control volume for the entire volume.
*     Net Source = Q + R*V 
*
*     Q(ID,JD) fixed source coefficient; output
*     R(ID,JD) linearized source coefficient; output
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input
*
************************************************************************
*
      PARAMETER(TREF= 293.15)
      REAL QV(ID,JD),RV(ID,JD)
      REAL DCCE(ID,JD),DCCN(ID,JD)
      REAL T(ID,JD),VOLP(ID,JD),RHO,GEY,BETA
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*     
      DO 20 J=JB,JE
      DO 10 I=IB,IE
        QV(I,J)= DCCE(I-1,J) - DCCE(I,J) + DCCN(I,J-1) - DCCN(I,J)
     C           -RHO*BETA*GEY*(T(I,J)-TREF)*VOLP(I,J) 
        RV(I,J)= 0.0
 10   CONTINUE
 20   CONTINUE
      RETURN
      END
      