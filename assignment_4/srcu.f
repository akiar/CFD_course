*            file srcu.f
*
******************************************************************
      SUBROUTINE SRCU(QU,RU, U,UHE,RHO,VISC,XE,ARO,AREP,IB,IE,ID)
*
*     Subroutine to calculate the net source of U in each interior
*     control volume:
*     Net Source = QU + RU*U 
*
*     QU(ID) fixed source coefficient; output
*     RU(ID) linearized source coefficient; output
*
*     IB,IE  first and last interior indices in i; input
*     ID,JD  array dimensions; input
*
*     DH     hydraulic diameter, Dh= 4A/P; computed
*     CF     drag coefficient; computed
*
******************************************************************
      REAL ...
      INTEGER ...
*
      RETURN
      END

