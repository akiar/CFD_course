*
*        file srct.f
******************************************************************
*
      SUBROUTINE SRCT(QT,RT, T,VOLP,ARO,HCONV,TINFC,IB,IE,ID)
*
*     Subroutine to calculate the net source of T in each interior
*     control volume/ unit volume.
*     Net Source = QT + RT*T (energy source/unit volume/cp )
*
*     REAL QT(ID) fixed source coefficient; output
*     REAL RT(ID) linearized source coefficient; output
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER ID,JD array dimensions; input
*
*     NOTE: For solution with radiation, temperature must be
*           absolute!!!
*
******************************************************************
*     Variable Declaration
      REAL T(ID),VOLP(ID),ARO(ID),    ! Inputs
     C     HCONV,TINFC,               ! Inputs
     C     QT(ID),RT(ID)              ! Outputs
      INTEGER IB,IE,ID,               ! Inputs
     C        I                       ! loop integer
******************************************************************
*     Source term calculation
*
      IF (HCONV /= 0) THEN
          DO 10 I=IB,IE
              QT(I) = HCONV * ARO(I) / TINFC
              RT(I) = -HCONV * ARO(I)
   10     CONTINUE
*
      ELSE
          DO 20 I = IB,IE
              QT(I) = 0
              RT(I) = 0
   20     CONTINUE
      END IF
      
      RETURN
      END
