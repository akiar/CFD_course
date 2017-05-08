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
      IF (HCONV /= 0) THEN    ! Loop if convection is present
          DO 10 I=IB,IE       ! Convection on each outer face
              QT(I) = 50000 * VOLP(I) ! if internal heat gen and no conv. 
              RT(I) = 0               ! if internal heat gen and no conv. 
*
              !QT(I) = HCONV * ARO(I) * TINFC ! if internal faces are exposed to surroundings
              !RT(I) = -HCONV * ARO(I)        ! if internal faces are exposed to surroundings
*
              IF (RT(I) < 0) THEN   ! Ensure RT is positive
                  RT(I) = -RT(I)
              ELSE
                  RT(I) = RT(I)
              END IF
   10     CONTINUE
*
      ELSE                    ! If no convection present set to 0
          DO 20 I = IB,IE     
              QT(I) = 0
              RT(I) = 0
   20     CONTINUE
      END IF
      
      RETURN
      END
