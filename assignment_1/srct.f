*
*        file srct.f
******************************************************************
*
      SUBROUTINE SRCT(QT,RT, T,VOLP,ARO,HCONV,TINFC,IB,IE,ID,EMIS)
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
************************************************************************
*     Variable Declaration
      PARAMETER (SBC=5.67E-8)        ! stefan-boltzmann constant for rad.
      REAL T(ID),VOLP(ID),ARO(ID),   ! Inputs
     C     HCONV,TINFC,              ! Inputs
     C     INTGEN,                   ! Internal heat generation
     C     QT(ID),RT(ID)             ! Outputs
      INTEGER IB,IE,ID,              ! Inputs
     C        I,LIN                  ! loop integer, linearization type
************************************************************************
*
*     Source term calculation: Loop over all CVs and calculate sources
*     --Set INTGEN to specified internal heat generation (W/m^3) 
*     --Set LIN to 1,2, or 3 for linearization technique
*     --Loops will adjust based on heat transfer methods present
*         -- EMIS = 0 / no radiation
*         -- HCONV = 0 / no convection
*
*--Linearization and heat generation flags
*
      INTGEN = 50000      ! set internal CV heat generation (W/m^3). QUESTION 3: INTGEN = 50000
      LIN = 1         ! choose linearization type: 1, 2, 3 (Newton-raphson) 
*
*--Begin loop over all CVs
*
      DO 10 I=IB,IE
*
*--Linearization 1
*
        IF (LIN == 1) THEN
          QT(I) = - EMIS * ARO(I) * SBC *(T(I)**4 - TINFC**4)
     C            + HCONV * ARO(I) * TINFC
     C            + INTGEN * VOLP(I)
          RT(I) = 0 !- HCONV * ARO(I) ! Question 3: RT(I) = 0 for only internal generation
*
*--Linearization 2
*
        ELSE IF  (LIN == 2) THEN
          QT(I) = EMIS * SBC * ARO(I) * (TINFC**4)
     C            + HCONV * ARO(I) * TINFC
     C            + INTGEN * VOLP(I)
          RT(I) = - EMIS * ARO(I) * SBC * (T(I)**3)
     C            - HCONV * ARO(I)
*     
*--Linearization 3, Newton - Raphson
*
        ELSE
          QT(I) = EMIS * ARO(I) * SBC * (3 * (T(I)**4) + TINFC**4)
     C            + HCONV * ARO(I) * TINFC
     C            + INTGEN * VOLP(I)
          RT(I) = - 4 * EMIS * ARO(I) * SBC *(T(I)**3)
     C            - HCONV * ARO(I)
*
        END IF
 10   CONTINUE

      RETURN
      END
