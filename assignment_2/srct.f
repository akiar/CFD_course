*     
*        file srct.f
************************************************************************************
*
      SUBROUTINE SRCT(QT,RT, T,VOLP,ARO,HCONV,TINFC,IB,IE,ID,EMIS,OMEG,
     C                DEOLD,TOLD)
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
************************************************************************************
*
*     Variable Declaration
      PARAMETER (SBC=5.67E-8)        ! stefan-boltzmann constant for rad.
      REAL T(ID),VOLP(ID),ARO(ID),   ! Inputs
     C     HCONV,TINFC,              ! Inputs
     C     INTGEN,                   ! Internal heat generation
     C     QT(ID),RT(ID),            ! Outputs
     C     OMEG,DEOLD(ID),TOLD(ID)
      INTEGER IB,IE,ID,              ! Inputs
     C        I
*     C        LIN! loop integer, linearization type
*
*--Linearization and heat generation flags
*
*     Set INTGEN to specified internal heat generation (W/m^3) 
*     Set LIN to 1,2, or 3 for linearization technique
*
      INTGEN = 0      ! set internal heat generation (W/m^3). Q3: INTGEN = 50000
*     LIN = 1         ! choose linearization type: 1, 2, 3 (Newton-raphson)
*
************************************************************************************
*
*     Source term calculation: Loop over all CVs and calculate sources
*     --Loops will adjust based on heat transfer methods present
*         -- EMIS = 0 in.dat / no radiation 
*         -- HCONV = 0 in.dat / no convection
*
*--Begin loop over all CVs
*
      DO 10 I=IB,IE
*     
*--Newton - Raphson Linearization
*
          QT(I) = EMIS * ARO(I) * SBC * (3 * (T(I)**4) + TINFC**4) !Fixed source
     C            + HCONV * ARO(I) * TINFC    !convection source 
     C            + INTGEN * VOLP(I)          !internal gen source
     C            + (1-OMEG)*(DEOLD(I)*(TOLD(I+1)-TOLD(I))
     C                        -DEOLD(I-1)*(TOLD(I)-TOLD(I-1)))
          RT(I) = - 4 * EMIS * ARO(I) * SBC *(T(I)**3)    !linearized source
     C            - HCONV * ARO(I)                        !convection
*
*
*--Linearization 1
*
*        IF (LIN == 1) THEN
*          QT(I) = - EMIS * ARO(I) * SBC *(T(I)**4 - TINFC**4) !Fixed source
*     C            + HCONV * ARO(I) * TINFC                    !Convection source
*     C            + INTGEN * VOLP(I)                          !Internal Gen source
*     C            + (1 - OMEG) * 
*          RT(I) = - HCONV * ARO(I) ! Linearized Source = 0, only CONV ! Q3: RT(I)=0 
*
*--Linearization 2
*
*        ELSE IF  (LIN == 2) THEN
*          QT(I) = EMIS * SBC * ARO(I) * (TINFC**4)    !Fixed source
*     C            + HCONV * ARO(I) * TINFC            !Convection source
*     C            + INTGEN * VOLP(I)                  !Internal Gen source
*          RT(I) = - EMIS * ARO(I) * SBC * (T(I)**3) !Linearized Source ! Q3: RT(I)=0 
*     C            - HCONV * ARO(I)                  ! Convection 
*
 10   CONTINUE
*
      RETURN
      END
