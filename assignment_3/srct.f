*     
*        file srct.f
************************************************************************************
*
      SUBROUTINE SRCT(QT,RT, T,VOLP,ARO,HCONV,TINFC,IB,IE,ID,EMIS,
     C                OMEG,DEOLD,TOLD,
     C                DCCE)
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
     C     OMEG,DEOLD(ID),TOLD(ID),
     C     DCCE(ID)                  !deffered corrections 
      INTEGER IB,IE,ID,              ! Inputs
     C        I
*
*--Heat generation flag
*
*     Set INTGEN to specified internal heat generation (W/m^3) 
*
      INTGEN = 0      ! set internal heat generation (W/m^3)
*
************************************************************************************
*
*     Source term calculation: Loop over all CVs and calculate sources
*     --Loops will adjust based on heat transfer methods present
*         -- EMIS = 0 in.dat / no radiation 
*         -- HCONV = 0 in.dat / no convection
*         -- OMEG will change transient discretiztion method, affecting QT
*
*--Begin loop over all CVs
*
      DO 10 I=IB,IE
*
*--Newton - Raphson Linearization
*
          QT(I) = EMIS * ARO(I) * SBC * (3 * (T(I)**4) + TINFC**4) !Fixed source
     C            + HCONV * ARO(I) * TINFC                    !convection source 
     C            + INTGEN * VOLP(I)                          !internal gen source
     C            + (1-OMEG)*(DEOLD(I)*(TOLD(I+1)-TOLD(I))    !Transient Term
     C                        -DEOLD(I-1)*(TOLD(I)-TOLD(I-1)))
     C            - DCCE(I) + DCCE(I-1)                       ! Deferred correction

          RT(I) = - 4 * EMIS * ARO(I) * SBC *(T(I)**3)    !linearized source
     C            - HCONV * ARO(I)                        !convection
*	                                                      !No trasient effect
*                                                         !no deferred correction
*
 10   CONTINUE
*
      RETURN
      END