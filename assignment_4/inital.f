*
*        file inital.f contains two subroutines: INITAL and SAVE
***********************************************************************
*
      SUBROUTINE INITAL(T,U,UHE,P, T0,U0,UHE0,P0,IRSI,IB,IE,ID)
*
*     Subroutine to set initial T fields by fixing values at T0
*
*     T(ID)   T array for interior and boundary nodes; output
*     T0      initial guess for uniform T-distribution; input
*     U       U array for interior and boundary nodes: output (velocity field)
*     U0      Initialized U value
*     UHE     Velocities evaluated at e face of each CV 
*     UHE0    Initial guess for UHE
*     P       Momentum field for interior and boundary nodes?
*     P0      Initial momentum guess? 
*     
*     IRSI    file number for restart data file; input
*     IB,IE   first and last interior indices in i; input
*     ID      array dimensions; input
*
***********************************************************************
*
      REAL T(ID),T0,
     C     U(ID),U0,
     C     UHE(ID),UHE0,
     C     P(ID),P0
      INTEGER IRSI,IB,IE,ID,I 
*
      DO 10 I=IB-1,IE+1
         T(I)= T0
         U(I) = U0            !Initialize flow field
         UHE(I) = UHE0
         P(I) = P0
 10   CONTINUE
*
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE SAVE(T,U,UHE,P, IRSO,IB,IE,ID)
*
*     Routine to save solution arrays in unformatted form to restart
*     code.
*          
*     T(ID)   T array to be saved;  input
*
*     IRSO    unit number of file for saved data; input
*     IB,IE   first and last interior indices in i; input
*     ID      array dimensions; input
*
************************************************************************
*
      REAL T(ID),
     C     U(ID),UHE(ID),P(ID)
      INTEGER IRSO,IB,IE,ID,I
*     
      WRITE(IRSO) (T(I),I=IB-1,IE+1)
      WRITE(IRSO) (U(I),I=IB-1,IE+1)
      WRITE(IRSO) (UHE(I),I=IB-1,IE+1)
      WRITE(IRSO) (P(I),I=IB-1,IE+1)
*     
      RETURN
      END
