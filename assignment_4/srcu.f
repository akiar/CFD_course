*            file srcu.f
*
******************************************************************
      SUBROUTINE SRCU(QU,RU, U,UHE,RHO,VISC,XE,ARO,AREP,IB,IE,ID,
     C                YNE,ZNE,DCCE)
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
      REAL QU(ID),RU(ID)              ! Output source terms
      REAL UHE(ID),U(ID),XE(ID),ARO(ID),AREP(ID)
      REAL RHO,VISC
      REAL DH,CF
      REAL RE,PERIM,YNE(ID),ZNE(ID),DCCE(ID)
      INTEGER IB,IE,ID
*
*      PRINT *, "QU RU PERIM DH RE CF"
      DO 10 I=IB,IE
        PERIM = 2*YNE(I)+2*ZNE(I)                 ! CHANGE BASED ON GEOMETERY OF THE DUCT
        DH = 4*AREP(I)/PERIM
        RE = RHO*DH*ABS(U(I))/VISC
        CF = 0 !(1.58*LOG(RE)-3.28)**(-2) ! QUESTION 3: 0
        QU(I) = DCCE(I-1) - DCCE(I) + 0.5*CF*RHO*ARO(I)*U(I)**2      !QUESTION 1 NO SOURCE
        RU(I) = - CF*RHO*ARO(I)*U(I)                    !QUESTION 1 NO SOURCE -ve
*        PRINT *, QU(I), RU(I), PERIM, DH, RE, CF
 10	  CONTINUE
      RETURN
      END