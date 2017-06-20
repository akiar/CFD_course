*            file srcu.f
*
******************************************************************
      SUBROUTINE SRCU(QU,RU, U,UHE,RHO,VISC,XE,ARO,AREP,IB,IE,ID,
     C                YNE,ZNE,UOLD)
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
      REAL RE,PERIM,YNE(ID),ZNE(ID)
      REAL UOLD(ID)
      INTEGER IB,IE,ID
*
      PRINT *, "PERIM DH RE CF QU RU"
      DO 10 I=IB,IE
*        PERIM = 2*YNE(I)+2*ZNE(I)                 ! CHANGE BASED ON GEOMETERY OF THE DUCT
*        DH = 4*AREP(I)/PERIM
*        RE = RHO*DH*U(I)/VISC
*        CF = (1.58*LOG(RE)-3.28)**(-2)
        QU(I) = 0 ! 0.5*CF*RHO*ARO(I)*UOLD(I)**2      !QUESTION 1 NO SOURCE
        RU(I) = 0 ! CF*RHO*UOLD(I)                    !QUESTION 1 NO SOURCE
*        PRINT *, PERIM, DH, RE, CF, QU(I), RU(I)
 10	  CONTINUE
      RETURN
      END

