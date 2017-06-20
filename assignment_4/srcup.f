*
*            file srcup.f
******************************************************************
      SUBROUTINE SRCUP(AUP,AUW,AUE,BU, DIEP,VOLP,IB,IE,ID)
*
*     Subroutine to add the pressure gradient source term to the 
*     momentum equation for each interior control volume:
*
*     AUP(2,2,ID) active coefficient for P node; output
*     AUW(2,2,ID) active coefficient for W node; output
*     AUE(2,2,ID) active coefficient for E node; output
*     BU(2,ID)  accumulated fixed source term; output
*
*     IB,IE     first and last interior indices in i; input
*     ID,JD     array dimensions; input
*
******************************************************************
*
      REAL AUP(2,2,ID),AUW(2,2,ID),AUE(2,2,ID),BU(2,ID)
      REAL DIEP(ID),VOLP(ID)
      INTEGER IB,IE,ID
      INTEGER I
*
      PRINT *, "AUP(2,1,I) AUW(2,1,I) AUE(2,1,I)"
      DO 10 I=IB,IE
        AUP(2,1,I) = 0.5*VOLP(I)*(1/DIEP(I-1)-1/DIEP(I))
        AUW(2,1,I) = 0.5*VOLP(I)/DIEP(I-1)
        AUE(2,1,I) = -0.5*VOLP(I)/DIEP(I)
        PRINT *, I,":",AUP(2,1,I), AUW(2,1,I), AUE(2,1,I)
 10   CONTINUE
*
      RETURN
      END