*     This file contains 2 subroutines: COEFCN and MASFLX
*
************************************************************************
*
      SUBROUTINE COEFCN(ACUW,ACUE,BC, RHO,AREP,IB,IE,ID)
*
*     Subroutine to calculate the coefficients of the continuity
*     equation for each control volume.
*
*     ACUW(ID)   coefficient for west face uhat; output
*     ACUE(ID)   coefficient for east face uhat; output
*     BC(ID)     mass source term; output
*
*     RHO        fluid density (kg/m^3); input
*     AREP(ID)   c.v. area of face at e point; input
*     IB,IE      first and last interior indices in i; input
*     ID         array dimension; input     
*
*     Notes: 1) The form of the continuity equation is:
*
*               ACUE*UHE(I)+ ACUW*UHE(I-1) + BC(I) = 0
*
***********************************************************************
*
*     Calculate coefficients for mass equation 
*     Equation of the form: ACUE(I)*UHE(I)+ACUW(I)*UHE(I-1)+BC(I)=0
*
      REAL ACUW(ID),ACUE(ID),BC(ID)
      REAL RHO,AREP(ID)
      INTEGER ID,IB,IE
      INTEGER I
*
      DO 10 I=IB,IE
        ACUE(I) = RHO*AREP(I)
        ACUW(I) = -RHO*AREP(I-1)
        BC(I) = 0
 10   CONTINUE    
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE MASFLX(ME, UHE,ACUW,ACUE,BC,IB,IE,ID)
*
*     Subroutine to calculate the mass flow for the 
*     normal advection fluxes through the east face of
*     each control volume. 
*
*     REAL ME(ID) normal mass flux for east face; output
*
*     UHE(ID)   x component of mass velocity at e point; input
*     ACUW(ID)  coefficient for west face uhat; input
*     ACUE(ID)  coefficient for east face uhat; input
*     BC(ID)    mass source term; input
*
*     IB,IE     first and last interior indices in i; input
*     ID        array dimension; input     
*
*     Notes: 1) ME is calculated for the e points on all control-
*               volume faces.
*
*            2) This subroutine must be preceded by a call to COEFCN.
*
***********************************************************************
*
      REAL ME(ID),UHE(ID),ACUW(ID),ACUE(ID),BC(ID),RHO
      INTEGER IB,IE,ID,I
*
      DO 20 I=IB,IE
        ME(I) = -ACUE(I)*UHE(I)
 20	  CONTINUE

      RETURN
      END
