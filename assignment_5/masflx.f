*     This file contains 2 subroutines: COEFCN and MASFLX
*
************************************************************************
*
      SUBROUTINE COEFCN(ACUW,ACUE,ACVS,ACVN,BC,
     C                  RHO,AREP,ARNP,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the coefficients of the continuity
*     equation for each control volume.
*
*     REAL ACUW(ID,JD) coefficient for west face uhat; output
*     REAL ACUE(ID,JD) coefficient for east face uhat; output
*     REAL ACVS(ID,JD) coefficient for south face vhat; output
*     REAL ACVN(ID,JD) coefficient for north face vhat; output
*     REAL BC(ID,JD)   mass source term; output
*
*     REAL RHO fluid density (kg/m^3); input
*     REAL AREP(JD) c.v. area of face at e point; input
*     REAL ARNP(ID) c.v. area of face at n point; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input     
*
*     Notes: 1) The form of the continuity equation is:
*     ACUE*UHE(I,J) + ACUW*UHE(I-1,J)
*                   + ACVN*VHN(I,J) + ACVS*VHN(I,J-1) + BC = 0
*
***********************************************************************
      REAL ACUW(ID,JD),ACUE(ID,JD),ACVS(ID,JD),ACVN(ID,JD)
      REAL BC(ID,JD),AREP(JD),ARNP(ID)
      REAL RHO
      INTEGER IB,IE,JB,JE,ID,JD
      INTEGER I,J
*     
      DO 20 J=JB,JE
      DO 10 I=IB,IE
        ACUE(I,J)=0.0+RHO*AREP(J)
        ACUW(I,J)=0.0-RHO*AREP(J)
        ACVN(I,J)=0.0+RHO*ARNP(I)
        ACVS(I,J)=0.0-RHO*ARNP(I)
        BC(I,J)=0.0
 10   CONTINUE
 20   CONTINUE
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE MASFLX(ME,MN, UHE,VHN,ACUW,ACUE,ACVS,ACVN,
     C                  IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the mass flow for the 
*     normal advection fluxes through the east and north faces of
*     each control volume. 
*
*     ME(ID,JD) normal mass flux for east face; output
*     MN(ID,JD) normal mass flux for north face; output
*
*     UHE(ID,JD) x component of mass velocity at e point; input
*     VHN(ID,JD) y component of mass velocity at n point; input
*     ACUW(ID,JD) coefficient for west face uhat; input
*     ACUE(ID,JD) coefficient for east face uhat; input
*     ACVS(ID,JD) coefficient for south face vhat; input
*     ACVN(ID,JD) coefficient for north face vhat; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input     
*
*     Notes: 1) ME and MN are calculated for the e and n points on
*            interior control volume and on the boundaries as required
*            to cover all interior control volume faces.
*
*            2) This subroutine must be preceded by a call to COEFCN.
*
***********************************************************************
      REAL ME(ID,JD),MN(ID,JD)
      REAL UHE(ID,JD),VHN(ID,JD)
      REAL ACUW(ID,JD),ACUE(ID,JD)
      REAL ACVS(ID,JD),ACVN(ID,JD)
      INTEGER IB,IE,JB,JE,I,J,ID,JD
*     
*  Bottom 'n' faces
*
      DO 10 I=IB,IE
        MN(I,JB-1)=0.0-ACVS(I,JB)*VHN(I,JB-1)
 10   CONTINUE
*
*  Interior 'e' and 'n' faces
      DO 30 J=JB,JE
        ME(IB-1,J)=0.0-ACUW(IB,J)*UHE(IB-1,J)
        DO 20 I=IB,IE
          ME(I,J)= ACUE(I,J)*UHE(I,J)
          MN(I,J)= ACVN(I,J)*VHN(I,J)
 20     CONTINUE
 30   CONTINUE
*
      RETURN
      END

