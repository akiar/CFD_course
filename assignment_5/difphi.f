*
*          file difphi.f
*
************************************************************************
*
      SUBROUTINE DIFPHI(DE,DN, GAMA,AREP,ARNP,DIEP,DJNP,
     C                  IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the diffusion coefficients for  
*     normal conduction fluxes through the east and north faces of
*     each control volume. These coefficients are used to calculate
*     the finite volume equations and to assess the relative role
*     of diffusion and advection.
*
*     REAL DE(ID,JD) diffusion coefficient for east face; output
*     REAL DN(ID,JD) diffusion coefficient for north face; output
*
*     REAL GAMA diffusivity * density (kg/s/m); input
*     REAL MUTE(ID,JD) turbulent viscosity at e faces; input
*     REAL MUTN(ID,JD) turbulent viscosity at n faces; input
*     REAL AREP(JD) c.v. area of face at e point; input
*     REAL ARNP(ID) c.v. area of face at n point; input
*     REAL DIEP(ID) distance from P to E through e; input
*     REAL DJNP(JD) distance from P to N through n; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input     
*
*     Notes: 1) DE and DN are calculated for the e and n points on
*            interior control volume and on the boundaries as required
*            to cover all interior control volume faces.
*
***********************************************************************
*
      REAL DE(ID,JD),DN(ID,JD)
      REAL AREP(JD),ARNP(ID)
      REAL DIEP(ID),DJNP(JD)
      REAL GAMA
      INTEGER IB,IE,JB,JE,I,J
*
*  Bottom 'n' faces
*
      DO 10 I=IB,IE
        J=JB-1
        DN(I,J)=GAMA*ARNP(I)/DJNP(J)
 10   CONTINUE
*
*  Interior 'e' and 'n' faces
*
      DO 30 J=JB,JE
        I=IB-1
        DE(I,J)=GAMA*AREP(J)/DIEP(I)
        DO 20 I=IB,IE
          DE(I,J)=GAMA*AREP(J)/DIEP(I)
          DN(I,J)=GAMA*ARNP(I)/DJNP(J)
 20     CONTINUE
 30   CONTINUE
*
      RETURN
      END
