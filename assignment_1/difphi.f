*
*        file difphi.f
************************************************************************
*
      SUBROUTINE DIFPHI(DE, GAMA,AREP,DIEP,IB,IE,ID)
*
*     Subroutine to calculate the diffusion coefficients for  
*     normal diffusion fluxes through the east face of
*     each control volume.
*
*     REAL DE(ID) diffusion coefficient for east face; output
*
*     REAL GAMA diffusivity * density (kg/s/m); input
*     REAL AREP(ID) c.v. area of face at e point; input
*     REAL DIEP(ID) grid distance in i at e point; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER ID array dimensions; input     
*
*     Notes: 1) DE is calculated for the e points on
*               interior control volume faces.
*
***********************************************************************
*     Initialize Variables
      REAL DE(ID),            ! diffusion coefficient for east face, output
     C     AREP(ID),DIEP(ID), ! CV area of face at e, grid distance in i at e, input
     C     GAMA               ! diffusivity * density, input
      INTEGER IB,IE,ID,       ! first and last interior indices, array dimensions
     C        I               ! integer for loops
***********************************************************************
*     East Diffusion coefficient calculation
*     
      DO 10 I=IB-1,IE                       ! Loop over all internal nodes
          DE(I) = GAMA * AREP(I) / DIEP(I)    ! Calculate De at each node and assign to I location
   10 CONTINUE
*
      RETURN
      END
