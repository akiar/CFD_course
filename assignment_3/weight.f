*
*    Note this file contains 2 subroutines, WEIGHT and PRFL
************************************************************************************
*
      SUBROUTINE WEIGHT(ALFAE, ME,DE,IB,IE,ID)
*
*     Subroutine to calculate the advection weighting factors
*     for estimating fluxes through the east faces of
*     each control volume.
*
*     REAL ALFAE(ID) advection factor for east face; output
*
*     REAL ME(ID)    normal mass flux for east face; input
*     REAL DE(ID)    diffusion coefficient for east face; input
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER ID     array dimensions; input
*
*     Notes:
*
*     1) It is best to calculate the weighting factors with 
*        a subroutine PRFL, to facilitate changing to profile used
*        to calculate alfa. This routine will make calls to
*        PRFL for each e point.
*     2) The boundary values of ALFAE are calculated outside of
*        the loop to ensure that the appropriate conditions are
*        carried into the domain.
*
************************************************************************************
*
      REAL ALFAE(ID)
      REAL ME(ID),DE(ID)
      INTEGER IB,IE,ID,IEM1,I
*
      ALFAE(IB-1)= 1      !Make sure physically correct
      DO 1 I=IB,IE-1              ! All internal faces
        CALL PRFL(ALFAE(I),ME(I),DE(I))
  1   CONTINUE
      ALFAE(IE)= -1       !Make sure physically correct
*
      RETURN
      END
*
************************************************************************************
*
      SUBROUTINE PRFL(ALFA,M,D)
*
*     Subroutine to calculate the weighting factors alfa and beta
*     at a point.
*
*     ALFA advection weighting factor; output
*
*     M    mass flux for face at point; input
*     D    diffusion coefficient for face at point; input
*
************************************************************************************
*
      REAL ALFA,M,D
*
*     Check the mass flow rate for alpha assignment
*     
      IF (M > 0) THEN     !Forward flow 
        ALFA = 1
      ELSEIF (M < 0) THEN !Reversed Flow
        ALFA = -1
      ELSE                ! No flow
        ALFA = 0
      ENDIF
*
      RETURN
      END