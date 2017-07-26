*
*     Note this file contains 2 subroutines, WEIGHT and PROFILE
*
*
************************************************************************
*
      SUBROUTINE WEIGHT(ALFAE,ALFAN, ME,MN,DE,DN,
     C                IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the advection weighting factors
*     for estimating fluxes through the east and north faces of
*     each control volume.
*
*     REAL ALFAE(ID,JD) advection factor for east face; output
*     REAL ALFAN(ID,JD) advection factor for north face; output
*
*     REAL ME(ID,JD) normal mass flux for east face; input
*     REAL MN(ID,JD) normal mass flux for north face; input
*     REAL DE(ID,JD) diffusion coefficient for east face; input
*     REAL DN(ID,JD) diffusion coefficient for north face; input
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input     
*
*     Notes: 1) It is best to calculate the weighting factors with 
*            a subroutine PRFL. This routine will make calls to
*            PRFL for each e and n point.
*
***********************************************************************
*
      REAL ALFAE(ID,JD),ALFAN(ID,JD)
      REAL ME(ID,JD),MN(ID,JD),DE(ID,JD),DN(ID,JD)
      INTEGER IB,IE,JB,JE,ID,JD
      INTEGER I,J
*
*  Interior north faces
*
      DO 10 I=IB,IE
         ALFAN(I,JB-1)=1.0
         DO 20 J=JB,JE-1
            CALL PRFL(ALFAN(I,J),MN(I,J),DN(I,J))
 20      CONTINUE
         ALFAN(I,JE)=-1.0
 10   CONTINUE
*
*  Interior east faces
*
      DO 30 J=JB,JE
         ALFAE(IB-1,J)=1.0
         DO 40 I=IB,IE-1
            CALL PRFL(ALFAE(I,J),ME(I,J),DE(I,J))      
 40      CONTINUE
         ALFAE(IE,J)=-1.0
 30   CONTINUE
*
      RETURN
      END
*
*     
***********************************************************************
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
*     Notes: 1) The weighting factor must be 1.0 or -1.0 for pure UDS
*
************************************************************************
*
      REAL ALFA,M,D,PE
*
      IF(M.GE.0.0) THEN
        ALFA= 1.0
      ELSE
        ALFA=-1.0
      ENDIF
*
      RETURN
      END
