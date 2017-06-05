*
*     file hoconv.f
**********************************************************************
*
      SUBROUTINE HOCONV(DCCE, ADVSCM,IB,IE,ID,DE,T,QT,RT,RHO,VOLP,DIEP
     C                  ME,ALFAE,UHE)
*
*     Routine to calculate deferred corrections on convection based
*     on the selection of a higher order convection scheme.
*
*     DCCE(ID)  Deferred Correction on Convection at the East face
*
**********************************************************************
*
      REAL DCCE(ID),DE(ID),T(ID),QT(ID),RT(ID),ME(ID),ALFAE(ID),UHE(ID)
      REAL THOS,PAC,PACP,PACE,RHO,VOLP(ID),DIEP(ID),BETA
      INTEGER ADVSCM,IB,IE,ID
*
*     -- SET ADVECTION SCHEME
*
      ADVSCM = 3
      BETA = 1
*
*----------------------------
* UDS Scheme
*----------------------------
*
      IF(ADVSCM == 1) THEN
*
        CALL NULL(DCCE, IB,IE,ID)
        PRINT *, "UDS"
*
*----------------------------
* CDS Scheme and QUICK
* UDS = UDS - HOS
*----------------------------
*
      ELSEIF(ADVSCM == 2) THEN
*
*     PAC Scheme - reduces to CDS
*
        PRINT *, "CDS"
        DO 10 I=IB,IE
*          PACP = (DE(I)*(T(I+1)-T(I))-DE(I-1)*(T(I)-T(I-1))
*     C              +QT(I)+RT(I)*T(I))/(RHO*VOLP(I))
*          PACE = (DE(I+1)*(T(I+2)-T(I+1))-DE(I)*(T(I+1)-T(I))
*     C               +QT(I+1)+RT(I+1)*T(I+1))/(RHO*VOLP(I))
*          PAC = 0.5 * (PACP + PACE)
          THOS = 0.5*(T(I)+T(I+1))
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*T(I))
          PRINT *, DCCE(I)
*
 10	    CONTINUE
*
      ELSEIF(ADVSCM == 3) THEN
*
*     QUICK scheme
*
        PRINT *, "QUICK"
        DO 20 I=IB,IE
*          THOS=DIEP(I)*DIEP(I)/4/(DIEP(I-1)*(DIEP(I)+DIEP(I-1)))*T(I-1)
*     C         +(DIEP(I-1)+DIEP(I)/2)*DIEP(I)/2/(DIEP(I-1)*DIEP(I))*T(I)
*     C         +(DIEP(I-1)+DIEP(I)/2)*DIEP(I)/2/((DIEP(I)+DIEP(I-1))
*     C                                                *DIEP(I))*T(I+1)
          THOS = -1.0/8.0*T(I-1) + 3.0/4.0*T(I) + 3.0/8.0*T(I+1)
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*T(I))
 20	    CONTINUE
      ENDIF
*
      RETURN
      END