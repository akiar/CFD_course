*
*     file hoconv.f
**********************************************************************
*
      SUBROUTINE HOCONV(DCCE, ADVSCM,IB,IE,ID,DE,T,QT,RT,RHO,
     C                  ME,ALFAE,XP,XE)
*
*     Routine to calculate deferred corrections on convection based
*     on the selection of a higher order convection scheme.
*
*     DCCE(ID)  Deferred Correction on Convection at the East face
*
**********************************************************************
*
      REAL DCCE(ID),DE(ID),T(ID),QT(ID),RT(ID),ME(ID),ALFAE(ID)
      REAL THOS,RHO,XP(ID),XE(ID),BETA,TE
      INTEGER ADVSCM,IB,IE,ID
*
*     -- SET BLENDING FACTOR
*
      BETA = 1.0
*
*----------------------------
* UDS Scheme
*----------------------------
*
      IF(ADVSCM == 1) THEN
*
        CALL NULL(DCCE, IB,IE,ID)
        DCCE(IB-1) = 0    ! Set beginning and end conditions
        DCCE(IE+1) = 0
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
        DCCE(IB-1) = 0
        DO 10 I=IB,IE
          TE = (1+ALFAE(I))*0.5*T(I) + (1-ALFAE(I))*0.5*T(I+1)
          THOS = 0.5*(T(I)+T(I+1))
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*TE)
 10	    CONTINUE
        DCCE(IE+1) = 0
*
      ELSEIF(ADVSCM == 3) THEN
*
*     QUICK scheme
*
        PRINT *, "QUICK"
        DCCE(IB-1) = 0
        DO 20 I=IB,IE
*          THOS = (XE(I)-XP(I))*(XE(I)-XP(I+1))/
*     C             (XP(I-1)-XP(I))/(XP(I-1)-XP(I+1))*T(I-1)
*     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I+1))/
*     C             (XP(I)-XP(I-1))/(XP(I)-XP(I+1))*T(I)
*     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I))/
*     C             (XP(I)-XP(I-1))/(XP(I+1)-XP(I))*T(I+1)
          THOS = -1.0/8.0*T(I-1) + 3.0/4.0*T(I) + 3.0/8.0*T(I+1)
          TE = (1+ALFAE(I))*0.5*T(I) + (1-ALFAE(I))*0.5*T(I+1)
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*TE)
 20	    CONTINUE
        DCCE(IE+1) = 0
      ENDIF
*
      RETURN
      END