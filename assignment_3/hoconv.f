*
*     file hoconv.f
************************************************************************************
*
      SUBROUTINE HOCONV(DCCE, ADVSCM,IB,IE,ID,DE,T,RHO,
     C                  ME,ALFAE,XP,XE)
*
*     Routine to calculate deferred corrections on convection based
*     on the selection of a higher order convection scheme.
*
*     DCCE(ID)  Deferred Correction on Convection at the East face
*
************************************************************************************
*
*     Define Variables
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
        DCCE(IB-1) = 0                    !Set external face
        DO 10 I=IB,IE-1                   !Only internal faces
          THOS = 0.5*(T(I)+T(I+1))        !CDS Temperature
          TE = T(I)*(1+ALFAE(I))/2 + T(I+1)*(1-ALFAE(I))/2
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*TE)
          print *, I, dcce(i)             !Check values
 10	    CONTINUE
        DCCE(IE) = 0                      !Set external face
*
      ELSEIF(ADVSCM == 3) THEN
*
*     QUICK scheme
*
        IF (ME(I) > 0) THEN
         PRINT *, "QUICK"
         DCCE(IB-1) = 0                    !Set external face
         DO 20 I=IB,IE-1                   ! ONly internal faces
          THOS = (XE(I)-XP(I))*(XE(I)-XP(I+1))/
     C             (XP(I-1)-XP(I))/(XP(I-1)-XP(I+1))*T(I-1)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I+1))/
     C             (XP(I)-XP(I-1))/(XP(I)-XP(I+1))*T(I)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I))/
     C             (XP(I+1)-XP(I-1))/(XP(I+1)-XP(I))*T(I+1)
          TE = T(I)*(1+ALFAE(I))/2 + T(I+1)*(1-ALFAE(I))/2
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*TE)
          print *, I, dcce(i)             !Check values
 20	     CONTINUE
         DCCE(IE) = 0                      !Set external face
        ELSE 
         PRINT *, "QUICK"
         DCCE(IB) = 0                    !Set external face
         DO 30 I=IB+1,IE                   ! ONly internal faces
          THOS = (XE(I)-XP(I))*(XE(I)-XP(I+1))/
     C             (XP(I-1)-XP(I))/(XP(I-1)-XP(I+1))*T(I-1)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I+1))/
     C             (XP(I)-XP(I-1))/(XP(I)-XP(I+1))*T(I)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I))/
     C             (XP(I+1)-XP(I-1))/(XP(I+1)-XP(I))*T(I+1)
          TE = T(I)*(1+ALFAE(I))/2 + T(I+1)*(1-ALFAE(I))/2
          DCCE(I) = BETA*(ME(I)*THOS-ME(I)*TE)
          print *, I, dcce(i)             !Check values
 30	     CONTINUE
         DCCE(IE+1) = 0                      !Set external face
        ENDIF 
      ENDIF
*
      RETURN
      END