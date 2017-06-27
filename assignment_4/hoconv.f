*
*     file hoconv.f
************************************************************************************
*
      SUBROUTINE HOCONV(DCCE, ADVSCM,IB,IE,ID,PHI,RHO,
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
      REAL DCCE(ID),DE(ID),PHI(ID),QT(ID),RT(ID),ME(ID),ALFAE(ID)
      REAL PHIHOS,RHO,XP(ID),XE(ID),BETA,PHIE
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
        DO 5 I=IB-1,IE+1
         PRINT *, I,DCCE(I)
 5      CONTINUE   
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
        PRINT *, "DCCE"
        DCCE(IB-1) = 0                    !Set external face
        print *, IB-1, dcce(IB-1)    
        DO 10 I=IB,IE-1                   !Only internal faces
          PHIHOS = 0.5*(PHI(I)+PHI(I+1))        !CDS Temperature
          PHIE = PHI(I)*(1+ALFAE(I))/2 + PHI(I+1)*(1-ALFAE(I))/2
          DCCE(I) = BETA*(ME(I)*PHIHOS-ME(I)*PHIE)
          print *, I, dcce(i)             !Check values
 10	    CONTINUE
        DCCE(IE) = 0                      !Set external face
        print *, IE, dcce(IE)  
*
      ELSEIF(ADVSCM == 3) THEN
*
*     QUICK scheme
*
        PRINT *, "QUICK"
        PRINT *, "DCCE"
        DCCE(IB-1) = 0                    !Set external face
        print *, IB-1, dcce(IB-1)
        DO 20 I=IB,IE-1                   ! Only internal faces
          PHIHOS = (XE(I)-XP(I))*(XE(I)-XP(I+1))/
     C             (XP(I-1)-XP(I))/(XP(I-1)-XP(I+1))*PHI(I-1)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I+1))/
     C             (XP(I)-XP(I-1))/(XP(I)-XP(I+1))*PHI(I)
     C          +(XE(I)-XP(I-1))*(XE(I)-XP(I))/
     C             (XP(I+1)-XP(I-1))/(XP(I+1)-XP(I))*PHI(I+1)
          PHIE = PHI(I)*(1+ALFAE(I))/2 + PHI(I+1)*(1-ALFAE(I))/2
          DCCE(I) = BETA*(ME(I)*PHIHOS-ME(I)*PHIE)
          print *, I, dcce(i)             !Check values
 20	    CONTINUE
        DCCE(IE) = 0                      !Set external face
        print *, IE, dcce(IE)
      ENDIF
*
      RETURN
      END