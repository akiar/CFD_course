*
*        file hoconv.f
************************************************************************
*
      SUBROUTINE HOCONV(DCCE,DCCN, ALFAE,ALFAN,ME,MN,DE,DN,
     C          XP,XNE,YP,YNE,DISE,DIEP,DISN,DJNP,RHO,VOLP,
     C          U,UHE,V,VHN,PHI,
     C          ADVSCM,BLEND,IB,IE,JB,JE,ID,JD,
     C          XE,YN)
*
*     Subroutine to calculate the deferred correction when using the 
*     CDS or QUICK convection schemes.
*
*     Variables:
*
*     DCCE,DCCN     deferred corrections for east and north faces
*     ALFAE,ALFAN   convection weights for east and north faces
*     ME,MN         mass fluxes for east and north faces
*     PHI           variable to which scheme is applied
*     ADVSCM        scheme selected
*     BLEND         blending function (set to 1 unless you have a 
*                   good reason not to)
*
***********************************************************************
*
      REAL DCCE(ID,JD),DCCN(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD),ME(ID,JD),MN(ID,JD)
      REAL DE(ID,JD),DN(ID,JD),XP(ID),XNE(ID),YP(JD),YNE(JD)
      REAL DISE(ID),DIEP(ID),DISN(JD),DJNP(JD),RHO,VOLP(ID,JD)
      REAL U(ID,JD),UHE(ID,JD),V(ID,JD),VHN(ID,JD)
      REAL PHI(ID,JD),BLEND
      INTEGER ADVSCM,IB,IE,JB,JE,ID,JD
*
      REAL PHIEHOS,PHIE,PHINHOS,PHIN,YN(JD),XE(ID)
*
*     UDS
*
      IF (ADVSCM == 1) THEN
       CALL NULLM(DCCE, IB,IE,JB,JE,ID,JD)
       CALL NULLM(DCCN, IB,IE,JB,JE,ID,JD)
*
*     QUICK
*
      ELSEIF (ADVSCM == 2) THEN
*
*      INTERIOR EAST FACES
*       
       DO 10 J=JB,JE
        DCCE(IB-1,J) = 0.0
        DO 20 I=IB,IE-1
         IF (ME(I,J) > 0) THEN 
          PHIEHOS = (XE(I)-XP(I))*(XE(I)-XP(I+1))/
     C              (XP(I-1)-XP(I))/(XP(I-1)-XP(I+1))*PHI(I-1,J)
     C             +(XE(I)-XP(I-1))*(XE(I)-XP(I+1))/
     C              (XP(I)-XP(I-1))/(XP(I)-XP(I+1))*PHI(I,J)
     C             +(XE(I)-XP(I-1))*(XE(I)-XP(I))/
     C              (XP(I+1)-XP(I-1))/(XP(I+1)-XP(I))*PHI(I+1,J)
          PHIE = PHI(I,J)*(1+ALFAE(I,J))/2+PHI(I+1,J)*(1-ALFAE(I,J))/2
          DCCE(I,J) = BLEND*ME(I,J)*(PHIEHOS-PHIE)
         ELSE 
          PHIEHOS = (XE(I)-XP(I+1))*(XE(I)-XP(I+2))/
     C              (XP(I)-XP(I+1))/(XP(I)-XP(I+2))*PHI(I,J)
     C             +(XE(I)-XP(I))*(XE(I)-XP(I+2))/
     C              (XP(I+1)-XP(I))/(XP(I+1)-XP(I+2))*PHI(I+1,J)
     C             +(XE(I)-XP(I))*(XE(I)-XP(I+1))/
     C              (XP(I+2)-XP(I))/(XP(I+2)-XP(I+1))*PHI(I+2,J)
          PHIE = PHI(I,J)*(1+ALFAE(I,J))/2+PHI(I+1,J)*(1-ALFAE(I,J))/2
          DCCE(I,J) = BLEND*ME(I,J)*(PHIEHOS-PHIE)
         ENDIF
   20   CONTINUE
        DCCE(IE,J)=0.0
   10  CONTINUE
*
*      INTERIOR NORTH FACES
*       
       DO 30 I=IB,IE
        DCCN(I,JB-1) = 0.0
        DO 40 J=JB,JE-1
         IF (MN(I,J) > 0) THEN 
          PHINHOS = (YN(J)-YP(J))*(YN(J)-YP(J+1))/
     C              (YP(J-1)-YP(J))/(YP(J-1)-YP(J+1))*PHI(I,J-1)
     C             +(YN(J)-YP(J-1))*(YN(J)-YP(J+1))/
     C              (YP(J)-YP(J-1))/(YP(J)-YP(J+1))*PHI(I,J)
     C             +(YN(J)-YP(J-1))*(YN(J)-YP(J))/
     C              (YP(J+1)-YP(J-1))/(YP(J+1)-YP(J))*PHI(I,J+1)
          PHIN = PHI(I,J)*(1+ALFAN(I,J))/2+PHI(I,J+1)*(1-ALFAN(I,J))/2
          DCCN(I,J) = BLEND*MN(I,J)*(PHINHOS-PHIN)
         ELSE 
          PHINHOS = (YN(J)-YP(J+1))*(YN(J)-YP(J+2))/
     C              (YP(J)-YP(J+1))/(YP(J)-YP(J+2))*PHI(I,J)
     C             +(YN(J)-YP(J))*(YN(J)-YP(J+2))/
     C              (YP(J+1)-YP(J))/(YP(J+1)-YP(J+2))*PHI(I,J+1)
     C             +(YN(J)-YP(J))*(YN(J)-YP(J+1))/
     C              (YP(J+2)-YP(J))/(YP(J+2)-YP(J+1))*PHI(I,J+2)
          PHIN = PHI(I,J)*(1+ALFAN(I,J))/2+PHI(I,J+1)*(1-ALFAN(I,J))/2
          DCCN(I,J) = BLEND*MN(I,J)*(PHINHOS-PHIN)
         ENDIF
   40   CONTINUE
        DCCN(I,JE)=0.0
   30  CONTINUE
*
*     CDS
*   
      ELSE 
*
*      INTERIOR EAST FACES
*       
       DO 50 J=JB,JE
        DCCE(IB-1,J) = 0.0
        DO 60 I=IB,IE-1
         PHIEHOS = 0.5*(PHI(I,J)+PHI(I+1,J))
         PHIE = PHI(I,J)*(1+ALFAE(I,J))/2+PHI(I+1,J)*(1-ALFAE(I,J))/2
         DCCE(I,J) = BLEND*ME(I,J)*(PHIEHOS-PHIE)
   60   CONTINUE
        DCCE(IE,J) = 0.0
   50  CONTINUE
*
*      INTERIOR NORTH FACES
*
       DO 70 I=IB,IE
        DCCN(I,JB-1) = 0.0
        DO 80 J=JB,JE-1
         PHINHOS = 0.5*(PHI(I,J)+PHI(I,J+1))
         PHIN = PHI(I,J)*(1+ALFAN(I,J))/2+PHI(I,J+1)*(1-ALFAN(I,J))/2
         DCCN(I,J) = BLEND*MN(I,J)*(PHINHOS-PHIN)
   80	CONTINUE
        DCCN(I,JE) = 0.0
   70  CONTINUE
      ENDIF
*
      RETURN
      END
