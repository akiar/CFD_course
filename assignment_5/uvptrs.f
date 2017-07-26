*
*        file uvptrs.f
************************************************************************
*
      SUBROUTINE UVPTRS(T,P,U,V,
     C            NORM,RSD,ARSD,RSDMAX,RSDOLD,I_MAX,J_MAX,
     C            AUP,AUW,AUE,AUS,AUN,BU,
     C            ATP,ATW,ATE,ATS,ATN,BT,
     C            IDATO,ITERMO,ISOTHM,IB,IE,JB,JE,N,ID,JD)
*
*  Routine to compute the residuals for the T,P,U,V equations
*
************************************************************************
*
      REAL T(ID,JD),P(ID,JD),U(ID,JD),V(ID,JD)
      REAL AUW(N,N,ID,JD),AUE(N,N,ID,JD),AUS(N,N,ID,JD)
      REAL AUN(N,N,ID,JD),AUP(N,N,ID,JD),BU(N,ID,JD)
      REAL ATP(ID,JD),ATW(ID,JD),ATE(ID,JD)
      REAL ATS(ID,JD),ATN(ID,JD),BT(ID,JD)
      REAL NORM(4),RSD(ID,JD),ARSD(4)
      REAL RSDMAX(4),RSDOLD(4),RATE
      INTEGER IDATO,ITERMO,ISOTHM,I_MAX,J_MAX,IB,IE,JB,JE,N,ID,JD
*
*  Find magnitude of average value in each solution field
*
      IF(ISOTHM.GT.0) CALL MEANA(NORM(1), T,IB,IE,JB,JE,ID,JD)
      NORM(2)= 1.0
      CALL MEANA(NORM(3), U,IB,IE,JB,JE,ID,JD)
      CALL MEANA(NORM(4), V,IB,IE,JB,JE,ID,JD)
*
*  Calculate and print residuals for each field
*
      IF(ISOTHM.GT.0) THEN
      CALL RESID(RSD,ARSD(1),RSDMAX(1),I_MAX,J_MAX,
     C           T,ATP,ATW,ATE,ATS,ATN,BT,
     C           IB,IE,JB,JE,ID,JD)
      RATE= (ARSD(1)/NORM(1))/RSDOLD(1)
      WRITE(IDATO,7221) 'T-energy ',RATE,ARSD(1)/NORM(1),
     C                  RSDMAX(1)/NORM(1),I_MAX,J_MAX
      WRITE(ITERMO,7221) 'T-energy ',RATE,ARSD(1)/NORM(1),
     C                  RSDMAX(1)/NORM(1),I_MAX,J_MAX
      ENDIF
*
      CALL RESIDM(RSD,ARSD(2),RSDMAX(2),I_MAX,J_MAX,
     C            P,U,V,AUP,AUW,AUE,AUS,AUN,BU,
     C            1,IB,IE,JB,JE,N,ID,JD)
      RATE= (ARSD(2)/NORM(2))/RSDOLD(2)
      WRITE(IDATO,7221) 'P-mass   ',RATE,ARSD(2)/NORM(2),
     C                  RSDMAX(2)/NORM(2),I_MAX,J_MAX
      WRITE(ITERMO,7221) 'P-mass   ',RATE,ARSD(2)/NORM(2),
     C                  RSDMAX(2)/NORM(2),I_MAX,J_MAX
*
      CALL RESIDM(RSD,ARSD(3),RSDMAX(3),I_MAX,J_MAX,
     C            P,U,V,AUP,AUW,AUE,AUS,AUN,BU,
     C            2,IB,IE,JB,JE,N,ID,JD)
      RATE= (ARSD(3)/NORM(3))/RSDOLD(3)
      WRITE(IDATO,7221) 'U-mom    ',RATE,ARSD(3)/NORM(3),
     C                  RSDMAX(3)/NORM(3),I_MAX,J_MAX
      WRITE(ITERMO,7221) 'U-mom    ',RATE,ARSD(3)/NORM(3),
     C                  RSDMAX(3)/NORM(3),I_MAX,J_MAX
*
      CALL RESIDM(RSD,ARSD(4),RSDMAX(4),I_MAX,J_MAX,
     C            P,U,V,AUP,AUW,AUE,AUS,AUN,BU,
     C            3,IB,IE,JB,JE,N,ID,JD)
      RATE= (ARSD(4)/NORM(4))/RSDOLD(4)
      WRITE(IDATO,7221) 'V-mom    ',RATE,ARSD(4)/NORM(4),
     C                  RSDMAX(4)/NORM(4),I_MAX,J_MAX
      WRITE(ITERMO,7221) 'V-mom    ',RATE,ARSD(4)/NORM(4),
     C                  RSDMAX(4)/NORM(4),I_MAX,J_MAX
      WRITE(IDATO,7222)
      WRITE(ITERMO,7222)
*
 7221 FORMAT(
     C'|',A9,' |',F6.2,' |'1PE9.1'|'1PE9.1'|  (',I3,',',I3,')   |')
 7222 FORMAT(
     C'+-----------------------------------------------------+',/)
*
      RETURN
      END
