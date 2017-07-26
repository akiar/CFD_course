*        file uvpcof.f
************************************************************************
*
      SUBROUTINE UVPCOF(DE,DN,ME,MN,ALFAE,ALFAN,QU,RU,QV,RV,
     C            DHUE,DHVN,
     C            DPDX,DPDY,AUP,AUW,AUE,AUS,AUN,BU,
     C            LVLCOF,IDATO,IB,IE,JB,JE,N,ID,JD)
*
*  Routine to print out coefficients for the mass-momentum equations.
*
************************************************************************
*
      REAL DE(ID,JD),DN(ID,JD),ME(ID,JD),MN(ID,JD)
      REAL QU(ID,JD),RU(ID,JD),QV(ID,JD),RV(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD)
      REAL DHUE(ID,JD),DHVN(ID,JD)
      REAL DPDX(ID,JD),DPDY(ID,JD)
      REAL AUW(N,N,ID,JD),AUE(N,N,ID,JD),AUS(N,N,ID,JD)
      REAL AUN(N,N,ID,JD),AUP(N,N,ID,JD),BU(N,ID,JD)
*
      INTEGER LVLCOF,IDATO,IB,IE,JB,JE,N,ID,JD
*
*  Produce output for LVLCOF = 4
*
      IF(LVLCOF.GE.4) THEN
	  WRITE(IDATO,6001)
 6001 FORMAT('      U,V,P Equation Transport Coefficients  ')
	  CALL OUT2D(DE,   'DE      ',IDATO,IB-1,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(DN,   'DN      ',IDATO,IB  ,IE  ,1,JB-1,JE  ,1,ID,JD)
	  CALL OUT2D(ME,   'ME      ',IDATO,IB-1,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(MN,   'MN      ',IDATO,IB  ,IE  ,1,JB-1,JE  ,1,ID,JD)
	  CALL OUT2D(QU,   'QU      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(RU,   'RU      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(QV,   'QV      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(RV,   'RV      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(ALFAE,'ALFAE   ',IDATO,IB-1,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(ALFAN,'ALFAN   ',IDATO,IB  ,IE  ,1,JB-1,JE  ,1,ID,JD)
	  CALL OUT2D(DPDX ,'DPDX    ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(DPDY ,'DPDY    ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(DHUE, 'DHUE    ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(DHVN, 'DHVN    ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      ENDIF
      IF(LVLCOF.GE.3) THEN
	   WRITE(IDATO,6010)
 6010  FORMAT('      U Equation Active Coefficients  ')
	   CALL WRT2D(AUP,AUW,AUE,AUS,AUN,BU,IDATO,
     C             2,2,IB,IE,JB,JE,N,ID,JD)
        WRITE(IDATO,6011)
 6011   FORMAT('      V Equation Active Coefficients  ')
        CALL WRT2D(AUP,AUW,AUE,AUS,AUN,BU,IDATO,
     C             3,3,IB,IE,JB,JE,N,ID,JD)
      ENDIF
*
      IF(LVLCOF.GE.3) THEN
	   WRITE(IDATO,6012)
 6012  FORMAT('      P Equation Active Coefficients  ')
	   CALL WRT2D(AUP,AUW,AUE,AUS,AUN,BU,IDATO,
     C             1,1,IB,IE,JB,JE,N,ID,JD)
      ENDIF
*
      RETURN
      END
