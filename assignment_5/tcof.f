*        file tcof.f
************************************************************************
*
      SUBROUTINE TCOF(DE,DN,ALFAE,ALFAN,QT,RT,
     C            ATP,ATW,ATE,ATS,ATN,BT,
     C            LVLCOF,IDATO,IB,IE,JB,JE,ID,JD)
*
*  Routine to print out coefficients for the energy equation.
*
************************************************************************
*
      REAL DE(ID,JD),DN(ID,JD),QT(ID,JD),RT(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD)
      REAL ATW(ID,JD),ATE(ID,JD),ATS(ID,JD)
      REAL ATN(ID,JD),ATP(ID,JD),BT(ID,JD)
*
      INTEGER LVLCOF,IDATO,IB,IE,JB,JE,ID,JD
*
*  Produce output for LVLCOF>0 
*
      IF(LVLCOF.GE.1) THEN
	  WRITE(IDATO,6001)
 6001 FORMAT('      T Equation Transport Coefficients  ')
      CALL OUT2D(DE,   'DE      ',IDATO,IB-1,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(DN,   'DN      ',IDATO,IB  ,IE  ,1,JB-1,JE  ,1,ID,JD)
	  CALL OUT2D(ALFAE,'ALFAE   ',IDATO,IB-1,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(ALFAN,'ALFAN   ',IDATO,IB  ,IE  ,1,JB-1,JE  ,1,ID,JD)
	  CALL OUT2D(QT,   'QT      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(RT,   'RT      ',IDATO,IB  ,IE  ,1,JB  ,JE  ,1,ID,JD)
	  CALL OUT2D(ATP,  'ATP     ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(ATW,  'ATW     ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(ATE,  'ATE     ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(ATS,  'ATS     ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(ATN,  'ATN     ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
	  CALL OUT2D(BT ,  'BT      ',IDATO,IB-1,IE+1,1,JB-1,JE+1,1,ID,JD)
      ENDIF
*
      RETURN
      END
