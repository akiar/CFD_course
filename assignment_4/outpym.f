*
*     outpy.f
************************************************************************
*
      SUBROUTINE OUTPYM(ID,IB,IE,DE,ATW,ATP,BT,P,U,XP,KNTOUT,TITLE)
*     
*     Print output of 1d analysis to a python formatted file for plotting
*     Formatted into columns by coefficient, rows are node number 
*     Used as input to assignment_1.py to plot temperature distributions
*     Format for astropy.table- Table.read()
*
************************************************************************
*
*     Variable declaration
*
      REAL P(ID),U(ID),DE(ID),ATW(ID),ATP(ID),BT(ID),XP(ID)
      CHARACTER*1 TITLE
      INTEGER I,IDATOP
*
************************************************************************
*
*     Open python formatted file
*
      IDATOP = 92
      OPEN(IDATOP,FILE="outpym.txt")
*
*     Write necessary output to outpy.txt
*
      IF (KNTOUT==0) THEN
          WRITE(IDATOP,*)'# I XP P U DE ATW ATP BT' ! Print headeR    
      ELSE
          DO 20 I=IB-1,IE+1                       ! Loop over all control volumes
              WRITE(IDATOP,7010) I,XP(I),P(I),U(I),
     C                           DE(I),ATW(I),ATP(I),BT(I)   !print information
 20       CONTINUE
 7010     FORMAT(' ',I5,' ',F12.4,' ',F12.4,      ! Set format
     C           ' ',F12.4,' ',F12.4,' ',F12.4,
     C           ' ',F12.4,' ',F12.4)
      ENDIF
      RETURN
      END
