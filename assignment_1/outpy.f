*
************************************************************************
*
      SUBROUTINE OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T)
*     
*     Print output of 1d to a python formatted file for plotting
*     Formatted into columns by coefficient, rows are node number 
*     Used as input to assignment_1.py to plot temperature distributions
*     Format for astropy.table- Table.read()
*
************************************************************************
*
*     Variable declaration
*
      REAL T(ID),DE(ID),ATW(ID),ATP(ID),BT(ID)
      INTEGER I,IDATOP
*
************************************************************************
*
*     Open python formatted file
*
      IDATOP = 91
      OPEN(IDATOP,FILE="outpy.txt")
*
*     Write necessary output to outpy.txt
*
      WRITE(IDATOP,*)'# I T DE ATW ATP BT'
      DO 10 I=IB-1,IE+1
          WRITE(IDATOP,7000) I, T(I), DE(I), ATW(I), ATP(I), BT(I)
 10   CONTINUE
 7000 FORMAT(' ',I5,' ',F7.3,
     C       ' ',F10.3,' ',F10.3,
     C       ' ',F10.3,' ',F10.3)
      RETURN
      END
