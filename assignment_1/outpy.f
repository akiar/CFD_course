*
************************************************************************
*
      SUBROUTINE OUTPY(ID,IB,IE,DE,ATW,ATP,BT,T)
*     
*     Print output of 1d to a python formatted file for plotting
*
*
************************************************************************
*      
      REAL T(ID),DE(ID),ATW(ID),ATP(ID),BT(ID)
      INTEGER I,IDATOP
      
      IDATOP = 90
      OPEN(IDATOP,FILE="outpy.txt")
*
      WRITE(IDATOP,*)'# I T DE ATW ATP BT'
      DO 10 I=IB-1,IE+1
          WRITE(IDATOP,7000) I, T(I), DE(I), ATW(I), ATP(I), BT(I) !ATE(I) ,
 10   CONTINUE
 7000 FORMAT(' ',I2,' ',F7.3,
     C       ' ',F10.3,' ',F10.3,
     C       ' ',F10.3,' ',F10.3)
      RETURN
      END
