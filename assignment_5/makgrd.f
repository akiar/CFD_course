*
*           file makgrd.f
********************************************************************
      SUBROUTINE MAKGRD(XNE,YNE, GRDXI,GRDX,GRDYI,GRDY,
     C                  IRATX1,IRATX2,IRATY1,IRATY2,
     C                  IBW,IEW,JBW,JEW,X_IB,X_IE,Y_JB,Y_JE,
     C                  IB,IE1,IE,JB,JE1,JE,ID,JD)
*
*     Subroutine to calculate the northeast corner locations of a
*     simple non-orthogonal grid. The grid is formed from a pin-jointed
*     structure with non-uniform properties.
*
*     XNE(ID,JD)   x location of northeast corners; output
*     YNE(ID,JD)   y location of northeast corners; output
*
*     GRDX,GRDY     grid lengths in 'x' and 'y' directions; input
*     GRDXI,GRDYI   position of interface in 'x' and 'y'; input
*     IRATX1,IRATY1 ratio of 1st to last c.v. in 'x' and 'y'; input
*     IRATX2,IRATY2 ratio of last to 1st c.v. in 'x' and 'y'; input
*     INTEGER IB    index of first interior volume in i direction;input
*     INTEGER IE    index of last interior volume in i direction; input
*     INTEGER JB    index of first interior volume in j direction; input
*     INTEGER JE    index of last interior volume in j direction; input
*     INTEGER ID,JD dimensions of arrays
*
*     Notes: 1) Corner locations must be calculated for volumes with 
*               indices running from IB-1 to IE and JB-1 to JE
*
********************************************************************
*
      REAL XNE(ID),YNE(JD)
      REAL GRDXI,GRDX,GRDYI,GRDY
      REAL IRATX1,IRATX2,IRATY1,IRATY2
      REAL SUMDX,SUMDY,DX1,DY1,GRADEX,GRADEY
      REAL X_IB,X_IE,Y_JB,Y_JE
      INTEGER IBW,IEW,JBW,JEW
      INTEGER IB,IE1,IE,JB,JE1,JE,ID,JD
      INTEGER NCVX,NCVY,I,J
*
*  Check for input errors
*
      IF(IE1.GT.IE.OR.JE1.GT.JE) GOTO 111
*
      NCVX= IE-IB+1
      NCVY= JE-JB+1
*
*  Set position of lower left corner
*
      XNE(IB-1)= 0.0
      YNE(JB-1)= 0.0
*
*  Construct 'x' direction
*
*  (a) from IB to IE1
*
      IF(IBW.EQ.1) THEN
        XNE(IB)= X_IB
      ENDIF
      NCVX= IE1-IB+1-IBW
      IF(NCVX.GT.1) THEN
        SUMDX= 0.0
        DO 1 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX1)/(NCVX-1))
  1     CONTINUE
        DX1= (GRDXI-X_IB)/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX1)
        XNE(IB+IBW)= XNE(IB+IBW-1)+DX1
        DO 2 I=2,NCVX
          XNE(I+IBW+1)= XNE(I+IBW)+(DX1-(I-1)*GRADEX)
  2     CONTINUE
      ELSE IF(NCVX.EQ.1) THEN
        XNE(IB)= GRDX
      ELSE
        GOTO 111
      ENDIF
*
*  (b) from IE1 to IE
*
      NCVX= IE-IE1-IEW
      IF(NCVX.LT.1) GOTO 3
        SUMDX= 0.0
        DO 4 I=1,NCVX
          SUMDX= SUMDX+(1.0-(I-1)*(1.0-1.0/IRATX2)/(NCVX-1))
  4     CONTINUE
        DX1= (GRDX-GRDXI-X_IB)/SUMDX
        GRADEX= DX1/(NCVX-1)*(1.0-1.0/IRATX2)
        XNE(IE1+1)= XNE(IE1)+DX1
        DO 5 I=2,NCVX
          XNE(IE1+I)= XNE(IE1+I-1)+(DX1-(I-1)*GRADEX)
  5     CONTINUE
        IF(IEW.EQ.1) THEN
          XNE(IE)= GRDX
        ENDIF
  3     XNE(IE+1)= XNE(IE)
*
*  Construct 'y' direction
*
*  (a) from JB to JE1
*
      IF(JBW.EQ.1) THEN
        YNE(JB)= Y_JB
      ENDIF
      NCVY= JE1-JB+1-JBW
      IF(NCVY.GT.1) THEN
        SUMDY= 0.0
        DO 6 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY1)/(NCVY-1))
  6     CONTINUE
        DY1= (GRDYI-Y_JB)/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY1)
        YNE(JB+JBW)= YNE(JB+JBW-1)+DY1
        DO 7 J=2,NCVY
          YNE(J+JBW+1)= YNE(J+JBW)+(DY1-(J-1)*GRADEY)
  7     CONTINUE
      ELSE IF(NCVY.EQ.1) THEN
        YNE(JB)= GRDY
      ELSE
        GOTO 111
      ENDIF
*
*  (b) from JE1 to JE
*
      NCVY= JE-JE1-JEW
      IF(NCVY.LT.1) GOTO 8
        SUMDY= 0.0
        DO 9 J=1,NCVY
          SUMDY= SUMDY+(1.0-(J-1)*(1.0-1.0/IRATY2)/(NCVY-1))
  9     CONTINUE
        DY1= (GRDY-GRDYI-Y_JE)/SUMDY
        GRADEY= DY1/(NCVY-1)*(1.0-1.0/IRATY2)
        YNE(JE1+1)= YNE(JE1)+DY1
        DO 10 J=2,NCVY
          YNE(JE1+J)= YNE(JE1+J-1)+(DY1-(J-1)*GRADEY)
 10     CONTINUE
        IF(JEW.EQ.1) THEN
          YNE(JE)= GRDY
        ENDIF
  8     YNE(JE+1)= YNE(JE)
*
      RETURN
*
 111  WRITE(*,*) 'ERROR IN makgrd.f'
      STOP
      END

