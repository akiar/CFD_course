*           file cargrd.f
***********************************************************************
*
      SUBROUTINE CARGRD(XNE,YNE, GRDXI,GRDX,GRDYI,GRDY,
     C            IRATX1,IRATX2,IRATY1,IRATY2,
     C            IBW,IEW,JBW,JEW,X_IB,X_IE,Y_JB,Y_JE,
     C            XP,YP,XE,YN,DIEP,DJNP,DISE,DISN,
     C            AREP,ARNP,VOLP,
     C            IDATO,LVLGEO,IB,IE1,IE,JB,JE1,JE,ID,JD)
*
*  Routine to organize the creation of a cartesian grid, and the 
*  writing of several grid files.
*
***********************************************************************
*
      REAL XNE(ID),YNE(JD)
      REAL GRDXI,GRDX,GRDYI,GRDY
      REAL IRATX1,IRATX2,IRATY1,IRATY2
      REAL XP(ID),YP(JD),XE(ID),YN(JD)
      REAL DIEP(ID),DJNP(JD),DISE(ID),DISN(JD)
      REAL AREP(JD),ARNP(ID),VOLP(ID,JD)
      REAL X_IB,X_IE,Y_JB,Y_JE
*
      INTEGER IBW,IEW,JBW,JEW
      INTEGER IB,IE1,IE,JB,JE1,JE,I,J
      INTEGER LVLGEO
*
      OPEN(UNIT=77,FILE='corn.grd')
      OPEN(UNIT=78,FILE='pont.grd')
      OPEN(UNIT=79,FILE='xp.dat')
      OPEN(UNIT=80,FILE='yp.dat')
      OPEN(UNIT=81,FILE='xne.dat')
      OPEN(UNIT=82,FILE='yne.dat')
*
*  Initialization for uniform Cartesian grid (MME9710 Assign 5)
*
      IE1=IE
      JE1=JE
      GRDXI= GRDX
      GRDYI= GRDY
      IRATX1= 1.0
      IRATX2= 1.0
      IRATY1= 1.0
      IRATY2= 1.0
      IBW=0
      IEW=0
      JBW=0
      JEW=0
      X_IB=0.0
      X_IE=0.0
      Y_JB=0.0
      Y_JE=0.0
*
*  Define integration points and c.v. corners; generate grid
*
      CALL MAKGRD(XNE,YNE, GRDXI,GRDX,GRDYI,GRDY,
     C            IRATX1,IRATX2,IRATY1,IRATY2,
     C            IBW,IEW,JBW,JEW,X_IB,X_IE,Y_JB,Y_JE,
     C            IB,IE1,IE,JB,JE1,JE,ID,JD)
      CALL GRDGEO(XP,YP,XE,YN,DIEP,DJNP,DISE,DISN,
     C            AREP,ARNP,VOLP,
     C            IBW,IEW,JBW,JEW,
     C            XNE,YNE,IB,IE,JB,JE,ID,JD) 
*  
*  Create coordinate lists for post-processing
*
      WRITE(79,*) 'XP'
      WRITE(81,*) 'XNE'
      WRITE(80,*) 'YP'
      WRITE(82,*) 'YNE'
      DO 61 I=IB-1,IE+1
	    WRITE(79,'(1PE13.5)') XP(I)
	    WRITE(81,'(1PE13.5)') XNE(I)
 61   CONTINUE
      DO 62 I=JB-1,JE+1
	    WRITE(80,'(1PE13.5)') YP(I)
	    WRITE(82,'(1PE13.5)') YNE(I)
 62   CONTINUE
*
*  Create 2 .grd files; one containing NE corners and one
*  containing nodal points.
*
      WRITE(77,'(3(I3,2X))') IE-IB+2, JE-JB+2,1
      DO 63 J=JB-1,JE
	    DO 64 I=IB-1,IE
	      WRITE(77,'(3(1PE13.5,2X))') XNE(I),YNE(J),0.0D+00
 64     CONTINUE
	    WRITE(77,*) ' '
 63   CONTINUE
*
      WRITE(78,'(3(I3,2X))') IE-IB+3, JE-JB+3,1
      DO 65 J=JB-1,JE+1
	    DO 66 I=IB-1,IE+1
	      WRITE(78,'(3(1PE13.5,2X))') XP(I),YP(J),0.0D+00
 66     CONTINUE
	    WRITE(78,*) ' '
 65   CONTINUE
*
*  Produce output for LVLGEO >= 1
*
      IF (LVLGEO .GE. 1) THEN
        CALL OUT1D(XNE ,'XNE     ',IDATO,IB-1,IE  ,1,ID)
        CALL OUT1D(YNE ,'YNE     ',IDATO,JB-1,JE  ,1,JD)
	    CALL OUT1D(XP  ,'XP      ',IDATO,IB-1,IE+1,1,ID)
	    CALL OUT1D(YP  ,'YP      ',IDATO,JB-1,JE+1,1,JD)
	    CALL OUT2D(VOLP,'VOLP    ',IDATO,IB ,IE ,1,JB ,JE ,1,ID,JD)
	    CALL OUT1D(DIEP,'DIEP    ',IDATO,IB-1,IE+1,1,ID)
	    CALL OUT1D(DJNP,'DJNP    ',IDATO,JB-1,JE+1,1,JD)
        CALL OUT1D(DISE,'DISE    ',IDATO,IB-1,IE+1,1,ID)
        CALL OUT1D(DISN,'DISN    ',IDATO,JB-1,JE+1,1,JD)
      ENDIF
*
      RETURN
      END
