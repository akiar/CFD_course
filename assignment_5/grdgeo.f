*
*      file grdgeo.f
************************************************************************
      SUBROUTINE GRDGEO(XP,YP,XE,YN,DIEP,DJNP,DISE,DISN,
     C                  AREP,ARNP,VOLP,
     C                  IBW,IEW,JBW,JEW,
     C                  XNE,YNE,IB,IE,JB,JE,ID,JD) 
*
*     Subroutine to calculate the geometrical properties for a simple
*     plane orthogonal grid. All properties are based on control
*     volumes obtained by joining the corner points with straight lines.
*
*     XP(ID) x location of P node for each c.v.; output
*     YP(JD) y location of P node for each c.v.; output
*     XE(ID) x location of e integration point;  output
*     YN(JD) y location of n integration point;  output
*
*     DIEP(ID) distance from P to E through e; output
*     DJNP(JD) distance from P to N through n; output
*     
*     DISE(ID) distance in i from P to e face; output
*     DISN(JD) distance in j from P to n face; output
*
*     AREP(JD) c.v. area of face at e point; output
*     ARNP(ID) c.v. area of face at n point; output
*     VOLP(ID,JD) volume of c.v. around P node; output
*
*     XNE(ID) x location of ne corner points; input
*     YNE(JD) y location of ne corner points; input
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input
*
*     Notes: 1) This is setup for an orthogonal grid.
*
*            2) The locations of the e and n integration points should
*               be calculated as the face midpoints and stored in arrays.
*
************************************************************************
*
      REAL XP(ID),YP(JD)
      REAL XE(ID),YN(JD)
      REAL DIEP(ID),DJNP(JD)
      REAL DISE(ID),DISN(JD)
      REAL AREP(JD),ARNP(ID),VOLP(ID,JD)
      REAL XNE(ID),YNE(JD)
      INTEGER IBW,IEW,JBW,JEW
      INTEGER IB,IE,JB,JE,IBM1,IEP1,JBM1,JEP1
      INTEGER I,J,ID,JD
*
      IBM1=IB-1
      IEP1=IE+1
      JBM1=JB-1
      JEP1=JE+1
*
*     Calculate locations of P nodes for interior and boundary c.v.
*
      XP(IBM1)=XNE(IBM1)
      DO 10 I=IB,IE
        XP(I)=0.5*(XNE(I)+XNE(I-1))
   10 CONTINUE
      XP(IEP1)=XNE(IE)
*
      YP(JBM1)=YNE(JBM1)
      DO 20 J=JB,JE
        YP(J)=0.5*(YNE(J)+YNE(J-1))
   20 CONTINUE
      YP(JEP1)=YNE(JE)
*
*  Compute grid geometry parameters
*
      DO 50 I=IBM1,IE
        DIEP(I)=XP(I+1)-XP(I)
 50   CONTINUE
*
      DO 60 J=JBM1,JE
        DJNP(J)=YP(J+1)-YP(J)
 60   CONTINUE
*
*  Calculate midface locations
*     
      DO 70 I=IBM1,IE
        XE(I)=XNE(I)
 70   CONTINUE
      XE(IEP1)=XNE(IE)
*
      DO 80 J=JBM1,JE
        YN(J)=YNE(J)
 80   CONTINUE
      YN(JEP1)=YNE(JE)
*
      DO 100 I=IBM1,IEP1
        DISE(I)=XE(I)-XP(I)
 100  CONTINUE
*
      DO 120 J=JBM1,JEP1
        DISN(J)=YN(J)-YP(J)
 120  CONTINUE
*
*  Calculate planar areas and volumes for Cartesian grid
*
      DO 200 I=IB,IE
        ARNP(I)=2.0*DISE(I)
 200  CONTINUE
      DO 210 J=JB,JE
        AREP(J)=2.0*DISN(J)
 210  CONTINUE
*     
      DO 230 J=JB,JE
      DO 220 I=IB,IE
        VOLP(I,J)=AREP(J)*ARNP(I)
 220  CONTINUE
 230  CONTINUE
*     
      RETURN
      END
