*
*          file gradp.f
************************************************************************
*
      SUBROUTINE GRADP(DPDX,DPDY, P,DIEP,DJNP,IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the pressure gradient at each nodal
*     location using the most current pressure field.
*
*     DPDX(ID,JD)  x-derivative of pressure at nodal location; output
*     DPDY(ID,JD)  y-derivative of pressure at nodal location; output
*
*     P(ID,JD)     current pressure field; input
*     DIEP(ID)  distance from P to E nodes; input
*     DJNP(JD)  distance from P to N nodes; input
*     IB,IE     first and last interior indices in i; input
*     JB,JE     first and last interior indices in j; input
*     ID,JD     array dimensions; input     
*
***********************************************************************
*
      REAL DPDX(ID,JD),DPDY(ID,JD),P(ID,JD),DIEP(ID),DJNP(JD)
      INTEGER IB,IE,JB,JE,ID,JD,I,J
*
*--Derivatives of P along JB-1 line
*
      J=JB-1
      DO 1 I=IB,IE
        DPDX(I,J)= 0.5*( (P(I+1,J)-P(I,J))/DIEP(I)
     C                  +(P(I,J)-P(I-1,J))/DIEP(I-1) )
        DPDY(I,J)= ( (P(I,J+1)-P(I,J))/DJNP(J) )
  1   CONTINUE
*
      DO 2 J=JB,JE
*
*--Derivatives of P at IB-1 node
*
        I=IB-1
        DPDX(I,J)= ( (P(I+1,J)-P(I,J))/DIEP(I) )
        DPDY(I,J)= 0.5*( (P(I,J+1)-P(I,J))/DJNP(J)
     C                  +(P(I,J)-P(I,J-1))/DJNP(J-1) )
*
*--Derivatives of P at interior nodes
*
        DO 3 I=IB,IE
          DPDX(I,J)= 0.5*( (P(I+1,J)-P(I,J))/DIEP(I)
     C                    +(P(I,J)-P(I-1,J))/DIEP(I-1) )
          DPDY(I,J)= 0.5*( (P(I,J+1)-P(I,J))/DJNP(J)
     C                    +(P(I,J)-P(I,J-1))/DJNP(J-1) )
  3     CONTINUE
*
*--Derivatives of P at IE+1 node
*
        I=IE+1
        DPDX(I,J)= ( (P(I,J)-P(I-1,J))/DIEP(I-1) )
        DPDY(I,J)= 0.5*( (P(I,J+1)-P(I,J))/DJNP(J)
     C                  +(P(I,J)-P(I,J-1))/DJNP(J-1) )
*
  2   CONTINUE
*
*--Derivatives of P along JE+1 line
*
      J=JE+1
      DO 4 I=IB,IE
        DPDX(I,J)= 0.5*( (P(I+1,J)-P(I,J))/DIEP(I)
     C                  +(P(I,J)-P(I-1,J))/DIEP(I-1) )
        DPDY(I,J)= ( (P(I,J)-P(I,J-1))/DJNP(J-1) )
  4   CONTINUE
*
      RETURN
      END

