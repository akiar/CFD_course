*
*             file adcont.f
*********************************************************************
*
      SUBROUTINE ADCONT(AP,AW,AE,AS,AN,B,
     C                  ACUW,ACUE,ACVS,ACVN,DHUE,DHVN,
     C                  DPDX,DPDY,RHO,
     C                  DIEP,DISE,DJNP,DISN,
     C                  IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to construct coefficients for all terms in the
*     conservation of mass equation.  The form of the mass equation
*     is:
*             Acue*Ue + Acuw*Uw + Acvn*Vn + Acvs*Vs = 0
*
*     The coefficients are stored in the same matrix arrays as the
*     other equations.  The form of those equations is:
*
*      [A]p*{}p = [A]w*{}w + [A]e*{}e + [A]s*{}s + [A]n*{}n
*                                                   + [b]p + [S]p
*
*     thus the terms on {}p are moved to the other side of the 
*     equation.  The order of the blocks is P,U,V.
*
*********************************************************************
*
      REAL AP(N,N,ID,JD),AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),B(N,ID,JD)
      REAL ACUW(ID,JD),ACUE(ID,JD),ACVS(ID,JD),ACVN(ID,JD)
      REAL DPDX(ID,JD),DPDY(ID,JD)
      REAL DHUE(ID,JD),DHVN(ID,JD),RHO
      REAL DIEP(ID),DISE(ID),DJNP(JD),DISN(JD)
      INTEGER IB,IE,JB,JE,I,J,N,ID,JD
*  
      DO 1 I=IB,IE
       DO 2 J=JB,JE
*
*--Add components for west face of control-volume
*
       IF(I.EQ.IB) THEN
	    AW(1,2,I,J)= AW(1,2,I,J)+ACUW(I,J)
       ELSE
        AP(1,1,I,J)= AP(1,1,I,J)+ACUW(I,J)*DHUE(I-1,J)/DIEP(I-1)
        AW(1,1,I,J)= AW(1,1,I,J)+ACUW(I,J)*DHUE(I-1,J)/DIEP(I-1)
*
        AP(1,2,I,J)= AP(1,2,I,J)-ACUW(I,J)*DISE(I-1)/DIEP(I-1)
        AW(1,2,I,J)= AW(1,2,I,J)+ACUW(I,J)*DISE(I)/DIEP(I-1)
*
        B(1,I,J)= B(1,I,J)+ACUW(I,J)*( DHUE(I-1,J)*
     C           (DISE(I)*DPDX(I-1,J)+DISE(I-1)*DPDX(I,J))/DIEP(I-1) )
       ENDIF
*
*--Add components for east face of control-volume
*
       IF(I.EQ.IE) THEN
	    AE(1,2,I,J)= AE(1,2,I,J)+ACUE(I,J)
       ELSE
	    AP(1,1,I,J)= AP(1,1,I,J)-ACUE(I,J)*DHUE(I,J)/DIEP(I)
	    AE(1,1,I,J)= AE(1,1,I,J)-ACUE(I,J)*DHUE(I,J)/DIEP(I)
*
        AP(1,2,I,J)= AP(1,2,I,J)-ACUE(I,J)*DISE(I+1)/DIEP(I)
        AE(1,2,I,J)= AE(1,2,I,J)+ACUE(I,J)*DISE(I)/DIEP(I)
*
        B(1,I,J)= B(1,I,J)+ACUE(I,J)*( DHUE(I,J)*
     C           (DISE(I+1)*DPDX(I,J)+DISE(I)*DPDX(I+1,J))/DIEP(I) )
       ENDIF
*
*--Add components for south face of control-volume
*
       IF(J.EQ.JB) THEN
	    AS(1,3,I,J)=AS(1,3,I,J)+ACVS(I,J)
       ELSE
	    AP(1,1,I,J)= AP(1,1,I,J)+ACVS(I,J)*DHVN(I,J-1)/DJNP(J-1)
	    AS(1,1,I,J)= AS(1,1,I,J)+ACVS(I,J)*DHVN(I,J-1)/DJNP(J-1)
*
        AP(1,3,I,J)= AP(1,3,I,J)-ACVS(I,J)*DISN(J-1)/DJNP(J-1)
        AS(1,3,I,J)= AS(1,3,I,J)+ACVS(I,J)*DISN(J)/DJNP(J-1)
*
        B(1,I,J)= B(1,I,J)+ACVS(I,J)*( DHVN(I,J-1)*
     C           (DISN(J)*DPDY(I,J-1)+DISN(J-1)*DPDY(I,J))/DJNP(J-1) )
       ENDIF
*
*--Add components for north face of control-volume
*
       IF(J.EQ.JE) THEN
	    AN(1,3,I,J)=AN(1,3,I,J)+ACVN(I,J)
       ELSE 
	    AP(1,1,I,J)= AP(1,1,I,J)-ACVN(I,J)*DHVN(I,J)/DJNP(J)
	    AN(1,1,I,J)= AN(1,1,I,J)-ACVN(I,J)*DHVN(I,J)/DJNP(J)
*
        AP(1,3,I,J)= AP(1,3,I,J)-ACVN(I,J)*DISN(J+1)/DJNP(J)
        AN(1,3,I,J)= AN(1,3,I,J)+ACVN(I,J)*DISN(J)/DJNP(J)
*
        B(1,I,J)= B(1,I,J)+ACVN(I,J)*( DHVN(I,J)*
     C           (DISN(J+1)*DPDY(I,J)+DISN(J)*DPDY(I,J+1))/DJNP(J) )
       ENDIF
*
  2    CONTINUE
  1   CONTINUE
*
      RETURN
      END
