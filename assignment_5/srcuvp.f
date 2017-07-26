*
*         file srcuvp.f
*********************************************************************
*
      SUBROUTINE SRCUVP(AP,AW,AE,AS,AN,B, 
     C                 VOLP,DIEP,DISE,DJNP,DISN,
     C                 IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to add the pressure source term in the U & V
*     momentum equations.  The order of the solver is P,U,V.
*
*********************************************************************
*
      REAL AP(N,N,ID,JD),AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),B(N,ID,JD)
      REAL VOLP(ID,JD),DIEP(ID),DISE(ID),DJNP(JD),DISN(JD)
      INTEGER IB,IE,JB,JE,I,J,N,ID,JD
*  
      DO 1 I=IB,IE
       DO 2 J=JB,JE
*
        AP(2,1,I,J)= 0.5*VOLP(I,J)*(1./DIEP(I-1)-1./DIEP(I))
        AW(2,1,I,J)= 0.5*VOLP(I,J)/DIEP(I-1)
        AE(2,1,I,J)=-0.5*VOLP(I,J)/DIEP(I)
        AS(2,1,I,J)= 0.0
        AN(2,1,I,J)= 0.0
*
	AP(3,1,I,J)= 0.5*VOLP(I,J)*(1./DJNP(J-1)-1./DJNP(J))
	AW(3,1,I,J)= 0.0
	AE(3,1,I,J)= 0.0
	AS(3,1,I,J)= 0.5*VOLP(I,J)/DJNP(J-1)
	AN(3,1,I,J)=-0.5*VOLP(I,J)/DJNP(J) 
*
  2    CONTINUE
  1   CONTINUE
*
      RETURN
      END
