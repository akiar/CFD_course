*
*     This file contains 2 subroutines: COEFF and COEFFM
*
***********************************************************************
*
      SUBROUTINE COEFF(AP,AW,AE,AS,AN,B,
     C                 ME,MN,DE,DN,Q,R,PHI,VOLP,
     C                 ALFAE,ALFAN,DIEP,DJNP,
     C                 RHO,DTIME,
     C                 IB,IE,JB,JE,ID,JD)
*
*     Subroutine to calculate the transport coefficients for the
*     variable PHI on a co-located orthogonal grid.  Used for
*     independently solved variables.
*
*     AP(ID,JD) active coefficient for P node; output
*     ASUM sum of active neighbour coefficients; output
*     AW(ID,JD) active coefficient for W node; output
*     AE(ID,JD) active coefficient for E node; output
*     AS(ID,JD) active coefficient for S node; output
*     AN(ID,JD) active coefficient for N node; output
*     B(ID,JD)  accumulated fixed source term; output
*
*     ME(ID,JD) mass flux for the east face;  input
*     MN(ID,JD) mass flux for the north face; input
*     DE(ID,JD) diffusion coefficient for east face; input
*     DN(ID,JD) diffusion coefficient for north face; input
*     Q(ID,JD) fixed source term; input
*     R(ID,JD) linearized source term; input
*     PHI(ID,JD) estimate of phi at old time; input
*     VOLP(ID,JD) c.v. volume; input
*     ALFAE(ID,JD) advection weight factor for east face; input
*     ALFAN(ID,JD) advection weight factor for north face: input
*     DIEP(JD) distance from P to E through e; input
*     DJNP(JD) distance from P to N through n; input
*     RHO  fluid density; input
*     DTIME  time step for implicit scheme; input
*
*     INTEGER IB,IE first and last interior indices in i; input
*     INTEGER JB,JE first and last interior indices in j; input
*     INTEGER ID,JD array dimensions; input     
*
*     Notes: 1) This subroutine should be proceeded by calls to
*                 NULL(B, ...
*                 DIFPHI(DE,DN, ...
*                 MASFLX(ME,MN,..
*                 WEIGHT(ALFAE,ALFAN,...
*                 SRC(Q,R, ...
*
*
**********************************************************************
*
      REAL ASUM,AW(ID,JD),AE(ID,JD),AS(ID,JD)
      REAL AN(ID,JD),AP(ID,JD),B(ID,JD),PHI(ID,JD)
      REAL ME(ID,JD),MN(ID,JD),DE(ID,JD),DN(ID,JD)
      REAL Q(ID,JD),R(ID,JD),VOLP(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD)
      REAL DIEP(ID),DJNP(JD)
      REAL RHO,DTIME
      INTEGER IB,IE,JB,JE,ID,JD
      INTEGER I,J
*     
      DO 20 J=JB,JE
        DO 10 I=IB,IE
*
          AW(I,J)=DE(I-1,J)
     C           +0.5*ME(I-1,J)*(1.0+ALFAE(I-1,J))
          AE(I,J)=DE(I,J)
     C           -0.5*ME(I  ,J)*(1.0-ALFAE(I  ,J))
          AS(I,J)=DN(I,J-1)
     C           +0.5*MN(I,J-1)*(1.0+ALFAN(I,J-1))
          AN(I,J)=DN(I,J)
     C           -0.5*MN(I,J  )*(1.0-ALFAN(I,J  ))
          ASUM=   AW(I,J)+AE(I,J)+AS(I,J)+AN(I,J)
          AP(I,J)=RHO*VOLP(I,J)/DTIME+ASUM-R(I,J)
*
          B(I,J)=B(I,J)+Q(I,J)
     C          +PHI(I,J)*RHO*VOLP(I,J)/DTIME
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE COEFFM(AP,AW,AE,AS,AN,B,
     C                  ME,MN,DE,DN,Q,R,PHI,VOLP,
     C                  ALFAE,ALFAN,DIEP,DJNP,
     C                  RHO,DTIME,
     C                  IV,IB,IE,JB,JE,N,ID,JD)
*
*     Subroutine to calculate the transport coefficients for the
*     variable PHI on a co-located non-orthogonal grid.  Used for
*     N simultaneous equations.
*
*     AP(N,N,ID,JD) active coefficient for P node; output
*     ASUM sum of active neighbour coefficients; output
*     AW(N,N,ID,JD) active coefficient for W node; output
*     AE(N,N,ID,JD) active coefficient for E node; output
*     AS(N,N,ID,JD) active coefficient for S node; output
*     AN(N,N,ID,JD) active coefficient for N node; output
*     B(N,ID,JD)    accumulated fixed source term; output
*
*     ME(ID,JD) mass flux for the east face;  input
*     MN(ID,JD) mass flux for the north face; input
*     DE(ID,JD) diffusion coefficient for east face; input
*     DN(ID,JD) diffusion coefficient for north face; input
*     Q(ID,JD) fixed source term; input
*     R(ID,JD) linearized source term; input
*     PHI(ID,JD) estimate of phi at old time; input
*     VOLP(ID,JD) c.v. volume; input
*     ALFAE(ID,JD) advection weight factor for east face; input
*     ALFAN(ID,JD) advection weight factor for north face: input
*     DIEP(JD) distance from P to E through e; input
*     DJNP(JD) distance from P to N through n; input
*     RHO  fluid density; input
*     DTIME  time step for implicit scheme; input
*
*     INTEGER IV      position in block; input
*     INTEGER IB,IE   first and last interior indices in i; input
*     INTEGER JB,JE   first and last interior indices in j; input
*     INTEGER N       number of variables in solution; input
*     INTEGER IJ      dimension of RS array; input
*     INTEGER ID,JD   grid-array dimensions; input     
*
*     Notes: 1) This subroutine should be proceeded by calls to
*                 NULL(B, ...
*                 DIFPHI(DE,DN, ...
*                 MASFLX(ME,MN,..
*                 WEIGHT(ALFAE,ALFAN,...
*                 SRC(Q,R, ...
*
**********************************************************************
*
      REAL ASUM,AW(N,N,ID,JD),AE(N,N,ID,JD)
      REAL AS(N,N,ID,JD),AN(N,N,ID,JD),AP(N,N,ID,JD)
      REAL B(N,ID,JD),PHI(ID,JD)
      REAL ME(ID,JD),MN(ID,JD),DE(ID,JD),DN(ID,JD)
      REAL Q(ID,JD),R(ID,JD),VOLP(ID,JD)
      REAL ALFAE(ID,JD),ALFAN(ID,JD)
      REAL DIEP(ID),DJNP(JD)
      REAL RHO,DTIME
      INTEGER IV,IB,IE,JB,JE,N,ID,JD
      INTEGER I,J
*
*
      DO 20 J=JB,JE
        DO 10 I=IB,IE
*
          AW(IV,IV,I,J)=DE(I-1,J)
     C                 +0.5*ME(I-1,J)*(1.0+ALFAE(I-1,J))
          AE(IV,IV,I,J)=DE(I  ,J)
     C                 -0.5*ME(I  ,J)*(1.0-ALFAE(I  ,J))
          AS(IV,IV,I,J)=DN(I,J-1)
     C                 +0.5*MN(I,J-1)*(1.0+ALFAN(I,J-1))
          AN(IV,IV,I,J)=DN(I,J)
     C                 -0.5*MN(I,J  )*(1.0-ALFAN(I,J  ))
          ASUM= AW(IV,IV,I,J)+AE(IV,IV,I,J)+AS(IV,IV,I,J)
     C              +AN(IV,IV,I,J)
          AP(IV,IV,I,J)=RHO*VOLP(I,J)/DTIME+ASUM-R(I,J)
*
          B(IV,I,J)=B(IV,I,J)+Q(I,J)
     C             +PHI(I,J)*RHO*VOLP(I,J)/DTIME
*
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END
