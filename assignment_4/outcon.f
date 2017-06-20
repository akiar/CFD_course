
      SUBROUTINE OUTCON(A, CVAR,IUNIT,IB,IE,ISKP,ID)
************************************************************************
*     Write output of active coefficients for momentum and conservation equations
************************************************************************
      REAL A(2,2,ID)
      INTEGER IUNIT,IB,IE,ISKP,ID
      INTEGER NWORD,IPNT(100),LI,LINE,I1,IN,IMAX
      CHARACTER*8 CVAR
      INTEGER I,J,K
*     
      PARAMETER (NWORD=7)
*
      DO 2000 J=1,2
        DO 200 K=1,2
          LI=0
          DO 10 I=IB,IE,ISKP
            LI=LI+1
            IPNT(LI)=I
 10       CONTINUE
          IMAX=LI
*
          WRITE(IUNIT,104)CVAR
          LINE=0
 30       LINE=LINE+1
          IN=LINE*NWORD
          I1=IN-NWORD+1
          IN=AMIN0(IN,IMAX)
          WRITE(IUNIT,100)(A(J,K,IPNT(I)),I=I1,IN)
          WRITE(IUNIT,101)(IPNT(I),I=I1,IN)
          IF(IMAX .GT. IN)GO TO 30
 100      FORMAT(' ',100(1PE11.3))
 101      FORMAT(' I-> ',3X,100(I4,7X))
 104      FORMAT(' ',/,' ',A8)
 200    CONTINUE
 2000 CONTINUE
      RETURN
      END
      