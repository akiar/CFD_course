*
*             file adcont.f
*********************************************************************
*
      SUBROUTINE ADCONT(AP,AW,AE,B,ACUW,ACUE,BC,DHUE,DPDX,
     C                  RHO,DIEP,DISE,IB,IE,ID)
*
*     Subroutine to construct coefficients for all terms in the
*     conservation of mass equation.  The form of the mass equation
*     is:
*                   ACUE*Ue + ACUW*Uw + BC = 0
*
*     The coefficients are stored in the first line of the 
*     coefficient blocks for mass and momentum.  The form of the 
*     equation is:
*
*               [A]p*{}p = [A]w*{}w + [A]e*{}e + [b]p
*
*********************************************************************
*
      REAL AP(2,2,ID),AW(2,2,ID),AE(2,2,ID),B(2,ID),ACUW(ID),ACUE(ID),
     C     BC(ID),DHUE(ID),DPDX(ID)
      REAL DIEP(ID),DISE(ID),RHO
      INTEGER IB,IE,ID,I
      PRINT *, "AP(1,1,I) AW(1,1,I) AE(1,1,I) 
     C          AP(1,2,I) AW(1,2,I) AE(1,2,I), B(1)"
*
*     BEGINNING NODE: IB
*
      AP(1,1,IB) = DHUE(IB)*ACUE(IB)/DIEP(IB)
      AP(1,2,IB) = 0.5*ACUE(IB)
      AW(1,1,IB) = 0
      AW(1,2,IB) = - ACUW(IB)
      AE(1,1,IB) = DHUE(IB)*ACUE(IB)/DIEP(IB)
      AE(1,2,IB) = - 0.5*ACUE(I)
      B(1,IB)    = - ACUE(IB)*(DHUE(IB)*(DPDX(IB)+DPDX(IB+1))/2)
     C             - ACUW(IB)*(DHUE(IB-1)*(DPDX(IB-1)+DPDX(IB))/2)
     C             - BC(IB)
      PRINT *, AP(1,1,IB), AW(1,1,IB), AE(1,1,IB), AP(1,2,IB), 
     C         AW(1,2,IB), AE(1,2,IB), B(1,IB)
*
*     INTERIOR NODES
*
      DO 10 I=IB+1,IE-1
        AP(1,1,I) = DHUE(I)*ACUE(I)/DIEP(I)
     C              - DHUE(I-1)*ACUW(I)/DIEP(I-1)
        AP(1,2,I) = 0.5*ACUE(I) + 0.5*ACUW(I)
        AW(1,1,I) = - DHUE(I-1)*ACUW(I)/DIEP(I-1)
        AW(1,2,I) = - 0.5*ACUW(I)
        AE(1,1,I) = DHUE(I)*ACUE(I)/DIEP(I)
        AE(1,2,I) = - 0.5*ACUE(I)
        B(1,I)    = - ACUE(I)*(DHUE(I)*(DPDX(I)+DPDX(I+1))/2)
     C              - ACUW(I)*(DHUE(I-1)*(DPDX(I-1)+DPDX(I))/2)
     C              - BC(I)
        PRINT *, AP(1,1,I), AW(1,1,I), AE(1,1,I), AP(1,2,I), 
     C           AW(1,2,I), AE(1,2,I), B(1,I)
 10	  CONTINUE
*
*     ENDING NODE: IE
*
      AP(1,1,IE) = -DHUE(IE-1)*ACUW(IE)/DIEP(IE-1)
      AW(1,1,IE) = -DHUE(IE-1)*ACUW(IE)/DIEP(IE-1)
      AE(1,1,IE) = 0
      AP(1,2,IE) = 0.5*ACUE(IE)
      AW(1,2,IE) = - 0.5*ACUW(IE)
      AE(1,2,IE) = - ACUE(IE)
      B(1,IE)    = - ACUE(IE)*(DHUE(IE)*(DPDX(IE)+DPDX(IE+1))/2)
     C             - ACUW(IE)*(DHUE(IE-1)*(DPDX(IE-1)+DPDX(IE))/2)
     C             - BC(IE)
      PRINT *, AP(1,1,IE), AW(1,1,IE), AE(1,1,IE), AP(1,2,IE), 
     C         AW(1,2,IE), AE(1,2,IE), B(1,IE)
      RETURN
      END
