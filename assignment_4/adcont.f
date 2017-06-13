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
      REAL ...
      INTEGER ...
*
      RETURN
      END
