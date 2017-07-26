*      This file contains 3 subroutines: FILDEF, INPUT and ECHO
*
*************************************************************************
*
      SUBROUTINE FILDEF(IDATI,IRSI,IDATO,IRSO,ITERMO)
*
*     Subroutine to define unit numbers and open data files.
*
*     INTEGER IDATI  unit number for file with input data; output
*     INTEGER IRSI   unit number for file with binary restart data; input
*     INTEGER IDATO  unit number for file with print output;output
*     INTEGER IRSO   unit number for file with binary output; output
*     INTEGER ITERMO unit number for output to terminal; output
*
*************************************************************************
*
      INTEGER IDATI,IRSI,IDATO,IRSO,ITERMO,ISOLV
*     
*      FILE UNIT NUMBERS
*
      IDATI=71
      IDATO=72
      ITERMO=6
      ISOLV=73
      IRSI=75
      IRSO=76
*
      OPEN(IDATI,FILE='in.dat')
      OPEN(IDATO,FILE='out.dat')
      OPEN(ISOLV,FILE='dlap.dat')
      OPEN(IRSI,FORM='UNFORMATTED',FILE='rsi.bin')
      OPEN(IRSO,FORM='UNFORMATTED',FILE='rso.bin')
      RETURN
      END
*
*
*************************************************************************
*
      SUBROUTINE INPUT(IB,IE,JB,JE,
     C            GRDX,GRDY,
     C            RHO,COND,CP,VISC,BETA,
     C            T0,U0,V0,P0,IRSTRT,ITIME,
     C            DTIME,PFTIME,KNTTM,KNTUVP,CRIT,
     C            LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM,IDATI)
*
*     Subroutine to read in input variable from a data file.
*
*     INTEGER IB,IE  first and last interior indices in i; output
*     INTEGER JB,JE  first and last interior indices in j; output
*     GRDX,GRDY      total dimension of domain in 'x' and 'y'; input
*     IRATX,IRATY    ratio if 1st to last c.v. in 'x' and 'y'; input
*     RHO            fluid density; output
*     COND           thermal conductivity; output
*     CP             specific heat at constant pressure; output
*     VISC           viscosity of fluid; output
*     BETA           expansion coefficient; output
*     T0             initial constant temperature; output
*     U0,V0          initial constant velocity; output
*     P0             initial constant pressure; output
*     INTEGER IRSTRT restart parameter; output
*     ITIME          initial real time; output
*     DTIME          time step ; output
*     INTEGER PFTIME printed frequency of output; output
*     INTEGER KNTTM  maximum number of time steps;  output
*     INTEGER KNTUVP maximum number of u-v-p iterations: output
*     CRIT           stopping criterion for inner loop; output
*     INTEGER LVLGEO parameter to control output of geometry; output
*     INTEGER LVLCOF parameter to control output of coef. ; output
*     INTEGER LVLMGD parameter to control output in solver; output
*
*     INTEGER IDATI unit number for file with input data; input
*
************************************************************************
*
      INTEGER IB,IE1,IE,JB,JE1,JE,IRSTRT,KNTTM,KNTUVP
      INTEGER LVLGEO,LVLCOF,LVLMGD,ADVSCM
      INTEGER IDATI,PFTIME
      REAL GRDXI,GRDX,GRDYI,GRDY
      REAL IRATX1,IRATX2,IRATY1,IRATY2
      REAL RHO,COND,CP,VISC,BETA,ITIME,DTIME,CRIT
      REAL U0,V0,P0,T0
*
      READ(IDATI,5000) IB,IE,JB,JE
      READ(IDATI,5010) GRDX,GRDY
      READ(IDATI,5020) RHO,COND,CP,VISC,BETA
      READ(IDATI,5030) T0,U0,V0,P0
      READ(IDATI,5040) IRSTRT,ITIME,DTIME,PFTIME,KNTTM,KNTUVP,CRIT
      READ(IDATI,5050) LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM
 5000 FORMAT(I5/I5/I5/I5)
 5010 FORMAT(E12.5/E12.5)
 5020 FORMAT(E12.5/E12.5/E12.5/E12.5/E12.5)
 5030 FORMAT(E12.5/E12.5/E12.5/E12.5)
 5040 FORMAT(I5/E12.5/E12.5/I5/I5/I5/E12.5)
 5050 FORMAT(I5/I5/I5/I5/I5)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE ECHO(IB,IE,JB,JE,
     C           GRDX,GRDY,
     C           RHO,COND,CP,VISC,BETA,
     C           T0,U0,V0,P0,IRSTRT,ITIME,
     C           DTIME,PFTIME,KNTTM,KNTUVP,CRIT,
     C           LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM,IDATO)
*
*     Subroutine to echo the input parameters read in from file.
*
*     INTEGER IDATO  unit number for file with formatted output; input
*
*     See INPUT for other variable definitions
*************************************************************************
*
      INTEGER IB,IE1,IE,JB,JE1,JE,IRSTRT,KNTTM,KNTUVP
      INTEGER LVLGEO,LVLCOF,LVLMGD,ADVSCM
      INTEGER IDATO,PFTIME
      REAL GRDXI,GRDX,GRDYI,GRDY
      REAL IRATX1,IRATX2,IRATY1,IRATY2
      REAL RHO,COND,CP,VISC,BETA,DTIME,ITIME,CRIT
      REAL U0,V0,P0,T0
*
      WRITE(IDATO,5990) 
      WRITE(IDATO,6000) IB,IE,JB,JE
      WRITE(IDATO,6010) GRDX,GRDY
      WRITE(IDATO,6020) RHO,COND,CP,VISC,BETA
      WRITE(IDATO,6030) T0,U0,V0,P0
      WRITE(IDATO,6040) IRSTRT,ITIME,DTIME,PFTIME,KNTTM,KNTUVP,CRIT
      WRITE(IDATO,6050) LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM
 5990 FORMAT(72('='),/
     C      ,'==',/
     C      ,'==   F L O W 2 D   Output File',/
     C      ,'==',/
     C      ,72('='),//)
 6000 FORMAT(' Grid Parameters:',//
     C      ,'              IB = ',I5,/
     C      ,'              IE = ',I5,/
     C      ,'              JB = ',I5,/
     C      ,'              JE = ',I5,/,' ')
 6010 FORMAT(' Geometry Parameters:',//
     C      ,'            GRDX = ',1PD13.5,/
     C      ,'            GRDY = ',1PD13.5,/,' ')
 6020 FORMAT(' Fluid Properties:',//
     C      ,'             RHO = ',1PD13.5,/
     C      ,'            COND = ',1PD13.5,/
     C      ,'              CP = ',1PD13.5,/
     C      ,'            VISC = ',1PD13.5,/
     C      ,'            BETA = ',1PD13.5,/,' ')
 6030 FORMAT(' Initial Conditions:',//
     C      ,'              T0 = ',1PD13.5,/
     C      ,'              U0 = ',1PD13.5,/
     C      ,'              V0 = ',1PD13.5,/
     C      ,'              P0 = ',1PD13.5,/,' ')
 6040 FORMAT(' Solution Parameters:',//
     C      ,'          IRSTRT = ',I5,/
     C      ,'           ITIME = ',1PD13.5,/
     C      ,'           DTIME = ',1PD13.5,/
     C      ,'          PFTIME = ',I5,/
     C      ,'           KNTTM = ',I5,/
     C      ,'          KNTUVP = ',I5,/
     C      ,'            CRIT = ',1PD13.5,/,' ')
 6050 FORMAT(' Output Parameters:',//
     C      ,'          LVLGEO = ',I5,/
     C      ,'          LVLCOF = ',I5,/
     C      ,'          LVLMGD = ',I5,/
     C      ,'          ADVSCM = ',I5,/
     C      ,'          ISOTHM = ',I5,//,' ')
      RETURN
      END
