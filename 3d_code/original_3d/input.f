*      This file contains 4 subroutines: FILDEF,FLDFIN, FLDFOT, and INPUT
*
*************************************************************************
*
      SUBROUTINE FILDEF(IDATI,IRSI,IRSO,ITERMO)
*
*     Subroutine to define unit numbers and open data files.
*
*     INTEGER IDATI  unit number for file with input data; output
*     INTEGER IRSI   unit number for file with binary restart data; input
*     INTEGER IRSO   unit number for file with binary output; output
*     INTEGER ITERMO unit number for output to terminal; output
*
*************************************************************************
*
      IMPLICIT NONE
      INTEGER IDATI,IRSI,IRSO,ITERMO
*     
*      FILE UNIT NUMBERS
*
      IDATI=71
      ITERMO=6
      IRSI=75
      IRSO=76
*
      OPEN(IDATI,FILE='in.dat')
      OPEN(IRSI,FORM='UNFORMATTED',FILE='rsi.bin')
      OPEN(IRSO,FORM='UNFORMATTED',FILE='rso.bin')
*
      RETURN
      END
*
*************************************************************************
*
      SUBROUTINE FLDFIN(IDATI,IRSI,ITERMO)
*
*     Subroutine to define unit numbers and open data files.
*
*     INTEGER IDATI  unit number for file with input data; output
*     INTEGER IRSI   unit number for file with binary restart data; input
*     INTEGER ITERMO unit number for output to terminal; output
*
*************************************************************************
*
      IMPLICIT NONE
      INTEGER IDATI,IRSI,ITERMO
*     
*      FILE UNIT NUMBERS
*
      IDATI=71
      IRSI=75
      ITERMO=6
*
      OPEN(IDATI,FILE='in.dat')
      OPEN(IRSI,FORM='UNFORMATTED',FILE='rsi.bin')
*      
      RETURN
      END
*
*************************************************************************
*
      SUBROUTINE FLDFOT(IRSO, FNAPP)
*
*     Subroutine to define unit numbers and open data files.
*
*     INTEGER IRSO   unit number for file with binary output; output
*
*************************************************************************
*
      IMPLICIT NONE
      INTEGER IRSO
      CHARACTER*5 FNAPP
*     
*      FILE UNIT NUMBERS
*
      IRSO=76
*
      OPEN(IRSO,FORM='UNFORMATTED',FILE='rso'//FNAPP//'.bin')
*
      RETURN
      END
*
*************************************************************************
*
      SUBROUTINE INPUT(IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C     KB,KE1,KE2,KE3,KE,GRDX1,GRDX2,GRDX3,GRDX,
     C     GRDY1,GRDY2,GRDY3,GRDY,GRDZ1,GRDZ2,GRDZ3,GRDZ,
     C     IRATX1,IRATX2,IRATX3,IRATX4,IRATY1,IRATY2,IRATY3,IRATY4,
     C     IRATZ1,IRATZ2,IRATZ3,IRATZ4,RHO,COND,CP,VISC,BETA,
     C     T0,U0,V0,W0,P0,UIN,DTMX,IRSTRT,ITIME,DTIME,PFTIME,KNTTM,
     C     KNTUVP,CRIT,TCRIT,LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM,
     C     POROUS,PCONV,EPS,PD,LD,VISCB,FORCH,PERM,CONDFE,CONDSX,
     C     CONDSY,CONDSZ,CSP,RHOSP,ASF,DEPSDZ,KS,IPB,IPE,JPB,JPE,
     C     KPB,KPE,IPB2,IPE2,JPB2,JPE2,KPB2,KPE2,SOLID,CONDS,CS,
     C     RHOS,ISB,ISE,JSB,JSE,KSB,KSE,ISB2,ISE2,JSB2,JSE2,
     C     KSB2,KSE2,FNAPP,IDATI)
*
*     Subroutine to read in input variables from a data file.
*
*     IB,IE  first and last interior indices in i; output
*     JB,JE  first and last interior indices in j; output
*     KB,KE  first and last interior indices in k; output
*     IE*,JE*,KE*   last indices of part domains; output
*     GRDX,GRDY,GRDZ  total dimensions of domain; output
*     GRDX*,GRDY*,GRDZ*  dimensions of part domains; output
*     IRATX*,IRATY*,  ratio if 1st to last CV in part domain; output
*     IRATZ*
*     RHO            fluid density; output
*     COND           thermal conductivity; output
*     CP             specific heat at constant pressure; output
*     VISC           viscosity of fluid; output
*     BETA           expansion coefficient; output
*     T0             initial constant temperature; output
*     U0,V0,W0       initial constant velocity; output
*     P0             initial constant pressure; output
*     UIN            u(x=0) for inlet velocity conditions; output
*     DTMX           heater to inlet temp. difference; output
*     IRSTRT         restart parameter; output
*     ITIME          initial real time; output
*     DTIME          time step ; output
*     PFTIME         printed frequency of output; output
*     KNTTM          maximum number of time steps;  output
*     KNTUVP         maximum number of u-v-p iterations: output
*     CRIT           stopping criterion for inner loop; output
*     TCRIT          stopping criterion for time loop; output
*     LVLGEO         parameter to control output of geometry; output
*     LVLCOF         parameter to control output of coef. ; output
*     LVLMGD         parameter to control output in solver; output
*     ADVSCM         parameter controlling advection scheme; output
*     ISOTHM         = 0, isothermal, = 1, solve energy; output
*     POROUS  = 0, no porous region, = 1, constant pptys, = 2, carbon 
*             foam; output
*     PCONV          Sets correlation, =0 AL, =1 POCO, =2 218
*     EPS            mean porosity of porous region; output
*     PD             pore diameter of porous medium; output
*     LD  ligament diameter of porous medium for POROUS = 1; output
*     VISCB  Brinkman effective viscosity; output     
*     PERM  permeability or coefficient for permeability correlation; 
*           output
*     FORCH  Forchheimer coefficient or coefficient for Forchheimer
*            correlation; output
*     CONDFE  fluid stagnant effective conductivity for POROUS = 1;
*             output
*     CONDSX,CONDSY,CONDSZ  components of effective consuctivity of 
*                           solid constituent for POROUS = 1; output
*     CSP  specific heat of solid constituent of porous region; output
*     RHOSP  density of solid constituent of porous region; output
*     ASF            specific surface are for POROUS = 1; output
*     DEPSDZ         porosity gradient in z-direction; input     
*     KS  material conductivity of solid constituent of foam; output     
*     IPB,IPE  first and last indices of porous region in i; output
*     JPB,JPE  first and last indices of porous region in j; output
*     KPB,KPE  first and last indices of porous region in k; output     
*     SOLID          = 0, no solid region, = 1, solid region; output
*     CONDS          thermal conductivity of pure solid region; output
*     CS             specific heat of pure solid region; output
*     RHOS           density of pure solid region; output
*     ISB,ISE  first and last indices of solid region in i; output
*     JSB,JSE  first and last indices of solid region in j; output
*     KSB,KSE  first and last indices of solid region in k; output
*     FNAPP          text to be appended to output filenames to 
*                    identify case; output
*     IDATI          unit number for file with input data; input
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,KB,KE1,KE2,KE3,KE
      INTEGER IRSTRT,KNTTM,KNTUVP,LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM
      INTEGER IDATI,PFTIME
      INTEGER POROUS,PCONV,IPB,IPE,JPB,JPE,KPB,KPE
      INTEGER IPB2,IPE2,JPB2,JPE2,KPB2,KPE2
      INTEGER SOLID,ISB,ISE,JSB,JSE,KSB,KSE
      INTEGER ISB2,ISE2,JSB2,JSE2,KSB2,KSE2
      REAL*8 GRDX1,GRDX2,GRDX3,GRDX,GRDY1,GRDY2,GRDY3,GRDY
      REAL*8 GRDZ1,GRDZ2,GRDZ3,GRDZ
      REAL*8 IRATX1,IRATX2,IRATX3,IRATX4,IRATY1,IRATY2,IRATY3,IRATY4
      REAL*8 IRATZ1,IRATZ2,IRATZ3,IRATZ4
      REAL*8 RHO,COND,CP,VISC,BETA,ITIME,DTIME,CRIT,TCRIT
      REAL*8 U0,V0,W0,P0,T0,UIN,DTMX
      REAL*8 EPS,PD,LD,VISCB,FORCH,PERM,CONDFE,CONDSX,CONDSY,CONDSZ
      REAL*8 CSP,RHOSP,ASF,CONDS,CS,RHOS,DEPSDZ,KS
      CHARACTER*5 FNAPP
*
      READ(IDATI,5000) IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C                 KB,KE1,KE2,KE3,KE
      READ(IDATI,5010) GRDX1,GRDX2,GRDX3,GRDX,GRDY1,GRDY2,GRDY3,GRDY,
     C                 GRDZ1,GRDZ2,GRDZ3,GRDZ
      READ(IDATI,5010) IRATX1,IRATX2,IRATX3,IRATX4,
     C                 IRATY1,IRATY2,IRATY3,IRATY4, 
     C                 IRATZ1,IRATZ2,IRATZ3,IRATZ4 
      READ(IDATI,5020) RHO,COND,CP,VISC,BETA
      READ(IDATI,5025) T0,U0,V0,W0,P0,UIN,DTMX
      READ(IDATI,5040) IRSTRT,ITIME,DTIME,PFTIME,KNTTM,KNTUVP,CRIT,
     C                 TCRIT
      READ(IDATI,5050) LVLGEO,LVLCOF,LVLMGD,ADVSCM,ISOTHM
      READ(IDATI,5060) POROUS,PCONV,EPS,PD,LD,VISCB,FORCH,PERM
      READ(IDATI,5070) CONDFE,CONDSX,CONDSY,CONDSZ,CSP,RHOSP,ASF
      READ(IDATI,5075) DEPSDZ,KS
      READ(IDATI,5080) IPB,IPE,JPB,JPE,KPB,KPE
      READ(IDATI,5080) IPB2,IPE2,JPB2,JPE2,KPB2,KPE2
      READ(IDATI,5090) SOLID,CONDS,CS,RHOS
      READ(IDATI,5080) ISB,ISE,JSB,JSE,KSB,KSE
      READ(IDATI,5080) ISB2,ISE2,JSB2,JSE2,KSB2,KSE2
      READ(IDATI,5100) FNAPP      
*      
 5000 FORMAT(I5/I5/I5/I5/I5/I5/I5/I5/I5/I5/I5/I5/I5/I5/I5)
 5010 FORMAT(E12.5/E12.5/E12.5/E12.5/E12.5/E12.5/E12.5/E12.5
     C       /E12.5/E12.5/E12.5/E12.5)
 5020 FORMAT(E12.5/E12.5/E12.5/E12.5/E12.5)
 5025 FORMAT(E12.5/E12.5/E12.5/E12.5/E12.5/E12.5/E12.5)
 5040 FORMAT(I5/E12.5/E12.5/I5/I5/I5/E12.5/E12.5)
 5050 FORMAT(I5/I5/I5/I5/I5)
 5060 FORMAT(I5/I5/E12.5/E12.5/E12.5/E12.5/E12.5/E12.5)
 5070 FORMAT(E12.5/E12.5/E12.5/E12.5/E12.5/E12.5/E12.5)
 5075 FORMAT(E12.5/E12.5)
 5080 FORMAT(I5/I5/I5/I5/I5/I5)
 5090 FORMAT(I5/E12.5/E12.5/E12.5)
 5100 FORMAT(A5) 
      RETURN
      END