*
*         file main.f
************************************************************************
*
*  Program MAIN solves 3-dimensional laminar flows on 3-dimensional 
*  cartesian grids.
*
*  The field equations are discretized using the finite-volume 
*  formulation of Patankar, but using collocated variable arrangement.
*  The default advection scheme is 1st order Pecelet weighting.
*
***********************************************************************
*
*============================
*  Declaration of Variables
*============================
*
      IMPLICIT NONE
      INTEGER ID,JD,KD,N,NT,NNB,ITER,A,B,C
      REAL*8 GEE,PI,BLEND,UREF,DBL1,DBL0,RLX,TIN
      PARAMETER(ID=82,JD=82,KD=22)
      PARAMETER(N=4)
      PARAMETER(NT=2)
      PARAMETER(NNB=6)
      PARAMETER(ITER=10)
      PARAMETER(GEE=9.80665)
      PARAMETER(PI=3.14159265)
      PARAMETER(UREF=1.0)
      PARAMETER(BLEND=1.0)
      PARAMETER(DBL1=1.0,DBL0=0.0)
      PARAMETER(RLX=0.0)
      PARAMETER(TIN=293.0) 
*  
*--Arrays for geometry description
*
      REAL*8 XNET(ID),YNET(JD),ZNET(KD)
      REAL*8 XP(ID),YP(JD),ZP(KD),XE(ID),YN(JD),ZT(KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DISE(ID),DISN(JD),DIST(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD),VOLP(ID,JD,KD)
      REAL*8 X_IB,X_IE,Y_JB,Y_JE,Z_KB,Z_KE
*
*--Arrays for active coefficients and source terms
*
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      REAL*8 DEF(ID,JD,KD),DNF(ID,JD,KD),DTF(ID,JD,KD)
      REAL*8 DES(ID,JD,KD),DNS(ID,JD,KD),DTS(ID,JD,KD)            
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 ALFAE(ID,JD,KD),ALFAN(ID,JD,KD),ALFAT(ID,JD,KD)
      REAL*8 QTF(ID,JD,KD),RTF(ID,JD,KD)
      REAL*8 QTS(ID,JD,KD),RTS(ID,JD,KD)      
      REAL*8 QU(ID,JD,KD),RU(ID,JD,KD),QV(ID,JD,KD),RV(ID,JD,KD)
      REAL*8 QW(ID,JD,KD),RW(ID,JD,KD),BC(ID,JD,KD)
*
      REAL*8 ACUW(ID,JD,KD),ACUE(ID,JD,KD),ACVS(ID,JD,KD)
      REAL*8 ACVN(ID,JD,KD),ACWT(ID,JD,KD),ACWB(ID,JD,KD)
*
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD)
      REAL*8 ATS(NT,NT,ID,JD,KD),ATN(NT,NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD)
      REAL*8 AUB(N,N,ID,JD,KD),AUT(N,N,ID,JD,KD)
      REAL*8 AUP(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
*
*--Arrays of active variables
*
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 P(ID,JD,KD),TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 DPDX(ID,JD,KD),DPDY(ID,JD,KD),DPDZ(ID,JD,KD)
      REAL*8 DPDXM(ID,JD,KD),DPDYM(ID,JD,KD),DPDZM(ID,JD,KD)
      REAL*8 DPDXP(ID,JD,KD),DPDYP(ID,JD,KD),DPDZP(ID,JD,KD)            
      REAL*8 UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD)
      REAL*8 TFOLD(ID,JD,KD),TSOLD(ID,JD,KD),UOLD(ID,JD,KD)
      REAL*8 VOLD(ID,JD,KD),WOLD(ID,JD,KD),POLD(ID,JD,KD)
      REAL*8 U0,V0,W0,P0,T0,UIN,DTMX
*
*--Arrays of coefficients for coupling mass-momentum equations
*
      REAL*8 DHUE(ID,JD,KD),DHVN(ID,JD,KD),DHWT(ID,JD,KD)
      REAL*8 VOLE(ID,JD,KD),VOLN(ID,JD,KD),VOLT(ID,JD,KD)
      REAL*8 AEAST(ID,JD,KD),ANORTH(ID,JD,KD),ATOP(ID,JD,KD)
*
*--Arrays for reporting convergence of equations
*
      REAL*8 NORM(6),RSD(ID,JD,KD),ARSD(6),MAXRSD
      REAL*8 RSDMAX(6),RSDOLD(6),CRIT
      REAL*8 TRMAX(6),MAXTR,TCRIT
*
*--Declaration of all other variables
*
      INTEGER I,J,K,L,IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE
      INTEGER KB,KE1,KE2,KE3,KE
      INTEGER IDATI,IRSI,IRSO,ITERMO
      INTEGER IRSTRT,PFTIME,KNTTM,KNTUVP
      INTEGER LVLGEO,LVLCOF,LVLMGD,ADVSCM
      INTEGER KNTOUT,KNTIN,PF,I_MAX,J_MAX,K_MAX
      INTEGER IBW,IEW,JBW,JEW,KBW,KEW,ISOTHM
*
      REAL*8 GRDX1,GRDX2,GRDX3,GRDX,GRDY1,GRDY2,GRDY3,GRDY
      REAL*8 GRDZ1,GRDZ2,GRDZ3,GRDZ
      REAL*8 IRATX1,IRATX2,IRATX3,IRATX4
      REAL*8 IRATY1,IRATY2,IRATY3,IRATY4
      REAL*8 IRATZ1,IRATZ2,IRATZ3,IRATZ4
      REAL*8 RHO,COND,CP,VISC,BETA
      REAL*8 ITIME,TTIME,DTIME
      REAL*8 WORK3(N,N),WORK4(N),WORK5(NT,NT),WORK6(NT)
      CHARACTER*8 CHTIME
      CHARACTER*5 FNAPP     
*
*	  Variables changed for updated temperature
*      REAL*8 RHO(ID,JD,KD),COND(ID,JD,KD),CP(ID,JD,KD),VISC(ID,JD,KD),BETA(ID,JD,KD)
*
      INTEGER POROUS,PCONV,IPB,IPE,JPB,JPE,KPB,KPE
      INTEGER IPB2,IPE2,JPB2,JPE2,KPB2,KPE2
      INTEGER SOLID,ISB,ISE,JSB,JSE,KSB,KSE
      INTEGER ISB2,ISE2,JSB2,JSE2,KSB2,KSE2
      REAL*8 EPS,PD,LD,VISCB,FORCH,PERM
      REAL*8 CONDFE,CONDSX,CONDSY,CONDSZ,CSP,RHOSP,ASF
      REAL*8 DEPSDZ,KS
      REAL*8 CONDS,CS,RHOS
*
*	  Variables changed for updated temperature 
*      REAL*8 VISCB(ID,JD,KD),CONDFE(ID,JD,KD)
*
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
      REAL*8 MUEFF(ID,JD,KD),PRSTY(ID,JD,KD),KEFF(ID,JD,KD)
      REAL*8 KEFFX(ID,JD,KD),KEFFY(ID,JD,KD),KEFFZ(ID,JD,KD)
      REAL*8 DENST(ID,JD,KD),CFLUID(ID,JD,KD),CSOLID(ID,JD,KD)
      REAL*8 DEQ(ID,JD,KD),CFORCH(ID,JD,KD),KPERM(ID,JD,KD)
      REAL*8 SPECSA(ID,JD,KD),KEFF0(ID,JD,KD),SLDTY(ID,JD,KD)
      REAL*8 PEAST(ID,JD,KD),PNORTH(ID,JD,KD),PTOP(ID,JD,KD)
      REAL*8 UNSTF(ID,JD,KD)
      REAL*8 CONTF(ID,JD,KD,NNB+1),DIFTF(ID,JD,KD,NNB+1)
      REAL*8 UNSTS(ID,JD,KD)
      REAL*8 CONTS(ID,JD,KD,NNB+1),DIFTS(ID,JD,KD,NNB+1)
      REAL*8 UNSU(ID,JD,KD)
      REAL*8 CONU(ID,JD,KD,NNB+1),DIFU(ID,JD,KD,NNB+1)   
      REAL*8 HSF(ID,JD,KD)         
      REAL*8 SGENT,SGEN(ID,JD,KD),SGENF(ID,JD,KD),SGENS(ID,JD,KD)
      REAL*8 RATIOT,RATIOS,RATIO(ID,JD,KD)
*
*     New variables
*
      REAL*8 HEATERTIME,FANTIME  !Time steps for heater and fan on/off
      INTEGER HEATER
*
*============================
*  Initialization and input
*============================
*
      HEATERTIME = 2400.0 !Heater off time 
      FANTIME = 1800.0    !Fan on time
      HEATER = 18         ! 60:14 80:18
*
*--Read input parameters
*
      CALL FLDFIN(IDATI,IRSI,ITERMO)
      CALL INPUT(IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
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
      CALL FLDFOT(IRSO, FNAPP)
*
*--Create grid
*
      CALL CARGRD(XNET,YNET,ZNET, GRDX1,GRDX2,GRDX3,GRDX,
     C            GRDY1,GRDY2,GRDY3,GRDY,GRDZ1,GRDZ2,GRDZ3,GRDZ,
     C            IRATX1,IRATX2,IRATX3,IRATX4,
     C            IRATY1,IRATY2,IRATY3,IRATY4,
     C            IRATZ1,IRATZ2,IRATZ3,IRATZ4,
     C            IBW,IEW,JBW,JEW,KBW,KEW,X_IB,X_IE,Y_JB,Y_JE,
     C            Z_KB,Z_KE,XP,YP,ZP,XE,YN,ZT,DIEP,DJNP,DKTP,
     C            DISE,DISN,DIST,AREP,ARNP,ARTP,VOLP,LVLGEO,
     C            POROUS,IPB,IPE,JPB,JPE,KPB,KPE,
     C            IPB2,IPE2,JPB2,JPE2,KPB2,KPE2,
     C            SOLID,ISB,ISE,JSB,JSE,KSB,KSE,
     C            ISB2,ISE2,JSB2,JSE2,KSB2,KSE2,
     C            CVTYPE,IB,IE1,IE2,IE3,IE,JB,JE1,JE2,JE3,JE,
     C            KB,KE1,KE2,KE3,KE,ID,JD,KD,NNB)
*
*--Initialize field variables
*
      CALL INITAL(TF,TS,P,U,V,W,UHE,VHN,WHT,
     C            T0,P0,U0,V0,W0,
     C            CVTYPE,IRSTRT,IRSI,
     C            IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*--Initialize residuals (for RATE calculations)
*
      DO 9 I=1,6
        ARSD(I)=   0.0
	    RSDOLD(I)= 1.0E+10
	    NORM(I)=   1.0E+10
   9  CONTINUE
*
*--Create needed property arrays. Note for POROUS = 1
*  LD is used in DEQ array
*
      IF(POROUS.EQ.2) THEN
       CALL CFPPTY(PRSTY,SLDTY,DEQ,CFORCH,KPERM,MUEFF,
     C             KEFF0,KEFFX,KEFFY,KEFFZ,SPECSA, 
     C             EPS,PD,DEPSDZ,VISC,VISCB,
     C             PERM,FORCH,KS,COND,CONDS,ZP,GRDZ,PI,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      ELSE
       CALL PPTYS(PRSTY, DBL1,EPS,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)   
       CALL PPTYS(SLDTY, DBL0,(DBL1-EPS),DBL1, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)     
       CALL PPTYS(KEFF0, COND,CONDFE,CONDS, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(KEFFX, COND,CONDSX,CONDS, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(KEFFY, COND,CONDSY,CONDS, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(KEFFZ, COND,CONDSZ,CONDS, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(DEQ, DBL0,LD,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(CFORCH, DBL0,FORCH,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(KPERM, DBL0,PERM,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL PPTYS(SPECSA, DBL0,ASF,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)    
       CALL PPTYS(MUEFF, VISC,VISCB/EPS,DBL0, 
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      ENDIF
*
      CALL PPTYS(CFLUID, CP,CP,CS, 
     C           CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)            
      CALL PPTYS(CSOLID, CP,CSP,CS, 
     C           CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)   
*
*--Create arrays of factors for unsteady, advective, and diffusive
*  terms
*      
      CALL FACTR(UNSTF,CONTF,DIFTF,
     C           UNSTS,CONTS,DIFTS,
     C           UNSU,CONU,DIFU,     
     C           PRSTY,CP,CSP,CS,MUEFF,VISCB,
     C           DISE,DIEP,DISN,DJNP,DIST,DKTP,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)     
*
*--Compute coefficients in conservation of mass equation (constant
*  for fixed-grid, incompressible flow)
*
      CALL COEFCN(ACUW,ACUE,ACVS,ACVN,ACWB,ACWT,BC,
     C            RHO,AREP,ARNP,ARTP,CVTYPE,
     C            IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*--Compute pressure gradient field from initial conditions
*
      CALL PINT(PEAST,PNORTH,PTOP, 
     C          P,DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C          ME,MN,MT,U,V,W,PRSTY,DE,DN,DT,MUEFF,
     C          AREP,ARNP,ARTP,CONU,BLEND,ALFAE,ALFAN,ALFAT,
     C          0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      CALL GRADP(DPDX,DPDXP,DPDXM,DPDY,DPDYP,DPDYM,
     C           DPDZ,DPDZP,DPDZM, P,PEAST,PNORTH,PTOP,
     C           DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*=================================
*  Begin outer, time loop
*=================================
*
      PF=0
      TTIME= ITIME
*
      DO 2000 KNTOUT=1,KNTTM
*
      PF=PF+1
      TTIME= TTIME+DTIME
*
*--Transfer variables to old values
*
      DO 12 K=KB-1,KE+1
       DO 11 J=JB-1,JE+1
        DO 10 I=IB-1,IE+1
         TFOLD(I,J,K)= TF(I,J,K)
         TSOLD(I,J,K)= TS(I,J,K)
         POLD(I,J,K)=  P(I,J,K)
	     UOLD(I,J,K)=  U(I,J,K)
	     VOLD(I,J,K)=  V(I,J,K)
         WOLD(I,J,K)=  W(I,J,K)
 10     CONTINUE
 11    CONTINUE
 12   CONTINUE
*
*----------------------------------
*  Begin inner, linearization loop
*----------------------------------
*     
      DO 1700 KNTIN=1,KNTUVP
*
      WRITE(ITERMO,7220) KNTOUT,KNTIN,TTIME
*
*--Compute active coefficients for T, wherein porous volumes
*  the first energy equation is for the fluid phase and the 
*  second is for the solid phase
*
      IF(ISOTHM.GT.0) THEN
        CALL NULLMM(BT, IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATP,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATW,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATE,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATS,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATN,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)        
        CALL NULLMN(ATB,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
        CALL NULLMN(ATT,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*     Update material properties with old fluid temperature field
*
*        CALL INTERP(RHO, 1,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE)
*        CALL INTERP(CP, 2,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE) 
*        CALL INTERP(MU, 3,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE) 
*        CALL INTERP(COND, 4,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE) 
*        CALL INTERP(CONDFE, 4,TF,ID,JD,KD,IE,IB,JB,JE,KB,KE) 
*        !Update BETA based on temperature 
*        DO 9000 I=IB,IE
*        DO 9001 J=JB,JE
*        DO 9002 K=KB,KE 
*         BETA(I,J,K) = 1/TF(I,J,K)
* 9002   CONTINUE    
* 9001   CONTINUE
* 9000   CONTINUE
*
*--Compute coefficients for fluid phase energy equation
*
        CALL CONV(HSF, U,V,W,PRSTY,PD,DEQ,
     C            RHO,CP,COND,VISC,
     C            PCONV,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL PPTYS(DENST, RHO,RHO,RHOS, 
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
        CALL ADDISP(KEFF, KEFF0,U,V,W,PRSTY,KPERM,PD,LD,
     C              RHO,CP,COND,VISC,
     C              CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL DIFPHI(DEF,DNF,DTF, KEFF,KEFF,KEFF,AREP,ARNP,ARTP,
     C              DIEP,DJNP,DKTP,DISE,DISN,DIST,PRSTY,
     C              0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL MASFLX(ME,MN,MT, UHE,VHN,WHT,ACUW,ACUE,ACVS,ACVN,
     C              ACWB,ACWT,IB,IE,JB,JE,KB,KE,ID,JD,KD)
        CALL WEIGHT(ALFAE,ALFAN,ALFAT, ME,MN,MT,
     C              DEF,DNF,DTF,CONTF,DIFTF,
     C              CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL HOCONV(DCCE,DCCN,DCCT, ALFAE,ALFAN,ALFAT,ME,MN,MT,
     C              DEF,DNF,DTF,KEFF,XP,XNET,YP,YNET,ZP,ZNET,
     C              DISE,DIEP,DISN,DJNP,DIST,DKTP,AREP,ARNP,ARTP,
     C              U,UHE,V,VHN,W,WHT,TF,CONTF,DIFTF,PRSTY,0,
     C              CVTYPE,ADVSCM,BLEND,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL SRCTF(QTF,RTF, DCCE,DCCN,DCCT,HSF,SPECSA,VOLP,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C             TTIME,
     C             HEATERTIME,HEATER)
        CALL COEFFM(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C              ME,MN,MT,DEF,DNF,DTF,QTF,RTF,TFOLD,TF,VOLP,
     C              ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C              DENST,DTIME,UNSTF,CONTF,DIFTF,
     C              1,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB)
*        
*--Compute coefficients for solid phase energy equation
*
        CALL PPTYS(DENST, RHO,RHOSP,RHOS,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)           
        CALL DIFPHI(DES,DNS,DTS, KEFFX,KEFFY,KEFFZ,AREP,ARNP,ARTP,
     C              DIEP,DJNP,DKTP,DISE,DISN,DIST,SLDTY,
     C              0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL SRCTS(QTS,RTS, HSF,SPECSA,VOLP,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C             TTIME,
     C             HEATERTIME,HEATER)
        CALL COEFFM(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C              ME,MN,MT,DES,DNS,DTS,QTS,RTS,TSOLD,TS,VOLP,
     C              ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C              DENST,DTIME,UNSTS,CONTS,DIFTS,
     C              2,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB)
*
*--Add cross-terms in energy equations
*
        CALL CRST(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C            HSF,SPECSA,VOLP,DEF,DNF,DTF,DES,DNS,DTS,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB) 
*
*--Set boundary conditions in T equation
*
        CALL BNDCT(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,TF,TS,
     C             XP,YP,ZP,DIEP,DJNP,DKTP,AREP,ARNP,ARTP,
     C             DEF,DNF,DTF,DES,DNS,DTS,DTMX,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB,
     C             HEATER)
      ENDIF
*
*--Compute active coefficients for U & V momentum
*
        CALL NULLMM(BU, IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUP,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUW,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUS,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUN,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUB,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
        CALL NULLMN(AUT,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
        CALL PPTYS(DENST, RHO,RHO,RHOS, 
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL DIFPHI(DE,DN,DT, MUEFF,MUEFF,MUEFF,AREP,ARNP,ARTP,
     C              DIEP,DJNP,DKTP,DISE,DISN,DIST,PRSTY,
     C              1,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)   
        CALL MASFLX(ME,MN,MT, UHE,VHN,WHT,ACUW,ACUE,ACVS,ACVN,
     C              ACWB,ACWT,IB,IE,JB,JE,KB,KE,ID,JD,KD)
        CALL WEIGHT(ALFAE,ALFAN,ALFAT, ME,MN,MT,
     C              DE,DN,DT,CONU,DIFU,
     C              CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL HOCONV(DCCE,DCCN,DCCT, ALFAE,ALFAN,ALFAT,ME,MN,MT,
     C              DE,DN,DT,MUEFF,XP,XNET,YP,YNET,ZP,ZNET,
     C              DISE,DIEP,DISN,DJNP,DIST,DKTP,AREP,ARNP,ARTP,
     C              U,UHE,V,VHN,W,WHT,U,CONU,DIFU,PRSTY,1,
     C              CVTYPE,ADVSCM,BLEND,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL SRCU(QU,RU, DCCE,DCCN,DCCT,VOLP,
     C            VISC,KPERM,RHO,CFORCH,PRSTY,U,V,W,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
        CALL HOCONV(DCCE,DCCN,DCCT, ALFAE,ALFAN,ALFAT,ME,MN,MT,
     C              DE,DN,DT,MUEFF,XP,XNET,YP,YNET,ZP,ZNET,
     C              DISE,DIEP,DISN,DJNP,DIST,DKTP,AREP,ARNP,ARTP,
     C              U,UHE,V,VHN,W,WHT,V,CONU,DIFU,PRSTY,1,
     C              CVTYPE,ADVSCM,BLEND,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)    
        CALL SRCV(QV,RV, DCCE,DCCN,DCCT,VOLP,
     C            VISC,KPERM,RHO,CFORCH,PRSTY,U,V,W,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,
     C            TF,BETA)
        CALL HOCONV(DCCE,DCCN,DCCT, ALFAE,ALFAN,ALFAT,ME,MN,MT,
     C              DE,DN,DT,MUEFF,XP,XNET,YP,YNET,ZP,ZNET,
     C              DISE,DIEP,DISN,DJNP,DIST,DKTP,AREP,ARNP,ARTP,
     C              U,UHE,V,VHN,W,WHT,W,CONU,DIFU,PRSTY,1,
     C              CVTYPE,ADVSCM,BLEND,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)    
        CALL SRCW(QW,RW, DCCE,DCCN,DCCT,
     C            TF,RHO,VOLP,GEE,BETA,
     C            VISC,KPERM,CFORCH,PRSTY,U,V,W,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)    
        CALL COEFFM(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              ME,MN,MT,DE,DN,DT,QU,RU,UOLD,U,VOLP,
     C              ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C              DENST,DTIME,UNSU,CONU,DIFU,
     C              2,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
        CALL COEFFM(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              ME,MN,MT,DE,DN,DT,QV,RV,VOLD,V,VOLP,
     C              ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C              DENST,DTIME,UNSU,CONU,DIFU,
     C              3,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
        CALL COEFFM(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              ME,MN,MT,DE,DN,DT,QW,RW,WOLD,W,VOLP,
     C              ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C              DENST,DTIME,UNSU,CONU,DIFU,
     C              4,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)      
*
*--Add the pressure terms to the momentum equations
*
        CALL SRCUVP(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              P,PEAST,PNORTH,PTOP,PRSTY,   
     C              VOLP,DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C              CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*--Set the boundary conditions in the momentum equations
*
        CALL BNDCU(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C             XP,YP,ZP,GRDZ,FORCH,PERM,EPS,VISC,RHO,UIN,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,U,
     C             TTIME,FANTIME)
        CALL BNDCV(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C             CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,
     C             UIN,TTIME,FANTIME)
        CALL BNDCW(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C             CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,
     C             TTIME,FANTIME)
*
*--Compute terms required to couple mass-momentum
*
        CALL DHAT(DHUE,DHVN,DHWT,VOLE,VOLN,VOLT,
     C            AEAST,ANORTH,ATOP, AUP,VOLP,PRSTY,RHO,
     C            DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*--Insert the mass coefficients into the coefficient blocks
*
        CALL ADCONT(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              ACUW,ACUE,ACVS,ACVN,ACWB,ACWT,
     C              AEAST,ANORTH,ATOP,DHUE,DHVN,DHWT,
     C              DPDXP,DPDXM,DPDYP,DPDYM,DPDZP,DPDZM,
     C              P,PEAST,PNORTH,PTOP,RHO,
     C              VOLP,VOLE,VOLN,VOLT,PRSTY,
     C              DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C              CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*--Set the boundary conditions in the P (mass) equation
*         
        CALL BNDCP(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C             P,DIEP,DJNP,DKTP,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,
     C             RHO,UIN)
*
*--Compute and report normalized residuals for energy, mass and momentum
*
        CALL UVPTRS(TF,TS,P,U,V,W,
     C              NORM,RSD,ARSD,RSDMAX,RSDOLD,I_MAX,J_MAX,K_MAX,
     C              AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C              ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C              CFLUID,CSOLID,CVTYPE,ITERMO,ISOTHM,
     C              IB,IE,JB,JE,KB,KE,N,NT,ID,JD,KD,NNB)
*
*--Check residuals of all field variables against convergence criterion
*
        IF(KNTIN.GT.1) THEN
         MAXRSD= ABS(RSDMAX(1)/NORM(1))
         DO 1710 I=2,6
	    IF(ABS(RSDMAX(I)/NORM(I)).GT.MAXRSD) 
     C                   MAXRSD=ABS(RSDMAX(I)/NORM(I))
 1710    CONTINUE
*
	 IF(MAXRSD.LT.CRIT) THEN
          WRITE(ITERMO,7216) KNTOUT
	  GOTO 1800
	 ELSEIF(MAXRSD.GT.1.0E+04) THEN
          WRITE(ITERMO,7215) KNTOUT
	  GOTO 9999
	 ENDIF
        ENDIF
*
*--Store current residual levels
*
	IF(ISOTHM.GT.0) THEN
	  RSDOLD(1)=ARSD(1)/NORM(1)
	  RSDOLD(2)=ARSD(2)/NORM(2)
	ENDIF
        DO 1720 I=3,6
          RSDOLD(I)=ARSD(I)/NORM(I)
 1720   CONTINUE
*
*--Solve for T using a direct solver
*     
        IF(ISOTHM.GT.0) THEN
         CALL SOLT(TF,TS,WORK5,WORK6,RLX,
     C             ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C             IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
	ENDIF
*
*--Solve P,U,V,W as a coupled set using a direct solver
* 
        CALL SOLUVP(P,U,V,W,UHE,VHN,WHT,RLX,
     C              AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,WORK3,WORK4,
     C              IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*--Update pressure gradient fieldand face velocities for 
*  the mass fluxes
*
      DO 1600 L=1,ITER
       CALL PINT(PEAST,PNORTH,PTOP, 
     C           P,DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C           ME,MN,MT,U,V,W,PRSTY,DE,DN,DT,MUEFF,
     C           AREP,ARNP,ARTP,CONU,BLEND,ALFAE,ALFAN,ALFAT,
     C           1,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
       CALL GRADP(DPDX,DPDXP,DPDXM,DPDY,DPDYP,DPDYM,
     C            DPDZ,DPDZP,DPDZM, P,PEAST,PNORTH,PTOP,
     C            DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL UHAT(UHE, U,P,PEAST,DHUE,VOLE,VOLP,
     C           DPDXP,DPDXM,AEAST,AUP,PRSTY,
     C           RHO,DIEP,DISE,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
       CALL VHAT(VHN, V,P,PNORTH,DHVN,VOLN,VOLP,
     C           DPDYP,DPDYM,ANORTH,AUP,PRSTY,
     C           RHO,DJNP,DISN,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
       CALL WHAT(WHT, W,P,PTOP,DHWT,VOLT,VOLP,
     C           DPDZP,DPDZM,ATOP,AUP,PRSTY,
     C           RHO,DKTP,DIST,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)     
       CALL MASFLX(ME,MN,MT, UHE,VHN,WHT,ACUW,ACUE,ACVS,ACVN,
     C             ACWB,ACWT,IB,IE,JB,JE,KB,KE,ID,JD,KD)
       CALL WEIGHT(ALFAE,ALFAN,ALFAT, ME,MN,MT,
     C             DE,DN,DT,CONU,DIFU,
     C             CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)     
 1600 CONTINUE   
*
*     Calculate heat transfer characteristics
*
*       CALL HTCHAR(ID,JD,KD,IB,IE,JB,JE,KB,KE,GRDX,YP,
*     C             TF,TS,TIN,
*     C             CONDFE,CONDSX)
*
*     Print to a TECplot formatted file tec.dat, python file py.txt
*
       CALL TECPLT(XP,YP,ZP,IB,IE,JB,JE,KB,KE,ID,JD,KD,
     C             U,V,W,P,TS,TF,KNTOUT)
       CALL PYPLOT(XP,YP,ZP,IB,IE,JB,JE,KB,KE,ID,JD,KD,
     C             U,V,W,P,TS,TF,KNTOUT)
       PRINT *, 'inside inner loop'
*
*----------------------------------
*  End inner, linearization loop
*----------------------------------
*
 1700 CONTINUE
*
*  Print solution for present time step
*
 1800 CONTINUE 
      PRINT *,'out of inner loop'
      IF((PF/PFTIME).EQ.1) THEN
       CALL SGENC(SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
     C            U,V,W,TF,TS,
     C            DE,DN,DT,DEF,DNF,DTF,DES,DNS,DTS,
     C            PRSTY,VISC,RHO,KPERM,CFORCH,MUEFF,
     C            KEFF,KEFFX,KEFFY,KEFFZ,HSF,SPECSA,
     C            DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C            AREP,ARNP,ARTP,VOLP,
     C            CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*       CALL SAVE(TF,TS,P,U,V,W,UHE,VHN,WHT,XP,YP,ZP,
*     C           SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
*     C           AREP,ME,CP,DIEP,DISE,
*     C           DTF,DTS,DNF,DNS,DEF,DES,DTMX,
*     C           FNAPP,IRSO,IB,IE,JB,JE,KB,KE,ID,JD,KD,
*     C           CVTYPE,KNTOUT,0,IPB,IPE,NNB)
       PF=0
      ENDIF
*
*  Check temporal convergence
*
      MAXTR=0.0
      DO 1820 I=1,6
       TRMAX(I)=0.0
 1820 CONTINUE
      IF(ISOTHM.GT.0) THEN
       CALL TREST(TRMAX, TF,TS,TFOLD,TSOLD,CVTYPE,
     C            1,1,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)  
      ENDIF
      CALL TRESU(TRMAX, P,U,V,W,POLD,UOLD,VOLD,WOLD,CVTYPE,
     C           1,1,0,0,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      WRITE(ITERMO,7230) TRMAX(1),TRMAX(2),TRMAX(3),
     C                   TRMAX(4),TRMAX(5),TRMAX(6)
      DO 1825 I=1,6
       IF(TRMAX(I).GT.MAXTR) THEN
        MAXTR = TRMAX(I)
       ENDIF
 1825 CONTINUE
      IF(MAXTR.LT.TCRIT .AND. MAXRSD.LT.CRIT) THEN
       GOTO 1840
      ENDIF
*
      IF(KNTOUT.GE.KNTTM) THEN
       WRITE(ITERMO,7213)      
       GOTO 1840
      ENDIF
*
*=================================
*  End outer, time loop
*=================================
*
 2000 CONTINUE
*
*--Write all field variables to an unformatted output file, rso.bin,
*  and to CVELO.dat for direct input to Python.
*
 1840 OPEN(UNIT=34,FILE='CVELO.dat')
        WRITE (34,*) '# XP YP ZP P U'
        B=2          !BJ: 81/101 PP: 11
        C=6           !BJ: 3 PP: 6
        DO 64 A=1,102  !BJ: 22/42 PP:102 
          WRITE (34,327)  XP(A),YP(B),ZP(C),P(A,B,C),U(A,B,C)   !BJ: XP(B),YP(A),ZP(C),P(B,A,C),U(B,A,C) PP: XP(A),YP(B),ZP(C),P(A,B,C),U(A,B,C)
  64    CONTINUE
  327  FORMAT(5(2X,1PE12.5))
        CALL FLUSH(34)
      CLOSE(UNIT=34)
*
       CALL SGENC(SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
     C           U,V,W,TF,TS,
     C           DE,DN,DT,DEF,DNF,DTF,DES,DNS,DTS,
     C           PRSTY,VISC,RHO,KPERM,CFORCH,MUEFF,
     C           KEFF,KEFFX,KEFFY,KEFFZ,HSF,SPECSA,
     C           DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C           AREP,ARNP,ARTP,VOLP,
     C           CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*      CALL SAVE(TF,TS,P,U,V,W,UHE,VHN,WHT,XP,YP,ZP,
*     C          SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
*     C          AREP,ME,CP,DIEP,DISE,
*     C          DTF,DTS,DNF,DNS,DEF,DES,DTMX,
*     C          FNAPP,IRSO,IB,IE,JB,JE,KB,KE,ID,JD,KD,
*     C          CVTYPE,KNTOUT,1,IPB,IPE,NNB)
*
*-----------------------------------------------------------------------
*  Format statements
*
 7213 FORMAT(/,'ITERATION LIMIT REACHED WITHOUT TIME CONVERGENCE')
 7214 FORMAT(//,'Time Step: ',I5,/)
 7215 FORMAT(/,'JOB STOPPED: MAXRSD EXCEEDS 1.0E+04 AT TIMESTEP ',I4,/)
 7216 FORMAT(/,'JOB CONVERGED AT TIMESTEP ',I4,/)
 7217 FORMAT(1PE13.5,2X,1PE13.5)
 7220 FORMAT(
     C '=======================================================',/
     C,'TIME STEP= 'I4'('I2')     SIMULATION TIME = '1PE13.5    ,/
     C,'-------------------------------------------------------',/  
     C,'|Equation  | Rate  | Ave RSD | Max RSD | Max Location |',/
     C,'+----------+-------+---------+---------+--------------+')
 7224 FORMAT(1PE11.4,3(1X,1PE11.4))
 7225 FORMAT(1PE11.4,4(1X,1PE11.4))
 7226 FORMAT(1PE11.4,5(1X,1PE11.4))
 7230 FORMAT(
     C '+--------------------------------+',/
     C,'| Maximum Temporal Residuals...  |',/
     C,'|    Tf: '1PE11.3,'             |',/
     C,'|    Ts: '1PE11.3,'             |',/
     C,'|    P:  '1PE11.3,'             |',/
     C,'|    U:  '1PE11.3,'             |',/
     C,'|    V:  '1PE11.3,'             |',/
     C,'|    W:  '1PE11.3,'             |',/
     C,'+--------------------------------+',/)
*
*--End of main
*
 9999 STOP
      END

