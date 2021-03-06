*
*   This file contains 5 subroutines: BNDCT,BNDCP,BNDCU,BNDCV and BNDCW
*
************************************************************************
*
      SUBROUTINE BNDCT(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,TF,TS,
     C                 XP,YP,ZP,DIEP,DJNP,DKTP,AREP,ARNP,ARTP,
     C                 DEF,DNF,DTF,DES,DNS,DTS,DTMX,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB,
     C                 HEATER)
*
*     Subroutine to put the boundary condition information for each
*     P at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
* 
*            2) A boundary condition must be set for pressure in 
*            all flow domains.  In enclosures, a position inside the
*            domain must be selected.
*
*     L = 1 Fluid, L = 2 Solid
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD)
      REAL*8 ATS(NT,NT,ID,JD,KD),ATN(NT,NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 XP(ID),YP(JD),ZP(KD),DIEP(ID),DJNP(JD),DKTP(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD)
      REAL*8 DEF(ID,JD,KD),DNF(ID,JD,KD),DTF(ID,JD,KD)
      REAL*8 DES(ID,JD,KD),DNS(ID,JD,KD),DTS(ID,JD,KD)
      REAL*8 DTMX,QWALL,QTOT,XNODE,TIN
      PARAMETER(QWALL = 10.0)
      PARAMETER (TIN = 293.0)                                 !set inlet temperature
      INTEGER IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
      INTEGER HEATER
*
      IBM1 = IB - 1
      JBM1 = JB - 1
      KBM1 = KB - 1
      IEP1 = IE + 1
      JEP1 = JE + 1
      KEP1 = KE + 1
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
*     West face boundary conditions
*
        I = IBM1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*        
        ELSEIF(CVTYPE(I,J,K,3).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = TIN
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 1.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*                 
        ELSE
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 1.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 1.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
        ENDIF
*
*     East face boundary conditions
*
        I = IEP1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*
        ELSEIF(CVTYPE(I,J,K,2).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 1.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 1.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*
        ELSE
         L = 1
         ATW(L,L,I,J,K) = 1.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 1.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
        ENDIF
  5    CONTINUE  
 10   CONTINUE
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
*     South face boundary conditions
*
        J = JBM1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*
        ELSEIF(CVTYPE(I,J,K,5).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 1.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 1.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0

         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*
        ELSE
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0 !Changed T1 Nov. 7 from 1.0 to 0.0 
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0
         ATP(L,L,I,J,K) = 1.0
*
         BT(L,I,J,K) = TIN    !Changed T1 Nov. 7 from 0.0 TO TIN 
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0 
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0
         ATP(L,L,I,J,K) = 1.0
*
         BT(L,I,J,K) = TIN
*
        ENDIF
*
*     North face boundary conditions
*
        J = JEP1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*        
        ELSEIF(CVTYPE(I,J,K,4).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 1.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 1.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*
        ELSE
         L = 1    !Fluid
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 1.0 !Changed T1 Nov. 7 from 0.0 TO 1.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0 
         BT(L,I,J,K) = 0.0    !Changed T1 Nov. 7 from TIN TO 0.0 
*
         L = 2    !Solid
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 1.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0

         ATP(L,L,I,J,K) = 1.0 
         BT(L,I,J,K) = 0.0 
        ENDIF
 15    CONTINUE
 20   CONTINUE
* 
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
*     Bottom face boundary conditions
*
        K = KBM1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*        
        ELSEIF(CVTYPE(I,J,K,7).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 1.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 1.0 
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*                 
        ELSE
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 1.0         
         ATP(L,L,I,J,K) = 1.0
*         BT(L,I,J,K) = TIN+DTMX
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 1.0         
         ATP(L,L,I,J,K) = 1.0
*         BT(L,I,J,K) = TIN+DTMX
         BT(L,I,J,K) = 0.0

        ENDIF
*
*    Top face boundary conditions
*
        K = KEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 0.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Fluid CV adjacent to boundary
*        
        ELSEIF(CVTYPE(I,J,K,6).EQ.0) THEN
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 1.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 1.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
*       Porous CV adjacent to boundary
*                 
        ELSE
         L = 1
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 1.0
         ATT(L,L,I,J,K) = 0.0         
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
*
         L = 2
         ATW(L,L,I,J,K) = 0.0
         ATE(L,L,I,J,K) = 0.0
         ATS(L,L,I,J,K) = 0.0
         ATN(L,L,I,J,K) = 0.0
         ATB(L,L,I,J,K) = 1.0
         ATT(L,L,I,J,K) = 0.0 
         ATP(L,L,I,J,K) = 1.0
         BT(L,I,J,K) = 0.0
        ENDIF
 25    CONTINUE
 30   CONTINUE
*
*     Set temperature at heater 
*
*      J = HEATER
*      L = 1
*      DO 50 I = IB,IE
*       DO 55 K = KB,KE
*        ATW(L,L,I,J,K) = 0.0
*        ATE(L,L,I,J,K) = 0.0
*        ATS(L,L,I,J,K) = 0.0
*        ATN(L,L,I,J,K) = 0.0
*        ATB(L,L,I,J,K) = 0.0
*        ATT(L,L,I,J,K) = 0.0
*        ATP(L,L,I,J,K) = 1.0
*        BT(L,I,J,K) = 293.0
* 55    CONTINUE
* 50   CONTINUE
*      I = IE/2 + 1
*      J = JE/2 + 1
*      K = KE/2 + 1
*      DO 50 M = 1,4
*        AUW(L,M,I,J,K) = 0.0
*        AUE(L,M,I,J,K) = 0.0
*        AUS(L,M,I,J,K) = 0.0
*        AUN(L,M,I,J,K) = 0.0
*        AUB(L,M,I,J,K) = 0.0
*        AUT(L,M,I,J,K) = 0.0
*        AUP(L,M,I,J,K) = 0.0
*        BU(L,I,J,K) = 0.0
* 50   CONTINUE
*      AUP(L,L,I,J,K) = 1.0
*      BU(L,I,J,K) = 0.0
*
*     Conditions in solid and fluid CVs
*
      DO 45 I=IB,IE
       DO 40 J=JB,JE
        DO 35 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          L = 1
          ATW(L,L,I,J,K) = 0.0
          ATE(L,L,I,J,K) = 0.0
          ATS(L,L,I,J,K) = 0.0
          ATN(L,L,I,J,K) = 0.0
          ATB(L,L,I,J,K) = 0.0
          ATT(L,L,I,J,K) = 0.0         
          ATP(L,L,I,J,K) = 1.0
          BT(L,I,J,K) = 0.0
         ELSEIF(CVTYPE(I,J,K,1).EQ.0) THEN
          L = 2
          ATW(L,L,I,J,K) = 0.0
          ATE(L,L,I,J,K) = 0.0
          ATS(L,L,I,J,K) = 0.0
          ATN(L,L,I,J,K) = 0.0
          ATB(L,L,I,J,K) = 0.0
          ATT(L,L,I,J,K) = 0.0         
          ATP(L,L,I,J,K) = 1.0
          BT(L,I,J,K) = 300
         ENDIF
 35     CONTINUE        
 40    CONTINUE
 45   CONTINUE  
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCP(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C                 P,DIEP,DJNP,DKTP,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,
     C                 RHO,UIN)
*
*     Subroutine to put the boundary condition information for each
*     P at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
* 
*            2) A boundary condition must be set for pressure in 
*            all flow domains.  In enclosures, a position inside the
*            domain must be selected.
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD),AUB(N,N,ID,JD,KD)
      REAL*8 AUT(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 P(ID,JD,KD),DIEP(ID),DJNP(JD),DKTP(KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,N,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L,M,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
*
      REAL*8 RHO(ID,JD,KD),UIN
*
      L = 1
      IBM1 = IB - 1
      JBM1 = JB - 1
      KBM1 = KB - 1
      IEP1 = IE + 1
      JEP1 = JE + 1
      KEP1 = KE + 1
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
*     West face boundary conditions
*
        I = IBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 1.0     !BJ, PP : 1
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = DIEP(I)*(P(IB,J,K)-P(IB+1,J,K))/DIEP(IB) 
        ENDIF
*
*     East face boundary conditions
*
        I = IEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 1.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = DIEP(IE)*(P(IE,J,K)-P(IE-1,J,K))/DIEP(IE-1) 
        ENDIF        
  5    CONTINUE
 10   CONTINUE
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
*     South face boundary conditions
*
        J = JBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 1.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = DJNP(J)*(P(I,JB,K)-P(I,JB+1,K))/DJNP(JB)
        ENDIF
*
*     North face boundary conditions
*
        J = JEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 101325.0 !ATMS P0 in Marco !DJNP(JE)*(P(I,JE,K)-P(I,JE-1,K))/DJNP(JE-1)
        ENDIF
*
 15    CONTINUE
 20   CONTINUE
*
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
*     Bottom face boundary conditions
*
        K = KBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 1.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = DKTP(K)*(P(I,J,KB)-P(I,J,KB+1))/DKTP(KB)   !Changed T1 Nov. 7 from 0.0 TO INTERPOLATE
*         BU(L,I,J,K) = 0.0
        ENDIF
*
*     Top face boundary conditions
*
        K = KEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 1.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = DKTP(KE)*(P(I,J,KE)-P(I,J,KE-1))/DKTP(KE-1)    !Changed T1 Nov. 7 from 0.0 TO INTERPOLATE 
*         BU(L,I,J,K) = 0.0
        ENDIF
*
 25    CONTINUE
 30   CONTINUE
*
*     Conditions in solid CVs
*
      DO 45 I=IB,IE
       DO 40 J=JB,JE
        DO 35 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          AUW(L,L,I,J,K) = 0.0
          AUE(L,L,I,J,K) = 0.0
          AUS(L,L,I,J,K) = 0.0
          AUN(L,L,I,J,K) = 0.0
          AUB(L,L,I,J,K) = 0.0
          AUT(L,L,I,J,K) = 0.0         
          AUP(L,L,I,J,K) = 1.0
          BU(L,I,J,K) = 0.0
         ENDIF
 35     CONTINUE
 40    CONTINUE
 45   CONTINUE
*
*  Pressure level for some point in the domain (be careful!)
*
*      I = IE/2 + 1
*      J = JE/2 + 1
*      K = KE/2 + 1
*      DO 50 M = 1,4
*        AUW(L,M,I,J,K) = 0.0
*        AUE(L,M,I,J,K) = 0.0
*        AUS(L,M,I,J,K) = 0.0
*        AUN(L,M,I,J,K) = 0.0
*        AUB(L,M,I,J,K) = 0.0
*        AUT(L,M,I,J,K) = 0.0
*        AUP(L,M,I,J,K) = 0.0
*        BU(L,I,J,K) = 0.0
* 50   CONTINUE
*      AUP(L,L,I,J,K) = 1.0
*      BU(L,I,J,K) = 0.0
*     
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCU(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C                 XP,YP,ZP,GRDZ,FORCH,PERM,EPS,VISC,RHO,UIN,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,U,
     C                 TTIME,FANTIME)
*
*     Subroutine to put the boundary condition information for U
*     at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD),AUB(N,N,ID,JD,KD)
      REAL*8 AUT(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 XP(ID),YP(JD),ZP(KD),GRDZ,FORCH,PERM,EPS,UIN
      REAL*8 UINF,H,DA,LAMBDA,A,B,C,D,ZSTAR,CPRIME,CARG,X0,TOL,UCAP
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
      REAL*8 U(ID,JD,KD),FANTIME,TTIME
*
      REAL*8 VISC(ID,JD,KD),RHO(ID,JD,KD)
*
*      UINF = UIN
       H=0.02     !PP: 0.02 BJ: 0.01
*
      L = 2
      IBM1 = IB - 1
      JBM1 = JB - 1
      KBM1 = KB - 1      
      IEP1 = IE + 1
      JEP1 = JE + 1
      KEP1 = KE + 1    
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
*     West face boundary conditions
*
        I = IBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 1.0 !Changed T1 Nov. 7 from 0.0 TO 1.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K)= 0.0
*
        ENDIF
*
*     East face boundary conditions
*
        I = IEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 1.0 !Changed T1 Nov. 7 from 0.0 TO 1.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
  5    CONTINUE
 10   CONTINUE
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
*     South face boundary conditions
*
        J = JBM1
*
*       Solid CV adjacent to boundary
*
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*
        ELSE  !Changed T1 Nov. 7: REMOVED FAN CONDITIONS, SETTING U=0 
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
        ENDIF
*
*     North face boundary conditions
*
        J = JEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE  !Changed T1 Nov. 7: REMOVED FAN CONDITIONS, SETTING U=0 
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0 !Nov. 7 T1: Changed 1.0 to 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
        ENDIF        
*
 15    CONTINUE
 20   CONTINUE

*
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
*     Bottom face boundary conditions
*
        K = KBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
*     Top face boundary conditions
*
        K = KEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF        
*
 25    CONTINUE
 30   CONTINUE 
*
*     Conditions in solid CVs
*
      DO 45 I=IB,IE
       DO 40 J=JB,JE
        DO 35 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          AUW(L,L,I,J,K) = 0.0
          AUE(L,L,I,J,K) = 0.0
          AUS(L,L,I,J,K) = 0.0
          AUN(L,L,I,J,K) = 0.0
          AUB(L,L,I,J,K) = 0.0
          AUT(L,L,I,J,K) = 0.0         
          AUP(L,L,I,J,K) = 1.0
          BU(L,I,J,K) = 0.0
         ENDIF
 35     CONTINUE
 40    CONTINUE
 45   CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCV(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,UIN,
     C                 TTIME,FANTIME)
*
*     Subroutine to put the boundary condition information for V
*     at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD),AUB(N,N,ID,JD,KD)
      REAL*8 AUT(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 UIN
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
      REAL*8 TTIME,FANTIME
*
      L = 3
      IBM1 = IB - 1
      JBM1 = JB - 1
      KBM1 = KB - 1      
      IEP1 = IE + 1
      JEP1 = JE + 1
      KEP1 = KE + 1      
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
*     West face boundary conditions
*
        I = IBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 1.0 !Nov. 7 T1: Changed 0.0 to 1.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
*     East face boundary conditions
*
        I = IEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 1.0 !Nov. 7 T1: Changed 0.0 to 1.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
  5    CONTINUE
 10   CONTINUE
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
*     South face boundary conditions
*
        J = JBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         IF (TTIME <=FANTIME) THEN    !Modelling fan condition
          BU(L,I,J,K) = 0.0
         ELSE
          BU(L,I,J,K) = UIN
         ENDIF 
        ENDIF
*
*     North face boundary conditions
*
        J = JEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 1.0 !Nov. 7 T1: Removed fan condition and set to 1.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF        
*
 15    CONTINUE
 20   CONTINUE
*
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
*     Bottom face boundary conditions
*
        K = KBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 1.0 !Nov. 7 T1: Changed 0.0 to 1.0     
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
*     Top face boundary conditions
*
        K = KEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 1.0 !Nov. 7 T1: Changed 0.0 to 1.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
 25    CONTINUE
 30   CONTINUE 
*
*     Conditions in solid CVs
*
      DO 45 I=IB,IE
       DO 40 J=JB,JE
        DO 35 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          AUW(L,L,I,J,K) = 0.0
          AUE(L,L,I,J,K) = 0.0
          AUS(L,L,I,J,K) = 0.0
          AUN(L,L,I,J,K) = 0.0
          AUB(L,L,I,J,K) = 0.0
          AUT(L,L,I,J,K) = 0.0         
          AUP(L,L,I,J,K) = 1.0
          BU(L,I,J,K) = 0.0
         ENDIF
 35     CONTINUE
 40    CONTINUE
 45   CONTINUE
*
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE BNDCW(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU, 
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,
     C                 TTIME,FANTIME)
*
*     Subroutine to put the boundary condition information for W
*     at each boundary node into the finite difference coefficients.
*
*     Notes: 1) This routine is restricted to grids with a Cartesian
*            index layout. The outline should be filled in for each
*            boundary.
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD),AUB(N,N,ID,JD,KD)
      REAL*8 AUT(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L,IBM1,IEP1,JBM1,JEP1,KBM1,KEP1
      REAL*8 FANTIME,TTIME
*
      L = 4
      IBM1 = IB - 1
      JBM1 = JB - 1
      KBM1 = KB - 1      
      IEP1 = IE + 1
      JEP1 = JE + 1
      KEP1 = KE + 1      
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
*     West face boundary conditions
*
        I = IBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,3).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0  
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0 
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
*     East face boundary conditions
*
        I = IEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,2).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
  5    CONTINUE
 10   CONTINUE
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
*     South face boundary conditions
*
        J = JBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,5).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE  !Nov. 7 T1: Removed fan condition, set to 0.0
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0 
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*        
        ENDIF
*
*     North face boundary conditions
*
        J = JEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,4).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE  !Nov. 7 T1: Removed fan condition, set to 0.0
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0 
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF        
*
 15    CONTINUE
 20   CONTINUE
*
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
*     Bottom face boundary conditions
*
        K = KBM1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,7).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 1.0 !Nov. 7 T1: Changed 0.0 to 1.0
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF
*
*     Top face boundary conditions
*
        K = KEP1
*
*       Solid CV adjacent to boundary
*        
        IF(CVTYPE(I,J,K,6).EQ.3) THEN
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0
         AUT(L,L,I,J,K) = 0.0         
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
*
*       Non-solid CV adjacent to boundary
*        
        ELSE
         AUW(L,L,I,J,K) = 0.0
         AUE(L,L,I,J,K) = 0.0
         AUS(L,L,I,J,K) = 0.0
         AUN(L,L,I,J,K) = 0.0
         AUB(L,L,I,J,K) = 0.0 !Nov. 7 T1: Changed 0.0 to 1.0
         AUT(L,L,I,J,K) = 0.0
         AUP(L,L,I,J,K) = 1.0
         BU(L,I,J,K) = 0.0
        ENDIF        
*
 25    CONTINUE
 30   CONTINUE 
*
*     Conditions in solid CVs
*
      DO 45 I=IB,IE
       DO 40 J=JB,JE
        DO 35 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.3) THEN
          AUW(L,L,I,J,K) = 0.0
          AUE(L,L,I,J,K) = 0.0
          AUS(L,L,I,J,K) = 0.0
          AUN(L,L,I,J,K) = 0.0
          AUB(L,L,I,J,K) = 0.0
          AUT(L,L,I,J,K) = 0.0         
          AUP(L,L,I,J,K) = 1.0
          BU(L,I,J,K) = 0.0
         ENDIF
 35     CONTINUE
 40    CONTINUE
 45   CONTINUE
*
      RETURN
      END
