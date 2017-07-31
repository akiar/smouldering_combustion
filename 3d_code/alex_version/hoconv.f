*
*        file hoconv.f
************************************************************************
*
      SUBROUTINE HOCONV(DCCE,DCCN,DCCT, ALFAE,ALFAN,ALFAT,ME,MN,MT,
     C    DE,DN,DT,GAMA,XP,XNET,YP,YNET,ZP,ZNET,
     C    DISE,DIEP,DISN,DJNP,DIST,DKTP,AREP,ARNP,ARTP,
     C    U,UHE,V,VHN,W,WHT,PHI,CFACT,DFACT,PRSTY,
     C    EQTYPE,CVTYPE,ADVSCM,BLEND,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the deferred correction when using the
*     CDS, interface corrected CDS, or interface corrected TVD-MUSCL
*     to augment the implicit Peclet weighted UDS scheme. Even for UDS
*     scheme, interface corrections are employed.
*
*     Subroutine to calculate deferred correction terms.
*
*           ADVSCM=2: CDS, corrected for porous interface
*           ADVSCM=3: TVD-MUSCL, corrected for porous interface
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 DCCE(ID,JD,KD),DCCN(ID,JD,KD),DCCT(ID,JD,KD)
      REAL*8 ALFAE(ID,JD,KD),ALFAN(ID,JD,KD),ALFAT(ID,JD,KD)
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD),GAMA(ID,JD,KD)
      REAL*8 XP(ID),XNET(ID),YP(JD),YNET(JD),ZP(KD),ZNET(KD)
      REAL*8 DISE(ID),DIEP(ID),DISN(JD),DJNP(JD),DIST(KD),DKTP(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD)
      REAL*8 U(ID,JD,KD),UHE(ID,JD,KD)
      REAL*8 V(ID,JD,KD),VHN(ID,JD,KD)
      REAL*8 W(ID,JD,KD),WHT(ID,JD,KD)      
      REAL*8 PHI(ID,JD,KD),BLEND,PRSTY(ID,JD,KD)
      REAL*8 CFACT(ID,JD,KD,NNB+1),DFACT(ID,JD,KD,NNB+1)
      REAL*8 RSUBE,RSUBN,RSUBT,PSIE,PSIN,PSIT,CSUBE,CSUBN,CSUBT
      REAL*8 CFACTP,G1,G2,G3,UHOS,UUDS,VHOS,VUDS,WHOS,WUDS
      REAL*8 LMTR,PHIMAX,PHIMIN,RP(ID,JD,KD),DELTAP,DELTAN,ECONST
      REAL*8 DISXMA,DISXMI,DISYMA,DISYMI,DISZMA,DISZMI
      REAL*8 PHIECD,PHINCD,PHITCD
      REAL*8 PSIE2(ID,JD,KD),PSIN2(ID,JD,KD),PSIT2(ID,JD,KD),PSIMIN
      INTEGER I2,J2,K2
      INTEGER ADVSCM,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),EQTYPE,I,J,K,IEM1,JEM1,KEM1
      PARAMETER(LMTR=100.0)
*

      CALL NULLM(DCCE, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(DCCN, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(DCCT, IB,IE,JB,JE,KB,KE,ID,JD,KD)      
*
*==========================================
*  Interface corrections
*==========================================
*
      IEM1 = IE - 1
      JEM1 = JE - 1
      KEM1 = KE - 1        
      DO 10 I=IB,IEM1
       DO 5 J=JB,JE
        DO 4 K=KB,KE
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,3).NE.3
     C     .AND. CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,3)) THEN
*         
*         Account for difference in form between energy equation
*         correction and momentum equation correction as per thesis;
*         EQTYPE = 0 corresponds to energy equation. Weightings come
*         from harmonic mean diffusion coefficients.
*
          IF(EQTYPE.EQ.0 .AND. CVTYPE(I,J,K,1).EQ.0) THEN
           CSUBE = DISE(I)*DE(I,J,K)/(PRSTY(I+1,J,K)*GAMA(I,J,K)
     C             *AREP(J,K))
          ELSE
           CSUBE = DISE(I)*DE(I,J,K)/(GAMA(I,J,K)*AREP(J,K))
          ENDIF
*         
*         Determine porous advection factor
*
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
           CFACTP = CFACT(I,J,K,3)
          ELSE
           CFACTP = CFACT(I,J,K,1)
          ENDIF
*         
*         Calculate interface correction as per thesis
*                    
          DCCE(I,J,K) = BLEND*ME(I,J,K)*((CFACTP*CSUBE
     C      -0.5*CFACT(I,J,K,3)*(1.0-ALFAE(I,J,K)))*PHI(I+1,J,K)
     C      +(CFACTP*(1.0-CSUBE)-0.5*CFACT(I,J,K,1)*(1.0
     C      +ALFAE(I,J,K)))*PHI(I,J,K))
         ENDIF
    4   CONTINUE
    5  CONTINUE
   10 CONTINUE          
*           
      DO 20 I=IB,IE
       DO 15 J=JB,JEM1
        DO 14 K=KB,KE
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,5).NE.3
     C     .AND. CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,5)) THEN
*         
*         Account for difference in form between energy equation 
*         correction and momentum equation correction as per thesis;
*         EQTYPE = 0 corresponds to energy equation. Weightings come 
*         from harmonic mean diffusion coefficients.
*               
          IF(EQTYPE.EQ.0 .AND. CVTYPE(I,J,K,1).EQ.0) THEN
           CSUBN = DISN(J)*DN(I,J,K)/(PRSTY(I,J+1,K)*GAMA(I,J,K)
     C             *ARNP(I,K))
          ELSE
           CSUBN = DISN(J)*DN(I,J,K)/(GAMA(I,J,K)*ARNP(I,K))
          ENDIF
*         
*         Determine porous advection factor
*          
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
           CFACTP = CFACT(I,J,K,5)
          ELSE
           CFACTP = CFACT(I,J,K,1)
          ENDIF          
*         
*         Calculate interface correction as per thesis
*          
          DCCN(I,J,K) = BLEND*MN(I,J,K)*((CFACTP*CSUBN
     C      -0.5*CFACT(I,J,K,5)*(1.0-ALFAN(I,J,K)))*PHI(I,J+1,K)
     C      +(CFACTP*(1.0-CSUBN)-0.5*CFACT(I,J,K,1)*(1.0
     C      +ALFAN(I,J,K)))*PHI(I,J,K))
         ENDIF
   14   CONTINUE
   15  CONTINUE
   20 CONTINUE
*           
      DO 30 I=IB,IE
       DO 25 J=JB,JE
        DO 24 K=KB,KEM1
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,7).NE.3
     C     .AND. CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,7)) THEN
*         
*         Account for difference in form between energy equation 
*         correction and momentum equation correction as per thesis;
*         EQTYPE = 0 corresponds to energy equation. Weightings come 
*         from harmonic mean diffusion coefficients.
*               
          IF(EQTYPE.EQ.0 .AND. CVTYPE(I,J,K,1).EQ.0) THEN
           CSUBT = DIST(K)*DT(I,J,K)/(PRSTY(I,J,K+1)*GAMA(I,J,K)
     C             *ARTP(I,J))
          ELSE
           CSUBT = DIST(K)*DT(I,J,K)/(GAMA(I,J,K)*ARTP(I,J))
          ENDIF
*         
*         Determine porous advection factor
*          
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
           CFACTP = CFACT(I,J,K,7)
          ELSE
           CFACTP = CFACT(I,J,K,1)
          ENDIF          
*         
*         Calculate interface correction as per thesis
*                    
          DCCT(I,J,K) = BLEND*MT(I,J,K)*((CFACTP*CSUBT
     C      -0.5*CFACT(I,J,K,7)*(1.0-ALFAT(I,J,K)))*PHI(I,J,K+1)
     C      +(CFACTP*(1.0-CSUBT)-0.5*CFACT(I,J,K,1)*(1.0
     C      +ALFAT(I,J,K)))*PHI(I,J,K))
         ENDIF
   24   CONTINUE
   25  CONTINUE
   30 CONTINUE
*
*==========================================
*  ADVSCM=2, CDS
*==========================================
*
      IF(ADVSCM.EQ.2) THEN
*       
        DO 40 I=IB,IEM1
         DO 35 J=JB,JE
          DO 34 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,3)) THEN
            DCCE(I,J,K) = BLEND*ME(I,J,K)*((DISE(I)/DIEP(I)
     C                    -0.5*(1.0-ALFAE(I,J,K)))
     C                    *CFACT(I,J,K,3)*PHI(I+1,J,K)
     C                    +(DISE(I+1)/DIEP(I)
     C                    -0.5*(1.0+ALFAE(I,J,K)))
     C                    *CFACT(I,J,K,1)*PHI(I,J,K))
           ENDIF
   34     CONTINUE
   35    CONTINUE
   40   CONTINUE
*
        DO 50 I=IB,IE
         DO 45 J=JB,JEM1
          DO 44 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,5)) THEN
            DCCN(I,J,K) = BLEND*MN(I,J,K)*((DISN(J)/DJNP(J)
     C                    -0.5*(1.0-ALFAN(I,J,K)))
     C                    *CFACT(I,J,K,5)*PHI(I,J+1,K)
     C                    +(DISN(J+1)/DJNP(J)
     C                    -0.5*(1.0+ALFAN(I,J,K)))
     C                    *CFACT(I,J,K,1)*PHI(I,J,K))
           ENDIF
   44     CONTINUE
   45    CONTINUE
   50   CONTINUE
*
        DO 60 I=IB,IE
         DO 55 J=JB,JE
          DO 54 K=KB,KEM1
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,7)) THEN
            DCCT(I,J,K) = BLEND*MT(I,J,K)*((DIST(K)/DKTP(K)
     C                    -0.5*(1.0-ALFAT(I,J,K)))
     C                    *CFACT(I,J,K,7)*PHI(I,J,K+1)
     C                    +(DIST(K+1)/DKTP(K)
     C                    -0.5*(1.0+ALFAT(I,J,K)))
     C                    *CFACT(I,J,K,1)*PHI(I,J,K))           
           ENDIF
   54     CONTINUE
   55    CONTINUE
   60   CONTINUE      
*
*================================================
*  ADVSCM=3, TVD-MUSCL
*================================================
*
      ELSEIF(ADVSCM.EQ.3) THEN
*
        DO 70 I=IB,IEM1
         DO 65 J=JB,JE
          DO 64 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,3)) THEN
            IF(ME(I,J,K).LT.0.0) THEN
             RSUBE=(CFACT(I,J,K,3)*PHI(I+1,J,K)
     C        -CFACT(I+2,J,K,1)*PHI(I+2,J,K))/(CFACT(I,J,K,1)
     C        *PHI(I,J,K)-CFACT(I,J,K,3)*PHI(I+1,J,K)+1.0E-20)
     C        *(DIEP(I)/DIEP(I+1))
             PSIE = (RSUBE+ABS(RSUBE))/(1.0+ABS(RSUBE))
             DCCE(I,J,K)=BLEND*ME(I,J,K)*0.5*(1.0-PSIE+ALFAE(I,J,K))
     C                   *(CFACT(I,J,K,3)*PHI(I+1,J,K)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ELSE
             RSUBE=(CFACT(I,J,K,1)*PHI(I,J,K)
     C        -CFACT(I,J,K,2)*PHI(I-1,J,K))/(CFACT(I,J,K,3)
     C        *PHI(I+1,J,K)-CFACT(I,J,K,1)*PHI(I,J,K)+1.0E-20)
     C        *(DIEP(I)/DIEP(I-1))
             PSIE = (RSUBE+ABS(RSUBE))/(1.0+ABS(RSUBE))
             DCCE(I,J,K)=BLEND*ME(I,J,K)*0.5*(PSIE-1.0+ALFAE(I,J,K))
     C                   *(CFACT(I,J,K,3)*PHI(I+1,J,K)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ENDIF
           ENDIF
   64     CONTINUE
   65    CONTINUE
   70   CONTINUE
*           
        DO 80 I=IB,IE
         DO 75 J=JB,JEM1
          DO 74 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,5)) THEN
            IF(MN(I,J,K).LT.0.0) THEN
             RSUBN=(CFACT(I,J,K,5)*PHI(I,J+1,K)
     C        -CFACT(I,J+2,K,1)*PHI(I,J+2,K))/(CFACT(I,J,K,1)
     C        *PHI(I,J,K)-CFACT(I,J,K,5)*PHI(I,J+1,K)+1.0E-20)
     C        *(DJNP(J)/DJNP(J+1))
             PSIN = (RSUBN+ABS(RSUBN))/(1.0+ABS(RSUBN))
             DCCN(I,J,K)=BLEND*MN(I,J,K)*0.5*(1.0-PSIN+ALFAN(I,J,K))
     C                   *(CFACT(I,J,K,5)*PHI(I,J+1,K)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ELSE
             RSUBN=(CFACT(I,J,K,1)*PHI(I,J,K)
     C        -CFACT(I,J,K,4)*PHI(I,J-1,K))/(CFACT(I,J,K,5)
     C        *PHI(I,J+1,K)-CFACT(I,J,K,1)*PHI(I,J,K)+1.0E-20)
     C        *(DJNP(J)/DJNP(J-1))
             PSIN = (RSUBN+ABS(RSUBN))/(1.0+ABS(RSUBN))
             DCCN(I,J,K)=BLEND*MN(I,J,K)*0.5*(PSIN-1.0+ALFAN(I,J,K))
     C                   *(CFACT(I,J,K,5)*PHI(I,J+1,K)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ENDIF
           ENDIF
   74     CONTINUE
   75    CONTINUE
   80   CONTINUE
*
        DO 90 I=IB,IE
         DO 85 J=JB,JE
          DO 84 K=KB,KEM1
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,7)) THEN
            IF(MT(I,J,K).LT.0.0) THEN
             RSUBT=(CFACT(I,J,K,7)*PHI(I,J,K+1)
     C        -CFACT(I,J,K+2,1)*PHI(I,J,K+2))/(CFACT(I,J,K,1)
     C        *PHI(I,J,K)-CFACT(I,J,K,7)*PHI(I,J,K+1)+1.0E-20)
     C        *(DKTP(K)/DKTP(K+1))
             PSIT = (RSUBT+ABS(RSUBT))/(1.0+ABS(RSUBT))
             DCCT(I,J,K)=BLEND*MT(I,J,K)*0.5*(1.0-PSIT+ALFAT(I,J,K))
     C                   *(CFACT(I,J,K,7)*PHI(I,J,K+1)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ELSE
             RSUBT=(CFACT(I,J,K,1)*PHI(I,J,K)
     C        -CFACT(I,J,K,6)*PHI(I,J,K-1))/(CFACT(I,J,K,7)
     C        *PHI(I,J,K+1)-CFACT(I,J,K,1)*PHI(I,J,K)+1.0E-20)
     C        *(DKTP(K)/DKTP(K-1))
             PSIT = (RSUBT+ABS(RSUBT))/(1.0+ABS(RSUBT))
             DCCT(I,J,K)=BLEND*MT(I,J,K)*0.5*(PSIT-1.0+ALFAT(I,J,K))
     C                   *(CFACT(I,J,K,7)*PHI(I,J,K+1)
     C                   -CFACT(I,J,K,1)*PHI(I,J,K))
            ENDIF
           ENDIF
   84     CONTINUE
   85    CONTINUE
   90   CONTINUE
*
*================================================
*  ADVSCM=4, QUICK SCHEME
*================================================          
       ELSEIF (ADVSCM.EQ.4) THEN
*
        DO 110 I=IB,IEM1
         DO 105 J=JB,JE
          DO 104 K=KB,KE
          IF (ME(I,J,K).GE.0.0) THEN
          G1=((XNET(I)-XP(I))*(XNET(I)-XP(I+1)))/
     C       ((XP(I-1)-XP(I))*(XP(I-1)-XP(I+1)))
          G2=((XNET(I)-XP(I-1))*(XNET(I)-XP(I+1)))/
     C     ((XP(I)-XP(I-1))*(XP(I)-XP(I+1)))
          G3=((XNET(I)-XP(I-1))*(XNET(I)-XP(I)))/
     C       ((XP(I+1)-XP(I-1))*(XP(I+1)-XP(I)))
          UHOS=G1*CFACT(I,J,K,2)*PHI(I-1,J,K)+
     C         G2*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         G3*CFACT(I,J,K,3)*PHI(I+1,J,K)
          UUDS=0.5*((1.0+ALFAE(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAE(I,J,K))*CFACT(I,J,K,3)*PHI(I+1,J,K))
          DCCE(I,J,K)=BLEND*ME(I,J,K)*(UHOS-UUDS)
          ELSE
          G1=((XNET(I)-XP(I+1))*(XNET(I)-XP(I)))/
     C       ((XP(I+2)-XP(I+1))*(XP(I+2)-XP(I)))
          G2=((XNET(I)-XP(I+2))*(XNET(I)-XP(I)))/
     C       ((XP(I+1)-XP(I+2))*(XP(I+1)-XP(I)))
          G3=((XNET(I)-XP(I+2))*(XNET(I)-XP(I+1)))/
     C       ((XP(I)-XP(I+2))*(XP(I)-XP(I+1)))
          UHOS=G1*CFACT(I+2,J,K,1)*PHI(I+2,J,K)+
     C         G2*CFACT(I+1,J,K,1)*PHI(I+1,J,K)+
     C         G3*CFACT(I,J,K,1)*PHI(I,J,K)
          UUDS=0.5*((1.0+ALFAE(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAE(I,J,K))*CFACT(I,J,K,3)*PHI(I+1,J,K))
          DCCE(I,J,K)=BLEND*ME(I,J,K)*(UHOS-UUDS)
          ENDIF
  104     CONTINUE
  105     CONTINUE
  110     CONTINUE    
*
        DO 120 I=IB,IE
         DO 115 J=JB,JEM1
          DO 114 K=KB,KE 
          If (MN(I,J,K).GE.0.0) THEN
          G1=((YNET(J)-YP(J))*(YNET(J)-YP(J+1)))/
     C       ((YP(J-1)-YP(J))*(YP(J-1)-YP(J+1)))
          G2=((YNET(J)-YP(J-1))*(YNET(J)-YP(J+1)))/
     C       ((YP(J)-YP(J-1))*(YP(J)-YP(J+1)))
          G3=((YNET(J)-YP(J-1))*(YNET(J)-YP(J)))/
     C       ((YP(J+1)-YP(J-1))*(YP(J+1)-YP(J)))
          VHOS=G1*CFACT(I,J,K,4)*PHI(I,J-1,K)+
     C         G2*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         G3*CFACT(I,J,K,5)*PHI(I,J+1,K)
          VUDS=0.5*((1.0+ALFAN(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAN(I,J,K))*CFACT(I,J,K,5)*PHI(I,J+1,K))
          DCCN(I,J,K)=BLEND*MN(I,J,K)*(VHOS-VUDS)
          ELSE
          G1=((YNET(J)-YP(J+1))*(YNET(J)-YP(J)))/
     C       ((YP(J+2)-YP(J+1))*(YP(J+2)-YP(J)))
          G2=((YNET(J)-YP(J+2))*(YNET(J)-YP(J)))/
     C       ((YP(J+1)-YP(J+2))*(YP(J+1)-YP(J)))
          G3=((YNET(J)-YP(J+2))*(YNET(J)-YP(J+1)))/
     C       ((YP(J)-YP(J+2))*(YP(J)-YP(J+1)))
          VHOS=G1*CFACT(I,J+2,K,1)*PHI(I,J+2,K)+
     C         G2*CFACT(I,J+1,K,1)*PHI(I,J+1,K)+
     C         G3*CFACT(I,J,K,1)*PHI(I,J,K)
          VUDS=0.5*((1.0+ALFAN(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAN(I,J,K))*CFACT(I,J,K,5)*PHI(I,J+1,K))
          DCCN(I,J,K)=BLEND*MN(I,J,K)*(VHOS-VUDS)
          ENDIF
  114     CONTINUE
  115     CONTINUE
  120     CONTINUE
*
        DO 130 I=IB,IE
         DO 125 J=JB,JE
          DO 124 K=KB,KEM1
          IF (MT(I,J,K).GE.0.0) THEN
          G1=((ZNET(K)-ZP(K))*(ZNET(K)-ZP(K+1)))/
     C       ((ZP(K-1)-ZP(K))*(ZP(K-1)-ZP(K+1)))
          G2=((ZNET(K)-ZP(K-1))*(ZNET(K)-ZP(K+1)))/
     C       ((ZP(K)-ZP(K-1))*(ZP(K)-ZP(K+1)))
          G3=((ZNET(K)-ZP(K-1))*(ZNET(K)-ZP(K)))/
     C       ((ZP(K+1)-ZP(K-1))*(ZP(K+1)-ZP(K)))
          WHOS=G1*CFACT(I,J,K,6)*PHI(I,J,K-1)+
     C         G2*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         G3*CFACT(I,J,K,7)*PHI(I,J,K+1)
          WUDS=0.5*((1.0+ALFAT(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAT(I,J,K))*CFACT(I,J,K,7)*PHI(I,J,K+1))
          DCCN(I,J,K)=BLEND*MT(I,J,K)*(WHOS-WUDS)
          ELSE
          G1=((ZNET(K)-ZP(K+1))*(ZNET(K)-ZP(K)))/
     C       ((ZP(K+2)-ZP(K+1))*(ZP(K+2)-ZP(K)))
          G2=((ZNET(K)-ZP(K+2))*(ZNET(K)-ZP(K)))/
     C       ((ZP(K+1)-ZP(K+2))*(ZP(K+1)-ZP(K)))
          G3=((ZNET(K)-ZP(K+2))*(ZNET(K)-ZP(K+1)))/
     C       ((ZP(K)-ZP(K+2))*(ZP(K)-ZP(K+1)))
          WHOS=G1*CFACT(I,J,K+2,1)*PHI(I,J,K+2)+
     C         G2*CFACT(I,J,K+1,1)*PHI(I,J,K+1)+
     C         G3*CFACT(I,J,K,1)*PHI(I,J,K)
          WUDS=0.5*((1.0+ALFAT(I,J,K))*CFACT(I,J,K,1)*PHI(I,J,K)+
     C         (1-ALFAT(I,J,K))*CFACT(I,J,K,7)*PHI(I,J,K+1))
          DCCN(I,J,K)=BLEND*MT(I,J,K)*(WHOS-WUDS)
          ENDIF
  124     CONTINUE
  125     CONTINUE
  130     CONTINUE
*================================================
*  ADVSCM=5, TVD VENKATAKRISHNAN LIMITER
*================================================ 
      ELSEIF (ADVSCM.EQ.5) THEN
*      PRINT *,'TVD-VNK'
*      PRINT *,'PHI(2,2,2)',PHI(2,2,2)
*      PRINT *,'PHI(1,2,2)',PHI(1,2,2)
*      PRINT *,'PHI(3,2,2)',PHI(3,2,2)      
        DO 140 I=IB,IEM1
         DO 145 J=JB,JE
          DO 144 K=KB,KE        
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,3)) THEN
*               PRINT *,'I',I
*               PRINT *,'J',J
*               PRINT *,'K',K
*               PRINT *,'PHI(I,2,2)',PHI(I,J,K)
*               PRINT *,'PHI(I-1,2,2)',PHI(I-1,J,K)
*               PRINT *,'PHI(I+1,2,2)',PHI(I+1,J,K)
*               PRINT *,'DIEP(I)',DIEP(I)
*               PRINT *,'DIEP(I-1)',DIEP(I-1) 
*               PRINT *,'DIEP(I+1)',DIEP(I+1)
*               PRINT *,'DISE(I)',DISE(I)
*               PRINT *,'DISE(I-1)',DISE(I-1) 
*               PRINT *,'DISE(I+1)',DISE(I+1)               
               PHIMAX=PHI(I-1,J,K)
               PHIMIN=PHI(I-1,J,K)
               DISXMA=DIEP(I-1)
               DISXMI=DIEP(I-1)
*               PRINT *,'PHIMAX',PHIMAX
*               PRINT *,'PHIMIN',PHIMIN               
               DO 143 I2=I-1,I+1
               IF (PHI(I2,J,K).GT.PHIMAX) THEN
                 PHIMAX=PHI(I2,J,K)*CFACT(I2,J,K,1)
                 DISXMA=DIEP(I2)
               ENDIF
               IF (PHI(I2,J,K).LT.PHIMIN) THEN
                 PHIMIN=PHI(I2,J,K)*CFACT(I2,J,K,1)
                 DISXMI=DIEP(I2)
               ENDIF            
  143          CONTINUE  
*               PRINT *,'PHIMAX',PHIMAX
*               PRINT *,'PHIMIN',PHIMIN
*               PRINT *,'DISMI',DISXMA
*               PRINT *,'DISMA',DISXMA               
               PHIECD=(DISE(I)/DIEP(I))*PHI(I+1,J,K)*CFACT(I,J,K,3)+
     C            (DISE(I+1)/DIEP(I))*PHI(I,J,K)*CFACT(I,J,K,1)
*               PRINT *,'PHIECD',PHIECD
*               PRINT *,'PHI',PHI(I,J,K) 
             IF (PHIECD.GT.PHI(I,J,K)) THEN
*             PRINT *,'C1'
             DELTAP=(PHIMAX-PHI(I,J,K)*CFACT(I,J,K,1))/DISXMA
*             PRINT *,'DELTAP',DELTAP
             DELTAN=(PHIECD-PHI(I,J,K)*CFACT(I,J,K,1))/DISE(I)
*             PRINT *,'DELTAN',DELTAN             
             ECONST=(LMTR*2*DISE(I))**3
*             PRINT *,'ECONST',ECONST              
             PSIE2(I,J,K)=(1/(DELTAN+1E-20))*(((DELTAP**2+ECONST**2)
     C            *DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))
*             PRINT *,'PSIE',PSIE                           
             ELSEIF(PHIECD.LT.PHI(I,J,K)) THEN
*             PRINT *,'C2'             
             DELTAP=(PHIMIN-PHI(I,J,K)*CFACT(I,J,K,1))/DISXMI
*             PRINT *,'DELTAP',DELTAP
             DELTAN=(PHIECD-PHI(I,J,K)*CFACT(I,J,K,1))/DISE(I)
*             PRINT *,'DELTAN',DELTAN               
             ECONST=(LMTR*2*DISE(I))**3
*             PRINT *,'ECONST',ECONST             
             PSIE2(I,J,K)=(1/DELTAN+1E-20)*(((DELTAP**2+ECONST**2)
     C            *DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))                          
             ELSE
             PSIE2(I,J,K)=1.0
*             PRINT *,'PSIE=1'
             ENDIF
              IF (PSIE2(I,J,K).GT.1.0) THEN
              PSIE2(I,J,K)=1.0
              ENDIF          
           ENDIF
  144     CONTINUE
  145    CONTINUE
  140   CONTINUE
        DO 148 I=IB,IEM1
         DO 147 J=JB,JE
          DO 146 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C     .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,3)) THEN          
            PSIMIN=PSIE2(I-1,J,K)
            IF (PSIMIN.GT.PSIE2(I,J,K))THEN
            PSIMIN=PSIE2(I,J,K)
            PSIE2(I,J,K)=PSIMIN            
            ENDIF
              DCCE(I,J,K)=BLEND*ME(I,J,K)*((1-PSIE2(I,J,K)
     C                   *(DISE(I)/DIEP(I))-0.5*(1+ALFAE(I,J,K)))
     C                   *PHI(I,J,K)*CFACT(I,J,K,1) +
     C                   (PSIE2(I,J,K)*(DISE(I)/DIEP(I))-
     C                   0.5*(1-ALFAE(I,J,K)))*
     C                   PHI(I+1,J,K)*CFACT(I,J,K,3))
           ENDIF   
  146     CONTINUE
  147    CONTINUE
  148   CONTINUE
        DO 150 I=IB,IE
         DO 155 J=JB,JEM1
          DO 154 K=KB,KE        
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,5)) THEN
               PHIMAX=PHI(I,J-1,K)
               PHIMIN=PHI(I,J-1,K)
               DISYMA=DJNP(J-1)
               DISYMI=DJNP(J-1)               
               DO 153 J2=J-1,J+1
               IF (PHI(I,J2,K).GT.PHIMAX) THEN
                 PHIMAX=PHI(I,J2,K)*CFACT(I,J2,K,1)
                 DISYMA=DJNP(J2)
               ENDIF
               IF (PHI(I,J2,K).LT.PHIMIN) THEN
                 PHIMIN=PHI(I,J2,K)*CFACT(I,J2,K,1)
                 DISYMI=DJNP(J2)
               ENDIF
  153          CONTINUE  

            PHINCD=(DISN(J)/DJNP(J))*PHI(I,J+1,K)*CFACT(I,J,K,5)+
     C            (DISN(J+1)/DJNP(J))*PHI(I,J,K)*CFACT(I,J,K,1)

             IF (PHINCD.GT.PHI(I,J,K)) THEN
             DELTAP=(PHIMAX-PHI(I,J,K)*CFACT(I,J,K,1))/DISYMA
             DELTAN=(PHINCD-PHI(I,J,K)*CFACT(I,J,K,1))/DISN(J)
             ECONST=(LMTR*2*DISN(J))**3
             PSIN2(I,J,K)=(1/DELTAN+1E-20)*(((DELTAP**2+ECONST**2)
     C            *DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))
             ELSEIF(PHINCD.LT.PHI(I,J,K)) THEN
             DELTAP=(PHIMIN-PHI(I,J,K)*CFACT(I,J,K,1))/DISYMI
             DELTAN=(PHIECD-PHI(I,J,K)*CFACT(I,J,K,1))/DISN(J)
             ECONST=(LMTR*2*DISN(J))**3
             PSIN2(I,J,K)=(1/DELTAN+1E-20)*(((DELTAP**2+ECONST**2)*
     C            DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))
             ELSE
             PSIN2(I,J,K)=1.0
             ENDIF
              IF (PSIN2(I,J,K).GT.1.0) THEN
              PSIN2(I,J,K)=1.0
              ENDIF 
           ENDIF
  154     CONTINUE
  155    CONTINUE
  150   CONTINUE
        DO 158 I=IB,IE
         DO 157 J=JB,JEM1
          DO 156 K=KB,KE
           IF(CVTYPE(I,J,K,1).NE.3 
     C     .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,5)) THEN
            PSIMIN=PSIN2(I-1,J,K)
            IF (PSIMIN.GT.PSIN2(I,J,K))THEN
            PSIMIN=PSIN2(I,J,K)
            PSIN2(I,J,K)=PSIMIN            
            ENDIF          
              DCCN(I,J,K)=BLEND*MN(I,J,K)*((1-PSIN2(I,J,K)*
     C                   (DISN(J)/DJNP(J))-0.5*(1+ALFAN(I,J,K)))
     C                   *PHI(I,J,K)*CFACT(I,J,K,1) +
     C                   (PSIN2(I,J,K)*(DISN(J)/DJNP(J))
     C                   -0.5*(1-ALFAN(I,J,K)))*        
     C                   PHI(I,J+1,K)*CFACT(I,J,K,5))
           ENDIF  
  156     CONTINUE
  157    CONTINUE
  158   CONTINUE

        DO 160 I=IB,IE
         DO 165 J=JB,JE
          DO 164 K=KB,KEM1        
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,7)) THEN
               PHIMAX=PHI(I,J,K-1)
               PHIMIN=PHI(I,J,K-1)
               DISZMA=DKTP(K-1)                
               DISZMI=DKTP(K-1)               
               DO 163 K2=K-1,K+1
               IF (PHI(I,J,K2).GT.PHIMAX) THEN
                 PHIMAX=PHI(I,J,K2)*CFACT(I,J,K2,1)
                 DISZMA=DKTP(K2)
               ENDIF
               IF (PHI(I,J,K2).LT.PHIMIN) THEN
                 PHIMIN=PHI(I,J,K2)*CFACT(I,J,K2,1)
                 DISZMI=DKTP(K2)
               ENDIF
  163          CONTINUE  

            PHITCD=(DIST(K)/DKTP(K))*PHI(I,J,K+1)*CFACT(I,J,K,7)+
     C            (DIST(K+1)/DKTP(K))*PHI(I,J,K)*CFACT(I,J,K,1)

             IF (PHITCD.GT.PHI(I,J,K)) THEN
             DELTAP=(PHIMAX-PHI(I,J,K)*CFACT(I,J,K,1))/DISZMA
             DELTAN=(PHITCD-PHI(I,J,K)*CFACT(I,J,K,1))/DIST(K)
             ECONST=(LMTR*2*DIST(K))**3
             PSIT2(I,J,K)=(1/DELTAN+1E-20)*(((DELTAP**2+ECONST**2)
     C            *DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))
             ELSEIF(PHITCD.LT.PHI(I,J,K)) THEN
             DELTAP=(PHIMIN-PHI(I,J,K)*CFACT(I,J,K,1))/DISZMI
             DELTAN=(PHITCD-PHI(I,J,K)*CFACT(I,J,K,1))/DIST(K)
             ECONST=(LMTR*2*DISN(J))**3
             PSIT2(I,J,K)=(1/DELTAN+1E-20)*(((DELTAP**2+ECONST**2)
     C            *DELTAN+2*DELTAN**2*DELTAP) / (DELTAP**2+2*DELTAN**2
     C            +DELTAN*DELTAP+ECONST**2))
             ELSE
             PSIT2(I,J,K)=1.0
             ENDIF
              IF (PSIT2(I,J,K).GT.1.0) THEN
              PSIT2(I,J,K)=1.0
              ENDIF 
           ENDIF
  164     CONTINUE
  165    CONTINUE
  160   CONTINUE
        DO 168 I=IB,IE
         DO 167 J=JB,JE
          DO 166 K=KB,KEM1
           IF(CVTYPE(I,J,K,1).NE.3 
     C       .AND. CVTYPE(I,J,K,1).EQ.CVTYPE(I,J,K,7)) THEN          
            PSIMIN=PSIT2(I-1,J,K)
            IF (PSIMIN.GT.PSIT2(I,J,K))THEN
            PSIMIN=PSIT2(I,J,K)
            PSIT2(I,J,K)=PSIMIN            
            ENDIF          
             DCCT(I,J,K)=BLEND*MN(I,J,K)*((1-PSIT2(I,J,K)*
     C                   (DIST(K)/DKTP(K))
     C                   -0.5*(1+ALFAT(I,J,K)))
     C                   *PHI(I,J,K)*CFACT(I,J,K,1) +
     C                   (PSIT2(I,J,K)*(DIST(K)/DKTP(K))-
     C                   0.5*(1-ALFAT(I,J,K)))*
     C                   PHI(I,J,K+1)*CFACT(I,J,K,7))
           ENDIF
  166     CONTINUE
  167    CONTINUE
  168   CONTINUE             
      ENDIF
* 
      RETURN
      END

