*     file crst.f
***********************************************************************
*
      SUBROUTINE CRST(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C                HSF,SPECSA,VOLP,DEF,DNF,DTF,DES,DNS,DTS,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB) 
*
*     Routine to insert additional need terms into LTNE energy equations
*     related to coupling between fluid and solid energy equations.
*
*     AT*  active coefficients for energy equations
*     BT  accumulated fixed source term for energy equations
*     HSF(ID,JD,KD)  interfacial exchange coefficient at porous 
*                    volumes; input
*     SPECSA(ID,JD,KD)  specific surface area array; input
*     VOLP(ID,JD,KD) CV volume; input
*     DEF,DNF,DTF  diffusion coefficients for fluid phase; input
*     DES,DNS,DTS  diffusion coefficients for solid phase; input
*     CVTYPE(ID,JD,KD,NNB+1)  indicates type of each CV; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  array dimensions; input
*     NNB  number of neighbors for interior nodes; input
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD)
      REAL*8 ATS(NT,NT,ID,JD,KD),ATN(NT,NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)      
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 HSF(ID,JD,KD),SPECSA(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 DEF(ID,JD,KD),DNF(ID,JD,KD),DTF(ID,JD,KD)
      REAL*8 DES(ID,JD,KD),DNS(ID,JD,KD),DTS(ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*
      DO 20 I=IB,IE
       DO 10 J=JB,JE
        DO 5 K=KB,KE
         IF(CVTYPE(I,J,K,1).EQ.2) THEN
*
*        Add surface heat exchange terms for porous volumes
*        
          ATP(1,2,I,J,K) = -1.0*HSF(I,J,K)*SPECSA(I,J,K)*VOLP(I,J,K)
          ATP(2,1,I,J,K) = -1.0*HSF(I,J,K)*SPECSA(I,J,K)*VOLP(I,J,K)
*
*        Move diffusion terms at interface adjacent porous 
*        volumes to correct locations in coefficient matrices
*        and add additional cross diffusion terms
*                 
          IF(CVTYPE(I,J,K,2).EQ.0) THEN
           ATW(2,1,I,J,K) = ATW(2,2,I,J,K)
           ATW(2,2,I,J,K) = 0.0
           ATP(1,1,I-1,J,K) = ATP(1,1,I-1,J,K)+DES(I-1,J,K)
           ATE(1,2,I-1,J,K) = ATE(1,2,I-1,J,K)+DES(I-1,J,K)
          ELSEIF(CVTYPE(I,J,K,2).EQ.3) THEN
           ATW(1,2,I,J,K) = ATW(1,1,I,J,K)
           ATW(1,1,I,J,K) = 0.0
           ATP(2,2,I-1,J,K) = ATP(2,2,I-1,J,K)+DEF(I-1,J,K)
           ATE(2,1,I-1,J,K) = ATE(2,1,I-1,J,K)+DEF(I-1,J,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,3).EQ.0) THEN
           ATE(2,1,I,J,K) = ATE(2,2,I,J,K)
           ATE(2,2,I,J,K) = 0.0
           ATP(1,1,I+1,J,K) = ATP(1,1,I+1,J,K)+DES(I,J,K)
           ATW(1,2,I+1,J,K) = ATW(1,2,I+1,J,K)+DES(I,J,K)
          ELSEIF(CVTYPE(I,J,K,3).EQ.3) THEN
           ATE(1,2,I,J,K) = ATE(1,1,I,J,K)
           ATE(1,1,I,J,K) = 0.0
           ATP(2,2,I+1,J,K) = ATP(2,2,I+1,J,K)+DEF(I,J,K)
           ATW(2,1,I+1,J,K) = ATW(2,1,I+1,J,K)+DEF(I,J,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,4).EQ.0) THEN
           ATS(2,1,I,J,K) = ATS(2,2,I,J,K)
           ATS(2,2,I,J,K) = 0.0
           ATP(1,1,I,J-1,K) = ATP(1,1,I,J-1,K)+DNS(I,J-1,K)
           ATN(1,2,I,J-1,K) = ATN(1,2,I,J-1,K)+DNS(I,J-1,K)
          ELSEIF(CVTYPE(I,J,K,4).EQ.3) THEN
           ATS(1,2,I,J,K) = ATS(1,1,I,J,K)
           ATS(1,1,I,J,K) = 0.0
           ATP(2,2,I,J-1,K) = ATP(2,2,I,J-1,K)+DNF(I,J-1,K)
           ATN(2,1,I,J-1,K) = ATN(2,1,I,J-1,K)+DNF(I,J-1,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,5).EQ.0) THEN
           ATN(2,1,I,J,K) = ATN(2,2,I,J,K)
           ATN(2,2,I,J,K) = 0.0
           ATP(1,1,I,J+1,K) = ATP(1,1,I,J+1,K)+DNS(I,J,K)
           ATS(1,2,I,J+1,K) = ATS(1,2,I,J+1,K)+DNS(I,J,K)
          ELSEIF(CVTYPE(I,J,K,5).EQ.3) THEN
           ATN(1,2,I,J,K) = ATN(1,1,I,J,K)
           ATN(1,1,I,J,K) = 0.0
           ATP(2,2,I,J+1,K) = ATP(2,2,I,J+1,K)+DNF(I,J,K)
           ATS(2,1,I,J+1,K) = ATS(2,1,I,J+1,K)+DNF(I,J,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,6).EQ.0) THEN
           ATB(2,1,I,J,K) = ATB(2,2,I,J,K)
           ATB(2,2,I,J,K) = 0.0
           ATP(1,1,I,J,K-1) = ATP(1,1,I,J,K-1)+DTS(I,J,K-1)
           ATT(1,2,I,J,K-1) = ATT(1,2,I,J,K-1)+DTS(I,J,K-1)
          ELSEIF(CVTYPE(I,J,K,6).EQ.3) THEN
           ATB(1,2,I,J,K) = ATB(1,1,I,J,K)
           ATB(1,1,I,J,K) = 0.0
           ATP(2,2,I,J,K-1) = ATP(2,2,I,J,K-1)+DTF(I,J,K-1)
           ATT(2,1,I,J,K-1) = ATT(2,1,I,J,K-1)+DTF(I,J,K-1)
          ENDIF
*
          IF(CVTYPE(I,J,K,7).EQ.0) THEN
           ATT(2,1,I,J,K) = ATT(2,2,I,J,K)
           ATT(2,2,I,J,K) = 0.0
           ATP(1,1,I,J,K+1) = ATP(1,1,I,J,K+1)+DTS(I,J,K)
           ATB(1,2,I,J,K+1) = ATB(1,2,I,J,K+1)+DTS(I,J,K)
          ELSEIF(CVTYPE(I,J,K,7).EQ.3) THEN
           ATT(1,2,I,J,K) = ATT(1,1,I,J,K)
           ATT(1,1,I,J,K) = 0.0
           ATP(2,2,I,J,K+1) = ATP(2,2,I,J,K+1)+DTF(I,J,K)
           ATB(2,1,I,J,K+1) = ATB(2,1,I,J,K+1)+DTF(I,J,K)
          ENDIF          
*
*        Move diffusion terms at solid/fluid interface
*        volumes to correct locations in coefficient matrices
* 
         ELSEIF(CVTYPE(I,J,K,1).EQ.0) THEN   
          IF(CVTYPE(I,J,K,2).EQ.3) THEN
           ATW(1,2,I,J,K) = ATW(1,1,I,J,K)
           ATW(1,1,I,J,K) = 0.0
           ATE(2,1,I-1,J,K) = ATE(2,2,I-1,J,K)
           ATE(2,2,I-1,J,K) = 0.0
          ENDIF
*
          IF(CVTYPE(I,J,K,3).EQ.3) THEN
           ATE(1,2,I,J,K) = ATE(1,1,I,J,K)
           ATE(1,1,I,J,K) = 0.0
           ATW(2,1,I+1,J,K) = ATW(2,2,I+1,J,K)
           ATW(2,2,I+1,J,K) = 0.0
          ENDIF
* 
          IF(CVTYPE(I,J,K,4).EQ.3) THEN
           ATS(1,2,I,J,K) = ATS(1,1,I,J,K)
           ATS(1,1,I,J,K) = 0.0
           ATN(2,1,I,J-1,K) = ATN(2,2,I,J-1,K)
           ATN(2,2,I,J-1,K) = 0.0
          ENDIF
* 
          IF(CVTYPE(I,J,K,5).EQ.3) THEN
           ATN(1,2,I,J,K) = ATN(1,1,I,J,K)
           ATN(1,1,I,J,K) = 0.0
           ATS(2,1,I,J+1,K) = ATS(2,2,I,J+1,K)
           ATS(2,2,I,J+1,K) = 0.0
          ENDIF
* 
          IF(CVTYPE(I,J,K,6).EQ.3) THEN
           ATB(1,2,I,J,K) = ATB(1,1,I,J,K)
           ATB(1,1,I,J,K) = 0.0
           ATT(2,1,I,J,K-1) = ATT(2,2,I,J,K-1)
           ATT(2,2,I,J,K-1) = 0.0
          ENDIF
* 
          IF(CVTYPE(I,J,K,7).EQ.3) THEN
           ATT(1,2,I,J,K) = ATT(1,1,I,J,K)
           ATT(1,1,I,J,K) = 0.0
           ATB(2,1,I,J,K+1) = ATB(2,2,I,J,K+1)
           ATB(2,2,I,J,K+1) = 0.0           
          ENDIF        
         ENDIF
    5   CONTINUE
   10  CONTINUE
   20 CONTINUE
*
      RETURN
      END