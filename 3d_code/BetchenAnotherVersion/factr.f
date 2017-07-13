*     file factr.f
***********************************************************************
*
      SUBROUTINE FACTR(UNSTF,CONTF,DIFTF,
     C                 UNSTS,CONTS,DIFTS,
     C                 UNSU,CONU,DIFU,     
     C                 PRSTY,CP,CSP,CS,MUEFF,VISCB,
     C                 DISE,DIEP,DISN,DJNP,DIST,DKTP,     
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB) 
*
*     Routine to create factor arrays.
*
*     UNSTF,UNSTS,UNSU  factor arrays for unsteady terms for fluid
*                       energy, solid energy, and momentum equations;
*                       output
*     CONTF,CONTS,CONU  factor arrays for advective terms for fluid
*                       energy, solid energy, and momentum equations;
*                       output
*     DIFTF,DIFTS,DIFU  factor arrays for diffusive terms for fluid
*                       energy, solid energy, and momentum equations;
*                       output
*
*     Notes:  For CON*(I,J,K,L):  L=1: P node
*                                 L=2: W node
*                                 L=3: E node
*                                 L=4: S node
*                                 L=5: N node
*                                 L=6: B node
*                                 L=7: T node
*             For DIF*(I,J,K,L):  L=2: w face
*                                 L=3: e face
*                                 L=4: s face
*                                 L=5: n face
*                                 L=6: b face
*                                 L=7: t face
*
*             Form of unsteady term is UNS*DENST*VOLP(PHI-PHIOLD)/DTIME
*             Form of advective term for e is ME*(CON(L=3)*FE*PHI_E
*                                                  +CON(L=3)*FP*PHI_P),
*             where FE and FP are the weightings of the advection scheme                                                  
*             Form of diffusive term for e is DIF(L=3)*DE*(PHI_E-PHI_P)
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 UNSTF(ID,JD,KD)
      REAL*8 CONTF(ID,JD,KD,NNB+1),DIFTF(ID,JD,KD,NNB+1)
      REAL*8 UNSTS(ID,JD,KD)
      REAL*8 CONTS(ID,JD,KD,NNB+1),DIFTS(ID,JD,KD,NNB+1)
      REAL*8 UNSU(ID,JD,KD)
      REAL*8 CONU(ID,JD,KD,NNB+1),DIFU(ID,JD,KD,NNB+1) 
      REAL*8 PRSTY(ID,JD,KD),CP,CSP,CS,MUEFF(ID,JD,KD),VISCB
      REAL*8 DISE(ID),DIEP(ID),DISN(JD),DJNP(JD),DIST(KD),DKTP(KD)
      REAL*8 GAMMA
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L
*
      DO 20 I=IB,IE      
       DO 15 J=JB,JE
        DO 10 K=KB,KE
*
         L=1
         IF(CVTYPE(I,J,K,L).EQ.0) THEN
          UNSTF(I,J,K) = CP
          UNSTS(I,J,K) = 0.0
          UNSU(I,J,K) = 1.0                    
          CONTF(I,J,K,L) = CP
          CONTS(I,J,K,L) = 0.0
          CONU(I,J,K,L) = 1.0                    
         ELSEIF(CVTYPE(I,J,K,L).EQ.2) THEN
          UNSTF(I,J,K) = PRSTY(I,J,K)*CP
          UNSTS(I,J,K) = (1.0-PRSTY(I,J,K))*CSP
          UNSU(I,J,K) = 1.0                    
          CONTF(I,J,K,L) = CP
          CONTS(I,J,K,L) = 0.0
          CONU(I,J,K,L) = 1.0/PRSTY(I,J,K)                    
         ELSEIF(CVTYPE(I,J,K,L).EQ.3) THEN
          UNSTF(I,J,K) = 0.0
          UNSTS(I,J,K) = CS
          UNSU(I,J,K) = 0.0                    
          CONTF(I,J,K,L) = 0.0
          CONTS(I,J,K,L) = 0.0
          CONU(I,J,K,L) = 0.0
         ENDIF
*
*
*
         DO 5 L=2,7
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
           DIFTF(I,J,K,L) = 1.0
           DIFTS(I,J,K,L) = 0.0           
           DIFU(I,J,K,L) = 1.0           
          ELSEIF(CVTYPE(I,J,K,1).EQ.2) THEN
           DIFTF(I,J,K,L) = 1.0
           DIFTS(I,J,K,L) = 1.0           
*
           IF(CVTYPE(I,J,K,L).EQ.2) THEN
            IF(L.EQ.2) THEN
             GAMMA = 1.0/((DISE(I)/DIEP(I-1))/MUEFF(I,J,K)
     C                   +(DISE(I-1)/DIEP(I-1))/MUEFF(I-1,J,K))
            ELSEIF(L.EQ.3) THEN
             GAMMA = 1.0/((DISE(I)/DIEP(I))/MUEFF(I,J,K)
     C                   +(DISE(I+1)/DIEP(I))/MUEFF(I+1,J,K))
            ELSEIF(L.EQ.4) THEN
             GAMMA = 1.0/((DISN(J)/DJNP(J-1))/MUEFF(I,J,K)
     C                   +(DISN(J-1)/DJNP(J-1))/MUEFF(I,J-1,K))
            ELSEIF(L.EQ.5) THEN
             GAMMA = 1.0/((DISN(J)/DJNP(J))/MUEFF(I,J,K)
     C                   +(DISN(J+1)/DJNP(J))/MUEFF(I,J+1,K))
            ELSEIF(L.EQ.6) THEN
             GAMMA = 1.0/((DIST(K)/DKTP(K-1))/MUEFF(I,J,K)
     C                   +(DIST(K-1)/DKTP(K-1))/MUEFF(I,J,K-1))
            ELSEIF(L.EQ.7) THEN
             GAMMA = 1.0/((DIST(K)/DKTP(K))/MUEFF(I,J,K)
     C                   +(DIST(K+1)/DKTP(K))/MUEFF(I,J,K+1))
            ENDIF           
            DIFU(I,J,K,L) = VISCB/GAMMA
           ELSE
            DIFU(I,J,K,L) = PRSTY(I,J,K)          
           ENDIF
*
          ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN          
           DIFTF(I,J,K,L) = 0.0
           DIFTS(I,J,K,L) = 1.0           
           DIFU(I,J,K,L) = 0.0          
          ENDIF
*          
*
          IF(CVTYPE(I,J,K,L).EQ.0) THEN
           CONTF(I,J,K,L) = CP
           CONTS(I,J,K,L) = 0.0
           CONU(I,J,K,L) = 1.0     
          ELSEIF(CVTYPE(I,J,K,L).EQ.1) THEN
           CONTF(I,J,K,L) = CONTF(I,J,K,1)
           CONTS(I,J,K,L) = CONTS(I,J,K,1)
           CONU(I,J,K,L) = CONU(I,J,K,1)
          ELSEIF(CVTYPE(I,J,K,L).EQ.2) THEN
           CONTF(I,J,K,L) = CP
           CONTS(I,J,K,L) = 0.0
*           
           IF(L.EQ.2) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I-1,J,K) 
           ELSEIF(L.EQ.3) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I+1,J,K) 
           ELSEIF(L.EQ.4) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I,J-1,K) 
           ELSEIF(L.EQ.5) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I,J+1,K) 
           ELSEIF(L.EQ.6) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I,J,K-1) 
           ELSEIF(L.EQ.7) THEN
            CONU(I,J,K,L) = 1.0/PRSTY(I,J,K+1) 
           ENDIF
*           
          ELSEIF(CVTYPE(I,J,K,L).EQ.3) THEN 
           CONTF(I,J,K,L) = 0.0
           CONTS(I,J,K,L) = 0.0
           CONU(I,J,K,L) = 0.0
          ENDIF          
    5    CONTINUE
   10   CONTINUE    
   15  CONTINUE
   20 CONTINUE
*
*     At boundaries
*
      DO 25 I=IB,IE
       DO 22 K=KB,KE
        J=JB-1
        CONTF(I,J,K,1) = CONTF(I,J+1,K,4)
        CONTF(I,J,K,5) = CONTF(I,J+1,K,1)
        CONTS(I,J,K,1) = CONTS(I,J+1,K,4)
        CONTS(I,J,K,5) = CONTS(I,J+1,K,1)
        CONU(I,J,K,1) = CONU(I,J+1,K,4)
        CONU(I,J,K,5) = CONU(I,J+1,K,1)                
*        
        J=JE+1
        CONTF(I,J,K,1) = CONTF(I,J-1,K,5)
        CONTF(I,J,K,4) = CONTF(I,J-1,K,1)
        CONTS(I,J,K,1) = CONTS(I,J-1,K,5)
        CONTS(I,J,K,4) = CONTS(I,J-1,K,1)
        CONU(I,J,K,1) = CONU(I,J-1,K,5)
        CONU(I,J,K,4) = CONU(I,J-1,K,1)                
   22  CONTINUE
   25 CONTINUE
*
      DO 30 J=JB,JE
       DO 27 K=KB,KE
        I=IB-1
        CONTF(I,J,K,1) = CONTF(I+1,J,K,2)
        CONTF(I,J,K,3) = CONTF(I+1,J,K,1)
        CONTS(I,J,K,1) = CONTS(I+1,J,K,2)
        CONTS(I,J,K,3) = CONTS(I+1,J,K,1)        
        CONU(I,J,K,1) = CONU(I+1,J,K,2)
        CONU(I,J,K,3) = CONU(I+1,J,K,1)        
*
        I=IE+1
        CONTF(I,J,K,1) = CONTF(I-1,J,K,3)
        CONTF(I,J,K,2) = CONTF(I-1,J,K,1)
        CONTS(I,J,K,1) = CONTS(I-1,J,K,3)
        CONTS(I,J,K,2) = CONTS(I-1,J,K,1)
        CONU(I,J,K,1) = CONU(I-1,J,K,3)
        CONU(I,J,K,2) = CONU(I-1,J,K,1)                
   27  CONTINUE
   30 CONTINUE
*   
      DO 40 I=IB,IE
       DO 35 J=JB,JE
        K=KB-1
        CONTF(I,J,K,1) = CONTF(I,J,K+1,6)
        CONTF(I,J,K,7) = CONTF(I,J,K+1,1)
        CONTS(I,J,K,1) = CONTS(I,J,K+1,6)
        CONTS(I,J,K,7) = CONTS(I,J,K+1,1)
        CONU(I,J,K,1) = CONU(I,J,K+1,6)
        CONU(I,J,K,7) = CONU(I,J,K+1,1)
*
        K=KE+1
        CONTF(I,J,K,1) = CONTF(I,J,K-1,7)
        CONTF(I,J,K,6) = CONTF(I,J,K-1,1)
        CONTS(I,J,K,1) = CONTS(I,J,K-1,7)
        CONTS(I,J,K,6) = CONTS(I,J,K-1,1)
        CONU(I,J,K,1) = CONU(I,J,K-1,7)
        CONU(I,J,K,6) = CONU(I,J,K-1,1)        
   35  CONTINUE
   40 CONTINUE
*   
      RETURN
      END