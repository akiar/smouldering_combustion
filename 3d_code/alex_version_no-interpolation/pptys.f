*     file pptys.f
***********************************************************************
*
      SUBROUTINE PPTYS(PPTY, PPTF,PPTP,PPTS,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Routine to create property arrays.
*
*     PPTY(ID,JD,KD)  property array; output
*     PPTF  property in fluid region; input
*     PPTP  property in porous region; input
*     PPTS  property in solid region; input
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
      REAL*8 PPTY(ID,JD,KD),PPTF,PPTP,PPTS
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
*     Interior volumes
*
      DO 20 I=IB,IE      
       DO 10 J=JB,JE
        DO 5 K=KB,KE
          IF(CVTYPE(I,J,K,1).EQ.0) THEN
            PPTY(I,J,K) = PPTF
          ELSEIF(CVTYPE(I,J,K,1).EQ.2) THEN
            PPTY(I,J,K) = PPTP
          ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            PPTY(I,J,K) = PPTS            
          ENDIF
    5   CONTINUE
   10  CONTINUE
   20 CONTINUE
*
*     Boundary nodes
*
      DO 30 I=IB,IE
       DO 25 K=KB,KE
        J = JB-1
        PPTY(I,J,K) = PPTY(I,J+1,K)
*
        J = JE+1        
        PPTY(I,J,K) = PPTY(I,J-1,K)
   25  CONTINUE
   30 CONTINUE
      DO 40 J=JB,JE
       DO 35 K=KB,KE
        I = IB-1
        PPTY(I,J,K) = PPTY(I+1,J,K)
*
        I = IE+1        
        PPTY(I,J,K) = PPTY(I-1,J,K)
   35  CONTINUE
   40 CONTINUE   
      DO 50 I=IB,IE
       DO 45 J=JB,JE
        K = KB-1
        PPTY(I,J,K) = PPTY(I,J,K+1)
*
        K = KE+1        
        PPTY(I,J,K) = PPTY(I,J,K-1)
   45  CONTINUE
   50 CONTINUE    
*
      RETURN
      END