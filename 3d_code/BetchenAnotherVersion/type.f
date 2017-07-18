*     file type.f
***********************************************************************
*
      SUBROUTINE TYPE(CVTYPE, POROUS,IPB,IPE,JPB,JPE,KPB,KPE,
     C                IPB2,IPE2,JPB2,JPE2,KPB2,KPE2,
     C                SOLID,ISB,ISE,JSB,JSE,KSB,KSE,
     C                ISB2,ISE2,JSB2,JSE2,KSB2,KSE2,
     C                IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*  Routine to create array CVTYPE, which indicates type of each CV
*  as follows:
*
*        CVTYPE(I,J,K,L) = 0: CV L contains clear fluid
*        CVTYPE(I,J,K,L) = 1: CV L on boundary
*        CVTYPE(I,J,K,L) = 2: CV L contains porous medium
*        CVTYPE(I,J,K,L) = 3: CV L contains solid
*  
*  Where:
*
*        L=1: P node
*        L=2: W node
*        L=3: E node
*        L=4: S node
*        L=5: N node
*        L=6: B node
*        L=7: T node
*
***********************************************************************
*
      IMPLICIT NONE
      INTEGER POROUS,IPB,IPE,JPB,JPE,KPB,KPE
      INTEGER IPB2,IPE2,JPB2,JPE2,KPB2,KPE2
      INTEGER SOLID,ISB,ISE,JSB,JSE,KSB,KSE
      INTEGER ISB2,ISE2,JSB2,JSE2,KSB2,KSE2
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K,L
*
      L = 1
      DO 3 I=IB,IE
       DO 2 J=JB,JE
        DO 1 K=KB,KE
         CVTYPE(I,J,K,L) = 0
    1   CONTINUE
    2  CONTINUE
    3 CONTINUE
*
*     Porous media control volumes
*
      IF(POROUS.NE.0) THEN
        DO 15 I=IPB,IPE
         DO 10 J=JPB,JPE
          DO 5 K=KPB,KPE
           CVTYPE(I,J,K,L) = 2
    5     CONTINUE
   10    CONTINUE
   15   CONTINUE
        DO 30 I=IPB2,IPE2
         DO 25 J=JPB2,JPE2
          DO 20 K=KPB2,KPE2
           CVTYPE(I,J,K,L) = 2
   20     CONTINUE
   25    CONTINUE
   30   CONTINUE
      ENDIF
*
*     Solid control volumes
*
      IF(SOLID.NE.0) THEN
        DO 45 I=ISB,ISE
         DO 40 J=JSB,JSE
          DO 35 K=KSB,KSE
           CVTYPE(I,J,K,L) = 3
   35     CONTINUE
   40    CONTINUE
   45   CONTINUE
        DO 60 I=ISB2,ISE2
         DO 55 J=JSB2,JSE2
          DO 50 K=KSB2,KSE2
           CVTYPE(I,J,K,L) = 3
   50     CONTINUE
   55    CONTINUE
   60   CONTINUE
      ENDIF
*
*     Bottom and top faces of domain
*
      DO 80 I=IB,IE
       DO 70 J=JB,JE
        CVTYPE(I,J,KB-1,1) = 1
        CVTYPE(I,J,KB-1,7) = CVTYPE(I,J,KB,1)
        CVTYPE(I,J,KE+1,1) = 1
        CVTYPE(I,J,KE+1,6) = CVTYPE(I,J,KE,1)   
   70  CONTINUE     
   80 CONTINUE
*
*     South and north faces of domain
*
      DO 100 I=IB,IE
       DO 90 K=KB,KE
        CVTYPE(I,JB-1,K,1) = 1
        CVTYPE(I,JB-1,K,5) = CVTYPE(I,JB,K,1)
        CVTYPE(I,JE+1,K,1) = 1
        CVTYPE(I,JE+1,K,4) = CVTYPE(I,JE,K,1)        
   90  CONTINUE
  100 CONTINUE
*
*     East and west faces of domain
*
      DO 120 J=JB,JE
       DO 110 K=KB,KE
        CVTYPE(IB-1,J,K,1) = 1
        CVTYPE(IB-1,J,K,3) = CVTYPE(IB,J,K,1)        
        CVTYPE(IE+1,J,K,1) = 1
        CVTYPE(IE+1,J,K,2) = CVTYPE(IE,J,K,1)        
  110  CONTINUE
  120 CONTINUE
*
*     Interior neighbor nodes
*
      DO 140 I=IB,IE
       DO 130 J=JB,JE
        DO 125 K=KB,KE
         L=2
         CVTYPE(I,J,K,L) = CVTYPE(I-1,J,K,1)
*
         L=3
         CVTYPE(I,J,K,L) = CVTYPE(I+1,J,K,1)
*
         L=4
         CVTYPE(I,J,K,L) = CVTYPE(I,J-1,K,1)
*
         L=5
         CVTYPE(I,J,K,L) = CVTYPE(I,J+1,K,1)
*
         L=6
         CVTYPE(I,J,K,L) = CVTYPE(I,J,K-1,1)
*
         L=7
         CVTYPE(I,J,K,L) = CVTYPE(I,J,K+1,1)     
  125   CONTINUE     
  130  CONTINUE
  140 CONTINUE  
*
      RETURN
      END