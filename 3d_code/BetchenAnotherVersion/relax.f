*     file relax.f
***********************************************************************
*
      SUBROUTINE RELAXT(ATP,BT, TF,TS,RLX,
     C                  IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*     Routine to add relaxation parameters to the energy equations.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD),RLX,RLXADD
      INTEGER IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,I,J,K
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
*
*       Fluid energy equation
*        
         RLXADD = ATP(1,1,I,J,K)*RLX
         ATP(1,1,I,J,K) = ATP(1,1,I,J,K)+RLXADD
         BT(1,I,J,K) = BT(1,I,J,K)+RLXADD*TF(I,J,K)
*
*       Solid energy equation
*        
         RLXADD = ATP(2,2,I,J,K)*RLX
         ATP(2,2,I,J,K) = ATP(2,2,I,J,K)+RLXADD
         BT(2,I,J,K) = BT(2,I,J,K)+RLXADD*TS(I,J,K)         
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE UNRLXT(ATP,BT, TF,TS,RLX,
     C                  IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*     Routine to remove relaxation parameters from the energy 
*     equations.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD),RLX,RLXSUB
      INTEGER IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,I,J,K
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
*
*       Fluid energy equation
*        
         RLXSUB = ATP(1,1,I,J,K)*RLX/(1.0+RLX)
         ATP(1,1,I,J,K) = ATP(1,1,I,J,K)-RLXSUB
         BT(1,I,J,K) = BT(1,I,J,K)-RLXSUB*TF(I,J,K)
*
*       Solid energy equation
*        
         RLXSUB = ATP(2,2,I,J,K)*RLX/(1.0+RLX)
         ATP(2,2,I,J,K) = ATP(2,2,I,J,K)-RLXSUB
         BT(2,I,J,K) = BT(2,I,J,K)-RLXSUB*TS(I,J,K)         
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE
*
      RETURN
      END      
* 
***********************************************************************
*
      SUBROUTINE  RELAXU(AUP,BU, U,V,W,RLX,
     C                   IB,IE,JB,JE,KB,KE,N,ID,JD,KD)

*
*     Routine to add relaxation parameters to the momentum equations.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 RLX,RLXADD
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,I,J,K
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
*
*       X-momentum equation
*        
         RLXADD = AUP(2,2,I,J,K)*RLX
         AUP(2,2,I,J,K) = AUP(2,2,I,J,K)+RLXADD
         BU(2,I,J,K) = BU(2,I,J,K)+RLXADD*U(I,J,K)       
*
*       Y-momentum equation
*        
         RLXADD = AUP(3,3,I,J,K)*RLX
         AUP(3,3,I,J,K) = AUP(3,3,I,J,K)+RLXADD
         BU(3,I,J,K) = BU(3,I,J,K)+RLXADD*V(I,J,K)   
*
*       Z-momentum equation
*        
         RLXADD = AUP(4,4,I,J,K)*RLX
         AUP(4,4,I,J,K) = AUP(4,4,I,J,K)+RLXADD
         BU(4,I,J,K) = BU(4,I,J,K)+RLXADD*W(I,J,K)                     
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE  UNRLXU(AUP,BU, U,V,W,RLX,
     C                   IB,IE,JB,JE,KB,KE,N,ID,JD,KD)

*
*     Routine to remove relaxation parameters from the momentum 
*     equations.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 AUP(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 RLX,RLXSUB
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,I,J,K
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
*
*       X-momentum equation
*        
         RLXSUB = AUP(2,2,I,J,K)*RLX/(1.0+RLX)
         AUP(2,2,I,J,K) = AUP(2,2,I,J,K)-RLXSUB
         BU(2,I,J,K) = BU(2,I,J,K)-RLXSUB*U(I,J,K)       
*
*       Y-momentum equation
*        
         RLXSUB = AUP(3,3,I,J,K)*RLX/(1.0+RLX)
         AUP(3,3,I,J,K) = AUP(3,3,I,J,K)-RLXSUB
         BU(3,I,J,K) = BU(3,I,J,K)-RLXSUB*V(I,J,K)   
*
*       Z-momentum equation
*        
         RLXSUB = AUP(4,4,I,J,K)*RLX/(1.0+RLX)
         AUP(4,4,I,J,K) = AUP(4,4,I,J,K)-RLXSUB
         BU(4,I,J,K) = BU(4,I,J,K)-RLXSUB*W(I,J,K)                     
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE
*
      RETURN
      END