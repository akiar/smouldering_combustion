*
*                    file vhat.f
**********************************************************************
*
      SUBROUTINE VHAT(VHN, V,P,PNORTH,DHVN,VOLN,VOLP,
     C                DPDYP,DPDYM,ANORTH,AUP,PRSTY,
     C                RHO,DJNP,DISN,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*     
*     Subroutine to calculate Vhn at the C.V. faces.
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 VHN(ID,JD,KD),V(ID,JD,KD)
      REAL*8 P(ID,JD,KD),PNORTH(ID,JD,KD),DHVN(ID,JD,KD)
      REAL*8 VOLN(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 DPDYP(ID,JD,KD),DPDYM(ID,JD,KD)
      REAL*8 ANORTH(ID,JD,KD),AUP(N,N,ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 DJNP(JD),DISN(JD)
      REAL*8 DPDYN,FACEPS
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*
      REAL*8 RHO(ID,JD,KD)
*
      DO 3 I=IB,IE
      DO 2 K=KB,KE
       VHN(I,JB-1,K)= V(I,JB-1,K)
       DO 1 J=JB,JE-1
        IF(CVTYPE(I,J,K,1).EQ.3 .OR. CVTYPE(I,J,K,5).EQ.3) THEN
         VHN(I,J,K)= 0.0
        ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,5)) THEN
         DPDYN= 0.5*(PRSTY(I,J,K)*VOLP(I,J,K)*(PNORTH(I,J,K)
     C          -P(I,J,K))/DISN(J)+PRSTY(I,J+1,K)*VOLP(I,J+1,K)
     C          *(P(I,J+1,K)-PNORTH(I,J,K))/DISN(J+1))/VOLN(I,J,K)
         VHN(I,J,K)= ((AUP(3,3,I,J,K)/ANORTH(I,J,K))*V(I,J,K)
     C    +(AUP(3,3,I,J+1,K)/ANORTH(I,J,K))*V(I,J+1,K))/2.0
     C    -DHVN(I,J,K)*(DPDYN-(PRSTY(I,J,K)*VOLP(I,J,K)*DPDYP(I,J,K)
     C    +PRSTY(I,J+1,K)*VOLP(I,J+1,K)*DPDYM(I,J+1,K))
     C    /(2.0*VOLN(I,J,K)))
        ELSE
         FACEPS= (PRSTY(I,J,K)*DISN(J+1)
     C           +PRSTY(I,J+1,K)*DISN(J))/DJNP(J)        
         VHN(I,J,K)= (V(I,J,K)*DISN(J+1)+V(I,J+1,K)*DISN(J))/DJNP(J)
     C    -DHVN(I,J,K)*(FACEPS*(P(I,J+1,K)-P(I,J,K))/DJNP(J)
     C    -(DPDYP(I,J,K)*PRSTY(I,J,K)*DISN(J+1)
     C    +DPDYM(I,J+1,K)*PRSTY(I,J+1,K)*DISN(J))/DJNP(J))
        ENDIF
  1    CONTINUE
       VHN(I,JE,K)= V(I,JE+1,K)
  2   CONTINUE
  3   CONTINUE
*
      RETURN
      END

