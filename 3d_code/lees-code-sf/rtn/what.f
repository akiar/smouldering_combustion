*
*                    file what.f
**********************************************************************
*
      SUBROUTINE WHAT(WHT, W,P,PTOP,DHWT,VOLT,VOLP,
     C                DPDZP,DPDZM,ATOP,AUP,PRSTY,
     C                RHO,DKTP,DIST,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*     
*     Subroutine to calculate Wht at the C.V. faces.
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 WHT(ID,JD,KD),W(ID,JD,KD)
      REAL*8 P(ID,JD,KD),PTOP(ID,JD,KD),DHWT(ID,JD,KD)
      REAL*8 VOLT(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 DPDZP(ID,JD,KD),DPDZM(ID,JD,KD)
      REAL*8 ATOP(ID,JD,KD),AUP(N,N,ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 RHO,DKTP(KD),DIST(KD)
      REAL*8 DPDZT,FACEPS
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*
      DO 3 I=IB,IE
      DO 2 J=JB,JE
       WHT(I,J,KB-1)= W(I,J,KB-1)
       DO 1 K=KB,KE-1
        IF(CVTYPE(I,J,K,1).EQ.3 .OR. CVTYPE(I,J,K,7).EQ.3) THEN
         WHT(I,J,K)= 0.0
        ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,7)) THEN
         DPDZT= 0.5*(PRSTY(I,J,K)*VOLP(I,J,K)*(PTOP(I,J,K)
     C          -P(I,J,K))/DIST(K)+PRSTY(I,J,K+1)*VOLP(I,J,K+1)
     C          *(P(I,J,K+1)-PTOP(I,J,K))/DIST(K+1))/VOLT(I,J,K)
         WHT(I,J,K)= ((AUP(4,4,I,J,K)/ATOP(I,J,K))*W(I,J,K)
     C    +(AUP(4,4,I,J,K+1)/ATOP(I,J,K))*W(I,J,K+1))/2.0
     C    -DHWT(I,J,K)*(DPDZT-(PRSTY(I,J,K)*VOLP(I,J,K)*DPDZP(I,J,K)
     C    +PRSTY(I,J,K+1)*VOLP(I,J,K+1)*DPDZM(I,J,K+1))
     C    /(2.0*VOLT(I,J,K)))
        ELSE
         FACEPS= (PRSTY(I,J,K)*DIST(K+1)
     C           +PRSTY(I,J,K+1)*DIST(K))/DKTP(K)        
         WHT(I,J,K)= (W(I,J,K)*DIST(K+1)+W(I,J,K+1)*DIST(K))/DKTP(K)
     C    -DHWT(I,J,K)*(FACEPS*(P(I,J,K+1)-P(I,J,K))/DKTP(K)
     C    -(DPDZP(I,J,K)*PRSTY(I,J,K)*DIST(K+1)
     C    +DPDZM(I,J,K+1)*PRSTY(I,J,K+1)*DIST(K))/DKTP(K))
        ENDIF
  1    CONTINUE
       WHT(I,J,KE)= W(I,J,KE+1)
  2   CONTINUE
  3   CONTINUE
*
      RETURN
      END

