*             file dhat.f
*     
************************************************************************
*
      SUBROUTINE DHAT(DHUE,DHVN,DHWT,VOLE,VOLN,VOLT,
     C                AEAST,ANORTH,ATOP, AUP,VOLP,PRSTY,RHO,
     C                DIEP,DISE,DJNP,DISN,DKTP,DIST,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*     Subroutine to calculate the diffusion coefficients for  
*     the special momentum equation required to compute the 
*     advecting velocity.See thesis for details of implementation.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 DHUE(ID,JD,KD),DHVN(ID,JD,KD),DHWT(ID,JD,KD)
      REAL*8 VOLE(ID,JD,KD),VOLN(ID,JD,KD),VOLT(ID,JD,KD)
      REAL*8 AEAST(ID,JD,KD),ANORTH(ID,JD,KD),ATOP(ID,JD,KD)
      REAL*8 AUP(N,N,ID,JD,KD),VOLP(ID,JD,KD),PRSTY(ID,JD,KD),RHO
      REAL*8 DIEP(ID),DISE(ID),DJNP(JD),DISN(JD),DKTP(KD),DIST(KD)
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*
*  Compute diffusion coefficients for i direction
*
      DO 25 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE-1
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,3).NE.3) THEN
          AEAST(I,J,K)= 0.5*(AUP(2,2,I,J,K)+AUP(2,2,I+1,J,K))
          VOLE(I,J,K)= 0.5*(VOLP(I,J,K)+VOLP(I+1,J,K))
          DHUE(I,J,K)= VOLE(I,J,K)/(AEAST(I,J,K)+1.E-20)
         ENDIF
  10    CONTINUE
  20   CONTINUE
  25  CONTINUE
*
*  Compute diffusion coefficients for j direction
*
      DO 40 I=IB,IE
       DO 35 J=JB,JE-1
        DO 30 K=KB,KE
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,5).NE.3) THEN
          ANORTH(I,J,K)= 0.5*(AUP(3,3,I,J,K)+AUP(3,3,I,J+1,K))
          VOLN(I,J,K)= 0.5*(VOLP(I,J,K)+VOLP(I,J+1,K))
          DHVN(I,J,K)= VOLN(I,J,K)/(ANORTH(I,J,K)+1.E-20)
         ENDIF
  30    CONTINUE
  35   CONTINUE
  40  CONTINUE
*
*  Compute diffusion coefficients for k direction
*
      DO 55 I=IB,IE
       DO 50 J=JB,JE
        DO 45 K=KB,KE-1
         IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,7).NE.3) THEN
          ATOP(I,J,K)= 0.5*(AUP(4,4,I,J,K)+AUP(4,4,I,J,K+1))
          VOLT(I,J,K)= 0.5*(VOLP(I,J,K)+VOLP(I,J,K+1))
          DHWT(I,J,K)= VOLT(I,J,K)/(ATOP(I,J,K)+1.E-20)
         ENDIF
  45    CONTINUE
  50   CONTINUE
  55  CONTINUE  
*
      RETURN
      END

