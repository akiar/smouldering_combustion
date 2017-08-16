*
*                    file uhat.f
**********************************************************************
*
      SUBROUTINE UHAT(UHE, U,P,PEAST,DHUE,VOLE,VOLP,
     C                DPDXP,DPDXM,AEAST,AUP,PRSTY,
     C                RHO,DIEP,DISE,
     C                CVTYPE,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*     
*     Subroutine to calculate Uhe at the C.V. faces.
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 UHE(ID,JD,KD),U(ID,JD,KD)
      REAL*8 P(ID,JD,KD),PEAST(ID,JD,KD),DHUE(ID,JD,KD)
      REAL*8 VOLE(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 DPDXP(ID,JD,KD),DPDXM(ID,JD,KD)
      REAL*8 AEAST(ID,JD,KD),AUP(N,N,ID,JD,KD),PRSTY(ID,JD,KD)
      REAL*8 DIEP(ID),DISE(ID)
      REAL*8 DPDXE,FACEPS
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1),I,J,K
*
      REAL*8 RHO(ID,JD,KD)
*
      DO 3 J=JB,JE
      DO 2 K=KB,KE
       UHE(IB-1,J,K)= U(IB-1,J,K)
       DO 1 I=IB,IE-1
        IF(CVTYPE(I,J,K,1).EQ.3 .OR. CVTYPE(I,J,K,3).EQ.3) THEN
         UHE(I,J,K)= 0.0
        ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,3)) THEN
         DPDXE= 0.5*(PRSTY(I,J,K)*VOLP(I,J,K)*(PEAST(I,J,K)
     C          -P(I,J,K))/DISE(I)+PRSTY(I+1,J,K)*VOLP(I+1,J,K)
     C          *(P(I+1,J,K)-PEAST(I,J,K))/DISE(I+1))/VOLE(I,J,K)
         UHE(I,J,K)= ((AUP(2,2,I,J,K)/AEAST(I,J,K))*U(I,J,K)
     C    +(AUP(2,2,I+1,J,K)/AEAST(I,J,K))*U(I+1,J,K))/2.0
     C    -DHUE(I,J,K)*(DPDXE-(PRSTY(I,J,K)*VOLP(I,J,K)*DPDXP(I,J,K)
     C    +PRSTY(I+1,J,K)*VOLP(I+1,J,K)*DPDXM(I+1,J,K))
     C    /(2.0*VOLE(I,J,K)))     
        ELSE
         FACEPS= (PRSTY(I,J,K)*DISE(I+1)
     C           +PRSTY(I+1,J,K)*DISE(I))/DIEP(I)
	 UHE(I,J,K)= (U(I,J,K)*DISE(I+1)+U(I+1,J,K)*DISE(I))/DIEP(I)
     C    -DHUE(I,J,K)*(FACEPS*(P(I+1,J,K)-P(I,J,K))/DIEP(I)
     C    -(DPDXP(I,J,K)*PRSTY(I,J,K)*DISE(I+1)
     C    +DPDXM(I+1,J,K)*PRSTY(I+1,J,K)*DISE(I))/DIEP(I))
        ENDIF
  1    CONTINUE
       UHE(IE,J,K)= U(IE+1,J,K)
  2   CONTINUE
  3   CONTINUE
*
      RETURN
      END

