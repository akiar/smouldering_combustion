*
*     This file contains 2 subroutines;  RESIDT and RESIDM
*
***********************************************************************
*
      SUBROUTINE RESIDT(RSD,AVRSD,RSDMAX,I_MAX,J_MAX,K_MAX,
     C                  TF,TS,AP,AW,AE,AS,AN,AB,AT,B,CARR,
     C                  IV,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*     Subroutine to calculate the residual at each interior c.v and
*     the average of the absolute residuals over all interior c.v.
*     Used for energy equations.
*
*     RSD(ID,JD,KD)  residual array for each interior c.v; output
*     AVRSD  average residual for all interior c.v.; output
*     RSDMAX  maximum residual; output
*     I_MAX,J_MAX,K_MAX  location of maximum residual; output
*     TF,TS            fluid and solid temperatures; input
*     AP(N,N,ID,JD,KD) active coefficient for P node; input
*     AW(N,N,ID,JD,KD) active coefficient for W node; input
*     AE(N,N,ID,JD,KD) active coefficient for E node; input
*     AS(N,N,ID,JD,KD) active coefficient for S node; input
*     AN(N,N,ID,JD,KD) active coefficient for N node; input
*     AB(N,N,ID,JD,KD) active coefficient for B node; input
*     AT(N,N,ID,JD,KD) active coefficient for T node; input
*     B(N,ID,JD,KD)    accumulated fixed source term; input
*     CARR(ID,JD,KD)   array of specific heats for normalization;
*                      input
*     IV        equation number to check; input
*     IB,IE     first and last interior indices in i; input
*     JB,JE     first and last interior indices in j; input
*     KB,KE     first and last interior indices in k; input
*     NT        number of simultaneously solved variables; input
*     ID,JD,KD  array dimensions; input
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 RSD(ID,JD,KD),AVRSD,RSDMAX,TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 AW(NT,NT,ID,JD,KD),AE(NT,NT,ID,JD,KD)
      REAL*8 AS(NT,NT,ID,JD,KD),AN(NT,NT,ID,JD,KD)
      REAL*8 AB(NT,NT,ID,JD,KD),AT(NT,NT,ID,JD,KD)
      REAL*8 AP(NT,NT,ID,JD,KD),B(NT,ID,JD,KD),CARR(ID,JD,KD)
      INTEGER IV,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD
      INTEGER I,J,K,NPTS,I_MAX,J_MAX,K_MAX
*
      AVRSD=0.0
      RSDMAX= 0.0
      NPTS=(IE-IB+1)*(JE-JB+1)*(KE-KB+1)
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
          RSD(I,J,K)=
     C      (AW(IV,1,I,J,K)*TF(I-1,J,K)+AE(IV,1,I,J,K)*TF(I+1,J,K)
     C      +AS(IV,1,I,J,K)*TF(I,J-1,K)+AN(IV,1,I,J,K)*TF(I,J+1,K)
     C      +AB(IV,1,I,J,K)*TF(I,J,K-1)+AT(IV,1,I,J,K)*TF(I,J,K+1)
     C      -AP(IV,1,I,J,K)*TF(I,J,K)
     C      +AW(IV,2,I,J,K)*TS(I-1,J,K)+AE(IV,2,I,J,K)*TS(I+1,J,K)
     C      +AS(IV,2,I,J,K)*TS(I,J-1,K)+AN(IV,2,I,J,K)*TS(I,J+1,K)
     C      +AB(IV,2,I,J,K)*TS(I,J,K-1)+AT(IV,2,I,J,K)*TS(I,J,K+1)
     C      -AP(IV,2,I,J,K)*TS(I,J,K)
     C      +B(IV,I,J,K))/CARR(I,J,K)
          IF(ABS(RSD(I,J,K)).GT.RSDMAX) THEN
            RSDMAX= ABS(RSD(I,J,K))
            I_MAX= I
            J_MAX= J
            K_MAX= K
          ENDIF
          AVRSD=AVRSD+ABS(RSD(I,J,K))
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      AVRSD=AVRSD/NPTS
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE RESIDM(RSD,AVRSD,RSDMAX,I_MAX,J_MAX,K_MAX,
     C                  P,U,V,W,AP,AW,AE,AS,AN,AB,AT,B,
     C                  IV,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*     Subroutine to calculate the residual at each interior c.v and
*     the average of the absolute residuals over all interior c.v.
*     Used for mass/momentum set.
*
*     RSD(ID,JD,KD)  residual array for each interior c.v; output
*     AVRSD  average residual for all interior c.v.; output
*     RSDMAX  maximum residual; output
*     I_MAX,J_MAX,K_MAX  location of maximum residual; output
*     P,U,V,W          fluid and solid temperatures; input
*     AP(N,N,ID,JD,KD) active coefficient for P node; input
*     AW(N,N,ID,JD,KD) active coefficient for W node; input
*     AE(N,N,ID,JD,KD) active coefficient for E node; input
*     AS(N,N,ID,JD,KD) active coefficient for S node; input
*     AN(N,N,ID,JD,KD) active coefficient for N node; input
*     AB(N,N,ID,JD,KD) active coefficient for B node; input
*     AT(N,N,ID,JD,KD) active coefficient for T node; input
*     B(N,ID,JD,KD)    accumulated fixed source term; input
*     IV        equation number to check; input
*     IB,IE     first and last interior indices in i; input
*     JB,JE     first and last interior indices in j; input
*     KB,KE     first and last interior indices in k; input
*     NT        number of simultaneously solved variables; input
*     ID,JD,KD  array dimensions; input
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 RSD(ID,JD,KD),AVRSD,RSDMAX
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD) 
      REAL*8 AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD)      
      REAL*8 AP(N,N,ID,JD,KD),B(N,ID,JD,KD)
      INTEGER IV,IB,IE,JB,JE,KB,KE,N,ID,JD,KD
      INTEGER I,J,K,NPTS,I_MAX,J_MAX,K_MAX
*
      AVRSD=0.0
      RSDMAX= 0.0
      NPTS=(IE-IB+1)*(JE-JB+1)*(KE-KB+1)
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
          RSD(I,J,K)=
     C      AW(IV,1,I,J,K)*P(I-1,J,K)+AE(IV,1,I,J,K)*P(I+1,J,K)
     C      +AS(IV,1,I,J,K)*P(I,J-1,K)+AN(IV,1,I,J,K)*P(I,J+1,K)
     C      +AB(IV,1,I,J,K)*P(I,J,K-1)+AT(IV,1,I,J,K)*P(I,J,K+1)     
     C      -AP(IV,1,I,J,K)*P(I,J,K)
     C      +AW(IV,2,I,J,K)*U(I-1,J,K)+AE(IV,2,I,J,K)*U(I+1,J,K)
     C      +AS(IV,2,I,J,K)*U(I,J-1,K)+AN(IV,2,I,J,K)*U(I,J+1,K)
     C      +AB(IV,2,I,J,K)*U(I,J,K-1)+AT(IV,2,I,J,K)*U(I,J,K+1)     
     C      -AP(IV,2,I,J,K)*U(I,J,K)
     C      +AW(IV,3,I,J,K)*V(I-1,J,K)+AE(IV,3,I,J,K)*V(I+1,J,K)
     C      +AS(IV,3,I,J,K)*V(I,J-1,K)+AN(IV,3,I,J,K)*V(I,J+1,K)
     C      +AB(IV,3,I,J,K)*V(I,J,K-1)+AT(IV,3,I,J,K)*V(I,J,K+1)     
     C      -AP(IV,3,I,J,K)*V(I,J,K)
     C      +AW(IV,4,I,J,K)*W(I-1,J,K)+AE(IV,4,I,J,K)*W(I+1,J,K)
     C      +AS(IV,4,I,J,K)*W(I,J-1,K)+AN(IV,4,I,J,K)*W(I,J+1,K)
     C      +AB(IV,4,I,J,K)*W(I,J,K-1)+AT(IV,4,I,J,K)*W(I,J,K+1)     
     C      -AP(IV,4,I,J,K)*W(I,J,K)
     C      +B(IV,I,J,K)
          IF(ABS(RSD(I,J,K)).GT.RSDMAX) THEN
            RSDMAX= ABS(RSD(I,J,K))
            I_MAX= I
            J_MAX= J
            K_MAX= K
          ENDIF
          AVRSD=AVRSD+ABS(RSD(I,J,K))
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
      AVRSD=AVRSD/NPTS
*
      RETURN
      END

