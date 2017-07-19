*
***********************************************************************
*
      SUBROUTINE COEFFM(AP,AW,AE,AS,AN,AB,AT,B,
     C                  ME,MN,MT,DE,DN,DT,Q,R,PHI,PHIN,VOLP,
     C                  ALFAE,ALFAN,ALFAT,DIEP,DJNP,DKTP,
     C                  DENST,DTIME,UFACT,CFACT,DFACT,
     C                  IV,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB)
*
*     Subroutine to calculate transport coefficients on a co-located 
*     orthogonal grid.
*
*     AP(N,N,ID,JD,KD) active coefficient for P node; output
*     AW(N,N,ID,JD,KD) active coefficient for W node; output
*     AE(N,N,ID,JD,KD) active coefficient for E node; output
*     AS(N,N,ID,JD,KD) active coefficient for S node; output
*     AN(N,N,ID,JD,KD) active coefficient for N node; output
*     AB(N,N,ID,JD,KD) active coefficient for B node; output
*     AT(N,N,ID,JD,KD) active coefficient for T node; output
*     B(N,ID,JD,KD)    accumulated fixed source term; output
*     ME(ID,JD,KD) mass flux for the east face;  input
*     MN(ID,JD,KD) mass flux for the north face; input
*     MT(ID,JD,KD) mass flux for the top face;  input
*     DE(ID,JD,KD) diffusion coefficient for east face; input
*     DN(ID,JD,KD) diffusion coefficient for north face; input
*     DT(ID,JD,KD) diffusion coefficient for top face; input
*     Q(ID,JD,KD) fixed source term; input
*     R(ID,JD,KD) linearized source term; input
*     PHI(ID,JD,KD) estimate of phi at old time; input
*     PHIN(ID,JD,KD) estimate of phi at new time from last 
*                    iteration; input
*     VOLP(ID,JD,KD) CV volume; input
*     ALFAE(ID,JD) advection weight factor for east face; input
*     ALFAN(ID,JD) advection weight factor for north face: input
*     ALFAT(ID,JD) advection weight factor for top face: input
*     DIEP(ID) distance from P to E through e; input
*     DJNP(JD) distance from P to N through n; input
*     DKTP(KD) distance from P to T through t; input
*     DENST(ID,JD,KD) array of densities for each CV; input
*     DTIME    time step for implicit scheme; input
*     UFACT,CFACT,DFACT factors for unsteady, advective, and 
*                       diffusive terms; input
*     IV       position in block; input
*     IB,IE    first and last interior indices in i; input
*     JB,JE    first and last interior indices in j; input
*     KB,KE    first and last interior indices in k; input
*     N        number of variables in solution; input
*     ID,JD,KD grid-array dimensions; input     
*     NNB  number of neighbors for interior nodes; input
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD),AP(N,N,ID,JD,KD)
      REAL*8 B(N,ID,JD,KD),ASUM,PHI(ID,JD,KD),PHIN(ID,JD,KD)
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      REAL*8 Q(ID,JD,KD),R(ID,JD,KD),VOLP(ID,JD,KD)
      REAL*8 ALFAE(ID,JD,KD),ALFAN(ID,JD,KD),ALFAT(ID,JD,KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DENST(ID,JD,KD),DTIME
      REAL*8 UFACT(ID,JD,KD)
      REAL*8 CFACT(ID,JD,KD,NNB+1),DFACT(ID,JD,KD,NNB+1)
      INTEGER IV,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,NNB
      INTEGER I,J,K
*
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
*
*        Calculate active coefficients of neighbors.
*
         AW(IV,IV,I,J,K)=DFACT(I,J,K,2)*DE(I-1,J,K)
     C                   +0.5*CFACT(I,J,K,2)*(ME(I-1,J,K)
     C                   +ABS(ME(I-1,J,K)*ALFAE(I-1,J,K)))
         AE(IV,IV,I,J,K)=DFACT(I,J,K,3)*DE(I,J,K)
     C                   -0.5*CFACT(I,J,K,3)*(ME(I,J,K)
     C                   -ABS(ME(I,J,K)*ALFAE(I,J,K)))
         AS(IV,IV,I,J,K)=DFACT(I,J,K,4)*DN(I,J-1,K)
     C                   +0.5*CFACT(I,J,K,4)*(MN(I,J-1,K)
     C                   +ABS(MN(I,J-1,K)*ALFAN(I,J-1,K)))
         AN(IV,IV,I,J,K)=DFACT(I,J,K,5)*DN(I,J,K)
     C                   -0.5*CFACT(I,J,K,5)*(MN(I,J,K)
     C                   -ABS(MN(I,J,K)*ALFAN(I,J,K)))
         AB(IV,IV,I,J,K)=DFACT(I,J,K,6)*DT(I,J,K-1)
     C                   +0.5*CFACT(I,J,K,6)*(MT(I,J,K-1)
     C                   +ABS(MT(I,J,K-1)*ALFAT(I,J,K-1)))
         AT(IV,IV,I,J,K)=DFACT(I,J,K,7)*DT(I,J,K)
     C                   -0.5*CFACT(I,J,K,7)*(MT(I,J,K)
     C                   -ABS(MT(I,J,K)*ALFAT(I,J,K)))
*
*        Calculate active coefficient of P volume.
*
         ASUM=DFACT(I,J,K,2)*DE(I-1,J,K)+DFACT(I,J,K,3)*DE(I,J,K)
     C      +DFACT(I,J,K,4)*DN(I,J-1,K)+DFACT(I,J,K,5)*DN(I,J,K)
     C      +DFACT(I,J,K,6)*DT(I,J,K-1)+DFACT(I,J,K,7)*DT(I,J,K)
     C      +0.5*CFACT(I,J,K,1)*(ME(I-1,J,K)
     C      +ABS(ME(I-1,J,K)*ALFAE(I-1,J,K))
     C      -ME(I,J,K)+ABS(ME(I,J,K)*ALFAE(I,J,K))
     C      +MN(I,J-1,K)+ABS(MN(I,J-1,K)*ALFAN(I,J-1,K)) 
     C      -MN(I,J,K)+ABS(MN(I,J,K)*ALFAN(I,J,K))
     C      +MT(I,J,K-1)+ABS(MT(I,J,K-1)*ALFAT(I,J,K-1)) 
     C      -MT(I,J,K)+ABS(MT(I,J,K)*ALFAT(I,J,K)))
         AP(IV,IV,I,J,K)=UFACT(I,J,K)*DENST(I,J,K)*VOLP(I,J,K)/DTIME
     C                 +ASUM-R(I,J,K)
*
*        Calculate fixed term.
*
         B(IV,I,J,K)=B(IV,I,J,K)+Q(I,J,K)
     C     +UFACT(I,J,K)*PHI(I,J,K)*DENST(I,J,K)*VOLP(I,J,K)/DTIME
*

 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
*
      RETURN
      END
