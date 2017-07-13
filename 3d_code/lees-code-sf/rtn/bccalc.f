*
*     This file contains 2 subroutines: BCCLCT and BCCLCU
*
**********************************************************************
      SUBROUTINE BCCLCT(TF,TS, AP,AW,AE,AS,AN,AB,AT,B,
     C                  IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*     Subroutine to update the phi field boundary values.  Used for
*     LTNE energy equations.
*
*     TF,TS updated estimate of temperature fields; output
*
*     AP(N,N,ID,JD,KD) active coefficient for P node; input
*     AW(N,N,ID,JD,KD) active coefficient for W node; input
*     AE(N,N,ID,JD,KD) active coefficient for E node; input
*     AS(N,N,ID,JD,KD) active coefficient for S node; input
*     AN(N,N,ID,JD,KD) active coefficient for N node; input
*     AB(N,N,ID,JD,KD) active coefficient for B node; input
*     AT(N,N,ID,JD,KD) active coefficient for T node; input
*     B(N,ID,JD,KD) accumulated fixed source term; input
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER KB,KE  first and last interior indices in k; input
*     INTEGER ID,JD,KD  array dimensions; input
*     INTEGER N      number of simultaneous variables; input
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD)
      REAL*8 AP(N,N,ID,JD,KD),B(N,ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,I,J,K,L,N,ID,JD,KD
*     
*  West and East faces of domain
*
      DO 10 J=JB,JE
       DO 5 K=KB,KE
*
        I=IB-1     
        TF(I,J,K)= (AE(1,1,I,J,K)*TF(I+1,J,K)
     C             +AE(1,2,I,J,K)*TS(I+1,J,K)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AE(2,1,I,J,K)*TF(I+1,J,K)
     C             +AE(2,2,I,J,K)*TS(I+1,J,K)
     C             +B(2,I,J,K))/AP(2,2,I,J,K)     
*
        I=IE+1
        TF(I,J,K)= (AW(1,1,I,J,K)*TF(I-1,J,K)
     C             +AW(1,2,I,J,K)*TS(I-1,J,K)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AW(2,1,I,J,K)*TF(I-1,J,K)
     C             +AW(2,2,I,J,K)*TS(I-1,J,K)
     C             +B(2,I,J,K))/AP(2,2,I,J,K) 
*
  5    CONTINUE       
 10   CONTINUE
*
*  South and North faces of domain
*     
      DO 20 I=IB,IE
       DO 15 K=KB,KE
*
        J=JB-1
        TF(I,J,K)= (AN(1,1,I,J,K)*TF(I,J+1,K)
     C             +AN(1,2,I,J,K)*TS(I,J+1,K)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AN(2,1,I,J,K)*TF(I,J+1,K)
     C             +AN(2,2,I,J,K)*TS(I,J+1,K)
     C             +B(2,I,J,K))/AP(2,2,I,J,K)     
*
        J=JE+1
        TF(I,J,K)= (AS(1,1,I,J,K)*TF(I,J-1,K)
     C             +AS(1,2,I,J,K)*TS(I,J-1,K)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AS(2,1,I,J,K)*TF(I,J-1,K)
     C             +AS(2,2,I,J,K)*TS(I,J-1,K)
     C             +B(2,I,J,K))/AP(2,2,I,J,K)     
*   
 15    CONTINUE     
 20   CONTINUE
*
*  Bottom and Top faces of domain
*     
      DO 30 I=IB,IE
       DO 25 J=JB,JE
*
        K=KB-1
        TF(I,J,K)= (AT(1,1,I,J,K)*TF(I,J,K+1)
     C             +AT(1,2,I,J,K)*TS(I,J,K+1)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AT(2,1,I,J,K)*TF(I,J,K+1)
     C             +AT(2,2,I,J,K)*TS(I,J,K+1)
     C             +B(2,I,J,K))/AP(2,2,I,J,K)     
*
        K=KE+1
        TF(I,J,K)= (AB(1,1,I,J,K)*TF(I,J,K-1)
     C             +AB(1,2,I,J,K)*TS(I,J,K-1)
     C             +B(1,I,J,K))/AP(1,1,I,J,K)
        TS(I,J,K)= (AB(2,1,I,J,K)*TF(I,J,K-1)
     C             +AB(2,2,I,J,K)*TS(I,J,K-1)
     C             +B(2,I,J,K))/AP(2,2,I,J,K)     
*   
 25    CONTINUE     
 30   CONTINUE 
*
*  Edges and corners of domain
*
      DO 40 K=KB,KE
       TF(IB-1,JB-1,K)= TF(IB-1,JB,K)
       TF(IE+1,JB-1,K)= TF(IE,JB-1,K)
       TF(IB-1,JE+1,K)= TF(IB-1,JE,K)
       TF(IE+1,JE+1,K)= TF(IE,JE+1,K)      
       TS(IB-1,JB-1,K)= TS(IB-1,JB,K)
       TS(IE+1,JB-1,K)= TS(IE,JB-1,K)
       TS(IB-1,JE+1,K)= TS(IB-1,JE,K)
       TS(IE+1,JE+1,K)= TS(IE,JE+1,K)
 40   CONTINUE
      DO 50 I=IB,IE
       TF(I,JB-1,KB-1)= TF(I,JB,KB-1)
       TF(I,JE+1,KB-1)= TF(I,JE,KB-1)
       TF(I,JB-1,KE+1)= TF(I,JB,KE+1)
       TF(I,JE+1,KE+1)= TF(I,JE,KE+1)      
       TS(I,JB-1,KB-1)= TS(I,JB,KB-1)
       TS(I,JE+1,KB-1)= TS(I,JE,KB-1)
       TS(I,JB-1,KE+1)= TS(I,JB,KE+1)
       TS(I,JE+1,KE+1)= TS(I,JE,KE+1)
 50   CONTINUE
      DO 60 J=JB,JE
       TF(IB-1,J,KB-1)= TF(IB-1,J,KB)
       TF(IE+1,J,KB-1)= TF(IE,J,KB-1)
       TF(IB-1,J,KE+1)= TF(IB-1,J,KE)
       TF(IE+1,J,KE+1)= TF(IE,J,KE+1)
       TS(IB-1,J,KB-1)= TS(IB-1,J,KB)
       TS(IE+1,J,KB-1)= TS(IE,J,KB-1)
       TS(IB-1,J,KE+1)= TS(IB-1,J,KE)
       TS(IE+1,J,KE+1)= TS(IE,J,KE+1)
 60   CONTINUE
*
      TF(IB-1,JB-1,KB-1)= TF(IB-1,JB,KB-1)
      TF(IE+1,JB-1,KB-1)= TF(IE,JB-1,KB-1)
      TF(IB-1,JE+1,KB-1)= TF(IB-1,JE,KB-1)
      TF(IE+1,JE+1,KB-1)= TF(IE+1,JB,KB-1)
      TF(IB-1,JB-1,KE+1)= TF(IB-1,JB,KE+1)
      TF(IE+1,JB-1,KE+1)= TF(IE,JB-1,KE+1)
      TF(IB-1,JE+1,KE+1)= TF(IB-1,JE,KE+1)
      TF(IE+1,JE+1,KE+1)= TF(IE+1,JB,KE+1)
      TS(IB-1,JB-1,KB-1)= TS(IB-1,JB,KB-1)
      TS(IE+1,JB-1,KB-1)= TS(IE,JB-1,KB-1)
      TS(IB-1,JE+1,KB-1)= TS(IB-1,JE,KB-1)
      TS(IE+1,JE+1,KB-1)= TS(IE+1,JB,KB-1)
      TS(IB-1,JB-1,KE+1)= TS(IB-1,JB,KE+1)
      TS(IE+1,JB-1,KE+1)= TS(IE,JB-1,KE+1)
      TS(IB-1,JE+1,KE+1)= TS(IB-1,JE,KE+1)
      TS(IE+1,JE+1,KE+1)= TS(IE+1,JB,KE+1)
*
      RETURN
      END
*      
**********************************************************************
      SUBROUTINE BCCLCU(PHI, AP,AW,AE,AS,AN,AB,AT,B,
     C                  IV,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*     Subroutine to update the phi field boundary values.  Used for
*     N simulataneously solved variables.
*
*     PHI(ID,JD,KD) updated estimate of phi field; output
*
*     AP(N,N,ID,JD,KD) active coefficient for P node; input
*     AW(N,N,ID,JD,KD) active coefficient for W node; input
*     AE(N,N,ID,JD,KD) active coefficient for E node; input
*     AS(N,N,ID,JD,KD) active coefficient for S node; input
*     AN(N,N,ID,JD,KD) active coefficient for N node; input
*     AB(N,N,ID,JD,KD) active coefficient for B node; input
*     AT(N,N,ID,JD,KD) active coefficient for T node; input
*     B(N,ID,JD,KD) accumulated fixed source term; input
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER KB,KE  first and last interior indices in k; input
*     INTEGER ID,JD,KD  array dimensions; input
*     INTEGER N      number of simultaneous variables; input
*
**********************************************************************
*
      IMPLICIT NONE
      REAL*8 AP(N,N,ID,JD,KD),AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD),B(N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD),PHI(ID,JD,KD)
      INTEGER IV,IB,IE,JB,JE,KB,KE,I,J,K,N,ID,JD,KD
*     
*  West and East faces of domain
*
      DO 15 K=KB,KE
       DO 10 J=JB,JE
*
        I=IB-1
        PHI(I,J,K)= (AE(IV,IV,I,J,K)*PHI(I+1,J,K)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
*
        I=IE+1
        PHI(I,J,K)= (AW(IV,IV,I,J,K)*PHI(I-1,J,K)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
 10    CONTINUE
 15   CONTINUE
*
*  South and North faces of domain
*     
      DO 25 K=KB,KE
       DO 20 I=IB,IE
*
        J=JB-1
        PHI(I,J,K)= (AN(IV,IV,I,J,K)*PHI(I,J+1,K)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
*
        J=JE+1
        PHI(I,J,K)= (AS(IV,IV,I,J,K)*PHI(I,J-1,K)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
 20    CONTINUE
 25   CONTINUE
*
*  Bottom and Top faces of domain
*     
      DO 35 J=JB,JE
       DO 30 I=IB,IE
*
        K=KB-1
        PHI(I,J,K)= (AT(IV,IV,I,J,K)*PHI(I,J,K+1)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
*
        K=KE+1
        PHI(I,J,K)= (AB(IV,IV,I,J,K)*PHI(I,J,K-1)
     C              +B(IV,I,J,K))/AP(IV,IV,I,J,K)
 30    CONTINUE
 35   CONTINUE
*
*  Edges and corners of domain
*
      DO 40 K=KB,KE
       PHI(IB-1,JB-1,K)= PHI(IB-1,JB,K)
       PHI(IE+1,JB-1,K)= PHI(IE,JB-1,K)
       PHI(IB-1,JE+1,K)= PHI(IB-1,JE,K)
       PHI(IE+1,JE+1,K)= PHI(IE,JE+1,K)
 40   CONTINUE
      DO 50 I=IB,IE
       PHI(I,JB-1,KB-1)= PHI(I,JB,KB-1)
       PHI(I,JE+1,KB-1)= PHI(I,JE,KB-1)
       PHI(I,JB-1,KE+1)= PHI(I,JB,KE+1)
       PHI(I,JE+1,KE+1)= PHI(I,JE,KE+1)
 50   CONTINUE
      DO 60 J=JB,JE
       PHI(IB-1,J,KB-1)= PHI(IB-1,J,KB)
       PHI(IE+1,J,KB-1)= PHI(IE,J,KB-1)
       PHI(IB-1,J,KE+1)= PHI(IB-1,J,KE)
       PHI(IE+1,J,KE+1)= PHI(IE,J,KE+1)
 60   CONTINUE
*
      PHI(IB-1,JB-1,KB-1)= PHI(IB-1,JB,KB-1)
      PHI(IE+1,JB-1,KB-1)= PHI(IE,JB-1,KB-1)
      PHI(IB-1,JE+1,KB-1)= PHI(IB-1,JE,KB-1)
      PHI(IE+1,JE+1,KB-1)= PHI(IE+1,JB,KB-1)
      PHI(IB-1,JB-1,KE+1)= PHI(IB-1,JB,KE+1)
      PHI(IE+1,JB-1,KE+1)= PHI(IE,JB-1,KE+1)
      PHI(IB-1,JE+1,KE+1)= PHI(IB-1,JE,KE+1)
      PHI(IE+1,JE+1,KE+1)= PHI(IE+1,JB,KE+1)
*
      RETURN
      END
