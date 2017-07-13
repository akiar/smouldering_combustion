*
*     This file contains 1 main subroutine: BCELIM.  The
*     additional routines are used in BCELIM for matrix algebra.
*
**********************************************************************
      SUBROUTINE BCELIM(AP,AW,AE,AS,AN,AB,AT,BP,A,B,
     C                  IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*     Subroutine to absorb the boundary conditions into the equations
*     at the first interior nodes.  Used for N simultaneous variables.
*
*     AP(N,N,ID,JD,KD) active coefficient for P node; output
*     AW(N,N,ID,JD,KD) active coefficient for W node; output
*     AE(N,N,ID,JD,KD) active coefficient for E node; output
*     AS(N,N,ID,JD,KD) active coefficient for S node; output
*     AN(N,N,ID,JD,KD) active coefficient for N node; output
*     AB(N,N,ID,JD,KD) active coefficient for B node; output
*     AT(N,N,ID,JD,KD) active coefficient for T node; output
*     B(N,ID,JD,KD)    accumulated fixed source term; output
*
*     INTEGER IB,IE  first and last interior indices in i; input
*     INTEGER JB,JE  first and last interior indices in j; input
*     INTEGER KB,KE  first and last interior indices in k; input
*     INTEGER ID,JD,KD  array dimensions; input
*     INTEGER N      number of simultaneous variables; input
*
*     NOTES: 1.) corner control-volumes are visited twice since
*                they are adjacent to two edges.
*
*            2.) for solutions of simultaneous variables, b.c.'s are
*                absorbed using matrix algebra.
*
*            3.) the indices in the active coefficients are:
*
*                  A(row, column, i-position, j-position, k-position)
*
*                  B(row, i-position, j-position, k-position)
*
*     IMPORTANT: Routine currently assumes that AP(N,N,IB-1,J,K),
*                AP(N,N,IE+1,J,K), AP(N,N,I,JB-1,K), AP(N,N,I,JE+1,K),
*                AP(N,N,I,J,KB-1), and AP(N,N,I,J,KE+1) are NxN 
*                identity matrices.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 AP(N,N,ID,JD,KD),AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD),BP(N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD)
      REAL*8 A(N,N),B(N)
      INTEGER IB,IE,JB,JE,KB,KE
      INTEGER I,J,K,N,ID,JD,KD
*     
*  West and East faces of domain
*
      DO 15 K=KB,KE
       DO 10 J=JB,JE
         I=IB
         CALL MMULT(A,  AW,AE,I,J,K,I-1,J,K,N,ID,JD,KD)
         CALL MVMULT(B, AW,BP,I,J,K,I-1,J,K,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AW, I,J,K,N,ID,JD,KD)
*
         I=IE
         CALL MMULT(A,  AE,AW,I,J,K,I+1,J,K,N,ID,JD,KD)
         CALL MVMULT(B, AE,BP,I,J,K,I+1,J,K,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AE, I,J,K,N,ID,JD,KD)
 10    CONTINUE
 15   CONTINUE
*
*  South and North faces of domain
*
      DO 25 K=KB,KE
       DO 20 I=IB,IE
         J=JB
         CALL MMULT(A,  AS,AN,I,J,K,I,J-1,K,N,ID,JD,KD)
         CALL MVMULT(B, AS,BP,I,J,K,I,J-1,K,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AS, I,J,K,N,ID,JD,KD)
*
         J=JE
         CALL MMULT(A,  AN,AS,I,J,K,I,J+1,K,N,ID,JD,KD)
         CALL MVMULT(B, AN,BP,I,J,K,I,J+1,K,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AN, I,J,K,N,ID,JD,KD)
 20    CONTINUE
 25   CONTINUE
*
*  Top and Bottom faces of domain
*
      DO 35 J=JB,JE
       DO 30 I=IB,IE
         K=KB
         CALL MMULT(A,  AB,AT,I,J,K,I,J,K-1,N,ID,JD,KD)
         CALL MVMULT(B, AB,BP,I,J,K,I,J,K-1,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AB, I,J,K,N,ID,JD,KD)
*
         K=KE
         CALL MMULT(A,  AT,AB,I,J,K,I,J,K+1,N,ID,JD,KD)
         CALL MVMULT(B, AT,BP,I,J,K,I,J,K+1,N,ID,JD,KD)
         CALL MSUBT(AP, AP,A,I,J,K,N,ID,JD,KD)
         CALL VPLUS(BP, BP,B,I,J,K,N,ID,JD,KD)
         CALL ZEROM(AT, I,J,K,N,ID,JD,KD)
 30    CONTINUE
 35   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE VPLUS(B3, B1,B2,I,J,K,N,ID,JD,KD)
*
*     Routine to add two vectors, B1 and B2. The result is stored in
*     vector B3 on output.  B2 has dimension B2(N).
*
************************************************************************
      IMPLICIT NONE
      REAL*8 B3(N,ID,JD,KD),B1(N,ID,JD,KD),B2(N)
      INTEGER H,I,J,K,N,ID,JD,KD
*
      DO 10 H=1,N
        B3(H,I,J,K)=B1(H,I,J,K)+B2(H)
 10   CONTINUE
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MSUBT(A3, A1,A2,I,J,K,N,ID,JD,KD)
*
*     Routine to subtract two matrices, A1 and A2. The result (A1-A2) 
*     is stored in matrix A3 on output.  A2 has dimension A3(N,N).
*
************************************************************************
      IMPLICIT NONE
      REAL*8 A3(N,N,ID,JD,KD),A1(N,N,ID,JD,KD),A2(N,N)
      INTEGER H,I,J,K,L,N,ID,JD,KD
*
      DO 20 H=1,N
        DO 10 L=1,N
          A3(H,L,I,J,K)= A1(H,L,I,J,K)-A2(H,L)
 10     CONTINUE
 20   CONTINUE
* 
      RETURN
      END
*
************************************************************************
      SUBROUTINE MMULT(A3, A1,A2,I1,J1,K1,I2,J2,K2,N,ID,JD,KD)
*
*     Routine to multiply two matrices, A1 and A2. The result (A1*A2)
*     is stored in A3 on output.  A3 has dimension A3(N,N).
*
*************************************************************************
      IMPLICIT NONE
      REAL*8 A3(N,N),A1(N,N,ID,JD,KD),A2(N,N,ID,JD,KD)
      INTEGER I1,J1,K1,I2,J2,K2,H,L,M,N,ID,JD,KD
*
      DO 30 H=1,N
        DO 20 L=1,N
          A3(H,L)= 0.0
          DO 10 M=1,N
            A3(H,L)= A1(H,M,I1,J1,K1)*A2(M,L,I2,J2,K2)+A3(H,L)
 10       CONTINUE
 20     CONTINUE
 30   CONTINUE
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE MVMULT(B2, A1,B1,I1,J1,K1,I2,J2,K2,N,ID,JD,KD)
*
*     Routine to multiply a vector, B1, by a matrix, A1. The result
*     is stored in the vector B2 on output.  B2 has dimension B2(N).
*
***********************************************************************
      IMPLICIT NONE
      REAL*8 B2(N),A1(N,N,ID,JD,KD),B1(N,ID,JD,KD)
      INTEGER I1,J1,K1,I2,J2,K2,L,M,N,ID,JD,KD
*
      DO 20 L=1,N
        B2(L)= 0.0
        DO 10 M=1,N
          B2(L)= A1(L,M,I1,J1,K1)*B1(M,I2,J2,K2)+B2(L)
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END
*
***********************************************************************
      SUBROUTINE ZEROM(A1, I,J,K,N,ID,JD,KD)
*
*     Routine to zero a matrix A1.
*
***********************************************************************
      IMPLICIT NONE
      REAL*8 A1(N,N,ID,JD,KD)
      INTEGER H,I,J,K,L,N,ID,JD,KD
*
      DO 20 H=1,N
        DO 10 L=1,N
          A1(H,L,I,J,K)= 0.0
 10     CONTINUE
 20   CONTINUE
*
      RETURN
      END

