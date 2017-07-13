*
*     file globm.f
*********************************************************************
*
      SUBROUTINE GLOBMS(GLOBM,RHS,IA,JA,
     C                  AP,AW,AE,AS,AN,AB,AT,B,
     C                  IB,IE,JB,JE,KB,KE,N,ID,JD,KD,
     C                  MG,MR,NELT,NEQ)
*
*
*     Subroutine to construct the coefficient array in the structure 
*     required by DLAP.
*
*     GLOBM(NELT)    global array in IA(),JA() format; output
*     RHS(NEQ)       right-hand-side vector; output
*     NELT           number of non-zero elements in GLOBM(); output
*     NEQ            order of global coefficient matrix; input
*     
*     AU#(N,N,ID,JD) NxN active coefficient block arrays
*     IB,IE          starting and ending indices in 'x' direction
*     JB,JE          starting and ending indices in 'y' direction
*     ID,JD          array dimensions of block coefficients
*
*     MG,MR          maximum array lengths of GLOBM() and RHS()     
*
*     NOTE:  The structure of the matrix is set up such that the
*            diagonal of each E+H-1 is the first element in the packed
*            global matrix.
*
*       Here is an example of the SLAP Triad storage format for a
*       5x5 Matrix.  The entries may appear in any order.
*
*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
*                              1  2  3  4  5  6  7  8  9 10 11
*       |11 12  0  0 15| GLOBM: 11 12 15 22 21 33 35 44 55 51 53
*       |21 22  0  0  0|    IA:  1  1  1  2  2  3  3  4  5  5  5
*       | 0  0 33  0 35|    JA:  1  2  3  2  1  3  5  4  5  1  3
*       | 0  0  0 44  0|
*       |51  0 53  0 55|
*
*                            Here; NEQ=5, NELT=11      
*           
*
*********************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION GLOBM(MG),RHS(MR)
      REAL*8 AP(N,N,ID,JD,KD),AW(N,N,ID,JD,KD),AE(N,N,ID,JD,KD)
      REAL*8 AS(N,N,ID,JD,KD),AN(N,N,ID,JD,KD),B(N,ID,JD,KD)
      REAL*8 AB(N,N,ID,JD,KD),AT(N,N,ID,JD,KD)
      INTEGER IA(MG),JA(MG)
      INTEGER IB,IE,JB,JE,KB,KE,E,F,H,L,N,ID,JD,KD
      INTEGER MG,MR,NELT,NEQ
      INTEGER W,Z,I,J,K
*
      W=N*(IE-IB+1)
      Z=(JE-JB+1)
*
*-------------------------------------------------------
*  [1] Form the IA(), JA() arrays
*-------------------------------------------------------
*
*  Initialize indices and counters
*
      I=IB
      J=JB
      K=KB
      E=1
      F=1
*
*  Bottom face control-volumes--------------------------------
*
      DO 21 H=1,N
        DO 22 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          F=F+1
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IF(JE .GT. JB) THEN
            IA(F)=E+H-1
            JA(F)=E+W+L-1
            GLOBM(F)= DBLE(-AN(H,L,I,J,K))
            F=F+1
          ENDIF
          IF(KE .GT. KB) THEN
            IA(F)=E+H-1
            JA(F)=E+W*Z+L-1
            GLOBM(F)= DBLE(-AT(H,L,I,J,K))
            F=F+1
          ENDIF
 22     CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 21   CONTINUE  
*
      DO 50 I=IB+1,IE-1
        E=E+N
        DO 51 H=1,N
          DO 52 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-N+L-1
            GLOBM(F+1)= DBLE(-AW(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E+N+L-1
            GLOBM(F+2)= DBLE(-AE(H,L,I,J,K))
            F=F+3
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 52       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 51     CONTINUE       
 50   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 56 H=1,N
          DO 57 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-N+L-1
            GLOBM(F+1)= DBLE(-AW(H,L,I,J,K))
            F=F+2
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 57       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))    
 56     CONTINUE       
      ENDIF
*------------------------------------------------------
      DO 70 J=JB+1,JE-1
*
        I=IB
        E=E+N
        DO 71 H=1,N
          DO 72 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W+L-1
            GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
            F=F+2
            IF(IE .GT. IB) THEN
              IA(F)=E+H-1
              JA(F)=E+N+L-1
              GLOBM(F)= DBLE(-AE(H,L,I,J,K))
              F=F+1
            ENDIF
            IA(F)=E+H-1
            JA(F)=E+W+L-1
            GLOBM(F)= DBLE(-AN(H,L,I,J,K))
            F=F+1
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 72       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 71     CONTINUE           
*
      DO 73 I=IB+1,IE-1
        E=E+N
        DO 74 H=1,N
          DO 75 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W+L-1
            GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E+N+L-1
            GLOBM(F+3)= DBLE(-AE(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+W+L-1
            GLOBM(F+4)= DBLE(-AN(H,L,I,J,K))
            F=F+5
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 75       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 74     CONTINUE 
 73   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 76 H=1,N
          DO 77 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W+L-1
            GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E+W+L-1
            GLOBM(F+3)= DBLE(-AN(H,L,I,J,K))
            F=F+4
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 77       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 76     CONTINUE 
      ENDIF
*
 70   CONTINUE
*-------------------------------------------------      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      E=E+N
      DO 78 H=1,N
        DO 79 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W+L-1
          GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
          F=F+2
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IF(KE .GT. KB) THEN
            IA(F)=E+H-1
            JA(F)=E+W*Z+L-1
            GLOBM(F)= DBLE(-AT(H,L,I,J,K))
            F=F+1
          ENDIF
 79     CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 78   CONTINUE 
*
      DO 80 I=IB+1,IE-1
        E=E+N
        DO 81 H=1,N
          DO 82 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W+L-1
            GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E+N+L-1
            GLOBM(F+3)= DBLE(-AE(H,L,I,J,K))
            F=F+4
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 82       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 81     CONTINUE 
 80   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 86 H=1,N
          DO 87 L=1,N
            IA(F)=E+H-1
            JA(F) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W+L-1
            GLOBM(F+1)= DBLE(-AS(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            F=F+3
            IF(KE .GT. KB) THEN
              IA(F)=E+H-1
              JA(F)=E+W*Z+L-1
              GLOBM(F)= DBLE(-AT(H,L,I,J,K))
              F=F+1
            ENDIF
 87       CONTINUE
            RHS(E+H-1)= DBLE( B(H,I,J,K))
 86     CONTINUE
      ENDIF
*
        ENDIF
*        PRINT *,'FIRST GLOBM DONE'
*
*  Interior control-volumes--------------------------
*
      DO 150 K=KB+1,KE-1
*
      I=IB
      J=JB
      E=E+N
      DO 88 H=1,N
        DO 89 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          F=F+2
          IF(IE .GT. IB) THEN
          IA(F)=E+H-1
          JA(F)=E+N+L-1
          GLOBM(F)= DBLE(-AE(H,L,I,J,K))
          F=F+1
          ENDIF
          IF(JE .GT. JB) THEN
          IA(F)=E+H-1
          JA(F)=E+W+L-1
          GLOBM(F)= DBLE(-AN(H,L,I,J,K))
          F=F+1
          ENDIF
          IA(F)=E+H-1
          JA(F)=E+W*Z+L-1
          GLOBM(F)= DBLE(-AT(H,L,I,J,K))
          F=F+1
 89     CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 88   CONTINUE

*
      DO 90 I=IB+1,IE-1
        E=E+N
        DO 91 H=1,N
          DO 92 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E+N+L-1
            GLOBM(F+3)= DBLE(-AE(H,L,I,J,K))
            F=F+4
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
            IA(F)=E+H-1
            JA(F)=E+W*Z+L-1
            GLOBM(F)= DBLE(-AT(H,L,I,J,K))
            F=F+1
 92       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 91     CONTINUE           
 90   CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 96 H=1,N
          DO 97 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            F=F+3
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
            IA(F)=E+H-1
            JA(F)=E+W*Z+L-1
            GLOBM(F)= DBLE(-AT(H,L,I,J,K))
            F=F+1
 97       CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 96     CONTINUE           
      ENDIF
*------------------------------------------------------
      DO 110 J=JB+1,JE-1
*
      I=IB
      E=E+N
      DO 98 H=1,N
        DO 99 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          IA(F+2)=E+H-1
          JA(F+2)=E-W+L-1
          GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
          F=F+3
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IA(F)=E+H-1
          JA(F)=E+W+L-1
          GLOBM(F)= DBLE(-AN(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E+W*Z+L-1
          GLOBM(F+1)= DBLE(-AT(H,L,I,J,K))
          F=F+2
 99     CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 98   CONTINUE         
*
      DO 100 I=IB+1,IE-1
        E=E+N
        DO 101 H=1,N
          DO 102 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+N+L-1
            GLOBM(F+4)= DBLE(-AE(H,L,I,J,K))
            IA(F+5)=E+H-1
            JA(F+5)=E+W+L-1
            GLOBM(F+5)= DBLE(-AN(H,L,I,J,K))
            IA(F+6)=E+H-1
            JA(F+6)=E+W*Z+L-1
            GLOBM(F+6)= DBLE(-AT(H,L,I,J,K))
            F=F+7
 102      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 101    CONTINUE           
 100  CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 106 H=1,N
          DO 107 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+W+L-1
            GLOBM(F+4)= DBLE(-AN(H,L,I,J,K))
            IA(F+5)=E+H-1
            JA(F+5)=E+W*Z+L-1
            GLOBM(F+5)= DBLE(-AT(H,L,I,J,K))
            F=F+6
 107      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 106    CONTINUE           
      ENDIF
      
*
 110  CONTINUE
*-----------------------------------------------      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      E=E+N
      DO 108 H=1,N
        DO 109 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          IA(F+2)=E+H-1
          JA(F+2)=E-W+L-1
          GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
          F=F+3
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IA(F)=E+H-1
          JA(F)=E+W*Z+L-1
          GLOBM(F)= DBLE(-AT(H,L,I,J,K))
          F=F+1
 109    CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 108  CONTINUE         
*
      DO 120 I=IB+1,IE-1
        E=E+N
        DO 121 H=1,N
          DO 122 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+N+L-1
            GLOBM(F+4)= DBLE(-AE(H,L,I,J,K))
            IA(F+5)=E+H-1
            JA(F+5)=E+W*Z+L-1
            GLOBM(F+5)= DBLE(-AT(H,L,I,J,K))
            F=F+6
 122      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 121    CONTINUE           
 120  CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 126 H=1,N
          DO 127 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+W*Z+L-1
            GLOBM(F+4)= DBLE(-AT(H,L,I,J,K))
            F=F+5
 127      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 126    CONTINUE           
      ENDIF
*
      ENDIF
      
*
 150   CONTINUE
*       PRINT *,'SECOND GLOBM DONE'

*
*  Top face control-volumes-------------------------------
*
      IF(KE .GT. KB) THEN
*
      I=IB
      J=JB
      K=KE
      E=E+N
      DO 128 H=1,N
        DO 129 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          F=F+2
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IF(JE .GT. JB) THEN
            IA(F)=E+H-1
            JA(F)=E+W+L-1
            GLOBM(F)= DBLE(-AN(H,L,I,J,K))
            F=F+1
          ENDIF
 129    CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 128  CONTINUE         
*
      DO 160 I=IB+1,IE-1
        E=E+N
        DO 161 H=1,N
          DO 162 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E+N+L-1
            GLOBM(F+3)= DBLE(-AE(H,L,I,J,K))
            F=F+4
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
 162      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 161    CONTINUE           
 160  CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 166 H=1,N
          DO 167 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-N+L-1
            GLOBM(F+2)= DBLE(-AW(H,L,I,J,K))
            F=F+3
            IF(JE .GT. JB) THEN
              IA(F)=E+H-1
              JA(F)=E+W+L-1
              GLOBM(F)= DBLE(-AN(H,L,I,J,K))
              F=F+1
            ENDIF
 167      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 166    CONTINUE           
      ENDIF
*-------------------------------------------
      DO 180 J=JB+1,JE-1
*
      I=IB
      E=E+N
      DO 168 H=1,N
        DO 169 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          IA(F+2)=E+H-1
          JA(F+2)=E-W+L-1
          GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
          F=F+3
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
          IA(F)=E+H-1
          JA(F)=E+W+L-1
          GLOBM(F)= DBLE(-AN(H,L,I,J,K))
          F=F+1
 169    CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 168  CONTINUE           
*
      DO 170 I=IB+1,IE-1
        E=E+N
        DO 171 H=1,N
          DO 172 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+N+L-1
            GLOBM(F+4)= DBLE(-AE(H,L,I,J,K))
            IA(F+5)=E+H-1
            JA(F+5)=E+W+L-1
            GLOBM(F+5)= DBLE(-AN(H,L,I,J,K))
            F=F+6
 172      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 171    CONTINUE           
 170  CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 176 H=1,N
          DO 177 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+W+L-1
            GLOBM(F+4)= DBLE(-AN(H,L,I,J,K))
            F=F+5
 177      CONTINUE
            RHS(E+H-1)= DBLE( B(H,I,J,K))
 176    CONTINUE           
      ENDIF
*
 180   CONTINUE
*----------------------------------------------------      
      IF(JE .GT. JB) THEN
*
      J=JE
      I=IB
      E=E+N
      DO 178 H=1,N
        DO 179 L=1,N
          IA(F)=E+H-1
          JA(F) = E+L-1
          GLOBM(F)= DBLE( AP(H,L,I,J,K))
          IA(F+1)=E+H-1
          JA(F+1)=E-W*Z+L-1
          GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
          IA(F+2)=E+H-1
          JA(F+2)=E-W+L-1
          GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
          F=F+3
          IF(IE .GT. IB) THEN
            IA(F)=E+H-1
            JA(F)=E+N+L-1
            GLOBM(F)= DBLE(-AE(H,L,I,J,K))
            F=F+1
          ENDIF
 179    CONTINUE
        RHS(E+H-1)= DBLE( B(H,I,J,K))
 178  CONTINUE         
*
      DO 190 I=IB+1,IE-1
        E=E+N
        DO 191 H=1,N
          DO 192 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            IA(F+4)=E+H-1
            JA(F+4)=E+N+L-1
            GLOBM(F+4)= DBLE(-AE(H,L,I,J,K))
            F=F+5
 192      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 191    CONTINUE           
 190  CONTINUE
*
      IF(IE .GT. IB) THEN
        I=IE
        E=E+N
        DO 196 H=1,N
          DO 197 L=1,N
            IA(F  )=E+H-1
            JA(F  ) = E+L-1
            GLOBM(F)= DBLE( AP(H,L,I,J,K))
            IA(F+1)=E+H-1
            JA(F+1)=E-W*Z+L-1
            GLOBM(F+1)= DBLE(-AB(H,L,I,J,K))
            IA(F+2)=E+H-1
            JA(F+2)=E-W+L-1
            GLOBM(F+2)= DBLE(-AS(H,L,I,J,K))
            IA(F+3)=E+H-1
            JA(F+3)=E-N+L-1
            GLOBM(F+3)= DBLE(-AW(H,L,I,J,K))
            F=F+4
 197      CONTINUE
          RHS(E+H-1)= DBLE( B(H,I,J,K))
 196    CONTINUE           
      ENDIF
*
      ENDIF
*      PRINT *,'THIRD GLOBM DONE'
*
      ENDIF
*
*  Indicate the number of non-zero elements
*      
        NELT=F-1
*
      RETURN
      END
*
