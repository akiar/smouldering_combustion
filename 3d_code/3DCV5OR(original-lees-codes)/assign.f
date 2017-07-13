*     
*     This file contains 2 subroutines: ASSGNT and ASSGNS
*
********************************************************************
      SUBROUTINE ASSGNT(TF,TS, SOLN,
     C                  IB,IE,JB,KB,KE,JE,NT,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of TF,TS from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      IMPLICIT NONE
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,NT,M,ID,JD,KD
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
          TF(I,J,K)= SOLN(E  )
          TS(I,J,K)= SOLN(E+1)
          E=E+NT
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END      
*      
********************************************************************
      SUBROUTINE ASSGNS(P,U,V,W, SOLN,
     C                  IB,IE,JB,JE,KB,KE,N,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of P,U,V,W from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      IMPLICIT NONE
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,N,M,ID,JD,KD
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 2 I=IB,IE
          P(I,J,K)= SOLN(E  )
          U(I,J,K)= SOLN(E+1)
          V(I,J,K)= SOLN(E+2)
          W(I,J,K)= SOLN(E+3)
          E=E+N
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END
