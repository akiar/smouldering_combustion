*     
*     This file contains ASSIGN AND GUESS subroutines
*
********************************************************************
      SUBROUTINE ASSGNT(TF,TS, SOLN,
     C                  IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of TF,TS from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      IMPLICIT NONE
      INTEGER IB,IE,JB,JE,KB,KE,NT,M,ID,JD,KD
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 SOLN(M)
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
      INTEGER IB,IE,JB,JE,KB,KE,N,M,ID,JD,KD
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
          P(I,J,K)= REAL(SOLN(E  ))
          U(I,J,K)= REAL(SOLN(E+1))
          V(I,J,K)= REAL(SOLN(E+2))
          W(I,J,K)= REAL(SOLN(E+3))
          E=E+N
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END
********************************************************************
      SUBROUTINE ASSGST(SOLN,TF,TS, 
     C                  IB,IE,JB,KB,KE,JE,NT,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of TF,TS to the solution
*     vector for PETSC* (assign solution vector for T)
********************************************************************
*
      IMPLICIT NONE
      INTEGER IB,IE,JB,JE,KB,KE,NT,M,ID,JD,KD
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 SOLN(M)
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
          SOLN(E  ) = TF(I,J,K)
          SOLN(E+1) = TS(I,J,K)
          E=E+NT
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END      
*      
********************************************************************
      SUBROUTINE ASSGSU(P,U,V,W, SOLN,
     C                  IB,IE,JB,JE,KB,KE,N,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of P,U,V,W to the 
*     solution vector for PETSC* (assign solution vector 
*     for P,U,V,W)
********************************************************************
*
      IMPLICIT NONE
      INTEGER IB,IE,JB,JE,KB,KE,N,ID,JD,KD,M
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 SOLN(M)
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
          SOLN(E  )= P(I,J,K)
          SOLN(E+1)= U(I,J,K)
          SOLN(E+2)= V(I,J,K)
          SOLN(E+3)= W(I,J,K) 
          E=E+N
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END
*******************************************************************
      SUBROUTINE ASSGNI(PHI, SOLN,IB,IE,JB,JE,KB,KE,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of PHI(ID,JD,KD) from  
*     the solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      IMPLICIT NONE
      REAL*8 PHI(ID,JD,KD)
      REAL*8 SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,I,J,K,E,M,ID,JD,KD
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
         PHI(I,J,K)= SOLN(E)
         E=E+1
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE ASSGNY(SOLN,PHI,IB,IE,JB,JE,KB,KE,ID,JD,KD,M)
*
*     Subroutine to assign Solution vector from prev solution for
*     solver
*
********************************************************************
*
      IMPLICIT NONE
      REAL*8 PHI(ID,JD,KD)
      REAL*8 SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,I,J,K,E,M,ID,JD,KD
*
      E=1
      DO 3 K=KB,KE
       DO 2 J=JB,JE
        DO 1 I=IB,IE
          SOLN(E)=PHI(I,J,K)
          E=E+1
  1     CONTINUE
  2    CONTINUE
  3   CONTINUE
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE IGUESS(SOLN, PHI,IB,IE,JB,JE,KB,KE,ID,JD,KD,M)
*
*     Subroutine to insert the nodal values of PHI(ID,JD,KD) into
*     solution vector SOLN(M) as an initial guess to solver.
*
********************************************************************
*
      REAL*8 PHI(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,I,J,K,E,M,ID,JD,KD
*
      E=0
      DO 3 K=KB,KE
        DO 2 J=JB,JE
          DO 1 I=IB,IE
            E=E+1
            SOLN(E)= DBLE(PHI(I,J,K))
  1       CONTINUE
  2     CONTINUE
  3   CONTINUE          
*
      RETURN
      END
*
********************************************************************
      SUBROUTINE SGUESS(SOLN,P,U,V,W,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of U,V,W,P from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,N,M,ID,JD,KD
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
        DO 2 J=JB,JE
          DO 1 I=IB,IE
            SOLN(E  )= DBLE(P(I,J,K))
            SOLN(E+1)= DBLE(U(I,J,K))
            SOLN(E+2)= DBLE(V(I,J,K))
            SOLN(E+3)= DBLE(W(I,J,K))
            E=E+N
  1       CONTINUE
  2     CONTINUE
  3   CONTINUE          
*
      RETURN
      END     
********************************************************************
      SUBROUTINE SGUESST(SOLN,TF,TS,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,M)
*
*     Subroutine to assign the nodal values of TF,TS from the 
*     solution vector of the direct solution SOLN(M).
*
********************************************************************
*
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      DOUBLE PRECISION SOLN(M)
      INTEGER IB,IE,JB,JE,KB,KE,N,M,ID,JD,KD
      INTEGER I,J,K,E
*
      E=1
      DO 3 K=KB,KE
        DO 2 J=JB,JE
          DO 1 I=IB,IE
            SOLN(E  )= DBLE(TF(I,J,K))
            SOLN(E+1)= DBLE(TS(I,J,K))
            E=E+N
  1       CONTINUE
  2     CONTINUE
  3   CONTINUE          
*
      RETURN
      END         
