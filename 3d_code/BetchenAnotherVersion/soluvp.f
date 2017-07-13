*       file soluvp.f for WATsitB
************************************************************************
*
      SUBROUTINE SOLUVP(P,U,V,W,UHE,VHN,WHT,RLX,
     C        AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,WORK3,WORK4,
     C        IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*  Routine to organize the solution of P,U,V,W
*
************************************************************************
*  
      IMPLICIT NONE
*
      REAL*8 AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD),AUS(N,N,ID,JD,KD)
      REAL*8 AUN(N,N,ID,JD,KD),AUP(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 AUB(N,N,ID,JD,KD),AUT(N,N,ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD),P(ID,JD,KD)
      REAL*8 UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD),RLX
      REAL*8 WORK3(N,N),WORK4(N)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,N,MR,MG
*
*  Declaration of storage arrays for solver; DLAP Triad format is used;
*  MR is the maximum number of equations, and MG should be 5*MR
*
      PARAMETER(MR=5000000)
      PARAMETER(MG=5*MR)
      DOUBLE PRECISION X(MR),B(MR)
      DOUBLE PRECISION A(MG)
      INTEGER IA(MG),JA(MG)
      INTEGER NEQ,NELT,ISYM,ITOL,ITMAX
      INTEGER ITER,IERR,IUNIT,LENW,LENIW,IWORK(MG+4*MR)
      DOUBLE PRECISION ERR,TOL,RWORK(MG+8*MR)  
*
*
**********************************
*  Executable statements
**********************************
*
*  Absorb boundary conditions to reduce the size of the global 
*  coefficient matrix
*
      CALL BCELIM(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,WORK3,WORK4,
     C            IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*  Calculate order of matrix for solver routines
*
      NEQ = N*(IE-IB+1)*(JE-JB+1)*(KE-KB+1)
*
*  Convert the coeff matrix to DLAP Triad format storing in GLOBM(NELT),
*  IA(NELT),JA(NELT), and store the B vector in RHS(NEQ)
*
*       PRINT *,'INSIDE SOLUVP'
      CALL GLOBMS(A,B,IA,JA,
     C            AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            IB,IE,JB,JE,KB,KE,N,ID,JD,KD,
     C            MG,MR,NELT,NEQ)
*       PRINT *,'OUT OF GLOBM'
*       PRINT*,'A', A(1)
*       PRINT*,'B', B(1)
*       PRINT*,'IA', IA(1)
*       PRINT*,'JA', JA(1)
*       PRINT*,'AP11', AUP(1,1,2,2,2)
*
*  Provide initial guess for DSDBCG solver. Current solution is passed 
*  as an initial guess in X(); solution to system is then also 
*  returned in X().
*
      CALL SGUESS(X, P,U,V,W,IB,IE,JB,JE,KB,KE,N,ID,JD,KD,MR)
*
*  Solve system using Incomplete LU BiConjugate Gradient Sparse Ax=b 
*  solver DSLUBC. Routine to solve a linear system  Ax = b  using the
*  BiConjugate Gradient  method  with  Incomplete LU decomposition 
*  preconditioning.

*
      ISYM= 0
      ITOL= 1
      TOL=  1.0D-12
      ITMAX=500
      IUNIT=73
      LENW= NELT+NEQ+8*NEQ
      LENIW=NELT+NEQ+4*NEQ+12
               PRINT *,'NELT IS', NELT
*
      CALL DSLUBC(NEQ,B,X,NELT,IA,JA,A,ISYM,ITOL,TOL,
     C            ITMAX,ITER,ERR,IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
*
*         PRINT *,'SOLUTION SOLVED'
*--Insert the new solution into the P(),U(),V() and W() array
*       
      CALL ASSGNS(P,U,V,W, X,
     C            IB,IE,JB,JE,KB,KE,N,ID,JD,KD,MR)
*
*  Update the boundary values
*
      CALL BCCLCU(P, AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            1,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      CALL BCCLCU(U, AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            2,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      CALL BCCLCU(V, AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            3,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      CALL BCCLCU(W, AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            4,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
      RETURN
      END

