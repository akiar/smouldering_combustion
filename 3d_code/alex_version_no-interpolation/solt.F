*       file solt.f for DLAP
************************************************************************
*
      SUBROUTINE SOLT(TF,TS,WORK5,WORK6,RLX,
     C        ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C        IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*  Routine to organize the solution of TF,TS
*
************************************************************************
*  
      IMPLICIT NONE
*
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD),ATS(NT,NT,ID,JD,KD)
      REAL*8 ATN(NT,NT,ID,JD,KD),ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 WORK5(NT,NT),WORK6(NT),RLX
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NT,MR,MG
*
*  Declaration of storage arrays for solver; DLAP Triad format is used;
*  MR is the maximum number of equations, and MG should be 5*MR
*
      PARAMETER(MR=1250000)
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
      CALL BCELIM(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,WORK5,WORK6,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*  Calculate order of matrix for solver routines
*
      NEQ = NT*(IE-IB+1)*(JE-JB+1)*(KE-KB+1)
*
*  Convert the coeff matrix to DLAP Triad format storing in GLOBM(NELT),
*  IA(NELT),JA(NELT), and store the B vector in RHS(NEQ)
*
*       PRINT *,'INSIDE SOLUVP'
      CALL GLOBMS(A,B,IA,JA,
     C            ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,
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
      CALL SGUESST(X,TF,TS,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,MR)
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
*--Insert the new solution into the TF() and TS() array
*       
      CALL ASSGNT(TF,TS,X,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD,MR)
*
*  Update the boundary values
*
      CALL BCCLCT(TF,TS, ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
      RETURN
      END

