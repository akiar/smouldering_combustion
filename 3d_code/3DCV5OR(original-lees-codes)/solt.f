*        file solt.f 
************************************************************************
*
      SUBROUTINE SOLT(TF,TS,WORK5,WORK6,RLX,
     C                ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C                IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*  Routine to organize the solution of T
*
************************************************************************
*  
      IMPLICIT NONE
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
*
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD)
      REAL*8 ATS(NT,NT,ID,JD,KD),ATN(NT,NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)      
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
*
      REAL*8 WORK5(NT,NT),WORK6(NT),RLX
*      
      INTEGER IB,IE,JB,JE,KB,KE,NT,ID,JD,KD
*
      INTEGER neqmax
      PARAMETER(neqmax=4000000)
      DOUBLE PRECISION SOL(neqmax)      
*
*  Absorb boundary conditions
*
      CALL BCELIM(ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,WORK5,WORK6,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
*  Construct the global coeff. matrix and RHS vector, and solve 
*  the system using WATSIT-B.  Solution is returned in 
*  SOL(neqmax).
*
      CALL WATSOLN(SOL, ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C             IB,IE,JB,JE,KB,KE,ID,JD,KD,NT)
*
      CALL ASSGNT(TF,TS, SOL,
     C            IB,IE,JB,KB,KE,JE,NT,ID,JD,KD,neqmax)
*
*  Update the boundary values
*
      CALL BCCLCT(TF,TS, ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C            IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
*
      RETURN
      END