*       file soluvp.f for WATsitB
************************************************************************
*
      SUBROUTINE SOLUVP(P,U,V,W,UHE,VHN,WHT,RLX,
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
*   
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD),P(ID,JD,KD)
      REAL*8 UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD)
*
      REAL*8 WORK3(N,N),WORK4(N),RLX
*
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,N
*
      INTEGER neqmax
      PARAMETER(neqmax=4000000)
      DOUBLE PRECISION SOL(neqmax)
*
*
**********************************
*  Executable statements
**********************************
*
*  Absorb boundary conditions
*
      CALL BCELIM(AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,WORK3,WORK4,
     C            IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
*
*  Construct the global coeff. matrix and RHS vector, and solve 
*  the system using WATSIT-B.  Solution is returned in 
*  SOL(neqmax).
*
      CALL WATSOLN(SOL, AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C             IB,IE,JB,JE,KB,KE,ID,JD,KD,N)
*
      CALL ASSGNS(P,U,V,W, SOL,
     C            IB,IE,JB,JE,KB,KE,N,ID,JD,KD,neqmax)
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

