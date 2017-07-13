*     This file contains 2 subroutines: COEFCN and MASFLX
*
************************************************************************
*
      SUBROUTINE COEFCN(ACUW,ACUE,ACVS,ACVN,ACWB,ACWT,BC,
     C                  RHO,AREP,ARNP,ARTP,CVTYPE,
     C                  IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the coefficients of the continuity
*     equation for each control volume.
*
*     ACUW(ID,JD,KD) coefficient for west face uhat; output
*     ACUE(ID,JD,KD) coefficient for east face uhat; output
*     ACVS(ID,JD,KD) coefficient for south face vhat; output
*     ACVN(ID,JD,KD) coefficient for north face vhat; output
*     ACWB(ID,JD,KD) coefficient for bottom face what; output
*     ACWT(ID,JD,KD) coefficient for top face what; output
*     BC(ID,JD,KD)   mass source term; output
*     RHO fluid density (kg/m^3); input
*     AREP(JD,KD) c.v. area of face at e point; input
*     ARNP(ID,KD) c.v. area of face at n point; input
*     ARTP(ID,JD) c.v. area of face at t point; input
*     CVTYPE(ID,JD,KD,NNB+1) indicates type of each CV; input
*     IB,IE first and last interior indices in i; input
*     JB,JE first and last interior indices in j; input
*     KB,KE first and last interior indices in k; input
*     ID,JD,KD array dimensions; input     
*     NNB number of neighbors for interior nodes; input
*
*     Notes: 1) The form of the continuity equation is:
*     ACUE*UHE(I,J,K) + ACUW*UHE(I-1,J,K) + ACVN*VHN(I,J,K) 
*     + ACVS*VHN(I,J-1,K) + ACWT*WHT(I,J,K) + ACWB*WHT(I,J,K-1)
*     + BC = 0
*
***********************************************************************
      IMPLICIT NONE
      REAL*8 ACUW(ID,JD,KD),ACUE(ID,JD,KD),ACVS(ID,JD,KD)
      REAL*8 ACVN(ID,JD,KD),ACWB(ID,JD,KD),ACWT(ID,JD,KD)
      REAL*8 BC(ID,JD,KD),AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD)
      REAL*8 RHO
      INTEGER CVTYPE(ID,JD,KD,NNB+1),IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB
      INTEGER I,J,K
*     
      CALL NULLM(ACUW, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(ACUE, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(ACVS, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(ACVN, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(ACWB, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(ACWT, IB,IE,JB,JE,KB,KE,ID,JD,KD)
      CALL NULLM(BC, IB,IE,JB,JE,KB,KE,ID,JD,KD)            
*     
      DO 30 K=KB,KE
       DO 20 J=JB,JE
        DO 10 I=IB,IE
         IF(CVTYPE(I,J,K,1).NE.3) THEN
          IF(CVTYPE(I,J,K,3).NE.3) THEN
           ACUE(I,J,K)=0.0+RHO*AREP(J,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,2).NE.3) THEN
           ACUW(I,J,K)=0.0-RHO*AREP(J,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,5).NE.3) THEN
           ACVN(I,J,K)=0.0+RHO*ARNP(I,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,4).NE.3) THEN
           ACVS(I,J,K)=0.0-RHO*ARNP(I,K)
          ENDIF
*
          IF(CVTYPE(I,J,K,7).NE.3) THEN
           ACWT(I,J,K)=0.0+RHO*ARTP(I,J)
          ENDIF
*
          IF(CVTYPE(I,J,K,6).NE.3) THEN
           ACWB(I,J,K)=0.0-RHO*ARTP(I,J)
          ENDIF
*         
          BC(I,J,K)=0.0
         ENDIF
 10     CONTINUE
 20    CONTINUE
 30   CONTINUE
*
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE MASFLX(ME,MN,MT, UHE,VHN,WHT,ACUW,ACUE,ACVS,ACVN,
     C                  ACWB,ACWT,IB,IE,JB,JE,KB,KE,ID,JD,KD)
*
*     Subroutine to calculate the mass flow for the 
*     normal advection fluxes through the faces of
*     each control volume. 
*
*     ME(ID,JD,KD) normal mass flux for east face; output
*     MN(ID,JD,KD) normal mass flux for north face; output
*     ME(ID,JD,KD) normal mass flux for top face; output
*
*     UHE(ID,JD,KD) x component of mass velocity at e point; input
*     VHN(ID,JD,KD) y component of mass velocity at n point; input
*     WHT(ID,JD,KD) z component of mass velocity at t point; input
*     ACUW(ID,JD,KD) coefficient for west face uhat; input
*     ACUE(ID,JD,KD) coefficient for east face uhat; input
*     ACVS(ID,JD,KD) coefficient for south face vhat; input
*     ACVN(ID,JD,KD) coefficient for north face vhat; input
*     ACWB(ID,JD,KD) coefficient for bottom face what; input
*     ACWT(ID,JD,KD) coefficient for top face what; input
*     IB,IE first and last interior indices in i; input
*     JB,JE first and last interior indices in j; input
*     KB,KE first and last interior indices in k; input
*     ID,JD,KD array dimensions; input     
*
***********************************************************************
      IMPLICIT NONE
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 UHE(ID,JD,KD),VHN(ID,JD,KD),WHT(ID,JD,KD)
      REAL*8 ACUW(ID,JD,KD),ACUE(ID,JD,KD)
      REAL*8 ACVS(ID,JD,KD),ACVN(ID,JD,KD)
      REAL*8 ACWB(ID,JD,KD),ACWT(ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,I,J,K,ID,JD,KD
*     
*  IB - 1 faces
*
      I = IB - 1
      DO 5 J=JB,JE
       DO 4 K=KB,KE
        ME(I,J,K)=0.0-ACUW(I+1,J,K)*UHE(I,J,K)
  4    CONTINUE
  5   CONTINUE
*     
*  JB - 1 faces
*
      J = JB - 1
      DO 10 I=IB,IE
       DO 9 K=KB,KE
        MN(I,J,K)=0.0-ACVS(I,J+1,K)*VHN(I,J,K)
  9    CONTINUE
 10   CONTINUE
*     
*  KB - 1 faces
*
      K = KB - 1
      DO 15 I=IB,IE
       DO 14 J=JB,JE
        MT(I,J,K)=0.0-ACWB(I,J,K+1)*WHT(I,J,K)
 14    CONTINUE
 15   CONTINUE
*
*  Interior faces
*
      DO 40 K=KB,KE
       DO 30 J=JB,JE
        DO 20 I=IB,IE
         ME(I,J,K)= ACUE(I,J,K)*UHE(I,J,K)
         MN(I,J,K)= ACVN(I,J,K)*VHN(I,J,K)
         MT(I,J,K)= ACWT(I,J,K)*WHT(I,J,K)
 20     CONTINUE
 30    CONTINUE
 40   CONTINUE
*
      RETURN
      END

