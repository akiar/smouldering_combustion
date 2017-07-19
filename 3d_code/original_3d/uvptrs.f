*
*        file uvptrs.f
************************************************************************
*
      SUBROUTINE UVPTRS(TF,TS,P,U,V,W,
     C            NORM,RSD,ARSD,RSDMAX,RSDOLD,I_MAX,J_MAX,K_MAX,
     C            AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,
     C            CFLUID,CSOLID,CVTYPE,ITERMO,ISOTHM,
     C            IB,IE,JB,JE,KB,KE,N,NT,ID,JD,KD,NNB)
*
*  Routine to compute the residuals for the T,P,U,V,W equations
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 AUW(N,N,ID,JD,KD),AUE(N,N,ID,JD,KD)
      REAL*8 AUS(N,N,ID,JD,KD),AUN(N,N,ID,JD,KD)
      REAL*8 AUB(N,N,ID,JD,KD),AUT(N,N,ID,JD,KD)      
      REAL*8 AUP(N,N,ID,JD,KD),BU(N,ID,JD,KD)
      REAL*8 ATW(NT,NT,ID,JD,KD),ATE(NT,NT,ID,JD,KD)
      REAL*8 ATS(NT,NT,ID,JD,KD),ATN(NT,NT,ID,JD,KD)
      REAL*8 ATB(NT,NT,ID,JD,KD),ATT(NT,NT,ID,JD,KD)      
      REAL*8 ATP(NT,NT,ID,JD,KD),BT(NT,ID,JD,KD)
      REAL*8 NORM(6),RSD(ID,JD,KD),ARSD(6),RSDMAX(6),RSDOLD(6),RATE
      REAL*8 CFLUID(ID,JD,KD),CSOLID(ID,JD,KD)
      INTEGER ITERMO,ISOTHM,I_MAX,J_MAX,K_MAX
      INTEGER IB,IE,JB,JE,KB,KE,N,NT,ID,JD,KD,NNB
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
*
*  Find magnitude of average value in each solution field
*
      IF(ISOTHM.GT.0) THEN
       CALL MEANA(NORM(1), 
     C            TF,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
       CALL MEANA(NORM(2), 
     C            TS,0,1,1,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)       
      ENDIF
      NORM(3)= 1.0
      CALL MEANA(NORM(4), 
     C           U,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      CALL MEANA(NORM(5), 
     C           V,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      CALL MEANA(NORM(6), 
     C           W,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
*
*  Calculate and print residuals for each field
*
      IF(ISOTHM.GT.0) THEN
*      
      CALL RESIDT(RSD,ARSD(1),RSDMAX(1),I_MAX,J_MAX,K_MAX,
     C            TF,TS,ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,CFLUID,
     C            1,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
      RATE= (ARSD(1)/NORM(1))/RSDOLD(1)
      WRITE(ITERMO,7221) 'Tf-enrgy  ',RATE,ARSD(1)/NORM(1),
     C                  RSDMAX(1)/NORM(1),I_MAX,J_MAX,K_MAX      
*
      CALL RESIDT(RSD,ARSD(2),RSDMAX(2),I_MAX,J_MAX,K_MAX,
     C            TF,TS,ATP,ATW,ATE,ATS,ATN,ATB,ATT,BT,CSOLID,
     C            2,IB,IE,JB,JE,KB,KE,NT,ID,JD,KD)
      RATE= (ARSD(2)/NORM(2))/RSDOLD(2)
      WRITE(ITERMO,7221) 'Ts-enrgy  ',RATE,ARSD(2)/NORM(2),
     C                  RSDMAX(2)/NORM(2),I_MAX,J_MAX,K_MAX  
*          
      ENDIF
*
      CALL RESIDM(RSD,ARSD(3),RSDMAX(3),I_MAX,J_MAX,K_MAX,
     C            P,U,V,W,AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            1,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      RATE= (ARSD(3)/NORM(3))/RSDOLD(3)
      WRITE(ITERMO,7221) 'P-mass    ',RATE,ARSD(3)/NORM(3),
     C                  RSDMAX(3)/NORM(3),I_MAX,J_MAX,K_MAX
*
      CALL RESIDM(RSD,ARSD(4),RSDMAX(4),I_MAX,J_MAX,K_MAX,
     C            P,U,V,W,AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            2,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      RATE= (ARSD(4)/NORM(4))/RSDOLD(4)
      WRITE(ITERMO,7221) 'U-mom     ',RATE,ARSD(4)/NORM(4),
     C                  RSDMAX(4)/NORM(4),I_MAX,J_MAX,K_MAX
*
      CALL RESIDM(RSD,ARSD(5),RSDMAX(5),I_MAX,J_MAX,K_MAX,
     C            P,U,V,W,AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            3,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      RATE= (ARSD(5)/NORM(5))/RSDOLD(5)
      WRITE(ITERMO,7221) 'V-mom     ',RATE,ARSD(5)/NORM(5),
     C                  RSDMAX(5)/NORM(5),I_MAX,J_MAX,K_MAX
*
      CALL RESIDM(RSD,ARSD(6),RSDMAX(6),I_MAX,J_MAX,K_MAX,
     C            P,U,V,W,AUP,AUW,AUE,AUS,AUN,AUB,AUT,BU,
     C            4,IB,IE,JB,JE,KB,KE,N,ID,JD,KD)
      RATE= (ARSD(6)/NORM(6))/RSDOLD(6)
      WRITE(ITERMO,7221) 'W-mom     ',RATE,ARSD(6)/NORM(6),
     C                  RSDMAX(6)/NORM(6),I_MAX,J_MAX,K_MAX     
      WRITE(ITERMO,7222)
*
 7221 FORMAT(
     C'|',A10,'|',F6.2,' |'1PE9.1'|'1PE9.1'| (',I3,',',I3,',',I3,')|')
 7222 FORMAT(
     C'+-----------------------------------------------------+',/)
*
      RETURN
      END
