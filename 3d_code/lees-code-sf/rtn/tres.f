*     file tres.f
***********************************************************************
*
      SUBROUTINE TREST(TRMAX, TF,TS,TFOLD,TSOLD,CVTYPE,
     C                 CI,CII,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB) 
*
*     Routine to evaluate temporal convergence of the energy equations.
*
*     TRMAX(6)  maximum percent change in each distribution over time 
*               step; output
*     TF,TS  present values of dependent variables; input
*     TFOLD,TSOLD  old values of dependent variables; input
*     C*  = 1, *th variable is checked for temporal convergence; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  grid-array dimensions; input 
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 TRMAX(6),DIFF,MEANTF,MEANTS
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 TFOLD(ID,JD,KD),TSOLD(ID,JD,KD)    
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CI,CII,I,J,K  
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
*
      CALL MEANA(MEANTF, 
     C           TF,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      CALL MEANA(MEANTS,
     C           TS,0,1,1,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
         CALL PDIFF(DIFF, TF(I,J,K),TFOLD(I,J,K),MEANTF) 
         IF(DIFF.GT.TRMAX(1) .AND. CI.EQ.1) THEN
          TRMAX(1) = DIFF
         ENDIF
*
         CALL PDIFF(DIFF, TS(I,J,K),TSOLD(I,J,K),MEANTS)
         IF(DIFF.GT.TRMAX(2) .AND. CII.EQ.1) THEN
          TRMAX(2) = DIFF
         ENDIF         
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE 
*
      RETURN
      END
*
***********************************************************************
*
      SUBROUTINE TRESU(TRMAX, P,U,V,W,POLD,UOLD,VOLD,WOLD,CVTYPE,
     C           CIII,CIV,CV,CVI,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB) 
*
*     Routine to evaluate temporal convergence of the mass and
*     momentum equations.
*
*     TRMAX(6)  maximum percent change in each distribution over time 
*               step; output
*     P,U,V,W  present values of dependent variables; input
*     POLD,UOLD,VOLD,WOLD  old values of dependent variables; input
*     C*  = 1, *th variable is checked for temporal convergence; input
*     IB,IE  first and last interior indices in i; input
*     JB,JE  first and last interior indices in j; input
*     KB,KE  first and last interior indices in k; input
*     ID,JD,KD  grid-array dimensions; input 
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 TRMAX(6),DIFF,MEANP,MEANU,MEANV,MEANW
      REAL*8 P(ID,JD,KD),U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 POLD(ID,JD,KD),UOLD(ID,JD,KD),VOLD(ID,JD,KD)
      REAL*8 WOLD(ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CIII,CIV,CV,CVI,I,J,K  
      INTEGER CVTYPE(ID,JD,KD,NNB+1)
*
      CALL MEANA(MEANP, 
     C           P,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
      CALL MEANA(MEANU, 
     C           U,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
      CALL MEANA(MEANV,
     C           V,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)      
      CALL MEANA(MEANW,
     C           W,1,1,0,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)                  
*      
      DO 30 I=IB,IE
       DO 20 J=JB,JE
        DO 10 K=KB,KE
         CALL PDIFF(DIFF, P(I,J,K),POLD(I,J,K),MEANP) 
         IF(DIFF.GT.TRMAX(3) .AND. CIII.EQ.1) THEN
          TRMAX(3) = DIFF
         ENDIF
*
         CALL PDIFF(DIFF, U(I,J,K),UOLD(I,J,K),MEANU)
         IF(DIFF.GT.TRMAX(4) .AND. CIV.EQ.1) THEN
          TRMAX(4) = DIFF
         ENDIF         
*
         CALL PDIFF(DIFF, V(I,J,K),VOLD(I,J,K),MEANV)
         IF(DIFF.GT.TRMAX(5) .AND. CV.EQ.1) THEN
          TRMAX(5) = DIFF
         ENDIF    
*
         CALL PDIFF(DIFF, W(I,J,K),WOLD(I,J,K),MEANW)
         IF(DIFF.GT.TRMAX(6) .AND. CVI.EQ.1) THEN
          TRMAX(6) = DIFF
         ENDIF                      
   10   CONTINUE
   20  CONTINUE
   30 CONTINUE
*
      RETURN
      END      
*
***********************************************************************
*
      SUBROUTINE PDIFF(DIFF, A,AOLD,MEANA)
*
*     Routine to evaluate percent difference, normalized to mean 
*     value of variable
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 DIFF,A,AOLD,MEANA
*
      DIFF = ABS((A-AOLD)/(MEANA+1.0E-20))
*
      RETURN
      END            