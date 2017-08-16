*
*     This file contains subroutine SGENC
************************************************************************
*
      SUBROUTINE SGENC(SGENT,SGEN,SGENF,SGENS,RATIOT,RATIOS,RATIO,
     C                 U,V,W,TF,TS,
     C                 DE,DN,DT,DEF,DNF,DTF,DES,DNS,DTS,
     C                 PRSTY,VISC,RHO,KPERM,CFORCH,MUEFF,
     C                 KEFF,KEFFX,KEFFY,KEFFZ,HSF,SPECSA,
     C                 DIEP,DJNP,DKTP,DISE,DISN,DIST,
     C                 AREP,ARNP,ARTP,VOLP,
     C                 CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the local and total rates of entropy 
*     generation.
*
*     SGENT   total rate of entropy generation; output    
*     SGEN    local volumetric rate of entropy generation; output    
*     SGENF   local volumetric rate of entropy generation in fluid 
*             constituent; output    
*     SGENS   local volumetric rate of entropy generation in solid 
*             constituent; output    
*     RATIOT  Bejan number; output
*     RATIOS  ratio of entropy generation by solid to total; output
*     RATIO   local ratio of entropy generation by fluid flow to 
*             entropy generation by heat transfer; output
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 SGENT,SGEN(ID,JD,KD),SGENF(ID,JD,KD),SGENS(ID,JD,KD)
      REAL*8 RATIOT,RATIOS,RATIO(ID,JD,KD)
      REAL*8 U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD)
      REAL*8 TF(ID,JD,KD),TS(ID,JD,KD)
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      REAL*8 DEF(ID,JD,KD),DNF(ID,JD,KD),DTF(ID,JD,KD)
      REAL*8 DES(ID,JD,KD),DNS(ID,JD,KD),DTS(ID,JD,KD)
      REAL*8 PRSTY(ID,JD,KD)
      REAL*8 KPERM(ID,JD,KD),CFORCH(ID,JD,KD)
      REAL*8 MUEFF(ID,JD,KD),KEFF(ID,JD,KD)
      REAL*8 KEFFX(ID,JD,KD),KEFFY(ID,JD,KD),KEFFZ(ID,JD,KD)
      REAL*8 HSF(ID,JD,KD),SPECSA(ID,JD,KD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DISE(ID),DISN(JD),DIST(KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD),VOLP(ID,JD,KD)
      REAL*8 DTFDX,DTFDY,DTFDZ,DTSDX,DTSDY,DTSDZ
      REAL*8 DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ
      REAL*8 VEL,DISS,RATIOB
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
*
      REAL*8 VISC(ID,JD,KD),RHO(ID,JD,KD)
*
      SGENT = 0.0
      RATIOT = 0.0
      RATIOB = 0.0
      RATIOS = 0.0
*
      DO 20 I=IB,IE
       DO 10 J=JB,JE
        DO 5 K=KB,KE
*
*        Cases for x derivatives:
*
*        Interface at east face  
*
         IF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,3)
     C      .AND. CVTYPE(I,J,K,3).NE.1) THEN
          DTFDX = 0.5*(TF(I,J,K)-TF(I-1,J,K))/DIEP(I-1)
          DTSDX = 0.5*(TS(I,J,K)-TS(I-1,J,K))/DIEP(I-1)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,3).EQ.0) THEN
            DTFDX = DTFDX+0.5*DEF(I,J,K)*(TF(I+1,J,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*AREP(J,K))
            DTSDX = DTSDX+0.5*DES(I,J,K)*(TF(I+1,J,K)
     C              -TS(I,J,K))/(KEFFX(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,3).EQ.3) THEN
            DTFDX = DTFDX+0.5*DEF(I,J,K)*(TS(I+1,J,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*AREP(J,K))
            DTSDX = DTSDX+0.5*DES(I,J,K)*(TS(I+1,J,K)
     C              -TS(I,J,K))/(KEFFX(I,J,K)*AREP(J,K))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,3).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDX = DTFDX+0.5*(DEF(I,J,K)*(TF(I+1,J,K)
     C              -TF(I,J,K))+DES(I,J,K)*(TS(I+1,J,K)
     C              -TF(I,J,K)))/(KEFF(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDX = DTFDX+0.5*(DEF(I,J,K)*(TF(I+1,J,K)
     C              -TS(I,J,K))+DES(I,J,K)*(TS(I+1,J,K)
     C              -TS(I,J,K)))/(KEFF(I,J,K)*AREP(J,K))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDX = DTFDX+0.5*DEF(I,J,K)*(TS(I+1,J,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDX = DTSDX+0.5*DES(I,J,K)*(TF(I+1,J,K)
     C              -TS(I,J,K))/(KEFF(I,J,K)*AREP(J,K))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDX = 0.5*((U(I,J,K)-U(I-1,J,K))/DIEP(I-1)
     C            +DE(I,J,K)*(U(I+1,J,K)-U(I,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
           DVDX = 0.5*((V(I,J,K)-V(I-1,J,K))/DIEP(I-1)
     C            +DE(I,J,K)*(V(I+1,J,K)-V(I,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
           DWDX = 0.5*((W(I,J,K)-W(I-1,J,K))/DIEP(I-1)
     C            +DE(I,J,K)*(W(I+1,J,K)-W(I,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
          ENDIF
*
*        Interface at west face  
*
         ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,2)
     C          .AND. CVTYPE(I,J,K,2).NE.1) THEN
          DTFDX = 0.5*(TF(I+1,J,K)-TF(I,J,K))/DIEP(I)
          DTSDX = 0.5*(TS(I+1,J,K)-TS(I,J,K))/DIEP(I)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,2).EQ.0) THEN
            DTFDX = DTFDX+0.5*DEF(I-1,J,K)*(TF(I,J,K)
     C              -TF(I-1,J,K))/(KEFF(I,J,K)*AREP(J,K))
            DTSDX = DTSDX+0.5*DES(I-1,J,K)*(TS(I,J,K)
     C              -TF(I-1,J,K))/(KEFFX(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,2).EQ.3) THEN
            DTFDX = DTFDX+0.5*DEF(I-1,J,K)*(TF(I,J,K)
     C              -TS(I-1,J,K))/(KEFF(I,J,K)*AREP(J,K))
            DTSDX = DTSDX+0.5*DES(I-1,J,K)*(TS(I,J,K)
     C              -TS(I-1,J,K))/(KEFFX(I,J,K)*AREP(J,K))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,2).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDX = DTFDX+0.5*(DEF(I-1,J,K)*(TF(I,J,K)
     C              -TF(I-1,J,K))+DES(I-1,J,K)*(TF(I,J,K)
     C              -TS(I-1,J,K)))/(KEFF(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDX = DTFDX+0.5*(DEF(I-1,J,K)*(TS(I,J,K)
     C              -TF(I-1,J,K))+DES(I-1,J,K)*(TS(I,J,K)
     C              -TS(I-1,J,K)))/(KEFF(I,J,K)*AREP(J,K))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDX = DTFDX+0.5*DEF(I-1,J,K)*(TF(I,J,K)
     C              -TS(I-1,J,K))/(KEFF(I,J,K)*AREP(J,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDX = DTSDX+0.5*DES(I-1,J,K)*(TS(I,J,K)
     C              -TF(I-1,J,K))/(KEFF(I,J,K)*AREP(J,K))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDX = 0.5*((U(I+1,J,K)-U(I,J,K))/DIEP(I)
     C            +DE(I-1,J,K)*(U(I,J,K)-U(I-1,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
           DVDX = 0.5*((V(I+1,J,K)-V(I,J,K))/DIEP(I)
     C            +DE(I-1,J,K)*(V(I,J,K)-V(I-1,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
           DWDX = 0.5*((W(I+1,J,K)-W(I,J,K))/DIEP(I)
     C            +DE(I-1,J,K)*(W(I,J,K)-W(I-1,J,K))
     C            /(MUEFF(I,J,K)*AREP(J,K)))
          ENDIF
*
*        No interface
*
         ELSE
          DTFDX = 0.5*((TF(I,J,K)-TF(I-1,J,K))/DIEP(I-1)
     C            +(TF(I+1,J,K)-TF(I,J,K))/DIEP(I))
          DTSDX = 0.5*((TS(I,J,K)-TS(I-1,J,K))/DIEP(I-1)
     C            +(TS(I+1,J,K)-TS(I,J,K))/DIEP(I))
          DUDX = 0.5*((U(I,J,K)-U(I-1,J,K))/DIEP(I-1)
     C           +(U(I+1,J,K)-U(I,J,K))/DIEP(I))
          DVDX = 0.5*((V(I,J,K)-V(I-1,J,K))/DIEP(I-1)
     C           +(V(I+1,J,K)-V(I,J,K))/DIEP(I))
          DWDX = 0.5*((W(I,J,K)-W(I-1,J,K))/DIEP(I-1)
     C           +(W(I+1,J,K)-W(I,J,K))/DIEP(I))
         ENDIF
*
*        Cases for y derivatives:
*
*        Interface at north face  
*
         IF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,5)
     C      .AND. CVTYPE(I,J,K,5).NE.1) THEN
          DTFDY = 0.5*(TF(I,J,K)-TF(I,J-1,K))/DJNP(J-1)
          DTSDY = 0.5*(TS(I,J,K)-TS(I,J-1,K))/DJNP(J-1)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,5).EQ.0) THEN
            DTFDY = DTFDY+0.5*DNF(I,J,K)*(TF(I,J+1,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARNP(I,K))
            DTSDY = DTSDY+0.5*DNS(I,J,K)*(TF(I,J+1,K)
     C              -TS(I,J,K))/(KEFFY(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,5).EQ.3) THEN
            DTFDY = DTFDY+0.5*DNF(I,J,K)*(TS(I,J+1,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARNP(I,K))
            DTSDY = DTSDY+0.5*DNS(I,J,K)*(TS(I,J+1,K)
     C              -TS(I,J,K))/(KEFFY(I,J,K)*ARNP(I,K))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,5).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDY = DTFDY+0.5*(DNF(I,J,K)*(TF(I,J+1,K)
     C              -TF(I,J,K))+DNS(I,J,K)*(TS(I,J+1,K)
     C              -TF(I,J,K)))/(KEFF(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDY = DTFDY+0.5*(DNF(I,J,K)*(TF(I,J+1,K)
     C              -TS(I,J,K))+DNS(I,J,K)*(TS(I,J+1,K)
     C              -TS(I,J,K)))/(KEFF(I,J,K)*ARNP(I,K))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDY = DTFDY+0.5*DNF(I,J,K)*(TS(I,J+1,K)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDY = DTSDY+0.5*DNS(I,J,K)*(TF(I,J+1,K)
     C              -TS(I,J,K))/(KEFF(I,J,K)*ARNP(I,K))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDY = 0.5*((U(I,J,K)-U(I,J-1,K))/DJNP(J-1)
     C            +DN(I,J,K)*(U(I,J+1,K)-U(I,J,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
           DVDY = 0.5*((V(I,J,K)-V(I,J-1,K))/DJNP(J-1)
     C            +DN(I,J,K)*(V(I,J+1,K)-V(I,J,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
           DWDY = 0.5*((W(I,J,K)-W(I,J-1,K))/DJNP(J-1)
     C            +DN(I,J,K)*(W(I,J+1,K)-W(I,J,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
          ENDIF
*
*        Interface at south face  
*
         ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,4)
     C      .AND. CVTYPE(I,J,K,4).NE.1) THEN
          DTFDY = 0.5*(TF(I+1,J,K)-TF(I,J,K))/DJNP(J)
          DTSDY = 0.5*(TS(I+1,J,K)-TS(I,J,K))/DJNP(J)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,4).EQ.0) THEN
            DTFDY = DTFDY+0.5*DNF(I,J-1,K)*(TF(I,J,K)
     C              -TF(I,J-1,K))/(KEFF(I,J,K)*ARNP(I,K))
            DTSDY = DTSDY+0.5*DNS(I,J-1,K)*(TS(I,J,K)
     C              -TF(I,J-1,K))/(KEFFY(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,4).EQ.3) THEN
            DTFDY = DTFDY+0.5*DNF(I,J-1,K)*(TF(I,J,K)
     C              -TS(I,J-1,K))/(KEFF(I,J,K)*ARNP(I,K))
            DTSDY = DTSDY+0.5*DNS(I,J-1,K)*(TS(I,J,K)
     C              -TS(I,J-1,K))/(KEFFY(I,J,K)*ARNP(I,K))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,4).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDY = DTFDY+0.5*(DNF(I,J-1,K)*(TF(I,J,K)
     C              -TF(I,J-1,K))+DNS(I,J-1,K)*(TF(I,J,K)
     C              -TS(I,J-1,K)))/(KEFF(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDY = DTFDY+0.5*(DNF(I,J-1,K)*(TS(I,J,K)
     C              -TF(I,J-1,K))+DNS(I,J-1,K)*(TS(I,J,K)
     C              -TS(I,J-1,K)))/(KEFF(I,J,K)*ARNP(I,K))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDY = DTFDY+0.5*DNF(I,J-1,K)*(TF(I,J,K)
     C              -TS(I,J-1,K))/(KEFF(I,J,K)*ARNP(I,K))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDY = DTSDY+0.5*DNS(I,J-1,K)*(TS(I,J,K)
     C              -TF(I,J-1,K))/(KEFF(I,J,K)*ARNP(I,K))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDY = 0.5*((U(I,J+1,K)-U(I,J,K))/DJNP(J)
     C            +DN(I,J-1,K)*(U(I,J,K)-U(I,J-1,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
           DVDY = 0.5*((V(I,J+1,K)-V(I,J,K))/DJNP(J)
     C            +DN(I,J-1,K)*(V(I,J,K)-V(I,J-1,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
           DWDY = 0.5*((W(I,J+1,K)-W(I,J,K))/DJNP(J)
     C            +DN(I,J-1,K)*(W(I,J,K)-W(I,J-1,K))
     C            /(MUEFF(I,J,K)*ARNP(I,K)))
          ENDIF
*
*        No interface
*
         ELSE
          DTFDY = 0.5*((TF(I,J,K)-TF(I,J-1,K))/DJNP(J-1)
     C            +(TF(I,J+1,K)-TF(I,J,K))/DJNP(J))
          DTSDY = 0.5*((TS(I,J,K)-TS(I,J-1,K))/DJNP(J-1)
     C            +(TS(I,J+1,K)-TS(I,J,K))/DJNP(J))
          DUDY = 0.5*((U(I,J,K)-U(I,J-1,K))/DJNP(J-1)
     C           +(U(I,J+1,K)-U(I,J,K))/DJNP(J))
          DVDY = 0.5*((V(I,J,K)-V(I,J-1,K))/DJNP(J-1)
     C           +(V(I,J+1,K)-V(I,J,K))/DJNP(J))
          DWDY = 0.5*((W(I,J,K)-W(I,J-1,K))/DJNP(J-1)
     C           +(W(I,J+1,K)-W(I,J,K))/DJNP(J))
         ENDIF
*
*        Cases for z derivatives:
*
*        Interface at top face  
*
         IF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,7)
     C      .AND. CVTYPE(I,J,K,7).NE.1) THEN
          DTFDZ = 0.5*(TF(I,J,K)-TF(I,J,K-1))/DKTP(K-1)
          DTSDZ = 0.5*(TS(I,J,K)-TS(I,J,K-1))/DKTP(K-1)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,7).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K)*(TF(I,J,K+1)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARTP(I,J))
            DTSDZ = DTSDZ+0.5*DTS(I,J,K)*(TF(I,J,K+1)
     C              -TS(I,J,K))/(KEFFZ(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,7).EQ.3) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K)*(TS(I,J,K+1)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARTP(I,J))
            DTSDZ = DTSDZ+0.5*DTS(I,J,K)*(TS(I,J,K+1)
     C              -TS(I,J,K))/(KEFFZ(I,J,K)*ARTP(I,J))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,7).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*(DTF(I,J,K)*(TF(I,J,K+1)
     C              -TF(I,J,K))+DTS(I,J,K)*(TS(I,J,K+1)
     C              -TF(I,J,K)))/(KEFF(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDZ = DTFDZ+0.5*(DTF(I,J,K)*(TF(I,J,K+1)
     C              -TS(I,J,K))+DTS(I,J,K)*(TS(I,J,K+1)
     C              -TS(I,J,K)))/(KEFF(I,J,K)*ARTP(I,J))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K)*(TS(I,J,K+1)
     C              -TF(I,J,K))/(KEFF(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDZ = DTSDZ+0.5*DTS(I,J,K)*(TF(I,J,K+1)
     C              -TS(I,J,K))/(KEFF(I,J,K)*ARTP(I,J))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDZ = 0.5*((U(I,J,K)-U(I,J,K-1))/DKTP(K-1)
     C            +DT(I,J,K)*(U(I,J,K+1)-U(I,J,K))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
           DVDZ = 0.5*((V(I,J,K)-V(I,J,K-1))/DKTP(K-1)
     C            +DT(I,J,K)*(V(I,J,K+1)-V(I,J,K))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
           DWDZ = 0.5*((W(I,J,K)-W(I,J,K-1))/DKTP(K-1)
     C            +DT(I,J,K)*(W(I,J,K+1)-W(I,J,K))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
          ENDIF
*
*        Interface at bottom face  
*
         ELSEIF(CVTYPE(I,J,K,1).NE.CVTYPE(I,J,K,6)
     C      .AND. CVTYPE(I,J,K,6).NE.1) THEN
          DTFDZ = 0.5*(TF(I,J,K+1)-TF(I,J,K))/DKTP(K)
          DTSDZ = 0.5*(TS(I,J,K+1)-TS(I,J,K))/DKTP(K)
          IF(CVTYPE(I,J,K,1).EQ.2) THEN
           IF(CVTYPE(I,J,K,6).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K-1)*(TF(I,J,K)
     C              -TF(I,J,K-1))/(KEFF(I,J,K)*ARTP(I,J))
            DTSDZ = DTSDZ+0.5*DTS(I,J,K-1)*(TS(I,J,K)
     C              -TF(I,J,K-1))/(KEFFZ(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,6).EQ.3) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K-1)*(TF(I,J,K)
     C              -TS(I,J,K-1))/(KEFF(I,J,K)*ARTP(I,J))
            DTSDZ = DTSDZ+0.5*DTS(I,J,K-1)*(TS(I,J,K)
     C              -TS(I,J,K-1))/(KEFFZ(I,J,K)*ARTP(I,J))
           ENDIF
          ELSEIF(CVTYPE(I,J,K,6).EQ.2) THEN
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*(DTF(I,J,K-1)*(TF(I,J,K)
     C              -TF(I,J,K-1))+DTS(I,J,K-1)*(TF(I,J,K)
     C              -TS(I,J,K-1)))/(KEFF(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDZ = DTFDZ+0.5*(DTF(I,J,K-1)*(TS(I,J,K)
     C              -TF(I,J,K-1))+DTS(I,J,K-1)*(TS(I,J,K)
     C              -TS(I,J,K-1)))/(KEFF(I,J,K)*ARTP(I,J))
           ENDIF
          ELSE
           IF(CVTYPE(I,J,K,1).EQ.0) THEN
            DTFDZ = DTFDZ+0.5*DTF(I,J,K-1)*(TF(I,J,K)
     C              -TS(I,J,K-1))/(KEFF(I,J,K)*ARTP(I,J))
           ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
            DTSDZ = DTSDZ+0.5*DTS(I,J,K-1)*(TS(I,J,K)
     C              -TF(I,J,K-1))/(KEFF(I,J,K)*ARTP(I,J))
           ENDIF
          ENDIF
          IF(CVTYPE(I,J,K,1).NE.3) THEN
           DUDZ = 0.5*((U(I,J,K+1)-U(I,J,K))/DKTP(K)
     C            +DT(I,J,K-1)*(U(I,J,K)-U(I,J,K-1))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
           DVDZ = 0.5*((V(I,J,K+1)-V(I,J,K))/DKTP(K)
     C            +DT(I,J,K-1)*(V(I,J,K)-V(I,J,K-1))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
           DWDZ = 0.5*((W(I,J,K+1)-W(I,J,K))/DKTP(K)
     C            +DT(I,J,K-1)*(W(I,J,K)-W(I,J,K-1))
     C            /(MUEFF(I,J,K)*ARTP(I,J)))
          ENDIF
*
*        No interface
*
         ELSE
          DTFDZ = 0.5*((TF(I,J,K)-TF(I,J,K-1))/DKTP(K-1)
     C            +(TF(I,J,K+1)-TF(I,J,K))/DKTP(K))
          DTSDZ = 0.5*((TS(I,J,K)-TS(I,J,K-1))/DKTP(K-1)
     C            +(TS(I,J,K+1)-TS(I,J,K))/DKTP(K))
          DUDZ = 0.5*((U(I,J,K)-U(I,J,K-1))/DKTP(K-1)
     C           +(U(I,J,K+1)-U(I,J,K))/DKTP(K))
          DVDZ = 0.5*((V(I,J,K)-V(I,J,K-1))/DKTP(K-1)
     C           +(V(I,J,K+1)-V(I,J,K))/DKTP(K))
          DWDZ = 0.5*((W(I,J,K)-W(I,J,K-1))/DKTP(K-1)
     C           +(W(I,J,K+1)-W(I,J,K))/DKTP(K))
         ENDIF
*
*        Calculation of entropy generation
*
         IF(CVTYPE(I,J,K,1).EQ.0) THEN
          DISS = VISC(I,J,K)*(2.0*(DUDX**2+DVDY**2+DWDZ**2)
     C           +(DVDX+DUDY)**2+(DWDY+DVDZ)**2+(DUDZ+DWDX)**2)
          SGENF(I,J,K) = KEFF(I,J,K)*(DTFDX**2+DTFDY**2+DTFDZ**2)
     C                   /(TF(I,J,K)**2+10.0**(-20))
     C                   +DISS/(TF(I,J,K)+10.0**(-20))
          SGENS(I,J,K) = 0.0
         ELSEIF(CVTYPE(I,J,K,1).EQ.2) THEN
          VEL = (U(I,J,K)**2+V(I,J,K)**2+W(I,J,K)**2)**0.5
          DISS = VISC(I,J,K)*(2.0*(DUDX**2+DVDY**2+DWDZ**2)
     C           +(DVDX+DUDY)**2+(DWDY+DVDZ)**2
     C           +(DUDZ+DWDX)**2)/PRSTY(I,J,K)+(VISC(I,J,K)/KPERM(I,J,K)
     C         +RHO(I,J,K)*CFORCH(I,J,K)*VEL/(KPERM(I,J,K)**0.5))*VEL**2
          SGENF(I,J,K) = KEFF(I,J,K)*(DTFDX**2+DTFDY**2+DTFDZ**2)
     C                   /(TF(I,J,K)**2+10.0**(-20))
     C                   +DISS/(TF(I,J,K)+10.0**(-20))
     C                   +HSF(I,J,K)*SPECSA(I,J,K)*((TS(I,J,K)
     C                   -TF(I,J,K))**2)/(TF(I,J,K)*TS(I,J,K)
     C                   +10.0**(-20))
          SGENS(I,J,K) = (KEFFX(I,J,K)*DTSDX**2
     C                   +KEFFY(I,J,K)*DTSDY**2
     C                   +KEFFZ(I,J,K)*DTSDZ**2)
     C                   /(TS(I,J,K)**2+10.0**(-20))
         ELSEIF(CVTYPE(I,J,K,1).EQ.3) THEN
          DISS = 0.0
          SGENF(I,J,K) = 0.0
          SGENS(I,J,K) = KEFF(I,J,K)*(DTSDX**2+DTSDY**2+DTSDZ**2)
     C                   /(TS(I,J,K)**2+10.0**(-20))
         ENDIF
         SGEN(I,J,K) = SGENF(I,J,K)+SGENS(I,J,K)
         RATIOS = RATIOS+SGENS(I,J,K)*VOLP(I,J,K)
         SGENT = SGENT+SGEN(I,J,K)*VOLP(I,J,K)
         RATIO(I,J,K) = DISS/(TF(I,J,K)*SGEN(I,J,K)-DISS
     C                  +10.0**(-20))
         RATIOT = RATIOT+DISS*VOLP(I,J,K)/(TF(I,J,K)+10.0**(-20))
         RATIOB = RATIOB+VOLP(I,J,K)*SGEN(I,J,K)
*
   5    CONTINUE
  10   CONTINUE
  20  CONTINUE
*
      RATIOT = RATIOT/RATIOB
      RATIOS = RATIOS/SGENT
*
      RETURN
      END