*
*     Note this file contains 2 subroutines, WEIGHT and PROFILE
************************************************************************
*
      SUBROUTINE WEIGHT(ALFAE,ALFAN,ALFAT, ME,MN,MT,
     C                  DE,DN,DT,CFACT,DFACT,
     C                  CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the advection weighting factors
*     for estimating fluxes through the east and north faces of
*     each control volume.
*
*     ALFAE(ID,JD,KD) advection factor for east face; output
*     ALFAN(ID,JD,KD) advection factor for north face; output
*     ALFAT(ID,JD,KD) advection factor for top face; output
*     ME(ID,JD,KD) normal mass flux for east face; input
*     MN(ID,JD,KD) normal mass flux for north face; input
*     MN(ID,JD,KD) normal mass flux for top face; input
*     DE(ID,JD,KD) diffusion coefficient for east face; input
*     DN(ID,JD,KD) diffusion coefficient for north face; input
*     DT(ID,JD,KD) diffusion coefficient for top face; input
*     CFACT,DFACT factors for advective and diffusive terms; input
*     CVTYPE(ID,JD,KD,NNB+1) indicates type of each CV; input
*     IB,IE first and last interior indices in i; input
*     JB,JE first and last interior indices in j; input
*     KB,KE first and last interior indices in k; input
*     ID,JD,KD array dimensions; input     
*     NNB number of neighbors for interior nodes; input
*
*     Notes: 1) It is best to calculate the weighting factors with 
*            a subroutine PRFL. This routine will make calls to
*            PRFL for each e and n point.
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 ALFAE(ID,JD,KD),ALFAN(ID,JD,KD),ALFAT(ID,JD,KD)
      REAL*8 CFACT(ID,JD,KD,NNB+1),DFACT(ID,JD,KD,NNB+1),CONV
      REAL*8 ME(ID,JD,KD),MN(ID,JD,KD),MT(ID,JD,KD)
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER I,J,K
      INTEGER JEM1,IEM1
*
      CALL NULLM(ALFAE, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)
      CALL NULLM(ALFAN, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)
      CALL NULLM(ALFAT, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)
*
      JEM1=JE-1
      IEM1=IE-1
*
      DO 20 I=IB,IE
       DO 15 K=KB,KE
        ALFAN(I,JB-1,K)=1.0
        DO 10 J=JB,JE-1
          IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,5).NE.3) THEN
           CONV=0.5*(CFACT(I,J,K,1)/DFACT(I,J,K,5)
     C               +CFACT(I,J,K,5)/DFACT(I,J+1,K,4))*MN(I,J,K)
           CALL PRFL(ALFAN(I,J,K),CONV,DN(I,J,K))
          ENDIF
 10     CONTINUE
        ALFAN(I,JE,K)=1.0
 15    CONTINUE
 20   CONTINUE
      DO 40 J=JB,JE
       DO 35 K=KB,KE
        ALFAE(IB-1,J,K)=1.0
        DO 30 I=IB,IE-1
          IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,3).NE.3) THEN
           CONV=0.5*(CFACT(I,J,K,1)/DFACT(I,J,K,3)
     C               +CFACT(I,J,K,3)/DFACT(I+1,J,K,2))*ME(I,J,K)
           CALL PRFL(ALFAE(I,J,K),CONV,DE(I,J,K))   
          ENDIF   
 30     CONTINUE
        ALFAE(IE,J,K)=1.0
 35    CONTINUE
 40   CONTINUE
      DO 60 I=IB,IE
       DO 55 J=JB,JE
        ALFAT(I,J,KB-1)=1.0
        DO 50 K=KB,KE-1
          IF(CVTYPE(I,J,K,1).NE.3 .AND. CVTYPE(I,J,K,7).NE.3) THEN
           CONV=0.5*(CFACT(I,J,K,1)/DFACT(I,J,K,7)
     C               +CFACT(I,J,K,7)/DFACT(I,J,K+1,6))*MT(I,J,K)
           CALL PRFL(ALFAT(I,J,K),CONV,DT(I,J,K))   
          ENDIF   
 50     CONTINUE
        ALFAT(I,J,KE)=1.0
 55    CONTINUE
 60   CONTINUE
*
      RETURN
      END
*     
***********************************************************************
*
      SUBROUTINE PRFL(ALFA,M,D)
*
*     Subroutine to calculate the weighting factors alfa and beta
*     at a point.
*
*     ALFA advection weighting factor; output
*
*     M    mass flux for face at point; input
*     D    diffusion coefficient for face at point; input
*
*     Notes: 1) The weighting factors should be based on the local
*               Peclet number.
*                    PE = M / ( D + 1.0E-20 )
*
*            2) Use the polynomial approximation to the exponential
*               weighting scheme.
*
*            3) The weighting factor varies from: -1.0 < ALFA < 1.0
*
************************************************************************
*
      IMPLICIT NONE
      REAL*8 ALFA,M,D,PE
*
      PE=M/(D+1.0E-20)
      ALFA=M/ABS(M+1.0E-20)*PE**2/(5.0+PE**2)
*
      RETURN
      END
