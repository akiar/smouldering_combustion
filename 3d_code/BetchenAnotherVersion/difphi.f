*
*          file difphi.f
*
************************************************************************
*
      SUBROUTINE DIFPHI(DE,DN,DT, GAMAX,GAMAY,GAMAZ,AREP,ARNP,ARTP,
     C                  DIEP,DJNP,DKTP,DISE,DISN,DIST,VFRAC,
     C                  EQTYPE,CVTYPE,IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB)
*
*     Subroutine to calculate the diffusion coefficients for  
*     normal conduction fluxes through the east and north faces of
*     each control volume. These coefficients are used to calculate
*     the finite volume equations and to assess the relative role
*     of diffusion and advection.
*
*     DE(ID,JD,KD) diffusion coefficient for east face; output
*     DN(ID,JD,KD) diffusion coefficient for north face; output
*     DT(ID,JD,KD) diffusion coefficient for top face; output
*     GAMAX(ID,JD,KD) diffusivity, multiplied by density or
*                      thermal inertia per volume, in x; input
*     GAMAY(ID,JD,KD) diffusivity, multiplied by density or
*                      thermal inertia per volume, in y; input                      
*     GAMAZ(ID,JD,KD) diffusivity, multiplied by density or
*                      thermal inertia per volume, in z; input
*     AREP(JD,KD) CV area of face at e point; input
*     ARNP(ID,KD) CV area of face at n point; input
*     ARTP(ID,JD) CV area of face at t point; input
*     DIEP(ID) distance from P to E through e; input
*     DJNP(JD) distance from P to N through n; input
*     DKTP(JD) distance from P to T through t; input
*     DISE(ID) distance in i from P to e face; input
*     DISN(JD) distance in j from P to n face; input
*     DIST(KD) distance in k from P to t face; input
*     VFRAC(ID,JD,KD) volume fraction of relevant constituent; input
*     CVTYPE(ID,JD,KD,NNB+1) control volume type; input
*     IB,IE first and last interior indices in i; input
*     JB,JE first and last interior indices in j; input
*     KB,KE first and last interior indices in j; input
*     ID,JD,KD array dimensions; input     
*     NNB number of neighbors for interior nodes; input
*
*     Note:  EQTYPE=0: Energy equation
*            EQTYPE=1: Momentum equation
*
***********************************************************************
*
      IMPLICIT NONE
      REAL*8 DE(ID,JD,KD),DN(ID,JD,KD),DT(ID,JD,KD)
      REAL*8 AREP(JD,KD),ARNP(ID,KD),ARTP(ID,JD)
      REAL*8 DIEP(ID),DJNP(JD),DKTP(KD),DISE(ID),DISN(JD),DIST(KD)
      REAL*8 GAMAX(ID,JD,KD),GAMAY(ID,JD,KD),GAMAZ(ID,JD,KD)
      REAL*8 VFRAC(ID,JD,KD)
      INTEGER IB,IE,JB,JE,KB,KE,ID,JD,KD,NNB,CVTYPE(ID,JD,KD,NNB+1)
      INTEGER EQTYPE,I,J,K
*
      CALL NULLM(DE, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)
      CALL NULLM(DN, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)            
      CALL NULLM(DT, IB-1,IE+1,JB-1,JE+1,KB-1,KE+1,ID,JD,KD)       
*
*  Diffusion coefficients for momentum equation
*
      IF(EQTYPE.EQ.1) THEN
*
*  Boundary faces
*
       DO 10 I=IB,IE
        DO 5 K=KB,KE
         J=JB-1
         DN(I,J,K)= GAMAY(I,J,K)*ARNP(I,K)/DJNP(J)       
  5     CONTINUE
 10    CONTINUE
       DO 20 J=JB,JE
        DO 15 K=KB,KE
         I=IB-1
         DE(I,J,K)= GAMAX(I,J,K)*AREP(J,K)/DIEP(I)
 15     CONTINUE
 20    CONTINUE
       DO 30 I=IB,IE
        DO 25 J=JB,JE
         K=KB-1
         DT(I,J,K)= GAMAZ(I,J,K)*ARTP(I,J)/DKTP(K)
 25     CONTINUE
 30    CONTINUE 
*
*  Interior faces
*
       DO 60 K=KB,KE
        DO 50 J=JB,JE
         DO 40 I=IB,IE
          IF(CVTYPE(I,J,K,1).EQ.3) THEN
*         
*          Solid CVs
* 
           IF(CVTYPE(I,J,K,3).EQ.0 .OR. CVTYPE(I,J,K,3).EQ.2) THEN
            DE(I,J,K)=GAMAX(I+1,J,K)*AREP(J,K)/DISE(I+1)
           ENDIF
*          
           IF(CVTYPE(I,J,K,5).EQ.0 .OR. CVTYPE(I,J,K,5).EQ.2) THEN           
            DN(I,J,K)=GAMAY(I,J+1,K)*ARNP(I,K)/DISN(J+1)
           ENDIF
*          
           IF(CVTYPE(I,J,K,7).EQ.0 .OR. CVTYPE(I,J,K,7).EQ.2) THEN           
            DT(I,J,K)=GAMAZ(I,J,K+1)*ARTP(I,J)/DIST(K+1)
           ENDIF           
*           
          ELSE
*
*          Clear fluid and porous CVs
*         
           IF(CVTYPE(I,J,K,3).EQ.3) THEN
            DE(I,J,K)=GAMAX(I,J,K)*AREP(J,K)/DISE(I)
           ELSE
            DE(I,J,K)=AREP(J,K)/(DISE(I)/GAMAX(I,J,K)
     C                         + DISE(I+1)/GAMAX(I+1,J,K))
           ENDIF
*          
           IF(CVTYPE(I,J,K,5).EQ.3) THEN
            DN(I,J,K)=GAMAY(I,J,K)*ARNP(I,K)/DISN(J)
           ELSE
            DN(I,J,K)=ARNP(I,K)/(DISN(J)/GAMAY(I,J,K)
     C                         + DISN(J+1)/GAMAY(I,J+1,K))
           ENDIF
*          
           IF(CVTYPE(I,J,K,7).EQ.3) THEN
            DT(I,J,K)=GAMAZ(I,J,K)*ARTP(I,J)/DIST(K)
           ELSE
            DT(I,J,K)=ARTP(I,J)/(DIST(K)/GAMAZ(I,J,K)
     C                         + DIST(K+1)/GAMAZ(I,J,K+1))
           ENDIF           
*          
          ENDIF
 40      CONTINUE
 50     CONTINUE
 60    CONTINUE
*
*  Boundary faces
*
       DO 62 I=IB,IE
        DO 61 K=KB,KE
         J=JE
         DN(I,J,K)= GAMAY(I,J+1,K)*ARNP(I,K)/DJNP(J)         
 61     CONTINUE
 62    CONTINUE
       DO 64 J=JB,JE
        DO 63 K=KB,KE
         I=IE
         DE(I,J,K)= GAMAX(I+1,J,K)*AREP(J,K)/DIEP(I)         
 63     CONTINUE
 64    CONTINUE
       DO 66 I=IB,IE
        DO 65 J=JB,JE
         K=KE
         DT(I,J,K)= GAMAZ(I,J,K+1)*ARTP(I,J)/DKTP(K)                  
 65     CONTINUE
 66    CONTINUE  
*
*  Diffusion coefficients for energy equations
*
      ELSEIF(EQTYPE.EQ.0) THEN
*
*  Boundary faces
*
       DO 70 I=IB,IE
        DO 67 K=KB,KE
         J=JB-1
         DN(I,J,K)= GAMAY(I,J,K)*ARNP(I,K)/DJNP(J)
 67     CONTINUE
 70    CONTINUE
       DO 80 J=JB,JE
        DO 75 K=KB,KE
         I=IB-1
         DE(I,J,K)= GAMAX(I,J,K)*AREP(J,K)/DIEP(I)
 75     CONTINUE
 80    CONTINUE
       DO 90 I=IB,IE
        DO 85 J=JB,JE
         K=KB-1
         DT(I,J,K)= GAMAZ(I,J,K)*ARTP(I,J)/DKTP(K)
 85     CONTINUE
 90    CONTINUE
 
*
*  Interior faces
*
       DO 120 K=KB,KE
        DO 110 J=JB,JE
         DO 100 I=IB,IE
*
           IF(CVTYPE(I,J,K,1).EQ.2 .AND. (CVTYPE(I,J,K,3).EQ.0
     C        .OR. CVTYPE(I,J,K,3).EQ.3)) THEN
            DE(I,J,K)=AREP(J,K)/(DISE(I)/GAMAX(I,J,K)
     C                + DISE(I+1)/(VFRAC(I,J,K)*GAMAX(I+1,J,K)))
           ELSEIF(CVTYPE(I,J,K,1).NE.2 
     C            .AND. CVTYPE(I,J,K,3).EQ.2) THEN
            DE(I,J,K)=AREP(J,K)/(DISE(I)/(VFRAC(I+1,J,K)
     C                *GAMAX(I,J,K)) + DISE(I+1)/GAMAX(I+1,J,K))          
           ELSE
            DE(I,J,K)=AREP(J,K)/(DISE(I)/GAMAX(I,J,K)
     C                         + DISE(I+1)/GAMAX(I+1,J,K))
           ENDIF
*
           IF(CVTYPE(I,J,K,1).EQ.2 .AND. (CVTYPE(I,J,K,5).EQ.0
     C        .OR. CVTYPE(I,J,K,5).EQ.3)) THEN
            DN(I,J,K)=ARNP(I,K)/(DISN(J)/GAMAY(I,J,K)
     C                + DISN(J+1)/(VFRAC(I,J,K)*GAMAY(I,J+1,K)))
           ELSEIF(CVTYPE(I,J,K,1).NE.2 
     C            .AND. CVTYPE(I,J,K,5).EQ.2) THEN
            DN(I,J,K)=ARNP(I,K)/(DISN(J)/(VFRAC(I,J+1,K)
     C                *GAMAY(I,J,K)) + DISN(J+1)/GAMAY(I,J+1,K))          
           ELSE
            DN(I,J,K)=ARNP(I,K)/(DISN(J)/GAMAY(I,J,K)
     C                         + DISN(J+1)/GAMAY(I,J+1,K))
           ENDIF
*
           IF(CVTYPE(I,J,K,1).EQ.2 .AND. (CVTYPE(I,J,K,7).EQ.0
     C        .OR. CVTYPE(I,J,K,7).EQ.3)) THEN
            DT(I,J,K)=ARTP(I,J)/(DIST(K)/GAMAZ(I,J,K)
     C                + DIST(K+1)/(VFRAC(I,J,K)*GAMAZ(I,J,K+1)))
           ELSEIF(CVTYPE(I,J,K,1).NE.2 
     C            .AND. CVTYPE(I,J,K,7).EQ.2) THEN
            DT(I,J,K)=ARTP(I,J)/(DIST(K)/(VFRAC(I,J,K+1)
     C                *GAMAZ(I,J,K)) + DIST(K+1)/GAMAZ(I,J,K+1))          
           ELSE
            DT(I,J,K)=ARTP(I,J)/(DIST(K)/GAMAZ(I,J,K)
     C                         + DIST(K+1)/GAMAZ(I,J,K+1))
           ENDIF
*
 100     CONTINUE
 110    CONTINUE
 120   CONTINUE      
*
*  Boundary faces
*
       DO 122 I=IB,IE
        DO 121 K=KB,KE
         J=JE
         DN(I,J,K)= GAMAY(I,J+1,K)*ARNP(I,K)/DJNP(J)         
 121    CONTINUE
 122   CONTINUE
       DO 124 J=JB,JE
        DO 123 K=KB,KE
         I=IE
         DE(I,J,K)= GAMAX(I+1,J,K)*AREP(J,K)/DIEP(I)         
 123    CONTINUE
 124   CONTINUE
       DO 126 I=IB,IE
        DO 125 J=JB,JE
         K=KE
         DT(I,J,K)= GAMAZ(I,J,K+1)*ARTP(I,J)/DKTP(K)                  
 125    CONTINUE
 126   CONTINUE      
      ENDIF 
*      
      RETURN
      END